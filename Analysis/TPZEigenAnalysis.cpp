//
// Created by Francisco Teixeira Orlandini on 11/17/17.

#include "TPZEigenAnalysis.h"
#include "TPZEigenSolver.h"
#include "TPZKrylovEigenSolver.h"
#include "TPZSpStructMatrix.h"
#include "pzysmp.h"
#include "pzsysmp.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "TPZMatGeneralisedEigenVal.h"
#include "TPZMaterial.h"
#include "pzcmesh.h"
#include "TPZBndCond.h" // para identificar boundary conditions

#include "pzlog.h"

#ifdef PZ_LOG
static TPZLogger logger("pz.analysis");
static TPZLogger loggerError("pz.analysis.error");
#endif

TPZEigenAnalysis::TPZEigenAnalysis()
: TPZRegisterClassId(&TPZEigenAnalysis::ClassId), TPZAnalysis()
{
}

TPZEigenAnalysis::TPZEigenAnalysis(TPZCompMesh *mesh,
                                   bool mustOptimizeBandwidth, std::ostream &out)
: TPZAnalysis(mesh, mustOptimizeBandwidth, out)
{
}

TPZEigenAnalysis::TPZEigenAnalysis(TPZAutoPointer<TPZCompMesh> mesh,
                                   bool mustOptimizeBandwidth,
                                   std::ostream &out)
: TPZAnalysis(mesh, mustOptimizeBandwidth, out)
{
}

template<class TVar>
TPZEigenSolver<TVar> &TPZEigenAnalysis::EigenSolver()
{
  const auto tmp = dynamic_cast<TPZEigenSolver<TVar>*>(fSolver);
  if (fSolver && !tmp) {
    PZError << __PRETTY_FUNCTION__;
    PZError << " incompatible Solver type! Aborting\n";
    DebugStop();
  }
  return *tmp;
}

void TPZEigenAnalysis::SetSolver(const TPZSolver &solver)
{
  if (fSolver) delete fSolver;
  auto *tmpState  = dynamic_cast<const TPZEigenSolver<STATE>*> (&solver);
  auto *tmpCState = dynamic_cast<const TPZEigenSolver<CSTATE>*>(&solver);
  if (tmpState && fSolType == EReal) {
    fSolver = (TPZEigenSolver<STATE>*) solver.Clone();
    return;
  } else if (tmpCState && fSolType == EComplex) {
    fSolver = (TPZEigenSolver<CSTATE>*) solver.Clone();
    return;
  }
  PZError << __PRETTY_FUNCTION__ << " Incompatible types!\nAborting...\n";
  DebugStop();
}

void TPZEigenAnalysis::Assemble()
{
  if (fSolType == EReal) AssembleT<STATE>();
  else                   AssembleT<CSTATE>();
}

template<class TVar>
void TPZEigenAnalysis::AssembleT()
{
  // evita perturbar pivôs no Pardiso
  auto PardisoConfig = [](TPZAutoPointer<TPZMatrix<TVar>> mat)
  {
    auto ssp = TPZAutoPointerDynamicCast<TPZSYsmpMatrix<TVar>>(mat);
    auto sp  = TPZAutoPointerDynamicCast<TPZFYsmpMatrix<TVar>>(mat);
    TPZPardisoSolver<TVar>* pardiso{nullptr};
    if (ssp) pardiso = &(ssp->GetPardisoControl());
    else if (sp) pardiso = &(sp->GetPardisoControl());
    if (pardiso) {
      auto param = pardiso->GetParam();
      param[0] = 0; // não perturbar pivôs
    }
  };

  TPZFMatrix<TVar> dummyRhs;
  if (!fCompMesh) {
    std::stringstream sout;
    sout << __PRETTY_FUNCTION__ << "\nNo computational mesh found!\n";
    #ifdef PZ_LOG
    LOGPZ_ERROR(logger, sout.str().c_str());
    #else
    std::cout << sout.str().c_str() << std::endl;
    #endif
    return;
  }

  if (!this->fStructMatrix) {
    #ifdef USING_MKL
    std::cout << "Setting default struct matrix: sparse(non-symmetric)\n";
    TPZSpStructMatrix<TVar> defaultMatrix(fCompMesh);
    #else
    std::cout << "Setting default struct matrix: skyline(non-symmetric)\n";
    TPZSkylineNSymStructMatrix<TVar> defaultMatrix(fCompMesh);
    #endif
    this->SetStructuralMatrix(defaultMatrix);
  }

  fStructMatrix->SetComputeRhs(false);

  { // garante que há solver
    auto *eigSolver = &(this->EigenSolver<TVar>());
    if (!eigSolver) {
      std::cout << "Setting default solver: Krylov\n";
      constexpr int nev{10};
      constexpr int dimKrylov{100};
      TPZKrylovEigenSolver<TVar> defaultSolver;
      defaultSolver.SetNEigenpairs(nev);
      defaultSolver.SetKrylovDim(dimKrylov);
      std::cout << "Setting nev: " << nev << "\n";
      std::cout << "Setting krylov dim: " << dimKrylov << "\n";
      this->SetSolver(defaultSolver);
    }
  }

  const auto sz = fStructMatrix->EquationFilter().NActiveEquations();
  auto &eigSolver = this->EigenSolver<TVar>();

  // --------- Sinaliza materiais para montar A (K) ----------
  if (eigSolver.IsGeneralised()) {
    auto &materialVec = fCompMesh->MaterialVec();
    bool foundGen = false;
    for (auto &&item : materialVec) {
      auto *mat = item.second;
      if (dynamic_cast<TPZBndCond*>(mat)) continue; // ignora BCs
      if (auto *eigmat = dynamic_cast<TPZMatGeneralisedEigenVal*>(mat)) {
        eigmat->SetMatrixA();
        foundGen = true;
      }
    }
    if (!foundGen) {
      PZError << __PRETTY_FUNCTION__
      << "\nNo material implements TPZMatGeneralisedEigenVal.\nAborting...\n";
      DebugStop();
    }
  }

  // --------- Monta Matrix A ----------
  auto matA = eigSolver.MatrixA();
  if (matA && matA->Rows() == sz) {
    matA->Zero();
    fStructMatrix->Assemble(*matA, dummyRhs, fGuiInterface);
  } else {
    TPZAutoPointer<TPZMatrix<TVar>> mat =
    dynamic_cast<TPZMatrix<TVar>*>(fStructMatrix->CreateAssemble(dummyRhs, fGuiInterface));
    eigSolver.SetMatrixA(mat);
  }
  PardisoConfig(eigSolver.MatrixA());
  matA = eigSolver.MatrixA(); // refresh
  if (matA) {
    //std::cout << "=== Stiffness Matrix (A) ===\n";
    //matA->Print("Matrix A", std::cout, EFormatted);
  }

  // --------- Sinaliza materiais para montar B (M) ----------
  if (eigSolver.IsGeneralised()) {
    auto &materialVec = fCompMesh->MaterialVec();
    for (auto &&item : materialVec) {
      auto *mat = item.second;
      if (dynamic_cast<TPZBndCond*>(mat)) continue; // ignora BCs
      if (auto *eigmat = dynamic_cast<TPZMatGeneralisedEigenVal*>(mat)) {
        eigmat->SetMatrixB();
      }
    }

    // --------- Monta Matrix B ----------
    auto matB = eigSolver.MatrixB();
    if (matB && matB->Rows() == sz) {
      matB->Zero();
      fStructMatrix->Assemble(*matB, dummyRhs, fGuiInterface);
    } else {
      TPZAutoPointer<TPZMatrix<TVar>> mat =
      dynamic_cast<TPZMatrix<TVar>*>(fStructMatrix->CreateAssemble(dummyRhs, fGuiInterface));
      eigSolver.SetMatrixB(mat);
    }
    PardisoConfig(eigSolver.MatrixB());
    matB = eigSolver.MatrixB(); // refresh
    if (matB) {
      //std::cout << "=== Mass Matrix (B) ===\n";
      //matB->Print("Matrix B", std::cout, EFormatted);
    }
  }
}

void TPZEigenAnalysis::Solve()
{
  if (fSolType == EReal) SolveT<STATE>();
  else                   SolveT<CSTATE>();
}

template<class TVar>
void TPZEigenAnalysis::SolveT()
{
  const auto nEq        = fCompMesh->NEquations();
  const auto nReducedEq = fStructMatrix->NReducedEquations();
  const bool isReduced  = nReducedEq != nEq;
  auto &eigSolver = this->EigenSolver<TVar>();

  fEigenvalues.Resize(nReducedEq);

  if (ComputeEigenvectors()) {
    TPZFMatrix<CTVar> *eigvectors = &fEigenvectors;
    if (isReduced) eigvectors = new TPZFMatrix<CTVar>;
    const auto nev = eigSolver.NEigenpairs();
    eigvectors->Redim(nReducedEq, nev);
    eigSolver.Solve(fEigenvalues, *eigvectors);
    if (isReduced) {
      fEigenvectors.Redim(nEq, nev);
      fStructMatrix->EquationFilter().Scatter(*eigvectors, fEigenvectors);
      delete eigvectors;
    }
  } else {
    eigSolver.Solve(fEigenvalues);
  }
}

TPZFMatrix<CSTATE>
TPZEigenAnalysis::GetEigenvectors() const
{
  if (!ComputeEigenvectors()) {
    std::cout << __PRETTY_FUNCTION__
    << "\nWARNING:No eigenvectors were computed\n";
  }
  return fEigenvectors;
}

TPZVec<CSTATE>
TPZEigenAnalysis::GetEigenvalues() const
{
  return fEigenvalues;
}

int TPZEigenAnalysis::ClassId() const
{
  return Hash("TPZEigenAnalysis") ^ (TPZAnalysis::ClassId() << 1);
}

void TPZEigenAnalysis::Write(TPZStream &buf, int withclassid) const
{
  TPZAnalysis::Write(buf, withclassid);
  buf.Write(fEigenvalues);
  fEigenvectors.Write(buf, withclassid);
  buf.Write(fCalcVectors);
}

void TPZEigenAnalysis::Read(TPZStream &buf, void *context)
{
  TPZAnalysis::Read(buf, context);
  buf.Read(fEigenvalues);
  fEigenvectors.Read(buf, context);
  buf.Read(fCalcVectors);
}

#define INSTANTIATE_TEMPLATES(TVar) \
template TPZEigenSolver<TVar> &TPZEigenAnalysis::EigenSolver<TVar>();

INSTANTIATE_TEMPLATES(STATE)
INSTANTIATE_TEMPLATES(CSTATE)
#undef INSTANTIATE_TEMPLATES

template class TPZRestoreClass<TPZEigenAnalysis>;
