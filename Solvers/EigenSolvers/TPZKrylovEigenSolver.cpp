#include "TPZKrylovEigenSolver.h"
#include "TPZLapackEigenSolver.h"
#include "TPZSimpleTimer.h"

#include <limits>
#include <cmath> // std::isfinite

template<class TVar>
void TPZKrylovEigenSolver<TVar>::SetTarget(TVar target)
{
  TPZEigenSolver<TVar>::SetTarget(target);
  fUserTarget = true;
  AdjustTargetST();
}

template<class TVar>
void TPZKrylovEigenSolver<TVar>::AdjustTargetST()
{
  auto st =
  dynamic_cast<TPZSTShiftOrigin<TVar>*>(this->SpectralTransform().operator->());
  if (st) {
    st->SetShift(this->Target());
  }
}

template<class TVar>
int TPZKrylovEigenSolver<TVar>::SolveImpl(TPZVec<CTVar> &w,
                                          TPZFMatrix<CTVar> &eigenVectors,
                                          bool computeVectors)
{
  if (this->NEigenpairs() < 1) SetNEigenpairs(1);

  #ifndef USING_LAPACK
  PZError << __PRETTY_FUNCTION__;
  PZError << "\nERROR: NeoPZ was not linked against LAPACK. Aborting...\n";
  DebugStop();
  #endif

  TPZSimpleTimer total("Arnoldi Solver");

  const int nRows = this->MatrixA()->Rows();

  if (fUserTarget) AdjustTargetST();

  TPZAutoPointer<TPZMatrix<TVar>> arnoldiMat{nullptr};
  auto st = this->SpectralTransform();
  const bool hasST = (bool)st;

  if (st) {
    TPZSimpleTimer calcMat("ST Calculating matrix");
    if (this->IsGeneralised()) {
      arnoldiMat = st->CalcMatrix(this->MatrixA(), this->MatrixB());
    } else {
      arnoldiMat = st->CalcMatrix(this->MatrixA());
    }
  } else {
    // sem ST, usamos diretamente A; se generalizado, decompomos B para aplicar B^{-1}A
    arnoldiMat = this->MatrixA();
    if (this->IsGeneralised() && this->MatrixB()) {
      TPZSimpleTimer binvert("Factorizing B");
      if (this->MatrixB()->IsSymmetric()) this->MatrixB()->Decompose_LDLt();
      else                                this->MatrixB()->Decompose_LU();
    }
  }

  const int nWanted = this->NEigenpairs();
  if (KrylovDim() == -1) SetKrylovDim(10 * nWanted);
  const int krylovDim = std::min(KrylovDim(), nRows);

  TPZManVector<TPZAutoPointer<TPZFMatrix<TVar>>, 20> qVecs;
  TPZFNMatrix<400, TVar> H(krylovDim, krylovDim, 0.);

  const bool success = ArnoldiIteration(*arnoldiMat, qVecs, H);
  if (!success) return -1;

  // Resolva o EVP de Hessenberg de dimensão krylovDim
  TPZFNMatrix<400, CTVar> lapackEV(krylovDim, krylovDim, 0.);
  TPZManVector<CTVar> w_all; // todos os Ritz valores (dim = krylovDim)

  const int lapackres = [&H, &w_all, &lapackEV]() {
    TPZSimpleTimer lapacktimer("Hessenberg EVP");
    TPZLapackEigenSolver<TVar> lapack;
    return lapack.SolveHessenbergEigenProblem(H, w_all, lapackEV);
  }();
  if (lapackres) return lapackres;

  // Transforme os autovalores de volta se usou ST
  if (st) st->TransformEigenvalues(w_all);

  // Ordene e selecione os nWanted
  TPZManVector<int, 20> indices;
  this->SortEigenvalues(w_all, indices);

  TPZManVector<CTVar> w_sorted_n(nWanted, 0.);
  for (int i = 0; i < nWanted; i++) {
    w_sorted_n[i] = w_all[indices[i]];
  }
  w = w_sorted_n;

  if (!computeVectors) return lapackres;

  // Reconstituir autovetores aproximados: V(:,i) = sum_j Q_j * y_j,
  // onde y é o autovetor de H correspondente (coluna de lapackEV)
  eigenVectors.Redim(nRows, nWanted);
  {
    TPZSimpleTimer evTimer("Computing eigenvectors");
    for (int i = 0; i < nWanted; i++) { // qual autovetor de A
      const int il = indices[i];        // coluna em lapackEV
      CTVar* ev = &eigenVectors.g(0, i);
      // zera a coluna i antes de acumular
      for (int k = 0; k < nRows; k++) ev[k] = CTVar(0);

      for (int j = 0; j < krylovDim; j++) { // qual vetor de Q
        const auto lev = lapackEV(j, il);
        TVar* q = &qVecs[j]->g(0, 0);
        for (int k = 0; k < nRows; k++) {
          ev[k] += lev * (*q++);
        }
      }
    }
  }

  return lapackres;
}

template<class TVar>
int TPZKrylovEigenSolver<TVar>::SolveEigenProblem(TPZVec<CTVar> &w, TPZFMatrix<CTVar> &eigenVectors)
{
  return SolveImpl(w, eigenVectors, true);
}

template<class TVar>
int TPZKrylovEigenSolver<TVar>::SolveEigenProblem(TPZVec<CTVar> &w)
{
  TPZFMatrix<CTVar> eigenVectors;
  return SolveImpl(w, eigenVectors, false);
}

template<class TVar>
int TPZKrylovEigenSolver<TVar>::SolveGeneralisedEigenProblem(TPZVec<CTVar> &w,
                                                             TPZFMatrix<CTVar> &eigenVectors)
{
  return SolveImpl(w, eigenVectors, true);
}

template<class TVar>
int TPZKrylovEigenSolver<TVar>::SolveGeneralisedEigenProblem(TPZVec<CTVar> &w)
{
  TPZFMatrix<CTVar> eigenVectors;
  return SolveImpl(w, eigenVectors, false);
}

template<class TVar>
bool TPZKrylovEigenSolver<TVar>::ArnoldiIteration(
  const TPZMatrix<TVar> &Aeff,
  TPZVec<TPZAutoPointer<TPZFMatrix<TVar>>> &Q,
  TPZFMatrix<TVar> &H)
{
  if (KrylovDim() < 2) {
    fKrylovDim = 10;
  }

  const int nRows = Aeff.Rows();
  const bool hasST = (bool)this->SpectralTransform();

  // Sem ST e problema generalizado: aplicamos o operador M = B^{-1}A
  // utilizando a fatoração previamente feita de B
  const bool useBInverse = (!hasST && this->fIsGeneralised && this->MatrixB());

  const int n = std::min(fKrylovDim, nRows);
  std::cout << "Calculating Krylov subspace of dimension " << n << '\n';

  H.Redim(n, n);
  H.Zero();
  Q.Resize(n, nullptr);

  // vetor inicial não-nulo
  if (fKrylovVector.Rows() != nRows || fKrylovVector.Cols() != 1 || Norm(fKrylovVector) == RTVar(0)) {
    fKrylovVector.Redim(nRows, 1);
    for (int i = 0; i < nRows; i++) fKrylovVector(i, 0) = TVar(1);
  }

  for (int i = 0; i < n; i++) Q[i] = new TPZFMatrix<TVar>;

  // inicializa primeiro vetor
  const RTVar v0norm = Norm(fKrylovVector);
  if (v0norm == RTVar(0)) return false;
  *(Q[0]) = fKrylovVector * (TVar)(RTVar(1) / v0norm);

  TPZSimpleTimer arnoldiIteration("ArnoldiIteration");
  const auto tol = Tolerance();

  // Função para aplicar o operador efetivo
  auto ApplyOp = [&](const TPZFMatrix<TVar>& v) -> TPZFMatrix<TVar> {
    TPZFMatrix<TVar> y = Aeff * v;
    if (useBInverse) {
      if (this->MatrixB()->IsSymmetric()) {
        // resolve B x = y usando a fatoração LDLt
        this->MatrixB()->Solve_LDLt(&y); // supõe-se solver in-place
      } else {
        // resolve B x = y usando a fatoração LU
        this->MatrixB()->Solve_LU(&y);   // supõe-se solver in-place
      }
    }
    return y;
  };

  for (int k = 1; k < n + 1; k++) {
    // w = M * q_{k-1}
    TPZFMatrix<TVar> w = ApplyOp(*(Q[k - 1]));

    RTVar normW{0};
    int restarts = 0;
    constexpr int kMaxRestarts = 12;

    while (true) {
      // zere a coluna k-1 de H antes de acumular (caso de recomeço)
      for (int j = 0; j < k; j++) H.PutVal(j, k - 1, TVar(0));

      // Gram-Schmidt com re-ortogonalização
      for (int pass = 0; pass < 2; pass++) {
        for (int j = k - 1; j >= 0; --j) {
          const auto& qj = *(Q[j]);
          const auto h = Dot(w, qj);
          H.PutVal(j, k - 1, H.GetVal(j, k - 1) + h);
          w -= qj * h;
        }
      }

      normW = Norm(w);

      // saiu? (norma válida e acima do tol) ou último passo
      const bool normOK = std::isfinite(static_cast<double>(normW)) && normW > tol;
      if (normOK || k == n) break;

      // evita laço infinito: declarou breakdown
      if (++restarts >= kMaxRestarts) { normW = RTVar(0); break; }

      // restart com vetor simples determinístico e normaliza
      w.Redim(nRows, 1);
      for (int i = 0; i < nRows; ++i) w(i, 0) = TVar((i & 1) ? 1.0 : -1.0);
      const RTVar nw = Norm(w);
      if (nw > RTVar(0)) w *= (TVar)(RTVar(1) / nw);
      else break;
    }

    if (k < n) {
      H.PutVal(k, k - 1, normW);

      if (std::isfinite(static_cast<double>(normW)) && normW > tol) {
        // caminho normal
        w *= (TVar)(RTVar(1) / normW);
        (*(Q[k])) = std::move(w);
      } else {
        // breakdown: injeta um vetor seguro ortogonal a Q
        TPZFMatrix<TVar> z(nRows, 1, 0.);
        z(k % nRows, 0) = TVar(1.0); // base canônica
        // ortogonaliza contra Q
        for (int j = k - 1; j >= 0; --j) {
          const auto& qj = *(Q[j]);
          z -= qj * Dot(z, qj);
        }
        RTVar nz = Norm(z);
        if (std::isfinite(static_cast<double>(nz)) && nz > tol) {
          z *= (TVar)(RTVar(1) / nz);
          (*(Q[k])) = std::move(z);
        } else {
          // deflação: mantém H(k,k-1)=0 e segue
        }
      }
    }
  } // for k

  return true;
}

// Instanciações explícitas
template class TPZKrylovEigenSolver<float>;
template class TPZKrylovEigenSolver<double>;
template class TPZKrylovEigenSolver<std::complex<float>>;
template class TPZKrylovEigenSolver<std::complex<double>>;
