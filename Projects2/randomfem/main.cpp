// main.cpp — Galerkin + Krylov (p=1) e VTK do campo KL (lognormal)

#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzgeoquad.h"
#include "pzintel.h"

#include "TPZEigenAnalysis.h"
#include "TPZKrylovEigenSolver.h"
#include "pzdoublestrmatriz.h"          // sua matriz estrutural dupla (K/M)
#include "pzskylstrmatrix.h"            // symmetric skyline matrix storage
#include "TPZMatKLKernel.h"

#include "Elasticity/TPZMatElastic2DMem.h"
#include "TPZVTKGeoMesh.h"
#include "pzstack.h"
#include "pzcompel.h"
#include "pzstepsolver.h"
#include "pztransfer.h"
#include "Poisson/TPZMatPoisson.h"
#include "TPZGenGrid2D.h"
#include "tpzgeoelrefpattern.h"
#include "pzgeoelbc.h"                  // << FIX: BC geométrico
#include <TPZGeoLinear.h>

#include <fstream>
#include <algorithm>
#include <cmath>
#include <random>
#include <vector>

// ---------------- geometria: 1 quad, BC: bottom (linha) = -1, ponto (nó 2) = -2 ----------------
static TPZAutoPointer<TPZGeoMesh>
CreateGeoMeshMathematicaLike(int matId, int nRef)
{

        TPZAutoPointer<TPZGeoMesh> gmesh = new TPZGeoMesh;
        gmesh->SetDimension(2);

        const double X[4][2] = {
                {-0.5,-0.5},{0.5,-0.5},{0.5,0.5},{-0.5, 0.5}
        };
        gmesh->NodeVec().Resize(4);
        for (int i=0;i<4;i++){
                TPZManVector<REAL,3> c(3,0.); c[0]=X[i][0]; c[1]=X[i][1];
                gmesh->NodeVec()[i].Initialize(c, *gmesh);
        }

        TPZManVector<int64_t,4> nodes(4);
        nodes[0]=0; nodes[1]=1; nodes[2]=2; nodes[3]=3;
        int64_t idx=0;
        gmesh->CreateGeoElement(EQuadrilateral, nodes, matId, idx);


        idx=1;
        TPZVec <long> TopoLine ( 2 ); TopoLine[0] = 0; TopoLine[1] = 1;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( idx, TopoLine, - 1, *gmesh );

        idx=2;
        TPZVec <long> node ( 1 ); node[0]=2;
        new TPZGeoElRefPattern< pzgeom::TPZGeoPoint> ( idx, node, - 2, *gmesh );//top

        gmesh->BuildConnectivity();

        for (int r=0;r<nRef;r++){
                const int nel = gmesh->NElements();
                TPZManVector<TPZGeoEl*> sub;
                for (int i=0;i<nel;i++) if (auto *eg = gmesh->ElementVec()[i]) eg->Divide(sub);
        }
        return gmesh;
}

// ---------------- CompMesh H1 com kernel + exata (opcional) ----------------
static TPZCompMesh*
BuildCompMesh(TPZAutoPointer<TPZGeoMesh> gmesh, int porder, int matId)
{
        auto *cmesh = new TPZCompMesh(gmesh);
        cmesh->SetDimModel(2);
        cmesh->SetDefaultOrder(porder);
        cmesh->SetAllCreateFunctionsContinuous();

        const REAL Lx=1., Ly=1.;
        auto KernelFn = [Lx,Ly](const TPZVec<REAL>& x, const TPZVec<REAL>& y)->STATE{
                return std::exp(-std::fabs(x[0]-y[0])/Lx - std::fabs(x[1]-y[1])/Ly);
        };

        auto *mat = new TPZMatKLKernel(matId, 2, KernelFn);
        mat->SetId(matId);

        // modo 0 analítico (opcional)
        constexpr double A = 1.15021, k = 1.30654;
        mat->SetExact([A,k](const TPZVec<REAL>& x, STATE& u, TPZFMatrix<STATE>& du){
                const double xx=x[0], yy=x[1];
                u = (STATE)(A*std::cos(k*xx)*std::cos(k*yy));
                du.Resize(2,1);
                du(0,0) = (STATE)(-A*k*std::sin(k*xx)*std::cos(k*yy));
                du(1,0) = (STATE)(-A*k*std::cos(k*xx)*std::sin(k*yy));
        });

        cmesh->InsertMaterialObject(mat);
        cmesh->AutoBuild();
        cmesh->AdjustBoundaryElements();
        cmesh->CleanUpUnconnectedNodes();
        return cmesh;
}

// ---------------- Malha elástica (usa seu material com memória) ----------------
static TPZCompMesh*
CreateCMeshElastic ( TPZAutoPointer<TPZGeoMesh> gmesh )
{
        TPZCompMesh *cmesh = new TPZCompMesh ( gmesh );
        cmesh->SetDimModel(2);
        cmesh->SetDefaultOrder ( TPZCompEl::GetgOrder() );
        cmesh->SetAllCreateFunctionsContinuous();

        // Ex.: TPZMatElastic2DMem<SeuMem> (id, planeStress,nu, E)
        auto *mat = new TPZMatElastic2DMem<TKLPointMem>(1, /*planeStress=*/true, /*nu=*/0., /*E=*/1.);


        TPZFMatrix<STATE> val1 ( 2,2,0. );
        TPZVec<STATE>     val2 ( 2,0. );

        // -1: bottom “clamped”
        val2[0]=1.; val2[1]=1.;
        TPZBndCond *bcclamp = mat->CreateBC ( mat, -1, 3, val1, val2 );

        // -2: carga 1 N no nó 2 (y)
        val2[0]=0.; val2[1]=1.;
        TPZBndCond *bcload  = mat->CreateBC ( mat, -2, 1, val1, val2 );

        cmesh->InsertMaterialObject ( mat );
        cmesh->InsertMaterialObject ( bcclamp );
        cmesh->InsertMaterialObject ( bcload );

        cmesh->AutoBuild();
        cmesh->AdjustBoundaryElements();
        cmesh->CleanUpUnconnectedNodes();
        return cmesh;
}

int main()
{
        static const int   kMatId  = 1;
        static const int   kPOrder = 1;
        static const int   kRef    = 4;      // n. de refinamentos
        // parâmetros do campo lognormal
        static const STATE kMu     = 1.0;    // média alvo do E
        static const STATE kCov    = 0.30;   // coeficiente de variação σ/μ
        static const int   kMmax   = 30;     // nº de modos KL usados
        static const uint64_t kSeed = 20250820ULL;

        // --- malhas ---
        TPZAutoPointer<TPZGeoMesh> gmesh   = CreateGeoMeshMathematicaLike(kMatId, kRef);
        TPZCompMesh* cmesh                 = BuildCompMesh(gmesh, kPOrder, kMatId);
        TPZCompMesh* cmeshelastic          = CreateCMeshElastic(gmesh);

        // --- análise de autovalores: Galerkin + Krylov ---
        TPZEigenAnalysis an(cmesh);

        pzdoublestrmatriz<STATE> sm(cmesh);
        sm.SetCAssembly(pzdoublestrmatriz<STATE>::ECAssembly::Galerkin);
        an.SetStructuralMatrix(sm);

        const int ndof = cmesh->NEquations();

        TPZKrylovEigenSolver<STATE> solver;
        solver.SetAsGeneralised(true);
        solver.SetEigenSorting(TPZEigenSort::AbsDescending);
        solver.SetNEigenpairs(ndof);
        solver.SetKrylovDim(ndof);
        solver.SetTolerance(1e-10);
        an.SetSolver(solver);

        an.Assemble();
        an.Solve();

        // ---- VTK do modo 0 (real) ----
        TPZStack<std::string> scalars, vectors;
        scalars.Push("Solution");       // mantenha apenas variáveis garantidas
        scalars.Push("ExactSolution");
        vectors.Push("Gradient");
        vectors.Push("ExactGradient");

        an.DefineGraphMesh(2, scalars, vectors, "kl_mode0.vtk"); // << FIX: define graph mesh
        an.PostProcessMode(0, 0, scalars, vectors, TPZEigenAnalysis::EEigPart::Real, true);

        // ---------------- Problema elástico (estático) ----------------
        TPZLinearAnalysis anelastic(cmeshelastic);
        TPZSkylineStructMatrix<STATE> smelastic(cmeshelastic);
        anelastic.SetStructuralMatrix(smelastic);

        TPZStepSolver<STATE> direct;
        direct.SetDirect(ELDLt);                // << FIX: define método direto
        anelastic.SetSolver(direct);

        anelastic.Assemble();
        anelastic.Solve();

        TPZStack<std::string> scalars2, vectors2;
        // Se seu TPZMatElastic2DMem expõe "E" como variável, mantenha; senão, remova.
        scalars2.Push("E");
        vectors2.Push("Displacement");

        anelastic.DefineGraphMesh(2, scalars2, vectors2, "elasticplate.vtk");
        anelastic.PostProcess(0);

        // ---------------- Campo KL → G ~ N(0, C) → lognormal E ----------------
        TPZVec<CSTATE> vals = an.GetEigenvalues();
        TPZFMatrix<CSTATE> Evecs = an.GetEigenvectors(); // colunas = autovetores

        const int M = std::min((int)vals.size(), kMmax);

        TPZFMatrix<STATE> PHI_nodal(ndof, M, 0.0);
        TPZVec<STATE>     lambdas(M, 0.0);
        for (int k = 0; k < M; ++k) {
                lambdas[k] = (STATE)vals[k].real();
                for (int i = 0; i < ndof; ++i) {
                        PHI_nodal(i, k) = (STATE)Evecs.GetVal(i, k).real();
                }
        }

        // Amostras alfa_k ~ N(0, λ_k)
        std::mt19937_64 rng(kSeed);
        std::normal_distribution<double> N01(0.0, 1.0);
        TPZFMatrix<STATE> Alpha(M, 1, 0.0);
        for (int k=0;k<M;k++) Alpha(k,0) = std::sqrt((double)lambdas[k]) * (STATE)N01(rng);

        TPZFMatrix<STATE> G_dof;
        PHI_nodal.Multiply(Alpha, G_dof); // [ndof x 1]

        // Var[G] por dof
        TPZVec<STATE> VarG_dof(ndof, 0.0);
        for (int i = 0; i < ndof; ++i) {
                double acc = 0.0;
                for (int k = 0; k < M; ++k) {
                        const double vik = (double)PHI_nodal.GetVal(i, k);
                        acc += (double)lambdas[k] * vik * vik;
                }
                VarG_dof[i] = (STATE)acc;
        }

        // mapeamento p/ lognormal (média kMu e COV=kCov)
        const STATE xi  = std::sqrt(std::log(1.0 + kCov*kCov));
        const STATE lam = std::log(kMu) - 0.5*xi*xi;

        TPZFMatrix<STATE> E_dof(ndof, 1, 0.0);
        for (int i = 0; i < ndof; ++i) {
                const double g    = (double)G_dof.GetVal(i, 0);
                const double invS = (VarG_dof[i] > 0.0) ? 1.0/std::sqrt((double)VarG_dof[i]) : 0.0;
                const double ghat = g * invS;
                E_dof(i, 0) = (STATE)std::exp((double)lam + (double)xi * ghat);
        }

        // VTK do campo KL (carrega como "Solution" na mesma cmesh do KL)
        cmesh->LoadSolution(E_dof);
        TPZStack<std::string> scal, vec;
        scal.Push("Solution");                 // se seu material expõe "KLField", troque aqui
        an.DefineGraphMesh(2, scal, vec, "KL_field.vtk");
        an.PostProcess(0);

        // limpeza
        delete cmeshelastic;
        delete cmesh;
        return 0;
}
