// main.cpp — sweep em nRef × {Galerkin,Nystrom} × {Lapack,Krylov}
// Gera results_all.csv com (asm, solver, tempos, erros) e VTK do modo 0
// somente no último nível (se habilitado).

#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzgeoquad.h"
#include "pzintel.h"

#include "TPZEigenAnalysis.h"
#include "TPZKrylovEigenSolver.h"
#include "TPZLapackEigenSolver.h"
#include "pzskylstrmatrix.h"

#include "TPZVTKGeoMesh.h"
#include "TPZMatKLKernel.h"
#include "pzdoublestrmatriz.h"
#include "TPZBndCond.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <string>
#include <vector>
#include <iomanip> // std::fixed, std::setprecision
// =================== CONFIGURAÇÕES EXPLÍCITAS =================== //
static const int   kMatId        = 1;
static const int   kPOrder       = 2;
static const int   kRef0         = 0;
static const int   kRefMax       = 3;
static const bool  kWriteVTKLast = true; // VTK só no último nível

using CAsm = pzdoublestrmatriz<STATE>::ECAssembly;
enum class ESolver { Lapack, Krylov };

// Quais montagens C comparar (ordem = ordem de relatório)
static const std::vector<CAsm> kAsmList    = { CAsm::Galerkin, CAsm::Nystrom };
// Quais solvers comparar
static const std::vector<ESolver> kSolverList = { ESolver::Lapack, ESolver::Krylov };

static const char* AsmName(CAsm a){
        switch (a){ case CAsm::Galerkin: return "Galerkin"; case CAsm::Nystrom: return "Nystrom"; }
        return "Unknown";
}
static const char* SolverName(ESolver s){ return s==ESolver::Lapack ? "Lapack" : "Krylov"; }
// ================================================================ //

// ------------- geometria 3x3 (+ refinos) -------------
static TPZAutoPointer<TPZGeoMesh>
CreateGeoMeshMathematicaLike(int matId, int nRef)
{
        TPZAutoPointer<TPZGeoMesh> gmesh = new TPZGeoMesh;
        gmesh->SetDimension(2);

        const double m = 1.0/6.0;
        const double X[16][2] = {
                {-0.5,-0.5},{-0.5,-m},{-0.5, m},{-0.5, 0.5},
                {-m ,-0.5},{-m ,-m},{-m , m},{-m , 0.5},
                { m ,-0.5},{ m ,-m},{ m , m},{ m , 0.5},
                { 0.5,-0.5},{ 0.5,-m},{ 0.5, m},{ 0.5, 0.5}
        };
        gmesh->NodeVec().Resize(16);
        for (int i=0;i<16;i++){
                TPZManVector<REAL,3> c(3,0.); c[0]=X[i][0]; c[1]=X[i][1];
                gmesh->NodeVec()[i].Initialize(c, *gmesh);
        }
        const int conn[9][4] = {
                { 1, 5, 6, 2}, { 2, 6, 7, 3}, { 3, 7, 8, 4},
                { 5, 9,10, 6}, { 6,10,11, 7}, { 7,11,12, 8},
                { 9,13,14,10}, {10,14,15,11}, {11,15,16,12}
        };
        for (int e=0;e<9;e++){
                TPZManVector<int64_t,4> nodes(4);
                nodes[0]=conn[e][0]-1; nodes[1]=conn[e][1]-1;
                nodes[2]=conn[e][2]-1; nodes[3]=conn[e][3]-1;
                int64_t idx;
                gmesh->CreateGeoElement(EQuadrilateral, nodes, matId, idx);
        }
        gmesh->BuildConnectivity();

        for (int r=0;r<nRef;r++){
                const int nel = gmesh->NElements();
                TPZManVector<TPZGeoEl*> sub;
                for (int i=0;i<nel;i++) if (auto *gel = gmesh->ElementVec()[i]) gel->Divide(sub);
        }
        return gmesh;
}

// ------------- CompMesh H1 com kernel + solução exata -------------
static TPZCompMesh*
BuildCompMesh(TPZAutoPointer<TPZGeoMesh> gmesh, int porder, int matId)
{
        auto *cmesh = new TPZCompMesh(gmesh);
        cmesh->SetDimModel(gmesh->Dimension());
        cmesh->SetDefaultOrder(porder);
        cmesh->SetAllCreateFunctionsContinuous();

        const REAL Lx=1., Ly=1.;
        auto KernelFn = [Lx,Ly](const TPZVec<REAL>& x, const TPZVec<REAL>& y)->STATE{
                const REAL dx=x[0]-y[0], dy=x[1]-y[1];
                return std::exp(-std::fabs(dx)/Lx - std::fabs(dy)/Ly);
        };

        auto *mat = new TPZMatKLKernel(matId, 2, KernelFn);
        mat->SetId(matId);

        // modo 0 analítico
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

// ------------- linha para relatório -------------
struct Row {
        std::string asmname, solver;
        int ref=0, nel=0, ndof=0;
        double asm_t=0.0, solve_t=0.0, l2=0.0, energy=0.0, eig0=0.0;
};



using CAsm = pzdoublestrmatriz<STATE>::ECAssembly;


int main()
{
        std::cout << "===== RUN (sem argc/argv) =====\n"
        << "p          = " << kPOrder << "\n"
        << "ref range  = [" << kRef0 << "," << kRefMax << "]\n\n";

        std::vector<Row> results;
        results.reserve(kAsmList.size()*kSolverList.size()*(kRefMax-kRef0+1));

        const std::string csvname = "results_all.csv";
        std::ofstream csv(csvname);
        if (!csv) { std::cerr << "ERRO: não consegui abrir " << csvname << "\n"; return 1; }
        // Formato consistente para o Python
        csv.setf(std::ios::fixed);
        csv << "asm,solver,ref,nel,ndof,asm_s,solve_s,L2,Energy,eig0\n";
        std::cout << ">> Escrevendo CSV em " << csvname << "\n";

        for (CAsm asmChoice : kAsmList) {                 // <<<<<< CHAVES EXPLÍCITAS
                for (ESolver sol : kSolverList)               // <<<<<< CHAVES EXPLÍCITAS
                {
                        std::cout << "\n==== Combo: C=" << AsmName(asmChoice)
                        << " | Solver=" << SolverName(sol) << " ====\n";

                        for (int ref=kRef0; ref<=kRefMax; ++ref)
                        {
                                std::cout << "\n--- Refinement level " << ref << " ---\n";

                                TPZAutoPointer<TPZGeoMesh> gmesh = CreateGeoMeshMathematicaLike(kMatId, ref);
                                std::cout << "GeoMesh: nel=" << gmesh->NElements()
                                << "  nnodes=" << gmesh->NNodes() << "\n";

                                TPZCompMesh* cmesh = BuildCompMesh(gmesh, kPOrder, kMatId);
                                std::cout << "CompMesh: nelem=" << cmesh->NElements()
                                << "  ndof="  << cmesh->NEquations() << "\n";

                                TPZEigenAnalysis an(cmesh);
                                pzdoublestrmatriz<STATE> sm(cmesh);
                                sm.SetCAssembly(asmChoice);           // Galerkin ou Nystrom
                                an.SetStructuralMatrix(sm);

                                const int ndof = cmesh->NEquations();
                                if (sol==ESolver::Lapack) {
                                        TPZLapackEigenSolver<STATE> solver;
                                        solver.SetAsGeneralised(true);
                                        solver.SetEigenSorting(TPZEigenSort::AbsDescending);
                                        solver.SetNEigenpairs(ndof);
                                        an.SetSolver(solver);
                                        std::cout << "Solver: LAPACK (neigs=" << ndof << ")\n";
                                } else {
                                        TPZKrylovEigenSolver<STATE> solver;
                                        solver.SetAsGeneralised(true);
                                        solver.SetEigenSorting(TPZEigenSort::AbsDescending);
                                        solver.SetNEigenpairs(std::min(ndof,150));
                                        solver.SetKrylovDim(10*solver.NEigenpairs());
                                        solver.SetTolerance(1e-10);
                                        an.SetSolver(solver);
                                        std::cout << "Solver: KRYLOV (neigs=" << std::min(ndof,150)
                                        << ", krylovDim=" << 10*solver.NEigenpairs() << ")\n";
                                }

                                Row row; row.asmname=AsmName(asmChoice); row.solver=SolverName(sol);
                                row.ref=ref; row.ndof=ndof; row.nel=cmesh->NElements();

                                // ---- montagem ----
                                std::cout << "Assembling (C+B)..." << std::flush;
                                auto t0 = std::chrono::steady_clock::now();
                                an.Assemble();
                                auto t1 = std::chrono::steady_clock::now();
                                row.asm_t = std::chrono::duration<double>(t1-t0).count();
                                std::cout << " done in " << std::scientific << row.asm_t << " s\n";

                                // ---- solve ----
                                std::cout << "Solving EVP..." << std::flush;
                                auto t2 = std::chrono::steady_clock::now();
                                an.Solve();
                                auto t3 = std::chrono::steady_clock::now();
                                row.solve_t = std::chrono::duration<double>(t3-t2).count();
                                std::cout << " done in " << std::scientific << row.solve_t << " s\n";

                                TPZVec<CSTATE> vals = an.GetEigenvalues();
                                row.eig0 = vals.size() ? vals[0].real() : 0.0;
                                std::cout << "  top eigenvalue (Re) = " << std::setprecision(8) << row.eig0 << "\n";

                                // ---- carregar modo 0 (VTK só no último nível, se pedido) ----
                                TPZStack<std::string> scalars, vectors;
                                scalars.Push("Solution");  scalars.Push("ExactSolution");
                                scalars.Push("Error");     scalars.Push("ErrorSquared");
                                vectors.Push("Gradient");  vectors.Push("ExactGradient");
                                vectors.Push("ErrorGrad");

                                if (!kWriteVTKLast || ref==kRefMax) {
                                        std::cout << "PostProcess (mode 0, normalize by integral)..." << std::flush;
                                        an.PostProcessMode(0, 0, scalars, vectors,
                                                           TPZEigenAnalysis::EEigPart::Real,
                                                           /*normalizeByIntegral=*/true);
                                        if (kWriteVTKLast && ref==kRefMax) std::cout << " (VTK: mode0.vtk)";
                                        std::cout << " ok\n";
                                } else {
                                        // Se o seu TPZEigenAnalysis tiver método “carregar sem VTK”, use-o aqui.
                                        an.PostProcessMode(0, 0, scalars, vectors,
                                                           TPZEigenAnalysis::EEigPart::Real,
                                                           /*normalizeByIntegral=*/true);
                                }

                                // ---- integra erros (materiais de volume) ----
                                std::cout << "Integrating errors..." << std::flush;
                                std::set<int> mats;
                                for (auto &it : cmesh->MaterialVec())
                                        if (it.second && dynamic_cast<TPZBndCond*>(it.second)==nullptr)
                                                mats.insert(it.first);

                                const STATE l2_sq     = an.Integrate("ErrorSquared", mats)[0];
                                const STATE h1semi_sq = an.Integrate("GradErrorSquared", mats)[0];
                                row.l2     = std::sqrt(l2_sq);
                                row.energy = std::sqrt(l2_sq + h1semi_sq);
                                std::cout << " L2=" << std::setprecision(6) << row.l2
                                << "  Energy=" << row.energy << "\n";

                                // ---- escreve e força flush no CSV ----
                                csv << row.asmname << "," << row.solver << ","
                                << row.ref << "," << row.nel << "," << row.ndof << ","
                                << std::setprecision(10) << row.asm_t << "," << row.solve_t << ","
                                << row.l2 << "," << row.energy << "," << row.eig0 << "\n";
                                csv.flush();

                                results.push_back(row);
                                delete cmesh; // gmesh (AutoPointer) é liberada ao final
                        } // ref
                }     // solver
        }         // asm

        csv.close();

        // ---- resumo no stdout ----
        std::cout << "\n============================================================\n";
        std::cout << "p = " << kPOrder << " | ref=[" << kRef0 << "," << kRefMax << "]\n";
        std::cout << "asm        solver   ref   nel    ndof     assemble[s]   solve[s]      ||e||_0      ||e||_E      eig0\n";
        std::cout << "----------------------------------------------------------------------------------------------------------------\n";
        for (auto &r : results){
                std::cout <<std::fixed<< std::left << std::setw(10) << r.asmname
                << std::setw(9)  << r.solver
                << std::right
                << std::setw(4)  << r.ref
                << std::setw(7)  << r.nel
                << std::setw(8)  << r.ndof
                << std::setw(13) << std::fixed << std::setprecision(5) << r.asm_t
                << std::setw(12) << r.solve_t
                << std::setw(14) << r.l2
                << std::setw(14) << r.energy
                << std::setw(14) << r.eig0 << "\n";
        }
        std::cout << "----------------------------------------------------------------------------------------------------------------\n";
        std::cout << "Resultados salvos em " << csvname << "\n";
        return 0;
}
