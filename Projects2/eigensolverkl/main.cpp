// main.cpp — sweep em nRef com prints, tempos e erros
// Uso:
//   ./prog --solver=lapack --p=2 --minref=0 --maxref=3
//   ./prog --solver=krylov --p=2 --minref=0 --maxref=3
// Flags extras: --vtk-last  (gera VTK só no último nível)

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


#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <string>
#include <vector>

// ---------------- util ----------------
template <class T>
static T read_int_arg(int argc, char** argv, const std::string& key, T defval){
        const std::string k = "--"+key+"=";
        for (int i=1;i<argc;i++){ std::string a=argv[i]; if (a.rfind(k,0)==0) return static_cast<T>(std::stoi(a.substr(k.size()))); }
        return defval;
}
static bool has_flag(int argc, char** argv, const std::string& flag){
        const std::string f="--"+flag; for (int i=1;i<argc;i++) if (f==argv[i]) return true; return false;
}

// -------- geometria 3x3 + nRef níveis (retorna AutoPointer!) --------
static TPZAutoPointer<TPZGeoMesh> CreateGeoMeshMathematicaLike(int matId, int nRef = 1)
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
                gmesh->NodeVec()[i].Initialize(c,*gmesh);
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
                gmesh->CreateGeoElement(EQuadrilateral,nodes,matId,idx);
        }
        gmesh->BuildConnectivity();

        for (int r=0;r<nRef;r++){
                const int nel = gmesh->NElements();
                TPZManVector<TPZGeoEl*> sub;
                for (int i=0;i<nel;i++) if (auto *gel = gmesh->ElementVec()[i]) gel->Divide(sub);
        }
        return gmesh;
}

// -------- CompMesh H1 com kernel + solução exata --------
static TPZCompMesh* BuildCompMesh(TPZAutoPointer<TPZGeoMesh> gmesh, int porder, int matId)
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

// -------- solver flag --------
enum class ESolver { Lapack, Krylov };
static ESolver parse_solver(int argc, char** argv){
        for (int i=1;i<argc;i++){ std::string a=argv[i];
                if (a=="--solver=lapack") return ESolver::Lapack;
                if (a=="--solver=krylov") return ESolver::Krylov;
        }
        return ESolver::Lapack;
}

// -------- linha de resultados --------
struct Row { int ref=0, nel=0, ndof=0; double asm_t=0.0, solve_t=0.0, l2=0.0, energy=0.0, eig0=0.0; };

int main(int argc, char** argv)
{
        const int matId  = 1;
        const int porder = read_int_arg(argc,argv,"p",1);
        const int ref0   = read_int_arg(argc,argv,"minref",0);
        const int refMax = read_int_arg(argc,argv,"maxref",3);
        const ESolver which = parse_solver(argc,argv);
        const bool vtkLast = has_flag(argc,argv,"vtk-last");

        std::cout << "===== RUN =====\n"
        << "solver       = " << (which==ESolver::Lapack?"LAPACK":"KRYLOV") << "\n"
        << "p            = " << porder << "\n"
        << "ref range    = [" << ref0 << "," << refMax << "]\n\n";

        std::vector<Row> results;

        for (int ref=ref0; ref<=refMax; ++ref)
        {
                std::cout << "\n=== Refinement level " << ref << " ===\n";

                // geo + comp
                TPZAutoPointer<TPZGeoMesh> gmesh = CreateGeoMeshMathematicaLike(matId, ref);
                std::cout << "GeoMesh: nel=" << gmesh->NElements()
                << "  nnodes=" << gmesh->NNodes() << std::endl;

                TPZCompMesh* cmesh = BuildCompMesh(gmesh, porder, matId);
                std::cout << "CompMesh: nelem=" << cmesh->NElements()
                << "  ndof=" << cmesh->NEquations() << std::endl;

                TPZEigenAnalysis an(cmesh);
                pzdoublestrmatriz<STATE> sm(cmesh);
                an.SetStructuralMatrix(sm);

                const int ndof = cmesh->NEquations();
                if (which==ESolver::Lapack) {
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

                Row row; row.ref=ref; row.ndof=ndof; row.nel=cmesh->NElements();

                // montagem
                std::cout << "Assembling (A+B)..." << std::flush;
                auto t0 = std::chrono::steady_clock::now();
                an.Assemble();
                auto t1 = std::chrono::steady_clock::now();
                row.asm_t = std::chrono::duration<double>(t1-t0).count();
                std::cout << " done in " << std::scientific << row.asm_t << " s\n";


                // solve
                std::cout << "Solving EVP..." << std::flush;
                auto t2 = std::chrono::steady_clock::now();
                an.Solve();
                auto t3 = std::chrono::steady_clock::now();
                row.solve_t = std::chrono::duration<double>(t3-t2).count();
                std::cout << " done in " << std::scientific << row.solve_t << " s\n";

                TPZVec<CSTATE> vals = an.GetEigenvalues();
                row.eig0 = vals.size() ? vals[0].real() : 0.0;
                std::cout << "  top eigenvalue (Re) = " << std::setprecision(8) << row.eig0 << "\n";

                // pós: carrega modo 0 e normaliza por integral
                std::cout << "Loading mode 0 + normalize by integral..." << std::flush;
                TPZStack<std::string> scalars, vectors;
                scalars.Push("Solution"); scalars.Push("ExactSolution"); scalars.Push("Error"); scalars.Push("ErrorSquared");
                vectors.Push("Gradient"); vectors.Push("ExactGradient"); vectors.Push("ErrorGrad");
                an.PostProcessMode(0, 0, scalars, vectors,
                                   TPZEigenAnalysis::EEigPart::Real,
                                   /*normalizeByIntegral=*/true);
                if (vtkLast && ref==refMax) std::cout << " VTK: mode0.vtk";
                std::cout << " ok\n";

                // integra erros (somente materiais de volume)
                std::cout << "Integrating errors..." << std::flush;
                std::set<int> mats;
                for (auto &it : cmesh->MaterialVec())if (it.second && dynamic_cast<TPZBndCond*>(it.second)==nullptr) mats.insert(it.first);

                const STATE l2_sq     = an.Integrate("ErrorSquared", mats)[0];
                const STATE h1semi_sq = an.Integrate("GradErrorSquared", mats)[0];
                row.l2     = std::sqrt(l2_sq);
                row.energy = std::sqrt(l2_sq + h1semi_sq);
                std::cout << " L2=" << std::setprecision(6) << row.l2
                << "  Energy=" << row.energy << "\n";

                results.push_back(row);

                // IMPORTANTE: destrua APENAS a compMesh.
                // gmesh é TPZAutoPointer; o último AutoPointer (no cmesh) libera a geoMesh.
                delete cmesh;
        }

        // tabela
        std::cout << "\n============================================================\n";
        std::cout << "p = " << porder << "  solver = " << (which==ESolver::Lapack?"LAPACK":"KRYLOV") << "\n";
        std::cout << "ref   nel    ndof        asm[s]     solve[s]        ||e||_0        ||e||_E        eig0\n";
        std::cout << "-------------------------------------------------------------------------------------------\n";
        std::ofstream csv("results.csv"); csv << "ref,nel,ndof,asm_s,solve_s,L2,Energy,eig0\n";
        for (auto &r : results){
                std::cout << std::setw(3) << r.ref
                << std::setw(7) << r.nel
                << std::setw(8) << r.ndof
                << std::setw(13) << std::scientific << std::setprecision(3) << r.asm_t
                << std::setw(12) << r.solve_t
                << std::setw(14) << r.l2
                << std::setw(14) << r.energy
                << std::setw(14) << r.eig0 << "\n";
                csv << r.ref << "," << r.nel << "," << r.ndof << ","
                << std::setprecision(10) << r.asm_t << "," << r.solve_t << ","
                << r.l2 << "," << r.energy << "," << r.eig0 << "\n";
        }
        std::cout << "-------------------------------------------------------------------------------------------\n";
        std::cout << "Resultados salvos em results.csv\n";
        return 0;
}
