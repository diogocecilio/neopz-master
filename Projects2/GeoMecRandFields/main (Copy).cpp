
#include <iostream>
#include <fstream>
#include <set>
#include <vector>
#include <limits>
#include <algorithm>
#include <cmath>
#include <string>

// ---------- NeoPZ core ----------
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzvec.h"
#include "pzmanvector.h"
#include "pzgeoel.h"

#include "pzcompel.h"
#include "pzconnect.h"
#include "pzblock.h"
#include "pzintel.h"
#include "pzgeoelside.h"

// ---------- NeoPZ geo + refs ----------
#include "pzgeotriangle.h"

#include "TPZVTKGeoMesh.h"

// ---------- NeoPZ material elastoplástico ----------
#include "TPZElasticResponse.h"
#include "TPZElastoPlasticMem.h"
#include "TPZPlasticStepPV.h"
#include "TPZYCMohrCoulombPV.h"
#include "TPZMatElastoPlastic2D.h"

// ---------- NeoPZ análise ----------

#include "pzstepsolver.h"
#include <iostream>
#include <fstream>
#include <thread>
#include <vector>
#include "TPZFileStream.h"
#include <TPZBFileStream.h>

#include <fstream>
#include <iostream>
#include <fstream>
#include <thread>
#include <vector>
#include <mutex>
#include <iostream>
#include <fstream>
#include <unistd.h>
#include <sys/wait.h>
#include <memory>

#include "pzskylstrmatrix.h"
#include "pzelastoplasticanalysis.h"
#include "pznonlinanalysis.h"
#include "TPZEigenSolver.h"
#include "TPZKrylovEigenSolver.h"
#include "TPZLapackEigenSolver.h" // ou outro solver concreto
#include "pzdoublestrmatriz.h"
#include "pzskylstrmatrix.h"
#include <TPZSSpStructMatrix.h> //symmetric sparse matrix storage
#include <pzskylstrmatrix.h> //symmetric skyline matrix storage
#include <pzstepsolver.h> //for TPZStepSolver
#include <TPZSimpleTimer.h>
#include "TPZPardisoSolver.h"
// main.cpp
#include <random>
#include <cmath>
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <set>

#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzgeoquad.h"
#include "TPZGeoLinear.h"
#include "tpzgeoelrefpattern.h"
#include "pzintel.h"
#include "pzstack.h"

#include "TPZEigenAnalysis.h"
#include "TPZKrylovEigenSolver.h"
#include "pzdoublestrmatriz.h"
#include "pzskylstrmatrix.h"

#include "TPZMatKLKernel.h"
#include "Elasticity/TPZMatElastic2DMem.h"
#include "Elasticity/TPZElasticMem.h"

#include "pzpostprocanalysis.h"
#include "pzstepsolver.h"
#include "TPZBFileStream.h"
#include "TPZVTKGeoMesh.h"
#include "Projection/TPZL2Projection.h"
typedef TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse> TPlasticMC;
typedef TPZMatElastoPlastic2D<TPlasticMC, TPZElastoPlasticMem>   plasticmat;

using namespace std;

// ------------------------------------------------------------
// Pós-processo
// ------------------------------------------------------------
void PostProcessVariables(TPZStack<std::string>& scal, TPZStack<std::string>& vec);

void CreatePostProcessingMesh(TPZCompMesh* cmesh,TPZPostProcAnalysis* pproc,int matid);

void PostElastoplastic(TPZCompMesh* cmesh,const std::string& vtkfile,int matid);

TPZGeoMesh*  TriGMesh(int ref);

plasticmat*  CreateMaterial(REAL young, REAL poisson, REAL coes, REAL atrito,TPZManVector<REAL,3> bodyforce,int planestrain = 1, int matid = 1);

TPZCompMesh* CreateCMesh(TPZGeoMesh* gmesh, int pOrder, plasticmat* mat);


void InitializeMemory(TPZCompMesh* cmesh, REAL coesion, REAL atrito);

void ComputeElementDeformation(TPZCompMesh* cmesh, TPZVec<REAL>& fPlasticDeformSqJ2);

void DivideElementsAbove(TPZCompMesh* cmesh, REAL refineaboveval, std::set<int64_t>& out_newels);

void PRefineElementsAbove (TPZCompMesh* cmesh, REAL refineaboveval, std::set<int64_t>& out_newels,int porder );

bool Hrefine(TPZCompMesh* cmesh, REAL coes, REAL atrito, REAL refineAboveVal);

bool HPrefine(TPZCompMesh* cmesh, REAL coes, REAL atrito,REAL refineAboveVal,int porder);

REAL UyAtNode(TPZCompMesh* cmesh, REAL x, REAL y);

bool RunAndAccept(TPZCompMesh* cmesh,
                  REAL coes, REAL atrito, REAL factor, int& iters_out,bool initmem=true);

REAL AutoRefine(TPZCompMesh* cmesh, REAL coes, REAL atrito, REAL refineAboveVal,
                int max_refines, REAL fs_start, REAL fs_step, REAL fs_max);


REAL BisectionFS(TPZCompMesh* cmesh, REAL coes, REAL atrito,
                 REAL lo, REAL hi, REAL tol_fs_rel, int max_bis,bool initmem=true);

REAL HybridBracketedFS(TPZCompMesh* cmesh, REAL coes, REAL atrito,
                       REAL lo, REAL hi, REAL tol_fs_rel, int max_it, bool initmem=true);


REAL Solve(TPZCompMesh* cmesh,REAL coes,REAL atrito,bool initmem=true);

enum class FieldModel { Gaussian, Lognormal };
struct FieldSpec {
        REAL mu_x   = 1.0;          // média alvo do campo X
        REAL cov_x  = 0.0;          // coeficiente de variação (sigma_x / mu_x)
        FieldModel model = FieldModel::Gaussian;
};

REAL SolveStochastic(REAL coes,REAL atrito,bool initmem,const std::vector<TPZCompMesh*>& sources,TPZCompMesh* target,const std::vector<FieldSpec>& specs,int idxE, int idxNu);


// KL paramétrico em Lx, Ly
TPZCompMesh * BuildCompMeshKL_Param(TPZGeoMesh *gmesh, int porder, int matId,
                                           REAL Lx, REAL Ly);

TPZFMatrix<STATE> BuildPhiSqrtLambda(TPZCompMesh* cmesh,const TPZFMatrix<CSTATE>& eigenvectors,const TPZVec<CSTATE>&eigenvalues,int M);



// ------------------------------------------------------------------
// ComputeFieldMulti: lê N malhas-fonte (cada uma já com uma realização "z(x)")
// e escreve em memória elástica do alvo, montando X_i(x) com mu_x e cov_x.
// - Gaussian:   X_i = mu_x[i] + (cov_x[i]*mu_x[i]) * z_i
// - Lognormal:  X_i = exp( mu_ln[i] + sigma_ln[i] * z_i )
//   onde: sigma_ln = sqrt( ln(1 + COV_x^2) ), mu_ln = ln(mu_x) - 0.5*sigma_ln^2
// idxE/idxNu escolhem quais campos alimentam E e nu.
// ------------------------------------------------------------------
void ComputeFieldMulti(const std::vector<TPZCompMesh*>& sources,
                              TPZCompMesh* target,
                              const std::vector<FieldSpec>& specs,
                              int idxE = 0, int idxNu = 1);

TPZFMatrix<REAL> ComputeTheta (int M, int samples, uint32_t seed=12345);

void BuildFields(int pOrder ,int ref,int kMatId,REAL Lx,REAL Ly,int M,int NsampGen);

void ApplyLoad(TPZCompMesh* cmesh,
               REAL coes, REAL atrito, TPZManVector<REAL> factors);
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <filesystem>

// Lê a última linha não vazia do CSV e retorna o último índice s (ou -1 se não houver dados)
static inline long LastSampleIndexFromCSV(const std::string& path) {
        std::ifstream in(path);
        if (!in) return -1;

        std::string line, last;
        while (std::getline(in, line)) {
                if (!line.empty()) last = line;
        }
        if (last.empty()) return -1;

        // Se tiver cabeçalho "s,FS", ignore
        if (last.find("s,") == 0 || last == "s,FS") return -1;

        // Parse "s,FS"
        std::istringstream iss(last);
        std::string tok_s;
        if (!std::getline(iss, tok_s, ',')) return -1;
        try {
                return std::stol(tok_s);
        } catch (...) {
                return -1;
        }
}
// int main()
// {
//         // 1) Malha geométrica
//         int ref = 0;
//         TPZGeoMesh* gmesh = TriGMesh(ref);
//
//         {
//                 std::ofstream vtk1("antes.vtk");
//                 TPZVTKGeoMesh::PrintGMeshVTK(gmesh, vtk1, true);
//                 std::cout << "antes.vtk escrito.\n";
//         }
//         // 2) Material
//         REAL young   = 20000.;
//         REAL poisson = 0.49;
//         // REAL coes    = 50.;
//         // REAL atrito  = 20. * M_PI / 180.;
//         REAL coes    = 10.;
//         REAL atrito  = 30. * M_PI / 180.;
//         TPZManVector<REAL,3> bodyforce(3,0.0);
//         bodyforce[1] = -20.0;
//
//         plasticmat* mat = CreateMaterial(young, poisson, coes, atrito, bodyforce);
//         int pOrder = 2;
//         TPZCompMesh* cmesh = CreateCMesh(gmesh, pOrder, mat);
//         mat->SetBodyForce(bodyforce);
//         InitializeMemory(cmesh, coes, atrito);
//         auto* body = dynamic_cast<plasticmat*>(cmesh->FindMaterial(1));
//         if (!body) { std::cerr << "Material id=1 não encontrado.\n"; return 1; }
//
//         using Clock = std::chrono::steady_clock;
//
//         auto t0 = Clock::now();
//         REAL FS = Solve(cmesh, coes, atrito,false);
//         auto t1 = Clock::now();
//
//         std::chrono::duration<double> secs = t1 - t0;
//         auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
//
//         std::cout << std::fixed << std::setprecision(3)
//         << "[Timing] Solve: " << secs.count() << " s  (" << ms << " ms)\n";
//
//         {
//                 std::ofstream vtk1("gmeshtri_refined_preGI.vtk");
//                 TPZVTKGeoMesh::PrintGMeshVTK(cmesh->Reference(), vtk1, true);
//                 std::cout << "[VTK] gmeshtri_refined_preGI.vtk escrito.\n";
//         }
//         string vtk2="post_plasticity.vtk";
//         int matid=1;
//         PostElastoplastic(cmesh,vtk2,matid);
//
//         cout <<"FS = "<<FS<<endl;
//         int nloads = 10;
//         const REAL FS_target = FS ;          // mantém seu “+0.1”
//         TPZManVector<REAL> factors(nloads+1);   // 0 .. nloads (inclusivo)
//
//         for (int i = 0; i <= nloads; ++i) {
//                 factors[i] = FS_target * REAL(i) / REAL(nloads); // 0, Δ, 2Δ, …, FS_target
//                 //cout<< factors[i] <<endl;
//         }
//
//
//         // ApplyLoad(cmesh,coes, atrito,  factors);
//
//         return 0;
// }

static long LastS(const std::string& path){
        std::ifstream in(path);
        if(!in) return -1;
        std::string line, last;
        while(std::getline(in,line)) if(!line.empty()) last = line;
        if(last.empty() || last.rfind("sample",0)==0) return -1; // cabeçalho
        auto p = last.find(','); if(p==std::string::npos) return -1;
        try { return std::stol(last.substr(0,p)); } catch(...) { return -1; }
}

int main()
{
        int pOrder =2;
        int refield=1;
        int ref=1;
        int kMatId=1;
        REAL Lx=20.;
        REAL Ly=2.;
        int M=100;
        int NsampGen=10000;
        if(false)
        {
                BuildFields( pOrder , refield, kMatId, Lx, Ly, M,NsampGen);
        }else
        {
                // malhas fontes (para mapear campo -> elasticidade)
                TPZGeoMesh*  gmesh0       =  TriGMesh(refield);
                TPZGeoMesh*  gmesh1        = TriGMesh(refield);

                TPZCompMesh* cmeshFieldCoes  = BuildCompMeshKL_Param(gmesh0, pOrder, kMatId, 1.0, 1.0);
                TPZCompMesh* cmeshFieldAtrito = BuildCompMeshKL_Param(gmesh1, pOrder, kMatId, 1.0, 1.0);

                // 2) Material
                REAL young   = 20000.;
                REAL poisson = 0.49;

                REAL coes    = 10.;
                REAL atrito  = 30. * M_PI / 180.;
                TPZManVector<REAL,3> bodyforce(3,0.0);
                bodyforce[1] = -20.0;

                // lê hhat
                std::ostringstream fin; fin.setf(std::ios::fixed);
                fin << "hhat_" << std::setprecision(6) << Lx << "_" << Ly << ".bin";
                TPZFMatrix<REAL> hhatE,hhatMU;
                {
                        TPZBFileStream in;
                        in.OpenRead(fin.str());
                        hhatE.Read(in, 0);
                        hhatMU.Read(in, 0);
                }
                const int NsampRun = std::min(500, (int)std::min(hhatE.Cols(), hhatMU.Cols()));


                // nome do CSV (mantém seu padrão)
                std::ostringstream fout; fout.setf(std::ios::fixed);
                fout << "mc_results_" << std::setprecision(6) << Lx << "_" << Ly << ".csv";
                const std::string csvname = fout.str();

                // retoma do último sample salvo
                const long s_start = LastS(csvname) + 1;

                // abre em append (escreve cabeçalho se arquivo novo)
                std::ofstream outcsv(csvname, std::ios::app);
                if(outcsv.tellp()==0) outcsv << "sample,FS\n";
                outcsv.setf(std::ios::fixed); outcsv << std::setprecision(10);

                std::cout << "solving mc " << fin.str()
                << " retomando em s=" << s_start
                << " até < " << NsampRun << std::endl;

                double sum=0.0, sum2=0.0; int nacc=0; // estatísticas desta sessão
                TPZFMatrix<REAL> colE(hhatE.Rows(),1), colNu(hhatMU.Rows(),1);

                // *** LOOP COM RETOMADA ***
                for (long s = s_start; s < NsampRun; ++s)
                {
                        cout << "========== IMC ========== "<< s <<endl;
                        TPZGeoMesh*  gmesh2  = TriGMesh(ref);
                        plasticmat*  mat     = CreateMaterial(young, poisson, coes, atrito, bodyforce);
                        int pOrder = 2;
                        TPZCompMesh* cmesh   = CreateCMesh(gmesh2, pOrder, mat);
                        mat->SetBodyForce(bodyforce);

                        auto* body = dynamic_cast<plasticmat*>(cmesh->FindMaterial(1));
                        if(!body){ std::cerr << "Material id=1 não encontrado.\n"; return 1; }

                        InitializeMemory(cmesh, coes, atrito);
                        for (int i=0;i<hhatE.Rows(); ++i) colE(i,0)  = hhatE(i,(int)s);
                        for (int i=0;i<hhatMU.Rows();++i) colNu(i,0) = hhatMU(i,(int)s);

                        cmeshFieldCoes->LoadSolution(colE);  cmeshFieldCoes->LoadReferences();
                        cmeshFieldAtrito->LoadSolution(colNu); cmeshFieldAtrito->LoadReferences();

                        std::vector<TPZCompMesh*> sources = { cmeshFieldCoes, cmeshFieldAtrito };
                        std::vector<FieldSpec> specs = {
                                {coes,  0.30, FieldModel::Lognormal},  // -> E
                                {atrito,0.20, FieldModel::Lognormal}   // -> nu
                        };
                        ComputeFieldMulti(sources, cmesh, specs, /*idxE=*/0, /*idxNu=*/1);

                        using Clock = std::chrono::steady_clock;
                        auto t0 = Clock::now();
                        REAL FS = SolveStochastic(coes, atrito, false, sources, cmesh, specs, 0, 1);
                        auto t1 = Clock::now();
                        auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
                        std::cout << std::fixed << std::setprecision(3)
                        << "[Timing] Solve: " << (ms/1000.0) << " s  (" << ms << " ms)\n";

                        // (opcionais)
                        // { std::ofstream vtk1("gmeshtri_refined_preGI.vtk");
                        //   TPZVTKGeoMesh::PrintGMeshVTK(cmesh->Reference(), vtk1, true); }
                        // PostElastoplastic(cmesh, "post_plasticity.vtk", /*matid=*/1);

                        if(FS <= 10.0){
                                outcsv << s << "," << FS << "\n";
                                outcsv.flush();          // grava imediatamente (seguro pra retomada)
                                sum  += FS; sum2 += (double)FS*(double)FS; ++nacc;
                        }

                        delete cmesh;
                        delete gmesh2;
                }
                outcsv.close();

                // estatísticas da sessão corrente (o que rodou agora)
                if(nacc>0){
                        const double mean = sum/nacc;
                        const double var  = std::max(0.0, sum2/nacc - mean*mean);
                        const double stdv = std::sqrt(var);
                        std::cout << "[Lx="<<Lx<<", Ly="<<Ly<<"] Nsamp(sessão)="<<nacc
                        << " mean="<<mean<<" std="<<stdv
                        << " -> salvo em " << csvname << "\n";
                } else {
                        std::cout << "Nenhuma amostra válida adicionada nesta sessão.\n";
                }

                delete cmeshFieldCoes;
                delete cmeshFieldAtrito;
                delete gmesh0;
                delete gmesh1;
        }


        return 0;
}


// ============================================================
// Implementações
// ============================================================
REAL SolveStochastic(REAL coes,REAL atrito,bool initmem,const std::vector<TPZCompMesh*>& sources,TPZCompMesh* target,const std::vector<FieldSpec>& specs,int idxE, int idxNu)
{
        REAL lo=0.25;
        REAL hi=20.;
        REAL tol_fs_rel=0.01;
        int max_bis = 30;
        REAL FS=5.;

        int porder=target->GetDefaultOrder();
        int iters_out;
        //FS=  BisectionFS(cmesh, coes, atrito, lo,  hi, tol_fs_rel ,  max_bis);
        //FS=1.77;
        //RunAndAccept( cmesh,coes,  atrito,  FS,  iters_out);
        //Hrefine(cmesh,coes, atrito,0.01);
        int maxref=4;
        for ( int iref=1; iref<=maxref; iref++ ) {

                int neq=target->NEquations();
                std::cout << "\n[solve] ===== Refinamento # "<< iref <<" ====="<<" neq = " <<neq << "\n";
                //FS=  BisectionFS(target, coes, atrito, lo,  hi, tol_fs_rel ,  max_bis,initmem);
                ComputeFieldMulti(sources, target, specs, /*idxE=*/0, /*idxNu=*/1);
                FS=HybridBracketedFS(target, coes, atrito, lo,  hi, tol_fs_rel ,  max_bis,initmem);
                if(iref==maxref)
                {
                        break;
                }
               // Hrefine(target,coes, atrito,0.01);

                HPrefine(target,coes, atrito,0.01,porder);

                porder+=1;
        }
        //para pos-processar
        RunAndAccept( target,coes,  atrito,  FS,  iters_out,initmem);
        return FS;
}
void BuildFields(int pOrder ,int ref,int kMatId,REAL Lx,REAL Ly,int M,int NsampGen)
{



        TPZGeoMesh* gmesh = TriGMesh(ref);
        TPZCompMesh* cmeshKL = BuildCompMeshKL_Param(gmesh, pOrder, kMatId, Lx, Ly);

        TPZEigenAnalysis an(cmeshKL,false);
        pzdoublestrmatriz<REAL> sm(cmeshKL);
        sm.SetCAssembly(pzdoublestrmatriz<REAL>::ECAssembly::Galerkin);
        an.SetStructuralMatrix(sm);

        TPZKrylovEigenSolver<STATE> esolver;
        esolver.SetAsGeneralised(true);
        esolver.SetEigenSorting(TPZEigenSort::AbsDescending);
        esolver.SetNEigenpairs(cmeshKL->NEquations());
        esolver.SetKrylovDim(cmeshKL->NEquations());
        esolver.SetTolerance(1e-10);
        an.SetSolver(esolver);

        an.Assemble();
        an.Solve();
        TPZFMatrix<CSTATE> evecs = an.Eigenvectors();
        TPZVec<CSTATE>     evals = an.Eigenvalues();

        TPZFMatrix<STATE> PHI   = BuildPhiSqrtLambda(cmeshKL, evecs, evals, M);



        // Duas famílias independentes para Coes e Phi (poderia ser correlacionado se desejado)
        TPZFMatrix<REAL>  THETAE  = ComputeTheta(M, NsampGen, /*seed*/12345);
        TPZFMatrix<REAL>  THETAMU = ComputeTheta(M, NsampGen, /*seed*/54321);

        TPZFMatrix<REAL>  hhatE,hhatMU; // (ndof x NsampGen)
        PHI.Multiply(THETAE,  hhatE);
        PHI.Multiply(THETAMU, hhatMU);

        // nome estável (fixed, 6 casas)
        std::ostringstream oss; oss.setf(std::ios::fixed);
        oss << "hhat_" << std::setprecision(6) << Lx << "_" << Ly << ".bin";
        TPZBFileStream out;
        out.OpenWrite(oss.str());
        hhatE.Write(out,0);
        hhatMU.Write(out,0);

        delete cmeshKL;
        delete gmesh;
        std::cout << "Gerado: " << oss.str() << "\n";

}
// ============================================================
// Implementações
// ============================================================
REAL Solve(TPZCompMesh* cmesh,REAL coes,REAL atrito,bool initmem)
{
        REAL lo=0.5;
        REAL hi=10.;
        REAL tol_fs_rel=0.01;
        int max_bis = 20;
        REAL FS=1.;

        int porder=cmesh->GetDefaultOrder();
        int iters_out;
        //FS=  BisectionFS(cmesh, coes, atrito, lo,  hi, tol_fs_rel ,  max_bis);
        //FS=1.77;
        //RunAndAccept( cmesh,coes,  atrito,  FS,  iters_out);
        //Hrefine(cmesh,coes, atrito,0.01);
        int maxref=4;
        for ( int iref=1; iref<=maxref; iref++ ) {

                int neq=cmesh->NEquations();
                std::cout << "\n[solve] ===== Refinamento # "<< iref <<" ====="<<" neq = " <<neq << "\n";
               //FS=  BisectionFS(cmesh, coes, atrito, lo,  hi, tol_fs_rel ,  max_bis,initmem);
                FS=HybridBracketedFS(cmesh, coes, atrito, lo,  hi, tol_fs_rel ,  max_bis,initmem);
                if(iref==maxref)
                {
                       break;
                }
                //Hrefine(cmesh,coes, atrito,0.01);

                HPrefine(cmesh,coes, atrito,0.01,porder);
                porder+=1;
        }
        //para pos-processar
        RunAndAccept( cmesh,coes,  atrito,  FS,  iters_out,initmem);
        return FS;
}
bool RunAndAccept(TPZCompMesh* cmesh,
                  REAL coes, REAL atrito, REAL factor, int& iters_out,bool initmem)
{
        auto* body = dynamic_cast<plasticmat*>(cmesh->FindMaterial(1));
        body->SetLoadFactor(factor);
        if(initmem)InitializeMemory(cmesh, coes, atrito);
        cmesh->Solution().Zero();


        TPZElastoPlasticAnalysis anal(cmesh, std::cout,TPZElastoPlasticAnalysis::ELineSearch::QuadraticArmijo);

        if(false)
        {
                TPZSSpStructMatrix<STATE> SSpStructMatrix ( cmesh );
                SSpStructMatrix.SetNumThreads(12);
                anal.SetStructuralMatrix(SSpStructMatrix);
                TPZPardisoSolver<REAL> *pardiso = new TPZPardisoSolver<REAL>;
                anal.SetSolver ( *pardiso );
        }else{
                TPZSkylineStructMatrix<STATE> matskl(cmesh);
                matskl.SetNumThreads(16);
                anal.SetStructuralMatrix(matskl);
                TPZStepSolver<STATE> step; step.SetDirect(ELDLt);
                anal.SetSolver(step);
        }

        int iters=30;
        bool ok = anal.FindRoot(iters_out);
        //bool ok = anal.IterativeProcess(std::cout, 1.e-3, iters, true, false, iters_out);
        if (!ok) return false;
        anal.AcceptSolution(1);
        //cmesh->LoadSolution(anal.CumulativeSolution());
        return true;
}
void ApplyLoad(TPZCompMesh* cmesh,
                  REAL coes, REAL atrito, TPZManVector<REAL> factors)
{

        // parâmetros de controle
        int nloads = factors.size();
        REAL FS_target = factors[nloads-1];
        REAL fator_atual = 0.0;
        REAL passo_base  = FS_target / REAL(nloads); // passo médio de referência

        cmesh->SetDefaultOrder(4);
        TPZElastoPlasticAnalysis anal(cmesh, std::cout,TPZElastoPlasticAnalysis::ELineSearch::Dicotomic);
        auto* body = dynamic_cast<plasticmat*>(cmesh->FindMaterial(1));
        body->ResetMemory();
        InitializeMemory(cmesh, coes, atrito);
        if(true)
        {
                TPZSSpStructMatrix<STATE> SSpStructMatrix ( cmesh );
                SSpStructMatrix.SetNumThreads(12);
                anal.SetStructuralMatrix(SSpStructMatrix);
                TPZPardisoSolver<REAL> *pardiso = new TPZPardisoSolver<REAL>;
                anal.SetSolver ( *pardiso );
        }else{
                TPZSkylineStructMatrix<STATE> matskl(cmesh);
                matskl.SetNumThreads(16);
                anal.SetStructuralMatrix(matskl);
                TPZStepSolver<STATE> step; step.SetDirect(ELDLt);
                anal.SetSolver(step);
        }
        const std::string csv_path = "loadsweep.csv";
        std::ofstream csv(csv_path);
        csv << "step,factor,uy,iters,ok\n";
        csv << std::setprecision(15) << std::scientific;
       // int nloads=factors.size();
        REAL x=30.;
        REAL y=40.;
        REAL uy=0.;
        int counter=0;
        int old_iters_out=0;
        int iters_out;
        while( counter<100)
        {
                body->SetLoadFactor(fator_atual);


                bool ok = anal.IterativeProcess(std::cout, 1.e-3, 1000, true, false, iters_out);
                if(old_iters_out<iters_out)
                {
                        old_iters_out=iters_out;
                }
                if(!ok)break;
                TPZFMatrix<REAL> tempsol=anal.Solution();
                anal.AcceptSolution(0);

                cmesh->LoadSolution(tempsol);
                uy += UyAtNode(cmesh,  x,  y);
                csv << counter << "," << -uy << "," << fator_atual << "," << iters_out << "," << 1 << "\n";
                std::cout << "uy = " << uy << "  factor = " << fator_atual
                << "  iters = " << iters_out <<  "  old_iters_out = " << old_iters_out << " counter =" << counter<< std::endl;




                REAL peso = 2.5 / std::max(1, old_iters_out);  // 1/iters
                REAL delta = passo_base * peso  ;      //

                fator_atual += delta;
                //if (fator_atual > FS_target) fator_atual = FS_target;

                if (fator_atual >= FS_target) old_iters_out=300;
                counter++;
        }

}

bool HPrefine(TPZCompMesh* cmesh, REAL coes, REAL atrito,REAL refineAboveVal,int porder)
{
        const int nels_before = cmesh->NElements();
        TPZVec<REAL> defel;
        ComputeElementDeformation(cmesh, defel);

        std::set<int64_t> novosp;
        PRefineElementsAbove (cmesh,refineAboveVal, novosp, porder );

        std::set<int64_t> novos;
        DivideElementsAbove(cmesh, refineAboveVal, novos);

        const int nels_after = cmesh->NElements();
        std::cout << "[HRefine] nels: " << nels_before << " -> " << nels_after
        << "  (refinados: " << (int)novos.size() << ")\n";

        std::cout << "[PRefine] nels: " << (int)novosp.size() << "\n";

        if (nels_after > nels_before) {
                InitializeMemory(cmesh, coes, atrito);
                return true;
        } else {
                std::cout << "[PreRefine] sem novos refinamentos; fim.\n";
                return false;
        }
}


bool Hrefine(TPZCompMesh* cmesh, REAL coes, REAL atrito,REAL refineAboveVal)
{
        const int nels_before = cmesh->NElements();
        TPZVec<REAL> defel;
        ComputeElementDeformation(cmesh, defel);

        std::set<int64_t> novos;
        DivideElementsAbove(cmesh, refineAboveVal, novos);

        const int nels_after = cmesh->NElements();
        std::cout << "[PreRefine] nels: " << nels_before << " -> " << nels_after
        << "  (refinados: " << (int)novos.size() << ")\n";

        if (nels_after > nels_before) {
                InitializeMemory(cmesh, coes, atrito);
                return true;
        } else {
                std::cout << "[PreRefine] sem novos refinamentos; fim.\n";
                return false;
        }
}

// Drop-in no lugar do BisectionFS: híbrido Secant + Bisection com salvaguardas
REAL HybridBracketedFS(TPZCompMesh* cmesh, REAL coes, REAL atrito,
                       REAL lo, REAL hi, REAL tol_fs_rel, int max_it,bool initmem)
{
        std::cout << "\n[bisect] ===== Início da bisseção de FS =====\n"
        << "[bisect] alvo: tol_rel=" << tol_fs_rel
        << "  max_bis=" << max_it
        << "  lo_in=" << lo << "  hi_in=" << hi << "\n";
        auto rel_gap = [](REAL a, REAL b){
                const REAL m = (REAL)0.5*(a+b);
                return (b-a)/std::max<REAL>(m,(REAL)1e-12);
        };
        if (hi < lo) std::swap(lo,hi);

        REAL flo = +1.; // +1 := converge
        REAL fhi = -1.; // -1 := falha
        // loop híbrido
        for (int k=0; k<max_it && rel_gap(lo,hi) > tol_fs_rel; ++k) {

                // tentativa secante regulada (Illinois)
                REAL fs = lo - (flo*(hi-lo))/(fhi - flo + 1.e-16);
                // salvaguarda: mantém dentro do bracket; se ruim, cai na bisseção
                if (!(fs>lo && fs<hi)) fs = 0.5*(lo+hi);

                int itmid=0;
                bool ok = RunAndAccept(cmesh, coes, atrito, fs, itmid,initmem);
                if (ok) {
                        lo = fs;
                        flo = +1.;
                        // Illinois: amortecer o lado que não muda de sinal
                        fhi *= 0.5;
                } else {
                        hi = fs;
                        fhi = -1.;
                        flo *= 0.5;
                }
                std::cout << "[bisect][it " << k << "] hi=" << hi <<" lo=" << lo<< "  -> " << (ok ? "OK" : "FAIL")<<endl;
        }
        std::cout << "[bisect] ===== FIM =====  FS*≈" << lo <<endl;
        return lo; // melhor piso convergente
}

REAL BisectionFS(TPZCompMesh* cmesh, REAL coes, REAL atrito,
                 REAL lo, REAL hi,
                 REAL tol_fs_rel, int max_bis,bool initmem)
{
        auto rel_gap = [](REAL a, REAL b){
                const REAL m = (REAL)0.5*(a+b);
                return (b - a) / std::max<REAL>(m, (REAL)1e-12);
        };

        if (hi < lo) std::swap(lo, hi);

        std::cout << "\n[bisect] ===== Início da bisseção de FS =====\n"
        << "[bisect] alvo: tol_rel=" << tol_fs_rel
        << "  max_bis=" << max_bis
        << "  lo_in=" << lo << "  hi_in=" << hi << "\n";


        // --- bisseção ---
        int k = 0;
        while (k < max_bis) {
                const REAL gap = rel_gap(lo, hi);
                if (gap <= tol_fs_rel) {
                        std::cout << "[bisect][stop] gap_rel=" << gap
                        << " <= tol_rel=" << tol_fs_rel
                        << "  it=" << k << "\n";
                        break;
                }

                const REAL mid = (REAL)0.5*(lo + hi);
                int it_mid = 0;
                const bool ok = RunAndAccept(cmesh, coes, atrito, mid, it_mid,initmem);

                std::cout << "[bisect][it " << k << "] mid=" << mid
                << "  gap_rel=" << gap
                << "  -> " << (ok ? "OK" : "FAIL")
                << " (iters=" << it_mid << ")  ";

                if (ok) {
                        lo = mid;
                        std::cout << "novo lo=" << lo << "\n";
                } else {
                        hi = mid;
                        std::cout << "novo hi=" << hi << "\n";
                }
                k++;
        }

        const REAL fs_star = lo;
        std::cout << "[bisect] ===== FIM =====  FS*≈" << fs_star
        << "  gap_rel_final=" << rel_gap(lo,hi)
        << "  it_bis=" << k << "\n";

        return fs_star; // melhor estimativa do FS* nesta malha
}


plasticmat* CreateMaterial(REAL young, REAL poisson, REAL coes, REAL atrito,
                           TPZManVector<REAL,3> bodyforce, int planestrain, int matid)
{
        TPZElasticResponse ER; ER.SetEngineeringData(young, poisson);

        TPlasticMC mc;
        mc.fYC.SetUp(atrito, atrito, coes, ER);
        mc.fER = ER;
        mc.SetStrengthReductionFactor(1.0);

        auto* material = new plasticmat(matid, planestrain);
        material->SetPlasticityModel(mc);
        material->SetId(matid);
        material->SetWhichLoadVector(0);
        material->SetLoadFactor(1.0);
        material->SetBodyForce(bodyforce);
        return material;
}

TPZCompMesh* CreateCMesh(TPZGeoMesh* gmesh, int pOrder, plasticmat* mat)
{
        TPZCompMesh* cmesh = new TPZCompMesh(gmesh);
        cmesh->SetDefaultOrder(pOrder);
        cmesh->SetDimModel(2);
        cmesh->InsertMaterialObject(mat);

        TPZFMatrix<STATE> val1(2,2,0.0);
        TPZManVector<STATE,2> val2(2,0.0);

        int dir = 3;
        val2[0]=1; val2[1]=1; auto* bc0 = mat->CreateBC(mat, -1, dir, val1, val2);
        val2[0]=1; val2[1]=0; auto* bc1 = mat->CreateBC(mat, -2, dir, val1, val2);
        val2[0]=1; val2[1]=0; auto* bc2 = mat->CreateBC(mat, -5, dir, val1, val2);

        cmesh->InsertMaterialObject(bc0);
        cmesh->InsertMaterialObject(bc1);
        cmesh->InsertMaterialObject(bc2);

        cmesh->SetAllCreateFunctionsContinuousWithMem();
        cmesh->AutoBuild();
        cmesh->AdjustBoundaryElements();
        cmesh->CleanUpUnconnectedNodes();
        return cmesh;
}

void InitializeMemory(TPZCompMesh* cmesh, REAL coesion, REAL atrito)
{
        auto* pMatWithMem2 = dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem>*>(cmesh->MaterialVec()[1]);
        if (!pMatWithMem2) { DebugStop(); }
        pMatWithMem2->SetUpdateMem(true);

        int nels = cmesh->NElements();
        for (int iel=0; iel<nels; iel++) {
                TPZCompEl* cel = cmesh->ElementVec()[iel];
                auto* intel = dynamic_cast<TPZInterpolationSpace*>(cel);
                if (!cel || !intel) continue;
                if (dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem>*>(intel->Material()) != pMatWithMem2
                        || intel->Material()->Id() < 1) continue;

                const TPZIntPoints& intpoints = intel->GetIntegrationRule();
                int nint = intpoints.NPoints();
                TPZManVector<REAL,3> q(2,0.0);

                TPZMaterialDataT<REAL> data;
                intel->InitMaterialData(data);
                data.fNeedsSol = true;

                for (int ip=0; ip<nint; ip++) {
                        REAL w; intpoints.Point(ip, q, w);
                        data.intLocPtIndex = ip;
                        intel->ComputeRequiredData(data, q);

                        int idx = data.intGlobPtIndex;
                        TPZElastoPlasticMem &mem = pMatWithMem2->MemItem(idx);
                        mem.m_elastoplastic_state.fmatprop.Resize(3);
                        mem.m_elastoplastic_state.fmatpropinit.Resize(3);

                        mem.m_elastoplastic_state.fmatpropinit[0] = coesion;
                        mem.m_elastoplastic_state.fmatpropinit[1] = atrito;
                        mem.m_elastoplastic_state.fmatpropinit[2] = atrito;

                        mem.m_elastoplastic_state.fmatprop[0] = coesion;
                        mem.m_elastoplastic_state.fmatprop[1] = atrito;
                        mem.m_elastoplastic_state.fmatprop[2] = atrito;
                }
        }
        pMatWithMem2->SetUpdateMem(false);
}

REAL UyAtNode(TPZCompMesh* cmesh, REAL x, REAL y)
{
        cmesh->LoadReferences();
        auto* gmesh = cmesh->Reference();
        if (!gmesh) return 0.0;

        TPZManVector<REAL,3> X(3,0.0); X[0]=x; X[1]=y;
        TPZManVector<REAL,3> qsi(3,0.0);
        int64_t elindex = 0;

        TPZGeoEl* gel = gmesh->FindElement(X, qsi, elindex,gmesh->Dimension());
        if (!gel || !gel->Reference()) return 0.0;

        auto* cel = gel->Reference();
        auto* intel = dynamic_cast<TPZInterpolationSpace*>(cel);
        if (!intel) return 0.0;

        int nn = gel->NNodes();
        int local = -1;
        for (int i=0;i<nn;i++){
                TPZManVector<REAL,3> co(3,0.0);
                gel->NodePtr(i)->GetCoordinates(co);
                if (std::fabs(co[0]-x) < 1e-10 && std::fabs(co[1]-y) < 1e-10){ local = i; break; }
        }
        if (local < 0) return 0.0;

        int ic = intel->ConnectIndex(local);
        if (ic < 0) return 0.0;

        TPZConnect &c = cmesh->ConnectVec()[ic];
        int64_t seq = c.SequenceNumber();
        TPZBlock &block = cmesh->Block();
        int pos = block.Position(seq);
        int ndof = block.Size(seq);

        TPZFMatrix<REAL> sol = cmesh->Solution();
        if (pos+1 >= sol.Rows() || ndof < 2) return 0.0;

        return sol(pos+1,0);
}

void ComputeElementDeformation(TPZCompMesh* cmesh, TPZVec<REAL>& fPlasticDeformSqJ2)
{
    const int64_t nelem = cmesh->NElements();
    fPlasticDeformSqJ2.resize(nelem);
    fPlasticDeformSqJ2.Fill(0.0);

    // 1 coluna para armazenar o indicador por elemento
    cmesh->ElementSolution().Redim(nelem, 1);

    for (int64_t el = 0; el < nelem; ++el) {
        TPZCompEl* cel = cmesh->ElementVec()[el];
        if (!cel) continue;
        // ignore contorno/interfaces
        if (cel->Dimension() != cmesh->Dimension()) continue;

        // material com memória
        auto* matmem = dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem>*>(cel->Material());
        if (!matmem) { fPlasticDeformSqJ2[el] = 0.0; continue; }

        TPZManVector<int64_t> memindices;
        cel->GetMemoryIndices(memindices);
        if (memindices.size()==0) { fPlasticDeformSqJ2[el] = 0.0; continue; }

        REAL vmax = 0.0;
        int  npts = 0;

        for (int64_t midx : memindices) {
            if (midx < 0) continue;
            // (opcional) checagem extra:
            // if (midx >= matmem->GetMemory().size()) continue;

            const auto& mem = matmem->MemItem(midx);
            const TPZTensor<REAL> epsp = mem.m_elastoplastic_state.EpsP();

            REAL J2 = epsp.J2();
            if (J2 < (REAL)0) J2 = 0;

            const REAL eqp = std::sqrt(J2); // medida equivalente (ajuste se quiser usar sqrt(2/3)*||dev||)
            vmax = std::max(vmax, eqp);
            ++npts;
        }

        // indicador escolhido: MÁXIMO nos IPs do elemento
        // (se quiser MÉDIA, acumule e divida por npts)
        fPlasticDeformSqJ2[el] = (npts > 0) ? vmax : (REAL)0.0;
    }

    // publica na coluna 0 do ElementSolution
    cmesh->SetElementSolution(0, fPlasticDeformSqJ2);
}


void DivideElementsAbove(TPZCompMesh* cmesh, REAL refineaboveval, std::set<int64_t>& out_newels)
{
        cmesh->LoadReferences();
        std::vector<int64_t> to_divide;
        const TPZFMatrix<STATE>& elsol = cmesh->ElementSolution();

        const int64_t ne0 = cmesh->NElements();
        for (int64_t el=0; el<ne0; el++) {
                TPZCompEl* cel = cmesh->ElementVec()[el];
                if (!cel) continue;
                auto* intel = dynamic_cast<TPZInterpolationSpace*>(cel);
                if (!intel) continue;
                if (!dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem>*>(cel->Material())) continue;

                if (elsol.Rows() > el && elsol.Cols() > 0) {
                        if (elsol(el,0) > refineaboveval) to_divide.push_back(el);
                }
        }

        for (auto el : to_divide) {
                TPZCompEl* cel = cmesh->ElementVec()[el]; if (!cel) continue;
                auto* intel = dynamic_cast<TPZInterpolationSpace*>(cel); if (!intel) continue;
                int p = intel->GetPreferredOrder();
                TPZStack<int64_t> sub; const int64_t idx = cel->Index();
                intel->Divide(idx, sub, /*create_boundary_elements=*/0);
                for (int i=0;i<sub.size();i++){
                        out_newels.insert(sub[i]);
                        auto* sc = cmesh->ElementVec()[sub[i]];
                        if (auto* si = dynamic_cast<TPZInterpolationSpace*>(sc)) si->SetPreferredOrder(p);
                }
        }

        // balance
        bool changed = true;
        while (changed){
                changed = false;
                std::set<int64_t> need;
                const int64_t ne = cmesh->NElements();
                for (int64_t el=0; el<ne; el++){
                        TPZCompEl* cel = cmesh->ElementVec()[el]; if (!cel) continue;
                        auto* intel = dynamic_cast<TPZInterpolationSpace*>(cel); if (!intel) continue;
                        TPZGeoEl* gel = cel->Reference(); if (!gel) continue;

                        const int ns = gel->NSides();
                        for (int s=0;s<ns;s++){
                                TPZGeoElSide gs(gel,s);
                                if (gs.Dimension() != gel->Dimension()-1) continue;
                                TPZCompElSide big = gs.LowerLevelCompElementList2(/*onlyintersect=*/1);
                                if (!big) continue;
                                TPZGeoElSide gbig(big.Reference());
                                if (gbig.Element()->Dimension() != gel->Dimension()) continue;
                                if (gel->Level() - gbig.Element()->Level() > 1){
                                        need.insert(big.Element()->Index());
                                }
                        }
                }
                for (auto el : need){
                        TPZCompEl* cel = cmesh->ElementVec()[el]; if (!cel) continue;
                        auto* intel = dynamic_cast<TPZInterpolationSpace*>(cel); if (!intel) continue;
                        if (!dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem>*>(cel->Material())) continue;
                        int p = intel->GetPreferredOrder();
                        TPZStack<int64_t> sub; const int64_t idx = cel->Index();
                        intel->Divide(idx, sub, /*create_boundary_elements=*/0);
                        for (int i=0;i<sub.size();i++){
                                out_newels.insert(sub[i]);
                                auto* sc = cmesh->ElementVec()[sub[i]];
                                if (auto* si = dynamic_cast<TPZInterpolationSpace*>(sc)) si->SetPreferredOrder(p);
                        }
                        changed = true;
                }
        }

        cmesh->AdjustBoundaryElements();
        cmesh->CleanUpUnconnectedNodes();
        cmesh->InitializeBlock();
        cmesh->ExpandSolution();

        // reset memória em todos os materiais com memória
        for (auto& kv : cmesh->MaterialVec()){
                TPZMaterial* m = kv.second;
                if (auto* mm = dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem>*>(m)){
                        mm->ResetMemory();
                }
        }
}

void PRefineElementsAbove(TPZCompMesh* cmesh,
                          REAL refineaboveval,
                          std::set<int64_t>& out_newels,
                          int porder)
{
        if (!cmesh) return;

        // 1) Garantir referências geométricas
        cmesh->LoadReferences();

        const TPZFMatrix<STATE>& elsol = cmesh->ElementSolution();

        // 2) Sanidade: checar se ElementSolution existe
        if (elsol.Rows() == 0 || elsol.Cols() == 0) {
                std::cerr << "[PRefineElementsAbove] ElementSolution vazio: "
                << "chame quem preenche (ex.: Analysis::ComputeElementSolution) antes.\n";
                return;
        }

        const int64_t nelem = cmesh->NElements();

        for (int64_t pos = 0; pos < nelem; ++pos) {
                TPZCompEl* cel = cmesh->ElementVec()[pos];
                if (!cel) continue;

                // só elementos com espaço de interpolação (descarta especiais / multiphysics sem interp)
                auto* intel = dynamic_cast<TPZInterpolationSpace*>(cel);
                if (!intel) continue;

                // filtra materiais: só continua se for um mat com memória elastoplástica
                auto* pMatWithMem =
                dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem>*>(cel->Material());
                if (!pMatWithMem) continue;

                // 3) usar o índice REAL do elemento (pode não coincidir com 'pos' se há buracos)
                const int64_t idx = cel->Index();
                if (idx < 0 || idx >= elsol.Rows()) continue;

                // 4) pega o escalar da 1ª coluna (ajuste se sua métrica estiver em outra coluna)
                const STATE val = elsol(idx, 0);
                if (val < refineaboveval) continue;

                //cout << "porder"<<porder <<endl;
                // 5) aplica p-refine: define a ordem preferida para o elemento
                intel->SetPreferredOrder(porder);

                out_newels.insert(idx);
        }

        // 6) Rebuild de conects/estrutura e vetor solução
        cmesh->AdjustBoundaryElements();
        cmesh->CleanUpUnconnectedNodes(); // remove connects órfãos
        cmesh->ExpandSolution();          // expande Solution conforme novas ordens
        cmesh->ComputeNodElCon();         // (opcional) recomputa incidências
        cmesh->InitializeBlock();         // re-inicializa blocagem conforme connects
}

// ------------------------------------------------------------
// Pós-processo
// ------------------------------------------------------------
void PostProcessVariables(TPZStack<std::string>& scal, TPZStack<std::string>& vec)
{
        scal.Push ( "POrder" );
        scal.Push ( "Atrito" );
        scal.Push ( "Coesion" );
        scal.Push ( "StrainPlasticJ2" );
        scal.Push ( "VolHardening" );
        vec.Push ( "Displacement" );
        vec.Push ( "ShearPlasticDeformation" );
        vec.Push ( "PlasticDeformation" );

}

void CreatePostProcessingMesh(TPZCompMesh* cmesh,TPZPostProcAnalysis* pproc,int matid)
{
        if (pproc->ReferenceCompMesh() != cmesh) {
                pproc->SetCompMesh(cmesh);
                TPZStack<std::string> scal, vec, all;
                PostProcessVariables(scal, vec);
                for (auto i=0; i<scal.size();  ++i) all.Push(scal[i]);
                for (auto i=0; i<vec.size();   ++i) all.Push(vec[i]);
                TPZVec<int> matids(1); matids[0] = matid;
                pproc->SetPostProcessVariables(matids, all);
                TPZFStructMatrix<REAL> str(pproc->Mesh());
                str.SetNumThreads(0);
                pproc->SetStructuralMatrix(str);
        }
        pproc->TransferSolution();
}

void PostElastoplastic(TPZCompMesh* cmesh,const std::string& vtkfile,int matid)
{
        TPZPostProcAnalysis pproc;
        CreatePostProcessingMesh(cmesh, &pproc, matid);
        TPZStack<std::string> scal, vec;
        PostProcessVariables(scal, vec);
        pproc.DefineGraphMesh(/*dim=*/2, scal, vec, vtkfile);
        pproc.PostProcess(2);
}

TPZGeoMesh* TriGMesh(int ref)
{
        TPZGeoMesh* gmesh = new TPZGeoMesh();
        gmesh->SetDimension(2);

        std::vector<std::vector<double>> co = {
                /*0*/{0,0},/*1*/{10,0},/*2*/{20,0},/*3*/{30,0},/*4*/{40,0},/*5*/{50,0},/*6*/{60,0},/*7*/{70,0},
                /*8*/{0,10},/*9*/{10,10},/*10*/{20,10},/*11*/{30,10},/*12*/{40,10},/*13*/{50,10},/*14*/{60,10},/*15*/{70,10},
                /*16*/{0,20},/*17*/{10,20},/*18*/{20,20},/*19*/{30,20},/*20*/{40,20},/*21*/{50,20},/*22*/{60,20},/*23*/{70,20},
                /*24*/{0,30},/*25*/{10,30},/*26*/{20,30},/*27*/{30,30},/*28*/{40,30},/*29*/{50,30},/*30*/{60,30},/*31*/{70,30},
                /*32*/{0,40},/*33*/{10,40},/*34*/{20,40},/*35*/{30,40}
        };

        std::vector<std::vector<int>> topol = {
                /*triangles*/
                {0,1,8},{1,9,8},{1,2,9},{2,10,9},{2,3,10},{3,11,10},{3,4,11},
                {4,12,11},{4,5,12},{5,13,12},{5,6,13},{6,14,13},{6,7,14},{7,15,14},
                {8,9,16},{9,17,16},{9,10,17},{10,18,17},{10,11,18},{11,19,18},{11,12,19},
                {12,20,19},{12,13,20},{13,21,20},{13,14,21},{14,22,21},{14,15,22},{15,23,22},
                {16,17,24},{17,25,24},{17,18,25},{18,26,25},{18,19,26},{19,27,26},{19,20,27},
                {20,28,27},{20,21,28},{21,29,28},{21,22,29},{22,30,29},{22,23,30},{23,31,30},
                {24,25,32},{25,33,32},{25,26,33},{26,34,33},{26,27,34},{27,35,34},{27,28,35},
                /*lines (BCs):*/
                {0,1},{1,2},{2,3},{3,4},{4,5},{5,6},{6,7},       // -1 bottom
                {7,15},{15,23},{23,31},                          // -2 right
                {31,30},{30,29},{29,28},                         // -3 top right
                {35,34},{34,33},{33,32},                         // -4 top left
                {32,24},{24,16},{16,8},{8,0},                    // -5 left
                {28,35}                                          // -6 ramp
        };

        gmesh->NodeVec().Resize(co.size());
        TPZVec<REAL> coord(2);
        for (int i = 0; i < (int)co.size(); i++) {
                coord[0] = co[i][0];
                coord[1] = co[i][1];
                gmesh->NodeVec()[i] = TPZGeoNode(i, coord, *gmesh);
        }

        TPZVec<long> topotri(3), topoline(2);
        for (int i = 0; i < (int)topol.size(); i++) {
                if (topol[i].size() == 3) {
                        topotri[0] = topol[i][0]; topotri[1] = topol[i][1]; topotri[2] = topol[i][2];
                        new TPZGeoElRefPattern<pzgeom::TPZGeoTriangle>(i, topotri, /*matid=*/1, *gmesh);
                } else {
                        topoline[0] = topol[i][0]; topoline[1] = topol[i][1];
                        REAL x0 = co[topoline[0]][0], y0 = co[topoline[0]][1];
                        REAL xf = co[topoline[1]][0], yf = co[topoline[1]][1];
                        REAL tol = 1.e-3, L = 70, h1 = 30, h2 = 10;

                        int bcid = 0;
                        if (std::fabs(y0-0) < tol && std::fabs(yf-0) < tol)           bcid = -1; // bottom
                        else if (std::fabs(x0-L)<tol && std::fabs(xf-L)<tol)          bcid = -2; // right
                        else if (std::fabs(y0-h1)<tol && std::fabs(yf-h1)<tol)        bcid = -3; // top right
                        else if (std::fabs(y0-(h1+h2))<tol && std::fabs(yf-(h1+h2))<tol) bcid = -4; // top left
                        else if (std::fabs(x0-0)<tol && std::fabs(xf-0)<tol)          bcid = -5; // left
                        else if (std::fabs(xf-x0)>tol && std::fabs(yf-y0)>tol)        bcid = -6; // ramp
                        else {
                                std::cout << "bc element not found.\n"; DebugStop();
                        }
                        new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(i, topoline, bcid, *gmesh);
                }
        }

        gmesh->BuildConnectivity();
        for (int d = 0; d < ref; d++) {
                int nel = gmesh->NElements();
                TPZManVector<TPZGeoEl*> sub;
                for (int iel = 0; iel < nel; iel++) {
                        gmesh->ElementVec()[iel]->Divide(sub);
                }
        }

        std::ofstream vtk("gmeshtri.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, vtk, true);
        return gmesh;
}

TPZFMatrix<REAL> ComputeTheta (int M, int samples, uint32_t seed)
{
        std::mt19937 gen(seed);
        std::normal_distribution<double> N01(0.,1.);
        TPZFMatrix<REAL>  THETA (M, samples, 0.);
        for (int j=0;j<samples;j++)
                for (int i=0;i<M;i++) THETA(i,j) = (REAL)N01(gen);

        return THETA;
}

// KL paramétrico em Lx, Ly
TPZCompMesh * BuildCompMeshKL_Param(TPZGeoMesh *gmesh, int porder, int matId,
                                           REAL Lx, REAL Ly)
{
        auto *cmesh = new TPZCompMesh(gmesh);
        cmesh->SetDimModel(2);
        cmesh->SetDefaultOrder(porder);
        cmesh->SetAllCreateFunctionsContinuous();

        auto *mat = new TPZMatKLKernel(matId, 2, Lx, Ly);
        mat->SetId(matId);
        cmesh->InsertMaterialObject(mat);
        cmesh->AutoBuild();
        cmesh->AdjustBoundaryElements();
        cmesh->CleanUpUnconnectedNodes();
        return cmesh;
}
// ------------------------------------------------------------
// PHI_sqrtLambda (ndof × nm). Mantido por compatibilidade.
// ------------------------------------------------------------
TPZFMatrix<STATE> BuildPhiSqrtLambda(
        TPZCompMesh* cmesh,
        const TPZFMatrix<CSTATE>& eigenvectors, // cols = modos
        const TPZVec<CSTATE>&     eigenvalues,  // tamanho >= nm
        int M)
{
        TPZFMatrix<STATE>soltemp=cmesh->Solution();
        const int ndof = cmesh->NEquations();
        const int nm   = std::min(M, (int)eigenvectors.Cols());
        //if (nm <= 0) return;

        TPZFMatrix<STATE> PHI_sqrtLambda(ndof, nm, (STATE)0);

        // materiais de volume (sem BC) para integração
        std::set<int> mats;
        for (auto &it : cmesh->MaterialVec()) {
                if (!it.second) continue;
                if (dynamic_cast<TPZBndCond*>(it.second)) continue;
                mats.insert(it.first);
        }

        TPZFMatrix<STATE> sol(ndof,1,(STATE)0);

        for (int k = 0; k < nm; ++k) {
                // 1) autovetor k (parte real)
                for (int i=0; i<ndof; i++) {
                        sol(i,0) = (STATE)eigenvectors.GetVal(i,k).real();
                }

                // 2) normalização na métrica de massa: ||phi||_M = sqrt( ∫ phi^2 dΩ )
                double scale = 1.0;

                cmesh->LoadSolution(sol);

                TPZVec<STATE> Ivec = cmesh->Integrate("SolutionSquared", mats);
                const double I = (Ivec.size() ? (double)Ivec[0] : 0.0);
                if (I > 1e-30) scale = 1.0/std::sqrt(I);

                // 3) √λ_k (parte real, truncada em 0)
                const double lamk    = std::max(0.0, (double)eigenvalues[k].real());
                const double rootlam = std::sqrt(lamk);

                // 4) coluna k := √λ_k * (phi_k * scale)
                for (int i=0; i<ndof; i++) {
                        PHI_sqrtLambda(i,k) = (STATE)(rootlam* (double)sol(i,0) * scale);
                }
        }

        cmesh->Solution().Zero();
        return PHI_sqrtLambda;
}
void ComputeFieldMulti(const std::vector<TPZCompMesh*>& sources,TPZCompMesh* target,
                              const std::vector<FieldSpec>& specs,int idxE, int idxNu)
{
        if (!target) DebugStop();
        if (sources.empty()) DebugStop();
        if ((int)specs.size() < (int)sources.size()) DebugStop();

        auto *targetmatwithmem =dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem> *> ( target->MaterialVec() [1] );
        if(!targetmatwithmem){ std::cout<<"material com memoria nao inicializado\n"; DebugStop(); }

        // Pré-cálculos por campo
        struct LN { REAL mu_ln=0., sigma_ln=0.; };
        std::vector<LN> ln(sources.size());
        std::vector<REAL> sigma_g(sources.size(), 0.0); // desvio para Gaussiano (sigma_x = cov_x*mu_x)

        for (size_t i=0;i<sources.size();++i){
                const auto &sp = specs[i];
                if (sp.model == FieldModel::Lognormal){
                        const REAL cov2 = sp.cov_x * sp.cov_x;
                        const REAL sigma_ln = std::sqrt(std::log(1.0 + cov2));      // σ_lnX
                        const REAL mu_ln    = std::log(std::max(sp.mu_x, (REAL)1e-14)) - 0.5*sigma_ln*sigma_ln; // μ_lnX
                        ln[i] = { mu_ln, sigma_ln };
                } else {
                        sigma_g[i] = sp.cov_x * sp.mu_x; // σ_x = COV_x * μ_x
                }
        }

        for (auto *sm : sources) { if(!sm) DebugStop(); sm->LoadReferences(); }
        target->LoadReferences();

        targetmatwithmem->SetUpdateMem(true);
        TPZAdmChunkVector<TPZElastoPlasticMem> &mem = *targetmatwithmem->GetMemory();

        const int nels = target->NElements();
        for (int iel=0; iel<nels; iel++)
        {

                auto *targetcel = target->ElementVec()[iel];

                if(!targetcel) continue;
                auto *targetintel = dynamic_cast<TPZInterpolationSpace*>(targetcel);
                if(!targetintel) continue;

                auto *intelmat = dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem>*>(targetintel->Material());
                if (intelmat != targetmatwithmem) continue;

                TPZMaterialDataT<STATE> datatarget;
                TPZIntPoints& intrule = targetintel->GetIntegrationRule();
                const int nip = intrule.NPoints();

                for (int ip=0; ip<nip; ip++)
                {
                        REAL w; (void)w;
                        TPZManVector<REAL,3> intpt(3,0.);
                        intrule.Point(ip, intpt, w);

                        datatarget.intLocPtIndex = ip;
                        datatarget.fNeedsSol = false;
                        targetintel->InitMaterialData(datatarget);
                        targetintel->ComputeRequiredData(datatarget, intpt);

                        // Avalia z_i(x) em cada fonte e transforma para X_i(x) com mu_x/cov_x
                        std::vector<REAL> X(sources.size(), 0.0);

                        for (size_t f=0; f<sources.size(); ++f)
                        {
                                TPZManVector<REAL,3> qsi(3,0.);
                                int64_t elidsrc;
                                TPZGeoEl *gelsrc = sources[f]->Reference()->FindElement(datatarget.x, qsi, elidsrc, sources[f]->Dimension());

                                if (!gelsrc || !gelsrc->Reference()) DebugStop();
                                auto *celsource = dynamic_cast<TPZInterpolationSpace*>(gelsrc->Reference());
                                if (!celsource) DebugStop();

                                TPZMaterialDataT<STATE> dsrc;
                                dsrc.intLocPtIndex = ip;
                                dsrc.fNeedsSol = true;
                                celsource->InitMaterialData(dsrc);
                                celsource->ComputeRequiredData(dsrc, qsi);

                                const REAL z = (REAL)dsrc.sol[0][0]; // realização padrão (campo KL normalizado)

                                if (specs[f].model == FieldModel::Lognormal) {
                                        X[f] = std::exp( ln[f].mu_ln + ln[f].sigma_ln * z );
                                } else {
                                        X[f] = specs[f].mu_x + sigma_g[f] * z; // Gaussiano com μₓ e COVₓ
                                }
                        }

                        // grava em memória: E <- X[idxE], nu <- X[idxNu]
                        const int indextarget = datatarget.intGlobPtIndex;
                        TPZElasticResponse ER;
                        const REAL coes  = (idxE  >=0 && idxE  < (int)X.size()) ? X[idxE]  : 1.0;
                        const REAL atrito = (idxNu >=0 && idxNu < (int)X.size()) ? X[idxNu] : 0.0;

                        // cout << "coes = "<< coes <<endl;
                        // cout << "atrito = "<< atrito <<endl;

                        mem[indextarget].m_elastoplastic_state.fmatprop.Resize ( 3 );
                        if ( !mem[indextarget].m_elastoplastic_state.fmatprop.size() ) {
                                cout << "deve-se inicializar corretamente o matprop"<<endl;
                                DebugStop();
                        }
                        mem[indextarget].m_elastoplastic_state.fmatprop.Resize ( 3 );
                        mem[indextarget].m_elastoplastic_state.fmatprop[0]=coes;
                        mem[indextarget].m_elastoplastic_state.fmatprop[1]=atrito;

                }

        }

        targetmatwithmem->SetUpdateMem(false);
}
