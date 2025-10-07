// // mc_sweep.cpp
// #include <random>
// #include <cmath>
// #include <string>
// #include <iostream>
// #include <fstream>
// #include <iomanip>
// #include <sstream>
// #include <set>
// #include <vector>
// #include <utility>
//
// // ---- NeoPZ ----
// #include "pzgmesh.h"
// #include "pzcmesh.h"
// #include "pzgeoquad.h"
// #include "TPZGeoLinear.h"
// #include "tpzgeoelrefpattern.h"
// #include "pzintel.h"
// #include "pzstack.h"
//
// #include "TPZEigenAnalysis.h"
// #include "TPZKrylovEigenSolver.h"
// #include "pzdoublestrmatriz.h"
// #include "pzskylstrmatrix.h"
//
// #include "TPZMatKLKernel.h"
// #include "Elasticity/TPZMatElastic2DMem.h"
// #include "Elasticity/TPZElasticMem.h"
//
// #include "pzpostprocanalysis.h"
// #include "pzstepsolver.h"
// #include "TPZBFileStream.h"
//
// // ------------------------------------------------------------------
// // Helpers (iguais/compatíveis ao seu setup)
// // ------------------------------------------------------------------
// static TPZGeoMesh * CreateGeoMeshMathematicaLike(int matId, int nRef)
// {
//     auto * gmesh = new TPZGeoMesh;
//     gmesh->SetDimension(2);
//
//     const double X[4][2] = { {-0.5,-0.5},{0.5,-0.5},{0.5,0.5},{-0.5,0.5} };
//     gmesh->NodeVec().Resize(4);
//     for (int i=0;i<4;i++){
//         TPZManVector<REAL,3> c(3,0.); c[0]=X[i][0]; c[1]=X[i][1];
//         gmesh->NodeVec()[i].Initialize(c, *gmesh);
//     }
//
//     TPZManVector<int64_t,4> nodes(4);
//     nodes[0]=0; nodes[1]=1; nodes[2]=2; nodes[3]=3;
//     int64_t idx=0;
//     gmesh->CreateGeoElement(EQuadrilateral, nodes, matId, idx);
//
//     idx=1; { TPZVec<long> L(2); L[0]=0; L[1]=1; new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(idx, L, -1, *gmesh); }
//     idx=2; { TPZVec<long> L(2); L[0]=2; L[1]=3; new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(idx, L, -2, *gmesh); }
//
//     gmesh->BuildConnectivity();
//
//     for (int r=0;r<nRef;r++){
//         const int nel = gmesh->NElements();
//         TPZManVector<TPZGeoEl*> sub;
//         for (int i=0;i<nel;i++) if (auto *eg = gmesh->ElementVec()[i]) eg->Divide(sub);
//     }
//     return gmesh;
// }
//
// static TPZCompMesh * CompMeshElastic(TPZGeoMesh * gmesh, int porder, int matId)
// {
//     auto *cmesh = new TPZCompMesh(gmesh);
//     cmesh->SetDimModel(2);
//     cmesh->SetDefaultOrder(porder);
//     cmesh->SetAllCreateFunctionsContinuousWithMem();
//
//     auto *mat = new TPZMatElastic2DMem<TPZElasticMem>(matId);
//     mat->SetId(matId);
//     TPZElasticResponse ER; ER.SetEngineeringData(/*E=*/1.0, /*nu=*/0.0);
//     mat->SetUpdateMem(true);
//     mat->SetElasticResponse(ER);
//     mat->SetUpdateMem(false);
//     cmesh->InsertMaterialObject(mat);
//
//     TPZFMatrix<STATE> val1 (2,2,0.);
//     TPZVec<STATE>     val2 (2,0.);
//
//     val2[0]=1.; val2[1]=1.; auto bcclamp = mat->CreateBC(mat,-1,3,val1,val2);
//     val2[0]=0.; val2[1]=1.; auto bcload  = mat->CreateBC(mat,-2,1,val1,val2);
//
//     cmesh->InsertMaterialObject(bcclamp);
//     cmesh->InsertMaterialObject(bcload);
//     cmesh->AutoBuild();
//     return cmesh;
// }
//
// // KL paramétrico em Lx, Ly
// static TPZCompMesh * BuildCompMeshKL_Param(TPZGeoMesh *gmesh, int porder, int matId,
//                                            REAL Lx, REAL Ly)
// {
//     auto *cmesh = new TPZCompMesh(gmesh);
//     cmesh->SetDimModel(2);
//     cmesh->SetDefaultOrder(porder);
//     cmesh->SetAllCreateFunctionsContinuous();
//
//     auto *mat = new TPZMatKLKernel(matId, 2, Lx, Ly);
//     mat->SetId(matId);
//     cmesh->InsertMaterialObject(mat);
//     cmesh->AutoBuild();
//     cmesh->AdjustBoundaryElements();
//     cmesh->CleanUpUnconnectedNodes();
//     return cmesh;
// }
//
// static TPZFMatrix<REAL> ComputeTheta (int M, int samples, uint32_t seed=12345)
// {
//     std::mt19937 gen(seed);
//     std::normal_distribution<double> N01(0.,1.);
//     TPZFMatrix<REAL>  THETA (M, samples, 0.);
//     for (int j=0;j<samples;j++)
//         for (int i=0;i<M;i++) THETA(i,j) = (REAL)N01(gen);
//         return THETA;
// }
//
// static TPZFMatrix<STATE> BuildPhiSqrtLambda(
//     TPZCompMesh* cmesh,
//     const TPZFMatrix<CSTATE>& evecs, // colunas
//     const TPZVec<CSTATE>&     evals,
//     int M)
// {
//     const int ndof = cmesh->NEquations();
//     const int nm   = std::min(M, (int)evecs.Cols());
//     TPZFMatrix<STATE> PHI(ndof, nm, (STATE)0);
//
//     std::set<int> mats;
//     for (auto &it : cmesh->MaterialVec()) {
//         if (!it.second) continue;
//         if (dynamic_cast<TPZBndCond*>(it.second)) continue;
//         mats.insert(it.first);
//     }
//
//     TPZFMatrix<STATE> sol(ndof,1,(STATE)0);
//     for (int k=0;k<nm;k++){
//         for (int i=0;i<ndof;i++) sol(i,0) = (STATE)evecs.GetVal(i,k).real();
//         cmesh->LoadSolution(sol);
//         TPZVec<STATE> Ivec = cmesh->Integrate("SolutionSquared", mats);
//         const double I = (Ivec.size()? (double)Ivec[0]:0.0);
//         const double scale = (I>1e-30)? 1.0/std::sqrt(I) : 1.0;
//         const double rootlam = std::sqrt(std::max(0.0,(double)evals[k].real()));
//         for (int i=0;i<ndof;i++) PHI(i,k) = (STATE)(rootlam * (double)sol(i,0) * scale);
//     }
//     cmesh->Solution().Zero();
//     return PHI;
// }
//
// // Campo -> memória elástica
// static void ComputeField(TPZVec<TPZCompMesh*> source,TPZCompMesh* target)
// {
//     auto *targetmatwithmem =
//     dynamic_cast<TPZMatWithMem<TPZElasticMem>*>(target->MaterialVec()[1]);
//     if(!targetmatwithmem){ std::cout<<"material com memoria nao inicializado\n"; DebugStop(); }
//     targetmatwithmem->SetUpdateMem(true);
//     TPZAdmChunkVector<TPZElasticMem> &mem = *targetmatwithmem->GetMemory();
//
//     const REAL mu    = 1.0;
//     const REAL sigma = 0.2;
//     const bool use_lognormal = false;
//     const REAL lambda_ln = std::log(mu) - 0.5*std::log(1.0 + (sigma*sigma)/(mu*mu));
//     const REAL xi_ln     = std::sqrt(std::log(1.0 + (sigma*sigma)/(mu*mu)));
//
//     const REAL mu2    = 0.2;
//     const REAL sigma2 = 0.2;
//     const REAL lambda_ln2 = std::log(mu2) - 0.5*std::log(1.0 + (sigma2*sigma2)/(mu2*mu2));
//     const REAL xi_ln2     = std::sqrt(std::log(1.0 + (sigma2*sigma2)/(mu2*mu2)));
//
//     const int nels = target->NElements();
//     for (int iel=0; iel<nels; iel++)
//     {
//         auto *targetcel = target->ElementVec()[iel];
//         if(!targetcel) continue;
//         auto *targetintel = dynamic_cast<TPZInterpolationSpace*>(targetcel);
//         if(!targetintel) continue;
//
//         auto *intelmat = dynamic_cast<TPZMatWithMem<TPZElasticMem>*>(targetintel->Material());
//         if (intelmat != targetmatwithmem) continue;
//
//         TPZMaterialDataT<STATE> datatarget,datasource,datasource2;
//
//         TPZIntPoints& intrule = targetintel->GetIntegrationRule();
//         const int nip = intrule.NPoints();
//
//         for (int ip=0; ip<nip; ip++)
//         {
//             REAL w;
//             TPZManVector<REAL,3> intpt(3,0.);
//             intrule.Point(ip, intpt, w);
//
//             datatarget.intLocPtIndex = ip;
//             datatarget.fNeedsSol = false;
//             targetintel->InitMaterialData(datatarget);
//             targetintel->ComputeRequiredData(datatarget, intpt);
//
//             TPZManVector<REAL,3> qsisource(2,0.),qsisource2(2,0.);
//             int64_t elidsrc;
//             TPZGeoEl *gelsource = source[0]->Reference()->FindElement(datatarget.x, qsisource, elidsrc, source[0]->Dimension());
//             TPZGeoEl *gelsource2 = source[1]->Reference()->FindElement(datatarget.x, qsisource, elidsrc, source[1]->Dimension());
//             if (!gelsource) DebugStop();
//
//             auto *celsource = gelsource->Reference();
//             auto *celsource2 = gelsource2->Reference();
//             if (!celsource) DebugStop();
//             auto *intelsource = dynamic_cast<TPZInterpolationSpace*>(celsource);
//             auto *intelsource2 = dynamic_cast<TPZInterpolationSpace*>(celsource2);
//             if (!intelsource) DebugStop();
//
//             datasource.intLocPtIndex = ip;
//             datasource.fNeedsSol = true;
//             intelsource->InitMaterialData(datasource);
//             intelsource->ComputeRequiredData(datasource, qsisource);
//
//             datasource2.intLocPtIndex = ip;
//             datasource2.fNeedsSol = true;
//             intelsource2->InitMaterialData(datasource2);
//             intelsource2->ComputeRequiredData(datasource2, qsisource2);
//
//             REAL H1,H2;
//             if (!use_lognormal) {
//                 H1 = mu + sigma * datasource.sol[0][0];
//                 H2 = mu2 + sigma2 * datasource.sol[0][0];
//             } else {
//                 H1 = std::exp(lambda_ln + xi_ln * datasource.sol[0][0]);
//                 H2 = std::exp(lambda_ln2 + xi_ln2 * datasource2.sol[0][0]);
//             }
//
//
//             const int indextarget = datatarget.intGlobPtIndex;
//             TPZElasticResponse ER;
//             ER.SetEngineeringData(/*E=*/H1, /*nu=*/H2);
//             mem[indextarget].m_ER = ER;
//         }
//     }
//     targetmatwithmem->SetUpdateMem(false);
// }
//
// static REAL UyAtPoint(TPZCompMesh* cmesh, REAL x, REAL y)
// {
//     TPZManVector<REAL,3> X(3,0.); X[0]=x; X[1]=y;
//     TPZManVector<REAL,3> qsi(3,0.); int64_t elid=0;
//     cmesh->LoadReferences();
//     TPZGeoEl* gel = cmesh->Reference()->FindElement(X,qsi,elid,2);
//     if(!gel || !gel->Reference()) DebugStop();
//     auto* cel = dynamic_cast<TPZInterpolationSpace*>(gel->Reference());
//     if(!cel) DebugStop();
//     int var = cel->Material()->VariableIndex("Displacement");
//     TPZVec<REAL> sol; cel->Solution(qsi, var, sol);
//     return sol[1];
// }
//
// static REAL SolveOnce_Uy(TPZCompMesh* cm, REAL x, REAL y)
// {
//     TPZLinearAnalysis an(cm);
//     TPZSkylineStructMatrix<STATE> str(cm);
//     an.SetStructuralMatrix(str);
//     TPZStepSolver<REAL> direct; direct.SetDirect(ELDLt);
//     an.SetSolver(direct);
//     an.Assemble();
//     an.Solve();
//     return UyAtPoint(cm,x,y);
// }
//
// // ------------------------------------------------------------
// // Sweep Lx,Ly -> gera hhat_<Lx>_<Ly>.bin (buildfields=1)
// // ou lê hhat e salva mc_results_<Lx>_<Ly>.csv (buildfields=0)
// // ------------------------------------------------------------
// int main()
// {
//     const int   kMatId   = 1;
//     const int   porder   = 2;
//     const int   nRef     = 3;
//     const int   M        = 20;
//     int         NsampGen = 30000;   // usado ao gerar hhat
//     const REAL  xprobe   = -0.5;
//     const REAL  yprobe   =  0.5;
//
//     std::vector<std::pair<REAL,REAL>> grid = {
//         {0.25,0.25},{0.25,0.50},{0.25,0.75},{0.25,1.0},
//         {0.50,0.25},{0.50,0.50},{0.50,0.75},{0.50,1.0},
//         {0.75,0.25},{0.75,0.50},{0.75,0.75},{0.75,1.0},
//         {1.00,0.25},{1.00,0.50},{1.00,0.75},{1.00,1.0}
//     };
//
//     const bool buildfields = false; // mude para false para rodar os MC e salvar CSV
//
//     if (buildfields)
//     {
//         for (auto [Lx,Ly] : grid)
//         {
//             TPZGeoMesh*  gmesh   = CreateGeoMeshMathematicaLike(kMatId, nRef);
//             TPZCompMesh* cmeshKL = BuildCompMeshKL_Param(gmesh, porder, kMatId, Lx, Ly);
//
//             TPZEigenAnalysis an(cmeshKL,false);
//             pzdoublestrmatriz<REAL> sm(cmeshKL);
//             sm.SetCAssembly(pzdoublestrmatriz<REAL>::ECAssembly::Galerkin);
//             an.SetStructuralMatrix(sm);
//
//             TPZKrylovEigenSolver<STATE> esolver;
//             esolver.SetAsGeneralised(true);
//             esolver.SetEigenSorting(TPZEigenSort::AbsDescending);
//             esolver.SetNEigenpairs(cmeshKL->NEquations());
//             esolver.SetKrylovDim(cmeshKL->NEquations());
//             esolver.SetTolerance(1e-10);
//             an.SetSolver(esolver);
//
//             an.Assemble();
//             an.Solve();
//             TPZFMatrix<CSTATE> evecs = an.Eigenvectors();
//             TPZVec<CSTATE>     evals = an.Eigenvalues();
//
//             TPZFMatrix<STATE> PHI   = BuildPhiSqrtLambda(cmeshKL, evecs, evals, M);
//             TPZFMatrix<REAL>  THETAE = ComputeTheta(M, NsampGen, /*seed*/12345);
//             TPZFMatrix<REAL>  THETAMU = ComputeTheta(M, NsampGen, /*seed*/12345);
//             TPZFMatrix<REAL>  hhatE,hhatMU; // (ndof x NsampGen)
//             PHI.Multiply(THETAE, hhatE);
//             PHI.Multiply(THETAMU, hhatMU);
//             // nome estável (fixed, 6 casas)
//             std::ostringstream oss; oss.setf(std::ios::fixed);
//             oss << "hhat_" << std::setprecision(6) << Lx << "_" << Ly << ".bin";
//             TPZBFileStream out;
//             out.OpenWrite(oss.str());
//             hhatE.Write(out,0);
//             hhatMU.Write(out,0);
//
//             delete cmeshKL;
//             delete gmesh;
//             std::cout << "Gerado: " << oss.str() << "\n";
//         }
//     }
//     else
//     {
//         for (auto [Lx,Ly] : grid)
//         {
//             // malhas (para mapear campo -> elasticidade)
//             TPZGeoMesh*  gmesh     = CreateGeoMeshMathematicaLike(kMatId, nRef);
//             TPZCompMesh* cmeshdummy1   = BuildCompMeshKL_Param(gmesh, porder, kMatId, 1.0, 1.0);
//             TPZCompMesh* cmeshdummy2   = BuildCompMeshKL_Param(gmesh, porder, kMatId, 1.0, 1.0);
//             TPZCompMesh* celastic  = CompMeshElastic(gmesh, porder, kMatId);
//
//             // lê hhat
//             std::ostringstream fin; fin.setf(std::ios::fixed);
//             fin << "hhat_" << std::setprecision(6) << Lx << "_" << Ly << ".bin";
//             TPZFMatrix<REAL> hhatE,hhatMU;
//             {
//                 TPZBFileStream in;
//                 in.OpenRead(fin.str());
//                 hhatE.Read(in, 0);
//                 hhatMU.Read(in, 0);
//             }
//
//             //hhatE.Print("E");
//            // hhatMU.Print("mu");
//
//             // quantas amostras usar nesta passada (no máx. hhat.Cols())
//             const int NsampRun = std::min(100, (int)hhatE.Cols());
//
//             // abre CSV compatível com postlx_ly.py: mc_results_<Lx>_<Ly>.csv
//             std::ostringstream fout; fout.setf(std::ios::fixed);
//             fout << "mc_results_" << std::setprecision(6) << Lx << "_" << Ly << ".csv";
//             std::ofstream outcsv(fout.str());
//             outcsv.setf(std::ios::fixed); outcsv << std::setprecision(10);
//             outcsv << "sample,uy\n";
//
//             std::cout << "solving mc "<<fin.str()<< " NsampRun = "<< NsampRun<<std::endl;
//             double sum=0.0, sum2=0.0;
//             for (int s=0; s<NsampRun; ++s)
//             {
//                 TPZFMatrix<REAL> colfieldE(hhatE.Rows(),1),colfieldMU(hhatE.Rows(),1);
//                 for (int i=0;i<hhatE.Rows();++i) colfieldE(i,0) = hhatE(i,s);
//                 for (int i=0;i<hhatMU.Rows();++i) colfieldMU(i,0) = hhatMU(i,s);
//
//                 TPZVec<TPZCompMesh*> colfieldvec(2);
//                 cmeshdummy1->LoadSolution(hhatE);
//                 cmeshdummy1->LoadReferences();
//                 colfieldvec[0]=cmeshdummy1;
//                 cmeshdummy2->LoadSolution(hhatMU);
//                 cmeshdummy2->LoadReferences();
//                 colfieldvec[1]=cmeshdummy2;
//
//                 ComputeField(colfieldvec, celastic);
//                 const REAL uy = SolveOnce_Uy(celastic, xprobe, yprobe);
//
//                 outcsv << s << "," << uy << "\n";
//                 sum  += uy;
//                 sum2 += (double)uy*(double)uy;
//             }
//             outcsv.close();
//
//             const double mean = sum/NsampRun;
//             const double var  = std::max(0.0, sum2/NsampRun - mean*mean);
//             const double stdv = std::sqrt(var);
//             std::cout << "[Lx="<<Lx<<", Ly="<<Ly<<"] Nsamp="<<NsampRun
//             << " mean="<<mean<<" std="<<stdv
//             << " -> salvo em " << fout.str() << "\n";
//
//             delete celastic;
//             delete cmeshdummy1;
//             delete cmeshdummy2;
//             delete gmesh;
//         }
//     }
//     return 0;
// }
// mc_sweep.cpp
#include <random>
#include <cmath>
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <set>
#include <vector>
#include <utility>

// ---- NeoPZ ----
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
#include "Plasticity/TPZElasticResponse.h"

#include "pzpostprocanalysis.h"
#include "pzstepsolver.h"
#include "TPZBFileStream.h"

// ------------------------------------------------------------------
// Helpers (iguais/compatíveis ao seu setup)
// ------------------------------------------------------------------
static TPZGeoMesh * CreateGeoMeshMathematicaLike(int matId, int nRef)
{
    auto * gmesh = new TPZGeoMesh;
    gmesh->SetDimension(2);

    const double X[4][2] = { {-0.5,-0.5},{0.5,-0.5},{0.5,0.5},{-0.5,0.5} };
    gmesh->NodeVec().Resize(4);
    for (int i=0;i<4;i++){
        TPZManVector<REAL,3> c(3,0.); c[0]=X[i][0]; c[1]=X[i][1];
        gmesh->NodeVec()[i].Initialize(c, *gmesh);
    }

    TPZManVector<int64_t,4> nodes(4);
    nodes[0]=0; nodes[1]=1; nodes[2]=2; nodes[3]=3;
    int64_t idx=0;
    gmesh->CreateGeoElement(EQuadrilateral, nodes, matId, idx);

    idx=1; { TPZVec<long> L(2); L[0]=0; L[1]=1; new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(idx, L, -1, *gmesh); }
    idx=2; { TPZVec<long> L(2); L[0]=2; L[1]=3; new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(idx, L, -2, *gmesh); }

    gmesh->BuildConnectivity();

    for (int r=0;r<nRef;r++){
        const int nel = gmesh->NElements();
        TPZManVector<TPZGeoEl*> sub;
        for (int i=0;i<nel;i++) if (auto *eg = gmesh->ElementVec()[i]) eg->Divide(sub);
    }
    return gmesh;
}

static TPZCompMesh * CompMeshElastic(TPZGeoMesh * gmesh, int porder, int matId)
{
    auto *cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(2);
    cmesh->SetDefaultOrder(porder);
    cmesh->SetAllCreateFunctionsContinuousWithMem();

    auto *mat = new TPZMatElastic2DMem<TPZElasticMem>(matId);
    mat->SetId(matId);
    TPZElasticResponse ER; ER.SetEngineeringData(/*E=*/1.0, /*nu=*/0.0);
    mat->SetUpdateMem(true);
    mat->SetElasticResponse(ER);
    mat->SetUpdateMem(false);
    cmesh->InsertMaterialObject(mat);

    TPZFMatrix<STATE> val1 (2,2,0.);
    TPZVec<STATE>     val2 (2,0.);

    val2[0]=1.; val2[1]=1.; auto bcclamp = mat->CreateBC(mat,-1,3,val1,val2);
    val2[0]=0.; val2[1]=1.; auto bcload  = mat->CreateBC(mat,-2,1,val1,val2);

    cmesh->InsertMaterialObject(bcclamp);
    cmesh->InsertMaterialObject(bcload);
    cmesh->AutoBuild();
    return cmesh;
}

// KL paramétrico em Lx, Ly
static TPZCompMesh * BuildCompMeshKL_Param(TPZGeoMesh *gmesh, int porder, int matId,
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

static TPZFMatrix<REAL> ComputeTheta (int M, int samples, uint32_t seed=12345)
{
    std::mt19937 gen(seed);
    std::normal_distribution<double> N01(0.,1.);
    TPZFMatrix<REAL>  THETA (M, samples, 0.);
    for (int j=0;j<samples;j++)
        for (int i=0;i<M;i++) THETA(i,j) = (REAL)N01(gen);
        return THETA;
}

static TPZFMatrix<STATE> BuildPhiSqrtLambda(
    TPZCompMesh* cmesh,
    const TPZFMatrix<CSTATE>& evecs, // colunas
    const TPZVec<CSTATE>&     evals,
    int M)
{
    const int ndof = cmesh->NEquations();
    const int nm   = std::min(M, (int)evecs.Cols());
    TPZFMatrix<STATE> PHI(ndof, nm, (STATE)0);

    std::set<int> mats;
    for (auto &it : cmesh->MaterialVec()) {
        if (!it.second) continue;
        if (dynamic_cast<TPZBndCond*>(it.second)) continue;
        mats.insert(it.first);
    }

    TPZFMatrix<STATE> sol(ndof,1,(STATE)0);
    for (int k=0;k<nm;k++){
        for (int i=0;i<ndof;i++) sol(i,0) = (STATE)evecs.GetVal(i,k).real();
        cmesh->LoadSolution(sol);
        TPZVec<STATE> Ivec = cmesh->Integrate("SolutionSquared", mats);
        const double I = (Ivec.size()? (double)Ivec[0]:0.0);
        const double scale = (I>1e-30)? 1.0/std::sqrt(I) : 1.0;
        const double rootlam = std::sqrt(std::max(0.0,(double)evals[k].real()));
        for (int i=0;i<ndof;i++) PHI(i,k) = (STATE)(rootlam * (double)sol(i,0) * scale);
    }
    cmesh->Solution().Zero();
    return PHI;
}

// ------------------------------------------------------------------
// Generalização: ComputeFieldMulti (N campos) -> escreve em memória elástica
// H[i] = (mu[i] + sigma[i]*z)  OU  lognormal: exp(lambda_ln + xi_ln*z)
// Escolhe quais índices entram como E e nu (idxE, idxNu).
// ------------------------------------------------------------------
// ------------------------------------------------------------------
// Especificação por campo i: média mu_x[i], COV_x[i], modelo (gauss/lognormal)
// ------------------------------------------------------------------
enum class FieldModel { Gaussian, Lognormal };

struct FieldSpec {
    REAL mu_x   = 1.0;          // média alvo do campo X
    REAL cov_x  = 0.0;          // coeficiente de variação (sigma_x / mu_x)
    FieldModel model = FieldModel::Gaussian;
};


// ------------------------------------------------------------------
// ComputeFieldMulti: lê N malhas-fonte (cada uma já com uma realização "z(x)")
// e escreve em memória elástica do alvo, montando X_i(x) com mu_x e cov_x.
// - Gaussian:   X_i = mu_x[i] + (cov_x[i]*mu_x[i]) * z_i
// - Lognormal:  X_i = exp( mu_ln[i] + sigma_ln[i] * z_i )
//   onde: sigma_ln = sqrt( ln(1 + COV_x^2) ), mu_ln = ln(mu_x) - 0.5*sigma_ln^2
// idxE/idxNu escolhem quais campos alimentam E e nu.
// ------------------------------------------------------------------
static void ComputeFieldMulti(const std::vector<TPZCompMesh*>& sources,
                              TPZCompMesh* target,
                              const std::vector<FieldSpec>& specs,
                              int idxE = 0, int idxNu = 1)
{
    if (!target) DebugStop();
    if (sources.empty()) DebugStop();
    if ((int)specs.size() < (int)sources.size()) DebugStop();

    auto *targetmatwithmem =
    dynamic_cast<TPZMatWithMem<TPZElasticMem>*>(target->MaterialVec()[1]);
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
    TPZAdmChunkVector<TPZElasticMem> &mem = *targetmatwithmem->GetMemory();

    const int nels = target->NElements();
    for (int iel=0; iel<nels; iel++)
    {

        auto *targetcel = target->ElementVec()[iel];

        if(!targetcel) continue;
        auto *targetintel = dynamic_cast<TPZInterpolationSpace*>(targetcel);
        if(!targetintel) continue;

        auto *intelmat = dynamic_cast<TPZMatWithMem<TPZElasticMem>*>(targetintel->Material());
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
            const REAL E  = (idxE  >=0 && idxE  < (int)X.size()) ? X[idxE]  : 1.0;
            const REAL nu = (idxNu >=0 && idxNu < (int)X.size()) ? X[idxNu] : 0.0;
            ER.SetEngineeringData(/*E=*/E, /*nu=*/nu);
            mem[indextarget].m_ER = ER;
            //std::cout << " E="<<ER.E() << " nu="<<ER.Poisson() << "\n";
        }

    }

    targetmatwithmem->SetUpdateMem(false);
}


// Consulta deslocamento Uy em ponto
static REAL UyAtPoint(TPZCompMesh* cmesh, REAL x, REAL y)
{
    TPZManVector<REAL,3> X(3,0.); X[0]=x; X[1]=y;
    TPZManVector<REAL,3> qsi(3,0.); int64_t elid=0;
    cmesh->LoadReferences();
    TPZGeoEl* gel = cmesh->Reference()->FindElement(X,qsi,elid,2);
    if(!gel || !gel->Reference()) DebugStop();
    auto* cel = dynamic_cast<TPZInterpolationSpace*>(gel->Reference());
    if(!cel) DebugStop();
    int var = cel->Material()->VariableIndex("Displacement");
    TPZVec<REAL> sol; cel->Solution(qsi, var, sol);
    return sol[1];
}

static REAL SolveOnce_Uy(TPZCompMesh* cm, REAL x, REAL y)
{
    TPZLinearAnalysis an(cm);
    TPZSkylineStructMatrix<STATE> str(cm);
    an.SetStructuralMatrix(str);
    TPZStepSolver<REAL> direct; direct.SetDirect(ELDLt);
    an.SetSolver(direct);
    an.Assemble();
    an.Solve();
    return UyAtPoint(cm,x,y);
}

// ------------------------------------------------------------
// Sweep Lx,Ly -> gera hhat_<Lx>_<Ly>.bin (buildfields=true)
// ou lê hhat e salva mc_results_<Lx>_<Ly>.csv (buildfields=false)
// ------------------------------------------------------------
int main()
{
    const int   kMatId   = 1;
    const int   porder   = 2;
    const int   nRef     = 3;
    const int   M        = 20;
    int         NsampGen = 30000;   // usado ao gerar hhat
    const REAL  xprobe   = -0.5;
    const REAL  yprobe   =  0.5;

    std::vector<std::pair<REAL,REAL>> grid = {
        {0.25,0.25},{0.25,0.50},{0.25,0.75},{0.25,1.0},
        {0.50,0.25},{0.50,0.50},{0.50,0.75},{0.50,1.0},
        {0.75,0.25},{0.75,0.50},{0.75,0.75},{0.75,1.0},
        {1.00,0.25},{1.00,0.50},{1.00,0.75},{1.00,1.0}
    };

    const bool buildfields = false; // true: gera hhat | false: roda MC e salva CSV

    if (buildfields)
    {
        for (auto [Lx,Ly] : grid)
        {
            TPZGeoMesh*  gmesh   = CreateGeoMeshMathematicaLike(kMatId, nRef);
            TPZCompMesh* cmeshKL = BuildCompMeshKL_Param(gmesh, porder, kMatId, Lx, Ly);

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

            // Duas famílias independentes para E e nu (poderia ser correlacionado se desejado)
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
    }
    else
    {
        for (auto [Lx,Ly] : grid)
        {
            // malhas fontes (para mapear campo -> elasticidade)
            TPZGeoMesh*  gmesh0       = CreateGeoMeshMathematicaLike(kMatId, nRef);
            TPZGeoMesh*  gmesh1        = CreateGeoMeshMathematicaLike(kMatId, nRef);
            TPZGeoMesh*  gmesh2        = CreateGeoMeshMathematicaLike(kMatId, nRef);
            TPZCompMesh* cmeshFieldE  = BuildCompMeshKL_Param(gmesh0, porder, kMatId, 1.0, 1.0);
            TPZCompMesh* cmeshFieldNu = BuildCompMeshKL_Param(gmesh1, porder, kMatId, 1.0, 1.0);
            TPZCompMesh* celastic     = CompMeshElastic(gmesh2, porder, kMatId);

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

           // hhatE.Print("E");
           // hhatMU.Print("mu");
            // quantas amostras usar (outra escolha possível: usar tudo)
            const int NsampRun = std::min(500, (int)std::min(hhatE.Cols(), hhatMU.Cols()));

            // abre CSV: mc_results_<Lx>_<Ly>.csv
            std::ostringstream fout; fout.setf(std::ios::fixed);
            fout << "mc_results_" << std::setprecision(6) << Lx << "_" << Ly << ".csv";
            std::ofstream outcsv(fout.str());
            outcsv.setf(std::ios::fixed); outcsv << std::setprecision(10);
            outcsv << "sample,uy\n";

            std::cout << "solving mc "<<fin.str()<< " NsampRun = "<< NsampRun<<std::endl;
            double sum=0.0, sum2=0.0;

            // buffers de coluna por amostra
            TPZFMatrix<REAL> colE(hhatE.Rows(),1), colNu(hhatMU.Rows(),1);

            for (int s=0; s<NsampRun; ++s)
            {
                for (int i=0;i<hhatE.Rows(); ++i)   colE(i,0)  = hhatE(i,s);
                for (int i=0;i<hhatMU.Rows();++i)   colNu(i,0) = hhatMU(i,s);

                //colE.Print("colE");
                // Carregar a *coluna da amostra* em cada malha fonte
                cmeshFieldE->LoadSolution(colE);
                cmeshFieldE->LoadReferences();
                cmeshFieldNu->LoadSolution(colNu);
                cmeshFieldNu->LoadReferences();

                // Vetor de fontes arbitrário (pode ter N>2 no futuro)
                std::vector<TPZCompMesh*> sources = { cmeshFieldE, cmeshFieldNu };

                // Especificações de cada campo (média, desvio, lognormal?)
                // Aqui: identidade (μ=0, σ=1, normal) pois hhat já é a realização.
                std::vector<FieldSpec> specs;
                specs.push_back( FieldSpec{ /*mu_x=*/1.0,  /*cov_x=*/0.20, FieldModel::Lognormal } );   // -> E
                specs.push_back( FieldSpec{ /*mu_x=*/0.25, /*cov_x=*/0.10, FieldModel::Lognormal } );  // -> nu

                // Escolha qual campo alimenta E e nu (índices em 'sources')
                ComputeFieldMulti(sources, celastic, specs, /*idxE=*/0, /*idxNu=*/1);
                // Após ComputeFieldMulti(...)
                // {
                //     celastic->LoadReferences();
                //     // pega o 1º elemento com material
                //     for (int el=0; el<celastic->NElements(); ++el){
                //         auto *cel = dynamic_cast<TPZInterpolationSpace*>(celastic->ElementVec()[el]);
                //         if(!cel) continue;
                //         auto *matmem = dynamic_cast<TPZMatWithMem<TPZElasticMem>*>(cel->Material());
                //         if(!matmem) continue;
                //         TPZIntPoints &intr = cel->GetIntegrationRule();
                //         TPZMaterialDataT<STATE> d; d.intLocPtIndex=0; cel->InitMaterialData(d);
                //         REAL w; (void)w;
                //         TPZManVector<REAL,3> intpt(3,0.);
                //         intr.Point(0, intpt, w);
                //         cel->ComputeRequiredData(d,intpt);
                //         int ig = d.intGlobPtIndex;
                //         const auto &mem = *matmem->GetMemory();
                //         auto ER = mem[ig].m_ER;
                //         std::cout << "sample " << s << " E="<<ER.E() << " nu="<<ER.Poisson() << "\n";
                //         break;
                //     }
                // }

                const REAL uy = SolveOnce_Uy(celastic, xprobe, yprobe);
                outcsv << s << "," << uy << "\n";
                sum  += uy;
                sum2 += (double)uy*(double)uy;
            }
            outcsv.close();

            const double mean = sum/NsampRun;
            const double var  = std::max(0.0, sum2/NsampRun - mean*mean);
            const double stdv = std::sqrt(var);
            std::cout << "[Lx="<<Lx<<", Ly="<<Ly<<"] Nsamp="<<NsampRun
            << " mean="<<mean<<" std="<<stdv
            << " -> salvo em " << fout.str() << "\n";

            delete celastic;
            delete cmeshFieldE;
            delete cmeshFieldNu;
            delete gmesh0;
            delete gmesh1;
            delete gmesh2;
        }
    }
    return 0;
}
