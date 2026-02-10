
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
//#include "TPZElasticResponse.h"
#include "Plasticity/TPZElastoPlasticMem.h"
#include "Plasticity/TPZPlasticStepPV.h"
#include "Plasticity/TPZYCMohrCoulombPV.h"
#include "Plasticity/TPZMatElastoPlastic2D.h"
#include "Plasticity/TPZElasticResponse.h"
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
#include "Plasticity/pzelastoplasticanalysis.h"
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
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <filesystem>
#include "Plasticity/TPZPlasticStepVoigt.h"
#include "Plasticity/TPZYCMohrCoulombPV2.h"

typedef TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse> TPlasticMC;
typedef TPZMatElastoPlastic2D<TPlasticMC, TPZElastoPlasticMem>   plasticmat;

typedef TPZPlasticStepVoigt<TPZYCMohrCoulombPV2, TPZElasticResponse> TPlasticStepVoigtMC;
typedef TPZMatElastoPlastic2D<TPlasticStepVoigtMC,TPZElastoPlasticMem> TMatElastoPlaticMC;

void RunImpactStudy_FS2();
REAL young   = 20000.;
REAL poisson = 0.49;

REAL coes    = 10.;
REAL atrito  = 30. * M_PI / 180.;
TPZManVector<REAL,3> bodyforce={0.,-20.,0.};


void RunDeterministic();

void RunImpactStudy_FS();
using namespace std;

// ------------------------------------------------------------
// Pós-processo
// ------------------------------------------------------------
void PostProcessVariables(TPZStack<std::string>& scal, TPZStack<std::string>& vec);

void CreatePostProcessingMesh(TPZCompMesh* cmesh,TPZPostProcAnalysis* pproc,int matid);

void PostElastoplastic(TPZCompMesh* cmesh,const std::string& vtkfile,int matid,int step,int dim);


TPZGeoMesh*  TriGMesh(int ref);

TPZCompMesh* CreateCMesh(TPZGeoMesh* gmesh, int pOrder, plasticmat* mat);

TPZCompMesh* CreateCMesh(TPZGeoMesh* gmesh, int pOrder);
//void InitializeMemory(TPZCompMesh* cmesh, REAL coesion, REAL atrito);

void ComputeElementDeformation(TPZCompMesh* cmesh, TPZVec<REAL>& fPlasticDeformSqJ2);

void ComputeElementDeformation2(TPZCompMesh* cmesh, TPZVec<REAL>& fIndicator);

void DivideElementsAbove(TPZCompMesh* cmesh, REAL refineaboveval, std::set<int64_t>& out_newels);
//void DivideElementsAbove(TPZCompMesh* cmesh, REAL refineaboveval,std::set<int64_t>& out_newels,int maxLevel=5);
void PRefineElementsAbove (TPZCompMesh* cmesh, REAL refineaboveval, std::set<int64_t>& out_newels,int porder );

bool Hrefine(TPZCompMesh* cmesh,REAL refineAboveVal);

bool HPrefine(TPZCompMesh* cmesh,REAL refineAboveVal,int porder);

REAL UyAtNode(TPZCompMesh* cmesh, REAL x, REAL y);

REAL FindFS_Bisection(TPZCompMesh* cmesh,
                      REAL lo, REAL hi,
                      REAL tol_fs_rel, int max_it,
                      int verbose,int loadmatid);

REAL Solve(TPZCompMesh* cmesh,REAL coes,REAL phi);

void ApplyLoad(TPZCompMesh* cmesh,
               REAL coes, REAL atrito, TPZManVector<REAL> factors);

bool RunAndAccept(TPZCompMesh* cmesh,REAL factor,int matid,bool post=false);
REAL IterativeProcessArcLength2(TPZElastoPlasticAnalysis &an,
                                int nsteps,
                                STATE lambda0,
                                STATE L0,
                                std::string vtkfile,int matid,STATE x,STATE y);

REAL SolveArc(TPZCompMesh* cmesh,int loadmatid,string vtkfile,int ref,STATE tol_fs_rel);
bool HPrefine(TPZCompMesh* cmesh,REAL refineAboveVal,int porder)
{
    const int nels_before = cmesh->NElements();
    TPZVec<REAL> defel;
    ComputeElementDeformation2(cmesh, defel);

    std::set<int64_t> novosp;
    PRefineElementsAbove (cmesh,refineAboveVal, novosp, porder );

    std::set<int64_t> novos;
    DivideElementsAbove(cmesh, refineAboveVal, novos);

    const int nels_after = cmesh->NElements();
    std::cout << "[HRefine] nels: " << nels_before << " -> " << nels_after
    << "  (refinados: " << (int)novos.size() << ")\n";

    std::cout << "[PRefine] nels: " << (int)novosp.size() << "\n";

    if (nels_after > nels_before) {
        return true;
    } else {
        std::cout << "[PreRefine] sem novos refinamentos; fim.\n";
        return false;
    }
}


bool Hrefine(TPZCompMesh* cmesh,REAL refineAboveVal)
{
    const int nels_before = cmesh->NElements();
    TPZVec<REAL> defel;
    ComputeElementDeformation2(cmesh, defel);

    std::set<int64_t> novos;
    DivideElementsAbove(cmesh, refineAboveVal, novos);

    const int nels_after = cmesh->NElements();
    std::cout << "[PreRefine] nels: " << nels_before << " -> " << nels_after
    << "  (refinados: " << (int)novos.size() << ")\n";

    if (nels_after > nels_before) {
        return true;
    } else {
        std::cout << "[PreRefine] sem novos refinamentos; fim.\n";
        return false;
    }
}

// Pequena utilidade para medir quão "fechado" está o bracket
static inline REAL RelGap(REAL a, REAL b) {
    const REAL m = (REAL)0.5*(a+b);
    return (b - a) / std::max<REAL>(m, (REAL)1e-12);
}


REAL FindFS_Bisection(TPZCompMesh* cmesh,
                      REAL lo, REAL hi,
                      REAL tol_fs_rel, int max_it,
                      int verbose,int loadmatid)
{
    if (hi < lo) std::swap(lo, hi);

    if (verbose) {
        std::cout << "\n[FS-Bisection] start"
        << "  lo=" << lo << "  hi=" << hi
        << "  tol_rel=" << tol_fs_rel
        << "  max_it=" << max_it <<std::endl;
    }

    int k = 0;
    while (k < max_it) {
        const REAL gap = RelGap(lo, hi);
        if (gap <= tol_fs_rel) {
            if (verbose) {
                std::cout << "[FS-Bisection] stop: gap=" << gap
                << " <= tol=" << tol_fs_rel
                << "  it=" << k <<std::endl;
            }
            break;
        }

        const REAL mid = (REAL)0.5*(lo + hi);
        int it_mid = 0;

        const bool ok = RunAndAccept(cmesh,  mid,loadmatid);

        if (verbose) {
            std::cout << "[FS-Bisection][it " << k << "] "
            << "mid=" << mid
            << " gap=" << gap
            << " -> " << (ok ? "OK" : "FAIL")
            <<std::endl;
        }

        if (ok) lo = mid; else hi = mid;
        ++k;
    }

    if (verbose) {
        std::cout << "[FS-Bisection] end  FS≈" << lo
        << "  gap_final=" << RelGap(lo,hi)
        << "  it=" << k <<std::endl;
    }
    return lo; // melhor piso convergente
}





TPZCompMesh* CreateCMesh(TPZGeoMesh* gmesh, int pOrder)
{
    TPZCompMesh* cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder);
    cmesh->SetDimModel(2);

    STATE phi=atrito;
    STATE psi=phi;
    STATE c =coes;
    TPZElasticResponse ER;
    ER.SetEngineeringData(young,poisson);
    auto mc = TPZYCMohrCoulombPV2( phi, psi, c,ER) ;
    TPlasticStepVoigtMC PlasticStepVoigt;
    PlasticStepVoigt.SetPlasticCriterion(mc);
    PlasticStepVoigt.SetElasticResponse(ER);


    int id=1;
    int planestrain=1;
    auto* material = new TMatElastoPlaticMC(id, planestrain);
    material->SetPlasticityModel(PlasticStepVoigt);
    material->SetId(id);
    material->SetBodyForce0(bodyforce);
    material->SetBodyForce(bodyforce);

    cmesh->InsertMaterialObject(material);

    TPZFMatrix<STATE> val1(2,2,0.0);
    TPZManVector<STATE,2> val2(2,0.0);

    int dir = 3;
    val2[0]=1; val2[1]=1; auto* bc0 = material->CreateBC(material, -1, dir, val1, val2);
    val2[0]=1; val2[1]=0; auto* bc1 = material->CreateBC(material, -2, dir, val1, val2);
    val2[0]=1; val2[1]=0; auto* bc2 = material->CreateBC(material, -5, dir, val1, val2);

    cmesh->InsertMaterialObject(bc0);
    cmesh->InsertMaterialObject(bc1);
    cmesh->InsertMaterialObject(bc2);

    cmesh->SetAllCreateFunctionsContinuousWithMem();
    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    return cmesh;
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

void ComputeElementDeformation2(TPZCompMesh* cmesh, TPZVec<REAL>& fIndicator)
{
    const int64_t nelem = cmesh->NElements();
    fIndicator.resize(nelem);
    fIndicator.Fill(0.0);

    // guardamos os 2 critérios separadamente
    TPZVec<REAL> vQ(nelem, 0.0);
    TPZVec<REAL> vP(nelem, 0.0);

    // 2 colunas: col0 = indicador combinado (para refino), col1 = J2, col2 = Umag (opcional)
    cmesh->ElementSolution().Redim(nelem, 3);

    REAL maxQ = 0.0;
    REAL maxP  = 0.0;

    for (int64_t el = 0; el < nelem; ++el) {
        TPZCompEl* cel = cmesh->ElementVec()[el];
        if (!cel) continue;
        if (cel->Dimension() != cmesh->Dimension()) continue;

        auto* matmem = dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem>*>(cel->Material());
        if (!matmem) continue;

        TPZManVector<int64_t> memindices;
        cel->GetMemoryIndices(memindices);
        if (memindices.size() == 0) continue;

        REAL vmaxQ = 0.0;
        REAL vmaxP  = 0.0;
        int  npts   = 0;

        for (int64_t midx : memindices) {
            if (midx < 0) continue;

            const auto& mem = matmem->MemItem(midx);

            const TPZTensor<REAL> epsp = mem.m_elastoplastic_state.EpsP();
            //const TPZTensor<REAL> sigma = mem.m_sigma;
            //REAL hard = mem.m_elastoplastic_state.;
            REAL J2 = epsp.J2();
            REAL I1=epsp.I1();
            if (J2 < (REAL)0) J2 = 0;
            const REAL eqq = std::sqrt(J2);


            const TPZVec<REAL> dispace = mem.m_u;
            const REAL umag = Norm(dispace);

            const REAL eqp = fabs(I1/3);

            //const REAL eqp = umag;
            vmaxQ = std::max(vmaxQ, eqq);
            vmaxP  = std::max(vmaxP,  eqp);

            ++npts;
        }

        if (npts > 0) {
            vQ[el] = vmaxQ;
            vP[el] = vmaxP;

            maxQ = std::max(maxQ, vmaxQ);
            maxP  = std::max(maxP,  vmaxP);
        }
    }

    // combina (com normalização) -> um único indicador p/ refino
    const REAL eps = 1e-12;
    const REAL invMaxQ = 1.0 / (maxQ + eps);
    const REAL invMaxP  = 1.0 / (maxP  + eps);

    // escolha da regra:
    // (1) produto (mais focado): ind = nj2 * nu
    // (2) soma ponderada:       ind = wj2*nj2 + wu*nu
    const REAL wj2 = 1.;
    const REAL wu  = 0.;

    for (int64_t el = 0; el < nelem; ++el) {
        const REAL nQ = vQ[el] * invMaxQ;
        const REAL nP  = vP[el] * invMaxP;

        //cout << "nQ = "<< nQ <<endl;
       // cout << "nP = "<< nP <<endl;
        // ---- TROQUE AQUI se quiser produto ----
        // fIndicator[el] = nj2 * nu;
        fIndicator[el] = wj2*nQ + wu*nP;

    }

    // publica o combinado na coluna 0 (o que seu HPrefine usa)
    cmesh->SetElementSolution(0, fIndicator);
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
            // REAL alpha = mem.m_elastoplastic_state.m_hardening;
            // const TPZTensor<REAL> sigma = mem.m_sigma;
            // const REAL J2sigma = sigma.J2();
            //vmax = std::max(vmax, eqp);
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

// void DivideElementsAbove(TPZCompMesh* cmesh, REAL refineaboveval,
//                          std::set<int64_t>& out_newels,
//                          int maxLevel /*ex: 3*/)
// {
//     cmesh->LoadReferences();
//     std::vector<int64_t> to_divide;
//     const TPZFMatrix<STATE>& elsol = cmesh->ElementSolution();
//
//     const int64_t ne0 = cmesh->NElements();
//     for (int64_t el=0; el<ne0; el++) {
//         TPZCompEl* cel = cmesh->ElementVec()[el];
//         if (!cel) continue;
//         auto* intel = dynamic_cast<TPZInterpolationSpace*>(cel);
//         if (!intel) continue;
//         if (!dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem>*>(cel->Material())) continue;
//
//         TPZGeoEl* gel = cel->Reference();
//         if (!gel) continue;
//
//         // >>> NOVO: bloqueio por nível (elemento e, indiretamente, seus filhos)
//         if (gel->Level() >= maxLevel) continue;
//
//         if (elsol.Rows() > el && elsol.Cols() > 0) {
//             if (elsol(el,0) > refineaboveval) to_divide.push_back(el);
//         }
//     }
//
//     for (auto el : to_divide) {
//         TPZCompEl* cel = cmesh->ElementVec()[el]; if (!cel) continue;
//         auto* intel = dynamic_cast<TPZInterpolationSpace*>(cel); if (!intel) continue;
//
//         TPZGeoEl* gel = cel->Reference(); if (!gel) continue;
//         // >>> NOVO: segurança extra
//         if (gel->Level() >= maxLevel) continue;
//
//         int p = intel->GetPreferredOrder();
//         TPZStack<int64_t> sub; const int64_t idx = cel->Index();
//         intel->Divide(idx, sub, /*create_boundary_elements=*/0);
//         for (int i=0;i<sub.size();i++){
//             out_newels.insert(sub[i]);
//             auto* sc = cmesh->ElementVec()[sub[i]];
//             if (auto* si = dynamic_cast<TPZInterpolationSpace*>(sc)) si->SetPreferredOrder(p);
//         }
//     }
//
//     // balance
//     bool changed = true;
//     while (changed){
//         changed = false;
//         std::set<int64_t> need;
//         const int64_t ne = cmesh->NElements();
//         for (int64_t el=0; el<ne; el++){
//             TPZCompEl* cel = cmesh->ElementVec()[el]; if (!cel) continue;
//             auto* intel = dynamic_cast<TPZInterpolationSpace*>(cel); if (!intel) continue;
//             TPZGeoEl* gel = cel->Reference(); if (!gel) continue;
//
//             const int ns = gel->NSides();
//             for (int s=0;s<ns;s++){
//                 TPZGeoElSide gs(gel,s);
//                 if (gs.Dimension() != gel->Dimension()-1) continue;
//                 TPZCompElSide big = gs.LowerLevelCompElementList2(/*onlyintersect=*/1);
//                 if (!big) continue;
//                 TPZGeoElSide gbig(big.Reference());
//                 if (gbig.Element()->Dimension() != gel->Dimension()) continue;
//                 if (gel->Level() - gbig.Element()->Level() > 1){
//                     need.insert(big.Element()->Index());
//                 }
//             }
//         }
//         for (auto el : need){
//             TPZCompEl* cel = cmesh->ElementVec()[el]; if (!cel) continue;
//             auto* intel = dynamic_cast<TPZInterpolationSpace*>(cel); if (!intel) continue;
//             if (!dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem>*>(cel->Material())) continue;
//
//             TPZGeoEl* gel = cel->Reference(); if (!gel) continue;
//             // >>> NOVO: não balancear refinando além do limite
//             if (gel->Level() >= maxLevel) continue;
//
//             int p = intel->GetPreferredOrder();
//             TPZStack<int64_t> sub; const int64_t idx = cel->Index();
//             intel->Divide(idx, sub, /*create_boundary_elements=*/0);
//             for (int i=0;i<sub.size();i++){
//                 out_newels.insert(sub[i]);
//                 auto* sc = cmesh->ElementVec()[sub[i]];
//                 if (auto* si = dynamic_cast<TPZInterpolationSpace*>(sc)) si->SetPreferredOrder(p);
//             }
//             changed = true;
//         }
//     }
//
//     cmesh->AdjustBoundaryElements();
//     cmesh->CleanUpUnconnectedNodes();
//     cmesh->InitializeBlock();
//     cmesh->ExpandSolution();
//
//     // reset memória em todos os materiais com memória
//     for (auto& kv : cmesh->MaterialVec()){
//         TPZMaterial* m = kv.second;
//         if (auto* mm = dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem>*>(m)){
//             mm->ResetMemory();
//         }
//     }
// }



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

bool RunAndAccept(TPZCompMesh* cmesh,REAL factor,int matid,bool post)
{
    auto* bodymat = dynamic_cast<TMatElastoPlaticMC*>(cmesh->FindMaterial(matid));
    auto* bcmat = dynamic_cast<TPZBndCondT<STATE>*>(cmesh->FindMaterial(matid));


    TPZManVector<REAL,3> f0;
    if(bodymat)
    {
        f0=bodymat->GetBodyForce();
        TPZManVector<REAL,3> fb=f0;
        fb[1]*=factor;
        //std::cout << "Factor = " << factor <<   "\n";
        bodymat->SetBodyForce(fb);
    }else
    {
        if(!bcmat)
        {
            std::cout << "material de contorno nao encontrado \n";
            DebugStop();
        }
        f0=bcmat->Val2();
        bcmat->Val2()[0] *= factor;
        bcmat->Val2()[1] *= factor;


    }

    cmesh->Solution().Zero();

    TPZElastoPlasticAnalysis anal(cmesh, std::cout,TPZElastoPlasticAnalysis::ELineSearch::Armijo);

    if(false)
    {
        TPZFStructMatrix<REAL> str(cmesh);
        anal.SetStructuralMatrix(str);
        TPZStepSolver<REAL> direct;
        direct.SetDirect(ELU);
        anal.SetSolver(direct);
    }else{
        TPZSkylineStructMatrix<STATE> matskl(cmesh);
        matskl.SetNumThreads(16);
        anal.SetStructuralMatrix(matskl);
        TPZStepSolver<STATE> step; step.SetDirect(ELDLt);
        anal.SetSolver(step);
    }
    bool ok = anal.NewtonRaphson(false);
    if(bodymat)
    {
        bodymat->SetBodyForce(f0);
    }else
    {
        bcmat->Val2()[0]=f0[0];
        bcmat->Val2()[1]=f0[1];
        bcmat->Val2()[2]=f0[2];
    }
    if (!ok) return false;
    if(post)
    {
        int dim=2;
        int matid=1;
        PostElastoplastic(cmesh,"post.vtk",matid,0,dim);
    }
    anal.AcceptSolution();
    return true;
}

REAL Solve(TPZCompMesh* cmesh,int loadmatid,string vtkfile,int ref,STATE tol_fs_rel)
{
    REAL lo=0.5;
    REAL hi=30.;
    int max_bis = 20;
    int verbose=1;
    REAL FS=1000.;
    REAL FSOLD=0.;

    int porder=cmesh->GetDefaultOrder();
    porder+=1;
    int iters_out;
    int maxref=ref;
    TPZStack<STATE> fsstack;
    for ( int iref=1; iref<maxref; iref++ ) {

        int neq=cmesh->NEquations();
        std::cout << "\n[solve] ===== Refinamento # "<< iref <<" ====="<<" neq = " <<neq <<std::endl;
        FSOLD=FS;
        FS=  FindFS_Bisection(cmesh, lo,  hi, tol_fs_rel ,  max_bis,verbose,loadmatid);


        fsstack.Push(FS);
        if(FSOLD<FS)
        {
            cout << "FSOLD<FS  "<< "FSOLD = " << FSOLD << " FS = "<< FS<<endl;
            cout  << "FS final = " << FSOLD<<endl;
            REAL resu,resf;
            RunAndAccept( cmesh,  FSOLD,loadmatid);
            return FSOLD;

        }
        //if(iref==maxref-1)break;
        //RunAndAccept( cmesh,  FS,loadmatid,false);
        HPrefine(cmesh,tol_fs_rel,porder+1);

        //porder+=1;
    }
    std::ofstream vtk("gmeshtrirefined.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(cmesh->Reference(), vtk, true);
    cout  << "FS final = " << FS<<endl;
    RunAndAccept( cmesh,  1.,loadmatid);
    return FS;

}
REAL FindFS_ArcLength(TPZCompMesh* cmesh,int loadmatid)
{

    auto* bodymat = dynamic_cast<TMatElastoPlaticMC*>(cmesh->FindMaterial(loadmatid));
    auto* bcmat = dynamic_cast<TPZBndCondT<STATE>*>(cmesh->FindMaterial(loadmatid));


    TPZManVector<REAL,3> f0;

    cmesh->Solution().Zero();

    TPZElastoPlasticAnalysis anal(cmesh, std::cout,TPZElastoPlasticAnalysis::ELineSearch::Armijo);

    TPZSkylineStructMatrix<STATE> matskl(cmesh);
    matskl.SetNumThreads(16);
    anal.SetStructuralMatrix(matskl);
    TPZStepSolver<STATE> step; step.SetDirect(ELDLt);
    anal.SetSolver(step);

    std::string vtkfile2="vtkfilearc.vtk";
    int nsteps=5;
    STATE lambda0=0.01;
    STATE L0=1.;
    STATE x=40.;
    STATE y=30.;
    REAL FS = IterativeProcessArcLength2(anal,nsteps,lambda0,L0,vtkfile2,loadmatid,x,y);

    return FS;
}

REAL SolveArc(TPZCompMesh* cmesh,int loadmatid,int ref,STATE tol_fs_rel)
{
    REAL FS=1000.;
    REAL FSOLD=0.;

    int porder=cmesh->GetDefaultOrder();

    int maxref=ref;
    TPZStack<STATE> fsstack;
    for ( int iref=1; iref<maxref; iref++ ) {

        int neq=cmesh->NEquations();
        std::cout << "\n[solve] ===== Refinamento # "<< iref <<" ====="<<" neq = " <<neq <<std::endl;
        FSOLD=FS;
        FS=  FindFS_ArcLength(cmesh,loadmatid);
        fsstack.Push(FS);
        REAL deltaFS  = fabs(FS-FSOLD);
        if(deltaFS<1.e-3 && iref>1)
        {
            cout << "deltaFS  "<< deltaFS  << "FSOLD = " << FSOLD << " FS = "<< FS <<endl;
            std::ofstream vtk("gmeshtrirefined.vtk");
            TPZVTKGeoMesh::PrintGMeshVTK(cmesh->Reference(), vtk, true);
            return FS;
        }else if(iref!=maxref-1){
            HPrefine(cmesh,tol_fs_rel,porder+1);
        }


    }

    cout << "FAILED TO CONVERGE FS IN HP REFINE " <<endl;

    return FS;

}

void PostProcessVariables(TPZStack<std::string>& scal, TPZStack<std::string>& vec)
{
    scal.Push ( "StrainPlasticJ2" );
    vec.Push ( "DisplacementDoF" );
    scal.Push ( "EBodyForce" );
    scal.Push ( "StressXX" );
    scal.Push ( "StressYY" );
    scal.Push ( "StressZZ" );
    scal.Push ( "StrainPlasticZZ" );
    scal.Push ( "StrainTotalZZ" );
    scal.Push ( "DamageVariable" );

     scal.Push ( "EOrder" );
    // scal.Push ( "POrder" );
    // scal.Push ( "Atrito" );
    // scal.Push ( "Coesion" );
    // scal.Push ( "VolHardening" );
    // vec.Push ( "ShearPlasticDeformation" );
    // vec.Push ( "PlasticDeformation" );

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

void PostElastoplastic(TPZCompMesh* cmesh,const std::string& vtkfile,int matid,int step,int dim)
{
    TPZPostProcAnalysis pproc;
    CreatePostProcessingMesh(cmesh, &pproc, matid);
    TPZStack<std::string> scal, vec;
    PostProcessVariables(scal, vec);
    pproc.DefineGraphMesh(/*dim=*/dim, scal, vec, vtkfile);
    pproc.SetStep(step);
    pproc.PostProcess(0);
}


void RunDeterministic()
{
    int pOrder = 2;
    int ref =0;
    TPZGeoMesh* gmesh = TriGMesh(ref);
    {
        std::ofstream vtk1("antes.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, vtk1, true);
        std::cout << "antes.vtk escrito.\n";
    }

    TPZCompMesh* cmesh = CreateCMesh(gmesh, pOrder);

    using Clock = std::chrono::steady_clock;

    auto t0 = Clock::now();
    std::string vtkfile2="Slope.vtk";
    int loadid=1;
    int iref=10;
    STATE tolref=0.2;
   // REAL FS = Solve(cmesh,loadid,vtkfile2,iref,tolref);

    REAL FS = SolveArc(cmesh,loadid,iref,tolref);

    auto t1 = Clock::now();

    std::chrono::duration<double> secs = t1 - t0;
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();

    std::cout << std::fixed << std::setprecision(3) << "[Timing] Solve: " << secs.count() << " s  (" << ms << " ms)\n";
    {
        std::ofstream vtk1("gmeshtri_refined_preGI.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(cmesh->Reference(), vtk1, true);
        std::cout << "[VTK] gmeshtri_refined_preGI.vtk escrito.\n";
    }

    int matid=1;
    int dim=2;
    std::string vtkfile3="SlopeX.vtk";
    PostElastoplastic(cmesh, vtkfile3, matid, /*out_step*/ 0 ,dim);

    //
    //
    // TPZElastoPlasticAnalysis an(cmesh, std::cout,TPZElastoPlasticAnalysis::ELineSearch::Armijo);
    // TPZSkylineStructMatrix<STATE> matskl(cmesh);
    // matskl.SetNumThreads(16);
    // an.SetStructuralMatrix(matskl);
    // TPZStepSolver<STATE> step;
    // step.SetDirect(ELDLt);
    // an.SetSolver(step);
    // int nsteps=20;
    // STATE lambda0=0.0001;
    // STATE L0=1.;
    // STATE x=30.;
    // STATE y=45;
    // std::string vtkfile4="SlopeARC.vtk";

}

int main()
{
    // RunStochastic(false);
    RunDeterministic();
   // RunImpactStudy_FS2();

    return 0;
}


// ============================================================
// ImpactStudy FS2 - versão corrigida e robusta (NeoPZ style)
//  - Seleciona top-M por indicador (ind por CompEl index)
//  - Converte cada candidato para coordenada global do centro do GeoEl
//  - Em cada caso: recria gmesh, acha o elemento por coordenada, refina 1-ring, calcula FS
//  - Faz ranking por |FSK-FS0|
//  - Pega o melhor caso e faz +1 nível de refinamento extra nesse mesmo patch
// ============================================================

#include <set>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <fstream>
#include <chrono>
#include <cmath>
#include <iostream>

// ------------------------------------------------------------
// Dependências que você já tem no projeto:
//  - TriGMesh(int ref0)
//  - CreateCMesh(TPZGeoMesh*, int pOrder)
//  - ComputeElementDeformation(TPZCompMesh*, TPZVec<REAL>& ind)
//  - FindFS_Bisection(TPZCompMesh*, REAL lo, REAL hi, REAL tol_rel, int max_it, int verbose, int loadid)
//  - Solve(TPZCompMesh*, int loadid, std::string vtkfile, int nref_solve, STATE tolref)
//  - TPZVTKGeoMesh::PrintGMeshVTK(...)
// ------------------------------------------------------------

struct ImpactResult {
    int icase = -1;
    TPZManVector<REAL,3> X = TPZManVector<REAL,3>(3,0.0); // coord global do candidato (centro)
    REAL FS = 0.0;
    REAL dFS = 0.0;
};

// --------- util: pega top-M índices por indicador ----------
static std::vector<int64_t> TopM_ByIndicator(const TPZVec<REAL>& ind, int M)
{
    std::vector<std::pair<REAL,int64_t>> v;
    v.reserve(ind.size());
    for (int64_t i = 0; i < (int64_t)ind.size(); ++i) v.push_back({ind[i], i});

    std::sort(v.begin(), v.end(),
              [](const std::pair<REAL,int64_t>& a, const std::pair<REAL,int64_t>& b)
              { return a.first > b.first; });

    std::vector<int64_t> out;
    out.reserve(std::min<int>(M, (int)v.size()));
    for (int i = 0; i < M && i < (int)v.size(); ++i) {
        if (v[i].first <= 0) break;
        out.push_back(v[i].second);
    }
    return out;
}

// --------- util: dado X, acha gel e monta patch (1-ring) ----
static std::set<int64_t> CollectPatchGeoEls(TPZGeoMesh* gmesh,
                                            const TPZManVector<REAL,3>& X,
                                            int dim)
{
    std::set<int64_t> patch;
    if(!gmesh) return patch;

    TPZManVector<REAL,3> qsi(2,0.0);
    int64_t elindex=0;
    TPZManVector<REAL,3> Xintternal=X;
   // TPZGeoEl* gel = gmesh->FindApproxElement(Xintternal, qsi, elindex, dim);
    TPZGeoEl* gel = gmesh->FindCloseElement(Xintternal, elindex, dim);
    //TPZGeoEl* gel = gmesh->FindElement(Xintternal, qsi, elindex, dim);
    if(!gel) return patch;

    patch.insert(gel->Index());

    const int ns = gel->NSides();
    for (int s = 0; s < ns; s++){
        TPZGeoElSide gs(gel, s);
        if(!gs.Element()) continue;

        TPZStack<TPZGeoElSide> neigh;
        gs.AllNeighbours(neigh);

        for (int i = 0; i < neigh.size(); i++){
            TPZGeoEl* ng = neigh[i].Element();
            if (!ng) continue;
            if (ng->Dimension() != gel->Dimension()) continue;
            patch.insert(ng->Index());
        }
    }
    return patch;
}

static std::set<int64_t> CollectPatchGeoEls_Rings(TPZGeoMesh* gmesh,
                                                  const TPZManVector<REAL,3>& X_in,
                                                  int dim,
                                                  int rings)
{
    std::set<int64_t> patch;
    if(!gmesh || rings < 0) return patch;

    // FindElement costuma exigir X não-const
    TPZManVector<REAL,3> X = X_in;

    TPZManVector<REAL,3> qsi(dim, 0.0);   // CRÍTICO: tamanho = dim
    int64_t elindex = 0;

    TPZGeoEl* gel0 = gmesh->FindElement(X, qsi, elindex, dim);
   // TPZGeoEl* gel0 = gmesh->FindCloseElement(X, elindex, dim);
    if(!gel0) return patch;

    // Sempre inclui o elemento central
    patch.insert(gel0->Index());

    if(rings == 0) return patch; // só o elemento que contém X

    // BFS por níveis (rings)
    std::set<int64_t> visited;
    visited.insert(gel0->Index());

    TPZStack<TPZGeoEl*> frontier;
    frontier.Push(gel0);

    for(int level = 0; level < rings; level++)
    {
        TPZStack<TPZGeoEl*> next_frontier;

        for(int ifr = 0; ifr < frontier.size(); ifr++)
        {
            TPZGeoEl* gel = frontier[ifr];
            if(!gel) continue;

            const int ns = gel->NSides();
            for(int s = 0; s < ns; s++)
            {
                TPZGeoElSide gs(gel, s);
                if(!gs.Element()) continue;

                TPZStack<TPZGeoElSide> neigh;
                gs.AllNeighbours(neigh);

                for(int i=0; i<neigh.size(); i++)
                {
                    TPZGeoEl* ng = neigh[i].Element();
                    if(!ng) continue;

                    // mantém só elementos da mesma dimensão "principal"
                    if(ng->Dimension() != gel0->Dimension()) continue;

                    const int64_t idx = ng->Index();
                    if(visited.insert(idx).second)
                    {
                        patch.insert(idx);
                        next_frontier.Push(ng);
                    }
                }
            }
        }

        frontier = next_frontier;
        if(frontier.size() == 0) break; // acabou o grafo
    }

    return patch;
}


// --------- util: refinar geometricamente um conjunto ---------
static void HRefineGeoSet_OneLevel(TPZGeoMesh* gmesh, const std::set<int64_t>& to_refine)
{
    if (!gmesh) return;

    // Um nível: divide apenas os elementos atuais (sem repetir no mesmo call)
    for (auto gelidx : to_refine) {
        if (gelidx < 0 || gelidx >= gmesh->ElementVec().NElements()) continue;

        TPZGeoEl* gel = gmesh->ElementVec()[gelidx];
        if (!gel) continue;
        if (gel->HasSubElement()) continue;

        TPZVec<TPZGeoEl*> sub;
        gel->Divide(sub);
    }
    gmesh->BuildConnectivity();
}
static void HRefineGeoSet_OneLevel_Limited(TPZGeoMesh* gmesh,
                                           const std::set<int64_t>& to_refine,
                                           int max_level = 2)
{
    if(!gmesh) return;

    for (auto gelidx : to_refine)
    {
        TPZGeoEl* gel = gmesh->ElementVec()[gelidx];
        if(!gel) continue;

        // já foi refinado "demais" -> não mexe
        if(gel->Level() >= max_level) continue;

        // se já tem filhos, não divide de novo (para 1-level)
        if(gel->HasSubElement()) continue;

        TPZVec<TPZGeoEl*> sub;
        gel->Divide(sub);
    }

    gmesh->BuildConnectivity();
}

// --------- util: cel -> coord global do centro do gel --------
// (corrige o erro: CenterPoint devolve qsi, não X global)
static std::vector<TPZManVector<REAL,3>>
MapTopCompElsToGeoCenterCoords(TPZCompMesh* cmesh,
                               const std::vector<int64_t>& top_cel_indices)
{
    std::vector<TPZManVector<REAL,3>> coords;
    if(!cmesh) return coords;

    auto &elvec = cmesh->ElementVec();
    coords.reserve(top_cel_indices.size());

    for (auto celidx : top_cel_indices) {

        if (celidx < 0 || celidx >= elvec.NElements()) continue;

        TPZCompEl* cel = elvec[celidx];
        if (!cel) continue;

        TPZGeoEl* gel = cel->Reference();
        if (!gel) continue;

        if (gel->Dimension() != cmesh->Dimension()) continue;

        TPZManVector<REAL,3> qsi(2,0.0);
        //const int side = gel->NSides() - 1; // interior do elemento (genérico)
        // const int side = 0; // interior do elemento (genérico)
        // gel->CenterPoint(side, qsi);

        qsi[0]=0.;
        qsi[1]=0.;
        TPZManVector<REAL,3> X(3,0.0);
        gel->X(qsi, X);

        coords.push_back(X);
    }

    return coords;
}
void RunImpactStudy_FS2()
{
    using Clock = std::chrono::steady_clock;

    const int pOrder    = 2;
    const int ref0      = 1;
    const int loadid    = 1;

    const REAL lo_fs    = 0.5;
    const REAL hi_fs    = 30.0;
    const REAL tol_rel  = 0.01;
    const int  max_bis  = 10;
    const int  verbose  = 0;

    const int M =10;

    std::vector<TPZManVector<REAL,3>> top_X;
    std::vector<TPZManVector<REAL,3>> best_history;  // <<< NOVO
    int rings = 2; // 1=vizinho, 2=vizinho do vizinho
    for (int ilevel=0; ilevel<10; ilevel++)
    {
        // ================= BASELINE (malha top acumulada) =================
        TPZGeoMesh* gmesh0 = TriGMesh(ref0);

        // reaplica TODOS os refinamentos vencedores anteriores
        for (auto &Xbest : best_history)
        {

            //auto patch_best = CollectPatchGeoEls_Rings(gmesh0, Xbest, /*dim=*/2, rings);
            auto patch_best = CollectPatchGeoEls(gmesh0, Xbest, /*dim=*/2);
            HRefineGeoSet_OneLevel(gmesh0, patch_best);
            //HRefineGeoSet_OneLevel_Limited(gmesh0, patch_best);
        }
        // >>>>>>>>>>>>>>> AQUI IMPRIME A MALHA ACUMULADA <<<<<<<<<<<<<<
        {
            std::string vtkname =
            "baseline_ilevel_" + std::to_string(ilevel) + ".vtk";
            std::ofstream vtk(vtkname.c_str());
            TPZVTKGeoMesh::PrintGMeshVTK(gmesh0, vtk, true);
            std::cout << "[ImpactStudy] escreveu " << vtkname << std::endl;
        }
        TPZCompMesh* cmesh0 = CreateCMesh(gmesh0, pOrder);

        // const REAL FS0 = FindFS_Bisection(
        //     cmesh0, lo_fs, hi_fs, tol_rel, max_bis, verbose, loadid);
        //
        //
        // const REAL FS0arc = FindFS_ArcLength(cmesh0,loadid);

/*
        std::cout << "\n[ImpactStudy] ilevel="<<ilevel
        << "  FS0 (baseline) = " << FS0  << "FS0 (ARC)"<<FS0arc << std::endl;*/



        const REAL u0arc = FindFS_ArcLength(cmesh0,loadid);
        std::cout << "\n[ImpactStudy] ilevel="<<ilevel
        << "  FS0 (baseline) = " << u0arc   << std::endl;
        // ================= CANDIDATOS =================
        TPZVec<REAL> ind;
        ComputeElementDeformation2(cmesh0, ind);

        auto top_cel = TopM_ByIndicator(ind, M);
        top_X = MapTopCompElsToGeoCenterCoords(cmesh0, top_cel);

        // ================= TESTAR IMPACTO =================
        std::vector<double> dFSVec;
        dFSVec.reserve(top_X.size());

        for (int icase=0; icase<(int)top_X.size(); icase++)
        {
            TPZGeoMesh* gmeshK = TriGMesh(ref0);

            // aplica baseline acumulado
            for (auto &Xbest : best_history)
            {
                auto patch_best = CollectPatchGeoEls(gmeshK, Xbest, /*dim=*/2);
               // auto patch_best = CollectPatchGeoEls_Rings(gmeshK, Xbest, /*dim=*/2, rings);
                HRefineGeoSet_OneLevel(gmeshK, patch_best);
                //HRefineGeoSet_OneLevel_Limited(gmeshK, patch_best);
            }

            // aplica patch do candidato
            auto patch_cand = CollectPatchGeoEls(gmeshK, top_X[icase], /*dim=*/2);
           // auto patch_cand = CollectPatchGeoEls_Rings(gmeshK, top_X[icase], /*dim=*/2, rings);
            HRefineGeoSet_OneLevel(gmeshK, patch_cand);
            //HRefineGeoSet_OneLevel_Limited(gmeshK, patch_cand);
            TPZCompMesh* cmeshK = CreateCMesh(gmeshK, pOrder);

            // const REAL FSK = FindFS_Bisection(
            //     cmeshK, lo_fs, hi_fs, tol_rel, max_bis, verbose, loadid);
            //
            // const REAL FSKarc = FindFS_ArcLength(cmeshK,loadid);

            const REAL uKarc = FindFS_ArcLength(cmeshK,loadid);

            const STATE dFS = std::abs(uKarc - u0arc);

            // const STATE dFS = std::abs(FSK - FS0);
            //
            // const STATE dFSKarc = std::abs(FSKarc - FS0);

            // std::cout << "\n[ImpactStudy] ilevel="<<ilevel
            // << " case="<<icase
            // << " X=("<<top_X[icase][0]<<","<<top_X[icase][1]<<","<<top_X[icase][2]<<")"
            // << " FSK="<<FSK
            // << " FSKarc="<<FSKarc
            // << " dFS="<<dFS
            // << " dFSKarc="<<dFSKarc
            // << std::endl;

            std::cout << "\n[ImpactStudy] ilevel="<<ilevel
            << " case="<<icase
            << " X=("<<top_X[icase][0]<<","<<top_X[icase][1]<<","<<top_X[icase][2]<<")"
            << " FSK="<<uKarc
            << " dFS="<<dFS
            << std::endl;

            dFSVec.push_back(dFS);

            delete cmeshK;
            delete gmeshK;
        }

        // ================= CAMPEÃO DO NÍVEL =================
        int indextop = (int)std::distance(
            dFSVec.begin(),
                                          std::max_element(dFSVec.begin(), dFSVec.end())
        );

        // auto rit = std::max_element(dFSVec.rbegin(), dFSVec.rend());
        // int indextop = (int)std::distance(dFSVec.begin(), rit.base()) - 1;
        std::cout << "[ImpactStudy] campeão do nível = "
        << indextop << std::endl;

        // guarda coordenada vencedora no histórico
        best_history.push_back(top_X[indextop]);

        delete cmesh0;
        delete gmesh0;
    }
}

#include <fstream>
#include <iomanip>
#include <vector> // <-- ADICIONE ISTO
// Eq. (4.123) – passo preditor (k = 1)
static STATE compute_dlambda0_riks(const TPZFMatrix<STATE>& dwb,
                                   const TPZFMatrix<STATE>& dw,
                                   STATE L)
{
    const STATE s    = Dot(dw, dwb);                 // Δu^T * dū
    const STATE ndwb = Norm(dwb);
    const STATE signum = (s > 0.0 ? -1.0 : 1.0);     // Souza Neto 4.123
    return signum * L / ndwb;
}

// Eqs. (4.116) + (4.118) – escolha da raiz para k > 1
static STATE compute_dlambda_riks(const TPZFMatrix<STATE>& dwb,
                                  const TPZFMatrix<STATE>& dws,
                                  const TPZFMatrix<STATE>& dw,
                                  STATE L, int& rootIdx)
{
    const STATE aa = Dot(dwb, dwb);
    TPZFMatrix<STATE> t = dw;
    t += dws;               // t = dw + dws
    const STATE bb = 2.0 * Dot(dwb, t);
    const STATE cc = Dot(t,t) - L*L;

    const STATE eps = 1e-14;
    if (aa < eps) {                                   // cai para linear
        rootIdx = 1;
        return (std::fabs(bb) > eps) ? (-cc/bb) : 0.0;
    }

    STATE disc = bb*bb - 4.0*aa*cc;
    if (disc < 0.0) disc = 0.0;                       // clamp numérico
    const STATE sq = std::sqrt(disc);

    const STATE dl1 = (-bb - sq) / (2.0*aa);          // “menor”
    const STATE dl2 = (-bb + sq) / (2.0*aa);          // “maior”

    auto score = [&](STATE dl)->STATE {
        TPZFMatrix<STATE> x = dw;                     // Δu^(k-1)
        TPZFMatrix<STATE> tmp = dwb; tmp *= dl;
        x += dws; x += tmp;                           // Δu^(k-1)+δu*+δλ dū
        return Dot(x, dw);                            // maximiza (4.118)
    };
    const STATE s1 = score(dl1), s2 = score(dl2);
    if (s1 > s2) { rootIdx = 1; return dl1; }
    else         { rootIdx = 2; return dl2; }
}
REAL IterativeProcessArcLength2(TPZElastoPlasticAnalysis &an,
                                int nsteps,
                                STATE lambda0,
                                STATE L0,
                                std::string vtkfile,int matid,STATE x,STATE y)
{
    // ======== ARQUIVOS DE SAÍDA (reuso dos nomes p/ Python) ========
    std::ofstream out_ld("arc_load_displacement.txt");
    out_ld << std::scientific << std::setprecision(15);
    out_ld << "# step lambda uy\n";

    std::ofstream out_res("arc_residuals.txt");
    out_res << std::scientific << std::setprecision(15);
    out_res << "# step iter lambda L(reserved) res_norm\n";

    auto cmesh = an.Mesh();
    cmesh->Solution().Zero();

    auto* bodymat = dynamic_cast<TMatElastoPlaticMC*>(cmesh->FindMaterial(matid));
    auto* bcmat = dynamic_cast<TPZBndCondT<STATE>*>(cmesh->FindMaterial(matid));

    auto* matmem = dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem>*>(cmesh->FindMaterial(1));
    matmem->ResetMemory();

    TPZManVector<REAL,3> f0;
    if(bodymat)
    {
        std::cout << "material de volume ENCONTRADO \n";
    }else
    {
        if(!bcmat)
        {
            std::cout << "material de contorno nao encontrado \n";
            DebugStop();
        }
        std::cout << "material de contorno ENCONTRADO \n";

    }
    REAL uy=0.;

    TPZFMatrix<STATE> FEXT;

    if(bodymat)
    {
        bodymat->SetBodyForce(bodymat->GetBodyForce0());
        an.Assemble();
        FEXT = an.Rhs();
        TPZManVector<REAL,3> f0={0,0,0};
        bodymat->SetBodyForce(f0);
    }else
    {
        if(!bcmat)
        {
            std::cout << "material de contorno nao encontrado \n";
            DebugStop();
        }
        an.Assemble();
        FEXT = an.Rhs();
        bcmat->Val2()[0]=0.;
        bcmat->Val2()[1]=0.;
        bcmat->Val2()[2]=0.;
    }

    cmesh->Solution().Zero();
    int dim = cmesh->Dimension();
    REAL lambda = lambda0;   // lambda corrente (pode ser alterado durante os passos)
    REAL lambdan = lambda;   // lambda aceito (para rollback seguro)
    REAL L = L0;             // comprimento de arco atual

    // Solução acumulada aceita (u_acc) e incremento acumulado do passo (dw)
    TPZFMatrix<STATE> u_acc = an.Solution();
    u_acc.Zero();

    int step = 0;
    const int max_cuts = 12;
    int cut_count = 0;

    while (step < nsteps) {

        // Estado de trabalho parte do aceito
        TPZFMatrix<STATE> u  = u_acc;
        bool okconv = false;
        REAL normR  = 1e30;
        TPZFMatrix<STATE> dw = an.Solution();
        dw.Zero();

        {
            const int  maxit_inner = 20;
            const REAL etol_inner  = 1e-8;

            TPZFMatrix<STATE> rhs, R, dws, dwb;

            int it = 0;
            int rootIdx = 0;

            // ----- LOG de resíduos deste step (armazenado, só imprime se convergir)
            struct ResEntry { int iter; STATE lambda; STATE Ldummy; STATE res; };
            std::vector<ResEntry> res_log; // <-- NOVO

            std::cout << " \n step = "<< step <<"\n";
            while (it < maxit_inner && normR > etol_inner)
            {
                an.LoadSolution(u);
                an.Assemble();
                rhs = an.Rhs();
                R = FEXT * lambda;
                R += rhs;

                an.Rhs() = R;
                an.Solve();
                dws = an.Solution();

                an.Rhs() = FEXT;
                an.Solve();
                dwb = an.Solution();

                REAL dl = (it == 0) ? compute_dlambda0_riks(dwb, dw, L)
                : compute_dlambda_riks(dwb, dws, dw, L, /*rootIdx*/ *(int[]){0});

                const REAL dl_max = 1;
                if (dl >  dl_max) dl =  dl_max;
                if (dl < -dl_max) dl = -dl_max;

                TPZFMatrix<STATE> dwtot = dwb*dl + dws;

                u      += dwtot;
                dw     += dwtot;
                lambda += dl;

                an.LoadSolution(u);
                an.Assemble();

                rhs   = an.Rhs();
                R     = FEXT * lambda;
                R    += rhs;
                normR = Norm(R);

                STATE normRFEXT = Norm(R)/Norm(FEXT);
                // ---- Guarda o residual desta iteração (só será impresso se o step convergir)
                res_log.push_back({it, lambda, L, normR}); // <-- NOVO

                cout << " iter = " << it
                << "  lambda = " << lambda
                << " dl = "     << dl
                << "  L = "     << L
                //<< " Norm(dw) = "  << Norm(dwtot)
                //<< " normRFEXT = "  << normRFEXT
                << " normR = "  << normR << std::endl;
                ++it;

                // if (normR > 1e5) {
                //     break;
                // }
            }

            okconv = (normR <= etol_inner);

            if (okconv) {
                // ---- imprime resíduos deste step (apenas em caso de CONVERGÊNCIA)
                for (const auto &r : res_log) {
                    out_res << step      << " "
                    << (r.iter+0) << " "
                    << r.lambda   << " "
                    << r.Ldummy   << " "
                    << r.res      << std::endl;
                }
            }
        }

        if (okconv) {
            // aceita o passo
            cmesh->LoadSolution(u);
            uy+=  UyAtNode(cmesh, x,y);
            out_ld << step   << " " << lambda << " " << uy << "\n";
            an.Assemble();
            an.AcceptSolution();
            u_acc   = an.Solution();

            // salva saída VTK (função externa no seu projeto)
            PostElastoplastic(cmesh, vtkfile, /*matid*/1, /*out_step*/ step , dim);

            // pronto para o próximo step
            ++step;
            cut_count = 0;

            // atualiza “lambda aceito” para rollback seguro
            if(fabs(lambda-lambdan)<1.e-3)break;
            lambdan = lambda;
        } else {
            // rollback e corta L
            an.LoadSolution(u_acc);
            lambda = lambdan;   // volta lambda aceito
            L *= 0.5;
            ++cut_count;

            if (cut_count > 12 || L < 1e-14) {
                std::cout << "[ArcLength] Falha em convergir no step " << step
                << " após " << cut_count << " cortes de L. Abortando ciclo.\n";
                break; // evita loop infinito do ciclo
            }

            // tenta novamente o MESMO step com L menor
            continue;
        }

    } // fim while step < nsteps

    an.AcceptSolution();
    return lambda;
}


// void RunImpactStudy_FS2()
// {
//     using Clock = std::chrono::steady_clock;
//
//     const int pOrder    = 2;
//     const int ref0      = 1;
//     const int loadid    = 1;
//
//     const REAL lo_fs    = 0.5;
//     const REAL hi_fs    = 30.0;
//     const REAL tol_rel  = 0.015;
//     const int  max_bis  = 20;
//     const int  verbose  = 0;
//
//     const int M = 4;
//
//     int indextop = 0; // precisa iniciar
//     std::vector<TPZManVector<REAL,3>> top_X;
//
//     for (int ilevel=0; ilevel<4; ilevel++)
//     {
//         // ---------------- baseline do nível (malha best acumulada) ----------------
//         TPZGeoMesh* gmesh0 = TriGMesh(ref0);
//
//         if (ilevel > 0)
//         {
//             cout << "indextop = "<<indextop<<" X=("<<top_X[indextop][0]<<","<<top_X[indextop][1]<<","<<top_X[indextop][2]<<")"<< std::endl;
//             auto patch_best = CollectPatchGeoEls(gmesh0, top_X[indextop], /*dim=*/2);
//             HRefineGeoSet_OneLevel(gmesh0, patch_best);
//         }
//
//
//         std::ofstream vtk(("impact_gmesh_" + std::to_string(indextop) +
//         "_ilevel_" + std::to_string(ilevel) + ".vtk").c_str());
//         TPZVTKGeoMesh::PrintGMeshVTK(gmesh0, vtk, true);
//
//         TPZCompMesh* cmesh0 = CreateCMesh(gmesh0, pOrder);
//
//         const REAL FS0 = FindFS_Bisection(cmesh0, lo_fs, hi_fs, tol_rel, max_bis, verbose, loadid);
//         std::cout << "\n[ImpactStudy] ilevel="<<ilevel<<"  FS0 (baseline) = " << FS0 << std::endl;
//
//         // ---------------- candidatos calculados NA malha baseline ----------------
//         TPZVec<REAL> ind;
//         ComputeElementDeformation(cmesh0, ind);
//
//         auto top_cel = TopM_ByIndicator(ind, M);
//
//         // se você já trocou para âncora (side=0), deixe sua função Map... retornar esse ponto
//         top_X = MapTopCompElsToGeoCenterCoords(cmesh0, top_cel);
//
//         // ---------------- impacto: comparar em cima do MESMO baseline ----------------
//         std::vector<double> dFSVec;
//         dFSVec.reserve(top_X.size());
//
//         for (int icase=0; icase<(int)top_X.size(); icase++)
//         {
//             TPZGeoMesh* gmeshK = TriGMesh(ref0);
//
//             // aplica o best acumulado também nos casos (senão dFS fica inconsistente)
//             if (ilevel > 0)
//             {
//                 auto patch_best = CollectPatchGeoEls(gmeshK, top_X[indextop], /*dim=*/2);
//                 HRefineGeoSet_OneLevel(gmeshK, patch_best);
//             }
//
//             // aplica o patch do candidato
//             auto patch_cand = CollectPatchGeoEls(gmeshK, top_X[icase], /*dim=*/2);
//             HRefineGeoSet_OneLevel(gmeshK, patch_cand);
//
//
//             TPZCompMesh* cmeshK = CreateCMesh(gmeshK, pOrder);
//
//             const REAL FSK = FindFS_Bisection(cmeshK, lo_fs, hi_fs, tol_rel, max_bis, verbose, loadid);
//             const STATE dFS = std::abs(FSK - FS0);
//
//             std::cout << "\n[ImpactStudy] ilevel="<<ilevel
//             << " case="<<icase
//             << " X=("<<top_X[icase][0]<<","<<top_X[icase][1]<<","<<top_X[icase][2]<<")"
//             << " FSK="<<FSK
//             << " dFS="<<dFS << std::endl;
//
//             dFSVec.push_back(dFS);
//
//             delete cmeshK;
//             delete gmeshK;
//         }
//
//         // campeão do nível
//         indextop = (int)std::distance(dFSVec.begin(),
//                                       std::max_element(dFSVec.begin(), dFSVec.end()));
//
//         cout << "indextop = "<< indextop << std::endl;
//         delete cmesh0;
//         delete gmesh0;
//     }
// }

// // ------------------------------------------------------------
// // Função principal do experimento
// // ------------------------------------------------------------
// void RunImpactStudy_FS2()
// {
//     using Clock = std::chrono::steady_clock;
//
//     const int pOrder    = 2;
//     const int ref0      = 1;
//     const int loadid    = 1;
//
//     // FS-bisection
//     const REAL lo_fs    = 0.5;
//     const REAL hi_fs    = 30.0;
//     const REAL tol_rel  = 0.015;
//     const int  max_bis  = 20;
//     const int  verbose  = 0;
//
//     // Solve interno (se você quiser usar)
//     const int nref_solve = 2;
//     const STATE tolref   = 0.015;
//
//     // --------------------------------------------------------
//     // 1) malha base + FS0
//     // --------------------------------------------------------
//     TPZGeoMesh* gmesh0 = TriGMesh(ref0);
//     TPZCompMesh* cmesh0 = CreateCMesh(gmesh0, pOrder);
//
//     auto t0 = Clock::now();
//
//     // Escolha 1: usar FindFS_Bisection diretamente (mais leve)
//     const REAL FS0 = FindFS_Bisection(cmesh0, lo_fs, hi_fs, tol_rel, max_bis, verbose, loadid);
//
//     // (se preferir usar seu Solve, troque a linha acima por:)
//     // const REAL FS0 = Solve(cmesh0, loadid, "base.vtk", nref_solve, tolref);
//
//     std::cout << "\n[ImpactStudy] FS0 (base) = " << FS0 << "\n";
//
//     // --------------------------------------------------------
//     // 2) escolher candidatos (top-M por indicador)
//     // --------------------------------------------------------
//     TPZVec<REAL> ind;
//     ComputeElementDeformation(cmesh0, ind);
//
//     const int M = 3;
//     auto top_cel = TopM_ByIndicator(ind, M);
//
//     REAL maxind = 0.0;
//     for (auto v : ind) maxind = std::max(maxind, v);
//     std::cout << "[ImpactStudy] max indicator = " << maxind << "\n";
//     std::cout << "[ImpactStudy] top_cel.size() = " << top_cel.size() << "\n";
//
//     // Coordenadas globais dos centros dos candidatos (na malha base)
//     auto top_X = MapTopCompElsToGeoCenterCoords(cmesh0, top_cel);
//     std::cout << "[ImpactStudy] candidatos (coords) = " << top_X.size() << "\n";
//
//     // --------------------------------------------------------
//     // 3) testar impacto por candidato
//     // --------------------------------------------------------
//     std::vector<ImpactResult> results;
//     results.reserve(top_X.size());
//
//     int icase = 0;
//     for (const auto& Xcand : top_X) {
//         icase++;
//
//         TPZGeoMesh* gmeshK = TriGMesh(ref0);
//
//         // patch 1-ring no elemento que contém Xcand
//         auto patch = CollectPatchGeoEls(gmeshK, Xcand, /*dim=*/2);
//         if (patch.empty()) {
//             delete gmeshK;
//             continue;
//         }
//
//         // 1 nível de refinamento no patch
//         HRefineGeoSet_OneLevel(gmeshK, patch);
//
//         // opcional: VTK do geo refinado
//         {
//             std::ofstream vtk(("impact_gmesh_" + std::to_string(icase) + ".vtk").c_str());
//             TPZVTKGeoMesh::PrintGMeshVTK(gmeshK, vtk, true);
//         }
//
//         TPZCompMesh* cmeshK = CreateCMesh(gmeshK, pOrder);
//
//         const REAL FSK = FindFS_Bisection(cmeshK, lo_fs, hi_fs, tol_rel, max_bis, verbose, loadid);
//         // ou: const REAL FSK = Solve(cmeshK, loadid, "impact_case_"+std::to_string(icase)+".vtk", nref_solve, tolref);
//
//         ImpactResult r;
//         r.icase = icase;
//         r.X = Xcand;
//         r.FS = FSK;
//         r.dFS = std::abs(FSK - FS0);
//         results.push_back(r);
//
//         std::cout << "\n [ImpactStudy] case " << icase
//         << "  X=(" << Xcand[0] << "," << Xcand[1] << "," << Xcand[2] << ")"
//         << "  FSK=" << FSK
//         << "  dFS=" << r.dFS << "\n";
//
//         delete cmeshK;
//         delete gmeshK;
//     }
//
//     if(results.empty()){
//         std::cout << "\n[ImpactStudy] Nenhum caso válido (patch vazio / FindElement falhou).\n";
//         delete cmesh0;
//         delete gmesh0;
//         return;
//     }
//
//     // --------------------------------------------------------
//     // 4) ranking final
//     // --------------------------------------------------------
//     std::sort(results.begin(), results.end(),
//               [](const ImpactResult& a, const ImpactResult& b){ return a.dFS > b.dFS; });
//
//     std::cout << "\n========== RANKING (maior impacto em FS) ==========\n";
//     std::cout << "rank  case   FS        |FS-FS0|      X\n";
//     for (int i = 0; i < (int)results.size(); i++) {
//         const auto& r = results[i];
//         std::cout << std::setw(4) << (i+1) << "  "
//         << std::setw(4) << r.icase << "  "
//         << std::setw(10) << std::setprecision(8) << r.FS << "  "
//         << std::setw(10) << std::setprecision(8) << r.dFS << "   "
//         << "(" << r.X[0] << "," << r.X[1] << "," << r.X[2] << ")\n";
//     }
//
//     // --------------------------------------------------------
//     // 5) pega o melhor e faz +1 nível de refinamento extra
//     // --------------------------------------------------------
//     const auto best = results.front();
//     std::cout << "\n[ImpactStudy] BEST case=" << best.icase
//     << "  dFS=" << best.dFS
//     << "  X=(" << best.X[0] << "," << best.X[1] << "," << best.X[2] << ")\n";
//
//     {
//         TPZGeoMesh* gmeshB = TriGMesh(ref0);
//
//         auto patchB = CollectPatchGeoEls(gmeshB, best.X, /*dim=*/2);
//         if(!patchB.empty()){
//             // primeiro nível (mesmo que no estudo)
//             HRefineGeoSet_OneLevel(gmeshB, patchB);
//             // segundo nível extra
//             HRefineGeoSet_OneLevel(gmeshB, patchB);
//
//             std::ofstream vtk(("impact_gmesh_" + std::to_string(best.icase) + ".vtk").c_str());
//             TPZVTKGeoMesh::PrintGMeshVTK(gmeshB, vtk, true);
//
//             TPZCompMesh* cmeshB = CreateCMesh(gmeshB, pOrder);
//             const REAL FSbest = FindFS_Bisection(cmeshB, lo_fs, hi_fs, tol_rel, max_bis, verbose, loadid);
//
//             std::cout << "[ImpactStudy] FS(best +1 level) = " << FSbest
//             << "  dFS=" << std::abs(FSbest - FS0) << "\n";
//
//             delete cmeshB;
//         } else {
//             std::cout << "[ImpactStudy] patch do BEST vazio no rerun.\n";
//         }
//         delete gmeshB;
//     }
//
//     auto t2 = Clock::now();
//     std::chrono::duration<double> secs = t2 - t0;
//     std::cout << "\n[ImpactStudy] total time = " << secs.count() << " s\n";
//
//     delete cmesh0;
//     delete gmesh0;
// }
