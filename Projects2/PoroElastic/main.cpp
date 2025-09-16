// ===== main.cpp =====
#include <pzgmesh.h>
#include <pzcmesh.h>
#include <pzmanvector.h>
#include <TPZBndCondT.h>
#include <TPZLinearAnalysis.h>
#include <pzskylstrmatrix.h>
#include <pzstepsolver.h>
#include <TPZVTKGeoMesh.h>
#include <pzbuildmultiphysicsmesh.h>
#include <pzinterpolationspace.h>

#include "Elasticity/TPZElasticity2D.h"
#include "Poisson/TPZMatPoisson.h"
#include "Elasticity/TPZMatPoroElastic2DMem.h"
#include "tpzgeoelrefpattern.h"
#include "pzmultiphysicselement.h"

#include <iostream>
#include <fstream>
#include <map>
#include <pzgraphmesh.h>
using std::cout; using std::endl;

class TPZSSpStructMatrix;
class TPZSSpStructMatrix;
class TPZSSpStructMatrix;
class TPZSSpStructMatrix;
// ---- parâmetros (equivalentes ao seu script) ----
static constexpr int   dim   = 2;
static constexpr int   order = 1;
static constexpr double tick = 1.0;   // (não usado aqui)
static constexpr double lh   = 0.1;
static constexpr double lv   = 1.0;

// materiais/BC ids
static constexpr int kMatVol   = 1;
static constexpr int kBC_Top   = -1;  // marker 1
static constexpr int kBC_Right = -2;  // marker 2
static constexpr int kBC_Left  = -3;  // marker 3
static constexpr int kBC_Bot   = -4;  // marker 4

// propriedades físicas
static constexpr double young = 3.e5;
static constexpr double nu    = 0.2;
static constexpr double alphaB= 1.0;
static constexpr double rhof  = 1000.0;
static constexpr double Se     = 1e-10;
static constexpr double mu     = 1e-2;
static constexpr double perm   = 1e-10;

// ===== helpers de debug (curtos) =====
static void PrintMeshSummary(const char* title, TPZCompMesh* cmesh) {
    cout << title
    << " nelem=" << cmesh->NElements()
    << " ncon="  << cmesh->NConnects()
    << " neq="   << cmesh->NEquations() << endl;
}
static void DumpGeoBCSummary(TPZGeoMesh* gmesh){
    const int gdim = gmesh->Dimension();
    std::map<int,int> byid; int nb=0;
    for (auto gel : gmesh->ElementVec()){
        if (!gel) continue;
        if (gel->Dimension()==gdim-1){ byid[gel->MaterialId()]++; nb++; }
    }
    cout << "[geo] boundary elems = " << nb << "\n";
    for (auto &kv: byid) cout << "  id " << kv.first << " : " << kv.second << "\n";
}

// ===== criação da malha geométrica: tirinha 1x10 de quads =====
static TPZGeoMesh* CreateGeoMeshStrip_Quad_1x10()
{
    auto *gmesh = new TPZGeoMesh();
    gmesh->SetDimension(dim);

    // nós: duas colunas (x=0 e x=0.1), y = 0:0.1:1
    TPZManVector<REAL,3> x(3,0.0);
    gmesh->NodeVec().Resize(22);
    int nid = 0;
    for (int i=0;i<=10;i++){ // coluna x=0
        x[0]=0.0; x[1]=0.1*i; gmesh->NodeVec()[nid] = TPZGeoNode(nid,x,*gmesh); nid++;
    }
    for (int i=0;i<=10;i++){ // coluna x=lh
        x[0]=lh; x[1]=0.1*i; gmesh->NodeVec()[nid] = TPZGeoNode(nid,x,*gmesh); nid++;
    }

    // quads (10 elementos) - conectividade 0-based: {i, 11+i, 12+i, i+1}
    for (int i=0;i<10;i++){
        TPZManVector<int64_t,4> nodes(4);
        nodes[0] = i;
        nodes[1] = 11+i;
        nodes[2] = 12+i;
        nodes[3] = i+1;
        int64_t elid;
        gmesh->CreateGeoElement(EQuadrilateral, nodes, kMatVol, elid);
    }

    // linhas de contorno
    // esquerda (x=0): 10 segmentos
    for (int i=0;i<10;i++){
        TPZManVector<int64_t,2> n = { i, i+1 };
        new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(n, kBC_Left, *gmesh);
    }
    // direita (x=lh)
    for (int i=11;i<21;i++){
        TPZManVector<int64_t,2> n = { i, i+1 };
        new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(n, kBC_Right, *gmesh);
    }
    // base (y=0)
    { TPZManVector<int64_t,2> n = { 0, 11 };
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(n, kBC_Bot, *gmesh); }
    // topo (y=lv)
    { TPZManVector<int64_t,2> n = { 10, 21 };
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(n, kBC_Top, *gmesh); }

    gmesh->BuildConnectivity();

    std::ofstream vtk("geom.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, vtk, false);
    return gmesh;
}

// ===== submalha H1(u) com BCs "dummy" (para criar comp-BCs) =====
static TPZCompMesh* CMeshElastic(TPZGeoMesh* gmesh, int pOrder)
{
    auto *cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(dim);
    cmesh->SetDefaultOrder(pOrder);
    cmesh->SetAllCreateFunctionsContinuous();

    auto *mat = new TPZElasticity2D(kMatVol, young,nu, 0., 0., 1); // E,nu não importam aqui
    cmesh->InsertMaterialObject(mat);

    // inserir materiais de BC na SUBMALHA (só para AutoBuild criar comp-BCs)
    TPZFMatrix<STATE> v10(2,2,0.); TPZManVector<STATE,2> v20(2,0.);

    cmesh->InsertMaterialObject(mat->CreateBC(mat, kBC_Top,   0, v10, v20)); // Neumann "dummy"
    cmesh->InsertMaterialObject(mat->CreateBC(mat, kBC_Bot,   0, v10, v20)); // Neumann "dummy"
    cmesh->InsertMaterialObject(mat->CreateBC(mat, kBC_Right,  0, v10, v20)); // Neumann "dummy"
    cmesh->InsertMaterialObject(mat->CreateBC(mat, kBC_Left,   0, v10, v20)); // Neumann "dummy"


    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    return cmesh;
}

// ===== submalha H1(p) com BCs "dummy" =====
static TPZCompMesh* CMeshPressure(TPZGeoMesh* gmesh, int pOrder)
{
    auto *cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(dim);
    cmesh->SetDefaultOrder(pOrder);
    cmesh->SetAllCreateFunctionsContinuous();

    auto *mat = new TPZMatPoisson<REAL>(kMatVol, dim);
    cmesh->InsertMaterialObject(mat);

    TPZFMatrix<STATE> v10(1,1,0.); TPZManVector<STATE,1> v20(1,0.);
     cmesh->InsertMaterialObject(mat->CreateBC(mat, kBC_Top,   0, v10, v20));
     cmesh->InsertMaterialObject(mat->CreateBC(mat, kBC_Right,0, v10, v20));
    cmesh->InsertMaterialObject(mat->CreateBC(mat, kBC_Left,  0, v10, v20));
     cmesh->InsertMaterialObject(mat->CreateBC(mat, kBC_Bot,   0, v10, v20));

    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    return cmesh;
}

// ===== cria a mista, insere material poroelástico e BCs "reais" =====
static TPZCompMesh* CreateMPhysWithMaterialsAndBCs(TPZGeoMesh* gmesh)
{
    auto *mphys = new TPZCompMesh(gmesh);
    mphys->SetDimModel(dim);
    mphys->SetAllCreateFunctionsMultiphysicElemWithMem();

    // material u–p (TPZMatPoroElastic2DMem)
    auto *mat = new TPZMatPoroElastic2DMem<TPZElasticMem>(kMatVol);
    mat->SetId(kMatVol);
    mat->SetUpdateMem(true);

    // propriedades
    mat->SetElasticity(young, nu);
    TPZElasticResponse ER;
    ER.SetEngineeringData(young, nu);
    mat->SetElasticResponse(ER);
    mat->SetAlpha(alphaB);
    mat->SetSe(Se);
    mat->SetPermeability(perm);
    mat->SetViscosity(mu);
    mat->SetRhoF(rhof);
    mat->SetGravity(0.0, -9.81);     // <<< gravidade NÃO zero
    mat->SetBodyForce(0.0, 0.0);
    mphys->InsertMaterialObject(mat);

    TPZFMatrix<STATE> v1(3,3,0.);
    TPZManVector<STATE,3> v2(3,0.);
    // Val2: máscara
    v2[0] = 0.;  v2[1] = 1.;
    auto *bcB = mat->CreateBC(mat, kBC_Bot, 4, v1, v2);
    mphys->InsertMaterialObject(bcB);

    v2[0] = 1.;  v2[1] = 0.;
    auto *bcL = mat->CreateBC(mat, kBC_Left, 4, v1, v2);
    mphys->InsertMaterialObject(bcL);

    v2[0] = 1.;  v2[1] = 0.;                 // (mx,my)=(1,0)
    auto *bcR = mat->CreateBC(mat, kBC_Right, 4, v1, v2); // type=4
    mphys->InsertMaterialObject(bcR);

    v2[0] = 0.0;                 // t_x
    v2[1] = 1000.;             // t_y = -σ0  (compressão; ajuste o sinal se seu eixo y/z for o oposto)
    v2[2] = 0.0;                // pD = 0 (drenado)
    auto *bcT = mat->CreateBC(mat, kBC_Top, 10, v1, v2);  // type=10
    mphys->InsertMaterialObject(bcT);

    mphys->AutoBuild();
    mphys->AdjustBoundaryElements();
    mphys->CleanUpUnconnectedNodes();
    return mphys;
}

// ===== inicializa memória elástica (E,nu) nos pontos de integração =====
static void InitializeMemory(TPZCompMesh* mphys, REAL E, REAL poisson)
{
    auto* pMatWithMem =
    dynamic_cast<TPZMatWithMem<TPZElasticMem>*>(mphys->MaterialVec()[kMatVol]);
    if (!pMatWithMem) return;

    mphys->LoadReferences();
    pMatWithMem->SetUpdateMem(true);
    TPZAdmChunkVector<TPZElasticMem> &mem = *pMatWithMem->GetMemory();

    const int nels = mphys->NElements();
    for (int iel=0; iel<nels; iel++) {
        auto *cel = mphys->ElementVec()[iel];
        auto *mpel = dynamic_cast<TPZMultiphysicsElement*>(cel);
        if (!mpel || !cel->Material() || cel->Material()->Id()!=kMatVol) continue;

        const TPZIntPoints& intpoints = mpel->GetIntegrationRule();
        const int nint = intpoints.NPoints();

        const int nsubs = (int)mpel->ElementVec().size();
        TPZVec<TPZMaterialDataT<STATE>> data; data.Resize(nsubs);
        TPZVec<TPZTransform<>> trvec;         trvec.Resize(nsubs);
        mpel->InitMaterialData(data);
        mpel->AffineTransform(trvec);

        TPZManVector<REAL,3> q(2,0.0);
        for (int ip=0; ip<nint; ip++) {
            REAL w; intpoints.Point(ip, q, w);
            for (int iv=0; iv<nsubs; ++iv) data[iv].intLocPtIndex = ip;
            mpel->ComputeRequiredData(q, trvec, data);
            const int64_t idx = data[0].intGlobPtIndex;
            if (idx >= mem.NElements()) mem.Resize(idx + 1);
            mem[idx].m_ER.SetEngineeringData(E, poisson);
        }
    }
    pMatWithMem->SetUpdateMem(false);
}

// ---- PoroLogger.hpp (ou num .cpp seu) -----------------------
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>

struct PoroLogger {
    std::ofstream csv;
    std::string   path;
    const char*   logo = "⟦ PZ·PORO ⟧";

    void OpenCSV(const std::string& p) {
        path = p;
        csv.open(path, std::ios::out | std::ios::trunc);
        // cabeçalho
        csv << "step,time,rhs_norm,rhs_u_norm,rhs_p_norm,"
        "delta_norm,delta_u_norm,delta_p_norm,"
        "p_min,p_mean,p_max,uy_top\n";
        csv.flush();
    }

    // imprime e escreve CSV em uma chamada
    void PrintAndCSV(int step, double t,
                     double rhs_norm, double rhsu_norm, double rhsp_norm,
                     double delta_norm, double deltau_norm, double deltap_norm,
                     double pmin, double pmean, double pmax,
                     double uy_top = 0.0)
    {
        std::cout << std::fixed << std::setprecision(6)
        << logo << " step=" << step
        << "  t=" << t
        << "  ||rhs||=" << rhs_norm
        << " (u=" << rhsu_norm << ", p=" << rhsp_norm << ")"
        << "  ||Δ||=" << delta_norm
        << "  ||Δu||=" << deltau_norm
        << "  ||Δp||=" << deltap_norm
        << "  p[min,mean,max]=[" << pmin << "," << pmean << "," << pmax << "]"
        << "  uy_top=" << uy_top
        << "\n";

        if (csv.is_open()) {
            csv << step << ','
            << std::setprecision(16) << t << ','
            << rhs_norm << ',' << rhsu_norm << ',' << rhsp_norm << ','
            << delta_norm << ',' << deltau_norm << ',' << deltap_norm << ','
            << pmin << ',' << pmean << ',' << pmax << ','
            << uy_top << '\n';
            csv.flush();
        }
    }
};

// utilitário para separar a norma do RHS entre u e p (assumindo ordenação [u|p])
inline void SplitRHS_UP(const TPZFMatrix<STATE>& rhs, int nequ_u,
                        double& rhsu_norm, double& rhsp_norm)
{
    double su=0., sp=0.;
    const int neq = rhs.Rows();
    for (int i=0; i<nequ_u; ++i) { const double v = rhs(i,0); su += v*v; }
    for (int i=nequ_u; i<neq;  ++i) { const double v = rhs(i,0); sp += v*v; }
    rhsu_norm = std::sqrt(su);
    rhsp_norm = std::sqrt(sp);
}


TPZVec<REAL>  UyAtPoint(TPZCompMesh* cmesh, REAL x, REAL y,int var)
{
    TPZGeoMesh* gmesh = cmesh->Reference();
    TPZManVector<REAL,3> X(3,0.); X[0]=x; X[1]=y;
    TPZManVector<REAL,3> qsi(3,0.);
    int64_t elid=0;

    cmesh->LoadReferences();          // liga cada TPZGeoEl ao seu TPZCompEl
    TPZGeoEl* gel = gmesh->FindElement(X, qsi, elid, /*dim=*/2);
    if(!gel || !gel->Reference()) DebugStop();

    auto* cel = dynamic_cast<TPZMultiphysicsElement*>(gel->Reference());
    if(!cel) DebugStop();

    TPZVec<REAL> sol;
    cel->Solution(qsi, var, sol);
    return sol; // componente y
}
#include "pzcreateapproxspace.h"  // se precisar do ApproxSpace()

static void AddBoundaryCompElsToMF(TPZCompMesh* mphys)
{
    TPZGeoMesh* gmesh = mphys->Reference();
    const int gdim = gmesh->Dimension();
    int64_t index = -1;

    // garanta que a mista sabe criar elementos multifísicos (com memória)
    mphys->SetAllCreateFunctionsMultiphysicElemWithMem();

    for (auto gel : gmesh->ElementVec()) {
        if (!gel) continue;
        if (gel->Dimension() != gdim-1) continue;       // só borda
        const int bcid = gel->MaterialId();
        if (!mphys->FindMaterial(bcid)) continue;       // só IDs com material BC inserido na mista

        // ✅ uma das duas linhas abaixo (conforme sua versão):
        mphys->CreateCompEl(gel);
         //mphys->ApproxSpace().CreateCompEl(gel, *mphys);
    }
}

#include "pzfstrmatrix.h"
#include "TPZBSpStructMatrix.h"
int main()
{
    // --- malhas ---
    TPZGeoMesh        *gmesh = CreateGeoMeshStrip_Quad_1x10();
    TPZCompMesh *cmesh_u    = CMeshElastic (gmesh, order);
    TPZCompMesh *cmesh_p    = CMeshPressure(gmesh, order);

    // mista com material e BCs reais
    TPZCompMesh *mphys = CreateMPhysWithMaterialsAndBCs(gmesh);

    // gerar elementos/ligação a partir das submalhas
    TPZVec<TPZCompMesh*> meshvec(2);
    meshvec[0]=cmesh_u;
    meshvec[1]=cmesh_p;
    TPZBuildMultiphysicsMesh::AddElements(meshvec, mphys);
    TPZBuildMultiphysicsMesh::AddConnects(meshvec, mphys);
    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, mphys);
    mphys->LoadReferences();

    // memória elástica
    InitializeMemory(mphys, young, nu);


    // pegue o ponteiro do material MF (precisa ser o seu poroelástico)
    auto *mat = dynamic_cast<TPZMatPoroElastic2DMem<TPZElasticMem>*>(
        mphys->MaterialVec()[kMatVol]);

    TPZLinearAnalysis an(mphys);
    TPZSkylineStructMatrix<REAL> str(mphys);
    an.SetStructuralMatrix(str);
    TPZStepSolver<REAL> direct; direct.SetDirect(ELDLt);
    an.SetSolver(direct);


    // solução total corrente
    TPZFMatrix<STATE> xtotal = mphys->Solution();
    xtotal.Zero();
    // --- antes do loop: CSV + “prev” para medir Δu, Δp ---
    std::ofstream csv("poroelog.csv");
    csv << "step,time,rhs_norm,rhs_u_norm,rhs_p_norm,"
    "delta_norm,delta_u_norm,delta_p_norm,"
    "p_min,p_mean,p_max\n";

    // projeta solução inicial para as submalhas
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphys);
    TPZFMatrix<STATE> u_prev = meshvec[0]->Solution();
    TPZFMatrix<STATE> p_prev = meshvec[1]->Solution();


    // ... antes do loop
    std::ofstream prof("poro_profile.csv");
    prof << "step,time,y,uy,p\n";

    int nsteps=5;
    double t = 0.0;
    REAL   dt = 1.;           // se quiser variável, atualize *depois* do passo.
    mat->SetTimeStep(dt);

    mat->UseMassScaledByInvDt(true); // incremental (S/Δt + H)

    for (int k=1; k<=nsteps; ++k) {

        mat->SetTimeStep(dt);
        // 1) monta e resolve pelo incremento Δx
        an.Assemble();
        an.Solve();
        TPZFMatrix<STATE> delta = an.Solution();

        // 2) acumula e carrega a solução total no multifísico
        xtotal += delta;
        mphys->LoadSolution(xtotal);

        // 3) **agora** projeta para submalhas e mede u^{n+1}, p^{n+1}
        TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphys);
        TPZFMatrix<STATE> u_now = meshvec[0]->Solution();
        TPZFMatrix<STATE> p_now = meshvec[1]->Solution();

        TPZFMatrix<STATE> du = u_now; du -= u_prev;
        TPZFMatrix<STATE> dp = p_now; dp -= p_prev;

        std::cout << std::fixed << std::setprecision(12)
        << "[PZ·PORO] step=" << k
        << "  t=" << t
        << "  dt=" << dt
        << "  ||Δ||="  << Norm(delta)
        << "  ||Δu||=" << Norm(du)
        << "  ||Δp||=" << Norm(dp) << "\n";

        // 4) perfil (x fixo, varrendo y) — grava antes de incrementar y
        REAL deltay = lv/10.;  // garanta que lv está definido (altura)
        REAL y = 0.;
        for (int idiv=0; idiv<=10; ++idiv) { // <= para incluir topo
            TPZVec<REAL> uy   = UyAtPoint(mphys, lh*0.5, y, 1);  // [ux, uy]
            TPZVec<REAL> pres = UyAtPoint(mphys, lh*0.5, y, 2);  // use função própria p/ pressão se tiver
            prof << k << "," << std::setprecision(6) << t << ","
            << y << "," << uy[1] << "," << pres[0] << "\n";
            y += deltay;
        }

        // 5) prepara próximo passo
        u_prev = std::move(u_now);
        p_prev = std::move(p_now);

        //dt = pow(1.1, k) - t;  // (rever se faz sentido p/ seu caso)
        t += dt;
        // se quiser dt variável, atualize aqui e só então aplique no material:

    }
    prof.flush();
    prof.close();

    return 0;

}
