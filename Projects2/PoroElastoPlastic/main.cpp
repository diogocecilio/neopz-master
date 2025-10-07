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
#include "DarcyFlow/TPZDarcyFlow.h"
#include "pzlog.h"
#include <log4cxx/basicconfigurator.h>
#include "Poisson/TPZMatPoisson.h"
#include "Elasticity/TPZMatPoroElastic2DMem.h"
#include "tpzgeoelrefpattern.h"
#include "pzmultiphysicselement.h"
#include "Elasticity/TPZMatElastic2DMem.h"
#include <iostream>
#include <fstream>
#include <map>
#include <pzgraphmesh.h>
#include <TPZSSpStructMatrix.h> //symmetric sparse matrix storage
#include <pzskylstrmatrix.h> //symmetric skyline matrix storage
#include "pzsfstrmatrix.h"
#include "TPZPardisoSolver.h"
#include <Projection/TPZL2Projection.h>
#include <TPZSSpStructMatrix.h>
#include "Plasticity/TPZMatPoroElastoPlastic2DMem.h"
//#include "TPZElasticResponse.h"
#include "Plasticity/TPZElastoPlasticMem.h"
#include "Plasticity/TPZPlasticStepPV.h"
#include "Plasticity/TPZYCMohrCoulombPV.h"
#include "Plasticity/TPZMatElastoPlastic2D.h"
#include "TPZMaterialDataT.h"
#include "Plasticity/pzelastoplasticanalysis.h"
#include "Plasticity/TPZElasticResponse.h"
#include "pzcompelwithmem.h"
#include "pzmultiphysicscompel.h"
#include "pzgeoquad.h"          // você tem GeoQuad no print
// #include "TPZGeoTriangle.h"  // se tiver triângulos
// #include "TPZGeoLinear.h"    // se tiver elementos 1D

using std::cout; using std::endl;



static constexpr int dim   =   2;
static constexpr int order   =   2;

static constexpr int kMatVol   =   1;
static constexpr int kBC_Top   =  -1;  // marker 1

static constexpr int kBC_Right =  -2;  // marker 2
static constexpr int kBC_Left  =  -3;  // marker 3
static constexpr int kBC_Bot   =  -4;  // marker 4
static constexpr int kBC_Bot2   =  -8;  // marker 4
static constexpr int kBC_Top2   = -5;  // marker 1
static constexpr int kBC_nodeleft   = -6;  // marker 1
static constexpr int kBC_noderigth   = -7;  // marker 1

static constexpr int bcidbottom = -1; // bottom
static constexpr int bcidright = -2; // right
static constexpr int bcidtopright = -3; // topright
static constexpr int bcidtopleft = -4; // top left
static constexpr int bcidleft = -5; // left
static constexpr int bcidramp = -6; // ramp
#ifdef PZ_LOG
// escolha um nome de categoria claro (use pontos para hierarquia)
static TPZLogger logger_poro_main("PoroPlastic_Main");
#endif


typedef TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse> TPlasticMC;
typedef TPZMatPoroElastoPlastic2DMem<TPlasticMC,TPZElastoPlasticMem> poroplasticmat;

// Guarde isto num header, p.ex. PoroMechParams.h
// Requer: typedef STATE já definido no seu projeto.

struct PoroMechParams {
    // --- parâmetros mecânicos ---
    STATE alpha   = 1.;      // [-]
    STATE young   = 20000;    // Pa
    STATE nu      = 0.49;      // [-]

    // derivados (preenchidos em Update())
    STATE G       = 0.0;      // Pa  (módulo de cisalhamento)
    STATE lambdaL = 0.0;      // Pa  (1º módulo de Lamé)

    // resistência/atrito (se usados)
    STATE cohesion = 10;   // Pa
    STATE phi      = 30.0 * M_PI / 180.0; // rad
    STATE psi      = 30.0 * M_PI / 180.0; // rad

    // pressão inicial
    STATE p0 = 0.0;           // Pa

    // --- fluxo ---
    STATE perm = 1.e-7;     // m^2
    STATE mu   = 1.0e-3;      // Pa·s
    STATE Se   = 0.0;         // 1/Pa (compressibilidade do fluido)
    STATE k_over_mu = 0.0;    // (m^2 / Pa·s) = perm / mu (preenchido em Update())

    // --- carregamentos volumétricos ---
    STATE fx =  0.0;          // N/m^3
    STATE fy = -20;        // N/m^3

    // --- densidade do fluido (se usado) ---
    STATE rhof = 0.0;         // kg/m^3

    // recalcula os derivados quando mudar E, nu, perm, mu
    void Update() {
        G        = young / (2.0 * (1.0 + nu));
        lambdaL  = young * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
        k_over_mu = perm / mu;
    }
};


static void InitializeMemory(TPZCompMesh* mphys)
{
    if (!mphys) return;

    auto* mat = dynamic_cast<poroplasticmat*>(mphys->FindMaterial(kMatVol));
    if (!mat) return;

    mphys->LoadReferences();


    mat->SetUpdateMem(true);

    auto memSP = mat->GetMemory(); // std::shared_ptr<TPZAdmChunkVector<TMEM>>
    // #ifdef PZ_LOG
    // if (logger_poro.isDebugEnabled()) {
    //     std::ostringstream oss;
    //     oss << "[InitializeMemory] BEGIN\n"
    //     << "CompMesh NElements = " << mphys->NElements() << "\n"
    //     << "Mem NElements (antes) = " << (memSP ? memSP->NElements() : -1) << "\n";
    //     LOGPZ_DEBUG(logger_poro, oss.str());
    // }
    // #endif

    const int nel = mphys->NElements();

    for (int iel = 0; iel < nel; ++iel)
    {
        TPZCompEl* cel = mphys->ElementVec()[iel];
        if (!cel) continue;


        auto* mpel = dynamic_cast<TPZMultiphysicsElement*>(cel);
        if (!mpel) continue;
        if (!cel->Material()) continue;
        if (cel->Material()->Id() != kMatVol) continue;

        const int nsubs = (int)mpel->ElementVec().size();
        if (nsubs == 0) continue;


        TPZVec<TPZMaterialDataT<STATE>> datavec(nsubs);
        TPZVec<TPZTransform<STATE>>     tr(nsubs);
        mpel->InitMaterialData(datavec);
        mpel->AffineTransform(tr);

        const TPZIntPoints& ir = mpel->GetIntegrationRule();
        const int nint = ir.NPoints();

        TPZManVector<REAL,3> q(dim, 0.0);

        for (int ip = 0; ip < nint; ++ip)
        {
            REAL w = 0.0;
            ir.Point(ip, q, w);

            mpel->ComputeRequiredData(q, tr, datavec);

            const int idx = mat->PushMemItem();

            auto& mem = (*memSP)[idx];

            PoroMechParams prm;
            mem.m_elastoplastic_state.fmatpropinit.Resize(3);
            mem.m_elastoplastic_state.fmatprop.Resize(3);
            mem.m_elastoplastic_state.fmatpropinit[0] = prm.cohesion; // unidades consistentes
            mem.m_elastoplastic_state.fmatpropinit[1] = prm.phi;      // rad
            mem.m_elastoplastic_state.fmatpropinit[2] = prm.psi;      // rad
            mem.m_elastoplastic_state.fmatprop = mem.m_elastoplastic_state.fmatpropinit;

            mem.m_elastoplastic_state.m_eps_t.Zero();
            mem.m_elastoplastic_state.m_eps_p.Zero();
            mem.m_sigma.Zero();

            mem.m_u.Resize(3, 0.0);
            mem.m_plastic_steps = 0;
            mem.m_phi = 0.0;

            mem.m_ER.SetEngineeringData(prm.young, prm.nu);

            mem.m_elastoplastic_state.fpressure = prm.p0;
            mem.m_elastoplastic_state.fdPorePressure.Resize(dim, 0.0);
            mem.m_elastoplastic_state.fSolU.Resize(dim, 0.0);
            mem.m_elastoplastic_state.fGradSolU.Zero();


            // #ifdef PZ_LOG
            // if (logger_poro.isDebugEnabled()) {
            //     std::ostringstream oss;
            //     oss << "[InitializeMemory] el=" << iel
            //     << " ip=" << ip
            //     << " idx=" << idx
            //     << " q=(" << q[0] << "," << q[1] << ") w=" << w << "\n";
            //     LOGPZ_DEBUG(logger_poro, oss.str());
            // }
            // #endif
        } // ip
    } // iel

    // #ifdef PZ_LOG
    // if (logger_poro.isDebugEnabled()) {
    //     std::ostringstream oss;
    //     oss << "Mem NElements (depois) = " << (memSP ? memSP->NElements() : -1) << "\n"
    //     << "[InitializeMemory] END";
    //     LOGPZ_DEBUG(logger_poro, oss.str());
    // }
    // #endif

    mat->SetUpdateMem(false);
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
            if (std::fabs(y0-0) < tol && std::fabs(yf-0) < tol)           bcid = bcidbottom; // bottom
            else if (std::fabs(x0-L)<tol && std::fabs(xf-L)<tol)          bcid = bcidright; // right
            else if (std::fabs(y0-h1)<tol && std::fabs(yf-h1)<tol)        bcid = bcidtopright; // top right
            else if (std::fabs(y0-(h1+h2))<tol && std::fabs(yf-(h1+h2))<tol) bcid = bcidtopleft; // top left
            else if (std::fabs(x0-0)<tol && std::fabs(xf-0)<tol)          bcid = bcidleft; // left
            else if (std::fabs(xf-x0)>tol && std::fabs(yf-y0)>tol)        bcid = bcidramp; // ramp
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

TPZGeoMesh* CreateSingleQuadMesh()
{

    REAL co[4][2] = {{0.,0.},{1,0.},{1.,1.},{0.,1.}};
    long indices[1][4] = {{0,1,2,3}};
    TPZGeoEl *elvec[1];
    TPZGeoMesh *gmesh = new TPZGeoMesh();
    gmesh->SetDimension ( 2 );
    long nnode = 4;
    long nod;
    for ( nod=0; nod<nnode; nod++ )
    {
        long nodind = gmesh->NodeVec().AllocateNewElement();
        TPZVec<REAL> coord ( 2 );
        coord[0] = co[nod][0];
        coord[1] = co[nod][1];
        gmesh->NodeVec() [nodind] = TPZGeoNode ( nod,coord,*gmesh );
    }

    long el;
    long nelem = 1;
    for ( el=0; el<nelem; el++ )
    {
        TPZVec<long> nodind ( 4 );
        for ( nod=0; nod<4; nod++ ) nodind[nod]=indices[el][nod];
        //    elvec[el] = new TPZGeoElQ2d(el,nodind,1);
        long index;
        elvec[el] = gmesh->CreateGeoElement ( EQuadrilateral,nodind,1,index );
    }

    TPZVec <long> TopoLine ( 2 );

    TopoLine[0] = 0;
    TopoLine[1] = 1;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( 1, TopoLine, kBC_Bot, *gmesh );

    TopoLine[0] = 1;
    TopoLine[1] = 2;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( 2, TopoLine, kBC_Right, *gmesh );

    TopoLine[0] = 2;
    TopoLine[1] = 3;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( 3, TopoLine, kBC_Top, *gmesh );

    TopoLine[0] = 3;
    TopoLine[1] = 0;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( 4, TopoLine, kBC_Left, *gmesh );

    TopoLine[0] = 2;
    TopoLine[1] = 3;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( 5, TopoLine, kBC_Top2, *gmesh );

    TPZVec <long> node ( 1 );
    node[0]=0;
    new TPZGeoElRefPattern< pzgeom::TPZGeoPoint> ( 6, node, kBC_nodeleft, *gmesh );//bottomrigth node
    gmesh->BuildConnectivity();

    node[0]=1;
    new TPZGeoElRefPattern< pzgeom::TPZGeoPoint> ( 7, node, kBC_noderigth, *gmesh );//bottomrigth node


    TopoLine[0] = 0;
    TopoLine[1] = 1;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( 8, TopoLine, kBC_Bot2, *gmesh );
    gmesh->BuildConnectivity();

    cout << "c" << endl;
    for ( int d = 0; d<2; d++ )
    {
        int nel = gmesh->NElements();
        TPZManVector<TPZGeoEl *> subels;
        for ( int iel = 0; iel<nel; iel++ )
        {
            TPZGeoEl *gel = gmesh->ElementVec() [iel];
            gel->Divide ( subels );
        }
    }
    // gmesh->BuildConnectivity();
    std::ofstream files ( "teste-mesh.vtk" );
    TPZVTKGeoMesh::PrintGMeshVTK ( gmesh,files,false );
    cout << "d" << endl;
    return gmesh;

}


static TPZCompMesh* CMeshElastic(TPZGeoMesh* gmesh){

    auto *cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(2);
    cmesh->SetDefaultOrder(order);
    cmesh->SetAllCreateFunctionsContinuousWithMem();

    auto *mat = new TPZMatElastic2DMem<TPZElasticMem>(1);
    mat->SetId(1);      // (mantido)

    PoroMechParams prm;
    TPZElasticResponse ER; ER.SetEngineeringData(prm.young, prm.nu);
    mat->SetUpdateMem(true);
    mat->SetElasticResponse(ER);
    mat->SetUpdateMem(false);

    cmesh->InsertMaterialObject(mat);
    TPZFMatrix<STATE> val1(2,2,0.);
    TPZVec<REAL> val2(2,0.);
    int dirichlet=0;
    // cmesh->InsertMaterialObject(mat->CreateBC(mat, kBC_Bot, dirichlet, val1, val2));
    // cmesh->InsertMaterialObject(mat->CreateBC(mat, kBC_Bot2, dirichlet, val1, val2));
    // cmesh->InsertMaterialObject(mat->CreateBC(mat, kBC_Top, dirichlet, val1, val2));
    // cmesh->InsertMaterialObject(mat->CreateBC(mat, kBC_Left, dirichlet, val1, val2));
    // cmesh->InsertMaterialObject(mat->CreateBC(mat, kBC_Right, dirichlet, val1, val2));
    // cmesh->InsertMaterialObject(mat->CreateBC(mat, kBC_Top2, dirichlet, val1, val2));
    // cmesh->InsertMaterialObject(mat->CreateBC(mat, kBC_nodeleft, dirichlet, val1, val2));
    // cmesh->InsertMaterialObject(mat->CreateBC(mat, kBC_noderigth, dirichlet, val1, val2));
    cmesh->InsertMaterialObject(mat->CreateBC(mat, bcidbottom, dirichlet, val1, val2));
    cmesh->InsertMaterialObject(mat->CreateBC(mat, bcidright, dirichlet, val1, val2));
    cmesh->InsertMaterialObject(mat->CreateBC(mat, bcidtopright, dirichlet, val1, val2));
    cmesh->InsertMaterialObject(mat->CreateBC(mat, bcidtopleft, dirichlet, val1, val2));
    cmesh->InsertMaterialObject(mat->CreateBC(mat, bcidleft, dirichlet, val1, val2));
    cmesh->InsertMaterialObject(mat->CreateBC(mat, bcidramp, dirichlet, val1, val2));
    cmesh->AutoBuild();
    return cmesh;
}

static TPZCompMesh* CMeshPressure(TPZGeoMesh* gmesh)
{
    auto *cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(dim);
    cmesh->SetDefaultOrder(order);
    cmesh->SetAllCreateFunctionsContinuous();

    auto *mat = new TPZDarcyFlow(kMatVol, dim);
    cmesh->InsertMaterialObject(mat);

    TPZFMatrix<STATE> val1(2,2,0.);
    TPZVec<REAL> val2(2,0.);
    int dirichlet=0;
    cmesh->InsertMaterialObject(mat->CreateBC(mat, bcidbottom, dirichlet, val1, val2));
    cmesh->InsertMaterialObject(mat->CreateBC(mat, bcidright, dirichlet, val1, val2));
    cmesh->InsertMaterialObject(mat->CreateBC(mat, bcidtopright, dirichlet, val1, val2));
    cmesh->InsertMaterialObject(mat->CreateBC(mat, bcidtopleft, dirichlet, val1, val2));
    cmesh->InsertMaterialObject(mat->CreateBC(mat, bcidleft, dirichlet, val1, val2));
    cmesh->InsertMaterialObject(mat->CreateBC(mat, bcidramp, dirichlet, val1, val2));
    // cmesh->InsertMaterialObject(mat->CreateBC(mat, kBC_Bot, dirichlet, val1, val2));
    // cmesh->InsertMaterialObject(mat->CreateBC(mat, kBC_Bot2, dirichlet, val1, val2));
    // cmesh->InsertMaterialObject(mat->CreateBC(mat, kBC_Top, dirichlet, val1, val2));
    // cmesh->InsertMaterialObject(mat->CreateBC(mat, kBC_Left, dirichlet, val1, val2));
    // cmesh->InsertMaterialObject(mat->CreateBC(mat, kBC_Right, dirichlet, val1, val2));
    // cmesh->InsertMaterialObject(mat->CreateBC(mat, kBC_Top2, dirichlet, val1, val2));
    // cmesh->InsertMaterialObject(mat->CreateBC(mat, kBC_nodeleft, dirichlet, val1, val2));
    // cmesh->InsertMaterialObject(mat->CreateBC(mat, kBC_noderigth, dirichlet, val1, val2));
    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    return cmesh;
}

#include "pzfstrmatrix.h"
#include <TPZCompMeshTools.h>
// ===== cria a mista, insere material poroelástico e BCs "reais" =====
static TPZCompMesh* CreateMPhysWithMaterialsAndBCs(TPZGeoMesh* gmesh)
{
    auto *mphys = new TPZCompMesh(gmesh);
    mphys->SetDimModel(dim);
    mphys->SetAllCreateFunctionsMultiphysicElemWithMem();

    auto mat = new poroplasticmat();
    mat->SetId(kMatVol);
    mat->SetUpdateMem(true);
    PoroMechParams prm;
    TPZElasticResponse ER;

    mat->SetElasticity(prm.young, prm.nu);
    mat->SetElasticResponse(ER);
    mat->SetAlpha(prm.alpha);
    mat->SetSe(prm.Se);
    mat->SetPermeability(prm.perm);
    mat->SetViscosity(prm.mu);
    mat->SetRhoF(prm.rhof);
    //mat->SetBodyForce(0.0, 0.0);
    mat->SetTimeStep(1.e-12);

    TPlasticMC mc;
    mc.fYC.SetUp(prm.phi, prm.psi, prm.cohesion, ER);
    mc.fER = ER;
    mc.SetStrengthReductionFactor(1.0);

    mat->SetPlasticModel(mc);




    mphys->InsertMaterialObject(mat);

    TPZFMatrix<STATE> v1(3,3,0.);
    TPZManVector<STATE,3> v2(3,0.);
    int dirdirichlet=3,pressure=2,newmannforca=1;


    // v2[0] = 1.0;
    // v2[1] = 0.0;
    // mphys->InsertMaterialObject(mat->CreateBC(mat, kBC_Left, dirdirichlet, v1, v2));//Direcional em x
    // v2[0] = 1.0;
    // v2[1] = 0.0;
    // mphys->InsertMaterialObject(mat->CreateBC(mat, kBC_Right, dirdirichlet, v1, v2));//Direcional em x
    //
    // v2[0] = 0.0;
    // v2[1] = 1.0;
    // mphys->InsertMaterialObject(mat->CreateBC(mat, kBC_Top, dirdirichlet, v1, v2));//base presa em y
    //
    //
    // v2[0] = 0.0;
    // v2[1] = 0.0;
    // v2[2] = 0.0;//pressao
    // mphys->InsertMaterialObject(mat->CreateBC(mat, kBC_Bot2, pressure, v1, v2));//pressao zero
    //
    // v2[0] = 0.0;
    // v2[1] = 1000.;
    // mphys->InsertMaterialObject(mat->CreateBC(mat, kBC_Bot, newmannforca, v1, v2));//tensao no top

    v2[0]=0.;
    v2[1]=1.;
    mphys->InsertMaterialObject(mat->CreateBC(mat, bcidbottom, dirdirichlet, v1, v2));
    v2[0]=1.;
    v2[1]=0.;
    mphys->InsertMaterialObject(mat->CreateBC(mat, bcidright, dirdirichlet, v1, v2));
    v2[0]=1.;
    v2[1]=0.;
    mphys->InsertMaterialObject(mat->CreateBC(mat, bcidleft, dirdirichlet,  v1, v2));

    v2[0]=0.;
    v2[1]=0.;
    v2[2]=0.;
    mphys->InsertMaterialObject(mat->CreateBC(mat, bcidtopright, pressure,  v1, v2));
    mphys->InsertMaterialObject(mat->CreateBC(mat, bcidtopleft, pressure,  v1, v2));
    mphys->InsertMaterialObject(mat->CreateBC(mat, bcidramp, pressure, v1, v2));

    mphys->AutoBuild();
    mphys->AdjustBoundaryElements();
    mphys->CleanUpUnconnectedNodes();
    mat->SetUpdateMem(false);


    return mphys;
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


#include <iostream>
#include <cmath>
#include <limits>
#include <algorithm>   // <-- precisa disso para std::min

static std::ofstream gnewtlog("poro_newton.log");

static inline REAL SafeNorm(const TPZFMatrix<STATE>& v) {
    REAL s = Norm(v);
    if (!std::isfinite(s)) s = std::numeric_limits<REAL>::infinity();
    return s;
}

static inline REAL BlockNorm(const TPZFMatrix<STATE>& v, int start, int n) {
    // v.Rows() é int64_t; padronize o tipo para evitar ambiguidade do std::min
    const int64_t rows = v.Rows();
    if (start < 0) return 0.;
    if ((int64_t)start >= rows) return 0.;

    const int64_t N = std::min<int64_t>((int64_t)n, rows - (int64_t)start);
    REAL s = 0.;
    for (int64_t k = 0; k < N; ++k) {
        const REAL a = v(start + (int)k, 0);
        s += a * a;
    }
    return std::sqrt(s);
}

bool FindRoot(TPZElastoPlasticAnalysis &an,int &iters,REAL &resu,REAL &resf,int ueqs)
{
    // estado e incremento
    TPZFMatrix<STATE> x(an.Solution()), dx(an.Solution());
    x.Zero(); dx.Zero();

    const REAL tol   = 1e-2;   // igual ao seu
    const int  n_it  = 3;
    const REAL EPS   = 1.e-30;

    if (!gnewtlog.is_open()) { /* nada */ }
    gnewtlog << std::setprecision(10) << std::scientific;

    // resíduo inicial
    an.AssembleResidual();
    REAL rhs_prev = SafeNorm(an.Rhs());
    if (!std::isfinite(rhs_prev)) rhs_prev = 1.0;
    if (rhs_prev < EPS) { iters = 0; resu = 0.; resf = 0.; return true; }

    gnewtlog << "[init] ||R|| = " << rhs_prev << "\n";

    REAL normdu  = 0.0;
    REAL normdu_u= 0.0, normdu_p=0.0;

    for (int it = 1; it <= n_it; ++it) {
        // monta tangente e resolve Δx
        an.Assemble();
        an.Solve();
        dx = an.Solution();
        // atualiza solução candidata
        x += dx;
        an.LoadSolution(x);

        // *** RECOMPUTA o resíduo para a solução ATUALIZADA ***
        an.AssembleResidual();

        // normas do incremento total e por blocos

        normdu   = SafeNorm(dx);
        normdu_u = BlockNorm(dx, 0, ueqs);
        normdu_p = BlockNorm(dx, ueqs, dx.Rows()-ueqs);

        // critério de convergência (resíduo relativo)
        if (normdu_u < tol ) {
            std::cout << "normrhs ="<< normdu_u<< " normdu_p ="<< normdu_p   <<endl;
            return true;
        }
        PoroMechParams prm;
        std::cout<< "[it " << it << "] "
        << " alpha=" << prm.alpha
        << "  ||Δx||=" << normdu
        << "  ||Δu||=" << normdu_u
        << "  ||Δp||=" << normdu_p
        << "\n";
        #ifdef PZ_LOG
        if (logger_poro_main.isDebugEnabled()) {
            std::ostringstream oss;
            oss<< "[it " << it << "] "
            << " alpha=" << prm.alpha
            << "  ||Δx||=" << normdu
            << "  ||Δu||=" << normdu_u
            << "  ||Δp||=" << normdu_p
            << "\n";
            LOGPZ_DEBUG(logger_poro_main, oss.str());
        }
        #endif
       // DebugStop();

    }

    // fallback: aceita se Δx ficou pequeno
    return (resu < tol);
}

#include <TPZSpStructMatrix.h>
#include "pzblockdiag.h"
#include "pzbdstrmatrix.h"
#include "pzblockdiag.h"
#include "pzsubcmesh.h"
// ====== header util (pode ficar no main.cpp mesmo) ==========================
// ---------- setup do solver (BiCGStab + Jacobi em bloco-diagonal) ----------
// ---------- setup do solver (BiCGStab + Jacobi em bloco-diagonal) ----------
static void ConfigureBiCGStabJacobi(TPZElastoPlasticAnalysis& an,
                                    int maxIt, STATE tol, int nThreads,
                                    TPZAutoPointer< TPZBlockDiagonal<STATE> >& Pblock)
{
    // matriz estrutural geral (não assume simetria)
    TPZSpStructMatrix<STATE> strmat(an.Mesh());
    strmat.SetNumThreads(nThreads);               // use 1 p/ reprodutibilidade
    an.SetStructuralMatrix(strmat);

    // cria A (valores virão do Assemble)
    TPZMatrix<STATE>* A = strmat.Create();

    // pré-condicionador: bloco-diagonal por connect (estrutura + valores)
    Pblock = new TPZBlockDiagonal<STATE>();
    TPZBlockDiagonalStructMatrix<STATE> bdiag(an.Mesh());
    bdiag.AssembleBlockDiagonal(*Pblock);

    TPZStepSolver<STATE> Prec;                    // Jacobi sobre o bloco
    Prec.SetMatrix(Pblock);
    Prec.SetJacobi(1, 0., 0);

    TPZStepSolver<STATE> Krylov;                  // BiCGStab(A, Prec)
    Krylov.SetBiCGStab(maxIt, Prec, tol, 0);
    Krylov.SetMatrix(A);

    an.SetSolver(Krylov);
    an.SetPrecond(Prec);
}

// ---------- atualizar valores do pré-condicionador após cada Assemble() ----
static void RefreshPreconditioner(TPZElastoPlasticAnalysis& an,
                                  TPZAutoPointer< TPZBlockDiagonal<STATE> > Pblock)
{
    if (!Pblock) return;
    TPZBlockDiagonalStructMatrix<STATE> bdiag(an.Mesh());
    bdiag.AssembleBlockDiagonal(*Pblock);         // reescreve os VALORES de P
}




// ---- helpers ----
static inline TPZAutoPointer<TPZMatrix<STATE>>
GetGlobalMatrix(TPZLinearAnalysis &an)
{
    auto *ms = dynamic_cast<TPZMatrixSolver<STATE>*>(an.Solver());
    return ms ? ms->Matrix() : TPZAutoPointer<TPZMatrix<STATE>>();
}
// Encapsula auditoria dos blocos K,Q,QT,H,S + solve denso de verificação
// Pré-condições:
//  - an.Solver() é TPZMatrixSolver<STATE>
//  - mat->SetWhichAssemble() com os enums {EK, EQ, EQT, EH, ES} já implementados no material
//  - fUpdateMem == false durante a checagem
static void AuditAssembleAndSolveOnce(TPZElastoPlasticAnalysis &an,
                                      poroplasticmat *mat,
                                      int ueqs, int peqs)
{
    const int neq = ueqs + peqs;

    auto ms = dynamic_cast<TPZMatrixSolver<STATE>*>(an.Solver());
    if(!ms){
        std::cerr << "[Audit] Solver não é TPZMatrixSolver<STATE>.\n";
        return;
    }

    auto assemble_block = [&](typename poroplasticmat::EWhichMatrix which)
    -> TPZFMatrix<STATE>
    {
        mat->SetWhichAssemble(which);
        an.Assemble();                              // monta SÓ o bloco pedido
        TPZFMatrix<STATE> M = *ms->Matrix();        // cópia densa do bloco na posição global
        return M;
    };

    // 1) Monta cada bloco separadamente (os sinais já vêm corretos do material)
    TPZFMatrix<STATE> EK  = assemble_block(mat->EK);
    TPZFMatrix<STATE> EQ  = assemble_block(mat->EQ);
    TPZFMatrix<STATE> EQT = assemble_block(mat->EQT);
    TPZFMatrix<STATE> EH  = assemble_block(mat->EH);
    TPZFMatrix<STATE> ES  = assemble_block(mat->ES);

    // 2) Soma para obter a tangente completa K (os blocos ocupam sub-blocos disjuntos no global)
    TPZFMatrix<STATE> K = EK;  K += EQ;  K += EQT;  K += EH;  K += ES;

    // 3) Monta o resíduo no estado ATUAL (sem atualizar memória!)
    an.AssembleResidual();
    TPZFMatrix<STATE> R = an.Rhs(); // tamanho neq x 1

    // 4) Particionamentos (opcional, só para logs)
    TPZFMatrix<STATE> fu(ueqs,1,0.), fp(peqs,1,0.);
    for(int i=0;i<ueqs;i++)     fu(i,0) = R(i,0);
    for(int i=0;i<peqs;i++)     fp(i,0) = R(ueqs+i,0);

    std::cout << "[Audit] ||fu|| = " << Norm(fu) << "\n";
    std::cout << "[Audit] ||fp|| = " << Norm(fp) << "\n";

    // 5) Resolve K dx = -R (cheque denso “de mesa”)
    TPZFMatrix<STATE> dx(R);
    dx *= -1.;
    TPZFMatrix<STATE> Kcopy(K);      // Solve_LU é destrutivo
    Kcopy.Solve_LU(&dx);

    // 6) Verificação: ||K dx + R|| deve ser ~0
    TPZFMatrix<STATE> Kdx;  K.Multiply(dx, Kdx);
    Kdx += R;
    std::cout << "[Audit] ||K dx + R|| = " << Norm(Kdx) << "\n";

    // 7) (Opcional) Particionar dx para inspeção
    TPZFMatrix<STATE> dx_u(ueqs,1,0.), dx_p(peqs,1,0.);
    for(int i=0;i<ueqs;i++) dx_u(i,0) = dx(i,0);
    for(int i=0;i<peqs;i++) dx_p(i,0) = dx(ueqs+i,0);

    std::cout << "[Audit] ||dx_u|| = " << Norm(dx_u) << "\n";
    std::cout << "[Audit] ||dx_p|| = " << Norm(dx_p) << "\n";

    // (Se quiser, pode carregar dx no estado para ver o pós-update, mas a auditoria termina aqui)
}


int main()
{
    const std::string configfile = "/home/diogo/projects/neopz-master-build-debug/Util/log4cxx.cfg";
    TPZLogger::InitializePZLOG(configfile);



    // 1) Um appender simples no console
    // -------------------- 1) MALHAS --------------------
    //TPZGeoMesh *gmesh   = CreateSingleQuadMesh();
   TPZGeoMesh *gmesh   = TriGMesh(0);

    TPZCompMesh *cmeshU = CMeshElastic (gmesh);
    TPZCompMesh *cmeshP = CMeshPressure(gmesh);

    // -------------------- 2) MULTIFÍSICA --------------------
    TPZCompMesh *mphys = CreateMPhysWithMaterialsAndBCs(gmesh);


    TPZVec<TPZCompMesh*> meshvec(2);
    meshvec[0] = cmeshU;  // u
    meshvec[1] = cmeshP;  // p
    TPZBuildMultiphysicsMesh::AddElements(meshvec, mphys);
    TPZBuildMultiphysicsMesh::AddConnects(meshvec, mphys);
    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, mphys);
    mphys->LoadReferences();
    //mphys->AutoBuild();



    mphys->SetDefaultOrder(order);
    TPZElastoPlasticAnalysis an(mphys, std::cout,TPZElastoPlasticAnalysis::ELineSearch::Dicotomic);
    auto* mat = dynamic_cast<poroplasticmat*>(mphys->FindMaterial(1));
    mat->ResetMemory();
    InitializeMemory(mphys);


    mat->SetUpdateMem(false);
    TPZFStructMatrix<REAL> str(mphys);
    an.SetStructuralMatrix(str);
    TPZStepSolver<REAL> direct;
    direct.SetDirect(ELU);
    an.SetSolver(direct);

    // TPZSkylineStructMatrix<REAL> str(mphys);
    // //TPZFStructMatrix<REAL> str(mphys);
    // an.SetStructuralMatrix(str);
    // TPZStepSolver<REAL> direct;
    // direct.SetDirect(ELDLt);
    // an.SetSolver(direct);

    // TPZAutoPointer< TPZBlockDiagonal<STATE> > Pblock;
    // ConfigureBiCGStabJacobi(an,  500,1.e-8,1, Pblock);
    // --- setup básico ---
    // ===================================================================
    //   Passo de depuração do sistema linear poro-mecânico
    // ===================================================================

    // --- parâmetros de corpo-força e contagens de equações ---
    PoroMechParams prm;
    mat->SetBodyForce(prm.fx, prm.fy);
    const int ueqs = cmeshU->NEquations();
    const int peqs = cmeshP->NEquations();
    const int neq  = ueqs + peqs;
    AuditAssembleAndSolveOnce(an,mat,ueqs, peqs);
    // const int ueqs = cmeshU->NEquations();
    // const int peqs = cmeshP->NEquations();
    // const int neq  = ueqs + peqs;
    //
    // // --- buffers auxiliares ---
    // TPZFMatrix<STATE> x(neq,1,0.), dx(neq,1,0.), rhs;
    // TPZFMatrix<STATE> fu(neq,1,0.), fp(neq,1,0.);
    // TPZFMatrix<STATE> solu(neq,1,0.), solp(neq,1,0.);
    //
    // auto ms = dynamic_cast<TPZMatrixSolver<STATE>*>(an.Solver());
    // mat->SetWhichAssemble(mat->EWhichMatrix::EK);
    // an.Assemble();
    // TPZFMatrix<STATE> EK = *ms->Matrix();
    // mat->SetWhichAssemble(mat->EWhichMatrix::EQ);
    // an.Assemble();
    // TPZFMatrix<STATE> EQ = *ms->Matrix();
    // mat->SetWhichAssemble(mat->EWhichMatrix::EQT);
    // an.Assemble();
    // TPZFMatrix<STATE> EQT = *ms->Matrix();
    // mat->SetWhichAssemble(mat->EWhichMatrix::EH);
    // an.Assemble();
    // TPZFMatrix<STATE> EH = *ms->Matrix();
    // mat->SetWhichAssemble(mat->EWhichMatrix::ES);
    // an.Assemble();
    // TPZFMatrix<STATE> ES = *ms->Matrix();
    //
    // //EH.Print(std::cout);
    //
    //
    //
    // rhs = an.Rhs();
    // dx = an.Solution();
    //
    // for (int i=0; i<ueqs; ++i) fu(i,0) = rhs(i,0);
    // for (int i=ueqs; i<neq; ++i) fp(i,0) = rhs(i,0);
    //
    // for (int i=0; i<ueqs; ++i) solu(i,0) = dx(i,0);
    // for (int i=ueqs; i<neq; ++i) solp(i,0) = dx(i,0);
    //
    //
    // TPZFMatrix<STATE> resp,temp,resup;
    // ES.Multiply(solp,resp);
    // EQ.Multiply(solu,temp);
    // //fq dt+S pn+Q un;
    // resp+=temp+fp;
    //
    // std::cout << "Norm resp =  "<< Norm(resup)<<std::endl;
    // resup=resp;
    // resup+=fu;
    //
    // EK+=EQ;
    // EK+=EQT;
    // EK+=EH;
    // EK+=ES;
    // TPZFMatrix<STATE> sol(resup);
    // EK.Solve_LU(&sol);
    //
    // std::cout << "Norm sol =  "<< Norm(sol)<<std::endl;


    // -------------------------------------------------------------------
    // 1) Monta K e R e resolve K dx = R  (primeira iteração de Newton)
    // -------------------------------------------------------------------
    // an.Assemble();
    // TPZFMatrix<STATE> R0 = an.Rhs();

    // ponteiro para a matriz tangente global
/*
    auto ms = dynamic_cast<TPZMatrixSolver<STATE>*>(an.Solver());
    mat->SetWhichAssemble(mat->EWhichMatrix::EK);
    an.Assemble();
    TPZFMatrix<STATE> EK = *ms->Matrix();
    mat->SetWhichAssemble(mat->EWhichMatrix::EQ);
    an.Assemble();
    TPZFMatrix<STATE> EQ = *ms->Matrix();
    mat->SetWhichAssemble(mat->EWhichMatrix::EQT);
    an.Assemble();
    TPZFMatrix<STATE> EQT = *ms->Matrix();
    mat->SetWhichAssemble(mat->EWhichMatrix::EH);
    an.Assemble();
    TPZFMatrix<STATE> EH = *ms->Matrix();
    mat->SetWhichAssemble(mat->EWhichMatrix::ES);
    an.Assemble();
    TPZFMatrix<STATE> ES = *ms->Matrix();

    //EH.Print(std::cout);



    rhs = an.Rhs();
    dx = an.Solution();

    for (int i=0; i<ueqs; ++i) fu(i,0) = rhs(i,0);
    for (int i=ueqs; i<neq; ++i) fp(i,0) = rhs(i,0);

    for (int i=0; i<ueqs; ++i) solu(i,0) = dx(i,0);
    for (int i=ueqs; i<neq; ++i) solp(i,0) = dx(i,0);


    TPZFMatrix<STATE> resp,temp,resup;
    ES.Multiply(solp,resp);
    EQ.Multiply(solu,temp);
    //fq dt+S pn+Q un;
    resp+=temp+fp;

    std::cout << "Norm resp =  "<< Norm(resup)<<std::endl;
    resup=resp;
    resup+=fu;

    EK+=EQ;
    EK+=EQT;
    EK+=EH;
    EK+=ES;
    TPZFMatrix<STATE> sol(resup);
    EK.Solve_LU(&sol);

    std::cout << "Norm sol =  "<< Norm(sol)<<std::endl;
*/

    // // resolve Δx
    // an.Solve();
    // dx = an.Solution();
    // //RefreshPreconditioner(an,Pblock);
    // // checa consistência: ||K dx − R|| deve ser ~0
    // TPZFMatrix<STATE> Kdx;
    // K.Multiply(dx, Kdx);
    // Kdx -= R0;
    // std::cout << "||K dx - R|| = " << Norm(Kdx) << "\n";
    //
    // // separa residual por blocos u e p
    // rhs = an.Rhs();
    // for (int i=0; i<ueqs; ++i) fu(i,0) = rhs(i,0);
    // for (int i=0; i<peqs; ++i) fp(i,0) = rhs(i+ueqs,0);
    //
    // std::cout << "Norm(fu)  = " << Norm(fu)  << "\n";
    // std::cout << "Norm(fp)  = " << Norm(fp)  << "\n";
    // std::cout << "Norm(dx)  = " << Norm(dx)  << "\n";
    //
    // // transfere para submalhas para inspecionar deslocamento e pressão
    // TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphys);
    // solu = meshvec[0]->Solution();
    // solp = meshvec[1]->Solution();
    //
    // std::cout << "Norm(solu)= " << Norm(solu) << "\n";
    // std::cout << "Norm(solp)= " << Norm(solp) << "\n";
    //
    // // -------------------------------------------------------------------
    // // 2) Atualiza a solução total x = x + dx e monta somente o resíduo
    // // -------------------------------------------------------------------
    // x += dx;
    // an.LoadSolution(x);
    //
    // an.Assemble(); // monta R(x^{k+1}) sem remontar K
    // rhs = an.Rhs();
    //
    // for (int i=0; i<ueqs; ++i) fu(i,0) = rhs(i,0);
    // for (int i=0; i<peqs; ++i) fp(i,0) = rhs(i+ueqs,0);
    //
    // std::cout << "[Residual-only] Norm(fu) = " << Norm(fu) << "\n";
    // std::cout << "[Residual-only] Norm(fp) = " << Norm(fp) << "\n";
    //
    // // -------------------------------------------------------------------
    // // 3) (opcional) Novo passo de Newton: remonta K e resolve outra vez
    // // -------------------------------------------------------------------
    // an.Assemble();
    // TPZFMatrix<STATE> R1 = an.Rhs();
    //
    // an.Solve();
    // dx = an.Solution();
    //
    // // verifica novamente ||K dx − R|| para a nova tangente
    // K = *ms->Matrix();
    // K.Multiply(dx, Kdx);
    // Kdx -= R1;
    // std::cout << "pos-update  ||K dx - R|| = " << Norm(Kdx) << "\n";
    //
    // // transfere a solução atualizada para as submalhas
    // TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphys);
    // solu = meshvec[0]->Solution();
    // solp = meshvec[1]->Solution();
    //
    // std::cout << "pos-update  Norm(dx)   = " << Norm(dx)   << "\n";
    // std::cout << "pos-update  Norm(solu) = " << Norm(solu) << "\n";
    // std::cout << "pos-update  Norm(solp) = " << Norm(solp) << "\n";
    //

    return 0;

    TPZStack<std::string> scal, vecs;
    vecs.Push("Displacement");
    //vecs.Push("PlasticStrain");
    //vecs.Push("ElasticStrain");
    scal.Push("Pressure");
    scal.Push("StrainPlasticJ2");
    //scal.Push("StrainPlasticI1");
    //scal.Push("Coesion");
    //scal.Push("Atrito");

    an.DefineGraphMesh(2, scal, vecs, "sol.vtk");


    std::ofstream out("/home/diogo/projects/neopz-master/Projects2/PoroElastic/kpz.csv");
    std::ofstream oute("/home/diogo/projects/neopz-master/Projects2/PoroElastic/fpz.csv");
    std::ofstream outs("/home/diogo/projects/neopz-master/Projects2/PoroElastic/solpz.csv");
    REAL scale=1.1;
    REAL t=0.;
    REAL dt;
    int g_time=0;




    for(int it=1;it<=1;it++)
    {

        dt = pow(scale,it) - t;

        std::cout<<"Time  = "<< t << " dt = "<< dt << " Normrhs =" << Norm(an.Rhs()) << std::endl;
        mat->SetTimeStep(dt);

        mat->SetUpdateMem(false);

        //an.Assemble();

        //bool ok;
         int iters=5,iters_out;
          STATE resu,resf;

         bool ok=FindRoot(an,iters,resu,resf,ueqs);
        //bool ok=an.FindRoot(iters,resu,resf);

        //bool ok = an.IterativeProcess(std::cout, 1.e-3, 100, true, false, iters_out);

        //PrintKWithSubmeshLabels(an,mphys,cmeshU,cmeshP);

        //an.Solve();


        std::cout<<"Time  = "<< t << " Normrhs =" << Norm(an.Rhs()) << " ok? "<< ok << std::endl;

        an.SetStep(it);

        an.PostProcess(0);

        mat->SetUpdateMem(false);

        t+=dt;

        g_time=t;
    }


    // int iters=30;
    // STATE resu,resf;
    // an.FindRoot(iters,resu,resf);
    //
    // //an.DicotomicLineSearch()
    //
    // TPZStack<std::string> scal, vecs;
    // vecs.Push("Displacement");
    // scal.Push("Pressure");
    //
    // an.DefineGraphMesh(2, scal, vecs, "sol.vtk");
    //
    // an.PostProcess(0);
    //
    // std::cout << " iters = "<<iters <<" resu = "<<resu << " resf = "<<resf << std::endl;
    // bool ok = an.IterativeProcess(std::cout, 1.e-3, 1000, true, false, iters_out);

    return 0;
}
void ApplyLoad(TPZCompMesh* cmesh,
               REAL coes, REAL atrito, TPZManVector<REAL> factors)
{

    // parâmetros de controle
    int nloads = factors.size();
    REAL FS_target = factors[nloads-1];
    REAL fator_atual = 0.0;
    REAL passo_base  = FS_target / REAL(nloads); // passo médio de referência

    cmesh->SetDefaultOrder(2);
    TPZElastoPlasticAnalysis anal(cmesh, std::cout,TPZElastoPlasticAnalysis::ELineSearch::Dicotomic);
    auto* body = dynamic_cast<poroplasticmat*>(cmesh->FindMaterial(1));
    body->ResetMemory();
    InitializeMemory(cmesh);
    if(true)
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
    PoroMechParams prm;
    while( counter<100)
    {
        body->SetBodyForce(prm.fx,prm.fy*fator_atual);


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
// rotula linhas de uma malha H1: "uxN/uyN" se ndof par, "pN" se ndof ímpar
static std::vector<std::string> LabelsFromSubmesh(TPZCompMesh* cmesh, char kind /*'u' ou 'p'*/){
    const int neq = cmesh->NEquations();
    std::vector<std::string> lab(neq, "");
    for (auto cel : cmesh->ElementVec()){
        if (!cel || !cel->Reference()) continue;
        auto *gel = cel->Reference();
        const int ncon = cel->NConnects();
        for (int ic=0; ic<ncon; ++ic){
            const int64_t cidx = cel->ConnectIndex(ic);
            if (cidx < 0) continue;
            TPZConnect &c = cmesh->ConnectVec()[cidx];
            int seq = c.SequenceNumber(); if (seq<0) continue;
            int pos = cmesh->Block().Position(seq), nd = c.NDof();

            // tenta pegar nó de canto correspondente
            int nc = gel->NCornerNodes();
            int which = std::min(ic, nc-1);
            int64_t node = gel->NodeIndex(which);

            if (kind=='u' && nd>=2){
                for(int k=0;k+1<nd;k+=2){
                    if (pos+k   < neq) lab[pos+k]   = "ux"+std::to_string(node);
                    if (pos+k+1 < neq) lab[pos+k+1] = "uy"+std::to_string(node);
                }
            } else if (kind=='p'){
                for(int k=0;k<nd;++k)
                    if (pos+k < neq) lab[pos+k] = "p"+std::to_string(node);
            }
        }
    }
    return lab;
}

// imprime as linhas da K do mphys com rótulos vindos das submalhas
static void PrintKWithSubmeshLabels(TPZLinearAnalysis& an,
                                    TPZCompMesh* mphys,
                                    TPZCompMesh* cmeshU,
                                    TPZCompMesh* cmeshP){
    auto labU = LabelsFromSubmesh(cmeshU,'u');
    auto labP = LabelsFromSubmesh(cmeshP,'p');
    const int off_p = cmeshU->NEquations();

    auto ms = dynamic_cast<TPZMatrixSolver<STATE>*>(an.Solver());
    auto K  = ms->Matrix();
    const int n = K->Rows();

    for (int i=0;i<n;++i){
        std::string tag = (i<off_p) ? labU[i]
        : labP[i-off_p];
        if (tag.empty()) tag = (i<off_p ? "u?" : "p?");
        std::cout << std::setw(8) << tag << " |";
        for (int j=0;j<n;++j) std::cout << " " << K->GetVal(i,j);
        std::cout << "\n";
    }
                                    }

