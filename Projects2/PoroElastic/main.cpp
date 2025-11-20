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
#include "Plasticity/TPZElasticResponse.h"
#include "Elasticity/TPZElasticity2D.h"
#include "DarcyFlow/TPZDarcyFlow.h"

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

using std::cout; using std::endl;



// ---- parâmetros (equivalentes ao seu script) ----
static constexpr int   dim   = 2;
static constexpr int   order = 2;
static constexpr double tick = 1.0;   // (não usado aqui)
static constexpr double lh   = 1.0;
static constexpr double lv   = 1.0;

// materiais/BC ids
static constexpr int kMatVol   =   1;
static constexpr int kBC_Top   =  -1;  // marker 1

static constexpr int kBC_Right =  -2;  // marker 2
static constexpr int kBC_Left  =  -3;  // marker 3
static constexpr int kBC_Bot   =  -4;  // marker 4
static constexpr int kBC_Bot2   =  -8;  // marker 4
static constexpr int kBC_Top2   = -5;  // marker 1
static constexpr int kBC_nodeleft   = -6;  // marker 1
static constexpr int kBC_noderigth   = -7;  // marker 1
double young = 300000.;
double nu    = 0.2;

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

TPZGeoMesh* CreateSingleQuadMesh()
{

    REAL co[4][2] = {{0.,0.},{lh,0.},{lh,lv},{0.,lv}};
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
    for ( int d = 0; d<3; d++ )
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
            mem[idx].fPorePressure=0.;
            TPZVec<REAL> DPorePressure(2,0.);
            mem[idx].fdPorePressure=DPorePressure;
            /**  displacements */
            TPZVec<REAL> fSolU(2,0.);

            /**gradient of u_n */
            TPZFMatrix<REAL> fGradSolU(2,2,0.);
            mem[idx].fSolU=fSolU;
            mem[idx].fGradSolU=fGradSolU;
        }
    }
    pMatWithMem->SetUpdateMem(false);
}

static TPZCompMesh* CMeshElastic(TPZGeoMesh* gmesh){

    auto *cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(2);
    cmesh->SetDefaultOrder(order);
    cmesh->SetAllCreateFunctionsContinuousWithMem();

    auto *mat = new TPZMatElastic2DMem<TPZElasticMem>(1);
    mat->SetId(1);      // (mantido)

    TPZElasticResponse ER; ER.SetEngineeringData(young, nu);
    mat->SetUpdateMem(true);
    mat->SetElasticResponse(ER);
    mat->SetUpdateMem(false);

    cmesh->InsertMaterialObject(mat);
    TPZFMatrix<STATE> val1(2,2,0.);
    TPZVec<REAL> val2(2,0.);
    int dirichlet=0;
    cmesh->InsertMaterialObject(mat->CreateBC(mat, kBC_Bot, dirichlet, val1, val2));
    cmesh->InsertMaterialObject(mat->CreateBC(mat, kBC_Bot2, dirichlet, val1, val2));
    cmesh->InsertMaterialObject(mat->CreateBC(mat, kBC_Top, dirichlet, val1, val2));
    cmesh->InsertMaterialObject(mat->CreateBC(mat, kBC_Left, dirichlet, val1, val2));
    cmesh->InsertMaterialObject(mat->CreateBC(mat, kBC_Right, dirichlet, val1, val2));
    cmesh->InsertMaterialObject(mat->CreateBC(mat, kBC_Top2, dirichlet, val1, val2));
    cmesh->InsertMaterialObject(mat->CreateBC(mat, kBC_nodeleft, dirichlet, val1, val2));
    cmesh->InsertMaterialObject(mat->CreateBC(mat, kBC_noderigth, dirichlet, val1, val2));
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
    cmesh->InsertMaterialObject(mat->CreateBC(mat, kBC_Bot, dirichlet, val1, val2));
    cmesh->InsertMaterialObject(mat->CreateBC(mat, kBC_Bot2, dirichlet, val1, val2));
    cmesh->InsertMaterialObject(mat->CreateBC(mat, kBC_Top, dirichlet, val1, val2));
    cmesh->InsertMaterialObject(mat->CreateBC(mat, kBC_Left, dirichlet, val1, val2));
    cmesh->InsertMaterialObject(mat->CreateBC(mat, kBC_Right, dirichlet, val1, val2));
    cmesh->InsertMaterialObject(mat->CreateBC(mat, kBC_Top2, dirichlet, val1, val2));
    cmesh->InsertMaterialObject(mat->CreateBC(mat, kBC_nodeleft, dirichlet, val1, val2));
    cmesh->InsertMaterialObject(mat->CreateBC(mat, kBC_noderigth, dirichlet, val1, val2));
    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    return cmesh;
}



static REAL g_E   =  young ;
static REAL g_nu  =  nu ;
static REAL g_perm= 1.e-10;
static REAL g_mu  = 1.e-2;
static REAL g_lv  = 1.;
static REAL g_time= 0.0;
static constexpr int KMAX = 20;

// ---- 1) VERSÃO ESCALAR (para mat->SetExact) ----
static void Exact(const TPZVec<REAL>& x,
                          STATE& p,                  // <- escalar
                          TPZFMatrix<STATE>& grad)   // gradiente (d/dx,d/dy,d/dz)
{
    REAL z = x[1];
    //REAL t = g_time;
    REAL L=g_lv;
    const REAL lambda = g_E*g_nu / ((1.+g_nu)*(1.-2.*g_nu));
    const REAL G      = g_E / (2.*(1.+g_nu));
    const REAL c = (lambda + 2.*G) * g_perm / g_mu;

    REAL sigma0=1000.;
    REAL gamma_val=1.;

    REAL sum = 0.0;
    REAL dpdz = 0.0;
    REAL dpdt = 0.0;

    for (int m = 0; m < 30; ++m) {
        REAL n = 2*m + 1;
        REAL coeff = 1.0 / n;
        REAL lambda = n * M_PI / (2.0 * L);
        REAL exponent = - (n*n * M_PI*M_PI * c * g_time) / (4.0 * L*L);
        REAL exp_term = std::exp(exponent);
        REAL sin_term = std::sin(lambda * z);
        REAL cos_term = std::cos(lambda * z);

        sum += coeff * exp_term * sin_term;

        // Derivatives
        dpdz += coeff * exp_term * lambda * cos_term;
        dpdt += coeff * (-n*n * M_PI*M_PI * c / (4.0 * L*L)) * exp_term * sin_term;
    }

    REAL factor = 4.0 * gamma_val * sigma0 / M_PI;
    p = factor * sum;

    grad.Resize(2,1); // [∂p/∂z, ∂p/∂t]
    grad(0,0) = 0;
    grad(1,0) = factor * dpdz;
}
// Exata do DESLOCAMENTO (u = [0, w(y,t)]), com gradiente espacial
// Eq. (6.40): Δw(z,t) = c_m γ σ0 { (L - z) - (8L/π²) Σ_{m=0}^∞ [ 1/(2m+1)² e^{-((2m+1)² π² c t)/(4L²)} cos((2m+1)π z/(2L)) ] }
static void ExactDisp(const TPZVec<REAL>& x,
                      TPZVec<STATE>& u,              // vetor deslocamento [ux, uy]
                      TPZFMatrix<STATE>& grad)       // grad(u): ∂u_j/∂x_i (i=row, j=col)
{
    const REAL y = x[1];      // coordenada vertical (z na equação)
    const REAL L = g_lv;

    // Constantes elásticas e hidráulicas
    const REAL lambda = g_E*g_nu / ((1.+g_nu)*(1.-2.*g_nu));
    const REAL G      = g_E / (2.*(1.+g_nu));

    // difusividade uniaxial c = (λ+2G) k / μ  [m²/s]
    const REAL c = (lambda + 2.*G) * g_perm / g_mu;

    // compressibilidade uniaxial c_m = 1/(λ+2G)  [1/Pa]
    const REAL cm = 1.0 / (lambda + 2.*G);

    const REAL gamma_val = 1;   // eficiência de carregamento γ
    const REAL sigma0    = 1000;  // amplitude de carga σ0

    // Série
    REAL series  = 0.0;  // termo com cos (...)
    REAL dseries = 0.0;  // para ∂w/∂y  (usa sin (...))

    for (int m = 0; m <= KMAX; ++m) {
        const REAL n     = 2*m + 1;
        const REAL alpha = n * M_PI / (2.0 * L);
        const REAL expo  = std::exp( -(n*n * M_PI*M_PI) * c * g_time / (4.0 * L * L) );
        series  += (1.0/(n*n)) * expo * std::cos(alpha * y);
        dseries += (1.0/n)     * expo * std::sin(alpha * y);
    }

    // w(y,t)
    const REAL w = cm * gamma_val * sigma0 *
                   ( (L - y) - (8.0*L/(M_PI*M_PI)) * series );

    // ∂w/∂y = c_m γ σ0 [ -1 + (4/π) Σ (1/(2m+1)) e^{...} sin(...) ]
    const REAL dwdY = cm * gamma_val * sigma0 * ( -1.0 + (4.0/M_PI) * dseries );

    // Saídas
    u.Resize(2);                 // 2D: [ux, uy]
    u[0] = 0.0;
    u[1] = w;

    grad.Redim(2,2);             // grad(u): rows= {x,y}, cols = {ux,uy}
    grad.Zero();
    // ∂ux/∂x = ∂ux/∂y = 0
    // ∂uy/∂x = 0, ∂uy/∂y = dwdY
    grad(1,1) = dwdY;
}


#include "pzfstrmatrix.h"
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
    // propriedades físicas

    double alphaB= 1.0;
    double rhof  = 1000.0;
    double Se     = 1e-10;
    double mu     = 1e-2;
    double perm   = 1e-10;
    mat->SetElasticity(young, nu);
    TPZElasticResponse ER;
    ER.SetEngineeringData(young, nu);
    mat->SetElasticResponse(ER);
    mat->SetAlpha(alphaB);
    mat->SetSe(Se);
    mat->SetPermeability(perm);
    mat->SetViscosity(mu);
    mat->SetRhoF(rhof);
    mat->SetGravity(0.0, 0.0);     // <<< gravidade NÃO zero
    mat->SetBodyForce(0.0, 0.0);
    REAL p0=0.;
    mat->SetMem(p0,ER);
    mphys->InsertMaterialObject(mat);

    TPZFMatrix<STATE> v1(3,3,0.);
    TPZManVector<STATE,3> v2(3,0.);



    v2[0] = 1.0;
    v2[1] = 0.0;
    mphys->InsertMaterialObject(mat->CreateBC(mat, kBC_Left, 3, v1, v2));//Direcional em x
    v2[0] = 1.0;
    v2[1] = 0.0;
    mphys->InsertMaterialObject(mat->CreateBC(mat, kBC_Right, 3, v1, v2));//Direcional em x

    v2[0] = 0.0;
    v2[1] = 1.0;
    mphys->InsertMaterialObject(mat->CreateBC(mat, kBC_Top, 3, v1, v2));//base presa em y


    v2[0] = 0.0;
    v2[1] = 0.0;
    v2[2] = 0.0;//pressao
    mphys->InsertMaterialObject(mat->CreateBC(mat, kBC_Bot2, 2, v1, v2));//pressao zero

    // v2[0] = 0.0;
    // v2[1] = 0.0;
    // v2[2] = 1000;//pressao
    // mphys->InsertMaterialObject(mat->CreateBC(mat, kBC_nodeleft, 2, v1, v2));//pressao zero

    v2[0] = 0.0;
    v2[1] = 1000.;
    mphys->InsertMaterialObject(mat->CreateBC(mat, kBC_Bot, 1, v1, v2));//tensao no top

    mat->SetExact(ExactDisp);
    mat->SetExact(Exact);
    mphys->AutoBuild();
    mphys->AdjustBoundaryElements();
    mphys->CleanUpUnconnectedNodes();

    return mphys;
}


int main()
{
    const std::string configfile = "/home/diogo/projects/neopz-master-build-debug/Util/log4cxx.cfg";
    //TPZLogger::InitializePZLOG(configfile);
    // -------------------- 1) MALHAS --------------------
    TPZGeoMesh *gmesh   = CreateSingleQuadMesh();

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

    InitializeMemory(mphys, young, nu);

    // pega o material para configurar dt e (opcional) histórico
    auto *mat = dynamic_cast<TPZMatPoroElastic2DMem<TPZElasticMem>*>(mphys->FindMaterial(kMatVol));
    if (!mat) { std::cerr << "Material poroelástico não encontrado.\n"; return 1; }



    // -------------------- 3) ANÁLISE --------------------
    TPZLinearAnalysis an(mphys,false);


   // // TPZSkylineStructMatrix<REAL> str(mphys);
    TPZFStructMatrix<REAL> str(mphys);
     an.SetStructuralMatrix(str);
     TPZStepSolver<REAL> direct;
     direct.SetDirect(ELU);
     an.SetSolver(direct);



     // TPZSSpStructMatrix<STATE> str ( mphys );
     // // an.SetStructuralMatrix ( str );
     // // TPZPardisoSolver<REAL> *pardiso = new TPZPardisoSolver<REAL>;
     // // an.SetSolver ( *pardiso );
     //
     // an.SetStructuralMatrix(str);
     // TPZStepSolver<REAL> direct;
     // direct.SetDirect(ELU);
     // an.SetSolver(direct);

    TPZStack<std::string> scal, vecs;
    vecs.Push("Displacement");
    vecs.Push("Flux");
    vecs.Push("ExactPressureGradiendSolution");
    vecs.Push("GradP");
    vecs.Push("ExactDisplacement");
    scal.Push("Pressure");
    scal.Push("ExactPressureSolution");

    an.DefineGraphMesh(2, scal, vecs, "sol.vtk");


    std::ofstream out("/home/diogo/projects/neopz-master/Projects2/PoroElastic/kpz.csv");
    std::ofstream oute("/home/diogo/projects/neopz-master/Projects2/PoroElastic/fpz.csv");
     std::ofstream outs("/home/diogo/projects/neopz-master/Projects2/PoroElastic/solpz.csv");
     REAL scale=1.1;
     REAL t=0.;
     REAL dt=1.e-12;
     g_time=t;

    for(int it=1;it<=1;it++)
    {
        std::cout<<"Time  = "<< t << " Normrhs =" << Norm(an.Rhs()) << std::endl;
        dt = pow(scale,it) - t;

        mat->SetTimeStep(dt);

        mat->SetUpdateMem(true);

        an.Assemble();

        an.Solve();



        an.SetStep(it);

        an.PostProcess(0);

        mat->SetUpdateMem(false);

        t+=dt;

        g_time=t;
    }

    return 0;
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


void PrintFMat(const TPZFMatrix<STATE>& mat,
               std::ostream& out,
               int decimals = 10,   // casas decimais
               char sep = ',')      // separador entre colunas
{
    out.setf(std::ios::fixed);                  // força formato fixo
    out << std::setprecision(decimals);

    const int r = mat.Rows();
    const int c = mat.Cols();

    for (int i = 0; i < r; ++i) {
        for (int j = 0; j < c; ++j) {
            out << mat.GetVal(i, j);
            if (j + 1 < c) out << sep;         // sem sep no fim da linha
        }
        out << '\n';
    }
}
