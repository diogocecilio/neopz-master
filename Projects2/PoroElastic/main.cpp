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
static constexpr double lh   = 0.1;
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

    REAL co[4][2] = {{0.,0.},{0.1,0.},{0.1,1.},{0.,1.}};
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


    TopoLine[0] = 0;
    TopoLine[1] = 1;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( 8, TopoLine, kBC_Bot2, *gmesh );

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
    gmesh->BuildConnectivity();

    cout << "c" << endl;
    for ( int d = 0; d<4; d++ )
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
static constexpr int KMAX = 10;

// ---- 1) VERSÃO ESCALAR (para mat->SetExact) ----
static void ExactP_scalar(const TPZVec<REAL>& x,
                          STATE& p,                  // <- escalar
                          TPZFMatrix<STATE>& grad)   // gradiente (d/dx,d/dy,d/dz)
{
    const REAL y = x[1];
    const REAL lambda = g_E*g_nu / ((1.+g_nu)*(1.-2.*g_nu));
    const REAL G      = g_E / (2.*(1.+g_nu));
    const REAL D = (lambda + 2.*G) * g_perm / g_mu;

    REAL val = 0.0, dpy = 0.0;
    for (int k=0; k<=KMAX; ++k) {
        const REAL m    = (M_PI*0.5) * (2*k+1);
        const REAL arg  = m * y / g_lv;
        const REAL expo = std::exp( -(m*m)*D*g_time/(g_lv*g_lv) );
        val += (2.0/m) * std::sin(arg) * expo;
        dpy += (2.0/m) * (m/g_lv) * std::cos(arg) * expo;
    }
    p = val*1000;
    grad.Redim(3,1); grad.Zero();
    grad(1,0) = dpy; // componente y
}

// ---- 2) ADAPTADOR VETORIAL (p/ TPZDummyFunction e BC/erros) ----
static void ExactP_vec(const TPZVec<REAL>& x,
                       TPZVec<STATE>& v,            // <- vetor tamanho 1
                       TPZFMatrix<STATE>& grad)
{
    STATE p;
    ExactP_scalar(x, p, grad);   // reutiliza cálculo
    v.Resize(1);
    v[0] = p;
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

    v2[0] = 0.0;
    v2[1] = 1.0;
    mphys->InsertMaterialObject(mat->CreateBC(mat, kBC_Top, 3, v1, v2));//Direcional em y
    v2[0] = 1.0;
    v2[1] = 0.0;
    mphys->InsertMaterialObject(mat->CreateBC(mat, kBC_Left, 3, v1, v2));//Direcional em x
    v2[0] = 1.0;
    v2[1] = 0.0;
    mphys->InsertMaterialObject(mat->CreateBC(mat, kBC_Right, 3, v1, v2));//Direcional em x
    v2[0] = 0.0;
    v2[1] = 1000.;
    mphys->InsertMaterialObject(mat->CreateBC(mat, kBC_Bot, 1, v1, v2));//tensao na base

    //mat->SetForcingFunctionBC(ExactP);
    mat->SetExact(ExactP_scalar);

    v2[0] = 0.0;
    v2[1] = 0.0;
    mphys->InsertMaterialObject(mat->CreateBC(mat, kBC_Bot2, 2, v1, v2));//pressao nula no topo

    mphys->AutoBuild();
    mphys->AdjustBoundaryElements();
    mphys->CleanUpUnconnectedNodes();

    return mphys;
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


int main()
{
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
   // TPZSkylineStructMatrix<REAL> str(mphys);
    TPZFStructMatrix<REAL> str(mphys);
    an.SetStructuralMatrix(str);
    TPZStepSolver<REAL> direct;
    direct.SetDirect(ELU);
    an.SetSolver(direct);


     // TPZSSpStructMatrix<STATE> SSpStructMatrix ( mphys );
     // an.SetStructuralMatrix ( SSpStructMatrix );
     // TPZPardisoSolver<REAL> *pardiso = new TPZPardisoSolver<REAL>;
     // an.SetSolver ( *pardiso );

     // an.SetStructuralMatrix(SSpStructMatrix);
     // TPZStepSolver<REAL> direct;
     // direct.SetDirect(ELU);
     // an.SetSolver(direct);

    TPZStack<std::string> scal, vecs;
    vecs.Push("Displacement");
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

    for(int it=1;it<=55;it++)
    {


        dt = pow(scale,it) - t;
        mat->SetTimeStep(dt);
        mat->SetUpdateMem(true);
        an.Assemble();
        mat->SetUpdateMem(false);
        auto *msolver = dynamic_cast<TPZMatrixSolver<STATE>*>(an.Solver());
        TPZFMatrix<STATE> K= *msolver->Matrix();
        //PrintKWithSubmeshLabels(an, mphys, meshvec[0], meshvec[1]);
        TPZFMatrix<STATE> f1=an.Rhs();
        ///f1.Print("f1");
        an.AssembleResidual();
        TPZFMatrix<STATE> f2=an.Rhs();
        //f2.Print("f2");
        f2+=f1;
        an.Rhs()=f2;
        //f2.Print("f2full");
        //std::cout << std::fixed << std::setprecision(20);
        //std::cout << "\nKglob (do solver) =\n";
        //K.Print("Kglob", std::cout, EFormatted);
        //PrintFMat(K,out);
        //PrintFMat(an.Rhs(),oute);
        //PrintFMat(an.Solution(),outs);


        an.Solve();
        //K.Print("KGlobal");
        //an.Rhs().Print("RHS");
        //an.Solution().Print("sol");


        // depois que mphys estiver montada (AddElements/AddConnects/TransferFromMeshes prontos):
       // K.Print("kglob ", std::cout, ECSV);
        // //f2.Print("b ", out, ECSV);

        TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphys);

        an.SetStep(it);
        an.PostProcess(0);


        mat->SetUpdateMem(false);

        //dt=0.1;
        t+=dt;
        g_time=t;
    }

    //an.Solve();

    // auto *msolver = dynamic_cast<TPZMatrixSolver<STATE>*>(an.Solver());
    //
    // auto K = msolver->Matrix();
    // std::cout << std::fixed << std::setprecision(10);
    // std::cout << "\nKglob (do solver) =\n";
    // //K->Print("Kglob", std::cout, EFormatted);
    // K->Print("Kglob");

//     // -------------------- 4) LOOP EM TEMPO (Backward-Euler) --------------------
//
//     const int varU = mat->VariableIndex("Displacement"); // [ux,uy]
//     const int varP = mat->VariableIndex("Pressure");     // p
//                 // ajuste
//     const int  nsteps = 5;          // ajuste
//
//     // gráfico (define uma vez)
//     TPZStack<std::string> scal, vecs;
//     scal.Push("Pressure");
//     vecs.Push("Displacement");
//     an.DefineGraphMesh(2, scal, vecs, "solution.vtk");
//
//     std::cout << std::scientific << std::setprecision(6);
//
//     std::ofstream prof("saida.csv");
//     prof << "step,time,y,uy,p\n";
//     REAL t=0.;
//     REAL dt =0.000001;
//     TPZFMatrix<STATE> u_n = meshvec[0]->Solution();
//     TPZFMatrix<STATE> p_n= meshvec[1]->Solution();
//     TPZFMatrix<STATE> u_n1 = meshvec[0]->Solution();
//     TPZFMatrix<STATE> p_n1= meshvec[1]->Solution();
//     for (int it = 1; it <= 2; ++it)
//     {
//
//         mat->SetTimeStep(dt);
//
//        // mat->SetUpdateMem(false);
//
//         an.Assemble();
//
//         //an.Rhs().Print("RHS");
//
//         an.Solve();
//
//         an.PostProcess(0);
//         u_n= u_n1;
//         p_n = p_n1;
//
//         // prints do passo
//         const TPZFMatrix<STATE> &sol = an.Solution();
//
//         mphys->LoadSolution(sol);
//
//
//         TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphys);
//         u_n1= meshvec[0]->Solution();
//         p_n1= meshvec[1]->Solution();
//
//         u_n-=u_n1;
//         p_n-=p_n1;
//         cout<<  std::scientific << std::setprecision(12)<< "time = "<< t << " dt = "<< dt << " normdu = " << Norm(u_n) << " normdp = " << Norm(p_n) <<endl;
//
//
//         mat->SetUpdateMem(true);
//         an.Assemble();          // <- percorre elementos/IPs, chama Contribute e UpdatePorePressure
//         mat->SetUpdateMem(false);
//
//
//         t+=dt;
//         dt = pow(1.1,it) - t;
//          // você verá 'solution_inc.vtk' com múltiplos steps
//     }

    return 0;
}

/*
int main2()
{
    // -------------------- 1) MALHAS --------------------
    TPZGeoMesh *gmesh   = CreateSingleQuadMesh();
    TPZCompMesh *cmeshU = CMeshElastic (gmesh);
    TPZCompMesh *cmeshP = CMeshPressure(gmesh);

    // -------------------- 2) MULTIFÍSICA --------------------
    TPZCompMesh *mphys = CreateMPhysWithMaterialsAndBCs(gmesh);

    // acopla submalhas → mista (sem AutoBuild)
    TPZVec<TPZCompMesh*> meshvec(2);
    meshvec[0] = cmeshU;  // u
    meshvec[1] = cmeshP;  // p
    TPZBuildMultiphysicsMesh::AddElements(meshvec, mphys);
    TPZBuildMultiphysicsMesh::AddConnects(meshvec, mphys);
    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, mphys);
    mphys->LoadReferences();
    AddBoundaryCompElsToMF(mphys); // cria comp-BCs na mista
    mphys->CleanUpUnconnectedNodes();

    // memória elástica (E,ν) nos IPs
    InitializeMemory(mphys, young, nu);

    // pega o material e índices de pós-processo
    auto *mat = dynamic_cast<TPZMatPoroElastic2DMem<TPZElasticMem>*>(mphys->FindMaterial(kMatVol));
    if (!mat) { std::cerr << "Material poroelástico não encontrado.\n"; return 1; }
    const int varU = mat->VariableIndex("Displacement"); // [ux,uy]
    const int varP = mat->VariableIndex("Pressure");     // p

    // -------------------- 3) ANÁLISE --------------------
    TPZLinearAnalysis an(mphys,false);
    TPZSkylineStructMatrix<REAL> str(mphys);
    an.SetStructuralMatrix(str);
    TPZStepSolver<REAL> direct; direct.SetDirect(ELDLt);
    an.SetSolver(direct);

    REAL dt=1;
    mat->SetTimeStep(dt);

    an.Assemble();
    an.Solve();

    auto *msolver = dynamic_cast<TPZMatrixSolver<STATE>*>(an.Solver());

    auto K = msolver->Matrix();
    std::cout << std::fixed << std::setprecision(10);
    std::cout << "\nKglob (do solver) =\n";
    K->Print("Kglob", std::cout, EFormatted);

//     */
// // ---- PoroLogger.hpp (ou num .cpp seu) -----------------------
// #include <fstream>
// #include <iomanip>
// #include <string>
// #include <cmath>
//
// struct PoroLogger {
//     std::ofstream csv;
//     std::string   path;
//     const char*   logo = "⟦ PZ·PORO ⟧";
//
//     void OpenCSV(const std::string& p) {
//         path = p;
//         csv.open(path, std::ios::out | std::ios::trunc);
//         // cabeçalho
//         csv << "step,time,rhs_norm,rhs_u_norm,rhs_p_norm,"
//         "delta_norm,delta_u_norm,delta_p_norm,"
//         "p_min,p_mean,p_max,uy_top\n";
//         csv.flush();
//     }
//
//     // imprime e escreve CSV em uma chamada
//     void PrintAndCSV(int step, double t,
//                      double rhs_norm, double rhsu_norm, double rhsp_norm,
//                      double delta_norm, double deltau_norm, double deltap_norm,
//                      double pmin, double pmean, double pmax,
//                      double uy_top = 0.0)
//     {
//         std::cout << std::fixed << std::setprecision(6)
//         << logo << " step=" << step
//         << "  t=" << t
//         << "  ||rhs||=" << rhs_norm
//         << " (u=" << rhsu_norm << ", p=" << rhsp_norm << ")"
//         << "  ||Δ||=" << delta_norm
//         << "  ||Δu||=" << deltau_norm
//         << "  ||Δp||=" << deltap_norm
//         << "  p[min,mean,max]=[" << pmin << "," << pmean << "," << pmax << "]"
//         << "  uy_top=" << uy_top
//         << "\n";
//
//         if (csv.is_open()) {
//             csv << step << ','
//             << std::setprecision(16) << t << ','
//             << rhs_norm << ',' << rhsu_norm << ',' << rhsp_norm << ','
//             << delta_norm << ',' << deltau_norm << ',' << deltap_norm << ','
//             << pmin << ',' << pmean << ',' << pmax << ','
//             << uy_top << '\n';
//             csv.flush();
//         }
//     }
// };
//
// // utilitário para separar a norma do RHS entre u e p (assumindo ordenação [u|p])
// inline void SplitRHS_UP(const TPZFMatrix<STATE>& rhs, int nequ_u,
//                         double& rhsu_norm, double& rhsp_norm)
// {
//     double su=0., sp=0.;
//     const int neq = rhs.Rows();
//     for (int i=0; i<nequ_u; ++i) { const double v = rhs(i,0); su += v*v; }
//     for (int i=nequ_u; i<neq;  ++i) { const double v = rhs(i,0); sp += v*v; }
//     rhsu_norm = std::sqrt(su);
//     rhsp_norm = std::sqrt(sp);
// }
//
//
// TPZVec<REAL>  UyAtPoint(TPZCompMesh* cmesh, REAL x, REAL y,int var)
// {
//     TPZGeoMesh* gmesh = cmesh->Reference();
//     TPZManVector<REAL,3> X(3,0.); X[0]=x; X[1]=y;
//     TPZManVector<REAL,3> qsi(3,0.);
//     int64_t elid=0;
//
//     cmesh->LoadReferences();          // liga cada TPZGeoEl ao seu TPZCompEl
//     TPZGeoEl* gel = gmesh->FindElement(X, qsi, elid, /*dim=*/2);
//     if(!gel || !gel->Reference()) DebugStop();
//
//     auto* cel = dynamic_cast<TPZMultiphysicsElement*>(gel->Reference());
//     if(!cel) DebugStop();
//
//     TPZVec<REAL> sol;
//     cel->Solution(qsi, var, sol);
//     return sol; // componente y
// }
// #include "pzcreateapproxspace.h"  // se precisar do ApproxSpace()

// static void AddBoundaryCompElsToMF(TPZCompMesh* mphys)
// {
//     TPZGeoMesh* gmesh = mphys->Reference();
//     const int gdim = gmesh->Dimension();
//     int64_t index = -1;
//
//     // garanta que a mista sabe criar elementos multifísicos (com memória)
//     mphys->SetAllCreateFunctionsMultiphysicElemWithMem();
//
//     for (auto gel : gmesh->ElementVec()) {
//         if (!gel) continue;
//         if (gel->Dimension() != gdim-1) continue;       // só borda
//         const int bcid = gel->MaterialId();
//         if (!mphys->FindMaterial(bcid)) continue;       // só IDs com material BC inserido na mista
//
//         // ✅ uma das duas linhas abaixo (conforme sua versão):
//         mphys->CreateCompEl(gel);
//          //mphys->ApproxSpace().CreateCompEl(gel, *mphys);
//     }
// }
