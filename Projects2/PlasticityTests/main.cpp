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
#include "pzgeoquad.h"
#include "TPZGeoLinear.h"
#include "pzgeotriangle.h"
#include "tpzquadrilateral.h"
#include "pzgnode.h"
#include "tpzarc3d.h"
#include "pzgeopoint.h"

#include "tpzcompmeshreferred.h"
#include "TPZRefPattern.h"
#include "tpzgeoelrefpattern.h"
#include "tpzgeoblend.h"

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
#include "Plasticity/TPZPlasticStepVoigt.h"
#include "Plasticity/TPZVonMises.h"
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
#include "Plasticity/TPZYCVonMisesPV.h"

#include "pzfstrmatrix.h"
using std::cout; using std::endl;



#ifdef PZ_LOG
static TPZLogger logger_plasticity("PlasticityTests");
#endif


typedef TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse> TPlasticMC;

typedef TPZPlasticStepVoigt<TPZYCVonMisesPV, TPZElasticResponse> TPlasticVM;



typedef TPZMatElastoPlastic<TPlasticMC,TPZElastoPlasticMem> plasticmat;

typedef TPZMatElastoPlastic2D<TPlasticVM,TPZElastoPlasticMem> plasticmatvm;


// struct MechParamsMorhCoulomb {
//     // --- parâmetros mecânicos ---
//     STATE young   = 2.e8;
//     STATE nu      = 0.2;      // [-]
//
//     // resistência/atrito (se usados)
//     STATE cohesion = 10000.;
//     STATE phi      = 30.0 * M_PI / 180.0; // rad
//     STATE psi      = 30.0 * M_PI / 180.0; // rad
//
//     // --- carregamentos volumétricos ---
//     STATE fx =  0.0;          //
//     STATE fy = 0.0;        //
//     int dim=3;
//     TPZManVector<STATE,3> BodyForce = {0.0, 0.0, 0.0};
//
// };

struct MechParamsMorhCoulomb {
    // --- parâmetros mecânicos (SI) ---
    STATE young   = 3.07e7;//kPa                 // ≈ 30.7 GPa (Ec ≈ 5600*sqrt(fck[MPa]))
    STATE nu      = 0.20;                    // [-]

    // resistência/atrito (calibrado p/ fck ≈ 30 MPa)
    STATE cohesion = 7.5e3;//kPa                  // ≈ 7.5 MPa (c = σc*(1 - sinφ)/(2 cosφ), φ≈37°)
    STATE phi      = 37.0 * M_PI / 180.0;    // rad
    STATE psi      = 10.0 * M_PI / 180.0;    // rad (não associado, dilatância menor que φ)

    // --- carregamentos volumétricos (peso próprio) ---
    STATE fx = 0.0;
    STATE fy =0.; //-24.;//kPa   // N/m³  (γ ≈ 24 kN/m³ → direção -y)
    int dim = 3;
    TPZManVector<STATE,3> BodyForce = {0.0, fy, 0.0};
};

struct MechParamsVonMises{
    // --- parâmetros mecânicos (SI) ---
    STATE young   = 210;//GPa
    STATE nu      = 0.30;

    STATE sigmay = 0.24;//GPa

    // --- carregamentos volumétricos (peso próprio) ---
    STATE fx = 0.0;
    STATE fy =0.;
    int dim = 3;
    TPZManVector<STATE,3> BodyForce = {0.0, 0.0, 0.0};
};

struct BoundaryIndexes {
    // --- parâmetros mecânicos ---
public:
    int bcbottom=-1;
    int bctop=-2;
    int bcright=-3;
    int bcback=-4;
    int bcleft=-5;
    int bcfront=-6;
    int bcnode0=-7;//{b2,0.,0.}
    int bcnode1=-8;//{b2,b1,0.}
    int bcnode2=-9;//{0.,b1,0}
    int bcnode3=-10;//{0.,0.,0.}
    int bcnode4=-11;//{b2,0.,h}

};
struct IndexesCylinder {
    // --- parâmetros mecânicos ---
public:
    int matvolume=1;
    int bcinner=-1;
    int bcouter=-2;
    int bcbottom=-3;
    int bctop=-4;;

};


TPZGeoMesh* CubeMesh()
{

    STATE b1=0.5,b2=0.5,h=1.;
    REAL co[8][3] = {{b2,0.,0.},{b2,b1,0.},{0.,b1,0},{0.,0.,0.},
                     {b2,0.,h},{b2,b1,h},{0.,b1,h},{0.,0.,h}};
    long indices[1][8] = {{0,1,2,3,4,5,6,7}};
    TPZGeoEl *elvec[1];
    TPZGeoMesh *gmesh = new TPZGeoMesh();
    gmesh->SetDimension ( 3 );
    long nnode = 8;
    long nod;
    for ( nod=0; nod<nnode; nod++ )
    {
        long nodind = gmesh->NodeVec().AllocateNewElement();
        TPZVec<REAL> coord ( 3 );
        coord[0] = co[nod][0];
        coord[1] = co[nod][1];
        coord[2] = co[nod][2];
        gmesh->NodeVec() [nodind] = TPZGeoNode ( nod,coord,*gmesh );
    }

    long el;
    long nelem = 1;
    long index=0;
    for ( el=0; el<nelem; el++ )
    {
        TPZVec<long> nodind ( nnode );
        for ( nod=0; nod<nnode; nod++ ) nodind[nod]=indices[el][nod];
        //    elvec[el] = new TPZGeoElQ2d(el,nodind,1);

        elvec[el] = gmesh->CreateGeoElement ( ECube,nodind,1,index );
    }

    TPZVec <long> TopoQuad ( 4 );

    BoundaryIndexes bcindexes;


    index++;
    TopoQuad[0] = 0;
    TopoQuad[1] = 1;
    TopoQuad[2] = 2;
    TopoQuad[3] = 3;
    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> ( index, TopoQuad, bcindexes.bcbottom, *gmesh );

    index++;
    TopoQuad[0] = 4;
    TopoQuad[1] = 5;
    TopoQuad[2] = 6;
    TopoQuad[3] = 7;
    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> ( index, TopoQuad, bcindexes.bctop, *gmesh );

    index++;
    TopoQuad[0] = 1;
    TopoQuad[1] = 2;
    TopoQuad[2] = 6;
    TopoQuad[3] = 5;
    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> ( index, TopoQuad, bcindexes.bcright, *gmesh );

    index++;
    TopoQuad[0] = 3;
    TopoQuad[1] = 2;
    TopoQuad[2] = 6;
    TopoQuad[3] = 7;
    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> ( index, TopoQuad, bcindexes.bcback, *gmesh );

    index++;
    TopoQuad[0] = 3;
    TopoQuad[1] = 7;
    TopoQuad[2] = 4;
    TopoQuad[3] = 0;
    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> ( index, TopoQuad, bcindexes.bcleft, *gmesh );

    index++;
    TopoQuad[0] = 0;
    TopoQuad[1] = 1;
    TopoQuad[2] = 5;
    TopoQuad[3] = 4;
    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> ( index, TopoQuad, bcindexes.bcfront, *gmesh );

    TPZVec <long> TopoNode ( 1 );
    index++;
    TopoNode[0]=0;
    new TPZGeoElRefPattern< pzgeom::TPZGeoPoint> ( index, TopoNode, bcindexes.bcnode0, *gmesh );
    index++;
    TopoNode[0]=1;
    new TPZGeoElRefPattern< pzgeom::TPZGeoPoint> ( index, TopoNode, bcindexes.bcnode1, *gmesh );
    index++;
    TopoNode[0]=2;
    new TPZGeoElRefPattern< pzgeom::TPZGeoPoint> ( index, TopoNode, bcindexes.bcnode2, *gmesh );
    index++;
    TopoNode[0]=3;
    new TPZGeoElRefPattern< pzgeom::TPZGeoPoint> ( index, TopoNode, bcindexes.bcnode3, *gmesh );
    index++;
    TopoNode[0]=4;
    new TPZGeoElRefPattern< pzgeom::TPZGeoPoint> ( index, TopoNode, bcindexes.bcnode4, *gmesh );
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
#include "pzgeotetrahedra.h"
#include "pzgeopyramid.h"
#include "TPZRefPatternTools.h"
#include "pzgeoelbc.h"

#include "pzgeoquad.h"
#include "TPZGeoLinear.h"
#include "pzgeotriangle.h"
#include "tpzgeoelrefpattern.h"
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "TPZVTKGeoMesh.h"
#include "tpzquadraticquad.h"
TPZGeoMesh* PressurizedCylinderMesh()
{

    IndexesCylinder cylindexes;
    STATE ri=100.;//mm
    STATE re=200.;//mm
    STATE h=20.;
    STATE theta=90*M_PI/180.;
    STATE s=sin(theta);
    STATE c=cos(theta);
    STATE s2=sin(theta/2.);
    STATE c2=cos(theta/2.);

    REAL co[6][2] = {

    {ri,0},{re,0},
    {ri*c,ri*s},{re*c,re*s},
    {ri*c2,ri*s2},{re*c2,re*s2}


    };

    TPZGeoEl *elvec[1];
    TPZGeoMesh *gmesh = new TPZGeoMesh();
    gmesh->SetDimension ( 2 );
    long nnode = 6;
    long nod;
    for ( nod=0; nod<nnode; nod++ )
    {
        long nodind = gmesh->NodeVec().AllocateNewElement();
        TPZVec<REAL> coord ( 2 );
        coord[0] = co[nod][0];
        coord[1] = co[nod][1];
        gmesh->NodeVec() [nodind] = TPZGeoNode ( nod,coord,*gmesh );
    }

    long id=0;
    TPZVec<long> TopolQuad(4);
    TopolQuad[0] = 0;
    TopolQuad[1] = 1;
    TopolQuad[2] = 3;
    TopolQuad[3] = 2;
    new TPZGeoElRefPattern< pzgeom::TPZGeoBlend< pzgeom::TPZGeoQuad> > (id,TopolQuad,cylindexes.matvolume,*gmesh);

    id++;
    TPZVec<long> TopolArc(3);
    TopolArc[0] = 0;
    TopolArc[1] = 2;
    TopolArc[2] = 4;
    new TPZGeoElRefPattern< pzgeom::TPZArc3D > (id,TopolArc,cylindexes.bcinner,*gmesh);

    id++;
    TopolArc[0] = 1;
    TopolArc[1] = 3;
    TopolArc[2] = 5;
    new TPZGeoElRefPattern< pzgeom::TPZArc3D > (id,TopolArc,cylindexes.bcouter,*gmesh);

    id++;
    TPZVec<long> TopolLine(2);
    TopolLine[0]=0;
    TopolLine[1]=1;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,cylindexes.bcbottom,*gmesh);

    id++;
    TopolLine[0]=3;
    TopolLine[1]=2;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,cylindexes.bctop,*gmesh);

    cout << "b" << endl;
    gmesh->BuildConnectivity();
    // cout << "c" << endl;
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

    std::ofstream files ( "teste-blend.vtk" );
    TPZVTKGeoMesh::PrintGMeshVTK ( gmesh,files,false );
    cout << "c" << endl;
    return gmesh;
}

static TPZCompMesh* CompMeshCyl(TPZGeoMesh* gmesh)
{
    MechParamsVonMises param;
    auto *mphys = new TPZCompMesh(gmesh);
    mphys->SetDimModel(2);
    mphys->SetAllCreateFunctionsContinuousWithMem();
    mphys->SetDefaultOrder(1);
    auto mat = new plasticmat(1);


     TPZElasticResponse ER;
     ER.SetEngineeringData(param.young, param.nu);

     //TPlasticVM vm;
     //vm.fYC.SetUp(param.sigmay, ER);
     //vm.fER = ER;
    //mat->SetPlasticityModel(vm);
     mat->SetId(1);
     //REIMPLEMENTAR
//     DebugStop();
     //body->SetLoadFactor ( factor );
    //mat->SetBodyForce(param.BodyForce);
     mphys->InsertMaterialObject(mat);//0
     mat->Print(std::cout);

    TPZFMatrix<STATE> v1(3,3,0.);
    TPZManVector<STATE,3> v2(3,0.);
    int dirdirichlet=3,pressure=5;

    IndexesCylinder bcindexes;
    v2[0]=0.;
    v2[1]=1.;
    mphys->InsertMaterialObject(mat->CreateBC(mat,bcindexes.bcbottom , dirdirichlet, v1, v2));

    v2[0]=1.;
    v2[1]=0.;
    mphys->InsertMaterialObject(mat->CreateBC(mat,bcindexes.bctop , dirdirichlet, v1, v2));

    v2[0]=0.19;
    v2[1]=0.;
    mphys->InsertMaterialObject(mat->CreateBC(mat,bcindexes.bcinner , pressure, v1, v2));

    mphys->AutoBuild();
    mphys->AdjustBoundaryElements();
    mphys->CleanUpUnconnectedNodes();
    mphys->Print(cout);
    return mphys;
}


void PostProcessVariables(TPZStack<std::string>& scal, TPZStack<std::string>& vec)
{
    //scal.Push ( "Yield" );
    //scal.Push ( "POrder" );
   // scal.Push ( "Atrito" );
    //scal.Push ( "Coesion" );
    //scal.Push ( "StrainPlasticJ2" );
    //scal.Push ( "VolHardening" );
    vec.Push ( "Displacement" );
    //vec.Push ( "ShearPlasticDeformation" );
    //vec.Push ( "PlasticDeformation" );
    //vec.Push ( "StrainElastic" );
    scal.Push ( "StrainElasticJ2" );

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

void PostElastoplastic(TPZCompMesh* cmesh,const std::string& vtkfile,int matid,int step)
{
    TPZPostProcAnalysis pproc;
    CreatePostProcessingMesh(cmesh, &pproc, matid);
    TPZStack<std::string> scal, vec;
    PostProcessVariables(scal, vec);
    pproc.DefineGraphMesh(/*dim=*/2, scal, vec, vtkfile);
    pproc.SetStep(step);
    pproc.PostProcess(0);
}
static TPZCompMesh* CompMesh(TPZGeoMesh* gmesh)
{
    MechParamsMorhCoulomb param;
    auto *mphys = new TPZCompMesh(gmesh);
    mphys->SetDimModel(3);
    mphys->SetAllCreateFunctionsContinuousWithMem();
    mphys->SetDefaultOrder(3);
    auto mat = new plasticmat(1);

    TPZElasticResponse ER;
    ER.SetEngineeringData(param.young, param.nu);



    TPlasticMC mc;
    mc.fYC.SetUp(param.phi, param.psi, param.cohesion, ER);
    STATE epspbar=0.;
    mc.fYC.SetEpsBar(epspbar);
    mc.fER = ER;
    mc.SetStrengthReductionFactor(1.0);
    mc.fER.Print(std::cout);

    mat->SetPlasticityModel(mc);
    mat->SetId(1);

    mat->SetBulkDensity(1.0);


    mphys->InsertMaterialObject(mat);//0
    TPZFMatrix<STATE> v1(3,3,0.);
    TPZManVector<STATE,3> v2(3,0.);
    int dirdirichlet=3,newmann=1;

    BoundaryIndexes bcindexes;
    v2[0]=0.;
    v2[1]=0.;
    v2[2]=1.;
    mphys->InsertMaterialObject(mat->CreateBC(mat,bcindexes.bcbottom , dirdirichlet, v1, v2));

    v2[0]=1.;
    v2[1]=0.;
    v2[2]=0.;
    mphys->InsertMaterialObject(mat->CreateBC(mat,bcindexes.bcfront , dirdirichlet, v1, v2));


    v2[0]=0.;
    v2[1]=1.;
    v2[2]=0.;
    mphys->InsertMaterialObject(mat->CreateBC(mat,bcindexes.bcright , dirdirichlet, v1, v2));//9


    v2[0]=0.;
    v2[1]=0.;
    v2[2]=-100000.;
    mphys->InsertMaterialObject(mat->CreateBC(mat,bcindexes.bctop , newmann, v1, v2));//13
    // v2[0]=0.;
    // v2[1]=0.;
    // v2[2]=1.;
    // mphys->InsertMaterialObject(mat->CreateBC(mat,bcindexes.bcquadbase0 , dirdirichlet, v1, v2));//5
    // mphys->InsertMaterialObject(mat->CreateBC(mat,bcindexes.bcquadbase1 , dirdirichlet, v1, v2));//6
    // mphys->InsertMaterialObject(mat->CreateBC(mat,bcindexes.bcquadbase2 , dirdirichlet, v1, v2));//7
    // mphys->InsertMaterialObject(mat->CreateBC(mat,bcindexes.bcquadbase3 , dirdirichlet, v1, v2));//8
    //
    //
    //
    // v2[0]=0.;
    // v2[1]=0.;
    // v2[2]=100.;
    // mphys->InsertMaterialObject(mat->CreateBC(mat,bcindexes.bcquadtopo0 , newmann, v1, v2));//9
    // mphys->InsertMaterialObject(mat->CreateBC(mat,bcindexes.bcquadtopo1 , newmann, v1, v2));//10
    // mphys->InsertMaterialObject(mat->CreateBC(mat,bcindexes.bcquadtopo2 , newmann, v1, v2));//11
    // mphys->InsertMaterialObject(mat->CreateBC(mat,bcindexes.bcquadtopo3 , newmann, v1, v2));//12
    //
    // v2[0]=1.;
    // v2[1]=1.;
    // v2[2]=1.;
    // mphys->InsertMaterialObject(mat->CreateBC(mat,bcindexes.bcnobase , dirdirichlet, v1, v2));//13

    mphys->AutoBuild();
    mphys->AdjustBoundaryElements();
    mphys->CleanUpUnconnectedNodes();
    //mphys->Print(cout);
    return mphys;
}
REAL UyAtNode(TPZCompMesh* cmesh, REAL x, REAL y,REAL z)
{
    cmesh->LoadReferences();
    auto* gmesh = cmesh->Reference();
    if (!gmesh) return 0.0;

    TPZManVector<REAL,3> X(3,0.0); X[0]=x; X[1]=y, X[2]=z;
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
        if (std::fabs(co[0]-x) < 1e-10 && std::fabs(co[1]-y) < 1e-10&& std::fabs(co[2]-z) < 1e-10){ local = i; break; }
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
    if (pos+2 >= sol.Rows() || ndof < 3) return 0.0;

    return sol(pos+2,0);
}
void ApplyLoad(TPZCompMesh* cmesh,TPZManVector<REAL> factors)
{
    BoundaryIndexes bcindexes;
    // parâmetros de controle
    int nloads = factors.size();

    TPZElastoPlasticAnalysis anal(cmesh, std::cout,TPZElastoPlasticAnalysis::ELineSearch::GoldenSection);

    auto* bcmat = dynamic_cast<TPZBndCondT<STATE>*>(cmesh->FindMaterial(bcindexes.bctop));

    if(!bcmat)
    {
        std::cout << "material do contorno nao encontrado"<<std::endl;
        DebugStop();

    }

    const REAL load0=bcmat->Val2()[2];

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
    REAL x=0.5;
    REAL y=0.5;
    REAL z=1.;
    REAL uy=0.;
    int counter=0;
    int old_iters_out=0;

    int matid=1;



     std::string vtkfile="applyload.vtk";

    for(int i =0;i< nloads;i++)
    {
         REAL fator_atual=factors[i];
        bcmat->Val2()[2]=load0*factors[i];
        cout << "Load step =" << i<<" factor =  "<<factors[i] <<endl;
        int iters_out;
        REAL resu,resf;
        //anal.FindRoot(iters_out,resu,resf);
         bool ok = anal.IterativeProcess(std::cout, 1.e-6,7, true, false, iters_out);

         TPZFMatrix<REAL> tempsol=anal.Solution();
         anal.AcceptSolution();

         //PostElastoplastic(cmesh,vtkfile,matid,i);
         cmesh->LoadSolution(tempsol);
         anal.LoadSolution(cmesh->Solution());
        // uy += UyAtNode(cmesh,  x,  y, z);
         csv << counter << "," << -uy << "," << fator_atual << "," << iters_out << "," << 1 << "\n";
         std::cout << "uy = " << uy << "  factor = " << fator_atual
         << "  iters = " << iters_out <<  "  old_iters_out = " << old_iters_out << " counter =" << counter<< std::endl;

    }

}

int main()
{
    const std::string configfile = "/home/diogo/projects/neopz-master-build-debug/Util/log4cxx.cfg";
    TPZLogger::InitializePZLOG(configfile);

    auto gmesh = PressurizedCylinderMesh();
    auto cmesh = CompMeshCyl(gmesh);

   // auto gmesh =CreateLCircleGeoMesh();
   //auto gmesh = CubeMesh();


   //  auto cmesh = CompMesh(gmesh);
   //
  //  InitializeMemory(cmesh);
     TPZElastoPlasticAnalysis an(cmesh, std::cout,TPZElastoPlasticAnalysis::ELineSearch::Armijo);
   //
   //  //mat->ResetMemory();
   //
     TPZFStructMatrix<STATE> str(cmesh);
     an.SetStructuralMatrix(str);
     TPZStepSolver<REAL> direct;
     direct.SetDirect(ELU);
     an.SetSolver(direct);
    int itersout;
    bool ok = an.IterativeProcess(std::cout, 1.e-6, 5, true, false, itersout);
   REAL resu,resf;
   int iter;
   // an.FindRoot(iter,resu,resf);
   //  int iters=3;
   //  REAL resu;
   //  REAL resf;
    //an.Assemble();
    //an.Rhs().Print("rhs");
    //an.Solve();
    //an.Solution().Print("sol");
    an.AcceptSolution();
    std::string vtkfile="load.vtk";

    PostElastoplastic(cmesh,vtkfile,1,0);
   //
   // // auto ms = dynamic_cast<TPZMatrixSolver<STATE>*>(an.Solver());
   //                            // monta SÓ o bloco pedido
   // // TPZFMatrix<STATE> M = *ms->Matrix();
   //  //M.Print("mat");
   //
   //
   //  // an.Assemble();
   //  // an.Solve();
   //  // an.AcceptSolution(0);
   //    int itersout;
   //    bool ok = an.IterativeProcess(std::cout, 1.e-6, iters, true, false, itersout);
   //    an.AcceptSolution();
   //    std::string vtkfile="load.vtk";
   //    PostElastoplastic(cmesh,vtkfile,1,0);
   //    return 0;
   //  // // //an.FindRoot(iters,resu,resf);
   //  //  REAL norm = an.AcceptSolution();
   //  //  cout << "norm = "<<norm <<endl;
   //   //std::string vtkfile="out.vtk";
   //   //int matid=1;
   //   //PostElastoplastic(cmesh,vtkfile,matid,0);
   //  // TPZElastoPlasticAnalysis an2(cmesh, std::cout,TPZElastoPlasticAnalysis::ELineSearch::Armijo);
   //  // BoundaryIndexes bcindexes;
   //  // auto* bcmat = dynamic_cast<TPZBndCondT<STATE>*>(cmesh->FindMaterial(bcindexes.bctop));
   //  // auto* mat = dynamic_cast<plasticmat*>(cmesh->FindMaterial(1));
   //  // TPZFStructMatrix<REAL> str2(cmesh);
   //  // an2.SetStructuralMatrix(str2);
   //  // TPZStepSolver<REAL> direct2;
   //  // direct2.SetDirect(ELU);
   //  // an2.SetSolver(direct2);
   //  //  int itersout2;
   //  //  //bcmat->Val2()[2]*=1/2.;
   //  //  mat->SetUpdateMem(true);
   //  // ok = an2.IterativeProcess(std::cout, 1.e-6, 5, true, false, itersout2);
   //
   //
   //  int nloads = 10;
   //  const REAL FS_target = 1. ;          // mantém seu “+0.1”
   //  TPZManVector<REAL> factors(nloads+1);   // 0 .. nloads (inclusivo)
   //  for (int i = 0; i <= nloads; ++i) {
   //      factors[i] = FS_target * REAL(i) / REAL(nloads); // 0, Δ, 2Δ, …, FS_target
   //      //cout<< factors[i] <<endl;
   //  }
   //  ApplyLoad( cmesh,factors);
   //
   //  //int itersout;
   //  //
    return 0;
}
struct BoundaryIndexes2 {
    // --- parâmetros mecânicos ---
public:

    int bcquadbase0=-1;
    int bcquadbase1=-2;
    int bcquadbase2=-3;
    int bcquadbase3=-4;

    int bcquadtopo0=-5;
    int bcquadtopo1=-6;
    int bcquadtopo2=-7;
    int bcquadtopo3=-8;

    int bcnobase=-9;


};



TPZGeoMesh* CubeMesh2()
{

    STATE b1=1.,b2=1.,h=1.;
    REAL co[18][3] = {

        {b1,0.,0.},{b1,b2/2,0.},{b1,b2,0},
        {b1/2,0.,0.},{b1/2,b2/2,0.},{b1/2,b2,0},
        {0,0.,0.},{0,b2/2,0.},{0,b2,0},

        {b1,0.,h},{b1,b2/2,h},{b1,b2,h},
        {b1/2,0.,h},{b1/2,b2/2,h},{b1/2,b2,h},
        {0,0.,h},{0,b2/2,h},{0,b2,h},


    };
    long indices[4][8] = {
        {0,1,4,3,9,10,13,12},
        {1,2,5,4,10,11,14,13},
        {4,5,8,7,13,14,17,16},
        {3,4,7,6,12,13,16,15}
    };
    TPZGeoEl *elvec[1];
    TPZGeoMesh *gmesh = new TPZGeoMesh();
    gmesh->SetDimension ( 3 );
    long nnode = 18;
    long nod;
    for ( nod=0; nod<nnode; nod++ )
    {
        long nodind = gmesh->NodeVec().AllocateNewElement();
        TPZVec<REAL> coord ( 3 );
        coord[0] = co[nod][0];
        coord[1] = co[nod][1];
        coord[2] = co[nod][2];
        gmesh->NodeVec() [nodind] = TPZGeoNode ( nod,coord,*gmesh );
    }

    long el;
    long nelem = 4;
    long index=0;
    for ( el=0; el<nelem; el++ )
    {
        TPZVec<long> nodind ( nnode );
        for ( nod=0; nod<nnode; nod++ ) nodind[nod]=indices[el][nod];
        //    elvec[el] = new TPZGeoElQ2d(el,nodind,1);

        elvec[el] = gmesh->CreateGeoElement ( ECube,nodind,1,index );
    }

    TPZVec <long> TopoQuad ( 4 );

    BoundaryIndexes2 bcindexes;


    index++;
    TopoQuad[0] = 0;
    TopoQuad[1] = 1;
    TopoQuad[2] = 4;
    TopoQuad[3] = 3;
    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> ( index, TopoQuad, bcindexes.bcquadbase0, *gmesh );

    index++;
    TopoQuad[0] = 1;
    TopoQuad[1] = 2;
    TopoQuad[2] = 5;
    TopoQuad[3] = 4;
    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> ( index, TopoQuad, bcindexes.bcquadbase1, *gmesh );

    index++;
    TopoQuad[0] = 4;
    TopoQuad[1] = 5;
    TopoQuad[2] = 8;
    TopoQuad[3] = 7;
    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> ( index, TopoQuad,bcindexes.bcquadbase2, *gmesh );

    index++;
    TopoQuad[0] = 3;
    TopoQuad[1] = 4;
    TopoQuad[2] = 7;
    TopoQuad[3] = 6;
    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> ( index, TopoQuad, bcindexes.bcquadbase3, *gmesh );

    index++;
    TopoQuad[0] = 9;
    TopoQuad[1] =10;
    TopoQuad[2] = 13;
    TopoQuad[3] = 12;
    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> ( index, TopoQuad,bcindexes.bcquadtopo0, *gmesh );

    index++;
    TopoQuad[0] = 10;
    TopoQuad[1] = 11;
    TopoQuad[2] = 14;
    TopoQuad[3] = 13;
    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> ( index, TopoQuad, bcindexes.bcquadtopo1, *gmesh );

    index++;
    TopoQuad[0] = 13;
    TopoQuad[1] = 14;
    TopoQuad[2] = 17;
    TopoQuad[3] = 16;
    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> ( index, TopoQuad, bcindexes.bcquadtopo2, *gmesh );

    index++;
    TopoQuad[0] = 12;
    TopoQuad[1] = 13;
    TopoQuad[2] = 16;
    TopoQuad[3] = 15;
    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> ( index, TopoQuad, bcindexes.bcquadtopo3, *gmesh );

    TPZVec <long> TopoNode ( 1 );
    index++;
    TopoNode[0]=4;
    new TPZGeoElRefPattern< pzgeom::TPZGeoPoint> ( index, TopoNode, bcindexes.bcnobase, *gmesh );

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

