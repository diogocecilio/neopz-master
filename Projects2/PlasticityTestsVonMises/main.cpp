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
#include "Plasticity/TPZYCVonMisesVoigt.h"

#include "pzfstrmatrix.h"

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

using std::cout; using std::endl;



#ifdef PZ_LOG
static TPZLogger logger_plasticity("PlasticityTests");
#endif

typedef TPZPlasticStepVoigt<TPZYCVonMisesVoigt, TPZElasticResponse> TPlasticStepVoigtVM;

typedef TPZMatElastoPlastic2D<TPlasticStepVoigtVM,TPZElastoPlasticMem> TMatElastoPlaticVoigtVM;

typedef TPZMatElastoPlastic<TPlasticStepVoigtVM,TPZElastoPlasticMem> TMatElastoPlaticVoigtVM3D;

void ApplyLoad(TPZCompMesh* cmesh,TPZManVector<REAL> factors,int loaddir,int indexbc);

void PostElastoplastic(TPZCompMesh* cmesh,const std::string& vtkfile,int matid,int step,int dim);

void PostProcessVariables(TPZStack<std::string>& scal, TPZStack<std::string>& vec);

void CreatePostProcessingMesh(TPZCompMesh* cmesh,TPZPostProcAnalysis* pproc,int matid);

REAL IterativeProcessArcLength(TPZElastoPlasticAnalysis &an,TPZCompMesh* cmesh,int loaddir,int indexbc,std::string vtkfile);

REAL IterativeProcessArcLength(TPZElastoPlasticAnalysis &an,
                               int loaddir,
                               int indexbc,
                               int nsteps,
                               int nloadcicles,
                               STATE bccondval,
                               STATE lambda0,
                               STATE L0,
                               std::string vtkfile);
struct MechParamsVonMises{
    // --- parâmetros mecânicos (SI) ---
    STATE young   = 210;//GPa
    STATE nu      = 0.30;

    STATE sigmay = 0.24;//GPa
    STATE H0=0.;

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
    int bctop=-4;

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
    for ( int d = 0; d<1; d++ )
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
    mphys->SetDefaultOrder(2);
    auto mat = new TMatElastoPlaticVoigtVM(1);


     TPZElasticResponse ER;
     ER.SetEngineeringData(210,0.3);

     TPZYCVonMisesVoigt vmyc;
     const STATE sigmaY0 = 0.24;
     const STATE Hiso    = 0.;
     vmyc.SetUp(sigmaY0,Hiso);


     TPlasticStepVoigtVM PlasticStepVoigt;



     PlasticStepVoigt.SetPlasticCriterion(vmyc);
     PlasticStepVoigt.SetElasticResponse(ER);



     mat->SetPlasticityModel(PlasticStepVoigt);
     mat->SetId(1);
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

    v2[0]=-0.19209;
    v2[1]=0.;
    mphys->InsertMaterialObject(mat->CreateBC(mat,bcindexes.bcinner , pressure, v1, v2));

    mphys->AutoBuild();
    mphys->AdjustBoundaryElements();
    mphys->CleanUpUnconnectedNodes();
    //mphys->Print(cout);
    return mphys;
}
static TPZCompMesh* CompMeshCube(TPZGeoMesh* gmesh)
{
    MechParamsVonMises param;
    auto *mphys = new TPZCompMesh(gmesh);
    mphys->SetDimModel(3);
    mphys->SetAllCreateFunctionsContinuousWithMem();
    mphys->SetDefaultOrder(2);
    auto mat = new TMatElastoPlaticVoigtVM3D(1);


    TPZElasticResponse ER;
    ER.SetEngineeringData(200000., 0.0);
    TPZYCVonMisesVoigt vmyc;

    const STATE sigmaY0 = 200.0;
    const STATE Hiso    = 10000.;
    vmyc.SetUp(sigmaY0,Hiso);

    TPlasticStepVoigtVM PlasticStepVoigt;


    PlasticStepVoigt.SetPlasticCriterion(vmyc);
    PlasticStepVoigt.SetElasticResponse(ER);


    PlasticStepVoigt.Print(std::cout);

    mat->SetPlasticityModel(PlasticStepVoigt);
    mat->SetId(1);
    mphys->InsertMaterialObject(mat);//0
    //mat->Print(std::cout);

    TPZFMatrix<STATE> v1(3,3,0.);
    TPZManVector<STATE,3> v2(3,0.);
    int dirdirichlet=3,pressure=5,newmann=1;

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
    mphys->InsertMaterialObject(mat->CreateBC(mat,bcindexes.bcright, dirdirichlet, v1, v2));

    v2[0]=0.;
    v2[1]=0.;
    v2[2]=0;
    mphys->InsertMaterialObject(mat->CreateBC(mat,bcindexes.bctop , newmann, v1, v2));
/*
    v2[0]=0.;
    v2[1]=0.;
    v2[2]=0.;
    mphys->InsertMaterialObject(mat->CreateBC(mat,bcindexes.bctop , 0, v1, v2));
*/

    mphys->AutoBuild();
    mphys->AdjustBoundaryElements();
    mphys->CleanUpUnconnectedNodes();
    //mphys->Print(cout);
    return mphys;
}


void SolveCyl()
{
    auto gmesh = PressurizedCylinderMesh();
    auto cmesh = CompMeshCyl(gmesh);

    TPZElastoPlasticAnalysis an(cmesh, std::cout,TPZElastoPlasticAnalysis::ELineSearch::Armijo);

    TPZFStructMatrix<STATE> str(cmesh);
    an.SetStructuralMatrix(str);
    TPZStepSolver<REAL> direct;
    direct.SetDirect(ELU);
    an.SetSolver(direct);
    int itersout;

    int nloads = 10;
    const REAL FS_target = 1. ;          // mantém seu “+0.1”
    TPZManVector<REAL> factors(nloads+1);   // 0 .. nloads (inclusivo)
    for (int i =0; i <= nloads; ++i) {
        factors[i] = FS_target * REAL(i) / REAL(nloads); // 0, Δ, 2Δ, …, FS_target
        cout<< factors[i] <<endl;
    }
    IndexesCylinder cylindexes;
    int loaddir=0;//direcao da pressao
    ApplyLoad( cmesh,factors,loaddir,cylindexes.bcinner);
}




void ApplyLoad(TPZCompMesh* cmesh,TPZManVector<REAL> factors,int loaddir,int indexbc)
{
    int dim=cmesh->Dimension();

    cout << "dimensao da malha computacional = " << dim << endl;

    // parâmetros de controle
    int nloads = factors.size();

    TPZElastoPlasticAnalysis anal(cmesh, std::cout,TPZElastoPlasticAnalysis::ELineSearch::Armijo);

    auto* bcmat = dynamic_cast<TPZBndCondT<STATE>*>(cmesh->FindMaterial(indexbc));

    if(!bcmat)
    {
        std::cout << "material do contorno nao encontrado"<<std::endl;
        DebugStop();

    }
    const REAL load0=bcmat->Val2()[loaddir];
    cout << " load0 = "<< load0 <<endl;
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


    int matid=1;
    std::string vtkfile="cyl.vtk";
    for(int i =0;i< nloads;i++)
    {

        bcmat->Val2()[loaddir]=load0*factors[i];
       // bcmat->Val2()[loaddir]=factors[i];
        cout << "Load step =" << i<<" factor =  "<<factors[i] <<endl;
        int iters_out;

        REAL resf,resuu;
        //bool ok = anal.FindRoot( iters_out,resf,resuu);
       bool ok = anal.IterativeProcess(std::cout, 1.e-6,100, true, false, iters_out);

        TPZFMatrix<REAL> tempsol=anal.Solution();
        bcmat->Val2()[loaddir]=0;
        anal.AcceptSolution(0);

        // cmesh->LoadSolution(anal.CumulativeSolution());
        PostElastoplastic(cmesh,vtkfile,matid,i,dim);

        //tempsol.Zero();
        //cmesh->LoadSolution(tempsol);
        //anal.LoadSolution();
    }

}
void SolveCube()
{
    auto gmesh = CubeMesh();
    auto cmesh = CompMeshCube(gmesh);



    TPZElastoPlasticAnalysis an(cmesh, std::cout,TPZElastoPlasticAnalysis::ELineSearch::Armijo);

    TPZFStructMatrix<STATE> str(cmesh);
    an.SetStructuralMatrix(str);
    TPZStepSolver<REAL> direct;
    direct.SetDirect(ELU);
    an.SetSolver(direct);

    //TPZManVector<REAL> factors={0,-0.0005,-0.0005,-0.0005,-0.0005,-0.0001,0,0,0,0.001,0.001,0.001};

    BoundaryIndexes cubeindexes;
    int loaddir=2;
    int nsteps=100;
    int nloadcicles=3;
    STATE bccondval=200.;
    STATE lambda0=0.0001;
    STATE L0=0.0003;
    std::string vtkfile="cubeal2.vtk";

    //IterativeProcessArcLength(an,cmesh, loaddir, cubeindexes.bctop,vtkfile);

    IterativeProcessArcLength(an,loaddir,cubeindexes.bctop,nsteps,nloadcicles,bccondval,lambda0,L0,vtkfile);

    //IterativeProcessArcLength(an ,cmesh, loaddir, cubeindexes.bctop,vtkfile);
    //TPZManVector<REAL> factors={200,150,100,50,0};
    //ApplyLoad( cmesh,factors,loaddir,cubeindexes.bctop);
}
int main()
{
    //const std::string configfile = "/home/diogo/projects/neopz-master-build-debug/Util/log4cxx.cfg";
    //TPZLogger::InitializePZLOG(configfile);

    //SolveCyl();

    SolveCube();

    return 0;
}
void PostProcessVariables(TPZStack<std::string>& scal, TPZStack<std::string>& vec)
{
   scal.Push ( "StrainPlasticJ2" );
    vec.Push ( "Displacement" );
    //scal.Push ( "StressXX" );
    //scal.Push ( "StressYY" );
    //scal.Push ( "StrainElasticJ2" );
    scal.Push ( "StressZZ" );
    //scal.Push ( "StrainPlasticXX" );
    //scal.Push ( "StrainPlasticYY" );
    scal.Push ( "StrainPlasticZZ" );
    scal.Push ( "StrainTotalZZ" );
    //scal.Push ( "StrainElasticZZ" );
    scal.Push ( "DamageVariable" );
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

REAL UyAtNode3D(TPZCompMesh* cmesh, REAL x, REAL y,REAL z)
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

// ===== utilidades simples =====
static inline STATE Dot(const TPZFMatrix<STATE>& a, const TPZFMatrix<STATE>& b) {
    TPZFMatrix<STATE> at,temp;
    a.Transpose(&at);
    at.Multiply(b,temp);
    if(temp.Rows()>1)DebugStop();
    STATE s=temp(0,0);
    return s;
}

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


REAL IterativeProcessArcLength(TPZElastoPlasticAnalysis &an,TPZCompMesh* cmesh,int loaddir,int indexbc,std::string vtkfile)
{

    auto* bcmat = dynamic_cast<TPZBndCondT<STATE>*>(cmesh->FindMaterial(indexbc));
    if (!bcmat) {
        std::cout << "[ArcLength] BndCond material não encontrado (indexbc="
                  << indexbc << ")\n";
        DebugStop();
    }
    std::cout << "\n========== CICLO DE CARGA " << "\n";
    int nsteps=100;
    // ===========================================================
    //           LOOP DE CICLOS DE CARGA (2 ciclos, por ex.)
    // ===========================================================
    for (int iloadcicle = 0; iloadcicle < 2; iloadcicle++) {

        std::cout << "\n========== CICLO DE CARGA " << (iloadcicle + 1)
                  << " ==========\n";

        int dim = cmesh->Dimension();
        const int maxit_inner = 20;
        const REAL etol_inner = 1.e-3;
        REAL lambda = 0.0001;
        REAL lambdan = lambda;
        REAL L = 0.0003;

        TPZFMatrix<STATE> u_acc = an.Solution(); u_acc.Zero();
        TPZFMatrix<STATE> dw = an.Solution(); dw.Zero();

        const STATE load0 = 200.0;
        bcmat->Val2()[loaddir] = load0;
        an.Assemble();
        TPZFMatrix<STATE> FEXT = an.Rhs();
        const REAL normFEXT = Norm(FEXT);
        bcmat->Val2()[loaddir] = 0.0;

        int step = 0;
        bool okconv = false;
        REAL diff = 1e9;

        // =======================================================
        //      Loop do método de comprimento de arco (Riks)
        // =======================================================
        const int max_cuts = 12;   // limite para reduzir L no mesmo step
        int cut_count = 0;

        while (step < nsteps) {

            // Estado de trabalho parte do aceito
            TPZFMatrix<STATE> u  = u_acc;
            TPZFMatrix<STATE> dw = an.Solution(); dw.Zero();   // incremento acumulado do step
            bool okconv = false;
            REAL normR  = 1e30;

            // ---------- NEWTON DE COMPRIMENTO DE ARCO ----------
            {
                const int  maxit_inner = 20;
                const REAL etol_inner  = 1e-3;

                TPZFMatrix<STATE> rhs, R, dws, dwb;

                int it = 0;
                while (it < maxit_inner && normR > etol_inner)
                {
                    // monta -FINT
                    bcmat->Val2()[loaddir] = 0.0;
                    an.LoadSolution(u);
                    an.Assemble();
                    rhs = an.Rhs(); // -FINT

                    // R = lambda * FEXT + (-FINT)
                    R = FEXT * lambda;
                    R += rhs;

                    // solve para dws
                    an.Rhs() = R;
                    an.Solve();
                    dws = an.Solution();

                    // solve para dwb (direção básica)
                    an.Rhs() = FEXT;
                    an.Solve();
                    dwb = an.Solution();

                    // atualiza dl pelo Riks
                    REAL dl = (it == 0) ? compute_dlambda0_riks(dwb, dw, L)
                    : compute_dlambda_riks(dwb, dws, dw, L, /*rootIdx*/ *(int[]){0});

                    const REAL dl_max = 0.1;
                    if (dl >  dl_max) dl =  dl_max;
                    if (dl < -dl_max) dl = -dl_max;

                    TPZFMatrix<STATE> dwtot = dwb*dl + dws;

                    u      += dwtot;
                    dw     += dwtot;
                    lambda += dl;

                    // remonta resíduo para checagem
                    an.LoadSolution(u);
                    bcmat->Val2()[loaddir] = 0.0;
                    an.Assemble();
                    rhs   = an.Rhs();          // -FINT
                    R     = FEXT * lambda;     // lambda*FEXT
                    R    += rhs;               // + (-FINT)
                    normR = Norm(R);

                    ++it;
                    if (normR > 1e5) break;   // guarda de divergência
                }

                okconv = (normR <= etol_inner);
            }
            // ---------- FIM NEWTON ----------

            if (okconv) {
                // aceita o passo
                an.AcceptSolution();
                u_acc   = an.Solution();
                PostElastoplastic(cmesh, vtkfile, /*matid*/1, /*out_step*/ step + iloadcicle*nsteps, cmesh->Dimension());

                // pronto para o próximo step
                ++step;
                cut_count = 0;
            } else {
                // rollback e corta L
                an.LoadSolution(u_acc);
                lambda = lambdan;   // volta lambda aceito
                L *= 0.5;
                ++cut_count;

                if (cut_count > max_cuts || L < 1e-14) {
                    std::cout << "[ArcLength] Falha em convergir no step " << step
                    << " após " << cut_count << " cortes de L. Abortando ciclo.\n";
                    break; // evita loop infinito
                }

                // tenta novamente o MESMO step com L menor
                continue;
            }

            // opcional: atualiza “lambda aceito” para rollback seguro
            lambdan = lambda;
        }

        // =======================================================
        //            DESCARREGAMENTO (unloading)
        // =======================================================
        std::cout << "Unloading ciclo " << (iloadcicle + 1) << std::endl;
        TPZManVector<REAL> factors = {
            lambda * load0,
            lambda * load0 * 0.8,
            lambda * load0 * 0.6,
            lambda * load0 * 0.4,
            lambda * load0 * 0.2,
            0.0
        };

        for (int i = 0; i < factors.size(); i++) {
            bcmat->Val2()[loaddir] = factors[i];
            std::cout << "Load step = " << i<< " factor = " << factors[i] << std::endl;
            int iters_out;
            bool ok = an.IterativeProcess(std::cout, 1.e-3, 30, true, false, iters_out);
            bcmat->Val2()[loaddir] = 0.0;
            an.AcceptSolution();
            PostElastoplastic(cmesh, vtkfile, 1, i + step + iloadcicle * nsteps, dim);
        }
    }

    an.AcceptSolution();
    return 0.0;
}

// Correções e documentação do método de Comprimento de Arco (Riks)
// Autor da correção: copilot (@copilot)
// Observação: este ficheiro contém apenas a função corrigida e documentada.
// Tip: adapte includes e nomes de tipos conforme o seu projeto (TPZElastoPlasticAnalysis, TPZFMatrix, TPZManVector, etc).

REAL IterativeProcessArcLength(TPZElastoPlasticAnalysis &an,
                               int loaddir,
                               int indexbc,
                               int nsteps,
                               int nloadcicles,
                               STATE bccondval,
                               STATE lambda0,
                               STATE L0,
                               std::string vtkfile)
{
    // Recupera a malha/condição de contorno
    auto cmesh = an.Mesh();
    auto* bcmat = dynamic_cast<TPZBndCondT<STATE>*>(cmesh->FindMaterial(indexbc));
    if (!bcmat) {
        std::cout << "[ArcLength] BndCond material não encontrado (indexbc="
        << indexbc << ")\n";
        DebugStop();
    }

    // ===========================================================
    //           LOOP DE CICLOS DE CARGA
    // ===========================================================
    for (int iloadcicle = 0; iloadcicle < nloadcicles; iloadcicle++) {

        std::cout << "\n========== CICLO DE CARGA " << (iloadcicle + 1)
        << " ==========\n";

        int dim = cmesh->Dimension();
        REAL lambda = lambda0;   // lambda corrente (pode ser alterado durante os passos)
        REAL lambdan = lambda;   // lambda aceito (para rollback seguro)
        REAL L = L0;             // comprimento de arco atual

        // Solução acumulada aceita (u_acc) e incremento acumulado do passo (dw)
        TPZFMatrix<STATE> u_acc = an.Solution();
        u_acc.Zero();

        // Monta vetor de cargas externas (FEXT) com a condição de contorno aplicada
        const STATE load0 = bccondval;
        bcmat->Val2()[loaddir] = load0;
        an.Assemble();
        TPZFMatrix<STATE> FEXT = an.Rhs();
        // remove o carregamento aplicado da matriz de cond. de contorno,
        // pois o termo externo será aplicado multiplicado por lambda
        bcmat->Val2()[loaddir] = 0.0;

        int step = 0;

        // =======================================================
        //      Loop do método de comprimento de arco (Riks)
        // =======================================================
        const int max_cuts = 12;   // limite para reduzir L no mesmo step
        int cut_count = 0;

        while (step < nsteps) {

            // Estado de trabalho parte do aceito
            TPZFMatrix<STATE> u  = u_acc;
            // NOTA: NÃO redeclarar dw aqui (shadowing). Usamos dw como incremento do passo.
            // bool de convergência e norma do resíduo
            bool okconv = false;
            REAL normR  = 1e30;
            TPZFMatrix<STATE> dw = an.Solution();
            dw.Zero();

            // ---------- NEWTON DE COMPRIMENTO DE ARCO ----------
            {
                const int  maxit_inner = 20;
                const REAL etol_inner  = 1e-3;

                TPZFMatrix<STATE> rhs, R, dws, dwb;

                int it = 0;
                // Mantemos um índice de "root" para funções auxiliares de escolha de raiz,
                // se a implementação de compute_dlambda_riks exigir.
                int rootIdx = 0;

                while (it < maxit_inner && normR > etol_inner)
                {
                    // monta -FINT (observação: usamos bcmat->Val2()[loaddir]=0 para evitar
                    // re-aplicar o BC na montagem do interno; o externo será lambda*FEXT)
                    bcmat->Val2()[loaddir] = 0.0;
                    an.LoadSolution(u);
                    an.Assemble();
                    rhs = an.Rhs(); // -FINT

                    // R = lambda * FEXT + (-FINT)
                    R = FEXT * lambda;
                    R += rhs;

                    // resolve para dws (incremento devido ao passo de Newton)
                    an.Rhs() = R;
                    an.Solve();
                    dws = an.Solution();

                    // resolve para dwb (direção básica / carregamento)
                    an.Rhs() = FEXT;
                    an.Solve();
                    dwb = an.Solution();

                    // atualiza dl pelo Riks
                    REAL dl = (it == 0) ? compute_dlambda0_riks(dwb, dw, L)
                    : compute_dlambda_riks(dwb, dws, dw, L, /*rootIdx*/ *(int[]){0});

                    const REAL dl_max = 0.1;
                    if (dl >  dl_max) dl =  dl_max;
                    if (dl < -dl_max) dl = -dl_max;

                    TPZFMatrix<STATE> dwtot = dwb*dl + dws;

                    u      += dwtot;
                    dw     += dwtot;
                    lambda += dl;

                    // remonta resíduo para checagem
                    an.LoadSolution(u);
                    bcmat->Val2()[loaddir] = 0.0;
                    an.Assemble();
                    rhs   = an.Rhs();          // -FINT
                    R     = FEXT * lambda;     // lambda*FEXT
                    R    += rhs;               // + (-FINT)
                    normR = Norm(R);

                    ++it;
                    // guarda de divergência: se a norma explode, aborta o inner
                    if (normR > 1e5) {
                        break;
                    }
                } // fim do while (Newton)

                okconv = (normR <= etol_inner);
            }
            // ---------- FIM NEWTON ----------

            if (okconv) {
                // aceita o passo
                an.AcceptSolution();
                u_acc   = an.Solution();

                // salva saída VTK (função externa no seu projeto)
                PostElastoplastic(cmesh, vtkfile, /*matid*/1, /*out_step*/ step + iloadcicle*nsteps, dim);

                // pronto para o próximo step
                ++step;
                cut_count = 0;

                // atualiza “lambda aceito” para rollback seguro
                lambdan = lambda;
            } else {
                // rollback e corta L
                an.LoadSolution(u_acc);
                lambda = lambdan;   // volta lambda aceito
                L *= 0.5;
                ++cut_count;

                if (cut_count > max_cuts || L < 1e-14) {
                    std::cout << "[ArcLength] Falha em convergir no step " << step
                    << " após " << cut_count << " cortes de L. Abortando ciclo.\n";
                    break; // evita loop infinito do ciclo
                }

                // tenta novamente o MESMO step com L menor
                continue;
            }

        } // fim while step < nsteps

        // =======================================================
        //            DESCARREGAMENTO (unloading)
        // =======================================================
        std::cout << "Unloading ciclo " << (iloadcicle + 1) << std::endl;
        TPZManVector<REAL> factors = {
            lambda * load0,
            lambda * load0 * 0.8,
            lambda * load0 * 0.6,
            lambda * load0 * 0.4,
            lambda * load0 * 0.2,
            0.0
        };

        // usar int para compatibilidade com TPZManVector::size() (depende da versão)
        int nFactors = static_cast<int>(factors.size());
        for (int i = 0; i < nFactors; i++) {
            bcmat->Val2()[loaddir] = factors[i];
            std::cout << "Load step = " << i << " factor = " << factors[i] << std::endl;
            int iters_out = 0;
            bool ok = an.IterativeProcess(std::cout, 1.e-3, 30, true, false, iters_out);
            bcmat->Val2()[loaddir] = 0.0;
            an.AcceptSolution();
            PostElastoplastic(cmesh, vtkfile, 1, i + step + iloadcicle * nsteps, dim);
        }
    } // fim ciclos de carga

    an.AcceptSolution();
    return 0.0;
}

// REAL IterativeProcessArcLength(TPZElastoPlasticAnalysis &an, TPZCompMesh* cmesh, int loaddir, int indexbc,std::string vtkfile)
// {
//
//
//     WriteHeaderCSV("saida_arc.csv");
//
//
//     auto* bcmat = dynamic_cast<TPZBndCondT<STATE>*>(cmesh->FindMaterial(indexbc));
//     if (!bcmat) {
//         std::cout << "[ArcLength] BndCond material não encontrado (indexbc=" << indexbc << ")\n";
//         DebugStop();
//     }
//     int globalstep=0;
//     bool okconv = false;
//     REAL lambda            = 0.0001;       // λ inicial
//     for(int iloadcicle=0;iloadcicle<2;iloadcicle++)
//     {
//     int dim =cmesh->Dimension();
//     // Parâmetros
//     const int  maxit_inner = 20;
//     const REAL etol_inner  = 1.e-3;
//     REAL lambdan=lambda;
//     REAL L                 = 0.0003;      // alvo de arco (norma do incremento total em u)
//
//     // Estado aceito acumulado
//     TPZFMatrix<STATE> u_acc = an.Solution(); u_acc.Zero();
//     // Incremento acumulado no passo corrente
//     TPZFMatrix<STATE> dw = an.Solution();   dw.Zero();
//
//     STATE load0=-200;
//     // ---- Monta FEXT com uma BC "unitária" (ou escala desejada) ----
//     bcmat->Val2()[loaddir] = load0;   // valor de referência para extrair o vetor de cargas
//     an.Assemble();
//     TPZFMatrix<STATE> FEXT = an.Rhs(); // isto é o "rhs" gerado pela BC -> aqui use como vetor de cargas
//     REAL normFEXT = Norm(FEXT);
//     std::cout << " Norm(FEXT) = " << normFEXT << std::endl;
//     bcmat->Val2()[loaddir] = 0.0;     // zera a BC mecânica para o loop de Newton
//     int matid=1;
//
//     // ---- Um único passo de arco (ajuste se quiser mais) ----
//     int  step   = 0;
//
//     STATE diff=10000;
//     STATE tollamb=0.001;
//     while (step < 100)
//     {
//         std::cout << "\n[Arc step " << (step+1) << "]  L=" << L << "  lambda=" << lambda << std::endl;
//
//         // Estado de trabalho (começa no aceito)
//         TPZFMatrix<STATE> u  = u_acc;
//         TPZFMatrix<STATE> du = an.Solution();
//         du.Zero(); // incremento total do passo
//
//         int  it    = 0;
//         REAL normR = 1e30;
//
//         TPZFMatrix<STATE> rhs, R, dws, dwb; // declara fora
//         dw.Zero();
//
//
//         an.LoadSolution(u);
//         while (it < maxit_inner && normR > etol_inner)
//         {
//             rhs.Zero(); R.Zero(); dws.Zero(); dwb.Zero();
//             bcmat->Val2()[loaddir] = 0.0;
//             an.Assemble();
//             rhs = an.Rhs();   // rhs = -FINT
//
//             TPZFMatrix<STATE> R = FEXT * lambda;
//             R += rhs;
//
//             an.Rhs() = R;
//             an.Solve();
//             dws = an.Solution();
//
//             an.Rhs() = FEXT;
//             an.Solve();
//             TPZFMatrix<STATE> dwb = an.Solution();
//             REAL normdws = Norm(dws);
//             REAL normdwb = Norm(dwb);
//
//             REAL dl = 0.0; int rootIdx = 0;
//             if (it == 0) {
//                 dl = compute_dlambda0_riks(dwb, dw, L);
//             } else {
//                 dl = compute_dlambda_riks(dwb, dws, dw, L, rootIdx);
//             }
//             const REAL dl_max = 0.1;        // ajuste conforme seu problema
//             if (dl >  dl_max) dl =  dl_max;
//             if (dl < -dl_max) dl = -dl_max;
//             TPZFMatrix<STATE> dwtot = dwb*dl + dws;
//             u      += dwtot;
//             dw     += dwtot;
//             lambda += dl;
//
//             REAL normdu = Norm(dwtot);
//
//             an.LoadSolution(u);
//
//
//             bcmat->Val2()[loaddir] = 0.0;
//             an.Assemble();
//             rhs   = an.Rhs();
//             R     = FEXT * lambda;
//             R += rhs;
//             normR = Norm(R);
//
//
//             // cout << "  [it " << it
//             // << "] diff=" << diff
//             // << " lambda=" << lambda
//             // << "] ||dws||=" << normdws
//             // << " ||dwb||=" << normdwb
//             // << " dl=" << dl
//             //
//             // << " ||dwtot||=" << normdu
//             // << " \n ||R||=" << normR
//             // //<< " ||FINT||=" << Norm(FINT)
//             // << " lambda||FEXT||=" << lambda*Norm(FEXT)
//             // << " ||FEXT||=" << Norm(FEXT)
//             // << endl;
//             it++;
//             if(normR>1.e5)break;
//         }
//
//
//         okconv = (normR <= etol_inner);
//
//         if(okconv)
//         {
//             std::cout << " (OK  ) — posprocessa, atualiza lambda e u \n";
//
//             const double uz = UyAtNode3D(cmesh, 0.5, 0.5, 1.0);
//             const double load = lambda * normFEXT;
//             AppendRowCSV("saida_arc.csv", globalstep, lambda, load, uz);
//
//             an.AcceptSolution();
//             PostElastoplastic(cmesh, vtkfile, matid, globalstep, dim);
//             diff=fabs(lambda-lambdan);
//             cout << "  [it " << it
//             << "] diff=" << diff <<"\n";
//             lambdan = lambda;
//             u_acc   = an.Solution();
//             const int ndesi = 5;
//             //L *= REAL(ndesi) / std::max(1, it);
//             //if (L > 1) L = 1;
//                              // refaz o MESMO step com L menor
//         }
//         else
//         {
//             std::cout << " (não convergiu) rollback e reduzir L\n" << std::endl;
//             lambda = lambdan;
//             an.LoadSolution(u_acc);
//             dw.Zero();
//             L *= 0.5;
//             continue;
//         }
//
//         step++;
//         globalstep++;
//     }
//     cout << "Unloading" <<endl;
//     TPZManVector<REAL> factors={lambda*load0,lambda*load0*0.8,lambda*load0*0.6,lambda*load0*0.4,lambda*load0*0.2,lambda*load0*0.0};
//     for(int i =0;i< factors.size();i++)
//     {
//         bcmat->Val2()[loaddir]=factors[i];
//         cout << "Load step =" << i<<" factor =  "<<factors[i] <<endl;
//         int iters_out;
//         bool ok = an.IterativeProcess(std::cout, 1.e-3,30, true, false, iters_out);
//         const double uz = UyAtNode3D(cmesh, 0.5, 0.5, 1.0);
//         const double load = factors[i];
//         AppendRowCSV("saida_arc.csv", globalstep, lambda, load, uz);
//         TPZFMatrix<REAL> tempsol=an.Solution();
//         bcmat->Val2()[loaddir]=0;
//         an.AcceptSolution();
//
//         PostElastoplastic(cmesh,vtkfile,matid,i+globalstep,dim);
//
//     }
//
//     }
//     an.AcceptSolution();
//     return okconv ? lambda : 0.;
// }


// REAL IterativeProcessArcLength(TPZCompMesh* cmesh, int loaddir, int indexbc)
// {
//     TPZElastoPlasticAnalysis an(
//         cmesh, std::cout,
//         TPZElastoPlasticAnalysis::ELineSearch::Dicotomic
//     );
//
//     auto* bcmat = dynamic_cast<TPZBndCondT<STATE>*>(cmesh->FindMaterial(indexbc));
//     if (!bcmat) {
//         std::cout << "[ArcLength] BndCond material não encontrado (indexbc=" << indexbc << ")\n";
//         DebugStop();
//     }
//
//     int dim =cmesh->Dimension();
//     // Parâmetros
//     const int  maxit_inner = 20;
//     const REAL etol_inner  = 1.e-3;
//     REAL lambda            = 0.5;       // λ inicial
//     REAL lambdan=lambda;
//     REAL L                 = 0.01;      // alvo de arco (norma do incremento total em u)
//
//     // Estado aceito acumulado
//     TPZFMatrix<STATE> u_acc = an.Solution(); u_acc.Zero();
//     // Incremento acumulado no passo corrente
//     TPZFMatrix<STATE> dw = an.Solution();   dw.Zero();
//
//     // ---- Monta FEXT com uma BC "unitária" (ou escala desejada) ----
//     bcmat->Val2()[loaddir] = 200;   // valor de referência para extrair o vetor de cargas
//     an.Assemble();
//     TPZFMatrix<STATE> FEXT = an.Rhs(); // isto é o "rhs" gerado pela BC -> aqui use como vetor de cargas
//     REAL normFEXT = Norm(FEXT);
//     std::cout << " Norm(FEXT) = " << normFEXT << std::endl;
//     bcmat->Val2()[loaddir] = 0.0;     // zera a BC mecânica para o loop de Newton
//     int matid=1;
//     std::string vtkfile="cubeal.vtk";
//     // ---- Um único passo de arco (ajuste se quiser mais) ----
//     int  step   = 0;
//     bool okconv = false;
//     while (step < 20)
//     {
//         std::cout << "\n[Arc step " << (step+1) << "]  L=" << L << "  lambda=" << lambda << std::endl;
//
//         // Estado de trabalho (começa no aceito)
//         TPZFMatrix<STATE> u  = u_acc;
//         TPZFMatrix<STATE> du = an.Solution();
//         du.Zero(); // incremento total do passo
//
//         int  it    = 0;
//         REAL normR = 1e30;
//
//         TPZFMatrix<STATE> rhs, R, dws, dwb; // declara fora
//         dw.Zero();
//
//
//         an.LoadSolution(u);
//         while (it < maxit_inner && normR > etol_inner)
//         {
//             rhs.Zero(); R.Zero(); dws.Zero(); dwb.Zero();
//             bcmat->Val2()[loaddir] = 0.0;
//             an.Assemble();
//             rhs = an.Rhs();   // rhs = -FINT
//
//             // --- resíduo físico: R = λ FEXT + rhs  (pois rhs = -FINT) ---
//             TPZFMatrix<STATE> R = FEXT * lambda;
//             R += rhs;
//
//             // --- PARTICULAR: K dws = R ---
//             an.Rhs() = R;
//             an.Solve();
//             dws = an.Solution();
//
//             // --- DIREÇÃO EM λ: K dwb = FEXT ---
//             an.Rhs() = FEXT;
//             an.Solve();
//             TPZFMatrix<STATE> dwb = an.Solution();
//             REAL normdws = Norm(dws);
//             REAL normdwb = Norm(dwb);
//             // ---------- cálculo do dl (com salvaguardas) ----------
//             REAL dl = 0.0; int rootIdx = 0;
//             if (it == 0) {
//                 // “curvature safeguard” no preditor
//                 static REAL sgn = 1.0;
//                 if (Dot(dwb, dw) < 0.0) sgn = -sgn;
//                 //dl = sgn * std::fabs(compute_dlambda0_riks(dwb, dw, L));
//                 dl = compute_dlambda0_riks(dwb, dw, L);
//             } else {
//                 dl = compute_dlambda_riks(dwb, dws, dw, L, rootIdx);
//             }
//
//             // clamp em dl
//             const REAL dl_max = 0.1;        // ajuste conforme seu problema
//             if (dl >  dl_max) dl =  dl_max;
//             if (dl < -dl_max) dl = -dl_max;
//             TPZFMatrix<STATE> dwtot = dwb*dl + dws;
//             u      += dwtot;
//             dw     += dwtot;
//             lambda += dl;
//
//             REAL normdu = Norm(dwtot);
//
//             an.LoadSolution(u);
//
//             // ---------- reavaliar resíduo ATUALIZADO (ESSENCIAL!) ----------
//             bcmat->Val2()[loaddir] = 0.0;
//             an.Assemble();
//             rhs   = an.Rhs();                 // -FINT(u atualizado)
//             R     = FEXT * lambda;   R += rhs;
//             normR = Norm(R);
//
//
//             cout << "  [it " << it
//             << "] ||dws||=" << normdws
//             << " ||dwb||=" << normdwb
//             << " dl=" << dl
//             << " lambda=" << lambda
//             << " ||dwtot||=" << normdu
//             << " \n ||R||=" << normR
//             //<< " ||FINT||=" << Norm(FINT)
//             << " lambda||FEXT||=" << lambda*Norm(FEXT)
//             << " ||FEXT||=" << Norm(FEXT)
//             << endl;
//             it++;
//             if(normR>1.e5)break;
//         }
//
//
//         okconv = (normR <= etol_inner);
//
//         if(okconv)
//         {
//             REAL uz = UyAtNode3D(cmesh, 0.5,0.5,1.);
//             std::cout << uz << " " << lambda * normFEXT << std::endl;
//             std::cout << " (convergiu)\n" << std::endl;
//
//             // -------- DUMMY NEWTON CHECK (não persiste) --------
//             // snapshot
//             TPZFMatrix<STATE> u_snap = an.Solution();
//             const REAL bc_old = bcmat->Val2()[loaddir];
//
//             // use exatamente a mesma escala que você usou para montar FEXT
//             const REAL load_ref = 200.0;                 // <- se FEXT veio de Val2()=100
//             bcmat->Val2()[loaddir] = load_ref * lambda;  // λ fixo no dummy
//
//             int iters_out = 0;
//             bool ok_dummy = an.IterativeProcess(std::cout, 1.e-3, 10,
//                                                 /*line-search*/ false,
//                                                 /*update_mem*/  false,                      // <— NÃO atualiza memória
//                                                 iters_out);
//
//             // restaura tudo (dummy não deixa rastro)
//             an.LoadSolution(u_snap);
//             bcmat->Val2()[loaddir] = bc_old;
//
//             // -------- FIM DUMMY --------
//
//             if(ok_dummy)
//             {
//                 std::cout << " (OK dummy ) — posprocessa, atualiza lambda e u \n";
//                 // agora sim, persiste o passo convergido
//                 an.AcceptSolution();                       // variante correta
//                 PostElastoplastic(cmesh, vtkfile, matid, step, dim);
//                 lambdan = lambda;                           // salva λ aceito
//                 u_acc   = an.Solution();                    // sincroniza com o que o analysis tem
//                 // (se quiser adaptar L aqui, ok)
//                 const int ndesi = 5;
//                 L *= REAL(ndesi) / std::max(1, it);
//                 if (L > 1) L = 1;
//                 //if (L > 1) L = 1;
//             }
//             else
//             {
//                 std::cout << " (dummy não convergiu) — rollback e reduzir L\n";
//                 // rollback para último estado aceito
//                 lambda = lambdan;
//                 an.LoadSolution(u_acc);
//                 dw.Zero();
//                 L *= 0.5;
//                 continue;                                   // refaz o MESMO step com L menor
//             }
//         }
//         else
//         {
//             std::cout << " (não convergiu) rollback e reduzir L\n" << std::endl;
//
//             // rollback para último aceito e repetir o step com L menor
//             lambda = lambdan;
//             an.LoadSolution(u_acc);
//             dw.Zero();
//             L *= 0.5;
//             continue;                                       // refaz o MESMO step
//         }
//
//
//
//
//         step++;
//     }
//     an.AcceptSolution();
//     return okconv ? lambda : 0.;
// }


// REAL IterativeProcessArcLength(TPZCompMesh* cmesh, int loaddir, int indexbc)
// {
//     TPZElastoPlasticAnalysis an(
//         cmesh, std::cout,
//         TPZElastoPlasticAnalysis::ELineSearch::Armijo
//     );
//
//     auto* bcmat = dynamic_cast<TPZBndCondT<STATE>*>(cmesh->FindMaterial(indexbc));
//     if (!bcmat) {
//         std::cout << "[ArcLength] BndCond material não encontrado (indexbc=" << indexbc << ")\n";
//         DebugStop();
//     }
//
//     int dim =cmesh->Dimension();
//     // Parâmetros
//     const int  maxit_inner = 10;
//     const REAL etol_inner  = 1.e-3;
//     REAL lambda            = 0.01;       // λ inicial
//     REAL lambdan=lambda;
//     REAL L                 = 1.e-4;      // alvo de arco (norma do incremento total em u)
//
//     // Estado aceito acumulado
//     TPZFMatrix<STATE> u_acc = an.Solution(); u_acc.Zero();
//     // Incremento acumulado no passo corrente
//     TPZFMatrix<STATE> dw = an.Solution();   dw.Zero();
//
//     // ---- Monta FEXT com uma BC "unitária" (ou escala desejada) ----
//     bcmat->Val2()[loaddir] = 100;   // valor de referência para extrair o vetor de cargas
//     an.Assemble();
//     TPZFMatrix<STATE> FEXT = an.Rhs(); // isto é o "rhs" gerado pela BC -> aqui use como vetor de cargas
//     REAL normFEXT = Norm(FEXT);
//     std::cout << " Norm(FEXT) = " << normFEXT << std::endl;
//     bcmat->Val2()[loaddir] = 0.0;     // zera a BC mecânica para o loop de Newton
//     int matid=1;
//     std::string vtkfile="cubeal.vtk";
//     // ---- Um único passo de arco (ajuste se quiser mais) ----
//     int  step   = 0;
//     bool okconv = false;
//     while (step < 20)
//     {
//         std::cout << "\n[Arc step " << (step+1) << "]  L=" << L << "  lambda=" << lambda << std::endl;
//
//         // Estado de trabalho (começa no aceito)
//         TPZFMatrix<STATE> u  = u_acc;
//         TPZFMatrix<STATE> du = an.Solution();
//         du.Zero(); // incremento total do passo
//
//         int  it    = 0;
//         REAL normR = 1e30;
//
//         dw.Zero();
//
//         while (it < maxit_inner && normR > etol_inner)
//         {
//             bcmat->Val2()[loaddir] = 0.0;
//             an.Assemble();
//             TPZFMatrix<STATE> rhs = an.Rhs();   // rhs = -FINT
//
//             // --- resíduo físico: R = λ FEXT + rhs  (pois rhs = -FINT) ---
//             TPZFMatrix<STATE> R = FEXT * lambda;
//             R += rhs;
//
//             // --- PARTICULAR: K dws = R ---
//             an.Rhs() = R;
//             an.Solve();
//             TPZFMatrix<STATE> dws = an.Solution();
//
//             // --- DIREÇÃO EM λ: K dwb = FEXT ---
//             an.Rhs() = FEXT;
//             an.Solve();
//             TPZFMatrix<STATE> dwb = an.Solution();
//
//             const REAL normdws = Norm(dws);
//             const REAL normdwb = Norm(dwb);
//
//             // ---------- coeficientes da restrição g(dl)=a dl^2 + b dl + c ----------
//             TPZFMatrix<STATE> t = dw;  t += dws;     // t = dw + dws
//             long double a_ld = (long double)Dot(dwb,dwb);
//             long double b_ld = 2.0L * (long double)Dot(dwb,t);
//             long double c_ld = (long double)Dot(t,t) - (long double)L*(long double)L;
//             auto g_eval = [&](REAL dl_try)->long double {
//                 long double d = (long double)dl_try;
//                 return a_ld*d*d + b_ld*d + c_ld;
//             };
//
//             // ---------- cálculo do dl (com salvaguardas) ----------
//             REAL dl = 0.0; int rootIdx = 0;
//             if (it == 0) {
//                 // “curvature safeguard” no preditor
//                 static REAL sgn = 1.0;
//                 if (Dot(dwb, dw) < 0.0) sgn = -sgn;
//                 dl = sgn * std::fabs(compute_dlambda0_riks(dwb, dw, L));
//             } else {
//                 dl = compute_dlambda_riks(dwb, dws, dw, L, rootIdx);
//             }
//
//             // clamp em dl
//             const REAL dl_max = 0.25;        // ajuste conforme seu problema
//             if (dl >  dl_max) dl =  dl_max;
//             if (dl < -dl_max) dl = -dl_max;
//
//             // // se g piorou demais, tente a outra raiz
//             // long double g_old = g_eval(dl);
//             // if (std::fabsl(g_old) > 0.5L*(long double)L*(long double)L) {
//             //     int dummy=0;
//             //     REAL dl_alt = compute_dlambda_riks(dwb, dws, dw, L, dummy); // recomputa (ou guarde o par)
//             //     long double g_alt = g_eval(dl_alt);
//             //     if (std::fabsl(g_alt) < std::fabsl(g_old)) dl = dl_alt;
//             // }
//
//             // ---------- passo candidato ----------
//             TPZFMatrix<STATE> dwtot = (dwb * dl);  dwtot += dws;
//
//             // ---------- backtracking no par (Δu, δλ) ----------
//             auto try_step = [&](REAL alpha, TPZFMatrix<STATE>& u_try, REAL &lambda_try){
//                 u_try = u; u_try += dwtot * alpha;
//                 lambda_try = lambda + alpha*dl;
//                 an.LoadSolution(u_try);
//                 bcmat->Val2()[loaddir] = 0.0;
//                 an.Assemble();
//                 TPZFMatrix<STATE> rhs_try = an.Rhs();           // -FINT(u_try)
//                 TPZFMatrix<STATE> R_try   = FEXT*lambda_try;
//                 R_try += rhs_try;
//                 return Norm(R_try);
//             };
//
//             REAL alpha = 1.0;
//             TPZFMatrix<STATE> u_try; REAL lambda_try = lambda;
//             REAL normR_curr = Norm(R); // resíduo antes do passo
//             REAL normR_test = try_step(alpha, u_try, lambda_try);
//
//             int bt = 0;
//             while (normR_test > 0.9*normR_curr && bt < 6) { // critério simples de Armijo
//                 alpha *= 0.5;
//                 normR_test = try_step(alpha, u_try, lambda_try);
//                 bt++;
//             }
//
//             // aceita o passo (talvez encurtado)
//             u = u_try;
//             lambda = lambda_try;
//             TPZFMatrix<STATE> dwtot_eff = dwtot * alpha;
//             du += dwtot_eff;
//             dw += dwtot_eff;
//
//             // ---------- reavaliar resíduo ATUALIZADO (ESSENCIAL!) ----------
//             bcmat->Val2()[loaddir] = 0.0;
//             an.Assemble();
//             rhs   = an.Rhs();                 // -FINT(u atualizado)
//             R     = FEXT * lambda;   R += rhs;
//             normR = Norm(R);
//
//             // ---------- logs úteis ----------
//             long double gp_ld = b_ld + 2.0L*a_ld*(long double)dl;
//             long double g_rel = std::fabsl(g_eval(dl)) / std::max<long double>((long double)L*(long double)L, 1e-30L);
//
//             std::cout << "  [it " << it
//             << "] ||dws||="   << normdws
//             << " ||dwb||="   << normdwb
//             << " dl="         << dl
//             << " alpha="      << alpha
//             << " lambda="     << lambda
//             << " ||dwtot||="  << Norm(dwtot_eff)
//             << "\n     ||R||=" << normR
//             << "  ||FINT||="  << Norm(rhs)
//             << "  lambda||FEXT||=" << lambda * normFEXT
//             << "\n     g="    << (STATE)g_eval(dl)
//             << "  g'="        << (STATE)gp_ld
//             << "  |g|/L^2="   << (STATE)g_rel
//             << std::endl;
//
//
//             if(normR>1.e5)break;
//             it++;
//
//         }
//
//
//         okconv = (normR <= etol_inner);
//
//         if(okconv)
//         {
//             REAL uz = UyAtNode3D(cmesh, 0.5,0.5,1.);
//             std::cout << uz << " " << lambda * normFEXT << std::endl;
//             std::cout <<" (convergiu)\n" <<std::endl;
//             an.AcceptSolution();
//             PostElastoplastic(cmesh,vtkfile,matid,step,dim);
//             lambdan=lambda;
//             u_acc = u;
//             const int ndesi = 5;
//             L *= REAL(ndesi) / std::max(1, it);
//             if (L > 0.01) L = 0.01;
//         }else{
//             std::cout <<" (não convergiu)\n" <<std::endl;
//             lambda=lambdan*0.9;
//             u.Zero();
//             an.LoadSolution(u);
//             L*=0.5;
//         }
//
//
//
//         step++;
//     }
//
//     return okconv ? lambda : 0.;
// }


// REAL IterativeProcessArcLength(TPZCompMesh* cmesh, int loaddir, int indexbc)
// {
//     TPZElastoPlasticAnalysis an(
//         cmesh, std::cout,
//         TPZElastoPlasticAnalysis::ELineSearch::Armijo
//     );
//
//     auto* bcmat = dynamic_cast<TPZBndCondT<STATE>*>(cmesh->FindMaterial(indexbc));
//     if (!bcmat) {
//         std::cout << "[ArcLength] BndCond material não encontrado (indexbc=" << indexbc << ")\n";
//         DebugStop();
//     }
//
//     int dim =cmesh->Dimension();
//     // Parâmetros
//     const int  maxit_inner = 10;
//     const REAL etol_inner  = 1.e-3;
//     REAL lambda            = 0.01;       // λ inicial
//     REAL lambdan=lambda;
//     REAL L                 = 1.e-4;      // alvo de arco (norma do incremento total em u)
//
//     // Estado aceito acumulado
//     TPZFMatrix<STATE> u_acc = an.Solution(); u_acc.Zero();
//     // Incremento acumulado no passo corrente
//     TPZFMatrix<STATE> dw = an.Solution();   dw.Zero();
//
//     // ---- Monta FEXT com uma BC "unitária" (ou escala desejada) ----
//     bcmat->Val2()[loaddir] = 100;   // valor de referência para extrair o vetor de cargas
//     an.Assemble();
//     TPZFMatrix<STATE> FEXT = an.Rhs(); // isto é o "rhs" gerado pela BC -> aqui use como vetor de cargas
//     REAL normFEXT = Norm(FEXT);
//     std::cout << " Norm(FEXT) = " << normFEXT << std::endl;
//     bcmat->Val2()[loaddir] = 0.0;     // zera a BC mecânica para o loop de Newton
//     int matid=1;
//     std::string vtkfile="cubeal.vtk";
//     // ---- Um único passo de arco (ajuste se quiser mais) ----
//     int  step   = 0;
//     bool okconv = false;
//     while (step < 20)
//     {
//         std::cout << "\n[Arc step " << (step+1) << "]  L=" << L << "  lambda=" << lambda << std::endl;
//
//         // Estado de trabalho (começa no aceito)
//         TPZFMatrix<STATE> u  = u_acc;
//         TPZFMatrix<STATE> du = an.Solution(); du.Zero(); // incremento total do passo
//
//         int  it    = 0;
//         REAL normR = 1e30;
//
//         dw.Zero();
//
//         while (it < maxit_inner && normR > etol_inner)
//         {
//             bcmat->Val2()[loaddir] = 0.0;
//             an.Assemble();
//             TPZFMatrix<STATE> rhs = an.Rhs();   // rhs = -FINT
//
//             // --- resíduo físico: R = λ FEXT + rhs  (pois rhs = -FINT) ---
//             TPZFMatrix<STATE> R = FEXT * lambda;
//             R += rhs;
//
//             // --- PARTICULAR: K dws = R ---
//             an.Rhs() = R;
//             an.Solve();
//             TPZFMatrix<STATE> dws = an.Solution();
//
//             // --- DIREÇÃO EM λ: K dwb = FEXT ---
//             an.Rhs() = FEXT;
//             an.Solve();
//             TPZFMatrix<STATE> dwb = an.Solution();
//
//             const REAL normdws = Norm(dws);
//             const REAL normdwb = Norm(dwb);
//
//             // ---------- coeficientes da restrição g(dl)=a dl^2 + b dl + c ----------
//             TPZFMatrix<STATE> t = dw;  t += dws;     // t = dw + dws
//             long double a_ld = (long double)Dot(dwb,dwb);
//             long double b_ld = 2.0L * (long double)Dot(dwb,t);
//             long double c_ld = (long double)Dot(t,t) - (long double)L*(long double)L;
//             auto g_eval = [&](REAL dl_try)->long double {
//                 long double d = (long double)dl_try;
//                 return a_ld*d*d + b_ld*d + c_ld;
//             };
//
//             // ---------- cálculo do dl (com salvaguardas) ----------
//             REAL dl = 0.0; int rootIdx = 0;
//             if (it == 0) {
//                 // “curvature safeguard” no preditor
//                 static REAL sgn = 1.0;
//                 if (Dot(dwb, dw) < 0.0) sgn = -sgn;
//                 dl = sgn * std::fabs(compute_dlambda0_riks(dwb, dw, L));
//             } else {
//                 dl = compute_dlambda_riks(dwb, dws, dw, L, rootIdx);
//             }
//
//             // clamp em dl
//             const REAL dl_max = 0.25;        // ajuste conforme seu problema
//             if (dl >  dl_max) dl =  dl_max;
//             if (dl < -dl_max) dl = -dl_max;
//
//             // se g piorou demais, tente a outra raiz
//             long double g_old = g_eval(dl);
//             if (std::fabsl(g_old) > 0.5L*(long double)L*(long double)L) {
//                 int dummy=0;
//                 REAL dl_alt = compute_dlambda_riks(dwb, dws, dw, L, dummy); // recomputa (ou guarde o par)
//                 long double g_alt = g_eval(dl_alt);
//                 if (std::fabsl(g_alt) < std::fabsl(g_old)) dl = dl_alt;
//             }
//
//             // ---------- passo candidato ----------
//             TPZFMatrix<STATE> dwtot = (dwb * dl);  dwtot += dws;
//
//             // ---------- backtracking no par (Δu, δλ) ----------
//             auto try_step = [&](REAL alpha, TPZFMatrix<STATE>& u_try, REAL &lambda_try){
//                 u_try = u; u_try += dwtot * alpha;
//                 lambda_try = lambda + alpha*dl;
//                 an.LoadSolution(u_try);
//                 bcmat->Val2()[loaddir] = 0.0;
//                 an.Assemble();
//                 TPZFMatrix<STATE> rhs_try = an.Rhs();           // -FINT(u_try)
//                 TPZFMatrix<STATE> R_try   = FEXT*lambda_try;
//                 R_try += rhs_try;
//                 return Norm(R_try);
//             };
//
//             REAL alpha = 1.0;
//             TPZFMatrix<STATE> u_try; REAL lambda_try = lambda;
//             REAL normR_curr = Norm(R); // resíduo antes do passo
//             REAL normR_test = try_step(alpha, u_try, lambda_try);
//
//             int bt = 0;
//             while (normR_test > 0.9*normR_curr && bt < 6) { // critério simples de Armijo
//                 alpha *= 0.5;
//                 normR_test = try_step(alpha, u_try, lambda_try);
//                 bt++;
//             }
//
//             // aceita o passo (talvez encurtado)
//             u = u_try;
//             lambda = lambda_try;
//             TPZFMatrix<STATE> dwtot_eff = dwtot * alpha;
//             du += dwtot_eff;
//             dw += dwtot_eff;
//
//             // ---------- reavaliar resíduo ATUALIZADO (ESSENCIAL!) ----------
//             bcmat->Val2()[loaddir] = 0.0;
//             an.Assemble();
//             rhs   = an.Rhs();                 // -FINT(u atualizado)
//             R     = FEXT * lambda;   R += rhs;
//             normR = Norm(R);
//
//             // ---------- logs úteis ----------
//             long double gp_ld = b_ld + 2.0L*a_ld*(long double)dl;
//             long double g_rel = std::fabsl(g_eval(dl)) / std::max<long double>((long double)L*(long double)L, 1e-30L);
//
//             std::cout << "  [it " << it
//             << "] ||dws||="   << normdws
//             << " ||dwb||="   << normdwb
//             << " dl="         << dl
//             << " alpha="      << alpha
//             << " lambda="     << lambda
//             << " ||dwtot||="  << Norm(dwtot_eff)
//             << "\n     ||R||=" << normR
//             << "  ||FINT||="  << Norm(rhs)
//             << "  lambda||FEXT||=" << lambda * normFEXT
//             << "\n     g="    << (STATE)g_eval(dl)
//             << "  g'="        << (STATE)gp_ld
//             << "  |g|/L^2="   << (STATE)g_rel
//             << std::endl;
//
//             it++;
//
//         }
//
//
//         okconv = (normR <= etol_inner);
//
//         if(okconv)
//         {
//             REAL uz = UyAtNode3D(cmesh, 0.5,0.5,1.);
//             std::cout << uz << " " << lambda * normFEXT << std::endl;
//              std::cout <<" (convergiu)\n" <<std::endl;
//             an.AcceptSolution(0);
//             PostElastoplastic(cmesh,vtkfile,matid,step,dim);
//             lambdan=lambda;
//             u_acc = u;
//         }else{
//             std::cout <<" (não convergiu)\n" <<std::endl;
//             lambda=lambdan;
//             u=u_acc;
//             an.LoadSolution(u);
//             // const REAL s = 0.5;
//             // FEXT *= s;
//             // lambda /= s;
//             // normFEXT = Norm(FEXT);
//         }
//
//         // Aceita o passo e atualiza acumulados
//
//
//
//
//         // Ajuste adaptativo simples do alvo de arco (opcional)
//         const int ndesi = 10;
//         //L *= REAL(ndesi) / std::max(1, it);
//         //if (L > 0.00000) L = 0.001;
//         // Adaptação simples de L pelo nº de iterações (trust-region)
//         auto adaptL = [&](REAL L_in, int it){
//             const REAL it_tar = 6.0; // alvo ~6 iterações
//             REAL scale = it_tar / std::max<REAL>(1., it);
//             scale = std::max<REAL>(0.5, std::min<REAL>(2.0, scale));
//             REAL L_out = L_in * scale;
//             // limites globais (ajuste ao seu problema)
//             L_out = std::max<REAL>(1e-6, std::min<REAL>(1e-2, L_out));
//             return L_out;
//         };
//         L = adaptL(L, it);
//
//
//         step++;
//     }
//
//     return okconv ? lambda : 0.;
// }

/*
static STATE compute_dlambda0_riks(const TPZFMatrix<STATE>& dwb,
                                   const TPZFMatrix<STATE>& dw,
                                   STATE L)
{
    const STATE ndwb = Norm(dwb);
    if (ndwb < (STATE)1e-30) return (STATE)0.0;
    STATE signum;
    if (Norm(dw) < (STATE)1e-30) signum = (STATE)1.0;     // 1º passo: seguir em frente
    else                         signum = (Dot(dw,dwb) >= 0.0 ? 1.0 : -1.0);
    return signum * (L / ndwb);
}

// Corretor (k>1): a δλ² + b δλ + c = 0, com DTANG = dwb e t = dw + dws
// a = ||dwb||²
// b = 2 (dwb · (dw + dws))
// c = ||dw + dws||² − L²
// escolha da raiz: maximize ( (dw + dws + δλ*dwb) · dw )
static STATE compute_dlambda_riks(const TPZFMatrix<STATE>& dwb,
                                  const TPZFMatrix<STATE>& dws,
                                  const TPZFMatrix<STATE>& dw,
                                  STATE L, int& rootIdx)
{
    const long double a = (long double)Dot(dwb,dwb);
    if (a < 1e-30L) { rootIdx = 1; return (STATE)0.0; }

    const STATE dot_t_dwb = Dot(dwb,dws) + Dot(dwb,dw);
    const long double b = 2.0L * (long double)dot_t_dwb;

    const STATE t2 = Dot(dws,dws) + Dot(dw,dw) + 2.0*Dot(dw,dws);
    const long double c = (long double)t2 - (long double)L*(long double)L;

    long double disc = b*b - 4.0L*a*c;
    if (disc < 0.0L) disc = 0.0L;
    const long double sq = sqrt(disc);

    const STATE dl1 = (STATE)((-b - sq) / (2.0L*a));
    const STATE dl2 = (STATE)((-b + sq) / (2.0L*a));

    auto score = [&](STATE dl)->STATE {
        // (dw + dws + dl*dwb) · dw  (Souza Neto 4.118)
        return Dot(dw,dw) + Dot(dws,dw) + dl*Dot(dwb,dw);
    };
    const STATE s1 = score(dl1), s2 = score(dl2);
    if (s1 > s2) { rootIdx = 1; return dl1; }
    else         { rootIdx = 2; return dl2; }
}
*/


// REAL IterativeProcessArcLength(TPZCompMesh* cmesh,int loaddir,int indexbc)
// {
//     TPZElastoPlasticAnalysis an(cmesh, std::cout,TPZElastoPlasticAnalysis::ELineSearch::Armijo);
//
//     auto* bcmat = dynamic_cast<TPZBndCondT<STATE>*>(cmesh->FindMaterial(indexbc));
//
//     int maxsteps=10;
//     int maxit_inner=10;
//     REAL etol_outer=0.01;
//     REAL etol_inner=1.e-6;
//     bool converged = false;
//     REAL lambda0=0.01;
//     REAL lambda = lambda0;
//     REAL L=0.00001;
//     REAL lambda_prev = std::numeric_limits<REAL>::infinity();
//
//     TPZFMatrix<STATE> u_acc = an.Solution();
//     u_acc.Zero();
//
//     TPZFMatrix<STATE> dw(an.Solution());
//     dw.Zero();
//
//     int step = 0;
//     REAL diff = 1e9;
//
//     const REAL load0=bcmat->Val2()[loaddir];
//
//     bcmat->Val2()[loaddir]=199;
//     TPZFMatrix<STATE> FEXT;
//     an.Assemble();
//     FEXT = an.Rhs();
//     cout << " Norm(FEXT) = " << Norm(FEXT)  <<std::endl;
//     bcmat->Val2()[loaddir]=0;
//
//     // ===== FIM DO CHECK =====
//     while (step < 1)
//     {
//         std::cout << "\n[Arc step " << (step+1) << "]  L=" << L<< "  lambda=" << lambda << std::endl;
//
//         TPZFMatrix<STATE> u = u_acc;
//
//         int it = 0;
//         REAL normdu = 1e30;
//         REAL normR = 1e30;
//         dw.Zero();
//
//         while (it < maxit_inner && normR > etol_inner)
//         {
//             TPZFMatrix<STATE> dwb,dws,FINT;
//
//             // --- monta K(u) e FINT(u) com BC zerado ---
//             bcmat->Val2()[loaddir]=0;
//             an.Assemble();
//             FINT=an.Rhs();
//
//             // --- resíduo físico: Rphys = FINT - lambda*FEXT ---
//            //TPZFMatrix<STATE> R = FINT - (FEXT*lambda);
//             TPZFMatrix<STATE> R =  (FEXT*lambda)-FINT;
//             // --- PARTICULAR: K dws = -Rphys (sinal corrigido) ---
//             R *= -1.0;
//             REAL n0 = Norm(R);
//             an.Rhs()=R;
//             an.Solve();
//             dws=an.Solution();
//
//             // --- DIREÇÃO λ: K dwb = FEXT ---
//             an.Rhs()=FEXT*(-1.0);
//             an.Solve();
//             dwb=an.Solution();
//
//             STATE normdwb,normdws;
//             normdwb=Norm(dwb);
//             normdws=Norm(dws);
//
//
//             REAL dl = 0.0;
//             if (it == 0) {
//                 dl = compute_dlambda0_riks(dwb, dw /*Δu^(k-1)*/, L);
//             } else {
//                 int rootIdx = 0;
//                 dl = compute_dlambda_riks(dwb, dws, dw /*Δu^(k-1)*/, L, rootIdx);
//             }
//
//             // // (opcional) freio simples para evitar saltos
//             // if (dl >  0.05) dl =  0.05;
//             // if (dl < -0.05) dl = -0.05;
//
//             TPZFMatrix<STATE> dwtot = dwb*dl + dws;
//             u      += dwtot;
//             dw     += dwtot;
//             lambda += dl;
//
//             normdu = Norm(dwtot);
//
//             an.LoadSolution(u);
//
//             an.LoadSolution(u);
//
//             // --- reavalia Rphys no estado ATUALIZADO (para medir normR corretamente) ---
//             bcmat->Val2()[loaddir]=0;
//             an.Assemble();
//             FINT = an.Rhs();
//             R    = (FEXT*lambda) - FINT;   // (mesma convenção que você usou no restante)
//             normR= Norm(R);
// /*
//             std::cout << "[ELASTIC CHECK] ||R0||="<< n0
//             << "  ||R1||="<< n1 << "  (esperado ~ 0 após 1 it)\n";*/
//
//             cout << "  [it " << it
//             << "] ||dws||=" << normdws
//             << " ||dwb||=" << normdwb
//             << " dl=" << dl
//             << " lambda=" << lambda
//             << " ||dwtot||=" << normdu
//             << " \n ||R||=" << normR
//             << " ||FINT||=" << Norm(FINT)
//             << " lambda||FEXT||=" << lambda*Norm(FEXT)
//             << " ||FEXT||=" << Norm(FEXT)
//             << endl;
//
//             it++;
//         }
//
//         bool converged = (normR <= etol_inner);
//         std::cout << (converged ? " (convergiu)\n" : " (não convergiu)\n");
//
//         an.AcceptSolution(0);
//         u_acc = u;
//         diff = std::fabs(lambda - lambda_prev);
//         lambda_prev = lambda;
//
//         const int ndesi = 10;
//         L *= REAL(ndesi) / std::max(1, it);
//         if (L > 0.1) L = 0.1;
//         step++;
//     }
//
//     return converged ? lambda : 0.;
// }
//
//
