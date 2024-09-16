#include <pzgmesh.h> //for TPZGeoMesh
#include <pzcmesh.h> //for TPZCompMesh
#include <TPZGeoMeshTools.h> //for TPZGeoMeshTools::CreateGeoMeshOnGrid
#include <MMeshType.h> //for MMeshType
#include <pzmanvector.h>//for TPZManVector
#include <TPZBndCond.h> //for TPZBndCond
#include <TPZLinearAnalysis.h> //for TPZLinearAnalysis
#include <TPZSSpStructMatrix.h> //symmetric sparse matrix storage
#include <pzskylstrmatrix.h> //symmetric skyline matrix storage
#include <pzstepsolver.h> //for TPZStepSolver
#include <TPZSimpleTimer.h>
#include "TPZPardisoSolver.h"
#include "tpzgeoelrefpattern.h"
#include "Elasticity/TPZElasticity2D.h"
#include "TPZVTKGeoMesh.h"
#include "tpzsparseblockdiagonalstructmatrix.h"
using namespace std;

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.adaptivity"));
static LoggerPtr loggerconv(Logger::getLogger("pz.adaptivity.conv"));
static LoggerPtr loggerpoint(Logger::getLogger("pz.adaptivity.points"));
#endif

#include <time.h>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <chrono>
#include "pzblockdiag.h"
#include "tpzsparseblockdiagonal.h"
using namespace std;
using namespace std;
using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::duration;
using std::chrono::milliseconds;
using std::chrono::seconds;

TPZGeoMesh *CreateGeoMesh();
TPZCompMesh *CreateMesh(TPZGeoMesh *gmesh);
void SolveSelfWeigthBar(string filename);
void SolveBendingClampedBeam(string filename);

TPZGeoMesh *CreateGeoMeshBending();
TPZCompMesh *CreateMeshBending(TPZGeoMesh *gmesh);
TPZVec<REAL> findnodalsol(TPZCompMesh *cmesh,TPZVec<REAL> coord) ;

// bi-dimensional problem for elasticity
int main() {
    string filename =  "selfweigth.vtk";
    SolveSelfWeigthBar(filename);

//     filename =  "clamped.vtk";
//     SolveBendingClampedBeam(filename);

}


void SolveSelfWeigthBar(string filename)
{
	// Creating geometric mesh
	TPZGeoMesh *gmesh = CreateGeoMesh();

	// Creating computational mesh (approximation space and materials)
	int p =6;
    TPZCompEl::SetgOrder(p);
    TPZCompMesh *cmesh = CreateMesh(gmesh);
	// Solving linear equations

	// Initial steps
	TPZLinearAnalysis an (cmesh,true);
    TPZLinearAnalysis an2 (cmesh,true);


    TPZSSpStructMatrix<STATE> SSpStructMatrix(cmesh);
	an.SetStructuralMatrix(SSpStructMatrix);

    TPZSkylineStructMatrix<STATE> SkylineStructMatrix(cmesh);
    an2.SetStructuralMatrix(SkylineStructMatrix);


	TPZStepSolver<REAL> *direct = new TPZStepSolver<REAL>;
    direct->SetDirect(ECholesky);

    TPZPardisoSolver<REAL> *pardiso = new TPZPardisoSolver<REAL>;

	an.SetSolver(*pardiso);
    an2.SetSolver(*direct);

    an.Assemble();
    an2.Assemble();

    int sz=an.Rhs().Rows();
    cout << "sz = " << sz <<endl;

    //an.Solve();



     TPZPardisoSolver<REAL>  * step1=dynamic_cast<TPZPardisoSolver<REAL> *> (an.Solver());

     TPZStepSolver<REAL>  * step2=dynamic_cast<TPZStepSolver<REAL> *> (an2.Solver());
//
     TPZAutoPointer<TPZMatrix<REAL> > stiffnessmatrix1 = step1->Matrix();

     TPZAutoPointer<TPZMatrix<REAL> > stiffnessmatrix2 = step2->Matrix();
//
     std::ofstream outmat1("stiffnessmatrixpardiso.dat");
     std::ofstream outmat2("stiffnessmatrixstepsolver.dat");
//
   //  stiffnessmatrix1->Print(outmat1);
    // stiffnessmatrix2->Print(outmat2);
//
//     an.Rhs().Print("load");

     std::cout << "start solving with TPZPardisoSolver"<< endl ;
     auto t1 = high_resolution_clock::now();

     an.Solve();
//
     auto t2 = high_resolution_clock::now();
     auto ms_int = duration_cast<seconds> ( t2 - t1 );
     std::cout << "tempo total Pardiso = "<<ms_int.count() << " s\n";

     std::cout << "start solving with TPZStepSolver"<< endl ;
     t1 = high_resolution_clock::now();

     an2.Solve();
//
     t2 = high_resolution_clock::now();
     ms_int = duration_cast<seconds> ( t2 - t1 );
     std::cout << "tempo total  StepSolver = "<<ms_int.count() << " s\n";

//
//
//     //TPZMatWithMem<TPZElastoPlasticMem> *pMatWithMem2 = dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem> *> ( fcmesh->MaterialVec() [1] );
    TPZFMatrix<REAL> sol = an.Solution();
    TPZVec<REAL> coord(2);
    //cout << "Displacement solution in coord x = "<< coord[0] << " y ="  << coord[1] <<endl;
    cout << "ux = "<< sol(0,0)<< " uy ="  << sol(1,0) <<endl;
    //sol.Print(cout);

  ///Calculating approximation error
  TPZManVector<REAL,3> error;
  std::ofstream anPostProcessFile("postprocess.txt");
  an.PostProcess(error,anPostProcessFile);
  ///vtk export
  TPZVec<std::string> scalarVars(1), vectorVars(1);
  vectorVars[0] = "displacement";
  scalarVars[0] = "SigmaY";
  an.DefineGraphMesh(2,scalarVars,vectorVars,filename);
  //constexpr int resolution{1};
  an.PostProcess(0);

}
void SolveBendingClampedBeam(string filename)
{
	// Creating geometric mesh
	TPZGeoMesh *gmesh = CreateGeoMeshBending();

	// Creating computational mesh (approximation space and materials)
	int p = 2;
    TPZCompEl::SetgOrder(p);
    TPZCompMesh *cmesh = CreateMeshBending(gmesh);
	// Solving linear equations
	// Initial steps
	TPZLinearAnalysis an (cmesh,false);
	TPZSkylineStructMatrix<STATE> strskyl(cmesh);
	an.SetStructuralMatrix(strskyl);
	// Solver (is your choose)
	TPZStepSolver<REAL> *direct = new TPZStepSolver<REAL>;
	direct->SetDirect(ECholesky);
	an.SetSolver(*direct);
	delete direct;
	direct = 0;

	an.Run();

    //cout << "\n SOLUTION "<< endl;
    TPZFMatrix<REAL> sol = an.Solution();
    TPZVec<REAL> coord(2);
    //coord[0]=0.;
    //coord[1]=0.;
    //TPZVec<REAL> sol = findnodalsol(cmesh,coord);

    cout << "Displacement solution in coord x = "<< coord[0] << " y ="  << coord[1] <<endl;
    cout << "ux = "<< sol(0,0)<< " uy ="  << sol(1,0) <<endl;

	// Post processing
	TPZManVector<std::string> scalarnames(3), vecnames(3);
	scalarnames[0] = "SigmaX";
	scalarnames[1] = "SigmaY";
	scalarnames[2] = "TauXY";


   //     if(!strcmp("NormalStress",name.c_str()))        return 23;
  //  if(!strcmp("ShearStress",name.c_str()))        return 24;
  //  if(!strcmp("NormalStrain",name.c_str()))        return 25;
  //  if(!strcmp("ShearStrain",name.c_str()))        return 26;
	vecnames[0] = "displacement";
    vecnames[1] = "Strain";
    vecnames[2] = "ShearStrain";
	//vecnames[1] = "";
	an.DefineGraphMesh(2,scalarnames,vecnames,filename);

	an.PostProcess(0);
}
TPZGeoMesh *CreateGeoMeshBending()
{
    REAL co[4][2] = {{0.,0.},{100,0},{100,20},{0,20}};
    long indices[1][4] = {{0,1,2,3}};
    TPZGeoEl *elvec[1];
    TPZGeoMesh *gmesh = new TPZGeoMesh();
    gmesh->SetDimension(2);
    long nnode = 4;
    long nod;
    for(nod=0; nod<nnode; nod++) {
        long nodind = gmesh->NodeVec().AllocateNewElement();
        TPZVec<REAL> coord(2);
        coord[0] = co[nod][0];
        coord[1] = co[nod][1];
        gmesh->NodeVec()[nodind] = TPZGeoNode(nod,coord,*gmesh);
    }

    long el;
    long nelem = 1;
    for(el=0; el<nelem; el++) {
        TPZVec<long> nodind(4);
        for(nod=0; nod<4; nod++) nodind[nod]=indices[el][nod];
        //    elvec[el] = new TPZGeoElQ2d(el,nodind,1);
        long index;
        elvec[el] = gmesh->CreateGeoElement(EQuadrilateral,nodind,1,index);
    }

    TPZVec <long> TopoLine ( 2 );

    //long index=2;
    TopoLine[0] = 1;
    TopoLine[1] = 2;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( 2, TopoLine, - 1, *gmesh );//clamped in right side

   TPZVec <long> node(1);
   node[0]=0;
   new TPZGeoElRefPattern< pzgeom::TPZGeoPoint> ( 3, node, - 2, *gmesh );//load node

    node[0]=1;
   new TPZGeoElRefPattern< pzgeom::TPZGeoPoint> ( 4, node, - 3, *gmesh );//bottomrigth node

    gmesh->BuildConnectivity();


    cout << "c" << endl;
    for ( int d = 0; d<5; d++ )
    {
        int nel = gmesh->NElements();
        TPZManVector<TPZGeoEl *> subels;
        for ( int iel = 0; iel<nel; iel++ )
        {
            TPZGeoEl *gel = gmesh->ElementVec() [iel];
            gel->Divide ( subels );
        }
    }

    std::ofstream files ( "teste-mesh.vtk" );
    TPZVTKGeoMesh::PrintGMeshVTK ( gmesh,files,false );
    cout << "d" << endl;
    return gmesh;

}

TPZCompMesh *CreateMeshBending(TPZGeoMesh *gmesh)
{
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(TPZCompEl::GetgOrder());

    //TPZElasticityMaterial(int id, REAL E, REAL nu, REAL fx, REAL fy, int planestress = 1);

    auto * mat = new TPZElasticity2D(1,21000000.,0.3,0.,0.);//selfweigth

    cmesh->SetDimModel(2);

   // TPZFMatrix<REAL> val1(2,2,0.),val2(2,1,0.);
    TPZFMatrix<STATE> val1(2,2,0.);
    TPZVec<STATE> val2(2,0.);
	//TPZMaterial *bcload,*bcclamp,*bcnode;

    val2[0]=1.;
    val2[1]=1.;
    TPZBndCond *bcclamp = mat->CreateBC(mat,-1,3,val1,val2);//clamped line restrictions

    val2[0]=0.;
    val2[1]=1.;
    TPZBndCond *bcnode = mat->CreateBC(mat,-3,3,val1,val2);//bottomrigth node restrictions

    val2[0]=0.;
    val2[1]=-1000.;
    TPZBndCond *bcload = mat->CreateBC(mat,-2,1,val1,val2);//-100 N in y direction node 4


    cmesh->InsertMaterialObject(mat);
	cmesh->InsertMaterialObject(bcclamp);
    cmesh->InsertMaterialObject(bcnode);
    cmesh->InsertMaterialObject(bcload);


	cmesh->SetAllCreateFunctionsContinuous();

    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();

    return cmesh;
}

TPZGeoMesh *CreateGeoMesh() {

    REAL co[4][2] = {{0.,0.},{1,0},{1,10},{0,10}};
    long indices[1][4] = {{0,1,2,3}};
    TPZGeoEl *elvec[1];
    TPZGeoMesh *gmesh = new TPZGeoMesh();
    long nnode = 4;
    long nod;
    for(nod=0; nod<nnode; nod++) {
        long nodind = gmesh->NodeVec().AllocateNewElement();
        TPZVec<REAL> coord(2);
        coord[0] = co[nod][0];
        coord[1] = co[nod][1];
        gmesh->NodeVec()[nodind] = TPZGeoNode(nod,coord,*gmesh);
    }

    long el;
    long nelem = 1;
    for(el=0; el<nelem; el++) {
        TPZVec<long> nodind(4);
        for(nod=0; nod<4; nod++) nodind[nod]=indices[el][nod];
        //    elvec[el] = new TPZGeoElQ2d(el,nodind,1);
        long index;
        elvec[el] = gmesh->CreateGeoElement(EQuadrilateral,nodind,1,index);
    }

    TPZVec <long> TopoLine ( 2 );

    TopoLine[0] = 0;
    TopoLine[1] = 1;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( 7, TopoLine, - 1, *gmesh );//bottom


    TopoLine[0] = 2;
    TopoLine[1] = 3;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( 8, TopoLine, - 2, *gmesh );//top

   TPZVec <long> node(1);
   node[0]=3;
   new TPZGeoElRefPattern< pzgeom::TPZGeoPoint> ( 9, node, - 3, *gmesh );//top
//     long index;
//     gmesh->CreateGeoElement(EPoint,3,-3,index);
//
//     TPZGeoElBC gbc3(elvec[0],6,-3);

    gmesh->BuildConnectivity();


    for(int d=0;d<5;d++) {
    int nel = gmesh->NElements();
    for (int iel=0; iel<nel; iel++) {
        TPZManVector<TPZGeoEl *> subels;
        TPZGeoEl *gel = gmesh->ElementVec()[iel];
        gel->Divide(subels);
    }
	}

    std::ofstream files ( "self-weigth-bar.vtk" );
    TPZVTKGeoMesh::PrintGMeshVTK ( gmesh,files,false );

    return gmesh;

	return gmesh;
}


TPZCompMesh *CreateMesh(TPZGeoMesh *gmesh) {

    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(TPZCompEl::GetgOrder());

    int id=1;
    REAL E=10;
    REAL nu=0.;
    REAL fx=0;
    REAL fy=-1;
    int planestress = 1;

    // Creating elasticity material
    auto * mat = new TPZElasticity2D(id,E,nu,fx,fy,planestress);//selfweigth

	// Creating four boundary condition
    mat->Print(std::cout);
    TPZFMatrix<STATE> val1(2,2,0.);
    TPZVec<STATE> val2(2,0.);

    val2[1]=1.;
    auto*bctop = mat->CreateBC(mat,-2,3,val1,val2);

    val2[0]=1.;
    auto*bcpoint = mat->CreateBC(mat,-3,3,val1,val2);

    //val2[0]=0.;
    //val2[1]=-1.;
    //bcload = mat->CreateBC(mat,-1,1,val1,val2);

    cmesh->InsertMaterialObject(mat);
	// Inserting boundary conditions into computational mesh
	cmesh->InsertMaterialObject(bctop);

    cmesh->InsertMaterialObject(bcpoint);

    //cmesh->InsertMaterialObject(bcload);

	cmesh->SetAllCreateFunctionsContinuous();

    cmesh->AutoBuild();
   // cmesh->AdjustBoundaryElements();
  //  cmesh->CleanUpUnconnectedNodes();

    return cmesh;
}


