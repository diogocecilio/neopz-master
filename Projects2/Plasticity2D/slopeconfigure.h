// slopeconfigure.h
#ifndef SLOPE_H // include guard
#define SLOPE_H

#include <cmath>
#include <set>
#include <iostream>
#include <fstream>
#include <string>
using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::duration;
using std::chrono::milliseconds;
using std::chrono::seconds;
    int bottombc_slope=-1;
    int rigthbc_slope = -2;
    int leftbc_slope =-3;
    int toprigthbc_slope = -4;
    int topleftbc_slope = -5;
    int rampbc_slope = -6;

typedef TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse> plasticmorh;
typedef   TPZMatElastoPlastic2D < TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse>, TPZElastoPlasticMem > plasticmat;
//typedef TPZNonLinearAnalysis anal;
typedef TPZElastoPlasticAnalysis anal;
using namespace std;

class Slope
{
public:

    TPZCompMesh *fCompMesh;
    Slope(TPZGeoMesh * gmesh,int porder);
    enum  SolverType { EStep, EPardiso};
    TPZCompMesh * CreateCMeshDarcy ( TPZGeoMesh *gmesh, int pOrder );

    void SlopeAnalysis(anal *analysis,int type,int numthreads);
    void ApplyGravityLoad(TPZManVector<REAL, 3> bodyforce);
    void PostPlasticity(string vtkd);
    void  CreatePostProcessingMesh (TPZPostProcAnalysis * PostProcess );
    void PostProcessVariables ( TPZStack<std::string> &scalNames, TPZStack<std::string> &vecNames );
    void InitializeSimulation(int porder,int ref);
    void Solve();
    TPZCompMesh * CreateCMeshElastoplastic ( TPZGeoMesh *gmesh, int pOrder );

protected:

    TPZGeoMesh *fGmesh;
    //TPZElastoPlasticAnalysis *fAnalysis;


};

void Slope::Solve()
{
    //TPZElastoPlasticAnalysis *analysis = new TPZElastoPlasticAnalysis( fCompMesh,std::cout );
    anal *analysis = new anal(fCompMesh,cout);

    TPZManVector<REAL,3> gravityload(3,0);
    gravityload[1]=-20.;
    ApplyGravityLoad(gravityload);

    int numthreads = 0;
    int type =Slope::EStep;

    SlopeAnalysis(analysis,type,numthreads);

 	REAL tol=1.e-3;
 	int numiter=10;
 	bool linesearch=true;
 	bool checkconv= false;
    bool convordiv;
    int iters;
    std::cout << "start solving"<< endl ;
    auto t1 = high_resolution_clock::now();
    //analysis->IterativeProcess(std::cout,tol,numiter,linesearch,checkconv);
    //analysis->IterativeProcessPrecomputedMatrix(std::cout,tol,numiter,linesearch);
    //analysis->IterativeProcess(std::cout,tol,numiter,linesearch,  checkconv,convordiv) ;
    int niter_update_jac=1;
    analysis->IterativeProcess(std::cout, tol,  numiter, niter_update_jac,  linesearch);
    auto t2 = high_resolution_clock::now();
    auto ms_int = duration_cast<milliseconds> ( t2 - t1 );
    std::cout << "tempo total = "<<ms_int.count() << " ms\n";

     string vtkd = "postlin.vtk";
     //PostPlasticity(vtkd);


     analysis->AcceptSolution();

     TPZLinearAnalysis * anal = new TPZLinearAnalysis(fCompMesh);
     //analysis->Solution().Print("sol");
    //TPZPostProcAnalysis * postprocdeter = new TPZPostProcAnalysis();
    //CreatePostProcessingMesh ( postprocdeter);
   // TPZNonLinearAnalysis * postprocdeter = new TPZNonLinearAnalysis(fCompMesh,cout);

    TPZVec<int> PostProcMatIds ( 1,1 );

    TPZStack<std::string> PostProcVars, scalNames, vecNames;

    PostProcessVariables ( scalNames, vecNames );

    //string vtkd = "postprocessdeter.vtk";
    anal->DefineGraphMesh ( 2,scalNames,vecNames,vtkd );

    anal->PostProcess ( 0 );

}

Slope::Slope( TPZGeoMesh * gmesh,int porder)
{
    fGmesh = gmesh;
    fCompMesh = CreateCMeshElastoplastic ( fGmesh, porder );

}

TPZCompMesh * Slope::CreateCMeshElastoplastic ( TPZGeoMesh *gmesh, int pOrder )
{
	    // Creating computational mesh:
    TPZCompMesh * cmesh = new TPZCompMesh ( gmesh );

    //cmesh->Print(std::cout);
	cmesh->SetDefaultOrder ( pOrder );

	int dim = 2 ;

	cmesh->SetDimModel ( dim );

	cmesh->SetAllCreateFunctionsContinuousWithMem();

	int matid=1;
	int planestrain=1;

	//auto * material = new plasticmat(matid,planestrain);

    REAL E=20000;
    REAL nu=0.49;
    TPZElasticResponse  elasticresponse;

    //REAL lambda = nu * E / ((1. + nu)*(1. - 2. * nu));
    //REAL mu = E / (2. * (1. + nu));
    elasticresponse.SetEngineeringData(E,nu);

    // Mohr Coulomb data
    REAL mc_cohesion    = 10.;//kPa
    REAL mc_phi         = 30.*M_PI/180;
    REAL mc_psi         = mc_phi;

    //elasticresponse.Print(std::cout);

    plasticmorh mohrcoulombplasticstep;

    mohrcoulombplasticstep.fYC.SetUp ( mc_phi, mc_psi, mc_cohesion, elasticresponse );

    mohrcoulombplasticstep.fER = elasticresponse;

    //mohrcoulombplasticstep.Print(std::cout);

    int PlaneStrain = 1;

    plasticmat * material = new plasticmat ( matid,PlaneStrain );

    material->SetPlasticityModel(mohrcoulombplasticstep);

    material->SetId(matid);


   // material->

    //material->Print(std::cout);

    cmesh->InsertMaterialObject ( material );

    //cmesh->Print(std::cout);

    // boundary condition
    TPZFMatrix<STATE>  val1 ( 2,2,0. );

    TPZManVector<STATE,2> val2 ( 2,0. );

	int directionaldirichlet = 3 ;

	val2[0]=1;
	val2[1]=1;
    auto * BCond0 = material->CreateBC ( material, bottombc_slope, directionaldirichlet, val1, val2 );

	val2[0]=1;
	val2[1]=0;
    auto * BCond1 = material->CreateBC ( material, rigthbc_slope, directionaldirichlet, val1, val2 );

	val2[0]=1;
	val2[1]=0;
	auto * BCond2 = material->CreateBC ( material, leftbc_slope, directionaldirichlet, val1, val2 );

	cmesh->InsertMaterialObject ( BCond0 );
    cmesh->InsertMaterialObject ( BCond1 );
	cmesh->InsertMaterialObject ( BCond2 );


    //Creating computational elements that manage the space of the mesh:
    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();

    return cmesh;
}


void Slope::ApplyGravityLoad(TPZManVector<REAL, 3> bodyforce)
{
    plasticmat * body= dynamic_cast<plasticmat *> ( fCompMesh->FindMaterial ( 1 ) );
    body->SetBodyForce ( bodyforce );

}

void  Slope::SlopeAnalysis(anal *analysis,int type,int numthreads)
{

    //TPZElastoPlasticAnalysis *analysis = new TPZElastoPlasticAnalysis();
    //analysis->SetCompMesh(fCompMesh,true);
switch(type) {
    case SolverType::EStep:
        {
			cout << "Solver called with TPZStepSolver\n";
            TPZSkylineStructMatrix<STATE> matskl(fCompMesh);
            matskl.SetNumThreads ( numthreads );
             analysis->SetStructuralMatrix ( matskl );
             TPZStepSolver<STATE> step;
             step.SetDirect ( ELDLt );
             analysis->SetSolver ( step );
             break;
        }
    case SolverType::EPardiso:
        {
            cout << "Solver called with TPZPardisoSolver\n";
            TPZSSpStructMatrix<STATE> SSpStructMatrix(fCompMesh);
            SSpStructMatrix.SetNumThreads ( numthreads );
            analysis->SetStructuralMatrix(SSpStructMatrix);
            TPZPardisoSolver<REAL> *pardiso = new TPZPardisoSolver<REAL>;
            analysis->SetSolver(*pardiso);
            break;
        }
    default:
        {
            cout << "Solver was not initialized properly\n";
            DebugStop();
        }
    }
}



void Slope::PostPlasticity(string vtkd)
{

   TPZPostProcAnalysis * postprocdeter = new TPZPostProcAnalysis();
    CreatePostProcessingMesh ( postprocdeter);
   // TPZNonLinearAnalysis * postprocdeter = new TPZNonLinearAnalysis(fCompMesh,cout);

    TPZVec<int> PostProcMatIds ( 1,1 );

    TPZStack<std::string> PostProcVars, scalNames, vecNames;

    PostProcessVariables ( scalNames, vecNames );

    //string vtkd = "postprocessdeter.vtk";
    postprocdeter->DefineGraphMesh ( 2,scalNames,vecNames,vtkd );

    postprocdeter->PostProcess ( 0 );

    delete postprocdeter;
}

void  Slope::CreatePostProcessingMesh (TPZPostProcAnalysis * PostProcess )
{
    if ( PostProcess->ReferenceCompMesh() != fCompMesh )
    {

        PostProcess->SetCompMesh ( fCompMesh );

        TPZVec<int> PostProcMatIds ( 1,1 );
        TPZStack<std::string> PostProcVars, scalNames, vecNames;
        PostProcessVariables ( scalNames, vecNames );

        for ( int i=0; i<scalNames.size(); i++ )
        {
            PostProcVars.Push ( scalNames[i] );
        }
        for ( int i=0; i<vecNames.size(); i++ )
        {
            PostProcVars.Push ( vecNames[i] );
        }
        //
        PostProcess->SetPostProcessVariables ( PostProcMatIds, PostProcVars );
        TPZFStructMatrix<REAL> structmatrix ( PostProcess->Mesh() );
        structmatrix.SetNumThreads ( 0 );
        PostProcess->SetStructuralMatrix ( structmatrix );
    }
    //
    //Chamar com o analysis e nao com o postanalysis pois tem o acumulo de sols
    PostProcess->TransferSolution();

}

void Slope::PostProcessVariables ( TPZStack<std::string> &scalNames, TPZStack<std::string> &vecNames )
{

     scalNames.Push ( "StressJ2" );
     vecNames.Push ( "Displacement" );


}

#endif /* SLOPE_H */
