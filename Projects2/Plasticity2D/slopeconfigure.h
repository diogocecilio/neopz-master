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


typedef TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse> plasticmorh;
typedef   TPZMatElastoPlastic2D < TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse>, TPZElastoPlasticMem > plasticmat;
//typedef TPZNonLinearAnalysis typedefanal;
typedef TPZElastoPlasticAnalysis typedefanal;
using namespace std;

class Slope
{
public:

    enum  SolverType { EStep, EPardiso};

    Slope( TPZGeoMesh * gmesh,int porder,REAL gammaagua, REAL gammasolo);

    TPZCompMesh * CreateCMeshElastoplastic ( TPZGeoMesh *gmesh, int pOrder );
    TPZCompMesh *CreateCompMesh(TPZGeoMesh *gmesh,int porder) ;
    typedefanal *  SlopeAnalysis(int type,int numthreads);

    void ApplyGravityLoad(TPZManVector<REAL, 3> bodyforce);

    void PostPlasticity(string vtkd);

    void  CreatePostProcessingMesh (TPZPostProcAnalysis * PostProcess );

    void PostProcessVariables ( TPZStack<std::string> &scalNames, TPZStack<std::string> &vecNames );

    void InitializeSimulation(int porder,int ref);

    void Solve();

    void LoadingRamp (  REAL factor );

    REAL ShearRed ( );


    void TransferSolutionFromShearRed ( TPZManVector<TPZCompMesh*,3> vecmesh,int isol,REAL FS );
public:
    TPZCompMesh *fCompMesh;
protected:

    REAL fgammaagua;
    REAL fgammasolo;
    TPZGeoMesh *fGmesh;
   // typedefanal *fAnalysis;


};

Slope::Slope( TPZGeoMesh * gmesh,int porder,REAL gammaagua, REAL gammasolo)
{
    fgammaagua=gammaagua;
    fgammasolo=gammasolo;
    fGmesh = gmesh;
    fCompMesh = CreateCMeshElastoplastic ( fGmesh, porder );
    //fCompMesh = CreateCompMesh ( fGmesh, porder );


}

void Slope::Solve()
{
    ShearRed();
    PostPlasticity("saidax.vtk");

}
void Slope::LoadingRamp (  REAL factor )
{
    plasticmat * body= dynamic_cast<plasticmat *> ( fCompMesh->FindMaterial ( 1 ) );
    //plasticmatcrisfield * body= dynamic_cast<plasticmatcrisfield *> ( cmesh->FindMaterial ( 1 ) );
    TPZManVector<REAL, 3> force ( 3,0. );

    force[1]= ( fgammaagua-fgammasolo );

    //force[1]=(-gammasolo);
    body->SetLoadFactor ( factor );
    body->SetBodyForce ( force );

}


REAL Slope::ShearRed ( )
{

    //fgmesh->ResetReference();
   // fcmesh->LoadReferences();
    LoadingRamp(1.);

    REAL FS=0.5,FSmax=5.,FSmin=0.,tol=0.01;
    int neq = fCompMesh->NEquations();

    TPZFMatrix<REAL> displace(neq,1),displace0(neq,1);

    int counterout = 0;

    plasticmat *material = dynamic_cast<plasticmat *> ( fCompMesh->MaterialVec() [1] );
    TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse> LEMC = material->GetPlasticity();
    TPZElasticResponse ER = LEMC.fER;

    REAL phi0 = LEMC.fYC.Phi();
    REAL cohesion0 = LEMC.fYC.Cohesion();
    REAL phi,psi,c;
    int numthreads = 10;
    int type =EStep;
    //int type =EPardiso;
    bool conv=false;

    do {

        fCompMesh->Solution().Zero();
        std::cout << "FS "<< FS <<  "| step = " << counterout  <<std::endl;
        typedefanal* anal =  SlopeAnalysis(type,numthreads);

        REAL norm = 1000.;
        REAL tol2 = 1.e-3;
        int NumIter = 30;
        bool linesearch = true;
        bool checkconv = false;
        int iters;



        chrono::steady_clock sc;
        auto start = sc.now();
//IterativeProcess2(std::ostream &out,REAL tol,int numiter, bool linesearch = false, bool checkconv = false);
        conv  =anal->IterativeProcess2(cout,tol2, NumIter,  linesearch,  checkconv,iters);
        //conv  = anal->IterativeProcess ( cout, tol2, NumIter,linesearch,checkconv,iters);
        cout << "iters = " << iters <<endl;
         //conv = anal->FindRoot ( );
        norm = Norm ( anal->Rhs() );
        std::cout << "conv "<< conv <<std::endl;

        //anal->AcceptSolution();


        auto end = sc.now();
        auto time_span = static_cast<chrono::duration<double>> ( end - start );

        std::cout << "safety factor "<< FS <<  "| step = " << counterout <<" | Rhs norm = " << norm  << " | IterativeProcess Time: " << time_span.count() << " seconds !!! " <<std::endl;

        if ( conv==false ) {

            FSmax = FS;
            FS = ( FSmin + FSmax ) / 2.;
        } else {

            FSmin = FS;
			FS = 1. / ( ( 1. / FSmin + 1. / FSmax ) / 2. );
        }

        c=cohesion0/FS;
        std::cout << "coes "<<  c <<std::endl;
        phi=atan ( tan ( phi0 ) /FS );
        psi=phi;
        //void SetUp(REAL Phi, REAL Psi, REAL c, TPZElasticResponse &ER)
        LEMC.fYC.SetUp ( phi, psi, c, ER );
        material->SetPlasticityModel ( LEMC );
        counterout++;
        if(( FSmax - FSmin ) / FS < tol)
        {
            anal->AcceptSolution();
        }
    }  while ( ( FSmax - FSmin ) / FS > tol || conv==false);


        std::cout << "final safety factor "<< FS <<std::endl;
        return ( FSmax + FSmin )/2;
}





TPZCompMesh *Slope::CreateCompMesh(TPZGeoMesh *gmesh,int porder) {


    unsigned int dim  = porder;
    const std::string name ( "ElastoPlastic COMP MESH Slope Problem " );

    TPZCompMesh * cmesh =  new TPZCompMesh ( gmesh );
    cmesh->SetName ( name );
    cmesh->SetDefaultOrder (porder );
    cmesh->SetDimModel ( dim );


    //Pu = (2+pi)c = 218.183

      // Mohr Coulomb data
    //REAL mc_cohesion    = 490.;//kPa
    REAL mc_cohesion    = 10.;//N/cm^2
    REAL mc_phi         = 30.*M_PI/180;
    REAL mc_psi         = mc_phi;

    /// ElastoPlastic Material using Mohr Coulomb
    // Elastic predictor
    TPZElasticResponse ER;
    REAL nu = 0.49;
    //REAL E = 10000000.;////kPa
    REAL E = 20000.;////N/cm^2

    plasticmorh plasticstep;
    //TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse> LEMC;
    //TPZPlasticStepVoigt<TPZMohrCoulombVoigt, TPZElasticResponse> LEMC;
    ER.SetEngineeringData(E,nu);
    //ER.SetUp ( E, nu );
    plasticstep.fER =ER;
    // LEMC.SetElasticResponse( ER );
    plasticstep.fYC.SetUp ( mc_phi, mc_psi, mc_cohesion, ER );

    int PlaneStrain = 1;

    plasticmat * material = new plasticmat ( 1,PlaneStrain );

    REAL factor;
    TPZManVector<REAL, 3> bodyforce ( 3,0. );
    factor=1.;
    bodyforce[1]=-20.;

   // material->SetPlasticity ( plasticstep );

    material->SetId ( 1 );

    material->SetWhichLoadVector ( 0 ); //option to compute the total internal force vecor fi=(Bt sigma+ N (b+gradu))

    material->SetLoadFactor ( factor );

    material->SetBodyForce ( bodyforce );

    cmesh->InsertMaterialObject(material);

    TPZFMatrix<STATE>  val1 ( 2,2,0. );

    TPZManVector<STATE,2> val2 ( 2,0. );

    val2[1]=1.;
    val2[0]=1.;
    auto * bcleft = material->CreateBC(material,-1,3,val1,val2);//bottom

    val2[0]=1.;
    val2[1]=0.;
    auto * bcbottom = material->CreateBC(material,-2,3,val1,val2);//rigth

    val2[0]=1.;
    val2[1]=0;
    auto *  bcrigth = material->CreateBC(material,-5,3,val1,val2);//left


	cmesh->InsertMaterialObject(bcleft);

    cmesh->InsertMaterialObject(bcbottom);

    cmesh->InsertMaterialObject(bcrigth);



	cmesh->SetAllCreateFunctionsContinuousWithMem();

    cmesh->AutoBuild();

    std::ofstream print("cmeshslopecphi.txt");
	cmesh->Print(print);

    return cmesh;
}

TPZCompMesh * Slope::CreateCMeshElastoplastic ( TPZGeoMesh *gmesh, int pOrder )
{
	    // Creating computational mesh:
    TPZCompMesh * cmesh = new TPZCompMesh ( gmesh );

    //cmesh->Print(std::cout);
	cmesh->SetDefaultOrder ( pOrder );

	int dim = 2 ;

	cmesh->SetDimModel ( dim );



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

    REAL factor;
    TPZManVector<REAL, 3> bodyforce ( 3,0. );
    factor=1.;
    bodyforce[1]=-20.;

   // material->SetPlasticity ( plasticstep );

    material->SetId ( 1 );

    material->SetWhichLoadVector ( 0 ); //option to compute the total internal force vecor fi=(Bt sigma+ N (b+gradu))

    material->SetLoadFactor ( factor );

    material->SetBodyForce ( bodyforce );
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
    auto * BCond0 = material->CreateBC ( material, -1, directionaldirichlet, val1, val2 );

	val2[0]=1;
	val2[1]=0;
    auto * BCond1 = material->CreateBC ( material, -2, directionaldirichlet, val1, val2 );

	val2[0]=1;
	val2[1]=0;
	auto * BCond2 = material->CreateBC ( material, -5, directionaldirichlet, val1, val2 );

	cmesh->InsertMaterialObject ( BCond0 );
    cmesh->InsertMaterialObject ( BCond1 );
	cmesh->InsertMaterialObject ( BCond2 );

    cmesh->SetAllCreateFunctionsContinuousWithMem();
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

typedefanal *  Slope::SlopeAnalysis(int type,int numthreads)
{

    typedefanal *analysis = new typedefanal(fCompMesh,cout);
   // analysis->SetCompMesh(fCompMesh,true);
switch(type) {
    case 0:
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
    case 1:
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
    return analysis;
}



void Slope::PostPlasticity(string vtkd)
{

   TPZPostProcAnalysis * postprocdeter = new TPZPostProcAnalysis();
    CreatePostProcessingMesh ( postprocdeter);


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
        TPZFStructMatrix<REAL> structmatrix ( PostProcess->Mesh() );
        PostProcess->SetStructuralMatrix ( structmatrix );
        PostProcess->SetPostProcessVariables ( PostProcMatIds, PostProcVars );


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
