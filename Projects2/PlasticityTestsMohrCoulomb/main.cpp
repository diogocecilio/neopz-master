
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
#include "Plasticity/TPZPlasticStep.h"
#include "Plasticity/TPZThermoForceA.h"
#include "Plasticity/TPZPlasticStepVoigt.h"
#include "Plasticity/TPZVonMises.h"
#include "Plasticity/TPZYCMohrCoulombPV.h"
#include "Plasticity/TPZMatElastoPlastic2D.h"

#include "Plasticity/TPZYCMohrCoulombPV2.h"

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
#include "TPZHWTools.h"
//#include "helpers.h"
#include "helpers2.h"
using namespace std;

void SolveFoot()
{
    // TPZGeoMesh* gmesh = ReadGiDMesh("/home/diogo/projects/neopz-master/Projects2/PlasticityTestsMohrCoulomb/mesh_footing-tri.msh", /*mat2D=*/1, /*mat1D=*/-1);
    TPZGeoMesh* gmesh = ReadGiDMesh("/home/diogo/projects/neopz-master/Projects2/PlasticityTestsMohrCoulomb/mesh_footing-tri.msh", /*mat2D=*/1, /*mat1D=*/-1);
    auto cmesh= CreateCMeshFoot(gmesh,2);
    std::string vtkfile2="Foot.vtk";
    int loadid=-4;
    int iref=5;
    STATE tolref=0.01;
    Solve(cmesh,loadid,vtkfile2,iref,tolref);


    TPZElastoPlasticAnalysis an(cmesh, std::cout,TPZElastoPlasticAnalysis::ELineSearch::Armijo);
    TPZSkylineStructMatrix<STATE> matskl(cmesh);
    matskl.SetNumThreads(16);
    an.SetStructuralMatrix(matskl);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    an.SetSolver(step);

    int nsteps=30;
    STATE lambda0=0.0001;
    STATE L0=0.01;
    STATE x=0.;
    STATE y=5;
    IterativeProcessArcLength2(an,nsteps,lambda0,L0,vtkfile2,loadid,x,y);
}

void SolveSlope()
{
    auto gmesh = TriGMesh(1);
    auto cmesh = CreateCMesh(gmesh,2);
    std::string vtkfile2="Slope.vtk";
    int loadid=1;
    int iref=5;
    STATE tolref=0.015;
    Solve(cmesh,loadid,vtkfile2,iref,tolref);


    TPZElastoPlasticAnalysis an(cmesh, std::cout,TPZElastoPlasticAnalysis::ELineSearch::Armijo);
    TPZSkylineStructMatrix<STATE> matskl(cmesh);
    matskl.SetNumThreads(16);
    an.SetStructuralMatrix(matskl);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    an.SetSolver(step);


    int nsteps=30;
    STATE lambda0=0.0001;
    STATE L0=0.3;
    STATE x=30.;
    STATE y=45;
    IterativeProcessArcLength2(an,nsteps,lambda0,L0,vtkfile2,loadid,x,y);
}

int main()
{

    std::cout << "HELLO WORLD"<<std::endl;
    SolveSlope();
    //SolveFoot();


   //  STATE phi=20*M_PI/180.;
   //  STATE psi=phi;
   //  STATE c =50.;
   //  TPZElasticResponse ER;
   //  //ER.SetEngineeringData(0.1e8,0.48);
   //  ER.SetEngineeringData(20000,0.49);
   //  auto mc = TPZYCMohrCoulombPV2( phi, psi, c,ER) ;
   //  TPlasticStepVoigtMC PlasticStepVoigt;
   //  PlasticStepVoigt.SetPlasticCriterion(mc);
   //  PlasticStepVoigt.SetElasticResponse(ER);
   //
   //
   //
   //  PlasticStepVoigt.SetPlasticCriterion(mc);
   //  PlasticStepVoigt.SetElasticResponse(ER);
   //  PlasticStepVoigt.Print(std::cout);
   //
   //
   //  PlasticStepVoigt.Print(std::cout);
   //  TPZFMatrix<REAL> Dep(6,6,0.);
   //  TPZTensor<STATE> epst,epsteste;
   //  TPZTensor<STATE> sigma;
   //  //TPZManVector<STATE,6> epstm={-7.1436465404952029e-5, 3.7632073167565092e-4, 0,-5.4913608546659984e-5 , 0, 1.1423294530637051e-4};
   //  TPZManVector<STATE,6> epstm = {3.6157562487484481e-2, (2.5111726131341262e-2), 0, -2.5800526093809846e-2, 0, 0};
   //  TPZManVector<STATE,6> epsp={0,0,0,0,0,0};
   //  //TPZPlasticState ps;
   //  for(int i=0;i<6;i++)PlasticStepVoigt.fN.EpsP()[i]=epsp[i];
   //  for(int i=0;i<6;i++)epst[i]=epstm[i];
   //
   //
   //  ER.ComputeStress(epst,sigma);
   //  std::cout << sigma <<std::endl;
   //  ER.ComputeStrain(sigma, epsteste);
   //  std::cout << epsteste <<std::endl;
   //  PlasticStepVoigt.ApplyStrainComputeSigma(epst,sigma,&Dep);
   //  std::cout << sigma <<std::endl;
   // Dep.Print("Dep");
   //


   //  TPZGeoMesh* gmesh = ReadGiDMeshSimple("/home/diogo/projects/neopz-master/Projects2/PlasticityTestsMohrCoulomb/mesh_footing.msh", /*mat2D=*/1, /*mat1D=*/-1);
   //  auto cmesh= CreateCMeshFoot(gmesh,2);
   //  //auto cmesh= CreateCMesh2(gmesh,2);
   //  //TPZManVector<REAL> factors={-0.000157143, -0.000357143, -0.000714286, -0.00107143, -0.00142857, -0.00178571, -0.00214286, -0.0025};
   // //  TPZManVector<REAL> factors={-1, -491, -981, -1471, -1961, -2451, -2941, -3431, -3921, -4411, \
   // //      -4901, -5391, -5881, -6371, -6861, -7351, -7841, -8331};
   // //  int loaddir=1;
   // //  int indexbc=-4;
   // //  std::string vtkfile ="postfooting.vtk";
   // // ApplyLoadFoot(cmesh,factors,loaddir, indexbc,vtkfile);
   //
   //  TPZElastoPlasticAnalysis an(cmesh, std::cout,TPZElastoPlasticAnalysis::ELineSearch::Armijo);
   //
   //
   //
   //  //TPZFStructMatrix<STATE> str(cmesh);
   //  TPZSkylineStructMatrix<STATE> str(cmesh);
   //  an.SetStructuralMatrix(str);
   //  TPZStepSolver<REAL> direct;
   //  direct.SetDirect(ELDLt);
   //  //direct.SetDirect(ELU);
   //  an.SetSolver(direct);
   //  int nsteps=30;
   //  STATE lambda0=1.31;
   //  STATE L0=0.01;
   //  std::string vtkfile2="IterativeProcessArcLength2.vtk";
   //  IterativeProcessArcLength2(an,nsteps,lambda0,L0,vtkfile2);

     // TPZGeoMesh* gmesh = ReadGiDMeshSimple("/home/diogo/projects/neopz-master/Projects2/PlasticityTestsMohrCoulomb/mesh_footing.msh", /*mat2D=*/1, /*mat1D=*/-1);
     // auto cmesh= CreateCMeshFoot(gmesh,2);
     //
     // TPZElastoPlasticAnalysis an(cmesh, std::cout,TPZElastoPlasticAnalysis::ELineSearch::Armijo);
     // TPZSkylineStructMatrix<STATE> str(cmesh);
     // an.SetStructuralMatrix(str);
     // TPZStepSolver<REAL> direct;
     // direct.SetDirect(ELDLt);
     // //direct.SetDirect(ELU);
     // an.SetSolver(direct);
     //
     // STATE L0=1;
     // std::string vtkfile2="BissectionFoot.vtk";
     // int loadid=-4;
     // Solve(cmesh,loadid,vtkfile2);




    // auto gmesh = TriGMesh(1);
    // auto cmesh = CreateCMesh(gmesh,2);
    // std::string vtkfile2="Slope.vtk";
    //  int loadid=1;
    //  Solve(cmesh,loadid,vtkfile2);
    //
    //
    // TPZElastoPlasticAnalysis an(cmesh, std::cout,TPZElastoPlasticAnalysis::ELineSearch::Armijo);
    // TPZSkylineStructMatrix<STATE> matskl(cmesh);
    // matskl.SetNumThreads(16);
    // an.SetStructuralMatrix(matskl);
    // TPZStepSolver<STATE> step;
    // step.SetDirect(ELDLt);
    // an.SetSolver(step);
    //
    //
    // int nsteps=30;
    // STATE lambda0=0.0001;
    // STATE L0=0.3;
    // IterativeProcessArcLength2(an,nsteps,lambda0,L0,vtkfile2,loadid);
/*
    //RunDeterministic();
     int porder=2;
     int ref0=0;
     int maxref=5;
  //   Solve( porder, ref0, maxref);

      TPZBFileStream in;
   //  //
    string fin = "gmesh.bin";
    TPZGeoMesh gmesh;
    in.OpenRead(fin);
    gmesh.Read(in, 0);
    //gmesh.Print(cout);

   // int porder=2;
   // auto *gmesh = TriGMesh(3);
    auto cmesh = CreateCMesh(&gmesh,porder);



    TPZElastoPlasticAnalysis an(cmesh, std::cout,TPZElastoPlasticAnalysis::ELineSearch::Armijo);

    //RunAndAccept(cmesh,0.1);
    RunAndAccept(cmesh,0.5);
    RunAndAccept(cmesh,0.6);
    RunAndAccept(cmesh,0.7);
    RunAndAccept(cmesh,0.8);


    // TPZFStructMatrix<STATE> str(cmesh);
    // an.SetStructuralMatrix(str);
    // TPZStepSolver<REAL> direct;
    // direct.SetDirect(ELU);
    // an.SetSolver(direct);
    //
    // int nsteps=10;
    // STATE lambda0=2;
    // STATE L0=0.01;
    // std::string vtkfile="slope-NewtonRaphson.vtk";
    //
    // auto* bodymat = dynamic_cast<TMatElastoPlaticMC*>(cmesh->FindMaterial(1));
    //
    // TPZManVector<REAL,3> fb=bodymat->GetBodyForce0();
    // fb[1]*=0.1;
    // bodymat->SetBodyForce(fb);
    //
    // an.NewtonRaphson();
    // //IterativeProcessArcLength(an,nsteps,lambda0,L0,vtkfile);*/


    return 0;
}



