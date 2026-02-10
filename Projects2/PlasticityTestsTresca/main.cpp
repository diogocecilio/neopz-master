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
#include "Plasticity/TPZYCTrescaVoigt.h"
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
#include "helper.h"
using std::cout;
using std::endl;


std::string meshfile = "/home/diogo/projects/neopz-master/Projects2/PlasticityTestsTresca/gmesh.msh";
#ifdef PZ_LOG
static TPZLogger logger_plasticity ( "PlasticityTests" );
#endif

TMatElastoPlaticVoigtVM* CreateMatVonMises(int id,
                                           STATE E,
                                           STATE nu,
                                           STATE sigy,
                                           STATE H)
{
        auto mat = new TMatElastoPlaticVoigtVM(id);

        TPZElasticResponse ER;
        ER.SetEngineeringData(E, nu);

        TPZYCVonMisesVoigt vmyc;
        const STATE sigmaY0 = sigy;
        const STATE Hiso    = H;
        vmyc.SetUp(sigmaY0, Hiso, ER);

        TPlasticStepVoigtVM PlasticStepVoigt;
        PlasticStepVoigt.SetPlasticCriterion(vmyc);
        PlasticStepVoigt.SetElasticResponse(ER);

        mat->SetPlasticityModel(PlasticStepVoigt);
        mat->SetId(id);

        return mat;
}

TMatElastoPlaticVoigtTresca* CreateMatTresca(int id,
                                             STATE E,
                                             STATE nu,
                                             STATE sigy,
                                             STATE H)
{
        auto mat = new TMatElastoPlaticVoigtTresca(id);

        TPZElasticResponse ER;
        ER.SetEngineeringData(E, nu);

        TPZYCTrescaVoigt vmyc;
        const STATE sigmaY0 = sigy;
        const STATE Hiso    = H;
        vmyc.SetUp(sigmaY0, Hiso, ER);

        TPlasticStepVoigtTresca PlasticStepVoigt;
        PlasticStepVoigt.SetPlasticCriterion(vmyc);
        PlasticStepVoigt.SetElasticResponse(ER);

        mat->SetPlasticityModel(PlasticStepVoigt);
        mat->SetId(id);

        return mat;
}

TMatElastoPlaticMC* CreateMatMohr(int id,
                                  STATE E,
                                  STATE nu,
                                  STATE coes,
                                  STATE phi)
{
        auto mat = new TMatElastoPlaticMC(id);

        TPZElasticResponse ER;
        ER.SetEngineeringData(E, nu);

        TPZYCMohrCoulombPV2 vmyc;
        vmyc.SetUp(phi, phi, coes, ER);

        TPlasticStepVoigtMC PlasticStepVoigt;
        PlasticStepVoigt.SetPlasticCriterion(vmyc);
        PlasticStepVoigt.SetElasticResponse(ER);

        mat->SetPlasticityModel(PlasticStepVoigt);
        mat->SetId(id);

        return mat;
}

template <class TMat>
void InsertFootLoadBCs(TMat* mat, TPZCompMesh* cmesh,int bctype)
{
        TPZFMatrix<STATE> val1(2, 2, 0.0);
        TPZManVector<STATE, 2> val2(2, 0.0);

        const int neumannDirichletType = 3; // teu "dir"

        // Apoio em -1 (fixo em x)
        val2[0] = 1.0;
        val2[1] = 0.0;
        auto* bc0 = mat->CreateBC(mat, -1, neumannDirichletType, val1, val2);

        // Apoio em -2 (fixo em y)
        val2[0] = 0.0;
        val2[1] = 1.0;
        auto* bc1 = mat->CreateBC(mat, -2, neumannDirichletType, val1, val2);

        // Apoio em -3 (fixo em x)
        val2[0] = 1.0;
        val2[1] = 0.0;
        auto* bc2 = mat->CreateBC(mat, -3, neumannDirichletType, val1, val2);

        //Carga em -4 (tipo 1, sem valor prescrito inicial)
        val2[0] = 0.0;
        val2[1] = 0.0;
        auto* bc4 = mat->CreateBC(mat, -4, bctype, val1, val2);

        // val1(1,1)=1.;
        // val2[0] = 0.0;
        // val2[1] = 0.0;
        // auto* bc4 = mat->CreateBC(mat, -4, 7, val1, val2);





        cmesh->InsertMaterialObject(bc0);
        cmesh->InsertMaterialObject(bc1);
        cmesh->InsertMaterialObject(bc2);
        cmesh->InsertMaterialObject(bc4);
}

enum class EMat {EMohr=1,EVonMises=2,ETresca=3};

TPZCompMesh* CreateCMeshFoot(TPZGeoMesh* gmesh,
                                 int pOrder,
                                 EMat type,int bctype)
{
        auto* cmesh = new TPZCompMesh(gmesh);
        cmesh->SetDefaultOrder(pOrder);
        cmesh->SetDimModel(2);

        const int   id   = 1;
        const STATE E    = 1.e7;
        const STATE nu   = 0.48;
        const STATE sigy = 848.7;
        const STATE H    = 0.0;
        const STATE phi  = 20.0 * M_PI / 180.0;
        const STATE coes = sigy / std::sqrt(3.0);

        switch (type)
        {
                case EMat::EMohr:
                {
                        auto* matMohr = CreateMatMohr(id, E, nu, coes, phi);
                        cmesh->InsertMaterialObject(matMohr);
                        InsertFootLoadBCs(matMohr, cmesh,bctype);
                        break;
                }
                case EMat::EVonMises:
                {
                        auto* matVM = CreateMatVonMises(id, E, nu, sigy, H);
                        cmesh->InsertMaterialObject(matVM);
                        InsertFootLoadBCs(matVM, cmesh,bctype);
                        break;
                }
                case EMat::ETresca:
                {
                        auto* matTresca = CreateMatTresca(id, E, nu, 2*coes, H);
                        cmesh->InsertMaterialObject(matTresca);
                        InsertFootLoadBCs(matTresca, cmesh,bctype);
                        break;
                }
                default:
                        DebugStop(); // ou lança exceção
        }

        cmesh->SetAllCreateFunctionsContinuousWithMem();
        cmesh->AutoBuild();
        cmesh->AdjustBoundaryElements();
        cmesh->CleanUpUnconnectedNodes();

        return cmesh;
}


void SolveFoot(EMat mattype)
{

        TPZGeoMesh* gmesh = ReadGiDMesh ( meshfile, /*mat2D=*/1, /*mat1D=*/-1 );
        int bctype=1;
        auto cmesh= CreateCMeshFoot ( gmesh,2,mattype,bctype );
        std::cout << "Equations = " << cmesh->NEquations() << std::endl;
        int indexfoot=-4;
        std::string namevtk="refine.vtk";
        int ref=4;
        STATE tol_fs_rel =0.0001;
        Solve ( cmesh,indexfoot,namevtk, ref, tol_fs_rel );
        TPZGeoMesh *gmesh2 = cmesh->Reference(); // sua rotina de geração
        WriteGeoMesh ( gmesh2, "malha_refinada.txt" );
}

void SolveFootApplyDisplacement(std::string namevtk,EMat mattype)
{
        TPZGeoMesh* gmesh = ReadGiDMesh ( meshfile, /*mat2D=*/1, /*mat1D=*/-1 );
        //TPZGeoMesh *gmesh = ReadGeoMesh ( "malha_refinada.txt" );
        int bctype=0;
        auto cmesh= CreateCMeshFoot ( gmesh,2 ,mattype, bctype);
        std::cout << "Equations = " << cmesh->NEquations() << std::endl;
       TPZManVector<REAL,14> factors= { -0.0001, -0.0001, -0.0001,  -0.0001,  -0.0001,  -0.0001,  -0.0001,  -0.0001,  -0.0001,  -0.0001,  -0.0001,  -0.0001,  -0.0001,  -0.0001,  -0.0001,  -0.0001,  -0.0001,  -0.0001,  -0.0001,  -0.0001};
       // TPZManVector<REAL,14> factors={-0.500000e-03, -0.500000e-03, -0.250000e-03, -0.250000e-03,-0.500000e-03,-0.500000e-03};
       //factors*=2;
        ApplyLoad2 ( cmesh,factors,1,-4,namevtk );

}
void SolveFootArcLength(std::string vtkfile,EMat type)
{

        TPZGeoMesh* gmesh3 = ReadGiDMesh ( meshfile, /*mat2D=*/1, /*mat1D=*/-1 );
        //TPZGeoMesh *gmesh3 = ReadGeoMesh ( "malha_refinada.txt" );
        int bctype=1;
        auto cmesh3= CreateCMeshFoot ( gmesh3,2 , type,bctype);
        std::cout << "Equations = " << cmesh3->NEquations() << std::endl;
        TPZElastoPlasticAnalysis anal ( cmesh3, std::cout,TPZElastoPlasticAnalysis::ELineSearch::Armijo );


        TPZSkylineStructMatrix<STATE> matskl ( cmesh3 );
        matskl.SetNumThreads ( 16 );
        anal.SetStructuralMatrix ( matskl );
        TPZStepSolver<STATE> step;
        step.SetDirect ( ELDLt );
        anal.SetSolver ( step );


        int itersout;
        int nsteps=100;
        int loaddir=1;
        STATE bccondval=-490;
        STATE lambda0=0.001;
        STATE L0=0.0001;

        REAL out = IterativeProcessArcLength ( anal,loaddir,-4,nsteps, bccondval,lambda0,L0,vtkfile );

}
int main()
{
        // SolveFoot();
        // std::string vtkfile="foot-arclength-mc.vtk";
        // SolveFootArcLength(vtkfile,EMat::EMohr);
        std::string namevtk="foot-diplacement-mc.vtk";
        SolveFootApplyDisplacement(namevtk,EMat::EMohr);

        // std::string vtkfile="foot-arclength-vm.vtk";
        // SolveFootArcLength(vtkfile,EMat::EVonMises);
        // std::string namevtk="foot-diplacement-mises.vtk";
        // SolveFootApplyDisplacement(namevtk,EMat::EVonMises);
//
        // std::string vtkfile="foot-arclength-tresca.vtk";
        // SolveFootArcLength(vtkfile,EMat::ETresca);
        // std::string namevtk="foot-diplacement-tresca.vtk";
        // SolveFootApplyDisplacement(namevtk,EMat::ETresca);
        //SolveFootApplyDisplacement();
        //SolveCyl();


        return 0;
}
