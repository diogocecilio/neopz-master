
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

using namespace std;

typedef TPZPlasticStepVoigt<TPZYCMohrCoulombPV2, TPZElasticResponse> TPlasticStepVoigtMC;

TPZManVector<REAL, 3> ComputePrincialVec(TPZTensor<STATE> &tensor)
{
    TPZTensor<REAL>::TPZDecomposed eigen_system;
    tensor.EigenSystem(eigen_system);
    return eigen_system.fEigenvalues;
}
#include "TPZHWTools.h"
int main()
{
    //subst={20000,0.49,20 Pi/180,50};
    // STATE phi=20*M_PI/180.;
    // STATE psi=phi;
    // STATE c =50.;
    // TPZElasticResponse ER;
    // ER.SetEngineeringData(20000.,0.49);

    STATE phi=20*M_PI/180.;
    STATE psi=phi;
    STATE c =490.;
    TPZElasticResponse ER;
    ER.SetEngineeringData(10000000.,0.48);
    auto mc = TPZYCMohrCoulombPV2( phi, psi, c,ER) ;


    TPZTensor<STATE>  sigmatr;
    TPZManVector<STATE,2> dlambda;
    STATE havarn=0.;
    STATE havarn1;
    TPZFNMatrix<9> Grad3x3;
    int myype;


    //TPZFMatrix<STATE> epst={3.6157562487484481e-2, 2.5111726131341262e-2, 0, -2.5800526093809846e-2, 0, 0};//main
   //TPZFMatrix<STATE> epst={2.1153582845979804e-2, -3.2774022989690031e-2, 0, -6.1244636739267791e-3, 0, 0};//rigth
   TPZFMatrix<STATE> epst={-7.1436465404952029e-5, 3.7632073167565092e-4, 0, -5.4913608546659984e-5, 0, 1.1423294530637051e-4};//left
   //TPZFMatrix<STATE> epst={-0.0103988, -0.0157616, -0.0259828, 0.0436286, -0.00650411, 0.0206444};//apex
    TPZFMatrix<STATE> epsp={0,0,0,0,0,0};
   // TPZFMatrix<STATE> epstr=epst-epsp;
    TPZTensor<STATE> epstr,sigprojtensor;
    for(int i=0;i<6;i++)epstr[i]=epst(i,0)-epsp(i,0);
    ER.ComputeStress(epstr,sigmatr);

    ER.ComputeStrain(sigmatr,epstr);

    cout << "epstr = "<< epstr <<endl;

    TPZManVector<STATE,3> sigprojvec,epstrvecout;
    TPZManVector<STATE,3> sigtrvec =ComputePrincialVec(sigmatr);
    TPZManVector<STATE,3> epstrvec =ComputePrincialVec(epstr);


    mc.ProjectSigma(sigtrvec,havarn,dlambda,sigprojvec,epstrvecout,Grad3x3,havarn1,myype);

    cout << " dlambda"<<dlambda <<endl;
    cout << " sigprojvec"<<sigprojvec <<endl;
    cout << " sigtrvec"<<sigtrvec <<endl;
    cout << " epstrvec"<<epstrvec <<endl;
    cout << " epstrvecout"<<epstrvecout <<endl;
    cout << " Grad3x3"<<Grad3x3 <<endl;
    TPZManVector<STATE,3> cyl(3);
    TPZHWTools::FromPrincipalToHWCyl(sigprojvec,cyl);
    cout << " cyl = "<<cyl <<endl;
/*
    TPZFMatrix<STATE> epst={3.6157562487484481e-2, 2.5111726131341262e-2, 0, -2.5800526093809846e-2, 0, 0};
    TPZFMatrix<STATE> epsp={0,0,0,0,0,0};*/


    TPlasticStepVoigtMC PlasticStepVoigt;
    PlasticStepVoigt.SetPlasticCriterion(mc);
    PlasticStepVoigt.SetElasticResponse(ER);
    // void TPZPlasticStepVoigt<YC,ER>::ApplyStrainComputeSigma(const TPZTensor<REAL>& epsTotal,
    //                                                          TPZTensor<REAL>& sigma,
    //                                                          TPZFMatrix<REAL>* tangent)
    TPZFMatrix<STATE> Dep;
    PlasticStepVoigt.ApplyStrainComputeSigma(epstr,sigprojtensor,&Dep);


    return 0;
}
