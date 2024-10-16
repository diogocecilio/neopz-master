#ifndef SLOPEANALYSIS_H
#define SLOPEANALYSIS_H
#include "tpzgeoelrefpattern.h"
#include "Plasticity/pzelastoplasticanalysis.h"
#include "Plasticity/TPZElasticResponse.h"
#include "Plasticity/pzelastoplasticanalysis.h"
#include "Plasticity/TPZYCMohrCoulombPV.h"
//#include "Plasticity/TPZMatElastoPlastic_impl.h"
#include "Plasticity/TPZMatElastoPlastic2D.h"
#include "Plasticity/TPZMatElastoPlastic.h"
#include "Plasticity/TPZPlasticStepPV.h"
#include <pzgmesh.h> //for TPZGeoMesh
#include <pzcmesh.h> //for TPZCompMesh
#include <time.h>

#include <cmath>
#include <set>
#include <iostream>
#include <fstream>
#include <ctime>
#include <ratio>
#include <chrono>
#include "pzinterpolationspace.h"
#include "pzskylstrmatrix.h"
#include <TPZSSpStructMatrix.h>
#include "TPZRandomFieldAnalysis.h"
#include "TPZPardisoSolver.h"
#include "TPZVTKGeoMesh.h"
#include <random>
#include "TPZSavable.h"
#include "pzfstrmatrix.h"
using namespace std;
typedef TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse> plasticmorh;
typedef   TPZMatElastoPlastic2D < TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse>, TPZElastoPlasticMem > plasticmat;
//typedef TPZNonLinearAnalysis typedefanal;
typedef TPZElastoPlasticAnalysis typedefanal;

class SlopeAnalysis{
public:
    /**
     * Default constructor
     */
    SlopeAnalysis();

    /**
     * Copy constructor
     *
     * @param other TODO
     */
    SlopeAnalysis(const SlopeAnalysis& other);

    /**
     *  constructor
     *
     * @param other TODO
     */
    SlopeAnalysis(  REAL gammaagua, REAL gammasolo,REAL coes,REAL atrito,int ref0, int porder);

    /**
     * Destructor
     */
    ~SlopeAnalysis();

    void ApplyGravityLoad(TPZManVector<REAL, 3> bodyforce);

    void LoadingRamp (  REAL factor );

    REAL ShearRed ( int maxcout,REAL FS0,REAL fstol );

    void TransferFieldsSolutionFrom ( int isol );

    void ShearReductionIntegrationPoints(REAL FS);

    void InitializeMemory();

    void InitializePOrderInRegion();

    void Write(TPZStream &buf, int withclassid) const ;

    void Read(TPZStream &buf, void *context) ;

    TPZGeoMesh * TriGMesh ( int ref );

    TPZCompMesh * CreateCMesh( TPZGeoMesh *gmesh, int pOrder, REAL coes,REAL atrito);

    TPZFMatrix<REAL>  GenerateRandomField (REAL mean, REAL cov ,int samples);

    void ManageFieldCretion();

    TPZElastoPlasticAnalysis  SetSlopeAnalysis(int type,int numthreads);

    REAL SolveSingleField(int ifield);

    void DivideElementsAbove(REAL refineaboveval, std::set<long> &elindices);

    void PRefineElementsAbove(REAL refineaboveval, int porder, std::set<long> &elindices);

    void ComputeElementDeformation();
public:
    void PostPlasticity(std::string vtkd);


    void CreatePostProcessingMesh (TPZPostProcAnalysis * PostProcess );

    void PostProcessVariables ( TPZStack<std::string> &scalNames, TPZStack<std::string> &vecNames );

    void SetFieldsData(TPZCompMesh  * CompMeshField,TPZFMatrix<REAL>SolutionValVec,TPZVec<REAL> meanvec,TPZVec<REAL> covvec, int samples)
    {
        fCompMeshField =CompMeshField;
        fSolutionValVec=SolutionValVec;
        fMeanvec =meanvec;
        fCovvec=covvec;
        fSamples=samples;

    }

    void SetFields(TPZVec<TPZFMatrix<REAL>> fields)
    {
        fFields=fields;
    }

    TPZVec<TPZFMatrix<REAL>> GetFields()
    {
        return fFields;
    }
    int ClassId() const;

    void IntegrateFieldOverARegion(int imc);

private:

    REAL fCohesion;
    REAL fAtrito;
    //int  fPorder;
    REAL fGammaW;
    REAL fGammaS;

    TPZCompMesh * fCompMeshField;
    TPZFMatrix<REAL> fSolutionValVec;

    TPZVec<TPZFMatrix<REAL>> fFields;

    TPZVec<REAL> fMeanvec;
    TPZVec<REAL> fCovvec;
    int fSamples;

    TPZCompMesh * fCompMesh;
    TPZGeoMesh * fGMesh;

    	TPZVec<REAL> fPlasticDeformSqJ2;
};

#endif // SLOPEANALYSIS_H
