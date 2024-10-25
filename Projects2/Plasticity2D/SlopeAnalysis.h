#ifndef SLOPEANALYSIS_H
#define SLOPEANALYSIS_H

#include "tpzgeoelrefpattern.h"
#include "Plasticity/pzelastoplasticanalysis.h"
#include "Plasticity/TPZElasticResponse.h"
#include "Plasticity/pzelastoplasticanalysis.h"
#include "Plasticity/TPZYCMohrCoulombPV.h"
#include "Plasticity/TPZMatElastoPlastic2D.h"
#include "Plasticity/TPZMatElastoPlastic.h"
#include "Plasticity/TPZPlasticStepPV.h"
#include <pzgmesh.h> // for TPZGeoMesh
#include <pzcmesh.h> // for TPZCompMesh
#include <time.h>
#include <numeric>
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
typedef TPZMatElastoPlastic2D <TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse>, TPZElastoPlasticMem> plasticmat;
typedef TPZElastoPlasticAnalysis typedefanal;

class SlopeAnalysis {
public:
    // Construtor padrão
    SlopeAnalysis();

    // Construtor de cópia
    SlopeAnalysis(const SlopeAnalysis &other);



    // Construtor customizado
    SlopeAnalysis(REAL gammaagua, REAL gammasolo, REAL coes, REAL atrito, int ref0, int porder,int therads,int solver);

    // Destrutor
    ~SlopeAnalysis();

    void ApplyGravityLoad(TPZManVector<REAL, 3> bodyforce);
    void LoadingRamp(REAL factor);
    REAL ShearRed(int maxcout, REAL FS0, REAL fstol);
    void TransferFieldsSolutionFrom(int isol);
    void ShearReductionIntegrationPoints(REAL FS);
    void InitializeMemory();
    REAL SolveDeterministic();
    void InitializePOrderInRegion();
    void Write(TPZStream &buf, int withclassid) const;
    void Read(TPZStream &buf, void *context);

    TPZGeoMesh* TriGMesh(int ref);
    TPZCompMesh* CreateCMesh(TPZGeoMesh *gmesh, int pOrder, REAL coes, REAL atrito);

    TPZFMatrix<REAL> GenerateRandomField(REAL mean, REAL cov, TPZFMatrix<REAL> valvec, TPZFMatrix<REAL> stdnormalsamples);
    TPZVec<TPZFMatrix<REAL>> GenerateRandomField2(REAL mean, REAL cov, TPZFMatrix<REAL> valvec);

    void ManageFieldCretion();
    void ManageFieldCretion(std::vector<int> fieldindexes);
    TPZElastoPlasticAnalysis SetSlopeAnalysis();
    REAL SolveSingleField(int ifield);

    void DivideElementsAbove(REAL refineaboveval, std::set<long> &elindices);
    void PRefineElementsAbove(REAL refineaboveval, int porder, std::set<long> &elindices);
    void ComputeElementDeformation();

    void PostPlasticity(std::string vtkd);
    void CreatePostProcessingMesh(TPZPostProcAnalysis *PostProcess);
    void PostProcessVariables(TPZStack<std::string> &scalNames, TPZStack<std::string> &vecNames);

    void SetFieldsData(TPZCompMesh *CompMeshField, TPZFMatrix<REAL> SolutionValVec, TPZVec<REAL> meanvec, TPZVec<REAL> covvec, int samples) {


        fSolutionValVec = SolutionValVec;
        fMeanvec = meanvec;
        fCovvec = covvec;
        fNSamples = samples;
        if ( fCompMeshField == CompMeshField ) {
                return;
        }else{
            fCompMeshField = CompMeshField;
        }
    }

    void SetFields(TPZVec<TPZFMatrix<REAL>> fields) {
        fFields = fields;
    }

    void SetFieldsSamples(TPZVec<TPZFMatrix<REAL>> fields) {
        fFieldSamples = fields;
    }

    TPZVec<TPZFMatrix<REAL>> GetFieldsSamples() {
        return fFieldSamples;
    }

    TPZVec<TPZFMatrix<REAL>> GetFields() {
        return fFields;
    }

    void SetHFields(TPZVec<TPZFMatrix<REAL>> fields) {
        fHFields = fields;
    }

    TPZVec<TPZFMatrix<REAL>> GetHFields() {
        return fHFields;
    }

    int ClassId() const;
    void IntegrateFieldOverARegion(int imc);
    bool FindCriticalMonteCarloSimulations(int imc);
    TPZFMatrix<REAL> CreateNormalStandardSamples();




private:
    REAL fCohesion;
    REAL fAtrito;
    REAL fGammaW;
    REAL fGammaS;

    TPZCompMesh *fCompMeshField; // Ponteiro para a malha computacional do campo
    TPZFMatrix<REAL> fSolutionValVec;

    TPZVec<TPZFMatrix<REAL>> fFields;
    TPZVec<TPZFMatrix<REAL>> fHFields;
    TPZVec<TPZFMatrix<REAL>> fFieldSamples;

    TPZVec<REAL> fMeanvec;
    TPZVec<REAL> fCovvec;
    int fNSamples;

    TPZCompMesh * fCompMesh;       // Usando std::unique_ptr
    TPZGeoMesh * fGMesh;           // Usando std::unique_ptr


    TPZVec<REAL> fPlasticDeformSqJ2;
    TPZVec<TPZFMatrix<REAL>> fPesos;
    int fRef0;
    int fPorder;
    int fNumThreads;
    int fSolver;
};

#endif // SLOPEANALYSIS_H
