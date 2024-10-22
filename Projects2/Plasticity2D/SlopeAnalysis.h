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
    SlopeAnalysis(const SlopeAnalysis & other);

    SlopeAnalysis& operator=(const SlopeAnalysis& other);

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

    void TransferFieldsSolutionFrom ( int isol);

    void ShearReductionIntegrationPoints(REAL FS);

    void InitializeMemory();

    void InitializePOrderInRegion();

    void Write(TPZStream &buf, int withclassid) const ;

    void Read(TPZStream &buf, void *context) ;

    TPZGeoMesh * TriGMesh ( int ref );

    TPZCompMesh * CreateCMesh( TPZGeoMesh *gmesh, int pOrder, REAL coes,REAL atrito);

    TPZFMatrix<REAL>  GenerateRandomField (REAL mean, REAL cov,TPZFMatrix<REAL> valvec, TPZFMatrix<REAL> stdnormalsamples);
    TPZVec<TPZFMatrix<REAL>>  GenerateRandomField2 (REAL mean, REAL cov, TPZFMatrix<REAL> valvec);

    void ManageFieldCretion();

    void ManageFieldCretion(std::vector<int> fieldindexes);



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
        fNSamples=samples;

    }

    void SetFields(TPZVec<TPZFMatrix<REAL>> fields)
    {
        fFields=fields;
    }

    void SetFieldsSamples(TPZVec<TPZFMatrix<REAL>> fields)
    {
        fFieldSamples=fields;
    }

    TPZVec<TPZFMatrix<REAL>> GetFieldsSamples()
    {
        return fFieldSamples;
    }

    TPZVec<TPZFMatrix<REAL>> GetFields()
    {
        return fFields;
    }

    void SetHFields(TPZVec<TPZFMatrix<REAL>> fields)
    {
        fHFields=fields;
    }

    TPZVec<TPZFMatrix<REAL>> GetHFields()
    {
        return fHFields;
    }
    int ClassId() const;

    void IntegrateFieldOverARegion(int imc);

    //para monte carlo amostragem por importancia
    bool FindCriticalMonteCarloSimulations(int imc);

    //cria matriz de tamanho MxNsamples e preenche com distribuicao normal padrao
    TPZFMatrix<REAL> CreateNormalStandardSamples( );

    REAL CreateNormalStandardSample( );
    //para monte carlo amostragem por importancia h(x)
    void ComputeH();

    std::vector<std::vector<int>> SelectCriticalIndexes();

   // std::vector<std::vector<int>> SelectCriticalIndexes2();
int CountCriticalFields();
    std::vector<int> GetIndex(std::vector<double> vetor,bool decrescente=true)
    {
        std::vector<int> indices(vetor.size());
        for (int i = 0; i < vetor.size(); ++i) {
            indices[i] = i;
        }

        std::sort(indices.begin(), indices.end(), [&](int a, int b) {
            if(decrescente)
            {
                return vetor[a] > vetor[b]; // Ordena de forma decrescente
            }else{
                return vetor[a] < vetor[b]; // Ordena de forma crescente
            }
        });

        return indices;
    }

    std::vector<std::pair<int, int>> encontrarRepetidos(const std::vector<int>& maiores) {
    std::unordered_map<int, int> repetidos;
    std::vector<std::pair<int, int>> resultado;

    // Contar as ocorrências de índices no vetor 'maiores'
    for (int i : maiores) {
        repetidos[i]++;
    }

    // Adicionar ao resultado apenas os índices que aparecem mais de uma vez
    for (const auto& entry : repetidos) {
        if (entry.second > 1) {
            resultado.push_back({entry.first, entry.second});
        }
    }

    return resultado;
}

private:

    REAL fCohesion;
    REAL fAtrito;
    //int  fPorder;
    REAL fGammaW;
    REAL fGammaS;

    TPZCompMesh * fCompMeshField;
    TPZFMatrix<REAL> fSolutionValVec;

    TPZVec<TPZFMatrix<REAL>> fFields;

    TPZVec<TPZFMatrix<REAL>> fHFields;

     TPZVec<TPZFMatrix<REAL>> fFieldSamples;

    TPZVec<REAL> fMeanvec;
    TPZVec<REAL> fCovvec;
    int fNSamples;

    TPZCompMesh * fCompMesh;
    TPZGeoMesh * fGMesh;

    TPZVec<REAL> fPlasticDeformSqJ2;

    TPZVec<TPZFMatrix<REAL>> fPesos;
};

#endif // SLOPEANALYSIS_H
