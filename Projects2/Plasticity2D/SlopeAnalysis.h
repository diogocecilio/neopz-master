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

    // Construtor de copia
    SlopeAnalysis(const SlopeAnalysis &other);

    // Construtor customizado
    SlopeAnalysis(REAL gammaagua, REAL gammasolo, REAL coes, REAL atrito, int ref0, int porder,int therads,int solver);

    // Destrutor: deleta os ponteiros
    ~SlopeAnalysis();

    //aplica a foca de volume
    void ApplyGravityLoad(TPZManVector<REAL, 3> bodyforce);

    //Aplica um factor multilplicador a forca de volume
    void LoadingRamp(REAL factor);

    //Metodo que reduz a resistencia ate que seja verificada a falha no contexto da elastoplasticidade incremental. Utiliza o Mohr coulomb.
    REAL ShearRed(int maxcout, REAL FS0, REAL fstol);

    //Transfere os valores de coesao e angulo de atrito da malha do campo estocastico para a malha elastoplastica
    void TransferFieldsSolutionFrom(int isol);

    //Reduz a resistencia dos parametros materiais nos pontos de integracao da malha
    void ShearReductionIntegrationPoints(REAL FS);

    //inicializa a memoria plastica da malha de elementos finitos com o valor medio dos parametros materiais (e.g. coesao angulo de atrio e permeabildiade)
    void InitializeMemory();

    //Resolve o prblema de plasticide incremental deterministico, ou seja, co o valor dos parametros materiais iguis a media em todos os pontos de integracao
    REAL SolveDeterministic();

    //Faz a integracao de um parametro material no ponto de integraçao em uma determinada regiao do dominio. Neste metodo deve-se especificar a regiao que se deseja integrar bem como o parametro material que se deseja conhecer o valor medio.
    void IntegrateFieldOverARegion(int imc);

    //cria  fNSamles com distribuicao normal padrao N~(0,1)
    TPZFMatrix<REAL> CreateNormalStandardSamples();

    //cria malha computacional geometrica estruturada com elementos triangulares
    TPZGeoMesh* TriGMesh(int ref);

    //cria malha computacional com material elastoplastico mohr-coulomb com memoria
    TPZCompMesh* CreateCMesh(TPZGeoMesh *gmesh, int pOrder, REAL coes, REAL atrito);

    //cria campo estocastico a partido dos auto valores e das auto funcoes de karhunen loeve
    TPZFMatrix<REAL> GenerateRandomField(REAL mean, REAL cov, TPZFMatrix<REAL> valvec, TPZFMatrix<REAL> stdnormalsamples);

    //cria campo estocastico a partir dos autovalores e autovetores mais a distribuicao N~(0,1) e os atributos dos camps estocasticos
    void ManageFieldCretion();

    //cria campo estocastico a partir dos autovalores e autovetores SELECIONADOS (fieldindexes) mais a distribuicao N~(0,1) e os atributos dos camps estocasticos
    void ManageFieldCretion(std::vector<int> fieldindexes);

    //especifica as propriedas daanalise, como solver
    TPZElastoPlasticAnalysis SetSlopeAnalysis();

    //resolve o problema fazendo a tranferecia da solucao da malha estocastica para a malha elastoplastica
    REAL SolveSingleField(int ifield);

    //faz o refinamento geometrico da malha a depender do parametro especificado. atualemte depende do valor de sqrt(j2)
    void DivideElementsAbove(REAL refineaboveval, std::set<long> &elindices);

    //aumenta o grau das funcoes de interpolacao da malha a depender do parametro especificado. atualemte depende do valor de sqrt(j2)
    void PRefineElementsAbove(REAL refineaboveval, int porder, std::set<long> &elindices);

    //calcula a deformacao nos elementos.metodo utilizado para refinamento da malha, ou seja, amalha precisa da deformacao para calcular sqrt(j2)
    void ComputeElementDeformation();

    //faz o pos processamento
    void PostPlasticity(std::string vtkd);

    //cria a malha de pos processamento
    void CreatePostProcessingMesh(TPZPostProcAnalysis *PostProcess);

    //especifica as variaveis a serem pos processadas
    void PostProcessVariables(TPZStack<std::string> &scalNames, TPZStack<std::string> &vecNames);

    //especifica os auto valores autovetores, a malha do campo estocastido e os valores da media e coeficiente de variacao para a geraoca do campo estocasito.
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



    void Write(TPZStream &buf, int withclassid) const;

    void Read(TPZStream &buf, void *context);


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
