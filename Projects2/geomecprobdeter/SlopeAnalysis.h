#ifndef SLOPEANALYSIS_H
#define SLOPEANALYSIS_H

// Inclui várias bibliotecas de análise e estruturas para elementos finitos e material elastoplástico
#include "tpzgeoelrefpattern.h"
#include "Plasticity/pzelastoplasticanalysis.h"
//#include "Plasticity/TPZElasticResponse.h"
//#include "Plasticity/TPZYCMohrCoulombPV.h"
//#include "Plasticity/TPZMatElastoPlastic2D.h"
//#include "Plasticity/TPZMatElastoPlastic.h"
//#include "Plasticity/TPZPlasticStepPV.h"
#include <pzgmesh.h> // para TPZGeoMesh, malha geométrica
#include <pzcmesh.h> // para TPZCompMesh, malha computacional
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
#include "TPZVTKGeoMesh.h"
#include <random>
#include "TPZSavable.h"
#include "pzfstrmatrix.h"
using namespace std;


class SlopeAnalysis {
public:
    // Construtor padrão
    SlopeAnalysis();

    // Construtor de cópia
    SlopeAnalysis(const SlopeAnalysis &other);

    SlopeAnalysis(TPZGeoMesh * gmesh, TPZCompMesh * cmesh,int nThreads=10,int solver=0);

    // Destrutor para liberar memória alocada dinamicamente
    ~SlopeAnalysis();

    /**
     * Aplica uma força gravitacional na malha computacional.
     * @param bodyforce Vetor de força de corpo para gravidade.
     */
    void ApplyGravityLoad(TPZManVector<REAL, 3> bodyforce);

    /**
     * Aplica um fator multiplicador à força gravitacional.
     * @param factor Fator multiplicador da força gravitacional.
     */
    void LoadingRamp(REAL factor);

    /**
     * Realiza a redução de resistência até atingir a falha, usando o critério de Mohr-Coulomb.
     * @param maxcount Número máximo de iterações para redução de resistência.
     * @param FS0 Fator de segurança inicial.
     * @param fstol Tolerância para o fator de segurança.
     * @return O fator de segurança encontrado.
     */
    REAL  computelamda ( TPZFMatrix<REAL>& dwb, TPZFMatrix<REAL>& dws, TPZFMatrix<REAL>& dw, REAL& l );
    TPZVec<REAL> computelamdacris ( TPZFMatrix<REAL>& dwb, TPZFMatrix<REAL>& dws, TPZFMatrix<REAL>& dw, REAL& l );
    REAL  computelamda0 ( TPZFMatrix<REAL>& dwb,  TPZFMatrix<REAL>& fext, REAL& l );
    REAL GravityIncrease ( );
    REAL ShearRedNoIntegrationPoints ( int maxcout,REAL FS0,REAL fstol );
    REAL ShearRed ( int maxcout,REAL FS0,REAL fstol );
    REAL IterativeProcessArcLength ( REAL tol,int numiter,REAL tol2,int numiter2,REAL l,REAL lambda0,bool &converge );

    REAL ArcLength(bool &conv);

    void ShearReductionIntegrationPoints ( REAL FS );

    void InitializeMemory ( REAL coesion, REAL atrito);

    //void InitializeMemory();
    /**
     * Resolve o problema elastoplástico determinístico usando valores médios dos parâmetros materiais.
     * @return O fator de segurança calculado.
     */
    REAL SolveDeterministic(bool IsSRM,REAL coes,REAL atrito);

    /**
     * Define as propriedades da análise, como o solver.
     * @return Objeto de análise elastoplástica configurado.
     */
    TPZElastoPlasticAnalysis SetSlopeAnalysis();

    /**
     * Refina a malha geométrica com base em um valor específico.
     * @param refineaboveval Valor para decidir se o refinamento é necessário.
     * @param elindices Índices dos elementos que precisam de refinamento.
     */
    void DivideElementsAbove(REAL refineaboveval, std::set<long> &elindices);

    /**
     * Aumenta o grau de polinômios de interpolação com base em um valor específico.
     * @param refineaboveval Valor para decidir se o refinamento é necessário.
     * @param porder Nova ordem de polinômio para refinamento.
     * @param elindices Índices dos elementos que precisam de refinamento.
     */
    void PRefineElementsAbove(REAL refineaboveval, int porder, std::set<long> &elindices);

    /**
     * Calcula a deformação nos elementos da malha para fins de refinamento.
     */
    void ComputeElementDeformation();

    /**
     * Realiza o pós-processamento e gera visualizações.
     * @param vtkd Caminho para salvar o arquivo de visualização (formato VTK).
     */
    void PostPlasticity(std::string vtkd);

    /**
     * Cria uma malha para pós-processamento.
     * @param PostProcess Ponteiro para o objeto de pós-processamento.
     */
    void CreatePostProcessingMesh(TPZPostProcAnalysis *PostProcess);

    /**
     * Especifica as variáveis para pós-processamento.
     * @param scalNames Nomes das variáveis escalares a serem processadas.
     * @param vecNames Nomes das variáveis vetoriais a serem processadas.
     */
    void PostProcessVariables(TPZStack<std::string> &scalNames, TPZStack<std::string> &vecNames);


    int ClassId() const; // Identificador de classe para serialização

private:

    TPZCompMesh *fCompMesh;   // Ponteiro para a malha computacional principal
    TPZGeoMesh *fGMesh;       // Ponteiro para a malha geométrica principal

    TPZVec<REAL> fPlasticDeformSqJ2; // Vetor de deformação plástica acumulada (sqrt(j2))
    TPZVec<TPZFMatrix<REAL>> fPesos; // Vetor de matrizes de pesos dos campos estocásticos
    int fRef0;             // Nível inicial de refinamento da malha
    int fPorder;           // Ordem dos polinômios de interpolação
    int fNumThreads;       // Número de threads para paralelização
    int fSolver;           // Tipo de solver configurado para a análise

    typedef TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse> plasticmorh;

    typedef TPZMatElastoPlastic2D <TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse>, TPZElastoPlasticMem> plasticmat;

    typedef TPZElastoPlasticAnalysis typedefanal;

};

#endif // SLOPEANALYSIS_H
