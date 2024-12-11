#ifndef SLOPEANALYSIS_H
#define SLOPEANALYSIS_H
#include "TPZFileStream.h"
#include <TPZBFileStream.h>
#include "TPZSavable.h"
// Inclui várias bibliotecas de análise e estruturas para elementos finitos e material elastoplástico
#include "tpzgeoelrefpattern.h"
#include "Plasticity/pzelastoplasticanalysis.h"
#include "Plasticity/TPZElasticResponse.h"
#include "Plasticity/TPZYCMohrCoulombPV.h"
#include "Plasticity/TPZMatElastoPlastic2D.h"
#include "Plasticity/TPZMatElastoPlastic.h"
#include "Plasticity/TPZPlasticStepPV.h"
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
#include "TPZRandomFieldAnalysis.h"
#include "TPZPardisoSolver.h"
#include "TPZVTKGeoMesh.h"
#include <random>
#include "TPZSavable.h"
#include "pzfstrmatrix.h"


#include "pzinterpolationspace.h"
#include "pzstack.h"
#include "pzcmesh.h"
#include "pzquad.h"
#include "TPZMaterial.h"
#include "TPZMatInterfaceSingleSpace.h"
#include "TPZMatInterfaceCombinedSpaces.h"
#include "TPZMatWithMem.h"
#include "pzelctemp.h"
#include "pzmultiphysicscompel.h"

using namespace std;


class SlopeAnalysis {
public:
    // Construtor padrão
    SlopeAnalysis();

    // Construtor de cópia
    SlopeAnalysis(const SlopeAnalysis &other);

    /**
     * Construtor customizado.
     * @param gammaagua Peso específico da água.
     * @param gammasolo Peso específico do solo.
     * @param coes Coesão do solo.
     * @param atrito Ângulo de atrito do solo.
     * @param ref0 Nível de refinamento inicial da malha.
     * @param porder Ordem do polinômio para funções de interpolação.
     * @param threads Número de threads para paralelização.
     * @param solver Tipo de solver a ser usado na análise.
     */
    SlopeAnalysis(REAL gammaagua, REAL gammasolo, REAL coes, REAL atrito, int ref0, int porder, int threads, int solver);

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
    REAL ShearRed(int maxcount, REAL FS0, REAL fstol);

    REAL ShearRedNoIntegrationPoints ( int maxcout,REAL FS0,REAL fstol );
    /**
     * Realiza o aumento gradativo da carga até atingir a falha, usando o critério de Mohr-Coulomb.
     * @return O fator de segurança encontrado.
     */
    REAL GravityIncrease ();

    REAL IterativeProcessArcLength ( REAL tol,int numiter,REAL tol2,int numiter2,REAL l,REAL lambda0,bool &converge );

    REAL  computelamda0 ( TPZFMatrix<REAL>& dwb,  TPZFMatrix<REAL>& dw, REAL& l );

    REAL  computelamda ( TPZFMatrix<REAL>& dwb, TPZFMatrix<REAL>& dws, TPZFMatrix<REAL>& dw, REAL& l );

    TPZVec<REAL>  computelamdacris ( TPZFMatrix<REAL>& dwb, TPZFMatrix<REAL>& dws, TPZFMatrix<REAL>& dw, REAL& l);


    REAL ArcLength(bool &conv);

    /**
     * Transfere valores de coesão e ângulo de atrito da malha de campo estocástico para a malha elastoplástica.
     * @param isol Índice da solução a ser transferida.
     */
    void TransferFieldsSolutionFrom(int isol);

    void TransferFieldsSolutionFrom ( TPZVec<TPZFMatrix<REAL>> sample );

    /**
     * Reduz a resistência dos parâmetros materiais nos pontos de integração.
     * @param FS Fator de segurança atual.
     */
    void ShearReductionIntegrationPoints(REAL FS);

    // Inicializa a memória plástica da malha de elementos finitos com valores médios dos parâmetros materiais
    void InitializeMemory();

    /**
     * Resolve o problema elastoplástico determinístico usando valores médios dos parâmetros materiais.
     * @return O fator de segurança calculado.
     */
    REAL SolveDeterministic(bool IsSRM);

    /**
     * Realiza a integração de um parâmetro material em uma região específica do domínio.
     * @param imc Índice da região onde a integração será realizada.
     */
    void IntegrateFieldOverARegion(int imc);

    /**
     * Realiza a integração de um parâmetro material em uma região específica em que houve plastificacoa.
     * @param refineaboveval valor a de deformacao plastica a cima do qual sera selecionado elemento para integrar field
     * @param imc Índice da região onde a integração será realizada.
     */
    void IntegrateFieldOverARegion ( REAL refineaboveval,int imc );

    void IntegrateFieldOverARegionB ( REAL refineaboveval,int imc );

    /**
     * Gera amostras normais padrão (N~(0,1)) para análise estocástica.
     * @return Matriz de amostras normais padrão.
     */

    TPZFMatrix<REAL> CreateNormalStandardSamples();


    std::vector<std::pair<int, double>> CrudeMonteCarlo(int a,int b);


    /**
     * Cria uma malha geométrica estruturada com elementos triangulares.
     * @param ref Nível de refinamento da malha.
     * @return Ponteiro para a malha geométrica criada.
     */
    TPZGeoMesh* TriGMesh(int ref);

    TPZGeoMesh * QuadGMesh(int ref);
    /**
     * Cria uma malha computacional elastoplástica com material Mohr-Coulomb.
     * @param gmesh Ponteiro para a malha geométrica.
     * @param pOrder Ordem do polinômio de interpolação.
     * @param coes Coesão do material.
     * @param atrito Ângulo de atrito do material.
     * @return Ponteiro para a malha computacional criada.
     */
    TPZCompMesh* CreateCMesh(TPZGeoMesh *gmesh, int pOrder, REAL coes, REAL atrito);

    /**
     * Gera um campo aleatório usando autovalores e autofunções de Karhunen-Loève.
     * @param mean Valor médio do campo.
     * @param cov Coeficiente de variação.
     * @param valvec Matriz de autovalores.
     * @param stdnormalsamples Amostras normais padrão.
     * @return Matriz representando o campo aleatório gerado.
     */
    TPZFMatrix<REAL> GenerateRandomField(REAL mean, REAL cov, TPZFMatrix<REAL> valvec, TPZFMatrix<REAL> stdnormalsamples);

    // Cria campos estocásticos com atributos de campo
    void ManageFieldCretion();
    void ManageFieldCretion(std::vector<int> fieldindexes);
    void ManageFieldCretion ( std::vector<std::vector<int>>  fieldindexes );
    void ManageFieldCretionSubSet();
    /**
     * Define as propriedades da análise, como o solver.
     * @return Objeto de análise elastoplástica configurado.
     */
    TPZElastoPlasticAnalysis SetSlopeAnalysis();

    /**
     * Resolve o problema usando a solução da malha estocástica.
     * @param ifield Índice do campo a ser resolvido.
     * @return O fator de segurança calculado para o campo especificado.
     */
    REAL SolveSingleField(int ifield);

    REAL SolveSingleField(TPZVec<TPZFMatrix<REAL>> sample );

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

    /**
     * Define dados para geração do campo estocástico.
     * @param CompMeshField Ponteiro para a malha do campo estocástico.
     * @param SolutionValVec Vetor de valores da solução.
     * @param meanvec Vetor de médias dos campos.
     * @param covvec Vetor de coeficientes de variação.
     * @param samples Número de amostras do campo.
     */
    void SetFieldsData(TPZCompMesh *CompMeshField, TPZFMatrix<REAL> SolutionValVec, TPZVec<REAL> meanvec, TPZVec<REAL> covvec, int samples);

    // Define e obtém os campos estocásticos e suas amostras
    void SetFields(TPZVec<TPZFMatrix<REAL>> fields);
    void SetFieldsSamples(TPZVec<TPZFMatrix<REAL>> fields);
    TPZVec<TPZFMatrix<REAL>> GetFieldsSamples();
    TPZVec<TPZFMatrix<REAL>> GetFields();

    int ClassId() const; // Identificador de classe para serialização

    void Write(TPZStream &buf, int withclassid) const; // salva a classe
    void Read(TPZStream &buf, void *context); // le a classe



    int GetM()
    {
        return fSolutionValVec.Cols();
    }

    TPZVec<TPZFMatrix<REAL>> GetFieldSamples()
    {
        return fFieldSamples;
    }

    void SetSubSetSamples(TPZVec<TPZFMatrix<REAL>> data)
    {
        fFieldSamplesSubSetFN.Push(data);

    }


    TPZVec<TPZFMatrix<REAL>> GetSubSetSamples(int isample)
    {
        return fFieldSamplesSubSetFN[isample];
    }

    void ResetSubSetSamples()
    {
        fFieldSamplesSubSetFN.resize(0);
    }

    TPZVec<TPZFMatrix<REAL>> GetIfield(int ifield)
    {
        TPZVec<TPZFMatrix<REAL>> field(2);
        field[0].Resize(GetM(),1);
        field[1].Resize(GetM(),1);
        for(int iM=0;iM<GetM();iM++)
        {
            field[0](iM,0)=fFieldSamples[0](iM,ifield);
            field[1](iM,0)=fFieldSamples[1](iM,ifield);
        }
        return field;
    }

private:
    REAL fCohesion;   // Coesão do solo
    REAL fAtrito;     // Ângulo de atrito do solo
    REAL fGammaW;     // Peso específico da água
    REAL fGammaS;     // Peso específico do solo

    TPZCompMesh *fCompMeshField; // Ponteiro para a malha computacional do campo
    TPZFMatrix<REAL> fSolutionValVec; // Matriz de valores da solução do campo

    TPZVec<TPZFMatrix<REAL>> fFields;       // Vetor de matrizes dos campos estocásticos
    TPZVec<TPZFMatrix<REAL>> fHFields;      // Vetor de matrizes para campos H
    TPZVec<TPZFMatrix<REAL>> fFieldSamples; // Amostras dos campos estocásticos

    TPZStack<TPZVec<TPZFMatrix<REAL>>> fFieldSamplesSubSetFN;

    TPZVec<REAL> fMeanvec;    // Vetor dos valores médios dos campos
    TPZVec<REAL> fCovvec;     // Vetor dos coeficientes de variação dos campos
    int fNSamples;            // Número de amostras geradas para o campo estocástico

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
