#ifndef ELASTOPLASTICANALYSIS_H
#define ELASTOPLASTICANALYSIS_H

#include "pznonlinanalysis.h"
#include "pzcompel.h"
#include "TPZGeoElement.h"
#include "pzfmatrix.h"
#include "pzvec.h"
#include "pzpostprocanalysis.h"
#include <iostream>
#include "tpzgeoelrefpattern.h"
#include "Plasticity/pzelastoplasticanalysis.h"
#include "Plasticity/TPZElasticResponse.h"
#include "Plasticity/TPZYCMohrCoulombPV.h"
#include "Plasticity/TPZMatElastoPlastic2D.h"
#include "Plasticity/TPZMatElastoPlastic.h"
#include "Plasticity/TPZPlasticStepPV.h"
#include <deque>
/**
 * @brief Análise não-linear elastoplástica (driver) para malhas do NeoPZ.
 *
 * Responsável por montar resíduo/tangente, executar o processo iterativo
 * (Newton/Krylov), globalizar com *line search* e aceitar a solução,
 * incluindo mecanismos de atualização de memória plástica e *post-process*.
 *
 * Convenções nesta classe:
 * - Métrica de mérito usada nos *line searches*: \f$\phi(\alpha)=\tfrac12\|R(W_n+\alpha\Delta W)\|_2^2\f$
 *   (salvo indicado).
 * - `AcceptSolution()` propaga o incremento para a solução acumulada e
 *   sinaliza materiais com memória para atualizar estados internos.
 */
class TPZElastoPlasticAnalysis : public TPZLinearAnalysis {

public:

	enum class ELineSearch {
		None = 0,
		Armijo,
		QuadraticArmijo,
		Dicotomic,
		GoldenSection,
		StrongWolfe,
		NonmonotoneArmijo
	};

	/** @name Construção/Destruição */
	///@{
	/** @brief Construtor principal.
	 *  @param mesh Malha computacional.
	 *  @param out  Saída de log/diagnóstico.
	 *  @param lsearch Tipo de Busca Unidimensional.
	 */

	// --- NOVO: tipos de line search ---

	/** Construtor principal. Agora aceita o tipo de busca (default: Dicotomic) */
	TPZElastoPlasticAnalysis(TPZCompMesh *mesh,
							 std::ostream &out,
						  ELineSearch lsearch = ELineSearch::Dicotomic);

	/** @brief Construtor default (sem malha). */
	TPZElastoPlasticAnalysis();

	/** @brief Destrutor. */
	virtual ~TPZElastoPlasticAnalysis();
	///@}

	/** @name Processo iterativo (nível alto) */
	///@{
	/**
	 * @brief Tenta localizar a raiz do resíduo \f$R(u)=0\f$.
	 * @param iters [out] Iterações efetivamente executadas.
	 * @return true se convergiu segundo critérios internos.
	 */
	bool FindRoot(int & iters,REAL &resu,REAL &resf);

	bool NewtonRaphson(bool verbose);
	bool NewtonRaphson(REAL tol,TPZStack<STATE> &outresF,TPZStack<STATE> &outresU);

	/**
	 * @brief Processo iterativo não-linear (tipo Newton) com *line search* opcional.
	 * @param out       Stream para log.
	 * @param tol       Tolerância alvo (ex.: norma relativa do resíduo).
	 * @param numiter   Máximo de iterações.
	 * @param linesearch Se true, habilita busca linear.
	 * @param checkconv Se true, checa convergência intermediária.
	 * @param iters     [out] Iterações realizadas.
	 * @return true se convergiu.
	 */
	bool IterativeProcess ( std::ostream &out, REAL tol, int numiter,
							bool linesearch, bool checkconv, int &iters );
	///@}

	/** @name Line searchs (globalização) */
	///@{
	/**
	 * @brief Line search de Armijo (backtracking).
	 *
	 * Procura \f$\alpha\in(0,1]\f$ tal que
	 * \f$\phi(\alpha)\le\phi(0)+c_1\alpha\phi'(0)\f$
	 * (decréscimo suficiente). Reduz o passo multiplicativamente
	 * até satisfazer Armijo ou esgotar tentativas.
	 *
	 * @param Wn     Estado atual.
	 * @param DeltaW Direção de busca (incremento de Newton/Krylov).
	 * @param NextW  [out] Estado atualizado \f$W_{n+1}=W_n+\alpha\Delta W\f$.
	 * @param tol    Passo mínimo permitido (limite inferior para \f$\alpha\f$).
	 * @param niter  Máximo de reduções.
	 * @return \f$\alpha\f$ aceito (0 se não houve avanço).
	 */
	REAL ArmijoLineSearch ( const TPZFMatrix<STATE> &Wn,
							TPZFMatrix<STATE> DeltaW,
						 TPZFMatrix<STATE> &NextW,
						 REAL tol, int niter );

	/**
	 * @brief Armijo com interpolação quadrática protegida.
	 *
	 * Usa modelo quadrático local com \f$\phi(0),\phi'(0),\phi(\alpha)\f$
	 * para propor novo passo entre \([0.1\alpha,0.5\alpha]\) (safeguard),
	 * reduzindo avaliações e aceitando passos maiores quando possível.
	 *
	 * @see ArmijoLineSearch
	 */
	REAL QuadraticArmijoLineSearch(const TPZFMatrix<STATE>& Wn,
								   TPZFMatrix<STATE> DeltaW,
								   TPZFMatrix<STATE>& NextW,
								   REAL tol, int niter);

	/**
	 * @brief Line search por busca dicotômica (dichotomous search).
	 *
	 * Minimiza \f$\phi(\alpha)\f$ em \([0,1]\) comparando dois pontos
	 * simétricos ao redor do meio do intervalo e descartando a metade pior
	 * até atingir a tolerância em \f$\alpha\f$ ou o limite de iterações.
	 */
	REAL DicotomicLineSearch(const TPZFMatrix<STATE>& Wn,
							 TPZFMatrix<STATE> DeltaW,
							 TPZFMatrix<STATE>& NextW,
							 REAL tol, int niter);

	REAL GoldenSectionLineSearch(const TPZFMatrix<STATE>& Wn,
														   TPZFMatrix<STATE> DeltaW,
														   TPZFMatrix<STATE>& NextW,
														   REAL tol, int niter);

	REAL StrongWolfeLineSearch(const TPZFMatrix<STATE>& Wn,
							   const TPZFMatrix<STATE>& d,
							   TPZFMatrix<STATE>& NextW,
							   REAL c1=1e-4, REAL c2=0.9,
							   REAL a0=1.0, int max_eval=20);

	REAL NonmonotoneArmijoGLL(const TPZFMatrix<STATE>& Wn,
														const TPZFMatrix<STATE>& d,
														TPZFMatrix<STATE>& NextW,
														std::deque<REAL>& phi_hist, // mantém últimas M φ
														REAL c1=1e-4, REAL beta=0.5,
														int M=5, REAL a0=1.0, int max_red=20);
	void SetLineSearch(ELineSearch kind) { fLineSearch = kind; }
	ELineSearch GetLineSearch() const { return fLineSearch; }
	///@}

	/** @name Carga/aceitação de solução */
	///@{
	/**
	 * @brief Carrega um vetor de solução na malha (e multiphysics, se houver).
	 * @param loadsol Solução a carregar (não-const por API do NeoPZ).
	 */
	void  LoadSolution ( TPZFMatrix<STATE> & loadsol );

	/**
	 * @brief Executa passos de carregamento em progressão geométrica numa BC.
	 *
	 * Atualiza valores de `bc.Val1/Val2` do início ao fim em `nsteps`, usando
	 * razão `PGRatio`. Em cada passo, roda o processo iterativo com tolerância
	 * `tol` e no máx. `numiter` iterações e pode pós-processar resultados.
	 *
	 * @param out        Stream para log.
	 * @param tol        Tolerância por passo.
	 * @param numiter    Máximo de iterações por passo.
	 * @param BCId       ID da condição de contorno a ser escalonada.
	 * @param nsteps     Nº de passos.
	 * @param PGRatio    Razão da PG.
	 * @param val1Begin  Valor inicial de Val1 (matriz NxN conforme material).
	 * @param val1End    Valor final de Val1.
	 * @param val2Begin  Valor inicial de Val2.
	 * @param val2End    Valor final de Val2.
	 * @param ppAnalysis (opcional) objeto de pós-processamento.
	 * @param res        (opcional) resolução de saída do pós-processamento.
	 */
	virtual void ManageIterativeProcess(std::ostream &out, REAL tol, int numiter,
										int BCId, int nsteps, REAL PGRatio,
										TPZFMatrix<REAL> & val1Begin, TPZFMatrix<REAL> & val1End,
										TPZFMatrix<REAL> & val2Begin, TPZFMatrix<REAL> & val2End,
										TPZPostProcAnalysis * ppAnalysis = NULL, int res = 0);

	/**
	 * @brief Aceita a solução atual: atualiza memória plástica e solução acumulada.
	 *
	 * Sinaliza materiais com memória para atualizar estados internos durante a
	 * montagem, monta o rhs para efetivar a atualização e então desativa a flag.
	 *
	 * @param ResetOutputDisplacements Se !=0, redefine o acumulado; caso contrário, soma o incremento.
	 * @return Norma do incremento aceito (pode ser usada para diagnósticos).
	 */
	virtual REAL AcceptSolution(const int ResetOutputDisplacements = 0);

	/** @brief Carrega a solução corrente em todas as malhas relevantes. */
	virtual void LoadSolution();
	///@}

	/** @name Acesso/solvers utilitários */
	///@{
	/** @brief Retorna a solução acumulada (cumulativa). */
	TPZFMatrix<REAL> &CumulativeSolution() { return fCumSol; }

	/** @brief Define pré-condicionador genérico. */
	void SetPrecond(TPZMatrixSolver<REAL> &precond);

	/** @brief Configura BiCGStab com tolerância e iterações. */
	void SetBiCGStab(int numiter, REAL tol);

	/** @brief Configura BiCGStab com Jacobi como pré-condicionador. */
	void SetBiCGStab_Jacobi(int numiter, REAL tol);

	/** @brief Define solver direto LU. */
	void SetLU();

	/** @brief Transfere solução para uma malha de pós-processamento. */
	void TransferSolution(TPZPostProcAnalysis & ppanalysis);

	/** @brief Transfere solução para a malha de pós-processamento interna. */
	void TransferSolution();
	///@}

	/** @name Restrições de não-penetração (util para contato simples) */
	///@{
	/**
	 * @brief Adiciona material com restrição de não penetração numa direção.
	 * @param matid     ID do material de contorno/contato.
	 * @param direction 0 = x, 1 = y (convenção adotada).
	 */
	void AddNoPenetration(int matid, int direction)
	{
		fMaterialIds.insert(std::pair<int,int>(matid, direction));
	}

	/** @brief Constrói a estrutura de equações a zerar com base nos materiais marcados. */
	void IdentifyEquationsToZero();

	/** @brief Obtém o vetor de índices de equações ativas (não zeradas). */
	void GetActiveEquations(TPZVec<int64_t> &activeEquations);
	///@}

protected:
	/** @name Suporte interno */
	///@{
	/**
	 * @brief Liga/desliga atualização de memória plástica nos materiais com memória.
	 * @param update 1 = atualizar memória nas próximas montagens; 0 = montagem convencional.
	 */
	void SetUpdateMem(int update);

	/** @brief Atualiza a matriz de pré-condicionador (bloco diagonal, etc.). */
	void UpdatePrecond();
	///@}

public:
	/** @name Ganchos de análise (para extensões) */
	///@{
	/** @brief Checa critérios de convergência e reporta no log. */
	void CheckConv(std::ostream &out, REAL range);

	/** @brief Monta a matriz tangente (ou combinação linear) para um caso. */
	virtual void ComputeTangent(TPZFMatrix<REAL> &tangent, TPZVec<REAL> &coefs, int icase);

	/** @brief Nº de casos (para estratégias multi-caso). */
	virtual int NumCases();

	/** @brief Monta o resíduo para um dado caso. */
	virtual void Residual(TPZFMatrix<REAL> &residual, int icase);

	/** @brief Ajusta *create functions* para elementos com memória. */
	static void SetAllCreateFunctionsWithMem(TPZCompMesh *cmesh);
	///@}

	/** @name Multiphysics (opcional) */
	///@{
	/** @brief Indica se há configuração multiphysics associada. */
	bool IsMultiPhysicsConfiguration()
	{
		// fixado como falso no momento
		return 0;
	}

	/** @brief Define estrutura multiphysics (malha acoplada e vetor de malhas base). */
	void SetMultiPhysics(TPZCompMesh *mphysics, TPZVec<TPZCompMesh *> &meshvec)
	{
		fMultiPhysics = mphysics;
		fMeshVec = meshvec;
	}

	/** @brief Reseta ponteiros/estruturas multiphysics. */
	void ResetMultiPhysics()
	{
		fMultiPhysics = 0;
		fMeshVec.Resize(0);
	}
	///@}

	/** @name Tipos auxiliares/material padrão */
	///@{
	/// Alias do material elastoplástico 2D usado frequentemente.
	typedef TPZMatElastoPlastic2D<
	TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse>, TPZElastoPlasticMem> plasticmat;
	///@}

	/** @name Utilidades de refinamento dirigido por indicador */
	///@{
	/** @brief Escala vetor de carga (ou BCs relevantes) por um fator. */
	void LoadingRamp ( REAL factor );

	/** @brief Marca e divide elementos cujo indicador excede um limiar. */
	void DivideElementsAbove(REAL refineaboveval, std::set<long> &elindices);

	/** @brief Aumenta a ordem p de elementos cujo indicador excede um limiar. */
	void PRefineElementsAbove(REAL refineaboveval, int porder, std::set<long> &elindices);

	/** @brief Recalcula indicadores de deformação plástica por elemento. */
	void ComputeElementDeformation();
	///@}

public:
	/** @name Pós-processamento */
	///@{
	/** @brief Gera VTK pós-processado com variáveis plásticas/estruturais. */
	void PostPlasticity(std::string vtkd);

	/** @brief Cria malha de pós-processamento (supermesh) associada. */
	void CreatePostProcessingMesh (TPZPostProcAnalysis * PostProcess );

	/** @brief Seleciona variáveis escalares/vetoriais de pós-processamento. */
	void PostProcessVariables ( TPZStack<std::string> &scalNames, TPZStack<std::string> &vecNames );
	///@}
	void ToggleUpdateMem(bool on) { this->SetUpdateMem(on ? 1 : 0); }
protected:
	/** @name Dados internos */
	///@{
	/** @brief Vetor de solução acumulada (somatória dos incrementos aceitos). */
	TPZFMatrix<REAL> fCumSol;

	/** @brief Pré-condicionador (propriedade não-donas do ponteiro). */
	TPZMatrixSolver<REAL> * fPrecond = nullptr;

	/** @brief Índices de equações a zerar (p.ex., restrições de contato). */
	std::set<int64_t> fEquationstoZero;

	/** @brief Mapeia materiais com *no-penetration* para direção 0(x)/1(y). */
	std::multimap<int,int> fMaterialIds;

	/** @brief Malha multiphysics (se configurada). */
	TPZCompMesh *fMultiPhysics = nullptr;

	/** @brief Malhas base associadas à configuração multiphysics. */
	TPZManVector<TPZCompMesh *,2> fMeshVec;

	ELineSearch fLineSearch = ELineSearch::Dicotomic;

	std::deque<REAL> fPhiHistory;

	///@}

	/**
	 * @name Fábricas de elementos com memória
	 * @brief Funções *create* para compor elementos computacionais com memória.
	 * @note Mantidas aqui (não como *static* templates) para evitar necessidade
	 *       de parâmetros *dummy* em chamadas genéricas.
	 */
	///@{
	static TPZCompEl * CreateCubeElWithMem(  TPZGeoEl *gel, TPZCompMesh &mesh, int64_t &index);
	static TPZCompEl * CreateLinearElWithMem(TPZGeoEl *gel, TPZCompMesh &mesh, int64_t &index);
	static TPZCompEl * CreatePointElWithMem( TPZGeoEl *gel, TPZCompMesh &mesh, int64_t &index);
	static TPZCompEl * CreatePrismElWithMem( TPZGeoEl *gel, TPZCompMesh &mesh, int64_t &index);
	static TPZCompEl * CreatePyramElWithMem( TPZGeoEl *gel, TPZCompMesh &mesh, int64_t &index);
	static TPZCompEl * CreateQuadElWithMem(  TPZGeoEl *gel, TPZCompMesh &mesh, int64_t &index);
	static TPZCompEl * CreateTetraElWithMem( TPZGeoEl *gel, TPZCompMesh &mesh, int64_t &index);
	static TPZCompEl * CreateTriangElWithMem(TPZGeoEl *gel, TPZCompMesh &mesh, int64_t &index);
	///@}
};

#endif
