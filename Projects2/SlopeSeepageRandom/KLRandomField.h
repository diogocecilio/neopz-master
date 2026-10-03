// KLRandomField.h
//
// Campos aleatórios gaussianos (e lognormais) discretizados pela expansão de Karhunen–Loève, com as rotinas de
// campo estocástico do NeoPZ (Projects2/GeoMecRandFields*, eigensolverkl):
//
//   * malha H1 contínua com TPZMatKLKernel (kernel C(x,y) = exp(-|x1-y1|/Lx - |x2-y2|/Ly), eq. 64 do artigo);
//   * C e a massa B montadas por pzdoublestrmatriz (pares de elementos, Galerkin);
//   * problema generalizado simétrico-definido C v = λ B v resolvido por TPZLapackEigenSolver (dsbgv, autovetores
//     reais B-ortonormais);
//   * Φ = [√λ_k φ_k] (φ_k normalizados em L2 por TPZCompMesh::Integrate("SolutionSquared"), como em
//     BuildPhiSqrtLambda);
//   * realização: H(x) = Σ_k √λ_k ξ_k φ_k(x), ξ_k ~ N(0,1) (eq. 61), avaliada nos pontos de interesse com
//     TPZInterpolationSpace::ComputeSolution (malha KL carregada com Φ ξ).
//
// O erro de truncamento médio em variância (Allaix & Carbone) é ε_M = 1 - Σ_{k<=M} λ_k / |Ω|; a variância pontual
// truncada v(x) = Σ_k λ_k φ_k(x)² < 1 é, por opção, compensada (H/√v), para que o CoV alvo seja respeitado.
//
// Pontos de interesse ("conjuntos de alvos") são registrados uma vez: para cada ponto guarda-se o elemento KL e
// as coordenadas paramétricas (FindElement com índice inicial válido + busca exaustiva de reserva), e v(x).
//
#ifndef KLRANDOMFIELD_H
#define KLRANDOMFIELD_H

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include "pzcmesh.h"
#include "pzgmesh.h"
#include "pzmanvector.h"

class TPZInterpolationSpace;

class TPZKLRandomField {
public:
    struct TOptions {
        REAL Lx = 20., Ly = 2.;          ///< distâncias de autocorrelação
        int porder = 2;                  ///< ordem da malha KL (2 em quadriláteros = 9 nós, como no artigo)
        int nModes = -1;                 ///< M (<= 0: todos os modos)
        REAL targetVarianceError = -1.;  ///< se > 0, usa o menor M com ε_M <= alvo
        bool normalizeVariance = true;   ///< divide H(x) por √v(x)
        /// arquivo binário com λ e Φ (vazio: sem cache). Gravação atômica (tmp + rename); o cabeçalho guarda
        /// versão, Lx, Ly, ordem, neq, número de elementos e de nós e assinatura da malha KL, |Ω| e o número de
        /// modos gravados; arquivo truncado, corrompido (FNV-1a) ou incompatível é recalculado
        std::string cacheFile;
    };

    /// gmesh: malha geométrica do domínio (a classe passa a ser dona dela)
    TPZKLRandomField(TPZGeoMesh *gmesh, const TOptions &opt);
    ~TPZKLRandomField();

    /// Monta e resolve o autoproblema (ou lê do cache) e calcula Φ
    void Compute();

    int NModes() const { return fM; }
    int NEquations() const;
    REAL DomainArea() const { return fArea; }
    const std::vector<REAL> &Eigenvalues() const { return fLambda; }
    /// ε_M = 1 - Σ_{k<M} λ_k / |Ω|
    REAL VarianceError(int M) const;

    /// Registra um conjunto de pontos; devolve o identificador do conjunto
    int AddTargetSet(const std::vector<TPZManVector<REAL, 3>> &points);
    int NTargets(int set) const { return (int)fSets[set].pts.size(); }
    /// v(x) nos pontos do conjunto (variância da expansão truncada)
    const std::vector<REAL> &PointVariance(int set) const { return fSets[set].var; }

    /// Variáveis ξ (M) da realização "sample" do campo "field": gerador mt19937_64 semeado com
    /// seed_seq{seed, sample, field}, reprodutível e independente da ordem em que as amostras são calculadas
    void Xi(uint64_t seed, int64_t sample, int field, std::vector<REAL> &xi) const;

    /// Campo gaussiano padrão (média 0, variância 1 após a normalização) nos pontos dos conjuntos dados, para
    /// vários campos de uma vez: values[f][set][i]. xi[f] tem M componentes.
    void Evaluate(const std::vector<std::vector<REAL>> &xi, const std::vector<int> &sets,
                  std::vector<std::vector<std::vector<REAL>>> &values);

    /// Transformação lognormal (eqs. 65-66): X = exp(μ' + σ' H), μ' = ln(μ/√(1+CoV²)), σ' = √ln(1+CoV²)
    static REAL Lognormal(REAL h, REAL mean, REAL cov);

    /// Escreve a realização H = Φ ξ (sem a normalização por √v) em VTK por DefineGraphMesh/PostProcess ("KLField")
    void WriteVTK(const std::vector<REAL> &xi, const std::string &file, int step = 0);

    TPZCompMesh *Mesh() { return fCMesh; }

private:
    struct TTargetSet {
        std::vector<TPZManVector<REAL, 3>> pts;
        std::vector<TPZInterpolationSpace *> cel;
        std::vector<TPZManVector<REAL, 3>> qsi;
        std::vector<REAL> var;
    };
    void BuildMesh();
    bool ReadCache();
    void WriteCache() const;
    /// M a partir de λ e das opções (nModes ou targetVarianceError)
    int ChooseModes() const;
    void LocatePoint(const TPZManVector<REAL, 3> &x, TPZInterpolationSpace *&cel, TPZManVector<REAL, 3> &qsi,
                     int64_t &start) const;

    TOptions fOpt;
    std::unique_ptr<TPZGeoMesh> fGMesh;
    TPZCompMesh *fCMesh = nullptr;
    REAL fArea = 0.;
    int fM = 0;
    std::vector<REAL> fLambda;   ///< todos os autovalores calculados (decrescentes)
    TPZFMatrix<STATE> fPhi;      ///< neq x M, colunas √λ_k φ_k
    std::vector<TTargetSet> fSets;
};

#endif
