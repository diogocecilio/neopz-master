// SlopeStability.h
//
// Análise elastoplástica do talude (deformação plana, tensões efetivas) com a estrutura nativa do NeoPZ:
//
//   * TPZCompMesh H1 vetorial com memória nos pontos de integração (SetAllCreateFunctionsContinuousWithMem);
//   * material TPZMatElastoPlastic2D<T> (T = TPZPlasticStepVoigt<TPZYCMohrCoulombPV2> ou TPZModifiedCamClay) com
//     os parâmetros de resistência de cada ponto em TPZPlasticState::fmatprop (c, φ, ... — campos aleatórios) e
//     o fator de redução de resistência nativo (SetStrengthReductionFactor);
//   * forças de corpo γ' (SetBodyForce) e forças de percolação -grad u, guardadas na memória de cada ponto
//     (TPZElastoPlasticMem::fdPorePressure = grad u, fPorePressure = u) e somadas por
//     TPZMatElastoPlastic2DSeepage, que só acrescenta esse termo ao Contribute nativo;
//   * TPZLinearAnalysis (renumeração de banda, skyline) com o laço de Newton e a aceitação do estado na memória
//     (SetUpdateMem + AssembleResidual), como em TPZElastoPlasticAnalysis::AcceptSolution;
//   * TPZPostProcAnalysis para a saída VTK das variáveis nos pontos de integração.
//
// Medidas de estabilidade:
//   Γ  (fator de estabilidade do artigo, eq. 54): multiplicador λ das cargas λ(γ' e_y - grad u) no colapso,
//       obtido por carregamento incremental com cortes de passo a partir do estado nulo (Mohr-Coulomb);
//   FS (fator de segurança, eq. 58): redução c/F, tan φ/F com as cargas reais (método de redução de resistência).
// O colapso é a não convergência do Newton (critério usual de Griffiths & Lane), com o passo dividido até a
// tolerância relativa dada.
//
#ifndef SLOPESTABILITY_H
#define SLOPESTABILITY_H

#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "Plasticity/TPZElasticResponse.h"
#include "Plasticity/TPZElastoPlasticMem.h"
#include "Plasticity/TPZMatElastoPlastic2D.h"
#include "Plasticity/TPZModifiedCamClay.h"
#include "Plasticity/TPZPlasticStepVoigt.h"
#include "Plasticity/TPZYCMohrCoulombPV2.h"
#include "SlopeGeometry.h"
#include "TPZLinearAnalysis.h"
#include "pzcmesh.h"

/// Mohr-Coulomb associado com retorno no espaço das tensões principais e tangente consistente em Voigt
/// (distorções de engenharia). O TPZPlasticStepPV<TPZYCMohrCoulombPV> não é usado: o seu operador tangente não
/// é consistente com a convenção de engenharia do TPZElasticResponse (o Newton estagna).
using TMohrCoulomb = TPZPlasticStepVoigt<TPZYCMohrCoulombPV2, TPZElasticResponse>;

/// TPZMatElastoPlastic2D + forças de percolação f_s = -λ_s grad u lidas da memória do ponto de integração.
/// Variáveis de pós-processamento extras: "Cohesion", "FrictionAngle" (graus), "ExcessPorePressure",
/// "SeepageForce", "PlasticStrainNorm".
template <class T>
class TPZMatElastoPlastic2DSeepage : public TPZMatElastoPlastic2D<T, TPZElastoPlasticMem> {
    using TBase = TPZMatElastoPlastic2D<T, TPZElastoPlasticMem>;

public:
    // índices < 100 (TPZCompEl::Solution trata var >= 100 como solução por elemento) e acima dos do TPZMatElastoPlastic
    enum { ECohesion = 80, EFriction, EExcessPorePressure, ESeepageForce, EPlasticStrainNorm };

    explicit TPZMatElastoPlastic2DSeepage(int id) : TBase(id, 1) {}
    TPZMaterial *NewMaterial() const override { return new TPZMatElastoPlastic2DSeepage<T>(*this); }
    std::string Name() const override { return "TPZMatElastoPlastic2DSeepage"; }

    void SetSeepageFactor(REAL f) { fSeepageFactor = f; }
    REAL SeepageFactor() const { return fSeepageFactor; }

    void Contribute(const TPZMaterialDataT<STATE> &data, REAL weight, TPZFMatrix<STATE> &ek,
                    TPZFMatrix<STATE> &ef) override {
        TBase::Contribute(data, weight, ek, ef);
        AddSeepage(data, weight, ef);
    }
    void Contribute(const TPZMaterialDataT<STATE> &data, REAL weight, TPZFMatrix<STATE> &ef) override {
        TBase::Contribute(data, weight, ef);
        AddSeepage(data, weight, ef);
    }

    int VariableIndex(const std::string &name) const override;
    int NSolutionVariables(int var) const override;
    void Solution(const TPZMaterialDataT<STATE> &data, int var, TPZVec<REAL> &Solout) override;

private:
    void AddSeepage(const TPZMaterialDataT<STATE> &data, REAL weight, TPZFMatrix<STATE> &ef) {
        if (fSeepageFactor == 0.) return;
        const TPZElastoPlasticMem &mem = this->MemItem(data.intGlobPtIndex);
        if (mem.fdPorePressure.size() < 2) return;
        const REAL fx = -fSeepageFactor * mem.fdPorePressure[0], fy = -fSeepageFactor * mem.fdPorePressure[1];
        const TPZFMatrix<REAL> &phi = data.phi;
        for (int in = 0; in < phi.Rows(); in++) {
            ef(2 * in + 0, 0) += weight * fx * phi(in, 0);
            ef(2 * in + 1, 0) += weight * fy * phi(in, 0);
        }
    }
    REAL fSeepageFactor = 0.;
};

/// Propriedades do solo (valores médios; os campos aleatórios entram por ponto)
struct TSoil {
    REAL gamma = 20.;        ///< peso específico saturado (ou natural, talude seco), kN/m³
    REAL gammaW = 10.;       ///< peso específico da água, kN/m³
    bool buoyant = true;     ///< força de corpo γ' = γ - γw (talude saturado, tensões efetivas); false: γ
    REAL E = 1.e5, nu = 0.3; ///< elasticidade (kPa)
    REAL c = 10., phiDeg = 30.;
    // Cam-Clay modificado (resistência de estado crítico equivalente: M(φ), p_t = c cot φ)
    REAL lambda = 0.10, kappa = 0.02, v0 = 2.0, OCR = 1.0;
    TPZModifiedCamClay::EStrengthMapping mapping = TPZModifiedCamClay::EPlaneStrain;
    REAL BodyForce() const { return buoyant ? gamma - gammaW : gamma; }
};

/// Modelos com as propriedades médias do solo (elasticidade, c, φ; Cam-Clay: M(φ), p_t = c cot φ, λ, κ, v0)
void SetupSoilModel(TMohrCoulomb &mc, const TSoil &soil);
void SetupSoilModel(TPZModifiedCamClay &mcc, const TSoil &soil);
/// Cam-Clay: com c = mp[0] e φ = mp[1] (rad), define p_c0 = OCR p_c(σ0) (mp[2]) e σ0 (mp[3..8]) a partir da
/// tensão inicial s (com p' >= 1 kPa de compressão; s é ajustada)
void CamClayInitialState(TPZVec<REAL> &mp, TPZTensor<REAL> &s, const TSoil &soil);

struct TSolverOptions {
    int porder = 2;
    REAL tol = 1.e-5;         ///< ||R_livre|| <= tol ||F_ref||
    int maxIter = 20;         ///< Newton: iterações por passo (não convergência = colapso)
    bool stagnation = true;   ///< interrompe o Newton sem redução de 1/2 do resíduo em 3 iterações
    REAL relTol = 5.e-3;      ///< precisão relativa do fator (Γ ou FS)
    REAL step0 = 0.25;        ///< passo inicial do fator
    REAL maxFactor = 20.;
    int verbose = 0;
};

struct TFactorResult {
    REAL factor = 0.;         ///< último valor convergido (limite inferior do colapso)
    REAL upper = 0.;          ///< primeiro valor sem convergência (limite superior)
    int steps = 0, iterations = 0, cuts = 0;
    bool bracketed = false;   ///< colapso encontrado (upper válido)
    std::string status;
};

/// Malha e análise elastoplástica do talude para o modelo T
template <class T>
class TSlopeFEM {
public:
    struct TPoint {
        int64_t gel;                 ///< índice do elemento geométrico
        TPZManVector<REAL, 3> qsi;   ///< coordenadas paramétricas
        TPZManVector<REAL, 3> x;     ///< coordenadas
    };

    TSlopeFEM(TPZGeoMesh *gmesh, const TSlopeGeometry &geo, const TSoil &soil, const TSolverOptions &opt);
    ~TSlopeFEM();

    int64_t NPoints() const { return (int64_t)fPoints.size(); }
    /// Pontos de integração na ordem dos índices de memória
    const std::vector<TPoint> &Points() const { return fPoints; }
    int64_t NEquations() const { return fCMesh->NEquations(); }

    /// Resistência por ponto (índice de memória): c (kPa) e φ (rad). Para o Cam-Clay também p_c0 e σ0
    /// (SetInitialStress) entram em fmatprop.
    void SetStrength(const std::vector<REAL> &c, const std::vector<REAL> &phi);
    void SetUniformStrength(REAL c, REAL phi);
    /// Excesso de poropressão e gradiente por ponto (vazios: sem percolação)
    void SetSeepage(const std::vector<REAL> &u, const std::vector<TPZManVector<REAL, 2>> &gradu);
    /// Tensão inicial por ponto (Cam-Clay): define σ0 e p_c0 = OCR p_c(σ0) com M e p_t locais
    void SetInitialStress(const std::vector<TPZTensor<REAL>> &sigma0);
    /// Tensões atuais (m_sigma) por ponto
    void Stresses(std::vector<TPZTensor<REAL>> &sigma) const;

    /// Zera deformações, deslocamentos e (MC) tensões; o Cam-Clay volta a σ0. Mantém fmatprop e a percolação.
    void ResetState();

    /// Γ: multiplicador das cargas λ(γ' e_y - grad u) no colapso, por carregamento incremental a partir do
    /// estado atual (λ0 = fator já aplicado)
    TFactorResult LoadFactor(REAL lambda0 = 0.);
    /// FS por redução de resistência: aplica as cargas reais (γ' e percolação; com F = F0) e aumenta F até o
    /// colapso
    TFactorResult StrengthReduction(REAL F0 = 0.5);

    /// Fixa cargas e fator de redução e resolve um passo a partir do estado aceito; aceita se convergir
    bool Solve(REAL lambdaGravity, REAL lambdaSeepage, REAL F, int &iterations);

    /// Indicador por elemento geométrico (índice): máximo, nos pontos de integração, de ||ε^p|| (increment =
    /// false) ou de ||Δε^p|| no último passo aceito (increment = true: o mecanismo ativo no colapso; o acúmulo de
    /// ε^p é dominado pela concentração de tensões no pé do talude)
    void PlasticIndicator(std::vector<REAL> &byGel, bool increment = true) const;
    /// Mecanismo de colapso por ponto de integração (índice de memória), no último passo aceito: ||Δε^p|| e
    /// ||Δu|| (deslocamento acumulado na memória). Vazios antes do primeiro passo aceito.
    void MechanismIndicators(std::vector<REAL> &depsp, std::vector<REAL> &du) const;

    /// Saída VTK (TPZPostProcAnalysis): tensões, deformação plástica, c, φ, u, forças de percolação
    void DefineVTK(const std::string &file);
    void WriteVTK(int step);

    TPZCompMesh *Mesh() { return fCMesh; }
    TPZMatElastoPlastic2DSeepage<T> *Material() { return fMat; }

private:
    void BuildPoints();
    void SetLoads(REAL lambdaGravity, REAL lambdaSeepage, REAL F);
    void Accept();
    void ComputeReference();
    TFactorResult Follow(const std::function<void(REAL)> &apply, REAL t0, REAL dt0, REAL tmax);

    TPZGeoMesh *fGMesh;
    TSlopeGeometry fGeo;
    TSoil fSoil;
    TSolverOptions fOpt;
    TPZCompMesh *fCMesh = nullptr;
    TPZMatElastoPlastic2DSeepage<T> *fMat = nullptr;
    std::unique_ptr<TPZLinearAnalysis> fAn;
    std::vector<TPoint> fPoints;
    std::vector<bool> fFree;
    REAL fFRef = 1.;
    REAL fLambdaG = 0., fLambdaS = 0., fF = 1.;                   ///< estado aceito
    REAL fLambdaGTarget = 0., fLambdaSTarget = 0., fFTarget = 1.; ///< alvo do passo em Follow
    std::vector<TPZTensor<REAL>> fSigma0;
    std::vector<TPZTensor<REAL>> fEpsPPrev;  ///< ε^p no início do último passo aceito
    std::vector<TPZManVector<REAL, 3>> fUPrev;  ///< deslocamento no início do último passo aceito
    class TPZPostProcAnalysis *fPost = nullptr;
};

#endif
