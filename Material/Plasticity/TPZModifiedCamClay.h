// TPZModifiedCamClay.h
//
// Modelo Cam-Clay modificado com integração implícita (return mapping) e operador
// tangente consistente, operando diretamente em Voigt 3D.
//
// Porte para o NeoPZ do módulo Python modified_cam_clay.py.
//
// Referência do algoritmo: de Souza Neto, Perić & Owen (2008), "Computational Methods
// for Plasticity", Seção 10.1, eqs. (10.1)-(10.22): sistema reduzido nas incógnitas
// {Δγ, α_{n+1}} resolvido por Newton-Raphson com a jacobiana (10.20).
//
// Extensões em relação ao texto do livro:
//   * endurecimento de estado crítico  p_c(α) = p_c0 exp(v0 α/(λ-κ))  (RS2, eq. 8.8, v constante);
//   * elasticidade "linear" (livro) ou "dependente da pressão":
//       p = p_ini exp(-(v0/κ) ε_v^e),  i.e.  K = -v0 p/κ  (RS2, eq. 8.7),
//     com G constante, ν constante (forma secante) ou ν constante hipoelástico (FLAC3D);
//   * tensão inicial σ0;
//   * operador tangente consistente 6x6.
//
// Convenções (as mesmas de TPZMatElastoPlastic / TPZElasticResponse neste ramo):
//   * tração positiva (p < 0 em compressão);
//   * Voigt {xx, xy, xz, yy, yz, zz} (a ordem de armazenamento de TPZTensor);
//   * deformações com distorção de engenharia (γ = 2ε) nas posições de cisalhamento;
//     tensões com as componentes tensoriais;
//   * α = -ε_v^p (deformação volumétrica plástica, compressão positiva), guardada em
//     TPZPlasticState::m_hardening;
//   * ε^e medida a partir do estado inicial (σ = σ0 quando ε^e = 0).
//
#ifndef TPZMODIFIEDCAMCLAY_H
#define TPZMODIFIEDCAMCLAY_H

#include <stdexcept>
#include <string>

#include "TPZPlasticBase.h"
#include "TPZPlasticCriterion.h"
#include "TPZPlasticState.h"
#include "TPZElasticResponse.h"
#include "TPZTensor.h"
#include "pzfmatrix.h"
#include "pzvec.h"

/// Superfície de escoamento e lei de endurecimento do Cam-Clay modificado
/// (eqs. 10.1-10.2 e 10.9 de de Souza Neto et al.).
///
///   f(p, q, a) = (p - p_t + a)²/b² + (q/M)² - a²,   b = 1 se p >= p_t - a, senão b = β
///   a(α) = (p_c(α) + p_t)/(1 + β),   p_c(α) = p_c0 exp(v0 α/(λ-κ))
class TPZYCModifiedCamClay : public TPZPlasticCriterion {
public:

    enum { NYield = 1 };

    TPZYCModifiedCamClay() = default;

    TPZYCModifiedCamClay(const TPZYCModifiedCamClay &other) = default;

    TPZYCModifiedCamClay &operator=(const TPZYCModifiedCamClay &other) = default;

    /// Define os parâmetros da superfície (pc0 > 0, pt >= 0: valores positivos)
    void SetUp(REAL M, REAL lambda, REAL kappa, REAL v0, REAL pc0, REAL pt = 0., REAL beta = 1.);

    REAL M() const { return fM; }
    REAL Lambda() const { return fLambda; }
    REAL Kappa() const { return fKappa; }
    REAL V0() const { return fV0; }
    REAL Pc0() const { return fPc0; }
    REAL Pt() const { return fPt; }
    REAL Beta() const { return fBeta; }

    /// p_c(α)
    REAL Pc(REAL alpha) const;

    /// a(α) (eq. 10.9) e H = da/dα (eq. 10.21)
    void Hardening(REAL alpha, REAL &a, REAL &H) const;

    /// fator b do ramo da elipse: 1 no lado "seco" (p >= p_t - a), β no lado "úmido"
    REAL BranchB(REAL p, REAL a) const { return (p >= fPt - a) ? 1. : fBeta; }

    /// f(p, q, a)  (p < 0 em compressão)
    REAL YieldValue(REAL p, REAL q, REAL a) const;

    // ----------------------------------------------------- TPZPlasticCriterion
    /// sigma: tensões principais; kprev: α. yield[0] = f(p, q, a(α))
    void YieldFunction(const TPZVec<STATE> &sigma, STATE kprev, TPZVec<STATE> &yield) const override;

    int GetNYield() const override { return NYield; }

    void SetLocalMatState(TPZPlasticState<REAL> &state) override {}

    TPZPlasticState<REAL> GetLocalMatState() override { return TPZPlasticState<REAL>(); }

    void ChangeLocalMatParameters(TPZPlasticState<REAL> &state, REAL factor) override {}

    int ClassId() const override;

    void Read(TPZStream &buf, void *context) override;

    void Write(TPZStream &buf, int withclassid) const override;

    void Print(std::ostream &out) const override;

private:
    REAL fM = 1.2;
    REAL fLambda = 0.077;
    REAL fKappa = 0.0066;
    REAL fV0 = 1.7;
    REAL fPc0 = 200.;
    REAL fPt = 0.;
    REAL fBeta = 1.;
};


/// Modelo Cam-Clay modificado: elasticidade + return mapping + tangente consistente.
///
/// Pode ser usado como argumento de template de TPZMatElastoPlastic / TPZMatElastoPlastic2D.
/// O estado (ε^p, α = m_hardening, ε_n = m_eps_t) é guardado em TPZPlasticState.
///
/// No modo hipoelástico (EHypoNu) a tensão do passo anterior σ_n é necessária: ela é lida do
/// argumento "sigma" de ApplyStrainComputeSigma na entrada (TPZMatElastoPlastic passa a
/// tensão guardada na memória do ponto de integração). Um tensor nulo é interpretado como
/// "estado inicial" e substituído por σ0 (o p = 0 é degenerado neste modelo).
class TPZModifiedCamClay : public TPZPlasticBase {
public:

    /// Lei elástica volumétrica
    enum EElasticity {
        ELinear = 0,            ///< K0 = v0 p0/κ constante (livro)
        EPressureDependent = 1  ///< p = p_ini exp(-(v0/κ) ε_v^e), K = -v0 p/κ (RS2, eq. 8.7)
    };

    /// Lei elástica desviadora
    enum EShear {
        EConstantG = 0,  ///< G constante
        EConstantNu = 1, ///< G = 3(1-2ν)/(2(1+ν)) K, forma secante s = 2G(p) e^e
        EHypoNu = 2      ///< hipoelástico: ds = 2 G_n de^e com G_n = 3(1-2ν)/(2(1+ν)) K(p_n) (FLAC3D)
    };

    /// Falha de convergência do return mapping local
    class ReturnMappingError : public std::runtime_error {
    public:
        using std::runtime_error::runtime_error;
    };

    /// Resultado do return mapping (equivalente ao dict retornado pelo Python)
    struct TResult {
        TPZTensor<REAL> stress;          ///< σ_{n+1}
        TPZFNMatrix<36, REAL> Dep;       ///< tangente consistente dσ/dε (6x6, Voigt, γ de engenharia)
        TPZTensor<REAL> elastic_strain;  ///< ε^e_{n+1} (engenharia)
        TPZTensor<REAL> plastic_strain;  ///< ε^p_{n+1} (engenharia)
        REAL alpha = 0.;                 ///< α_{n+1}
        REAL dgamma = 0.;                ///< Δγ
        bool plastic = false;
        int iterations = 0;
        REAL b = 1.;                     ///< ramo da elipse
        REAL p = 0., q = 0., pc = 0.;
        TResult() : Dep(6, 6, 0.) {}
    };

    TPZModifiedCamClay();

    TPZModifiedCamClay(const TPZModifiedCamClay &other);

    TPZModifiedCamClay &operator=(const TPZModifiedCamClay &other);

    ~TPZModifiedCamClay() override;

    /// Monta os parâmetros (equivalente a mcc_parameters do Python).
    ///
    /// p0, pc0: valores positivos (compressão), como na Tabela 8.1 do RS2.
    /// v0 <= 0  ->  v0 = N - λ ln(pc0) + κ ln(pc0/p0)  (linha virgem + linha de descarregamento).
    /// A tensão inicial é isotrópica, σ0 = -p0 I; use SetInitialStress para outra.
    void SetUp(REAL M, REAL lambda, REAL kappa, REAL N, REAL v0, REAL pc0, REAL p0,
               REAL pt = 0., REAL beta = 1., EElasticity elasticity = ELinear,
               EShear shear = EConstantNu, REAL G = 20000., REAL nu = 0.3);

    /// Redefine σ0 (tração positiva) e recalcula as grandezas derivadas (p_ini, K0, G0, e0)
    void SetInitialStress(const TPZTensor<REAL> &sigma0);

    /// Substitui K0 = -v0 p_ini/κ por um módulo volumétrico dado (só para ELinear, p.ex. uma
    /// variante com E e ν fixos). K0 <= 0 restaura o valor calculado.
    void SetLinearBulkModulus(REAL K0);

    /// Tolerância e número máximo de iterações do Newton local
    void SetLocalTolerance(REAL tol, int maxit = 50) { fTol = tol; fMaxIt = maxit; }

    /// Relação entre M e o ângulo de atrito de estado crítico, usada pelos parâmetros por ponto
    enum EStrengthMapping {
        ETriaxialCompression = 0, ///< M = 6 sinφ/(3 - sinφ) (compressão triaxial)
        EPlaneStrain = 1          ///< M = √3 sinφ (estado crítico em deformação plana com fluxo associado, θ = 0)
    };
    static REAL MFromFriction(REAL phi, EStrengthMapping mapping);
    void SetStrengthMapping(EStrengthMapping mapping) { fStrengthMapping = mapping; }
    EStrengthMapping StrengthMapping() const { return fStrengthMapping; }

    /// Fator de redução de resistência F (método de redução de resistência, como TPZPlasticStepPV):
    /// só atua nos parâmetros por ponto (φ_r = atan(tan φ / F), c_r = c/F)
    void SetStrengthReductionFactor(REAL F) { fReductionFactor = F; }
    REAL StrengthReductionFactor() const { return fReductionFactor; }

    /// Parâmetros por ponto de integração lidos de TPZPlasticState::fmatprop (o mesmo mecanismo do
    /// TPZYCMohrCoulombPV, usado nos campos aleatórios): quando fmatprop.size() >= 3 e fmatprop[2] > 0,
    ///   fmatprop[0] = c, [1] = φ (rad), [2] = p_c0 (> 0), [3..8] = σ0 (opcional; Voigt, tração positiva),
    /// e o modelo usa M = M(φ_r) e p_t = c cot φ = c_r cot φ_r (invariante pela redução de resistência), mantendo
    /// λ, κ, v0, β e a elasticidade. Chamado no início de ApplyStrainComputeSigma.
    void ApplyLocalProperties();

    // ------------------------------------------------------------- acesso
    const TPZYCModifiedCamClay &YC() const { return fYC; }
    REAL M() const { return fYC.M(); }
    REAL Lambda() const { return fYC.Lambda(); }
    REAL Kappa() const { return fYC.Kappa(); }
    REAL NIso() const { return fNIso; }
    REAL V0() const { return fYC.V0(); }
    REAL Pc0() const { return fYC.Pc0(); }
    REAL P0() const { return fP0; }
    REAL Pt() const { return fYC.Pt(); }
    REAL Beta() const { return fYC.Beta(); }
    REAL ShearModulusParameter() const { return fG; }
    REAL Nu() const { return fNu; }
    EElasticity Elasticity() const { return fElasticity; }
    EShear Shear() const { return fShear; }
    const TPZTensor<REAL> &InitialStress() const { return fSigma0; }
    REAL PIni() const { return fPIni; }  ///< p inicial (negativo em compressão)
    REAL K0() const { return fK0; }      ///< K no estado inicial
    REAL G0() const { return fG0; }      ///< G no estado inicial

    // ------------------------------------------------------------- lei constitutiva
    /// {p, q} de uma tensão
    static void Invariants(const TPZTensor<REAL> &sigma, REAL &p, REAL &q);

    /// p(ε_v^e), K = dp/dε_v^e e dK/dε_v^e
    void Pressure(REAL eev, REAL &p, REAL &K, REAL &dK) const;

    /// G(ε_v^e) e dG/dε_v^e (formas secantes; não usado no modo hipoelástico)
    void ShearModulus(REAL eev, REAL &G, REAL &dG) const;

    /// p_c(α)
    REAL Pc(REAL alpha) const { return fYC.Pc(alpha); }

    /// Return mapping implícito.
    ///
    /// eps: ε_{n+1} total; epsp_n: ε^p_n; alpha_n: α_n (Voigt, engenharia).
    /// sig_n, eps_n: tensão e deformação total do passo anterior (só para EHypoNu;
    /// nullptr -> σ0 e 0). Lança ReturnMappingError se o Newton local não convergir.
    void ReturnMapping(const TPZTensor<REAL> &eps, const TPZTensor<REAL> &epsp_n, REAL alpha_n,
                       TResult &res, const TPZTensor<REAL> *sig_n = nullptr,
                       const TPZTensor<REAL> *eps_n = nullptr) const;

    // ------------------------------------------------------------- TPZPlasticBase
    int ClassId() const override;
    void Write(TPZStream &buf, int withclassid) const override;
    void Read(TPZStream &buf, void *context) override;
    void Print(std::ostream &out) const override;
    const char *Name() const override { return "TPZModifiedCamClay"; }

    void ApplyStrain(const TPZTensor<REAL> &epsTotal) override;

    /// σ e (opcionalmente) a tangente consistente para ε_total; atualiza o estado.
    /// No modo EHypoNu, "sigma" deve conter σ_n na entrada (ver descrição da classe).
    void ApplyStrainComputeSigma(const TPZTensor<REAL> &epsTotal, TPZTensor<REAL> &sigma,
                                 TPZFMatrix<REAL> *tangent = nullptr) override;

    void ApplyStrainComputeDep(const TPZTensor<REAL> &epsTotal, TPZTensor<REAL> &sigma,
                               TPZFMatrix<REAL> &Dep) override;

    /// Inversa: deformação total que produz a tensão dada (Newton com a tangente consistente, a
    /// partir de ε_n). Com amolecimento (lado seco) a inversa pode não ser única: o Newton pode
    /// convergir para a solução elástica se a tensão estiver dentro da superfície do passo anterior.
    void ApplyLoad(const TPZTensor<REAL> &sigma, TPZTensor<REAL> &epsTotal) override;

    void SetState(const TPZPlasticState<REAL> &state) override { fN = state; }
    TPZPlasticState<REAL> GetState() const override { return fN; }

    TPZPlasticCriterion &GetYC() override { return fYC; }

    /// phi[0] = f(p, q, a(α_n))/a² com σ obtida da deformação ELÁSTICA dada (forma secante)
    void Phi(const TPZTensor<REAL> &epsElastic, TPZVec<REAL> &phi) const override;

    int IntegrationSteps() const override { return fLastIterations; }

    /// A elasticidade é a do próprio modelo: ER é guardada mas não é usada no cálculo
    void SetElasticResponse(TPZElasticResponse &ER) override { fER = ER; }

    /// Resposta elástica linear equivalente ao estado inicial (K0, G0)
    TPZElasticResponse GetElasticResponse() const override;

    /// Estado atual
    TPZPlasticState<REAL> fN;

private:
    void UpdateDerived();

    TPZYCModifiedCamClay fYC;

    REAL fNIso = 1.788;  ///< N: volume específico da NCL em p = 1
    REAL fP0 = 200.;     ///< p0 (positivo) usado em SetUp
    REAL fG = 20000.;    ///< G (modo EConstantG)
    REAL fNu = 0.3;      ///< ν (modos EConstantNu e EHypoNu)
    EElasticity fElasticity = ELinear;
    EShear fShear = EConstantNu;
    TPZTensor<REAL> fSigma0;

    REAL fK0User = 0.;   ///< K0 imposto (ELinear); <= 0 -> calculado

    // grandezas derivadas
    REAL fPIni = -200.;
    REAL fK0 = 0.;
    REAL fGFac = 0.;
    REAL fG0 = 0.;
    TPZTensor<REAL> fE0; ///< desviador inicial / 2G0 (componentes tensoriais)

    REAL fTol = 1.e-11;
    int fMaxIt = 50;

    EStrengthMapping fStrengthMapping = ETriaxialCompression;
    REAL fReductionFactor = 1.;
    int fLastIterations = 0;

    TPZElasticResponse fER;
};

#endif // TPZMODIFIEDCAMCLAY_H
