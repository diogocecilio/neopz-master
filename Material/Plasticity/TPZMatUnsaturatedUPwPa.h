/**
 * @file TPZMatUnsaturatedUPwPa.h
 * @brief Coupled three-phase material for unsaturated soils: deformable skeleton (u), water (p_w) and air (p_a),
 * elastoplastic effective (Bishop) stress, plane strain.
 */

#ifndef TPZMATUNSATURATEDUPWPA_H
#define TPZMATUNSATURATEDUPWPA_H

#include "TPZMatBase.h"
#include "TPZMatCombinedSpaces.h"
#include "TPZMatWithMem.h"
#include "TPZMaterialDataT.h"
#include "TPZBndCondT.h"
#include "TPZElastoPlasticMem.h"
#include "TPZUnsaturatedHydraulics.h"
#include "TPZTensor.h"
#include "pzfmatrix.h"
#include <functional>

class TPZCompMesh;

/**
 * @brief Memory of an integration point of TPZMatUnsaturatedUPwPa: the elastoplastic state (TPZElastoPlasticMem) and the
 * pressures of the last converged step, at the point and at the vertices of the element (lumped storage)
 */
class TPZUnsaturatedMem : public TPZElastoPlasticMem {
public:
    TPZUnsaturatedMem() : TPZElastoPlasticMem() {
        for (int i = 0; i < 3; ++i) fPwV[i] = fPaV[i] = 0.;
    }
    TPZUnsaturatedMem(const TPZUnsaturatedMem &cp) = default;
    TPZUnsaturatedMem &operator=(const TPZUnsaturatedMem &cp) = default;

    /** @brief p_w and p_a at the point (last converged step) */
    REAL fPw = 0., fPa = 0.;
    /** @brief water saturation at the point (last converged step, output) */
    REAL fS = 1.;
    /** @brief p_w and p_a at the three vertices of the (P1) pressure element, last converged step */
    REAL fPwV[3], fPaV[3];

    int ClassId() const override;
    void Write(TPZStream &buf, int withclassid) const override;
    void Read(TPZStream &buf, void *context) override;
    void Print(std::ostream &out = std::cout) const override;
    friend std::ostream &operator<<(std::ostream &out, const TPZUnsaturatedMem &m) {
        m.Print(out);
        return out;
    }
};

/**
 * @ingroup matcombinedspaces
 * @brief u - p_w - p_a material with a Mohr-Coulomb (or any TPZPlasticStepVoigt-like) effective stress.
 *
 * Formulation: [F] "Formulação acoplada u-p_w-p_a com Mohr-Coulomb associado"
 * (docs/biot_nao_saturado/formulacao/formulacao.tex, Fases 1 e 2); equation numbers below are those of [F]. [F]
 * derives from [N] the notes "Formulação poroelástica não saturada: esqueleto deformável, água e ar"
 * (docs/biot_nao_saturado/biot_nao_saturado_agua_ar3.tex) and [A] Lira Cecílio (2026), J. Eng. Math. 157:10.
 *
 * Unknowns (TOTAL values in the mesh solution), Eq. (1.1): displacement u (datavec[0], H1, quadratic), water
 * pressure p_w (datavec[1], H1, linear) and air pressure p_a (datavec[2], H1, linear); gauge pressures positive in
 * compression, stresses positive in tension, strains in Voigt form with engineering shear (F Sect. 1). Units: kN, m, s.
 *
 * Field equations of [F] and where they enter:
 *  - kinematics and elastic constants: Eqs. (2.1)-(2.3);
 *  - suction, saturation and retention capacity: Eq. (3.1); retention curves (3.2)-(3.5) in TPZUnsaturatedHydraulics;
 *  - Bishop equivalent pressure p_B: Eq. (3.6) with chi = S_w (default, EChiSaturation); the options chi = S_e and
 *    a constant chi (tan(phi^b)/tan(phi'), Fredlund et al. 1978, as in Rahardjo et al. 2007) are kept (F Sect. 3);
 *    tangent Bishop coefficients chi_w, chi_a: Eq. (3.7);
 *  - equilibrium and effective stress: Eq. (4.1), with the weight factor lambda_g and the pressure factor lambda_p
 *    of F Sect. 4 (both 1 in the transient analysis; see SetGravityFactor and SetGravityScalesSuction);
 *  - Darcy's law with the conductivity tensor K and the mobilities: Eqs. (5.1)-(5.2); relative permeabilities
 *    (5.3)-(5.5) in TPZUnsaturatedHydraulics;
 *  - fluid compressibility and the air density-gradient term: Eq. (6.1);
 *  - storage and coupling coefficients: Eqs. (7.1)-(7.3); mass balances: Eq. (8.1);
 *  - effective stress: the plastic step T (TPZPlasticStepVoigt with the Mohr-Coulomb PV2 of [A], F Sect. 10),
 *    with the strength reduction of F Sect. 10.7.
 *
 * Discretization ([F] Sects. 11-20):
 *  - weak form: Eqs. (11.1)-(11.2); Taylor-Hood triangles T6/T3 (12.1)-(12.4), u in the hierarchical quadratic basis
 *    of NeoPZ; quadrature order set by SetIntegrationOrder (the test projects use the 3-point rule of [N] and [A]);
 *  - storage condensed at the vertices of the P1 pressure element, the operator <.> of Eq. (13.2) with the diagonal
 *    derivative (13.5) (SetLumpedStorage; without condensation the terms are integrated at the points);
 *  - water storage in the secant form (14.4), which conserves the water mass exactly (14.5); air storage (15.2) with
 *    S_a at the previous step, which conserves the air mass exactly (15.4)-(15.5);
 *  - backward-Euler step: Eqs. (16.4)-(16.6) (Theta of (16.6) is the variable X of the code); residuals R_u, R_w, R_a
 *    of Eqs. (16.10)-(16.12), the mass balances multiplied by Delta t; ef = -R (NeoPZ convention, table of F Sect. 16);
 *  - tangent: the full consistent Jacobian dR/dX, Eqs. (17.1)-(17.14) (full Newton), checked by the Taylor test
 *    (17.16); it is not symmetric (F Sect. 18), hence a non-symmetric LU in the analysis;
 *  - Newton iteration, Dirichlet elimination and convergence criterion (19.1)-(19.3): TPZUnsatAnalysis;
 *  - factors of safety by strength reduction (20.1) and gravity increase (20.2) with the pore pressures of an
 *    instant frozen: SetFrozenPressures(true) assembles only the momentum rows (R_u, J_uu and the coupling J_up,
 *    Eqs. (17.8)-(17.9)); the analysis eliminates the pressure equations, so only J_uu is solved (F Sects. 17, 20).
 *
 * Simplifying hypotheses (the items Dk are those of the table of discrepancies, F Sect. 21):
 *  - small strains, quasi-static loading, constant temperature, two immiscible fluid phases without phase change, Darcy
 *    flow with relative permeabilities and Bishop effective stress [N Sect. 1] (F Sect. 1);
 *  - elastoplastic skeleton (associated Mohr-Coulomb of [A]) instead of the elastic skeleton of the notes (F Sect. 1);
 *  - plane strain (eps_xz = eps_yz = eps_zz = 0; the plastic step works with the 3D tensor and returns sigma_zz);
 *  - Biot coefficient alpha and grain modulus K_s constant (D2); the tests use alpha = 1 and incompressible grains
 *    (eta = 0);
 *  - porosity constant in all the coefficients, n = n0: the change of the pore volume enters only through the
 *    coupling term, and there is no geometric stiffness G_uu (D3);
 *  - conductivity tensor K constant (isotropic k_s I or anisotropic), no Kozeny-Carman dependence on the porosity (D4);
 *  - unit weight constant by default; optionally that of the mixture (1-n) gamma_s + n S_w gamma_w, without the weight
 *    of the air and with n constant (D5);
 *  - unit weight of the air constant in the gravity term of the air flux (D6);
 *  - the water density-gradient term c_w q_w . grad p_w of Eq. (6.1) is omitted, as the notes allow (D10);
 *  - water storage in the secant form instead of the rate form (14.3) of the notes (D8); condensed storage evaluated
 *    with the vertex pressures instead of diag(int C N_i) (D13);
 *  - default hydraulic laws of TPZUnsaturatedHydraulics (Fredlund-Xing retention, k_rw = S_e^4, Brooks-Corey air)
 *    instead of the van Genuchten-Mualem of the notes, which stays available, as the Gardner law (D7);
 *  - lower bounds k_rw >= 1e-12 and k_ra >= 1e-4 (zero derivative at the bound, D7); the air pressure is an unknown in
 *    the whole domain, also where S_a -> 0 (F Sects. 5 and 8);
 *  - associated Mohr-Coulomb only (psi = phi, D14); the strength reduction divides c, tan(phi) and tan(psi) by F (D15);
 *  - the prescribed fluxes (EFluxW, EFluxA) are positive entering the domain (D20);
 *  - in the analyses: initial state computed by the model (geostatic step with frozen pressures, displacements then
 *    reset keeping the effective stress, ResetDisplacementReference) instead of given sigma_0 and p_B0 (D16); factors
 *    of safety with the pore pressures of an instant frozen and, in the gravity increase, the suction not scaled
 *    (D18); Newton convergence on the momentum residual and the pressure correction only (Eq. (19.3), D19).
 *
 * The Dirichlet conditions are imposed by the analysis (equation filter, TPZUnsatAnalysis); ContributeBC assembles
 * only tractions (ENeumannU) and fluid fluxes (EFluxW, EFluxA).
 *
 * @tparam T plastic step with ApplyStrainComputeSigma(eps, sigma, tangent) in total strain (engineering shear
 * strains), e.g. TPZPlasticStepVoigt<TPZYCMohrCoulombPV2, TPZElasticResponse>
 * @tparam TMEM memory type (TPZUnsaturatedMem)
 */
template <class T, class TMEM = TPZUnsaturatedMem>
class TPZMatUnsaturatedUPwPa : public TPZMatBase<STATE, TPZMatCombinedSpacesT<STATE>, TPZMatWithMem<TMEM>> {
    using TBase = TPZMatBase<STATE, TPZMatCombinedSpacesT<STATE>, TPZMatWithMem<TMEM>>;

public:
    /** @brief Boundary condition types (the Dirichlet conditions are imposed by the analysis) */
    enum EBCType {
        EDirichletU = 0,   ///< displacement prescribed (imposed by the analysis; no contribution)
        ENeumannU = 1,     ///< traction Val2[0..1]
        EDirichletPw = 2,  ///< water pressure prescribed (imposed by the analysis; no contribution)
        EDirichletPa = 3,  ///< air pressure prescribed (imposed by the analysis; no contribution)
        EFluxW = 4,        ///< inflow of water Val2[0] (m/s, volume per unit area, positive entering)
        EFluxA = 5         ///< inflow of air Val2[0] (m/s)
    };

    /**
     * @brief Bishop parameter: \f$\chi = S_w\f$, constant, or \f$\chi = S_e = (S_w - S_{wr})/(1 - S_{wr})\f$.
     * The constant \f$\chi\f$ (\f$\tan\phi^b = \chi\tan\phi'\f$) acts on the suction only: a positive
     * \f$p_w - p_a\f$ enters with \f$\chi = 1\f$ (Terzaghi), as in the extended Mohr-Coulomb criterion.
     */
    enum EBishop { EChiSaturation = 0, EChiConstant = 1, EChiEffective = 2 };

    /** @brief Post-processing variables */
    enum ESolutionVar {
        EDisplacement = 1,   ///< ("Displacement")
        EWaterPressure = 2,  ///< ("WaterPressure")
        EAirPressure = 3,    ///< ("AirPressure")
        ESuction = 4,        ///< \f$p_a - p_w\f$ ("Suction")
        ESaturation = 5,     ///< \f$S_w\f$ ("Saturation")
        EKrw = 6,            ///< \f$k_{rw}\f$ ("RelativePermeabilityWater")
        EBishopPressure = 7, ///< \f$p_B\f$ ("BishopPressure")
        EDisplacementX = 8,  ///< ("DisplacementX")
        EDisplacementY = 9,  ///< ("DisplacementY")
        // variables of the integration points (memory of the last converged step, zero at other points; meant to be
        // projected with TPZPostProcAnalysis); names of TPZMatElastoPlastic, stresses positive in tension
        EStressXX = 10,         ///< effective stress \f$\sigma'_{xx}\f$ ("StressXX")
        EStressYY = 11,         ///< ("StressYY")
        EStressZZ = 12,         ///< ("StressZZ")
        EStressXY = 13,         ///< ("StressXY")
        EStrainPlastic = 14,    ///< plastic strain XX, YY, ZZ, XY, XZ, YZ, engineering shear ("StrainPlastic")
        EStrainPlasticI1 = 15,  ///< volumetric plastic strain ("StrainPlasticI1")
        EStrainPlasticJ2 = 16,  ///< second invariant of the deviatoric plastic strain tensor ("StrainPlasticJ2")
        EFailureType = 17,      ///< type of the last plastic response, 0 elastic ("FailureType")
        EHardening = 18         ///< hardening variable of the plastic step ("Hardening")
    };

    TPZMatUnsaturatedUPwPa();
    explicit TPZMatUnsaturatedUPwPa(int id);
    TPZMatUnsaturatedUPwPa(const TPZMatUnsaturatedUPwPa &cp) = default;
    virtual ~TPZMatUnsaturatedUPwPa() = default;

    /** @name Set up */
    /** @{ */
    /** @brief Plastic model; it is also used to build the default memory item */
    void SetPlasticModel(const T &model);
    T &PlasticModel() { return fPlasticity; }
    const T &PlasticModel() const { return fPlasticity; }

    void SetHydraulics(const TPZUnsaturatedHydraulics &h) { fHydro = h; }
    const TPZUnsaturatedHydraulics &Hydraulics() const { return fHydro; }

    /** @brief Porosity n */
    void SetPorosity(REAL n) { fPorosity = n; }
    /** @brief Biot coefficient and bulk modulus of the grains (Ks <= 0: incompressible grains, eta = 0) */
    void SetBiot(REAL alpha, REAL Ks) { fAlpha = alpha; fKs = Ks; }
    /** @brief Bulk modulus of the water (kPa) */
    void SetWaterBulkModulus(REAL Kw) { fKw = Kw; }
    /** @brief Saturated hydraulic conductivity k_s (m/s, isotropic), unit weight of the water and mu_w / mu_a */
    void SetConductivity(REAL ks, REAL gammaw, REAL viscosityRatio) {
        fKx = fKy = ks;
        fKxy = 0.;
        fGammaW = gammaw;
        fMuRatio = viscosityRatio;
    }
    /**
     * @brief Anisotropic saturated conductivity (m/s): principal values k_1 (direction at the angle theta, radians,
     * counterclockwise from the x axis) and k_2 (normal to it), K = R diag(k_1, k_2) R^T. The relative
     * permeabilities multiply the tensor; the air uses the same intrinsic permeability tensor (times mu_w / mu_a).
     * theta = 0: K = diag(k_x, k_y) with k_x = k_1, k_y = k_2
     */
    void SetAnisotropicConductivity(REAL k1, REAL k2, REAL theta = 0.) {
        const REAL c = std::cos(theta), s = std::sin(theta);
        fKx = k1 * c * c + k2 * s * s;
        fKy = k1 * s * s + k2 * c * c;
        fKxy = (k1 - k2) * c * s;
    }
    /** @brief Unit weight of the air (kN/m3, constant) and atmospheric pressure (kPa) */
    void SetAir(REAL gammaa, REAL patm) { fGammaA = gammaa; fPatm = patm; }
    /** @brief Constant unit weight of the soil (kN/m3) */
    void SetUnitWeight(REAL gamma) { fGamma = gamma; fMixture = false; }
    /** @brief Unit weight of the mixture (1-n) gamma_s + n S_w gamma_w; gamma_s is the unit weight of the grains */
    void SetUnitWeightMixture(REAL gammas) { fGammaS = gammas; fMixture = true; }
    /** @brief Unit vector of gravity (default (0, -1, 0)) */
    void SetGravityDirection(const TPZVec<REAL> &g) { for (int i = 0; i < 3; ++i) fG[i] = i < g.size() ? g[i] : 0.; }
    /** @brief Bishop parameter */
    void SetBishop(EBishop type, REAL chi = 1.) { fBishop = type; fChi = chi; }
    /** @brief Lumped (default) or consistent storage */
    void SetLumpedStorage(bool lumped) { fLumped = lumped; }
    /**
     * @brief Hydro-mechanical coupling of the mass balances (default true). With false the terms
     * \f$S_\pi(\alpha\Delta\varepsilon_v + \eta\Delta p_B)\f$ are dropped: the seepage is that of an uncoupled analysis
     * (SEEP/W-like, rigid skeleton), the momentum balance still uses the pressures (one-way coupling)
     */
    void SetHydroMechanicalCoupling(bool coupled) { fCoupled = coupled; }
    /** @brief Order of the integration rule (0: NeoPZ default); must be set before the memory is created */
    void SetIntegrationOrder(int order) { fIntegrationOrder = order; }
    /** @} */

    /** @name State of the analysis */
    /** @{ */
    void SetTimeStep(REAL dt) { fTimeStep = dt; }
    REAL TimeStep() const { return fTimeStep; }
    /** @brief Factor of the weight of the soil and of the pore pressures in the momentum balance */
    void SetGravityFactor(REAL lambda) { fGravityFactor = lambda; }
    REAL GravityFactor() const { return fGravityFactor; }
    /**
     * @brief If true (default) the gravity factor multiplies the equivalent pressure p_B of the momentum balance
     * (geostatic initialization with frozen pressures); if false it multiplies only its positive part, the suction
     * being kept (gravity increase: the suction is not gravitational, and scaling it with the weight would make the
     * strength grow with lambda; the positive pore pressures come from the weight of the water)
     */
    void SetGravityScalesSuction(bool flag) { fScaleSuction = flag; }
    /** @brief Only the momentum balance (pore pressures prescribed by the analysis) */
    void SetFrozenPressures(bool frozen) { fFrozen = frozen; }
    bool FrozenPressures() const { return fFrozen; }
    /** @brief When true, Contribute assembles only the weight of the soil (normalization of the momentum residual) */
    void SetAssembleGravityOnly(bool flag) { fGravityOnly = flag; }
    /** @} */

    REAL Porosity() const { return fPorosity; }
    REAL ConductivityX() const { return fKx; }
    REAL ConductivityY() const { return fKy; }
    REAL ConductivityXY() const { return fKxy; }
    REAL GammaW() const { return fGammaW; }


    /**
     * @brief Resets the reference of the displacements keeping the effective stress: at every integration point
     * eps_p <- eps_p - eps_t and eps_t <- 0 (total-strain plastic step); the analysis then sets u = 0
     */
    void ResetDisplacementReference();

    /** @brief Number of integration points with plastic response in the last update of the memory */
    int64_t NPlasticPoints() const;

    /** @name TPZMaterial interface */
    /** @{ */
    int Dimension() const override { return 2; }
    int NStateVariables() const override { return 2; }
    int IntegrationRuleOrder(const TPZVec<int> &elPMaxOrder) const override;
    void FillDataRequirements(TPZVec<TPZMaterialDataT<STATE>> &datavec) const override;
    void FillBoundaryConditionDataRequirements(int type, TPZVec<TPZMaterialDataT<STATE>> &datavec) const override;
    void Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek,
                    TPZFMatrix<STATE> &ef) override;
    void Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ef) override;
    void ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek,
                      TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc) override;
    void ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ef,
                      TPZBndCondT<STATE> &bc) override;
    int VariableIndex(const std::string &name) const override;
    int NSolutionVariables(int var) const override;
    void Solution(const TPZVec<TPZMaterialDataT<STATE>> &datavec, int var, TPZVec<STATE> &sol) override;
    TPZMaterial *NewMaterial() const override { return new TPZMatUnsaturatedUPwPa(*this); }
    std::string Name() const override { return "TPZMatUnsaturatedUPwPa"; }
    void Print(std::ostream &out = std::cout) const override;
    int ClassId() const override;
    void Write(TPZStream &buf, int withclassid) const override;
    void Read(TPZStream &buf, void *context) override;
    /** @} */

    /** @brief Equivalent pressure and its derivatives chi_w = dp_B/dp_w, chi_a = dp_B/dp_a */
    void Bishop(REAL pw, REAL pa, REAL S, REAL Gamma, REAL &pB, REAL &chiw, REAL &chia) const;
    /** @brief Unit weight of the soil and its derivative with respect to p_w */
    void UnitWeight(REAL S, REAL Gamma, REAL &gamma, REAL &dgamma) const;

protected:
    void ContributeInternal(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> *ek,
                            TPZFMatrix<STATE> &ef);
    void ContributeBCInternal(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> *ek,
                              TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc);
    /** @brief Global gradients (2 x n) of shape functions given in the axes of the element */
    static void GlobalGradients(const TPZMaterialDataT<STATE> &data, TPZFMatrix<REAL> &grad);
    T fPlasticity;
    TPZUnsaturatedHydraulics fHydro;
    REAL fPorosity = 0.45, fAlpha = 1., fKs = 0., fKw = 2.2e6;
    REAL fKx = 1.e-5, fKy = 1.e-5, fKxy = 0., fGammaW = 9.81, fMuRatio = 1.e-3 / 1.8e-5;
    REAL fGammaA = 1.2 * 9.81e-3, fPatm = 101.325;
    REAL fGamma = 20., fGammaS = 26.;
    bool fMixture = false;
    REAL fG[3] = {0., -1., 0.};
    EBishop fBishop = EChiSaturation;
    REAL fChi = 1.;
    bool fLumped = true;
    int fIntegrationOrder = 0;
    REAL fTimeStep = 0., fGravityFactor = 1.;
    bool fFrozen = false;
    bool fGravityOnly = false;
    bool fScaleSuction = true;
    bool fCoupled = true;
};

#endif
