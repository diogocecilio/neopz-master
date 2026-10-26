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
 * @brief u - p_w - p_a material of the notes "Formulação poroelástica não saturada: esqueleto deformável, água e ar"
 * (Sects. 1-13), with a Mohr-Coulomb (or any TPZPlasticStepVoigt-like) effective stress.
 *
 * Unknowns (TOTAL values in the mesh solution): displacement u (datavec[0], H1, quadratic), water pressure p_w
 * (datavec[1], H1, linear) and air pressure p_a (datavec[2], H1, linear), gauge pressures positive in compression
 * (kPa); stresses positive in tension. Units: kN, m, s.
 *
 * Effective (Bishop) stress and equivalent pressure:
 * \f[ \sigma = \sigma' - \alpha\,p_B\,I, \qquad p_B = \chi p_w + (1-\chi) p_a, \f]
 * with \f$\chi = S_w\f$ (notes) or a constant \f$\chi = \tan\phi^b/\tan\phi'\f$ (shear strength of Fredlund et al.
 * 1978 used by Rahardjo et al. 2007; \f$\phi^b = \phi'\f$ gives \f$\chi = 1\f$).
 *
 * Backward-Euler residuals of the step (the mass balances are multiplied by \f$\Delta t\f$; \f$\Delta(\cdot)\f$ is the
 * increment in the step, \f$X = \alpha\Delta\varepsilon_v + \eta\Delta p_B\f$, \f$\eta=(\alpha-n)/K_s\f$):
 * \f[ R_u = \int B^T\sigma' - \lambda_g\int\alpha\,m^TB\,p_B - \lambda_g\int N_u^T\gamma\,\mathbf{g} - f_t, \f]
 * \f[ R_w = \int N_p\Big[S_w X + \frac{nS_w}{K_w}\Delta p_w\Big] + \langle n\Delta S_w\rangle
 *     + \Delta t\int\nabla N_p\cdot\lambda_w(\nabla p_w - \gamma_w\mathbf{g}) - \Delta t f_w, \f]
 * \f[ R_a = \int N_p\Big[S_a X - \Delta t\,c_a\lambda_a(\nabla p_a - \gamma_a\mathbf{g})\cdot\nabla p_a\Big]
 *     + \Big\langle n\Delta S_a + \frac{nS_a}{p_a^{abs}}\Delta p_a\Big\rangle
 *     + \Delta t\int\nabla N_p\cdot\lambda_a(\nabla p_a - \gamma_a\mathbf{g}) - \Delta t f_a, \f]
 * where \f$\lambda_w = \mathbf{K} k_{rw}/\gamma_w\f$, \f$\lambda_a = (\mathbf{K}/\gamma_w)(\mu_w/\mu_a)k_{ra}\f$ with the
 * conductivity tensor \f$\mathbf{K} = R(\theta)\,\mathrm{diag}(k_1, k_2)\,R(\theta)^T\f$ (isotropic: \f$k_1 = k_2 = k_s\f$),
 * \f$c_a = 1/p_a^{abs}\f$ (ideal isothermal gas), \f$\mathbf{g}\f$ the unit vector of gravity, \f$\lambda_g\f$ the gravity
 * factor (1 in the transient analysis; the geostatic initialization scales the weight of the soil and the pore
 * pressures together, the gravity-increase factor of safety the weight and the positive pressures, see SetGravityScalesSuction) and \f$\langle\cdot\rangle\f$ the storage terms, lumped at the
 * vertices of the P1 pressure element (\f$\int N_i\,d\Omega\f$ times the term evaluated with the vertex pressures,
 * conservative secant form \f$n(S^{n+1}-S^n)\f$) or, if lumping is disabled, integrated at the points. The porosity
 * is constant in the coefficients (first-order closure of Sect. 9 with \f$\alpha = 1\f$, incompressible grains).
 * The unit weight is constant (Rahardjo et al. 2007) or that of the mixture
 * \f$\gamma = (1-n)\gamma_s + nS_w\gamma_w\f$. The tangent is the full consistent Jacobian of these residuals.
 *
 * With SetFrozenPressures(true) only the momentum balance is assembled (factor of safety with the pore pressures of
 * an instant; the pressure equations must then be prescribed by the analysis). The Dirichlet conditions are
 * imposed by the analysis (equation filter); ContributeBC assembles only tractions (ENeumannU) and fluid
 * fluxes (EFluxW, EFluxA).
 *
 * As usual in NeoPZ the element vector is ef = -R and the element matrix is the tangent dR/dx.
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

    /** @brief Bishop parameter: \f$\chi = S_w\f$ or constant */
    enum EBishop { EChiSaturation = 0, EChiConstant = 1 };

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
        EDisplacementY = 9   ///< ("DisplacementY")
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

protected:
    void ContributeInternal(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> *ek,
                            TPZFMatrix<STATE> &ef);
    void ContributeBCInternal(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> *ek,
                              TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc);
    /** @brief Global gradients (2 x n) of shape functions given in the axes of the element */
    static void GlobalGradients(const TPZMaterialDataT<STATE> &data, TPZFMatrix<REAL> &grad);
    /** @brief Equivalent pressure and its derivatives chi_w = dp_B/dp_w, chi_a = dp_B/dp_a */
    void Bishop(REAL pw, REAL pa, REAL S, REAL Gamma, REAL &pB, REAL &chiw, REAL &chia) const;
    /** @brief Unit weight of the soil and its derivative with respect to p_w */
    void UnitWeight(REAL S, REAL Gamma, REAL &gamma, REAL &dgamma) const;

    T fPlasticity;
    TPZUnsaturatedHydraulics fHydro;
    REAL fPorosity = 0.45, fAlpha = 1., fKs = 0., fKw = 2.2e6;
    REAL fKx = 1.e-5, fKy = 1.e-5, fKxy = 0., fGammaW = 9.81, fMuRatio = 1.e-3 / 1.8e-5;
    REAL fGammaA = 1.2 * 9.81e-3, fPatm = 101.325;
    REAL fGamma = 20., fGammaS = 26.;
    bool fMixture = false;
    REAL fG[3] = {0., -1., 0.};
    EBishop fBishop = EChiConstant;
    REAL fChi = 1.;
    bool fLumped = true;
    int fIntegrationOrder = 0;
    REAL fTimeStep = 0., fGravityFactor = 1.;
    bool fFrozen = false;
    bool fGravityOnly = false;
    bool fScaleSuction = true;
};

#endif
