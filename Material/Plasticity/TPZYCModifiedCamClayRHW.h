/**
 * @file TPZYCModifiedCamClayRHW.h
 * @brief Modified Cam-Clay yield criterion and closest-point projection in rotated
 * Haigh-Westergaard (RHW) space.
 *
 * Implements Sect. 4 of "Return mapping for Modified Cam-Clay plasticity in rotated
 * Haigh-Westergaard space with consistent tangent operator and coupled u-p consolidation"
 * (D. Lira Cecilio). It is the C++ counterpart of the routines HardeningCC, PhiCC, ResCC,
 * JacCC, dResdTrialCC, ProjectHWCC and GradCC of camclay-perf.m (Table 9 of the article).
 *
 * Conventions: tension positive, \f$p = I_1/3\f$ (negative in compression),
 * \f$\xi = I_1/\sqrt{3}\f$, \f$\rho = \sqrt{2 J_2}\f$, \f$q = \sqrt{3/2}\,\rho\f$.
 */

#ifndef TPZYCMODIFIEDCAMCLAYRHW_H
#define TPZYCMODIFIEDCAMCLAYRHW_H

#include "pzreal.h"
#include "pzvec.h"
#include "pzfmatrix.h"
#include "TPZPlasticState.h"
#include "TPZPlasticCriterion.h"

/**
 * @ingroup material
 * @brief Modified Cam-Clay (MCC) yield criterion with the local return mapping of the article
 * written in rotated Haigh-Westergaard coordinates.
 *
 * The yield function, in the form of de Souza Neto et al., eq. (11)-(12) of the article, is
 * \f[ \Phi(\xi,\rho,a) = \frac{1}{b^2}\left(\frac{\xi}{\sqrt3}-p_t+a\right)^2
 *     + \frac{3\rho^2}{2M^2} - a^2, \qquad a = \frac{p_c+p_t}{1+\omega}, \f]
 * with \f$b=1\f$ in the supercritical region (\f$\bar p\ge0\f$) and \f$b=\omega\f$ in the
 * subcritical region (\f$\bar p<0\f$), \f$\bar p = p - p_t + a\f$. The preconsolidation pressure
 * follows the exponential hardening law (13)
 * \f$ p_c = p_{c,n}\exp\left(\frac{v_0\,\Delta\alpha}{\lambda-\kappa}\right)\f$.
 *
 * Because \f$\Phi\f$ does not depend on the Lode angle the projection is radial in the deviatoric
 * plane and the local problem (18) has four unknowns \f$X=[\xi,\rho,\Delta\alpha,\Delta\gamma]\f$.
 * The class solves (18) by Newton's method with the analytical Jacobian (A.1) and computes the
 * Jacobian of the projection in principal stresses, eq. (22), by implicit differentiation (21).
 *
 * The criterion does not hold the elastic moduli: the trial data (trial RHW coordinates, shear
 * modulus of the step, specific volume and preconsolidation pressure of the previous step) are
 * passed in a TPZYCModifiedCamClayRHW::TTrial structure built by the plastic step
 * (TPZPlasticStepModifiedCamClay).
 */
class TPZYCModifiedCamClayRHW : public TPZPlasticCriterion {
public:

    enum {
        NYield = 1 ///< number of yield functions
    };

    /** @brief Volumetric elastic law used in the residual R1 of (18) */
    enum EVolumetricLaw {
        ELinear = 0, ///< linear elasticity, \f$p = p_n + K_0\,\Delta\varepsilon^e_v\f$
        EPorous = 1  ///< porous (logarithmic) elasticity, \f$p = p_n\exp(-v_0\Delta\varepsilon^e_v/\kappa)\f$
    };

    /** @brief Integration of the porous law during the plastic correction (Sect. 6.1, Table 3) */
    enum EPorousIntegration {
        EExact = 0, ///< exact integral of the porous law over the step (this work)
        EFrozen = 1 ///< bulk modulus frozen at its trial value (Sanei et al. 2020), comparison only
    };

    /** @brief Data of the elastic trial state needed by the local projection */
    struct TTrial {
        REAL fXiTr = 0.;  ///< hydrostatic RHW coordinate of the trial state, \f$\xi_{tr}=\sqrt3\,p_{tr}\f$
        REAL fRhoTr = 0.; ///< deviatoric radius of the trial state, \f$\rho_{tr}\f$
        REAL fG = 0.;     ///< shear modulus of the step (kept fixed during the step)
        REAL fV0 = 1.;    ///< specific volume \f$v_0\f$ of the integration point
        REAL fPcn = 0.;   ///< preconsolidation pressure of the last converged step \f$p_{c,n}\f$
        REAL fB = 1.;     ///< shape parameter \f$b\f$ of the branch (1 or \f$\omega\f$), fixed during the iterations
    };

    /** @brief Default constructor (M = 1, lambda = 0.2, kappa = 0.05, pt = 0, omega = 1, porous law) */
    TPZYCModifiedCamClayRHW();

    /** @brief Copy constructor */
    TPZYCModifiedCamClayRHW(const TPZYCModifiedCamClayRHW &cp) = default;

    /** @brief Assignment operator */
    TPZYCModifiedCamClayRHW &operator=(const TPZYCModifiedCamClayRHW &cp) = default;

    virtual ~TPZYCModifiedCamClayRHW() = default;

    /**
     * @brief Sets the material parameters (Table 1 of the article)
     * @param M slope of the critical state line
     * @param lambda slope of the normal compression line in the v-ln p' plane
     * @param kappa slope of the swelling line in the v-ln p' plane
     * @param pt tensile strength \f$p_t\f$ (0 in all the examples)
     * @param omega shape parameter of the subcritical region (1 in all the examples)
     */
    void SetUp(REAL M, REAL lambda, REAL kappa, REAL pt = 0., REAL omega = 1.);

    /** @brief Selects the volumetric elastic law of the residual R1 (linear needs the bulk modulus K0) */
    void SetVolumetricLaw(EVolumetricLaw law, REAL K0 = 0.) { fVolumetricLaw = law; fK0 = K0; }

    /** @brief Selects exact or frozen integration of the porous law (only for the comparison of Sect. 6.1) */
    void SetPorousIntegration(EPorousIntegration integ) { fPorousIntegration = integ; }

    /** @brief Tolerance and maximum number of iterations of the local Newton method */
    void SetNewtonParameters(REAL tol, int maxit) { fNewtonTol = tol; fMaxNewton = maxit; }

    /** @name Access to the parameters */
    /** @{ */
    REAL M() const { return fM; }
    REAL Lambda() const { return fLambda; }
    REAL Kappa() const { return fKappa; }
    REAL Pt() const { return fPt; }
    REAL Omega() const { return fOmega; }
    REAL K0() const { return fK0; }
    EVolumetricLaw VolumetricLaw() const { return fVolumetricLaw; }
    EPorousIntegration PorousIntegration() const { return fPorousIntegration; }
    /** @} */

    /**
     * @brief Hardening law (13): semi-axis a, hardening modulus H = da/d(dal) and preconsolidation pressure
     * @param pcn preconsolidation pressure of the last converged step
     * @param dal hardening increment \f$\Delta\alpha = -\Delta\varepsilon^p_v\f$
     * @param v0 specific volume
     * @param[out] a semi-axis of the ellipse at the end of the step
     * @param[out] H \f$da/d\Delta\alpha\f$
     * @param[out] pc preconsolidation pressure at the end of the step
     */
    void Hardening(REAL pcn, REAL dal, REAL v0, REAL &a, REAL &H, REAL &pc) const;

    /**
     * @brief Yield function (11) in terms of the mean stress and the deviatoric radius
     * @param p mean stress (tension positive)
     * @param rho deviatoric radius \f$\rho=\sqrt{2J_2}\f$
     * @param a semi-axis of the ellipse
     * @param b shape parameter (1 or omega)
     */
    REAL PhiCC(REAL p, REAL rho, REAL a, REAL b) const;

    /** @brief Shape parameter b of the region of the state (mean stress p, semi-axis a) */
    REAL BFromP(REAL p, REAL a) const { return (p - fPt + a >= 0.) ? 1. : fOmega; }

    /**
     * @brief Residual (18) of the local problem, scaled by \f$a_n\f$ and \f$a_n^2\f$ (ResCC)
     * @param X unknowns \f$[\xi,\rho,\Delta\alpha,\Delta\gamma]\f$
     * @param trial trial state data
     * @param[out] R residual vector (4)
     */
    void Residual(const TPZVec<REAL> &X, const TTrial &trial, TPZVec<REAL> &R) const;

    /**
     * @brief Jacobian (A.1) of the residual and derivative of the residual with respect to the
     * trial coordinates \f$Y=[\xi_{tr},\rho_{tr}]\f$, eq. (21) (JacCC and dResdTrialCC)
     * @param X unknowns
     * @param trial trial state data
     * @param[out] J 4x4 Jacobian \f$\partial R/\partial X\f$
     * @param[out] dRdY 4x2 matrix \f$\partial R/\partial Y\f$
     */
    void Jacobian(const TPZVec<REAL> &X, const TTrial &trial, TPZFMatrix<REAL> &J, TPZFMatrix<REAL> &dRdY) const;

    /**
     * @brief Local Newton iterations of Sect. 4.3 (ProjectHWCC) starting from \f$X^{(0)}=[\xi_{tr},\rho_{tr},0,0]\f$
     * @param trial trial state data (the branch b is taken from trial.fB)
     * @param[out] X converged unknowns
     * @param[out] niter number of Newton corrections
     * @return true if \f$\|R\|\f$ reached the tolerance
     */
    bool ProjectHW(const TTrial &trial, TPZVec<REAL> &X, int &niter) const;

    /**
     * @brief Jacobian of the projection in principal stresses, eq. (22) (GradCC)
     * @param X converged unknowns
     * @param trial trial state data
     * @param n unit deviatoric direction of the trial state in principal stresses (zero if isotropic)
     * @param isotropic true when the trial state is isotropic (\f$\rho_{tr}=0\f$)
     * @param[out] Dproj 3x3 matrix \f$D_{ij}=\partial\sigma^{proj}_i/\partial\sigma^{tr}_j\f$
     */
    void GradProjection(const TPZVec<REAL> &X, const TTrial &trial, const TPZVec<REAL> &n, bool isotropic,
                        TPZFMatrix<REAL> &Dproj) const;

    /**
     * @brief Yield function evaluated at principal stresses
     * @param sigma principal stresses
     * @param kprev preconsolidation pressure \f$p_c\f$
     * @param[out] yield value of \f$\Phi\f$ (size 1)
     */
    void YieldFunction(const TPZVec<STATE> &sigma, STATE kprev, TPZVec<STATE> &yield) const override;

    int GetNYield() const override { return as_integer(NYield); }

    /** @brief Local material parameters are not used by this criterion (no-op) */
    void SetLocalMatState(TPZPlasticState<REAL> &state) override {}

    /** @brief Local material parameters are not used by this criterion */
    TPZPlasticState<REAL> GetLocalMatState() override { return TPZPlasticState<REAL>(); }

    /** @brief Strength reduction is not defined for this criterion (no-op) */
    void ChangeLocalMatParameters(TPZPlasticState<REAL> &state, REAL factor) override {}

    void Print(std::ostream &out) const override;

    int ClassId() const override;

    void Write(TPZStream &buf, int withclassid) const override;

    void Read(TPZStream &buf, void *context) override;

protected:
    REAL fM;      ///< slope of the critical state line
    REAL fLambda; ///< slope of the normal compression line
    REAL fKappa;  ///< slope of the swelling line
    REAL fPt;     ///< tensile strength
    REAL fOmega;  ///< shape parameter of the subcritical region (key "beta" in the WL code)
    REAL fK0;     ///< bulk modulus of the linear volumetric law
    EVolumetricLaw fVolumetricLaw;         ///< volumetric elastic law
    EPorousIntegration fPorousIntegration; ///< exact or frozen porous law (comparison of Sect. 6.1)
    REAL fNewtonTol; ///< tolerance of the local Newton method on \f$\|R\|\f$ (1e-12 as in the Python code)
    int fMaxNewton;  ///< maximum number of local Newton iterations
};

#endif
