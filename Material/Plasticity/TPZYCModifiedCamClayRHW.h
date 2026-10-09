/**
 * @file TPZYCModifiedCamClayRHW.h
 * @brief Modified Cam-Clay yield criterion and closest-point projection in rotated
 * Haigh-Westergaard (RHW) space.
 *
 * Implements the local problem of "Two-variable closest-point projection of the Modified Cam-Clay model in
 * rotated Haigh-Westergaard space with consistent tangent operator" (D. Lira Cecilio): the reduced system in
 * the angle of the meridian ellipse and the hardening increment (ProjectReduced), and, as a cross-check, the
 * four-equation system of the Wolfram Language routines HardeningCC, PhiCC, ResCC, JacCC, dResdTrialCC,
 * ProjectHWCC and GradCC of camclay-perf.m (ProjectHW).
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
 * The yield function, in the form of de Souza Neto et al., is
 * \f[ \Phi(\xi,\rho,a) = \frac{1}{b^2}\left(\frac{\xi}{\sqrt3}-p_t+a\right)^2
 *     + \frac{3\rho^2}{2M^2} - a^2, \qquad a = \frac{p_c+p_t}{1+\omega}, \f]
 * with \f$b=1\f$ in the supercritical region (\f$\bar p\ge0\f$) and \f$b=\omega\f$ in the
 * subcritical region (\f$\bar p<0\f$), \f$\bar p = p - p_t + a\f$. The preconsolidation pressure
 * follows the exponential hardening law
 * \f$ p_c = p_{c,n}\exp\left(\frac{v_0\,\Delta\alpha}{\lambda-\kappa}\right)\f$.
 *
 * The local problem is the closest-point projection of the trial state onto the end-of-step surface in the
 * complementary-energy metric of the RHW space. The distance is stationary in the Lode angle at
 * \f$\beta=\beta_{tr}\f$ (radial return in the deviatoric plane), and with the meridian ellipse parametrized by
 * the angle \f$\theta\f$,
 * \f[ \bar p = -a\,b\cos\theta,\qquad \rho = \sqrt{2/3}\,M a\sin\theta,\qquad
 *     \xi = \sqrt3\,(p_t - a + \bar p), \f]
 * (\f$\theta=0\f$ at the compression apex, \f$\theta=\pi/2\f$ on the critical state line) the problem reduces to
 * two unknowns \f$X=[\theta,\Delta\alpha]\f$ and two equations: the stationarity of the distance along the
 * surface for the end-of-step value of \f$a\f$ (E1) and the volumetric elastic law, which defines the hardening
 * increment \f$\Delta\alpha=-\Delta\varepsilon^p_v\f$ (E2). The plastic multiplier is not an unknown; it is
 * recovered from \f$\Delta\gamma=M^2(\rho_{tr}/\rho-1)/(6G)\f$ (reduced solver, the default, ProjectReduced).
 * The four-equation system in \f$[\xi,\rho,\Delta\alpha,\Delta\gamma]\f$ of the earlier versions (ProjectHW)
 * is kept as a cross-check (SetLocalSolver(EFull)); both have the same solution and the same consistent tangent.
 * The Jacobian of the projection in principal stresses, \f$\partial\sigma^{proj}/\partial\sigma^{tr}\f$, is
 * obtained by implicit differentiation of the converged local system (GradProjectionReduced, GradProjection).
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

    /** @brief Volumetric elastic law of the elastic predictor and of the local problem */
    enum EVolumetricLaw {
        ELinear = 0, ///< linear elasticity, \f$p = p_n + K_0\,\Delta\varepsilon^e_v\f$
        EPorous = 1  ///< porous (logarithmic) elasticity, \f$p = p_n\exp(-v_0\Delta\varepsilon^e_v/\kappa)\f$
    };

    /** @brief Integration of the porous law during the plastic correction (comparison with the frozen bulk modulus) */
    enum EPorousIntegration {
        EExact = 0, ///< exact integral of the porous law over the step (this work)
        EFrozen = 1 ///< bulk modulus frozen at its trial value (Sanei et al. 2020), comparison only
    };

    /** @brief System solved by the local Newton iterations */
    enum ELocalSolver {
        EReduced = 0, ///< two unknowns \f$[\theta,\Delta\alpha]\f$, distance stationary along the surface (default)
        EFull = 1     ///< four unknowns \f$[\xi,\rho,\Delta\alpha,\Delta\gamma]\f$ (cross-check)
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

    /** @brief Default constructor (M = 1, lambda = 0.2, kappa = 0.05, pt = 0, omega = 1, porous law, reduced solver) */
    TPZYCModifiedCamClayRHW();

    /** @brief Copy constructor */
    TPZYCModifiedCamClayRHW(const TPZYCModifiedCamClayRHW &cp) = default;

    /** @brief Assignment operator */
    TPZYCModifiedCamClayRHW &operator=(const TPZYCModifiedCamClayRHW &cp) = default;

    virtual ~TPZYCModifiedCamClayRHW() = default;

    /**
     * @brief Sets the material parameters
     * @param M slope of the critical state line
     * @param lambda slope of the normal compression line in the v-ln p' plane
     * @param kappa slope of the swelling line in the v-ln p' plane
     * @param pt tensile strength \f$p_t\f$ (0 in all the examples)
     * @param omega shape parameter of the subcritical region (1 in all the examples)
     */
    void SetUp(REAL M, REAL lambda, REAL kappa, REAL pt = 0., REAL omega = 1.);

    /** @brief Selects the volumetric elastic law (linear needs the bulk modulus K0) */
    void SetVolumetricLaw(EVolumetricLaw law, REAL K0 = 0.) { fVolumetricLaw = law; fK0 = K0; }

    /** @brief Selects exact or frozen integration of the porous law (frozen only for the comparison with Sanei et al. 2020) */
    void SetPorousIntegration(EPorousIntegration integ) { fPorousIntegration = integ; }

    /** @brief Tolerance and maximum number of iterations of the local Newton method */
    void SetNewtonParameters(REAL tol, int maxit) { fNewtonTol = tol; fMaxNewton = maxit; }

    /** @brief Selects the local system: reduced (two unknowns, default) or full (four unknowns, cross-check) */
    void SetLocalSolver(ELocalSolver solver) { fLocalSolver = solver; }

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
    ELocalSolver LocalSolver() const { return fLocalSolver; }
    REAL NewtonTolerance() const { return fNewtonTol; }
    /** @} */

    /**
     * @brief Hardening law: semi-axis a, hardening modulus H = da/d(dal) and preconsolidation pressure
     * @param pcn preconsolidation pressure of the last converged step
     * @param dal hardening increment \f$\Delta\alpha = -\Delta\varepsilon^p_v\f$
     * @param v0 specific volume
     * @param[out] a semi-axis of the ellipse at the end of the step
     * @param[out] H \f$da/d\Delta\alpha\f$
     * @param[out] pc preconsolidation pressure at the end of the step
     */
    void Hardening(REAL pcn, REAL dal, REAL v0, REAL &a, REAL &H, REAL &pc) const;

    /**
     * @brief Yield function in terms of the mean stress and the deviatoric radius
     * @param p mean stress (tension positive)
     * @param rho deviatoric radius \f$\rho=\sqrt{2J_2}\f$
     * @param a semi-axis of the ellipse
     * @param b shape parameter (1 or omega)
     */
    REAL PhiCC(REAL p, REAL rho, REAL a, REAL b) const;

    /** @brief Shape parameter b of the region of the state (mean stress p, semi-axis a) */
    REAL BFromP(REAL p, REAL a) const { return (p - fPt + a >= 0.) ? 1. : fOmega; }

    /**
     * @brief Residual of the full local problem (four equations: volumetric elastic law, deviatoric elastic law
     * with the flow rule, hardening rule, consistency), scaled by \f$a_n\f$ and \f$a_n^2\f$ (ResCC)
     * @param X unknowns \f$[\xi,\rho,\Delta\alpha,\Delta\gamma]\f$
     * @param trial trial state data
     * @param[out] R residual vector (4)
     */
    void Residual(const TPZVec<REAL> &X, const TTrial &trial, TPZVec<REAL> &R) const;

    /**
     * @brief Jacobian of the full residual and derivative of the residual with respect to the
     * trial coordinates \f$Y=[\xi_{tr},\rho_{tr}]\f$ (JacCC and dResdTrialCC)
     * @param X unknowns
     * @param trial trial state data
     * @param[out] J 4x4 Jacobian \f$\partial R/\partial X\f$
     * @param[out] dRdY 4x2 matrix \f$\partial R/\partial Y\f$
     */
    void Jacobian(const TPZVec<REAL> &X, const TTrial &trial, TPZFMatrix<REAL> &J, TPZFMatrix<REAL> &dRdY) const;

    /**
     * @brief Local Newton iterations of the full system (ProjectHWCC) starting from \f$X^{(0)}=[\xi_{tr},\rho_{tr},0,0]\f$
     * @param trial trial state data (the branch b is taken from trial.fB)
     * @param[out] X converged unknowns
     * @param[out] niter number of Newton corrections
     * @return true if \f$\|R\|\f$ reached the tolerance
     */
    bool ProjectHW(const TTrial &trial, TPZVec<REAL> &X, int &niter) const;

    /**
     * @brief Jacobian of the projection in principal stresses from the full system (GradCC)
     * @param X converged unknowns of the full system
     * @param trial trial state data
     * @param n unit deviatoric direction of the trial state in principal stresses (zero if isotropic)
     * @param isotropic true when the trial state is isotropic (\f$\rho_{tr}=0\f$)
     * @param[out] Dproj 3x3 matrix \f$D_{ij}=\partial\sigma^{proj}_i/\partial\sigma^{tr}_j\f$
     * @return false if the local Jacobian is singular (Dproj is not computed)
     */
    bool GradProjection(const TPZVec<REAL> &X, const TTrial &trial, const TPZVec<REAL> &n, bool isotropic,
                        TPZFMatrix<REAL> &Dproj) const;

    /** @name Reduced local problem (two unknowns) */
    /** @{ */

    /**
     * @brief Hardening increment given by the volumetric elastic law, \f$\Delta\alpha_e(\xi)\f$, and its derivatives
     *
     * The volumetric elastic law integrated over the step relates the projected and the trial hydrostatic
     * coordinates through the plastic volumetric strain \f$-\Delta\alpha\f$:
     * \f$\Delta\alpha_e=-(\kappa/v_0)\ln(\xi/\xi_{tr})\f$ (porous law, exact), \f$(\kappa/v_0)(1-\xi/\xi_{tr})\f$
     * (porous law, bulk modulus frozen at the trial value) or \f$(\xi-\xi_{tr})/(\sqrt3K_0)\f$ (linear law).
     * \f$\Delta\alpha_e/\sqrt3\f$ is the derivative with respect to \f$\xi\f$ of the volumetric part of the
     * distance (the Bregman divergence of the complementary energy of the elastic law).
     * @param xi projected hydrostatic coordinate
     * @param xitr trial hydrostatic coordinate
     * @param v0 specific volume
     * @param[out] dxi \f$\partial\Delta\alpha_e/\partial\xi\f$
     * @param[out] dxitr \f$\partial\Delta\alpha_e/\partial\xi_{tr}\f$
     */
    REAL ElasticHardeningIncrement(REAL xi, REAL xitr, REAL v0, REAL &dxi, REAL &dxitr) const;

    /**
     * @brief Point of the meridian ellipse of semi-axis a at the angle theta, and its derivatives
     * @param theta angle on the ellipse (0 at the compression apex, pi/2 on the critical state line)
     * @param a semi-axis
     * @param b shape parameter of the branch
     * @param[out] xi hydrostatic coordinate \f$\sqrt3(p_t-a-ab\cos\theta)\f$
     * @param[out] rho deviatoric radius \f$\sqrt{2/3}Ma\sin\theta\f$
     * @param[out] dxi derivatives \f$[\partial\xi/\partial\theta,\ \partial\xi/\partial a]\f$
     * @param[out] drho derivatives \f$[\partial\rho/\partial\theta,\ \partial\rho/\partial a]\f$
     */
    void SurfacePoint(REAL theta, REAL a, REAL b, REAL &xi, REAL &rho, REAL dxi[2], REAL drho[2]) const;

    /**
     * @brief Residual of the reduced local problem
     * @param X2 unknowns \f$[\theta,\Delta\alpha]\f$
     * @param trial trial state data
     * @param[out] E residual (2): E1 = stationarity of the distance along the surface, divided by the semi-axis
     * (a strain); E2 = \f$\Delta\alpha-\Delta\alpha_e(\xi)\f$ (volumetric elastic law)
     */
    void ResidualReduced(const TPZVec<REAL> &X2, const TTrial &trial, TPZVec<REAL> &E) const;

    /**
     * @brief Jacobian of the reduced residual and its derivative with respect to the trial coordinates
     * @param X2 unknowns \f$[\theta,\Delta\alpha]\f$
     * @param trial trial state data
     * @param[out] J 2x2 Jacobian \f$\partial E/\partial X\f$
     * @param[out] dEdY 2x2 matrix \f$\partial E/\partial Y\f$, \f$Y=[\xi_{tr},\rho_{tr}]\f$
     */
    void JacobianReduced(const TPZVec<REAL> &X2, const TTrial &trial, TPZFMatrix<REAL> &J, TPZFMatrix<REAL> &dEdY) const;

    /**
     * @brief Newton iterations of the reduced problem, with backtracking on the norm of the residual, from the
     * angle of the trial state on the ellipse of \f$a_n\f$ and \f$\Delta\alpha=0\f$
     * @param trial trial state data (the branch b is taken from trial.fB)
     * @param[out] X2 converged unknowns \f$[\theta,\Delta\alpha]\f$
     * @param[out] niter number of Newton corrections
     * @return true if \f$\|E\|\f$ reached the tolerance
     */
    bool ProjectReduced(const TTrial &trial, TPZVec<REAL> &X2, int &niter) const;

    /**
     * @brief Quantities of the four-unknown system from the reduced solution:
     * \f$X=[\xi,\rho,\Delta\alpha,\Delta\gamma]\f$ with \f$\Delta\gamma=M^2(\rho_{tr}/\rho-1)/(6G)\f$
     * (or \f$-\Delta\alpha b^2/(2\bar p)\f$ at the apexes, where \f$\rho=0\f$)
     */
    void ReducedToFull(const TPZVec<REAL> &X2, const TTrial &trial, TPZVec<REAL> &X) const;

    /**
     * @brief Jacobian of the projection in principal stresses from the reduced system: implicit differentiation
     * \f$\partial X/\partial Y=-J^{-1}\partial E/\partial Y\f$ and the chain rule through \f$\xi(\theta,a)\f$,
     * \f$\rho(\theta,a)\f$ (same arguments and result as GradProjection)
     */
    bool GradProjectionReduced(const TPZVec<REAL> &X2, const TTrial &trial, const TPZVec<REAL> &n, bool isotropic,
                               TPZFMatrix<REAL> &Dproj) const;

    /**
     * @brief Local projection with the solver selected by SetLocalSolver
     * @param trial trial state data
     * @param[out] X \f$[\xi,\rho,\Delta\alpha,\Delta\gamma]\f$ of the converged state (ReducedToFull for the reduced solver)
     * @param[out] X2 \f$[\theta,\Delta\alpha]\f$ (reduced solver only; otherwise resized to zero)
     * @param[out] niter number of Newton corrections
     * @return true if the iterations converged
     */
    bool Project(const TTrial &trial, TPZVec<REAL> &X, TPZVec<REAL> &X2, int &niter) const;

    /** @brief Jacobian of the projection for the solver selected by SetLocalSolver (GradProjectionReduced or GradProjection) */
    bool ProjectionJacobian(const TPZVec<REAL> &X, const TPZVec<REAL> &X2, const TTrial &trial, const TPZVec<REAL> &n,
                            bool isotropic, TPZFMatrix<REAL> &Dproj) const;
    /** @} */

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
    REAL fNewtonTol; ///< tolerance of the local Newton method on the norm of the residual (1e-12)
    int fMaxNewton;  ///< maximum number of local Newton iterations
    ELocalSolver fLocalSolver; ///< reduced (two unknowns) or full (four unknowns) local system

    /**
     * @brief Jacobian of the projection in principal stresses from the derivatives of the projected coordinates
     * @param dxi \f$[\partial\xi/\partial\xi_{tr},\ \partial\xi/\partial\rho_{tr}]\f$
     * @param drho \f$[\partial\rho/\partial\xi_{tr},\ \partial\rho/\partial\rho_{tr}]\f$
     * @param f ratio \f$\rho/\rho_{tr}\f$ (scaling of the deviatoric components normal to n)
     * @param n unit deviatoric direction of the trial state (zero if isotropic)
     * @param[out] Dproj 3x3 matrix
     */
    static void AssembleDproj(const REAL dxi[2], const REAL drho[2], REAL f, const TPZVec<REAL> &n, TPZFMatrix<REAL> &Dproj);
};

#endif
