/**
 * @file IntegrationSchemes.h
 * @brief Sect. 6.2, Table 4 and Fig. 6 of the article: accuracy and work of the return mapping of this work and
 * of the rival integration schemes of the Modified Cam-Clay model (implicit backward-Euler variants and explicit
 * adaptive Runge-Kutta schemes) in material-point tests with closed-form references. C++ port of rivais.py and
 * of the function rivais() of gen_data.py of the Python code of the article.
 */
#pragma once

#include "MCCPaperTools.h"
#include "IntegrationSchemesReference.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

/**
 * @ingroup mccpaper
 * @brief Comparison of the return mapping of this work with the integration schemes of the recent literature
 * (Sect. 6.2, Table 4 and Fig. 6 of the article).
 *
 * The schemes integrate the same Modified Cam-Clay model (\f$\omega=1\f$, \f$p_t=0\f$) with the same elastic,
 * flow and hardening laws: porous volumetric law \f$p=p_n\exp(-v_0\Delta\varepsilon^e_v/\kappa)\f$, shear modulus
 * \f$G=rK\f$ with \f$r=3(1-2\nu)/(2(1+\nu))\f$ and the hardening law \f$p_c=p_{c,n}\exp(v_0\Delta\alpha/(\lambda-\kappa))\f$
 * (\f$v_0\f$ constant). They are a faithful C++ transcription of rivais.py:
 *
 * 1. Implicit (backward-Euler, BE) schemes written in tensors with the radial return (BackwardEuler, be_tensor):
 *    unknowns \f$X=(\xi,\rho,\Delta\alpha,\Delta\gamma)\f$ and the residuals
 *    \f[ R_1=(\xi-\xi_{tr}e_1)/a_n,\quad R_2=(\rho(1+6G\Delta\gamma/M^2)-\rho_{tr}(G))/a_n,\quad
 *        R_3=\Delta\alpha+2\Delta\gamma\bar p,\quad R_4=(\bar p^2+\tfrac32\rho^2/M^2-a^2)/a_n^2, \f]
 *    \f$\bar p=\xi/\sqrt3-p_t+a\f$, \f$\rho_{tr}(G)=\|s_n+2G\Delta e\|\f$, with the options
 *    - volumetric law EExactK (\f$e_1=\exp(-v_0\Delta\alpha/\kappa)\f$, exact integral, this work) or EFrozenK
 *      (\f$e_1=1-v_0\Delta\alpha/\kappa\f$, bulk modulus frozen at its trial value, Sanei et al. 2020);
 *    - shear modulus EShearAtPn (\f$G=rK(p_n)\f$ kept in the step, this work) or ESecantShear (secant
 *      \f$\bar G=r\bar K\f$, \f$\bar K=p_n(\exp(-v_0x/\kappa)-1)/x\f$, \f$x=\Delta\varepsilon_v+\Delta\alpha\f$: the
 *      implicit systems of Zhou et al. 2022, Lu et al. 2023, Krabbenhoft and Lyamin 2012 and Bui et al. 2026 at
 *      convergence).
 *    With (EExactK, EShearAtPn) the scheme is the return mapping of this work; it is computed also with the library
 *    class TPZPlasticStepModifiedCamClay (LibraryUpdate, spectral form in rotated Haigh-Westergaard space), which
 *    must give the same stresses (to round-off) and the same numbers of local iterations; the frozen variant is
 *    cross-checked with the option TPZYCModifiedCamClayRHW::EFrozen of the library.
 *    BackwardEulerSubstep (be_substep) halves the increment recursively only when the single step fails, as
 *    Bui et al. (2026) do; it is used with the secant shear modulus.
 * 2. Explicit adaptive Runge-Kutta schemes (RungeKutta, rk_update; Sloan et al. 2001, Xie et al. 2026): embedded
 *    pairs ME2(1) and RKDP5(4) (TTableau), exact integration of the elastic part, Newton search (with bisection
 *    safeguard) of the intersection with the yield surface, one drift correction at the end of the increment;
 *    step control of Sloan et al. (the sub-step grows by up to 1.1 after an accepted one, ESloan) or of Xie et al.
 *    (the rest of the increment after an accepted sub-step, EXie); model EConstantV (constant specific volume, as in
 *    the article) or EVoidRatio (void ratio updated, as in Xie et al.). Internal convention of the class TRKModel:
 *    compression positive, tensor components (11, 22, 33, 12, 23, 31).
 * 3. Material point drivers: UndrainedPath (\f$\Delta\varepsilon_{xx}=\Delta\varepsilon_{yy}=-\Delta\varepsilon_{zz}/2\f$,
 *    strain control) and DrainedPath (\f$\sigma_r=-p'_0\f$, \f$\varepsilon_{zz}\f$ prescribed: mixed control, the lateral
 *    strain of each increment is found by a Newton iteration with a finite-difference derivative and a bisection
 *    safeguard; the strain path is a straight line within the increment); closed forms: UndrainedClosed (eq. (B.7)
 *    and the elastic shear strain (B.9) with constant \f$\nu\f$) and mcc::TriaxialDrainedClosed (Appendix B.1).
 *
 * Tests (rivais() of gen_data.py; the work is the sum, over the increments, of the local Newton iterations (BE) or
 * of the evaluations of the elastoplastic operator in the Runge-Kutta stages (RK) of the converged updates):
 *  - test B of Xie et al. (2026) (RunXieTestB): undrained, normally consolidated, \f$M=0.896\f$, \f$\lambda=0.240\f$,
 *    \f$\kappa=0.045\f$, \f$v_0=2.27\f$ (\f$e_0=1.27\f$), \f$\nu=0.2\f$, \f$p'_0=p'_{c0}=689.02\f$ kPa, up to
 *    \f$\varepsilon_a=10\%\f$; relative error of the stress against the closed form at \f$\varepsilon_q=10\%\f$;
 *    BE with 1 to 1024 increments, RK in one increment with tolerances \f$10^{-1}\f$ to \f$10^{-8}\f$ (Sloan control)
 *    and the validation of the transcription with the control and model of Xie et al. (tolerances \f$10^{-1}\f$ to
 *    \f$10^{-5}\f$; Fig. 3 of Xie et al.);
 *  - example 1 of Krabbenhoft and Lyamin (2012) (RunKLUndrained, RunKLDrained):
 *    \f$M=3\sin\phi/\sqrt{3+\sin^2\phi}\f$ with \f$\phi=24^\circ\f$ (0.686), \f$\lambda=0.07\f$, \f$\kappa=0.008\f$,
 *    \f$v_0=1.5\f$, \f$\nu=0.3\f$, \f$p'_0=200\f$ kPa; undrained with OCR = 1 (to \f$\varepsilon_a=2.5\%\f$) and
 *    OCR = 10 (to 8%), errors of \f$p'\f$ and \f$q\f$, BE with 5 to 1000 increments and RK in 10 increments with
 *    tolerances \f$10^{-2}\f$ to \f$10^{-7}\f$; drained with OCR = 1 and 10 to \f$\varepsilon_a=25\%\f$, errors of
 *    \f$q\f$ and \f$\varepsilon_v\f$ and largest \f$q\f$, BE with 10 to 500 increments and RK (tolerance \f$10^{-4}\f$)
 *    with 10 to 100 increments; RK with tolerance \f$10^{-5}\f$ is also run (statement of Sect. 6.2, not in the
 *    pickled results).
 *
 * The numbers are compared with the results of the Python code (data_rivais.pkl, IntegrationSchemesReference.h):
 * the work counts are identical in all the 183 results and the values agree to round-off (largest differences about
 * \f$10^{-11}\f$ kPa), with one exception: the drained test with OCR = 10 and the secant shear modulus in 10
 * increments. In its first increment the Newton iterations of the drained driver pass through lateral strains where
 * the single secant step either fails or converges, after 20 to 50 local iterations, to spurious states (p' = 0 or
 * another branch), depending on round-off (RunSecantSweep), so that \f$\sigma_r(\varepsilon_r)+p'_0\f$ of the
 * sub-stepped update jumps between branches. Here the Newton iterations converge to an equilibrated state
 * (\f$|\sigma_r+p'_0|=3\times10^{-14}\f$ kPa at \f$\varepsilon_r=0.016479\f$); in the Python run they diverged, and
 * the bisection safeguard closed in on a jump of the function at \f$\varepsilon_r=0.0164014\f$ (from -2.56 kPa on
 * the sub-stepped branch to +202 kPa on a spurious p' = 0 branch), where it stopped after its 200 halvings and
 * accepted a state with \f$\sigma_r+p'_0=-2.56\f$ kPa, out of equilibrium. q at the end differs by 0.04 kPa and the
 * largest q by 2.2 kPa (that of the unbalanced state of the Python run). The driver records, for every drained run,
 * the increments solved by the bisection safeguard, those whose bisection did not converge and the largest
 * \f$|\sigma_r+p'_0|\f$ of the accepted increments (TPath): all the drained increments of this program are in
 * equilibrium to the tolerance \f$10^{-10}p'_0\f$.
 *
 * Two details of the transcription: the loading criterion of the Runge-Kutta schemes,
 * \f$a:\mathbb{D}:\Delta\varepsilon\f$, is evaluated as \f$a:(\lambda\,\mathrm{tr}\Delta\varepsilon\,m+2G\Delta\varepsilon)\f$,
 * which is exactly zero for the isochoric first increment from the isotropic state on the surface of the undrained
 * NC tests, as the matrix product of the Python code (the increment then starts with an elastic part of about
 * \f$3\times10^{-7}\f$ of the increment); and Brent's method is the transcription of scipy's brentq, so the
 * closed-form references coincide with those of the Python code to the last digit.
 *
 * There is no finite element mesh in this example: every test is a single material point; the methods follow
 * the sequence of the finite element examples (material set-up, solution, post-processing).
 */
class IntegrationSchemes {
public:
    /** @brief Symmetric tensor in the Voigt order of TPZTensor (XX, XY, XZ, YY, YZ, ZZ), tension positive; strains
     * with engineering shear components (convention of camclay_hw.py) */
    typedef std::array<REAL, 6> TVoigt;

    /** @brief Indices of the Voigt components */
    enum EVoigt { iXX = 0, iXY = 1, iXZ = 2, iYY = 3, iYZ = 4, iZZ = 5 };

    /** @brief Failure of a local projection (ProjectionError of camclay_hw.py) */
    class TProjectionError : public std::runtime_error {
    public:
        /** @brief Exception with the reason of the failure */
        explicit TProjectionError(const std::string &what) : std::runtime_error(what) {}
    };

    /** @brief Failure of a material point driver (RuntimeError of rivais.py) */
    class TDriverError : public std::runtime_error {
    public:
        /** @brief Exception with the reason of the failure */
        explicit TDriverError(const std::string &what) : std::runtime_error(what) {}
    };

    /** @brief Parameters of the Modified Cam-Clay model with porous elasticity (cc_parameters of camclay_hw.py) */
    struct TMaterial {
        REAL fM = 1.;       ///< slope of the critical state line
        REAL fLambda = 0.2; ///< slope of the normal compression line
        REAL fKappa = 0.05; ///< slope of the swelling line
        REAL fV0 = 2.;      ///< specific volume (constant)
        REAL fNu = 0.3;     ///< Poisson ratio of the shear modulus G = r K
        REAL fG = 0.;       ///< constant shear modulus if positive (option 'G'); otherwise G = r K
        REAL fPt = 0.;      ///< tensile strength
        REAL fOmega = 1.;   ///< shape parameter of the subcritical region (only 1 is supported)
        /** @brief Ratio \f$r=G/K=3(1-2\nu)/(2(1+\nu))\f$ (gfac) */
        REAL GFactor() const { return 3. * (1. - 2. * fNu) / (2. * (1. + fNu)); }
        /** @brief True with a constant shear modulus */
        bool ConstantG() const { return fG > 0.; }
        /** @brief Hardening law: a, \f$H=da/d\Delta\alpha\f$ and \f$p_c=p_{c,n}\exp(v_0\Delta\alpha/(\lambda-\kappa))\f$ */
        void Hardening(REAL pcn, REAL dal, REAL &a, REAL &H, REAL &pc) const {
            const REAL c = fV0 / (fLambda - fKappa);
            pc = pcn * std::exp(c * dal);
            a = (pc + fPt) / (1. + fOmega);
            H = c * pc / (1. + fOmega);
        }
        /** @brief Yield function \f$\Phi=(p-p_t+a)^2/b^2+\tfrac32\rho^2/M^2-a^2\f$ (p tension positive) */
        REAL Phi(REAL p, REAL rho, REAL a, REAL b) const {
            const REAL t = p - fPt + a;
            return (t * t) / (b * b) + 1.5 * (rho * rho) / (fM * fM) - a * a;
        }
    };

    /** @brief Volumetric law in the plastic correction of the backward-Euler schemes */
    enum EVolumetric {
        EExactK = 0, ///< exact integral of the porous law (this work)
        EFrozenK = 1 ///< bulk modulus frozen at its trial value (Sanei et al. 2020)
    };

    /** @brief Shear modulus in the plastic correction of the backward-Euler schemes */
    enum EShearModulus {
        EShearAtPn = 0,  ///< G = r K(p_n), kept in the step (this work)
        ESecantShear = 1 ///< secant G = r K_secant (Zhou 2022, Lu 2023, Krabbenhoft and Lyamin 2012, Bui 2026)
    };

    /** @brief Step size control of the Runge-Kutta schemes */
    enum EStepControl {
        ESloan = 0, ///< the sub-step grows by up to 1.1 after an accepted sub-step (Sloan et al. 2001)
        EXie = 1    ///< the next sub-step is the rest of the increment, eq. (36) of Xie et al. (2026)
    };

    /** @brief Specific volume in the Runge-Kutta model */
    enum ERKModel {
        EConstantV = 0, ///< constant specific volume (as in the article), model 'v0'
        EVoidRatio = 1  ///< void ratio updated with the volumetric strain (as in Xie et al.), model 'xie'
    };

    /** @brief Result of one stress update */
    struct TUpdate {
        TVoigt fSigma{};     ///< updated stress
        REAL fPc = 0.;       ///< updated preconsolidation pressure
        REAL fDal = 0.;      ///< plastic volumetric compaction \f$\Delta\alpha=-\Delta\varepsilon^p_v\f$ (BE)
        REAL fE = 0.;        ///< void ratio at the end of the increment (RK)
        int fWork = 0;       ///< local Newton iterations (BE) or evaluations of the elastoplastic operator (RK)
        int fSubsteps = 1;   ///< number of sub-steps (BackwardEulerSubstep)
        int fAttempts = 0;   ///< attempted sub-steps (RK)
        int fRejections = 0; ///< rejected sub-steps (RK)
    };

    /** @brief Stress update (sigma_n, strain increment, p_c,n) -> TUpdate; throws TProjectionError on failure */
    typedef std::function<TUpdate(const TVoigt &sign, const TVoigt &deps, REAL pcn)> TUpdateFunction;

    /** @brief Butcher tableau of an embedded explicit Runge-Kutta pair */
    struct TTableau {
        std::string fName;                ///< ME2(1) or RKDP5(4)
        std::vector<std::vector<REAL>> fA; ///< coefficients a_ij (lower triangular)
        std::vector<REAL> fB;             ///< weights of the solution that is propagated
        std::vector<REAL> fBhat;          ///< weights of the embedded solution (error estimate)
        int fOrder = 1;                   ///< order used in the step size control
    };

    /** @brief Path of a material point test */
    struct TPath {
        std::vector<std::array<REAL, 4>> fRows; ///< (eps_a, p', q, eps_v), compression positive
        long fWork = 0;                         ///< total work of the converged updates
        int fSubsteps = 0;                      ///< total number of sub-steps (= increments without sub-stepping)
        int fSubsteppedIncrements = 0;          ///< increments that needed sub-stepping
        long fAttempts = 0;                     ///< RK: attempted sub-steps of the converged updates
        long fRejections = 0;                   ///< RK: rejected sub-steps of the converged updates
        int fBisectionIncrements = 0;           ///< drained: increments solved by the bisection safeguard
        int fUnbalancedIncrements = 0;          ///< drained: bisections stopped without |sigma_r + p'0| < tol p'0
        REAL fMaxResidual = 0.;                 ///< drained: largest |sigma_r + p'0| of the accepted increments (kPa)
    };

    /**
     * @brief One run of Table 4 / Fig. 6: a scheme in a test with a number of increments (and a tolerance)
     *
     * The values fV follow data_rivais.pkl: test B of Xie et al.: (relative error of the stress, -, -); undrained
     * tests: (p' - p'_exact, q - q_exact, -) at the end; drained tests: (q - q_exact, eps_v - eps_v,exact, largest
     * q along the path) at eps_a = 25%.
     */
    struct TRecord {
        std::string fTest;     ///< xieB, kl_undrained_ocr1, kl_undrained_ocr10, kl_drained_ocr1, kl_drained_ocr10
        std::string fScheme;   ///< exact_n, exact_secant, frozen_n, library_exact, library_frozen, ME2(1), RKDP5(4)
        std::string fControl;  ///< RK step control (sloan, xie), "-" for BE
        std::string fModel;    ///< RK model (v0, xie), "-" for BE
        int fN = 0;            ///< number of increments
        REAL fTol = 0.;        ///< RK tolerance (0 for BE)
        bool fConverged = false; ///< false if the test has no solution
        REAL fV[3] = {std::numeric_limits<REAL>::quiet_NaN(), std::numeric_limits<REAL>::quiet_NaN(),
                      std::numeric_limits<REAL>::quiet_NaN()}; ///< values (see above)
        long fWork = -1;       ///< work (-1 if no solution)
        int fSubsteps = 0;     ///< BE: total number of sub-steps
        int fSubsteppedIncrements = 0; ///< BE: increments that needed sub-stepping
        long fAttempts = 0;    ///< RK: attempted sub-steps
        long fRejections = 0;  ///< RK: rejected sub-steps
        int fBisectionIncrements = 0;  ///< drained: increments solved by the bisection safeguard of the driver
        int fUnbalancedIncrements = 0; ///< drained: increments whose bisection stopped without equilibrium
        REAL fMaxResidual = std::numeric_limits<REAL>::quiet_NaN(); ///< drained: largest |sigma_r + p'0| (kPa)
        REAL fPEnd = std::numeric_limits<REAL>::quiet_NaN();   ///< p' at the end
        REAL fQEnd = std::numeric_limits<REAL>::quiet_NaN();   ///< q at the end
        REAL fEpsVEnd = std::numeric_limits<REAL>::quiet_NaN(); ///< eps_v at the end (drained tests)
        REAL fTime = 0.;       ///< computing time (s)
        std::string fMessage;  ///< reason of the failure
        TPath fPath;           ///< path of the test (kept for the cross-check of the library)
    };

    /** @brief Closed-form reference of a test */
    struct TReference {
        std::string fTest;   ///< name of the test
        REAL fEpsA = 0.;     ///< axial strain at which the errors are measured
        REAL fEta = std::numeric_limits<REAL>::quiet_NaN();   ///< stress ratio of the reference state (undrained)
        REAL fP = std::numeric_limits<REAL>::quiet_NaN();     ///< p' of the reference state (undrained)
        REAL fQ = 0.;        ///< q of the reference state
        REAL fEpsV = 0.;     ///< eps_v of the reference state (drained; 0 in the undrained tests)
        REAL fQPeak = std::numeric_limits<REAL>::quiet_NaN(); ///< largest q of the closed form (drained)
    };

    /** @name Materials of the tests */
    /** @{ */
    /** @brief Clay of test B of Xie et al. (2026) */
    static TMaterial XieMaterial();
    /** @brief Clay of example 1 of Krabbenhoft and Lyamin (2012) */
    static TMaterial KLMaterial();
    REAL fXieP0 = 689.02; ///< initial p' = p'_c of test B (kPa)
    REAL fXieE0 = 1.27;   ///< initial void ratio of test B
    REAL fKLP0 = 200.;    ///< initial p' of the Krabbenhoft and Lyamin tests (kPa)
    REAL fKLE0 = 0.5;     ///< initial void ratio of the Krabbenhoft and Lyamin tests
    /** @} */

    /** @name 1. Backward-Euler schemes (be_tensor, be_substep) */
    /** @{ */
    /**
     * @brief Backward-Euler update written in tensors with the radial return (be_tensor of rivais.py)
     * @param P material (omega = 1, porous law)
     * @param sign converged stress
     * @param deps strain increment (engineering shear components)
     * @param pcn converged preconsolidation pressure
     * @param volumetric exact or frozen porous law
     * @param shear shear modulus at p_n or secant
     * @param tol tolerance on the norm of the scaled residual
     * @param maxit largest number of Newton iterations
     * @return stress, p_c, Delta alpha and the number of local Newton iterations (0 for elastic steps)
     */
    static TUpdate BackwardEuler(const TMaterial &P, const TVoigt &sign, const TVoigt &deps, REAL pcn,
                                 EVolumetric volumetric, EShearModulus shear, REAL tol = 1e-12, int maxit = 50);

    /**
     * @brief Backward Euler with sub-stepping only when the single step fails (be_substep, strategy of Bui et al.
     * 2026): the increment is halved recursively until the projection converges (at most maxdepth levels)
     * @return the update with the sum of the iterations of the sub-steps and the number of sub-steps
     */
    static TUpdate BackwardEulerSubstep(const TMaterial &P, const TVoigt &sign, const TVoigt &deps, REAL pcn,
                                        EVolumetric volumetric, EShearModulus shear, int maxdepth = 12, int depth = 0);

    /** @brief Secant bulk modulus \f$\bar K=p_n(\exp(-cx)-1)/x\f$ and its derivative in x (series for small cx) */
    static void SecantBulkModulus(REAL pn, REAL c, REAL x, REAL &k, REAL &dk);
    /** @} */

    /** @name 2. Explicit Runge-Kutta schemes (RKModel, rk_update) */
    /** @{ */
    /** @brief Modified Euler pair ME2(1) */
    static TTableau ME21();
    /** @brief Runge-Kutta-Dormand-Prince pair RKDP5(4) of Sloan et al. (2001) */
    static TTableau RKDP54();

    /**
     * @brief Modified Cam-Clay of Xie et al. (2026), eqs. 7-11, in the internal convention of the Runge-Kutta
     * schemes (compression positive, tensor components 11, 22, 33, 12, 23, 31; contraction with the weights
     * (1, 1, 1, 2, 2, 2)); state X = (sigma, p_x = p_c, e)
     */
    class TRKModel {
    public:
        /** @brief State: six stress components, preconsolidation pressure and void ratio */
        typedef std::array<REAL, 8> TState;
        /** @brief Vector of six internal components */
        typedef std::array<REAL, 6> TVec6;
        /**
         * @brief Model of a material
         * @param P material
         * @param model constant specific volume or updated void ratio
         * @param tolf tolerance on the yield function (intersection and drift correction)
         */
        TRKModel(const TMaterial &P, ERKModel model, REAL tolf)
            : fM2(P.fM * P.fM), fLambda(P.fLambda), fKappa(P.fKappa), fR(P.GFactor()), fConstantG(P.ConstantG()),
              fGc(P.fG), fModel(model), fTolF(tolf) {}
        /** @brief Specific volume v = 1 + e */
        REAL V(const TState &X) const { return 1. + X[7]; }
        /** @brief Mean stress p, deviatoric stress q and deviator s of the stress of X */
        void Invariants(const TState &X, REAL &p, REAL &q, TVec6 &s) const;
        /** @brief Elastic matrix at the mean stress p (K = v p / kappa, G = r K or constant) */
        void Elastic(REAL p, const TState &X, REAL D[6][6]) const;
        /** @brief Yield function \f$f=M^2p(p-p_x)+q^2\f$ */
        REAL F(const TState &X) const;
        /** @brief Mean stress, gradient a = df/dsigma and \f$c_p=(\lambda-\kappa)/v\f$ */
        void Grads(const TState &X, REAL &p, TVec6 &a, REAL &cp) const;
        /** @brief Elastoplastic operator (8 x 6) of the state rates; counts one evaluation */
        void Cep(const TState &X, REAL C[8][6]);
        /** @brief Exact elastic update with the fraction beta of the strain increment de (eqs. 23, 24 and 28) */
        TState ElasticUpdate(const TState &X, REAL beta, const TVec6 &de) const;
        /** @brief Derivative of the yield function along the elastic path, df/dbeta */
        REAL DfDbeta(const TState &X, const TVec6 &de) const;
        /** @brief Drift correction back to the yield surface (at most maxit corrections) */
        TState Drift(TState X, int maxit = 50) const;
        /** @brief One sub-step with the pair: propagated and embedded solutions */
        void Stages(const TState &X, const TVec6 &de, const TTableau &tab, TState &Xs, TState &Xh);
        /** @brief Number of evaluations of the elastoplastic operator */
        int NEval() const { return fNEval; }
        /** @brief Tolerance on the yield function */
        REAL TolF() const { return fTolF; }

    private:
        REAL fM2;        ///< M^2
        REAL fLambda;    ///< lambda
        REAL fKappa;     ///< kappa
        REAL fR;         ///< G / K
        bool fConstantG; ///< constant shear modulus
        REAL fGc;        ///< constant shear modulus
        ERKModel fModel; ///< model of the specific volume
        REAL fTolF;      ///< tolerance on the yield function
        int fNEval = 0;  ///< evaluations of the elastoplastic operator
    };

    /**
     * @brief One strain increment with an adaptive explicit Runge-Kutta scheme (rk_update of rivais.py)
     * @param P material
     * @param sign converged stress
     * @param deps strain increment (engineering shear components)
     * @param pcn converged preconsolidation pressure
     * @param e0 void ratio at the start of the increment
     * @param tab embedded pair
     * @param tol tolerance on the relative local error
     * @param control step size control
     * @param model model of the specific volume
     * @param tolf tolerance on the yield function
     * @param maxatt largest number of attempted sub-steps
     * @return stress, p_c, void ratio, evaluations of the elastoplastic operator, attempts and rejections
     */
    static TUpdate RungeKutta(const TMaterial &P, const TVoigt &sign, const TVoigt &deps, REAL pcn, REAL e0,
                              const TTableau &tab, REAL tol, EStepControl control, ERKModel model, REAL tolf = 1e-5,
                              int maxatt = 200000);
    /** @} */

    /** @name 3. Material point drivers and closed forms */
    /** @{ */
    /**
     * @brief Undrained test with eps_v = 0 imposed (undrained_path): n equal increments
     * \f$\Delta\varepsilon_{xx}=\Delta\varepsilon_{yy}=-\Delta\varepsilon_{zz}/2=\varepsilon_{a,max}/(2n)\f$ (strain control)
     */
    static TPath UndrainedPath(const TUpdateFunction &update, REAL p0, REAL pc0, REAL eamax, int n);

    /**
     * @brief Drained test with sigma_r = -p0 constant and eps_zz prescribed (drained_path): mixed control; Newton
     * iteration on the lateral strain with a finite-difference derivative and, if it fails, bisection in an
     * interval where sigma_r + p0 changes sign. The work is that of the converged update of each increment.
     * Throws TProjectionError or TDriverError when the test has no solution.
     */
    static TPath DrainedPath(const TUpdateFunction &update, REAL p0, REAL pc0, REAL eamax, int n, REAL tol = 1e-10);

    /**
     * @brief Closed-form undrained path with constant Poisson ratio (undrained_closed; Appendix B.2, eq. (B.7)
     * for p'(eta) and eq. (B.9) for the elastic shear strain)
     * @param p0 initial p'
     * @param pc0 initial p'_c
     * @param P material
     * @param eta stress ratio q/p'
     * @param[out] p p'
     * @param[out] q q
     * @param[out] epsq shear strain \f$\varepsilon_q=\varepsilon^p_q+\varepsilon^e_q\f$ (= eps_a when eps_v = 0)
     */
    static void UndrainedClosed(REAL p0, REAL pc0, const TMaterial &P, REAL eta, REAL &p, REAL &q, REAL &epsq);

    /** @brief Brent's method of scipy.optimize.brentq (same iterations; xtol, rtol = 4 eps, 100 iterations) */
    static REAL Brent(const std::function<REAL(REAL)> &f, REAL xa, REAL xb, REAL xtol,
                      REAL rtol = 4. * std::numeric_limits<REAL>::epsilon(), int maxiter = 100);

    /** @brief p' and q of a stress (pq of camclay_hw.py), compression positive */
    static void PQ(const TVoigt &sig, REAL &p, REAL &q);

    /** @brief Axisymmetric compression stress (tension positive) with the given p' and q */
    static TVoigt Axisymmetric(REAL p, REAL q);
    /** @} */

    /** @name Stress updates used by the drivers */
    /** @{ */
    /** @brief Backward Euler (with sub-stepping when substep is true) */
    static TUpdateFunction BEUpdate(const TMaterial &P, EVolumetric volumetric, EShearModulus shear, bool substep);
    /** @brief Runge-Kutta, with the void ratio e0 at the start of every increment (as in gen_data.py) */
    static TUpdateFunction RKUpdate(const TMaterial &P, REAL e0, const TTableau &tab, REAL tol,
                                    EStepControl control = ESloan, ERKModel model = EConstantV);
    /**
     * @brief The library class TPZPlasticStepModifiedCamClay (mcc::ApplyStrain from the converged state with
     * eps_n = 0 and eps = deps): the return mapping of this work in spectral form; the work is
     * LastNewtonIterations()
     */
    static TUpdateFunction LibraryUpdate(const TMaterial &P, TPZYCModifiedCamClayRHW::EPorousIntegration integration);
    /** @} */

    /** @name Tests */
    /** @{ */
    /** @brief Test B of Xie et al. (2026): BE with 1 to 1024 increments, RK in one increment (Fig. 6a, Table 4a) */
    void RunXieTestB();
    /** @brief Undrained test of Krabbenhoft and Lyamin, OCR = 1 (eps_a = 2.5%) or 10 (8%) (Fig. 6b, Table 4b) */
    void RunKLUndrained(int ocr, REAL eamax);
    /** @brief Drained test of Krabbenhoft and Lyamin, OCR = 1 or 10, eps_a = 25% (Fig. 6c, Table 4c) */
    void RunKLDrained(int ocr);
    /**
     * @brief Robustness of the single backward-Euler step in the first increment of the drained test with OCR = 10
     * and 10 increments (Delta eps_a = 2.5%): sweep of the lateral strain eps_r (2001 points in [0.010, 0.030]) with
     * this work and with the secant shear modulus without sub-stepping; counts the failures and the spurious
     * converged states (p' = 0)
     */
    void RunSecantSweep();
    /**
     * @brief Error of the return mapping of this work against the tolerance of its local Newton iterations
     * (test B of Xie et al. with 1, 4, 16 and 64 increments, tolerances 1e-2 to 1e-14 on the scaled residual):
     * the error of a backward-Euler step is its truncation error, which a tighter tolerance does not reduce
     */
    void RunToleranceStudy();
    /** @brief All the tests, the post-processing and the comparisons */
    void RunAll();
    /** @} */

    /** @name Post-processing */
    /** @{ */
    /** @brief Writes the CSV files (results, Table 4, Fig. 6, references, library check, Python comparison) */
    void PostProcess() const;
    /** @brief Prints Table 4 with the values of the article (v0.6) */
    void PrintTable4() const;
    /** @brief Prints the numbers quoted in the text of Sect. 6.2 */
    void PrintTextNumbers() const;
    /** @brief Prints the comparison of the library class with be_tensor (exact and frozen) */
    void PrintLibraryCheck() const;
    /** @brief Prints the largest differences with data_rivais.pkl */
    void PrintPythonComparison() const;
    /** @} */

    /** @brief Results of all the runs */
    const std::vector<TRecord> &Records() const { return fRecords; }

    /** @brief Finds a run; nullptr if absent */
    const TRecord *Find(const std::string &test, const std::string &scheme, int n, REAL tol = 0.,
                        const std::string &control = "") const;

private:
    /** @brief Runs a test with an update function, catching the failures; stores and returns the record */
    TRecord &Store(TRecord rec);
    /** @brief Runs one drained or undrained path and fills the record */
    void RunPath(TRecord &rec, const std::function<TPath()> &path, const TReference &ref);
    /** @brief Reference of a test */
    const TReference &Reference(const std::string &test) const;
    /** @brief Number formatted with 17 significant digits ("nan" for NaN) */
    static std::string Num(REAL v);
    /** @brief Writes a CSV file with string cells */
    static void WriteTable(const std::string &file, const std::vector<std::string> &header,
                           const std::vector<std::vector<std::string>> &rows);
    /** @brief The value of the error quantity drawn in Fig. 6 for a record */
    static REAL FigureError(const TRecord &rec);
    /** @brief Solution of a 4 x 4 system by Gaussian elimination with partial pivoting (numpy.linalg.solve) */
    static bool Solve4(std::array<std::array<REAL, 4>, 4> A, std::array<REAL, 4> &b);

    /** @brief One point of RunSecantSweep */
    struct TSweepPoint {
        REAL fEpsR = 0.;      ///< lateral strain of the increment
        int fScheme = 0;      ///< 0 this work (exact K, G(p_n)), 1 secant G (single step)
        bool fConverged = false; ///< the local Newton iterations converged
        int fIterations = 0;  ///< local Newton iterations of a converged step (0 if the step failed after 50)
        REAL fP = 0.;         ///< p' of the converged state (NaN in the CSV file if the step failed)
        REAL fQ = 0.;         ///< q of the converged state (NaN in the CSV file if the step failed)
        REAL fSigmaR = 0.;    ///< sigma_r + p'_0 of the converged state (zero at the solution of the increment)
    };
    std::vector<TSweepPoint> fSweep;     ///< results of RunSecantSweep
    /** @brief One run of RunToleranceStudy: (n, tolerance, converged, relative error, local iterations) */
    struct TTolerancePoint {
        int fN = 0;              ///< number of increments
        REAL fTol = 0.;          ///< tolerance of the local Newton iterations
        bool fConverged = false; ///< the test has a solution
        REAL fError = 0.;        ///< relative error of the stress at eps_a = 10%
        long fWork = 0;          ///< local Newton iterations
    };
    std::vector<TTolerancePoint> fTolerance; ///< results of RunToleranceStudy
    std::vector<TRecord> fRecords;       ///< all the runs
    std::vector<TReference> fReferences; ///< closed-form references of the tests
    REAL fTotalTime = 0.;                ///< total computing time (s)
};

// ============================================================================================ materials and helpers

inline IntegrationSchemes::TMaterial IntegrationSchemes::XieMaterial() {
    TMaterial P;
    P.fM = 0.896;
    P.fLambda = 0.240;
    P.fKappa = 0.045;
    P.fV0 = 2.27;
    P.fNu = 0.2;
    return P;
}

inline IntegrationSchemes::TMaterial IntegrationSchemes::KLMaterial() {
    // M = 3 sin(phi) / sqrt(3 + sin^2(phi)), phi = 24 degrees (np.radians(24) = 24 * (pi / 180))
    const REAL s = std::sin(24. * (M_PI / 180.));
    TMaterial P;
    P.fM = 3. * s / std::sqrt(3. + s * s);
    P.fLambda = 0.07;
    P.fKappa = 0.008;
    P.fV0 = 1.5;
    P.fNu = 0.3;
    return P;
}

inline void IntegrationSchemes::PQ(const TVoigt &sig, REAL &p, REAL &q) {
    const REAL pm = (sig[iXX] + sig[iYY] + sig[iZZ]) / 3.;
    // s = voigt_to_cart(sig) - p I, (s * s).sum() over the 9 entries (numpy pairwise sum of 9 values)
    const REAL c[9] = {sig[iXX] - pm, sig[iXY], sig[iXZ], sig[iXY], sig[iYY] - pm, sig[iYZ],
                       sig[iXZ], sig[iYZ], sig[iZZ] - pm};
    REAL r[9];
    for (int i = 0; i < 9; ++i) r[i] = c[i] * c[i];
    const REAL sum = (((r[0] + r[1]) + (r[2] + r[3])) + ((r[4] + r[5]) + (r[6] + r[7]))) + r[8];
    p = -pm;
    q = std::sqrt(1.5 * sum);
}

inline IntegrationSchemes::TVoigt IntegrationSchemes::Axisymmetric(REAL p, REAL q) {
    TVoigt s{};
    s[iZZ] = -(p + 2. * q / 3.);
    s[iXX] = -(p - q / 3.);
    s[iYY] = -(p - q / 3.);
    return s;
}

inline bool IntegrationSchemes::Solve4(std::array<std::array<REAL, 4>, 4> A, std::array<REAL, 4> &b) {
    for (int j = 0; j < 4; ++j) {
        int ip = j;
        REAL amax = std::fabs(A[j][j]);
        for (int i = j + 1; i < 4; ++i)
            if (std::fabs(A[i][j]) > amax) {
                amax = std::fabs(A[i][j]);
                ip = i;
            }
        if (A[ip][j] == 0.) return false; // exactly singular (LinAlgError of numpy)
        if (ip != j) {
            std::swap(A[ip], A[j]);
            std::swap(b[ip], b[j]);
        }
        const REAL rp = 1. / A[j][j];
        for (int i = j + 1; i < 4; ++i) {
            A[i][j] *= rp;
            for (int k = j + 1; k < 4; ++k) A[i][k] -= A[i][j] * A[j][k];
            b[i] -= A[i][j] * b[j];
        }
    }
    for (int i = 3; i >= 0; --i) {
        REAL s = b[i];
        for (int k = i + 1; k < 4; ++k) s -= A[i][k] * b[k];
        b[i] = s / A[i][i];
    }
    return true;
}

// ============================================================================================ 1. backward Euler

inline void IntegrationSchemes::SecantBulkModulus(REAL pn, REAL c, REAL x, REAL &k, REAL &dk) {
    if (std::fabs(c * x) < 1e-4) {
        const REAL y = c * x;
        k = -c * pn * (1. - y / 2. + (y * y) / 6. - std::pow(y, 3.) / 24. + std::pow(y, 4.) / 120.);
        dk = -c * pn * (-c / 2. + c * y / 3. - c * (y * y) / 8. + c * std::pow(y, 3.) / 30.);
        return;
    }
    const REAL e = std::exp(-c * x);
    k = pn * (e - 1.) / x;
    dk = pn * (-c * x * e - (e - 1.)) / (x * x);
}

inline IntegrationSchemes::TUpdate IntegrationSchemes::BackwardEuler(const TMaterial &P, const TVoigt &sign,
                                                                     const TVoigt &deps, REAL pcn,
                                                                     EVolumetric volumetric, EShearModulus shear,
                                                                     REAL tol, int maxit) {
    if (P.fOmega != 1.) throw std::invalid_argument("BackwardEuler: only omega = 1 is implemented");
    static const REAL MV[6] = {1., 0., 0., 1., 0., 1.};
    static const REAL SHF[6] = {1., 0.5, 0.5, 1., 0.5, 1.}; // engineering -> tensor shear components
    static const REAL WT[6] = {1., 2., 2., 1., 2., 1.};     // contraction of tensors stored in Voigt form
    const REAL sq3 = std::sqrt(3.);
    const REAL c = P.fV0 / P.fKappa;
    const REAL M2 = P.fM * P.fM;
    const REAL pn = (sign[iXX] + sign[iYY] + sign[iZZ]) / 3.;
    const REAL dev = deps[iXX] + deps[iYY] + deps[iZZ];
    TVoigt sn, de;
    for (int i = 0; i < 6; ++i) {
        sn[i] = sign[i] - pn * MV[i];
        de[i] = (deps[i] - dev / 3. * MV[i]) * SHF[i];
    }
    const REAL Kn = -c * pn;
    const REAL Gn = P.GFactor() * Kn;
    // shear modulus G(x) of the correction and dG/dx, x = Delta eps_v + Delta alpha
    auto gfun = [&](REAL x, REAL &G, REAL &dG) {
        if (P.ConstantG()) {
            G = P.fG;
            dG = 0.;
        } else if (shear == EShearAtPn) {
            G = Gn;
            dG = 0.;
        } else {
            REAL k, dk;
            SecantBulkModulus(pn, c, x, k, dk);
            G = P.GFactor() * k;
            dG = P.GFactor() * dk;
        }
    };
    auto dott = [](const TVoigt &a, const TVoigt &b) {
        REAL s = 0.;
        for (int i = 0; i < 6; ++i) s += WT[i] * a[i] * b[i];
        return s;
    };
    const REAL A = dott(sn, sn), B = dott(sn, de), C = dott(de, de);
    // trial deviatoric norm rho_tr(G) = ||s_n + 2 G de|| and its derivative in G
    auto rhotr = [&](REAL G, REAL &r, REAL &dr) {
        r = std::sqrt(std::max(A + 4. * G * B + 4. * G * G * C, 0.));
        dr = r > 0. ? (2. * B + 4. * G * C) / r : 0.;
    };
    const REAL ptr = pn * std::exp(-c * dev);
    const REAL xitr = sq3 * ptr;
    REAL G0, dG0, rt0, drt0;
    gfun(dev, G0, dG0);
    rhotr(G0, rt0, drt0);
    REAL an, Hn, pcdum;
    P.Hardening(pcn, 0., an, Hn, pcdum);
    TUpdate res;
    if (P.Phi(ptr, rt0, an, 1.) <= 1e-11 * (an * an)) { // elastic step
        for (int i = 0; i < 6; ++i) res.fSigma[i] = ptr * MV[i] + sn[i] + 2. * G0 * de[i];
        res.fPc = pcn;
        res.fWork = 0;
        return res;
    }
    std::array<REAL, 4> X = {xitr, rt0, 0., 0.};
    int it = 0;
    for (it = 0; it <= maxit; ++it) {
        const REAL xi = X[0], rho = X[1], dal = X[2], dg = X[3];
        REAL a, H, pc;
        P.Hardening(pcn, dal, a, H, pc);
        const REAL pbar = xi / sq3 - P.fPt + a;
        REAL G, dG, rt, drt;
        gfun(dev + dal, G, dG);
        rhotr(G, rt, drt);
        const REAL e1 = volumetric == EExactK ? std::exp(-c * dal) : 1. - c * dal;
        std::array<REAL, 4> R = {(xi - xitr * e1) / an, (rho * (1. + 6. * G * dg / M2) - rt) / an,
                                 dal + 2. * dg * pbar,
                                 ((pbar * pbar) + 1.5 * (rho * rho) / M2 - (a * a)) / (an * an)};
        const REAL norm = std::sqrt(R[0] * R[0] + R[1] * R[1] + R[2] * R[2] + R[3] * R[3]);
        if (norm < tol) break; // a NaN residual never converges: the iterations go on until maxit, as in Python
        if (it == maxit) throw TProjectionError("backward Euler (tensors) did not converge");
        const REAL dexp = volumetric == EExactK ? std::exp(-c * dal) : 1.;
        std::array<std::array<REAL, 4>, 4> J;
        J[0] = {1. / an, 0. / an, c * xitr * dexp / an, 0. / an};
        J[1] = {0. / an, (1. + 6. * G * dg / M2) / an, ((6. * rho * dg / M2 - drt) * dG) / an, (6. * G * rho / M2) / an};
        J[2] = {2. * dg / sq3, 0., 1. + 2. * dg * H, 2. * pbar};
        const REAL an2 = an * an;
        J[3] = {(2. * pbar / sq3) / an2, (3. * rho / M2) / an2, (2. * pbar * H - 2. * a * H) / an2, 0. / an2};
        if (!Solve4(J, R)) throw TProjectionError("backward Euler (tensors): singular Jacobian");
        for (int i = 0; i < 4; ++i) X[i] -= R[i];
    }
    const REAL xi = X[0], rho = X[1], dal = X[2];
    REAL G, dG, rt, drt;
    gfun(dev + dal, G, dG);
    rhotr(G, rt, drt);
    for (int i = 0; i < 6; ++i) {
        const REAL s = rt > 0. ? (rho / rt) * (sn[i] + 2. * G * de[i]) : 0.;
        res.fSigma[i] = xi / sq3 * MV[i] + s;
    }
    REAL a, H, pc;
    P.Hardening(pcn, dal, a, H, pc);
    res.fPc = pc;
    res.fDal = dal;
    res.fWork = it;
    return res;
}

inline IntegrationSchemes::TUpdate IntegrationSchemes::BackwardEulerSubstep(const TMaterial &P, const TVoigt &sign,
                                                                            const TVoigt &deps, REAL pcn,
                                                                            EVolumetric volumetric,
                                                                            EShearModulus shear, int maxdepth,
                                                                            int depth) {
    try {
        TUpdate u = BackwardEuler(P, sign, deps, pcn, volumetric, shear);
        u.fSubsteps = 1;
        return u;
    } catch (TProjectionError &) {
        if (depth == maxdepth) throw;
    }
    TVoigt half;
    for (int i = 0; i < 6; ++i) half[i] = deps[i] / 2.;
    const TUpdate u1 = BackwardEulerSubstep(P, sign, half, pcn, volumetric, shear, maxdepth, depth + 1);
    TUpdate u2 = BackwardEulerSubstep(P, u1.fSigma, half, u1.fPc, volumetric, shear, maxdepth, depth + 1);
    u2.fWork += u1.fWork;
    u2.fSubsteps += u1.fSubsteps;
    return u2;
}

// ============================================================================================ 2. Runge-Kutta

inline IntegrationSchemes::TTableau IntegrationSchemes::ME21() {
    TTableau t;
    t.fName = "ME2(1)";
    t.fA = {{0., 0.}, {1., 0.}};
    t.fB = {0.5, 0.5};
    t.fBhat = {1.0, 0.0};
    t.fOrder = 2;
    return t;
}

inline IntegrationSchemes::TTableau IntegrationSchemes::RKDP54() {
    TTableau t;
    t.fName = "RKDP5(4)";
    t.fA.assign(6, std::vector<REAL>(6, 0.));
    t.fA[1][0] = 1. / 5.;
    t.fA[2][0] = 3. / 40.;
    t.fA[2][1] = 9. / 40.;
    t.fA[3][0] = 3. / 10.;
    t.fA[3][1] = -9. / 10.;
    t.fA[3][2] = 6. / 5.;
    t.fA[4][0] = 226. / 729.;
    t.fA[4][1] = -25. / 27.;
    t.fA[4][2] = 880. / 729.;
    t.fA[4][3] = 55. / 729.;
    t.fA[5][0] = -181. / 270.;
    t.fA[5][1] = 5. / 2.;
    t.fA[5][2] = -266. / 297.;
    t.fA[5][3] = -91. / 27.;
    t.fA[5][4] = 189. / 55.;
    t.fB = {19. / 216., 0., 1000. / 2079., -125. / 216., 81. / 88., 5. / 56.};
    t.fBhat = {31. / 540., 0., 190. / 297., -145. / 108., 351. / 220., 1. / 20.};
    t.fOrder = 5;
    return t;
}

namespace mccschemes {
/** @brief m = (1, 1, 1, 0, 0, 0) in the internal convention of the Runge-Kutta schemes */
static const REAL kM[6] = {1., 1., 1., 0., 0., 0.};
/** @brief Weights of the contraction of tensors stored as (11, 22, 33, 12, 23, 31) */
static const REAL kW[6] = {1., 1., 1., 2., 2., 2.};
/** @brief Voigt index (XX, XY, XZ, YY, YZ, ZZ) of each internal component (11, 22, 33, 12, 23, 31) */
static const int kIdx[6] = {0, 3, 5, 1, 4, 2};
} // namespace mccschemes

inline void IntegrationSchemes::TRKModel::Invariants(const TState &X, REAL &p, REAL &q, TVec6 &s) const {
    using namespace mccschemes;
    p = (X[0] + X[1] + X[2]) / 3.;
    REAL sum = 0.;
    for (int i = 0; i < 6; ++i) {
        s[i] = X[i] - p * kM[i];
        sum += kW[i] * s[i] * s[i];
    }
    q = std::sqrt(1.5 * sum);
}

inline void IntegrationSchemes::TRKModel::Elastic(REAL p, const TState &X, REAL D[6][6]) const {
    using namespace mccschemes;
    const REAL K = V(X) * p / fKappa;
    const REAL G = fConstantG ? fGc : fR * K;
    const REAL lam = K - 2. * G / 3.;
    for (int i = 0; i < 6; ++i)
        for (int j = 0; j < 6; ++j) D[i][j] = lam * (kM[i] * kM[j]) + (i == j ? (i < 3 ? 2. * G : G) : 0.);
}

inline REAL IntegrationSchemes::TRKModel::F(const TState &X) const {
    REAL p, q;
    TVec6 s;
    Invariants(X, p, q, s);
    return fM2 * p * (p - X[6]) + q * q;
}

inline void IntegrationSchemes::TRKModel::Grads(const TState &X, REAL &p, TVec6 &a, REAL &cp) const {
    using namespace mccschemes;
    REAL q;
    TVec6 s;
    Invariants(X, p, q, s);
    const REAL px = X[6];
    for (int i = 0; i < 6; ++i) a[i] = fM2 * (2. * p - px) / 3. * kM[i] + 3. * s[i];
    cp = (fLambda - fKappa) / V(X);
}

inline void IntegrationSchemes::TRKModel::Cep(const TState &X, REAL C[8][6]) {
    using namespace mccschemes;
    fNEval++;
    REAL p, cp;
    TVec6 a;
    Grads(X, p, a, cp);
    const REAL px = X[6];
    REAL D[6][6];
    Elastic(p, X, D);
    TVec6 Da, aD;
    for (int i = 0; i < 6; ++i) {
        REAL s = 0.;
        for (int j = 0; j < 6; ++j) s += D[i][j] * (kW[j] * a[j]);
        Da[i] = s;
    }
    for (int j = 0; j < 6; ++j) {
        REAL s = 0.;
        for (int i = 0; i < 6; ++i) s += (a[i] * kW[i]) * D[i][j];
        aD[j] = s;
    }
    const REAL h = (px / cp) * fM2 * (2. * p - px); // dpx/deps_v^p * dg/dp
    REAL den = 0.;
    for (int i = 0; i < 6; ++i) den += a[i] * (kW[i] * Da[i]);
    den += fM2 * p * h;
    for (int i = 0; i < 6; ++i)
        for (int j = 0; j < 6; ++j) C[i][j] = D[i][j] - Da[i] * aD[j] / den;
    for (int j = 0; j < 6; ++j) {
        C[6][j] = h * aD[j] / den;
        C[7][j] = fModel == EVoidRatio ? -V(X) * kM[j] : 0.;
    }
}

inline IntegrationSchemes::TRKModel::TState IntegrationSchemes::TRKModel::ElasticUpdate(const TState &X, REAL beta,
                                                                                        const TVec6 &de) const {
    using namespace mccschemes;
    REAL p, q;
    TVec6 s;
    Invariants(X, p, q, s);
    const REAL px = X[6], e = X[7];
    const REAL dev = de[0] + de[1] + de[2];
    const REAL cK = V(X) / fKappa;
    const REAL pnew = p * std::exp(cK * beta * dev);
    TState out;
    for (int i = 0; i < 6; ++i) {
        const REAL dee = de[i] - dev / 3. * kM[i];
        REAL snew;
        if (fConstantG) snew = s[i] + 2. * fGc * beta * dee;
        else if (std::fabs(dev) > 1e-14) snew = s[i] + 2. * fR * (pnew - p) / dev * dee;
        else snew = s[i] + 2. * fR * cK * p * beta * dee;
        out[i] = snew + pnew * kM[i];
    }
    out[6] = px;
    out[7] = fModel == EVoidRatio ? e - (1. + e) * (1. - std::exp(-beta * dev)) : e;
    return out;
}

inline REAL IntegrationSchemes::TRKModel::DfDbeta(const TState &X, const TVec6 &de) const {
    using namespace mccschemes;
    REAL p, cp, pp, q;
    TVec6 a, s;
    Grads(X, p, a, cp);
    Invariants(X, pp, q, s);
    const REAL dev = de[0] + de[1] + de[2];
    const REAL cK = V(X) / fKappa;
    const REAL G2 = 2. * (fConstantG ? fGc : fR * cK * p);
    REAL sum = 0.;
    for (int i = 0; i < 6; ++i) {
        const REAL dee = de[i] - dev / 3. * kM[i];
        sum += kW[i] * s[i] * (G2 * dee);
    }
    return fM2 * (2. * p - X[6]) * cK * p * dev + 3. * sum;
}

inline IntegrationSchemes::TRKModel::TState IntegrationSchemes::TRKModel::Drift(TState X, int maxit) const {
    using namespace mccschemes;
    for (int k = 0; k < maxit; ++k) {
        const REAL fv = F(X);
        if (std::fabs(fv) <= fTolF) return X;
        REAL p, cp;
        TVec6 a;
        Grads(X, p, a, cp);
        const REAL px = X[6];
        REAL D[6][6];
        Elastic(p, X, D);
        TState R;
        for (int i = 0; i < 6; ++i) {
            REAL s = 0.;
            for (int j = 0; j < 6; ++j) s += D[i][j] * (kW[j] * a[j]);
            R[i] = s;
        }
        R[6] = -(px / cp) * fM2 * (2. * p - px);
        R[7] = 0.;
        REAL den = 0.;
        for (int i = 0; i < 6; ++i) den += a[i] * kW[i] * R[i];
        den += -fM2 * p * R[6];
        const REAL fac = fv / den;
        for (int i = 0; i < 8; ++i) X[i] -= fac * R[i];
    }
    return X;
}

inline void IntegrationSchemes::TRKModel::Stages(const TState &X, const TVec6 &de, const TTableau &tab, TState &Xs,
                                                 TState &Xh) {
    using namespace mccschemes;
    const int ns = (int)tab.fB.size();
    std::vector<TState> ks(ns);
    TVec6 wde;
    for (int j = 0; j < 6; ++j) wde[j] = kW[j] * de[j];
    REAL C[8][6];
    for (int i = 0; i < ns; ++i) {
        TState acc{};
        for (int j = 0; j < i; ++j)
            for (int l = 0; l < 8; ++l) acc[l] += tab.fA[i][j] * ks[j][l];
        TState Xi;
        for (int l = 0; l < 8; ++l) Xi[l] = X[l] + acc[l];
        Cep(Xi, C);
        for (int l = 0; l < 8; ++l) {
            REAL s = 0.;
            for (int j = 0; j < 6; ++j) s += C[l][j] * wde[j];
            ks[i][l] = s;
        }
    }
    TState sb{}, sh{};
    for (int i = 0; i < ns; ++i)
        for (int l = 0; l < 8; ++l) {
            sb[l] += tab.fB[i] * ks[i][l];
            sh[l] += tab.fBhat[i] * ks[i][l];
        }
    for (int l = 0; l < 8; ++l) {
        Xs[l] = X[l] + sb[l];
        Xh[l] = X[l] + sh[l];
    }
}

inline IntegrationSchemes::TUpdate IntegrationSchemes::RungeKutta(const TMaterial &P, const TVoigt &sign,
                                                                  const TVoigt &deps, REAL pcn, REAL e0,
                                                                  const TTableau &tab, REAL tol,
                                                                  EStepControl control, ERKModel model, REAL tolf,
                                                                  int maxatt) {
    using namespace mccschemes;
    typedef TRKModel::TState TState;
    TRKModel mod(P, model, tolf);
    TState X;
    TRKModel::TVec6 de;
    for (int k = 0; k < 6; ++k) {
        X[k] = -sign[kIdx[k]];
        de[k] = -deps[kIdx[k]] * (k >= 3 ? 0.5 : 1.);
    }
    X[6] = pcn;
    X[7] = e0;
    auto external = [](const TState &Y, TUpdate &u) {
        for (int k = 0; k < 6; ++k) u.fSigma[kIdx[k]] = -Y[k];
        u.fPc = Y[6];
        u.fE = Y[7];
    };
    TUpdate res;
    const TState Xtr = mod.ElasticUpdate(X, 1., de);
    if (mod.F(Xtr) < tolf) { // elastic increment
        external(Xtr, res);
        res.fPc = pcn;
        res.fWork = 0;
        return res;
    }
    const REAL fn = mod.F(X);
    {
        // loading criterion a : D : de, with D : de = (K - 2G/3) tr(de) m + 2G de in closed form. For an isochoric
        // increment from an isotropic state on the surface (test B of Xie et al., undrained NC tests) the criterion is
        // zero and the increment is treated as not loading (an elastic part of about 3e-7 of the increment is found
        // before the surface is reached), as in the Python code, where the matrix product gives exactly zero there.
        REAL p, cp;
        TRKModel::TVec6 a;
        mod.Grads(X, p, a, cp);
        const REAL K = mod.V(X) * p / P.fKappa;
        const REAL G = P.ConstantG() ? P.fG : P.GFactor() * K;
        const REAL dev = de[0] + de[1] + de[2];
        REAL load = 0.;
        for (int i = 0; i < 6; ++i) load += a[i] * kW[i] * ((K - 2. * G / 3.) * dev * kM[i] + 2. * G * de[i]);
        const bool loading = load > 0.;
        REAL beta = 0.;
        if (fn < -tolf || !loading) {
            // elastic part up to the yield surface (Newton in beta with a bisection safeguard)
            REAL lo = 0., hi = 1.;
            if (fn >= -tolf) { // on the surface, unloading first and then reloading: bracket by sampling
                std::vector<REAL> bs(100), fs(100);
                for (int i = 1; i <= 100; ++i) bs[i - 1] = i == 100 ? 1. : i * (1. / 100.); // np.linspace(0, 1, 101)[1:]
                for (int i = 0; i < 100; ++i) fs[i] = mod.F(mod.ElasticUpdate(X, bs[i], de));
                int k = -1;
                for (int i = 0; i < 100; ++i)
                    if (fs[i] > 0.) {
                        k = i;
                        break;
                    }
                if (k < 0) throw TDriverError("Runge-Kutta: no intersection with the yield surface");
                lo = k > 0 ? bs[k - 1] : 0.;
                hi = bs[k];
            }
            beta = fn >= -tolf ? lo + (hi - lo) * 0.5 : fn / (fn - mod.F(Xtr));
            for (int it = 0; it < 100; ++it) {
                const TState Xb = mod.ElasticUpdate(X, beta, de);
                const REAL fb = mod.F(Xb);
                if (std::fabs(fb) <= tolf) break;
                if (fb > 0.) hi = beta;
                else lo = beta;
                const REAL nb = beta - fb / mod.DfDbeta(Xb, de);
                beta = (lo < nb && nb < hi) ? nb : 0.5 * (lo + hi);
            }
            X = mod.ElasticUpdate(X, beta, de);
        }
        TRKModel::TVec6 rem;
        for (int k = 0; k < 6; ++k) rem[k] = (1. - beta) * de[k];
        REAL T = 0., dT = 1.;
        bool failed = false;
        int natt = 0, nrej = 0;
        while (T < 1. - 1e-14) {
            if (natt > maxatt) throw TDriverError("Runge-Kutta: largest number of attempts");
            TRKModel::TVec6 d;
            for (int k = 0; k < 6; ++k) d[k] = dT * rem[k];
            TState Xs, Xh;
            mod.Stages(X, d, tab, Xs, Xh);
            natt++;
            REAL n1 = 0., n2 = 0.;
            for (int l = 0; l < 8; ++l) {
                n1 += (Xs[l] - Xh[l]) * (Xs[l] - Xh[l]);
                n2 += Xs[l] * Xs[l];
            }
            const REAL r2 = std::sqrt(n1) / std::sqrt(n2);
            if (r2 > tol) {
                const REAL q = std::max(0.1, 0.9 * std::pow(tol / r2, 1. / tab.fOrder));
                dT = q * dT;
                nrej++;
                failed = true;
                continue;
            }
            X = Xs;
            T += dT;
            if (control == EXie) {
                dT = 1. - T; // eq. (36): the rest of the increment
            } else {
                REAL q = r2 > 0. ? std::min(1.1, 0.9 * std::pow(tol / r2, 1. / tab.fOrder)) : 1.1;
                if (failed) q = std::min(q, 1.0);
                failed = false;
                dT = std::min(q * dT, 1. - T);
            }
        }
        X = mod.Drift(X);
        external(X, res);
        res.fWork = mod.NEval();
        res.fAttempts = natt;
        res.fRejections = nrej;
    }
    return res;
}

// ============================================================================================ 3. drivers

inline IntegrationSchemes::TPath IntegrationSchemes::UndrainedPath(const TUpdateFunction &update, REAL p0, REAL pc0,
                                                                   REAL eamax, int n) {
    TPath path;
    TVoigt sig{};
    sig[iXX] = sig[iYY] = sig[iZZ] = -p0;
    REAL pc = pc0;
    TVoigt eps{}, d{}, inc;
    d[iXX] = d[iYY] = 0.5;
    d[iZZ] = -1.0;
    for (int i = 0; i < 6; ++i) inc[i] = eamax / n * d[i];
    path.fRows.push_back({0., p0, 0., 0.});
    for (int k = 0; k < n; ++k) {
        const TUpdate u = update(sig, inc, pc);
        sig = u.fSigma;
        pc = u.fPc;
        for (int i = 0; i < 6; ++i) eps[i] = eps[i] + inc[i];
        path.fWork += u.fWork;
        path.fSubsteps += u.fSubsteps;
        if (u.fSubsteps > 1) path.fSubsteppedIncrements++;
        path.fAttempts += u.fAttempts;
        path.fRejections += u.fRejections;
        REAL p, q;
        PQ(sig, p, q);
        path.fRows.push_back({-eps[iZZ], p, q, -(eps[iXX] + eps[iYY] + eps[iZZ])});
    }
    return path;
}

inline IntegrationSchemes::TPath IntegrationSchemes::DrainedPath(const TUpdateFunction &update, REAL p0, REAL pc0,
                                                                 REAL eamax, int n, REAL tol) {
    const REAL nan = std::numeric_limits<REAL>::quiet_NaN();
    TPath path;
    TVoigt sig{};
    sig[iXX] = sig[iYY] = sig[iZZ] = -p0;
    REAL pc = pc0;
    path.fRows.push_back({0., p0, 0., 0.});
    const REAL dea = eamax / n;
    REAL er = 0.3 * dea;
    TVoigt eps{};
    struct TData {
        bool fOk = false;
        TUpdate fU;
        TVoigt fD{};
    };
    for (int k = 0; k < n; ++k) {
        // sigma_r + p0 for the lateral strain x (throws TProjectionError)
        auto sr = [&](REAL x, TData &data) {
            TVoigt d{};
            d[iZZ] = -dea;
            d[iXX] = d[iYY] = x;
            data.fU = update(sig, d, pc);
            data.fD = d;
            data.fOk = true;
            return data.fU.fSigma[iXX] + p0;
        };
        TData data;
        try { // Newton with a finite-difference derivative
            REAL r0 = sr(er, data);
            bool conv = false;
            for (int it = 0; it < 40; ++it) {
                if (std::fabs(r0) < tol * p0) {
                    conv = true;
                    break;
                }
                const REAL h = 1e-7 * std::max(1e-3, std::fabs(er));
                TData tmp;
                const REAL r1 = sr(er + h, tmp);
                er = er - r0 * h / (r1 - r0);
                r0 = sr(er, data);
            }
            if (!conv) data.fOk = false;
        } catch (TProjectionError &) {
            data.fOk = false;
        }
        if (!data.fOk) { // safeguard: bisection in an interval where sigma_r + p0 changes sign
            auto g = [&](REAL x) {
                TData tmp;
                try {
                    return sr(x, tmp);
                } catch (TProjectionError &) {
                    return nan;
                }
            };
            const int nx = 81;
            std::vector<REAL> xs(nx), gs(nx);
            const REAL start = -2. * dea, stop = 2. * dea, step = (stop - start) / (nx - 1);
            for (int i = 0; i < nx; ++i) xs[i] = i == nx - 1 ? stop : i * step + start; // np.linspace
            for (int i = 0; i < nx; ++i) gs[i] = g(xs[i]);
            auto sign = [](REAL v) { return v > 0. ? 1. : (v < 0. ? -1. : 0.); };
            int ok = -1;
            for (int i = 0; i + 1 < nx; ++i)
                if (std::isfinite(gs[i]) && std::isfinite(gs[i + 1]) && sign(gs[i]) != sign(gs[i + 1])) {
                    ok = i;
                    break;
                }
            if (ok < 0) throw TDriverError("drained driver did not converge");
            REAL a = xs[ok], b = xs[ok + 1], c = a;
            bool balanced = false;
            for (int it = 0; it < 200; ++it) {
                c = 0.5 * (a + b);
                const REAL gc = g(c);
                if (!std::isfinite(gc)) throw TDriverError("drained driver: projection failed in the interval");
                if (std::fabs(gc) < tol * p0) {
                    balanced = true;
                    break;
                }
                if (sign(gc) == sign(g(a))) a = c;
                else b = c;
            }
            // as in rivais.py, the last midpoint is accepted even if the bisection stopped without equilibrium (it
            // closes in on a jump of sigma_r + p0 when the update changes branch); such increments are counted
            path.fBisectionIncrements++;
            if (!balanced) path.fUnbalancedIncrements++;
            er = c;
            sr(er, data); // may throw TProjectionError (no solution)
        }
        path.fMaxResidual = std::max(path.fMaxResidual, std::fabs(data.fU.fSigma[iXX] + p0));
        sig = data.fU.fSigma;
        pc = data.fU.fPc;
        for (int i = 0; i < 6; ++i) eps[i] = eps[i] + data.fD[i];
        path.fWork += data.fU.fWork;
        path.fSubsteps += data.fU.fSubsteps;
        if (data.fU.fSubsteps > 1) path.fSubsteppedIncrements++;
        path.fAttempts += data.fU.fAttempts;
        path.fRejections += data.fU.fRejections;
        REAL p, q;
        PQ(sig, p, q);
        path.fRows.push_back({-eps[iZZ], p, q, -(eps[iXX] + eps[iYY] + eps[iZZ])});
    }
    return path;
}

inline void IntegrationSchemes::UndrainedClosed(REAL p0, REAL pc0, const TMaterial &P, REAL eta, REAL &p, REAL &q,
                                                REAL &epsq) {
    const REAL M = P.fM, lam = P.fLambda, kap = P.fKappa, v0 = P.fV0;
    const REAL Lam = (lam - kap) / lam;
    const REAL R = pc0 / p0;
    const REAL etay = M * std::sqrt(std::max(R - 1., 0.));
    p = p0 * std::pow(((M * M) + (eta * eta)) / ((M * M) * R), -Lam);
    q = eta * p;
    auto g = [M](REAL x) { return 0.5 * std::log(std::fabs((M + x) / (M - x))) - std::atan(x / M); };
    const REAL eqp = 2. * Lam * kap / (M * v0) * (g(eta) - g(etay));
    REAL eqe;
    if (P.ConstantG()) {
        eqe = q / (3. * P.fG);
    } else { // d eps_q^e = kappa dq / (3 r v0 p'), with q = eta p' and ln p'(eta)
        const REAL r = P.GFactor();
        const REAL el = etay * kap / (3. * r * v0); // elastic part, p' = p0: q_y / (3 G0)
        const REAL dat = std::atan(eta / M) - std::atan(etay / M);
        eqe = el + kap / (3. * r * v0) * ((eta - etay) - 2. * Lam * ((eta - etay) - M * dat));
    }
    epsq = eqp + eqe;
}

inline REAL IntegrationSchemes::Brent(const std::function<REAL(REAL)> &f, REAL xa, REAL xb, REAL xtol, REAL rtol,
                                      int maxiter) {
    // transcription of scipy/optimize/Zeros/brentq.c
    REAL xpre = xa, xcur = xb, xblk = 0., fblk = 0., spre = 0., scur = 0.;
    REAL fpre = f(xpre), fcur = f(xcur);
    if (fpre == 0.) return xpre;
    if (fcur == 0.) return xcur;
    if (std::signbit(fpre) == std::signbit(fcur)) throw TDriverError("Brent: f(a) and f(b) must have different signs");
    for (int i = 0; i < maxiter; ++i) {
        if (fpre != 0. && fcur != 0. && std::signbit(fpre) != std::signbit(fcur)) {
            xblk = xpre;
            fblk = fpre;
            spre = scur = xcur - xpre;
        }
        if (std::fabs(fblk) < std::fabs(fcur)) {
            xpre = xcur;
            xcur = xblk;
            xblk = xpre;
            fpre = fcur;
            fcur = fblk;
            fblk = fpre;
        }
        const REAL delta = (xtol + rtol * std::fabs(xcur)) / 2.;
        const REAL sbis = (xblk - xcur) / 2.;
        if (fcur == 0. || std::fabs(sbis) < delta) return xcur;
        if (std::fabs(spre) > delta && std::fabs(fcur) < std::fabs(fpre)) {
            REAL stry;
            if (xpre == xblk) { // interpolate
                stry = -fcur * (xcur - xpre) / (fcur - fpre);
            } else { // extrapolate
                const REAL dpre = (fpre - fcur) / (xpre - xcur);
                const REAL dblk = (fblk - fcur) / (xblk - xcur);
                stry = -fcur * (fblk * dblk - fpre * dpre) / (dblk * dpre * (fblk - fpre));
            }
            if (2. * std::fabs(stry) < std::min(std::fabs(spre), 3. * std::fabs(sbis) - delta)) { // good short step
                spre = scur;
                scur = stry;
            } else { // bisect
                spre = sbis;
                scur = sbis;
            }
        } else { // bisect
            spre = sbis;
            scur = sbis;
        }
        xpre = xcur;
        fpre = fcur;
        if (std::fabs(scur) > delta) xcur += scur;
        else xcur += (sbis > 0. ? delta : -delta);
        fcur = f(xcur);
    }
    throw TDriverError("Brent: no convergence");
}

// ============================================================================================ update functions

inline IntegrationSchemes::TUpdateFunction IntegrationSchemes::BEUpdate(const TMaterial &P, EVolumetric volumetric,
                                                                        EShearModulus shear, bool substep) {
    return [P, volumetric, shear, substep](const TVoigt &s, const TVoigt &d, REAL pc) {
        return substep ? BackwardEulerSubstep(P, s, d, pc, volumetric, shear)
                       : BackwardEuler(P, s, d, pc, volumetric, shear);
    };
}

inline IntegrationSchemes::TUpdateFunction IntegrationSchemes::RKUpdate(const TMaterial &P, REAL e0,
                                                                        const TTableau &tab, REAL tol,
                                                                        EStepControl control, ERKModel model) {
    return [P, e0, tab, tol, control, model](const TVoigt &s, const TVoigt &d, REAL pc) {
        return RungeKutta(P, s, d, pc, e0, tab, tol, control, model);
    };
}

inline IntegrationSchemes::TUpdateFunction
IntegrationSchemes::LibraryUpdate(const TMaterial &P, TPZYCModifiedCamClayRHW::EPorousIntegration integration) {
    mcc::TPlastic model;
    model.SetModifiedCamClay(P.fM, P.fLambda, P.fKappa, P.fPt, P.fOmega);
    model.SetPorousElasticity();
    if (P.ConstantG()) model.SetConstantShearModulus(P.fG);
    else model.SetPoissonRatio(P.fNu);
    model.SetDefaultSpecificVolume(P.fV0);
    model.SetPorousIntegration(integration);
    const REAL v0 = P.fV0;
    return [model, v0](const TVoigt &s, const TVoigt &d, REAL pc) {
        TPZTensor<REAL> sign, eps, epsn, sigma;
        for (int i = 0; i < 6; ++i) {
            sign[i] = s[i];
            eps[i] = d[i]; // total strain = increment, from eps_n = 0: deps = eps - eps_n is exact
        }
        const mcc::TPointState st(sign, pc, v0);
        TPZFNMatrix<36, REAL> D(6, 6, 0.);
        REAL pcnew = pc;
        int type = 0;
        mcc::TLocalStats stats;
        if (!mcc::ApplyStrain(model, epsn, st, eps, sigma, D, pcnew, type, &stats))
            throw TProjectionError("TPZPlasticStepModifiedCamClay: projection failed");
        TUpdate u;
        for (int i = 0; i < 6; ++i) u.fSigma[i] = sigma[i];
        u.fPc = pcnew;
        u.fWork = stats.fIts;
        return u;
    };
}

// ============================================================================================ tests

inline const IntegrationSchemes::TReference &IntegrationSchemes::Reference(const std::string &test) const {
    for (auto &r : fReferences)
        if (r.fTest == test) return r;
    throw std::logic_error("no reference for test " + test);
}

inline IntegrationSchemes::TRecord &IntegrationSchemes::Store(TRecord rec) {
    fRecords.push_back(std::move(rec));
    return fRecords.back();
}

inline void IntegrationSchemes::RunPath(TRecord &rec, const std::function<TPath()> &path, const TReference &ref) {
    const auto t0 = std::chrono::steady_clock::now();
    try {
        rec.fPath = path();
        rec.fConverged = true;
    } catch (TProjectionError &e) {
        rec.fConverged = false;
        rec.fMessage = e.what();
    } catch (TDriverError &e) {
        rec.fConverged = false;
        rec.fMessage = e.what();
    }
    rec.fTime = std::chrono::duration<REAL>(std::chrono::steady_clock::now() - t0).count();
    if (!rec.fConverged) return;
    const TPath &p = rec.fPath;
    const auto &last = p.fRows.back();
    rec.fWork = p.fWork;
    rec.fSubsteps = p.fSubsteps;
    rec.fSubsteppedIncrements = p.fSubsteppedIncrements;
    rec.fAttempts = p.fAttempts;
    rec.fRejections = p.fRejections;
    rec.fBisectionIncrements = p.fBisectionIncrements;
    rec.fUnbalancedIncrements = p.fUnbalancedIncrements;
    if (rec.fTest.rfind("kl_drained", 0) == 0) rec.fMaxResidual = p.fMaxResidual;
    rec.fPEnd = last[1];
    rec.fQEnd = last[2];
    rec.fEpsVEnd = last[3];
    if (rec.fTest.rfind("kl_undrained", 0) == 0) {
        rec.fV[0] = last[1] - ref.fP;
        rec.fV[1] = last[2] - ref.fQ;
    } else if (rec.fTest.rfind("kl_drained", 0) == 0) {
        REAL qmax = -std::numeric_limits<REAL>::infinity();
        for (auto &r : p.fRows) qmax = std::max(qmax, r[2]);
        rec.fV[0] = last[2] - ref.fQ;
        rec.fV[1] = last[3] - ref.fEpsV;
        rec.fV[2] = qmax;
    } else { // test B of Xie et al.: relative error of the stress (axisymmetric state with the final p', q)
        const TVoigt s = Axisymmetric(last[1], last[2]), sref = Axisymmetric(ref.fP, ref.fQ);
        REAL n1 = 0., n2 = 0.;
        for (int i = 0; i < 6; ++i) {
            n1 += (s[i] - sref[i]) * (s[i] - sref[i]);
            n2 += sref[i] * sref[i];
        }
        rec.fV[0] = std::sqrt(n1) / std::sqrt(n2);
    }
}

inline void IntegrationSchemes::RunXieTestB() {
    const TMaterial P = XieMaterial();
    const REAL p0 = fXieP0, e0 = fXieE0;
    // reference: closed form at eps_q = eps_a = 10% (Brent on eta, as brentq of gen_data.py)
    TReference ref;
    ref.fTest = "xieB";
    ref.fEpsA = 0.1;
    ref.fEta = Brent([&](REAL x) {
        REAL p, q, eq;
        UndrainedClosed(p0, p0, P, x, p, q, eq);
        return eq - 0.1;
    }, 1e-12, P.fM - 1e-13, 1e-15);
    REAL eq;
    UndrainedClosed(p0, p0, P, ref.fEta, ref.fP, ref.fQ, eq);
    fReferences.push_back(ref);
    // implicit schemes, 1 to 1024 increments
    struct TBE {
        const char *fName;
        std::function<TUpdateFunction()> fUpdate;
    };
    const std::vector<TBE> bes = {
        {"exact_n", [&]() { return BEUpdate(P, EExactK, EShearAtPn, false); }},
        {"exact_secant", [&]() { return BEUpdate(P, EExactK, ESecantShear, true); }},
        {"frozen_n", [&]() { return BEUpdate(P, EFrozenK, EShearAtPn, false); }},
        {"library_exact", [&]() { return LibraryUpdate(P, TPZYCModifiedCamClayRHW::EExact); }},
        {"library_frozen", [&]() { return LibraryUpdate(P, TPZYCModifiedCamClayRHW::EFrozen); }}};
    for (auto &be : bes)
        for (int n = 1; n <= 1024; n *= 2) {
            TRecord rec;
            rec.fTest = ref.fTest;
            rec.fScheme = be.fName;
            rec.fControl = rec.fModel = "-";
            rec.fN = n;
            const TUpdateFunction up = be.fUpdate();
            RunPath(rec, [&]() { return UndrainedPath(up, p0, p0, 0.1, n); }, ref);
            Store(rec);
        }
    // explicit schemes: one increment of 10%, Sloan control (tolerances 1e-1 to 1e-8) and the control and model of
    // Xie et al. (validation, tolerances 1e-1 to 1e-5)
    TVoigt d1{};
    d1[iZZ] = -0.1;
    d1[iXX] = d1[iYY] = 0.05;
    TVoigt s0{};
    s0[iXX] = s0[iYY] = s0[iZZ] = -p0;
    const TVoigt sref = Axisymmetric(ref.fP, ref.fQ);
    for (const TTableau &tab : {ME21(), RKDP54()})
        for (REAL tol : {1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8})
            for (EStepControl control : {ESloan, EXie}) {
                if (control == EXie && tol < 1e-5) continue;
                TRecord rec;
                rec.fTest = ref.fTest;
                rec.fScheme = tab.fName;
                rec.fControl = control == ESloan ? "sloan" : "xie";
                rec.fModel = control == ESloan ? "v0" : "xie";
                rec.fN = 1;
                rec.fTol = tol;
                const auto t0 = std::chrono::steady_clock::now();
                try {
                    const TUpdate u = RungeKutta(P, s0, d1, p0, e0, tab, tol, control,
                                                 control == ESloan ? EConstantV : EVoidRatio);
                    REAL n1 = 0., n2 = 0.;
                    for (int i = 0; i < 6; ++i) {
                        n1 += (u.fSigma[i] - sref[i]) * (u.fSigma[i] - sref[i]);
                        n2 += sref[i] * sref[i];
                    }
                    rec.fConverged = true;
                    rec.fV[0] = std::sqrt(n1) / std::sqrt(n2);
                    rec.fWork = u.fWork;
                    rec.fAttempts = u.fAttempts;
                    rec.fRejections = u.fRejections;
                    PQ(u.fSigma, rec.fPEnd, rec.fQEnd);
                } catch (std::runtime_error &e) {
                    rec.fConverged = false;
                    rec.fMessage = e.what();
                }
                rec.fTime = std::chrono::duration<REAL>(std::chrono::steady_clock::now() - t0).count();
                Store(rec);
            }
}

inline void IntegrationSchemes::RunKLUndrained(int ocr, REAL eamax) {
    const TMaterial P = KLMaterial();
    const REAL p0 = fKLP0, e0 = fKLE0, pc0 = ocr * p0;
    TReference ref;
    ref.fTest = "kl_undrained_ocr" + std::to_string(ocr);
    ref.fEpsA = eamax;
    const REAL etay = P.fM * std::sqrt(std::max(REAL(ocr - 1), 0.));
    const REAL lo = ocr < 2 ? 1e-12 : P.fM + 1e-13, hi = ocr < 2 ? P.fM - 1e-13 : etay - 1e-13;
    ref.fEta = Brent([&](REAL x) {
        REAL p, q, eq;
        UndrainedClosed(p0, pc0, P, x, p, q, eq);
        return eq - eamax;
    }, lo, hi, 1e-15);
    REAL eq;
    UndrainedClosed(p0, pc0, P, ref.fEta, ref.fP, ref.fQ, eq);
    fReferences.push_back(ref);
    const std::vector<std::pair<const char *, TUpdateFunction>> bes = {
        {"exact_n", BEUpdate(P, EExactK, EShearAtPn, false)},
        {"exact_secant", BEUpdate(P, EExactK, ESecantShear, true)},
        {"frozen_n", BEUpdate(P, EFrozenK, EShearAtPn, false)},
        {"library_exact", LibraryUpdate(P, TPZYCModifiedCamClayRHW::EExact)},
        {"library_frozen", LibraryUpdate(P, TPZYCModifiedCamClayRHW::EFrozen)}};
    for (auto &be : bes)
        for (int n : {5, 10, 20, 50, 100, 200, 500, 1000}) {
            TRecord rec;
            rec.fTest = ref.fTest;
            rec.fScheme = be.first;
            rec.fControl = rec.fModel = "-";
            rec.fN = n;
            RunPath(rec, [&]() { return UndrainedPath(be.second, p0, pc0, eamax, n); }, ref);
            Store(rec);
        }
    for (const TTableau &tab : {ME21(), RKDP54()})
        for (REAL tol : {1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7}) {
            TRecord rec;
            rec.fTest = ref.fTest;
            rec.fScheme = tab.fName;
            rec.fControl = "sloan";
            rec.fModel = "v0";
            rec.fN = 10;
            rec.fTol = tol;
            const TUpdateFunction up = RKUpdate(P, e0, tab, tol);
            RunPath(rec, [&]() { return UndrainedPath(up, p0, pc0, eamax, 10); }, ref);
            Store(rec);
        }
}

inline void IntegrationSchemes::RunKLDrained(int ocr) {
    const TMaterial P = KLMaterial();
    const REAL p0 = fKLP0, e0 = fKLE0, pc0 = ocr * p0, eamax = 0.25;
    TReference ref;
    ref.fTest = "kl_drained_ocr" + std::to_string(ocr);
    ref.fEpsA = eamax;
    // closed form of Appendix B.1 with 20000 points (triaxial_closed), interpolated at eps_a = 25% (numpy.interp)
    const auto ana = mcc::TriaxialDrainedClosed(p0, pc0, P.fV0, P.fM, P.fLambda, P.fKappa, 0., P.fNu, 20000);
    ref.fQ = mcc::Interpolate(ana, eamax, 2);
    ref.fEpsV = mcc::Interpolate(ana, eamax, 3);
    ref.fQPeak = -std::numeric_limits<REAL>::infinity();
    for (auto &r : ana) ref.fQPeak = std::max(ref.fQPeak, r[2]);
    fReferences.push_back(ref);
    const std::vector<std::pair<const char *, TUpdateFunction>> bes = {
        {"exact_n", BEUpdate(P, EExactK, EShearAtPn, false)},
        {"exact_secant", BEUpdate(P, EExactK, ESecantShear, true)},
        {"frozen_n", BEUpdate(P, EFrozenK, EShearAtPn, false)},
        {"library_exact", LibraryUpdate(P, TPZYCModifiedCamClayRHW::EExact)},
        {"library_frozen", LibraryUpdate(P, TPZYCModifiedCamClayRHW::EFrozen)}};
    for (auto &be : bes)
        for (int n : {10, 20, 50, 100, 200, 500}) {
            TRecord rec;
            rec.fTest = ref.fTest;
            rec.fScheme = be.first;
            rec.fControl = rec.fModel = "-";
            rec.fN = n;
            RunPath(rec, [&]() { return DrainedPath(be.second, p0, pc0, eamax, n); }, ref);
            Store(rec);
        }
    // RK with tolerance 1e-4 (data of the article) and 1e-5 (statement of Sect. 6.2: the same errors)
    for (REAL tol : {1e-4, 1e-5})
        for (const TTableau &tab : {ME21(), RKDP54()})
            for (int n : {10, 20, 50, 100}) {
                TRecord rec;
                rec.fTest = ref.fTest;
                rec.fScheme = tab.fName;
                rec.fControl = "sloan";
                rec.fModel = "v0";
                rec.fN = n;
                rec.fTol = tol;
                const TUpdateFunction up = RKUpdate(P, e0, tab, tol);
                RunPath(rec, [&]() { return DrainedPath(up, p0, pc0, eamax, n); }, ref);
                Store(rec);
            }
}

inline void IntegrationSchemes::RunSecantSweep() {
    const TMaterial P = KLMaterial();
    const REAL p0 = fKLP0, pc0 = 10. * p0, dea = 0.25 / 10.;
    TVoigt s0{};
    s0[iXX] = s0[iYY] = s0[iZZ] = -p0;
    fSweep.clear();
    for (int i = 0; i <= 2000; ++i) {
        const REAL x = 0.010 + i * 1e-5;
        TVoigt d{};
        d[iZZ] = -dea;
        d[iXX] = d[iYY] = x;
        for (int scheme = 0; scheme < 2; ++scheme) {
            TSweepPoint pt;
            pt.fEpsR = x;
            pt.fScheme = scheme;
            try {
                const TUpdate u = BackwardEuler(P, s0, d, pc0, EExactK, scheme ? ESecantShear : EShearAtPn);
                pt.fConverged = true;
                pt.fIterations = u.fWork;
                PQ(u.fSigma, pt.fP, pt.fQ);
                pt.fSigmaR = u.fSigma[iXX] + p0;
            } catch (TProjectionError &) {
                pt.fConverged = false;
            }
            fSweep.push_back(pt);
        }
    }
}

inline void IntegrationSchemes::RunToleranceStudy() {
    const TMaterial P = XieMaterial();
    const TReference &ref = Reference("xieB");
    const TVoigt sref = Axisymmetric(ref.fP, ref.fQ);
    fTolerance.clear();
    for (int n : {1, 4, 16, 64})
        for (REAL tol : {1e-2, 1e-4, 1e-6, 1e-8, 1e-10, 1e-12, 1e-14}) {
            TTolerancePoint pt;
            pt.fN = n;
            pt.fTol = tol;
            const TUpdateFunction up = [P, tol](const TVoigt &s, const TVoigt &d, REAL pc) {
                return BackwardEuler(P, s, d, pc, EExactK, EShearAtPn, tol);
            };
            try {
                const TPath path = UndrainedPath(up, fXieP0, fXieP0, 0.1, n);
                const TVoigt sv = Axisymmetric(path.fRows.back()[1], path.fRows.back()[2]);
                REAL n1 = 0., n2 = 0.;
                for (int i = 0; i < 6; ++i) {
                    n1 += (sv[i] - sref[i]) * (sv[i] - sref[i]);
                    n2 += sref[i] * sref[i];
                }
                pt.fConverged = true;
                pt.fError = std::sqrt(n1) / std::sqrt(n2);
                pt.fWork = path.fWork;
            } catch (TProjectionError &) {
                pt.fConverged = false;
            }
            fTolerance.push_back(pt);
        }
}

inline const IntegrationSchemes::TRecord *IntegrationSchemes::Find(const std::string &test, const std::string &scheme,
                                                                   int n, REAL tol,
                                                                   const std::string &control) const {
    for (auto &r : fRecords) {
        if (r.fTest != test || r.fScheme != scheme || r.fN != n) continue;
        if (!control.empty() && r.fControl != control) continue;
        if (tol > 0. ? std::fabs(r.fTol - tol) > 1e-6 * tol : r.fTol != 0.) continue;
        return &r;
    }
    return nullptr;
}

// ============================================================================================ post-processing

inline std::string IntegrationSchemes::Num(REAL v) {
    if (std::isnan(v)) return "nan";
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%.17g", v);
    return buf;
}

inline void IntegrationSchemes::WriteTable(const std::string &file, const std::vector<std::string> &header,
                                           const std::vector<std::vector<std::string>> &rows) {
    std::ofstream out(file);
    for (size_t i = 0; i < header.size(); ++i) out << header[i] << (i + 1 < header.size() ? "," : "\n");
    for (auto &r : rows)
        for (size_t i = 0; i < r.size(); ++i) out << r[i] << (i + 1 < r.size() ? "," : "\n");
}

inline REAL IntegrationSchemes::FigureError(const TRecord &rec) {
    if (!rec.fConverged) return std::numeric_limits<REAL>::quiet_NaN();
    return rec.fTest == "xieB" ? rec.fV[0] : std::fabs(rec.fV[0]);
}

inline void IntegrationSchemes::PostProcess() const {
    // 1. all the runs
    std::vector<std::vector<std::string>> rows;
    for (auto &r : fRecords) {
        const bool rk = r.fControl != "-";
        rows.push_back({r.fTest, r.fScheme, rk ? "RK" : (r.fScheme.rfind("library", 0) == 0 ? "library" : "BE"),
                        r.fControl, r.fModel, std::to_string(r.fN), Num(r.fTol), std::to_string(int(r.fConverged)),
                        Num(r.fV[0]), Num(r.fV[1]), Num(r.fV[2]), std::to_string(r.fWork),
                        std::to_string(r.fSubsteps), std::to_string(r.fSubsteppedIncrements),
                        std::to_string(r.fAttempts), std::to_string(r.fRejections), Num(r.fPEnd), Num(r.fQEnd),
                        Num(r.fEpsVEnd), Num(r.fTime), std::to_string(r.fBisectionIncrements),
                        std::to_string(r.fUnbalancedIncrements), r.fConverged ? Num(r.fMaxResidual) : "nan"});
    }
    WriteTable("schemes_results.csv",
               {"test", "scheme", "family", "control", "model", "n", "tol", "converged", "v0", "v1", "v2", "work",
                "substeps", "substepped_increments", "rk_attempts", "rk_rejections", "p_end", "q_end", "eps_v_end",
                "time_s", "bisection_increments", "unbalanced_increments", "max_abs_sigma_r_plus_p0"},
               rows);
    // 2. Fig. 6: (panel, scheme, parameter, work, error)
    rows.clear();
    const std::vector<std::pair<std::string, std::string>> panels = {
        {"a", "xieB"}, {"b", "kl_undrained_ocr10"}, {"c", "kl_drained_ocr1"}};
    for (auto &pan : panels)
        for (auto &r : fRecords) {
            if (r.fTest != pan.second || r.fControl == "xie") continue;
            if (r.fScheme.rfind("library", 0) == 0) continue;
            if (pan.first == "c" && r.fControl == "sloan" && std::fabs(r.fTol - 1e-4) > 1e-12) continue;
            rows.push_back({pan.first, r.fTest, r.fScheme, r.fControl == "-" ? "BE" : "RK", std::to_string(r.fN),
                            Num(r.fTol), std::to_string(int(r.fConverged)), std::to_string(r.fWork),
                            Num(FigureError(r))});
        }
    WriteTable("schemes_fig05.csv", {"panel", "test", "scheme", "family", "n", "tol", "converged", "work", "error"},
               rows);
    // 3. Table 4
    rows.clear();
    const std::vector<std::pair<std::string, std::string>> t4 = {{"exact_n", "This work: BE, exact K, G(p_n)"},
                                                                 {"exact_secant", "BE, exact K, secant G"},
                                                                 {"frozen_n", "BE, frozen K [26]"},
                                                                 {"ME2(1)", "RK ME2(1) [14]"},
                                                                 {"RKDP5(4)", "RK RKDP5(4) [16]"},
                                                                 {"library_exact", "TPZPlasticStepModifiedCamClay, exact K"},
                                                                 {"library_frozen", "TPZPlasticStepModifiedCamClay, frozen K"}};
    for (auto &s : t4) {
        const bool rk = s.first[0] == 'M' || s.first[0] == 'R';
        const TRecord *a = rk ? Find("xieB", s.first, 1, 1e-4, "sloan") : Find("xieB", s.first, 1);
        const TRecord *b = rk ? Find("kl_undrained_ocr10", s.first, 10, 1e-4) : Find("kl_undrained_ocr10", s.first, 10);
        const TRecord *c = rk ? Find("kl_drained_ocr1", s.first, 10, 1e-4) : Find("kl_drained_ocr1", s.first, 10);
        auto v = [](const TRecord *r, int i, REAL scale = 1.) {
            return (r && r->fConverged) ? Num(scale * r->fV[i]) : std::string("nan");
        };
        auto w = [](const TRecord *r) { return (r && r->fConverged) ? std::to_string(r->fWork) : std::string("nan"); };
        rows.push_back({s.first, "\"" + s.second + "\"", v(a, 0), w(a), v(b, 0), v(b, 1), w(b), v(c, 0),
                        v(c, 1, 100.), w(c), std::to_string(int(c && c->fConverged))});
    }
    WriteTable("schemes_table4.csv",
               {"scheme", "label", "a_xieB_error", "a_work", "b_undrained_ocr10_dp", "b_undrained_ocr10_dq", "b_work",
                "c_drained_nc_dq", "c_drained_nc_deps_v_percent", "c_work", "c_converged"},
               rows);
    // 4. references
    rows.clear();
    for (auto &r : fReferences)
        rows.push_back({r.fTest, Num(r.fEpsA), Num(r.fEta), Num(r.fP), Num(r.fQ), Num(r.fEpsV), Num(r.fQPeak)});
    WriteTable("schemes_references.csv", {"test", "eps_a", "eta", "p_eff", "q", "eps_v", "q_peak"}, rows);
    // 5. library against be_tensor (exact and frozen)
    rows.clear();
    for (auto &r : fRecords) {
        const bool lib = r.fScheme.rfind("library", 0) == 0;
        if (!lib) continue;
        const TRecord *t = Find(r.fTest, r.fScheme == "library_exact" ? "exact_n" : "frozen_n", r.fN);
        if (!t) continue;
        REAL dmax = std::numeric_limits<REAL>::quiet_NaN();
        if (r.fConverged && t->fConverged && r.fPath.fRows.size() == t->fPath.fRows.size()) {
            dmax = 0.;
            for (size_t i = 0; i < r.fPath.fRows.size(); ++i)
                for (int c = 1; c < 3; ++c)
                    dmax = std::max(dmax, std::fabs(r.fPath.fRows[i][c] - t->fPath.fRows[i][c]));
        }
        rows.push_back({r.fTest, r.fScheme, t->fScheme, std::to_string(r.fN), std::to_string(int(r.fConverged)),
                        std::to_string(int(t->fConverged)), Num(r.fPEnd), Num(t->fPEnd), Num(r.fQEnd), Num(t->fQEnd),
                        Num(dmax), std::to_string(r.fWork), std::to_string(t->fWork)});
    }
    WriteTable("schemes_library_check.csv",
               {"test", "library", "be_tensor", "n", "converged_library", "converged_be_tensor", "p_end_library",
                "p_end_be_tensor", "q_end_library", "q_end_be_tensor", "max_abs_diff_p_q_along_path", "work_library",
                "work_be_tensor"},
               rows);
    // 6. comparison with data_rivais.pkl
    rows.clear();
    for (int i = 0; i < gIntegrationSchemesPythonSize; ++i) {
        const TIntegrationSchemesPythonRow &py = gIntegrationSchemesPython[i];
        std::vector<std::string> row = {py.fTest, py.fScheme, py.fControl, std::to_string(py.fN), Num(py.fTol)};
        REAL cpp[3] = {std::numeric_limits<REAL>::quiet_NaN(), std::numeric_limits<REAL>::quiet_NaN(),
                       std::numeric_limits<REAL>::quiet_NaN()};
        long work = -1;
        int conv = 0;
        if (std::string(py.fScheme) == "reference") {
            const TReference &ref = Reference(py.fTest);
            const bool drained = std::string(py.fTest).rfind("kl_drained", 0) == 0;
            cpp[0] = drained ? ref.fQ : ref.fP;
            cpp[1] = drained ? ref.fEpsV : ref.fQ;
            cpp[2] = drained ? ref.fQPeak : std::numeric_limits<REAL>::quiet_NaN();
            work = 0;
            conv = 1;
        } else if (const TRecord *r = Find(py.fTest, py.fScheme, py.fN, py.fTol,
                                           std::string(py.fControl) == "-" ? "" : py.fControl)) {
            conv = r->fConverged;
            work = r->fWork;
            for (int k = 0; k < 3; ++k) cpp[k] = r->fV[k];
        } else {
            conv = -1;
        }
        row.push_back(std::to_string(py.fConverged));
        row.push_back(std::to_string(conv));
        for (int k = 0; k < 3; ++k) {
            row.push_back(Num(py.fV[k]));
            row.push_back(Num(cpp[k]));
            row.push_back(Num(cpp[k] - py.fV[k]));
        }
        row.push_back(std::to_string(py.fWork));
        row.push_back(std::to_string(work));
        rows.push_back(row);
    }
    WriteTable("schemes_python_comparison.csv",
               {"test", "scheme", "control", "n", "tol", "converged_python", "converged_cpp", "v0_python", "v0_cpp",
                "v0_diff", "v1_python", "v1_cpp", "v1_diff", "v2_python", "v2_cpp", "v2_diff", "work_python",
                "work_cpp"},
               rows);
    // 7. single step in the first increment of the drained test with OCR = 10 (RunSecantSweep)
    rows.clear();
    for (auto &pt : fSweep)
        rows.push_back({Num(pt.fEpsR), pt.fScheme ? "exact_secant_single_step" : "exact_n",
                        std::to_string(int(pt.fConverged)), std::to_string(pt.fIterations),
                        std::to_string(int(pt.fConverged && pt.fP < 1e-6 * fKLP0)),
                        pt.fConverged ? Num(pt.fP) : "nan", pt.fConverged ? Num(pt.fQ) : "nan",
                        pt.fConverged ? Num(pt.fSigmaR) : "nan"});
    WriteTable("schemes_secant_first_increment_sweep.csv",
               {"eps_r", "scheme", "converged", "iterations", "spurious_p_zero", "p_eff", "q", "sigma_r_plus_p0"}, rows);
    // 8. error of this work against the tolerance of the local Newton iterations (RunToleranceStudy)
    rows.clear();
    for (auto &pt : fTolerance)
        rows.push_back({std::to_string(pt.fN), Num(pt.fTol), std::to_string(int(pt.fConverged)),
                        pt.fConverged ? Num(pt.fError) : std::string("nan"), std::to_string(pt.fWork)});
    WriteTable("schemes_newton_tolerance.csv", {"n", "newton_tol", "converged", "error", "work"}, rows);
    // 9. paths of the Krabbenhoft and Lyamin tests with 10 increments (and of test B with 1 to 1024 increments)
    for (const std::string test : {"xieB", "kl_undrained_ocr1", "kl_undrained_ocr10", "kl_drained_ocr1",
                                   "kl_drained_ocr10"}) {
        rows.clear();
        for (auto &r : fRecords) {
            if (r.fTest != test || !r.fConverged || r.fControl == "xie" || r.fScheme.rfind("library", 0) == 0)
                continue;
            for (size_t i = 0; i < r.fPath.fRows.size(); ++i) {
                const auto &w = r.fPath.fRows[i];
                rows.push_back({r.fScheme, std::to_string(r.fN), Num(r.fTol), std::to_string(i), Num(w[0]), Num(w[1]),
                                Num(w[2]), Num(w[3])});
            }
        }
        WriteTable("schemes_paths_" + test + ".csv",
                   {"scheme", "n", "tol", "increment", "eps_a", "p_eff", "q", "eps_v"}, rows);
    }
}

inline void IntegrationSchemes::PrintTable4() const {
    std::cout << "\nTable 4: accuracy and work at a material point. Values: this work [article v0.6]\n";
    std::cout << "  (a) test B of Xie et al., one increment of eps_a = 10%: relative error of the stress, work\n";
    std::cout << "  (b) K&L example 1, undrained, OCR = 10, eps_a = 8% in 10 increments: p' - p'_exact (kPa), work\n";
    std::cout << "  (c) K&L example 1, drained, NC, eps_a = 25% in 10 increments: q - q_exact (kPa), "
                 "eps_v - eps_v,exact (%), work\n";
    struct TArt {
        const char *fScheme, *fLabel, *fA, *fB, *fC;
    };
    const TArt art[] = {
        {"exact_n", "This work: BE, exact K, G(pn)", "8.04e-02 9", "-1.70 56", "-3.59 -0.09 84"},
        {"exact_secant", "BE, exact K, secant G", "8.12e-02 9", "-1.61 56", "-3.58 -0.09 85"},
        {"frozen_n", "BE, frozen K [26]", "4.90e-02 9", "-28.74 56", "no convergence"},
        {"ME2(1)", "RK ME2(1) [14]", "3.25e-05 288", "-0.044 252", "-2.01 -0.05 734"},
        {"RKDP5(4)", "RK RKDP5(4) [16]", "2.35e-07 150", "-0.010 72", "-2.01 -0.05 432"},
        {"library_exact", "TPZPlasticStepModifiedCamClay", "(= this work)", "", ""},
        {"library_frozen", "TPZPlasticStepModifiedCamClay, EFrozen", "(= frozen K)", "", ""}};
    for (auto &a : art) {
        const bool rk = a.fScheme[0] == 'M' || a.fScheme[0] == 'R';
        const TRecord *ra = rk ? Find("xieB", a.fScheme, 1, 1e-4, "sloan") : Find("xieB", a.fScheme, 1);
        const TRecord *rb = rk ? Find("kl_undrained_ocr10", a.fScheme, 10, 1e-4) : Find("kl_undrained_ocr10", a.fScheme, 10);
        const TRecord *rc = rk ? Find("kl_drained_ocr1", a.fScheme, 10, 1e-4) : Find("kl_drained_ocr1", a.fScheme, 10);
        std::ostringstream sa, sb, sc;
        sa << std::scientific << std::setprecision(2) << ra->fV[0] << " " << ra->fWork;
        sb << std::fixed << std::setprecision(rk ? 3 : 2) << rb->fV[0] << " " << rb->fWork;
        if (rc->fConverged)
            sc << std::fixed << std::setprecision(2) << rc->fV[0] << " " << 100. * rc->fV[1] << " " << rc->fWork;
        else
            sc << "no convergence";
        std::cout << "  " << std::left << std::setw(40) << a.fLabel << std::setw(30)
                  << (sa.str() + " [" + a.fA + "]") << std::setw(26) << (sb.str() + " [" + a.fB + "]")
                  << sc.str() << " [" << a.fC << "]" << std::right << "\n";
    }
}

inline void IntegrationSchemes::PrintTextNumbers() const {
    auto rec = [&](const std::string &t, const std::string &s, int n, REAL tol = 0., const std::string &c = "") {
        const TRecord *r = Find(t, s, n, tol, c);
        if (!r) throw std::logic_error("missing run " + t + " " + s);
        return *r;
    };
    std::cout << "\nNumbers quoted in the text of Sect. 6.2 (this work [article v0.6]):\n" << std::setprecision(3);
    const TReference &xb = Reference("xieB");
    std::cout << "  test B, closed form at eps_q = 10%: p' = " << std::fixed << std::setprecision(4) << xb.fP
              << " kPa, q = " << xb.fQ << " kPa (eta = " << std::setprecision(6) << xb.fEta << ")\n";
    std::cout << std::scientific << std::setprecision(2);
    std::cout << "  control of Xie et al., tolerance 1e-4: ME2(1) " << rec("xieB", "ME2(1)", 1, 1e-4, "xie").fV[0]
              << " [4.30e-05], RKDP5(4) " << rec("xieB", "RKDP5(4)", 1, 1e-4, "xie").fV[0] << " [2.87e-06]\n";
    std::cout << "  one increment: this work " << rec("xieB", "exact_n", 1).fV[0] << " [8.0e-02], secant G "
              << rec("xieB", "exact_secant", 1).fV[0] << " [8.1e-02], frozen K " << rec("xieB", "frozen_n", 1).fV[0]
              << " [4.9e-02], " << rec("xieB", "exact_n", 1).fWork << " [9] Newton iterations\n";
    std::cout << "  128 increments: " << rec("xieB", "exact_n", 128).fV[0] << " [2.3e-04] in "
              << rec("xieB", "exact_n", 128).fWork << " [511] iterations\n";
    std::cout << "  RKDP5(4): " << rec("xieB", "RKDP5(4)", 1, 1e-3, "sloan").fV[0] << " [2.2e-06] with "
              << rec("xieB", "RKDP5(4)", 1, 1e-3, "sloan").fWork << " [126] evaluations, "
              << rec("xieB", "RKDP5(4)", 1, 1e-4, "sloan").fV[0] << " [2.3e-07] with "
              << rec("xieB", "RKDP5(4)", 1, 1e-4, "sloan").fWork << " [150]; ME2(1): "
              << rec("xieB", "ME2(1)", 1, 1e-4, "sloan").fV[0] << " [3.2e-05] with "
              << rec("xieB", "ME2(1)", 1, 1e-4, "sloan").fWork << " [288]\n";
    std::cout << std::fixed << std::setprecision(2);
    const std::string dn = "kl_drained_ocr1";
    std::cout << "  drained NC, 10 increments, error in q: ME2(1) " << rec(dn, "ME2(1)", 10, 1e-4).fV[0] << " / "
              << rec(dn, "ME2(1)", 10, 1e-5).fV[0] << ", RKDP5(4) " << rec(dn, "RKDP5(4)", 10, 1e-4).fV[0] << " / "
              << rec(dn, "RKDP5(4)", 10, 1e-5).fV[0] << " (tol 1e-4 / 1e-5) [2.0]; BE " << rec(dn, "exact_n", 10).fV[0]
              << " [3.6] kPa\n";
    std::cout << "  drained NC, 100 increments: RKDP5(4) " << rec(dn, "RKDP5(4)", 100, 1e-4).fV[0] << " [0.05] kPa with "
              << rec(dn, "RKDP5(4)", 100, 1e-4).fWork << " [666] evaluations; BE " << rec(dn, "exact_n", 100).fV[0]
              << " [0.37] kPa with " << rec(dn, "exact_n", 100).fWork << " [520] iterations\n";
    const std::string un = "kl_undrained_ocr10";
    std::cout << "  undrained OCR 10, error in p': frozen K, 10 increments " << rec(un, "frozen_n", 10).fV[0]
              << " [28.7], exact " << rec(un, "exact_n", 10).fV[0] << " [1.7]; frozen K, 1000 increments "
              << rec(un, "frozen_n", 1000).fV[0] << " [0.38], exact with 50 " << rec(un, "exact_n", 50).fV[0]
              << " [0.11] kPa\n";
    const TRecord f10 = rec(dn, "frozen_n", 10), f20 = rec(dn, "frozen_n", 20);
    const TReference &rd = Reference(dn);
    std::cout << "  drained NC, frozen K: 10 increments " << (f10.fConverged ? "converged" : "no convergence")
              << " [no convergence]; 20 increments eps_v error " << 100. * f20.fV[1] << " [3.87] percentage points ("
              << 100. * (rd.fEpsV + f20.fV[1]) << "% [7.83%] against " << 100. * rd.fEpsV << "% [3.96%])\n";
    // secant shear modulus against this work: ratio of the errors (primary quantity of each test)
    REAL ncmax = 0., ocmin = 1e300, ocmax = -1e300;
    for (auto &r : fRecords) {
        if (r.fScheme != "exact_secant" || !r.fConverged) continue;
        const TRecord *t = Find(r.fTest, "exact_n", r.fN);
        if (!t || !t->fConverged) continue;
        const REAL ratio = std::fabs(r.fV[0]) / std::fabs(t->fV[0]);
        if (r.fTest == "xieB" || r.fTest.find("ocr1") + 4 == r.fTest.size()) ncmax = std::max(ncmax, std::fabs(ratio - 1.));
        else {
            ocmin = std::min(ocmin, 1. - ratio);
            ocmax = std::max(ocmax, 1. - ratio);
        }
    }
    std::cout << "  secant G against this work: NC tests, errors changed by at most " << 100. * ncmax
              << "% [less than 3%]; OC tests, errors reduced by " << 100. * ocmin << "% to " << 100. * ocmax
              << "% [2 to 20%]\n";
    const TRecord s10 = rec("kl_drained_ocr10", "exact_secant", 10);
    std::cout << "  secant G with sub-stepping, drained OCR 10, 10 increments: " << s10.fSubsteppedIncrements
              << " increment(s) needed sub-stepping, " << s10.fSubsteps << " sub-steps in all [first increment, two "
              << "sub-steps]\n";
    for (int scheme = 0; scheme < 2; ++scheme) {
        int npts = 0, nconv = 0, nspur = 0, itmin = 1000, itmax = 0;
        for (auto &pt : fSweep) {
            if (pt.fScheme != scheme) continue;
            npts++;
            if (!pt.fConverged) continue;
            if (pt.fP < 1e-6 * fKLP0) nspur++;
            else nconv++;
            itmin = std::min(itmin, pt.fIterations);
            itmax = std::max(itmax, pt.fIterations);
        }
        std::cout << "  first increment of the drained test, OCR 10, 10 increments, single step for " << npts
                  << " lateral strains in [0.010, 0.030]: " << (scheme ? "secant G" : "this work") << " converged at "
                  << nconv << ", spurious state p' = 0 at " << nspur << ", failed at " << npts - nconv - nspur
                  << " (" << itmin << " to " << itmax << " local iterations)\n";
    }
    std::cout << "  test B, this work, tolerance of the local Newton iterations (1e-2 to 1e-14): relative error / work\n";
    for (int n : {1, 4, 16, 64}) {
        std::cout << "    n = " << std::setw(2) << n << ":";
        for (auto &pt : fTolerance)
            if (pt.fN == n) {
                std::ostringstream o;
                if (pt.fConverged) o << std::scientific << std::setprecision(3) << pt.fError << " / " << pt.fWork;
                else o << "no convergence";
                std::cout << "  " << o.str();
            }
        std::cout << "\n";
    }
    int failures = 0, total = 0;
    for (auto &r : fRecords)
        if (r.fScheme == "exact_n" || r.fScheme == "library_exact") {
            total++;
            failures += !r.fConverged || r.fSubsteps != r.fN;
        }
    std::cout << "  this work (be_tensor and library): " << total - failures << " of " << total
              << " tests converged with a single step per increment [all]\n";
    for (auto &r : fRecords)
        if (r.fScheme == "exact_secant" && r.fSubsteppedIncrements > 0)
            std::cout << "  secant G: sub-stepping in " << r.fTest << ", n = " << r.fN << ": "
                      << r.fSubsteppedIncrements << " increment(s), " << r.fSubsteps << " sub-steps\n";
    // equilibrium of the drained tests (mixed control): every accepted increment must satisfy |sigma_r + p'0| <
    // 1e-10 p'0; the bisection safeguard of the driver accepts its last midpoint even without equilibrium
    int ndrained = 0, nbis = 0, nunbal = 0;
    REAL resmax = 0.;
    std::string bisruns;
    for (auto &r : fRecords) {
        if (r.fTest.rfind("kl_drained", 0) != 0 || !r.fConverged) continue;
        ndrained++;
        nbis += r.fBisectionIncrements;
        nunbal += r.fUnbalancedIncrements;
        resmax = std::max(resmax, r.fMaxResidual);
        if (r.fBisectionIncrements > 0)
            bisruns += " " + r.fTest + " " + r.fScheme + " n = " + std::to_string(r.fN) + ";";
    }
    std::cout << "  drained tests: " << ndrained << " runs with a solution, " << nbis
              << " increment(s) solved by the bisection safeguard," << (bisruns.empty() ? "" : bisruns) << " "
              << nunbal << " accepted without equilibrium; largest |sigma_r + p'0| of the accepted increments "
              << std::scientific << std::setprecision(2) << resmax << " kPa (tolerance " << 1e-10 * fKLP0
              << " kPa)\n";
    std::cout << std::defaultfloat << std::setprecision(6);
}

inline void IntegrationSchemes::PrintLibraryCheck() const {
    std::cout << "\nThe return mapping of this work in the library (TPZPlasticStepModifiedCamClay, spectral form) "
                 "against be_tensor:\n";
    for (const std::string lib : {"library_exact", "library_frozen"}) {
        const std::string be = lib == "library_exact" ? "exact_n" : "frozen_n";
        int nruns = 0, workdiff = 0, convdiff = 0;
        REAL dmax = 0., rel = 0.;
        for (auto &r : fRecords) {
            if (r.fScheme != lib) continue;
            const TRecord *t = Find(r.fTest, be, r.fN);
            nruns++;
            if (r.fConverged != t->fConverged) {
                convdiff++;
                continue;
            }
            if (!r.fConverged) continue;
            if (r.fWork != t->fWork) workdiff++;
            for (size_t i = 0; i < r.fPath.fRows.size(); ++i)
                for (int c = 1; c < 3; ++c) {
                    const REAL d = std::fabs(r.fPath.fRows[i][c] - t->fPath.fRows[i][c]);
                    dmax = std::max(dmax, d);
                    rel = std::max(rel, d / std::max(1., std::fabs(t->fPath.fRows[i][c])));
                }
        }
        std::cout << "  " << lib << " vs " << be << ": " << nruns << " runs, different convergence " << convdiff
                  << ", different work " << workdiff << ", largest |difference| of p', q along the paths "
                  << std::scientific << std::setprecision(2) << dmax << " kPa (relative " << rel << ")"
                  << std::defaultfloat << std::setprecision(6) << "\n";
    }
}

inline void IntegrationSchemes::PrintPythonComparison() const {
    std::cout << "\nComparison with the Python code (data_rivais.pkl, " << gIntegrationSchemesPythonSize
              << " results):\n";
    struct TStat {
        int fN = 0, fConv = 0, fWork = 0, fMissing = 0;
        REAL fMax[3] = {0., 0., 0.};
        std::string fWhere[3];
    };
    std::map<std::string, TStat> stats;
    std::vector<std::string> workdiff;
    for (int i = 0; i < gIntegrationSchemesPythonSize; ++i) {
        const TIntegrationSchemesPythonRow &py = gIntegrationSchemesPython[i];
        const bool isref = std::string(py.fScheme) == "reference";
        const bool rk = std::string(py.fControl) != "-";
        const std::string group = std::string(py.fTest) + (isref ? " reference" : (rk ? " RK" : " BE"));
        TStat &s = stats[group];
        s.fN++;
        REAL cpp[3];
        bool conv = true;
        long work = 0;
        if (isref) {
            const TReference &ref = Reference(py.fTest);
            const bool drained = std::string(py.fTest).rfind("kl_drained", 0) == 0;
            cpp[0] = drained ? ref.fQ : ref.fP;
            cpp[1] = drained ? ref.fEpsV : ref.fQ;
            cpp[2] = drained ? ref.fQPeak : std::numeric_limits<REAL>::quiet_NaN();
        } else {
            const TRecord *r = Find(py.fTest, py.fScheme, py.fN, py.fTol, rk ? py.fControl : "");
            if (!r) {
                s.fMissing++;
                continue;
            }
            conv = r->fConverged;
            work = r->fWork;
            for (int k = 0; k < 3; ++k) cpp[k] = r->fV[k];
        }
        if (int(conv) != py.fConverged) {
            s.fConv++;
            continue;
        }
        if (!conv) continue;
        if (!isref && work != py.fWork) {
            s.fWork++;
            std::ostringstream o;
            o << py.fTest << " " << py.fScheme << " " << py.fControl << " n = " << py.fN << " tol = " << py.fTol
              << ": " << work << " (Python " << py.fWork << ")";
            workdiff.push_back(o.str());
        }
        for (int k = 0; k < 3; ++k) {
            if (std::isnan(py.fV[k])) continue;
            const REAL d = std::fabs(cpp[k] - py.fV[k]);
            if (d > s.fMax[k] || std::isnan(d) || s.fWhere[k].empty()) {
                s.fMax[k] = d;
                std::ostringstream o;
                o << py.fScheme << " n=" << py.fN;
                if (py.fTol > 0.) o << " tol=" << py.fTol;
                s.fWhere[k] = o.str();
            }
        }
    }
    std::cout << "  largest |C++ - Python| per test and family (value 0, 1, 2 as in TRecord::fV), "
                 "results with different convergence / work / missing\n";
    for (auto &kv : stats) {
        const TStat &s = kv.second;
        std::cout << "  " << std::left << std::setw(30) << kv.first << std::right << std::setw(4) << s.fN << " results:";
        for (int k = 0; k < 3; ++k)
            if (!s.fWhere[k].empty())
                std::cout << " v" << k << " " << std::scientific << std::setprecision(2) << s.fMax[k] << " ("
                          << s.fWhere[k] << ")";
        std::cout << std::defaultfloat << std::setprecision(6) << "; " << s.fConv << " / " << s.fWork << " / "
                  << s.fMissing << "\n";
    }
    if (workdiff.empty()) std::cout << "  the work counts are identical in all the results\n";
    // results whose values differ by more than 1e-9 (relative to max(|value|, 1)): only the drained test with
    // OCR = 10 and the secant G with 10 increments, in which the single secant step of the first increment fails or
    // converges to spurious states depending on round-off (see RunSecantSweep and the README): the C++ run converges
    // to an equilibrated first increment, whereas the bisection of the Python run stopped at a jump of
    // sigma_r + p'0 and accepted a state with sigma_r + p'0 = -2.56 kPa
    for (int i = 0; i < gIntegrationSchemesPythonSize; ++i) {
        const TIntegrationSchemesPythonRow &py = gIntegrationSchemesPython[i];
        if (std::string(py.fScheme) == "reference" || !py.fConverged) continue;
        const bool rk = std::string(py.fControl) != "-";
        const TRecord *r = Find(py.fTest, py.fScheme, py.fN, py.fTol, rk ? py.fControl : "");
        if (!r || !r->fConverged) continue;
        for (int k = 0; k < 3; ++k) {
            if (std::isnan(py.fV[k])) continue;
            const REAL scale = std::max(std::fabs(py.fV[k]), REAL(1.)) * (std::string(py.fTest) == "xieB" ? 1e-3 : 1.);
            if (std::fabs(r->fV[k] - py.fV[k]) > 1e-9 * scale)
                std::cout << "  not reproduced to round-off: " << py.fTest << " " << py.fScheme << " n = " << py.fN
                          << ", value " << k << ": " << std::setprecision(10) << r->fV[k] << " (Python " << py.fV[k]
                          << ")" << std::setprecision(6) << "\n";
        }
    }
    for (auto &w : workdiff) std::cout << "  different work: " << w << "\n";
    if (const TRecord *r = Find("kl_drained_ocr10", "exact_secant", 10))
        std::cout << "  (kl_drained_ocr10 exact_secant n = 10: here every increment is in equilibrium, largest "
                     "|sigma_r + p'0| = "
                  << std::scientific << std::setprecision(1) << r->fMaxResidual << std::defaultfloat
                  << std::setprecision(6)
                  << " kPa; the first increment of the Python run ended its bisection at a jump of sigma_r + p'0, "
                     "with sigma_r + p'0 = -2.56 kPa, out of equilibrium: see the README)\n";
}

inline void IntegrationSchemes::RunAll() {
    const auto t0 = std::chrono::steady_clock::now();
    std::cout << "Integration schemes of the Modified Cam-Clay model at a material point (Sect. 6.2, Table 4, Fig. 6)\n";
    const TMaterial X = XieMaterial(), K = KLMaterial();
    std::cout << std::setprecision(10) << "test B of Xie et al.: M = " << X.fM << ", lambda = " << X.fLambda
              << ", kappa = " << X.fKappa << ", v0 = " << X.fV0 << ", nu = " << X.fNu << ", p'0 = p'c0 = " << fXieP0
              << " kPa\nK&L example 1: M = " << K.fM << ", lambda = " << K.fLambda << ", kappa = " << K.fKappa
              << ", v0 = " << K.fV0 << ", nu = " << K.fNu << ", p'0 = " << fKLP0 << " kPa\n"
              << std::setprecision(6);
    RunXieTestB();
    RunKLUndrained(1, 0.025);
    RunKLUndrained(10, 0.08);
    RunKLDrained(1);
    RunKLDrained(10);
    RunSecantSweep();
    RunToleranceStudy();
    fTotalTime = std::chrono::duration<REAL>(std::chrono::steady_clock::now() - t0).count();
    std::cout << "closed-form references:\n" << std::setprecision(12);
    for (auto &r : fReferences)
        std::cout << "  " << std::left << std::setw(20) << r.fTest << std::right << " eps_a = " << r.fEpsA
                  << ": eta = " << r.fEta << ", p' = " << r.fP << ", q = " << r.fQ << ", eps_v = " << r.fEpsV
                  << ", peak q = " << r.fQPeak << "\n";
    std::cout << std::setprecision(6);
    PrintTable4();
    PrintTextNumbers();
    PrintLibraryCheck();
    PrintPythonComparison();
    PostProcess();
    std::map<std::string, REAL> times;
    for (auto &r : fRecords) times[r.fScheme] += r.fTime;
    std::cout << "\ncomputing time per scheme (s):";
    for (auto &kv : times) std::cout << " " << kv.first << " " << std::setprecision(3) << kv.second;
    std::cout << std::setprecision(6) << "\n" << fRecords.size() << " runs in " << fTotalTime << " s\n";
}
