/**
 * @file YieldSurfaceProjection.h
 * @brief Figs. 1 and 2 of the article: the Modified Cam-Clay yield surface in principal stresses and in rotated
 * Haigh-Westergaard (RHW) space, and the closest-point projection in the meridian plane.
 *
 * C++ counterpart of fig_surface.py (functions fig_surface, cpp_linear, ellipse and fig_meridian) of the Python
 * transcription of the article "Return mapping for Modified Cam-Clay plasticity in rotated Haigh-Westergaard
 * space with consistent tangent operator and coupled u-p consolidation" (D. Lira Cecilio).
 */
#pragma once

#include "MCCPaperTools.h"
#include "TPZHWTools.h"
#include "TPZVTKGeoMesh.h"
#include <algorithm>
#include <array>
#include <random>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

/**
 * @ingroup mccpaper
 * @brief Yield surface of the Modified Cam-Clay model (Fig. 1) and closest-point projection in the meridian
 * plane (Fig. 2).
 *
 * This example has no finite element problem; it exercises the two classes that implement the constitutive
 * model of the article and keeps the structure of the other examples (data of the problem, geometric mesh,
 * "solution" and post-processing):
 *  - Fig. 1: the surface (12),
 *    \f$\Phi(\xi,\rho,a)=\frac{1}{b^2}\left(\frac{\xi}{\sqrt3}-p_t+a\right)^2+\frac{3\rho^2}{2M^2}-a^2=0\f$,
 *    for \f$M=1\f$, \f$p'_c=100\f$ kPa, \f$p_t=0\f$, \f$\omega=1\f$ (\f$a=50\f$ kPa) is sampled on the grid of
 *    fig_surface.py (121 values of \f$\xi\in[\sqrt3p_t,-\sqrt3p'_c]\f$ and 97 values of the Lode angle
 *    \f$\beta\in[0,2\pi]\f$). The principal stresses and the RHW coordinates
 *    \f$\sigma^*=\{\xi,\rho\cos\beta,\rho\sin\beta\}\f$ are obtained with the native TPZHWTools, the yield
 *    function is evaluated at every point with TPZYCModifiedCamClayRHW::YieldFunction (\f$\Phi\approx0\f$), and the
 *    surface is written as a closed TPZGeoMesh (triangles at the two apexes, quadrilaterals elsewhere, line
 *    elements on the critical state circle) exported with TPZVTKGeoMesh.
 *  - Fig. 2: closest-point projection with linear elasticity (\f$K=6\f$ MPa, \f$G=3\f$ MPa), \f$M=1\f$,
 *    \f$p_{c,n}=200\f$ kPa, \f$\lambda=0.2\f$, \f$\kappa=0.05\f$, \f$v_0=2\f$
 *    (\f$v_0/(\lambda-\kappa)=13.3\f$) for the trial states \f$(p',q)=(190,150)\f$ kPa (subcritical) and
 *    \f$(55,128)\f$ kPa (supercritical). The local problem (18) is solved with
 *    TPZYCModifiedCamClayRHW::ProjectHW and the result is verified through the complete stress update
 *    TPZPlasticStepModifiedCamClay::ApplyStrainComputeSigma (strain driven), first for the triaxial trial stress
 *    (projected state, \f$p_c\f$ and consistent tangent compared with camclay_hw.apply_strain) and then for 12
 *    Lode angles with the same trial invariants in a rotated frame (radial return).
 *
 * Sign conventions: the constitutive classes use tension positive, \f$p=I_1/3\f$, \f$\xi=\sqrt3\,p\f$,
 * \f$\rho=\sqrt{2J_2}\f$; the figures and the printout use the soil mechanics quantities \f$p'=-p\f$ and
 * \f$q=\sqrt{3/2}\,\rho\f$ (kPa).
 */
class YieldSurfaceProjection {
public:
    /** @brief Material ids of the geometric mesh of the yield surface (Fig. 1) */
    enum {
        ESubcritical = 1,   ///< part of the surface with \f$\bar p<0\f$ (\f$p'>p'_c/2\f$, compaction)
        ESupercritical = 2, ///< part of the surface with \f$\bar p\ge0\f$ (\f$p'\le p'_c/2\f$, dilation)
        ECriticalState = 3  ///< line elements on the critical state circle \f$\bar p=0\f$
    };

    /** @brief Space in which the yield surface mesh is drawn */
    enum ESpace {
        EPrincipalCompression = 0, ///< principal stresses with compression axes \f$(-\sigma_1,-\sigma_2,-\sigma_3)\f$ (Fig. 1a)
        ERHW = 1                   ///< rotated Haigh-Westergaard space \f$(\xi,\rho\cos\beta,\rho\sin\beta)\f$ (Fig. 1b)
    };

    /** @name Data of Fig. 1 */
    /** @{ */
    REAL fM1 = 1.;     ///< slope of the critical state line
    REAL fPc1 = 100.;  ///< preconsolidation pressure \f$p'_c\f$ (kPa)
    REAL fPt1 = 0.;    ///< tensile strength \f$p_t\f$ (kPa)
    REAL fOmega1 = 1.; ///< shape parameter \f$\omega\f$ of the subcritical region
    int fNXi = 121;    ///< number of values of \f$\xi\f$ of the grid (as in fig_surface.py)
    int fNBeta = 97;   ///< number of values of the Lode angle of the grid, both ends included
    /** @} */

    /** @name Data of Fig. 2 (cpp_linear of fig_surface.py) */
    /** @{ */
    REAL fK = 6000.;     ///< bulk modulus of the linear volumetric law (kPa)
    REAL fG = 3000.;     ///< constant shear modulus (kPa)
    REAL fM2 = 1.;       ///< slope of the critical state line
    REAL fLambda = 0.2;  ///< slope of the normal compression line
    REAL fKappa = 0.05;  ///< slope of the swelling line
    REAL fV0 = 2.;       ///< specific volume, \f$v_0/(\lambda-\kappa)=13.3\f$
    REAL fPcn = 200.;    ///< preconsolidation pressure at the start of the step \f$p_{c,n}\f$ (kPa)
    REAL fPStart = 100.; ///< isotropic converged state \f$p'_n\f$ of the strain driven stress update (kPa)
    int fNCurve = 400;   ///< number of points of the curves of Fig. 2 (ellipses and energy contour)
    int fNLode = 12;     ///< number of Lode angles \f$\beta=2\pi k/n\f$ of the radial return check
    /** @} */

    /** @brief Trial state of Fig. 2 and the reference values of the Python code (fig_surface.py, camclay_hw.py) */
    struct TMeridianCase {
        std::string fName;  ///< suffix of the output files
        std::string fTitle; ///< title of the panel of Fig. 2
        REAL fPTrial = 0.;  ///< trial mean effective stress \f$p'_{tr}\f$ (kPa)
        REAL fQTrial = 0.;  ///< trial deviatoric stress \f$q_{tr}\f$ (kPa)
        REAL fPRef = 0.;    ///< Python (cpp_linear): projected \f$p'\f$
        REAL fQRef = 0.;    ///< Python (cpp_linear): projected \f$q\f$
        REAL fARef = 0.;    ///< Python (cpp_linear): semi-axis \f$a\f$ at the end of the step
        REAL fDalRef = 0.;  ///< Python (cpp_linear): \f$\Delta\alpha\f$
        REAL fDgRef = 0.;   ///< Python (cpp_linear): \f$\Delta\gamma\f$
        int fItRef = 0;     ///< Python (cpp_linear and apply_strain): local Newton iterations
        REAL fPcRef = 0.;   ///< Python (apply_strain): \f$p_c\f$ at the end of the step
        REAL fAArticle = 0.; ///< value of \f$a\f$ quoted in the caption of Fig. 2
        /** @brief Python consistent tangent (apply_strain): D(xx,xx), D(xx,yy), D(xx,zz), D(zz,xx), D(zz,zz), D(xy,xy) */
        std::array<REAL, 6> fDRef{};
    };

    /** @brief Solution of the local problem in the meridian plane (reduced problem, with the full system as a check) */
    struct TMeridianResult {
        REAL fP = 0.;   ///< projected \f$p'\f$ (kPa)
        REAL fQ = 0.;   ///< projected \f$q\f$ (kPa)
        REAL fA = 0.;   ///< semi-axis at the end of the step (kPa)
        REAL fAn = 0.;  ///< semi-axis at the start of the step (kPa)
        REAL fDal = 0.; ///< hardening increment \f$\Delta\alpha=-\Delta\varepsilon^p_v\f$
        REAL fDg = 0.;  ///< plastic multiplier \f$\Delta\gamma\f$ recovered from the reduced solution
        REAL fTheta = 0.; ///< angle of the projected state on the meridian ellipse (0 at the compression apex)
        int fIter = 0;  ///< number of local Newton corrections of the reduced problem (ProjectReduced)
        bool fConverged = false; ///< true if ProjectReduced converged
        int fIterFull = 0;       ///< number of Newton corrections of the four-unknown system (ProjectHW)
        bool fConvergedFull = false; ///< true if ProjectHW converged
        REAL fDpFull = 0.;  ///< \f$|p'-p'_{full}|\f$ between the two solvers (kPa)
        REAL fDqFull = 0.;  ///< \f$|q-q_{full}|\f$ (kPa)
        REAL fDaFull = 0.;  ///< \f$|a-a_{full}|\f$ (kPa)
        REAL fDgFull = 0.;  ///< \f$|\Delta\gamma-\Delta\gamma_{full}|\f$
    };

    /** @brief Summary of the cross-check of the two local solvers on random trial states */
    struct TSolverCheck {
        int fNStates = 0;      ///< number of trial states outside the surface
        int fFailReduced = 0;  ///< failures of the reduced solver
        int fFailFull = 0;     ///< failures of the full solver
        REAL fDp = 0.;         ///< largest \f$|p'-p'_{full}|\f$ (kPa)
        REAL fDq = 0.;         ///< largest \f$|q-q_{full}|\f$ (kPa)
        REAL fDa = 0.;         ///< largest \f$|a-a_{full}|\f$ (kPa)
        REAL fDproj = 0.;      ///< largest entry of \f$|D_{proj}-D_{proj,full}|\f$
        REAL fDprojP = 0.;     ///< trial p' of the state of the largest difference in \f$D_{proj}\f$ (kPa)
        REAL fDprojQ = 0.;     ///< trial q of that state (kPa)
        int fNegDgReduced = 0; ///< states where the reduced solution has \f$\Delta\gamma<0\f$ (rejected by the stress update)
        int fNegDgFull = 0;    ///< states where the full solution has \f$\Delta\gamma<0\f$
        int fNotMinimum = 0;   ///< states where the distance of the reduced solution exceeds that of the full one by more than 1e-8
        REAL fMeanIterReduced = 0.; ///< mean Newton corrections of the reduced solver
        REAL fMeanIterFull = 0.;    ///< mean Newton corrections of the full solver
        int fMaxIterReduced = 0;    ///< largest number of corrections of the reduced solver
        int fMaxIterFull = 0;       ///< largest number of corrections of the full solver
    };

    /** @brief Result of the complete stress update for one trial state */
    struct TStressUpdateResult {
        TPZTensor<REAL> fSigmaTrial; ///< trial stress computed by TPZPlasticStepModifiedCamClay::TrialStress
        TPZTensor<REAL> fSigma;      ///< updated stress
        TPZFNMatrix<36, REAL> fDep;  ///< consistent tangent (engineering strains, Voigt XX, XY, XZ, YY, YZ, ZZ)
        REAL fPc = 0.;               ///< preconsolidation pressure at the end of the step
        int fType = 0;               ///< 0 elastic, 1 subcritical, 2 supercritical
        int fIter = 0;               ///< local Newton iterations
        bool fFailed = false;        ///< true if the local projection failed
    };

    /** @brief Largest deviations of the stress update over the Lode angles of a trial state (radial return check) */
    struct TLodeCheck {
        int fNTests = 0;    ///< number of trial states (Lode angles)
        int fFailures = 0;  ///< number of failed local projections
        REAL fDp = 0.;      ///< \f$\max|p'-p'_{ref}|\f$ (kPa)
        REAL fDq = 0.;      ///< \f$\max|q-q_{ref}|\f$ (kPa)
        REAL fDpc = 0.;     ///< \f$\max|p_c-p_{c,ref}|\f$ (kPa)
        REAL fDn = 0.;      ///< \f$\max\|n-n_{tr}\|\f$, unit deviatoric directions after and before the return
    };

    /** @name Fig. 1: yield surface */
    /** @{ */

    /** @brief Yield criterion of Fig. 1 (\f$\lambda\f$ and \f$\kappa\f$ do not affect the surface) */
    TPZYCModifiedCamClayRHW CreateSurfaceCriterion() const;

    /** @brief Semi-axis of the surface of Fig. 1, \f$a=(p'_c+p_t)/(1+\omega)\f$ (hardening law (13) with \f$\Delta\alpha=0\f$) */
    REAL SurfaceSemiAxis(const TPZYCModifiedCamClayRHW &yc) const;

    /** @brief Hydrostatic coordinate of the i-th section of the grid, \f$\xi_i=\sqrt3\,[p_t+(-p'_c-p_t)\,i/(n_\xi-1)]\f$ */
    REAL GridXi(int i) const;

    /** @brief Lode angle of the j-th meridian of the grid, \f$\beta_j=2\pi j/(n_\beta-1)\f$ */
    REAL GridBeta(int j) const;

    /**
     * @brief Deviatoric radius of the surface at the section \f$\xi\f$ (rho_s of fig_surface.py),
     * \f$\rho_s(\xi)=\sqrt{2/3}\,M\sqrt{\max(a^2-\bar p^2/b^2,0)}\f$ with \f$\bar p=\xi/\sqrt3-p_t+a\f$
     */
    REAL SurfaceRadius(const TPZYCModifiedCamClayRHW &yc, REAL a, REAL xi) const;

    /**
     * @brief Geometric mesh of the yield surface of Fig. 1 (a closed surface of revolution about the hydrostatic
     * axis): triangles at the apexes, quadrilaterals elsewhere (material ESubcritical or ESupercritical according to
     * the sign of \f$\bar p\f$ at the centre of the cell) and line elements on the critical state circle
     * (ECriticalState). The nodes are the points of the grid of fig_surface.py, with the coincident points merged
     * (one node at each apex, meridian \f$\beta=2\pi\f$ identified with \f$\beta=0\f$), so that the mesh is closed.
     * @param space principal stresses (compression axes) or RHW coordinates
     * @param[out] elData deviatoric stress \f$q\f$ at the centre of each element (kPa), written in the VTK file
     */
    TPZGeoMesh *CreateGeoMesh(ESpace space, TPZVec<REAL> &elData) const;

    /** @brief Fig. 1: CSV of the surface and of the critical state circle, VTK of the surface and the printout */
    void RunSurface();
    /** @} */

    /** @name Fig. 2: closest-point projection in the meridian plane */
    /** @{ */

    /** @brief Yield criterion of Fig. 2 with the linear volumetric law (\f$K_0=K\f$) */
    TPZYCModifiedCamClayRHW CreateMeridianCriterion() const;

    /** @brief Plastic step of Fig. 2 (linear volumetric law, constant G, \f$p_c=p_{c,n}\f$, \f$v_0=2\f$) */
    TPZPlasticStepModifiedCamClay CreatePlasticStep() const;

    /**
     * @brief The two trial states of Fig. 2 with the reference values of the article (caption) and of the Python
     * code: cpp_linear of fig_surface.py (local problem) and camclay_hw.apply_strain from \f$\sigma_n=-p'_nI\f$ with the
     * strain that gives the same trial stress (\f$p_c\f$ and consistent tangent)
     */
    std::vector<TMeridianCase> MeridianCases() const;

    /**
     * @brief Local problem for the trial invariants of a case (cpp_linear of fig_surface.py):
     * \f$\xi_{tr}=-\sqrt3\,p'_{tr}\f$, \f$\rho_{tr}=\sqrt{2/3}\,q_{tr}\f$, \f$b=1\f$, solved with the reduced
     * problem (two unknowns) and, as a check, with the four-unknown system
     */
    TMeridianResult ProjectMeridian(const TMeridianCase &c) const;

    /**
     * @brief Cross-check of the reduced and of the full local solvers on random trial states outside the surface
     * (both regions, including states close to the apexes and to the critical state line): projected state,
     * semi-axis, Jacobian of the projection and number of Newton corrections
     * @param yc criterion (the volumetric law and the elastic moduli of the trial are those of the arguments)
     * @param G shear modulus
     * @param v0 specific volume
     * @param pcn preconsolidation pressure
     * @param n number of trial states
     * @param seed seed of the random generator
     */
    static TSolverCheck CrossCheckSolvers(const TPZYCModifiedCamClayRHW &yc, REAL G, REAL v0, REAL pcn, int n,
                                          unsigned seed);

    /**
     * @brief Complete stress update (Algorithm 1) reaching the trial stress sigmaTrial from the isotropic converged
     * state \f$\sigma_n=-p'_nI\f$ with \f$\varepsilon_n=0\f$: the total strain is
     * \f$\varepsilon=\mathbb{C}^{-1}(\sigma_{tr}-\sigma_n)\f$ (TPZElasticResponse::ComputeStrain with the pair E,
     * nu of K and G), so that the elastic predictor of TPZPlasticStepModifiedCamClay returns sigmaTrial
     */
    TStressUpdateResult StressUpdate(const TPZTensor<REAL> &sigmaTrial) const;

    /**
     * @brief Radial return check: stress update of the trial states with the invariants of the case, Lode angles
     * \f$\beta_k=2\pi k/n\f$ and principal directions rotated by an arbitrary rotation; the projected \f$p'\f$,
     * \f$q\f$ and \f$p_c\f$ must be those of the meridian projection and the deviatoric direction must not change
     */
    TLodeCheck LodeAngleCheck(const TMeridianCase &c) const;

    /**
     * @brief Writes the curves of a panel of Fig. 2 (fig_meridian of fig_surface.py): ellipses at the start and at
     * the end of the step, energy-norm contour centred at the trial state through the projected state, and the
     * points (trial, projected, tip of the flow direction arrow)
     */
    void WriteMeridianFiles(const TMeridianCase &c, const TMeridianResult &r) const;

    /** @brief Fig. 2: both trial states, printout side by side with the reference values and CSV files */
    void RunProjection();
    /** @} */

    /** @brief Runs Figs. 1 and 2 */
    void RunAll();

    /** @name Auxiliary functions */
    /** @{ */

    /** @brief Relative difference \f$|x-r|/\max(|r|,10^{-300})\f$ */
    static REAL RelDiff(REAL x, REAL r);

    /** @brief Value i (0..n-1) of numpy.linspace(start, stop, n), computed as numpy does (the last value is exactly stop) */
    static REAL Linspace(REAL start, REAL stop, int n, int i);

    /** @brief Returns x with the sign of a zero removed (avoids "-0" in the CSV files) */
    static REAL NoNegativeZero(REAL x);

    /**
     * @brief Principal stresses of the HW cylindrical coordinates (TPZHWTools) and value of the yield function
     * (TPZYCModifiedCamClayRHW::YieldFunction) at them
     * @param yc yield criterion
     * @param pc preconsolidation pressure
     * @param cyl HW cylindrical coordinates \f$(\xi,\rho,\beta)\f$
     * @param[out] sigma principal stresses (tension positive)
     * @return \f$\Phi\f$
     */
    static REAL PhiAtHW(const TPZYCModifiedCamClayRHW &yc, REAL pc, const TPZVec<REAL> &cyl, TPZVec<REAL> &sigma);

    /** @brief Rotation matrix \f$R=R_z(a)R_y(b)R_x(c)\f$ */
    static TPZFNMatrix<9, REAL> Rotation(REAL a, REAL b, REAL c);

    /** @brief Symmetric tensor with principal values s and principal directions in the columns of R, \f$R\,\mathrm{diag}(s)R^T\f$ */
    static TPZTensor<REAL> TensorFromPrincipal(const TPZVec<REAL> &s, const TPZFMatrix<REAL> &R);

    /** @brief Unit deviatoric direction \f$n=s/\|s\|\f$ of a stress tensor (zero for an isotropic tensor) */
    static TPZTensor<REAL> DeviatoricDirection(const TPZTensor<REAL> &sig);
    /** @} */
};

// --------------------------------------------------------------------------------------------- Fig. 1

inline TPZYCModifiedCamClayRHW YieldSurfaceProjection::CreateSurfaceCriterion() const {
    TPZYCModifiedCamClayRHW yc;
    yc.SetUp(fM1, fLambda, fKappa, fPt1, fOmega1);
    return yc;
}

inline REAL YieldSurfaceProjection::SurfaceSemiAxis(const TPZYCModifiedCamClayRHW &yc) const {
    REAL a, H, pc;
    yc.Hardening(fPc1, 0., fV0, a, H, pc);
    return a;
}

inline REAL YieldSurfaceProjection::GridXi(int i) const {
    return std::sqrt(3.) * Linspace(fPt1, -fPc1, fNXi, i);
}

inline REAL YieldSurfaceProjection::GridBeta(int j) const {
    return Linspace(0., 2. * M_PI, fNBeta, j);
}

inline REAL YieldSurfaceProjection::SurfaceRadius(const TPZYCModifiedCamClayRHW &yc, REAL a, REAL xi) const {
    const REAL p = xi / std::sqrt(3.);
    const REAL pbar = p - yc.Pt() + a;
    const REAL b = yc.BFromP(p, a);
    return std::sqrt(2. / 3.) * yc.M() * std::sqrt(std::max(a * a - pbar * pbar / (b * b), REAL(0.)));
}

inline TPZGeoMesh *YieldSurfaceProjection::CreateGeoMesh(ESpace space, TPZVec<REAL> &elData) const {
    const TPZYCModifiedCamClayRHW yc = CreateSurfaceCriterion();
    const REAL a = SurfaceSemiAxis(yc);
    const int nring = fNXi - 2;    // sections with rho > 0 (the first and the last are the apexes)
    const int nmer = fNBeta - 1;   // distinct meridians (beta = 2 pi coincides with beta = 0)
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    gmesh->SetDimension(2);
    gmesh->NodeVec().Resize(2 + int64_t(nring) * nmer);
    // node coordinates from the HW cylindrical coordinates (xi, rho, beta) with the native TPZHWTools
    auto setnode = [&](int64_t index, REAL xi, REAL rho, REAL beta) {
        TPZManVector<REAL, 3> cyl = {xi, rho, beta}, x(3, 0.);
        if (space == ERHW) {
            TPZHWTools::FromHWCylToHWCart(cyl, x);
        } else {
            TPZHWTools::FromHWCylToPrincipal(cyl, x);
            for (int k = 0; k < 3; ++k) x[k] = -x[k];
        }
        gmesh->NodeVec()[index].Initialize(x, *gmesh);
    };
    const int64_t apex0 = 0, apex1 = 1;
    setnode(apex0, GridXi(0), 0., 0.);
    setnode(apex1, GridXi(fNXi - 1), 0., 0.);
    auto ring = [nmer](int i, int j) { return int64_t(2) + int64_t(i - 1) * nmer + (j % nmer); };
    for (int i = 1; i <= nring; ++i) {
        const REAL xi = GridXi(i);
        const REAL rho = SurfaceRadius(yc, a, xi);
        for (int j = 0; j < nmer; ++j) setnode(ring(i, j), xi, rho, GridBeta(j));
    }
    // elements: material from the sign of pbar at the centre of the cell, data = q at the centre of the cell
    std::vector<REAL> data;
    auto region = [&](int i, REAL &q) {
        const REAL xim = 0.5 * (GridXi(i) + GridXi(i + 1));
        q = std::sqrt(1.5) * SurfaceRadius(yc, a, xim);
        return (xim / std::sqrt(3.) - yc.Pt() + a >= 0.) ? ESupercritical : ESubcritical;
    };
    for (int i = 0; i < fNXi - 1; ++i) {
        REAL q;
        const int matid = region(i, q);
        for (int j = 0; j < nmer; ++j) {
            if (i == 0 || i == fNXi - 2) {
                TPZManVector<int64_t, 3> tri(3);
                if (i == 0) {
                    tri[0] = apex0; tri[1] = ring(1, j); tri[2] = ring(1, j + 1);
                } else {
                    tri[0] = ring(i, j); tri[1] = apex1; tri[2] = ring(i, j + 1);
                }
                new TPZGeoElRefPattern<pzgeom::TPZGeoTriangle>(tri, matid, *gmesh);
            } else {
                TPZManVector<int64_t, 4> quad(4);
                quad[0] = ring(i, j);
                quad[1] = ring(i + 1, j);
                quad[2] = ring(i + 1, j + 1);
                quad[3] = ring(i, j + 1);
                new TPZGeoElRefPattern<pzgeom::TPZGeoQuad>(quad, matid, *gmesh);
            }
            data.push_back(q);
        }
    }
    // critical state circle: the section of maximum rho, pbar = 0 (index (n-1)/2 of the grid when pt = 0)
    int ics = 1;
    for (int i = 1; i <= nring; ++i)
        if (SurfaceRadius(yc, a, GridXi(i)) > SurfaceRadius(yc, a, GridXi(ics))) ics = i;
    for (int j = 0; j < nmer; ++j) {
        TPZManVector<int64_t, 2> line(2);
        line[0] = ring(ics, j);
        line[1] = ring(ics, j + 1);
        new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(line, ECriticalState, *gmesh);
        data.push_back(std::sqrt(1.5) * SurfaceRadius(yc, a, GridXi(ics)));
    }
    gmesh->BuildConnectivity();
    elData.Resize(data.size());
    for (size_t k = 0; k < data.size(); ++k) elData[k] = data[k];
    return gmesh;
}

inline void YieldSurfaceProjection::RunSurface() {
    const REAL sq3 = std::sqrt(3.);
    const TPZYCModifiedCamClayRHW yc = CreateSurfaceCriterion();
    const REAL a = SurfaceSemiAxis(yc); // a = (pc + pt)/(1 + omega)
    const REAL pc = fPc1;

    // grid of fig_surface.py: principal stresses (HW formula), RHW coordinates and Phi
    std::vector<std::vector<REAL>> rows;
    REAL ximax = -1.e300, ximin = 1.e300, rhomax = 0., xirhomax = 0., phimax = 0.;
    for (int j = 0; j < fNBeta; ++j) {
        const REAL beta = GridBeta(j);
        for (int i = 0; i < fNXi; ++i) {
            const REAL xi = GridXi(i);
            const REAL rho = SurfaceRadius(yc, a, xi);
            TPZManVector<REAL, 3> cyl = {xi, rho, beta}, s(3), x(3);
            const REAL phi = PhiAtHW(yc, pc, cyl, s);
            TPZHWTools::FromHWCylToHWCart(cyl, x);
            phimax = std::max(phimax, std::fabs(phi) / (a * a));
            ximax = std::max(ximax, xi);
            ximin = std::min(ximin, xi);
            if (rho > rhomax) {
                rhomax = rho;
                xirhomax = xi;
            }
            rows.push_back({xi, beta, rho, NoNegativeZero(-xi / sq3), std::sqrt(1.5) * rho, s[0], s[1], s[2], x[0],
                            x[1], x[2], phi / (a * a)});
        }
    }
    mcc::WriteCSV("fig1_surface.csv", {"xi", "beta", "rho", "p_eff", "q", "sigma1", "sigma2", "sigma3", "sstar1",
                                       "sstar2", "sstar3", "phi_over_a2"}, rows);
    // critical state circle (xi_cs = sqrt3 (pt - a), as in fig_surface.py)
    const REAL xics = sq3 * (yc.Pt() - a);
    const REAL rcs = SurfaceRadius(yc, a, xics);
    rows.clear();
    for (int j = 0; j < fNBeta; ++j) {
        TPZManVector<REAL, 3> cyl = {xics, rcs, GridBeta(j)}, s(3), x(3);
        TPZHWTools::FromHWCylToPrincipal(cyl, s);
        TPZHWTools::FromHWCylToHWCart(cyl, x);
        rows.push_back({GridBeta(j), s[0], s[1], s[2], x[0], x[1], x[2]});
    }
    mcc::WriteCSV("fig1_critical_state_circle.csv", {"beta", "sigma1", "sigma2", "sigma3", "sstar1", "sstar2", "sstar3"},
                  rows);
    // VTK of the surface (geometric mesh written with the native TPZVTKGeoMesh)
    int64_t nel = 0, nnodes = 0;
    for (ESpace space : {EPrincipalCompression, ERHW}) {
        TPZManVector<REAL> elData;
        TPZGeoMesh *gmesh = CreateGeoMesh(space, elData);
        std::ofstream out(space == ERHW ? "fig1_surface_rhw.vtk" : "fig1_surface_principal.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out, elData);
        nel = gmesh->NElements();
        nnodes = gmesh->NNodes();
        delete gmesh;
    }

    // printout: characteristic values of the caption of Fig. 1 and checks of Phi with TPZYCModifiedCamClayRHW
    auto row = [](const std::string &name, REAL val, const std::string &ref) {
        std::cout << "  " << std::left << std::setw(46) << name << std::right << std::setw(14) << val
                  << std::setw(14) << ref << "\n";
    };
    std::cout << std::setprecision(6);
    std::cout << "\n=== Fig. 1: MCC yield surface (12), M = " << fM1 << ", p'c = " << fPc1 << " kPa, pt = " << fPt1
              << ", omega = " << fOmega1 << " ===\n";
    std::cout << "  " << std::left << std::setw(46) << "quantity" << std::right << std::setw(14) << "this code"
              << std::setw(14) << "article" << "\n";
    row("a = (p'c + pt)/(1 + omega) (kPa)", a, "50");
    row("length sqrt3 p'c (kPa), from the grid", ximax - ximin, "173");
    row("width 2 rho_max (kPa), from the grid", 2. * rhomax, "82");
    row("rho_max (kPa), from the grid", rhomax, "40.8");
    row("rho_max = sqrt(2/3) M a (kPa), closed form", std::sqrt(2. / 3.) * yc.M() * a, "40.8");
    row("p' of the section of rho_max (kPa)", -xirhomax / sq3, "p'c/2 = 50");
    row("q on the critical state circle (kPa)", std::sqrt(1.5) * rcs, "M a = 50");
    std::cout << "  Phi/a^2 (TPZYCModifiedCamClayRHW::YieldFunction) at points of the surface and inside it:\n";
    // points of the surface: (p'/p'c, Lode angle); the first and the last are the apexes
    const REAL surfacePoints[5][2] = {{0., 0.}, {0.25, M_PI / 6.}, {0.5, M_PI / 3.}, {0.8, 2.}, {1., 4.}};
    TPZManVector<REAL, 3> s(3);
    for (auto &pt : surfacePoints) {
        const REAL xi = -sq3 * (pt[0] * fPc1);
        const REAL rho = SurfaceRadius(yc, a, xi);
        TPZManVector<REAL, 3> cyl = {xi, rho, pt[1]};
        const REAL phi = PhiAtHW(yc, pc, cyl, s);
        std::ostringstream name;
        name << "surface: p' = " << pt[0] * fPc1 << ", q = " << std::sqrt(1.5) * rho << ", beta = " << pt[1];
        row(name.str(), phi / (a * a), "0");
    }
    {
        TPZManVector<REAL, 3> cyl = {sq3 * (yc.Pt() - a), 0., 0.};
        std::ostringstream name;
        name << "centre of the ellipse: p' = " << a - yc.Pt() << ", q = 0";
        row(name.str(), PhiAtHW(yc, pc, cyl, s) / (a * a), "-1");
    }
    std::cout << "  max |Phi|/a^2 on the " << fNXi << " x " << fNBeta << " grid of fig_surface.py = " << phimax << "\n";
    std::cout << "  files: fig1_surface.csv (" << fNXi * fNBeta << " points), fig1_critical_state_circle.csv, "
              << "fig1_surface_principal.vtk, fig1_surface_rhw.vtk (" << nel << " elements, " << nnodes << " nodes)\n";
}

// --------------------------------------------------------------------------------------------- Fig. 2

inline TPZYCModifiedCamClayRHW YieldSurfaceProjection::CreateMeridianCriterion() const {
    TPZYCModifiedCamClayRHW yc;
    yc.SetUp(fM2, fLambda, fKappa);
    yc.SetVolumetricLaw(TPZYCModifiedCamClayRHW::ELinear, fK);
    return yc;
}

inline TPZPlasticStepModifiedCamClay YieldSurfaceProjection::CreatePlasticStep() const {
    mcc::TPlastic model;
    model.SetModifiedCamClay(fM2, fLambda, fKappa);
    model.SetLinearVolumetric(fK);
    model.SetConstantShearModulus(fG);
    model.SetDefaultSpecificVolume(fV0);
    TPZPlasticState<REAL> st = model.GetState();
    st.m_eps_t.Zero();
    st.m_hardening = fPcn;
    model.SetState(st);
    return model;
}

inline std::vector<YieldSurfaceProjection::TMeridianCase> YieldSurfaceProjection::MeridianCases() const {
    std::vector<TMeridianCase> cases(2);
    cases[0].fName = "a_subcritical";
    cases[0].fTitle = "(a) subcritical region: compaction and hardening";
    cases[0].fPTrial = 190.;
    cases[0].fQTrial = 150.;
    cases[0].fPRef = 163.66180702920812;
    cases[0].fQRef = 88.99523828096234;
    cases[0].fARef = 106.02760701083203;
    cases[0].fDalRef = 0.004389698828465314;
    cases[0].fDgRef = 3.808241310771821e-05;
    cases[0].fItRef = 5;
    cases[0].fPcRef = 212.05521402166406;
    cases[0].fAArticle = 106.0;
    cases[0].fDRef = {6602.247265001305, 3042.437733762812, 2574.2963992271384, 2901.6525291634684,
                      2439.9231856983192, 1779.9047656192472};
    cases[1].fName = "b_supercritical";
    cases[1].fTitle = "(b) supercritical region: dilation and softening";
    cases[1].fPTrial = 55.;
    cases[1].fQTrial = 128.;
    cases[1].fPRef = 63.92816707217939;
    cases[1].fQRef = 91.91110414682453;
    cases[1].fARef = 98.0355153663691;
    cases[1].fDalRef = -0.0014880278453632317;
    cases[1].fDgRef = 2.181388937844696e-05;
    cases[1].fItRef = 5;
    cases[1].fPcRef = 196.0710307327382;
    cases[1].fAArticle = 98.0;
    cases[1].fDRef = {5364.935358600966, 1056.6023517185656, 4690.6376235275475, 4952.366698431286,
                      6782.816183652566, 2154.166503441201};
    return cases;
}

inline YieldSurfaceProjection::TMeridianResult YieldSurfaceProjection::ProjectMeridian(const TMeridianCase &c) const {
    const TPZYCModifiedCamClayRHW yc = CreateMeridianCriterion();
    TPZYCModifiedCamClayRHW::TTrial trial;
    trial.fXiTr = std::sqrt(3.) * (-c.fPTrial);
    trial.fRhoTr = std::sqrt(2. / 3.) * c.fQTrial;
    trial.fG = fG;
    trial.fV0 = fV0;
    trial.fPcn = fPcn;
    trial.fB = 1.;
    TMeridianResult r;
    TPZManVector<REAL, 4> X(4, 0.), Xfull(4, 0.);
    TPZManVector<REAL, 2> X2(2, 0.);
    r.fConverged = yc.ProjectReduced(trial, X2, r.fIter);
    yc.ReducedToFull(X2, trial, X);
    REAL H, pc;
    yc.Hardening(fPcn, X[2], fV0, r.fA, H, pc);
    yc.Hardening(fPcn, 0., fV0, r.fAn, H, pc);
    r.fP = -X[0] / std::sqrt(3.);
    r.fQ = X[1] * std::sqrt(1.5);
    r.fDal = X[2];
    r.fDg = X[3];
    r.fTheta = X2[0];
    // the same problem with the four-unknown system
    r.fConvergedFull = yc.ProjectHW(trial, Xfull, r.fIterFull);
    REAL afull;
    yc.Hardening(fPcn, Xfull[2], fV0, afull, H, pc);
    r.fDpFull = std::fabs(r.fP + Xfull[0] / std::sqrt(3.));
    r.fDqFull = std::fabs(r.fQ - Xfull[1] * std::sqrt(1.5));
    r.fDaFull = std::fabs(r.fA - afull);
    r.fDgFull = std::fabs(r.fDg - Xfull[3]);
    return r;
}

inline YieldSurfaceProjection::TSolverCheck
YieldSurfaceProjection::CrossCheckSolvers(const TPZYCModifiedCamClayRHW &yc, REAL G, REAL v0, REAL pcn, int n,
                                          unsigned seed) {
    const REAL sq3 = std::sqrt(3.);
    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<REAL> uni(0., 1.);
    REAL an, H, pc;
    yc.Hardening(pcn, 0., v0, an, H, pc);
    TSolverCheck chk;
    REAL sumRed = 0., sumFull = 0.;
    int nboth = 0;
    while (chk.fNStates < n) {
        // log-uniform p' in [0.02 p_c,n, 20 p_c,n] (compression) and q in [1e-6 a_n, 20 a_n]
        const REAL p = pcn * std::pow(10., -1.7 + 3. * uni(rng));
        const REAL q = an * std::pow(10., -6. + 7.3 * uni(rng));
        if (yc.PhiCC(-p, std::sqrt(2. / 3.) * q, an, yc.BFromP(-p, an)) <= 0.) continue;
        chk.fNStates++;
        TPZYCModifiedCamClayRHW::TTrial trial;
        trial.fXiTr = -sq3 * p;
        trial.fRhoTr = std::sqrt(2. / 3.) * q;
        trial.fG = G;
        trial.fV0 = v0;
        trial.fPcn = pcn;
        trial.fB = yc.BFromP(-p, an);
        TPZManVector<REAL, 4> X(4, 0.), Xf(4, 0.);
        TPZManVector<REAL, 2> X2(2, 0.);
        int itr = 0, itf = 0;
        const bool okr = yc.ProjectReduced(trial, X2, itr);
        const bool okf = yc.ProjectHW(trial, Xf, itf);
        if (!okr) chk.fFailReduced++;
        if (!okf) chk.fFailFull++;
        if (!okr || !okf) continue;
        yc.ReducedToFull(X2, trial, X);
        if (X[3] < -1.e-14) chk.fNegDgReduced++;
        if (Xf[3] < -1.e-14) chk.fNegDgFull++;
        // the four-unknown Newton iterations may converge to a stationary point of the distance that is not its
        // minimum (negative multiplier): only the admissible solutions are compared
        if (X[3] < -1.e-14 || Xf[3] < -1.e-14) continue;
        nboth++;
        sumRed += itr;
        sumFull += itf;
        chk.fMaxIterReduced = std::max(chk.fMaxIterReduced, itr);
        chk.fMaxIterFull = std::max(chk.fMaxIterFull, itf);
        REAL a, af;
        yc.Hardening(pcn, X[2], v0, a, H, pc);
        yc.Hardening(pcn, Xf[2], v0, af, H, pc);
        chk.fDp = std::max(chk.fDp, std::fabs(X[0] - Xf[0]) / sq3);
        chk.fDq = std::max(chk.fDq, std::fabs(X[1] - Xf[1]) * std::sqrt(1.5));
        chk.fDa = std::max(chk.fDa, std::fabs(a - af));
        // squared distance (linear law) or Bregman divergence of the porous law, both solutions
        auto dist2 = [&](const TPZVec<REAL> &Z) {
            REAL dv;
            if (yc.VolumetricLaw() == TPZYCModifiedCamClayRHW::ELinear) {
                dv = (Z[0] - trial.fXiTr) * (Z[0] - trial.fXiTr) / (6. * yc.K0());
            } else {
                dv = -yc.Kappa() / (sq3 * v0) * (Z[0] * std::log(Z[0] / trial.fXiTr) - Z[0] + trial.fXiTr);
            }
            return dv + (Z[1] - trial.fRhoTr) * (Z[1] - trial.fRhoTr) / (4. * G);
        };
        if (dist2(X) > dist2(Xf) * (1. + 1.e-8) + 1.e-14) chk.fNotMinimum++;
        // Jacobian of the projection for a triaxial direction n (sigma_zz axial)
        TPZManVector<REAL, 3> nvec(3, 0.);
        const bool isotropic = !(trial.fRhoTr > 1.e-14 * std::max(1., p));
        if (!isotropic) {
            nvec[0] = nvec[1] = 1. / std::sqrt(6.);
            nvec[2] = -2. / std::sqrt(6.);
        }
        TPZFNMatrix<9, REAL> Dr(3, 3, 0.), Df(3, 3, 0.);
        if (yc.GradProjectionReduced(X2, trial, nvec, isotropic, Dr) && yc.GradProjection(Xf, trial, nvec, isotropic, Df)) {
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    if (std::fabs(Dr(i, j) - Df(i, j)) > chk.fDproj) {
                        chk.fDproj = std::fabs(Dr(i, j) - Df(i, j));
                        chk.fDprojP = p;
                        chk.fDprojQ = q;
                    }
                }
            }
        }
    }
    if (nboth > 0) {
        chk.fMeanIterReduced = sumRed / nboth;
        chk.fMeanIterFull = sumFull / nboth;
    }
    return chk;
}

inline YieldSurfaceProjection::TStressUpdateResult
YieldSurfaceProjection::StressUpdate(const TPZTensor<REAL> &sigmaTrial) const {
    TPZPlasticStepModifiedCamClay model = CreatePlasticStep();
    const TPZTensor<REAL> sigman = mcc::IsotropicTensor(-fPStart);
    // strain increment that brings sigma_n to sigma_tr with the linear elastic law (E, nu equivalent to K, G)
    TPZElasticResponse er;
    er.SetEngineeringData(9. * fK * fG / (3. * fK + fG), (3. * fK - 2. * fG) / (2. * (3. * fK + fG)));
    TPZTensor<REAL> eps;
    er.ComputeStrain(sigmaTrial - sigman, eps);
    TStressUpdateResult r;
    REAL Ktr, G;
    model.TrialStress(eps, sigman, model.SpecificVolume(), r.fSigmaTrial, Ktr, G);
    r.fSigma = sigman;
    model.ApplyStrainComputeSigma(eps, r.fSigma, &r.fDep);
    r.fFailed = model.LastProjectionFailed();
    r.fPc = model.GetState().m_hardening;
    r.fType = model.GetState().m_m_type;
    r.fIter = model.LastNewtonIterations();
    return r;
}

inline YieldSurfaceProjection::TLodeCheck YieldSurfaceProjection::LodeAngleCheck(const TMeridianCase &c) const {
    const TPZFNMatrix<9, REAL> R = Rotation(0.3, 0.7, 1.1);
    TLodeCheck chk;
    for (int k = 0; k < fNLode; ++k) {
        TPZManVector<REAL, 3> cyl = {-std::sqrt(3.) * c.fPTrial, std::sqrt(2. / 3.) * c.fQTrial, 2. * M_PI * k / fNLode};
        TPZManVector<REAL, 3> s(3);
        TPZHWTools::FromHWCylToPrincipal(cyl, s);
        const TStressUpdateResult u = StressUpdate(TensorFromPrincipal(s, R));
        chk.fNTests++;
        if (u.fFailed) {
            chk.fFailures++;
            continue;
        }
        const TPZTensor<REAL> ntr = DeviatoricDirection(u.fSigmaTrial), n = DeviatoricDirection(u.fSigma);
        chk.fDp = std::max(chk.fDp, std::fabs(mcc::MeanEffectiveStress(u.fSigma) - c.fPRef));
        chk.fDq = std::max(chk.fDq, std::fabs(mcc::DeviatoricStress(u.fSigma) - c.fQRef));
        chk.fDpc = std::max(chk.fDpc, std::fabs(u.fPc - c.fPcRef));
        chk.fDn = std::max(chk.fDn, (n - ntr).Norm());
    }
    return chk;
}

inline void YieldSurfaceProjection::WriteMeridianFiles(const TMeridianCase &c, const TMeridianResult &r) const {
    // ellipses p' in [0, 2a] (ellipse of fig_surface.py) and the energy-norm contour through the projected state,
    // (p' - p'_tr)^2/K + (q - q_tr)^2/(3G) = d^2
    const REAL d2 = (c.fPTrial - r.fP) * (c.fPTrial - r.fP) / fK + (c.fQTrial - r.fQ) * (c.fQTrial - r.fQ) / (3. * fG);
    std::vector<std::vector<REAL>> rows;
    for (int k = 0; k < fNCurve; ++k) {
        const REAL t = Linspace(0., M_PI, fNCurve, k);
        const REAL th = Linspace(0., 2. * M_PI, fNCurve, k);
        std::vector<REAL> rw;
        for (REAL a : {r.fAn, r.fA}) {
            const REAL p = -a + a * std::cos(t);
            rw.push_back(NoNegativeZero(-p));
            rw.push_back(fM2 * std::sqrt(std::max(a * a - (p + a) * (p + a), REAL(0.))));
        }
        rw.push_back(c.fPTrial + std::sqrt(d2 * fK) * std::cos(th));
        rw.push_back(c.fQTrial + std::sqrt(3. * fG * d2) * std::sin(th));
        rows.push_back(rw);
    }
    mcc::WriteCSV("fig2" + c.fName + "_curves.csv", {"p_start", "q_start", "p_end", "q_end", "p_energy", "q_energy"},
                  rows);
    // points: 0 trial, 1 projected, 2 tip of the arrow of the flow direction (length 28 kPa, as in the figure)
    REAL n0 = 2. * (r.fP - r.fA), n1 = 2. * r.fQ / (fM2 * fM2);
    const REAL nn = std::sqrt(n0 * n0 + n1 * n1);
    n0 /= nn;
    n1 /= nn;
    mcc::WriteCSV("fig2" + c.fName + "_points.csv", {"point", "p", "q"},
                  {{0., c.fPTrial, c.fQTrial}, {1., r.fP, r.fQ}, {2., r.fP + 28. * n0, r.fQ + 28. * n1}});
}

inline void YieldSurfaceProjection::RunProjection() {
    std::cout << "\n=== Fig. 2: closest-point projection in the meridian plane, linear elasticity K = " << fK
              << " kPa, G = " << fG << " kPa, M = " << fM2 << ", p_c,n = " << fPcn
              << " kPa, v0/(lambda - kappa) = " << std::setprecision(3) << fV0 / (fLambda - fKappa) << " ===\n";
    // critical state line q = M p' (two points, as in fig_meridian)
    mcc::WriteCSV("fig2_critical_state_line.csv", {"p", "q"}, {{0., 0.}, {240., fM2 * 240.}});

    auto row = [](const std::string &name, REAL val, REAL ref, const std::string &art) {
        std::cout << "  " << std::left << std::setw(22) << name << std::right << std::setprecision(15)
                  << std::setw(24) << val << std::setw(24) << ref << std::setprecision(1) << std::scientific
                  << std::setw(11) << RelDiff(val, ref) << std::defaultfloat << std::setw(10) << art << "\n";
    };
    for (const TMeridianCase &c : MeridianCases()) {
        const TMeridianResult r = ProjectMeridian(c);
        WriteMeridianFiles(c, r);
        std::cout << std::defaultfloat << std::setprecision(6) << "\n" << c.fTitle << ": trial p' = " << c.fPTrial
                  << " kPa, q = " << c.fQTrial << " kPa\n";
        std::cout << "  reduced local problem (theta, Delta alpha), TPZYCModifiedCamClayRHW::ProjectReduced"
                  << (r.fConverged ? "" : " NOT CONVERGED") << "\n";
        std::cout << "  " << std::left << std::setw(22) << "quantity" << std::right << std::setw(24) << "this code"
                  << std::setw(24) << "Python" << std::setw(11) << "rel.diff" << std::setw(10) << "article" << "\n";
        std::ostringstream art;
        art << std::fixed << std::setprecision(1) << c.fAArticle;
        row("p' (kPa)", r.fP, c.fPRef, "");
        row("q (kPa)", r.fQ, c.fQRef, "");
        row("a (kPa)", r.fA, c.fARef, art.str());
        row("a_n (kPa)", r.fAn, 100., "100");
        row("Delta alpha", r.fDal, c.fDalRef, c.fDalRef > 0 ? "> 0" : "< 0");
        row("Delta gamma", r.fDg, c.fDgRef, "");
        std::cout << "  " << std::left << std::setw(22) << "theta (rad)" << std::right << std::setprecision(15)
                  << std::setw(24) << r.fTheta << std::defaultfloat << "\n";
        std::cout << "  " << std::left << std::setw(22) << "Newton iterations" << std::right << std::setw(24) << r.fIter
                  << std::setw(24) << c.fItRef << "   (Python: four unknowns)\n";
        std::cout << std::scientific << std::setprecision(1) << "  four-unknown system (ProjectHW): " << r.fIterFull
                  << " iterations" << (r.fConvergedFull ? "" : " NOT CONVERGED") << ", |p' - p'_full| = " << r.fDpFull
                  << " kPa, |q - q_full| = " << r.fDqFull << " kPa, |a - a_full| = " << r.fDaFull
                  << " kPa, |Delta gamma - Delta gamma_full| = " << r.fDgFull << std::defaultfloat << "\n";
        std::cout << "  surface " << (r.fA > r.fAn ? "expands" : "contracts") << ": a = " << std::fixed
                  << std::setprecision(1) << r.fA << (r.fA > r.fAn ? " > " : " < ") << "a_n = " << std::setprecision(0)
                  << r.fAn << std::defaultfloat << " (article: a = " << art.str() << " kPa)\n";
        // tangency of the energy contour and the end-of-step surface at the projected state
        const REAL ny0 = 2. * (r.fP - r.fA), ny1 = 2. * r.fQ / (fM2 * fM2);
        const REAL ne0 = (c.fPTrial - r.fP) / fK, ne1 = (c.fQTrial - r.fQ) / (3. * fG);
        const REAL sine = (ny0 * ne1 - ny1 * ne0) / (std::hypot(ny0, ny1) * std::hypot(ne0, ne1));
        std::cout << std::setprecision(2) << std::scientific
                  << "  energy-norm contour tangent to Phi = 0 at the projected state: sin(angle between normals) = "
                  << sine << std::defaultfloat << "\n";

        // complete stress update with the triaxial trial stress (sigma_zz axial); p' and q are compared with the
        // local projection (camclay_hw.apply_strain gives the same values to the last bit)
        TPZTensor<REAL> sigtr;
        sigtr.XX() = sigtr.YY() = -(c.fPTrial - c.fQTrial / 3.);
        sigtr.ZZ() = -(c.fPTrial + 2. * c.fQTrial / 3.);
        const TStressUpdateResult su = StressUpdate(sigtr);
        std::cout << "  stress update TPZPlasticStepModifiedCamClay::ApplyStrainComputeSigma (from sigma_n = -"
                  << std::setprecision(4) << fPStart << " I kPa, eps_n = 0):\n";
        std::cout << std::setprecision(15) << "    trial p' = " << mcc::MeanEffectiveStress(su.fSigmaTrial)
                  << ", q = " << mcc::DeviatoricStress(su.fSigmaTrial) << (su.fFailed ? "  PROJECTION FAILED" : "")
                  << "\n";
        row("  p' (kPa)", mcc::MeanEffectiveStress(su.fSigma), c.fPRef, "");
        row("  q (kPa)", mcc::DeviatoricStress(su.fSigma), c.fQRef, "");
        row("  p_c (kPa)", su.fPc, c.fPcRef, "");
        std::cout << "    type " << su.fType << " (" << (su.fType == 1 ? "subcritical" : su.fType == 2 ? "supercritical" : "elastic")
                  << "), local iterations " << su.fIter << " (Python " << c.fItRef << ")\n";
        const int idx[6][2] = {{_XX_, _XX_}, {_XX_, _YY_}, {_XX_, _ZZ_}, {_ZZ_, _XX_}, {_ZZ_, _ZZ_}, {_XY_, _XY_}};
        const char *names[6] = {"D(xx,xx)", "D(xx,yy)", "D(xx,zz)", "D(zz,xx)", "D(zz,zz)", "D(xy,xy)"};
        std::cout << "    consistent tangent (kPa), non-symmetric:\n";
        for (int k = 0; k < 6; ++k) row(std::string("    ") + names[k], su.fDep(idx[k][0], idx[k][1]), c.fDRef[k], "");

        // invariance with respect to the Lode angle and to the orientation of the principal axes
        const TLodeCheck chk = LodeAngleCheck(c);
        std::cout << std::scientific << std::setprecision(1) << "    " << chk.fNTests
                  << " Lode angles beta = 2 pi k/" << fNLode << " in a rotated frame (" << chk.fFailures
                  << " failures): max |p' - p'_ref| = " << chk.fDp << " kPa, max |q - q_ref| = " << chk.fDq
                  << " kPa, max |p_c - p_c,ref| = " << chk.fDpc << " kPa, max |n - n_tr| = " << chk.fDn
                  << std::defaultfloat << "\n";
        std::cout << "  files: fig2" << c.fName << "_curves.csv, fig2" << c.fName << "_points.csv\n";
    }
    std::cout << "  file: fig2_critical_state_line.csv\n";

    // cross-check of the two local solvers on random trial states: linear elasticity of Fig. 2 and the porous
    // elasticity of the Abaqus clay (M = 1, lambda = 0.174, kappa = 0.026, v0 = 2.08, nu = 0.3 at p'_n = 100 kPa)
    std::cout << "\n  cross-check of the reduced (theta, Delta alpha) and of the four-unknown local solvers on random "
                 "trial states outside the surface\n";
    auto report = [](const std::string &name, const TSolverCheck &chk) {
        std::cout << std::scientific << std::setprecision(1) << "  " << name << ": " << chk.fNStates << " states, failures "
                  << chk.fFailReduced << " (reduced) " << chk.fFailFull << " (full); max |p' - p'_full| = " << chk.fDp
                  << " kPa, |q - q_full| = " << chk.fDq << " kPa, |a - a_full| = " << chk.fDa
                  << " kPa, |D_proj - D_proj,full| = " << chk.fDproj << " (at p' = " << std::defaultfloat
                  << std::setprecision(4) << chk.fDprojP << ", q = " << chk.fDprojQ << ")" << std::setprecision(3)
                  << "; Newton corrections: reduced " << chk.fMeanIterReduced << " (max " << chk.fMaxIterReduced
                  << "), full " << chk.fMeanIterFull << " (max " << chk.fMaxIterFull << "); solutions with Delta gamma < 0: "
                  << chk.fNegDgReduced << " (reduced), " << chk.fNegDgFull << " (full); reduced solution farther than the full one: "
                  << chk.fNotMinimum << "\n";
    };
    report("linear elasticity, K = 6 MPa, G = 3 MPa, p_c,n = 200 kPa", CrossCheckSolvers(CreateMeridianCriterion(), fG, fV0, fPcn, 2000, 2026));
    TPZYCModifiedCamClayRHW ycp;
    ycp.SetUp(1., 0.174, 0.026);
    const REAL v0p = 2.08, Kp = v0p * 100. / 0.026, Gp = 3. * (1. - 0.6) / (2. * 1.3) * Kp;
    report("porous elasticity, Abaqus clay, p'_n = 100 kPa, p_c,n = 116.6 kPa", CrossCheckSolvers(ycp, Gp, v0p, 116.6, 2000, 2026));
}

inline void YieldSurfaceProjection::RunAll() {
    RunSurface();
    RunProjection();
}

// --------------------------------------------------------------------------------------------- auxiliary functions

inline REAL YieldSurfaceProjection::RelDiff(REAL x, REAL r) {
    return std::fabs(x - r) / std::max(std::fabs(r), REAL(1.e-300));
}

inline REAL YieldSurfaceProjection::Linspace(REAL start, REAL stop, int n, int i) {
    if (i == n - 1) return stop;
    return i * ((stop - start) / (n - 1)) + start;
}

inline REAL YieldSurfaceProjection::NoNegativeZero(REAL x) { return x + 0.; }

inline REAL YieldSurfaceProjection::PhiAtHW(const TPZYCModifiedCamClayRHW &yc, REAL pc, const TPZVec<REAL> &cyl,
                                            TPZVec<REAL> &sigma) {
    TPZHWTools::FromHWCylToPrincipal(cyl, sigma);
    TPZManVector<STATE, 1> phi(1, 0.);
    yc.YieldFunction(sigma, pc, phi);
    return phi[0];
}

inline TPZFNMatrix<9, REAL> YieldSurfaceProjection::Rotation(REAL a, REAL b, REAL c) {
    TPZFNMatrix<9, REAL> Rz(3, 3, 0.), Ry(3, 3, 0.), Rx(3, 3, 0.), tmp, R;
    Rz(0, 0) = std::cos(a); Rz(0, 1) = -std::sin(a); Rz(1, 0) = std::sin(a); Rz(1, 1) = std::cos(a); Rz(2, 2) = 1.;
    Ry(0, 0) = std::cos(b); Ry(0, 2) = std::sin(b); Ry(2, 0) = -std::sin(b); Ry(2, 2) = std::cos(b); Ry(1, 1) = 1.;
    Rx(1, 1) = std::cos(c); Rx(1, 2) = -std::sin(c); Rx(2, 1) = std::sin(c); Rx(2, 2) = std::cos(c); Rx(0, 0) = 1.;
    Rz.Multiply(Ry, tmp);
    tmp.Multiply(Rx, R);
    return R;
}

inline TPZTensor<REAL> YieldSurfaceProjection::TensorFromPrincipal(const TPZVec<REAL> &s, const TPZFMatrix<REAL> &R) {
    TPZFNMatrix<9, REAL> D(3, 3, 0.), RD, Rt, RDRt;
    for (int i = 0; i < 3; ++i) D(i, i) = s[i];
    R.Multiply(D, RD);
    R.Transpose(&Rt);
    RD.Multiply(Rt, RDRt);
    TPZTensor<REAL> sig;
    for (int r = 0; r < 3; ++r)
        for (int c = r; c < 3; ++c) sig(r, c) = 0.5 * (RDRt(r, c) + RDRt(c, r));
    return sig;
}

inline TPZTensor<REAL> YieldSurfaceProjection::DeviatoricDirection(const TPZTensor<REAL> &sig) {
    TPZTensor<REAL> s;
    sig.S(s);
    const REAL norm = s.Norm();
    if (norm > 0.) s *= (1. / norm);
    return s;
}
