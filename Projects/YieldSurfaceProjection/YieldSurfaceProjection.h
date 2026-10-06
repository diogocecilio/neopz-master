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
#include <iostream>
#include <string>

/**
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
 *    TPZPlasticStepModifiedCamClay::ApplyStrainComputeSigma (strain driven, for 13 orientations and Lode angles
 *    of the same trial invariants) and against the consistent tangent of the Python code.
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
    /** @} */

    /** @brief Trial state of Fig. 2 and the reference values of the Python code (fig_surface.py, camclay_hw.py) */
    struct TMeridianCase {
        std::string fName;  ///< prefix of the output files
        std::string fTitle; ///< title of the panel of Fig. 2
        REAL fPTrial = 0.;  ///< trial mean effective stress \f$p'_{tr}\f$ (kPa)
        REAL fQTrial = 0.;  ///< trial deviatoric stress \f$q_{tr}\f$ (kPa)
        REAL fPRef = 0.;    ///< Python: projected \f$p'\f$
        REAL fQRef = 0.;    ///< Python: projected \f$q\f$
        REAL fARef = 0.;    ///< Python: semi-axis \f$a\f$ at the end of the step
        REAL fDalRef = 0.;  ///< Python: \f$\Delta\alpha\f$
        REAL fDgRef = 0.;   ///< Python: \f$\Delta\gamma\f$
        int fItRef = 0;     ///< Python: local Newton iterations
        REAL fPcRef = 0.;   ///< Python: \f$p_c\f$ returned by apply_strain
        REAL fAArticle = 0.; ///< value of \f$a\f$ quoted in the caption of Fig. 2
        /** @brief Python consistent tangent (apply_strain): D(xx,xx), D(xx,yy), D(xx,zz), D(zz,xx), D(zz,zz), D(xy,xy) */
        std::array<REAL, 6> fDRef{};
    };

    /** @brief Solution of the local problem (18) in the meridian plane */
    struct TMeridianResult {
        REAL fP = 0.;   ///< projected \f$p'\f$ (kPa)
        REAL fQ = 0.;   ///< projected \f$q\f$ (kPa)
        REAL fA = 0.;   ///< semi-axis at the end of the step (kPa)
        REAL fAn = 0.;  ///< semi-axis at the start of the step (kPa)
        REAL fDal = 0.; ///< hardening increment \f$\Delta\alpha=-\Delta\varepsilon^p_v\f$
        REAL fDg = 0.;  ///< plastic multiplier \f$\Delta\gamma\f$
        int fIter = 0;  ///< number of local Newton corrections
        bool fConverged = false; ///< true if ProjectHW converged
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

    /** @name Fig. 1: yield surface */
    /** @{ */

    /** @brief Yield criterion of Fig. 1 (\f$\lambda\f$ and \f$\kappa\f$ do not affect the surface) */
    TPZYCModifiedCamClayRHW CreateSurfaceCriterion() const;

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
     * @brief Local problem (18) for the trial invariants of a case (cpp_linear of fig_surface.py):
     * \f$\xi_{tr}=-\sqrt3\,p'_{tr}\f$, \f$\rho_{tr}=\sqrt{2/3}\,q_{tr}\f$, \f$b=1\f$
     */
    TMeridianResult ProjectMeridian(const TMeridianCase &c) const;

    /**
     * @brief Complete stress update (Algorithm 1) reaching the trial stress sigmaTrial from the isotropic converged
     * state \f$\sigma_n=-p'_nI\f$ with \f$\varepsilon_n=0\f$: the total strain is
     * \f$\varepsilon=\mathbb{C}^{-1}(\sigma_{tr}-\sigma_n)\f$ (TPZElasticResponse::ComputeStrain with the pair E,
     * nu of K and G), so that the elastic predictor of TPZPlasticStepModifiedCamClay returns sigmaTrial
     */
    TStressUpdateResult StressUpdate(const TPZTensor<REAL> &sigmaTrial) const;

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
    void Run();
};

// --------------------------------------------------------------------------------------------- helpers

namespace ysp {
/** @brief Relative difference \f$|x-r|/\max(|r|,10^{-300})\f$ */
inline REAL RelDiff(REAL x, REAL r) { return std::fabs(x - r) / std::max(std::fabs(r), REAL(1.e-300)); }

/** @brief Values 0..n-1 of numpy.linspace(start, stop, n) (the last value is exactly stop) */
inline REAL Linspace(REAL start, REAL stop, int n, int i) {
    if (i == n - 1) return stop;
    return i * ((stop - start) / (n - 1)) + start;
}

/** @brief Returns x with the sign of a zero removed (avoids "-0" in the CSV files) */
inline REAL NoNegativeZero(REAL x) { return x + 0.; }

/** @brief Rotation matrix \f$R=R_z(a)R_y(b)R_x(c)\f$ */
inline TPZFNMatrix<9, REAL> Rotation(REAL a, REAL b, REAL c) {
    TPZFNMatrix<9, REAL> Rz(3, 3, 0.), Ry(3, 3, 0.), Rx(3, 3, 0.), tmp, R;
    Rz(0, 0) = std::cos(a); Rz(0, 1) = -std::sin(a); Rz(1, 0) = std::sin(a); Rz(1, 1) = std::cos(a); Rz(2, 2) = 1.;
    Ry(0, 0) = std::cos(b); Ry(0, 2) = std::sin(b); Ry(2, 0) = -std::sin(b); Ry(2, 2) = std::cos(b); Ry(1, 1) = 1.;
    Rx(1, 1) = std::cos(c); Rx(1, 2) = -std::sin(c); Rx(2, 1) = std::sin(c); Rx(2, 2) = std::cos(c); Rx(0, 0) = 1.;
    Rz.Multiply(Ry, tmp);
    tmp.Multiply(Rx, R);
    return R;
}

/** @brief Symmetric tensor with principal values s and principal directions in the columns of R, \f$R\,\mathrm{diag}(s)R^T\f$ */
inline TPZTensor<REAL> TensorFromPrincipal(const TPZVec<REAL> &s, const TPZFMatrix<REAL> &R) {
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

/** @brief Unit deviatoric direction \f$n=s/\|s\|\f$ of a stress tensor (zero for an isotropic tensor) */
inline TPZTensor<REAL> DeviatoricDirection(const TPZTensor<REAL> &sig) {
    TPZTensor<REAL> s;
    sig.S(s);
    const REAL norm = s.Norm();
    if (norm > 0.) s *= (1. / norm);
    return s;
}
} // namespace ysp

// --------------------------------------------------------------------------------------------- Fig. 1

inline TPZYCModifiedCamClayRHW YieldSurfaceProjection::CreateSurfaceCriterion() const {
    TPZYCModifiedCamClayRHW yc;
    yc.SetUp(fM1, 0.2, 0.05, fPt1, fOmega1);
    return yc;
}

inline REAL YieldSurfaceProjection::GridXi(int i) const {
    return std::sqrt(3.) * ysp::Linspace(fPt1, -fPc1, fNXi, i);
}

inline REAL YieldSurfaceProjection::GridBeta(int j) const {
    return ysp::Linspace(0., 2. * M_PI, fNBeta, j);
}

inline REAL YieldSurfaceProjection::SurfaceRadius(const TPZYCModifiedCamClayRHW &yc, REAL a, REAL xi) const {
    const REAL p = xi / std::sqrt(3.);
    const REAL pbar = p - yc.Pt() + a;
    const REAL b = yc.BFromP(p, a);
    return std::sqrt(2. / 3.) * yc.M() * std::sqrt(std::max(a * a - pbar * pbar / (b * b), REAL(0.)));
}

inline TPZGeoMesh *YieldSurfaceProjection::CreateGeoMesh(ESpace space, TPZVec<REAL> &elData) const {
    const TPZYCModifiedCamClayRHW yc = CreateSurfaceCriterion();
    REAL a, H, pc;
    yc.Hardening(fPc1, 0., 2., a, H, pc);
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
    REAL a, H, pc;
    yc.Hardening(fPc1, 0., 2., a, H, pc); // a = (pc + pt)/(1 + omega)

    // grid of fig_surface.py: principal stresses (HW formula), RHW coordinates and Phi
    std::vector<std::vector<REAL>> rows;
    REAL ximax = -1.e300, ximin = 1.e300, rhomax = 0., xirhomax = 0., phimax = 0.;
    for (int j = 0; j < fNBeta; ++j) {
        const REAL beta = GridBeta(j);
        for (int i = 0; i < fNXi; ++i) {
            const REAL xi = GridXi(i);
            const REAL rho = SurfaceRadius(yc, a, xi);
            TPZManVector<REAL, 3> cyl = {xi, rho, beta}, s(3), x(3);
            TPZHWTools::FromHWCylToPrincipal(cyl, s);
            TPZHWTools::FromHWCylToHWCart(cyl, x);
            TPZManVector<STATE, 1> phi(1, 0.);
            yc.YieldFunction(s, pc, phi);
            phimax = std::max(phimax, std::fabs(phi[0]) / (a * a));
            ximax = std::max(ximax, xi);
            ximin = std::min(ximin, xi);
            if (rho > rhomax) {
                rhomax = rho;
                xirhomax = xi;
            }
            rows.push_back({xi, beta, rho, ysp::NoNegativeZero(-xi / sq3), std::sqrt(1.5) * rho, s[0], s[1], s[2], x[0], x[1], x[2],
                            phi[0] / (a * a)});
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
    row("q on the critical state circle = M a (kPa)", std::sqrt(1.5) * rcs, "50");
    std::cout << "  Phi/a^2 (TPZYCModifiedCamClayRHW::YieldFunction) at points of the surface and inside it:\n";
    const REAL pts[6][2] = {{0., 0.}, {25., M_PI / 6.}, {50., M_PI / 3.}, {80., 2.}, {100., 4.}, {50., -1.}};
    for (auto &pt : pts) {
        const REAL pp = pt[0];
        const bool center = pt[1] < 0.;
        const REAL xi = -sq3 * pp;
        const REAL rho = center ? 0. : SurfaceRadius(yc, a, xi);
        TPZManVector<REAL, 3> cyl = {xi, rho, center ? 0. : pt[1]}, s(3);
        TPZHWTools::FromHWCylToPrincipal(cyl, s);
        TPZManVector<STATE, 1> phi(1, 0.);
        yc.YieldFunction(s, pc, phi);
        std::ostringstream name;
        name << (center ? "centre of the ellipse" : "surface") << ": p' = " << pp << ", q = " << std::sqrt(1.5) * rho;
        if (!center) name << ", beta = " << pt[1];
        row(name.str(), phi[0] / (a * a), center ? "-1" : "0");
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
    TPZManVector<REAL, 4> X(4, 0.);
    r.fConverged = yc.ProjectHW(trial, X, r.fIter);
    REAL H, pc;
    yc.Hardening(fPcn, X[2], fV0, r.fA, H, pc);
    yc.Hardening(fPcn, 0., fV0, r.fAn, H, pc);
    r.fP = -X[0] / std::sqrt(3.);
    r.fQ = X[1] * std::sqrt(1.5);
    r.fDal = X[2];
    r.fDg = X[3];
    return r;
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

inline void YieldSurfaceProjection::WriteMeridianFiles(const TMeridianCase &c, const TMeridianResult &r) const {
    // ellipses p' in [0, 2a] (ellipse of fig_surface.py) and the energy-norm contour through the projected state,
    // (p' - p'_tr)^2/K + (q - q_tr)^2/(3G) = d^2
    const REAL d2 = (c.fPTrial - r.fP) * (c.fPTrial - r.fP) / fK + (c.fQTrial - r.fQ) * (c.fQTrial - r.fQ) / (3. * fG);
    std::vector<std::vector<REAL>> rows;
    for (int k = 0; k < fNCurve; ++k) {
        const REAL t = ysp::Linspace(0., M_PI, fNCurve, k);
        const REAL th = ysp::Linspace(0., 2. * M_PI, fNCurve, k);
        std::vector<REAL> rw;
        for (REAL a : {r.fAn, r.fA}) {
            const REAL p = -a + a * std::cos(t);
            rw.push_back(ysp::NoNegativeZero(-p));
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
    std::vector<TMeridianCase> cases(2);
    // reference values: fig_surface.py (cpp_linear) and camclay_hw.apply_strain with the same trial stress
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

    std::cout << "\n=== Fig. 2: closest-point projection in the meridian plane, linear elasticity K = " << fK
              << " kPa, G = " << fG << " kPa, M = " << fM2 << ", p_c,n = " << fPcn
              << " kPa, v0/(lambda - kappa) = " << std::setprecision(3) << fV0 / (fLambda - fKappa) << " ===\n";
    // critical state line q = M p' (two points, as in fig_meridian)
    mcc::WriteCSV("fig2_critical_state_line.csv", {"p", "q"}, {{0., 0.}, {240., fM2 * 240.}});

    for (auto &c : cases) {
        const TMeridianResult r = ProjectMeridian(c);
        WriteMeridianFiles(c, r);
        std::cout << std::defaultfloat << std::setprecision(6) << "\n" << c.fTitle << ": trial p' = " << c.fPTrial << " kPa, q = " << c.fQTrial << " kPa\n";
        std::cout << "  local problem (18), TPZYCModifiedCamClayRHW::ProjectHW" << (r.fConverged ? "" : " NOT CONVERGED")
                  << "\n";
        std::cout << "  " << std::left << std::setw(22) << "quantity" << std::right << std::setw(24) << "this code"
                  << std::setw(24) << "Python" << std::setw(11) << "rel.diff" << std::setw(10) << "article" << "\n";
        auto row = [](const std::string &name, REAL val, REAL ref, const std::string &art) {
            std::cout << "  " << std::left << std::setw(22) << name << std::right << std::setprecision(15)
                      << std::setw(24) << val << std::setw(24) << ref << std::setprecision(1) << std::scientific
                      << std::setw(11) << ysp::RelDiff(val, ref) << std::defaultfloat << std::setw(10) << art << "\n";
        };
        std::ostringstream art;
        art << std::fixed << std::setprecision(1) << c.fAArticle;
        row("p' (kPa)", r.fP, c.fPRef, "");
        row("q (kPa)", r.fQ, c.fQRef, "");
        row("a (kPa)", r.fA, c.fARef, art.str());
        row("a_n (kPa)", r.fAn, 100., "100");
        row("Delta alpha", r.fDal, c.fDalRef, c.fDalRef > 0 ? "> 0" : "< 0");
        row("Delta gamma", r.fDg, c.fDgRef, "");
        std::cout << "  " << std::left << std::setw(22) << "Newton iterations" << std::right << std::setw(24) << r.fIter
                  << std::setw(24) << c.fItRef << "\n";
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

        // complete stress update: triaxial trial stress (sigma_zz axial) and 12 Lode angles in a rotated frame
        TPZTensor<REAL> sigtr;
        sigtr.XX() = sigtr.YY() = -(c.fPTrial - c.fQTrial / 3.);
        sigtr.ZZ() = -(c.fPTrial + 2. * c.fQTrial / 3.);
        const TStressUpdateResult su = StressUpdate(sigtr);
        std::cout << std::setprecision(15) << "  stress update TPZPlasticStepModifiedCamClay::ApplyStrainComputeSigma "
                  << "(from sigma_n = -" << std::setprecision(4) << fPStart << " I kPa, eps_n = 0):\n";
        std::cout << std::setprecision(15) << "    trial p' = " << mcc::MeanEffectiveStress(su.fSigmaTrial)
                  << ", q = " << mcc::DeviatoricStress(su.fSigmaTrial) << (su.fFailed ? "  PROJECTION FAILED" : "")
                  << "\n";
        row("  p' (kPa)", mcc::MeanEffectiveStress(su.fSigma), c.fPRef, "");
        row("  q (kPa)", mcc::DeviatoricStress(su.fSigma), c.fQRef, "");
        row("  p_c (kPa)", su.fPc, c.fPcRef, "");
        std::cout << "    type " << su.fType << " (" << (su.fType == 1 ? "subcritical" : su.fType == 2 ? "supercritical" : "elastic")
                  << "), local iterations " << su.fIter << "\n";
        const int idx[6][2] = {{_XX_, _XX_}, {_XX_, _YY_}, {_XX_, _ZZ_}, {_ZZ_, _XX_}, {_ZZ_, _ZZ_}, {_XY_, _XY_}};
        const char *names[6] = {"D(xx,xx)", "D(xx,yy)", "D(xx,zz)", "D(zz,xx)", "D(zz,zz)", "D(xy,xy)"};
        std::cout << "    consistent tangent (kPa), non-symmetric:\n";
        for (int k = 0; k < 6; ++k) row(std::string("    ") + names[k], su.fDep(idx[k][0], idx[k][1]), c.fDRef[k], "");

        // invariance with respect to the Lode angle and to the orientation of the principal axes
        const TPZFNMatrix<9, REAL> R = ysp::Rotation(0.3, 0.7, 1.1);
        REAL dp = 0., dq = 0., dn = 0.;
        int ntests = 0;
        for (int k = 0; k < 12; ++k) {
            TPZManVector<REAL, 3> cyl = {-std::sqrt(3.) * c.fPTrial, std::sqrt(2. / 3.) * c.fQTrial, k * M_PI / 6.};
            TPZManVector<REAL, 3> s(3);
            TPZHWTools::FromHWCylToPrincipal(cyl, s);
            const TPZTensor<REAL> st = ysp::TensorFromPrincipal(s, R);
            const TStressUpdateResult u = StressUpdate(st);
            const TPZTensor<REAL> ntr = ysp::DeviatoricDirection(u.fSigmaTrial), n = ysp::DeviatoricDirection(u.fSigma);
            dp = std::max(dp, std::fabs(mcc::MeanEffectiveStress(u.fSigma) - c.fPRef));
            dq = std::max(dq, std::fabs(mcc::DeviatoricStress(u.fSigma) - c.fQRef));
            dn = std::max(dn, (n - ntr).Norm());
            if (u.fFailed) dp = 1.e300;
            ntests++;
        }
        std::cout << std::scientific << std::setprecision(1) << "    " << ntests
                  << " Lode angles beta = k pi/6 in a rotated frame: max |p' - p'_ref| = " << dp
                  << " kPa, max |q - q_ref| = " << dq << " kPa, max |n - n_tr| = " << dn << std::defaultfloat << "\n";
        std::cout << "  files: fig2" << c.fName << "_curves.csv, fig2" << c.fName << "_points.csv\n";
    }
    std::cout << "  file: fig2_critical_state_line.csv\n";
}

inline void YieldSurfaceProjection::Run() {
    RunSurface();
    RunProjection();
}
