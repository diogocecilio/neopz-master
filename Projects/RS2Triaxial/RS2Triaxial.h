/**
 * @file RS2Triaxial.h
 * @brief Sect. 6.1 of the article: drained triaxial tests of the RS2 manual (Rocscience, Sect. 8.7) at a
 * material point, Fig. 4 and Table 2.
 *
 * Mirrors the function rs2() of gen_data.py (camclay_hw.triaxial_point and camclay_hw.triaxial_closed).
 */
#pragma once

#include "MCCPaperTools.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "pzstepsolver.h"

#include <array>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <vector>

/**
 * @ingroup mccpaper
 * @brief Drained triaxial compression tests of the RS2 manual with the Modified Cam-Clay return mapping
 * in rotated Haigh-Westergaard space (TPZPlasticStepModifiedCamClay).
 *
 * Material (Table 1): \f$M=1.2\f$, \f$\lambda=0.077\f$, \f$\kappa=0.0066\f$, \f$v_0=1.70\f$, porous elasticity
 * with constant Poisson's ratio \f$\nu=0.3\f$ or constant shear modulus \f$G=20\f$ MPa, \f$p_t=0\f$, \f$\omega=1\f$.
 * The cell pressure is kept equal to the initial isotropic effective stress \f$p'_0\f$ and the axial strain is
 * increased to \f$\varepsilon_a=20\%\f$. Four cases (figures of the RS2 manual):
 *  - Fig8.5: normally consolidated (NC), \f$p'_0=p'_{c0}=200\f$ kPa, constant \f$\nu\f$;
 *  - Fig8.6: NC, \f$p'_0=p'_{c0}=200\f$ kPa, constant \f$G\f$;
 *  - Fig8.7: OCR = 2, \f$p'_0=100\f$, \f$p'_{c0}=200\f$ kPa, constant \f$\nu\f$;
 *  - Fig8.8: OCR = 5, \f$p'_0=100\f$, \f$p'_{c0}=500\f$ kPa, constant \f$\nu\f$ (yield in the supercritical region).
 *
 * For each case the example computes
 *  -# the closed-form solution of Appendix B.1 (mcc::TriaxialDrainedClosed with npts = 600, as gen_data.py:
 *     41 points on the elastic branch and two sets of 600 stress ratios on the plastic branch, uniform and
 *     clustered near \f$M\f$);
 *  -# the material point solution (mcc::TriaxialDrained, a transcription of TriaxialPointCC: \f$\varepsilon_{zz}\f$
 *     prescribed in equal increments, lateral strains from a Newton iteration on
 *     \f$\sigma'_{xx}=\sigma'_{yy}=-p'_0\f$ with the xx-yy block of the consistent tangent) with 100, 200, 400, 800
 *     and 1600 increments, and the error \f$q_{num}-q_{exact}\f$ at \f$\varepsilon_a=20\%\f$ (Table 2);
 *  -# the largest difference along the path between the 400-increment solution and the closed form,
 *     \f$\max_k|q_k-q_{closed}(\varepsilon_{a,k})|\f$, with \f$q_{closed}\f$ linearly interpolated at the numerical
 *     axial strains (numpy.interp, as in gen_data.py);
 *  -# for OCR = 5, the peak of \f$q\f$, the value at 20% and the maximum compression \f$\varepsilon_v\f$;
 *  -# a finite element check: the 400-increment test repeated with one Hex20-Hex8 u-p element on the unit cube
 *     (2 x 2 x 2 Gauss points; u_x = 0 on x = 0, u_y = 0 on y = 0, u_z = 0 on z = 0, cell pressure p'0 on the
 *     faces x = 1 and y = 1, vertical displacement of the top z = 1 controlled, pore pressure prescribed as zero
 *     at the eight vertices), assembled and solved with the native NeoPZ classes (geometric mesh, atomic and
 *     multiphysics computational meshes, TPZPoroElastoPlasticUPAnalysis with TPZSkylineNSymStructMatrix and
 *     TPZStepSolver). The strain field is homogeneous, so the element must reproduce the material point; the
 *     same element, with the boundary conditions of each test, is used for the FLAC3D tests (FLAC3DTriaxial).
 *
 * Two definitions of \f$q_{exact}\f$ at \f$\varepsilon_a=20\%\f$ are reported:
 *  - interpolated: linear interpolation of the 600-point closed-form table, as gen_data.py and the reference
 *    numbers of the Python transcription;
 *  - exact: the stress ratio \f$\eta\f$ with \f$\varepsilon_a(\eta)=0.2\f$ is found by bisection on the
 *    parametric closed form (ClosedFormAtRatio). The differences are below \f$10^{-3}\f$ kPa; with this value the
 *    errors agree with Table 2 of the article in 19 of the 20 entries to the printed digits (the exception,
 *    NC with constant G and 800 increments, is -0.117497 kPa, at the rounding boundary of the article's -0.118)
 *    and the four values of \f$q_{exact}\f$ agree, whereas the interpolated value differs from the article in 8
 *    entries and in \f$q_{exact}\f$ of NC with constant G (387.879).
 *
 * Output files (working directory), with \<case\> = nc_nu, nc_g, ocr2, ocr5:
 *  - rs2_\<case\>_closed.csv: closed form up to 20% (eps_a, p', q, eps_v, eps_q, sigma_a), dashed lines of Fig. 4;
 *  - rs2_\<case\>_n400.csv: material point with 400 increments (same columns, plus the closed form interpolated at
 *    eps_a and the difference), solid lines of Fig. 4;
 *  - rs2_table2.csv: Table 2 (both definitions of q_exact; the last row, increments = 0, holds q_exact);
 *  - rs2_\<case\>_fe.csv, rs2_\<case\>_fe.scal_vec.0.vtk, rs2_\<case\>_fe_gauss.vtk: finite element check;
 *  - rs2_fe_check.csv: summary of the finite element check (q at 20%, largest differences from the material
 *    point, global iterations, wall time);
 *  - rs2_mesh_*.csv: geometric mesh of the element (mcc::WriteMeshCSV), for the figure of the model.
 *
 * Stresses in kPa, compression positive in the printed and written invariants (p', q, eps_v, eps_a).
 */
class RS2Triaxial {
public:
    /** @brief Number of entries (increment counts) of the convergence study (Table 2) */
    static constexpr int NIncrements = 5;

    /**
     * @brief Material and boundary ids of the single element finite element check: faces x = 0, x = 1, y = 0,
     * y = 1, z = 0, z = 1 (displacement and cell pressure) and the drained faces (pore pressure, coincident with
     * the six faces); the same ids as FLAC3DTriaxial
     */
    enum { EMatId = 1, EX0 = -21, EX1 = -22, EY0 = -23, EY1 = -24, EZ0 = -25, EZ1 = -26, EPDrained = -30 };

    /**
     * @brief Reference values of a test (Python transcription and article); "400" refers to the default number
     * of increments of the curves of Fig. 4 (fNFigure)
     */
    struct TReference {
        REAL fQExactInterp = 0.;                ///< Python: q_exact(0.2) by numpy.interp on the closed form
        std::array<REAL, NIncrements> fErr{};   ///< Python: q_num - q_exact(interp) at eps_a = 0.2
        REAL fQEnd400 = 0.;                     ///< Python: q at eps_a = 0.2 with 400 increments
        REAL fEvEnd400 = 0.;                    ///< Python: eps_v at eps_a = 0.2 with 400 increments
        REAL fMaxDiff400 = 0.;                  ///< Python: max |q - q_closed| along the path, 400 increments
        REAL fEaMaxDiff400 = 0.;                ///< Python: eps_a of the largest difference
        std::array<REAL, NIncrements> fErrArticle{}; ///< article, Table 2 (3 decimals)
        REAL fQExactArticle = 0.;               ///< article, Table 2
        /// article: max |q - q_closed| with 400 increments (text of Sect. 6.1 or Table 3; NaN if not reported)
        REAL fMaxDiffArticle = std::numeric_limits<REAL>::quiet_NaN();
        REAL fLocalItsMean = 0.;                ///< Python: mean local Newton iterations per projection (400 increments)
        int fLocalItsMax = 0;                   ///< Python: maximum local Newton iterations (400 increments)
        /// article: mean local Newton iterations with 400 increments (Table 3; NaN if not reported)
        REAL fLocalItsArticle = std::numeric_limits<REAL>::quiet_NaN();
    };

    /** @brief Definition of a test */
    struct TCase {
        std::string fName;  ///< figure of the RS2 manual (Fig8.5 ... Fig8.8)
        std::string fTag;   ///< tag of the output files
        std::string fLabel; ///< label of Fig. 4 / Table 2
        REAL fP0 = 200.;    ///< initial isotropic effective stress p'0 (kPa), also the cell pressure
        REAL fPc0 = 200.;   ///< initial preconsolidation pressure p'c0 (kPa)
        bool fConstantG = false; ///< constant shear modulus G (true) or constant Poisson's ratio (false)
        TReference fRef;    ///< reference values
    };

    /** @brief Results of a test */
    struct TResult {
        std::vector<std::array<REAL, 6>> fClosed;                 ///< closed form (eps_a, p', q, eps_v, eps_q, sigma_a)
        std::map<int, std::vector<std::array<REAL, 6>>> fPoint;   ///< material point solutions by number of increments
        REAL fQExactInterp = 0.;                                  ///< q_exact(0.2), interpolated
        REAL fQExactRoot = 0.;                                    ///< q_exact(0.2), exact (bisection on eta)
        std::array<REAL, NIncrements> fErrInterp{};               ///< q_num - q_exact(interp)
        std::array<REAL, NIncrements> fErrRoot{};                 ///< q_num - q_exact(exact)
        REAL fMaxDiff400 = 0., fEaMaxDiff400 = 0.;                ///< largest |q - q_closed| along the path (400)
        mcc::TLocalStats fStats400;                               ///< local Newton iterations (400 increments)
        std::vector<std::array<REAL, 6>> fFE;                     ///< FE check (eps_a, p', q, eps_v, eps_q, evaluations)
        REAL fFEMaxDiff = 0.;                                     ///< max |q_FE - q_point| along the path (400)
        REAL fFEMaxDiffP = 0.;                                    ///< max |p'_FE - p'_point| along the path
        REAL fFEMaxDiffEv = 0.;                                   ///< max |eps_v,FE - eps_v,point| along the path
        REAL fFEMeanEvaluations = 0.;                             ///< mean global residual evaluations per increment
        int fFEMaxEvaluations = 0;                                ///< largest number of evaluations of an increment
        int64_t fFEGlobalIterations = 0;                          ///< all the global iterations (NGlobalIterations)
        int64_t fFEBisections = 0;                                ///< bisected increments (NBisections)
        REAL fFEWallTime = 0.;                                    ///< wall time of the incremental solution (s)
        REAL fFESpreadQ = 0.;                                     ///< largest spread of q over the 8 points (kPa)
        int64_t fFEEquations = 0;                                 ///< equations of the multiphysics mesh
    };

    /** @name Material parameters (Table 1) */
    /** @{ */
    REAL fM = 1.2, fLambda = 0.077, fKappa = 0.0066, fV0 = 1.70, fG = 2.e4, fNu = 0.3;
    /** @} */
    /** @brief Final axial strain of the tests */
    REAL fEaMax = 0.2;
    /** @brief Number of points of each branch of the closed form (npts of triaxial_closed) */
    int fNClosed = 600;
    /** @brief Numbers of increments of the convergence study */
    std::array<int, NIncrements> fIncrements = {100, 200, 400, 800, 1600};
    /** @brief Number of increments of the curves of Fig. 4 and of the finite element check */
    int fNFigure = 400;
    /** @brief Runs the finite element check */
    bool fFECheck = true;

    /** @brief The four tests with their reference values */
    std::vector<TCase> Cases() const;

    /** @brief Plastic model of a test: MCC with porous elasticity and the shear law of the case */
    mcc::TPlastic CreateMaterial(const TCase &tcase) const;

    /** @name Closed form (Appendix B.1) */
    /** @{ */
    /**
     * @brief Yield point of the drained path \f$q=3(p'-p'_0)\f$:
     * \f$(9+M^2)p'^2-(18p'_0+M^2p'_{c0})p'+9p'^2_0=0\f$
     * @param[out] py mean effective stress at yield
     * @return stress ratio at yield \f$\eta_y=q_y/p'_y\f$
     */
    REAL YieldRatio(const TCase &tcase, REAL &py) const;

    /**
     * @brief Closed-form state on the plastic branch for the stress ratio eta, eqs. (B.1)-(B.6):
     * \f$p'=3p'_0/(3-\eta)\f$, \f$q=\eta p'\f$, \f$p'_c=p'(1+\eta^2/M^2)\f$,
     * \f$\varepsilon_v=\varepsilon_v^e+\frac{\lambda-\kappa}{v_0}\ln\frac{p'_c}{p'_{c0}}\f$,
     * \f$\varepsilon_q=\varepsilon_q^e+\frac{\lambda-\kappa}{v_0}(F(\eta)-F(\eta_y))\f$,
     * \f$\varepsilon_a=\varepsilon_q+\varepsilon_v/3\f$ (same expressions of mcc::TriaxialDrainedClosed)
     * @return (eps_a, p', q, eps_v, eps_q)
     */
    std::array<REAL, 5> ClosedFormAtRatio(const TCase &tcase, REAL eta) const;

    /**
     * @brief Deviatoric stress of the closed form at the axial strain ea, without interpolation: bisection on
     * \f$\eta\f$ between \f$\eta_y\f$ and \f$M\f$ (increasing in the subcritical region, decreasing in the
     * supercritical region) for \f$\varepsilon_a(\eta)=\varepsilon_a\f$
     * @param ea axial strain beyond the yield point (NaN is returned in the elastic branch)
     */
    REAL ExactDeviatoricStress(const TCase &tcase, REAL ea) const;
    /** @} */

    /** @name Finite element check: one Hex20-Hex8 element */
    /** @{ */
    /**
     * @brief Geometric mesh: the unit cube (one trilinear hexahedron) with a boundary quadrilateral for the
     * displacement condition and a coincident one for the pore pressure on each face
     */
    TPZGeoMesh *CreateGeoMesh();

    /** @brief Displacement (Hex20), pore pressure (Hex8) and multiphysics meshes, materials and initial state */
    TPZMultiphysicsCompMesh *CreateCompMesh(TPZGeoMesh *gmesh, const TCase &tcase, mcc::TPoroMaterial *&mat);

    /**
     * @brief Solves the test with the element: analysis, structural matrix (TPZSkylineNSymStructMatrix), direct
     * solver (TPZStepSolver, ELU), nsteps increments of the top displacement and post-processing (CSV, VTK)
     * @param tcase test
     * @param nsteps number of increments; res.fPoint must hold the material point solution with nsteps increments
     * @param res results of the test: fFE, fFEMaxDiff and fFEMeanEvaluations are filled
     */
    void RunFiniteElement(const TCase &tcase, int nsteps, TResult &res);
    /** @} */

    /** @brief Runs a test: closed form, material point solutions, finite element check and CSV files */
    TResult Run(const TCase &tcase);

    /** @brief Runs the four tests and prints Table 2 and the values of Fig. 4 next to the reference values */
    void RunAll();

private:
    /** @brief Function F of eq. (B.5) */
    REAL F(REAL x) const {
        const REAL M = fM;
        return (1. / M) * std::log(std::fabs((M + x) / (M - x))) - (2. / M) * std::atan(x / M) -
               std::log(std::fabs(M - x)) / (3. - M) - std::log(M + x) / (3. + M) + 6. * std::log(3. - x) / (9. - M * M);
    }
    /** @brief Writes a table of fixed-size rows as CSV */
    template <size_t N>
    static void WriteTable(const std::string &file, const std::vector<std::string> &header,
                           const std::vector<std::array<REAL, N>> &tab) {
        std::vector<std::vector<REAL>> rows;
        for (auto &r : tab) rows.emplace_back(r.begin(), r.end());
        mcc::WriteCSV(file, header, rows);
    }
};

inline std::vector<RS2Triaxial::TCase> RS2Triaxial::Cases() const {
    const REAL nan = std::numeric_limits<REAL>::quiet_NaN();
    std::vector<TCase> cases(4);
    // Reference values, in the order of TReference:
    //   Python (gen_data.py rs2, RS2 lines of the reference numbers): q_exact (interpolated), errors for
    //   100...1600 increments, q and eps_v at 20% (400 increments), max |q - q_closed| and its eps_a;
    //   article: Table 2 errors and q_exact, max |q - q_closed| (Sect. 6.1; Table 3 for constant G);
    //   Python: mean and maximum local iterations (400 increments); article: mean local iterations (Table 3).
    cases[0] = {"Fig8.5", "nc_nu", "NC, constant nu", 200., 200., false,
                {388.345014514804,
                 {-0.9271152607461772, -0.4653471852926714, -0.23319763506560776, -0.11660994508389422,
                  -0.05815190967706485},
                 388.1118168797384, 0.05055284332005455, 1.38335978861, 0.0085,
                 {-0.928, -0.466, -0.234, -0.117, -0.059}, 388.345, 1.4,
                 4.070781893004115, 5, nan}};
    cases[1] = {"Fig8.6", "nc_g", "NC, constant G", 200., 200., true,
                {387.87931751567714,
                 {-0.9322888920397645, -0.4676679620955042, -0.234278829244829, -0.11715511929116929,
                  -0.05845798283939985},
                 387.6450386864323, 0.05050173562793553, 1.25602047411, 0.011,
                 {-0.933, -0.468, -0.235, -0.118, -0.059}, 387.880, 1.26,
                 4.051070840197694, 5, 4.05}};
    cases[2] = {"Fig8.7", "ocr2", "OCR = 2", 100., 200., false,
                {196.86572443881326,
                 {-0.20723535233290136, -0.10393816753074248, -0.05187258970593689, -0.02585923859714967,
                  -0.01285719431635357},
                 196.81385184910732, 0.022449480426678547, 2.7352857887, 0.003,
                 {-0.207, -0.104, -0.052, -0.026, -0.013}, 196.866, nan,
                 4., 4, nan}};
    cases[3] = {"Fig8.8", "ocr5", "OCR = 5", 100., 500., false,
                {202.87765357015488,
                 {0.19603904867420852, 0.09867007181884446, 0.04912937943888096, 0.024558549624657644,
                  0.012244299472143894},
                 202.92678294959376, -0.014181972004958021, 8.46654345612, 0.0065,
                 {0.196, 0.099, 0.049, 0.025, 0.012}, 202.878, 8.5,
                 4., 4, nan}};
    return cases;
}

inline mcc::TPlastic RS2Triaxial::CreateMaterial(const TCase &tcase) const {
    mcc::TPlastic model;
    model.SetModifiedCamClay(fM, fLambda, fKappa);
    model.SetPorousElasticity();
    if (tcase.fConstantG)
        model.SetConstantShearModulus(fG);
    else
        model.SetPoissonRatio(fNu);
    model.SetDefaultSpecificVolume(fV0);
    return model;
}

inline REAL RS2Triaxial::YieldRatio(const TCase &tcase, REAL &py) const {
    const REAL p0 = tcase.fP0, pc0 = tcase.fPc0, M = fM;
    const REAL A = 9. + M * M, B = -(18. * p0 + M * M * pc0), C = 9. * p0 * p0;
    const REAL disc = std::sqrt(B * B - 4. * A * C);
    py = 1e300;
    for (REAL root : {(-B + disc) / (2. * A), (-B - disc) / (2. * A)})
        if (root >= p0 - 1e-9) py = std::min(py, root);
    return 3. * (py - p0) / py;
}

inline std::array<REAL, 5> RS2Triaxial::ClosedFormAtRatio(const TCase &tcase, REAL eta) const {
    REAL py;
    const REAL etay = YieldRatio(tcase, py);
    const REAL p0 = tcase.fP0;
    const REAL pp = 3. * p0 / (3. - eta), q = eta * pp, pc = pp * (1. + eta * eta / (fM * fM));
    const REAL r = 3. * (1. - 2. * fNu) / (2. * (1. + fNu));
    REAL ev = fKappa / fV0 * std::log(pp / p0);
    REAL eq = tcase.fConstantG ? q / (3. * fG) : fKappa / (r * fV0) * std::log(pp / p0);
    ev += (fLambda - fKappa) / fV0 * std::log(pc / tcase.fPc0);
    eq += (fLambda - fKappa) / fV0 * (F(eta) - F(etay));
    return {eq + ev / 3., pp, q, ev, eq};
}

inline REAL RS2Triaxial::ExactDeviatoricStress(const TCase &tcase, REAL ea) const {
    REAL py;
    const REAL etay = YieldRatio(tcase, py);
    // eps_a grows monotonically from the yield point (eta = eta_y) to infinity (eta -> M)
    REAL a = etay, b = fM - (fM - etay) * 1e-15;
    if (ClosedFormAtRatio(tcase, a)[0] > ea) {
        std::cerr << "RS2Triaxial::ExactDeviatoricStress: eps_a = " << ea << " is in the elastic branch\n";
        return std::numeric_limits<REAL>::quiet_NaN();
    }
    for (int it = 0; it < 200 && a != b; ++it) {
        const REAL c = 0.5 * (a + b);
        if (c == a || c == b) break;
        if (ClosedFormAtRatio(tcase, c)[0] < ea)
            a = c;
        else
            b = c;
    }
    return ClosedFormAtRatio(tcase, 0.5 * (a + b))[2];
}

inline TPZGeoMesh *RS2Triaxial::CreateGeoMesh() {
    // one trilinear hexahedron on [0,1]^3 (the serendipity displacement space needs no geometric mid-edge nodes);
    // each face receives a quadrilateral for the displacement condition and a coincident one for the pore pressure
    return mcc::CreateUnitCubeMesh(EMatId, [](const std::array<std::array<REAL, 3>, 4> &X) {
        std::vector<int> ids;
        const int planes[3][2] = {{EX0, EX1}, {EY0, EY1}, {EZ0, EZ1}};
        for (int axis = 0; axis < 3; ++axis)
            for (int side = 0; side < 2; ++side)
                if (mcc::FaceOnPlane(X, axis, REAL(side))) ids.push_back(planes[axis][side]);
        ids.push_back(EPDrained);
        return ids;
    });
}

inline TPZMultiphysicsCompMesh *RS2Triaxial::CreateCompMesh(TPZGeoMesh *gmesh, const TCase &tcase,
                                                            mcc::TPoroMaterial *&mat) {
    const std::set<int> bcids = {EX0, EX1, EY0, EY1, EZ0, EZ1, EPDrained};
    TPZCompMesh *cmeshU = mcc::CreateDisplacementMesh(gmesh, 3, EMatId, bcids); // serendipity Hex20
    TPZCompMesh *cmeshP = mcc::CreatePressureMesh(gmesh, 3, EMatId, bcids);     // trilinear Hex8

    TPZMultiphysicsCompMesh *mphys = new TPZMultiphysicsCompMesh(gmesh);
    mphys->SetDimModel(3);
    mat = new mcc::TPoroMaterial(EMatId, TPZMatPoroElastoPlasticUPBase::EThreeDimensional);
    mat->SetPlasticModel(CreateMaterial(tcase));
    mat->SetBiot(1., 0.);       // drained: incompressible constituents, pore pressure prescribed as zero
    mat->SetPermeability(0.);
    mat->SetIntegrationOrder(3); // 2 x 2 x 2 Gauss points
    mphys->InsertMaterialObject(mat);

    using B = TPZMatPoroElastoPlasticUPBase;
    TPZFNMatrix<9, STATE> val1(3, 3, 0.);
    TPZManVector<STATE, 3> val2(3, 0.);
    auto directional = [&](int id, int comp) {
        val1.Zero();
        val1(comp, comp) = 1.;
        val2.Fill(0.);
        mphys->InsertMaterialObject(mat->CreateBC(mat, id, B::EDirichletUDirectional, val1, val2));
    };
    directional(EX0, 0); // symmetry planes
    directional(EY0, 1);
    directional(EZ0, 2);
    directional(EZ1, 2); // top: u_z controlled by the analysis
    // faces x = 1 and y = 1: cell pressure p'0 (total = effective stress, the pore pressure is zero)
    for (int k = 0; k < 2; ++k) {
        val1.Zero();
        val2.Fill(0.);
        val2[k] = -tcase.fP0;
        mphys->InsertMaterialObject(mat->CreateBC(mat, k == 0 ? EX1 : EY1, B::ENeumannU, val1, val2));
    }
    // drained: p = 0 on the whole boundary (the eight vertices)
    TPZManVector<STATE, 3> zero(1, 0.);
    val1.Zero();
    mphys->InsertMaterialObject(mat->CreateBC(mat, EPDrained, B::EDirichletP, val1, zero));
    mcc::BuildMultiphysics(mphys, cmeshU, cmeshP, EMatId, bcids);

    // initial state of the integration points: isotropic stress p'0, preconsolidation p'c0, specific volume v0
    const REAL v0 = fV0;
    mat->InitializeMemory(mphys, [&](const TPZVec<REAL> &, TPZElastoPlasticMem &mem) {
        mem.m_sigma = mcc::IsotropicTensor(-tcase.fP0);
        mem.m_elastoplastic_state.m_hardening = tcase.fPc0;
        mem.m_elastoplastic_state.fmatprop.Resize(1, v0);
        mem.m_elastoplastic_state.fmatprop[0] = v0;
        mem.m_elastoplastic_state.fpressure = 0.;
    });
    return mphys;
}

inline void RS2Triaxial::RunFiniteElement(const TCase &tcase, int nsteps, TResult &res) {
    TPZGeoMesh *gmesh = CreateGeoMesh();
    mcc::TPoroMaterial *mat = nullptr;
    TPZMultiphysicsCompMesh *mphys = CreateCompMesh(gmesh, tcase, mat);
    res.fFEEquations = mphys->NEquations();

    // analysis: non-symmetric skyline matrix and LU decomposition (no renumbering, see TPZPoroElastoPlasticUPAnalysis)
    mcc::TAnalysis analysis(mphys, mat);
    TPZSkylineNSymStructMatrix<STATE> skyl(mphys);
    skyl.SetNumThreads(0);
    analysis.SetStructuralMatrix(skyl);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELU);
    analysis.SetSolver(step);
    analysis.SetControlledDisplacement(EZ1, 2);
    analysis.SetPredictor(true);

    // increments of the top displacement (height 1 m: u_z = -eps_a), constant cell pressure (load factor 1)
    std::vector<mcc::TAnalysis::TLoadState> steps;
    for (int k = 1; k <= nsteps; ++k) steps.emplace_back(0., 1., -fEaMax * k / nsteps);

    // XX and YY are the lateral components, ZZ the axial one; the state is homogeneous, so the first integration
    // point represents the element (0. - x avoids writing -0 for the initial state); the spread of q over the
    // points measures the deviation from homogeneity
    res.fFE.clear();
    res.fFESpreadQ = 0.;
    res.fFEMaxEvaluations = 0;
    size_t nlog = 0;
    auto monitor = [&](int, const mcc::TAnalysis::TLoadState &s) {
        const auto gps = mcc::GaussPoints(mat, mphys);
        const auto &g = gps[0];
        const REAL q = mcc::DeviatoricStress(g.fSigma);
        for (auto &gi : gps) res.fFESpreadQ = std::max(res.fFESpreadQ, std::fabs(mcc::DeviatoricStress(gi.fSigma) - q));
        int nev = 0;
        const auto &log = analysis.StepLog();
        for (; nlog < log.size(); ++nlog) nev += int(log[nlog].fResiduals.size());
        res.fFEMaxEvaluations = std::max(res.fFEMaxEvaluations, nev);
        res.fFE.push_back({0. - s.fUc, mcc::MeanEffectiveStress(g.fSigma), q, 0. - g.fEps.I1(),
                           2. / 3. * std::fabs(g.fEps.ZZ() - g.fEps.XX()), REAL(nev)});
    };
    analysis.ResetCounters();
    const auto start = std::chrono::steady_clock::now();
    if (!analysis.Run(steps, monitor)) std::cerr << "RS2Triaxial: the finite element solution failed\n";
    res.fFEWallTime = std::chrono::duration<REAL>(std::chrono::steady_clock::now() - start).count();
    res.fFEMeanEvaluations = mcc::MeanEvaluations(analysis.StepLog());
    res.fFEGlobalIterations = analysis.NGlobalIterations();
    res.fFEBisections = analysis.NBisections();

    // post-processing: integration points and nodal fields (VTK) and the history (CSV)
    const std::string base = "rs2_" + tcase.fTag + "_fe";
    mcc::WriteGaussPointsVTK(mat, mphys, base + "_gauss.vtk");
    mcc::WriteNodalVTK(analysis, 3, base + ".vtk", 0);

    // comparison with the material point solution with the same increments
    const auto &mp = res.fPoint.at(nsteps);
    res.fFEMaxDiff = res.fFEMaxDiffP = res.fFEMaxDiffEv = 0.;
    std::vector<std::vector<REAL>> rows;
    for (size_t k = 0; k < res.fFE.size() && k < mp.size(); ++k) {
        const auto &f = res.fFE[k];
        res.fFEMaxDiff = std::max(res.fFEMaxDiff, std::fabs(f[2] - mp[k][2]));
        res.fFEMaxDiffP = std::max(res.fFEMaxDiffP, std::fabs(f[1] - mp[k][1]));
        res.fFEMaxDiffEv = std::max(res.fFEMaxDiffEv, std::fabs(f[3] - mp[k][3]));
        rows.push_back({f[0], f[1], f[2], f[3], f[4], mp[k][1], mp[k][2], mp[k][3], f[2] - mp[k][2], f[5]});
    }
    mcc::WriteCSV(base + ".csv",
                  {"eps_a", "p_eff", "q", "eps_v", "eps_q", "p_point", "q_point", "eps_v_point", "q_minus_q_point",
                   "evaluations"},
                  rows);
    mcc::DeleteMeshes(mphys);
}

inline RS2Triaxial::TResult RS2Triaxial::Run(const TCase &tcase) {
    TResult res;
    const mcc::TPlastic model = CreateMaterial(tcase);
    const REAL G = tcase.fConstantG ? fG : 0.;

    // closed form (Appendix B.1) and the value of q at eps_a = 20%
    res.fClosed = mcc::TriaxialDrainedClosed(tcase.fP0, tcase.fPc0, fV0, fM, fLambda, fKappa, G, fNu, fNClosed);
    res.fQExactInterp = mcc::Interpolate(res.fClosed, fEaMax, 2);
    res.fQExactRoot = ExactDeviatoricStress(tcase, fEaMax);

    // material point solutions (TriaxialPointCC) and the errors at eps_a = 20% (Table 2)
    for (int i = 0; i < NIncrements; ++i) {
        const int n = fIncrements[i];
        mcc::TLocalStats stats;
        auto path = mcc::TriaxialDrained(model, tcase.fP0, tcase.fPc0, fV0, fEaMax, n, &stats);
        if (path.empty()) {
            std::cerr << "RS2Triaxial: projection failure, " << tcase.fName << ", " << n << " increments\n";
            res.fErrInterp[i] = res.fErrRoot[i] = std::numeric_limits<REAL>::quiet_NaN();
            continue;
        }
        res.fErrInterp[i] = path.back()[2] - res.fQExactInterp;
        res.fErrRoot[i] = path.back()[2] - res.fQExactRoot;
        if (n == fNFigure) res.fStats400 = stats;
        res.fPoint[n] = std::move(path);
    }

    // largest difference along the path: closed form interpolated at the numerical axial strains (numpy.interp)
    std::vector<std::vector<REAL>> rows;
    const auto &num = res.fPoint.at(fNFigure);
    for (auto &r : num) {
        const REAL qc = mcc::Interpolate(res.fClosed, r[0], 2), evc = mcc::Interpolate(res.fClosed, r[0], 3);
        if (std::fabs(r[2] - qc) > res.fMaxDiff400) {
            res.fMaxDiff400 = std::fabs(r[2] - qc);
            res.fEaMaxDiff400 = r[0];
        }
        rows.push_back({r[0], r[1], r[2], r[3], r[4], r[5], qc, r[2] - qc, evc});
    }

    // curves of Fig. 4: material point (fNFigure increments) and closed form up to eps_a = 20%
    const std::string base = "rs2_" + tcase.fTag;
    mcc::WriteCSV(base + "_n" + std::to_string(fNFigure) + ".csv",
                  {"eps_a", "p_eff", "q", "eps_v", "eps_q", "sigma_a", "q_closed", "q_minus_q_closed", "eps_v_closed"},
                  rows);
    std::vector<std::array<REAL, 6>> closed;
    for (auto &r : res.fClosed)
        if (r[0] <= fEaMax + 1e-9) closed.push_back(r);
    WriteTable(base + "_closed.csv", {"eps_a", "p_eff", "q", "eps_v", "eps_q", "sigma_a"}, closed);

    if (fFECheck) RunFiniteElement(tcase, fNFigure, res);
    return res;
}

inline void RS2Triaxial::RunAll() {
    const std::vector<TCase> cases = Cases();
    std::vector<TResult> results;
    for (auto &c : cases) results.push_back(Run(c));

    std::cout << "RS2 drained triaxial tests at a material point (Sect. 6.1, Fig. 4, Table 2)\n"
              << "M = " << fM << ", lambda = " << fLambda << ", kappa = " << fKappa << ", v0 = " << fV0
              << ", nu = " << fNu << " or G = " << fG / 1000. << " MPa; cell pressure p'0, eps_a up to "
              << 100. * fEaMax << "%\n\n";

    auto header = [&](const char *first, int width) {
        std::cout << std::setw(12) << std::left << first << std::right;
        for (auto &c : cases) std::cout << std::setw(width) << c.fLabel;
        std::cout << "\n";
    };
    auto cell = [](REAL a, REAL b, int prec) {
        std::ostringstream s;
        s << std::fixed << std::setprecision(prec) << a << " | " << b;
        return s.str();
    };
    // value of the article as printed there, or "-" if the article does not report it
    auto article = [](REAL a) {
        std::ostringstream s;
        if (std::isnan(a))
            s << "-";
        else
            s << a;
        return s.str();
    };

    std::cout << "Table 2 (a): q_num - q_exact at eps_a = 20% (kPa), q_exact interpolated on the closed form "
                 "(npts = 600)\n             [this code | Python gen_data.py]\n";
    header("Increments", 28);
    for (int i = 0; i < NIncrements; ++i) {
        std::cout << std::setw(12) << std::left << fIncrements[i] << std::right;
        for (size_t j = 0; j < cases.size(); ++j)
            std::cout << std::setw(28) << cell(results[j].fErrInterp[i], cases[j].fRef.fErr[i], 8);
        std::cout << "\n";
    }
    std::cout << std::setw(12) << std::left << "q_exact" << std::right;
    for (size_t j = 0; j < cases.size(); ++j)
        std::cout << std::setw(28) << cell(results[j].fQExactInterp, cases[j].fRef.fQExactInterp, 7);
    std::cout << "\n\nTable 2 (b): the same with q_exact solved exactly on the closed form (bisection on eta)\n"
                 "             [this code | article, Table 2]\n";
    header("Increments", 22);
    for (int i = 0; i < NIncrements; ++i) {
        std::cout << std::setw(12) << std::left << fIncrements[i] << std::right;
        for (size_t j = 0; j < cases.size(); ++j)
            std::cout << std::setw(22) << cell(results[j].fErrRoot[i], cases[j].fRef.fErrArticle[i], 3);
        std::cout << "\n";
    }
    std::cout << std::setw(12) << std::left << "q_exact" << std::right;
    for (size_t j = 0; j < cases.size(); ++j)
        std::cout << std::setw(22) << cell(results[j].fQExactRoot, cases[j].fRef.fQExactArticle, 3);
    std::cout << "\n";
    // entries that do not round to the article's value: print them with more digits
    for (size_t j = 0; j < cases.size(); ++j)
        for (int i = 0; i < NIncrements; ++i)
            if (std::fabs(results[j].fErrRoot[i] - cases[j].fRef.fErrArticle[i]) > 5.e-4)
                std::cout << "note: " << cases[j].fLabel << ", " << fIncrements[i] << " increments: " << std::fixed
                          << std::setprecision(6) << results[j].fErrRoot[i]
                          << " (rounding boundary of the article's value " << std::setprecision(3)
                          << cases[j].fRef.fErrArticle[i] << ")\n";
    REAL rel100 = 0.;
    for (size_t j = 0; j < cases.size(); ++j)
        rel100 = std::max(rel100, std::fabs(results[j].fErrRoot[0]) / results[j].fQExactRoot);
    std::cout << "largest error at eps_a = 20% with 100 increments: " << std::fixed << std::setprecision(2)
              << 100. * rel100 << "% (article: 0.24%)\n\n";

    std::cout << "Fig. 4, " << fNFigure
              << " increments [this code | Python | article: Sect. 6.1, Table 3 for constant G]\n";
    for (size_t j = 0; j < cases.size(); ++j) {
        const auto &r = results[j];
        const auto &ref = cases[j].fRef;
        const auto &e = r.fPoint.at(fNFigure).back();
        std::cout << "  " << std::setw(17) << std::left << cases[j].fLabel << std::right << "q(20%) = "
                  << cell(e[2], ref.fQEnd400, 6) << " kPa, eps_v(20%) = " << cell(e[3], ref.fEvEnd400, 8) << "\n"
                  << std::setw(19) << "" << "max|q - q_closed| along the path = " << cell(r.fMaxDiff400, ref.fMaxDiff400, 4)
                  << " | " << article(ref.fMaxDiffArticle) << " kPa at eps_a = " << std::setprecision(3)
                  << 100. * r.fEaMaxDiff400 << " | " << 100. * ref.fEaMaxDiff400 << " %\n";
    }

    // OCR = 5 (last case): peak, softening and dilation; reference values of the Python transcription
    // (400 increments, closed form with npts = 600) and of the text of Sect. 6.1 of the article
    {
        const REAL pyPeak = 292.95216553962547, pyPeakEa = 0.007, pyClosedPeak = 293.3863425,
                   pyClosedPeakEa = 0.00662003, pyQ20 = 202.92678294959376, pyEvMax = 0.00257737248965,
                   pyEvMaxEa = 0.007, pyEv20 = -0.014181972004958021;
        const REAL artPeak = 293.0, artPeakEa = 0.0070, artClosedPeak = 293.4, artClosedPeakEa = 0.0066,
                   artQ20 = 202.9, artEvMax = 0.0026;
        const auto &r = results.back();
        const auto &num = r.fPoint.at(fNFigure);
        size_t kq = 0, kv = 0, cq = 0;
        for (size_t k = 0; k < num.size(); ++k) {
            if (num[k][2] > num[kq][2]) kq = k;
            if (num[k][3] > num[kv][3]) kv = k;
        }
        for (size_t k = 0; k < r.fClosed.size(); ++k)
            if (r.fClosed[k][2] > r.fClosed[cq][2]) cq = k;
        auto cell3 = [](REAL a, REAL b, REAL c, int prec, int precart) {
            std::ostringstream s;
            s << std::fixed << std::setprecision(prec) << a << " | " << b << " | " << std::setprecision(precart) << c;
            return s.str();
        };
        std::cout << "\n" << cases.back().fLabel << ", " << fNFigure << " increments [this code | Python | article]\n"
                  << "  peak q            = " << cell3(num[kq][2], pyPeak, artPeak, 4, 1) << " kPa at eps_a = "
                  << cell3(100. * num[kq][0], 100. * pyPeakEa, 100. * artPeakEa, 3, 2) << " %\n"
                  << "  closed-form peak  = " << cell3(r.fClosed[cq][2], pyClosedPeak, artClosedPeak, 4, 1)
                  << " kPa at eps_a = " << cell3(100. * r.fClosed[cq][0], 100. * pyClosedPeakEa, 100. * artClosedPeakEa, 3, 2)
                  << " %\n"
                  << "  q(20%)            = " << cell3(num.back()[2], pyQ20, artQ20, 4, 1) << " kPa\n"
                  << "  max. compression  = " << cell3(100. * num[kv][3], 100. * pyEvMax, 100. * artEvMax, 4, 2)
                  << " % at eps_a = " << cell(100. * num[kv][0], 100. * pyEvMaxEa, 3) << " %; eps_v(20%) = "
                  << cell(100. * num.back()[3], 100. * pyEv20, 4) << " % (dilation)\n";
    }

    std::cout << "\nLocal Newton iterations per projection (" << fNFigure
              << " increments): mean [this code | Python | article, Table 3], max [this code | Python]\n";
    for (size_t j = 0; j < cases.size(); ++j)
        std::cout << "  " << std::setw(17) << std::left << cases[j].fLabel << std::right
                  << cell(results[j].fStats400.Mean(), cases[j].fRef.fLocalItsMean, 4) << " | "
                  << article(cases[j].fRef.fLocalItsArticle) << ", " << results[j].fStats400.fMax << " | "
                  << cases[j].fRef.fLocalItsMax << " (" << results[j].fStats400.fCalls << " projections)\n";

    if (fFECheck) {
        std::cout << "\nFinite element check: one Hex20-Hex8 element (2 x 2 x 2 points), " << fNFigure
                  << " increments [FE | material point]\n";
        std::vector<std::vector<REAL>> fe;
        for (size_t j = 0; j < cases.size(); ++j) {
            const auto &r = results[j];
            const auto &mp = r.fPoint.at(fNFigure).back();
            std::cout << "  " << std::setw(17) << std::left << cases[j].fLabel << std::right << "q(20%) = "
                      << cell(r.fFE.back()[2], mp[2], 9) << " kPa, max|q_FE - q_point| = " << std::scientific
                      << std::setprecision(2) << r.fFEMaxDiff << ", max|p'_FE - p'_point| = " << r.fFEMaxDiffP
                      << " kPa, max|eps_v,FE - eps_v,point| = " << r.fFEMaxDiffEv << ", spread of q over the points "
                      << r.fFESpreadQ << " kPa" << std::fixed << "\n" << std::setw(19) << ""
                      << "global evaluations per increment = " << std::setprecision(4) << r.fFEMeanEvaluations
                      << " (max " << r.fFEMaxEvaluations << "), total " << r.fFEGlobalIterations << ", bisections "
                      << r.fFEBisections << ", " << std::setprecision(3) << r.fFEWallTime << " s, " << r.fFEEquations
                      << " equations\n";
            fe.push_back({REAL(j), REAL(fNFigure), r.fFE.back()[2], mp[2], r.fFE.back()[1], mp[1], r.fFE.back()[3],
                          mp[3], r.fFEMaxDiff, r.fFEMaxDiffP, r.fFEMaxDiffEv, r.fFESpreadQ, r.fFEMeanEvaluations,
                          REAL(r.fFEMaxEvaluations), REAL(r.fFEGlobalIterations), REAL(r.fFEBisections),
                          r.fFEWallTime, REAL(r.fFEEquations)});
        }
        // summary of the check (case index in the order nc_nu, nc_g, ocr2, ocr5) and the mesh of the element
        mcc::WriteCSV("rs2_fe_check.csv",
                      {"case", "increments", "q_fe", "q_point", "p_fe", "p_point", "eps_v_fe", "eps_v_point",
                       "max_diff_q", "max_diff_p", "max_diff_eps_v", "spread_q", "mean_evaluations", "max_evaluations",
                       "global_iterations", "bisections", "wall_time_s", "equations"},
                      fe);
        TPZGeoMesh *gmesh = CreateGeoMesh();
        mcc::WriteMeshCSV(gmesh, "rs2_mesh");
        delete gmesh;
    }
    std::cout << "\nFiles: rs2_<case>_closed.csv, rs2_<case>_n" << fNFigure
              << ".csv (Fig. 4), rs2_table2.csv (Table 2), rs2_<case>_fe.csv, rs2_fe_check.csv, rs2_mesh_*.csv and VTK "
                 "files of the FE check\n";

    // Table 2 as CSV
    std::vector<std::vector<REAL>> rows;
    std::vector<std::string> head = {"increments"};
    for (auto &c : cases) {
        head.push_back("err_interp_" + c.fTag);
        head.push_back("err_exact_" + c.fTag);
    }
    for (int i = 0; i < NIncrements; ++i) {
        std::vector<REAL> row = {REAL(fIncrements[i])};
        for (auto &r : results) {
            row.push_back(r.fErrInterp[i]);
            row.push_back(r.fErrRoot[i]);
        }
        rows.push_back(row);
    }
    std::vector<REAL> qrow = {0.};
    for (auto &r : results) {
        qrow.push_back(r.fQExactInterp);
        qrow.push_back(r.fQExactRoot);
    }
    rows.push_back(qrow); // last row (increments = 0): q_exact
    mcc::WriteCSV("rs2_table2.csv", head, rows);
}
