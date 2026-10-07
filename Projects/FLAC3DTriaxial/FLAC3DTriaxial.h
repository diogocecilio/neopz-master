/**
 * @file FLAC3DTriaxial.h
 * @brief Sect. 6.3 of the article: drained and undrained triaxial tests of the FLAC3D verification problem
 * with a single Hex20-Hex8 u-p element (model in Fig. 4, results in Fig. 7 and Table 5), and the FLAC3D column of
 * Table 10 (Sect. 6.7): global iterations with five tangent operators.
 */
#pragma once

#include "MCCPaperTools.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "pzstepsolver.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <set>
#include <sstream>
#include <string>
#include <vector>

/**
 * @ingroup mccpaper
 * @brief Triaxial tests of the FLAC3D verification problem on a single three-dimensional u-p element.
 *
 * The FLAC3D zone, a cube of 1 m, is modelled by one Hex20-Hex8 element (quadratic serendipity displacement,
 * trilinear pore pressure) with 2 x 2 x 2 integration points. Boundary conditions: u_x = 0 on the face x = 0,
 * u_y = 0 on y = 0 and u_z = 0 on z = 0 (symmetry planes), total cell pressure p'0 = 5 kPa on the faces x = 1
 * and y = 1, and vertical displacement of the top z = 1 controlled (eps_a = -u_z). The initial state is the
 * isotropic effective stress p'0 with p'c0 = R p'0 and the specific volume on the normal compression line.
 *  - Drained tests: pore pressure prescribed as zero at the eight vertices (all faces drained), eps_a up to 50% in
 *    500 increments.
 *  - Undrained tests: no flow (k = 0, Dt = 0, no drained face), pore fluid with K_w = 2e4 kPa (M_B = K_w/n,
 *    n = (v0 - 1)/v0), eps_a up to 10% in 400 increments.
 *
 * Material (Table 1): M = 1.02, lambda = 0.2, kappa = 0.05, v_lambda = 3.32, G = 250 kPa, porous law. The state
 * is homogeneous, so the element reproduces the axisymmetric Q8-Q4 element of the previous version of the
 * article and the Python code (gen_data.py itasca); the spread of p' and q over the eight points is reported.
 *
 * The closed-form solutions of Appendix B (drained: B.1 with constant G; undrained: constant volume, B.7) and
 * the critical state values of Table 5 are evaluated at the final axial strain, without interpolation
 * (bisection on the stress ratio).
 *
 * Sect. 6.7 (Table 10): the drained test with R = 1.6 to eps_a = 5% in 50 increments is repeated with the
 * consistent tangent D, central differences, the symmetric part (D + D^T)/2, the continuum tangent and D^T
 * (TPZPlasticStepModifiedCamClay::SetTangentMode); the work counters of TPZPoroElastoPlasticUPAnalysis
 * (NGlobalIterations, NBisections) and the wall time are reported.
 *
 * The class follows the structure of the NeoPZ examples (footing.h): geometric mesh, computational meshes
 * (displacement, pore pressure and multiphysics), analysis with structural matrix and solver, incremental
 * solution and post-processing (CSV files of the figures and tables, VTK files).
 */
class FLAC3DTriaxial {
public:
    /** @brief Tangent operator returned by the stress update (Sect. 6.7) */
    typedef TPZPlasticStepModifiedCamClay::ETangentMode ETangentMode;

    /**
     * @brief Material and boundary ids: faces x = 0, x = 1, y = 0, y = 1, z = 0, z = 1 (displacement and cell
     * pressure) and the drained faces (pore pressure, coincident with the six faces in the drained tests)
     */
    enum { EMatId = 1, EX0 = -21, EX1 = -22, EY0 = -23, EY1 = -24, EZ0 = -25, EZ1 = -26, EPDrained = -30 };

    /** @brief Definition of a test */
    struct TCase {
        REAL fR = 1.6;            ///< overconsolidation ratio p'c0/p'0
        bool fDrained = true;     ///< drained or undrained test
        REAL fEaMax = 0.5;        ///< final axial strain
        int fNSteps = 500;        ///< number of increments
        ETangentMode fTangent = TPZPlasticStepModifiedCamClay::EConsistentTangent; ///< tangent operator
        std::string fName = "drained_R1.6"; ///< name of the output files
        bool fOutput = true;      ///< write the history (CSV) and the VTK files
    };

    /** @brief Columns of a row of the history */
    enum { EEa = 0, EP = 1, EQ = 2, EV = 3, EU = 4, EEvaluations = 5 };

    /** @brief History and work of a test */
    struct TResult {
        /// rows (eps_a, p', q, v, mean pore pressure of the vertices, residual evaluations of the increment)
        std::vector<std::array<REAL, 6>> fHistory;
        REAL fV0 = 0.;                    ///< initial specific volume v0 (on the normal compression line)
        REAL fPc0 = 0.;                   ///< initial preconsolidation pressure p'c0 = R p'0 (kPa)
        REAL fMeanEvaluations = 0.;       ///< mean residual evaluations per converged increment (step log)
        int fMaxEvaluations = 0;          ///< largest number of evaluations of an increment
        int64_t fNGlobalIterations = 0;   ///< all the global iterations, failed attempts included
        int64_t fNBisections = 0;         ///< bisected increments
        REAL fWallTime = 0.;              ///< wall time of the incremental solution (s)
        bool fCompleted = false;          ///< all the increments converged
        REAL fSpreadP = 0.;               ///< largest spread of p' over the integration points (kPa)
        REAL fSpreadQ = 0.;               ///< largest spread of q over the integration points (kPa)
        REAL fSpreadU = 0.;               ///< largest spread of the pore pressure over the vertices (kPa)
        int64_t fNEquations = 0;          ///< equations of the multiphysics mesh
        int fNPoints = 0;                 ///< integration points
    };

    /** @brief Values of Table 5 for a test: this work, closed form, FLAC3D and critical state */
    struct TTable5 {
        std::array<REAL, 3> fThis{};     ///< this work: (p', q, v) drained or (p', q, u) undrained
        std::array<REAL, 3> fClosed{};   ///< closed form at the same axial strain (same quantities)
        std::array<REAL, 3> fFLAC{};     ///< FLAC3D (Table 5 of the article)
        std::array<REAL, 3> fCritical{}; ///< critical state reached by the test
    };

    /** @name Material parameters (Table 1) and loading */
    /** @{ */
    REAL fM = 1.02;      ///< slope of the critical state line M
    REAL fLambda = 0.2;  ///< slope of the normal compression line lambda
    REAL fKappa = 0.05;  ///< slope of the swelling lines kappa
    REAL fVLambda = 3.32; ///< specific volume of the normal compression line at p' = 1 kPa
    REAL fG = 250.;      ///< shear modulus G (kPa, constant)
    REAL fP0 = 5.;       ///< initial isotropic effective stress and cell pressure p'0 (kPa)
    REAL fKw = 2.e4;     ///< bulk modulus of the pore fluid K_w (kPa, undrained tests)
    /** @} */

    /** @brief Specific volume of the initial state on the normal compression line */
    REAL V0(REAL R) const { return mcc::SpecificVolumeNCL(fVLambda, fLambda, fKappa, R * fP0, fP0); }

    /** @brief Geometric mesh: the unit cube with the boundary quadrilaterals (and the drained faces) */
    TPZGeoMesh *CreateGeoMesh(bool drained);

    /** @brief Computational meshes of displacement (Hex20), pore pressure (Hex8) and the multiphysics mesh */
    TPZMultiphysicsCompMesh *CreateCompMesh(TPZGeoMesh *gmesh, const TCase &tcase, mcc::TPoroMaterial *&mat);

    /** @brief Runs a test: analysis, structural matrix, solver, increments and post-processing */
    TResult Run(const TCase &tcase);

    /** @name Closed-form solutions (Appendix B) and critical state */
    /** @{ */
    /**
     * @brief Drained test (B.1)-(B.6) with constant G at the axial strain ea: elastic branch (q = 3(p' - p'0),
     * bisection on p') or plastic branch (bisection on the stress ratio eta between eta_y and M)
     * @return (p', q, v = v0 (1 - eps_v), eps_v)
     */
    std::array<REAL, 4> DrainedClosed(REAL R, REAL ea) const;
    /** @brief Axial strain, p', q and eps_v of the drained closed form on the plastic branch at the stress ratio eta */
    std::array<REAL, 4> DrainedClosedAtRatio(REAL R, REAL eta) const;
    /** @brief Yield point of the drained path: returns eta_y, py is p' at yield */
    REAL DrainedYieldRatio(REAL R, REAL &py) const;
    /**
     * @brief Undrained test (constant volume, B.7) at the axial strain ea:
     * \f$p'=p'_0((M^2+\eta^2)/(M^2R))^{-\Lambda}\f$, \f$\varepsilon_a=q/(3G)+2\Lambda\kappa/(Mv_0)(g(\eta)-g(\eta_y))\f$
     * with \f$g(x)=\frac12\ln|(M+x)/(M-x)|-\arctan(x/M)\f$, \f$\eta_y=M\sqrt{R-1}\f$ (the expressions of figs.py)
     * @return (p', q, u = q/3 + p'0 - p')
     */
    std::array<REAL, 3> UndrainedClosed(REAL R, REAL ea) const;
    /** @brief Axial strain, p', q and u of the undrained closed form on the plastic branch at the stress ratio eta */
    std::array<REAL, 4> UndrainedClosedAtRatio(REAL R, REAL eta) const;
    /**
     * @brief Critical state reached by the test (last column of Table 5): drained, \f$p'=3p'_0/(3-M)\f$, \f$q=Mp'\f$,
     * \f$v=\Gamma-\lambda\ln p'\f$ with \f$\Gamma=v_\lambda-(\lambda-\kappa)\ln 2\f$; undrained, \f$p'=p'_0(2/R)^{-\Lambda}\f$,
     * \f$q=Mp'\f$, \f$u=q/3+p'_0-p'\f$
     */
    std::array<REAL, 3> CriticalState(REAL R, bool drained) const;
    /**
     * @brief Peak of the closed form (largest q and its axial strain) in the supercritical region, by dense
     * sampling of the stress ratio; NaN in the subcritical region, where q grows monotonically
     */
    std::pair<REAL, REAL> ClosedPeak(REAL R, bool drained) const;
    /** @} */

    /**
     * @brief Writes the closed-form curves of Fig. 7 up to the final axial strain of the test (dashed lines):
     * drained, mcc::TriaxialDrainedClosed with npts = 600 (columns eps_a, p_eff, q, v); undrained, 40 points on the
     * elastic branch and the stress ratios of fig_itasca of figs.py (npts = 600, clustered near M; columns eps_a,
     * p_eff, q, u)
     */
    void WriteClosedForm(const TCase &tcase, const std::string &file) const;

    /** @brief Runs the four tests of Table 5 (Fig. 7) and writes the tables */
    void RunTests();

    /** @brief Runs the drained test with R = 1.6 in 50 increments with the five tangent operators (Table 10) */
    void RunTangents();

    /** @brief Writes the geometric meshes of the drained and undrained tests (model of Fig. 4) */
    void WriteMeshes();

    /** @brief Runs everything */
    void RunAll();

private:
    /** @brief Writes a CSV file with a header and rows of text cells */
    static void WriteTextCSV(const std::string &file, const std::string &header,
                             const std::vector<std::string> &rows) {
        std::ofstream out(file);
        out << header << "\n";
        for (auto &r : rows) out << r << "\n";
    }
    /** @brief Bisection for f(x) = 0 on [a, b] with f(a) < 0 < f(b) or f(a) > 0 > f(b) */
    template <class TFunc>
    static REAL Bisection(TFunc f, REAL a, REAL b) {
        const REAL fa = f(a);
        for (int it = 0; it < 300; ++it) {
            const REAL c = 0.5 * (a + b);
            if (c == a || c == b) break;
            if ((f(c) < 0.) == (fa < 0.))
                a = c;
            else
                b = c;
        }
        return 0.5 * (a + b);
    }
    /** @brief Function F of eq. (B.5) (drained closed form) */
    REAL F(REAL x) const {
        const REAL M = fM;
        return (1. / M) * std::log(std::fabs((M + x) / (M - x))) - (2. / M) * std::atan(x / M) -
               std::log(std::fabs(M - x)) / (3. - M) - std::log(M + x) / (3. + M) + 6. * std::log(3. - x) / (9. - M * M);
    }
    /** @brief Function g of the undrained closed form */
    REAL g(REAL x) const { return 0.5 * std::log(std::fabs((fM + x) / (fM - x))) - std::atan(x / fM); }
};

inline TPZGeoMesh *FLAC3DTriaxial::CreateGeoMesh(bool drained) {
    // one trilinear hexahedron on [0,1]^3 (the serendipity displacement space needs no geometric mid-edge nodes);
    // each face receives a quadrilateral for the displacement condition and, in the drained tests, a coincident
    // quadrilateral for the pore pressure
    return mcc::CreateUnitCubeMesh(EMatId, [drained](const std::array<std::array<REAL, 3>, 4> &X) {
        std::vector<int> ids;
        const int planes[3][2] = {{EX0, EX1}, {EY0, EY1}, {EZ0, EZ1}};
        for (int axis = 0; axis < 3; ++axis)
            for (int side = 0; side < 2; ++side)
                if (mcc::FaceOnPlane(X, axis, REAL(side))) ids.push_back(planes[axis][side]);
        if (drained) ids.push_back(EPDrained);
        return ids;
    });
}

inline TPZMultiphysicsCompMesh *FLAC3DTriaxial::CreateCompMesh(TPZGeoMesh *gmesh, const TCase &tcase,
                                                              mcc::TPoroMaterial *&mat) {
    std::set<int> bcids = {EX0, EX1, EY0, EY1, EZ0, EZ1};
    if (tcase.fDrained) bcids.insert(EPDrained);
    TPZCompMesh *cmeshU = mcc::CreateDisplacementMesh(gmesh, 3, EMatId, bcids); // serendipity Hex20
    TPZCompMesh *cmeshP = mcc::CreatePressureMesh(gmesh, 3, EMatId, bcids);     // trilinear Hex8

    const REAL pc0 = tcase.fR * fP0;
    const REAL v0 = V0(tcase.fR);
    const REAL n0 = (v0 - 1.) / v0;

    TPZMultiphysicsCompMesh *mphys = new TPZMultiphysicsCompMesh(gmesh);
    mphys->SetDimModel(3);
    mat = new mcc::TPoroMaterial(EMatId, TPZMatPoroElastoPlasticUPBase::EThreeDimensional);
    mcc::TPlastic model;
    model.SetModifiedCamClay(fM, fLambda, fKappa);
    model.SetConstantShearModulus(fG);
    model.SetTangentMode(tcase.fTangent);
    mat->SetPlasticModel(model);
    mat->SetBiot(1., tcase.fDrained ? 0. : n0 / fKw); // 1/M_B = n/K_w in the undrained tests
    mat->SetPermeability(0.);                         // no flow (and Dt = 0 in all the increments)
    mat->SetIntegrationOrder(3);                      // 2 x 2 x 2 Gauss points
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
    // faces x = 1 and y = 1: total cell pressure (traction -p'0 n)
    for (int k = 0; k < 2; ++k) {
        val1.Zero();
        val2.Fill(0.);
        val2[k] = -fP0;
        mphys->InsertMaterialObject(mat->CreateBC(mat, k == 0 ? EX1 : EY1, B::ENeumannU, val1, val2));
    }
    if (tcase.fDrained) {
        TPZManVector<STATE, 3> zero(1, 0.);
        val1.Zero();
        mphys->InsertMaterialObject(mat->CreateBC(mat, EPDrained, B::EDirichletP, val1, zero));
    }
    mcc::BuildMultiphysics(mphys, cmeshU, cmeshP, EMatId, bcids);

    // initial state: isotropic effective stress p'0, preconsolidation pc0, specific volume v0, no pore pressure
    mat->InitializeMemory(mphys, [&](const TPZVec<REAL> &, TPZElastoPlasticMem &mem) {
        mem.m_sigma = mcc::IsotropicTensor(-fP0);
        mem.m_elastoplastic_state.m_hardening = pc0;
        mem.m_elastoplastic_state.fmatprop.Resize(1, v0);
        mem.m_elastoplastic_state.fmatprop[0] = v0;
        mem.m_elastoplastic_state.fpressure = 0.;
    });
    return mphys;
}

inline FLAC3DTriaxial::TResult FLAC3DTriaxial::Run(const TCase &tcase) {
    TPZGeoMesh *gmesh = CreateGeoMesh(tcase.fDrained);
    mcc::TPoroMaterial *mat = nullptr;
    TPZMultiphysicsCompMesh *mphys = CreateCompMesh(gmesh, tcase, mat);

    TResult res;
    res.fPc0 = tcase.fR * fP0;
    res.fV0 = V0(tcase.fR);
    res.fNEquations = mphys->NEquations();

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

    // increments of the top displacement (height 1 m: u_z = -eps_a), constant cell pressure, Dt = 0
    std::vector<mcc::TAnalysis::TLoadState> steps;
    for (int k = 1; k <= tcase.fNSteps; ++k) steps.emplace_back(0., 1., -tcase.fEaMax * k / tcase.fNSteps);

    // the state is homogeneous: the first integration point represents the element; the spread of p' and q over
    // the points and of the pore pressure over the vertices measures the deviation from homogeneity
    const std::set<int64_t> vertices = analysis.NodesOfMaterials({EMatId});
    size_t nlog = 0;
    auto monitor = [&](int, const mcc::TAnalysis::TLoadState &s) {
        const auto gps = mcc::GaussPoints(mat, mphys);
        const auto &g = gps[0];
        const REAL p = mcc::MeanEffectiveStress(g.fSigma), q = mcc::DeviatoricStress(g.fSigma);
        for (auto &gi : gps) {
            res.fSpreadP = std::max(res.fSpreadP, std::fabs(mcc::MeanEffectiveStress(gi.fSigma) - p));
            res.fSpreadQ = std::max(res.fSpreadQ, std::fabs(mcc::DeviatoricStress(gi.fSigma) - q));
        }
        res.fNPoints = int(gps.size());
        REAL pmean = 0., pmin = 1e300, pmax = -1e300;
        for (auto n : vertices) {
            const REAL pw = analysis.NodalValue(n, 1, 0);
            pmean += pw;
            pmin = std::min(pmin, pw);
            pmax = std::max(pmax, pw);
        }
        res.fSpreadU = std::max(res.fSpreadU, pmax - pmin);
        // residual evaluations of the increment (sum over its sub-increments if it was bisected)
        int nev = 0;
        const auto &log = analysis.StepLog();
        for (; nlog < log.size(); ++nlog) nev += int(log[nlog].fResiduals.size());
        res.fMaxEvaluations = std::max(res.fMaxEvaluations, nev);
        res.fHistory.push_back({0. - s.fUc, p, q, res.fV0 * (1. + g.fEps.I1()), pmean / vertices.size(), REAL(nev)});
    };
    analysis.ResetCounters();
    const auto start = std::chrono::steady_clock::now();
    res.fCompleted = analysis.Run(steps, monitor);
    res.fWallTime = std::chrono::duration<REAL>(std::chrono::steady_clock::now() - start).count();
    if (!res.fCompleted) std::cerr << "FLAC3DTriaxial: test " << tcase.fName << " stopped\n";
    res.fMeanEvaluations = mcc::MeanEvaluations(analysis.StepLog());
    res.fNGlobalIterations = analysis.NGlobalIterations();
    res.fNBisections = analysis.NBisections();

    if (tcase.fOutput) {
        mcc::WriteGaussPointsVTK(mat, mphys, "flac3d_" + tcase.fName + "_gauss.vtk");
        mcc::WriteNodalVTK(analysis, 3, "flac3d_" + tcase.fName + ".vtk", 0);
        std::vector<std::vector<REAL>> rows;
        for (auto &h : res.fHistory) rows.emplace_back(h.begin(), h.end());
        mcc::WriteCSV("flac3d_" + tcase.fName + ".csv", {"eps_a", "p_eff", "q", "v", "u", "evaluations"}, rows);
    }
    mcc::DeleteMeshes(mphys);
    return res;
}

inline REAL FLAC3DTriaxial::DrainedYieldRatio(REAL R, REAL &py) const {
    const REAL p0 = fP0, pc0 = R * fP0, M = fM;
    const REAL A = 9. + M * M, B = -(18. * p0 + M * M * pc0), C = 9. * p0 * p0;
    const REAL disc = std::sqrt(B * B - 4. * A * C);
    py = 1e300;
    for (REAL root : {(-B + disc) / (2. * A), (-B - disc) / (2. * A)})
        if (root >= p0 - 1e-9) py = std::min(py, root);
    return 3. * (py - p0) / py;
}

inline std::array<REAL, 4> FLAC3DTriaxial::DrainedClosedAtRatio(REAL R, REAL eta) const {
    REAL py;
    const REAL etay = DrainedYieldRatio(R, py);
    const REAL p0 = fP0, v0 = V0(R);
    const REAL pp = 3. * p0 / (3. - eta), q = eta * pp, pc = pp * (1. + eta * eta / (fM * fM));
    const REAL ev = fKappa / v0 * std::log(pp / p0) + (fLambda - fKappa) / v0 * std::log(pc / (R * p0));
    const REAL eq = q / (3. * fG) + (fLambda - fKappa) / v0 * (F(eta) - F(etay));
    return {eq + ev / 3., pp, q, ev};
}

inline std::array<REAL, 4> FLAC3DTriaxial::DrainedClosed(REAL R, REAL ea) const {
    REAL py;
    const REAL etay = DrainedYieldRatio(R, py);
    const REAL p0 = fP0, v0 = V0(R);
    auto elastic = [&](REAL pp) {
        const REAL q = 3. * (pp - p0), ev = fKappa / v0 * std::log(pp / p0);
        return std::array<REAL, 4>{q / (3. * fG) + ev / 3., pp, q, ev};
    };
    std::array<REAL, 4> s;
    if (elastic(py)[0] >= ea) {
        const REAL pp = Bisection([&](REAL x) { return elastic(x)[0] - ea; }, p0, py);
        s = elastic(pp);
    } else {
        // eps_a increases from the yield point (eta_y) to infinity (eta -> M) in both regions
        const REAL eta = Bisection([&](REAL x) { return DrainedClosedAtRatio(R, x)[0] - ea; }, etay,
                                   fM - (fM - etay) * 1e-15);
        s = DrainedClosedAtRatio(R, eta);
    }
    return {s[1], s[2], v0 * (1. - s[3]), s[3]};
}

inline std::array<REAL, 4> FLAC3DTriaxial::UndrainedClosedAtRatio(REAL R, REAL eta) const {
    const REAL Lam = (fLambda - fKappa) / fLambda, etay = fM * std::sqrt(R - 1.), v0 = V0(R);
    const REAL p = mcc::UndrainedClosedP(fP0, R, fM, fLambda, fKappa, eta), q = eta * p;
    return {q / (3. * fG) + 2. * Lam * fKappa / (fM * v0) * (g(eta) - g(etay)), p, q, q / 3. + fP0 - p};
}

inline std::array<REAL, 3> FLAC3DTriaxial::UndrainedClosed(REAL R, REAL ea) const {
    const REAL etay = fM * std::sqrt(R - 1.), qy = etay * fP0;
    if (qy / (3. * fG) >= ea) {
        const REAL q = 3. * fG * ea; // elastic: p' constant, eps_q = eps_a
        return {fP0, q, q / 3.};
    }
    const REAL eta = Bisection([&](REAL x) { return UndrainedClosedAtRatio(R, x)[0] - ea; }, etay,
                               fM - (fM - etay) * 1e-15);
    const auto s = UndrainedClosedAtRatio(R, eta);
    return {s[1], s[2], s[3]};
}

inline std::array<REAL, 3> FLAC3DTriaxial::CriticalState(REAL R, bool drained) const {
    if (drained) {
        const REAL p = 3. * fP0 / (3. - fM), Gamma = fVLambda - (fLambda - fKappa) * std::log(2.);
        return {p, fM * p, Gamma - fLambda * std::log(p)};
    }
    const REAL p = mcc::UndrainedClosedP(fP0, R, fM, fLambda, fKappa, fM), q = fM * p;
    return {p, q, q / 3. + fP0 - p};
}

inline std::pair<REAL, REAL> FLAC3DTriaxial::ClosedPeak(REAL R, bool drained) const {
    // the largest q is at the yield point or on the plastic branch: sample the stress ratio between eta_y and M;
    // in the subcritical region (eta_y < M) q grows monotonically to the critical state and there is no peak
    REAL py;
    const REAL etay = drained ? DrainedYieldRatio(R, py) : fM * std::sqrt(R - 1.);
    if (etay <= fM) return {std::numeric_limits<REAL>::quiet_NaN(), std::numeric_limits<REAL>::quiet_NaN()};
    std::pair<REAL, REAL> best(-1., 0.);
    const int n = 200000;
    for (int i = 0; i < n; ++i) {
        const REAL eta = etay + (fM - etay) * REAL(i) / n;
        const auto s = drained ? DrainedClosedAtRatio(R, eta) : UndrainedClosedAtRatio(R, eta);
        if (s[2] > best.first) best = {s[2], s[0]};
    }
    return best;
}

inline void FLAC3DTriaxial::WriteClosedForm(const TCase &tcase, const std::string &file) const {
    const REAL R = tcase.fR, v0 = V0(R);
    std::vector<std::vector<REAL>> rows;
    if (tcase.fDrained) {
        for (auto &r : mcc::TriaxialDrainedClosed(fP0, R * fP0, v0, fM, fLambda, fKappa, fG, 0., 600))
            if (r[0] <= tcase.fEaMax + 1e-12) rows.push_back({r[0], r[1], r[2], v0 * (1. - r[3])});
        mcc::WriteCSV(file, {"eps_a", "p_eff", "q", "v"}, rows);
        return;
    }
    const REAL etay = fM * std::sqrt(R - 1.), qy = etay * fP0;
    for (int i = 0; i < 40; ++i) {
        const REAL q = qy * i / 39.;
        rows.push_back({q / (3. * fG), fP0, q, q / 3.});
    }
    const int npts = 600;
    std::vector<REAL> etas;
    for (int i = 1; i < npts; ++i) etas.push_back(etay + (fM - etay) * REAL(i) / npts);
    for (int i = 1; i <= npts; ++i) etas.push_back(fM - (fM - etay) * std::exp(-14. * REAL(i) / npts));
    std::sort(etas.begin(), etas.end());
    etas.erase(std::unique(etas.begin(), etas.end()), etas.end());
    if (etay > fM) std::reverse(etas.begin(), etas.end());
    for (REAL eta : etas) {
        const auto s = UndrainedClosedAtRatio(R, eta);
        if (s[0] <= tcase.fEaMax + 1e-12) rows.push_back({s[0], s[1], s[2], s[3]});
    }
    mcc::WriteCSV(file, {"eps_a", "p_eff", "q", "u"}, rows);
}

inline void FLAC3DTriaxial::WriteMeshes() {
    for (bool drained : {true, false}) {
        TPZGeoMesh *gmesh = CreateGeoMesh(drained);
        mcc::WriteMeshCSV(gmesh, std::string("flac3d_mesh_") + (drained ? "drained" : "undrained"));
        delete gmesh;
    }
}

inline void FLAC3DTriaxial::RunTests() {
    std::cout << std::setprecision(6);
    std::cout << "FLAC3D triaxial tests (Sect. 6.3, Table 5, Fig. 7): one Hex20-Hex8 element, 2 x 2 x 2 points\n";
    std::vector<TCase> cases(4);
    cases[0] = {1.6, true, 0.5, 500, TPZPlasticStepModifiedCamClay::EConsistentTangent, "drained_R1.6"};
    cases[1] = {8.0, true, 0.5, 500, TPZPlasticStepModifiedCamClay::EConsistentTangent, "drained_R8"};
    cases[2] = {1.6, false, 0.1, 400, TPZPlasticStepModifiedCamClay::EConsistentTangent, "undrained_R1.6"};
    cases[3] = {8.0, false, 0.1, 400, TPZPlasticStepModifiedCamClay::EConsistentTangent, "undrained_R8"};
    // reference values: final states of the Python code (gen_data.py itasca, axisymmetric Q8-Q4 element, which
    // gave the numbers of v0.6) and of FLAC3D (Table 5); evaluations per increment of the Python code
    const REAL python[4][4] = {{7.573292329890411, 7.719876989665036, 2.8111971280068224, 0.0},
                               {7.5838818484038555, 7.751645545163063, 2.8105109616874806, 0.0},
                               {4.234378307859034, 4.31852679256991, 2.927399341327785, 2.205130622997605},
                               {14.047913734943828, 14.422435155178112, 2.6865536965569574, -4.240435349884442}};
    const REAL pyEvaluations[4] = {2.228, 2.34, 2.0025, 2.0025};
    const REAL flac[4][3] = {{7.573, 7.718, 2.811}, {7.583, 7.747, 2.811}, {4.234, 4.312, 2.203}, {14.05, 14.42, -4.241}};
    std::vector<std::string> summary, table5;
    std::vector<std::string> diffs;
    REAL maxClosedDrained = 0., maxClosedUndrained = 0., maxFLAC = 0.;
    std::string maxClosedUndrainedWhere;
    for (int ic = 0; ic < 4; ++ic) {
        const TCase &c = cases[ic];
        const TResult r = Run(c);
        WriteClosedForm(c, "flac3d_" + c.fName + "_closed.csv");
        const auto &f = r.fHistory.back();
        size_t kp = 0;
        for (size_t k = 0; k < r.fHistory.size(); ++k)
            if (r.fHistory[k][EQ] > r.fHistory[kp][EQ]) kp = k;
        const auto peak = ClosedPeak(c.fR, c.fDrained);
        TTable5 t;
        t.fThis = {f[EP], f[EQ], c.fDrained ? f[EV] : f[EU]};
        if (c.fDrained) {
            const auto cl = DrainedClosed(c.fR, f[EEa]);
            t.fClosed = {cl[0], cl[1], cl[2]};
        } else {
            t.fClosed = UndrainedClosed(c.fR, f[EEa]);
        }
        t.fFLAC = {flac[ic][0], flac[ic][1], flac[ic][2]};
        t.fCritical = CriticalState(c.fR, c.fDrained);
        std::ostringstream closedpeak;
        if (std::isnan(peak.first))
            closedpeak << "monotonic";
        else
            closedpeak << "closed form " << peak.first << " at " << peak.second;
        std::cout << c.fName << ": v0 = " << r.fV0 << " pc0 = " << r.fPc0 << " | eps_a = " << f[EEa] << " p' = " << f[EP]
                  << " q = " << f[EQ] << " v = " << f[EV] << " u = " << f[EU] << " | peak q = " << r.fHistory[kp][EQ]
                  << " at eps_a = " << r.fHistory[kp][EEa] << " (" << closedpeak.str()
                  << ") | evaluations per increment = " << r.fMeanEvaluations << " (max " << r.fMaxEvaluations
                  << "; Python " << pyEvaluations[ic] << "), total " << r.fNGlobalIterations << ", bisections "
                  << r.fNBisections << ", " << r.fWallTime << " s\n";
        std::cout << "   difference from the Python Q8-Q4 final state (p', q, v, u): " << std::scientific
                  << std::setprecision(2) << f[EP] - python[ic][0] << " " << f[EQ] - python[ic][1] << " "
                  << f[EV] - python[ic][2] << " " << f[EU] - python[ic][3] << std::defaultfloat << std::setprecision(6)
                  << "; spread over the " << r.fNPoints << " points: p' " << r.fSpreadP << ", q " << r.fSpreadQ
                  << ", u over the vertices " << r.fSpreadU << " kPa\n";
        const char *names[2][3] = {{"p_eff", "q", "u"}, {"p_eff", "q", "v"}};
        for (int k = 0; k < 3; ++k) {
            const REAL dc = 100. * std::fabs(t.fClosed[k] - t.fThis[k]) / std::fabs(t.fThis[k]);
            const REAL df = 100. * std::fabs(t.fFLAC[k] - t.fThis[k]) / std::fabs(t.fThis[k]);
            if (c.fDrained) maxClosedDrained = std::max(maxClosedDrained, dc);
            else if (dc > maxClosedUndrained) {
                maxClosedUndrained = dc;
                maxClosedUndrainedWhere = c.fName + " " + names[c.fDrained][k];
            }
            maxFLAC = std::max(maxFLAC, df);
            std::ostringstream s;
            s << std::setprecision(10) << c.fName << "," << names[c.fDrained][k] << "," << t.fThis[k] << ","
              << t.fClosed[k] << "," << t.fFLAC[k] << "," << t.fCritical[k] << "," << dc << "," << df;
            table5.push_back(s.str());
        }
        std::ostringstream s;
        s << std::setprecision(10) << c.fName << "," << c.fR << "," << int(c.fDrained) << "," << c.fNSteps << ","
          << f[EEa] << "," << r.fV0 << "," << r.fPc0 << "," << f[EP] << "," << f[EQ] << "," << f[EV] << "," << f[EU]
          << "," << f[EQ] / f[EP] << "," << r.fHistory[kp][EQ] << "," << r.fHistory[kp][EEa] << "," << peak.first << ","
          << peak.second << "," << r.fMeanEvaluations << "," << r.fMaxEvaluations << "," << r.fNGlobalIterations << ","
          << r.fNBisections << "," << r.fWallTime << "," << r.fSpreadP << "," << r.fSpreadQ << "," << r.fSpreadU << ","
          << r.fNPoints << "," << r.fNEquations << "," << pyEvaluations[ic] << "," << f[EP] - python[ic][0] << ","
          << f[EQ] - python[ic][1] << "," << f[EV] - python[ic][2] << "," << f[EU] - python[ic][3] << "," << fM << ","
          << fLambda << "," << fKappa << "," << fVLambda << "," << fG << "," << fP0 << "," << fKw;
        summary.push_back(s.str());
    }
    std::cout << "Table 5: largest difference closed form - this work: drained " << maxClosedDrained
              << "% (v0.6: 0.02%), undrained " << maxClosedUndrained << "% in " << maxClosedUndrainedWhere
              << " (v0.6: 0.6% in u for R = 8); FLAC3D - this work: " << maxFLAC << "% (v0.6: 0.2%)\n";
    WriteTextCSV("flac3d_table5.csv",
                 "test,quantity,this_work,closed_form,flac3d,critical_state,diff_closed_percent,diff_flac3d_percent",
                 table5);
    WriteTextCSV("flac3d_summary.csv",
                 "test,R,drained,nsteps,eps_a_end,v0,pc0,p_end,q_end,v_end,u_end,eta_end,q_peak,eps_a_peak,"
                 "q_peak_closed,eps_a_peak_closed,mean_evaluations,max_evaluations,global_iterations,bisections,"
                 "wall_time_s,spread_p,spread_q,spread_u,points,equations,python_mean_evaluations,diff_python_p,"
                 "diff_python_q,diff_python_v,diff_python_u,M,lambda,kappa,v_lambda,G,p0,Kw",
                 summary);
}

inline void FLAC3DTriaxial::RunTangents() {
    using P = TPZPlasticStepModifiedCamClay;
    std::cout << "\nSect. 6.7, Table 10 (FLAC3D column): drained test, R = 1.6, 50 increments to eps_a = 5%\n";
    const ETangentMode modes[5] = {P::EConsistentTangent, P::EFiniteDifferenceTangent, P::ESymmetricTangent,
                                   P::EContinuumTangent, P::ETransposedTangent};
    const REAL python[5] = {3.02, 3.02, 3.02, 5.48, 3.02};
    std::vector<TResult> res;
    std::vector<std::string> rows;
    for (int m = 0; m < 5; ++m) {
        TCase c = {1.6, true, 0.05, 50, modes[m], std::string("drained_R1.6_50steps_") + P::TangentModeName(modes[m])};
        c.fOutput = (m == 0);
        res.push_back(Run(c));
        const TResult &r = res.back();
        const auto &f = r.fHistory.back();
        std::cout << "  " << std::setw(4) << P::TangentModeName(modes[m]) << ": " << std::fixed << std::setprecision(2)
                  << r.fMeanEvaluations << " evaluations per increment (max " << r.fMaxEvaluations << "), total "
                  << r.fNGlobalIterations << ", bisections " << r.fNBisections << ", " << std::setprecision(4)
                  << r.fWallTime << " s; final q = " << std::setprecision(10) << f[EQ] << std::defaultfloat
                  << std::setprecision(6) << " kPa (v0.6, Python: " << python[m] << ")\n";
        std::ostringstream s;
        s << std::setprecision(10) << P::TangentModeName(modes[m]) << "," << c.fNSteps << "," << r.fMeanEvaluations
          << "," << r.fMaxEvaluations << "," << r.fNGlobalIterations << "," << r.fNBisections << "," << r.fWallTime << ","
          << int(r.fCompleted) << "," << f[EP] << "," << f[EQ] << "," << f[EV] << "," << python[m];
        rows.push_back(s.str());
    }
    REAL dq = 0.;
    for (auto &r : res) dq = std::max(dq, std::fabs(r.fHistory.back()[EQ] - res[0].fHistory.back()[EQ]));
    std::cout << "  largest difference of the final q with respect to D: " << dq << " kPa\n";
    WriteTextCSV("flac3d_table10.csv",
                 "tangent,nsteps,mean_evaluations,max_evaluations,global_iterations,bisections,wall_time_s,completed,"
                 "p_end,q_end,v_end,python_mean_evaluations",
                 rows);
    // evaluations of each increment with each operator
    std::vector<std::vector<REAL>> ev;
    for (size_t k = 1; k < res[0].fHistory.size(); ++k) {
        ev.push_back({REAL(k), res[0].fHistory[k][EEa]});
        for (auto &r : res) ev.back().push_back(k < r.fHistory.size() ? r.fHistory[k][EEvaluations] : 0.);
    }
    std::vector<std::string> header = {"increment", "eps_a"};
    for (auto m : modes) header.push_back(P::TangentModeName(m));
    mcc::WriteCSV("flac3d_tangents_evaluations.csv", header, ev);
}

inline void FLAC3DTriaxial::RunAll() {
    const auto start = std::chrono::steady_clock::now();
    WriteMeshes();
    RunTests();
    RunTangents();
    std::cout << "\nFiles: flac3d_<test>.csv and flac3d_<test>_closed.csv (Fig. 7), flac3d_table5.csv, "
                 "flac3d_summary.csv, flac3d_table10.csv, flac3d_tangents_evaluations.csv, "
                 "flac3d_mesh_<drained|undrained>_*.csv (Fig. 4), VTK files\n"
              << "total run time " << std::chrono::duration<REAL>(std::chrono::steady_clock::now() - start).count()
              << " s\n";
}
