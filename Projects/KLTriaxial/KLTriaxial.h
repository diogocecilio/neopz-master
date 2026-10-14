/**
 * @file KLTriaxial.h
 * @brief Sect. 7.1 of the article (v1.0): triaxial tests of example 1 of Krabbenhoft and Lyamin (2012), drained
 * and undrained, normally consolidated and with OCR = 10, with one Hex20-Hex8 u-p element against the closed forms
 * of Appendix B and the material point solution, and the convergence of the error with the number of increments.
 */
#pragma once

#include "MCCPaperTools.h"
#include "IntegrationSchemes.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "pzstepsolver.h"

#include <array>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <set>
#include <string>
#include <vector>

/**
 * @ingroup mccpaper
 * @brief Triaxial compression tests of example 1 of Krabbenhoft and Lyamin (2012) with the return mapping of this
 * work (TPZPlasticStepModifiedCamClay) in one Hex20-Hex8 u-p element.
 *
 * Material: \f$M=3\sin\phi/\sqrt{3+\sin^2\phi}\f$ with \f$\phi=24^\circ\f$ (0.6858), \f$\lambda=0.07\f$,
 * \f$\kappa=0.008\f$, \f$v_0=1.5\f$, porous elasticity with constant Poisson's ratio \f$\nu=0.3\f$, \f$p_t=0\f$,
 * \f$\omega=1\f$, initial isotropic effective stress \f$p'_0=200\f$ kPa (IntegrationSchemes::KLMaterial, the clay
 * of the material-point comparison of the integration schemes). Four tests:
 *  - drained, NC (\f$p'_{c0}=200\f$ kPa) and OCR = 10 (\f$p'_{c0}=2000\f$ kPa), cell pressure \f$p'_0\f$, axial
 *    strain to 25% (the OCR = 10 clay yields in the supercritical region: peak of q, softening and dilation);
 *  - undrained, NC to \f$\varepsilon_a=2.5\%\f$ and OCR = 10 to 8% (strain control, \f$\varepsilon_v=0\f$).
 *
 * For each test the example computes
 *  -# the closed form of Appendix B: the drained path (mcc::TriaxialDrainedClosed, 20000 points) and the undrained
 *     path (IntegrationSchemes::UndrainedClosed, eq. (B.7) with the elastic shear strain (B.9)), with the reference
 *     values at the final axial strain (interpolation of the drained table; Brent's method on the stress ratio of
 *     the undrained path) and the peak of q of the drained tests;
 *  -# the material point solution (mcc::TriaxialDrained with the xx-yy block of the consistent tangent for the
 *     lateral strain, mcc::TriaxialUndrained) with 10, 100 and 1000 increments;
 *  -# the finite element solution with one Hex20-Hex8 element on the unit cube (2 x 2 x 2 Gauss points): u_x = 0 on
 *     x = 0, u_y = 0 on y = 0, u_z = 0 on z = 0, total traction \f$-p'_0\f$ on the faces x = 1 and y = 1, vertical
 *     displacement of the top z = 1 controlled by the analysis; drained tests: pore pressure prescribed as zero at
 *     the eight vertices (permeability zero, incompressible constituents); undrained tests: no drainage, the pore
 *     pressure is the Lagrange multiplier of the incompressibility (\f$\Delta t=0\f$); the same increments as the
 *     material point. The deformation is homogeneous, so the element reproduces the material point to round-off,
 *     and the error of both against the closed form is the error of the integration scheme with the given
 *     increments (first order, Sect. 5).
 *
 * Output files (working directory), with \<case\> = drained_nc, drained_ocr10, undrained_nc, undrained_ocr10:
 *  - kl_\<case\>_closed.csv: closed form (eps_a, p', q, eps_v or excess pore pressure u);
 *  - kl_\<case\>_fe_n\<n\>.csv: finite element history (eps_a, p', q, eps_v, u, evaluations of the increment) with
 *    the material point columns (p', q, eps_v or u) and the differences;
 *  - kl_table.csv: errors at the end of the tests against the number of increments (finite element and material
 *    point), largest difference between them along the path, global evaluations per increment, bisections, wall
 *    time, local Newton iterations (mean and largest) of the material point;
 *  - kl_fe.vtk files of the last state of the 100-increment runs.
 *
 * Stresses in kPa, compression positive in the printed and written invariants.
 */
class KLTriaxial {
public:
    /** @brief Material and boundary ids of the single element (as RS2Triaxial) */
    enum { EMatId = 1, EX0 = -21, EX1 = -22, EY0 = -23, EY1 = -24, EZ0 = -25, EZ1 = -26, EPDrained = -30 };

    /** @brief Definition of a test */
    struct TCase {
        std::string fTag;   ///< tag of the output files
        std::string fLabel; ///< label of the figure
        bool fDrained;      ///< drained (cell pressure constant) or undrained (eps_v = 0)
        int fOCR;           ///< overconsolidation ratio p'_c0 / p'_0
        REAL fEaMax;        ///< final axial strain
    };

    /** @brief One finite element (and material point) run with n increments */
    struct TRun {
        int fN = 0;              ///< increments
        bool fConverged = false; ///< the finite element solution completed
        /** @brief rows (eps_a, p', q, eps_v, u, evaluations of the increment), the initial state first */
        std::vector<std::array<REAL, 6>> fHistory;
        std::vector<std::array<REAL, 4>> fPoint; ///< material point (eps_a, p', q, eps_v or u)
        REAL fErrQ = 0.;         ///< q_FE - q_exact at the end
        REAL fErrP = 0.;         ///< p'_FE - p'_exact at the end
        REAL fErrEv = 0.;        ///< eps_v,FE - eps_v,exact at the end (drained)
        REAL fPointErrQ = 0.;    ///< q_point - q_exact at the end
        REAL fMaxDiffPoint = 0.; ///< largest |q_FE - q_point| along the path
        REAL fMeanEval = 0.;     ///< mean global evaluations per increment
        int fMaxEval = 0;        ///< largest number of evaluations of an increment
        int64_t fGlobalIts = 0;  ///< global iterations (failed attempts included)
        int64_t fBisections = 0; ///< bisected increments
        REAL fWall = 0.;         ///< wall time of the incremental solution (s)
        REAL fSpreadQ = 0.;      ///< largest spread of q over the 8 points (kPa)
        REAL fQPeak = 0.;        ///< largest q of the finite element path
        REAL fEaPeak = 0.;       ///< axial strain of the largest q
        REAL fUEnd = 0.;         ///< nodal pore pressure at the end (undrained)
        mcc::TLocalStats fStats; ///< local Newton iterations of the material point run
    };

    /** @brief Results of a test */
    struct TResult {
        TCase fCase;
        std::vector<std::array<REAL, 4>> fClosed; ///< closed form (eps_a, p', q, eps_v or u)
        REAL fQRef = 0., fPRef = 0., fEvRef = 0.; ///< closed form at the end of the test
        REAL fQPeakRef = 0., fEaPeakRef = 0.;     ///< peak of q of the closed form (drained)
        std::vector<TRun> fRuns;                  ///< runs with the numbers of increments
    };

    /** @brief Initial isotropic effective stress p'0 (kPa), also the cell pressure */
    REAL fP0 = 200.;
    /** @brief Numbers of increments of the convergence study */
    std::vector<int> fIncrements = {10, 100, 1000};
    /** @brief Number of increments of the curves of the figure (VTK files written) */
    int fNFigure = 100;
    /** @brief Runs the finite element solutions (false: material point only) */
    bool fFECheck = true;
    /** @brief Tangent operator returned to the global iterations (TPZPlasticStepModifiedCamClay::SetTangentMode) */
    TPZPlasticStepModifiedCamClay::ETangentMode fTangent = TPZPlasticStepModifiedCamClay::EConsistentTangent;
    /** @brief Number of increments of the comparison of the tangent operators (RunTangents) */
    int fNTangents = 100;

    /** @brief The four tests */
    std::vector<TCase> Cases() const;
    /** @brief Clay of example 1 of Krabbenhoft and Lyamin (IntegrationSchemes::KLMaterial) */
    static IntegrationSchemes::TMaterial Material() { return IntegrationSchemes::KLMaterial(); }
    /** @brief Plastic model of the library with the parameters of Material() */
    mcc::TPlastic CreateMaterial() const;
    /** @brief Closed form of a test and its reference values (Appendix B) */
    void ClosedForm(const TCase &c, TResult &res) const;
    /** @brief Material point solution with n increments (eps_a, p', q, eps_v or u) */
    std::vector<std::array<REAL, 4>> MaterialPoint(const TCase &c, int n, mcc::TLocalStats *stats) const;
    /** @brief Geometric mesh of the element (the drained faces exist only in the drained tests) */
    TPZGeoMesh *CreateGeoMesh(bool drained) const;
    /** @brief Displacement (Hex20), pore pressure (Hex8) and multiphysics meshes, materials and initial state */
    TPZMultiphysicsCompMesh *CreateCompMesh(TPZGeoMesh *gmesh, const TCase &c, mcc::TPoroMaterial *&mat) const;
    /** @brief Finite element solution with n increments (the material point columns of run must be filled) */
    void RunFiniteElement(const TCase &c, TRun &run) const;
    /** @brief Runs a test: closed form, material point and finite element solutions, CSV files */
    TResult Run(const TCase &c) const;
    /** @brief Runs the four tests, prints the table and writes kl_table.csv */
    void RunAll();
    /**
     * @brief Global iterations of the drained tests (NC and OCR = 10, fNTangents increments) with the tangent
     * operators D, central differences, (D + D^T)/2, the continuum tangent and D^T (Table 8 of the article, v1.1);
     * writes kl_tangents.csv
     */
    void RunTangents();
};

inline std::vector<KLTriaxial::TCase> KLTriaxial::Cases() const {
    return {{"drained_nc", "drained, NC", true, 1, 0.25},
            {"drained_ocr10", "drained, OCR = 10", true, 10, 0.25},
            {"undrained_nc", "undrained, NC", false, 1, 0.025},
            {"undrained_ocr10", "undrained, OCR = 10", false, 10, 0.08}};
}

inline mcc::TPlastic KLTriaxial::CreateMaterial() const {
    const IntegrationSchemes::TMaterial P = Material();
    mcc::TPlastic model;
    model.SetModifiedCamClay(P.fM, P.fLambda, P.fKappa);
    model.SetPorousElasticity();
    model.SetPoissonRatio(P.fNu);
    model.SetDefaultSpecificVolume(P.fV0);
    model.SetTangentMode(fTangent);
    return model;
}

inline void KLTriaxial::ClosedForm(const TCase &c, TResult &res) const {
    const IntegrationSchemes::TMaterial P = Material();
    const REAL p0 = fP0, pc0 = c.fOCR * p0;
    res.fClosed.clear();
    if (c.fDrained) {
        const auto ana = mcc::TriaxialDrainedClosed(p0, pc0, P.fV0, P.fM, P.fLambda, P.fKappa, 0., P.fNu, 20000);
        res.fQRef = mcc::Interpolate(ana, c.fEaMax, 2);
        res.fPRef = mcc::Interpolate(ana, c.fEaMax, 1);
        res.fEvRef = mcc::Interpolate(ana, c.fEaMax, 3);
        res.fQPeakRef = -1.;
        for (auto &r : ana) {
            if (r[0] <= c.fEaMax + 1e-12) res.fClosed.push_back({r[0], r[1], r[2], r[3]});
            if (r[2] > res.fQPeakRef) {
                res.fQPeakRef = r[2];
                res.fEaPeakRef = r[0];
            }
        }
        return;
    }
    // undrained: elastic branch with p' = p'0 up to eta_y = M sqrt(R - 1), then eq. (B.7) with the elastic shear
    // strain (B.9); the stress ratio runs from eta_y to M (increasing for NC, decreasing for OCR = 10)
    const REAL M = P.fM, R = pc0 / p0;
    const REAL etay = M * std::sqrt(std::max(R - 1., 0.));
    const REAL G0 = P.GFactor() * P.fV0 * p0 / P.fKappa;
    auto row = [&](REAL eta) {
        REAL p, q, eq;
        IntegrationSchemes::UndrainedClosed(p0, pc0, P, eta, p, q, eq);
        return std::array<REAL, 4>{eq, p, q, p0 + q / 3. - p}; // u = total mean stress - p'
    };
    const int ne = 200;
    if (c.fOCR > 1)
        for (int k = 0; k < ne; ++k) {
            const REAL q = etay * p0 * k / ne;
            res.fClosed.push_back({q / (3. * G0), p0, q, q / 3.});
        }
    const int np = 4000;
    const REAL lo = c.fOCR > 1 ? M * (1. + 1e-7) : 1e-9, hi = c.fOCR > 1 ? etay : M * (1. - 1e-7);
    std::vector<std::array<REAL, 4>> plastic;
    for (int k = 0; k <= np; ++k) {
        // geometric spacing towards the critical state (eq grows without bound there)
        const REAL t = REAL(k) / np;
        const REAL eta = c.fOCR > 1 ? hi - (hi - lo) * (1. - std::pow(1. - t, 4.)) : lo + (hi - lo) * (1. - std::pow(1. - t, 4.));
        plastic.push_back(row(eta));
    }
    std::sort(plastic.begin(), plastic.end(), [](const std::array<REAL, 4> &a, const std::array<REAL, 4> &b) {
        return a[0] < b[0];
    });
    for (auto &r : plastic)
        if (r[0] <= c.fEaMax + 1e-12) res.fClosed.push_back(r);
    const REAL eta = IntegrationSchemes::Brent([&](REAL x) { return row(x)[0] - c.fEaMax; },
                                               c.fOCR > 1 ? hi : lo, c.fOCR > 1 ? lo : hi, 1e-15);
    const auto ref = row(eta);
    res.fPRef = ref[1];
    res.fQRef = ref[2];
    res.fEvRef = 0.;
    res.fClosed.push_back({c.fEaMax, ref[1], ref[2], ref[3]});
    res.fQPeakRef = std::numeric_limits<REAL>::quiet_NaN();
    res.fEaPeakRef = std::numeric_limits<REAL>::quiet_NaN();
}

inline std::vector<std::array<REAL, 4>> KLTriaxial::MaterialPoint(const TCase &c, int n, mcc::TLocalStats *stats) const {
    const mcc::TPlastic model = CreateMaterial();
    const IntegrationSchemes::TMaterial P = Material();
    const REAL p0 = fP0, pc0 = c.fOCR * p0;
    std::vector<std::array<REAL, 4>> out;
    if (c.fDrained) {
        for (auto &r : mcc::TriaxialDrained(model, p0, pc0, P.fV0, c.fEaMax, n, stats)) out.push_back({r[0], r[1], r[2], r[3]});
    } else {
        for (auto &r : mcc::TriaxialUndrained(model, p0, pc0, P.fV0, c.fEaMax, n, stats))
            out.push_back({r[0], r[1], r[2], p0 + r[2] / 3. - r[1]});
    }
    return out;
}

inline TPZGeoMesh *KLTriaxial::CreateGeoMesh(bool drained) const {
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

inline TPZMultiphysicsCompMesh *KLTriaxial::CreateCompMesh(TPZGeoMesh *gmesh, const TCase &c,
                                                           mcc::TPoroMaterial *&mat) const {
    std::set<int> bcids = {EX0, EX1, EY0, EY1, EZ0, EZ1};
    if (c.fDrained) bcids.insert(EPDrained);
    TPZCompMesh *cmeshU = mcc::CreateDisplacementMesh(gmesh, 3, EMatId, bcids); // serendipity Hex20
    TPZCompMesh *cmeshP = mcc::CreatePressureMesh(gmesh, 3, EMatId, bcids);     // trilinear Hex8

    TPZMultiphysicsCompMesh *mphys = new TPZMultiphysicsCompMesh(gmesh);
    mphys->SetDimModel(3);
    mat = new mcc::TPoroMaterial(EMatId, TPZMatPoroElastoPlasticUPBase::EThreeDimensional);
    mat->SetPlasticModel(CreateMaterial());
    mat->SetBiot(1., 0.); // incompressible constituents
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
    // faces x = 1 and y = 1: cell pressure p'0 (total stress)
    for (int k = 0; k < 2; ++k) {
        val1.Zero();
        val2.Fill(0.);
        val2[k] = -fP0;
        mphys->InsertMaterialObject(mat->CreateBC(mat, k == 0 ? EX1 : EY1, B::ENeumannU, val1, val2));
    }
    if (c.fDrained) { // p = 0 on the whole boundary (the eight vertices)
        TPZManVector<STATE, 3> zero(1, 0.);
        val1.Zero();
        mphys->InsertMaterialObject(mat->CreateBC(mat, EPDrained, B::EDirichletP, val1, zero));
    }
    mcc::BuildMultiphysics(mphys, cmeshU, cmeshP, EMatId, bcids);

    const REAL v0 = Material().fV0, pc0 = c.fOCR * fP0, p0 = fP0;
    mat->InitializeMemory(mphys, [&](const TPZVec<REAL> &, TPZElastoPlasticMem &mem) {
        mem.m_sigma = mcc::IsotropicTensor(-p0);
        mem.m_elastoplastic_state.m_hardening = pc0;
        mem.m_elastoplastic_state.fmatprop.Resize(1, v0);
        mem.m_elastoplastic_state.fmatprop[0] = v0;
        mem.m_elastoplastic_state.fpressure = 0.;
    });
    return mphys;
}

inline void KLTriaxial::RunFiniteElement(const TCase &c, TRun &run) const {
    TPZGeoMesh *gmesh = CreateGeoMesh(c.fDrained);
    mcc::TPoroMaterial *mat = nullptr;
    TPZMultiphysicsCompMesh *mphys = CreateCompMesh(gmesh, c, mat);

    mcc::TAnalysis analysis(mphys, mat);
    TPZSkylineNSymStructMatrix<STATE> skyl(mphys);
    skyl.SetNumThreads(0);
    analysis.SetStructuralMatrix(skyl);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELU);
    analysis.SetSolver(step);
    analysis.SetControlledDisplacement(EZ1, 2);
    analysis.SetPredictor(true);

    // increments of the top displacement (height 1 m: u_z = -eps_a), constant cell pressure, no time (dt = 0)
    std::vector<mcc::TAnalysis::TLoadState> steps;
    for (int k = 1; k <= run.fN; ++k) steps.emplace_back(0., 1., -c.fEaMax * k / run.fN);

    const std::set<int64_t> pnodes = analysis.NodesOfMaterials({EMatId});
    const int64_t pnode = *pnodes.begin();
    run.fHistory.clear();
    run.fSpreadQ = 0.;
    run.fMaxEval = 0;
    size_t nlog = 0;
    auto monitor = [&](int, const mcc::TAnalysis::TLoadState &s) {
        const auto gps = mcc::GaussPoints(mat, mphys);
        const auto &g = gps[0];
        const REAL q = mcc::DeviatoricStress(g.fSigma), p = mcc::MeanEffectiveStress(g.fSigma);
        for (auto &gi : gps) run.fSpreadQ = std::max(run.fSpreadQ, std::fabs(mcc::DeviatoricStress(gi.fSigma) - q));
        int nev = 0;
        const auto &log = analysis.StepLog();
        for (; nlog < log.size(); ++nlog) nev += int(log[nlog].fResiduals.size());
        run.fMaxEval = std::max(run.fMaxEval, nev);
        run.fUEnd = analysis.NodalValue(pnode, 1, 0);
        // excess pore pressure from the stresses: total mean stress p'0 + q/3 minus p' (homogeneous state)
        run.fHistory.push_back({0. - s.fUc, p, q, 0. - g.fEps.I1(), fP0 + q / 3. - p, REAL(nev)});
    };
    analysis.ResetCounters();
    const auto start = std::chrono::steady_clock::now();
    run.fConverged = analysis.Run(steps, monitor);
    run.fWall = std::chrono::duration<REAL>(std::chrono::steady_clock::now() - start).count();
    run.fMeanEval = mcc::MeanEvaluations(analysis.StepLog());
    run.fGlobalIts = analysis.NGlobalIterations();
    run.fBisections = analysis.NBisections();
    if (!run.fConverged) std::cerr << "KLTriaxial: the finite element solution of " << c.fTag << " with " << run.fN
                                   << " increments stopped before the end\n";
    if (run.fN == fNFigure) {
        mcc::WriteGaussPointsVTK(mat, mphys, "kl_" + c.fTag + "_fe_gauss.vtk");
        mcc::WriteNodalVTK(analysis, 3, "kl_" + c.fTag + "_fe.vtk", 0);
    }
    mcc::DeleteMeshes(mphys);
}

inline KLTriaxial::TResult KLTriaxial::Run(const TCase &c) const {
    TResult res;
    res.fCase = c;
    ClosedForm(c, res);
    {
        std::vector<std::vector<REAL>> rows;
        for (auto &r : res.fClosed) rows.push_back({r[0], r[1], r[2], r[3]});
        mcc::WriteCSV("kl_" + c.fTag + "_closed.csv", {"eps_a", "p_eff", "q", c.fDrained ? "eps_v" : "u"}, rows);
    }
    for (int n : fIncrements) {
        TRun run;
        run.fN = n;
        run.fPoint = MaterialPoint(c, n, &run.fStats);
        if (run.fPoint.empty()) {
            std::cerr << "KLTriaxial: the material point solution of " << c.fTag << " with " << n << " increments failed\n";
            res.fRuns.push_back(run);
            continue;
        }
        run.fPointErrQ = run.fPoint.back()[2] - res.fQRef;
        if (fFECheck) RunFiniteElement(c, run);
        if (!run.fHistory.empty()) {
            const auto &e = run.fHistory.back();
            run.fErrQ = e[2] - res.fQRef;
            run.fErrP = e[1] - res.fPRef;
            run.fErrEv = e[3] - res.fEvRef;
            run.fQPeak = -1.;
            for (auto &h : run.fHistory)
                if (h[2] > run.fQPeak) {
                    run.fQPeak = h[2];
                    run.fEaPeak = h[0];
                }
            std::vector<std::vector<REAL>> rows;
            for (size_t k = 0; k < run.fHistory.size() && k < run.fPoint.size(); ++k) {
                const auto &f = run.fHistory[k];
                const auto &m = run.fPoint[k];
                run.fMaxDiffPoint = std::max(run.fMaxDiffPoint, std::fabs(f[2] - m[2]));
                rows.push_back({f[0], f[1], f[2], f[3], f[4], f[5], m[1], m[2], m[3], f[2] - m[2]});
            }
            mcc::WriteCSV("kl_" + c.fTag + "_fe_n" + std::to_string(n) + ".csv",
                          {"eps_a", "p_eff", "q", "eps_v", "u", "evaluations", "p_point", "q_point",
                           c.fDrained ? "eps_v_point" : "u_point", "q_minus_q_point"},
                          rows);
        }
        res.fRuns.push_back(run);
    }
    return res;
}

inline void KLTriaxial::RunAll() {
    const IntegrationSchemes::TMaterial P = Material();
    std::cout << "Triaxial tests of example 1 of Krabbenhoft and Lyamin (2012) with one Hex20-Hex8 element (Sect. 7.1)\n"
              << std::setprecision(6) << "M = " << P.fM << ", lambda = " << P.fLambda << ", kappa = " << P.fKappa
              << ", v0 = " << P.fV0 << ", nu = " << P.fNu << ", p'0 = " << fP0 << " kPa\n";
    std::vector<std::vector<std::string>> table;
    for (const TCase &c : Cases()) {
        const TResult r = Run(c);
        std::cout << "\n" << c.fLabel << " (p'_c0 = " << c.fOCR * fP0 << " kPa, eps_a = " << 100. * c.fEaMax
                  << "%): closed form at the end p' = " << std::fixed << std::setprecision(4) << r.fPRef
                  << ", q = " << r.fQRef << " kPa";
        if (c.fDrained)
            std::cout << ", eps_v = " << std::setprecision(6) << r.fEvRef << "; peak q = " << std::setprecision(3)
                      << r.fQPeakRef << " kPa at eps_a = " << 100. * r.fEaPeakRef << "%";
        else
            std::cout << ", u = " << fP0 + r.fQRef / 3. - r.fPRef << " kPa";
        std::cout << std::defaultfloat << "\n";
        std::cout << "  increments  q_FE - q_exact  q_point - q_exact  " << (c.fDrained ? "eps_v - exact " : "p' - p'_exact ")
                  << " max|q_FE - q_point|  spread q  eval/inc (max)  global its  bisections  time (s)  local its mean (max)\n";
        for (auto &run : r.fRuns) {
            std::cout << "  " << std::setw(10) << run.fN << std::scientific << std::setprecision(3) << std::setw(16)
                      << run.fErrQ << std::setw(19) << run.fPointErrQ << std::setw(15)
                      << (c.fDrained ? run.fErrEv : run.fErrP) << std::setw(21) << run.fMaxDiffPoint << std::setw(10)
                      << run.fSpreadQ << std::fixed << std::setprecision(2) << std::setw(10) << run.fMeanEval << " ("
                      << run.fMaxEval << ")" << std::setw(12) << run.fGlobalIts << std::setw(12) << run.fBisections
                      << std::setw(10) << run.fWall << std::setw(12) << run.fStats.Mean() << " (" << run.fStats.fMax
                      << ")" << (run.fConverged ? "" : "  FE STOPPED") << std::defaultfloat << "\n";
            if (!c.fDrained && !run.fHistory.empty())
                std::cout << "              nodal pore pressure at the end " << std::fixed << std::setprecision(4)
                          << run.fUEnd << " kPa, from the stresses " << run.fHistory.back()[4] << " kPa"
                          << std::defaultfloat << "\n";
            if (c.fDrained && !run.fHistory.empty())
                std::cout << "              peak q of the element " << std::fixed << std::setprecision(3) << run.fQPeak
                          << " kPa at eps_a = " << 100. * run.fEaPeak << "% (closed form " << r.fQPeakRef << " at "
                          << 100. * r.fEaPeakRef << "%)" << std::defaultfloat << "\n";
            std::ostringstream o;
            o << std::setprecision(17);
            auto num = [&o](REAL v) {
                o.str("");
                o << v;
                return o.str();
            };
            table.push_back({c.fTag, std::to_string(run.fN), std::to_string(int(run.fConverged)), num(run.fErrQ),
                             num(run.fPointErrQ), num(run.fErrP), num(run.fErrEv), num(run.fMaxDiffPoint),
                             num(run.fSpreadQ), num(run.fMeanEval), std::to_string(run.fMaxEval),
                             std::to_string(run.fGlobalIts), std::to_string(run.fBisections), num(run.fWall),
                             num(run.fStats.Mean()), std::to_string(run.fStats.fMax), num(r.fQRef), num(r.fPRef),
                             num(r.fEvRef), num(run.fQPeak), num(run.fEaPeak), num(r.fQPeakRef), num(r.fEaPeakRef)});
        }
    }
    std::ofstream out("kl_table.csv");
    out << "case,increments,converged,err_q_fe,err_q_point,err_p_fe,err_eps_v_fe,max_diff_q_fe_point,spread_q,"
           "mean_evaluations,max_evaluations,global_iterations,bisections,wall_time_s,local_its_mean,local_its_max,"
           "q_exact,p_exact,eps_v_exact,q_peak_fe,eps_a_peak_fe,q_peak_exact,eps_a_peak_exact\n";
    for (auto &row : table) {
        for (size_t i = 0; i < row.size(); ++i) out << row[i] << (i + 1 < row.size() ? "," : "\n");
    }
    std::cout << "\nFiles: kl_<case>_closed.csv, kl_<case>_fe_n<n>.csv, kl_table.csv and the VTK files of the "
              << fNFigure << "-increment runs\n";
}

inline void KLTriaxial::RunTangents() {
    using P = TPZPlasticStepModifiedCamClay;
    const P::ETangentMode modes[5] = {P::EConsistentTangent, P::EFiniteDifferenceTangent, P::ESymmetricTangent,
                                      P::EContinuumTangent, P::ETransposedTangent};
    std::cout << "\nGlobal iterations of the drained tests with the tangent operators (" << fNTangents
              << " increments, one Hex20-Hex8 element)\n";
    std::ofstream out("kl_tangents.csv");
    out << "case,tangent,increments,converged,mean_evaluations,max_evaluations,global_iterations,bisections,wall_time_s,"
           "q_end,p_end\n";
    const P::ETangentMode keep = fTangent;
    for (const TCase &c : Cases()) {
        if (!c.fDrained) continue;
        REAL q0 = 0.;
        for (int m = 0; m < 5; ++m) {
            fTangent = modes[m];
            TRun run;
            run.fN = fNTangents;
            RunFiniteElement(c, run);
            const auto &e = run.fHistory.back();
            if (m == 0) q0 = e[2];
            std::cout << "  " << std::setw(16) << c.fTag << std::setw(6) << P::TangentModeName(modes[m]) << ": "
                      << std::fixed << std::setprecision(2) << run.fMeanEval << " evaluations per increment (max "
                      << run.fMaxEval << "), global iterations " << run.fGlobalIts << ", bisections "
                      << run.fBisections << ", " << std::setprecision(3) << run.fWall << " s; q_end - q_end(D) = "
                      << std::scientific << std::setprecision(2) << e[2] - q0 << std::defaultfloat
                      << (run.fConverged ? "" : "  STOPPED") << "\n";
            out << std::setprecision(17) << c.fTag << "," << P::TangentModeName(modes[m]) << "," << run.fN << ","
                << int(run.fConverged) << "," << run.fMeanEval << "," << run.fMaxEval << "," << run.fGlobalIts << ","
                << run.fBisections << "," << run.fWall << "," << e[2] << "," << e[1] << "\n";
        }
    }
    fTangent = keep;
}
