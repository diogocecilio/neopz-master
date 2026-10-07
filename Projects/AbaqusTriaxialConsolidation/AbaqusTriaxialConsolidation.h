/**
 * @file AbaqusTriaxialConsolidation.h
 * @brief Sect. 6.5 of the article: consolidation of a triaxial specimen (Abaqus benchmark 1.15.2) with the
 * three-dimensional Hex20-Hex8 model (Figs. 9 to 11 and Table 7), and the global Newton iterations of
 * Sect. 6.7 in the same benchmark (Table 9 and the Abaqus columns of Table 10).
 */
#pragma once

#include "MCCPaperTools.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "pzstepsolver.h"
#include <chrono>
#include <cstdint>
#include <ctime>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <tuple>

/**
 * @ingroup mccpaper
 * @brief Abaqus benchmark 1.15.2: drained triaxial test with displacement control on a cylindrical clay
 * specimen, modelled in three dimensions (functions abaqus3d() of gen_data3d.py, abaqus_mp() and
 * abaqus_states() of gen_data.py, and the Abaqus part of tangentes() of gen_data.py, there with the
 * axisymmetric model).
 *
 * By symmetry a quarter of the upper half of the specimen is modelled (H = 60 mm, R = 20 mm), with the
 * mixed u-p element Hex20-Hex8 (serendipity quadratic displacements, trilinear pore pressures) and the
 * quadratic geometry of mcc::CreateQuarterCylinderMesh: 2 x 2 x 4 divisions, i.e. 48 elements, 321 nodes
 * and 95 pore pressure nodes, or the refined mesh 4 x 4 x 8 (384 elements) of the softening study.
 * Material: M = 1, lambda = 0.174, kappa = 0.026, v0 = 2.08 (e0 = 1.08), porous elasticity with nu = 0.3,
 * p'c0 = 116.6 kPa; mobility k = 1.728e-4 m/day / gamma_w = 2e-10 m^2/(kPa s). The specimen starts from
 * the isotropic effective stress p'0 (100 kPa in the benchmark) in equilibrium with the cell pressure,
 * which is kept on the curved lateral face as the radial traction -p'0 (x, y, 0)/r; the top platen
 * (z = H) moves down at a constant rate up to delta/H = 0.6 in 400 days, with free drainage (p_w = 0)
 * at the top. The mid-plane z = 0 is impermeable with u_z = 0, and the symmetry planes x = 0 and y = 0
 * have u_x = 0 and u_y = 0. With the smooth platen only u_z is prescribed on the top, with the rough
 * platen also u_x = u_y = 0.
 *
 * Monitored quantities: effective stress at the point A (x = 5 mm, y = 0, z = 7.5 mm; r = 5 mm,
 * z = 7.5 mm as in the benchmark), interpolated from the integration points of the element that
 * contains it (mcc::StressAtPoint). In the axisymmetric mesh of the benchmark A is the centroid of an
 * element; in the quarter of cylinder it lies on the edge x = 5 mm, y = 0 of the element next to the axis
 * (the first element that contains it, as in the Python code), so that its stress is extrapolated from the
 * integration points (linearly with 2 x 2 x 2 points, quadratically with 3 x 3 x 3 points); the smallest,
 * mean and largest q at the integration points of that element are also recorded. Further: average axial
 * stress on the platen (reaction of the top divided by the area pi R^2/4) and the global deviatoric stress
 * sigma_a - p'0; largest excess pore pressure; volumetric strain from the displacement of the outer top node
 * (exact only for a homogeneous deformation: smooth platen with p'0 = 100 kPa).
 */
class AbaqusTriaxialConsolidation {
public:
    /** @brief Material and boundary ids: x = 0, lateral face, y = 0, mid-plane z = 0, platen z = H (displacement
     * and pore pressure conditions) */
    enum { EMatId = 1, EX0 = -21, ELateral = -22, EY0 = -23, EZ0 = -25, EZH = -26, EPZH = -36 };

    /** @brief Tangent operator of the global iterations (Sect. 6.7, Table 10) */
    typedef TPZPlasticStepModifiedCamClay::ETangentMode ETangentMode;

    /** @brief Definition of a finite element run */
    struct TConfig {
        bool fRough = false;  ///< rough platen (u_x = u_y = 0 on the top)
        int fIntegration = 3; ///< integration order: 3 -> 2 x 2 x 2 (reduced), 4 -> 3 x 3 x 3 (full) points
        /** @brief operator returned by the stress update (consistent tangent D by default) */
        ETangentMode fTangent = TPZPlasticStepModifiedCamClay::EConsistentTangent;
        int fNSteps = 150;    ///< number of increments
        REAL fP0 = 100.;      ///< initial isotropic effective stress = cell pressure (kPa)
        int fNc = 2;          ///< divisions of the central square and along each outer arc
        int fNr = 2;          ///< radial divisions of the outer blocks
        int fNz = 4;          ///< divisions along the height
        /** @brief VTK file series: every fVTKStride-th state is written (0: no series for this run) */
        int fVTKStride = 1;
        /**
         * @brief Renumber the equations with the bandwidth optimization of NeoPZ (TPZSloanRenumbering) before the u-p
         * analysis is built. By default (false) the pressure equations follow the displacement equations, as
         * TPZPoroElastoPlasticUPAnalysis requires for undrained steps; with true they are interleaved with them,
         * which is safe here because every increment has a positive time step (the pressure block -dt H is
         * negative definite), and the skyline LU of the refined mesh costs 4.4 times less (Softening).
         */
        bool fOptimizeBandwidth = false;
        /** @brief tolerance of the normalized residual of the global iterations (Sect. 5.4; 1e-8 in the article) */
        REAL fTolerance = 1.e-8;
        std::string fName = "smooth_2x2x2"; ///< name of the run (files abaqus_\<name\>*)
    };

    /** @brief Columns of the history of a run (file abaqus_\<name\>.csv) */
    enum EHistory {
        EDH = 0, EPA, EQA, ESigmaA, EQPlaten, EMaxPw, EEpsV, EEvaluations, EQElMin, EQElMean, EQElMax, ENHistory
    };

    /** @brief Results of a run */
    struct TResult {
        /** @brief rows (delta/H, p' at A, q at A, axial stress on the platen, sigma_a - p'0, max |p_w|, eps_v,
         * evaluations of the residual in the increment, smallest, mean and largest q at the integration points
         * of the element that contains A), the initial state first */
        std::vector<std::array<REAL, ENHistory>> fHistory;
        std::vector<mcc::TAnalysis::TStepLog> fLog; ///< convergence records of the converged (sub)increments
        REAL fMeanEvaluations = 0.;                 ///< mean evaluations per converged (sub)increment
        int fMaxEvaluations = 0;                    ///< largest number of evaluations of a (sub)increment
        int64_t fNGlobalIterations = 0;             ///< all global iterations, failed attempts included
        int64_t fNBisections = 0;                   ///< bisected (sub)increments
        /**
         * @brief wall time of the solution (s): meshes and increments (assembly, stress updates, linear
         * solutions), without the monitored quantities (fMonitorTime) and the VTK series (fVTKTime)
         */
        REAL fWallTime = 0.;
        /** @brief wall time of the monitored quantities (s): stress at A, reaction of the platen, nodal values */
        REAL fMonitorTime = 0.;
        REAL fVTKTime = 0.;                         ///< wall time of the VTK file series (s)
        /**
         * @brief processor time of the solution (s, std::clock), without the monitored quantities and the VTK
         * series: the solution runs on one thread, so that this is the wall time it would take on an idle machine
         * (the wall time grows with the load of the machine)
         */
        REAL fCPUTime = 0.;
        int fNElements = 0, fNNodes = 0, fNPressureNodes = 0;
        /** @brief degrees of freedom of the multiphysics mesh (TPZCompMesh::NEquations, Dirichlet ones included) */
        int64_t fNEquations = 0;
        int64_t fNFreeEquations = 0;                ///< equations solved (without the Dirichlet conditions)
        bool fCompleted = true;                     ///< false if the analysis stopped before the end
    };

    REAL fR = 0.020, fH = 0.060, fM = 1., fLambda = 0.174, fKappa = 0.026, fV0 = 2.08, fNu = 0.3, fPc0 = 116.6;
    REAL fPerm = 1.728e-4 / 86400. / 10.; ///< mobility k (m^2/(kPa s))
    REAL fTend = 34.56e6;                 ///< 400 days
    REAL fDHend = 0.6;                    ///< final delta/H
    REAL fXA = 0.005, fZA = 0.0075;       ///< point A (x = r, y = 0, z)
    REAL fToleranceSensitivity = 1.e-9;   ///< tolerance of the part "tolerance" (ToleranceSensitivity)

    /**
     * @brief Write the VTK file series of the finite element runs (mcc::TVTKSeries) in vtk/\<run name\>, with
     * the series time delta/H (the command line argument "novtk" disables it)
     */
    bool fWriteVTK = true;

    /** @brief The plastic model of the benchmark with the given tangent operator */
    mcc::TPlastic Model(ETangentMode tangent = TPZPlasticStepModifiedCamClay::EConsistentTangent) const;

    /** @brief Geometric mesh: quarter of cylinder with fNc x fNr x fNz divisions and its boundary faces */
    TPZGeoMesh *CreateGeoMesh(const TConfig &cfg);

    /** @brief Computational meshes, material, boundary conditions and initial state */
    TPZMultiphysicsCompMesh *CreateCompMesh(TPZGeoMesh *gmesh, const TConfig &cfg, mcc::TPoroMaterial *&mat);

    /**
     * @brief Solves a configuration and writes its post-processing files: abaqus_\<name\>.csv (history),
     * abaqus_\<name\>_profile.csv (final displacement of the outer generatrix y = 0, r = R), the VTK files of the
     * final state and, with fWriteVTK and cfg.fVTKStride > 0, the VTK file series in vtk/\<name\>
     * (abaqus_\<name\>_nodal.vtk.series, abaqus_\<name\>_intpoints.vtk.series, abaqus_\<name\>_gausspoints.vtk.series
     * and abaqus_\<name\>_states.csv with delta/H, the time and the platen displacement of each state)
     */
    TResult Run(const TConfig &cfg);

    /** @brief Fig. 9: CSV files of the geometric meshes (coarse 2 x 2 x 4 and refined 4 x 4 x 8, mcc::WriteMeshCSV) */
    void WriteMeshes();

    /** @brief Fig. 10: drained tests at a material point from p'0 = 100 and 20 kPa (600 increments) and closed form */
    void MaterialPointStates();

    /** @brief Table 7 (column delta delta/H = 0.02): material point solution with 30, 150 and 3000 increments */
    void MaterialPointIncrements();

    /**
     * @brief Fig. 11 and Table 7: smooth platen with 2 x 2 x 2 points, rough platen with 2 x 2 x 2 (reduced) and
     * 3 x 3 x 3 (full) points, 150 increments; comparison with the digitized Abaqus curves
     */
    void FiniteElementModels();

    /**
     * @brief Sect. 6.7, Tables 9 and 10: rough platen with 2 x 2 x 2 points and the operators D, central
     * differences, (D + D^T)/2, continuum tangent and D^T (no VTK series; the times of Table 10 are the processor
     * times of the solutions, TResult::fCPUTime, which do not depend on the load of the machine)
     */
    void TangentOperators();

    /**
     * @brief Sect. 6.7, sensitivity of Table 10 to the tolerance of the global iterations: the five operators of
     * TangentOperators with the tolerance fToleranceSensitivity (abaqus_table10_tolerance.csv). Since
     * \f$\|f_{ext}\| < 1\f$ kN, the normalized residual \f$\|R\|/\max(\|f_{ext}\|,1)\f$ is the norm of the residual
     * in kN, and the residuals of the quarter model are about ten times smaller than those of the axisymmetric
     * model of the previous version of the article (whole ring, 37 nodes; 1.30e-2 against 1.30e-1 in the
     * predictor of the first increment): the tolerance 1e-9 here is the 1e-8 of that model.
     */
    void ToleranceSensitivity();

    /** @brief Sect. 6.5: smooth platen with 600 increments from p'0 = 100 and 20 kPa (mesh 2 x 2 x 4) */
    void States();

    /**
     * @brief Sect. 6.5: softening state p'0 = 20 kPa with the refined mesh 4 x 4 x 8 (smooth platen, 600 increments,
     * 6576 degrees of freedom, 5727 free; the equations are renumbered, TConfig::fOptimizeBandwidth, so that the
     * run takes about 25 min instead of about 2.3 h with the default numbering)
     */
    void Softening();

    /**
     * @brief Runs the parts given by name ("mesh", "mp", "fe", "tangents", "tolerance", "states", "softening"; "all"
     * runs everything)
     */
    void RunAll(const std::vector<std::string> &parts);

    /** @brief Digitized Abaqus curve q at A against delta/H (reference/abaqus_1_15_2_digitalizado.json, "qd") */
    static const std::vector<std::array<REAL, 2>> &AbaqusCurve(bool rough);

protected:
    /** @brief Geometric node closest to x */
    static int64_t ClosestNode(TPZGeoMesh *gmesh, const TPZVec<REAL> &x);
    /** @brief Value of the history column c at delta/H = d (linear interpolation) */
    static REAL At(const TResult &r, REAL d, int c) { return mcc::Interpolate(r.fHistory, d, c); }
    /** @brief Largest value of the history column c and its delta/H */
    static std::pair<REAL, REAL> Peak(const TResult &r, int c);
    /**
     * @brief Root mean square and largest difference of the column c of a history (abscissa in column 0)
     * with respect to a digitized curve, and the abscissa of the largest difference
     */
    template <size_t N>
    static std::array<REAL, 3> CompareWithCurve(const std::vector<std::array<REAL, N>> &h, int c,
                                                const std::vector<std::array<REAL, 2>> &ref);
    /** @brief Writes one row per run with the counters and the final values (summary files of the parts) */
    void WriteSummary(const std::string &file, const std::vector<std::pair<TConfig, TResult>> &runs) const;
    /**
     * @brief Solves smooth-platen configurations and compares them with the material point (same increments):
     * files abaqus_\<file\>_summary.csv and abaqus_\<file\>_numbers.csv
     */
    void SmoothRuns(const std::vector<TConfig> &configs, const std::string &file);
    /** @brief Writes a list of named numbers (quantity, value, description) */
    static void WriteNumbers(const std::string &file,
                             const std::vector<std::tuple<std::string, REAL, std::string>> &numbers);
};

inline mcc::TPlastic AbaqusTriaxialConsolidation::Model(ETangentMode tangent) const {
    mcc::TPlastic model;
    model.SetModifiedCamClay(fM, fLambda, fKappa);
    model.SetPoissonRatio(fNu);
    model.SetPorousElasticity();
    model.SetDefaultSpecificVolume(fV0);
    model.SetTangentMode(tangent);
    return model;
}

inline const std::vector<std::array<REAL, 2>> &AbaqusTriaxialConsolidation::AbaqusCurve(bool rough) {
    // Abaqus Benchmarks Manual 2016, 1.15.2, Figs. 1.15.2-2 and -3, point A (digitized; uncertainty 1-2 kPa)
    static const std::vector<std::array<REAL, 2>> smooth = {
        {0.0297, 59.68},  {0.0791, 89.96},  {0.1301, 109.94}, {0.1695, 122.81}, {0.2189, 131.09}, {0.27, 137.23},
        {0.3192, 141.09}, {0.3687, 144.17}, {0.4081, 145.85}, {0.4593, 147.25}, {0.5086, 148.08}};
    static const std::vector<std::array<REAL, 2>> rougher = {
        {0.0304, 59.69},  {0.08, 91.89},    {0.1294, 111.84}, {0.1688, 125.81}, {0.2201, 134.73}, {0.2697, 140.85},
        {0.3191, 144.98}, {0.3702, 147.78}, {0.4097, 150.08}, {0.4593, 152.19}, {0.5096, 153.09}};
    return rough ? rougher : smooth;
}

inline TPZGeoMesh *AbaqusTriaxialConsolidation::CreateGeoMesh(const TConfig &cfg) {
    const REAL R = fR, H = fH;
    auto marker = [R, H](const std::array<std::array<REAL, 3>, 4> &X) {
        const REAL tol = 1e-9 * R;
        if (mcc::FaceOnPlane(X, 0, 0., tol)) return std::vector<int>{EX0};
        if (mcc::FaceOnPlane(X, 1, 0., tol)) return std::vector<int>{EY0};
        if (mcc::FaceOnPlane(X, 2, 0., tol)) return std::vector<int>{EZ0};
        if (mcc::FaceOnPlane(X, 2, H, tol)) return std::vector<int>{EZH, EPZH};
        return std::vector<int>{ELateral};
    };
    return mcc::CreateQuarterCylinderMesh(fR, fH, cfg.fNc, cfg.fNr, cfg.fNz, EMatId, marker, ELateral);
}

inline TPZMultiphysicsCompMesh *AbaqusTriaxialConsolidation::CreateCompMesh(TPZGeoMesh *gmesh, const TConfig &cfg,
                                                                           mcc::TPoroMaterial *&mat) {
    const int dim = 3;
    const std::set<int> bcids = {EX0, ELateral, EY0, EZ0, EZH, EPZH};
    TPZCompMesh *cmeshU = mcc::CreateDisplacementMesh(gmesh, dim, EMatId, bcids);
    TPZCompMesh *cmeshP = mcc::CreatePressureMesh(gmesh, dim, EMatId, bcids);

    TPZMultiphysicsCompMesh *mphys = new TPZMultiphysicsCompMesh(gmesh);
    mphys->SetDimModel(dim);
    mat = new mcc::TPoroMaterial(EMatId, TPZMatPoroElastoPlasticUPBase::EThreeDimensional);
    mat->SetPlasticModel(Model(cfg.fTangent));
    mat->SetBiot(1., 0.);
    mat->SetPermeability(fPerm);
    mat->SetIntegrationOrder(cfg.fIntegration);
    mphys->InsertMaterialObject(mat);

    using B = TPZMatPoroElastoPlasticUPBase;
    TPZFNMatrix<9, STATE> val1(dim, dim, 0.);
    TPZManVector<STATE, 3> val2(dim, 0.), zero(1, 0.);
    auto directional = [&](int id, std::initializer_list<int> comps) {
        val1.Zero();
        for (int c : comps) val1(c, c) = 1.;
        mphys->InsertMaterialObject(mat->CreateBC(mat, id, B::EDirichletUDirectional, val1, val2));
    };
    const REAL p0 = cfg.fP0;
    directional(EX0, {0});                      // symmetry plane x = 0
    directional(EY0, {1});                      // symmetry plane y = 0
    directional(EZ0, {2});                      // mid-plane
    if (cfg.fRough) directional(EZH, {0, 1, 2}); // rough platen: u_x = u_y = 0, u_z controlled
    else directional(EZH, {2});                 // smooth platen: u_z controlled
    val1.Zero();
    auto *lateral = mat->CreateBC(mat, ELateral, B::ENeumannU, val1, val2);
    // radial cell pressure on the curved face: t = -p0 (x, y, 0)/r
    lateral->SetForcingFunctionBC([p0](const TPZVec<REAL> &x, TPZVec<STATE> &val, TPZFMatrix<STATE> &) {
        const REAL r = std::hypot(x[0], x[1]);
        val[0] = -p0 * x[0] / r;
        val[1] = -p0 * x[1] / r;
        val[2] = 0.;
    });
    mphys->InsertMaterialObject(lateral);
    mphys->InsertMaterialObject(mat->CreateBC(mat, EPZH, B::EDirichletP, val1, zero)); // drained platen
    mcc::BuildMultiphysics(mphys, cmeshU, cmeshP, EMatId, bcids);
    const REAL pc0 = fPc0, v0 = fV0;
    mat->InitializeMemory(mphys, [p0, pc0, v0](const TPZVec<REAL> &, TPZElastoPlasticMem &mem) {
        mem.m_sigma = mcc::IsotropicTensor(-p0);
        mem.m_elastoplastic_state.m_hardening = pc0;
        mem.m_elastoplastic_state.fmatprop.Resize(1, v0);
        mem.m_elastoplastic_state.fmatprop[0] = v0;
        mem.m_elastoplastic_state.fpressure = 0.;
    });
    return mphys;
}

inline int64_t AbaqusTriaxialConsolidation::ClosestNode(TPZGeoMesh *gmesh, const TPZVec<REAL> &x) {
    int64_t best = -1;
    REAL dbest = 1e300;
    for (int64_t n = 0; n < gmesh->NNodes(); ++n) {
        TPZManVector<REAL, 3> co(3);
        gmesh->NodeVec()[n].GetCoordinates(co);
        REAL d = 0.;
        for (int k = 0; k < x.size(); ++k) d += (co[k] - x[k]) * (co[k] - x[k]);
        if (d < dbest) {
            dbest = d;
            best = n;
        }
    }
    return best;
}

inline std::pair<REAL, REAL> AbaqusTriaxialConsolidation::Peak(const TResult &r, int c) {
    std::pair<REAL, REAL> pk(-1e300, 0.);
    for (auto &h : r.fHistory)
        if (h[c] > pk.first) pk = {h[c], h[EDH]};
    return pk;
}

template <size_t N>
inline std::array<REAL, 3> AbaqusTriaxialConsolidation::CompareWithCurve(const std::vector<std::array<REAL, N>> &h,
                                                                         int c,
                                                                         const std::vector<std::array<REAL, 2>> &ref) {
    REAL sum = 0., dmax = 0., at = 0.;
    for (auto &pt : ref) {
        const REAL d = mcc::Interpolate(h, pt[0], c) - pt[1];
        sum += d * d;
        if (std::fabs(d) > dmax) {
            dmax = std::fabs(d);
            at = pt[0];
        }
    }
    return {std::sqrt(sum / ref.size()), dmax, at};
}

inline AbaqusTriaxialConsolidation::TResult AbaqusTriaxialConsolidation::Run(const TConfig &cfg) {
    const auto start = std::chrono::steady_clock::now();
    const std::clock_t cstart = std::clock();
    REAL cpuExcluded = 0.; // processor time of the monitored quantities and of the VTK series
    TPZGeoMesh *gmesh = CreateGeoMesh(cfg);
    mcc::TPoroMaterial *mat = nullptr;
    TPZMultiphysicsCompMesh *mphys = CreateCompMesh(gmesh, cfg, mat);

    if (cfg.fOptimizeBandwidth) {
        TPZLinearAnalysis renumber(mphys, true); // permutes the sequence numbers of the connects of mphys
    }
    mcc::TAnalysis analysis(mphys, mat);
    TPZSkylineNSymStructMatrix<STATE> skyl(mphys);
    skyl.SetNumThreads(0);
    analysis.SetStructuralMatrix(skyl);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELU);
    analysis.SetSolver(step);
    analysis.SetControlledDisplacement(EZH, 2);
    analysis.SetPredictor(true);
    analysis.SetNewtonParameters(cfg.fTolerance, 25, 8);
    analysis.ResetCounters();

    TResult res;
    res.fNNodes = gmesh->NNodes();
    res.fNPressureNodes = analysis.NodesOfMaterials({EMatId}).size();
    for (int64_t i = 0; i < gmesh->NElements(); ++i)
        if (gmesh->Element(i) && gmesh->Element(i)->MaterialId() == EMatId) res.fNElements++;
    res.fNEquations = mphys->NEquations();
    analysis.IdentifyEquations();
    for (bool c : analysis.ConstrainedEquations()) res.fNFreeEquations += !c;
    TPZManVector<REAL, 3> corner(3, 0.), pointA(3, 0.);
    corner[0] = fR;
    corner[2] = fH;
    const int64_t iR = ClosestNode(gmesh, corner); // outer top node on y = 0 (u_x = u_r)
    pointA[0] = fXA;
    pointA[2] = fZA;
    // element that contains A (the first one, as in mcc::StressAtPoint): A lies on its edge x = 5 mm, y = 0
    TPZGeoEl *gelA = nullptr;
    {
        TPZManVector<REAL, 3> qsi(3, 0.);
        gelA = mcc::LocatePoint(gmesh, EMatId, pointA, qsi);
        if (!gelA) DebugStop();
    }
    const REAL area = M_PI * fR * fR / 4.;
    std::vector<mcc::TAnalysis::TLoadState> steps;
    for (int k = 1; k <= cfg.fNSteps; ++k)
        steps.emplace_back(fTend * k / cfg.fNSteps, 1., -fDHend * fH * k / cfg.fNSteps);
    const std::string prefix = "abaqus_" + cfg.fName;
    // VTK file series (series time delta/H)
    std::unique_ptr<mcc::TVTKSeries> vtk;
    if (fWriteVTK && cfg.fVTKStride > 0)
        vtk = std::make_unique<mcc::TVTKSeries>(mphys, mat, "vtk/" + cfg.fName, prefix, "delta_H",
                                                std::vector<std::string>{"t", "uc"});
    size_t nlog = 0;
    const std::set<int64_t> pnodes = analysis.NodesOfMaterials({EMatId});
    auto monitor = [&](int istep, const mcc::TAnalysis::TLoadState &s) {
        const auto tm = std::chrono::steady_clock::now();
        const std::clock_t cm = std::clock();
        const TPZTensor<REAL> sig = mcc::StressAtPoint(mat, mphys, pointA);
        const REAL dH = fDHend * s.fTime / fTend;
        const REAL sa = (s.fTime > 0.) ? -analysis.Reaction({EZH}, 2) / area : cfg.fP0;
        REAL pmax = 0.;
        for (auto n : pnodes) pmax = std::max(pmax, std::fabs(analysis.NodalValue(n, 1, 0)));
        // exact only for a homogeneous deformation (smooth platen with hardening, p'0 = 100 kPa)
        const REAL ev = dH - 2. * analysis.NodalValue(iR, 0, 0) / fR;
        // evaluations of the residual in this increment (its converged sub-increments)
        const auto &log = analysis.StepLog();
        REAL nev = 0.;
        for (size_t i = nlog; i < log.size(); ++i) nev += log[i].fResiduals.size();
        nlog = log.size();
        // q at the integration points of the element that contains A: the spread shows the oscillation of the
        // stresses within the element (volumetric locking of the full integration)
        REAL qmin = 1e300, qmax = -1e300, qsum = 0.;
        int nq = 0;
        for (auto &g : mcc::GaussPoints(mat, mphys)) {
            if (mphys->Element(g.fElement)->Reference() != gelA) continue;
            const REAL q = mcc::DeviatoricStress(g.fSigma);
            qmin = std::min(qmin, q);
            qmax = std::max(qmax, q);
            qsum += q;
            nq++;
        }
        res.fHistory.push_back({dH, mcc::MeanEffectiveStress(sig), mcc::DeviatoricStress(sig), sa, sa - cfg.fP0, pmax,
                                ev, nev, qmin, qsum / std::max(nq, 1), qmax});
        res.fMonitorTime += std::chrono::duration<REAL>(std::chrono::steady_clock::now() - tm).count();
        if (vtk && (istep % cfg.fVTKStride == 0 || istep == cfg.fNSteps)) {
            const auto t0 = std::chrono::steady_clock::now();
            vtk->Write(dH, {s.fTime, s.fUc});
            res.fVTKTime += std::chrono::duration<REAL>(std::chrono::steady_clock::now() - t0).count();
        }
        cpuExcluded += REAL(std::clock() - cm) / CLOCKS_PER_SEC;
    };
    res.fCompleted = analysis.Run(steps, monitor);
    res.fWallTime = std::chrono::duration<REAL>(std::chrono::steady_clock::now() - start).count() - res.fVTKTime -
                    res.fMonitorTime;
    res.fCPUTime = REAL(std::clock() - cstart) / CLOCKS_PER_SEC - cpuExcluded;
    if (!res.fCompleted) std::cout << "  " << cfg.fName << ": the analysis stopped before the end\n";
    res.fLog = analysis.StepLog();
    res.fMeanEvaluations = mcc::MeanEvaluations(res.fLog);
    for (auto &l : res.fLog) res.fMaxEvaluations = std::max(res.fMaxEvaluations, int(l.fResiduals.size()));
    res.fNGlobalIterations = analysis.NGlobalIterations();
    res.fNBisections = analysis.NBisections();

    std::vector<std::vector<REAL>> rows;
    for (auto &h : res.fHistory) rows.push_back(std::vector<REAL>(h.begin(), h.end()));
    mcc::WriteCSV(prefix + ".csv",
                  {"delta_H", "p_A", "q_A", "sigma_a_platen", "q_platen", "max_pw", "eps_v", "evaluations", "q_elA_min",
                   "q_elA_mean", "q_elA_max"},
                  rows);
    // final displacement of the outer generatrix (y = 0, r = R), four points per element height: the bulging of
    // the specimen (displacement interpolated by the element that contains each point)
    {
        mphys->LoadReferences();
        const int varu = mat->VariableIndex("Displacement");
        std::vector<std::vector<REAL>> prof;
        for (int i = 0; i <= 4 * cfg.fNz; ++i) {
            TPZManVector<REAL, 3> x(3, 0.), qsi(3, 0.);
            x[0] = fR;
            x[2] = fH * i / (4. * cfg.fNz);
            TPZGeoEl *gel = mcc::LocatePoint(gmesh, EMatId, x, qsi);
            if (!gel || !gel->Reference()) continue;
            TPZManVector<STATE, 3> u(3, 0.);
            gel->Reference()->Solution(qsi, varu, u);
            prof.push_back({x[2], u[0], u[2]});
        }
        mcc::WriteCSV(prefix + "_profile.csv", {"z", "u_r", "u_z"}, prof);
    }
    mcc::WriteNodalVTK(analysis, 3, prefix + ".vtk", 0);
    mcc::WriteGaussPointsVTK(mat, mphys, prefix + "_gauss.vtk");
    vtk.reset(); // the post-processing meshes refer to the meshes of the run
    mcc::DeleteMeshes(mphys);
    return res;
}

inline void AbaqusTriaxialConsolidation::WriteSummary(const std::string &file,
                                                      const std::vector<std::pair<TConfig, TResult>> &runs) const {
    std::ofstream out(file);
    out << "run,platen,points,tangent,tolerance,nsteps,p0,nc,nr,nz,elements,nodes,pressure_nodes,equations,"
           "free_equations,completed,mean_evaluations,max_evaluations,global_iterations,bisections,wall_time_s,"
           "cpu_time_s,monitor_time_s,vtk_time_s,delta_H_end,p_A_end,q_A_end,sigma_a_end,q_platen_end,max_pw_end,"
           "eps_v_end,q_A_max,"
           "delta_H_q_A_max,q_platen_max,delta_H_q_platen_max,q_elA_min_end,q_elA_mean_end,q_elA_max_end\n";
    out << std::setprecision(12);
    for (auto &rc : runs) {
        const TConfig &c = rc.first;
        const TResult &r = rc.second;
        const auto &e = r.fHistory.back();
        const auto pa = Peak(r, EQA), pp = Peak(r, EQPlaten);
        const int np = c.fIntegration == 3 ? 8 : 27;
        out << c.fName << "," << (c.fRough ? "rough" : "smooth") << "," << np << ","
            << TPZPlasticStepModifiedCamClay::TangentModeName(c.fTangent) << "," << c.fTolerance << "," << c.fNSteps
            << "," << c.fP0 << "," << c.fNc << "," << c.fNr << "," << c.fNz << "," << r.fNElements << ","
            << r.fNNodes << "," << r.fNPressureNodes << "," << r.fNEquations << "," << r.fNFreeEquations << ","
            << int(r.fCompleted) << "," << r.fMeanEvaluations << "," << r.fMaxEvaluations << ","
            << r.fNGlobalIterations << "," << r.fNBisections << "," << r.fWallTime << "," << r.fCPUTime << ","
            << r.fMonitorTime << "," << r.fVTKTime << "," << e[EDH] << "," << e[EPA] << "," << e[EQA] << ","
            << e[ESigmaA] << "," << e[EQPlaten] << "," << e[EMaxPw] << "," << e[EEpsV] << "," << pa.first << ","
            << pa.second << "," << pp.first << "," << pp.second << "," << e[EQElMin] << "," << e[EQElMean] << ","
            << e[EQElMax] << "\n";
    }
}

inline void AbaqusTriaxialConsolidation::WriteNumbers(
    const std::string &file, const std::vector<std::tuple<std::string, REAL, std::string>> &numbers) {
    std::ofstream out(file);
    out << "quantity,value,description\n" << std::setprecision(12);
    for (auto &n : numbers) out << std::get<0>(n) << "," << std::get<1>(n) << ",\"" << std::get<2>(n) << "\"\n";
}

inline void AbaqusTriaxialConsolidation::WriteMeshes() {
    std::cout << "\nFig. 9: geometric meshes of the quarter of the specimen (Hex20 with quadratic geometry)\n";
    for (auto n : {std::array<int, 3>{2, 2, 4}, std::array<int, 3>{4, 4, 8}}) {
        TConfig cfg;
        cfg.fNc = n[0];
        cfg.fNr = n[1];
        cfg.fNz = n[2];
        TPZGeoMesh *gmesh = CreateGeoMesh(cfg);
        const std::string prefix =
            "abaqus_mesh_" + std::to_string(n[0]) + "x" + std::to_string(n[1]) + "x" + std::to_string(n[2]);
        mcc::WriteMeshCSV(gmesh, prefix);
        int nel = 0;
        for (int64_t i = 0; i < gmesh->NElements(); ++i)
            if (gmesh->Element(i) && gmesh->Element(i)->MaterialId() == EMatId) nel++;
        std::cout << "  " << n[0] << " x " << n[1] << " x " << n[2] << ": " << nel << " elements, " << gmesh->NNodes()
                  << " nodes (" << prefix << "_*.csv)\n";
        delete gmesh;
    }
}

inline void AbaqusTriaxialConsolidation::MaterialPointStates() {
    std::cout << "\nFig. 10: drained test at a material point, p'c0 = 116.6 kPa (600 increments) and closed form\n";
    for (REAL p0 : {100., 20.}) {
        mcc::TLocalStats stats;
        auto mp = mcc::TriaxialDrained(Model(), p0, fPc0, fV0, fDHend, 600, &stats);
        auto cf = mcc::TriaxialDrainedClosed(p0, fPc0, fV0, fM, fLambda, fKappa, 0., fNu, 600);
        REAL qpk = 0., eapk = 0., dq = 0., dev = 0.;
        for (auto &r : mp) {
            if (r[2] > qpk) { qpk = r[2]; eapk = r[0]; }
            dq = std::max(dq, std::fabs(r[2] - mcc::Interpolate(cf, r[0], 2)));
            dev = std::max(dev, std::fabs(r[3] - mcc::Interpolate(cf, r[0], 3)));
        }
        REAL qcf = 0.;
        for (auto &r : cf) qcf = std::max(qcf, r[2]);
        const auto &e = mp.back();
        std::cout << "  p'0 = " << p0 << " (R = " << fPc0 / p0 << "): peak q = " << qpk << " at eps_1 = " << eapk
                  << " (closed form " << qcf << "); end p' = " << e[1] << " q = " << e[2] << " eps_v = " << e[3]
                  << "; max |q - closed| = " << dq << " kPa, max |eps_v - closed| = " << dev << "\n";
        std::vector<std::vector<REAL>> rows, rowsc;
        for (auto &r : mp) rows.push_back({r[0], r[1], r[2], r[3]});
        for (auto &r : cf) rowsc.push_back({r[0], r[1], r[2], r[3]});
        const std::string tag = std::to_string(int(p0));
        mcc::WriteCSV("abaqus_material_point_p0_" + tag + ".csv", {"eps_1", "p", "q", "eps_v"}, rows);
        mcc::WriteCSV("abaqus_closed_form_p0_" + tag + ".csv", {"eps_1", "p", "q", "eps_v"}, rowsc);
    }
    std::cout << "  (article: p'0 = 100 kPa reaches the critical state p' = q = 150 kPa with eps_v = 7.2% at 60%;\n"
                 "   p'0 = 20 kPa peaks at q = 54.7 kPa (eps_1 = 2.0%) and softens to q = 30 kPa, eps_v = -4.2%;\n"
                 "   agreement with the closed form within 1 kPa in q and 2e-4 in eps_v)\n";
}

inline void AbaqusTriaxialConsolidation::MaterialPointIncrements() {
    std::cout << "\nMaterial point with the state of the benchmark (smooth platen = homogeneous problem):\n";
    const REAL ref[3] = {149.262, 149.503, 149.553};
    int i = 0;
    for (int n : {30, 150, 3000}) {
        auto mp = mcc::TriaxialDrained(Model(), 100., fPc0, fV0, fDHend, n);
        std::cout << "  " << n << " increments: q(0.6) = " << mp.back()[2] << " kPa (Python " << ref[i++] << ")\n";
        std::vector<std::vector<REAL>> rows;
        for (auto &r : mp) rows.push_back({r[0], r[1], r[2], r[3]});
        mcc::WriteCSV("abaqus_material_point_" + std::to_string(n) + ".csv", {"delta_H", "p", "q", "eps_v"}, rows);
    }
}

inline void AbaqusTriaxialConsolidation::FiniteElementModels() {
    std::cout << "\nQuarter of the specimen, 48 Hex20-Hex8 elements, 150 increments: q at A (kPa), Table 7\n";
    TConfig smooth, rough, full;
    smooth.fName = "smooth_2x2x2";
    rough.fRough = true;
    rough.fName = "rough_2x2x2";
    full.fRough = true;
    full.fIntegration = 4;
    full.fName = "rough_3x3x3";
    std::vector<std::pair<TConfig, TResult>> runs;
    for (auto *c : {&smooth, &rough, &full}) {
        runs.emplace_back(*c, Run(*c));
        const TResult &r = runs.back().second;
        std::cout << "  " << c->fName << ": " << r.fNElements << " elements, " << r.fNNodes << " nodes, "
                  << r.fNPressureNodes << " pressure nodes, " << r.fNEquations << " degrees of freedom ("
                  << r.fNFreeEquations << " free); end q_A = "
                  << r.fHistory.back()[EQA] << " kPa; " << r.fMeanEvaluations << " evaluations per increment (max "
                  << r.fMaxEvaluations << "), " << r.fNBisections << " bisections, " << r.fWallTime << " s (processor "
                  << r.fCPUTime << " s)\n";
    }
    WriteSummary("abaqus_fe_summary.csv", runs);
    const TResult &rs = runs[0].second, &rr = runs[1].second, &rf = runs[2].second;
    auto mp30 = mcc::TriaxialDrained(Model(), 100., fPc0, fV0, fDHend, 30);
    auto mp150 = mcc::TriaxialDrained(Model(), 100., fPc0, fV0, fDHend, 150);

    // Table 7: q at A at the abscissas of the digitized Abaqus curves (every other digitized point) and at 0.6
    const auto &as = AbaqusCurve(false), &ar = AbaqusCurve(true);
    const REAL nominal[7] = {0.03, 0.13, 0.22, 0.32, 0.41, 0.51, 0.60};
    const REAL v06[7][6] = {{59.7, 60.2, 56.7, 59.7, 61.6, 61.6},       {109.9, 114.0, 110.1, 111.8, 116.8, 117.4},
                            {131.1, 133.8, 131.0, 134.7, 137.4, 138.7}, {141.1, 143.5, 141.9, 145.0, 147.2, 147.4},
                            {145.8, 147.1, 146.2, 150.1, 151.6, 149.0}, {148.1, 148.8, 148.4, 153.1, 154.3, 147.7},
                            {-1, 149.5, 149.3, -1, 155.6, 145.5}};
    const REAL nan = std::numeric_limits<REAL>::quiet_NaN();
    std::vector<std::vector<REAL>> t7;
    std::cout << "  d/H  | Abaqus smooth | smooth 2x2x2 | mat. point 0.02 | Abaqus rough | rough 2x2x2 | rough 3x3x3"
                 "   (v0.6 axisymmetric values in parentheses)\n" << std::fixed;
    for (int i = 0; i < 7; ++i) {
        const REAL ds = i < 6 ? as[2 * i][0] : fDHend, dr = i < 6 ? ar[2 * i][0] : fDHend;
        const REAL abs_ = i < 6 ? as[2 * i][1] : nan, abr = i < 6 ? ar[2 * i][1] : nan;
        const REAL vals[6] = {abs_, At(rs, ds, EQA), mcc::Interpolate(mp30, ds, 2), abr, At(rr, dr, EQA),
                              At(rf, dr, EQA)};
        t7.push_back({nominal[i], ds, dr, vals[0], vals[1], vals[2], vals[3], vals[4], vals[5]});
        std::cout << "  " << std::setprecision(2) << nominal[i] << std::setprecision(1);
        for (int k = 0; k < 6; ++k) {
            if (std::isnan(vals[k])) std::cout << " | " << std::setw(13) << "-";
            else std::cout << " | " << std::setw(6) << vals[k] << " (" << std::setw(5) << v06[i][k] << ")";
        }
        std::cout << "\n";
    }
    std::cout << std::defaultfloat << std::setprecision(6);
    mcc::WriteCSV("abaqus_table7.csv",
                  {"delta_H", "delta_H_smooth", "delta_H_rough", "abaqus_smooth", "smooth_2x2x2", "material_point_30",
                   "abaqus_rough", "rough_2x2x2", "rough_3x3x3"},
                  t7);

    // comparison with the digitized Abaqus curves and between the integration rules
    const auto cs = CompareWithCurve(rs.fHistory, EQA, as), cm30 = CompareWithCurve(mp30, 2, as),
               cm150 = CompareWithCurve(mp150, 2, as), cr = CompareWithCurve(rr.fHistory, EQA, ar),
               cf = CompareWithCurve(rf.fHistory, EQA, ar);
    REAL dsm = 0., dsm_at = 0., dsm01 = 0.;
    for (size_t k = 0; k < rs.fHistory.size() && k < mp150.size(); ++k) {
        const REAL d = std::fabs(rs.fHistory[k][EQA] - mp150[k][2]);
        if (d > dsm) { dsm = d; dsm_at = rs.fHistory[k][EDH]; }
        if (rs.fHistory[k][EDH] >= 0.1 - 1e-12) dsm01 = std::max(dsm01, d);
    }
    // full against reduced integration at A: the full integration is first slightly stiffer (q above), then the
    // locking makes it fall below the reduced one
    REAL depart = nan, d33max = 0., excess = 0., excess_at = 0., cross = nan, below = nan;
    for (size_t k = 0; k < rr.fHistory.size() && k < rf.fHistory.size(); ++k) {
        const REAL ds = rf.fHistory[k][EQA] - rr.fHistory[k][EQA], d = std::fabs(ds), x = rr.fHistory[k][EDH];
        if (std::isnan(depart) && d > 1.) depart = x;
        d33max = std::max(d33max, d);
        if (ds > excess) { excess = ds; excess_at = x; }
        if (std::isnan(cross) && excess > 0. && ds <= 0.) cross = x;
        if (std::isnan(below) && ds < -1.) below = x;
    }
    const auto pkf = Peak(rf, EQA);
    // spread of q at the integration points of the element that contains A (largest along the test)
    auto spread = [](const TResult &r) {
        REAL s = 0.;
        for (auto &h : r.fHistory) s = std::max(s, h[EQElMax] - h[EQElMin]);
        return s;
    };
    const auto &er = rr.fHistory.back(), &ef = rf.fHistory.back();
    std::vector<std::tuple<std::string, REAL, std::string>> nums = {
        {"q_A_end_smooth_2x2x2", rs.fHistory.back()[EQA], "q at A at delta/H = 0.6, smooth platen (v0.6: 149.50)"},
        {"q_A_end_material_point_150", mp150.back()[2], "material point, 150 increments (homogeneous solution)"},
        {"q_A_end_material_point_30", mp30.back()[2], "material point, 30 increments (delta delta/H = 0.02)"},
        {"max_diff_smooth_vs_material_point_150", dsm, "largest |q_A(FE smooth) - q(material point)| along the test"},
        {"delta_H_max_diff_smooth_vs_material_point", dsm_at, "delta/H of that largest difference"},
        {"max_diff_smooth_vs_material_point_150_after_0.1", dsm01, "same, for delta/H >= 0.1"},
        {"rms_smooth_2x2x2_vs_abaqus", cs[0], "RMS difference at the 11 digitized points (v0.6 2D: 2.3)"},
        {"max_diff_smooth_2x2x2_vs_abaqus", cs[1], "largest difference at the digitized points (v0.6: 4 kPa)"},
        {"delta_H_max_diff_smooth_vs_abaqus", cs[2], "abscissa of the largest difference (v0.6: 0.13)"},
        {"rms_material_point_30_vs_abaqus", cm30[0], "RMS difference, material point 30 increments (v0.6: 1.1)"},
        {"rms_material_point_150_vs_abaqus", cm150[0], "RMS difference, material point 150 increments"},
        {"rms_rough_2x2x2_vs_abaqus", cr[0], "RMS difference, rough platen 2x2x2"},
        {"max_diff_rough_2x2x2_vs_abaqus", cr[1], "largest difference, rough platen 2x2x2"},
        {"rms_rough_3x3x3_vs_abaqus", cf[0], "RMS difference, rough platen 3x3x3"},
        {"q_A_end_rough_2x2x2", rr.fHistory.back()[EQA], "q at A at the end, rough 2x2x2 (v0.6 3D: 155.94; 2D 155.56)"},
        {"q_A_0.51_rough_2x2x2", At(rr, ar[10][0], EQA), "q at A at the last digitized point 0.5096 (Abaqus 153.09)"},
        {"q_A_end_rough_3x3x3", rf.fHistory.back()[EQA], "q at A at the end, rough 3x3x3 (v0.6 2D 3x3: 145.5)"},
        {"q_A_max_rough_3x3x3", pkf.first, "largest q at A, rough 3x3x3 (v0.6 2D 3x3: 149.0)"},
        {"delta_H_q_A_max_rough_3x3x3", pkf.second, "delta/H of the largest q at A, rough 3x3x3 (v0.6 2D: 0.40)"},
        {"delta_H_departure_3x3x3", depart, "first delta/H with |q_A(3x3x3) - q_A(2x2x2)| > 1 kPa (v0.6 2D: 0.17)"},
        {"max_diff_3x3x3_vs_2x2x2", d33max, "largest |q_A(3x3x3) - q_A(2x2x2)|"},
        {"max_excess_3x3x3_over_2x2x2", excess, "largest q_A(3x3x3) - q_A(2x2x2) (full integration stiffer)"},
        {"delta_H_max_excess_3x3x3", excess_at, "delta/H of that largest excess (v0.6 2D: 1.31 kPa at 0.228)"},
        {"delta_H_crossing_3x3x3", cross, "first delta/H with q_A(3x3x3) <= q_A(2x2x2) after it (v0.6 2D: 0.332)"},
        {"delta_H_below_1kPa_3x3x3", below, "first delta/H with q_A(3x3x3) < q_A(2x2x2) - 1 kPa (v0.6 2D: 0.364)"},
        {"q_elA_min_end_rough_2x2x2", er[EQElMin], "smallest q at the 8 points of the element that contains A, end"},
        {"q_elA_mean_end_rough_2x2x2", er[EQElMean], "mean q at the 8 points of the element that contains A, end"},
        {"q_elA_max_end_rough_2x2x2", er[EQElMax], "largest q at the 8 points of the element that contains A, end"},
        {"q_elA_min_end_rough_3x3x3", ef[EQElMin], "smallest q at the 27 points of the element that contains A, end"},
        {"q_elA_mean_end_rough_3x3x3", ef[EQElMean], "mean q at the 27 points of the element that contains A, end"},
        {"q_elA_max_end_rough_3x3x3", ef[EQElMax], "largest q at the 27 points of the element that contains A, end"},
        {"q_elA_spread_max_rough_2x2x2", spread(rr), "largest (max - min) of q in the element of A along the test"},
        {"q_elA_spread_max_rough_3x3x3", spread(rf), "largest (max - min) of q in the element of A along the test"},
        {"sigma_a_end_smooth_2x2x2", rs.fHistory.back()[ESigmaA], "average axial stress on the platen (v0.6 3D: 249.5)"},
        {"sigma_a_end_rough_2x2x2", rr.fHistory.back()[ESigmaA], "average axial stress on the platen (v0.6: 250.8)"},
        {"sigma_a_end_rough_3x3x3", rf.fHistory.back()[ESigmaA], "average axial stress on the platen (v0.6 2D: 251.0)"},
        {"max_pw_end_smooth_2x2x2", rs.fHistory.back()[EMaxPw], "largest excess pore pressure at the end (kPa)"},
        {"max_pw_end_rough_2x2x2", rr.fHistory.back()[EMaxPw], "largest excess pore pressure at the end (v0.6: 4e-4)"},
        {"max_pw_end_rough_3x3x3", rf.fHistory.back()[EMaxPw], "largest excess pore pressure at the end (kPa)"},
        {"evaluations_smooth_2x2x2", rs.fMeanEvaluations, "mean evaluations per increment (v0.6 3D: 2.48)"},
        {"evaluations_rough_2x2x2", rr.fMeanEvaluations, "mean evaluations per increment (v0.6 3D: 2.73)"},
        {"evaluations_rough_3x3x3", rf.fMeanEvaluations, "mean evaluations per increment (v0.6 2D 3x3: 2.97)"},
        {"bisections_total", REAL(rs.fNBisections + rr.fNBisections + rf.fNBisections), "bisections of the 3 runs"}};
    WriteNumbers("abaqus_fe_numbers.csv", nums);
    std::cout << "  smooth: end q_A = " << rs.fHistory.back()[EQA] << " (material point 150 increments "
              << mp150.back()[2] << "; largest difference " << dsm << " kPa at delta/H = " << dsm_at << ", " << dsm01
              << " kPa beyond 0.1); RMS difference to Abaqus " << cs[0] << " kPa (largest " << cs[1]
              << " at " << cs[2] << "), material point 30 increments " << cm30[0] << " kPa (v0.6: 2.3 and 1.1)\n";
    std::cout << "  rough 2x2x2: end q_A = " << rr.fHistory.back()[EQA] << " (Python 3D 155.94), sigma_a = "
              << rr.fHistory.back()[ESigmaA] << " (250.8), RMS to Abaqus " << cr[0] << " kPa\n";
    std::cout << "  rough 3x3x3: max q_A = " << pkf.first << " at delta/H = " << pkf.second << ", end "
              << rf.fHistory.back()[EQA] << ", departs from 2x2x2 by more than 1 kPa at delta/H = " << depart
              << " (v0.6 2D: 149.0 at 0.40, 145.5, 0.17): first above it (largest excess " << excess << " kPa at "
              << excess_at << "), below it from " << cross << " and by more than 1 kPa from " << below
              << " (2D: 1.31 at 0.228, 0.332, 0.364);"
              << " sigma_a = " << rf.fHistory.back()[ESigmaA] << " (2D 251.0)\n";
    std::cout << "  q at the integration points of the element that contains A at the end (A is on its edge x = 5 mm,"
                 " y = 0): 2x2x2 " << er[EQElMin] << " to " << er[EQElMax] << " (mean " << er[EQElMean] << "), 3x3x3 "
              << ef[EQElMin] << " to " << ef[EQElMax] << " (mean " << ef[EQElMean] << ")\n";
}

inline void AbaqusTriaxialConsolidation::TangentOperators() {
    std::cout << "\nSect. 6.7: rough platen (2 x 2 x 2 points, 150 increments) with five tangent operators\n";
    using P = TPZPlasticStepModifiedCamClay;
    const ETangentMode modes[5] = {P::EConsistentTangent, P::EFiniteDifferenceTangent, P::ESymmetricTangent,
                                   P::EContinuumTangent, P::ETransposedTangent};
    std::vector<std::pair<TConfig, TResult>> runs;
    for (auto mode : modes) {
        TConfig cfg;
        cfg.fRough = true;
        cfg.fTangent = mode;
        cfg.fVTKStride = 0;
        cfg.fName = std::string("rough_2x2x2_tangent_") + P::TangentModeName(mode);
        runs.emplace_back(cfg, Run(cfg));
        const TResult &r = runs.back().second;
        std::cout << "  " << std::setw(4) << P::TangentModeName(mode) << ": " << r.fMeanEvaluations << " evaluations per"
                  << " increment (max " << r.fMaxEvaluations << "), total " << r.fNGlobalIterations << ", bisections "
                  << r.fNBisections << ", " << r.fWallTime << " s (processor " << r.fCPUTime << " s), end q_A = "
                  << std::setprecision(10)
                  << r.fHistory.back()[EQA] << std::setprecision(6) << (r.fCompleted ? "" : " (stopped)") << "\n";
    }
    WriteSummary("abaqus_table10.csv", runs);
    const TResult &rd = runs[0].second, &rt = runs[4].second;
    REAL dq = 0.;
    for (auto &rc : runs)
        if (rc.second.fCompleted) dq = std::max(dq, std::fabs(rc.second.fHistory.back()[EQA] - rd.fHistory.back()[EQA]));
    std::cout << "  largest difference of the final q_A with respect to D: " << dq << " kPa (v0.6 2D: 1e-5)\n";
    std::cout << "  (v0.6, axisymmetric 2 x 2 model, Python: D 3.03 (5) 454; fd 3.03 (5) 454; sym 7.11 (9) 1066;"
                 " cont 11.47 (16) 1720; DT 12.23 (16) 1835)\n";
    // evaluations of each increment with each operator
    {
        std::vector<std::vector<REAL>> rows;
        for (size_t k = 1; k < rd.fHistory.size(); ++k) {
            rows.push_back({REAL(k), rd.fHistory[k][EDH]});
            for (auto &rc : runs)
                rows.back().push_back(k < rc.second.fHistory.size() ? rc.second.fHistory[k][EEvaluations] : 0.);
        }
        std::vector<std::string> header = {"increment", "delta_H"};
        for (auto mode : modes) header.push_back(P::TangentModeName(mode));
        mcc::WriteCSV("abaqus_tangents_evaluations.csv", header, rows);
    }
    // Table 9: normalized residual of the iterations of four increments with D and D^T
    std::cout << "Table 9: normalized residual r = |R|/max(|f_ext|,1) in four increments (D | D^T)\n";
    const int incs[4] = {1, 41, 76, 150};
    std::vector<std::vector<REAL>> t9;
    auto entry = [this](const TResult &r, int inc) -> const mcc::TAnalysis::TStepLog * {
        // last converged (sub)increment that ends at the end of the increment inc
        const REAL tend = fTend * inc / 150.;
        for (auto &l : r.fLog)
            if (std::fabs(l.fState.fTime - tend) <= 1e-9 * fTend) return &l;
        return nullptr;
    };
    for (int inc : incs) {
        const auto *a = entry(rd, inc), *b = entry(rt, inc);
        if (!a || !b) continue;
        std::cout << "  increment " << inc << " (delta/H = " << -a->fState.fUc / fH << "):\n";
        for (size_t it = 0; it < std::max(a->fResiduals.size(), b->fResiduals.size()); ++it) {
            std::cout << "    " << std::setw(2) << it + 1 << "  " << std::scientific << std::setprecision(2);
            if (it < a->fResiduals.size()) std::cout << std::setw(10) << a->fResiduals[it];
            else std::cout << std::setw(10) << " ";
            if (it < b->fResiduals.size()) std::cout << "  " << std::setw(10) << b->fResiduals[it];
            std::cout << std::defaultfloat << std::setprecision(6) << "\n";
        }
        for (int which = 0; which < 2; ++which) {
            const auto *l = which ? b : a;
            for (size_t it = 0; it < l->fResiduals.size(); ++it)
                t9.push_back({REAL(inc), -l->fState.fUc / fH, REAL(which), REAL(it + 1), l->fResiduals[it]});
        }
    }
    mcc::WriteCSV("abaqus_table9.csv", {"increment", "delta_H", "transposed", "iteration", "residual"}, t9);
}

inline void AbaqusTriaxialConsolidation::ToleranceSensitivity() {
    std::cout << "\nSect. 6.7: Table 10 with the tolerance " << fToleranceSensitivity
              << " (rough platen, 2 x 2 x 2 points, 150 increments)\n";
    using P = TPZPlasticStepModifiedCamClay;
    const ETangentMode modes[5] = {P::EConsistentTangent, P::EFiniteDifferenceTangent, P::ESymmetricTangent,
                                   P::EContinuumTangent, P::ETransposedTangent};
    std::vector<std::pair<TConfig, TResult>> runs;
    for (auto mode : modes) {
        TConfig cfg;
        cfg.fRough = true;
        cfg.fTangent = mode;
        cfg.fTolerance = fToleranceSensitivity;
        cfg.fVTKStride = 0;
        cfg.fName = std::string("rough_2x2x2_tolerance_tangent_") + P::TangentModeName(mode);
        runs.emplace_back(cfg, Run(cfg));
        const TResult &r = runs.back().second;
        std::cout << "  " << std::setw(4) << P::TangentModeName(mode) << ": " << r.fMeanEvaluations << " evaluations per"
                  << " increment (max " << r.fMaxEvaluations << "), total " << r.fNGlobalIterations << ", bisections "
                  << r.fNBisections << ", ratio to D " << r.fMeanEvaluations / runs[0].second.fMeanEvaluations << ", "
                  << r.fWallTime << " s (processor " << r.fCPUTime << " s), end q_A = " << std::setprecision(10)
                  << r.fHistory.back()[EQA]
                  << std::setprecision(6) << (r.fCompleted ? "" : " (stopped)") << "\n";
    }
    WriteSummary("abaqus_table10_tolerance.csv", runs);
    std::cout << "  (v0.6, axisymmetric 2 x 2 model with the tolerance 1e-8, Python: D 3.03, sym 7.11, cont 11.47,"
                 " DT 12.23 evaluations per increment, ratios to D 2.35, 3.79 and 4.04)\n";
}

inline void AbaqusTriaxialConsolidation::SmoothRuns(const std::vector<TConfig> &configs, const std::string &file) {
    std::vector<std::pair<TConfig, TResult>> runs;
    std::vector<std::tuple<std::string, REAL, std::string>> nums;
    for (const TConfig &c : configs) {
        runs.emplace_back(c, Run(c));
        const TResult &r = runs.back().second;
        const auto mp = mcc::TriaxialDrained(Model(), c.fP0, fPc0, fV0, fDHend, c.fNSteps);
        const auto pk = Peak(r, EQA), pkp = Peak(r, EQPlaten);
        REAL dmp = 0.;
        for (size_t k = 0; k < r.fHistory.size() && k < mp.size(); ++k)
            dmp = std::max(dmp, std::fabs(r.fHistory[k][EQA] - mp[k][2]));
        const auto &e = r.fHistory.back();
        std::cout << "  " << c.fName << " (" << r.fNElements << " elements, " << r.fNEquations << " degrees of freedom, "
                  << r.fNFreeEquations << " free): peak q_A = "
                  << pk.first << " at delta/H = " << pk.second << "; end p'_A = " << e[EPA] << " q_A = " << e[EQA]
                  << " eps_v = " << e[EEpsV] << "; global q (platen) peak " << pkp.first << " at " << pkp.second
                  << ", end " << e[EQPlaten] << " (material point " << mp.back()[2] << "); largest |q_A - q_mp| = "
                  << dmp << "; " << r.fMeanEvaluations << " evaluations per increment (max " << r.fMaxEvaluations
                  << "), " << r.fNBisections << " bisections, " << r.fWallTime << " s (processor " << r.fCPUTime
                  << " s)\n";
        const std::string t = c.fName;
        nums.push_back({t + "_q_A_peak", pk.first, "largest q at A"});
        nums.push_back({t + "_delta_H_q_A_peak", pk.second, "delta/H of the largest q at A"});
        nums.push_back({t + "_p_A_end", e[EPA], "p' at A at delta/H = 0.6"});
        nums.push_back({t + "_q_A_end", e[EQA], "q at A at delta/H = 0.6"});
        nums.push_back({t + "_eps_v_end", e[EEpsV],
                        "eps_v from the outer top node at delta/H = 0.6 (exact only for a homogeneous deformation)"});
        nums.push_back({t + "_q_platen_peak", pkp.first, "largest global q = sigma_a - p'0 (platen force)"});
        nums.push_back({t + "_delta_H_q_platen_peak", pkp.second, "delta/H of the largest global q"});
        nums.push_back({t + "_q_platen_end", e[EQPlaten], "global q = sigma_a - p'0 at delta/H = 0.6"});
        nums.push_back({t + "_q_material_point_end", mp.back()[2], "material point (same increments) at 0.6"});
        nums.push_back({t + "_max_diff_q_A_vs_material_point", dmp, "largest |q_A - q(material point)|"});
        nums.push_back({t + "_evaluations", r.fMeanEvaluations, "mean evaluations per increment"});
        nums.push_back({t + "_max_evaluations", REAL(r.fMaxEvaluations), "largest evaluations of an increment"});
        nums.push_back({t + "_bisections", REAL(r.fNBisections), "bisections"});
        nums.push_back({t + "_wall_time_s", r.fWallTime, "wall time of the run (s)"});
        nums.push_back({t + "_cpu_time_s", r.fCPUTime, "processor time of the run (s)"});
    }
    WriteSummary("abaqus_" + file + "_summary.csv", runs);
    WriteNumbers("abaqus_" + file + "_numbers.csv", nums);
}

inline void AbaqusTriaxialConsolidation::States() {
    std::cout << "\nFinite element solutions of the two initial states (smooth platen, 600 increments, 48 elements)\n";
    TConfig c100, c20;
    c100.fNSteps = c20.fNSteps = 600;
    c100.fVTKStride = c20.fVTKStride = 4; // 151 states
    c20.fP0 = 20.;
    c100.fName = "smooth_600_p0_100";
    c20.fName = "smooth_600_p0_20";
    SmoothRuns({c100, c20}, "states");
    std::cout << "  (v0.6, axisymmetric 2 x 4 mesh: p'0 = 100: peak 149.54 at 0.6, end p' = 149.85, eps_v = 0.0723, 2.17"
                 " evaluations;\n   p'0 = 20: peak 54.61 at 0.021, end p' = 36.99 q = 37.18 eps_v = -0.681, 2.09"
                 " evaluations; global q at the end 33.4 against 30.1 at the material point)\n";
}

inline void AbaqusTriaxialConsolidation::Softening() {
    std::cout << "\nSoftening state p'0 = 20 kPa with the refined mesh 4 x 4 x 8 (smooth platen, 600 increments)\n";
    TConfig c;
    c.fNSteps = 600;
    c.fVTKStride = 12; // 51 states
    c.fP0 = 20.;
    c.fNc = c.fNr = 4;
    c.fNz = 8;
    c.fOptimizeBandwidth = true; // 5727 free equations: Sloan renumbering, about 4 times faster (results unchanged)
    c.fName = "smooth_600_p0_20_mesh4x4x8";
    SmoothRuns({c}, "softening");
    std::cout << "  (v0.6, axisymmetric 4 x 8 mesh: global q at the end 33.6 against 33.4 with 2 x 4 and 30.1 at the"
                 " material point)\n";
}

inline void AbaqusTriaxialConsolidation::RunAll(const std::vector<std::string> &parts) {
    std::cout << std::setprecision(6);
    const std::vector<std::string> known = {"all", "mesh", "mp", "fe", "tangents", "tolerance", "states", "softening"};
    for (auto &p : parts)
        if (std::find(known.begin(), known.end(), p) == known.end())
            std::cout << "unknown part \"" << p << "\" (parts: mesh mp fe tangents tolerance states softening all)\n";
    auto want = [&parts](const std::string &p) {
        return std::find(parts.begin(), parts.end(), p) != parts.end() ||
               std::find(parts.begin(), parts.end(), "all") != parts.end();
    };
    if (want("mesh")) WriteMeshes();
    if (want("mp")) {
        MaterialPointStates();
        MaterialPointIncrements();
    }
    if (want("fe")) FiniteElementModels();
    if (want("tangents")) TangentOperators();
    if (want("tolerance")) ToleranceSensitivity();
    if (want("states")) States();
    if (want("softening")) Softening();
}
