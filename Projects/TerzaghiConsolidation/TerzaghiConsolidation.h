/**
 * @file TerzaghiConsolidation.h
 * @brief Sect. 6.4 of the article: Terzaghi's consolidation of an elastic column with 1 x 1 x 10 Hex20-Hex8
 * elements (Fig. 8 and Table 6).
 */
#pragma once

#include "MCCPaperTools.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "pzstepsolver.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>

/**
 * @ingroup mccpaper
 * @brief Terzaghi's consolidation of a linear elastic column (functions terzaghi() of gen_data.py and the
 * Hex20-Hex8 model of gen_data3d.py).
 *
 * Column [0,1] x [0,1] x [0,H] with H = 10 m, E = 1e4 kPa, nu = 0.25, incompressible constituents
 * (alpha_B = 1, 1/M_B = 0) and mobility k = 1e-6 m^2/(kPa s), so that c_v = k E_oed = 0.012 m^2/s. The column
 * is discretized with 1 x 1 x 10 Hex20-Hex8 elements (quadratic serendipity displacement, trilinear pore
 * pressure) with 3 x 3 x 3 Gauss points. The base z = 0 is fixed (u_z = 0) and impermeable, the four lateral
 * faces have zero normal displacement and are impermeable, and the load q = 10 kPa is applied at t = 0 on the
 * drained top z = H (p_w = 0). The load is applied in an undrained step (Dt = 0), followed by 102 time steps: the
 * 101 times with 20 steps per decade from T = c_v t/H^2 = 1e-5 to 1, which include T = 0.001, 0.01 and 0.1, and
 * T = 0.5 (Table 6).
 *
 * Monitored: the settlement of the top vertex of the edge x = y = 0 and the pore pressure at the eleven
 * vertices of this edge (z = 0, 1, ..., 10 m). The model reproduces the plane strain Q8-Q4 column of the
 * previous version of the article (gen_data.py terzaghi) and its Hex20-Hex8 counterpart (gen_data3d.py): the
 * strain field is one-dimensional, so both discretizations give the same nodal values.
 */
class TerzaghiConsolidation {
public:
    /** @brief Material and boundary ids: faces x = 0, x = 1, y = 0, y = 1, z = 0, z = H and the drained top */
    enum { EMatId = 1, EX0 = -21, EX1 = -22, EY0 = -23, EY1 = -24, EZ0 = -25, EZ1 = -26, EPZ1 = -36 };

    /** @brief Results of the analysis: history of the monitored nodes, work counters and checks of the solution */
    struct TResult {
        /// rows (t, settlement of the top, p_w at the vertices of x = y = 0 sorted by height), initial state first
        std::vector<std::vector<REAL>> fHistory;
        std::vector<REAL> fHeights;   ///< heights z of the monitored pressure vertices
        REAL fMeanEvaluations = 0.;   ///< mean residual evaluations per increment
        int64_t fNGlobalIterations = 0; ///< all the global iterations (NGlobalIterations)
        int64_t fNBisections = 0;     ///< bisected increments
        REAL fWallTime = 0.;          ///< wall time of the incremental solution (s)
        REAL fSpreadW = 0.;           ///< largest difference of the settlement over the four vertices of the top (m)
        REAL fSpreadP = 0.;           ///< largest difference of p_w between the four vertical edges at the same height (kPa)
        int64_t fNEquations = 0;      ///< equations of the multiphysics mesh
        int fNElements = 0;           ///< volume elements
        int fNPoints = 0;             ///< integration points
    };

    /** @name Data of the problem */
    /** @{ */
    REAL fE = 1.e4;  ///< Young's modulus E (kPa)
    REAL fNu = 0.25; ///< Poisson's ratio nu
    REAL fK = 1.e-6; ///< mobility k (m^2/(kPa s))
    REAL fQ = 10.;   ///< load q on the top (kPa)
    REAL fH = 10.;   ///< height H of the column (m)
    int fNz = 10; ///< number of elements along the height
    /** @} */

    /** @brief Oedometric modulus \f$E_{oed}=E(1-\nu)/[(1+\nu)(1-2\nu)]\f$ */
    REAL Eoed() const { return fE * (1. - fNu) / ((1. + fNu) * (1. - 2. * fNu)); }
    /** @brief Consolidation coefficient \f$c_v = k E_{oed}\f$ */
    REAL Cv() const { return fK * Eoed(); }

    /** @brief Geometric mesh of the column: 1 x 1 x fNz trilinear hexahedra with the boundary quadrilaterals */
    TPZGeoMesh *CreateGeoMesh();

    /** @brief Displacement (Hex20), pore pressure (Hex8) and multiphysics meshes, material and boundary conditions */
    TPZMultiphysicsCompMesh *CreateCompMesh(TPZGeoMesh *gmesh, mcc::TPoroMaterial *&mat);

    /** @brief Times of the analysis (the undrained step at t = 0 is added by Run) */
    std::vector<REAL> Times() const;

    /** @brief Solves the model and post-processes it */
    TResult Run();

    /** @brief Runs the model, prints Table 6 and writes the CSV files of Fig. 8 and Table 6 */
    void RunAll();
};

inline TPZGeoMesh *TerzaghiConsolidation::CreateGeoMesh() {
    const REAL H = fH;
    return mcc::CreateBoxMesh(1., 1., H, 1, 1, fNz, EMatId, [H](const std::array<std::array<REAL, 3>, 4> &X) {
        const REAL tol = 1e-9 * H;
        if (mcc::FaceOnPlane(X, 0, 0., tol)) return std::vector<int>{EX0};
        if (mcc::FaceOnPlane(X, 0, 1., tol)) return std::vector<int>{EX1};
        if (mcc::FaceOnPlane(X, 1, 0., tol)) return std::vector<int>{EY0};
        if (mcc::FaceOnPlane(X, 1, 1., tol)) return std::vector<int>{EY1};
        if (mcc::FaceOnPlane(X, 2, 0., tol)) return std::vector<int>{EZ0};
        return std::vector<int>{EZ1, EPZ1}; // top: load and drainage (coincident quadrilaterals)
    });
}

inline TPZMultiphysicsCompMesh *TerzaghiConsolidation::CreateCompMesh(TPZGeoMesh *gmesh, mcc::TPoroMaterial *&mat) {
    const std::set<int> bcids = {EX0, EX1, EY0, EY1, EZ0, EZ1, EPZ1};
    TPZCompMesh *cmeshU = mcc::CreateDisplacementMesh(gmesh, 3, EMatId, bcids); // serendipity Hex20
    TPZCompMesh *cmeshP = mcc::CreatePressureMesh(gmesh, 3, EMatId, bcids);     // trilinear Hex8

    TPZMultiphysicsCompMesh *mphys = new TPZMultiphysicsCompMesh(gmesh);
    mphys->SetDimModel(3);
    mat = new mcc::TPoroMaterial(EMatId, TPZMatPoroElastoPlasticUPBase::EThreeDimensional);
    mcc::TPlastic model;
    model.SetLinearElastic(fE, fNu);
    mat->SetPlasticModel(model);
    mat->SetBiot(1., 0.);       // incompressible constituents
    mat->SetPermeability(fK);
    mat->SetIntegrationOrder(4); // 3 x 3 x 3 Gauss points
    mphys->InsertMaterialObject(mat);

    using B = TPZMatPoroElastoPlasticUPBase;
    TPZFNMatrix<9, STATE> val1(3, 3, 0.);
    TPZManVector<STATE, 3> val2(3, 0.);
    auto directional = [&](int id, int comp) {
        val1.Zero();
        val1(comp, comp) = 1.;
        mphys->InsertMaterialObject(mat->CreateBC(mat, id, B::EDirichletUDirectional, val1, val2));
    };
    directional(EX0, 0); // lateral faces: zero normal displacement (impermeable: no pressure condition)
    directional(EX1, 0);
    directional(EY0, 1);
    directional(EY1, 1);
    directional(EZ0, 2); // base: u_z = 0 (impermeable)
    // top: load q (traction -q e_z, scaled by the load factor) and drained face p_w = 0
    TPZManVector<STATE, 3> load(3, 0.), zero(1, 0.);
    load[2] = -fQ;
    val1.Zero();
    mphys->InsertMaterialObject(mat->CreateBC(mat, EZ1, B::ENeumannU, val1, load));
    mphys->InsertMaterialObject(mat->CreateBC(mat, EPZ1, B::EDirichletP, val1, zero));
    mcc::BuildMultiphysics(mphys, cmeshU, cmeshP, EMatId, bcids);
    // initial state: no effective stress, no pore pressure (the model is elastic)
    mat->InitializeMemory(mphys, [](const TPZVec<REAL> &, TPZElastoPlasticMem &mem) {
        mem.m_sigma.Zero();
        mem.m_elastoplastic_state.fpressure = 0.;
    });
    return mphys;
}

inline std::vector<REAL> TerzaghiConsolidation::Times() const {
    const REAL scale = fH * fH / Cv();
    std::vector<REAL> all;
    for (int i = 0; i <= 100; ++i) all.push_back(std::pow(10., -5. + 5. * i / 100.) * scale);
    for (REAL T : {0.001, 0.01, 0.1, 0.5}) all.push_back(T * scale);
    std::sort(all.begin(), all.end());
    std::vector<REAL> times;
    for (REAL t : all)
        if (times.empty() || std::fabs(t - times.back()) >= 1e-9 * times.back()) times.push_back(t);
    return times;
}

inline TerzaghiConsolidation::TResult TerzaghiConsolidation::Run() {
    TPZGeoMesh *gmesh = CreateGeoMesh();
    mcc::WriteMeshCSV(gmesh, "terzaghi_mesh");
    mcc::TPoroMaterial *mat = nullptr;
    TPZMultiphysicsCompMesh *mphys = CreateCompMesh(gmesh, mat);

    TResult res;
    res.fNEquations = mphys->NEquations();
    mcc::TAnalysis analysis(mphys, mat);
    TPZSkylineNSymStructMatrix<STATE> skyl(mphys);
    skyl.SetNumThreads(0);
    analysis.SetStructuralMatrix(skyl);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELU);
    analysis.SetSolver(step);
    analysis.SetPredictor(false);

    // monitored nodes: the vertices of the vertical edge x = y = 0 (top vertex: settlement); the other three
    // vertical edges and the nodes of the top are used to check that the solution is one-dimensional
    int64_t topnode = -1;
    std::vector<std::pair<REAL, int64_t>> line;
    std::map<std::pair<long, long>, std::vector<int64_t>> edges; // (x, y) of a vertical edge -> vertices
    std::vector<int64_t> top;
    for (int64_t n = 0; n < gmesh->NNodes(); ++n) {
        TPZManVector<REAL, 3> x(3);
        gmesh->NodeVec()[n].GetCoordinates(x);
        edges[{std::lround(x[0]), std::lround(x[1])}].push_back(n);
        if (std::fabs(x[2] - fH) < 1e-9 * fH) top.push_back(n);
        if (std::fabs(x[0]) > 1e-12 || std::fabs(x[1]) > 1e-12) continue;
        line.push_back({x[2], n});
        if (std::fabs(x[2] - fH) < 1e-9 * fH) topnode = n;
    }
    std::sort(line.begin(), line.end());
    for (auto &l : line) res.fHeights.push_back(l.first);
    for (int64_t iel = 0; iel < gmesh->NElements(); ++iel)
        if (gmesh->Element(iel)->MaterialId() == EMatId) res.fNElements++;

    std::vector<mcc::TAnalysis::TLoadState> steps;
    steps.emplace_back(0., 1., 0.); // undrained step: the load is applied with Dt = 0
    for (REAL t : Times()) steps.emplace_back(t, 1., 0.);
    const REAL scale = fH * fH / Cv();
    auto monitor = [&](int istep, const mcc::TAnalysis::TLoadState &s) {
        std::vector<REAL> row = {s.fTime, -analysis.NodalValue(topnode, 0, 2)};
        for (auto &l : line) row.push_back(analysis.NodalValue(l.second, 1, 0));
        res.fHistory.push_back(row);
        // one-dimensional solution: same settlement at all the nodes of the top, same p_w on the four edges
        for (auto n : top) res.fSpreadW = std::max(res.fSpreadW, std::fabs(-analysis.NodalValue(n, 0, 2) - row[1]));
        for (auto &e : edges)
            for (auto n : e.second) {
                if (analysis.NodeEquation(n, 1, 0) < 0) continue; // no pressure connect
                TPZManVector<REAL, 3> x(3);
                gmesh->NodeVec()[n].GetCoordinates(x);
                for (size_t j = 0; j < line.size(); ++j)
                    if (std::fabs(line[j].first - x[2]) < 1e-9 * fH)
                        res.fSpreadP = std::max(res.fSpreadP, std::fabs(analysis.NodalValue(n, 1, 0) - row[2 + j]));
            }
        for (REAL T : {0.001, 0.01, 0.1, 0.5})
            if (std::fabs(s.fTime - T * scale) < 1e-9 * s.fTime) mcc::WriteNodalVTK(analysis, 3, "terzaghi.vtk", istep);
    };
    analysis.ResetCounters();
    const auto start = std::chrono::steady_clock::now();
    if (!analysis.Run(steps, monitor)) std::cerr << "TerzaghiConsolidation: the analysis stopped\n";
    res.fWallTime = std::chrono::duration<REAL>(std::chrono::steady_clock::now() - start).count();
    res.fMeanEvaluations = mcc::MeanEvaluations(analysis.StepLog());
    res.fNGlobalIterations = analysis.NGlobalIterations();
    res.fNBisections = analysis.NBisections();
    res.fNPoints = int(mcc::GaussPoints(mat, mphys).size());

    std::vector<std::string> header = {"t", "settlement"};
    for (REAL h : res.fHeights) header.push_back("p_z" + std::to_string(int(std::lround(h))));
    mcc::WriteCSV("terzaghi_history.csv", header, res.fHistory);
    mcc::DeleteMeshes(mphys);
    return res;
}

inline void TerzaghiConsolidation::RunAll() {
    const auto start = std::chrono::steady_clock::now();
    std::cout << std::setprecision(6);
    const REAL scale = fH * fH / Cv();
    std::cout << "Terzaghi consolidation (Sect. 6.4, Fig. 8, Table 6): 1 x 1 x " << fNz
              << " Hex20-Hex8 elements, 3 x 3 x 3 points; c_v = " << Cv() << " m2/s, E_oed = " << Eoed()
              << " kPa, w_inf = qH/E_oed = " << 1e3 * fQ * fH / Eoed() << " mm\n";
    const TResult r = Run();
    std::cout << r.fNElements << " elements, " << r.fNPoints << " integration points, " << r.fNEquations
              << " equations, " << r.fHistory.size() - 1 << " increments (undrained step and "
              << r.fHistory.size() - 2 << " time steps)\n";
    // Table 6: errors of the normalized pore pressure at the vertices and settlement of the top
    const REAL refErr[4] = {0.075, 0.011, 0.007, 0.015}, refS[4] = {0.366, 0.954, 2.959, 6.289},
               refEx[4] = {0.297, 0.940, 2.974, 6.366};
    std::cout << "Table 6                T:     0.001      0.01       0.1       0.5\n";
    std::vector<std::vector<REAL>> table, isochrones;
    const REAL Ts[4] = {0.001, 0.01, 0.1, 0.5};
    std::array<REAL, 4> err, errz, s, sex;
    for (int i = 0; i < 4; ++i) {
        const REAL t = Ts[i] * scale;
        size_t best = 0;
        for (size_t k = 0; k < r.fHistory.size(); ++k)
            if (std::fabs(r.fHistory[k][0] - t) < std::fabs(r.fHistory[best][0] - t)) best = k;
        const auto &row = r.fHistory[best];
        REAL U = 0.;
        err[i] = 0.;
        errz[i] = 0.;
        for (size_t j = 0; j < r.fHeights.size(); ++j) {
            const REAL z = r.fHeights[j];
            const REAL pex = mcc::TerzaghiPressure((fH - z) / fH, Ts[i], U);
            if (std::fabs(row[2 + j] / fQ - pex) > err[i]) {
                err[i] = std::fabs(row[2 + j] / fQ - pex);
                errz[i] = z;
            }
            isochrones.push_back({Ts[i], z, row[2 + j] / fQ, pex});
        }
        s[i] = row[1];
        sex[i] = fQ * fH / Eoed() * U;
        table.push_back({Ts[i], REAL(best), err[i], errz[i], 1e3 * s[i], 1e3 * sex[i],
                         100. * (s[i] - sex[i]) / sex[i]});
    }
    auto line = [](const char *name, const std::array<REAL, 4> &v, REAL f, const REAL *ref) {
        std::cout << std::left << std::setw(22) << name << std::right;
        for (int i = 0; i < 4; ++i) std::cout << std::setw(10) << std::fixed << std::setprecision(3) << f * v[i];
        std::cout << "   (v0.6:";
        for (int i = 0; i < 4; ++i) std::cout << " " << ref[i];
        std::cout << ")\n" << std::defaultfloat;
    };
    line("max|pw/q - exact|", err, 1., refErr);
    line("settlement (mm)", s, 1e3, refS);
    line("settlement exact (mm)", sex, 1e3, refEx);
    // Python code of v0.6 (gen_data.py terzaghi, plane strain Q8-Q4 column): errors and settlements (mm)
    const REAL pyErr[4] = {0.07542243437519269, 0.011285986117618774, 0.007075274691403122, 0.015360723540644772},
               pyS[4] = {0.3661850953043024, 0.9539622253702426, 2.958902496264679, 6.288605301075821};
    REAL dpy = 0.;
    for (int i = 0; i < 4; ++i) dpy = std::max({dpy, std::fabs(err[i] - pyErr[i]), std::fabs(1e3 * s[i] - pyS[i])});
    std::cout << "largest difference from the Python Q8-Q4 values of v0.6 (errors and settlements in mm): "
              << std::scientific << std::setprecision(2) << dpy << std::defaultfloat << std::setprecision(6) << "\n";
    std::cout << "largest error at z = " << errz[0] << ", " << errz[1] << ", " << errz[2] << ", " << errz[3]
              << " m; settlement at T = 0.5: " << std::setprecision(3) << 100. * (s[3] - sex[3]) / sex[3]
              << "% of the exact value (v0.6: 1.2% below)\n" << std::setprecision(6);
    // undrained step: oscillation next to the drained face (v0.6: 1.27 at y = 9 m and 0.93 at y = 8 m)
    const auto &u = r.fHistory[1];
    std::cout << "undrained step: p_w/q at z = 9 m: " << u[2 + 9] / fQ << " (v0.6 1.27), z = 8 m: " << u[2 + 8] / fQ
              << " (0.93); settlement " << 1e3 * u[1] << " mm\n";
    std::cout << "one-dimensional solution: settlement spread over the top vertices " << r.fSpreadW
              << " m, pore pressure spread between the four vertical edges " << r.fSpreadP << " kPa\n";
    std::cout << "evaluations of the residual per increment: " << r.fMeanEvaluations << ", total "
              << r.fNGlobalIterations << ", bisections " << r.fNBisections << "; solution " << r.fWallTime << " s\n";
    // degree of consolidation U(T) (Fig. 8c)
    std::vector<std::vector<REAL>> degree;
    for (auto &row : r.fHistory) {
        const REAL T = row[0] / scale;
        REAL U = 0.;
        if (T > 0.) mcc::TerzaghiPressure(0., T, U);
        degree.push_back({T, row[1] * Eoed() / (fQ * fH), U});
    }
    // exact series solution (B.9) for the lines of Fig. 8b and 8c: isochrones at 201 heights and U at 301 values of T
    std::vector<std::vector<REAL>> iso, deg;
    for (REAL T : Ts)
        for (int k = 0; k <= 200; ++k) {
            const REAL z = fH * k / 200.;
            REAL U;
            iso.push_back({T, z, mcc::TerzaghiPressure((fH - z) / fH, T, U)});
        }
    for (int k = 0; k <= 300; ++k) {
        const REAL T = std::pow(10., -5. + 5. * k / 300.);
        REAL U;
        mcc::TerzaghiPressure(0., T, U);
        deg.push_back({T, U});
    }
    mcc::WriteCSV("terzaghi_exact_isochrones.csv", {"T", "z", "pw_over_q"}, iso);
    mcc::WriteCSV("terzaghi_exact_degree.csv", {"T", "U"}, deg);
    mcc::WriteCSV("terzaghi_table6.csv",
                  {"T", "increment", "max_err_pw", "z_max_err", "settlement_mm", "settlement_exact_mm",
                   "settlement_diff_percent"},
                  table);
    mcc::WriteCSV("terzaghi_isochrones.csv", {"T", "z", "pw_over_q", "exact"}, isochrones);
    mcc::WriteCSV("terzaghi_degree.csv", {"T", "U_numerical", "U_exact"}, degree);
    mcc::WriteCSV("terzaghi_summary.csv",
                  {"H", "E", "nu", "k", "q", "cv", "Eoed", "w_inf", "elements", "points", "equations", "increments",
                   "mean_evaluations", "global_iterations", "bisections", "wall_time_s", "spread_settlement",
                   "spread_pw", "pw_z9_undrained", "pw_z8_undrained", "settlement_undrained"},
                  {{fH, fE, fNu, fK, fQ, Cv(), Eoed(), fQ * fH / Eoed(), REAL(r.fNElements), REAL(r.fNPoints),
                    REAL(r.fNEquations), REAL(r.fHistory.size() - 1), r.fMeanEvaluations, REAL(r.fNGlobalIterations),
                    REAL(r.fNBisections), r.fWallTime, r.fSpreadW, r.fSpreadP, u[2 + 9] / fQ, u[2 + 8] / fQ, u[1]}});
    std::cout << "Files: terzaghi_history.csv, terzaghi_isochrones.csv, terzaghi_degree.csv, terzaghi_exact_isochrones.csv, "
                 "terzaghi_exact_degree.csv, terzaghi_table6.csv, terzaghi_summary.csv, terzaghi_mesh_*.csv, "
                 "terzaghi.scal_vec.<step>.vtk\ntotal run time "
              << std::chrono::duration<REAL>(std::chrono::steady_clock::now() - start).count() << " s\n";
}
