/**
 * @file TerzaghiConsolidation.h
 * @brief Sect. 6.3 of the article: Terzaghi's consolidation of an elastic column with Q8-Q4 elements in
 * plane strain and with Hex20-Hex8 elements (Fig. 6 and Table 5).
 */
#pragma once

#include "MCCPaperTools.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "pzstepsolver.h"
#include <iostream>

/**
 * @ingroup mccpaper
 * @brief Terzaghi's consolidation of a linear elastic column (functions terzaghi() of gen_data.py and the
 * Hex20-Hex8 model of gen_data3d.py).
 *
 * Column of height H = 10 m and width 1 m, E = 1e4 kPa, nu = 0.25, incompressible constituents
 * (alpha_B = 1, 1/M_B = 0) and mobility k = 1e-6 m^2/(kPa s), so that c_v = k E_oed = 0.012 m^2/s.
 * The base is fixed and impermeable, the lateral faces are impermeable and restrained horizontally, and
 * the load q = 10 kPa is applied at t = 0 on the drained top. The load is applied in an undrained step
 * (Dt = 0), followed by 101 time steps with 20 steps per decade from T = c_v t/H^2 = 1e-5 to 1 (plus the
 * times T = 0.001, 0.01, 0.1 and 0.5). The 2D model has 1 x 10 Q8-Q4 elements (3 x 3 Gauss points),
 * the 3D model 1 x 1 x 10 Hex20-Hex8 elements (3 x 3 x 3 points).
 */
class TerzaghiConsolidation {
public:
    /** @brief Material and boundary ids (2D: bottom, right, top, left; 3D: x=0, x=1, y=0, y=1, z=0, z=H) */
    enum { EMatId = 1, EBottom = -1, ERight = -2, ETop = -3, ELeft = -4, EPTop = -13,
           EX0 = -21, EX1 = -22, EY0 = -23, EY1 = -24, EZ0 = -25, EZ1 = -26, EPZ1 = -36 };

    /** @brief Results of a model: rows (t, settlement of the top, pore pressures at the vertices of x = 0 sorted by height) */
    struct TResult {
        std::vector<std::vector<REAL>> fHistory;
        std::vector<REAL> fHeights; ///< heights of the monitored pressure vertices
        REAL fMeanEvaluations = 0.;
    };

    REAL fE = 1.e4, fNu = 0.25, fK = 1.e-6, fQ = 10., fH = 10.;

    /** @brief Oedometric modulus \f$E_{oed}=E(1-\nu)/[(1+\nu)(1-2\nu)]\f$ */
    REAL Eoed() const { return fE * (1. - fNu) / ((1. + fNu) * (1. - 2. * fNu)); }
    /** @brief Consolidation coefficient \f$c_v = k E_{oed}\f$ */
    REAL Cv() const { return fK * Eoed(); }

    /** @brief Geometric mesh of the column: 1 x 10 quadrilaterals (dim = 2) or 1 x 1 x 10 hexahedra (dim = 3) */
    TPZGeoMesh *CreateGeoMesh(int dim);

    /** @brief Displacement, pore pressure and multiphysics meshes, material and boundary conditions */
    TPZMultiphysicsCompMesh *CreateCompMesh(TPZGeoMesh *gmesh, int dim, mcc::TPoroMaterial *&mat);

    /** @brief Times of the analysis (undrained step at t = 0 followed by the consolidation steps) */
    std::vector<REAL> Times() const;

    /** @brief Solves the model of dimension dim and post-processes it */
    TResult Run(int dim);

    /** @brief Runs the 2D and 3D models and prints Table 5 */
    void RunAll();
};

inline TPZGeoMesh *TerzaghiConsolidation::CreateGeoMesh(int dim) {
    if (dim == 2) {
        return mcc::CreateRectangleMesh(0., 0., 1., fH, 1, 10, EMatId, [](int side, const TPZVec<REAL> &) {
            const int ids[4] = {EBottom, ERight, ETop, ELeft};
            std::vector<int> out = {ids[side]};
            if (side == 2) out.push_back(EPTop);
            return out;
        });
    }
    const REAL H = fH;
    return mcc::CreateBoxMesh(1., 1., H, 1, 1, 10, EMatId, [H](const std::array<std::array<REAL, 3>, 4> &X) {
        auto all = [&X, H](int c, REAL v) {
            for (auto &p : X)
                if (std::fabs(p[c] - v) > 1e-9 * H) return false;
            return true;
        };
        if (all(0, 0.)) return std::vector<int>{EX0};
        if (all(0, 1.)) return std::vector<int>{EX1};
        if (all(1, 0.)) return std::vector<int>{EY0};
        if (all(1, 1.)) return std::vector<int>{EY1};
        if (all(2, 0.)) return std::vector<int>{EZ0};
        return std::vector<int>{EZ1, EPZ1};
    });
}

inline TPZMultiphysicsCompMesh *TerzaghiConsolidation::CreateCompMesh(TPZGeoMesh *gmesh, int dim,
                                                                      mcc::TPoroMaterial *&mat) {
    std::set<int> bcids = (dim == 2) ? std::set<int>{EBottom, ERight, ETop, ELeft, EPTop}
                                     : std::set<int>{EX0, EX1, EY0, EY1, EZ0, EZ1, EPZ1};
    TPZCompMesh *cmeshU = mcc::CreateDisplacementMesh(gmesh, dim, EMatId, bcids);
    TPZCompMesh *cmeshP = mcc::CreatePressureMesh(gmesh, dim, EMatId, bcids);

    TPZMultiphysicsCompMesh *mphys = new TPZMultiphysicsCompMesh(gmesh);
    mphys->SetDimModel(dim);
    mat = new mcc::TPoroMaterial(EMatId, dim == 2 ? TPZMatPoroElastoPlasticUPBase::EPlaneStrain
                                                  : TPZMatPoroElastoPlasticUPBase::EThreeDimensional);
    mcc::TPlastic model;
    model.SetLinearElastic(fE, fNu);
    mat->SetPlasticModel(model);
    mat->SetBiot(1., 0.);
    mat->SetPermeability(fK);
    mat->SetIntegrationOrder(4); // 3 x 3 (x 3) Gauss points
    mphys->InsertMaterialObject(mat);

    using B = TPZMatPoroElastoPlasticUPBase;
    TPZFNMatrix<9, STATE> val1(dim, dim, 0.);
    TPZManVector<STATE, 3> val2(dim, 0.);
    auto directional = [&](int id, int comp) {
        val1.Zero();
        val1(comp, comp) = 1.;
        mphys->InsertMaterialObject(mat->CreateBC(mat, id, B::EDirichletUDirectional, val1, val2));
    };
    TPZManVector<STATE, 3> load(dim, 0.), zero(1, 0.);
    load[dim - 1] = -fQ;
    val1.Zero();
    if (dim == 2) {
        directional(ELeft, 0);
        directional(ERight, 0);
        directional(EBottom, 1);
        mphys->InsertMaterialObject(mat->CreateBC(mat, ETop, B::ENeumannU, val1, load));
        mphys->InsertMaterialObject(mat->CreateBC(mat, EPTop, B::EDirichletP, val1, zero));
    } else {
        directional(EX0, 0);
        directional(EX1, 0);
        directional(EY0, 1);
        directional(EY1, 1);
        directional(EZ0, 2);
        val1.Zero();
        mphys->InsertMaterialObject(mat->CreateBC(mat, EZ1, B::ENeumannU, val1, load));
        mphys->InsertMaterialObject(mat->CreateBC(mat, EPZ1, B::EDirichletP, val1, zero));
    }
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

inline TerzaghiConsolidation::TResult TerzaghiConsolidation::Run(int dim) {
    TPZGeoMesh *gmesh = CreateGeoMesh(dim);
    mcc::TPoroMaterial *mat = nullptr;
    TPZMultiphysicsCompMesh *mphys = CreateCompMesh(gmesh, dim, mat);

    mcc::TAnalysis analysis(mphys, mat);
    TPZSkylineNSymStructMatrix<STATE> skyl(mphys);
    skyl.SetNumThreads(0);
    analysis.SetStructuralMatrix(skyl);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELU);
    analysis.SetSolver(step);
    analysis.SetPredictor(false);

    // monitored nodes: top node on the axis x = 0 (and y = 0 in 3D) and the vertices along this line
    TResult res;
    int64_t topnode = -1;
    std::vector<std::pair<REAL, int64_t>> line;
    for (int64_t n = 0; n < gmesh->NNodes(); ++n) {
        TPZManVector<REAL, 3> x(3);
        gmesh->NodeVec()[n].GetCoordinates(x);
        const bool onaxis = std::fabs(x[0]) < 1e-12 && (dim == 2 || std::fabs(x[1]) < 1e-12);
        if (!onaxis) continue;
        const REAL h = x[dim - 1];
        line.push_back({h, n});
        if (std::fabs(h - fH) < 1e-9) topnode = n;
    }
    std::sort(line.begin(), line.end());
    for (auto &l : line) res.fHeights.push_back(l.first);

    std::vector<mcc::TAnalysis::TLoadState> steps;
    steps.emplace_back(0., 1., 0.);
    for (REAL t : Times()) steps.emplace_back(t, 1., 0.);
    const std::string prefix = (dim == 2) ? "terzaghi_q8q4" : "terzaghi_hex20hex8";
    const REAL scale = fH * fH / Cv();
    auto monitor = [&](int istep, const mcc::TAnalysis::TLoadState &s) {
        std::vector<REAL> row = {s.fTime, -analysis.NodalValue(topnode, 0, dim - 1)};
        for (auto &l : line) row.push_back(analysis.NodalValue(l.second, 1, 0));
        res.fHistory.push_back(row);
        for (REAL T : {0.001, 0.01, 0.1, 0.5})
            if (std::fabs(s.fTime - T * scale) < 1e-9 * s.fTime) mcc::WriteNodalVTK(analysis, dim, prefix + ".vtk", istep);
    };
    analysis.Run(steps, monitor);
    res.fMeanEvaluations = mcc::MeanEvaluations(analysis.StepLog());

    std::vector<std::string> header = {"t", "settlement"};
    for (REAL h : res.fHeights) header.push_back("p_y" + std::to_string(int(std::lround(h))));
    mcc::WriteCSV(prefix + "_history.csv", header, res.fHistory);
    mcc::DeleteMeshes(mphys);
    return res;
}

inline void TerzaghiConsolidation::RunAll() {
    std::cout << std::setprecision(6);
    const REAL scale = fH * fH / Cv();
    std::cout << "Terzaghi consolidation (Sect. 6.3): c_v = " << Cv() << " m2/s, E_oed = " << Eoed()
              << " kPa, w_inf = qH/E_oed = " << 1e3 * fQ * fH / Eoed() << " mm\n";
    TResult r2 = Run(2);
    TResult r3 = Run(3);
    // Table 5: errors of the normalized pore pressure at the vertices and settlement of the top
    const REAL refErr[4] = {0.075, 0.011, 0.007, 0.015}, refS[4] = {0.366, 0.954, 2.959, 6.289},
               refEx[4] = {0.297, 0.940, 2.974, 6.366};
    std::cout << "Table 5 (Q8-Q4)        T:     0.001      0.01       0.1       0.5\n";
    std::vector<std::vector<REAL>> table, fig6;
    REAL Ts[4] = {0.001, 0.01, 0.1, 0.5};
    std::array<REAL, 4> err, s, sex;
    for (int i = 0; i < 4; ++i) {
        const REAL t = Ts[i] * scale;
        size_t best = 0;
        for (size_t k = 0; k < r2.fHistory.size(); ++k)
            if (std::fabs(r2.fHistory[k][0] - t) < std::fabs(r2.fHistory[best][0] - t)) best = k;
        const auto &row = r2.fHistory[best];
        REAL U = 0.;
        err[i] = 0.;
        for (size_t j = 0; j < r2.fHeights.size(); ++j) {
            const REAL y = r2.fHeights[j];
            const REAL pex = mcc::TerzaghiPressure((fH - y) / fH, Ts[i], U);
            err[i] = std::max(err[i], std::fabs(row[2 + j] / fQ - pex));
            fig6.push_back({Ts[i], y, row[2 + j] / fQ, pex});
        }
        s[i] = row[1];
        sex[i] = fQ * fH / Eoed() * U;
        table.push_back({Ts[i], err[i], 1e3 * s[i], 1e3 * sex[i]});
    }
    auto line = [](const char *name, const std::array<REAL, 4> &v, REAL f, const REAL *ref) {
        std::cout << std::left << std::setw(22) << name << std::right;
        for (int i = 0; i < 4; ++i) std::cout << std::setw(10) << std::fixed << std::setprecision(3) << f * v[i];
        std::cout << "   (article:";
        for (int i = 0; i < 4; ++i) std::cout << " " << ref[i];
        std::cout << ")\n" << std::defaultfloat;
    };
    line("max|pw/q - exact|", err, 1., refErr);
    line("settlement (mm)", s, 1e3, refS);
    line("settlement exact (mm)", sex, 1e3, refEx);
    // undrained step: oscillation next to the drained face (article: 1.27 at y = 9 m and 0.93 at y = 8 m)
    const auto &u = r2.fHistory[1];
    std::cout << "undrained step: p_w/q at y = 9 m: " << u[2 + 9] / fQ << " (article 1.27), y = 8 m: " << u[2 + 8] / fQ
              << " (0.93); settlement " << 1e3 * u[1] << " mm\n";
    // the Hex20-Hex8 model gives the same results (article: to 5e-11)
    REAL diff = 0.;
    for (size_t k = 0; k < r2.fHistory.size() && k < r3.fHistory.size(); ++k)
        for (size_t j = 1; j < r2.fHistory[k].size(); ++j) diff = std::max(diff, std::fabs(r2.fHistory[k][j] - r3.fHistory[k][j]));
    std::cout << "largest difference between the Hex20-Hex8 and Q8-Q4 histories: " << diff << " (article 5e-11)\n";
    std::cout << "evaluations of the residual per increment: Q8-Q4 " << r2.fMeanEvaluations << ", Hex20-Hex8 "
              << r3.fMeanEvaluations << "\n";
    // degree of consolidation U(T) (Fig. 6b)
    std::vector<std::vector<REAL>> fig6b;
    for (auto &row : r2.fHistory) {
        const REAL T = row[0] / scale;
        REAL U = 0.;
        if (T > 0.) mcc::TerzaghiPressure(0., T, U);
        fig6b.push_back({T, row[1] * Eoed() / (fQ * fH), U});
    }
    mcc::WriteCSV("terzaghi_table5.csv", {"T", "max_err_pw", "settlement_mm", "settlement_exact_mm"}, table);
    mcc::WriteCSV("terzaghi_fig6a.csv", {"T", "y", "pw_over_q", "exact"}, fig6);
    mcc::WriteCSV("terzaghi_fig6b.csv", {"T", "U_numerical", "U_exact"}, fig6b);
}
