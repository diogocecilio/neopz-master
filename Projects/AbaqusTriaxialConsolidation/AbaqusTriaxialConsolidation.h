/**
 * @file AbaqusTriaxialConsolidation.h
 * @brief Sects. 6.4 and 6.6 of the article: consolidation of a triaxial specimen (Abaqus benchmark 1.15.2),
 * Figs. 7 to 9 and Tables 6 and 8.
 */
#pragma once

#include "MCCPaperTools.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "pzstepsolver.h"
#include <iostream>
#include <memory>
#include <string>

/**
 * @ingroup mccpaper
 * @brief Abaqus benchmark 1.15.2: drained triaxial test with displacement control on a cylindrical clay
 * specimen (functions abaqus(), abaqus_mp() and abaqus_states() of gen_data.py and abaqus3d() of
 * gen_data3d.py).
 *
 * By symmetry the upper half of the specimen is modelled (H = 60 mm, R = 20 mm). Material: M = 1,
 * lambda = 0.174, kappa = 0.026, v0 = 2.08 (e0 = 1.08), porous elasticity with nu = 0.3,
 * p'c0 = 116.6 kPa; mobility k = 1.728e-4 m/day / gamma_w = 2e-10 m^2/(kPa s). The specimen starts from
 * the isotropic effective stress p'0 (100 kPa in the benchmark) in equilibrium with the cell pressure
 * kept on the lateral face; the top platen moves down at a constant rate up to delta/H = 0.6 in 400 days,
 * with free drainage at the top; the mid-plane is impermeable with u_z = 0. With the smooth platen only
 * u_z is prescribed on the top, with the rough platen also u_r = 0 (u_x = u_y = 0 in 3D).
 *
 * Models: axisymmetric 2 x 4 Q8-Q4 elements with 2 x 2 (reduced) or 3 x 3 (full) integration, a 4 x 8 mesh
 * for the softening state, and a quarter of the specimen with 48 Hex20-Hex8 elements (2 x 2 x 2 points).
 * Monitored quantities: effective stress at the point A (r = 5 mm, z = 7.5 mm, centroid of the element
 * next to the axis and the mid-plane, interpolated from the integration points), average axial stress on
 * the platen (reaction of the top divided by the area), largest excess pore pressure and volumetric
 * strain from the displacements.
 */
class AbaqusTriaxialConsolidation {
public:
    /** @brief Material and boundary ids (2D: bottom = mid-plane, right = lateral face, top = platen, left = axis;
     * 3D: x = 0, lateral, y = 0, z = 0, z = H) */
    enum { EMatId = 1, EBottom = -1, ERight = -2, ETop = -3, ELeft = -4, EPTop = -13,
           EX0 = -21, ELateral = -22, EY0 = -23, EZ0 = -25, EZH = -26, EPZH = -36 };

    /** @brief Definition of a finite element run */
    struct TConfig {
        bool fRough = false;      ///< rough platen (u_r = 0 on the top)
        int fIntegration = 3;     ///< integration order: 3 -> 2 x 2 (x 2) points, 4 -> 3 x 3 points
        bool fTransposed = false; ///< use the transpose of the consistent tangent (Sect. 6.6)
        int fNSteps = 150;        ///< number of increments
        REAL fP0 = 100.;          ///< initial isotropic effective stress = cell pressure
        int fNr = 2, fNz = 4;     ///< 2D mesh (elements in r and z)
        int fDim = 2;             ///< 2 axisymmetric, 3 quarter of cylinder
        std::string fName = "smooth_2x2";
    };

    /** @brief History of a run: rows (delta/H, p' at A, q at A, axial stress of the platen, max |p_w|, eps_v) */
    struct TResult {
        std::vector<std::array<REAL, 6>> fHistory;
        std::vector<TPZPoroElastoPlasticUPAnalysis::TStepLog> fLog;
        REAL fMeanEvaluations = 0.;
        int fNNodes = 0, fNPressureNodes = 0;
    };

    REAL fR = 0.020, fH = 0.060, fM = 1., fLambda = 0.174, fKappa = 0.026, fV0 = 2.08, fNu = 0.3, fPc0 = 116.6;
    REAL fPerm = 1.728e-4 / 86400. / 10.; ///< mobility k (m^2/(kPa s))
    REAL fTend = 34.56e6;                 ///< 400 days
    REAL fDHend = 0.6;                    ///< final delta/H

    /**
     * @brief Write the VTK file series of every increment of each finite element run (mcc::TVTKSeries) in
     * vtk/\<run name\>, with the series time delta/H (the command line argument "novtk" disables it)
     */
    bool fWriteVTK = true;

    /** @brief The plastic model of the benchmark */
    mcc::TPlastic Model(bool transposed = false) const;

    /** @brief Geometric mesh: axisymmetric rectangle (dim = 2) or quarter of cylinder (dim = 3) */
    TPZGeoMesh *CreateGeoMesh(const TConfig &cfg);

    /** @brief Computational meshes, material and boundary conditions */
    TPZMultiphysicsCompMesh *CreateCompMesh(TPZGeoMesh *gmesh, const TConfig &cfg, mcc::TPoroMaterial *&mat);

    /**
     * @brief Solves a configuration and writes its post-processing files: abaqus_\<name\>.csv, the VTK files of
     * the final state and, with fWriteVTK, the VTK file series of every increment in vtk/\<name\>
     * (abaqus_\<name\>_nodal.vtk.series, abaqus_\<name\>_intpoints.vtk.series, abaqus_\<name\>_gausspoints.vtk.series
     * and abaqus_\<name\>_states.csv with the time and the platen displacement of each state)
     */
    TResult Run(const TConfig &cfg);

    /** @brief Fig. 8: drained tests at a material point from p'0 = 100 and 20 kPa (600 increments) and closed form */
    void MaterialPointStates();

    /** @brief Table 6 (column delta delta/H = 0.02): material point solution with 30, 150 and 3000 increments */
    void MaterialPointIncrements();

    /** @brief Fig. 9, Table 6 and Table 8: 2D models with smooth and rough platens, full and reduced integration */
    void AxisymmetricModels();

    /** @brief Fig. 8 and Sect. 6.4: finite element solutions of the two initial states (600 increments) and mesh 4 x 8 */
    void States();

    /** @brief Fig. 9b and Table 6: quarter of the specimen with Hex20-Hex8 elements */
    void ThreeDimensionalModels();

    /** @brief Runs the parts given by name ("mp", "axi", "states", "3d"; "all" runs everything) */
    void RunAll(const std::vector<std::string> &parts);

protected:
    /** @brief Geometric node closest to x */
    static int64_t ClosestNode(TPZGeoMesh *gmesh, const TPZVec<REAL> &x);
    /** @brief Value of the history column c at delta/H = d (linear interpolation) */
    static REAL At(const TResult &r, REAL d, int c);
};

inline mcc::TPlastic AbaqusTriaxialConsolidation::Model(bool transposed) const {
    mcc::TPlastic model;
    model.SetModifiedCamClay(fM, fLambda, fKappa);
    model.SetPoissonRatio(fNu);
    model.SetPorousElasticity();
    model.SetDefaultSpecificVolume(fV0);
    model.SetTransposedTangent(transposed);
    return model;
}

inline TPZGeoMesh *AbaqusTriaxialConsolidation::CreateGeoMesh(const TConfig &cfg) {
    if (cfg.fDim == 2) {
        return mcc::CreateRectangleMesh(0., 0., fR, fH, cfg.fNr, cfg.fNz, EMatId, [](int side, const TPZVec<REAL> &) {
            const int ids[4] = {EBottom, ERight, ETop, ELeft};
            std::vector<int> out = {ids[side]};
            if (side == 2) out.push_back(EPTop);
            return out;
        });
    }
    const REAL R = fR, H = fH;
    auto marker = [R, H](const std::array<std::array<REAL, 3>, 4> &X) {
        auto all = [&X, R](int c, REAL v) {
            for (auto &p : X)
                if (std::fabs(p[c] - v) > 1e-9 * R) return false;
            return true;
        };
        if (all(0, 0.)) return std::vector<int>{EX0};
        if (all(1, 0.)) return std::vector<int>{EY0};
        if (all(2, 0.)) return std::vector<int>{EZ0};
        if (all(2, H)) return std::vector<int>{EZH, EPZH};
        return std::vector<int>{ELateral};
    };
    return mcc::CreateQuarterCylinderMesh(fR, fH, 2, 2, 4, EMatId, marker, ELateral);
}

inline TPZMultiphysicsCompMesh *AbaqusTriaxialConsolidation::CreateCompMesh(TPZGeoMesh *gmesh, const TConfig &cfg,
                                                                           mcc::TPoroMaterial *&mat) {
    const int dim = cfg.fDim;
    std::set<int> bcids = (dim == 2) ? std::set<int>{EBottom, ERight, ETop, ELeft, EPTop}
                                     : std::set<int>{EX0, ELateral, EY0, EZ0, EZH, EPZH};
    TPZCompMesh *cmeshU = mcc::CreateDisplacementMesh(gmesh, dim, EMatId, bcids);
    TPZCompMesh *cmeshP = mcc::CreatePressureMesh(gmesh, dim, EMatId, bcids);

    TPZMultiphysicsCompMesh *mphys = new TPZMultiphysicsCompMesh(gmesh);
    mphys->SetDimModel(dim);
    mat = new mcc::TPoroMaterial(EMatId, dim == 2 ? TPZMatPoroElastoPlasticUPBase::EAxisymmetric
                                                  : TPZMatPoroElastoPlasticUPBase::EThreeDimensional);
    mat->SetPlasticModel(Model(cfg.fTransposed));
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
    if (dim == 2) {
        directional(ELeft, {0});                                   // axis: u_r = 0
        directional(EBottom, {1});                                 // mid-plane: u_z = 0
        if (cfg.fRough) directional(ETop, {0, 1});                 // rough platen: u_r = 0, u_z controlled
        else directional(ETop, {1});                               // smooth platen: u_z controlled
        val1.Zero();
        TPZManVector<STATE, 3> cell(2, 0.);
        cell[0] = -p0;                                             // cell pressure on the lateral face
        mphys->InsertMaterialObject(mat->CreateBC(mat, ERight, B::ENeumannU, val1, cell));
        mphys->InsertMaterialObject(mat->CreateBC(mat, EPTop, B::EDirichletP, val1, zero));
    } else {
        directional(EX0, {0});                                     // symmetry plane x = 0
        directional(EY0, {1});                                     // symmetry plane y = 0
        directional(EZ0, {2});                                     // mid-plane
        if (cfg.fRough) directional(EZH, {0, 1, 2});
        else directional(EZH, {2});
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
        mphys->InsertMaterialObject(mat->CreateBC(mat, EPZH, B::EDirichletP, val1, zero));
    }
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

inline REAL AbaqusTriaxialConsolidation::At(const TResult &r, REAL d, int c) {
    return mcc::Interpolate(r.fHistory, d, c);
}

inline AbaqusTriaxialConsolidation::TResult AbaqusTriaxialConsolidation::Run(const TConfig &cfg) {
    TPZGeoMesh *gmesh = CreateGeoMesh(cfg);
    mcc::TPoroMaterial *mat = nullptr;
    TPZMultiphysicsCompMesh *mphys = CreateCompMesh(gmesh, cfg, mat);
    const int dim = cfg.fDim;

    mcc::TAnalysis analysis(mphys, mat);
    TPZSkylineNSymStructMatrix<STATE> skyl(mphys);
    skyl.SetNumThreads(0);
    analysis.SetStructuralMatrix(skyl);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELU);
    analysis.SetSolver(step);
    analysis.SetControlledDisplacement(dim == 2 ? ETop : EZH, dim - 1);
    analysis.SetPredictor(true);

    TResult res;
    res.fNNodes = gmesh->NNodes();
    res.fNPressureNodes = analysis.NodesOfMaterials({EMatId}).size();
    TPZManVector<REAL, 3> corner(3, 0.), pointA(3, 0.);
    corner[0] = fR;
    corner[dim - 1] = fH;
    const int64_t iR = ClosestNode(gmesh, corner); // outer top corner (u_r gives the volumetric strain in 2D)
    pointA[0] = 0.005;
    pointA[dim - 1] = 0.0075;
    const REAL area = (dim == 2) ? M_PI * fR * fR : M_PI * fR * fR / 4.;
    std::vector<mcc::TAnalysis::TLoadState> steps;
    for (int k = 1; k <= cfg.fNSteps; ++k)
        steps.emplace_back(fTend * k / cfg.fNSteps, 1., -fDHend * fH * k / cfg.fNSteps);
    const std::string prefix = "abaqus_" + cfg.fName;
    // VTK file series of every increment (series time delta/H)
    std::unique_ptr<mcc::TVTKSeries> vtk;
    if (fWriteVTK)
        vtk = std::make_unique<mcc::TVTKSeries>(mphys, mat, "vtk/" + cfg.fName, prefix, "delta_H",
                                                std::vector<std::string>{"t", "uc"});
    auto monitor = [&](int istep, const mcc::TAnalysis::TLoadState &s) {
        const TPZTensor<REAL> sig = mcc::StressAtPoint(mat, mphys, pointA);
        const REAL dH = fDHend * s.fTime / fTend;
        const REAL sa = (s.fTime > 0.) ? -analysis.Reaction({dim == 2 ? ETop : EZH}, dim - 1) / area : cfg.fP0;
        REAL pmax = 0.;
        for (auto n : analysis.NodesOfMaterials({EMatId})) pmax = std::max(pmax, std::fabs(analysis.NodalValue(n, 1, 0)));
        const REAL ev = dH - 2. * analysis.NodalValue(iR, 0, 0) / fR; // exact with the smooth platen (homogeneous)
        res.fHistory.push_back({dH, mcc::MeanEffectiveStress(sig), mcc::DeviatoricStress(sig), sa, pmax, ev});
        if (vtk) vtk->Write(dH, {s.fTime, s.fUc});
        if (istep == cfg.fNSteps) {
            mcc::WriteNodalVTK(analysis, dim, prefix + ".vtk", 0);
            mcc::WriteGaussPointsVTK(mat, mphys, prefix + "_gauss.vtk");
        }
    };
    if (!analysis.Run(steps, monitor)) std::cout << "  " << cfg.fName << ": the analysis stopped before the end\n";
    res.fLog = analysis.StepLog();
    res.fMeanEvaluations = mcc::MeanEvaluations(res.fLog);
    std::vector<std::vector<REAL>> rows;
    for (auto &h : res.fHistory) rows.push_back({h[0], h[1], h[2], h[3], h[4], h[5]});
    mcc::WriteCSV(prefix + ".csv", {"delta_H", "p_A", "q_A", "sigma_a_platen", "max_pw", "eps_v"}, rows);
    vtk.reset(); // the post-processing meshes refer to the meshes of the run
    mcc::DeleteMeshes(mphys);
    return res;
}

inline void AbaqusTriaxialConsolidation::MaterialPointStates() {
    std::cout << "\nFig. 8: drained test at a material point, p'c0 = 116.6 kPa (600 increments) and closed form\n";
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

inline void AbaqusTriaxialConsolidation::AxisymmetricModels() {
    std::cout << "\nAxisymmetric models (2 x 4 Q8-Q4, 150 increments): q at A (kPa), Table 6\n";
    TConfig smooth{false, 3, false, 150, 100., 2, 4, 2, "smooth_2x2"};
    TConfig rough{true, 3, false, 150, 100., 2, 4, 2, "rough_2x2"};
    TConfig roughfull{true, 4, false, 150, 100., 2, 4, 2, "rough_3x3"};
    TResult rs = Run(smooth), rr = Run(rough), rf = Run(roughfull);
    auto mp = mcc::TriaxialDrained(Model(), 100., fPc0, fV0, fDHend, 30);
    // rows of Table 6: abscissas of the digitized Abaqus curves (smooth and rough platens)
    const REAL dH[7] = {0.03, 0.13, 0.22, 0.32, 0.41, 0.51, 0.60};
    const REAL dHs[7] = {0.0297, 0.1301, 0.2189, 0.3192, 0.4081, 0.5086, 0.60};
    const REAL dHr[7] = {0.0304, 0.1294, 0.2201, 0.3191, 0.4097, 0.5096, 0.60};
    const REAL aS[7] = {59.7, 109.9, 131.1, 141.1, 145.8, 148.1, -1}, tS[7] = {60.2, 114.0, 133.8, 143.5, 147.1, 148.8, 149.5},
               tMP[7] = {56.7, 110.1, 131.0, 141.9, 146.2, 148.4, 149.3}, aR[7] = {59.7, 111.8, 134.7, 145.0, 150.1, 153.1, -1},
               t22[7] = {61.6, 116.8, 137.4, 147.2, 151.6, 154.3, 155.6}, t33[7] = {61.6, 117.4, 138.7, 147.4, 149.0, 147.7, 145.5};
    std::cout << " (q interpolated at the abscissas of the digitized Abaqus curves, as in Table 6)\n";
    std::cout << " d/H | Abaqus smooth | smooth (article) | mat. point 0.02 (article) | Abaqus rough | rough 2x2 (article) | rough 3x3 (article)\n";
    std::cout << std::fixed << std::setprecision(1);
    for (int i = 0; i < 7; ++i) {
        std::cout << std::setprecision(2) << std::setw(4) << dH[i] << std::setprecision(1) << " | " << std::setw(13) << (aS[i] > 0 ? std::to_string(aS[i]).substr(0, 5) : "-")
                  << " | " << std::setw(7) << At(rs, dHs[i], 2) << " (" << tS[i] << ")"
                  << " | " << std::setw(7) << mcc::Interpolate(mp, dHs[i], 2) << " (" << tMP[i] << ")"
                  << "           | " << std::setw(12) << (aR[i] > 0 ? std::to_string(aR[i]).substr(0, 5) : "-")
                  << " | " << std::setw(7) << At(rr, dHr[i], 2) << " (" << t22[i] << ")"
                  << "     | " << std::setw(7) << At(rf, dHr[i], 2) << " (" << t33[i] << ")\n";
    }
    std::cout << std::defaultfloat << std::setprecision(6);
    std::cout << "  end: smooth q_A = " << rs.fHistory.back()[2] << " (149.50), rough 2x2 q_A = " << rr.fHistory.back()[2]
              << " (155.56), rough 3x3 q_A = " << rf.fHistory.back()[2] << " (145.5, max 149.0 at 0.40)\n";
    std::cout << "  average axial stress on the platen at the end: rough 2x2 " << rr.fHistory.back()[3] << " (250.8), 3x3 "
              << rf.fHistory.back()[3] << " (251.0); largest excess pore pressure " << rr.fHistory.back()[4]
              << " kPa (4e-4)\n";
    std::cout << "  evaluations of the residual per increment: smooth " << rs.fMeanEvaluations << " (2.65), rough 2x2 "
              << rr.fMeanEvaluations << " (3.03), rough 3x3 " << rf.fMeanEvaluations << " (2.97)\n";
    REAL qmax = 0., dmax = 0.;
    for (auto &h : rf.fHistory)
        if (h[2] > qmax) { qmax = h[2]; dmax = h[0]; }
    std::cout << "  full integration: maximum q_A = " << qmax << " at delta/H = " << dmax << " (149.0 at 0.40, volumetric locking)\n";

    // Sect. 6.6 and Table 8: transposed tangent with the rough platen
    TConfig transposed = rough;
    transposed.fTransposed = true;
    transposed.fName = "rough_2x2_transposed";
    TResult rt = Run(transposed);
    std::cout << "\nSect. 6.6: rough platen, evaluations per increment with D = " << rr.fMeanEvaluations << " (3.0) and D^T = "
              << rt.fMeanEvaluations << " (12.2); end q_A with D^T = " << rt.fHistory.back()[2] << "\n";
    std::cout << "Table 8: normalized residual r = |R|/max(|f_ext|,1) in four increments (D | D^T)\n";
    const int incs[4] = {0, 40, 75, 149};
    for (int k : incs) {
        if (k >= (int)rr.fLog.size() || k >= (int)rt.fLog.size()) continue;
        std::cout << "  increment " << k + 1 << " (delta/H = " << -rr.fLog[k].fState.fUc / fH << "):\n";
        const auto &a = rr.fLog[k].fResiduals, &b = rt.fLog[k].fResiduals;
        for (size_t it = 0; it < std::max(a.size(), b.size()); ++it) {
            std::cout << "    " << std::setw(2) << it + 1 << "  " << std::scientific << std::setprecision(2);
            if (it < a.size()) std::cout << std::setw(10) << a[it]; else std::cout << std::setw(10) << " ";
            if (it < b.size()) std::cout << "  " << std::setw(10) << b[it];
            std::cout << std::defaultfloat << std::setprecision(6) << "\n";
        }
    }
    std::vector<std::vector<REAL>> t8;
    for (int k : incs)
        for (int which = 0; which < 2; ++which) {
            const auto &log = which ? rt.fLog : rr.fLog;
            if (k >= (int)log.size()) continue;
            for (size_t it = 0; it < log[k].fResiduals.size(); ++it)
                t8.push_back({REAL(k + 1), REAL(which), REAL(it + 1), log[k].fResiduals[it]});
        }
    mcc::WriteCSV("abaqus_table8.csv", {"increment", "transposed", "iteration", "residual"}, t8);
}

inline void AbaqusTriaxialConsolidation::States() {
    std::cout << "\nFinite element solutions of the two initial states (smooth platen, 600 increments)\n";
    for (REAL p0 : {100., 20.}) {
        TConfig cfg{false, 3, false, 600, p0, 2, 4, 2, "smooth_600_p0_" + std::to_string(int(p0))};
        TResult r = Run(cfg);
        REAL qpk = 0., dpk = 0.;
        for (auto &h : r.fHistory)
            if (h[2] > qpk) { qpk = h[2]; dpk = h[0]; }
        const auto &e = r.fHistory.back();
        std::cout << "  p'0 = " << p0 << ": peak q_A = " << qpk << " at delta/H = " << dpk << "; end p' = " << e[1]
                  << " q = " << e[2] << " eps_v = " << e[5] << "; evaluations per increment " << r.fMeanEvaluations << "\n";
    }
    std::cout << "  (Python: p'0 = 100: peak 149.54 at 0.6, end p' = 149.85 eps_v = 0.0723, 2.17 evaluations;\n"
                 "   p'0 = 20: peak 54.61 at 0.021, end p' = 36.99 q = 37.18 eps_v = -0.681, 2.09 evaluations)\n";
    std::cout << "Softening state p'0 = 20 kPa: global deviatoric stress from the platen force (sigma_a - p'0)\n";
    const REAL refpk[2] = {54.614, 54.614}, refend[2] = {33.425, 33.635};
    int i = 0;
    for (auto n : {std::pair<int, int>(2, 4), std::pair<int, int>(4, 8)}) {
        TConfig cfg{false, 3, false, 600, 20., n.first, n.second, 2,
                    "smooth_600_p0_20_mesh" + std::to_string(n.first) + "x" + std::to_string(n.second)};
        TResult r = Run(cfg);
        REAL pk = 0.;
        for (auto &h : r.fHistory) pk = std::max(pk, h[3] - 20.);
        std::cout << "  mesh " << n.first << " x " << n.second << ": peak " << pk << " (" << refpk[i] << "), end "
                  << r.fHistory.back()[3] - 20. << " (" << refend[i] << "); material point 30.1 kPa\n";
        ++i;
    }
}

inline void AbaqusTriaxialConsolidation::ThreeDimensionalModels() {
    std::cout << "\nQuarter of the specimen with 48 Hex20-Hex8 elements (2 x 2 x 2 points, 150 increments)\n";
    const REAL dH[7] = {0.03, 0.13, 0.22, 0.32, 0.41, 0.51, 0.60};
    const REAL dHr[7] = {0.0304, 0.1294, 0.2201, 0.3191, 0.4097, 0.5096, 0.60}; // abscissas of Table 6
    const REAL t3d[7] = {61.5, 116.4, 136.9, 146.9, 151.6, 154.5, 155.9};
    for (bool rough : {false, true}) {
        TConfig cfg{rough, 3, false, 150, 100., 2, 4, 3, rough ? "3d_rough" : "3d_smooth"};
        TResult r = Run(cfg);
        std::cout << "  " << (rough ? "rough" : "smooth") << ": " << r.fNNodes << " nodes (321), " << r.fNPressureNodes
                  << " pressure nodes (95); evaluations per increment " << r.fMeanEvaluations
                  << (rough ? " (2.73)" : " (2.48)") << "\n    d/H:  ";
        for (REAL d : dH) std::cout << std::setw(8) << d;
        std::cout << "\n    q_A:  " << std::fixed << std::setprecision(2);
        for (REAL d : dHr) std::cout << std::setw(8) << At(r, d, 2);
        if (rough) {
            std::cout << "\n    (art.)";
            for (REAL v : t3d) std::cout << std::setw(8) << v;
        }
        std::cout << std::defaultfloat << std::setprecision(6) << "\n    end: q_A = " << r.fHistory.back()[2]
                  << (rough ? " (155.94)" : " (149.50)") << ", average axial stress " << r.fHistory.back()[3]
                  << (rough ? " (250.8)" : " (249.5)") << ", max |p_w| = " << r.fHistory.back()[4] << "\n";
    }
}

inline void AbaqusTriaxialConsolidation::RunAll(const std::vector<std::string> &parts) {
    std::cout << std::setprecision(6);
    auto want = [&parts](const std::string &p) {
        return std::find(parts.begin(), parts.end(), p) != parts.end() ||
               std::find(parts.begin(), parts.end(), "all") != parts.end();
    };
    if (want("mp")) {
        MaterialPointStates();
        MaterialPointIncrements();
    }
    if (want("axi")) AxisymmetricModels();
    if (want("states")) States();
    if (want("3d")) ThreeDimensionalModels();
}
