/**
 * @file FLAC3DTriaxial.h
 * @brief Sect. 6.2 of the article: drained and undrained triaxial tests of the FLAC3D verification
 * problem with a single axisymmetric Q8-Q4 element (Fig. 5 and Table 4).
 */
#pragma once

#include "MCCPaperTools.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "pzstepsolver.h"
#include <iostream>

/**
 * @ingroup mccpaper
 * @brief Triaxial tests of the FLAC3D verification problem on a single u-p element.
 *
 * A 1 m x 1 m axisymmetric Q8-Q4 element with 2 x 2 integration points: radial displacement fixed on
 * the axis (left), vertical displacement fixed at the base, total cell pressure p'0 = 5 kPa on the
 * lateral face and vertical displacement of the top controlled. Drained tests: pore pressure prescribed
 * as zero at the four vertices, eps_a up to 50% in 500 increments. Undrained tests: no flow (k = 0,
 * Dt = 0), pore fluid with Kw = 2e4 kPa (M_B = Kw/n), eps_a up to 10% in 400 increments.
 * Material: M = 1.02, lambda = 0.2, kappa = 0.05, v_lambda = 3.32, G = 250 kPa (Table 1).
 *
 * The class follows the structure of the NeoPZ examples: geometric mesh, computational meshes
 * (displacement, pore pressure and multiphysics), analysis with structural matrix and solver,
 * incremental solution and post-processing.
 */
class FLAC3DTriaxial {
public:
    /** @brief Material and boundary ids */
    enum { EMatId = 1, EBottom = -1, ERight = -2, ETop = -3, ELeft = -4, EPBottom = -11, EPRight = -12,
           EPTop = -13, EPLeft = -14 };

    /** @brief Definition of a test */
    struct TCase {
        REAL fR = 1.6;           ///< overconsolidation ratio p'c0/p'0
        bool fDrained = true;    ///< drained or undrained test
        REAL fEaMax = 0.5;       ///< final axial strain
        int fNSteps = 500;       ///< number of increments
        bool fTransposed = false; ///< use the transpose of the consistent tangent (Sect. 6.6)
        std::string fName = "drained_R1.6";
    };

    /** @brief History of a test: rows (eps_a, p', q, v, mean pore pressure) */
    struct TResult {
        std::vector<std::array<REAL, 5>> fHistory;
        REAL fMeanEvaluations = 0.;
        REAL fV0 = 0., fPc0 = 0.;
    };

    /** @brief Material parameters (Table 1) */
    REAL fM = 1.02, fLambda = 0.2, fKappa = 0.05, fVLambda = 3.32, fG = 250., fP0 = 5., fKw = 2.e4;

    /** @brief Geometric mesh: one quadrilateral with the boundary lines (coincident lines for the pore pressure) */
    TPZGeoMesh *CreateGeoMesh(bool drained);

    /** @brief Computational meshes of displacement, pore pressure and the multiphysics mesh */
    TPZMultiphysicsCompMesh *CreateCompMesh(TPZGeoMesh *gmesh, const TCase &tcase, mcc::TPoroMaterial *&mat);

    /** @brief Runs a test: analysis, structural matrix, solver, increments and post-processing */
    TResult Run(const TCase &tcase);

    /** @brief Runs the four tests of Table 4 and the comparison of the transposed tangent */
    void RunAll();
};

inline TPZGeoMesh *FLAC3DTriaxial::CreateGeoMesh(bool drained) {
    return mcc::CreateRectangleMesh(0., 0., 1., 1., 1, 1, EMatId, [drained](int side, const TPZVec<REAL> &) {
        const int uid[4] = {EBottom, ERight, ETop, ELeft};
        const int pid[4] = {EPBottom, EPRight, EPTop, EPLeft};
        std::vector<int> ids = {uid[side]};
        if (drained) ids.push_back(pid[side]);
        return ids;
    });
}

inline TPZMultiphysicsCompMesh *FLAC3DTriaxial::CreateCompMesh(TPZGeoMesh *gmesh, const TCase &tcase,
                                                              mcc::TPoroMaterial *&mat) {
    std::set<int> bcids = {EBottom, ERight, ETop, ELeft};
    if (tcase.fDrained) bcids.insert({EPBottom, EPRight, EPTop, EPLeft});
    TPZCompMesh *cmeshU = mcc::CreateDisplacementMesh(gmesh, 2, EMatId, bcids);
    TPZCompMesh *cmeshP = mcc::CreatePressureMesh(gmesh, 2, EMatId, bcids);

    const REAL pc0 = tcase.fR * fP0;
    const REAL v0 = mcc::SpecificVolumeNCL(fVLambda, fLambda, fKappa, pc0, fP0);
    const REAL n0 = (v0 - 1.) / v0;

    TPZMultiphysicsCompMesh *mphys = new TPZMultiphysicsCompMesh(gmesh);
    mphys->SetDimModel(2);
    mat = new mcc::TPoroMaterial(EMatId, TPZMatPoroElastoPlasticUPBase::EAxisymmetric);
    mcc::TPlastic model;
    model.SetModifiedCamClay(fM, fLambda, fKappa);
    model.SetConstantShearModulus(fG);
    model.SetTransposedTangent(tcase.fTransposed);
    mat->SetPlasticModel(model);
    mat->SetBiot(1., tcase.fDrained ? 0. : n0 / fKw);
    mat->SetPermeability(0.);
    mat->SetIntegrationOrder(3); // 2 x 2 Gauss points
    mphys->InsertMaterialObject(mat);

    TPZFNMatrix<4, STATE> val1(2, 2, 0.);
    TPZManVector<STATE, 3> val2(2, 0.);
    // axis: u_r = 0
    val1(0, 0) = 1.;
    mphys->InsertMaterialObject(mat->CreateBC(mat, ELeft, TPZMatPoroElastoPlasticUPBase::EDirichletUDirectional, val1, val2));
    // base: u_z = 0
    val1.Zero();
    val1(1, 1) = 1.;
    mphys->InsertMaterialObject(mat->CreateBC(mat, EBottom, TPZMatPoroElastoPlasticUPBase::EDirichletUDirectional, val1, val2));
    // top: u_z controlled
    mphys->InsertMaterialObject(mat->CreateBC(mat, ETop, TPZMatPoroElastoPlasticUPBase::EDirichletUDirectional, val1, val2));
    // lateral face: total cell pressure
    val1.Zero();
    val2[0] = -fP0;
    mphys->InsertMaterialObject(mat->CreateBC(mat, ERight, TPZMatPoroElastoPlasticUPBase::ENeumannU, val1, val2));
    if (tcase.fDrained) {
        TPZManVector<STATE, 3> zero(1, 0.);
        for (int id : {EPBottom, EPRight, EPTop, EPLeft})
            mphys->InsertMaterialObject(mat->CreateBC(mat, id, TPZMatPoroElastoPlasticUPBase::EDirichletP, val1, zero));
    }
    mcc::BuildMultiphysics(mphys, cmeshU, cmeshP, EMatId, bcids);

    // initial state: isotropic effective stress p'0, preconsolidation pc0, specific volume v0
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
    res.fV0 = mcc::SpecificVolumeNCL(fVLambda, fLambda, fKappa, res.fPc0, fP0);

    // analysis: non-symmetric skyline matrix and LU decomposition (no renumbering, see TPZPoroElastoPlasticUPAnalysis)
    mcc::TAnalysis analysis(mphys, mat);
    TPZSkylineNSymStructMatrix<STATE> skyl(mphys);
    skyl.SetNumThreads(0);
    analysis.SetStructuralMatrix(skyl);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELU);
    analysis.SetSolver(step);
    analysis.SetControlledDisplacement(ETop, 1);
    analysis.SetPredictor(true);

    std::vector<mcc::TAnalysis::TLoadState> steps;
    for (int k = 1; k <= tcase.fNSteps; ++k) steps.emplace_back(0., 1., -tcase.fEaMax * k / tcase.fNSteps);

    // top-right vertex (r = 1, z = 1): its vertical displacement is the controlled value
    auto monitor = [&](int, const mcc::TAnalysis::TLoadState &s) {
        auto gps = mcc::GaussPoints(mat, mphys);
        const auto &g = gps[0];
        REAL pmean = 0.;
        int np = 0;
        for (auto &n : analysis.NodesOfMaterials({EMatId})) {
            pmean += analysis.NodalValue(n, 1, 0);
            np++;
        }
        res.fHistory.push_back({-s.fUc, mcc::MeanEffectiveStress(g.fSigma), mcc::DeviatoricStress(g.fSigma),
                                res.fV0 * (1. + g.fEps.I1()), pmean / np});
    };
    analysis.Run(steps, monitor);
    res.fMeanEvaluations = mcc::MeanEvaluations(analysis.StepLog());

    mcc::WriteGaussPointsVTK(mat, mphys, "flac3d_" + tcase.fName + "_gauss.vtk");
    mcc::WriteNodalVTK(analysis, 2, "flac3d_" + tcase.fName + ".vtk", 0);
    std::vector<std::vector<REAL>> rows;
    for (auto &h : res.fHistory) rows.push_back({h[0], h[1], h[2], h[3], h[4]});
    mcc::WriteCSV("flac3d_" + tcase.fName + ".csv", {"eps_a", "p_eff", "q", "v", "u"}, rows);

    mcc::DeleteMeshes(mphys);
    return res;
}

inline void FLAC3DTriaxial::RunAll() {
    std::cout << std::setprecision(6);
    std::cout << "FLAC3D triaxial tests (Sect. 6.2, Table 4): final states (p', q in kPa, v, u in kPa)\n";
    std::vector<TCase> cases(4);
    cases[0] = {1.6, true, 0.5, 500, false, "drained_R1.6"};
    cases[1] = {8.0, true, 0.5, 500, false, "drained_R8"};
    cases[2] = {1.6, false, 0.1, 400, false, "undrained_R1.6"};
    cases[3] = {8.0, false, 0.1, 400, false, "undrained_R8"};
    for (auto &c : cases) {
        TResult r = Run(c);
        const auto &f = r.fHistory.back();
        REAL qmax = 0., eapeak = 0.;
        for (auto &h : r.fHistory)
            if (h[2] > qmax) { qmax = h[2]; eapeak = h[0]; }
        std::cout << c.fName << ": v0 = " << r.fV0 << " pc0 = " << r.fPc0 << " | eps_a = " << f[0] << " p' = " << f[1]
                  << " q = " << f[2] << " v = " << f[3] << " u = " << f[4] << " | peak q = " << qmax
                  << " at eps_a = " << eapeak << " | evaluations per increment = " << r.fMeanEvaluations << std::endl;
    }
    // Sect. 6.6: in the homogeneous test the transposition of the tangent has no effect
    TCase consistent = {1.6, true, 0.05, 50, false, "drained_R1.6_50steps"};
    TCase transposed = {1.6, true, 0.05, 50, true, "drained_R1.6_50steps_transposed"};
    TResult rc = Run(consistent), rt = Run(transposed);
    std::cout << "R = 1.6, 50 increments to 5%: evaluations per increment, consistent D = " << rc.fMeanEvaluations
              << ", transposed D^T = " << rt.fMeanEvaluations << std::endl;
}
