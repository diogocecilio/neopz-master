// Steady seepage after the rapid drawdown (Ceron et al. 2025, Eq. 20-21) on a slope mesh of SlopeGeometry.h:
//   -div(K grad u) = 0, K = diag(k_h, k_v), u = excess pore pressure (p = u + gamma_w y_paper = u - gamma_w y);
//   Dirichlet on the ground surface (toe -3, crest -4, face -6): u = gamma_w max(y, -h_w) (NeoPZ y up), i.e.
//   u = 0 on the crest, u = gamma_w y on the face above the water level and u = -gamma_w h_w below it and on the
//   toe ground; zero flux on the base and the sides (-1, -2, -5). H1 (TPZAnisotropicDarcy), order 2 by default.
// Hydraulic functional J(u) = 1/2 int grad u . K grad u (and J / (k_h H^2 gamma_w^2)) and the check p >= 0.
#ifndef SEEPAGEFE_H
#define SEEPAGEFE_H

#include "AnisotropicDarcy.h"
#include "SeepageForceField.h"
#include "SlopeGeometry.h"

#include "TPZLinearAnalysis.h"
#include "pzcmesh.h"
#include "pzskylstrmatrix.h"
#include "pzstepsolver.h"

#include <chrono>
#include <functional>
#include <map>
#include <memory>

namespace slope {

struct Hydraulics {
    REAL kh = 1., kv = 1.; ///< hydraulic conductivities (only alpha = kh / kv matters for the forces)
    REAL gammaw = 9.81;    ///< unit weight of water (kN/m^3)
    int order = 2;         ///< H1 order of u (1 or 2: PoreField interpolates P2)
};

/// Dirichlet data of Eq. 21 in NeoPZ coordinates, valid on the whole ground surface
inline REAL DrawdownExcessPressure(const SlopeGeometry &g, REAL gammaw, REAL y) {
    return gammaw * std::max(y, -g.hw);
}

using ScalarBC = std::function<REAL(const TPZVec<REAL> &x)>;

/// H1 mesh of the seepage problem: Dirichlet u = ud(x) on dirichlet ids, Neumann g(x) = (K grad u) . n on the
/// neumann ids (no element: zero flux)
inline TPZCompMesh *CreateSeepageCMesh(TPZGeoMesh *gmesh, const Hydraulics &hy, const std::map<int, ScalarBC> &dirichlet,
                                       const std::map<int, ScalarBC> &neumann) {
    auto *cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(2);
    cmesh->SetDefaultOrder(hy.order);
    cmesh->SetAllCreateFunctionsContinuous();
    auto *mat = new TPZAnisotropicDarcy(ESoil, hy.kh, hy.kv);
    cmesh->InsertMaterialObject(mat);
    TPZFMatrix<STATE> val1(1, 1, 0.);
    TPZManVector<STATE, 1> val2(1, 0.);
    for (int type = 0; type < 2; type++)
        for (auto &bcf : type == 0 ? dirichlet : neumann) {
            auto *bc = mat->CreateBC(mat, bcf.first, type == 0 ? TPZAnisotropicDarcy::EDirichlet : TPZAnisotropicDarcy::ENeumann,
                                     val1, val2);
            const ScalarBC fn = bcf.second;
            bc->SetForcingFunctionBC([fn](const TPZVec<REAL> &x, TPZVec<STATE> &v, TPZFMatrix<STATE> &) { v[0] = fn(x); });
            cmesh->InsertMaterialObject(bc);
        }
    cmesh->AutoBuild();
    return cmesh;
}

/// Assembles and solves (skyline Cholesky, bandwidth renumbering); the solution is left in cmesh
inline int64_t SolveSeepage(TPZCompMesh *cmesh) {
    TPZLinearAnalysis an(cmesh, true);
    TPZSkylineStructMatrix<STATE> skyl(cmesh);
    an.SetStructuralMatrix(skyl);
    TPZStepSolver<STATE> step;
    step.SetDirect(ECholesky);
    an.SetSolver(step);
    an.Run();
    return cmesh->NEquations();
}

/// J = 1/2 int grad u . K grad u (edge-midpoint rule: exact for the linear gradient of a P2 field)
inline REAL HydraulicFunctional(const PoreField &pf, const Hydraulics &hy) {
    static const REAL mid[3][2] = {{0.5, 0.}, {0.5, 0.5}, {0., 0.5}};
    REAL J = 0.;
    for (const PoreField::Tri &t : pf.Triangles())
        for (auto &m : mid) {
            REAL u, g[2];
            PoreField::EvaluateTri(t, m[0], m[1], u, g);
            J += 0.5 * (hy.kh * g[0] * g[0] + hy.kv * g[1] * g[1]) * t.area / 3.;
        }
    return J;
}

/// min of the total pore pressure p = u - gamma_w y over the P2 nodes and the centroid of every triangle
inline REAL MinTotalPressure(const PoreField &pf, REAL gammaw, REAL xmin[2]) {
    static const REAL pts[7][2] = {{0., 0.}, {1., 0.}, {0., 1.}, {0.5, 0.}, {0.5, 0.5}, {0., 0.5}, {1. / 3., 1. / 3.}};
    REAL pmin = 1.e300;
    for (const PoreField::Tri &t : pf.Triangles()) {
        const REAL det = t.inv[0] * t.inv[3] - t.inv[1] * t.inv[2];
        for (auto &q : pts) {
            REAL u, g[2];
            PoreField::EvaluateTri(t, q[0], q[1], u, g);
            const REAL x = t.x0[0] + (t.inv[3] * q[0] - t.inv[1] * q[1]) / det;
            const REAL y = t.x0[1] + (-t.inv[2] * q[0] + t.inv[0] * q[1]) / det;
            const REAL p = u - gammaw * y;
            if (p < pmin) pmin = p, xmin[0] = x, xmin[1] = y;
        }
    }
    return pmin;
}

struct SeepageResult {
    std::shared_ptr<const PoreField> field;
    int64_t neq = 0, nel = 0;
    REAL J = 0., Jnorm = 0.; ///< J and J / (k_h H^2 gamma_w^2)
    REAL pmin = 0., xpmin[2] = {0., 0.};
    double seconds = 0.;
};

/// Drawdown seepage on gmesh (the geometric mesh is not modified; cmesh is deleted)
inline SeepageResult DrawdownSeepage(TPZGeoMesh *gmesh, const SlopeGeometry &g, const Hydraulics &hy) {
    const auto t0 = std::chrono::steady_clock::now();
    const ScalarBC ud = [g, hy](const TPZVec<REAL> &x) { return DrawdownExcessPressure(g, hy.gammaw, x[1]); };
    TPZCompMesh *cmesh = CreateSeepageCMesh(gmesh, hy, {{EToe, ud}, {ECrest, ud}, {EFace, ud}}, {});
    SeepageResult r;
    r.neq = SolveSeepage(cmesh);
    auto pf = std::make_shared<PoreField>(cmesh, ESoil, TPZAnisotropicDarcy::EExcessPorePressure);
    r.nel = pf->Triangles().size();
    r.J = HydraulicFunctional(*pf, hy);
    r.Jnorm = r.J / (hy.kh * g.H * g.H * hy.gammaw * hy.gammaw);
    r.pmin = MinTotalPressure(*pf, hy.gammaw, r.xpmin);
    r.field = pf;
    delete cmesh;
    r.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    return r;
}

} // namespace slope

#endif
