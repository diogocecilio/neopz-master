// Steady seepage after the rapid drawdown (Ceron et al. 2025, Eq. 20-21) on a slope mesh of SlopeGeometry.h:
//   -div(K grad u) = 0, K = diag(k_h, k_v), u = excess pore pressure (p = u + gamma_w y_paper = u - gamma_w y);
//   Dirichlet on the ground surface (toe -3, crest -4, face -6): u = gamma_w max(y, -h_w) (NeoPZ y up), i.e.
//   u = 0 on the crest, u = gamma_w y on the face above the water level and u = -gamma_w h_w below it and on the
//   toe ground. Far sides of the truncated half-space (base -1, right -2, left -5): zero flux, u = 0 or
//   u = -gamma_w h_w per side (FarSides, presets of scripts/fe_seepage.py). Default zero_lb: u = 0 (far field still
//   hydrostatic at the crest level) on the left side and the base, zero flux on the right side; on the paper-like box
//   (50 H / 10 H / 30 H) it reproduces the dashed curves of Fig. 5 within 0.22 % (alpha = 1, 2, 4, 10, beta = 15..90),
//   while the fully impermeable box gives a J 3-9 % lower. The iso-lines of u then leave the slope and end on the
//   right side, as in Fig. 4b. H1 (TPZAnisotropicDarcy), order 2 by default.
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

#include <array>
#include <chrono>
#include <functional>
#include <map>
#include <memory>
#include <set>
#include <string>

namespace slope {

/// Condition on a far side of the hydraulic box: zero flux, u = 0 (hydrostatic at the crest level) or
/// u = -gamma_w h_w (hydrostatic at the final water level)
enum class EFarBC { ENoFlow, EZero, EToeLevel };

struct FarSides { ///< default: preset zero_lb
    EFarBC left = EFarBC::EZero, bottom = EFarBC::EZero, right = EFarBC::ENoFlow;

    /// presets of scripts/fe_seepage.py (BC_PRESETS); false if the name is unknown
    bool SetPreset(const std::string &name) {
        const EFarBC N = EFarBC::ENoFlow, Z = EFarBC::EZero, T = EFarBC::EToeLevel;
        static const std::map<std::string, std::array<EFarBC, 3>> presets = {
            {"impermeable", {N, N, N}}, {"zero_lb", {Z, Z, N}},  {"zero_b", {N, Z, N}},
            {"zero_l", {Z, N, N}},      {"zero_lbr", {Z, Z, Z}}, {"toe_r", {Z, N, T}}};
        auto it = presets.find(name);
        if (it == presets.end()) return false;
        left = it->second[0], bottom = it->second[1], right = it->second[2];
        return true;
    }
    static bool Parse(const std::string &s, EFarBC &bc) {
        if (s == "noflow") bc = EFarBC::ENoFlow;
        else if (s == "zero") bc = EFarBC::EZero;
        else if (s == "toe") bc = EFarBC::EToeLevel;
        else return false;
        return true;
    }
    static const char *Name(EFarBC bc) {
        return bc == EFarBC::ENoFlow ? "noflow" : (bc == EFarBC::EZero ? "zero" : "toe");
    }
};

struct Hydraulics {
    REAL kh = 1., kv = 1.; ///< hydraulic conductivities (only alpha = kh / kv matters for the forces)
    REAL gammaw = 9.81;    ///< unit weight of water (kN/m^3)
    int order = 2;         ///< H1 order of u (1 or 2: PoreField interpolates P2)
    FarSides far;          ///< conditions on the left side, the base and the right side
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

/// Assembles and solves (skyline Cholesky, bandwidth renumbering); the solution is left in cmesh. vtk: optional
/// file of u, -grad u (seepage force) and the Darcy velocity (resolution 1)
inline int64_t SolveSeepage(TPZCompMesh *cmesh, const std::string &vtk = "") {
    TPZLinearAnalysis an(cmesh, true);
    TPZSkylineStructMatrix<STATE> skyl(cmesh);
    an.SetStructuralMatrix(skyl);
    TPZStepSolver<STATE> step;
    step.SetDirect(ECholesky);
    an.SetSolver(step);
    an.Run();
    if (!vtk.empty()) {
        TPZManVector<std::string, 1> scal = {"ExcessPorePressure"};
        TPZManVector<std::string, 2> vec = {"SeepageForce", "DarcyVelocity"};
        an.DefineGraphMesh(2, std::set<int>{ESoil}, scal, vec, vtk);
        an.PostProcess(1);
    }
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

/// Drawdown seepage on gmesh (the geometric mesh is not modified; cmesh is deleted). Far sides with u = 0 or
/// u = -gamma_w h_w are Dirichlet (penalty); where two Dirichlet sides meeting at a box corner carry different data
/// (among the presets only the top right corner of zero_lbr; with per-side overrides also e.g. base zero + right toe)
/// the two penalties average them at the corner node, whereas fe_seepage.py lets the surface win.
inline SeepageResult DrawdownSeepage(TPZGeoMesh *gmesh, const SlopeGeometry &g, const Hydraulics &hy,
                                     const std::string &vtk = "") {
    const auto t0 = std::chrono::steady_clock::now();
    const ScalarBC ud = [g, hy](const TPZVec<REAL> &x) { return DrawdownExcessPressure(g, hy.gammaw, x[1]); };
    std::map<int, ScalarBC> dirichlet = {{EToe, ud}, {ECrest, ud}, {EFace, ud}};
    const std::pair<int, EFarBC> far[3] = {{ELeft, hy.far.left}, {EBase, hy.far.bottom}, {ERight, hy.far.right}};
    for (auto &side : far) {
        if (side.second == EFarBC::ENoFlow) continue;
        const REAL v = side.second == EFarBC::EZero ? 0. : -hy.gammaw * g.hw;
        dirichlet[side.first] = [v](const TPZVec<REAL> &) { return v; };
    }
    TPZCompMesh *cmesh = CreateSeepageCMesh(gmesh, hy, dirichlet, {});
    SeepageResult r;
    r.neq = SolveSeepage(cmesh, vtk);
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
