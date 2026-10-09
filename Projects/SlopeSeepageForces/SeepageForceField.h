// Seepage force field f = -grad u (kN/m^3, NeoPZ coordinates, y up) at any point of the soil, for the FEM stability
// (forcing function evaluated at every integration point of every assembly, from several threads) and for the
// limit analysis (to be ported from scripts/limit_analysis.py).
//
//  - ForceField: generic interface; FE fields (PoreField below), the analytical K^-1 v'_opt field
//    (AnalyticalSeepage.h, port of scripts/analytical_seepage.py) and the dry/no-seepage case all plug in through it.
//  - PoreField: generalization of PoreField of Projects/SlopeDrawdown/main.cpp (P1 pressure, linear search) to the
//    excess pore pressure u of the seepage FE solution: per straight triangle, the quadratic Lagrange interpolant of
//    the FE solution at the 3 vertices and 3 edge midpoints (exact for order <= 2, so grad u is the exact P2
//    gradient), and point location with a quadtree of the triangle bounding boxes. Read-only after construction
//    (thread safe).
#ifndef SEEPAGEFORCEFIELD_H
#define SEEPAGEFORCEFIELD_H

#include "pzcmesh.h"
#include "pzgeoel.h"
#include "pzinterpolationspace.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <memory>
#include <vector>

namespace slope {

/// f(x) = seepage force density (kN/m^3) at x = (x, y, z) in NeoPZ coordinates; zero outside the soil
using ForceField = std::function<void(const TPZVec<REAL> &x, REAL f[2])>;

inline ForceField NoSeepage() {
    return [](const TPZVec<REAL> &, REAL f[2]) { f[0] = f[1] = 0.; };
}

/// Adapter for a field written in the paper coordinates (x right, y down), as the Python reference fields
/// (field.force(x, y) of scripts/): fpaper(x, y_paper, f_paper) -> ForceField in NeoPZ coordinates
inline ForceField FromPaperCoordinates(const std::function<void(REAL x, REAL yp, REAL fp[2])> &fpaper) {
    return [fpaper](const TPZVec<REAL> &x, REAL f[2]) {
        REAL fp[2];
        fpaper(x[0], -x[1], fp);
        f[0] = fp[0], f[1] = -fp[1];
    };
}

class PoreField {
public:
    struct Tri {
        REAL x0[2], inv[4]; ///< x = x0 + J xi, inv = J^-1 (row major)
        REAL u[6];          ///< u at the vertices 0, 1, 2 and at the midpoints of the edges 01, 12, 20
        REAL area;
    };

    PoreField() = default;

    /// From the solution of a 2D H1 mesh (material matid, order <= 2, straight triangles); var = index of the
    /// scalar solution variable of the material
    PoreField(TPZCompMesh *cmesh, int matid, int var) {
        static const REAL qsi[6][2] = {{0., 0.}, {1., 0.}, {0., 1.}, {0.5, 0.}, {0.5, 0.5}, {0., 0.5}};
        for (int64_t el = 0; el < cmesh->NElements(); el++) {
            auto *intel = dynamic_cast<TPZInterpolationSpace *>(cmesh->ElementVec()[el]);
            if (!intel || !intel->Reference() || intel->Reference()->MaterialId() != matid) continue;
            TPZGeoEl *gel = intel->Reference();
            if (gel->Type() != ETriangle || !gel->IsLinearMapping() || intel->MaxOrder() > 2) DebugStop();
            REAL x[3][2];
            for (int k = 0; k < 3; k++)
                for (int d = 0; d < 2; d++) x[k][d] = gel->NodePtr(k)->Coord(d);
            Tri t;
            SetGeometry(t, x);
            for (int k = 0; k < 6; k++) {
                TPZManVector<REAL, 3> xi = {qsi[k][0], qsi[k][1], 0.};
                TPZManVector<STATE, 3> sol(1, 0.);
                intel->Solution(xi, var, sol);
                t.u[k] = sol[0];
            }
            fTri.push_back(t);
        }
        BuildTree();
    }

    /// From nodal data: triangles as vertex coordinates and the 6 values (vertices, then edge midpoints)
    PoreField(const std::vector<std::array<REAL, 6>> &xy, const std::vector<std::array<REAL, 6>> &u) {
        for (size_t i = 0; i < xy.size(); i++) {
            const REAL x[3][2] = {{xy[i][0], xy[i][1]}, {xy[i][2], xy[i][3]}, {xy[i][4], xy[i][5]}};
            Tri t;
            SetGeometry(t, x);
            std::copy(u[i].begin(), u[i].end(), t.u);
            fTri.push_back(t);
        }
        BuildTree();
    }

    bool Empty() const { return fTri.empty(); }
    const std::vector<Tri> &Triangles() const { return fTri; }

    /// Quadratic interpolant on triangle t at the parametric point (xi, eta): u and grad u (x, y)
    static void EvaluateTri(const Tri &t, REAL xi, REAL eta, REAL &u, REAL grad[2]) {
        const REAL l0 = 1. - xi - eta, l1 = xi, l2 = eta;
        const REAL N[6] = {l0 * (2. * l0 - 1.), l1 * (2. * l1 - 1.), l2 * (2. * l2 - 1.),
                           4. * l0 * l1,        4. * l1 * l2,        4. * l2 * l0};
        // derivatives with respect to xi and eta (dl0 = (-1, -1), dl1 = (1, 0), dl2 = (0, 1))
        const REAL Nx[6] = {-(4. * l0 - 1.), 4. * l1 - 1., 0., 4. * (l0 - l1), 4. * l2, -4. * l2};
        const REAL Ne[6] = {-(4. * l0 - 1.), 0., 4. * l2 - 1., -4. * l1, 4. * l1, 4. * (l0 - l2)};
        REAL ux = 0., ue = 0.;
        u = 0.;
        for (int k = 0; k < 6; k++) {
            u += t.u[k] * N[k];
            ux += t.u[k] * Nx[k];
            ue += t.u[k] * Ne[k];
        }
        grad[0] = t.inv[0] * ux + t.inv[2] * ue; // J^-T grad_xi
        grad[1] = t.inv[1] * ux + t.inv[3] * ue;
    }

    /// u and grad u at x (false outside the mesh)
    bool Evaluate(const TPZVec<REAL> &x, REAL &u, REAL grad[2]) const {
        REAL xi, eta;
        const int i = Locate(x[0], x[1], xi, eta);
        if (i < 0) return false;
        EvaluateTri(fTri[i], xi, eta, u, grad);
        return true;
    }

    /// seepage force f = -grad u, zero outside the mesh
    void Force(const TPZVec<REAL> &x, REAL f[2]) const {
        REAL u, g[2];
        if (!Evaluate(x, u, g)) g[0] = g[1] = 0.;
        f[0] = -g[0], f[1] = -g[1];
    }

    ForceField AsForceField(std::shared_ptr<const PoreField> self) const {
        return [self](const TPZVec<REAL> &x, REAL f[2]) { self->Force(x, f); };
    }

    /// index of the triangle that contains (x, y) and its parametric coordinates (-1 outside)
    int Locate(REAL x, REAL y, REAL &xi, REAL &eta) const {
        if (fNodes.empty()) return -1;
        const Node *n = &fNodes[0];
        if (x < n->box[0] || x > n->box[2] || y < n->box[1] || y > n->box[3]) return -1;
        while (n->child >= 0) {
            const REAL xm = 0.5 * (n->box[0] + n->box[2]), ym = 0.5 * (n->box[1] + n->box[3]);
            n = &fNodes[n->child + (x > xm ? 1 : 0) + (y > ym ? 2 : 0)];
        }
        int best = -1;
        REAL bestm = -1.e300, bxi = 0., beta = 0.;
        for (int k = n->first; k < n->first + n->count; k++) {
            const Tri &t = fTri[fItems[k]];
            const REAL dx = x - t.x0[0], dy = y - t.x0[1];
            const REAL a = t.inv[0] * dx + t.inv[1] * dy, b = t.inv[2] * dx + t.inv[3] * dy;
            const REAL m = std::min({a, b, 1. - a - b});
            if (m >= 0.) {
                xi = a, eta = b;
                return fItems[k];
            }
            if (m > bestm) bestm = m, best = fItems[k], bxi = a, beta = b;
        }
        if (best >= 0 && bestm > -1.e-9) { // on an edge, within rounding
            xi = bxi, eta = beta;
            return best;
        }
        return -1;
    }

private:
    struct Node {
        REAL box[4]; ///< xmin, ymin, xmax, ymax
        int child = -1, first = 0, count = 0;
    };
    std::vector<Tri> fTri;
    std::vector<Node> fNodes;
    std::vector<int> fItems; ///< triangle indices of the leaves (contiguous per leaf)

    static void SetGeometry(Tri &t, const REAL x[3][2]) {
        const REAL a = x[1][0] - x[0][0], b = x[2][0] - x[0][0], c = x[1][1] - x[0][1], d = x[2][1] - x[0][1];
        const REAL det = a * d - b * c;
        if (!(std::fabs(det) > 0.)) DebugStop();
        t.x0[0] = x[0][0], t.x0[1] = x[0][1];
        t.inv[0] = d / det, t.inv[1] = -b / det, t.inv[2] = -c / det, t.inv[3] = a / det;
        t.area = 0.5 * std::fabs(det);
    }

    void TriBox(const Tri &t, REAL box[4]) const {
        // vertices from x0 and J = inv^-1
        const REAL det = t.inv[0] * t.inv[3] - t.inv[1] * t.inv[2];
        const REAL J[4] = {t.inv[3] / det, -t.inv[1] / det, -t.inv[2] / det, t.inv[0] / det};
        const REAL xs[3] = {t.x0[0], t.x0[0] + J[0], t.x0[0] + J[1]};
        const REAL ys[3] = {t.x0[1], t.x0[1] + J[2], t.x0[1] + J[3]};
        box[0] = std::min({xs[0], xs[1], xs[2]}), box[2] = std::max({xs[0], xs[1], xs[2]});
        box[1] = std::min({ys[0], ys[1], ys[2]}), box[3] = std::max({ys[0], ys[1], ys[2]});
    }

    /// quadtree: a leaf holds the triangles whose bounding box overlaps it; split while it has more than 8
    void BuildTree() {
        fNodes.clear();
        fItems.clear();
        if (fTri.empty()) return;
        std::vector<std::array<REAL, 4>> boxes(fTri.size());
        Node root;
        root.box[0] = root.box[1] = 1.e300, root.box[2] = root.box[3] = -1.e300;
        for (size_t i = 0; i < fTri.size(); i++) {
            TriBox(fTri[i], boxes[i].data());
            root.box[0] = std::min(root.box[0], boxes[i][0]), root.box[1] = std::min(root.box[1], boxes[i][1]);
            root.box[2] = std::max(root.box[2], boxes[i][2]), root.box[3] = std::max(root.box[3], boxes[i][3]);
        }
        const REAL pad = 1.e-9 * std::max(root.box[2] - root.box[0], root.box[3] - root.box[1]);
        root.box[0] -= pad, root.box[1] -= pad, root.box[2] += pad, root.box[3] += pad;
        std::vector<int> all(fTri.size());
        for (size_t i = 0; i < all.size(); i++) all[i] = int(i);
        fNodes.push_back(root);
        Build(0, all, boxes, 0);
    }

    void Build(int node, const std::vector<int> &items, const std::vector<std::array<REAL, 4>> &boxes, int depth) {
        if (items.size() <= 8 || depth >= 24) {
            fNodes[node].first = int(fItems.size());
            fNodes[node].count = int(items.size());
            fItems.insert(fItems.end(), items.begin(), items.end());
            return;
        }
        const Node n = fNodes[node];
        const REAL xm = 0.5 * (n.box[0] + n.box[2]), ym = 0.5 * (n.box[1] + n.box[3]);
        const int child = int(fNodes.size());
        fNodes[node].child = child;
        for (int c = 0; c < 4; c++) {
            Node s;
            s.box[0] = (c & 1) ? xm : n.box[0], s.box[2] = (c & 1) ? n.box[2] : xm;
            s.box[1] = (c & 2) ? ym : n.box[1], s.box[3] = (c & 2) ? n.box[3] : ym;
            fNodes.push_back(s);
        }
        for (int c = 0; c < 4; c++) {
            std::vector<int> sub;
            const Node s = fNodes[child + c];
            for (int i : items)
                if (boxes[i][0] <= s.box[2] && boxes[i][2] >= s.box[0] && boxes[i][1] <= s.box[3] && boxes[i][3] >= s.box[1])
                    sub.push_back(i);
            Build(child + c, sub, boxes, depth + 1);
        }
    }
};

} // namespace slope

#endif
