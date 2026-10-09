// Parametric slope of Ceron et al. (2025) in NeoPZ coordinates (y up): crest edge O = (0, 0), crest y = 0 for
// x <= 0, face O -> T at beta degrees, toe T = (H / tan(beta), -H), toe ground y = -H, water-line point
// W = (h_w / tan(beta), -h_w) on the face (level after the drawdown). Box: 'left' to the left of O, 'right' to the
// right of T, 'depth' below T. Paper coordinates (y down): x_paper = x, y_paper = -y.
//
// Geometric meshes of linear triangles (the field evaluators of SeepageForceField.h assume straight triangles):
//  - CreateSlopeGMesh: Delaunay refinement (DelaunayMesher.h) with a size function graded towards O, W and T
//    (h0) and the face (hs), growing by 'grade' per unit distance up to hmax; nodes at O, W and T; optional
//    uniform refinement. Used for the large hydraulic domain (paper Fig. 4) and for the stability domain.
//  - SlopeMohrCoulombGMesh: TriGMesh of Projects/SlopeMohrCoulomb (70 x 40 m, H = 10 m, beta = 45 deg) moved to
//    the convention above, for the regression against Projects/SlopeDrawdown.
// Boundary ids (as in TriGMesh): -1 base, -2 right side, -3 toe ground, -4 crest, -5 left side, -6 face; soil 1.
#ifndef SLOPEGEOMETRY_H
#define SLOPEGEOMETRY_H

#include "../SlopeMohrCoulomb/SlopeModel.h"
#include "DelaunayMesher.h"

#include "TPZGeoLinear.h"
#include "pzgeotriangle.h"
#include "pzgmesh.h"
#include "tpzgeoelrefpattern.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <map>
#include <vector>

namespace slope {

enum EBoundary { EBase = -1, ERight = -2, EToe = -3, ECrest = -4, ELeft = -5, EFace = -6 };
constexpr int ESoil = 1;

struct SlopeGeometry {
    REAL H = 10.;    ///< slope height (m)
    REAL beta = 45.; ///< face angle (deg), 15 <= beta <= 90
    REAL hw = 10.;   ///< water level after the drawdown, below the crest (m), 0 <= hw <= H
    REAL left = 30., right = 30., depth = 30.; ///< box extents (m): left of O, right of T, below T

    REAL Cot() const { return beta >= 90. ? 0. : 1. / std::tan(beta * M_PI / 180.); }
    REAL XT() const { return H * Cot(); }
    Pt O() const { return {0., 0.}; }
    Pt T() const { return {XT(), -H}; }
    Pt W() const { return {hw * Cot(), -hw}; }

    /// counterclockwise boundary polygon and the boundary id of each edge (i -> i + 1)
    std::vector<Pt> Polygon(std::vector<int> &marker) const {
        const REAL xT = XT();
        std::vector<Pt> p = {{-left, -H - depth}, {xT + right, -H - depth}, {xT + right, -H}, T()};
        marker = {EBase, ERight, EToe, EFace};
        if (hw > 1.e-9 * H && hw < H * (1. - 1.e-9)) {
            p.push_back(W());
            marker.push_back(EFace);
        }
        p.push_back(O());
        marker.push_back(ECrest);
        p.push_back({-left, 0.});
        marker.push_back(ELeft);
        return p;
    }

    void Check() const {
        if (H <= 0. || beta < 1. || beta > 90. || hw < 0. || hw > H || left <= 0. || right <= 0. || depth <= 0.) {
            std::cerr << "SlopeGeometry: invalid data (H " << H << ", beta " << beta << ", hw " << hw << ")\n";
            DebugStop();
        }
    }
};

/// Box of the stability analysis: a H + H / tan(beta) on each side (left of O, right of T, below T). With a = 2 it
/// is the 70 x 40 m box of SlopeMohrCoulomb for H = 10 m, beta = 45 deg; the extra H / tan(beta) keeps the box
/// proportional to the length of the face, which sets the size of the failure mechanism of flat slopes.
inline void StabilityExtents(SlopeGeometry &g, REAL a) {
    g.left = g.right = g.depth = a * g.H + g.XT();
}

/// Target element size (in units of H): h0 at O, W and T, hs along the face, + grade * distance, at most hmax
struct MeshSize {
    REAL h0 = 1. / 40., hs = 1. / 16., grade = 0.15, hmax = 2.;
    REAL minAngle = 28.; ///< degrees (Delaunay refinement)

    REAL operator()(const SlopeGeometry &g, REAL x, REAL y) const {
        const Pt O = g.O(), T = g.T(), W = g.W();
        REAL h = hmax;
        for (const Pt &a : {O, T, W}) h = std::min(h, h0 + grade * std::hypot(x - a.x, y - a.y) / g.H);
        const REAL dx = T.x - O.x, dy = T.y - O.y;
        const REAL t = std::clamp(((x - O.x) * dx + (y - O.y) * dy) / (dx * dx + dy * dy), 0., 1.);
        h = std::min(h, hs + grade * std::hypot(x - O.x - t * dx, y - O.y - t * dy) / g.H);
        return h * g.H;
    }
};

/// Uniform refinement of all the leaf elements (volume and boundary)
inline void UniformRefine(TPZGeoMesh *gmesh, int nref) {
    for (int d = 0; d < nref; d++) {
        const int64_t nel = gmesh->NElements();
        TPZManVector<TPZGeoEl *> sub;
        for (int64_t i = 0; i < nel; i++) {
            TPZGeoEl *gel = gmesh->Element(i);
            if (gel && !gel->HasSubElement()) gel->Divide(sub);
        }
    }
}

struct MeshStats {
    int64_t nodes = 0, triangles = 0;
    REAL minAngle = 180., hmin = 1.e300, hmax = 0.;
};

/// Linear triangles and boundary lines of a polygon meshed by DelaunayMesher
inline TPZGeoMesh *CreateSlopeGMesh(const SlopeGeometry &g, const MeshSize &ms, int uniformRef = 0,
                                    MeshStats *stats = nullptr) {
    g.Check();
    std::vector<int> marker;
    const std::vector<Pt> poly = g.Polygon(marker);
    auto size = [&g, &ms](double x, double y) { return ms(g, x, y); };
    const DelaunayMesher::Result r = DelaunayMesher::Mesh(poly, marker, size, ms.minAngle);
    auto *gmesh = new TPZGeoMesh();
    gmesh->SetDimension(2);
    gmesh->NodeVec().Resize(r.nodes.size());
    for (size_t i = 0; i < r.nodes.size(); i++) {
        TPZManVector<REAL, 3> x = {r.nodes[i].x, r.nodes[i].y, 0.};
        gmesh->NodeVec()[i] = TPZGeoNode(i, x, *gmesh);
    }
    for (auto &t : r.triangles) {
        TPZManVector<int64_t, 3> nodes = {t[0], t[1], t[2]};
        new TPZGeoElRefPattern<pzgeom::TPZGeoTriangle>(nodes, ESoil, *gmesh);
    }
    for (auto &e : r.edges) {
        TPZManVector<int64_t, 2> nodes = {e[0], e[1]};
        new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodes, e[2], *gmesh);
    }
    gmesh->BuildConnectivity();
    UniformRefine(gmesh, uniformRef);
    if (stats) {
        *stats = MeshStats();
        stats->nodes = r.nodes.size();
        stats->triangles = r.triangles.size();
        for (auto &t : r.triangles)
            for (int i = 0; i < 3; i++) {
                const Pt &p = r.nodes[t[i]], &q = r.nodes[t[(i + 1) % 3]], &s = r.nodes[t[(i + 2) % 3]];
                const REAL l = std::hypot(q.x - p.x, q.y - p.y);
                const REAL c = ((q.x - p.x) * (s.x - p.x) + (q.y - p.y) * (s.y - p.y)) / (l * std::hypot(s.x - p.x, s.y - p.y));
                stats->minAngle = std::min<REAL>(stats->minAngle, std::acos(std::clamp<REAL>(c, -1., 1.)) * 180. / M_PI);
                stats->hmin = std::min(stats->hmin, l), stats->hmax = std::max(stats->hmax, l);
            }
    }
    return gmesh;
}

/// Consistency of a slope mesh with its polygon: sum of the leaf triangle areas and total length of the leaf boundary
/// lines of each id, relative to the exact values (round-off when the mesh is valid); every leaf triangle must be
/// counterclockwise with a positive area, and the mesh must be conforming: every triangle edge (pair of node
/// indices) is shared by two triangles or by one triangle and one boundary line (otherwise the error is 1)
inline REAL CheckSlopeGMesh(TPZGeoMesh *gmesh, const SlopeGeometry &g, bool verbose = false) {
    std::vector<int> marker;
    const std::vector<Pt> poly = g.Polygon(marker);
    REAL area = 0.;
    std::map<int, REAL> len, lenex;
    for (size_t i = 0; i < poly.size(); i++) {
        const Pt &a = poly[i], &b = poly[(i + 1) % poly.size()];
        area += 0.5 * (a.x * b.y - b.x * a.y);
        lenex[marker[i]] += std::hypot(b.x - a.x, b.y - a.y);
    }
    REAL sum = 0., minArea = 1.e300;
    std::map<std::pair<int64_t, int64_t>, std::array<int, 2>> edges; ///< node pair -> (triangles, boundary lines)
    auto edge = [&edges](int64_t a, int64_t b) -> std::array<int, 2> & { return edges[{std::min(a, b), std::max(a, b)}]; };
    for (int64_t i = 0; i < gmesh->NElements(); i++) {
        TPZGeoEl *gel = gmesh->Element(i);
        if (!gel || gel->HasSubElement()) continue;
        REAL x[3][2];
        for (int k = 0; k < gel->NCornerNodes(); k++)
            for (int d = 0; d < 2; d++) x[k][d] = gel->NodePtr(k)->Coord(d);
        if (gel->Dimension() == 2) {
            const REAL a = 0.5 * ((x[1][0] - x[0][0]) * (x[2][1] - x[0][1]) - (x[1][1] - x[0][1]) * (x[2][0] - x[0][0]));
            sum += a;
            minArea = std::min(minArea, a);
            for (int k = 0; k < 3; k++) edge(gel->NodeIndex(k), gel->NodeIndex((k + 1) % 3))[0]++;
        } else if (gel->Dimension() == 1) {
            len[gel->MaterialId()] += std::hypot(x[1][0] - x[0][0], x[1][1] - x[0][1]);
            edge(gel->NodeIndex(0), gel->NodeIndex(1))[1]++;
        }
    }
    int64_t nonconforming = 0;
    for (auto &e : edges)
        if (!((e.second[0] == 2 && e.second[1] == 0) || (e.second[0] == 1 && e.second[1] == 1))) nonconforming++;
    REAL err = std::fabs(sum - area) / area;
    for (auto &l : lenex) err = std::max(err, std::fabs(len[l.first] - l.second) / l.second);
    if (minArea <= 0. || nonconforming > 0) err = std::max<REAL>(err, 1.);
    if (verbose)
        std::cout << "  mesh check: area " << sum << " (exact " << area << "), min triangle area " << minArea
                  << ", non-conforming edges " << nonconforming
                  << ", max relative error of the area and of the boundary lengths " << err << "\n";
    return err;
}

/// TriGMesh(ref) of SlopeMohrCoulomb (base y = 0, crest y = 40, O = (30, 40)) translated to O = (0, 0): the slope of
/// SlopeDrawdown, H = 10 m, beta = 45 deg, 30 m left of O, right of T and below T
inline TPZGeoMesh *SlopeMohrCoulombGMesh(int ref) {
    TPZGeoMesh *gmesh = TriGMesh(0);
    for (int64_t i = 0; i < gmesh->NNodes(); i++) {
        TPZGeoNode &n = gmesh->NodeVec()[i];
        n.SetCoord(0, n.Coord(0) - 30.);
        n.SetCoord(1, n.Coord(1) - 40.);
    }
    UniformRefine(gmesh, ref); // after the translation: the sub-element nodes are created at the new coordinates
    return gmesh;
}

inline SlopeGeometry SlopeMohrCoulombGeometry(REAL hw) {
    SlopeGeometry g;
    g.H = 10., g.beta = 45., g.hw = hw, g.left = g.right = g.depth = 30.;
    return g;
}

} // namespace slope

#endif
