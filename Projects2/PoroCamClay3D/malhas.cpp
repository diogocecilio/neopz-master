// malhas.cpp — ver malhas.h

#include "malhas.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <map>
#include <vector>

#include "pzgeoel.h"
#include "pzgnode.h"
#include "tpzgeoelrefpattern.h"
#include "tpzquadraticcube.h"
#include "tpzquadraticquad.h"

namespace {

using V3 = std::array<REAL, 3>;
using Cell = std::array<int64_t, 8>;

// arestas do hexaedro na ordem dos nós 8-19 do TPZQuadraticCube
const int kEdges[12][2] = {{0, 1}, {1, 2}, {2, 3}, {3, 0}, {0, 4}, {1, 5},
                           {2, 6}, {3, 7}, {4, 5}, {5, 6}, {6, 7}, {7, 4}};
// faces do hexaedro (vértices locais)
const int kFaces[6][4] = {{0, 3, 2, 1}, {4, 5, 6, 7}, {0, 1, 5, 4}, {1, 2, 6, 5}, {2, 3, 7, 6}, {3, 0, 4, 7}};

/// Monta o TPZGeoMesh a partir dos vértices e das células (8 vértices); com quadratic = true cria os nós de
/// meio de aresta (adjustMid pode movê-los, p.ex. para um arco). As faces que aparecem em uma única célula
/// são faces de contorno, com o marcador dado por classify(vértices da face).
TPZGeoMesh *Build(std::vector<V3> coords, const std::vector<Cell> &cells, bool quadratic,
                  const std::function<void(int64_t, int64_t, V3 &)> &adjustMid,
                  const std::function<int(const std::vector<V3> &)> &classify) {
    std::map<std::pair<int64_t, int64_t>, int64_t> edgeNode;
    auto mid = [&](int64_t a, int64_t b) {
        const auto key = std::make_pair(std::min(a, b), std::max(a, b));
        auto it = edgeNode.find(key);
        if (it != edgeNode.end()) return it->second;
        V3 p;
        for (int i = 0; i < 3; i++) p[i] = 0.5 * (coords[a][i] + coords[b][i]);
        if (adjustMid) adjustMid(a, b, p);
        coords.push_back(p);
        edgeNode[key] = int64_t(coords.size()) - 1;
        return int64_t(coords.size()) - 1;
    };
    std::vector<std::vector<int64_t>> elems;
    for (const Cell &c : cells) {
        std::vector<int64_t> nodes(c.begin(), c.end());
        if (quadratic)
            for (auto &e : kEdges) nodes.push_back(mid(c[e[0]], c[e[1]]));
        elems.push_back(nodes);
    }
    // faces de contorno
    std::map<std::array<int64_t, 4>, std::vector<std::array<int64_t, 4>>> count;
    std::vector<std::array<int64_t, 4>> order;
    for (const Cell &c : cells)
        for (auto &f : kFaces) {
            std::array<int64_t, 4> v = {c[f[0]], c[f[1]], c[f[2]], c[f[3]]}, key = v;
            std::sort(key.begin(), key.end());
            if (!count.count(key)) order.push_back(key);
            count[key].push_back(v);
        }
    std::vector<std::pair<std::vector<int64_t>, int>> faces;
    for (auto &key : order) {
        const auto &lst = count[key];
        if (lst.size() != 1) continue;
        const auto &f = lst[0];
        std::vector<V3> X = {coords[f[0]], coords[f[1]], coords[f[2]], coords[f[3]]};
        std::vector<int64_t> nodes(f.begin(), f.end());
        if (quadratic)
            for (int i = 0; i < 4; i++) {
                const int64_t a = f[i], b = f[(i + 1) % 4];
                nodes.push_back(edgeNode.at(std::make_pair(std::min(a, b), std::max(a, b))));
            }
        faces.push_back({nodes, classify(X)});
    }
    // malha do NeoPZ
    auto *gmesh = new TPZGeoMesh;
    gmesh->SetDimension(3);
    gmesh->NodeVec().Resize(int64_t(coords.size()));
    for (size_t i = 0; i < coords.size(); i++) {
        TPZManVector<REAL, 3> x = {coords[i][0], coords[i][1], coords[i][2]};
        gmesh->NodeVec()[i].Initialize(x, *gmesh);
    }
    int64_t index;
    for (auto &e : elems) {
        TPZManVector<int64_t, 20> nodes(int64_t(e.size()));
        for (size_t i = 0; i < e.size(); i++) nodes[i] = e[i];
        if (quadratic) new TPZGeoElRefPattern<pzgeom::TPZQuadraticCube>(nodes, kMatVolume, *gmesh);
        else gmesh->CreateGeoElement(ECube, nodes, kMatVolume, index);
    }
    for (auto &f : faces) {
        TPZManVector<int64_t, 8> nodes(int64_t(f.first.size()));
        for (size_t i = 0; i < f.first.size(); i++) nodes[i] = f.first[i];
        if (quadratic) new TPZGeoElRefPattern<pzgeom::TPZQuadraticQuad>(nodes, f.second, *gmesh);
        else gmesh->CreateGeoElement(EQuadrilateral, nodes, f.second, index);
    }
    gmesh->BuildConnectivity();
    return gmesh;
}

bool All(const std::vector<V3> &X, const std::function<bool(const V3 &)> &crit) {
    return std::all_of(X.begin(), X.end(), crit);
}

} // namespace

TPZGeoMesh *BoxMesh(const std::array<REAL, 3> &L, const std::array<int, 3> &n, bool quadratic) {
    const int nx = n[0], ny = n[1], nz = n[2];
    auto idx = [&](int i, int j, int k) { return int64_t(k) * (ny + 1) * (nx + 1) + int64_t(j) * (nx + 1) + i; };
    std::vector<V3> coords;
    for (int k = 0; k <= nz; k++)
        for (int j = 0; j <= ny; j++)
            for (int i = 0; i <= nx; i++) coords.push_back({L[0] * i / nx, L[1] * j / ny, L[2] * k / nz});
    std::vector<Cell> cells;
    for (int k = 0; k < nz; k++)
        for (int j = 0; j < ny; j++)
            for (int i = 0; i < nx; i++)
                cells.push_back({idx(i, j, k), idx(i + 1, j, k), idx(i + 1, j + 1, k), idx(i, j + 1, k),
                                 idx(i, j, k + 1), idx(i + 1, j, k + 1), idx(i + 1, j + 1, k + 1), idx(i, j + 1, k + 1)});
    const REAL tol = 1.e-9 * std::max({L[0], L[1], L[2]});
    auto classify = [&](const std::vector<V3> &X) {
        if (All(X, [&](const V3 &p) { return std::fabs(p[0]) < tol; })) return 1;
        if (All(X, [&](const V3 &p) { return std::fabs(p[0] - L[0]) < tol; })) return 2;
        if (All(X, [&](const V3 &p) { return std::fabs(p[1]) < tol; })) return 3;
        if (All(X, [&](const V3 &p) { return std::fabs(p[1] - L[1]) < tol; })) return 4;
        if (All(X, [&](const V3 &p) { return std::fabs(p[2]) < tol; })) return 5;
        return 6;
    };
    return Build(coords, cells, quadratic, nullptr, classify);
}

TPZGeoMesh *QuarterCylinderMesh(REAL r, REAL h, int nc, int nr, int nz, bool quadratic, REAL aratio) {
    // ---- quarto de disco "O-grid" (porte de quarter_disk do Python)
    const REAL a = aratio * r;
    std::map<std::pair<long long, long long>, int64_t> index;
    std::vector<std::array<REAL, 2>> pts;
    auto node = [&](REAL x, REAL y) {
        const auto key = std::make_pair(std::llround(x / r * 1.e9), std::llround(y / r * 1.e9));
        auto it = index.find(key);
        if (it != index.end()) return it->second;
        pts.push_back({x, y});
        index[key] = int64_t(pts.size()) - 1;
        return int64_t(pts.size()) - 1;
    };
    std::vector<std::array<int64_t, 4>> quads;
    std::vector<std::vector<int64_t>> S(nc + 1, std::vector<int64_t>(nc + 1));
    for (int i = 0; i <= nc; i++)
        for (int j = 0; j <= nc; j++) S[i][j] = node(a * i / nc, a * j / nc);
    for (int i = 0; i < nc; i++)
        for (int j = 0; j < nc; j++) quads.push_back({S[i][j], S[i + 1][j], S[i + 1][j + 1], S[i][j + 1]});
    for (int block = 0; block < 2; block++) {
        std::vector<std::vector<int64_t>> G;
        for (int j = 0; j <= nc; j++) {
            REAL inner[2], th;
            if (block == 0) {  // leste
                inner[0] = a; inner[1] = a * j / nc; th = 0.25 * M_PI * j / nc;
            } else {           // norte
                inner[0] = a * j / nc; inner[1] = a; th = 0.5 * M_PI - 0.25 * M_PI * j / nc;
            }
            const REAL outer[2] = {r * std::cos(th), r * std::sin(th)};
            std::vector<int64_t> row;
            for (int k = 0; k <= nr; k++)
                row.push_back(node(inner[0] + REAL(k) / nr * (outer[0] - inner[0]),
                                   inner[1] + REAL(k) / nr * (outer[1] - inner[1])));
            G.push_back(row);
        }
        for (int j = 0; j < nc; j++)
            for (int k = 0; k < nr; k++) quads.push_back({G[j][k], G[j][k + 1], G[j + 1][k + 1], G[j + 1][k]});
    }
    for (auto &q : quads) {  // orientação anti-horária (normal +z)
        REAL area = 0.;
        for (int i = 0; i < 4; i++) {
            const auto &p = pts[q[i]], &pn = pts[q[(i + 1) % 4]];
            area += 0.5 * (p[0] * pn[1] - pn[0] * p[1]);
        }
        if (area < 0.) std::swap(q[1], q[3]);
    }
    // ---- extrusão em z
    const int64_t n2 = int64_t(pts.size());
    std::vector<V3> coords;
    for (int k = 0; k <= nz; k++)
        for (auto &p : pts) coords.push_back({p[0], p[1], h * k / nz});
    std::vector<Cell> cells;
    for (int k = 0; k < nz; k++)
        for (auto &q : quads)
            cells.push_back({k * n2 + q[0], k * n2 + q[1], k * n2 + q[2], k * n2 + q[3],
                             (k + 1) * n2 + q[0], (k + 1) * n2 + q[1], (k + 1) * n2 + q[2], (k + 1) * n2 + q[3]});
    const REAL tol = 1.e-9 * r;
    const std::vector<V3> &cref = coords;
    auto onLateral = [&](int64_t i) { return std::fabs(std::hypot(cref[i][0], cref[i][1]) - r) < tol; };
    auto adjust = [&](int64_t i, int64_t j, V3 &p) {  // aresta na superfície lateral: nó no arco
        if (onLateral(i) && onLateral(j)) {
            const REAL s = r / std::hypot(p[0], p[1]);
            p[0] *= s;
            p[1] *= s;
        }
    };
    auto classify = [&](const std::vector<V3> &X) {
        if (All(X, [&](const V3 &p) { return std::fabs(p[0]) < tol; })) return int(EX0);
        if (All(X, [&](const V3 &p) { return std::fabs(p[1]) < tol; })) return int(EY0);
        if (All(X, [&](const V3 &p) { return std::fabs(p[2]) < tol; })) return int(EBase);
        if (All(X, [&](const V3 &p) { return std::fabs(p[2] - h) < tol; })) return int(ETopo);
        return int(ELateral);
    };
    return Build(coords, cells, quadratic, adjust, classify);
}
