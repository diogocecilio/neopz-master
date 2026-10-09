// Quality triangular meshes of a simple polygon with a size function: conforming Delaunay triangulation and
// Delaunay refinement (Ruppert / Shewchuk). The boundary is first divided according to the size function; the
// subsegments are kept as edges (a circumcenter that encroaches a subsegment, i.e. lies inside its diametral circle,
// splits the subsegment instead), so every triangle lies entirely inside or outside the polygon. Triangles inside
// the polygon are refined while their smallest angle is below minAngle or their size exceeds size(centroid).
// Self-contained (no NeoPZ dependency); used by SlopeGeometry.h. The input polygon must have no angle below 60 deg
// (the slope domains have 90 deg, 180 - beta and 180 + beta).
#ifndef SLOPE_DELAUNAYMESHER_H
#define SLOPE_DELAUNAYMESHER_H

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <deque>
#include <functional>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace slope {

struct Pt {
    double x, y;
};

class DelaunayMesher {
public:
    struct Result {
        std::vector<Pt> nodes;
        std::vector<std::array<int, 3>> triangles; ///< counterclockwise
        std::vector<std::array<int, 3>> edges;     ///< boundary edges {a, b, marker}, oriented as the polygon
    };

    /// polygon: counterclockwise; marker[i] belongs to the edge polygon[i] -> polygon[i + 1];
    /// size(x, y): target edge length; minAngle in degrees (Ruppert terminates up to ~ 30)
    static Result Mesh(const std::vector<Pt> &polygon, const std::vector<int> &marker,
                       const std::function<double(double, double)> &size, double minAngle = 28.,
                       int64_t maxNodes = 4000000) {
        DelaunayMesher m(polygon, marker, size, minAngle, maxNodes);
        m.Run();
        return m.Extract();
    }

private:
    struct Tri {
        int v[3]; ///< counterclockwise vertices
        int n[3]; ///< n[i]: neighbour across the edge opposite to v[i] (-1: none)
        bool alive, inside;
    };
    struct Seg {
        int a, b, marker; ///< oriented as the polygon
    };

    std::vector<Pt> fPoly;
    std::vector<int> fMarker;
    std::function<double(double, double)> fSize;
    double fRatio; ///< circumradius / shortest edge limit = 1 / (2 sin(minAngle))
    int64_t fMaxNodes;
    std::vector<Pt> fP;
    std::vector<Tri> fT;
    std::vector<int> fVTri; ///< one triangle incident to each vertex
    std::unordered_map<uint64_t, Seg> fSeg;
    std::vector<int> fStamp;
    int fCurStamp = 0, fHint = 0;
    uint32_t fRand = 12345u;
    std::deque<int> fQueue;

    DelaunayMesher(const std::vector<Pt> &polygon, const std::vector<int> &marker,
                   const std::function<double(double, double)> &size, double minAngle, int64_t maxNodes)
        : fPoly(polygon), fMarker(marker), fSize(size), fMaxNodes(maxNodes) {
        if (polygon.size() < 3 || marker.size() != polygon.size()) throw std::invalid_argument("DelaunayMesher: polygon");
        fRatio = 1. / (2. * std::sin(minAngle * M_PI / 180.));
    }

    static uint64_t Key(int a, int b) {
        if (a > b) std::swap(a, b);
        return (uint64_t(uint32_t(a)) << 32) | uint32_t(b);
    }
    static double Orient(const Pt &a, const Pt &b, const Pt &c) {
        return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
    }
    /// > 0 when d lies inside the circumcircle of the counterclockwise triangle abc
    static double InCircle(const Pt &a, const Pt &b, const Pt &c, const Pt &d) {
        const double adx = a.x - d.x, ady = a.y - d.y, bdx = b.x - d.x, bdy = b.y - d.y, cdx = c.x - d.x,
                     cdy = c.y - d.y;
        const double ad = adx * adx + ady * ady, bd = bdx * bdx + bdy * bdy, cd = cdx * cdx + cdy * cdy;
        return adx * (bdy * cd - bd * cdy) - ady * (bdx * cd - bd * cdx) + ad * (bdx * cdy - bdy * cdx);
    }
    static Pt Circumcenter(const Pt &a, const Pt &b, const Pt &c) {
        const double bx = b.x - a.x, by = b.y - a.y, cx = c.x - a.x, cy = c.y - a.y;
        const double d = 2. * (bx * cy - by * cx);
        const double b2 = bx * bx + by * by, c2 = cx * cx + cy * cy;
        return {a.x + (cy * b2 - by * c2) / d, a.y + (bx * c2 - cx * b2) / d};
    }
    bool InsidePolygon(const Pt &p) const {
        bool in = false;
        const size_t n = fPoly.size();
        for (size_t i = 0, j = n - 1; i < n; j = i++) {
            const Pt &a = fPoly[i], &b = fPoly[j];
            if ((a.y > p.y) != (b.y > p.y) && p.x < (b.x - a.x) * (p.y - a.y) / (b.y - a.y) + a.x) in = !in;
        }
        return in;
    }
    uint32_t Rand() {
        fRand = fRand * 1664525u + 1013904223u;
        return fRand >> 8;
    }

    int NewTri(int a, int b, int c) {
        Tri t;
        t.v[0] = a, t.v[1] = b, t.v[2] = c;
        t.n[0] = t.n[1] = t.n[2] = -1;
        t.alive = true;
        const Pt g = {(fP[a].x + fP[b].x + fP[c].x) / 3., (fP[a].y + fP[b].y + fP[c].y) / 3.};
        t.inside = InsidePolygon(g);
        fT.push_back(t);
        fStamp.push_back(0);
        const int id = int(fT.size()) - 1;
        fVTri[a] = fVTri[b] = fVTri[c] = id;
        return id;
    }

    /// triangle containing p (stochastic visibility walk; brute force if the walk fails)
    int Locate(const Pt &p) {
        int t = (fHint >= 0 && fHint < int(fT.size()) && fT[fHint].alive) ? fHint : -1;
        if (t < 0)
            for (int i = int(fT.size()) - 1; i >= 0; i--)
                if (fT[i].alive) { t = i; break; }
        for (size_t step = 0; step < 4 * fT.size() + 100; step++) {
            const Tri &tr = fT[t];
            const int s = int(Rand() % 3);
            int next = -1;
            for (int k = 0; k < 3; k++) {
                const int i = (s + k) % 3;
                if (Orient(fP[tr.v[(i + 1) % 3]], fP[tr.v[(i + 2) % 3]], p) < 0.) { next = tr.n[i]; break; }
            }
            if (next == -1) {
                bool in = true;
                for (int i = 0; i < 3; i++) in = in && Orient(fP[tr.v[(i + 1) % 3]], fP[tr.v[(i + 2) % 3]], p) >= 0.;
                if (in) return t;
                break; // outside the hull (cannot happen inside the bounding box)
            }
            t = next;
        }
        int best = -1;
        double bestv = -1.e300;
        for (int i = 0; i < int(fT.size()); i++) { // fallback: least violated triangle
            if (!fT[i].alive) continue;
            double m = 1.e300;
            for (int k = 0; k < 3; k++) m = std::min(m, Orient(fP[fT[i].v[(k + 1) % 3]], fP[fT[i].v[(k + 2) % 3]], p));
            if (m > bestv) { bestv = m; best = i; }
        }
        return best;
    }

    /// Bowyer-Watson insertion of the vertex iv; seeds: triangles that must be in the cavity; the subsegment
    /// 'cross' (key, 0 = none) may be crossed by the cavity (it is being split at iv)
    void Insert(int iv, const std::vector<int> &seeds, uint64_t cross) {
        const Pt p = fP[iv];
        std::vector<int> cav;
        auto grow = [&]() {
            ++fCurStamp;
            cav.clear();
            for (int s : seeds) { fStamp[s] = fCurStamp; cav.push_back(s); }
            for (size_t k = 0; k < cav.size(); k++) {
                const Tri &t = fT[cav[k]];
                for (int i = 0; i < 3; i++) {
                    const int nb = t.n[i];
                    if (nb < 0 || fStamp[nb] == fCurStamp || fStamp[nb] == -fCurStamp) continue;
                    const uint64_t key = Key(t.v[(i + 1) % 3], t.v[(i + 2) % 3]);
                    if (key != cross && fSeg.count(key)) continue;
                    const Tri &u = fT[nb];
                    if (InCircle(fP[u.v[0]], fP[u.v[1]], fP[u.v[2]], p) > 0.) {
                        fStamp[nb] = fCurStamp;
                        cav.push_back(nb);
                    }
                }
            }
        };
        grow();
        // repair: the cavity must be star-shaped from p and keep all its vertices on its boundary
        for (int iter = 0;; iter++) {
            if (iter > 1000) throw std::runtime_error("DelaunayMesher: cavity repair failed");
            int bad = -1;
            std::unordered_map<int, int> onBoundary;
            for (int c : cav) {
                const Tri &t = fT[c];
                for (int i = 0; i < 3; i++) {
                    const int nb = t.n[i];
                    if (nb >= 0 && fStamp[nb] == fCurStamp) continue;
                    const int a = t.v[(i + 1) % 3], b = t.v[(i + 2) % 3];
                    onBoundary[a]++;
                    if (Orient(fP[a], fP[b], p) <= 0. && std::find(seeds.begin(), seeds.end(), c) == seeds.end())
                        bad = c;
                }
            }
            if (bad < 0)
                for (int c : cav) {
                    if (std::find(seeds.begin(), seeds.end(), c) != seeds.end()) continue;
                    for (int i = 0; i < 3 && bad < 0; i++)
                        if (!onBoundary.count(fT[c].v[i]) || onBoundary[fT[c].v[i]] != 1) bad = c;
                    if (bad >= 0) break;
                }
            if (bad < 0) break;
            // drop 'bad' and keep the part of the cavity connected to the seeds
            fStamp[bad] = -fCurStamp;
            std::vector<int> keep;
            const int st = fCurStamp;
            ++fCurStamp;
            for (int s : seeds) { fStamp[s] = fCurStamp; keep.push_back(s); }
            for (size_t k = 0; k < keep.size(); k++) {
                const Tri &t = fT[keep[k]];
                for (int i = 0; i < 3; i++) {
                    const int nb = t.n[i];
                    if (nb < 0 || fStamp[nb] != st) continue;
                    const uint64_t key = Key(t.v[(i + 1) % 3], t.v[(i + 2) % 3]);
                    if (key != cross && fSeg.count(key)) continue;
                    fStamp[nb] = fCurStamp;
                    keep.push_back(nb);
                }
            }
            fStamp[bad] = -fCurStamp; // never re-enter
            cav = keep;
        }
        // retriangulate: one triangle per boundary edge
        struct BE {
            int a, b, outer, old;
        };
        std::vector<BE> be;
        for (int c : cav) {
            const Tri &t = fT[c];
            for (int i = 0; i < 3; i++) {
                const int nb = t.n[i];
                if (nb >= 0 && fStamp[nb] == fCurStamp) continue;
                be.push_back({t.v[(i + 1) % 3], t.v[(i + 2) % 3], nb, c});
            }
        }
        for (int c : cav) fT[c].alive = false;
        std::unordered_map<int, int> startAt, endAt;
        std::vector<int> created;
        for (const BE &e : be) {
            const int id = NewTri(e.a, e.b, iv);
            Tri &t = fT[id];
            t.n[2] = e.outer;
            if (e.outer >= 0) {
                Tri &o = fT[e.outer];
                for (int j = 0; j < 3; j++)
                    if (o.n[j] == e.old) o.n[j] = id;
            }
            startAt[e.a] = id;
            endAt[e.b] = id;
            created.push_back(id);
        }
        for (int id : created) {
            Tri &t = fT[id];
            t.n[0] = startAt.at(t.v[1]); // across (b, p)
            t.n[1] = endAt.at(t.v[0]);   // across (p, a)
        }
        fHint = created.back();
        for (int id : created) fQueue.push_back(id);
    }

    int AddPoint(const Pt &p) {
        if (int64_t(fP.size()) >= fMaxNodes) throw std::runtime_error("DelaunayMesher: too many nodes");
        fP.push_back(p);
        fVTri.push_back(-1);
        return int(fP.size()) - 1;
    }

    /// the two triangles sharing the edge (a, b) ({-1, -1} when it is not an edge)
    std::array<int, 2> EdgeTris(int a, int b) {
        std::array<int, 2> r = {-1, -1};
        int t = fVTri[a];
        auto scan = [&](int tri) {
            const Tri &tr = fT[tri];
            for (int i = 0; i < 3; i++) {
                const int x = tr.v[(i + 1) % 3], y = tr.v[(i + 2) % 3];
                if ((x == a && y == b) || (x == b && y == a)) return i;
            }
            return -1;
        };
        if (t >= 0 && fT[t].alive) { // rotate around a
            for (int k = 0; k < 64 && t >= 0; k++) {
                const int i = scan(t);
                if (i >= 0) { r[0] = t; r[1] = fT[t].n[i]; return r; }
                const Tri &tr = fT[t];
                int ka = 0;
                while (tr.v[ka] != a) ka++;
                t = tr.n[(ka + 2) % 3]; // across (a, v[ka + 1])
            }
        }
        for (int i = 0; i < int(fT.size()); i++) { // fallback
            if (!fT[i].alive) continue;
            const int j = scan(i);
            if (j >= 0) { r[0] = i; r[1] = fT[i].n[j]; return r; }
        }
        return r;
    }

    /// splits the subsegment key at its midpoint
    void SplitSegment(uint64_t key) {
        const Seg s = fSeg.at(key);
        const Pt m = {0.5 * (fP[s.a].x + fP[s.b].x), 0.5 * (fP[s.a].y + fP[s.b].y)};
        const int iv = AddPoint(m);
        const std::array<int, 2> et = EdgeTris(s.a, s.b);
        std::vector<int> seeds;
        for (int t : et)
            if (t >= 0) seeds.push_back(t);
        if (seeds.empty()) seeds.push_back(Locate(m)); // recovery phase: not an edge yet
        fSeg.erase(key);
        Insert(iv, seeds, key);
        fSeg[Key(s.a, iv)] = {s.a, iv, s.marker};
        fSeg[Key(iv, s.b)] = {iv, s.b, s.marker};
    }

    bool Encroached(const Seg &s, const Pt &q) const {
        const Pt &a = fP[s.a], &b = fP[s.b];
        const double l2 = (b.x - a.x) * (b.x - a.x) + (b.y - a.y) * (b.y - a.y);
        return (a.x - q.x) * (b.x - q.x) + (a.y - q.y) * (b.y - q.y) < -1.e-12 * l2;
    }

    /// is the subsegment encroached by the apex of one of its triangles?
    bool SegmentEncroached(const Seg &s) {
        const std::array<int, 2> et = EdgeTris(s.a, s.b);
        for (int t : et) {
            if (t < 0) continue;
            for (int i = 0; i < 3; i++) {
                const int v = fT[t].v[i];
                if (v >= 4 && v != s.a && v != s.b && Encroached(s, fP[v])) return true;
            }
        }
        return false;
    }

    bool Bad(const Tri &t) const {
        const Pt &a = fP[t.v[0]], &b = fP[t.v[1]], &c = fP[t.v[2]];
        const double l[3] = {std::hypot(b.x - c.x, b.y - c.y), std::hypot(c.x - a.x, c.y - a.y),
                             std::hypot(a.x - b.x, a.y - b.y)};
        const double area = 0.5 * Orient(a, b, c);
        const double R = l[0] * l[1] * l[2] / (4. * area);
        const double lmin = std::min({l[0], l[1], l[2]});
        if (R > fRatio * lmin * (1. + 1.e-9)) return true;
        const Pt g = {(a.x + b.x + c.x) / 3., (a.y + b.y + c.y) / 3.};
        return R * std::sqrt(3.) > fSize(g.x, g.y); // R sqrt(3) = edge of the equilateral triangle
    }

    void Run() {
        // bounding box and two enclosing triangles (vertices 0..3)
        double xmin = 1.e300, xmax = -1.e300, ymin = 1.e300, ymax = -1.e300;
        for (const Pt &p : fPoly) {
            xmin = std::min(xmin, p.x), xmax = std::max(xmax, p.x);
            ymin = std::min(ymin, p.y), ymax = std::max(ymax, p.y);
        }
        const double d = std::max(xmax - xmin, ymax - ymin);
        AddPoint({xmin - 3. * d, ymin - 3. * d});
        AddPoint({xmax + 3. * d, ymin - 3. * d});
        AddPoint({xmax + 3. * d, ymax + 3. * d});
        AddPoint({xmin - 3. * d, ymax + 3. * d});
        const int t0 = NewTri(0, 1, 2), t1 = NewTri(0, 2, 3);
        fT[t0].n[1] = t1; // across (2, 0)
        fT[t1].n[2] = t0; // across (0, 2)
        // boundary points: the polygon edges divided according to the size function (equidistributed 1/h)
        std::vector<int> pid(fPoly.size());
        for (size_t i = 0; i < fPoly.size(); i++) pid[i] = AddPoint(fPoly[i]);
        std::vector<std::array<int, 3>> segs;
        for (size_t i = 0; i < fPoly.size(); i++) {
            const Pt a = fPoly[i], b = fPoly[(i + 1) % fPoly.size()];
            const double L = std::hypot(b.x - a.x, b.y - a.y);
            const int ns = 2000;
            std::vector<double> cum(ns + 1, 0.);
            for (int k = 0; k < ns; k++) {
                const double s = (k + 0.5) / ns;
                cum[k + 1] = cum[k] + L / ns / fSize(a.x + s * (b.x - a.x), a.y + s * (b.y - a.y));
            }
            const int n = std::max(1, int(std::ceil(cum[ns] - 1.e-9)));
            int prev = pid[i];
            for (int k = 1; k < n; k++) {
                const double target = cum[ns] * k / n;
                const int j = int(std::upper_bound(cum.begin(), cum.end(), target) - cum.begin()) - 1;
                const double s = (j + (target - cum[j]) / (cum[j + 1] - cum[j])) / ns;
                const int q = AddPoint({a.x + s * (b.x - a.x), a.y + s * (b.y - a.y)});
                segs.push_back({prev, q, fMarker[i]});
                prev = q;
            }
            segs.push_back({prev, pid[(i + 1) % fPoly.size()], fMarker[i]});
        }
        for (int iv = 4; iv < int(fP.size()); iv++) Insert(iv, {Locate(fP[iv])}, 0);
        fQueue.clear();
        // subsegments: recover the missing ones and split the encroached ones
        for (auto &s : segs) fSeg[Key(s[0], s[1])] = {s[0], s[1], s[2]};
        for (bool changed = true; changed;) {
            changed = false;
            std::vector<uint64_t> keys;
            for (auto &kv : fSeg) keys.push_back(kv.first);
            std::sort(keys.begin(), keys.end());
            for (uint64_t k : keys) {
                if (!fSeg.count(k)) continue;
                const Seg s = fSeg.at(k);
                const std::array<int, 2> et = EdgeTris(s.a, s.b);
                if (et[0] < 0 || SegmentEncroached(s)) { SplitSegment(k); changed = true; }
            }
        }
        // Delaunay refinement of the triangles inside the polygon
        fQueue.clear();
        for (int i = 0; i < int(fT.size()); i++)
            if (fT[i].alive && fT[i].inside) fQueue.push_back(i);
        while (!fQueue.empty()) {
            const int t = fQueue.front();
            fQueue.pop_front();
            if (!fT[t].alive || !fT[t].inside || !Bad(fT[t])) continue;
            const Tri tr = fT[t];
            const Pt c = Circumcenter(fP[tr.v[0]], fP[tr.v[1]], fP[tr.v[2]]);
            std::vector<uint64_t> enc;
            for (auto &kv : fSeg)
                if (Encroached(kv.second, c)) enc.push_back(kv.first);
            if (!enc.empty()) {
                for (uint64_t k : enc)
                    if (fSeg.count(k)) SplitSegment(k);
                // the halves may be encroached by existing vertices
                for (bool changed = true; changed;) {
                    changed = false;
                    std::vector<uint64_t> keys;
                    for (auto &kv : fSeg) keys.push_back(kv.first);
                    for (uint64_t k : keys)
                        if (fSeg.count(k) && SegmentEncroached(fSeg.at(k))) { SplitSegment(k); changed = true; }
                }
                if (fT[t].alive) fQueue.push_back(t);
                continue;
            }
            const int loc = Locate(c);
            if (loc < 0 || !fT[loc].inside) continue; // rounding: circumcenter outside without encroachment
            const int iv = AddPoint(c);
            Insert(iv, {loc}, 0);
        }
    }

    Result Extract() const {
        Result r;
        std::vector<int> map(fP.size(), -1);
        for (const Tri &t : fT) {
            if (!t.alive || !t.inside) continue;
            std::array<int, 3> tri;
            for (int i = 0; i < 3; i++) {
                if (map[t.v[i]] < 0) {
                    map[t.v[i]] = int(r.nodes.size());
                    r.nodes.push_back(fP[t.v[i]]);
                }
                tri[i] = map[t.v[i]];
            }
            r.triangles.push_back(tri);
        }
        for (auto &kv : fSeg) {
            const Seg &s = kv.second;
            if (map[s.a] < 0 || map[s.b] < 0) throw std::runtime_error("DelaunayMesher: lost boundary edge");
            r.edges.push_back({map[s.a], map[s.b], s.marker});
        }
        std::sort(r.edges.begin(), r.edges.end());
        return r;
    }
};

} // namespace slope

#endif
