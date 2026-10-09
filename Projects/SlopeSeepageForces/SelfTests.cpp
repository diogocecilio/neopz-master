// Self-tests of the components (see main.cpp, commands check and verify): manufactured solutions of the anisotropic
// seepage solver and point location; Eq. 21, far-side data, boundary ids and thread safety of the solved field; the
// load vector of the stability problem; the limit analysis (closed forms, quadratures, port of limit_analysis.py,
// reference values); the analytical field (Python reference, Fig. 5, Eq. 40, admissibility). 'check' also runs the
// self-tests of the drivers, CheckFigures (LimitAnalysisCommands.cpp) and CheckFEMBatch (FEMCommands.cpp).
#include "Commands.h"
#include "FEMStability.h"
#include "FigureCases.h"

#include <array>
#include <random>
#include <thread>

namespace slope {

namespace {

// ------------------------------------------------------------------------------------------------------------------
// seepage solver and field
// ------------------------------------------------------------------------------------------------------------------

/// Manufactured solutions: u_ex linear and quadratic with div(K grad u_ex) = 0 (kv x^2 - kh y^2 and xy are
/// K-harmonic only for K = diag(kh, kv) in (x, y): a swapped orientation of K fails), Dirichlet on the ground surface,
/// exact flux (K grad u_ex) . n on the base and the sides. err[case][0..2]: relative errors of u, grad u and J at the
/// P2 nodes and centroids; with locate, also the point-location test (1e6 random points). Max error returned.
REAL VerifyManufactured(const Problem &p, int href, bool locate, REAL err[2][3]) {
    const REAL kh = p.hy.kh, kv = p.hy.kv, H = p.hyd.H;
    REAL worst = 0.;
    struct Exact {
        const char *name;
        REAL c[6]; ///< u = c0 + c1 x + c2 y + c3 (kv x^2 - kh y^2) + c4 x y (c5 unused)
    };
    const Exact cases[2] = {{"linear", {3., 2., -1.5, 0., 0., 0.}},
                            {"quadratic K-harmonic", {3., 2., -1.5, 0.4 / H, 0.3 / H, 0.}}};
    PrintGeometry("hydraulic domain", p.hyd);
    std::cout << "K = diag(" << kh << ", " << kv << "), order " << p.hy.order << " (exact flux on the far sides)\n";
    for (int ic = 0; ic < 2; ic++) {
        const Exact &e = cases[ic];
        err[ic][0] = err[ic][1] = err[ic][2] = 0.;
        if (p.hy.order == 1 && e.c[3] != 0.) continue;
        auto uex = [e, kh, kv](REAL x, REAL y) {
            return e.c[0] + e.c[1] * x + e.c[2] * y + e.c[3] * (kv * x * x - kh * y * y) + e.c[4] * x * y;
        };
        auto gex = [e, kh, kv](REAL x, REAL y, REAL g[2]) {
            g[0] = e.c[1] + 2. * e.c[3] * kv * x + e.c[4] * y;
            g[1] = e.c[2] - 2. * e.c[3] * kh * y + e.c[4] * x;
        };
        auto flux = [gex, kh, kv](REAL nx, REAL ny) {
            return ScalarBC([gex, kh, kv, nx, ny](const TPZVec<REAL> &x) {
                REAL g[2];
                gex(x[0], x[1], g);
                return kh * g[0] * nx + kv * g[1] * ny;
            });
        };
        const ScalarBC ud = [uex](const TPZVec<REAL> &x) { return uex(x[0], x[1]); };
        TPZGeoMesh *gmesh = HydraulicMesh(p, href);
        const auto t0 = std::chrono::steady_clock::now();
        TPZCompMesh *cmesh = CreateSeepageCMesh(gmesh, p.hy, {{EToe, ud}, {ECrest, ud}, {EFace, ud}},
                                                {{EBase, flux(0., -1.)}, {ERight, flux(1., 0.)}, {ELeft, flux(-1., 0.)}});
        const int64_t neq = SolveSeepage(cmesh);
        PoreField pf(cmesh, ESoil, TPZAnisotropicDarcy::EExcessPorePressure);
        const double dt = SecondsSince(t0);
        // errors at the P2 nodes and at the centroids; exact J by the same (exact) edge-midpoint rule
        static const REAL pts[7][2] = {{0., 0.}, {1., 0.}, {0., 1.}, {0.5, 0.}, {0.5, 0.5}, {0., 0.5}, {1. / 3., 1. / 3.}};
        REAL eu = 0., eg = 0., umax = 0., gmax = 0., Jex = 0.;
        for (const PoreField::Tri &t : pf.Triangles()) {
            const REAL det = t.inv[0] * t.inv[3] - t.inv[1] * t.inv[2];
            for (int q = 0; q < 7; q++) {
                const REAL x = t.x0[0] + (t.inv[3] * pts[q][0] - t.inv[1] * pts[q][1]) / det;
                const REAL y = t.x0[1] + (-t.inv[2] * pts[q][0] + t.inv[0] * pts[q][1]) / det;
                REAL u, g[2], ge[2];
                PoreField::EvaluateTri(t, pts[q][0], pts[q][1], u, g);
                gex(x, y, ge);
                eu = std::max(eu, std::fabs(u - uex(x, y)));
                eg = std::max(eg, std::hypot(g[0] - ge[0], g[1] - ge[1]));
                umax = std::max(umax, std::fabs(uex(x, y)));
                gmax = std::max(gmax, std::hypot(ge[0], ge[1]));
                if (q >= 3 && q < 6) Jex += 0.5 * (kh * ge[0] * ge[0] + kv * ge[1] * ge[1]) * t.area / 3.;
            }
        }
        const REAL J = HydraulicFunctional(pf, p.hy);
        err[ic][0] = eu / umax, err[ic][1] = eg / gmax, err[ic][2] = std::fabs(J - Jex) / Jex;
        worst = std::max({worst, err[ic][0], err[ic][1], err[ic][2]});
        std::cout << std::scientific << std::setprecision(3) << e.name << ": " << neq << " equations, " << pf.Triangles().size()
                  << " triangles, max|u - u_ex| / max|u_ex| = " << eu / umax << ", max|grad e| / max|grad u_ex| = "
                  << eg / gmax << ", |J - J_ex| / J_ex = " << std::fabs(J - Jex) / Jex << std::defaultfloat
                  << std::setprecision(6) << "  (" << dt << " s)\n";
        if (!locate) {
            delete cmesh;
            delete gmesh;
            continue;
        }
        // point location (the evaluator of the seepage forces): random points of the bounding box of the soil and
        // of a band of 0.5 H above the ground surface; inside the soil Evaluate must succeed and reproduce u_ex
        const SlopeGeometry &g = p.hyd;
        std::mt19937 gen(11);
        std::uniform_real_distribution<REAL> X(-g.left, g.XT() + g.right), Y(-g.H - g.depth, 0.5 * g.H);
        const int npts = 1000000;
        int64_t nin = 0, nout = 0, miss = 0, spurious = 0;
        REAL elu = 0., elg = 0.;
        const auto t1 = std::chrono::steady_clock::now();
        for (int k = 0; k < npts; k++) {
            TPZManVector<REAL, 3> x = {X(gen), Y(gen), 0.};
            const REAL ys = x[0] <= 0. ? 0. : (x[0] >= g.XT() ? -g.H : -x[0] * g.H / g.XT()); // ground surface
            REAL u, gr[2], ge[2];
            const bool found = pf.Evaluate(x, u, gr);
            const REAL dist = std::min(std::fabs(x[1] - ys), std::fabs(x[0]) + std::fabs(x[0] - g.XT()));
            if (dist < 1.e-9 * g.H) continue; // on the ground surface within rounding
            if (x[1] < ys) {
                nin++;
                if (!found) { miss++; continue; }
                gex(x[0], x[1], ge);
                elu = std::max(elu, std::fabs(u - uex(x[0], x[1])) / umax);
                elg = std::max(elg, std::hypot(gr[0] - ge[0], gr[1] - ge[1]) / gmax);
            } else {
                nout++;
                if (found) spurious++;
            }
        }
        const double tl = SecondsSince(t1);
        std::cout << std::scientific << std::setprecision(3) << "  point location: " << nin << " points in the soil ("
                  << miss << " not found), " << nout << " outside (" << spurious << " found); max error u " << elu
                  << ", grad u " << elg << std::defaultfloat << std::setprecision(4) << "; " << 1.e6 * tl / npts
                  << " us per point\n";
        if (miss > 0 || spurious > 0) worst = std::max<REAL>(worst, 1.);
        worst = std::max({worst, elu, elg});
        delete cmesh;
        delete gmesh;
    }
    return worst;
}

/// manufactured solutions (alpha = 5, P2) exact in the discrete space, and a transposed K detected
void CheckManufactured(Checker &ck, const Problem &p0) {
    std::cout << "manufactured solutions (alpha = 5, P2): exact in the discrete space\n";
    Problem p = p0;
    p.hmesh = p.smesh = "gen";
    p.SetSlope(30., 0.6);
    p.hy.kh = 5. * p.hy.kv;
    REAL err[2][3];
    const REAL worst = VerifyManufactured(p, 0, false, err);
    ck.Expect("    max relative error of u, grad u and J (linear and quadratic K-harmonic)", worst, 1.e-9);
    Problem q = p; // the same with K transposed: the quadratic solution is not harmonic, the error must be large
    q.hy.kh = p.hy.kv, q.hy.kv = p.hy.kh;
    std::cout << "    (K = diag(kv, kh) against the solution K-harmonic for diag(kh, kv): must fail)\n";
    TPZGeoMesh *gm = HydraulicMesh(q, 0);
    const REAL kh = p.hy.kh, kv = p.hy.kv, H = p.hyd.H;
    auto uex = [kh, kv, H](const TPZVec<REAL> &x) { return 0.4 / H * (kv * x[0] * x[0] - kh * x[1] * x[1]); };
    TPZCompMesh *cm = CreateSeepageCMesh(gm, q.hy, {{EToe, uex}, {ECrest, uex}, {EFace, uex}, {ELeft, uex}, {EBase, uex}, {ERight, uex}}, {});
    SolveSeepage(cm);
    PoreField pf(cm, ESoil, TPZAnisotropicDarcy::EExcessPorePressure);
    REAL e = 0., m = 0.;
    for (const PoreField::Tri &t : pf.Triangles()) {
        REAL u, gr[2];
        PoreField::EvaluateTri(t, 1. / 3., 1. / 3., u, gr);
        const REAL det = t.inv[0] * t.inv[3] - t.inv[1] * t.inv[2];
        TPZManVector<REAL, 3> x = {t.x0[0] + (t.inv[3] - t.inv[1]) / (3. * det), t.x0[1] + (-t.inv[2] + t.inv[0]) / (3. * det), 0.};
        e = std::max(e, std::fabs(u - uex(x))), m = std::max(m, std::fabs(uex(x)));
    }
    delete cm;
    delete gm;
    ck.Expect("    transposed K detected: 1e-3 / (max|u - u_ex| / max|u_ex|)", 1.e-3 / (e / m), 1.);
}

/// Dirichlet data of Eq. 21 written independently in the paper coordinates (x right, y_p down) at a point of the
/// ground surface: u = 0 on the crest, u = -gamma_w y_p on the face above the water level (y_p <= h_w), u = -gamma_w h_w
/// on the face below it and on the toe ground
REAL Eq21Paper(const SlopeGeometry &g, REAL gammaw, REAL xp, REAL yp) {
    const REAL tol = 1.e-12 * g.H;
    if (yp <= tol && xp <= tol) return 0.;              // crest
    if (yp >= g.H - tol) return -gammaw * g.hw;         // toe ground
    return yp <= g.hw ? -gammaw * yp : -gammaw * g.hw;  // face
}

/// is the segment a-b on the part of the boundary with the given id?
bool OnBoundaryPart(const SlopeGeometry &g, int id, const REAL a[2], const REAL b[2]) {
    const REAL tol = 1.e-9 * g.H, xT = g.XT();
    for (const REAL *q : {a, b}) {
        const REAL x = q[0], y = q[1];
        bool ok = false;
        switch (id) {
        case EBase: ok = std::fabs(y + g.H + g.depth) < tol; break;
        case ERight: ok = std::fabs(x - xT - g.right) < tol; break;
        case EToe: ok = std::fabs(y + g.H) < tol && x >= xT - tol; break;
        case ECrest: ok = std::fabs(y) < tol && x <= tol; break;
        case ELeft: ok = std::fabs(x + g.left) < tol; break;
        case EFace: { // distance to the line O-T and between O and T
            const REAL L = std::hypot(xT, g.H);
            ok = std::fabs(x * g.H + y * xT) / L < tol && y <= tol && y >= -g.H - tol;
            break;
        }
        default: ok = false;
        }
        if (!ok) return false;
    }
    return true;
}

/// Eq. 21 and the far-side conditions on the solved field, boundary ids, p >= 0, f = -grad u (finite differences),
/// point location on vertices, edges and outside, continuity of the P2 field, thread safety of the evaluator
void CheckSeepageField(Checker &ck, Problem p, REAL beta, REAL hwr, REAL alpha, const std::string &preset) {
    p.hmesh = p.smesh = "gen";
    p.SetSlope(beta, hwr);
    p.hy.kh = alpha * p.hy.kv;
    if (!p.hy.far.SetPreset(preset)) DebugStop();
    const SlopeGeometry &g = p.hyd;
    const REAL gw = p.hy.gammaw, scale = gw * g.H;
    std::cout << "seepage field: beta " << beta << ", h_w / H " << hwr << ", alpha " << alpha << ", far sides " << preset
              << "\n";
    TPZGeoMesh *gmesh = HydraulicMesh(p, 0);
    ck.Expect("    mesh: area, boundary lengths, conformity", CheckSlopeGMesh(gmesh, g), 1.e-10);
    SeepageResult r = DrawdownSeepage(gmesh, g, p.hy);
    const PoreField &pf = *r.field;
    // boundary ids and boundary data
    int64_t misplaced = 0, nmiss = 0, nlines = 0;
    REAL esurf = 0., efar = 0.;
    std::map<int, int64_t> count;
    for (int64_t i = 0; i < gmesh->NElements(); i++) {
        TPZGeoEl *gel = gmesh->Element(i);
        if (!gel || gel->HasSubElement() || gel->Dimension() != 1) continue;
        nlines++;
        const int id = gel->MaterialId();
        count[id]++;
        REAL a[2], b[2];
        for (int d = 0; d < 2; d++) a[d] = gel->NodePtr(0)->Coord(d), b[d] = gel->NodePtr(1)->Coord(d);
        if (!OnBoundaryPart(g, id, a, b)) misplaced++;
        const EFarBC far = id == ELeft ? p.hy.far.left : (id == EBase ? p.hy.far.bottom : p.hy.far.right);
        for (REAL t : {0., 0.25, 0.5, 0.75, 1.}) {
            TPZManVector<REAL, 3> x = {a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1]), 0.};
            REAL u, gr[2];
            if (!pf.Evaluate(x, u, gr)) {
                nmiss++;
                continue;
            }
            if (id == EToe || id == ECrest || id == EFace) esurf = std::max(esurf, std::fabs(u - Eq21Paper(g, gw, x[0], -x[1])) / scale);
            else if (far != EFarBC::ENoFlow) efar = std::max(efar, std::fabs(u - (far == EFarBC::EZero ? 0. : -gw * g.hw)) / scale);
        }
    }
    ck.Expect("    boundary lines off their part of the boundary (ids -1 .. -6)", REAL(misplaced), 0.);
    ck.Expect("    boundary ids missing (each of -1 .. -6 present)", REAL(6 - count.size()), 0.);
    ck.Expect("    boundary points not located", REAL(nmiss), 0.);
    ck.Expect("    Eq. 21 on crest, face and toe ground: max|u - u_D| / (gamma_w H)", esurf, 1.e-8);
    ck.Expect("    far-side Dirichlet data: max|u - u_D| / (gamma_w H)", efar, 1.e-8);
    ck.Expect("    total pore pressure p = u - gamma_w y >= 0: -min p / (gamma_w H)", -r.pmin / scale, 1.e-9);
    // f = -grad u: central differences of u inside random triangles; the adapter of paper-coordinate fields
    std::mt19937 gen(7);
    std::uniform_real_distribution<REAL> U(0., 1.);
    const auto &tris = pf.Triangles();
    REAL efd = 0., gmax = 0.;
    for (int k = 0; k < 2000; k++) {
        const PoreField::Tri &t = tris[size_t(U(gen) * tris.size()) % tris.size()];
        REAL a = U(gen), b = U(gen);
        if (a + b > 1.) a = 1. - a, b = 1. - b;
        const REAL l1 = 0.1 + 0.7 * a, l2 = 0.1 + 0.7 * b; // all barycentric coordinates >= 0.1
        const REAL det = t.inv[0] * t.inv[3] - t.inv[1] * t.inv[2];
        const REAL x = t.x0[0] + (t.inv[3] * l1 - t.inv[1] * l2) / det, y = t.x0[1] + (-t.inv[2] * l1 + t.inv[0] * l2) / det;
        const REAL h = 1.e-4 * std::sqrt(t.area);
        REAL f[2], up = 0., um = 0., gr[2], d[2];
        pf.Force(TPZManVector<REAL, 3>{x, y, 0.}, f);
        for (int dir = 0; dir < 2; dir++) {
            TPZManVector<REAL, 3> xp = {x + (dir == 0 ? h : 0.), y + (dir == 1 ? h : 0.), 0.};
            TPZManVector<REAL, 3> xm = {x - (dir == 0 ? h : 0.), y - (dir == 1 ? h : 0.), 0.};
            const bool found = pf.Evaluate(xp, up, gr) && pf.Evaluate(xm, um, gr);
            d[dir] = found ? (up - um) / (2. * h) : 1.e300; // not found: the check fails
        }
        efd = std::max(efd, std::hypot(f[0] + d[0], f[1] + d[1]));
        gmax = std::max(gmax, std::hypot(d[0], d[1]));
    }
    ck.Expect("    f = -grad u (central differences of u): max|f + grad_h u| / max|grad u|", efd / gmax, 1.e-5);
    {
        // u_paper = 2 + 3 x - 5 y_p: f_paper = (-3, 5); in NeoPZ (y = -y_p) f = (-3, -5)
        const ForceField fp = FromPaperCoordinates([](REAL, REAL, REAL fpp[2]) { fpp[0] = -3., fpp[1] = 5.; });
        REAL f[2];
        fp(TPZManVector<REAL, 3>{1., -2., 0.}, f);
        ck.Expect("    FromPaperCoordinates: |f - (-3, -5)|", std::hypot(f[0] + 3., f[1] + 5.), 0.);
    }
    // location of the vertices, edge midpoints and centroid of every triangle: found, and the same P2 value from
    // whichever triangle the search returns (continuity across the edges)
    {
        static const REAL qp[7][2] = {{0., 0.}, {1., 0.}, {0., 1.}, {0.5, 0.}, {0.5, 0.5}, {0., 0.5}, {1. / 3., 1. / 3.}};
        int64_t miss = 0;
        REAL ejump = 0.;
        for (const PoreField::Tri &t : tris) {
            const REAL det = t.inv[0] * t.inv[3] - t.inv[1] * t.inv[2];
            for (auto &q : qp) {
                const REAL x = t.x0[0] + (t.inv[3] * q[0] - t.inv[1] * q[1]) / det;
                const REAL y = t.x0[1] + (-t.inv[2] * q[0] + t.inv[0] * q[1]) / det;
                REAL uo, go[2], u, gr[2];
                PoreField::EvaluateTri(t, q[0], q[1], uo, go);
                if (!pf.Evaluate(TPZManVector<REAL, 3>{x, y, 0.}, u, gr)) {
                    miss++;
                    continue;
                }
                ejump = std::max(ejump, std::fabs(u - uo) / scale);
            }
        }
        ck.Expect("    vertices / edge midpoints / centroids not located", REAL(miss), 0.);
        ck.Expect("    continuity of u at shared vertices and edge midpoints / (gamma_w H)", ejump, 1.e-11);
    }
    // points just outside: above the ground surface and beyond the box, by 1e-6 H
    {
        std::uniform_real_distribution<REAL> X(-g.left, g.XT() + g.right), Y(-g.H - g.depth, 0.);
        const REAL d = 1.e-6 * g.H;
        int64_t spurious = 0, n = 0;
        for (int k = 0; k < 20000; k++, n += 5) {
            const REAL x = X(gen), y = Y(gen);
            REAL ys = x <= 0. ? 0. : (x >= g.XT() ? -g.H : -x * g.H / g.XT()), u, gr[2];
            if (g.XT() > 0. && x > 0. && x < g.XT()) ys += d * std::hypot(g.XT(), g.H) / g.XT(); // normal offset d
            else ys += d;
            const REAL out[5][2] = {{x, ys}, {-g.left - d, y}, {g.XT() + g.right + d, y}, {x, -g.H - g.depth - d}, {x, 1.e3 * g.H}};
            for (auto &q : out) spurious += pf.Evaluate(TPZManVector<REAL, 3>{q[0], q[1], 0.}, u, gr);
        }
        if (g.XT() == 0.) { // vertical face: points just to the right of it, between the toe and the crest
            for (int k = 0; k < 20000; k++, n++) {
                REAL u, gr[2];
                spurious += pf.Evaluate(TPZManVector<REAL, 3>{d, -g.H * (1.e-3 + (1. - 2.e-3) * U(gen)), 0.}, u, gr);
            }
        }
        ck.Expect("    points 1e-6 H outside the soil found (of " + std::to_string(n) + ")", REAL(spurious), 0.);
    }
    // thread safety: the evaluator called from 4 threads gives bitwise the serial values
    {
        std::uniform_real_distribution<REAL> X(-g.left, g.XT() + g.right), Y(-g.H - g.depth, 0.);
        const int npts = 400000, nth = 4;
        std::vector<REAL> xs(2 * npts), f1(2 * npts), f4(2 * npts);
        for (int k = 0; k < npts; k++) xs[2 * k] = X(gen), xs[2 * k + 1] = Y(gen);
        const ForceField ff = r.field->AsForceField(r.field);
        for (int k = 0; k < npts; k++) ff(TPZManVector<REAL, 3>{xs[2 * k], xs[2 * k + 1], 0.}, &f1[2 * k]);
        std::vector<std::thread> th;
        for (int t = 0; t < nth; t++)
            th.emplace_back([&, t]() {
                for (int k = t; k < npts; k += nth) ff(TPZManVector<REAL, 3>{xs[2 * k], xs[2 * k + 1], 0.}, &f4[2 * k]);
            });
        for (auto &t : th) t.join();
        int64_t diff = 0;
        for (int k = 0; k < 2 * npts; k++) diff += f1[k] != f4[k];
        ck.Expect("    force field from 4 threads: values differing from the serial ones", REAL(diff), 0.);
    }
    delete gmesh;
}

// ------------------------------------------------------------------------------------------------------------------
// load vector of the stability problem (FEMStability.h)
// ------------------------------------------------------------------------------------------------------------------

/// int_boundary u n ds over the counterclockwise polygon of g (outward normal), m subintervals per edge, 3 Gauss
/// points each (exact for a P2 field when the subintervals contain no element vertex); nmiss: points not located
void BoundaryIntegralUN(const PoreField &pf, const SlopeGeometry &g, int m, REAL I[2], int64_t &nmiss) {
    static const REAL gp[3] = {-std::sqrt(0.6), 0., std::sqrt(0.6)}, gw[3] = {5. / 9., 8. / 9., 5. / 9.};
    std::vector<int> marker;
    const std::vector<Pt> poly = g.Polygon(marker);
    I[0] = I[1] = 0.;
    nmiss = 0;
    for (size_t i = 0; i < poly.size(); i++) {
        const Pt &a = poly[i], &b = poly[(i + 1) % poly.size()];
        const REAL dx = b.x - a.x, dy = b.y - a.y; // n ds = (dy, -dx) dt
        for (int k = 0; k < m; k++)
            for (int q = 0; q < 3; q++) {
                const REAL t = (k + 0.5 * (1. + gp[q])) / m, w = 0.5 * gw[q] / m;
                TPZManVector<REAL, 3> x = {a.x + t * dx, a.y + t * dy, 0.};
                REAL u, gr[2];
                if (!pf.Evaluate(x, u, gr)) {
                    nmiss++;
                    continue;
                }
                I[0] += w * u * dy, I[1] -= w * u * dx;
            }
    }
}

REAL PolygonArea(const SlopeGeometry &g) {
    std::vector<int> marker;
    const std::vector<Pt> poly = g.Polygon(marker);
    REAL A = 0.;
    for (size_t i = 0; i < poly.size(); i++) {
        const Pt &a = poly[i], &b = poly[(i + 1) % poly.size()];
        A += 0.5 * (a.x * b.y - b.x * a.y);
    }
    return A;
}

REAL RelDiff(const TPZFMatrix<STATE> &a, const TPZFMatrix<STATE> &b) {
    TPZFMatrix<STATE> d(a);
    d -= b;
    return Norm(d) / Norm(a);
}

/// Load vector of the stability problem: resultant against the divergence theorem int_Omega -grad u = -int_dOmega u n,
/// lambda scales gamma' AND f, forms u and p identical, threaded assembly equal to the serial one
void CheckLoadVector(Checker &ck, Problem p, REAL beta, REAL hwr, REAL alpha, bool trig, REAL tolDiv) {
    p.hmesh = p.smesh = trig ? "trig" : "gen";
    if (trig) p.hyd.H = 10.;
    p.SetSlope(trig ? 45. : beta, hwr);
    p.hy.kh = alpha * p.hy.kv;
    const SlopeGeometry &g = p.stab;
    std::cout << "load vector: meshes " << p.hmesh << ", H " << g.H << ", beta " << g.beta << ", h_w / H " << hwr
              << ", alpha " << alpha << ", hydraulic order " << p.hy.order << "\n";
    const SeepageResult r = SolveDrawdown(p, 0);
    TPZGeoMesh *smesh = StabilityMesh(p);
    ck.Expect("    stability mesh: area, boundary lengths, conformity", CheckSlopeGMesh(smesh, g), 1.e-10);
    const TMCVoigt model = ModelVoigt(p.soil);
    const REAL gs = p.soil.gamma, gw = p.hy.gammaw, gp = gs - gw, A = PolygonArea(g);
    {
        TPZCompMesh *cmesh = CreateCMesh(smesh, 2, model, p.soil);
        int64_t nneg, ntot, nout;
        REAL pmin;
        NegativePressurePoints(cmesh, *r.field, gw, nneg, ntot, pmin, nout);
        delete cmesh;
        ck.Expect("    integration points of the stability mesh outside the hydraulic mesh", REAL(nout), 0.);
    }
    const ForceField fu = r.field->AsForceField(r.field), fp = TotalPressureForce(r.field, gw, false);
    REAL R1[2], Rd[2], I[2];
    const TPZFMatrix<STATE> F1 = LoadVector(smesh, model, p.soil, fu, gp, 1., 0, R1);
    const TPZFMatrix<STATE> F27 = LoadVector(smesh, model, p.soil, fu, gp, 2.7);
    const TPZFMatrix<STATE> Fp = LoadVector(smesh, model, p.soil, fp, gs);
    const TPZFMatrix<STATE> F4 = LoadVector(smesh, model, p.soil, fu, gp, 1., 4);
    const TPZFMatrix<STATE> Fd = LoadVector(smesh, model, p.soil, NoSeepage(), gp, 1., 0, Rd);
    int64_t nmiss;
    BoundaryIntegralUN(*r.field, g, 672, I, nmiss);
    const REAL Rex[2] = {-I[0], -I[1] - gp * A}; // int (f + gamma' g) with int f = -int u n ds
    const REAL Rn = std::hypot(Rex[0], Rex[1]);
    std::cout << std::setprecision(10) << "    resultant at lambda = 1: load vector (" << R1[0] << ", " << R1[1]
              << "), divergence theorem (" << Rex[0] << ", " << Rex[1] << "); seepage part (" << R1[0] - Rd[0] << ", "
              << R1[1] - Rd[1] << ") kN/m" << std::setprecision(6) << "\n";
    ck.Expect("    boundary points of the stability domain not located in the hydraulic mesh", REAL(nmiss), 0.);
    ck.Expect("    gravity only: |R + gamma' A e_y| / (gamma' A)", std::hypot(Rd[0], Rd[1] + gp * A) / (gp * A), 1.e-12);
    ck.Expect("    resultant vs -int u n ds - gamma' A e_y: relative difference",
              std::hypot(R1[0] - Rex[0], R1[1] - Rex[1]) / Rn, tolDiv);
    TPZFMatrix<STATE> F1s(F1);
    F1s *= 2.7;
    ck.Expect("    lambda = 2.7 scales the whole load (gamma' and f): |F(2.7) - 2.7 F(1)| / |F(2.7)|", RelDiff(F27, F1s), 1.e-14);
    ck.Expect("    forms u and p: |F_u - F_p| / |F_u|", RelDiff(F1, Fp), 1.e-13);
    ck.Expect("    assembly on 4 threads vs serial: |F_4 - F_1| / |F_1|", RelDiff(F1, F4), 1.e-13);
    delete smesh;
}

// ------------------------------------------------------------------------------------------------------------------
// kinematic limit analysis (LimitAnalysis.h, port of scripts/limit_analysis.py)
// ------------------------------------------------------------------------------------------------------------------

/// random admissible mechanisms of a class (uniform in the search box)
std::vector<la::Mechanism> RandomMechanisms(const la::Problem &prob, int kind, int n, std::mt19937_64 &gen) {
    std::uniform_real_distribution<REAL> U(0., 1.);
    REAL lb[3], ub[3];
    prob.Bounds(kind, lb, ub);
    std::vector<la::Mechanism> out;
    for (int it = 0; it < 4000000 && int(out.size()) < n; it++) {
        REAL x[3];
        for (int d = 0; d < 3; d++) x[d] = lb[d] + U(gen) * (ub[d] - lb[d]);
        if (prob.Feasible(x, kind)) out.push_back(prob.Mech(x));
    }
    return out;
}

/// int_Omega (Cx - x) dA (= P_gamma / gamma') by the shoelace formulas of the polygon A-O-(T)-B + n-chord spiral,
/// Richardson-extrapolated in the number of chords (independent of the closed forms and of the polar quadrature)
REAL PolygonMoment(const la::Mechanism &m, const la::Geometry &g, int n) {
    auto moment = [&](int ns) {
        std::vector<std::array<REAL, 2>> P = {{m.Ax, 0.}, {0., 0.}};
        if (m.II) P.push_back({g.xtoe, g.H});
        P.push_back({m.Bx, m.By});
        for (int i = 1; i < ns - 1; i++) {
            const REAL th = m.t2 + (m.t1 - m.t2) * i / (ns - 1.), r = m.R(th, g.k);
            P.push_back({m.Cx - r * std::cos(th), -m.Cy + r * std::sin(th)});
        }
        REAL A = 0., Mx = 0.;
        for (size_t i = 0; i < P.size(); i++) {
            const auto &a = P[i], &b = P[(i + 1) % P.size()];
            const REAL cr = a[0] * b[1] - b[0] * a[1];
            A += 0.5 * cr, Mx += (a[0] + b[0]) * cr / 6.;
        }
        if (A < 0.) A = -A, Mx = -Mx;
        return m.Cx * A - Mx;
    };
    const REAL a = moment(n), b = moment(2 * n - 1);
    return b + (b - a) / 3.;
}

/// limit analysis: closed forms against the polar quadrature and an independent polygon integration, P_mr, domain
/// quadrature against the boundary formula, values of scripts/limit_analysis.py at fixed mechanisms (port), sampled
/// admissibility, dry stability numbers and Fig. 8 h_w = 0 ends against the Python reference and the paper, scale
/// invariance, uniform field, threads, FE field (boundary vs domain, Gamma vs the Python reference)
void CheckLimitAnalysis(Checker &ck, const Problem &p0) {
    std::cout << "limit analysis (LimitAnalysis.h)\n";
    std::mt19937_64 gen(12345);
    const la::Quadrature xfine = la::Quadrature::Level("xfine");
    { // closed forms vs quadrature vs polygon; P_u domain vs boundary for a smooth potential; sampled admissibility
        REAL eq = 0., epol = 0., e53 = 0., epot = 0., above = 0.;
        int n = 0;
        const REAL a = 0.7, b = 0.3, cc = 1.1, d = -0.4; // u = a x^2 y - b y^3 + c x y + d x (paper coordinates)
        auto up = [&](REAL x, REAL y) { return a * x * x * y - b * y * y * y + cc * x * y + d * x; };
        auto fp = [&](REAL x, REAL y, REAL f[2]) { f[0] = -(2 * a * x * y + cc * y + d), f[1] = -(a * x * x - 3 * b * y * y + cc * x); };
        for (REAL beta : {30., 60., 90.})
            for (REAL phi : {0., 20., 32.})
                for (int kind : {1, 2}) {
                    const la::Problem prob(beta, 1., 1., phi, 1., 0., la::Seepage());
                    for (const la::Mechanism &m : RandomMechanisms(prob, kind, 25, gen)) {
                        if (m.rh > 20.) continue; // huge near-degenerate mechanisms: P_gamma ~ 1e-4 of the f terms
                        n++;
                        const REAL pc = la::PgammaClosed(m, prob.geo, 1.);
                        eq = std::max(eq, std::fabs(prob.PgammaQuadrature(m, xfine) - pc) / std::fabs(pc));
                        epol = std::max(epol, std::fabs(PolygonMoment(m, prob.geo, 2001) - pc) / std::fabs(pc));
                        if (kind == 1 && beta < 90.) {
                            REAL f[4];
                            la::FTerms(m, prob.geo, f);
                            e53 = std::max(e53, std::fabs(la::F3Paper(m, prob.geo) - f[2]) / std::fabs(f[2]));
                        }
                        const REAL pd = la::DomainPower(m, prob.geo, xfine, fp, {}), pb = la::BoundaryPower(m, prob.geo, up, 0.37, 48, 10);
                        epot = std::max(epot, std::fabs(pd - pb) / std::fabs(pd));
                        for (int i = 1; i < 400; i++) { // spiral below the ground surface (paper y >= y_surface)
                            const REAL th = m.t1 + (m.t2 - m.t1) * i / 400., r = m.R(th, prob.geo.k);
                            const REAL x = m.Cx - r * std::cos(th), y = -m.Cy + r * std::sin(th);
                            const REAL ys = x <= 0. ? 0. : (x >= prob.geo.xtoe ? 1. : x * prob.geo.sb / prob.geo.cb);
                            above = std::max(above, ys - y);
                        }
                    }
                }
        std::cout << "    " << n << " random admissible mechanisms (beta 30/60/90, phi 0/20/32, I and II, r_h <= 20 H)\n";
        ck.Expect("    P_gamma: closed form f1 - f2 - f3 (- f4) vs polar quadrature (xfine), max rel", eq, 1.e-9);
        ck.Expect("    P_gamma: closed form vs polygon A-O-(T)-B + spiral chords (shoelace), max rel", epol, 1.e-8);
        ck.Expect("    f3: Eq. 53 as printed vs the sin form (I, beta < 90), max rel", e53, 1.e-9);
        ck.Expect("    P_u of u = 0.7 x^2 y - 0.3 y^3 + 1.1 x y - 0.4 x: domain (xfine) vs boundary formula, max rel", epot, 1.e-9);
        ck.Expect("    sampled spiral points above the ground surface: max height / H", above, 1.e-9);
    }
    { // P_mr: Eq. 48 vs Gauss quadrature of c r^2 dtheta (= c cos(phi) |U| ds), continuity at phi -> 0
        REAL e = 0.;
        std::vector<REAL> x, w;
        la::Gauss01(40, x, w);
        for (REAL phi : {0., 30.}) {
            const la::Geometry g(1., 60., phi);
            const la::Mechanism m = la::MakeMechanism(g, 0.5, 1.6, 1.);
            REAL q = 0.;
            for (size_t i = 0; i < x.size(); i++) {
                const REAL r = m.R(m.t1 + (m.t2 - m.t1) * x[i], g.k);
                q += (m.t2 - m.t1) * w[i] * r * r;
            }
            e = std::max(e, std::fabs(la::Pmr(m, g, 1.) - q) / q);
        }
        const la::Geometry g0(1., 60., 0.), g1(1., 60., 1.e-7);
        const REAL p0v = la::Pmr(la::MakeMechanism(g0, 0.5, 1.6, 1.), g0, 1.), p1v = la::Pmr(la::MakeMechanism(g1, 0.5, 1.6, 1.), g1, 1.);
        ck.Expect("    P_mr: Eq. 48 vs quadrature of c r^2 dtheta along the spiral (phi = 0, 30), max rel", e, 1.e-13);
        ck.Expect("    P_mr continuity: |P_mr(phi = 1e-7 deg) - P_mr(phi = 0)| / P_mr", std::fabs(p1v - p0v) / p0v, 1.e-8);
    }
    { // port: values of scripts/limit_analysis.py at fixed mechanisms (printed with 17 digits)
        struct Ref {
            REAL beta, phi, x[3];
            const char *what;
            REAL value;
        };
        const Ref refs[] = {
            {60., 30., {0.62, 1.75, 1.}, "P_mr", 1.3593651407732807},
            {60., 30., {0.62, 1.75, 1.}, "Pg", 0.07737538847106358},
            {60., 30., {0.62, 1.75, 1.}, "poly coarse", -0.077399121888196115},
            {60., 30., {0.62, 1.75, 1.}, "circ fine", 0.21673871318454482},
            {90., 0., {0.8, 1.3, 0.6}, "circ xfine", 0.64175807184453182},
            {90., 0., {0.8, 1.3, 0.6}, "circ2 coarse", 0.45795732405646439},
            {45., 20., {0.9, 2.2, 1.8}, "Pg", -0.005294751116120211},
            {45., 20., {0.9, 2.2, 1.8}, "circ fine", 0.15376831158339679},
            {45., 20., {0.9, 2.2, 1.8}, "circ2 fine", -0.77577934565328532},
            {45., 20., {0.9, 2.2, 1.8}, "boundary 24", -1.0023548344041124},
            {30., 24.7, {1., 2.1, 0.5}, "circ2 xfine", 0.0050952195358424301},
            {35., 32., {1., 2.05, 0.2}, "Pg coarse", -0.00033676399952368134}};
        auto circ = [](REAL qx, REAL qy, REAL R) { // f = (0.3, 1.7) inside the circle, (-0.4, 0.9) outside
            return [=](REAL x, REAL y, REAL f[2]) {
                const bool in = std::hypot(x - qx, y - qy) < R;
                f[0] = in ? 0.3 : -0.4, f[1] = in ? 1.7 : 0.9;
            };
        };
        const REAL a = 0.7, b = 0.3, cc = 1.1, d = -0.4;
        auto up = [&](REAL x, REAL y) { return a * x * x * y - b * y * y * y + cc * x * y + d * x; };
        auto fp = [&](REAL x, REAL y, REAL f[2]) { f[0] = -(2 * a * x * y + cc * y + d), f[1] = -(a * x * x - 3 * b * y * y + cc * x); };
        REAL e = 0.;
        for (const Ref &r : refs) {
            const la::Geometry g(1., r.beta, r.phi);
            const la::Mechanism m = la::MakeMechanism(g, r.x[0], r.x[1], r.x[2]);
            const std::string w = r.what;
            const std::string lev = w.find("coarse") != std::string::npos ? "coarse" : (w.find("xfine") != std::string::npos ? "xfine" : "fine");
            const la::Quadrature q = la::Quadrature::Level(lev);
            REAL v = 0.;
            if (w == "P_mr") v = la::Pmr(m, g, 1.);
            else if (w == "Pg") v = la::PgammaClosed(m, g, 1.);
            else if (w == "Pg coarse") v = la::DomainPower(m, g, q, [](REAL, REAL, REAL f[2]) { f[0] = 0., f[1] = 1.; }, {});
            else if (w.rfind("poly", 0) == 0) v = la::DomainPower(m, g, q, fp, {});
            else if (w.rfind("circ2", 0) == 0) v = la::DomainPower(m, g, q, circ(0., 0., 0.5), {{0., 0., 0.5, true}, {0., 0., 0.15, false}});
            else if (w.rfind("circ", 0) == 0) v = la::DomainPower(m, g, q, circ(0.1, 0.4, 0.7), {{0.1, 0.4, 0.7, true}});
            else if (w == "boundary 24") v = la::BoundaryPower(m, g, up, 0.37, 24, 8);
            e = std::max(e, std::fabs(v - r.value) / std::fabs(r.value));
        }
        ck.Expect("    port: limit_analysis.py values at fixed mechanisms (P_mr, P_gamma, P_u domain with circle jumps, "
                  "boundary formula), max rel",
                  e, 1.e-12);
    }
    { // dry stability numbers N = gamma H_c / c (gamma = c = H = 1, gamma_w = 0) vs limit_analysis.py
        // (results/limit_analysis/dry_stability_numbers.csv, 6 decimals) and the vertical cut of Chen (1975), 3.83
        const REAL ref[5][3] = {{90., 0., 3.831337}, {60., 30., 16.035178}, {45., 20., 16.160944}, {30., 25., 119.920313}, {45., 0., 5.525981}};
        REAL e = 0., N90 = 0.;
        bool dmax = false;
        for (auto &r : ref) {
            const la::Result res = la::StabilityFactor(la::Problem(r[0], 1., 1., r[1], 1., 0., la::Seepage()));
            e = std::max(e, std::fabs(res.Gamma / r[2] - 1.));
            if (r[0] == 90.) N90 = res.Gamma;
            if (r[0] == 45. && r[1] == 0.) dmax = res.dmaxReached && res.kind == 2;
        }
        ck.Expect("    dry N = gamma H_c / c (beta/phi 90/0, 60/30, 45/20, 30/25, 45/0) vs limit_analysis.py: max rel", e, 2.e-7);
        ck.Expect("    vertical cut, phi = 0: |N - 3.83| / 3.83 (Chen 1975)", std::fabs(N90 - 3.83) / 3.83, 1.e-3);
        ck.Expect("    beta = 45, phi = 0: optimum = mechanism II at d = d_max (flagged)", dmax ? 0. : 1., 0.);
        const la::Result inf = la::StabilityFactor(la::Problem(30., 10., 6., 32., 18., 9.8, la::Seepage()));
        ck.Expect("    beta = 30 < phi = 32 without seepage: no mechanism with P_ext > 0 (Gamma = inf)", inf.found ? 1. : 0., 0.);
    }
    { // Fig. 8 at h_w = 0 (buoyant, f = 0), (c, phi) of Table 1 swapped between the panels, gamma_w = 9.8 (SPEC)
        struct Panel {
            REAL beta, c, phi, python, paper;
        };
        const Panel panels[4] = {{30., 11.7, 24.7, 156.52744, 156.555}, {60., 11.7, 24.7, 17.94946, 17.949},
                                 {35., 6., 32., 228.82857, 229.319}, {60., 6., 32., 12.98603, 12.986}};
        REAL ep = 0., epap = 0.;
        for (const Panel &q : panels) {
            const la::Result r = la::StabilityFactor(la::Problem(q.beta, 10., q.c, q.phi, 18., 9.8, la::Seepage()));
            ep = std::max(ep, std::fabs(r.Hcrit / q.python - 1.)), epap = std::max(epap, std::fabs(r.Hcrit / q.paper - 1.));
            std::cout << std::setprecision(8) << "    Fig. 8, h_w = 0, beta " << q.beta << ", c " << q.c << ", phi " << q.phi
                      << ": H_crit " << r.Hcrit << " m (Python " << q.python << ", paper " << q.paper << ")\n";
        }
        ck.Expect("    Fig. 8 h_w = 0 ends: H_crit vs limit_analysis.py (fig8_hw0_buoyant.csv), max rel", ep, 1.e-6);
        ck.Expect("    Fig. 8 h_w = 0 ends: H_crit vs the paper (data/paper_fig8.csv), max rel", epap, 2.5e-3);
    }
    { // scale invariance and threads
        const la::Result r5 = la::StabilityFactor(la::Problem(60., 5., 6., 32., 18., 9.81, la::Seepage()));
        const la::Result r10 = la::StabilityFactor(la::Problem(60., 10., 6., 32., 18., 9.81, la::Seepage()));
        ck.Expect("    scale invariance (dry): |Gamma(H = 5) / Gamma(H = 10) - 2|", std::fabs(r5.Gamma / r10.Gamma - 2.), 1.e-9);
        ck.Expect("    Gamma(H = 10) vs limit_analysis.py 1.30018873 (brute-force verified)", std::fabs(r10.Gamma / 1.30018873 - 1.), 1.e-8);
    }
    { // uniform field f = (0, g) (paper, downwards): Gamma(gamma', f) = Gamma(gamma' + g, no field), both P_u methods
        la::Seepage s;
        const REAL g = 4.;
        s.force = [g](const TPZVec<REAL> &, REAL f[2]) { f[0] = 0., f[1] = -g; }; // NeoPZ: y up
        s.u = [g](const TPZVec<REAL> &x, REAL &u) {
            u = g * x[1]; // f = -grad u = (0, -g) in NeoPZ coordinates
            return true;
        };
        s.gradient = true;
        la::Settings st;
        st.seeds = {0};
        const la::Result rd = la::StabilityFactor(la::Problem(60., 5., 6., 32., 18., 9.81, s, "domain"), st);
        const la::Result rb = la::StabilityFactor(la::Problem(60., 5., 6., 32., 18., 9.81, s, "boundary"), st);
        const la::Result rg = la::StabilityFactor(la::Problem(60., 5., 6., 32., 18. + g, 9.81, la::Seepage()), st);
        ck.Expect("    uniform field (0, g): |Gamma(domain P_u) / Gamma(gamma' + g) - 1|", std::fabs(rd.Gamma / rg.Gamma - 1.), 1.e-9);
        ck.Expect("    uniform field (0, g): |Gamma(boundary P_u) / Gamma(gamma' + g) - 1|", std::fabs(rb.Gamma / rg.Gamma - 1.), 1.e-9);
        la::Settings s1 = st, s4 = st;
        s1.seeds = s4.seeds = {0, 1, 2};
        s1.nThreads = 1, s4.nThreads = 4;
        const la::Problem pd(60., 5., 6., 32., 18., 9.81, s, "domain");
        const la::Result t1 = la::StabilityFactor(pd, s1), t4 = la::StabilityFactor(pd, s4);
        ck.Expect("    runs on 4 threads vs serial: Gamma and x bitwise equal (differences)",
                  REAL((t1.Gamma != t4.Gamma) + (t1.x[0] != t4.x[0]) + (t1.x[1] != t4.x[1]) + (t1.x[2] != t4.x[2])), 0.);
    }
    { // FE field of Fig. 9 (beta = 60, h_w = H = 5 m, alpha = 1, box 50 / 10 / 30 m, zero_lb): P_u boundary vs domain
        // at the Python optimum, and Gamma vs scripts/reproduce_python.py (results/reproduce_python/cache.jsonl, FE_box_m)
        Problem p = p0;
        p.hmesh = p.smesh = "gen";
        p.hyd.H = 5.;
        p.hext[0] = 10., p.hext[1] = 2., p.hext[2] = 6.;
        p.hy.kh = p.hy.kv, p.hy.gammaw = 9.81;
        if (!p.hy.far.SetPreset("zero_lb")) DebugStop();
        p.SetSlope(60., 1.);
        const LAField lf = MakeLAField(p, "fe", 0);
        const la::Problem prob(60., 5., 10., 30., 20., 9.81, lf.seep);
        const REAL xpy[3] = {0.6191997989032978, 1.749558364700103, 1.0};
        const la::Mechanism m = prob.Mech(xpy);
        const REAL pb = prob.Pu(m, la::Quadrature::Level("fine"), la::PuMethod::EBoundary);
        const REAL pdom = prob.Pu(m, la::Quadrature::Level("ref"), la::PuMethod::EDomain);
        std::cout << std::setprecision(10) << "    FE field (Fig. 9, beta 60, box 50/10/30 m): P_u at the Python optimum: boundary "
                  << pb << ", domain (ref) " << pdom << ", Python 256.1165173\n";
        ck.Expect("    FE field: P_u domain quadrature (ref) vs boundary formula, rel", std::fabs(pdom / pb - 1.), 2.e-5);
        ck.Expect("    FE field: P_u (C++ FE field) vs Python (fe_seepage.py field) at the same mechanism, rel",
                  std::fabs(pb / 256.1165173225788 - 1.), 1.e-3);
        la::Settings st;
        st.seeds = {0};
        const la::Result r = la::StabilityFactor(prob, st);
        std::cout << "    Gamma " << r.Gamma << " (Python 0.95704370), " << r.seconds << " s\n";
        ck.Expect("    FE field: Gamma vs limit_analysis.py with fe_seepage.py (0.95704370), rel", std::fabs(r.Gamma / 0.9570437 - 1.), 2.e-3);
    }
}

// ------------------------------------------------------------------------------------------------------------------
// analytical seepage field K^-1 v'_opt (AnalyticalSeepage.h)
// ------------------------------------------------------------------------------------------------------------------

/// a case of the reference file of scripts/analytical_seepage_reference.py: data, scalars of the optimum and f at points
struct AnalyticalRefCase {
    REAL beta = 0., alpha = 1., H = 1., hw = 0., gammaw = 9.81, lm = 10., m = 0., Cn = 0., Dn = 0., F = 0., Jn = 0.;
    bool degenerate = false;
    std::vector<std::array<REAL, 4>> pts; ///< x, y, fx, fy in the paper coordinates (y down)
};

std::vector<AnalyticalRefCase> ReadAnalyticalReference(const std::string &file) {
    std::vector<AnalyticalRefCase> cases;
    std::ifstream in(file);
    std::string line;
    std::vector<REAL> v;
    while (std::getline(in, line)) {
        if (line.size() < 2 || line[0] == '#') continue;
        v.clear();
        const char *s = line.c_str() + 2;
        for (char *end = nullptr;; s = end + (*end == ',')) {
            v.push_back(strtod(s, &end));
            if (end == s || *end == '\0') break;
        }
        if (line[0] == 'C' && v.size() == 13) {
            AnalyticalRefCase c;
            c.beta = v[1], c.alpha = v[2], c.H = v[3], c.hw = v[4], c.gammaw = v[5], c.lm = v[6], c.degenerate = v[7] != 0.;
            c.m = v[8], c.Cn = v[9], c.Dn = v[10], c.F = v[11], c.Jn = v[12];
            cases.push_back(c);
        } else if (line[0] == 'P' && v.size() == 5 && !cases.empty() && int64_t(v[0]) == int64_t(cases.size()) - 1) {
            cases.back().pts.push_back({v[1], v[2], v[3], v[4]});
        } else {
            std::cerr << "unexpected line in " << file << ": " << line.substr(0, 80) << "\n";
            return {};
        }
    }
    return cases;
}

/// C++ against Python on one reference case. The force is evaluated through the ForceField (NeoPZ coordinates) and
/// converted back, so that the coordinate conversion is part of the comparison. Excluded: points within 1e-9 R_k of
/// the discontinuity circles r = R_w, R, R_e and within 1e-10 H of the ground surface (zone or soil test decided by
/// round-off)
struct AnalyticalRefDiff {
    int64_t nforce = 0, nzero = 0, nexcl = 0; ///< compared points with f_py != 0 and f_py = 0, excluded points
    REAL frel = 0.;      ///< max |f - f_py| / |f_py| (f_py != 0)
    REAL worst[3] = {0., 0., 0.}; ///< where: r / R_w (r / R if h_w = 0), theta (deg), |f_py| / gamma_w
    REAL fzero = 0.;     ///< max |f| / gamma_w where f_py = 0
    REAL dm = 0., dF = 0., dJ = 0.; ///< |m - m_py| and the relative differences of F and J*
    bool sameFlag = true;           ///< same degenerate flag
    double build = 0., eval = 0.;   ///< construction time (s) and evaluation time per point (s)
    REAL m = 0.;
};

AnalyticalRefDiff CompareAnalyticalReference(const AnalyticalRefCase &c) {
    AnalyticalRefDiff d;
    auto t0 = std::chrono::steady_clock::now();
    auto an = std::make_shared<const AnalyticalSeepage>(c.beta, c.H, c.hw, c.alpha, 1., c.gammaw, c.lm);
    d.build = SecondsSince(t0);
    const ForceField f = AnalyticalForceField(an);
    d.m = an->M();
    d.dm = std::fabs(an->M() - c.m);
    d.dF = std::fabs(an->F() / c.F - 1.);
    d.dJ = c.Jn != 0. ? std::fabs(an->JstarNormalized() / c.Jn - 1.) : std::fabs(an->JstarNormalized());
    d.sameFlag = an->Degenerate() == c.degenerate;
    const REAL eps = 1.e-10 * c.H;
    std::vector<char> skip(c.pts.size(), 0);
    for (size_t i = 0; i < c.pts.size(); i++) {
        const REAL x = c.pts[i][0], yp = c.pts[i][1], r = std::sqrt(x * x + yp * yp);
        for (REAL Rk : {an->Rw(), an->R(), an->Re()})
            if (Rk > 0. && std::fabs(r - Rk) <= 1.e-9 * Rk) skip[i] = 1;
        const bool in = an->InSoilPaper(x, yp);
        for (int k = 0; k < 4; k++)
            if (an->InSoilPaper(x + (k & 1 ? eps : -eps), yp + (k & 2 ? eps : -eps)) != in) skip[i] = 1;
    }
    std::vector<std::array<REAL, 2>> fc(c.pts.size());
    t0 = std::chrono::steady_clock::now();
    for (size_t i = 0; i < c.pts.size(); i++) {
        TPZManVector<REAL, 3> X = {c.pts[i][0], -c.pts[i][1], 0.}; // paper -> NeoPZ
        f(X, fc[i].data());
    }
    d.eval = SecondsSince(t0) / std::max<size_t>(c.pts.size(), 1);
    for (size_t i = 0; i < c.pts.size(); i++) {
        if (skip[i]) {
            d.nexcl++;
            continue;
        }
        const REAL fx = fc[i][0], fy = -fc[i][1]; // NeoPZ -> paper
        const REAL diff = std::hypot(fx - c.pts[i][2], fy - c.pts[i][3]), ref = std::hypot(c.pts[i][2], c.pts[i][3]);
        if (ref > 0.) {
            d.nforce++;
            if (diff / ref > d.frel) {
                const REAL x = c.pts[i][0], yp = c.pts[i][1];
                d.frel = diff / ref;
                d.worst[0] = std::sqrt(x * x + yp * yp) / (an->Rw() > 0. ? an->Rw() : an->R());
                d.worst[1] = std::atan2(yp, -x) * 180. / M_PI, d.worst[2] = ref / c.gammaw;
            }
        } else {
            d.nzero++;
            d.fzero = std::max(d.fzero, std::hypot(fx, fy) / c.gammaw);
        }
    }
    return d;
}

/// table of CompareAnalyticalReference over all the cases of a reference file; returns the worst values
void AnalyticalReferenceReport(const std::vector<AnalyticalRefCase> &cases, AnalyticalRefDiff &worst, bool print) {
    worst = AnalyticalRefDiff();
    if (print)
        std::cout << "case  beta  alpha  h_w/H  H  gamma_w  m (C++)  |m - m_py|  degenerate  |dF|/F  |dJ*|/|J*|  points: f != 0 "
                     "/ f = 0 / excluded  max |df|/|f_py| at (r/R_w, theta deg; |f_py|/gw)  max |f|/gw (f_py = 0)  build (ms)  eval (ns/pt)\n";
    for (size_t k = 0; k < cases.size(); k++) {
        const AnalyticalRefCase &c = cases[k];
        const AnalyticalRefDiff d = CompareAnalyticalReference(c);
        if (print)
            std::cout << std::setprecision(6) << k << "  " << c.beta << "  " << c.alpha << "  " << c.hw / c.H << "  " << c.H << "  "
                      << c.gammaw << "  " << std::setprecision(10) << d.m << std::setprecision(3) << "  " << d.dm << "  "
                      << (c.degenerate ? "yes" : "no") << (d.sameFlag ? "" : " (C++ differs)") << "  " << d.dF << "  " << d.dJ << "  "
                      << d.nforce << " / " << d.nzero << " / " << d.nexcl << "  " << d.frel << " at (" << d.worst[0] << ", "
                      << d.worst[1] << "; " << d.worst[2] << ")  " << d.fzero << "  "
                      << 1.e3 * d.build << "  " << 1.e9 * d.eval << std::endl;
        worst.nforce += d.nforce, worst.nzero += d.nzero, worst.nexcl += d.nexcl;
        worst.frel = std::max(worst.frel, d.frel), worst.fzero = std::max(worst.fzero, d.fzero);
        worst.dm = std::max(worst.dm, d.dm), worst.dF = std::max(worst.dF, d.dF), worst.dJ = std::max(worst.dJ, d.dJ);
        worst.sameFlag = worst.sameFlag && d.sameFlag;
        worst.build = std::max(worst.build, d.build), worst.eval = std::max(worst.eval, d.eval);
    }
}

/// J*(v) of Eq. 23 by quadrature for v = K f, with f given by a ForceField in NeoPZ coordinates (K = diag(k_h, k_v)):
/// 1/2 int f.K.f over the soil inside r < R_e (polar coordinates about O split at R_w and R, 64-point Gauss in theta,
/// adaptive Gauss-Kronrod in r; zone 1 with r = R_w t^p, p = 1/m, against the r^(m-1) singularity; zone 3 with
/// r = sqrt(H^2 + X^2)) + int u^d (K f).n ds on the face and the toe ground (u^d of Eq. 21 in NeoPZ coordinates,
/// u = 0 on the crest, Cartesian outward normals, f taken 1e-11 H inside). Independent of the closed form (Eq. 40);
/// meant for 0.3 <= m or the degenerate case m = 0 (for small m > 0 the zone-1 energy converges too slowly at O).
REAL JstarByQuadrature(const ForceField &f, const AnalyticalSeepage &a) {
    const REAL kh = a.Kh(), kv = kh / a.Alpha(), gw = a.Gammaw(), H = a.H(), hw = a.Hw();
    const REAL beta = a.BetaDeg() * M_PI / 180., sb = std::sin(beta), cb = std::cos(beta), Theta = M_PI - beta;
    std::vector<REAL> xg, wg;
    la::Gauss01(64, xg, wg);
    auto ring = [&](REAL r, REAL thmax) { // int_0^thmax 1/2 f.K.f r dtheta
        REAL s = 0.;
        for (size_t i = 0; i < xg.size(); i++) {
            const REAL th = thmax * xg[i];
            TPZManVector<REAL, 3> X = {-r * std::cos(th), -r * std::sin(th), 0.};
            REAL F[2];
            f(X, F);
            s += wg[i] * 0.5 * (kh * F[0] * F[0] + kv * F[1] * F[1]);
        }
        return thmax * s * r;
    };
    const REAL eps = 1.e-11 * H, n[2] = {sb, cb}; // outward normal of the face
    auto faceFlux = [&](REAL rho) {                // (K f).n just inside the face, at the distance rho from O
        TPZManVector<REAL, 3> X = {rho * cb - eps * n[0], -rho * sb - eps * n[1], 0.};
        REAL F[2];
        f(X, F);
        return kh * F[0] * n[0] + kv * F[1] * n[1];
    };
    const REAL Rw = a.Rw(), R = a.R(), Re = a.Re(), m = a.M(), tol = 1.e-12;
    REAL J = 0.;
    if (Rw > 0.) {
        const REAL p = m > 0. && m < 1. ? 1. / m : 1.;
        J += numerics::AdaptiveGaussKronrod(
            [&](REAL t) {
                const REAL r = Rw * std::pow(t, p), jac = Rw * p * std::pow(t, p - 1.);
                return (ring(r, Theta) + gw * (-r * sb) * faceFlux(r)) * jac; // u^d = gamma_w y on the face above the water
            },
            0., 1., 0., tol);
    }
    if (R > Rw)
        J += numerics::AdaptiveGaussKronrod([&](REAL r) { return ring(r, Theta) - gw * hw * faceFlux(r); }, Rw, R, 0., tol);
    const REAL X0 = a.XToe(), X1 = std::sqrt(Re * Re - H * H);
    J += numerics::AdaptiveGaussKronrod(
        [&](REAL X) {
            const REAL r = std::hypot(H, X);
            return ring(r, M_PI - std::atan2(H, X)) * X / r;
        },
        X0, X1, 0., tol);
    J += numerics::AdaptiveGaussKronrod( // toe ground: u^d = -gamma_w h_w, outward normal (0, 1)
        [&](REAL X) {
            TPZManVector<REAL, 3> P = {X, -H - eps, 0.};
            REAL F[2];
            f(P, F);
            return -gw * hw * kv * F[1];
        },
        X0, X1, 0., tol);
    return J;
}

/// self-tests of the analytical field
void CheckAnalytical(Checker &ck, const std::string &refFile, const std::string &fig5File) {
    std::cout << "analytical field K^-1 v'_opt (AnalyticalSeepage.h)\n";
    const auto tStart = std::chrono::steady_clock::now();
    // 1) against the Python reference (scripts/analytical_seepage_reference.py: 13 cases, 300 points each)
    const std::vector<AnalyticalRefCase> cases = ReadAnalyticalReference(refFile);
    ck.Expect("    reference cases read from " + refFile.substr(refFile.find_last_of('/') + 1) + " (missing ones)",
              REAL(std::max<int64_t>(13 - int64_t(cases.size()), 0)), 0.);
    if (!cases.empty()) {
        AnalyticalRefDiff w;
        AnalyticalReferenceReport(cases, w, false);
        std::cout << "    " << cases.size() << " cases (beta 15..90, alpha 1..10, h_w/H 0..1, degenerate and near-threshold), "
                  << w.nforce << " points with f != 0, " << w.nzero << " with f = 0, " << w.nexcl << " excluded\n";
        ck.Expect("    vs Python: max |f - f_py| / |f_py| (NeoPZ ForceField converted back to the paper coordinates)", w.frel, 1.e-6);
        ck.Expect("    vs Python: max |f| / gamma_w where f_py = 0 (outside the soil, r >= R_e, h_w = 0)", w.fzero, 0.);
        ck.Expect("    vs Python: max |m - m_py|", w.dm, 1.e-8);
        ck.Expect("    vs Python: max relative difference of F = h2e^2 / (sqrt C + sqrt D)^2 and of J*",
                  std::max(w.dF, w.dJ), 1.e-10);
        ck.Expect("    vs Python: cases with a different degenerate flag", w.sameFlag ? 0. : 1., 0.);
    }
    // 2) Fig. 5 solid curves (lower edge of the bands of data/fig5_vector_fill_polygons.csv), h_w = H, L_m = 10 H
    {
        const auto paper = ReadFig5(fig5File, true);
        REAL worst = 0.;
        int n = 0;
        for (const auto &kv : paper) {
            const AnalyticalSeepage an(kv.first.second, 1., 1., kv.first.first);
            worst = std::max(worst, std::fabs(-an.JstarNormalized() / kv.second - 1.));
            n++;
        }
        std::cout << "    Fig. 5: " << n << " points (alpha 1, 2, 4, 10; beta 15..90 every 5 deg)\n";
        ck.Expect("    Fig. 5 solid: points read (missing of 64)", REAL(std::max(64 - n, 0)), 0.);
        ck.Expect("    Fig. 5 solid: max |-J*/(kh H^2 gw^2) / paper - 1|", worst, 5.e-4);
    }
    // 3) degenerate optimum m -> 0 exactly where Phi(m) = A + m (A - P) + O(m^2) increases, A > P with
    //    P = int_0^Theta (int_0^t d)^2 / c dt (thresholds 80.761 / 81.101 / 88.398 deg for alpha = 1 / 2 / 4)
    {
        int bad = 0, ndeg = 0, ncases = 0;
        std::vector<std::pair<REAL, REAL>> ab;
        for (REAL alpha : {1., 2., 4., 5., 10.})
            for (REAL beta = 15.; beta <= 90.; beta += 5.) ab.push_back({alpha, beta});
        for (auto t : {std::make_pair(1., 80.71), std::make_pair(1., 80.81), std::make_pair(2., 81.05), std::make_pair(2., 81.15),
                       std::make_pair(4., 88.35), std::make_pair(4., 88.45)})
            ab.push_back(t);
        for (auto &t : ab) {
            const REAL alpha = t.first, beta = t.second, Theta = M_PI - beta * M_PI / 180.;
            const AnalyticalSeepage an(beta, 1., 1., alpha);
            const REAL P = numerics::AdaptiveGaussKronrod(
                [alpha](REAL s) {
                    const REAL Id = (alpha + 1.) * s / 2. + (alpha - 1.) * std::sin(2. * s) / 4.; // int_0^s d
                    return Id * Id / (std::cos(s) * std::cos(s) + alpha * std::sin(s) * std::sin(s));
                },
                0., Theta, 0., 1.e-13);
            if (an.Degenerate() != (an.A() > P)) bad++;
            ndeg += an.Degenerate(), ncases++;
        }
        std::cout << "    degenerate optimum: " << ndeg << " of " << ncases << " cases\n";
        ck.Expect("    degenerate flag different from the criterion A > P (count)", REAL(bad), 0.);
    }
    // 4) J* by quadrature of the field (through the NeoPZ ForceField) against the closed form of Eq. 40
    for (auto c : {std::array<REAL, 3>{30., 1., 1.}, std::array<REAL, 3>{60., 4., 0.6}, std::array<REAL, 3>{45., 10., 0.3},
                   std::array<REAL, 3>{20., 2., 0.5}, std::array<REAL, 3>{90., 10., 1.}, std::array<REAL, 3>{85., 1., 1.},
                   std::array<REAL, 3>{90., 1., 0.4}}) {
        auto an = std::make_shared<const AnalyticalSeepage>(c[0], 5., c[2] * 5., c[1], 2., 9.81);
        const REAL Jq = JstarByQuadrature(AnalyticalForceField(an), *an);
        std::ostringstream os;
        os << "    J* by quadrature vs Eq. 40, beta " << c[0] << " alpha " << c[1] << " h_w/H " << c[2] << " (m " << std::setprecision(4)
           << an->M() << "): relative difference";
        ck.Expect(os.str(), std::fabs(Jq / an->Jstar() - 1.), 1.e-9);
    }
    // 5) admissibility: div v = div(K f) = 0 (central differences, NeoPZ coordinates) at random points of zones 1-3,
    //    normal flux v.e_r continuous across r = R_w and R (zero on both sides), h2 interpolation and natural BC
    {
        std::mt19937_64 rng(7);
        std::uniform_real_distribution<REAL> U(0., 1.);
        REAL divmax = 0., jump = 0., interp = 0., ident = 0.;
        for (auto c : {std::array<REAL, 3>{30., 1., 1.}, std::array<REAL, 3>{60., 4., 0.6}, std::array<REAL, 3>{45., 10., 0.3},
                       std::array<REAL, 3>{90., 2., 1.}, std::array<REAL, 3>{20., 2., 0.5}, std::array<REAL, 3>{88., 4., 0.5},
                       std::array<REAL, 3>{85., 1., 0.7}}) {
            auto an = std::make_shared<const AnalyticalSeepage>(c[0], 1., c[2], c[1], 1., 9.81);
            const ForceField f = AnalyticalForceField(an);
            const REAL kh = an->Kh(), kv = kh / an->Alpha();
            auto v = [&](REAL x, REAL y, REAL V[2]) {
                TPZManVector<REAL, 3> X = {x, y, 0.};
                REAL F[2];
                f(X, F);
                V[0] = kh * F[0], V[1] = kv * F[1];
            };
            for (int i = 0; i < 2000; i++) {
                const REAL r = an->Rw() * 1.e-3 * std::pow(an->Re() / (an->Rw() * 1.e-3), U(rng));
                const REAL thmax = r < an->R() ? an->Theta() : M_PI - std::asin(an->H() / r);
                const REAL th = thmax * (0.01 + 0.98 * U(rng));
                bool nearCircle = false;
                for (REAL Rk : {an->Rw(), an->R(), an->Re()})
                    if (std::fabs(r - Rk) < 1.e-3 * r) nearCircle = true;
                const REAL x = -r * std::cos(th), y = -r * std::sin(th), h = 1.e-6 * r; // B(r) ~ sqrt(r - H) near R = H
                if (nearCircle || std::fabs(y + an->H()) < 1.e-3 * r) continue;
                REAL V[2], a[2], b[2], cc[2], d[2];
                v(x, y, V), v(x + h, y, a), v(x - h, y, b), v(x, y + h, cc), v(x, y - h, d);
                const REAL div = (a[0] - b[0] + cc[1] - d[1]) / (2. * h);
                divmax = std::max(divmax, std::fabs(div) / (std::hypot(V[0], V[1]) / r));
            }
            for (int i = 1; i < 50; i++) { // v_r on both sides of R_w and R
                const REAL th = an->Theta() * i / 50.;
                REAL vr1, vt1, vr2, vt2;
                for (REAL Rk : {an->Rw(), an->R()}) {
                    an->PolarVelocity(Rk * (1. - 1.e-12), th, vr1, vt1);
                    an->PolarVelocity(Rk * (1. + 1.e-12), th, vr2, vt2);
                    jump = std::max(jump, std::fabs(vr1 - vr2) / (std::fabs(vt1) + std::fabs(vt2)));
                }
            }
            if (!an->Degenerate()) {
                REAL y[4];
                for (int i = 0; i < 7; i++) { // h2 at angles between the nodes against a direct integration
                    const REAL th = an->Theta() * (0.0371 + 0.137 * i);
                    an->Shoot(an->M(), y, 1.e-13, th);
                    REAL h2, dh2;
                    an->H2(th, h2, dh2);
                    const REAL s = std::sin(th), c2 = std::cos(th) * std::cos(th), cth = c2 + an->Alpha() * s * s;
                    interp = std::max({interp, std::fabs(h2 / y[0] - 1.), std::fabs(dh2 - an->M() * y[1] / cth) / std::fabs(y[1] / cth)});
                }
                // natural condition at the face, an identity at the optimum: c h2' h2 (Theta) = sqrt C (sqrt C + sqrt D)
                REAL h2, dh2;
                an->H2(an->Theta(), h2, dh2);
                const REAL ce = std::cos(an->Theta()) * std::cos(an->Theta()) + an->Alpha() * std::pow(std::sin(an->Theta()), 2);
                ident = std::max(ident, std::fabs(ce * dh2 * h2 - std::sqrt(an->C()) * (std::sqrt(an->C()) + std::sqrt(an->D()))) / an->C());
            }
        }
        ck.Expect("    max |div(K f)| / (|K f| / r) by central differences (zones 1-3)", divmax, 1.e-6);
        ck.Expect("    jump of v.e_r across r = R_w and R / |v_theta|", jump, 1.e-9);
        ck.Expect("    h2, h2' interpolated vs integrated at angles between the nodes (relative)", interp, 1.e-10);
        ck.Expect("    natural condition at the face c h2' h2 = sqrt C (sqrt C + sqrt D) (relative to C)", ident, 1.e-9);
    }
    // 6) zero outside the soil and for r >= R_e; h_w = 0; evaluator from 4 threads bit-identical to the serial one
    {
        auto an = std::make_shared<const AnalyticalSeepage>(60., 5., 3., 4., 1., 9.81);
        const ForceField f = AnalyticalForceField(an);
        std::mt19937_64 rng(11);
        std::uniform_real_distribution<REAL> U(0., 1.);
        REAL fout = 0.;
        int nin = 0;
        const AnalyticalSeepage dry(60., 5., 0., 4.);
        for (int i = 0; i < 20000; i++) {
            const REAL r = an->Re() * (0.01 + 2. * U(rng)), th = 2. * M_PI * U(rng);
            TPZManVector<REAL, 3> X = {r * std::cos(th), r * std::sin(th), 0.};
            REAL F[2], F0[2];
            f(X, F);
            dry.Force(X, F0);
            fout = std::max({fout, std::fabs(F0[0]), std::fabs(F0[1])});
            const bool soil = an->InSoilPaper(X[0], -X[1]) && r < an->Re();
            if (!soil) fout = std::max({fout, std::fabs(F[0]), std::fabs(F[1])});
            else nin += std::hypot(F[0], F[1]) > 0.;
        }
        ck.Expect("    f = 0 outside the soil, for r >= R_e and for h_w = 0 (max |f|)", fout, 0.);
        ck.Expect("    f != 0 inside the soil for r < R_e (points found)", nin > 1000 ? 0. : 1., 0.);
        const int np = 100000, nt = 4;
        std::vector<REAL> xs(2 * np), f1(2 * np), f4(2 * np);
        for (int i = 0; i < np; i++) xs[2 * i] = -3. + 9. * U(rng), xs[2 * i + 1] = -8. * U(rng);
        auto run = [&](std::vector<REAL> &out, int k0, int k1) {
            for (int i = k0; i < k1; i++) {
                TPZManVector<REAL, 3> X = {xs[2 * i], xs[2 * i + 1], 0.};
                f(X, &out[2 * i]);
            }
        };
        run(f1, 0, np);
        std::vector<std::thread> th;
        for (int t = 0; t < nt; t++) th.emplace_back(run, std::ref(f4), t * np / nt, (t + 1) * np / nt);
        for (auto &t : th) t.join();
        REAL diff = 0.;
        for (int i = 0; i < 2 * np; i++) diff = std::max(diff, std::fabs(f1[i] - f4[i]));
        ck.Expect("    evaluator from 4 threads vs serial, 1e5 points (max |difference|)", diff, 0.);
    }
    std::cout << "    (analytical field tests: " << SecondsSince(tStart) << " s)\n";
}

} // namespace

int AnalyticalReferenceTable(const std::string &ref) {
    auto t0 = std::chrono::steady_clock::now();
    const std::vector<AnalyticalRefCase> cases = ReadAnalyticalReference(ref);
    int64_t n = 0;
    for (auto &c : cases) n += int64_t(c.pts.size());
    std::cout << "reference " << ref << ": " << cases.size() << " cases, " << n << " points (read in " << SecondsSince(t0) << " s)\n";
    if (cases.empty()) return 1;
    AnalyticalRefDiff w;
    AnalyticalReferenceReport(cases, w, true);
    std::cout << std::setprecision(3) << "all: " << w.nforce << " points with f != 0 (max |f - f_py| / |f_py| " << w.frel << "), "
              << w.nzero << " with f = 0 (max |f| / gw " << w.fzero << "), " << w.nexcl << " excluded; max |m - m_py| " << w.dm
              << ", max relative difference of F " << w.dF << " and of J* " << w.dJ << "; degenerate flags "
              << (w.sameFlag ? "equal" : "DIFFERENT") << "\n";
    return w.frel <= 1.e-6 && w.fzero == 0. && w.sameFlag ? 0 : 1;
}

int CmdVerify(Options &o) {
    if (!o.Has("alpha")) std::cout << "(verify: alpha = 5 by default)\n";
    const REAL alpha = o.Num("alpha", 5.);
    Problem p = ReadProblem(o);
    p.hy.kh = alpha * p.hy.kv;
    const int href = o.IntList("href", {0})[0];
    o.CheckUnused();
    REAL err[2][3];
    const REAL worst = VerifyManufactured(p, href, true, err);
    std::cout << "max relative error " << worst << "\n";
    return worst < 1.e-9 ? 0 : 1;
}

/// check: the component self-tests above (no nonlinear analysis), then those of the figure drivers and of fembatch
/// (tiny limit analyses and elastoplastic runs); about 20-30 s
int CmdCheck(Options &o) {
    Problem p0 = ReadProblem(o);
    o.CheckUnused();
    const auto t0 = std::chrono::steady_clock::now();
    Checker ck;
    CheckManufactured(ck, p0);
    CheckSeepageField(ck, p0, 45., 1., 1., "zero_lb");
    CheckSeepageField(ck, p0, 30., 0.5, 4., "zero_lb");
    CheckSeepageField(ck, p0, 60., 0.25, 2., "toe_r");
    CheckSeepageField(ck, p0, 90., 0.6, 10., "impermeable");
    CheckSeepageField(ck, p0, 15., 1., 5., "zero_b");
    // same mesh for u (P2) and the stability (P2, exact integration of the linear f): divergence theorem to round-off
    CheckLoadVector(ck, p0, 45., 1., 1., true, 1.e-11);
    CheckLoadVector(ck, p0, 45., 0.5, 4., true, 1.e-11);
    // generated meshes: f is discontinuous inside the stability elements (quadrature error only)
    CheckLoadVector(ck, p0, 30., 1., 5., false, 2.e-3);
    CheckLoadVector(ck, p0, 75., 0.4, 1., false, 2.e-3);
    CheckLimitAnalysis(ck, p0);
    CheckAnalytical(ck, ProjectDir() + "data/analytical_seepage_reference.csv", ProjectDir() + "data/fig5_vector_fill_polygons.csv");
    CheckFigures(ck, p0);
    CheckFEMBatch(ck, p0);
    std::cout << "\ncheck: " << ck.total - ck.failed << " of " << ck.total << " passed, " << SecondsSince(t0) << " s\n";
    return ck.failed == 0 ? 0 : 1;
}

} // namespace slope
