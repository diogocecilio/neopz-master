// Slope stability under the seepage forces of a rapid drawdown (Ceron, Cecilio, Linn & Maghous, IJNAMG 2025,
// doi:10.1002/nag.3993, Section 4.3), built on Projects/SlopeDrawdown and Projects/SlopeMohrCoulomb:
//  1. parametric slope meshes (SlopeGeometry.h, DelaunayMesher.h): large hydraulic domain and stability domain;
//  2. steady anisotropic seepage of the excess pore pressure u after the drawdown (SeepageFE.h,
//     AnisotropicDarcy.h): J(u_FE) = 1/2 int grad u . K grad u;
//  3. seepage force field f = -grad u_FE at any point (SeepageForceField.h);
//  4. FEM stability: gravity increase of SlopeAnalysis.h with b = lambda (gamma' g + f) (FEMStability.h);
//     Gamma_FEM = lambda_crit, H_crit = lambda_crit H.
// The limit analysis and the analytical field K^-1 v'_opt (scripts/limit_analysis.py, scripts/analytical_seepage.py)
// are to be ported later: they plug in through slope::ForceField.
//
// Usage: SlopeSeepageForces <command> [key=value ...]
//   mesh     meshes, quality and consistency (vtk=1 writes them; sweep=1 checks beta = 15..90, h_w / H = 0..1)
//   verify   manufactured solutions of the anisotropic seepage solver (exact for P2) and point location
//   seepage  drawdown seepage: J, J / (k_h H^2 gamma_w^2), min p, for href=<list> uniform refinements
//   fig5     J / (k_h H^2 gamma_w^2) for alphas=<list> betas=<list> against the dashed curves of the paper's Fig. 5
//   probe    u and f = -grad u at the points of pts=<file> (lines x_paper,y_paper), paper coordinates
//   fs       gravity-increase factor of the stability domain with the seepage forces (or dry)
//   check    self-tests (about 6 s, exit code 1 on failure): manufactured solutions and orientation of K; Eq. 21
//            and the far-side data on the solved field, boundary ids, p >= 0, f = -grad u, point location on
//            vertices / edges / just outside, continuity of the P2 field, evaluator from 4 threads; load vector of
//            the stability problem: resultant = -int u n ds - gamma' A e_y (divergence theorem), lambda scales
//            gamma' and f, forms u and p identical, threaded assembly = serial
// Options (default):
//   slope/soil: H=5 beta=45 hw=1 (h_w / H) gamma=20 gammaw=9.81 c=10 phi=30 E=20000 nu=0.3
//   seepage:    alpha=1 (k_h / k_v) kv=1 horder=2 hbc=zero_lb|impermeable|zero_b|zero_l|zero_lbr|toe_r (far sides,
//               presets of scripts/fe_seepage.py; zero_lb: u = 0 left and base, see SeepageFE.h) hbcleft= hbcbottom=
//               hbcright= (noflow|zero|toe, override of the preset)
//   hydraulic mesh: hmesh=gen|trig hleft=50 hright=10 hdepth=30 (units of H, from O, T, T) hh0=0.025 hhs=0.0625
//               hgrade=0.15 hhmax=2 (sizes in units of H) href=0 (list for seepage, e.g. href=0,1,2)
//   stability mesh: smesh=gen|trig sa=2 (extents sa H + H / tan(beta)) sleft= sright= sdepth= (override, units of
//               H) sh0=0.25 shs=0.25 sgrade=0.25 shmax=1 sref=0
//   fs:         water=seepage|dry form=u|p|p+ (b = lambda (gamma' g - grad u) | lambda (gamma_sat g - grad p) |
//               lambda (gamma_sat g - grad max(p, 0)) as SlopeDrawdown) nref=3 srm=0 vtk=<prefix>
//               mark=0.1 (refinement of the elements with sqrt(J2(eps_p)) >= mark * max at collapse)
//               maxnewton=100 tolfs=0.002 (driver: Newton iterations, relative step of the continuation;
//               maxnewton=30 is the setting of SlopeMohrCoulomb / SlopeDrawdown, see FEMStability.h)
//               checkforms=1 (compares the load vectors of the forms u and p)
//   fig5:       alphas=1,2,4,10 betas=15,30,45,60,75,90 paper=<csv> (default data/fig5_vector_fill_polygons.csv)
//   output:     seepage: csv=<file> (u and f on a grid around the slope, paper coordinates), vtk=<file> (u, -grad u
//               and the Darcy velocity on the hydraulic mesh, NeoPZ coordinates), both for the last href
// mesh=trig: TriGMesh(1 + ref) of SlopeMohrCoulomb (H = 10 m, beta = 45 deg, 70 x 40 m), the mesh of SlopeDrawdown.
#include "FEMStability.h"
#include "SeepageFE.h"

#include "TPZVTKGeoMesh.h"

#include <chrono>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

using namespace slope;

/// key=value options; unknown keys are an error
class Options {
    std::map<std::string, std::string> fKV;
    std::set<std::string> fUsed;

public:
    Options(int argc, char *argv[], int first) {
        for (int i = first; i < argc; i++) {
            const char *eq = strchr(argv[i], '=');
            if (!eq) {
                std::cerr << "option without '=': " << argv[i] << "\n";
                exit(1);
            }
            fKV[std::string(argv[i], size_t(eq - argv[i]))] = eq + 1;
        }
    }
    bool Has(const std::string &k) const { return fKV.count(k) > 0; }
    std::string Str(const std::string &k, const std::string &def) {
        fUsed.insert(k);
        return Has(k) ? fKV.at(k) : def;
    }
    double Num(const std::string &k, double def) {
        fUsed.insert(k);
        return Has(k) ? atof(fKV.at(k).c_str()) : def;
    }
    std::vector<double> NumList(const std::string &k, const std::vector<double> &def) {
        fUsed.insert(k);
        if (!Has(k)) return def;
        std::vector<double> v;
        std::stringstream ss(fKV.at(k));
        std::string item;
        while (std::getline(ss, item, ',')) v.push_back(atof(item.c_str()));
        if (v.empty()) {
            std::cerr << "empty list for option " << k << "\n";
            exit(1);
        }
        return v;
    }
    std::vector<int> IntList(const std::string &k, const std::vector<int> &def) {
        std::vector<int> v;
        for (double d : NumList(k, std::vector<double>(def.begin(), def.end()))) v.push_back(int(d));
        return v;
    }
    void CheckUnused() const {
        for (auto &kv : fKV)
            if (!fUsed.count(kv.first)) {
                std::cerr << "unknown option " << kv.first << "\n";
                exit(1);
            }
    }
};

struct Problem {
    SlopeGeometry hyd, stab; ///< hydraulic and stability domains
    MeshSize hsize, ssize;
    std::string hmesh, smesh;
    Hydraulics hy;
    Soil soil;
    int sref = 0;
    REAL hext[3] = {50., 10., 30.}; ///< hydraulic box / H: left of O, right of T, below T
    REAL sa = 2.;                   ///< stability box: sa H + H / tan(beta) on each side ...
    REAL sext[3] = {-1., -1., -1.}; ///< ... unless overridden (/ H, > 0): left of O, right of T, below T

    /// both domains for the slope angle beta (deg) and the water level h_w = hwr H
    void SetSlope(REAL beta, REAL hwr) {
        if ((hmesh == "trig" || smesh == "trig") && (std::fabs(beta - 45.) > 1.e-12 || std::fabs(hyd.H - 10.) > 1.e-12)) {
            std::cerr << "mesh=trig is the slope of SlopeMohrCoulomb: H = 10, beta = 45 only (asked: H " << hyd.H << ", beta "
                      << beta << ")\n";
            exit(1);
        }
        SlopeGeometry g;
        g.H = hyd.H, g.beta = beta, g.hw = hwr * g.H;
        hyd = g;
        hyd.left = hext[0] * g.H, hyd.right = hext[1] * g.H, hyd.depth = hext[2] * g.H;
        if (hmesh == "trig") hyd = SlopeMohrCoulombGeometry(g.hw);
        stab = g;
        StabilityExtents(stab, sa);
        REAL *ext[3] = {&stab.left, &stab.right, &stab.depth};
        for (int i = 0; i < 3; i++)
            if (sext[i] > 0.) *ext[i] = sext[i] * g.H;
        if (smesh == "trig") stab = SlopeMohrCoulombGeometry(g.hw);
    }
};

Problem ReadProblem(Options &o) {
    Problem p;
    p.hyd.H = o.Num("H", 5.);
    const REAL beta = o.Num("beta", 45.), hwr = o.Num("hw", 1.);
    p.hmesh = o.Str("hmesh", "gen");
    p.smesh = o.Str("smesh", "gen");
    for (const std::string &m : {p.hmesh, p.smesh})
        if (m != "gen" && m != "trig") {
            std::cerr << "hmesh and smesh must be gen or trig\n";
            exit(1);
        }
    if (p.hmesh == "trig" || p.smesh == "trig") {
        if (std::fabs(p.hyd.H - 10.) > 1.e-12 || std::fabs(beta - 45.) > 1.e-12) {
            std::cerr << "mesh=trig is the slope of SlopeMohrCoulomb: use H=10 beta=45\n";
            exit(1);
        }
    }
    const char *hkeys[3] = {"hleft", "hright", "hdepth"}, *skeys[3] = {"sleft", "sright", "sdepth"};
    for (int i = 0; i < 3; i++) {
        p.hext[i] = o.Num(hkeys[i], p.hext[i]);
        p.sext[i] = o.Num(skeys[i], -1.);
    }
    p.sa = o.Num("sa", 2.);
    p.SetSlope(beta, hwr);
    p.hsize.h0 = o.Num("hh0", 0.025);
    p.hsize.hs = o.Num("hhs", 0.0625);
    p.hsize.grade = o.Num("hgrade", 0.15);
    p.hsize.hmax = o.Num("hhmax", 2.);
    p.ssize.h0 = o.Num("sh0", 0.25);
    p.ssize.hs = o.Num("shs", 0.25);
    p.ssize.grade = o.Num("sgrade", 0.25);
    p.ssize.hmax = o.Num("shmax", 1.);
    p.sref = int(o.Num("sref", 0));
    p.hy.kv = o.Num("kv", 1.);
    p.hy.kh = o.Num("alpha", 1.) * p.hy.kv;
    p.hy.gammaw = o.Num("gammaw", 9.81);
    p.hy.order = int(o.Num("horder", 2));
    if (p.hy.order < 1 || p.hy.order > 2) {
        std::cerr << "horder must be 1 or 2\n";
        exit(1);
    }
    const std::string preset = o.Str("hbc", "zero_lb"); // see SeepageFE.h
    if (!p.hy.far.SetPreset(preset)) {
        std::cerr << "unknown hbc preset " << preset << " (impermeable, zero_lb, zero_b, zero_l, zero_lbr, toe_r)\n";
        exit(1);
    }
    const std::pair<const char *, EFarBC *> sides[3] = {
        {"hbcleft", &p.hy.far.left}, {"hbcbottom", &p.hy.far.bottom}, {"hbcright", &p.hy.far.right}};
    for (auto &side : sides)
        if (o.Has(side.first) && !FarSides::Parse(o.Str(side.first, ""), *side.second)) {
            std::cerr << side.first << " must be noflow, zero or toe\n";
            exit(1);
        }
    p.soil.gamma = o.Num("gamma", 20.);
    p.soil.c = o.Num("c", 10.);
    p.soil.phi = o.Num("phi", 30.) * M_PI / 180.;
    p.soil.E = o.Num("E", 20000.);
    p.soil.nu = o.Num("nu", 0.3);
    return p;
}

TPZGeoMesh *HydraulicMesh(const Problem &p, int ref, MeshStats *st = nullptr) {
    if (p.hmesh == "trig") return SlopeMohrCoulombGMesh(1 + ref);
    if (p.hmesh != "gen") {
        std::cerr << "hmesh must be gen or trig\n";
        exit(1);
    }
    return CreateSlopeGMesh(p.hyd, p.hsize, ref, st);
}

TPZGeoMesh *StabilityMesh(const Problem &p, MeshStats *st = nullptr) {
    if (p.smesh == "trig") return SlopeMohrCoulombGMesh(1 + p.sref);
    if (p.smesh != "gen") {
        std::cerr << "smesh must be gen or trig\n";
        exit(1);
    }
    return CreateSlopeGMesh(p.stab, p.ssize, p.sref, st);
}

void PrintGeometry(const char *name, const SlopeGeometry &g) {
    std::cout << name << ": H " << g.H << " m, beta " << g.beta << " deg, h_w " << g.hw << " m, box x in [" << -g.left
              << ", " << g.XT() + g.right << "], y in [" << -g.H - g.depth << ", 0]\n";
}

void PrintHydraulics(const Hydraulics &hy, const std::string &mesh) {
    std::cout << "K = diag(" << hy.kh << ", " << hy.kv << "), alpha " << hy.kh / hy.kv << ", gamma_w " << hy.gammaw
              << ", order " << hy.order << ", mesh " << mesh << "; far sides: left " << FarSides::Name(hy.far.left)
              << ", base " << FarSides::Name(hy.far.bottom) << ", right " << FarSides::Name(hy.far.right) << "\n";
}

int64_t CountLeafTriangles(TPZGeoMesh *gmesh) {
    int64_t n = 0;
    for (int64_t i = 0; i < gmesh->NElements(); i++) {
        TPZGeoEl *gel = gmesh->Element(i);
        if (gel && !gel->HasSubElement() && gel->Dimension() == 2) n++;
    }
    return n;
}

/// mesh sweep=1: both generated domains for beta = 15 .. 90 deg and h_w / H = 0 .. 1, consistency and quality
int MeshSweep(Problem p0) {
    p0.hmesh = p0.smesh = "gen";
    const auto t0 = std::chrono::steady_clock::now();
    REAL worst = 0., minAngle = 180.;
    int64_t maxTri = 0, n = 0;
    for (REAL beta = 15.; beta <= 90.; beta += 5.)
        for (REAL hw : {0., 0.1, 0.25, 0.5, 0.75, 0.9, 1.}) {
            Problem p = p0;
            p.SetSlope(beta, hw);
            for (int k = 0; k < 2; k++) {
                MeshStats st;
                const SlopeGeometry &g = k == 0 ? p.hyd : p.stab;
                TPZGeoMesh *gmesh = CreateSlopeGMesh(g, k == 0 ? p.hsize : p.ssize, 0, &st);
                const REAL err = CheckSlopeGMesh(gmesh, g);
                if (err > 1.e-10)
                    std::cout << "  mesh check failed: beta " << beta << " h_w/H " << hw << " domain " << k << " error " << err << "\n";
                worst = std::max(worst, err), minAngle = std::min(minAngle, st.minAngle);
                maxTri = std::max(maxTri, st.triangles), n++;
                delete gmesh;
            }
        }
    std::cout << n << " meshes: max relative error of the area / boundary lengths " << worst << ", min angle " << minAngle
              << " deg, max triangles " << maxTri << ", "
              << std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count() << " s\n";
    return worst < 1.e-10 ? 0 : 1;
}

int CmdMesh(Options &o) {
    Problem p = ReadProblem(o);
    const bool vtk = o.Num("vtk", 0) != 0.;
    const int href = o.IntList("href", {0})[0];
    const bool sweep = o.Num("sweep", 0) != 0.;
    o.CheckUnused();
    if (sweep) return MeshSweep(p);
    for (int k = 0; k < 2; k++) {
        MeshStats st;
        const auto t0 = std::chrono::steady_clock::now();
        TPZGeoMesh *gmesh = k == 0 ? HydraulicMesh(p, href, &st) : StabilityMesh(p, &st);
        const double dt = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        PrintGeometry(k == 0 ? "hydraulic domain" : "stability domain", k == 0 ? p.hyd : p.stab);
        std::cout << "  leaf triangles " << CountLeafTriangles(gmesh);
        if (st.triangles) std::cout << " (base mesh: " << st.nodes << " nodes, " << st.triangles << " triangles, min angle "
                                    << st.minAngle << " deg, edges " << st.hmin << " .. " << st.hmax << " m)";
        std::cout << ", " << dt << " s\n";
        CheckSlopeGMesh(gmesh, k == 0 ? p.hyd : p.stab, true);
        if (vtk) {
            std::ofstream f(k == 0 ? "mesh_hydraulic.vtk" : "mesh_stability.vtk");
            TPZVTKGeoMesh::PrintGMeshVTK(gmesh, f, true);
        }
        delete gmesh;
    }
    return 0;
}

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
        const double dt = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
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
        const double tl = std::chrono::duration<double>(std::chrono::steady_clock::now() - t1).count();
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

/// u and f = -grad u on a grid around the slope, paper coordinates (y down): x, y, u, fx, fy
void WriteGridCSV(const std::string &file, const PoreField &pf, const SlopeGeometry &g) {
    std::ofstream f(file);
    f << "# x_paper,y_paper,u_kPa,fx_paper,fy_paper (f = -grad u, kN/m^3); H " << g.H << " beta " << g.beta << " hw " << g.hw
      << "\n";
    f << std::setprecision(10);
    const REAL d = g.H / 20.;
    for (REAL x = -2. * g.H; x <= g.XT() + 2. * g.H + 1.e-9; x += d)
        for (REAL yp = 0.; yp <= 2. * g.H + 1.e-9; yp += d) {
            TPZManVector<REAL, 3> X = {x, -yp, 0.};
            REAL u, gr[2];
            if (!pf.Evaluate(X, u, gr)) continue;
            f << x << "," << yp << "," << u << "," << -gr[0] << "," << gr[1] << "\n"; // fy_paper = -fy = grad_y u
        }
}

int CmdSeepage(Options &o) {
    Problem p = ReadProblem(o);
    const std::vector<int> refs = o.IntList("href", {0});
    const std::string csv = o.Str("csv", ""), vtk = o.Str("vtk", "");
    o.CheckUnused();
    PrintGeometry("hydraulic domain", p.hyd);
    PrintHydraulics(p.hy, p.hmesh);
    std::vector<REAL> J;
    std::cout << "href  triangles  equations  J  J/(kh H^2 gw^2)  min p (kPa) at (x, y)  time (s)\n";
    for (int ref : refs) {
        TPZGeoMesh *gmesh = HydraulicMesh(p, ref);
        SeepageResult r = DrawdownSeepage(gmesh, p.hyd, p.hy, ref == refs.back() ? vtk : "");
        std::cout << std::setprecision(10) << ref << "  " << r.nel << "  " << r.neq << "  " << r.J << "  " << r.Jnorm
                  << std::setprecision(4) << "  " << r.pmin << " at (" << r.xpmin[0] << ", " << r.xpmin[1] << ")  " << r.seconds
                  << "\n";
        J.push_back(r.Jnorm);
        if (!csv.empty() && ref == refs.back()) WriteGridCSV(csv, *r.field, p.hyd);
        delete gmesh;
    }
    if (J.size() >= 3) { // nested uniform refinements: observed order in h and Richardson extrapolation
        std::cout << std::setprecision(8);
        for (size_t k = 2; k < J.size(); k++) {
            const REAL d1 = J[k - 1] - J[k - 2], d2 = J[k] - J[k - 1];
            const REAL q = std::log(std::fabs(d1 / d2)) / std::log(2.);
            const REAL Jinf = J[k] + d2 / (std::pow(2., q) - 1.);
            std::cout << "href " << refs[k - 2] << ".." << refs[k] << ": order " << q << ", extrapolated J/(kh H^2 gw^2) " << Jinf
                      << ", error of the last " << (J[k] - Jinf) / Jinf << "\n";
        }
    }
    return 0;
}

/// u and f = -grad u at points given in paper coordinates (file: lines x_paper,y_paper; '#' comments); output in
/// paper coordinates: x, y, u, fx, fy (fy_paper = -fy)
int CmdProbe(Options &o) {
    Problem p = ReadProblem(o);
    const std::string file = o.Str("pts", "");
    const int href = o.IntList("href", {0})[0];
    o.CheckUnused();
    std::ifstream in(file);
    if (!in) {
        std::cerr << "probe: cannot read pts=" << file << "\n";
        return 1;
    }
    TPZGeoMesh *gmesh = HydraulicMesh(p, href);
    SeepageResult r = DrawdownSeepage(gmesh, p.hyd, p.hy);
    delete gmesh;
    PrintGeometry("hydraulic domain", p.hyd);
    PrintHydraulics(p.hy, p.hmesh);
    std::cout << std::setprecision(10) << "J/(kh H^2 gw^2) " << r.Jnorm << "\nx_paper,y_paper,u,fx_paper,fy_paper\n";
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        REAL xp, yp;
        if (sscanf(line.c_str(), "%lf,%lf", &xp, &yp) != 2) continue;
        TPZManVector<REAL, 3> x = {xp, -yp, 0.};
        REAL u, g[2];
        if (!r.field->Evaluate(x, u, g)) {
            std::cout << xp << "," << yp << ",nan,0,0\n";
            continue;
        }
        std::cout << xp << "," << yp << "," << u << "," << -g[0] << "," << g[1] << "\n";
    }
    return 0;
}

/// Fig. 5 of the paper (vector data extracted from the PDF): alpha, beta -> J(u'_FE) / (k_h H^2 gamma_w^2), dashed
std::map<std::pair<int, int>, REAL> ReadFig5(const std::string &file) {
    std::map<std::pair<int, int>, REAL> m;
    std::ifstream f(file);
    std::string line;
    while (std::getline(f, line)) {
        if (line.empty() || line[0] == '#') continue;
        REAL a, b, solid, dashed;
        if (sscanf(line.c_str(), "%lf,%lf,%lf,%lf", &a, &b, &solid, &dashed) == 4)
            m[{int(std::lround(a)), int(std::lround(b))}] = dashed;
    }
    return m;
}

/// J(u_FE) / (k_h H^2 gamma_w^2) versus beta and alpha (paper Fig. 5, dashed curves) on the hydraulic domain
int CmdFig5(Options &o) {
    Problem p = ReadProblem(o);
    const std::vector<double> alphas = o.NumList("alphas", {1., 2., 4., 10.});
    const std::vector<double> betas = o.NumList("betas", {15., 30., 45., 60., 75., 90.});
    const int href = o.IntList("href", {0})[0];
    const std::string src(__FILE__);
    const std::string file = o.Str("paper", src.substr(0, src.find_last_of('/') + 1) + "data/fig5_vector_fill_polygons.csv");
    const REAL hwr = p.hyd.hw / p.hyd.H;
    o.CheckUnused();
    const auto paper = ReadFig5(file);
    if (paper.empty()) std::cout << "(paper data not found: " << file << ")\n";
    PrintGeometry("hydraulic domain", p.hyd);
    PrintHydraulics(p.hy, p.hmesh);
    std::cout << "h_w / H " << hwr << ", href " << href << "\nalpha  beta  triangles  equations  J/(kh H^2 gw^2)  paper (dashed)"
              << "  difference  min p (kPa)  time (s)\n";
    for (double alpha : alphas)
        for (double beta : betas) {
            p.SetSlope(beta, hwr);
            p.hy.kh = alpha * p.hy.kv;
            TPZGeoMesh *gmesh = HydraulicMesh(p, href);
            SeepageResult r = DrawdownSeepage(gmesh, p.hyd, p.hy);
            delete gmesh;
            auto it = paper.find({int(std::lround(alpha)), int(std::lround(beta))});
            std::cout << std::setprecision(6) << alpha << "  " << beta << "  " << r.nel << "  " << r.neq << "  " << r.Jnorm << "  ";
            if (it != paper.end()) std::cout << it->second << "  " << std::setprecision(3) << 100. * (r.Jnorm / it->second - 1.) << "%";
            else std::cout << "-  -";
            std::cout << std::setprecision(3) << "  " << r.pmin << "  " << r.seconds << "\n";
        }
    return 0;
}

int CmdFS(Options &o) {
    Problem p = ReadProblem(o);
    const std::string water = o.Str("water", "seepage"), form = o.Str("form", "u");
    const int nref = int(o.Num("nref", 3));
    const bool srm = o.Num("srm", 0) != 0.;
    const std::string vtk = o.Str("vtk", "");
    const int href = o.IntList("href", {0})[0];
    DriverSettings ds;
    ds.maxNewton = int(o.Num("maxnewton", ds.maxNewton));
    ds.tolFS = o.Num("tolfs", ds.tolFS);
    ds.markFrac = o.Num("mark", ds.markFrac);
    const bool checkForms = o.Num("checkforms", 0) != 0.;
    o.CheckUnused();
    const auto t0 = std::chrono::steady_clock::now();
    ForceField f = NoSeepage();
    REAL gammaRef = p.soil.gamma;
    std::shared_ptr<const PoreField> field;
    if (water == "seepage") {
        TPZGeoMesh *hmesh = HydraulicMesh(p, href);
        SeepageResult r = DrawdownSeepage(hmesh, p.hyd, p.hy);
        delete hmesh;
        PrintGeometry("hydraulic domain", p.hyd);
        PrintHydraulics(p.hy, p.hmesh);
        std::cout << "seepage: " << r.neq << " equations, J/(kh H^2 gw^2) " << r.Jnorm << ", min p " << r.pmin << " kPa, "
                  << r.seconds << " s\n";
        field = r.field;
        if (form == "u") {
            f = r.field->AsForceField(r.field);
            gammaRef = p.soil.gamma - p.hy.gammaw; // b = lambda (gamma' g - grad u)
        } else if (form == "p" || form == "p+") {
            f = TotalPressureForce(r.field, p.hy.gammaw, form == "p+"); // b = lambda (gamma_sat g - grad p[+])
        } else {
            std::cerr << "form must be u, p or p+\n";
            return 1;
        }
    } else if (water != "dry") {
        std::cerr << "water must be seepage or dry\n";
        return 1;
    }
    MeshStats st;
    TPZGeoMesh *smesh = StabilityMesh(p, &st);
    PrintGeometry("stability domain", p.stab);
    std::cout << "stability mesh " << p.smesh << ": " << CountLeafTriangles(smesh) << " triangles; water " << water << ", form "
              << form << ", gamma_ref " << gammaRef << "; Newton cap " << ds.maxNewton << ", continuation tol " << ds.tolFS << "\n";
    const TMCVoigt model = ModelVoigt(p.soil);
    if (field) { // where p < 0 the forms u / p and p+ differ
        TPZCompMesh *cmesh = CreateCMesh(smesh, 2, model, p.soil);
        int64_t nneg, ntot, nout;
        REAL pmin;
        NegativePressurePoints(cmesh, *field, p.hy.gammaw, nneg, ntot, pmin, nout);
        std::cout << "total pore pressure at the " << ntot << " integration points of the initial stability mesh: " << nneg
                  << " with p < 0, min p " << pmin << " kPa\n";
        delete cmesh;
        if (nout > 0) { // the seepage force would silently be zero there
            std::cerr << nout << " integration points of the stability mesh lie outside the hydraulic mesh: enlarge the "
                      << "hydraulic box (hleft, hright, hdepth) or reduce the stability box (sa, sleft, sright, sdepth)\n";
            return 1;
        }
        if (checkForms) { // b = lambda (gamma' g - grad u) and lambda (gamma_sat g - grad p) give the same load vector
            const TPZFMatrix<STATE> Fu = LoadVector(smesh, model, p.soil, field->AsForceField(field), p.soil.gamma - p.hy.gammaw);
            const TPZFMatrix<STATE> Fp = LoadVector(smesh, model, p.soil, TotalPressureForce(field, p.hy.gammaw, false), p.soil.gamma);
            const TPZFMatrix<STATE> Fd = LoadVector(smesh, model, p.soil, NoSeepage(), p.soil.gamma - p.hy.gammaw);
            TPZFMatrix<STATE> d(Fu);
            d -= Fp;
            TPZFMatrix<STATE> s(Fu);
            s -= Fd;
            std::cout << "load vectors at lambda = 1: |F_u| " << Norm(Fu) << ", |F_u - F_p| / |F_u| " << Norm(d) / Norm(Fu)
                      << ", seepage part |F_u - F(gamma' g)| / |F_u| " << Norm(s) / Norm(Fu) << "\n";
        }
    }
    std::vector<FSCycle> cyc = GravityIncreaseFS(smesh, model, p.soil, f, gammaRef, nref, srm, vtk, ds);
    delete smesh;
    const double total = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    std::cout << "\ncycle  equations  lambda_GI (= Gamma_FEM)  H_crit (m)  gamma H_crit / c" << (srm ? "  FS_SRM" : "")
              << "  time (s)  plastic zone / H (mark; 1 %): x_min  x_max - x_T  y_min - y_T\n";
    const REAL H = p.stab.H;
    for (const FSCycle &c : cyc) {
        std::cout << "  " << c.cycle << "  " << c.neq << "  " << c.gi << "  " << c.gi * H << "  " << p.soil.gamma * c.gi * H / p.soil.c;
        if (srm) std::cout << "  " << c.srm;
        std::cout << "  " << c.seconds << "  " << c.zone[0] / H << "  " << (c.zone[2] - p.stab.XT()) / H << "  "
                  << (c.zone[1] + H) / H << ";  " << c.zone1[0] / H << "  " << (c.zone1[2] - p.stab.XT()) / H << "  "
                  << (c.zone1[1] + H) / H << "\n";
    }
    std::cout << "box / H: x_min " << -p.stab.left / H << ", x_max - x_T " << p.stab.right / H << ", y_min - y_T "
              << -p.stab.depth / H << "\n";
    std::cout << "total time " << total << " s\n";
    return 0;
}

// ------------------------------------------------------------------------------------------------------------------
// check: focused self-tests of the conventions and of the coupling (about 1 min, no nonlinear analysis)
// ------------------------------------------------------------------------------------------------------------------

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

struct Checker {
    int total = 0, failed = 0;
    /// value <= tol passes
    void Expect(const std::string &what, REAL value, REAL tol) {
        total++;
        const bool ok = value <= tol && std::isfinite(value);
        if (!ok) failed++;
        std::cout << (ok ? "  [ok]   " : "  [FAIL] ") << what << ": " << std::setprecision(3) << value << " (tol " << tol
                  << ")" << std::setprecision(6) << std::endl;
    }
};

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
    TPZGeoMesh *hmesh = HydraulicMesh(p, 0);
    SeepageResult r = DrawdownSeepage(hmesh, p.hyd, p.hy);
    delete hmesh;
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

int CmdCheck(Options &o) {
    Problem p0 = ReadProblem(o);
    o.CheckUnused();
    const auto t0 = std::chrono::steady_clock::now();
    Checker ck;
    std::cout << "manufactured solutions (alpha = 5, P2): exact in the discrete space\n";
    {
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
    std::cout << "\ncheck: " << ck.total - ck.failed << " of " << ck.total << " passed, "
              << std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count() << " s\n";
    return ck.failed == 0 ? 0 : 1;
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        std::cout << "usage: SlopeSeepageForces mesh|verify|seepage|fig5|probe|fs|check [key=value ...] (see main.cpp)\n";
        return 1;
    }
    Options o(argc, argv, 2);
    const std::string cmd = argv[1];
    if (cmd == "mesh") return CmdMesh(o);
    if (cmd == "verify") return CmdVerify(o);
    if (cmd == "seepage") return CmdSeepage(o);
    if (cmd == "fig5") return CmdFig5(o);
    if (cmd == "probe") return CmdProbe(o);
    if (cmd == "fs") return CmdFS(o);
    if (cmd == "check") return CmdCheck(o);
    std::cerr << "unknown command " << cmd << "\n";
    return 1;
}
