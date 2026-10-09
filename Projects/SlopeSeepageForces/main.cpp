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
//   mesh     meshes and their quality (vtk=1 writes them)
//   verify   manufactured solutions of the anisotropic seepage solver (exact for P2)
//   seepage  drawdown seepage: J, J / (k_h H^2 gamma_w^2), min p, for href=<list> uniform refinements
//   fs       gravity-increase factor of the stability domain with the seepage forces (or dry)
// Options (default):
//   slope/soil: H=5 beta=45 hw=1 (h_w / H) gamma=20 gammaw=9.81 c=10 phi=30 E=20000 nu=0.3
//   seepage:    alpha=1 (k_h / k_v) kv=1 horder=2
//   hydraulic mesh: hmesh=gen|trig hleft=50 hright=10 hdepth=30 (units of H, from O, T, T) hh0=0.025 hhs=0.0625
//               hgrade=0.15 hhmax=2 (sizes in units of H) href=0 (list for seepage, e.g. href=0,1,2)
//   stability mesh: smesh=gen|trig sa=2 (extents sa H + H / tan(beta)) sleft= sright= sdepth= (override, units of
//               H) sh0=0.25 shs=0.25 sgrade=0.25 shmax=1 sref=0
//   fs:         water=seepage|dry form=u|p nref=3 srm=0 vtk=<prefix>
//   output:     csv=<file> (seepage: u and f on a grid around the slope, paper coordinates)
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
#include <set>
#include <sstream>
#include <string>
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
    std::vector<int> IntList(const std::string &k, const std::vector<int> &def) {
        fUsed.insert(k);
        if (!Has(k)) return def;
        std::vector<int> v;
        std::stringstream ss(fKV.at(k));
        std::string item;
        while (std::getline(ss, item, ',')) v.push_back(atoi(item.c_str()));
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
};

Problem ReadProblem(Options &o) {
    Problem p;
    SlopeGeometry g;
    g.H = o.Num("H", 5.);
    g.beta = o.Num("beta", 45.);
    g.hw = o.Num("hw", 1.) * g.H;
    p.hmesh = o.Str("hmesh", "gen");
    p.smesh = o.Str("smesh", "gen");
    if (p.hmesh == "trig" || p.smesh == "trig") {
        if (std::fabs(g.H - 10.) > 1.e-12 || std::fabs(g.beta - 45.) > 1.e-12) {
            std::cerr << "mesh=trig is the slope of SlopeMohrCoulomb: use H=10 beta=45\n";
            exit(1);
        }
    }
    p.hyd = g;
    p.hyd.left = o.Num("hleft", 50.) * g.H;
    p.hyd.right = o.Num("hright", 10.) * g.H;
    p.hyd.depth = o.Num("hdepth", 30.) * g.H;
    if (p.hmesh == "trig") p.hyd = SlopeMohrCoulombGeometry(g.hw);
    p.stab = g;
    StabilityExtents(p.stab, o.Num("sa", 2.));
    if (o.Has("sleft")) p.stab.left = o.Num("sleft", 0.) * g.H;
    if (o.Has("sright")) p.stab.right = o.Num("sright", 0.) * g.H;
    if (o.Has("sdepth")) p.stab.depth = o.Num("sdepth", 0.) * g.H;
    if (p.smesh == "trig") p.stab = SlopeMohrCoulombGeometry(g.hw);
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

int64_t CountLeafTriangles(TPZGeoMesh *gmesh) {
    int64_t n = 0;
    for (int64_t i = 0; i < gmesh->NElements(); i++) {
        TPZGeoEl *gel = gmesh->Element(i);
        if (gel && !gel->HasSubElement() && gel->Dimension() == 2) n++;
    }
    return n;
}

int CmdMesh(Options &o) {
    Problem p = ReadProblem(o);
    const bool vtk = o.Num("vtk", 0) != 0.;
    const int href = o.IntList("href", {0})[0];
    o.CheckUnused();
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
        if (vtk) {
            std::ofstream f(k == 0 ? "mesh_hydraulic.vtk" : "mesh_stability.vtk");
            TPZVTKGeoMesh::PrintGMeshVTK(gmesh, f, true);
        }
        delete gmesh;
    }
    return 0;
}

/// Manufactured solutions: u_ex linear and quadratic with div(K grad u_ex) = 0 (kv x^2 - kh y^2 and xy are
/// K-harmonic), Dirichlet on the ground surface, exact flux (K grad u_ex) . n on the base and the sides
int CmdVerify(Options &o) {
    if (!o.Has("alpha")) std::cout << "(verify: alpha = 5 by default)\n";
    const REAL alpha = o.Num("alpha", 5.);
    Problem p = ReadProblem(o);
    p.hy.kh = alpha * p.hy.kv;
    const int href = o.IntList("href", {0})[0];
    o.CheckUnused();
    const REAL kh = p.hy.kh, kv = p.hy.kv, H = p.hyd.H;
    struct Exact {
        const char *name;
        REAL c[6]; ///< u = c0 + c1 x + c2 y + c3 (kv x^2 - kh y^2) + c4 x y (c5 unused)
    };
    const Exact cases[2] = {{"linear", {3., 2., -1.5, 0., 0., 0.}},
                            {"quadratic K-harmonic", {3., 2., -1.5, 0.4 / H, 0.3 / H, 0.}}};
    PrintGeometry("hydraulic domain", p.hyd);
    std::cout << "K = diag(" << kh << ", " << kv << "), order " << p.hy.order << "\n";
    for (const Exact &e : cases) {
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
        std::cout << std::scientific << std::setprecision(3) << e.name << ": " << neq << " equations, " << pf.Triangles().size()
                  << " triangles, max|u - u_ex| / max|u_ex| = " << eu / umax << ", max|grad e| / max|grad u_ex| = "
                  << eg / gmax << ", |J - J_ex| / J_ex = " << std::fabs(J - Jex) / Jex << std::defaultfloat
                  << std::setprecision(6) << "  (" << dt << " s)\n";
        delete cmesh;
        delete gmesh;
    }
    return 0;
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
    const std::string csv = o.Str("csv", "");
    o.CheckUnused();
    PrintGeometry("hydraulic domain", p.hyd);
    std::cout << "K = diag(" << p.hy.kh << ", " << p.hy.kv << "), alpha " << p.hy.kh / p.hy.kv << ", gamma_w " << p.hy.gammaw
              << ", order " << p.hy.order << ", mesh " << p.hmesh << "\n";
    std::vector<REAL> J;
    std::cout << "href  triangles  equations  J  J/(kh H^2 gw^2)  min p (kPa) at (x, y)  time (s)\n";
    for (int ref : refs) {
        TPZGeoMesh *gmesh = HydraulicMesh(p, ref);
        SeepageResult r = DrawdownSeepage(gmesh, p.hyd, p.hy);
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

int CmdFS(Options &o) {
    Problem p = ReadProblem(o);
    const std::string water = o.Str("water", "seepage"), form = o.Str("form", "u");
    const int nref = int(o.Num("nref", 3));
    const bool srm = o.Num("srm", 0) != 0.;
    const std::string vtk = o.Str("vtk", "");
    const int href = o.IntList("href", {0})[0];
    o.CheckUnused();
    const auto t0 = std::chrono::steady_clock::now();
    ForceField f = NoSeepage();
    REAL gammaRef = p.soil.gamma;
    if (water == "seepage") {
        TPZGeoMesh *hmesh = HydraulicMesh(p, href);
        SeepageResult r = DrawdownSeepage(hmesh, p.hyd, p.hy);
        delete hmesh;
        PrintGeometry("hydraulic domain", p.hyd);
        std::cout << "seepage: alpha " << p.hy.kh / p.hy.kv << ", " << r.neq << " equations, J/(kh H^2 gw^2) " << r.Jnorm
                  << ", min p " << r.pmin << " kPa, " << r.seconds << " s\n";
        if (form == "u") {
            f = r.field->AsForceField(r.field);
            gammaRef = p.soil.gamma - p.hy.gammaw; // b = lambda (gamma' g - grad u)
        } else if (form == "p") {
            f = TotalPressureForce(r.field, p.hy.gammaw); // b = lambda (gamma_sat g - grad p+)
        } else {
            std::cerr << "form must be u or p\n";
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
              << form << ", gamma_ref " << gammaRef << "\n";
    const TMCVoigt model = ModelVoigt(p.soil);
    std::vector<FSCycle> cyc = GravityIncreaseFS(smesh, model, p.soil, f, gammaRef, nref, srm, vtk);
    delete smesh;
    const double total = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    std::cout << "\ncycle  equations  lambda_GI (= Gamma_FEM)  H_crit (m)  gamma H_crit / c" << (srm ? "  FS_SRM" : "")
              << "  time (s)  plastic zone / H: x_min  x_max - x_T  y_min\n";
    const REAL H = p.stab.H;
    for (const FSCycle &c : cyc) {
        std::cout << "  " << c.cycle << "  " << c.neq << "  " << c.gi << "  " << c.gi * H << "  " << p.soil.gamma * c.gi * H / p.soil.c;
        if (srm) std::cout << "  " << c.srm;
        std::cout << "  " << c.seconds << "  " << c.zone[0] / H << "  " << (c.zone[2] - p.stab.XT()) / H << "  " << c.zone[1] / H
                  << "\n";
    }
    std::cout << "box / H: x_min " << -p.stab.left / H << ", x_max - x_T " << p.stab.right / H << ", y_min "
              << (-p.stab.H - p.stab.depth) / H << "\n";
    std::cout << "total time " << total << " s\n";
    return 0;
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        std::cout << "usage: SlopeSeepageForces mesh|verify|seepage|fs [key=value ...] (see main.cpp)\n";
        return 1;
    }
    Options o(argc, argv, 2);
    const std::string cmd = argv[1];
    if (cmd == "mesh") return CmdMesh(o);
    if (cmd == "verify") return CmdVerify(o);
    if (cmd == "seepage") return CmdSeepage(o);
    if (cmd == "fs") return CmdFS(o);
    std::cerr << "unknown command " << cmd << "\n";
    return 1;
}
