// Slope stability under the seepage forces of a rapid drawdown (Ceron, Cecilio, Linn & Maghous, IJNAMG 2025,
// doi:10.1002/nag.3993, Section 4.3), built on Projects/SlopeDrawdown and Projects/SlopeMohrCoulomb:
//  1. parametric slope meshes (SlopeGeometry.h, DelaunayMesher.h): large hydraulic domain and stability domain;
//  2. steady anisotropic seepage of the excess pore pressure u after the drawdown (SeepageFE.h,
//     AnisotropicDarcy.h): J(u_FE) = 1/2 int grad u . K grad u;
//  3. seepage force field f = -grad u_FE at any point (SeepageForceField.h);
//  4. FEM stability: gravity increase of SlopeAnalysis.h with b = lambda (gamma' g + f) (FEMStability.h);
//     Gamma_FEM = lambda_crit, H_crit = lambda_crit H.
//  5. kinematic limit analysis (LimitAnalysis.h): upper bound Gamma = min P_mr / (P_gamma + P_u) over rotational
//     log-spiral mechanisms with the seepage forces of any slope::ForceField (P_u by the boundary formula for the FE
//     field f = -grad u_h); H_crit = Gamma H.
//  6. semi-analytical optimal field K^-1 v'_opt of the paper (AnalyticalSeepage.h, port of
//     scripts/analytical_seepage.py: Eqs. 29-40 with the corrected Eq. 31 and the degenerate m -> 0 optimum of steep
//     slopes), a slope::ForceField like the FE field.
//
// Usage: SlopeSeepageForces <command> [key=value ...]
//   mesh     meshes, quality and consistency (vtk=1 writes them; sweep=1 checks beta = 15..90, h_w / H = 0..1)
//   verify   manufactured solutions of the anisotropic seepage solver (exact for P2) and point location
//   seepage  drawdown seepage: J, J / (k_h H^2 gamma_w^2), min p, for href=<list> uniform refinements
//   fig5     J / (k_h H^2 gamma_w^2) for alphas=<list> betas=<list> against the dashed curves of the paper's Fig. 5,
//            and -J*(v'_opt) / (k_h H^2 gamma_w^2) of the analytical field against the solid curves
//   analytical  analytical field K^-1 v'_opt (AnalyticalSeepage.h): m, C, D, F, J* (Eq. 40) and the force at sample
//            points; ref=<file> compares with the Python reference of scripts/analytical_seepage_reference.py
//   probe    u and f = -grad u at the points of pts=<file> (lines x_paper,y_paper), paper coordinates
//   fs       gravity-increase factor of the stability domain with the seepage forces (FE or analytical field, or dry)
//   check    self-tests (about 6 s, exit code 1 on failure): manufactured solutions and orientation of K; Eq. 21
//            and the far-side data on the solved field, boundary ids, p >= 0, f = -grad u, point location on
//            vertices / edges / just outside, continuity of the P2 field, evaluator from 4 threads; load vector of
//            the stability problem: resultant = -int u n ds - gamma' A e_y (divergence theorem), lambda scales
//            gamma' and f, forms u and p identical, threaded assembly = serial; limit analysis (LimitAnalysis.h): closed
//            forms vs quadrature and polygon, domain vs boundary P_u, limit_analysis.py values at fixed mechanisms, dry
//            stability numbers, Fig. 8 h_w = 0 ends, scale invariance, uniform field, threads, FE field vs Python;
//            analytical field (AnalyticalSeepage.h): Python reference (data/analytical_seepage_reference.csv), Fig. 5
//            solid curves, degenerate optimum vs A > P, J* by quadrature vs Eq. 40, div v = 0, f = 0 outside, threads;
//            figure drivers: resumable csv (keys, cut rows), hydraulic box fixed in metres, Fig. 8 soil sets, tiny
//            fig8 / fig9 runs repeated (the second must skip every case); fembatch: extrapolations, stability box
//            inside the hydraulic box, similarity (load vectors at H and 4 H, lambda H of tiny runs), tiny runs repeated
//   la       kinematic limit analysis (LimitAnalysis.h, port of scripts/limit_analysis.py): Gamma = min P_mr / (P_gamma
//            + P_u) over log-spiral mechanisms I (B on the face) and II (B on the toe ground), H_crit = Gamma H;
//            prints Gamma, the mechanism (theta1, theta2, eta or d/H; A, B, C), P_mr, P_gamma, P_u
//   labatch  la for every line of cases=<file> (key=value options per line), rows appended to out=<csv>; cases
//            already in out are skipped (resumable)
//   fig8     paper Fig. 8: H_crit = Gamma(H) H versus h_w / H (alpha = 1, H = 1 m), London (beta = 30, 60) and Israeli
//            (35, 60) clay panels, curves vopt (K^-1 v'_opt) and FE (-grad u'_FE, box 50 / 10 / 30 m), limit analysis
//   fig9     paper Fig. 9: Gamma versus beta (H = 5 m, h_w = H) for alphas=1,5,10, curves vopt and FE (box in metres)
//            fig5 out=, fig8 and fig9 append one row per case to a resumable csv (default results/cpp/fig*.csv;
//            cases already there with the same settings are skipped) and log the progress to the .log next to it
//   fembatch FEM gravity-increase factor (FEMStability.h) of cases of Fig. 9 (fig=9: alphas=1,5,10 x betas=30,45,60,
//            75,90 x curves FE, vopt) or Fig. 8 (fig=8: London30/60, Israeli35/60 x hws=0,0.2,0.5,1, h_w = 0 once),
//            an independent check of the limit-analysis curves: each row holds lambda of every refinement cycle, the
//            extrapolations to h -> 0 and the limit analysis of the same case; resumable csv results/cpp/fem_fig9.csv
//            / fem_fig8.csv (scripts/run_fem_batch.sh)
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
//               water=fe (= seepage) | none (h_w = 0: f = 0 with gamma' = gamma - gamma_w) hboxm=<left>,<right>,<depth>
//               (hydraulic box in metres, as la; overrides hleft, hright, hdepth)
//   fig5:     alphas=1,2,4,10 betas=15,30,45,60,75,90 paper=<csv> (default data/fig5_vector_fill_polygons.csv)
//               fe=1 (0: analytical curve only) lm=10 (L_m / H of the analytical field)
//   analytical: lm=10 (L_m / H, R_e) pts=<file> (lines x_paper,y_paper: zone and f, paper coordinates) csv=<file>
//               (grid) ref=<file> (Python reference, e.g. data/analytical_seepage_reference.csv) bench=<n> (timing)
//   fs water=analytical: K^-1 v'_opt with gamma' = gamma - gamma_w (lm=10)
//   fig5 out=<csv>: production grid alphas=1,2,4,10 betas=15,22.5,..,90 href=1 (resumable csv, see fig8 / fig9 block)
//   fig8:       soil=swapped|table1 (Table 1 (c, phi) pairs exchanged between the panels | as printed) gammaw=9.8 H=1
//               panels=London30,London60,Israeli35,Israeli60 hws=0,0.05,0.1,0.2,..,1 curves=vopt,FE out=<csv>
//               (default results/cpp/fig8.csv, fig8_table1.csv for soil=table1)
//   fig9:       H=5 c=10 phi=30 gamma=20 gammaw=9.81 hw=1 alphas=1,5,10 betas=15,20,..,90 (with 37.5, 52.5, 67.5, 82.5)
//               curves=vopt,FE out=<csv> (default results/cpp/fig9.csv)
//   fig8, fig9: limit analysis as la (seeds=0,1,2 np=40 niter=150 pools=25 qsearch=coarse qfinal=fine polish=1
//               dmax=10 mech=I,II threads=2 lm=10) with hboxm=50,10,30 (hydraulic box in metres) href=1 hbc=zero_lb
//   output:     seepage: csv=<file> (u and f on a grid around the slope, paper coordinates), vtk=<file> (u, -grad u
//               and the Darcy velocity on the hydraulic mesh, NeoPZ coordinates), both for the last href
//   la:         water=fe|analytical|none|dry (FE field -grad u_h | K^-1 v'_opt (lm=10) | f = 0 with gamma' = gamma -
//               gamma_w (h_w = 0) | f = 0 and gamma_w = 0) hboxm=<left>,<right>,<depth> (hydraulic box in metres, e.g. 50,10,30; overrides hleft..)
//               href=0 pu=auto|domain|boundary (P_u: auto = boundary formula for the FE field) mech=I,II seeds=0,1,2
//               np=40 niter=150 pools=25 (PSO; pools: further initial pools of 4 np admissible mechanisms when
//               fewer than np / 2 have P_ext > 0; 0 = limit_analysis.py) qsearch=coarse qfinal=fine (coarse|medium|fine|xfine|ref|dense) polish=1
//               dmax=10 (d / H of mechanism II) threads=4 (runs class x seed in parallel) verbose=0 out=<csv> (append
//               a row; A, B, C in paper coordinates) x=theta1,theta2,s (no optimisation: rates of work of this mechanism by every rule)
//   fembatch:   fig=9|8; FEM nref=3 mark=0.05 sa=2 (stability box sa H + H / tan(beta), capped at the hydraulic box)
//               sh0=0.25 shs=0.25 sgrade=0.25 shmax=1 maxnewton=100 tolfs=0.002; limit analysis and seepage field as
//               fig8 / fig9 (href=1 also for the FE field of the FEM); cases: fig=9 H=5 c=10 phi=30 gamma=20
//               gammaw=9.81 hw=1 alphas= betas= curves=; fig=8 soil=swapped gammaw=9.8 panels= hws= curves= (H_ref =
//               H_crit of the limit analysis, 3 digits, box 50 / 10 / 30 H_ref); out=<csv>
// mesh=trig: TriGMesh(1 + ref) of SlopeMohrCoulomb (H = 10 m, beta = 45 deg, 70 x 40 m), the mesh of SlopeDrawdown.
#include "AnalyticalSeepage.h"
#include "FEMStability.h"
#include "LimitAnalysis.h"
#include "SeepageFE.h"

#include "TPZVTKGeoMesh.h"

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <filesystem>
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
/// (solid = true: -J*(v'_opt) / (k_h H^2 gamma_w^2), the solid curves = lower edges of the bands)
std::map<std::pair<int, int>, REAL> ReadFig5(const std::string &file, bool solidCurve = false) {
    std::map<std::pair<int, int>, REAL> m;
    std::ifstream f(file);
    std::string line;
    while (std::getline(f, line)) {
        if (line.empty() || line[0] == '#') continue;
        REAL a, b, solid, dashed;
        if (sscanf(line.c_str(), "%lf,%lf,%lf,%lf", &a, &b, &solid, &dashed) == 4)
            m[{int(std::lround(a)), int(std::lround(b))}] = solidCurve ? solid : dashed;
    }
    return m;
}

/// the analytical field K^-1 v'_opt (AnalyticalSeepage.h) of the slope and water level of p.hyd (alpha = k_h / k_v),
/// L_m = lm H
std::shared_ptr<const AnalyticalSeepage> MakeAnalyticalField(const Problem &p, REAL lm) {
    return std::make_shared<const AnalyticalSeepage>(p.hyd.beta, p.hyd.H, p.hyd.hw, p.hy.kh / p.hy.kv, p.hy.kh, p.hy.gammaw, lm);
}

int CmdFig5CSV(Options &o); // fig5 out=<csv>: production run into a resumable csv (fig8 / fig9 block below)

/// J(u_FE) / (k_h H^2 gamma_w^2) versus beta and alpha (paper Fig. 5, dashed curves) on the hydraulic domain, and
/// -J*(v'_opt) / (k_h H^2 gamma_w^2) of the analytical field (AnalyticalSeepage.h, L_m = lm H; solid curves)
int CmdFig5(Options &o) {
    if (o.Has("out")) return CmdFig5CSV(o);
    Problem p = ReadProblem(o);
    const std::vector<double> alphas = o.NumList("alphas", {1., 2., 4., 10.});
    const std::vector<double> betas = o.NumList("betas", {15., 30., 45., 60., 75., 90.});
    const int href = o.IntList("href", {0})[0];
    const std::string src(__FILE__);
    const std::string file = o.Str("paper", src.substr(0, src.find_last_of('/') + 1) + "data/fig5_vector_fill_polygons.csv");
    const bool fe = o.Num("fe", 1.) != 0.;
    const REAL lm = o.Num("lm", 10.);
    const REAL hwr = p.hyd.hw / p.hyd.H;
    o.CheckUnused();
    const auto paper = ReadFig5(file), paperSolid = ReadFig5(file, true);
    if (paper.empty()) std::cout << "(paper data not found: " << file << ")\n";
    if (fe) {
        PrintGeometry("hydraulic domain", p.hyd);
        PrintHydraulics(p.hy, p.hmesh);
    }
    std::cout << "h_w / H " << hwr << ", href " << href << ", analytical field: L_m / H " << lm << "\nalpha  beta"
              << (fe ? "  triangles  equations  J/(kh H^2 gw^2)  paper (dashed)  difference  min p (kPa)  time (s)" : "")
              << "  -J*(v'_opt)/(kh H^2 gw^2)  paper (solid)  difference  m\n";
    for (double alpha : alphas)
        for (double beta : betas) {
            p.SetSlope(beta, hwr);
            p.hy.kh = alpha * p.hy.kv;
            SeepageResult r;
            if (fe) { // solve before printing the row (the solver writes to stdout)
                TPZGeoMesh *gmesh = HydraulicMesh(p, href);
                r = DrawdownSeepage(gmesh, p.hyd, p.hy);
                delete gmesh;
            }
            std::cout << std::setprecision(6) << alpha << "  " << beta;
            if (fe) {
                auto it = paper.find({int(std::lround(alpha)), int(std::lround(beta))});
                std::cout << "  " << r.nel << "  " << r.neq << "  " << r.Jnorm << "  ";
                if (it != paper.end()) std::cout << it->second << "  " << std::setprecision(3) << 100. * (r.Jnorm / it->second - 1.) << "%";
                else std::cout << "-  -";
                std::cout << std::setprecision(3) << "  " << r.pmin << "  " << r.seconds;
            }
            const auto an = MakeAnalyticalField(p, lm);
            const REAL Jn = -an->JstarNormalized();
            auto it = paperSolid.find({int(std::lround(alpha)), int(std::lround(beta))});
            std::cout << std::setprecision(6) << "  " << Jn << "  ";
            if (it != paperSolid.end()) std::cout << it->second << "  " << std::setprecision(3) << 100. * (Jn / it->second - 1.) << "%";
            else std::cout << "-  -";
            std::cout << std::setprecision(6) << "  " << an->M() << (an->Degenerate() ? " (degenerate)" : "") << "\n";
        }
    return 0;
}

int CmdFS(Options &o) {
    Problem p = ReadProblem(o);
    { // hboxm=<left>,<right>,<depth>: hydraulic box in metres (as la / fig8 / fig9; overrides hleft, hright, hdepth)
        const std::vector<double> boxm = o.NumList("hboxm", {});
        if (!boxm.empty() && boxm.size() != 3) {
            std::cerr << "fs: hboxm=<left>,<right>,<depth> (m)\n";
            return 1;
        }
        for (size_t i = 0; i < boxm.size(); i++) p.hext[i] = boxm[i] / p.hyd.H;
        if (!boxm.empty()) p.SetSlope(p.hyd.beta, p.hyd.hw / p.hyd.H);
    }
    std::string water = o.Str("water", "seepage");
    if (water == "fe") water = "seepage"; // the name used by la
    const std::string form = o.Str("form", "u");
    const int nref = int(o.Num("nref", 3));
    const bool srm = o.Num("srm", 0) != 0.;
    const std::string vtk = o.Str("vtk", "");
    const int href = o.IntList("href", {0})[0];
    DriverSettings ds;
    ds.maxNewton = int(o.Num("maxnewton", ds.maxNewton));
    ds.tolFS = o.Num("tolfs", ds.tolFS);
    ds.markFrac = o.Num("mark", ds.markFrac);
    const bool checkForms = o.Num("checkforms", 0) != 0.;
    const REAL lm = o.Num("lm", 10.); // water=analytical: L_m / H
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
    } else if (water == "analytical") { // K^-1 v'_opt (AnalyticalSeepage.h): b = lambda (gamma' g + f)
        if (form != "u") {
            std::cerr << "water=analytical: form must be u\n";
            return 1;
        }
        const auto an = MakeAnalyticalField(p, lm);
        std::cout << "analytical field: " << an->Summary() << "\n";
        f = AnalyticalForceField(an);
        gammaRef = p.soil.gamma - p.hy.gammaw;
    } else if (water == "none") { // h_w = 0 (submerged slope, as la water=none): f = 0 with the buoyant gamma'
        gammaRef = p.soil.gamma - p.hy.gammaw;
    } else if (water != "dry") {
        std::cerr << "water must be seepage (fe), analytical, none or dry\n";
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
    const FSExtrapolation ex = ExtrapolateCycles(cyc); // FEMStability.h (NaN: not enough cycles)
    std::cout << "lambda extrapolated to h -> 0 from the last cycles: order 1 " << ex.h1 << ", observed order " << ex.order
              << " -> " << ex.richardson << ", linear in 1 / sqrt(neq) " << ex.sqrtneq << "\n";
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

// ------------------------------------------------------------------------------------------------------------------
// la, labatch: kinematic limit analysis (LimitAnalysis.h, port of scripts/limit_analysis.py) and its self-tests
// ------------------------------------------------------------------------------------------------------------------

/// seepage input of the limit analysis: water = fe (FE field -grad u_h on the hydraulic domain of p, href uniform
/// refinements; P_u by the boundary formula unless pu=domain), analytical (K^-1 v'_opt of AnalyticalSeepage.h, L_m =
/// lm H; P_u by the domain quadrature with the rays split on its discontinuity circles r = R_w, R, R_e about O and on
/// circles R_w 0.3^k graded towards the singular point O, as limit_analysis.py), none (h_w = 0: f = 0 with the
/// buoyant gamma') or dry (f = 0 and gamma_w = 0)
struct LAField {
    la::Seepage seep;
    std::string info;
    REAL Jnorm = 0.;
    double seconds = 0.;
};

LAField MakeLAField(const Problem &p, const std::string &water, int href, REAL lm = 10.) {
    LAField lf;
    if (water == "fe") {
        TPZGeoMesh *gmesh = HydraulicMesh(p, href);
        SeepageResult r = DrawdownSeepage(gmesh, p.hyd, p.hy);
        delete gmesh;
        lf.seep = la::FESeepage(r.field, p.hyd.hw);
        lf.Jnorm = r.Jnorm, lf.seconds = r.seconds;
        std::ostringstream s;
        s << "FE field: hydraulic box x in [" << -p.hyd.left << ", " << p.hyd.XT() + p.hyd.right << "] m, y >= "
          << -p.hyd.H - p.hyd.depth << " m, far sides left " << FarSides::Name(p.hy.far.left) << " / base "
          << FarSides::Name(p.hy.far.bottom) << " / right " << FarSides::Name(p.hy.far.right) << ", order " << p.hy.order
          << ", href " << href << ": " << r.neq << " equations, J/(kh H^2 gw^2) " << std::setprecision(8) << r.Jnorm
          << ", " << std::setprecision(3) << r.seconds << " s";
        lf.info = s.str();
    } else if (water == "analytical") {
        const auto t0 = std::chrono::steady_clock::now();
        const auto an = MakeAnalyticalField(p, lm);
        lf.seep.force = AnalyticalForceField(an);
        for (REAL R : {an->Rw(), an->R(), an->Re()}) {
            bool dup = !(R > 0.);
            for (const la::Circle &c : lf.seep.circles) dup = dup || c.R == R;
            if (!dup) lf.seep.circles.push_back({0., 0., R, true});
        }
        if (an->Rw() > 0.) lf.seep.grade = {0., 0., an->Rw(), false}, lf.seep.gradeRatio = 0.3;
        lf.Jnorm = -an->JstarNormalized();
        lf.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        lf.info = "analytical field K^-1 v'_opt: " + an->Summary();
    } else if (water == "none" || water == "dry") {
        lf.info = water == "none" ? "no seepage (h_w = 0): f = 0, buoyant gamma' = gamma - gamma_w" : "dry: f = 0, gamma_w = 0";
    } else {
        std::cerr << "la: water must be fe, analytical, none or dry\n";
        exit(1);
    }
    return lf;
}

/// one limit analysis with the options of 'la' (see the header of this file); id: case label of the csv row
int RunLA(Options &o, const std::string &id, int defThreads, const std::string &csvOut = "") {
    Problem p = ReadProblem(o);
    const REAL beta = p.hyd.beta, H = p.hyd.H, hwr = p.hyd.hw / H;
    const std::vector<double> boxm = o.NumList("hboxm", {}); // hydraulic box in metres (left of O, right of T, below T)
    if (!boxm.empty()) {
        if (boxm.size() != 3) {
            std::cerr << "la: hboxm=<left>,<right>,<depth> (m)\n";
            return 1;
        }
        for (int i = 0; i < 3; i++) p.hext[i] = boxm[i] / H;
        p.SetSlope(beta, hwr);
    }
    const REAL phiDeg = o.Num("phi", 30.), c = p.soil.c, gamma = p.soil.gamma; // phi in degrees, as given
    const std::string water = o.Str("water", hwr > 0. ? "fe" : "none");
    const int href = o.IntList("href", {0})[0];
    const REAL lm = o.Num("lm", 10.); // water=analytical: L_m / H
    const std::string pu = o.Str("pu", "auto"), mech = o.Str("mech", "I,II"), out = csvOut.empty() ? o.Str("out", "") : csvOut;
    la::Settings st;
    st.kinds.clear();
    if (mech == "I" || mech == "I,II" || mech == "II,I") st.kinds.push_back(1);
    if (mech == "II" || mech == "I,II" || mech == "II,I") st.kinds.push_back(2);
    if (st.kinds.empty()) {
        std::cerr << "la: mech must be I, II or I,II\n";
        return 1;
    }
    st.seeds.clear();
    for (int s : o.IntList("seeds", {0, 1, 2})) st.seeds.push_back(uint64_t(s));
    st.nParticles = int(o.Num("np", 40));
    st.nIter = int(o.Num("niter", 150));
    st.quadSearch = o.Str("qsearch", "coarse");
    st.quadFinal = o.Str("qfinal", "fine");
    st.polish = o.Num("polish", 1) != 0.;
    st.extraPools = int(o.Num("pools", 25));
    st.nThreads = int(o.Num("threads", defThreads));
    const REAL dmax = o.Num("dmax", 10.);
    const std::vector<double> xgiven = o.NumList("x", {}); // theta1,theta2,s: evaluate this mechanism only
    const bool verbose = o.Num("verbose", 0) != 0.;
    o.CheckUnused();
    if (!la::Quadrature::Level(st.quadSearch).Valid() || !la::Quadrature::Level(st.quadFinal).Valid()) {
        std::cerr << "la: qsearch / qfinal must be coarse, medium, fine, xfine, ref or dense\n";
        return 1;
    }
    const REAL gammaw = water == "dry" ? 0. : p.hy.gammaw;
    LAField lf = MakeLAField(p, water, href, lm);
    la::Problem prob(beta, H, c, phiDeg, gamma, gammaw, lf.seep, pu, dmax);
    std::cout << std::setprecision(10) << "limit analysis: H " << H << " m, beta " << beta << " deg, h_w / H " << hwr
              << ", c " << c << " kPa, phi " << phiDeg << " deg, gamma " << gamma << ", gamma_w " << gammaw
              << " (gamma' " << prob.gammap << "), alpha " << p.hy.kh / p.hy.kv << "; water " << water << ", P_u "
              << la::PuMethodName(prob.pu) << "\n";
    if (!lf.info.empty()) std::cout << lf.info << "\n";
    if (!xgiven.empty()) { // one mechanism: rates of work by every available rule
        if (xgiven.size() != 3) {
            std::cerr << "la: x=theta1,theta2,s\n";
            return 1;
        }
        const la::Mechanism m = la::MakeMechanism(prob.geo, xgiven[0], xgiven[1], xgiven[2]);
        std::cout << std::setprecision(12) << "mechanism " << (m.II ? "II" : "I") << " x = (" << xgiven[0] << ", " << xgiven[1]
                  << ", " << xgiven[2] << "): admissible " << m.ok << ", r0 " << m.r0 << ", C (paper) = (" << m.Cx << ", "
                  << -m.Cy << "), A_x " << m.Ax << ", B (paper) = (" << m.Bx << ", " << m.By << ")\n";
        if (!m.ok) return 1;
        const la::Quadrature qf = la::Quadrature::Level(st.quadFinal);
        const la::Powers P = prob.Evaluate(m, qf);
        std::cout << std::setprecision(15) << "P_mr " << P.Pmr << "\nP_gamma closed form " << P.Pgamma << ", quadrature ("
                  << qf.name << ") " << prob.PgammaQuadrature(m, qf) << "\n";
        if (prob.pu != la::PuMethod::ENone) {
            std::cout << "P_u (" << la::PuMethodName(prob.pu) << ", " << qf.name << ") " << P.Pu << "\n";
            if (lf.seep.u) std::cout << "P_u boundary formula " << prob.Pu(m, qf, la::PuMethod::EBoundary) << "\n";
            for (const char *lev : {"fine", "xfine", "ref"})
                std::cout << "P_u domain quadrature " << lev << " " << prob.Pu(m, la::Quadrature::Level(lev), la::PuMethod::EDomain) << "\n";
        }
        std::cout << "Gamma = P_mr / (P_gamma + P_u) = " << P.Pmr / (P.Pgamma + P.Pu) << "\n";
        return 0;
    }
    const la::Result r = la::StabilityFactor(prob, st);
    if (verbose || !r.found) {
        std::cout << "runs (class, seed: PSO -> Nelder-Mead, evaluations, s):\n";
        for (const la::Run &run : r.runs)
            std::cout << "  " << (run.kind == 1 ? "I " : "II") << " " << run.seed << ": " << std::setprecision(10)
                      << (run.found ? run.GammaPSO : INFINITY) << " -> " << (run.found ? run.Gamma : INFINITY) << "  x = ("
                      << run.x[0] << ", " << run.x[1] << ", " << run.x[2] << ")  " << run.nEval << "  "
                      << std::setprecision(3) << run.seconds << "\n";
    }
    const la::Mechanism &m = r.mech;
    if (!r.found) std::cout << "no admissible mechanism with P_gamma + P_u > 0: Gamma = inf\n";
    else {
        std::cout << std::setprecision(10) << "Gamma = " << r.Gamma << "  (H_crit = Gamma H = " << r.Hcrit
                  << " m, gamma H_crit / c = " << gamma * r.Hcrit / c << "), mechanism " << (r.kind == 1 ? "I" : "II")
                  << "\n  theta1 " << r.x[0] << ", theta2 " << r.x[1] << (r.kind == 1 ? ", eta " : ", d / H ")
                  << (r.kind == 1 ? r.x[2] : r.x[2] - 1.) << ", theta_O " << m.thO << "\n  A = (" << r.A[0] << ", "
                  << r.A[1] << "), B = (" << r.B[0] << ", " << r.B[1] << "), C = (" << r.C[0] << ", " << r.C[1]
                  << ") (paper coordinates, y down; NeoPZ: y = -y_paper)\n  r0 " << m.r0 << ", r_h " << m.rh << ", L = OA "
                  << m.L << "\n  P_mr " << r.P.Pmr << ", P_gamma " << r.P.Pgamma << ", P_u " << r.P.Pu
                  << " (per unit omega, " << st.quadFinal << ")\n  best by class: I " << r.bestByClass[1] << ", II "
                  << r.bestByClass[2] << "; seed spread " << std::setprecision(2) << r.seedSpread;
        for (const std::string &b : r.atBound) std::cout << "; at search bound " << b;
        if (r.dmaxReached) std::cout << "; d = d_max (deeper mechanisms may lower Gamma)";
        std::cout << "\n";
    }
    std::cout << std::setprecision(4) << "evaluations " << r.nEval << ", " << r.seconds << " s (field " << lf.seconds << " s)\n";
    if (!out.empty()) { // csv row (header if the file is new)
        std::ifstream test(out);
        const bool isNew = !test.good() || test.peek() == std::ifstream::traits_type::eof();
        test.close();
        std::ofstream f(out, std::ios::app);
        if (isNew)
            f << "id,H,beta_deg,hw_over_H,alpha,c,phi_deg,gamma,gamma_w,water,pu,Gamma,Hcrit,mechanism,theta1,theta2,s,"
                 "Ax,Bx,By,Cx,Cy,r0,rh,L,P_mr,P_gamma,P_u,seed_spread,n_eval,Jnorm,t_field,t_la\n";
        f << std::setprecision(12) << '"' << id << "\"," << H << "," << beta << "," << hwr << "," << p.hy.kh / p.hy.kv << ","
          << c << "," << phiDeg << "," << gamma << "," << gammaw << "," << water << "," << la::PuMethodName(prob.pu) << ",";
        if (r.found)
            f << r.Gamma << "," << r.Hcrit << "," << (r.kind == 1 ? "I" : "II") << "," << r.x[0] << "," << r.x[1] << "," << r.x[2]
              << "," << r.A[0] << "," << r.B[0] << "," << r.B[1] << "," << r.C[0] << "," << r.C[1] << "," << m.r0 << "," << m.rh
              << "," << m.L << "," << r.P.Pmr << "," << r.P.Pgamma << "," << r.P.Pu << "," << r.seedSpread;
        else f << "inf,inf,,,,,,,,,,,,,,,,";
        f << "," << r.nEval << "," << lf.Jnorm << "," << lf.seconds << "," << r.seconds << "\n";
    }
    return 0;
}

int CmdLA(Options &o, int argc, char *argv[]) {
    std::string id;
    for (int i = 2; i < argc; i++) id += (i > 2 ? " " : "") + std::string(argv[i]);
    return RunLA(o, id, int(std::min(4u, std::max(1u, std::thread::hardware_concurrency()))));
}

/// labatch cases=<file> out=<csv> [threads=]: one 'la' case per line of the file (key=value options; '#' comments);
/// cases whose line is already in the first column of out are skipped (resumable: rows are appended one by one)
int CmdLABatch(Options &o) {
    const std::string cases = o.Str("cases", ""), out = o.Str("out", "");
    const int threads = int(o.Num("threads", std::min(4u, std::max(1u, std::thread::hardware_concurrency()))));
    o.CheckUnused();
    std::ifstream in(cases);
    if (!in || out.empty()) {
        std::cerr << "labatch: cases=<file> out=<csv>\n";
        return 1;
    }
    std::set<std::string> done;
    {
        std::ifstream prev(out);
        std::string line;
        while (std::getline(prev, line))
            if (line.size() > 1 && line[0] == '"') done.insert(line.substr(1, line.find('"', 1) - 1));
    }
    std::string line;
    int n = 0, skipped = 0, failed = 0;
    while (std::getline(in, line)) {
        const size_t a = line.find_first_not_of(" \t"), b = line.find_last_not_of(" \t\r");
        if (a == std::string::npos || line[a] == '#') continue;
        const std::string id = line.substr(a, b - a + 1);
        if (done.count(id)) {
            skipped++;
            continue;
        }
        std::vector<std::string> tok;
        std::istringstream ss(id);
        for (std::string t; ss >> t;) tok.push_back(t);
        std::vector<char *> argv = {nullptr};
        for (std::string &t : tok) argv.push_back(&t[0]);
        Options oc(int(argv.size()), argv.data(), 1);
        std::cout << "\n=== case: " << id << std::endl;
        if (RunLA(oc, id, threads, out) != 0) failed++;
        n++;
    }
    std::cout << "\nlabatch: " << n << " cases run, " << skipped << " already in " << out << ", " << failed << " failed\n";
    return failed == 0 ? 0 : 1;
}

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
// analytical seepage field K^-1 v'_opt (AnalyticalSeepage.h): comparison with the Python reference
// (scripts/analytical_seepage_reference.py), J* by quadrature, self-tests and the command 'analytical'
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
    d.build = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
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
    d.eval = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count() / std::max<size_t>(c.pts.size(), 1);
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

/// Gauss-Legendre rule on [-1, 1] (Newton on the Legendre polynomial)
void GaussLegendreRule(int n, std::vector<REAL> &x, std::vector<REAL> &w) {
    x.resize(n), w.resize(n);
    for (int i = 0; i < n; i++) {
        REAL z = std::cos(M_PI * (i + 0.75) / (n + 0.5)), dp = 1.;
        for (int it = 0; it < 100; it++) {
            REAL p0 = 1., p1 = z;
            for (int k = 2; k <= n; k++) {
                const REAL p2 = ((2. * k - 1.) * z * p1 - (k - 1.) * p0) / k;
                p0 = p1, p1 = p2;
            }
            dp = n * (z * p1 - p0) / (z * z - 1.);
            const REAL dz = p1 / dp;
            z -= dz;
            if (std::fabs(dz) < 1.e-16) break;
        }
        x[i] = z, w[i] = 2. / ((1. - z * z) * dp * dp);
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
    GaussLegendreRule(64, xg, wg);
    auto ring = [&](REAL r, REAL thmax) { // int_0^thmax 1/2 f.K.f r dtheta
        REAL s = 0.;
        for (size_t i = 0; i < xg.size(); i++) {
            const REAL th = 0.5 * thmax * (xg[i] + 1.);
            TPZManVector<REAL, 3> X = {-r * std::cos(th), -r * std::sin(th), 0.};
            REAL F[2];
            f(X, F);
            s += wg[i] * 0.5 * (kh * F[0] * F[0] + kv * F[1] * F[1]);
        }
        return 0.5 * thmax * s * r;
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

/// self-tests of the analytical field (part of 'check')
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
    std::cout << "    (analytical field tests: " << std::chrono::duration<double>(std::chrono::steady_clock::now() - tStart).count()
              << " s)\n";
}

/// analytical K^-1 v'_opt: scalars of the optimum, J*, the force at sample points; pts=<file> (x_paper,y_paper lines),
/// csv=<file> (grid), ref=<file> (comparison with scripts/analytical_seepage_reference.py), bench=<n> (timing)
int CmdAnalytical(Options &o) {
    Problem p = ReadProblem(o);
    const REAL lm = o.Num("lm", 10.);
    const std::string pts = o.Str("pts", ""), csv = o.Str("csv", ""), ref = o.Str("ref", "");
    const int64_t bench = int64_t(o.Num("bench", 0.));
    o.CheckUnused();
    if (!ref.empty()) {
        auto t0 = std::chrono::steady_clock::now();
        const std::vector<AnalyticalRefCase> cases = ReadAnalyticalReference(ref);
        int64_t n = 0;
        for (auto &c : cases) n += int64_t(c.pts.size());
        std::cout << "reference " << ref << ": " << cases.size() << " cases, " << n << " points (read in "
                  << std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count() << " s)\n";
        if (cases.empty()) return 1;
        AnalyticalRefDiff w;
        AnalyticalReferenceReport(cases, w, true);
        std::cout << std::setprecision(3) << "all: " << w.nforce << " points with f != 0 (max |f - f_py| / |f_py| " << w.frel << "), "
                  << w.nzero << " with f = 0 (max |f| / gw " << w.fzero << "), " << w.nexcl << " excluded; max |m - m_py| " << w.dm
                  << ", max relative difference of F " << w.dF << " and of J* " << w.dJ << "; degenerate flags "
                  << (w.sameFlag ? "equal" : "DIFFERENT") << "\n";
        return w.frel <= 1.e-6 && w.fzero == 0. && w.sameFlag ? 0 : 1;
    }
    const SlopeGeometry &g = p.hyd;
    auto t0 = std::chrono::steady_clock::now();
    auto an = MakeAnalyticalField(p, lm);
    const double tb = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    const AnalyticalSeepage &a = *an;
    std::cout << std::setprecision(10) << "slope: H " << g.H << " m, beta " << g.beta << " deg, h_w " << g.hw << " m (h_w/H "
              << g.hw / g.H << "), alpha = kh/kv " << a.Alpha() << ", kh " << a.Kh() << ", gamma_w " << a.Gammaw() << ", L_m "
              << a.Lm() << " m\n";
    std::cout << "R_w " << a.Rw() << ", R " << a.R() << ", R_e " << a.Re() << " (m), Theta = pi - beta " << a.Theta()
              << ", A (Eq. 34) " << a.A() << ", int_R^Re 4 dr / (B r) " << a.I3() << "\n";
    std::cout << "h2 problem: m = sqrt(C/D) " << a.M() << (a.Degenerate() && a.Hw() > 0. ? " (degenerate optimum m -> 0)" : "")
              << ", C " << a.C() << ", D " << a.D() << ", h2(Theta) " << a.H2e() << ", C/h2e^2 " << a.C() / (a.H2e() * a.H2e())
              << ", D/h2e^2 " << a.D() / (a.H2e() * a.H2e()) << ", Phi_min " << a.Phi() << ", F " << a.F() << ", 1/A " << 1. / a.A()
              << ", a' " << a.APrime() << "\n";
    REAL J[3];
    a.JstarParts(J);
    std::cout << "J*(v'_opt) (Eq. 40) " << a.Jstar() << " = zone 1 " << J[0] << " + zone 2 " << J[1] << " + zone 3 " << J[2]
              << "; J*/(kh H^2 gw^2) " << a.JstarNormalized() << "\n";
    std::cout << std::setprecision(4) << "construction " << 1.e3 * tb << " ms (" << a.ShootingCount() << " integrations of h2)\n";
    // force at sample points (those of scripts/analytical_seepage.py, check 6), paper coordinates in units of H
    const REAL H = g.H, xT = a.XToe() / H, tb_ = g.beta < 90. ? std::tan(g.beta * M_PI / 180.) : 0.;
    const REAL sp[11][2] = {{-1., 0.25}, {-0.5, 0.5}, {-0.2, 0.2}, {-0.05, 0.6}, {0.1, g.beta < 90. ? 0.12 + 0.1 * tb_ : 0.5},
                            {0.5 * xT, 0.55}, {0.9 * xT, 1.}, {xT + 0.2, 1.1}, {xT + 1., 1.05}, {0., 1.5}, {-2., 2.}};
    std::cout << "force f = K^-1 v'_opt, paper coordinates (y down; angle from +x towards +y = down):\n"
              << "   x/H    y/H  zone  |f|/gw  angle (deg)  fx/gw  fy/gw\n";
    for (auto &q : sp) {
        const REAL x = q[0] * H, yp = q[1] * H;
        if (!a.InSoilPaper(x, yp)) continue;
        REAL f[2];
        a.ForcePaper(x, yp, f);
        std::cout << std::fixed << std::setprecision(2) << "  " << std::setw(5) << q[0] << "  " << std::setw(5) << q[1] << "  "
                  << a.Zone(x, yp) << std::setprecision(4) << "  " << std::hypot(f[0], f[1]) / a.Gammaw() << "  " << std::setprecision(1)
                  << std::setw(6) << std::atan2(f[1], f[0]) * 180. / M_PI << std::setprecision(4) << "  " << f[0] / a.Gammaw()
                  << "  " << f[1] / a.Gammaw() << std::defaultfloat << "\n";
    }
    if (!pts.empty()) { // x_paper,y_paper -> x, y, zone, fx, fy (paper coordinates)
        std::ifstream in(pts);
        std::string line;
        std::cout << "# x_paper,y_paper,zone,fx_paper,fy_paper\n" << std::setprecision(17);
        while (std::getline(in, line)) {
            REAL x, yp;
            if (line.empty() || line[0] == '#' || sscanf(line.c_str(), "%lf,%lf", &x, &yp) != 2) continue;
            REAL f[2];
            a.ForcePaper(x, yp, f);
            std::cout << x << "," << yp << "," << a.Zone(x, yp) << "," << f[0] << "," << f[1] << "\n";
        }
    }
    if (!csv.empty()) { // grid around the slope, as WriteGridCSV
        std::ofstream out(csv);
        out << "# x_paper,y_paper,zone,fx_paper,fy_paper (f = K^-1 v'_opt, kN/m^3); H " << g.H << " beta " << g.beta << " hw "
            << g.hw << " alpha " << a.Alpha() << "\n"
            << std::setprecision(10);
        const REAL d = H / 20.;
        for (REAL x = -2. * H; x <= a.XToe() + 2. * H + 1.e-9; x += d)
            for (REAL yp = 0.; yp <= 2. * H + 1.e-9; yp += d) {
                if (!a.InSoilPaper(x, yp)) continue;
                REAL f[2];
                a.ForcePaper(x, yp, f);
                out << x << "," << yp << "," << a.Zone(x, yp) << "," << f[0] << "," << f[1] << "\n";
            }
        std::cout << "grid written to " << csv << "\n";
    }
    if (bench > 0) { // evaluation through the ForceField (NeoPZ coordinates) at random points of a box around the slope
        const ForceField f = AnalyticalForceField(an);
        std::mt19937_64 rng(3);
        std::uniform_real_distribution<REAL> ux(-3. * H, a.XToe() + 3. * H), uy(-4. * H, 0.);
        std::vector<REAL> xs(2 * bench);
        for (int64_t i = 0; i < bench; i++) xs[2 * i] = ux(rng), xs[2 * i + 1] = uy(rng);
        REAL s = 0.;
        t0 = std::chrono::steady_clock::now();
        for (int64_t i = 0; i < bench; i++) {
            TPZManVector<REAL, 3> X = {xs[2 * i], xs[2 * i + 1], 0.};
            REAL F[2];
            f(X, F);
            s += F[0] + F[1];
        }
        const double t = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        std::cout << std::setprecision(3) << "evaluation: " << 1.e9 * t / bench << " ns per point (ForceField, " << bench
                  << " points in x in [-3H, x_T + 3H], y in [-4H, 0]; checksum " << s << ")\n";
    }
    return 0;
}

// ------------------------------------------------------------------------------------------------------------------
// fig5 out=, fig8, fig9: production runs of the paper's Figs. 5, 8 and 9 (limit analysis with both seepage fields)
// into resumable CSVs (default results/cpp/): the columns of data/paper_fig*.csv, then the mechanism, run information
// and a settings string; a progress log next to each CSV (.log)
// ------------------------------------------------------------------------------------------------------------------

/// directory of this source file (data/ and results/ are looked up there)
std::string ProjectDir() {
    const std::string src(__FILE__);
    return src.substr(0, src.find_last_of('/') + 1);
}

/// number of the figure CSVs: prec significant digits, "inf" / "-inf" / "nan" when not finite
std::string CsvNum(REAL v, int prec = 10) {
    if (!std::isfinite(v)) return v > 0. ? "inf" : (v < 0. ? "-inf" : "nan");
    std::ostringstream s;
    s << std::setprecision(prec) << v;
    return s.str();
}

/// Resumable CSV: the header is written when the file is new or empty (an existing file must have the same header);
/// a row is identified by its key columns (indices into the header, the settings string being one of them). Rows are
/// appended and flushed one at a time, so a run killed at any point loses at most the case in progress; a row cut by
/// the kill (wrong number of columns) is ignored, and a missing final newline is restored before the next append.
class ResumableCSV {
    std::string fFile;
    std::vector<size_t> fKey;
    size_t fNCols = 0;
    std::set<std::string> fDone;
    bool fNeedNewline = false;
    int fIgnored = 0;

public:
    static std::vector<std::string> Split(const std::string &line) {
        std::vector<std::string> v;
        std::string item;
        std::istringstream ss(line);
        while (std::getline(ss, item, ',')) v.push_back(item);
        if (!line.empty() && line.back() == ',') v.push_back("");
        return v;
    }
    std::string Key(const std::vector<std::string> &cols) const {
        std::string k;
        for (size_t i : fKey) k += (i < cols.size() ? cols[i] : std::string("?")) + "|";
        return k;
    }
    ResumableCSV(const std::string &file, const std::string &header, const std::vector<size_t> &key)
        : fFile(file), fKey(key), fNCols(Split(header).size()) {
        const std::filesystem::path dir = std::filesystem::path(file).parent_path();
        if (!dir.empty()) std::filesystem::create_directories(dir);
        bool haveHeader = false;
        {
            std::ifstream in(file);
            std::string line;
            while (std::getline(in, line)) {
                if (!line.empty() && line.back() == '\r') line.pop_back();
                if (line.empty()) continue;
                if (!haveHeader) {
                    if (line != header) {
                        std::cerr << file << ": the header differs from this version's\n  " << header << "\n(use another out=)\n";
                        exit(1);
                    }
                    haveHeader = true;
                    continue;
                }
                const std::vector<std::string> cols = Split(line);
                if (cols.size() == fNCols) fDone.insert(Key(cols));
                else fIgnored++;
            }
        }
        if (!haveHeader) std::ofstream(file, std::ios::trunc) << header << "\n";
        else { // the last row may have been cut by a kill: the next one starts on a new line
            std::ifstream f(file, std::ios::binary | std::ios::ate);
            const std::streamoff size = f.tellg();
            if (size > 0) {
                f.seekg(size - 1);
                char c = '\n';
                f.get(c);
                fNeedNewline = c != '\n';
            }
        }
    }
    const std::string &File() const { return fFile; }
    int Ignored() const { return fIgnored; }
    size_t NDone() const { return fDone.size(); }
    bool Done(const std::vector<std::string> &cols) const { return fDone.count(Key(cols)) > 0; }
    void Append(const std::vector<std::string> &cols) {
        if (cols.size() != fNCols) {
            std::cerr << "ResumableCSV: " << cols.size() << " columns for a header of " << fNCols << "\n";
            DebugStop();
        }
        std::string line = fNeedNewline ? "\n" : "";
        for (size_t i = 0; i < cols.size(); i++) line += (i ? "," : "") + cols[i];
        std::ofstream f(fFile, std::ios::app);
        f << line << "\n";
        f.flush();
        fNeedNewline = false;
        fDone.insert(Key(cols));
    }
};

/// progress log: each message to stdout and, with a UTC time stamp, appended to a file
struct ProgressLog {
    std::string file;
    void operator()(const std::string &msg) const {
        const std::time_t t = std::time(nullptr);
        char stamp[32];
        std::strftime(stamp, sizeof stamp, "%Y-%m-%d %H:%M:%S", std::gmtime(&t));
        std::cout << msg << std::endl;
        if (!file.empty()) std::ofstream(file, std::ios::app) << stamp << "  " << msg << "\n";
    }
};

/// the .log file next to a .csv
std::string LogFileOf(const std::string &csv) {
    const size_t n = csv.size();
    return (n > 4 && csv.compare(n - 4, 4, ".csv") == 0 ? csv.substr(0, n - 4) : csv) + ".log";
}

/// limit-analysis settings of the figure drivers (options and defaults as for 'la', but hydraulic mesh href = 1)
struct FigureSettings {
    la::Settings st;
    int href = 1;
    REAL lm = 10., dmax = 10.;
    std::vector<double> boxm = {50., 10., 30.}; ///< hydraulic box in metres: left of O, right of T, below T

    /// settings string of the CSV rows (no commas): what changes a result besides the row's own columns
    std::string Describe(const Problem &p) const {
        std::ostringstream s;
        s << std::setprecision(10) << "box_m=" << boxm[0] << "/" << boxm[1] << "/" << boxm[2] << ";hbc="
          << FarSides::Name(p.hy.far.left) << "/" << FarSides::Name(p.hy.far.bottom) << "/" << FarSides::Name(p.hy.far.right)
          << ";horder=" << p.hy.order << ";href=" << href << ";hmesh=" << p.hmesh << "/" << p.hsize.h0 << "/" << p.hsize.hs << "/"
          << p.hsize.grade << "/" << p.hsize.hmax << ";lm=" << lm << ";mech=";
        for (size_t i = 0; i < st.kinds.size(); i++) s << (i ? "+" : "") << (st.kinds[i] == 1 ? "I" : "II");
        s << ";seeds=";
        for (size_t i = 0; i < st.seeds.size(); i++) s << (i ? "+" : "") << st.seeds[i];
        s << ";np=" << st.nParticles << ";niter=" << st.nIter << ";pools=" << st.extraPools << ";q=" << st.quadSearch << "/"
          << st.quadFinal << ";polish=" << st.polish << ";dmax=" << dmax;
        return s.str();
    }
};

FigureSettings ReadFigureSettings(Options &o) {
    FigureSettings fs;
    fs.boxm = o.NumList("hboxm", fs.boxm);
    fs.href = o.IntList("href", {1})[0];
    fs.lm = o.Num("lm", 10.);
    fs.dmax = o.Num("dmax", 10.);
    const std::string mech = o.Str("mech", "I,II");
    fs.st.kinds.clear();
    if (mech == "I" || mech == "I,II" || mech == "II,I") fs.st.kinds.push_back(1);
    if (mech == "II" || mech == "I,II" || mech == "II,I") fs.st.kinds.push_back(2);
    fs.st.seeds.clear();
    for (int s : o.IntList("seeds", {0, 1, 2})) fs.st.seeds.push_back(uint64_t(s));
    fs.st.nParticles = int(o.Num("np", 40));
    fs.st.nIter = int(o.Num("niter", 150));
    fs.st.quadSearch = o.Str("qsearch", "coarse");
    fs.st.quadFinal = o.Str("qfinal", "fine");
    fs.st.polish = o.Num("polish", 1) != 0.;
    fs.st.extraPools = int(o.Num("pools", 25));
    fs.st.nThreads = int(o.Num("threads", 2));
    if (fs.boxm.size() != 3 || fs.st.kinds.empty() || !la::Quadrature::Level(fs.st.quadSearch).Valid() ||
        !la::Quadrature::Level(fs.st.quadFinal).Valid()) {
        std::cerr << "figure settings: hboxm=<left>,<right>,<depth> (m), mech=I|II|I,II, qsearch/qfinal = coarse, medium, "
                     "fine, xfine, ref or dense\n";
        exit(1);
    }
    return fs;
}

/// the problem of one figure case: slope height H (m), angle beta (deg), water level h_w = hwr H, hydraulic box fixed
/// in metres (boxm: left of O, right of T, below T) whatever H (the box of the paper's Fig. 4, see SPEC / README)
Problem FigureProblem(const Problem &base, REAL H, REAL beta, REAL hwr, const std::vector<double> &boxm) {
    Problem p = base;
    p.hyd.H = H;
    for (int i = 0; i < 3; i++) p.hext[i] = boxm[i] / H;
    p.SetSlope(beta, hwr);
    return p;
}

/// one limit analysis of a figure: water = fe (-grad u'_FE), analytical (K^-1 v'_opt) or none (h_w = 0: f = 0 with
/// the buoyant gamma' = gamma - gamma_w); phi in degrees
struct FigureCase {
    la::Result r;
    LAField lf;
};

FigureCase RunFigureCase(const Problem &p, REAL phiDeg, const std::string &water, const FigureSettings &fs) {
    FigureCase fc;
    fc.lf = MakeLAField(p, water, fs.href, fs.lm);
    const la::Problem prob(p.hyd.beta, p.hyd.H, p.soil.c, phiDeg, p.soil.gamma, p.hy.gammaw, fc.lf.seep, "auto", fs.dmax);
    fc.r = la::StabilityFactor(prob, fs.st);
    return fc;
}

/// mechanism and run columns of the figure CSVs (after the paper's columns; A, B, C in paper coordinates, metres)
const std::string kFigureMechHeader = "mechanism,theta1,theta2,eta_or_dH,L_over_H,Ax,Bx,By,Cx,Cy,r0,P_mr,P_gamma,P_u,"
                                      "seed_spread,at_bound,n_eval,Jnorm,t_field_s,t_la_s";
constexpr size_t kFigureMechCols = 20;

std::vector<std::string> FigureMechColumns(const FigureCase &fc, REAL H) {
    const la::Result &r = fc.r;
    std::vector<std::string> c(16, "");
    if (r.found) {
        std::string atb;
        for (const std::string &b : r.atBound) atb += (atb.empty() ? "" : ";") + b;
        if (r.dmaxReached) atb += (atb.empty() ? "" : ";") + std::string("d=dmax");
        c = {r.kind == 1 ? "I" : "II", CsvNum(r.x[0]), CsvNum(r.x[1]), CsvNum(r.kind == 1 ? r.x[2] : r.x[2] - 1.),
             CsvNum(r.mech.L / H), CsvNum(r.A[0]), CsvNum(r.B[0]), CsvNum(r.B[1]),
             CsvNum(r.C[0]), CsvNum(r.C[1]), CsvNum(r.mech.r0), CsvNum(r.P.Pmr), CsvNum(r.P.Pgamma), CsvNum(r.P.Pu),
             CsvNum(r.seedSpread, 3), atb};
    }
    c.push_back(std::to_string(r.nEval));
    c.push_back(CsvNum(fc.lf.Jnorm));
    c.push_back(CsvNum(fc.lf.seconds, 4));
    c.push_back(CsvNum(r.seconds, 4));
    return c;
}

/// one-line summary of a figure case for the progress log
std::string FigureCaseSummary(const FigureCase &fc) {
    const la::Result &r = fc.r;
    std::ostringstream s;
    if (!r.found) s << "Gamma = inf (no mechanism with P_gamma + P_u > 0)";
    else
        s << std::setprecision(7) << "Gamma = " << r.Gamma << ", H_crit = " << r.Hcrit << " m, mechanism "
          << (r.kind == 1 ? "I eta = " : "II d/H = ") << std::setprecision(4) << (r.kind == 1 ? r.x[2] : r.x[2] - 1.)
          << ", seed spread " << std::setprecision(2) << r.seedSpread;
    for (const std::string &b : r.atBound) s << ", " << b;
    if (r.found && r.dmaxReached) s << ", d = d_max";
    s << std::setprecision(3) << " (" << fc.lf.seconds + r.seconds << " s)";
    return s.str();
}

std::vector<std::string> SplitList(const std::string &s) {
    std::vector<std::string> v;
    std::string item;
    std::istringstream ss(s);
    while (std::getline(ss, item, ','))
        if (!item.empty()) v.push_back(item);
    return v;
}

/// the curves of a figure: vopt (K^-1 v'_opt, water=analytical) and FE (-grad u'_FE, water=fe)
std::vector<std::string> FigureCurves(Options &o) {
    const std::vector<std::string> cv = SplitList(o.Str("curves", "vopt,FE"));
    for (const std::string &c : cv)
        if (c != "vopt" && c != "FE") {
            std::cerr << "curves must be vopt, FE or vopt,FE\n";
            exit(1);
        }
    if (cv.empty()) exit(1);
    return cv;
}

/// soil parameters of a Fig. 8 panel ("London" or "Israeli"): Table 1 (soil=table1) or with the (c, phi) pairs
/// exchanged between the panels (soil=swapped, the set that reproduces the h_w = 0 ends of Fig. 8); gamma = 18
void Fig8Soil(const std::string &panel, const std::string &set, REAL &c, REAL &phiDeg, REAL &gamma) {
    const bool london = panel == "London";
    const bool table = set == "table1";
    gamma = 18.;
    if (london == table) c = 6., phiDeg = 32.; // London clay of Table 1
    else c = 11.7, phiDeg = 24.7;              // Israeli clay of Table 1
}

/// fig8: H_crit versus h_w / H of the paper's Fig. 8 (alpha = 1): London clay panel (beta = 30, 60 deg) and Israeli
/// clay panel (beta = 35, 60 deg), curves vopt (K^-1 v'_opt) and FE (-grad u'_FE, hydraulic box fixed in metres,
/// default 50 / 10 / 30 m with u = 0 on the left side and the base), and their common h_w = 0 end (f = 0 with the
/// buoyant gamma', both rows from one run). H_crit = Gamma(H) H (similarity) with the reference height H = 1 m, at
/// which the box is the paper's 50 / 10 / 30 H. soil=swapped (default): the (c, phi) pairs of Table 1 exchanged
/// between the panels; soil=table1: as printed. gammaw=9.8 by default (it reproduces the h_w = 0 ends, see README).
int CmdFig8(Options &o) {
    for (const char *k : {"beta", "hw", "c", "phi", "gamma"})
        if (o.Has(k)) {
            std::cerr << "fig8: " << k << " is set by the panels (panels=, hws=, soil=)\n";
            return 1;
        }
    const auto t0 = std::chrono::steady_clock::now();
    Problem base = ReadProblem(o);
    const REAL H = o.Num("H", 1.);
    base.hy.gammaw = o.Num("gammaw", 9.8);
    const std::string soil = o.Str("soil", "swapped");
    const std::vector<std::string> panels = SplitList(o.Str("panels", "London30,London60,Israeli35,Israeli60"));
    const std::vector<double> hws = o.NumList("hws", {0., 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.});
    const std::vector<std::string> curves = FigureCurves(o);
    const FigureSettings fs = ReadFigureSettings(o);
    const std::string out = o.Str("out", ProjectDir() + "results/cpp/" + (soil == "table1" ? "fig8_table1.csv" : "fig8.csv"));
    o.CheckUnused();
    if (soil != "swapped" && soil != "table1") {
        std::cerr << "fig8: soil must be swapped or table1\n";
        return 1;
    }
    std::vector<std::pair<std::string, REAL>> pan;
    for (const std::string &s : panels) {
        const std::string name = s.rfind("London", 0) == 0 ? "London" : (s.rfind("Israeli", 0) == 0 ? "Israeli" : "");
        if (name.empty() || s.size() == name.size()) {
            std::cerr << "fig8: panels are London<beta> or Israeli<beta>, e.g. London30,Israeli35\n";
            return 1;
        }
        pan.push_back({name, atof(s.c_str() + name.size())});
    }
    const REAL alpha = base.hy.kh / base.hy.kv, gw = base.hy.gammaw;
    const std::string settings = fs.Describe(base) + ";soil=" + soil;
    const std::string header = "soil,beta_deg,curve,hw_over_H,Hcrit_m,Gamma,H_ref,c,phi_deg,gamma,gamma_w,alpha,water," +
                               kFigureMechHeader + ",settings";
    const size_t nCols = 13 + kFigureMechCols + 1;
    ResumableCSV csv(out, header, {0, 1, 2, 3, 6, 7, 8, 9, 10, 11, nCols - 1});
    const ProgressLog log{LogFileOf(out)};
    std::ostringstream intro;
    intro << "fig8: " << pan.size() << " panels x " << hws.size() << " h_w / H x " << curves.size() << " curves, soil " << soil
          << ", gamma_w " << gw << ", alpha " << alpha << ", H_ref " << H << " m; " << settings << "; out " << out << " ("
          << csv.NDone() << " rows already there" << (csv.Ignored() ? ", " + std::to_string(csv.Ignored()) + " cut rows ignored" : "")
          << ")";
    log(intro.str());
    int run = 0, skipped = 0, inf = 0;
    for (const auto &pb : pan)
        for (double hw : hws) {
            REAL c, phi, gamma;
            Fig8Soil(pb.first, soil, c, phi, gamma);
            Problem p = FigureProblem(base, H, pb.second, hw, fs.boxm);
            p.soil.c = c, p.soil.gamma = gamma, p.soil.phi = phi * M_PI / 180.;
            // h_w = 0: one run without seepage gives the common end of both curves
            std::vector<std::pair<std::vector<std::string>, std::string>> todo; // (curves of the run, water)
            if (hw > 0.)
                for (const std::string &cv : curves) todo.push_back({{cv}, cv == "vopt" ? "analytical" : "fe"});
            else todo.push_back({curves, "none"});
            for (auto &job : todo) {
                std::vector<std::vector<std::string>> rows;
                bool all = true;
                for (const std::string &cv : job.first) {
                    std::vector<std::string> row = {pb.first, CsvNum(pb.second), cv, CsvNum(hw), "", "", CsvNum(H), CsvNum(c),
                                                    CsvNum(phi), CsvNum(gamma), CsvNum(gw), CsvNum(alpha), job.second};
                    row.resize(nCols - 1, "");
                    row.push_back(settings);
                    all = all && csv.Done(row);
                    rows.push_back(row);
                }
                std::ostringstream id;
                id << pb.first << " beta=" << pb.second << " hw/H=" << hw << " " << (hw > 0. ? job.first[0] : "none");
                if (all) {
                    skipped++;
                    continue;
                }
                const FigureCase fc = RunFigureCase(p, phi, job.second, fs);
                const std::vector<std::string> mc = FigureMechColumns(fc, H);
                for (auto &row : rows) {
                    if (csv.Done(row)) continue;
                    row[4] = CsvNum(fc.r.Hcrit), row[5] = CsvNum(fc.r.Gamma);
                    std::copy(mc.begin(), mc.end(), row.begin() + 13);
                    csv.Append(row);
                }
                run++;
                inf += !fc.r.found;
                log("  [fig8 " + std::to_string(run) + "] " + id.str() + ": " + FigureCaseSummary(fc));
            }
        }
    std::ostringstream end;
    end << "fig8: " << run << " limit analyses run (" << inf << " with Gamma = inf), " << skipped << " already in " << out
        << ", " << std::setprecision(4) << std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count() << " s";
    log(end.str());
    return 0;
}

/// fig9: stability factor Gamma versus beta of the paper's Fig. 9 (H = 5 m, c = 10 kPa, phi = 30 deg, gamma = 20,
/// h_w = H, gamma_w = 9.81 by default), alphas=1,5,10, curves vopt (K^-1 v'_opt) and FE (-grad u'_FE, hydraulic box
/// fixed in metres, default 50 / 10 / 30 m = 10 / 2 / 6 H, u = 0 on the left side and the base); betas as in
/// scripts/reproduce_python.py (every 5 deg from 15 to 90 plus 37.5, 52.5, 67.5, 82.5)
int CmdFig9(Options &o) {
    for (const char *k : {"beta", "alpha"})
        if (o.Has(k)) {
            std::cerr << "fig9: use betas= and alphas=\n";
            return 1;
        }
    const auto t0 = std::chrono::steady_clock::now();
    Problem base = ReadProblem(o); // defaults H = 5, c = 10, phi = 30, gamma = 20, gamma_w = 9.81, h_w / H = 1
    const REAL phi = o.Num("phi", 30.), H = base.hyd.H, hwr = base.hyd.hw / base.hyd.H;
    const std::vector<double> alphas = o.NumList("alphas", {1., 5., 10.});
    const std::vector<double> betas = o.NumList("betas", {15., 20., 25., 30., 35., 37.5, 40., 45., 50., 52.5, 55., 60., 65.,
                                                          67.5, 70., 75., 80., 82.5, 85., 90.});
    const std::vector<std::string> curves = FigureCurves(o);
    const FigureSettings fs = ReadFigureSettings(o);
    const std::string out = o.Str("out", ProjectDir() + "results/cpp/fig9.csv");
    o.CheckUnused();
    const REAL c = base.soil.c, gamma = base.soil.gamma, gw = base.hy.gammaw;
    const std::string settings = fs.Describe(base);
    const std::string header = "alpha,beta_deg,curve,Gamma,Hcrit_m,H,c,phi_deg,gamma,gamma_w,hw_over_H,water," +
                               kFigureMechHeader + ",settings";
    const size_t nCols = 12 + kFigureMechCols + 1;
    ResumableCSV csv(out, header, {0, 1, 2, 5, 6, 7, 8, 9, 10, nCols - 1});
    const ProgressLog log{LogFileOf(out)};
    std::ostringstream intro;
    intro << "fig9: " << alphas.size() << " alphas x " << betas.size() << " betas x " << curves.size() << " curves, H " << H
          << " m, c " << c << ", phi " << phi << ", gamma " << gamma << ", gamma_w " << gw << ", h_w / H " << hwr << "; "
          << settings << "; out " << out << " (" << csv.NDone() << " rows already there"
          << (csv.Ignored() ? ", " + std::to_string(csv.Ignored()) + " cut rows ignored" : "") << ")";
    log(intro.str());
    int run = 0, skipped = 0, inf = 0;
    for (double alpha : alphas)
        for (double beta : betas)
            for (const std::string &cv : curves) {
                const std::string water = hwr > 0. ? (cv == "vopt" ? "analytical" : "fe") : "none";
                std::vector<std::string> row = {CsvNum(alpha), CsvNum(beta), cv, "", "", CsvNum(H), CsvNum(c), CsvNum(phi),
                                                CsvNum(gamma), CsvNum(gw), CsvNum(hwr), water};
                row.resize(nCols - 1, "");
                row.push_back(settings);
                if (csv.Done(row)) {
                    skipped++;
                    continue;
                }
                Problem p = FigureProblem(base, H, beta, hwr, fs.boxm);
                p.hy.kh = alpha * p.hy.kv;
                const FigureCase fc = RunFigureCase(p, phi, water, fs);
                row[3] = CsvNum(fc.r.Gamma), row[4] = CsvNum(fc.r.Hcrit);
                const std::vector<std::string> mc = FigureMechColumns(fc, H);
                std::copy(mc.begin(), mc.end(), row.begin() + 12);
                csv.Append(row);
                run++;
                inf += !fc.r.found;
                std::ostringstream id;
                id << "alpha=" << alpha << " beta=" << beta << " " << cv;
                log("  [fig9 " + std::to_string(run) + "] " + id.str() + ": " + FigureCaseSummary(fc));
            }
    std::ostringstream end;
    end << "fig9: " << run << " limit analyses run (" << inf << " with Gamma = inf), " << skipped << " already in " << out
        << ", " << std::setprecision(4) << std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count() << " s";
    log(end.str());
    return 0;
}

/// fig5 out=<csv>: the Fig. 5 functionals (h_w = H) J(u'_FE) / (k_h H^2 gamma_w^2) (curve J_FE) and
/// -J*(v'_opt) / (k_h H^2 gamma_w^2) (curve minus_Jstar_opt) on the grid of scripts/reproduce_python.py (betas 15 to 90
/// every 7.5 deg, href=1 by default; hydraulic box in units of H, the paper's 50 / 10 / 30 m at H = 1 m), appended to
/// a resumable CSV: the columns of data/paper_fig5.csv, then m and the degenerate flag of the analytical optimum, the
/// number of equations of the FE solution, the time and the settings
int CmdFig5CSV(Options &o) {
    const auto t0 = std::chrono::steady_clock::now();
    Problem p = ReadProblem(o);
    const std::vector<double> alphas = o.NumList("alphas", {1., 2., 4., 10.});
    std::vector<double> def;
    for (int k = 0; k <= 10; k++) def.push_back(15. + 7.5 * k);
    const std::vector<double> betas = o.NumList("betas", def);
    const int href = o.IntList("href", {1})[0];
    const bool fe = o.Num("fe", 1.) != 0.;
    const REAL lm = o.Num("lm", 10.);
    const std::string out = o.Str("out", "");
    o.Str("paper", ""); // the comparison with the paper is made by scripts/plot_results.py
    o.CheckUnused();
    if (out.empty()) {
        std::cerr << "fig5: out=<csv>\n";
        return 1;
    }
    const REAL hwr = p.hyd.hw / p.hyd.H;
    std::ostringstream st;
    st << std::setprecision(10) << "hw=" << hwr << ";box_H=" << p.hext[0] << "/" << p.hext[1] << "/" << p.hext[2] << ";hbc="
       << FarSides::Name(p.hy.far.left) << "/" << FarSides::Name(p.hy.far.bottom) << "/" << FarSides::Name(p.hy.far.right)
       << ";horder=" << p.hy.order << ";href=" << href << ";hmesh=" << p.hmesh << "/" << p.hsize.h0 << "/" << p.hsize.hs << "/"
       << p.hsize.grade << "/" << p.hsize.hmax << ";lm=" << lm;
    const std::string settings = st.str();
    ResumableCSV csv(out, "alpha,beta_deg,curve,value,m,degenerate,neq,time_s,settings", {0, 1, 2, 8});
    const ProgressLog log{LogFileOf(out)};
    log("fig5: " + std::to_string(alphas.size()) + " alphas x " + std::to_string(betas.size()) + " betas; " + settings +
        "; out " + out + " (" + std::to_string(csv.NDone()) + " rows already there)");
    int run = 0, skipped = 0;
    for (double alpha : alphas)
        for (double beta : betas) {
            p.SetSlope(beta, hwr);
            p.hy.kh = alpha * p.hy.kv;
            std::vector<std::string> rowA = {CsvNum(alpha), CsvNum(beta), "minus_Jstar_opt", "", "", "", "", "", settings};
            std::vector<std::string> rowF = {CsvNum(alpha), CsvNum(beta), "J_FE", "", "", "", "", "", settings};
            std::ostringstream msg;
            msg << "  [fig5] alpha=" << alpha << " beta=" << beta << std::setprecision(8);
            bool any = false;
            if (!csv.Done(rowA)) {
                const auto ta = std::chrono::steady_clock::now();
                const auto an = MakeAnalyticalField(p, lm);
                rowA[3] = CsvNum(-an->JstarNormalized()), rowA[4] = CsvNum(an->M()), rowA[5] = an->Degenerate() ? "1" : "0";
                rowA[7] = CsvNum(std::chrono::duration<double>(std::chrono::steady_clock::now() - ta).count(), 4);
                csv.Append(rowA);
                msg << ": -J*/(kh H^2 gw^2) = " << rowA[3] << " (m = " << rowA[4] << (an->Degenerate() ? ", degenerate" : "") << ")";
                any = true;
            }
            if (fe && !csv.Done(rowF)) {
                TPZGeoMesh *gmesh = HydraulicMesh(p, href);
                const SeepageResult r = DrawdownSeepage(gmesh, p.hyd, p.hy);
                delete gmesh;
                rowF[3] = CsvNum(r.Jnorm), rowF[6] = std::to_string(r.neq), rowF[7] = CsvNum(r.seconds, 4);
                csv.Append(rowF);
                msg << "; J_FE/(kh H^2 gw^2) = " << rowF[3] << " (" << r.neq << " equations, " << std::setprecision(3) << r.seconds << " s)";
                any = true;
            }
            if (any) run++, log(msg.str());
            else skipped++;
        }
    std::ostringstream end;
    end << "fig5: " << run << " cases run, " << skipped << " already in " << out << ", " << std::setprecision(4)
        << std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count() << " s";
    log(end.str());
    return 0;
}

// ------------------------------------------------------------------------------------------------------------------
// fembatch: FEM gravity-increase stability factor (FEMStability.h: elastoplastic Mohr-Coulomb, refinement cycles of
// the plastic zone) of cases of the paper's Figs. 8 and 9, an independent check of the limit-analysis curves; each
// row also holds the limit analysis of the same case (figure settings, the same seepage field). Resumable CSVs
// results/cpp/fem_fig9.csv and fem_fig8.csv (one row per case; cases already there with the same settings are
// skipped), progress in the .log next to each CSV.
// ------------------------------------------------------------------------------------------------------------------

/// FEM settings of fembatch; defaults = production settings of the convergence study (README, section (g))
struct FEMBatchSettings {
    DriverSettings ds;
    int nref = 3;    ///< refinement cycles of the plastic zone (nref + 1 gravity-increase analyses)
    REAL sa = 2.;    ///< stability box: sa H + H / tan(beta) left of O, right of T and below T, capped at the hydraulic box
    MeshSize ssize;  ///< initial stability mesh (units of H)
    FEMBatchSettings() {
        ds.markFrac = 0.05;
        ssize.h0 = ssize.hs = 0.25, ssize.grade = 0.25, ssize.hmax = 1.;
    }
    /// settings string of the CSV rows (no commas)
    std::string Describe() const {
        std::ostringstream s;
        s << std::setprecision(10) << "fem:nref=" << nref << ";mark=" << ds.markFrac << ";smesh=" << ssize.h0 << "/" << ssize.hs
          << "/" << ssize.grade << "/" << ssize.hmax << ";sa=" << sa << ";maxnewton=" << ds.maxNewton << ";tolfs=" << ds.tolFS;
        return s.str();
    }
};

FEMBatchSettings ReadFEMBatchSettings(Options &o) {
    FEMBatchSettings fs;
    fs.nref = int(o.Num("nref", fs.nref));
    fs.sa = o.Num("sa", fs.sa);
    fs.ds.markFrac = o.Num("mark", fs.ds.markFrac);
    fs.ds.maxNewton = int(o.Num("maxnewton", fs.ds.maxNewton));
    fs.ds.tolFS = o.Num("tolfs", fs.ds.tolFS);
    fs.ssize.h0 = o.Num("sh0", fs.ssize.h0);
    fs.ssize.hs = o.Num("shs", fs.ssize.hs);
    fs.ssize.grade = o.Num("sgrade", fs.ssize.grade);
    fs.ssize.hmax = o.Num("shmax", fs.ssize.hmax);
    if (fs.nref < 0 || fs.sa <= 0. || fs.ds.markFrac <= 0. || fs.ds.markFrac >= 1. || fs.ds.maxNewton < 1 || fs.ds.tolFS <= 0.) {
        std::cerr << "fembatch: nref >= 0, sa > 0, 0 < mark < 1, maxnewton >= 1, tolfs > 0\n";
        exit(1);
    }
    return fs;
}

/// stability box of a figure case: sa H + H / tan(beta) left of O, right of T and below T (as fs with sa), each side
/// capped at the hydraulic box of p (the FE field is zero outside it; the vopt curve uses the same domain)
void FEMStabilityBox(Problem &p, const FEMBatchSettings &fs) {
    p.ssize = fs.ssize;
    p.smesh = "gen";
    p.sref = 0;
    p.stab = p.hyd;
    StabilityExtents(p.stab, fs.sa);
    p.stab.left = std::min(p.stab.left, p.hyd.left);
    p.stab.right = std::min(p.stab.right, p.hyd.right);
    p.stab.depth = std::min(p.stab.depth, p.hyd.depth);
}

/// one FEM gravity-increase analysis of a figure case: water = fe (-grad u'_FE on the hydraulic box of p, href
/// uniform refinements), analytical (K^-1 v'_opt, L_m = lm H) or none (h_w = 0: f = 0); b = lambda (gamma' g + f)
/// with gamma' = gamma - gamma_w in the three cases. nout > 0: integration points of the stability mesh outside the
/// hydraulic mesh (not run)
struct FEMCase {
    std::vector<FSCycle> cyc;
    FSExtrapolation ex;
    int64_t nout = 0;
    double tField = 0., tFEM = 0.;
    REAL Lambda() const { return cyc.empty() ? NAN : cyc.back().gi; } ///< lambda of the last cycle
    /// Gamma_FEM: order-1 extrapolation of the last two cycles to h -> 0 (the last lambda when there is one cycle)
    REAL Estimate() const { return cyc.size() >= 2 ? ex.h1 : Lambda(); }
};

FEMCase RunFEMCase(const Problem &p, const std::string &water, const FEMBatchSettings &fs, int href, REAL lm) {
    FEMCase fc;
    auto t0 = std::chrono::steady_clock::now();
    ForceField f = NoSeepage();
    std::shared_ptr<const PoreField> field;
    if (water == "fe") {
        TPZGeoMesh *hmesh = HydraulicMesh(p, href);
        const SeepageResult r = DrawdownSeepage(hmesh, p.hyd, p.hy);
        delete hmesh;
        field = r.field;
        f = r.field->AsForceField(r.field);
    } else if (water == "analytical") f = AnalyticalForceField(MakeAnalyticalField(p, lm));
    else if (water != "none") DebugStop();
    fc.tField = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    t0 = std::chrono::steady_clock::now();
    TPZGeoMesh *smesh = StabilityMesh(p);
    const TMCVoigt model = ModelVoigt(p.soil);
    if (field) {
        TPZCompMesh *cmesh = CreateCMesh(smesh, 2, model, p.soil);
        int64_t nneg, ntot;
        REAL pmin;
        NegativePressurePoints(cmesh, *field, p.hy.gammaw, nneg, ntot, pmin, fc.nout);
        delete cmesh;
    }
    if (fc.nout == 0) {
        fc.cyc = GravityIncreaseFS(smesh, model, p.soil, f, p.soil.gamma - p.hy.gammaw, fs.nref, false, "", fs.ds);
        fc.ex = ExtrapolateCycles(fc.cyc);
    }
    delete smesh;
    fc.tFEM = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    return fc;
}

/// FEM columns of the fembatch CSVs (after the case columns): cycles, equations of the last one, lambda and the
/// equations of every cycle (';' lists), extrapolations, extent of the 1 % plastic zone of the last cycle (units of
/// H: x_min / H, (x_max - x_T) / H, (y_min - y_T) / H), stability box (left of O, right of T, below T, / H), times
const std::string kFEMHeader = "ncycles,neq,lambda_cycles,neq_cycles,lambda_richardson,order_richardson,lambda_sqrtneq,"
                               "zone1_H,stab_box_H,t_field_s,t_fem_s,t_la_s";
constexpr size_t kFEMCols = 12;

std::vector<std::string> FEMColumns(const FEMCase &fc, const Problem &p, double tLA) {
    const REAL H = p.stab.H;
    std::string lam, neq, zone, box;
    for (size_t i = 0; i < fc.cyc.size(); i++) {
        lam += (i ? ";" : "") + CsvNum(fc.cyc[i].gi, 8);
        neq += (i ? ";" : "") + std::to_string(fc.cyc[i].neq);
    }
    if (!fc.cyc.empty()) {
        const FSCycle &c = fc.cyc.back();
        zone = CsvNum(c.zone1[0] / H, 4) + ";" + CsvNum((c.zone1[2] - p.stab.XT()) / H, 4) + ";" + CsvNum((c.zone1[1] + H) / H, 4);
    }
    box = CsvNum(p.stab.left / H, 6) + ";" + CsvNum(p.stab.right / H, 6) + ";" + CsvNum(p.stab.depth / H, 6);
    return {std::to_string(fc.cyc.size()), fc.cyc.empty() ? "" : std::to_string(fc.cyc.back().neq), lam, neq,
            CsvNum(fc.ex.richardson), CsvNum(fc.ex.order, 4), CsvNum(fc.ex.sqrtneq), zone, box, CsvNum(fc.tField, 4),
            CsvNum(fc.tFEM, 5), CsvNum(tLA, 4)};
}

/// fembatch fig=9: Fig. 9 data (H = 5 m, c = 10, phi = 30, gamma = 20, gamma_w = 9.81, h_w = H), alphas=1,5,10
/// betas=30,45,60,75,90 curves=FE,vopt -> results/cpp/fem_fig9.csv: Gamma_FEM = order-1 extrapolation of the last two
/// cycles to h -> 0 (FEMCase::Estimate), the lambda of the last cycle and of every cycle, the other extrapolations, and
/// Gamma_LA of the limit analysis of the same case (FigureSettings, threads=2)
int FEMBatchFig9(Options &o, const FEMBatchSettings &fem, const FigureSettings &la, const std::string &outDefault) {
    Problem base = ReadProblem(o);
    const REAL phi = o.Num("phi", 30.), H = base.hyd.H, hwr = base.hyd.hw / base.hyd.H;
    const std::vector<double> alphas = o.NumList("alphas", {1., 5., 10.});
    const std::vector<double> betas = o.NumList("betas", {30., 45., 60., 75., 90.});
    const std::vector<std::string> curves = FigureCurves(o);
    const std::string out = o.Str("out", outDefault);
    o.CheckUnused();
    const auto t0 = std::chrono::steady_clock::now();
    const REAL c = base.soil.c, gamma = base.soil.gamma, gw = base.hy.gammaw;
    const std::string settings = fem.Describe() + "|la:" + la.Describe(base);
    const std::string header = "alpha,beta_deg,curve,Gamma_FEM,Gamma_FEM_last,Gamma_LA,FEM_over_LA_minus_1,H,c,phi_deg,gamma,"
                               "gamma_w,hw_over_H,water," + kFEMHeader + ",settings";
    const size_t nCols = 14 + kFEMCols + 1;
    ResumableCSV csv(out, header, {0, 1, 2, 7, 8, 9, 10, 11, 12, nCols - 1});
    const ProgressLog log{LogFileOf(out)};
    std::ostringstream intro;
    intro << "fembatch fig=9: " << alphas.size() << " alphas x " << betas.size() << " betas x " << curves.size() << " curves, H "
          << H << " m, c " << c << ", phi " << phi << ", gamma " << gamma << ", gamma_w " << gw << ", h_w / H " << hwr << "; "
          << settings << "; out " << out << " (" << csv.NDone() << " rows already there)";
    log(intro.str());
    int run = 0, skipped = 0, failed = 0;
    for (double alpha : alphas)
        for (double beta : betas)
            for (const std::string &cv : curves) {
                const std::string water = hwr > 0. ? (cv == "vopt" ? "analytical" : "fe") : "none";
                std::vector<std::string> row = {CsvNum(alpha), CsvNum(beta), cv, "", "", "", "", CsvNum(H), CsvNum(c), CsvNum(phi),
                                                CsvNum(gamma), CsvNum(gw), CsvNum(hwr), water};
                row.resize(nCols - 1, "");
                row.push_back(settings);
                if (csv.Done(row)) {
                    skipped++;
                    continue;
                }
                Problem p = FigureProblem(base, H, beta, hwr, la.boxm);
                p.hy.kh = alpha * p.hy.kv;
                FEMStabilityBox(p, fem);
                std::ostringstream id;
                id << "alpha=" << alpha << " beta=" << beta << " " << cv;
                log("  [fem fig9] " + id.str() + ": limit analysis, then " + std::to_string(fem.nref + 1) + " FEM cycles");
                const FigureCase lac = RunFigureCase(p, phi, water, la);
                const FEMCase fc = RunFEMCase(p, water, fem, la.href, la.lm);
                if (fc.nout > 0 || fc.cyc.empty()) {
                    log("  [fem fig9] " + id.str() + ": " + std::to_string(fc.nout) +
                        " integration points of the stability mesh outside the hydraulic mesh, not run");
                    failed++;
                    continue;
                }
                row[3] = CsvNum(fc.Estimate()), row[4] = CsvNum(fc.Lambda()), row[5] = CsvNum(lac.r.Gamma);
                row[6] = CsvNum(fc.Estimate() / lac.r.Gamma - 1., 6);
                const std::vector<std::string> mc = FEMColumns(fc, p, lac.lf.seconds + lac.r.seconds);
                std::copy(mc.begin(), mc.end(), row.begin() + 14);
                csv.Append(row);
                run++;
                std::ostringstream s;
                s << std::setprecision(6) << "  [fem fig9 " << run << "] " << id.str() << ": Gamma_FEM = " << fc.Estimate()
                  << " (h -> 0; last cycle " << fc.Lambda() << ", " << fc.cyc.back().neq << " equations), Gamma_LA = "
                  << lac.r.Gamma << " (" << std::setprecision(4) << fc.tField + fc.tFEM << " s)";
                log(s.str());
            }
    std::ostringstream end;
    end << "fembatch fig=9: " << run << " cases run, " << skipped << " already in " << out << ", " << failed << " failed, "
        << std::setprecision(4) << std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count() << " s";
    log(end.str());
    return failed == 0 ? 0 : 1;
}

/// fembatch fig=8: Fig. 8 panels (alpha = 1, gamma = 18, (c, phi) of Table 1 swapped between the panels by default,
/// gamma_w = 9.8) panels=London30,London60,Israeli35,Israeli60 hws=0,0.2,0.5,1 curves=FE,vopt (h_w = 0: one run
/// without seepage, f = 0 with gamma', rows for both curves) -> results/cpp/fem_fig8.csv (H_crit FEM = Estimate()
/// H_ref, the last cycle, the limit analysis H_crit). Each case is computed at
/// the reference height H_ref = H_crit of its limit analysis (3 significant digits; hydraulic box 50 / 10 / 30 H_ref,
/// the paper's box at H = 1 m by similarity), so that lambda ~ 1 (the continuation of SlopeAnalysis starts with
/// steps of 0.5 and stops at 100); H_crit = lambda H_ref (exact similarity, see the check)
int FEMBatchFig8(Options &o, const FEMBatchSettings &fem, const FigureSettings &la, const std::string &outDefault) {
    for (const char *k : {"beta", "hw", "c", "phi", "gamma", "H"})
        if (o.Has(k)) {
            std::cerr << "fembatch fig=8: " << k << " is set by the panels (panels=, hws=, soil=)\n";
            return 1;
        }
    Problem base = ReadProblem(o);
    base.hy.gammaw = o.Num("gammaw", 9.8);
    const std::string soil = o.Str("soil", "swapped");
    const std::vector<std::string> panels = SplitList(o.Str("panels", "London30,London60,Israeli35,Israeli60"));
    const std::vector<double> hws = o.NumList("hws", {0., 0.2, 0.5, 1.});
    const std::vector<std::string> curves = FigureCurves(o);
    const std::string out = o.Str("out", outDefault);
    o.CheckUnused();
    if (soil != "swapped" && soil != "table1") {
        std::cerr << "fembatch fig=8: soil must be swapped or table1\n";
        return 1;
    }
    std::vector<std::pair<std::string, REAL>> pan;
    for (const std::string &s : panels) {
        const std::string name = s.rfind("London", 0) == 0 ? "London" : (s.rfind("Israeli", 0) == 0 ? "Israeli" : "");
        if (name.empty() || s.size() == name.size()) {
            std::cerr << "fembatch fig=8: panels are London<beta> or Israeli<beta>\n";
            return 1;
        }
        pan.push_back({name, atof(s.c_str() + name.size())});
    }
    const auto t0 = std::chrono::steady_clock::now();
    const REAL alpha = base.hy.kh / base.hy.kv, gw = base.hy.gammaw;
    // the hydraulic box of Fig. 8 is the paper's 50 / 10 / 30 m at H = 1 m: 50 / 10 / 30 H_ref (la.boxm / 1 m)
    const std::string settings = fem.Describe() + "|la:" + la.Describe(base) + ";soil=" + soil + ";Href=LA3";
    const std::string header = "soil,beta_deg,curve,hw_over_H,Hcrit_FEM_m,Hcrit_FEM_last_m,Hcrit_LA_m,FEM_over_LA_minus_1,H_ref,"
                               "lambda_last,c,phi_deg,gamma,gamma_w,alpha,water," + kFEMHeader + ",settings";
    const size_t nCols = 16 + kFEMCols + 1;
    ResumableCSV csv(out, header, {0, 1, 2, 3, 10, 11, 12, 13, 14, nCols - 1});
    const ProgressLog log{LogFileOf(out)};
    std::ostringstream intro;
    intro << "fembatch fig=8: " << pan.size() << " panels x " << hws.size() << " h_w / H x " << curves.size() << " curves, soil "
          << soil << ", gamma_w " << gw << ", alpha " << alpha << "; " << settings << "; out " << out << " (" << csv.NDone()
          << " rows already there)";
    log(intro.str());
    int run = 0, skipped = 0, failed = 0;
    for (const auto &pb : pan)
        for (double hw : hws) {
            REAL c, phi, gamma;
            Fig8Soil(pb.first, soil, c, phi, gamma);
            std::vector<std::pair<std::vector<std::string>, std::string>> todo; // (curves of the run, water)
            if (hw > 0.)
                for (const std::string &cv : curves) todo.push_back({{cv}, cv == "vopt" ? "analytical" : "fe"});
            else todo.push_back({curves, "none"});
            for (auto &job : todo) {
                std::vector<std::vector<std::string>> rows;
                bool all = true;
                for (const std::string &cv : job.first) {
                    std::vector<std::string> row = {pb.first, CsvNum(pb.second), cv, CsvNum(hw), "", "", "", "", "", "",
                                                    CsvNum(c), CsvNum(phi), CsvNum(gamma), CsvNum(gw), CsvNum(alpha), job.second};
                    row.resize(nCols - 1, "");
                    row.push_back(settings);
                    all = all && csv.Done(row);
                    rows.push_back(row);
                }
                if (all) {
                    skipped++;
                    continue;
                }
                std::ostringstream id;
                id << pb.first << " beta=" << pb.second << " hw/H=" << hw << " " << (hw > 0. ? job.first[0] : "none");
                // limit analysis at H = 1 m (the paper's scale, box 50 / 10 / 30 m), then the FEM at H_ref = H_crit
                Problem p1 = FigureProblem(base, 1., pb.second, hw, la.boxm);
                p1.soil.c = c, p1.soil.gamma = gamma, p1.soil.phi = phi * M_PI / 180.;
                log("  [fem fig8] " + id.str() + ": limit analysis, then " + std::to_string(fem.nref + 1) + " FEM cycles");
                const FigureCase lac = RunFigureCase(p1, phi, job.second, la);
                REAL Href = 1.;
                if (lac.r.found && std::isfinite(lac.r.Hcrit) && lac.r.Hcrit > 0.) {
                    const REAL e = std::pow(10., std::floor(std::log10(lac.r.Hcrit)) - 2.);
                    Href = std::round(lac.r.Hcrit / e) * e;
                }
                std::vector<double> boxH = la.boxm;
                for (double &b : boxH) b *= Href; // la.boxm metres at H = 1 m -> the same box in units of H at H_ref
                Problem p = FigureProblem(base, Href, pb.second, hw, boxH);
                p.soil = p1.soil;
                FEMStabilityBox(p, fem);
                const FEMCase fc = RunFEMCase(p, job.second, fem, la.href, la.lm);
                if (fc.nout > 0 || fc.cyc.empty()) {
                    log("  [fem fig8] " + id.str() + ": " + std::to_string(fc.nout) +
                        " integration points of the stability mesh outside the hydraulic mesh, not run");
                    failed++;
                    continue;
                }
                const REAL Hc = fc.Estimate() * Href;
                const std::vector<std::string> mc = FEMColumns(fc, p, lac.lf.seconds + lac.r.seconds);
                for (auto &row : rows) {
                    if (csv.Done(row)) continue;
                    row[4] = CsvNum(Hc), row[5] = CsvNum(fc.Lambda() * Href), row[6] = CsvNum(lac.r.Hcrit);
                    row[7] = CsvNum(Hc / lac.r.Hcrit - 1., 6), row[8] = CsvNum(Href), row[9] = CsvNum(fc.Lambda());
                    std::copy(mc.begin(), mc.end(), row.begin() + 16);
                    csv.Append(row);
                }
                run++;
                std::ostringstream s;
                s << std::setprecision(6) << "  [fem fig8 " << run << "] " << id.str() << ": H_crit FEM = " << Hc
                  << " m (h -> 0; H_ref " << Href << " m, last cycle lambda " << fc.Lambda() << ", " << fc.cyc.back().neq
                  << " equations), LA = " << lac.r.Hcrit << " m (" << std::setprecision(4) << fc.tField + fc.tFEM << " s)";
                log(s.str());
            }
        }
    std::ostringstream end;
    end << "fembatch fig=8: " << run << " cases run, " << skipped << " already in " << out << ", " << failed << " failed, "
        << std::setprecision(4) << std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count() << " s";
    log(end.str());
    return failed == 0 ? 0 : 1;
}

/// fembatch fig=9|8 [FEM options nref= mark= sa= sh0= shs= sgrade= shmax= maxnewton= tolfs=] [limit-analysis
/// options of fig8 / fig9, href= also for the FE field of the FEM] [case options of FEMBatchFig9 / FEMBatchFig8] out=<csv>
int CmdFEMBatch(Options &o) {
    const std::string fig = o.Str("fig", "");
    for (const char *k : {"sleft", "sright", "sdepth", "sref", "smesh", "beta", "alpha"})
        if (o.Has(k)) {
            std::cerr << "fembatch: " << k << " is not an option here (stability box: sa=, capped at the hydraulic box; cases: "
                      << "betas=, alphas= or panels=, hws=)\n";
            return 1;
        }
    const FEMBatchSettings fem = ReadFEMBatchSettings(o);
    const FigureSettings la = ReadFigureSettings(o);
    if (fig == "9") return FEMBatchFig9(o, fem, la, ProjectDir() + "results/cpp/fem_fig9.csv");
    if (fig == "8") return FEMBatchFig8(o, fem, la, ProjectDir() + "results/cpp/fem_fig8.csv");
    std::cerr << "fembatch: fig=9 or fig=8\n";
    return 1;
}

/// self-tests of fembatch (part of 'check'): extrapolations exact on model sequences; stability box inside the
/// hydraulic box for every production case; similarity behind H_ref of fig=8 (geometrically similar problems: equal
/// numbers of equations and load vectors F(H2, lambda H1 / H2) = (H2 / H1) F(H1, lambda) for the FE field, the
/// analytical field and f = 0, and equal lambda H of tiny nonlinear runs); tiny fig=9 / fig=8 runs repeated (the
/// second must skip every case)
void CheckFEMBatch(Checker &ck, const Problem &p0) {
    std::cout << "FEM batch (fembatch)\n";
    const auto t0 = std::chrono::steady_clock::now();
    {
        struct C { int64_t neq; REAL gi; };
        REAL err = 0.;
        for (REAL p : {1., 1.5}) { // lambda_k = 0.9 + 0.4 2^-pk
            std::vector<C> c;
            for (int k = 0; k < 4; k++) c.push_back({int64_t(1000) << (2 * k), 0.9 + 0.4 * std::pow(2., -p * k)});
            const FSExtrapolation e = ExtrapolateCycles(c);
            err = std::max({err, std::fabs(e.richardson - 0.9), std::fabs(e.order - p)});
            if (p == 1.) err = std::max({err, std::fabs(e.h1 - 0.9), std::fabs(e.sqrtneq - 0.9)}); // h ~ neq^-1/2
        }
        std::vector<C> c = {{100, 1.}, {400, 0.9}, {1600, 0.85}, {6400, 0.84}}; // order 2.3: NaN only if not decreasing
        const FSExtrapolation e = ExtrapolateCycles(c);
        std::vector<C> c2 = {{100, 1.}, {400, 0.75}, {1600, 0.5}};                // equal differences: no order
        const FSExtrapolation e2 = ExtrapolateCycles(c2);
        err = std::max(err, REAL(!std::isfinite(e.richardson) + std::isfinite(e2.richardson) + !std::isfinite(e2.h1)));
        ck.Expect("    extrapolations (order 1 and 1.5 sequences, h ~ neq^-1/2) exact", err, 1.e-12);
    }
    Problem base = p0;
    base.hmesh = base.smesh = "gen";
    const FEMBatchSettings def;
    {
        // over the production cases: max of (stability extent - hydraulic extent) / H (<= 0: inside), and the error of
        // each side against the rule sa H + H / tan(beta) capped at the hydraulic box (Fig. 9: right side 10 m from T)
        REAL outside = -1.e300, rule = 0.;
        auto side = [&](const Problem &p) {
            const REAL H = p.hyd.H, ext = (def.sa + 1. / std::tan(p.hyd.beta * M_PI / 180.)) * H;
            const REAL st[3] = {p.stab.left, p.stab.right, p.stab.depth}, hy[3] = {p.hyd.left, p.hyd.right, p.hyd.depth};
            for (int i = 0; i < 3; i++) {
                outside = std::max(outside, (st[i] - hy[i]) / H);
                rule = std::max(rule, std::fabs(st[i] - std::min(ext, hy[i])) / H);
            }
        };
        for (REAL beta : {30., 45., 60., 75., 90.}) {
            Problem p = FigureProblem(base, 5., beta, 1., {50., 10., 30.});
            FEMStabilityBox(p, def);
            side(p);
            rule = std::max(rule, std::fabs(p.stab.right - 10.) / 5.);
        }
        for (REAL beta : {30., 35., 60.})
            for (REAL hw : {0.2, 0.5, 1.}) {
                Problem p = FigureProblem(base, 37., beta, hw, {50. * 37., 10. * 37., 30. * 37.});
                FEMStabilityBox(p, def);
                side(p);
            }
        ck.Expect("    stability box inside the hydraulic box and as sa H + H / tan(beta) capped (Fig. 9: right side at 10 m "
                  "from T; Fig. 8)", std::max<REAL>(outside, 0.) + rule, 1.e-12);
    }
    { // similarity: at H2 = 4 H1 every coordinate scales exactly (power of 2), so the meshes are the same up to the scale
        FEMBatchSettings fs;
        fs.ssize.h0 = fs.ssize.hs = 0.5, fs.ssize.hmax = 2.;
        REAL errF = 0.;
        int64_t bad = 0;
        for (const std::string water : {"fe", "analytical", "none"}) {
            TPZFMatrix<STATE> F[2];
            const REAL Hs[2] = {1.5, 6.}, lambda1 = 2.;
            for (int i = 0; i < 2; i++) {
                const REAL H = Hs[i];
                Problem p = FigureProblem(base, H, 60., water == "none" ? 0. : 0.5, {50. * H, 10. * H, 30. * H});
                p.hsize.hmax = 4.;
                FEMStabilityBox(p, fs);
                ForceField f = NoSeepage();
                if (water == "fe") {
                    TPZGeoMesh *hm = HydraulicMesh(p, 0);
                    const SeepageResult r = DrawdownSeepage(hm, p.hyd, p.hy);
                    delete hm;
                    f = r.field->AsForceField(r.field);
                } else if (water == "analytical") f = AnalyticalForceField(MakeAnalyticalField(p, 10.));
                TPZGeoMesh *sm = StabilityMesh(p);
                F[i] = LoadVector(sm, ModelVoigt(p.soil), p.soil, f, p.soil.gamma - p.hy.gammaw, lambda1 * Hs[0] / H);
                delete sm;
            }
            if (F[0].Rows() != F[1].Rows() || F[0].Rows() == 0) {
                bad++;
                continue;
            }
            TPZFMatrix<STATE> s(F[0]), d(F[1]);
            s *= Hs[1] / Hs[0];
            d -= s;
            errF = std::max(errF, Norm(d) / Norm(s));
        }
        ck.Expect("    similarity: F(4 H, lambda / 4) = 4 F(H, lambda), FE / analytical / no field (+ mesh mismatches)",
                  errF + REAL(bad), 1.e-10);
        // the nonlinear problem: lambda H of tiny runs (FE field, h_w = H / 2) at H and 4 H equal within the
        // continuation tolerance (the trial sequences differ: the first step is 0.5 in lambda)
        fs.nref = 0;
        fs.ds.tolFS = 0.02, fs.ds.maxNewton = 20;
        REAL lh[2];
        const REAL Hs[2] = {2., 8.};
        for (int i = 0; i < 2; i++) {
            Problem p = FigureProblem(base, Hs[i], 60., 0.5, {50. * Hs[i], 10. * Hs[i], 30. * Hs[i]});
            p.hsize.hmax = 4.;
            p.soil.c = 6., p.soil.phi = 32. * M_PI / 180., p.soil.gamma = 18.;
            FEMStabilityBox(p, fs);
            std::streambuf *old = std::cout.rdbuf();
            std::ostringstream sink;
            std::cout.rdbuf(sink.rdbuf());
            const FEMCase fc = RunFEMCase(p, "fe", fs, 0, 10.);
            std::cout.rdbuf(old);
            lh[i] = fc.Lambda() * Hs[i];
        }
        std::cout << "    lambda H at H = 2 and 8 m: " << std::setprecision(8) << lh[0] << ", " << lh[1] << std::setprecision(6) << "\n";
        ck.Expect("    similarity of the elastoplastic collapse (coarse, tolfs 0.02): |lambda H (8 m) / lambda H (2 m) - 1|",
                  std::fabs(lh[1] / lh[0] - 1.), 2. * fs.ds.tolFS);
    }
    { // tiny fembatch runs, each repeated: the second must skip every case and leave the csv unchanged
        const std::filesystem::path dir = std::filesystem::temp_directory_path() /
            ("SlopeSeepageForces_femcheck_" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
        std::filesystem::create_directories(dir);
        auto runCmd = [](std::vector<std::string> args) {
            std::vector<char *> argv = {nullptr};
            for (std::string &a : args) argv.push_back(&a[0]);
            Options o(int(argv.size()), argv.data(), 1);
            std::streambuf *old = std::cout.rdbuf();
            std::ostringstream sink;
            std::cout.rdbuf(sink.rdbuf());
            const int rc = CmdFEMBatch(o);
            std::cout.rdbuf(old);
            return rc;
        };
        auto lines = [](const std::string &f) {
            std::ifstream in(f);
            std::vector<std::string> v;
            for (std::string l; std::getline(in, l);) v.push_back(l);
            return v;
        };
        const std::string f9 = (dir / "fem_fig9.csv").string(), f8 = (dir / "fem_fig8.csv").string();
        const std::vector<std::string> coarse = {"sh0=0.5", "shs=0.5", "shmax=2", "tolfs=0.05", "maxnewton=20", "seeds=0",
                                                 "np=12", "niter=15", "threads=1"};
        std::vector<std::string> a9 = {"fig=9", "alphas=1", "betas=90", "curves=vopt", "nref=1", "out=" + f9};
        std::vector<std::string> a8 = {"fig=8", "panels=Israeli60", "hws=0", "curves=vopt,FE", "nref=0", "out=" + f8};
        a9.insert(a9.end(), coarse.begin(), coarse.end());
        a8.insert(a8.end(), coarse.begin(), coarse.end());
        int rc = runCmd(a9) + runCmd(a8);
        const auto l9 = lines(f9), l8 = lines(f8);
        rc += runCmd(a9) + runCmd(a8);
        const bool same = lines(f9) == l9 && lines(f8) == l8;
        // fig=9: header + 1 row (2 cycles); fig=8: header + the h_w = 0 rows of both curves from one run
        bool shape = l9.size() == 2 && l8.size() == 3, rows = false;
        if (shape) {
            const auto r9 = ResumableCSV::Split(l9[1]), r1 = ResumableCSV::Split(l8[1]), r2 = ResumableCSV::Split(l8[2]);
            rows = r1.size() == r2.size() && r1[2] == "vopt" && r2[2] == "FE" && r1[4] == r2[4] && r1[15] == "none" &&
                   r1[16] == "1" && r1[4] == r1[5] && r9[14] == "2" && std::isfinite(atof(r1[4].c_str())) &&
                   std::isfinite(atof(r9[3].c_str())) && std::isfinite(atof(r9[5].c_str())) && r9[13] == "analytical" &&
                   r9[3] != r9[4];
        }
        ck.Expect("    fembatch fig=9 / fig=8 tiny runs repeated: second run skips all, h_w = 0 rows shared (mismatches)",
                  REAL(rc != 0) + REAL(!same) + REAL(!shape) + REAL(!rows), 0.);
        std::error_code ec;
        std::filesystem::remove_all(dir, ec);
    }
    std::cout << "    (fembatch checks: " << std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count() << " s)\n";
}

/// self-tests of the figure drivers: resumable CSV (cut row, missing newline, keys), hydraulic box fixed in metres,
/// Fig. 8 soil sets, and a tiny fig9 / fig8 run repeated (the second run must skip every case)
void CheckFigures(Checker &ck, const Problem &p0) {
    std::cout << "figure drivers (fig5 out=, fig8, fig9)\n";
    const std::filesystem::path dir = std::filesystem::temp_directory_path() /
                                      ("SlopeSeepageForces_check_" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    std::filesystem::create_directories(dir);
    {
        const std::string file = (dir / "r.csv").string();
        {
            ResumableCSV a(file, "k,v,s", {0, 2});
            a.Append({"1", "10", "x"});
            a.Append({"2", "20", "x"});
        }
        std::ofstream(file, std::ios::app) << "3,3"; // a row cut by a kill (no newline, too few columns)
        ResumableCSV b(file, "k,v,s", {0, 2});
        const bool keys = b.Done({"1", "", "x"}) && b.Done({"2", "", "x"}) && !b.Done({"1", "", "y"}) && !b.Done({"3", "", "x"});
        b.Append({"3", "30", "x"});
        ResumableCSV c(file, "k,v,s", {0, 2});
        std::ifstream in(file);
        std::string line;
        int nl = 0;
        while (std::getline(in, line)) nl++;
        ck.Expect("    resumable csv: keys, cut row ignored, next row on a new line (mismatches)",
                  REAL(!keys) + REAL(b.Ignored() != 1) + REAL(c.Ignored() != 1) + REAL(c.NDone() != 3) + REAL(nl != 5), 0.);
    }
    {
        Problem base = p0;
        base.hmesh = base.smesh = "gen";
        REAL err = 0.;
        for (REAL H : {1., 5., 10.}) {
            const Problem p = FigureProblem(base, H, 60., 0.5, {50., 10., 30.});
            err = std::max({err, std::fabs(p.hyd.left - 50.), std::fabs(p.hyd.right - 10.), std::fabs(p.hyd.depth - 30.),
                            std::fabs(p.hyd.H - H), std::fabs(p.hyd.hw - 0.5 * H), std::fabs(p.hyd.beta - 60.)});
        }
        ck.Expect("    hydraulic box fixed in metres (50 / 10 / 30 m for H = 1, 5, 10 m)", err, 1.e-12);
        REAL c, phi, g, d = 0.;
        Fig8Soil("London", "swapped", c, phi, g), d += std::fabs(c - 11.7) + std::fabs(phi - 24.7) + std::fabs(g - 18.);
        Fig8Soil("Israeli", "swapped", c, phi, g), d += std::fabs(c - 6.) + std::fabs(phi - 32.);
        Fig8Soil("London", "table1", c, phi, g), d += std::fabs(c - 6.) + std::fabs(phi - 32.);
        Fig8Soil("Israeli", "table1", c, phi, g), d += std::fabs(c - 11.7) + std::fabs(phi - 24.7);
        ck.Expect("    Fig. 8 soil sets (swapped / Table 1)", d, 1.e-12);
    }
    { // tiny runs, each repeated: the second must skip every case and leave the csv unchanged
        auto runCmd = [](int (*cmd)(Options &), std::vector<std::string> args) {
            std::vector<char *> argv = {nullptr};
            for (std::string &a : args) argv.push_back(&a[0]);
            Options o(int(argv.size()), argv.data(), 1);
            std::streambuf *old = std::cout.rdbuf();
            std::ostringstream sink;
            std::cout.rdbuf(sink.rdbuf());
            const int rc = cmd(o);
            std::cout.rdbuf(old);
            return rc;
        };
        auto lines = [](const std::string &f) {
            std::ifstream in(f);
            std::vector<std::string> v;
            for (std::string l; std::getline(in, l);) v.push_back(l);
            return v;
        };
        const std::string f9 = (dir / "fig9.csv").string(), f8 = (dir / "fig8.csv").string();
        const std::vector<std::string> a9 = {"alphas=1", "betas=60,90", "curves=vopt", "seeds=0", "np=12", "niter=15", "threads=1",
                                             "out=" + f9};
        const std::vector<std::string> a8 = {"panels=London60", "hws=0,0.5", "curves=vopt,FE", "seeds=0", "np=12", "niter=15",
                                             "threads=1", "href=0", "hhmax=4", "out=" + f8};
        int rc = runCmd(CmdFig9, a9) + runCmd(CmdFig8, a8);
        const auto l9 = lines(f9), l8 = lines(f8);
        rc += runCmd(CmdFig9, a9) + runCmd(CmdFig8, a8);
        const bool same = lines(f9) == l9 && lines(f8) == l8;
        // fig9: header + 2 rows; fig8: header + h_w = 0 (one run, rows vopt and FE with the same value) + 2 rows
        bool shape = l9.size() == 3 && l8.size() == 5, ends = false;
        if (l8.size() == 5) {
            const auto r1 = ResumableCSV::Split(l8[1]), r2 = ResumableCSV::Split(l8[2]);
            ends = r1.size() == r2.size() && r1[2] == "vopt" && r2[2] == "FE" && r1[4] == r2[4] && r1[12] == "none" &&
                   std::isfinite(atof(r1[4].c_str()));
        }
        ck.Expect("    fig8 / fig9 tiny runs repeated: second run skips all, h_w = 0 rows shared (mismatches)",
                  REAL(rc != 0) + REAL(!same) + REAL(!shape) + REAL(!ends), 0.);
    }
    std::error_code ec;
    std::filesystem::remove_all(dir, ec);
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
    CheckLimitAnalysis(ck, p0);
    {
        const std::string src(__FILE__), data = src.substr(0, src.find_last_of('/') + 1) + "data/";
        CheckAnalytical(ck, data + "analytical_seepage_reference.csv", data + "fig5_vector_fill_polygons.csv");
    }
    CheckFigures(ck, p0);
    CheckFEMBatch(ck, p0);
    std::cout << "\ncheck: " << ck.total - ck.failed << " of " << ck.total << " passed, "
              << std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count() << " s\n";
    return ck.failed == 0 ? 0 : 1;
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        std::cout << "usage: SlopeSeepageForces mesh|verify|seepage|fig5|fig8|fig9|analytical|probe|fs|check|la|labatch|fembatch [key=value ...] (see main.cpp)\n";
        return 1;
    }
    Options o(argc, argv, 2);
    const std::string cmd = argv[1];
    if (cmd == "mesh") return CmdMesh(o);
    if (cmd == "verify") return CmdVerify(o);
    if (cmd == "seepage") return CmdSeepage(o);
    if (cmd == "fig5") return CmdFig5(o);
    if (cmd == "analytical") return CmdAnalytical(o);
    if (cmd == "probe") return CmdProbe(o);
    if (cmd == "fs") return CmdFS(o);
    if (cmd == "check") return CmdCheck(o);
    if (cmd == "la") return CmdLA(o, argc, argv);
    if (cmd == "labatch") return CmdLABatch(o);
    if (cmd == "fig8") return CmdFig8(o);
    if (cmd == "fig9") return CmdFig9(o);
    if (cmd == "fembatch") return CmdFEMBatch(o);
    std::cerr << "unknown command " << cmd << "\n";
    return 1;
}
