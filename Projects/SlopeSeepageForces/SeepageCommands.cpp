// Commands of the hydraulic part (see main.cpp): mesh (slope meshes), seepage (FE solution of the drawdown seepage),
// probe (u and f at given points), fig5 (normalized hydraulic functionals of the paper's Fig. 5, as a table or into a
// resumable CSV) and analytical (the field K^-1 v'_opt of AnalyticalSeepage.h).
#include "Commands.h"
#include "ResumableCSV.h"

#include "TPZVTKGeoMesh.h"

#include <random>

namespace slope {

namespace {

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
              << " deg, max triangles " << maxTri << ", " << SecondsSince(t0) << " s\n";
    return worst < 1.e-10 ? 0 : 1;
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

/// fig5 out=<csv>: the Fig. 5 functionals (h_w = H) J(u'_FE) / (k_h H^2 gamma_w^2) (curve J_FE) and
/// -J*(v'_opt) / (k_h H^2 gamma_w^2) (curve minus_Jstar_opt) on the grid of scripts/reproduce_python.py (betas 15 to 90
/// every 7.5 deg, href=1 by default; hydraulic box in units of H, the paper's 50 / 10 / 30 m at H = 1 m), appended to
/// a resumable CSV: the columns of data/paper_fig5.csv, then m and the degenerate flag of the analytical optimum, the
/// number of equations of the FE solution, the time and the settings
int Fig5CSV(Options &o) {
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
                rowA[7] = CsvNum(SecondsSince(ta), 4);
                csv.Append(rowA);
                msg << ": -J*/(kh H^2 gw^2) = " << rowA[3] << " (m = " << rowA[4] << (an->Degenerate() ? ", degenerate" : "") << ")";
                any = true;
            }
            if (fe && !csv.Done(rowF)) {
                const SeepageResult r = SolveDrawdown(p, href);
                rowF[3] = CsvNum(r.Jnorm), rowF[6] = std::to_string(r.neq), rowF[7] = CsvNum(r.seconds, 4);
                csv.Append(rowF);
                msg << "; J_FE/(kh H^2 gw^2) = " << rowF[3] << " (" << r.neq << " equations, " << std::setprecision(3) << r.seconds << " s)";
                any = true;
            }
            if (any) run++, log(msg.str());
            else skipped++;
        }
    std::ostringstream end;
    end << "fig5: " << run << " cases run, " << skipped << " already in " << out << ", " << std::setprecision(4) << SecondsSince(t0)
        << " s";
    log(end.str());
    return 0;
}

} // namespace

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
        const double dt = SecondsSince(t0);
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
        const SeepageResult r = SolveDrawdown(p, ref, ref == refs.back() ? vtk : "");
        std::cout << std::setprecision(10) << ref << "  " << r.nel << "  " << r.neq << "  " << r.J << "  " << r.Jnorm
                  << std::setprecision(4) << "  " << r.pmin << " at (" << r.xpmin[0] << ", " << r.xpmin[1] << ")  " << r.seconds
                  << "\n";
        J.push_back(r.Jnorm);
        if (!csv.empty() && ref == refs.back()) WriteGridCSV(csv, *r.field, p.hyd);
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
    const SeepageResult r = SolveDrawdown(p, href);
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

std::map<std::pair<int, int>, REAL> ReadFig5(const std::string &file, bool solidCurve) {
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

/// J(u_FE) / (k_h H^2 gamma_w^2) versus beta and alpha (paper Fig. 5, dashed curves) on the hydraulic domain, and
/// -J*(v'_opt) / (k_h H^2 gamma_w^2) of the analytical field (AnalyticalSeepage.h, L_m = lm H; solid curves);
/// out=<csv>: production run into a resumable CSV (Fig5CSV)
int CmdFig5(Options &o) {
    if (o.Has("out")) return Fig5CSV(o);
    Problem p = ReadProblem(o);
    const std::vector<double> alphas = o.NumList("alphas", {1., 2., 4., 10.});
    const std::vector<double> betas = o.NumList("betas", {15., 30., 45., 60., 75., 90.});
    const int href = o.IntList("href", {0})[0];
    const std::string file = o.Str("paper", ProjectDir() + "data/fig5_vector_fill_polygons.csv");
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
            if (fe) r = SolveDrawdown(p, href); // solve before printing the row (the solver writes to stdout)
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

/// analytical K^-1 v'_opt: scalars of the optimum, J*, the force at sample points; pts=<file> (x_paper,y_paper lines),
/// csv=<file> (grid), ref=<file> (comparison with scripts/analytical_seepage_reference.py), bench=<n> (timing)
int CmdAnalytical(Options &o) {
    Problem p = ReadProblem(o);
    const REAL lm = o.Num("lm", 10.);
    const std::string pts = o.Str("pts", ""), csv = o.Str("csv", ""), ref = o.Str("ref", "");
    const int64_t bench = int64_t(o.Num("bench", 0.));
    o.CheckUnused();
    if (!ref.empty()) return AnalyticalReferenceTable(ref);
    const SlopeGeometry &g = p.hyd;
    auto t0 = std::chrono::steady_clock::now();
    auto an = MakeAnalyticalField(p, lm);
    const double tb = SecondsSince(t0);
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
        const double t = SecondsSince(t0);
        std::cout << std::setprecision(3) << "evaluation: " << 1.e9 * t / bench << " ns per point (ForceField, " << bench
                  << " points in x in [-3H, x_T + 3H], y in [-4H, 0]; checksum " << s << ")\n";
    }
    return 0;
}

} // namespace slope
