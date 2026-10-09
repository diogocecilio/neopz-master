// Commands of the kinematic limit analysis (LimitAnalysis.h, see main.cpp): la (one case), labatch (a file of cases)
// and the production of the paper's Figs. 8 and 9 (fig8, fig9) into resumable CSVs (default results/cpp/): the
// columns of data/paper_fig*.csv, then the mechanism, run information and a settings string; a progress log next to
// each CSV (.log). Also the self-tests of the figure drivers (CheckFigures, part of 'check').
#include "Commands.h"
#include "FigureCases.h"

#include <thread>

namespace slope {

namespace {

/// one limit analysis with the options of 'la' (see main.cpp); id: case label of the csv row
int RunLA(Options &o, const std::string &id, int defThreads, const std::string &csvOut = "") {
    Problem p = ReadProblem(o);
    const REAL beta = p.hyd.beta, H = p.hyd.H, hwr = p.hyd.hw / H;
    if (!ReadHydraulicBoxMetres(o, p, "la")) return 1;
    const REAL phiDeg = o.Num("phi", 30.), c = p.soil.c, gamma = p.soil.gamma; // phi in degrees, as given
    const std::string water = o.Str("water", hwr > 0. ? "fe" : "none");
    const int href = o.IntList("href", {0})[0];
    const REAL lm = o.Num("lm", 10.); // water=analytical: L_m / H
    const std::string pu = o.Str("pu", "auto"), out = csvOut.empty() ? o.Str("out", "") : csvOut;
    la::Settings st;
    REAL dmax;
    if (!ReadLASettings(o, "la", defThreads, st, dmax)) return 1;
    const std::vector<double> xgiven = o.NumList("x", {}); // theta1,theta2,s: evaluate this mechanism only
    const bool verbose = o.Num("verbose", 0) != 0.;
    o.CheckUnused();
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

int DefaultThreads() { return int(std::min(4u, std::max(1u, std::thread::hardware_concurrency()))); }

} // namespace

int CmdLA(Options &o, int argc, char *argv[]) {
    std::string id;
    for (int i = 2; i < argc; i++) id += (i > 2 ? " " : "") + std::string(argv[i]);
    return RunLA(o, id, DefaultThreads());
}

/// labatch cases=<file> out=<csv> [threads=]: one 'la' case per line of the file (key=value options; '#' comments);
/// cases whose line is already in the first column of out are skipped (resumable: rows are appended one by one)
int CmdLABatch(Options &o) {
    const std::string cases = o.Str("cases", ""), out = o.Str("out", "");
    const int threads = int(o.Num("threads", DefaultThreads()));
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

/// fig8: H_crit versus h_w / H of the paper's Fig. 8 (alpha = 1): London clay panel (beta = 30, 60 deg) and Israeli
/// clay panel (beta = 35, 60 deg), curves vopt (K^-1 v'_opt) and FE (-grad u'_FE, hydraulic box fixed in metres,
/// default 50 / 10 / 30 m with u = 0 on the left side and the base), and their common h_w = 0 end (f = 0 with the
/// buoyant gamma', both rows from one run). H_crit = Gamma(H) H (similarity) with the reference height H = 1 m, at
/// which the box is the paper's 50 / 10 / 30 H. soil=swapped (default): the (c, phi) pairs of Table 1 exchanged
/// between the panels; soil=table1: as printed. gammaw=9.8 by default (it reproduces the h_w = 0 ends, see README).
int CmdFig8(Options &o) {
    if (o.Forbidden({"beta", "hw", "c", "phi", "gamma"}, "fig8", "is set by the panels (panels=, hws=, soil=)")) return 1;
    const auto t0 = std::chrono::steady_clock::now();
    Problem base = ReadProblem(o);
    const REAL H = o.Num("H", 1.);
    base.hy.gammaw = o.Num("gammaw", 9.8);
    const std::string soil = o.Str("soil", "swapped");
    std::vector<std::pair<std::string, REAL>> pan;
    if (!ReadFig8Panels(o, "fig8", pan)) return 1;
    const std::vector<double> hws = o.NumList("hws", {0., 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.});
    const std::vector<std::string> curves = FigureCurves(o);
    const FigureSettings fs = ReadFigureSettings(o, "fig8");
    const std::string out = o.Str("out", ProjectDir() + "results/cpp/" + (soil == "table1" ? "fig8_table1.csv" : "fig8.csv"));
    o.CheckUnused();
    if (soil != "swapped" && soil != "table1") {
        std::cerr << "fig8: soil must be swapped or table1\n";
        return 1;
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
            for (auto &job : Fig8Runs(hw, curves)) {
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
                if (all) {
                    skipped++;
                    continue;
                }
                std::ostringstream id;
                id << pb.first << " beta=" << pb.second << " hw/H=" << hw << " " << (hw > 0. ? job.first[0] : "none");
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
        << ", " << std::setprecision(4) << SecondsSince(t0) << " s";
    log(end.str());
    return 0;
}

/// fig9: stability factor Gamma versus beta of the paper's Fig. 9 (H = 5 m, c = 10 kPa, phi = 30 deg, gamma = 20,
/// h_w = H, gamma_w = 9.81 by default), alphas=1,5,10, curves vopt (K^-1 v'_opt) and FE (-grad u'_FE, hydraulic box
/// fixed in metres, default 50 / 10 / 30 m = 10 / 2 / 6 H, u = 0 on the left side and the base); betas as in
/// scripts/reproduce_python.py (every 5 deg from 15 to 90 plus 37.5, 52.5, 67.5, 82.5)
int CmdFig9(Options &o) {
    if (o.Forbidden({"beta", "alpha"}, "fig9", "is not an option here (betas=, alphas=)")) return 1;
    const auto t0 = std::chrono::steady_clock::now();
    Problem base = ReadProblem(o); // defaults H = 5, c = 10, phi = 30, gamma = 20, gamma_w = 9.81, h_w / H = 1
    const REAL phi = o.Num("phi", 30.), H = base.hyd.H, hwr = base.hyd.hw / base.hyd.H;
    const std::vector<double> alphas = o.NumList("alphas", {1., 5., 10.});
    const std::vector<double> betas = o.NumList("betas", {15., 20., 25., 30., 35., 37.5, 40., 45., 50., 52.5, 55., 60., 65.,
                                                          67.5, 70., 75., 80., 82.5, 85., 90.});
    const std::vector<std::string> curves = FigureCurves(o);
    const FigureSettings fs = ReadFigureSettings(o, "fig9");
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
                const std::string water = CurveWater(cv, hwr);
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
        << ", " << std::setprecision(4) << SecondsSince(t0) << " s";
    log(end.str());
    return 0;
}

/// self-tests of the figure drivers: resumable CSV (cut row, missing newline, keys), hydraulic box fixed in metres,
/// Fig. 8 soil sets, and a tiny fig9 / fig8 run repeated (the second run must skip every case)
void CheckFigures(Checker &ck, const Problem &p0) {
    std::cout << "figure drivers (fig5 out=, fig8, fig9)\n";
    const std::filesystem::path dir = ScratchDir("SlopeSeepageForces_check_");
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
        const int nl = int(ReadLines(file).size());
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
        const std::string f9 = (dir / "fig9.csv").string(), f8 = (dir / "fig8.csv").string();
        const std::vector<std::string> a9 = {"alphas=1", "betas=60,90", "curves=vopt", "seeds=0", "np=12", "niter=15", "threads=1",
                                             "out=" + f9};
        const std::vector<std::string> a8 = {"panels=London60", "hws=0,0.5", "curves=vopt,FE", "seeds=0", "np=12", "niter=15",
                                             "threads=1", "href=0", "hhmax=4", "out=" + f8};
        int rc = RunQuiet(CmdFig9, a9) + RunQuiet(CmdFig8, a8);
        const auto l9 = ReadLines(f9), l8 = ReadLines(f8);
        rc += RunQuiet(CmdFig9, a9) + RunQuiet(CmdFig8, a8);
        const bool same = ReadLines(f9) == l9 && ReadLines(f8) == l8;
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

} // namespace slope
