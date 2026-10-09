// Commands of the FEM stability (FEMStability.h: gravity increase of SlopeAnalysis.h with the seepage forces, refinement
// cycles of the plastic zone), see main.cpp: fs (one analysis, any field and form of the load) and fembatch (the cases
// of the paper's Figs. 9 and 8, an independent check of the limit-analysis curves: each row also holds the limit
// analysis of the same case with the same seepage field; resumable CSVs results/cpp/fem_fig9.csv and fem_fig8.csv,
// progress in the .log next to each CSV). Also the self-tests of fembatch (CheckFEMBatch, part of 'check').
#include "Commands.h"
#include "FEMStability.h"
#include "FigureCases.h"

namespace slope {

int CmdFS(Options &o) {
    Problem p = ReadProblem(o);
    if (!ReadHydraulicBoxMetres(o, p, "fs")) return 1;
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
        const SeepageResult r = SolveDrawdown(p, href);
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
    const double total = SecondsSince(t0);
    std::cout << "\ncycle  equations  lambda_GI  H_crit (m)  gamma H_crit / c" << (srm ? "  FS_SRM" : "")
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
// fembatch
// ------------------------------------------------------------------------------------------------------------------

namespace {

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
    /// settings string of the CSV rows (no commas); the elastic constants of the soil (they change the path, in
    /// theory not the collapse factor) are listed only when they differ from the defaults E = 20000, nu = 0.3, which
    /// keeps the keys of the rows computed before they were listed
    std::string Describe(const Soil &soil) const {
        std::ostringstream s;
        s << std::setprecision(10) << "fem:nref=" << nref << ";mark=" << ds.markFrac << ";smesh=" << ssize.h0 << "/" << ssize.hs
          << "/" << ssize.grade << "/" << ssize.hmax << ";sa=" << sa << ";maxnewton=" << ds.maxNewton << ";tolfs=" << ds.tolFS;
        if (soil.E != 20000. || soil.nu != 0.3) s << ";E=" << soil.E << ";nu=" << soil.nu;
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

/// seepage force of a FEM case: water = fe (-grad u'_FE on the hydraulic box of p, href uniform refinements; the
/// field is returned in fe), analytical (K^-1 v'_opt, L_m = lm H) or none (f = 0)
ForceField FEMForceField(const Problem &p, const std::string &water, int href, REAL lm, std::shared_ptr<const PoreField> &fe) {
    fe.reset();
    if (water == "fe") {
        fe = SolveDrawdown(p, href).field;
        return fe->AsForceField(fe);
    }
    if (water == "analytical") return AnalyticalForceField(MakeAnalyticalField(p, lm));
    if (water != "none") DebugStop();
    return NoSeepage();
}

/// one FEM gravity-increase analysis of a figure case (FEMForceField); b = lambda (gamma' g + f) with gamma' = gamma -
/// gamma_w for the three fields. nout > 0: integration points of the stability mesh outside the hydraulic mesh (not run)
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
    std::shared_ptr<const PoreField> field;
    const ForceField f = FEMForceField(p, water, href, lm, field);
    fc.tField = SecondsSince(t0);
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
    fc.tFEM = SecondsSince(t0);
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
    const std::string settings = fem.Describe(base.soil) + "|la:" + la.Describe(base);
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
                const std::string water = CurveWater(cv, hwr);
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
                row[6] = CsvNum(lac.r.found ? fc.Estimate() / lac.r.Gamma - 1. : NAN, 6);
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
        << std::setprecision(4) << SecondsSince(t0) << " s";
    log(end.str());
    return failed == 0 ? 0 : 1;
}

/// fembatch fig=8: Fig. 8 panels (alpha = 1, gamma = 18, (c, phi) of Table 1 swapped between the panels by default,
/// gamma_w = 9.8) panels=London30,London60,Israeli35,Israeli60 hws=0,0.2,0.5,1 curves=FE,vopt (h_w = 0: one run
/// without seepage, f = 0 with gamma', rows for both curves) -> results/cpp/fem_fig8.csv (H_crit FEM = Estimate()
/// H_ref, the last cycle, the limit analysis H_crit). The limit analysis runs at H = 1 m (the paper's scale, box
/// la.boxm = 50 / 10 / 30 m); the FEM at the reference height H_ref = H_crit of that limit analysis (3 significant
/// digits; hydraulic box 50 / 10 / 30 H_ref, the same problem by similarity), so that lambda ~ 1 (the continuation of
/// SlopeAnalysis starts with steps of 0.5 and stops at 100); H_crit = lambda H_ref (exact similarity, see the check)
int FEMBatchFig8(Options &o, const FEMBatchSettings &fem, const FigureSettings &la, const std::string &outDefault) {
    if (o.Forbidden({"beta", "hw", "c", "phi", "gamma", "H"}, "fembatch fig=8", "is set by the panels (panels=, hws=, soil=)"))
        return 1;
    Problem base = ReadProblem(o);
    base.hy.gammaw = o.Num("gammaw", 9.8);
    const std::string soil = o.Str("soil", "swapped");
    std::vector<std::pair<std::string, REAL>> pan;
    if (!ReadFig8Panels(o, "fembatch fig=8", pan)) return 1;
    const std::vector<double> hws = o.NumList("hws", {0., 0.2, 0.5, 1.});
    const std::vector<std::string> curves = FigureCurves(o);
    const std::string out = o.Str("out", outDefault);
    o.CheckUnused();
    if (soil != "swapped" && soil != "table1") {
        std::cerr << "fembatch fig=8: soil must be swapped or table1\n";
        return 1;
    }
    const auto t0 = std::chrono::steady_clock::now();
    const REAL alpha = base.hy.kh / base.hy.kv, gw = base.hy.gammaw;
    const std::string settings = fem.Describe(base.soil) + "|la:" + la.Describe(base) + ";soil=" + soil + ";Href=LA3";
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
            for (auto &job : Fig8Runs(hw, curves)) {
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
                    row[7] = CsvNum(lac.r.found ? Hc / lac.r.Hcrit - 1. : NAN, 6), row[8] = CsvNum(Href), row[9] = CsvNum(fc.Lambda());
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
        << std::setprecision(4) << SecondsSince(t0) << " s";
    log(end.str());
    return failed == 0 ? 0 : 1;
}

} // namespace

/// fembatch fig=9|8 [FEM options nref= mark= sa= sh0= shs= sgrade= shmax= maxnewton= tolfs=] [limit-analysis
/// options of fig8 / fig9, href= also for the FE field of the FEM] [case options of FEMBatchFig9 / FEMBatchFig8] out=<csv>
int CmdFEMBatch(Options &o) {
    const std::string fig = o.Str("fig", "");
    if (o.Forbidden({"sleft", "sright", "sdepth", "sref", "smesh", "beta", "alpha"}, "fembatch",
                    "is not an option here (stability box: sa=, capped at the hydraulic box; cases: betas=, alphas= or "
                    "panels=, hws=)"))
        return 1;
    const FEMBatchSettings fem = ReadFEMBatchSettings(o);
    const FigureSettings la = ReadFigureSettings(o, "fembatch");
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
                std::shared_ptr<const PoreField> field;
                const ForceField f = FEMForceField(p, water, 0, 10., field);
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
        const std::filesystem::path dir = ScratchDir("SlopeSeepageForces_femcheck_");
        const std::string f9 = (dir / "fem_fig9.csv").string(), f8 = (dir / "fem_fig8.csv").string();
        const std::vector<std::string> coarse = {"sh0=0.5", "shs=0.5", "shmax=2", "tolfs=0.05", "maxnewton=20", "seeds=0",
                                                 "np=12", "niter=15", "threads=1"};
        std::vector<std::string> a9 = {"fig=9", "alphas=1", "betas=90", "curves=vopt", "nref=1", "out=" + f9};
        std::vector<std::string> a8 = {"fig=8", "panels=Israeli60", "hws=0", "curves=vopt,FE", "nref=0", "out=" + f8};
        a9.insert(a9.end(), coarse.begin(), coarse.end());
        a8.insert(a8.end(), coarse.begin(), coarse.end());
        int rc = RunQuiet(CmdFEMBatch, a9) + RunQuiet(CmdFEMBatch, a8);
        const auto l9 = ReadLines(f9), l8 = ReadLines(f8);
        rc += RunQuiet(CmdFEMBatch, a9) + RunQuiet(CmdFEMBatch, a8);
        const bool same = ReadLines(f9) == l9 && ReadLines(f8) == l8;
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
    std::cout << "    (fembatch checks: " << SecondsSince(t0) << " s)\n";
}

} // namespace slope
