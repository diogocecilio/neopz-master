// Limit-analysis cases, shared by la, labatch, fig8, fig9 and fembatch: the seepage input of the limit analysis
// (FE field, analytical field, none, dry), its optimiser settings, the problem of a case of the paper's figures
// (hydraulic box fixed in metres), the Fig. 8 panels and soils, and the CSV columns of the optimal mechanism.
#ifndef SLOPESEEPAGE_FIGURECASES_H
#define SLOPESEEPAGE_FIGURECASES_H

#include "Commands.h"
#include "LimitAnalysis.h"
#include "ResumableCSV.h"

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <string>
#include <utility>
#include <vector>

namespace slope {

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

inline LAField MakeLAField(const Problem &p, const std::string &water, int href, REAL lm = 10.) {
    LAField lf;
    if (water == "fe") {
        const SeepageResult r = SolveDrawdown(p, href);
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
        lf.seconds = SecondsSince(t0);
        lf.info = "analytical field K^-1 v'_opt: " + an->Summary();
    } else if (water == "none" || water == "dry") {
        lf.info = water == "none" ? "no seepage (h_w = 0): f = 0, buoyant gamma' = gamma - gamma_w" : "dry: f = 0, gamma_w = 0";
    } else {
        std::cerr << "la: water must be fe, analytical, none or dry\n";
        exit(1);
    }
    return lf;
}

/// optimiser options of la, fig8, fig9 and fembatch: mech=I,II seeds=0,1,2 np=40 niter=150 qsearch=coarse qfinal=fine
/// polish=1 pools=25 threads=<defThreads> dmax=10; false (with a message) if invalid
inline bool ReadLASettings(Options &o, const std::string &cmd, int defThreads, la::Settings &st, REAL &dmax) {
    const std::string mech = o.Str("mech", "I,II");
    st.kinds.clear();
    if (mech == "I" || mech == "I,II" || mech == "II,I") st.kinds.push_back(1);
    if (mech == "II" || mech == "I,II" || mech == "II,I") st.kinds.push_back(2);
    st.seeds.clear();
    for (int s : o.IntList("seeds", {0, 1, 2})) st.seeds.push_back(uint64_t(s));
    st.nParticles = int(o.Num("np", 40));
    st.nIter = int(o.Num("niter", 150));
    st.quadSearch = o.Str("qsearch", "coarse");
    st.quadFinal = o.Str("qfinal", "fine");
    st.polish = o.Num("polish", 1) != 0.;
    st.extraPools = int(o.Num("pools", 25));
    st.nThreads = int(o.Num("threads", defThreads));
    dmax = o.Num("dmax", 10.);
    if (st.kinds.empty() || !la::Quadrature::Level(st.quadSearch).Valid() || !la::Quadrature::Level(st.quadFinal).Valid()) {
        std::cerr << cmd << ": mech=I|II|I,II, qsearch / qfinal = coarse, medium, fine, xfine, ref or dense\n";
        return false;
    }
    return true;
}

/// limit-analysis settings of the figure drivers (options and defaults of 'la', but threads=2 and hydraulic mesh
/// href = 1) and the hydraulic box in metres
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

/// the figure settings; the hydraulic box is hboxm= (hleft, hright, hdepth would be overridden: an error)
inline FigureSettings ReadFigureSettings(Options &o, const std::string &cmd) {
    if (o.Forbidden({"hleft", "hright", "hdepth"}, cmd, "is not an option here: the hydraulic box is hboxm=<left>,<right>,<depth> (m)"))
        exit(1);
    FigureSettings fs;
    fs.boxm = o.NumList("hboxm", fs.boxm);
    fs.href = o.IntList("href", {1})[0];
    fs.lm = o.Num("lm", 10.);
    const bool ok = ReadLASettings(o, cmd, 2, fs.st, fs.dmax);
    if (fs.boxm.size() != 3) std::cerr << cmd << ": hboxm=<left>,<right>,<depth> (m)\n";
    if (!ok || fs.boxm.size() != 3) exit(1);
    return fs;
}

/// the problem of one figure case: slope height H (m), angle beta (deg), water level h_w = hwr H, hydraulic box fixed
/// in metres (boxm: left of O, right of T, below T) whatever H (the box of the paper's Fig. 4, see README)
inline Problem FigureProblem(const Problem &base, REAL H, REAL beta, REAL hwr, const std::vector<double> &boxm) {
    Problem p = base;
    p.hyd.H = H;
    for (int i = 0; i < 3; i++) p.hext[i] = boxm[i] / H;
    p.SetSlope(beta, hwr);
    return p;
}

/// the curves of a figure: vopt (K^-1 v'_opt, water=analytical) and FE (-grad u'_FE, water=fe)
inline std::vector<std::string> FigureCurves(Options &o) {
    const std::vector<std::string> cv = SplitList(o.Str("curves", "vopt,FE"));
    for (const std::string &c : cv)
        if (c != "vopt" && c != "FE") {
            std::cerr << "curves must be vopt, FE or vopt,FE\n";
            exit(1);
        }
    if (cv.empty()) exit(1);
    return cv;
}

/// seepage input of a curve at the water level h_w / H = hwr (none: h_w = 0, f = 0 with the buoyant gamma')
inline std::string CurveWater(const std::string &curve, REAL hwr) {
    return hwr > 0. ? (curve == "vopt" ? "analytical" : "fe") : "none";
}

/// one limit analysis of a figure: water = fe (-grad u'_FE), analytical (K^-1 v'_opt) or none (h_w = 0: f = 0 with
/// the buoyant gamma' = gamma - gamma_w); phi in degrees
struct FigureCase {
    la::Result r;
    LAField lf;
};

inline FigureCase RunFigureCase(const Problem &p, REAL phiDeg, const std::string &water, const FigureSettings &fs) {
    FigureCase fc;
    fc.lf = MakeLAField(p, water, fs.href, fs.lm);
    const la::Problem prob(p.hyd.beta, p.hyd.H, p.soil.c, phiDeg, p.soil.gamma, p.hy.gammaw, fc.lf.seep, "auto", fs.dmax);
    fc.r = la::StabilityFactor(prob, fs.st);
    return fc;
}

/// mechanism and run columns of the figure CSVs (after the paper's columns; A, B, C in paper coordinates, metres)
inline const std::string kFigureMechHeader = "mechanism,theta1,theta2,eta_or_dH,L_over_H,Ax,Bx,By,Cx,Cy,r0,P_mr,P_gamma,P_u,"
                                             "seed_spread,at_bound,n_eval,Jnorm,t_field_s,t_la_s";
constexpr size_t kFigureMechCols = 20;

inline std::vector<std::string> FigureMechColumns(const FigureCase &fc, REAL H) {
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
inline std::string FigureCaseSummary(const FigureCase &fc) {
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

// ------------------------------------------------------------------------------------------------------------------
// Fig. 8: London clay panel (beta = 30, 60 deg) and Israeli clay panel (beta = 35, 60 deg), alpha = 1, gamma = 18
// ------------------------------------------------------------------------------------------------------------------

/// soil parameters of a Fig. 8 panel ("London" or "Israeli"): Table 1 (soil=table1) or with the (c, phi) pairs
/// exchanged between the panels (soil=swapped, the set that reproduces the h_w = 0 ends of Fig. 8); gamma = 18
inline void Fig8Soil(const std::string &panel, const std::string &set, REAL &c, REAL &phiDeg, REAL &gamma) {
    const bool london = panel == "London";
    const bool table = set == "table1";
    gamma = 18.;
    if (london == table) c = 6., phiDeg = 32.; // London clay of Table 1
    else c = 11.7, phiDeg = 24.7;              // Israeli clay of Table 1
}

/// panels=London30,London60,Israeli35,Israeli60 (default): (soil, beta in degrees); false (with a message) if invalid
inline bool ReadFig8Panels(Options &o, const std::string &cmd, std::vector<std::pair<std::string, REAL>> &pan) {
    pan.clear();
    for (const std::string &s : SplitList(o.Str("panels", "London30,London60,Israeli35,Israeli60"))) {
        const std::string name = s.rfind("London", 0) == 0 ? "London" : (s.rfind("Israeli", 0) == 0 ? "Israeli" : "");
        char *end = nullptr;
        const REAL beta = name.empty() ? 0. : std::strtod(s.c_str() + name.size(), &end);
        if (name.empty() || end == s.c_str() + name.size() || *end != '\0') {
            std::cerr << cmd << ": panels are London<beta> or Israeli<beta>, e.g. London30,Israeli35\n";
            return false;
        }
        pan.push_back({name, beta});
    }
    return true;
}

/// the runs of a point h_w / H = hwr of Fig. 8: (curves, water) of each run; at h_w = 0 one run without seepage
/// (f = 0 with the buoyant gamma') gives the common end of both curves
inline std::vector<std::pair<std::vector<std::string>, std::string>> Fig8Runs(REAL hwr, const std::vector<std::string> &curves) {
    std::vector<std::pair<std::vector<std::string>, std::string>> runs;
    if (hwr > 0.)
        for (const std::string &cv : curves) runs.push_back({{cv}, CurveWater(cv, hwr)});
    else runs.push_back({curves, "none"});
    return runs;
}

} // namespace slope

#endif
