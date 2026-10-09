// Command-line layer of SlopeSeepageForces (the commands and their options are documented in main.cpp): key=value
// options, the slope problem built from them (hydraulic and stability domains and meshes, hydraulics, soil), helpers
// shared by the commands, the self-test counter, and the declarations of the commands (SeepageCommands.cpp,
// LimitAnalysisCommands.cpp, FEMCommands.cpp, SelfTests.cpp).
// Coordinates: NeoPZ (y up, crest edge O at the origin, see SlopeGeometry.h); files and tables in "paper coordinates"
// use x_paper = x, y_paper = -y.
#ifndef SLOPESEEPAGE_COMMANDS_H
#define SLOPESEEPAGE_COMMANDS_H

#include "AnalyticalSeepage.h"
#include "SeepageFE.h"

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <initializer_list>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace slope {

/// key=value options; unknown keys and values that are not numbers (numeric options) are an error
class Options {
    std::map<std::string, std::string> fKV;
    std::set<std::string> fUsed;

    static double ToNumber(const std::string &k, const std::string &s) {
        char *end = nullptr;
        const double v = std::strtod(s.c_str(), &end);
        if (end == s.c_str() || *end != '\0') {
            std::cerr << "option " << k << ": '" << s << "' is not a number\n";
            exit(1);
        }
        return v;
    }

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
        return Has(k) ? ToNumber(k, fKV.at(k)) : def;
    }
    std::vector<double> NumList(const std::string &k, const std::vector<double> &def) {
        fUsed.insert(k);
        if (!Has(k)) return def;
        std::vector<double> v;
        std::stringstream ss(fKV.at(k));
        std::string item;
        while (std::getline(ss, item, ',')) v.push_back(ToNumber(k, item));
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
    /// true (with the message "<cmd>: <key> <why>") if one of the keys is given: options the command sets itself
    bool Forbidden(std::initializer_list<const char *> keys, const std::string &cmd, const std::string &why) const {
        for (const char *k : keys)
            if (Has(k)) {
                std::cerr << cmd << ": " << k << " " << why << "\n";
                return true;
            }
        return false;
    }
};

/// The slope problem of the options common to all commands (main.cpp): hydraulic domain (large box of the seepage
/// solution) and stability domain, their mesh sizes, the hydraulics and the soil
struct Problem {
    SlopeGeometry hyd, stab; ///< hydraulic and stability domains
    MeshSize hsize, ssize;
    std::string hmesh, smesh; ///< gen (CreateSlopeGMesh) or trig (TriGMesh of SlopeMohrCoulomb)
    Hydraulics hy;
    Soil soil;
    int sref = 0;
    REAL hext[3] = {50., 10., 30.}; ///< hydraulic box / H: left of O, right of T, below T
    REAL sa = 2.;                   ///< stability box: sa H + H / tan(beta) on each side ...
    REAL sext[3] = {-1., -1., -1.}; ///< ... unless overridden (/ H, > 0): left of O, right of T, below T

    /// both domains for the slope angle beta (deg) and the water level h_w = hwr H (H = hyd.H)
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
    /// hydraulic box fixed in metres (left of O, right of T, below T), whatever H
    void SetHydraulicBoxMetres(const std::vector<double> &boxm) {
        for (int i = 0; i < 3; i++) hext[i] = boxm[i] / hyd.H;
        SetSlope(hyd.beta, hyd.hw / hyd.H);
    }
};

inline Problem ReadProblem(Options &o) {
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

/// hboxm=<left>,<right>,<depth> (fs, la): hydraulic box in metres, overrides hleft, hright, hdepth; false (with a
/// message) if malformed
inline bool ReadHydraulicBoxMetres(Options &o, Problem &p, const std::string &cmd) {
    const std::vector<double> boxm = o.NumList("hboxm", {});
    if (boxm.empty()) return true;
    if (boxm.size() != 3) {
        std::cerr << cmd << ": hboxm=<left>,<right>,<depth> (m)\n";
        return false;
    }
    p.SetHydraulicBoxMetres(boxm);
    return true;
}

inline TPZGeoMesh *HydraulicMesh(const Problem &p, int ref, MeshStats *st = nullptr) {
    return p.hmesh == "trig" ? SlopeMohrCoulombGMesh(1 + ref) : CreateSlopeGMesh(p.hyd, p.hsize, ref, st);
}

inline TPZGeoMesh *StabilityMesh(const Problem &p, MeshStats *st = nullptr) {
    return p.smesh == "trig" ? SlopeMohrCoulombGMesh(1 + p.sref) : CreateSlopeGMesh(p.stab, p.ssize, p.sref, st);
}

/// drawdown seepage (SeepageFE.h) on the hydraulic mesh of p with href uniform refinements; vtk: optional output
inline SeepageResult SolveDrawdown(const Problem &p, int href, const std::string &vtk = "") {
    TPZGeoMesh *gmesh = HydraulicMesh(p, href);
    SeepageResult r = DrawdownSeepage(gmesh, p.hyd, p.hy, vtk);
    delete gmesh;
    return r;
}

/// the analytical field K^-1 v'_opt (AnalyticalSeepage.h) of the slope and water level of p.hyd (alpha = k_h / k_v),
/// L_m = lm H
inline std::shared_ptr<const AnalyticalSeepage> MakeAnalyticalField(const Problem &p, REAL lm) {
    return std::make_shared<const AnalyticalSeepage>(p.hyd.beta, p.hyd.H, p.hyd.hw, p.hy.kh / p.hy.kv, p.hy.kh, p.hy.gammaw, lm);
}

inline void PrintGeometry(const char *name, const SlopeGeometry &g) {
    std::cout << name << ": H " << g.H << " m, beta " << g.beta << " deg, h_w " << g.hw << " m, box x in [" << -g.left
              << ", " << g.XT() + g.right << "], y in [" << -g.H - g.depth << ", 0]\n";
}

inline void PrintHydraulics(const Hydraulics &hy, const std::string &mesh) {
    std::cout << "K = diag(" << hy.kh << ", " << hy.kv << "), alpha " << hy.kh / hy.kv << ", gamma_w " << hy.gammaw
              << ", order " << hy.order << ", mesh " << mesh << "; far sides: left " << FarSides::Name(hy.far.left)
              << ", base " << FarSides::Name(hy.far.bottom) << ", right " << FarSides::Name(hy.far.right) << "\n";
}

inline int64_t CountLeafTriangles(TPZGeoMesh *gmesh) {
    int64_t n = 0;
    for (int64_t i = 0; i < gmesh->NElements(); i++) {
        TPZGeoEl *gel = gmesh->Element(i);
        if (gel && !gel->HasSubElement() && gel->Dimension() == 2) n++;
    }
    return n;
}

/// directory of the project sources (data/ and results/ are looked up there)
inline std::string ProjectDir() {
    const std::string src(__FILE__);
    return src.substr(0, src.find_last_of('/') + 1);
}

inline double SecondsSince(std::chrono::steady_clock::time_point t0) {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
}

// ------------------------------------------------------------------------------------------------------------------
// self-tests ('check')
// ------------------------------------------------------------------------------------------------------------------

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

/// runs a command with the options args (key=value), its standard output discarded; returns its exit code
inline int RunQuiet(int (*cmd)(Options &), std::vector<std::string> args) {
    std::vector<char *> argv = {nullptr};
    for (std::string &a : args) argv.push_back(&a[0]);
    Options o(int(argv.size()), argv.data(), 1);
    std::streambuf *old = std::cout.rdbuf();
    std::ostringstream sink;
    std::cout.rdbuf(sink.rdbuf());
    const int rc = cmd(o);
    std::cout.rdbuf(old);
    return rc;
}

inline std::vector<std::string> ReadLines(const std::string &file) {
    std::ifstream in(file);
    std::vector<std::string> v;
    for (std::string l; std::getline(in, l);) v.push_back(l);
    return v;
}

/// a new empty directory in the system temporary directory (the caller removes it)
inline std::filesystem::path ScratchDir(const std::string &prefix) {
    const std::filesystem::path dir = std::filesystem::temp_directory_path() /
        (prefix + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    std::filesystem::create_directories(dir);
    return dir;
}

// ------------------------------------------------------------------------------------------------------------------
// commands (dispatched by main.cpp) and the functions shared between their translation units
// ------------------------------------------------------------------------------------------------------------------

// SeepageCommands.cpp: meshes, seepage solution, Fig. 5 functionals, analytical field
int CmdMesh(Options &o);
int CmdSeepage(Options &o);
int CmdProbe(Options &o);
int CmdFig5(Options &o);
int CmdAnalytical(Options &o);
/// Fig. 5 of the paper (data/fig5_vector_fill_polygons.csv): (alpha, beta) -> J(u'_FE) / (k_h H^2 gamma_w^2), the
/// dashed curves (solidCurve: -J*(v'_opt) / (k_h H^2 gamma_w^2), the solid curves = lower edges of the bands)
std::map<std::pair<int, int>, REAL> ReadFig5(const std::string &file, bool solidCurve = false);

// LimitAnalysisCommands.cpp: kinematic limit analysis and the production of Figs. 8 and 9
int CmdLA(Options &o, int argc, char *argv[]);
int CmdLABatch(Options &o);
int CmdFig8(Options &o);
int CmdFig9(Options &o);
void CheckFigures(Checker &ck, const Problem &p0); ///< self-tests of fig5 out=, fig8, fig9

// FEMCommands.cpp: FEM gravity increase (FEMStability.h)
int CmdFS(Options &o);
int CmdFEMBatch(Options &o);
void CheckFEMBatch(Checker &ck, const Problem &p0); ///< self-tests of fembatch

// SelfTests.cpp
int CmdCheck(Options &o);
int CmdVerify(Options &o);
/// analytical ref=<file>: C++ field against the Python reference of scripts/analytical_seepage_reference.py
int AnalyticalReferenceTable(const std::string &file);

} // namespace slope

#endif
