// Homogeneous unsaturated soil slope under rainfall: reproduction of Rahardjo, Ong, Rezaur & Leong (2007), "Factors
// controlling instability of homogeneous soil slopes under rainfall", J. Geotech. Geoenviron. Eng. 133(12), with the
// coupled u - p_w - p_a model TPZMatUnsaturatedUPwPa (Mohr-Coulomb in the Bishop stress) - port of caso.py (talude3f).
//
// One run = one case of the article: soil (f10, f50, f100), slope angle alpha, height Hs, depth Hw of the water table
// below the toe, rainfall intensity I_r (mm/h) during t_rain hours, total time t_total. The factor of safety at the
// instants t_fs is computed as in the article (an independent stability analysis with the p_w and p_a of the instant
// frozen): strength reduction (SRM, from the geostatic state with the strength increased, F0 = 0.6, so that FS < 1 is
// found) and gravity increase (GI, from the stress-free state: the weight of the soil and the positive pore pressures
// are multiplied by lambda, the suction is kept).
// If the coupled transient loses equilibrium with the real strength (FS < 1), it continues with the strength
// increased (F0) only to follow the pressures, and the time of the rupture is recorded.
//
// Anisotropic conductivity: principal values kx (along the direction at theta degrees, counterclockwise from the x
// axis) and ky (normal to it).
//
// Usage: SlopeRainfall [solo=f50] [alpha=45] [hs=10] [hw=5] [I=36] [kx=k_s] [ky=k_s] [theta=0] [kramin=1e-4] [psi=phi (deg; SlopeRainfallPV)] [smax=75] [train=24] [ttotal=48] [tfs=3,6,...] [passivo=0] [acoplado=1]
//                      [gi=1] [fs=1] [vtk=0] [mesh=<file .msh>] [out=<folder>] [verbose=1]
//                      [threads=0 (assembly threads)]
// Other soils and slopes (the values of solo are the defaults):
//   c= (kPa) phi= (deg) gamma= (kN/m3) E= (kPa) nu= ks= (m/s) poro= (porosity)
//   chi=1|<value>|sw|se       Bishop parameter: constant (on the suction; tan phi^b = chi tan phi'), S_w or S_e
//   swcc=fx:a,n,m | vg:alpha,n,swr   retention curve (Fredlund-Xing, a in kPa; van Genuchten, alpha in 1/kPa)
//   krw=pot:p | mualem        relative permeability of the water
//   camadas=<file>            one line per soil layer (material id of the mesh): "id key=value ..." (same keys)
//   na=x:y,x:y,...            water table (polyline; replaces hw and the inclined water table of the article)
//   lateral=1                 hydrostatic p_w below the water table on the lateral boundaries (0: closed)
//   basena=0                  1: head of the water table prescribed on the base
// Rain:
//   chuva=<file>              hyetograph, rows "t0,t1,I" (h, h, mm/h), replaces I and train
//   padrao=uniforme|adiantado|atrasado   linear hyetographs of mean I during train (Tang et al. 2018)
//   evap=0 (mm/day) psimin=  potential evaporation when it does not rain, limited by the suction psimin (kPa)
//   regime=<q mm/h>           initial state: steady infiltration q (pseudo-transient up to the steady state)
//   alagado=0                 1: surface ponded (p_w = 0) while it rains (I -> infinity)
//   dtmax=<h>                 time steps growing (x 1.25) up to dtmax (default: the steps of the article)
//   dtchuva=0.25              time step during the rain (h)
//   mecanica=1                0 (with acoplado=0): the transient solves only the seepage (displacements fixed; the
//                             FS analyses solve the mechanics from the stress-free state with the pressures frozen)
// Outputs in <out>: fs.csv, hist.csv, snaps_pw.csv, snaps_pa.csv, vertices.csv, triangles.csv, outline.csv,
//                   summary.txt (and vtk/ if vtk=1). Figures and animations: scripts/*.py and mma/*.wl.
#include "UnsatSlope.h"
#include "TPZVTKGeoMesh.h"

#include <array>
#include <chrono>
#include <cstring>
#include <iomanip>
#include <limits>
#include <sys/stat.h>

#ifndef SLOPERAINFALL_SOURCE_DIR
#define SLOPERAINFALL_SOURCE_DIR "."
#endif

static const REAL HOUR = 3600.;
static const REAL F0 = 0.6;                    // strength increase of the geostatic state of the FS analysis
static const REAL BETA = 7. * M_PI / 180.;     // inclination of the initial water table (Fig. 2)

struct Case {
    std::string solo = "f50";
    REAL alpha = 45., hs = 10., hw = 5., I = 36., train = 24., ttotal = 48.;
    REAL kx = 0., ky = 0., theta = 0.; ///< principal conductivities (m/s; <= 0: k_s of the soil), angle of k_x (deg)
    REAL kramin = 1.e-4;                ///< minimum relative permeability of the air
    REAL psi = -1.;                     ///< dilatancy angle (deg); < 0: psi = phi
    REAL smax = 75.;                    ///< limit of the initial suction (kPa)
    std::vector<REAL> tfs = {3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 36, 48};
    bool passivo = false, gi = true, fs = true, vtk = false, acoplado = true;
    std::string mesh, out = "out";
    int verbose = 1, threads = 0;
    /// soil options (key, value) applied over the preset of solo, in order; layers file
    std::vector<std::pair<std::string, std::string>> soil;
    std::string camadas;
    std::vector<std::pair<REAL, REAL>> na; ///< water table polyline
    bool lateral = true, basena = false;
    std::string chuva, padrao = "uniforme";
    REAL evap = 0., psimin = 1.e30, regime = -1., dtmax = -1., dtchuva = 0.25;
    bool alagado = false;
    bool mecanica = true; ///< false (uncoupled only): displacements fixed in the transient, only the seepage is solved
};

/// Geometry of Fig. 2: crest at y = 3 Hs from x = 0 to xb, toe at y = 2 Hs from xc to L
struct Geometry {
    REAL hs, xb, xc, L;
    Geometry(REAL h, REAL alpha) : hs(h), xb(3 * h), xc(3 * h + h / std::tan(alpha * M_PI / 180.)), L(xc + 3 * h) {}
    REAL WaterTable(REAL hw, REAL x) const { return (2 * hs - hw) - std::tan(BETA) * (x - xc); }
};

static std::vector<REAL> ParseList(const char *s) {
    std::vector<REAL> v;
    std::stringstream ss(s);
    std::string tok;
    while (std::getline(ss, tok, ',')) v.push_back(std::atof(tok.c_str()));
    return v;
}

/// Soil option key=value (also the keys of the lines of the layers file)
static bool SoilOption(UnsatParams &p, const std::string &k, const std::string &v) {
    const REAL d = std::atof(v.c_str());
    if (k == "c") p.c = d;
    else if (k == "phi") {
        const bool assoc = std::fabs(p.psi - p.phi) < 1.e-12;
        p.phi = d * M_PI / 180.;
        if (assoc) p.psi = p.phi;
    } else if (k == "psi") p.psi = d * M_PI / 180.;
    else if (k == "gamma") p.gamma = d;
    else if (k == "E") p.E = d;
    else if (k == "nu") p.nu = d;
    else if (k == "ks") p.ks = d;
    else if (k == "poro") p.n = d;
    else if (k == "kx") p.kx = d;
    else if (k == "ky") p.ky = d;
    else if (k == "theta") p.theta = d;
    else if (k == "chi") {
        if (v == "sw") p.bishop = TUnsat::EChiSaturation;
        else if (v == "se") p.bishop = TUnsat::EChiEffective;
        else p.bishop = TUnsat::EChiConstant, p.chi = d;
    } else if (k == "swcc") {
        const std::vector<REAL> a = ParseList(v.c_str() + 3);
        if (v.rfind("fx:", 0) == 0 && a.size() == 3) {
            p.hydro.fRetention = TPZUnsaturatedHydraulics::EFredlundXing;
            p.hydro.fA = a[0], p.hydro.fN = a[1], p.hydro.fM = a[2];
        } else if (v.rfind("vg:", 0) == 0 && a.size() == 3) {
            p.hydro.fRetention = TPZUnsaturatedHydraulics::EVanGenuchten;
            p.hydro.fAVG = a[0], p.hydro.fBVG = a[1], p.hydro.fSwr = a[2];
        } else return false;
    } else if (k == "krw") {
        if (v == "mualem") p.hydro.fKrw = TPZUnsaturatedHydraulics::EMualem;
        else if (v.rfind("pot:", 0) == 0) p.hydro.fKrw = TPZUnsaturatedHydraulics::EPowerLaw, p.hydro.fPower = std::atof(v.c_str() + 4);
        else return false;
    } else return false;
    return true;
}

/// Collapse search of nsat3.collapse_search: follows xs; at the first x without convergence, bisection to the
/// tolerance; the failed x is tried again from the last converged state (increment <= tol)
static REAL CollapseSearch(const std::function<bool(REAL)> &solveAt, REAL x0, const std::vector<REAL> &xs, REAL tol,
                           REAL xmax, bool verbose) {
    REAL xok = x0, xfail = std::numeric_limits<REAL>::infinity();
    while (true) {
        REAL x;
        if (std::isfinite(xfail) && xfail - xok > tol) x = 0.5 * (xok + xfail);
        else if (std::isfinite(xfail)) x = xfail;
        else {
            x = xok + (xs.back() - xs[xs.size() - 2]);
            for (REAL v : xs)
                if (v > xok + 1.e-12) {
                    x = v;
                    break;
                }
        }
        if (x > xmax) return xok;
        const bool conv = solveAt(x);
        if (verbose) std::cout << "     x = " << std::setprecision(5) << x << (conv ? " converged" : " no convergence") << "\n";
        if (conv) {
            if (x == xfail) xfail = std::numeric_limits<REAL>::infinity();
            xok = x;
        } else if (std::isfinite(xfail) && x == xfail) break;
        else xfail = x;
    }
    return xok;
}

/// Rainfall in time: constant I during train (the article), linear patterns, or a hyetograph read from a file
class Rainfall {
public:
    Rainfall(const Case &c) : fI(c.I), fTrain(c.train), fPattern(c.padrao) {
        if (!c.chuva.empty()) {
            std::ifstream in(c.chuva);
            if (!in) {
                std::cout << "rainfall file " << c.chuva << " not found\n";
                DebugStop();
            }
            std::string line;
            while (std::getline(in, line)) {
                if (line.empty() || line[0] == '#') continue;
                const std::vector<REAL> v = ParseList(line.c_str());
                if (v.size() >= 3 && v[1] > v[0]) fBlocks.push_back({v[0], v[1], v[2]});
            }
        } else if (fPattern != "uniforme" && fPattern != "adiantado" && fPattern != "atrasado") {
            std::cout << "unknown rainfall pattern " << fPattern << "\n";
            DebugStop();
        }
    }
    bool Series() const { return !fBlocks.empty(); }
    /// mean intensity (mm/h) in [t0, t1] (h)
    REAL Mean(REAL t0, REAL t1) const {
        if (t1 <= t0) return 0.;
        REAL vol = 0.;
        if (Series()) {
            for (auto &b : fBlocks) vol += b[2] * std::max(REAL(0.), std::min(t1, b[1]) - std::max(t0, b[0]));
        } else {
            const REAL a = std::max(REAL(0.), t0), z = std::min(t1, fTrain);
            if (z > a) {
                // integral of I(t) on [a, z]
                if (fPattern == "uniforme") vol = fI * (z - a);
                else if (fPattern == "atrasado") vol = fI / fTrain * (z * z - a * a);
                else vol = 2. * fI * (z - a) - fI / fTrain * (z * z - a * a);
            }
        }
        return vol / (t1 - t0);
    }
    /// end of the rain period that contains t (or t if it does not rain at t)
    bool Raining(REAL t) const {
        if (!Series()) return t < fTrain - 1.e-9;
        for (auto &b : fBlocks)
            if (t >= b[0] - 1.e-9 && t < b[1] - 1.e-9 && b[2] > 0.) return true;
        return false;
    }
    /// instants where the intensity changes
    std::vector<REAL> Breaks() const {
        std::vector<REAL> b;
        if (Series())
            for (auto &x : fBlocks) b.push_back(x[0]), b.push_back(x[1]);
        else b.push_back(fTrain);
        return b;
    }
    /// end of the last rain before t (-infinity if none)
    REAL LastRainEnd(REAL t) const {
        REAL e = -1.e30;
        if (!Series()) return t >= fTrain - 1.e-9 ? fTrain : e;
        for (auto &b : fBlocks)
            if (b[2] > 0. && b[1] <= t + 1.e-9) e = std::max(e, b[1]);
        return e;
    }

private:
    REAL fI, fTrain;
    std::string fPattern;
    std::vector<std::array<REAL, 3>> fBlocks;
};

class RainfallSlope {
public:
    RainfallSlope(const Case &c) : fCase(c), fGeo(c.hs, c.alpha), fRain(c) {
        if (c.solo == "f10") fPar.hydro.fA = 10., fPar.ks = 1.e-4;
        else if (c.solo == "f50") fPar.hydro.fA = 50., fPar.ks = 1.e-5;
        else if (c.solo == "f100") fPar.hydro.fA = 100., fPar.ks = 1.e-6;
        else DebugStop();
        fPar.kx = c.kx;
        fPar.ky = c.ky;
        fPar.theta = c.theta;
        fPar.hydro.fKraMin = c.kramin;
        if (c.psi >= 0.) fPar.psi = c.psi * M_PI / 180.;
        for (auto &o : c.soil)
            if (!SoilOption(fPar, o.first, o.second)) {
                std::cout << "invalid soil option " << o.first << "=" << o.second << "\n";
                DebugStop();
            }
        if (!c.camadas.empty()) ReadLayers(c.camadas);
        std::string mesh = c.mesh;
        if (mesh.empty()) {
            std::ostringstream s;
            s << SLOPERAINFALL_SOURCE_DIR << "/meshes/talude_H" << c.hs << "_a" << c.alpha << ".msh";
            mesh = s.str();
        }
        fGMesh = ReadMsh(mesh);
        fXmin = 1.e30, fXmax = -1.e30;
        for (int64_t i = 0; i < fGMesh->NNodes(); ++i) {
            fXmin = std::min(fXmin, fGMesh->NodeVec()[i].Coord(0));
            fXmax = std::max(fXmax, fGMesh->NodeVec()[i].Coord(0));
        }
        fMPhys = UnsatMesh(fGMesh, fPar, fLayers);
        fAn = new UnsatAnalysis(fMPhys, c.threads);
        fAn->SetVerbose(c.verbose > 2 ? 2 : 0);
        fAn->SetHydroMechanicalCoupling(c.acoplado);
        SetBoundary();
    }

    ~RainfallSlope() {
        delete fAn;
        TPZManVector<TPZCompMesh *, 3> atomic(fMPhys->MeshVector());
        delete fMPhys;
        for (auto *m : atomic) delete m;
        delete fGMesh;
    }

    void Run();

private:
    /// parameters of a material id (layer)
    const UnsatParams &Params(int id) const {
        auto it = fLayers.find(id);
        return it == fLayers.end() ? fPar : it->second;
    }

    void ReadLayers(const std::string &file) {
        std::ifstream in(file);
        if (!in) {
            std::cout << "layers file " << file << " not found\n";
            DebugStop();
        }
        std::string line;
        while (std::getline(in, line)) {
            if (line.empty() || line[0] == '#') continue;
            std::istringstream ss(line);
            int id;
            if (!(ss >> id)) continue;
            UnsatParams p = fPar;
            std::string kv;
            while (ss >> kv) {
                const auto e = kv.find('=');
                if (e == std::string::npos || !SoilOption(p, kv.substr(0, e), kv.substr(e + 1))) {
                    std::cout << "invalid layer option " << kv << "\n";
                    DebugStop();
                }
            }
            fLayers[id] = p;
        }
    }

    /// elevation of the water table at x
    REAL WaterTable(REAL x) const {
        const auto &na = fCase.na;
        if (na.empty()) return fGeo.WaterTable(fCase.hw, x);
        if (x <= na.front().first) return na.front().second;
        for (size_t i = 1; i < na.size(); ++i)
            if (x <= na[i].first) {
                const REAL s = (x - na[i - 1].first) / (na[i].first - na[i - 1].first);
                return na[i - 1].second + s * (na[i].second - na[i - 1].second);
            }
        return na.back().second;
    }

    void SetBoundary() {
        std::set<int64_t> surf, left, right, base;
        for (int64_t i = 0; i < fGMesh->NElements(); ++i) {
            TPZGeoEl *gel = fGMesh->Element(i);
            if (!gel || gel->Dimension() != 1) continue;
            const int id = gel->MaterialId();
            const int64_t a = gel->NodeIndex(0), b = gel->NodeIndex(1);
            if (id == ECrest || id == EFace || id == EToe) {
                surf.insert(a);
                surf.insert(b);
                // rain on the horizontal projection, half to each vertex
                const REAL dx = std::fabs(gel->NodePtr(1)->Coord(0) - gel->NodePtr(0)->Coord(0));
                fSupply[a] += dx / 2.;
                fSupply[b] += dx / 2.;
            } else if (id == ELeft) left.insert(a), left.insert(b);
            else if (id == ERight) right.insert(a), right.insert(b);
            else if (id == EBase) base.insert(a), base.insert(b);
        }
        fAtm.assign(surf.begin(), surf.end());
        fFixed = fAn->DisplacementConstraints();
        const REAL gw = fPar.gw;
        auto hydrostatic = [&](const std::set<int64_t> &nodes) {
            for (auto n : nodes) {
                const TPZGeoNode &nd = fGMesh->NodeVec()[n];
                const REAL h = WaterTable(nd.Coord(0)), y = nd.Coord(1);
                if (y <= h + 1.e-6) fFixed[fAn->WaterEq().at(n)] = gw * (h - y);
            }
        };
        if (fCase.lateral) hydrostatic(left), hydrostatic(right);
        if (fCase.basena) hydrostatic(base);
        if (fCase.passivo)
            for (auto &e : fAn->AirEq()) fFixed[e.second] = 0.;
        else
            for (auto n : surf) fFixed[fAn->AirEq().at(n)] = 0.;
    }

    /// Initial pore pressures: hydrostatic below the water table, suction limited to smax, p_a = 0
    void InitialPressures(TPZFMatrix<STATE> &x) {
        for (auto &e : fAn->WaterEq()) {
            const TPZGeoNode &nd = fGMesh->NodeVec()[e.first];
            x(e.second, 0) = std::max(fPar.gw * (WaterTable(nd.Coord(0)) - nd.Coord(1)), -fCase.smax);
        }
        for (auto &e : fAn->AirEq()) x(e.second, 0) = 0.;
    }

    /// Fixed equations of a frozen analysis: displacements and all the pressures at their values in x
    std::map<int64_t, REAL> FrozenFixed(const TPZFMatrix<STATE> &x) {
        std::map<int64_t, REAL> fx = fAn->DisplacementConstraints();
        for (int64_t i = 0; i < fAn->NEq(); ++i)
            if (fAn->IsPressure(i)) fx[i] = x.GetVal(i, 0);
        return fx;
    }

    /// Weight of the soil applied in 4 increments with the pressures (scaled together) frozen, strength factor F
    bool Geostatic(REAL F) {
        fAn->SetStrengthReductionFactor(F);
        const auto fx = FrozenFixed(fAn->X());
        for (REAL lam : {0.25, 0.5, 0.75, 1.0}) {
            fAn->SetGravityFactor(lam);
            if (!fAn->Solve(0., fx, {}, true).conv) return false;
            fAn->Accept();
        }
        return true;
    }

    /// Resets the memory to the default (stress free) state and the displacements to zero
    void ResetState(const TPZFMatrix<STATE> &x) {
        fAn->ResetMemory();
        TPZFMatrix<STATE> &X = fAn->X();
        X = x;
        for (int64_t i = 0; i < fAn->NEq(); ++i)
            if (!fAn->IsPressure(i)) X(i, 0) = 0.;
        fAn->Restore();
    }

    /// Factors of safety (SRM, GI) with the pressures of x frozen; the transient state is restored at the end
    std::pair<REAL, REAL> FactorOfSafety(const TPZFMatrix<STATE> &x, REAL prevSRM, REAL prevGI) {
        const auto saved = fAn->SaveMemory();
        const TPZFMatrix<STATE> xsaved(fAn->X());
        const REAL Fsaved = fStrength;
#ifdef SLOPERAINFALL_PV
        // the dilatancy acts in the coupled transient; the factor of safety is computed with psi = phi (associative),
        // the same method of the closed-form model (the strength reduction with psi < phi loses the convergence of
        // Newton before the collapse)
        std::vector<TMC> transientModels;
        for (auto *m : fAn->Materials()) {
            transientModels.push_back(m->PlasticModel());
            const UnsatParams &p = Params(m->Id());
            TPZElasticResponse ER = m->PlasticModel().fER;
            m->PlasticModel().fYC.SetUp(p.phi, p.phi, p.c, ER);
        }
#endif
        auto around = [](REAL prev, REAL lo) {
            std::vector<REAL> l;
            for (REAL d : {-0.25, -0.08, -0.02, 0.03, 0.1, 0.4})
                if (prev + d > lo) l.push_back(std::round((prev + d) * 1.e4) / 1.e4);
            return l;
        };
        const bool v = fCase.verbose > 1;
        REAL fsSRM = -1., fsGI = -1.;
        // strength reduction from the geostatic state with F0
        ResetState(x);
        if (Geostatic(F0)) {
            const auto fx = FrozenFixed(fAn->X());
            fAn->SetGravityFactor(1.);
            std::vector<REAL> xs = prevSRM < 0 ? std::vector<REAL>{0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.5, 3.0, 3.5, 4.0}
                                               : around(prevSRM, F0);
            if (xs.size() < 2) xs = {F0 + 0.1, F0 + 0.2};
            if (v) std::cout << "   SRM\n";
            fsSRM = CollapseSearch([&](REAL F) {
                fAn->SetStrengthReductionFactor(F);
                const bool ok = fAn->Solve(0., fx, {}, true).conv;
                if (ok) fAn->Accept();
                else fAn->Restore();
                return ok;
            }, F0, xs, 0.01, 20., v);
        }
        // gravity increase from the stress-free state, real strength
        if (fCase.gi) {
            ResetState(x);
            fAn->SetStrengthReductionFactor(1.);
            fAn->SetGravityScalesSuction(false);
            const auto fx = FrozenFixed(fAn->X());
            // increments of at most 0.5 up to the neighbourhood of the previous factor
            std::vector<REAL> xs;
            const REAL near = prevGI < 0 ? 1.e30 : prevGI - 0.25;
            for (REAL g = 0.5; g < near - 1.e-9 && (prevGI > 0 || g <= 4.); g += 0.5) xs.push_back(g);
            for (REAL g : (prevGI < 0 ? std::vector<REAL>{4.5, 5.0, 6.0, 7.0, 8.0, 10.0} : around(prevGI, 0.)))
                if (xs.empty() || g > xs.back()) xs.push_back(g);
            if (v) std::cout << "   GI\n";
            fsGI = CollapseSearch([&](REAL lam) {
                fAn->SetGravityFactor(lam);
                const bool ok = fAn->Solve(0., fx, {}, true).conv;
                if (ok) fAn->Accept();
                else fAn->Restore();
                return ok;
            }, 0., xs, 0.01, 20., v);
            fAn->SetGravityScalesSuction(true);
        }
        // transient state back
        fAn->RestoreMemory(saved);
        fAn->X() = xsaved;
        fAn->Restore();
        fAn->SetGravityFactor(1.);
#ifdef SLOPERAINFALL_PV
        for (size_t i = 0; i < fAn->Materials().size(); ++i) fAn->Materials()[i]->PlasticModel() = transientModels[i];
#endif
        fAn->SetStrengthReductionFactor(Fsaved);
        return {fsSRM, fsGI};
    }

    /**
     * One time step with the unilateral conditions of the surface: rain (rate > 0, m/s on the horizontal projection:
     * active set of ponded vertices, p_w <= 0, inflow <= rain) or evaporation (rate < 0: active set of dry vertices,
     * p_w >= -psimin, outflow <= evaporation). ponded = true: the whole surface ponded (p_w = 0).
     */
    bool UnilateralStep(REAL dt, REAL rate, bool pondedAll = false) {
        std::set<int64_t> ponded = fPonded, dry = rate < 0. ? fDry : std::set<int64_t>{};
        if (pondedAll) ponded.insert(fAtm.begin(), fAtm.end());
        bool conv = false;
        for (int pass = 0; pass < 12; ++pass) {
            std::map<int64_t, REAL> fixed(fFixed), fq;
            for (auto n : ponded) fixed[fAn->WaterEq().at(n)] = 0.;
            for (auto n : dry) fixed[fAn->WaterEq().at(n)] = -fCase.psimin;
            for (auto n : fAtm)
                if (!ponded.count(n) && !dry.count(n) && rate != 0.) fq[fAn->WaterEq().at(n)] = rate * fSupply[n];
            const auto r = fAn->Solve(dt, fixed, fq, false);
            fIter += r.iter;
            conv = r.conv;
            if (!conv && rate > 0. && ponded.size() < fAtm.size()) {
                // rain much larger than the conductivity on a dry soil: the flux condition gives very large pressures
                // before the active set is updated; start again with the whole surface ponded (p_w = 0) and let the
                // active set release the vertices whose inflow would exceed the rain
                fAn->Restore();
                ponded.insert(fAtm.begin(), fAtm.end());
                continue;
            }
            if (!conv) break;
            fRw = fAn->FullResidual();
            if (pondedAll) break;
            const TPZFMatrix<STATE> &x = fMPhys->Solution();
            std::set<int64_t> add, drop, adddry, dropdry;
            for (auto n : fAtm) {
                if (ponded.count(n) || dry.count(n)) continue;
                const REAL pw = x.GetVal(fAn->WaterEq().at(n), 0);
                if (pw > 1.e-6) add.insert(n);
                else if (rate < 0. && pw < -fCase.psimin - 1.e-6) adddry.insert(n);
            }
            const REAL qin = std::max(rate, REAL(0.));
            for (auto n : ponded)
                if (fRw.GetVal(fAn->WaterEq().at(n), 0) / dt > qin * fSupply[n] * (1. + 1.e-9) + 1.e-14) drop.insert(n);
            for (auto n : dry)
                if (fRw.GetVal(fAn->WaterEq().at(n), 0) / dt < rate * fSupply[n] * (1. + 1.e-9) - 1.e-14) dropdry.insert(n);
            if (add.empty() && drop.empty() && adddry.empty() && dropdry.empty()) break;
            for (auto n : add) ponded.insert(n);
            for (auto n : drop) ponded.erase(n);
            for (auto n : adddry) dry.insert(n);
            for (auto n : dropdry) dry.erase(n);
        }
        if (!conv) {
            fAn->Restore();
            return false;
        }
        fAn->Accept();
        fPonded = pondedAll ? std::set<int64_t>(fAtm.begin(), fAtm.end()) : ponded;
        fDry = dry;
        return true;
    }

    /// Step with bisection of the time step (avancar of caso.py)
    bool Advance(REAL dt, REAL rate, bool pondedAll, int depth, int maxdepth) {
        if (UnilateralStep(dt, rate, pondedAll)) return true;
        if (depth >= maxdepth) return false;
        if (fCase.verbose > 1) std::cout << "    bisection: dt = " << dt / 2. / HOUR << " h\n";
        return Advance(dt / 2., rate, pondedAll, depth + 1, maxdepth) && Advance(dt / 2., rate, pondedAll, depth + 1, maxdepth);
    }

    /// surface flux (m/s, positive: rain) of the interval [t0, t1] (h): mean rain, or the evaporation if dry
    REAL SurfaceRate(REAL t0, REAL t1) const {
        const REAL r = fRain.Mean(t0, t1);
        if (r > 0.) return r / 1000. / HOUR;
        return fCase.evap > 0. ? -fCase.evap / 24. / 1000. / HOUR : 0.;
    }

    /// Time steps (h): 15 min during the rain, 1 h in the next 24 h, 4 h afterwards (the article); with dtmax,
    /// steps growing 1.25x up to dtmax; aligned with the changes of the rain and the instants of the FS
    std::vector<std::pair<REAL, bool>> TimeSteps() const {
        std::vector<std::pair<REAL, bool>> steps;
        std::vector<REAL> breaks = fRain.Breaks();
        if (fRain.Series() || fCase.dtmax > 0.)
            for (REAL tf : fCase.tfs) breaks.push_back(tf);
        std::sort(breaks.begin(), breaks.end());
        REAL dtprev = 0.;
        for (REAL t = 0.; t < fCase.ttotal - 1.e-9;) {
            const bool rain = fRain.Raining(t);
            REAL dt;
            if (fCase.dtmax > 0.) dt = dtprev <= 0. ? 0.25 : std::min(1.25 * dtprev, fCase.dtmax);
            else dt = rain ? fCase.dtchuva : (t < fRain.LastRainEnd(t) + 24. - 1.e-9 ? 1. : 4.);
            dt = std::min(dt, fCase.ttotal - t);
            for (REAL b : breaks)
                if (b > t + 1.e-9) {
                    dt = std::min(dt, b - t);
                    break;
                }
            steps.push_back({dt, rain});
            dtprev = dt;
            t += dt;
        }
        return steps;
    }

    /// Steady infiltration q (mm/h) from the current state: pseudo-time steps doubling from 24 h until the largest
    /// change of p_w is below 1e-3 kPa
    bool SteadyInfiltration(REAL q) {
        const REAL rate = q / 1000. / HOUR;
        REAL dt = 24.;
        for (int k = 0; k < 60; ++k) {
            const TPZFMatrix<STATE> x0(fAn->X());
            if (!Advance(dt * HOUR, rate, false, 0, 6)) return false;
            REAL dp = 0.;
            for (auto &e : fAn->WaterEq()) dp = std::max(dp, std::fabs(fAn->X().GetVal(e.second, 0) - x0.GetVal(e.second, 0)));
            if (fCase.verbose) std::cout << fCase.out << ": steady infiltration, dt = " << dt << " h, max |dp_w| = " << dp << " kPa\n";
            if (dp < 1.e-3 && dt > 1.e3) return true;
            dt = std::min(dt * 2., 1.e6);
        }
        return true;
    }

    void WriteMeshFiles();
    void Snapshot(REAL th);
    void WriteVTK(int k);

    Case fCase;
    Geometry fGeo;
    Rainfall fRain;
    UnsatParams fPar;
    std::map<int, UnsatParams> fLayers;
    TPZGeoMesh *fGMesh = nullptr;
    TPZMultiphysicsCompMesh *fMPhys = nullptr;
    UnsatAnalysis *fAn = nullptr;
    REAL fXmin = 0., fXmax = 0.;
    std::map<int64_t, REAL> fFixed;
    std::map<int64_t, REAL> fSupply;
    std::vector<int64_t> fAtm;
    std::set<int64_t> fPonded, fDry;
    TPZFMatrix<STATE> fRw;
    REAL fStrength = 1.;
    int64_t fIter = 0;
    std::ofstream fSnapW, fSnapA;
    std::vector<int64_t> fNodes;
    bool fVTKDefined = false;
};

void RainfallSlope::WriteMeshFiles() {
    std::ofstream v(fCase.out + "/vertices.csv"), t(fCase.out + "/triangles.csv"), o(fCase.out + "/outline.csv");
    std::map<int64_t, int64_t> local;
    for (auto &e : fAn->WaterEq()) {
        local[e.first] = fNodes.size();
        fNodes.push_back(e.first);
    }
    v << "x,y\n" << std::setprecision(10);
    for (auto n : fNodes) v << fGMesh->NodeVec()[n].Coord(0) << "," << fGMesh->NodeVec()[n].Coord(1) << "\n";
    t << "v0,v1,v2\n";
    for (int64_t i = 0; i < fGMesh->NElements(); ++i) {
        TPZGeoEl *gel = fGMesh->Element(i);
        if (!gel || gel->Dimension() != 2) continue;
        t << local[gel->NodeIndex(0)] << "," << local[gel->NodeIndex(1)] << "," << local[gel->NodeIndex(2)] << "\n";
    }
    // outline: the boundary edges chained from the lowest leftmost vertex, corners only
    std::map<int64_t, std::vector<int64_t>> adj;
    for (int64_t i = 0; i < fGMesh->NElements(); ++i) {
        TPZGeoEl *gel = fGMesh->Element(i);
        if (!gel || gel->Dimension() != 1) continue;
        adj[gel->NodeIndex(0)].push_back(gel->NodeIndex(1));
        adj[gel->NodeIndex(1)].push_back(gel->NodeIndex(0));
    }
    auto X = [&](int64_t n) { return fGMesh->NodeVec()[n].Coord(0); };
    auto Y = [&](int64_t n) { return fGMesh->NodeVec()[n].Coord(1); };
    int64_t start = adj.begin()->first;
    for (auto &a : adj)
        if (X(a.first) < X(start) - 1.e-9 || (std::fabs(X(a.first) - X(start)) < 1.e-9 && Y(a.first) < Y(start))) start = a.first;
    std::vector<int64_t> loop = {start};
    int64_t prev = -1, cur = start;
    while (true) {
        int64_t next = -1;
        for (auto n : adj[cur])
            if (n != prev) {
                next = n;
                break;
            }
        if (next < 0 || next == start) break;
        loop.push_back(next);
        prev = cur;
        cur = next;
        if (loop.size() > adj.size()) break;
    }
    o << "x,y\n";
    const size_t m = loop.size();
    for (size_t i = 0; i < m; ++i) {
        const int64_t a = loop[(i + m - 1) % m], b = loop[i], c = loop[(i + 1) % m];
        const REAL cr = (X(b) - X(a)) * (Y(c) - Y(b)) - (Y(b) - Y(a)) * (X(c) - X(b));
        if (std::fabs(cr) > 1.e-8) o << X(b) << "," << Y(b) << "\n";
    }
    fSnapW.open(fCase.out + "/snaps_pw.csv");
    fSnapA.open(fCase.out + "/snaps_pa.csv");
    fSnapW << std::setprecision(7);
    fSnapA << std::setprecision(7);
}

void RainfallSlope::Snapshot(REAL th) {
    const TPZFMatrix<STATE> &x = fAn->X();
    fSnapW << th;
    fSnapA << th;
    for (auto n : fNodes) {
        fSnapW << "," << fAn->Nodal(x, n, 1);
        fSnapA << "," << fAn->Nodal(x, n, 2);
    }
    fSnapW << "\n";
    fSnapA << "\n";
    fSnapW.flush();
    fSnapA.flush();
}

void RainfallSlope::WriteVTK(int k) {
    if (!fCase.vtk) return;
    if (!fVTKDefined) {
        mkdir((fCase.out + "/vtk").c_str(), 0755);
        TPZStack<std::string> scal, vec;
        for (auto s : {"WaterPressure", "AirPressure", "Suction", "Saturation", "BishopPressure"}) scal.Push(s);
        vec.Push("Displacement");
        fAn->DefineGraphMesh(2, scal, vec, fCase.out + "/vtk/rain.vtk");
        fVTKDefined = true;
    }
    fAn->SetStep(k);
    fAn->PostProcess(0, 2);
}

void RainfallSlope::Run() {
    const auto t0 = std::chrono::steady_clock::now();
    auto elapsed = [&]() { return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count(); };
    mkdir(fCase.out.c_str(), 0755);
    WriteMeshFiles();
    std::ofstream fs(fCase.out + "/fs.csv"), hist(fCase.out + "/hist.csv");
    fs << "t_h,fs_srm,fs_gi,ponded,pa_max\n" << std::setprecision(6);
    hist << "t_h,dt_h,raining,rain_m3s,infil_m3s,ponded,pw_max,pa_max,iterations\n" << std::setprecision(8);

    // initial state: hydrostatic pressures with the suction limited, weight of the soil, u = 0
    TPZFMatrix<STATE> &X = fAn->X();
    InitialPressures(X);
    for (auto &f : fFixed) X(f.first, 0) = f.second;
    fAn->Restore();
    // weight with the real strength; if it does not converge (slope at the limit), weight with the strength
    // increased (F0) and then the strength reduced back to the real one; if that fails, the slope is unstable at
    // t = 0 and the transient follows the pressures with F0
    REAL tFail = -1.;
    // uncoupled seepage: the deformation does not change the pressures, the transient mechanics is solved with the
    // strength increased (F = 0.25) only to keep it in equilibrium; the stability is given by the FS analyses
    if (!fCase.acoplado) fStrength = 0.25;
    if (!Geostatic(fStrength)) {
        fAn->Restore();
        ResetState(TPZFMatrix<STATE>(fAn->X()));
        if (!Geostatic(F0)) {
            std::cout << "geostatic state did not converge\n";
            return;
        }
        fAn->SetGravityFactor(1.);
        const auto fx = FrozenFixed(fAn->X());
        bool ok = true;
        for (REAL F : {0.7, 0.8, 0.9, 0.95, 1.0}) {
            fAn->SetStrengthReductionFactor(F);
            ok = fAn->Solve(0., fx, {}, true).conv;
            if (!ok) break;
            fAn->Accept();
        }
        if (!ok) {
            tFail = 0.;
            fStrength = F0;
            fAn->Restore();
            fAn->SetStrengthReductionFactor(F0);
            std::cout << fCase.out << ": no equilibrium with the real strength at t = 0 (FS < 1)\n";
        }
    }
    fAn->ResetDisplacements();
    fAn->SetGravityFactor(1.);
    fAn->Accept();
    // uncoupled seepage without the transient mechanics: all the displacement equations fixed
    if (!fCase.acoplado && !fCase.mecanica)
        for (int64_t i = 0; i < fAn->NEq(); ++i)
            if (!fAn->IsPressure(i)) fFixed[i] = X.GetVal(i, 0);
    // initial steady infiltration (Pradhan et al. 2022: 3.3e-8 m/s for five years)
    if (fCase.regime >= 0.) {
        if (!SteadyInfiltration(fCase.regime)) std::cout << fCase.out << ": steady infiltration did not converge\n";
        fAn->ResetDisplacements();
        fAn->Accept();
    }
    auto pmax = [&](int space) {
        REAL m = -1.e300;
        for (auto &e : (space == 1 ? fAn->WaterEq() : fAn->AirEq())) m = std::max(m, X.GetVal(e.second, 0));
        return m;
    };
    REAL fsS = -1., fsG = -1.;
    std::vector<std::pair<REAL, std::pair<REAL, REAL>>> fsHist;
    auto computeFS = [&](REAL th) {
        if (!fCase.fs) return;
        const auto t1 = elapsed();
        const TPZFMatrix<STATE> xc(X);
        auto f = FactorOfSafety(xc, fsS, fsG);
        fsS = f.first;
        fsG = f.second;
        fsHist.push_back({th, f});
        fs << th << "," << fsS << "," << fsG << "," << fPonded.size() << "," << pmax(2) << "\n";
        fs.flush();
        if (fCase.verbose)
            std::cout << fCase.out << ": t = " << std::setw(6) << th << " h  FS(SRM) = " << std::setprecision(4) << fsS
                      << "  FS(GI) = " << fsG << "  p_a max = " << pmax(2) << " kPa  (" << int(elapsed() - t1) << " s; assembly "
                      << int(fAn->TimeAssemble()) << " s, solver " << int(fAn->TimeSolve())
                      << " s, total " << std::setprecision(3) << elapsed() / 60. << " min)" << std::endl;
    };
    Snapshot(0.);
    WriteVTK(0);
    computeFS(0.);

    const std::vector<std::pair<REAL, bool>> steps = TimeSteps();
    REAL t = 0.;
    int k = 1;
    for (auto &st : steps) {
        const REAL dt = st.first * HOUR;
        const REAL th0 = t / HOUR;
        const REAL rate = SurfaceRate(th0, th0 + st.first);
        const bool pondedAll = fCase.alagado && st.second;
        int64_t it0 = fIter;
        // a failure with 2 bisection levels may be numerical (Newton, active set of the rain): the step is tried
        // again with 5 levels and the real strength before it is taken as the loss of equilibrium of the slope
        if (!Advance(dt, rate, pondedAll, 0, tFail >= 0. ? 5 : 2) && (tFail >= 0. || !Advance(dt, rate, pondedAll, 0, 5))) {
            if (tFail >= 0.) {
                std::cout << fCase.out << ": time step did not converge at t = " << t / HOUR << " h\n";
                break;
            }
            tFail = t / HOUR;
            if (fCase.acoplado) { // uncoupled: the mechanics already uses the increased strength (0.25)
                fStrength = F0;
                fAn->SetStrengthReductionFactor(F0);
            }
            std::cout << fCase.out << ": transient without equilibrium after t = " << tFail << " h"
                      << (fCase.acoplado ? " (FS < 1)" : " (seepage or mechanics without convergence)") << "\n";
            if (!Advance(dt, rate, pondedAll, 0, 5)) {
                std::cout << fCase.out << ": time step did not converge at t = " << t / HOUR << " h\n";
                break;
            }
        }
        t += dt;
        const REAL th = t / HOUR;
        // infiltration: rain on the flux vertices plus the inflow at the ponded ones
        REAL rain = 0., infil = 0.;
        for (auto n : fAtm) {
            rain += std::max(rate, REAL(0.)) * fSupply[n];
            if (!fPonded.count(n) && !fDry.count(n)) infil += rate * fSupply[n];
            else if (fRw.Rows()) infil += fRw.GetVal(fAn->WaterEq().at(n), 0) / dt;
        }
        hist << th << "," << st.first << "," << st.second << "," << rain << "," << infil << "," << fPonded.size()
             << "," << pmax(1) << "," << pmax(2) << "," << fIter - it0 << "\n";
        const bool integer = std::fabs(th - std::round(th)) < 1.e-6;
        if (integer) Snapshot(th);
        for (REAL tf : fCase.tfs)
            if (std::fabs(th - tf) < 1.e-6) {
                WriteVTK(k++);
                computeFS(th);
            }
    }
    // summary
    std::ofstream sum(fCase.out + "/summary.txt");
    REAL minS = 1.e300, minG = 1.e300, tminS = 0., tminG = 0.;
    for (auto &h : fsHist) {
        if (h.second.first < minS) minS = h.second.first, tminS = h.first;
        if (h.second.second < minG) minG = h.second.second, tminG = h.first;
    }
    sum << "solo " << fCase.solo << "\nalpha " << fCase.alpha << "\nhs " << fCase.hs << "\nhw " << fCase.hw << "\nI "
         << fCase.I << "\nkx " << (fPar.kx > 0 ? fPar.kx : fPar.ks) << "\nky "
        << (fPar.ky > 0 ? fPar.ky : fPar.ks) << "\ntheta " << fPar.theta << "\nkramin " << fCase.kramin << "\npsi " << fPar.psi * 180. / M_PI << "\nsmax " << fCase.smax << "\ntrain " << fCase.train << "\nttotal " << fCase.ttotal << "\npassivo " << fCase.passivo << "\nacoplado " << fCase.acoplado
        << "\nelements " << fGMesh->NElements() << "\nequations " << fAn->NEq() << "\nfs_ini_srm "
        << (fsHist.empty() ? -1. : fsHist[0].second.first) << "\nfs_ini_gi "
        << (fsHist.empty() ? -1. : fsHist[0].second.second) << "\nfs_min_srm " << minS << "\nt_min_srm " << tminS
        << "\nfs_min_gi " << minG << "\nt_min_gi " << tminG << "\nt_fail " << tFail << "\nnewton_iterations " << fIter
        << "\nminutes " << elapsed() / 60. << "\n";
    // options other than the defaults
    for (auto &o : fCase.soil) sum << "opt_" << o.first << " " << o.second << "\n";
    if (!fCase.camadas.empty()) sum << "camadas " << fCase.camadas << "\n";
    if (!fCase.chuva.empty()) sum << "chuva " << fCase.chuva << "\n";
    if (fCase.padrao != "uniforme") sum << "padrao " << fCase.padrao << "\n";
    if (fCase.evap > 0.) sum << "evap " << fCase.evap << "\npsimin " << fCase.psimin << "\n";
    if (fCase.regime >= 0.) sum << "regime " << fCase.regime << "\n";
    if (fCase.alagado) sum << "alagado 1\n";
    if (fCase.dtmax > 0.) sum << "dtmax " << fCase.dtmax << "\n";
    if (fCase.dtchuva != 0.25) sum << "dtchuva " << fCase.dtchuva << "\n";
    if (!fCase.mecanica) sum << "mecanica 0\n";
    if (!fCase.na.empty()) sum << "na_pontos " << fCase.na.size() << "\n";
    if (!fCase.lateral) sum << "lateral 0\n";
    if (fCase.basena) sum << "basena 1\n";
    std::cout << "== " << fCase.out << " finished in " << elapsed() / 60. << " min (FS min SRM " << minS << " at "
              << tminS << " h, GI " << minG << ")\n";
}

int main(int argc, char *argv[]) {
    Case c;
    for (int i = 1; i < argc; ++i) {
        const char *a = argv[i];
        auto is = [&](const char *k) { return !std::strncmp(a, k, std::strlen(k)); };
        auto val = [&](const char *k) { return a + std::strlen(k); };
        if (is("solo=")) c.solo = val("solo=");
        else if (is("alpha=")) c.alpha = std::atof(val("alpha="));
        else if (is("hs=")) c.hs = std::atof(val("hs="));
        else if (is("hw=")) c.hw = std::atof(val("hw="));
        else if (is("I=")) c.I = std::atof(val("I="));
        else if (is("kx=")) c.kx = std::atof(val("kx="));
        else if (is("ky=")) c.ky = std::atof(val("ky="));
        else if (is("theta=")) c.theta = std::atof(val("theta="));
        else if (is("kramin=")) c.kramin = std::atof(val("kramin="));
        else if (is("psi=")) c.psi = std::atof(val("psi="));
        else if (is("smax=")) c.smax = std::atof(val("smax="));
        else if (is("train=")) c.train = std::atof(val("train="));
        else if (is("ttotal=")) c.ttotal = std::atof(val("ttotal="));
        else if (is("tfs=")) c.tfs = ParseList(val("tfs="));
        else if (is("passivo=")) c.passivo = std::atoi(val("passivo="));
        else if (is("acoplado=")) c.acoplado = std::atoi(val("acoplado="));
        else if (is("gi=")) c.gi = std::atoi(val("gi="));
        else if (is("fs=")) c.fs = std::atoi(val("fs="));
        else if (is("vtk=")) c.vtk = std::atoi(val("vtk="));
        else if (is("mesh=")) c.mesh = val("mesh=");
        else if (is("out=")) c.out = val("out=");
        else if (is("verbose=")) c.verbose = std::atoi(val("verbose="));
        else if (is("threads=")) c.threads = std::atoi(val("threads="));
        else if (is("camadas=")) c.camadas = val("camadas=");
        else if (is("na=")) {
            std::stringstream ss(val("na="));
            std::string tok;
            while (std::getline(ss, tok, ',')) {
                const auto e = tok.find(':');
                if (e == std::string::npos) {
                    std::cout << "invalid water table point " << tok << "\n";
                    return 1;
                }
                c.na.push_back({std::atof(tok.substr(0, e).c_str()), std::atof(tok.substr(e + 1).c_str())});
            }
        } else if (is("lateral=")) c.lateral = std::atoi(val("lateral="));
        else if (is("basena=")) c.basena = std::atoi(val("basena="));
        else if (is("chuva=")) c.chuva = val("chuva=");
        else if (is("padrao=")) c.padrao = val("padrao=");
        else if (is("evap=")) c.evap = std::atof(val("evap="));
        else if (is("psimin=")) c.psimin = std::atof(val("psimin="));
        else if (is("regime=")) c.regime = std::atof(val("regime="));
        else if (is("alagado=")) c.alagado = std::atoi(val("alagado="));
        else if (is("dtmax=")) c.dtmax = std::atof(val("dtmax="));
        else if (is("dtchuva=")) c.dtchuva = std::atof(val("dtchuva="));
        else if (is("mecanica=")) c.mecanica = std::atoi(val("mecanica="));
        else {
            const char *e = std::strchr(a, '=');
            UnsatParams test;
            if (e && SoilOption(test, std::string(a, e - a), e + 1)) c.soil.push_back({std::string(a, e - a), e + 1});
            else {
                std::cout << "unknown argument " << a << "\n";
                return 1;
            }
        }
    }
    RainfallSlope slope(c);
    slope.Run();
    return 0;
}
