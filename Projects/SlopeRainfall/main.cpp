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
// Usage: SlopeRainfall [solo=f50] [alpha=45] [hs=10] [hw=5] [I=36] [kx=k_s] [ky=k_s] [theta=0] [train=24] [ttotal=48] [tfs=3,6,...] [passivo=0] [acoplado=1]
//                      [gi=1] [fs=1] [vtk=0] [mesh=<file .msh>] [out=<folder>] [verbose=1]
//                      [threads=0 (assembly threads)]
// Outputs in <out>: fs.csv, hist.csv, snaps_pw.csv, snaps_pa.csv, vertices.csv, triangles.csv, outline.csv,
//                   summary.txt (and vtk/ if vtk=1). Figures and animations: scripts/*.py and mma/*.wl.
#include "UnsatSlope.h"
#include "TPZVTKGeoMesh.h"

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
static const REAL SMAX = 75.;                  // limit of the initial suction (kPa)
static const REAL BETA = 7. * M_PI / 180.;     // inclination of the initial water table (Fig. 2)

struct Case {
    std::string solo = "f50";
    REAL alpha = 45., hs = 10., hw = 5., I = 36., train = 24., ttotal = 48.;
    REAL kx = 0., ky = 0., theta = 0.; ///< principal conductivities (m/s; <= 0: k_s of the soil), angle of k_x (deg)
    std::vector<REAL> tfs = {3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 36, 48};
    bool passivo = false, gi = true, fs = true, vtk = false, acoplado = true;
    std::string mesh, out = "out";
    int verbose = 1, threads = 0;
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

class RainfallSlope {
public:
    RainfallSlope(const Case &c) : fCase(c), fGeo(c.hs, c.alpha) {
        if (c.solo == "f10") fPar.hydro.fA = 10., fPar.ks = 1.e-4;
        else if (c.solo == "f50") fPar.hydro.fA = 50., fPar.ks = 1.e-5;
        else if (c.solo == "f100") fPar.hydro.fA = 100., fPar.ks = 1.e-6;
        else DebugStop();
        fPar.kx = c.kx;
        fPar.ky = c.ky;
        fPar.theta = c.theta;
        std::string mesh = c.mesh;
        if (mesh.empty()) {
            std::ostringstream s;
            s << SLOPERAINFALL_SOURCE_DIR << "/meshes/talude_H" << c.hs << "_a" << c.alpha << ".msh";
            mesh = s.str();
        }
        fGMesh = ReadMsh(mesh);
        fMPhys = UnsatMesh(fGMesh, fPar);
        fAn = new UnsatAnalysis(fMPhys, c.threads);
        fAn->SetVerbose(c.verbose > 2 ? 2 : 0);
        fMat = fAn->Material();
        fMat->SetHydroMechanicalCoupling(c.acoplado);
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
    void SetBoundary() {
        std::set<int64_t> surf, left, right;
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
        }
        fAtm.assign(surf.begin(), surf.end());
        fFixed = fAn->DisplacementConstraints();
        const REAL gw = fPar.gw;
        const REAL hl = fGeo.WaterTable(fCase.hw, 0.), hr = fGeo.WaterTable(fCase.hw, fGeo.L);
        for (auto n : left) {
            const REAL y = fGMesh->NodeVec()[n].Coord(1);
            if (y <= hl + 1.e-6) fFixed[fAn->WaterEq().at(n)] = gw * (hl - y);
        }
        for (auto n : right) {
            const REAL y = fGMesh->NodeVec()[n].Coord(1);
            if (y <= hr + 1.e-6) fFixed[fAn->WaterEq().at(n)] = gw * (hr - y);
        }
        if (fCase.passivo)
            for (auto &e : fAn->AirEq()) fFixed[e.second] = 0.;
        else
            for (auto n : surf) fFixed[fAn->AirEq().at(n)] = 0.;
    }

    /// Initial pore pressures: hydrostatic below the water table, suction limited to SMAX, p_a = 0
    void InitialPressures(TPZFMatrix<STATE> &x) {
        for (auto &e : fAn->WaterEq()) {
            const TPZGeoNode &nd = fGMesh->NodeVec()[e.first];
            x(e.second, 0) = std::max(fPar.gw * (fGeo.WaterTable(fCase.hw, nd.Coord(0)) - nd.Coord(1)), -SMAX);
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
        fMat->PlasticModel().SetStrengthReductionFactor(F);
        const auto fx = FrozenFixed(fAn->X());
        for (REAL lam : {0.25, 0.5, 0.75, 1.0}) {
            fMat->SetGravityFactor(lam);
            if (!fAn->Solve(0., fx, {}, true).conv) return false;
            fAn->Accept();
        }
        return true;
    }

    /// Resets the memory to the default (stress free) state and the displacements to zero
    void ResetState(const TPZFMatrix<STATE> &x) {
        auto mem = fMat->GetMemory();
        const TPZUnsaturatedMem def = fMat->GetDefaultMemory();
        for (int64_t i = 0; i < mem->NElements(); ++i) (*mem)[i] = def;
        TPZFMatrix<STATE> &X = fAn->X();
        X = x;
        for (int64_t i = 0; i < fAn->NEq(); ++i)
            if (!fAn->IsPressure(i)) X(i, 0) = 0.;
        fAn->Restore();
    }

    /// Factors of safety (SRM, GI) with the pressures of x frozen; the transient state is restored at the end
    std::pair<REAL, REAL> FactorOfSafety(const TPZFMatrix<STATE> &x, REAL prevSRM, REAL prevGI) {
        const TPZAdmChunkVector<TPZUnsaturatedMem> saved(*fMat->GetMemory());
        const TPZFMatrix<STATE> xsaved(fAn->X());
        const REAL Fsaved = fStrength;
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
            fMat->SetGravityFactor(1.);
            std::vector<REAL> xs = prevSRM < 0 ? std::vector<REAL>{0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.5, 3.0, 3.5, 4.0}
                                               : around(prevSRM, F0);
            if (xs.size() < 2) xs = {F0 + 0.1, F0 + 0.2};
            if (v) std::cout << "   SRM\n";
            fsSRM = CollapseSearch([&](REAL F) {
                fMat->PlasticModel().SetStrengthReductionFactor(F);
                const bool ok = fAn->Solve(0., fx, {}, true).conv;
                if (ok) fAn->Accept();
                else fAn->Restore();
                return ok;
            }, F0, xs, 0.01, 20., v);
        }
        // gravity increase from the stress-free state, real strength
        if (fCase.gi) {
            ResetState(x);
            fMat->PlasticModel().SetStrengthReductionFactor(1.);
            fMat->SetGravityScalesSuction(false);
            const auto fx = FrozenFixed(fAn->X());
            // increments of at most 0.5 up to the neighbourhood of the previous factor
            std::vector<REAL> xs;
            const REAL near = prevGI < 0 ? 1.e30 : prevGI - 0.25;
            for (REAL g = 0.5; g < near - 1.e-9 && (prevGI > 0 || g <= 4.); g += 0.5) xs.push_back(g);
            for (REAL g : (prevGI < 0 ? std::vector<REAL>{4.5, 5.0, 6.0, 7.0, 8.0, 10.0} : around(prevGI, 0.)))
                if (xs.empty() || g > xs.back()) xs.push_back(g);
            if (v) std::cout << "   GI\n";
            fsGI = CollapseSearch([&](REAL lam) {
                fMat->SetGravityFactor(lam);
                const bool ok = fAn->Solve(0., fx, {}, true).conv;
                if (ok) fAn->Accept();
                else fAn->Restore();
                return ok;
            }, 0., xs, 0.01, 20., v);
            fMat->SetGravityScalesSuction(true);
        }
        // transient state back
        *fMat->GetMemory() = saved;
        fAn->X() = xsaved;
        fAn->Restore();
        fMat->SetGravityFactor(1.);
        fMat->PlasticModel().SetStrengthReductionFactor(Fsaved);
        return {fsSRM, fsGI};
    }

    /// One time step with the unilateral rainfall condition (active set of ponded surface vertices)
    bool UnilateralStep(REAL dt, bool raining) {
        const REAL rate = raining ? fCase.I / 1000. / HOUR : 0.;
        std::set<int64_t> ponded = fPonded;
        bool conv = false;
        for (int pass = 0; pass < 12; ++pass) {
            std::map<int64_t, REAL> fixed(fFixed), fq;
            for (auto n : ponded) fixed[fAn->WaterEq().at(n)] = 0.;
            for (auto n : fAtm)
                if (!ponded.count(n) && rate > 0.) fq[fAn->WaterEq().at(n)] = rate * fSupply[n];
            const auto r = fAn->Solve(dt, fixed, fq, false);
            fIter += r.iter;
            conv = r.conv;
            if (!conv) break;
            const TPZFMatrix<STATE> &x = fMPhys->Solution();
            std::set<int64_t> add, drop;
            for (auto n : fAtm)
                if (!ponded.count(n) && x.GetVal(fAn->WaterEq().at(n), 0) > 1.e-6) add.insert(n);
            fRw = fAn->FullResidual();
            for (auto n : ponded)
                if (fRw.GetVal(fAn->WaterEq().at(n), 0) / dt > rate * fSupply[n] * (1. + 1.e-9) + 1.e-14) drop.insert(n);
            if (add.empty() && drop.empty()) break;
            for (auto n : add) ponded.insert(n);
            for (auto n : drop) ponded.erase(n);
        }
        if (!conv) {
            fAn->Restore();
            return false;
        }
        fAn->Accept();
        fPonded = ponded;
        return true;
    }

    /// Step with bisection of the time step (avancar of caso.py)
    bool Advance(REAL dt, bool raining, int depth, int maxdepth) {
        if (UnilateralStep(dt, raining)) return true;
        if (depth >= maxdepth) return false;
        if (fCase.verbose > 1) std::cout << "    bisection: dt = " << dt / 2. / HOUR << " h\n";
        return Advance(dt / 2., raining, depth + 1, maxdepth) && Advance(dt / 2., raining, depth + 1, maxdepth);
    }

    void WriteMeshFiles();
    void Snapshot(REAL th);
    void WriteVTK(int k);

    Case fCase;
    Geometry fGeo;
    UnsatParams fPar;
    TPZGeoMesh *fGMesh = nullptr;
    TPZMultiphysicsCompMesh *fMPhys = nullptr;
    UnsatAnalysis *fAn = nullptr;
    TUnsat *fMat = nullptr;
    std::map<int64_t, REAL> fFixed;
    std::map<int64_t, REAL> fSupply;
    std::vector<int64_t> fAtm;
    std::set<int64_t> fPonded;
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
        if (!gel || gel->MaterialId() != 1) continue;
        t << local[gel->NodeIndex(0)] << "," << local[gel->NodeIndex(1)] << "," << local[gel->NodeIndex(2)] << "\n";
    }
    const REAL h = fGeo.hs;
    o << "x,y\n" << 0. << "," << 0. << "\n" << fGeo.L << "," << 0. << "\n" << fGeo.L << "," << 2 * h << "\n"
      << fGeo.xc << "," << 2 * h << "\n" << fGeo.xb << "," << 3 * h << "\n" << 0. << "," << 3 * h << "\n";
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
    if (!Geostatic(1.)) {
        std::cout << "geostatic state did not converge\n";
        return;
    }
    fAn->ResetDisplacements();
    fMat->SetGravityFactor(1.);
    fAn->Accept();
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

    // time steps: 15 min during the rain, 1 h in the next 24 h, 4 h afterwards
    std::vector<std::pair<REAL, bool>> steps;
    for (REAL t = 0.; t < fCase.ttotal - 1.e-9;) {
        const bool rain = t < fCase.train - 1.e-9;
        REAL dt = rain ? 0.25 : (t < fCase.train + 24. - 1.e-9 ? 1. : 4.);
        dt = std::min(dt, fCase.ttotal - t);
        if (rain) dt = std::min(dt, fCase.train - t);
        steps.push_back({dt, rain});
        t += dt;
    }
    REAL t = 0.;
    REAL tFail = -1.;
    int k = 1;
    for (auto &st : steps) {
        const REAL dt = st.first * HOUR;
        int64_t it0 = fIter;
        if (!Advance(dt, st.second, 0, tFail >= 0. ? 5 : 2)) {
            if (tFail >= 0.) {
                std::cout << fCase.out << ": time step did not converge at t = " << t / HOUR << " h\n";
                break;
            }
            tFail = t / HOUR;
            fStrength = F0;
            fMat->PlasticModel().SetStrengthReductionFactor(F0);
            std::cout << fCase.out << ": transient without equilibrium after t = " << tFail << " h (FS < 1)\n";
            if (!Advance(dt, st.second, 0, 5)) {
                std::cout << fCase.out << ": time step did not converge at t = " << t / HOUR << " h\n";
                break;
            }
        }
        t += dt;
        const REAL th = t / HOUR;
        // infiltration: rain on the flux vertices plus the inflow at the ponded ones
        REAL rain = 0., infil = 0.;
        const REAL rate = st.second ? fCase.I / 1000. / HOUR : 0.;
        for (auto n : fAtm) {
            rain += rate * fSupply[n];
            if (!fPonded.count(n)) infil += rate * fSupply[n];
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
         << fCase.I << "\nkx " << (fCase.kx > 0 ? fCase.kx : fPar.ks) << "\nky "
        << (fCase.ky > 0 ? fCase.ky : fPar.ks) << "\ntheta " << fCase.theta << "\ntrain " << fCase.train << "\nttotal " << fCase.ttotal << "\npassivo " << fCase.passivo << "\nacoplado " << fCase.acoplado
        << "\nelements " << fGMesh->NElements() << "\nequations " << fAn->NEq() << "\nfs_ini_srm "
        << (fsHist.empty() ? -1. : fsHist[0].second.first) << "\nfs_ini_gi "
        << (fsHist.empty() ? -1. : fsHist[0].second.second) << "\nfs_min_srm " << minS << "\nt_min_srm " << tminS
        << "\nfs_min_gi " << minG << "\nt_min_gi " << tminG << "\nt_fail " << tFail << "\nnewton_iterations " << fIter
        << "\nminutes " << elapsed() / 60. << "\n";
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
        else {
            std::cout << "unknown argument " << a << "\n";
            return 1;
        }
    }
    RainfallSlope slope(c);
    slope.Run();
    return 0;
}
