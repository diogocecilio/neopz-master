// main.cpp — SlopeSeepageRandom
//
// Reprodução, por elementos finitos elastoplásticos no NeoPZ, de
//   M. Vargas Ceron, D. L. Cecílio, R. V. Linn, S. Maghous, "Stability Analysis of Slope Subjected to Seepage
//   Forces Considering Spatial Variability of Soil Properties", Int J Numer Anal Methods Geomech 49 (2025)
//   2459-2491,
// com Mohr-Coulomb (TPZPlasticStepVoigt<TPZYCMohrCoulombPV2>) e Cam-Clay modificado (TPZModifiedCamClay), campos
// aleatórios de c, φ e kv pela expansão de Karhunen-Loève (TPZMatKLKernel + pzdoublestrmatriz + LAPACK) e
// forças de percolação do problema de Darcy desacoplado (TPZDarcyFlow).
//
// Uso: SlopeSeepageRandom <comando> [opções]   (ver Usage())
//
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <type_traits>
#include <string>
#include <vector>

#include "CoupledDrawdown.h"
#include "KLRandomField.h"
#include "SeepageProblem.h"
#include "SlopeGeometry.h"
#include "SlopeStability.h"
#include "TPZVTKGeoMesh.h"
#include "pzgeoel.h"

namespace {

using Clock = std::chrono::steady_clock;
double Seconds(Clock::time_point t0) { return std::chrono::duration<double>(Clock::now() - t0).count(); }

/// Opções de linha de comando no formato chave=valor
struct TArgs {
    std::map<std::string, std::string> kv;
    TArgs(int argc, char **argv, int first) {
        for (int i = first; i < argc; i++) {
            std::string a = argv[i];
            auto p = a.find('=');
            if (p == std::string::npos) kv[a] = "1";
            else kv[a.substr(0, p)] = a.substr(p + 1);
        }
    }
    REAL Get(const std::string &k, REAL def) const {
        auto it = kv.find(k);
        return it == kv.end() ? def : std::atof(it->second.c_str());
    }
    int GetI(const std::string &k, int def) const {
        auto it = kv.find(k);
        return it == kv.end() ? def : std::atoi(it->second.c_str());
    }
    std::string GetS(const std::string &k, const std::string &def) const {
        auto it = kv.find(k);
        return it == kv.end() ? def : it->second;
    }
};

/// Gradiente do excesso de poropressão do problema de Darcy nos pontos de integração da malha mecânica
template <class T>
void TransferSeepage(TSeepageProblem &seep, TSlopeFEM<T> &fem) {
    const auto &pts = fem.Points();
    std::vector<REAL> u(pts.size(), 0.);
    std::vector<TPZManVector<REAL, 2>> g(pts.size(), TPZManVector<REAL, 2>(2, 0.));
    for (size_t i = 0; i < pts.size(); i++) {
        if (pts[i].gel < 0) continue;
        seep.Evaluate(pts[i].gel, pts[i].qsi, u[i], g[i]);
    }
    fem.SetSeepage(u, g);
}

struct TCase {
    std::string name;
    TSlopeGeometry geo;
    TSoil soil;
    bool seepage = false;
    REAL hw = 5.;
    REAL alpha = 1.;
};

TCase MakeCase(const std::string &name, REAL h) {
    TCase c;
    c.name = name;
    if (name == "cho_coesivo") {  // Cho (2010) / seção 5.3.1: cu = 23 kPa, φu = 0, γ = 20, 2:1, H = 5 m
        c.geo = TSlopeGeometry::Cho2H1V(h);
        c.soil.c = 23.;
        c.soil.phiDeg = 0.;
        c.soil.gamma = 20.;
        c.soil.buoyant = false;
    } else if (name == "cho_cphi") {
        // seção 5.3.2 / Cho (2010): c = 10 kPa, φ = 30°, γ = 20, 1:1 (seco). O FS = 1.204 de Cho (e Γ = 1.777 do
        // artigo) corresponde a H = 10 m: com H = 5 m o próprio Bishop simplificado dá FS = 1.61 (ver README).
        c.geo = TSlopeGeometry::Cho1H1V(h);
        c.geo.H = 10.;
        c.soil.c = 10.;
        c.soil.phiDeg = 30.;
        c.soil.gamma = 20.;
        c.soil.buoyant = false;
    } else {  // "percolacao": seção 6 (Tabela 2): c = 10, φ = 30°, γ = 20, β = 45°, H = hw = 5 m, α = 1
        c.geo = TSlopeGeometry::Cho1H1V(h);
        c.soil.c = 10.;
        c.soil.phiDeg = 30.;
        c.soil.gamma = 20.;
        c.soil.buoyant = true;
        c.seepage = true;
    }
    return c;
}

void AdaptMesh(const TCase &cs, TPZGeoMesh *gmesh, const TArgs &args);

/// Tensão efetiva geostática: análise elástica (Mohr-Coulomb com c muito alta) com a força de corpo (γ') e sem
/// percolação, na mesma malha e ordem; devolve σ'0 nos pontos dados (que devem coincidir com os da análise)
template <class TPoints>
std::vector<TPZTensor<REAL>> GeostaticStress(const TCase &cs, TPZGeoMesh *gmesh, const TSolverOptions &opt,
                                             const TPoints &target) {
    TSoil elastic = cs.soil;
    elastic.c = 1.e6;
    TSlopeFEM<TMohrCoulomb> fe(gmesh, cs.geo, elastic, opt);
    fe.ResetState();
    int its = 0;
    if (!fe.Solve(1., 0., 1., its)) throw std::runtime_error("análise geostática não convergiu");
    std::vector<TPZTensor<REAL>> s;
    fe.Stresses(s);
    const auto &pa = fe.Points();
    if (pa.size() != target.size()) throw std::runtime_error("pontos de integração diferentes (geostática)");
    for (size_t i = 0; i < pa.size(); i++)
        if (target[i].gel >= 0 && std::hypot(pa[i].x[0] - target[i].x[0], pa[i].x[1] - target[i].x[1]) > 1.e-9)
            throw std::runtime_error("pontos de integração diferentes (geostática)");
    return s;
}

/// Tensão inicial do Cam-Clay (Mohr-Coulomb: nada a fazer)
template <class T>
void GeostaticStress(const TCase &cs, TPZGeoMesh *gmesh, const TSolverOptions &opt, TSlopeFEM<T> &target) {
    if constexpr (std::is_same_v<T, TPZModifiedCamClay>) {
        target.SetInitialStress(GeostaticStress(cs, gmesh, opt, target.Points()));
    } else {
        (void)cs;
        (void)gmesh;
        (void)opt;
        (void)target;
    }
}

/// Estatísticas de um campo lognormal (média, CoV; CoV = 0: determinístico)
struct TFieldSpec {
    REAL mean = 1., cov = 0.;
};

template <class T>
void RunMonteCarlo(const TCase &cs, const TArgs &args, const std::string &modelName) {
    TSolverOptions opt;
    opt.porder = args.GetI("p", 2);
    opt.verbose = args.GetI("verbose", 0);
    opt.relTol = args.Get("reltol", 5.e-3);
    opt.maxIter = args.GetI("maxit", 20);
    opt.stagnation = args.GetI("estagnacao", 1) != 0;
    const std::string medida = args.GetS("medida", "gamma");  // gamma (fator de carga) ou fs (redução)
    const int64_t n = args.GetI("n", 100), first = args.GetI("inicio", 0);
    const uint64_t seed = (uint64_t)args.GetI("seed", 2025);
    TFieldSpec fc{cs.soil.c, args.Get("covc", 0.3)}, fphi{cs.soil.phiDeg, args.Get("covphi", 0.1)},
        fk{1., args.Get("covk", cs.seepage ? 0.6 : 0.)};
    const REAL Lx = args.Get("Lx", 20.), Ly = args.Get("Ly", 2.);

    std::unique_ptr<TPZGeoMesh> gmesh(cs.geo.CreateGeoMesh());
    AdaptMesh(cs, gmesh.get(), args);
    TSlopeFEM<T> fem(gmesh.get(), cs.geo, cs.soil, opt);
    GeostaticStress(cs, gmesh.get(), opt, fem);
    std::cout << "[MC] malha mecânica: " << fem.NEquations() << " equações, " << fem.NPoints()
              << " pontos de integração\n";

    // campo aleatório: malha KL própria (quadriláteros de 9 nós, h_KL)
    TSlopeGeometry geoKL = cs.geo;
    geoKL.h = args.Get("hkl", std::min(Ly / 2., 1.));
    geoKL.triangles = false;
    TPZKLRandomField::TOptions klopt;
    klopt.Lx = Lx;
    klopt.Ly = Ly;
    klopt.porder = 2;
    klopt.nModes = args.GetI("M", -1);
    klopt.targetVarianceError = args.Get("epsM", -1.);
    klopt.normalizeVariance = args.GetI("normvar", 1) != 0;
    {
        std::ostringstream f;
        f << "kl_H" << cs.geo.H << "_b" << std::lround(cs.geo.betaDeg * 100) << "_Lx" << Lx << "_Ly" << Ly << "_h"
          << geoKL.h << ".bin";
        klopt.cacheFile = args.GetS("klcache", f.str());
    }
    auto t0 = Clock::now();
    TPZKLRandomField kl(geoKL.CreateGeoMesh(), klopt);
    kl.Compute();
    std::cout << "[KL] " << kl.NEquations() << " equações, M = " << kl.NModes() << ", ε_M = "
              << kl.VarianceError(kl.NModes()) << " (" << Seconds(t0) << " s)\n";

    // pontos alvo: pontos de integração da malha mecânica (c, φ) e centróides dos elementos (kv)
    std::vector<TPZManVector<REAL, 3>> ipts;
    for (auto &p : fem.Points()) ipts.push_back(p.gel >= 0 ? p.x : TPZManVector<REAL, 3>(3, 0.));
    const int setIP = kl.AddTargetSet(ipts);
    std::vector<int64_t> elems;
    std::vector<TPZManVector<REAL, 3>> cpts;
    for (TPZGeoEl *gel : gmesh->ElementVec()) {
        if (!gel || gel->Dimension() != 2 || gel->HasSubElement()) continue;
        TPZManVector<REAL, 3> qc(2, 0.), xc(3, 0.);
        gel->CenterPoint(gel->NSides() - 1, qc);
        gel->X(qc, xc);
        elems.push_back(gel->Index());
        cpts.push_back(xc);
    }
    const int setEl = kl.AddTargetSet(cpts);
    {
        REAL vmin = 1., vmean = 0.;
        for (REAL v : kl.PointVariance(setIP)) {
            vmin = std::min(vmin, v);
            vmean += v / kl.PointVariance(setIP).size();
        }
        std::cout << "[KL] variância truncada nos pontos de integração: média " << vmean << ", mínima " << vmin
                  << (klopt.normalizeVariance ? " (compensada)" : " (não compensada)") << "\n";
    }

    std::unique_ptr<TSeepageProblem> seep;
    if (cs.seepage) {
        TSeepageProblem::TParams sp;
        sp.hw = cs.hw;
        sp.alpha = cs.alpha;
        sp.gammaW = cs.soil.gammaW;
        seep = std::make_unique<TSeepageProblem>(gmesh.get(), cs.geo, sp);
        if (fk.cov == 0.) {
            seep->Solve();
            TransferSeepage(*seep, fem);
        }
    }

    const std::string csv = args.GetS("saida", "mc_" + cs.name + "_" + modelName + ".csv");
    const bool exists = std::ifstream(csv).good();
    std::ofstream out(csv, std::ios::app);
    if (!exists) out << "amostra,fator,limite_superior,status,passos,iteracoes,cortes,tempo_s,c_medio,phi_medio,kv_medio\n";
    out << std::setprecision(8);
    const int vtkEvery = args.GetI("vtk", 0);
    REAL sum = 0., sum2 = 0.;
    int64_t nf = 0, nok = 0;
    std::cout << "[MC] " << cs.name << " (" << modelName << "), medida = " << medida << ", amostras " << first << ".."
              << first + n - 1 << ", CoV(c) = " << fc.cov << ", CoV(phi) = " << fphi.cov << ", CoV(kv) = " << fk.cov
              << ", Lx = " << Lx << ", Ly = " << Ly << " -> " << csv << "\n";
    std::vector<REAL> c(fem.NPoints()), phi(fem.NPoints());
    std::vector<REAL> kvEl(gmesh->NElements(), 1.);
    for (int64_t s = first; s < first + n; s++) {
        auto ts = Clock::now();
        std::vector<std::vector<REAL>> xi(3);
        for (int f = 0; f < 3; f++) kl.Xi(seed, s, f, xi[f]);
        std::vector<std::vector<std::vector<REAL>>> H;
        kl.Evaluate(xi, {setIP, setEl}, H);
        REAL cm = 0., pm = 0., km = 0.;
        int64_t np = 0;
        for (int64_t i = 0; i < fem.NPoints(); i++) {
            if (fem.Points()[i].gel < 0) continue;
            c[i] = TPZKLRandomField::Lognormal(H[0][setIP][i], fc.mean, fc.cov);
            const REAL phideg = TPZKLRandomField::Lognormal(H[1][setIP][i], fphi.mean, fphi.cov);
            phi[i] = phideg * M_PI / 180.;
            cm += c[i];
            pm += phideg;
            np++;
        }
        cm /= np;
        pm /= np;
        fem.SetStrength(c, phi);
        if (seep && fk.cov > 0.) {
            for (size_t e = 0; e < elems.size(); e++) {
                kvEl[elems[e]] = TPZKLRandomField::Lognormal(H[2][setEl][e], fk.mean, fk.cov);
                km += kvEl[elems[e]] / elems.size();
            }
            seep->SetElementPermeability(kvEl);
            seep->Solve();
            TransferSeepage(*seep, fem);
        } else {
            km = 1.;
        }
        fem.ResetState();
        TFactorResult r = (medida == "fs") ? fem.StrengthReduction(args.Get("F0", 0.5)) : fem.LoadFactor(0.);
        const double dt = Seconds(ts);
        out << s << "," << r.factor << "," << r.upper << "," << (r.status.empty() ? "ok" : r.status) << "," << r.steps
            << "," << r.iterations << "," << r.cuts << "," << dt << "," << cm << "," << pm << "," << km << "\n";
        out.flush();
        sum += r.factor;
        sum2 += r.factor * r.factor;
        nok++;
        if (r.factor < 1.) nf++;
        if (vtkEvery > 0 && (s - first) % vtkEvery == 0) {
            const std::string base = cs.name + "_" + modelName + "_amostra" + std::to_string(s);
            fem.DefineVTK(base + ".vtk");
            fem.WriteVTK(0);
            if (seep) {
                seep->DefineVTK(base + "_darcy.vtk");
                seep->WriteVTK(0);
            }
        }
        if (opt.verbose || (s - first) % 10 == 0)
            std::cout << "  amostra " << s << ": fator = " << r.factor << " (" << dt << " s), Pf parcial = "
                      << REAL(nf) / nok << "\n";
    }
    const REAL mean = sum / std::max<int64_t>(nok, 1);
    const REAL sd = std::sqrt(std::max(REAL(0.), sum2 / std::max<int64_t>(nok, 1) - mean * mean));
    std::cout << "[MC] " << nok << " amostras: média " << mean << ", desvio " << sd << ", Pf = " << REAL(nf) / nok
              << " (tempo " << Seconds(t0) << " s)\n";
}

/// Adapta o TPZGeoMesh ao mecanismo de colapso do problema com propriedades médias (Mohr-Coulomb): em cada nível
/// resolve Γ (ou FS), marca os elementos com ||ε^p|| > frac max e "camadas" de vizinhos e os divide.
void AdaptMesh(const TCase &cs, TPZGeoMesh *gmesh, const TArgs &args) {
    const int levels = args.GetI("adapt", 0);
    if (levels <= 0) return;
    const REAL frac = args.Get("frac", 0.02);
    const int layers = args.GetI("camadas", 1);
    const bool fs = args.GetS("medida", "gamma") == "fs";
    TSolverOptions opt;
    opt.porder = args.GetI("p", 2);
    opt.relTol = args.Get("reltol", 5.e-3);
    opt.maxIter = args.GetI("maxit", 20);
    opt.stagnation = args.GetI("estagnacao", 1) != 0;
    for (int l = 0; l < levels; l++) {
        auto t0 = Clock::now();
        TSlopeFEM<TMohrCoulomb> fem(gmesh, cs.geo, cs.soil, opt);
        std::unique_ptr<TSeepageProblem> seep;
        if (cs.seepage) {
            TSeepageProblem::TParams sp;
            sp.hw = cs.hw;
            sp.alpha = cs.alpha;
            sp.gammaW = cs.soil.gammaW;
            seep = std::make_unique<TSeepageProblem>(gmesh, cs.geo, sp);
            seep->Solve();
            TransferSeepage(*seep, fem);
        }
        fem.ResetState();
        TFactorResult r = fs ? fem.StrengthReduction(args.Get("F0", 0.5)) : fem.LoadFactor(0.);
        std::vector<REAL> ind;
        fem.PlasticIndicator(ind, args.GetI("incremento", 1) != 0);
        REAL vmax = 0.;
        for (REAL v : ind) vmax = std::max(vmax, v);
        std::vector<int64_t> marked;
        for (TPZGeoEl *gel : gmesh->ElementVec())
            if (gel && gel->Dimension() == 2 && !gel->HasSubElement() && ind[gel->Index()] > frac * vmax)
                marked.push_back(gel->Index());
        const size_t nmarked = marked.size();
        TSlopeGeometry::Grow(gmesh, marked, layers);
        std::cout << "[adapt] nível " << l << ": " << fem.NEquations() << " equações, " << (fs ? "FS" : "Gamma")
                  << " = " << r.factor << "; dividindo " << nmarked << " + " << marked.size() - nmarked
                  << " elementos (" << Seconds(t0) << " s)\n";
        TSlopeGeometry::Refine(gmesh, marked);
    }
}

void Usage(const char *prog) {
    std::cout << "uso: " << prog << " <comando> [chave=valor ...]\n"
              << "  det caso=cho_coesivo|cho_cphi|percolacao modelo=mc|mcc h=0.5 p=2 hw=5 alpha=1 beta=45 H= c= phi= gam=\n"
              << "      adapt=0 frac=0.02 camadas=1 gamma=1 fs=1 vtk=0 (Cam-Clay: lambda= kappa= v0= OCR=)\n"
              << "      Lc= Lt= Hb= (dimensões do domínio: crista, pé, base abaixo do pé)\n"
              << "  mc  caso=... modelo=mc|mcc medida=gamma|fs n=100 inicio=0 seed=2025 covc=0.3 covphi=0.1 covk=0.6\n"
              << "      Lx=20 Ly=2 hkl=1 M=-1 epsM=-1 normvar=1 saida=arquivo.csv vtk=0\n"
              << "  rebaixamento caso=percolacao modelo=mc|mcc Td=0.1 k=1e-5 Se=0 hw=5 alpha=1 tempos=0.01,0.1,1,3\n"
              << "      fs=1 gamma=0 saida=rebaixamento.csv vtk=0 (u-p acoplado; T = c_v t / H^2)\n";
}

/// Lista de números separados por vírgula
std::vector<REAL> ParseList(const std::string &str) {
    std::vector<REAL> v;
    std::stringstream ss(str);
    std::string item;
    while (std::getline(ss, item, ','))
        if (!item.empty()) v.push_back(std::atof(item.c_str()));
    return v;
}

/// Rebaixamento rápido/lento com o problema u-p acoplado (TPZMatPoroElastoPlastic3DMem em deformação plana):
/// em tempos escolhidos a poropressão é congelada e transferida à análise de estabilidade (Γ e/ou FS), como no
/// artigo, mas com p(x, t) do adensamento em vez do fluxo estacionário desacoplado.
template <class T>
void RunDrawdown(const TCase &cs, const TArgs &args, const std::string &modelName) {
    TSolverOptions opt;
    opt.porder = args.GetI("p", 2);
    opt.verbose = args.GetI("verbose", 0);
    opt.relTol = args.Get("reltol", 5.e-3);
    opt.maxIter = args.GetI("maxit", 20);
    opt.stagnation = args.GetI("estagnacao", 1) != 0;
    const bool doFS = args.GetI("fs", 1) != 0, doGamma = args.GetI("gamma", 0) != 0;
    std::unique_ptr<TPZGeoMesh> gmesh(cs.geo.CreateGeoMesh());
    AdaptMesh(cs, gmesh.get(), args);

    typename TCoupledDrawdown<T>::TParams par;
    par.hw = cs.hw;
    par.alpha = cs.alpha;
    par.kv = args.Get("k", 1.e-5);  // artigo: k_v/γw = 1e-6 m⁴/(kN s)
    par.Se = args.Get("Se", 0.);
    par.Td = args.Get("Td", 0.1);
    par.porderU = opt.porder;
    par.nDrawdown = args.GetI("nrebaixamento", 10);
    par.growth = args.Get("crescimento", 1.5);
    par.tol = args.Get("tolup", 1.e-7);
    par.verbose = args.GetI("verbose", 0);

    TSlopeFEM<T> fem(gmesh.get(), cs.geo, cs.soil, opt);
    GeostaticStress(cs, gmesh.get(), opt, fem);
    TCoupledDrawdown<T> cd(gmesh.get(), cs.geo, cs.soil, par);
    if constexpr (std::is_same_v<T, TPZModifiedCamClay>) cd.SetInitialStress(GeostaticStress(cs, gmesh.get(), opt, cd.Points()));
    std::cout << "\n=== rebaixamento acoplado (" << modelName << "): " << cs.geo.Describe() << "\n"
              << "    h_w = " << cs.hw << " m em T_d = " << par.Td << " (t_d = " << par.Td * cd.TimeScale() / 3600.
              << " h), k_v = " << par.kv << " m/s, k_h/k_v = " << par.alpha << ", c_v = " << cd.Cv()
              << " m2/s, H2/c_v = " << cd.TimeScale() / 3600. << " h\n"
              << "    u-p: " << cd.NEquations() << " equacoes; estabilidade: " << fem.NEquations() << " equacoes\n";

    const std::string out = args.GetS("saida", "rebaixamento_" + modelName + ".csv");
    std::ofstream csv(out);
    csv << "Td,T,t_h,zw,u_A,umax_desloc,pontos_plasticos,FS,FS_sup,FS_status,Gamma,Gamma_sup,Gamma_status,tempo_s\n";
    const bool vtk = args.GetI("vtk", 0) != 0;
    if (vtk) cd.DefineVTK(cs.name + "_" + modelName + "_Td" + args.GetS("Td", "0.1"));
    // ponto A: abaixo da borda da crista, a meia altura do talude
    TPZManVector<REAL, 3> xA = {cs.geo.Lc, cs.geo.Hb + 0.5 * cs.geo.H, 0.};
    int vtkStep = 0;
    auto evaluate = [&](const std::string &label) {
        auto t1 = Clock::now();
        cd.TransferSeepage(fem);
        TFactorResult rf, rg;
        if (doFS) {
            fem.ResetState();
            rf = fem.StrengthReduction(args.Get("F0", 0.5));
        }
        if (doGamma) {
            fem.ResetState();
            rg = fem.LoadFactor(0.);
        }
        const REAL Tnow = cd.Time();
        const REAL uA = cd.PorePressureAt(xA) - cs.soil.gammaW * (cs.geo.D() - xA[1]);
        std::cout << "    " << std::setw(10) << label << " T = " << std::setw(9) << Tnow << "  z_w = " << std::setw(6)
                  << cd.WaterLevel() << "  u_A = " << std::setw(8) << uA << "  |du|max = " << std::setw(10)
                  << cd.MaxDisplacementIncrement() << "  plast. " << std::setw(5) << cd.NPlasticPoints();
        if (doFS) std::cout << "  FS = " << rf.factor << " (" << rf.status << ")";
        if (doGamma) std::cout << "  Gamma = " << rg.factor << " (" << rg.status << ")";
        std::cout << "  [" << Seconds(t1) << " s]" << std::endl;
        csv << par.Td << "," << Tnow << "," << Tnow * cd.TimeScale() / 3600. << "," << cd.WaterLevel() << "," << uA << ","
            << cd.MaxDisplacementIncrement() << "," << cd.NPlasticPoints() << "," << rf.factor << "," << rf.upper
            << "," << rf.status << "," << rg.factor << "," << rg.upper << "," << rg.status << "," << Seconds(t1)
            << std::endl;
        if (vtk) cd.WriteVTK(vtkStep++);
    };

    auto t0 = Clock::now();
    // referência do artigo: fluxo estacionário desacoplado (Darcy) com as mesmas condições de contorno finais
    {
        TSeepageProblem::TParams sp;
        sp.hw = cs.hw;
        sp.alpha = cs.alpha;
        sp.gammaW = cs.soil.gammaW;
        TSeepageProblem seep(gmesh.get(), cs.geo, sp);
        seep.Solve();
        TransferSeepage(seep, fem);
        TFactorResult rf, rg;
        if (doFS) {
            fem.ResetState();
            rf = fem.StrengthReduction(args.Get("F0", 0.5));
        }
        if (doGamma) {
            fem.ResetState();
            rg = fem.LoadFactor(0.);
        }
        std::cout << "    estacionario desacoplado (artigo):";
        if (doFS) std::cout << " FS = " << rf.factor << " (" << rf.status << ")";
        if (doGamma) std::cout << " Gamma = " << rg.factor << " (" << rg.status << ")";
        std::cout << std::endl;
        csv << par.Td << ",inf,inf," << cs.geo.D() - cs.hw << ",nan,nan,0," << rf.factor << "," << rf.upper << ","
            << rf.status << "," << rg.factor << "," << rg.upper << "," << rg.status << ",0" << std::endl;
    }
    if (!cd.Initialize()) {
        std::cout << "    equilibrio inicial nao convergiu\n";
        return;
    }
    evaluate("inicial");
    std::vector<REAL> times;
    for (REAL f : {0.25, 0.5, 0.75, 1.}) times.push_back(f * par.Td);
    for (REAL t : ParseList(args.GetS("tempos", "0.01,0.03,0.1,0.3,1,3"))) times.push_back(par.Td + t);
    for (REAL Tout : times) {
        if (!cd.AdvanceTo(Tout, args.Get("dtmax", 0.5))) {
            std::cout << "    colapso no adensamento acoplado em T = " << cd.Time() << " (z_w = " << cd.WaterLevel()
                      << ")\n";
            csv << par.Td << "," << cd.Time() << "," << cd.Time() * cd.TimeScale() / 3600. << "," << cd.WaterLevel()
                << ",nan,nan,0,nan,nan,colapso_acoplado,nan,nan,colapso_acoplado,0" << std::endl;
            break;
        }
        evaluate(Tout <= par.Td * (1. + 1.e-9) ? "rebaixando" : "dissipacao");
    }
    std::cout << "    tempo total " << Seconds(t0) << " s\n";
}

template <class T>
void RunDeterministic(const TCase &cs, const TArgs &args, const std::string &modelName) {
    TSolverOptions opt;
    opt.porder = args.GetI("p", 2);
    opt.verbose = args.GetI("verbose", 0);
    opt.relTol = args.Get("reltol", 5.e-3);
    opt.maxIter = args.GetI("maxit", 20);
    opt.stagnation = args.GetI("estagnacao", 1) != 0;
    std::unique_ptr<TPZGeoMesh> gmesh(cs.geo.CreateGeoMesh());
    AdaptMesh(cs, gmesh.get(), args);
    std::cout << "\n=== " << cs.name << " (" << modelName << "): " << cs.geo.Describe() << "\n";
    std::cout << "    c = " << cs.soil.c << " kPa, phi = " << cs.soil.phiDeg << " graus, gamma = " << cs.soil.gamma
              << (cs.soil.buoyant ? " (forca de corpo gamma')" : "") << (cs.seepage ? ", percolacao hw = " : "")
              << (cs.seepage ? std::to_string(cs.hw) : std::string()) << "\n";
    auto t0 = Clock::now();
    TSlopeFEM<T> fem(gmesh.get(), cs.geo, cs.soil, opt);
    GeostaticStress(cs, gmesh.get(), opt, fem);
    std::cout << "    " << fem.NEquations() << " equacoes, " << fem.NPoints() << " pontos de integracao\n";
    if (args.GetI("vtk", 0)) {  // malha (elementos computacionais = folhas da malha adaptada)
        std::ofstream f(cs.name + "_malha.vtk");
        TPZVTKGeoMesh::PrintCMeshVTK(fem.Mesh(), f, true);
    }
    std::unique_ptr<TSeepageProblem> seep;
    if (cs.seepage) {
        TSeepageProblem::TParams sp;
        sp.hw = cs.hw;
        sp.alpha = cs.alpha;
        sp.gammaW = cs.soil.gammaW;
        seep = std::make_unique<TSeepageProblem>(gmesh.get(), cs.geo, sp);
        seep->Solve();
        TransferSeepage(*seep, fem);
        if (args.GetI("vtk", 0)) {
            seep->DefineVTK(cs.name + "_darcy.vtk");
            seep->WriteVTK(0);
        }
    }
    if (args.GetI("gamma", 1)) {
        fem.ResetState();
        auto t1 = Clock::now();
        TFactorResult r = fem.LoadFactor(0.);
        std::cout << "    Gamma (fator de carga) = " << r.factor << " (colapso em " << r.upper << "; " << r.steps
                  << " passos, " << r.iterations << " it., " << r.cuts << " cortes, " << Seconds(t1) << " s) "
                  << r.status << "\n";
        if (args.GetI("vtk", 0)) {
            fem.DefineVTK(cs.name + "_" + modelName + "_gamma.vtk");
            fem.WriteVTK(0);
        }
    }
    if (args.GetI("fs", 1)) {
        fem.ResetState();
        auto t1 = Clock::now();
        TFactorResult r = fem.StrengthReduction(args.Get("F0", 0.5));
        std::cout << "    FS (reducao de resistencia) = " << r.factor << " (colapso em " << r.upper << "; "
                  << r.steps << " passos, " << r.iterations << " it., " << r.cuts << " cortes, " << Seconds(t1)
                  << " s) " << r.status << "\n";
        if (args.GetI("vtk", 0)) {
            fem.DefineVTK(cs.name + "_" + modelName + "_fs.vtk");
            fem.WriteVTK(0);
        }
    }
    std::cout << "    tempo total " << Seconds(t0) << " s\n";
}

} // namespace

int main(int argc, char **argv) {
    if (argc < 2) {
        Usage(argv[0]);
        return 1;
    }
    const std::string cmd = argv[1];
    TArgs args(argc, argv, 2);
    if (cmd == "det" || cmd == "mc" || cmd == "rebaixamento") {
        TCase cs = MakeCase(args.GetS("caso", "percolacao"), args.Get("h", 0.5));
        cs.hw = args.Get("hw", cs.hw);
        cs.alpha = args.Get("alpha", cs.alpha);
        if (args.kv.count("beta")) cs.geo.betaDeg = args.Get("beta", 45.);
        if (args.kv.count("H")) {  // escala a geometria com H (crista, pé e base proporcionais)
            const REAL s = args.Get("H", cs.geo.H) / cs.geo.H;
            cs.geo.H *= s;
            cs.geo.Lc *= s;
            cs.geo.Lt *= s;
            cs.geo.Hb *= s;
            cs.geo.h *= s;
            if (!args.kv.count("hw")) cs.hw *= s;
        }
        // dimensões do domínio (crista, pé e base abaixo do pé), em m
        if (args.kv.count("Lc")) cs.geo.Lc = args.Get("Lc", cs.geo.Lc);
        if (args.kv.count("Lt")) cs.geo.Lt = args.Get("Lt", cs.geo.Lt);
        if (args.kv.count("Hb")) cs.geo.Hb = args.Get("Hb", cs.geo.Hb);
        if (args.kv.count("c")) cs.soil.c = args.Get("c", cs.soil.c);
        if (args.kv.count("phi")) cs.soil.phiDeg = args.Get("phi", cs.soil.phiDeg);
        if (args.kv.count("nu")) cs.soil.nu = args.Get("nu", cs.soil.nu);
        cs.soil.gamma = args.Get("gam", cs.soil.gamma);
        cs.soil.lambda = args.Get("lambda", cs.soil.lambda);
        cs.soil.kappa = args.Get("kappa", cs.soil.kappa);
        cs.soil.v0 = args.Get("v0", cs.soil.v0);
        cs.soil.OCR = args.Get("OCR", cs.soil.OCR);
        if (args.GetS("mapeamento", "deformacao_plana") == "triaxial")
            cs.soil.mapping = TPZModifiedCamClay::ETriaxialCompression;
        if (args.kv.count("E")) cs.soil.E = args.Get("E", cs.soil.E);
        cs.geo.triangles = args.GetI("tri", 0) != 0;
        const std::string model = args.GetS("modelo", "mc");
        if (cmd == "det") {
            if (model == "mcc") RunDeterministic<TPZModifiedCamClay>(cs, args, model);
            else RunDeterministic<TMohrCoulomb>(cs, args, model);
        } else if (cmd == "rebaixamento") {
            if (model == "mcc") RunDrawdown<TPZModifiedCamClay>(cs, args, model);
            else RunDrawdown<TMohrCoulomb>(cs, args, model);
        } else {
            if (model == "mcc") RunMonteCarlo<TPZModifiedCamClay>(cs, args, model);
            else RunMonteCarlo<TMohrCoulomb>(cs, args, model);
        }
        return 0;
    }
    Usage(argv[0]);
    return 1;
}
