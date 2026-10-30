// Testes 0 e 4 do plano: Jacobiano do modelo completo (com plasticidade), Terzaghi, Mandel e Liakopoulos.
#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "UnsatFS.h"
#include "UnsatProblem.h"
#include "UnsatReference.h"
#include "UnsatTaylor.h"

#include <fstream>
#include <iostream>
#include <sstream>

namespace {
/// saturado e elastico, sem peso (Terzaghi e Mandel)
UnsatParams SaturatedElastic() {
    UnsatParams p;
    p.E = 1.e4, p.nu = 0.3, p.c = 1.e8, p.gamma = 0., p.g = {0., 0., 0.};
    p.n = 0.3, p.Kw = 2.2e6, p.ks = 1.e-6;
    p.hydro.fAVG = 1.e-8; // S_w = 1 para qualquer succao pequena
    return p;
}
REAL Constrained(const UnsatParams &p) { return p.E * (1. - p.nu) / ((1. + p.nu) * (1. - 2. * p.nu)); }
REAL BiotModulus(const UnsatParams &p) { return p.Kw / p.n; } // 1/M = n/K_w (graos incompressiveis)
/// coeficiente de adensamento c = (k/gamma_w) M (lambda + 2G) / (lambda + 2G + alpha^2 M)
REAL Consolidation(const UnsatParams &p) {
    const REAL M = BiotModulus(p), L = Constrained(p);
    return p.ks / p.gw * M * L / (L + p.alpha * p.alpha * M);
}
} // namespace

TEST_CASE("Teste 0: Jacobiano do modelo completo (u, p_w, p_a, plasticidade) com ordem 2", "[acoplado][jacobiano]") {
    UnsatParams p;
    p.E = 2.e4, p.nu = 0.1, p.c = 1., p.phi = 30. * M_PI / 180., p.gamma = 20.;
    p.hydro.fAVG = 0.1, p.hydro.fBVG = 2.5, p.hydro.fSwr = 0.05;
    const std::map<int, std::array<bool, 2>> fix = {{EBase, {true, true}}, {ELeft, {true, false}}, {ERight, {true, false}}};
    // variantes: chi = S_w; chi constante; chi = S_w com o peso da mistura (d gamma / d p), graos compressiveis
    // (eta > 0, alpha = 1 - K/K_s) e armazenamento integrado nos pontos (sem condensacao); saturado (a_VG = 1e-8:
    // S_w = 1 e Gamma = 0 com as succoes do teste)
    struct Variant {
        const char *name;
        TUnsat::EBishop bishop;
        bool mixture, lumped;
        REAL Ks, avg;
    };
    const Variant variants[] = {{"chi S_w", TUnsat::EChiSaturation, false, true, 0., 0.1},
                                {"chi constante", TUnsat::EChiConstant, false, true, 0., 0.1},
                                {"chi S_w, mistura, eta > 0, consistente", TUnsat::EChiSaturation, true, false, 5.e4, 0.1},
                                {"saturado", TUnsat::EChiSaturation, false, true, 0., 1.e-8}};
    for (const Variant &v : variants) {
        p.bishop = v.bishop, p.mixture = v.mixture, p.lumped = v.lumped, p.Ks = v.Ks, p.hydro.fAVG = v.avg;
        p.alpha = v.Ks > 0. ? 1. - p.E / (3. * (1. - 2. * p.nu)) / v.Ks : 1.;
        UnsatProblem pr(RectangleMesh(1., 2., 2, 8), {{1, p}}, fix);
        auto &an = pr.A();
        TPZFMatrix<STATE> &X = an.X();
        for (auto &e : an.WaterEq()) X(e.second, 0) = -12. - 4. * pr.Y(e.first);
        for (auto &e : an.AirEq()) X(e.second, 0) = 1.;
        an.Restore();
        REQUIRE(Geostatic(an, 1.)); // peso com pressoes congeladas: pontos plasticos (K0 baixo)
        const int64_t npl = an.Materials()[0]->NPlasticPoints();
        // ponto-base fora do estado convergido: no estado convergido os pontos plasticos estao sobre a superficie
        // (Phi = 0) e a resposta nao e derivavel (carga ou descarga); um incremento finito de deslocamento os tira dela
        TPZFMatrix<STATE> x(an.X());
        unsigned k = 0;
        std::map<int64_t, REAL> fixed = an.DisplacementConstraints();
        for (int64_t i = 0; i < an.NEq(); ++i)
            if (!an.IsPressure(i) && !fixed.count(i)) x(i, 0) += 2.e-4 * std::sin(0.37 * (++k));
        for (auto &e : an.WaterEq()) x(e.second, 0) += 2. * std::sin(1.3 * (++k));
        for (auto &e : an.AirEq()) x(e.second, 0) += 1. * std::cos(0.7 * (++k));
        for (auto &n : an.BoundaryNodes(ETop)) fixed[an.AirEq().at(n.first)] = x(an.AirEq().at(n.first), 0);
        const TaylorResult r = TaylorTest(an, x, 600., false, fixed, 1.e-5, 1., 6, 1.e-3);
        std::cout << "Jacobiano completo (" << v.name << "): pontos plasticos " << npl << ", ordem mediana " << r.median
                  << ", minima " << r.min << "\n";
        CHECK(npl > 0);
        CHECK(r.median > 1.9);
        CHECK(r.median < 2.1);
    }
}

TEST_CASE("Teste 4a: Terzaghi no limite saturado contra a serie", "[acoplado]") {
    const UnsatParams p = SaturatedElastic();
    const REAL H = 1., q = 10.;
    UnsatProblem pr(RectangleMesh(0.1, H, 1, 40), {{1, p}},
                    {{EBase, {true, true}}, {ELeft, {true, false}}, {ERight, {true, false}}}, {{ETop, 0., -q}});
    auto &an = pr.A();
    an.SetForceReference(q * 0.1);
    std::map<int64_t, REAL> fixed = an.DisplacementConstraints();
    for (auto &e : an.AirEq()) fixed[e.second] = 0.;
    for (auto &n : an.BoundaryNodes(ETop)) fixed[an.WaterEq().at(n.first)] = 0.;
    const REAL c = Consolidation(p), M = BiotModulus(p), L = Constrained(p);
    const REAL p0 = p.alpha * M / (L + p.alpha * p.alpha * M) * q; // resposta nao drenada
    REQUIRE(an.Solve(1.e-6, fixed, {}, false).conv);
    an.Accept();
    REAL T = 0., worst = 0.;
    const REAL dT = 1.e-3;
    for (REAL Tcheck : {0.05, 0.1, 0.2, 0.5}) {
        while (T < Tcheck - 1.e-12) {
            REQUIRE(an.Solve(dT * H * H / c, fixed, {}, false).conv);
            an.Accept();
            T += dT;
        }
        REAL e = 0.;
        for (auto &w : an.WaterEq()) {
            const REAL z = (H - pr.Y(w.first)) / H;
            e = std::max(e, std::fabs(an.X()(w.second, 0) / p0 - unsatref::Terzaghi(z, T)));
        }
        std::cout << "Terzaghi: T = " << Tcheck << ", erro max de p/p0 = " << e << "\n";
        worst = std::max(worst, e);
    }
    CHECK(worst < 0.01);
}

TEST_CASE("Teste 4b: Mandel (placa rigida lisa) contra a solucao analitica", "[acoplado]") {
    const UnsatParams p = SaturatedElastic();
    const REAL a = 1., b = 1., sigma0 = -10.;
    UnsatProblem pr(StructuredMesh(Join({Graded(0., 0.7, 10), Graded(0.7, a, 10, 0.4)}), [b](REAL) { return b; },
                                   Graded(0., 1., 8)),
                    {{1, p}}, {{ELeft, {true, false}}, {EBase, {false, true}}, {ETop, {false, true}}});
    auto &an = pr.A();
    an.SetForceReference(std::fabs(sigma0) * a);
    std::map<int64_t, REAL> fixed = an.DisplacementConstraints();
    for (auto &e : an.AirEq()) fixed[e.second] = 0.;
    for (auto &n : an.BoundaryNodes(ERight)) fixed[an.WaterEq().at(n.first)] = 0.;
    const auto topy = an.ConstrainedEqs(ETop, 1);
    REQUIRE(!topy.empty());
    // base hierarquica do NeoPZ: deslocamento uniforme v da placa = v nos vertices e zero nas funcoes de aresta
    std::set<int64_t> topv;
    for (auto &n : an.BoundaryNodes(ETop)) topv.insert(an.DispEq().at(n.first) + 1);
    auto plate = [&](std::map<int64_t, REAL> &f, REAL vt) {
        for (auto i : topy) f[i] = topv.count(i) ? vt : 0.;
    };
    const REAL M = BiotModulus(p), L = Constrained(p);
    const REAL G = p.E / (2. * (1. + p.nu)), K = p.E / (3. * (1. - 2. * p.nu));
    const REAL B = p.alpha * M / (K + p.alpha * p.alpha * M);
    const REAL nuu = (3. * p.nu + p.alpha * B * (1. - 2. * p.nu)) / (3. - p.alpha * B * (1. - 2. * p.nu));
    const REAL c = Consolidation(p), p0 = std::fabs(sigma0) * B * (1. + nuu) / 3.;
    // mesmo coeficiente com o modulo volumetrico drenado K no lugar do confinado lambda + 2G
    const REAL cKdr = p.ks / p.gw / (1. / M + p.alpha * p.alpha / K);
    std::cout << "Mandel: c = " << c << " m2/s (modulo confinado), " << cKdr << " m2/s (K_dr volumetrico)\n";
    (void)G, (void)L;
    // passo com o deslocamento v da placa tal que a forca vertical total seja sigma0 a (problema linear: duas
    // solucoes dao v exato)
    REAL v = 0.;
    auto step = [&](REAL dt) {
        auto force = [&](REAL vt) {
            auto f = fixed;
            plate(f, vt);
            if (!an.Solve(dt, f, {}, false).conv) return std::numeric_limits<REAL>::quiet_NaN();
            const TPZFMatrix<STATE> r = an.FullResidual();
            REAL s = 0.;
            for (auto i : topv) s += r.GetVal(i, 0); // particao da unidade: so os vertices
            an.Restore();
            return s;
        };
        const REAL v1 = v, v2 = v - 1.e-4;
        const REAL f1 = force(v1), f2 = force(v2);
        // soma dos residuos nas equacoes verticais da placa = forca vertical total = sigma0 a (compressao < 0)
        v = v1 + (sigma0 * a - f1) * (v2 - v1) / (f2 - f1);
        auto f = fixed;
        plate(f, v);
        REQUIRE(an.Solve(dt, f, {}, false).conv);
        const TPZFMatrix<STATE> r = an.FullResidual();
        REAL s = 0.;
        for (auto i : topv) s += r.GetVal(i, 0); // particao da unidade: so os vertices
        REQUIRE(s == Approx(sigma0 * a).epsilon(1.e-6));
        an.Accept();
    };
    step(1.e-6);
    REAL pc0 = 0.;
    for (auto &w : an.WaterEq())
        if (pr.X(w.first) < 1.e-9) pc0 = an.X()(w.second, 0);
    REAL T = 0., worst = 0., pcmax = pc0;
    const REAL dT = 2.e-3;
    for (REAL Tcheck : {0.01, 0.05, 0.1, 0.3, 0.6}) {
        while (T < Tcheck - 1.e-12) {
            step(dT * a * a / c);
            T += dT;
            for (auto &w : an.WaterEq())
                if (pr.X(w.first) < 1.e-9) pcmax = std::max(pcmax, an.X()(w.second, 0));
        }
        REAL e = 0., ec = 0., ek = 0.;
        for (auto &w : an.WaterEq()) {
            const REAL ref = unsatref::MandelPressure(pr.X(w.first), T * a * a / c, a, sigma0, B, p.nu, nuu, c);
            const REAL d = std::fabs(an.X()(w.second, 0) - ref) / p0;
            e = std::max(e, d);
            if (pr.X(w.first) < 1.e-9) {
                ec = std::max(ec, d);
                // coeficiente de Haghighat et al. (2021, arXiv:2110.03049, p. 25) com K_dr volumetrico (p. 5)
                const REAL refk = unsatref::MandelPressure(0., T * a * a / c, a, sigma0, B, p.nu, nuu, cKdr);
                ek = std::fabs(an.X()(w.second, 0) - refk) / p0;
            }
        }
        std::cout << "Mandel: ct/a^2 = " << Tcheck << ", erro de p/p0 no centro " << ec << ", maximo no dominio " << e
                  << "; no centro com c(K_dr volumetrico) " << ek << "\n";
        worst = std::max(worst, ec); // criterio do plano: pressao no centro
    }
    std::cout << "Mandel: p0 = " << p0 << " kPa, p no centro apos o carregamento " << pc0 << ", maximo " << pcmax
              << " (efeito Mandel-Cryer)\n";
    CHECK(pcmax > 1.02 * pc0);
    CHECK(worst < 0.02);
}

namespace {
struct LiakData {
    std::vector<std::array<REAL, 2>> q;  ///< t, q (1e-6 m/s)
    std::vector<std::array<REAL, 3>> pw; ///< t, z, p_w
};
LiakData ReadData(const std::string &flux, const std::string &press) {
    LiakData d;
    std::ifstream f(flux), g(press);
    std::string l;
    while (std::getline(f, l)) {
        if (l.empty() || l[0] == '#') continue;
        std::replace(l.begin(), l.end(), ',', ' ');
        std::istringstream s(l);
        REAL t, q;
        s >> t >> q;
        d.q.push_back({t, q});
    }
    while (std::getline(g, l)) {
        if (l.empty() || l[0] == '#') continue;
        std::replace(l.begin(), l.end(), ',', ' ');
        std::istringstream s(l);
        REAL t, z, v;
        s >> t >> z >> v;
        d.pw.push_back({t, z, v});
    }
    return d;
}
LiakData ReadCallari(const std::string &file) {
    LiakData d;
    std::ifstream f(file);
    std::string l;
    while (std::getline(f, l)) {
        if (l.empty() || l[0] == '#') continue;
        std::replace(l.begin(), l.end(), ',', ' ');
        std::istringstream s(l);
        std::string k;
        s >> k;
        if (k == "Q") {
            REAL t, q;
            s >> t >> q;
            d.q.push_back({t, q});
        } else {
            REAL t, z, v;
            s >> t >> z >> v;
            d.pw.push_back({t, z, v});
        }
    }
    return d;
}
REAL Interp(const std::vector<std::array<REAL, 2>> &c, REAL t) {
    if (t <= c.front()[0]) return c.front()[1];
    for (size_t i = 1; i < c.size(); ++i)
        if (t <= c[i][0]) return c[i - 1][1] + (c[i][1] - c[i - 1][1]) * (t - c[i - 1][0]) / (c[i][0] - c[i - 1][0]);
    return c.back()[1];
}

/// Liakopoulos (Callari e Abati 2009, Tabela 1, Eq. 58): curvas de vazao (t, q) e p_w por (t, z)
struct LiakResult {
    std::vector<std::array<REAL, 2>> q;
    std::map<std::pair<int, int>, REAL> pw; ///< (t s, z cm) -> p_w
};
LiakResult RunLiakopoulos(bool activeAir) {
    UnsatParams p;
    p.E = 1300., p.nu = 0.4, p.c = 1.e5, p.gamma = 16.7, p.n = 0.2975, p.Kw = 2.e6;
    p.ks = 4.5e-7 * 9.81; // k_w = 4,5e-7 m2/(kPa s)
    p.hydro.fRetention = TPZUnsaturatedHydraulics::ELiakopoulos;
    p.hydro.fLC = 3.79010e-4, p.hydro.fLB = 2.4279, p.hydro.fSwr = 0.2;
    p.hydro.fKrw = TPZUnsaturatedHydraulics::ELiakopoulosW, p.hydro.fLD = 2.207, p.hydro.fLE = 1.0121;
    p.hydro.fKra = TPZUnsaturatedHydraulics::EBrooksCorey, p.hydro.fLambdaBC = 3.;
    const REAL W = 0.1;
    UnsatProblem pr(RectangleMesh(W, 1., 2, 40), {{1, p}},
                    {{EBase, {true, true}}, {ELeft, {true, false}}, {ERight, {true, false}}});
    auto &an = pr.A();
    REQUIRE(Geostatic(an, 1.)); // p_w = p_a = 0 congeladas
    an.ResetDisplacements();
    an.Accept();
    std::map<int64_t, REAL> fixed = an.DisplacementConstraints();
    const auto base = an.BoundaryNodes(EBase);
    for (auto &n : base) fixed[an.WaterEq().at(n.first)] = 0.;
    if (activeAir) for (int id : {EBase, ETop}) for (auto &n : an.BoundaryNodes(id)) fixed[an.AirEq().at(n.first)] = 0.;
    else for (auto &e : an.AirEq()) fixed[e.second] = 0.;
    LiakResult res;
    const REAL dt = 10.;
    for (int s = 1; s <= 720; ++s) {
        REQUIRE(an.Solve(dt, fixed, {}, false).conv);
        const TPZFMatrix<STATE> r = an.FullResidual();
        REAL in = 0.;
        for (auto &n : base) in += r.GetVal(an.WaterEq().at(n.first), 0);
        an.Accept();
        const REAL t = s * dt;
        res.q.push_back({t, -in / dt / W * 1.e6});
        if (s % 60 == 0 && s <= 180) {
            std::map<int, std::pair<REAL, int>> acc;
            for (auto &w : an.WaterEq()) {
                const int zc = int(std::lround(pr.Y(w.first) * 100.));
                acc[zc].first += an.X()(w.second, 0), acc[zc].second++;
            }
            for (auto &z : acc) res.pw[{int(t), z.first}] = z.second.first / z.second.second;
        }
    }
    return res;
}
} // namespace

TEST_CASE("Teste 4c: Liakopoulos com ar passivo contra os dados de Callari e Abati (2009, Fig. 3)", "[acoplado]") {
    const std::string dir = std::string(UNSAT_SOURCE_DIR) + "/data/";
    const LiakData exp = ReadData(dir + "liakopoulos_vazao_exp.csv", dir + "liakopoulos_pw_exp.csv");
    const LiakData cal = ReadCallari(dir + "liakopoulos_callari_calc.csv");
    const LiakResult r = RunLiakopoulos(false);
    REAL eq = 0., eqc = 0., ep = 0., epc = 0.;
    int nq = 0, np = 0;
    for (auto &d : exp.q) {
        if (d[0] < 10.) continue;
        eq += std::pow(Interp(r.q, d[0]) - d[1], 2), eqc += std::pow(Interp(cal.q, d[0]) - d[1], 2), ++nq;
    }
    for (auto &d : exp.pw) {
        const auto it = r.pw.find({int(d[0]), int(std::lround(d[1] * 100.))});
        REQUIRE(it != r.pw.end());
        REAL c = 0.;
        for (auto &k : cal.pw)
            if (k[0] == d[0] && std::fabs(k[1] - d[1]) < 1.e-9) c = k[2];
        ep += std::pow(it->second - d[2], 2), epc += std::pow(c - d[2], 2), ++np;
    }
    eq = std::sqrt(eq / nq), eqc = std::sqrt(eqc / nq), ep = std::sqrt(ep / np), epc = std::sqrt(epc / np);
    std::cout << "Liakopoulos, ar passivo: RMS da vazao " << eq << "e-6 m/s (Callari " << eqc << "), RMS de p_w " << ep
              << " kPa (Callari " << epc << ")\n";
    for (int t : {600, 1200, 1800}) {
        std::cout << "  t = " << t << " s, p_w (kPa) em z = 0,8/0,5/0,2 m: " << r.pw.at({t, 80}) << " / "
                  << r.pw.at({t, 50}) << " / " << r.pw.at({t, 20}) << "\n";
    }
    CHECK(eq <= 1.5 * eqc);
    CHECK(ep <= 1.5 * epc);
}

TEST_CASE("Teste 4d: Liakopoulos com ar ativo (comparacao qualitativa)", "[acoplado][.lento]") {
    const LiakResult pas = RunLiakopoulos(false), act = RunLiakopoulos(true);
    for (REAL t : {300., 1200., 3600., 7200.})
        std::cout << "t = " << t << " s: vazao passivo " << Interp(pas.q, t) << ", ativo " << Interp(act.q, t)
                  << " (1e-6 m/s)\n";
    CHECK(act.q.size() == pas.q.size());
}


