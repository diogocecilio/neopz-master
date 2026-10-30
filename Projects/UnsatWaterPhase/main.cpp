// Teste 2 do plano: agua com ar passivo (p_a = 0 em todo o dominio) e esqueleto rigido (deslocamentos fixos).
// A equacao resolvida e a de Richards na forma secante condensada (documento, Secoes 13 e 14).
#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "UnsatProblem.h"
#include "UnsatReference.h"
#include "UnsatTaylor.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>

namespace {
const std::map<int, std::array<bool, 2>> kFix = {{EBase, {true, true}}, {ELeft, {true, false}}, {ERight, {true, false}}};

/// rigid skeleton and passive air: all displacement and air equations fixed (air at zero)
std::map<int64_t, REAL> RigidPassive(UnsatAnalysis &an) {
    std::map<int64_t, REAL> f = an.AllDisplacements(an.X());
    for (auto &e : an.AirEq()) f[e.second] = 0.;
    return f;
}
} // namespace

TEST_CASE("Teste 2.0: Jacobiano da fase agua (esqueleto rigido, ar passivo) com ordem 2", "[agua][jacobiano]") {
    UnsatParams p;
    p.hydro.fAVG = 0.1, p.hydro.fBVG = 2., p.hydro.fSwr = 0.05;
    UnsatProblem pr(RectangleMesh(0.2, 1., 2, 10), {{1, p}}, kFix);
    auto &an = pr.A();
    TPZFMatrix<STATE> &X = an.X();
    for (auto &e : an.WaterEq()) X(e.second, 0) = -10. - 8. * pr.Y(e.first);
    an.Restore();
    an.Accept(); // estado convergido (memoria) com p_w de X
    TPZFMatrix<STATE> x(X);
    unsigned k = 0;
    for (auto &e : an.WaterEq()) x(e.second, 0) += 3. * std::sin(1.7 * (++k));
    const auto fixed = RigidPassive(an);
    for (bool lumped : {true, false}) {
        for (auto *m : an.Materials()) m->SetLumpedStorage(lumped);
        const TaylorResult r = TaylorTest(an, x, 600., false, fixed, 0., 1.);
        std::cout << "Jacobiano da agua (condensado " << lumped << "): ordem mediana " << r.median << ", minima " << r.min << "\n";
        CHECK(r.min > 1.9);
        CHECK(r.min < 2.1);
    }
}

TEST_CASE("Teste 2.1: perfil hidrostatico e solucao exata (fluxo nulo, sem variacao)", "[agua]") {
    UnsatParams p;
    UnsatProblem pr(RectangleMesh(0.2, 2., 2, 20), {{1, p}}, kFix);
    auto &an = pr.A();
    TPZFMatrix<STATE> &X = an.X();
    for (auto &e : an.WaterEq()) X(e.second, 0) = -p.gw * pr.Y(e.first); // nivel d'agua na base
    an.Restore();
    an.Accept();
    auto fixed = RigidPassive(an);
    for (auto &n : an.BoundaryNodes(EBase)) fixed[an.WaterEq().at(n.first)] = 0.;
    REAL dmax = 0.;
    for (int s = 0; s < 5; ++s) {
        REQUIRE(an.Solve(3600., fixed, {}, false).conv);
        an.Accept();
    }
    for (auto &e : an.WaterEq()) dmax = std::max(dmax, std::fabs(an.X()(e.second, 0) + p.gw * pr.Y(e.first)));
    std::cout << "hidrostatico: max |p_w - p_w0| = " << dmax << " kPa\n";
    CHECK(dmax < 1.e-8);
}

TEST_CASE("Teste 2.2: infiltracao permanente em solo de Gardner contra a solucao fechada", "[agua]") {
    UnsatParams p;
    p.hydro.fRetention = TPZUnsaturatedHydraulics::EGardner;
    p.hydro.fAG = 0.2, p.hydro.fSwr = 0.;
    p.hydro.fKrw = TPZUnsaturatedHydraulics::EPowerLaw, p.hydro.fPower = 1.; // k = k_s exp(-a p_c)
    p.ks = 1.e-5;
    const REAL H = 3., W = 0.2, q = 0.3 * p.ks;
    REAL err[2];
    const int ny[2] = {30, 60};
    for (int m = 0; m < 2; ++m) {
        UnsatProblem pr(RectangleMesh(W, H, 1, ny[m]), {{1, p}}, kFix);
        auto &an = pr.A();
        TPZFMatrix<STATE> &X = an.X();
        for (auto &e : an.WaterEq()) X(e.second, 0) = -p.gw * pr.Y(e.first);
        an.Restore();
        an.Accept();
        auto fixed = RigidPassive(an);
        for (auto &n : an.BoundaryNodes(EBase)) fixed[an.WaterEq().at(n.first)] = 0.;
        std::map<int64_t, REAL> fq;
        for (auto &n : an.BoundaryNodes(ETop)) fq[an.WaterEq().at(n.first)] = q * n.second;
        REAL dt = 600., change = 1.;
        for (int s = 0; s < 200 && change > 1.e-9; ++s) {
            const TPZFMatrix<STATE> x0(an.X());
            REQUIRE(an.Solve(dt, fixed, fq, false).conv);
            an.Accept();
            change = 0.;
            for (auto &e : an.WaterEq()) change = std::max(change, std::fabs(an.X()(e.second, 0) - x0.GetVal(e.second, 0)));
            dt *= 1.5;
        }
        REQUIRE(change < 1.e-9);
        err[m] = 0.;
        REAL smax = 0.;
        for (auto &e : an.WaterEq()) {
            const REAL ref = unsatref::GardnerSteadySuction(pr.Y(e.first), p.hydro.fAG, p.gw, q, p.ks);
            err[m] = std::max(err[m], std::fabs(-an.X()(e.second, 0) - ref));
            smax = std::max(smax, ref);
        }
        std::cout << "Gardner: ny = " << ny[m] << ", erro max " << err[m] << " kPa, succao max " << smax << " kPa\n";
        CHECK(err[m] < 0.01 * smax);
    }
    INFO("razao dos erros (malha / 2) " << err[0] / err[1]);
    CHECK(err[0] / err[1] > 3.); // segunda ordem nos nos
}

TEST_CASE("Teste 2.3: balanco de massa de agua na infiltracao transiente (forma secante)", "[agua]") {
    UnsatParams p;
    p.hydro.fAVG = 0.1, p.hydro.fBVG = 2.5, p.hydro.fSwr = 0.05;
    p.Kw = 1.e15; // so o termo de saturacao armazena agua
    UnsatProblem pr(RectangleMesh(0.2, 1., 2, 20), {{1, p}}, kFix);
    auto &an = pr.A();
    TPZFMatrix<STATE> &X = an.X();
    for (auto &e : an.WaterEq()) X(e.second, 0) = -30.;
    an.Restore();
    an.Accept();
    const auto fixed = RigidPassive(an); // coluna fechada: so entra a chuva no topo
    std::map<int64_t, REAL> fq;
    const REAL q = 2.e-6;
    REAL Lw = 0.;
    for (auto &n : an.BoundaryNodes(ETop)) fq[an.WaterEq().at(n.first)] = q * n.second, Lw += n.second;
    REAL Vw0, Va0, Vw, Va;
    pr.StoredVolumes(an.X(), p.hydro, p.n, Vw0, Va0);
    REAL t = 0.;
    for (REAL dt : {60., 120., 300., 600., 1200., 1800., 3600., 3600., 3600.}) {
        REQUIRE(an.Solve(dt, fixed, fq, false).conv);
        an.Accept();
        t += dt;
    }
    pr.StoredVolumes(an.X(), p.hydro, p.n, Vw, Va);
    const REAL in = q * Lw * t;
    std::cout << "balanco de agua: armazenada " << Vw - Vw0 << ", entrada " << in << ", erro relativo "
              << std::fabs((Vw - Vw0) - in) / in << "\n";
    CHECK(std::fabs((Vw - Vw0) - in) < 1.e-6 * in);
}

namespace {
/// elevation where the nodal profile p_w(z) of the left edge crosses pstar (front moving down: p_w decreases with depth)
REAL FrontElevation(UnsatProblem &pr, REAL pstar) {
    auto &an = pr.A();
    std::vector<std::pair<REAL, REAL>> zp; // (z, p_w) on x = 0
    for (auto &e : an.WaterEq())
        if (pr.X(e.first) < 1.e-9) zp.push_back({pr.Y(e.first), an.X()(e.second, 0)});
    std::sort(zp.begin(), zp.end());
    for (size_t i = 1; i < zp.size(); ++i) {
        const REAL p0 = zp[i - 1].second, p1 = zp[i].second;
        if ((p0 - pstar) * (p1 - pstar) <= 0. && p0 != p1)
            return zp[i - 1].first + (pstar - p0) * (zp[i].first - zp[i - 1].first) / (p1 - p0);
    }
    return -1.;
}

/// Callari e Abati (2009, Secao 4.2): coluna rigida de 0,30 m, p_w = -100 kPa inicial, -7,5 kPa no topo e -100 kPa na
/// base a partir de t = 0 (Fig. 5, p. 7), Tabela 2 (p. 8), 1000 passos de 1 s e depois passos de 10 s (p. 9); ny
/// elementos na altura e passos multiplicados por dtf. Devolve o perfil em t = 21600 s pelas elevacoes da frente.
std::map<int, REAL> Philip(int ny, REAL dtf) {
    UnsatParams p;
    p.n = 0.368, p.Kw = 1.e12;         // liquido incompressivel
    p.ks = 0.922e-5 * p.gw;            // k_w^sat = 0,922e-5 m2/(kPa s)
    p.hydro.fRetention = TPZUnsaturatedHydraulics::EVanGenuchten;
    p.hydro.fAVG = 0.335, p.hydro.fBVG = 2.01, p.hydro.fSwr = 0.277;
    p.hydro.fKrw = TPZUnsaturatedHydraulics::EMualem;
    UnsatProblem pr(RectangleMesh(0.03, 0.30, 1, ny), {{1, p}}, kFix);
    auto &an = pr.A();
    TPZFMatrix<STATE> &X = an.X();
    for (auto &e : an.WaterEq()) X(e.second, 0) = -100.;
    an.Restore();
    an.Accept();
    auto fixed = RigidPassive(an);
    for (auto &n : an.BoundaryNodes(ETop)) fixed[an.WaterEq().at(n.first)] = -7.5;
    for (auto &n : an.BoundaryNodes(EBase)) fixed[an.WaterEq().at(n.first)] = -100.;
    REAL t = 0.;
    while (t < 21600. - 1.e-9) {
        const REAL dt = (t < 1000. - 1.e-9 ? 1. : 10.) * dtf;
        REQUIRE(an.Solve(dt, fixed, {}, false).conv);
        an.Accept();
        t += dt;
    }
    std::map<int, REAL> zf;
    for (int ps : {-30, -50, -70, -90}) zf[ps] = FrontElevation(pr, ps);
    return zf;
}

/// elevation of the digitized semi-analytical curve at p_w = pstar (data sorted by elevation, from the top)
REAL ReferenceElevation(const std::vector<std::array<REAL, 2>> &c, REAL pstar) {
    for (size_t i = 1; i < c.size(); ++i) {
        const REAL p0 = c[i - 1][0], p1 = c[i][0];
        if ((p0 - pstar) * (p1 - pstar) <= 0. && p0 != p1)
            return c[i - 1][1] + (pstar - p0) * (c[i][1] - c[i - 1][1]) / (p1 - p0);
    }
    return -1.;
}
} // namespace

namespace {
std::vector<std::array<REAL, 2>> PhilipReference() { // (p_w, z), digitalizado da Fig. 7 de Callari e Abati (2009)
    std::vector<std::array<REAL, 2>> ref;
    std::ifstream f(std::string(UNSAT_SOURCE_DIR) + "/data/philip_callari_fig7.csv");
    std::string l;
    while (std::getline(f, l)) {
        if (l.empty() || l[0] == '#') continue;
        std::replace(l.begin(), l.end(), ',', ' ');
        std::istringstream s(l);
        REAL pw, z;
        s >> pw >> z;
        ref.push_back({pw, z});
    }
    return ref;
}
const REAL kPhilipH = 0.30 / 60.; // um elemento da malha fina de Callari e Abati (2 x 6 x 60)
} // namespace

TEST_CASE("Teste 2.4: infiltracao de Philip (Callari e Abati 2009, Secao 4.2) contra a solucao semianalitica", "[agua]") {
    const auto ref = PhilipReference();
    REQUIRE(ref.size() > 20);
    const auto zf = Philip(60, 1.);
    REAL err = 0.;
    for (auto &e : zf) {
        const REAL zr = ReferenceElevation(ref, e.first);
        REQUIRE(zr > 0.);
        std::cout << "Philip, t = 21600 s, p_w = " << e.first << " kPa: z = " << e.second << " m, semianalitica " << zr
                  << " m\n";
        err = std::max(err, std::fabs(e.second - zr));
    }
    std::cout << "Philip: maior distancia vertical a solucao semianalitica " << err << " m (elemento " << kPhilipH << " m)\n";
    CHECK(err < kPhilipH);
}

TEST_CASE("Teste 2.4 (lento): Philip estavel com a malha refinada e o passo dividido por 2", "[agua][.lento]") {
    const auto base = Philip(60, 1.), fine = Philip(120, 0.5);
    REAL d = 0.;
    for (auto &e : base) {
        std::cout << "Philip, p_w = " << e.first << " kPa: z = " << e.second << " m (malha 60), " << fine.at(e.first)
                  << " m (malha 120, dt/2)\n";
        d = std::max(d, std::fabs(e.second - fine.at(e.first)));
    }
    std::cout << "Philip: mudanca com a malha e o passo refinados " << d << " m\n";
    CHECK(d < kPhilipH);
}
