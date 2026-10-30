// Teste 3 do plano: agua e ar com esqueleto rigido (deslocamentos fixos). Armazenamento do ar com S_a^n
// (documento, Secao "Armazenamento do ar").
#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "UnsatProblem.h"
#include "UnsatTaylor.h"

#include <iostream>

namespace {
const std::map<int, std::array<bool, 2>> kFix = {{EBase, {true, true}}, {ELeft, {true, false}}, {ERight, {true, false}}};

UnsatParams Soil() {
    UnsatParams p;
    p.hydro.fAVG = 0.1, p.hydro.fBVG = 2.5, p.hydro.fSwr = 0.05;
    p.hydro.fKraMin = 1.e-4;
    p.ks = 1.e-5;
    return p;
}
} // namespace

TEST_CASE("Teste 3.0: Jacobiano agua-ar (esqueleto rigido) com ordem 2", "[ar][jacobiano]") {
    UnsatParams p = Soil();
    UnsatProblem pr(RectangleMesh(0.2, 1., 2, 10), {{1, p}}, kFix);
    auto &an = pr.A();
    TPZFMatrix<STATE> &X = an.X();
    for (auto &e : an.WaterEq()) X(e.second, 0) = -12. - 6. * pr.Y(e.first);
    for (auto &e : an.AirEq()) X(e.second, 0) = 2. + pr.Y(e.first);
    an.Restore();
    an.Accept();
    TPZFMatrix<STATE> x(X);
    unsigned k = 0;
    for (auto &e : an.WaterEq()) x(e.second, 0) += 2. * std::sin(1.3 * (++k));
    for (auto &e : an.AirEq()) x(e.second, 0) += 1.5 * std::cos(0.7 * (++k));
    std::map<int64_t, REAL> fixed = an.AllDisplacements(X);
    for (auto &n : an.BoundaryNodes(ETop)) fixed[an.AirEq().at(n.first)] = x(an.AirEq().at(n.first), 0);
    for (bool lumped : {true, false}) {
        for (auto *m : an.Materials()) m->SetLumpedStorage(lumped);
        const TaylorResult r = TaylorTest(an, x, 600., false, fixed, 0., 1.);
        std::cout << "Jacobiano agua-ar (condensado " << lumped << "): ordem mediana " << r.median << ", minima " << r.min << "\n";
        CHECK(r.min > 1.9);
        CHECK(r.min < 2.1);
    }
}

TEST_CASE("Teste 3a: ar confinado (Boyle): p_abs V_a constante com S_a^n", "[ar]") {
    UnsatParams p = Soil();
    p.Kw = 1.e15;
    UnsatProblem pr(RectangleMesh(0.1, 1., 1, 20), {{1, p}}, kFix);
    auto &an = pr.A();
    TPZFMatrix<STATE> &X = an.X();
    for (auto &e : an.WaterEq()) X(e.second, 0) = -20. - p.gw * pr.Y(e.first);
    an.Restore();
    an.Accept();
    const auto fixed = an.AllDisplacements(an.X()); // ar fechado em todo o contorno
    REAL Vw0, Va0, M0;
    pr.StoredVolumes(an.X(), p.hydro, p.n, Vw0, Va0, &M0, p.patm);
    // agua injetada pela base: 30% do volume inicial de ar em 10 h
    REAL Lb = 0.;
    for (auto &n : an.BoundaryNodes(EBase)) Lb += n.second;
    const REAL q = 0.3 * Va0 / (10. * 3600.) / Lb;
    std::map<int64_t, REAL> fq;
    for (auto &n : an.BoundaryNodes(EBase)) fq[an.WaterEq().at(n.first)] = q * n.second;
    REAL worst = 0., pamax = 0.;
    for (int s = 0; s < 40; ++s) {
        REQUIRE(an.Solve(900., fixed, fq, false).conv);
        an.Accept();
        REAL Vw, Va, M;
        pr.StoredVolumes(an.X(), p.hydro, p.n, Vw, Va, &M, p.patm);
        worst = std::max(worst, std::fabs(M - M0) / M0);
        for (auto &e : an.AirEq()) pamax = std::max(pamax, an.X()(e.second, 0));
    }
    REAL Vw, Va, M;
    pr.StoredVolumes(an.X(), p.hydro, p.n, Vw, Va, &M, p.patm);
    std::cout << "Boyle: V_a/V_a0 = " << Va / Va0 << ", p_a max = " << pamax << " kPa, max |dM|/M0 = " << worst << "\n";
    CHECK(Va / Va0 < 0.75);
    CHECK(pamax > 30.);
    CHECK(worst < 1.e-4);
}

TEST_CASE("Teste 3b: limite de ar passivo quando a mobilidade do ar cresce", "[ar]") {
    REAL tolp = 1.e-6;
    auto run = [&tolp](REAL factor, bool passive) {
        UnsatParams p = Soil();
        p.muRatio *= factor;
        // o ar passivo (p_a = 0) supoe ar sem peso: com gamma_a > 0 e p_a = 0 no topo e na base, p_a interno fica da
        // ordem de gamma_a H para qualquer mobilidade e o limite nao e zero
        p.ga = 0.;
        UnsatProblem pr(RectangleMesh(0.1, 1., 1, 20), {{1, p}}, kFix);
        auto &an = pr.A();
        an.SetTolerances(1.e-7, tolp, 40, 25);
        TPZFMatrix<STATE> &X = an.X();
        for (auto &e : an.WaterEq()) X(e.second, 0) = -25. - p.gw * pr.Y(e.first);
        an.Restore();
        an.Accept();
        auto fixed = an.AllDisplacements(an.X());
        if (passive) for (auto &e : an.AirEq()) fixed[e.second] = 0.;
        else for (int id : {ETop, EBase}) for (auto &n : an.BoundaryNodes(id)) fixed[an.AirEq().at(n.first)] = 0.;
        std::map<int64_t, REAL> fq;
        for (auto &n : an.BoundaryNodes(ETop)) fq[an.WaterEq().at(n.first)] = 0.5 * p.ks * n.second;
        for (int s = 0; s < 12; ++s) {
            if (!an.Solve(600., fixed, fq, false).conv) return std::map<int64_t, REAL>();
            an.Accept();
        }
        std::map<int64_t, REAL> pw;
        for (auto &e : an.WaterEq()) pw[e.first] = an.X()(e.second, 0);
        return pw;
    };
    if (const char *e = std::getenv("TOLP")) tolp = std::atof(e);
    const auto pas = run(1., true);
    REQUIRE(!pas.empty());
    REAL lo = 1.e30, hi = -1.e30;
    for (auto &v : pas) lo = std::min(lo, v.second), hi = std::max(hi, v.second);
    std::vector<REAL> diff;
    for (REAL f : {1., 10., 100.}) {
        const auto act = run(f, false);
        REQUIRE(!act.empty());
        REAL d = 0.;
        for (auto &v : act) d = std::max(d, std::fabs(v.second - pas.at(v.first)));
        diff.push_back(d);
        std::cout << "mu_w/mu_a x " << f << ": max |p_w(ativo) - p_w(passivo)| = " << d << " kPa (variacao de p_w "
                  << hi - lo << " kPa)\n";
    }
    CHECK(diff[1] < diff[0]);
    CHECK(diff[2] < diff[1]);
    CHECK(diff[2] < 0.01 * (hi - lo));
}

TEST_CASE("Teste 3c: balanco de ar e de agua com ar livre no topo", "[ar]") {
    UnsatParams p = Soil();
    p.Kw = 1.e15;
    UnsatProblem pr(RectangleMesh(0.1, 1., 1, 20), {{1, p}}, kFix);
    auto &an = pr.A();
    TPZFMatrix<STATE> &X = an.X();
    for (auto &e : an.WaterEq()) X(e.second, 0) = -20. - p.gw * pr.Y(e.first);
    an.Restore();
    an.Accept();
    auto fixed = an.AllDisplacements(an.X());
    const auto top = an.BoundaryNodes(ETop);
    for (auto &n : top) fixed[an.AirEq().at(n.first)] = 0.;
    std::map<int64_t, REAL> fq;
    REAL Lt = 0.;
    for (auto &n : top) fq[an.WaterEq().at(n.first)] = 0.8 * p.ks * n.second, Lt += n.second;
    REAL Vw0, Va0, M0;
    pr.StoredVolumes(an.X(), p.hydro, p.n, Vw0, Va0, &M0, p.patm);
    REAL airOut = 0., waterIn = 0.;
    for (int s = 0; s < 24; ++s) {
        const REAL dt = 600.;
        REQUIRE(an.Solve(dt, fixed, fq, false).conv);
        const TPZFMatrix<STATE> r = an.FullResidual();
        for (auto &n : top) airOut -= r.GetVal(an.AirEq().at(n.first), 0) * p.patm; // massa (x R T) que sai
        waterIn += 0.8 * p.ks * Lt * dt;
        an.Accept();
    }
    REAL Vw, Va, M;
    pr.StoredVolumes(an.X(), p.hydro, p.n, Vw, Va, &M, p.patm);
    const REAL eAir = std::fabs((M0 - M) - airOut) / M0, eW = std::fabs((Vw - Vw0) - waterIn) / waterIn;
    std::cout << "balanco: ar que saiu " << airOut / M0 << " M0, erro de ar " << eAir << " M0; erro de agua " << eW
              << " da entrada\n";
    CHECK(airOut > 0.);
    CHECK(eW < 1.e-6);
    CHECK(eAir < 1.e-3);
}
