// Teste 1 do plano (SRM): talude seco (p_w = p_a = 0 congeladas) do projeto SlopeMohrCoulomb: altura 10 m, 45 graus,
// c = 10 kPa, phi = 30 graus, gamma = 20 kN/m3. Referencias: Bishop simplificado 1,206 (SRM) e o FS do modelo puramente
// mecanico, 1,203 no ciclo 5 de refinamento adaptativo (README de Projects/SlopeMohrCoulomb, branch claude/relaxed-volta-1r4gog).
#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "UnsatFS.h"
#include "UnsatProblem.h"

#include <chrono>
#include <iostream>

namespace {
UnsatParams DrySoil() {
    UnsatParams p;
    p.E = 2.e4, p.nu = 0.3, p.c = 10., p.phi = 30. * M_PI / 180., p.gamma = 20.;
    return p;
}
const std::map<int, std::array<bool, 2>> kFix = {{EBase, {true, true}}, {ELeft, {true, false}}, {ERight, {true, false}}};
const REAL kBishopSRM = 1.206;

REAL SRM(int nf, REAL F0, REAL tol, double *seconds = nullptr) {
    UnsatProblem pr(SlopeMeshFace(nf), {{1, DrySoil()}}, kFix);
    auto t0 = std::chrono::steady_clock::now();
    const FSResult r = FactorSRM(pr.A(), pr.A().X(), F0, tol);
    if (seconds) *seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    std::cout << "SRM: malha da face nf = " << nf << " (" << pr.A().NEq() << " equacoes), F0 = " << F0 << ": FS = " << r.fs
              << " (" << (seconds ? *seconds : 0.) << " s)\n";
    return r.fs;
}
} // namespace

TEST_CASE("Teste 1a: SRM do talude seco diminui com o refinamento e fica a 2% de Bishop", "[srm]") {
    double s;
    const REAL f20 = SRM(20, 0.6, 0.002, &s), f40 = SRM(40, 0.6, 0.002, &s);
    CHECK(f40 < f20);
    CHECK(std::fabs(f40 - kBishopSRM) <= 0.02 * kBishopSRM);
}

TEST_CASE("Teste 1a: SRM independente do estado de partida F0", "[srm]") {
    std::vector<REAL> fs;
    double s;
    for (REAL F0 : {0.5, 0.6, 0.8}) fs.push_back(SRM(20, F0, 0.002, &s));
    const REAL lo = *std::min_element(fs.begin(), fs.end()), hi = *std::max_element(fs.begin(), fs.end());
    CHECK(hi - lo <= 0.01);
}

TEST_CASE("Teste 1a (lento): malha nf = 80 a 1% de Bishop", "[srm][.lento]") {
    double s;
    CHECK(std::fabs(SRM(80, 0.6, 0.002, &s) - kBishopSRM) <= 0.01 * kBishopSRM);
}
