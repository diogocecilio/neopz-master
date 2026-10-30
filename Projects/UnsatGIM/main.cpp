// Teste 1 do plano (GIM): talude seco do projeto SlopeMohrCoulomb. Referencias: Bishop simplificado 1,806 (aumento de
// gravidade) e 1,793 do modelo puramente mecanico no ciclo 5 (README de Projects/SlopeMohrCoulomb, branch
// claude/relaxed-volta-1r4gog). Tambem: GIM >= SRM e a succao (nao escalada) aumenta o fator.
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
const REAL kBishopGIM = 1.806;
} // namespace

namespace {
REAL GIM(int nf, REAL tol) {
    UnsatProblem pr(SlopeMeshFace(nf), {{1, DrySoil()}}, kFix);
    auto t0 = std::chrono::steady_clock::now();
    const REAL fs = FactorGIM(pr.A(), pr.A().X(), tol).fs;
    std::cout << "GIM: malha da face nf = " << nf << " (" << pr.A().NEq() << " equacoes): FS = " << fs << " ("
              << std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count() << " s)\n";
    return fs;
}
} // namespace

TEST_CASE("Teste 1b: GIM do talude seco diminui com o refinamento e fica a 5% de Bishop", "[gim]") {
    const REAL f20 = GIM(20, 0.002), f40 = GIM(40, 0.002);
    CHECK(f40 < f20);
    CHECK(std::fabs(f40 - kBishopGIM) <= 0.05 * kBishopGIM);
}

TEST_CASE("Teste 1b (lento): malha nf = 80 a 2% de Bishop", "[gim][.lento]") {
    CHECK(std::fabs(GIM(80, 0.002) - kBishopGIM) <= 0.02 * kBishopGIM);
}

TEST_CASE("Teste 1b: GIM >= SRM e succao nao escalada aumenta o fator", "[gim]") {
    UnsatParams p = DrySoil();
    UnsatProblem pr(SlopeMeshFace(20), {{1, p}}, kFix);
    auto &an = pr.A();
    const REAL gim = FactorGIM(an, an.X(), 0.005).fs;
    const REAL srm = FactorSRM(an, an.X(), 0.6, 0.005).fs;
    // succao uniforme de 20 kPa (p_w = -20, p_a = 0): p_B = S_w p_w < 0, nao escalada pelo fator de gravidade
    TPZFMatrix<STATE> x(an.X());
    for (auto &e : an.WaterEq()) x(e.second, 0) = -20.;
    const REAL gims = FactorGIM(an, x, 0.005).fs;
    std::cout << "nf = 20: GIM seco " << gim << ", SRM seco " << srm << ", GIM com succao de 20 kPa " << gims << "\n";
    CHECK(gim >= srm);
    CHECK(gims > gim + 0.05);
}
