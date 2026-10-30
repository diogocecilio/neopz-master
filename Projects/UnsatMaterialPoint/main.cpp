// Ponto material: leis hidraulicas, pressao de Bishop (chi = S_w), Mohr-Coulomb PV2 do artigo (Tabelas 4-6,
// projecao de ponto mais proximo, Jacobiano) e tangente consistente (teste de Taylor, Eqs. 63-65 do artigo).
// Documento: docs/biot_nao_saturado/formulacao/formulacao.tex, Secoes 3, 5 e 10.
#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "UnsatModel.h"

#include <random>

using H = TPZUnsaturatedHydraulics;

namespace {
H Hydro(H::ERetention r, H::EWaterPermeability kw, H::EAirPermeability ka) {
    H h;
    h.fRetention = r, h.fKrw = kw, h.fKra = ka;
    h.fAVG = 0.05, h.fBVG = 2.3, h.fSwr = 0.1, h.fAG = 0.08, h.fA = 50., h.fN = 1.5, h.fM = 0.8;
    h.fLC = 3.79e-4, h.fLB = 2.4279, h.fPower = 4.;
    h.fKrwMin = 1.e-30, h.fKraMin = 1.e-30;
    return h;
}
} // namespace

TEST_CASE("Leis hidraulicas: Gamma = -dS/dp_c e dk_r/dS por diferencas centrais", "[ponto][hidraulica]") {
    for (auto r : {H::EFredlundXing, H::EVanGenuchten, H::ELiakopoulos, H::EGardner})
        for (auto kw : {H::EPowerLaw, H::EMualem, H::ELiakopoulosW})
            for (auto ka : {H::EBrooksCorey, H::EMualemAir}) {
                const H h = Hydro(r, kw, ka);
                for (REAL pc : {0.5, 2., 5., 12., 30.}) {
                    REAL S, G, Sp, Gp, Sm, Gm;
                    h.Retention(pc, S, G);
                    if (r == H::ELiakopoulos && S <= h.fSwr + 1.e-9) continue;
                    const REAL dp = 1.e-6 * pc;
                    h.Retention(pc + dp, Sp, Gp);
                    h.Retention(pc - dp, Sm, Gm);
                    INFO("retencao " << r << " p_c = " << pc);
                    CHECK(G == Approx(-(Sp - Sm) / (2. * dp)).epsilon(1.e-6));
                    REAL kr, dkr, krp, dk1, krm, dk2;
                    const REAL dS = 1.e-7;
                    h.WaterPermeability(S, kr, dkr);
                    h.WaterPermeability(S + dS, krp, dk1);
                    h.WaterPermeability(S - dS, krm, dk2);
                    if (kw == H::ELiakopoulosW && kr < 0.) continue; // fora do intervalo de validade da lei
                    INFO("k_rw " << kw);
                    CHECK(dkr == Approx((krp - krm) / (2. * dS)).epsilon(1.e-5).margin(1.e-9));
                    h.AirPermeability(S, kr, dkr);
                    h.AirPermeability(S + dS, krp, dk1);
                    h.AirPermeability(S - dS, krm, dk2);
                    INFO("k_ra " << ka);
                    CHECK(dkr == Approx((krp - krm) / (2. * dS)).epsilon(1.e-5).margin(1.e-9));
                }
            }
}

TEST_CASE("Leis hidraulicas: van Genuchten e Mualem iguais as notas (Secao 3 e 5)", "[ponto][hidraulica]") {
    const H h = Hydro(H::EVanGenuchten, H::EMualem, H::EMualemAir);
    const REAL pc = 7.3, m = 1. - 1. / h.fBVG;
    const REAL Se = std::pow(1. + std::pow(h.fAVG * pc, h.fBVG), -m);
    REAL S, G;
    h.Retention(pc, S, G);
    CHECK(S == Approx(h.fSwr + (1. - h.fSwr) * Se).epsilon(1.e-14));
    REAL kr, dkr;
    h.WaterPermeability(S, kr, dkr);
    CHECK(kr == Approx(std::sqrt(Se) * std::pow(1. - std::pow(1. - std::pow(Se, 1. / m), m), 2)).epsilon(1.e-12));
    h.AirPermeability(S, kr, dkr);
    CHECK(kr == Approx(std::sqrt(1. - Se) * std::pow(1. - std::pow(Se, 1. / m), 2. * m)).epsilon(1.e-12));
    REAL Ssat;
    h.Retention(-3., Ssat, G);
    CHECK(Ssat == 1.);
    CHECK(G == 0.);
}

TEST_CASE("Bishop com chi = S_w (padrao): p_B, chi_w = S_w - Gamma p_c e chi_a = S_a + Gamma p_c", "[ponto][bishop]") {
    UnsatParams p;
    TUnsat mat(1); // padrao do material: chi = S_w
    mat.SetHydraulics(p.hydro);
    for (auto pr : {std::pair<REAL, REAL>{-12., 3.}, {-40., 0.}, {5., 2.}, {-1., -0.5}}) {
        const REAL pw = pr.first, pa = pr.second;
        auto ev = [&](REAL w, REAL a, REAL &pB, REAL &cw, REAL &ca) {
            const auto hs = p.hydro.Evaluate(w, a);
            mat.Bishop(w, a, hs.S, hs.Gamma, pB, cw, ca);
            return hs;
        };
        REAL pB, cw, ca, p1, p2, d1, d2;
        const auto hs = ev(pw, pa, pB, cw, ca);
        CHECK(pB == Approx(hs.S * pw + (1. - hs.S) * pa).epsilon(1.e-14));
        CHECK(cw + ca == Approx(1.).epsilon(1.e-14));
        CHECK(cw == Approx(hs.S - hs.Gamma * (pa - pw)).epsilon(1.e-12));
        const REAL h = 1.e-6;
        ev(pw + h, pa, p1, d1, d2);
        ev(pw - h, pa, p2, d1, d2);
        CHECK(cw == Approx((p1 - p2) / (2. * h)).epsilon(1.e-6));
        ev(pw, pa + h, p1, d1, d2);
        ev(pw, pa - h, p2, d1, d2);
        CHECK(ca == Approx((p1 - p2) / (2. * h)).epsilon(1.e-6).margin(1.e-9));
    }
}

// ------------------------------------------------------------------ Mohr-Coulomb PV2 x artigo
namespace {
/// Coeficientes A_ij e denominador das Tabelas 4-6 do artigo (kind 0 face principal, 1 aresta direita, 2 esquerda)
void ArticleA(int kind, REAL G, REAL K, REAL c, REAL phi, REAL A[3][4], REAL &D) {
    const REAL s = std::sin(phi), cs = std::cos(phi);
    if (kind == 0) {
        D = 2. * (3. * G + (G + 3. * K) * s * s);
        const REAL t[3][4] = {{-3. * G * (s - 1.), (G + 3. * K) * (s - 1.) * s, 0., 2. * c * cs * (3. * K * s + G * (3. + s))},
                              {G * (s * s - 3.) - 6. * K * s * s, -3. * K * s * (1. + s) - G * (3. + (s - 2.) * s),
                               6. * G + 2. * (G + 3. * K) * s * s, (6. * K - 4. * G) * c * cs * s},
                              {3. * G * (1. + s), -(G + 3. * K) * s * (1. + s), 0., 2. * c * cs * (G * (s - 3.) + 3. * K * s)}};
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 4; ++j) A[i][j] = t[i][j];
    } else if (kind == 1) {
        D = 9. * G + 6. * G * s + G * s * s + 12. * K * s * s;
        const REAL a = G * (3. + s) - 6. * K * s, b = G * (3. + s) + 12. * K * s, d = G * (3. + s) + 3. * K * s;
        const REAL r1[4] = {-0.5 * (s - 1.) * a, 0.5 * (s - 1.) * b, -(s - 1.) * d, 4. * c * cs * d};
        const REAL r2[4] = {0.5 * (1. + s) * a, -0.5 * (1. + s) * b, (1. + s) * d, -2. * c * cs * a};
        for (int j = 0; j < 4; ++j) A[0][j] = r1[j], A[1][j] = r2[j], A[2][j] = r2[j];
    } else {
        D = 9. * G - 6. * G * s + G * s * s + 12. * K * s * s;
        const REAL a = G * (s - 3.) - 6. * K * s, d = G * (s - 3.) + 3. * K * s;
        const REAL r1[4] = {0.5 * (s - 1.) * a, -0.5 * (s - 1.) * a, (s - 1.) * d, -2. * c * cs * a};
        const REAL r3[4] = {-0.5 * (1. + s) * a, 0.5 * (1. + s) * a, -(1. + s) * d, 4. * c * cs * d};
        for (int j = 0; j < 4; ++j) A[0][j] = r1[j], A[1][j] = r1[j], A[2][j] = r3[j];
    }
}

/// Projecao de ponto mais proximo na norma de energia por enumeracao de conjuntos ativos dos seis planos (KKT)
std::array<REAL, 3> ClosestPoint(REAL G, REAL K, REAL c, REAL phi, const std::array<REAL, 3> &s) {
    const REAL sn = std::sin(phi), r = 2. * c * std::cos(phi);
    const int pairs[6][2] = {{0, 2}, {1, 2}, {1, 0}, {2, 0}, {2, 1}, {0, 1}};
    REAL N[6][3] = {};
    for (int k = 0; k < 6; ++k) N[k][pairs[k][0]] += 1. + sn, N[k][pairs[k][1]] += sn - 1.;
    const REAL a = K + 4. * G / 3., b = K - 2. * G / 3.;
    auto Dmul = [&](const REAL *v, REAL *o) {
        for (int i = 0; i < 3; ++i) o[i] = b * (v[0] + v[1] + v[2]) + (a - b) * v[i];
    };
    auto dot = [](const REAL *u, const REAL *v) { return u[0] * v[0] + u[1] * v[1] + u[2] * v[2]; };
    const REAL scale = 1. + std::max({std::fabs(s[0]), std::fabs(s[1]), std::fabs(s[2]), c});
    bool feas = true;
    for (auto &n : N) feas = feas && dot(n, s.data()) - r <= 1.e-12 * scale;
    if (feas) return s;
    for (int k = 0; k < 6; ++k) { // one plane
        REAL Dn[3];
        Dmul(N[k], Dn);
        const REAL lam = (dot(N[k], s.data()) - r) / dot(N[k], Dn);
        std::array<REAL, 3> sp = {s[0] - lam * Dn[0], s[1] - lam * Dn[1], s[2] - lam * Dn[2]};
        bool ok = lam >= -1.e-12 * scale;
        for (auto &n : N) ok = ok && dot(n, sp.data()) - r <= 1.e-10 * scale;
        if (ok) return sp;
    }
    for (int k = 0; k < 6; ++k) // two planes
        for (int l = k + 1; l < 6; ++l) {
            REAL Dk[3], Dl[3];
            Dmul(N[k], Dk), Dmul(N[l], Dl);
            const REAL m11 = dot(N[k], Dk), m12 = dot(N[k], Dl), m22 = dot(N[l], Dl);
            const REAL det = m11 * m22 - m12 * m12;
            if (std::fabs(det) < 1.e-12 * m11 * m22) continue;
            const REAL r1 = dot(N[k], s.data()) - r, r2 = dot(N[l], s.data()) - r;
            const REAL l1 = (m22 * r1 - m12 * r2) / det, l2 = (m11 * r2 - m12 * r1) / det;
            std::array<REAL, 3> sp;
            for (int i = 0; i < 3; ++i) sp[i] = s[i] - l1 * Dk[i] - l2 * Dl[i];
            bool ok = l1 >= -1.e-12 * scale && l2 >= -1.e-12 * scale;
            for (auto &n : N) ok = ok && dot(n, sp.data()) - r <= 1.e-10 * scale;
            if (ok) return sp;
        }
    const REAL apex = c / std::tan(phi);
    return {apex, apex, apex};
}
} // namespace

TEST_CASE("PV2: expressoes = Tabelas 4-6 e Eqs. 57-58 do artigo; regiao = ponto mais proximo", "[ponto][mc]") {
    std::mt19937 gen(2026);
    std::uniform_real_distribution<REAL> U(0., 1.);
    std::normal_distribution<REAL> N01(0., 1.);
    int count[5] = {0, 0, 0, 0, 0};
    REAL eart = 0., ejac = 0., ecpp = 0.;
    for (int t = 0; t < 20000; ++t) {
        const REAL E = std::pow(10., 3. + 4. * U(gen)), nu = 0.49 * U(gen);
        const REAL phi = (5. + 40. * U(gen)) * M_PI / 180., c = std::pow(10., 3. * U(gen));
        TPZElasticResponse ER;
        ER.SetEngineeringData(E, nu);
        const REAL G = ER.G(), K = ER.K();
        TPZYCMohrCoulombPV2 yc(phi, phi, c, ER);
        const REAL sc = c * (0.5 + 19.5 * U(gen)), shift = 5. * c * N01(gen);
        std::array<REAL, 3> s = {N01(gen) * sc + shift, N01(gen) * sc + shift, N01(gen) * sc + shift};
        std::sort(s.begin(), s.end(), std::greater<REAL>());
        const REAL ref = std::max({std::fabs(s[0]), std::fabs(s[2]), c});
        // (a) formulas of the three candidates, regardless of admissibility
        for (int kind = 0; kind < 3; ++kind) {
            TPZManVector<STATE, 3> st = {s[0], s[1], s[2]}, sp(3), eps(3);
            TPZManVector<STATE, 2> dl;
            TPZFNMatrix<9> Gr(3, 3, 0.);
            STATE al = 0., al1 = 0.;
            if (kind == 0) yc.ComputeLambdaSigmaMainPlane(st, al, dl, sp, eps, Gr, al1);
            else if (kind == 1) yc.ComputeLambdaSigmaRigth(st, al, dl, sp, eps, Gr, al1);
            else yc.ComputeLambdaSigmaLeft(st, al, dl, sp, eps, Gr, al1);
            REAL A[3][4], D;
            ArticleA(kind, G, K, c, phi, A, D);
            for (int i = 0; i < 3; ++i) {
                const REAL sa = ((s[0] + s[2]) * A[i][0] + (s[0] - s[2]) * A[i][1] + (s[0] + s[1]) * A[i][2] + A[i][3]) / D;
                eart = std::max(eart, std::fabs(sa - sp[i]) / ref);
                const REAL J[3] = {(A[i][0] + A[i][1] + A[i][2]) / D, A[i][2] / D, (A[i][0] - A[i][1]) / D};
                for (int j = 0; j < 3; ++j) ejac = std::max(ejac, std::fabs(J[j] - Gr(i, j)));
            }
        }
        // (b) region chosen by the code = closest point in the energy norm
        TPZManVector<STATE, 3> st = {s[0], s[1], s[2]}, sp(3), eps(3);
        TPZManVector<STATE, 2> dl;
        TPZFNMatrix<9> Gr(3, 3, 0.);
        STATE al = 0., al1 = 0.;
        int type = 0;
        yc.ProjectSigma(st, al, dl, sp, eps, Gr, al1, type);
        count[type == -1 ? 4 : type]++;
        const auto cp = ClosestPoint(G, K, c, phi, s);
        for (int i = 0; i < 3; ++i) ecpp = std::max(ecpp, std::fabs(cp[i] - sp[i]) / ref);
    }
    INFO("elastico " << count[0] << " face " << count[1] << " direita " << count[2] << " esquerda " << count[3]
                     << " apice " << count[4]);
    CHECK(eart < 1.e-12);
    CHECK(ejac < 1.e-12);
    CHECK(ecpp < 1.e-12);
    for (int k = 0; k < 5; ++k) CHECK(count[k] > 1000);
}

TEST_CASE("PV2: tangente consistente com ordem 2 no teste de Taylor; apice e elastico exatos", "[ponto][mc]") {
    UnsatParams p;
    p.E = 2.e4, p.nu = 0.49, p.c = 10., p.phi = 30. * M_PI / 180.;
    const TMC model = MohrCoulomb(p);
    std::mt19937 gen(7);
    std::uniform_real_distribution<REAL> U(-1., 1.);
    const REAL cases[][3] = {{4.e-3, 1.e-3, -5.e-3}, {5.e-3, -2.5e-3, -2.5e-3}, {2.5e-3, 2.5e-3, -5.e-3},
                             {1.e-3, 1.e-3, 1.e-3},  {3.e-3, -1.e-3, -2.e-3},   {2.e-4, 1.e-4, -3.e-4}};
    const int expected[] = {1, 2, 3, -1, 1, 0};
    for (int test = 0; test < 6; ++test) {
        TPZFNMatrix<9, REAL> Q(3, 3), A(3, 3); // rotacao aleatoria
        for (int i = 0; i < 9; ++i) A(i / 3, i % 3) = U(gen);
        for (int i = 0; i < 3; ++i) {
            for (int k = 0; k < i; ++k) {
                REAL d = 0.;
                for (int j = 0; j < 3; ++j) d += A(j, i) * Q(j, k);
                for (int j = 0; j < 3; ++j) A(j, i) -= d * Q(j, k);
            }
            REAL n = 0.;
            for (int j = 0; j < 3; ++j) n += A(j, i) * A(j, i);
            for (int j = 0; j < 3; ++j) Q(j, i) = A(j, i) / std::sqrt(n);
        }
        const REAL shift = test == 3 ? 0. : -2.e-5;
        TPZTensor<REAL> eps, deps;
        const int ij[6][2] = {{0, 0}, {0, 1}, {0, 2}, {1, 1}, {1, 2}, {2, 2}};
        for (int v = 0; v < 6; ++v) {
            REAL e = 0.;
            for (int k = 0; k < 3; ++k) e += Q(ij[v][0], k) * (cases[test][k] + shift) * Q(ij[v][1], k);
            eps[v] = (ij[v][0] == ij[v][1]) ? e : 2. * e; // distorcao de engenharia
            deps[v] = 1.e-3 * U(gen);
        }
        auto stress = [&](const TPZTensor<REAL> &e, TPZFMatrix<REAL> *D) {
            TMC m(model);
            m.SetState(TPZPlasticState<REAL>());
            TPZTensor<REAL> sig;
            m.ApplyStrainComputeSigma(e, sig, D);
            return std::make_pair(sig, m.GetState().m_m_type);
        };
        TPZFNMatrix<36, REAL> D(6, 6, 0.);
        auto [s0, type] = stress(eps, &D);
        INFO("caso " << test);
        CHECK(type == expected[test]);
        REAL err[2], alpha[2] = {1.e-3, 2.e-3};
        for (int k = 0; k < 2; ++k) {
            TPZTensor<REAL> e(eps);
            e.Add(deps, alpha[k]);
            TPZTensor<REAL> s1 = stress(e, nullptr).first;
            REAL n2 = 0.;
            for (int i = 0; i < 6; ++i) {
                REAL lin = 0.;
                for (int j = 0; j < 6; ++j) lin += D(i, j) * deps[j] * alpha[k];
                n2 += std::pow(s1[i] - s0[i] - lin, 2);
            }
            err[k] = std::sqrt(n2);
        }
        REAL asym = 0., dmax = 0.;
        for (int i = 0; i < 6; ++i)
            for (int j = 0; j < 6; ++j)
                asym = std::max(asym, std::fabs(D(i, j) - D(j, i))), dmax = std::max(dmax, std::fabs(D(i, j)));
        if (dmax > 0.) CHECK(asym / dmax < 1.e-12);
        if (expected[test] == 0 || expected[test] == -1) CHECK(err[0] < 1.e-10 * (Norm(s0) + 1.));
        else {
            const REAL order = std::log(err[1] / err[0]) / std::log(alpha[1] / alpha[0]);
            INFO("ordem " << order);
            CHECK(order > 1.85);
            CHECK(order < 2.15);
        }
    }
}

TEST_CASE("Reducao de resistencia: c/F, atan(tan phi/F), regra associada mantida", "[ponto][mc]") {
    UnsatParams p;
    TMC model = MohrCoulomb(p);
    for (REAL F : {0.6, 1., 1.37, 2.5}) {
        model.SetStrengthReductionFactor(F);
        const TPZYCMohrCoulombPV2 yc = model.LocalCriterion();
        CHECK(yc.Cohesion() == Approx(p.c / F).epsilon(1.e-14));
        CHECK(yc.Phi() == Approx(std::atan(std::tan(p.phi) / F)).epsilon(1.e-14));
        TPZPlasticState<REAL> st = const_cast<TPZYCMohrCoulombPV2 &>(yc).GetLocalMatState();
        CHECK(st.fmatprop[2] == Approx(st.fmatprop[1]).epsilon(1.e-14)); // psi = phi
    }
}
