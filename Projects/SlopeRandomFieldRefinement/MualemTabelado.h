// MualemTabelado.h - k_r de Mualem (van Genuchten, m = 1 - 1/n) tabelado contra S_w para o ETabulated do
// TPZUnsaturatedHydraulics, com um teto em S_e = 1 - eps (opcoes vg_tab= e vg_tab_psi= do SlopeRainMonteCarlo).
// Por que: com n < 2, dk_r/dpsi ~ psi^(n - 2) explode perto da saturacao, e o EMualem (que limita S_e a 1 - 1e-10) da
// k_r(1) = 0,984 com a inclinacao enorme junto de psi = 0. A tabela segue Mualem ate psi_c (S_e = 1 - eps) e fecha ate
// (S_w = 1, k_r = 1) de um destes jeitos (t = S_w - S_c, L_c = ln k_r(psi_c)):
//   reta:    ln k_r linear em S_w de (S_c, L_c) a (1, 0);
//   c1:      ln k_r = L_c + g_c t + A t^2 em 50 pontos, com g_c = d ln k_r / dS_w de Mualem em S_c (derivada continua em
//            S_c) e A = (-L_c - g_c eps') / eps'^2, eps' = 1 - S_c (chega em (1, 0));
//   ippisch: k_r / k_c (Ippisch et al. 2006, com a entrada de ar em psi_c): L_c sai de todos os pontos e o fecho e
//            (1, 0), isto e, k_r = 1 de S_c ate a saturacao.
// Acima do ultimo S_w o ETabulated da 1 com derivada nula, e Gamma = 0 em p_c <= 0: dlw = lw0 k g Gamma e continuo na
// saturacao. Os nos sao N + 1 valores de x = a psi logaritmicos de 1e6 a x_c, com
//   y = x^n, S_e = (1 + y)^(-m), S_w = S_wr + (1 - S_wr) S_e, ln k_r = -(m/2) log1p(y) + 2 log(-expm1(-m log1p(1/y)))
// (sem cancelamento perto da saturacao): a tabela so depende de n, S_wr e eps; a so muda psi_c. Nao usar o fecho do
// kfredlund (montecarlo.cpp, TabelaFredlund): ele deixa k_r < 1 na saturacao
#pragma once

#include "Plasticity/TPZUnsaturatedHydraulics.h"

#include <algorithm>
#include <cmath>
#include <vector>

enum class EFechoVG { Reta, C1, Ippisch };

inline const char *NomeFechoVG(EFechoVG f) { return f == EFechoVG::Reta ? "reta" : f == EFechoVG::C1 ? "c1" : "ippisch"; }

struct TabelaVG {
    std::vector<REAL> S, lnK; // S_w (estritamente crescente) e ln k_r: o fTabS e o fTabLnK do ETabulated
    // eps = 1 - S_e(psi_c); psi_c (kPa); kC: k_r de Mualem em psi_c (no ippisch, o divisor); erroMax: maior |k_tab / k - 1|
    // em [psi_c, 1e3 kPa] (2000 succoes logaritmicas), contra Mualem (reta e c1) ou Mualem / k_c (ippisch)
    REAL eps = 0., psiC = 0., kC = 0., erroMax = 0.;
};

namespace mualemtab {
// S_e, ln k_r e d ln k_r / dS_e de Mualem em y = (a psi)^n
inline REAL Se(REAL y, REAL m) { return std::exp(-m * std::log1p(y)); }
inline REAL LnK(REAL y, REAL m) { return -0.5 * m * std::log1p(y) + 2. * std::log(-std::expm1(-m * std::log1p(1. / y))); }
inline REAL DLnKdSe(REAL y, REAL m)
{
    // ln k_r = ln(S_e) / 2 + 2 ln f, f = 1 - (1 - q)^m, q = S_e^(1/m) = 1 / (1 + y): df/dS_e = (1 - q)^(m - 1) q / S_e
    const REAL q = 1. / (1. + y), f = -std::expm1(-m * std::log1p(1. / y));
    const REAL umqm1 = std::exp(-(m - 1.) * std::log1p(1. / y)); // (1 - q)^(m - 1), 1 - q = y / (1 + y)
    return (0.5 + 2. * umqm1 * q / f) / Se(y, m);
}
// k_r da tabela em S_w: a mesma interpolacao do ETabulated (ln k linear entre os pontos, constante fora deles)
inline REAL Interpola(const TabelaVG &t, REAL S)
{
    const size_t n = t.S.size();
    if (n < 2 || S <= t.S[0] || S >= t.S[n - 1]) return std::exp(n == 0 ? 0. : S <= t.S[0] ? t.lnK[0] : t.lnK[n - 1]);
    const size_t i = std::upper_bound(t.S.begin(), t.S.end(), S) - t.S.begin() - 1;
    const REAL g = (t.lnK[i + 1] - t.lnK[i]) / (t.S[i + 1] - t.S[i]);
    return std::exp(t.lnK[i] + g * (S - t.S[i]));
}
} // namespace mualemtab

// 1 - S_e(psi_c) da van Genuchten (a em 1/kPa, psi_c em kPa): o eps de vg_tab_psi=
inline REAL EpsDePsi(REAL a, REAL n, REAL psiC)
{
    const REAL m = 1. - 1. / n;
    return -std::expm1(-m * std::log1p(std::pow(a * psiC, n)));
}

// A tabela de ln k_r contra S_w (n > 1, 0 <= S_wr < 1, a em 1/kPa, 0 < eps < 1; N intervalos de x de 1e6 a x_c)
inline TabelaVG MualemTabelado(REAL n, REAL swr, REAL a, REAL eps, EFechoVG f, int N = 400)
{
    using namespace mualemtab;
    const REAL m = 1. - 1. / n;
    TabelaVG t;
    t.eps = eps;
    const REAL xc = std::pow(std::expm1(-std::log1p(-eps) / m), 1. / n); // ((1 - eps)^(-1/m) - 1)^(1/n)
    t.psiC = xc / a;
    const REAL l0 = std::log(1.e6), l1 = std::log(xc);
    for (int i = 0; i <= N; i++) { // do seco para o umido: S_w crescente (so os estritamente crescentes ficam)
        const REAL x = i == N ? xc : std::exp(l0 + (l1 - l0) * i / N), y = std::pow(x, n);
        const REAL S = swr + (1. - swr) * Se(y, m);
        if (!t.S.empty() && !(S > t.S.back())) continue;
        t.S.push_back(S), t.lnK.push_back(LnK(y, m));
    }
    const REAL Sc = t.S.back(), Lc = t.lnK.back();
    t.kC = std::exp(Lc);
    if (f == EFechoVG::C1) {
        const REAL gc = DLnKdSe(std::pow(xc, n), m) / (1. - swr); // d ln k_r / dS_w em S_c
        const REAL ep = 1. - Sc, A = (-Lc - gc * ep) / (ep * ep);
        for (int k = 1; k <= 50; k++) {
            const REAL s = ep * k / 50.;
            t.S.push_back(k == 50 ? 1. : Sc + s), t.lnK.push_back(k == 50 ? 0. : Lc + gc * s + A * s * s);
        }
    } else {
        if (f == EFechoVG::Ippisch)
            for (REAL &l : t.lnK) l -= Lc;
        t.S.push_back(1.), t.lnK.push_back(0.);
    }
    // erro da interpolacao em [psi_c, 1e3 kPa]
    const REAL p0 = std::log(t.psiC), p1 = std::log(std::max(REAL(1.e3), t.psiC));
    for (int i = 0; i < 2000; i++) {
        const REAL y = std::pow(a * std::exp(p0 + (p1 - p0) * i / 1999.), n);
        const REAL k = std::exp(LnK(y, m) - (f == EFechoVG::Ippisch ? Lc : 0.));
        t.erroMax = std::max(t.erroMax, std::abs(Interpola(t, swr + (1. - swr) * Se(y, m)) / k - 1.));
    }
    return t;
}

// O material passa a usar a tabela: fKrw = ETabulated, fTabS e fTabLnK (a retencao continua a van Genuchten)
inline void UsaTabela(TPZUnsaturatedHydraulics &h, const TabelaVG &t)
{
    h.fKrw = TPZUnsaturatedHydraulics::ETabulated;
    h.fTabS = t.S, h.fTabLnK = t.lnK;
}
