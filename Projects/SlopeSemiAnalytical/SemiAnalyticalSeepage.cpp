#include "SemiAnalyticalSeepage.h"

#include "SlopeSolver.h" // Print

#include <algorithm>
#include <cmath>

SemiAnalyticalSeepage::SemiAnalyticalSeepage(const slope::Data &d, REAL Lm)
    : fGw(d.gammaw), fHw(d.hw), fH(d.H), fBeta(d.beta), fOmega(M_PI - d.beta), fRw(d.hw / std::sin(d.beta)),
      fR(d.H / std::sin(d.beta)), fRe(std::hypot(d.H, Lm + d.H / std::tan(d.beta))), fXO(d.xCrista), fYO(d.yCrista)
{
    // m em (0, 1): Sigma (1 - m) = omega (1 + m) / 2 (m = 0 e a raiz trivial), por bisseccao
    auto g = [this](REAL m) {
        const REAL a = std::sqrt(m);
        return std::sinh(2. * a * fOmega) / (4. * a) * (1. - m) - fOmega * (1. + m) / 2.;
    };
    REAL lo = 1.e-4, hi = 1. - 1.e-12;
    if (!(g(lo) > 0. && g(hi) < 0.)) DebugStop(); // sem raiz: face quase vertical (omega < sqrt(3))
    for (int i = 0; i < 200; i++) {
        const REAL mid = 0.5 * (lo + hi);
        (g(mid) > 0. ? lo : hi) = mid;
    }
    fM = 0.5 * (lo + hi);
    fA = std::sqrt(fM);
    const REAL Sigma = std::sinh(2. * fA * fOmega) / (4. * fA);
    fC = fM * (Sigma - fOmega / 2.); // int h2'^2 = m int sinh^2
    fD = Sigma + fOmega / 2.;        // int h2^2 = int cosh^2
    fK0 = fGw * std::sin(fBeta) * H2(fOmega) / (fD - fC);
}

REAL SemiAnalyticalSeepage::ThetaMax(REAL r) const { return r < fR ? fOmega : M_PI - std::asin(fH / r); }

std::array<REAL, 2> SemiAnalyticalSeepage::ForcePolar(REAL r, REAL t) const
{
    if (r >= fRe) return {0., 0.};
    if (r < fRw) {
        const REAL q = std::pow(r / fRw, fM - 1.);
        return {-fK0 * (1. - q) * DH2(t), fK0 * (1. - fM * q) * H2(t)};
    }
    return {0., fGw * fHw / (ThetaMax(r) * r)}; // arcos: a perda gamma_w h_w distribuida no comprimento do arco
}

std::array<REAL, 2> SemiAnalyticalSeepage::Force(const TPZVec<REAL> &x) const
{
    const REAL dx = x[0] - fXO, dy = x[1] - fYO, r = std::hypot(dx, dy);
    if (r < 1.e-12 || r >= fRe) return {0., 0.};
    const REAL t = std::min(std::max(std::atan2(-dy, -dx), REAL(0)), ThetaMax(r));
    const auto f = ForcePolar(r, t);
    const REAL c = std::cos(t), s = std::sin(t);
    return {-f[0] * c + f[1] * s, -f[0] * s - f[1] * c}; // f_r e_r + f_theta e_theta
}

REAL SemiAnalyticalSeepage::JStarClosed() const
{
    const int n = 4000; // int_R^Re 4 / (B r) dr, B = 2 theta_max (alpha = 1), Simpson em ln r
    const REAL a = std::log(fR), b = std::log(fRe), h = (b - a) / n;
    REAL I = 0.;
    for (int i = 0; i <= n; i++) {
        const REAL r = std::exp(a + i * h), w = (i == 0 || i == n) ? 1. : (i % 2 ? 4. : 2.);
        I += w * 2. / ThetaMax(r); // (4 / (B r)) dr = (2 / theta_max) d(ln r)
    }
    I *= h / 3.;
    // J*/k = -(gamma_w h_w)^2 / 4 [zona 1 + zona 2 + zona 3]; cada zona contribui com -1/2 int |f|^2 < 0, entao a zona 3
    // entra somada (no texto extraido do PDF da Eq. 40 ela aparece subtraida)
    const REAL t1 = H2(fOmega) / (std::sqrt(fC) + std::sqrt(fD));
    return -fGw * fGw * fHw * fHw / 4. * (t1 * t1 - 2. / fOmega * std::log(fHw / fH) + I);
}

REAL SemiAnalyticalSeepage::JStarNumeric() const
{
    auto simpson = [](int n, REAL a, REAL b, const std::function<REAL(REAL)> &f) {
        const REAL h = (b - a) / n;
        REAL s = 0.;
        for (int i = 0; i <= n; i++) s += ((i == 0 || i == n) ? 1. : (i % 2 ? 4. : 2.)) * f(a + i * h);
        return s * h / 3.;
    };
    // |f|^2 r integrado em theta de 0 a theta_max(r)
    auto anel = [&](REAL r) {
        return simpson(400, 0., ThetaMax(r), [&](REAL t) {
            const auto f = ForcePolar(r, t);
            return (f[0] * f[0] + f[1] * f[1]) * r;
        });
    };
    // zona 1 com r = R_w s^2 (a singularidade r^(2m-1) em O fica suave)
    const REAL z1 = simpson(2000, 0., 1., [&](REAL s) { return s > 0. ? anel(fRw * s * s) * 2. * fRw * s : 0.; });
    const REAL z2 = fRw < fR ? simpson(400, fRw, fR, anel) : 0.;
    const REAL z3 = simpson(2000, std::log(fR), std::log(fRe), [&](REAL u) { return anel(std::exp(u)) * std::exp(u); });
    return -0.5 * (z1 + z2 + z3);
}

void SemiAnalyticalSeepage::Print() const
{
    using slope::Print;
    Print("  campo semi-analitico (alpha = 1): O = (%.2f, %.2f), omega = %.4f rad", fXO, fYO, fOmega);
    Print("    R_w = %.4f m, R = %.4f m, R_e = %.4f m (L_m = 10 H)", fRw, fR, fRe);
    Print("    zona 1: m = sqrt(C/D) = %.6f, h2(omega) = %.6f, C = %.6f, D = %.6f, K0 = %.6f kN/m3", fM, H2(fOmega), fC,
          fD, fK0);
    const REAL jc = JStarClosed(), jn = JStarNumeric();
    Print("    J*/k: Eq. (40) %.6f, integrado %.6f (diferenca %.1e)", jc, jn, std::fabs(jc - jn) / std::fabs(jc));
    REAL div = 0.; // div f = (1/r) d(r f_r)/dr + (1/r) d f_theta / d theta: deve ser 0
    for (REAL s : {0.2, 0.5, 0.8})
        for (REAL t : {0.3, 1.2, 2.0}) {
            const REAL r = s * fRw, h = 1.e-6 * r, k = 1.e-6;
            const auto f = ForcePolar(r, t), fp = ForcePolar(r + h, t), fm = ForcePolar(r - h, t);
            const auto gp = ForcePolar(r, t + k), gm = ForcePolar(r, t - k);
            const REAL d = ((r + h) * fp[0] - (r - h) * fm[0]) / (2. * h * r) + (gp[1] - gm[1]) / (2. * k * r);
            div = std::max(div, std::fabs(d) * r / std::hypot(f[0], f[1]));
        }
    Print("    div f na zona 1 (relativo): %.1e", div);
}
