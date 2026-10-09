// Semi-analytical optimal filtration velocity v'_opt of Ceron, Cecilio, Linn & Maghous (IJNAMG 2025), Sect. 4.1
// (Eqs. 29-40), and the seepage force f = K^-1 . v'_opt as a slope::ForceField: C++ port of
// scripts/analytical_seepage.py (derivation, with the corrections of the printed Eqs. 31 and 39, in
// scripts/analytical_seepage_derivation.md).
//
// Polar coordinates of the paper's Fig. 3 around the crest edge O, written in the paper coordinates (x right, y DOWN):
// x = -r cos(theta), y = r sin(theta), theta = atan2(y, -x) from the crest (theta = 0) to the face
// (theta = Theta = pi - beta), e_r = (-cos, sin), e_theta = (sin, cos). Class of velocities (Eq. 29):
//   zone 1, r < R_w      : v = -h1(r) h2'(theta) e_r + (r h1)' h2(theta) e_theta, h1(R_w) = 0
//   zone 2, R_w <= r < R : v = h3 e_theta, h3 = k_h gamma_w h_w / (A r)                    (Eqs. 32, 34)
//   zone 3, R <= r < R_e : v = h4 e_theta, h4 = 2 k_h gamma_w h_w / (B(r) r)               (Eqs. 33, 35)
//   v = 0 for r >= R_e and outside the soil (and everywhere for h_w = 0)
// with R_w = h_w / sin(beta), R = H / sin(beta), R_e = sqrt(H^2 + (L_m + H / tan(beta))^2), L_m = 10 H by default.
// Zone 1 at the optimum (Eq. 31 corrected: prefactor 1 / (D - C), exponent and coefficient sqrt(C/D) in the e_theta
// term), with s = r / R_w, m = sqrt(C / D) and E(s) = (1 - s^(m-1)) / (1 - m) (-> ln s for m = 1):
//   v_r = -a' E h2'(theta),  v_theta = a' (E + s^(m-1)) h2(theta),  a' = k_h gamma_w sin(beta) h2(Theta) / (D (1 + m)),
// h2 the solution of (c h2')' = m d h2, h2(0) = 1, h2'(0) = 0 (Eqs. 37-38; the printed Eq. 39 only normalizes h2),
// c = cos^2 + alpha sin^2, d = sin^2 + alpha cos^2, C = int_0^Theta c h2'^2, D = int_0^Theta d h2^2 (Eq. 36). The
// optimal m minimizes Phi(m) = (1 + 1/m) c h2' / h2 (Theta), whose minimum is (sqrt C + sqrt D)^2 / h2(Theta)^2 at the
// fixed point m = sqrt(C/D). Steep slopes (Phi increasing from Phi(0+) = A; beta > 80.76 deg for alpha = 1, 81.10 deg
// for alpha = 2, 88.40 deg for alpha = 4, none for alpha >= 5): degenerate optimum, the limit m -> 0 of the class,
// h2 = 1 and the purely tangential zone-1 field v = k_h gamma_w sin(beta) / A e_theta.
//   J*(v'_opt) = -(k_h h_w^2 gamma_w^2 / 4) [h2(Theta)^2 / (sqrt C + sqrt D)^2 + (2 / A) ln(H / h_w)
//                                            + int_R^R_e 4 dr / (B r)]                        (Eq. 40)
// f = K^-1 v'_opt with K^-1 = diag(1 / k_h, alpha / k_h) in (x, y), alpha = k_h / k_v: f does not depend on k_h.
//
// Numerics: h2 by an adaptive Dormand-Prince 5(4) integrator (rtol 1e-13) of the scaled state (h2, w = c h2' / m,
// int w^2 / c, int d h2^2), well conditioned for m -> 0; m as in the Python (scan of Phi on 81 points of ln m in
// [ln 1e-7, ln 30], Brent minimization around the smallest value, then Brent root of the fixed point m = sqrt(C/D);
// the smallest value at m = 1e-7 means the degenerate optimum); h2 and h2' tabulated on 2049 uniform angles
// (integrator stopped at each one) and interpolated by cubic Hermite polynomials (with h2' and h2'' = (m d h2 -
// c' h2') / c from the ODE); the zone-3 integral of Eq. 40 by adaptive Gauss-Kronrod (7-15). Construction takes a
// few ms; evaluation is O(1) (40-80 ns per point) and the object is read-only after construction (thread safe).
// Coordinates: Force(x, f) takes NeoPZ coordinates (y up, O at the origin, y_paper = -y) and returns f with y up, as
// the other slope::ForceField; ForcePaper / VelocityPaper / PolarVelocity use the paper coordinates.
#ifndef ANALYTICALSEEPAGE_H
#define ANALYTICALSEEPAGE_H

#include "SeepageForceField.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace slope {

namespace numerics {

/// Adaptive Dormand-Prince 5(4) integration of y' = f(t, y, dy) (N components) from t0 to t1, y updated in place.
/// Local error control (RMS over the components) error <= atol + rtol max(|y_old|, |y_new|), as scipy's RK45.
/// h: in, the step to try first (<= 0: automatic); out, the step proposed by the controller (kept across output
/// points: a step shortened to land on t1 does not reduce it). Returns the number of steps, -1 on failure.
template <int N, class Rhs>
int DormandPrince54(const Rhs &f, REAL t0, REAL t1, REAL y[N], REAL rtol, REAL atol, REAL &h, int maxSteps = 1000000) {
    static constexpr REAL c2 = 1. / 5., c3 = 3. / 10., c4 = 4. / 5., c5 = 8. / 9.;
    static constexpr REAL a21 = 1. / 5., a31 = 3. / 40., a32 = 9. / 40., a41 = 44. / 45., a42 = -56. / 15., a43 = 32. / 9.,
                          a51 = 19372. / 6561., a52 = -25360. / 2187., a53 = 64448. / 6561., a54 = -212. / 729.,
                          a61 = 9017. / 3168., a62 = -355. / 33., a63 = 46732. / 5247., a64 = 49. / 176.,
                          a65 = -5103. / 18656.;
    static constexpr REAL b1 = 35. / 384., b3 = 500. / 1113., b4 = 125. / 192., b5 = -2187. / 6784., b6 = 11. / 84.;
    // error weights b - b* (5th minus embedded 4th order)
    static constexpr REAL e1 = 71. / 57600., e3 = -71. / 16695., e4 = 71. / 1920., e5 = -17253. / 339200.,
                          e6 = 22. / 525., e7 = -1. / 40.;
    const REAL span = t1 - t0;
    if (!(span > 0.)) return 0;
    if (!(h > 0.)) h = 1.e-3 * span;
    REAL k1[N], k2[N], k3[N], k4[N], k5[N], k6[N], k7[N], yt[N], yn[N];
    REAL t = t0;
    f(t, y, k1);
    int steps = 0;
    while (t < t1) {
        if (++steps > maxSteps) return -1;
        const bool last = t + h >= t1 - 1.e-12 * span;
        const REAL hs = last ? t1 - t : h;
        for (int i = 0; i < N; i++) yt[i] = y[i] + hs * a21 * k1[i];
        f(t + c2 * hs, yt, k2);
        for (int i = 0; i < N; i++) yt[i] = y[i] + hs * (a31 * k1[i] + a32 * k2[i]);
        f(t + c3 * hs, yt, k3);
        for (int i = 0; i < N; i++) yt[i] = y[i] + hs * (a41 * k1[i] + a42 * k2[i] + a43 * k3[i]);
        f(t + c4 * hs, yt, k4);
        for (int i = 0; i < N; i++) yt[i] = y[i] + hs * (a51 * k1[i] + a52 * k2[i] + a53 * k3[i] + a54 * k4[i]);
        f(t + c5 * hs, yt, k5);
        for (int i = 0; i < N; i++)
            yt[i] = y[i] + hs * (a61 * k1[i] + a62 * k2[i] + a63 * k3[i] + a64 * k4[i] + a65 * k5[i]);
        f(t + hs, yt, k6);
        for (int i = 0; i < N; i++) yn[i] = y[i] + hs * (b1 * k1[i] + b3 * k3[i] + b4 * k4[i] + b5 * k5[i] + b6 * k6[i]);
        f(t + hs, yn, k7);
        REAL err = 0.;
        for (int i = 0; i < N; i++) {
            const REAL sc = atol + rtol * std::max(std::fabs(y[i]), std::fabs(yn[i]));
            const REAL ei = hs * (e1 * k1[i] + e3 * k3[i] + e4 * k4[i] + e5 * k5[i] + e6 * k6[i] + e7 * k7[i]) / sc;
            err += ei * ei;
        }
        err = std::sqrt(err / N);
        if (!std::isfinite(err)) {
            h = 0.2 * hs;
            continue;
        }
        const REAL fac = err > 0. ? std::min(5., std::max(0.2, 0.9 * std::pow(err, -0.2))) : 5.;
        if (err <= 1.) {
            t = last ? t1 : t + hs;
            for (int i = 0; i < N; i++) y[i] = yn[i], k1[i] = k7[i]; // first same as last
            if (!last || hs >= h) h = hs * fac; // a step shortened to land on t1 keeps the previous proposal
        } else {
            h = hs * fac;
        }
    }
    return steps;
}

/// Adaptive Gauss-Kronrod (7-15 points, QUADPACK nodes) quadrature of f on [a, b]: the interval with the largest
/// error estimate |K15 - G7| is bisected until the sum of the estimates <= max(epsabs, epsrel |I|)
template <class Fun>
REAL AdaptiveGaussKronrod(const Fun &f, REAL a, REAL b, REAL epsabs, REAL epsrel, int maxIntervals = 2000) {
    static const REAL xgk[8] = {0.991455371120812639206854697526329, 0.949107912342758524526189684047851,
                                0.864864423359769072789712788640926, 0.741531185599394439863864773280788,
                                0.586087235467691130294144845693013, 0.405845151377397166906606412076961,
                                0.207784955007898467600689403773245, 0.};
    static const REAL wgk[8] = {0.022935322010529224963732008058970, 0.063092092629978553290700663189204,
                                0.104790010322250183839876322541518, 0.140653259715525918745189590510238,
                                0.169004726639267902826583426598550, 0.190350578064785409913256402421014,
                                0.204432940075298892414161999234649, 0.209482141084727828012999174891714};
    static const REAL wg[4] = {0.129484966168869693270611432679082, 0.279705391489276667901467771423780,
                               0.381830050505118944950369775488975, 0.417959183673469387755102040816327};
    struct Piece {
        REAL a, b, I, err;
    };
    auto rule = [&](REAL lo, REAL hi) {
        const REAL c = 0.5 * (lo + hi), hl = 0.5 * (hi - lo);
        const REAL fc = f(c);
        REAL K = fc * wgk[7], G = fc * wg[3];
        for (int j = 0; j < 7; j++) {
            const REAL s = f(c - hl * xgk[j]) + f(c + hl * xgk[j]);
            K += wgk[j] * s;
            if (j % 2 == 1) G += wg[j / 2] * s;
        }
        return Piece{lo, hi, K * hl, std::fabs((K - G) * hl)};
    };
    std::vector<Piece> pieces = {rule(a, b)};
    for (;;) {
        REAL I = 0., E = 0.;
        size_t worst = 0;
        for (size_t i = 0; i < pieces.size(); i++) {
            I += pieces[i].I, E += pieces[i].err;
            if (pieces[i].err > pieces[worst].err) worst = i;
        }
        if (E <= std::max(epsabs, epsrel * std::fabs(I)) || int(pieces.size()) >= maxIntervals) return I;
        const Piece w = pieces[worst];
        const REAL mid = 0.5 * (w.a + w.b);
        if (!(mid > w.a && mid < w.b)) return I; // interval at round-off size
        pieces[worst] = rule(w.a, mid);
        pieces.push_back(rule(mid, w.b));
    }
}

/// Bounded minimization of f on [x1, x2] by Brent's method (port of scipy.optimize.fminbound: golden section and
/// parabolic interpolation; stops when the bracket is within 2 (sqrt(eps) |x| + xatol / 3))
template <class Fun>
REAL BrentMinimize(const Fun &func, REAL x1, REAL x2, REAL xatol, int maxfun = 500) {
    const REAL sqrtEps = std::sqrt(2.2e-16), goldenMean = 0.5 * (3. - std::sqrt(5.));
    REAL a = x1, b = x2;
    REAL fulc = a + goldenMean * (b - a), nfc = fulc, xf = fulc;
    REAL rat = 0., e = 0., x = xf, fx = func(x);
    int num = 1;
    REAL ffulc = fx, fnfc = fx, xm = 0.5 * (a + b);
    REAL tol1 = sqrtEps * std::fabs(xf) + xatol / 3., tol2 = 2. * tol1;
    while (std::fabs(xf - xm) > tol2 - 0.5 * (b - a)) {
        bool golden = true;
        if (std::fabs(e) > tol1) { // parabolic fit
            golden = false;
            REAL r = (xf - nfc) * (fx - ffulc), q = (xf - fulc) * (fx - fnfc);
            REAL p = (xf - fulc) * q - (xf - nfc) * r;
            q = 2. * (q - r);
            if (q > 0.) p = -p;
            q = std::fabs(q);
            r = e, e = rat;
            if (std::fabs(p) < std::fabs(0.5 * q * r) && p > q * (a - xf) && p < q * (b - xf)) {
                rat = p / q;
                x = xf + rat;
                if (x - a < tol2 || b - x < tol2) rat = tol1 * ((xm - xf) >= 0. ? 1. : -1.);
            } else {
                golden = true;
            }
        }
        if (golden) {
            e = xf >= xm ? a - xf : b - xf;
            rat = goldenMean * e;
        }
        const REAL si = rat >= 0. ? 1. : -1.;
        x = xf + si * std::max(std::fabs(rat), tol1);
        const REAL fu = func(x);
        num++;
        if (fu <= fx) {
            if (x >= xf) a = xf;
            else b = xf;
            fulc = nfc, ffulc = fnfc;
            nfc = xf, fnfc = fx;
            xf = x, fx = fu;
        } else {
            if (x < xf) a = x;
            else b = x;
            if (fu <= fnfc || nfc == xf) {
                fulc = nfc, ffulc = fnfc;
                nfc = x, fnfc = fu;
            } else if (fu <= ffulc || fulc == xf || fulc == nfc) {
                fulc = x, ffulc = fu;
            }
        }
        xm = 0.5 * (a + b);
        tol1 = sqrtEps * std::fabs(xf) + xatol / 3.;
        tol2 = 2. * tol1;
        if (num >= maxfun) break;
    }
    return xf;
}

/// Root of f in [xa, xb] (f(xa) f(xb) < 0) by Brent's method (port of scipy.optimize.brentq), converged when the
/// bracket is below xtol + rtol |x|
template <class Fun>
REAL BrentRoot(const Fun &f, REAL xa, REAL xb, REAL xtol, REAL rtol, int maxiter = 200) {
    REAL xpre = xa, xcur = xb, xblk = 0., fpre = f(xpre), fcur = f(xcur), fblk = 0., spre = 0., scur = 0.;
    if (fpre == 0.) return xpre;
    if (fcur == 0.) return xcur;
    if (std::signbit(fpre) == std::signbit(fcur)) throw std::runtime_error("BrentRoot: no sign change");
    for (int i = 0; i < maxiter; i++) {
        if (fpre != 0. && fcur != 0. && std::signbit(fpre) != std::signbit(fcur)) {
            xblk = xpre, fblk = fpre;
            spre = scur = xcur - xpre;
        }
        if (std::fabs(fblk) < std::fabs(fcur)) {
            xpre = xcur, xcur = xblk, xblk = xpre;
            fpre = fcur, fcur = fblk, fblk = fpre;
        }
        const REAL delta = 0.5 * (xtol + rtol * std::fabs(xcur)), sbis = 0.5 * (xblk - xcur);
        if (fcur == 0. || std::fabs(sbis) < delta) return xcur;
        if (std::fabs(spre) > delta && std::fabs(fcur) < std::fabs(fpre)) {
            REAL stry;
            if (xpre == xblk) { // interpolate
                stry = -fcur * (xcur - xpre) / (fcur - fpre);
            } else { // extrapolate
                const REAL dpre = (fpre - fcur) / (xpre - xcur), dblk = (fblk - fcur) / (xblk - xcur);
                stry = -fcur * (fblk * dblk - fpre * dpre) / (dblk * dpre * (fblk - fpre));
            }
            if (2. * std::fabs(stry) < std::min(std::fabs(spre), 3. * std::fabs(sbis) - delta)) {
                spre = scur, scur = stry; // good short step
            } else {
                spre = sbis, scur = sbis; // bisect
            }
        } else {
            spre = sbis, scur = sbis;
        }
        xpre = xcur, fpre = fcur;
        xcur += std::fabs(scur) > delta ? scur : (sbis > 0. ? delta : -delta);
        fcur = f(xcur);
    }
    return xcur;
}

} // namespace numerics

class AnalyticalSeepage {
public:
    /// beta (deg) in (0, 90]; 0 <= h_w <= H (m); alpha = k_h / k_v >= 1; k_h only scales v and J*; L_m = LmOverH H
    AnalyticalSeepage(REAL betaDeg, REAL H, REAL hw, REAL alpha, REAL kh = 1., REAL gammaw = 9.81, REAL LmOverH = 10.)
        : fBetaDeg(betaDeg), fH(H), fHw(std::min(hw, H)), fAlpha(alpha), fKh(kh), fGammaw(gammaw), fLm(LmOverH * H) {
        const REAL beta = betaDeg * M_PI / 180.;
        if (!(beta > 0. && beta <= 0.5 * M_PI + 1.e-12)) throw std::invalid_argument("AnalyticalSeepage: beta must be in (0, 90] deg");
        if (!(hw >= 0. && hw <= H * (1. + 1.e-12))) throw std::invalid_argument("AnalyticalSeepage: hw must be in [0, H]");
        if (!(alpha >= 1.)) throw std::invalid_argument("AnalyticalSeepage: alpha = k_h / k_v must be >= 1");
        if (!(H > 0. && kh > 0. && LmOverH > 0.)) throw std::invalid_argument("AnalyticalSeepage: H, k_h, L_m must be > 0");
        fSb = std::sin(beta), fCb = std::cos(beta);
        fTanFinite = std::fabs(fCb) >= 1.e-15;
        fTb = fTanFinite ? fSb / fCb : 0.;
        fTheta = M_PI - beta;
        fXToe = fH * fCb / fSb;
        fR = fH / fSb;
        fRw = fHw / fSb;
        fRe = std::hypot(fH, fLm + fXToe);
        fA = (fAlpha + 1.) * fTheta / 2. - (fAlpha - 1.) * std::sin(2. * beta) / 4.; // Eq. 34
        // zone-3 integral of Eq. 40, int_R^Re 4 dr / (B r), with r = sqrt(H^2 + X^2) (X = abscissa on the toe ground):
        // dr / r = X dX / r^2, and B written with atan2 is smooth up to X = 0 (beta = 90 deg)
        const REAL a = fAlpha, HH = fH;
        auto i3 = [a, HH](REAL X) {
            const REAL r2 = HH * HH + X * X;
            const REAL B = (a + 1.) * (M_PI - std::atan2(HH, X)) - (a - 1.) * HH * X / r2;
            return 4. * X / (B * r2);
        };
        fI3 = fRe > fR ? numerics::AdaptiveGaussKronrod(i3, fXToe, fLm + fXToe, 1.e-15, 1.e-14) : 0.;
        if (fHw > 0.) {
            SolveH2();
        } else { // no drawdown: v = 0 (scalars as the Python: degenerate values)
            fM = 0., fC = 0., fD = fA, fH2e = 1., fPhi = fA, fDegenerate = true;
        }
        fF = fH2e * fH2e / std::pow(std::sqrt(fC) + std::sqrt(fD), 2);
        fAp = fHw > 0. ? fKh * fGammaw * fSb * fH2e / (fD * (1. + fM)) : 0.;
        fUnitM = std::fabs(fM - 1.) < 1.e-12;
    }

    // ---------------------------------------------------------------- data and results
    REAL BetaDeg() const { return fBetaDeg; }
    REAL H() const { return fH; }
    REAL Hw() const { return fHw; }
    REAL Alpha() const { return fAlpha; }
    REAL Kh() const { return fKh; }
    REAL Gammaw() const { return fGammaw; }
    REAL Lm() const { return fLm; }
    REAL Theta() const { return fTheta; } ///< opening of zones 1 and 2, pi - beta
    REAL XToe() const { return fXToe; }   ///< H / tan(beta)
    REAL Rw() const { return fRw; }
    REAL R() const { return fR; }
    REAL Re() const { return fRe; }
    REAL A() const { return fA; }   ///< Eq. 34
    REAL I3() const { return fI3; } ///< int_R^Re 4 dr / (B r)
    REAL M() const { return fM; }   ///< m = sqrt(C / D) (0: degenerate optimum)
    REAL C() const { return fC; }   ///< Eq. 36 with h2(0) = 1
    REAL D() const { return fD; }
    REAL H2e() const { return fH2e; } ///< h2(Theta)
    REAL Phi() const { return fPhi; } ///< min Phi = 1 / F
    REAL F() const { return fF; }     ///< h2(Theta)^2 / (sqrt C + sqrt D)^2, J1* = -k_h gamma_w^2 h_w^2 F / 4
    REAL APrime() const { return fAp; }
    bool Degenerate() const { return fDegenerate; }
    int ShootingCount() const { return fShootings; }

    /// Eq. 35: B(r) = 2 int_0^{pi - arcsin(H / r)} d(theta) dtheta, r >= H
    REAL B(REAL r) const {
        const REAL X = std::sqrt(std::max(r * r - fH * fH, 0.));
        return (fAlpha + 1.) * (M_PI - std::atan2(fH, X)) - (fAlpha - 1.) * fH * X / (r * r);
    }

    /// contributions of the zones 1, 2, 3 to J*(v'_opt) (Eq. 40)
    void JstarParts(REAL J[3]) const {
        J[0] = J[1] = J[2] = 0.;
        if (fHw <= 0.) return;
        const REAL pre = -fKh * fGammaw * fGammaw * fHw * fHw / 4.;
        J[0] = pre * fF, J[1] = pre * 2. / fA * std::log(fH / fHw) + 0., J[2] = pre * fI3; // + 0: no -0 for h_w = H
    }
    REAL Jstar() const {
        REAL J[3];
        JstarParts(J);
        return J[0] + J[1] + J[2];
    }
    /// J*(v'_opt) / (k_h H^2 gamma_w^2) (negative; the paper's Fig. 5 plots minus this, solid curves)
    REAL JstarNormalized() const { return Jstar() / (fKh * fH * fH * fGammaw * fGammaw); }

    // ---------------------------------------------------------------- h2 problem (independent of h_w)
    /// integrates (c h2')' = m d h2, h2(0) = 1, h2'(0) = 0 from 0 to thetaEnd (default Theta):
    /// y = {h2, w = c h2' / m, int w^2 / c, int d h2^2}
    void Shoot(REAL m, REAL y[4], REAL rtol = 1.e-13, REAL thetaEnd = -1.) const {
        y[0] = 1., y[1] = y[2] = y[3] = 0.;
        REAL h = 0.;
        if (numerics::DormandPrince54<4>(H2Rhs{fAlpha, m}, 0., thetaEnd < 0. ? fTheta : thetaEnd, y, rtol, 1.e-15 * fA, h) < 0)
            throw std::runtime_error("AnalyticalSeepage: h2 integration failed");
    }
    /// reduced objective Phi(m) = (1 + 1/m) c h2'(Theta) / h2(Theta) = (1 + m) w / h2 (Theta); Phi(0+) = A
    REAL PhiOfM(REAL m) const {
        if (m <= 0.) return fA;
        REAL y[4];
        Shoot(m, y);
        return (1. + m) * y[1] / y[0];
    }
    /// m - sqrt(C(m) / D(m)) (zero at the optimum), C = m^2 int w^2 / c
    REAL FixedPointResidual(REAL m) const {
        REAL y[4];
        Shoot(m, y);
        return m - m * std::sqrt(y[2] / y[3]);
    }
    /// h2(theta) and h2'(theta), theta clipped to [0, Theta]
    void H2(REAL theta, REAL &h2, REAL &dh2) const {
        if (fDegenerate) {
            h2 = 1., dh2 = 0.;
            return;
        }
        const REAL t = std::min(std::max(theta, 0.), fTheta) * fInvDT;
        const int n = int(fTab.size()) - 1;
        const int i = std::min(int(t), n - 1);
        const REAL u = t - i, u1 = 1. - u;
        const REAL H00 = (1. + 2. * u) * u1 * u1, H10 = u * u1 * u1, H01 = u * u * (3. - 2. * u), H11 = -u * u * u1;
        const H2Node &p = fTab[i], &q = fTab[i + 1];
        h2 = H00 * p.h + H01 * q.h + fDT * (H10 * p.dh + H11 * q.dh);
        dh2 = H00 * p.dh + H01 * q.dh + fDT * (H10 * p.d2h + H11 * q.d2h);
    }

    // ---------------------------------------------------------------- fields
    /// soil region in the paper coordinates: crest y = 0 (x <= 0), face, toe ground y = H
    bool InSoilPaper(REAL x, REAL yp) const {
        if (x <= 0.) return yp >= 0.;
        if (x < fXToe) return yp >= (fTanFinite ? x * fTb : 0.);
        return yp >= fH;
    }
    /// 0 outside the soil, 1, 2, 3 (r < R_w, R_w <= r < R, R <= r < R_e), 4 (r >= R_e); paper coordinates
    int Zone(REAL x, REAL yp) const {
        if (!InSoilPaper(x, yp)) return 0;
        const REAL r = std::sqrt(x * x + yp * yp);
        return r < fRw ? 1 : (r < fR ? 2 : (r < fRe ? 3 : 4));
    }
    /// (v_r, v_theta) of v'_opt at the polar point (r, theta) assumed in the soil (zero for r = 0 or r >= R_e)
    void PolarVelocity(REAL r, REAL theta, REAL &vr, REAL &vt) const {
        vr = vt = 0.;
        if (fHw <= 0. || !(r > 0.) || r >= fRe) return;
        if (r < fRw) {
            if (fDegenerate) { // m -> 0: E + s^(m-1) = 1, h2 = 1, h2' = 0
                vt = fAp;
                return;
            }
            REAL h2, dh2;
            H2(theta, h2, dh2);
            const REAL ls = std::log(r / fRw), em = (fM - 1.) * ls;
            const REAL E = fUnitM ? ls : std::expm1(em) / (fM - 1.);
            vr = -fAp * E * dh2;
            vt = fAp * (E + std::exp(em)) * h2;
        } else if (r < fR) {
            vt = fKh * fGammaw * fHw / (fA * r); // Eq. 32
        } else {
            vt = 2. * fKh * fGammaw * fHw / (B(r) * r); // Eq. 33
        }
    }
    /// Darcy velocity v'_opt in the paper coordinates (x right, y down); 0 outside the soil and for r >= R_e
    void VelocityPaper(REAL x, REAL yp, REAL v[2]) const {
        v[0] = v[1] = 0.;
        if (fHw <= 0. || !InSoilPaper(x, yp)) return;
        const REAL r = std::sqrt(x * x + yp * yp);
        if (!(r > 0.) || r >= fRe) return;
        REAL vr, vt;
        if (r < fRw) PolarVelocity(r, std::atan2(yp, -x), vr, vt);
        else PolarVelocity(r, 0., vr, vt); // zones 2 and 3 do not depend on theta
        const REAL ct = -x / r, st = yp / r;
        v[0] = -vr * ct + vt * st;
        v[1] = vr * st + vt * ct;
    }
    /// seepage force f = K^-1 v'_opt (kN/m^3) in the paper coordinates (x right, y down)
    void ForcePaper(REAL x, REAL yp, REAL f[2]) const {
        REAL v[2];
        VelocityPaper(x, yp, v);
        f[0] = v[0] / fKh, f[1] = fAlpha * v[1] / fKh;
    }
    /// seepage force in the NeoPZ coordinates (y up, y_paper = -y)
    void Force(REAL x, REAL y, REAL f[2]) const {
        REAL fp[2];
        ForcePaper(x, -y, fp);
        f[0] = fp[0], f[1] = -fp[1];
    }
    void Force(const TPZVec<REAL> &x, REAL f[2]) const { Force(x[0], x[1], f); }

    ForceField AsForceField(std::shared_ptr<const AnalyticalSeepage> self) const {
        return [self](const TPZVec<REAL> &x, REAL f[2]) { self->Force(x[0], x[1], f); };
    }

    std::string Summary() const {
        char buf[512];
        snprintf(buf, sizeof(buf),
                 "beta %.4g alpha %.4g h_w/H %.4g: m = sqrt(C/D) %.10g, C/h2e^2 %.10g, D/h2e^2 %.10g, F %.10g, 1/A %.10g, "
                 "I3 %.10g, -J*/(kh H^2 gw^2) %.10g%s",
                 fBetaDeg, fAlpha, fHw / fH, fM, fC / (fH2e * fH2e), fD / (fH2e * fH2e), fF, 1. / fA, fI3, -JstarNormalized(),
                 fDegenerate && fHw > 0. ? " [degenerate m -> 0: purely tangential zone-1 field]" : "");
        return buf;
    }

private:
    struct H2Rhs { ///< scaled h2 system, y = {h2, w = c h2' / m, int w^2 / c, int d h2^2}
        REAL alpha, m;
        void operator()(REAL t, const REAL y[4], REAL dy[4]) const {
            const REAL s = std::sin(t), s2 = s * s, c2 = 1. - s2;
            const REAL c = c2 + alpha * s2, d = s2 + alpha * c2;
            dy[0] = m * y[1] / c;
            dy[1] = d * y[0];
            dy[2] = y[1] * y[1] / c;
            dy[3] = d * y[0] * y[0];
        }
    };
    struct H2Node {
        REAL h, dh, d2h; ///< h2, h2', h2''
    };

    REAL fBetaDeg, fH, fHw, fAlpha, fKh, fGammaw, fLm;
    REAL fSb = 0., fCb = 0., fTb = 0., fTheta = 0., fXToe = 0., fR = 0., fRw = 0., fRe = 0., fA = 0., fI3 = 0.;
    bool fTanFinite = true, fDegenerate = false, fUnitM = false;
    REAL fM = 0., fC = 0., fD = 0., fH2e = 1., fPhi = 0., fF = 0., fAp = 0.;
    std::vector<H2Node> fTab;
    REAL fDT = 0., fInvDT = 0.;
    int fShootings = 0; ///< number of integrations of the h2 problem done by the construction

    void SolveH2() {
        // 1) scan Phi(m) on a log grid; 2) bounded Brent minimization around the smallest value; 3) sharpen with the
        //    root of the fixed-point equation m = sqrt(C/D) (Phi is flat at its minimum), as scripts/analytical_seepage.py
        const int ng = 81;
        const REAL t0 = std::log(1.e-7), t1 = std::log(30.);
        std::vector<REAL> tg(ng), ph(ng);
        int imin = 0;
        for (int k = 0; k < ng; k++) {
            tg[k] = t0 + (t1 - t0) * k / (ng - 1);
            ph[k] = PhiOfM(std::exp(tg[k]));
            fShootings++;
            if (ph[k] < ph[imin]) imin = k;
        }
        fDegenerate = imin == 0; // Phi increases from Phi(0+) = A: the infimum is the limit m -> 0 of the class
        if (fDegenerate) {
            fM = 0., fC = 0., fD = fA, fH2e = 1., fPhi = fA;
            return;
        }
        const REAL lo = tg[std::max(imin - 1, 0)], hi = tg[std::min(imin + 1, ng - 1)];
        REAL m = std::exp(numerics::BrentMinimize([this](REAL t) { return fShootings++, PhiOfM(std::exp(t)); }, lo, hi, 1.e-11));
        const REAL mlo = std::exp(lo), mhi = std::exp(hi);
        fShootings += 2;
        if (FixedPointResidual(mlo) * FixedPointResidual(mhi) < 0.)
            m = numerics::BrentRoot([this](REAL mm) { return fShootings++, FixedPointResidual(mm); }, mlo, mhi, 1.e-14, 1.e-14);
        fM = m;
        // tabulation of h2, h2', h2'' on a uniform grid (the integrator stops at every node)
        const int n = 2049;
        fTab.resize(n);
        fDT = fTheta / (n - 1), fInvDT = 1. / fDT;
        REAL y[4] = {1., 0., 0., 0.}, h = 0.;
        const H2Rhs rhs{fAlpha, m};
        auto node = [&](int i) {
            const REAL t = i * fDT, s = std::sin(t), s2 = s * s, c2 = 1. - s2;
            const REAL c = c2 + fAlpha * s2, d = s2 + fAlpha * c2, dc = (fAlpha - 1.) * std::sin(2. * t);
            const REAL dh = m * y[1] / c;
            fTab[i] = {y[0], dh, (m * d * y[0] - dc * dh) / c};
        };
        node(0);
        for (int i = 1; i < n; i++) {
            const REAL ta = (i - 1) * fDT, tb = i == n - 1 ? fTheta : i * fDT;
            if (numerics::DormandPrince54<4>(rhs, ta, tb, y, 1.e-13, 1.e-15 * fA, h) < 0)
                throw std::runtime_error("AnalyticalSeepage: h2 tabulation failed");
            node(i);
        }
        fShootings++;
        fC = m * m * y[2], fD = y[3], fH2e = y[0];
        fPhi = (1. + m) * y[1] / y[0];
    }
};

/// the analytical seepage force K^-1 v'_opt as a ForceField (NeoPZ coordinates)
inline ForceField AnalyticalForceField(std::shared_ptr<const AnalyticalSeepage> field) {
    return field->AsForceField(field);
}

} // namespace slope

#endif
