// Reference solutions of the verification tests (formulation document, Sect. "Plano de verificação").
#ifndef UNSATREFERENCE_H
#define UNSATREFERENCE_H

#include <cmath>
#include <vector>

namespace unsatref {

/**
 * Terzaghi: excess pore pressure p/q at the distance z from the drained face divided by the thickness H, time
 * factor T = c t / H^2 (series with 200 terms, as mcc::TerzaghiPressure of the former Projects/Common).
 */
inline double Terzaghi(double z, double T) {
    double p = 0.;
    for (int m = 0; m < 200; ++m) {
        const double mu = (2. * m + 1.) * M_PI / 2.;
        p += 2. / mu * std::sin(mu * z) * std::exp(-mu * mu * T);
    }
    return p;
}

/**
 * Steady vertical flow in a soil with k = k_s exp(-a p_c) (Gardner), water table at z = 0, downward flux q (m/s,
 * infiltration > 0): suction p_c(z) = -(1/a) ln[(1 - q/k_s) exp(-a gamma_w z) + q/k_s]. Obtained from Darcy,
 * k (d psi/dz + 1) = q with psi = -p_c/gamma_w; it is Eq. (9) of Li et al. (2024, p. 6) with their Q = -q.
 */
inline double GardnerSteadySuction(double z, double a, double gw, double q, double ks) {
    return -std::log((1. - q / ks) * std::exp(-a * gw * z) + q / ks) / a;
}

/// Positive roots of tan(x) = c x, one in each interval ((i-1) pi, (i-1/2) pi), i = 1..n (Mandel, c > 1)
inline std::vector<double> MandelRoots(double c, int n) {
    std::vector<double> r;
    for (int i = 1; i <= n; ++i) {
        // f = sin x - c x cos x has opposite signs at the ends of the interval
        double lo = (i - 1) * M_PI + 1.e-9, hi = (i - 0.5) * M_PI - 1.e-12;
        auto f = [c](double x) { return std::sin(x) - c * x * std::cos(x); };
        for (int k = 0; k < 200; ++k) {
            const double m = 0.5 * (lo + hi);
            if ((f(lo) < 0.) == (f(m) < 0.)) lo = m;
            else hi = m;
        }
        r.push_back(0.5 * (lo + hi));
    }
    return r;
}

/**
 * Mandel: pore pressure at the distance x from the centre, half-width a, under the total vertical force 2F per unit
 * length (F/a = |sigma_0|), consolidation coefficient c. Series of Appendix A.1 of Haghighat et al. (2021, p. 25):
 * p = 2|sigma_0| B (1 + nu_u)/3 sum sin(al)/(al - sin(al) cos(al)) (cos(al x/a) - cos(al)) exp(-al^2 c t/a^2).
 */
inline double MandelPressure(double x, double t, double a, double sigma0, double B, double nu, double nuu,
                             double c, int nterms = 100) {
    static thread_local std::vector<double> roots;
    static thread_local double croot = -1.;
    const double cc = (1. - nu) / (nuu - nu);
    if (cc != croot || (int)roots.size() < nterms) roots = MandelRoots(cc, nterms), croot = cc;
    double s = 0.;
    for (int i = 0; i < nterms; ++i) {
        const double al = roots[i];
        s += std::sin(al) / (al - std::sin(al) * std::cos(al)) * (std::cos(al * x / a) - std::cos(al)) *
             std::exp(-al * al * c * t / (a * a));
    }
    return 2. * std::fabs(sigma0) * B * (1. + nuu) / 3. * s;
}

} // namespace unsatref

#endif
