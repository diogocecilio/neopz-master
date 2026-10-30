// Taylor test of the global Jacobian (test 0 of the plan): E(alpha) = |R(x + alpha d) - R(x) - alpha J d| on the free
// equations should decrease as alpha^2 (formulation document, Sect. "Plano de verificação").
#ifndef UNSATTAYLOR_H
#define UNSATTAYLOR_H

#include "UnsatAnalysis.h"

#include <algorithm>
#include <random>

struct TaylorResult {
    std::vector<REAL> orders; ///< slope between consecutive amplitudes, per direction
    REAL median = 0., min = 0.;
};

/// Random directions d (scaled per field: scaleU for displacements, scaleP for pressures; zero on the fixed
/// equations); amplitudes alpha and alpha/2
inline TaylorResult TaylorTest(UnsatAnalysis &an, const TPZFMatrix<STATE> &x, REAL dt, bool frozen,
                               const std::map<int64_t, REAL> &fixed, REAL scaleU, REAL scaleP, int ndir = 5,
                               REAL alpha = 1.e-2, unsigned seed = 7) {
    std::mt19937 gen(seed);
    std::uniform_real_distribution<REAL> U(-1., 1.);
    TaylorResult res;
    const int64_t neq = an.NEq();
    const TPZFMatrix<STATE> R0 = an.Residual(x, dt, frozen);
    for (int k = 0; k < ndir; ++k) {
        TPZFMatrix<STATE> d(neq, 1, 0.);
        for (int64_t i = 0; i < neq; ++i)
            if (!fixed.count(i)) d(i, 0) = U(gen) * (an.IsPressure(i) ? scaleP : scaleU);
        const TPZFMatrix<STATE> Jd = an.TangentTimes(x, dt, frozen, d);
        REAL E[2];
        for (int s = 0; s < 2; ++s) {
            const REAL a = alpha / (s == 0 ? 1. : 2.);
            TPZFMatrix<STATE> xa(x);
            for (int64_t i = 0; i < neq; ++i) xa(i, 0) += a * d.GetVal(i, 0);
            const TPZFMatrix<STATE> Ra = an.Residual(xa, dt, frozen);
            REAL e2 = 0.;
            for (int64_t i = 0; i < neq; ++i) {
                if (fixed.count(i)) continue;
                if (frozen && an.IsPressure(i)) continue;
                e2 += std::pow(Ra.GetVal(i, 0) - R0.GetVal(i, 0) - a * Jd.GetVal(i, 0), 2);
            }
            E[s] = std::sqrt(e2);
        }
        res.orders.push_back(std::log(E[0] / E[1]) / std::log(2.));
    }
    an.LoadSolution(x);
    std::vector<REAL> o(res.orders);
    std::sort(o.begin(), o.end());
    res.median = o[o.size() / 2];
    res.min = o.front();
    return res;
}

#endif
