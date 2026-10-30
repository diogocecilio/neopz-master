// Factors of safety with the pore pressures of an instant frozen (only the momentum balance is solved):
// strength reduction (SRM) from a geostatic state with the strength increased (F0), and gravity increase (GIM)
// from the stress-free state with the real strength (port of Projects/SlopeRainfall/main.cpp, branch
// claude/relaxed-volta-1r4gog). Formulation document, Sect. "Fator de segurança: SRM e GIM".
#ifndef UNSATFS_H
#define UNSATFS_H

#include "UnsatAnalysis.h"

#include <functional>
#include <limits>

/// Follows xs; at the first x without convergence, bisection down to tol; the failed x is tried again from the last
/// converged state. Returns the last x that converged.
inline REAL CollapseSearch(const std::function<bool(REAL)> &solveAt, REAL x0, const std::vector<REAL> &xs, REAL tol,
                           REAL xmax) {
    REAL xok = x0, xfail = std::numeric_limits<REAL>::infinity();
    while (true) {
        REAL x;
        if (std::isfinite(xfail) && xfail - xok > tol) x = 0.5 * (xok + xfail);
        else if (std::isfinite(xfail)) x = xfail;
        else {
            x = xok + (xs.back() - xs[xs.size() - 2]);
            for (REAL v : xs)
                if (v > xok + 1.e-12) {
                    x = v;
                    break;
                }
        }
        if (x > xmax) return xok;
        const bool conv = solveAt(x);
        if (conv) {
            if (x == xfail) xfail = std::numeric_limits<REAL>::infinity();
            xok = x;
        } else if (std::isfinite(xfail) && x == xfail) break;
        else xfail = x;
    }
    return xok;
}

/// Displacement constraints and all the pressures at their values in x
inline std::map<int64_t, REAL> FrozenFixed(UnsatAnalysis &an, const TPZFMatrix<STATE> &x) {
    std::map<int64_t, REAL> fx = an.DisplacementConstraints();
    for (int64_t i = 0; i < an.NEq(); ++i)
        if (an.IsPressure(i)) fx[i] = x.GetVal(i, 0);
    return fx;
}

/// Memory back to the stress-free state, displacements zero, pressures of x
inline void ResetState(UnsatAnalysis &an, const TPZFMatrix<STATE> &x) {
    an.ResetMemory();
    TPZFMatrix<STATE> &X = an.X();
    X = x;
    for (int64_t i = 0; i < an.NEq(); ++i)
        if (!an.IsPressure(i)) X(i, 0) = 0.;
    an.Restore();
}

/// Weight (and the equivalent pressure, scaled together) in increments of 1/n with frozen pressures, factor F
inline bool Geostatic(UnsatAnalysis &an, REAL F, int n = 4) {
    an.SetStrengthReductionFactor(F);
    an.SetGravityScalesSuction(true);
    const auto fx = FrozenFixed(an, an.X());
    for (int k = 1; k <= n; ++k) {
        an.SetGravityFactor(REAL(k) / n);
        if (!an.Solve(0., fx, {}, true).conv) return false;
        an.Accept();
    }
    return true;
}

struct FSResult {
    REAL fs = -1.;
    int solves = 0;
};

/// SRM with the pressures of x: geostatic with F0 from the stress-free state, then F increased (list xs) and bisection
inline FSResult FactorSRM(UnsatAnalysis &an, const TPZFMatrix<STATE> &x, REAL F0 = 0.6, REAL tol = 0.01,
                          std::vector<REAL> xs = {0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.5, 3.0, 3.5, 4.0}) {
    FSResult r;
    ResetState(an, x);
    if (!Geostatic(an, F0)) return r;
    const auto fx = FrozenFixed(an, an.X());
    an.SetGravityFactor(1.);
    std::vector<REAL> l;
    for (REAL v : xs)
        if (v > F0 + 1.e-9) l.push_back(v);
    r.fs = CollapseSearch([&](REAL F) {
        r.solves++;
        an.SetStrengthReductionFactor(F);
        const bool ok = an.Solve(0., fx, {}, true).conv;
        if (ok) an.Accept();
        else an.Restore();
        return ok;
    }, F0, l, tol, 20.);
    an.SetStrengthReductionFactor(1.);
    return r;
}

/// GIM with the pressures of x: from the stress-free state, real strength, lambda on the weight and on the positive
/// equivalent pressure (the suction keeps its value)
inline FSResult FactorGIM(UnsatAnalysis &an, const TPZFMatrix<STATE> &x, REAL tol = 0.01,
                          std::vector<REAL> xs = {0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 10.0}) {
    FSResult r;
    ResetState(an, x);
    an.SetStrengthReductionFactor(1.);
    an.SetGravityScalesSuction(false);
    const auto fx = FrozenFixed(an, an.X());
    r.fs = CollapseSearch([&](REAL lam) {
        r.solves++;
        an.SetGravityFactor(lam);
        const bool ok = an.Solve(0., fx, {}, true).conv;
        if (ok) an.Accept();
        else an.Restore();
        return ok;
    }, 0., xs, tol, 20.);
    an.SetGravityScalesSuction(true);
    an.SetGravityFactor(1.);
    return r;
}

#endif
