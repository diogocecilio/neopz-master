/**
 * @file TPZYCModifiedCamClayRHW.cpp
 * @brief Implementation of the Modified Cam-Clay closest-point projection in rotated
 * Haigh-Westergaard space (see TPZYCModifiedCamClayRHW.h).
 */

#include "TPZYCModifiedCamClayRHW.h"
#include "TPZHash.h"
#include "TPZStream.h"
#include "pzerror.h"
#include <cmath>
#include <algorithm>

namespace {
/**
 * @brief Solves the small dense system A x = B (B with several columns) by Gaussian elimination
 * with partial pivoting. A and B are overwritten. Returns false for a (numerically) singular A.
 * @note TPZFMatrix::Solve_LU does not pivot unless the library is built with LAPACK, and the
 * consistency row of the local Jacobian (A.1) has a zero diagonal entry.
 */
bool SolvePivoting(TPZFMatrix<REAL> &A, TPZFMatrix<REAL> &B) {
    const int n = A.Rows(), nrhs = B.Cols();
    REAL scale = 0.;
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) scale = std::max(scale, std::fabs(A(i, j)));
    if (scale == 0. || !std::isfinite(scale)) return false;
    for (int k = 0; k < n; ++k) {
        int piv = k;
        for (int i = k + 1; i < n; ++i)
            if (std::fabs(A(i, k)) > std::fabs(A(piv, k))) piv = i;
        if (std::fabs(A(piv, k)) < 1.e-14 * scale) return false;
        if (piv != k) {
            for (int j = 0; j < n; ++j) std::swap(A(k, j), A(piv, j));
            for (int j = 0; j < nrhs; ++j) std::swap(B(k, j), B(piv, j));
        }
        for (int i = k + 1; i < n; ++i) {
            const REAL f = A(i, k) / A(k, k);
            if (f == 0.) continue;
            for (int j = k; j < n; ++j) A(i, j) -= f * A(k, j);
            for (int j = 0; j < nrhs; ++j) B(i, j) -= f * B(k, j);
        }
    }
    for (int k = n - 1; k >= 0; --k) {
        for (int j = 0; j < nrhs; ++j) {
            REAL s = B(k, j);
            for (int i = k + 1; i < n; ++i) s -= A(k, i) * B(i, j);
            B(k, j) = s / A(k, k);
        }
    }
    return true;
}
} // namespace

TPZYCModifiedCamClayRHW::TPZYCModifiedCamClayRHW()
    : fM(1.), fLambda(0.2), fKappa(0.05), fPt(0.), fOmega(1.), fK0(0.), fVolumetricLaw(EPorous),
      fPorousIntegration(EExact), fNewtonTol(1.e-12), fMaxNewton(50) {
}

void TPZYCModifiedCamClayRHW::SetUp(REAL M, REAL lambda, REAL kappa, REAL pt, REAL omega) {
    if (lambda <= kappa || kappa <= 0. || M <= 0. || omega <= 0.) {
        std::cout << __PRETTY_FUNCTION__ << " invalid parameters M = " << M << " lambda = " << lambda
                  << " kappa = " << kappa << " omega = " << omega << std::endl;
        DebugStop();
    }
    fM = M;
    fLambda = lambda;
    fKappa = kappa;
    fPt = pt;
    fOmega = omega;
}

void TPZYCModifiedCamClayRHW::Hardening(REAL pcn, REAL dal, REAL v0, REAL &a, REAL &H, REAL &pc) const {
    const REAL c = v0 / (fLambda - fKappa);
    pc = pcn * std::exp(c * dal);
    a = (pc + fPt) / (1. + fOmega);
    H = c * pc / (1. + fOmega);
}

REAL TPZYCModifiedCamClayRHW::PhiCC(REAL p, REAL rho, REAL a, REAL b) const {
    const REAL pbar = p - fPt + a;
    return pbar * pbar / (b * b) + 1.5 * rho * rho / (fM * fM) - a * a;
}

void TPZYCModifiedCamClayRHW::Residual(const TPZVec<REAL> &X, const TTrial &trial, TPZVec<REAL> &R) const {
    const REAL sq3 = std::sqrt(3.);
    const REAL xi = X[0], rho = X[1], dal = X[2], dg = X[3];
    const REAL b = trial.fB;
    REAL an, a, H, pc;
    Hardening(trial.fPcn, 0., trial.fV0, an, H, pc);
    Hardening(trial.fPcn, dal, trial.fV0, a, H, pc);
    const REAL pbar = xi / sq3 - fPt + a;
    const REAL M2 = fM * fM;
    R.Resize(4);
    if (fVolumetricLaw == ELinear) {
        R[0] = (xi - trial.fXiTr - sq3 * fK0 * dal) / an;
    } else if (fPorousIntegration == EFrozen) {
        R[0] = (xi - trial.fXiTr * (1. - trial.fV0 * dal / fKappa)) / an;
    } else {
        R[0] = (xi - trial.fXiTr * std::exp(-trial.fV0 * dal / fKappa)) / an;
    }
    R[1] = (rho - trial.fRhoTr + 6. * trial.fG * dg * rho / M2) / an;
    R[2] = dal + dg * 2. * pbar / (b * b);
    R[3] = (pbar * pbar / (b * b) + 1.5 * rho * rho / M2 - a * a) / (an * an);
}

void TPZYCModifiedCamClayRHW::Jacobian(const TPZVec<REAL> &X, const TTrial &trial, TPZFMatrix<REAL> &J,
                                       TPZFMatrix<REAL> &dRdY) const {
    const REAL sq3 = std::sqrt(3.);
    const REAL xi = X[0], rho = X[1], dal = X[2], dg = X[3];
    const REAL b = trial.fB, b2 = b * b;
    REAL an, a, H, pc;
    Hardening(trial.fPcn, 0., trial.fV0, an, H, pc);
    Hardening(trial.fPcn, dal, trial.fV0, a, H, pc);
    const REAL pbar = xi / sq3 - fPt + a;
    const REAL M2 = fM * fM;
    const REAL G = trial.fG;
    J.Redim(4, 4);
    dRdY.Redim(4, 2);
    REAL e1 = 1.;
    // row 1: volumetric elastic law
    if (fVolumetricLaw == ELinear) {
        J(0, 0) = 1.;
        J(0, 2) = -sq3 * fK0;
        e1 = 1.;
    } else if (fPorousIntegration == EFrozen) {
        J(0, 0) = 1.;
        J(0, 2) = trial.fXiTr * trial.fV0 / fKappa;
        e1 = 1. - trial.fV0 * dal / fKappa;
    } else {
        e1 = std::exp(-trial.fV0 * dal / fKappa);
        J(0, 0) = 1.;
        J(0, 2) = trial.fXiTr * trial.fV0 / fKappa * e1;
    }
    J(0, 0) /= an;
    J(0, 2) /= an;
    // row 2: deviatoric elastic law with the flow rule
    J(1, 1) = (1. + 6. * G * dg / M2) / an;
    J(1, 3) = 6. * G * rho / M2 / an;
    // row 3: hardening rule with the flow rule
    J(2, 0) = dg * 2. / b2 / sq3;
    J(2, 2) = 1. + dg * 2. / b2 * H;
    J(2, 3) = 2. * pbar / b2;
    // row 4: consistency condition
    J(3, 0) = 2. * pbar / b2 / sq3 / (an * an);
    J(3, 1) = 3. * rho / M2 / (an * an);
    J(3, 2) = (2. * pbar / b2 * H - 2. * a * H) / (an * an);
    // sensitivity with respect to the trial coordinates: only R1 depends on xitr and R2 on rhotr
    dRdY(0, 0) = -e1 / an;
    dRdY(1, 1) = -1. / an;
}

bool TPZYCModifiedCamClayRHW::ProjectHW(const TTrial &trial, TPZVec<REAL> &X, int &niter) const {
    X.Resize(4);
    X[0] = trial.fXiTr;
    X[1] = trial.fRhoTr;
    X[2] = 0.;
    X[3] = 0.;
    TPZManVector<REAL, 4> R(4, 0.);
    TPZFNMatrix<16, REAL> J(4, 4, 0.), dRdY(4, 2, 0.);
    TPZFNMatrix<4, REAL> dX(4, 1, 0.);
    for (niter = 0; niter < fMaxNewton; ++niter) {
        Residual(X, trial, R);
        REAL norm = 0.;
        for (int i = 0; i < 4; ++i) norm += R[i] * R[i];
        norm = std::sqrt(norm);
        if (!std::isfinite(norm)) return false;
        if (norm < fNewtonTol) return true;
        Jacobian(X, trial, J, dRdY);
        for (int i = 0; i < 4; ++i) dX(i, 0) = R[i];
        if (!SolvePivoting(J, dX)) return false;
        for (int i = 0; i < 4; ++i) X[i] -= dX(i, 0);
    }
    return false;
}

void TPZYCModifiedCamClayRHW::GradProjection(const TPZVec<REAL> &X, const TTrial &trial, const TPZVec<REAL> &n,
                                             bool isotropic, TPZFMatrix<REAL> &Dproj) const {
    const REAL sq3 = std::sqrt(3.);
    TPZFNMatrix<16, REAL> J(4, 4, 0.), dRdY(4, 2, 0.);
    Jacobian(X, trial, J, dRdY);
    // dX/dY = -J^{-1} dR/dY
    TPZFNMatrix<8, REAL> dXdY(dRdY);
    if (!SolvePivoting(J, dXdY)) {
        std::cout << __PRETTY_FUNCTION__ << " singular local Jacobian" << std::endl;
        DebugStop();
    }
    dXdY *= -1.;
    const REAL f = isotropic ? 1. / (1. + 6. * trial.fG * X[3] / (fM * fM)) : X[1] / trial.fRhoTr;
    Dproj.Redim(3, 3);
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            const REAL delta = (i == j) ? 1. : 0.;
            Dproj(i, j) = dXdY(0, 0) / 3. + dXdY(0, 1) / sq3 * n[j] + dXdY(1, 0) / sq3 * n[i] +
                          dXdY(1, 1) * n[i] * n[j] + f * (delta - 1. / 3. - n[i] * n[j]);
        }
    }
}

void TPZYCModifiedCamClayRHW::YieldFunction(const TPZVec<STATE> &sigma, STATE kprev, TPZVec<STATE> &yield) const {
    const REAL p = (sigma[0] + sigma[1] + sigma[2]) / 3.;
    REAL rho2 = 0.;
    for (int i = 0; i < 3; ++i) rho2 += (sigma[i] - p) * (sigma[i] - p);
    const REAL a = (kprev + fPt) / (1. + fOmega);
    yield.Resize(NYield);
    yield[0] = PhiCC(p, std::sqrt(rho2), a, BFromP(p, a));
}

void TPZYCModifiedCamClayRHW::Print(std::ostream &out) const {
    out << "TPZYCModifiedCamClayRHW: M = " << fM << " lambda = " << fLambda << " kappa = " << fKappa
        << " pt = " << fPt << " omega = " << fOmega
        << " volumetric law = " << (fVolumetricLaw == ELinear ? "linear (K0 = " : "porous (K0 = ") << fK0 << ")"
        << " porous integration = " << (fPorousIntegration == EExact ? "exact" : "frozen")
        << " Newton tol = " << fNewtonTol << " max it = " << fMaxNewton << std::endl;
}

int TPZYCModifiedCamClayRHW::ClassId() const {
    return Hash("TPZYCModifiedCamClayRHW");
}

void TPZYCModifiedCamClayRHW::Write(TPZStream &buf, int withclassid) const {
    buf.Write(&fM);
    buf.Write(&fLambda);
    buf.Write(&fKappa);
    buf.Write(&fPt);
    buf.Write(&fOmega);
    buf.Write(&fK0);
    int law = fVolumetricLaw, integ = fPorousIntegration;
    buf.Write(&law);
    buf.Write(&integ);
    buf.Write(&fNewtonTol);
    buf.Write(&fMaxNewton);
}

void TPZYCModifiedCamClayRHW::Read(TPZStream &buf, void *context) {
    buf.Read(&fM);
    buf.Read(&fLambda);
    buf.Read(&fKappa);
    buf.Read(&fPt);
    buf.Read(&fOmega);
    buf.Read(&fK0);
    int law, integ;
    buf.Read(&law);
    buf.Read(&integ);
    fVolumetricLaw = EVolumetricLaw(law);
    fPorousIntegration = EPorousIntegration(integ);
    buf.Read(&fNewtonTol);
    buf.Read(&fMaxNewton);
}
