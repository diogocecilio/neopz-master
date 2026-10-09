/**
 * @file TPZYCModifiedCamClayRHW.cpp
 * @brief Implementation of the Modified Cam-Clay closest-point projection in rotated
 * Haigh-Westergaard space (see TPZYCModifiedCamClayRHW.h): reduced local problem in the angle of the meridian
 * ellipse and the hardening increment (default), and the four-unknown system kept as a cross-check.
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
 * consistency row of the Jacobian of the four-unknown system has a zero diagonal entry.
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
      fPorousIntegration(EExact), fNewtonTol(1.e-12), fMaxNewton(50), fLocalSolver(EReduced) {
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

bool TPZYCModifiedCamClayRHW::GradProjection(const TPZVec<REAL> &X, const TTrial &trial, const TPZVec<REAL> &n,
                                             bool isotropic, TPZFMatrix<REAL> &Dproj) const {
    TPZFNMatrix<16, REAL> J(4, 4, 0.), dRdY(4, 2, 0.);
    Jacobian(X, trial, J, dRdY);
    // dX/dY = -J^{-1} dR/dY
    TPZFNMatrix<8, REAL> dXdY(dRdY);
    if (!SolvePivoting(J, dXdY)) return false;
    dXdY *= -1.;
    const REAL f = isotropic ? 1. / (1. + 6. * trial.fG * X[3] / (fM * fM)) : X[1] / trial.fRhoTr;
    const REAL dxi[2] = {dXdY(0, 0), dXdY(0, 1)}, drho[2] = {dXdY(1, 0), dXdY(1, 1)};
    AssembleDproj(dxi, drho, f, n, Dproj);
    return true;
}

void TPZYCModifiedCamClayRHW::AssembleDproj(const REAL dxi[2], const REAL drho[2], REAL f, const TPZVec<REAL> &n,
                                            TPZFMatrix<REAL> &Dproj) {
    const REAL sq3 = std::sqrt(3.);
    Dproj.Redim(3, 3);
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            const REAL delta = (i == j) ? 1. : 0.;
            Dproj(i, j) = dxi[0] / 3. + dxi[1] / sq3 * n[j] + drho[0] / sq3 * n[i] + drho[1] * n[i] * n[j] +
                          f * (delta - 1. / 3. - n[i] * n[j]);
        }
    }
}

// ------------------------------------------------------------------------------------ reduced local problem

REAL TPZYCModifiedCamClayRHW::ElasticHardeningIncrement(REAL xi, REAL xitr, REAL v0, REAL &dxi, REAL &dxitr) const {
    const REAL sq3 = std::sqrt(3.);
    if (fVolumetricLaw == ELinear) {
        dxi = 1. / (sq3 * fK0);
        dxitr = -dxi;
        return (xi - xitr) / (sq3 * fK0);
    }
    const REAL c = fKappa / v0;
    if (fPorousIntegration == EFrozen) {
        dxi = -c / xitr;
        dxitr = c * xi / (xitr * xitr);
        return c * (1. - xi / xitr);
    }
    dxi = -c / xi;
    dxitr = c / xitr;
    return -c * std::log(xi / xitr);
}

void TPZYCModifiedCamClayRHW::SurfacePoint(REAL theta, REAL a, REAL b, REAL &xi, REAL &rho, REAL dxi[2],
                                           REAL drho[2]) const {
    const REAL sq3 = std::sqrt(3.), s23 = std::sqrt(2. / 3.);
    const REAL st = std::sin(theta), ct = std::cos(theta);
    xi = sq3 * (fPt - a - a * b * ct);
    rho = s23 * fM * a * st;
    dxi[0] = sq3 * a * b * st;
    dxi[1] = -sq3 * (1. + b * ct);
    drho[0] = s23 * fM * a * ct;
    drho[1] = s23 * fM * st;
}

void TPZYCModifiedCamClayRHW::ResidualReduced(const TPZVec<REAL> &X2, const TTrial &trial, TPZVec<REAL> &E) const {
    const REAL s23 = std::sqrt(2. / 3.);
    const REAL theta = X2[0], dal = X2[1], b = trial.fB, G = trial.fG;
    REAL a, H, pc;
    Hardening(trial.fPcn, dal, trial.fV0, a, H, pc);
    REAL xi, rho, dxi[2], drho[2];
    SurfacePoint(theta, a, b, xi, rho, dxi, drho);
    REAL de_xi, de_xitr;
    const REAL de = ElasticHardeningIncrement(xi, trial.fXiTr, trial.fV0, de_xi, de_xitr);
    E.Resize(2);
    // E1 = (1/a) d(d^2/2)/dtheta = (de/sqrt3)(dxi/dtheta)/a + (rho - rhotr)/(2G) (drho/dtheta)/a
    E[0] = de * b * std::sin(theta) + (rho - trial.fRhoTr) * s23 * fM * std::cos(theta) / (2. * G);
    // E2: volumetric elastic law
    E[1] = dal - de;
}

void TPZYCModifiedCamClayRHW::JacobianReduced(const TPZVec<REAL> &X2, const TTrial &trial, TPZFMatrix<REAL> &J,
                                              TPZFMatrix<REAL> &dEdY) const {
    const REAL s23 = std::sqrt(2. / 3.);
    const REAL theta = X2[0], dal = X2[1], b = trial.fB, G = trial.fG;
    const REAL st = std::sin(theta), ct = std::cos(theta);
    REAL a, H, pc;
    Hardening(trial.fPcn, dal, trial.fV0, a, H, pc);
    REAL xi, rho, dxi[2], drho[2];
    SurfacePoint(theta, a, b, xi, rho, dxi, drho);
    REAL de_xi, de_xitr;
    const REAL de = ElasticHardeningIncrement(xi, trial.fXiTr, trial.fV0, de_xi, de_xitr);
    J.Redim(2, 2);
    dEdY.Redim(2, 2);
    J(0, 0) = de_xi * dxi[0] * b * st + de * b * ct +
              (drho[0] * s23 * fM * ct - (rho - trial.fRhoTr) * s23 * fM * st) / (2. * G);
    J(0, 1) = H * (de_xi * dxi[1] * b * st + drho[1] * s23 * fM * ct / (2. * G));
    J(1, 0) = -de_xi * dxi[0];
    J(1, 1) = 1. - de_xi * dxi[1] * H;
    dEdY(0, 0) = de_xitr * b * st;
    dEdY(0, 1) = -s23 * fM * ct / (2. * G);
    dEdY(1, 0) = -de_xitr;
    dEdY(1, 1) = 0.;
}

bool TPZYCModifiedCamClayRHW::ProjectReduced(const TTrial &trial, TPZVec<REAL> &X2, int &niter) const {
    const REAL sq3 = std::sqrt(3.);
    REAL an, H, pc;
    Hardening(trial.fPcn, 0., trial.fV0, an, H, pc);
    const REAL pbtr = trial.fXiTr / sq3 - fPt + an;
    // start: angle of the trial state seen from the centre of the ellipse of a_n, no hardening
    X2.Resize(2);
    X2[0] = std::atan2(trial.fRhoTr * std::sqrt(1.5) / (fM * an), -pbtr / (an * trial.fB));
    X2[1] = 0.;
    TPZManVector<REAL, 2> E(2, 0.), En(2, 0.), Xn(2, 0.);
    TPZFNMatrix<4, REAL> J(2, 2, 0.), dEdY(2, 2, 0.);
    ResidualReduced(X2, trial, E);
    REAL norm = std::hypot(E[0], E[1]);
    for (niter = 0; niter < fMaxNewton; ++niter) {
        if (!std::isfinite(norm)) return false;
        if (norm < fNewtonTol) return true;
        JacobianReduced(X2, trial, J, dEdY);
        const REAL det = J(0, 0) * J(1, 1) - J(0, 1) * J(1, 0);
        const REAL scale = std::max(std::max(std::fabs(J(0, 0)), std::fabs(J(0, 1))),
                                    std::max(std::fabs(J(1, 0)), std::fabs(J(1, 1))));
        if (!(std::fabs(det) > 1.e-14 * scale * scale)) return false;
        const REAL d0 = (J(1, 1) * E[0] - J(0, 1) * E[1]) / det;
        const REAL d1 = (-J(1, 0) * E[0] + J(0, 0) * E[1]) / det;
        // backtracking: the step is halved until the norm of the residual decreases and theta stays in [0, pi]
        bool accepted = false;
        REAL lam = 1.;
        for (int k = 0; k < 40; ++k, lam *= 0.5) {
            Xn[0] = X2[0] - lam * d0;
            Xn[1] = X2[1] - lam * d1;
            if (Xn[0] < 0. || Xn[0] > M_PI) continue;
            ResidualReduced(Xn, trial, En);
            const REAL nn = std::hypot(En[0], En[1]);
            if (std::isfinite(nn) && nn < norm) {
                accepted = true;
                norm = nn;
                break;
            }
        }
        if (!accepted) return false;
        X2[0] = Xn[0];
        X2[1] = Xn[1];
        E[0] = En[0];
        E[1] = En[1];
    }
    return false;
}

void TPZYCModifiedCamClayRHW::ReducedToFull(const TPZVec<REAL> &X2, const TTrial &trial, TPZVec<REAL> &X) const {
    const REAL sq3 = std::sqrt(3.);
    const REAL theta = X2[0], dal = X2[1], b = trial.fB;
    REAL a, H, pc;
    Hardening(trial.fPcn, dal, trial.fV0, a, H, pc);
    REAL xi, rho, dxi[2], drho[2];
    SurfacePoint(theta, a, b, xi, rho, dxi, drho);
    X.Resize(4);
    X[0] = xi;
    X[1] = rho;
    X[2] = dal;
    // plastic multiplier: from the deviatoric law, rho (1 + 6 G dg/M^2) = rho_tr, where the surface is not close
    // to the apexes (sin theta >= |cos theta|), and from the hardening rule, dal = -2 dg pbar/b^2, near them
    const REAL st = std::sin(theta), ct = std::cos(theta);
    if (st >= std::fabs(ct) && trial.fRhoTr > 0. && rho > 0.) {
        X[3] = fM * fM * (trial.fRhoTr / rho - 1.) / (6. * trial.fG);
    } else {
        const REAL pbar = xi / sq3 - fPt + a;
        X[3] = -dal * b * b / (2. * pbar);
    }
}

bool TPZYCModifiedCamClayRHW::GradProjectionReduced(const TPZVec<REAL> &X2, const TTrial &trial, const TPZVec<REAL> &n,
                                                    bool isotropic, TPZFMatrix<REAL> &Dproj) const {
    TPZFNMatrix<4, REAL> J(2, 2, 0.), dEdY(2, 2, 0.);
    JacobianReduced(X2, trial, J, dEdY);
    const REAL det = J(0, 0) * J(1, 1) - J(0, 1) * J(1, 0);
    const REAL scale = std::max(std::max(std::fabs(J(0, 0)), std::fabs(J(0, 1))),
                                std::max(std::fabs(J(1, 0)), std::fabs(J(1, 1))));
    if (!(std::fabs(det) > 1.e-14 * scale * scale)) return false;
    // dX/dY = -J^{-1} dE/dY, X = [theta, dal], Y = [xi_tr, rho_tr]
    REAL dXdY[2][2];
    for (int j = 0; j < 2; ++j) {
        dXdY[0][j] = -(J(1, 1) * dEdY(0, j) - J(0, 1) * dEdY(1, j)) / det;
        dXdY[1][j] = -(-J(1, 0) * dEdY(0, j) + J(0, 0) * dEdY(1, j)) / det;
    }
    // chain rule through xi(theta, a) and rho(theta, a), a = a(dal)
    REAL a, H, pc;
    Hardening(trial.fPcn, X2[1], trial.fV0, a, H, pc);
    REAL xi, rho, dxi[2], drho[2];
    SurfacePoint(X2[0], a, trial.fB, xi, rho, dxi, drho);
    REAL dxiY[2], drhoY[2];
    for (int j = 0; j < 2; ++j) {
        dxiY[j] = dxi[0] * dXdY[0][j] + dxi[1] * H * dXdY[1][j];
        drhoY[j] = drho[0] * dXdY[0][j] + drho[1] * H * dXdY[1][j];
    }
    TPZManVector<REAL, 4> X(4, 0.);
    ReducedToFull(X2, trial, X);
    const REAL f = isotropic ? 1. / (1. + 6. * trial.fG * X[3] / (fM * fM)) : rho / trial.fRhoTr;
    AssembleDproj(dxiY, drhoY, f, n, Dproj);
    return true;
}

bool TPZYCModifiedCamClayRHW::Project(const TTrial &trial, TPZVec<REAL> &X, TPZVec<REAL> &X2, int &niter) const {
    if (fLocalSolver == EFull) {
        X2.Resize(0);
        return ProjectHW(trial, X, niter);
    }
    if (!ProjectReduced(trial, X2, niter)) return false;
    ReducedToFull(X2, trial, X);
    return true;
}

bool TPZYCModifiedCamClayRHW::ProjectionJacobian(const TPZVec<REAL> &X, const TPZVec<REAL> &X2, const TTrial &trial,
                                                 const TPZVec<REAL> &n, bool isotropic, TPZFMatrix<REAL> &Dproj) const {
    if (fLocalSolver == EFull) return GradProjection(X, trial, n, isotropic, Dproj);
    return GradProjectionReduced(X2, trial, n, isotropic, Dproj);
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
        << " local solver = " << (fLocalSolver == EReduced ? "reduced (theta, dal)" : "full (xi, rho, dal, dg)")
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
    int solver = fLocalSolver;
    buf.Write(&solver);
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
    int solver;
    buf.Read(&solver);
    fLocalSolver = ELocalSolver(solver);
}
