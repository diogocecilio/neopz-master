/**
 * @file TPZPlasticStepModifiedCamClay.cpp
 * @brief Implementation of the Modified Cam-Clay stress update (see TPZPlasticStepModifiedCamClay.h).
 */

#include "TPZPlasticStepModifiedCamClay.h"
#include "TPZHash.h"
#include "TPZStream.h"
#include "pzerror.h"
#include <cmath>
#include <algorithm>

namespace {

/** @brief Voigt vector (XX, XY, XZ, YY, YZ, ZZ) of the symmetric tensor a b^T + b a^T scaled by 1/2 (stress type) */
inline void SymVoigt(const TPZManVector<REAL, 3> &a, const TPZManVector<REAL, 3> &b, REAL v[6]) {
    v[_XX_] = a[0] * b[0];
    v[_XY_] = 0.5 * (a[0] * b[1] + a[1] * b[0]);
    v[_XZ_] = 0.5 * (a[0] * b[2] + a[2] * b[0]);
    v[_YY_] = a[1] * b[1];
    v[_YZ_] = 0.5 * (a[1] * b[2] + a[2] * b[1]);
    v[_ZZ_] = a[2] * b[2];
}

/** @brief Cyclic Jacobi rotations for a symmetric 3x3 matrix: eigenvalues on the diagonal, eigenvectors in the columns of V */
void Jacobi3(REAL A[3][3], REAL V[3][3]) {
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) V[i][j] = (i == j) ? 1. : 0.;
    REAL norm = 0.;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) norm += A[i][j] * A[i][j];
    for (int sweep = 0; sweep < 60; ++sweep) {
        const REAL off = A[0][1] * A[0][1] + A[0][2] * A[0][2] + A[1][2] * A[1][2];
        if (off <= 1.e-32 * norm || off == 0.) break;
        const int pq[3][2] = {{0, 1}, {0, 2}, {1, 2}};
        for (auto &e : pq) {
            const int p = e[0], q = e[1];
            if (A[p][q] == 0.) continue;
            const REAL theta = (A[q][q] - A[p][p]) / (2. * A[p][q]);
            const REAL t = (theta >= 0. ? 1. : -1.) / (std::fabs(theta) + std::sqrt(theta * theta + 1.));
            const REAL c = 1. / std::sqrt(t * t + 1.), s = t * c;
            for (int k = 0; k < 3; ++k) { // A <- A J
                const REAL akp = A[k][p], akq = A[k][q];
                A[k][p] = c * akp - s * akq;
                A[k][q] = s * akp + c * akq;
            }
            for (int k = 0; k < 3; ++k) { // A <- J^T A
                const REAL apk = A[p][k], aqk = A[q][k];
                A[p][k] = c * apk - s * aqk;
                A[q][k] = s * apk + c * aqk;
            }
            for (int k = 0; k < 3; ++k) { // V <- V J
                const REAL vkp = V[k][p], vkq = V[k][q];
                V[k][p] = c * vkp - s * vkq;
                V[k][q] = s * vkp + c * vkq;
            }
        }
    }
}

} // namespace

TPZPlasticStepModifiedCamClay::TPZPlasticStepModifiedCamClay()
    : fYC(), fER(), fModel(EModifiedCamClay), fShear(EPoisson), fG(0.), fNu(0.3), fV0Default(2.),
      fTangentMode(EConsistentTangent), fFDStep(1.e-7), fN(), fSigma(), fLastNewtonIterations(0), fFailed(false) {
}

const char *TPZPlasticStepModifiedCamClay::TangentModeName(ETangentMode mode) {
    switch (mode) {
    case EConsistentTangent:
        return "D";
    case ETransposedTangent:
        return "DT";
    case ESymmetricTangent:
        return "sym";
    case EContinuumTangent:
        return "cont";
    case EFiniteDifferenceTangent:
        return "fd";
    }
    return "unknown";
}

void TPZPlasticStepModifiedCamClay::SetModifiedCamClay(REAL M, REAL lambda, REAL kappa, REAL pt, REAL omega) {
    fModel = EModifiedCamClay;
    fYC.SetUp(M, lambda, kappa, pt, omega);
}

void TPZPlasticStepModifiedCamClay::SetFiniteDifferenceStep(REAL h) {
    if (!(h > 0.) || !std::isfinite(h)) {
        std::cout << __PRETTY_FUNCTION__ << " invalid step h = " << h << " (it must be positive and finite)"
                  << std::endl;
        DebugStop();
    }
    fFDStep = h;
}

void TPZPlasticStepModifiedCamClay::SetLinearElastic(REAL E, REAL nu) {
    fModel = ELinearElastic;
    fER.SetEngineeringData(E, nu);
}

REAL TPZPlasticStepModifiedCamClay::SpecificVolume() const {
    return (fN.fmatprop.size() > 0) ? fN.fmatprop[0] : fV0Default;
}

void TPZPlasticStepModifiedCamClay::SetSpecificVolume(REAL v0) {
    if (fN.fmatprop.size() < 1) fN.fmatprop.Resize(1, 0.);
    fN.fmatprop[0] = v0;
}

void TPZPlasticStepModifiedCamClay::ElasticOperator(REAL K, REAL G, TPZFMatrix<REAL> &C) {
    C.Redim(6, 6);
    const int diag[3] = {_XX_, _YY_, _ZZ_};
    for (int i : diag) {
        for (int j : diag) C(i, j) = K - 2. * G / 3.;
        C(i, i) = K + 4. * G / 3.;
    }
    C(_XY_, _XY_) = G;
    C(_XZ_, _XZ_) = G;
    C(_YZ_, _YZ_) = G;
}

void TPZPlasticStepModifiedCamClay::EigenSystem(const TPZTensor<REAL> &sigma, TPZManVector<REAL, 3> &eigval,
                                                TPZManVector<TPZManVector<REAL, 3>, 3> &eigvec) {
    eigval.Resize(3);
    eigvec.Resize(3);
    REAL scale = 0.;
    for (int i = 0; i < 6; ++i) scale = std::max(scale, std::fabs(sigma[i]));
    bool ok = true;
    if (scale > 0.) {
        TPZTensor<REAL>::TPZDecomposed dec;
        sigma.EigenSystem(dec);
        for (int i = 0; i < 3; ++i) {
            eigval[i] = dec.fEigenvalues[i];
            eigvec[i] = dec.fEigenvectors[i];
            if (eigvec[i].size() != 3) ok = false;
        }
        // orthonormality and reconstruction checks
        for (int i = 0; ok && i < 3; ++i) {
            for (int j = i; j < 3; ++j) {
                REAL dot = 0.;
                for (int k = 0; k < 3; ++k) dot += eigvec[i][k] * eigvec[j][k];
                if (std::fabs(dot - (i == j ? 1. : 0.)) > 1.e-13) ok = false;
            }
        }
        for (int r = 0; ok && r < 3; ++r) {
            for (int c = r; c < 3; ++c) {
                REAL val = 0.;
                for (int i = 0; i < 3; ++i) val += eigval[i] * eigvec[i][r] * eigvec[i][c];
                if (std::fabs(val - sigma(r, c)) > 1.e-13 * scale) ok = false;
            }
        }
        for (int i = 0; ok && i < 2; ++i)
            if (eigval[i] < eigval[i + 1]) ok = false;
    } else {
        ok = false;
    }
    if (ok) return;
    // fall back: cyclic Jacobi rotations, robust for repeated eigenvalues
    REAL A[3][3], V[3][3];
    for (int r = 0; r < 3; ++r)
        for (int c = 0; c < 3; ++c) A[r][c] = sigma(r, c);
    Jacobi3(A, V);
    int order[3] = {0, 1, 2};
    std::stable_sort(order, order + 3, [&A](int a, int b) { return A[a][a] > A[b][b]; });
    for (int i = 0; i < 3; ++i) {
        eigval[i] = A[order[i]][order[i]];
        eigvec[i].Resize(3);
        for (int k = 0; k < 3; ++k) eigvec[i][k] = V[k][order[i]];
    }
}

void TPZPlasticStepModifiedCamClay::TrialStress(const TPZTensor<REAL> &deps, const TPZTensor<REAL> &sigman, REAL v0,
                                                TPZTensor<REAL> &sigtr, REAL &Ktr, REAL &G) const {
    const REAL pn = sigman.I1() / 3.;
    const REAL dev = deps.I1();
    const REAL kappa = fYC.Kappa();
    const bool linear = fYC.VolumetricLaw() == TPZYCModifiedCamClayRHW::ELinear;
    const REAL Kn = linear ? fYC.K0() : -v0 * pn / kappa;
    G = (fShear == EConstantG) ? fG : 3. * (1. - 2. * fNu) / (2. * (1. + fNu)) * Kn;
    REAL ptr;
    if (linear) {
        ptr = pn + fYC.K0() * dev;
        Ktr = fYC.K0();
    } else {
        ptr = pn * std::exp(-v0 * dev / kappa);
        Ktr = -v0 * ptr / kappa;
    }
    // sigma_tr = s_n + 2 G de + p_tr I (engineering shear strains: sigma_xy += G gamma_xy)
    sigtr = sigman;
    sigtr.Add(deps, 2. * G);
    sigtr.XY() -= G * deps.XY();
    sigtr.XZ() -= G * deps.XZ();
    sigtr.YZ() -= G * deps.YZ();
    const REAL shift = ptr - pn - 2. * G * dev / 3.;
    sigtr.XX() += shift;
    sigtr.YY() += shift;
    sigtr.ZZ() += shift;
}

void TPZPlasticStepModifiedCamClay::ComputedDep(const TPZFMatrix<REAL> &Dproj,
                                                const TPZVec<TPZManVector<REAL, 3>> &eigvec, REAL Ktr, REAL G,
                                                REAL ratio, TPZFMatrix<REAL> &Dep) {
    TPZFNMatrix<36, REAL> C(6, 6, 0.);
    ElasticOperator(Ktr, G, C);
    // projection vectors: v_ii (stress type) and v_jj (strain type, engineering shear)
    REAL vsig[3][6], veps[3][6];
    for (int i = 0; i < 3; ++i) {
        SymVoigt(eigvec[i], eigvec[i], vsig[i]);
        for (int k = 0; k < 6; ++k) veps[i][k] = vsig[i][k];
        veps[i][_XY_] *= 2.;
        veps[i][_XZ_] *= 2.;
        veps[i][_YZ_] *= 2.;
    }
    // rotational correction: kappa_ij s_ij s_ij^T for i < j, eq. (9), with
    // kappa_ij = (sigma_i - sigma_j)/(eps_tr,i - eps_tr,j) = 2 G rho/rho_tr for every pair
    REAL sij[3][6];
    const REAL kap = 2. * G * ratio;
    const int pairs[3][2] = {{0, 1}, {0, 2}, {1, 2}};
    for (int ip = 0; ip < 3; ++ip) SymVoigt(eigvec[pairs[ip][0]], eigvec[pairs[ip][1]], sij[ip]);
    Dep.Redim(6, 6);
    // column icol of the operator = response to the unit engineering strain e_icol
    for (int icol = 0; icol < 6; ++icol) {
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                // (v_jj^eps)^T C e_icol
                REAL proj = 0.;
                for (int k = 0; k < 6; ++k) proj += veps[j][k] * C.GetVal(k, icol);
                const REAL fac = Dproj.GetVal(i, j) * proj;
                for (int l = 0; l < 6; ++l) Dep(l, icol) += fac * vsig[i][l];
            }
        }
        for (int ip = 0; ip < 3; ++ip) {
            const REAL fac = 2. * kap * sij[ip][icol];
            for (int l = 0; l < 6; ++l) Dep(l, icol) += fac * sij[ip][l];
        }
    }
}

void TPZPlasticStepModifiedCamClay::ContinuumTangent(const TPZTensor<REAL> &sigma, REAL pc, REAL v0,
                                                     TPZFMatrix<REAL> &Dc) const {
    // moduli at the state (continuum_tangent of camclay_hw.py)
    const REAL p = sigma.I1() / 3.;
    const bool linear = fYC.VolumetricLaw() == TPZYCModifiedCamClayRHW::ELinear;
    const REAL K = linear ? fYC.K0() : -v0 * p / fYC.Kappa();
    const REAL G = (fShear == EConstantG) ? fG : 3. * (1. - 2. * fNu) / (2. * (1. + fNu)) * K;
    TPZFNMatrix<36, REAL> C(6, 6, 0.);
    ElasticOperator(K, G, C);
    // gradient of the yield function in engineering form (shear components doubled)
    const REAL omega = fYC.Omega(), pt = fYC.Pt(), M2 = fYC.M() * fYC.M();
    const REAL a = (pc + pt) / (1. + omega);
    const REAL pbar = p - pt + a;
    const REAL b = (pbar >= 0.) ? 1. : omega;
    const REAL b2 = b * b;
    REAL n[6];
    for (int k = 0; k < 6; ++k) {
        const bool diag = (k == _XX_ || k == _YY_ || k == _ZZ_);
        const REAL sk = diag ? sigma[k] - p : sigma[k];
        n[k] = (diag ? 2. * pbar / (3. * b2) : 0.) + 3. * sk / M2;
        if (!diag) n[k] *= 2.;
    }
    // hardening: dPhi/da and da/dgamma (dal = -2 dgamma pbar/b^2)
    const REAL H = v0 * pc / ((fYC.Lambda() - fYC.Kappa()) * (1. + omega));
    const REAL dphida = 2. * pbar / b2 - 2. * a;
    const REAL dadg = -H * 2. * pbar / b2;
    REAL Cn[6], nC[6], nCn = 0.;
    for (int l = 0; l < 6; ++l) {
        Cn[l] = 0.;
        nC[l] = 0.;
        for (int k = 0; k < 6; ++k) {
            Cn[l] += C.GetVal(l, k) * n[k];
            nC[l] += n[k] * C.GetVal(k, l);
        }
    }
    for (int k = 0; k < 6; ++k) nCn += n[k] * Cn[k];
    const REAL den = nCn - dphida * dadg;
    Dc.Redim(6, 6);
    for (int l = 0; l < 6; ++l)
        for (int k = 0; k < 6; ++k) Dc(l, k) = C.GetVal(l, k) - Cn[l] * nC[k] / den;
}

bool TPZPlasticStepModifiedCamClay::FiniteDifferenceTangent(const TPZTensor<REAL> &epsTotal,
                                                            const TPZTensor<REAL> &sigman, TPZFMatrix<REAL> &Dfd) {
    // the object is restored after the 12 updates
    const TPZPlasticState<REAL> N(fN);
    const TPZElasticResponse ER(fER);
    const TPZTensor<REAL> sig(fSigma);
    const int nit = fLastNewtonIterations;
    const bool failed = fFailed;
    const REAL h = fFDStep;
    Dfd.Redim(6, 6);
    bool ok = true;
    for (int k = 0; k < 6 && ok; ++k) {
        TPZTensor<REAL> sp, sm;
        for (int sign = 0; sign < 2 && ok; ++sign) {
            TPZTensor<REAL> eps(epsTotal);
            eps[k] += (sign == 0) ? h : -h;
            TPZTensor<REAL> &s = (sign == 0) ? sp : sm;
            s = sigman;
            fN = N;
            StressUpdate(eps, s, nullptr);
            if (fFailed) ok = false;
        }
        if (!ok) break;
        for (int l = 0; l < 6; ++l) Dfd(l, k) = (sp[l] - sm[l]) / (2. * h);
    }
    fN = N;
    fER = ER;
    fSigma = sig;
    fLastNewtonIterations = nit;
    fFailed = failed;
    return ok;
}

void TPZPlasticStepModifiedCamClay::ApplyStrainComputeSigma(const TPZTensor<REAL> &epsTotal, TPZTensor<REAL> &sigma,
                                                            TPZFMatrix<REAL> *tangent) {
    if (!tangent || fModel == ELinearElastic || fTangentMode == EConsistentTangent) {
        StressUpdate(epsTotal, sigma, tangent);
        return;
    }
    // converged state, stress and elastic response, for the finite differences and the failure case
    const TPZPlasticState<REAL> N0(fN);
    const TPZElasticResponse ER0(fER);
    const TPZTensor<REAL> sigma0(fSigma);
    const TPZTensor<REAL> sigman(sigma);
    StressUpdate(epsTotal, sigma, tangent);
    if (fFailed) return;
    TPZFMatrix<REAL> &D = *tangent;
    switch (fTangentMode) {
    case EConsistentTangent:
        break;
    case ETransposedTangent: {
        TPZFNMatrix<36, REAL> Dt;
        D.Transpose(&Dt);
        D = Dt;
        break;
    }
    case ESymmetricTangent: {
        TPZFNMatrix<36, REAL> Ds(6, 6, 0.);
        for (int l = 0; l < 6; ++l)
            for (int k = 0; k < 6; ++k) Ds(l, k) = 0.5 * (D.GetVal(l, k) + D.GetVal(k, l));
        D = Ds;
        break;
    }
    case EContinuumTangent:
        if (fN.m_m_type != 0) ContinuumTangent(sigma, fN.m_hardening, SpecificVolume(), D);
        break;
    case EFiniteDifferenceTangent: {
        // the extra updates start from the converged state of the point
        const TPZPlasticState<REAL> N1(fN);
        fN = N0;
        TPZFNMatrix<36, REAL> Dfd(6, 6, 0.);
        const bool ok = FiniteDifferenceTangent(epsTotal, sigman, Dfd);
        fN = N1;
        if (ok) {
            D = Dfd;
            break;
        }
        // failure of an extra update: the call fails as a whole (state, stress and elastic response unchanged)
        fN = N0;
        fER = ER0;
        fSigma = sigma0;
        fFailed = true;
        sigma = sigman;
        TPZTensor<REAL> sigtr;
        REAL Ktr, G;
        TrialStress(epsTotal - N0.m_eps_t, sigman, SpecificVolume(), sigtr, Ktr, G);
        ElasticOperator(Ktr, G, D);
        break;
    }
    }
}

void TPZPlasticStepModifiedCamClay::StressUpdate(const TPZTensor<REAL> &epsTotal, TPZTensor<REAL> &sigma,
                                                 TPZFMatrix<REAL> *tangent) {
    fFailed = false;
    fLastNewtonIterations = 0;
    const TPZTensor<REAL> sigman(sigma);
    const TPZTensor<REAL> deps = epsTotal - fN.m_eps_t;
    if (tangent) tangent->Redim(6, 6);

    if (fModel == ELinearElastic) {
        TPZFNMatrix<36, REAL> De(6, 6, 0.);
        fER.De(De);
        for (int l = 0; l < 6; ++l) {
            REAL val = 0.;
            for (int k = 0; k < 6; ++k) val += De(l, k) * deps[k];
            sigma[l] = sigman[l] + val;
        }
        if (tangent) *tangent = De;
        fN.m_eps_t = epsTotal;
        fN.m_m_type = 0;
        fSigma = sigma;
        return;
    }

    const REAL sq3 = std::sqrt(3.);
    const REAL v0 = SpecificVolume();
    const REAL pcn = fN.m_hardening;

    // 1. elastic trial state
    TPZTensor<REAL> sigtr;
    REAL Ktr, G;
    TrialStress(deps, sigman, v0, sigtr, Ktr, G);
    // equivalent pair (E, nu) of (K_tr, G), stored in fER when the step succeeds
    const REAL Eq = 9. * Ktr * G / (3. * Ktr + G);
    const REAL nuq = (3. * Ktr - 2. * G) / (2. * (3. * Ktr + G));
    const bool validER = std::isfinite(Eq) && std::isfinite(nuq) && Eq > 0.;

    // 2. spectral decomposition of the trial stress
    TPZManVector<REAL, 3> strhw(3, 0.);
    TPZManVector<TPZManVector<REAL, 3>, 3> vecs(3);
    EigenSystem(sigtr, strhw, vecs);
    const REAL ptr = (strhw[0] + strhw[1] + strhw[2]) / 3.;
    TPZManVector<REAL, 3> s(3, 0.);
    REAL rhotr = 0.;
    for (int i = 0; i < 3; ++i) {
        s[i] = strhw[i] - ptr;
        rhotr += s[i] * s[i];
    }
    rhotr = std::sqrt(rhotr);

    TPZYCModifiedCamClayRHW::TTrial trial;
    trial.fXiTr = sq3 * ptr;
    trial.fRhoTr = rhotr;
    trial.fG = G;
    trial.fV0 = v0;
    trial.fPcn = pcn;
    REAL an, H, pc;
    fYC.Hardening(pcn, 0., v0, an, H, pc);
    const REAL omega = fYC.Omega();
    const REAL btr = fYC.BFromP(ptr, an);

    // 3. elastic step
    if (fYC.PhiCC(ptr, rhotr, an, btr) <= 1.e-11 * an * an) {
        sigma = sigtr;
        if (tangent) ElasticOperator(Ktr, G, *tangent);
        if (validER) fER.SetEngineeringData(Eq, nuq);
        fN.m_eps_t = epsTotal;
        fN.m_m_type = 0;
        fSigma = sigma;
        return;
    }

    // 4. local Newton iterations, with the region check for omega != 1
    TPZManVector<REAL, 3> blist;
    if (omega == 1.) {
        blist.Resize(1, 1.);
    } else {
        blist.Resize(2);
        blist[0] = btr;
        blist[1] = (btr == 1.) ? omega : 1.;
    }
    bool found = false;
    TPZManVector<REAL, 4> X(4, 0.);
    REAL pbar = 0.;
    for (int ib = 0; ib < blist.size() && !found; ++ib) {
        trial.fB = blist[ib];
        int niter = 0;
        if (!fYC.ProjectHW(trial, X, niter)) continue;
        REAL a;
        fYC.Hardening(pcn, X[2], v0, a, H, pc);
        pbar = X[0] / sq3 - fYC.Pt() + a;
        const bool region = (omega == 1.) || (trial.fB == 1. && pbar >= 0.) || (trial.fB != 1. && pbar < 0.);
        if (region && X[3] >= -1.e-14) {
            found = true;
            fLastNewtonIterations = niter;
        }
    }
    if (!found) {
        fFailed = true;
        sigma = sigman;
        if (tangent) ElasticOperator(Ktr, G, *tangent);
        return;
    }

    // 5. projected principal stresses
    const bool isotropic = !(rhotr > 1.e-14 * std::max(1., std::fabs(ptr)));
    TPZManVector<REAL, 3> n(3, 0.), sigproj(3, 0.);
    for (int i = 0; i < 3; ++i) {
        n[i] = isotropic ? 0. : s[i] / rhotr;
        sigproj[i] = X[0] / sq3 + X[1] * n[i];
    }
    sigma.Zero();
    for (int i = 0; i < 3; ++i) {
        for (int r = 0; r < 3; ++r) {
            for (int c = r; c < 3; ++c) sigma(r, c) += sigproj[i] * vecs[i][r] * vecs[i][c];
        }
    }

    // 6-7. Jacobian of the projection and consistent tangent
    if (tangent) {
        TPZFNMatrix<9, REAL> Dproj(3, 3, 0.);
        if (!fYC.GradProjection(X, trial, n, isotropic, Dproj)) {
            fFailed = true;
            sigma = sigman;
            ElasticOperator(Ktr, G, *tangent);
            return;
        }
        // rho/rho_tr from the second local residual, rho (1 + 6 G dgamma/M^2) = rho_tr (exact also for
        // isotropic trial states, where the quotient of the principal differences is undefined)
        const REAL M = fYC.M();
        const REAL ratio = 1. / (1. + 6. * G * X[3] / (M * M));
        ComputedDep(Dproj, vecs, Ktr, G, ratio, *tangent);
    }

    if (validER) fER.SetEngineeringData(Eq, nuq);
    fYC.Hardening(pcn, X[2], v0, an, H, pc);
    fN.m_hardening = pc;
    fN.m_eps_t = epsTotal;
    fN.m_m_type = (pbar < 0.) ? 1 : 2;
    fSigma = sigma;
}

void TPZPlasticStepModifiedCamClay::ApplyStrainComputeDep(const TPZTensor<REAL> &epsTotal, TPZTensor<REAL> &sigma,
                                                          TPZFMatrix<REAL> &Dep) {
    ApplyStrainComputeSigma(epsTotal, sigma, &Dep);
}

void TPZPlasticStepModifiedCamClay::ApplyStrain(const TPZTensor<REAL> &epsTotal) {
    std::cout << __PRETTY_FUNCTION__ << ": the model is incremental, use ApplyStrainComputeSigma with the "
              << "converged stress" << std::endl;
    DebugStop();
}

void TPZPlasticStepModifiedCamClay::ApplyLoad(const TPZTensor<REAL> &sigma, TPZTensor<REAL> &epsTotal) {
    std::cout << __PRETTY_FUNCTION__ << " not implemented" << std::endl;
    DebugStop();
}

void TPZPlasticStepModifiedCamClay::Phi(const TPZTensor<REAL> &epsTotal, TPZVec<REAL> &phi) const {
    TPZManVector<REAL, 3> eigval(3);
    TPZManVector<TPZManVector<REAL, 3>, 3> eigvec(3);
    EigenSystem(fSigma, eigval, eigvec);
    fYC.YieldFunction(eigval, fN.m_hardening, phi);
}

void TPZPlasticStepModifiedCamClay::Print(std::ostream &out) const {
    out << Name() << "\n model = " << (fModel == EModifiedCamClay ? "Modified Cam-Clay" : "linear elastic")
        << "\n shear = " << (fShear == EConstantG ? "constant G = " : "Poisson nu = ")
        << (fShear == EConstantG ? fG : fNu) << "\n default v0 = " << fV0Default
        << "\n tangent mode = " << TangentModeName(fTangentMode) << " (finite difference step " << fFDStep << ")\n";
    fYC.Print(out);
    fER.Print(out);
    fN.Print(out);
}

int TPZPlasticStepModifiedCamClay::ClassId() const {
    return Hash("TPZPlasticStepModifiedCamClay") ^ TPZPlasticBase::ClassId() << 1;
}

void TPZPlasticStepModifiedCamClay::Write(TPZStream &buf, int withclassid) const {
    fYC.Write(buf, withclassid);
    fER.Write(buf, withclassid);
    int model = fModel, shear = fShear, mode = fTangentMode;
    buf.Write(&model);
    buf.Write(&shear);
    buf.Write(&fG);
    buf.Write(&fNu);
    buf.Write(&fV0Default);
    buf.Write(&mode);
    buf.Write(&fFDStep);
    fN.Write(buf, withclassid);
    fSigma.Write(buf, withclassid);
}

void TPZPlasticStepModifiedCamClay::Read(TPZStream &buf, void *context) {
    fYC.Read(buf, context);
    fER.Read(buf, context);
    int model, shear, mode;
    buf.Read(&model);
    buf.Read(&shear);
    buf.Read(&fG);
    buf.Read(&fNu);
    buf.Read(&fV0Default);
    buf.Read(&mode);
    buf.Read(&fFDStep);
    fModel = EModel(model);
    fShear = EShear(shear);
    fTangentMode = ETangentMode(mode);
    fN.Read(buf, context);
    fSigma.Read(buf, context);
}
