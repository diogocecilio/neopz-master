// TPZModifiedCamClay.cpp
//
// Porte para o NeoPZ de modified_cam_clay.py (return mapping implícito do Cam-Clay
// modificado, de Souza Neto, Perić & Owen 2008, Seção 10.1). Ver TPZModifiedCamClay.h.

#include "TPZModifiedCamClay.h"

#include <cmath>
#include <sstream>
#include <iomanip>

#include "TPZStream.h"
#include "pzerror.h"

namespace {

// Voigt {xx, xy, xz, yy, yz, zz}
const REAL kM[6] = {1., 0., 0., 1., 0., 1.};          // identidade
const REAL kShearFac[6] = {1., .5, .5, 1., .5, 1.};   // engenharia -> componentes tensoriais
const REAL kSqrt23 = std::sqrt(2. / 3.);

/// Norma de um tensor simétrico guardado em Voigt (componentes tensoriais)
inline REAL TNorm(const REAL t[6]) {
    return std::sqrt(t[_XX_] * t[_XX_] + t[_YY_] * t[_YY_] + t[_ZZ_] * t[_ZZ_] +
                     2. * (t[_XY_] * t[_XY_] + t[_XZ_] * t[_XZ_] + t[_YZ_] * t[_YZ_]));
}

/// I_DEV ε: deformação em Voigt (engenharia) -> desviador (componentes tensoriais)
inline void IDev(const REAL e[6], REAL out[6]) {
    const REAL tr = e[_XX_] + e[_YY_] + e[_ZZ_];
    for (int i = 0; i < 6; i++) out[i] = kShearFac[i] * e[i] - kM[i] * tr / 3.;
}

/// Entrada (i, j) da matriz I_DEV = diag(SHEAR_FAC) - m mᵀ/3
inline REAL IDevIJ(int i, int j) {
    return (i == j ? kShearFac[i] : 0.) - kM[i] * kM[j] / 3.;
}

inline void ToArray(const TPZTensor<REAL> &t, REAL a[6]) {
    for (int i = 0; i < 6; i++) a[i] = t[i];
}

inline void FromArray(const REAL a[6], TPZTensor<REAL> &t) {
    for (int i = 0; i < 6; i++) t[i] = a[i];
}

/// Resolve A x = b (n x n, eliminação de Gauss com pivotamento parcial). Retorna false se singular.
bool SolveDense(int n, REAL *A, REAL *b) {
    for (int k = 0; k < n; k++) {
        int piv = k;
        REAL vmax = std::fabs(A[k * n + k]);
        for (int i = k + 1; i < n; i++) {
            if (std::fabs(A[i * n + k]) > vmax) { vmax = std::fabs(A[i * n + k]); piv = i; }
        }
        if (vmax == 0.) return false;
        if (piv != k) {
            for (int j = 0; j < n; j++) std::swap(A[k * n + j], A[piv * n + j]);
            std::swap(b[k], b[piv]);
        }
        for (int i = k + 1; i < n; i++) {
            const REAL fac = A[i * n + k] / A[k * n + k];
            for (int j = k; j < n; j++) A[i * n + j] -= fac * A[k * n + j];
            b[i] -= fac * b[k];
        }
    }
    for (int k = n - 1; k >= 0; k--) {
        REAL s = b[k];
        for (int j = k + 1; j < n; j++) s -= A[k * n + j] * b[j];
        b[k] = s / A[k * n + k];
    }
    return true;
}

} // namespace


// =====================================================================================
// TPZYCModifiedCamClay
// =====================================================================================

void TPZYCModifiedCamClay::SetUp(REAL M, REAL lambda, REAL kappa, REAL v0, REAL pc0, REAL pt, REAL beta) {
    fM = M;
    fLambda = lambda;
    fKappa = kappa;
    fV0 = v0;
    fPc0 = pc0;
    fPt = pt;
    fBeta = beta;
}

REAL TPZYCModifiedCamClay::Pc(REAL alpha) const {
    return fPc0 * std::exp(fV0 * alpha / (fLambda - fKappa));
}

void TPZYCModifiedCamClay::Hardening(REAL alpha, REAL &a, REAL &H) const {
    const REAL c = fV0 / (fLambda - fKappa);
    const REAL pc = fPc0 * std::exp(c * alpha);
    a = (pc + fPt) / (1. + fBeta);
    H = c * pc / (1. + fBeta);
}

REAL TPZYCModifiedCamClay::YieldValue(REAL p, REAL q, REAL a) const {
    const REAL b = BranchB(p, a);
    const REAL pb = p - fPt + a;
    return pb * pb / (b * b) + (q / fM) * (q / fM) - a * a;
}

void TPZYCModifiedCamClay::YieldFunction(const TPZVec<STATE> &sigma, STATE kprev, TPZVec<STATE> &yield) const {
    const REAL p = (sigma[0] + sigma[1] + sigma[2]) / 3.;
    const REAL q = std::sqrt(0.5 * ((sigma[0] - sigma[1]) * (sigma[0] - sigma[1]) +
                                    (sigma[1] - sigma[2]) * (sigma[1] - sigma[2]) +
                                    (sigma[2] - sigma[0]) * (sigma[2] - sigma[0])));
    REAL a, H;
    Hardening(kprev, a, H);
    yield.Resize(NYield);
    yield[0] = YieldValue(p, q, a);
}

int TPZYCModifiedCamClay::ClassId() const {
    return Hash("TPZYCModifiedCamClay");
}

void TPZYCModifiedCamClay::Read(TPZStream &buf, void *context) {
    buf.Read(&fM);
    buf.Read(&fLambda);
    buf.Read(&fKappa);
    buf.Read(&fV0);
    buf.Read(&fPc0);
    buf.Read(&fPt);
    buf.Read(&fBeta);
}

void TPZYCModifiedCamClay::Write(TPZStream &buf, int withclassid) const {
    buf.Write(&fM);
    buf.Write(&fLambda);
    buf.Write(&fKappa);
    buf.Write(&fV0);
    buf.Write(&fPc0);
    buf.Write(&fPt);
    buf.Write(&fBeta);
}

void TPZYCModifiedCamClay::Print(std::ostream &out) const {
    out << "\nTPZYCModifiedCamClay:"
        << "\n M = " << fM << "  lambda = " << fLambda << "  kappa = " << fKappa
        << "\n v0 = " << fV0 << "  pc0 = " << fPc0 << "  pt = " << fPt << "  beta = " << fBeta;
}


// =====================================================================================
// TPZModifiedCamClay
// =====================================================================================

TPZModifiedCamClay::TPZModifiedCamClay() : TPZPlasticBase(), fN() {
    // mesmos valores padrão de mcc_parameters (Tabela 8.1 do RS2)
    SetUp(1.2, 0.077, 0.0066, 1.788, -1., 200., 200.);
}

TPZModifiedCamClay::TPZModifiedCamClay(const TPZModifiedCamClay &other) = default;

TPZModifiedCamClay &TPZModifiedCamClay::operator=(const TPZModifiedCamClay &other) = default;

TPZModifiedCamClay::~TPZModifiedCamClay() = default;

void TPZModifiedCamClay::SetUp(REAL M, REAL lambda, REAL kappa, REAL N, REAL v0, REAL pc0, REAL p0,
                               REAL pt, REAL beta, EElasticity elasticity, EShear shear, REAL G, REAL nu) {
    if (lambda <= kappa || kappa <= 0. || M <= 0. || pc0 <= 0. || p0 <= 0. || beta <= 0.) {
        PZError << "TPZModifiedCamClay::SetUp: parâmetros inválidos (exige lambda > kappa > 0, M > 0, "
                   "pc0 > 0, p0 > 0, beta > 0)\n";
        DebugStop();
    }
    if (v0 <= 0.) {
        v0 = N - lambda * std::log(pc0) + kappa * std::log(pc0 / p0);
    }
    fYC.SetUp(M, lambda, kappa, v0, pc0, pt, beta);
    fNIso = N;
    fP0 = p0;
    fG = G;
    fNu = nu;
    fElasticity = elasticity;
    fShear = shear;
    fK0User = 0.;
    fSigma0.Zero();
    fSigma0[_XX_] = fSigma0[_YY_] = fSigma0[_ZZ_] = -p0;
    UpdateDerived();
    fN.CleanUp();
}

void TPZModifiedCamClay::SetInitialStress(const TPZTensor<REAL> &sigma0) {
    fSigma0 = sigma0;
    UpdateDerived();
}

void TPZModifiedCamClay::SetLinearBulkModulus(REAL K0) {
    if (K0 > 0. && fElasticity != ELinear) {
        PZError << "TPZModifiedCamClay::SetLinearBulkModulus: só faz sentido com elasticidade ELinear\n";
        DebugStop();
    }
    fK0User = K0;
    UpdateDerived();
}

void TPZModifiedCamClay::UpdateDerived() {
    fPIni = (fSigma0[_XX_] + fSigma0[_YY_] + fSigma0[_ZZ_]) / 3.;
    if (fPIni >= 0.) {
        PZError << "TPZModifiedCamClay: a tensão inicial deve ser de compressão (p_ini < 0)\n";
        DebugStop();
    }
    fK0 = (fElasticity == ELinear && fK0User > 0.) ? fK0User : -fYC.V0() * fPIni / fYC.Kappa();
    fGFac = 3. * (1. - 2. * fNu) / (2. * (1. + fNu));
    fG0 = (fShear == EConstantG) ? fG : fGFac * fK0;
    for (int i = 0; i < 6; i++) fE0[i] = (fSigma0[i] - fPIni * kM[i]) / (2. * fG0);
}

void TPZModifiedCamClay::Invariants(const TPZTensor<REAL> &sigma, REAL &p, REAL &q) {
    REAL s[6];
    ToArray(sigma, s);
    p = (s[_XX_] + s[_YY_] + s[_ZZ_]) / 3.;
    for (int i = 0; i < 6; i++) s[i] -= p * kM[i];
    q = std::sqrt(1.5) * TNorm(s);
}

void TPZModifiedCamClay::Pressure(REAL eev, REAL &p, REAL &K, REAL &dK) const {
    if (fElasticity == ELinear) {
        p = fPIni + fK0 * eev;
        K = fK0;
        dK = 0.;
        return;
    }
    const REAL c = fYC.V0() / fYC.Kappa();
    p = fPIni * std::exp(-c * eev);
    K = -c * p;
    dK = -c * K;
}

void TPZModifiedCamClay::ShearModulus(REAL eev, REAL &G, REAL &dG) const {
    if (fShear == EConstantG) {
        G = fG;
        dG = 0.;
        return;
    }
    if (fElasticity == ELinear) {
        G = fGFac * fK0;
        dG = 0.;
        return;
    }
    REAL p, K, dK;
    Pressure(eev, p, K, dK);
    G = fGFac * K;
    dG = fGFac * dK;
}

void TPZModifiedCamClay::ReturnMapping(const TPZTensor<REAL> &epsTensor, const TPZTensor<REAL> &epspTensor,
                                       REAL alpha_n, TResult &res, const TPZTensor<REAL> *sig_n,
                                       const TPZTensor<REAL> *eps_n) const {
    const REAL M = fYC.M();
    const REAL M2 = M * M;
    const REAL pt = fYC.Pt();
    const REAL beta = fYC.Beta();
    const bool hypo = (fShear == EHypoNu);

    REAL eps[6], epsp_n[6], e0[6];
    ToArray(epsTensor, eps);
    ToArray(epspTensor, epsp_n);
    ToArray(fE0, e0);

    // ---- estado tentativa elástico (10.10)
    REAL epse_tr[6], ehat[6];
    for (int i = 0; i < 6; i++) epse_tr[i] = eps[i] - epsp_n[i];
    const REAL xv = epse_tr[_XX_] + epse_tr[_YY_] + epse_tr[_ZZ_]; // ε_v^{e,trial}
    REAL Gh = 0.;
    if (hypo) {
        // G congelado no início do passo; s^trial = s_n + 2 G_n Δe  ->  ê = s_n/(2G_n) + Δe
        REAL sn[6], en[6] = {0., 0., 0., 0., 0., 0.};
        ToArray(sig_n ? *sig_n : fSigma0, sn);
        if (sig_n && eps_n) ToArray(*eps_n, en);
        const REAL p_n = (sn[_XX_] + sn[_YY_] + sn[_ZZ_]) / 3.;
        const REAL Kn = (fElasticity == ELinear) ? fK0 : -(fYC.V0() / fYC.Kappa()) * p_n;
        Gh = fGFac * Kn;
        REAL deps[6], de[6];
        for (int i = 0; i < 6; i++) deps[i] = eps[i] - en[i];
        IDev(deps, de);
        for (int i = 0; i < 6; i++) ehat[i] = (sn[i] - p_n * kM[i]) / (2. * Gh) + de[i];
    } else {
        IDev(epse_tr, ehat);                             // ε_d^{e,trial}
        for (int i = 0; i < 6; i++) ehat[i] += e0[i];    // (+ desviador inicial)
    }
    const REAL nrm = TNorm(ehat);
    const REAL eq_tr = kSqrt23 * nrm;
    REAL n[6];
    for (int i = 0; i < 6; i++) n[i] = (nrm > 1.e-14) ? ehat[i] / nrm : 0.;

    REAL p_tr, K_tr, dK_tr, G_tr, dG_tr;
    Pressure(xv, p_tr, K_tr, dK_tr);
    if (hypo) {
        G_tr = Gh;
        dG_tr = 0.;
    } else {
        ShearModulus(xv, G_tr, dG_tr);
    }
    const REAL q_tr = 3. * G_tr * eq_tr;
    REAL a_n, H_n;
    fYC.Hardening(alpha_n, a_n, H_n);

    res.Dep.Redim(6, 6);
    if (fYC.YieldValue(p_tr, q_tr, a_n) <= fTol * a_n * a_n) {
        // ---- passo elástico
        REAL sig[6];
        for (int i = 0; i < 6; i++) sig[i] = 2. * G_tr * ehat[i] + p_tr * kM[i];
        FromArray(sig, res.stress);
        for (int i = 0; i < 6; i++) {
            for (int j = 0; j < 6; j++) {
                res.Dep(i, j) = K_tr * kM[i] * kM[j] + 2. * G_tr * IDevIJ(i, j) + 2. * dG_tr * ehat[i] * kM[j];
            }
        }
        if (hypo) {
            FromArray(epse_tr, res.elastic_strain);
            res.plastic_strain = epspTensor;
        } else {
            REAL epse[6], epsp[6];
            for (int i = 0; i < 6; i++) {
                epse[i] = (ehat[i] - e0[i] + xv / 3. * kM[i]) / kShearFac[i];
                epsp[i] = eps[i] - epse[i];
            }
            FromArray(epse, res.elastic_strain);
            FromArray(epsp, res.plastic_strain);
        }
        res.alpha = alpha_n;
        res.dgamma = 0.;
        res.plastic = false;
        res.iterations = 0;
        res.b = 1.;
        res.p = p_tr;
        res.q = q_tr;
        res.pc = fYC.Pc(alpha_n);
        return;
    }

    // ---- corretor plástico: Newton-Raphson no sistema reduzido (10.17)
    REAL blist[2] = {1., beta};
    int nb = 2;
    if (beta == 1.) {
        nb = 1;
    } else if (p_tr < pt - a_n) {
        blist[0] = beta;
        blist[1] = 1.;
    }
    bool done = false;
    REAL b = 1., dg = 0., al = alpha_n;
    REAL ee = 0., p = 0., K = 0., dK = 0., G = 0., dG = 0., f = 1., q = 0., a = 0., H = 0., pb = 0.;
    int it = 0;
    for (int ib = 0; ib < nb && !done; ib++) {
        b = blist[ib];
        const REAL b2 = b * b;
        dg = 0.;
        al = alpha_n;
        bool conv = false;
        for (it = 1; it <= fMaxIt; it++) {
            ee = xv + al - alpha_n;                   // ε_v^e_{n+1}
            Pressure(ee, p, K, dK);                   // p(α)  (10.16)
            if (hypo) {
                G = Gh;
                dG = 0.;
            } else {
                ShearModulus(ee, G, dG);
            }
            f = M2 / (M2 + 6. * G * dg);
            q = 3. * G * f * eq_tr;                   // q(Δγ)  (10.15)
            fYC.Hardening(al, a, H);
            pb = p - pt + a;                          // p̄  (10.22)
            const REAL R1 = pb * pb / b2 + (q / M) * (q / M) - a * a;
            const REAL R2 = al - alpha_n + dg * 2. * pb / b2;
            if (std::fabs(R1) <= fTol * a * a && std::fabs(R2) <= fTol) {
                conv = true;
                break;
            }
            // jacobiana (10.20)
            const REAL J00 = -12. * G * f * q * q / (M2 * M2);
            const REAL J01 = 2. * pb / b2 * (K + H) + 2. * q / M2 * (q * f / G) * dG - 2. * a * H;
            const REAL J10 = 2. * pb / b2;
            const REAL J11 = 1. + 2. * dg / b2 * (K + H);
            const REAL det = J00 * J11 - J01 * J10;
            if (det == 0. || !std::isfinite(det)) break;
            const REAL d0 = (-R1 * J11 + R2 * J01) / det;
            const REAL d1 = (-R2 * J00 + R1 * J10) / det;
            dg += d0;
            al += d1;
            if (!std::isfinite(dg) || !std::isfinite(al)) break;
        }
        if (it > fMaxIt) it = fMaxIt;
        if (conv && (beta == 1. || (b == 1. && pb >= 0.) || (b != 1. && pb < 0.))) {
            done = true;
        }
    }
    if (!done || dg < -1.e-12) {
        std::stringstream sout;
        sout << "TPZModifiedCamClay::ReturnMapping: return mapping não convergiu (eps = {";
        sout << std::setprecision(10);
        for (int i = 0; i < 6; i++) sout << eps[i] << (i < 5 ? ", " : "})");
        throw ReturnMappingError(sout.str());
    }
    const REAL b2 = b * b;

    // ---- atualização (10.14), (10.18)
    REAL sig[6], epse[6], epsp[6];
    for (int i = 0; i < 6; i++) sig[i] = 2. * G * f * ehat[i] + p * kM[i];
    if (hypo) {
        for (int i = 0; i < 6; i++) {
            const REAL depsp = ((1. - f) * ehat[i] - (al - alpha_n) / 3. * kM[i]) / kShearFac[i];
            epsp[i] = epsp_n[i] + depsp;
            epse[i] = eps[i] - epsp[i];
        }
    } else {
        for (int i = 0; i < 6; i++) {
            epse[i] = (f * ehat[i] - e0[i] + ee / 3. * kM[i]) / kShearFac[i];
            epsp[i] = eps[i] - epse[i];
        }
    }

    // ---- tangente consistente
    const REAL dqdG = q * f / G;
    const REAL dqdg = -6. * G * f * q / M2;
    const REAL J00 = -12. * G * f * q * q / (M2 * M2);
    const REAL J01 = 2. * pb / b2 * (K + H) + 2. * q / M2 * dqdG * dG - 2. * a * H;
    const REAL J10 = 2. * pb / b2;
    const REAL J11 = 1. + 2. * dg / b2 * (K + H);
    const REAL det = J00 * J11 - J01 * J10;
    const REAL dR1dx = 2. * pb / b2 * K + 2. * q / M2 * dqdG * dG;  // ∂R1/∂ε_v^trial
    const REAL dR1deq = 2. * q / M2 * 3. * G * f;                   // ∂R1/∂ε_q^trial
    const REAL dR2dx = 2. * dg / b2 * K;                            // ∂R2/∂ε_v^trial
    REAL deq_de[6], ddg_de[6], dal_de[6], dq_de[6], dp_de[6];
    for (int j = 0; j < 6; j++) {
        deq_de[j] = kSqrt23 * n[j];
        const REAL rhs0 = dR1dx * kM[j] + dR1deq * deq_de[j];
        const REAL rhs1 = dR2dx * kM[j];
        ddg_de[j] = -(J11 * rhs0 - J01 * rhs1) / det;
        dal_de[j] = -(-J10 * rhs0 + J00 * rhs1) / det;
        dq_de[j] = dqdg * ddg_de[j] + dqdG * dG * (kM[j] + dal_de[j]) + 3. * G * f * deq_de[j];
        dp_de[j] = K * (kM[j] + dal_de[j]);
    }
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            res.Dep(i, j) = 2. * G * f * (IDevIJ(i, j) - n[i] * n[j]) + kSqrt23 * n[i] * dq_de[j] + kM[i] * dp_de[j];
        }
    }
    FromArray(sig, res.stress);
    FromArray(epse, res.elastic_strain);
    FromArray(epsp, res.plastic_strain);
    res.alpha = al;
    res.dgamma = dg;
    res.plastic = true;
    res.iterations = it;
    res.b = b;
    res.p = p;
    res.q = q;
    res.pc = fYC.Pc(al);
}

// ------------------------------------------------------------------------- TPZPlasticBase

int TPZModifiedCamClay::ClassId() const {
    return Hash("TPZModifiedCamClay") ^ TPZPlasticBase::ClassId() << 1;
}

void TPZModifiedCamClay::Write(TPZStream &buf, int withclassid) const {
    fYC.Write(buf, withclassid);
    buf.Write(&fNIso);
    buf.Write(&fP0);
    buf.Write(&fG);
    buf.Write(&fNu);
    int elasticity = fElasticity, shear = fShear;
    buf.Write(&elasticity);
    buf.Write(&shear);
    fSigma0.Write(buf, withclassid);
    buf.Write(&fK0User);
    buf.Write(&fTol);
    buf.Write(&fMaxIt);
    int mapping = fStrengthMapping;
    buf.Write(&mapping);
    buf.Write(&fReductionFactor);
    fN.Write(buf, withclassid);
}

void TPZModifiedCamClay::Read(TPZStream &buf, void *context) {
    fYC.Read(buf, context);
    buf.Read(&fNIso);
    buf.Read(&fP0);
    buf.Read(&fG);
    buf.Read(&fNu);
    int elasticity, shear;
    buf.Read(&elasticity);
    buf.Read(&shear);
    fElasticity = EElasticity(elasticity);
    fShear = EShear(shear);
    fSigma0.Read(buf, context);
    buf.Read(&fK0User);
    buf.Read(&fTol);
    buf.Read(&fMaxIt);
    int mapping;
    buf.Read(&mapping);
    fStrengthMapping = EStrengthMapping(mapping);
    buf.Read(&fReductionFactor);
    fN.Read(buf, context);
    UpdateDerived();
}

void TPZModifiedCamClay::Print(std::ostream &out) const {
    static const char *elname[] = {"linear", "pressure_dependent"};
    static const char *shname[] = {"constant_G", "constant_nu", "hypo_nu"};
    out << "\n" << Name();
    fYC.Print(out);
    out << "\n N = " << fNIso << "  p0 = " << fP0
        << "\n elasticity = " << elname[fElasticity] << "  shear = " << shname[fShear]
        << "  G = " << fG << "  nu = " << fNu
        << "\n sigma0 = " << fSigma0
        << "\n p_ini = " << fPIni << "  K0 = " << fK0 << "  G0 = " << fG0
        << "\n tol = " << fTol << "  maxit = " << fMaxIt
        << "\n fN = ";
    fN.Print(out);
}

void TPZModifiedCamClay::ApplyStrain(const TPZTensor<REAL> &epsTotal) {
    TPZTensor<REAL> sigma;
    ApplyStrainComputeSigma(epsTotal, sigma, nullptr);
}

REAL TPZModifiedCamClay::MFromFriction(REAL phi, EStrengthMapping mapping) {
    const REAL s = std::sin(phi);
    if (mapping == EPlaneStrain) return std::sqrt(3.) * s;
    return 6. * s / (3. - s);
}

void TPZModifiedCamClay::ApplyLocalProperties() {
    const TPZVec<REAL> &mp = fN.fmatprop;
    if (mp.size() < 3 || !(mp[2] > 0.)) return;
    const REAL c = mp[0], phi = mp[1];
    const REAL tphi = std::tan(phi);
    if (!(tphi > 1.e-6)) {
        PZError << "TPZModifiedCamClay::ApplyLocalProperties: φ deve ser positivo (M = M(φ) > 0)\n";
        DebugStop();
    }
    const REAL phir = std::atan(tphi / fReductionFactor);
    const REAL M = MFromFriction(phir, fStrengthMapping);
    const REAL pt = std::max(c, REAL(0.)) / tphi;
    fYC.SetUp(M, fYC.Lambda(), fYC.Kappa(), fYC.V0(), mp[2], pt, fYC.Beta());
    if (mp.size() >= 9) {
        for (int i = 0; i < 6; i++) fSigma0[i] = mp[3 + i];
        UpdateDerived();
    }
}

void TPZModifiedCamClay::ApplyStrainComputeSigma(const TPZTensor<REAL> &epsTotal, TPZTensor<REAL> &sigma,
                                                 TPZFMatrix<REAL> *tangent) {
    ApplyLocalProperties();
    TResult res;
    if (fShear == EHypoNu) {
        bool zero = true;
        for (int i = 0; i < 6; i++) {
            if (sigma[i] != 0.) zero = false;
        }
        const TPZTensor<REAL> sig_n = zero ? fSigma0 : sigma;
        ReturnMapping(epsTotal, fN.m_eps_p, fN.m_hardening, res, &sig_n, &fN.m_eps_t);
    } else {
        ReturnMapping(epsTotal, fN.m_eps_p, fN.m_hardening, res);
    }
    sigma = res.stress;
    if (tangent) {
        tangent->Redim(6, 6);
        for (int i = 0; i < 6; i++) {
            for (int j = 0; j < 6; j++) (*tangent)(i, j) = res.Dep(i, j);
        }
    }
    fN.m_eps_t = epsTotal;
    fN.m_eps_p = res.plastic_strain;
    fN.m_hardening = res.alpha;
    fN.m_m_type = res.plastic ? 1 : 0;
    fLastIterations = res.iterations;
}

void TPZModifiedCamClay::ApplyStrainComputeDep(const TPZTensor<REAL> &epsTotal, TPZTensor<REAL> &sigma,
                                               TPZFMatrix<REAL> &Dep) {
    ApplyStrainComputeSigma(epsTotal, sigma, &Dep);
}

void TPZModifiedCamClay::ApplyLoad(const TPZTensor<REAL> &sigma, TPZTensor<REAL> &epsTotal) {
    const TPZPlasticState<REAL> state_n = fN;
    // no modo hipoelástico σ_n não faz parte do estado: usa-se σ0 (exato apenas a partir do estado inicial)
    const TPZTensor<REAL> sig_n(fSigma0);
    const REAL alpha_n = state_n.m_hardening;
    REAL scale = 0.;
    for (int i = 0; i < 6; i++) scale = std::max(scale, std::fabs(sigma[i]));
    scale = std::max(scale, std::fabs(fPIni));
    TPZTensor<REAL> eps(state_n.m_eps_t);
    TResult res;
    bool conv = false;
    for (int it = 0; it < 100; it++) {
        if (fShear == EHypoNu) {
            ReturnMapping(eps, state_n.m_eps_p, alpha_n, res, &sig_n, &state_n.m_eps_t);
        } else {
            ReturnMapping(eps, state_n.m_eps_p, alpha_n, res);
        }
        REAL A[36], r[6], rmax = 0.;
        for (int i = 0; i < 6; i++) {
            r[i] = sigma[i] - res.stress[i];
            rmax = std::max(rmax, std::fabs(r[i]));
            for (int j = 0; j < 6; j++) A[i * 6 + j] = res.Dep(i, j);
        }
        if (rmax <= 1.e-12 * scale) {
            conv = true;
            break;
        }
        if (!SolveDense(6, A, r)) break;
        for (int i = 0; i < 6; i++) eps[i] += r[i];
    }
    if (!conv) {
        PZError << "TPZModifiedCamClay::ApplyLoad: Newton não convergiu (tensão inatingível?)\n";
    }
    epsTotal = eps;
    fN.m_eps_t = eps;
    fN.m_eps_p = res.plastic_strain;
    fN.m_hardening = res.alpha;
    fN.m_m_type = res.plastic ? 1 : 0;
    fLastIterations = res.iterations;
}

void TPZModifiedCamClay::Phi(const TPZTensor<REAL> &epsElastic, TPZVec<REAL> &phi) const {
    REAL ee[6], de[6];
    ToArray(epsElastic, ee);
    const REAL xv = ee[_XX_] + ee[_YY_] + ee[_ZZ_];
    REAL p, K, dK, G, dG;
    Pressure(xv, p, K, dK);
    if (fShear == EHypoNu) {
        G = fGFac * ((fElasticity == ELinear) ? fK0 : K);
    } else {
        ShearModulus(xv, G, dG);
    }
    IDev(ee, de);
    REAL s[6];
    for (int i = 0; i < 6; i++) s[i] = 2. * G * (de[i] + fE0[i]);
    const REAL q = std::sqrt(1.5) * TNorm(s);
    REAL a, H;
    fYC.Hardening(fN.m_hardening, a, H);
    phi.Resize(1);
    phi[0] = fYC.YieldValue(p, q, a) / (a * a);
}

TPZElasticResponse TPZModifiedCamClay::GetElasticResponse() const {
    TPZElasticResponse ER;
    const REAL K = fK0, G = fG0;
    ER.SetEngineeringData(9. * K * G / (3. * K + G), (3. * K - 2. * G) / (2. * (3. * K + G)));
    return ER;
}
