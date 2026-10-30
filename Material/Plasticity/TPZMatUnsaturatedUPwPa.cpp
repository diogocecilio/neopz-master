/**
 * @file TPZMatUnsaturatedUPwPa.cpp
 * @brief Implementation of the u - p_w - p_a material for unsaturated soils (see TPZMatUnsaturatedUPwPa.h).
 */

#include "TPZMatUnsaturatedUPwPa.h"
#include "TPZPlasticStepVoigt.h"
#include "TPZYCMohrCoulombPV2.h"
#include "TPZPlasticStepPV.h"
#include "TPZYCMohrCoulombPV.h"
#include "TPZElasticResponse.h"
#include "TPZHash.h"
#include "TPZStream.h"
#include "pzerror.h"
#include <cmath>
#include <map>

// ---------------------------------------------------------------- memory
int TPZUnsaturatedMem::ClassId() const { return Hash("TPZUnsaturatedMem") ^ TPZElastoPlasticMem::ClassId() << 1; }

void TPZUnsaturatedMem::Write(TPZStream &buf, int withclassid) const {
    TPZElastoPlasticMem::Write(buf, withclassid);
    buf.Write(&fPw);
    buf.Write(&fPa);
    buf.Write(&fS);
    buf.Write(fPwV, 3);
    buf.Write(fPaV, 3);
}

void TPZUnsaturatedMem::Read(TPZStream &buf, void *context) {
    TPZElastoPlasticMem::Read(buf, context);
    buf.Read(&fPw);
    buf.Read(&fPa);
    buf.Read(&fS);
    buf.Read(fPwV, 3);
    buf.Read(fPaV, 3);
}

void TPZUnsaturatedMem::Print(std::ostream &out) const {
    TPZElastoPlasticMem::Print(out);
    out << " p_w = " << fPw << " p_a = " << fPa << " S_w = " << fS << " vertices p_w = {" << fPwV[0] << ", "
        << fPwV[1] << ", " << fPwV[2] << "} p_a = {" << fPaV[0] << ", " << fPaV[1] << ", " << fPaV[2] << "}\n";
}

// ---------------------------------------------------------------- material
template <class T, class TMEM>
TPZMatUnsaturatedUPwPa<T, TMEM>::TPZMatUnsaturatedUPwPa() : TBase(), fPlasticity() {}

template <class T, class TMEM>
TPZMatUnsaturatedUPwPa<T, TMEM>::TPZMatUnsaturatedUPwPa(int id) : TBase(id), fPlasticity() {}

template <class T, class TMEM>
void TPZMatUnsaturatedUPwPa<T, TMEM>::SetPlasticModel(const T &model) {
    fPlasticity = model;
    TMEM mem;
    mem.m_elastoplastic_state = fPlasticity.GetState();
    mem.m_sigma.Zero();
    this->SetDefaultMem(mem);
}

template <class T, class TMEM>
int TPZMatUnsaturatedUPwPa<T, TMEM>::IntegrationRuleOrder(const TPZVec<int> &elPMaxOrder) const {
    if (fIntegrationOrder > 0) return fIntegrationOrder;
    return TPZMatCombinedSpacesT<STATE>::IntegrationRuleOrder(elPMaxOrder);
}

template <class T, class TMEM>
void TPZMatUnsaturatedUPwPa<T, TMEM>::FillDataRequirements(TPZVec<TPZMaterialDataT<STATE>> &datavec) const {
    for (int i = 0; i < datavec.size(); ++i) {
        datavec[i].SetAllRequirements(false);
        datavec[i].fNeedsSol = true;
    }
}

template <class T, class TMEM>
void TPZMatUnsaturatedUPwPa<T, TMEM>::FillBoundaryConditionDataRequirements(
    int type, TPZVec<TPZMaterialDataT<STATE>> &datavec) const {
    for (int i = 0; i < datavec.size(); ++i) {
        datavec[i].SetAllRequirements(false);
        datavec[i].fNeedsSol = true;
    }
}

template <class T, class TMEM>
void TPZMatUnsaturatedUPwPa<T, TMEM>::GlobalGradients(const TPZMaterialDataT<STATE> &data, TPZFMatrix<REAL> &grad) {
    const int nshape = data.dphix.Cols();
    const int dimel = data.dphix.Rows();
    grad.Redim(2, nshape);
    for (int a = 0; a < nshape; ++a)
        for (int i = 0; i < 2; ++i) {
            REAL val = 0.;
            for (int k = 0; k < dimel; ++k) val += data.axes.GetVal(k, i) * data.dphix.GetVal(k, a);
            grad(i, a) = val;
        }
}

template <class T, class TMEM>
void TPZMatUnsaturatedUPwPa<T, TMEM>::Bishop(REAL pw, REAL pa, REAL S, REAL Gamma, REAL &pB, REAL &chiw,
                                             REAL &chia) const {
    if (fBishop == EChiConstant) {
        // phi^b = atan(chi tan phi') on the suction; saturated (pw > pa): chi = 1
        const REAL chi = pw > pa ? 1. : fChi;
        pB = chi * pw + (1. - chi) * pa;
        chiw = chi;
        chia = 1. - chi;
    } else if (fBishop == EChiEffective) {
        REAL Se, dSe;
        fHydro.Effective(S, Se, dSe);
        if (Se <= 0. || Se >= 1.) dSe = 0.;
        pB = Se * pw + (1. - Se) * pa;
        chiw = Se + dSe * Gamma * (pw - pa);
        chia = (1. - Se) - dSe * Gamma * (pw - pa);
    } else {
        pB = S * pw + (1. - S) * pa;
        chiw = S + Gamma * (pw - pa); // chi_w = S_w - Gamma p_c
        chia = (1. - S) - Gamma * (pw - pa);
    }
}

template <class T, class TMEM>
void TPZMatUnsaturatedUPwPa<T, TMEM>::UnitWeight(REAL S, REAL Gamma, REAL &gamma, REAL &dgamma) const {
    if (fMixture) {
        gamma = (1. - fPorosity) * fGammaS + fPorosity * S * fGammaW;
        dgamma = fPorosity * fGammaW * Gamma;
    } else {
        gamma = fGamma;
        dgamma = 0.;
    }
}

namespace {
/** @brief Coefficients c of the P1 basis with sum c_i phi_i = p and sum c_i grad phi_i = grad p (3 x 3 system) */
void VertexValues(const TPZFMatrix<REAL> &phi, const TPZFMatrix<REAL> &grad, REAL p, const REAL gp[2], REAL c[3]) {
    REAL A[3][3], b[3] = {p, gp[0], gp[1]};
    for (int j = 0; j < 3; ++j) {
        A[0][j] = phi.GetVal(j, 0);
        A[1][j] = grad.GetVal(0, j);
        A[2][j] = grad.GetVal(1, j);
    }
    const REAL det = A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1]) - A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0]) +
                     A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]);
    for (int k = 0; k < 3; ++k) {
        REAL M[3][3];
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) M[i][j] = (j == k) ? b[i] : A[i][j];
        c[k] = (M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1]) - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0]) +
                M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0])) /
               det;
    }
}
} // namespace

template <class T, class TMEM>
void TPZMatUnsaturatedUPwPa<T, TMEM>::Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight,
                                                 TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) {
    ContributeInternal(datavec, weight, &ek, ef);
}

template <class T, class TMEM>
void TPZMatUnsaturatedUPwPa<T, TMEM>::Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight,
                                                 TPZFMatrix<STATE> &ef) {
    ContributeInternal(datavec, weight, nullptr, ef);
}

template <class T, class TMEM>
void TPZMatUnsaturatedUPwPa<T, TMEM>::ContributeInternal(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight,
                                                         TPZFMatrix<STATE> *ek, TPZFMatrix<STATE> &ef) {
    const TPZMaterialDataT<STATE> &dU = datavec[0], &dW = datavec[1], &dA = datavec[2];
    const int nu = dU.phi.Rows(), np = dW.phi.Rows();
    if (dA.phi.Rows() != np) DebugStop(); // p_w and p_a in the same space
    const int nequ = 2 * nu, ow = nequ, oa = nequ + np;
    const REAL dv = weight;
    if (fGravityOnly) {
        const TPZUnsaturatedHydraulics::TState h0 = fHydro.Evaluate(dW.sol[0][0], dA.sol[0][0]);
        REAL gam, dgam;
        UnitWeight(h0.S, h0.Gamma, gam, dgam);
        for (int a = 0; a < nu; ++a)
            for (int d = 0; d < 2; ++d) ef(2 * a + d, 0) += dv * fGravityFactor * dU.phi.GetVal(a, 0) * gam * fG[d];
        return;
    }

    TPZFNMatrix<40, REAL> gU, gP;
    GlobalGradients(dU, gU);
    GlobalGradients(dW, gP);
    // strain (engineering shear) and B operator columns: mB = divergence of each displacement function
    TPZTensor<REAL> eps;
    {
        const TPZFMatrix<STATE> &ds = dU.dsol[0];
        REAL du[2][2];
        for (int c = 0; c < 2; ++c)
            for (int i = 0; i < 2; ++i) {
                REAL v = 0.;
                for (int k = 0; k < ds.Rows(); ++k) v += dU.axes.GetVal(k, i) * ds.GetVal(k, c);
                du[c][i] = v;
            }
        eps.XX() = du[0][0];
        eps.YY() = du[1][1];
        eps.XY() = du[0][1] + du[1][0];
    }
    auto gradOf = [](const TPZMaterialDataT<STATE> &d, REAL g[2]) {
        const TPZFMatrix<STATE> &ds = d.dsol[0];
        for (int i = 0; i < 2; ++i) {
            REAL v = 0.;
            for (int k = 0; k < ds.Rows(); ++k) v += d.axes.GetVal(k, i) * ds.GetVal(k, 0);
            g[i] = v;
        }
    };
    const REAL pw = dW.sol[0][0], pa = dA.sol[0][0];
    REAL gpw[2], gpa[2];
    gradOf(dW, gpw);
    gradOf(dA, gpa);

    // memory: converged state
    const int64_t ip = dU.intGlobPtIndex;
    if (ip < 0) DebugStop();
    TMEM &mem = this->MemItem(ip);
    T model(fPlasticity);
    model.SetState(mem.m_elastoplastic_state);
    TPZTensor<REAL> sigma(mem.m_sigma);
    TPZFNMatrix<36, REAL> Dep(6, 6, 0.);
    model.ApplyStrainComputeSigma(eps, sigma, ek ? &Dep : nullptr);
    const REAL evn = mem.m_elastoplastic_state.m_eps_t.I1();
    const REAL pwn = mem.fPw, pan = mem.fPa;

    // hydraulic state at the point
    const TPZUnsaturatedHydraulics::TState hs = fHydro.Evaluate(pw, pa);
    const REAL S = hs.S, Gm = hs.Gamma, Sa = 1. - S;
    REAL pB, chiw, chia;
    Bishop(pw, pa, S, Gm, pB, chiw, chia);
    REAL gam, dgam;
    UnitWeight(S, Gm, gam, dgam);
    // factor of p_B in the momentum balance: lambda_g, or 1 on the suction (pB < 0) if the suction is not scaled
    const REAL lg = fGravityFactor, lp = (fScaleSuction || pB > 0.) ? fGravityFactor : 1.;

    REAL vw[3] = {0., 0., 0.}, va[3] = {0., 0., 0.};
    const bool lump = fLumped && np == 3;
    if (lump) {
        VertexValues(dW.phi, gP, pw, gpw, vw);
        VertexValues(dA.phi, gP, pa, gpa, va);
    }
    if (this->fUpdateMem) {
        mem.m_sigma = sigma;
        mem.m_elastoplastic_state = model.GetState();
        mem.fPw = pw;
        mem.fPa = pa;
        mem.fS = S;
        if (lump)
            for (int i = 0; i < 3; ++i) {
                mem.fPwV[i] = vw[i];
                mem.fPaV[i] = va[i];
            }
    }

    // ---------------- momentum balance
    TPZFNMatrix<40, REAL> mB(nequ, 1, 0.);
    for (int a = 0; a < nu; ++a) {
        mB(2 * a, 0) = gU.GetVal(0, a);
        mB(2 * a + 1, 0) = gU.GetVal(1, a);
    }
    for (int a = 0; a < nu; ++a) {
        const REAL ga0 = gU.GetVal(0, a), ga1 = gU.GetVal(1, a), Na = dU.phi.GetVal(a, 0);
        // B^T sigma' (Voigt XX, XY, XZ, YY, YZ, ZZ; engineering shear)
        const REAL fx = ga0 * sigma.XX() + ga1 * sigma.XY();
        const REAL fy = ga1 * sigma.YY() + ga0 * sigma.XY();
        ef(2 * a, 0) -= dv * (fx - lp * fAlpha * ga0 * pB - lg * Na * gam * fG[0]);
        ef(2 * a + 1, 0) -= dv * (fy - lp * fAlpha * ga1 * pB - lg * Na * gam * fG[1]);
    }
    if (ek) {
        // B (6 x nequ)
        TPZFNMatrix<240, REAL> B(6, nequ, 0.), DB(6, nequ, 0.);
        for (int a = 0; a < nu; ++a) {
            B(_XX_, 2 * a) = gU.GetVal(0, a);
            B(_XY_, 2 * a) = gU.GetVal(1, a);
            B(_XY_, 2 * a + 1) = gU.GetVal(0, a);
            B(_YY_, 2 * a + 1) = gU.GetVal(1, a);
        }
        for (int l = 0; l < 6; ++l)
            for (int J = 0; J < nequ; ++J) {
                REAL v = 0.;
                for (int m = 0; m < 6; ++m) v += Dep.GetVal(l, m) * B.GetVal(m, J);
                DB(l, J) = v;
            }
        for (int I = 0; I < nequ; ++I) {
            const REAL NI = dU.phi.GetVal(I / 2, 0), gI = fG[I % 2];
            for (int J = 0; J < nequ; ++J) {
                REAL v = 0.;
                for (int l = 0; l < 6; ++l) v += B.GetVal(l, I) * DB.GetVal(l, J);
                (*ek)(I, J) += dv * v;
            }
            for (int j = 0; j < np; ++j) {
                const REAL Nj = dW.phi.GetVal(j, 0);
                (*ek)(I, ow + j) += dv * (-lp * fAlpha * mB(I, 0) * chiw - lg * NI * gI * dgam) * Nj;
                (*ek)(I, oa + j) += dv * (-lp * fAlpha * mB(I, 0) * chia + lg * NI * gI * dgam) * Nj;
            }
        }
    }
    if (fFrozen) return;

    // ---------------- mass balances of the water and of the air
    const REAL n = fPorosity, dt = fTimeStep;
    // uncoupled seepage (SetHydroMechanicalCoupling(false)): the deformation does not enter the mass balances
    const REAL alc = fCoupled ? fAlpha : 0.;
    const REAL eta = (fCoupled && fKs > 0.) ? (fAlpha - n) / fKs : 0.;
    const REAL dEv = eps.I1() - evn;
    const REAL dpw = pw - pwn, dpa = pa - pan;
    REAL pBn = 0., chin_w, chin_a;
    if (eta != 0.) {
        const TPZUnsaturatedHydraulics::TState hn = fHydro.Evaluate(pwn, pan);
        Bishop(pwn, pan, hn.S, hn.Gamma, pBn, chin_w, chin_a);
    }
    const REAL dpB = (eta != 0.) ? pB - pBn : 0.;
    const REAL X = alc * dEv + eta * dpB;
    // conductivity tensor K = [k_xx k_xy; k_xy k_yy]; mobilities lambda = K k_r / gamma_w (times mu_w/mu_a)
    const REAL kxx = fKx, kyy = fKy, kxy = fKxy;
    const REAL lw0 = 1. / fGammaW, la0 = lw0 * fMuRatio;
    const REAL lw = lw0 * hs.krw, la = la0 * hs.kra;
    const REAL dlw = lw0 * hs.dkrw * Gm, dla = la0 * hs.dkra * Gm; // d/dp_w (d/dp_a = -d/dp_w)
    const REAL hw[2] = {gpw[0] - fGammaW * fG[0], gpw[1] - fGammaW * fG[1]};
    const REAL ha[2] = {gpa[0] - fGammaA * fG[0], gpa[1] - fGammaA * fG[1]};
    const REAL Pabs = pa + fPatm, cag = 1. / Pabs;
    const REAL Khw[2] = {kxx * hw[0] + kxy * hw[1], kxy * hw[0] + kyy * hw[1]};
    const REAL Kha[2] = {kxx * ha[0] + kxy * ha[1], kxy * ha[0] + kyy * ha[1]};
    // K (h_a + grad p_a): derivative of the advective term of the air with respect to grad p_a
    const REAL Kag[2] = {kxx * (ha[0] + gpa[0]) + kxy * (ha[1] + gpa[1]), kxy * (ha[0] + gpa[0]) + kyy * (ha[1] + gpa[1])};
    const REAL adv = Kha[0] * gpa[0] + Kha[1] * gpa[1]; // -q_a . grad p_a / lambda_a

    // saturation of the last converged step at the point (storage integrated at the points, no lumping)
    const REAL Sn = lump ? 0. : fHydro.Evaluate(pwn, pan).S;
    for (int j = 0; j < np; ++j) {
        const REAL Nj = dW.phi.GetVal(j, 0), gj0 = gP.GetVal(0, j), gj1 = gP.GetVal(1, j);
        REAL Rw = Nj * (S * X + n * S / fKw * dpw) + dt * lw * (gj0 * Khw[0] + gj1 * Khw[1]);
        REAL Ra = Nj * (Sa * X - dt * cag * la * adv) + dt * la * (gj0 * Kha[0] + gj1 * Kha[1]);
        if (!lump) {
            Rw += Nj * n * (S - Sn);
            // air: n Delta S_a + n S_a^n Delta p_a / p_a^abs (mass of air n S_a rho_a, rho_a ~ p_a^abs, divided by
            // rho_a^{n+1}: exactly conservative for uniform p_a)
            Ra += Nj * (-n * (S - Sn) + n * (1. - Sn) * cag * dpa);
        } else {
            // storage at the vertex j with the vertex pressures: (int N_j) n (S(p_j) - S(p_j^n)) ...
            const TPZUnsaturatedHydraulics::TState hv = fHydro.Evaluate(vw[j], va[j]);
            REAL Svn, Gvn;
            fHydro.Retention(mem.fPaV[j] - mem.fPwV[j], Svn, Gvn);
            const REAL Pv = va[j] + fPatm, dpav = va[j] - mem.fPaV[j];
            Rw += Nj * n * (hv.S - Svn);
            Ra += Nj * (-n * (hv.S - Svn) + n * (1. - Svn) / Pv * dpav);
            if (ek) {
                const REAL Gv = hv.Gamma, M = dv * Nj * n;
                (*ek)(ow + j, ow + j) += M * Gv;
                (*ek)(ow + j, oa + j) -= M * Gv;
                (*ek)(oa + j, ow + j) += M * (-Gv);
                (*ek)(oa + j, oa + j) += M * (Gv + (1. - Svn) / Pv - (1. - Svn) * dpav / (Pv * Pv));
            }
        }
        ef(ow + j, 0) -= dv * Rw;
        ef(oa + j, 0) -= dv * Ra;
        if (!ek) continue;
        for (int I = 0; I < nequ; ++I) {
            (*ek)(ow + j, I) += dv * Nj * S * alc * mB(I, 0);
            (*ek)(oa + j, I) += dv * Nj * Sa * alc * mB(I, 0);
        }
        const REAL hwj = gj0 * Khw[0] + gj1 * Khw[1], haj = gj0 * Kha[0] + gj1 * Kha[1];
        for (int k = 0; k < np; ++k) {
            const REAL Nk = dW.phi.GetVal(k, 0), gk0 = gP.GetVal(0, k), gk1 = gP.GetVal(1, k);
            const REAL gjk = gj0 * (kxx * gk0 + kxy * gk1) + gj1 * (kxy * gk0 + kyy * gk1);
            REAL kww = Nj * (Gm * X + S * eta * chiw + n / fKw * (Gm * dpw + S)) * Nk + dt * (lw * gjk + dlw * hwj * Nk);
            REAL kwa = Nj * (-Gm * X + S * eta * chia - n / fKw * Gm * dpw) * Nk - dt * dlw * hwj * Nk;
            REAL kaw = Nj * (-Gm * X + Sa * eta * chiw - dt * cag * dla * adv) * Nk + dt * dla * haj * Nk;
            REAL kaa = Nj * (Gm * X + Sa * eta * chia + dt * cag * cag * la * adv + dt * cag * dla * adv) * Nk -
                       dt * Nj * cag * la * (Kag[0] * gk0 + Kag[1] * gk1) +
                       dt * (la * gjk - dla * haj * Nk);
            if (!lump) {
                kww += Nj * n * Gm * Nk;
                kwa -= Nj * n * Gm * Nk;
                kaw += Nj * (-n * Gm) * Nk;
                kaa += Nj * (n * Gm + n * (1. - Sn) * (cag - dpa * cag * cag)) * Nk;
            }
            (*ek)(ow + j, ow + k) += dv * kww;
            (*ek)(ow + j, oa + k) += dv * kwa;
            (*ek)(oa + j, ow + k) += dv * kaw;
            (*ek)(oa + j, oa + k) += dv * kaa;
        }
    }
}

template <class T, class TMEM>
void TPZMatUnsaturatedUPwPa<T, TMEM>::ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight,
                                                   TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef,
                                                   TPZBndCondT<STATE> &bc) {
    ContributeBCInternal(datavec, weight, &ek, ef, bc);
}

template <class T, class TMEM>
void TPZMatUnsaturatedUPwPa<T, TMEM>::ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight,
                                                   TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc) {
    ContributeBCInternal(datavec, weight, nullptr, ef, bc);
}

template <class T, class TMEM>
void TPZMatUnsaturatedUPwPa<T, TMEM>::ContributeBCInternal(const TPZVec<TPZMaterialDataT<STATE>> &datavec,
                                                           REAL weight, TPZFMatrix<STATE> *ek, TPZFMatrix<STATE> &ef,
                                                           TPZBndCondT<STATE> &bc) {
    const TPZMaterialDataT<STATE> &dU = datavec[0], &dW = datavec[1];
    const int nu = dU.phi.Rows(), np = dW.phi.Rows();
    const int nequ = 2 * nu;
    TPZManVector<STATE, 3> v2(3, 0.);
    for (int i = 0; i < bc.Val2().size() && i < 3; ++i) v2[i] = bc.Val2()[i];
    if (bc.HasForcingFunctionBC()) {
        TPZFNMatrix<9, STATE> v1(3, 3, 0.);
        bc.ForcingFunctionBC()(dU.x, v2, v1);
    }
    switch (bc.Type()) {
    case ENeumannU:
        for (int a = 0; a < nu; ++a)
            for (int d = 0; d < 2; ++d) ef(2 * a + d, 0) += weight * v2[d] * dU.phi.GetVal(a, 0);
        break;
    case EFluxW:
    case EFluxA: {
        if (fFrozen) break;
        const int off = nequ + (bc.Type() == EFluxW ? 0 : np);
        for (int j = 0; j < np; ++j) ef(off + j, 0) += weight * fTimeStep * v2[0] * dW.phi.GetVal(j, 0);
        break;
    }
    case EDirichletU:
    case EDirichletPw:
    case EDirichletPa:
        break; // imposed by the analysis
    default:
        std::cout << __PRETTY_FUNCTION__ << " unknown boundary condition type " << bc.Type() << std::endl;
        DebugStop();
    }
}

template <class T, class TMEM>
void TPZMatUnsaturatedUPwPa<T, TMEM>::ResetDisplacementReference() {
    auto mem = this->GetMemory();
    for (int64_t i = 0; i < mem->NElements(); ++i) {
        TMEM &m = (*mem)[i];
        m.m_elastoplastic_state.m_eps_p -= m.m_elastoplastic_state.m_eps_t;
        m.m_elastoplastic_state.m_eps_t.Zero();
    }
}

template <class T, class TMEM>
int64_t TPZMatUnsaturatedUPwPa<T, TMEM>::NPlasticPoints() const {
    auto mem = const_cast<TPZMatUnsaturatedUPwPa *>(this)->GetMemory();
    int64_t count = 0;
    for (int64_t i = 0; i < mem->NElements(); ++i)
        if ((*mem)[i].m_elastoplastic_state.m_m_type != 0) count++;
    return count;
}

template <class T, class TMEM>
int TPZMatUnsaturatedUPwPa<T, TMEM>::VariableIndex(const std::string &name) const {
    static const std::map<std::string, int> names = {
        {"Displacement", EDisplacement},   {"WaterPressure", EWaterPressure}, {"AirPressure", EAirPressure},
        {"Suction", ESuction},             {"Saturation", ESaturation},       {"RelativePermeabilityWater", EKrw},
        {"BishopPressure", EBishopPressure}, {"DisplacementX", EDisplacementX}, {"DisplacementY", EDisplacementY}};
    auto it = names.find(name);
    if (it != names.end()) return it->second;
    return TBase::VariableIndex(name);
}

template <class T, class TMEM>
int TPZMatUnsaturatedUPwPa<T, TMEM>::NSolutionVariables(int var) const {
    if (var == EDisplacement) return 3;
    if (var >= EWaterPressure && var <= EDisplacementY) return 1;
    return TBase::NSolutionVariables(var);
}

template <class T, class TMEM>
void TPZMatUnsaturatedUPwPa<T, TMEM>::Solution(const TPZVec<TPZMaterialDataT<STATE>> &datavec, int var,
                                               TPZVec<STATE> &sol) {
    const TPZVec<STATE> &u = datavec[0].sol[0];
    const REAL pw = datavec[1].sol[0][0], pa = datavec[2].sol[0][0];
    sol.Resize(NSolutionVariables(var));
    sol.Fill(0.);
    const TPZUnsaturatedHydraulics::TState hs = fHydro.Evaluate(pw, pa);
    switch (var) {
    case EDisplacement:
        sol[0] = u[0];
        sol[1] = u[1];
        break;
    case EDisplacementX: sol[0] = u[0]; break;
    case EDisplacementY: sol[0] = u[1]; break;
    case EWaterPressure: sol[0] = pw; break;
    case EAirPressure: sol[0] = pa; break;
    case ESuction: sol[0] = pa - pw; break;
    case ESaturation: sol[0] = hs.S; break;
    case EKrw: sol[0] = hs.krw; break;
    case EBishopPressure: {
        REAL pB, cw, ca;
        Bishop(pw, pa, hs.S, hs.Gamma, pB, cw, ca);
        sol[0] = pB;
        break;
    }
    default:
        break;
    }
}

template <class T, class TMEM>
void TPZMatUnsaturatedUPwPa<T, TMEM>::Print(std::ostream &out) const {
    out << Name() << " id = " << this->Id() << "\n  n = " << fPorosity << " alpha = " << fAlpha << " Ks = " << fKs
        << " Kw = " << fKw << "\n  k_x = " << fKx << " k_y = " << fKy << " k_xy = " << fKxy << " gamma_w = " << fGammaW << " mu_w/mu_a = " << fMuRatio
        << " gamma_a = " << fGammaA << " p_atm = " << fPatm << "\n  unit weight "
        << (fMixture ? "mixture, gamma_s = " : "constant ") << (fMixture ? fGammaS : fGamma) << "\n  chi "
        << (fBishop == EChiSaturation ? "= S_w" : (fBishop == EChiEffective ? "= S_e" : "constant "))
        << (fBishop == EChiConstant ? fChi : 0.)
        << "\n  storage " << (fLumped ? "lumped" : "consistent") << " dt = " << fTimeStep
        << " gravity factor = " << fGravityFactor << (fFrozen ? " (frozen pressures)" : "") << "\n";
    fHydro.Print(out);
    fPlasticity.Print(out);
}

template <class T, class TMEM>
int TPZMatUnsaturatedUPwPa<T, TMEM>::ClassId() const {
    return Hash("TPZMatUnsaturatedUPwPa") ^ TBase::ClassId() << 1 ^ fPlasticity.ClassId() << 2;
}

template <class T, class TMEM>
void TPZMatUnsaturatedUPwPa<T, TMEM>::Write(TPZStream &buf, int withclassid) const {
    TBase::Write(buf, withclassid);
    REAL vals[] = {fPorosity, fAlpha, fKs, fKw, fKx, fKy, fKxy, fGammaW, fMuRatio, fGammaA, fPatm, fGamma, fGammaS,
                   fG[0], fG[1], fG[2], fChi, fTimeStep, fGravityFactor};
    buf.Write(vals, sizeof(vals) / sizeof(REAL));
    int ints[] = {fMixture, int(fBishop), fLumped, fIntegrationOrder, fFrozen, fCoupled};
    buf.Write(ints, 6);
    fPlasticity.Write(buf, withclassid);
}

template <class T, class TMEM>
void TPZMatUnsaturatedUPwPa<T, TMEM>::Read(TPZStream &buf, void *context) {
    TBase::Read(buf, context);
    REAL vals[19];
    buf.Read(vals, 19);
    REAL *p[] = {&fPorosity, &fAlpha, &fKs, &fKw, &fKx, &fKy, &fKxy, &fGammaW, &fMuRatio, &fGammaA, &fPatm, &fGamma, &fGammaS,
                 &fG[0], &fG[1], &fG[2], &fChi, &fTimeStep, &fGravityFactor};
    for (int i = 0; i < 19; ++i) *p[i] = vals[i];
    int ints[6];
    buf.Read(ints, 6);
    fMixture = ints[0];
    fBishop = EBishop(ints[1]);
    fLumped = ints[2];
    fIntegrationOrder = ints[3];
    fFrozen = ints[4];
    fCoupled = ints[5];
    fPlasticity.Read(buf, context);
}

template class TPZMatUnsaturatedUPwPa<TPZPlasticStepVoigt<TPZYCMohrCoulombPV2, TPZElasticResponse>, TPZUnsaturatedMem>;
// legacy iterative Mohr-Coulomb (non-associative: psi != phi)
template class TPZMatUnsaturatedUPwPa<TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse>, TPZUnsaturatedMem>;
