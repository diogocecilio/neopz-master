// TPZMatPoroElastoPlastic3DMem.cpp — ver TPZMatPoroElastoPlastic3DMem.h

#include "TPZMatPoroElastoPlastic3DMem.h"

#include <array>
#include <vector>

#include "TPZBndCondT.h"
#include "TPZStream.h"
#include "pzerror.h"
#include "TPZModifiedCamClay.h"
#include "TPZPlasticStepVoigt.h"
#include "TPZYCMohrCoulombPV2.h"
#include "TPZElasticResponse.h"

namespace {

/// Componentes não nulas da coluna (3a + i) de B (Voigt {xx, xy, xz, yy, yz, zz}, distorções de engenharia)
inline void BColumn(int i, const REAL dN[3], int rows[3], REAL vals[3]) {
    switch (i) {
        case 0:
            rows[0] = _XX_; vals[0] = dN[0];
            rows[1] = _XY_; vals[1] = dN[1];
            rows[2] = _XZ_; vals[2] = dN[2];
            break;
        case 1:
            rows[0] = _YY_; vals[0] = dN[1];
            rows[1] = _XY_; vals[1] = dN[0];
            rows[2] = _YZ_; vals[2] = dN[2];
            break;
        default:
            rows[0] = _ZZ_; vals[0] = dN[2];
            rows[1] = _XZ_; vals[1] = dN[0];
            rows[2] = _YZ_; vals[2] = dN[1];
            break;
    }
}

/// Derivadas nas coordenadas globais: out = axesᵀ in (in: dim x n, nas direções de axes)
void ToXYZ(const TPZFMatrix<REAL> &axes, const TPZFMatrix<REAL> &in, TPZFMatrix<REAL> &out) {
    axes.Multiply(in, out, 1);
}

} // namespace

template <class T, class TMEM>
TPZMatPoroElastoPlastic3DMem<T, TMEM>::TPZMatPoroElastoPlastic3DMem() : TBase() {}

template <class T, class TMEM>
TPZMatPoroElastoPlastic3DMem<T, TMEM>::TPZMatPoroElastoPlastic3DMem(int matid, int dim) : TBase(matid), fDim(dim) {
    if (dim != 2 && dim != 3) {
        PZError << Name() << ": dimensão " << dim << " inválida (2: deformação plana, 3: sólido)\n";
        DebugStop();
    }
}

template <class T, class TMEM>
int TPZMatPoroElastoPlastic3DMem<T, TMEM>::IntegrationRuleOrder(const TPZVec<int> &elPMaxOrder) const {
    if (fIntegrationOrder > 0) return fIntegrationOrder;
    return TPZMatCombinedSpacesT<STATE>::IntegrationRuleOrder(elPMaxOrder);
}

template <class T, class TMEM>
void TPZMatPoroElastoPlastic3DMem<T, TMEM>::FillDataRequirements(TPZVec<TPZMaterialDataT<STATE>> &datavec) const {
    for (auto &d : datavec) {
        d.SetAllRequirements(false);
        d.fNeedsSol = true;
    }
}

template <class T, class TMEM>
void TPZMatPoroElastoPlastic3DMem<T, TMEM>::FillBoundaryConditionDataRequirements(
    int type, TPZVec<TPZMaterialDataT<STATE>> &datavec) const {
    for (auto &d : datavec) {
        d.SetAllRequirements(false);
        d.fNeedsSol = true;
        d.fNeedsNormal = (type == ENormalPressure);
    }
}

// ---------------------------------------------------------------------------------------------------
template <class T, class TMEM>
void TPZMatPoroElastoPlastic3DMem<T, TMEM>::ContributeInternal(const TPZVec<TPZMaterialDataT<STATE>> &datavec,
                                                               REAL weight, TPZFMatrix<STATE> *ek,
                                                               TPZFMatrix<STATE> &ef) {
    if (datavec.size() != 2) {
        PZError << Name() << "::Contribute: datavec deve ter 2 espaços (u, p)\n";
        DebugStop();
    }
    const TPZMaterialDataT<STATE> &dU = datavec[0], &dP = datavec[1];
    const int dim = fDim;
    const int nU = dU.phi.Rows(), nP = dP.phi.Rows(), offp = dim * nU;

    if (fMode == EExternalForces) {
        for (int a = 0; a < nU; a++)
            for (int c = 0; c < dim; c++) ef(dim * a + c, 0) += dU.phi.GetVal(a, 0) * fBody[c] * weight;
        return;
    }

    TPZFNMatrix<180, REAL> dphiU, dphiP;
    ToXYZ(dU.axes, dU.dphix, dphiU);
    ToXYZ(dP.axes, dP.dphix, dphiP);
    TPZFNMatrix<9, REAL> D, gp;
    ToXYZ(dU.axes, dU.dsol[0], D);    // D(i,j) = ∂u_j/∂x_i (3 x dim)
    ToXYZ(dP.axes, dP.dsol[0], gp);   // ∇p (3 x 1)
    TPZTensor<REAL> eps;              // deformação plana: ε_zz = ε_xz = ε_yz = 0
    eps[_XX_] = D(0, 0);
    eps[_YY_] = D(1, 1);
    eps[_XY_] = D(1, 0) + D(0, 1);
    if (dim == 3) {
        eps[_ZZ_] = D(2, 2);
        eps[_XZ_] = D(2, 0) + D(0, 2);
        eps[_YZ_] = D(2, 1) + D(1, 2);
    }
    const REAL p = dP.sol[0][0];

    // ---- lei constitutiva a partir do estado convergido (memória)
    const int64_t igp = dU.intGlobPtIndex;
    if (igp < 0) {
        PZError << Name() << "::Contribute: ponto de integração sem memória (malha de u deve ter memória)\n";
        DebugStop();
    }
    TMEM &mem = this->MemItem(igp);
    T model(fPlasticModel);
    if (fModelUpdate) fModelUpdate(dU.x, model);
    model.SetState(mem.m_elastoplastic_state);
    TPZTensor<REAL> sigma(mem.m_sigma);  // σ'_n na entrada (leis hipoelásticas)
    TPZFNMatrix<36, REAL> Dep(6, 6, 0.);
    model.ApplyStrainComputeSigma(eps, sigma, ek ? &Dep : nullptr);
    const TPZTensor<REAL> &epsn = mem.m_elastoplastic_state.m_eps_t;
    const REAL evn = epsn[_XX_] + epsn[_YY_] + epsn[_ZZ_];
    const REAL ev = eps[_XX_] + eps[_YY_] + eps[_ZZ_];
    const REAL pn = mem.fPorePressure;
    const REAL mob[3] = {fK[0] / fMu, fK[1] / fMu, fK[2] / fMu};

    // ---- resíduo (ef = -R)
    for (int a = 0; a < nU; a++) {
        const REAL d0 = dphiU(0, a), d1 = dphiU(1, a), d2 = dphiU(2, a), N = dU.phi.GetVal(a, 0);
        const REAL f[3] = {sigma[_XX_] * d0 + sigma[_XY_] * d1 + sigma[_XZ_] * d2,
                           sigma[_XY_] * d0 + sigma[_YY_] * d1 + sigma[_YZ_] * d2,
                           sigma[_XZ_] * d0 + sigma[_YZ_] * d1 + sigma[_ZZ_] * d2};
        const REAL dN[3] = {d0, d1, d2};
        for (int c = 0; c < dim; c++) ef(dim * a + c, 0) += (-f[c] + fAlpha * p * dN[c] + N * fBody[c]) * weight;
    }
    for (int b = 0; b < nP; b++) {
        REAL r = dP.phi.GetVal(b, 0) * (fAlpha * (ev - evn) + fSe * (p - pn));
        if (fFlow)
            for (int c = 0; c < dim; c++) r += fTimeStep * mob[c] * dphiP(c, b) * (gp(c, 0) - fRhoF * fG[c]);
        ef(offp + b, 0) -= r * weight;
    }

    // ---- jacobiana
    if (ek) {
        const int nd = dim * nU;
        std::vector<std::array<REAL, 6>> DB(nd);
        std::vector<std::array<int, 3>> brow(nd);
        std::vector<std::array<REAL, 3>> bval(nd);
        for (int a = 0; a < nU; a++) {
            const REAL dN[3] = {dphiU(0, a), dphiU(1, a), dim == 3 ? dphiU(2, a) : 0.};
            for (int i = 0; i < dim; i++) {
                const int col = dim * a + i;
                BColumn(i, dN, brow[col].data(), bval[col].data());
                for (int rr = 0; rr < 6; rr++) {
                    REAL v = 0.;
                    for (int k = 0; k < 3; k++) v += Dep(rr, brow[col][k]) * bval[col][k];
                    DB[col][rr] = v;
                }
            }
        }
        for (int row = 0; row < nd; row++)
            for (int col = 0; col < nd; col++) {
                REAL v = 0.;
                for (int k = 0; k < 3; k++) v += bval[row][k] * DB[col][brow[row][k]];
                (*ek)(row, col) += v * weight;
            }
        for (int a = 0; a < nU; a++)
            for (int c = 0; c < dim; c++)
                for (int b = 0; b < nP; b++) {
                    const REAL q = fAlpha * dphiU(c, a) * dP.phi.GetVal(b, 0) * weight;
                    (*ek)(dim * a + c, offp + b) -= q;
                    (*ek)(offp + b, dim * a + c) += q;
                }
        for (int b = 0; b < nP; b++)
            for (int c = 0; c < nP; c++) {
                REAL v = fSe * dP.phi.GetVal(b, 0) * dP.phi.GetVal(c, 0);
                if (fFlow)
                    for (int k = 0; k < dim; k++) v += fTimeStep * mob[k] * dphiP(k, b) * dphiP(k, c);
                (*ek)(offp + b, offp + c) += v * weight;
            }
    }

    // ---- atualização da memória (passo aceito)
    if (this->GetUpdateMem()) {
        mem.m_elastoplastic_state = model.GetState();
        mem.m_sigma = sigma;
        mem.fPorePressure = p;
        mem.m_plastic_steps = model.IntegrationSteps();
        mem.m_u.Resize(dim);
        for (int c = 0; c < dim; c++) mem.m_u[c] = dU.sol[0][c];
    }
}

template <class T, class TMEM>
void TPZMatPoroElastoPlastic3DMem<T, TMEM>::Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight,
                                                       TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) {
    ContributeInternal(datavec, weight, &ek, ef);
}

template <class T, class TMEM>
void TPZMatPoroElastoPlastic3DMem<T, TMEM>::Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight,
                                                       TPZFMatrix<STATE> &ef) {
    ContributeInternal(datavec, weight, nullptr, ef);
}

// ---------------------------------------------------------------------------------------------------
template <class T, class TMEM>
void TPZMatPoroElastoPlastic3DMem<T, TMEM>::ContributeBCInternal(const TPZVec<TPZMaterialDataT<STATE>> &datavec,
                                                                 REAL weight, TPZFMatrix<STATE> *ek,
                                                                 TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc) {
    const TPZMaterialDataT<STATE> &dU = datavec[0], &dP = datavec[1];
    const int dim = fDim;
    const int nU = dU.phi.Rows(), nP = dP.phi.Rows(), offp = dim * nU;
    // valores da condição: Val2 ou a função do ponto (rhsVal com 4 componentes {v_x, v_y, v_z, p})
    TPZManVector<STATE, 4> v2(bc.Val2());
    if (bc.HasForcingFunctionBC()) {
        v2.Resize(4);
        v2.Fill(0.);
        TPZFNMatrix<16, STATE> v1f(4, 4, 0.);
        bc.ForcingFunctionBC()(dU.x, v2, v1f);
    }
    const TPZFMatrix<STATE> &v1 = bc.Val1();
    const bool penalty = (fMode == EFull);
    const REAL big = BigNumber();
    REAL u[3] = {0., 0., 0.};
    if (dU.sol.size() && dU.sol[0].size() >= dim)
        for (int c = 0; c < dim; c++) u[c] = dU.sol[0][c];
    const REAL p = (dP.sol.size() && dP.sol[0].size()) ? dP.sol[0][0] : 0.;

    auto traction = [&](const REAL t[3]) {
        for (int a = 0; a < nU; a++)
            for (int c = 0; c < dim; c++) ef(dim * a + c, 0) += t[c] * dU.phi.GetVal(a, 0) * weight;
    };
    auto penaltyU = [&](const bool mask[3], const REAL val[3]) {
        if (!penalty) return;
        for (int a = 0; a < nU; a++)
            for (int c = 0; c < dim; c++) {
                if (!mask[c]) continue;
                ef(dim * a + c, 0) += big * (val[c] - u[c]) * dU.phi.GetVal(a, 0) * weight;
                if (ek)
                    for (int b = 0; b < nU; b++)
                        (*ek)(dim * a + c, dim * b + c) += big * dU.phi.GetVal(a, 0) * dU.phi.GetVal(b, 0) * weight;
            }
    };
    auto penaltyP = [&](REAL val) {
        if (!penalty) return;
        for (int a = 0; a < nP; a++) {
            ef(offp + a, 0) += big * (val - p) * dP.phi.GetVal(a, 0) * weight;
            if (ek)
                for (int b = 0; b < nP; b++)
                    (*ek)(offp + a, offp + b) += big * dP.phi.GetVal(a, 0) * dP.phi.GetVal(b, 0) * weight;
        }
    };
    auto needP = [&]() {
        if (v2.size() < 4) {
            PZError << Name() << "::ContributeBC: o tipo " << bc.Type() << " exige Val2 com 4 componentes\n";
            DebugStop();
        }
    };
    const REAL vv[3] = {v2[0], v2.size() > 1 ? v2[1] : 0., v2.size() > 2 ? v2[2] : 0.};
    const REAL zero[3] = {0., 0., 0.};
    switch (bc.Type()) {
        case EDirichletU: {
            const bool all[3] = {true, true, true};
            penaltyU(all, vv);
        } break;
        case ETraction:
            traction(vv);
            break;
        case EDirichletP:
            needP();
            penaltyP(v2[3]);
            break;
        case EDirectionalNullU: {
            const bool mask[3] = {vv[0] != 0., vv[1] != 0., vv[2] != 0.};
            penaltyU(mask, zero);
        } break;
        case ENormalPressure: {
            const REAL t[3] = {v2[0] * dU.normal[0], v2[0] * dU.normal[1], v2[0] * dU.normal[2]};
            traction(t);
        } break;
        case EDirectionalU: {
            const bool mask[3] = {v1.GetVal(0, 0) != 0., v1.GetVal(1, 1) != 0., dim == 3 && v1.GetVal(2, 2) != 0.};
            penaltyU(mask, vv);
        } break;
        case ETractionDirichletP:
            needP();
            traction(vv);
            penaltyP(v2[3]);
            break;
        case EDirectionalUDirichletP: {
            needP();
            const bool mask[3] = {v1.GetVal(0, 0) != 0., v1.GetVal(1, 1) != 0., dim == 3 && v1.GetVal(2, 2) != 0.};
            penaltyU(mask, vv);
            penaltyP(v2[3]);
        } break;
        default:
            PZError << Name() << "::ContributeBC: tipo de condição de contorno desconhecido " << bc.Type() << "\n";
            DebugStop();
    }
}

template <class T, class TMEM>
void TPZMatPoroElastoPlastic3DMem<T, TMEM>::ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight,
                                                         TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef,
                                                         TPZBndCondT<STATE> &bc) {
    ContributeBCInternal(datavec, weight, &ek, ef, bc);
}

template <class T, class TMEM>
void TPZMatPoroElastoPlastic3DMem<T, TMEM>::ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight,
                                                         TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc) {
    ContributeBCInternal(datavec, weight, nullptr, ef, bc);
}

// ---------------------------------------------------------------------------------------------------
template <class T, class TMEM>
int TPZMatPoroElastoPlastic3DMem<T, TMEM>::VariableIndex(const std::string &name) const {
    if (name == "Displacement") return EDisplacement;
    if (name == "Pressure") return EPressure;
    if (name == "ExcessPressure") return EExcessPressure;
    if (name == "Flux") return EFlux;
    if (name == "PressureGradient") return EPressureGradient;
    return TBase::VariableIndex(name);
}

template <class T, class TMEM>
int TPZMatPoroElastoPlastic3DMem<T, TMEM>::NSolutionVariables(int var) const {
    switch (var) {
        case EDisplacement: return 3;
        case EPressure: return 1;
        case EExcessPressure: return 1;
        case EFlux: return 3;
        case EPressureGradient: return 3;
        default: return TBase::NSolutionVariables(var);
    }
}

template <class T, class TMEM>
void TPZMatPoroElastoPlastic3DMem<T, TMEM>::Solution(const TPZVec<TPZMaterialDataT<STATE>> &datavec, int var,
                                                     TPZVec<STATE> &Solout) {
    const TPZMaterialDataT<STATE> &dU = datavec[0], &dP = datavec[1];
    TPZFNMatrix<3, REAL> gp(3, 1, 0.);
    if (dP.dsol.size()) ToXYZ(dP.axes, dP.dsol[0], gp);
    switch (var) {
        case EDisplacement:
            Solout.Resize(3);
            Solout.Fill(0.);
            for (int c = 0; c < fDim; c++) Solout[c] = dU.sol[0][c];
            break;
        case EPressure:
            Solout.Resize(1);
            Solout[0] = dP.sol[0][0];
            break;
        case EExcessPressure:
            Solout.Resize(1);
            Solout[0] = dP.sol[0][0] - (fHydro ? fHydro(dP.x) : 0.);
            break;
        case EFlux:
            Solout.Resize(3);
            for (int c = 0; c < 3; c++) Solout[c] = c < fDim ? -(fK[c] / fMu) * (gp(c, 0) - fRhoF * fG[c]) : 0.;
            break;
        case EPressureGradient:
            Solout.Resize(3);
            for (int c = 0; c < 3; c++) Solout[c] = gp(c, 0);
            break;
        default:
            Solout.Resize(0);
    }
}

// ---------------------------------------------------------------------------------------------------
template <class T, class TMEM>
int TPZMatPoroElastoPlastic3DMem<T, TMEM>::ClassId() const {
    return Hash("TPZMatPoroElastoPlastic3DMem") ^ TBase::ClassId() << 1;
}

template <class T, class TMEM>
void TPZMatPoroElastoPlastic3DMem<T, TMEM>::Write(TPZStream &buf, int withclassid) const {
    TBase::Write(buf, withclassid);
    fPlasticModel.Write(buf, withclassid);
    buf.Write(&fDim);
    buf.Write(&fAlpha);
    buf.Write(&fSe);
    buf.Write(fK, 3);
    buf.Write(&fMu);
    buf.Write(&fRhoF);
    buf.Write(fG, 3);
    buf.Write(fBody, 3);
    buf.Write(&fTimeStep);
    int flow = fFlow, order = fIntegrationOrder;
    buf.Write(&flow);
    buf.Write(&order);
}

template <class T, class TMEM>
void TPZMatPoroElastoPlastic3DMem<T, TMEM>::Read(TPZStream &buf, void *context) {
    TBase::Read(buf, context);
    fPlasticModel.Read(buf, context);
    buf.Read(&fDim);
    buf.Read(&fAlpha);
    buf.Read(&fSe);
    buf.Read(fK, 3);
    buf.Read(&fMu);
    buf.Read(&fRhoF);
    buf.Read(fG, 3);
    buf.Read(fBody, 3);
    buf.Read(&fTimeStep);
    int flow, order;
    buf.Read(&flow);
    buf.Read(&order);
    fFlow = flow;
    fIntegrationOrder = order;
}

template <class T, class TMEM>
void TPZMatPoroElastoPlastic3DMem<T, TMEM>::Print(std::ostream &out) const {
    out << Name() << " id " << this->Id() << " dim " << fDim << "\n alpha = " << fAlpha << " Se = " << fSe
        << " k = (" << fK[0] << ", " << fK[1] << ", " << fK[2] << ") mu = " << fMu << " rhof = " << fRhoF << "\n g = (" << fG[0] << ", " << fG[1] << ", " << fG[2]
        << ") body = (" << fBody[0] << ", " << fBody[1] << ", " << fBody[2] << ")\n dt = " << fTimeStep
        << " flow = " << fFlow << " integration order = " << fIntegrationOrder << "\n";
    fPlasticModel.Print(out);
}

template class TPZMatPoroElastoPlastic3DMem<TPZModifiedCamClay, TPZElastoPlasticMem>;
template class TPZMatPoroElastoPlastic3DMem<TPZPlasticStepVoigt<TPZYCMohrCoulombPV2, TPZElasticResponse>,
                                            TPZElastoPlasticMem>;
