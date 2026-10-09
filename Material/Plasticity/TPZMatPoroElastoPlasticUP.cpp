/**
 * @file TPZMatPoroElastoPlasticUP.cpp
 * @brief Implementation of the coupled u-p poro-elastoplastic material (see TPZMatPoroElastoPlasticUP.h).
 */

#include "TPZMatPoroElastoPlasticUP.h"
#include "TPZPlasticStepModifiedCamClay.h"
#include "TPZPlasticStepVoigt.h"
#include "TPZYCMohrCoulombPV2.h"
#include "TPZElasticResponse.h"
#include "pzcmesh.h"
#include "pzcompel.h"
#include "pzmultiphysicselement.h"
#include "pzinterpolationspace.h"
#include "pzgeoel.h"
#include "pzquad.h"
#include "TPZHash.h"
#include "TPZStream.h"
#include "pzerror.h"
#include <algorithm>
#include <cmath>
#include <map>
#include <memory>
#include <type_traits>

template <class T, class TMEM>
TPZMatPoroElastoPlasticUP<T, TMEM>::TPZMatPoroElastoPlasticUP()
    : TBase(), TPZMatPoroElastoPlasticUPBase(), fPlasticity(), fBodyForce(3, 0.), fFluidWeight(3, 0.) {
    fKinematics = EPlaneStrain;
}

template <class T, class TMEM>
TPZMatPoroElastoPlasticUP<T, TMEM>::TPZMatPoroElastoPlasticUP(int id, EKinematics kinematics)
    : TBase(id), TPZMatPoroElastoPlasticUPBase(), fPlasticity(), fBodyForce(3, 0.), fFluidWeight(3, 0.) {
    fKinematics = kinematics;
}

template <class T, class TMEM>
TPZMatPoroElastoPlasticUP<T, TMEM>::TPZMatPoroElastoPlasticUP(const TPZMatPoroElastoPlasticUP &cp)
    : TBase(cp), TPZMatPoroElastoPlasticUPBase(cp), fPlasticity(cp.fPlasticity), fAlpha(cp.fAlpha),
      fInvBiotModulus(cp.fInvBiotModulus), fPermeability(cp.fPermeability), fBodyForce(cp.fBodyForce),
      fFluidWeight(cp.fFluidWeight), fIntegrationOrder(cp.fIntegrationOrder) {
}

template <class T, class TMEM>
void TPZMatPoroElastoPlasticUP<T, TMEM>::SetPlasticModel(const T &model) {
    fPlasticity = model;
    TMEM mem;
    mem.m_elastoplastic_state = fPlasticity.GetState();
    mem.m_sigma.Zero();
    this->SetDefaultMem(mem);
}

template <class T, class TMEM>
void TPZMatPoroElastoPlasticUP<T, TMEM>::SetBodyForce(const TPZVec<REAL> &b) {
    fBodyForce.Fill(0.);
    for (int i = 0; i < b.size() && i < 3; ++i) fBodyForce[i] = b[i];
}

template <class T, class TMEM>
void TPZMatPoroElastoPlasticUP<T, TMEM>::SetFluidWeight(const TPZVec<REAL> &rhowg) {
    fFluidWeight.Fill(0.);
    for (int i = 0; i < rhowg.size() && i < 3; ++i) fFluidWeight[i] = rhowg[i];
}

template <class T, class TMEM>
int TPZMatPoroElastoPlasticUP<T, TMEM>::IntegrationRuleOrder(const TPZVec<int> &elPMaxOrder) const {
    if (fIntegrationOrder > 0) return fIntegrationOrder;
    return TPZMatCombinedSpacesT<STATE>::IntegrationRuleOrder(elPMaxOrder);
}

template <class T, class TMEM>
void TPZMatPoroElastoPlasticUP<T, TMEM>::FillDataRequirements(TPZVec<TPZMaterialDataT<STATE>> &datavec) const {
    for (int i = 0; i < datavec.size(); ++i) {
        datavec[i].SetAllRequirements(false);
        datavec[i].fNeedsSol = true;
    }
}

template <class T, class TMEM>
void TPZMatPoroElastoPlasticUP<T, TMEM>::FillBoundaryConditionDataRequirements(
    int type, TPZVec<TPZMaterialDataT<STATE>> &datavec) const {
    for (int i = 0; i < datavec.size(); ++i) {
        datavec[i].SetAllRequirements(false);
        datavec[i].fNeedsSol = true;
    }
}

template <class T, class TMEM>
void TPZMatPoroElastoPlasticUP<T, TMEM>::GlobalGradients(const TPZMaterialDataT<STATE> &data,
                                                         TPZFMatrix<REAL> &grad) const {
    const int dim = DimensionU();
    const int nshape = data.dphix.Cols();
    const int dimel = data.dphix.Rows();
    grad.Redim(dim, nshape);
    for (int a = 0; a < nshape; ++a) {
        for (int i = 0; i < dim; ++i) {
            REAL val = 0.;
            for (int k = 0; k < dimel; ++k) val += data.axes.GetVal(k, i) * data.dphix.GetVal(k, a);
            grad(i, a) = val;
        }
    }
}

template <class T, class TMEM>
void TPZMatPoroElastoPlasticUP<T, TMEM>::ComputeB(const TPZFMatrix<REAL> &phi, const TPZFMatrix<REAL> &G,
                                                  const TPZVec<REAL> &x, TPZFMatrix<REAL> &B) const {
    const int n = phi.Rows();
    if (fKinematics == EThreeDimensional) {
        B.Redim(6, 3 * n);
        for (int a = 0; a < n; ++a) {
            const int c = 3 * a;
            B(_XX_, c) = G.GetVal(0, a);
            B(_XY_, c) = G.GetVal(1, a);
            B(_XY_, c + 1) = G.GetVal(0, a);
            B(_XZ_, c) = G.GetVal(2, a);
            B(_XZ_, c + 2) = G.GetVal(0, a);
            B(_YY_, c + 1) = G.GetVal(1, a);
            B(_YZ_, c + 1) = G.GetVal(2, a);
            B(_YZ_, c + 2) = G.GetVal(1, a);
            B(_ZZ_, c + 2) = G.GetVal(2, a);
        }
    } else {
        B.Redim(6, 2 * n);
        for (int a = 0; a < n; ++a) {
            const int c = 2 * a;
            B(_XX_, c) = G.GetVal(0, a);
            B(_XY_, c) = G.GetVal(1, a);
            B(_XY_, c + 1) = G.GetVal(0, a);
            B(_YY_, c + 1) = G.GetVal(1, a);
            if (fKinematics == EAxisymmetric) B(_ZZ_, c) = phi.GetVal(a, 0) / x[0];
        }
    }
}

template <class T, class TMEM>
void TPZMatPoroElastoPlasticUP<T, TMEM>::ComputeN(const TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &N) const {
    const int dim = DimensionU();
    const int n = phi.Rows();
    N.Redim(dim * n, dim);
    for (int a = 0; a < n; ++a)
        for (int d = 0; d < dim; ++d) N(dim * a + d, d) = phi.GetVal(a, 0);
}

namespace {
/** @brief ek(row0 + i, col0 + j) += scale blk(i, j): adds a block to the element matrix (or vector) */
inline void AddBlock(TPZFMatrix<STATE> &ek, int64_t row0, int64_t col0, const TPZFMatrix<REAL> &blk, REAL scale) {
    const int64_t nr = blk.Rows(), nc = blk.Cols();
    for (int64_t j = 0; j < nc; ++j)
        for (int64_t i = 0; i < nr; ++i) ek(row0 + i, col0 + j) += scale * blk.GetVal(i, j);
}
} // namespace

template <class T, class TMEM>
void TPZMatPoroElastoPlasticUP<T, TMEM>::Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight,
                                                    TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) {
    ContributeInternal(datavec, weight, &ek, ef);
}

template <class T, class TMEM>
void TPZMatPoroElastoPlasticUP<T, TMEM>::Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight,
                                                    TPZFMatrix<STATE> &ef) {
    ContributeInternal(datavec, weight, nullptr, ef);
}

template <class T, class TMEM>
void TPZMatPoroElastoPlasticUP<T, TMEM>::ContributeInternal(const TPZVec<TPZMaterialDataT<STATE>> &datavec,
                                                            REAL weight, TPZFMatrix<STATE> *ek,
                                                            TPZFMatrix<STATE> &ef) {
    const TPZMaterialDataT<STATE> &dataU = datavec[0];
    const TPZMaterialDataT<STATE> &dataP = datavec[1];
    const int dim = DimensionU();
    const int nu = dataU.phi.Rows();
    const int np = dataP.phi.Rows();
    const int nequ = dim * nu;
    const TPZVec<REAL> &x = dataU.x;
    const bool axi = (fKinematics == EAxisymmetric);
    const REAL dvol = weight * (axi ? 2. * M_PI * x[0] : 1.);

    // displacement shape functions N_u^T (nequ x dim, as BuildBN of TPZMatElastoPlastic), body force b (dim x 1)
    // and f_b = N_u^T b (nequ x 1)
    TPZFNMatrix<180, REAL> NuT(nequ, dim);
    TPZFNMatrix<3, REAL> b(dim, 1);
    TPZFNMatrix<60, REAL> fb(nequ, 1);
    ComputeN(dataU.phi, NuT);
    for (int d = 0; d < dim; ++d) b(d, 0) = fBodyForce[d];
    NuT.Multiply(b, fb);                        // (nequ x dim) (dim x 1)

    // external forces only: ef_u += dvol f_b (used for the normalization of the residual)
    if (fExternalOnly) {
        AddBlock(ef, 0, 0, fb, dvol);
        return;
    }

    // strain-displacement operator B (6 x nequ, engineering shear strains, hoop strain in axisymmetry),
    // pressure shape functions N_p (1 x np; N_p^T is dataP.phi itself) and their global gradients G_p (dim x np)
    TPZFNMatrix<360, REAL> B(6, nequ), Bt(nequ, 6);
    TPZFNMatrix<60, REAL> gradU(dim, nu), Np(1, np), Gp(dim, np), Gpt(np, dim);
    GlobalGradients(dataU, gradU);
    GlobalGradients(dataP, Gp);
    ComputeB(dataU.phi, gradU, x, B);
    B.Transpose(&Bt);
    dataP.phi.Transpose(&Np);
    Gp.Transpose(&Gpt);
    const TPZFMatrix<REAL> &Npt = dataP.phi;

    // total strain at the point (engineering components), from the current solution
    TPZFNMatrix<9, REAL> du(3, 3, 0.);
    {
        const TPZFMatrix<STATE> &dsol = dataU.dsol[0];
        const int dimel = dsol.Rows();
        for (int i = 0; i < dim; ++i) {
            for (int c = 0; c < dim; ++c) {
                REAL val = 0.;
                for (int k = 0; k < dimel; ++k) val += dataU.axes.GetVal(k, i) * dsol.GetVal(k, c);
                du(c, i) = val; // du(c,i) = d u_c / d x_i
            }
        }
    }
    // engineering shear strains (convention of the strain tensors of the plastic models)
    TPZTensor<REAL> eps;
    eps.XX() = du(0, 0);
    eps.YY() = du(1, 1);
    eps.XY() = du(0, 1) + du(1, 0);
    if (fKinematics == EThreeDimensional) {
        eps.ZZ() = du(2, 2);
        eps.XZ() = du(0, 2) + du(2, 0);
        eps.YZ() = du(1, 2) + du(2, 1);
    } else if (axi) {
        eps.ZZ() = dataU.sol[0][0] / x[0];
    }

    // pore pressure and its gradient
    const REAL p = dataP.sol[0][0];
    TPZManVector<REAL, 3> gradp(3, 0.);
    {
        const TPZFMatrix<STATE> &dsolp = dataP.dsol[0];
        const int dimel = dsolp.Rows();
        for (int i = 0; i < dim; ++i) {
            REAL val = 0.;
            for (int k = 0; k < dimel; ++k) val += dataP.axes.GetVal(k, i) * dsolp.GetVal(k, 0);
            gradp[i] = val;
        }
    }

    // stress update from the converged state of the integration point
    const int64_t gp = dataU.intGlobPtIndex;
    if (gp < 0) DebugStop();
    TMEM &mem = this->MemItem(gp);
    T model(fPlasticity);
    model.SetState(mem.m_elastoplastic_state);
    TPZTensor<REAL> sigma(mem.m_sigma);
    TPZFNMatrix<36, REAL> Dep(6, 6, 0.);
    model.ApplyStrainComputeSigma(eps, sigma, ek ? &Dep : nullptr);
    if (model.LastProjectionFailed()) fNFailed++;
    const REAL trepsn = mem.m_elastoplastic_state.m_eps_t.I1();
    const REAL pn = mem.m_elastoplastic_state.fpressure;
    if (this->fUpdateMem) {
        mem.m_sigma = sigma;
        mem.m_elastoplastic_state = model.GetState();
        mem.m_elastoplastic_state.fpressure = p;
    }

    // Voigt vectors (order of TPZTensor: XX XY XZ YY YZ ZZ): effective stress sigma' and m = (1 0 0 1 0 1)^T;
    // B^T m (nequ x 1) is the divergence of the displacement shape functions (hoop strain included)
    TPZFNMatrix<6, REAL> sig(6, 1), m(6, 1, 0.);
    sigma.CopyTo(sig);
    m(_XX_, 0) = m(_YY_, 0) = m(_ZZ_, 0) = 1.;
    TPZFNMatrix<60, REAL> Btm(nequ, 1), mB(1, nequ);
    Bt.Multiply(m, Btm);                        // B^T m  (nequ x 6) (6 x 1)
    Btm.Transpose(&mB);                         // m^T B  (1 x nequ)

    // momentum balance: R_u = B^T sigma' - alpha p B^T m - N_u^T b;  ef_u -= dvol R_u
    TPZFNMatrix<60, REAL> Ru(nequ, 1), Qp(Btm);
    Bt.Multiply(sig, Ru);                       // B^T sigma'  (nequ x 6) (6 x 1)
    Qp *= fAlpha * p;                           // alpha p B^T m
    Ru -= Qp;
    Ru -= fb;                                   // N_u^T b
    AddBlock(ef, 0, 0, Ru, -dvol);

    // mass balance: R_p = N_p^T [alpha (tr eps - tr eps_n) + (p - p_n)/M] + Dt k G_p^T (grad p - rho_w g);
    // ef_p -= dvol R_p
    TPZFNMatrix<3, REAL> flow(dim, 1);          // grad p - rho_w g
    for (int i = 0; i < dim; ++i) flow(i, 0) = gradp[i] - fFluidWeight[i];
    TPZFNMatrix<60, REAL> Rp(Npt), q(np, 1);
    Rp *= fAlpha * (eps.I1() - trepsn) + fInvBiotModulus * (p - pn); // N_p^T [alpha (tr eps - tr eps_n) + (p - p_n)/M]
    Gpt.Multiply(flow, q);                      // G_p^T (grad p - rho_w g)  (np x dim) (dim x 1)
    q *= fTimeStep * fPermeability;
    Rp += q;
    AddBlock(ef, nequ, 0, Rp, -dvol);

    if (!ek) return;
    // tangent: ek += dvol [[B^T Dep B, -Q], [Q^T, S + Dt H]], Q = alpha B^T m N_p
    TPZFNMatrix<360, REAL> DB(6, nequ);
    TPZFNMatrix<3600, REAL> Kuu(nequ, nequ);
    Dep.Multiply(B, DB);                        // Dep B      (6 x 6) (6 x nequ)
    Bt.Multiply(DB, Kuu);                       // B^T Dep B  (nequ x 6) (6 x nequ)
    AddBlock(*ek, 0, 0, Kuu, dvol);             // ek_uu += dvol B^T Dep B

    TPZFNMatrix<480, REAL> Q(nequ, np), Qt(np, nequ);
    Btm.Multiply(Np, Q);                        // B^T m N_p    (nequ x 1) (1 x np)
    Npt.Multiply(mB, Qt);                       // N_p^T m^T B  (np x 1) (1 x nequ)
    AddBlock(*ek, 0, nequ, Q, -fAlpha * dvol);  // ek_up -= Q
    AddBlock(*ek, nequ, 0, Qt, fAlpha * dvol);  // ek_pu += Q^T

    TPZFNMatrix<64, REAL> S(np, np), H(np, np);
    Npt.Multiply(Np, S);                        // N_p^T N_p  (np x 1) (1 x np)
    S *= fInvBiotModulus;
    Gpt.Multiply(Gp, H);                        // G_p^T G_p  (np x dim) (dim x np)
    H *= fTimeStep * fPermeability;
    S += H;                                     // S + Dt H
    AddBlock(*ek, nequ, nequ, S, dvol);         // ek_pp += dvol (S + Dt H)
}

template <class T, class TMEM>
void TPZMatPoroElastoPlasticUP<T, TMEM>::ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight,
                                                      TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef,
                                                      TPZBndCondT<STATE> &bc) {
    ContributeBCInternal(datavec, weight, &ek, ef, bc);
}

template <class T, class TMEM>
void TPZMatPoroElastoPlasticUP<T, TMEM>::ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight,
                                                      TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc) {
    ContributeBCInternal(datavec, weight, nullptr, ef, bc);
}

template <class T, class TMEM>
void TPZMatPoroElastoPlasticUP<T, TMEM>::ContributeBCInternal(const TPZVec<TPZMaterialDataT<STATE>> &datavec,
                                                              REAL weight, TPZFMatrix<STATE> *ek,
                                                              TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc) {
    const TPZMaterialDataT<STATE> &dataU = datavec[0];
    const TPZMaterialDataT<STATE> &dataP = datavec[1];
    const int dim = DimensionU();
    const int nu = dataU.phi.Rows();
    const int np = dataP.phi.Rows();
    const int nequ = dim * nu;
    const int type = bc.Type();

    TPZManVector<STATE, 3> v2(bc.Val2());
    if (v2.size() < 3) {
        const int n0 = v2.size();
        v2.Resize(3);
        for (int i = n0; i < 3; ++i) v2[i] = 0.;
    }
    if (bc.HasForcingFunctionBC()) {
        TPZFNMatrix<9, STATE> v1(3, 3, 0.);
        bc.ForcingFunctionBC()(dataU.x, v2, v1);
    }

    switch (type) {
    case ENeumannU:
    case ENeumannUFixed: {
        const REAL factor = (type == ENeumannU) ? fLoadFactor : 1.;
        const REAL w = weight * (fKinematics == EAxisymmetric ? 2. * M_PI * dataU.x[0] : 1.);
        for (int a = 0; a < nu; ++a)
            for (int d = 0; d < dim; ++d) ef(a * dim + d, 0) += w * factor * v2[d] * dataU.phi.GetVal(a, 0);
        break;
    }
    case EDirichletU:
    case EDirichletUDirectional: {
        if (fExternalOnly) break;
        for (int d = 0; d < dim; ++d) {
            if (type == EDirichletUDirectional && bc.Val1().GetVal(d, d) == 0.) continue;
            const REAL u = dataU.sol[0][d];
            for (int a = 0; a < nu; ++a) {
                ef(a * dim + d, 0) += fBig * (v2[d] - u) * dataU.phi.GetVal(a, 0) * weight;
                if (!ek) continue;
                for (int b = 0; b < nu; ++b)
                    (*ek)(a * dim + d, b * dim + d) += fBig * dataU.phi.GetVal(a, 0) * dataU.phi.GetVal(b, 0) * weight;
            }
        }
        break;
    }
    case EDirichletP: {
        if (fExternalOnly) break;
        const REAL p = dataP.sol[0][0];
        for (int i = 0; i < np; ++i) {
            ef(nequ + i, 0) += fBig * (v2[0] - p) * dataP.phi.GetVal(i, 0) * weight;
            if (!ek) continue;
            for (int j = 0; j < np; ++j)
                (*ek)(nequ + i, nequ + j) += fBig * dataP.phi.GetVal(i, 0) * dataP.phi.GetVal(j, 0) * weight;
        }
        break;
    }
    default:
        std::cout << __PRETTY_FUNCTION__ << " unknown boundary condition type " << type << std::endl;
        DebugStop();
    }
}

namespace {
/**
 * @brief Eigenvalues of a symmetric tensor (cyclic Jacobi rotations, accurate also for repeated
 * eigenvalues), sorted in decreasing order
 */
void SortedEigenvalues(const TPZTensor<REAL> &t, REAL eig[3]) {
    REAL a[3][3] = {{t.XX(), t.XY(), t.XZ()}, {t.XY(), t.YY(), t.YZ()}, {t.XZ(), t.YZ(), t.ZZ()}};
    for (int sweep = 0; sweep < 50; ++sweep) {
        const REAL off = std::fabs(a[0][1]) + std::fabs(a[0][2]) + std::fabs(a[1][2]);
        const REAL diag = std::fabs(a[0][0]) + std::fabs(a[1][1]) + std::fabs(a[2][2]);
        if (off <= 1.e-15 * diag || off < 1.e-300) break;
        for (int p = 0; p < 2; ++p) {
            for (int q = p + 1; q < 3; ++q) {
                if (a[p][q] == 0.) continue;
                // rotation that annihilates a[p][q]
                const REAL theta = (a[q][q] - a[p][p]) / (2. * a[p][q]);
                const REAL tt = (theta >= 0. ? 1. : -1.) / (std::fabs(theta) + std::sqrt(theta * theta + 1.));
                const REAL c = 1. / std::sqrt(tt * tt + 1.), sn = tt * c;
                for (int k = 0; k < 3; ++k) { // columns p and q
                    const REAL akp = a[k][p], akq = a[k][q];
                    a[k][p] = c * akp - sn * akq;
                    a[k][q] = sn * akp + c * akq;
                }
                for (int k = 0; k < 3; ++k) { // rows p and q
                    const REAL apk = a[p][k], aqk = a[q][k];
                    a[p][k] = c * apk - sn * aqk;
                    a[q][k] = sn * apk + c * aqk;
                }
            }
        }
    }
    for (int i = 0; i < 3; ++i) eig[i] = a[i][i];
    std::sort(eig, eig + 3, [](REAL x, REAL y) { return x > y; });
}

/** @brief Components of a tensor by rows (XX XY XZ, XY YY YZ, XZ YZ ZZ) */
void TensorRows(const TPZTensor<REAL> &t, TPZVec<STATE> &sol) {
    sol.Resize(9);
    sol[0] = t.XX(); sol[1] = t.XY(); sol[2] = t.XZ();
    sol[3] = t.XY(); sol[4] = t.YY(); sol[5] = t.YZ();
    sol[6] = t.XZ(); sol[7] = t.YZ(); sol[8] = t.ZZ();
}
} // namespace

template <class T, class TMEM>
int TPZMatPoroElastoPlasticUP<T, TMEM>::VariableIndex(const std::string &name) const {
    static const std::map<std::string, int> names = {
        {"Displacement", EDisplacement},
        {"PorePressure", EPorePressure},
        {"Pressure", EPorePressure},
        {"DisplacementX", EDisplacementX},
        {"DisplacementY", EDisplacementY},
        {"DisplacementZ", EDisplacementZ},
        {"MeanEffectiveStress", EMeanEffectiveStress},
        {"DeviatoricStress", EDeviatoricStress},
        {"EffectiveStressXX", EEffectiveStressXX},
        {"EffectiveStressYY", EEffectiveStressYY},
        {"EffectiveStressZZ", EEffectiveStressZZ},
        {"EffectiveStressXY", EEffectiveStressXY},
        {"EffectiveStressXZ", EEffectiveStressXZ},
        {"EffectiveStressYZ", EEffectiveStressYZ},
        {"PrincipalEffectiveStress", EPrincipalEffectiveStress},
        {"PreconsolidationPressure", EPreconsolidationPressure},
        {"PlasticType", EPlasticType},
        {"VolumetricStrain", EVolumetricStrain},
        {"SpecificVolume", ESpecificVolume},
        {"TotalStressXX", ETotalStressXX},
        {"TotalStressYY", ETotalStressYY},
        {"TotalStressZZ", ETotalStressZZ},
        {"TotalStressXY", ETotalStressXY},
        {"TotalStressXZ", ETotalStressXZ},
        {"TotalStressYZ", ETotalStressYZ},
        {"EffectiveStress", EEffectiveStress},
        {"TotalStress", ETotalStress}};
    auto it = names.find(name);
    if (it != names.end()) return it->second;
    return TBase::VariableIndex(name);
}

template <class T, class TMEM>
int TPZMatPoroElastoPlasticUP<T, TMEM>::NSolutionVariables(int var) const {
    switch (var) {
    case EDisplacement:
    case EPrincipalEffectiveStress:
        return 3;
    case EEffectiveStress:
    case ETotalStress:
        return 9;
    default:
        if (var >= EPorePressure && var <= ETotalStressYZ) return 1;
        return TBase::NSolutionVariables(var);
    }
}

template <class T, class TMEM>
void TPZMatPoroElastoPlasticUP<T, TMEM>::Solution(const TPZVec<TPZMaterialDataT<STATE>> &datavec, int var,
                                                  TPZVec<STATE> &sol) {
    const int dim = DimensionU();
    const TPZVec<STATE> &u = datavec[0].sol[0];
    switch (var) {
    case EDisplacement:
        sol.Resize(3);
        for (int i = 0; i < 3; ++i) sol[i] = (i < dim) ? u[i] : 0.;
        return;
    case EPorePressure:
        sol.Resize(1);
        sol[0] = datavec[1].sol[0][0];
        return;
    case EDisplacementX:
    case EDisplacementY:
    case EDisplacementZ: {
        const int c = var - EDisplacementX;
        sol.Resize(1);
        sol[0] = (c < dim) ? u[c] : 0.;
        return;
    }
    default:
        break;
    }
    if (var < EMeanEffectiveStress || var > ETotalStress) {
        sol.Resize(NSolutionVariables(var));
        sol.Fill(0.);
        return;
    }
    // variables of the integration points: state of the last converged step in the memory
    sol.Resize(NSolutionVariables(var));
    sol.Fill(0.);
    const int64_t gp = datavec[0].intGlobPtIndex;
    if (gp < 0 || gp >= int64_t(this->GetMemory()->NElements())) return;
    const TMEM &mem = this->MemItem(gp);
    const TPZTensor<REAL> &sig = mem.m_sigma;
    const auto &state = mem.m_elastoplastic_state;
    // total stress with the pore pressure of the point (compression positive)
    TPZTensor<REAL> total(sig);
    for (int i : {_XX_, _YY_, _ZZ_}) total[i] -= fAlpha * state.fpressure;
    // Voigt index of TPZTensor of the components XX, YY, ZZ, XY, XZ, YZ
    static const int comp[6] = {_XX_, _YY_, _ZZ_, _XY_, _XZ_, _YZ_};
    switch (var) {
    case EMeanEffectiveStress:
        sol[0] = -sig.I1() / 3.;
        break;
    case EDeviatoricStress:
        sol[0] = std::sqrt(3. * std::max(REAL(0.), sig.J2()));
        break;
    case EPrincipalEffectiveStress: {
        REAL eig[3];
        SortedEigenvalues(sig, eig);
        for (int i = 0; i < 3; ++i) sol[i] = eig[i];
        break;
    }
    case EPreconsolidationPressure:
        sol[0] = state.m_hardening;
        break;
    case EPlasticType:
        sol[0] = state.m_m_type;
        break;
    case EVolumetricStrain:
        sol[0] = state.m_eps_t.I1();
        break;
    case ESpecificVolume:
        sol[0] = state.fmatprop.size() ? state.fmatprop[0] : 0.;
        break;
    case EEffectiveStress:
        TensorRows(sig, sol);
        break;
    case ETotalStress:
        TensorRows(total, sol);
        break;
    default:
        if (var >= EEffectiveStressXX && var <= EEffectiveStressYZ) sol[0] = sig[comp[var - EEffectiveStressXX]];
        else if (var >= ETotalStressXX && var <= ETotalStressYZ) sol[0] = total[comp[var - ETotalStressXX]];
        break;
    }
}

template <class T, class TMEM>
void TPZMatPoroElastoPlasticUP<T, TMEM>::ForEachIntegrationPoint(
    TPZCompMesh *mesh, const std::function<void(TPZCompEl *, int, const TPZVec<REAL> &, const TPZVec<REAL> &, REAL,
                                                TMEM &)> &f) {
    const int64_t nel = mesh->NElements();
    for (int64_t iel = 0; iel < nel; ++iel) {
        TPZCompEl *cel = mesh->Element(iel);
        if (!cel || !cel->Material() || cel->Material()->Id() != this->Id()) continue;
        TPZMultiphysicsElement *mfel = dynamic_cast<TPZMultiphysicsElement *>(cel);
        if (!mfel) continue;
        TPZGeoEl *gel = cel->Reference();
        const int dim = gel->Dimension();
        // same integration rule of TPZMultiphysicsCompEl::CalcStiff
        TPZManVector<int, 4> ordervec;
        for (int iref = 0; iref < mfel->NMeshes(); ++iref) {
            TPZInterpolationSpace *msp = dynamic_cast<TPZInterpolationSpace *>(mfel->Element(iref));
            if (!msp) continue;
            ordervec.Resize(ordervec.size() + 1);
            ordervec[ordervec.size() - 1] = msp->MaxOrder();
        }
        const int order = IntegrationRuleOrder(ordervec);
        std::unique_ptr<TPZIntPoints> intrule(gel->CreateSideIntegrationRule(gel->NSides() - 1, order));
        TPZManVector<int, 4> intorder(dim, order);
        intrule->SetOrder(intorder);
        TPZManVector<int64_t, 64> indices;
        cel->GetMemoryIndices(indices);
        const int npts = intrule->NPoints();
        if (indices.size() != npts) {
            std::cout << __PRETTY_FUNCTION__ << " element " << iel << " has " << indices.size()
                      << " memory items and " << npts << " integration points" << std::endl;
            DebugStop();
        }
        TPZManVector<REAL, 3> qsi(dim, 0.), x(3, 0.);
        TPZFNMatrix<9, REAL> jac, axes, jacinv;
        REAL detjac, w;
        for (int ip = 0; ip < npts; ++ip) {
            intrule->Point(ip, qsi, w);
            gel->Jacobian(qsi, jac, axes, detjac, jacinv);
            gel->X(qsi, x);
            f(cel, ip, x, qsi, w * std::fabs(detjac), this->MemItem(indices[ip]));
        }
    }
}

template <class T, class TMEM>
void TPZMatPoroElastoPlasticUP<T, TMEM>::InitializeMemory(
    TPZCompMesh *mesh, const std::function<void(const TPZVec<REAL> &, TMEM &)> &init) {
    TMEM def = this->GetDefaultMemory();
    ForEachIntegrationPoint(mesh, [&](TPZCompEl *, int, const TPZVec<REAL> &x, const TPZVec<REAL> &, REAL, TMEM &mem) {
        mem = def;
        init(x, mem);
        // a total-strain plastic step (TPZPlasticStepVoigt) does not read m_sigma: the initial effective stress
        // becomes the eigenstrain eps_p = eps_t - De^-1 sigma0, so that De (eps - eps_p) = sigma0 at eps = eps_t
        if constexpr (!std::is_same_v<T, TPZPlasticStepModifiedCamClay>) {
            TPZTensor<REAL> epse;
            fPlasticity.GetElasticResponse().ComputeStrain(mem.m_sigma, epse);
            mem.m_elastoplastic_state.m_eps_p = mem.m_elastoplastic_state.m_eps_t - epse;
        }
    });
}

template <class T, class TMEM>
TPZMaterial *TPZMatPoroElastoPlasticUP<T, TMEM>::NewMaterial() const {
    return new TPZMatPoroElastoPlasticUP<T, TMEM>(*this);
}

template <class T, class TMEM>
void TPZMatPoroElastoPlasticUP<T, TMEM>::Print(std::ostream &out) const {
    out << Name() << " id = " << this->Id() << "\n kinematics = "
        << (fKinematics == EPlaneStrain ? "plane strain" : fKinematics == EAxisymmetric ? "axisymmetric" : "3D")
        << "\n alpha = " << fAlpha << " 1/M = " << fInvBiotModulus << " k = " << fPermeability
        << "\n body force = " << fBodyForce << " fluid weight = " << fFluidWeight
        << "\n time step = " << fTimeStep << " load factor = " << fLoadFactor
        << " integration order = " << fIntegrationOrder << " big = " << fBig << "\n";
    fPlasticity.Print(out);
}

template <class T, class TMEM>
int TPZMatPoroElastoPlasticUP<T, TMEM>::ClassId() const {
    return Hash("TPZMatPoroElastoPlasticUP") ^ TBase::ClassId() << 1 ^ fPlasticity.ClassId() << 2;
}

template <class T, class TMEM>
void TPZMatPoroElastoPlasticUP<T, TMEM>::Write(TPZStream &buf, int withclassid) const {
    TBase::Write(buf, withclassid);
    int kin = fKinematics;
    buf.Write(&kin);
    buf.Write(&fTimeStep);
    buf.Write(&fLoadFactor);
    buf.Write(&fBig);
    buf.Write(&fAlpha);
    buf.Write(&fInvBiotModulus);
    buf.Write(&fPermeability);
    buf.Write(fBodyForce);
    buf.Write(fFluidWeight);
    buf.Write(&fIntegrationOrder);
    fPlasticity.Write(buf, withclassid);
}

template <class T, class TMEM>
void TPZMatPoroElastoPlasticUP<T, TMEM>::Read(TPZStream &buf, void *context) {
    TBase::Read(buf, context);
    int kin;
    buf.Read(&kin);
    fKinematics = EKinematics(kin);
    buf.Read(&fTimeStep);
    buf.Read(&fLoadFactor);
    buf.Read(&fBig);
    buf.Read(&fAlpha);
    buf.Read(&fInvBiotModulus);
    buf.Read(&fPermeability);
    buf.Read(fBodyForce);
    buf.Read(fFluidWeight);
    buf.Read(&fIntegrationOrder);
    fPlasticity.Read(buf, context);
}

template class TPZMatPoroElastoPlasticUP<TPZPlasticStepModifiedCamClay, TPZElastoPlasticMem>;
template class TPZRestoreClass<TPZMatPoroElastoPlasticUP<TPZPlasticStepModifiedCamClay, TPZElastoPlasticMem>>;
// Mohr-Coulomb (closed-form RHW projection, paper model): total-strain update, the stress of the memory is not used
template class TPZMatPoroElastoPlasticUP<TPZPlasticStepVoigt<TPZYCMohrCoulombPV2, TPZElasticResponse>, TPZElastoPlasticMem>;
template class TPZRestoreClass<TPZMatPoroElastoPlasticUP<TPZPlasticStepVoigt<TPZYCMohrCoulombPV2, TPZElasticResponse>, TPZElastoPlasticMem>>;
