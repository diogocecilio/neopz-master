/**
 * @file TPZMatPoroElastoPlasticUP.cpp
 * @brief Implementation of the coupled u-p poro-elastoplastic material (see TPZMatPoroElastoPlasticUP.h).
 */

#include "TPZMatPoroElastoPlasticUP.h"
#include "TPZPlasticStepModifiedCamClay.h"
#include "pzcmesh.h"
#include "pzcompel.h"
#include "pzmultiphysicselement.h"
#include "pzinterpolationspace.h"
#include "pzgeoel.h"
#include "pzquad.h"
#include "TPZHash.h"
#include "TPZStream.h"
#include "pzerror.h"
#include <cmath>

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

    // external forces only: body forces f_b (used for the normalization of the residual)
    if (fExternalOnly) {
        for (int a = 0; a < nu; ++a) {
            for (int d = 0; d < dim; ++d) ef(a * dim + d, 0) += dvol * dataU.phi.GetVal(a, 0) * fBodyForce[d];
        }
        return;
    }

    TPZFNMatrix<60, REAL> gradU, gradP, B;
    GlobalGradients(dataU, gradU);
    GlobalGradients(dataP, gradP);
    ComputeB(dataU.phi, gradU, x, B);

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

    // m^T B for each displacement dof (divergence including the hoop strain)
    TPZFNMatrix<60, REAL> mB(nequ, 1, 0.);
    for (int I = 0; I < nequ; ++I) mB(I, 0) = B(_XX_, I) + B(_YY_, I) + B(_ZZ_, I);

    // residual of the momentum balance: R_u = B^T sigma' - alpha m^T B p N_p - N_u^T b
    for (int a = 0; a < nu; ++a) {
        for (int d = 0; d < dim; ++d) {
            const int I = a * dim + d;
            REAL fint = 0.;
            for (int l = 0; l < 6; ++l) fint += B(l, I) * sigma[l];
            const REAL Ru = fint - fAlpha * mB(I, 0) * p - dataU.phi.GetVal(a, 0) * fBodyForce[d];
            ef(I, 0) -= dvol * Ru;
        }
    }
    // residual of the mass balance: R_p = N_p alpha (tr eps - tr eps_n) + N_p S (p - p_n) + Dt k gradN_p (grad p - rho_w g)
    const REAL treps = eps.I1();
    for (int j = 0; j < np; ++j) {
        const REAL Np = dataP.phi.GetVal(j, 0);
        REAL flow = 0.;
        for (int i = 0; i < dim; ++i) flow += gradP.GetVal(i, j) * (gradp[i] - fFluidWeight[i]);
        const REAL Rp = Np * fAlpha * (treps - trepsn) + Np * fInvBiotModulus * (p - pn) +
                        fTimeStep * fPermeability * flow;
        ef(nequ + j, 0) -= dvol * Rp;
    }

    if (!ek) return;
    // K_T = B^T D B
    TPZFNMatrix<360, REAL> DB(6, nequ, 0.);
    for (int l = 0; l < 6; ++l)
        for (int J = 0; J < nequ; ++J) {
            REAL val = 0.;
            for (int m = 0; m < 6; ++m) val += Dep.GetVal(l, m) * B.GetVal(m, J);
            DB(l, J) = val;
        }
    for (int I = 0; I < nequ; ++I) {
        for (int J = 0; J < nequ; ++J) {
            REAL val = 0.;
            for (int l = 0; l < 6; ++l) val += B.GetVal(l, I) * DB.GetVal(l, J);
            (*ek)(I, J) += dvol * val;
        }
        // -Q and Q^T
        for (int j = 0; j < np; ++j) {
            const REAL q = fAlpha * mB(I, 0) * dataP.phi.GetVal(j, 0) * dvol;
            (*ek)(I, nequ + j) -= q;
            (*ek)(nequ + j, I) += q;
        }
    }
    // S + Dt H
    for (int i = 0; i < np; ++i) {
        for (int j = 0; j < np; ++j) {
            REAL h = 0.;
            for (int k = 0; k < dim; ++k) h += gradP.GetVal(k, i) * gradP.GetVal(k, j);
            (*ek)(nequ + i, nequ + j) += dvol * (fInvBiotModulus * dataP.phi.GetVal(i, 0) * dataP.phi.GetVal(j, 0) +
                                                 fTimeStep * fPermeability * h);
        }
    }
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

template <class T, class TMEM>
int TPZMatPoroElastoPlasticUP<T, TMEM>::VariableIndex(const std::string &name) const {
    if (name == "Displacement") return EDisplacement;
    if (name == "PorePressure" || name == "Pressure") return EPorePressure;
    if (name == "DisplacementX") return EDisplacementX;
    if (name == "DisplacementY") return EDisplacementY;
    if (name == "DisplacementZ") return EDisplacementZ;
    return TBase::VariableIndex(name);
}

template <class T, class TMEM>
int TPZMatPoroElastoPlasticUP<T, TMEM>::NSolutionVariables(int var) const {
    switch (var) {
    case EDisplacement:
        return 3;
    case EPorePressure:
    case EDisplacementX:
    case EDisplacementY:
    case EDisplacementZ:
        return 1;
    default:
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
        break;
    case EPorePressure:
        sol.Resize(1);
        sol[0] = datavec[1].sol[0][0];
        break;
    case EDisplacementX:
    case EDisplacementY:
    case EDisplacementZ: {
        const int c = var - EDisplacementX;
        sol.Resize(1);
        sol[0] = (c < dim) ? u[c] : 0.;
        break;
    }
    default:
        sol.Resize(NSolutionVariables(var));
        sol.Fill(0.);
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
        TPZIntPoints *intrule = gel->CreateSideIntegrationRule(gel->NSides() - 1, order);
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
        delete intrule;
    }
}

template <class T, class TMEM>
void TPZMatPoroElastoPlasticUP<T, TMEM>::InitializeMemory(
    TPZCompMesh *mesh, const std::function<void(const TPZVec<REAL> &, TMEM &)> &init) {
    TMEM def = this->GetDefaultMemory();
    ForEachIntegrationPoint(mesh, [&](TPZCompEl *, int, const TPZVec<REAL> &x, const TPZVec<REAL> &, REAL, TMEM &mem) {
        mem = def;
        init(x, mem);
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
