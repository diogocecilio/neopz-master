// Steady anisotropic Darcy flow in H1 (Ceron et al. 2025, Eq. 20-21): unknown u (excess pore pressure, kPa),
//     div v = s,   v = -K grad u,   K symmetric positive definite (constant),
// weak form  int grad w . K grad u = int w s + int_GN w g,  g = (K grad u) . n = inflow -v . n.
// TPZDarcyFlow of NeoPZ is isotropic only. Boundary conditions: type 0 Dirichlet u = val2[0] (penalty), type 1
// Neumann g = val2[0]; both read ForcingFunctionBC when it is set. Post-processing variables: ExcessPorePressure,
// Gradient, DarcyVelocity (-K grad u) and SeepageForce (-grad u, kN/m^3). Plane problem (dimension 2).
#ifndef ANISOTROPICDARCY_H
#define ANISOTROPICDARCY_H

#include "TPZBndCondT.h"
#include "TPZMatBase.h"
#include "TPZMatSingleSpace.h"
#include "TPZMaterialDataT.h"
#include "pzerror.h"

#include <cstring>
#include <string>

class TPZAnisotropicDarcy : public TPZMatBase<STATE, TPZMatSingleSpaceT<STATE>> {
    using TBase = TPZMatBase<STATE, TPZMatSingleSpaceT<STATE>>;

public:
    enum EBC { EDirichlet = 0, ENeumann = 1 };
    enum EVar { ENone = 0, EExcessPorePressure = 1, EGradient = 2, EDarcyVelocity = 3, ESeepageForce = 4 };

    TPZAnisotropicDarcy() = default;
    /// K = diag(kh, kv) in (x, y): horizontal and vertical hydraulic conductivities
    TPZAnisotropicDarcy(int id, REAL kh, REAL kv) : TBase(id) { SetPermeability(kh, 0., kv); }

    void SetPermeability(REAL kxx, REAL kxy, REAL kyy) {
        if (kxx <= 0. || kyy <= 0. || kxx * kyy - kxy * kxy <= 0.) DebugStop();
        fK[0][0] = kxx, fK[0][1] = fK[1][0] = kxy, fK[1][1] = kyy;
        // Dirichlet penalty relative to the conductivity (the default fBigNumber of TPZMaterial is ~ 7e16)
        SetBigNumber(1.e12 * std::max(kxx, kyy));
    }
    REAL K(int i, int j) const { return fK[i][j]; }

    std::string Name() const override { return "TPZAnisotropicDarcy"; }
    int Dimension() const override { return 2; }
    int NStateVariables() const override { return 1; }

    /// grad w in (x, y) from the derivatives in the element axes
    static void GradXY(const TPZFMatrix<REAL> &axes, const TPZFMatrix<REAL> &dax, int i, REAL g[2]) {
        g[0] = axes.GetVal(0, 0) * dax.GetVal(0, i) + axes.GetVal(1, 0) * dax.GetVal(1, i);
        g[1] = axes.GetVal(0, 1) * dax.GetVal(0, i) + axes.GetVal(1, 1) * dax.GetVal(1, i);
    }

    void Contribute(const TPZMaterialDataT<STATE> &data, REAL weight, TPZFMatrix<STATE> &ek,
                    TPZFMatrix<STATE> &ef) override {
        const int nphi = data.phi.Rows();
        STATE s = 0.;
        if (this->HasForcingFunction()) {
            TPZManVector<STATE, 1> res(1, 0.);
            this->fForcingFunction(data.x, res);
            s = res[0];
        }
        TPZManVector<REAL, 60> gx(nphi), gy(nphi);
        for (int i = 0; i < nphi; i++) {
            REAL g[2];
            GradXY(data.axes, data.dphix, i, g);
            gx[i] = g[0], gy[i] = g[1];
        }
        for (int i = 0; i < nphi; i++) {
            const REAL kgx = fK[0][0] * gx[i] + fK[0][1] * gy[i], kgy = fK[1][0] * gx[i] + fK[1][1] * gy[i];
            ef(i, 0) += weight * s * data.phi.GetVal(i, 0);
            for (int j = 0; j < nphi; j++) ek(i, j) += weight * (kgx * gx[j] + kgy * gy[j]);
        }
    }

    void ContributeBC(const TPZMaterialDataT<STATE> &data, REAL weight, TPZFMatrix<STATE> &ek,
                      TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc) override {
        const int nphi = data.phi.Rows();
        STATE v2 = bc.Val2()[0];
        if (bc.HasForcingFunctionBC()) {
            TPZManVector<STATE, 1> rhs(1, 0.);
            TPZFNMatrix<1, STATE> m(1, 1, 0.);
            bc.ForcingFunctionBC()(data.x, rhs, m);
            v2 = rhs[0];
        }
        switch (bc.Type()) {
        case EDirichlet: {
            const REAL big = this->BigNumber();
            for (int i = 0; i < nphi; i++) {
                ef(i, 0) += big * v2 * data.phi.GetVal(i, 0) * weight;
                for (int j = 0; j < nphi; j++) ek(i, j) += big * data.phi.GetVal(i, 0) * data.phi.GetVal(j, 0) * weight;
            }
            break;
        }
        case ENeumann:
            for (int i = 0; i < nphi; i++) ef(i, 0) += v2 * data.phi.GetVal(i, 0) * weight;
            break;
        default:
            PZError << "TPZAnisotropicDarcy: boundary condition type " << bc.Type() << " not implemented\n";
            DebugStop();
        }
    }

    int VariableIndex(const std::string &name) const override {
        if (name == "ExcessPorePressure") return EExcessPorePressure;
        if (name == "Gradient") return EGradient;
        if (name == "DarcyVelocity") return EDarcyVelocity;
        if (name == "SeepageForce") return ESeepageForce;
        return TBase::VariableIndex(name);
    }

    int NSolutionVariables(int var) const override {
        if (var == EExcessPorePressure) return 1;
        if (var == EGradient || var == EDarcyVelocity || var == ESeepageForce) return 3; // 3D vectors for VTK
        return TBase::NSolutionVariables(var);
    }

    void Solution(const TPZMaterialDataT<STATE> &data, int var, TPZVec<STATE> &sol) override {
        const TPZFMatrix<STATE> &dsol = data.dsol[0];
        if (var == EExcessPorePressure) {
            sol[0] = data.sol[0][0];
            return;
        }
        REAL g[2];
        GradXY(data.axes, dsol, 0, g);
        sol.Resize(3);
        sol[2] = 0.;
        switch (var) {
        case EGradient: sol[0] = g[0], sol[1] = g[1]; break;
        case EDarcyVelocity:
            sol[0] = -(fK[0][0] * g[0] + fK[0][1] * g[1]);
            sol[1] = -(fK[1][0] * g[0] + fK[1][1] * g[1]);
            break;
        case ESeepageForce: sol[0] = -g[0], sol[1] = -g[1]; break;
        default: DebugStop();
        }
    }

    void GetSolDimensions(uint64_t &u_len, uint64_t &du_row, uint64_t &du_col) const override {
        u_len = 1, du_row = 2, du_col = 1;
    }

    TPZMaterial *NewMaterial() const override { return new TPZAnisotropicDarcy(*this); }

    int ClassId() const override { return Hash("TPZAnisotropicDarcy") ^ TBase::ClassId() << 1; }

private:
    REAL fK[2][2] = {{1., 0.}, {0., 1.}};
};

#endif
