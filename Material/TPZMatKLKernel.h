#pragma once
#include "TPZMaterial.h"
#include "TPZMatSingleSpace.h"
#include "TPZMatGeneralisedEigenVal.h"
#include "pzinterpolationspace.h"
#include <functional>
#include <ostream>
#include "TPZMatBase.h"
#include "TPZMatErrorSingleSpace.h"
#include "TPZMatLoadCases.h"

class TPZMatKLKernel
: public TPZMatBase<STATE, TPZMatSingleSpaceT<STATE>, TPZMatGeneralisedEigenVal>
{
public:
    using TBase = TPZMatBase<STATE, TPZMatSingleSpaceT<STATE>, TPZMatGeneralisedEigenVal>;
    using KernelFn = std::function<STATE(const TPZVec<REAL>&, const TPZVec<REAL>&)>;
    enum class Target { A, B }; // A ≡ C (covariância), B ≡ massa

    // ---- tipos de kernel persistíveis ----
    enum class EKernelKind : int { Unknown = 0, ExpSeparable = 1, GaussSeparable = 2, User = 999 };

    // ---- construtores PM-safe ----
    TPZMatKLKernel();                                   // default seguro (dim=2, ExpSeparable Lx=1, Ly=0.1)
    TPZMatKLKernel(int matid, int dim);                 // default com id
    TPZMatKLKernel(int matid, int dim, REAL Lx, REAL Ly); // ExpSeparable com parâmetros

    // ---- TPZMatGeneralisedEigenVal ----
    void SetMatrixA() override { fTarget = Target::A; }
    void SetMatrixB() override { fTarget = Target::B; }
    Target CurrentTarget() const { return fTarget; }

    // ---- TPZMaterial (mínimo) ----
    int NStateVariables() const override { return 1; }
    int Dimension() const override { return fDim; }
    void Print(std::ostream &out) const override;

    // ---- MASSA L² (B) ----
    void Contribute(const TPZMaterialDataT<STATE>&,
                    REAL,
                    TPZFMatrix<STATE>&,
                    TPZFMatrix<STATE>&) override {}
    void ContributeBC(const TPZMaterialDataT<STATE>&,
                      REAL,
                      TPZFMatrix<STATE>&,
                      TPZFMatrix<STATE>&,
                      TPZBndCondT<STATE>&) override {}
    void FillDataRequirements(TPZMaterialData& data) const override;

    // ---- blocos locais usados pela struct matrix ----
    void CalcStiffNystrom(TPZInterpolationSpace* elx,
                          TPZInterpolationSpace* ely,
                          TPZElementMatrixT<STATE> &ce) const;
    void CalcStiffGalerkin(TPZInterpolationSpace* elx,
                           TPZInterpolationSpace* ely,
                           TPZElementMatrixT<STATE>& ce) const;
    void CalcStiffMass(TPZInterpolationSpace* el,
                       TPZElementMatrixT<STATE> &be,
                       int qmass) const;

    int  VariableIndex(const std::string &name) const override;
    int  NSolutionVariables(int var) const override;
    void Solution(const TPZMaterialDataT<STATE> &data, int var,
                  TPZVec<STATE> &Solout) override;

    enum EVars { ESolution=1, EGradient=2, EExact=3, EExactGrad=4, EError=5, EErrorGrad=6 };

    using TExact = std::function<void(const TPZVec<REAL>&, STATE&, TPZFMatrix<STATE>&)>;
    void SetExact(TExact f) { fExact = std::move(f); }

    // ---- IO / RTTI ----
    int  ClassId() const override;
    void Write(TPZStream &buf, int withclassid) const override; // grava dim, kind, params (não grava std::function)
    void Read (TPZStream &buf, void *context) override;         // lê dim, kind, params e chama RebuildKernel()
    bool HasForcingFunction() const override;

    // ---- helpers de kernel ----
    // custom (não serializado): use após Load() se precisar de forma arbitrária
    void SetKernel(KernelFn ker) { fKind = EKernelKind::User; fKernel = std::move(ker); }

    // persistíveis (serializados como tipo+parâmetros)
    void SetExpKernel(REAL Lx, REAL Ly);     // f(x,y)=exp(-|x0-y0|/Lx - |x1-y1|/Ly)
    void SetGaussKernel(REAL Lx, REAL Ly);   // f(x,y)=exp(-((dx/Lx)^2+(dy/Ly)^2))

    EKernelKind KernelKind()   const { return fKind; }
    const TPZManVector<REAL,3>& KernelParams() const { return fParams; }



private:
    void RebuildKernel(); // reconstrói fKernel a partir de fKind/fParams

    // estado
    KernelFn   fKernel;                 // não serializado diretamente
    int        fDim   = 2;
    Target     fTarget = Target::A;
    TExact     fExact;

    // descrição persistível do kernel
    EKernelKind          fKind   = EKernelKind::ExpSeparable;
    TPZManVector<REAL,3> fParams = TPZManVector<REAL,3>(3, 1.0); // [Lx,Ly,Lz(opc.)]
};
