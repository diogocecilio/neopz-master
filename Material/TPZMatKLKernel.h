#pragma once
#include "TPZMaterial.h"
#include "TPZMatSingleSpace.h"
#include "TPZMatGeneralisedEigenVal.h"
#include "pzinterpolationspace.h"
#include <functional>
#include <ostream>
#include "TPZMatBase.h"
#include "TPZMatSingleSpace.h"
#include "TPZMatErrorSingleSpace.h"
#include "TPZMatLoadCases.h"
// fwd: versão *templated* do ElementMatrix é a que costuma existir hoje no PZ
//template<class TVar> class TPZElementMatrixT;

class TPZMatKLKernel :     public TPZMatBase<STATE,TPZMatSingleSpaceT<STATE>,TPZMatGeneralisedEigenVal>{
public:
    using TBase = TPZMatBase<STATE,
    TPZMatSingleSpaceT<STATE>,
    TPZMatGeneralisedEigenVal>;
    using KernelFn = std::function<STATE(const TPZVec<REAL>&, const TPZVec<REAL>&)>;
    enum class Target { A, B }; // A ≡ C (covariância), B ≡ massa

    TPZMatKLKernel(int matid, int dim, KernelFn ker/*, int qx, int qy*/);

    // ---- TPZMatGeneralisedEigenVal ----
    void SetMatrixA() override { fTarget = Target::A; }
    void SetMatrixB() override { fTarget = Target::B; }
    Target CurrentTarget() const { return fTarget; }

    // ---- TPZMaterial (mínimo) ----
    int NStateVariables() const  override { return 1; }
    int Dimension() const  override{ return fDim; }
    void Print(std::ostream &out) const override;

    // ---- MASSA L² (B) ----
    void Contribute(const TPZMaterialDataT<STATE>& data,
                    REAL weight,
                    TPZFMatrix<STATE>& ek,
                    TPZFMatrix<STATE>& ef) override {}


                    void ContributeBC(const TPZMaterialDataT<STATE>&,
                                      REAL,
                                      TPZFMatrix<STATE>&,
                                      TPZFMatrix<STATE>&,
                                      TPZBndCondT<STATE>&) override {}
    void FillDataRequirements(TPZMaterialData&data)const override;

    // ---- blocos locais usados pela struct matrix ----
    void CalcStiffNystrom(TPZInterpolationSpace* elx,
                          TPZInterpolationSpace* ely,
                          TPZElementMatrixT<STATE> &ce) const;

    void CalcStiffMass(TPZInterpolationSpace* el,
                       TPZElementMatrixT<STATE> &be,
                       int qmass) const;


    int  VariableIndex(const std::string &name) const override;
    int  NSolutionVariables(int var) const override;
    void Solution(const TPZMaterialDataT<STATE> &data, int var,
                TPZVec<STATE> &Solout) override;

    enum EVars { ESolution=1, EGradient=2, EExact=3, EExactGrad=4,
               EError=5, EErrorGrad=6 };

    using TExact = std::function<void(const TPZVec<REAL>&, STATE&, TPZFMatrix<STATE>&)>;
    void SetExact(TExact f) { fExact = std::move(f); }


    // ---- IO / RTTI ----
    int ClassId() const override;
    void Write(TPZStream &buf, int withclassid) const override;
    void Read(TPZStream &buf, void *context) override;
    bool HasForcingFunction() const override;


    // ---- helpers ----
    void SetKernel(KernelFn ker) { fKernel = std::move(ker); }
    //void SetOrders(int qx, int qy) { fQx = qx; fQy = qy; }

private:
    KernelFn fKernel;
    int fDim = 2;
    //int fQx  = 2;   // ordem para integração em x (C)
    //int fQy  = 2;   // ordem para integração em y (C)
    Target fTarget = Target::A;
    TExact fExact;
};
