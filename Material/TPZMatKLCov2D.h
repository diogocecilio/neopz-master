#pragma once

#include <functional>
#include <vector>
#include <string>
#include <algorithm>

#include "TPZMatBase.h"
#include "TPZMatSingleSpace.h"

#include "TPZMatErrorSingleSpace.h"

#include "TPZMatLoadCases.h"

#include "pzfmatrix.h"
#include "pzvec.h"
#include "pzmanvector.h"

class TPZCompMesh;
class TPZInterpolationSpace;
template<class TVar> class TPZMaterialDataT;
template<class TVar> class TPZBndCondT;

/// Material L² + utilitários Nyström p/ KL (C = Φᵀ W K W Φ)
class TPZMatKLCov2D
: public TPZMatBase<
      STATE,
      TPZMatSingleSpaceT<STATE>,
      TPZMatErrorSingleSpace<STATE>,
      TPZMatLoadCases<STATE>>
{
public:
    using TBase = TPZMatBase<
        STATE,
        TPZMatSingleSpaceT<STATE>,
        TPZMatErrorSingleSpace<STATE>,
        TPZMatLoadCases<STATE>>;

    using KernelType = std::function<STATE(const TPZVec<REAL>&,
                                           const TPZVec<REAL>&)>;

    // ---- ctors ----
    TPZMatKLCov2D(int id, int dim = 2);
    TPZMatKLCov2D() : TPZMatKLCov2D(-1,2) {}

    // ---- metadata ----
    std::string Name() const override { return "TPZMatKLCov2D"; }
    int  Dimension() const override   { return fDim; }
    int  NStateVariables() const override { return 1; }
    //bool IsMatImpl() override         { return true; }

    // ---- MASSA L² (B) ----
    void Contribute(const TPZMaterialDataT<STATE>& data,
                    REAL weight,
                    TPZFMatrix<STATE>& ek,
                    TPZFMatrix<STATE>& ef) override;


    void ContributeBC(const TPZMaterialDataT<STATE>&,
                      REAL,
                      TPZFMatrix<STATE>&,
                      TPZFMatrix<STATE>&,
                      TPZBndCondT<STATE>&) override {}

    // ---- pós-processo ----
    int  VariableIndex(const std::string& name) const override;
    int  NSolutionVariables(int var) const override;
    void Solution(const TPZMaterialDataT<STATE>& data,int var, TPZVec<STATE>& sol) override;
    void Errors(const TPZMaterialDataT<STATE>& data,TPZVec<double>& values) override;

    // ---- serialização ----
    int  ClassId() const override;
    void Read (TPZStream& buf, void* context) override;
    void Write(TPZStream& buf, int withclassid) const override;

    // ---- kernel ----
    void SetKernelParams(REAL Lx, REAL Ly, STATE sigma2 = 1.) { fLx=Lx; fLy=Ly; fSigma2=sigma2; }
    void SetKernel(KernelType ker) { fKernel = std::move(ker); }
    static STATE ExpKernel(const TPZVec<REAL>& x, const TPZVec<REAL>& y,
                           REAL Lx, REAL Ly, STATE sigma2);

    // ---- Nyström explícito (fora do core) ----
    static void BuildB_Nystrom(TPZCompMesh& cmesh, int qorder, TPZFMatrix<STATE>& B);
    static void BuildC_Nystrom(TPZCompMesh& cmesh, int qorder,
                               const KernelType& ker, TPZFMatrix<STATE>& C);

private:
    struct GPEntry {
        TPZManVector<REAL,3> x;     // coord global
        double w = 0.;              // peso*|J|
        std::vector<int64_t> dof;   // índices globais
        std::vector<double>  phi;   // valores das shapes
    };

    static void BuildGlobalGPCatalog(TPZCompMesh& cmesh, int qorder, std::vector< TPZMatKLCov2D::GPEntry >& GP, int64_t& neq_out);

    static void ApplyPhiColumn(int64_t j, const std::vector<GPEntry>& GP,
                               std::vector<double>& t); // t = Φ e_j
    static void PhiT_times_vec(const std::vector<GPEntry>& GP,
                               const std::vector<double>& s,
                               TPZFMatrix<STATE>& y);   // y = Φᵀ s
    static void KernelApply(const std::vector<GPEntry>& GP,
                            const KernelType& ker,
                            const std::vector<double>& r,
                            std::vector<double>& s);   // s = K r

private:
    int   fDim     = 2;
    REAL  fLx      = 1.;
    REAL  fLy      = 1.;
    STATE fSigma2  = 1.;
    KernelType fKernel; // default = ExpKernel(Lx,Ly,σ²)
};
