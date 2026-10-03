// TPZMatPoroElastoPlastic3DMem.h
//
// Material multifísico u-p (Biot) 3D ou em deformação plana (dim = 2) com lei elastoplástica para a tensão
// efetiva e memória nos pontos de integração (versão com Newton e tangente consistente do
// TPZMatPoroElastoPlastic2DMem).
//
//   datavec[0] = u (H1 vetorial, dim componentes), datavec[1] = p (H1 escalar)
//
// Em deformação plana ε_zz = ε_xz = ε_yz = 0 e a lei constitutiva (3D, Voigt) dá σ'_zz.
//
// Formulação total (a solução guarda u e p totais; a memória guarda o estado convergido do passo n),
// Euler implícito, com o resíduo e a jacobiana do método de Newton:
//
//   R_u = ∫ Bᵀσ'(ε) - α ∫ Bᵀm N_p p - ∫ N_uᵀ b                        (- forças de contorno)
//   R_p = ∫ N_pᵀ [α (ε_v - ε_v,n) + S (p - p_n)] + Δt ∫ ∇N_pᵀ (K/μ) (∇p - ρ_f g⃗)
//   J   = [[∫ Bᵀ D_ep B, -α ∫ Bᵀm N_pᵀ], [α ∫ N_p mᵀB, S ∫ N_p N_pᵀ + Δt ∫ ∇N_pᵀ (K/μ) ∇N_p]]
//
// K = diag(k_x, k_y, k_z) (SetPermeability com um ou três valores).
//
// ek = J e ef = -R (o NeoPZ resolve ek Δx = ef). σ' vem de T::ApplyStrainComputeSigma com o estado da
// memória (ε_n, ε^p_n, α_n) e σ'_n na entrada (necessária nas leis hipoelásticas). O fluxo (H, f_g) só
// entra com SetFlow(true) (etapas drenadas / adensamento); com SetFlow(false) o passo é não drenado.
// Convenções: tração positiva, p > 0 em compressão, σ = σ' - α p I, Voigt do TPZTensor com distorções
// de engenharia. O índice de memória é o do elemento atômico de u (datavec[0].intGlobPtIndex).
//
// Condições de contorno (Val2 com 4 componentes {v_x, v_y, v_z, p}, também em 2D; com SetForcingFunctionBC
// os valores vêm da função, avaliada em cada ponto, com rhsVal de 4 componentes — p.ex. o nível d'água de um
// rebaixamento que varia no tempo):
//   0  u = v (todas as componentes)                 1  tração t = v
//   2  p = Val2[3]                                  3  u_i = 0 nas direções com v_i != 0
//   5  pressão normal: t = Val2[0] n                 6  u_i = v_i nas direções com Val1(i,i) != 0
//   12 tração t = v e p = Val2[3]                   16 tipo 6 e p = Val2[3]
// As condições de Dirichlet são impostas por penalidade (BigNumber()).
//
#ifndef TPZMATPOROELASTOPLASTIC3DMEM_H
#define TPZMATPOROELASTOPLASTIC3DMEM_H

#include <functional>
#include <string>

#include "TPZMatBase.h"
#include "TPZMatCombinedSpaces.h"
#include "TPZMatWithMem.h"
#include "TPZMaterialDataT.h"
#include "TPZBndCondT.h"
#include "TPZElastoPlasticMem.h"
#include "TPZTensor.h"
#include "pzfmatrix.h"

template <class T, class TMEM = TPZElastoPlasticMem>
class TPZMatPoroElastoPlastic3DMem
    : public TPZMatBase<STATE, TPZMatCombinedSpacesT<STATE>, TPZMatWithMem<TMEM>> {
    using TBase = TPZMatBase<STATE, TPZMatCombinedSpacesT<STATE>, TPZMatWithMem<TMEM>>;

public:
    /// O que Contribute monta
    enum EAssembleMode {
        EFull = 0,            ///< resíduo e jacobiana completos
        EExternalForces = 1,  ///< só as forças externas (forças de corpo e cargas de contorno) em ef
        ENoPenalty = 2        ///< resíduo completo sem as penalidades de Dirichlet (reações nos gdl restritos)
    };

    enum EBCType {
        EDirichletU = 0, ETraction = 1, EDirichletP = 2, EDirectionalNullU = 3, ENormalPressure = 5,
        EDirectionalU = 6, ETractionDirichletP = 12, EDirectionalUDirichletP = 16
    };

    enum ESolutionVar {
        EDisplacement = 1, EPressure = 2, EExcessPressure = 3, EFlux = 4, EPressureGradient = 5
    };

    TPZMatPoroElastoPlastic3DMem();
    /// dim = 3 (sólido) ou 2 (deformação plana)
    explicit TPZMatPoroElastoPlastic3DMem(int matid, int dim = 3);

    // ---------------------------------------------------------------- dados
    void SetPlasticModel(const T &model) { fPlasticModel = model; }
    T &GetPlasticModel() { return fPlasticModel; }
    /// Ajusta o modelo em cada ponto (p.ex. tensão inicial e parâmetros que variam com a profundidade)
    void SetModelUpdate(const std::function<void(const TPZVec<REAL> &x, T &model)> &f) { fModelUpdate = f; }
    void SetAlpha(STATE a) { fAlpha = a; }                       ///< coeficiente de Biot
    void SetSe(STATE se) { fSe = se; }                           ///< armazenamento 1/M
    void SetPermeability(STATE k) { fK[0] = fK[1] = fK[2] = k; }
    /// Permeabilidade ortótropa nos eixos globais
    void SetPermeability(STATE kx, STATE ky, STATE kz) { fK[0] = kx; fK[1] = ky; fK[2] = kz; }
    void SetViscosity(STATE mu) { fMu = mu; }
    void SetRhoF(STATE rhof) { fRhoF = rhof; }
    void SetGravity(STATE gx, STATE gy, STATE gz) { fG[0] = gx; fG[1] = gy; fG[2] = gz; }
    void SetBodyForce(STATE fx, STATE fy, STATE fz) { fBody[0] = fx; fBody[1] = fy; fBody[2] = fz; }
    void SetTimeStep(STATE dt) { fTimeStep = dt; }
    void SetFlow(bool flow) { fFlow = flow; }
    void SetAssembleMode(EAssembleMode m) { fMode = m; }
    /// Ordem da regra de integração (<= 0: padrão, 2 p_max)
    void SetIntegrationOrder(int order) { fIntegrationOrder = order; }
    /// Poropressão hidrostática (para a variável ExcessPressure)
    void SetHydrostatic(const std::function<REAL(const TPZVec<REAL> &x)> &f) { fHydro = f; }
    static constexpr REAL BigNumber() { return 1.e16; }

    // ---------------------------------------------------------------- TPZMaterial
    std::string Name() const override { return "TPZMatPoroElastoPlastic3DMem"; }
    int Dimension() const override { return fDim; }
    int NStateVariables() const override { return fDim + 1; }
    TPZMaterial *NewMaterial() const override { return new TPZMatPoroElastoPlastic3DMem<T, TMEM>(*this); }
    int IntegrationRuleOrder(const TPZVec<int> &elPMaxOrder) const override;

    void FillDataRequirements(TPZVec<TPZMaterialDataT<STATE>> &datavec) const override;
    void FillBoundaryConditionDataRequirements(int type, TPZVec<TPZMaterialDataT<STATE>> &datavec) const override;

    void Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek,
                    TPZFMatrix<STATE> &ef) override;
    void Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ef) override;
    void ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek,
                      TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc) override;
    void ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ef,
                      TPZBndCondT<STATE> &bc) override;

    int VariableIndex(const std::string &name) const override;
    int NSolutionVariables(int var) const override;
    void Solution(const TPZVec<TPZMaterialDataT<STATE>> &datavec, int var, TPZVec<STATE> &Solout) override;

    int ClassId() const override;
    void Write(TPZStream &buf, int withclassid) const override;
    void Read(TPZStream &buf, void *context) override;
    void Print(std::ostream &out = std::cout) const override;

private:
    void ContributeInternal(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> *ek,
                            TPZFMatrix<STATE> &ef);
    void ContributeBCInternal(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> *ek,
                              TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc);

    T fPlasticModel;
    std::function<void(const TPZVec<REAL> &, T &)> fModelUpdate;
    std::function<REAL(const TPZVec<REAL> &)> fHydro;
    int fDim = 3;
    STATE fAlpha = 1., fSe = 0., fMu = 1., fRhoF = 0.;
    STATE fK[3] = {0., 0., 0.};
    STATE fG[3] = {0., 0., 0.};
    STATE fBody[3] = {0., 0., 0.};
    STATE fTimeStep = 1.;
    bool fFlow = false;
    EAssembleMode fMode = EFull;
    int fIntegrationOrder = -1;
};

#endif
