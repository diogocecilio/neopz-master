// TPZPlasticStepVoigt.h
#ifndef TPZPLASTICSTEPVOIGT_H
#define TPZPLASTICSTEPVOIGT_H

#include "TPZPlasticBase.h"
#include "pzvec.h"
#include <memory>
#include "TPZYCVonMisesPV.h"
#include "TPZElasticResponse.h"
/// Classe constitutiva elasto-plástica em Voigt (3D) [11,22,33,13,23,12],
/// com critério de escoamento genérico YC e resposta elástica ER,
/// no mesmo espírito do TPZPlasticStepPV (mas operando em Voigt).
template <class YC, class ER>
class TPZPlasticStepVoigt : public TPZPlasticBase
{
public:
    using YieldCriterion = YC;
    using ElasticResponse = ER;

    // ===================== Construtores / Dtor =====================
    TPZPlasticStepVoigt(const YC& yc, const ER& er);
    TPZPlasticStepVoigt();                               // construtor padrão
    TPZPlasticStepVoigt(const TPZPlasticStepVoigt& other); // construtor de cópia
    ~TPZPlasticStepVoigt() override;                     // destrutor

    // ===================== Identificação / IO =====================
    int ClassId() const override;
    void Write(TPZStream& buf, int withclassid) const override;
    void Read(TPZStream& buf, void* context) override;
    void Print(std::ostream& out) const override;
    const char* Name() const override;

public:
    // --- exigidas pela TPZPlasticBase ---
    void ApplyStrain(const TPZTensor<REAL>& epsTotal) override;
    void ApplyStrainComputeSigma(const TPZTensor<REAL>& epsTotal,
                                 TPZTensor<REAL>& sigma,
                                 TPZFMatrix<REAL>* tangent = nullptr) override;
    void ApplyStrainComputeDep(const TPZTensor<REAL>& epsTotal,
                               TPZTensor<REAL>& sigma,
                               TPZFMatrix<REAL>& Dep) override;
    void ApplyLoad(const TPZTensor<REAL>& sigma, TPZTensor<REAL>& epsTotal) override;

    void SetState(const TPZPlasticState<REAL>& state) override;
    TPZPlasticState<REAL> GetState() const override;

    void Phi(const TPZTensor<REAL>& epsTotal, TPZVec<REAL>& phi) const override;

    void SetElasticResponse(TPZElasticResponse& ERin) override;
    TPZElasticResponse GetElasticResponse() const override;
    TPZPlasticCriterion& GetYC() override;


protected:
    ER   fER;                 // resposta elástica (p.ex., armazena K,G ou E,nu)
    YC   fYC;                 // critério de escoamento (deve derivar de TPZPlasticCriterion)
    TPZPlasticState<REAL> fN;
};

#endif // TPZPLASTICSTEPVOIGT_H
