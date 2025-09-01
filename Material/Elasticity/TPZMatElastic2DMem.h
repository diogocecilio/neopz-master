#pragma once

#include "TPZMatBase.h"
#include "TPZMatSingleSpace.h"
#include "TPZMatWithMem.h"
#include "TPZBndCondT.h"
#include "TPZElasticMem.h"
#include "TPZHash.h"
//#include "TPZFNMatrix.h"

/// Material elástico 2D com memória (E, nu) em cada ponto de integração.
/// TMEM deve expor um membro `TPZElasticResponse m_ER;` do qual se obtém E() e Poisson().
template <class TMEM = TPZElasticMem>
class TPZMatElastic2DMem
: public TPZMatBase<STATE,
                    TPZMatSingleSpaceT<STATE>,
                    TPZMatWithMem<TMEM>>
{
    //using TSpace   = TPZMatSingleSpaceT<STATE>;
    //using TWithMem = TPZMatWithMem<TMEM>;
    //using TBase    = TPZMatBase<STATE, TSpace, TWithMem>;

    using TBase = TPZMatBase<STATE,TPZMatSingleSpaceT<STATE>,TPZMatWithMem<TMEM>>;

public:
    enum class EPlaneType { PlaneStress, PlaneStrain };

    explicit TPZMatElastic2DMem(int matid,
                                EPlaneType plane = EPlaneType::PlaneStress);

    // ---------- configuração ----------
    void SetPlaneType(EPlaneType p){ fPlane = p; }
    EPlaneType PlaneType() const { return fPlane; }

    /// Valores fallback (usados apenas se o IP ainda não tem memória escrita)
    void SetElasticityFallback(REAL E, REAL nu){ fE_fallback = E; fNu_fallback = nu; }

    // ---------- overrides obrigatórios ----------
    std::string Name() const override { return "TPZMatElastic2DMem"; }
    int Dimension() const override { return 2; }
    int NStateVariables() const override { return 2; } // ux, uy

    void FillDataRequirements(TPZMaterialData &data) const override;
    void FillBoundaryConditionDataRequirements(int type, TPZMaterialData &data) const override;

    void Contribute(const TPZMaterialDataT<STATE> &data,
                    REAL weight,
                    TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef) override;

    void Contribute(const TPZMaterialDataT<STATE> &data,
                    REAL weight, TPZFMatrix<REAL> &ef) override;

    void ContributeBC(const TPZMaterialDataT<STATE> &data, REAL weight,
                      TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef,
                      TPZBndCondT<STATE> &bc) override;

    void ContributeBC(const TPZMaterialDataT<STATE> &data, REAL weight,
                      TPZFMatrix<REAL> &ef, TPZBndCondT<STATE> &bc) override;

    void SetElasticResponse(const TPZElasticResponse &ER);

    void SetElasticity(REAL E, REAL nu) {
        TPZElasticResponse er;
        er.SetEngineeringData(E, nu);   // ajuste o setter se sua branch tiver outro nome
        SetElasticResponse(er);

    }

     void Print(std::ostream & out = std::cout) const override;

    // ---------- pós-processo ----------
    enum ESolutionVar {
        EDisplacement = 1, // (ux,uy,0)
        EStrain,           // (exx, eyy, gxy) com gxy = 2*exy
        EStress,           // (sxx, syy, sxy)
        EYoung,            // E do IP
        EPoisson,          // nu do IP
        EPOrder
    };

    int  VariableIndex(const std::string &name) const override;
    int  NSolutionVariables(int var) const override;
    void Solution(const TPZMaterialDataT<STATE> &data, int var, TPZVec<REAL> &Solout) override;

    // ---------- persistência ----------
    int  ClassId() const override;
    void Write(TPZStream &buf, int withclassid) const override;
    void Read (TPZStream &buf, void *ctx) override;

private:
    // estado plano
    EPlaneType fPlane = EPlaneType::PlaneStress;

    // fallback
    REAL fE_fallback  = 1.0;
    REAL fNu_fallback = 0.0;

    // util: constitutiva 2D (Voigt xx,yy,xy) a partir de E,nu
    void BuildConstitutiveMatrix(REAL E, REAL nu, TPZFMatrix<REAL> &D) const;

    // util: monta B (3 x 2n) a partir de dphix (global) — clássico
    void BuildBMatrix(const TPZFMatrix<REAL> &dphix, TPZFMatrix<REAL> &B) const;

    // util: lê E, nu da memória (ou fallback)
    void ERFromMem(const TPZMaterialDataT<STATE> &data, REAL &E, REAL &nu) const;

    // atalhos para TPZMatWithMem
    const TPZMatWithMem<TMEM>* WithMem() const { return dynamic_cast<const TPZMatWithMem<TMEM>*>(this); }
    TPZMatWithMem<TMEM>*       WithMem()       { return dynamic_cast<TPZMatWithMem<TMEM>*>(this); }
};

// ---------------------- Instanciação explícita mais comum ----------------------
extern template class TPZMatElastic2DMem<TPZElasticMem>;
