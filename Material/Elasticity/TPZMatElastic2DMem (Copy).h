#pragma once
#include <ostream>
#include "pzvec.h"
#include "pzmatrix.h"
#include "pzmanvector.h"
#include "TPZMatWithMem.h"
#include "TPZMatBase.h"
#include "TPZMatSingleSpace.h"
#include "TPZBndCondT.h"

// Memória mínima por ponto de integração
struct TKLPointMem {
    STATE fE = 0.;
    // persistência
    int ClassId() const { return Hash("TKLPointMem"); }
    void Write(TPZStream& buf, int withclassid) const { buf.Write(&fE,1); }
    void Read (TPZStream& buf, void*)              { buf.Read(&fE,1); }
    // debug
    void Print(std::ostream& out) const { out << "E=" << fE; }
};
// para o TPZMatWithMem::Print
inline std::ostream& operator<<(std::ostream& os, const TKLPointMem& m){
    m.Print(os); return os;
}

/// Material elástico 2D isotrópico com memória: E(x_qp) vem de mem.fE
template<class TMEM = TKLPointMem>
class TPZMatElastic2DMem :
    public TPZMatBase<STATE,
                      TPZMatSingleSpaceT<STATE>,
                      TPZMatWithMem<TMEM>>
{
public:
    using TBase = TPZMatBase<STATE, TPZMatSingleSpaceT<STATE>, TPZMatWithMem<TMEM>>;

    using TBase::FillDataRequirements;
    using TBase::FillBoundaryConditionDataRequirements;

    // as suas versões (assinatura igual à da base) – agora realmente override
    //void FillDataRequirements(TPZMaterialDataT<STATE>& data) const override;
    //void FillBoundaryConditionDataRequirements(int type,TPZMaterialDataT<STATE>& data) const override;

    TPZMatElastic2DMem(int matid, bool planeStress, STATE nu, STATE Eref=(STATE)1.)
    : TBase(matid), fNu(nu), fEref(Eref), fPlaneStress(planeStress) {}

    // parâmetros
    void SetPoisson(STATE nu){ fNu = nu; }
    void SetPlaneStress(bool ps){ fPlaneStress = ps; }
    void SetEref(STATE Eref){ fEref = Eref; }
    //void SetBodyForce(STATE fx, STATE fy) { fBody[0]=fx; fBody[1]=fy; }

    // básicos
    int NStateVariables() const override { return 2; }
    int Dimension()        const override { return 2; }

    // montagem
    void Contribute(const TPZMaterialDataT<STATE>& data,
                    REAL weight, TPZFMatrix<STATE>& ek,
                    TPZFMatrix<STATE>& ef) override;

    void ContributeBC(const TPZMaterialDataT<STATE>& data,
                      REAL weight, TPZFMatrix<STATE>& ek,
                      TPZFMatrix<STATE>& ef,
                      TPZBndCondT<STATE>& bc) override;

    // requerimentos (definidos no .cpp para evitar tipo incompleto)
    void FillDataRequirements(TPZMaterialDataT<STATE>& data) const;
    void FillBoundaryConditionDataRequirements(int, TPZMaterialDataT<STATE>& data) const;

    // pós-processo
    enum EVarIds { EVar_Displacement=1, EVar_Strain, EVar_Stress, EVar_E };
    int VariableIndex(const std::string& name) const override {
        if(name=="Displacement") return EVar_Displacement;
        if(name=="Strain")       return EVar_Strain;
        if(name=="Stress")       return EVar_Stress;
        if(name=="E")            return EVar_E;
        return -1;
    }
    int NSolutionVariables(int var) const override {
        switch(var){
            case EVar_Displacement: return 2;
            case EVar_Strain:       return 3;
            case EVar_Stress:       return 3;
            case EVar_E:            return 1;
        }
        return 0;
    }
    void Solution(const TPZMaterialDataT<STATE>& data,
                  int var, TPZVec<STATE>& solout) override;

    // persistência
    int ClassId() const override { return Hash("TPZMatElastic2DMem") ^ TBase::ClassId() << 1; }
    void Write(TPZStream& buf, int withclassid) const override {
        TBase::Write(buf, withclassid);
        buf.Write(&fNu,1);
        buf.Write(&fEref,1);
        buf.Write((bool)fPlaneStress);
    }
    void Read(TPZStream& buf, void* ctx) override {
        TBase::Read(buf, ctx);
        buf.Read(&fNu,1);
        buf.Read(&fEref,1);
        bool ps=true; buf.Read(ps); fPlaneStress = ps;
    }

private:
    // matriz constitutiva 3x3
    void ElasticD(STATE E, STATE nu, TPZFNMatrix<9,STATE>& D) const;
    // pega E do mem (ou Eref)
    STATE EfromMem(const TPZMaterialDataT<STATE>& data) ;

private:
    STATE fNu = 0.3;
    STATE fEref = (STATE)1.;
    bool  fPlaneStress = true;
};
