#pragma once

#include "TPZMatBase.h"
#include "TPZMatCombinedSpaces.h"     // TPZMatCombinedSpacesT<STATE>
#include "TPZMatWithMem.h"
#include "TPZMaterialDataT.h"
#include "TPZBndCondT.h"
#include "TPZElasticMem.h"
//#include "TPZElasticResponse.h"
#include "TPZHash.h"
#include "pzreal.h"
#include "pzfmatrix.h"
#include "pzvec.h"
#include "pzfunction.h"
#include <string>
#include <iostream>
#ifdef PZ_LOG
// escolha um nome de categoria claro (use pontos para hierarquia)
static TPZLogger logger_poro_elastic("materials.PoroElastic2DMem");
#endif
enum class EPlaneType { PlaneStress=0, PlaneStrain=1 };

/// Poroelástico 2D (PRIMAL: u, p) com memória elástica (E, nu).
/// datavec[0] = u (vetor H1), datavec[1] = p (escalar H1).
template <class TMEM = TPZElasticMem>
class TPZMatPoroElastic2DMem
: public TPZMatBase<STATE, TPZMatCombinedSpacesT<STATE>, TPZMatWithMem<TMEM>>
{
    using TBase = TPZMatBase<STATE, TPZMatCombinedSpacesT<STATE>, TPZMatWithMem<TMEM>>;
public:
    explicit TPZMatPoroElastic2DMem(int matid,
                                    EPlaneType plane = EPlaneType::PlaneStress);

    // ---------- setup ----------
    void SetElasticResponse(const TPZElasticResponse &ER);
    void SetElasticity(STATE E, STATE nu){ TPZElasticResponse er; er.SetEngineeringData(E,nu); SetElasticResponse(er); }
    void SetAlpha(STATE a)                 { falpha = a; }
    void SetSe(STATE se)                   { fSe = se; }
    void SetPermeability(STATE k)          { fk = k; }
    void SetViscosity(STATE mu)            { fmu = mu; }
    void SetRhoF(STATE rhof)               { frhof = rhof; }
    void SetGravity(STATE gx, STATE gy)    { fG[0]=gx; fG[1]=gy; }
    //void SetGradP0(STATE dpdx, STATE dpdy) { fGradP0[0]=dpdx; fGradP0[1]=dpdy; }
    void SetBodyForce(STATE fx, STATE fy)  { fBody[0]=fx; fBody[1]=fy; }
    void SetTimeStep(STATE dt)             { fTimeStep = dt; }
    void UseMassScaledByInvDt(bool on)     { fMassInvDt = on; } // true: (1/dt)S; false: S + dt H
    void SetForcingFunctionP(TPZAutoPointer<TPZFunction<STATE>> f){ fForcingP = f; }

    // ---------- TPZMaterial básicos ----------
    std::string Name() const override { return "TPZMatPoroElastic2DMem"; }
    int Dimension() const override { return 2; }
    int NStateVariables() const override { return 1; } // irrelevante em multiphysics

// requisitos (mantêm sem const, pois você seta flags dentro)
void FillDataRequirements(TPZVec<TPZMaterialDataT<STATE>> &datavec) const override;
void FillBoundaryConditionDataRequirements(int, TPZVec<TPZMaterialDataT<STATE>> &datavec) const override;

// *** IMPORTANTE: datavec é const nestas três abaixo ***
void Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec,
                REAL weight,
                TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override;

void Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec,REAL weight, TPZFMatrix<STATE> &ef) override;

void ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec,
                  REAL weight,
                  TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef,
                  TPZBndCondT<STATE> &bc) override;

void ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec,REAL weight,TPZFMatrix<STATE> &ef,TPZBndCondT<STATE> &bc) override;


    // ---------- pós-processo ----------
    enum ESolutionVar { EDisplacement=1, EPressure=2, EFlux=3, EStrain=4, EStress=5, EYoung=6, EPoisson=7, EPOrder=8 ,EExactPressure=9,ExactPressureGradiendSolution=10,EGradP=11,EExactDisplacement=12};
    using TExact = std::function<void(const TPZVec<REAL>&, STATE&, TPZFMatrix<STATE>&)>;
    using TExactVec = std::function<void(const TPZVec<REAL>&, TPZVec<STATE>&, TPZFMatrix<STATE>&)>;
    void SetExact(TExact f) { fExact = std::move(f); }
    void SetExact(TExactVec f) { fExactVec = std::move(f); }
    int  VariableIndex(const std::string &name) const override;
    int  NSolutionVariables(int var) const override;
    void Solution(const TPZVec<TPZMaterialDataT<STATE>> &datavec,
                  int var, TPZVec<REAL> &Solout) override;//- persistência ----------
    int  ClassId() const override;
    void Write(TPZStream &buf, int withclassid) const override;
    void Read (TPZStream &buf, void *ctx) override;
    void Print(std::ostream & out = std::cout) const override;
    void SetHOnly(bool on)          { fHOnly = on; }
    bool HOnly() const              { return fHOnly; }

    void ComputePorePressure(const TPZMaterialDataT<STATE> & data, REAL & Pp, TPZVec<REAL> & dPp);
    void UpdatePorePressure(const TPZMaterialDataT<STATE> & data);
    void GetPrevPorePressure(const TPZMaterialDataT<STATE>& data,STATE &P0, TPZVec<STATE> &gradP0) const;
    void SetMem(const REAL Pp,const TPZElasticResponse ER);
    void UpdateMemory(const TPZVec< TPZMaterialDataT< STATE > >& datavec);
private:
    // helpers
    void ERFromMem(const TPZMaterialDataT<STATE> &data, STATE &E, STATE &nu) const;
    void BuildConstitutiveMatrix(STATE E, STATE nu, TPZFMatrix<STATE> &D) const; // 3x3
    void BuildBMatrix(const TPZFMatrix<STATE> &dphix, TPZFMatrix<STATE> &B) const; // 3 x 2n
    void BuildBu(const TPZFMatrix<STATE>& dphiU, TPZFMatrix<STATE>& Bu);
    void BuildBp(const TPZFMatrix<STATE>& dphiP, TPZFMatrix<STATE>& Bp);
    void BuildNpCol(const TPZFMatrix<STATE>& phiP, TPZFMatrix<STATE>& Np);


private:
    // plano + fallback
    EPlaneType fPlane = EPlaneType::PlaneStress;
    STATE fE_fallback = 1.0, fNu_fallback = 0.3;

    // poro
    STATE falpha = 1.0;
    STATE fSe    = 1.0;
    STATE fk     = 1.0;
    STATE fmu    = 1.0;
    STATE frhof  = 1.0;

    // ambiente
    STATE fG[2]      = {0.0, -9.81};
    STATE fBody[2]   = {0.0, 0.0};
    STATE fTimeStep  = 1.0;
    bool  fMassInvDt = true;

    bool  fHOnly = false;
    TPZAutoPointer<TPZFunction<STATE>> fForcingP = nullptr;

    // atalhos para TPZMatWithMem
    const TPZMatWithMem<TMEM>* WithMem() const { return dynamic_cast<const TPZMatWithMem<TMEM>*>(this); }
    TPZMatWithMem<TMEM>*       WithMem()       { return dynamic_cast<TPZMatWithMem<TMEM>*>(this); }
    TExact     fExact;
    TExactVec     fExactVec;
};




// instância explícita comum
extern template class TPZMatPoroElastic2DMem<TPZElasticMem>;
