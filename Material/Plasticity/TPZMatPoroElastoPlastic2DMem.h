#pragma once

#include "TPZMatBase.h"
#include "TPZMatCombinedSpaces.h"     // TPZMatCombinedSpacesT<STATE>
#include "TPZMatWithMem.h"
#include "TPZMaterialDataT.h"
#include "TPZBndCondT.h"

#include "TPZHash.h"
#include "pzreal.h"
#include "pzfmatrix.h"
#include "pzvec.h"
#include "pzfunction.h"
#include <string>
#include <iostream>
#include "TPZMatElastoPlastic.h"

#include "TPZMatSingleSpace.h"
#include "TPZMatWithMem.h"
#include "TPZMatErrorSingleSpace.h"

#include "pzporoelastoplasticmem.h"
#include "TPZElasticResponse.h"
#ifdef PZ_LOG
// escolha um nome de categoria claro (use pontos para hierarquia)
static TPZLogger logger_poro("materials.PoroPlastic2DMem");
#endif

#include "TPZPlasticStepPV.h"
#include "TPZYCMohrCoulombPV.h"
#include "TPZMatElastoPlastic2D.h"
#include "TPZPorousElastoPlasticMem.h"
#include "TPZTensor.h"
//Implementa material poroelastoplastico multifisico
template <class T, class TMEM = TPZPorousElastoPlasticMem>
class  TPZMatPoroElastoPlastic2DMem :
public TPZMatBase<STATE,TPZMatCombinedSpacesT<STATE>,TPZMatWithMem<TMEM>,TPZMatErrorSingleSpace<STATE>>
{
	using TBase = TPZMatBase<STATE,TPZMatCombinedSpacesT<STATE>,TPZMatWithMem<TMEM>,TPZMatErrorSingleSpace<STATE>>;

public:

	enum ESolutionVar { EDisplacement=1,EPressure=2,ElasticStrain=3,EPlasticStrain=4,EStrainPlasticJ2=5,EStrainPlasticI1=6,ECoesion=7,EAtrito=8};

	enum EWhichMatrix { EK=1,EQ=2,EQT=3,ES=4,EH=5};

	TPZMatPoroElastoPlastic2DMem();

	TPZMatPoroElastoPlastic2DMem(int id);

	TPZMatPoroElastoPlastic2DMem(const TPZMatPoroElastoPlastic2DMem<T,TMEM> &cp);

	virtual ~TPZMatPoroElastoPlastic2DMem();

	virtual void Print(std::ostream & out = std::cout) const override;

	virtual int  VariableIndex(const std::string &name) const override;

	virtual int  NSolutionVariables(int var) const override;

	virtual void Solution(const TPZVec<TPZMaterialDataT<STATE>> &datavec,int var, TPZVec<REAL> &Sol) override;

	void FillDataRequirements(TPZVec<TPZMaterialDataT<STATE>> &datavec) const override;

	void FillBoundaryConditionDataRequirements(int, TPZVec<TPZMaterialDataT<STATE>> &datavec) const override;

	virtual void Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec,REAL weight,TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override;

	virtual void Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec,REAL weight, TPZFMatrix<STATE> &ef) override;

	virtual void ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec,REAL weight,TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc) override;

	void ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec,REAL weight,TPZFMatrix<STATE> &ef,TPZBndCondT<STATE> &bc) override;

	void Errors(const TPZMaterialDataT<STATE>&data,
				TPZVec<REAL> &values) override;

	virtual void ComputeDeltaStrainVector(const TPZMaterialDataT<STATE> & data, TPZFMatrix<REAL> &DeltaStrain);

	virtual void ApplyDeltaStrainComputeDep(const TPZMaterialDataT<STATE> & data, TPZFMatrix<REAL> & DeltaStrain,
											TPZFMatrix<REAL> & Stress, TPZFMatrix<REAL> & Dep);

	virtual void ApplyDeltaStrain(const TPZMaterialDataT<STATE> & data, TPZFMatrix<REAL> & DeltaStrain,
								  TPZFMatrix<REAL> & Stress);

	void BuildBu(const TPZFMatrix<STATE>& dphiU, TPZFMatrix<STATE>& Bu);

	void BuildNu(const TPZFMatrix<STATE>& phiU, TPZFMatrix<STATE>& Nu);

	void BuildBp(const TPZFMatrix<STATE>& dphiP, TPZFMatrix<STATE>& Bp);


	void SetPlasticModel(T & plasticmodel);


	void SetBulkDensity(REAL & RhoB);

	void SetPorousElasticity(TPZElasticResponse & PER);

	TPZElasticResponse GetPorousElasticity(TPZElasticResponse & PER);

	virtual T & GetPlasticModel();

	virtual TPZMaterial * NewMaterial() const override;

	virtual int ClassId() const override;

	virtual void Write(TPZStream &buf, int withclassid) const override;

	virtual void Read(TPZStream &buf, void *context) override;

	virtual std::string Name() const override;

	virtual int Dimension() const override { return 2; }

	virtual int NStateVariables() const override { return 2; }

	void SetElasticResponse(const TPZElasticResponse &ER);
	void SetElasticity(STATE E, STATE nu){ TPZElasticResponse er; er.SetEngineeringData(E,nu); SetElasticResponse(er); }
	void SetAlpha(STATE a)                 { falpha = a; }
	void SetSe(STATE se)                   { fSe = se; }
	void SetPermeability(STATE k)          { fk = k; }
	void SetViscosity(STATE mu)            { fmu = mu; }
	void SetRhoF(STATE rhof)               { frhof = rhof; }
	void SetBodyForce(STATE fx, STATE fy)  { fBody[0]=fx; fBody[1]=fy; }
	void SetTimeStep(STATE dt)             { fTimeStep = dt; }
	// --- DEBUG helpers (mínimos) -----------------------------------------------

	void SetWhichAssemble(EWhichMatrix Whichassemble)
	{
		fWhichassemble=Whichassemble;
	}

	EWhichMatrix GetWhichAssemble()
	{
		return fWhichassemble;
	}


protected:

	EWhichMatrix fWhichassemble=EK;

    int fPlaneStrain;

	REAL fFactor;

	T fPlasticityModel;

	TPZElasticResponse m_PER;

	const TPZMatWithMem<TMEM>* WithMem() const { return dynamic_cast<const TPZMatWithMem<TMEM>*>(this); }
	TPZMatWithMem<TMEM>*       WithMem()       { return dynamic_cast<TPZMatWithMem<TMEM>*>(this); }

	STATE falpha = 1.0;
	STATE fSe    = 1.0;
	STATE fk     = 1.0;
	STATE fmu    = 1.0;
	STATE frhof  = 1.0;

	// ambiente
	STATE fBody[2]   = {0.0, 0.0};
	STATE fTimeStep  = 1.0;

	bool fConsistentMatrix=false;
};




// instância explícita comum
//extern template class TPZMatPoroElastoPlastic2DMem<TPZPlasticStepPV<TPZYCMohrCoulombPV,TPZElasticResponse>,TPZPorousElastoPlasticMem>;
