/*
 *  FEMPZ
 *
 *  Created by Nathan Shauer on 5/4/13.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef TPZYCMOHRCOULOMBPV2_H
#define TPZYCMOHRCOULOMBPV2_H

#include "pzlog.h"
#include "TPZTensor.h"
#include "pzvec_extras.h"
#include "TPZPlasticState.h"
#include "TPZElasticResponse.h"
#include "TPZPlasticCriterion.h"

#ifdef PZ_LOG
static TPZLogger loggerMohrCoulombPV2("pz.plasticity.mohrcoulombpv");
#endif

class TPZYCMohrCoulombPV2 : public TPZPlasticCriterion {
public:

    enum {
        NYield = 3
    };

private:

    REAL fPhi;
    REAL fPsi;
    REAL fc;
    TPZElasticResponse fER;

protected:
    REAL fEpsPlasticBar;

public:



    TPZYCMohrCoulombPV2();

    TPZYCMohrCoulombPV2(REAL Phi, REAL Psi, REAL c, TPZElasticResponse &ER);

    TPZYCMohrCoulombPV2(const TPZYCMohrCoulombPV2 &cp);

    void SetUp(REAL Phi, REAL Psi, REAL c, TPZElasticResponse &ER) {
        fPhi = Phi;
        fPsi = Psi;
        fc = c;
        fER = ER;
    }
    TPZYCMohrCoulombPV2 & operator=(const TPZYCMohrCoulombPV2 &cp);

    virtual int ClassId() const override;

    void Read(TPZStream& buf, void* context) override;

    void Write(TPZStream& buf, int withclassid) const override;
    virtual void YieldFunction(const TPZVec<STATE> &sigma, STATE kprev, TPZVec<STATE> &yield) const override;

    virtual int GetNYield() const override;

    virtual void SetLocalMatState ( TPZPlasticState<REAL> & state )override;

    virtual TPZPlasticState<REAL> GetLocalMatState (  )override;

    virtual void ChangeLocalMatParameters( TPZPlasticState<REAL> & state ,REAL factor)override;


    virtual void Print(std::ostream &out) const {
        std::cout << __PRETTY_FUNCTION__ << " Should not be called, please check children classes." << std::endl;
        //        DebugStop();
    }

    TPZTensor<STATE> ComputeN(const TPZTensor<STATE> stresstensor)const;

    TPZFMatrix<STATE> GetNdSigma(const TPZTensor<STATE>& sigma) const;

    STATE ProjectSigma(const TPZTensor<STATE> & sigmatr,  TPZTensor<STATE> & sigmaproj, TPZElasticResponse &ER,STATE &havarn,STATE &havarn1, int & m_type);

    STATE ProjectSigma(TPZManVector<STATE,3> &sigtr,STATE &alphan,TPZManVector<STATE,2> &dlambda,TPZManVector<STATE,3> &sigpr,TPZManVector<STATE,3> &epstr,TPZFNMatrix<9> &Grad3x3,STATE &alphan1,int & m_type);

    STATE DSigmaYDepsbar(STATE hardeningvar)  const { return  0;  }



    bool ComputeLambdaSigmaMainPlane(TPZManVector<STATE,3> &sigtr,STATE &alphan,TPZManVector<STATE,2> &dlambda,TPZManVector<STATE,3> &sigpr,TPZManVector<STATE,3> &epstr,TPZFNMatrix<9> &Grad3x3,STATE &alphan1);

    bool ComputeLambdaSigmaLeft(TPZManVector<STATE,3> &sigtr,STATE &alphan,TPZManVector<STATE,2> &dlambda,TPZManVector<STATE,3> &sigpr,TPZManVector<STATE,3> &epstr,TPZFNMatrix<9> &Grad3x3,STATE &alphan1);

    bool ComputeLambdaSigmaRigth(TPZManVector<STATE,3> &sigtr,STATE &alphan,TPZManVector<STATE,2> &dlambda,TPZManVector<STATE,3> &sigpr,TPZManVector<STATE,3> &epstr,TPZFNMatrix<9> &Grad3x3,STATE &alphan1);

    bool ReturnMapApex(TPZManVector<STATE,3> &sigtr,STATE &alphan,TPZManVector<STATE,2> &dlambda,TPZManVector<STATE,3> &sigpr,TPZManVector<STATE,3> &epstr,TPZFNMatrix<9> &Grad3x3,STATE &alphan1);


};


#endif //TPZYCMOHRCOULOMBPV_H
