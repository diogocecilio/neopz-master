

#ifndef TPZYCTRESCAVoigt_H
#define TPZYCTRESCAVoigt_H

#include "pzlog.h"
#include "TPZTensor.h"
#include "pzvec_extras.h"
#include "TPZPlasticState.h"
#include "TPZElasticResponse.h"
#include "TPZPlasticCriterion.h"
#include "TPZHWTools.h"


class TPZYCTrescaVoigt : public TPZPlasticCriterion {
private:
    STATE fSigmaY0 = 0.0;
    STATE fH0 = 0.0;
    STATE fHard=0.;
    TPZElasticResponse fER;

public:

    enum {
        NYield = 1
    };

    TPZYCTrescaVoigt();

    TPZYCTrescaVoigt(const TPZYCTrescaVoigt &cp);

    virtual void SetLocalMatState ( TPZPlasticState<REAL> & state )override;

    virtual TPZPlasticState<REAL> GetLocalMatState (  )override;

    virtual void ChangeLocalMatParameters( TPZPlasticState<REAL> & state ,REAL factor) override;


    virtual int ClassId() const override;

    void Read(TPZStream& buf, void* context) override;

    void Write(TPZStream& buf, int withclassid) const override;

    /**
     * @brief Print Method
     */
    virtual void Print(std::ostream &out) const override;


    STATE ProjectSigma(const TPZTensor<STATE> & sigmatr,  TPZTensor<STATE> & sigmaproj, TPZElasticResponse &ER,STATE &havarn,STATE &havarn1, int & m_type);

    STATE ProjectSigmaDep(const TPZTensor<STATE> & sigmatr,  TPZTensor<STATE> & sigmaproj, TPZElasticResponse &ER,STATE &havarn,STATE &havarn1, int & m_type);

    STATE ProjectSigma(TPZManVector<STATE,3> &sigtr,STATE &alphan,TPZManVector<STATE,2> &dlambda,TPZManVector<STATE,3> &sigpr,TPZManVector<STATE,3> &epstr,TPZFNMatrix<9> &Grad3x3,STATE &alphan1,int & m_type);


    bool ComputeLambdaSigmaMainPlane(TPZManVector<STATE,3> &sigtr,STATE &alphan,TPZManVector<STATE,2> &dlambda,TPZManVector<STATE,3> &sigpr,TPZManVector<STATE,3> &epstr,TPZFNMatrix<9> &Grad3x3,STATE &alphan1);

    bool ComputeLambdaSigmaLeft(TPZManVector<STATE,3> &sigtr,STATE &alphan,TPZManVector<STATE,2> &dlambda,TPZManVector<STATE,3> &sigpr,TPZManVector<STATE,3> &epstr,TPZFNMatrix<9> &Grad3x3,STATE &alphan1);

    bool ComputeLambdaSigmaRigth(TPZManVector<STATE,3> &sigtr,STATE &alphan,TPZManVector<STATE,2> &dlambda,TPZManVector<STATE,3> &sigpr,TPZManVector<STATE,3> &epstr,TPZFNMatrix<9> &Grad3x3,STATE &alphan1);
    // STATE ComputeGamma(const TPZTensor<STATE>sig, const TPZFMatrix<STATE> elasticmat)const;

    void SetHardening(STATE H)      { fH0 = H; }

    void SetYieldStress(STATE sy)   { fSigmaY0 = sy; }


    void SetUp(STATE sigmaY0, STATE Hiso,TPZElasticResponse &ER);


    // acesso seguro
    STATE H()      const { return fH0 ;}

   //
     STATE SigmaY(STATE hardeningvar) const { return  fSigmaY0+hardeningvar*fH0;  }
     STATE DSigmaYDepsbar(STATE hardeningvar)  const { return  fH0;  }

    // STATE SigmaY(STATE hardeningvar) const { return   (278.51775588600316 + 107.03078301291825*(1 - pow(exp(1),-450.3386387920765*hardeningvar)) + 1122.6997583510627*hardeningvar - 1996.235137824671*pow(hardeningvar,2)) ;}
    // STATE DSigmaYDepsbar(STATE hardeningvar) const {return (1122.6997583510627 + 48200.097130887705/pow(2.718281828459045,450.3386387920765*hardeningvar) - 3992.470275649342*hardeningvar);}

    /**
     Evaluates the yield criterion

     @param sig_vec principal stress
     @param alpha internal damage variable
     @param phi yield criterion function
     */
    void Phi(TPZTensor<STATE> sig, STATE alpha, TPZVec<STATE> &phi)const;


    virtual void YieldFunction(const TPZVec<STATE>& sigma, STATE kprev, TPZVec<STATE>& yield) const override{
        DebugStop();
    }

    virtual int GetNYield() const override{
        return as_integer(NYield);
    }

};


#endif //TPZYCVonMisesPV
