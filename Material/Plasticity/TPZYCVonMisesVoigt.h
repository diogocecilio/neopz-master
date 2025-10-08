

#ifndef TPZYCVONMISESVoigt_H
#define TPZYCVONMISESVoigt_H

#include "pzlog.h"
#include "TPZTensor.h"
#include "pzvec_extras.h"
#include "TPZPlasticState.h"
#include "TPZElasticResponse.h"
#include "TPZPlasticCriterion.h"
#include "TPZHWTools.h"
#ifdef PZ_LOG
static TPZLogger loggerVonMIsesVoigt("pz.plasticity.vonmisespv");
#endif

class TPZYCVonMisesVoigt : public TPZPlasticCriterion {
private:

    REAL fSigmaY0;

    std::function<STATE(STATE)>  fSigmaY;

    std::function<STATE(STATE)> fH;


public:

    enum {
        NYield = 1
    };

    TPZYCVonMisesVoigt();

    TPZYCVonMisesVoigt(STATE sigmaY0, STATE Hiso );

    TPZYCVonMisesVoigt(const TPZYCVonMisesVoigt &cp);

    virtual void SetLocalMatState ( TPZPlasticState<REAL> & state )override;

    virtual TPZPlasticState<REAL> GetLocalMatState (  )override;

    virtual void ChangeLocalMatParameters( TPZPlasticState<REAL> & state ,REAL factor) override;

    /**
     * @brief Operator =
     */
    TPZYCVonMisesVoigt & operator=(const TPZYCVonMisesVoigt &cp);

    virtual int ClassId() const override;

    void Read(TPZStream& buf, void* context) override;

    void Write(TPZStream& buf, int withclassid) const override;

    /**
     * @brief Print Method
     */
    virtual void Print(std::ostream &out) const override;


    void ProjectSigma(const TPZTensor<STATE> & sigmatr, STATE k_prev, TPZTensor<STATE> & sigmaproj, STATE &k_proj, int & m_type);


    TPZTensor<STATE> ComputeN(const TPZTensor<STATE> stresstensor)const;

    TPZFMatrix<STATE> GetNdSigma(const TPZTensor<STATE>& sigma) const;

    STATE ComputeGamma(const TPZTensor<STATE>sig, const TPZFMatrix<STATE> elasticmat)const;

    void SetHardening(std::function<STATE(STATE)>  H)      { fH = std::move(H); }
    void SetYieldStressLaw(std::function<STATE(STATE)>  sy)   { fSigmaY = std::move(sy); }

    void SetUp(STATE sigmaY0, STATE Hiso);

    // acesso seguro
    STATE H(STATE kappa)      const { return fH ? fH(kappa) : STATE(0); }
    STATE SigmaY(STATE kappa) const { return fSigmaY ? fSigmaY(kappa) : fSigmaY0; }

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
