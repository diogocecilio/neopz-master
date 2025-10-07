

#ifndef TPZYCVONMISESPV_H
#define TPZYCVONMISESPV_H

#include "pzlog.h"
#include "TPZTensor.h"
#include "pzvec_extras.h"
#include "TPZPlasticState.h"
#include "TPZElasticResponse.h"
#include "TPZPlasticCriterion.h"
#include "TPZHWTools.h"
#ifdef PZ_LOG
static TPZLogger loggerVonMIsesPV("pz.plasticity.vonmisespv");
#endif

class TPZYCVonMisesPV : public TPZPlasticCriterion {
private:

    REAL fYieldStress;


public:

    enum {
        NYield = 1
    };

    TPZYCVonMisesPV();

    TPZYCVonMisesPV(REAL yieldstress, TPZElasticResponse &ER);

    TPZYCVonMisesPV(const TPZYCVonMisesPV &cp);

    void SetUp(REAL yieldstress, TPZElasticResponse &ER);

    virtual void SetLocalMatState ( TPZPlasticState<REAL> & state )override;

    virtual TPZPlasticState<REAL> GetLocalMatState (  )override;

    virtual void ChangeLocalMatParameters( TPZPlasticState<REAL> & state ,REAL factor) override;

    /**
     * @brief Operator =
     */
    TPZYCVonMisesPV & operator=(const TPZYCVonMisesPV &cp);

    virtual int ClassId() const override;

    void Read(TPZStream& buf, void* context) override;

    void Write(TPZStream& buf, int withclassid) const override;


    /**
     * @brief Print Method
     */
    virtual void Print(std::ostream &out) const override;

    void SetElasticResponse(const TPZElasticResponse &ER) {
        //fER = ER;
    }

    virtual TPZElasticResponse GetElasticResponse() const {
        //return fER;
    }
    /**
     * @brief sigma = lambda Tr(E)I + 2 mu E
     */
    template<class T>
    TPZVec<T> SigmaElastPV(const TPZVec<T> &deform) const;



    void ProjectSigma(const TPZTensor<STATE> & epst,const TPZTensor<STATE> & epsp,STATE k_prev);


    /**
     Evaluates the yield criterion

     @param sig_vec principal stress
     @param alpha internal damage variable
     @param phi yield criterion function
     */
    void Phi(TPZVec<STATE> sig_vec, STATE alpha, TPZVec<STATE> &phi)const;


    virtual void YieldFunction(const TPZVec<STATE>& sigma, STATE kprev, TPZVec<STATE>& yield) const override{
        Phi(sigma, kprev, yield);
    }

    virtual int GetNYield() const override{
        return as_integer(NYield);
    }


};


#endif //TPZYCVonMisesPV
