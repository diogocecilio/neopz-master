

#ifndef TPZYCVONMISESVoigt_H
#define TPZYCVONMISESVoigt_H

#include "pzlog.h"
#include "TPZTensor.h"
#include "pzvec_extras.h"
#include "TPZPlasticState.h"
#include "TPZElasticResponse.h"
#include "TPZPlasticCriterion.h"
#include "TPZHWTools.h"
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#ifdef PZ_LOG
static TPZLogger loggerVonMIsesVoigt("pz.plasticity.vonmisespv");
#endif

class TPZYCVonMisesVoigt : public TPZPlasticCriterion {
private:
    STATE fSigmaY0 = 0.0;
    STATE fH0 = 0.0;


public:

    enum {
        NYield = 1
    };

    TPZYCVonMisesVoigt();

    TPZYCVonMisesVoigt(const TPZYCVonMisesVoigt &cp);

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


    void ProjectSigma(const TPZTensor<STATE> & sigmatr,  TPZTensor<STATE> & sigmaproj, STATE &cumhardenig, int & m_type);


    TPZTensor<STATE> ComputeN(const TPZTensor<STATE> stresstensor)const;

    TPZFMatrix<STATE> GetNdSigma(const TPZTensor<STATE>& sigma) const;

    // STATE ComputeGamma(const TPZTensor<STATE>sig, const TPZFMatrix<STATE> elasticmat)const;

    void SetHardening(STATE H)      { fH0 = H; }

    void SetYieldStress(STATE sy)   { fSigmaY0 = sy; }

    STATE UpdateHardeningVar(const TPZTensor<STATE>sig,const TPZElasticResponse& ER, STATE &hardeningvar);

    void SetUp(STATE sigmaY0, STATE Hiso);

    // acesso seguro
    STATE H()      const { return fH0 ;}

    STATE SigmaY(STATE hardeningvar) const { return   278.51775588600316 + 107.03078301291825*(1 - pow(exp(1),-450.3386387920765*hardeningvar)) + 1122.6997583510627*hardeningvar - 1996.235137824671*pow(hardeningvar,2) ;}
    //STATE SigmaY(STATE hardeningvar) const { return  fSigmaY0+hardeningvar*fH0;  }

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

    inline void Decomp3x3(const TPZTensor<STATE>& sig,
                          TPZTensor<STATE>::TPZDecomposed &out)
    {
        Eigen::Matrix3d A;
        A << sig.XX(), sig.XY(), sig.XZ(),
        sig.XY(), sig.YY(), sig.YZ(),
        sig.XZ(), sig.YZ(), sig.ZZ();
        Eigen::ComplexEigenSolver<Eigen::MatrixXd> ces;

        ces.compute ( A );
        //Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(A);
        // garanta tamanhos antes de escrever
        out.fEigenvalues.Resize(3);
        out.fEigenvectors.Resize(3);
        for (int k = 0; k < 3; ++k) out.fEigenvectors[k].Resize(3);

        for (int k = 0; k < 3; k++) {
            out.fEigenvalues[k]    = ces.eigenvalues()(k).real();
            out.fEigenvectors[k][0]= ces.eigenvectors()(0,k).real();
            out.fEigenvectors[k][1]= ces.eigenvectors()(1,k).real();
            out.fEigenvectors[k][2]= ces.eigenvectors()(2,k).real();
        }
    }
};


#endif //TPZYCVonMisesPV
