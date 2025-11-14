// TPZPlasticStepVoigt.h
#ifndef TPZPLASTICSTEPVOIGT_H
#define TPZPLASTICSTEPVOIGT_H

#include "TPZPlasticBase.h"
#include "pzvec.h"
#include <memory>
#include "TPZYCVonMisesVoigt.h"
#include "TPZYCMohrCoulombPV2.h"
#include "TPZElasticResponse.h"
/// Classe constitutiva elasto-plástica em Voigt (3D) [11,12,13,22,23,33],
/// com critério de escoamento genérico YC e resposta elástica ER,
/// no mesmo espírito do TPZPlasticStepPV (mas operando em Voigt).
template <class YC, class ER>
class TPZPlasticStepVoigt : public TPZPlasticBase
{
public:
    //using YieldCriterion = YC;
    //using ElasticResponse = ER;

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


    // (A) Forte/Tipada: barata e infalível
    void SetPlasticCriterion(const YC& pc);

    // (B) Polimórfica: aceita a classe base, verifica compatibilidade em runtime
    void SetPlasticCriterion(const TPZPlasticCriterion& pc_base);


    TPZTensor<STATE> FromFMatToTensor(TPZFMatrix<STATE> mat);

    // Calcula a matriz tangente consistente pela fórmula:
    // gamma = yield / (a^T Ce a)
    // Q = I + gamma Ce dadsig
    // R = Q^{-1} Ce
    // Dep = R - ( (R a) ⊗ (R a) ) / (a^T R a)
    void ConsistentTangent(const TPZTensor<STATE>& sigmatr,const TPZTensor<STATE>& sigmapr,STATE gamma,TPZFMatrix<STATE>& Dep) const;

    void ConsistentTangent(TPZManVector<STATE,3>& sigtrial, TPZManVector<STATE,3>& sigproj,TPZManVector<STATE,3>&epstrial, TPZFNMatrix<9> &Grad3x3,TPZManVector<TPZManVector<STATE,3>,3>&eigenvetors, TPZFNMatrix<36>& Dep) const;
    TPZTensor<STATE> MultiplyMatrixTensor(const TPZFMatrix<STATE> &A, const TPZTensor<STATE> &S);

    TPZManVector<REAL, 3> ComputePrincialVal(TPZTensor<STATE> &tensor)
    {
        TPZTensor<REAL>::TPZDecomposed eigen_system;
        tensor.EigenSystem(eigen_system);
        return eigen_system.fEigenvalues;
    }
    TPZManVector<TPZManVector<REAL,3>,3> ComputePrincialVec(TPZTensor<STATE> &tensor)
    {
        TPZTensor<REAL>::TPZDecomposed eigen_system;
        tensor.EigenSystem(eigen_system);
        return eigen_system.fEigenvectors;
    }

    TPZFMatrix< STATE > TensorProduct(const TPZManVector< STATE > a, const TPZManVector< STATE > b)const
    {
        int sza=a.size();
        int szb=a.size();
        TPZFMatrix<STATE> out(sza,szb,0.);
        for(int i=0;i<sza;i++)
        {
            for(int j=0;j<szb;j++)
            {
               out(i,j)= a[i]*b[j];
            }

        }
        return out;
    }
    TPZFMatrix< STATE > TensorProduct(const TPZFNMatrix<6> a, const TPZFNMatrix<6> b)const
    {
        int sza=6;
        int szb=6;
        TPZFMatrix<STATE> out(sza,szb,0.);
        for(int i=0;i<sza;i++)
        {
            for(int j=0;j<szb;j++)
            {
                out(i,j)= a[i]*b[j];
            }

        }
        return out;
    }

    TPZFNMatrix<9> EBasisGrad(int k) const
    {
        TPZFNMatrix<9>Base;
        switch (k) {
            case _XX_:
                Base={{1., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
                break;
            case _XY_:
                 Base={{0., 1./2., 0.}, {1./2., 0., 0.}, {0., 0., 0.}};
                break;
            case _XZ_:
                Base={{0., 0., 1./2.}, {0., 0., 0.}, {1./2., 0., 0.}};
                break;
            case _YY_:
                Base={{0., 0., 0.}, {0., 1., 0.}, {0., 0., 0.}};
                break;
            case _YZ_:
                Base={{0., 0., 0.}, {0., 0., 1./2.}, {0., 1./2., 0.}};
                break;
            case _ZZ_:
                Base={{0., 0., 0.}, {0., 0., 0.}, {0., 0., 1.}};
                break;

            default:
                DebugStop();
                break; // Optional for the last case/default
        }

        return Base;

    }
    TPZFNMatrix<9> EBasis(int k) const
    {
        TPZFNMatrix<9>Base;
        switch (k) {
            case _XX_:
                Base={{1., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
                break;
            case _XY_:
                Base={{0., 1., 0.}, {1., 0., 0.}, {0., 0., 0.}};
                break;
            case _XZ_:
                Base={{0., 0., 1.}, {0., 0., 0.}, {1., 0., 0.}};
                break;
            case _YY_:
                Base={{0., 0., 0.}, {0., 1., 0.}, {0., 0., 0.}};
                break;
            case _YZ_:
                Base={{0., 0., 0.}, {0., 0., 1.}, {0., 1., 0.}};
                break;
            case _ZZ_:
                Base={{0., 0., 0.}, {0., 0., 0.}, {0., 0., 1.}};
                break;

            default:
                DebugStop();
                break; // Optional for the last case/default
        }

        return Base;

    }
    TPZFNMatrix<6> FormCartToVoigt(TPZFNMatrix<9> cart)const
    {
        if(cart.Rows()!=3||cart.Cols()!=3)
        {
            DebugStop();
        }
        TPZFNMatrix<6> voigt = {{cart(0, 0)}, {cart(0, 1)},{cart(0, 2)},{cart(1, 1)}, {cart(1, 2)},{cart(2, 2)}};

        return voigt;
    }

    TPZFNMatrix<6> FormCartToVoigtGrad(TPZFNMatrix<9> cart)const
    {
        if(cart.Rows()!=3||cart.Cols()!=3)
        {
            DebugStop();
        }
       TPZFNMatrix<6> voigt = {{cart(0, 0)}, {2*cart(0, 1)},{2*cart(0, 2)},{cart(1, 1)}, {2*cart(1, 2)},{cart(2, 2)}};

        return voigt;
    }
    TPZPlasticState<REAL> fN;
protected:
    ER   fER;                 // resposta elástica (p.ex., armazena K,G ou E,nu)
    YC   fYC;                 // critério de escoamento (deve derivar de TPZPlasticCriterion)
    //TPZPlasticState<REAL> fN;
};





#endif // TPZPLASTICSTEPVOIGT_H
