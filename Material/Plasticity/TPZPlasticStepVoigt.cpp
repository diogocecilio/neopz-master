// TPZPlasticStepVoigt.cpp
#include "TPZPlasticStepVoigt.h"
#include "TPZPlasticBase.h"
#include "TPZPlasticCriterion.h"

#include "pzvec.h"

// ===================== Construtores / Dtor =====================
template <class YC, class ER>
TPZPlasticStepVoigt<YC,ER>::TPZPlasticStepVoigt(const YC& yc, const ER& er)
: fER(er), fYC(yc)
{
    fN.CleanUp();
}

// --- Construtor padrão ---
template <class YC, class ER>
TPZPlasticStepVoigt<YC,ER>::TPZPlasticStepVoigt()
: TPZPlasticBase()        // copia/aciona estado base, se houver
, fER()                   // ER default-constructed
, fYC()                   // YC default-constructed

{
    fN.CleanUp();
    // nada além do init de membros
}

// --- Construtor de cópia ---
template <class YC, class ER>
TPZPlasticStepVoigt<YC,ER>::TPZPlasticStepVoigt(const TPZPlasticStepVoigt& other)
: TPZPlasticBase(other)   // copia parte base
, fER(other.fER)
, fYC(other.fYC)
{
    fN.CleanUp();
    // nada extra
}

// --- Destrutor ---
template <class YC, class ER>
TPZPlasticStepVoigt<YC,ER>::~TPZPlasticStepVoigt()
{
    // sem recursos especiais para liberar
}
template <class YC, class ER>
int TPZPlasticStepVoigt<YC,ER>::ClassId() const
{
    return 0;
}

template <class YC, class ER>
void TPZPlasticStepVoigt<YC,ER>::Write(TPZStream& buf, int withclassid) const
{
}

template <class YC, class ER>
void TPZPlasticStepVoigt<YC,ER>::Read(TPZStream& buf, void* context)
{
}

template <class YC, class ER>
void TPZPlasticStepVoigt<YC,ER>::Print(std::ostream& out) const
{

        out << "\n" << this->Name();
        out << "\n YC_t:";
        fYC.Print(out);
        out << "\n ER_t:";
        fER.Print(out);
        out << "\nTPZPlasticStepPV Internal members:";
        out << "\n fN = "; // PlasticState
        fN.Print(out);
}

template <class YC, class ER>
const char* TPZPlasticStepVoigt<YC,ER>::Name() const
{
    return "TPZPlasticStepVoigt";
}

// --- ApplyStrain: stub (apenas grava deformação total no estado, se desejar) ---
template<class YC,class ER>
void TPZPlasticStepVoigt<YC,ER>::ApplyStrain(const TPZTensor<REAL>& epsTotal)
{
    fN.m_eps_t = epsTotal; // se o TPZPlasticState tiver esse campo; ajuste conforme seu struct
}

// // --- ComputeSigma: stub (retorna zero) ---
template<class YC,class ER>
void TPZPlasticStepVoigt<YC,ER>::ApplyStrainComputeSigma(const TPZTensor<REAL>& epsTotal,
                                                         TPZTensor<REAL>& sigma,
                                                         TPZFMatrix<REAL>* tangent)
{

    TPZTensor<REAL>sigtrtensor;
    TPZTensor<REAL> eps_e_trial = epsTotal - fN.m_eps_p;
    fER.ComputeStress(eps_e_trial,sigtrtensor);

    int type;
    TPZFMatrix<STATE> Ce,Dep;
    fER.De(Ce);

    STATE hvarnew;
    STATE gamma=fYC.ProjectSigma(sigtrtensor,sigma,fER,fN.m_hardening,hvarnew,type);
    fN.m_hardening = hvarnew;
    fN.m_m_type = type;

    if(type==1)//plastico
    {
        this->ConsistentTangent(sigtrtensor,sigma,gamma,Dep);

    }else{
        Dep=Ce;
    }
    if (tangent) {
        *tangent = Dep;
    }
    // Reconstruction of sigmaprTensor
    TPZTensor<REAL> eps_e_Np1;
    fER.ComputeStrain(sigma, eps_e_Np1);
    fN.m_eps_t = epsTotal;
    fN.m_eps_p = epsTotal - eps_e_Np1;

}

template<class YC, class ER>
void TPZPlasticStepVoigt<YC,ER>::ConsistentTangent(const TPZTensor<STATE>& sigmatr,const TPZTensor<STATE>& sigmapr,STATE gamma, TPZFMatrix<STATE>& Dep) const
{

    TPZTensor<STATE> Nvec = fYC.ComputeN(sigmapr);


    TPZFMatrix<STATE> dadsig = fYC.GetNdSigma(sigmapr); // 6x6 na sua ordem

    TPZFMatrix<REAL>  Ce,invCe;
    fER.De(Ce) ;
     Ce(_XY_,_XY_)/=2.;
     Ce(_XZ_,_XZ_)/=2.;
     Ce(_YZ_,_YZ_)/=2.;

    //Q=(IdentityMatrix[6]+gamma Ce . dadsigg);

    TPZFMatrix<STATE> Cedadsig,Q(6,6,0.),Qinv,R;
    Q.Identity();
    Ce.Multiply(dadsig,Cedadsig);
    //Cedadsig.Print("Cedadsig");
    Cedadsig*=gamma;
    Q+=Cedadsig;

    Q.Inverse(Qinv,ELU);

    Qinv.Multiply(Ce,R);

    TPZFMatrix<STATE> Noriginal(6,1,0.),RN,RNt,temp,tempreal;

    Nvec.CopyTo(Noriginal);

    R.Multiply(Noriginal,RN);

    RN.Transpose(&RNt);

    RN.Multiply(RNt,temp);

    RNt.Multiply(Noriginal,tempreal);

    temp*=1./tempreal(0,0);

    Dep=R;

    Dep-=temp;

   // Dep= R-1/(asol . R . asol) Outer[Times,R . asol,R . asol];

}


// --- ComputeDep: stub (mesma ideia, Dep zerada) ---
template<class YC,class ER>
void TPZPlasticStepVoigt<YC,ER>::ApplyStrainComputeDep(const TPZTensor<REAL>& epsTotal,
                                                       TPZTensor<REAL>& sigma,
                                                       TPZFMatrix<REAL>& Dep)
{
    (void)epsTotal;
    sigma.Zero();
    Dep.Redim(6,6); Dep.Zero();
}

// --- ApplyLoad: stub (inverso constitutivo a implementar) ---
template<class YC,class ER>
void TPZPlasticStepVoigt<YC,ER>::ApplyLoad(const TPZTensor<REAL>& sigma, TPZTensor<REAL>& epsTotal)
{
    (void)sigma;
    epsTotal.Zero(); // TODO
}

// --- Set/GetState exigidos pela base ---
template<class YC,class ER>
void TPZPlasticStepVoigt<YC,ER>::SetState(const TPZPlasticState<REAL>& state)
{
    fN = state;
}

template<class YC,class ER>
TPZPlasticState<REAL> TPZPlasticStepVoigt<YC,ER>::GetState() const
{
    return fN;
}
template<class YC, class ER>
inline TPZPlasticCriterion& TPZPlasticStepVoigt<YC,ER>::GetYC()
{
    return static_cast<TPZPlasticCriterion&>(fYC);
}
// --- Phi: stub (devolva um vetor com o(s) valor(es) de f; aqui 1 componente = 0) ---
template<class YC,class ER>
void TPZPlasticStepVoigt<YC,ER>::Phi(const TPZTensor<REAL>& epsTotal, TPZVec<REAL>& phi) const
{
    (void)epsTotal;
    phi.Resize(1); phi[0] = 0.; // TODO: usar fYC
}

// --- ElasticResponse: a base exige TPZElasticResponse exatamente ---
template<class YC,class ER>
void TPZPlasticStepVoigt<YC,ER>::SetElasticResponse(TPZElasticResponse& ERin)
{
    fER = ERin;
}

// (A) Forte/Tipada: barata e infalível
template<class YC,class ER>
void TPZPlasticStepVoigt<YC,ER>::SetPlasticCriterion(const YC& pc)
{
    fYC = pc;
}

// (B) Polimórfica: aceita a classe base, verifica compatibilidade em runtime
template<class YC,class ER>
void TPZPlasticStepVoigt<YC,ER>::SetPlasticCriterion(const TPZPlasticCriterion& pc_base)
{
    // tenta cast direto
    if (const auto* ycp = dynamic_cast<const YC*>(&pc_base)) {
        fYC = *ycp;
        return;
    }
    #ifdef PZDEBUG
    PZError << "TPZPlasticStepVoigt::SetPlasticCriterion: "
    "critério incompatível com YC do template.\n";
    DebugStop();
    #else
    throw std::invalid_argument(
        "TPZPlasticStepVoigt::SetPlasticCriterion: tipo incompatível");
    #endif
}


template<class YC,class ER>
TPZElasticResponse TPZPlasticStepVoigt<YC,ER>::GetElasticResponse() const
{
    return fER;
}

template<class YC,class ER>
TPZTensor<STATE> TPZPlasticStepVoigt<YC,ER>::FromFMatToTensor(TPZFMatrix<STATE> mat)
{
    if(mat.Cols()>1)
    {
        PZError << "FromFMatToTensor: matriz incompativel com tranformacao para TPZTensor \n";
        DebugStop();
    }
    TPZTensor<STATE> localt;
    localt.XX()=mat(_XX_,0);
    localt.YY()=mat(_YY_,0);
    localt.ZZ()=mat(_ZZ_,0);
    localt.XY()=mat(_XY_,0);
    localt.XZ()=mat(_XZ_,0);
    localt.YZ()=mat(_YZ_,0);
    return localt;
}

template class TPZPlasticStepVoigt<TPZYCVonMisesVoigt, TPZElasticResponse>;
