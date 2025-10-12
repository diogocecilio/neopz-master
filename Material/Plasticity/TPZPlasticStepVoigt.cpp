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
    fER.ComputeStress(epsTotal,sigtrtensor);

    int type;
    TPZFMatrix<STATE> Ce,Dep;
    fER.De(Ce);
    REAL hardening=0.;
    //std::cout<< "epsTotal= "<< epsTotal << std::endl;
    STATE gamma=fYC.UpdateHardeningVar(sigtrtensor, fER,hardening);
    fN.m_hardening=hardening;
    //std::cout<< "fN.m_hardening DEPOIS = "<< fN.m_hardening << std::endl;
    fYC.ProjectSigma(sigtrtensor,sigma,fN.m_hardening,type);

    fN.m_m_type = type;


    if(type==1)//plastico
    {
         //std::cout<< " ==================== = \n";
         //std::cout<< "sigtrtensor  ="<<sigtrtensor << std::endl;
         //std::cout<< "sigma = "<< sigma << std::endl;
        //std::cout<< "gamma = "<< gamma << std::endl;
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
// #define DEBUG_PLASTIC
// #ifdef DEBUG_PLASTIC
// #include <iomanip>
// #include <iostream>
// #endif
//
// template<class YC,class ER>
// void TPZPlasticStepVoigt<YC,ER>::ApplyStrainComputeSigma(const TPZTensor<REAL>& epsTotal,
//                                                          TPZTensor<REAL>& sigma,
//                                                          TPZFMatrix<REAL>* tangent)
// {
//     #ifdef DEBUG_PLASTIC
//     std::cout.setf(std::ios::scientific);
//     std::cout << std::setprecision(6);
//     std::cout << "\n=== ApplyStrainComputeSigma ===\n";
//     std::cout << "epsTotal = [" << epsTotal.XX() << ", " << epsTotal.YY() << ", " << epsTotal.ZZ()
//     << ", " << epsTotal.XY() << ", " << epsTotal.YZ() << ", " << epsTotal.XZ() << "]\n";
//     #endif
//
//     // 1) Preditor elástico (trial)
//     TPZTensor<REAL> sigtrtensor;
//     fER.ComputeStress(epsTotal, sigtrtensor);
//
//     #ifdef DEBUG_PLASTIC
//     std::cout << "sig_trial = [" << sigtrtensor.XX() << ", " << sigtrtensor.YY() << ", " << sigtrtensor.ZZ()
//     << ", " << sigtrtensor.XY() << ", " << sigtrtensor.YZ() << ", " << sigtrtensor.XZ() << "]\n";
//     // I1, desviador, J2 e sigma_eq do trial
//     REAL I1tr = sigtrtensor.I1();
//     TPZTensor<REAL> str = sigtrtensor;
//     const REAL ptr = I1tr/3.0;
//     str.XX() -= ptr; str.YY() -= ptr; str.ZZ() -= ptr;
//     REAL ss_tr =  str.XX()*str.XX() + str.YY()*str.YY() + str.ZZ()*str.ZZ()
//     + 2.0*( str.XY()*str.XY() + str.YZ()*str.YZ() + str.XZ()*str.XZ() );
//     REAL J2tr  = 0.5*ss_tr;
//     REAL seqtr = std::sqrt(3.0*J2tr);
//     std::cout << "I1_trial = " << I1tr << "  J2_trial = " << J2tr << "  seq_trial = " << seqtr << "\n";
//     #endif
//
//     // 2) Tensor elástico e hardening
//     int type = -1;
//     TPZFMatrix<STATE> Ce, Dep;
//     fER.De(Ce);
//
//     #ifdef DEBUG_PLASTIC
//     if (Ce.Rows()>=6 && Ce.Cols()>=6) {
//         std::cout << "Ce diag: C11=" << Ce(0,0) << " C22=" << Ce(1,1)
//         << " C33=" << Ce(2,2) << " C44=" << Ce(3,3)
//         << " C55=" << Ce(4,4) << " C66=" << Ce(5,5) << "\n";
//     }
//     std::cout << "hard_var(before) = " << fN.m_hardening << "\n";
//     #endif
//
//     REAL hardening = 0.;
//     STATE gamma = fYC.UpdateHardeningVar(sigtrtensor, Ce, hardening);
//     fN.m_hardening += hardening;
//
//     #ifdef DEBUG_PLASTIC
//     std::cout << "UpdateHardeningVar: gamma=" << gamma
//     << "  hardening_inc=" << hardening
//     << "  hard_var(after)=" << fN.m_hardening << "\n";
//     #endif
//
//     // 3) Projeção (retorno) e tipo (0 elástico, 1 plástico)
//     fYC.ProjectSigma(sigtrtensor, sigma, fN.m_hardening, type);
//     fN.m_m_type = type;
//
//     #ifdef DEBUG_PLASTIC
//     std::cout << "ProjectSigma: type=" << type
//     << "  sigma_corr = [" << sigma.XX() << ", " << sigma.YY() << ", " << sigma.ZZ()
//     << ", " << sigma.XY() << ", " << sigma.YZ() << ", " << sigma.XZ() << "]\n";
//     // checar consistência no fim (f≈0)
//     REAL I1c = sigma.I1();
//     TPZTensor<REAL> sc = sigma;
//     const REAL pc = I1c/3.0;
//     sc.XX() -= pc; sc.YY() -= pc; sc.ZZ() -= pc;
//     REAL ss_c =  sc.XX()*sc.XX() + sc.YY()*sc.YY() + sc.ZZ()*sc.ZZ()
//     + 2.0*( sc.XY()*sc.XY() + sc.YZ()*sc.YZ() + sc.XZ()*sc.XZ() );
//     REAL J2c  = 0.5*ss_c;
//     REAL seqc = std::sqrt(3.0*J2c);
//     std::cout << "I1_corr = " << I1c << "  J2_corr = " << J2c << "  seq_corr = " << seqc << "\n";
//     // Se existir fYC.Phi(σ,hard), você pode conferir f_trial e f_n1 aqui.
//     #endif
//
//     // 4) Tangente consistente
//     if (type == 1) {
//         this->ConsistentTangent(sigtrtensor, sigma, gamma, Dep);
//     } else {
//         Dep = Ce;
//     }
//     if (tangent) { *tangent = Dep; }
//
//     // 5) Atualização de strains elástica e plástica
//     TPZTensor<REAL> eps_e_Np1;
//     fER.ComputeStrain(sigma, eps_e_Np1);   // inverte a lei elástica
//     fN.m_eps_t = epsTotal;
//     fN.m_eps_p = epsTotal - eps_e_Np1;
//
//     #ifdef DEBUG_PLASTIC
//     std::cout << "eps_e(n+1) = [" << eps_e_Np1.XX() << ", " << eps_e_Np1.YY() << ", " << eps_e_Np1.ZZ()
//     << ", " << eps_e_Np1.XY() << ", " << eps_e_Np1.YZ() << ", " << eps_e_Np1.XZ() << "]\n";
//     std::cout << "eps_p(n+1) = [" << fN.m_eps_p.XX() << ", " << fN.m_eps_p.YY() << ", " << fN.m_eps_p.ZZ()
//     << ", " << fN.m_eps_p.XY() << ", " << fN.m_eps_p.YZ() << ", " << fN.m_eps_p.XZ() << "]\n";
//     std::cout << "tr(eps_p) = " << (fN.m_eps_p.XX()+fN.m_eps_p.YY()+fN.m_eps_p.ZZ())
//     << " (≈0 em J2)\n";
//     std::cout << "=== end ApplyStrainComputeSigma ===\n";
//     #endif
// }

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

    // STATE gammafYC.UpdateHardeningVar(sigmapr, Ce,fN.m_hardening);
    //STATE gamma = fYC.ComputeGamma(sigmatr,Ce );


    // Ce.Print("Ce");
    // dadsig.Print("dadsig");
    // std::cout << "gamma = " << gamma <<std::endl;

    //Q=(IdentityMatrix[6]+gamma Ce . dadsigg);

    TPZFMatrix<STATE> Cedadsig,Q(6,6,0.),Qinv,R;
    Q.Identity();
    Ce.Multiply(dadsig,Cedadsig);
    //Cedadsig.Print("Cedadsig");
    Cedadsig*=gamma;
    Q+=Cedadsig;

    Q.Inverse(Qinv,ELU);

    Qinv.Multiply(Ce,R);

   // R.Print("R");

    //std::cout << "a" <<std::endl;
    TPZFMatrix<STATE> Noriginal(6,1,0.),RN,RNt,temp,tempreal;

    Nvec.CopyTo(Noriginal);

    //std::cout << "b" <<std::endl;
    //Noriginal.Print("Noriginal");

    R.Multiply(Noriginal,RN);

   // RN.Print("RN");

    RN.Transpose(&RNt);

     //RNt.Print("RNt");

    RN.Multiply(RNt,temp);

    //temp.Print("temp");

    RNt.Multiply(Noriginal,tempreal);

    //tempreal.Print("tempreal");

    temp*=1./tempreal(0,0);

    Dep=R;

    Dep-=temp;

    //Dep.Print("Dep");
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
