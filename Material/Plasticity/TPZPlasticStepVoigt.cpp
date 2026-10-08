// TPZPlasticStepVoigt.cpp
#include "TPZPlasticStepVoigt.h"
#include "TPZPlasticBase.h"
#include "TPZPlasticCriterion.h"

#include "pzvec.h"

namespace {
/// true if the criterion stores its own elastic response (it must follow the one of the step)
template <class Y, class = void> struct HasSetER : std::false_type {};
template <class Y>
struct HasSetER<Y, std::void_t<decltype(std::declval<Y &>().SetElasticResponse(std::declval<const TPZElasticResponse &>()))>>
    : std::true_type {};
}

// ===================== Construtores / Dtor =====================
template <class YC, class ER>
TPZPlasticStepVoigt<YC,ER>::TPZPlasticStepVoigt(const YC& yc, const ER& er)
: fER(er), fYC(yc)
{
    fN.CleanUp();
    if constexpr (HasSetER<YC>::value) fYC.SetElasticResponse(er); // one elastic response for predictor and projection
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
, fN(other.fN)            // o estado faz parte do objeto (como no operator= implícito)
, fER(other.fER)
, fYC(other.fYC)
, fReductionFactor(other.fReductionFactor)
{
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
    return Hash("TPZPlasticStepVoigt") ^ YC().ClassId() << 1 ^ ER().ClassId() << 2;
}

template <class YC, class ER>
void TPZPlasticStepVoigt<YC,ER>::Write(TPZStream& buf, int withclassid) const
{
    fYC.Write(buf, withclassid);
    fER.Write(buf, withclassid);
    fN.Write(buf, withclassid);
    buf.Write(&fReductionFactor);
}

template <class YC, class ER>
void TPZPlasticStepVoigt<YC,ER>::Read(TPZStream& buf, void* context)
{
    fYC.Read(buf, context);
    fER.Read(buf, context);
    fN.Read(buf, context);
    buf.Read(&fReductionFactor);
}

template <class YC, class ER>
YC TPZPlasticStepVoigt<YC,ER>::LocalCriterion() const
{
    YC yc(fYC);
    // point properties {c, phi, psi}; an all-zero vector is a placeholder (e.g. written by post-processing)
    const bool props = fN.fmatprop.size() >= 3 && (fN.fmatprop[0] != 0. || fN.fmatprop[1] != 0.);
    if (props || fReductionFactor != 1.) {
        TPZPlasticState<REAL> state(fN);
        if (props) yc.SetLocalMatState(state);
        if (fReductionFactor != 1.) yc.ChangeLocalMatParameters(state, fReductionFactor);
    }
    return yc;
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
    TPZTensor<REAL> sigma;
    ApplyStrainComputeSigma(epsTotal, sigma); // updates eps_t, eps_p and the hardening
}

// // // --- ComputeSigma: stub (retorna zero) ---
// template<class YC,class ER>
// void TPZPlasticStepVoigt<YC,ER>::ApplyStrainComputeSigma(const TPZTensor<REAL>& epsTotal,
//                                                          TPZTensor<REAL>& sigma,
//                                                          TPZFMatrix<REAL>* tangent)
// {
//
//     TPZTensor<REAL>sigtrtensor,sigtrtensor2;
//     TPZTensor<REAL> eps_e_trial = epsTotal - fN.m_eps_p;
//     fER.ComputeStress(eps_e_trial,sigtrtensor);
//     TPZFMatrix<STATE> Cmat;
//     fER.De(Cmat);
//
//     int type;
//     TPZFMatrix<STATE> Dep;
//
//
//     STATE hvarnew;
//     STATE gamma=fYC.ProjectSigma(sigtrtensor,sigma,fER,fN.m_hardening,hvarnew,type);
//
//     TPZManVector<STATE,3> sigtrvec =ComputePrincialVal(sigtrtensor);
//     TPZManVector<STATE,3> epstrvecout,sigprojvec;
//     TPZManVector<STATE,2> dlambda;
//     TPZFNMatrix<9> Grad3x3(3,3);
//
//
//     fYC.ProjectSigma(sigtrvec,fN.m_hardening,dlambda,sigprojvec,epstrvecout,Grad3x3,hvarnew,type);
//
//     fN.m_hardening = hvarnew;
//     fN.m_m_type = type;
//
//     if(type==1)//plastico
//     {
//         this->ConsistentTangent(sigtrtensor,sigma,gamma,Dep);
//         TPZManVector<TPZManVector<REAL,3>,3> eigenvetors = ComputePrincialVec(sigtrtensor);
//         ConsistentTangent(sigtrvec,sigprojvec,epstrvecout,Grad3x3,eigenvetors);
//
//     }else{
//         Dep=Cmat;
//     }
//     //Dep=Ce;
//     if (tangent) {
//         *tangent = Dep;
//     }
//     // Reconstruction of sigmaprTensor
//     TPZTensor<REAL> eps_e_Np1;
//     fER.ComputeStrain(sigma, eps_e_Np1);
//     fN.m_eps_t = epsTotal;
//     fN.m_eps_p = epsTotal - eps_e_Np1;
//
// }
// // --- ComputeSigma: stub (retorna zero) ---
template<class YC,class ER>
void TPZPlasticStepVoigt<YC,ER>::ApplyStrainComputeSigma(const TPZTensor<REAL>& epsTotal,
                                                         TPZTensor<REAL>& sigma,
                                                         TPZFMatrix<REAL>* tangent)
{
    // elastic predictor: Voigt strains with engineering shear (paper Eq. 22 and 26)
    TPZTensor<REAL> sigtrtensor;
    fER.ComputeStress(epsTotal - fN.m_eps_p, sigtrtensor);

    TPZTensor<REAL>::TPZDecomposed eigen_system;
    sigtrtensor.EigenSystem(eigen_system);
    TPZManVector<STATE,3> sigtrvec = eigen_system.fEigenvalues, epstrvecout, sigprojvec;
    TPZManVector<STATE,2> dlambda;
    TPZFNMatrix<9> Grad3x3(3,3,0.);
    STATE hvarnew = fN.m_hardening;
    int type = 0;
    YC yc(LocalCriterion());
    yc.ProjectSigma(sigtrvec,fN.m_hardening,dlambda,sigprojvec,epstrvecout,Grad3x3,hvarnew,type);
    fN.m_hardening = hvarnew;
    fN.m_m_type = type;

    if (type == 0) {
        sigma = sigtrtensor; // elastic step: no spectral reconstruction round-off
        if (tangent) fER.De(*tangent);
    } else {
        TPZManVector<TPZManVector<REAL,3>,3> eigenvetors = eigen_system.fEigenvectors;
        // isotropy: principal directions are kept (paper Eq. 62); adding only the plastic correction keeps the
        // eigen-decomposition round-off away from the elastic part
        for (int i = 0; i < 3; i++) eigen_system.fEigenvalues[i] = sigprojvec[i] - sigtrvec[i];
        sigma = sigtrtensor + TPZTensor<REAL>(eigen_system);
        if (tangent) {
            TPZFNMatrix<36> Dep;
            ConsistentTangent(sigtrvec,sigprojvec,epstrvecout,Grad3x3,eigenvetors,Dep);
            *tangent = Dep;
        }
    }
    TPZTensor<REAL> eps_e_Np1;
    fER.ComputeStrain(sigma, eps_e_Np1);
    fN.m_eps_t = epsTotal;
    fN.m_eps_p = epsTotal - eps_e_Np1;
}

template<class YC,class ER>
void TPZPlasticStepVoigt<YC,ER>::ConsistentTangent(TPZManVector<STATE,3>& sigtrial, TPZManVector<STATE,3>& sigproj,TPZManVector<STATE,3>&epstrial, TPZFNMatrix<9> &Grad3x3,TPZManVector<TPZManVector<STATE,3>,3>&eigenvetors, TPZFNMatrix<36>& Dep) const
{

    TPZFNMatrix<36> gradpart(6,6,0.);
    TPZFNMatrix<36> rotationpart(6,6,0.);
    TPZFNMatrix<36> Cmat(6,6,0.),tempmat0;
    fER.De(Cmat);
    STATE G= fER.G();

    for( int icol=0;icol<6;icol++)
    {

        TPZFNMatrix<9> temprot(3,3,0.);
         TPZFNMatrix<9> deltaE = EBasisGrad(icol);
         for(int i=0;i<3;i++)
         {
             for(int j=0;j<3;j++)
             {
                 TPZFNMatrix<6> prodii= FormCartToVoigt(TensorProduct(eigenvetors[i],eigenvetors[i]));
                 TPZFNMatrix<6> prodjj= FormCartToVoigtGrad(TensorProduct(eigenvetors[j],eigenvetors[j]));
                 TPZFNMatrix<9> tempmat=TensorProduct(prodii,prodjj);
                 tempmat*=Grad3x3(i,j);
                 tempmat.Multiply(Cmat,tempmat0);
                 TPZFNMatrix<6> prodeltaEVoigth= FormCartToVoigtGrad(deltaE);
                 tempmat0.Multiply(prodeltaEVoigth,tempmat);

                 for(int irow=0;irow<6;irow++)
                 {
                     gradpart(irow,icol)+=tempmat[irow];
                 }
                 if(i<=j)continue;
                 STATE depstr = (epstrial[i] - epstrial[j]);
                 STATE dsigproj = (sigproj[i] - sigproj[j]);
                 //std::cout<< "tempmat = "<<tempmat<<std::endl;
                 STATE fac=0.;
                 if(fabs(depstr) < 1.e-12)
                 {
                      fac = G*(Grad3x3(i, i) -Grad3x3(i, j) - Grad3x3(j, i) + Grad3x3(j, j));
                 }else{
                     fac = dsigproj/depstr;
                }
                TPZFNMatrix<3> vecj = {{eigenvetors[j][0],eigenvetors[j][1],eigenvetors[j][2]}};
                TPZFNMatrix<3> veci = {{eigenvetors[i][0]},{eigenvetors[i][1]},{eigenvetors[i][2]}};
                TPZFMatrix<STATE> temp1,temp2;
                deltaE.Multiply(veci,temp1);
                // temp1.Print("temp1");
                // vecj.Print("vecj");
                vecj.Multiply(temp1,temp2);
                //temp2.Print("temp2");
                TPZFNMatrix<9>tempmat2=TensorProduct(eigenvetors[i],eigenvetors[j])+TensorProduct(eigenvetors[j],eigenvetors[i]);
                tempmat2*=fac;
                tempmat2*=temp2(0,0);
                temprot+=tempmat2;
                //std::cout <<" temprot = " <<temprot << std::endl;

             }
        }
        for(int irow=0;irow<6;irow++)
        {
            rotationpart(irow,icol)+=FormCartToVoigt(temprot)[irow];
        }
    }
    // std::cout <<" rotationpart = " <<rotationpart << std::endl;
    // std::cout <<" gradpart = " <<gradpart << std::endl;
    Dep=gradpart+rotationpart;
  //  std::cout <<" Dep = " <<Dep << std::endl;

}
// template<class YC,class ER>
// void TPZPlasticStepVoigt<YC,ER>::ConsistentTangent(TPZManVector<STATE,3>& sigtrial, TPZManVector<STATE,3>& sigproj,TPZManVector<STATE,3>&epstrial, TPZFNMatrix<9> &Grad3x3,TPZManVector<TPZManVector<STATE,3>,3>&eigenvetors, TPZFNMatrix<36>& Dep) const
// {
//
//     TPZFNMatrix<36> A(6,6,0.);
//     TPZFNMatrix<36> R(6,6,0.);
//     TPZFNMatrix<36> Cmat(6,6,0.);
//     fER.De(Cmat);
//     STATE G= fER.G();
//     //std::cout<<"Cmat" << Cmat <<std::endl;
//
//     for( int icol=0;icol<6;icol++)
//     {
//
//         TPZFNMatrix<36> tempmat0;
//         TPZFNMatrix<6> ColA(6,1,0.);
//         TPZFNMatrix<6> ColR(6,1,0.);
//         TPZFNMatrix<9> deltaE = EBasis(icol);
//         for(int i=0;i<3;i++)
//         {
//             for(int j=0;j<3;j++)
//             {
//                 TPZFNMatrix<6> prodii= FormCartToVoigt(TensorProduct(eigenvetors[i],eigenvetors[i]));
//                 TPZFNMatrix<6> prodjj= FormCartToVoigt(TensorProduct(eigenvetors[j],eigenvetors[j]));
//                 TPZFNMatrix<36> tempmat=TensorProduct(prodii,prodjj),tempmat00;
//                 //if(icol==5)std::cout<<"Grad3x3(i,j)" << Grad3x3(i,j) <<std::endl;
//                 //if(icol==5)std::cout<< "tempmat" << tempmat <<std::endl;
//                 tempmat*=Grad3x3(i,j);
//                 //if(icol==5)std::cout<< "tempmat*=Grad3x3(i,j)" << tempmat <<std::endl;
//                 tempmat.Multiply(Cmat,tempmat00);
//
//                 //if(icol==5)std::cout<<"Cmat" << Cmat <<std::endl;
//                 //std::cout<< "tempmat" << tempmat <<std::endl;
//                 TPZFNMatrix<6> prodeltaEVoigth= FormCartToVoigt(deltaE);
//                 //if(icol==5)std::cout<<"tempmatafter" << tempmat00 <<std::endl;
//                 //if(icol==5)std::cout<<"prodeltaEVoigth" << prodeltaEVoigth <<std::endl;
//                 tempmat00.Multiply(prodeltaEVoigth,tempmat0);
//                 //if(icol==5)std::cout<< "tempmat0" << tempmat0 <<std::endl;
//
//                 for(int irow=0;irow<6;irow++)
//                 {
//                     ColA(irow,0)+=tempmat0(irow,0);
//                 }
//                 //std::cout<< "ColA" << ColA <<std::endl;
//                 if(j<=i)continue;
//                 STATE depstr = (epstrial[i] - epstrial[j]);
//                 STATE dsigproj = (sigproj[i] - sigproj[j]);
//                 //std::cout<< "tempmat = "<<tempmat<<std::endl;
//                 STATE fac=0.;
//                 if(fabs(depstr) < 1.e-12)
//                 {
//                     fac = G*(Grad3x3(i, i) -Grad3x3(i, j) - Grad3x3(j, i) + Grad3x3(j, j));
//                 }else{
//                     fac = dsigproj/depstr;
//                 }
//                 TPZFNMatrix<9>tempmat2a=TensorProduct(eigenvetors[i],eigenvetors[j]);
//                 //std::cout<< "tempmat2a" << tempmat2a <<std::endl;
//                 TPZFNMatrix<9>tempmat2b=TensorProduct(eigenvetors[j],eigenvetors[i]);
//                // std::cout<< "tempmat2b" << tempmat2b <<std::endl;
//                 TPZFNMatrix<9> tempmat2 = 0.5*(tempmat2a+tempmat2b);
//                 //std::cout<< "tempmat2" << tempmat2 <<std::endl;
//                 TPZFNMatrix<6> sij= FormCartToVoigt(tempmat2);
//                 TPZFNMatrix<36>tempmat3=TensorProduct(sij,sij);
//                 //std::cout<< "tempmat3" << tempmat3 <<std::endl;
//                 tempmat3*=2.*fac;
//                 TPZFNMatrix<6> tempmat4;
//
//                 tempmat3.Multiply(prodeltaEVoigth,tempmat4);
//                 for(int irow=0;irow<6;irow++)
//                 {
//                     ColR(irow,0)+=tempmat4(irow,0);
//                 }
//                 //std::cout<< "ColR" << ColR <<std::endl;
//             }
//         }
//         for(int irow=0;irow<6;irow++)
//         {
//             A(irow,icol)+=ColA(irow,0);
//             R(irow,icol)+=ColR(irow,0);;
//         }
//     }
//      //std::cout <<" R = " <<R << std::endl;
//      //std::cout <<" A = " <<A << std::endl;
//     Dep=A+R;
//     //  std::cout <<" Dep = " <<Dep << std::endl;
//
// }




// --- ComputeDep: stub (mesma ideia, Dep zerada) ---
template<class YC,class ER>
void TPZPlasticStepVoigt<YC,ER>::ApplyStrainComputeDep(const TPZTensor<REAL>& epsTotal,
                                                       TPZTensor<REAL>& sigma,
                                                       TPZFMatrix<REAL>& Dep)
{
    Dep.Redim(6,6);
    ApplyStrainComputeSigma(epsTotal, sigma, &Dep);
}

// --- ApplyLoad: stub (inverso constitutivo a implementar) ---
template<class YC,class ER>
void TPZPlasticStepVoigt<YC,ER>::ApplyLoad(const TPZTensor<REAL>& sigma, TPZTensor<REAL>& epsTotal)
{
    (void)sigma;
    epsTotal.Zero();
    DebugStop(); // not implemented
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
void TPZPlasticStepVoigt<YC,ER>::Phi(const TPZTensor<REAL>& epsElastic, TPZVec<REAL>& phi) const
{
    TPZTensor<REAL> sigma;
    fER.ComputeStress(epsElastic, sigma);
    TPZTensor<REAL>::TPZDecomposed eigen_system;
    sigma.EigenSystem(eigen_system);
    LocalCriterion().YieldFunction(eigen_system.fEigenvalues, fN.m_hardening, phi);
}

// --- ElasticResponse: a base exige TPZElasticResponse exatamente ---
template<class YC,class ER>
void TPZPlasticStepVoigt<YC,ER>::SetElasticResponse(TPZElasticResponse& ERin)
{
    fER = ERin;
    if constexpr (HasSetER<YC>::value) fYC.SetElasticResponse(ERin); // the projection uses the criterion's G, K
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



template class TPZPlasticStepVoigt<TPZYCMohrCoulombPV2, TPZElasticResponse>;
template class TPZPlasticStepVoigt<TPZYCVonMisesVoigt, TPZElasticResponse>;
template class TPZPlasticStepVoigt<TPZYCTrescaVoigt, TPZElasticResponse>;
