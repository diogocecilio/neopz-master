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

    TPZTensor<REAL>sigtrtensor,sigtrtensor2;

    TPZTensor<REAL> eps_e_trial = epsTotal - fN.m_eps_p;


    fER.ComputeStress(eps_e_trial,sigtrtensor);

   // std::cout<< "strial = " << sigtrtensor << std::endl;

    TPZFMatrix<STATE> Cmat;
    fER.De(Cmat);

    int type;
    //TPZFMatrix<STATE> Dep;
    TPZFNMatrix<36> Dep;
    STATE hvarnew;

    TPZTensor<REAL>::TPZDecomposed eigen_system;
    sigtrtensor.EigenSystem(eigen_system);

    TPZManVector<STATE,3> sigtrvec =eigen_system.fEigenvalues;
    TPZManVector<STATE,3> epstrvecout,sigprojvec;
    TPZManVector<STATE,2> dlambda;
    TPZFNMatrix<9> Grad3x3(3,3);


    fYC.ProjectSigma(sigtrvec,fN.m_hardening,dlambda,sigprojvec,epstrvecout,Grad3x3,hvarnew,type);

    fN.m_hardening = hvarnew;
    fN.m_m_type = type;

    if(type==1)//plastico
    {
        TPZManVector<TPZManVector<REAL,3>,3> eigenvetors = eigen_system.fEigenvectors;

        ConsistentTangent(sigtrvec,sigprojvec,epstrvecout,Grad3x3,eigenvetors,Dep);

    }else{
        Dep=Cmat;
    }
    //Dep=Ce;
    if (tangent) {
        *tangent = Dep;
    }
    // Reconstruction of sigmaprTensor
    // Reconstruction of sigmaprTensor
    eigen_system.fEigenvalues = sigprojvec; // Under the assumption of isotropic material eigen vectors remain unaltered
    sigma = TPZTensor<REAL>(eigen_system);
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



template<class YC, class ER>
void TPZPlasticStepVoigt<YC,ER>::ConsistentTangent(const TPZTensor<STATE>& sigmatr,const TPZTensor<STATE>& sigmapr,STATE gamma, TPZFMatrix<STATE>& Dep) const
{

    TPZTensor<STATE> Nvec = fYC.ComputeN(sigmapr);


    TPZFMatrix<STATE> dadsig = fYC.GetNdSigma(sigmapr); // 6x6 na sua ordem

    TPZFMatrix<REAL>  Ce,invCe;
    fER.De(Ce) ;
    // Ce(_XY_,_XY_)/=2.;
    // Ce(_XZ_,_XZ_)/=2.;
    // Ce(_YZ_,_YZ_)/=2.;

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

    temp*=1./(tempreal(0,0)+fYC.DSigmaYDepsbar(fN.m_hardening));

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



template class TPZPlasticStepVoigt<TPZYCMohrCoulombPV2, TPZElasticResponse>;
template class TPZPlasticStepVoigt<TPZYCVonMisesVoigt, TPZElasticResponse>;
