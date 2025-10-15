#include "TPZYCVonMisesVoigt.h"
#include "pzerror.h"
#include "pzlog.h"
#include "pzvec.h"
#include "pzmanvector.h"
#include "TPZStream.h"
#include <cmath>
#include <iomanip>
#include "TPZKrylovEigenSolver.h"
#include "TPZEigenSolver.h"

TPZYCVonMisesVoigt::TPZYCVonMisesVoigt()
: fSigmaY0(0.0), fH0(0.0)
{

}


TPZYCVonMisesVoigt::TPZYCVonMisesVoigt(const TPZYCVonMisesVoigt& cp)
: TPZPlasticCriterion(cp)
, fSigmaY0(cp.fSigmaY0)
, fH0(cp.fH0)

{

}

void TPZYCVonMisesVoigt::SetUp(STATE sigmaY0, STATE Hiso)
{
    fSigmaY0 = sigmaY0;
    fH0=Hiso;
}

void TPZYCVonMisesVoigt::SetLocalMatState(TPZPlasticState<REAL> & /*state*/)
{
     DebugStop();
}

TPZPlasticState<REAL> TPZYCVonMisesVoigt::GetLocalMatState()
{
    TPZPlasticState<REAL> st;
    return st;
}


void TPZYCVonMisesVoigt::ChangeLocalMatParameters(TPZPlasticState<REAL> & /*state*/, REAL /*factor*/)
{
}

int TPZYCVonMisesVoigt::ClassId() const
{
        return Hash("TPZYCVonMisesVoigt") ;
}

void TPZYCVonMisesVoigt::Read(TPZStream& buf, void* /*context*/)
{
    buf.Read(&fSigmaY0,1);
}

void TPZYCVonMisesVoigt::Write(TPZStream& buf, int /*withclassid*/) const
{
    buf.Write(&fSigmaY0,1);
}

void TPZYCVonMisesVoigt::Print(std::ostream &out) const
{
    out << "----- TPZYCVonMisesVoigt -----\n";
    out << "Yield stress (sigma_y): " << std::setprecision(12) << fSigmaY0 << "\n";
    out << "NYield = " << NYield << "\n";
    out << "---------------------------\n";
}


void TPZYCVonMisesVoigt::Phi(TPZTensor<STATE>sig, STATE alpha, TPZVec<STATE> &phi) const
{
    phi.resize(1);
    STATE j2 =sig.J2();
    STATE q=sqrt(3.* j2);
    STATE f=q-SigmaY(alpha);
    phi[0]=f;
}



STATE TPZYCVonMisesVoigt::ProjectSigma(const TPZTensor<STATE> & sigmatr,  TPZTensor<STATE> & sigmaproj,TPZElasticResponse &ER, STATE &havarn,STATE &havarn1, int & m_type)
{

    fHard=havarn;

    TPZTensor<STATE> S;
    sigmatr.S(S);
    STATE j2 =sigmatr.J2();
    STATE q=sqrt(3.* j2);
    STATE sigy  = SigmaY(havarn);
   // STATE sigy  = fSigmaY0;
    STATE f=q-sigy;

    if(f<0.)
    {
        sigmaproj=sigmatr;
        //sigmaproj+=0.;
        havarn1=havarn;
        m_type=0;//elastic
        return 0.;

    }else
    {
        STATE H=DSigmaYDepsbar(havarn);
        STATE G=ER.G();
        STATE gamma =  (q-sigy)/(3*G+H);
        havarn1=fHard+gamma;
        STATE sigyn1  = SigmaY(havarn);

        TPZTensor<REAL>::TPZDecomposed sig_eigen_system;
        //Decomp3x3( sigmatr,sig_eigen_system);

        sigmatr.EigenSystem(sig_eigen_system);

        TPZManVector<STATE,3> sigprincipal=sig_eigen_system.fEigenvalues;
        TPZManVector<TPZManVector<STATE,3>,3> eigenvectors=sig_eigen_system.fEigenvectors;
        STATE sig1=sigprincipal[0];
        STATE sig2=sigprincipal[1];
        STATE sig3=sigprincipal[2];
        STATE betaproj=atan((sqrt(3.)* (-sig2+sig3))/(-2.* sig1+sig2+sig3));
        //STATE betaproj = (1.0/3.0) *atan2( sqrt(3.0)*(sig2 - sig3), (2.0*sig1 - sig2 - sig3) );

        STATE xiproj=sigmatr.I1()/sqrt(3.);

        STATE rhoproj=sqrt(2./3.)*sigyn1 ;

        TPZManVector<STATE,3> HWCylCoords(3),HWCart(3);
        HWCylCoords[0]=xiproj;
        HWCylCoords[1]=rhoproj;
        HWCylCoords[2]=betaproj;


        TPZHWTools::FromHWCylToPrincipal(HWCylCoords,HWCart);
        sig_eigen_system.fEigenvalues=HWCart;

        TPZTensor<STATE>sigmaproj2(sig_eigen_system);
        sigmaproj=sigmaproj2;

        m_type=1;//plastic

        return gamma;
    }

}

TPZTensor<STATE> TPZYCVonMisesVoigt::ComputeN(const TPZTensor<STATE> stresstensor)const
{
    STATE j2 =stresstensor.J2();
    STATE temp=0.;
    //std::cout << "j2  = "<<j2 <<std::endl;
    if(j2<1.e-12)
    {
        j2=1.e-12;
    }
    temp=sqrt(3.)/(2.* sqrt(j2));

    //std::cout << "temp  = "<<temp <<std::endl;
    TPZTensor<STATE> S;
    stresstensor.S(S);
    S*=temp;
    //std::cout << "S.Norm()" << S.Norm() << std::endl ;
    return S;
}

TPZFMatrix<STATE> TPZYCVonMisesVoigt::GetNdSigma(const TPZTensor<STATE>& sigma) const
{
    TPZFMatrix<STATE> dnds(6,6,0.0);

    // --- Deviador trial S e J2 vindos do TPZTensor ---
    TPZTensor<STATE> S;
    sigma.S(S);                 // S = deviador(sigma)
    const STATE J2 = sigma.J2(); // J2 = 1/2 S:S

    // Proteção numérica
    const STATE eps = 1e-20;
    const STATE J2s = std::max(J2, eps);
    const STATE rootJ2 = std::sqrt(J2s);

    // ===== Projetor deviadorico P na SUA ordem de Voigt =====
    // Você especificou: P = (1/3) * [[2,-1,-1,0,0,0],[-1,2,-1,0,0,0],[-1,-1,2,0,0,0],[0,0,0,6,0,0],[0,0,0,0,6,0],[0,0,0,0,0,6]]
    // OBS: isso equivale a bloco normal (2/3,-1/3,...) e cisalhantes = 2 na diagonal,
    // só que REORDENADO para a convenção [_XX_, _XY_, _XZ_, _YY_, _YZ_, _ZZ_].
    TPZFMatrix<STATE> P(6,6,0.0);
    // bloco "normal" (XX,YY,ZZ) — posições (_XX_, _YY_, _ZZ_) = (0,3,5)
    P(_XX_, _XX_) = 2.0/3.0;  P(_XX_, _YY_) = -1.0/3.0; P(_XX_, _ZZ_) = -1.0/3.0;
    P(_YY_, _XX_) = -1.0/3.0; P(_YY_, _YY_) =  2.0/3.0; P(_YY_, _ZZ_) = -1.0/3.0;
    P(_ZZ_, _XX_) = -1.0/3.0; P(_ZZ_, _YY_) = -1.0/3.0; P(_ZZ_, _ZZ_) =  2.0/3.0;
    // cisalhantes na diagonal (na SUA ordem: XY, XZ, YZ em 1,2,4) → valor 2
    P(_XY_, _XY_) = 2.0;
    P(_XZ_, _XZ_) = 2.0;
    P(_YZ_, _YZ_) = 2.0;

    // ===== vetor s (deviador) em VOIGT na SUA ordem =====
    TPZFMatrix<STATE> svec(6,1,0.0); // coluna 6x1
    svec(_XX_,0) = S.XX();
    svec(_XY_,0) = S.XY(); // atenção: você pediu _XY_ = 1
    svec(_XZ_,0) = S.XZ(); // _XZ_ = 2
    svec(_YY_,0) = S.YY(); // _YY_ = 3
    svec(_YZ_,0) = S.YZ(); // _YZ_ = 4
    svec(_ZZ_,0) = S.ZZ(); // _ZZ_ = 5

    // ===== outer product (s ⊗ s) usando TPZFMatrix =====
    TPZFMatrix<STATE> sT;          // 1x6
    svec.Transpose(&sT);           // sT = s^T
    TPZFMatrix<STATE> s_outer(6,6,0.0);
    svec.Multiply(sT, s_outer);    // s_outer = s * s^T

    // ===== coeficientes =====
    const STATE c1 = std::sqrt(3.0) / (2.0 * rootJ2);
    const STATE c2 = std::sqrt(3.0) / (4.0 * std::pow(J2s, 1.5));

    // d n / d sigma = c1 * P  -  c2 * (s ⊗ s)
    dnds = P;
    dnds *= c1;
    dnds -= c2 * s_outer;

    return dnds;
}

// STATE TPZYCVonMisesVoigt::ComputeGamma(const TPZTensor<STATE>sig, const TPZFMatrix<STATE> elasticmat)const
// {
//     STATE j2 =sig.J2();
//     STATE q=sqrt(3.* j2);
//     STATE f=q-fSigmaY0;
//     TPZTensor<STATE> Nvec = ComputeN(sig);
//     //std::cout << "Nvec  = "<<Nvec <<std::endl;
//     TPZFMatrix<STATE> NFmat(6,1,0.),tempsol,NFmatT,sol;
//     Nvec.CopyTo(NFmat);
//     elasticmat.Multiply(NFmat,tempsol);
//     //tempsol.Print("tempsol");
//     NFmat.Transpose(&NFmatT);
//     NFmatT.Multiply(tempsol,sol);
//     //sol.Print("sol");
//  return f/sol(0,0);
//
// }

