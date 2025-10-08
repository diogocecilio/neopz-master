#include "TPZYCVonMisesVoigt.h"
#include "pzerror.h"
#include "pzlog.h"
#include "pzvec.h"
#include "pzmanvector.h"
#include "TPZStream.h"
#include <cmath>
#include <iomanip>



// --------- util local: obtém λ e μ a partir de E e ν ----------
static inline void LameFromER(const TPZElasticResponse &ER, REAL &lambda, REAL &mu)
{
    // Ajuste os nomes abaixo conforme sua versão de TPZElasticResponse
    const REAL E  = ER.E();          // ou ER.E()
    const REAL nu = ER.Poisson();   // ou ER.Nu()

    mu     = E/(2.0*(1.0 + nu));
    lambda = (E*nu)/((1.0 + nu)*(1.0 - 2.0*nu));
}

// ======================================================================
//                              CTORs / DTORs
// ======================================================================

TPZYCVonMisesVoigt::TPZYCVonMisesVoigt() :
fSigmaY0(0.0)//, fER()
{
    // vazio
}

TPZYCVonMisesVoigt::TPZYCVonMisesVoigt(STATE sigmaY0, STATE Hiso)
: fSigmaY0(sigmaY0)
{

    // σy(κ) = σy0 + Hiso * κ
    fSigmaY = [=](STATE kappa) { return sigmaY0 + Hiso * kappa; };

    // H(κ) ≡ Hiso (constante)
    fH = [=](STATE /*kappa*/) { return Hiso; };
}
TPZYCVonMisesVoigt::TPZYCVonMisesVoigt(const TPZYCVonMisesVoigt& cp)
: TPZPlasticCriterion(cp)
, fSigmaY0(cp.fSigmaY0)
, fSigmaY(cp.fSigmaY)
, fH(cp.fH)
{
}

TPZYCVonMisesVoigt & TPZYCVonMisesVoigt::operator=(const TPZYCVonMisesVoigt &cp)
{
    if(this != &cp)
    {
        TPZPlasticCriterion::operator=(cp);
        fSigmaY0 = cp.fSigmaY0;
        //fER = cp.fER;
    }
    return *this;
}

// ======================================================================
//                         Estado/local e parâmetros
// ======================================================================

void TPZYCVonMisesVoigt::SetUp(STATE sigmaY0, STATE Hiso)
{
    fSigmaY0     = sigmaY0;

    // σy = σy0 + Hiso * κ
    fSigmaY = [=](STATE kappa){ return sigmaY0 + Hiso * kappa; };
    // H(κ) ≡ Hiso
    fH      = [=](STATE){ return Hiso; };
}

void TPZYCVonMisesVoigt::SetLocalMatState(TPZPlasticState<REAL> & /*state*/)
{
    // Caso sua formulação precise guardar/alterar estado interno local,
    // implemente aqui. Por ora, sem estado interno específico.
     DebugStop(); // remova o DebugStop do header para compilar.
}

TPZPlasticState<REAL> TPZYCVonMisesVoigt::GetLocalMatState()
{
    // Retorne o estado interno local, se houver. Aqui devolvemos um default.
    TPZPlasticState<REAL> st;
    return st;
}


void TPZYCVonMisesVoigt::ChangeLocalMatParameters(TPZPlasticState<REAL> & /*state*/, REAL /*factor*/)
{
    // Atualize parâmetros locais (ex.: hardening), se necessário.
    // Mantido como stub intencional.
}

int TPZYCVonMisesVoigt::ClassId() const
{
    // Ajuste se sua infraestrutura exigir um ID fixo/Hash específico.
    // Retornar um valor estável é suficiente para muitos casos.
        return Hash("TPZYCVonMisesVoigt") ;
}

void TPZYCVonMisesVoigt::Read(TPZStream& buf, void* /*context*/)
{
    buf.Read(&fSigmaY0,1);
    //fER.Read(buf,nullptr);
}

void TPZYCVonMisesVoigt::Write(TPZStream& buf, int /*withclassid*/) const
{
    buf.Write(&fSigmaY0,1);
    //fER.Write(buf,0);
}

void TPZYCVonMisesVoigt::Print(std::ostream &out) const
{
    out << "----- TPZYCVonMisesVoigt -----\n";
    out << "Yield stress (sigma_y): " << std::setprecision(12) << fSigmaY0 << "\n";
    //out << "Elastic response (E, nu): E=" << fER.E() << "  nu=" << fER.Poisson() << "\n";
    out << "NYield = " << NYield << "\n";
    out << "---------------------------\n";
}

// ======================================================================
//                 Phi (função de escoamento) e helpers elásticos
// ======================================================================

void TPZYCVonMisesVoigt::Phi(TPZTensor<STATE>sig, STATE alpha, TPZVec<STATE> &phi) const
{
    phi.resize(1);
    STATE j2 =sig.J2();
    STATE q=sqrt(3.* j2);
    STATE f=q-fSigmaY0;
    phi[0]=f;
}



void TPZYCVonMisesVoigt::ProjectSigma(const TPZTensor<STATE> & sigmatr, STATE k_prev, TPZTensor<STATE> & sigmaproj, STATE &k_proj, int & m_type)
{


    STATE p =sigmatr.I1()/3.;
    TPZTensor<STATE> S;
    sigmatr.S(S);
    STATE j2 =sigmatr.J2();
    STATE q=sqrt(3.* j2);
    STATE f=q-fSigmaY0;
  //  std::cout<< "p = "<< p << std::endl;
  //  std::cout<< "S = "<< S << std::endl;
  //  std::cout<< "j2 = "<< j2 << std::endl;
   // std::cout<< "q = "<< q << std::endl;
  //  std::cout<< "f = "<< f << std::endl;
    if(f<0.)
    {
        sigmaproj=sigmatr;
        k_proj=k_prev;
        m_type=0;//elastic

    }else
    {

        TPZTensor<REAL>::TPZDecomposed sig_eigen_system;
        sigmatr.EigenSystem(sig_eigen_system);
        TPZManVector<STATE,3> sigprincipal=sig_eigen_system.fEigenvalues;
        TPZManVector<TPZManVector<STATE,3>,3> eigenvectors=sig_eigen_system.fEigenvectors;
        STATE sig1=sigprincipal[0];
        STATE sig2=sigprincipal[1];
        STATE sig3=sigprincipal[2];
        STATE betaproj=atan((sqrt(3.)* (-sig2+sig3))/(-2.* sig1+sig2+sig3));

        STATE xiproj=sigmatr.I1()/sqrt(3.);

        STATE rhoproj=sqrt(2./3.)*fSigmaY0 ;

        TPZManVector<STATE,3> HWCylCoords(3),HWCart(3);
        HWCylCoords[0]=xiproj;
        HWCylCoords[1]=rhoproj;
        HWCylCoords[2]=betaproj;


        TPZHWTools::FromHWCylToPrincipal(HWCylCoords,HWCart);

       // std::cout<< " sig1 = "<< sig1 << " sig2 = "<< sig2 <<" sig3 = "<< sig3 <<std::endl;
      //  std::cout<< " betasol = "<< betaproj <<std::endl;
      //  std::cout<< " eigenvectors = "<< eigenvectors <<std::endl;
     //  std::cout<< " HWCart = "<< HWCart <<std::endl;
        TPZFMatrix<REAL> v1,v2,v3,v1t,v2t,v3t,temp1,temp2,temp3;
        v1.CopyFrom(eigenvectors[0]);
        v2.CopyFrom(eigenvectors[1]);
        v3.CopyFrom(eigenvectors[2]);
        v1.Transpose(&v1t);
        v2.Transpose(&v2t);
        v3.Transpose(&v3t);
        v1.Multiply(v1t,temp1);
        v2.Multiply(v2t,temp2);
        v3.Multiply(v3t,temp3);

        temp1*=HWCart[0];
        temp2*=HWCart[1];
        temp3*=HWCart[2];
        temp1+=temp2;
        temp1+=temp3;

       // temp1.Print("sigma projected recontructed");

        sigmaproj.XX()=temp1(0,0);sigmaproj.XY()=temp1(0,1);sigmaproj.XZ()=temp1(0,2);
        sigmaproj.XY()=temp1(1,0);sigmaproj.YY()=temp1(1,1);sigmaproj.YZ()=temp1(1,2);
        sigmaproj.XZ()=temp1(2,0);sigmaproj.YZ()=temp1(2,1);sigmaproj.ZZ()=temp1(2,2);

        k_proj=0.;
    //    std::cout<< " sigmaproj = "<< sigmaproj <<std::endl;

        m_type=1;//plastic
    }


}
TPZTensor<STATE> TPZYCVonMisesVoigt::ComputeN(const TPZTensor<STATE> stresstensor)const
{
    STATE j2 =stresstensor.J2();
    STATE temp=sqrt(3.)/(2.* sqrt(j2));
    TPZTensor<STATE> S;
    stresstensor.S(S);
    S*=temp;
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

STATE TPZYCVonMisesVoigt::ComputeGamma(const TPZTensor<STATE>sig, const TPZFMatrix<STATE> elasticmat)const
{
    STATE j2 =sig.J2();
    STATE q=sqrt(3.* j2);
    STATE f=q-fSigmaY0;
    TPZTensor<STATE> Nvec = ComputeN(sig);
    //std::cout << "Nvec  = "<<Nvec <<std::endl;
    TPZFMatrix<STATE> NFmat(6,1,0.),tempsol,NFmatT,sol;
    Nvec.CopyTo(NFmat);
    elasticmat.Multiply(NFmat,tempsol);
    //tempsol.Print("tempsol");
    NFmat.Transpose(&NFmatT);
    NFmatT.Multiply(tempsol,sol);
    //sol.Print("sol");
 return f/sol(0,0);

}

