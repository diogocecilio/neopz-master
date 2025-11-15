#include "TPZYCTrescaVoigt.h"
#include "pzerror.h"
#include "pzlog.h"
#include "pzvec.h"
#include "pzmanvector.h"
#include "TPZStream.h"
#include <cmath>
#include <iomanip>
#include "TPZKrylovEigenSolver.h"
#include "TPZEigenSolver.h"

TPZYCTrescaVoigt::TPZYCTrescaVoigt()
: fSigmaY0(0.0), fH0(0.0),fER()
{

}


TPZYCTrescaVoigt::TPZYCTrescaVoigt(const TPZYCTrescaVoigt& cp)
: TPZPlasticCriterion(cp)
, fSigmaY0(cp.fSigmaY0)
, fH0(cp.fH0),fER(cp.fER)

{

}

void TPZYCTrescaVoigt::SetUp(STATE sigmaY0, STATE Hiso,TPZElasticResponse &ER)
{
    fSigmaY0 = sigmaY0;
    fH0=Hiso;
    fER=ER;
}

void TPZYCTrescaVoigt::SetLocalMatState(TPZPlasticState<REAL> & /*state*/)
{
     DebugStop();
}

TPZPlasticState<REAL> TPZYCTrescaVoigt::GetLocalMatState()
{
    TPZPlasticState<REAL> st;
    return st;
}


void TPZYCTrescaVoigt::ChangeLocalMatParameters(TPZPlasticState<REAL> & /*state*/, REAL /*factor*/)
{
}

int TPZYCTrescaVoigt::ClassId() const
{
        return Hash("TPZYCVonMisesVoigt") ;
}

void TPZYCTrescaVoigt::Read(TPZStream& buf, void* /*context*/)
{
    buf.Read(&fSigmaY0,1);
}

void TPZYCTrescaVoigt::Write(TPZStream& buf, int /*withclassid*/) const
{
    buf.Write(&fSigmaY0,1);
}

void TPZYCTrescaVoigt::Print(std::ostream &out) const
{
    out << "----- TPZYCTrescaVoigt -----\n";
    out << "Yield stress (sigma_y): " << std::setprecision(12) << fSigmaY0 << "\n";
    out << "NYield = " << NYield << "\n";
    out << "---------------------------\n";
}


void TPZYCTrescaVoigt::Phi(TPZTensor<STATE>sig, STATE alpha, TPZVec<STATE> &phi) const
{
    DebugStop();
    phi.resize(1);
    STATE j2 =sig.J2();
    STATE q=sqrt(3.* j2);
    STATE f=q-SigmaY(alpha);
    phi[0]=f;
}



STATE TPZYCTrescaVoigt::ProjectSigma(const TPZTensor<STATE> & sigmatr,  TPZTensor<STATE> & sigmaproj,TPZElasticResponse &ER, STATE &havarn,STATE &havarn1, int & m_type)
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

        STATE I13=sigmatr.I1()/3.;

        STATE s2igy3=2*sigy/3.;

        TPZManVector<STATE,3> HWCart={I13+s2igy3*cos(betaproj),I13+s2igy3*cos(betaproj-2*M_PI/3.),I13+s2igy3*cos(betaproj+2*M_PI/3.)};


        sig_eigen_system.fEigenvalues=HWCart;

        TPZTensor<STATE>sigmaproj2(sig_eigen_system);
        sigmaproj=sigmaproj2;

        m_type=1;//plastic

        return gamma;
    }

}

STATE TPZYCTrescaVoigt::ProjectSigma(TPZManVector<STATE,3> &sigtr,STATE &alphan,TPZManVector<STATE,2> &dlambda,TPZManVector<STATE,3> &sigpr,TPZManVector<STATE,3> &epstr,TPZFNMatrix<9> &Grad3x3,STATE &alphan1,int & m_type)
{
    STATE G = fER.G();
    STATE K = fER.K();
    fHard=alphan;
    STATE s1=sigtr[0],s2=sigtr[1],s3=sigtr[2];
    STATE j2 =1./3. *(s1*s1 + s2*s2 - s2 *s3 + s3*s3 - s1 *(s2 + s3));
    STATE q=sqrt(3.* j2);
    STATE sigy  = SigmaY(alphan);
    // STATE sigy  = fSigmaY0;
    STATE f=(s1-s3)-sigy;

    STATE phimain=(s1-s3)-sigy;
    if(phimain<0)
    {
        //Elastico
        sigpr=sigtr;
        alphan=alphan1;
        m_type=0;
        return 0.;
    }


    m_type=1;
    bool check=ComputeLambdaSigmaMainPlane(sigtr,alphan,dlambda,sigpr,epstr,Grad3x3,alphan1);

    if(check)
    {
        return check;
    }
    STATE valcheck= sigtr[0]-2. * sigtr[1]+sigtr[2];
    //valcheck = strhw[[1]] + strhw[[3]] - 2 strhw[[2]];
    if(valcheck>0)
    {
        //std::cout<< "Rigth"<<std::endl;
        check=ComputeLambdaSigmaRigth(sigtr,alphan,dlambda,sigpr,epstr,Grad3x3,alphan1);
        //std::cout<< "Depois do Rigth"<<std::endl;
        if(check)
        {
            return check;
        }

    }else{

        check=ComputeLambdaSigmaLeft(sigtr,alphan,dlambda,sigpr,epstr,Grad3x3,alphan1);
        if(check)
        {
            return check;
        }
    }

    DebugStop();
    return check;
}

bool TPZYCTrescaVoigt::ComputeLambdaSigmaMainPlane(TPZManVector<STATE,3> &sigtr,STATE &alphan,TPZManVector<STATE,2> &dlambda,TPZManVector<STATE,3> &sigpr,TPZManVector<STATE,3> &epstr,TPZFNMatrix<9> &Grad3x3,STATE &alphan1)
{
    STATE G = fER.G();
    STATE K = fER.K();
    STATE s1=sigtr[0],s2=sigtr[1],s3=sigtr[2];
    STATE sigy=SigmaY(alphan);


    dlambda.Resize(2);
    dlambda[0]=0;
    dlambda[1]=0;

    sigpr={1./2. *(s1+s3+sigy),s2,1./2. *(s1+s3-sigy)};


    Grad3x3(0,0)=0.5;
    Grad3x3(0,1)=0;
    Grad3x3(0,2)=0.5;

    Grad3x3(1,0)=0;
    Grad3x3(1,1)=1;
    Grad3x3(1,2)=0;

    Grad3x3(2,0)=0.5;
    Grad3x3(2,1)=0;
    Grad3x3(2,2)=0.5;

    epstr={-0.16666666666666666*(-2*s1 + s2 + s3)/G + (s1 + s2 + s3)/(9.*K),-0.16666666666666666*(s1 - 2*s2 + s3)/G + (s1 + s2 + s3)/(9.*K),-0.16666666666666666*(s1 + s2 - 2*s3)/G + (s1 + s2 + s3)/(9.*K)};
    alphan1 = alphan+(dlambda[0] + dlambda[1]) ;

    STATE sp1=sigpr[0];
    STATE sp2=sigpr[1];
    STATE sp3=sigpr[2];

    if((sp1 > sp2 || IsZero(sp1 - sp2)) && (sp2 > sp3 || IsZero(sp2 - sp3))){
        return true;
    }else{
        return false;
    }
}


bool TPZYCTrescaVoigt::ComputeLambdaSigmaLeft(TPZManVector<STATE,3> &sigtr,STATE &alphan,TPZManVector<STATE,2> &dlambda,TPZManVector<STATE,3> &sigpr,TPZManVector<STATE,3> &epstr,TPZFNMatrix<9> &Grad3x3,STATE &alphan1)
{
    STATE G = fER.G();
    STATE K = fER.K();
    STATE s1=sigtr[0],s2=sigtr[1],s3=sigtr[2];
    STATE sigy=SigmaY(alphan);


    dlambda.Resize(2);
    dlambda[0]=0;
    dlambda[1]=0;

    sigpr={1./3.* (s1+s2+s3+sigy),1./3.*(s1+s2+s3+sigy),1./3.*(s1+s2+s3-2 *sigy)};


    Grad3x3(0,0)=1./3.;
    Grad3x3(0,1)=1./3.;
    Grad3x3(0,2)=1./3.;

    Grad3x3(1,0)=1./3.;
    Grad3x3(1,1)=1./3.;
    Grad3x3(1,2)=1./3.;

    Grad3x3(2,0)=1./3.;
    Grad3x3(2,1)=1./3.;
    Grad3x3(2,2)=1./3.;

    epstr={-0.16666666666666666*(-2*s1 + s2 + s3)/G + (s1 + s2 + s3)/(9.*K),-0.16666666666666666*(s1 - 2*s2 + s3)/G + (s1 + s2 + s3)/(9.*K),-0.16666666666666666*(s1 + s2 - 2*s3)/G + (s1 + s2 + s3)/(9.*K)};
    alphan1 = alphan+(dlambda[0] + dlambda[1]) ;

    STATE sp1=sigpr[0];
    STATE sp2=sigpr[1];
    STATE sp3=sigpr[2];

    if((sp1 > sp2 || IsZero(sp1 - sp2)) && (sp2 > sp3 || IsZero(sp2 - sp3))){
        return true;
    }else{
        return false;
    }
}

bool TPZYCTrescaVoigt::ComputeLambdaSigmaRigth(TPZManVector<STATE,3> &sigtr,STATE &alphan,TPZManVector<STATE,2> &dlambda,TPZManVector<STATE,3> &sigpr,TPZManVector<STATE,3> &epstr,TPZFNMatrix<9> &Grad3x3,STATE &alphan1)
{
    STATE G = fER.G();
    STATE K = fER.K();
    STATE s1=sigtr[0],s2=sigtr[1],s3=sigtr[2];
    STATE sigy=SigmaY(alphan);


    dlambda.Resize(2);
    dlambda[0]=0;
    dlambda[1]=0;

    sigpr={1./3.* (s1+s2+s3+2 *sigy),1./3.* (s1+s2+s3-sigy),1./3.* (s1+s2+s3-sigy)};


        Grad3x3(0,0)=1./3.;
        Grad3x3(0,1)=1./3.;
        Grad3x3(0,2)=1./3.;

        Grad3x3(1,0)=1./3.;
        Grad3x3(1,1)=1./3.;
        Grad3x3(1,2)=1./3.;

        Grad3x3(2,0)=1./3.;
        Grad3x3(2,1)=1./3.;
        Grad3x3(2,2)=1./3.;

        epstr={-0.16666666666666666*(-2*s1 + s2 + s3)/G + (s1 + s2 + s3)/(9.*K),-0.16666666666666666*(s1 - 2*s2 + s3)/G + (s1 + s2 + s3)/(9.*K),-0.16666666666666666*(s1 + s2 - 2*s3)/G + (s1 + s2 + s3)/(9.*K)};
        alphan1 = alphan+(dlambda[0] + dlambda[1]) ;

        STATE sp1=sigpr[0];
        STATE sp2=sigpr[1];
        STATE sp3=sigpr[2];

        if((sp1 > sp2 || IsZero(sp1 - sp2)) && (sp2 > sp3 || IsZero(sp2 - sp3))){
            return true;
        }else{
            return false;
        }

}

