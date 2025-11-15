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
: fSigmaY0(0.0), fH0(0.0),fER()
{

}


TPZYCVonMisesVoigt::TPZYCVonMisesVoigt(const TPZYCVonMisesVoigt& cp)
: TPZPlasticCriterion(cp)
, fSigmaY0(cp.fSigmaY0)
, fH0(cp.fH0),fER(cp.fER)

{

}

void TPZYCVonMisesVoigt::SetUp(STATE sigmaY0, STATE Hiso,TPZElasticResponse &ER)
{
    fSigmaY0 = sigmaY0;
    fH0=Hiso;
    fER=ER;
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

STATE TPZYCVonMisesVoigt::ProjectSigma(TPZManVector<STATE,3> &sigtr,STATE &alphan,TPZManVector<STATE,2> &dlambda,TPZManVector<STATE,3> &sigpr,TPZManVector<STATE,3> &epstr,TPZFNMatrix<9> &Grad3x3,STATE &alphan1,int & m_type)
{
    STATE G = fER.G();
    STATE K = fER.K();
    fHard=alphan;
    STATE s1=sigtr[0],s2=sigtr[1],s3=sigtr[2];
    STATE j2 =1./3. *(s1*s1 + s2*s2 - s2 *s3 + s3*s3 - s1 *(s2 + s3));
    STATE q=sqrt(3.* j2);
    STATE sigy  = SigmaY(alphan);
    // STATE sigy  = fSigmaY0;
    STATE f=q-sigy;

    if(f<0.)
    {
        sigpr=sigtr;
        //sigmaproj+=0.;
        alphan1=alphan;
        m_type=0;//elastic
        return 0.;

    }else
    {
        STATE H=DSigmaYDepsbar(alphan);

        STATE gamma =  (q-sigy)/(3*G+H);
        dlambda.Resize(1);
        dlambda[0]=gamma;
        alphan1=fHard+gamma;
        STATE sigyn1  = SigmaY(alphan);


        STATE betaproj=atan((sqrt(3.)* (-s2+s3))/(-2.* s1+s2+s3));

        //STATE betaproj = (1.0/3.0) *atan2( sqrt(3.0)*(sig2 - sig3), (2.0*sig1 - sig2 - sig3) );

        STATE I13=(s1+s2+s3)/3.;

        STATE s2igy3=2*sigyn1/3.;

        TPZManVector<STATE,3> HWCart={I13+s2igy3*cos(betaproj),I13+s2igy3*cos(betaproj-2*M_PI/3.),I13+s2igy3*cos(betaproj+2*M_PI/3.)};

        sigpr=HWCart;
        Grad3x3.Zero();
        STATE tempscal=2*sigyn1/(3*sqrt(3*j2));
        TPZManVector<STATE,3> a={sin(betaproj),sin(betaproj-2.*M_PI/3.),sin(betaproj+2.*M_PI/3.)};
        int sza=a.size();
        for(int i=0;i<sza;i++)
        {
            for(int j=0;j<sza;j++)
            {
                Grad3x3(i,j)= a[i]*a[j]*tempscal;
            }

        }
        Grad3x3+=1./3.;

        epstr={-0.16666666666666666*(-2*s1 + s2 + s3)/G + (s1 + s2 + s3)/(9.*K),-0.16666666666666666*(s1 - 2*s2 + s3)/G + (s1 + s2 + s3)/(9.*K),-0.16666666666666666*(s1 + s2 - 2*s3)/G + (s1 + s2 + s3)/(9.*K)};

        alphan1 = alphan+dlambda[0];
        m_type=1;//plastic

        return gamma;
    }
}



