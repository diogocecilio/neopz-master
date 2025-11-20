#include "TPZYCMohrCoulombPV2.h"

TPZYCMohrCoulombPV2::TPZYCMohrCoulombPV2() : fPhi(0.), fPsi(0.), fc(0.), fER(), fEpsPlasticBar(0.) {

}

TPZYCMohrCoulombPV2::TPZYCMohrCoulombPV2(REAL Phi, REAL Psi, REAL c, TPZElasticResponse &ER) : fPhi(Phi), fPsi(Psi), fc(c), fER(ER), fEpsPlasticBar(0.) {

}

TPZYCMohrCoulombPV2::TPZYCMohrCoulombPV2(const TPZYCMohrCoulombPV2 &cp)
{
    TPZYCMohrCoulombPV2::operator=(cp);
}

TPZYCMohrCoulombPV2 & TPZYCMohrCoulombPV2::operator=(const TPZYCMohrCoulombPV2 &cp) {
    fPhi = cp.fPhi;
    fPsi = cp.fPsi;
    fc = cp.fc;
    fEpsPlasticBar = cp.fEpsPlasticBar;
    fER = cp.fER;
    return *this;
}



int TPZYCMohrCoulombPV2::ClassId() const{
    return Hash("TPZYCMohrCoulombPV2");
}

void TPZYCMohrCoulombPV2::Read(TPZStream& buf, void* context) { //ok
    buf.Read(&fPhi);
    buf.Read(&fPsi);
    buf.Read(&fc);
    buf.Read(&fEpsPlasticBar);
    fER.Read(buf, context);
}

void TPZYCMohrCoulombPV2::Write(TPZStream& buf, int withclassid) const { //ok
    buf.Write(&fPhi);
    buf.Write(&fPsi);
    buf.Write(&fc);
    buf.Write(&fEpsPlasticBar);
    fER.Write(buf, withclassid);
}
void TPZYCMohrCoulombPV2::YieldFunction(const TPZVec<STATE> &sigma, STATE kprev, TPZVec<STATE> &yield) const
{

}

int TPZYCMohrCoulombPV2::GetNYield() const
{
    return 3;
}

void TPZYCMohrCoulombPV2::SetLocalMatState ( TPZPlasticState<REAL> & state )
{

}

TPZPlasticState<REAL>  TPZYCMohrCoulombPV2::GetLocalMatState (  )
{

}

void  TPZYCMohrCoulombPV2::ChangeLocalMatParameters( TPZPlasticState<REAL> & state ,REAL factor)
{

}

TPZTensor<STATE> TPZYCMohrCoulombPV2::ComputeN(const TPZTensor<STATE> stresstensor)const
{

}

TPZFMatrix<STATE> TPZYCMohrCoulombPV2::GetNdSigma(const TPZTensor<STATE>& sigma) const
{

}

bool TPZYCMohrCoulombPV2::ComputeLambdaSigmaMainPlane(TPZManVector<STATE,3> &sigtr,STATE &alphan,TPZManVector<STATE,2> &dlambda,TPZManVector<STATE,3> &sigpr,TPZManVector<STATE,3> &epstr,TPZFNMatrix<9> &Grad3x3,STATE &alphan1)
{
    STATE G = fER.G();
    STATE K = fER.K();
    STATE s1=sigtr[0],s2=sigtr[1],s3=sigtr[2];
    STATE phi=fPhi;
    STATE c=fc;

    STATE cosphi=cos(phi);
    STATE sinphi=sin(phi);
    dlambda.Resize(1);
    dlambda[0]=(0.75*(s1 - 1.*s3 - 2.*c*cosphi + (s1 + s3)*sinphi))/(3.*G + (G + 3.*K)*pow(sinphi,2));

    sigpr={(2.*c*cosphi*(3.*G + (G + 3.*K)*sinphi) + (-1. + sinphi)*(-3.*G*(s1 + s3) + (G + 3.*K)*(s1 - 1.*s3)*sinphi))/(6.*G + 2.*(G + 3.*K)*pow(sinphi,2)),
        (6.*G*s2 + (2.*G - 3.*K)*(s1 - 1.*s3 - 2.*c*cosphi)*sinphi + (-3.*K*(s1 - 2.*s2 + s3) + 2.*G*(s1 + s2 + s3))*pow(sinphi,2))/
        (6.*G + 2.*(G + 3.*K)*pow(sinphi,2)),(2.*c*cosphi*(-3.*G + (G + 3.*K)*sinphi) - 1.*(1. + sinphi)*(-3.*G*(s1 + s3) + (G + 3.*K)*(s1 - 1.*s3)*sinphi))/
        (6.*G + 2.*(G + 3.*K)*pow(sinphi,2))};

    Grad3x3.Resize(3,3);
    Grad3x3(0,0)=((-1. + sinphi)*(-3.*G + (G + 3.*K)*sinphi))/(6.*G + 2.*(G + 3.*K)*pow(sinphi,2));
    Grad3x3(0,1)=0.;
    Grad3x3(0,2)=1/((-2.*sinphi)/(-1. + sinphi) + (6.*G)/(3.*G + (G + 3.*K)*sinphi));

    Grad3x3(1,0)=((2.*G - 3.*K)*sinphi*(1. + sinphi))/(6.*G + 2.*(G + 3.*K)*pow(sinphi,2));
    Grad3x3(1,1)=1.;
    Grad3x3(1,2)=((2.*G - 3.*K)*(-1. + sinphi)*sinphi)/(6.*G + 2.*(G + 3.*K)*pow(sinphi,2));

    Grad3x3(2,0)=1/(-2. + 2./(1. + sinphi) - (6.*G)/(-3.*G + (G + 3.*K)*sinphi));
    Grad3x3(2,1)=0.;
    Grad3x3(2,2)=((1. + sinphi)*(3.*G + (G + 3.*K)*sinphi))/(6.*G + 2.*(G + 3.*K)*pow(sinphi,2));

    epstr={-0.16666666666666666*(-2*s1 + s2 + s3)/G + (s1 + s2 + s3)/(9.*K),-0.16666666666666666*(s1 - 2*s2 + s3)/G + (s1 + s2 + s3)/(9.*K),-0.16666666666666666*(s1 + s2 - 2*s3)/G + (s1 + s2 + s3)/(9.*K)};

    alphan1=alphan+ dlambda[0] * 2. * cosphi;
    STATE sp1=sigpr[0];
    STATE sp2=sigpr[1];
    STATE sp3=sigpr[2];

    if((sp1 > sp2 || IsZero(sp1 - sp2)) && (sp2 > sp3 || IsZero(sp2 - sp3))){
        return true;
    }else{
        return false;
    }
}


bool TPZYCMohrCoulombPV2::ComputeLambdaSigmaLeft(TPZManVector<STATE,3> &sigtr,STATE &alphan,TPZManVector<STATE,2> &dlambda,TPZManVector<STATE,3> &sigpr,TPZManVector<STATE,3> &epstr,TPZFNMatrix<9> &Grad3x3,STATE &alphan1)
{
    STATE G = fER.G();
    STATE K = fER.K();
    STATE s1=sigtr[0],s2=sigtr[1],s3=sigtr[2];
    STATE phi=fPhi;
    STATE c=fc;

    STATE cosphi=cos(phi);
    STATE sinphi=sin(phi);
    dlambda.Resize(2);
    dlambda[0]=0.25*((s1 - 1.*s2)/(G + G*sinphi) + (3.*(s1 + s2 - 2.*s3 - 4.*c*cosphi + (s1 + s2 + 2.*s3)*sinphi))/(9.*G - 6.*G*sinphi + (G + 12.*K)*pow(sinphi,2)));

    dlambda[1]=0.25*((-1.*s1 + s2)/(G + G*sinphi) + (3.*(s1 + s2 - 2.*s3 - 4.*c*cosphi + (s1 + s2 + 2.*s3)*sinphi))/(9.*G - 6.*G*sinphi + (G + 12.*K)*pow(sinphi,2)));

    sigpr={(2.*c*cosphi*(3.*G - 1.*(G - 6.*K)*sinphi) + (-1. + sinphi)*(-3.*G*(s1 + s2 + s3) + (3.*K*(s1 + s2 - 2.*s3) + G*(s1 + s2 + s3))*sinphi))/(9.*G - 6.*G*sinphi + (G + 12.*K)*pow(sinphi,2)),(2.*c*cosphi*(3.*G - 1.*(G - 6.*K)*sinphi) +
        (-1. + sinphi)*(-3.*G*(s1 + s2 + s3) + (3.*K*(s1 + s2 - 2.*s3) + G*(s1 + s2 + s3))*sinphi))/(9.*G - 6.*G*sinphi + (G + 12.*K)*pow(sinphi,2)),(4.*c*G*cosphi*(-3. + sinphi) - 1.*(1. + sinphi)*(-3.*G*(s1 + s2 + s3) + (3.*K*(s1 + s2 - 2.*s3) + G*(s1 + s2 + s3))*sinphi) + 6.*c*K*sin(2.*phi))/(9.*G - 6.*G*sinphi + (G + 12.*K)*pow(sinphi,2))};


    Grad3x3(0,0)=((-1. + sinphi)*(-3.*G + (G + 3.*K)*sinphi))/(9.*G - 6.*G*sinphi + (G + 12.*K)*pow(sinphi,2));
    Grad3x3(0,1)=((-1. + sinphi)*(-3.*G + (G + 3.*K)*sinphi))/(9.*G - 6.*G*sinphi + (G + 12.*K)*pow(sinphi,2));
    Grad3x3(0,2)=((-1. + sinphi)*(-3.*G + (G - 6.*K)*sinphi))/(9.*G - 6.*G*sinphi + (G + 12.*K)*pow(sinphi,2));

    Grad3x3(1,0)=((-1. + sinphi)*(-3.*G + (G + 3.*K)*sinphi))/(9.*G - 6.*G*sinphi + (G + 12.*K)*pow(sinphi,2));
    Grad3x3(1,1)=((-1. + sinphi)*(-3.*G + (G + 3.*K)*sinphi))/(9.*G - 6.*G*sinphi + (G + 12.*K)*pow(sinphi,2));
    Grad3x3(1,2)=((-1. + sinphi)*(-3.*G + (G - 6.*K)*sinphi))/(9.*G - 6.*G*sinphi + (G + 12.*K)*pow(sinphi,2));

    Grad3x3(2,0)=(-1.*(1. + sinphi)*(-3.*G + (G + 3.*K)*sinphi))/(9.*G - 6.*G*sinphi + (G + 12.*K)*pow(sinphi,2));
    Grad3x3(2,1)=(-1.*(1. + sinphi)*(-3.*G + (G + 3.*K)*sinphi))/(9.*G - 6.*G*sinphi + (G + 12.*K)*pow(sinphi,2));
    Grad3x3(2,2)=(-1.*(1. + sinphi)*(-3.*G + (G - 6.*K)*sinphi))/(9.*G - 6.*G*sinphi + (G + 12.*K)*pow(sinphi,2));

    epstr={-0.16666666666666666*(-2*s1 + s2 + s3)/G + (s1 + s2 + s3)/(9.*K),-0.16666666666666666*(s1 - 2*s2 + s3)/G + (s1 + s2 + s3)/(9.*K),-0.16666666666666666*(s1 + s2 - 2*s3)/G + (s1 + s2 + s3)/(9.*K)};
    alphan1 = alphan+(dlambda[0] + dlambda[1]) * 2. * cosphi;

    STATE sp1=sigpr[0];
    STATE sp2=sigpr[1];
    STATE sp3=sigpr[2];

    if((sp1 > sp2 || IsZero(sp1 - sp2)) && (sp2 > sp3 || IsZero(sp2 - sp3))){
        return true;
    }else{
        return false;
    }

}

bool TPZYCMohrCoulombPV2::ComputeLambdaSigmaRigth(TPZManVector<STATE,3> &sigtr,STATE &alphan,TPZManVector<STATE,2> &dlambda,TPZManVector<STATE,3> &sigpr,TPZManVector<STATE,3> &epstr,TPZFNMatrix<9> &Grad3x3,STATE &alphan1)
{
    STATE G = fER.G();
    STATE K = fER.K();
    STATE s1=sigtr[0],s2=sigtr[1],s3=sigtr[2];
    STATE phi=fPhi;
    STATE c=fc;

    STATE cosphi=cos(phi);
    STATE sinphi=sin(phi);
    dlambda.Resize(2);
    dlambda[0]=0.25*((s2 - 1.*s3)/(G - 1.*G*sinphi) + (3.*(2.*s1 - 1.*s2 - 1.*s3 - 4.*c*cosphi + (2.*s1 + s2 + s3)*sinphi))/
    (9.*G + 6.*G*sinphi + (G + 12.*K)*pow(sinphi,2)));

    dlambda[1]=0.25*((s2 - 1.*s3)/(G*(-1. + sinphi)) + (3.*(2.*s1 - 1.*s2 - 1.*s3 - 4.*c*cosphi + (2.*s1 + s2 + s3)*sinphi))/
    (9.*G + 6.*G*sinphi + (G + 12.*K)*pow(sinphi,2)));

    sigpr={(4.*c*G*cosphi*(3. + sinphi) - 1.*(-1. + sinphi)*(3.*G*(s1 + s2 + s3) + (3.*K*(-2.*s1 + s2 + s3) + G*(s1 + s2 + s3))*sinphi) + 6.*c*K*sin(2.*phi))/
        (9.*G + 6.*G*sinphi + (G + 12.*K)*pow(sinphi,2)),(-2.*c*G*cosphi*(3. + sinphi) +
        (1. + sinphi)*(3.*G*(s1 + s2 + s3) + (3.*K*(-2.*s1 + s2 + s3) + G*(s1 + s2 + s3))*sinphi) + 6.*c*K*sin(2.*phi))/
        (9.*G + 6.*G*sinphi + (G + 12.*K)*pow(sinphi,2)),(-2.*c*G*cosphi*(3. + sinphi) +
        (1. + sinphi)*(3.*G*(s1 + s2 + s3) + (3.*K*(-2.*s1 + s2 + s3) + G*(s1 + s2 + s3))*sinphi) + 6.*c*K*sin(2.*phi))/
        (9.*G + 6.*G*sinphi + (G + 12.*K)*pow(sinphi,2))};


    Grad3x3(0,0)=(-1.*(-1. + sinphi)*(3.*G + (G - 6.*K)*sinphi))/(9.*G + 6.*G*sinphi + (G + 12.*K)*pow(sinphi,2));
    Grad3x3(0,1)=(-1.*(-1. + sinphi)*(3.*G + (G + 3.*K)*sinphi))/(9.*G + 6.*G*sinphi + (G + 12.*K)*pow(sinphi,2));
    Grad3x3(0,2)=(-1.*(-1. + sinphi)*(3.*G + (G + 3.*K)*sinphi))/(9.*G + 6.*G*sinphi + (G + 12.*K)*pow(sinphi,2));

    Grad3x3(1,0)=((1. + sinphi)*(3.*G + (G - 6.*K)*sinphi))/(9.*G + 6.*G*sinphi + (G + 12.*K)*pow(sinphi,2));
    Grad3x3(1,1)=((1. + sinphi)*(3.*G + (G + 3.*K)*sinphi))/(9.*G + 6.*G*sinphi + (G + 12.*K)*pow(sinphi,2));
    Grad3x3(1,2)=((1. + sinphi)*(3.*G + (G + 3.*K)*sinphi))/(9.*G + 6.*G*sinphi + (G + 12.*K)*pow(sinphi,2));

    Grad3x3(2,0)=((1. + sinphi)*(3.*G + (G - 6.*K)*sinphi))/(9.*G + 6.*G*sinphi + (G + 12.*K)*pow(sinphi,2));
    Grad3x3(2,1)=((1. + sinphi)*(3.*G + (G + 3.*K)*sinphi))/(9.*G + 6.*G*sinphi + (G + 12.*K)*pow(sinphi,2));
    Grad3x3(2,2)=((1. + sinphi)*(3.*G + (G + 3.*K)*sinphi))/(9.*G + 6.*G*sinphi + (G + 12.*K)*pow(sinphi,2));

    epstr={-0.16666666666666666*(-2*s1 + s2 + s3)/G + (s1 + s2 + s3)/(9.*K),-0.16666666666666666*(s1 - 2*s2 + s3)/G + (s1 + s2 + s3)/(9.*K),-0.16666666666666666*(s1 + s2 - 2*s3)/G + (s1 + s2 + s3)/(9.*K)};
    alphan1 = alphan+(dlambda[0] + dlambda[1]) * 2. * cosphi;

    STATE sp1=sigpr[0];
    STATE sp2=sigpr[1];
    STATE sp3=sigpr[2];

    if((sp1 > sp2 || IsZero(sp1 - sp2)) && (sp2 > sp3 || IsZero(sp2 - sp3))){
        return true;
    }else{
        return false;
    }

}

bool TPZYCMohrCoulombPV2::ReturnMapApex(TPZManVector<STATE,3> &sigtr,STATE &alphan,TPZManVector<STATE,2> &dlambda,TPZManVector<STATE,3> &sigpr,TPZManVector<STATE,3> &epstr,TPZFNMatrix<9> &Grad3x3,STATE &alphan1)
{
    STATE c=fc;
    STATE phi=fPhi;
    STATE cotphi=1/tan(phi);
    STATE s1=sigtr[0],s2=sigtr[1],s3=sigtr[2];
    sigpr={c*cotphi,c*cotphi,c*cotphi};

    STATE G = fER.G();

    STATE K = fER.K();

    Grad3x3.Zero();
    Grad3x3+=-10e-12;

    epstr={-0.16666666666666666*(-2*s1 + s2 + s3)/G + (s1 + s2 + s3)/(9.*K),-0.16666666666666666*(s1 - 2*s2 + s3)/G + (s1 + s2 + s3)/(9.*K),-0.16666666666666666*(s1 + s2 - 2*s3)/G + (s1 + s2 + s3)/(9.*K)};
    alphan1=alphan;

    return true;
}


STATE TPZYCMohrCoulombPV2::ProjectSigma(const TPZTensor<STATE> & sigmatr,  TPZTensor<STATE> & sigmaproj, TPZElasticResponse &ER,STATE &havarn,STATE &havarn1,int & m_type)
{

}
STATE TPZYCMohrCoulombPV2::ProjectSigma(TPZManVector<STATE,3> &sigtr,STATE &alphan,TPZManVector<STATE,2> &dlambda,TPZManVector<STATE,3> &sigpr,TPZManVector<STATE,3> &epstr,TPZFNMatrix<9> &Grad3x3,STATE &alphan1,int & m_type)
{
    STATE c=fc;
    STATE phi=fPhi;
    STATE sinphi=sin(phi);
    STATE cosphi=cos(phi);


    STATE phimain=(sigtr[0]-sigtr[2])+(sigtr[0]+sigtr[2])*sin(phi)-2.* c* cos(phi);
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
    STATE valcheck= (1-sinphi)*sigtr[0]-2. * sigtr[1]+(1+sinphi)*sigtr[2];

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
    check= ReturnMapApex(sigtr,alphan,dlambda,sigpr,epstr,Grad3x3,alphan1);
    return check;

}



































