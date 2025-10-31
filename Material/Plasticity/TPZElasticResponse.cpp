//
//  TPZElasticResponse.h
//  pz
//
//  Created by Erick Slis Raggio Santos on 04/07/2009.
//

#include "TPZElasticResponse.h"

int TPZElasticResponse::ClassId() const{
    return Hash("TPZElasticResponse");
}

TPZElasticResponse::TPZElasticResponse() : m_lambda(0.), m_mu(0.) {
    m_epsilon_star.Zero();
    m_sigma_star.Zero();
}

TPZElasticResponse::TPZElasticResponse(const TPZElasticResponse & other) {
    m_E       = other.m_E;
    m_nu            = other.m_nu;
    m_lambda        = other.m_lambda;
    m_mu            = other.m_mu;
    m_epsilon_star      = other.m_epsilon_star;
    m_sigma_star    = other.m_sigma_star;
}

TPZElasticResponse & TPZElasticResponse::operator=(const TPZElasticResponse & other) {
    m_E       = other.m_E;
    m_nu            = other.m_nu;
    m_lambda        = other.m_lambda;
    m_mu            = other.m_mu;
    m_epsilon_star      = other.m_epsilon_star;
    return *this;
}


void TPZElasticResponse::Write(TPZStream& buf, int withclassid) const { //ok
    buf.Write(&m_lambda);
    buf.Write(&m_mu);
    m_epsilon_star.Write(buf,withclassid);
    m_sigma_star.Write(buf,withclassid);
}

void TPZElasticResponse::Read(TPZStream& buf, void* context) { //ok
    buf.Read(&m_lambda);
    buf.Read(&m_mu);
    m_epsilon_star.Read(buf, context);
    m_sigma_star.Read(buf, context);
}


const char * TPZElasticResponse::Name() const {
    return "TPZElasticResponse";
}

void TPZElasticResponse::Print(std::ostream & out) const {
    out << this->Name();
    out << "\n Young = " << E();
    out << "\n Poisson = " << Poisson();
    out << "\n m_lambda = " << m_lambda;
    out << "\n m_mu = " << m_mu;
    m_epsilon_star.Print(out);
    m_sigma_star.Print(out);
}

void TPZElasticResponse::De(TPZFMatrix<STATE> & DeMat) const {
    DeMat.Redim(6,6);
    DeMat.Zero();

    const STATE nu = Poisson();
    const STATE E  = this->E();
    const STATE factor = E/((1.+nu)*(1.-2.*nu));

    // parte normal
    DeMat(_XX_, _XX_) = 1.-nu;  DeMat(_XX_, _YY_) = nu;     DeMat(_XX_, _ZZ_) = nu;
    DeMat(_YY_, _XX_) = nu;     DeMat(_YY_, _YY_) = 1.-nu;  DeMat(_YY_, _ZZ_) = nu;
    DeMat(_ZZ_, _XX_) = nu;     DeMat(_ZZ_, _YY_) = nu;     DeMat(_ZZ_, _ZZ_) = 1.-nu;


    DeMat(_XY_, _XY_) = (1.-2.*nu)/2.;
    DeMat(_XZ_, _XZ_) = (1.-2.*nu)/2.;
    DeMat(_YZ_, _YZ_) = (1.-2.*nu)/2.;

    DeMat *= factor;
}


void TPZElasticResponse::InverseDe(TPZFMatrix<STATE> & DeMat) const {
    DeMat.Redim(6,6);
    DeMat.Zero();

    const STATE nu = Poisson();
    const STATE E  = this->E();
    const STATE factor = 1./E;

    // parte normal
    DeMat(_XX_, _XX_) = 1.;  DeMat(_XX_, _YY_) = -nu;     DeMat(_XX_, _ZZ_) = -nu;
    DeMat(_YY_, _XX_) = -nu;     DeMat(_YY_, _YY_) = 1.;  DeMat(_YY_, _ZZ_) = -nu;
    DeMat(_ZZ_, _XX_) = -nu;     DeMat(_ZZ_, _YY_) = -nu;     DeMat(_ZZ_, _ZZ_) = 1.;


    DeMat(_XY_, _XY_) =2.*(1.+nu);
    DeMat(_XZ_, _XZ_) =2.*(1.+nu);
    DeMat(_YZ_, _YZ_) =2.*(1.+nu);

    DeMat *= factor;
}
//template<class T>
void TPZElasticResponse::ComputeStrain(const TPZTensor<STATE> & sigma, TPZTensor<STATE> & epsilon) const
{

    TPZFMatrix<STATE> InvCmat,cpsigma(6,1,0.),temp;

    for(int i=0;i<6;i++)cpsigma(i,0)=sigma[i];

    InverseDe(InvCmat);

    InvCmat.Multiply(cpsigma,temp);

    epsilon.CopyFrom(temp);

}


void TPZElasticResponse::ComputeStress(const TPZTensor<STATE> & epsilon, TPZTensor<STATE> & sigma) const {

    TPZFMatrix<STATE> Cmat,cpeps(6,1,0.),temp;

    for(int i=0;i<6;i++)cpeps(i,0)=epsilon[i];

    De(Cmat);

    Cmat.Multiply(cpeps,temp);

    sigma.CopyFrom(temp);


}
void TPZElasticResponse::SetEngineeringData(REAL Eyoung, REAL Poisson) {
    m_E = Eyoung;
    m_nu = Poisson;
    SetLameData();
}

void TPZElasticResponse::SetLameData() {
     m_lambda = m_nu * m_E / ((1. + m_nu)*(1. - 2. * m_nu));
     m_mu = m_E / (2. * (1. + m_nu));
}

REAL TPZElasticResponse::Lambda() const {
    return m_lambda;
}

REAL TPZElasticResponse::K() const {
    return m_lambda + 2. * m_mu / 3.;
}

REAL TPZElasticResponse::Mu() const {
    return m_mu;
}

REAL TPZElasticResponse::G() const {
    return Mu();
}

REAL TPZElasticResponse::E() const {
    REAL E = m_mu * (3. * m_lambda + 2. * m_mu) / (m_lambda + m_mu);
    return E;
}

REAL TPZElasticResponse::Poisson() const {
    REAL poisson = m_lambda / (2. * (m_lambda + m_mu));
    return poisson;
}

void TPZElasticResponse::SetReferenceStrainData(TPZTensor<REAL> & eps_star){
    m_epsilon_star = eps_star;
}

TPZTensor<REAL> & TPZElasticResponse::ReferenceStrainData(){
    return m_epsilon_star;
}

void TPZElasticResponse::SetReferenceStressData(TPZTensor<REAL> & sigma_star){
    m_sigma_star = sigma_star;
}

TPZTensor<REAL> & TPZElasticResponse::ReferenceStressData(){
    return m_sigma_star;
}

