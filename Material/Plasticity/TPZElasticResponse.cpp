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
    m_lambda        = other.m_lambda;
    m_mu            = other.m_mu;
    m_epsilon_star      = other.m_epsilon_star;
    m_sigma_star    = other.m_sigma_star;
}

TPZElasticResponse & TPZElasticResponse::operator=(const TPZElasticResponse & other) {
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

void TPZElasticResponse::De(TPZFMatrix<REAL> & De) {
    REAL Mu2 = 2 * m_mu;
    
    De.Redim(6,6);
    De.Zero();
    
    De(_XX_, _XX_) += m_lambda;
    De(_XX_, _YY_) += m_lambda;
    De(_XX_, _ZZ_) += m_lambda;
    De(_YY_, _XX_) += m_lambda;
    De(_YY_, _YY_) += m_lambda;
    De(_YY_, _ZZ_) += m_lambda;
    De(_ZZ_, _XX_) += m_lambda;
    De(_ZZ_, _YY_) += m_lambda;
    De(_ZZ_, _ZZ_) += m_lambda;
    
    int i;
    for (i = 0; i < 6; i++)De(i, i) += Mu2;
}
void TPZElasticResponse::CMatrix(TPZFMatrix<REAL> & Cmat)
{
    Cmat.Redim(6,6);
    Cmat.Zero();

    const STATE a = (4.0*G())/3.0 + K();
    const STATE b = (-2.0*G())/3.0 + K();

    // parte volumétrica + deviadora (diagonal 3x3 e simétricos)
    Cmat(0,0) = a;  Cmat(0,1) = b;  Cmat(0,2) = b;
    Cmat(1,0) = b;  Cmat(1,1) = a;  Cmat(1,2) = b;
    Cmat(2,0) = b;  Cmat(2,1) = b;  Cmat(2,2) = a;

    // cisalhamentos (engenharia)
    Cmat(3,3) = G();  // xy
    Cmat(4,4) = G();  // yz
    Cmat(5,5) = G();  // zx
}

void TPZElasticResponse::InvCMatrix(TPZFMatrix<STATE>& S)
{
    S.Redim(6,6);
    S.Zero();

    const STATE c = (G() + 3.0*K()) / (9.0 * G() * K());
    const STATE d = (-1.0 / (6.0*G())) + (1.0 / (9.0*K()));

    // bloco 3x3 superior-esquerdo (simétrico)
    S(0,0) = c;  S(0,1) = d;  S(0,2) = d;
    S(1,0) = d;  S(1,1) = c;  S(1,2) = d;
    S(2,0) = d;  S(2,1) = d;  S(2,2) = c;

    // cisalhamentos (engenharia)
    S(3,3) = 1.0/G();  // xy
    S(4,4) = 1.0/G();  // yz
    S(5,5) = 1.0/G();  // zx
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

