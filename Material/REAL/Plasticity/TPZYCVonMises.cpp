// $Id: TPZYCVonMises.cpp,v 1.3 2008-03-08 03:12:52 erick Exp $

#include "TPZYCVonMises.h"


TPZYCVonMises::TPZYCVonMises() 
{

}

TPZYCVonMises::~TPZYCVonMises()
{

}

TPZYCVonMises::TPZYCVonMises(const TPZYCVonMises & copy)
{
	
}


void TPZYCVonMises::SetElasticResponse(TPZElasticResponse &ER) {

}

TPZElasticResponse TPZYCVonMises::GetElasticResponse() const {
	return fElasticResponse;
}

void TPZYCVonMises::Read(TPZStream &buf) {

}

void TPZYCVonMises::Write(TPZStream &buf) const {

}

TPZElasticResponse TPZYCVonMises::GetElasticResponse() {
	return fElasticResponse;
}



void TPZYCVonMises::FromHWCylToPrincipal(const TPZVec<STATE> &HWCylCoords, TPZVec<STATE> &HWCart) {

	HWCart[0] = (1. / sqrt(3.)) * HWCylCoords[0] + sqrt(2. / 3.) * HWCylCoords[1] * cos(HWCylCoords[2]);
	HWCart[1] = (1. / sqrt(3.)) * HWCylCoords[0] + sqrt(2. / 3.) * HWCylCoords[1] * cos(HWCylCoords[2] - (2. * M_PI / 3.));
	HWCart[2] = (1. / sqrt(3.)) * HWCylCoords[0] + sqrt(2. / 3.) * HWCylCoords[1] * cos(HWCylCoords[2] + (2. * M_PI / 3.));
}

void TPZYCVonMises::FromHWCylToHWCart(const TPZVec<STATE> &HWCylCoords, TPZVec<STATE> &cart) {
	cart[0] = HWCylCoords[0];
	cart[1] = HWCylCoords[1] * cos(HWCylCoords[2]);
	cart[2] = HWCylCoords[1] * sin(HWCylCoords[2]);

}

void TPZYCVonMises::FromPrincipalToHWCart(const TPZVec<STATE> &PrincipalCoords, TPZVec<STATE> &HWCart) {
	TPZFNMatrix<9, STATE> Rot(3, 3, 0.), temp(3, 1, 0.), cart(3, 1, 0.);
	HWCart.Resize(3, 0.);
	temp(0, 0) = PrincipalCoords[0];
	temp(1, 0) = PrincipalCoords[1];
	temp(2, 0) = PrincipalCoords[2];
	GetRotMatrix(Rot);
	Rot.Multiply(temp, cart);
	HWCart[0] = cart(0, 0);
	HWCart[1] = cart(1, 0);
	HWCart[2] = cart(2, 0);

}

void TPZYCVonMises::FromPrincipalToHWCyl(const TPZVec<STATE> &PrincipalCoords, TPZVec<STATE> &HWCyl) {
	TPZFNMatrix<9, STATE> Rot(3, 3, 0.), temp(3, 1, 0.), cart(3, 1, 0.);
	temp(0, 0) = PrincipalCoords[0];
	temp(1, 0) = PrincipalCoords[1];
	temp(2, 0) = PrincipalCoords[2];
	GetRotMatrix(Rot);
	Rot.Multiply(temp, cart);
	HWCyl[0] = cart(0, 0);
	HWCyl[1] = sqrt(cart(1, 0) * cart(1, 0) + cart(2, 0) * cart(2, 0));
	HWCyl[2] = atan2(cart(2, 0), cart(1, 0));
	//    HWCyl[2]=atan(cart(2,0)/cart(1,0));
}


void TPZYCVonMises::GetRotMatrix(TPZFMatrix<STATE> &Rot) {
	Rot.Resize(3, 3);
	Rot(0, 0) = 1. / sqrt(3.);
	Rot(0, 1) = 1. / sqrt(3.);
	Rot(0, 2) = 1. / sqrt(3.);
	Rot(1, 0) = sqrt(2. / 3.);
	Rot(1, 1) = -1. / sqrt(6.);
	Rot(1, 2) = -1. / sqrt(6.);
	Rot(2, 0) = 0;
	Rot(2, 1) = 1. / sqrt(2.);
	Rot(2, 2) = -1. / sqrt(2.);
}

void TPZYCVonMises::ComputeI1(TPZVec<STATE> stress, STATE &I1)const {
	STATE sig1, sig2, sig3;
	sig1 = stress[0];
	sig2 = stress[1];
	sig3 = stress[2];
	I1 = sig1 + sig2 + sig3;

}

void TPZYCVonMises::ComputeJ2(TPZVec<STATE> stress, STATE &J2)const {
	STATE sig1, sig2, sig3;
	sig1 = stress[0];
	sig2 = stress[1];
	sig3 = stress[2];
	J2 = (2. * sig1 * sig2 + pow(sig1 + (-sig1 - sig2 - sig3) / 3., 2.) +
		pow(sig2 + (-sig1 - sig2 - sig3) / 3., 2.) + 2 * sig1 * sig3 + 2. * sig2 * sig3 +
		pow((-sig1 - sig2 - sig3) / 3. + sig3, 2.)) / 2.;
}

void TPZYCVonMises::ApplyStrainComputeElasticStress(TPZVec<STATE> &strain, TPZVec<STATE> &stress)const {
	STATE sig1, sig2, sig3, s1, s2, s3;
	sig1 = strain[0];
	sig2 = strain[1];
	sig3 = strain[2];

	s1 = sig1 - (1. / 3.)*(sig1 + sig2 + sig3);
	s2 = sig2 - (1. / 3.)*(sig1 + sig2 + sig3);
	s3 = sig3 - (1. / 3.)*(sig1 + sig2 + sig3);

	stress[0] = s1 * (2 * fG) + fK * (sig1 + sig2 + sig3);
	stress[1] = s2 * (2 * fG) + fK * (sig1 + sig2 + sig3);
	stress[2] = s3 * (2 * fG) + fK * (sig1 + sig2 + sig3);
}

void TPZYCVonMises::ApplyStressComputeElasticStrain(TPZVec<STATE> &stress, TPZVec<STATE> &strain)const {
	STATE sig1, sig2, sig3, s1, s2, s3;
	sig1 = stress[0];
	sig2 = stress[1];
	sig3 = stress[2];

	s1 = sig1 - (1. / 3.)*(sig1 + sig2 + sig3);
	s2 = sig2 - (1. / 3.)*(sig1 + sig2 + sig3);
	s3 = sig3 - (1. / 3.)*(sig1 + sig2 + sig3);

	strain[0] = s1 / (2. * fG) + (sig1 + sig2 + sig3) / (9. * fK);
	strain[1] = s2 / (2. * fG) + (sig1 + sig2 + sig3) / (9. * fK);
	strain[2] = s3 / (2. * fG) + (sig1 + sig2 + sig3) / (9. * fK);

}

/**
* Imposes the specified strain tensor and returns the correspondent stress state.
*
* @param[in] epsTotal Imposed total strain tensor
* @param[out] sigma Resultant stress
*/
void TPZYCVonMises::ApplyStrainComputeSigma(TPZVec<STATE> &epst, TPZVec<STATE> &epsp, STATE & kprev, TPZVec<STATE> &epspnext, TPZVec<STATE> &stressnext, STATE & knext) const {


}

void TPZYCVonMises::ProjectSigma(const TPZVec<STATE> &sigtrial, STATE kprev, TPZVec<STATE> &sigproj, STATE &kproj) const {
	

}


void TPZYCVonMises::ProjectSigmaDep(const TPZVec<STATE> &sigtrial, STATE kprev, TPZVec<STATE> &sigproj, STATE &kproj, TPZFMatrix<STATE> &GradSigma) const {
	
}
