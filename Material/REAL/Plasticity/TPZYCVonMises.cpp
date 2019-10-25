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
	TPZYCVonMises::operator=(copy);
}


void TPZYCVonMises::SetElasticResponse(TPZElasticResponse &ER) {
	fElasticResponse = ER;
}

TPZElasticResponse TPZYCVonMises::GetElasticResponse() const {
	return fElasticResponse;
}

void TPZYCVonMises::Read(TPZStream &buf) {

}

void TPZYCVonMises::Write(TPZStream &buf) const {

}

void TPZYCVonMises::Phi(TPZVec<STATE> sigvec, STATE alpha, TPZVec<STATE> &phi)const
{
	DebugStop();
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

void TPZYCVonMises::GetPMatrix(TPZFMatrix<STATE> &P) {
	P.Resize(6, 6);
	P(0, 0) = 2. / 3.; P(0, 1) = -1. / 3; P(0, 2) = -1. / 3; P(0, 3) = 0.; P(0, 4) = 0.; P(0, 5) = 0.;
	P(1, 0) = -1. / 3.; P(1, 1) = 2. / 3; P(1, 2) = -1. / 3; P(1, 3) = 0.; P(1, 4) = 0.; P(1, 5) = 0.;
	P(2, 0) = -1. / 3.; P(2, 1) = -1. / 3; P(2, 2) = 2. / 3; P(2, 3) = 0.; P(2, 4) = 0.; P(2, 5) = 0.;
	P(3, 0) = 0.; P(3, 1) = 0.; P(3, 2) = 0.; P(3, 3) = 2.; P(3, 4) = 0.; P(3, 5) = 0.;
	P(4, 0) = 0.; P(4, 1) = 0.; P(4, 2) = 0.; P(4, 3) = 0.; P(4, 4) = 2.; P(4, 5) = 0.;
	P(5, 0) = 0.; P(5, 1) = 0.; P(5, 2) = 0.; P(5, 3) = 0.; P(5, 4) = 0.; P(5, 5) = 2.;

}

//void TPZYCVonMises::GetPMatrix(TPZFMatrix<STATE> &P) {
//	P.Resize(6, 6);
//	P(0, 0) = 2. / 3.; P(0, 1) = 0.; P(0, 2) = 0.; P(0, 3) = -1. / 3.; P(0, 4) = 0.; P(0, 5) = -1. / 3.;
//	P(1, 0) = 0.; P(1, 1) = 2.; P(1, 2) = 0.; P(1 ,3) = 0.; P(1, 4) = 0.; P(1, 5) =0.;//XY
//	P(2, 0) = 0.; P(2, 1) = 0.; P(2, 2) = 2.; P(2, 3) = 0.; P(2, 4) = 0.; P(2, 5) = 0.;//XZ
//	P(3, 0) = -1. / 3.; P(3, 1) = 0.; P(3, 2) = 0.; P(3, 3) = 2./3.; P(3, 4) =0.; P(3, 5) = -1. / 3.;
//	P(4, 0) = 0.; P(4, 1) = 0.; P(4, 2) = 2.; P(4, 3) = 0.; P(4, 4) = 0.; P(4, 5) = 0.;//YZ
//	P(5, 0) = -1. / 3.; P(5, 1) = 0.; P(5, 2) = 0; P(5, 3) = -1. / 3.; P(5, 4) = 0.; P(5, 5) = 2. / 3.;
//	
//	//#define _XX_ 0
//	//#define _XY_ 1
//	//#define _XZ_ 2
//	//#define _YY_ 3
//	//#define _YZ_ 4
//	//#define _ZZ_ 5
//}

//void TPZYCVonMises::GetCMatrix(TPZFMatrix<STATE> &C) {
//	C.Resize(6, 6);
//	STATE val1 = 4. * fG / 3. + fK;
//	STATE val2 = -2. * fG / 3. + fK;
//	STATE val3 = fG;
//	C(0, 0) = val1; C(0, 1) =   0.; C(0, 2) =  0.; C(0, 3) =val2; C(0, 4) =  0.; C(0, 5) = val2;//XX
//	C(1, 0) =   0.; C(1, 1) = val3; C(1, 2) =  0.; C(1, 3) =   0.; C(1, 4) = 0.; C(1, 5) =   0.;//XY
//	C(2, 0) =   0.; C(2, 1) =   0.; C(2, 2) =val3; C(2, 3) =   0.; C(2, 4) = 0.; C(2, 5) =   0.;//XZ
//	C(3, 0) = val2; C(3, 1) =   0.; C(3, 2) =  0.; C(3, 3) = val1; C(3, 4) = 0.; C(3, 5) = val2;//YY
//	C(4, 0) =   0.; C(4, 1) =   0.; C(4, 2) =  0.; C(4, 3) =   0.; C(4, 4) = val3; C(4, 5) =   0.;//YZ
//	C(5, 0) = val2; C(5, 1) =   0.; C(5, 2) =  0.; C(5, 3) = val2; C(5, 4) = 0.; C(5, 5) = val3;//ZZ
//}

void TPZYCVonMises::GetCMatrix(TPZFMatrix<STATE> &C) {
	C.Resize(6, 6);
	STATE val1 = 4. * fG / 3. + fK;
	STATE val2 = -2. * fG / 3. + fK;
	STATE val3 = fG;
	C(0, 0) = val1; C(0, 1) = val2; C(0, 2) = val2; C(0, 3) = 0.; C(0, 4) = 0.; C(0, 5) = 0.;//XX
	C(1, 0) = val2; C(1, 1) = val1; C(1, 2) = val2; C(1, 3) = 0.; C(1, 4) = 0.; C(1, 5) = 0.;//XY
	C(2, 0) = val2; C(2, 1) = val2; C(2, 2) = val1; C(2, 3) = 0.; C(2, 4) = 0.; C(2, 5) = 0.;//XZ
	C(3, 0) = 0.; C(3, 1) = 0.; C(3, 2) = 0.; C(3, 3) = val3; C(3, 4) = 0.; C(3, 5) = 0.;//YY
	C(4, 0) = 0.; C(4, 1) = 0.; C(4, 2) = 0.; C(4, 3) = 0.; C(4, 4) = val3; C(4, 5) = 0.;//YZ
	C(5, 0) = 0.; C(5, 1) = 0.; C(5, 2) = 0.; C(5, 3) = 0.; C(5, 4) = 0.; C(5, 5) = val3;//ZZ
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
	//J2 = (2. * sig1 * sig2 + pow(sig1 + (-sig1 - sig2 - sig3) / 3., 2.) +
	//	pow(sig2 + (-sig1 - sig2 - sig3) / 3., 2.) + 2 * sig1 * sig3 + 2. * sig2 * sig3 +
	//	pow((-sig1 - sig2 - sig3) / 3. + sig3, 2.)) / 2.;
	J2 = (pow(sig1, 2) + pow(sig2, 2) - sig2*sig3 + pow(sig3, 2) - sig1*(sig2 + sig3)) / 3.;
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
	
	STATE I1, J2;

	ComputeI1(sigtrial, I1);

	ComputeJ2(sigtrial, J2);

	STATE xitrial = I1 / sqrt(3.);

	STATE rhotrial = sqrt(2.* J2);

	STATE xisol = xitrial;

	STATE rhosol = sqrt(2. / 3.)*fsigy;

	STATE betasol;

	STATE denom = (-2.*sigtrial[0] + sigtrial[1] + sigtrial[2]);

	if (fabs(denom) < 1.e-6) {
		betasol = 0.;
	}
	else {
		betasol = atan((sqrt(3.)*(-sigtrial[1] + sigtrial[2])) / denom);
	}

	TPZVec<STATE> sigprojcyl(3, 0.);

	sigprojcyl[0] = xisol;

	sigprojcyl[1] = rhosol;

	sigprojcyl[2] = betasol;

	sigproj.resize(3);

	FromHWCylToPrincipal(sigprojcyl, sigproj);
}


void TPZYCVonMises::ProjectSigmaDep(const TPZVec<STATE> &sigtrial, STATE kprev, TPZVec<STATE> &sigproj, STATE &kproj, TPZFMatrix<STATE> &GradSigma) const {

	STATE I1, J2;

	ComputeI1(sigtrial, I1);

	ComputeJ2(sigtrial, J2);

	STATE xitrial = I1 / sqrt(3.);

	STATE rhotrial = sqrt(2.* J2);

	STATE xisol = xitrial;

	STATE rhosol = sqrt(2. / 3.)*fsigy;

	STATE betasol;

	STATE denom = (-2.*sigtrial[0] + sigtrial[1] + sigtrial[2]);

	if (fabs(denom) < 1.e-6){
		 betasol = 0.;
	}
	else {
		 betasol = atan((sqrt(3.)*(-sigtrial[1] + sigtrial[2])) / denom);
	}

	TPZVec<STATE> sigprojcyl(3,0.);

	sigprojcyl[0] = xisol;

	sigprojcyl[1] = rhosol;

	sigprojcyl[2] = betasol;

	sigproj.resize(3);

	FromHWCylToPrincipal(sigprojcyl, sigproj);
}

void TPZYCVonMises::N(TPZTensor<STATE> sigma, TPZTensor<STATE> &asol)
{
	sigma.S(asol);
	STATE J2 = sigma.J2();
	if (fabs(J2) < 1.e-6) {
		asol = 0.;
	}
	else {
		asol *= sqrt(3.) / (2.*sqrt(J2));
	}
	
	asol.XY() *= 2.; asol.XZ() *= 2.; asol.YZ() *= 2.;
}

void TPZYCVonMises::dadsig(TPZTensor<STATE> sigma, TPZFMatrix<STATE> &dadsigmat)
{
	STATE J2 = sigma.J2();

	TPZFMatrix<STATE> P,temp;

	TPZTensor<STATE> Sdev;

	sigma.S(Sdev);

	Sdev.XY() *= 2.; Sdev.XZ() *= 2.; Sdev.YZ() *= 2.;

	temp = Sdev.ProdT2(Sdev);

	GetPMatrix(P);

	P *= sqrt(3.) / (2.*sqrt(J2));

	temp *= sqrt(3.) / (4.*pow(J2, 3. / 2.));

	P -= temp;

	dadsigmat = P;

}

void TPZYCVonMises::ComputeDep(TPZTensor<STATE> sigma, TPZTensor<STATE> epsTr, TPZTensor<STATE> epsElaNp1, TPZFMatrix<REAL> &Dep)
{

	TPZFNMatrix<36> dSigDe(6, 6, 0.);

	TPZTensor<REAL> asol, diff;

	N(sigma, asol);

	diff = epsTr;

	diff -= epsElaNp1;

	STATE norm = diff.Norm();

	STATE gamma;
	if (fabs(asol.Norm()) < 1.e-6) {
		 gamma = 0.;
	}
	else {
		 gamma = norm / asol.Norm();
	}



	TPZFMatrix<STATE> dadsigmat,Q,C,I(6,6,0.),temp,invQ,R;

	dadsig(sigma, dadsigmat);

	GetCMatrix(C);

	C.Multiply(dadsigmat,temp);

	temp *= gamma;

	I.Identity();

	I += temp;

	Q = I;

	Q.Inverse(invQ, ELU);

	invQ.Multiply(C, R);

	TPZFMatrix<STATE> asolcopy = asol.FromTensorToStandardOrder();

	temp.Zero();

	R.Multiply(asolcopy, temp);

	STATE tempscalar = 1./Dot(asolcopy, temp);

	TPZFMatrix<STATE> T;

	TPZTensor<STATE> temptensor;

	temptensor = temptensor.FromStandardToTensor(temp);

	T = temptensor.ProdT2(temptensor);

	T *= tempscalar;

	dSigDe = R;

	dSigDe -= T;

	TransForm(dSigDe, Dep);

	//std::cout << " - - C - - " << std::endl;
	//C.Print(std::cout);
	//std::cout << " - - Q - - " << std::endl;
	//Q.Print(std::cout);
	//std::cout << " - - Q^-1 - - " << std::endl;
	//invQ.Print(std::cout);
	//std::cout << " - - R - - " << std::endl;
	//R.Print(std::cout);
	//std::cout << " - - asolcopy - - " << std::endl;
	//asolcopy.Print(std::cout);
	//std::cout << " - - R.asol - - " << std::endl;
	//temp.Print(std::cout);
	//std::cout << " - - R.asol(X)R.asol - - " << std::endl;
	//T.Print(std::cout);
	//std::cout << " - - Dep - - " << std::endl;
	//Dep.Print(std::cout);


}
void TPZYCVonMises::TransForm(TPZFMatrix<REAL> &A, TPZFMatrix<REAL> &B)
{
	B.Resize(6, 6);
	int xxa=0, yya=1, zza=2,  xza=3, yza=4, xya = 5;
	int xxb=0, yyb=3, zzb=5,  xyb=1, xzb=2, yzb = 4;
	B(xxb, xxb) = A(xxa, xxa);
	B(xxb, yyb) = A(xxa, yya);
	B(xxb, zzb) = A(xxa, zza);
	B(xxb, xyb) = A(xxa, xya);
	B(xxb, xzb) = A(xxa, xza);
	B(xxb, yzb) = A(xxa, yza);

	B(yyb, xxb) = A(yya, xxa);
	B(yyb, yyb) = A(yya, yya);
	B(yyb, zzb) = A(yya, zza);
	B(yyb, xyb) = A(yya, xya);
	B(yyb, xzb) = A(yya, xza);
	B(yyb, yzb) = A(yya, yza);


	B(zzb, xxb) = A(zza, xxa);
	B(zzb, yyb) = A(zza, yya);
	B(zzb, zzb) = A(zza, zza);
	B(zzb, xyb) = A(zza, xya);
	B(zzb, xzb) = A(zza, xza);
	B(zzb, yzb) = A(zza, yza);


	B(xyb, xxb) = A(xya, xxa);
	B(xyb, yyb) = A(xya, yya);
	B(xyb, zzb) = A(xya, zza);
	B(xyb, xyb) = A(xya, xya);
	B(xyb, xzb) = A(xya, xza);
	B(xyb, yzb) = A(xya, yza);


	B(xzb, xxb) = A(xza, xxa);
	B(xzb, yyb) = A(xza, yya);
	B(xzb, zzb) = A(xza, zza);
	B(xzb, xyb) = A(xza, xya);
	B(xzb, xzb) = A(xza, xza);
	B(xzb, yzb) = A(xza, yza);

	B(yzb, xxb) = A(yza, xxa);
	B(yzb, yyb) = A(yza, yya);
	B(yzb, zzb) = A(yza, zza);
	B(yzb, xyb) = A(yza, xya);
	B(yzb, xzb) = A(yza, xza);
	B(yzb, yzb) = A(yza, yza);

}
