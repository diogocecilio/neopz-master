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

//void TPZYCVonMises::GetPMatrix(TPZFMatrix<STATE> &P) {
//	P.Resize(6, 6);
//	P(0, 0) = 2. / 3.; P(0, 1) = -1. / 3; P(0, 2) = -1. / 3; P(0, 3) = 0.; P(0, 4) = 0.; P(0, 5) = 0.;
//	P(1, 0) = -1. / 3.; P(1, 1) = 2. / 3; P(1, 2) = -1. / 3; P(1, 3) = 0.; P(1, 4) = 0.; P(1, 5) = 0.;
//	P(2, 0) = -1. / 3.; P(2, 1) = -1. / 3; P(2, 2) = 2. / 3; P(2, 3) = 0.; P(2, 4) = 0.; P(2, 5) = 0.;
//	P(3, 0) = 0.; P(3, 1) = 0.; P(3, 2) = 0.; P(3, 3) = 2.; P(3, 4) = 0.; P(3, 5) = 0.;
//	P(4, 0) = 0.; P(4, 1) = 0.; P(4, 2) = 0.; P(4, 3) = 0.; P(4, 4) = 2.; P(4, 5) = 0.;
//	P(5, 0) = 0.; P(5, 1) = 0.; P(5, 2) = 0.; P(5, 3) = 0.; P(5, 4) = 0.; P(5, 5) = 2.;
//
//}

void TPZYCVonMises::GetPMatrix(TPZFMatrix<STATE> &P) {
	P.Resize(6, 6);
	P(0, 0) = 2. / 3.; P(0, 1) = 0.; P(0, 2) = 0.; P(0, 3) = -1. / 3.; P(0, 4) = 0.; P(0, 5) = -1. / 3.;
	P(1, 0) = 0.; P(1, 1) = 2.; P(1, 2) = 0.; P(1 ,3) = 0.; P(1, 4) = 0.; P(1, 5) =0.;//XY
	P(2, 0) = 0.; P(2, 1) = 0.; P(2, 2) = 2.; P(2, 3) = 0.; P(2, 4) = 0.; P(2, 5) = 0.;//XZ
	P(3, 0) = -1. / 3.; P(3, 1) = 0.; P(3, 2) = 0.; P(3, 3) = 2./3.; P(3, 4) =0.; P(3, 5) = -1. / 3.;
	P(4, 0) = 0.; P(4, 1) = 0.; P(4, 2) = 2.; P(4, 3) = 0.; P(4, 4) = 0.; P(4, 5) = 0.;//YZ
	P(5, 0) = -1. / 3.; P(5, 1) = 0.; P(5, 2) = 0; P(5, 3) = -1. / 3.; P(5, 4) = 0.; P(5, 5) = 2. / 3.;
	
	//#define _XX_ 0
	//#define _XY_ 1
	//#define _XZ_ 2
	//#define _YY_ 3
	//#define _YZ_ 4
	//#define _ZZ_ 5
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
	//proj = HW[F1HWCylVonMises[{xisol, rho, betasol}]] //. subst2;
	//		sigprojvoigth =
	//		FromCartToVoigt[
	//			Sum[proj[[i]] Outer[Times, vec[[i]], vec[[i]]], { i, 1, 3 }]];
	//	asol = ax[sigprojvoigth] //. subst2;
	//		(*dadsigg = dadsigmax[sigprojvoigth]//.subst2;*)
	//			epse = invCe.sigprojvoigth;
	//	gamma = Norm[epsetrial - epse] / Norm[asol];

	//	(*dadsigg = dadsigmax[sigma]//.subst2;
	//		CT = Ce.((IdentityMatrix[6] - (Outer[Times, asol, asol].Ce) / (asol.Ce.asol)));
	//	T = (IdentityMatrix[6] - gamma dadsigg.Ce);
	//	Dep = CT.T; *)

	//		dadsigg = dadsigmax[sigprojvoigth] //. subst2;
	//		Q = (IdentityMatrix[6] + gamma Ce.dadsigg);
	//	invQ = Inverse[Q];
	//	R = invQ.Ce;
	//	Dep = R - 1 / (asol.R.asol) Outer[Times, R.asol, R.asol];
}

void TPZYCVonMises::N(TPZTensor<STATE> sigma, TPZTensor<STATE> &asol)
{
	sigma.S(asol);
	STATE J2 = sigma.J2();
	asol *= sqrt(3.) / (2.*sqrt(J2));
	asol.XY() *= 2.; asol.XZ() *= 2.; asol.YZ() *= 2.;
}

void TPZYCVonMises::dadsig(TPZTensor<STATE> sigma, TPZFMatrix<STATE> &dadsigmat)
{
	STATE J2 = sigma.J2();
	TPZFMatrix<STATE> P,temp;
	TPZTensor<STATE> Sdev;
	sigma.S(Sdev);
	Sdev.XY() *= 2.; Sdev.XZ() *= 2.; Sdev.YZ() *= 2.;
	Sdev.Print(std::cout);

	Sdev.ProdT(Sdev, temp);
	temp.Print(std::cout);

	GetPMatrix(P);
	P *= sqrt(3.) / (2.*sqrt(J2));
	//TPZFMatrix<REAL> TPZTensor<T>::ProdT( TPZTensor<T> t2)

	temp *= sqrt(3.) / (4.*pow(J2, 3. / 2.));
	P -= temp;
	dadsigmat = P;
	dadsigmat.Print(std::cout);
	//dadsigmax[sigma_] : =
		//P Sqrt[3] / (2 Sqrt[ComputeJ2[sigma]]) -
		//Sqrt[3] / (4 ComputeJ2[sigma] ^ (3 / 2)) Outer[Times, ComputeS[sigma],
		//ComputeS[sigma]]

}

void TPZYCVonMises::ComputeDep(TPZTensor<REAL>::TPZDecomposed DecompSig, TPZTensor<REAL>::TPZDecomposed  DecompEps, TPZTensor<REAL> sigprojvoigt, TPZFMatrix<REAL> &Dep)
{

	TPZTensor<REAL> asol,strainproj,straintrial,diff;//Flow vector
	sigprojvoigt.S(asol);
	STATE J2 = sigprojvoigt.J2();
	asol *= sqrt(3.) / (2.*sqrt(J2));
	asol.XY() *=  2.;asol.XZ() *= 2.;asol.YZ() *= 2.;
	diff = straintrial ;
	diff -= strainproj;
	STATE norm = diff.Norm();
	STATE gamma = norm / asol.Norm();

}
/**
* @brief It computes z = beta * y + alpha * opt(this)*x but z and x can not overlap in memory.
* @param x Is x on the above operation
* @param y Is y on the above operation
* @param z Is z on the above operation
* @param alpha Is alpha on the above operation
* @param beta Is beta on the above operation
* @param opt Indicates if is Transpose or not
*/
//template <class TVar>
//void TPZFMatrix<TVar>::MultAdd(const TPZFMatrix<TVar> &x, const TPZFMatrix<TVar> &y, TPZFMatrix<TVar> &z,
//	const TVar alpha, const TVar beta, const int opt) const {
//ax[sigma_] := Sqrt[3]/(2 Sqrt[ComputeJ2[sigma]]) ComputeS[sigma]
//dadsigmax[sigma_] : = P Sqrt[3] / (2 Sqrt[ComputeJ2[sigma]]) -
// - Sqrt[3] / (4 ComputeJ2[sigma] ^ (3 / 2)) Outer[Times, ComputeS[sigma], ComputeS[sigma]]
//#define _XX_ 0
//#define _XY_ 1
//#define _XZ_ 2
//#define _YY_ 3
//#define _YZ_ 4
//#define _ZZ_ 5

