/**
 * @file
 */

#ifndef TPZYCVONMISES_H
#define TPZYCVONMISES_H

#include "TPZTensor.h"
#include "pzfmatrix.h"

#include "pzlog.h"

#include <iostream>
#include "TPZTensor.h"
#include "TPZElasticResponse.h"
#include "TPZPlasticState.h"

/**
 * @brief Implementa  a plastificacao do criterio de Von Mises
 */
class TPZYCVonMises {
    



//#endif //TPZYCVONMISES_H

public:

	enum {
		NYield = 2
	};

	/// Constructor, with all parameters which define the Sandler DiMaggio model
	TPZYCVonMises();
	/// Copy constructor
	TPZYCVonMises(const TPZYCVonMises & copy);


	TPZYCVonMises & operator=(const TPZYCVonMises & source) {
		
		fE = source.fE;
		fnu = source.fnu;
		fElasticResponse = source.fElasticResponse;

		return *this;
	}

	/// Desctructor
	virtual ~TPZYCVonMises();


	TPZElasticResponse GetElasticResponse();

	void SetElasticResponse(TPZElasticResponse &ER);

	virtual TPZElasticResponse GetElasticResponse() const;

	void YieldFunction(const TPZVec<STATE> &sigma, STATE kprev, TPZVec<STATE> &yield) const;

	void Read(TPZStream &buf);

	void Write(TPZStream &buf) const;

	virtual void Print(std::ostream & out) const
	{
		out << "\n Von Mises: Method not implemented!\n";

	}

private:
	/// The function which defines the plastic surface
	static void GetRotMatrix(TPZFMatrix<STATE> &Rot);
	/// Transform from Haigh Westergaard cylindrical coordinates to Haigh Westergaard cartesian coordinates
	static void FromHWCylToHWCart(const TPZVec<STATE> &HWCylCoords, TPZVec<STATE> &Cart);
	//TPZManVector<STATE> FromHWCartToHWCyl(TPZManVector<STATE>&HWCartCoords);
	/// Transform from eigenvalues to HW Cylindrical coordinates
	static void FromPrincipalToHWCyl(const TPZVec<STATE> &PrincipalCoords, TPZVec<STATE> &HWCyl);
	/// Transform from eigenvalues to HW cartesian coordinates
	static void FromPrincipalToHWCart(const TPZVec<STATE> &PrincipalCoords, TPZVec<STATE> &HWCart);
	/// Transform from HW Cylindrical coordinates to eigenvalues
	static void FromHWCylToPrincipal(const TPZVec<STATE> &HWCylCoords, TPZVec<STATE> &PrincipalCoords);

public:

	void Phi(TPZVec<REAL> sigma, STATE alpha, TPZVec<STATE> &phi)const;
	
	//    void ApplyStrainComputeSigma(const TPZVec<STATE> &eps, STATE kprev, TPZVec<STATE> &sigma,STATE &kproj) const;
	void ApplyStrainComputeSigma(TPZVec<STATE> &epst, TPZVec<STATE> &epsp, STATE & kprev, TPZVec<STATE> &epspnext, TPZVec<STATE> &stressnext, STATE & knext) const;



	void ComputeI1(TPZVec<STATE> stress, STATE &I1)const;
	void ComputeJ2(TPZVec<STATE> stress, STATE &J2)const;
	void ApplyStrainComputeElasticStress(TPZVec<STATE> &strain, TPZVec<STATE> &stress)const;
	void ApplyStressComputeElasticStrain(TPZVec<STATE> &stress, TPZVec<STATE> &strain)const;
	void ProjectSigma(const TPZVec<STATE> &sigmatrial, STATE kprev, TPZVec<STATE> &sigmaproj, STATE &kproj) const;

	void ProjectSigmaDep(const TPZVec<STATE> &sigmatrial, STATE kprev, TPZVec<STATE> &sigmaproj, STATE &kproj, TPZFMatrix<STATE> &GradSigma) const;


public:


private:
	STATE fE, fnu,fK,fG; //,fk0;

																   //    bool fIsonCap;
	TPZElasticResponse fElasticResponse;


};


#endif //TPZYCVONMISES_H
