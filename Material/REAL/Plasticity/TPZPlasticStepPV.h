/**
 * @file
 */

#ifndef TPZPlasticStepPV_H
#define TPZPlasticStepPV_H


#include "TPZTensor.h"
#include "pzreal.h"
#include "pzfmatrix.h"
#include "TPZPlasticState.h"
#include "TPZPlasticIntegrMem.h"
#include "TPZPlasticStep.h"
#include "pzlog.h"
#include "pzstepsolver.h"
#include "TPZElasticResponse.h"

#include <set>
#include <ostream>

// Metodos para deixar o prog mais "encapsulado"
REAL NormVecOfMat(TPZFNMatrix <9> mat);
REAL InnerVecOfMat(TPZFMatrix<REAL> &m1, TPZFMatrix<REAL> &m2);
TPZFMatrix<REAL> ProdT(TPZManVector<REAL,3> &v1, TPZManVector<REAL,3> &v2);
TPZFNMatrix <6> FromMatToVoight(TPZFNMatrix <9> mat);

/*
 
 enum EElastoPlastic
 {
 EAuto = 0,
 EForceElastic = 1,
 EForcePlastic = 2
 };
 */

class TPZYCBase
{
public:

	virtual	~TPZYCBase() {};


	virtual TPZElasticResponse GetElasticResponse()=0;
	virtual void SetElasticResponse(TPZElasticResponse &ER) = 0;
	virtual void SetUp(REAL K, REAL G, REAL sigy) = 0;
	virtual TPZElasticResponse GetElasticResponse() const = 0;
	virtual void Phi(TPZVec<STATE> sigvec, STATE alpha, TPZVec<STATE> &phi)const = 0;
	virtual void N(TPZTensor<STATE> sigma, TPZTensor<STATE> &asol) = 0;
	virtual void dadsig(TPZTensor<STATE> sigma, TPZFMatrix<STATE> &dadsigmat) = 0;
	virtual void Read(TPZStream &buf) = 0;
	virtual void Write(TPZStream &buf) const = 0;
	virtual virtual void Print(std::ostream & out) = 0;

	virtual void ApplyStrainComputeSigma(TPZVec<STATE> &epst, TPZVec<STATE> &epsp, STATE & kprev, TPZVec<STATE> &epspnext, TPZVec<STATE> &stressnext, STATE & knext) const = 0;
	virtual void ComputeI1(TPZVec<STATE> stress, STATE &I1)const = 0;
	virtual void ComputeJ2(TPZVec<STATE> stress, STATE &J2)const = 0;
	virtual void ApplyStrainComputeElasticStress(TPZVec<STATE> &strain, TPZVec<STATE> &stress)const = 0;
	virtual void ApplyStressComputeElasticStrain(TPZVec<STATE> &stress, TPZVec<STATE> &strain)const = 0;
	virtual void ProjectSigma(const TPZVec<STATE> &sigmatrial, STATE kprev, TPZVec<STATE> &sigmaproj, STATE &kproj) const = 0;
	virtual void ProjectSigmaDep(const TPZVec<STATE> &sigmatrial, STATE kprev, TPZVec<STATE> &sigmaproj, STATE &kproj, TPZFMatrix<STATE> &GradSigma) const = 0;
	virtual void ComputeDep(TPZTensor<REAL>::TPZDecomposed DecompSig, TPZTensor<REAL>::TPZDecomposed  DecompEps, TPZTensor<REAL> sigprojvoigt, TPZFMatrix<REAL> &Dep) = 0;


	/// The function which defines the plastic surface
	inline static void GetPMatrix(TPZFMatrix<STATE> &P) {
		P.Resize(6, 6);
		P(0, 0) = 2. / 3.; P(0, 1) = 0.; P(0, 2) = 0.; P(0, 3) = -1. / 3.; P(0, 4) = 0.; P(0, 5) = -1. / 3.;//XX
		P(1, 0) = 0.; P(1, 1) = 2.; P(1, 2) = 0.; P(1, 3) = 0.; P(1, 4) = 0.; P(1, 5) = 0.;//XY
		P(2, 0) = 0.; P(2, 1) = 0.; P(2, 2) = 2.; P(2, 3) = 0.; P(2, 4) = 0.; P(2, 5) = 0.;//XZ
		P(3, 0) = -1. / 3.; P(3, 1) = 0.; P(3, 2) = 0.; P(3, 3) = 2. / 3.; P(3, 4) = 0.; P(3, 5) = -1. / 3.;//YY
		P(4, 0) = 0.; P(4, 1) = 0.; P(4, 2) = 2.; P(4, 3) = 0.; P(4, 4) = 0.; P(4, 5) = 0.;//YZ
		P(5, 0) = -1. / 3.; P(5, 1) = 0.; P(5, 2) = 0; P(5, 3) = -1. / 3.; P(5, 4) = 0.; P(5, 5) = 2. / 3.;//ZZ

	}


	/// The function which defines the plastic surface
	inline static void GetRotMatrix(TPZFMatrix<STATE> &Rot) {
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
	/// Transform from Haigh Westergaard cylindrical coordinates to Haigh Westergaard cartesian coordinates
	inline static void FromHWCylToHWCart(const TPZVec<STATE> &HWCylCoords, TPZVec<STATE> &cart) {
		cart[0] = HWCylCoords[0];
		cart[1] = HWCylCoords[1] * cos(HWCylCoords[2]);
		cart[2] = HWCylCoords[1] * sin(HWCylCoords[2]);
	}
	//TPZManVector<STATE> FromHWCartToHWCyl(TPZManVector<STATE>&HWCartCoords);
	/// Transform from eigenvalues to HW Cylindrical coordinates
	inline static void FromPrincipalToHWCyl(const TPZVec<STATE> &PrincipalCoords, TPZVec<STATE> &HWCyl) {
		TPZFNMatrix<9, STATE> Rot(3, 3, 0.), temp(3, 1, 0.), cart(3, 1, 0.);
		temp(0, 0) = PrincipalCoords[0];
		temp(1, 0) = PrincipalCoords[1];
		temp(2, 0) = PrincipalCoords[2];
		GetRotMatrix(Rot);
		Rot.Multiply(temp, cart);
		HWCyl[0] = cart(0, 0);
		HWCyl[1] = sqrt(cart(1, 0) * cart(1, 0) + cart(2, 0) * cart(2, 0));
		HWCyl[2] = atan2(cart(2, 0), cart(1, 0));
	}
	/// Transform from eigenvalues to HW cartesian coordinates
	inline static void FromPrincipalToHWCart(const TPZVec<STATE> &PrincipalCoords, TPZVec<STATE> &HWCart) {
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
	/// Transform from HW Cylindrical coordinates to eigenvalues
	inline static void  FromHWCylToPrincipal(const TPZVec<STATE> &HWCylCoords, TPZVec<STATE> &HWCart) {
		HWCart[0] = (1. / sqrt(3.)) * HWCylCoords[0] + sqrt(2. / 3.) * HWCylCoords[1] * cos(HWCylCoords[2]);
		HWCart[1] = (1. / sqrt(3.)) * HWCylCoords[0] + sqrt(2. / 3.) * HWCylCoords[1] * cos(HWCylCoords[2] - (2. * M_PI / 3.));
		HWCart[2] = (1. / sqrt(3.)) * HWCylCoords[0] + sqrt(2. / 3.) * HWCylCoords[1] * cos(HWCylCoords[2] + (2. * M_PI / 3.));
	}

};


/**
 * @brief Classe que efetua avanco de um passo de plastificacao utilizando o metodo de Newton
 */
template <class YC_t, class ER_t>
class TPZPlasticStepPV : public TPZPlasticBase
{
public:

    /**
     * @brief Constructor which Initialize the plastic material damage variable only
     *
     * @param[in] alpha damage variable
     */

  TPZPlasticStepPV(REAL alpha=0.):fYC(), fER(), fResTol(1.e-12), fMaxNewton(30), fN()
	{ 
        fN.fAlpha = alpha;
    }

    /**
     * @brief Copy Constructor
     *
     * @param[in] source of copy
     */
	TPZPlasticStepPV(const TPZPlasticStepPV & source)
	{
        fYC = source.fYC;
        fER = source.fER;
        fResTol = source.fResTol;
        fMaxNewton = source.fMaxNewton;
        fN = source.fN;
    }

    /**
     * @brief Operator =
     *
     * @param[in] source of copy
     */
	TPZPlasticStepPV & operator=(const TPZPlasticStepPV & source)
	{
        fYC = source.fYC;
        fER = source.fER;
        fResTol = source.fResTol;
        fMaxNewton = source.fMaxNewton;
        fN = source.fN;

        return *this;
    }

    /**
     * @brief Name of the class ina string
     */
	virtual const char * Name() const
	{
        return "TPZPlasticStepPV";
    }

	virtual void Print(std::ostream & out) const
	{
        out << "\n" << this->Name();
        out << "\n YC_t:";
        //fYC.Print(out); FAZER O PRINT
        out << "\n ER_t:";
        fER.Print(out);
        out << "\nTPZPlasticStepPV Internal members:";
        out << "\n fResTol = " << fResTol;
        out << "\n fMaxNewton = " << fMaxNewton;
        out << "\n fN = "; // PlasticState
        fN.Print(out);
    }

    typedef YC_t fNYields;

    /**
     * Imposes the specified strain tensor, evaluating the plastic integration if necessary.
     *
     * @param[in] epsTotal Imposed total strain tensor
     */
    virtual void ApplyStrain(const TPZTensor<REAL> &epsTotal);

    /**
     * Imposes the specified strain tensor and returns the correspondent stress state.
     *
     * @param[in] epsTotal Imposed total strain tensor
     * @param[out] sigma Resultant stress
     */
    virtual void ApplyStrainComputeSigma(const TPZTensor<REAL> &epsTotal, TPZTensor<REAL> &sigma);



    //virtual void ApplySigmaComputeStrain(const TPZTensor<REAL> &sigma, TPZTensor<REAL> &epsTotal);

    /**
     * Imposes the specified strain tensor and returns the corresp. stress state and tangent
     * stiffness matrix.
     *
     * @param[in] epsTotal Imposed total strain tensor
     * @param[out] sigma Resultant stress
     * @param[out] Dep Incremental constitutive relation
     */
    virtual void ApplyStrainComputeDep(const TPZTensor<REAL> &epsTotal, TPZTensor<REAL> &sigma, TPZFMatrix<REAL> &Dep);




    /**
     * Attempts to compute an epsTotal value in order to reach an imposed stress state sigma.
     * This method should be used only for test purposes because it isn't fully robust. Some
     * materials, specially those perfectly plastic and with softening, may fail when applying
     * the Newton Method on ProcessLoad.
     *
     * @param[in] sigma stress tensor
     * @param[out] epsTotal deformation tensor
     */
    virtual void ApplyLoad(const TPZTensor<REAL> & sigma, TPZTensor<REAL> &epsTotal);

    virtual TPZPlasticState<REAL> GetState() const;
    /**
     * @brief Return the value of the yield functions for the given strain
     * @param[in] epsTotal strain tensor (total strain)
     * @param[out] phi vector of yield functions
     */
    virtual void Phi(const TPZTensor<REAL> &epsTotal, TPZVec<REAL> &phi) const;

    virtual void SetElasticResponse(TPZElasticResponse &ER);

    virtual TPZElasticResponse GetElasticResponse() const
    {
        return fER;
    }

    /**
     * @brief Update the damage values
     * @param[in] state Plastic state proposed
     */
    virtual void SetState(const TPZPlasticState<REAL> &state);


    //void CopyFromFNMatrixToTensor(TPZFNMatrix<6> FNM,TPZTensor<STATE> &copy);
    void CopyFromTensorToFMatrix(TPZTensor<STATE> tensor, TPZFMatrix<STATE> &copy);


    //void CopyFromTensorToFNMatrix(TPZTensor<STATE> tensor,TPZFNMatrix<6> &copy);
    void CopyFromFMatrixToTensor(TPZFMatrix<STATE> FNM, TPZTensor<STATE> &copy);




    virtual void Read(TPZStream &buf);

    virtual void Write(TPZStream &buf) const;


    /**
     * Does the TaylorCheck of the tangent matrix
     *
     * @param[in] epsTotal Imposed total strain tensor
     * @param[out] sigma Resultant stress
     * @param[out] Dep Incremental constitutive relation
     */
    void TaylorCheck(TPZTensor<REAL> &EpsIni, TPZTensor<REAL> &deps, REAL kprev, TPZVec<REAL> &conv);


    REAL ComputeNFromTaylorCheck(REAL alpha1, REAL alpha2, TPZFMatrix<REAL> &error1Mat, TPZFMatrix<REAL> &error2Mat);



public:

	void ComputeDep(TPZTensor<REAL>::TPZDecomposed DecompSig, TPZTensor<REAL>::TPZDecomposed  DecompEps, TPZManVector<REAL, 3> sigprvec,TPZFMatrix<REAL> &Dep);

	void F1HWCylVonMises();

	//void SetUp(REAL Phi, REAL Psi, REAL c, TPZElasticResponse &ER) {
	//	fPhi = Phi;
	//	fPsi = Psi;
	//	fc = c;
	//	fER = ER;
	//}

	void SetResidualTolerance(STATE tol)
	{
        fResTol = tol;
    }

    void ResetPlasticMem()
    {
        //fPlasticMem.Resize(0);
    }

    //virtual void Write(TPZStream &buf) const;

    //virtual void Read(TPZStream &buf);
public:

    /** @brief Object which represents the yield criterium */
    YC_t fYC;

    /** @brief Object representing the elastic response */
    ER_t fER;

protected:

    /** @brief Residual tolerance accepted in the plastic loop processes */
    REAL fResTol;

    /** @brief Maximum number of Newton interations allowed in the nonlinear solvers */
    int fMaxNewton; // COLOCAR = 30 (sugestao do erick!)



public:

    /** @brief Plastic State Variables (EpsT, EpsP, Alpha) at the current time step */
    TPZPlasticState<STATE> fN;

    int fYield;


};

#endif //TPZPlasticStepPV_H
