/**
 * @file pzmaterial.h
 * @brief Header file for abstract class KLMaterial.\n
 * It implements the weak statement of the differential equation within the PZ environment.
 */

#ifndef TPZKarhunenLoeveMat_H
#define TPZKarhunenLoeveMat_H

#include "pzreal.h"
#include "pzvec.h"
#include "pzmatrix.h"
#include "pzfmatrix.h"
#include "pzadmchunk.h"
#include "tpzautopointer.h"
#include "pzfunction.h"
#include "pzcompel.h"
#include "TPZMaterial.h"
#include <iostream>
#include <string>
#include "TPZMaterialDataT.h"
#include "TPZMatBase.h"
#include "TPZMatSingleSpace.h"
#include "TPZMatErrorSingleSpace.h"
#include "TPZMatLoadCases.h"
//class TPZBndCond;

//class TPZMaterialData;
//class TPZIntPoints;

/**
 * @brief Data of one element at the points of an integration rule, used in the double integral of the
 * Karhunen-Loeve matrix C (see TPZInterpolationSpace::KLPointData and TPZKarhunenLoeveMat::ContributeC).
 */
struct TPZKLPointData {
    /** @brief Coordinates of the points (3 x npoints) */
    TPZFMatrix<REAL> fX;
    /** @brief Weights of the rule multiplied by |det J| (npoints) */
    TPZManVector<REAL,27> fW;
    /** @brief Shape functions at the points (nshape x npoints) */
    TPZFMatrix<REAL> fPhi;
};

/**
 * @ingroup material
 * @brief This abstract class defines the behaviour which each derived class needs to implement
 */
/**
 * Classes derived from the KLMaterial class implement the weak statement of the differential equation
 * within the PZ environment \n
 * It is noteworthy to observe that this definition does not depend on the definition of the interpolation space \n
 * KLMaterial objects also need to implement the interface for post processing the results
 */
class  TPZKarhunenLoeveMat : public TPZMatBase<STATE,
                                          TPZMatSingleSpaceT<STATE>,
                                          TPZMatErrorSingleSpace<STATE>,
                                          TPZMatLoadCases<STATE>>
{
	    using TBase = TPZMatBase<STATE,
                             TPZMatSingleSpaceT<STATE>,
                             TPZMatErrorSingleSpace<STATE>,
                             TPZMatLoadCases<STATE>>;
    enum SOLUTIONVARS{ENone = -1,
	  // Strain
	  EVEC = 0,
	  EVEC1 = 1,
	  EVEC2 = 2,
	  EVEC3 = 3,
	  EVEC4 = 4,
	  EVEC5 = 5,
	  EVECSQR =6
};
private:
    int fId;
    
protected:
   

//	int fExpasionOrder;
	
	REAL fLx;
	
	REAL fLy;
	
	REAL fLz;
	
	int ftype;
	
	int fDim;
    
public:
    /** @brief Big number to penalization method, used for Dirichlet conditions */
    static REAL gBigNumber;
    
    /** @brief Creates a material object and inserts it in the vector of material pointers of the mesh. */
	/** Upon return vectorindex contains the index of the material object within the vector */
    TPZKarhunenLoeveMat(int id,REAL Lx,REAL Ly,REAL Lz,int dim, int type);
	
	TPZKarhunenLoeveMat(int id);
    
    /** @brief Default constructor */
    TPZKarhunenLoeveMat();
    
    /** @brief Creates a material object based on the referred object and inserts it in the vector of material pointers of the mesh. */
	/** Upon return vectorindex contains the index of the material object within the vector */
    TPZKarhunenLoeveMat(const TPZKarhunenLoeveMat &mat);
    /** @brief Default destructor */
    virtual ~TPZKarhunenLoeveMat();
    
    void FillDataRequirements(TPZMaterialData &data) const override;

    void FillBoundaryConditionDataRequirements(int type, TPZMaterialData &data) const override;

    /** @brief Returns the name of the material */
    std::string Name()  const override { return "Karhunen Loeve material"; }
    
    /** @brief Returns the integrable dimension of the material (fDim, 2 unless given in the constructor or SetDim) */
	int Dimension() const  override { return fDim;}
    
    
    void SetDim(int dim)
	{
		fDim=dim;
	}
    
    /** @brief Returns the number of state variables associated with the material */
    int NStateVariables() const override;
    
    /** @brief Returns the number of components which form the flux function */
    virtual int NFluxes() {return 0;}

    virtual TPZMaterial * NewMaterial() const override
    {
        TPZKarhunenLoeveMat *mat = new TPZKarhunenLoeveMat(fId,fLx,fLy,fLz,fDim,ftype);
        return mat;
    }


    
    /** @brief Prints out the data associated with the material */
    void Print(std::ostream & out = std::cout) const override
	{
		out << "KLMaterial" << std::endl;
		out << "Lx = " << fLx << std::endl;
		out << "Ly = " << fLy <<std::endl;
		out << "Lz = " << fLz <<std::endl;
		out << "dim = " << fDim <<std::endl;
		out << "Type = " << ftype <<std::endl;
	}
    
    /** @brief Returns the variable index associated with the name */
    int VariableIndex(const std::string &name) const override;
    
    /** 
	 * @brief Returns the number of variables associated with the variable indexed by var. 
	 * @param var Index variable into the solution, is obtained by calling VariableIndex
	 */
    int NSolutionVariables(int var) const override;
    
    void Solution(const TPZMaterialDataT<STATE> &data, int var,
                  TPZVec<STATE> &Solout) override;
    /** @brief Computes the value of the flux function to be used by ZZ error estimator */
    virtual void Flux(TPZVec<REAL> &x, TPZVec<STATE> &Sol,
                      TPZFMatrix<STATE> &DSol, TPZFMatrix<REAL> &axes,
                      TPZVec<STATE> &flux) {}
    
    virtual void ContributeB(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek);
	
	virtual void ContributeC(TPZMaterialData &data1,TPZMaterialData &data2, REAL w1,REAL w2, TPZFMatrix<STATE> &ek);

	/**
	 * @brief Computes ek(i,j) = sum_p sum_q w1_p w2_q phi1_i(x_p) G(x_p,x_q) phi2_j(x_q), the double integral of the
	 * matrix C between two elements, from their point data (ek is resized to nshape1 x nshape2)
	 */
	void ContributeC(const TPZKLPointData &data1, const TPZKLPointData &data2, TPZFMatrix<STATE> &ek) const;

	/** @brief Autocorrelation function G(x1,x2), without the variance (type 1 to 4, see the implementation) */
	REAL AutocorrelationFunc ( const TPZVec<REAL> &x1, const TPZVec<REAL> &x2 ) const;

	/** @brief Calculates the element stiffness matrix */
	void Contribute(const TPZMaterialDataT<STATE> &data, STATE weight,
                    TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef) override;

    /** @brief Applies the element boundary conditions */
	void ContributeBC(const TPZMaterialDataT<STATE> &data,STATE weight,
                      TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,
                      TPZBndCondT<STATE> &bc) override;

	void Errors(const TPZMaterialDataT<STATE> &data,
                TPZVec<STATE> &values) override;

    
    /** @brief Reads data of the material from a istream (file data)*/
    virtual void SetData(std::istream &data)
	{
		// before, this method called itself (infinite recursion)
		PZError << __PRETTY_FUNCTION__ << " not implemented\n";
		DebugStop();
	}
    
    /** @brief To return a numerical flux type to apply over the interfaces of the elements */
    virtual int FluxType() { return 2; }
	
    virtual void ContributeErrors(TPZMaterialData &data,
                                  REAL weight,
                                  TPZVec<REAL> &nk,
                                  int &errorid){
        PZError << "Error at " << __PRETTY_FUNCTION__ << " - Method not implemented\n";
    }
    
    /** @brief Sets fLinearContext attribute */
    void SetLinearContext(bool IsLinear)
	{
		// before, this method called itself (infinite recursion)
		PZError << __PRETTY_FUNCTION__ << " not implemented\n";
		DebugStop();
	}
    
    /** @{
     * @name Save and Load methods
     */
    
    /** @brief Unique identifier for serialization purposes */
    int ClassId() const override;

// 	int GetExpansioOrder()
// 	{
// 		return fExpasionOrder;
// 	}



     //virtual bool IsMatImpl() = 0;
};

/** @brief Extern variable - Vector of force values */
extern TPZVec< void(*) (const TPZVec<REAL> &, TPZVec<STATE>& ) > GFORCINGVEC;

#endif

