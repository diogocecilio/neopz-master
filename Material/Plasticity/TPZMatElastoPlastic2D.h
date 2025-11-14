/**
 * @file
 */

#ifndef PZELASTOPLASTIC2D_H
#define PZELASTOPLASTIC2D_H

#include "TPZMaterial.h"
#include "TPZMatWithMem.h"
#include "TPZElastoPlasticMem.h"
#include "pzporoelastoplasticmem.h"
#include "TPZMatElastoPlastic.h"
#include "TPZMaterial.h"

/**
 * Implements an elastoplastic material and uses the memory feature to store the damage variables
 * This material works only together with the Plasticity Library.
 */
#ifdef PZ_LOG
// escolha um nome de categoria claro (use pontos para hierarquia)
static TPZLogger loggerplastic2d("materials.plastic2d");
#endif
template <class T, class TMEM = TPZElastoPlasticMem>
class  TPZMatElastoPlastic2D : public TPZMatElastoPlastic<T,TMEM> //, TPZMatWithMem<TMEM>
{
public:
	
	//enum SOLUTIONVARS{ENone = -1};
	/**
	 * Default constructor
	 */
	TPZMatElastoPlastic2D();		
	
	/** Creates a material object and inserts it in the vector of
	 *  material pointers of the mesh. Upon return vectorindex
	 *  contains the index of the material object within the
	 *  vector
	 */
	TPZMatElastoPlastic2D(int id ,  int PlaneStrainOrPlaneStress=1);
	
	/** Creates a material object based on the referred object and
	 *  inserts it in the vector of material pointers of the mesh.
	 *  Upon return vectorindex contains the index of the material
	 *  object within the vector
	 */
	TPZMatElastoPlastic2D(const TPZMatElastoPlastic2D<T,TMEM> &mat);
	
	virtual ~TPZMatElastoPlastic2D();
	
	/** returns the name of the material*/
	virtual std::string Name() const override;
	
	/**returns the integrable dimension of the material*/
	virtual int Dimension() const override { return 2; }
	
	/** returns the number of state variables associated with the material*/
	virtual int NStateVariables() const override { return 2; }


	void GetSolDimensions(uint64_t &u_len,
												uint64_t &du_row,
												uint64_t &du_col) const override
	{u_len = 2; du_row = 2; du_col=2;}
	
    /** @brief Prints out the data associated with the material */
    virtual void Print(std::ostream &out) const override;

	/** print out the data associated with the material*/
	virtual void Print(std::ostream &out, const int memory) const override;
	
	
	/**returns the solution associated with the var index based on
	 * the finite element approximation*/
	virtual void Solution(const TPZMaterialDataT<STATE> &data, int var, TPZVec<REAL> &Solout) override;
	
	/**
	 * It computes a contribution to the stiffness matrix and load vector at one integration point.
	 * @param data [in] stores all input data
	 * @param weight [in] is the weight of the integration rule
	 * @param ek [out] is the stiffness matrix
	 * @param ef [out] is the load vector
	 */
	virtual void Contribute(const TPZMaterialDataT<STATE> &data, REAL weight, TPZFMatrix<REAL> &ef) override;
	
	/**
	 * It computes a contribution to the stiffness matrix and load vector at one integration point.
	 * @param data [in] stores all input data
	 * @param weight [in] is the weight of the integration rule
	 * @param ek [out] is the stiffness matrix
	 * @param ef [out] is the load vector
	 */
	virtual void Contribute(const TPZMaterialDataT<STATE> &data, REAL weight, TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef) override;
	
    /**
     * This method defines which parameters need to be initialized in order to compute the contribution of the boundary condition
     */
    virtual void FillBoundaryConditionDataRequirements(int type,TPZMaterialData &data) const override;
    
	/**
	 * It computes a contribution to the stiffness matrix and load vector at one BC integration point.
	 * @param data [in] stores all input data
	 * @param weight [in] is the weight of the integration rule
	 * @param ek [out] is the stiffness matrix
	 * @param ef [out] is the load vector
	 * @param bc [in] is the boundary condition material
	 */
	virtual void ContributeBC(const TPZMaterialDataT<STATE> &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc) override;

    /**
     * It computes a contribution to the stiffness matrix and load vector at one BC integration point.
     * @param data [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     * @param bc [in] is the boundary condition material
     */
    virtual void ContributeBC(const TPZMaterialDataT<STATE> &data, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc) override;

	
	/** Evaluates the Strain vector based on an available DSol (solution derivatives set) vector.
	 * @param DeltaStrain [out]
	 * @param data [in]
	 */
	virtual void ComputeDeltaStrainVector(const TPZMaterialDataT<STATE> & data, TPZFMatrix<REAL> &DeltaStrain);
    
	
	/** Calls the plasticity template aggregate applyStrainComputeDep method
	 *  @param data [in]
	 *  @param DeltaStrain [in]
	 *  @param Stress [out]
	 *  @param Dep [out]
	 */
	virtual void ApplyDeltaStrainComputeDep(const TPZMaterialDataT<STATE> & data, TPZFMatrix<REAL> & DeltaStrain,
									TPZFMatrix<REAL> & Stress, TPZFMatrix<REAL> & Dep);
	
	/** Calls the plasticity template aggregate applyStrainComputeDep method
	 *  @param data [in]
	 *  @param DeltaStrain [in]
	 *  @param Stress [out]
	 */
	virtual void ApplyDeltaStrain(const TPZMaterialDataT<STATE> & data, TPZFMatrix<REAL> & DeltaStrain,
									TPZFMatrix<REAL> & Stress);
	
	
	/**To create another material of the same type*/
	virtual TPZMaterial * NewMaterial() const override;
	
	
	/**
	 * Unique identifier for serialization purposes
	 */
	public:
virtual int ClassId() const override;

	
	/**
	 * Save the element data to a stream
	 */
	virtual void Write(TPZStream &buf, int withclassid) const override;
	
	/**
	 * Read the element data from a stream
	 */
	virtual void Read(TPZStream &buf, void *context) override;
    


	void BuildConstitutiveMatrix(STATE E, STATE nu, TPZFMatrix<STATE> &D) const {
		D.Redim(3,3); D.Zero();
		const STATE mu = E/(2.0*(1.0+nu));
		if (fPlaneStrain==false){
			const STATE c = E/(1.0 - nu*nu);
			D(0,0)=c; D(0,1)=c*nu; D(1,0)=c*nu; D(1,1)=c; D(2,2)=mu;
		} else { // PlaneStrain
			const STATE c = E/((1.0+nu)*(1.0-2.0*nu));
			D(0,0)=c*(1.0-nu); D(0,1)=c*nu; D(1,0)=c*nu; D(1,1)=c*(1.0-nu);
			D(2,2)=c*(1.0-2.0*nu)/2.0; // = mu
		}
	}

	// void BuildBu(const TPZFMatrix<STATE>& dphiU, TPZFMatrix<STATE>& Bu)
	// {
	// 	const int nU = dphiU.Cols();
	// 	Bu.Redim(3, 2*nU); Bu.Zero();
	// 	for (int a=0; a<nU; ++a) {
	// 		const STATE dNdx = dphiU(0,a), dNdy = dphiU(1,a);
	// 		const int iu = 2*a, iv = 2*a+1;
	// 		Bu(0,iu) = dNdx;          // exx = dudx
	// 		Bu(1,iv) = dNdy;          // eyy = dvdy
	// 		Bu(2,iu) = dNdy;          // gxy = dudy
	// 		Bu(2,iv) = dNdx;          // gxy = dvdx
	// 	}
	// }


	void BuildBu(const TPZFMatrix<STATE>& dphiU, TPZFMatrix<STATE>& Bu)
	{
		const int nU = dphiU.Cols();
		Bu.Redim(3, 2*nU);
		Bu.Zero();

		for (int a = 0; a < nU; ++a) {
			const STATE dNdx = dphiU(0,a);
			const STATE dNdy = dphiU(1,a);
			const int iu = 2*a;       // DOF u
			const int iv = 2*a + 1;   // DOF v

			// exx
			Bu(0, iu) = dNdx;

			// eyy
			Bu(1, iv) = dNdy;

			// exy (tensorial)
			Bu(2, iu) =  dNdy;
			Bu(2, iv) =  dNdx;

		}
	}
protected:
	
	
	int fPlaneStrain;
    
};

template <class T, class TMEM>
int TPZMatElastoPlastic2D<T, TMEM>::ClassId() const{
    return Hash("TPZMatElastoPlastic2D") ^ TPZMatElastoPlastic<T,TMEM>::ClassId() << 1;
}

#endif
