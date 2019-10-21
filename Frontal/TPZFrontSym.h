/**
 * @file
 * @brief Contains the TPZFrontSym class which implements decomposition process of the frontal matrix (case symmetric).
 */
template<class TVar>
class TPZEqnArray;

#ifndef TPZFRONTSYM_H
#define TPZFRONTSYM_H

#ifdef USING_ATLAS
extern "C"{
#include <cblas.h>
};
#endif

#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <fstream>

#include "pzmatrix.h"
#include "pzstack.h"
#include "pzvec.h"
#include "TPZFront.h"
#include "TPZFileEqnStorage.h"
#include "TPZStackEqnStorage.h"

/** 
 * The Front matrix itself.
 * It is controled by TPZFrontMatrix.\n
 * TPZFrontSym is a symmetrical matrix. It uses a Cholesky decomposition scheme.
 */
/**
 * @brief Abstract class implements storage and decomposition process of the frontal matrix involving simmetry characteristics. \ref frontal "Frontal"
 * @ingroup frontal
 */
template <class TVar>
class TPZFrontSym : public TPZFront<TVar> {
public:
	/** @brief Returns its type*/
	std::string GetMatrixType();
	
    /** Static main used for testing */
	static void main();
	/** @brief Simple destructor */
    ~TPZFrontSym();
    /** @brief Simple constructor */
    TPZFrontSym();
    
    TPZFrontSym(const TPZFrontSym<TVar> &cp) : TPZFront<TVar>(cp)
    {
    }
    /** @brief Constructor with a initial size parameter */
	TPZFrontSym(long GlobalSize);
	
    /// Set the decomposition type
    virtual void SetDecomposeType(DecomposeType dectype)
    {
        if (dectype == ECholesky || dectype == ELDLt) {
            this->fDecomposeType = dectype;
        }
        else
        {
            DebugStop();
        }
    }
    

    /**
     * @brief Decompose these equations and put the result in eqnarray \n
     * Default decompose method is Cholesky
	 * @param mineq index of equations to be decomposed
	 * @param maxeq index of equations to be decomposed
	 * @param result result of decomposition
     */
    void DecomposeEquations(long mineq, long maxeq, TPZEqnArray<TVar> & result);
	
    /**
     * @brief Decompose these equations in a symbolic way and store freed indexes in fFree 
	 * @param mineq Initial equation index
	 * @param maxeq Final equation index
     */
    void SymbolicDecomposeEquations(long mineq, long maxeq);
	
	/** @brief Add a contribution of a stiffness matrix using the indexes to compute the frontwidth */
	void SymbolicAddKel(TPZVec < long > & destinationindex);
	
    /** @brief Compress data structure */
    void Compress();
	
	/** @brief Expand the front matrix */
	void Expand(int largefrontsize);
	
    /** @brief Returns ith, jth element of matrix. \f$ (sourceindex[i],sourceindex[j]) \f$ */
	TVar & Element(long i, long j){
		if(i>j){
			long i_temp=i;
			i=j;
			j=i_temp;
		}
		return this->fData[(j*(j+1))/2+i];
	}
    /** @brief Returns ith, jth element of matrix. \f$ (sourceindex[i],sourceindex[j]) \f$ */
    const TVar & Element(long i, long j) const {
        if(i>j){
            long i_temp=i;
            i=j;
            j=i_temp;
        }
        return this->fData[(j*(j+1))/2+i];
    }
    /** @brief Add a contribution of a stiffness matrix*/
    void AddKel(TPZFMatrix<TVar> &elmat, TPZVec<long> &destinationindex);
	
    /**@brief Add a contribution of a stiffness matrix*/
    virtual void AddKel(TPZFMatrix<TVar> &elmat, TPZVec<long> &sourceindex,  TPZVec<long> &destinationindex);
	
	/** @brief Reorders the elements of the frontmatrix into the full matrix */
	virtual void ExtractFrontMatrix(TPZFMatrix<TVar> &front);
	
private:    

	TVar & Element4JGreatEqualI(long i, long j){
#ifdef PZDEBUG
    if(i>j){
      DebugStop();
    }
#endif
	  return this->fData[(j*(j+1))/2+i];
  }
	
	
	/** 
	 * @brief Decomposes ieq equation and add the result to EqnArray 
	 * @param ieq Index of equation to be decomposed 
	 * @param eqnarray EqnArray to store resulting members
	 */
    void DecomposeOneEquation(long ieq, TPZEqnArray<TVar> &eqnarray);
	
    /**
     * @brief Sets the global equation as freed, allowing the space 
	 */
	/** 
	 * Used by this equation to be used by future assembly processes 
     */
    void FreeGlobal(long global);
    /** @brief return a local index corresponding to a global equation number */
    int Local(long global);
public:
    /** @brief Returns the number of free equations */
	virtual long NFree();
    /** @brief Resets data structure */
	void Reset(long GlobalSize=0);
    /** @brief Allocates data for Front */
	void AllocData();
	
	/** @brief Prints TPZFront data */
	void Print(const char *name, std::ostream& out=std::cout) const;
	void PrintGlobal(const char *name, std::ostream& out = std::cout);
	
	/** @brief Returns decomposition type*/
	DecomposeType GetDecomposeType() const;
	
	/** @brief Does the tensor product betweem two vectors in the positions dependent of ithread*/
	virtual void TensorProductIJ(int ithread, typename TPZFront<TVar>::STensorProductMTData *data);
	
    /** @link dependency */
    /*#  TPZFileEqnStorage lnkTPZFileEqnStorage; */
	
    /** @link dependency */
    /*#  TPZStackEqnStorage lnkTPZStackEqnStorage; */
};

#endif //TPZFRONTSYM_H
