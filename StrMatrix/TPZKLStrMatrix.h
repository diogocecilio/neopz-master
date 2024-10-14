/**
 * @brief Implements TPZKLStrMatrix Structural Matrices. \ref structural "Structural Matrix"
 * @ingroup structural
 */

#ifndef TPZKLStrMatrix_H
#define TPZKLStrMatrix_H
#include "TPZElementMatrixT.h"
#include "TPZStructMatrix.h"
#include <Eigen/Core>

using namespace Eigen;
#include <Eigen/Dense>
 
using namespace std;
using namespace Eigen;
class TPZCompMesh;
template<class TVar>
class TPZFMatrix;
template<class TVar>
class TPZMatrix;


#include "TPZStructMatrixT.h"
#include "pzstrmatrixor.h"

/**
 * @brief Implements a full symmetric structural matrix
 * using TPZSFMatrix as a storage format.
 * @ingroup structural
 */
template<class TVar=STATE, class TPar=TPZStructMatrixOR<TVar>>
class TPZKLStrMatrix : public TPZStructMatrixT<TVar>,
                                  public TPar
{
public:    
	
    using TPZStructMatrixT<TVar>::TPZStructMatrixT;
    TPZMatrix<TVar>* Create() override;

    TPZStructMatrix * Clone() override;

    int ClassId() const override;

    void Read(TPZStream& buf, void* context) override;

    void Write(TPZStream& buf, int withclassid) const override;

    friend TPZPersistenceManager;
	
	virtual void AssembleC (TPZFMatrix<REAL> &C);
	
	virtual void AssembleB(TPZFMatrix<REAL> &B);

	void GetDestIndex(int iel, int nshape, TPZManVector<long> &fSourceIndex, TPZManVector<long> &fDestinationIndex);
	
	void FromEigen ( MatrixXd eigenmat, TPZFMatrix<REAL>  &pzmat )
    {
		int rows = eigenmat.rows();
		int cols = eigenmat.cols();
		pzmat.Resize(rows,cols);
		for(int irow=0;irow<rows;irow++)
		{
			for(int icol=0;icol<cols;icol++)
			{
				pzmat(irow,icol)=eigenmat(irow,icol);
			}
		}
	
    }

	void ToEigen ( TPZFMatrix<REAL>  pzmat,MatrixXd &eigenmat )
    {
		TPZFMatrix<REAL> intpz(pzmat);
		int rows = pzmat.Rows();
		int cols = pzmat.Cols();
		eigenmat.resize(rows,cols);
		for(int irow=0;irow<rows;irow++)
		{
			for(int icol=0;icol<cols;icol++)
			{
				eigenmat(irow,icol)=intpz(irow,icol);
			}
		}
	
    }
	
private:
	
	TPZFMatrix<REAL> fB;
	TPZFMatrix<REAL> fC;
	

};

#endif //TPZFSTRUCTMATRIX_H
