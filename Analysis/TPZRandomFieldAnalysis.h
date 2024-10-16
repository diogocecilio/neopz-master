//
// Created by Diogo Cecilio on 12/10/24.
//

#ifndef PZ_TPZEIGENANALYSIS_H
#define PZ_TPZEIGENANALYSIS_H
#include "TPZAnalysis.h"     //For TPZAnalysis
#include "tpzautopointer.h" //For TPZAutoPointer
#include "pzmatrix.h"       //For TPZFMatrix
#include "TPZLapackEigenSolver.h"
#include "TPZEigenSolver.h"
#include "pzinterpolationspace.h"
#include "TPZElementMatrixT.h"
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;
class TPZRandomFieldAnalysis : public TPZAnalysis{
public:

    TPZRandomFieldAnalysis();

    TPZRandomFieldAnalysis(TPZCompMesh *mesh, bool mustOptimizeBandwidth = true, std::ostream &out = std::cout);

    /** @brief Gets the eigensolver */
    template<class TVar>
    TPZEigenSolver<TVar> &EigenSolver();

    /** @brief Set the solver
    @note In this function it will be checked if the solver is a TPZEigenSolver*/
    void SetSolver(const TPZSolver &solver) override;

    virtual void Assemble()override;

	virtual void Solve()override;

    void AssembleC (TPZFMatrix<double> &C);

    void AssembleB(TPZFMatrix<double> &B);

	void GetDestIndex(int iel, int nshape, TPZManVector<long> &fSourceIndex, TPZManVector<long> &fDestinationIndex);

    REAL IntegrateSolution ( int varid );

    TPZFMatrix<REAL> NormalizeSolution();

    void LoadRealSol(TPZFMatrix<CSTATE> sol);

    void FromEigen ( MatrixXd eigenmat, TPZFMatrix<REAL>  &pzmat );

    void ToEigen ( TPZFMatrix<REAL>  pzmat,MatrixXd &eigenmat );

    int ClassId() const override;

    void Write(TPZStream &buf, int withclassid) const override;

    void Read(TPZStream &buf, void *context) override;

    void PrintField(TPZFMatrix<REAL> matfield,string filename);

    void SetNEigenpairs(int n)
    {
        fNValues=n;
    }



    TPZFMatrix<REAL> GetSolutionValVec()
    {
        return fSolutionValVec;
    }

protected:

    int fNValues;
    TPZVec<REAL> fEigenvalues;
    TPZFMatrix<REAL> fEigenvetors;
    TPZFMatrix<REAL> fSolutionValVec;



};




#endif //PZ_TPZEIGENANALYSIS_H

