// SPDX-FileCopyrightText: 2024 <copyright holder> <email>
// SPDX-License-Identifier: Apache-2.0

#ifndef TPZKLINTERPOLATIONSPACE_H
#define TPZKLINTERPOLATIONSPACE_H

#include "pzinterpolationspace.h"
#include "TPZMaterialDataT.h"
#include "TPZElementMatrixT.h"
#include "pzcompel.h"
#include "TPZKarhunenLoeveMat.h"
#include "pzinterpolationspace.h"
#include "TPZMaterial.h"
#include "pzcompel.h"
class TPZMaterialData;
template<class TVar>
class TPZMaterialDataT;
template<class TVar>
class TPZTransfer;
/**
 * @todo write docs
 */
class TPZKLInterpolationSpace : public TPZInterpolationSpace
{
public:


    /**Default constructor*/
    TPZKLInterpolationSpace();

    /**Copy constructor*/
    TPZKLInterpolationSpace(const TPZKLInterpolationSpace& other);

    /**Destructor*/
    ~TPZKLInterpolationSpace();

    /** @brief Puts a copy of the element in the referred mesh */
	TPZKLInterpolationSpace(TPZCompMesh &mesh, const TPZKLInterpolationSpace &copy);

	/** @brief Puts a copy of the element in the patch mesh */
	TPZKLInterpolationSpace(TPZCompMesh &mesh, const TPZKLInterpolationSpace &copy, std::map<int64_t,int64_t> &gl2lcElMap);

	/** Inserts the element within the data structure of the mesh */
	TPZKLInterpolationSpace(TPZCompMesh &mesh, TPZGeoEl *gel);

	virtual void CalcStiffC(TPZCompEl *jel, TPZElementMatrixT<REAL> &ce);


	virtual void CalcStiffB(TPZElementMatrixT<REAL> &be);

protected:
    int fPreferredOrder;

};

#endif // TPZKLINTERPOLATIONSPACE_H
