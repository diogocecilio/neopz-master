// //$Id: pzpostprocanalysis.h,v 1.5 2010-10-20 18:41:37 diogo Exp $
// #ifndef PZPOSTPROCANALYSIS_H
// #define PZPOSTPROCANALYSIS_H
//
// #include "TPZLinearAnalysis.h"
// #include "pzpostprocanalysis.h"
// #include "pzpostprocmat.h"
// #include "pzcompelpostproc.h"
// #include "pzcmesh.h"
// #include "pzgmesh.h"
// #include "pzvec.h"
// #include "tpzautopointer.h"
//
// #include "pzstring.h"
// //#include "pzelastoplasticanalysis.h"
// #include "pzcreateapproxspace.h"
// #include "pzmultiphysicselement.h"
//
// #include <map>
// #include <set>
// #include <stdio.h>
// #include "pzlog.h"
//
// #include "pzintel.h"
//
// #include "pzrefpoint.h"
// #include "pzgeopoint.h"
// #include "pzshapepoint.h"
// #include "tpzpoint.h"
//
// #include "pzshapelinear.h"
// #include "TPZGeoLinear.h"
// #include "TPZRefLinear.h"
// #include "tpzline.h"
//
// #include "pzshapetriang.h"
// #include "pzreftriangle.h"
// #include "pzgeotriangle.h"
// #include "tpztriangle.h"
//
// #include "pzrefquad.h"
// #include "pzshapequad.h"
// #include "pzgeoquad.h"
// #include "tpzquadrilateral.h"
//
// #include "pzshapeprism.h"
// #include "pzrefprism.h"
// #include "pzgeoprism.h"
// #include "tpzprism.h"
//
// #include "pzshapetetra.h"
// #include "pzreftetrahedra.h"
// #include "pzgeotetrahedra.h"
// #include "tpztetrahedron.h"
//
// #include "pzshapepiram.h"
// #include "pzrefpyram.h"
// #include "pzgeopyramid.h"
// #include "tpzpyramid.h"
//
// #include "TPZGeoCube.h"
// #include "pzshapecube.h"
// #include "TPZRefCube.h"
// #include "tpzcube.h"
// #include "pzelctemp.h"
// #include "tpzcompmeshreferred.h"
// /**
//  * The Post Processing Analysis makes the interface among the computational
//  * analysis and itself, being also in charge of getting the current solution
//  * from the main analysis.
//  */
//
// class TPZPostProcAnalysis : public TPZLinearAnalysis {
//
// public:
//
// TPZPostProcAnalysis(TPZCompMesh * pRef);
//
// TPZPostProcAnalysis();
//
//     TPZPostProcAnalysis(const TPZPostProcAnalysis &copy);
//
//     TPZPostProcAnalysis &operator=(const TPZPostProcAnalysis &copy);
//
// virtual ~TPZPostProcAnalysis();
//
//     /// Set the computational mesh we are going to post process
//     void SetCompMesh(TPZCompMesh *pRef);
//
//     TPZCompMesh *ReferenceCompMesh()
//     {
//         return fpMainMesh;
//     }
// /**
//  *	Assemble() blank implementation in order to avoid its usage. In such an Analysis
//  * class the Assemble() method is useless.
//  */
// virtual  void Assemble();
//
// /**
//  *	Solve() blank implementation in order to avoid its usage. In such an Analysis
//  * class the Solve() method is useless.
//  */
//
// virtual void Solve();
//
// /**
//  * TransferSolution is in charge of transferring the solution from the base analysis/mesh
//  * to the current post processing mesh, solving for the dof to provide a LSM extrapolation.
//  */
// void TransferSolution();
//
// /**
//  * Informs the material IDs and the variable names that are to be postprocessed, if
//  * matching. Should be called right after the class instantiation
//  */
// void SetPostProcessVariables(TPZVec<int> & matIds, TPZVec<std::string> &varNames);
//
// static void SetAllCreateFunctionsPostProc(TPZCompMesh *cmesh);
// static void SetAllCreateFunctionsContinuous();
// 		void AutoBuildDisc();
//
//     /** @brief Returns the unique identifier for reading/writing objects to streams */
// 	virtual int ClassId() const;
// 	/** @brief Save the element data to a stream */
// 	/** @brief Save the element data to a stream */
// 	void Write(TPZStream &buf, int withclassid) const override;
//
// 	/** @brief Read the element data from a stream */
// 	void Read(TPZStream &buf, void *context) override;
//
// // 	virtual void Print( const std::string &name, std::ostream &out)
// // 	{
// // 		Print(name,out);
// // 	}
//
//
//
// protected:
//
// 	TPZCompMesh * fpMainMesh;
//
//     TPZVec<TPZCompEl *> fReferredElements;
//
// /**
//  * TPZCompElPostProc<TCOMPEL> creation function setup
//  */
//
// public:
//
//
// static TPZCompEl * CreatePostProcDisc(  TPZGeoEl *gel, TPZCompMesh &mesh);
//
// static TPZCompEl * CreatePointEl( TPZGeoEl *gel, TPZCompMesh &mesh);
// static TPZCompEl * CreateLinearEl( TPZGeoEl *gel, TPZCompMesh &mesh);
// static TPZCompEl * CreateQuadEl( TPZGeoEl *gel, TPZCompMesh &mesh);
// static TPZCompEl * CreateTriangleEl( TPZGeoEl *gel, TPZCompMesh &mesh);
// static TPZCompEl * CreateCubeEl( TPZGeoEl *gel, TPZCompMesh &mesh);
// static TPZCompEl * CreatePyramEl( TPZGeoEl *gel, TPZCompMesh &mesh);
// static TPZCompEl * CreateTetraEl( TPZGeoEl *gel, TPZCompMesh &mesh);
// static TPZCompEl * CreatePrismEl( TPZGeoEl *gel, TPZCompMesh &mesh);
//
// };
//
// #endif
//$Id: pzpostprocanalysis.h,v 1.5 2010-10-20 18:41:37 diogo Exp $
#ifndef PZPOSTPROCANALYSIS_H
#define PZPOSTPROCANALYSIS_H

#include "TPZLinearAnalysis.h"
#include "pzcompel.h"
#include "TPZGeoElement.h"
#include "pzfmatrix.h"
#include "pzvec.h"

#include <iostream>
#include <string>


/**
 * The Post Processing Analysis makes the interface among the computational
 * analysis and itself, being also in charge of getting the current solution
 * from the main analysis.
 *
 * Typical use (variables stored at the integration points of a material with memory):
 * @code
 * TPZPostProcAnalysis post;
 * post.SetCompMesh(cmesh);                         // mesh of the simulation
 * post.SetPostProcessVariables(matids, varnames);  // creates the discontinuous post-processing mesh
 * TPZFStructMatrix<STATE> sm(post.Mesh());
 * sm.SetNumThreads(0);
 * post.SetStructuralMatrix(sm);
 * post.DefineGraphMesh(dim, scalnames, vecnames, "file.vtk");
 * // after each converged step of the simulation:
 * post.TransferSolution();                         // projection of the current values
 * post.PostProcess(0);                             // file.scal_vec.<step>.vtk
 * @endcode
 *
 * The post-processing mesh has one discontinuous element (TPZCompElPostProc) for each element of the
 * simulation mesh with a post-processed material. It refers to the element of the simulation mesh and
 * uses a clone of its integration rule. The referred element may be
 *  - an interpolation space (TPZInterpolationSpace) with a single-space material, or
 *  - a multiphysics element (TPZMultiphysicsElement, e.g. of a TPZMultiphysicsCompMesh built with
 *    BuildMultiphysicsSpaceWithMemory) with a combined-space material (TPZMatCombinedSpacesT, possibly
 *    TPZMatWithMem): the variables are then evaluated by the material from the vector of data of the
 *    atomic spaces, with the memory index of each integration point.
 *
 * The order of the post-processing element is the preferred order of the referred element (the highest
 * order of the atomic spaces for a multiphysics element), reduced if needed so that the number of shape
 * functions does not exceed the number of integration points (otherwise the element-wise projection is
 * singular). With n x n (x n) Gauss points and a preferred order of at least n-1, the order n-1 is used and
 * the projection is the Lagrange interpolation of the integration point values.
 *
 * SetPostProcessVariables and TransferSolution change the references of the geometric mesh: when
 * SetPostProcessVariables returns, the geometric elements have no reference (the post-processing elements
 * are created discontinuous); when TransferSolution returns, they refer to the elements of the simulation
 * mesh. A caller that relies on the references (e.g. of another mesh) must restore them.
 */

class TPZPostProcAnalysis : public TPZLinearAnalysis {

public:

TPZPostProcAnalysis(TPZCompMesh * pRef);

TPZPostProcAnalysis();

    TPZPostProcAnalysis(const TPZPostProcAnalysis &copy);

    TPZPostProcAnalysis &operator=(const TPZPostProcAnalysis &copy);

    virtual ~TPZPostProcAnalysis();

    /// Set the computational mesh we are going to post process
    void SetCompMesh(TPZCompMesh *pRef, bool mustOptimizeBandwidth = false) override;

    TPZCompMesh *ReferenceCompMesh()
    {
        return fpMainMesh;
    }
/**
 *	Assemble() blank implementation in order to avoid its usage. In such an Analysis
 * class the Assemble() method is useless.
 */
virtual  void Assemble() override;

/**
 *	Solve() blank implementation in order to avoid its usage. In such an Analysis
 * class the Solve() method is useless.
 */

virtual void Solve() override;

/**
 * TransferSolution is in charge of transferring the solution from the base analysis/mesh
 * to the current post processing mesh, solving for the dof to provide a LSM extrapolation.
 * It does not change the simulation mesh nor the memory of its material; on return the geometric
 * elements refer to the elements of the simulation mesh.
 */
void TransferSolution();

/**
 * Informs the material IDs and the variable names that are to be postprocessed, if
 * matching. Should be called right after the class instantiation
 */
void SetPostProcessVariables(TPZVec<int> & matIds, TPZVec<std::string> &varNames);

static void SetAllCreateFunctionsPostProc(TPZCompMesh *cmesh);
static void SetAllCreateFunctionsContinuous();
    /**
     * @brief Creates the discontinuous post-processing elements: one TPZCompElPostProc for each element of the
     * simulation mesh with a post-processed material, referring to it (interpolation space or multiphysics
     * element), with a clone of its integration rule and the order described in the class documentation.
     * For a multiphysics element with memory it verifies that the cloned rule is the rule of the assembly and of
     * the memory.
     */
		void AutoBuildDisc();

    /** @brief Returns the unique identifier for reading/writing objects to streams */
	public:
int ClassId() const override;

	/** @brief Save the element data to a stream */
	void Write(TPZStream &buf, int withclassid) const override;

	/** @brief Read the element data from a stream */
	void Read(TPZStream &buf, void *context) override;


protected:

	TPZCompMesh * fpMainMesh;

    TPZVec<TPZCompEl *> fReferredElements;

/**
 * TPZCompElPostProc<TCOMPEL> creation function setup
 */

public:


static TPZCompEl * CreatePostProcDisc(  TPZGeoEl *gel, TPZCompMesh &mesh);

static TPZCompEl * CreatePointEl( TPZGeoEl *gel, TPZCompMesh &mesh);
static TPZCompEl * CreateLinearEl( TPZGeoEl *gel, TPZCompMesh &mesh);
static TPZCompEl * CreateQuadEl( TPZGeoEl *gel, TPZCompMesh &mesh);
static TPZCompEl * CreateTriangleEl( TPZGeoEl *gel, TPZCompMesh &mesh);
static TPZCompEl * CreateCubeEl( TPZGeoEl *gel, TPZCompMesh &mesh);
static TPZCompEl * CreatePyramEl( TPZGeoEl *gel, TPZCompMesh &mesh);
static TPZCompEl * CreateTetraEl( TPZGeoEl *gel, TPZCompMesh &mesh);
static TPZCompEl * CreatePrismEl( TPZGeoEl *gel, TPZCompMesh &mesh);

};

#endif
