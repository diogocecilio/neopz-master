/* Generated by Together */

#ifndef TPZADAPTMESH_H
#define TPZADAPTMESH_H

#include "pzcmesh.h"
#include "pzcclonemesh.h"
#include "pzvec.h"

#include <pthread.h>

class TPZInterpolatedElement;
class TPZTransform;
template<class T, class V>
class TPZAvlMap;
class TPZOneDRef;

/// interface to generate adapted meshes
class TPZAdaptMesh {
 public:
  
  void RemoveCloneBC(TPZCompMesh *mesh);

  /**
   * Simple constructor
   */
  TPZAdaptMesh();    
  
  /**
   * Simple destructor
   */
  ~TPZAdaptMesh();
  
  /**
   * Defines the computational reference mesh
   */
  void SetCompMesh(TPZCompMesh * mesh);

  /**
   * Defines the maximum p order of an element
   */
  void SetMaxP(int maxp);
  
  /**
   * Public interface to get the optmally refined mesh 
   * @param error: returns the estimated error
   * @param truerror: returns the true error if analitical solution is provided
   * @param ervec: estimated element error for original mesh element vector
   * @param f: analitical solution
   * @param truervec: real element error at each orginal mesh element
   * @param effect: error estimator effectivity
   * @param use_trueerror: evaluates the error throgh the analitical solution provided by f
   */
  TPZCompMesh * GetAdaptedMesh(REAL &error,
			       REAL &truerror,
			       TPZVec<REAL> &ervec, 
			       void (*f)(TPZVec<REAL> &loc,TPZVec<REAL> &val,TPZFMatrix &deriv),
			       TPZVec<REAL> &truervec, 
			       TPZVec<REAL> &effect,
			       int use_trueerror = 0);

  /**
   * ??
   */
  static void DeleteElements(TPZCompMesh *mesh);

  REAL UseTrueError(TPZInterpolatedElement *coarse, void (*f)(TPZVec<REAL> &loc, TPZVec<REAL> &val, TPZFMatrix &deriv));


  //multi-threading control variables...
  static int fThreads_in_use;

 private:
  static int fNClones_to_Analyse;

 protected:
  
  
  /**
   * Retrieves the geometric reference elements to create the patches 
   */
  void GetReferenceElements();
  
  /**
   * Builds the patch of all reference elements. 
   * The patches are stored into patch vectors
   */   
  void BuildReferencePatch();
  
  /**
   * Fill the vector of clone meshes
   */
  void CreateClones();
  
  /**
   * Sorts the elements by the error vector vec, returning permutation vector
   */
  void Sort(TPZVec<REAL> &vec, TPZVec<int> &perm);
  
  /**
   * Sort
   */
  void HeapSort(TPZVec<REAL> &sol, TPZVec<int> &perm);
  
  /**
   * Sorts the errvec returning the ordering indexes in perm param.
   * errpercent is the percentual of the error that must be considered in returning minimum error
   */
  REAL SortMinError (TPZVec<REAL> errvec, TPZVec<int> perm, REAL errpercent);

  /**
   * Creates an adpted computational mesh based on original mesh and in a hp refinement pattern also
   * @param mesh: original mesh
   * @param gelstack: h refinement pattern given by a list of an adapted geometric elements
   * @param porders: p refinement pattern for each element of gelstack
   */
  TPZCompMesh* CreateCompMesh (TPZCompMesh *mesh,TPZVec<TPZGeoEl *> &gelstack,TPZVec<int> &porders);

  /**
   * Verifies if one clone, specified by its index, must be analysed \
   * This method only be called when the true solution is available and the \
   * option usetrueerror in void  GetAdaptedMesh is seted to 1.
   * @param clindex index of the clone to be verified
   * @param minerror minimum error to the clone be analysed
   * @param ervec vector containing the treu error 
   */
  int HasTrueError(int clindex, REAL &minerror, TPZVec<REAL> &ervec);


 private:   
  
  static TPZInterpolatedElement * LargeElement(TPZInterpolatedElement *cint);
  
  /**
   * Computational reference mesh
   */
  TPZCompMesh *fReference;
  
  /**
   * Geometric reference elements vector
   */
  std::set< TPZGeoEl * > fGeoRef;
  
  /**
   * Patches vector
   */
  TPZStack < TPZGeoEl * > fPatch;
  
  /**
   * Maps the start position of each patch into patches vector
   */
  TPZStack < int > fPatchIndex;
  
  /**
   * Element error vector
   */
  TPZStack < REAL > fElementError;

  /**
   * True Element error vector
   */
  TPZVec < REAL > fTrueErrorVec;
  
  /**
   * Clone meshes vector
   */
  TPZStack<TPZCompCloneMesh *> fCloneMesh;
  
  /**
   * Refined clone meshes
   */
  TPZStack <TPZCompMesh *> fFineCloneMesh;

  /**
   * Indexes of the clones that must be analysed
   */
  TPZStack <int> fClonestoAnalyse;
  
  /** 
   * Delete temporary clone meshes from memory
   */
  void CleanUp();

  /**
   * Mesh Error void -- to be used in multi-threading
   */
  static void  * MeshError (void *t);

  void (*fExact)(TPZVec<REAL> &loc, TPZVec<REAL> &result, TPZFMatrix &deriv);

  /*
   * Maximum p order of an element
   */
  int fMaxP;

  static  pthread_mutex_t fLock_clindex;
  static  pthread_cond_t fSignal_free;
  
};
#endif //TPZADAPTMESH_H
