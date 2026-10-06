/**
 * @file
 * @brief Contains the declaration of TPZCompElPostProc class
 */

#ifndef PZCOMPELPOSTPROC_H
#define PZCOMPELPOSTPROC_H

class TPZMaterialData;

//#include "pzreferredcompel.h"
#include "pzinterpolationspace.h"
#include "tpzautopointer.h"
#include "TPZMaterial.h"
#include "TPZMaterialDataT.h"
#include "TPZElementMatrixT.h"
#include "pzstack.h"
#include "pzcmesh.h"
#include "pzquad.h"
#include <cmath>
#include "pzlog.h"
#include "pzpostprocmat.h"
#include "pzvec.h"
#include "pzreal.h"
#include "pzmatrix.h"
#include "pzfmatrix.h"

#include "TPZShapeDisc.h"


#include "pzmultiphysicselement.h"
#include "TPZMatSingleSpace.h"
#include "TPZMatCombinedSpaces.h"



/**
 * @brief Data shared by all the post-processing elements: the element of the simulation mesh whose
 * solution (usually stored at the integration points) is projected on the post-processing element.
 *
 * The referred element is either
 *  - a TPZInterpolationSpace whose material is a single-space material (TPZMatSingleSpaceT), e.g. the
 *    elastoplastic materials with memory (TPZMatElastoPlastic2D inside TPZCompElWithMem), or
 *  - a TPZMultiphysicsElement (e.g. TPZCompElWithMem<TPZMultiphysicsCompEl>) whose material is a
 *    combined-space material (TPZMatCombinedSpacesT), e.g. a TPZMatWithMem material of a
 *    TPZMultiphysicsCompMesh. The variables are then evaluated by the material from the vector of
 *    material data of the atomic spaces, built by the multiphysics element itself.
 */
struct TPZCompElPostProcBase
{
    /// element of the simulation mesh (TPZInterpolationSpace or TPZMultiphysicsElement)
    TPZCompEl *fReferredElement = 0;
};

/**
 * @brief This class implements the TPZCompEl structure to enable copying the solution
 * of the referred compEl at the integration points to itself and interpolating it inside the element
 * @since May, 1 2009
 *
 * The post-processing element uses a clone of the integration rule of the referred element. In
 * CalcResidual the variables of the referred material are evaluated at each integration point (for a
 * material with memory, from the memory item of the point) and projected on the shape functions of this
 * element by a least squares fit computed element by element (TPZPostProcMat::Contribute). When the
 * number of shape functions equals the number of points (order n-1 with n x n or n x n x n Gauss points,
 * see TPZPostProcAnalysis::AutoBuildDisc) the projection is the Lagrange interpolation of the values at
 * the points, i.e. the usual extrapolation of the integration point values to the element.
 */
template <class TCOMPEL >
class TPZCompElPostProc : public TCOMPEL, public TPZCompElPostProcBase
{
    
public:
    
    TPZCompElPostProc();
    
    virtual ~TPZCompElPostProc();
    
    TPZCompElPostProc(TPZCompMesh &mesh, TPZGeoEl *gel);
    
    TPZCompElPostProc(TPZCompMesh &mesh, const TPZCompElPostProc<TCOMPEL> &copy);
    
    /**
     * @brief Used to generate patch mesh... generates a map of connect index from
     * global mesh to clone mesh
     */
    TPZCompElPostProc(TPZCompMesh &mesh,
                      const TPZCompElPostProc<TCOMPEL> &copy,
                      std::map<int64_t,int64_t> & gl2lcConMap,
                      std::map<int64_t,int64_t> & gl2lcElMap);
    
    /** @brief Initializes the shape function type in order to allow non ill-conditioned L2 Transfer matrix */
    void InitializeShapeFunctions();
    
    /** @brief Element of the simulation mesh whose solution is post-processed */
    TPZCompEl *ReferredElement()
    {
        if(!fReferredElement) DebugStop();
        return fReferredElement;
    }

    /**
     * @brief Sets the element of the simulation mesh whose solution is post-processed: a
     * TPZInterpolationSpace (single-space material) or a TPZMultiphysicsElement (combined-space material)
     */
    void SetReferredElement(TPZCompEl *cel)
    {
        if(!dynamic_cast<TPZInterpolationSpace *>(cel) && !dynamic_cast<TPZMultiphysicsElement *>(cel)) DebugStop();
        fReferredElement = cel;
    }
    
    virtual TPZCompEl *Clone(TPZCompMesh &mesh) const override;
    
    /**
     * @brief Create a copy of the given element. The clone copy have the connect indexes
     * mapped to the local clone connects by the given map
     * @param mesh Patch clone mesh
     * @param gl2lcConMap map the connects indexes from global element (original) to the local copy.
     * @param gl2lcElMap map the indexes of the elements between the original element and the patch element
     */
    virtual TPZCompEl *ClonePatchEl(TPZCompMesh &mesh,std::map<int64_t,int64_t> & gl2lcConMap,std::map<int64_t,int64_t>&gl2lcElMap) const override;
    
    /**
     * @brief Prints element data
     * @param out indicates the device where the data will be printed
     */
    virtual void Print(std::ostream & out = std::cout) const override
    {
        out << __PRETTY_FUNCTION__ << " calling print from superclass\n";
        TCOMPEL::Print(out);
    }
    
    void ComputeRequiredData(TPZMaterialDataT<STATE> &data, TPZVec<REAL> &qsi) override;
    void ComputeRequiredData(TPZMaterialDataT<CSTATE> &data, TPZVec<REAL> &qsi) override{
        DebugStop();
    }
    
    /**
     * @brief The CalcResidual reimplementation is in charge of extracting the data from the
     * referred compEl at the integration points and building a minimum square residual
     * method to extrapolate this information throughout the element subdomain.
     * The final ef vector shall be copied onto the solution vector, as it represents
     * the shape functions multipliers of the extrapolation functions.
     *
     * If the referred element is a TPZMultiphysicsElement, the vector of material data is built by the
     * multiphysics element as in TPZMultiphysicsCompEl::CalcStiff (InitMaterialData, AffineTransform and
     * ComputeRequiredData at each point), the indices intLocPtIndex and intGlobPtIndex of all the
     * entries are set from the memory indices of the referred element (GetMemoryIndices) and the
     * variables are evaluated by TPZMatCombinedSpacesT::Solution(datavec, var, sol).
     */
    virtual void CalcResidual(TPZElementMatrixT<STATE> &ef) override{
        CalcResidualInternal(ef);
    }
    
    /**
     * @brief Null implementation of the CalcStiff in order to ensure it wouldn't produce
     * any valid system of equations.
     */
    void CalcStiff(TPZElementMatrixT<STATE> &ek, TPZElementMatrixT<STATE> &ef) override;
    
    /** @brief Compare some fields of 2 TPZMaterialData and return true if these do match. */
    bool dataequal(TPZMaterialData &d1,TPZMaterialData &d2);
    
   
    
    /**
     * @brief Reimplemented in order to ensure the shape functions are computed with local support
     * and thus benefits from the orthogonality properties of the legendre polynomials defined
     * by the method InitializeShapeFunctios().
     */
    /** @brief Compute shape functions based on master element in the classical FEM manner. */
    virtual void ComputeShape(TPZVec<REAL> &intpoint, TPZVec<REAL> &X, TPZFMatrix<REAL> &jacobian, TPZFMatrix<REAL> &axes,
                              REAL &detjac, TPZFMatrix<REAL> &jacinv, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi, TPZFMatrix<REAL> &dphidx) override;
    
   
    int ClassId() const override;
    
    /** @brief Save the element data to a stream */
    void Write(TPZStream &buf, int withclassid) const override;
    
    /** @brief Read the element data from a stream */
    void Read(TPZStream &buf, void *context) override;
protected:
    template<class TVar>
    void CalcResidualInternal(TPZElementMatrixT<TVar> &ef);
    
};

template<class TCOMPEL>
inline TPZCompElPostProc<TCOMPEL>::TPZCompElPostProc() : TCOMPEL() {
    TPZCompElPostProc<TCOMPEL>::InitializeShapeFunctions();
}

template<class TCOMPEL>
inline TPZCompElPostProc<TCOMPEL>::~TPZCompElPostProc() {
    
}

template<class TCOMPEL>
inline TPZCompElPostProc<TCOMPEL>::TPZCompElPostProc(TPZCompMesh &mesh, TPZGeoEl *gel) :
TCOMPEL(mesh, gel){
    TPZCompElPostProc<TCOMPEL>::InitializeShapeFunctions();
}

template<class TCOMPEL>
inline TPZCompElPostProc<TCOMPEL>::TPZCompElPostProc(TPZCompMesh &mesh, const TPZCompElPostProc<TCOMPEL> &copy) :
TCOMPEL(mesh, copy) {
    TPZCompElPostProc<TCOMPEL>::InitializeShapeFunctions();
}

template<class TCOMPEL>
inline TPZCompElPostProc<TCOMPEL>::TPZCompElPostProc(TPZCompMesh &mesh,
                                                     const TPZCompElPostProc<TCOMPEL> &copy,
                                                     std::map<int64_t,int64_t> & gl2lcConMap,
                                                     std::map<int64_t,int64_t> & gl2lcElMap):
TCOMPEL(mesh,copy,gl2lcConMap,gl2lcElMap)
{
    TPZCompElPostProc<TCOMPEL>::InitializeShapeFunctions();
}

template <class TCOMPEL>
inline TPZCompEl * TPZCompElPostProc<TCOMPEL>::Clone(TPZCompMesh &mesh) const{
    return new TPZCompElPostProc<TCOMPEL> (mesh, *this);
}

template <class TCOMPEL>
inline TPZCompEl * TPZCompElPostProc<TCOMPEL>::ClonePatchEl(TPZCompMesh &mesh,std::map<int64_t,int64_t> & gl2lcConMap,std::map<int64_t,int64_t>&gl2lcElMap)const{
    return new TPZCompElPostProc<TCOMPEL> (mesh, *this, gl2lcConMap, gl2lcElMap);
}

template <class TCOMPEL>
inline void TPZCompElPostProc<TCOMPEL>::InitializeShapeFunctions(){
    //TPZReferredCompEl<TCOMPEL>::fShapefunctionType = pzshape::TPZShapeDisc::ETensorial;//pzshape::TPZShapeDisc::EOrdemTotal;
    // an orthogonal (or one closest possible to) polynomial function is very important
    // here to ensure the L2 solution transfer matrix isn't ill-conditioned
    //TPZReferredCompEl<TCOMPEL>::SetOrthogonalFunction(pzshape::TPZShapeDisc::ChebyshevWithoutScale);
    //   TPZReferredCompEl<TCOMPEL>::SetOrthogonalFunction(pzshape::TPZShapeDisc::LegendreWithoutScale);
}

template <class TCOMPEL>
inline void TPZCompElPostProc<TCOMPEL>::ComputeRequiredData(TPZMaterialDataT<STATE> &data,
                                                            TPZVec<REAL> &qsi){
    TCOMPEL::ComputeRequiredData(data, qsi);
}



/**
 * @brief write the element data to a stream
 */
template <class TCOMPEL>
int TPZCompElPostProc<TCOMPEL>::ClassId() const{
    return Hash("TPZCompElPostProc") ^ TCOMPEL::ClassId() << 1;
}


/**
 * @brief write the element data to a stream
 */
template <class TCOMPEL>
inline void TPZCompElPostProc<TCOMPEL>::Write(TPZStream &buf, int withclassid) const
{
    TCOMPEL::Write(buf,withclassid);
}

/**
 * @brief Read the element data from a stream
 */
template <class TCOMPEL>
inline void TPZCompElPostProc<TCOMPEL>::Read(TPZStream &buf, void *context)
{
    TCOMPEL::Read(buf,context);
}




template <class TCOMPEL>
template<class TVar>
inline void TPZCompElPostProc<TCOMPEL>::CalcResidualInternal(TPZElementMatrixT<TVar> &ef)
{
    ef.Reset();
    
    this->InitializeElementMatrix(ef);
    
    TPZCompEl * pCompElRef = ReferredElement();
    
    TPZInterpolationSpace * pIntSpRef = dynamic_cast<TPZInterpolationSpace *>(pCompElRef);
    
    TPZMultiphysicsElement *pMultiRef = dynamic_cast<TPZMultiphysicsElement *>(pCompElRef);
    
    TPZPostProcMat * pPostProcMat = dynamic_cast<TPZPostProcMat *>(this->Material());
    
    if(!pPostProcMat) DebugStop();
    
    // material of a referred TPZInterpolationSpace
    auto* pMaterialRef =
        dynamic_cast<TPZMatSingleSpaceT<STATE>*>(pCompElRef->Material());
    // material of a referred TPZMultiphysicsElement
    auto* pMatCombinedRef =
        dynamic_cast<TPZMatCombinedSpacesT<STATE>*>(pCompElRef->Material());
    if((pMultiRef && !pMatCombinedRef) || (!pIntSpRef && !pMultiRef))
    {
        PZError << "Error at " << __PRETTY_FUNCTION__ << " the referred element must be a TPZInterpolationSpace or "
                << "a TPZMultiphysicsElement with a combined-space material\n";
        DebugStop();
    }
    
    
    if (this->NConnects() == 0) return;///boundary discontinuous elements have this characteristic
    
    int64_t numeq = ef.fMat.Rows();
    TPZFNMatrix<1,STATE> efTemp(numeq,1,0.);
    
    TPZMaterialDataT<STATE> data, dataRef;
    this->InitMaterialData(data);
    if(pIntSpRef)
    {
        pIntSpRef->InitMaterialData(dataRef);
        dataRef.fNeedsSol = true;
        dataRef.p = pIntSpRef->MaxOrder();
    }
    data.p    = this->MaxOrder();
    int dim = this->Reference()->Dimension();
    TPZManVector<REAL,3> intpoint(dim,0.);
    TPZManVector<REAL,3> intpointRef(dim,0);
    REAL weight = 0.;
    REAL weightRef = 0;
  
    TPZIntPoints* pFIntegrationRule = &this->GetIntegrationRule();
    if (!pFIntegrationRule) {
        PZError << "Error at " << __PRETTY_FUNCTION__ << " Integration rule must be used from TPZCompEl. \n";
        DebugStop();
    }
    
    const TPZIntPoints &intrule    = this->GetIntegrationRule();
    const TPZIntPoints &intruleRef = pCompElRef->GetIntegrationRule();
    
    int intrulepoints = intrule.NPoints();
    int intrulepointsRef = intruleRef.NPoints();
    if(intrulepoints != intrulepointsRef)
    {
        PZError << "Error at " << __PRETTY_FUNCTION__ << " Referred CompEl with different number of integration points\n";
        return;
    }
    
    // Referred multiphysics element: vector of material data of the atomic spaces, built as in
    // TPZMultiphysicsCompEl::CalcStiff, affine transforms to the atomic elements and memory indices
    TPZManVector<TPZMaterialDataT<STATE>,4> datavecRef;
    TPZManVector<TPZTransform<>,4> trvecRef;
    TPZManVector<int64_t,64> memIndices;
    if(pMultiRef)
    {
        const int64_t nref = pMultiRef->NMeshes();
        datavecRef.resize(nref);
        // InitMaterialData of the atomic elements and FillDataRequirements of the material
        pMultiRef->InitMaterialData(datavecRef);
        for(int64_t iref = 0; iref < nref; iref++)
        {
            TPZInterpolationSpace *msp = dynamic_cast<TPZInterpolationSpace *>(pMultiRef->Element(iref));
            if(!msp) continue;
            datavecRef[iref].p = msp->MaxOrder();
            datavecRef[iref].fNeedsSol = true;
        }
        pMultiRef->AffineTransform(trvecRef);
        // memory items of the integration points (empty if the element has no memory)
        pCompElRef->GetMemoryIndices(memIndices);
        if(memIndices.size() && memIndices.size() != intrulepoints)
        {
            PZError << "Error at " << __PRETTY_FUNCTION__ << " the referred element has " << memIndices.size()
                    << " memory items and " << intrulepoints << " integration points\n";
            DebugStop();
        }
    }
    
    int nshape = this->NShapeF();
    TPZFNMatrix<10,STATE> ekTemp(nshape, nshape, 0.);
    
    TPZManVector<int,10> varIndex;
    pPostProcMat->GetPostProcessVarIndexList(varIndex);
    
    int stackedVarSize = pPostProcMat->NStateVariables();
    TPZVec<STATE> Sol;
    
    for(int int_ind = 0; int_ind < intrulepoints; ++int_ind)
    {
        intrule.   Point(int_ind,intpoint,   weight);
        intruleRef.Point(int_ind,intpointRef,weightRef);
        data   .intLocPtIndex = int_ind;
        dataRef.intLocPtIndex = int_ind;
        this->ComputeRequiredData(data,    intpoint);
        
        /*pIntSpRef ->ComputeShape(intpointRef, dataRef.x, dataRef.jacobian,
         dataRef.axes, dataRef.detjac, dataRef.jacinv,
         dataRef.phi, dataRef.dphix); */
        dataRef.intLocPtIndex = int_ind;
        if(pIntSpRef)
        {
            pIntSpRef->ComputeShape(intpointRef,dataRef);
            pIntSpRef->ComputeRequiredData(dataRef, intpointRef);
            weightRef *= fabs(dataRef.detjac);
        }
        if(pMultiRef)
        {
            // local index of the point in the rule of the referred element (the cloned rule has the same points)
            const int64_t nref = datavecRef.size();
            for(int64_t iref = 0; iref < nref; iref++) datavecRef[iref].intLocPtIndex = int_ind;
            pMultiRef->ComputeRequiredData(intpointRef, trvecRef, datavecRef);
            // global index of the memory item of the point (-1 without memory)
            const int64_t globindex = memIndices.size() ? memIndices[int_ind] : -1;
            for(int64_t iref = 0; iref < nref; iref++) datavecRef[iref].intGlobPtIndex = globindex;
        }
        weight    *= fabs(data.detjac);
        
//        if(pIntSpRef)
//        {
//            if(!dataequal(data,dataRef)){
//                PZError << "Error at " << __PRETTY_FUNCTION__ << " this and Referred CompEl TPZMaterialData(s) do not match\n";
//                ef.Reset();
//                return;
//            }
//        }
        data.sol[0].Resize(stackedVarSize,0.);
        int64_t index = 0;
        // stacking the solutions to post process.
#ifdef PZ_LOG
        TPZLogger pzcompelpostproclogger("pz.mesh.TPZCompElPostProc");
        if(pIntSpRef && pzcompelpostproclogger.isDebugEnabled())
        {
            std::stringstream sout;
            sout << "Integration point " << int_ind << " x = " << dataRef.x << " GradSol = " << dataRef.dsol[0] ;
            LOGPZ_DEBUG(pzcompelpostproclogger, sout.str())
        }
#endif
        int n_var_indexes = varIndex.NElements();
        for(int var_ind = 0; var_ind < n_var_indexes; var_ind++)
        {
            int variableindex = varIndex[var_ind];
            int nsolvars = pCompElRef->Material()->NSolutionVariables(variableindex);
            Sol.Resize(nsolvars);

            // diferenca entre variavel de interpolacao e variavel de elemento
            if (variableindex < 99) {
                if(pMultiRef) pMatCombinedRef->Solution(datavecRef, variableindex, Sol);
                else pMaterialRef->Solution(dataRef, variableindex, Sol);
            }
            else {
                pCompElRef->Solution(intpointRef, variableindex, Sol);
            }
            if(pMultiRef && Sol.size() < nsolvars)
            {
                PZError << "Error at " << __PRETTY_FUNCTION__ << " variable " << variableindex << " returned "
                        << Sol.size() << " values instead of " << nsolvars << "\n";
                DebugStop();
            }

#ifdef PZ_LOG
            if(pzcompelpostproclogger.isDebugEnabled())
            {
                std::stringstream sout;
                std::string varname;
                pPostProcMat->GetPostProcVarName(var_ind, varname);
                sout << varname << " -value- " << Sol;
                LOGPZ_DEBUG(pzcompelpostproclogger, sout.str())
            }
#endif
            for(int i = 0; i <nsolvars; i++) data.sol[0][index+i] = Sol[i];
            index += nsolvars;
        }
        
        pPostProcMat->Contribute(data,weight,ekTemp,efTemp);
        
    }//loop over integration points
    
    if(pMultiRef)
    {
        // as TPZMultiphysicsCompEl::CleanupMaterialData
        for(int64_t iref = 0; iref < datavecRef.size(); iref++)
        {
            TPZInterpolationSpace *msp = dynamic_cast<TPZInterpolationSpace *>(pMultiRef->Element(iref));
            if(msp) msp->CleanupMaterialData(datavecRef[iref]);
        }
    }
    
    TPZFNMatrix<90,STATE> ekCopy(ekTemp);
    
    TPZFNMatrix<10,STATE> rhsTemp(nshape, 1, 0.);
    for(int i_st = 0; i_st < stackedVarSize; i_st++)
    {
        
        efTemp.GetSub(i_st*nshape, 0, nshape, 1, rhsTemp);
        
        TPZFNMatrix<9,STATE> rhsCopy(rhsTemp);
//        int status = ekTemp.Solve_Cholesky(&(rhsTemp));
        int status = ekTemp.Solve_LU(&(rhsTemp));
#ifdef PZDEBUG
        {
            TPZFNMatrix<9,STATE> result(nshape,1,0.);
            ekCopy.MultAdd(rhsTemp, rhsCopy, result, 1., -1.);
            REAL invRes = Norm(result);
            if(!status ){
                PZError << "Error at " << __PRETTY_FUNCTION__ << " Unable to solve the transference linear system\n";
                ef.Reset();
                return;
            }
            
            if(invRes > 1.e-7)
            {
                PZError << "Error at " << __PRETTY_FUNCTION__
                << " Transference linear system solved with residual norm = "
                << invRes << " at " << i_st << " export variable\n";
            }
        }
#endif
        for(int i_sh = 0; i_sh < nshape; i_sh++)
            ef.fMat(i_sh * stackedVarSize + i_st, 0) = rhsTemp(i_sh);
    }

#ifdef PZ_LOG2
    {
        std::stringstream sout;
        sout << "Element index " << this->fIndex << std::endl;
        ef.fMat.Print("Post Processed ",sout);
        LOGPZ_DEBUG(pzcompelpostproclogger, sout.str())
    }
#endif
    
}

template <class TCOMPEL>
inline void TPZCompElPostProc<TCOMPEL>::CalcStiff(TPZElementMatrixT<STATE> &ek, TPZElementMatrixT<STATE> &ef){
    PZError << "\nTPZCompElPostProc<TCOMPEL>::CalcStiff() Should never be called!!!\n";
    return;
}

template <class TCOMPEL>
inline bool TPZCompElPostProc<TCOMPEL>::dataequal(TPZMaterialData &d1,TPZMaterialData &d2)
{
    const REAL SMALLNUMBER = 1.e-8;
    int i;
    if(d1.p!=d2.p)
    {
        DebugStop();
        return 0;
    }
    REAL res = 0;
    int dim = d1.x.NElements();
    int64_t nshape = d1.phi.Rows();
    int64_t nshape2 = d2.phi.Rows();
    if(dim != d2.x.NElements() || nshape!= nshape2)
    {
        DebugStop();
        return 0; // dimensions and number of integration points shall match
    }
    
    for(i = 0; i < dim; i++)res += pow(d1.x[i]-d2.x[i],(REAL)2.0); // integration points must be at the same locations
    if(res > SMALLNUMBER)
    {
        DebugStop();
        return 0;
    }
    return 1;
}



template <class TCOMPEL>
inline void TPZCompElPostProc<TCOMPEL>::ComputeShape(TPZVec<REAL> &intpoint, TPZVec<REAL> &X,
                                                     TPZFMatrix<REAL> &jacobian, TPZFMatrix<REAL> &axes,
                                                     REAL &detjac, TPZFMatrix<REAL> &jacinv,
                                                     TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi, TPZFMatrix<REAL> &dphix){
    TPZGeoEl * ref = this->Reference();
    if (!ref){
        PZError << "\nERROR AT " << __PRETTY_FUNCTION__ << " - this->Reference() == NULL\n";
        return;
    }//if
    ref->Jacobian( intpoint, jacobian, axes, detjac , jacinv);
    
    ref->X(intpoint, X);
    //  this->Shape(intpoint,intpoint,phi,dphix);
    this->Shape(intpoint,phi,dphi);
    //this->Shape(intpoint,X,phi,dphix);
    
    //    ///axes is identity in discontinuous elements
    //    axes.Resize(dphix.Rows(), dphix.Rows());
    //    axes.Identity();
}

#endif





