// //$Id: pzpostprocanalysis.cpp,v 1.10 2010-11-23 18:58:35 diogo Exp $
// #include "pzpostprocanalysis.h"
// #include <map>
// #include <set>
// #include <stdio.h>
// #include "pzlog.h"
//
// #ifdef LOG4CXX
// static LoggerPtr PPAnalysisLogger(Logger::getLogger("pz.analysis.postproc"));
// #endif
//
// using namespace std;
//
// TPZPostProcAnalysis::TPZPostProcAnalysis() :TPZLinearAnalysis(),fpMainMesh(NULL)
// {
// }
//
// TPZPostProcAnalysis::TPZPostProcAnalysis(TPZCompMesh * pRef):TPZLinearAnalysis(), fpMainMesh(pRef)
// {
//
//     SetCompMesh(pRef);
//
// }
//
// TPZPostProcAnalysis::TPZPostProcAnalysis(const TPZPostProcAnalysis &copy) : TPZLinearAnalysis(), fpMainMesh(0)
// {
//
// }
//
// TPZPostProcAnalysis &TPZPostProcAnalysis::operator=(const TPZPostProcAnalysis &copy)
// {
//     SetCompMesh(0);
//     return *this;
// }
//
// TPZPostProcAnalysis::~TPZPostProcAnalysis()
// {
//     if (fCompMesh) {
//         delete fCompMesh;
//     }
//
// }
//
// /// Set the computational mesh we are going to post process
// void TPZPostProcAnalysis::SetCompMesh(TPZCompMesh *pRef)
// {
//     // the postprocess mesh already exists, do nothing
//     if (fpMainMesh == pRef) {
//         return;
//     }
//
//     if (fCompMesh) {
//         delete fCompMesh;
//         fCompMesh = 0;
//         TPZLinearAnalysis::CleanUp();
//     }
//
//     fpMainMesh = pRef;
//
//     if (!pRef) {
//         return;
//     }
//
//     TPZCompMesh* pcMainMesh = fpMainMesh;
//
//     TPZGeoMesh * pgmesh = pcMainMesh->Reference();
//
//     TPZCompMesh * pcPostProcMesh = new TPZCompMesh(pgmesh);
//
//     fCompMesh = pcPostProcMesh;
//
//     TPZPostProcAnalysis::SetAllCreateFunctionsPostProc(pcPostProcMesh);
//
// }
//
//
// void TPZPostProcAnalysis::SetPostProcessVariables(TPZVec<int> & matIds, TPZVec<std::string> &varNames)
// {
//     int nMat, matNumber;
// 	TPZCompMesh * pcMainMesh = fpMainMesh;
//
// 	TPZCompMesh * pcPostProcMesh = (this->Mesh());
//
//     if (!pcPostProcMesh) {
//         DebugStop();
//     }
//
// 	nMat = matIds.NElements();
// 	for(int i = 0; i < nMat; i++)
// 	{
// 		TPZMaterial * pmat = pcMainMesh->FindMaterial(matIds[i]);
// 		if(!pmat)
// 		{
// 			PZError << "Error at " << __PRETTY_FUNCTION__ << " TPZPostProcAnalysis::SetPostProcessVariables() material Id " << matIds[i] << " not found in original mesh!\n";
// 			continue;
// 		}
//
// 		TPZPostProcMat * pPostProcMat = new TPZPostProcMat(matIds[i]);
//
// 		pPostProcMat->SetPostProcessVarIndexList(varNames,pmat);
//
// 		matNumber = pcPostProcMesh->InsertMaterialObject(pPostProcMat);
//
// 	}
//
// 	AutoBuildDisc();
//
//     pcPostProcMesh->ExpandSolution();
// }
//
// void TPZPostProcAnalysis::AutoBuildDisc()
// {
//     	TPZAdmChunkVector<TPZGeoEl *> &elvec = Mesh()->Reference()->ElementVec();
// 	int64_t i, nelem = elvec.NElements();
// 	int neltocreate = 0;
//
//     // build a data structure indicating which geometric elements will be post processed
//     fpMainMesh->LoadReferences();
//     std::map<TPZGeoEl *,TPZCompEl *> geltocreate;
//     TPZCompMesh * pcPostProcMesh = this->Mesh();
//     for (i=0; i<nelem; i++) {
//         TPZGeoEl * gel = elvec[i];
//         if (!gel) {
//             continue;
//         }
//
//         if (gel->HasSubElement()) {
//             //continue;
//         }
//
//         TPZMaterial * mat = pcPostProcMesh->FindMaterial(gel->MaterialId());
//         if(!mat)
//         {
//             //continue;
//         }
//
//         if (gel->Reference()) {
//             geltocreate[elvec[i]] = gel->Reference();
//         }
//     }
//     Mesh()->Reference()->ResetReference();
//     Mesh()->LoadReferences();
//     neltocreate = geltocreate.size();
//
// 	std::set<int> matnotfound;
// 	int nbl = Mesh()->Block().NBlocks();
// 	if(neltocreate > nbl) Mesh()->Block().SetNBlocks(neltocreate);
// 	Mesh()->Block().SetNBlocks(nbl);
//
//     std::map<TPZGeoEl *, TPZCompEl *>::iterator it;
//     for (it=geltocreate.begin(); it!= geltocreate.end(); it++)
//     {
// 		TPZGeoEl *gel = it->first;
// 		if(!gel) continue;
//         int matid = gel->MaterialId();
//         TPZMaterial * mat = Mesh()->FindMaterial(matid);
//         if(!mat)
//         {
//             matnotfound.insert(matid);
//             continue;
//         }
//         TPZCompEl *cel = Mesh()->CreateCompEl(gel);
//         TPZCompElPostProcBase *celpost = dynamic_cast<TPZCompElPostProcBase *>(cel);
//         if(!celpost) DebugStop();
//         TPZCompEl *celref = it->second;
//         int nc = cel->NConnects();
//         int ncref = celref->NConnects();
//         TPZInterpolationSpace *celspace = dynamic_cast<TPZInterpolationSpace *>(cel);
//         TPZInterpolationSpace *celrefspace = dynamic_cast<TPZInterpolationSpace *>(celref);
//         int porder;
//         if (!celrefspace) {
//             TPZMultiphysicsElement *celrefmf = dynamic_cast<TPZMultiphysicsElement *>(celref);
//             if (celrefmf){
//                 celrefspace = dynamic_cast<TPZInterpolationSpace *>(celrefmf->Element(0));
//             } else {
//                 DebugStop();
//             }
//         }
//         celpost->fReferredElement = celrefspace;
//
//         if (celrefspace) {
//             porder = celrefspace->GetPreferredOrder();
//         } else {
//             DebugStop();
//         }
//
//         celspace->SetPreferredOrder(porder);
//         for (int ic=0; ic<nc; ic++) {
//             cel->Connect(ic).SetOrder(porder,cel->ConnectIndex(ic));
//             int nshape = celspace->NConnectShapeF(ic,porder);
//             cel->Connect(ic).SetNShape(nshape);
//         }
//
//         TPZIntPoints &intrule = celspace->GetIntegrationRule();
//         const TPZIntPoints &intruleref = celref->GetIntegrationRule();
//         TPZIntPoints * cloned_rule = intruleref.Clone();
//         cel->SetIntegrationRule(cloned_rule);
//
// #ifdef PZDEBUG
//         if (cel->GetIntegrationRule().NPoints() != intruleref.NPoints()) {
//             DebugStop();
//         }
// #endif
//         // this is why the mesh will be discontinuous!!
//         gel->ResetReference();
//
// 	}
//
//     // we changed the properties of the connects
//     // now synchronize the connect properties with the block sizes
// 	int64_t nc= Mesh()->NConnects();
//     for (int64_t ic=0; ic<nc; ic++) {
//         TPZConnect &c = Mesh()->ConnectVec()[ic];
//         int blsize = c.NShape()*c.NState();
//         int64_t seqnum = c.SequenceNumber();
//         Mesh()->Block().Set(seqnum, blsize);
//     }
// 	Mesh()->InitializeBlock();
// #ifdef PZ_LOG
//   //  if(PPAnalysisLogger.isDebugEnabled())
//     {
//         std::stringstream sout;
//         Mesh()->Print(sout);
//        // LOGPZ_DEBUG(PPAnalysisLogger, sout.str())
//     }
// #endif
//
// //#ifdef PZDEBUG
// 	if(matnotfound.size())
// 	{
// 		std::cout << "Post-processing mesh was created without these materials: ";
// 		std::set<int>::iterator it;
// 		for(it = matnotfound.begin(); it!= matnotfound.end(); it++)
// 		{
// 			std::cout << *it << " ";
// 		}
// 		std::cout << std::endl;
//      //   DebugStop();
// 	}
// //#endif
// // 	TPZAdmChunkVector<TPZGeoEl *> &elvec = Mesh()->Reference()->ElementVec();
// // 	long i, nelem = elvec.NElements();
// // 	int neltocreate = 0;
// // 	long index;
// //     // build a data structure indicating which geometric elements will be post processed
// //     fpMainMesh->LoadReferences();
// //     std::map<TPZGeoEl *,TPZCompEl *> geltocreate;
// //     for (i=0; i<nelem; i++) {
// //         if (!elvec[i]) {
// //             continue;
// //         }
// //         if (elvec[i]->Reference()) {
// //             geltocreate[elvec[i]] = elvec[i]->Reference();
// //         }
// //     }
// //     Mesh()->Reference()->ResetReference();
// //     Mesh()->LoadReferences();
// //     neltocreate = geltocreate.size();
// //
// // 	std::set<int> matnotfound;
// // 	int nbl = Mesh()->Block().NBlocks();
// // 	if(neltocreate > nbl) Mesh()->Block().SetNBlocks(neltocreate);
// // 	Mesh()->Block().SetNBlocks(nbl);
// //
// //     std::map<TPZGeoEl *, TPZCompEl *>::iterator it;
// //     for (it=geltocreate.begin(); it!= geltocreate.end(); it++)
// //     {
// // 		TPZGeoEl *gel = it->first;
// // 		if(!gel) continue;
// //         int matid = gel->MaterialId();
// //         TPZMaterial * mat = Mesh()->FindMaterial(matid);
// //         if(!mat)
// //         {
// //             matnotfound.insert(matid);
// //             continue;
// //         }
// //         int printing = 0;
// //         if (printing) {
// //             gel->Print(cout);
// //         }
// //
// //
// //         TPZCompEl *cel =  Mesh()->CreateCompEl(gel);
// //         TPZCompEl *celref = it->second;
// //         int nc = cel->NConnects();
// //         int ncref = celref->NConnects();
// //         if (nc != ncref) {
// //             DebugStop();
// //         }
// //         TPZInterpolationSpace *celspace = dynamic_cast<TPZInterpolationSpace *>(cel);
// //         TPZInterpolationSpace *celrefspace = dynamic_cast<TPZInterpolationSpace *>(celref);
// //         int porder = celrefspace->GetPreferredOrder();
// // //        if (porder != 2) {
// // //            std::cout << "I should stop porder = " << porder << std::endl;
// // //        }
// //         celspace->SetPreferredOrder(porder);
// //         for (int ic=0; ic<nc; ic++) {
// //             cel->Connect(ic).SetOrder(porder,cel->ConnectIndex(ic));
// //             int nshape = celspace->NConnectShapeF(ic,porder);
// //             cel->Connect(ic).SetNShape(nshape);
// //         }
// //         TPZIntPoints &intrule = celspace->GetIntegrationRule();
// //         TPZVec<int> intorder(gel->Dimension(),0);
// //         const TPZIntPoints &intruleref = celrefspace->GetIntegrationRule();
// //         intruleref.GetOrder(intorder);
// //         intrule.SetOrder(intorder);
// // #ifdef DEBUG
// //         if (intrule.NPoints() != intruleref.NPoints()) {
// //             DebugStop();
// //         }
// // #endif
// //         gel->ResetReference();
// //
// // 	}
// //
// //     // we changed the properties of the connects
// //     // now synchronize the connect properties with the block sizes
// // 	long nc= Mesh()->NConnects();
// //     for (long ic=0; ic<nc; ic++) {
// //         TPZConnect &c = Mesh()->ConnectVec()[ic];
// //         int blsize = c.NShape()*c.NState();
// //         long seqnum = c.SequenceNumber();
// //         Mesh()->Block().Set(seqnum, blsize);
// //     }
// // 	Mesh()->InitializeBlock();
// // #ifdef LOG4CXX
// //     if(PPAnalysisLogger->isDebugEnabled())
// //     {
// //         std::stringstream sout;
// //         Mesh()->Print(sout);
// //         LOGPZ_DEBUG(PPAnalysisLogger, sout.str())
// //     }
// // #endif
// // 	if(matnotfound.size())
// // 	{
// // 		std::cout << "Malha post proc was created without these materials ";
// // 		std::set<int>::iterator it;
// // 		for(it = matnotfound.begin(); it!= matnotfound.end(); it++)
// // 		{
// // 			std::cout << *it << " ";
// // 		}
// // 		std::cout << std::endl;
// // 	}
//
//
// }
//
// void TPZPostProcAnalysis::Assemble()
// {
//    PZError << "Error at " << __PRETTY_FUNCTION__ << " TPZPostProcAnalysis::Assemble() should never be called\n";
// }
//
// void TPZPostProcAnalysis::Solve(){
//    PZError << "Error at " << __PRETTY_FUNCTION__ << " TPZPostProcAnalysis::Solve() should never be called\n";
// }
//
// void TPZPostProcAnalysis::TransferSolution()
// {
//  // this is where we compute the projection of the post processed variables
//     TPZLinearAnalysis::AssembleResidual();
//     fSolution = Rhs();
//     TPZLinearAnalysis::LoadSolution();
//
//     TPZCompMesh *compmeshPostProcess = (Mesh());
//     if (!compmeshPostProcess) {
//         DebugStop();
//     }
//     // fpMainMesh is the mesh with the actual finite element approximation, but probably stored at
//     // integration points
//     TPZCompMesh *solmesh = fpMainMesh;
//
//     fpMainMesh->Reference()->ResetReference();
//     fpMainMesh->LoadReferences();
//     //In case the post processing computed element solutions
//
//     // copy the values from the post processing mesh to the finite element mesh
//     TPZFMatrix<STATE> &comprefElSol = compmeshPostProcess->ElementSolution();
//     // solmesh if the finite element simulation mesh
//     const TPZFMatrix<STATE> &solmeshElSol = solmesh->ElementSolution();
//
//
//     int64_t numelsol = solmesh->ElementSolution().Rows();
//     int64_t nelem = compmeshPostProcess->NElements();
//     compmeshPostProcess->ElementSolution().Redim(nelem, numelsol);
//     if (numelsol)//caso tenha solucao no elemento
//     {
//      //   DebugStop();
//         for (int64_t el=0; el<nelem; el++) {
//             TPZCompEl *celpost = compmeshPostProcess->Element(el);
//             TPZGeoEl *gel = celpost->Reference();
//             // we dont acount for condensed elements submeshes etc
//             if(!gel) DebugStop();
//             TPZCompEl *cel = gel->Reference();
//             if (!cel) {
//                 DebugStop();
//             }
//             int64_t index = cel->Index();
//             // we copy from the simulation mesh to the post processing mesh
//             for (int64_t isol=0; isol<numelsol; isol++) {
//                 comprefElSol(el,isol) = solmeshElSol.Get(index,isol);
//             }
//         }
//     }
//
// }
//
//
// void TPZPostProcAnalysis::SetAllCreateFunctionsPostProc(TPZCompMesh *cmesh)
// {
//
//     TPZManVector<TCreateFunction,10> functions(8);
//
//     functions[EPoint] = &TPZPostProcAnalysis::CreatePointEl;
//     functions[EOned] = TPZPostProcAnalysis::CreateLinearEl;
//     functions[EQuadrilateral] = TPZPostProcAnalysis::CreateQuadEl;
//     functions[ETriangle] = TPZPostProcAnalysis::CreateTriangleEl;
//     functions[EPrisma] = TPZPostProcAnalysis::CreatePrismEl;
//     functions[ETetraedro] = TPZPostProcAnalysis::CreateTetraEl;
//     functions[EPiramide] = TPZPostProcAnalysis::CreatePyramEl;
//     functions[ECube] = TPZPostProcAnalysis::CreateCubeEl;
//     cmesh->ApproxSpace().SetCreateFunctions(functions);
// }
//
//
// #include "TPZCompElH1.h"
//
// using namespace pzshape;
//
// template class TPZCompElPostProc< TPZCompElH1<TPZShapePoint> >;
// template class TPZCompElPostProc< TPZCompElH1<TPZShapeLinear> >;
// template class TPZCompElPostProc< TPZCompElH1<TPZShapeQuad> >;
// template class TPZCompElPostProc< TPZCompElH1<TPZShapeTriang> >;
// template class TPZCompElPostProc< TPZCompElH1<TPZShapeCube> >;
// template class TPZCompElPostProc< TPZCompElH1<TPZShapePrism> >;
// template class TPZCompElPostProc< TPZCompElH1<TPZShapePiram> >;
// template class TPZCompElPostProc< TPZCompElH1<TPZShapeTetra> >;
// template class TPZCompElPostProc< TPZCompElDisc >;
//
// template class TPZRestoreClass<TPZCompElPostProc< TPZCompElH1<TPZShapePoint> >>;
// template class TPZRestoreClass<TPZCompElPostProc< TPZCompElH1<TPZShapeLinear> >>;
// template class TPZRestoreClass<TPZCompElPostProc< TPZCompElH1<TPZShapeQuad> >>;
// template class TPZRestoreClass<TPZCompElPostProc< TPZCompElH1<TPZShapeTriang> >>;
// template class TPZRestoreClass<TPZCompElPostProc< TPZCompElH1<TPZShapeCube> >>;
// template class TPZRestoreClass<TPZCompElPostProc< TPZCompElH1<TPZShapePrism> >>;
// template class TPZRestoreClass<TPZCompElPostProc< TPZCompElH1<TPZShapePiram> >>;
// template class TPZRestoreClass<TPZCompElPostProc< TPZCompElH1<TPZShapeTetra> >>;
// template class TPZRestoreClass<TPZCompElPostProc< TPZCompElDisc >>;
//
// TPZCompEl *TPZPostProcAnalysis::CreatePointEl(TPZGeoEl *gel,TPZCompMesh &mesh) {
// 	if(!gel->Reference() && gel->NumInterfaces() == 0)
// 		return new TPZCompElPostProc< TPZCompElH1<TPZShapePoint> >(mesh,gel);
// 	return NULL;
// }
// TPZCompEl *TPZPostProcAnalysis::CreateLinearEl(TPZGeoEl *gel,TPZCompMesh &mesh) {
// 	if(!gel->Reference() && gel->NumInterfaces() == 0)
// 		return new TPZCompElPostProc<TPZCompElH1<TPZShapeLinear> >(mesh,gel);
// 	return NULL;
// }
// TPZCompEl *TPZPostProcAnalysis::CreateQuadEl(TPZGeoEl *gel,TPZCompMesh &mesh) {
// 	if(!gel->Reference() && gel->NumInterfaces() == 0)
// 		return new TPZCompElPostProc<TPZCompElH1<TPZShapeQuad> >(mesh,gel);
// 	return NULL;
// }
// TPZCompEl *TPZPostProcAnalysis::CreateTriangleEl(TPZGeoEl *gel,TPZCompMesh &mesh) {
// 	if(!gel->Reference() && gel->NumInterfaces() == 0)
// 		return new TPZCompElPostProc<TPZCompElH1<TPZShapeTriang> >(mesh,gel);
// 	return NULL;
// }
// TPZCompEl *TPZPostProcAnalysis::CreateCubeEl(TPZGeoEl *gel,TPZCompMesh &mesh) {
// 	if(!gel->Reference() && gel->NumInterfaces() == 0)
// 		return new TPZCompElPostProc<TPZCompElH1<TPZShapeCube> >(mesh,gel);
// 	return NULL;
// }
// TPZCompEl *TPZPostProcAnalysis::CreatePrismEl(TPZGeoEl *gel,TPZCompMesh &mesh) {
// 	if(!gel->Reference() && gel->NumInterfaces() == 0)
// 		return new TPZCompElPostProc< TPZCompElH1<TPZShapePrism> >(mesh,gel);
// 	return NULL;
// }
// TPZCompEl *TPZPostProcAnalysis::CreatePyramEl(TPZGeoEl *gel,TPZCompMesh &mesh) {
// 	if(!gel->Reference() && gel->NumInterfaces() == 0)
// 		return new TPZCompElPostProc<TPZCompElH1<TPZShapePiram> >(mesh,gel);
// 	return NULL;
// }
// TPZCompEl *TPZPostProcAnalysis::CreateTetraEl(TPZGeoEl *gel,TPZCompMesh &mesh) {
// 	if(!gel->Reference() && gel->NumInterfaces() == 0)
// 		return new TPZCompElPostProc<TPZCompElH1<TPZShapeTetra> >(mesh,gel);
// 	return NULL;
// }
//
//
// TPZCompEl * TPZPostProcAnalysis::CreatePostProcDisc(TPZGeoEl *gel, TPZCompMesh &mesh)
// {
// 	return new TPZCompElPostProc< TPZCompElDisc > (mesh,gel);
// }
//
// /** @brief Returns the unique identifier for reading/writing objects to streams */
// int TPZPostProcAnalysis::ClassId() const{
//     return Hash("TPZPostProcAnalysis") ^ TPZLinearAnalysis::ClassId() << 1;
// }
// /** @brief Save the element data to a stream */
// void TPZPostProcAnalysis::Write(TPZStream &buf, int withclassid) const
// {
//     TPZLinearAnalysis::Write(buf, withclassid);
//     TPZPersistenceManager::WritePointer(fpMainMesh, &buf);
// }
//
// /** @brief Read the element data from a stream */
// void TPZPostProcAnalysis::Read(TPZStream &buf, void *context)
// {
//     TPZLinearAnalysis::Read(buf, context);
//     fpMainMesh = dynamic_cast<TPZCompMesh*>(TPZPersistenceManager::GetInstance(&buf));
// }
//$Id: pzpostprocanalysis.cpp,v 1.10 2010-11-23 18:58:35 diogo Exp $
#include "TPZLinearAnalysis.h"
#include "pzpostprocanalysis.h"
#include "pzpostprocmat.h"
#include "pzcompelpostproc.h"
#include "pzcmesh.h"
#include "pzgmesh.h"
#include "pzvec.h"
#include "tpzautopointer.h"

#include "pzstring.h"
//#include "pzelastoplasticanalysis.h"
#include "pzcreateapproxspace.h"
#include "pzmultiphysicselement.h"

#include <map>
#include <set>
#include <stdio.h>
#include "pzlog.h"

#include "pzintel.h"

#include "pzrefpoint.h"
#include "pzgeopoint.h"
#include "pzshapepoint.h"
#include "tpzpoint.h"

#include "pzshapelinear.h"
#include "TPZGeoLinear.h"
#include "TPZRefLinear.h"
#include "tpzline.h"

#include "pzshapetriang.h"
#include "pzreftriangle.h"
#include "pzgeotriangle.h"
#include "tpztriangle.h"

#include "pzrefquad.h"
#include "pzshapequad.h"
#include "pzgeoquad.h"
#include "tpzquadrilateral.h"

#include "pzshapeprism.h"
#include "pzrefprism.h"
#include "pzgeoprism.h"
#include "tpzprism.h"

#include "pzshapetetra.h"
#include "pzreftetrahedra.h"
#include "pzgeotetrahedra.h"
#include "tpztetrahedron.h"

#include "pzshapepiram.h"
#include "pzrefpyram.h"
#include "pzgeopyramid.h"
#include "tpzpyramid.h"

#include "TPZGeoCube.h"
#include "pzshapecube.h"
#include "TPZRefCube.h"
#include "tpzcube.h"
#include "pzelctemp.h"
#include "TPZMatCombinedSpaces.h"

#include <algorithm>
#include <cmath>
#include <memory>
#include <vector>

#ifdef PZ_LOG
static TPZLogger PPAnalysisLogger ( "pz.analysis.postproc" );
#endif

using namespace std;

TPZPostProcAnalysis::TPZPostProcAnalysis() : TPZRegisterClassId ( &TPZPostProcAnalysis::ClassId ),
        fpMainMesh ( NULL )
{
}

TPZPostProcAnalysis::TPZPostProcAnalysis ( TPZCompMesh * pRef ) :TPZRegisterClassId ( &TPZPostProcAnalysis::ClassId ),
        TPZLinearAnalysis(), fpMainMesh ( pRef )
{

        SetCompMesh ( pRef );

}

TPZPostProcAnalysis::TPZPostProcAnalysis ( const TPZPostProcAnalysis &copy ) : TPZRegisterClassId ( &TPZPostProcAnalysis::ClassId ),
        TPZLinearAnalysis ( copy ), fpMainMesh ( 0 )
{

}

TPZPostProcAnalysis &TPZPostProcAnalysis::operator= ( const TPZPostProcAnalysis &copy )
{
        SetCompMesh ( 0 );
        return *this;
}

TPZPostProcAnalysis::~TPZPostProcAnalysis()
{
        if ( fCompMesh ) {
                delete fCompMesh;
        }

}

/// Set the computational mesh we are going to post process
void TPZPostProcAnalysis::SetCompMesh ( TPZCompMesh *pRef, bool mustOptimizeBandwidth )
{
        // the postprocess mesh already exists, do nothing
        if ( fpMainMesh == pRef ) {
                return;
        }

        if ( fCompMesh ) {
                delete fCompMesh;
                fCompMesh = 0;
                TPZLinearAnalysis::CleanUp();
        }

        fpMainMesh = pRef;

        if ( !pRef ) {
                return;
        }

        TPZCompMesh* pcMainMesh = fpMainMesh;

        TPZGeoMesh * pgmesh = pcMainMesh->Reference();

        TPZCompMesh * pcPostProcMesh = new TPZCompMesh ( pgmesh );

        fCompMesh = pcPostProcMesh;

        TPZPostProcAnalysis::SetAllCreateFunctionsPostProc ( pcPostProcMesh );

}


void TPZPostProcAnalysis::SetPostProcessVariables ( TPZVec<int> & matIds, TPZVec<std::string> &varNames )
{

        int nMat, matNumber;
        TPZCompMesh * pcMainMesh = fpMainMesh;

        TPZCompMesh * pcPostProcMesh = ( this->Mesh() );

        if ( !pcPostProcMesh ) {
                DebugStop();
        }

        nMat = matIds.NElements();
        for ( int i = 0; i < nMat; i++ ) {
                TPZMaterial * pmat = pcMainMesh->FindMaterial ( matIds[i] );
                if ( !pmat ) {
                        PZError << "Error at " << __PRETTY_FUNCTION__ << " TPZPostProcAnalysis::SetPostProcessVariables() material Id " << matIds[i] << " not found in original mesh!\n";
                        continue;
                }

                TPZPostProcMat * pPostProcMat = new TPZPostProcMat ( matIds[i] );

                pPostProcMat->SetPostProcessVarIndexList ( varNames,pmat );

                matNumber = pcPostProcMesh->InsertMaterialObject ( pPostProcMat );

        }

        AutoBuildDisc();

        pcPostProcMesh->ExpandSolution();
}

/**
 * @brief Verifies that the integration rule of a multiphysics element with memory (GetIntegrationRule, cloned
 * by the post-processing element) is the rule of its assembly and of its memory
 *
 * TPZMultiphysicsCompEl::CalcStiff and CalcResidual integrate with the rule of the interior side of the
 * geometric element whose order is chosen by the material (IntegrationRuleOrder of the maximum orders of the
 * atomic spaces), and the material reads the memory item intGlobPtIndex = memory index of the point
 * intLocPtIndex of that rule. The post-processing must therefore use the same points in the same order;
 * the memory indices must be one per point. An element without memory imposes no condition: its variables
 * are computed from the solution, at the points of any rule.
 */
static void CheckMultiphysicsIntegrationRule ( TPZMultiphysicsElement *celmf )
{
        TPZManVector<int64_t,64> memindices;
        celmf->GetMemoryIndices ( memindices );
        if ( memindices.size() == 0 ) return;
        auto *mat = dynamic_cast<TPZMatCombinedSpaces *> ( celmf->Material() );
        TPZGeoEl *gel = celmf->Reference();
        if ( !mat || !gel ) {
                PZError << "Error at " << __PRETTY_FUNCTION__ << " multiphysics element without combined-space material\n";
                DebugStop();
        }
        TPZManVector<int,4> ordervec;
        for ( int64_t iref = 0; iref < celmf->NMeshes(); iref++ ) {
                TPZInterpolationSpace *msp = dynamic_cast<TPZInterpolationSpace *> ( celmf->Element ( iref ) );
                if ( !msp ) continue;
                ordervec.Resize ( ordervec.size() +1 );
                ordervec[ordervec.size()-1] = msp->MaxOrder();
        }
        const int order = mat->IntegrationRuleOrder ( ordervec );
        std::unique_ptr<TPZIntPoints> assemblyrule ( gel->CreateSideIntegrationRule ( gel->NSides()-1, order ) );
        TPZManVector<int,3> orderdim ( gel->Dimension(), order );
        assemblyrule->SetOrder ( orderdim );
        const TPZIntPoints &elrule = celmf->GetIntegrationRule();
        bool same = assemblyrule->NPoints() == elrule.NPoints();
        TPZManVector<REAL,3> pa ( gel->Dimension(),0. ), pb ( gel->Dimension(),0. );
        REAL wa, wb;
        for ( int ip = 0; same && ip < elrule.NPoints(); ip++ ) {
                assemblyrule->Point ( ip, pa, wa );
                elrule.Point ( ip, pb, wb );
                REAL diff = std::fabs ( wa-wb );
                for ( int d = 0; d < gel->Dimension(); d++ ) diff += std::fabs ( pa[d]-pb[d] );
                if ( diff > 1.e-12 ) same = false;
        }
        if ( !same || memindices.size() != elrule.NPoints() ) {
                PZError << "Error at " << __PRETTY_FUNCTION__ << " element " << celmf->Index()
                        << ": the integration rule of the multiphysics element (" << elrule.NPoints()
                        << " points) is not the rule of its assembly (" << assemblyrule->NPoints()
                        << " points) or of its memory (" << memindices.size() << " items)\n";
                DebugStop();
        }
}

void TPZPostProcAnalysis::AutoBuildDisc()
{
        TPZAdmChunkVector<TPZGeoEl *> &elvec = Mesh()->Reference()->ElementVec();
        int64_t i, nelem = elvec.NElements();
        int neltocreate = 0;

        // build a data structure indicating which geometric elements will be post processed
        // (in the order of the geometric mesh, so that the post-processing elements and the graphical
        // output follow that order)
        fpMainMesh->LoadReferences();
        std::vector<std::pair<TPZGeoEl *,TPZCompEl *> > geltocreate;
        TPZCompMesh * pcPostProcMesh = this->Mesh();
        for ( i=0; i<nelem; i++ ) {
                TPZGeoEl * gel = elvec[i];
                if ( !gel ) {
                        continue;
                }

                if ( gel->HasSubElement() ) {
                        continue;
                }

                TPZMaterial * mat = pcPostProcMesh->FindMaterial ( gel->MaterialId() );
                if ( !mat ) {
                        continue;
                }

                if ( gel->Reference() ) {
                        geltocreate.push_back ( std::make_pair ( elvec[i], gel->Reference() ) );
                }
        }
        Mesh()->Reference()->ResetReference();
        Mesh()->LoadReferences();
        neltocreate = geltocreate.size();

        std::set<int> matnotfound;
        int nbl = Mesh()->Block().NBlocks();
        if ( neltocreate > nbl ) Mesh()->Block().SetNBlocks ( neltocreate );
        Mesh()->Block().SetNBlocks ( nbl );

        for ( auto it=geltocreate.begin(); it!= geltocreate.end(); it++ ) {
                TPZGeoEl *gel = it->first;
                if ( !gel ) continue;
                int matid = gel->MaterialId();
                TPZMaterial * mat = Mesh()->FindMaterial ( matid );
                if ( !mat ) {
                        matnotfound.insert ( matid );
                        continue;
                }
                TPZCompEl *cel = Mesh()->CreateCompEl ( gel );
                TPZCompElPostProcBase *celpost = dynamic_cast<TPZCompElPostProcBase *> ( cel );
                if ( !celpost ) DebugStop();
                TPZCompEl *celref = it->second;
                int nc = cel->NConnects();
                TPZInterpolationSpace *celspace = dynamic_cast<TPZInterpolationSpace *> ( cel );
                TPZInterpolationSpace *celrefspace = dynamic_cast<TPZInterpolationSpace *> ( celref );
                TPZMultiphysicsElement *celrefmf = dynamic_cast<TPZMultiphysicsElement *> ( celref );
                int porder = -1;
                if ( celrefspace ) {
                        porder = celrefspace->GetPreferredOrder();
                } else if ( celrefmf ) {
                        // the multiphysics element itself is referred: its material (a combined-space material,
                        // possibly with memory) evaluates the variables from the data of all the atomic spaces.
                        // The order of the projection starts from the highest order of the atomic spaces.
                        for ( int64_t iref = 0; iref < celrefmf->NMeshes(); iref++ ) {
                                TPZInterpolationSpace *msp = dynamic_cast<TPZInterpolationSpace *> ( celrefmf->Element ( iref ) );
                                if ( msp ) porder = std::max ( porder, msp->GetPreferredOrder() );
                        }
                        CheckMultiphysicsIntegrationRule ( celrefmf );
                } else {
                        DebugStop();
                }
                if ( porder < 1 ) DebugStop();
                celpost->fReferredElement = celref;

                const TPZIntPoints &intruleref = celref->GetIntegrationRule();
                // The element-wise L2 projection (TPZCompElPostProc::CalcResidual) is singular when the element
                // has more shape functions than integration points (e.g. quadratic elements with 2 x 2 or
                // 2 x 2 x 2 points): the order is reduced until the number of shape functions does not exceed
                // the number of points. With n x n (x n) Gauss points the order n-1 gives the full tensor
                // product space, and the projection is the Lagrange interpolation of the values at the points.
                auto nshapeorder = [&] ( int order ) {
                        int nshape = 0;
                        for ( int ic=0; ic<nc; ic++ ) nshape += celspace->NConnectShapeF ( ic,order );
                        return nshape;
                };
                while ( porder > 1 && nshapeorder ( porder ) > intruleref.NPoints() ) porder--;

                celspace->SetPreferredOrder ( porder );
                for ( int ic=0; ic<nc; ic++ ) {
                        cel->Connect ( ic ).SetOrder ( porder,cel->ConnectIndex ( ic ) );
                        int nshape = celspace->NConnectShapeF ( ic,porder );
                        cel->Connect ( ic ).SetNShape ( nshape );
                }

                TPZIntPoints * cloned_rule = intruleref.Clone();
                cel->SetIntegrationRule ( cloned_rule );

#ifdef PZDEBUG
                if ( cel->GetIntegrationRule().NPoints() != intruleref.NPoints() ) {
                        DebugStop();
                }
#endif
                // this is why the mesh will be discontinuous!!
                gel->ResetReference();

        }

        // we changed the properties of the connects
        // now synchronize the connect properties with the block sizes
        int64_t nc= Mesh()->NConnects();
        for ( int64_t ic=0; ic<nc; ic++ ) {
                TPZConnect &c = Mesh()->ConnectVec() [ic];
                int blsize = c.NShape() *c.NState();
                int64_t seqnum = c.SequenceNumber();
                Mesh()->Block().Set ( seqnum, blsize );
        }
        Mesh()->InitializeBlock();
#ifdef PZ_LOG
        if ( PPAnalysisLogger.isDebugEnabled() ) {
                std::stringstream sout;
                Mesh()->Print ( sout );
                LOGPZ_DEBUG ( PPAnalysisLogger, sout.str() )
        }
#endif

#ifdef PZDEBUG
        if ( matnotfound.size() ) {
                std::cout << "Post-processing mesh was created without these materials: ";
                std::set<int>::iterator it;
                for ( it = matnotfound.begin(); it!= matnotfound.end(); it++ ) {
                        std::cout << *it << " ";
                }
                std::cout << std::endl;
                DebugStop();
        }
#endif

}

void TPZPostProcAnalysis::Assemble()
{
        PZError << "Error at " << __PRETTY_FUNCTION__ << " TPZPostProcAnalysis::Assemble() should never be called\n";
}

void TPZPostProcAnalysis::Solve()
{
        PZError << "Error at " << __PRETTY_FUNCTION__ << " TPZPostProcAnalysis::Solve() should never be called\n";
}

void TPZPostProcAnalysis::TransferSolution()
{

        // this is where we compute the projection of the post processed variables
        TPZLinearAnalysis::AssembleResidual();
        fSolution = Rhs();
        TPZLinearAnalysis::LoadSolution();

        TPZCompMesh *compmeshPostProcess = ( Mesh() );
        if ( !compmeshPostProcess ) {
                DebugStop();
        }
        // fpMainMesh is the mesh with the actual finite element approximation, but probably stored at
        // integration points
        TPZCompMesh *solmesh = fpMainMesh;
        fpMainMesh->Reference()->ResetReference();
        fpMainMesh->LoadReferences();
        //In case the post processing computed element solutions
        // copy the values from the post processing mesh to the finite element mesh
        TPZFMatrix<STATE> &comprefElSol = compmeshPostProcess->ElementSolution();
        // solmesh if the finite element simulation mesh
        const TPZFMatrix<STATE> &solmeshElSol = solmesh->ElementSolution();
        int64_t numelsol = solmesh->ElementSolution().Cols();
        int64_t nelem = compmeshPostProcess->NElements();
        compmeshPostProcess->ElementSolution().Redim ( nelem, numelsol );
//         if ( numelsol ) {
//                 for ( int64_t el=0; el<nelem; el++ ) {
//                         TPZCompEl *celpost = compmeshPostProcess->Element ( el );
//                         TPZGeoEl *gel = celpost->Reference();
//                         // we dont acount for condensed elements submeshes etc
//                         if ( !gel ) DebugStop();
//                         TPZCompEl *cel = gel->Reference();
//                         if ( !cel ) {
//                                 DebugStop();
//                         }
//                         int64_t index = cel->Index();
//                         // we copy from the simulation mesh to the post processing mesh
//                         for ( int64_t isol=0; isol<numelsol; isol++ ) {
//                                 comprefElSol ( el,isol ) = solmeshElSol.Get ( index,isol );
//                         }
//                 }
//         }
}


void TPZPostProcAnalysis::SetAllCreateFunctionsPostProc ( TPZCompMesh *cmesh )
{

        TPZManVector<TCreateFunction,10> functions ( 8 );

        functions[EPoint] = &TPZPostProcAnalysis::CreatePointEl;
        functions[EOned] = TPZPostProcAnalysis::CreateLinearEl;
        functions[EQuadrilateral] = TPZPostProcAnalysis::CreateQuadEl;
        functions[ETriangle] = TPZPostProcAnalysis::CreateTriangleEl;
        functions[EPrisma] = TPZPostProcAnalysis::CreatePrismEl;
        functions[ETetraedro] = TPZPostProcAnalysis::CreateTetraEl;
        functions[EPiramide] = TPZPostProcAnalysis::CreatePyramEl;
        functions[ECube] = TPZPostProcAnalysis::CreateCubeEl;
        cmesh->ApproxSpace().SetCreateFunctions ( functions );
}


#include "TPZCompElH1.h"

using namespace pzshape;

template class TPZCompElPostProc< TPZCompElH1<TPZShapePoint> >;
template class TPZCompElPostProc< TPZCompElH1<TPZShapeLinear> >;
template class TPZCompElPostProc< TPZCompElH1<TPZShapeQuad> >;
template class TPZCompElPostProc< TPZCompElH1<TPZShapeTriang> >;
template class TPZCompElPostProc< TPZCompElH1<TPZShapeCube> >;
template class TPZCompElPostProc< TPZCompElH1<TPZShapePrism> >;
template class TPZCompElPostProc< TPZCompElH1<TPZShapePiram> >;
template class TPZCompElPostProc< TPZCompElH1<TPZShapeTetra> >;
template class TPZCompElPostProc< TPZCompElDisc >;

template class TPZRestoreClass<TPZCompElPostProc< TPZCompElH1<TPZShapePoint> >>;
template class TPZRestoreClass<TPZCompElPostProc< TPZCompElH1<TPZShapeLinear> >>;
template class TPZRestoreClass<TPZCompElPostProc< TPZCompElH1<TPZShapeQuad> >>;
template class TPZRestoreClass<TPZCompElPostProc< TPZCompElH1<TPZShapeTriang> >>;
template class TPZRestoreClass<TPZCompElPostProc< TPZCompElH1<TPZShapeCube> >>;
template class TPZRestoreClass<TPZCompElPostProc< TPZCompElH1<TPZShapePrism> >>;
template class TPZRestoreClass<TPZCompElPostProc< TPZCompElH1<TPZShapePiram> >>;
template class TPZRestoreClass<TPZCompElPostProc< TPZCompElH1<TPZShapeTetra> >>;
template class TPZRestoreClass<TPZCompElPostProc< TPZCompElDisc >>;

TPZCompEl *TPZPostProcAnalysis::CreatePointEl ( TPZGeoEl *gel,TPZCompMesh &mesh )
{
        if ( !gel->Reference() && gel->NumInterfaces() == 0 )
                return new TPZCompElPostProc< TPZCompElH1<TPZShapePoint> > ( mesh,gel );
        return NULL;
}
TPZCompEl *TPZPostProcAnalysis::CreateLinearEl ( TPZGeoEl *gel,TPZCompMesh &mesh )
{
        if ( !gel->Reference() && gel->NumInterfaces() == 0 )
                return new TPZCompElPostProc<TPZCompElH1<TPZShapeLinear> > ( mesh,gel );
        return NULL;
}
TPZCompEl *TPZPostProcAnalysis::CreateQuadEl ( TPZGeoEl *gel,TPZCompMesh &mesh )
{
        if ( !gel->Reference() && gel->NumInterfaces() == 0 )
                return new TPZCompElPostProc<TPZCompElH1<TPZShapeQuad> > ( mesh,gel );
        return NULL;
}
TPZCompEl *TPZPostProcAnalysis::CreateTriangleEl ( TPZGeoEl *gel,TPZCompMesh &mesh )
{
        if ( !gel->Reference() && gel->NumInterfaces() == 0 )
                return new TPZCompElPostProc<TPZCompElH1<TPZShapeTriang> > ( mesh,gel );
        return NULL;
}
TPZCompEl *TPZPostProcAnalysis::CreateCubeEl ( TPZGeoEl *gel,TPZCompMesh &mesh )
{
        if ( !gel->Reference() && gel->NumInterfaces() == 0 )
                return new TPZCompElPostProc<TPZCompElH1<TPZShapeCube> > ( mesh,gel );
        return NULL;
}
TPZCompEl *TPZPostProcAnalysis::CreatePrismEl ( TPZGeoEl *gel,TPZCompMesh &mesh )
{
        if ( !gel->Reference() && gel->NumInterfaces() == 0 )
                return new TPZCompElPostProc< TPZCompElH1<TPZShapePrism> > ( mesh,gel );
        return NULL;
}
TPZCompEl *TPZPostProcAnalysis::CreatePyramEl ( TPZGeoEl *gel,TPZCompMesh &mesh )
{
        if ( !gel->Reference() && gel->NumInterfaces() == 0 )
                return new TPZCompElPostProc<TPZCompElH1<TPZShapePiram> > ( mesh,gel );
        return NULL;
}
TPZCompEl *TPZPostProcAnalysis::CreateTetraEl ( TPZGeoEl *gel,TPZCompMesh &mesh )
{
        if ( !gel->Reference() && gel->NumInterfaces() == 0 )
                return new TPZCompElPostProc<TPZCompElH1<TPZShapeTetra> > ( mesh,gel );
        return NULL;
}


TPZCompEl * TPZPostProcAnalysis::CreatePostProcDisc ( TPZGeoEl *gel, TPZCompMesh &mesh )
{
        return new TPZCompElPostProc< TPZCompElDisc > ( mesh,gel );
}

/** @brief Returns the unique identifier for reading/writing objects to streams */
int TPZPostProcAnalysis::ClassId() const
{
        return Hash ( "TPZPostProcAnalysis" ) ^ TPZLinearAnalysis::ClassId() << 1;
}
/** @brief Save the element data to a stream */
void TPZPostProcAnalysis::Write ( TPZStream &buf, int withclassid ) const
{
        TPZLinearAnalysis::Write ( buf, withclassid );
        TPZPersistenceManager::WritePointer ( fpMainMesh, &buf );
}

/** @brief Read the element data from a stream */
void TPZPostProcAnalysis::Read ( TPZStream &buf, void *context )
{
        TPZLinearAnalysis::Read ( buf, context );
        fpMainMesh = dynamic_cast<TPZCompMesh*> ( TPZPersistenceManager::GetInstance ( &buf ) );
}
