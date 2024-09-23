//$Id: pzpostprocanalysis.cpp,v 1.10 2010-11-23 18:58:35 diogo Exp $
#include "pzpostprocanalysis.h"
#include <map>
#include <set>
#include <stdio.h>
#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr PPAnalysisLogger(Logger::getLogger("pz.analysis.postproc"));
#endif

using namespace std;

TPZPostProcAnalysis::TPZPostProcAnalysis() :TPZLinearAnalysis(),fpMainMesh(NULL)
{
}

TPZPostProcAnalysis::TPZPostProcAnalysis(TPZCompMesh * pRef):TPZLinearAnalysis(), fpMainMesh(pRef)
{

    SetCompMesh(pRef);

}

TPZPostProcAnalysis::TPZPostProcAnalysis(const TPZPostProcAnalysis &copy) : TPZLinearAnalysis(), fpMainMesh(0)
{

}

TPZPostProcAnalysis &TPZPostProcAnalysis::operator=(const TPZPostProcAnalysis &copy)
{
    SetCompMesh(0);
    return *this;
}

TPZPostProcAnalysis::~TPZPostProcAnalysis()
{
    if (fCompMesh) {
        delete fCompMesh;
    }

}

/// Set the computational mesh we are going to post process
void TPZPostProcAnalysis::SetCompMesh(TPZCompMesh *pRef)
{
    // the postprocess mesh already exists, do nothing
    if (fpMainMesh == pRef) {
        return;
    }

    if (fCompMesh) {
        std::cout << "PostProcAnalysis deleting the mesh " << (void *) fCompMesh << std::endl;
        delete fCompMesh;
        fCompMesh = 0;
        TPZAnalysis::CleanUp();
    }

    fpMainMesh = pRef;

    if (!pRef) {
        return;
    }

    TPZCompMesh* pcMainMesh = fpMainMesh;

    TPZGeoMesh * pgmesh = pcMainMesh->Reference();

    // TPZPostProcAnalysis::SetAllCreateFunctionsPostProc();


    TPZCompMeshReferred * pcPostProcMesh = new TPZCompMeshReferred(pgmesh);

    fCompMesh = pcPostProcMesh;

    TPZPostProcAnalysis::SetAllCreateFunctionsPostProc(pcPostProcMesh);

   // fpMainMesh->Solution().Print("SSSSSSSSSSS");

}


void TPZPostProcAnalysis::SetPostProcessVariables(TPZVec<int> & matIds, TPZVec<std::string> &varNames)
{
	//int j;
    int nMat, matNumber;


	TPZCompMesh * pcMainMesh = fpMainMesh;

	//TPZGeoMesh * pgmesh = pcMainMesh->Reference();

	TPZCompMeshReferred * pcPostProcMesh = dynamic_cast<TPZCompMeshReferred *>(this->Mesh());

    if (!pcPostProcMesh) {
        DebugStop();
    }

    if (pcPostProcMesh->ReferredMesh() == pcMainMesh) {
        return;
    }

    /*
	TPZStack<int> avlMatIds;
	long nel = pgmesh->NElements(), i;
	for(i = 0; i < nel; i++)
	{
		int matId = pgmesh->ElementVec()[i]->MaterialId();
		int isMatPostProc = 0;
		int isMatAvl = 0;
		j = 0;
		nMat = matIds.NElements();
		while(j < nMat && !isMatPostProc)
		{
			if(matId == matIds[j])isMatPostProc = 1;
			j++;
		}

		if(!isMatPostProc)
		{
			nMat = avlMatIds.NElements();
			j = 0;
			while(j < nMat && !isMatAvl)
			{
				if(matId == avlMatIds[j])isMatAvl = 1;
				j++;
			}

			if(!isMatAvl)
			{
				avlMatIds.Push(matId);
			}
		}
	}
	*/
	nMat = matIds.NElements();
	for(int i = 0; i < nMat; i++)
	{
		TPZMaterial * pmat = pcMainMesh->FindMaterial(matIds[i]);
		if(!pmat)
		{
			PZError << "Error at " << __PRETTY_FUNCTION__ << " TPZPostProcAnalysis::SetPostProcessVariables() material Id " << matIds[i] << " not found in original mesh!\n";
			continue;
		}

		TPZPostProcMat * pPostProcMat = new TPZPostProcMat(matIds[i]);

		pPostProcMat->SetPostProcessVarIndexList(varNames,pmat);

		matNumber = pcPostProcMesh->InsertMaterialObject(pPostProcMat);

	}

	//pcPostProcMesh->AutoBuild();
	//pcPostProcMesh->AutoBuildDisc();

	AutoBuildDisc();

	pcPostProcMesh->LoadReferred(pcMainMesh);
}

void TPZPostProcAnalysis::AutoBuildDisc()
{
    	TPZAdmChunkVector<TPZGeoEl *> &elvec = Mesh()->Reference()->ElementVec();
	int64_t i, nelem = elvec.NElements();
	int neltocreate = 0;

    // build a data structure indicating which geometric elements will be post processed
    fpMainMesh->LoadReferences();
    std::map<TPZGeoEl *,TPZCompEl *> geltocreate;
    TPZCompMesh * pcPostProcMesh = this->Mesh();
    for (i=0; i<nelem; i++) {
        TPZGeoEl * gel = elvec[i];
        if (!gel) {
            continue;
        }

        if (gel->HasSubElement()) {
            continue;
        }

        TPZMaterial * mat = pcPostProcMesh->FindMaterial(gel->MaterialId());
        if(!mat)
        {
            continue;
        }

        if (gel->Reference()) {
            geltocreate[elvec[i]] = gel->Reference();
        }
    }
    Mesh()->Reference()->ResetReference();
    Mesh()->LoadReferences();
    neltocreate = geltocreate.size();

	std::set<int> matnotfound;
	int nbl = Mesh()->Block().NBlocks();
	if(neltocreate > nbl) Mesh()->Block().SetNBlocks(neltocreate);
	Mesh()->Block().SetNBlocks(nbl);

    std::map<TPZGeoEl *, TPZCompEl *>::iterator it;
    for (it=geltocreate.begin(); it!= geltocreate.end(); it++)
    {
		TPZGeoEl *gel = it->first;
		if(!gel) continue;
        int matid = gel->MaterialId();
        TPZMaterial * mat = Mesh()->FindMaterial(matid);
        if(!mat)
        {
            matnotfound.insert(matid);
            continue;
        }
        TPZCompEl *cel = Mesh()->CreateCompEl(gel);
        TPZCompElPostProcBase *celpost = dynamic_cast<TPZCompElPostProcBase *>(cel);
        if(!celpost) DebugStop();
        TPZCompEl *celref = it->second;
        int nc = cel->NConnects();
        int ncref = celref->NConnects();
        TPZInterpolationSpace *celspace = dynamic_cast<TPZInterpolationSpace *>(cel);
        TPZInterpolationSpace *celrefspace = dynamic_cast<TPZInterpolationSpace *>(celref);
        int porder;
        if (!celrefspace) {
            TPZMultiphysicsElement *celrefmf = dynamic_cast<TPZMultiphysicsElement *>(celref);
            if (celrefmf){
                celrefspace = dynamic_cast<TPZInterpolationSpace *>(celrefmf->Element(0));
            } else {
                DebugStop();
            }
        }
        celpost->fReferredElement = celrefspace;

        if (celrefspace) {
            porder = celrefspace->GetPreferredOrder();
        } else {
            DebugStop();
        }

        celspace->SetPreferredOrder(porder);
        for (int ic=0; ic<nc; ic++) {
            cel->Connect(ic).SetOrder(porder,cel->ConnectIndex(ic));
            int nshape = celspace->NConnectShapeF(ic,porder);
            cel->Connect(ic).SetNShape(nshape);
        }

        TPZIntPoints &intrule = celspace->GetIntegrationRule();
        const TPZIntPoints &intruleref = celref->GetIntegrationRule();
        TPZIntPoints * cloned_rule = intruleref.Clone();
        cel->SetIntegrationRule(cloned_rule);

#ifdef PZDEBUG
        if (cel->GetIntegrationRule().NPoints() != intruleref.NPoints()) {
            DebugStop();
        }
#endif
        // this is why the mesh will be discontinuous!!
        gel->ResetReference();

	}

    // we changed the properties of the connects
    // now synchronize the connect properties with the block sizes
	int64_t nc= Mesh()->NConnects();
    for (int64_t ic=0; ic<nc; ic++) {
        TPZConnect &c = Mesh()->ConnectVec()[ic];
        int blsize = c.NShape()*c.NState();
        int64_t seqnum = c.SequenceNumber();
        Mesh()->Block().Set(seqnum, blsize);
    }
	Mesh()->InitializeBlock();
#ifdef PZ_LOG
  //  if(PPAnalysisLogger.isDebugEnabled())
    {
        std::stringstream sout;
        Mesh()->Print(sout);
       // LOGPZ_DEBUG(PPAnalysisLogger, sout.str())
    }
#endif

#ifdef PZDEBUG
	if(matnotfound.size())
	{
		std::cout << "Post-processing mesh was created without these materials: ";
		std::set<int>::iterator it;
		for(it = matnotfound.begin(); it!= matnotfound.end(); it++)
		{
			std::cout << *it << " ";
		}
		std::cout << std::endl;
        DebugStop();
	}
#endif
// 	TPZAdmChunkVector<TPZGeoEl *> &elvec = Mesh()->Reference()->ElementVec();
// 	long i, nelem = elvec.NElements();
// 	int neltocreate = 0;
// 	long index;
//     // build a data structure indicating which geometric elements will be post processed
//     fpMainMesh->LoadReferences();
//     std::map<TPZGeoEl *,TPZCompEl *> geltocreate;
//     for (i=0; i<nelem; i++) {
//         if (!elvec[i]) {
//             continue;
//         }
//         if (elvec[i]->Reference()) {
//             geltocreate[elvec[i]] = elvec[i]->Reference();
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
//         int printing = 0;
//         if (printing) {
//             gel->Print(cout);
//         }
//
//
//         TPZCompEl *cel =  Mesh()->CreateCompEl(gel);
//         TPZCompEl *celref = it->second;
//         int nc = cel->NConnects();
//         int ncref = celref->NConnects();
//         if (nc != ncref) {
//             DebugStop();
//         }
//         TPZInterpolationSpace *celspace = dynamic_cast<TPZInterpolationSpace *>(cel);
//         TPZInterpolationSpace *celrefspace = dynamic_cast<TPZInterpolationSpace *>(celref);
//         int porder = celrefspace->GetPreferredOrder();
// //        if (porder != 2) {
// //            std::cout << "I should stop porder = " << porder << std::endl;
// //        }
//         celspace->SetPreferredOrder(porder);
//         for (int ic=0; ic<nc; ic++) {
//             cel->Connect(ic).SetOrder(porder,cel->ConnectIndex(ic));
//             int nshape = celspace->NConnectShapeF(ic,porder);
//             cel->Connect(ic).SetNShape(nshape);
//         }
//         TPZIntPoints &intrule = celspace->GetIntegrationRule();
//         TPZVec<int> intorder(gel->Dimension(),0);
//         const TPZIntPoints &intruleref = celrefspace->GetIntegrationRule();
//         intruleref.GetOrder(intorder);
//         intrule.SetOrder(intorder);
// #ifdef DEBUG
//         if (intrule.NPoints() != intruleref.NPoints()) {
//             DebugStop();
//         }
// #endif
//         gel->ResetReference();
//
// 	}
//
//     // we changed the properties of the connects
//     // now synchronize the connect properties with the block sizes
// 	long nc= Mesh()->NConnects();
//     for (long ic=0; ic<nc; ic++) {
//         TPZConnect &c = Mesh()->ConnectVec()[ic];
//         int blsize = c.NShape()*c.NState();
//         long seqnum = c.SequenceNumber();
//         Mesh()->Block().Set(seqnum, blsize);
//     }
// 	Mesh()->InitializeBlock();
// #ifdef LOG4CXX
//     if(PPAnalysisLogger->isDebugEnabled())
//     {
//         std::stringstream sout;
//         Mesh()->Print(sout);
//         LOGPZ_DEBUG(PPAnalysisLogger, sout.str())
//     }
// #endif
// 	if(matnotfound.size())
// 	{
// 		std::cout << "Malha post proc was created without these materials ";
// 		std::set<int>::iterator it;
// 		for(it = matnotfound.begin(); it!= matnotfound.end(); it++)
// 		{
// 			std::cout << *it << " ";
// 		}
// 		std::cout << std::endl;
// 	}


}

void TPZPostProcAnalysis::Assemble()
{
   PZError << "Error at " << __PRETTY_FUNCTION__ << " TPZPostProcAnalysis::Assemble() should never be called\n";
}

void TPZPostProcAnalysis::Solve(){
   PZError << "Error at " << __PRETTY_FUNCTION__ << " TPZPostProcAnalysis::Solve() should never be called\n";
}

void TPZPostProcAnalysis::TransferSolution()
{
//     TPZLinearAnalysis::AssembleResidual();
//     fSolution = Rhs();
//     TPZLinearAnalysis::LoadSolution();
//
//     TPZCompMeshReferred *compref = dynamic_cast<TPZCompMeshReferred *>(Mesh());
//     if (!compref) {
//         DebugStop();
//     }
//     TPZCompMesh *solmesh = fpMainMesh;
//     long numelsol = solmesh->ElementSolution().Cols();
//     long nelem = compref->NElements();
//     compref->ElementSolution().Redim(nelem, numelsol);
//     if (numelsol)
//     {
//         for (long el=0; el<nelem; el++) {
//             TPZCompEl *cel = compref->ReferredEl(el);
//             if (!cel) {
//                 continue;
//             }
//             long index = cel->Index();
//             for (long isol=0; isol<numelsol; isol++) {
//                // compref->ElementSolution()(el,isol) = solmesh->ElementSolution()(index,isol);
//             }
//         }
//     }
//     TPZLinearAnalysis::AssembleResidual();
//     fSolution = Rhs();
//     TPZLinearAnalysis::LoadSolution();
//
//     TPZCompMeshReferred *compref = dynamic_cast<TPZCompMeshReferred *>(Mesh());
//     if (!compref) {
//         DebugStop();
//     }
//     TPZCompMesh *solmesh = fpMainMesh;
//
//     TPZFMatrix<STATE> &comprefElSol = compref->ElementSolution();
//     const TPZFMatrix<STATE> &solmeshElSol = solmesh->ElementSolution();
//
//
//     comprefElSol.Print("comprefElSol");
//     solmeshElSol.Print("solmeshElSol");
//     long numelsol = solmeshElSol.Cols();
//     long nelem = compref->NElements();
//     comprefElSol.Redim(nelem, numelsol);
//     if (numelsol)
//     {
//         for (long el=0; el<nelem; el++) {
//             TPZCompEl *cel = compref->ReferredEl(el);
//             if (!cel) {
//                 continue;
//             }
//             long index = cel->Index();
//             for (long isol=0; isol<numelsol; isol++) {
//
//                 comprefElSol(el,isol) = solmeshElSol.Get(index,isol);
//                 //compref->ElementSolution()(el,isol) = solmesh->ElementSolution()(index,isol);
//             }
//         }
//     }

        // this is where we compute the projection of the post processed variables
    TPZLinearAnalysis::AssembleResidual();
    fSolution = Rhs();
    TPZLinearAnalysis::LoadSolution();

    TPZCompMesh *compmeshPostProcess = (Mesh());
    if (!compmeshPostProcess) {
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
    compmeshPostProcess->ElementSolution().Redim(nelem, numelsol);
    if (numelsol)
    {
        for (int64_t el=0; el<nelem; el++) {
            TPZCompEl *celpost = compmeshPostProcess->Element(el);
            TPZGeoEl *gel = celpost->Reference();
            // we dont acount for condensed elements submeshes etc
            if(!gel) DebugStop();
            TPZCompEl *cel = gel->Reference();
            if (!cel) {
                DebugStop();
            }
            int64_t index = cel->Index();
            // we copy from the simulation mesh to the post processing mesh
            for (int64_t isol=0; isol<numelsol; isol++) {
                comprefElSol(el,isol) = solmeshElSol.Get(index,isol);
            }
        }
    }
}


void TPZPostProcAnalysis::SetAllCreateFunctionsPostProc(TPZCompMesh *cmesh)
{

    TPZManVector<TCreateFunction,10> functions(8);

    functions[EPoint] = &TPZPostProcAnalysis::CreatePointEl;
    functions[EOned] = TPZPostProcAnalysis::CreateLinearEl;
    functions[EQuadrilateral] = TPZPostProcAnalysis::CreateQuadEl;
    functions[ETriangle] = TPZPostProcAnalysis::CreateTriangleEl;
    functions[EPrisma] = TPZPostProcAnalysis::CreatePrismEl;
    functions[ETetraedro] = TPZPostProcAnalysis::CreateTetraEl;
    functions[EPiramide] = TPZPostProcAnalysis::CreatePyramEl;
    functions[ECube] = TPZPostProcAnalysis::CreateCubeEl;
    cmesh->ApproxSpace().SetCreateFunctions(functions);
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

TPZCompEl *TPZPostProcAnalysis::CreatePointEl(TPZGeoEl *gel,TPZCompMesh &mesh) {
	if(!gel->Reference() && gel->NumInterfaces() == 0)
		return new TPZCompElPostProc< TPZCompElH1<TPZShapePoint> >(mesh,gel);
	return NULL;
}
TPZCompEl *TPZPostProcAnalysis::CreateLinearEl(TPZGeoEl *gel,TPZCompMesh &mesh) {
	if(!gel->Reference() && gel->NumInterfaces() == 0)
		return new TPZCompElPostProc<TPZCompElH1<TPZShapeLinear> >(mesh,gel);
	return NULL;
}
TPZCompEl *TPZPostProcAnalysis::CreateQuadEl(TPZGeoEl *gel,TPZCompMesh &mesh) {
	if(!gel->Reference() && gel->NumInterfaces() == 0)
		return new TPZCompElPostProc<TPZCompElH1<TPZShapeQuad> >(mesh,gel);
	return NULL;
}
TPZCompEl *TPZPostProcAnalysis::CreateTriangleEl(TPZGeoEl *gel,TPZCompMesh &mesh) {
	if(!gel->Reference() && gel->NumInterfaces() == 0)
		return new TPZCompElPostProc<TPZCompElH1<TPZShapeTriang> >(mesh,gel);
	return NULL;
}
TPZCompEl *TPZPostProcAnalysis::CreateCubeEl(TPZGeoEl *gel,TPZCompMesh &mesh) {
	if(!gel->Reference() && gel->NumInterfaces() == 0)
		return new TPZCompElPostProc<TPZCompElH1<TPZShapeCube> >(mesh,gel);
	return NULL;
}
TPZCompEl *TPZPostProcAnalysis::CreatePrismEl(TPZGeoEl *gel,TPZCompMesh &mesh) {
	if(!gel->Reference() && gel->NumInterfaces() == 0)
		return new TPZCompElPostProc< TPZCompElH1<TPZShapePrism> >(mesh,gel);
	return NULL;
}
TPZCompEl *TPZPostProcAnalysis::CreatePyramEl(TPZGeoEl *gel,TPZCompMesh &mesh) {
	if(!gel->Reference() && gel->NumInterfaces() == 0)
		return new TPZCompElPostProc<TPZCompElH1<TPZShapePiram> >(mesh,gel);
	return NULL;
}
TPZCompEl *TPZPostProcAnalysis::CreateTetraEl(TPZGeoEl *gel,TPZCompMesh &mesh) {
	if(!gel->Reference() && gel->NumInterfaces() == 0)
		return new TPZCompElPostProc<TPZCompElH1<TPZShapeTetra> >(mesh,gel);
	return NULL;
}


TPZCompEl * TPZPostProcAnalysis::CreatePostProcDisc(TPZGeoEl *gel, TPZCompMesh &mesh)
{
	return new TPZCompElPostProc< TPZCompElDisc > (mesh,gel);
}

/** @brief Returns the unique identifier for reading/writing objects to streams */
int TPZPostProcAnalysis::ClassId() const{
    return Hash("TPZPostProcAnalysis") ^ TPZLinearAnalysis::ClassId() << 1;
}
/** @brief Save the element data to a stream */
void TPZPostProcAnalysis::Write(TPZStream &buf, int withclassid) const
{
    TPZLinearAnalysis::Write(buf, withclassid);
    TPZPersistenceManager::WritePointer(fpMainMesh, &buf);
}

/** @brief Read the element data from a stream */
void TPZPostProcAnalysis::Read(TPZStream &buf, void *context)
{
    TPZLinearAnalysis::Read(buf, context);
    fpMainMesh = dynamic_cast<TPZCompMesh*>(TPZPersistenceManager::GetInstance(&buf));
}
