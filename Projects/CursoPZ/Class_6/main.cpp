
#include "pzvec.h"
#include "pzadmchunk.h"
#include "pzcmesh.h"
#include "pzvec_extras.h"
#include "pzdebug.h"
#include "pzcheckgeom.h"

#include "pzgeoel.h"
#include "pzgnode.h"
#include "pzgeoelside.h"
#include "pzgeoelbc.h"

#include "pzintel.h"
#include "pzcompel.h"

#include "pzmatrix.h"

#include "pzanalysis.h"
#include "pzfstrmatrix.h"
#include "pzskylstrmatrix.h"
#include "pzstepsolver.h"

#include "pzmaterial.h"
#include "pzbndcond.h"
#include "pzelasmat.h"
#include "pzplaca.h"
#include "pzmat2dlin.h"
#include "pzmathyperelastic.h"
#include "pzmattest3d.h"
#include "pzmatplaca2.h"

#include "pzfunction.h"
#include "TPZVTKGeoMesh.h"
#include "tpzgeoelrefpattern.h"
#include "pzlog.h"
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzcompel.h"
#include "pzgeoelside.h"
#include "TPZGeoLinear.h"
#include "pzgeopoint.h"
#include "tpzgeoblend.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.adaptivity"));
static LoggerPtr loggerconv(Logger::getLogger("pz.adaptivity.conv"));
static LoggerPtr loggerpoint(Logger::getLogger("pz.adaptivity.points"));
#endif

#include <time.h>
#include <stdio.h>
#include <fstream>

TPZGeoMesh *CreateGeoMesh();
TPZCompMesh *CreateMesh(TPZGeoMesh *gmesh);
TPZGeoMesh *CreateGeoMesh2();
TPZCompMesh *CreateMesh2(TPZGeoMesh *gmesh);

void studysol( TPZCompMesh * cmesh, TPZAnalysis  an);
REAL ProcessingError(TPZAnalysis &analysis, TPZVec<REAL> &ervec, TPZVec<REAL> &ervecbyel, TPZVec<REAL> &gradervecbyel, REAL &MinErrorByElement, REAL &MinGradErrorByElement);
// bi-dimensional problem for elasticity
int main() {
	
#ifdef LOG4CXX
	InitializePZLOG();
#endif
	
	// Creating geometric mesh
	TPZGeoMesh *gmesh = CreateGeoMesh2();

	// Creating computational mesh (approximation space and materials)
	int p = 1;
    TPZCompEl::SetgOrder(p);
	//gmesh->Print();
    TPZCompMesh *cmesh = CreateMesh2(gmesh);
	// Solving linear equations
	// Initial steps
	TPZAnalysis an (cmesh);
	TPZSkylineStructMatrix strskyl(cmesh);
	an.SetStructuralMatrix(strskyl);
	// Solver (is your choose) 
	TPZStepSolver<STATE> *direct = new TPZStepSolver<STATE>;
	direct->SetDirect(ECholesky);
	an.SetSolver(*direct);
	delete direct;
	direct = 0;

	an.Run();

	//studysol(cmesh, an);
	//an.Rhs().Print("RHS");
	//an.Mesh()->;
	an.Solution().Print("SOL");
	// Post processing
	TPZManVector<std::string> scalarnames(3), vecnames(1);
	scalarnames[0] = "SigmaX";
	scalarnames[1] = "SigmaY";
	scalarnames[2] = "Pressure";
	vecnames[0] = "displacement";
	//vecnames[1] = "";
	an.DefineGraphMesh(2,scalarnames,vecnames,"ElasticitySolutions.vtk");

	an.PostProcess(0);
	std::cout << "HELLO " << std::endl;
	system("Pause");
	return 0;
}

void studysol(TPZCompMesh * cmesh, TPZAnalysis  an)
{
	TPZVec<REAL> psi(3, 0.);
	TPZVec<STATE> Laplacian(1);
	TPZVec<REAL> center(3, 0.);
	TPZVec<long> subels;
	long MaxHLevel =  2;
	
	long nels = cmesh->ElementVec().NElements();

	REAL MinErrorByElement, MinGradErrorByElement;
	TPZVec<REAL> gradervecbyel,ervec, ervecbyel,ErrorVec(100, 0.0);
	ervecbyel.Resize(0);
	gradervecbyel.Resize(0);
	REAL MaxErrorByElement = ProcessingError(an, ervec, ervecbyel, gradervecbyel, MinErrorByElement, MinGradErrorByElement);

	REAL factorGrad = 1./3.;
	//REAL MaxGrad = factorGrad*gradervecbyel[nels] + (1. - factorGrad)*MinGrad;
	REAL MaxGrad = factorGrad*gradervecbyel[nels];
	for (long el = 0; el < nels; el++)
	{
		TPZCompEl* cel = cmesh->ElementVec()[el];
		cel->Reference()->CenterPoint(cel->Reference()->NSides() - 1, psi);
		cel->Reference()->X(psi, center);

		long index = cel->Index();
		long level = cel->Reference()->Level();

		// Applying hp refinement depends on high gradient and high laplacian value, and depends on computed error by element
		if (gradervecbyel[el] > MaxGrad && level < MaxHLevel) {
			cel->Divide(index, subels);
			cel = NULL;
			level++;
		}
	}
}

/**
* Get Global L2 Error for solution and the L2 error for each element.
* Return the maxime L2 error by elements. Also return in MinErrorByElement argument the minime L2 error for all elements of the mesh.
*/

REAL ProcessingError(TPZAnalysis &analysis, TPZVec<REAL> &ervec, TPZVec<REAL> &ervecbyel, TPZVec<REAL> &gradervecbyel, REAL &MinErrorByElement, REAL &MinGradErrorByElement) {
	long neq = analysis.Mesh()->NEquations();
	int ModelDimension = 2;

	TPZVec<REAL> ux(neq);
	TPZVec<REAL> sigx(neq);
	TPZManVector<REAL, 10> totalerror(10, 0.);
	analysis.Mesh()->LoadSolution(analysis.Solution());

	TPZAdmChunkVector<TPZCompEl *> elvec = analysis.Mesh()->ElementVec();
	TPZManVector<REAL, 10> errors(10);
	errors.Fill(0.0);
	long i, nel = elvec.NElements();
	ervecbyel.Resize(nel, 0.0);
	// The last position will be store the maxime value of the gradient errors
	gradervecbyel.Resize(nel + 1, 0.0);
	REAL maxError = 0.0;
	MinErrorByElement = 1000.0;
	MinGradErrorByElement = 10000.0;

	/** Computing error for all elements with same dimension of the model */
	for (i = 0L; i<nel; i++) {
		TPZCompEl *el = (TPZCompEl *)elvec[i];
		if (!el || el->Dimension() != ModelDimension) continue;
		if (el) {
			errors.Fill(0.0);
			el->EvaluateError(analysis.fExact, errors, 0);
			int nerrors = errors.NElements();
			totalerror.resize(nerrors);
			for (int ier = 0; ier < nerrors; ier++)
			{
				totalerror[ier] += errors[ier] * errors[ier];
			}
			// L2 error for each element
			if (i % 100 == 0)
			{
				std::cout << "Computed " << i << " elements from " << nel << " total error " << sqrt(totalerror[1]) << std::endl;
			}
			ervecbyel[i] = sqrt(errors[1] * errors[1]);
			gradervecbyel[i] = sqrt(errors[2] * errors[2]);
			if (gradervecbyel[i] > gradervecbyel[nel])
				gradervecbyel[nel] = gradervecbyel[i];
			if (gradervecbyel[i] < MinGradErrorByElement)
				MinGradErrorByElement = gradervecbyel[i];
			// The computed error by current element is compared with max and min values to return
			if (ervecbyel[i] > maxError)
				maxError = ervecbyel[i];
			else if (ervecbyel[i] < MinErrorByElement)
				MinErrorByElement = ervecbyel[i];
		}
	}

	int nerrors = errors.NElements();
	ervec.Resize(nerrors);
	ervec.Fill(-1.0);

	// Returns the square of the calculated errors.
	for (i = 0; i<nerrors; i++)
		ervec[i] = sqrt(totalerror[i]);
	return maxError;
}



//*******Shell to deforming************
TPZGeoMesh *CreateGeoMesh2() {


	TPZGeoMesh *gmesh = new TPZGeoMesh();
	gmesh->NodeVec().Resize(4);
	TPZVec<REAL> coord(2);
	coord[0] = 0.;
	coord[1] = 0.;
	gmesh->NodeVec()[0] = TPZGeoNode(0, coord, *gmesh);
	coord[0] = 8.;
	coord[1] = 0.;
	gmesh->NodeVec()[1] = TPZGeoNode(1, coord, *gmesh);
	coord[0] = 8.;
	coord[1] = 4.;
	gmesh->NodeVec()[2] = TPZGeoNode(2, coord, *gmesh);
	coord[0] = 0.;
	coord[1] = 4.;
	gmesh->NodeVec()[3] = TPZGeoNode(3, coord, *gmesh);

	TPZVec <long> TopoQuad(4);
	TopoQuad[0] = 0;
	TopoQuad[1] = 1;
	TopoQuad[2] = 2;
	TopoQuad[3] = 3;
	long index;
	//gmesh->CreateGeoElement(EQuadrilateral, TopoQuad, 1, index);
	//new TPZGeoElRefPattern< pzgeom::TPZGeoQuad>(elid, TopolQuad, matid, *gmesh);
	new TPZGeoElRefPattern< pzgeom::TPZGeoQuad>(0, TopoQuad, 1, *gmesh);

	//gmesh->ElementVec()[index];
	
	TPZVec <long> TopoLine(2);
	TopoLine[0] = TopoQuad[0];
	TopoLine[1] = TopoQuad[1];
	new TPZGeoElRefPattern< pzgeom::TPZGeoLinear>(1, TopoLine, -2, *gmesh);
	//gmesh->CreateGeoElement(EOned, TopoLine, -1, index);
	//gmesh->ElementVec()[index];

	TopoLine[0] = TopoQuad[1];
	TopoLine[1] = TopoQuad[2];
	//gmesh->CreateGeoElement(EOned, TopoLine, -2, index);
	new TPZGeoElRefPattern< pzgeom::TPZGeoLinear>(2, TopoLine, -3, *gmesh);
	//gmesh->ElementVec()[index];

	TopoLine[0] = TopoQuad[2];
	TopoLine[1] = TopoQuad[3];
	//gmesh->CreateGeoElement(EOned, TopoLine,-3, index);
	new TPZGeoElRefPattern< pzgeom::TPZGeoLinear>(3, TopoLine, -4, *gmesh);
	//gmesh->ElementVec()[index];

	TopoLine[0] = TopoQuad[3];
	TopoLine[1] = TopoQuad[0];
	new TPZGeoElRefPattern< pzgeom::TPZGeoLinear>(4, TopoLine, -5, *gmesh);
	//gmesh->CreateGeoElement(EOned, TopoLine, -4, index);
	//gmesh->ElementVec()[index];

	 gmesh->BuildConnectivity();

	for (int d = 0; d<0; d++) {
		int nel = gmesh->NElements();
		for (int iel = 0; iel<nel; iel++) {
			TPZManVector<TPZGeoEl *> subels;
			TPZGeoEl *gel = gmesh->ElementVec()[iel];
			gel->Divide(subels);
		}
	}

	return gmesh;
}

/** FORCING FUNCTION */
void ForcingFunction2(const TPZVec<REAL> &x, TPZVec<STATE> &sol, TPZFMatrix<STATE> &dsol)
{
	std::cout << "x = " << x[0] << std::endl;
	std::cout << "y = " << x[1] << std::endl;
	sol[1] = 100.*sin(M_PI*x[0]/16.);
	//sol[1] = 10000.;
	dsol.Zero();
}



//*************************************
//*******L Shape Quadrilateral*********
//*************************************
TPZCompMesh *CreateMesh2(TPZGeoMesh *gmesh) {

	TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
	// cmesh->SetDefaultOrder(TPZCompEl::GetgOrder());
	//cmesh->SetPord

	// Creating elasticity material
	TPZMaterial * mat = new TPZElasticityMaterial(1, 30000000, 0.3, 0., 0.,1);

	// Creating four boundary condition
	TPZFMatrix<STATE> val1(2, 2, 0.), val2(2, 1, 0.);
	TPZMaterial *bcBottom, *bcRight,*bcTop;

	
	val2(1, 0) = 1.;
	bcBottom = mat->CreateBC(mat, -2, 3, val1, val2);

	val2(0, 0) = 1.;
	val2(1, 0) = 0.;
	bcRight = mat->CreateBC(mat, -3, 3, val1, val2);

	val2(1, 0) = 0.;
	val2(0, 0) = 0.;
	bcTop = mat->CreateBC(mat, -4, 1, val1, val2);

	bcTop->SetForcingFunction(new TPZDummyFunction<STATE>(ForcingFunction2));

	cmesh->InsertMaterialObject(mat);
	cmesh->InsertMaterialObject(bcBottom);
	cmesh->InsertMaterialObject(bcRight);
	cmesh->InsertMaterialObject(bcTop);

	cmesh->SetAllCreateFunctionsContinuous();

	cmesh->AutoBuild();
	cmesh->AdjustBoundaryElements();
	cmesh->CleanUpUnconnectedNodes();

#ifdef LOG4CXX
	if (logger->isDebugEnabled())
	{
		std::stringstream sout;
		cmesh->Print(sout);
		LOGPZ_DEBUG(logger, sout.str())
	}
#endif
	return cmesh;
}



//*******Shell to deforming************
TPZGeoMesh *CreateGeoMesh() {

    REAL co[8][2] = {{0.,0.},{0.,-1.},{1.,-1.},{1.,0.},{1.,1.},{0.,1.},{-1.,1.},{-1.,0.}};
    long indices[3][4] = {{0,1,2,3},{0,3,4,5},{0,5,6,7}};
    TPZGeoEl *elvec[3];
    TPZGeoMesh *gmesh = new TPZGeoMesh();
    long nnode = 8;
    long nod;
    for(nod=0; nod<nnode; nod++) {
        long nodind = gmesh->NodeVec().AllocateNewElement();
        TPZVec<REAL> coord(2);
        coord[0] = co[nod][0];
        coord[1] = co[nod][1];
        gmesh->NodeVec()[nodind] = TPZGeoNode(nod,coord,*gmesh);
    }
    
    long el;
    long nelem = 3;
    for(el=0; el<nelem; el++) {
        TPZVec<long> nodind(4);
        for(nod=0; nod<4; nod++) nodind[nod]=indices[el][nod];
        //    elvec[el] = new TPZGeoElQ2d(el,nodind,1);
        long index;
        elvec[el] = gmesh->CreateGeoElement(EQuadrilateral,nodind,1,index);
    }
    
	//TPZVec <long> TopoLine(2);


	//TopoLine[0] = 0;
	//TopoLine[1] = 1;
	//long index = -1;
	//gmesh->CreateGeoElement(EOned, TopoLine, 1, index);

	//TopoLine[0] = 1;
	//TopoLine[1] = 2;
	//index = -2;
	//gmesh->CreateGeoElement(EOned, TopoLine, 1, index);

	//TopoLine[0] = 2;
	//TopoLine[1] = 3;
	//index = -3;
	//gmesh->CreateGeoElement(EOned, TopoLine, 1, index);

	//TopoLine[0] = 3;
	//TopoLine[1] = 4;
	//index = -4;
	//gmesh->CreateGeoElement(EOned, TopoLine, 1, index);

	//TopoLine[0] = 4;
	//TopoLine[1] = 5;
	//index = -5;
	//gmesh->CreateGeoElement(EOned, TopoLine, 1, index);

	//TopoLine[0] = 5;
	//TopoLine[1] = 6;
	//index = -6;
	//gmesh->CreateGeoElement(EOned, TopoLine, 1, index);

	//TopoLine[0] = 6;
	//TopoLine[1] = 7;
	//index = -7;
	//gmesh->CreateGeoElement(EOned, TopoLine, 1,index);


	//TopoLine[0] = 7;
	//TopoLine[1] = 0;
	//index = -8;
	//gmesh->CreateGeoElement(EOned, TopoLine, 1, index);


    gmesh->BuildConnectivity();
	gmesh->Print();
    TPZVec<TPZGeoEl *> sub; 
	/** @brief Creates a geometric element along side of el.  */
	/** The new geometric element is inserted in mesh and a pointer to it is stored here. */
	//TPZGeoElBC(TPZGeoEl *el, int side, int matid);

    // bc -1 -> Dirichlet
    TPZGeoElBC gbc1(elvec[0],4,-1);

    // bc -2 -> Neumann at the bottom y==-1
    TPZGeoElBC gbc2(elvec[0],5,-2);


    // bc -3 -> Neumann at the right x==1
    TPZGeoElBC gbc3(elvec[0],6,-3);


    // bc -3 -> Neumann at the right x==1
    TPZGeoElBC gbc4(elvec[1],5,-3);
    
    // bc -4 -> Neumann at the top y==1
    TPZGeoElBC gbc5(elvec[1],6,-4);
    
    // bc -4 -> Neumann at the top y==1
    TPZGeoElBC gbc6(elvec[2],5,-4);
    
    // bc -5 -> Neumann at the left x==-1
    TPZGeoElBC gbc7(elvec[2],6,-5);
    
    // bc -6 -> Homogeneous Neumann
    TPZGeoElBC gbc8(elvec[2],7,-6);


    for(int d=0;d<2;d++) {
    int nel = gmesh->NElements();
    for (int iel=0; iel<nel; iel++) {
        TPZManVector<TPZGeoEl *> subels;
        TPZGeoEl *gel = gmesh->ElementVec()[iel];
        gel->Divide(subels);
    }
	}
	//gmesh->Print();


	//TPZVec <long> TopoLine(2);
	//int bottonmat = 2;
	//elid = 4;
	//TopoLine[0] = 1;
	//TopoLine[1] = 2;
	//new TPZGeoElRefPattern< pzgeom::TPZGeoLinear>(elid, TopoLine, bottonmat, *gmesh);
	
	gmesh->Print();

	return gmesh;
}



/** FORCING FUNCTION */
void ForcingFunction(const TPZVec<REAL> &x, TPZVec<STATE> &sol,  TPZFMatrix<STATE> &dsol)
{
	sol[1] =100.*sin(M_PI*x[0]);
	dsol.Zero();
}

//*************************************
//*******L Shape Quadrilateral*********
//*************************************
TPZCompMesh *CreateMesh(TPZGeoMesh *gmesh) {
    
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(TPZCompEl::GetgOrder());

    // Creating elasticity material
    TPZMaterial * mat = new TPZElasticityMaterial(1,2000000000.,0.3,0.,0.);

	// Creating four boundary condition
    TPZFMatrix<STATE> val1(2,2,0.),val2(2,1,0.);
	TPZMaterial *bcBottom, *bcRight,*bcdisplacement;
	//val1(1,1) = 1000000.;

	/** @} */

	/** @brief Creates an object TPZBndCond derived of TPZMaterial*/
	//virtual TPZBndCond *CreateBC(TPZMaterial *reference, int id, int typ, TPZFMatrix<STATE> &val1,
		//TPZFMatrix<STATE> &val2);

	//val2(1,0) = 100.;
   // bcBottom = mat->CreateBC(mat,-2,1,val1,val2);

	val2(1,0) = -10000.;
    bcRight = mat->CreateBC(mat,-2,1,val1,val2);
	//bcRight->SetfBCForcingFunction(new TPZDummyFunction<STATE>(ForcingFunction));
	//bcRight->SetForcingFunction(new TPZDummyFunction<STATE>(ForcingFunction));
	//bc1->SetForcingFunction(new TPZDummyFunction<STATE>(ExactSolLaplaceBC));


	val2(1, 0) = 0.;
	val2(0, 0) = 0.;
	val1(0, 0) = 10e12;
	val1(1, 1) = 10e12;
	bcdisplacement = mat->CreateBC(mat, -5, 0, val1, val2);

    cmesh->InsertMaterialObject(mat);
	// Inserting boundary conditions into computational mesh
	//cmesh->InsertMaterialObject(bcBottom);
	cmesh->InsertMaterialObject(bcRight);
	cmesh->InsertMaterialObject(bcdisplacement);
	cmesh->SetAllCreateFunctionsContinuous();
    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    
#ifdef LOG4CXX
    if (logger->isDebugEnabled())
    {
        std::stringstream sout;
        cmesh->Print(sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    return cmesh;
}

