#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>

#ifndef REFPATTERNDIR
#define REFPATTERNDIR "/home/santos/PZlib_teste/install/include/refpatterns"
#endif

#ifndef REALdouble
#define REALdouble
#endif

#ifndef STATEdouble
#define STATEdouble
#endif

#include <TPZRefPatternDataBase.h>
#include <TPZRefPattern.h>
#include <TPZRefPatternTools.h>

#include <pzreal.h>
#include <pzsave.h>
#include <pzgmesh.h>
#include <pzvec.h>
#include <pzeltype.h>

#include <tpzchangeel.h>
#include <TPZGeoElement.h>
#include <pzreftriangle.h>
#include <pzgeotriangle.h>
#include <tpzgeoelrefpattern.h>
#include <TPZVTKGeoMesh.h>

int main(int argc, char *argv[]) {

	//TPZErrorHandler::init(argc, argv);

	TPZGeoMesh *gmesh = new TPZGeoMesh();
	int nnodes = 4;
	const int mat = 1;
	const int reftype = 1; //1 is refpatern
	const int hmax = 1;
	long index;
	TPZAutoPointer<TPZRefPattern> refp = NULL;
	TPZManVector<REAL, 3> coord(3, 0.);
	TPZManVector<long, 3> tria(3, 0);
	TPZVec<TPZGeoEl *> sons;
	gRefDBase.InitializeUniformRefPattern(ETriangle);

	gmesh->NodeVec().Resize(nnodes);
	coord[0] = 0.; coord[1] = 0.; //nó 0
	gmesh->NodeVec()[0].SetCoord(coord); gmesh->NodeVec()[0].SetNodeId(0);
	coord[0] = 10.; coord[1] = 0.; //nó 1
	gmesh->NodeVec()[1].SetCoord(coord); gmesh->NodeVec()[1].SetNodeId(1);
	coord[0] = 10.; coord[1] = 10.; //nó 2
	gmesh->NodeVec()[2].SetCoord(coord); gmesh->NodeVec()[2].SetNodeId(2);
	coord[0] = 0.; coord[1] = 10.; //nó 3
	gmesh->NodeVec()[3].SetCoord(coord); gmesh->NodeVec()[3].SetNodeId(3);

	tria[0] = 3; tria[1] = 1; tria[2] = 2; //element 0
	gmesh->CreateGeoElement(ETriangle, tria, mat, index, reftype);
	gmesh->ElementVec()[index]->SetId(index);
	tria[0] = 0; tria[1] = 1; tria[2] = 3; //element 1
	gmesh->CreateGeoElement(ETriangle, tria, mat, index, reftype);
	gmesh->ElementVec()[index]->SetId(index);

	gmesh->BuildConnectivity();

	std::ofstream file0("C:/neopz-cmake/Projects/CursoPZ/ExemploSimples2D/mesh0.txt");     gmesh->Print(file0);
	std::ofstream filevtk0("C:/neopz-cmake/Projects/CursoPZ/ExemploSimples2D/mesh0.vtk");  TPZVTKGeoMesh::PrintGMeshVTK(gmesh, filevtk0);

	/* Refinement process */
	for (int l = 0; l<hmax; l++) {
		const long nelem = gmesh->NElements();
		for (long i = 0; i<nelem; i++) {
			if (gmesh->Element(i)->HasSubElement()) continue;
			gmesh->Element(i)->Divide(sons); sons.clear();
		}
	}
	if (reftype == 0) {
		refp = TPZRefPatternTools::PerfectMatchRefPattern(gmesh->Element(1));
		if (refp) {
			gmesh->Element(1)->SetRefPattern(refp);
			gmesh->Element(1)->Divide(sons); sons.clear();
		}
		else {
			DebugStop();
		}
	}
	gmesh->BuildConnectivity();

	std::ofstream file1("C:/neopz-cmake/Projects/CursoPZ/ExemploSimples2D/mesh1.txt");     gmesh->Print(file1);
	std::ofstream filevtk1("C:/neopz-cmake/Projects/CursoPZ/ExemploSimples2D/mesh1.vtk");  TPZVTKGeoMesh::PrintGMeshVTK(gmesh, filevtk1);

	/* Delete process */
	index = 1;//element to be deleted
	gmesh->Element(index)->GetHigherSubElements(sons);
	gmesh->Element(index)->ResetSubElements();
	for (int i = 0; i<sons.size(); i++) {
		gmesh->DeleteElement(sons[i], sons[i]->Index());
	}
	sons.clear();
	gmesh->BuildConnectivity();

	std::ofstream file2("C:/neopz-cmake/Projects/CursoPZ/ExemploSimples2D/mesh2.txt");     gmesh->Print(file2);
	std::ofstream filevtk2("C:/neopz-cmake/Projects/CursoPZ/ExemploSimples2D/mesh2.vtk");  TPZVTKGeoMesh::PrintGMeshVTK(gmesh, filevtk2);

	/* Refine again*/
	if (1)
	{ //Divide again, uniform
		gmesh->Element(index)->Divide(sons); sons.clear();
	}
	else
	{ //Divide again, not uniform
		if (gmesh->Element(index)->HasSubElement()) DebugStop();
		refp = TPZRefPatternTools::PerfectMatchRefPattern(gmesh->Element(index));
		if (refp) {
			gmesh->Element(index)->SetRefPattern(refp);
			gmesh->Element(index)->Divide(sons); sons.clear();
		}
		else {
			DebugStop();
		}
	}
	gmesh->BuildConnectivity();

	std::ofstream file3("C:/neopz-cmake/Projects/CursoPZ/ExemploSimples2D/mesh3.txt");     gmesh->Print(file3);
	std::ofstream filevtk3("C:/neopz-cmake/Projects/CursoPZ/ExemploSimples2D/mesh3.vtk");  TPZVTKGeoMesh::PrintGMeshVTK(gmesh, filevtk3);
	std::cout << "HELLO " << std::endl;
	system("Pause");
	return 0;
}
