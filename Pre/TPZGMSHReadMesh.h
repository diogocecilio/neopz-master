/**
 * @file
 * @brief Contains the TPZGMSHReadMesh class which manages the manipulation of geometric meshes.
 */

#ifndef TPZGMSH_HPP
#define TPZGMSH_HPP

class TPZGeoMesh;
class TPZGeoEl;
class TPZGeoElSide;
#include "pzstack.h"
#include <fstream>
/**
 * @ingroup pre
 * @brief Manages the manipulation of geometric meshes. \ref pre "Getting Data"
 * generated by the GMSH in order to be used for the NeoPZ
 */
/** 
 * The archive generated for the GMSH must contain volume elements 
 * and contour elements, \n
 * The data of these elements must be: \n
 * - number of nodes \n
 * - co-ordinated of the nodes \n
 * - number of elements \n
 * - incidences of the elements
 */
class TPZGMSHReadMesh {
	
	/**
	 * @brief archive generated for the GMSH to be used (interpreted) inside of the NeoPZ
	 */
	std::ifstream fInGSMHGeoMesh;
	
	/**
	 * @brief defined geometric mesh in the NeoPZ
	 */
	TPZGeoMesh *fGeoMesh;
	
public:
	
	TPZGMSHReadMesh(TPZGeoMesh *gmesh);
	//TPZGMSHReadMesh(char *meshGMSH,TPZGeoMesh *gmesh);
	
	~TPZGMSHReadMesh(){};
	
	/*
	 * @brief Readings of meshes 2D, the contour and volume elements are returned in  
	 * elemlist and elembclist respectively 
	 *
	 * For a rectangular mesh the CC are: \n
	 * 1: for it wing \n
	 * 2: left lateral contour \n
	 * 3: inferior contour \n
	 * 4: right lateral contour \n
	 * 5: superior contour \n
	 * the plan is with axle X for right and axle Y for top
	 */
	void ReadMesh2D(const char *meshfile,TPZStack<TPZGeoEl *> &elemlist,TPZStack<TPZGeoElSide> &elembclist);
	void ReadMesh2D2(char *meshfile,TPZStack<TPZGeoEl *> &elemlist,TPZStack<TPZGeoElSide> &elembclist);
	/*
	 * @brief Readings of meshes 3D, the contour and volume elements are returned in  
	 * elemlist and elembclist respectively 
	 * 
	 * For one mesh of tetrahedrons the CC are \n
	 * 1: for it wing \n
	 * 2: left lateral contour \n
	 * 3: inferior contour \n
	 * 4: right lateral contour \n
	 * 5: superior contour \n
	 * 6: posterior contour (Z = 0) \n
	 * 7: previous contour (Z > 0) \n
	 * taking as base that the domain is with axle X for right, axle Y for top and axle -Z for the deep one
	 */
	void ReadMesh3D(char *meshfile,TPZStack<TPZGeoEl *> &elemlist,TPZStack<TPZGeoElSide> &elembclist);
	
	/**
	 * @brief Rearranges the nodal numeration given by the GMSH of sequential form
	 */
	void Resequence(TPZStack<long> &Indexes,const char *meshfile);
	
	/*
	 * @brief Prints in the exit defined for out the characteristics of the geometric mesh 
	 * created by the NeoPZ based on the mesh generated for the GMSH
	 */
	void PrintGeoMesh(std::ostream &out = std::cout);
	
};

#endif
