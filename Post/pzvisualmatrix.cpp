/**
 * @file
 * @brief Contains the implementation of the VisualMatrix for DataExplorer and VTK packages. 
 */

#include "pzvisualmatrix.h"
#include "pzfmatrix.h"

using namespace std;

/** This function creates adequated file that allow to visualization of the value of a matrix passed as parameter. \n
 *  Depends on disponible visualization package 
 */
template<class TVar>
void VisualMatrix(TPZFMatrix<TVar> & matrix, const std::string &outfilename)
{
	int posdx = outfilename.rfind(".dx");
	int posvtk = outfilename.rfind(".vtk");
	int filelength = outfilename.size();
	if(filelength-posdx == 3) {
		VisualMatrixDX(matrix,outfilename);
	} else if(filelength-posvtk == 4) {
		VisualMatrixVTK(matrix,outfilename);
	} else {
		cout << "visualmatrix output was not created\n";
	}
}	

/** This function creates a Data Explorer file that allow to visualization of the value of a matrix passed as parameter */
template<class TVar>
void VisualMatrixDX(TPZFMatrix<TVar> & matrix, const std::string &outfilename)
{
	const int nelx = matrix.Cols();
	const int nely = matrix.Rows();
	const int neltotal = nelx * nely;
	int i,j;
	ofstream out(outfilename.c_str());
	out << "# Graphical Visualization of Matrix." << endl;
	out << "# Positions as the indexes of the matrix, beginning by column." << endl;
	out << "# The number of elements in x direction correspond to the number of the columns of the matrix." << endl;
	out << "# The number of elements in y direction correspond to the number of the rows of the matrix." << endl;
	
	out  << "object 1 class gridpositions counts " << nelx+1 << " " << nely +1 << endl;
	out << "origin 0. 0." << endl;
	out << "delta 1. 0." << endl;
	out << "delta 0. 1." << endl;
	out << "attribute \"dep\" string \"positions\"" << endl;
	out << endl;
	
	
	out << "object 2 class gridconnections counts " << nelx+1 << " " << nely +1 << endl;
	
 	out << "attribute \"element type\" string \"quads\"" << endl;
	out << "attribute \"ref\" string \"positions\"" << endl;
	
	out.precision(5);
	out  << "object 3 class array type float rank 0 items " << neltotal << " data follows" << endl;
	for (i = 0; i < nelx; i++) {
		for(j=0; j< nely ; j++) out << matrix(i,j) << endl;
	}
	out << "attribute \"dep\" string \"connections\" " << endl;
	out << endl;
	
	out << "object 4 class field" << endl;
	out << "component \"data\" value 3" << endl;
	out << "component \"positions\" value 1" << endl;
	out << "component \"connections\" value 2" << endl;
	out << "attribute \"name\" string \"Matrix\"" << endl;
	
	out << endl;
	out << "end" << endl;
	
	out.close();
	
	cout << "Data Explorer file " << outfilename << " was created with success!\n";
}

/** This function creates a Visualization Tool Kit (VTK) file that allow to visualization of the value of a matrix passed as parameter */
template<class TVar>
void VisualMatrixVTK(TPZFMatrix<TVar> & matrix, const std::string &outfilename)
{
	const int nelx = matrix.Cols();
	const int nely = matrix.Rows();
	const int neltotal = nelx * nely;
	int i,j;
	ofstream out(outfilename.c_str());
	out << "# vtk DataFile Version 3.0\n";
	out << "Generated by PZ\n";
	out << "ASCII\n";
	out << "DATASET RECTILINEAR_GRID\n";
	out << "DIMENSIONS " << (nelx+1) << " " <<  (nely+1) << " 1\n";
	out << "X_COORDINATES " << nelx+1 << " float\n";
	for (i=0; i<=nelx; i++) {
		out << i << " ";
	}
	out << std::endl;
	out << "Y_COORDINATES " << nely+1 << " float\n";
	for (j=0; j<=nely; j++) {
		out << j << " ";
	}
	out << std::endl;
	out << "Z_COORDINATES " << 1 << " float\n0.\n";
	out << "CELL_DATA " << nelx*nely << std::endl;
	out << "SCALARS mat_value float 1\n";
	out << "LOOKUP_TABLE default\n";
	const TVar *elem = 0;
    if (neltotal) {
        elem = &matrix(0,0);
    }
	for (i=0; i<neltotal; i++) {
		out << *(elem+i) << std::endl;
	}
}

template void VisualMatrix<float>(TPZFMatrix<float> & matrix, const std::string &outfilename);
template void VisualMatrix<double>(TPZFMatrix<double> & matrix, const std::string &outfilename);
template void VisualMatrix<long double>(TPZFMatrix<long double> & matrix, const std::string &outfilename);

template void VisualMatrix<std::complex<float> >(TPZFMatrix<std::complex<float> > & matrix, const std::string &outfilename);
template void VisualMatrix<std::complex<double> >(TPZFMatrix<std::complex<double> > & matrix, const std::string &outfilename);
template void VisualMatrix<std::complex<long double> >(TPZFMatrix<std::complex<long double> > & matrix, const std::string &outfilename);
