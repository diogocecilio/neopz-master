//
// Created by Diogo Cecílio on 10/12/21.
//

#pragma once

#include <cmath>
#include <set>

#include <iostream>
#include <fstream>
#include <string>
#include "pzgmesh.h"
#include "pzstack.h"
#include "TPZVTKGeoMesh.h"
//#include "TPZDarcyMaterial.h"

#include <pzgeoel.h>
#include "pzgeoelbc.h"
#include "pzfmatrix.h"
#include "pzbstrmatrix.h"
#include <TPZGeoElement.h>
#include "TPZVTKGeoMesh.h"
#include "pzbuildmultiphysicsmesh.h"
#include "TPZInterfaceEl.h"
#include "TPZMultiphysicsInterfaceEl.h"
#include "pzfstrmatrix.h"
#include "pzskylstrmatrix.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "pzstepsolver.h"
#include "TPZGeoLinear.h"
#include "tpzgeoelrefpattern.h"
#include "TPZParFrontStructMatrix.h"
#include <math.h>
#include <cmath>
#include <iostream>
#include <chrono>

using namespace std;

class readgidmesh
{
public:
    /**
     * @brief Default constructor
     */
    readgidmesh();

    /**
     * @brief Class constructor
     * @param [in] mesh
     */
    readgidmesh ( string file );
    ~readgidmesh();




    template <class T>
    std::vector<T> str_vec( std::vector<std::string> &vs );


    void ReadMesh ( );

    void ReadMesh2 ( std::vector<std::vector<int>>& topol, std::vector<std::vector<double>> &coords,std::string ffile);
	

    
    void  FindIds ( TPZVec<double> constcoorddata,TPZVec<int> constcoord, std::vector<int>& ids )
    {
		REAL tol = 1.e-12;
		//constcoorddata vector containig info about face to search id. It must contain any coodinate locate in the in the face
        TPZFMatrix<REAL> elcoords;
		std::vector<double> elcoodsvec;
        int nels = fallcoords.size();
        GetElCoords (  0, elcoords );
        int nnodes = elcoords.Rows();
        int sum=0;
        //constcoord.size() = 1 face
        //constcoord.size() = 2 linha
        //constcoord.size() = 3 pontos

        std::vector<int> dirs;
        for ( int iconst=0; iconst<constcoord.size(); iconst++ ) {
            sum+=constcoord[iconst];
            if ( constcoord[iconst]==1 ) {
                dirs.push_back ( iconst );
            }

        }
        for ( int iel = 0; iel < nels; iel++ ) {
            GetElCoords (  iel, elcoords );

            for ( int inode = 0; inode < nnodes; inode++ ) {

                if ( sum==1 ) {
                    if ( fabs ( elcoords(inode,dirs[0]) - constcoorddata[dirs[0]] ) <tol ) {
                        ids.push_back ( fmeshtopology(iel,inode) );
                    }
                } else if ( sum==2 ) {
                    if ( fabs ( elcoords(inode,dirs[0]) - constcoorddata[dirs[0]] ) <tol && abs ( elcoords(inode,dirs[1]) - constcoorddata[dirs[1]] ) <tol ) {
                        ids.push_back (  fmeshtopology(iel,inode)  );
                    }
                } else if ( sum==3 ) {
                    if ( fabs ( elcoords(inode,dirs[0]) - constcoorddata[dirs[0]] ) <tol && abs ( elcoords(inode,dirs[1]) - constcoorddata[dirs[1]] ) <tol && abs ( elcoords(inode,dirs[2]) - constcoorddata[dirs[2]] ) <tol ) {
                        ids.push_back (  fmeshtopology(iel,inode)  );
                    }
                }


            }


        }

        sort ( ids.begin(), ids.end() );
        ids.erase ( unique ( ids.begin(), ids.end() ), ids.end() );
    }
    void  GetElCoords ( int el, TPZFMatrix<double>  & elcoords )
    {
        elcoords.Resize ( fallcoords[el].size(), 3 );
		
        for ( int j = 0; j < fallcoords[el].size(); j++ ) {
            double x = fallcoords[el][j][0];
            double y = fallcoords[el][j][1];
            double z = fallcoords[el][j][2];
            elcoords(j,0) = x;
            elcoords(j,1) = y;
            elcoords(j,2) = z;
        }
    }
    
	void  GetElCoords ( int el, vector<double>  & elcoords )
    {
        elcoords.resize ( 3 );
		
        for ( int j = 0; j < fallcoords[el].size(); j++ ) {
            double x = fallcoords[el][j][0];
            double y = fallcoords[el][j][1];
            double z = fallcoords[el][j][2];
            elcoords[0] = x;
            elcoords[1] = y;
            elcoords[2] = z;
        }
    }
	
	inline TPZFMatrix<int> GetTopology()
	{
		return fmeshtopology;
	}
	inline TPZFMatrix<double> GetCoords()
	{
		return fmeshcoords;
	}
	std::vector<std::vector< std::vector<double > > > GetAllCoords()
	{
		return fallcoords;
	}

    std::vector<std::vector<int>>   LineTopology ( std::vector<int> ids )
    {
        int  k = 0;


        std::vector<std::vector<int>> vg;
        for ( int j = 0; j < ids.size() / (fOrder+1); j++ ) {
            std::vector<int> v;
            for ( int i = 0; i < fOrder + 1; i++ ) {
                v.push_back ( ids[i + k] );
            }
            vg.push_back ( v );
            k += fOrder;
        }
        return vg;
    }

void  ToMatInt ( std::vector<std::vector<int>> in, TPZFMatrix<int> & out )
{
    int  rows = in.size();
    int cols = in[0].size();
    out.Resize ( rows, cols );
    for ( int i = 0; i < rows; i++ ) for ( int j = 0; j < cols; j++ ) out(i,j)= in[i][j];
}
	
	
void FindIdsInPath ( const TPZFMatrix<REAL>& path, std::vector<int>& idpath,REAL delta )
{
    TPZFMatrix<REAL> elcoords;
    int nels = fallcoords.size();
    GetElCoords (  0, elcoords );
    int nnodes = elcoords.Rows();
    for ( int iel = 0; iel < nels; iel++ ) {
        GetElCoords ( iel, elcoords );
        for ( int inode = 0; inode < nnodes; inode++ ) {
            REAL x = elcoords(inode,0);
            REAL y = elcoords(inode,1);

            for ( int ipath = 0; ipath < path.Rows(); ipath++ ) {
                REAL copathx = path.Get(ipath,0);
                REAL copathy = path.Get(ipath,1);

                if ( fabs ( x - copathx ) < delta && fabs ( y - copathy ) <delta) {
                    idpath.push_back ( fmeshtopology(iel,inode) );
                    ipath = path.Rows();
                }
            }
        }

    }

    sort ( idpath.begin(), idpath.end() );
    idpath.erase ( unique ( idpath.begin(), idpath.end() ), idpath.end() );

}

void Line ( TPZManVector<REAL,2> a, TPZManVector<REAL,2> b, int ndivs, TPZFMatrix<REAL>& path )
{
    double x0 = a[0];
    double xf = b[0];

    double y0 = a[1];
    double yf = b[1];

    double dx = ( xf - x0 ) / ndivs;
    double dy = ( yf - y0 ) / ndivs;

    path.Resize( ndivs, 2 );

    for ( int idiv = 0; idiv < ndivs; idiv++ ) {
        path(idiv,0) = x0 + idiv * dx;
        path(idiv,1) = y0 + idiv * dy;
    }

}


void FindElements(TPZVec<double> constcoorddata,TPZVec<int> constcoord, std::vector<int>& ids)
	{
		REAL tol = 1.e-12;
		//constcoorddata vector containig info about face to search id. It must contain any coodinate locate in the in the face
        TPZFMatrix<REAL> elcoords;
		std::vector<double> elcoodsvec;
        int nels = fallcoords.size();
        GetElCoords (  0, elcoords );

        int sum=0;

        std::vector<int> dirs;
        for ( int iconst=0; iconst<constcoord.size(); iconst++ ) {
            sum+=constcoord[iconst];
            if ( constcoord[iconst]==1 ) {
                dirs.push_back ( iconst );
            }

        }
        for ( int iel = 0; iel < nels; iel++ ) {
            GetElCoords (  iel, elcoords );
			int nnodes = elcoords.Rows();
			if(nnodes==4)continue;//warning this only works for 3D tetrahedrom and triangle
            for ( int inode = 0; inode < nnodes; inode++ ) {
                if ( sum!=1 )DebugStop();
				if ( fabs ( elcoords(0,dirs[0]) - constcoorddata[dirs[0]] ) <tol  && fabs ( elcoords(1,dirs[0]) - constcoorddata[dirs[0]] ) <tol && fabs ( elcoords(2,dirs[0]) - constcoorddata[dirs[0]] ) <tol) {
					ids.push_back ( iel );
				}
            }


        }

        sort ( ids.begin(), ids.end() );
        ids.erase ( unique ( ids.begin(), ids.end() ), ids.end() );
	}



template <class T>
void PrintVecVec(std::vector<std::vector<T>>& data)
{
    int rows= data.size();
    for(int i=0;i<rows;i++)
    {
        int cols=data[i].size();
        for(int j=0;j<cols;j++)
        {
            cout << data[i][j] << " ";
        }
        cout << endl;
    }

}
private:

string ffile;
TPZFMatrix<int> fmeshtopology;
TPZFMatrix<double> fmeshcoords;
std::vector<std::vector< std::vector<double > > > fallcoords;
int fOrder;
};


