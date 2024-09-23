#include "readgidmesh.h"

#include <cmath>
#include <set>
#include <iostream>
#include <fstream>
#include <string>
#include "pzgmesh.h"
#include "pzstack.h"
#include "TPZVTKGeoMesh.h"
#include "TPZAnalysis.h"
#include <TPZLinearAnalysis.h> //for TPZLinearAnalysis
#include <TPZSSpStructMatrix.h>
#include "TPZBndCond.h"
#include <pzgeoel.h>
#include "pzgeoelbc.h"
#include "pzfmatrix.h"
#include "pzbstrmatrix.h"

#include <TPZGeoElement.h>

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
#include "TPZSSpStructMatrix.h"
#include "TPZGmshReader.h"
#include <pzgmesh.h> //for TPZGeoMesh
#include <pzcmesh.h> //for TPZCompMesh
#include <TPZGeoMeshTools.h> //for TPZGeoMeshTools::CreateGeoMeshOnGrid
#include <MMeshType.h> //for MMeshType
#include <pzmanvector.h>//for TPZManVector
#include <Poisson/TPZMatPoisson.h> //for TPZMatPoisson
#include <TPZBndCond.h> //for TPZBndCond
#include <TPZLinearAnalysis.h> //for TPZLinearAnalysis
#include <TPZSSpStructMatrix.h> //symmetric sparse matrix storage
#include <pzskylstrmatrix.h> //symmetric skyline matrix storage
#include <pzstepsolver.h> //for TPZStepSolver
#include <TPZSimpleTimer.h>
#include "TPZBndCondT.h"
#include "DarcyFlow/TPZDarcyFlow.h"
#include <pzlog.h>
#include "TPZPardisoSolver.h"
#include "Plasticity/TPZElasticResponse.h"
#include "Plasticity/pzelastoplasticanalysis.h"
#include "Plasticity/TPZYCMohrCoulombPV.h"
#include "Plasticity/TPZMatElastoPlastic_impl.h"
#include "Plasticity/TPZMatElastoPlastic2D.h"
#include "Plasticity/TPZMatElastoPlastic.h"
#include "readgidmesh.h"
#include "slopeconfigure.h"
using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::duration;
using std::chrono::milliseconds;
using std::chrono::seconds;
using namespace std;

class TPZMaterial;
TPZCompMesh * CreateCMeshElastoplastic ( TPZGeoMesh *gmesh, int pOrder );
TPZGeoMesh *  CreateGMeshSlope ( int ref );
TPZGeoMesh * CreateGMesh ( int ref,string file,std::vector<double> coordbc );

#ifdef PZ_LOG
static TPZLogger teste("logmain");
#endif

int main()
{
#ifdef PZ_LOG
    TPZLogger::InitializePZLOG();
#endif

    int ref=2;
    string file;
    file="/home/diogo/projects/neopz-master/Projects2/Plasticity2D/tri-struc-v2.msh";
    std::vector<double> coordbc(3);
    //coordbc[0]=30.;coordbc[1]=5.;coordbc[2]=5.;
    coordbc[0]=75.;coordbc[1]=30.;coordbc[2]=10.;
    //coordbc[0]=75.;coordbc[1]=30.;coordbc[2]=5.;
    //TPZGeoMesh *gmesh = CreateGMesh ( ref,file,coordbc );
    TPZGeoMesh *gmesh =   CreateGMeshSlope (  ref );

    int porder=3;
    REAL gammaagua=0.;
    REAL gammasolo=20.;
    Slope*SlopeManager = new Slope(gmesh,porder,gammaagua,gammasolo);

    //SlopeManager->Solve( );
    SlopeManager->ShearRed ( );
    //string vtk="file.vtk";
    //SlopeManager->PostPlasticity(vtk);

    return 0;
}




TPZGeoMesh *  CreateGMeshSlope ( int ref )
{
    const std::string name ( "Darcy Flow Slope" );

    TPZGeoMesh *gmesh  =  new TPZGeoMesh();

    gmesh->SetName ( name );

    gmesh->SetDimension ( 2 );

    TPZVec<REAL> coord ( 2 );

vector<vector<double>> co= {
    {0., 0.}, {75., 0.}, {75., 30.},{45., 30.},
    {35., 40.}, {0.,40.},{35./3., 40.},{2 * 35/3., 40.},
    {30., 40.},{30., 30.}, {60.,30.},{2* 35./3.,2* 35/3.},
    {45., 2* 35/3.},{35./3., 35/3.}, {60., 35./3.}
};
vector<vector<int>> topol = {
    {0,  1,  14, 13},{1,  2,  10, 14}, {14, 10, 3,  12},
    {13, 14, 12, 11},{11, 12, 3,  9}, {9,  3,  4,  8},
    {11, 9,  8,  7},{13, 11, 7, 6},{0, 13,  6, 5}
};

    gmesh->NodeVec().Resize ( co.size() );

    for ( int inode=0; inode<co.size(); inode++ ) {
        coord[0] = co[inode][0];
        coord[1] = co[inode][1];
        gmesh->NodeVec() [inode] = TPZGeoNode ( inode, coord, *gmesh );
    }
    TPZVec <long> TopoQuad ( 4 );
    for ( int iel=0; iel<topol.size(); iel++ ) {
        TopoQuad[0] = topol[iel][0];
        TopoQuad[1] = topol[iel][1];
        TopoQuad[2] =	topol[iel][2];
        TopoQuad[3] = topol[iel][3];
        new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> ( iel, TopoQuad, 1,*gmesh );
    }




    int id = topol.size();
    TPZVec <long> TopoLine ( 2 );
    TopoLine[0] = 0;
    TopoLine[1] = 1;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -1, *gmesh );//bottom

    id++;
    TopoLine[0] = 1;
    TopoLine[1] = 2;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -2, *gmesh );//rigth

    id++;
    TopoLine[0] = 2;
    TopoLine[1] = 3;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -3, *gmesh );//top-rigth

    id++;
    TopoLine[0] = 4;
    TopoLine[1] = 5;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -4, *gmesh ); //top-left

    id++;
    TopoLine[0] = 3;
    TopoLine[1] = 4;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -6, *gmesh ); //ramp

    id++;
    TopoLine[0] = 5;
    TopoLine[1] = 0;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -5, *gmesh ); //left





    gmesh->BuildConnectivity();
    for ( int d = 0; d<ref; d++ ) {
        int nel = gmesh->NElements();
        TPZManVector<TPZGeoEl *> subels;
        for ( int iel = 0; iel<nel; iel++ ) {
            TPZGeoEl *gel = gmesh->ElementVec() [iel];
            gel->Divide ( subels );
        }
    }

    string meshref = "gmesh.vtk";
    std::ofstream files ( meshref );
    TPZVTKGeoMesh::PrintGMeshVTK ( gmesh,files,true );
    return gmesh;
}

TPZGeoMesh * CreateGMesh ( int ref,string file,std::vector<double> coordbc )
{


    REAL L=coordbc[0];
    REAL h1=coordbc[1];
    REAL h2=coordbc[2];

    TPZGeoMesh *gmesh  =  new TPZGeoMesh();

    gmesh->SetDimension ( 2 );

    readgidmesh read;

    std::vector<std::vector<int>> meshtopol;
    std::vector<std::vector<double>> meshcoords;

    read.ReadMesh2 ( meshtopol,meshcoords,file );


    cout << "a" << endl;
    int ncoords = meshcoords.size();
    gmesh->NodeVec().Resize ( ncoords );

    TPZVec<REAL> coord ( 2 );
    for ( int inode=0; inode<ncoords; inode++ )
    {
        coord[0] = meshcoords[inode][1];
        coord[1] = meshcoords[inode][2];
        gmesh->NodeVec() [inode] = TPZGeoNode ( inode, coord, *gmesh );
    }

    int sz = meshtopol.size();
    TPZVec <long> TopoQuad ( 4 );
    TPZVec <long> TopoTri ( 3 );
    TPZVec <long> TopoLine ( 2 );
    for ( int iel=0; iel<meshtopol.size(); iel++ )
    {
        if ( meshtopol[iel].size() ==4 )
        {
            TopoTri[0] =meshtopol[iel][1]-1;
            TopoTri[1] =meshtopol[iel][2]-1;
            TopoTri[2] =meshtopol[iel][3]-1;
            new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> ( iel, TopoTri, 1,*gmesh );

        }
        else if ( meshtopol[iel].size() ==5 )
        {

            TopoQuad[0] =meshtopol[iel][1]-1;
            TopoQuad[1] =meshtopol[iel][2]-1;
            TopoQuad[2] =meshtopol[iel][3]-1;
            TopoQuad[3] =meshtopol[iel][4]-1;
            new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> ( iel, TopoQuad, 1,*gmesh );

        }
        else if ( meshtopol[iel].size() ==3 )
        {
            TopoLine[0] = meshtopol[iel][1]-1;
            TopoLine[1] = meshtopol[iel][2]-1;
            REAL x0 = meshcoords[TopoLine[0]][1];
            REAL y0 = meshcoords[TopoLine[0]][2];
            REAL xf = meshcoords[TopoLine[1]][1];
            REAL yf = meshcoords[TopoLine[1]][2];
            REAL tol=1.e-3;
            if ( ( fabs ( ( y0-0 ) ) <tol && fabs ( ( yf-0 ) ) <tol ) )
            {
                new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( iel, TopoLine, -1, *gmesh );//bottom
            }
            else if ( ( fabs ( ( x0-L ) ) <tol && fabs ( ( xf-L ) ) <tol ) )
            {
                new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( iel, TopoLine, -2, *gmesh );//rigth
            }
            else if ( ( fabs ( ( y0-h1 ) ) <tol && fabs ( ( yf-h1 ) ) <tol ) )
            {
                new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( iel, TopoLine, -3, *gmesh );//toprigth
            }
            else if ( ( fabs ( ( y0-(h1+h2) ) ) <tol && fabs ( ( yf-(h1+h2) ) ) <tol ) )
            {
                new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( iel, TopoLine, -4, *gmesh );//topleft
            }
            else if ( ( fabs ( ( x0-0 ) ) <tol && fabs ( ( xf-0 ) ) <tol ) )
            {
                new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( iel, TopoLine, -5, *gmesh );//left
            }
            else if ( ( fabs ( ( xf-x0 ) ) >tol && fabs ( ( yf-y0 ) ) >tol ) )
            {
                new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( iel, TopoLine, -6, *gmesh );//ramp
            }
            else
            {
                cout<< "bc element not found."<<endl;
                cout<< "x0 = " << x0 << " y0 = "<< y0 << endl;
                cout<< "xf = " << xf << " yf = "<< yf << endl;
                DebugStop();
            }

        }

    }


    gmesh->BuildConnectivity();
    cout << "c" << endl;
    for ( int d = 0; d<ref; d++ )
    {
        int nel = gmesh->NElements();
        TPZManVector<TPZGeoEl *> subels;
        for ( int iel = 0; iel<nel; iel++ )
        {
            TPZGeoEl *gel = gmesh->ElementVec() [iel];
            gel->Divide ( subels );
        }
    }



    std::ofstream files ( "teste-mesh.vtk" );
    TPZVTKGeoMesh::PrintGMeshVTK ( gmesh,files,false );
    cout << "d" << endl;
    return gmesh;


}
