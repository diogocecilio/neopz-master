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
#include "pzeuleranalysis.h"
#include "pzerror.h"
#include "TPZCompElDisc.h"
#include "pzfstrmatrix.h"
#include "TPZParFrontStructMatrix.h"
#include "TPZFrontNonSym.h"
#include "TPZBSpStructMatrix.h"
#include "TPZElementMatrixT.h"
#include "pzbdstrmatrix.h"
#include "pzelmat.h"
#include <time.h>
#include "pzlog.h"
#include "TPZFileStream.h"
#include <TPZBFileStream.h>

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

TPZCompMesh * CreateCMeshField ( TPZGeoMesh *gmesh);

void TransferSolutionFrom ( TPZFMatrix<REAL> dataexter,TPZCompMesh *cmesh);
int FindElement ( TPZCompMesh * cmesh,TPZManVector<REAL,3>&vecx, TPZManVector<REAL,3>&vecxi );
void SetSol(TPZInterpolationSpace * intel,TPZMaterialDataT<REAL>& data,TPZFMatrix<REAL> dataexter);
//1 Ler os dados e colocar em matrizes cheias
// x1 y1 z1 dado1 dado2 ... dado n
// x2 y2 z2 dado1 dado2 ... dado n
// xn yn zn dado1 dado2 ... dado n
//opcao 1 - Percorrer as matrizes de dados buscando as coordenadas correspondentes à malha atual considerando um percentual de erro
//opcao 2 - Utilziar a mesma malha empregada para gerar os dados para interpolar a solucao

#ifdef PZ_LOG
static TPZLogger teste ( "logmain" );
#endif
//#include <bits/stdc++.h>
using namespace std;
int main()
{
#ifdef PZ_LOG
        TPZLogger::InitializePZLOG();
#endif

        int ref0=2;
        string file;
        file="/home/diogo/projects/neopz-master/Projects2/Plasticity2D/tri-struc-v2.msh";
        std::vector<double> coordbc ( 3 );
        //coordbc[0]=30.;coordbc[1]=5.;coordbc[2]=5.;
        coordbc[0]=75.;
        coordbc[1]=30.;
        coordbc[2]=10.;
        //coordbc[0]=75.;coordbc[1]=30.;coordbc[2]=5.;
        TPZGeoMesh *gmesh = CreateGMesh ( ref0+1,file,coordbc );
       // TPZGeoMesh *gmesh =   CreateGMeshSlope ( ref0 );

        //TPZCompMesh * cmesh = new TPZCompMesh(gmesh);

        int ref1=2;
        int porder=2;
        REAL gammaagua=0.;
        REAL gammasolo=20.;
        Slope*SlopeManager = new Slope ( gmesh,porder,ref1,gammaagua,gammasolo );

        TPZGeoMesh *gmesh2 = CreateGMesh ( ref0,file,coordbc );
        TPZCompMesh * cmeshfield = CreateCMeshField ( gmesh2);
        SlopeManager->fCompMeshField=cmeshfield;



//
        auto t1 = high_resolution_clock::now();

        SlopeManager->Solve( );
        auto t2 = high_resolution_clock::now();
        auto ms_int = duration_cast<seconds> ( t2 - t1 );
        //std::cout << "tempo total iterative process = "<<ms_int.count() << " s\n";
        //conv  = anal->IterativeProcess ( cout, tol2, NumIter,linesearch,checkconv,iters);
        cout <<" tempo total da simulacao   = "<<ms_int.count() << " s " << endl;
         //conv = anal->FindRoot ( );



        //SlopeManager->ShearRed();

        //SlopeManager->PostPlasticity("saida0.vtk");
//     cout << "Refining.."<<endl;
//     std::set<long> elindices;
//     for (int iref=1; iref<=0; iref++ )
//     {
//         Slope*SlopeManagertemp = new Slope(gmesh,porder,gammaagua,gammasolo);
//         SlopeManagertemp->ShearRed();
//         SlopeManagertemp->ComputeElementDeformation();
//         SlopeManagertemp->DivideElementsAbove ( 0.01,elindices );
//         gmesh=SlopeManagertemp->fCompMesh->Reference();
//     }
//
//     Slope*SlopeManager = new Slope(gmesh,porder,gammaagua,gammasolo);
//
//     cout << "Solving again.."<<endl;
//     //fCompMesh = CreateCMeshElastoplastic ( fGmesh, fPorder );
//     SlopeManager->ShearRed();
//     SlopeManager->PostPlasticity("saida3.vtk");




//
//     TPZCompMesh  *cmesh = SlopeManager->fCompMesh;
//     TPZManVector<REAL,3> vecx(3);
//     TPZManVector<REAL,3> vecxi(2);
//     vecx[0]=40.5;
//     vecx[1]=29.7;
//     vecx[2]=0.;
//     FindElement(cmesh,vecx,vecxi);




//     {
//         TPZBFileStream save;
//         save.OpenWrite("Config1-0.bin");
//         SlopeManager->Write(save);
//     }
//
//     {
//         TPZBFileStream read;
//         read.OpenRead("Config1-0.bin");
//         SlopeManager->Read(read);
//     }

        return 0;
}


TPZCompMesh * CreateCMeshField ( TPZGeoMesh *gmesh )
{
        string filename="/home/diogo/projects/neopz-master/Projects2/Plasticity2D/coesao3.dat";
        readgidmesh read;
        TPZFMatrix<REAL> pzdata  = read.ReadData(filename);
        int dim=2;
        TPZCompMesh * cmesh = new TPZCompMesh ( gmesh );
        cmesh->SetDefaultOrder ( 1);
        cmesh->SetDimModel ( dim );
        int matid=1;
        TPZDarcyFlow* material = new TPZDarcyFlow (matid,  dim);
        material->SetId(matid);
        cmesh->InsertMaterialObject ( material );
        cmesh->SetAllCreateFunctionsContinuous();
        cmesh->AutoBuild();

        int eqs = cmesh->NEquations();
        //cout << "eqs = " << eqs << endl;
        //cout << "pzdata.Rows() = " << pzdata.Rows() << endl;
        //cout << "pzdata.Cols() = " << pzdata.Cols() << endl;

        cmesh->LoadSolution(pzdata);
       // TransferSolutionFrom ( pzdata,cmesh);
        return cmesh;
}

void TransferSolutionFrom ( TPZFMatrix<REAL> dataexter,TPZCompMesh *cmesh)
{
    int rows=dataexter.Rows();
    int cols=dataexter.Cols();

    int nels =  cmesh->NElements();

    TPZGeoMesh *gmesh = cmesh->Reference();

    REAL tol=1.e-6;

    for(int iel=0;iel<nels;iel++)
    {
            TPZCompEl * cel =  cmesh->Element(iel);

            TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *> ( cel );

            TPZMaterialDataT<REAL> data;
            intel->InitMaterialData ( data );
            data.fNeedsSol = true;

            int numbersol= dataexter.Cols();
            data.sol.resize ( numbersol );
            for ( int is = 0; is<numbersol; is++ ){
                    data.sol[is].Resize ( 1 );
                    data.sol[is].Fill ( 0. );
        }


            TPZGeoEl * gel =cel->Reference();

            TPZVec<int64_t> nodeinds;
            gel->GetNodeIndices(nodeinds);

            int nnodes=nodeinds.size();

            int nconnects = intel->NConnects();

             cout << "nodeinds[0]  = "<<nodeinds[0]<< endl;
            cout << "intel->ConnectIndex(0)  = "<<intel->ConnectIndex(0)<< endl;

            for(int inode=0;inode<nconnects;inode++)
            {
                TPZConnect *df = &intel->Connect(inode);
                TPZGeoNode node = gmesh->NodeVec()[inode];

                TPZManVector<REAL,3> co(3);
                node.GetCoordinates(co);
                for(int irow=0;irow<rows;irow++)
                {
                        REAL xdata = dataexter(irow,0);
                        REAL ydata = dataexter(irow,1);
                        REAL zdata = dataexter(irow,2);
                        if(fabs(xdata-co[0])<tol&&fabs(ydata-co[1])<tol&&fabs(zdata-co[2])<tol)
                        {

                                for(int icol=0;icol<cols;icol++)
                                {
                                        data.sol[icol][0]=dataexter(irow,icol);
                                }

                                intel->LoadSolution();

                        }

                }

            }

    }


}

// void TransferSolutionFrom ( TPZFMatrix<REAL> dataexter,TPZCompMesh *cmesh)
// {
//     std::vector<std::vector<double>> datastd;
//
//
//     TPZGeoMesh *gmesh = cmesh->Reference();
//
//     TPZAdmChunkVector<TPZGeoNode> & nodevec =  gmesh->NodeVec();
//
//     int nnodes = nodevec.NElements();
//
//     for(int inode=0;inode<nnodes;inode++)
//     {
//         TPZManVector<REAL,3> co(3),xi(2);
//         nodevec[inode].GetCoordinates(co);
//         int el = FindElement ( cmesh,co,xi );
//
//         //cout << "el  = " << el <<endl;
//         TPZCompEl *cel =  cmesh->ElementVec()[el];
//         if(!cel)
//         {
//             continue;
//         }
//         TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *> ( cel );
//         if(!intel)
//         {
//             continue;
//         }
//         TPZMaterialDataT<REAL> data;
//         intel->InitMaterialData ( data );
//         data.fNeedsSol = true;
//
//         data.intLocPtIndex = inode;
//         intel->ComputeRequiredData ( data, xi );
//
//         //intel->ComputeSolution(xi,data);
//
//
//
//         REAL x,y,z,solu;
//         x=data.x[0];
//         y=data.x[1];
//         z=data.x[2];
//
//
//         REAL xext = dataexter(inode,0);
//         REAL yext = dataexter(inode,1);
//         REAL zext = dataexter(inode,2);
//         if(fabs(xext-x)>1.e-3||fabs(yext-y)>1.e-3||fabs(zext-z)>1.e-3)
//         {
//             cout << "as coordendas dos dados internos nao batem com a mablha interna."<<endl;
//             DebugStop();
//         }
//
//        // SetSol(intel,data,dataexter);
//        // intel->LoadSolution();
//     }
//
// }
void SetSol(TPZInterpolationSpace * intel,TPZMaterialDataT<REAL>& data,TPZFMatrix<REAL> dataexter)
{
        const TPZFMatrix<REAL> &phi = data.phi;

	const int nstate = intel->Material()->NStateVariables();
	const int ncon = intel->NConnects();
	TPZBlock &block = intel->Mesh()->Block();
	//TPZFMatrix<REAL> &MeshSol = intel->Mesh()->Solution();
    //const int64_t numbersol = MeshSol.Cols();

        int numbersol= dataexter.Cols()-3;
        data.sol.resize(numbersol);
    for (int is = 0; is<numbersol; is++) {
        data.sol[is].Resize(1);
        data.sol[is].Fill(0.);
    }

        	int64_t iv = 0;
	for(int in=0; in<ncon; in++) {
		TPZConnect *df = &intel->Connect(in);
		const int64_t dfseq = df->SequenceNumber();
		const int dfvar = block.Size(dfseq);
		const int64_t pos = block.Position(dfseq);
		for(int jn=0; jn<dfvar; jn++) {
            for (int64_t is=0; is<numbersol-3; is++) {
                data.sol[is][iv%nstate] +=
                    (REAL)phi.Get(iv/nstate,0)*dataexter(pos+jn,is+3);
            }
			iv++;
		}
	}
}
int FindElement ( TPZCompMesh * cmesh,TPZManVector<REAL,3>&vecx, TPZManVector<REAL,3>&vecxi )
{
        int nels = cmesh->NElements();

        REAL tol=1.e-5;

        for ( int iel=0; iel<nels; iel++ ) {

                TPZCompEl * cel =  cmesh->ElementVec() [iel];
                if ( !cel ) {
                        continue;
                }
                TPZGeoEl * gel = cel->Reference();
                if ( !gel ) {
                        continue;
                }
                bool find = gel->ComputeXInverse ( vecx,vecxi,tol );
                if ( find ) {
                        //cout << "ponto encontrado no elemento id "<<iel <<endl;
                        //cout << "coordenadas no elemento mestre = "<< vecxi << endl;
                        return iel;
                } else {
                        //cout << " ponto NAO encontrado! " <<endl;
                        //DebugStop();
                }

        }

            cout << " int KLAnalysis::FindElement  ponto NAO encontrado! " <<endl;
            DebugStop();

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
        for ( int inode=0; inode<ncoords; inode++ ) {
                coord[0] = meshcoords[inode][1];
                coord[1] = meshcoords[inode][2];
                gmesh->NodeVec() [inode] = TPZGeoNode ( inode, coord, *gmesh );
        }

        int sz = meshtopol.size();
        TPZVec <long> TopoQuad ( 4 );
        TPZVec <long> TopoTri ( 3 );
        TPZVec <long> TopoLine ( 2 );
        for ( int iel=0; iel<meshtopol.size(); iel++ ) {
                if ( meshtopol[iel].size() ==4 ) {
                        TopoTri[0] =meshtopol[iel][1]-1;
                        TopoTri[1] =meshtopol[iel][2]-1;
                        TopoTri[2] =meshtopol[iel][3]-1;
                        new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> ( iel, TopoTri, 1,*gmesh );

                } else if ( meshtopol[iel].size() ==5 ) {

                        TopoQuad[0] =meshtopol[iel][1]-1;
                        TopoQuad[1] =meshtopol[iel][2]-1;
                        TopoQuad[2] =meshtopol[iel][3]-1;
                        TopoQuad[3] =meshtopol[iel][4]-1;
                        new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> ( iel, TopoQuad, 1,*gmesh );

                } else if ( meshtopol[iel].size() ==3 ) {
                        TopoLine[0] = meshtopol[iel][1]-1;
                        TopoLine[1] = meshtopol[iel][2]-1;
                        REAL x0 = meshcoords[TopoLine[0]][1];
                        REAL y0 = meshcoords[TopoLine[0]][2];
                        REAL xf = meshcoords[TopoLine[1]][1];
                        REAL yf = meshcoords[TopoLine[1]][2];
                        REAL tol=1.e-3;
                        if ( ( fabs ( ( y0-0 ) ) <tol && fabs ( ( yf-0 ) ) <tol ) ) {
                                new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( iel, TopoLine, -1, *gmesh );//bottom
                        } else if ( ( fabs ( ( x0-L ) ) <tol && fabs ( ( xf-L ) ) <tol ) ) {
                                new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( iel, TopoLine, -2, *gmesh );//rigth
                        } else if ( ( fabs ( ( y0-h1 ) ) <tol && fabs ( ( yf-h1 ) ) <tol ) ) {
                                new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( iel, TopoLine, -3, *gmesh );//toprigth
                        } else if ( ( fabs ( ( y0- ( h1+h2 ) ) ) <tol && fabs ( ( yf- ( h1+h2 ) ) ) <tol ) ) {
                                new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( iel, TopoLine, -4, *gmesh );//topleft
                        } else if ( ( fabs ( ( x0-0 ) ) <tol && fabs ( ( xf-0 ) ) <tol ) ) {
                                new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( iel, TopoLine, -5, *gmesh );//left
                        } else if ( ( fabs ( ( xf-x0 ) ) >tol && fabs ( ( yf-y0 ) ) >tol ) ) {
                                new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( iel, TopoLine, -6, *gmesh );//ramp
                        } else {
                                cout<< "bc element not found."<<endl;
                                cout<< "x0 = " << x0 << " y0 = "<< y0 << endl;
                                cout<< "xf = " << xf << " yf = "<< yf << endl;
                                DebugStop();
                        }

                }

        }


        gmesh->BuildConnectivity();
        cout << "c" << endl;
        for ( int d = 0; d<ref; d++ ) {
                int nel = gmesh->NElements();
                TPZManVector<TPZGeoEl *> subels;
                for ( int iel = 0; iel<nel; iel++ ) {
                        TPZGeoEl *gel = gmesh->ElementVec() [iel];
                        gel->Divide ( subels );
                }
        }



        std::ofstream files ( "teste-mesh.vtk" );
        TPZVTKGeoMesh::PrintGMeshVTK ( gmesh,files,false );
        cout << "d" << endl;
        return gmesh;


}
