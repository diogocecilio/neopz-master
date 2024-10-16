
#include "slopeconfigure.h"
#include "SlopeAnalysis.h"
#include "RandomFieldAnalysis.h"
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

TPZCompMesh * CreateCMeshField (  TPZGeoMesh *gmesh,string filename );
void SolveDeter();
void TransferSolutionFrom ( TPZFMatrix<REAL> dataexter,TPZCompMesh *cmesh);
int FindElement ( TPZCompMesh * cmesh,TPZManVector<REAL,3>&vecx, TPZManVector<REAL,3>&vecxi );
void SetSol(TPZInterpolationSpace * intel,TPZMaterialDataT<REAL>& data,TPZFMatrix<REAL> dataexter);
void MonteCarlo();
TPZGeoMesh * TriGMesh(int ref,string file);

TPZGeoMesh *  CreateGMeshSlopeMista ( int ref );
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
int main2()
{
#ifdef PZ_LOG
        TPZLogger::InitializePZLOG();
#endif

        SlopeAnalysis  anal;// =new  slopeanalysis();
        //anal->CreateFields();
        //SolveDeter();
       // MonteCarlo();
        return 0;
}

void MonteCarloFields()
{

}

void MonteCarlo()
{

        string filenamecoes="/home/diogo/projects/neopz-master/Projects2/Plasticity2D/coesao2tri.dat";
        string filenameatrito="/home/diogo/projects/neopz-master/Projects2/Plasticity2D/atrito2tri.dat";

        int ref0=2;//original ref
        int ref1=5;//adaptive ref

        std::vector<double> coordbc ( 3 );
        coordbc[0]=75.;
        coordbc[1]=30.;
        coordbc[2]=10.;

        int porder=1;
        REAL gammaagua=0.;
        REAL gammasolo=20.;
        REAL coes=10.;
        REAL atrito=30.*M_PI/180;

        for ( int imc=750; imc<1000; imc++ )
        {
                TPZGeoMesh *gmesh = TriGMesh(ref0,"malhaplastica");
                Slope*SlopeManager = new Slope ( gmesh,porder,ref1,gammaagua,gammasolo,coes,atrito );
                TPZGeoMesh *gmesh2 = TriGMesh(ref0,"malhacoes");
                TPZGeoMesh *gmesh3 = TriGMesh(ref0,"malhaatrito");

                TPZCompMesh * cmeshfieldatrito= CreateCMeshField ( gmesh2,filenameatrito);
                TPZCompMesh * cmeshfieldcoes = CreateCMeshField ( gmesh3,filenamecoes);

                TPZVec<TPZCompMesh*> cmeshvec(2);

                cmeshvec[0]=cmeshfieldcoes;
                cmeshvec[1]=cmeshfieldatrito;

                SlopeManager->SetCompMeshField(cmeshvec);

                string saidafs = "post/fs";
                string vtk = "postvtk/saidamontecarloplasticity";
                auto var=to_string ( imc );
                saidafs+=var;
                vtk+=var;
                saidafs+=".dat";
                vtk+=".vtk";
                ofstream out ( saidafs );

                std::cout << "imc = " <<  imc << std::endl;

                auto t1 = high_resolution_clock::now();
                REAL FS = SlopeManager->Solve( imc);
                auto t2 = high_resolution_clock::now();
                auto ms_int = duration_cast<seconds> ( t2 - t1 );
                cout <<" tempo total da simulacao   = "<<ms_int.count() << " s " << endl;

                cout << " post processing..." <<endl;
                SlopeManager->PostPlasticity(vtk);
                out << FS << endl;
                delete gmesh;
                delete gmesh2;
                delete gmesh3;
                delete SlopeManager;

        }
}

//Cria uma malha com material generico com um unico grau de liberdade para armazenar a solucao do field nos nos;
//A malha que le deve ser identica a que gerou o field
TPZCompMesh * CreateCMeshField ( TPZGeoMesh *gmesh,string filename )
{
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
        cmesh->LoadSolution(pzdata);
        return cmesh;
}


void SolveDeter()
{

        REAL gammaagua=0.;
        REAL gammasolo=20.;
        REAL coes=10.;
        REAL atrito=30.*M_PI/180;

        std::set<long> elindices,elindices2;
        TPZGeoMesh *gmeshtest = TriGMesh(  2,"teste" );
        Slope*SlopeManagertest = new Slope ( gmeshtest,1,1,gammaagua,gammasolo,coes,atrito );
        SlopeManagertest->ShearRed(20,0.5,0.01);
        SlopeManagertest->ComputeElementDeformation();
        SlopeManagertest->PRefineElementsAbove(0.01, 3,elindices2);
        SlopeManagertest->DivideElementsAbove ( 0.01,elindices );
        SlopeManagertest->InitializeMemory();
        SlopeManagertest->ShearRed(20,0.5,0.01);
        string meshref = "refinidemesh-grid.vtk";
        std::ofstream files ( meshref );
        TPZVTKGeoMesh::PrintGMeshVTK (SlopeManagertest->fCompMesh->Reference(),files,true );
        string vtk0 = "benchmark.vtk";
        SlopeManagertest->PostPlasticity(vtk0);
}

TPZGeoMesh * TriGMesh(int ref,string file)
{
        const std::string name ( file);

        TPZGeoMesh *gmesh  =  new TPZGeoMesh();

        gmesh->SetName ( name );

        gmesh->SetDimension ( 2 );

        TPZVec<REAL> coord ( 2 );

        vector<vector<double>> co= {
                /*0*/{0,0},/*1*/{10,0},/*2*/{20,0},/*3*/{30,0},/*4*/{40,0},/*5*/{50,0},/*6*/{60,0},/*7*/{70,0},
                /*8*/{0,10},/*9*/{10,10},/*10*/{20,10},/*11*/{30,10},/*12*/{40,10},/*13*/{50,10},/*14*/{60,10},/*15*/{70,10},
                /*16*/{0,20},/*17*/{10,20},/*18*/{20,20},/*19*/{30,20},/*20*/{40,20},/*21*/{50,20},/*22*/{60,20},/*23*/{70,20},
                /*24*/{0,30},/*25*/{10,30},/*26*/{20,30},/*27*/{30,30},/*28*/{40,30},/*29*/{50,30},/*30*/{60,30},/*31*/{70,30},
               /*32*/{0,40},/*33*/{10,40},/*34*/{20,40},/*35*/{30,40}
        };
        vector<vector<int>> topol = {

                /*0*/{0,1,8},/*1*/{1,9,8},/*2*/{1,2,9},/*3*/{2,10,9},/*4*/{2,3,10},/*5*/{3,11,10},/*6*/{3,4,11},
                /*7*/{4,12,11},/*8*/{4,5,12}/*9*/,{5,13,12},/*10*/{5,6,13},/*11*/{6,14,13},/*12*/{6,7,14},/*13*/{7,15,14},

                /*14*/{8,9,16},/*15*/{9,17,16},/*16*/{9,10,17},/*17*/{10,18,17},/*18*/{10,11,18},/*19*/{11,19,18},/*20*/{11,12,19},
                /*21*/{12,20,19},/*22*/{12,13,20}/*23*/,{13,21,20},/*24*/{13,14,21},/*25*/{14,22,21},/*26*/{14,15,22},/*27*/{15,23,22},

                /*28*/{16,17,24},/*29*/{17,25,24},/*30*/{17,18,25},/*31*/{18,26,25},/*32*/{18,19,26},/*33*/{19,27,26},/*34*/{19,20,27},
                /*35*/{20,28,27},/*36*/{20,21,28}/*37*/,{21,29,28},/*38*/{21,22,29},/*39*/{22,30,29},/*40*/{22,23,30},/*41*/{23,31,30},

                /*42*/{24,25,32},/*43*/{25,33,32},/*44*/{25,26,33},/*45*/{26,34,33},/*46*/{26,27,34},/*47*/{27,35,34},/*48*/{27,28,35},

                {0,1},{1,2},{2,3},{3,4},{4,5},{5,6},{6,7},/*-1 bottom*/

                {7,15},{15,23},{23,31},/*-2 right*/

                {31,30},{30,29},{29,28},/*-3 top right*/

                {35,34},{34,33},{33,32},/*-4 top left*/

                {32,24},{24,16},{16,8},{8,0},/*-5 left*/

                {28,35}/*-6 ramp*/


        };

        gmesh->NodeVec().Resize ( co.size() );

        for ( int inode=0; inode<co.size(); inode++ ) {
                coord[0] = co[inode][0];
                coord[1] = co[inode][1];
                gmesh->NodeVec() [inode] = TPZGeoNode ( inode, coord, *gmesh );
        }
        TPZVec <long> topotri( 3 );
        TPZVec <long> TopoLine( 2 );
        for ( int iel=0; iel<topol.size(); iel++ ) {
                if(topol[iel].size()==3)
                {
                        topotri[0] = topol[iel][0];
                        topotri[1] = topol[iel][1];
                        topotri[2] = topol[iel][2];
                        new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> ( iel, topotri, 1,*gmesh );
                }else if ( topol[iel].size() ==2 ) {

                        TopoLine[0] = topol[iel][0];
                        TopoLine[1] = topol[iel][1];
                        REAL x0 = co[TopoLine[0]][0];
                        REAL y0 = co[TopoLine[0]][1];
                        REAL xf = co[TopoLine[1]][0];
                        REAL yf = co[TopoLine[1]][1];
                        REAL tol=1.e-3;
                        REAL L=70;
                        REAL h1=30;
                        REAL h2=10;
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
        for ( int d = 0; d<ref; d++ ) {
                int nel = gmesh->NElements();
                TPZManVector<TPZGeoEl *> subels;
                for ( int iel = 0; iel<nel; iel++ ) {
                        TPZGeoEl *gel = gmesh->ElementVec() [iel];
                        gel->Divide ( subels );
                }
        }

        string meshref = "gmeshtri.vtk";
        std::ofstream files ( meshref );
        TPZVTKGeoMesh::PrintGMeshVTK ( gmesh,files,true );
        return gmesh;
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



TPZGeoMesh *  CreateGMeshSlopeMista ( int ref )
{
        const std::string name ( "Malha Mista" );

        TPZGeoMesh *gmesh  =  new TPZGeoMesh();

        gmesh->SetName ( name );

        gmesh->SetDimension ( 2 );

        TPZVec<REAL> coord ( 2 );

        vector<vector<double>> co= {
                {0., 0.}, {35., 0.}, {45.,0.},{75., 0.},
                {75., 30.}, {45.,30.},{35., 30.},{0., 30.},
                {35., 40.},{0., 40.}
        };
        vector<vector<int>> topol = {
                {0,  1,  6, 7},{1,  2,  5, 6}, {2, 3, 4,  5},
                {5, 6, 8},{6, 8, 9,  7}
        };

        gmesh->NodeVec().Resize ( co.size() );

        for ( int inode=0; inode<co.size(); inode++ ) {
                coord[0] = co[inode][0];
                coord[1] = co[inode][1];
                gmesh->NodeVec() [inode] = TPZGeoNode ( inode, coord, *gmesh );
        }
        TPZVec <long> TopoQuad ( 4 ),topoltri( 3 );

        int el=0;
        TopoQuad[0] = topol[el][0];
        TopoQuad[1] = topol[el][1];
        TopoQuad[2] = topol[el][2];
        TopoQuad[3] = topol[el][3];
        new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> ( el , TopoQuad, 1,*gmesh );

        el++;
        TopoQuad[0] = topol[el][0];
        TopoQuad[1] = topol[el][1];
        TopoQuad[2] = topol[el][2];
        TopoQuad[3] = topol[el][3];
        new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> ( el , TopoQuad, 1,*gmesh );

        el++;
        TopoQuad[0] = topol[el][0];
        TopoQuad[1] = topol[el][1];
        TopoQuad[2] = topol[el][2];
        TopoQuad[3] = topol[el][3];
        new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> ( el , TopoQuad, 1,*gmesh );

        el++;
        topoltri[0] = topol[el][0];
        topoltri[1] = topol[el][1];
        topoltri[2] = topol[el][2];
        new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> ( el , topoltri, 1,*gmesh );

        el++;
        TopoQuad[0] = topol[el][0];
        TopoQuad[1] = topol[el][1];
        TopoQuad[2] = topol[el][2];
        TopoQuad[3] = topol[el][3];
        new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> ( el , TopoQuad, 1,*gmesh );




        int id = 4;
        TPZVec <long> TopoLine ( 2 );
        TopoLine[0] = 0;
        TopoLine[1] = 1;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -1, *gmesh );//bottom

        id++;
        TopoLine[0] = 1;
        TopoLine[1] = 2;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -1, *gmesh );//bottom

        id++;
        TopoLine[0] = 2;
        TopoLine[1] = 3;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -1, *gmesh );//bottom

        id++;
        TopoLine[0] = 3;
        TopoLine[1] = 4;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -2, *gmesh ); //rigth

        id++;
        TopoLine[0] = 4;
        TopoLine[1] = 5;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -3, *gmesh ); //toprigth

        id++;
        TopoLine[0] = 8;
        TopoLine[1] = 9;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -4, *gmesh ); //topleft

        id++;
        TopoLine[0] = 9;
        TopoLine[1] = 7;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -5, *gmesh ); //left

        id++;
        TopoLine[0] = 7;
        TopoLine[1] = 0;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -5, *gmesh ); //left

        id++;
        TopoLine[0] = 5;
        TopoLine[1] = 8;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -6, *gmesh ); //ramp






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
