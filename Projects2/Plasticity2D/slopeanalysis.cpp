// SPDX-FileCopyrightText: 2024 <copyright holder> <email>
// SPDX-License-Identifier: Apache-2.0

#include "slopeanalysis.h"
slopeanalysis::slopeanalysis( ):
Slope()
{

}
slopeanalysis::slopeanalysis( TPZGeoMesh * gmesh,int porder,int ref, REAL gammaagua, REAL gammasolo,REAL coes,REAL atrito):
Slope( gmesh, porder, ref, gammaagua, gammasolo, coes,atrito)
{

}

slopeanalysis::slopeanalysis ( const slopeanalysis& other ):Slope(other)
{

}

slopeanalysis::~slopeanalysis()
{

}

void slopeanalysis::CreateFields()
{
        int ref0=0;
        TPZGeoMesh * gmesh =  TriGMesh(ref0);
        TPZCompMesh* fCompMeshSlope = CreateCompMeshKL ( gmesh );
        TPZRandomFieldAnalysis * fFieldAnalysis = new TPZRandomFieldAnalysis ( fCompMeshSlope );

        TPZManVector<std::string> scalarnames ( 7 ), vecnames;
        scalarnames[0] = "vec";
        scalarnames[1] = "vec1";
        scalarnames[2] = "vec2";
        scalarnames[3] = "vec3";
        scalarnames[4] = "vec4";
        scalarnames[5] = "vec5";
        scalarnames[6] = "vecsqr";

        int Startfrom=0;
        if ( Startfrom ==0 ) {

                fFieldAnalysis->SetNEigenpairs ( 250 );
                fFieldAnalysis->Solve();

                TPZBFileStream save;
                save.OpenWrite ( "Config1-0.bin" );
                fFieldAnalysis->Write ( save,fFieldAnalysis->ClassId() );

        }

        if ( Startfrom ==1 ) {

                TPZBFileStream read;
                read.OpenRead ( "Config1-0.bin" );
                fFieldAnalysis->Read ( read,0 );

                fFieldAnalysis->LoadSolution();
                fFieldAnalysis->DefineGraphMesh ( 2,scalarnames,vecnames,"filename2.vtk" );
                fFieldAnalysis->PostProcess ( 0 );

                TPZVec<REAL> meanvec(2);
                meanvec[0]=10;
                meanvec[1]=30 *M_PI/180.;
                TPZVec<REAL> covvec(2);
                covvec[0]=0.3;
                covvec[1]=0.2;
                int samples=1000;
                fFieldAnalysis->SetFieldsData( meanvec,covvec,  samples);

                fFieldAnalysis->ManageFieldCretion();

                TPZBFileStream save;
                save.OpenWrite ( "Config2-0.bin" );
                fFieldAnalysis->Write ( save,fFieldAnalysis->ClassId() );
        }

        if ( Startfrom ==2 ) {
                //TPZRandomFieldAnalysis * dymmyanalysis = new TPZRandomFieldAnalysis(cmesh);

                TPZBFileStream read;
                read.OpenRead ( "Config2-0.bin" );
                fFieldAnalysis->Read ( read,0 );

                TPZVec<TPZFMatrix<REAL>> fields;

                fields = fFieldAnalysis->GetFields();

                fFieldAnalysis->LoadSolution ( fields[1] );
                fFieldAnalysis->DefineGraphMesh ( 2,scalarnames,vecnames,"fields[1].vtk" );
                fFieldAnalysis->PostProcess ( 0 );
        }

}

TPZCompMesh* slopeanalysis::CreateCompMeshKL ( TPZGeoMesh * gmesh )
{

        int id=1;
        int dim = gmesh->Dimension();
        TPZCompMesh * cmesh = new TPZCompMesh ( gmesh );
        REAL lz=1.;
        REAL Lx=20;
        REAL Ly=2;
        REAL Lz=1;
        int type=3;
        int porder=1;
        //TPZKarhunenLoeveMat(int id,REAL Lx,REAL Ly,REAL Lz,int dim, int type,int expansionorder);
        TPZKarhunenLoeveMat * mat = new TPZKarhunenLoeveMat ( id,Lx,Ly,Lz,dim,type );

        cmesh->SetDefaultOrder ( porder );
        cmesh->SetDimModel ( dim );
        cmesh->InsertMaterialObject ( mat );
        cmesh->SetAllCreateFunctionsContinuous();
        cmesh->AutoBuild();

        return cmesh;
}


TPZGeoMesh * slopeanalysis::TriGMesh ( int ref )
{

        TPZGeoMesh *gmesh  =  new TPZGeoMesh();


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
        TPZVec <long> topotri ( 3 );
        TPZVec <long> TopoLine ( 2 );
        for ( int iel=0; iel<topol.size(); iel++ ) {
                if ( topol[iel].size() ==3 ) {
                        topotri[0] = topol[iel][0];
                        topotri[1] = topol[iel][1];
                        topotri[2] = topol[iel][2];
                        new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> ( iel, topotri, 1,*gmesh );
                } else if ( topol[iel].size() ==2 ) {

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
