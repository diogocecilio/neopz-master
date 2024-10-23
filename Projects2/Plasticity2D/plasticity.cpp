
#include "slopeconfigure.h"
#include "SlopeAnalysis.h"
#include <fstream>
void SolveSlope ( int Startfrom );
TPZGeoMesh * TriGMesh ( int ref );
TPZCompMesh* CreateCompMeshKL ( TPZGeoMesh * gmesh,int porder,REAL Lx, REAL Ly, REAL Lz, int id,int type );
TPZCompMesh * CreateCMeshElastoplastic ( TPZGeoMesh *gmesh, int pOrder, REAL coes,REAL atrito );

void saveVector(const std::vector<double>& vec, const std::string& filename) ;

void readVector( std::vector<double>& vec,const std::string& filename);

std::vector<int> GetIndex ( std::vector<double> vetor );

void SolveSlopeIS ( int Startfrom );
int main()
{

        int Startfrom =2;
        SolveSlopeIS ( Startfrom );


        return 0;
}

void SolveSlopeIS ( int Startfrom )
{
        //create random analysis
        int ref=3;
        TPZGeoMesh * gmesh =  TriGMesh ( ref );
        int porder=1;
        REAL Lx=20.;
        REAL Ly=2;
        REAL Lz=1.;
        int id=1;
        int type=3;
        TPZCompMesh * cmesh = CreateCompMeshKL ( gmesh, porder, Lx,  Ly,  Lz,  id, type );
        TPZRandomFieldAnalysis * randonanalysis = new TPZRandomFieldAnalysis ( cmesh );
        TPZManVector<std::string> scalarnames = {"vec","vec1","vec2","vec3","vec4"}, vecnames;

        //create slope analysis
        int ref0slope=1;
        int porderslope=1;
        REAL gammaagua=0.;
        REAL gammasolo=20.;
        REAL coes=10.;
        REAL atrito=30*M_PI/180.;
        REAL fsfail=1.2;
        REAL coesh=coes/fsfail;
        REAL atritoh= atan ( tan ( atrito ) /fsfail );

        SlopeAnalysis  * slopeanalysisf =  new SlopeAnalysis ( gammaagua,gammasolo,coes,atrito,ref0slope,porderslope );

        SlopeAnalysis  * slopeanalysish =  new SlopeAnalysis ( gammaagua,gammasolo,coesh,atritoh,ref0slope,porderslope );

        int solvertype=0;
        int numthreads=12;
 //        slopeanalysish->SetSlopeAnalysis ( solvertype,numthreads );
//
 //        slopeanalysish->SolveDeterministic();
//
//         slopeanalysish->PostPlasticity("deterministic.vtk");

        //solve generalized eigenvalu problem
        if ( Startfrom ==0 ) {
                randonanalysis->SetNEigenpairs ( 500 );
                //randonanalysis->Assemble();
                randonanalysis->Solve();

                //save sqrt(lambda)*phi
                TPZBFileStream save;
                save.OpenWrite ( "Config1-0.bin" );
                randonanalysis->Write ( save,randonanalysis->ClassId() );
                //randonanalysis->LoadSolution();
                randonanalysis->DefineGraphMesh ( 2,scalarnames,vecnames,"filename2Assemble.vtk" );
                randonanalysis->PostProcess ( 0 );
        }


        if ( Startfrom >0 ) {
                //read sqrt(lambda)*phi
                TPZBFileStream read;
                read.OpenRead ( "Config1-0.bin" );
                randonanalysis->Read ( read,0 );

                //seting fied data
                TPZVec<REAL> meanvec ( 2 );
                meanvec[0]=coes;
                meanvec[1]=atrito;
                TPZVec<REAL> covvec ( 2 );
                covvec[0]=0.3;
                covvec[1]=0.2;
                int samples=1000;
                slopeanalysisf->SetFieldsData ( cmesh,randonanalysis->GetSolutionValVec(), meanvec,covvec,  samples );

                meanvec[0]=coesh;
                meanvec[1]=atritoh;
                covvec[0]=0.3;
                covvec[1]=0.2;
                slopeanalysish->SetFieldsData ( cmesh,randonanalysis->GetSolutionValVec(), meanvec,covvec,  samples );

                //crate random fields
                if ( Startfrom==1 ) {
                        //SensivityAnalysisOnEigenvalues ( slopeanalysis);
                        TPZFMatrix<REAL>samples1 = slopeanalysisf->CreateNormalStandardSamples();
                        TPZFMatrix<REAL>samples2 = slopeanalysisf->CreateNormalStandardSamples();
                        TPZVec<TPZFMatrix<REAL>> samples ( 2 );
                        samples[0]=samples1;
                        samples[1]=samples2;
                        slopeanalysisf->SetFieldsSamples ( samples );
                        slopeanalysish->SetFieldsSamples ( samples );

                        slopeanalysisf->ManageFieldCretion();
                        slopeanalysish->ManageFieldCretion();



                        TPZBFileStream save,saveh;
                        save.OpenWrite ( "Config2-0.bin" );
                        saveh.OpenWrite ( "Configh-2-0.bin" );
                        slopeanalysisf->Write ( save,slopeanalysisf->ClassId() );
                        slopeanalysish->Write ( saveh,slopeanalysish->ClassId()+1 );

                } else { //Startfrom>1 solve monte carlo
                        TPZBFileStream read,readh;
                        read.OpenRead ( "Config2-0.bin" );
                        readh.OpenRead ( "Configh-2-0.bin" );
                        slopeanalysisf->Read ( read,0 );
                        slopeanalysish->Read ( readh,0 );

                        //slopeanalysish->TransferFieldsSolutionFrom ( 0 );
                        //slopeanalysish->PostPlasticity ( "fieldtesteH111111111.vtk" );
                        //ofstream out3 ( "saidafsex750-1000.txt" );

                        TPZStack<REAL> data;
                        TPZBFileStream savestack;
                        savestack.OpenWrite ( "stack.txt" );
                        data.Push(1000);
                        savestack.Write(data);

                        std::vector<double> vecout,vecout2,vecin,vecin2;
                        //vecout.push_back(1000);
                        //saveVector(vecout,"out.bin") ;

                        //readVector(vecin,"out.bin") ;

                        //cout<< vecin[0] << endl;
                        //if(true)
                        //{
                                readVector(vecin,"FSf.bin");
                                readVector(vecin2,"FSh.bin");
                                ofstream out ( "fsf.txt" );
                                ofstream out2 ( "fsh.txt" );
                                for(int i=0;i<vecin.size();i++)out<< vecin[i] <<endl;
                                for(int i=0;i<vecin2.size();i++)out2<< vecin2[i] <<endl;

                        //}

                        for ( int imc=0; imc<=1000; imc++ ) {

                                SlopeAnalysis* slopeanalysisf1 = new SlopeAnalysis ( gammaagua,gammasolo,coes,atrito,ref0slope,porderslope );
                                slopeanalysisf1->SetFieldsData ( cmesh,randonanalysis->GetSolutionValVec(), meanvec,covvec,  samples );
                                slopeanalysisf1->SetFields ( slopeanalysisf->GetFields() );
                                slopeanalysisf1->SetSlopeAnalysis ( solvertype,numthreads );
                                slopeanalysisf1->SetFieldsSamples ( slopeanalysisf->GetFieldsSamples() );

                                SlopeAnalysis* slopeanalysish1 = new SlopeAnalysis ( gammaagua,gammasolo,coesh,atritoh,ref0slope,porderslope  );
                                slopeanalysish1->SetFieldsData ( cmesh,randonanalysis->GetSolutionValVec(), meanvec,covvec,  samples );
                                slopeanalysish1->SetFields ( slopeanalysish->GetFields() );
                                slopeanalysish1->SetSlopeAnalysis ( solvertype,numthreads );
                                slopeanalysish1->SetFieldsSamples ( slopeanalysish->GetFieldsSamples() );

                                string saidafs = "post/fs";
                                string vtk = "postvtk/postplasticityvtk";
                                auto var=to_string ( imc );
                                saidafs+=var;
                                vtk+=var;
                                auto saidafs2=saidafs;
                                saidafs+=".dat";
                                saidafs2+=".dat";
                                vtk+=".vtk";



                                ofstream out ( saidafs );
                                ofstream out2 ( saidafs2 );

                                REAL FSf = slopeanalysisf1->SolveSingleField ( imc );

                                REAL FSh = slopeanalysish1->SolveSingleField ( imc );

                                //slopeanalysisf1->PostPlasticity ( vtk );

                                //slopeanalysish1->PostPlasticity ( vtk2 );
                                //out  << FSf << endl;
                                //out2 <<  FSh <<  endl;

                                //if(vecout.size()>vecin.size())
                                //{
                                        vecout.push_back(FSf);
                                        vecout2.push_back(FSh);
                                        saveVector(vecout,"FSf.bin");
                                        saveVector(vecout,"FSh.bin");
                                //}

                        }
                }

        }



}

void SolveSlope ( int Startfrom )
{
        //create random analysis
        int ref=3;
        TPZGeoMesh * gmesh =  TriGMesh ( ref );
        int porder=1;
        REAL Lx=20.;
        REAL Ly=2;
        REAL Lz=1.;
        int id=1;
        int type=3;
        TPZCompMesh * cmesh = CreateCompMeshKL ( gmesh, porder, Lx,  Ly,  Lz,  id, type );
        TPZRandomFieldAnalysis * randonanalysis = new TPZRandomFieldAnalysis ( cmesh );
        TPZManVector<std::string> scalarnames = {"vec","vec1","vec2","vec3","vec4"}, vecnames;

        //create slope analysis
        int ref0slope=3;
        int porderslope=1;
        REAL gammaagua=0.;
        REAL gammasolo=20.;
        //REAL coes=50.;
       // REAL atrito=20*M_PI/180.;
        REAL fsfail=2.23;
        REAL coes=50./fsfail;
        REAL atrito= atan ( tan ( 20*M_PI/180 ) /fsfail );

        SlopeAnalysis  * slopeanalysis =  new SlopeAnalysis ( gammaagua,gammasolo,coes,atrito,ref0slope,porderslope );

        int solvertype=0;
        int numthreads=12;
        slopeanalysis->SetSlopeAnalysis ( solvertype,numthreads );

        slopeanalysis->SolveDeterministic();

        slopeanalysis->PostPlasticity("deterministic.vtk");
        return;

        //solve generalized eigenvalu problem
        if ( Startfrom ==0 ) {
                randonanalysis->SetNEigenpairs ( 500 );
                //randonanalysis->Assemble();
                randonanalysis->Solve();

                //save sqrt(lambda)*phi
                TPZBFileStream save;
                save.OpenWrite ( "Config1-0.bin" );
                randonanalysis->Write ( save,randonanalysis->ClassId() );
                //randonanalysis->LoadSolution();
                randonanalysis->DefineGraphMesh ( 2,scalarnames,vecnames,"filename2Assemble.vtk" );
                randonanalysis->PostProcess ( 0 );
        }


        if ( Startfrom >0 ) {
                //read sqrt(lambda)*phi
                TPZBFileStream read;
                read.OpenRead ( "Config1-0.bin" );
                randonanalysis->Read ( read,0 );

                //seting fied data
                TPZVec<REAL> meanvec ( 2 );
                meanvec[0]=50.;
                meanvec[1]=20. *M_PI/180.;
                TPZVec<REAL> covvec ( 2 );
                covvec[0]=0.3;
                covvec[1]=0.2;
                int samples=1000;

                slopeanalysis->SetFieldsData ( cmesh,randonanalysis->GetSolutionValVec(), meanvec,covvec,  samples );

                //crate random fields
                if ( Startfrom==1 ) {
                        //SensivityAnalysisOnEigenvalues ( slopeanalysis);
                        //TPZFMatrix<REAL>samples1 = slopeanalysis->CreateNormalStandardSamples();
                        //TPZFMatrix<REAL>samples2 = slopeanalysis->CreateNormalStandardSamples();
                        //TPZVec<TPZFMatrix<REAL>> samples ( 2 );
                        //samples[0]=samples1;
                        //samples[1]=samples2;
                        //slopeanalysis->SetFieldsSamples ( samples );
                        //std::vector<int> data = {2,3,4};
                        //slopeanalysis->ManageFieldCretion (  );
                        //slopeanalysis->TransferFieldsSolutionFrom ( 0 );
                        //slopeanalysis->PostPlasticity ( "fieldteste.vtk" );
                        //slopeanalysis->ManageFieldCretion();
                        //slopeanalysis->ComputeH();
                        TPZBFileStream save;
                        save.OpenWrite ( "Config2-0.bin" );
                        slopeanalysis->Write ( save,slopeanalysis->ClassId() );
                } else { //Startfrom>1 solve monte carlo
                        TPZBFileStream read;
                        read.OpenRead ( "Config2-0.bin" );
                        slopeanalysis->Read ( read,0 );



                        //slopeanalysis->SelectCriticalIndexes();
                        //slopeanalysis->ComputeH();
                        ofstream out3 ( "saidafs.txt" );
                        for ( int imc=0; imc<=100; imc++ ) {

                                SlopeAnalysis* slopeanalysis2 = new SlopeAnalysis ( gammaagua,gammasolo,coes,atrito,ref0slope,porderslope ); // Outra cópia independente
                                slopeanalysis2->SetFieldsData ( cmesh,randonanalysis->GetSolutionValVec(), meanvec,covvec,  samples );
                                slopeanalysis2->SetFields ( slopeanalysis->GetFields() );
                                slopeanalysis2->SetSlopeAnalysis ( solvertype,numthreads );
                                slopeanalysis2->SetFieldsSamples ( slopeanalysis->GetFieldsSamples() );

                                string saidafs = "post/fsh";
                                string vtk = "postvtk/postplasticityvtk";
                                auto var=to_string ( imc );
                                saidafs+=var;
                                vtk+=var;
                                saidafs+=".dat";

                                vtk+=".vtk";
                                ofstream out ( saidafs );

                                std::cout << "imc = " <<  imc << std::endl;
                                auto t1 = high_resolution_clock::now();
                                REAL FS2 = slopeanalysis2->SolveSingleField ( imc );
                                auto t2 = high_resolution_clock::now();
                                auto ms_int = duration_cast<seconds> ( t2 - t1 );
                                cout <<" tempo total da simulacao   = "<<ms_int.count() << " s " << endl;
                                cout << " post processing..." <<endl;
                                slopeanalysis2->PostPlasticity ( vtk );
                                //out << FS << endl;
                                out3<<"imc = "<< imc <<endl;
                                out3<<"FSF = "<< FS2 <<endl;
                                //delete slopeanalysis2;
                                //delete slopeanalysis3;

                        }
                }

        }



}



std::vector<int> GetIndex ( std::vector<double> vetor )
{
        std::vector<int> indices ( vetor.size() );
        for ( int i = 0; i < vetor.size(); ++i ) {
                indices[i] = i;
        }

        std::sort ( indices.begin(), indices.end(), [&] ( int a, int b ) {
                return vetor[a] > vetor[b]; // Ordena de forma decrescente
        } );

        return indices;

}

TPZCompMesh* CreateCompMeshKL ( TPZGeoMesh * gmesh,int porder,REAL Lx, REAL Ly, REAL Lz, int id,int type )
{

        int dim = gmesh->Dimension();
        TPZCompMesh * cmesh = new TPZCompMesh ( gmesh );
        TPZKarhunenLoeveMat * mat = new TPZKarhunenLoeveMat ( id,Lx,Ly,Lz,dim,type );
        cmesh->SetDefaultOrder ( porder );
        cmesh->SetDimModel ( dim );
        cmesh->InsertMaterialObject ( mat );
        cmesh->SetAllCreateFunctionsContinuous();
        cmesh->AutoBuild();
        return cmesh;
}

TPZGeoMesh * TriGMesh ( int ref )
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

        string meshref = "gmeshtri3.vtk";
        std::ofstream files ( meshref );
        TPZVTKGeoMesh::PrintGMeshVTK ( gmesh,files,true );
        return gmesh;
}

void saveVector(const std::vector<double>& vec, const std::string& filename) {
    // Open the file in binary mode
    std::ofstream outFile(filename, std::ios::binary);

    // Check if the file is open
    if (!outFile.is_open()) {
        std::cerr << "Error: Could not open file for writing." << std::endl;
        return;
    }

    // Write the size of the vector first
    size_t size = vec.size();
    outFile.write(reinterpret_cast<const char*>(&size), sizeof(size));

    // Write the vector data
    outFile.write(reinterpret_cast<const char*>(vec.data()), size * sizeof(double));

    // Close the file
    outFile.close();
}



void readVector( std::vector<double>& vec,const std::string& filename)
{
    // Open the file in binary mode
    std::ifstream inFile(filename, std::ios::binary);

    // Check if the file is open
    if (!inFile.is_open()) {
        std::cerr << "Error: Could not open file for reading." << std::endl;
    }

    // Read the size of the vector first
    size_t size;
    inFile.read(reinterpret_cast<char*>(&size), sizeof(size));

    // Resize the vector to hold the dat

    vec.resize(size);
    // Read the vector data
    inFile.read(reinterpret_cast<char*>(vec.data()), size * sizeof(double));

    // Close the file
    inFile.close();

}
