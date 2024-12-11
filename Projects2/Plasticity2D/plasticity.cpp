#include <iostream>
#include <fstream>
#include <thread>
#include <vector>
#include "TPZFileStream.h"
#include <TPZBFileStream.h>
#include "SlopeAnalysis.h"
#include <fstream>
#include <iostream>
#include <fstream>
#include <thread>
#include <vector>
#include <mutex>
#include <iostream>
#include <fstream>
#include <unistd.h>
#include <sys/wait.h>
//std::mutex mtx; // Mutex para proteger a escrita nos arquivos
void RunParallelSlopeAnalysis ( int imc_start, int imc_end, int num_processes, SlopeAnalysis* slopeanalysisf, SlopeAnalysis* slopeanalysish );
void RunParallelSlopeAnalysis ( int imc_start, int imc_end, int num_processes, SlopeAnalysis* slopeanalysis);
void SolveSlope ( int imc_start, int imc_end, SlopeAnalysis* slopeanalysisf, SlopeAnalysis* slopeanalysish );
void SolveSlope ( int imc_start, int imc_end, SlopeAnalysis* slopeanalysis );
void SolveSlope ( int Startfrom );
TPZGeoMesh * TriGMesh ( int ref );
TPZCompMesh* CreateCompMeshKL ( TPZGeoMesh * gmesh,int porder,REAL Lx, REAL Ly, REAL Lz, int id,int type );
TPZCompMesh * CreateCMeshElastoplastic ( TPZGeoMesh *gmesh, int pOrder, REAL coes,REAL atrito );
TPZCompMesh CreateCompMeshKL2 ( TPZGeoMesh  gmesh,int porder,REAL Lx, REAL Ly, REAL Lz, int id,int type );

void saveVector ( const std::vector<double>& vec, const std::string& filename, bool append = true ) ;

void readVector ( std::vector<double>& vec,const std::string& filename );

std::vector<int> GetIndex ( std::vector<double> vetor );

void SolveSlopeIS ( int Startfrom );

void ManageStartFrom(int Startfrom);

REAL func(REAL theta,REAL cov,REAL mean);
void PrintBases(TPZVec<TPZFMatrix<REAL>> allsamples, std::vector<int> indexes, TPZVec<REAL>covvec, TPZVec<REAL>meanvec);
//void PrintBasesMathematica(TPZVec<TPZFMatrix<REAL>> allsamples, std::vector<int> indexes, TPZVec<REAL>covvec, TPZVec<REAL>meanvec);
void PrintBasesMathematica(const TPZVec<TPZFMatrix<REAL>>& allsamples, const std::vector<int>& indexes, const TPZVec<REAL>& covvec, const TPZVec<REAL>& meanvec);
int main()
{

        int Startfrom =2;
        ManageStartFrom ( Startfrom );


        return 0;
}

REAL func(REAL theta,REAL cov,REAL mean)
{
        REAL xi = sqrt ( log ( 1. + cov*cov ) );
        REAL lambda = log ( mean ) - 0.5*xi * xi;
        //return exp ( lambda + xi *theta );
        return  exp(xi *theta) ;
}



void ManageStartFrom(int Startfrom)
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
        int solvertype=0;
        int numthreads=10;
        int ref0slope=3;
        int porderslope=1;
        REAL gammaagua=0.;
        REAL gammasolo=20.;
        REAL coes=10.;
        REAL atrito=30.*M_PI/180.;

        SlopeAnalysis  * slopeanalysis =  new SlopeAnalysis ( gammaagua,gammasolo,coes,atrito,ref0slope,porderslope,numthreads,solvertype );

       // bool issrm=true;
       // slopeanalysis->SolveDeterministic(issrm);
       // std::string saidavtk2 = "postdeter.vtk";
       // slopeanalysis->PostPlasticity ( saidavtk2 );

     //   return;
        if ( Startfrom ==0 ) {
                randonanalysis->SetNEigenpairs ( 1500 );
                //randonanalysis->Assemble();
                randonanalysis->Solve();

                //save sqrt(lambda)*phi
                TPZBFileStream save;
                save.OpenWrite ( "Config0.bin" );
                randonanalysis->Write ( save,randonanalysis->ClassId() );
                //randonanalysis->LoadSolution();
                randonanalysis->DefineGraphMesh ( 2,scalarnames,vecnames,"filename2Assemble.vtk" );
                randonanalysis->PostProcess ( 0 );
        }


        if ( Startfrom >0 ) {
                //read sqrt(lambda)*phi
                TPZBFileStream read;
                read.OpenRead ( "Config0.bin" );
                randonanalysis->Read ( read,0 );

                //seting fied data
                TPZVec<REAL> meanvec ( 2 );
                meanvec[0]=coes;
                meanvec[1]=atrito;
                TPZVec<REAL> covvec ( 2 );
                covvec[0]=0.3;
                covvec[1]=0.2;
                int samples=1000;
                slopeanalysis->SetFieldsData ( cmesh,randonanalysis->GetSolutionValVec(), meanvec,covvec,  samples );


                if ( Startfrom==1 ) {


                        slopeanalysis->ManageFieldCretion();
                        TPZBFileStream save;
                        save.OpenWrite ( "Config1.bin" );
                        slopeanalysis->Write ( save,slopeanalysis->ClassId() );


                } else { //Startfrom>1 solve monte carlo
                        TPZBFileStream read;
                        read.OpenRead ( "Config1.bin" );
                        slopeanalysis->Read ( read,0 );

                       // std::ofstream posprocfs ( "posprocfs.txt" );


                        int n=10;
                        REAL p0=0.5;
                      //  slopeanalysis->SubSet(n,  p0);

                        for(int imc=74;imc<75;imc++)
                        {

                                cout << "simulacao de monte carlo numero = "<<imc <<endl;
                                SlopeAnalysis* analysis = new SlopeAnalysis ( *slopeanalysis );
                                //REAL fs = analysis->SolveSingleField ( imc );

                                //analysis->MetropolisHastings(1,74);

                               // posprocfs << imc << " "<<fs << endl;
                                //std::string saidavtk2 = "postvtkarclength/saidavtk" + std::to_string ( imc ) + ".vtk";
                                //analysis->PostPlasticity ( saidavtk2 );
                                //std::string saidafs2 = "postarclength/fs" + std::to_string ( imc ) + ".dat";
                                //std::ofstream out2 ( saidafs2 );
                                //out2 << fs << std::endl;

                               // delete analysis;
                        }


                }


        }


        cout << "EXIT SUCESS"<<endl;
}



/*
void ManageStartFrom(int Startfrom)
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
        int solvertype=0;
        int numthreads=10;
        int ref0slope=3;
        int porderslope=1;
        REAL gammaagua=0.;
        REAL gammasolo=20.;
        REAL coes=10.;
        REAL atrito=30.*M_PI/180.;

        SlopeAnalysis  * slopeanalysisf =  new SlopeAnalysis ( gammaagua,gammasolo,coes,atrito,ref0slope,porderslope,numthreads,solvertype );

       // slopeanalysisf->SolveDeterministic();
        if ( Startfrom ==0 ) {
                randonanalysis->SetNEigenpairs ( 1500 );
                //randonanalysis->Assemble();
                randonanalysis->Solve();

                //save sqrt(lambda)*phi
                TPZBFileStream save;
                save.OpenWrite ( "Config1-0-fs13.bin" );
                randonanalysis->Write ( save,randonanalysis->ClassId() );
                //randonanalysis->LoadSolution();
                randonanalysis->DefineGraphMesh ( 2,scalarnames,vecnames,"filename2Assemble.vtk" );
                randonanalysis->PostProcess ( 0 );
        }


        if ( Startfrom >0 ) {
                //read sqrt(lambda)*phi
                TPZBFileStream read;
                read.OpenRead ( "Config1-0-fs13.bin" );
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


                if ( Startfrom==1 ) {

                        TPZFMatrix<REAL>samples1 = slopeanalysisf->CreateNormalStandardSamples();
                        TPZFMatrix<REAL>samples2 = slopeanalysisf->CreateNormalStandardSamples();
                        TPZVec<TPZFMatrix<REAL>> samples ( 2 );
                        samples[0]=samples1;
                        samples[1]=samples2;
                        samples1.Print("samples");
                        samples2.Print("samples");
                        slopeanalysisf->SetFieldsSamples ( samples );
                        slopeanalysisf->ManageFieldCretion();
                        TPZBFileStream save,saveh;
                        save.OpenWrite ( "Config2-0-fs13.bin" );
                        saveh.OpenWrite ( "Configh-2-0-fs13.bin" );
                        slopeanalysisf->Write ( save,slopeanalysisf->ClassId() );


                } else { //Startfrom>1 solve monte carlo
                        TPZBFileStream read,readh;
                        read.OpenRead ( "Config2-0-fs13.bin" );
                        readh.OpenRead ( "Configh-2-0-fs13.bin" );
                        slopeanalysisf->Read ( read,0 );

//                         int imc_start = 0;
//                         int imc_end = 20;
//
//                         TPZVec<TPZFMatrix<REAL>> allsamples =slopeanalysisf->GetFieldsSamples();
//                        std::vector<int> indexes={28,28,51,110,182,187,201,214,266,307,322,324,348,388,394,402,416,417,485,503,555,598,630,650,660,683,
//                                698,720,746,765,845,870,874,923,928,940,949,955,969,985,992,996
//                 };
//                         std::vector<int> indexes={52,74,81,142,177,204,243,259,261,277,281,308,330,356,384,446,458,504,553,590,592,607,625,637,676,715,739,773,
// 778,783,785,789,793,894,900,905,926,936,959};
                        //PrintBases(allsamples, indexes, covvec, meanvec);
//                       // PrintBasesMathematica(allsamples, indexes, covvec, meanvec);
//                        return;
                        //std::vector<int> critical={4,1446,889,25,646,1099,1491,558,1042,413,1117};
//std::vector<std::vector<int>> critical={{4,1446,889,25,646,1099,1491,558,1042,413},{4,2,1327,135,237,25,22,572,751,1447}};//min
//std::vector<std::vector<int>> critical={{1334,1174,133,254,1207,1005,217,1395,1309,105},{1280,295,138,1245,430,14,1312,1071,1473,1206}};//max
  //                       std::vector<std::vector<int>> critical={
 //                       {1334,1174,133,254,1207,1005,217,1395,1309,105,4,1446,889,25,646,1099,1491,558,1042,413},
//                         {1280,295,138,1245,430,14,1312,1071,1473,1206,4,2,1327,135,237,25,22,572,751,1447}};//max+min
//                         std::vector<std::vector<int>> critical={
//                         {4,22,849,576,461,889,1143,242,1124,13},
//                          {1379,1372,937,1171,510,471,866,274,349,149}};//min sigma
//                          std::vector<std::vector<int>> critical={
//                         {896,503,96,1334,1281,69,1174,1311,690,1156},
//                          {1078,1198,1177,1203,14,664,854,1161,518,975}};//max sigma
 //                       slopeanalysisf->ManageFieldCretion(critical);
                       // std::vector<std::vector<double>> failedfields1,failedfields2;
                       // std::string posprocfs = "posprocfs.txt";
                        std::ofstream posprocfs ( "posprocfs.txt" );

                        // std::vector<int> critical={30,50,35,70,9,21,76,97,11,19,49,4,99,73,58,54,29,67};
                         std::vector<int> critical={4, 9, 11, 19, 21, 29, 30, 35, 49, 54, 58, 70, 73, 76, 97, 99, 111, 126, 127, 128, 132, 135, 136, 137, 140, 167, 176, 179, 185, 186, 187, 193, 194, 198, 206, 207, 215, 218, 219, 224, 233, 242, 247, 252, 253, 266, 269, 271, 275, 284, 287, 288, 294, 297, 299, 302, 307, 319, 323, 326, 336, 337, 345, 346, 350, 357, 369, 372, 382, 391, 401, 402, 413, 425, 431, 438, 439, 452, 462, 465, 470, 484, 487, 488, 491, 507, 520, 539, 541, 544, 551, 553, 556, 558, 564, 566, 570, 574, 585, 593, 596, 602, 604, 609, 611, 612, 614, 622, 625, 631, 633, 644, 646, 652, 661, 662, 664, 665, 678, 684, 685, 688, 696, 713, 714, 716, 718, 720, 722, 733, 742, 744, 746, 753, 754, 758, 765, 767, 775, 777, 783, 796, 797, 802, 805, 816, 818, 820, 824, 827, 838, 839, 840, 855, 858, 860, 865, 871, 878, 884, 886, 888, 897, 903, 904, 909, 910, 914, 915, 919, 920, 922, 923, 927, 928, 943, 952, 954, 958, 963, 969, 973, 989, 993
};

                       // for(int imc=0;imc<critical.size();imc++)
                        for(int imc=0;imc<1000;imc++)
                        {
                                //int imcc=critical[imc];
                                cout << "simulacao de monte carlo numero = "<<imc <<endl;
                                SlopeAnalysis* analysis = new SlopeAnalysis ( *slopeanalysisf );
                               /// analysis->IntegrateFieldOverARegionB (0, imc );
                                REAL fs = analysis->SolveSingleField ( imc );
                                posprocfs << imc << " "<<fs << endl;
                                //slopeanalysisf->IntegrateFieldOverARegion(0.02,imc);
                                std::string saidavtk2 = "postvtkselected/saidavtk" + std::to_string ( imc ) + ".vtk";
                                analysis->PostPlasticity ( saidavtk2 );
                                std::string saidafs2 = "postselected/fs" + std::to_string ( imc ) + ".dat";
                                std::ofstream out2 ( saidafs2 );
                                out2 << fs << std::endl;

                                delete analysis;
                        }


                }


        }


        cout << "total"<<endl;
}*/
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
        int solvertype=0;
        int numthreads=10;
        int ref0slope=3;
        int porderslope=1;
        REAL gammaagua=0.;
        REAL gammasolo=20.;
        REAL coes=10.;
        REAL atrito=30.*M_PI/180.;
        REAL fsfail=1.2;
        REAL coesh=coes/fsfail;
        REAL atritoh= atan ( tan ( atrito ) /fsfail );

        SlopeAnalysis  * slopeanalysisf =  new SlopeAnalysis ( gammaagua,gammasolo,coes,atrito,ref0slope,porderslope,numthreads,solvertype );

       // slopeanalysisf->SolveDeterministic();

//return;
        SlopeAnalysis  * slopeanalysish =  new SlopeAnalysis ( gammaagua,gammasolo,coesh,atritoh,ref0slope,porderslope,numthreads,solvertype );

//         slopeanalysish->SolveDeterministic();
//         return;
        if ( Startfrom ==0 ) {
                randonanalysis->SetNEigenpairs ( 1653 );
                //randonanalysis->Assemble();
                randonanalysis->Solve();

                //save sqrt(lambda)*phi
                TPZBFileStream save;
                save.OpenWrite ( "Config1-0-fs13.bin" );
                randonanalysis->Write ( save,randonanalysis->ClassId() );
                //randonanalysis->LoadSolution();
                randonanalysis->DefineGraphMesh ( 2,scalarnames,vecnames,"filename2Assemble.vtk" );
                randonanalysis->PostProcess ( 0 );
        }


        if ( Startfrom >0 ) {
                //read sqrt(lambda)*phi
                TPZBFileStream read;
                read.OpenRead ( "Config1-0-fs13.bin" );
                randonanalysis->Read ( read,0 );

                //seting fied data
                TPZVec<REAL> meanvec ( 2 );
                meanvec[0]=coes;
                meanvec[1]=atrito;
                TPZVec<REAL> covvec ( 2 );
                covvec[0]=0.3;
                covvec[1]=0.2;
                int samples=10000;
                slopeanalysisf->SetFieldsData ( cmesh,randonanalysis->GetSolutionValVec(), meanvec,covvec,  samples );

                TPZVec<REAL> meanvech ( 2 );
                meanvech[0]=coesh;
                meanvech[1]=atritoh;
                TPZVec<REAL> covvech ( 2 );
                covvech[0]=0.3;
                covvech[1]=0.2;
                slopeanalysish->SetFieldsData ( cmesh,randonanalysis->GetSolutionValVec(), meanvech,covvech,  samples );

                //crate random fields
                if ( Startfrom==1 ) {

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
                        save.OpenWrite ( "Config2-0-fs13.bin" );
                        saveh.OpenWrite ( "Configh-2-0-fs13.bin" );
                        slopeanalysisf->Write ( save,slopeanalysisf->ClassId() );
                        slopeanalysish->Write ( saveh,slopeanalysish->ClassId()+1 );

                } else { //Startfrom>1 solve monte carlo
                        TPZBFileStream read,readh;
                        read.OpenRead ( "Config2-0-fs13.bin" );
                        readh.OpenRead ( "Configh-2-0-fs13.bin" );
                        slopeanalysisf->Read ( read,0 );
                        slopeanalysish->Read ( readh,0 );

                        int imc_start = 100;
                        int imc_end = 200;
                        int num_processes =2;  // Dividir para 4 processos

                       // SolveSlope (  imc_start,  imc_end,slopeanalysisf );
                        // Chama a função para executar a análise paralela entre imc_start e imc_end
                       //RunParallelSlopeAnalysis ( imc_start, imc_end, num_processes, slopeanalysisf, slopeanalysish );
                       //SolveSlope(imc_start, imc_end, slopeanalysisf, slopeanalysish);
                       // RunParallelSlopeAnalysis ( imc_start, imc_end, num_processes, slopeanalysisf );




                }


        }


        cout << "total"<<endl;

}

void RunParallelSlopeAnalysis ( int imc_start, int imc_end, int num_processes, SlopeAnalysis* analysis )
{
        int total_imc = imc_end - imc_start;
        int imc_per_process = total_imc / num_processes;

        for ( int i = 0; i < num_processes; ++i ) {
                pid_t pid = fork();  // Cria um novo processo

                if ( pid == 0 ) { // Processo filho
                        int local_imc_start = imc_start + i * imc_per_process;
                        int local_imc_end = local_imc_start + imc_per_process;

                        if ( i == num_processes - 1 ) {
                                local_imc_end = imc_end;  // O último processo vai até o fim
                        }

                        // Chama a função para resolver o intervalo de imc
                        SolveSlope ( local_imc_start, local_imc_end, analysis );

                        _exit ( 0 ); // Termina o processo filho quando completar o intervalo
                } else if ( pid > 0 ) {
                        // Processo pai continua e cria outro filho
                        continue;
                } else {
                        std::cerr << "Erro ao criar processo!" << std::endl;
                        return;
                }
        }

        // Processo pai aguarda todos os filhos finalizarem
        for ( int i = 0; i < num_processes; ++i ) {
                int status;
                wait ( &status ); // Espera pelo término de cada processo filho
        }
}

void RunParallelSlopeAnalysis ( int imc_start, int imc_end, int num_processes, SlopeAnalysis* slopeanalysisf, SlopeAnalysis* slopeanalysish )
{
        int total_imc = imc_end - imc_start;
        int imc_per_process = total_imc / num_processes;

        for ( int i = 0; i < num_processes; ++i ) {
                pid_t pid = fork();  // Cria um novo processo

                if ( pid == 0 ) { // Processo filho
                        int local_imc_start = imc_start + i * imc_per_process;
                        int local_imc_end = local_imc_start + imc_per_process;

                        if ( i == num_processes - 1 ) {
                                local_imc_end = imc_end;  // O último processo vai até o fim
                        }

                        // Chama a função para resolver o intervalo de imc
                        SolveSlope ( local_imc_start, local_imc_end, slopeanalysisf, slopeanalysish );

                        _exit ( 0 ); // Termina o processo filho quando completar o intervalo
                } else if ( pid > 0 ) {
                        // Processo pai continua e cria outro filho
                        continue;
                } else {
                        std::cerr << "Erro ao criar processo!" << std::endl;
                        return;
                }
        }

        // Processo pai aguarda todos os filhos finalizarem
        for ( int i = 0; i < num_processes; ++i ) {
                int status;
                wait ( &status ); // Espera pelo término de cada processo filho
        }
}

void SolveSlope ( int imc_start, int imc_end, SlopeAnalysis* slopeanalysis )
{
         for ( int imc = imc_start; imc < imc_end; ++imc ) {
                std::cout << "imc = " << imc << " (PID: " << getpid() << ")" << std::endl;

                SlopeAnalysis* analysis = new SlopeAnalysis ( *slopeanalysis );

                REAL fs = analysis->SolveSingleField ( imc );

                analysis->IntegrateFieldOverARegion(0.02,imc);

                std::string saidavtk2 = "postvtkx/saidavtk" + std::to_string ( imc ) + "h.vtk";

                analysis->PostPlasticity ( saidavtk2 );

                std::string saidafs2 = "postx/fs" + std::to_string ( imc ) + "h.dat";

                std::ofstream out2 ( saidafs2 );

                out2 << fs << std::endl;

                delete analysis;
        }
}
void SolveSlope ( int imc_start, int imc_end, SlopeAnalysis* slopeanalysisf, SlopeAnalysis* slopeanalysish )
{
        for ( int imc = imc_start; imc < imc_end; ++imc ) {
                std::cout << "imc = " << imc << " (PID: " << getpid() << ")" << std::endl;

                SlopeAnalysis* slopeanalysisf1 = new SlopeAnalysis ( *slopeanalysisf );
                SlopeAnalysis* slopeanalysish1 = new SlopeAnalysis ( *slopeanalysish );

                REAL fsf = slopeanalysisf1->SolveSingleField ( imc );
                REAL fsh = slopeanalysish1->SolveSingleField ( imc );

                //std::string saidavtk = "postvtk-fs13/saidavtk" + std::to_string ( imc ) + ".vtk";
                //std::string saidavtk2 = "postvtk-fs13/saidavtk" + std::to_string ( imc ) + "h.vtk";


                //slopeanalysisf1->PostPlasticity ( saidavtk );

                //slopeanalysish1->PostPlasticity ( saidavtk2 );

                std::string saidafs = "postx/fs" + std::to_string ( imc ) + ".dat";
                std::string saidafs2 = "postx/fs" + std::to_string ( imc ) + "h.dat";
                std::ofstream out ( saidafs );
                std::ofstream out2 ( saidafs2 );

                out << fsf << std::endl;
                out2 << fsh << std::endl;

                delete slopeanalysisf1;
                delete slopeanalysish1;
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

TPZCompMesh CreateCompMeshKL2 ( TPZGeoMesh  gmesh,int porder,REAL Lx, REAL Ly, REAL Lz, int id,int type )
{

        int dim = gmesh.Dimension();
        TPZCompMesh  cmesh =  TPZCompMesh ( &gmesh );
        TPZKarhunenLoeveMat * mat = new TPZKarhunenLoeveMat ( id,Lx,Ly,Lz,dim,type );
        cmesh.SetDefaultOrder ( porder );
        cmesh.SetDimModel ( dim );
        cmesh.InsertMaterialObject ( mat );
        cmesh.SetAllCreateFunctionsContinuous();
        cmesh.AutoBuild();
        return cmesh;
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

void saveVector ( const std::vector<double>& vec, const std::string& filename, bool append )
{
        // Open the file in binary mode, append if needed
        std::ios_base::openmode mode = std::ios::binary;
        if ( append ) {
                mode |= std::ios::app;  // Use append mode
        }

        std::ofstream outFile ( filename, mode );

        // Check if the file is open
        if ( !outFile.is_open() ) {
                std::cerr << "Error: Could not open file for writing/appending." << std::endl;
                return;
        }

        // Write the size of the vector first
        size_t size = vec.size();
        outFile.write ( reinterpret_cast<const char*> ( &size ), sizeof ( size ) );

        // Write the vector data
        outFile.write ( reinterpret_cast<const char*> ( vec.data() ), size * sizeof ( double ) );

        // Close the file
        outFile.close();
}



void readVector ( std::vector<double>& vec, const std::string& filename )
{
        // Open the file in binary mode
        std::ifstream inFile ( filename, std::ios::binary );

        // Check if the file is open
        if ( !inFile.is_open() ) {
                std::cerr << "Error: Could not open file for reading." << std::endl;
                return;
        }

        // Read the size of the vector first
        size_t size;
        inFile.read ( reinterpret_cast<char*> ( &size ), sizeof ( size ) );

        // Resize the vector to hold the data
        vec.resize ( size );

        // Read the vector data
        inFile.read ( reinterpret_cast<char*> ( vec.data() ), size * sizeof ( double ) );

        // Close the file
        inFile.close();
}
void PrintBases(TPZVec<TPZFMatrix<REAL>> allsamples, std::vector<int> indexes, TPZVec<REAL>covvec, TPZVec<REAL>meanvec)
{
        std::ofstream datacoes ( "basesc.dat" );
        std::ofstream dataatrito ( "basesa.dat" );

        for(int irow=0;irow<allsamples[0].Rows();irow++)
        {
                datacoes<< " base "<< irow ;
                dataatrito<< " base "<< irow;
                for(int ifailed=0;ifailed<indexes.size();ifailed++)
                {
                        datacoes<<  "  "<< func(allsamples[0](irow,indexes[ifailed]),covvec[0],meanvec[0]);
                        dataatrito<< "  "<< func(allsamples[1](irow,indexes[ifailed]),covvec[1],meanvec[1]);
                }
//                 for(int ifailed=0;ifailed<allsamples[0].Cols();ifailed++){
//                         datacoes<<  "  "<< func(allsamples[0](irow,ifailed),covvec[0],meanvec[0]);
//                         dataatrito<< "  "<< func(allsamples[1](irow,ifailed),covvec[1],meanvec[1]);
//                 }
                datacoes<< endl;
                dataatrito<< endl;

        }

}


// Função para exportar as bases no formato Mathematica
void PrintBasesMathematica(const TPZVec<TPZFMatrix<REAL>>& allsamples, const std::vector<int>& indexes, const TPZVec<REAL>& covvec, const TPZVec<REAL>& meanvec) {
    std::ofstream arquivo("dados_bases.nb");
    if (!arquivo) {
        std::cerr << "Erro ao abrir o arquivo para escrita!" << std::endl;
        return;
    }

    arquivo << "(*Dados no formato especificado*)\n";
    arquivo << "dados = {";

    for (int irow = 0; irow < allsamples[0].Rows(); irow++) {
        arquivo << "{\"base " << irow << "\", ";

        for (size_t i = 0; i < indexes.size(); i++) {
            double valorCoesao = func(allsamples[0](irow, indexes[i]), covvec[0], meanvec[0]);
            double valorAtrito = func(allsamples[1](irow, indexes[i]), covvec[1], meanvec[1]);

            arquivo << std::fixed << std::setprecision(6) << valorCoesao;

            if (i < indexes.size() - 1) {
                arquivo << ", ";
            }
        }

        arquivo << "}";
        if (irow < allsamples[0].Rows() - 1) {
            arquivo << ", ";
        }
    }

    arquivo << "};\n";
    std::cout << "Dados exportados com sucesso para 'dados_bases.nb'" << std::endl;
}
/*
void PrintBasesMathematica(TPZVec<TPZFMatrix<REAL>> allsamples, std::vector<int> indexes, TPZVec<REAL> covvec, TPZVec<REAL> meanvec)
{
    // Abre os arquivos de saída para coesão e atrito
    std::ofstream datacoes("basesdacoesao.nb");
    std::ofstream dataatrito("basesdoatrito.nb");

    // Início da estrutura em formato de lista para Mathematica
    datacoes << "{";
    dataatrito << "{";

    // Loop sobre as linhas das amostras
    for (int irow = 0; irow < allsamples[0].Rows(); irow++)
    {
        datacoes << "{";
        dataatrito << "{";

        // Loop sobre os índices falhos (indexes)
        for (int ifailed = 0; ifailed < indexes.size(); ifailed++)
        {
            // Calcula os valores normalizados usando a função `func`
            REAL coesaoValue = func(allsamples[0](irow, indexes[ifailed]), covvec[0], meanvec[0]);
            REAL atritoValue = func(allsamples[1](irow, indexes[ifailed]), covvec[1], meanvec[1]);

            // Escreve os valores no arquivo com vírgula até o último elemento
            datacoes << coesaoValue;
            dataatrito << atritoValue;

            if (ifailed < indexes.size() - 1)
            {
                datacoes << ", ";
                dataatrito << ", ";
            }
        }

        datacoes << "}";
        dataatrito << "}";

        // Adiciona uma vírgula entre as linhas, exceto na última
        if (irow < allsamples[0].Rows() - 1)
        {
            datacoes << ", ";
            dataatrito << ", ";
        }
    }

    // Finaliza a estrutura de lista em formato Mathematica
    datacoes << "};";
    dataatrito << "};";

    // Fecha os arquivos de saída
    datacoes.close();
    dataatrito.close();
}*/

// void PrintBasesMathematica(TPZVec<TPZFMatrix<REAL>> allsamples, std::vector<int> indexes, TPZVec<REAL>covvec, TPZVec<REAL>meanvec)
// {
//         std::ofstream datacoes ( "basesdacoesao.nb" );
//         std::ofstream dataatrito ( "basesdoatrito.nb" );
//
//         for(int irow=0;irow<allsamples[0].Rows();irow++)
//         {
//
//                 datacoes << "{";
//                 dataatrito << "{";
//                 for(int ifailed=0;ifailed<indexes.size();ifailed++)
//                 {
//                         if(ifailed<indexes.size()-1)
//                         {
//                                 datacoes<<  func(allsamples[0](irow,indexes[ifailed]),covvec[0],meanvec[0])<< ",";
//                                 dataatrito<< func(allsamples[1](irow,indexes[ifailed]),covvec[1],meanvec[1]) << ",";
//                         }else{
//                                                               datacoes<<  func(allsamples[0](irow,indexes[ifailed]),covvec[0],meanvec[0])<< " ";
//                                 dataatrito<< func(allsamples[1](irow,indexes[ifailed]),covvec[1],meanvec[1])<< " ";
//                         }
//                 }
//                 if(irow<allsamples[0].Rows()-1)
//                 {
//                         datacoes<< "},";
//                         dataatrito<< "},";
//                 }else{
//                         datacoes<< "};";
//                         dataatrito<< "};";
//                 }
//
//         }
//
// }
