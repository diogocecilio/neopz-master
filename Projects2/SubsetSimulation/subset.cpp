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

void ManageStartFrom(int Startfrom);

void SubSet(int n, REAL p0,SlopeAnalysis* slopeanalysis );

TPZStack<TPZVec<TPZFMatrix<REAL>>> MetropolisHastings(int nnewsamples,TPZVec<TPZFMatrix<REAL>> seedinit,SlopeAnalysis* analysis, std::vector<double> &fsvec, REAL b);

std::vector<std::pair<int, double>> CrudeMonteCarlo(int a,int b, SlopeAnalysis* slopeanalysis);

void Write(TPZStream &buf, int withclassid) ;

void Read(TPZStream &buf, void *context);

vector<size_t> sort_indexes(const vector<double> &v);

std::pair<std::vector<double>,std::vector<int>> Sort(std::vector<double> vec,int chop);

std::pair<std::vector<double>,std::vector<int>> Sort(std::vector<double> vec,int chop)
{
        std::pair<std::vector<double>,std::vector<int>> oudata;
        int counter=0;
        for (auto i: sort_indexes(vec))
        {
                if(counter>=chop)
                {
                        oudata.first.push_back(vec[i]);
                        oudata.second.push_back(i);
                }
                counter++;
        }
        return oudata;
}

int main()
{


        int Startfrom =2;
        ManageStartFrom ( Startfrom );

/*
        std::vector<double> fsvec ={3.,8.,1.};
        std::pair<std::vector<double>,std::vector<int>> out;

        out=Sort(fsvec,0);



        for (size_t i = 0; i < out.first.size(); ++i)
        {
                cout << "First (valor): " << out.first[i] << ", Second (índice): " << out.second[i] << endl;
        }
        cout << "HELLO WORLD" <<endl;

TPZVec<TPZFMatrix<REAL>> vecmat(2),vecmat2(2),vecmat3(2);
TPZFMatrix<REAL> a,b,c,d,e,f;
a.AutoFill(3,3,1);
b.AutoFill(2,3,1);
c.AutoFill(3,4,1);
d.AutoFill(4,5,1);
e.AutoFill(1,3,1);
f.AutoFill(2,2,1);

vecmat[0]=a;
vecmat[1]=b;
vecmat2[0]=c;
vecmat2[1]=d;
vecmat3[0]=e;
vecmat3[1]=f;


TPZStack<TPZVec<TPZFMatrix<REAL>>> stackvecmat;


stackvecmat.Push(vecmat);

stackvecmat.Push(vecmat2);

stackvecmat.Push(vecmat3);

                int szstack = stackvecmat.size();

                cout << "szstack = "<< szstack << " out.first.size() ="<<out.first.size() << endl;

                for (int i = 0; i < stackvecmat.size(); i++)
                {
                        cout << "First : " << out.first[i] << ", Second  : " << out.second[i] << endl;
                        //armazena as novas amostras geradas em ordem decrescente do fator de seguranca
                        stackvecmat[out.second[i]][0].Print(cout);
                }*/
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
        int porderslope=2;
        REAL gammaagua=0.;
        REAL gammasolo=20.;
        REAL coes=15.;
        REAL atrito=30.*M_PI/180.;

        SlopeAnalysis  * slopeanalysis =  new SlopeAnalysis ( gammaagua,gammasolo,coes,atrito,ref0slope,porderslope,numthreads,solvertype );

//        bool issrm=false;
//        slopeanalysis->SolveDeterministic(issrm);
//        std::string saidavtk2 = "postdeter.vtk";
//        slopeanalysis->PostPlasticity ( saidavtk2 );
//
//        return;

        if ( Startfrom ==0 ) {
                randonanalysis->SetNEigenpairs ( 1500 );
                //randonanalysis->Assemble();
                randonanalysis->Solve();

                //save sqrt(lambda)*phi
                TPZBFileStream save;
                save.OpenWrite ( "Configlowpf.bin" );
                randonanalysis->Write ( save,randonanalysis->ClassId() );
                //randonanalysis->LoadSolution();
                randonanalysis->DefineGraphMesh ( 2,scalarnames,vecnames,"filename2Assemble.vtk" );
                randonanalysis->PostProcess ( 0 );
        }


        if ( Startfrom >0 ) {
                cout << "AQ"<<endl;
                //read sqrt(lambda)*phi
                TPZBFileStream read;
                read.OpenRead ( "Configlowpf.bin" );
                randonanalysis->Read ( read,0 );
                cout << "A2"<<endl;
                //seting fied data
                TPZVec<REAL> meanvec ( 2 );
                meanvec[0]=coes;
                meanvec[1]=atrito;
                TPZVec<REAL> covvec ( 2 );
                covvec[0]=0.3;
                covvec[1]=0.2;
                int samples=10000;
                slopeanalysis->SetFieldsData ( cmesh,randonanalysis->GetSolutionValVec(), meanvec,covvec,  samples );


                if ( Startfrom==1 ) {


                        slopeanalysis->ManageFieldCretion();
                        TPZBFileStream save;
                        save.OpenWrite ( "Configlowpf1.bin" );
                        slopeanalysis->Write ( save,slopeanalysis->ClassId() );


                } else { //Startfrom>1 solve monte carlo


                        TPZBFileStream read;
                        read.OpenRead ( "Configlowpf1.bin" );
                        slopeanalysis->Read ( read,0 );
                        //CrudeMonteCarlo(1000,10000,slopeanalysis);
                        cout << "A3"<<endl;
                        int n=10;
                        REAL p0=0.1;
                        SubSet(n, p0,slopeanalysis );
                }


        }


        cout << "EXIT SUCESS"<<endl;
}


vector<size_t> sort_indexes(const vector<double> &v) {
  vector<size_t> idx(v.size());
  iota(idx.begin(), idx.end(), 0);
  stable_sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return v[i1] > v[i2];});

  return idx;
}

std::ofstream posprocfs ( "posprocfs.txt" );
std::ofstream posprocfs2 ( "posprocfs2.txt" );
std::ofstream posprocfs3 ( "posprocfs3.txt" );


void SubSet(int n, REAL p0,SlopeAnalysis* slopeanalysis )
{


        //simulacao inicial de monte carlo com n amostras

//                 std::vector<std::pair<int, double>> fsdata ={
//                         {0,1.65563},{14,1.55228},{91,1.44918},{29,1.43568},{98,1.41855},
//                         {31,1.39284},{18,1.37799},{55,1.36194},{94,1.36182},{84,1.35785},{70,1.33855},{100,1.33394},{44,1.33035},{16,1.32551},{51,1.32047},{23,1.31908},{43,1.31907},{87,1.31888},{46,1.31668},{54,1.31108},{38,1.30393},{89,1.30203},{19,1.2896},{2,1.28696},{49,1.286},{82,1.2821},{58,1.27589},{52,1.27565},{71,1.27106},{35,1.2644},{56,1.26138},{25,1.25786},{1,1.25218},{92,1.2511},{3,1.24322},{50,1.24257},{24,1.23259},{95,1.23068},{41,1.2301},{33,1.23},{99,1.22459},{34,1.22171},{21,1.20724},{39,1.20113},{73,1.20106},{97,1.19908},{10,1.19815},{7,1.19322},{60,1.18987},{88,1.18959},{5,1.1831},{75,1.18262},{22,1.18223},{45,1.17178},{48,1.16989},{76,1.16828},{77,1.16828},{27,1.16365},{67,1.16152},{4,1.16039},{57,1.15895},{64,1.15819},{13,1.14604},{62,1.14531},{28,1.14335},{9,1.1431},{26,1.14245},{78,1.13907},{40,1.13674},{32,1.13039},{59,1.12995},{66,1.12423},{20,1.12217},{42,1.11966},{68,1.1133},{37,1.11108},{15,1.10585},{47,1.10535},{86,1.0959},{83,1.09571},{8,1.09043},{63,1.08922},{12,1.07339},{6,1.07149},{36,1.06869},{61,1.0668},{11,1.06353},{69,1.06345},{65,1.04853},{80,1.04738},{96,1.03642},{53,1.03314},{72,1.01948},{17,1.0192},{81,1.01457},{30,1.00811},{90,0.98553},{93,0.98541},{79,0.9845},{74,0.97595},{85,0.95036}};
//                 std::vector<std::pair<int, double>> fsdata =
//                 {
//                         {14,1.55228},{91,1.44918},{29,1.43568},{98,1.41855},{31,1.39284},
//                         {18,1.37799},{55,1.36194},{94,1.36182},{84,1.35785},{70,1.33855}
//                 };
                std::vector<std::pair<int, double>> fsdata ={
                        {54,1.89556},{23,1.84837},{8,1.75914},{7,1.73126},{30,1.73126},
                        {57,1.73126},{71,1.73126},{94,1.73126},{48,1.67853},{22,1.67062},
                        {24,1.65457},{38,1.64212},{61,1.64212},{2,1.63141},{35,1.61408},
                        {74,1.60421},{34,1.59933},{73,1.59933},{13,1.58809},{59,1.58809},
                        {98,1.58809},{9,1.57856},{17,1.57856},{52,1.57856},{64,1.57379},
                        {81,1.5528},{39,1.548},{78,1.53674},{93,1.53674},{42,1.52705},
                        {92,1.52705},{28,1.51172},{88,1.51172},{90,1.5033},{77,1.48952},
                        {43,1.47729},{27,1.46797},{50,1.45626},{65,1.44718},{69,1.44718},
                        {3,1.43937},{19,1.43937},{36,1.43937},{56,1.43937},{97,1.43937},
                        {80,1.43546},{10,1.42607},{46,1.42607},{75,1.42607},{0,1.41789},
                        {31,1.41789},{53,1.41789},{55,1.41789},{85,1.41789},{86,1.41383},
                        {58,1.41181},{21,1.41079},{6,1.4045},{20,1.4045},{26,1.4045},
                        {32,1.4045},{95,1.4045},{11,1.39658},{15,1.39658},{67,1.39262},
                        {5,1.38325},{14,1.38325},{76,1.38325},{4,1.37519},{51,1.36771},
                        {45,1.36184},{63,1.3538},{82,1.3538},{83,1.34978},{29,1.33417},
                        {1,1.32018},{89,1.32018},{70,1.31887},{18,1.31417},{25,1.31417},
                        {87,1.31417},{49,1.30909},{12,1.30067},{47,1.30067},{41,1.29306},
                        {60,1.2873},{79,1.2873},{91,1.28248},{66,1.27454},{62,1.26982},
                        {16,1.26746},{37,1.26746},{99,1.2619},{33,1.25475},{68,1.2492},
                        {40,1.24443},{44,1.24443},{84,1.2027},{72,1.15331},{96,1.15013}
                };//Da simulacao de monte carlo

        int levels=4;//numero de niveis da simulacao subset

        n=fsdata.size();
        //numeroro de cadeias de markov
        int nc = p0*n;

        //numero de amostras por cadeia
        int ns = 1/p0;
        cout << "p0 = "<< p0 <<" n = "<< n  << " nc = " << nc  << " ns = " << ns << endl;

        slopeanalysis->ResetSubSetSamples();

        //especificando os campos inicias da simulacao subset baseado nos menores valores de fs. Comecando do valor n-nc=10-5=5 ate n
        for ( int i=n-nc; i<n; i++ ) {
                cout << " fsdata[n-nc+1].first " << fsdata[i].first << endl;
                TPZVec<TPZFMatrix<REAL>> samples=slopeanalysis->GetIfield ( fsdata[i].first );
                slopeanalysis->SetSubSetSamples ( samples );

        }


        //Salvando
        TPZBFileStream save;
        save.OpenWrite ( "teste.bin" );
        slopeanalysis->Write ( save,slopeanalysis->ClassId() );
        REAL b = fsdata[n-nc].second;
        cout << "estimando a probabilidade de falha incial" << endl;
        cout<< b <<" " <<  p0     << endl;
        posprocfs3<< fsdata[n-nc].second<<" " <<  pow ( p0,0 ) *  nc/n << endl;


//         SlopeAnalysis* analysist = new SlopeAnalysis ( *slopeanalysis );
//                 //le os dados dos atributos salvos
//                 TPZBFileStream read;
//                 read.OpenRead ( "teste.bin" );
//                 analysist->Read ( read,0 );
//
//                 return;



        for ( int j=1; j<=levels; j++ ) { //levels

                //cria uma copia da analise para evitar modificaoes na estrutura interna(dados de pontos de integracao), copiando apenas atributos;
                SlopeAnalysis* analysis = new SlopeAnalysis ( *slopeanalysis );

                //le os dados dos atributos salvos
                TPZBFileStream read;
                read.OpenRead ( "teste.bin" );
                analysis->Read ( read,0 );

                //cria vetor para empilhar os valores dos fatores de seguraca calculados
                std::vector<double> fsvec;
                //cria vetor para empilhar os valores dos campos estocasticos correspondentes aos fatores de seguranca
                TPZStack<TPZVec<TPZFMatrix<REAL>>> outsamplesfull;


                //faz um loop sobre as cadeias de markov
                for ( int i=0; i<nc; i++ ) {

                        cout << " nc  = "<< i << endl;

                        //pega o campo selecionado na simualcao anterior com menores fs como semente
                        TPZVec<TPZFMatrix<REAL>> seed = analysis->GetSubSetSamples ( i );
                        //seed[0].Print(cout);
                        //return;
                        TPZStack<TPZVec<TPZFMatrix<REAL>>> outsamples;
                        //cada semente gera ns novos campos armazenados no outsamples utilizando o algoritmo de metropolis hastings
                        outsamples=  MetropolisHastings ( ns,seed,analysis,fsvec,b );
                        for ( int ins=0; ins<outsamples.size(); ins++ ) {
                                outsamplesfull.Push ( outsamples[ins] );
                        }

                        //outsamples[0][0].Print(cout);
                }

                cout << "\n fatores de seguranca =  " <<endl;
                for ( int ii=0; ii<fsvec.size(); ii++ ) cout <<fsvec[ii] << endl;
                //descarta os campos anteriores
                analysis->ResetSubSetSamples();

                posprocfs2 << "\n level =  "<< j <<endl;


                std::string saidafs2 = "postsubset/fs" + std::to_string ( j ) + ".dat";
                std::ofstream out2 ( saidafs2 );
                for ( int ii=0; ii<fsvec.size(); ii++ ) out2 <<fsvec[ii] << endl;

                //cria um pair para armazenar os fatores de seguranca e os indexes
                std::pair<std::vector<double>,std::vector<int>> out;

                //ordena em ordem decrescente e pega somente os ultimos n-nc=5 ate n=10 valores
                out=Sort ( fsvec,n-nc );

                cout << "szstack = "<< outsamplesfull.size() << " out.first.size() ="<< out.first.size() << endl;
                for ( int i = 0; i < out.first.size(); i++ ) {
                        cout << "First : " << out.first[i] << ", Second  : " << out.second[i] << endl;
                        analysis->SetSubSetSamples ( outsamplesfull[out.second[i]] );
                        posprocfs2 << out.first[i] <<endl;
                }

                std::string saidafs3 = "postsubset/pf" + std::to_string ( j ) + ".dat";
                std::ofstream out3 ( saidafs3 );
                out3 << b <<" " <<  pow ( p0,j ) *  nc/n << endl;
                b=out.first[0];
                cout << "estimando a probabilidade de falha" << endl;
                //estimando a probabilidade de falha
                cout<< b <<" " <<  pow ( p0,j ) *  nc/n << endl;
                posprocfs3<< b <<" " <<  pow ( p0,j ) *  nc/n << endl;

                TPZBFileStream save;
                save.OpenWrite ( "teste.bin" );
                analysis->Write ( save,analysis->ClassId() );
                delete analysis;

        }


}
TPZStack<TPZVec<TPZFMatrix<REAL>>> MetropolisHastings(int nnewsamples,TPZVec<TPZFMatrix<REAL>> seedinit,SlopeAnalysis* analysis,std::vector<double> &fsvec, REAL b)
{


        SlopeAnalysis* analysis3 = new SlopeAnalysis ( *analysis );

        REAL fsx0 =analysis3->SolveSingleField ( seedinit );

        delete analysis3;

        int M = analysis->GetM();

        std::normal_distribution<double> distribution ( 0., 1. );

        std::uniform_real_distribution<double> distribution2 ( -0.5, 0.5);

        std::uniform_real_distribution<double> distributionunif ( 0, 1. );

        TPZStack<TPZVec<TPZFMatrix<REAL>>> outsamples;

         for(int ins=0;ins<nnewsamples;ins++)
        {
                cout<< "-----------------------ins = "<< ins <<endl;
                TPZVec<TPZFMatrix<REAL>> newfield(2);
                newfield[0].Resize(M,1);
                newfield[1].Resize(M,1);

                 for ( int irdvar = 0; irdvar < M; irdvar++ )
                 {
                        std::random_device rd{};
                        std::mt19937 generator{ rd() };
                        std::random_device rd2{};
                        std::mt19937 generator2{ rd2() };
                        REAL xic = distribution ( generator );
                        newfield[0] ( irdvar,0 ) = xic+seedinit[0]( irdvar,0);
                        xic = distribution ( generator2 );
                        newfield[1] ( irdvar,0 ) = xic+seedinit[1]( irdvar,0);

                }

                SlopeAnalysis* analysis2 = new SlopeAnalysis ( *analysis );

                REAL fsx1 =analysis2->SolveSingleField(newfield );

                REAL alpha = min(fsx1/fsx0, 1.);

                //(*Aceitação ou rejeição*)
                std::random_device rd2{};

                std::mt19937 generator2{ rd2() };

                REAL u = distributionunif(generator2);


                if(u<alpha)
                {
                        if(fsx1<b)
                        {
                                cout<< "Aceita com "<<" b = " << b  << " fsx1 = " <<fsx1 << " fsx0 = " << fsx0  << " alpha " << alpha << " u = " << u << endl;
                                fsx0=fsx1;//aceita fs novo
                                seedinit=newfield;//atualiza campo novo
                                outsamples.Push(newfield);//aceita campo novo
                        }
                        else
                        {
                                cout<< "Rejeita 0 com "<<" b = " << b  << " fsx1 = " <<fsx1  << " alpha " << alpha << " u = " << u << endl;
                                fsx1=fsx0;//mantem fs antigo
                                outsamples.Push(seedinit);//descarta amostra e pega campo antigo
                        }

                }
                else
                {
                         cout<< "Rejeita 1 com "<<" b = " << b  << " fsx1 = " <<fsx1 << " fsx0 = " << fsx0   << " alpha " << alpha << " u = " << u << endl;
                        fsx1=fsx0;//mantem fs antigo
                        outsamples.Push(seedinit);//descarta amostra e pega campo antigo
                }

                fsvec.push_back(fsx1);

                delete analysis2;

        }

        return outsamples;
}
std::vector<std::pair<int, double>> CrudeMonteCarlo(int a,int b, SlopeAnalysis* slopeanalysis)
{
        std::vector<std::pair<int, double>> fsvec;
        for(int imc=a;imc<b;imc++)
        {
                cout << "imc  = "<<imc<<endl;
                SlopeAnalysis* analysis = new SlopeAnalysis ( *slopeanalysis );
                REAL fs = analysis->SolveSingleField ( imc );
                fsvec.emplace_back(imc, fs);
                posprocfs << imc << " "<<fs << endl;
//                 std::string saidavtk2 = "postvtk/saidavtk" + std::to_string ( imc ) + ".vtk";
//                 analysis->PostPlasticity ( saidavtk2 );
                std::string saidafs2 = "post/fs" + std::to_string ( imc ) + ".dat";
                std::ofstream out2 ( saidafs2 );
                out2 << fs << std::endl;
                delete analysis;
        }
        return fsvec;
}

void Write ( TPZStream &buf, int withclassid )
{
     //   fSolutionValVec.Write ( buf,withclassid );
    //    fFieldSamples[0].Write ( buf,withclassid );
    //    fFieldSamples[1].Write ( buf,withclassid );
    //    fFields[0].Write ( buf,withclassid );
     //   fFields[1].Write ( buf,withclassid );
//fHFields[0].Write ( buf,withclassid );
//fHFields[1].Write ( buf,withclassid );
//fPesos[0].Write ( buf,withclassid );
//fPesos[1].Write ( buf,withclassid );
}

void Read ( TPZStream &buf, void *context )
{

       // fSolutionValVec.Read ( buf,context );
      //  fFieldSamples.resize ( 2 );
      //  fFieldSamples[0].Read ( buf,context );
      //  fFieldSamples[1].Read ( buf,context );
     //   fFields.resize ( 2 );
      //  fFields[0].Read ( buf,context );
      //  fFields[1].Read ( buf,context );
//fHFields.resize ( 2 );
//fHFields[0].Read ( buf,context );
//fHFields[1].Read ( buf,context );
//fPesos.resize ( 2 );
        //fPesos[0].Read ( buf,context );
        //fPesos[1].Read ( buf,context );
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
