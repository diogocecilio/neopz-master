// SPDX-FileCopyrightText: 2024 <diogo cecilio> <diogo.cecilio@ufrgs.br>
// SPDX-License-Identifier: Apache-2.0

#include "SubsetMonteCarlo.h"



SubsetMonteCarlo::SubsetMonteCarlo() : fCurrentConfig(), fSequence(), fSlopeAnalysis(),fp0(),fNsamples()
{

}

SubsetMonteCarlo::SubsetMonteCarlo(SlopeAnalysis * analysis, REAL p0, int samples): fCurrentConfig(), fSequence(), fSlopeAnalysis(analysis),fp0(p0),fNsamples(samples)
{

}

SubsetMonteCarlo::SubsetMonteCarlo(const SubsetMonteCarlo &copy) : fCurrentConfig(copy.fCurrentConfig), fSequence(copy.fSequence),fSlopeAnalysis(copy.fSlopeAnalysis),fp0(copy.fp0),fNsamples(copy.fNsamples)
{

}

SubsetMonteCarlo &SubsetMonteCarlo::operator=(const SubsetMonteCarlo &copy)
{
    if (this == &copy) {
        return *this;
    }
    fCurrentConfig = copy.fCurrentConfig;
    fSequence = copy.fSequence;
    fSlopeAnalysis= copy.fSlopeAnalysis;
    fp0=copy.fp0;
    fNsamples=copy.fNsamples;
    return *this;
}

SubsetMonteCarlo::~SubsetMonteCarlo()
{

    delete fSlopeAnalysis;
}

/// write the object on the stream
void SubsetMonteCarlo::Write(TPZStream &buf, int withclassid) const
{
    fCurrentConfig.Write(buf,withclassid);
    int seqsize = fSequence.size();
    buf.Write(&seqsize);
    //std::list<SubsetMonteCarlo::TConfig>::iterator it;
    auto it = fSequence.begin();
    for (it = fSequence.begin(); it != fSequence.end(); it++)
    {
        it->Write(buf,withclassid);
    }

}

/// read the object from the stream
void SubsetMonteCarlo::Read(TPZStream &buf, void *context)
{
    fCurrentConfig.Read(buf,context);
    int seqsize;
    buf.Read(&seqsize);
    for (int i=0; i<seqsize; i++) {
        SubsetMonteCarlo::TConfig config;
        config.Read(buf,context);
        fSequence.push_back(config);
    }

}




SubsetMonteCarlo::TConfig::TConfig(): fSimulateFields(),fHistoryLog()
{

}

SubsetMonteCarlo::TConfig::TConfig(const TConfig &conf) : fSimulateFields(conf.fSimulateFields),fHistoryLog(conf.fHistoryLog)
{

}

SubsetMonteCarlo::TConfig::~TConfig()
{
        // Destrutor, limpando a memória alocada

}

SubsetMonteCarlo::TConfig &SubsetMonteCarlo::TConfig::operator=(const SubsetMonteCarlo::TConfig &copy)
{
    if (this == &copy) {
        return *this;
    }
    fSimulateFields = copy.fSimulateFields;
    fHistoryLog=copy.fHistoryLog;
    return *this;
}

/// Write the data to the output stream
void SubsetMonteCarlo::TConfig::Write(TPZStream &buf, int withclassid) const
{
        //buf.Write(&fHistoryLog);
        int nc = fSimulateFields.size();
        int sz=fSimulateFields[0].second.size();

		buf.Write(&nc);
        buf.Write(&sz);

        cout << "Write"<<endl;
        cout << "nc = "<< nc << " sz = "<< sz <<endl;
		for(int c=0; c<nc; c++)
        {
                buf.Write(&fSimulateFields[c].first);
                for(int j=0;j<sz;j++)
                {
                        fSimulateFields[c].second[j].Write(buf,withclassid);
                }

        }

}

/// Read the data from the input stream
void SubsetMonteCarlo::TConfig::Read(TPZStream &buf, void *context)
{

        //buf.Read(&fHistoryLog);
        int nc,sz;
		buf.Read(&nc);
        buf.Read(&sz);

        fSimulateFields.resize(nc);

        cout << "read"<<endl;
        cout << "nc = "<< nc << " sz = "<< sz <<endl;

		for(int c=0; c<nc; c++)
        {
                buf.Read(&fSimulateFields[c].first);
                fSimulateFields[c].second.resize(sz);
                for(int j=0;j<sz;j++)
                {
                        fSimulateFields[c].second[j].Read ( buf,context );
                }
        }

}

void SubsetMonteCarlo::SaveConfig(std::stringstream &strout)
{
    fCurrentConfig.fHistoryLog = strout.str();
    fSequence.push_back(fCurrentConfig);
}




//#include "matplotlibcpp.h"

void SubsetMonteCarlo::SubSet( )
{

    int level=1;
    if(level==0)
    {
        ExecuteInitialMonteCarloSimulation(0,100);
        TPZBFileStream save;
        save.OpenWrite("SubSetInitialMonteCarloConfig.bin");
        Write(save,ClassId());
    }

    if(level==1)
    {
        TPZBFileStream read2;
        read2.OpenRead ( "SubSetInitialMonteCarloConfig.bin" );
        Read ( read2,0 );

        TPZStack<std::pair<REAL,TPZVec<TPZFMatrix<REAL>>>> copy = fCurrentConfig.fSimulateFields;

        fCurrentConfig.fSimulateFields = LevelLoop(copy);
        std::string name =  "Level"+ std::to_string ( 1 );
        std::stringstream strout;
        strout << name;
        SaveConfig(strout);
        TPZBFileStream save;
        name+=".bin";
        save.OpenWrite(name);
        Write(save,ClassId());
    }

    if(level==2)
    {
        TPZBFileStream read2;
        read2.OpenRead ( "Level1.bin" );
        Read ( read2,0 );

        TPZStack<std::pair<REAL,TPZVec<TPZFMatrix<REAL>>>> copy = fCurrentConfig.fSimulateFields;

        fCurrentConfig.fSimulateFields = LevelLoop(copy);
        std::string name =  "Level"+ std::to_string ( 2 );
        std::stringstream strout;
        strout << name;
        SaveConfig(strout);
        TPZBFileStream save;
        name+=".bin";
        save.OpenWrite(name);
        Write(save,ClassId());
    }

    if(level==3)
    {
        TPZBFileStream read2;
        read2.OpenRead ( "Level2.bin" );
        Read ( read2,0 );

        TPZStack<std::pair<REAL,TPZVec<TPZFMatrix<REAL>>>> copy = fCurrentConfig.fSimulateFields;

        fCurrentConfig.fSimulateFields = LevelLoop(copy);
        std::string name =  "Level"+ std::to_string ( 3 );
        std::stringstream strout;
        strout << name;
        SaveConfig(strout);
        TPZBFileStream save;
        name+=".bin";
        save.OpenWrite(name);
        Write(save,ClassId());
    }
    if(level==4)
    {
        cout<< "ads"<<endl;
        TPZBFileStream read2;
        read2.OpenRead ( "Level3.bin" );
        Read ( read2,0 );

        TPZStack<TPZFMatrix<REAL>> oudata;

        cout<< " -- level 0 -- "<<endl;
        TConfig *conf0 =  GetConfig (0);
        REAL prod=fp0;
        TPZFMatrix<REAL> pfdata=  ComputePf(conf0->fSimulateFields,prod);

        oudata.Push(pfdata);
        cout<< " -- level 1 -- "<<endl;
        TConfig *conf1 =  GetConfig (1);

        prod=pow(fp0,2);
        pfdata= ComputePf(conf1->fSimulateFields,prod);

        oudata.Push(pfdata);
        cout<< " -- level 2 -- "<<endl;
        TConfig *conf2 =  GetConfig (2);

        prod=pow(fp0,3);
        pfdata= ComputePf(conf2->fSimulateFields,prod);

        oudata.Push(pfdata);
        cout<< " -- level 3 -- "<<endl;
        TConfig *conf3 =  GetConfig (3);

        prod=pow(fp0,4);
        pfdata= ComputePf(conf3->fSimulateFields,prod);
        oudata.Push(pfdata);

        std::ofstream posprocfs ( "posprocfs.txt" );
        for(int  i=0;i<oudata.size();i++)
        {
            for(int j=0;j<oudata[i].Rows();j++)
            {
                posprocfs << oudata[i](j,0) << " " << oudata[i](j,1) <<endl;
            }
        }
    }

}
TPZStack<std::pair<REAL,TPZVec<TPZFMatrix<REAL>>>>  SubsetMonteCarlo::LevelLoop(TPZStack<std::pair<REAL,TPZVec<TPZFMatrix<REAL>>>> copy )
{

        TPZStack<std::pair<REAL,TPZVec<TPZFMatrix<REAL>>>> newSimulateFields;
        int n= copy.size();

        SortData ( copy );


        REAL p0=fp0;

        REAL nc = p0*n;
        REAL ns = 1/p0;
        cout << "p0 = "<< p0 << endl;
        cout << "n = "<< n << endl;
        cout << "fNsamples = "<< fNsamples << endl;
        cout << "nc = "<< nc << endl;
        cout << "ns = "<< ns<<endl;
        REAL b =copy[n-nc].first;


        TPZStack<std::pair<REAL,TPZVec<TPZFMatrix<REAL>>>> outsamplesfull;

        for ( int i=0; i<nc; i++ ) {
                std::pair<REAL,TPZVec<TPZFMatrix<REAL>>> temp;
                TPZStack<std::pair<REAL,TPZVec<TPZFMatrix<REAL>>>> outsamples;
                outsamples=  MetropolisHastings ( ns, copy[n-nc+i],b );
                for ( int ins=0; ins<outsamples.size(); ins++ ) {
                        outsamplesfull.Push ( outsamples[ins] );
                }
        }

        return outsamplesfull;
}

TPZStack<std::pair<REAL,TPZVec<TPZFMatrix<REAL>>>> SubsetMonteCarlo::MetropolisHastings(int nnewsamples,std::pair<REAL,TPZVec<TPZFMatrix<REAL>>> seedinit,REAL b)
{


        TPZStack<std::pair<REAL,TPZVec<TPZFMatrix<REAL>>>> copy = fCurrentConfig.fSimulateFields;

        SortData(copy);

        int M = fSlopeAnalysis->GetM();

        cout << "MetropolisHastings with proposal distribution ~N(0.,1) << " << endl;
        std::normal_distribution<double> distribution ( 0., 1. );

        std::uniform_real_distribution<double> distribution2 ( -0.5, 0.5);

        std::uniform_real_distribution<double> distributionunif ( 0, 1. );

        TPZStack<std::pair<REAL,TPZVec<TPZFMatrix<REAL>>>> outsamples;

        REAL fsx0 = seedinit.first;

         for(int ins=0;ins<nnewsamples;ins++)
        {
                std::pair<REAL,TPZVec<TPZFMatrix<REAL>>> couple;
                cout<< "ins = "<< ins <<endl;
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
                        newfield[0] ( irdvar,0 ) = xic+seedinit.second[0]( irdvar,0);
                        xic = distribution ( generator2 );
                        newfield[1] ( irdvar,0 ) = xic+seedinit.second[1]( irdvar,0);

                }

                SlopeAnalysis* analysis2 = new SlopeAnalysis ( *fSlopeAnalysis );

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
                                seedinit.second=newfield;//atualiza campo novo
                                couple.first=fsx0;
                                couple.second = newfield;
                                outsamples.Push(couple);//aceita campo novo
                       }
                        else
                       {
                               cout<< "Rejeita 0 com "<<" b = " << b  << " fsx1 = " <<fsx1  << " alpha " << alpha << " u = " << u << endl;
                               fsx1=fsx0;//mantem fs antigo
                               couple.first=fsx1;
                               couple.second = seedinit.second;
                               outsamples.Push(couple);//descarta amostra e pega campo antigo
                       }

                }
                else
                {
                         cout<< "Rejeita 1 com "<<" b = " << b  << " fsx1 = " <<fsx1 << " fsx0 = " << fsx0   << " alpha " << alpha << " u = " << u << endl;
                        fsx1=fsx0;//mantem fs antigo
                        couple.first=fsx1;
                        couple.second = seedinit.second;
                        outsamples.Push(couple);//descarta amostra e pega campo antigo


                }

                delete analysis2;
        }

        return outsamples;
}
void SubsetMonteCarlo::ExecuteInitialMonteCarloSimulation ( int a, int b )
{
        std::stringstream strout;

        std::vector<std::pair<int, double>> fsvec;
        for(int imc=a;imc<b;imc++)
        {

                cout << "imc  = "<<imc<<endl;

                SlopeAnalysis* analysis = new SlopeAnalysis ( *fSlopeAnalysis );

                TPZVec<TPZFMatrix<REAL>> field = analysis->GetIfield (imc );

                std::pair<REAL ,TPZVec<TPZFMatrix<REAL>> > pairdata;


                REAL fs = analysis->SolveSingleField ( field );

                pairdata.first=fs;
                pairdata.second=field;


                fCurrentConfig.fSimulateFields.Push(pairdata);

                delete analysis;
        }
        strout << "Initial Monte Carlo Simulation";
        SaveConfig(strout);

}
void SubsetMonteCarlo::SortData(TPZStack<std::pair<REAL,TPZVec<TPZFMatrix<REAL>>>> & data)
{
    int sizeofstack=data.size();
    for(int i=0;i<sizeofstack;i++)
    {
        for(int j=0;j<sizeofstack-1;j++)
        {
            if(data[j].first<data[j+1].first)
            {
                REAL temp = data[j].first;
                TPZVec<TPZFMatrix<REAL>> temp2 = data[j].second;

                data[j].first=data[j+1].first;
                data[j].second=data[j+1].second;

                data[j+1].first=temp;
                data[j+1].second=temp2;
            }
        }
    }

//      for(int i=0;i<sizeofstack;i++)
//     {
//         cout<<  data[i].first << endl;
//     }

}
int SubsetMonteCarlo::TConfig::ClassId() const
{
        return Hash ( "SubsetMonteCarloTConfig" );
}

int SubsetMonteCarlo::ClassId() const
{
        return Hash ( "SubsetMonteCarlo" );
}
