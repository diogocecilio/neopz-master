// SPDX-FileCopyrightText: 2024 <diogo cecilio> <diogo.cecilio@ufrgs.br>
// SPDX-License-Identifier: Apache-2.0

#ifndef SUBSETMONTECARLO_H
#define SUBSETMONTECARLO_H

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
//std::ofstream pfout ( "pfsubset.txt" );
class SubsetMonteCarlo
{

    struct TConfig
    {
            TConfig();

            ~TConfig();

            TConfig(const TConfig &copy);

            TConfig &operator=(const TConfig &copy);

            void Write(TPZStream &buf, int withclassid) const;

            void Read(TPZStream &buf, void *context);

            int ClassId() const; // Identificador de classe para serialização

            //guarda os campos simulados e seus respectivos FSs
            TPZStack<std::pair<REAL,TPZVec<TPZFMatrix<REAL>>>> fSimulateFields;

            std::string fHistoryLog;

    };

public:
    SubsetMonteCarlo();

    SubsetMonteCarlo(SlopeAnalysis * analysis, REAL p0, int samples);

    ~SubsetMonteCarlo();

    SubsetMonteCarlo(const SubsetMonteCarlo &copy);

    SubsetMonteCarlo &operator=(const SubsetMonteCarlo &copy);

    void SetSlopeAnalysis(SlopeAnalysis * analysis)
    {
        fSlopeAnalysis=analysis;
    }

    void Write(TPZStream &buf, int withclassid) const;

    void Read(TPZStream &buf, void *context);

    void SortData(TPZStack<std::pair<REAL,TPZVec<TPZFMatrix<REAL>>>> &data);

    //TPZStack<TPZVec<TPZFMatrix<REAL>>> MetropolisHastings(int nnewsamples,TPZVec<TPZFMatrix<REAL>> seedinit,std::vector<double> &fsvec, REAL b);

    TPZStack<std::pair<REAL,TPZVec<TPZFMatrix<REAL>>>> MetropolisHastings(int nnewsamples,std::pair<REAL,TPZVec<TPZFMatrix<REAL>>> seedinit,REAL b);

    TPZStack<std::pair<REAL,TPZVec<TPZFMatrix<REAL>>>>  LevelLoop( TPZStack<std::pair<REAL,TPZVec<TPZFMatrix<REAL>>>>);

    void ExecuteInitialMonteCarloSimulation ( int a, int b );

    void SaveConfig(std::stringstream &strout);


    void SubSet();

    int ClassId() const;
    TConfig * GetCurrentConfig ()
    {
        return &fCurrentConfig;
    }

    void PopConfiguration()
    {
        if (this->fSequence.size() ==0) {
            return;
        }
        fSequence.pop_back();
        fCurrentConfig = *(this->fSequence.rbegin());
    }

    /// Access method
    TConfig * GetConfig (int index) {
        if (index < 0 || index >= fSequence.size())
            DebugStop();

        list<SubsetMonteCarlo::TConfig>::iterator inte;
        int i=0;
        for (inte=fSequence.begin(); inte!=fSequence.end(); ++inte, i++)
        {
            if (i == index)
            {
                return &(*inte);
            }
        }
        DebugStop();
        return 0;
    }

    /// Return size of config list
    int GetConfigListSize () {
        return fSequence.size();
    }

    void PostProcessPf(string filename);

   TPZFMatrix<REAL> ComputePf(TPZStack<std::pair<REAL,TPZVec<TPZFMatrix<REAL>>>> copy,REAL &prod0,int index)
    {

        SortData(copy);
        int n= copy.size();

        REAL p0=fp0;
        REAL nc = p0*n;
        TPZFMatrix<REAL> mat(n-nc,2);

        for(int i=0;i<n-nc;i++)
        {

            REAL val1=copy[i+1].first;
            REAL val2=copy[i].first;
            REAL count1=1.;
            REAL count2=1.;
                for (int j =0;j<n;j++)
                {
                    if (val1 >= copy[j].first)
                    {
                        count1++;
                    }
                    if (val2 >= copy[j].first)
                    {
                        count2++;
                    }
                }


            prod0*=count1/count2;
           // cout <<"val1 = " << val1 <<" val2 = " << val2 << " count1 = " << count1 << " count2 = " << count2 <<endl;
            cout << val1 << " " << prod0 <<endl;
                mat(i,0)=val1;
                mat(i,1)=prod0;
        }
        return mat;
    }


private:

    /// The object with the current configuration
    TConfig fCurrentConfig;

    /// The list of all previous configurations
    std::list<TConfig> fSequence;

    SlopeAnalysis * fSlopeAnalysis;

    REAL fp0;

    REAL fNsamples;

};

#endif // SUBSETMONTECARLO_H
