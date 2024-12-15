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

            REAL fp0;

            REAL fNsamples;




    };

public:
    SubsetMonteCarlo();

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

    TPZStack<std::pair<REAL,TPZVec<TPZFMatrix<REAL>>>> MetropolisHastings(int nnewsamples,TPZVec<TPZFMatrix<REAL>> seedinit,REAL b);

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


    void ComputePf(TPZStack<std::pair<REAL,TPZVec<TPZFMatrix<REAL>>>> copy,REAL prod0)
    {
        SortData(copy);
        for(int i=0;i<copy.size()-1;i++)
        {

            REAL val1=copy[i+1].first;
            REAL val2=copy[i].first;
            REAL count1=1.;
            REAL count2=1.;
                for (int j =0;j<copy.size();j++)
                {
                    if (val1 <= copy[j].first)
                    {
                        count1++;
                    }
                    if (val2 <= copy[j].first)
                    {
                        count2++;
                    }
                }
            prod0*=count2/count1;
            //cout <<"val1 = " << val1 <<" val2 = " << val2 << " count1/count2 = " << count1/count2 <<endl;
            cout << val1 << " " << prod0 <<endl;

        }
    }


private:

    /// The object with the current configuration
    TConfig fCurrentConfig;

    /// The list of all previous configurations
    std::list<TConfig> fSequence;

    SlopeAnalysis * fSlopeAnalysis;

};

#endif // SUBSETMONTECARLO_H
