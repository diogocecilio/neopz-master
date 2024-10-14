// SPDX-FileCopyrightText: 2024 <copyright holder> <email>
// SPDX-License-Identifier: Apache-2.0

#ifndef SLOPEANALYSIS_H
#define SLOPEANALYSIS_H

#include "slopeconfigure.h"
#include <time.h>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <chrono>
#include <fstream>
#include "TPZRandomFieldAnalysis.h"
#include "TPZKarhunenLoeveMat.h"
#include "pzcmesh.h"
#include "tpzgeoelrefpattern.h"
#include "TPZVTKGeoMesh.h"
#include "TPZFileStream.h"
#include "TPZBFileStream.h"
/**
 * @todo write docs
 */
class slopeanalysis: Slope
{

   // friend Slope;
public:
    /**
     * Default constructor
     */
    slopeanalysis();

    slopeanalysis( TPZGeoMesh * gmesh,int porder,int ref, REAL gammaagua, REAL gammasolo,REAL coes,REAL atrito);

    /**
     * Copy constructor
     *
     * @param other TODO
     */
    slopeanalysis(const slopeanalysis& other);

    /**
     * Destructor
     */
    ~slopeanalysis();

    void CreateFields();

    TPZCompMesh* CreateCompMeshKL ( TPZGeoMesh * gmesh );

    TPZGeoMesh * TriGMesh(int ref);

protected:

    TPZCompMesh*fCompMeshSlope;
    TPZRandomFieldAnalysis*fFieldAnalysis;

};

#endif // SLOPEANALYSIS_H
