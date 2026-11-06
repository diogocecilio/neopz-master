// Incremental analysis of the u - p_w - p_a model: the class moved to the library (Analysis/TPZUnsatAnalysis.h);
// the name UnsatAnalysis is kept for the test projects.
#ifndef UNSATANALYSIS_H
#define UNSATANALYSIS_H

#include "UnsatModel.h"
#include "TPZUnsatAnalysis.h"
#include "TPZEquationFilter.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "TPZSpStructMatrix.h"
#include "pzstepsolver.h"

#include <chrono>
#include <iostream>

using UnsatAnalysis = TPZUnsatAnalysis;

#endif
