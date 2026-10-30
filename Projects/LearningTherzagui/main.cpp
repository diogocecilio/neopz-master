// LearningTherzagui: coluna de Terzaghi resolvida passo a passo (projeto de estudo).
// Por enquanto so os includes; o programa sera escrito aqui.

// geometria
#include "TPZGenGrid2D.h"
#include "TPZGeoLinear.h"
#include "pzgeoel.h"
#include "pzgeoquad.h"
#include "pzgeotriangle.h"
#include "pzgmesh.h"
#include "tpzgeoelrefpattern.h"

// malhas computacionais
#include "TPZBndCondT.h"
#include "TPZMultiphysicsCompMesh.h"
#include "TPZNullMaterial.h"
#include "pzcmesh.h"

// material u - p_w - p_a e Mohr-Coulomb PV2 do artigo
#include "Plasticity/TPZElasticResponse.h"
#include "Plasticity/TPZMatUnsaturatedUPwPa.h"
#include "Plasticity/TPZPlasticStepVoigt.h"
#include "Plasticity/TPZUnsaturatedHydraulics.h"
#include "Plasticity/TPZYCMohrCoulombPV2.h"

// analise e solver
#include "TPZLinearAnalysis.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "TPZSpStructMatrix.h"
#include "pzstepsolver.h"

#include <cmath>
#include <fstream>
#include <iostream>

int main() {
    return 0;
}
