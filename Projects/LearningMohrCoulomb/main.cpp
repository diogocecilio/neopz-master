// LearningMohrCoulomb: verificacao do Mohr-Coulomb (PV2) em deformacao plana, passo a passo
// (projeto de estudo): capacidade de carga de uma sapata corrida e estabilidade de um talude.
// Por enquanto so os includes; o programa sera escrito aqui.

// geometria
#include "TPZGenGrid2D.h"
#include "TPZGeoLinear.h"
#include "pzgeoel.h"
#include "pzgeoelbc.h"
#include "pzgeoquad.h"
#include "pzgeotriangle.h"
#include "pzgmesh.h"
#include "tpzgeoelrefpattern.h"

// malhas computacionais
#include "TPZBndCondT.h"
#include "TPZMultiphysicsCompMesh.h"
#include "TPZNullMaterial.h"
#include "pzcmesh.h"

// material u - p com tensao efetiva elastoplastica e Mohr-Coulomb PV2 do artigo
#include "Plasticity/TPZElasticResponse.h"
#include "Plasticity/TPZMatPoroElastoPlasticUP.h"
#include "Plasticity/TPZPlasticStepVoigt.h"
#include "Plasticity/TPZYCMohrCoulombPV2.h"

// analise incremental (Newton, bissecao, deslocamento controlado e reacoes) e solver
#include "Plasticity/TPZPoroElastoPlasticUPAnalysis.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "pzstepsolver.h"

// pos-processamento
#include "TPZVTKGeoMesh.h"

#include <array>
#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>

int main() {
    return 0;
}
