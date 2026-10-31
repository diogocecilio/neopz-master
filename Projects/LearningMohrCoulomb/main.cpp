// LearningMohrCoulomb: verificacao do Mohr-Coulomb (PV2) em deformacao plana, passo a passo
// (projeto de estudo): capacidade de carga de uma sapata corrida e estabilidade de um talude.
// Por enquanto so os includes; o programa sera escrito aqui.

// geometria
#include "TPZGmshReader.h"
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

// ids de material: grupos fisicos das malhas do Gmsh (meshes/talude.geo e meshes/sapata.geo)
enum EIds {
    matSolo = 1,       // solo (dominio)
    bcBase = -1,       // base
    bcDireita = -2,    // lateral direita
    bcSuperficie = -3, // superficie livre do terreno
    bcEsquerda = -4,   // lateral esquerda (talude) ou eixo de simetria (sapata)
    bcSapata = -5      // contato da sapata
};

// Le uma malha do Gmsh (.msh 4.1) e associa os nomes dos grupos fisicos aos ids de material.
//   domains    : grupos de dimensao 2 (dominio)
//   boundaries : grupos de dimensao 1 (contornos)
// O TPZGmshReader (formato 4) cria os elementos com o numero do grupo fisico do Gmsh, sem aplicar
// o mapa de nomes; por isso os ids sao trocados aqui pelo nome do grupo.
TPZGeoMesh *ReadGmshMesh(const std::string &file, const std::map<std::string, int> &domains,
                         const std::map<std::string, int> &boundaries)
{
    TPZGmshReader reader;
    TPZGeoMesh *gmesh = reader.GeometricGmshMesh(file);
    auto &tagName = reader.GetDimPhysicalTagName(); // [dimensao][numero do Gmsh] -> nome
    for (TPZGeoEl *gel : gmesh->ElementVec()) {
        if (!gel) continue;
        const int dim = gel->Dimension();
        const auto &ids = (dim == 2) ? domains : boundaries;
        auto name = tagName[dim].find(gel->MaterialId());
        auto id = (name == tagName[dim].end()) ? ids.end() : ids.find(name->second);
        if (dim == 0 || id == ids.end()) {
            std::cout << "ReadGmshMesh: grupo fisico " << gel->MaterialId() << " de dimensao " << dim
                      << " sem id de material" << std::endl;
            DebugStop();
        }
        gel->SetMaterialId(id->second);
    }
    gmesh->BuildConnectivity();
    return gmesh;
}

// Talude: 75 m x 40/30 m, h = 10 m, 45 graus
TPZGeoMesh *ReadSlopeMesh()
{
    return ReadGmshMesh(std::string(MESH_DIR) + "/talude.msh", {{"Solo", matSolo}},
                        {{"Base", bcBase}, {"Direita", bcDireita}, {"Superficie", bcSuperficie}, {"Esquerda", bcEsquerda}});
}

// Sapata corrida: 5 m x 5 m, meia sapata de 0,5 m com simetria em x = 0
TPZGeoMesh *ReadFootingMesh()
{
    return ReadGmshMesh(std::string(MESH_DIR) + "/sapata.msh", {{"Solo", matSolo}},
                        {{"Base", bcBase}, {"Direita", bcDireita}, {"Superficie", bcSuperficie},
                         {"Simetria", bcEsquerda}, {"Sapata", bcSapata}});
}

int main() {
    for (auto &[name, read] : std::map<std::string, std::function<TPZGeoMesh *()>>{{"talude", ReadSlopeMesh},
                                                                                     {"sapata", ReadFootingMesh}}) {
        TPZGeoMesh *gmesh = read();
        std::map<int, int> count;
        for (TPZGeoEl *gel : gmesh->ElementVec())
            if (gel) count[gel->MaterialId()]++;
        std::cout << name << ": " << gmesh->NNodes() << " nos";
        for (auto &[id, n] : count) std::cout << ", id " << id << ": " << n;
        std::cout << std::endl;
        std::ofstream out(name + ".vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out, true);
        delete gmesh;
    }
    return 0;
}
