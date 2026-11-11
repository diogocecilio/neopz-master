// Malhas do talude: leitura do gmsh, malha plana (sem hierarquia) e refinamentos uniforme, Rivara (h) e de ordem (p)
#ifndef SLOPEMESH_H
#define SLOPEMESH_H

#include "pzgmesh.h"

#include <array>
#include <map>
#include <set>
#include <string>
#include <vector>

namespace slope {

// Ids de material: grupos fisicos do .geo e gemeos das laterais (poropressao da malha u-p)
enum EIds {
    matSolo = 1,       // Physical Surface("Solo")
    bcBase = -1,       // Physical Curve("Base")
    bcDireita = -2,    // Physical Curve("Direita")
    bcSuperficie = -3, // Physical Curve("Superficie"): terreno do pe, face e crista
    bcEsquerda = -4,   // Physical Curve("Esquerda")
    bcDireitaP = -12,  // gemeos (DuplicateBoundaryBC)
    bcEsquerdaP = -14
};

// Le o .msh pelos nomes dos grupos fisicos; para se algum elemento ficar sem id conhecido
TPZGeoMesh *ReadGMesh(const std::string &file);

// Gemeos das laterais: o original leva u_x = 0 e o gemeo a poropressao. Chamar antes das malhas computacionais
void DuplicateBoundaryBC(TPZGeoMesh *gmesh);

// Malha sem hierarquia: cada nivel vira uma TPZGeoMesh nova e conforme (as analises nao aceitam nos pendurados)
struct FlatMesh {
    std::vector<std::array<REAL, 2>> nodes;
    std::vector<std::array<int64_t, 3>> tri;       // triangulos do solo, na ordem dos elementos da TPZGeoMesh
    std::map<std::pair<int64_t, int64_t>, int> bc; // aresta de contorno (nos em ordem crescente) -> id
    std::vector<int> order;                        // ordem de u de cada triangulo
};

FlatMesh FromGeoMesh(TPZGeoMesh *gmesh, int order = 2);
TPZGeoMesh *ToGeoMesh(const FlatMesh &m); // triangulos primeiro: o elemento i e m.tri[i]
FlatMesh RefineUniform(const FlatMesh &in); // 4 triangulos semelhantes; os filhos herdam a ordem
// Bisseccao pela aresta mais longa (Rivara): marcados divididos duas vezes (h/2), vizinhos so para a conformidade
FlatMesh RefineMarked(const FlatMesh &in, const std::set<int64_t> &marked);

// VTK dos triangulos com escalares por celula e por no (nome -> valores)
void WriteFlatVTK(const FlatMesh &m, const std::map<std::string, std::vector<REAL>> &cellData,
                  const std::map<std::string, std::vector<REAL>> &pointData, const std::string &file);

} // namespace slope

#endif
