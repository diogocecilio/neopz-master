// Malha em camadas que seguem o terreno. Colunas verticais nas abscissas xs; cada coluna vai da base (y = 0) ate a
// superficie y_s(x) (a poligonal do terreno) e e dividida nas mesmas fracoes da altura, eta (0 na base, 1 na superficie):
// as linhas das camadas acompanham o topo, a face e o pe. Cada celula e um quadrilatero ou dois triangulos (pela diagonal
// menor; num empate, alternando como num xadrez). Contornos: elementos lineares no sentido anti-horario, como na
// TriGMesh (base da esquerda para a direita, lateral direita para cima, superficie da direita para a esquerda, lateral
// esquerda para baixo), com um id por trecho da superficie. Malha plana, sem hierarquia; com triangulos, serve de base
// para a RefinedGeoMesh
#pragma once

#include "pzgmesh.h"

#include <array>
#include <vector>

struct LayerGMeshData {
    std::vector<std::array<REAL, 2>> superficie; // vertices do terreno (x, y), da esquerda para a direita
    std::vector<int> idSuperficie;               // id de cada trecho da superficie (um a menos que os vertices)
    int idMaterial = 0, idBase = 0, idEsquerda = 0, idDireita = 0;
    std::vector<REAL> xs;  // abscissas das colunas, crescentes, da primeira a ultima abscissa da superficie e passando
                           // por todos os vertices dela
    std::vector<REAL> eta; // fracoes da altura das linhas das camadas, crescentes: 0 (base), ..., 1 (superficie)
    bool quadrilateros = false;
};

TPZGeoMesh *LayerGMesh(const LayerGMeshData &d);

// Passos que somam L, comecando em h (o lado fino) e crescendo na razao r ate hmax. O ultimo passo raramente fecha a
// conta: todos sao esticados ou encolhidos na mesma proporcao (com ou sem o ultimo, o que mexer menos)
std::vector<REAL> GeometricSteps(REAL L, REAL h, REAL r, REAL hmax);
