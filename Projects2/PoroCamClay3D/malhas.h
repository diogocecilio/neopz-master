// malhas.h — malhas geométricas do NeoPZ (TPZGeoMesh) para os exemplos u-p com Cam-Clay
#ifndef POROCAMCLAY_MALHAS_H
#define POROCAMCLAY_MALHAS_H

#include <array>
#include "pzgmesh.h"

/// Material id dos elementos 3D
constexpr int kMatVolume = 100;

/// Caixa [0,Lx]x[0,Ly]x[0,Lz] com nx x ny x nz hexaedros: TPZGeoCube (8 nós) ou, com quadratic = true,
/// TPZQuadraticCube (20 nós). Faces de contorno (TPZGeoQuad / TPZQuadraticQuad) com material id:
/// 1 x = 0, 2 x = Lx, 3 y = 0, 4 y = Ly, 5 z = 0, 6 z = Lz  (marcadores do box_mesh do Python).
TPZGeoMesh *BoxMesh(const std::array<REAL, 3> &L, const std::array<int, 3> &n, bool quadratic);

/// Marcadores das faces do quarto de cilindro
enum ECylinderFace { EX0 = 1, ELateral = 2, EY0 = 3, EBase = 5, ETopo = 6 };

/// Quarto de cilindro 0 <= z <= h (simetria em x = 0 e y = 0) com malha "O-grid": quadrado central
/// [0, a]² (nc x nc, a = aratio r) e dois blocos externos (nc ao longo do arco de 45°, nr na direção
/// radial), nz camadas. Com quadratic = true os nós de meio de aresta da superfície lateral ficam no arco.
TPZGeoMesh *QuarterCylinderMesh(REAL r, REAL h, int nc, int nr, int nz, bool quadratic, REAL aratio = 0.5);

#endif
