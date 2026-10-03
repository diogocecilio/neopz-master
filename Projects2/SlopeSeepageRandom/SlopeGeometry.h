// SlopeGeometry.h
//
// Geometria do talude dos exemplos de Vargas Ceron, Cecílio, Linn & Maghous (IJNAMG 2025) e de Cho (2010):
//
//            crista (y = D)
//   (0,D) +-----------------+ (Lc, D)
//         |                  \  face (inclinação β)
//         |                   \
//         |    bloco B         \ (Lc+Ls, Hb)      pé (y = Hb)
//   (0,Hb)+--------------------+--------------------+ (W, Hb)
//         |                    bloco A               |
//   (0,0) +----------------------------------------- + (W, 0)
//
//   H = altura do talude, Ls = H/tan β, Lc = largura da crista, Lt = largura do pé,
//   Hb = espessura abaixo do pé, D = Hb + H, W = Lc + Ls + Lt.
//
// A malha é estruturada (quadriláteros de 4 nós ou triângulos), com dois blocos que compartilham os nós em
// y = Hb; no bloco B as colunas se afunilam até a face. O mesmo TPZGeoMesh serve às malhas de KL (campos
// aleatórios), de fluxo (Darcy) e mecânica (elastoplástica).
//
#ifndef SLOPEGEOMETRY_H
#define SLOPEGEOMETRY_H

#include <string>

#include "pzgmesh.h"

/// Identificadores de material do TPZGeoMesh
enum ESlopeIds {
    ESoil = 1,       ///< solo
    EBottom = -1,    ///< base (y = 0)
    ELeft = -2,      ///< lado esquerdo (x = 0)
    ERight = -3,     ///< lado direito (x = W)
    ECrest = -4,     ///< crista (y = D)
    EFace = -5,      ///< face do talude
    EToe = -6        ///< superfície do pé (y = Hb, x > Lc + Ls)
};

struct TSlopeGeometry {
    REAL H = 5.;          ///< altura (m)
    REAL betaDeg = 45.;   ///< inclinação (graus)
    REAL Lc = 10.;        ///< largura da crista (m)
    REAL Lt = 10.;        ///< largura à direita do pé (m)
    REAL Hb = 5.;         ///< espessura abaixo do pé (m)
    REAL h = 0.5;         ///< tamanho alvo dos elementos (m)
    bool triangles = false;

    REAL Ls() const;               ///< projeção horizontal da face
    REAL D() const { return Hb + H; }
    REAL W() const { return Lc + Ls() + Lt; }
    /// profundidade abaixo da crista, y' = D - y (eixo Oy do artigo, para baixo, origem na crista)
    REAL Depth(REAL y) const { return D() - y; }

    /// Monta o TPZGeoMesh (elementos de material ESoil e de contorno ESlopeIds)
    TPZGeoMesh *CreateGeoMesh() const;

    std::string Describe() const;

    /// Geometrias do artigo: Cho (2010) coesivo 2:1 e c-φ 1:1 (e o talude de referência com percolação)
    static TSlopeGeometry Cho2H1V(REAL h = 0.5);
    static TSlopeGeometry Cho1H1V(REAL h = 0.5);
};

#endif
