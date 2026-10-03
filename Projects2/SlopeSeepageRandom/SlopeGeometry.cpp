// SlopeGeometry.cpp — ver SlopeGeometry.h

#include "SlopeGeometry.h"

#include <array>
#include <cmath>
#include <sstream>
#include <vector>

#include "pzgeoel.h"
#include "pzgnode.h"

REAL TSlopeGeometry::Ls() const {
    if (betaDeg >= 90. - 1.e-9) return 0.;
    return H / std::tan(betaDeg * M_PI / 180.);
}

TSlopeGeometry TSlopeGeometry::Cho2H1V(REAL h) {
    TSlopeGeometry g;
    g.H = 5.;
    g.betaDeg = std::atan(0.5) * 180. / M_PI;  // 2:1 (≈ 26.57°)
    g.Lc = 10.;
    g.Lt = 10.;
    g.Hb = 5.;
    g.h = h;
    return g;
}

TSlopeGeometry TSlopeGeometry::Cho1H1V(REAL h) {
    TSlopeGeometry g;
    g.H = 5.;
    g.betaDeg = 45.;
    g.Lc = 10.;
    g.Lt = 10.;
    g.Hb = 5.;
    g.h = h;
    return g;
}

std::string TSlopeGeometry::Describe() const {
    std::ostringstream s;
    s << "H = " << H << " m, beta = " << betaDeg << " graus, crista " << Lc << " m, pe " << Lt << " m, base " << Hb
      << " m abaixo do pe (dominio " << W() << " x " << D() << " m), h = " << h << " m, "
      << (triangles ? "triangulos" : "quadrilateros");
    return s.str();
}

TPZGeoMesh *TSlopeGeometry::CreateGeoMesh() const {
    const REAL ls = Ls();
    auto ndiv = [this](REAL len) { return std::max(1, (int)std::lround(len / h)); };
    const int nxB = ndiv(Lc + ls);  // colunas do bloco B (e da parte esquerda do bloco A)
    const int nxT = ndiv(Lt);       // colunas à direita do pé
    const int nyA = ndiv(Hb);
    const int nyB = ndiv(H);
    const int nxA = nxB + nxT;

    auto *gmesh = new TPZGeoMesh;
    gmesh->SetDimension(2);
    std::vector<std::array<REAL, 2>> xy;
    auto addNode = [&](REAL x, REAL y) {
        xy.push_back({x, y});
        return (int64_t)xy.size() - 1;
    };
    // bloco A: (nxA+1) x (nyA+1)
    std::vector<std::vector<int64_t>> nA(nxA + 1, std::vector<int64_t>(nyA + 1));
    for (int j = 0; j <= nyA; j++) {
        const REAL y = Hb * j / nyA;
        for (int i = 0; i <= nxA; i++) {
            const REAL x = (i <= nxB) ? (Lc + ls) * i / nxB : (Lc + ls) + Lt * (i - nxB) / nxT;
            nA[i][j] = addNode(x, y);
        }
    }
    // bloco B: (nxB+1) x (nyB+1), linha j = 0 compartilhada com o topo do bloco A
    std::vector<std::vector<int64_t>> nB(nxB + 1, std::vector<int64_t>(nyB + 1));
    for (int i = 0; i <= nxB; i++) nB[i][0] = nA[i][nyA];
    for (int j = 1; j <= nyB; j++) {
        const REAL t = REAL(j) / nyB;
        const REAL xr = Lc + ls * (1. - t);
        for (int i = 0; i <= nxB; i++) nB[i][j] = addNode(xr * i / nxB, Hb + t * H);
    }
    gmesh->NodeVec().Resize(xy.size());
    for (size_t k = 0; k < xy.size(); k++) {
        TPZManVector<REAL, 3> co = {xy[k][0], xy[k][1], 0.};
        gmesh->NodeVec()[k].Initialize(co, *gmesh);
    }

    auto quad = [&](int64_t a, int64_t b, int64_t c, int64_t dd) {
        int64_t index;
        if (!triangles) {
            TPZManVector<int64_t, 4> top = {a, b, c, dd};
            gmesh->CreateGeoElement(EQuadrilateral, top, ESoil, index);
        } else {
            // diagonal alternada para não privilegiar uma direção
            TPZManVector<int64_t, 3> t1 = {a, b, c}, t2 = {a, c, dd};
            gmesh->CreateGeoElement(ETriangle, t1, ESoil, index);
            gmesh->CreateGeoElement(ETriangle, t2, ESoil, index);
        }
    };
    for (int j = 0; j < nyA; j++)
        for (int i = 0; i < nxA; i++) quad(nA[i][j], nA[i + 1][j], nA[i + 1][j + 1], nA[i][j + 1]);
    for (int j = 0; j < nyB; j++)
        for (int i = 0; i < nxB; i++) quad(nB[i][j], nB[i + 1][j], nB[i + 1][j + 1], nB[i][j + 1]);

    auto line = [&](int64_t a, int64_t b, int id) {
        int64_t index;
        TPZManVector<int64_t, 2> top = {a, b};
        gmesh->CreateGeoElement(EOned, top, id, index);
    };
    for (int i = 0; i < nxA; i++) line(nA[i][0], nA[i + 1][0], EBottom);
    for (int j = 0; j < nyA; j++) line(nA[nxA][j], nA[nxA][j + 1], ERight);
    for (int i = nxA; i > nxB; i--) line(nA[i][nyA], nA[i - 1][nyA], EToe);
    for (int j = 0; j < nyB; j++) line(nB[nxB][j], nB[nxB][j + 1], EFace);
    for (int i = nxB; i > 0; i--) line(nB[i][nyB], nB[i - 1][nyB], ECrest);
    for (int j = nyB; j > 0; j--) line(nB[0][j], nB[0][j - 1], ELeft);
    for (int j = nyA; j > 0; j--) line(nA[0][j], nA[0][j - 1], ELeft);
    gmesh->BuildConnectivity();
    return gmesh;
}
