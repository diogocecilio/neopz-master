#include "SlopeFields.h"

#include <algorithm>
#include <cmath>

namespace slope {

namespace {
// Retangulo dos nos: {x0, y0, Dx, Dy}
std::array<REAL, 4> Box(const FlatMesh &m)
{
    REAL x0 = 1.e30, x1 = -1.e30, y0 = 1.e30, y1 = -1.e30;
    for (const auto &p : m.nodes) x0 = std::min(x0, p[0]), x1 = std::max(x1, p[0]), y0 = std::min(y0, p[1]), y1 = std::max(y1, p[1]);
    return {x0, y0, x1 - x0, y1 - y0};
}
klsep::KL2D MakeKL(const FieldStats &s, const FlatMesh &m)
{
    const auto b = Box(m);
    return klsep::KL2D(b[0], b[1], b[2], b[3], s.Lx, s.Ly, 1., s.M);
}
} // namespace

RandomFields::RandomFields(const FieldStats &s, const FlatMesh &m)
    : fStats(s), fKL(MakeKL(s, m)), fLc(s.c, s.covC), fLphi(s.phiDeg, s.covPhi), fLk(1., s.covK)
{
    fColX.assign(fKL.mx.lambda.size(), -1);
    fColY.assign(fKL.my.lambda.size(), -1);
    for (int k = 0; k < fKL.M(); k++) {
        if (fColX[fKL.ix[k]] < 0) fColX[fKL.ix[k]] = fUsedX.size(), fUsedX.push_back(fKL.ix[k]);
        if (fColY[fKL.iy[k]] < 0) fColY[fKL.iy[k]] = fUsedY.size(), fUsedY.push_back(fKL.iy[k]);
    }
    for (auto &xi : fXi) xi.assign(fKL.M(), 0.); // sem Sample: Z = 0, a mediana de cada campo
}

void RandomFields::Sample(uint64_t seed, int64_t realization)
{
    for (int campo = 0; campo < 3; campo++)
        fXi[campo] = klsep::Normal(seed * 0x9E3779B97F4A7C15ULL + 3ULL * realization + campo).Vetor(fKL.M());
}

FieldValues RandomFields::At(REAL x, REAL y) const
{
    // modos 1D no ponto (medidos do centro do retangulo) e Z = sum_k sqrt(lambda_k) xi_k f_k(x) g_k(y) dos tres campos
    std::vector<REAL> f(fUsedX.size()), g(fUsedY.size());
    for (size_t i = 0; i < fUsedX.size(); i++) f[i] = fKL.mx.Avalia(fUsedX[i], x - (fKL.x0 + 0.5 * fKL.Dx));
    for (size_t j = 0; j < fUsedY.size(); j++) g[j] = fKL.my.Avalia(fUsedY[j], y - (fKL.y0 + 0.5 * fKL.Dy));
    REAL z[3] = {0., 0., 0.};
    for (int k = 0; k < fKL.M(); k++) {
        const REAL w = fKL.raiz[k] * f[fColX[fKL.ix[k]]] * g[fColY[fKL.iy[k]]];
        for (int campo = 0; campo < 3; campo++) z[campo] += w * fXi[campo][k];
    }
    return {fLc(z[0]), fLphi(z[1]) * M_PI / 180., fLk(z[2])};
}

} // namespace slope
