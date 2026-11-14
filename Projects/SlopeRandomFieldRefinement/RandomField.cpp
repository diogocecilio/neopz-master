#include "RandomField.h"

RandomField::RandomField(double x0, double y0, double Dx, double Dy, double Lx, double Ly, int M,
                         const std::vector<Stats> &campos)
    : fKL(x0, y0, Dx, Dy, Lx, Ly, 1., M)
{
    for (const Stats &s : campos) fLog.emplace_back(s.media, s.cov);
    // modos 1D que aparecem em algum termo: so eles sao avaliados em cada ponto
    fColX.assign(fKL.mx.lambda.size(), -1);
    fColY.assign(fKL.my.lambda.size(), -1);
    for (int k = 0; k < fKL.M(); k++) {
        if (fColX[fKL.ix[k]] < 0) fColX[fKL.ix[k]] = fUsadosX.size(), fUsadosX.push_back(fKL.ix[k]);
        if (fColY[fKL.iy[k]] < 0) fColY[fKL.iy[k]] = fUsadosY.size(), fUsadosY.push_back(fKL.iy[k]);
    }
    fXi.assign(campos.size(), std::vector<double>(fKL.M(), 0.)); // sem Sample: Z = 0, a mediana de cada campo
}

void RandomField::Sample(uint64_t semente, int64_t realizacao)
{
    const uint64_t n = fXi.size();
    for (uint64_t campo = 0; campo < n; campo++)
        fXi[campo] = klsep::Normal(semente * 0x9E3779B97F4A7C15ULL + n * realizacao + campo).Vetor(fKL.M());
}

std::vector<double> RandomField::At(double x, double y) const
{
    // modos 1D no ponto (medidos a partir do centro do retangulo) e Z = sum_k sqrt(lambda_k) xi_k f_k(x) g_k(y) de cada
    // campo; depois a transformacao lognormal
    std::vector<double> f(fUsadosX.size()), g(fUsadosY.size());
    for (size_t i = 0; i < fUsadosX.size(); i++) f[i] = fKL.mx.Avalia(fUsadosX[i], x - (fKL.x0 + 0.5 * fKL.Dx));
    for (size_t j = 0; j < fUsadosY.size(); j++) g[j] = fKL.my.Avalia(fUsadosY[j], y - (fKL.y0 + 0.5 * fKL.Dy));
    std::vector<double> z(fXi.size(), 0.);
    for (int k = 0; k < fKL.M(); k++) {
        const double w = fKL.raiz[k] * f[fColX[fKL.ix[k]]] * g[fColY[fKL.iy[k]]];
        for (size_t c = 0; c < fXi.size(); c++) z[c] += w * fXi[c][k];
    }
    for (size_t c = 0; c < fXi.size(); c++) z[c] = fLog[c](z[c]);
    return z;
}
