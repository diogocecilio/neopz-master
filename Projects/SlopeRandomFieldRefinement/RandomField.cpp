#include "RandomField.h"

#include <algorithm>
#include <cmath>

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
    for (uint64_t campo = 0; campo < n; campo++) {
        fXi[campo] = klsep::Normal(semente * 0x9E3779B97F4A7C15ULL + n * realizacao + campo).Vetor(fKL.M());
        if (!fMu.empty())
            for (int k = 0; k < fKL.M(); k++) fXi[campo][k] += fMu[campo][k];
    }
}

std::vector<double> RandomField::Direction(const std::vector<std::array<double, 3>> &pontos) const
{
    std::vector<double> a(fKL.M(), 0.), f(fUsadosX.size()), g(fUsadosY.size());
    double soma = 0.;
    for (const auto &p : pontos) {
        for (size_t i = 0; i < fUsadosX.size(); i++) f[i] = fKL.mx.Avalia(fUsadosX[i], p[0] - (fKL.x0 + 0.5 * fKL.Dx));
        for (size_t j = 0; j < fUsadosY.size(); j++) g[j] = fKL.my.Avalia(fUsadosY[j], p[1] - (fKL.y0 + 0.5 * fKL.Dy));
        for (int k = 0; k < fKL.M(); k++) a[k] += p[2] * fKL.raiz[k] * f[fColX[fKL.ix[k]]] * g[fColY[fKL.iy[k]]];
        soma += p[2];
    }
    double n2 = 0.;
    for (double &v : a) v /= soma, n2 += v * v;
    if (!(n2 > 0.)) return {};
    for (double &v : a) v /= std::sqrt(n2);
    return a;
}

void RandomField::SetShift(const std::vector<double> &d, const std::vector<double> &s)
{
    fMu.clear();
    if (s.empty()) return;
    fMu.assign(fXi.size(), std::vector<double>(fKL.M(), 0.));
    for (size_t c = 0; c < fXi.size() && c < s.size(); c++)
        for (int k = 0; k < fKL.M(); k++) fMu[c][k] = s[c] * d[k];
}

void RandomField::SetShift(const std::vector<std::vector<double>> &mu)
{
    fMu.assign(fXi.size(), std::vector<double>(fKL.M(), 0.));
    for (size_t c = 0; c < fXi.size() && c < mu.size(); c++)
        for (int k = 0; k < fKL.M() && k < (int)mu[c].size(); k++) fMu[c][k] = mu[c][k];
}

void RandomField::SetXi(const std::vector<std::vector<double>> &xi)
{
    for (size_t c = 0; c < fXi.size(); c++) {
        std::fill(fXi[c].begin(), fXi[c].end(), 0.);
        if (c < xi.size())
            for (int k = 0; k < fKL.M() && k < (int)xi[c].size(); k++) fXi[c][k] = xi[c][k];
    }
}

double RandomField::LogWeight() const
{
    double lw = 0.;
    for (size_t c = 0; c < fMu.size(); c++)
        for (int k = 0; k < fKL.M(); k++) lw += -fXi[c][k] * fMu[c][k] + 0.5 * fMu[c][k] * fMu[c][k];
    return lw;
}

std::vector<double> RandomField::Projection(const std::vector<double> &d) const
{
    std::vector<double> p(fXi.size(), 0.);
    if (d.empty()) return p;
    for (size_t c = 0; c < fXi.size(); c++)
        for (int k = 0; k < fKL.M(); k++) p[c] += fXi[c][k] * d[k];
    return p;
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
    if (fCorrB >= 0) z[fCorrB] = fRho * z[fCorrA] + std::sqrt(1. - fRho * fRho) * z[fCorrB];
    for (size_t c = 0; c < fXi.size(); c++) z[c] = fLog[c](z[c]);
    return z;
}

void RandomField::SetCrossCorrelation(int a, int b, double rho)
{
    if (a < 0 || b < 0 || a == b || a >= NCampos() || b >= NCampos() || !(std::abs(rho) < 1.)) {
        fCorrA = fCorrB = -1, fRho = 0.;
        return;
    }
    fCorrA = a, fCorrB = b, fRho = rho;
}
