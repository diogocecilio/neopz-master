#include "RandomField.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <stdexcept>

RandomField::RandomField(double x0, double y0, double Dx, double Dy, double Lx, double Ly, int M,
                         const std::vector<Stats> &campos)
    : fKL(x0, y0, Dx, Dy, Lx, Ly, 1., M), fM(M), fLx(Lx), fLy(Ly)
{
    fM = fKL.M();
    for (const Stats &s : campos) fLog.emplace_back(s.media, s.cov);
    // modos 1D que aparecem em algum termo: so eles sao avaliados em cada ponto
    fColX.assign(fKL.mx.lambda.size(), -1);
    fColY.assign(fKL.my.lambda.size(), -1);
    for (int k = 0; k < fM; k++) {
        if (fColX[fKL.ix[k]] < 0) fColX[fKL.ix[k]] = fUsadosX.size(), fUsadosX.push_back(fKL.ix[k]);
        if (fColY[fKL.iy[k]] < 0) fColY[fKL.iy[k]] = fUsadosY.size(), fUsadosY.push_back(fKL.iy[k]);
    }
    fXi.assign(campos.size(), std::vector<double>(fM, 0.)); // sem Sample: Z = 0, a mediana de cada campo
}

void RandomField::Sample(uint64_t semente, int64_t realizacao)
{
    const uint64_t n = fXi.size();
    for (uint64_t campo = 0; campo < n; campo++) {
        fXi[campo] = klsep::Normal(semente * 0x9E3779B97F4A7C15ULL + n * realizacao + campo).Vetor(fM);
        if (!fMu.empty())
            for (int k = 0; k < fM; k++) fXi[campo][k] += fMu[campo][k];
    }
    fNysSujo = true;
}

std::vector<double> RandomField::Direction(const std::vector<std::array<double, 3>> &pontos) const
{
    std::vector<double> a(fM, 0.);
    double soma = 0.;
    for (const auto &p : pontos) {
        const std::vector<double> m = ModeValues(p[0], p[1], fM);
        for (int k = 0; k < fM; k++) a[k] += p[2] * m[k];
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
    fMu.assign(fXi.size(), std::vector<double>(fM, 0.));
    for (size_t c = 0; c < fXi.size() && c < s.size(); c++)
        for (int k = 0; k < fM; k++) fMu[c][k] = s[c] * d[k];
}

void RandomField::SetShift(const std::vector<std::vector<double>> &mu)
{
    fMu.assign(fXi.size(), std::vector<double>(fM, 0.));
    for (size_t c = 0; c < fXi.size() && c < mu.size(); c++)
        for (int k = 0; k < fM && k < (int)mu[c].size(); k++) fMu[c][k] = mu[c][k];
}

void RandomField::SetXi(const std::vector<std::vector<double>> &xi)
{
    for (size_t c = 0; c < fXi.size(); c++) {
        std::fill(fXi[c].begin(), fXi[c].end(), 0.);
        if (c < xi.size())
            for (int k = 0; k < fM && k < (int)xi[c].size(); k++) fXi[c][k] = xi[c][k];
    }
    fNysSujo = true;
}

double RandomField::LogWeight() const
{
    double lw = 0.;
    for (size_t c = 0; c < fMu.size(); c++)
        for (int k = 0; k < fM; k++) lw += -fXi[c][k] * fMu[c][k] + 0.5 * fMu[c][k] * fMu[c][k];
    return lw;
}

std::vector<double> RandomField::Projection(const std::vector<double> &d) const
{
    std::vector<double> p(fXi.size(), 0.);
    if (d.empty()) return p;
    for (size_t c = 0; c < fXi.size(); c++)
        for (int k = 0; k < fM; k++) p[c] += fXi[c][k] * d[k];
    return p;
}

std::vector<double> RandomField::At(double x, double y) const
{
    if (Nystrom()) {
        const size_t n = fNysNos.size();
        if (fNysSujo) {
            fNysA.assign(fXi.size(), std::vector<double>(n, 0.));
            for (size_t c = 0; c < fXi.size(); c++)
                for (size_t j = 0; j < n; j++) {
                    double v = 0.;
                    for (int k = 0; k < fM; k++) v += fXi[c][k] * fNysF[j][k] / std::sqrt(fNysLam[k]);
                    fNysA[c][j] = fNysNos[j][2] * v;
                }
            fNysSujo = false;
        }
        std::vector<double> z(fXi.size(), 0.);
        for (size_t j = 0; j < n; j++) {
            const double cxy = std::exp(-std::abs(x - fNysNos[j][0]) / fLx - std::abs(y - fNysNos[j][1]) / fLy);
            for (size_t c = 0; c < fXi.size(); c++) z[c] += cxy * fNysA[c][j];
        }
        if (fCorrB >= 0) z[fCorrB] = fRho * z[fCorrA] + std::sqrt(1. - fRho * fRho) * z[fCorrB];
        for (size_t c = 0; c < fXi.size(); c++) z[c] = fLog[c](z[c]);
        return z;
    }
    // modos 1D no ponto (medidos a partir do centro do retangulo) e Z = sum_k sqrt(lambda_k) xi_k f_k(x) g_k(y) de cada
    // campo; depois a transformacao lognormal
    std::vector<double> f(fUsadosX.size()), g(fUsadosY.size());
    for (size_t i = 0; i < fUsadosX.size(); i++) f[i] = fKL.mx.Avalia(fUsadosX[i], x - (fKL.x0 + 0.5 * fKL.Dx));
    for (size_t j = 0; j < fUsadosY.size(); j++) g[j] = fKL.my.Avalia(fUsadosY[j], y - (fKL.y0 + 0.5 * fKL.Dy));
    std::vector<double> z(fXi.size(), 0.);
    for (int k = 0; k < fM; k++) {
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

std::vector<double> RandomField::ModeValues(double x, double y, int K) const
{
    K = std::min(K, fM);
    std::vector<double> m(K, 0.);
    if (Nystrom()) {
        for (size_t j = 0; j < fNysNos.size(); j++) {
            const double cw =
                fNysNos[j][2] * std::exp(-std::abs(x - fNysNos[j][0]) / fLx - std::abs(y - fNysNos[j][1]) / fLy);
            for (int k = 0; k < K; k++) m[k] += cw * fNysF[j][k];
        }
        for (int k = 0; k < K; k++) m[k] /= std::sqrt(fNysLam[k]); // sqrt(lambda) f(x) = sum w C f_j / sqrt(lambda)
        return m;
    }
    std::vector<double> f(fUsadosX.size()), g(fUsadosY.size());
    for (size_t i = 0; i < fUsadosX.size(); i++) f[i] = fKL.mx.Avalia(fUsadosX[i], x - (fKL.x0 + 0.5 * fKL.Dx));
    for (size_t j = 0; j < fUsadosY.size(); j++) g[j] = fKL.my.Avalia(fUsadosY[j], y - (fKL.y0 + 0.5 * fKL.Dy));
    for (int k = 0; k < K; k++) m[k] = fKL.raiz[k] * f[fColX[fKL.ix[k]]] * g[fColY[fKL.iy[k]]];
    return m;
}

void RandomField::UseNystrom(const std::string &arquivo, int M)
{
    std::ifstream f(arquivo);
    size_t n = 0;
    int mArq = 0;
    if (!(f >> n >> mArq) || n == 0 || mArq <= 0) throw std::runtime_error("base de Nystrom ilegivel: " + arquivo);
    fNysNos.assign(n, {0., 0., 0.});
    double area = 0.;
    for (auto &p : fNysNos) f >> p[0] >> p[1] >> p[2], area += p[2];
    std::vector<double> lam(mArq);
    for (double &l : lam) f >> l;
    fM = std::min(M, mArq);
    fNysLam.assign(lam.begin(), lam.begin() + fM);
    fNysF.assign(n, std::vector<double>(fM, 0.));
    std::vector<double> linha(mArq);
    for (size_t j = 0; j < n; j++) {
        for (double &v : linha) f >> v;
        std::copy(linha.begin(), linha.begin() + fM, fNysF[j].begin());
    }
    if (!f) throw std::runtime_error("base de Nystrom incompleta: " + arquivo);
    double soma = 0.;
    for (double l : fNysLam) soma += l;
    fNysFrac = soma / area;
    for (auto &x : fXi) x.assign(fM, 0.);
    fMu.clear();
    fNysSujo = true;
}
