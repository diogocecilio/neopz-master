#include "RandomField.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <map>
#include <sstream>
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
            fGQ.clear(); // as somas da grade sao da realizacao anterior
        }
        std::vector<double> z(fXi.size(), 0.);
        if (fGrade) { // as somas acumuladas na grade: O(1) por ponto
            if (fGQ.size() != fXi.size()) PrepareGrid();
            for (size_t c = 0; c < fXi.size(); c++) z[c] = GridValue(int(c), x, y);
        } else
            for (size_t j = 0; j < n; j++) {
                const double cxy = Kernel(x - fNysNos[j][0], y - fNysNos[j][1]);
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
            const double cw = fNysNos[j][2] * Kernel(x - fNysNos[j][0], y - fNysNos[j][1]);
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

double RandomField::Kernel(double dx, double dy) const
{
    if (fNucleo == EGauss) { // eixos girados de fAngulo
        const double u = dx * std::cos(fAngulo) + dy * std::sin(fAngulo), v = -dx * std::sin(fAngulo) + dy * std::cos(fAngulo);
        return std::exp(-(u / fLx) * (u / fLx) - (v / fLy) * (v / fLy));
    }
    return std::exp(-std::abs(dx) / fLx - std::abs(dy) / fLy);
}

std::string RandomField::KernelDescription() const
{
    std::ostringstream s;
    if (fNucleo == EGauss)
        s << "gauss lh=" << fLx << " lv=" << fLy << " angulo=" << fAngulo * 180. / M_PI;
    else
        s << "exp lh=" << fLx << " lv=" << fLy;
    return s.str();
}

void RandomField::UseNystrom(const std::string &arquivo, int M, bool rapido)
{
    std::ifstream f(arquivo);
    size_t n = 0;
    int mArq = 0;
    // 1a linha: n M e, nas bases novas, chave=valor do nucleo e da grade
    std::string linha1;
    std::getline(f, linha1);
    std::istringstream l1(linha1);
    if (!(l1 >> n >> mArq) || n == 0 || mArq <= 0) throw std::runtime_error("base de Nystrom ilegivel: " + arquivo);
    std::map<std::string, std::string> cab;
    for (std::string t; l1 >> t;) {
        const size_t eq = t.find('=');
        if (eq != std::string::npos) cab[t.substr(0, eq)] = t.substr(eq + 1);
    }
    if (cab.count("nucleo")) {
        if (cab["nucleo"] != "exp" && cab["nucleo"] != "gauss")
            throw std::runtime_error("base de Nystrom com nucleo desconhecido (" + cab["nucleo"] + "): " + arquivo);
        if (!cab.count("lh") || !cab.count("lv")) throw std::runtime_error("base de Nystrom sem lh e lv: " + arquivo);
        fNucleo = cab["nucleo"] == "gauss" ? EGauss : EExp;
        fLx = std::stod(cab["lh"]), fLy = std::stod(cab["lv"]);
        fAngulo = cab.count("angulo") ? std::stod(cab["angulo"]) * M_PI / 180. : 0.;
    }
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
    fGQ.clear();
    // a grade regular dos nos (so o nucleo exponencial, e so se pedida): o passo do cabecalho ou o menor intervalo entre
    // as coordenadas; cada no tem de cair num centro de celula (tolerancia de 1e-4 da celula: as bases antigas tem 6 casas)
    fGrade = false;
    if (!rapido || fNucleo != EExp) return;
    double xmin = 1e300, xmax = -1e300, ymin = 1e300, ymax = -1e300;
    for (const auto &p : fNysNos)
        xmin = std::min(xmin, p[0]), xmax = std::max(xmax, p[0]), ymin = std::min(ymin, p[1]), ymax = std::max(ymax, p[1]);
    auto passo = [&](int d, const char *chave) {
        if (cab.count(chave)) return std::stod(cab[chave]);
        std::vector<double> v;
        for (const auto &p : fNysNos) v.push_back(p[d]);
        std::sort(v.begin(), v.end());
        double h = 1e300;
        for (size_t i = 1; i < v.size(); i++)
            if (v[i] - v[i - 1] > 1e-9) h = std::min(h, v[i] - v[i - 1]);
        return h;
    };
    const double hx = passo(0, "dx"), hy = passo(1, "dy");
    if (!(hx > 0. && hx < 1e299 && hy > 0. && hy < 1e299)) return;
    const int nx = int(std::lround((xmax - xmin) / hx)) + 1, ny = int(std::lround((ymax - ymin) / hy)) + 1;
    std::vector<int64_t> celula(n);
    std::vector<char> usada(size_t(nx) * ny, 0);
    for (size_t j = 0; j < n; j++) {
        const long i = std::lround((fNysNos[j][0] - xmin) / hx), k = std::lround((fNysNos[j][1] - ymin) / hy);
        if (i < 0 || i >= nx || k < 0 || k >= ny || std::abs(fNysNos[j][0] - (xmin + i * hx)) > 1e-4 * hx ||
            std::abs(fNysNos[j][1] - (ymin + k * hy)) > 1e-4 * hy || usada[size_t(k) * nx + i])
            return; // nao e uma grade regular: fica a soma direta
        celula[j] = int64_t(k) * nx + i;
        usada[celula[j]] = 1;
    }
    fGrade = true;
    fGNx = nx, fGNy = ny, fGx0 = xmin, fGy0 = ymin, fGdx = hx, fGdy = hy;
    fGCelula = celula;
}

// Z(x) = sum_j exp(-|x - X_j|/Lx - |y - Y_j|/Ly) a_j com os nos nos centros da grade (X_i = fGx0 + i dx, Y_k = fGy0 + k dy):
// os nos abaixo de y somam e^{-(y - Y_kl)/Ly} Pf[kl] e os de cima e^{-(Y_kl+1 - y)/Ly} Pb[kl + 1], com
//   Pf[k] = A[k] + ry Pf[k - 1],  Pb[k] = A[k] + ry Pb[k + 1],  ry = e^{-dy/Ly}   (em cada coluna i),
// e o mesmo em x sobre Pf e Pb (rx = e^{-dx/Lx}): Q[0] = frente(Pf), Q[1] = tras(Pf), Q[2] = frente(Pb), Q[3] = tras(Pb)
void RandomField::PrepareGrid() const
{
    const int nx = fGNx, ny = fGNy;
    const size_t nc = size_t(nx) * ny;
    const double rx = std::exp(-fGdx / fLx), ry = std::exp(-fGdy / fLy);
    fGQ.assign(fXi.size(), {});
    std::vector<double> A(nc), Pf(nc), Pb(nc);
    for (size_t c = 0; c < fXi.size(); c++) {
        std::fill(A.begin(), A.end(), 0.);
        for (size_t j = 0; j < fGCelula.size(); j++) A[fGCelula[j]] = fNysA[c][j];
        for (int i = 0; i < nx; i++) {
            for (int k = 0; k < ny; k++) Pf[size_t(k) * nx + i] = A[size_t(k) * nx + i] + (k > 0 ? ry * Pf[size_t(k - 1) * nx + i] : 0.);
            for (int k = ny - 1; k >= 0; k--)
                Pb[size_t(k) * nx + i] = A[size_t(k) * nx + i] + (k < ny - 1 ? ry * Pb[size_t(k + 1) * nx + i] : 0.);
        }
        auto &Q = fGQ[c];
        for (auto &q : Q) q.assign(nc, 0.);
        for (int k = 0; k < ny; k++) {
            const size_t r = size_t(k) * nx;
            for (int i = 0; i < nx; i++) {
                Q[0][r + i] = Pf[r + i] + (i > 0 ? rx * Q[0][r + i - 1] : 0.);
                Q[2][r + i] = Pb[r + i] + (i > 0 ? rx * Q[2][r + i - 1] : 0.);
            }
            for (int i = nx - 1; i >= 0; i--) {
                Q[1][r + i] = Pf[r + i] + (i < nx - 1 ? rx * Q[1][r + i + 1] : 0.);
                Q[3][r + i] = Pb[r + i] + (i < nx - 1 ? rx * Q[3][r + i + 1] : 0.);
            }
        }
    }
}

double RandomField::GridValue(int c, double x, double y) const
{
    const int nx = fGNx, ny = fGNy;
    const int il = std::clamp(int(std::floor((x - fGx0) / fGdx)), -1, nx - 1); // ultima coluna com X <= x (-1: nenhuma)
    const int kl = std::clamp(int(std::floor((y - fGy0) / fGdy)), -1, ny - 1);
    const double exf = il >= 0 ? std::exp(-(x - (fGx0 + il * fGdx)) / fLx) : 0.;
    const double exb = il + 1 < nx ? std::exp(-((fGx0 + (il + 1) * fGdx) - x) / fLx) : 0.;
    const double eyf = kl >= 0 ? std::exp(-(y - (fGy0 + kl * fGdy)) / fLy) : 0.;
    const double eyb = kl + 1 < ny ? std::exp(-((fGy0 + (kl + 1) * fGdy) - y) / fLy) : 0.;
    const auto &Q = fGQ[c];
    double z = 0.;
    if (kl >= 0) {
        const size_t r = size_t(kl) * nx;
        if (il >= 0) z += eyf * exf * Q[0][r + il];
        if (il + 1 < nx) z += eyf * exb * Q[1][r + il + 1];
    }
    if (kl + 1 < ny) {
        const size_t r = size_t(kl + 1) * nx;
        if (il >= 0) z += eyb * exf * Q[2][r + il];
        if (il + 1 < nx) z += eyb * exb * Q[3][r + il + 1];
    }
    return z;
}

std::vector<std::vector<double>> RandomField::ProjectModes(const std::vector<std::array<double, 2>> &pontos,
                                                           const std::vector<std::vector<double>> &pesos, int K) const
{
    K = std::min(K, fM);
    std::vector<std::vector<double>> out(pesos.size(), std::vector<double>(K, 0.));
    if (!Nystrom()) { // analitica: ponto a ponto
        for (size_t p = 0; p < pontos.size(); p++) {
            const std::vector<double> m = ModeValues(pontos[p][0], pontos[p][1], K);
            for (size_t a = 0; a < pesos.size(); a++)
                for (int k = 0; k < K; k++) out[a][k] += pesos[a][p] * m[k];
        }
        return out;
    }
    const size_t n = fNysNos.size();
    std::vector<std::vector<double>> b(pesos.size(), std::vector<double>(n, 0.)); // b_j = sum_p a_p C(x_p, x_j)
    for (size_t p = 0; p < pontos.size(); p++)
        for (size_t j = 0; j < n; j++) {
            const double cxy = Kernel(pontos[p][0] - fNysNos[j][0], pontos[p][1] - fNysNos[j][1]);
            for (size_t a = 0; a < pesos.size(); a++) b[a][j] += pesos[a][p] * cxy;
        }
    for (size_t a = 0; a < pesos.size(); a++) {
        for (size_t j = 0; j < n; j++) {
            const double cw = fNysNos[j][2] * b[a][j];
            for (int k = 0; k < K; k++) out[a][k] += cw * fNysF[j][k];
        }
        for (int k = 0; k < K; k++) out[a][k] /= std::sqrt(fNysLam[k]);
    }
    return out;
}
