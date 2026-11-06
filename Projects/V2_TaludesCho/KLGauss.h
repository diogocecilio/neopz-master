// KLGauss.h - expansao de Karhunen-Loeve separavel em 2D com autopares 1D numericos (metodo de Nystrom), para o nucleo
// gaussiano (exponencial quadratico, rho2 de Liao et al. 2025, eq. 1)
//     C(x1, x2) = exp(-((x1 - x2) / hx)^2 - ((z1 - z2) / hz)^2) = exp(-(dx / hx)^2) exp(-(dz / hz)^2)
// e para o exponencial separavel (rho1), cujos autopares 1D tem forma fechada (Comum/KLSeparavel.h). Os dois nucleos sao
// produtos de nucleos 1D, entao os autopares 2D sao os produtos dos autopares 1D (como em KLSeparavel.h).
//
// Nystrom (Atkinson 1997; Betz, Papaioannou e Straub 2014): a equacao integral de Fredholm de segunda especie
//     int_{-a}^{a} k(s - t) f(t) dt = lambda f(s)
// discretizada com a quadratura de Gauss-Legendre (pontos t_j, pesos w_j) vira o problema simetrico
//     B v = lambda v,  B_ij = sqrt(w_i) k(t_i - t_j) sqrt(w_j),  f(t_j) = v_j / sqrt(w_j)
// (resolvido por Jacobi ciclico), com sum_j w_j f(t_j)^2 = 1; fora dos pontos a autofuncao vem da propria equacao
// integral (interpolacao de Nystrom): f(s) = (1 / lambda) sum_j w_j k(s - t_j) f(t_j). Para o nucleo gaussiano, analitico,
// a convergencia em n e exponencial. Sem dependencia do NeoPZ.
#pragma once
#include "KLSeparavel.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <memory>
#include <stdexcept>
#include <vector>

namespace v2kl {

/// Pontos e pesos de Gauss-Legendre em [-1, 1] (n pontos): Newton nas raizes do polinomio de Legendre P_n
inline void GaussLegendre(int n, std::vector<double> &x, std::vector<double> &w) {
    if (n < 1) throw std::invalid_argument("GaussLegendre: n < 1");
    x.assign(n, 0.);
    w.assign(n, 0.);
    auto legendre = [n](double z, double &pn, double &dpn) { // P_n(z) e P_n'(z)
        double p1 = 1., p2 = 0.;
        for (int j = 1; j <= n; ++j) {
            const double p3 = p2;
            p2 = p1;
            p1 = ((2. * j - 1.) * z * p2 - (j - 1.) * p3) / j;
        }
        pn = p1;
        dpn = n * (z * p1 - p2) / (z * z - 1.);
    };
    for (int i = 0; i < (n + 1) / 2; ++i) {
        double z = std::cos(M_PI * (i + 0.75) / (n + 0.5)), pn = 0., dpn = 1.;
        for (int it = 0; it < 100; ++it) {
            legendre(z, pn, dpn);
            const double dz = pn / dpn;
            z -= dz;
            if (std::fabs(dz) < 1.e-16) break;
        }
        legendre(z, pn, dpn);
        x[i] = -z;
        x[n - 1 - i] = z;
        w[i] = w[n - 1 - i] = 2. / ((1. - z * z) * dpn * dpn);
    }
}

/// Autovalores d e autovetores (colunas de V, n x n por linhas) da matriz simetrica A (n x n por linhas), por Jacobi
/// ciclico (rotacoes de Givens que zeram A_pq; Press et al., Numerical Recipes, sec. 11.1)
inline void JacobiSimetrico(std::vector<double> A, int n, std::vector<double> &d, std::vector<double> &V) {
    V.assign((size_t)n * n, 0.);
    for (int i = 0; i < n; ++i) V[(size_t)i * n + i] = 1.;
    for (int varredura = 0; varredura < 60; ++varredura) {
        double fora = 0., diag = 0.;
        for (int i = 0; i < n; ++i) {
            diag += A[(size_t)i * n + i] * A[(size_t)i * n + i];
            for (int j = i + 1; j < n; ++j) fora += A[(size_t)i * n + j] * A[(size_t)i * n + j];
        }
        if (fora <= 1.e-30 * diag) break;
        for (int p = 0; p < n - 1; ++p)
            for (int q = p + 1; q < n; ++q) {
                const double apq = A[(size_t)p * n + q];
                if (std::fabs(apq) < 1.e-300) continue;
                const double theta = (A[(size_t)q * n + q] - A[(size_t)p * n + p]) / (2. * apq);
                const double t = (theta >= 0. ? 1. : -1.) / (std::fabs(theta) + std::sqrt(theta * theta + 1.));
                const double c = 1. / std::sqrt(t * t + 1.), s = t * c;
                for (int k = 0; k < n; ++k) { // A P
                    const double akp = A[(size_t)k * n + p], akq = A[(size_t)k * n + q];
                    A[(size_t)k * n + p] = c * akp - s * akq;
                    A[(size_t)k * n + q] = s * akp + c * akq;
                }
                for (int k = 0; k < n; ++k) { // P^T (A P)
                    const double apk = A[(size_t)p * n + k], aqk = A[(size_t)q * n + k];
                    A[(size_t)p * n + k] = c * apk - s * aqk;
                    A[(size_t)q * n + k] = s * apk + c * aqk;
                }
                for (int k = 0; k < n; ++k) { // V P
                    const double vkp = V[(size_t)k * n + p], vkq = V[(size_t)k * n + q];
                    V[(size_t)k * n + p] = c * vkp - s * vkq;
                    V[(size_t)k * n + q] = s * vkp + c * vkq;
                }
            }
    }
    d.resize(n);
    for (int i = 0; i < n; ++i) d[i] = A[(size_t)i * n + i];
}

/// Autopares 1D de um nucleo estacionario num intervalo de comprimento D centrado na origem ([-D/2, D/2]), em ordem
/// decrescente de lambda; f(i, s): a autofuncao i (norma L2 unitaria no intervalo), s medido a partir do centro
struct Base1D {
    double comprimento = 0.;
    std::vector<double> lambda;
    std::function<double(int, double)> f;
    int N() const { return (int)lambda.size(); }
};

/// Base do nucleo exponencial exp(-|s - t| / L): forma fechada de Comum/KLSeparavel.h (n autopares)
inline Base1D BaseExponencial(double comprimento, double L, int n) {
    auto m = std::make_shared<klsep::Modos1D>(klsep::CalculaModos1D(comprimento, L, n));
    Base1D b;
    b.comprimento = comprimento;
    b.lambda = m->lambda;
    b.f = [m](int i, double s) { return m->Avalia(i, s); };
    return b;
}

/// Base de Nystrom do nucleo k(r) (k(0) = 1) com nq pontos de Gauss-Legendre; guarda os modos com
/// lambda > limiar * lambda_1 (os outros sao ruido de arredondamento: a interpolacao divide por lambda)
inline Base1D BaseNystrom(double comprimento, std::function<double(double)> k, int nq, double limiar = 1.e-11) {
    if (comprimento <= 0. || nq < 2) throw std::invalid_argument("BaseNystrom: argumentos invalidos");
    const double a = 0.5 * comprimento;
    std::vector<double> t, w;
    GaussLegendre(nq, t, w);
    for (int j = 0; j < nq; ++j) { t[j] *= a; w[j] *= a; }
    std::vector<double> B((size_t)nq * nq);
    for (int i = 0; i < nq; ++i)
        for (int j = 0; j < nq; ++j) B[(size_t)i * nq + j] = std::sqrt(w[i] * w[j]) * k(t[i] - t[j]);
    std::vector<double> d, V;
    JacobiSimetrico(B, nq, d, V);
    std::vector<int> ordem(nq);
    for (int i = 0; i < nq; ++i) ordem[i] = i;
    std::sort(ordem.begin(), ordem.end(), [&](int p, int q) { return d[p] > d[q]; });
    // g[m][j] = w_j f_m(t_j) / lambda_m: a interpolacao de Nystrom e f_m(s) = sum_j k(s - t_j) g[m][j]
    auto dados = std::make_shared<std::pair<std::vector<double>, std::vector<std::vector<double>>>>();
    dados->first = t;
    Base1D b;
    b.comprimento = comprimento;
    for (int m = 0; m < nq; ++m) {
        const int col = ordem[m];
        const double lam = d[col];
        if (!(lam > limiar * d[ordem[0]])) break;
        std::vector<double> g(nq);
        for (int j = 0; j < nq; ++j) g[j] = w[j] * (V[(size_t)j * nq + col] / std::sqrt(w[j])) / lam;
        b.lambda.push_back(lam);
        dados->second.push_back(std::move(g));
    }
    b.f = [dados, k](int i, double s) {
        const auto &tt = dados->first;
        const auto &g = dados->second[i];
        double v = 0.;
        for (size_t j = 0; j < tt.size(); ++j) v += k(s - tt[j]) * g[j];
        return v;
    };
    return b;
}

/// Base do nucleo gaussiano exp(-((s - t) / h)^2) (rho2 de Liao et al. 2025); nq = 0: escolhe pelo comprimento
inline Base1D BaseGaussiana(double comprimento, double h, int nq = 0) {
    if (h <= 0.) throw std::invalid_argument("BaseGaussiana: h <= 0");
    if (nq <= 0) nq = std::min(400, std::max(80, (int)std::ceil(12. * comprimento / h)));
    return BaseNystrom(comprimento, [h](double r) { return std::exp(-(r / h) * (r / h)); }, nq);
}

/// Expansao 2D separavel truncada: os M maiores produtos lambda_i^x lambda_j^y, ate a fracao de variancia pedida
/// (media no retangulo: sum lambda_k / (Dx Dy)) ou Mmax termos, o que vier primeiro. Mesma escolha de KL2D de
/// KLSeparavel.h, com bases 1D quaisquer.
struct KLSep2D {
    double x0 = 0., y0 = 0., Dx = 0., Dy = 0.;
    Base1D bx, by;
    std::vector<int> ix, iy;  ///< modo 1D em x e em y do termo k
    std::vector<double> raiz; ///< sqrt(lambda_k)
    double fracao = 0.;       ///< fracao de variancia retida

    int M() const { return (int)raiz.size(); }
    double Lambda(int k) const { return raiz[k] * raiz[k]; }

    KLSep2D(double x0_, double y0_, double Dx_, double Dy_, Base1D bx_, Base1D by_, double fracaoAlvo, int Mmax)
        : x0(x0_), y0(y0_), Dx(Dx_), Dy(Dy_), bx(std::move(bx_)), by(std::move(by_)) {
        struct T { double lam; int i, j; };
        std::vector<T> cand;
        for (int i = 0; i < bx.N(); ++i)
            for (int j = 0; j < by.N() && (long)(i + 1) * (j + 1) <= Mmax; ++j)
                cand.push_back({bx.lambda[i] * by.lambda[j], i, j});
        std::sort(cand.begin(), cand.end(), [](const T &p, const T &q) { return p.lam > q.lam; });
        const double total = Dx * Dy;
        double soma = 0.;
        for (const auto &t : cand) {
            if ((int)raiz.size() >= Mmax || soma / total >= fracaoAlvo) break;
            ix.push_back(t.i);
            iy.push_back(t.j);
            raiz.push_back(std::sqrt(t.lam));
            soma += t.lam;
        }
        fracao = soma / total;
    }

    /// Valores dos modos 1D usados nos pontos fixos (os pontos de integracao)
    struct Pontos {
        int np = 0;
        std::vector<int> mapx, mapy;               ///< modo 1D -> coluna (-1 se nao usado)
        std::vector<std::vector<double>> fx, gy;   ///< [coluna][ponto]
        std::vector<double> var;                   ///< variancia truncada em cada ponto
    };

    Pontos PreparaPontos(const std::vector<double> &x, const std::vector<double> &y) const {
        Pontos P;
        P.np = (int)x.size();
        P.mapx.assign(bx.N(), -1);
        P.mapy.assign(by.N(), -1);
        for (int k = 0; k < M(); ++k) {
            if (P.mapx[ix[k]] < 0) { P.mapx[ix[k]] = (int)P.fx.size(); P.fx.emplace_back(P.np); }
            if (P.mapy[iy[k]] < 0) { P.mapy[iy[k]] = (int)P.gy.size(); P.gy.emplace_back(P.np); }
        }
        const double xc = x0 + 0.5 * Dx, yc = y0 + 0.5 * Dy;
        for (int i = 0; i < bx.N(); ++i)
            if (P.mapx[i] >= 0)
                for (int p = 0; p < P.np; ++p) P.fx[P.mapx[i]][p] = bx.f(i, x[p] - xc);
        for (int j = 0; j < by.N(); ++j)
            if (P.mapy[j] >= 0)
                for (int p = 0; p < P.np; ++p) P.gy[P.mapy[j]][p] = by.f(j, y[p] - yc);
        P.var.assign(P.np, 0.);
        for (int k = 0; k < M(); ++k) {
            const auto &f = P.fx[P.mapx[ix[k]]], &g = P.gy[P.mapy[iy[k]]];
            for (int p = 0; p < P.np; ++p) P.var[p] += Lambda(k) * f[p] * f[p] * g[p] * g[p];
        }
        return P;
    }

    /// Campo gaussiano padrao (truncado) nos pontos, para as normais padrao xi (M valores)
    std::vector<double> Campo(const Pontos &P, const std::vector<double> &xi) const {
        if ((int)xi.size() != M()) throw std::invalid_argument("KLSep2D::Campo: xi deve ter M valores");
        std::vector<double> z(P.np, 0.);
        for (int k = 0; k < M(); ++k) {
            const double s = raiz[k] * xi[k];
            const auto &f = P.fx[P.mapx[ix[k]]], &g = P.gy[P.mapy[iy[k]]];
            for (int p = 0; p < P.np; ++p) z[p] += s * f[p] * g[p];
        }
        return z;
    }

    /// Covariancia da expansao truncada entre os pontos p e q de P: sum_k lambda_k phi_k(p) phi_k(q)
    double Covariancia(const Pontos &P, int p, int q) const {
        double c = 0.;
        for (int k = 0; k < M(); ++k) {
            const auto &f = P.fx[P.mapx[ix[k]]], &g = P.gy[P.mapy[iy[k]]];
            c += Lambda(k) * f[p] * g[p] * f[q] * g[q];
        }
        return c;
    }
};

/// O nucleo pedido ("exp": rho1, exponencial; "gauss": rho2, gaussiano) no retangulo [x0, x0 + Dx] x [y0, y0 + Dy]
inline KLSep2D CriaKL(const std::string &nucleo, double x0, double y0, double Dx, double Dy, double hx, double hy,
                      double fracaoAlvo, int Mmax) {
    if (nucleo == "exp")
        return KLSep2D(x0, y0, Dx, Dy, BaseExponencial(Dx, hx, std::max(1, Mmax)), BaseExponencial(Dy, hy, std::max(1, Mmax)),
                       fracaoAlvo, Mmax);
    if (nucleo == "gauss")
        return KLSep2D(x0, y0, Dx, Dy, BaseGaussiana(Dx, hx), BaseGaussiana(Dy, hy), fracaoAlvo, Mmax);
    throw std::invalid_argument("CriaKL: nucleo deve ser exp ou gauss");
}

/// O nucleo exato (sem truncamento) entre dois pontos
inline double NucleoExato(const std::string &nucleo, double dx, double dy, double hx, double hy) {
    if (nucleo == "exp") return std::exp(-std::fabs(dx) / hx - std::fabs(dy) / hy);
    return std::exp(-(dx / hx) * (dx / hx) - (dy / hy) * (dy / hy));
}

} // namespace v2kl
