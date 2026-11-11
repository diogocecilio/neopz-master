// KLSeparavel.h (copia do SlopeRandomField, sem mudancas) - expansao de Karhunen-Loeve de um campo gaussiano estacionario de media 0 e variancia 1 com covariancia
// exponencial separavel
//     C(x1, x2) = exp(-|x1 - x2| / Lx - |y1 - y2| / Ly)
// (Vargas Ceron, Cecilio, Linn e Maghous 2025, eq. 64) num retangulo [x0, x0 + Dx] x [y0, y0 + Dy] que contem o dominio.
// Como o nucleo e o produto de dois nucleos exponenciais 1D, os autopares 2D sao os produtos dos autopares 1D, que tem
// forma fechada (Ghanem e Spanos 1991, sec. 2.3.3; e a abordagem de Cho 2010 para o retangulo):
//   nucleo exp(-c |s - t|), c = 1/L, em [-a, a]:
//     pares:   omega tan(omega a) = c,   f(s) = cos(omega s) / sqrt(a + sin(2 omega a) / (2 omega))
//     impares: omega + c tan(omega a) = 0, f(s) = sin(omega s) / sqrt(a - sin(2 omega a) / (2 omega))
//     lambda = 2 c / (omega^2 + c^2); a soma de todos os lambda e 2a (o traco do nucleo).
// O campo truncado e Z(x) = sum_k sqrt(lambda_k) xi_k f_k(x), com os M maiores produtos lambda_i^x lambda_j^y; a fracao de
// variancia retida (media no retangulo) e sum lambda_k / (Dx Dy). Sem dependencia do NeoPZ.
#pragma once
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <vector>

namespace klsep {

/// Autopares do nucleo exp(-|s - t| / L) em [-a, a], em ordem decrescente de lambda
struct Modos1D {
    double a = 0., L = 0.;
    std::vector<double> lambda, omega;
    std::vector<int> par; ///< 1: cosseno (par), 0: seno (impar)
    std::vector<double> norma; ///< 1/sqrt(integral de cos^2 ou sin^2)

    double Avalia(int i, double s) const { // s medido a partir do centro do intervalo
        return norma[i] * (par[i] ? std::cos(omega[i] * s) : std::sin(omega[i] * s));
    }
};

namespace detalhe {
// raiz de g em (lo, hi) por bissecao; g muda de sinal no intervalo e e monotona nele
template <class G> double Bissecao(G g, double lo, double hi) {
    double glo = g(lo);
    for (int it = 0; it < 200; ++it) {
        const double m = 0.5 * (lo + hi), gm = g(m);
        if ((gm > 0.) == (glo > 0.)) { lo = m; glo = gm; } else hi = m;
        if (hi - lo < 1.e-15 * std::max(1., hi)) break;
    }
    return 0.5 * (lo + hi);
}
} // namespace detalhe

/// Os primeiros n autopares 1D (n/2 pares e n/2 impares, aproximadamente), ordenados por lambda decrescente
inline Modos1D CalculaModos1D(double comprimento, double L, int n) {
    if (comprimento <= 0. || L <= 0. || n <= 0) throw std::invalid_argument("CalculaModos1D: argumentos invalidos");
    Modos1D m;
    m.a = 0.5 * comprimento;
    m.L = L;
    const double a = m.a, c = 1. / L, pi = M_PI, eps = 1.e-12;
    struct Par { double lam, w; int par; };
    std::vector<Par> v;
    for (int k = 1; (int)v.size() < n; ++k) {
        // par: omega tan(omega a) - c = 0 em ((k-1) pi, (k-1/2) pi) / a
        {
            const double lo = ((k - 1) * pi + eps) / a, hi = ((k - 0.5) * pi - eps) / a;
            const double w = detalhe::Bissecao([&](double x) { return x * std::tan(x * a) - c; }, lo, hi);
            v.push_back({2. * c / (w * w + c * c), w, 1});
        }
        if ((int)v.size() >= n) break;
        // impar: omega + c tan(omega a) = 0 em ((k-1/2) pi, k pi) / a
        {
            const double lo = ((k - 0.5) * pi + eps) / a, hi = (k * pi - eps) / a;
            const double w = detalhe::Bissecao([&](double x) { return x + c * std::tan(x * a); }, lo, hi);
            v.push_back({2. * c / (w * w + c * c), w, 0});
        }
    }
    std::sort(v.begin(), v.end(), [](const Par &p, const Par &q) { return p.lam > q.lam; });
    for (const auto &p : v) {
        m.lambda.push_back(p.lam);
        m.omega.push_back(p.w);
        m.par.push_back(p.par);
        const double s2 = std::sin(2. * p.w * a) / (2. * p.w);
        m.norma.push_back(1. / std::sqrt(p.par ? a + s2 : a - s2));
    }
    return m;
}

/// Expansao 2D truncada: os M maiores produtos lambda_i^x lambda_j^y, ate atingir a fracao de variancia pedida
struct KL2D {
    double x0 = 0., y0 = 0., Dx = 0., Dy = 0.;
    Modos1D mx, my;
    std::vector<int> ix, iy;    ///< modo 1D em x e em y do termo k
    std::vector<double> raiz;   ///< sqrt(lambda_k)
    double fracao = 0.;         ///< fracao de variancia retida (media no retangulo)

    int M() const { return (int)raiz.size(); }

    /// fracaoAlvo: para ao atingir essa fracao; Mmax: limite de termos (o que vier primeiro)
    KL2D(double x0_, double y0_, double Dx_, double Dy_, double Lx, double Ly, double fracaoAlvo, int Mmax)
        : x0(x0_), y0(y0_), Dx(Dx_), Dy(Dy_) {
        // se (i, j) esta entre os M maiores produtos, os (i+1)(j+1) pares (i' <= i, j' <= j) tambem estao, pois os
        // lambda 1D estao em ordem decrescente: basta considerar (i+1)(j+1) <= M, cerca de M ln M candidatos
        const int n1 = std::max(1, Mmax);
        mx = CalculaModos1D(Dx, Lx, n1);
        my = CalculaModos1D(Dy, Ly, n1);
        struct T { double lam; int i, j; };
        std::vector<T> cand;
        for (int i = 0; i < n1; ++i)
            for (int j = 0; j < n1 && (long)(i + 1) * (j + 1) <= Mmax; ++j)
                cand.push_back({mx.lambda[i] * my.lambda[j], i, j});
        std::sort(cand.begin(), cand.end(), [](const T &p, const T &q) { return p.lam > q.lam; });
        const double total = Dx * Dy; // traco do nucleo 2D no retangulo
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

    /// Prepara a avaliacao em pontos fixos (os pontos de integracao): guarda f_i(x_p) e g_j(y_p) dos modos usados
    struct Pontos {
        int np = 0;
        std::vector<int> mapx, mapy;           ///< modo 1D -> coluna em fx/gy (-1 se nao usado)
        std::vector<std::vector<double>> fx, gy; ///< [coluna][ponto]
        std::vector<double> var;               ///< variancia truncada em cada ponto (sum lambda_k f_k^2)
    };

    Pontos PreparaPontos(const std::vector<double> &x, const std::vector<double> &y) const {
        Pontos P;
        P.np = (int)x.size();
        P.mapx.assign(mx.lambda.size(), -1);
        P.mapy.assign(my.lambda.size(), -1);
        for (int k = 0; k < M(); ++k) {
            if (P.mapx[ix[k]] < 0) { P.mapx[ix[k]] = (int)P.fx.size(); P.fx.emplace_back(P.np); }
            if (P.mapy[iy[k]] < 0) { P.mapy[iy[k]] = (int)P.gy.size(); P.gy.emplace_back(P.np); }
        }
        for (size_t i = 0; i < P.mapx.size(); ++i)
            if (P.mapx[i] >= 0)
                for (int p = 0; p < P.np; ++p) P.fx[P.mapx[i]][p] = mx.Avalia((int)i, x[p] - (x0 + 0.5 * Dx));
        for (size_t j = 0; j < P.mapy.size(); ++j)
            if (P.mapy[j] >= 0)
                for (int p = 0; p < P.np; ++p) P.gy[P.mapy[j]][p] = my.Avalia((int)j, y[p] - (y0 + 0.5 * Dy));
        P.var.assign(P.np, 0.);
        for (int k = 0; k < M(); ++k) {
            const auto &f = P.fx[P.mapx[ix[k]]], &g = P.gy[P.mapy[iy[k]]];
            const double l = raiz[k] * raiz[k];
            for (int p = 0; p < P.np; ++p) P.var[p] += l * f[p] * f[p] * g[p] * g[p];
        }
        return P;
    }

    /// Campo gaussiano padrao Z nos pontos, para as variaveis normais padrao xi (tamanho M)
    std::vector<double> Campo(const Pontos &P, const std::vector<double> &xi) const {
        if ((int)xi.size() != M()) throw std::invalid_argument("KL2D::Campo: xi deve ter M valores");
        std::vector<double> z(P.np, 0.);
        for (int k = 0; k < M(); ++k) {
            const double s = raiz[k] * xi[k];
            const auto &f = P.fx[P.mapx[ix[k]]], &g = P.gy[P.mapy[iy[k]]];
            for (int p = 0; p < P.np; ++p) z[p] += s * f[p] * g[p];
        }
        return z;
    }
};

/// Media e CoV de uma variavel lognormal -> parametros da normal de ln H (Vargas Ceron et al. 2025, eqs. 65-66)
struct Lognormal {
    double mu_ln = 0., sigma_ln = 0.;
    Lognormal(double media, double cov) {
        sigma_ln = std::sqrt(std::log(1. + cov * cov));
        mu_ln = std::log(media / std::sqrt(1. + cov * cov));
    }
    double operator()(double z) const { return std::exp(mu_ln + sigma_ln * z); }
};

/// Gerador reprodutivel em qualquer compilador: xoshiro256** semeado por splitmix64, normal por Box-Muller
class Normal {
    uint64_t s[4];
    bool tem = false;
    double guarda = 0.;
    static uint64_t SplitMix(uint64_t &x) {
        uint64_t z = (x += 0x9e3779b97f4a7c15ULL);
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
        z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
        return z ^ (z >> 31);
    }
    static uint64_t Rot(uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }
    uint64_t Proximo() {
        const uint64_t r = Rot(s[1] * 5, 7) * 9, t = s[1] << 17;
        s[2] ^= s[0]; s[3] ^= s[1]; s[1] ^= s[2]; s[0] ^= s[3]; s[2] ^= t; s[3] = Rot(s[3], 45);
        return r;
    }
    double Uniforme() { return ((Proximo() >> 11) + 0.5) * (1.0 / 9007199254740992.0); } // em (0, 1)
public:
    explicit Normal(uint64_t semente) { for (auto &v : s) v = SplitMix(semente); }
    double operator()() {
        if (tem) { tem = false; return guarda; }
        const double u1 = Uniforme(), u2 = Uniforme(), r = std::sqrt(-2. * std::log(u1)), t = 2. * M_PI * u2;
        guarda = r * std::sin(t);
        tem = true;
        return r * std::cos(t);
    }
    std::vector<double> Vetor(int n) { std::vector<double> v(n); for (auto &x : v) x = (*this)(); return v; }
};

} // namespace klsep
