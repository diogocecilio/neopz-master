// klarco.h - V6: o campo de Karhunen-Loeve de Comum/KLSeparavel.h ao longo de um arco: a media de cada termo no arco,
// as variancias exatas do campo truncado (no ponto, entre dois pontos e da media no arco) e o Monte Carlo.
//
// Campo truncado (como os drivers o usam: KL2D::Campo, sem renormalizar): Z_M(P) = sum_k sqrt(lambda_k) xi_k phi_k(P),
// phi_k(x, y) = f_i(x) g_j(y). A media no arco e X_M = sum_k sqrt(lambda_k) xi_k F_k, F_k = (1 / T) int phi_k ds, logo
//     Var[X_M] = sum_k lambda_k F_k^2        (cresce com M e tende a Gamma^2 da Eq. 14 pelo teorema de Mercer),
//     Var[Z_M(P)] = sum_k lambda_k phi_k(P)^2  (a media no retangulo e a fracao retida sum lambda_k / (Dx Dy)).
// O Monte Carlo sorteia xi ~ N(0, 1) com o gerador de KLSeparavel.h e mede a media no arco de duas formas identicas
// em aritmetica exata: (a) KL2D::Campo nos pontos de quadratura do arco e a soma ponderada (o caminho dos drivers) e
// (b) sum_k sqrt(lambda_k) xi_k F_k (rapida: M operacoes por realizacao). O teste confere (a) = (b) e usa (b) nas
// amostras grandes.
#pragma once
#include "arco.h"

#include "KLSeparavel.h"

#include <cstdint>
#include <memory>
#include <vector>

namespace v6 {

/// Pontos de quadratura ao longo do arco (Gauss-Legendre composta, paineis de comprimento <= h, ng pontos cada); os
/// pesos somam 1, entao sum_p w_p z(P_p) e a media de z no arco
struct QuadArco {
    std::vector<double> l, x, y, w;
    int N() const { return (int)w.size(); }
};

inline QuadArco CriaQuadArco(const Arco &A, double h, int ng = 8) {
    std::vector<double> gx, gw;
    GaussLegendre(ng, gx, gw);
    const int np = std::max(1, (int)std::ceil(A.T / h));
    const double dl = A.T / np;
    QuadArco q;
    for (int i = 0; i < np; ++i)
        for (int j = 0; j < ng; ++j) {
            const double l = dl * (i + 0.5 * (gx[j] + 1.));
            double x, y;
            A.Ponto(l, x, y);
            q.l.push_back(l);
            q.x.push_back(x);
            q.y.push_back(y);
            q.w.push_back(0.5 * dl * gw[j] / A.T);
        }
    return q;
}

/// Maior frequencia (rad/m) de um termo da KL ao longo de qualquer curva: max_k (omega_x(i_k) + omega_y(j_k))
inline double FrequenciaMaxima(const klsep::KL2D &kl) {
    double w = 0.;
    for (int k = 0; k < kl.M(); ++k) w = std::max(w, kl.mx.omega[kl.ix[k]] + kl.my.omega[kl.iy[k]]);
    return w;
}

/// Quadratura do arco fina o bastante para todos os termos da KL: omega_max h <= 5 com 8 pontos por painel (o erro
/// de Gauss-Legendre em exp(i omega s) num painel fica < (omega h / 2)^16 / 16! ~ 1e-7)
inline QuadArco CriaQuadArcoParaKL(const Arco &A, const klsep::KL2D &kl, double hmax = 0.5) {
    const double w = FrequenciaMaxima(kl);
    return CriaQuadArco(A, std::min(hmax, 5. / std::max(w, 1.e-12)), 8);
}

/// F_k = sum_p w_p f_{i_k}(x_p) g_{j_k}(y_p): a media de cada termo da KL no arco (sem o fator sqrt(lambda_k))
inline std::vector<double> MediasDosTermos(const klsep::KL2D &kl, const QuadArco &q) {
    const int M = kl.M();
    std::vector<int> colx(kl.mx.lambda.size(), -1), coly(kl.my.lambda.size(), -1), modx, mody;
    for (int k = 0; k < M; ++k) {
        if (colx[kl.ix[k]] < 0) { colx[kl.ix[k]] = (int)modx.size(); modx.push_back(kl.ix[k]); }
        if (coly[kl.iy[k]] < 0) { coly[kl.iy[k]] = (int)mody.size(); mody.push_back(kl.iy[k]); }
    }
    std::vector<int> cx(M), cy(M);
    for (int k = 0; k < M; ++k) { cx[k] = colx[kl.ix[k]]; cy[k] = coly[kl.iy[k]]; }
    const double xc = kl.x0 + 0.5 * kl.Dx, yc = kl.y0 + 0.5 * kl.Dy;
    std::vector<double> F(M, 0.), fx(modx.size()), gy(mody.size());
    for (int p = 0; p < q.N(); ++p) {
        for (size_t i = 0; i < modx.size(); ++i) fx[i] = kl.mx.Avalia(modx[i], q.x[p] - xc);
        for (size_t j = 0; j < mody.size(); ++j) gy[j] = q.w[p] * kl.my.Avalia(mody[j], q.y[p] - yc);
        const double *pfx = fx.data(), *pgy = gy.data();
        for (int k = 0; k < M; ++k) F[k] += pfx[cx[k]] * pgy[cy[k]];
    }
    return F;
}

/// Somas acumuladas de lambda_k F_k^2: V[m] = Var[X] com os m + 1 primeiros termos (os termos estao em ordem
/// decrescente de lambda, entao V[M' - 1] e a KL de M' termos)
inline std::vector<double> VarianciaMediaAcumulada(const klsep::KL2D &kl, const std::vector<double> &F) {
    std::vector<double> V(F.size());
    double s = 0.;
    for (size_t k = 0; k < F.size(); ++k) {
        s += kl.raiz[k] * kl.raiz[k] * F[k] * F[k];
        V[k] = s;
    }
    return V;
}

/// phi_k(P) para os pontos dados: [ponto][k]
inline std::vector<std::vector<double>> TermosNosPontos(const klsep::KL2D &kl, const std::vector<double> &x,
                                                        const std::vector<double> &y) {
    const double xc = kl.x0 + 0.5 * kl.Dx, yc = kl.y0 + 0.5 * kl.Dy;
    std::vector<std::vector<double>> phi(x.size(), std::vector<double>(kl.M()));
    for (size_t p = 0; p < x.size(); ++p)
        for (int k = 0; k < kl.M(); ++k)
            phi[p][k] = kl.mx.Avalia(kl.ix[k], x[p] - xc) * kl.my.Avalia(kl.iy[k], y[p] - yc);
    return phi;
}

/// Covariancia exata do campo truncado com os m primeiros termos: sum_{k<m} lambda_k phi_k(P) phi_k(Q)
inline double CovarianciaTruncada(const klsep::KL2D &kl, const std::vector<double> &phiP,
                                  const std::vector<double> &phiQ, int m = -1) {
    if (m < 0 || m > kl.M()) m = kl.M();
    double s = 0.;
    for (int k = 0; k < m; ++k) s += kl.raiz[k] * kl.raiz[k] * phiP[k] * phiQ[k];
    return s;
}

/// Media e variancia amostrais (Welford) e a correlacao amostral de pares de sondas
struct Momentos {
    long n = 0;
    double media = 0., m2 = 0.;
    void Soma(double v) {
        ++n;
        const double d = v - media;
        media += d / n;
        m2 += d * (v - media);
    }
    double Var() const { return n > 1 ? m2 / (n - 1) : 0.; }
    double DP() const { return std::sqrt(Var()); }
};

struct CoMomento { // correlacao de Pearson amostral de (u, v)
    long n = 0;
    double mu = 0., mv = 0., cuu = 0., cvv = 0., cuv = 0.;
    void Soma(double u, double v) {
        ++n;
        const double du = u - mu, dv = v - mv;
        mu += du / n;
        mv += dv / n;
        cuu += du * (u - mu);
        cvv += dv * (v - mv);
        cuv += du * (v - mv);
    }
    double Rho() const { return cuv / std::sqrt(cuu * cvv); }
};

/// Resultado do Monte Carlo de um arco
struct ResultadoMC {
    long N = 0;
    Momentos media;               ///< X_M (media no arco) por sum_k sqrt(lambda_k) xi_k F_k
    std::vector<Momentos> sonda;  ///< Z_M em cada sonda (KL2D::Campo)
    std::vector<CoMomento> par;   ///< pares de sondas (correlacao)
    double maxDifCampo = 0.;      ///< max |X por Campo nos pontos do arco - X por F_k| nas realizacoes conferidas
    long nConferidas = 0;
    double dpExato = 0.;          ///< sqrt(sum lambda_k F_k^2): o desvio exato da media no arco da KL truncada
    double s3 = 0., s4 = 0.;      ///< somas de (X / dpExato)^3 e ^4 (a media exata e 0)
    long nFora = 0;               ///< realizacoes com |X| > 1,959964 dpExato (5% numa normal)
    // com a media e o desvio exatos (nao os amostrais) os erros padrao numa normal sao sqrt(E[u^6] / N) = sqrt(15 / N)
    // para a assimetria e sqrt((E[u^8] - 9) / N) = sqrt(96 / N) para o excesso de curtose
    double Assimetria() const { return s3 / N; }
    double ExcessoCurtose() const { return s4 / N - 3.; }
};

/// Monte Carlo: N realizacoes (semente fixa); sondas: pontos onde se guarda o campo (Campo) e pares (i, j) de sondas
/// para a correlacao; nConferir realizacoes tambem calculam a media no arco por Campo nos pontos de quadratura
inline ResultadoMC MonteCarloArco(const klsep::KL2D &kl, const std::vector<double> &F, const QuadArco &q, long N,
                                  uint64_t semente, const std::vector<double> &sx, const std::vector<double> &sy,
                                  const std::vector<std::pair<int, int>> &pares, long nConferir = 0) {
    ResultadoMC R;
    R.N = N;
    R.sonda.resize(sx.size());
    R.par.resize(pares.size());
    klsep::Normal rng(semente);
    const auto Psondas = kl.PreparaPontos(sx, sy);
    std::unique_ptr<klsep::KL2D::Pontos> Parco;
    if (nConferir > 0) Parco.reset(new klsep::KL2D::Pontos(kl.PreparaPontos(q.x, q.y)));
    const int M = kl.M();
    std::vector<double> sf(M);
    double v = 0.;
    for (int k = 0; k < M; ++k) {
        sf[k] = kl.raiz[k] * F[k];
        v += sf[k] * sf[k];
    }
    R.dpExato = std::sqrt(v);
    for (long n = 0; n < N; ++n) {
        const std::vector<double> xi = rng.Vetor(M);
        double X = 0.;
        for (int k = 0; k < M; ++k) X += sf[k] * xi[k];
        R.media.Soma(X);
        const double u = X / R.dpExato;
        R.s3 += u * u * u;
        R.s4 += u * u * u * u;
        if (std::fabs(u) > 1.959964) ++R.nFora;
        if (!sx.empty()) {
            const auto z = kl.Campo(Psondas, xi);
            for (size_t i = 0; i < z.size(); ++i) R.sonda[i].Soma(z[i]);
            for (size_t i = 0; i < pares.size(); ++i) R.par[i].Soma(z[pares[i].first], z[pares[i].second]);
        }
        if (n < nConferir) {
            const auto z = kl.Campo(*Parco, xi);
            double Xc = 0.;
            for (int p = 0; p < q.N(); ++p) Xc += q.w[p] * z[p];
            R.maxDifCampo = std::max(R.maxDifCampo, std::fabs(Xc - X));
            ++R.nConferidas;
        }
    }
    return R;
}

/// Erro padrao do desvio padrao amostral de uma normal (N grande): sd / sqrt(2 (N - 1))
inline double EPdoDP(double sd, long N) { return sd / std::sqrt(2. * (N - 1.)); }
/// Erro padrao da variancia amostral de uma normal: var sqrt(2 / (N - 1))
inline double EPdaVar(double var, long N) { return var * std::sqrt(2. / (N - 1.)); }
/// Erro padrao da correlacao amostral (N grande): (1 - rho^2) / sqrt(N)
inline double EPdaCorrelacao(double rho, long N) { return (1. - rho * rho) / std::sqrt((double)N); }

} // namespace v6
