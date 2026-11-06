// arco.h - V6: geometria do talude (Fig. 9) e do arco circular pelo pe (Eqs. 4 a 12) de Yuan et al. (2021), e o fator
// de reducao de variancia da media do campo ao longo do arco (Eqs. 3, 13, 14 e 15), por quadratura de Gauss-Legendre.
//
// Convencao da correlacao. O artigo (Eq. 13) usa rho = exp[-2 (|tau1| / delta1 + |tau2| / delta2)], delta = distancia de
// correlacao (a "escala de flutuacao" de Vanmarcke). O Comum/KLSeparavel.h usa exp(-|tau_x| / Lx - |tau_y| / Ly). Logo
// Lx = delta1 / 2 e Ly = delta2 / 2 (delta1 = 3 m e delta2 = 1 m dao Lx = 1,5 m e Ly = 0,5 m). A secao 4 do artigo (KL)
// escreve exp(-c |x1 - x2|) com c = 1 / delta, o que corresponderia a L = delta: a ambiguidade de fator 2 do artigo.
// Aqui a Eq. 13 (a das Tabelas 5 e 6) e a referencia; a outra leitura aparece so como sensibilidade.
//
// Fator de variancia: Gamma^2 = Var[X_L] / sigma^2 = (1 / T^2) int_0^T int_0^T rho_L(xi, eta) dxi deta (Eq. 3). Com
// s = eta - xi e u = (eta + xi) / 2 (jacobiano 1) e a simetria em s:
//     Gamma^2 = (2 / T^2) int_0^T ds int_{s/2}^{T - s/2} du rho(dx, dy),
//     dx = 2 r sin(s / 2r) cos(u / r + theta),  dy = 2 r sin(s / 2r) sin(u / r + theta)        (Eq. 14 / Eq. 15)
// (a forma fechada da Eq. 15, com o sinal corrigido: o artigo imprime cos(... - theta)). O integrando e suave em s > 0
// e tem bicos em u so onde cos(u / r + theta) ou sin(u / r + theta) se anulam (corda vertical ou horizontal): a
// quadratura de Gauss-Legendre composta quebra os paineis nesses pontos e converge exponencialmente.
#pragma once
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace v6 {

constexpr double kPi = 3.14159265358979323846;

/// Pontos e pesos de Gauss-Legendre em [-1, 1] (n pontos): Newton nas raizes do polinomio de Legendre P_n
inline void GaussLegendre(int n, std::vector<double> &x, std::vector<double> &w) {
    if (n < 1) throw std::invalid_argument("GaussLegendre: n < 1");
    x.assign(n, 0.);
    w.assign(n, 0.);
    for (int i = 0; i < (n + 1) / 2; ++i) {
        double z = std::cos(kPi * (i + 0.75) / (n + 0.5)), p1 = 1., p2 = 0.;
        for (int it = 0; it < 100; ++it) {
            p1 = 1.;
            p2 = 0.;
            for (int j = 1; j <= n; ++j) {
                const double p3 = p2;
                p2 = p1;
                p1 = ((2. * j - 1.) * z * p2 - (j - 1.) * p3) / j;
            }
            const double dp = n * (z * p1 - p2) / (z * z - 1.), dz = p1 / dp;
            z -= dz;
            if (std::fabs(dz) < 1.e-15) break;
        }
        p1 = 1.;
        p2 = 0.;
        for (int j = 1; j <= n; ++j) {
            const double p3 = p2;
            p2 = p1;
            p1 = ((2. * j - 1.) * z * p2 - (j - 1.) * p3) / j;
        }
        const double dp = n * (z * p1 - p2) / (z * z - 1.);
        x[i] = -z;
        x[n - 1 - i] = z;
        w[i] = w[n - 1 - i] = 2. / ((1. - z * z) * dp * dp);
    }
}

/// Talude da Fig. 9 (pe na origem): altura H, inclinacao beta (graus); terreno a esquerda do pe com 1,5 H, crista com
/// 2,5 H a direita da borda, fundacao com H abaixo do pe (altura total 2H). O retangulo [X0, X1] x [Y0, Y1] e o dominio
/// do campo aleatorio (o mesmo de um RFEM do talude).
struct Talude {
    double H = 10., beta = 40.;
    double XCrista() const { return H / std::tan(beta * kPi / 180.); }
    double X0() const { return -1.5 * H; }
    double X1() const { return XCrista() + 2.5 * H; }
    double Y0() const { return -H; }
    double Y1() const { return H; }
    double Superficie(double x) const { return x <= 0. ? 0. : (x >= XCrista() ? H : x * std::tan(beta * kPi / 180.)); }
};

/// Arco OA (Figs. 3 a 5): circulo de centro C(a, b) que passa pelo pe O = (0, 0), percorrido de O no sentido
/// anti-horario em torno de C; l e o comprimento de arco a partir de O (Eq. 12)
struct Arco {
    double a = 0., b = 0.;   ///< centro C
    double r = 0.;           ///< raio |OC| (Eq. 4)
    double T = 0.;           ///< comprimento (Eq. 8)
    double alpha = 0.;       ///< angulo central (Eq. 7)
    double theta = 0.;       ///< direcao da tangente em O (Eq. 11: -arctan(a / b), para b > 0)
    double M = 0., N = 0.;   ///< r cos(theta), r sin(theta) (Eq. 12)
    double c = 0., d = 0.;   ///< ponto final A = (c, d)

    /// Eqs. 4 a 8 e 11: o arco do pe ate o ponto A = (c, H) da crista (Eq. 6: a intersecao da direita com y = H)
    static Arco PeloPe(double a, double b, double H) {
        Arco A = Base(a, b);
        const double disc = a * a + b * b - (H - b) * (H - b);
        if (disc < 0.) throw std::invalid_argument("Arco::PeloPe: o circulo nao alcanca a cota H");
        A.c = std::sqrt(disc) + a; // Eq. 6
        A.d = H;
        const double corda = std::hypot(A.c, A.d);
        // Eq. 7 (vale para alpha <= pi: o centro do lado esquerdo da corda OA, como nas Figs. 3 e 4)
        if (a * A.d - b * A.c > 0.) throw std::invalid_argument("Arco::PeloPe: centro do lado errado da corda (alpha > pi)");
        A.alpha = 2. * std::asin(std::min(1., corda / (2. * A.r)));
        A.T = A.r * A.alpha; // Eq. 8
        return A;
    }
    /// arco de comprimento T a partir de O (os testes usam raios enormes para obter um segmento reto)
    static Arco ComComprimento(double a, double b, double T) {
        Arco A = Base(a, b);
        A.T = T;
        A.alpha = T / A.r;
        A.Ponto(T, A.c, A.d);
        return A;
    }
    /// Eq. 12: o ponto a distancia l de O ao longo do arco
    void Ponto(double l, double &x, double &y) const {
        const double g = l / r, sg = std::sin(g), umc = 2. * std::sin(0.5 * g) * std::sin(0.5 * g); // 1 - cos(g)
        x = M * sg - N * umc;
        y = N * sg + M * umc;
    }

private:
    static Arco Base(double a, double b) {
        Arco A;
        A.a = a;
        A.b = b;
        A.r = std::hypot(a, b); // Eq. 4
        if (A.r <= 0.) throw std::invalid_argument("Arco: centro na origem");
        A.theta = std::atan2(-a, b); // = -arctan(a / b) para b > 0 (Eq. 11)
        A.M = A.r * std::cos(A.theta);
        A.N = A.r * std::sin(A.theta);
        return A;
    }
};

/// O arco esta dentro do solo do talude (abaixo da superficie e acima da base do dominio) nos n pontos de teste?
inline bool ArcoDentro(const Arco &A, const Talude &t, int n = 400, double tol = 1.e-9) {
    for (int i = 1; i < n; ++i) {
        double x, y;
        A.Ponto(A.T * i / n, x, y);
        if (y > t.Superficie(x) + tol * t.H || y < t.Y0() || x < t.X0() || x > t.X1()) return false;
    }
    return true;
}

/// Correlacao separavel exp(-|dx| / Lx - |dy| / Ly) (a de KLSeparavel.h)
struct Correlacao {
    double Lx = 1.5, Ly = 0.5;
    double operator()(double dx, double dy) const { return std::exp(-std::fabs(dx) / Lx - std::fabs(dy) / Ly); }
    /// Eq. 13 do artigo: exp[-2 (|tau1| / delta1 + |tau2| / delta2)]  ->  L = delta / 2
    static Correlacao Eq13(double delta1, double delta2) { return {0.5 * delta1, 0.5 * delta2}; }
    /// a outra leitura (secao 4 do artigo: exp(-c |x1 - x2|), c = 1 / delta)  ->  L = delta
    static Correlacao CIgualUmSobreDelta(double delta1, double delta2) { return {delta1, delta2}; }
};

/// Opcoes da quadratura (Gauss-Legendre composta)
struct OpcoesQuad {
    int ng = 8;          ///< pontos por painel
    double fator = 1.;   ///< < 1 refina os paineis (teste de convergencia)
};

namespace detalhe {
/// (2 / T^2) int_0^smax ds int_{s/2}^{T-s/2} du f(s, u), com paineis de largura hs em s e hu em u, quebrados nos
/// pontos de u em "quebras" (os bicos do integrando)
template <class F>
double IntegralSU(double T, double smax, double hs, double hu, const std::vector<double> &quebras, int ng, F f) {
    std::vector<double> gx, gw;
    GaussLegendre(ng, gx, gw);
    const int ns = std::max(1, (int)std::ceil(smax / hs));
    const double ds = smax / ns;
    double total = 0.;
    std::vector<double> pts;
    for (int is = 0; is < ns; ++is)
        for (int i = 0; i < ng; ++i) {
            const double s = ds * (is + 0.5 * (gx[i] + 1.)), ws = 0.5 * ds * gw[i];
            const double u0 = 0.5 * s, u1 = T - 0.5 * s;
            if (u1 <= u0) continue;
            pts.assign(1, u0);
            for (double q : quebras)
                if (q > u0 && q < u1) pts.push_back(q);
            pts.push_back(u1);
            std::sort(pts.begin(), pts.end());
            double inner = 0.;
            for (size_t k = 0; k + 1 < pts.size(); ++k) {
                const double len = pts[k + 1] - pts[k];
                const int nu = std::max(1, (int)std::ceil(len / hu));
                const double du = len / nu;
                for (int iu = 0; iu < nu; ++iu)
                    for (int j = 0; j < ng; ++j) inner += 0.5 * du * gw[j] * f(s, pts[k] + du * (iu + 0.5 * (gx[j] + 1.)));
            }
            total += ws * inner;
        }
    return 2. * total / (T * T);
}
/// abscissas u em [0, T] onde u / r + desloc = k pi / 2 (os bicos de |cos| e |sin| do angulo u / r + desloc)
inline void Quebras(double T, double r, double desloc, std::vector<double> &q) {
    const double g0 = desloc, g1 = T / r + desloc;
    for (int k = (int)std::floor(g0 / (0.5 * kPi)) - 1; k <= (int)std::ceil(g1 / (0.5 * kPi)) + 1; ++k) {
        const double u = r * (0.5 * kPi * k - desloc);
        if (u > 0. && u < T) q.push_back(u);
    }
}
/// s maximo: a partir dele rho <= exp(-corda / Lmax) < e^-46 (a corda e 2 r sin(s / 2r) >= 2 s / pi para s <= pi r)
inline double SMax(const Arco &A, double Lmax) {
    const double alvo = 46. * Lmax;
    if (alvo >= 2. * A.r) return A.T;
    return std::min(A.T, 2. * A.r * std::asin(alvo / (2. * A.r)));
}
} // namespace detalhe

/// Gamma^2 = Var[X_L] / sigma^2 do arco (Eq. 3 com a Eq. 14; = Eq. 15 com o sinal corrigido)
inline double FatorVariancia(const Arco &A, const Correlacao &rho, const OpcoesQuad &op = {}) {
    // paineis: em s, menos de meio decaimento-e por painel; em u, 1/16 do raio (a direcao da corda varia com u / r)
    const double Lmax = std::max(rho.Lx, rho.Ly);
    const double hs = op.fator * 0.5 / (1. / rho.Lx + 1. / rho.Ly);
    const double hu = op.fator * std::min(A.T, A.r) / 16.;
    std::vector<double> q;
    detalhe::Quebras(A.T, A.r, A.theta, q);
    const double r = A.r, th = A.theta, ax = 2. * r / rho.Lx, ay = 2. * r / rho.Ly;
    return detalhe::IntegralSU(A.T, detalhe::SMax(A, Lmax), hs, hu, q, op.ng, [&](double s, double u) {
        const double S = std::sin(0.5 * s / r), g = u / r + th;
        return std::exp(-S * (ax * std::fabs(std::cos(g)) + ay * std::fabs(std::sin(g))));
    });
}

/// A ultima linha da Eq. 15 como impressa: exp[-4r (|sin((eta - xi)/2r) cos((eta + xi)/2r - theta)| / delta1 +
/// |sin((eta - xi)/2r) sin((eta + xi)/2r + theta)| / delta2)] (com 4r / delta = 2r / L); so para medir o erro de sinal
inline double FatorVarianciaEq15Impressa(const Arco &A, const Correlacao &rho, const OpcoesQuad &op = {}) {
    const double hs = op.fator * 0.5 / (1. / rho.Lx + 1. / rho.Ly), hu = op.fator * std::min(A.T, A.r) / 16.;
    std::vector<double> q;
    detalhe::Quebras(A.T, A.r, A.theta, q);  // bicos de sin(u / r + theta)
    detalhe::Quebras(A.T, A.r, -A.theta, q); // bicos de cos(u / r - theta)
    const double r = A.r, th = A.theta, ax = 2. * r / rho.Lx, ay = 2. * r / rho.Ly;
    // com o sinal trocado o expoente deixa de ser a distancia real entre os pontos e o corte em s de FatorVariancia
    // nao vale: integra-se s ate T
    return detalhe::IntegralSU(A.T, A.T, hs, hu, q, op.ng, [&](double s, double u) {
        const double S = std::sin(0.5 * s / r), g = u / r;
        return std::exp(-S * (ax * std::fabs(std::cos(g - th)) + ay * std::fabs(std::sin(g + th))));
    });
}

/// Forma fechada para um segmento reto de comprimento T com correlacao exp(-|s| / Le) ao longo dele (Vanmarcke 1983):
/// Gamma^2 = 2 (Le / T)^2 (T / Le - 1 + exp(-T / Le)). Para o nucleo separavel num segmento na direcao phi,
/// Le = 1 / (|cos phi| / Lx + |sin phi| / Ly).
inline double FatorVarianciaSegmento(double T, double Le) {
    const double t = T / Le;
    return 2. / (t * t) * (t - 1. + std::exp(-t));
}

/// Forca bruta: regra do ponto medio em n x n celulas de (xi, eta) com as coordenadas da Eq. 12 (teste independente)
inline double FatorVarianciaPontoMedio(const Arco &A, const Correlacao &rho, int n) {
    std::vector<double> x(n), y(n);
    for (int i = 0; i < n; ++i) A.Ponto(A.T * (i + 0.5) / n, x[i], y[i]);
    double s = 0.;
    for (int i = 0; i < n; ++i) {
        double si = 0.;
        for (int j = 0; j < n; ++j) si += rho(x[i] - x[j], y[i] - y[j]);
        s += si;
    }
    return s / ((double)n * n);
}

} // namespace v6
