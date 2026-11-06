// bishop.h - V6: o circulo critico pelo pe de um talude seco homogeneo c-phi pelo metodo de Bishop simplificado, para
// reconstruir os arcos das Tabelas 5 e 6 (o artigo ajusta um circulo aos nos de maior deformacao plastica do SRM no
// ABAQUS e nao imprime centro nem raio; so a Fig. 15 mostra um). E uma reconstrucao NOSSA: o peso especifico nao e dado
// no artigo e e calibrado para que o FS minimo do talude H = 50 m, beta = 40 graus seja o da Tabela 7 (0,9938).
//
// Bishop simplificado (sem agua): FS = sum[(c b + W tan phi) / m_alpha] / sum W sin(alpha),
// m_alpha = cos(alpha) + sin(alpha) tan(phi) / FS, com fatias verticais de x = 0 (o pe) ate a saida na crista e
// sin(alpha) = (x - a) / r (o talude da Fig. 9 desce para a esquerda: a massa gira no sentido anti-horario em torno de C).
#pragma once
#include "arco.h"

#include <limits>

namespace v6 {

/// Solo homogeneo (Tabela 2: medias de c e phi); gamma nao e dado no artigo
struct SoloMC {
    double c = 36.9;     ///< kPa
    double phi = 22.6;   ///< graus
    double gamma = 18.;  ///< kN/m3 (suposto; calibrado em CalibraGamma)
};

/// FS de Bishop simplificado do circulo de centro (a, b) pelo pe, da saida A = (xA, H) na crista ate o pe; infinito se
/// o circulo nao sai na crista ou se a saida fica na metade de cima do circulo (superficie em balanco)
inline double FSBishop(const Talude &t, const SoloMC &s, double a, double b, int nf = 400) {
    const double inf = std::numeric_limits<double>::infinity(), H = t.H, r = std::hypot(a, b);
    if (b <= H) return inf;
    const double disc = r * r - (H - b) * (H - b);
    if (disc <= 0.) return inf;
    const double xA = a + std::sqrt(disc);
    if (xA <= t.XCrista()) return inf;
    const double tf = std::tan(s.phi * kPi / 180.), dx = xA / nf;
    std::vector<double> W(nf), sa(nf), ca(nf);
    double mot = 0.;
    for (int i = 0; i < nf; ++i) {
        const double x = dx * (i + 0.5), yb = b - std::sqrt(std::max(0., r * r - (x - a) * (x - a)));
        const double h = std::max(0., t.Superficie(x) - yb);
        W[i] = s.gamma * h * dx;
        sa[i] = (x - a) / r;
        ca[i] = std::sqrt(std::max(0., 1. - sa[i] * sa[i]));
        mot += W[i] * sa[i];
    }
    if (mot <= 0.) return inf;
    double F = 1.;
    for (int it = 0; it < 200; ++it) {
        double res = 0.;
        for (int i = 0; i < nf; ++i) {
            const double m = ca[i] + sa[i] * tf / F;
            if (m < 0.05) return inf; // m_alpha pequeno: fatia de base muito inclinada para tras (circulo absurdo)
            res += (s.c * dx + W[i] * tf) / m;
        }
        const double Fn = res / mot;
        if (std::fabs(Fn - F) < 1.e-12 * Fn) { F = Fn; break; }
        F = Fn;
    }
    return F;
}

/// Circulo critico pelo pe: centro (a, b) e o FS minimo
struct CirculoCritico {
    double a = 0., b = 0., FS = 0.;
};

/// Busca do centro: grade grossa em a em [-3H, 2H], b em (H, 8H], depois grades locais 21 x 21 que encolhem 4 vezes
inline CirculoCritico BuscaCirculoCritico(const Talude &t, const SoloMC &s, int nf = 400) {
    const double H = t.H;
    CirculoCritico best{0., 0., std::numeric_limits<double>::infinity()};
    const int n = 60;
    for (int i = 0; i <= n; ++i)
        for (int j = 1; j <= n; ++j) {
            const double a = -3. * H + 5. * H * i / n, b = H + 7. * H * j / n, F = FSBishop(t, s, a, b, nf);
            if (F < best.FS) best = {a, b, F};
        }
    double ha = 5. * H / n, hb = 7. * H / n;
    for (int nivel = 0; nivel < 8; ++nivel) {
        const CirculoCritico c0 = best;
        for (int i = -10; i <= 10; ++i)
            for (int j = -10; j <= 10; ++j) {
                const double a = c0.a + ha * i / 10., b = c0.b + hb * j / 10., F = FSBishop(t, s, a, b, nf);
                if (F < best.FS) best = {a, b, F};
            }
        ha /= 4.;
        hb /= 4.;
    }
    return best;
}

/// gamma tal que o FS minimo do talude t seja FSalvo (secante em ln gamma; o FS cai com gamma)
inline double CalibraGamma(const Talude &t, SoloMC s, double FSalvo, int nf = 400) {
    auto f = [&](double lg) { s.gamma = std::exp(lg); return BuscaCirculoCritico(t, s, nf).FS - FSalvo; };
    double x0 = std::log(15.), x1 = std::log(20.), f0 = f(x0), f1 = f(x1);
    for (int it = 0; it < 30 && std::fabs(f1) > 1.e-9; ++it) {
        const double x2 = x1 - f1 * (x1 - x0) / (f1 - f0);
        x0 = x1;
        f0 = f1;
        x1 = x2;
        f1 = f(x1);
    }
    return std::exp(x1);
}

} // namespace v6
