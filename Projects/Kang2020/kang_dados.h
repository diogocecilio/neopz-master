// Os dados de Kang, Cho, Kim e Go (2020), "Effects of two-phase flow of water and air on shallow slope failures induced
// by rainfall: insights from slope stability assessment at a regional scale", Water 12(3):812, doi:10.3390/w12030812.
// Unidades: kN, m, s, kPa (a chuva em mm/h). O que o artigo nao diz explicitamente esta marcado como HIPOTESE ou LIDO
// (lido das figuras).
#pragma once

#include "Plasticity/TPZUnsaturatedHydraulics.h"

#include <array>
#include <vector>

namespace kang {

constexpr REAL gammaW = 9.81; // peso especifico da agua (kN/m3)
constexpr REAL H = 2.;        // espessura do solo (m), medida na VERTICAL: na Fig. 8, com a coluna saturada, Pw no fundo
                              // e ~19 kPa = gamma_w x 2 m (com 2 m normais ao talude de 30 graus seriam 22,6 kPa)
constexpr REAL pw0 = -20.;    // Pw inicial (kPa), constante com a profundidade; Pa inicial = 0 (Secao 2.3.3)
constexpr REAL duracao = 17.; // horas de chuva simuladas: das 16:00 de 26/07/2011 as 09:00 de 27/07 (Secao 3.2)
// o esqueleto: o artigo usa o FLAC acoplado mas nao da os modulos. HIPOTESE: os do Rahardjo (E = 20 MPa, nu = 0,3)
constexpr REAL E = 2.e4;  // kPa
constexpr REAL nu = 0.3;

// Tabela 2: a agua e o ar
constexpr REAL muRatio = 56.;   // mu_w / mu_a
constexpr REAL Kw = 2.e6;       // modulo volumetrico da agua (kPa; 2e9 Pa)
constexpr REAL rhoA = 1.25;     // massa especifica do ar (kg/m3)
constexpr REAL Ka = 100.;       // modulo volumetrico do ar (kPa; 1e5 Pa): o gas ideal isotermico a pressao atmosferica

/// Um lado do morro: as duas curvas de van Genuchten da Fig. 4 (norte: zonas 1 e 2; sul: zonas 3, 4 e 5)
struct Lado {
    const char *nome;
    REAL invP0;  // 1/P0 (1/kPa), Tabela 1
    REAL m;      // o "a" do artigo (Eq. 1): S_e = [1 + (psi/P0)^(1/(1-m))]^(-m), Tabela 1
    REAL thetaS; // teor de umidade volumetrico saturado = porosidade (LIDO da Fig. 4a)
    REAL thetaR; // teor residual (LIDO da Fig. 4a); S_rw = thetaR / thetaS
    REAL Srw() const { return thetaR / thetaS; }
    /// As relacoes hidraulicas na biblioteca: van Genuchten (Eq. 1), Mualem da agua (Eq. 4, b = 0,5) e do ar (Eq. 5,
    /// c = 0,5, Lenhard e Parker), as mesmas formulas do artigo
    TPZUnsaturatedHydraulics Hidraulica() const {
        TPZUnsaturatedHydraulics h;
        h.fRetention = TPZUnsaturatedHydraulics::EVanGenuchten;
        h.fAVG = invP0;         // a_VG = 1/P0
        h.fBVG = 1. / (1. - m); // n_VG, com m = 1 - 1/n_VG
        h.fSwr = Srw();
        h.fKrw = TPZUnsaturatedHydraulics::EMualem;
        h.fKra = TPZUnsaturatedHydraulics::EMualemAir;
        return h;
    }
};
inline const Lado &Norte() {
    static const Lado l{"norte", 0.051, 0.35, 0.50, 0.19};
    return l;
}
inline const Lado &Sul() {
    static const Lado l{"sul", 0.038, 0.63, 0.51, 0.183};
    return l;
}

/// Uma zona da Tabela 1, com os angulos dos taludes rompidos (legenda da Fig. 9) e o FS inicial da Fig. 9
struct Zona {
    int id;
    REAL phi;    // angulo de atrito (graus)
    REAL c;      // coesao (kPa)
    REAL gammaT; // peso especifico total (kN/m3)
    REAL gammaD; // peso especifico seco (kN/m3)
    REAL ks;     // condutividade hidraulica saturada (m/s)
    bool norte;  // zonas 1 e 2: a curva do norte; 3, 4 e 5: a do sul
    std::array<REAL, 3> beta;    // angulos maximo, medio e minimo (graus)
    std::array<REAL, 3> fs0Fig9; // FS em t = 0 na Fig. 9 para os tres angulos (LIDO, +-0,03)
    const Lado &L() const { return norte ? Norte() : Sul(); }
};
inline const std::vector<Zona> &Zonas() {
    static const std::vector<Zona> z = {
        {1, 25.3, 9.6, 18.1, 15.0, 7.15e-6, true, {38.8, 30.7, 26.5}, {1.57, 1.87, 2.13}},
        {2, 28.5, 5.8, 17.7, 14.9, 3.87e-6, true, {37.4, 31.0, 23.3}, {1.53, 1.80, 2.33}},
        {3, 37.6, 7.7, 17.0, 14.0, 1.80e-6, false, {40.5, 31.2, 23.3}, {2.10, 2.62, 3.42}},
        {4, 30.9, 7.6, 17.3, 14.5, 9.70e-6, false, {40.5, 31.8, 24.8}, {1.70, 2.07, 2.58}},
        {5, 28.2, 6.3, 18.2, 15.0, 3.69e-6, false, {40.3, 31.1, 24.7}, {1.47, 1.82, 2.25}}};
    return z;
}

/// A chuva horaria da Fig. 2 (estacao de Namhyun), DIGITALIZADA: mm/h na hora i, i = 0 das 16:00 as 17:00 de 26/07/2011
/// ate i = 16 das 08:00 as 09:00 de 27/07 (a barra do rotulo HH:00 e a chuva da hora anterior). Soma 381,5 mm; a curva
/// acumulada da figura da ~390 mm as 09:00 (as barras pequenas de ~0,5 mm/h nao foram lidas)
inline const std::vector<REAL> &Chuva() {
    static const std::vector<REAL> c = {7.7,  14.7, 24.3, 3.3,  2.6,  18.2, 0.0,  0.0, 0.0,
                                        20.9, 17.0, 27.5, 19.1, 27.3, 32.3, 86.3, 80.2};
    return c;
}

} // namespace kang
