// Transiente da chuva com o material nao saturado, nas hipoteses do SEEP/W de Rahardjo et al. (2007): esqueleto rigido
// (u = 0), ar passivo (p_a = 0) e so a percolacao da agua, por Euler implicito, com o TPZUnsatAnalysis do neopz.
// Contornos: p_w hidrostatica nas laterais abaixo do nivel d'agua inicial; base e laterais acima do nivel sem fluxo;
// na superficie, o fluxo da chuva ou, onde a agua empocaria (p_w > 0), a face de percolacao (p_w = 0).
// E o RainTransient do Common (validado contra o Rahardjo), com os ids das laterais e da superficie como parametros
#pragma once

#include "TPZUnsatAnalysis.h"

#include <cmath>
#include <map>
#include <set>
#include <utility>
#include <vector>

// Solo, chuva e estado inicial (kPa, m, s; tempos em horas). Padrao: Rahardjo et al. (2007), Serie A, solo f50,-5 sob
// I_r = 36 mm/h durante 24 h, com o nivel d'agua inicial 5 m abaixo do pe, subindo 7 graus para dentro do macico
struct RainData {
    // solo: retencao de Fredlund-Xing com C(psi) = 1 (a em kPa) e k_w = k_s Theta^power (Leong e Rahardjo 1997)
    REAL ks = 1.e-5;               // condutividade hidraulica saturada (m/s)
    REAL n = 0.45;                 // porosidade
    REAL aFX = 50., nFX = 1., mFX = 1., power = 4.;
    REAL phib = 26. * M_PI / 180.; // phi^b (resistencia da succao): chi = tan(phi^b) / tan(phi')
    REAL Kw = 2.2e6;               // modulo volumetrico da agua (kPa)
    // chuva
    REAL intensity = 36.; // I_r (mm/h)
    REAL duration = 24.;  // horas de chuva
    REAL dt = 0.1;        // passo de tempo (h)
    // estado inicial: hidrostatico abaixo do nivel d'agua (reta que sobe para a esquerda), succao ate suctionMax
    REAL xToe = 40., yToe = 30.;   // pe do talude (o da TriGMesh)
    REAL Hw = 5.;                  // profundidade do nivel abaixo do pe (m)
    REAL betaW = 7. * M_PI / 180.; // inclinacao do nivel
    REAL suctionMax = 75.;         // kPa
    REAL WaterLevel(REAL x) const { return yToe - Hw + (xToe - x) * std::tan(betaW); }
    REAL Intensity() const { return intensity / 3.6e6; } // m/s
};

class RainTransient {
public:
    // an: analise da malha do transiente, com u e p_a nulos; laterais: ids das laterais (p_w hidrostatica abaixo do
    // nivel inicial); superficie: ids da superficie (recebem a chuva)
    RainTransient(TPZUnsatAnalysis &an, const RainData &r, REAL gammaw, const std::set<int> &laterais,
                  const std::set<int> &superficie);

    // t = 0: p_w = gamma_w (nivel - y), com a succao limitada a suctionMax
    void InitialState();
    // Passos de dt ate o instante t (h); false se um passo nao convergiu (o estado fica no ultimo passo aceito)
    bool AdvanceTo(REAL t);

    REAL Time() const { return fTime / 3600.; } // h
    int64_t Steps() const { return fSteps; }
    int64_t NewtonIterations() const { return fNewton; } // soma de todos os passos
    size_t FaceNodes() const { return fFace.size(); }    // nos da superficie na face de percolacao agora
    std::pair<REAL, REAL> PressureRange() const;         // menor e maior p_w (kPa)
    // p_w de cada no da malha geometrica (pelo indice do no; 0 nos que nao tem p_w): o que vai para a outra malha
    std::vector<REAL> NodalWaterPressure() const;

private:
    // Um passo com a chuva media q (m/s): resolve com uma face tentativa e a corrige ate ela se repetir
    bool Step(REAL dt, REAL q);

    TPZUnsatAnalysis &fAn;
    RainData fRain;
    REAL fGammaW;
    std::map<int64_t, REAL> fFixed;   // equacoes sempre impostas: u = 0, p_a = 0 e p_w das laterais
    std::map<int64_t, REAL> fSurface; // no da superficie -> comprimento de influencia (m), que recebe a chuva
    std::set<int64_t> fFace;          // nos da face de percolacao no ultimo passo aceito
    REAL fTime = 0.;                  // s
    int64_t fSteps = 0, fNewton = 0;
};
