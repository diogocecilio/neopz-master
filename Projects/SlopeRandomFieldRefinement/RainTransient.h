// Transiente da chuva com o material nao saturado, nas hipoteses do SEEP/W de Rahardjo et al. (2007): esqueleto rigido
// (u = 0), ar passivo (p_a = 0) e so a percolacao da agua, por Euler implicito, com o TPZUnsatAnalysis do neopz.
// Contornos: p_w hidrostatica nas laterais abaixo do nivel d'agua inicial (ou de uma carga total dada); base e laterais
// acima do nivel sem fluxo; na superficie, o fluxo da chuva ou, onde a agua empocaria (p_w > 0), a face de percolacao
// (p_w = 0). A chuva entra por metro de superficie inclinada (normal a face: o "unit flux" do GeoStudio; o padrao) ou
// por metro de projecao horizontal (chuva vertical: rainVertical).
// E o RainTransient do Common (validado contra o Rahardjo), com os ids das laterais e da superficie como parametros
#pragma once

#include "TPZUnsatAnalysis.h"

#include <cmath>
#include <map>
#include <set>
#include <string>
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
    REAL dtMin = 0.;      // > 0: um passo que nao converge e refeito em duas metades, ate dtMin (h); 0: para (o antigo)
    int lineSearch = 4;   // metades do passo de Newton na busca linear (a do TPZUnsatAnalysis::Solve)
    bool rainVertical = false; // false: q = I por metro inclinado (normal a face); true: por metro de projecao horizontal
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
    // nivel inicial); superficie: ids da superficie (recebem a chuva); cargas: as laterais com uma carga total h dada (m,
    // cota do nivel: p_w = gamma_w (h - y) abaixo de h) no lugar do nivel inicial
    RainTransient(TPZUnsatAnalysis &an, const RainData &r, REAL gammaw, const std::set<int> &laterais,
                  const std::set<int> &superficie, const std::map<int, REAL> &cargas = {});

    // t = 0: p_w = gamma_w (nivel - y), com a succao limitada a suctionMax
    void InitialState();
    // Passos de dt ate o instante t (h); false se um passo nao convergiu (o estado fica no ultimo passo aceito)
    bool AdvanceTo(REAL t);
    // Percolacao permanente com a infiltracao q (m/s, na convencao da chuva: por metro inclinado ou da projecao
    // horizontal), a partir do estado atual: passos de Euler implicito com dt dobrando de dt0 (h) ate 1e8 h, e metade
    // quando um passo nao converge; para quando a maior mudanca de p_w num passo com dt >= 1e4 h fica abaixo de tol
    // (kPa). Nesses passos o armazenamento e desprezivel: o passo e a solucao permanente. O relogio volta a 0 (a chuva
    // comeca depois). false se nao estacionou (o estado fica no ultimo passo aceito); mudanca: a ultima maior mudanca
    bool SteadyState(REAL q, REAL dt0 = 0.1, REAL tol = 1.e-6, REAL *mudanca = nullptr);
    int64_t SteadySteps() const { return fSteadySteps; }

    REAL Time() const { return fTime / 3600.; } // h
    int64_t Steps() const { return fSteps; }
    int64_t NewtonIterations() const { return fNewton; } // soma de todos os passos
    size_t FaceNodes() const { return fFace.size(); }    // nos da superficie na face de percolacao agora
    const std::string &LastFailure() const { return fFalha; } // por que o ultimo passo nao convergiu (Newton ou face)
    std::pair<REAL, REAL> PressureRange() const;         // menor e maior p_w (kPa)
    // p_w de cada no da malha geometrica (pelo indice do no; 0 nos que nao tem p_w): o que vai para a outra malha
    std::vector<REAL> NodalWaterPressure() const;

private:
    // Um passo com a chuva media q (m/s): resolve com uma face tentativa e a corrige ate ela se repetir
    bool Step(REAL dt, REAL q);
    // Um passo de dt a partir de fTime com a chuva media do intervalo (e, com dtMin > 0, as metades se nao convergir)
    bool Advance(REAL dt);

    TPZUnsatAnalysis &fAn;
    RainData fRain;
    REAL fGammaW;
    std::map<int64_t, REAL> fFixed;   // equacoes sempre impostas: u = 0, p_a = 0 e p_w das laterais
    std::map<int64_t, REAL> fSurface; // no da superficie -> comprimento de influencia (m), que recebe a chuva
    std::set<int64_t> fFace;          // nos da face de percolacao no ultimo passo aceito
    REAL fTime = 0.;                  // s
    int64_t fSteps = 0, fNewton = 0, fSteadySteps = 0;
    std::string fFalha;
};
