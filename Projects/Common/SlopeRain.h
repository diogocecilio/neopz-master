// Transiente da chuva de Rahardjo et al. (2007) com o material nao saturado, nas hipoteses do SEEP/W do artigo: esqueleto
// rigido (u = 0), ar passivo (p_a = 0) e so a percolacao da agua, por Euler implicito. Contornos (Fig. 2): p_w hidrostatica
// nas laterais abaixo do nivel inicial; base e laterais acima do nivel sem fluxo; na superficie, o fluxo da chuva ou, onde a
// agua empocaria (p_w > 0) ou sairia do talude, a face de percolacao (p_w = 0)
#ifndef SLOPERAIN_H
#define SLOPERAIN_H

#include "SlopeModel.h"

#include "TPZUnsatAnalysis.h"

#include <set>
#include <utility>

namespace slope {

class RainTransient {
public:
    // an: analise da malha de CreateRainCompMesh, com u e p_a ainda nulos
    RainTransient(TPZUnsatAnalysis &an, const RainData &r, REAL gammaw);

    // t = 0: p_w = gamma_w (nivel - y), com a succao limitada a suctionMax
    void InitialState();
    // Passos de dt ate o instante t (h); false se um passo nao convergiu (o estado fica no ultimo passo aceito)
    bool AdvanceTo(REAL t);

    REAL Time() const { return fTime / 3600.; } // h
    int64_t Steps() const { return fSteps; }
    int64_t NewtonIterations() const { return fNewton; } // soma de todos os passos
    size_t FaceNodes() const { return fFace.size(); }    // nos da superficie na face de percolacao agora
    std::pair<REAL, REAL> PressureRange() const;         // menor e maior p_w (kPa)

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

} // namespace slope

#endif
