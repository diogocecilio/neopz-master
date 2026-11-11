// Campo de velocidade de filtracao semi-analitico de Vargas Ceron et al. (2025), Sec. 4.1: o melhor campo da familia
// (29) pelo minimo do funcional J*(v) (Eq. 23), permeabilidade isotropica (alpha = 1). A forca de percolacao e
// -grad u = K^-1 v = v / k, independente de k. Polares em O (borda da crista): theta = 0 ao longo da crista e
// theta = omega = pi - beta ao longo da face; e_r = (-cos, -sin), e_theta = (sin, -cos), y para cima.
//   r < R_w = h_w / sin(beta):  f1 = K0 [-(1 - s^(m-1)) h2'(theta) e_r + (1 - m s^(m-1)) h2(theta) e_theta], s = r/R_w,
//                               h2 = cosh(sqrt(m) theta), m = sqrt(C/D), C = int h2'^2, D = int h2^2, ou seja
//                               Sigma (1 - m) = omega (1 + m) / 2 com Sigma = sinh(2 sqrt(m) omega) / (4 sqrt(m)),
//                               K0 = gamma_w sin(beta) h2(omega) / (D - C)                          (Eqs. 29, 31, 36-39)
//   R_w <= r < R = H/sin(beta): f2 = gamma_w h_w / (omega r) e_theta                                (Eq. 32)
//   R <= r < R_e:               f3 = gamma_w h_w / (theta_max(r) r) e_theta, theta_max = pi - asin(H/r)  (Eq. 33)
//   r >= R_e:                   f = 0;  R_e = sqrt(H^2 + (L_m + H / tan(beta))^2)
#ifndef SEMIANALYTICALSEEPAGE_H
#define SEMIANALYTICALSEEPAGE_H

#include "SlopeModel.h"

#include <array>

class SemiAnalyticalSeepage {
public:
    SemiAnalyticalSeepage(const slope::Data &d, REAL Lm);
    std::array<REAL, 2> Force(const TPZVec<REAL> &x) const; // -grad u (kN/m3) num ponto da malha
    std::array<REAL, 2> ForcePolar(REAL r, REAL theta) const; // (f_r, f_theta)
    REAL JStarClosed() const;  // J*/k da Eq. (40)
    REAL JStarNumeric() const; // -1/2 int |f|^2 (o minimo de J*/k), integrado no dominio do escoamento
    void Print() const;        // parametros e conferencias (J* e div f = 0)

private:
    REAL ThetaMax(REAL r) const; // angulo do fim do arco de raio r (face ou terreno do pe)
    REAL H2(REAL t) const { return std::cosh(fA * t); }
    REAL DH2(REAL t) const { return fA * std::sinh(fA * t); }
    REAL fGw, fHw, fH, fBeta, fOmega, fRw, fR, fRe, fXO, fYO;
    REAL fM = 0., fA = 0., fC = 0., fD = 0., fK0 = 0.; // m, sqrt(m), C, D e K0 da zona 1
};

#endif
