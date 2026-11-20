"""Campo de velocidades de filtracao semianalitico v'_opt (Secao 4.1 do
artigo, Eqs. 29-40), usado como forca de percolacao -grad u ~ K^-1 . v'_opt
na abordagem cinematica. So permeabilidade isotropica (alfa = 1).

Coordenadas polares em O (aresta da crista), como na Fig. 3 e no resto do
codigo: x = -r cos(theta), y = r sin(theta) (y para baixo); theta = 0 e a
crista, theta = pi - beta e a face. e_r = (-cos, sin), e_theta = (sin, cos).

Classe (Eq. 29), com R_w = h_w/sin(beta), R = H/sin(beta) e
R_e = sqrt(H^2 + (L_m + H/tan(beta))^2), L_m = 10 H:
  zona 1 (r < R_w)     : v = -h1(r) h2'(theta) e_r + (r h1)'(r) h2(theta) e_theta
  zona 2 (R_w <= r < R): v = h3(r) e_theta
  zona 3 (R <= r < R_e): v = h4(r) e_theta  (theta ate o terreno do pe)
  r >= R_e             : v = 0.
Funcional (Eq. 23) com u^d = 0 na crista, -gamma_w y na face exposta e
-gamma_w h_w na face submersa e no terreno do pe.

Deducao (alfa = 1, K = k I), conferida com as Eqs. 32, 33, 34, 35 e 40:
  * h3 = k gamma_w h_w/(A r), h4 = 2 k gamma_w h_w/(B r)  (Eqs. 32-35);
  * zona 1, com w = r h1: para h2 dado, a equacao de Euler-Lagrange em w e
    D (r w')' - C w/r = q,  q = k gamma_w sin(beta) h2(pi - beta),
    C = int h2'^2, D = int h2^2 (Eq. 36 com alfa = 1), w(0) = w(R_w) = 0:
        w = q r/(D - C) (1 - (r/R_w)^(m-1)),  m = sqrt(C/D);
    para h1 dado, a de h2 e h2'' = (b/a) h2 com h2'(0) = 0 (Eqs. 37-38),
    entao h2 = cosh(kappa theta) e o otimo conjunto e o minimo em kappa de
        J1(kappa) = -k gamma_w^2 h_w^2/4 (h2(pi-beta)/(sqrt C + sqrt D))^2,
    que e o primeiro termo da Eq. 40. O minimo em kappa faz o papel da
    Eq. 39.
  * Componentes na zona 1 (h1 = w/r):
        v_r     = -q/(D-C) (1 - (r/R_w)^(m-1)) h2'(theta)
        v_theta =  q/(D-C) (1 - m (r/R_w)^(m-1)) h2(theta)
    A Eq. 31 impressa traz (r/R_w)^(C/D - 1) C/D no termo de e_theta e o
    fator 1/(C - D) na frente; com isso div v != 0 e o fluxo na face exposta
    entraria no talude. A forma acima e a que sai da classe (29) e do
    funcional (23); o primeiro termo da Eq. 40 confere com ela.
"""
from __future__ import annotations

import math

import numpy as np
from scipy import optimize

from cinematico import GAMMA_W


class CampoVopt:
    """c, phi uniformes e forca de percolacao -grad u = v'_opt/k (alfa = 1).

    Interface de campo de cinematico.potencias: grad_u(x, y) devolve
    grad u = -v'/k. Nao ha potencial u (v'_opt nao e um gradiente), entao a
    conferencia de P_u pela forma de contorno nao se aplica."""
    uniforme = True
    tem_potencial = False
    tem_u = True

    def __init__(self, c, phi, gamma_ef, H, beta, hw=None, Lm=None, alfa=1.0):
        if alfa != 1.0:
            raise NotImplementedError("so alfa = 1 (permeabilidade isotropica)")
        self.c0, self.phi0, self.gamma_ef = float(c), float(phi), float(gamma_ef)
        self.phi_ref = self.phi0
        self.H, self.beta = float(H), float(beta)
        self.hw = self.H if hw is None else float(hw)
        Lm = 10.0 * self.H if Lm is None else float(Lm)
        sb = math.sin(self.beta)
        self.Th = math.pi - self.beta
        self.Rw = self.hw / sb
        self.R = self.H / sb
        self.Re = math.hypot(self.H, Lm + self.H / math.tan(self.beta))
        self.A = math.pi - self.beta          # A(1) = (alfa+1)(pi-beta)/2 (Eq. 34)
        Th = self.Th

        def CD(kap):
            if kap < 1e-8:
                return kap * kap * Th ** 3 / 3.0, Th          # limite kappa -> 0
            s2 = math.sinh(2 * kap * Th) / (4 * kap)
            return kap * kap * (s2 - Th / 2), s2 + Th / 2

        def Jn(lk):        # J1 normalizado: -(h2(Th)/(sqrt C + sqrt D))^2, h2 = cosh
            kap = math.exp(lk)
            C, D = CD(kap)
            return -(math.cosh(kap * Th) / (math.sqrt(C) + math.sqrt(D))) ** 2
        r = optimize.minimize_scalar(Jn, bounds=(math.log(1e-4), math.log(20.0)),
                                     method="bounded", options=dict(xatol=1e-10))
        self.kappa = math.exp(r.x)
        self.C, self.D = CD(self.kappa)
        self.m = math.sqrt(self.C / self.D)
        self.J1n = r.fun
        # forca de percolacao por k: q/k = gamma_w sin(beta) h2(Th)
        self.p = GAMMA_W * sb * math.cosh(self.kappa * Th) / (self.D - self.C)

    def B(self, r):
        r = np.asarray(r, float)
        return 2.0 * (math.pi - np.arcsin(np.minimum(self.H / r, 1.0)))   # alfa = 1

    def J_normalizado(self):
        """J*(v'_opt)/(k gamma_w^2 H^2) (Eq. 40), para comparar com a Fig. 5."""
        from scipy import integrate
        I3 = integrate.quad(lambda r: 4.0 / (self.B(r) * r), self.R, self.Re, limit=200)[0]
        J = -(self.hw ** 2) / 4.0 * (-self.J1n - 2.0 / self.A * math.log(self.hw / self.H) + I3)
        return J / self.H ** 2

    def c(self, x, y):
        return np.full(np.shape(x), self.c0)

    def phi(self, x, y):
        return np.full(np.shape(x), self.phi0)

    def phi_max_arco(self, xs, ys, flecha=0.0):
        return self.phi0

    def v_sobre_k(self, x, y):
        """v'_opt/k em coordenadas cartesianas (x, y para baixo)."""
        x = np.asarray(x, float)
        y = np.asarray(y, float)
        r = np.hypot(x, y)
        th = np.arctan2(y, -x)
        vr = np.zeros_like(r)
        vt = np.zeros_like(r)
        Th, kap = self.Th, self.kappa
        z1 = (r > 0) & (r < self.Rw) & (th >= 0) & (th <= Th)
        rr = r[z1] / self.Rw
        e = rr ** (self.m - 1.0)
        vr[z1] = -self.p * (1.0 - e) * kap * np.sinh(kap * th[z1])
        vt[z1] = self.p * (1.0 - self.m * e) * np.cosh(kap * th[z1])
        z2 = (r >= self.Rw) & (r < self.R) & (th >= 0) & (th <= Th)
        vt[z2] = GAMMA_W * self.hw / (self.A * r[z2])
        z3 = (r >= self.R) & (r < self.Re) & (th >= 0)
        if np.any(z3):
            thg = math.pi - np.arcsin(np.minimum(self.H / r[z3], 1.0))
            ok = th[z3] <= thg
            idx = np.nonzero(z3)[0][ok]
            vt[idx] = 2.0 * GAMMA_W * self.hw / (self.B(r[idx]) * r[idx])
        s, co = np.sin(th), np.cos(th)
        vx = -vr * co + vt * s
        vy = vr * s + vt * co
        return vx, vy

    def grad_u(self, x, y):
        vx, vy = self.v_sobre_k(x, y)
        return -vx, -vy

    def div_v(self, x, y, h=1e-6):
        """div(v/k) por diferencas finitas (conferencia: deve ser ~0)."""
        a = self.v_sobre_k(x + h, y)[0] - self.v_sobre_k(x - h, y)[0]
        b = self.v_sobre_k(x, y + h)[1] - self.v_sobre_k(x, y - h)[1]
        return (a + b) / (2 * h)
