"""Referência 1D da dessaturação da coluna de areia (Liakopoulos, 1965), com os dados de
Callari e Abati (2009), Seção 4.1.

Modelo (o mesmo do código, reduzido ao caso oedométrico):
  - ar passivo (p_a = 0), Bishop com chi = S_w, logo p_B = S_w p_w;
  - esqueleto elástico em oedômetro: com peso constante e topo livre, a tensão total vertical
    não muda, e eps_v = alpha p_B / M, com M o módulo oedométrico;
  - equação da água com armazenamento secante concentrado nos nós (n dS_w), compressibilidade
    da água e acoplamento S_w alpha d(eps_v); fluxo de Darcy com k_rw(S_w) (média dos nós);
  - Euler implícito e Newton (jacobiana tridiagonal por diferenças finitas).
Variantes: saturado (S_w = 1, k_rw = 1), coussy (dp_f = S_w dp_w, Eq. 47 do artigo) e
taxa (n C_s dp_w no lugar de n dS_w, Eq. 37 do artigo).
Unidades: kN, m, s, kPa. y = 0 na base, y = H no topo.

Uso:  python coluna_liakopoulos_1d.py            (imprime a tabela de referência)
"""
import numpy as np
from scipy.linalg import solve_banded

# ---------------- dados (Callari e Abati 2009, Tabela 1 e Eq. 58) ----------------
E, nu, alpha = 1300.0, 0.4, 1.0          # kPa, -, -   (grãos incompressíveis: 1/N = 0)
n, Kw, gw = 0.2975, 2.0e6, 9.81          # -, kPa, kN/m3
ks = gw*4.5e-7                           # m/s  (k_w^sat = 4,5e-7 m2/(kPa s) = kappa/mu_w)
cL = 1.9722e-11*1000.0**2.4279           # S_w = 1 - cL pc^bL, pc em kPa (no artigo, pc em Pa)
bL, dL, eL = 2.4279, 2.207, 1.0121       # k_rw = 1 - dL (1 - S_w)^eL
M = E*(1 - nu)/((1 + nu)*(1 - 2*nu))     # módulo oedométrico, kPa
H = 1.0                                  # altura da coluna, m


def leis(saturado):
    if saturado:
        return (lambda p: np.ones_like(p)), (lambda s: np.ones_like(s))
    Sw = lambda p: np.where(p < 0, 1 - cL*np.abs(np.minimum(p, 0))**bL, 1.0)
    krw = lambda s: np.maximum(1 - dL*np.clip(1 - s, 0, None)**eL, 1e-12)
    return Sw, krw


def passo_referencia(t):
    return 1.0 if t < 60 else (5.0 if t < 600 else (10.0 if t < 2000 else 20.0))


def passo_artigo(t):
    """Passos sugeridos pelos marcadores das Figs. 3a e 4a do artigo."""
    return 10.0 if t < 300 else (100.0 if t < 2400 else 1000.0)


def resolve(N=201, saturado=False, coussy=False, taxa=False, dt_fn=passo_referencia,
            t_fim=7500.0, t_perfis=()):
    """Retorna y, histórico {t, q_b (m/s), w_topo (m)} e perfis {t: (p_w, u_y)}."""
    Sw, krw = leis(saturado)
    y = np.linspace(0.0, H, N); h = y[1] - y[0]
    V = np.full(N, h); V[0] = V[-1] = h/2
    p = np.zeros(N); S = Sw(p); eps = np.zeros(N); t = 0.0
    hist = {'t': [0.0], 'q': [ks], 'w': [0.0]}; perfis = {}

    def residuo(pn, p_, Sn, epsn, dt):
        S_ = Sw(p_)
        eps_ = epsn + alpha*S_*(p_ - pn)/M if coussy else alpha*S_*p_/M
        if taxa and not saturado:
            dS = cL*bL*np.maximum(-p_, 0.0)**(bL - 1)*(p_ < 0)*(p_ - pn)
        else:
            dS = S_ - Sn
        R = V*(n*dS + n*S_*(p_ - pn)/Kw + S_*alpha*(eps_ - epsn))
        kr = krw(S_); q = -(ks/gw)*0.5*(kr[:-1] + kr[1:])*((p_[1:] - p_[:-1])/h + gw)
        R[:-1] += dt*q; R[1:] -= dt*q      # q: fluxo para cima em cada elemento
        return R, q, eps_

    while t < t_fim - 1e-9:
        dt = min(dt_fn(t), t_fim - t)
        pn, Sn, epsn = p.copy(), S.copy(), eps.copy(); p_ = p.copy()
        for _ in range(50):
            R, q, _ = residuo(pn, p_, Sn, epsn, dt); R[0] = 0.0     # p_w = 0 na base
            ab = np.zeros((3, N)); d = 1e-7
            for c in range(3):
                pp = p_.copy(); idx = np.arange(c, N, 3); pp[idx] += d
                dR = (residuo(pn, pp, Sn, epsn, dt)[0] - R)/d
                for j in idx:
                    for i in (j - 1, j, j + 1):
                        if 0 <= i < N:
                            ab[1 + i - j, j] = dR[i]
            ab[:, 0] = 0.0; ab[1, 0] = 1.0; ab[0, 1] = 0.0
            dp = solve_banded((1, 1), ab, -R); p_ += dp
            if np.max(np.abs(dp)) < 1e-10:
                break
        else:
            raise RuntimeError(f"Newton não convergiu em t = {t}")
        R, q, eps = residuo(pn, p_, Sn, epsn, dt)
        p, S, t = p_, Sw(p_), t + dt
        hist['t'].append(t); hist['q'].append(-q[0]); hist['w'].append(np.trapezoid(eps, y))
        for tp in t_perfis:
            if abs(t - tp) < 1e-6:
                u = np.concatenate([[0.0], np.cumsum(0.5*(eps[1:] + eps[:-1])*h)])
                perfis[tp] = (p.copy(), u)
    return y, {k: np.array(v) for k, v in hist.items()}, perfis


def saturado_analitico(t, kmax=4000):
    """Terzaghi com excesso inicial gamma_w*y, drenagem na base e topo impermeável."""
    cv = (ks/gw)/(1/M + n/Kw); T = cv*t/H**2
    m = (2*np.arange(kmax) + 1)*np.pi/2; s = (-1.0)**np.arange(kmax)
    w = gw*H**2/(2*M) - (2*gw*H**2/M)*np.sum(s*np.exp(-m**2*T)/m**3)
    q = 2*ks*np.sum(s*np.exp(-m**2*T)/m) if t > 0 else ks
    return -w, q


if __name__ == '__main__':
    T = (300, 600, 1200, 1800, 3600, 6400, 7500)
    y, h, perfis = resolve(t_perfis=T)
    print(f"M = {M:.1f} kPa, k_s = {ks:.4e} m/s, c_L = {cL:.4e} kPa^-bL")
    print("   t (s)  q_b (1e-6 m/s)  p_w(0,5 m)  p_w(topo)  w (mm)  | saturada: q_b    w (mm)")
    for t in T:
        p, _ = perfis[t]; wa, qa = saturado_analitico(t)
        print(f"  {t:6d}  {1e6*np.interp(t, h['t'], h['q']):10.3f}  {np.interp(0.5, y, p):10.2f}"
              f"  {p[-1]:9.2f}  {1000*np.interp(t, h['t'], h['w']):7.3f}  |  {1e6*qa:10.3f}  {1000*wa:7.3f}")
