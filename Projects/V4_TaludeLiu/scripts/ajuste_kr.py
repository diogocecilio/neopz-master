#!/usr/bin/env python3
# ajuste_kr.py: a funcao de permeabilidade de Fredlund, Xing e Huang (1994) (FXH, eq. 9 de Liu et al. 2024) com a curva
# de retencao de Fredlund e Xing (eqs. 7-8: a = 5 kPa, n = 2, m = 1, theta_s = 0.4, psi_r = 3000 kPa) e o ajuste da lei
# de potencia k_r = S^p da biblioteca (S da curva de Fredlund-Xing sem C(psi), como no TPZUnsaturatedHydraulics).
# psi_aev (o limite inferior da integral da eq. 9) nao e dado no artigo: e estimado pela construcao grafica usual
# (tangente no ponto de inflexao da curva theta x log psi, interceptada com theta = theta_s) e o script mostra a
# sensibilidade a ele. O ajuste de p e por minimos quadrados de log10(k_r) entre psi_aev e PSI_MAX (padrao 20 kPa, a
# faixa de succao do talude). Tambem a succao do patamar de infiltracao permanente (gradiente unitario: k_r(psi) = q/k_s)
# para q/k_s = 0.01, 0.05 e 0.1. Uso: python3 ajuste_kr.py [PSI_MAX] [PSI_AEV]
import sys
import numpy as np

a, n, m, ts, psir = 5.0, 2.0, 1.0, 0.4, 3000.0

def C(psi):
    return 1.0 - np.log(1.0 + psi / psir) / np.log(1.0 + 1.0e6 / psir)

def theta(psi, corr=True):
    t = ts / np.log(np.e + (psi / a) ** n) ** m
    return t * C(psi) if corr else t

def dtheta(psi, corr=True, h=1e-6):
    return (theta(psi * (1 + h), corr) - theta(psi * (1 - h), corr)) / (2 * h * psi)

def aev_tangente():
    # inflexao de theta(log10 psi): maximo de |d theta / d log10 psi|
    lp = np.linspace(-2, 3, 200001)
    psi = 10 ** lp
    d = np.gradient(theta(psi), lp)
    i = np.argmin(d)
    # tangente: theta = theta_i + d_i (lp - lp_i) = theta_s  ->  lp = lp_i + (theta_s - theta_i)/d_i
    return 10 ** (lp[i] + (ts - theta(psi[i])) / d[i]), psi[i]

def kr_fxh(psis, aev, npts=20000):
    b = np.log(1.0e6)
    def integral(lo, tref):
        y = np.linspace(lo, b, npts)
        ey = np.exp(y)
        f = (theta(ey) - tref) / ey * dtheta(ey)
        return np.trapz(f, y) if hasattr(np, 'trapz') else np.trapezoid(f, y)
    den = integral(np.log(aev), ts)
    out = []
    for p in np.atleast_1d(psis):
        out.append(1.0 if p <= aev else integral(np.log(p), theta(p)) / den)
    return np.array(out)

def S_fx(psi):  # retencao da biblioteca: C(psi) = 1, S = theta/theta_s
    return 1.0 / np.log(np.e + (psi / a) ** n) ** m

def patamar(kfun, alvo, lo=1e-3, hi=1e3):
    for _ in range(200):
        mid = np.sqrt(lo * hi)
        if kfun(mid) > alvo: lo = mid
        else: hi = mid
    return np.sqrt(lo * hi)

psimax = float(sys.argv[1]) if len(sys.argv) > 1 else 20.0
aev_t, psi_inf = aev_tangente()
aevs = [float(sys.argv[2])] if len(sys.argv) > 2 else [aev_t, 1.0, 0.1, a]
print(f"psi_aev pela tangente no ponto de inflexao (psi = {psi_inf:.3f} kPa): {aev_t:.3f} kPa")
print(f"C(psi) em 20 kPa: {C(20.):.5f}")
for aev in aevs:
    psi = np.geomspace(max(aev, 1e-3) * 1.0001, psimax, 300)
    kr = kr_fxh(psi, aev)
    lS = np.log10(S_fx(psi))
    p = np.sum(lS * np.log10(kr)) / np.sum(lS * lS)  # minimos quadrados de log10 kr = p log10 S (sem intercepto)
    print(f"\npsi_aev = {aev:.3f} kPa: p ajustado = {p:.3f} (log10 kr, {psi[0]:.2f} a {psimax:.0f} kPa)")
    for qk in (0.01, 0.05, 0.1):
        pf = patamar(lambda s: kr_fxh([s], aev)[0], qk, max(aev, 1e-3), 1e3)
        pp = patamar(lambda s: S_fx(s) ** p, qk)
        print(f"  q/ks = {qk:4.2f}: patamar FXH {pf:6.2f} kPa, lei de potencia p = {p:.2f}: {pp:6.2f} kPa")
    print("  psi (kPa)   kr FXH      S^p        S_FX")
    for s in (1, 2, 3, 5, 7, 10, 12.5, 15, 20):
        if s <= aev: continue
        print(f"  {s:8.2f}  {kr_fxh([s], aev)[0]:.4e}  {S_fx(s) ** p:.4e}  {S_fx(s):.4f}")
