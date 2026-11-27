#!/usr/bin/env python3
"""Richards 1D (so agua, p_a = 0), coluna vertical, z para baixo, silte do Cho (2016). Newton com busca linear.
FV nos nos (no 0 na superficie, armazenamento concentrado no no, como o FD do FLAC).
Topo: chuva 3 k_s ate empocar, depois h = 0. Base: h = h0 (-30 kPa).
K na interface: arit (media de K), satm (K da media de S_e: zona do FLAC), gauss (2 pontos de Gauss de h linear: MEF).
Uso: r1d.py dz media [L] [dtmax_h]"""
import sys, numpy as np
from scipy.linalg import solve_banded
P0, a, b, por, ks = 4.905, 0.336, 0.5, 0.5, 3.e-6 * 3600.
gw = 9.81; nvg = 1. / (1. - a); hP0 = P0 / gw; h0 = -30. / gw
def Se(h): return (1. + (np.maximum(-h, 0.) / hP0) ** nvg) ** (-a)
def Kr(s):
    s = np.clip(s, 1e-14, 1.); return s ** b * (1. - np.maximum(1. - s ** (1. / a), 0.) ** a) ** 2
EPS = 1e-3  # fecho linear de k_r entre S_e = 1 - EPS e 1 (como o MualemTabelado do NeoPZ): derivada finita
def Krf(s):
    s = np.asarray(s, dtype=float)
    k1 = Kr(1. - EPS)
    return np.where(s > 1. - EPS, k1 + (1. - k1) * (s - (1. - EPS)) / EPS, Kr(s))
def K(h): return ks * Krf(Se(h))
def th(h): return por * Se(h)
dz = float(sys.argv[1]); media = sys.argv[2]; L = float(sys.argv[3]) if len(sys.argv) > 3 else 10.
dtmax = float(sys.argv[4]) if len(sys.argv) > 4 else 0.05
N = int(round(L / dz)) + 1; z = np.linspace(0., L, N); dz = z[1] - z[0]
V = np.full(N, dz); V[0] = V[-1] = dz / 2.; qr = 3. * ks
gp = 0.5 / np.sqrt(3.)
def Kf2(hi, hj):
    if media == 'arit': return 0.5 * (K(hi) + K(hj))
    if media == 'satm': return ks * Krf(0.5 * (Se(hi) + Se(hj)))
    if media == 'upw': return np.where(hi - hj + dz > 0., K(hi), K(hj))  # o no de montante (H = h - z, z para baixo)
    if media == 'max': return np.maximum(K(hi), K(hj))
    if media == 'sup': return K(hi)  # montante fixo: o no de cima (o fluxo e todo para baixo na infiltracao)
    if media == 'geo': return np.sqrt(K(hi) * K(hj))
    return 0.5 * (K((0.5 + gp) * hi + (0.5 - gp) * hj) + K((0.5 - gp) * hi + (0.5 + gp) * hj))
def residuo(h, thn, dtt, emp):
    Kf = Kf2(h[:-1], h[1:]); q = Kf * ((h[:-1] - h[1:]) / dz + 1.)
    R = V * (th(h) - thn) / dtt
    R[:-1] += q; R[1:] -= q
    if emp: R[0] = h[0]
    else: R[0] -= qr
    R[-1] = h[-1] - h0
    return R
def jac(h, thn, dtt, emp):
    e = 1e-7
    hp = h + e
    Kf = Kf2(h[:-1], h[1:]); g = (h[:-1] - h[1:]) / dz + 1.
    dKi = (Kf2(hp[:-1], h[1:]) - Kf) / e; dKj = (Kf2(h[:-1], hp[1:]) - Kf) / e
    dqi = Kf / dz + dKi * g; dqj = -Kf / dz + dKj * g
    diag = V * (th(hp) - th(h)) / e / dtt
    diag[:-1] += dqi; diag[1:] -= dqj
    up = np.zeros(N); lo = np.zeros(N)
    up[1:] = dqj          # dR_i/dh_{i+1}
    lo[:-1] = -dqi        # dR_{i+1}/dh_i
    if emp: diag[0] = 1.; up[1] = 0.
    diag[-1] = 1.; lo[-2] = 0.
    return np.vstack([up, diag, lo])
def passo(h, dtt, emp):
    thn = th(h); hk = h.copy(); R = residuo(hk, thn, dtt, emp); nr0 = np.abs(R).max()
    for it in range(30):
        dh = solve_banded((1, 1), jac(hk, thn, dtt, emp), -R)
        lam = 1.; nr = np.abs(R).max()
        while lam > 1e-3:
            hn = hk + lam * dh; Rn = residuo(hn, thn, dtt, emp)
            if np.abs(Rn).max() < nr * (1 - 1e-4 * lam) or np.abs(Rn).max() < 1e-13: break
            lam /= 2.
        hk, R = hn, Rn
        if np.abs(R).max() < 1e-9 and np.abs(lam * dh).max() < 1e-6: return hk, it
    return None, it
h = np.full(N, h0); t, dt = 0., 1e-4
saidas = [2, 4, 8, 16, 32, 64, 96, 128, 160]; iout = 0; pond = None; Iacc = 0.; res = []
def frente(h):
    i = np.where(h >= -15. / gw)[0]
    if len(i) == 0: return 0.
    j = i.max()
    return z[j] + (h[j] + 15. / gw) / (h[j] - h[j + 1]) * dz if j + 1 < N else z[-1]
while iout < len(saidas):
    if saidas[iout] - t < 1e-7:
        t = saidas[iout]; res.append((saidas[iout], frente(h), Iacc)); iout += 1; continue
    tn = min(t + dt, saidas[iout]); dtt = tn - t
    hk, it = passo(h, dtt, pond is not None)
    if hk is not None and pond is None and hk[0] > 0.:
        hk, it = passo(h, dtt, True)
        if hk is not None: pond = tn
    if hk is None:
        dt /= 2.
        if dt < 1e-9: raise SystemExit(f"dt pequeno demais em t={t}")
        continue
    Iacc += np.sum(V * (th(hk) - th(h))); h = hk; t = tn
    dt = min(dt * (1.3 if it < 5 else 0.8), dtmax)
print(f"dz={dz:.4f} {media:5s} empoca={pond:.2f}h " + " ".join(f"{T}h:{f:.2f}" for T, f, I in res) + "  I(32h)=" + f"{[I for T,f,I in res if T==32][0]:.3f}", flush=True)
