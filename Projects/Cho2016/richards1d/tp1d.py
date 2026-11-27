#!/usr/bin/env python3
"""tp1d.py: fluxo bifasico agua-ar 1D (coluna vertical, z para baixo) do Cho (2016), independente do NeoPZ.
Incognitas por no: p_w e p_a (kPa). Volumes finitos nos nos (armazenamento concentrado no no, como o FD do FLAC):
  agua: n dS_w/dt + n S_w/K_w dp_w/dt + dq_w/dz = 0,  q_w = (k_s k_rw / g_w) ((p_w,i - p_w,i+1)/dz + g_w)
  ar:   n dS_a/dt + n S_a/K_a dp_a/dt + dq_a/dz = 0,  q_a = (k_s mu k_ra / g_w) ((p_a,i - p_a,i+1)/dz + g_a)
  K_a = p_atm + p_a (gas ideal, como o NeoPZ; o FLAC usa 100 kPa constante), S_w = S_e(p_a - p_w) de van Genuchten.
  k_rw de Mualem com o fecho linear entre S_e = 1 - 1e-3 e 1 (como o MualemTabelado); k_ra de Lenhard-Parker com piso.
Topo: p_a = 0; a chuva 3 k_s ate empocar, depois p_w = 0. Base: sem fluxo de agua e de ar (L = 10 m, como a coluna).
k_r na interface: gauss (2 pontos de Gauss de p_w e p_a lineares: como o MEF), arit (media), sup (no de montante de
cada fase, pelo sinal do gradiente de potencial: o que reproduz o FLAC no sp).
Uso: tp1d.py dz media [kramin=1e-4] [mu=56] [L=10] [dtmax=0.05] [solo=silte]
Saida: uma linha com a frente (p_w = -15 kPa) e p_a no fundo em 16, 32, 64, 96, 128 e 160 h"""
import sys, numpy as np
from scipy.linalg import solve_banded

args = dict(a.split('=') for a in sys.argv[3:])
silte = args.get('solo', 'silte') != 'areia'
P0, a, por, ks = (4.905, 0.336, 0.5, 3.e-6 * 3600.) if silte else (1.5, 0.445, 0.35, 1.e-5 * 3600.)
b = 0.5; gw = 9.81; ga = 1.25 * 9.81e-3; patm = 101.325; Kw = 2.e6
nvg = 1. / (1. - a)
kramin = float(args.get('kramin', 1e-4)); mu = float(args.get('mu', 56.)); L = float(args.get('L', 10.))
dtmax = float(args.get('dtmax', 0.05)); EPS = 1e-3
dz = float(sys.argv[1]); media = sys.argv[2]
N = int(round(L / dz)) + 1; z = np.linspace(0., L, N); dz = z[1] - z[0]
V = np.full(N, dz); V[0] = V[-1] = dz / 2.
qr = 3. * ks


def Se(s): return (1. + (np.maximum(s, 0.) / P0) ** nvg) ** (-a)


def Krw0(S):
    S = np.clip(S, 1e-14, 1.); return S ** b * (1. - np.maximum(1. - S ** (1. / a), 0.) ** a) ** 2


K1 = Krw0(1. - EPS)


def Krw(S): return np.where(S > 1. - EPS, K1 + (1. - K1) * (S - (1. - EPS)) / EPS, Krw0(S))


def Kra(S):
    S = np.clip(S, 0., 1.)
    return np.maximum((1. - S) ** 0.5 * np.maximum(1. - S ** (1. / a), 0.) ** (2. * a), kramin)


gp = 0.5 / np.sqrt(3.)


def faces(pw, pa):
    """k_rw e k_ra nas interfaces i+1/2"""
    if media == 'gauss':
        kw = ka = 0.
        for c in (0.5 + gp, 0.5 - gp):
            S = Se(c * (pa[:-1] - pw[:-1]) + (1 - c) * (pa[1:] - pw[1:]))
            kw = kw + 0.5 * Krw(S); ka = ka + 0.5 * Kra(S)
        return kw, ka
    S = Se(pa - pw)
    if media == 'arit':
        return 0.5 * (Krw(S[:-1]) + Krw(S[1:])), 0.5 * (Kra(S[:-1]) + Kra(S[1:]))
    # sup: montante pelo potencial de cada fase (z para baixo: potencial p - g z; desce se (p_i - p_i+1)/dz + g > 0)
    dw = (pw[:-1] - pw[1:]) / dz + gw; da = (pa[:-1] - pa[1:]) / dz + ga
    kw = np.where(dw > 0., Krw(S[:-1]), Krw(S[1:])); ka = np.where(da > 0., Kra(S[:-1]), Kra(S[1:]))
    return kw, ka


def residuo(x, xn, dtt, emp):
    pw, pa = x[0::2], x[1::2]; pwn, pan = xn[0::2], xn[1::2]
    S, Sn = Se(pa - pw), Se(pan - pwn)
    kw, ka = faces(pw, pa)
    qw = ks * kw / gw * ((pw[:-1] - pw[1:]) / dz + gw)
    qa = ks * mu * ka / gw * ((pa[:-1] - pa[1:]) / dz + ga)
    Rw = V * por * ((S - Sn) + S * (pw - pwn) / Kw) / dtt
    Ra = V * por * ((Sn - S) + (1. - S) * (pa - pan) / (patm + pa)) / dtt
    Rw[:-1] += qw; Rw[1:] -= qw; Ra[:-1] += qa; Ra[1:] -= qa
    if emp: Rw[0] = pw[0]
    else: Rw[0] -= qr
    Ra[0] = pa[0]
    R = np.empty(2 * N); R[0::2] = Rw; R[1::2] = Ra
    return R


def jacobiano(x, xn, dtt, emp, R0):
    """jacobiano em banda (3, 3) por diferencas com 6 cores (colunas que afetam linhas disjuntas)"""
    n = 2 * N; ab = np.zeros((7, n)); e = 1e-7
    for c in range(6):
        cols = np.arange(c, n, 6)
        xp = x.copy(); xp[cols] += e
        dR = (residuo(xp, xn, dtt, emp) - R0) / e
        for j in cols:
            lo, hi = max(0, j - 3), min(n, j + 4)
            ab[3 + np.arange(lo, hi) - j, j] = dR[lo:hi]
    return ab


def passo(x, dtt, emp):
    xn = x.copy(); xk = x.copy(); R = residuo(xk, xn, dtt, emp)
    for it in range(25):
        dx = solve_banded((3, 3), jacobiano(xk, xn, dtt, emp, R), -R)
        nr = np.abs(R).max(); lam = 1.
        while True:
            xt = xk + lam * dx; Rt = residuo(xt, xn, dtt, emp)
            if np.abs(Rt).max() < nr * (1. - 1e-4 * lam) or lam < 1e-3: break
            lam /= 2.
        xk, R = xt, Rt
        if not np.all(np.isfinite(xk)): return None, it
        if np.abs(R).max() < 1e-9 and np.abs(lam * dx).max() < 1e-5: return xk, it
    return None, it


x = np.zeros(2 * N); x[0::2] = -30.
t, dt = 0., 1e-4; saidas = [16, 32, 64, 96, 128, 160]; iout = 0; pond = None; res = []


def frente(pw):
    i = np.where(pw >= -15.)[0]
    if len(i) == 0: return 0.
    j = i.max()
    return z[j] + (pw[j] + 15.) / (pw[j] - pw[j + 1]) * dz if j + 1 < N else z[-1]


while iout < len(saidas):
    if saidas[iout] - t < 1e-7:
        t = saidas[iout]; res.append((t, frente(x[0::2]), x[1::2].max())); iout += 1; continue
    tn = min(t + dt, saidas[iout]); dtt = tn - t
    xk, it = passo(x, dtt, pond is not None)
    if xk is not None and pond is None and xk[0] > 0.:
        xk, it = passo(x, dtt, True)
        if xk is not None: pond = tn
    if xk is None:
        dt /= 2.
        if dt < 1e-9: print(f"dz={dz:.4f} {media} dt pequeno demais em t={t:.3f} " +
                            " ".join(f"{T:g}h:{f:.2f}/{p:.2f}" for T, f, p in res), flush=True); sys.exit(1)
        continue
    x = xk; t = tn; dt = min(dt * (1.3 if it < 5 else 0.8), dtmax)
print(f"dz={dz:.4f} {media:5s} kramin={kramin:g} mu={mu:g} empoca={pond:.2f}h frente/pa: " +
      " ".join(f"{T:g}h:{f:.2f}/{p:.2f}" for T, f, p in res), flush=True)
