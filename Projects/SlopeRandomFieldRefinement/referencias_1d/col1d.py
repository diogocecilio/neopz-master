# col1d.py
import numpy as np
from scipy.linalg import solve_banded
GW = 9.81
class Solo:   # VG-Mualem; psi_e_kpa > 0: Ippisch et al. 2006 (k_r normalizado em psi_e, saturado acima)
    def __init__(s, ts=0.358, tr=3.58e-4, a_kpa=0.149, n=1.289, ks=3e-6, psi_e_kpa=0.0, Ss=1.6e-6):
        s.ts, s.tr, s.n, s.m = ts, tr, n, 1 - 1 / n
        s.a = a_kpa * GW; s.ks = ks; s.he = psi_e_kpa / GW; s.Ss = Ss
        s.Sc = (1 + (s.a * s.he) ** s.n) ** (-s.m) if s.he > 0 else 1.0
        s.krc = s._krm(s.Sc) if s.he > 0 else 1.0
    def _krm(s, Se):
        Se = np.clip(Se, 1e-300, 1.0); f = 1 - (1 - Se ** (1 / s.m)) ** s.m
        return np.sqrt(Se) * f * f
    def Se(s, h):
        psi = np.maximum(-h, 0.0); Sv = (1 + (s.a * psi) ** s.n) ** (-s.m)
        return np.where(psi > s.he, Sv / s.Sc, 1.0)
    def theta(s, h): return s.tr + (s.ts - s.tr) * s.Se(h)
    def C(s, h):
        psi = np.maximum(-h, 0.0); x = (s.a * psi) ** s.n
        d = -s.m * s.n * s.a ** s.n * psi ** (s.n - 1) * (1 + x) ** (-s.m - 1) / s.Sc
        return np.where(psi > s.he, -(s.ts - s.tr) * d, 0.0)
    def kr(s, h):
        psi = np.maximum(-h, 0.0); Sv = (1 + (s.a * psi) ** s.n) ** (-s.m)
        return np.where(psi > s.he, np.minimum(s._krm(Sv) / s.krc, 1.0), 1.0)
    def dkr(s, h):   # dk_r/dh analitica: ~ psi^(n-2) perto de psi = 0
        psi = np.maximum(-h, 1e-300); y = (s.a * psi) ** s.n; Sv = (1 + y) ** (-s.m)
        q = 1 / (1 + y); omq = y / (1 + y); f = 1 - omq ** s.m
        dSv = -s.m * s.n * s.a ** s.n * psi ** (s.n - 1) * (1 + y) ** (-s.m - 1)
        with np.errstate(all='ignore'):
            d = -(0.5 / np.sqrt(Sv) * f * f + np.sqrt(Sv) * 2 * f * omq ** (s.m - 1) * q / Sv) * dSv / s.krc
        d = np.nan_to_num(d, nan=0.0, posinf=1e30, neginf=-1e30)
        return np.where(-h > s.he, d, 0.0)

class SoloTab(Solo):  # k_r de Mualem tabelado em S_e (ln k_r linear), ultimo ponto em 1 - eps
    def __init__(s, eps=1e-3, **kw):
        super().__init__(**kw)
        u = np.logspace(np.log10(eps), 0, 200)[::-1]
        St = np.concatenate([1 - u, [1.0]]); St[0] = max(St[0], 1e-6)
        s.St = St; s.lnk = np.log(np.maximum(s._krm(np.minimum(St, 1 - 1e-15)), 1e-300)); s.lnk[-1] = 0.0
        s.g = np.diff(s.lnk) / np.diff(s.St)
    def kr(s, h): return np.exp(np.interp(s.Se(h), s.St, s.lnk))
    def dkr(s, h):
        Se = s.Se(h); i = np.clip(np.searchsorted(s.St, Se) - 1, 0, len(s.g) - 1)
        return np.where(Se < 1.0, s.kr(h) * s.g[i] * s.C(h) / (s.ts - s.tr), 0.0)

class Coluna:   # y = elevacao (0 na base); h = carga de pressao (m); km = k_s(y)/k_s por celula
    def __init__(s, solo, H=2.0, N=200, kmult=None, r=3e-6, h0=-20 / GW):
        s.so, s.H, s.N, s.dz = solo, H, N, H / N
        s.y = (np.arange(N) + 0.5) * s.dz
        s.km = np.ones(N) if kmult is None else np.asarray(kmult)
        s.r, s.h, s.t, s.face = r, np.full(N, h0), 0.0, False
        s.stats = dict(it=0, fail=0)
    def K(s, h): return s.so.ks * s.km * s.so.kr(h)
    def solve(s, hn, dt, face, metodo, tol, maxit):
        so, dz, N = s.so, s.dz, s.N; h = hn.copy()
        for it in range(1, maxit + 1):
            K = s.K(h); Kf = 0.5 * (K[:-1] + K[1:]); g = (h[1:] - h[:-1]) / dz + 1.0
            R = (so.theta(h) - so.theta(hn) + so.Ss * (h - hn)) * dz / dt
            R[:-1] += -Kf * g; R[1:] -= -Kf * g
            gtop = -h[-1] / (dz / 2) + 1.0
            R[-1] -= K[-1] * gtop if face else s.r
            main = so.C(h) * dz / dt + so.Ss * dz / dt; up = np.zeros(N); lo = np.zeros(N); a = Kf / dz
            main[:-1] += a; up[1:] -= a; main[1:] += a; lo[:-1] -= a
            if face: main[-1] += K[-1] / (dz / 2)
            if metodo == 'newton':
                dK = so.ks * s.km * so.dkr(h)
                main[:-1] += -0.5 * dK[:-1] * g; up[1:] += -0.5 * dK[1:] * g
                main[1:] += 0.5 * dK[1:] * g; lo[:-1] += 0.5 * dK[:-1] * g
                if face: main[-1] -= dK[-1] * gtop
            dh = solve_banded((1, 1), np.vstack([up, main, lo]), -R)
            if not np.all(np.isfinite(dh)): return None, it
            h = h + dh
            if np.abs(dh).max() * GW < tol: return h, it      # |dp| em kPa, como o tol_p do NeoPZ
        return None, maxit
    def step(s, dt, metodo='newton', tol=1e-6, maxit=40):  # troca fluxo <-> face como o RainTransient
        hn = s.h.copy(); face = s.face; its = 0
        for k in range(40):
            h, it = s.solve(hn, dt, face, metodo, tol, maxit); its += it
            if h is None: return False, its
            K = s.K(h)
            if face: nova = not (K[-1] * (-h[-1] / (s.dz / 2) + 1.0) > s.r * (1 + 1e-9))
            else: nova = h[-1] + (s.dz / 2) * (s.r / K[-1] - 1.0) > 1e-6 / GW
            if nova == face:
                s.h, s.face, s.t = h, face, s.t + dt; return True, its
            face = nova
        return False, its
    def advance_to(s, T, dt=360.0, metodo='newton', dtmin=1.0, **kw):
        while s.t < T - 1e-9:
            pilha = [min(dt, T - s.t)]
            while pilha:
                d = pilha.pop(); ok, its = s.step(d, metodo, **kw); s.stats['it'] += its
                if not ok:
                    s.stats['fail'] += 1
                    if d / 2 < dtmin: raise RuntimeError(f'nao convergiu em t={s.t/3600:.4f} h')
                    pilha += [d / 2, d / 2]

def fs_perfil(col, c=5.0, phi=32.0, beta=40.0, gd=16.0):   # Eqs. 5, 7 e 9 de Cho (2014)
    so = col.so; b = np.radians(beta); tp = np.tan(np.radians(phi))
    z = col.H - col.y; o = np.argsort(z)
    gt = gd + so.theta(col.h)[o] * GW
    W = np.cumsum(gt * col.dz) - 0.5 * gt * col.dz
    uw = GW * col.h[o]; sig = np.where(uw < 0, so.Se(col.h)[o] * uw, uw)
    FS = (c + (W * np.cos(b) ** 2 - sig) * tp) / (W * np.sin(b) * np.cos(b))
    return z[o], FS, col.h[o]
