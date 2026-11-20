#!/usr/bin/env python3
"""campos_cho.py: os campos c' e phi' do Cho2010 (Projects/Cho2010/cho2010.cpp) em Python, IGUAIS aos do MEF dados os xi.

O que o C++ faz (RandomField.cpp, KLSeparavel.h) e o que esta copia reproduz:
  * gerador: realizacao n, campo f (0 = c', 1 = phi'): klsep::Normal(semente * 0x9E3779B97F4A7C15 + 2 n + f).Vetor(M),
    xoshiro256** semeado por splitmix64 e Box-Muller (copia de Projects/Cho2010/diagnostico_ai.py, que confere com o
    log_peso das ISMC e com o campos.csv do Cho2010). Na ISMC: xi = sorteio + mu (mu so nos K primeiros termos).
  * campo gaussiano Z de cada campo:
      - 'analitico' (o MC direto do MEF): KL analitica do nucleo exp(-|dx|/20 - |dy|/2) no retangulo [0,30] x [0,15],
        os M = 2000 maiores produtos lambda_x lambda_y; Z = sum_k sqrt(lambda_k) xi_k f_k(x) g_k(y). Aqui na forma
        tensorial (grade): Z = Fx A Gy^T, com A[a, b] = soma de sqrt(lambda_k) xi_k dos termos k = (a, b).
      - 'nystrom' (as ISMC e o line sampling do MEF): base de Nystrom de 2000 modos do talude (15 000 nos, pesos w),
        Z(x) = sum_j C(x, x_j) w_j sum_k xi_k F_jk / sqrt(lambda_k). Aqui pre-calculado uma vez nos nos de solo da grade:
        Psi = C(grade, nos) W F Lambda^(-1/2) (float32, gravado em .npy) e Z = Psi xi.
  * correlacao cruzada: Z_phi <- r Z_c + sqrt(1 - r^2) Z_phi (depois da soma, antes da lognormal);
  * lognormal: v = exp(mu_ln + sigma_ln Z), sigma_ln = sqrt(ln(1 + CoV^2)), mu_ln = ln(media / sqrt(1 + CoV^2));
    c' (10 kPa, CoV 0,3), phi' (30 graus, CoV 0,2; lognormal em graus).

Coordenadas: o Cho2010 usa X para a direita (crista a esquerda, X em [0, 30]) e Y para cima (base Y = 0, crista Y = 15).
A analise limite usa as do Ceron: origem na aresta da crista (10, 15), x = X - 10 (para fora do talude), y = 15 - Y
(para baixo). A grade da analise limite e x = -10 + h i, y = h j (i = 0..30/h, j = 0..15/h)."""
from __future__ import annotations

import math
import os

import numpy as np

# ---------------------------------------------------------------------------------------------------------------- dados
MEDIA_C, COV_C = 10.0, 0.3
MEDIA_PHI, COV_PHI = 30.0, 0.2
LH, LV = 20.0, 2.0
DX, DY = 30.0, 15.0
M_PADRAO = 2000
X_CRISTA, Y_CRISTA = 10.0, 15.0        # aresta da crista nas coordenadas do Cho2010


def lognormal(media, cov, z):
    s = math.sqrt(math.log(1.0 + cov * cov))
    return np.exp(math.log(media / math.sqrt(1.0 + cov * cov)) + s * z)


# ------------------------------------------------- gerador (copia do klsep::Normal e do RandomField::Sample do C++)
def _u(x):
    return np.uint64(x)


def _rot(x, k):
    return (x << _u(k)) | (x >> _u(64 - k))


def normais(sementes, n):
    """os n primeiros normais de cada fluxo (sementes: uint64), (len(sementes), n)"""
    with np.errstate(over="ignore"):
        x = np.asarray(sementes).astype(np.uint64).copy()
        s = []
        for _ in range(4):
            x = x + _u(0x9E3779B97F4A7C15)
            z = (x ^ (x >> _u(30))) * _u(0xBF58476D1CE4E5B9)
            z = (z ^ (z >> _u(27))) * _u(0x94D049BB133111EB)
            s.append(z ^ (z >> _u(31)))
        s0, s1, s2, s3 = s
        out = np.empty((len(x), n))
        for j in range(0, n, 2):
            u = []
            for _ in range(2):
                r = _rot(s1 * _u(5), 7) * _u(9)
                t = s1 << _u(17)
                s2 = s2 ^ s0
                s3 = s3 ^ s1
                s1 = s1 ^ s2
                s0 = s0 ^ s3
                s2 = s2 ^ t
                s3 = _rot(s3, 45)
                u.append(((r >> _u(11)).astype(np.float64) + 0.5) * (1.0 / 9007199254740992.0))
            raio, ang = np.sqrt(-2.0 * np.log(u[0])), 2.0 * np.pi * u[1]
            out[:, j] = raio * np.cos(ang)
            if j + 1 < n:
                out[:, j + 1] = raio * np.sin(ang)
    return out


def xi_cpp(semente, reais, n=M_PADRAO, ncampos=2):
    """os xi do Cho2010 (sem deslocamento) das realizacoes dadas: (len(reais), ncampos, n)"""
    base = (int(semente) * 0x9E3779B97F4A7C15) & ((1 << 64) - 1)
    res = np.empty((len(reais), ncampos, n))
    for c in range(ncampos):
        sem = np.array([(base + ncampos * int(r) + c) & ((1 << 64) - 1) for r in reais], dtype=np.uint64)
        res[:, c, :] = normais(sem, n)
    return res


def le_desloc(arq, M=M_PADRAO):
    """mu (2, M) de um ustar.txt (2K numeros: os K de c' e os K de phi') ou de um form.json (chave 'u')"""
    if arq.endswith(".json"):
        import json
        u = np.array(json.load(open(arq))["u"], float)
    else:
        u = np.loadtxt(arq).ravel()
    K = len(u) // 2
    mu = np.zeros((2, M))
    mu[0, :K], mu[1, :K] = u[:K], u[K:]
    return mu


# ------------------------------------------------------------------- KL analitica (copia do klsep::KL2D do C++)
def _bissecao(g, lo, hi):
    glo = g(lo)
    for _ in range(200):
        m = 0.5 * (lo + hi)
        gm = g(m)
        if (gm > 0.) == (glo > 0.):
            lo, glo = m, gm
        else:
            hi = m
        if hi - lo < 1.e-15 * max(1., hi):
            break
    return 0.5 * (lo + hi)


def _modos1d(comp, L, n):
    a, c, eps = 0.5 * comp, 1.0 / L, 1e-12
    v, k = [], 1
    while len(v) < n:
        w = _bissecao(lambda x: x * math.tan(x * a) - c, ((k - 1) * math.pi + eps) / a, ((k - 0.5) * math.pi - eps) / a)
        v.append((2 * c / (w * w + c * c), w, 1))
        if len(v) >= n:
            break
        w = _bissecao(lambda x: x + c * math.tan(x * a), ((k - 0.5) * math.pi + eps) / a, (k * math.pi - eps) / a)
        v.append((2 * c / (w * w + c * c), w, 0))
        k += 1
    v.sort(key=lambda p: -p[0])
    om, par = np.array([p[1] for p in v]), np.array([p[2] for p in v])
    s2 = np.sin(2 * om * a) / (2 * om)
    return dict(lam=np.array([p[0] for p in v]), om=om, par=par, norma=1 / np.sqrt(np.where(par == 1, a + s2, a - s2)))


class KLAnalitica:
    """o retangulo [0, 30] x [0, 15] do Cho2010 (coordenadas do Cho), l_h = 20, l_v = 2, M termos"""

    def __init__(self, M=M_PADRAO, Dx=DX, Dy=DY, Lx=LH, Ly=LV):
        self.Dx, self.Dy, self.M = Dx, Dy, M
        mx, my = _modos1d(Dx, Lx, M), _modos1d(Dy, Ly, M)
        cand = sorted(((mx["lam"][i] * my["lam"][j], i, j) for i in range(M) for j in range(M)
                       if (i + 1) * (j + 1) <= M), key=lambda t: -t[0])[:M]
        ix, iy = np.array([t[1] for t in cand]), np.array([t[2] for t in cand])
        self.raiz = np.sqrt(np.array([t[0] for t in cand]))
        self.fracao = float((self.raiz ** 2).sum() / (Dx * Dy))
        self.ux, self.ia = np.unique(ix, return_inverse=True)      # modos 1D usados e a coluna de cada termo
        self.uy, self.ib = np.unique(iy, return_inverse=True)
        self.mx, self.my = mx, my

    @staticmethod
    def _f(m, idx, s):
        s = np.asarray(s, float)[:, None]
        om, par, nor = m["om"][idx][None, :], m["par"][idx][None, :], m["norma"][idx][None, :]
        return nor * np.where(par == 1, np.cos(om * s), np.sin(om * s))

    def Fx(self, X):
        return self._f(self.mx, self.ux, np.asarray(X) - 0.5 * self.Dx)

    def Gy(self, Y):
        return self._f(self.my, self.uy, np.asarray(Y) - 0.5 * self.Dy)

    def coef(self, xi):
        """A (nx_usados, ny_usados) de um vetor xi (M,)"""
        A = np.zeros((len(self.ux), len(self.uy)))
        np.add.at(A, (self.ia, self.ib), self.raiz * xi[:self.M])
        return A

    def psi(self, X, Y):
        """sqrt(lambda_k) f_k(X) g_k(Y) nos pontos: (npontos, M) (para testes e projecoes)"""
        return self.raiz[None, :] * self.Fx(X)[:, self.ia] * self.Gy(Y)[:, self.ib]


# ------------------------------------------------------------------------------------------------- base de Nystrom
def le_nystrom(arq):
    """le a base de Nystrom (texto do nystrom_base.py) e guarda um .npz ao lado (a leitura do texto leva ~1 min)"""
    npz = os.path.splitext(arq)[0] + ".npz"
    if os.path.exists(npz) and os.path.getmtime(npz) >= os.path.getmtime(arq):
        d = np.load(npz)
        return {k: d[k] for k in d.files}
    with open(arq) as f:
        n, M = map(int, f.readline().split())
        nos = np.array([list(map(float, f.readline().split())) for _ in range(n)])
        lam = np.array(list(map(float, f.readline().split())))
        F = np.loadtxt(f, max_rows=n).reshape(n, M)
    d = dict(x=nos[:, 0], y=nos[:, 1], w=nos[:, 2], lam=lam, F=F)
    try:
        tmp = f"{npz}.{os.uname().nodename}.{os.getpid()}.tmp.npz"
        np.savez(tmp, **d)
        os.replace(tmp, npz)
    except OSError:
        pass
    return d


def psi_nystrom(base, X, Y, M=M_PADRAO, bloco=2000, dtype=np.float32):
    """Psi = C(pontos, nos) W F Lambda^(-1/2): (npontos, M). Z(pontos) = Psi xi."""
    B = (base["F"][:, :M] / np.sqrt(base["lam"][:M])[None, :]) * base["w"][:, None]
    X, Y = np.asarray(X, float).ravel(), np.asarray(Y, float).ravel()
    out = np.empty((len(X), M), dtype)
    for a in range(0, len(X), bloco):
        b = min(a + bloco, len(X))
        C = np.exp(-np.abs(X[a:b, None] - base["x"][None, :]) / LH - np.abs(Y[a:b, None] - base["y"][None, :]) / LV)
        out[a:b] = C @ B
    return out


# ------------------------------------------------------------------------------------------------- campos na grade
class CamposCho:
    """c' (kPa) e phi' (rad) na grade da analise limite (coordenadas do Ceron) dados os xi (2, M).

    campo = 'analitico' ou 'nystrom' (base = o arquivo da base de Nystrom; cache = .npy de Psi, refeito se faltar)."""

    def __init__(self, r, campo="analitico", h=0.025, M=M_PADRAO, base=None, cache=None):
        self.r, self.campo, self.h, self.M = float(r), campo, float(h), M
        nx, ny = int(round(DX / h)) + 1, int(round(DY / h)) + 1
        self.x = -X_CRISTA + h * np.arange(nx)            # coordenadas do Ceron (x para fora, y para baixo)
        self.y = h * np.arange(ny)
        self.X, self.Y = self.x + X_CRISTA, Y_CRISTA - self.y   # coordenadas do Cho2010
        self.shape = (nx, ny)
        if campo == "analitico":
            self.kl = KLAnalitica(M)
            self.Fx, self.Gy = self.kl.Fx(self.X), self.kl.Gy(self.Y)
            self.fracao = self.kl.fracao
        elif campo == "nystrom":
            # so os nos da grade no solo ou a menos de uma celula da superficie (os outros nunca sao usados)
            xs, ys = np.meshgrid(self.x, self.y, indexing="ij")
            ysup = np.where(xs <= 0, 0.0, np.where(xs < 10.0, xs, 10.0))
            self.mask = (ys >= ysup - 1.5 * h).ravel()
            if cache and os.path.exists(cache):
                self.psi = np.load(cache, mmap_mode="r")
                if self.psi.shape != (int(self.mask.sum()), M):
                    raise ValueError(f"{cache}: forma {self.psi.shape}, esperada {(int(self.mask.sum()), M)}")
                fr = os.path.splitext(cache)[0] + "_fracao.txt"
                self.fracao = float(open(fr).read()) if os.path.exists(fr) else float("nan")
            else:
                if base is None:
                    raise ValueError("campo nystrom precisa de base= (ou do cache .npy)")
                b = le_nystrom(base) if isinstance(base, str) else base
                Xg, Yg = (xs + X_CRISTA).ravel()[self.mask], (Y_CRISTA - ys).ravel()[self.mask]
                self.psi = psi_nystrom(b, Xg, Yg, M)
                self.fracao = float(b["lam"][:M].sum() / b["w"].sum())
                if cache:
                    os.makedirs(os.path.dirname(os.path.abspath(cache)), exist_ok=True)
                    # nome temporario com no + pid (dois jobs em nos diferentes podem calcular ao mesmo tempo); a
                    # fracao e gravada antes do .npy, entao quem acha o .npy acha tambem a fracao
                    tmp = cache + f".{os.uname().nodename}.{os.getpid()}.tmp.npy"
                    fr = os.path.splitext(cache)[0] + "_fracao.txt"
                    open(fr + f".{os.uname().nodename}.{os.getpid()}", "w").write(repr(self.fracao))
                    os.replace(fr + f".{os.uname().nodename}.{os.getpid()}", fr)
                    np.save(tmp, self.psi)
                    os.replace(tmp, cache)
                    self.psi = np.load(cache, mmap_mode="r")
        else:
            raise ValueError("campo = analitico ou nystrom")

    def Z(self, xi):
        """os dois campos gaussianos (sem a correlacao cruzada) na grade: (2, nx, ny)"""
        xi = np.asarray(xi, float).reshape(2, -1)
        if self.campo == "analitico":
            return np.stack([self.Fx @ self.kl.coef(xi[f]) @ self.Gy.T for f in range(2)])
        z = np.zeros((2, self.mask.size))
        z[:, self.mask] = (np.asarray(self.psi) @ xi[:, :self.M].T.astype(self.psi.dtype)).T
        return z.reshape(2, *self.shape)

    def campos(self, xi):
        """(c em kPa, phi em rad) na grade, (nx, ny) cada"""
        z = self.Z(xi)
        zc, zf = z[0], self.r * z[0] + math.sqrt(1.0 - self.r ** 2) * z[1]
        return lognormal(MEDIA_C, COV_C, zc), np.radians(lognormal(MEDIA_PHI, COV_PHI, zf))


def campos_em_pontos(xi, X, Y, r, campo="analitico", base=None, M=M_PADRAO):
    """c (kPa) e phi (graus) em pontos soltos (coordenadas do Cho2010): o teste contra o campos.csv do C++"""
    xi = np.asarray(xi, float).reshape(2, -1)
    if campo == "analitico":
        P = KLAnalitica(M).psi(X, Y)
    else:
        P = psi_nystrom(le_nystrom(base) if isinstance(base, str) else base, X, Y, M, dtype=np.float64)
    zc, zf = P @ xi[0, :M], P @ xi[1, :M]
    zf = r * zc + math.sqrt(1 - r * r) * zf
    return lognormal(MEDIA_C, COV_C, zc), lognormal(MEDIA_PHI, COV_PHI, zf)


_PC = {}


def _psi_bloco(ab):
    a, b = ab
    out = np.lib.format.open_memmap(_PC["tmp"], mode="r+")
    out[a:b] = psi_nystrom(_PC["base"], _PC["X"][a:b], _PC["Y"][a:b], _PC["M"])
    out.flush()
    return b - a


def construir_cache(base, cache, h=0.025, M=M_PADRAO, procs=1, bloco=20000):
    """Psi de Nystrom na grade direto num .npy (memmap), em procs processos (fork): o mesmo do CamposCho(cache=...),
    para quando o BLAS do numpy e de uma thread so (o 'blas' de referencia do sistema)"""
    import multiprocessing as mp
    nx, ny = int(round(DX / h)) + 1, int(round(DY / h)) + 1
    x, y = -X_CRISTA + h * np.arange(nx), h * np.arange(ny)
    xs, ys = np.meshgrid(x, y, indexing="ij")
    ysup = np.where(xs <= 0, 0.0, np.where(xs < 10.0, xs, 10.0))
    mask = (ys >= ysup - 1.5 * h).ravel()                 # o mesmo criterio do CamposCho
    b = le_nystrom(base)
    no = f"{os.uname().nodename}.{os.getpid()}"
    tmp = cache + f".{no}.tmp.npy"
    os.makedirs(os.path.dirname(os.path.abspath(cache)), exist_ok=True)
    n = int(mask.sum())
    np.lib.format.open_memmap(tmp, mode="w+", dtype=np.float32, shape=(n, M)).flush()
    _PC.update(base=b, X=(xs + X_CRISTA).ravel()[mask], Y=(Y_CRISTA - ys).ravel()[mask], M=M, tmp=tmp)
    with mp.get_context("fork").Pool(procs) as pool:
        feitos = 0
        for k in pool.imap_unordered(_psi_bloco, [(a, min(a + bloco, n)) for a in range(0, n, bloco)]):
            feitos += k
            print(f"  Psi: {feitos}/{n} nos", flush=True)
    fr = os.path.splitext(cache)[0] + "_fracao.txt"
    open(fr + "." + no, "w").write(repr(float(b["lam"][:M].sum() / b["w"].sum())))
    os.replace(fr + "." + no, fr)
    os.replace(tmp, cache)


if __name__ == "__main__":
    # python3 -I campos_cho.py cache base=<base_nys2000.txt> cache=<psi.npy> [h=0.025] [procs=1]: so pre-calcula Psi de
    # Nystrom (procs > 1: em processos, para BLAS de uma thread; com OpenBLAS de varias threads procs=1 basta)
    import sys
    import time
    kv = dict(x.split("=", 1) for x in sys.argv[2:] if "=" in x)
    if len(sys.argv) < 2 or sys.argv[1] != "cache" or "base" not in kv or "cache" not in kv:
        sys.exit("uso: python3 -I campos_cho.py cache base=<base_nys2000.txt> cache=<psi.npy> [h=0.025]")
    t0 = time.perf_counter()
    if int(kv.get("procs", 1)) > 1 and not os.path.exists(kv["cache"]):
        construir_cache(kv["base"], kv["cache"], float(kv.get("h", 0.025)), procs=int(kv["procs"]))
    C = CamposCho(0.0, "nystrom", float(kv.get("h", 0.025)), base=kv["base"], cache=kv["cache"])
    print(f"Psi {C.psi.shape} em {kv['cache']} (variancia retida {C.fracao:.4f}), {time.perf_counter() - t0:.0f} s")
