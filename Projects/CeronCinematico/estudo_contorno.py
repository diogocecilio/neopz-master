#!/usr/bin/env python3
"""estudo_contorno.py: o Gamma deterministico do Ceron (2025) com percolacao PERMANENTE saturada (Eqs. 19-21),
variando o dominio e as condicoes de contorno da agua, para achar o que leva ao 1,336 do artigo (Fig. 24).

Malha propria (estruturada mapeada, P1), c = 10, phi = 30, gamma = 20, H = 5, beta = 45, hw = H, k isotropico.
O Gamma e o do metodo cinematico (cinematico.minimizar), que no passo 1 coincidiu com o GIM extrapolado (+0,3%).

    python3 -I estudo_contorno.py            -> tabela no terminal e dados/estudo_contorno.csv
"""
import csv
import math
import os
import sys
import time

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cinematico as k  # noqa: E402

GW = k.GAMMA_W
H, BETA = 5.0, math.radians(45.0)


def graduada(a, b, h_fino, h_grosso, foco):
    """Pontos de a a b, espacamento h_fino perto de foco (lista de pontos) e crescendo ate h_grosso."""
    xs = [a]
    while xs[-1] < b - 1e-9:
        d = min(abs(xs[-1] - f) for f in foco)
        h = min(h_grosso, h_fino + 0.25 * d)
        xs.append(min(b, xs[-1] + h))
    return np.array(xs)


def malha(Le, Ld, D, h=0.125):
    """Dominio: crista de -Le a 0, face de 0 a xT, terreno do pe de xT a xT + Ld, base em y = H + D (y para baixo)."""
    tal = k.Talude(H=H, beta=BETA, x_esq=-Le, x_dir=H / math.tan(BETA) + Ld, y_base=H + D)
    xT = tal.xT
    xs = np.unique(np.r_[graduada(-Le, xT + Ld, h, 2.0, [0.0, xT]), 0.0, xT])
    eta = graduada(0.0, 1.0, h / (H + D), 2.0 / (H + D), [0.0])   # fino perto da superficie
    nx, ny = len(xs), len(eta)
    ys = tal.y_sup(xs)
    X = np.repeat(xs[:, None], ny, 1)
    Y = ys[:, None] + (tal.y_base - ys[:, None]) * eta[None, :]
    xy = np.c_[X.ravel(), Y.ravel()]
    idx = np.arange(nx * ny).reshape(nx, ny)
    a, b, c, d = idx[:-1, :-1].ravel(), idx[1:, :-1].ravel(), idx[1:, 1:].ravel(), idx[:-1, 1:].ravel()
    tri = np.r_[np.c_[a, b, c], np.c_[a, c, d]]
    return tal, xy, tri


def permanente(xy, tri, tal, esq, dir_, base, hw=H):
    """div(grad u) = 0; Dirichlet na crista (0), na face (-gw y ate hw, -gw hw abaixo) e no terreno do pe (-gw hw).
    esq / dir_ / base: None (sem fluxo) ou o valor de u ('0' ou 'hw' = -gw hw)."""
    n = len(xy)
    P = xy[tri]
    B = np.stack([P[:, 1] - P[:, 0], P[:, 2] - P[:, 0]], axis=2)
    A = 0.5 * np.abs(np.linalg.det(B))
    G = np.einsum("kji,jl->kil", np.linalg.inv(B), np.array([[-1.0, 1, 0], [-1.0, 0, 1]]))
    Ke = np.einsum("k,kia,kib->kab", A, G, G)
    K = sp.csr_matrix((Ke.ravel(), (np.repeat(tri, 3, 1).ravel(), np.tile(tri, (1, 3)).ravel())), shape=(n, n))
    x, y = xy[:, 0], xy[:, 1]
    t = 1e-7
    ud = np.full(n, np.nan)
    val = {"0": 0.0, "hw": -GW * hw}
    if base is not None:
        ud[np.abs(y - tal.y_base) < t] = val[base]
    if esq is not None:
        ud[np.abs(x - tal.x_esq) < t] = val[esq]
    if dir_ is not None:
        ud[np.abs(x - tal.x_dir) < t] = val[dir_]
    sup = np.abs(y - tal.y_sup(x)) < t
    ud[sup & (x <= t)] = 0.0
    face = sup & (x > t) & (x < tal.xT - t)
    ud[face] = np.where(y[face] <= hw, -GW * y[face], -GW * hw)
    ud[sup & (x >= tal.xT - t)] = -GW * hw
    D = ~np.isnan(ud)
    u = np.where(D, ud, 0.0)
    F = ~D
    u[F] = spl.spsolve(K[F][:, F].tocsc(), -K[F][:, D] @ u[D])
    return u


DOMINIOS = {  # (Le, Ld, D) em m
    "notebook (3,5H; 3H; 3H)": (17.5, 15.0, 15.0),
    "nosso transiente (7H; 6H; 6H)": (35.0, 30.0, 30.0),
    "artigo Fig.4 aprox. (13H; 2H; 8H)": (65.0, 10.0, 40.0),
}
CONTORNOS = {  # (esq, dir, base)
    "A: so crista/face/pe (resto sem fluxo; notebook cel. 1)": (None, None, None),
    "B: + esq u=0, dir u=-gw hw (notebook cel. 9)": ("0", "hw", None),
    "C: tudo Dirichlet (dOmega_v vazio): esq 0, base 0, dir -gw hw": ("0", "hw", "0"),
    "D: esq 0, base 0, dir sem fluxo": ("0", None, "0"),
}


def main():
    saida = os.path.join(os.path.dirname(os.path.abspath(__file__)), "dados", "estudo_contorno.csv")
    linhas = []
    print(f"{'dominio':36s} {'contorno':62s} {'nos':>6s} {'Gamma':>7s} {'Gamma/1,336':>11s}  mec", flush=True)
    for nd, (Le, Ld, D) in DOMINIOS.items():
        tal, xy, tri = malha(Le, Ld, D)
        n = len(xy)
        for nc, (e, d, b) in CONTORNOS.items():
            t0 = time.perf_counter()
            u = permanente(xy, tri, tal, e, d, b)
            campo = k.CampoMalha(xy, tri, np.full(n, 10.0), np.full(n, math.radians(30.0)), u, 20.0 - GW)
            r = k.minimizar(campo, tal)
            m = r.mecanismo
            print(f"{nd:36s} {nc:62s} {n:6d} {r.Gamma:7.4f} {r.Gamma / 1.336:11.3f}  {m.classe} xA={m.xA:.2f} "
                  f"({time.perf_counter() - t0:.0f} s)", flush=True)
            linhas.append({"dominio": nd, "contorno": nc, "nos": n, "Gamma": r.Gamma, "classe": m.classe, "xA": m.xA})
    with open(saida, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(linhas[0]))
        w.writeheader()
        w.writerows(linhas)
    print("salvo em", saida)


if __name__ == "__main__":
    main()
