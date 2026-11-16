#!/usr/bin/env python3
"""nystrom_base.py: a base de Karhunen-Loeve de Nystrom do talude do Cho (2010) para o RandomField (opcao nystrom= do
Cho2010).

O kernel e o do RandomField: C = exp(-|dx|/lh - |dy|/lv) (exponencial separavel). A quadratura e a regra do ponto medio
numa grade de celulas dx x dy: os nos sao os centros das celulas que caem dentro do talude, e o peso de cada no e a area
da celula. Os autopares saem de W^(1/2) C W^(1/2) v = lambda v, com f_k = v_k / sqrt(w) nos nos (sum_j w_j f_k^2 = 1). O
RandomField interpola f_k em qualquer ponto pela formula de Nystrom: f_k(x) = (1/lambda_k) sum_j w_j C(x, x_j) f_k(x_j).

Uso: python3 -I nystrom_base.py <saida.txt> [M=200] [lh=20] [lv=2] [dx=0.5] [dy=0.25]
Arquivo: 'n M' e, depois, n linhas 'x y w', 1 linha com os M autovalores e n linhas com os M valores f_k(x_j)."""
import sys

import numpy as np
from scipy.linalg import eigh

TALUDE = [(0, 0), (30, 0), (30, 5), (20, 5), (10, 15), (0, 15)]  # o cho2010.cpp (kTalude)


def dentro(x, y):
    """ponto dentro do poligono do talude (raio)"""
    c = False
    n = len(TALUDE)
    for i in range(n):
        x1, y1 = TALUDE[i]
        x2, y2 = TALUDE[(i + 1) % n]
        if (y1 > y) != (y2 > y) and x < x1 + (y - y1) * (x2 - x1) / (y2 - y1):
            c = not c
    return c


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    p = {"M": 200, "lh": 20., "lv": 2., "dx": 0.5, "dy": 0.25}
    for a in sys.argv[2:]:
        k, v = a.split("=", 1)
        p[k] = type(p[k])(v)
    xs = np.arange(p["dx"] / 2, 30, p["dx"])
    ys = np.arange(p["dy"] / 2, 15, p["dy"])
    pts = np.array([(x, y) for y in ys for x in xs if dentro(x, y)])
    n = len(pts)
    w = np.full(n, p["dx"] * p["dy"])
    C = np.exp(-np.abs(pts[:, 0, None] - pts[None, :, 0]) / p["lh"] - np.abs(pts[:, 1, None] - pts[None, :, 1]) / p["lv"])
    s = np.sqrt(w)
    K = C * s[:, None] * s[None, :]
    M = min(p["M"], n)
    lam, v = eigh(K, subset_by_index=[n - M, n - 1])
    ordem = np.argsort(lam)[::-1]
    lam, v = lam[ordem], v[:, ordem]
    f = v / s[:, None]
    area = w.sum()
    print(f"{n} nos (area {area:.1f} m2, a do talude e 300), M = {M}")
    for m in (10, 20, 50, 100, 200):
        if m <= M:
            print(f"  variancia retida com M = {m:3d}: {lam[:m].sum() / area:.3f}")
    with open(sys.argv[1], "w") as arq:
        arq.write(f"{n} {M}\n")
        for (x, y), wi in zip(pts, w):
            arq.write(f"{x:.6f} {y:.6f} {wi:.8f}\n")
        arq.write(" ".join(f"{x:.12e}" for x in lam) + "\n")
        for j in range(n):
            arq.write(" ".join(f"{x:.10e}" for x in f[j]) + "\n")


if __name__ == "__main__":
    main()
