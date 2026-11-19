#!/usr/bin/env python3
"""nystrom_base.py: a base de Karhunen-Loeve de Nystrom do talude do Cho (2010) para o RandomField (opcao nystrom= do
Cho2010).

O kernel e o do RandomField: C = exp(-|dx|/lh - |dy|/lv) (exponencial separavel). A quadratura e a regra do ponto medio
numa grade de celulas dx x dy: os nos sao os centros das celulas que caem dentro do talude, e o peso de cada no e a area
da celula. Os autopares saem de W^(1/2) C W^(1/2) v = lambda v, com f_k = v_k / sqrt(w) nos nos (sum_j w_j f_k^2 = 1). O
RandomField interpola f_k em qualquer ponto pela formula de Nystrom: f_k(x) = (1/lambda_k) sum_j w_j C(x, x_j) f_k(x_j).

Uso: python3 -I nystrom_base.py <saida.txt> [M=200] [lh=20] [lv=2] [dx=0.5] [dy=0.25] [copias=] [nos=] [autopares=]
  copias=20,100    tambem grava <saida>_M20.txt e <saida>_M100.txt com os 20 e os 100 primeiros modos (os mesmos
                   autovetores, com os mesmos sinais: os M primeiros termos de uma valem para a outra). O FORM com poucos
                   termos le o arquivo pequeno; a amostragem com todos os termos, o grande
  nos=<arquivo>    so grava os nos (x y por linha) e sai: para calcular os autopares fora (base grande, ver abaixo)
  autopares=<pref> le os autopares de <pref>.lam (M autovalores) e <pref>.f (n x M valores f_k(x_j), linha por no), em
                   Real64, no lugar do eigh

Base grande (milhares de nos e de modos): o eigh daqui pode ser lento (LAPACK sem threads). O caminho e:
  python3 -I nystrom_base.py base.txt M=2000 dx=0.2 dy=0.1 nos=base.nos
  wolframscript -file nystrom_autopares.wl base.nos 2000 20 2 0.02 base      # o Eigensystem do Mathematica (MKL)
  python3 -I nystrom_base.py base.txt M=2000 dx=0.2 dy=0.1 autopares=base copias=100
A grade de 0,2 x 0,1 m (15 000 nos) retem 99,7% da variancia com 2000 modos; a de 0,5 x 0,25 m (2400 nos) erra os
autovalores acima de ~50 modos (lambda_100 23% e lambda_1000 57% contra a de 0,2 x 0,1).
Arquivo: 'n M' e, depois, n linhas 'x y w', 1 linha com os M autovalores e n linhas com os M valores f_k(x_j)."""
import sys

import numpy as np

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


def grava(arq, pts, w, lam, f):
    """a base no formato que o RandomField::UseNystrom le"""
    n, M = f.shape
    with open(arq, "w") as a:
        a.write(f"{n} {M}\n")
        np.savetxt(a, np.column_stack([pts, w]), fmt=["%.6f", "%.6f", "%.8f"])
        a.write(" ".join(f"{x:.12e}" for x in lam) + "\n")
        np.savetxt(a, f, fmt="%.10e")


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    p = {"M": 200, "lh": 20., "lv": 2., "dx": 0.5, "dy": 0.25, "copias": "", "nos": "", "autopares": ""}
    for a in sys.argv[2:]:
        k, v = a.split("=", 1)
        if k not in p:
            sys.exit(f"opcao desconhecida: {a}")
        p[k] = type(p[k])(v)
    xs = np.arange(p["dx"] / 2, 30, p["dx"])
    ys = np.arange(p["dy"] / 2, 15, p["dy"])
    pts = np.array([(x, y) for y in ys for x in xs if dentro(x, y)])
    n = len(pts)
    w = np.full(n, p["dx"] * p["dy"])
    if p["nos"]:
        np.savetxt(p["nos"], pts, fmt="%.6f")
        print(f"{n} nos em {p['nos']}")
        return
    M = min(p["M"], n)
    if p["autopares"]:
        lam = np.fromfile(p["autopares"] + ".lam", dtype="<f8")
        f = np.fromfile(p["autopares"] + ".f", dtype="<f8")
        if len(lam) < M or f.size != n * len(lam):
            sys.exit(f"autopares de {p['autopares']}: {len(lam)} autovalores e {f.size} valores, esperados {M} e n x M"
                     f" com n = {n}")
        f = f.reshape(n, len(lam))[:, :M]
        lam = lam[:M]
        if np.any(np.diff(lam) > 0):
            sys.exit("os autovalores nao estao em ordem decrescente")
    else:
        from scipy.linalg import eigh
        C = np.exp(-np.abs(pts[:, 0, None] - pts[None, :, 0]) / p["lh"] - np.abs(pts[:, 1, None] - pts[None, :, 1]) / p["lv"])
        s = np.sqrt(w)
        K = C * s[:, None] * s[None, :]
        lam, v = eigh(K, subset_by_index=[n - M, n - 1])
        ordem = np.argsort(lam)[::-1]
        lam, v = lam[ordem], v[:, ordem]
        f = v / s[:, None]
    area = w.sum()
    print(f"{n} nos (area {area:.1f} m2, a do talude e 300), M = {M}")
    print(f"  normalizacao: max |sum_j w_j f_k^2 - 1| = {np.abs((w[:, None] * f * f).sum(axis=0) - 1).max():.2e}")
    for m in (10, 20, 50, 100, 200, 500, 1000, 2000):
        if m <= M:
            print(f"  variancia retida com M = {m:4d}: {lam[:m].sum() / area:.4f}")
    grava(sys.argv[1], pts, w, lam, f)
    base = sys.argv[1][:-4] if sys.argv[1].endswith(".txt") else sys.argv[1]
    for m in (int(c) for c in p["copias"].split(",") if c):
        if m < M:
            grava(f"{base}_M{m}.txt", pts, w, lam[:m], f[:, :m])
            print(f"  copia com os {m} primeiros modos: {base}_M{m}.txt")


if __name__ == "__main__":
    main()
