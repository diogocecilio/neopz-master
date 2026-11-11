#!/usr/bin/env python3
"""Forcas de percolacao -grad u nos pontos de integracao: MEF (rebaixforca) e semi-analitico, lado a lado.

Uso: python3 figuras_forca.py [pasta] [material: up | unsat]   (padrao: pasta atual, up)
Le gim_rebaixforca_ceron_<mat>_nivel0_forca.csv e gim_semianalitico_ceron_<mat>_nivel0_forca.csv (x, y, fx, fy) e a
banda do ultimo nivel de cada um (_nivel<k>_pontos.csv); grava forca_percolacao_<mat>.png.
"""
import csv
import glob
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

pasta = sys.argv[1] if len(sys.argv) > 1 else "."
mat = sys.argv[2] if len(sys.argv) > 2 else "up"
casos = [("rebaixforca", "MEF: $-$grad u da percolacao do MEF"), ("semianalitico", "campo semi-analitico (artigo, Sec. 4.1)")]

fig, eixos = plt.subplots(len(casos), 1, figsize=(10, 4.4 * len(casos)))
fmax = 0.
dados = []
for nome, _ in casos:
    p = np.loadtxt(os.path.join(pasta, f"gim_{nome}_ceron_{mat}_nivel0_forca.csv"), delimiter=",", skiprows=1)
    niveis = sorted(glob.glob(os.path.join(pasta, f"gim_{nome}_ceron_{mat}_nivel*_pontos.csv")))
    b = np.loadtxt(niveis[-1], delimiter=",", skiprows=1)
    b = b[b[:, 2] >= 0.1 * b[:, 2].max()]
    with open(os.path.join(pasta, f"gim_{nome}_ceron_{mat}.csv")) as f:
        gamma = [float(r["gamma"]) for r in csv.DictReader(f)][-1]
    dados.append((p, b, gamma))
    fmax = max(fmax, np.percentile(np.hypot(p[:, 2], p[:, 3]), 99))

for ax, (nome, titulo), (p, b, gamma) in zip(eixos, casos, dados):
    mag = np.hypot(p[:, 2], p[:, 3])
    sel = np.arange(0, len(p), 2)  # metade dos pontos, para as setas nao se sobreporem
    q = ax.quiver(p[sel, 0], p[sel, 1], p[sel, 2] / np.maximum(mag[sel], 1e-12), p[sel, 3] / np.maximum(mag[sel], 1e-12),
                  np.minimum(mag[sel], fmax), cmap="viridis", scale=60, width=0.0022, clim=(0, fmax))
    ax.plot([0, 10, 15, 25, 25, 0, 0], [10, 10, 5, 5, 0, 0, 10], "k-", lw=1)
    ax.scatter(b[:, 0], b[:, 1], s=2, c="tab:red", label="banda do ultimo nivel")
    ax.set_aspect("equal")
    ax.set_xlim(-0.5, 25.5)
    ax.set_ylim(-0.5, 10.5)
    ax.set_title(f"{titulo}: $\\Gamma$ = " + f"{gamma:.3f}".replace(".", ","), fontsize=10)
    ax.set_ylabel("y (m)")
    ax.legend(loc="lower left", fontsize=8)
    fig.colorbar(q, ax=ax, label="|$-$grad u| (kN/m$^3$)")
eixos[-1].set_xlabel("x (m)")
fig.suptitle("Forcas de percolacao apos o rebaixamento total (direcao: setas; intensidade: cor)", fontsize=11)
fig.tight_layout()
saida = os.path.join(pasta, f"forca_percolacao_{mat}.png")
fig.savefig(saida, dpi=130)
print("gravado", saida)
