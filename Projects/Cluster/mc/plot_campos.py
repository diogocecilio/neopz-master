#!/usr/bin/env python3
"""plot_campos.py <campos.csv> [saida.png]: os campos c', phi' e k de uma realizacao (rodada com CAMPOS_CSV=1), nos pontos
de integracao, como na Figura 12 de Vargas Ceron et al. (2025). Requer numpy e matplotlib."""
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np

arq = sys.argv[1] if len(sys.argv) > 1 else "campos.csv"
saida = sys.argv[2] if len(sys.argv) > 2 else arq.replace(".csv", ".png")
d = np.genfromtxt(arq, delimiter=",", names=True)
x, y = d["x"], d["y"]
tri = mtri.Triangulation(x, y)
# so os triangulos dentro do talude: a triangulacao de Delaunay dos pontos cobre o fecho convexo, que inclui o vazio
# acima da face do talude; descarta os triangulos com lado maior que 3 m (os elementos da malha tem ate 3 m)
xt, yt = x[tri.triangles], y[tri.triangles]
lados = np.max(np.hypot(xt - np.roll(xt, 1, axis=1), yt - np.roll(yt, 1, axis=1)), axis=1)
tri.set_mask(lados > 3.0)
fig, ax = plt.subplots(3, 1, figsize=(9, 10), constrained_layout=True)
for a, (campo, titulo, cmap, log) in zip(ax, [("c", "coesao c' (kPa)", "viridis", False),
                                            ("phi_deg", "angulo de atrito phi' (graus)", "plasma", False),
                                            ("k_mult", "k / k_s (multiplicador da permeabilidade)", "cividis", True)]):
    v = np.log10(d[campo]) if log else d[campo]
    cf = a.tricontourf(tri, v, levels=30, cmap=cmap)
    cb = fig.colorbar(cf, ax=a)
    if log:
        cb.set_label("log10")
    a.set_title(f"{titulo}: media {np.mean(d[campo]):.3g}, CoV {np.std(d[campo]) / np.mean(d[campo]):.2f}")
    a.set_aspect("equal")
    a.set_xlabel("x (m)")
    a.set_ylabel("y (m)")
fig.savefig(saida, dpi=110)
print(saida)
