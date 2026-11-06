#!/usr/bin/env python3
"""plot_campos.py <campos.csv> [saida.png]: os campos aleatorios de uma realizacao (o campos.csv que o programa grava com
CAMPOS_CSV=1, ou com SO_CAMPOS=1 nos que tem esse modo), nos pontos de integracao: os que houver entre c (c' ou c_u),
phi_deg e k_mult (k/k_s, em escala log), como na Figura 12 de Vargas Ceron et al. (2025). Requer numpy e matplotlib."""
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np

if len(sys.argv) < 2:
    print(__doc__)
    sys.exit(2)
arq = sys.argv[1]
saida = sys.argv[2] if len(sys.argv) > 2 else arq.replace(".csv", ".png")
d = np.genfromtxt(arq, delimiter=",", names=True)
x, y = d["x"], d["y"]
tri = mtri.Triangulation(x, y)
# so os triangulos dentro do dominio: a triangulacao de Delaunay dos pontos cobre o fecho convexo, que num talude inclui
# o vazio acima da face; descarta os triangulos com o lado maior que 3 vezes a mediana
xt, yt = x[tri.triangles], y[tri.triangles]
lados = np.max(np.hypot(xt - np.roll(xt, 1, axis=1), yt - np.roll(yt, 1, axis=1)), axis=1)
tri.set_mask(lados > 3.0 * np.median(lados))
campos = [(c, t, m, lg) for (c, t, m, lg) in [("c", "coesao (kPa)", "viridis", False),
                                              ("phi_deg", "angulo de atrito (graus)", "plasma", False),
                                              ("k_mult", "k / k_s (multiplicador da permeabilidade)", "cividis", True)]
          if c in d.dtype.names]
if not campos:
    sys.exit(f"{arq}: sem as colunas c, phi_deg ou k_mult")
fig, ax = plt.subplots(len(campos), 1, figsize=(9, 3.4 * len(campos)), constrained_layout=True, squeeze=False)
for a, (campo, titulo, cmap, log) in zip(ax[:, 0], campos):
    v = np.log10(d[campo]) if log else d[campo]
    if np.ptp(v) > 0:
        cf = a.tricontourf(tri, v, levels=30, cmap=cmap)
        cb = fig.colorbar(cf, ax=a)
        if log:
            cb.set_label("log10")
    else:
        a.scatter(x, y, s=2, color="0.6")
    m = np.mean(d[campo])
    a.set_title(f"{titulo}: media {m:.3g}, CoV {np.std(d[campo]) / m if m else 0:.2f}", fontsize=10)
    a.set_aspect("equal")
    a.set_xlabel("x (m)")
    a.set_ylabel("y (m)")
fig.savefig(saida, dpi=110)
print(saida)
