#!/usr/bin/env python3
# plot_campos.py <saida.png> <titulo> <campos.csv> [<campos.csv> ...]: desenha c (e phi, se variar) nos pontos de
# integracao de cada campos.csv (SO_CAMPOS=1 ou CAMPOS_CSV=1 do V2_TaludesCho), um painel por arquivo e variavel, com a
# escala de cores da Fig. 4 de Liao et al. (2025) para c_u (14 a 28 kPa) e da Fig. 9 para c (6 a 14) e phi (20 a 40).
import csv
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri


def le(arq):
    x, y, c, p = [], [], [], []
    with open(arq, newline="") as f:
        for l in csv.DictReader(f):
            x.append(float(l["x"])), y.append(float(l["y"])), c.append(float(l["c"])), p.append(float(l["phi_deg"]))
    return x, y, c, p


saida, titulo, arqs = sys.argv[1], sys.argv[2], sys.argv[3:]
paineis = []
for a in arqs:
    x, y, c, p = le(a)
    nome = a.split("/")[-2] if "/" in a else a
    paineis.append((nome, "c (kPa)", x, y, c))
    if max(p) - min(p) > 1e-6:
        paineis.append((nome, "phi (graus)", x, y, p))
fig, axs = plt.subplots(len(paineis), 1, figsize=(7, 2.6 * len(paineis)), squeeze=False)
for ax, (nome, rot, x, y, v) in zip(axs[:, 0], paineis):
    t = mtri.Triangulation(x, y)
    # remove os triangulos de Delaunay por cima da face (fora do talude)
    zc = max(y)
    sup = lambda xx: zc if xx <= 10 else (5. if xx >= 20 else zc + (5. - zc) * (xx - 10.) / 10.)
    import numpy as np
    xm, ym = np.array(x)[t.triangles].mean(axis=1), np.array(y)[t.triangles].mean(axis=1)
    t.set_mask(ym > np.array([sup(v) for v in xm]) + 1e-6)
    if rot.startswith("c"):
        lim = (14., 28.) if zc <= 10.5 else (6., 14.)
    else:
        lim = (20., 40.)
    h = ax.tricontourf(t, v, levels=np.linspace(lim[0], lim[1], 15), cmap="Blues", extend="both")
    fig.colorbar(h, ax=ax, label=rot)
    ax.set_aspect("equal")
    ax.set_title(f"{nome}: {rot}, de {min(v):.1f} a {max(v):.1f}", fontsize=9)
    ax.set_xlabel("x (m)")
    ax.set_ylabel("z (m)")
fig.suptitle(titulo, fontsize=10)
fig.tight_layout()
fig.savefig(saida, dpi=110)
print(saida)
