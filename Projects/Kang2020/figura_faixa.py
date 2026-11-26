#!/usr/bin/env python3
"""figura_faixa.py <malha.csv> <saida sem extensao>: a malha da faixa do FLAC (Passo 6) com os contornos da Fig. 5."""
import csv
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.collections import LineCollection, PolyCollection  # noqa: E402

TINTA, TINTA2, SUPERF, AZUL, SOLO = "#0b0b0b", "#5c5c58", "#fcfcfb", "#2a78d6", "#efe6d8"
arq, saida = sys.argv[1], sys.argv[2]
els = [r for r in csv.DictReader(l for l in open(arq) if not l.startswith("#"))]
tri = [[(float(r[f"x{i}"]), float(r[f"y{i}"])) for i in range(3)] for r in els if r["id"] == "1"]
lin = {k: [[(float(r["x0"]), float(r["y0"])), (float(r["x1"]), float(r["y1"]))] for r in els if r["id"] == k]
       for k in ("-1", "-2", "-5", "-6")}
plt.rcParams.update({"font.size": 8, "savefig.facecolor": SUPERF, "figure.facecolor": SUPERF})
fig, ax = plt.subplots(figsize=(7.4, 4.4))
ax.add_collection(PolyCollection(tri, facecolors=SOLO, edgecolors=TINTA2, linewidths=0.15))
ax.add_collection(LineCollection(lin["-6"], colors=AZUL, linewidths=2.0))
ax.add_collection(LineCollection(lin["-1"], colors=TINTA, linewidths=2.0))
ax.add_collection(LineCollection(lin["-2"] + lin["-5"], colors=TINTA, linewidths=1.2))
xs = [x for t in tri for x, _ in t]
ys = [y for t in tri for _, y in t]
ax.plot([5, 5], [5 * (max(ys) - 2) / 10, 5 * (max(ys) - 2) / 10 + 2], color="#e34948", lw=1.4)
ax.annotate("perfil do FS\n(x = 5 m)", (5, 5 * (max(ys) - 2) / 10 + 1), xytext=(5.6, 1.2), fontsize=7,
            color="#e34948", arrowprops=dict(arrowstyle="-", color="#e34948", lw=0.7))
ax.text(0.3, max(ys) - 0.4, "topo: $q_w$ = I(t) normal à superfície, ou $p_w$ = 0 sob lâmina;\n$p_a$ = 0 no modelo TP",
        fontsize=7.5, color=TINTA, va="top")
ax.text(9.9, 0.3, "laterais: $u_x$ = 0\nbase: $u_x$ = $u_y$ = 0\nsem fluxo nas laterais e na base", fontsize=7.5,
        color=TINTA2, ha="right")
ax.set_aspect("equal")
ax.set_xlim(min(xs) - 0.2, max(xs) + 0.2)
ax.set_ylim(min(ys) - 0.3, max(ys) + 0.3)
ax.set_xlabel("x (m)")
ax.set_ylabel("y (m)")
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
ax.set_title(f"A malha do FLAC (Kang et al. 2020): 40 x 20 zonas, 4 triângulos por zona ({len(tri)})",
             fontsize=8.5, loc="left", color=TINTA)
fig.tight_layout()
fig.savefig(saida + ".pdf")
fig.savefig(saida + ".png", dpi=150)
print("pronto:", saida + ".pdf")
