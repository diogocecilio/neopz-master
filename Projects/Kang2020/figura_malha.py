#!/usr/bin/env python3
"""figura_malha.py <malha.csv> <saida sem extensao>: a malha da coluna do Kang2020 (gravada pelo driver com malha=) em
tres paineis: (a) o talude infinito do artigo e a coluna vertical que o representa; (b) a coluna inteira em escala real;
(c) o topo ampliado, com os triangulos, os vertices (p_w, p_a e u), os nos do meio das arestas (so u, P2) e os
contornos."""
import csv
import math
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyArrowPatch, Polygon  # noqa: E402

TINTA, TINTA2, GRADE, SUPERF = "#0b0b0b", "#5c5c58", "#e4e3df", "#fcfcfb"
AZUL, LARANJA, SOLO = "#2a78d6", "#eb6834", "#efe6d8"

arq, saida = sys.argv[1], sys.argv[2]
cab = open(arq).readline()
els = [r for r in csv.DictReader(l for l in open(arq) if not l.startswith("#"))]
tri = [[(float(r[f"x{i}"]), float(r[f"y{i}"])) for i in range(3)] for r in els if r["id"] == "1"]
H = max(y for t in tri for _, y in t)
W = max(x for t in tri for x, _ in t)
ny = sum(1 for r in els if r["id"] == "-2")
neq = cab.split(",")[1].split()[0]

plt.rcParams.update({"font.size": 8, "savefig.facecolor": SUPERF, "figure.facecolor": SUPERF})
fig = plt.figure(figsize=(7.4, 5.0))
axA = fig.add_axes([0.01, 0.10, 0.43, 0.78])
axB = fig.add_axes([0.52, 0.12, 0.05, 0.74])
axC = fig.add_axes([0.64, 0.24, 0.21, 0.56])

# (a) o talude infinito: uma faixa inclinada de solo com 2 m na vertical sobre a base impermeavel; a coluna vertical
b = math.radians(30.)
L = 7.
sup = [(0., 0.), (L * math.cos(b), L * math.sin(b))]
poly = [sup[0], sup[1], (sup[1][0], sup[1][1] - H), (sup[0][0], sup[0][1] - H)]
axA.add_patch(Polygon(poly, closed=True, fc=SOLO, ec=TINTA2, lw=0.8))
axA.plot([poly[3][0] - 0.3, poly[2][0] + 0.3], [poly[3][1] - 0.3 * math.tan(b), poly[2][1] + 0.3 * math.tan(b)],
         color=TINTA, lw=1.6)  # a base impermeavel
for k in range(9):  # hachura da base
    x = poly[3][0] + k * (L * math.cos(b)) / 8
    y = poly[3][1] + (x - poly[3][0]) * math.tan(b)
    axA.plot([x, x - 0.25], [y, y - 0.25], color=TINTA2, lw=0.6)
xc = 0.5 * L * math.cos(b)
yc = xc * math.tan(b)
axA.add_patch(Polygon([(xc - 0.12, yc - 0.12 * math.tan(b)), (xc + 0.12, yc + 0.12 * math.tan(b)),
                       (xc + 0.12, yc + 0.12 * math.tan(b) - H), (xc - 0.12, yc - 0.12 * math.tan(b) - H)],
                      closed=True, fc=AZUL, alpha=0.35, ec=AZUL, lw=1.0))
for k in range(7):  # a chuva
    x = 0.6 + k * 0.85
    y = x * math.tan(b)
    axA.add_patch(FancyArrowPatch((x, y + 1.3), (x, y + 0.25), arrowstyle="-|>", mutation_scale=7, color=AZUL, lw=0.8))
axA.annotate("chuva horária\n(Fig. 2 do artigo)", (0.4, 3.6), fontsize=7, color=TINTA2)
axA.annotate("a coluna do modelo\n(2 m na vertical)", (xc + 0.05, yc - 1.4), xytext=(xc + 0.9, -2.4), fontsize=7,
             color=TINTA, arrowprops=dict(arrowstyle="-", color=TINTA2, lw=0.7))
xb = L * math.cos(b) - 1.6
axA.text(xb, xb * math.tan(b) - H - 0.55, "base impermeável", fontsize=7, color=TINTA2, rotation=30,
         rotation_mode="anchor")
axA.text(0.0, -H - 1.0, "β = 23–40° (só na Eq. 14)", fontsize=7, color=TINTA2)
axA.set_xlim(-0.3, L * math.cos(b) + 0.4)
axA.set_ylim(-H - 1.3, L * math.sin(b) + 1.6)
axA.set_aspect("equal")
axA.axis("off")
axA.set_title("(a) talude infinito de Kang et al. (2020)", fontsize=8.5, loc="left", color=TINTA)

# (b) a coluna inteira, em escala real, e o trecho ampliado em (c)
zoom = 4 * H / ny  # as 4 celulas do topo
for t in tri:
    axB.add_patch(Polygon(t, closed=True, fc=SOLO, ec=TINTA2, lw=0.25))
axB.add_patch(Polygon([(-0.004, H - zoom), (W + 0.004, H - zoom), (W + 0.004, H + 0.004), (-0.004, H + 0.004)],
                      closed=True, fill=False, ec=AZUL, lw=1.0))
axB.set_xlim(-0.01, W + 0.01)
axB.set_ylim(-0.02, H + 0.02)
axB.set_aspect("equal")
axB.set_xticks([])
axB.set_xlabel(f"{W * 100:g} cm", fontsize=6.5, color=TINTA2)
axB.set_yticks([0, 0.5, 1, 1.5, 2])
axB.set_yticklabels(["2", "1,5", "1", "0,5", "0"], fontsize=6.5)
axB.set_ylabel("profundidade (m)", fontsize=7.5, color=TINTA)
for s in ("top", "right"):
    axB.spines[s].set_visible(False)
axB.set_title("(b) coluna", fontsize=8.5, loc="left", color=TINTA)
axB.annotate("base: u = 0,\nsem fluxo", (W, 0), xytext=(W + 0.03, 0.04), ha="left", fontsize=6.5, color=TINTA2,
             annotation_clip=False)

# (c) o topo ampliado: triangulos, vertices e nos do meio das arestas
for t in tri:
    if min(y for _, y in t) < H - zoom - 1e-9:
        continue
    axC.add_patch(Polygon(t, closed=True, fc=SOLO, ec=TINTA2, lw=0.8))
vert = {(round(x, 6), round(y, 6)) for t in tri for x, y in t if y >= H - zoom - 1e-9}
meio = set()
for t in tri:
    if min(y for _, y in t) < H - zoom - 1e-9:
        continue
    for i in range(3):
        (x1, y1), (x2, y2) = t[i], t[(i + 1) % 3]
        meio.add((round(0.5 * (x1 + x2), 6), round(0.5 * (y1 + y2), 6)))
axC.plot([x for x, _ in meio], [y for _, y in meio], "o", ms=3, color=LARANJA, mec=SUPERF, mew=0.5, zorder=3,
         label="nó do meio da aresta: u (P2)")
axC.plot([x for x, _ in vert], [y for _, y in vert], "o", ms=5.5, color=AZUL, mec=SUPERF, mew=0.7, zorder=4,
         label="vértice: u, $p_w$, $p_a$ (P1)")
axC.plot([0, W], [H, H], color=AZUL, lw=2.2, zorder=2)
for k in range(3):
    x = W * (0.2 + 0.3 * k)
    axC.add_patch(FancyArrowPatch((x, H + 0.03), (x, H + 0.004), arrowstyle="-|>", mutation_scale=8, color=AZUL,
                                  lw=0.9))
axC.text(-0.004, H + 0.034, "topo: chuva q(t), ou $p_w$ = 0 sob lâmina;\n$p_a$ = 0 no modelo TP", fontsize=7,
         color=TINTA, va="bottom", clip_on=False)
fig.text(0.865, 0.52, "laterais:\n$u_x$ = 0,\nsem fluxo\nde água\ne de ar", fontsize=7, color=TINTA2, va="center")
axC.set_xlim(-0.004, W + 0.004)
axC.set_ylim(H - zoom - 0.004, H + 0.032)
axC.set_aspect("equal")
axC.set_xticks([0, W / 2, W])
axC.set_xticklabels(["0", f"{W / 2:g}".replace(".", ","), f"{W:g}".replace(".", ",")], fontsize=6.5)
yt = [H - k * H / ny for k in range(5)]
axC.set_yticks(yt)
axC.set_yticklabels([f"{H - y:.2f}".replace(".", ",") for y in yt], fontsize=6.5)
for s in ("top", "right"):
    axC.spines[s].set_visible(False)
axC.set_title("(c) topo ampliado", fontsize=8.5, loc="left", color=TINTA, pad=30)
for s_ in ("left",):
    axC.spines[s_].set_bounds(H - zoom, H)
axC.legend(fontsize=6.5, frameon=False, loc="upper left", bbox_to_anchor=(-0.25, -0.07), ncol=1)
fig.savefig(saida + ".pdf")
fig.savefig(saida + ".png", dpi=150)
print("pronto:", saida + ".pdf")
