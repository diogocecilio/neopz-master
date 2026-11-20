#!/usr/bin/env python3
"""fig4b_ceron.py: tenta reproduzir os contornos do excesso de poropressao da Fig. 4b do Ceron (2025), por cima da
propria figura do artigo (pagina 8 rasterizada a 300 dpi), para descobrir o dominio e o contorno da agua.

Medidas tiradas da figura (pixels): crista em x = 795, pe em x = 810 (H = 15 px, 45 graus), moldura de x = 24 a 963
e de y = 106 a 583 -> dominio ~ 50H atras da crista, 10H a frente do pe e 30H abaixo do pe. Os contornos da figura
sao -9,81 ... -0,98 (fracoes 1,0 ... 0,1 de gw hw, com hw = H = 1).
    python3 -I fig4b_ceron.py <pagina8.png>   -> figs/fig4b_reproducao.png/.pdf"""
import math, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
from PIL import Image
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cinematico as k
import estudo_contorno as e

pag = Image.open(sys.argv[1]).convert("L")
W, Hp = pag.size
x0, y0 = int(W * 0.51), int(Hp * 0.19)
img = np.array(pag)[y0:int(Hp * 0.42), x0:W - int(W * 0.09)]
FX0, FX1, FY0, FY1 = 24, 963, 106, 583          # moldura do painel (b) no recorte
Le, Ld, D = 50.0, 10.0, 30.0                       # em H
CASOS = {"B: esq u=0, dir u=-gw hw, base sem fluxo": ("0", "hw", None),
         "C: tudo Dirichlet (esq 0, base 0, dir -gw hw)": ("0", "hw", "0"),
         "D: esq 0, base 0, dir sem fluxo": ("0", None, "0"),
         "E: esq 0, base e dir sem fluxo": ("0", None, None)}
niv = np.sort(-9.81 * np.arange(1, 10) / 10.0)
fig, ax = plt.subplots(2, 2, figsize=(15, 9), layout="constrained")
res = []
for a, (nome, cc) in zip(ax.ravel(), CASOS.items()):
    # contornos (H = 1, como a figura)
    e.H = 1.0
    tal, xy, tri = e.malha(Le, Ld, D, 0.05)
    u = e.permanente(xy, tri, tal, *cc, hw=1.0)
    # Gamma do caso do artigo (H = 5 m, mesmo dominio em H)
    e.H = 5.0
    tal5, xy5, tri5 = e.malha(5 * Le, 5 * Ld, 5 * D, 0.0625)
    u5 = e.permanente(xy5, tri5, tal5, *cc, hw=5.0)
    n5 = len(xy5)
    G = k.minimizar(k.CampoMalha(xy5, tri5, np.full(n5, 10.), np.full(n5, math.radians(30)), u5, 20 - e.GW), tal5).Gamma
    res.append((nome, G))
    # mapa pixels <- coordenadas: x de -Le a 1+Ld em FX0..FX1, y (para baixo) de 0 a 1+D em FY0..FY1
    sx = (FX1 - FX0) / (Le + 1 + Ld); sy = (FY1 - FY0) / (1 + D)
    a.imshow(img, cmap="gray", extent=(-Le - FX0 / sx, -Le + (img.shape[1] - FX0) / sx,
                                       -(img.shape[0] - FY0) / sy, FY0 / sy), alpha=0.6)
    T = Triangulation(xy[:, 0], -xy[:, 1], tri)
    a.tricontour(T, u, levels=niv, colors="r", linewidths=1.2)
    a.set_xlim(-Le - 1, 1 + Ld + 1); a.set_ylim(-1 - D - 1, 2)
    a.set_title(f"{nome}\nGamma (H = 5 m, mesmo dominio) = {G:.3f}  (artigo 1,336)", fontsize=11)
fig.suptitle("Fig. 4b do Ceron (cinza) x os nossos contornos (vermelho), u = -0,1 ... -0,9 gw hw; "
             "dominio 50H x 10H x 30H", fontsize=12)
for ext in ("png", "pdf"):
    fig.savefig(os.path.join("figs", "fig4b_reproducao." + ext), dpi=120)
for r in res:
    print(f"{r[0]:50s} Gamma = {r[1]:.4f}")
