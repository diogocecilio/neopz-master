#!/usr/bin/env python3
"""fig4c_ceron.py: as forcas de percolacao -grad u (Fig. 4c do Ceron 2025) com o contorno que reproduz a Fig. 4b
(lateral esquerda e base u = 0, lateral direita sem fluxo; dominio 14H x 2,8H x 8,4H, o que da Gamma = 1,337), por
cima da figura do artigo (pagina 8 a 300 dpi) e sozinhas.
    python3 -I fig4c_ceron.py <pagina8.png>   -> figs/fig4c_reproducao.png/.pdf"""
import os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
from PIL import Image
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import estudo_contorno as e

pag = Image.open(sys.argv[1]).convert("L")
W, Hp = pag.size
img = np.array(pag)[int(Hp * 0.37):int(Hp * 0.53), int(W * 0.51):W - int(W * 0.07)]
PX_CRISTA, PY_CRISTA, S = 445.0, 129.0, 262.5        # crista e H em pixels no recorte (face a 45 graus)

e.H = 1.0
tal, xy, tri = e.malha(14.0, 2.8, 8.4, 0.02)
u = e.permanente(xy, tri, tal, "0", None, "0", hw=1.0)       # excesso de poropressao (gw hw = 9,81)
T = Triangulation(xy[:, 0], -xy[:, 1], tri)                   # y para cima
P = xy[tri]
M = np.concatenate([np.ones((len(tri), 3, 1)), P[:, :, :1], -P[:, :, 1:]], axis=2)
coef = np.linalg.solve(M, u[tri][..., None])[..., 0]                           # u = a + bx x + by y (y para cima)
fx, fy = -coef[:, 1], -coef[:, 2]                             # -grad u (kN/m3)
finder = T.get_trifinder()
gx, gy = np.meshgrid(np.arange(-1.6, 2.1, 0.07), np.arange(-1.5, 0.0, 0.07))
it = finder(gx, gy)
ok = it >= 0
X, Y = gx[ok], gy[ok]
U, V = fx[it[ok]], fy[it[ok]]

fig, ax = plt.subplots(1, 2, figsize=(16, 4.6), layout="constrained")
ext = (-PX_CRISTA / S, (img.shape[1] - PX_CRISTA) / S, -(img.shape[0] - PY_CRISTA) / S, PY_CRISTA / S)
ax[0].imshow(img, cmap="gray", extent=ext)
esp = (np.round(X / 0.07) % 2 == 0) & (np.round(Y / 0.07) % 2 == 0)     # metade das setas no painel sobreposto
ax[0].quiver(X[esp], Y[esp], U[esp], V[esp], color="r", alpha=0.6, angles="xy", scale_units="xy", scale=110,
             width=0.0018, headwidth=4)
ax[1].quiver(X, Y, U, V, color="k", angles="xy", scale_units="xy", scale=110, width=0.0016, headwidth=4)
for a in ax:
    a.plot([-1.7, 0, 1, 2.2], [0, 0, -1, -1], "k", lw=1)
    a.set_xlim(-1.65, 2.2); a.set_ylim(-1.5, 0.15); a.set_aspect("equal")
    a.set_xlabel("x / H (origem na crista)"); a.set_ylabel("cota / H")
ax[0].set_title("Fig. 4c do Ceron (cinza) x a nossa -grad u (vermelho)")
ax[1].set_title("a nossa -grad u: esq. e base u = 0, dir. sem fluxo, dominio 14H x 2,8H x 8,4H")
for ext_ in ("png", "pdf"):
    fig.savefig(os.path.join("figs", "fig4c_reproducao." + ext_), dpi=120)
