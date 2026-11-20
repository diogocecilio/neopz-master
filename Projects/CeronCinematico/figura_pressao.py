#!/usr/bin/env python3
"""figura_pressao.py: os contornos do excesso de poropressao u (kPa) da agua permanente usados no estudo do contorno:
(a) a celula do notebook darcy-simulation-slope-v2.nb (dominio 44H x 30H x 39H; u = 0 na crista e na lateral esquerda,
face com a coluna d'agua, -gw hw no pe e na lateral direita, base sem fluxo) e (b) o contorno literal do artigo
(tudo Dirichlet, nosso dominio 7H x 6H x 6H). Em cima o dominio inteiro, embaixo o zoom perto do talude.
    python3 -I figura_pressao.py   -> figs/pressao_permanente.pdf/.png"""
import math, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cinematico as k
import estudo_contorno as e

H = 5.0
CASOS = [("(a) notebook (celula colada): base sem fluxo", (44 * H, 30 * H, 39 * H), ("0", "hw", None)),
         ("(b) artigo (dOmega_v vazio): tudo Dirichlet\nnosso dominio 7H x 6H x 6H", (7 * H, 6 * H, 6 * H), ("0", "hw", "0")),
         ("(c) como a Fig. 4b: esq. e base u = 0, dir. sem fluxo\ndominio 13H x 2H x 8H", (13 * H, 2 * H, 8 * H), ("0", None, "0"))]
plt.rcParams.update({"font.size": 10})
fig, ax = plt.subplots(2, 3, figsize=(16, 7.5), layout="constrained")
niveis = np.arange(-50, 0.1, 5.0)
for j, (tit, (Le, Ld, D), cc) in enumerate(CASOS):
    tal, xy, tri = e.malha(Le, Ld, D, 0.0625)
    u = e.permanente(xy, tri, tal, *cc)
    n = len(xy)
    G = k.minimizar(k.CampoMalha(xy, tri, np.full(n, 10.), np.full(n, math.radians(30)), u, 20 - e.GW), tal)
    T = Triangulation(xy[:, 0], -xy[:, 1], tri)      # y para cima no desenho
    for i in range(2):
        a = ax[i, j]
        cs = a.tricontourf(T, u, levels=niveis, cmap="viridis", extend="both")
        a.tricontour(T, u, levels=niveis, colors="k", linewidths=0.4)
        # contorno do dominio, cor por condicao
        xs = [tal.x_esq, 0, tal.xT, tal.x_dir]
        a.plot([tal.x_esq, 0], [0, 0], "r", lw=3)
        a.plot([0, tal.xT], [0, -H], "r", lw=3)
        a.plot([tal.xT, tal.x_dir], [-H, -H], "m", lw=3)
        a.plot([tal.x_esq] * 2, [0, -tal.y_base], "r" if cc[0] == "0" else "0.6", lw=3, ls="-" if cc[0] else ":")
        a.plot([tal.x_dir] * 2, [-H, -tal.y_base], "m" if cc[1] == "hw" else "0.5", lw=3, ls="-" if cc[1] else "--")
        a.plot([tal.x_esq, tal.x_dir], [-tal.y_base] * 2, "r" if cc[2] == "0" else "0.5", lw=3,
               ls="-" if cc[2] else "--")
        a.set_aspect("equal")
        if i == 1:
            m = G.mecanismo
            a.set_xlim(-4 * H, tal.xT + 4 * H); a.set_ylim(-4 * H, 1)
            mx, my = m.pontos()
            a.plot(mx, -my, "w-", lw=2, label=f"mecanismo critico, Gamma = {G.Gamma:.3f}")
            a.legend(loc="lower left", fontsize=9)
            a.set_xlabel("x (m), origem na crista")
        else:
            a.set_title(tit + f"\nGamma cinematico = {G.Gamma:.3f} (artigo 1,336)")
        a.set_ylabel("cota (m)")
cb = fig.colorbar(cs, ax=ax, shrink=0.8, label="excesso de poropressao u = p - gw*y (kPa)")
fig.supxlabel( "bordas: vermelho u = 0 / coluna d'agua na face; magenta u = -gw*hw; cinza tracejado sem fluxo. "
         "H = 5 m, beta = 45, hw = H, k isotropico", fontsize=10)
for ext in ("pdf", "png"):
    fig.savefig(os.path.join("figs", "pressao_permanente." + ext), dpi=130, bbox_inches="tight")
