#!/usr/bin/env python3
"""figura_pz_contorno_ceron.py: os contornos do excesso de poropressao u = p_w - gw (35 - y) do PZ com a opcao
contorno=ceron (os contornos das Eqs. 19-21 do Ceron: crista p_w = 0, face e terreno do pe com o nivel final, lateral
esquerda e base u = 0, lateral direita sem fluxo) no rebaixamento instantaneo (trebaix = 1e-4 h), deterministico:
agua saturada em t = 10 h e nao saturada em t = 10 h e t = 100 h. Duas linhas: o dominio inteiro e o zoom perto do
talude (x de -20 a 25 m a partir da crista). No titulo: o Gamma do GIM do PZ (u-pw-pa, 3 niveis), o u-p (so no sat,
se houver) e o Gamma cinematico (cinematico.py) com as mesmas pressoes. A referencia: o permanente P1 do
estudo_contorno.py (caso D) na malha da exportacao, com o erro max |u_PZ - u_perm|.
    python3 -I figura_pz_contorno_ceron.py [pasta das rodadas]   (padrao ~/projects/runs/ceron_contorno)
    -> figs/pz_contorno_ceron.png/.pdf e dados/pz_contorno_ceron.csv"""
import csv
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cinematico as k  # noqa: E402
import estudo_contorno as ec  # noqa: E402

base = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else os.path.expanduser("~/projects/runs/ceron_contorno"))
aqui = os.path.dirname(os.path.abspath(__file__))
CASOS = [("sat_ceron", "agua saturada, t = 10 h", "up_sat_ceron.csv"),
         ("unsat_ceron_t10", "agua nao saturada, t = 10 h", None),
         ("unsat_ceron_t100", "agua nao saturada, t = 100 h", None)]
plt.rcParams.update({"font.size": 11})


def gim(arq):
    with open(arq) as f:
        return float(next(csv.DictReader(f))["FS_min"])


def gim_up(arq):
    with open(arq) as f:
        r = list(csv.DictReader(f))
    return float(r[-1]["Gamma"])


niv = np.arange(-50, 0.1, 5.0)
fig, ax = plt.subplots(2, 3, figsize=(17, 6.6), layout="constrained", height_ratios=[1, 1.2])
linhas = []
cs = None
for j, (s, nome, arq_up) in enumerate(CASOS):
    pasta = os.path.join(base, s + "_real-1")
    if not os.path.isdir(pasta):
        for i in range(2):
            ax[i, j].set_title(f"{nome}: falta rodar")
            ax[i, j].axis("off")
        continue
    tal, campo = k.ler_exportacao(pasta)
    V = np.loadtxt(os.path.join(pasta, "vertices.csv"), delimiter=",", skiprows=1)
    tri = np.loadtxt(os.path.join(pasta, "triangles.csv"), delimiter=",", skiprows=1).astype(int)
    u = campo.un
    g_pz = gim(os.path.join(base, s + "_0_-1.csv"))
    g_up = gim_up(os.path.join(base, arq_up)) if arq_up and os.path.exists(os.path.join(base, arq_up)) else None
    g_cin = k.minimizar(campo, tal).Gamma
    # referencia: permanente saturado (caso D) na mesma malha, no sistema do artigo (y para baixo a partir da crista)
    xy = np.c_[V[:, 1] - 35.0, 35.0 - V[:, 2]]
    u_perm = ec.permanente(xy, tri, tal, "0", None, "0")
    erro = float(np.max(np.abs(u - u_perm)))
    linhas.append({"caso": s, "Gamma_GIM_upwpa": g_pz, "Gamma_GIM_up": g_up if g_up else "", "Gamma_cin": g_cin,
                   "pw_min": campo.pw.min(), "pw_max": campo.pw.max(), "u_min": u.min(), "u_max": u.max(),
                   "erro_max_u_permanente": erro})
    T = Triangulation(V[:, 1] - 35.0, V[:, 2] - 35.0, tri)
    tit = f"{nome}\nGIM u-pw-pa {g_pz:.3f}"
    if g_up:
        tit += f" | GIM u-p {g_up:.3f}"
    tit += f" | cinematico {g_cin:.3f}"
    for i in range(2):
        a = ax[i, j]
        cs = a.tricontourf(T, u, levels=niv, cmap="viridis", extend="both")
        a.tricontour(T, u, levels=niv, colors="k", linewidths=0.4)
        a.plot([-35, 0, 5, 35], [0, 0, -5, -5], "r", lw=1.5)
        a.set_aspect("equal")
        if i == 0:
            a.set_xlim(-35, 35); a.set_ylim(-35, 1)
            a.set_title(tit, fontsize=11)
        else:
            a.set_xlim(-20, 25); a.set_ylim(-20, 1)
            a.set_title(f"zoom; max |u - u permanente| = {erro:.2f} kPa", fontsize=10)
            a.set_xlabel("x (m), origem na crista")
        if j == 0:
            a.set_ylabel("cota (m), origem na crista")
fig.colorbar(cs, ax=ax, shrink=0.8, label="excesso de poropressao u = p_w - gw (35 - y) (kPa)")
fig.suptitle("PZ com contorno=ceron (Eqs. 19-21), rebaixamento instantaneo (H = 5 m, hw = H), deterministico, "
             "3 niveis; o artigo da 1,336", fontsize=12)
for ext in ("png", "pdf"):
    fig.savefig(os.path.join(aqui, "figs", "pz_contorno_ceron." + ext), dpi=110)
if linhas:
    with open(os.path.join(aqui, "dados", "pz_contorno_ceron.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(linhas[0]))
        w.writeheader()
        w.writerows(linhas)
for r in linhas:
    print(r)
