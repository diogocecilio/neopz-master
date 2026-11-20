#!/usr/bin/env python3
"""figura_instantaneo.py: os contornos do excesso de poropressao u = p_w - gw (35 - y) do PZ no rebaixamento INSTANTANEO
(trebaix = 1e-4 h) do talude do Ceron, deterministico, com a agua saturada (agua=sat, a do u-p) e a nao saturada
(agua=unsat, a do u-pw-pa), em t = 0+, 0,1, 1 e 10 h. No titulo: o Gamma do GIM do PZ (u-pw-pa e u-p, 3 niveis) e o
Gamma cinematico (cinematico.py) com as mesmas pressoes.
    python3 -I figura_instantaneo.py <pasta das rodadas>   (~/projects/runs/ceron_instantaneo; rodar.sh)
    -> figs/instantaneo_pz.png/.pdf e dados/instantaneo_pz.csv"""
import csv, glob, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cinematico as k

base = os.path.abspath(sys.argv[1])
aqui = os.path.dirname(os.path.abspath(__file__))
TEMPOS = ["0.0001", "0.1", "1", "10"]
ROTULO = {"0.0001": "t = 0+ (0,36 s)", "0.1": "t = 0,1 h", "1": "t = 1 h", "10": "t = 10 h"}
AGUAS = [("sat", "agua saturada (u-p)"), ("unsat", "agua nao saturada (u-pw-pa)")]


def gim(arq):
    with open(arq) as f:
        return float(next(csv.DictReader(f))["FS_min"])


def gim_up(arq):
    with open(arq) as f:
        r = list(csv.DictReader(f))
    return float(r[-1]["Gamma"])


niv = np.arange(-50, 0.1, 5.0)
fig, ax = plt.subplots(2, 4, figsize=(18, 6.8), layout="constrained")
linhas = []
for i, (agua, nome) in enumerate(AGUAS):
    for j, t in enumerate(TEMPOS):
        a = ax[i, j]
        s = f"{agua}_t{t}"
        pasta = os.path.join(base, s + "_real-1")
        if not os.path.isdir(pasta):
            a.set_title(f"{nome}\n{ROTULO[t]}: falta rodar"); a.axis("off"); continue
        tal, campo = k.ler_exportacao(pasta)
        V = np.loadtxt(os.path.join(pasta, "vertices.csv"), delimiter=",", skiprows=1)
        tri = np.loadtxt(os.path.join(pasta, "triangles.csv"), delimiter=",", skiprows=1).astype(int)
        u = campo.un
        g_pz = gim(os.path.join(base, s + "_0_-1.csv"))
        arq_up = os.path.join(base, f"up_{s}.csv")
        g_up = gim_up(arq_up) if os.path.exists(arq_up) else float("nan")
        g_cin = k.minimizar(campo, tal).Gamma
        linhas.append({"agua": agua, "t_h": t, "Gamma_GIM_upwpa": g_pz, "Gamma_GIM_up": g_up, "Gamma_cin": g_cin,
                       "u_min": u.min(), "u_max": u.max()})
        T = Triangulation(V[:, 1] - 35.0, V[:, 2] - 35.0, tri)
        cs = a.tricontourf(T, u, levels=niv, cmap="viridis", extend="both")
        a.tricontour(T, u, levels=niv, colors="k", linewidths=0.4)
        a.plot([-20, 0, 5, 25], [0, 0, -5, -5], "r", lw=1.5)
        a.set_xlim(-20, 25); a.set_ylim(-20, 1); a.set_aspect("equal")
        a.set_title(f"{nome}, {ROTULO[t]}\nGIM u-pw-pa {g_pz:.3f} | GIM u-p {g_up:.3f} | cinematico {g_cin:.3f}",
                    fontsize=10)
        if i == 1: a.set_xlabel("x (m), origem na crista")
        if j == 0: a.set_ylabel("cota (m)")
fig.colorbar(cs, ax=ax, shrink=0.8, label="excesso de poropressao u = p_w - gw (35 - y) (kPa)")
fig.suptitle("PZ, rebaixamento instantaneo (H = 5 m, hw = H), deterministico; o artigo (agua permanente) da 1,336",
             fontsize=12)
for ext in ("png", "pdf"):
    fig.savefig(os.path.join(aqui, "figs", "instantaneo_pz." + ext), dpi=110)
with open(os.path.join(aqui, "dados", "instantaneo_pz.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(linhas[0])); w.writeheader(); w.writerows(linhas)
for r in linhas:
    print(r)
