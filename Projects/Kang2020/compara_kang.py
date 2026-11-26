#!/usr/bin/env python3
"""compara_kang.py <pasta dos CSV> [prefixo=kang]: as Figs. 7, 8 e 9 de Kang et al. (2020) com as nossas contas.
  - fig9_nosso.png: o FS minimo (Eq. 14) contra o tempo nas 5 zonas, nos tres angulos (max, medio, min), SP (tracejado)
    e TP (cheio) quando existirem, com o FS do angulo medio LIDO da Fig. 9 do artigo (marcadores; data/fig9_lido.csv);
  - fig7_nosso.png (SP) e fig8_nosso.png (TP): os perfis de p_w, (p_a) e succao a 0, 3, 6, 9, 12, 15 e 17 h."""
import csv
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

AQUI = os.path.dirname(os.path.abspath(__file__))
pasta = sys.argv[1] if len(sys.argv) > 1 else "."
pre = sys.argv[2] if len(sys.argv) > 2 else "kang"
TINTA, TINTA2, GRADE = "#0b0b0b", "#5c5c58", "#e4e3df"
COR_ANG = ["#e34948", "#2a78d6", "#1baf7a"]  # max, medio, min (vermelho, azul, verde)
HORAS = [0, 3, 6, 9, 12, 15, 17]
COR_T = plt.cm.viridis([i / (len(HORAS) - 1) for i in range(len(HORAS))])


def le(arq):
    return list(csv.DictReader(open(arq))) if os.path.exists(arq) else []


def estilo(ax):
    ax.grid(True, color=GRADE, lw=0.7)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.tick_params(labelsize=7, colors=TINTA2)


lido = [r for r in csv.reader(open(os.path.join(AQUI, "data", "fig9_lido.csv"))) if r and not r[0].startswith("#")]
lido = [dict(zip(lido[0], r)) for r in lido[1:]]

# ---- Fig. 9
fig, axs = plt.subplots(3, 2, figsize=(7.4, 8.4))
axs = axs.ravel()
for z in range(1, 6):
    ax = axs[z - 1]
    estilo(ax)
    for mod, ls in (("sp", (0, (4, 2))), ("tp", "-")):
        fs = le(os.path.join(pasta, f"{pre}_zona{z}_{mod}_fs.csv"))
        if not fs:
            continue
        t = [float(r["t_h"]) for r in fs]
        for k, (col, nome) in enumerate(zip(("FS_max", "FS_medio", "FS_min"), ("max", "médio", "min"))):
            ax.plot(t, [float(r[col]) for r in fs], color=COR_ANG[k], ls=ls, lw=1.3,
                    label=f"{mod.upper()} {nome} (nosso)")
        L = [r for r in lido if int(r["zona"]) == z and r["modelo"] == mod]
        ax.plot([float(r["t_h"]) for r in L], [float(r["FS"]) for r in L], "o" if mod == "sp" else "s", ms=4,
                mfc="none", mec=TINTA, mew=0.9, label=f"{mod.upper()} médio (Kang, lido)")
    ax.axhline(1, color=TINTA2, lw=0.8, ls=(0, (1, 1.5)))
    ax.set_title(f"Zona {z}", fontsize=9, loc="left", color=TINTA)
    ax.set_xlim(0, 17)
    ax.set_xlabel("t (h)", fontsize=8)
    ax.set_ylabel("FS mínimo", fontsize=8)
h, l = axs[0].get_legend_handles_labels()
axs[5].axis("off")
axs[5].legend(h, l, fontsize=7, frameon=False, loc="center")
fig.tight_layout()
fig.savefig(os.path.join(pasta, "fig9_nosso.png"), dpi=130)

# ---- Figs. 7 (SP) e 8 (TP): os perfis
for mod, nome, cols in (("sp", "fig7_nosso.png", ("pw_kPa", "succao_kPa")),
                        ("tp", "fig8_nosso.png", ("pw_kPa", "pa_kPa", "succao_kPa"))):
    if not any(os.path.exists(os.path.join(pasta, f"{pre}_zona{z}_{mod}_perfis.csv")) for z in range(1, 6)):
        continue
    fig, axs = plt.subplots(5, len(cols), figsize=(2.5 * len(cols) + 0.8, 11), squeeze=False)
    for z in range(1, 6):
        P = le(os.path.join(pasta, f"{pre}_zona{z}_{mod}_perfis.csv"))
        for j, c in enumerate(cols):
            ax = axs[z - 1][j]
            estilo(ax)
            for hh, cor in zip(HORAS, COR_T):
                R = [r for r in P if int(float(r["t_h"])) == hh]
                if R:
                    ax.plot([float(r[c]) for r in R], [float(r["z_m"]) for r in R], color=cor, lw=1.2, label=f"{hh} h")
            ax.set_ylim(2, 0)
            ax.set_xlim(-20 if c == "pw_kPa" else 0, 20)
            if z == 1:
                ax.set_title({"pw_kPa": "p_w (kPa)", "pa_kPa": "p_a (kPa)", "succao_kPa": "sucção (kPa)"}[c], fontsize=9)
            if j == 0:
                ax.set_ylabel(f"Zona {z}\nprofundidade (m)", fontsize=8)
    axs[0][-1].legend(fontsize=7, frameon=False, loc="lower right")
    fig.tight_layout()
    fig.savefig(os.path.join(pasta, nome), dpi=130)
print("pronto:", ", ".join(f for f in ("fig9_nosso.png", "fig7_nosso.png", "fig8_nosso.png")
                           if os.path.exists(os.path.join(pasta, f))))
