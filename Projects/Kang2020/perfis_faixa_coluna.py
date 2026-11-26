#!/usr/bin/env python3
"""perfis_faixa_coluna.py <pasta da coluna> <pasta da faixa> <modelo sp|tp> <saida.png>: os perfis de p_w, (p_a) e
succao a 0, 3, 6, 9, 12, 15 e 17 h nas 5 zonas, na coluna (linhas cheias; kang_zona<z>_<m>_perfis.csv) e na faixa do
FLAC, na vertical do meio (tracejadas; faixa_zona<z>_<m>_perfis.csv)."""
import csv
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

col, fx, mod, saida = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
TINTA, TINTA2, GRADE = "#0b0b0b", "#5c5c58", "#e4e3df"
HORAS = [0, 3, 6, 9, 12, 15, 17]
COR = plt.cm.viridis([i / (len(HORAS) - 1) for i in range(len(HORAS))])
CAMPOS = ("pw_kPa", "pa_kPa", "succao_kPa") if mod == "tp" else ("pw_kPa", "succao_kPa")
TIT = {"pw_kPa": "$p_w$ (kPa)", "pa_kPa": "$p_a$ (kPa)", "succao_kPa": "sucção (kPa)"}


def le(arq):
    return list(csv.DictReader(open(arq))) if os.path.exists(arq) else []


fig, axs = plt.subplots(5, len(CAMPOS), figsize=(2.6 * len(CAMPOS) + 0.6, 11.2), squeeze=False)
for z in range(1, 6):
    C = le(os.path.join(col, f"kang_zona{z}_{mod}_perfis.csv"))
    F = le(os.path.join(fx, f"faixa_zona{z}_{mod}_perfis.csv"))
    for j, c in enumerate(CAMPOS):
        ax = axs[z - 1][j]
        ax.grid(True, color=GRADE, lw=0.7)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        ax.tick_params(labelsize=7, colors=TINTA2)
        for hh, cor in zip(HORAS, COR):
            for R, ls, lw in ((C, "-", 1.3), (F, (0, (3, 1.6)), 1.3)):
                P = [r for r in R if int(float(r["t_h"])) == hh]
                if P:
                    ax.plot([float(r[c]) for r in P], [float(r["z_m"]) for r in P], color=cor, ls=ls, lw=lw)
        ax.set_ylim(2, 0)
        ax.set_xlim(-20 if c == "pw_kPa" else 0, 20)
        if z == 1:
            ax.set_title(TIT[c], fontsize=9, color=TINTA)
        if j == 0:
            ax.set_ylabel(f"Zona {z}\nprofundidade (m)", fontsize=8, color=TINTA)
alca = [Line2D([], [], color=cor, lw=1.6) for cor in COR] + [Line2D([], [], color=TINTA2, lw=1.3),
                                                            Line2D([], [], color=TINTA2, lw=1.3, ls=(0, (3, 1.6)))]
rot = [f"{h} h" for h in HORAS] + ["coluna", "faixa do FLAC (x = 5 m)"]
fig.legend(alca, rot, loc="upper center", ncol=9, fontsize=7, frameon=False, bbox_to_anchor=(0.5, 1.0),
           handlelength=2.2, columnspacing=1.0)
fig.tight_layout(rect=(0, 0, 1, 0.975))
fig.savefig(saida, dpi=130)
print("pronto:", saida)
