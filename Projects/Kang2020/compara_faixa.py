#!/usr/bin/env python3
"""compara_faixa.py <pasta da faixa> <pasta da coluna> <saida.png>: o FS no angulo medio de cada zona, na faixa do FLAC
(Passo 6, FS_beta), na coluna vertical (Passo 4) e no artigo (lido da Fig. 9), SP e TP."""
import csv
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

AQUI = os.path.dirname(os.path.abspath(__file__))
fx, col, saida = sys.argv[1], sys.argv[2], sys.argv[3]
TINTA, TINTA2, GRADE = "#0b0b0b", "#5c5c58", "#e4e3df"
COR = {"sp": "#e34948", "tp": "#2a78d6"}
lido = [r for r in csv.reader(open(os.path.join(AQUI, "data", "fig9_lido.csv"))) if r and not r[0].startswith("#")]
lido = [dict(zip(lido[0], r)) for r in lido[1:]]


def le(arq, col_):
    if not os.path.exists(arq):
        return [], []
    R = list(csv.DictReader(open(arq)))
    return [float(r["t_h"]) for r in R], [float(r[col_]) for r in R]


fig, axs = plt.subplots(3, 2, figsize=(7.4, 8.4))
axs = axs.ravel()
for z in range(1, 6):
    ax = axs[z - 1]
    ax.grid(True, color=GRADE, lw=0.7)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for m in ("sp", "tp"):
        t, f = le(os.path.join(fx, f"faixa_zona{z}_{m}_fs.csv"), "FS_beta")
        ax.plot(t, f, color=COR[m], lw=1.6, label=f"{m.upper()} faixa (malha do FLAC)")
        t, f = le(os.path.join(col, f"kang_zona{z}_{m}_fs.csv"), "FS_medio")
        ax.plot(t, f, color=COR[m], lw=1.0, ls=(0, (4, 2)), label=f"{m.upper()} coluna")
        L = [r for r in lido if int(r["zona"]) == z and r["modelo"] == m]
        ax.plot([float(r["t_h"]) for r in L], [float(r["FS"]) for r in L], "o" if m == "sp" else "s", ms=4,
                mfc="none", mec=COR[m], mew=1.0, label=f"{m.upper()} Kang (lido)")
    ax.axhline(1, color=TINTA2, lw=0.8, ls=(0, (1, 1.5)))
    ax.set_title(f"Zona {z}", fontsize=9, loc="left", color=TINTA)
    ax.set_xlim(0, 17)
    ax.set_xlabel("t (h)", fontsize=8)
    ax.set_ylabel("FS no ângulo médio", fontsize=8)
    ax.tick_params(labelsize=7)
h, l = axs[0].get_legend_handles_labels()
axs[5].axis("off")
axs[5].legend(h, l, fontsize=7.5, frameon=False, loc="center")
fig.tight_layout()
fig.savefig(saida, dpi=130)
print("pronto:", saida)
