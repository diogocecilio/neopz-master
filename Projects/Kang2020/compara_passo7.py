#!/usr/bin/env python3
"""compara_passo7.py <pasta do passo7> <pasta do passo4> <saida.png>: o FS no angulo medio de cada zona na coluna normal
ao talude (Passo 7, normal_zona<z>_<m>_fs.csv, FS_beta), na faixa do FLAC com a chuva por metro horizontal (Passo 7,
faixah_zona<z>_<m>_fs.csv, FS_beta), na coluna vertical (Passo 4, kang_zona<z>_<m>_fs.csv, FS_medio) e no artigo (lido
da Fig. 9), sem ar (SP) e com ar ativo (TP). Imprime tambem a tabela do FS final e a maior diferenca entre a coluna
normal e a faixa ao longo das 17 h."""
import csv
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

AQUI = os.path.dirname(os.path.abspath(__file__))
p7, p4, saida = sys.argv[1], sys.argv[2], sys.argv[3]
TINTA, TINTA2, GRADE = "#0b0b0b", "#5c5c58", "#e4e3df"
COR = {"sp": "#e34948", "tp": "#2a78d6"}
lido = [r for r in csv.reader(open(os.path.join(AQUI, "data", "fig9_lido.csv"))) if r and not r[0].startswith("#")]
lido = [dict(zip(lido[0], r)) for r in lido[1:]]


def le(arq, col):
    if not os.path.exists(arq):
        return [], []
    R = list(csv.DictReader(open(arq)))
    return [float(r["t_h"]) for r in R], [float(r[col]) for r in R]


fig, axs = plt.subplots(3, 2, figsize=(7.4, 8.4))
axs = axs.ravel()
print("zona modelo | FS final: normal  faixa(h)  coluna  Kang | max |normal - faixa(h)| em 0-17 h")
for z in range(1, 6):
    ax = axs[z - 1]
    ax.grid(True, color=GRADE, lw=0.7)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for m in ("sp", "tp"):
        tn, fn = le(os.path.join(p7, f"normal_zona{z}_{m}_fs.csv"), "FS_beta")
        tf, ff = le(os.path.join(p7, f"faixah_zona{z}_{m}_fs.csv"), "FS_beta")
        tc, fc = le(os.path.join(p4, f"kang_zona{z}_{m}_fs.csv"), "FS_medio")
        ax.plot(tn, fn, color=COR[m], lw=1.8, label=f"{m.upper()} coluna normal")
        ax.plot(tf, ff, color=COR[m], lw=1.2, ls=(0, (4, 2)), label=f"{m.upper()} faixa, chuva por m horizontal")
        ax.plot(tc, fc, color=COR[m], lw=0.8, ls=(0, (1, 1.2)), label=f"{m.upper()} coluna vertical")
        L = [r for r in lido if int(r["zona"]) == z and r["modelo"] == m]
        ax.plot([float(r["t_h"]) for r in L], [float(r["FS"]) for r in L], "o" if m == "sp" else "s", ms=4,
                mfc="none", mec=COR[m], mew=1.0, label=f"{m.upper()} Kang (lido)")
        kf = [float(r["FS"]) for r in L if float(r["t_h"]) == 17.]
        dmax = max((abs(a - b) for a, b in zip(fn, ff)), default=float("nan")) if fn and ff and len(fn) == len(ff) \
            else float("nan")
        print(f"  {z}  {m}   |  {fn[-1] if fn else float('nan'):6.3f}  {ff[-1] if ff else float('nan'):6.3f}  "
              f"{fc[-1] if fc else float('nan'):6.3f}  {kf[0] if kf else float('nan'):5.2f} | {dmax:.3f}")
    ax.axhline(1, color=TINTA2, lw=0.8, ls=(0, (1, 1.5)))
    ax.set_title(f"Zona {z}", fontsize=9, loc="left", color=TINTA)
    ax.set_xlim(0, 17)
    ax.set_xlabel("t (h)", fontsize=8)
    ax.set_ylabel("FS no ângulo médio", fontsize=8)
    ax.tick_params(labelsize=7)
h, l = axs[0].get_legend_handles_labels()
axs[5].axis("off")
axs[5].legend(h, l, fontsize=7, frameon=False, loc="center")
fig.tight_layout()
fig.savefig(saida, dpi=130)
print("pronto:", saida)
