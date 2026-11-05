#!/usr/bin/env python3
"""plot_mc.py <pasta da combinacao> [pasta deterministica] [saida.png]: FS(t) de cada realizacao (cinza), o caso
deterministico (vermelho) e a media; histograma do FS minimo; P_f(t) = P(FS < 1). Le mc_fs.csv e mc_resumo.csv de
rodar_mc.sh. Requer numpy e matplotlib."""
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

D = sys.argv[1]
det = sys.argv[2] if len(sys.argv) > 2 and os.path.isdir(sys.argv[2]) else None
saida = sys.argv[-1] if sys.argv[-1].endswith(".png") else os.path.join(D, "mc.png")

fs = np.genfromtxt(os.path.join(D, "mc_fs.csv"), delimiter=",", names=True)
res = np.genfromtxt(os.path.join(D, "mc_resumo.csv"), delimiter=",", names=True, dtype=None, encoding=None)
reais = np.unique(fs["real"])
ts = np.unique(fs["t_h"])
fig, ax = plt.subplots(1, 3, figsize=(15, 4.2), constrained_layout=True)
for r in reais:
    s = fs[fs["real"] == r]
    o = np.argsort(s["t_h"])
    ax[0].plot(s["t_h"][o], s["FS_SRM"][o], color="0.6", lw=0.8)
media = [np.mean(fs["FS_SRM"][(fs["t_h"] == t) & (fs["FS_SRM"] > 0)]) for t in ts]
ax[0].plot(ts, media, "k-o", lw=2, ms=4, label=f"media de {len(reais)} realizacoes")
if det:
    d = np.genfromtxt(os.path.join(det, "mc_fs.csv"), delimiter=",", names=True)
    o = np.argsort(d["t_h"])
    ax[0].plot(d["t_h"][o], d["FS_SRM"][o], "r-s", lw=2, ms=4, label="deterministico (medias)")
ax[0].axhline(1, color="k", ls="--", lw=0.8)
ax[0].axvspan(0, 24, color="tab:blue", alpha=0.08, label="chuva")
ax[0].set_xlabel("t (h)")
ax[0].set_ylabel("FS (SRM)")
ax[0].legend(fontsize=8)
fsmin = np.atleast_1d(res["FS_min"])
ax[1].hist(fsmin[fsmin > 0], bins=max(5, int(np.sqrt(len(fsmin)))), color="0.5", edgecolor="k")
ax[1].axvline(1, color="k", ls="--", lw=0.8)
ax[1].set_xlabel("FS minimo de cada realizacao")
ax[1].set_ylabel("realizacoes")
pf = [np.mean(fs["FS_SRM"][fs["t_h"] == t] < 1) for t in ts]
ax[2].plot(ts, pf, "k-o", ms=4)
ax[2].set_xlabel("t (h)")
ax[2].set_ylabel("P_f = P(FS < 1)")
ax[2].set_ylim(bottom=0)
fig.suptitle(os.path.basename(os.path.normpath(D)))
fig.savefig(saida, dpi=110)
print(saida)
