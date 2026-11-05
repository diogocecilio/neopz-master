#!/usr/bin/env python3
"""plot_mc.py <pasta da combinacao> [pasta deterministica] [saida.png]: junta as realizacoes concluidas (junta_mc.sh) e
desenha FS(t) de cada realizacao (cinza), o caso deterministico (vermelho) e a media; o histograma do FS minimo;
P_f(t) = P(FS < 1). Em <saida>_correlacoes.png: o FS minimo contra as medias dos campos no solo e na banda de ruptura
(com a correlacao r de Pearson). Le mc_fs.csv e mc_resumo.csv. Requer numpy e matplotlib."""
import os
import subprocess
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

D = sys.argv[1]
det = sys.argv[2] if len(sys.argv) > 2 and os.path.isdir(sys.argv[2]) else None
saida = sys.argv[-1] if sys.argv[-1].endswith(".png") else os.path.join(D, "mc.png")
junta = os.path.join(os.path.dirname(os.path.abspath(__file__)), "junta_mc.sh")
for d in [D] + ([det] if det else []):
    if os.path.exists(junta):
        subprocess.run([junta, d], check=False, stdout=subprocess.DEVNULL)


def le(arq):
    return np.genfromtxt(arq, delimiter=",", names=True, dtype=None, encoding=None)


fs = le(os.path.join(D, "mc_fs.csv"))
res = le(os.path.join(D, "mc_resumo.csv"))
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
    d = le(os.path.join(det, "mc_fs.csv"))
    o = np.argsort(np.atleast_1d(d["t_h"]))
    ax[0].plot(np.atleast_1d(d["t_h"])[o], np.atleast_1d(d["FS_SRM"])[o], "r-s", lw=2, ms=4, label="deterministico (medias)")
ax[0].axhline(1, color="k", ls="--", lw=0.8)
ax[0].axvspan(0, 24, color="tab:blue", alpha=0.08, label="chuva")
ax[0].set_xlabel("t (h)")
ax[0].set_ylabel("FS (SRM)")
ax[0].legend(fontsize=8)
fsmin = np.atleast_1d(res["FS_min"]).astype(float)
ax[1].hist(fsmin[fsmin > 0], bins=max(5, int(np.sqrt(len(fsmin)))), color="0.5", edgecolor="k")
# LIMITE_FS (padrao 1): o limiar de P_f (1.01 nas variantes acopladas com F_TRANS = 1, em que o FS nao cai abaixo de 1)
LIM = float(os.environ.get("LIMITE_FS", "1"))
ax[1].axvline(LIM, color="k", ls="--", lw=0.8)
ax[1].set_xlabel("FS minimo de cada realizacao")
ax[1].set_ylabel("realizacoes")
pf = [np.mean(fs["FS_SRM"][fs["t_h"] == t] < LIM) for t in ts]
ax[2].plot(ts, pf, "k-o", ms=4)
ax[2].set_xlabel("t (h)")
ax[2].set_ylabel(f"P_f = P(FS < {LIM:g})")
ax[2].set_ylim(bottom=0)
fig.suptitle(os.path.basename(os.path.normpath(D)))
fig.savefig(saida, dpi=110)
print(saida)

# FS minimo contra as propriedades medias: no solo todo e na banda de ruptura do instante do FS minimo
pares = [("c_med", "c' medio no solo (kPa)"), ("phi_med", "phi' medio no solo (graus)"), ("k_mult_med", "k/k_s medio no solo"),
         ("banda_c", "c' medio na banda (kPa)"), ("banda_phi_deg", "phi' medio na banda (graus)"),
         ("banda_k_mult", "k/k_s medio na banda")]
pares = [p for p in pares if p[0] in (res.dtype.names or ())]
if pares and len(fsmin) >= 3:
    fig, ax = plt.subplots(2, 3, figsize=(13, 7.5), constrained_layout=True)
    for a, (col, rotulo) in zip(ax.flat, pares):
        x = np.atleast_1d(res[col]).astype(float)
        ok = np.isfinite(x) & (fsmin > 0)
        a.scatter(x[ok], fsmin[ok], s=12, color="k")
        a.axhline(1, color="k", ls="--", lw=0.8)
        r = np.corrcoef(x[ok], fsmin[ok])[0, 1] if ok.sum() >= 3 else np.nan
        a.set_title(f"r = {r:.2f}")
        a.set_xlabel(rotulo)
        a.set_ylabel("FS minimo")
    for a in ax.flat[len(pares):]:
        a.axis("off")
    fig.suptitle(os.path.basename(os.path.normpath(D)))
    s2 = saida.replace(".png", "_correlacoes.png")
    fig.savefig(s2, dpi=110)
    print(s2)
