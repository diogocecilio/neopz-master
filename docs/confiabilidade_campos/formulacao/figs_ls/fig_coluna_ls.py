# Figura do LS na coluna de Cho (2014), 18 h (24/11/2026): direcoes do piloto, truncada e em estagios, contra o MC em Python.
import csv, math, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import norm
D = "/home/diogo/projects/runs/coluna_ls_piloto/"
runs = [("ls_p500_T18", "piloto de 500 (6 falhas)", "tab:gray"),
        ("ls_p5000_T18_rec", "piloto de 5000", "tab:blue"),
        ("ls_p5000tr10_T18", "piloto de 5000, modos 0-9", "tab:green"),
        ("ls_est2_50a199_T18", "em estagios: 2o estagio", "tab:red")]
pf_mc, se_mc = 0.0148, 0.0017077
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
cpl = []
for pasta, rot, cor in runs:
    r = list(csv.DictReader(open(D + pasta + "/linhas.csv")))
    r.sort(key=lambda x: int(x["linha"]))
    p = np.array([float(x["p"]) for x in r]); an = np.array([float(x["analises"]) for x in r])
    k = np.arange(1, len(p) + 1)
    ax[0].plot(np.cumsum(an), np.cumsum(p) / k, color=cor, label=rot)
    cpl.append((rot, p.std(ddof=1) / p.mean(), cor))
ax[0].axhline(pf_mc, color="k", lw=1)
ax[0].axhspan(pf_mc - 1.96 * se_mc, pf_mc + 1.96 * se_mc, color="0.85", zorder=0, label="MC em Python (5000), IC 95%")
ax[0].set_xscale("log"); ax[0].set_ylim(0, 0.05)
ax[0].set_xlabel("análises do SlopeRainColumn"); ax[0].set_ylabel("$P_f$ (18 h)")
ax[0].set_title("(a) a estimativa ao longo das linhas"); ax[0].legend(fontsize=8, loc="upper right")
ax[1].bar(range(len(cpl)), [c[1] for c in cpl], color=[c[2] for c in cpl])
ax[1].set_xticks(range(len(cpl))); ax[1].set_xticklabels([c[0] for c in cpl], rotation=15, ha="right", fontsize=8)
for i, c in enumerate(cpl):
    ax[1].text(i, c[1] + 0.05, f"{c[1]:.2f}", ha="center", fontsize=9)
ax[1].set_ylabel("CoV de uma linha"); ax[1].set_title("(b) o que a direção melhora")
fig.tight_layout(); fig.savefig("figs_ls/ls_coluna_cho.pdf")
print([(c[0], round(c[1], 3)) for c in cpl])
