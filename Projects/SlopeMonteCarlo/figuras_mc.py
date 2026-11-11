#!/usr/bin/env python3
"""Figuras do Monte Carlo (os mc_<config>_*.csv de uma pasta), como a Fig. 25 de Vargas Ceron et al. (2025):
  (a) histograma de Gamma com o ajuste lognormal deste Monte Carlo e o lognormal do artigo (media 1,353, CoV 23,54 %);
  (b) P_f = n_f / S nas S primeiras realizacoes (em ordem de numero), com a faixa de +-2 desvios (Eq. 68) e os 11,50 %
      do artigo.
Uso: figuras_mc.py [pasta = .] [config = u0h3p2]   ->   <pasta>/mc_<config>.png
"""
import glob
import math
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

pasta = sys.argv[1] if len(sys.argv) > 1 else "."
config = sys.argv[2] if len(sys.argv) > 2 else "u0h3p2"

gamma = {}
for arq in glob.glob(os.path.join(pasta, f"mc_{config}_*.csv")):
    dados = np.genfromtxt(arq, delimiter=",", names=True)
    for r, g in zip(np.atleast_1d(dados["real"]), np.atleast_1d(dados["gamma"])):
        gamma[int(r)] = float(g)
det = gamma.pop(-1, None)
reais = sorted(gamma)
g = np.array([gamma[r] for r in reais])
if g.size == 0:
    sys.exit(f"nenhuma realizacao em {pasta}/mc_{config}_*.csv")

AZUL, LARANJA, TINTA, TINTA2 = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e"
plt.rcParams.update({"font.size": 10, "axes.edgecolor": TINTA2, "axes.labelcolor": TINTA, "xtick.color": TINTA2,
                     "ytick.color": TINTA2, "axes.spines.top": False, "axes.spines.right": False})
fig, (a, b) = plt.subplots(1, 2, figsize=(11, 4.2), layout="constrained")

# (a) histograma e lognormais
pos = g[g > 0]
mu, s = np.mean(np.log(pos)), np.std(np.log(pos))  # maxima verossimilhanca
cov_c = 0.2354
s_c = math.sqrt(math.log(1 + cov_c**2))
mu_c = math.log(1.353) - 0.5 * s_c**2
x = np.linspace(0.3, 3.0, 400)  # Gamma > 3 sao raros: ficam fora do grafico, mas entram nas contas
a.hist(g[g <= 3.0], bins=np.arange(0.3, 3.0001, 0.05), density=True, color="#b7d3f6", edgecolor="#fcfcfb", linewidth=0.8, label=f"este ({g.size} realizacoes)")
a.plot(x, stats.lognorm.pdf(x, s, scale=math.exp(mu)), color=AZUL, lw=2,
       label=f"lognormal deste MC: media {math.exp(mu + s * s / 2):.3f}, CoV {100 * math.sqrt(math.exp(s * s) - 1):.1f}%")
a.plot(x, stats.lognorm.pdf(x, s_c, scale=math.exp(mu_c)), color=LARANJA, lw=2, ls="--",
       label="lognormal do artigo: media 1,353, CoV 23,5%")
a.axvline(1.0, ymax=0.7, color=TINTA2, lw=1, ls=":")
if det is not None:
    a.axvline(det, ymax=0.7, color=TINTA, lw=1, label=f"deterministico (medias): {det:.3f}")
a.set_xlabel("fator de estabilidade $\\Gamma$")
a.set_ylabel("densidade de probabilidade")
a.set_title("(a) distribuicao de $\\Gamma$", loc="left", color=TINTA)
a.set_xlim(0.3, 3.0)
a.set_ylim(0, a.get_ylim()[1] * 1.35)  # espaco para a legenda acima das barras
a.legend(frameon=False, fontsize=8, loc="upper right")
a.grid(axis="y", color="#e6e5e1", lw=0.6)

# (b) convergencia de P_f
S = np.arange(1, g.size + 1)
pf = np.cumsum(g < 1.0) / S
desvio = np.sqrt(np.clip(pf * (1 - pf), 0, None) / S)
b.fill_between(S, 100 * (pf - 2 * desvio), 100 * (pf + 2 * desvio), color="#cde2fb", lw=0, label="$\\pm 2$ desvios (Eq. 68)")
b.plot(S, 100 * pf, color=AZUL, lw=2, label=f"este: $P_f$ = {100 * pf[-1]:.2f}% com S = {g.size}")
b.axhline(11.50, color=LARANJA, lw=2, ls="--", label="artigo: $P_f$ = 11,50% (S = 10000)")
b.set_xlabel("numero de realizacoes S")
b.set_ylabel("$P_f$ (%)")
b.set_ylim(0, max(25, 100 * pf[min(len(pf) - 1, 20):].max() * 1.2))
b.set_title("(b) convergencia da probabilidade de ruptura", loc="left", color=TINTA)
b.legend(frameon=False, fontsize=8.5, loc="upper right")
b.grid(axis="y", color="#e6e5e1", lw=0.6)

saida = os.path.join(pasta, f"mc_{config}.png")
fig.savefig(saida, dpi=150)
print(f"{saida}: S = {g.size}, P_f = {100 * pf[-1]:.2f}%, lognormal media {math.exp(mu + s * s / 2):.4f}")
