"""Graficos do LearningSlopeWater: FS (SRM) e poropressao nos pontos de controle ao longo do tempo.

Uso (na pasta onde o programa gravou os CSVs):
    python3 plot_fs_pressao.py [fs_tempo.csv] [pressao_tempo.csv] [saida.png]
"""
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

fs_file = sys.argv[1] if len(sys.argv) > 1 else "fs_tempo.csv"
p_file = sys.argv[2] if len(sys.argv) > 2 else "pressao_tempo.csv"
out = sys.argv[3] if len(sys.argv) > 3 else "fs_pressao.png"

fs = pd.read_csv(fs_file)
pw = pd.read_csv(p_file)
estilo = {"acoplado": "-", "desacoplado": "--"}

fig, ax = plt.subplots(1, 2, figsize=(13, 5))
for caso, d in fs.groupby("caso"):
    t = d.t_h.clip(lower=1e-3)  # t = 0 no eixo logaritmico
    ax[0].semilogx(t, d.FS, estilo[caso], marker="o", label=caso)
ax[0].set_xlabel("tempo (h)  [t = 0 desenhado em 1e-3 h]")
ax[0].set_ylabel("FS (SRM)")
ax[0].set_title("Fator de seguranca ao longo do tempo")
ax[0].grid(alpha=0.3, which="both")
ax[0].legend()

cores = {}
for (caso, ponto), d in pw.groupby(["caso", "ponto"]):
    cor = cores.setdefault(ponto, f"C{len(cores)}")
    t = d.t_h.clip(lower=1e-3)
    ax[1].semilogx(t, d.pw_kPa, estilo[caso], color=cor,
                   label=f"{ponto} ({d.x.iloc[0]:.0f}, {d.y.iloc[0]:.0f})" if caso == "acoplado" else None)
ax[1].plot([], [], "k-", label="acoplado")
ax[1].plot([], [], "k--", label="desacoplado")
ax[1].set_xlabel("tempo (h)")
ax[1].set_ylabel("p_w (kPa)")
ax[1].set_title("Poropressao nos pontos de controle")
ax[1].grid(alpha=0.3, which="both")
ax[1].legend(fontsize=8)
fig.tight_layout()
fig.savefig(out, dpi=120)
print("grafico:", out)
