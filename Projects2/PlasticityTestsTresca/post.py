import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

# -------------------------------------------------------------------
# Caminho do CSV
# -------------------------------------------------------------------
#csv_path = "/home/diogo/projects/neopz-master-build-release/Projects2/PlasticityTestsTresca/loadsweep.csv"
csv_path = "/home/diogo/projects/neopz-master-build-release/Projects2/PlasticityTestsTresca/loadsweepmises.csv"
#loadsweep.csv
#loadsweeptresca.csv
#loadsweepmises.csv
print("Working dir  :", os.getcwd())
print("CSV filename :", csv_path)
print("Existe CSV?  :", os.path.exists(csv_path))

df = pd.read_csv(csv_path)

print("\nColunas encontradas:", list(df.columns))
print("\nPrimeiras linhas do CSV:")
print(df.head())
print("\nShape (linhas, colunas):", df.shape)

u           = df["u"].values
intStressYY = df["intStressYY"].values

# -------------------------------------------------------------------
# Configurações de figura (estilo “livro”, tudo em preto)
# -------------------------------------------------------------------
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype']  = 42
plt.rcParams['font.size']    = 10

# --- parâmetros físicos ---
B  = 1.0      # largura da sapata (ajuste se precisar)
c0 = 490.0    # coesão em kPa

# integral * 2 = pressão no topo (intStressYY < 0 → compressão)
pressure = -2.0 * intStressYY        # P > 0

# recalque positivo e normalizado
u_norm = np.abs(u) / B               # u/B >= 0

# pressão normalizada
P_norm = pressure / c0               # P/c

slipline_limit = 5.14

# aspect ratio parecido com a figura do livro
fig, ax = plt.subplots(figsize=(3.6, 2.7))   # largura x altura (polegadas)

# -------------------------------------------------------------------
# Curva EF em preto
# -------------------------------------------------------------------
ax.plot(
    u_norm, P_norm,
    linestyle='-',
    marker='D',
    markersize=4,
    linewidth=1.0,
    color='black',
    markerfacecolor='white',
    markeredgecolor='black',
    label="finite element results"
)

# Linha limite (Prandtl) em preto tracejado
ax.axhline(
    slipline_limit,
    linestyle='--',
    linewidth=1.0,
    color='black',
    label=f"Prandtl solution: {slipline_limit:.1f}"
)

# Rótulos dos eixos
ax.set_xlabel(r"normalised settlement, $u/B$")
ax.set_ylabel(r"normalised pressure, $P/c$")

# Grade leve em cinza
ax.grid(True, linestyle=':', linewidth=0.5, color='0.7', alpha=0.9)

# Limites (ajuste fino se quiser igual ao livro)
#ax.set_xlim(0.0, 0.003)   # para ficar igual à figura de referência
#ax.set_ylim(0.0, 17.0)

# Eixos em preto
for spine in ax.spines.values():
    spine.set_color('black')

ax.tick_params(axis='both', colors='black')

# Legenda sem moldura, posição semelhante
ax.legend(loc="lower right", frameon=False)

fig.tight_layout()
plt.savefig("u_vs_Pnorm_tresca.pdf")
plt.savefig("u_vs_Pnorm_tresca.png", dpi=300)
plt.show()
