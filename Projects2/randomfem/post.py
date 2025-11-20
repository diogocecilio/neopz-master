#!/usr/bin/env python3
# post_pdf_clean.py — análise dos resultados mc_results_*.csv
import glob, re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

# procurar arquivos
files = sorted(glob.glob("mc_results_*.csv"))
if not files:
    print("Nenhum arquivo mc_results_*.csv encontrado!")
    exit(0)

pat = re.compile(r"mc_results_([0-9.]+)_([0-9.]+)\.csv")

stats = []
data_map = {}

for f in files:
    m = pat.search(f)
    if not m:
        continue
    Lx, Ly = float(m.group(1)), float(m.group(2))

    try:
        df = pd.read_csv(f)
    except pd.errors.EmptyDataError:
        print(f"[WARN] arquivo vazio ignorado: {f}")
        continue
    except pd.errors.ParserError:
        try:
            df = pd.read_csv(f, sep=";")
        except Exception as e:
            print(f"[WARN] não consegui ler {f}: {e}")
            continue

    if df.empty or "uy" not in df.columns:
        print(f"[WARN] sem coluna 'uy' em {f}, ignorando.")
        continue

    uy = df["uy"].values
    mean = np.mean(uy)
    std  = np.std(uy, ddof=1)
    q25, q50, q75 = np.percentile(uy, [25,50,75])
    stats.append((Lx, Ly, mean, std, q25, q50, q75, uy.min(), uy.max()))
    data_map[(Lx,Ly)] = uy

stats_df = pd.DataFrame(stats,
    columns=["Lx","Ly","mean","std","q25","median","q75","min","max"])
print("\nResumo estatístico:")
print(stats_df)

# --------------------------------------------------
# Subfiguras: PDF + histograma + caixa de estatísticas
n_cases = len(data_map)
ncols = 4
nrows = int(np.ceil(n_cases / ncols))

fig, axes = plt.subplots(nrows, ncols, figsize=(4*ncols, 3*nrows))
axes = axes.ravel()

for idx, ((Lx,Ly), uy) in enumerate(sorted(data_map.items())):
    ax = axes[idx]

    # histograma normalizado (densidade)
    ax.hist(uy, bins=20, density=True, alpha=0.4,
            color="gray", edgecolor="black")

    # PDF via KDE
    kde = gaussian_kde(uy)
    xgrid = np.linspace(uy.min(), uy.max(), 200)
    pdf = kde(xgrid)
    ax.plot(xgrid, pdf, color="steelblue", lw=2)

    # Estatísticas
    mean = uy.mean()
    std  = uy.std(ddof=1)
    q25, q50, q75 = np.percentile(uy, [25,50,75])
    umin, umax = uy.min(), uy.max()

    # Texto formatado em várias linhas
    textstr = (
        f"Lx={Lx:.2f}, Ly={Ly:.2f}\n"
        f"μ={mean:.3f}, σ={std:.3f}\n"
        f"Q25={q25:.3f}, Med={q50:.3f}, Q75={q75:.3f}\n"
        f"min={umin:.3f}, max={umax:.3f}"
    )

    ax.text(0.98, 0.95, textstr, transform=ax.transAxes,
            fontsize=8, va="top", ha="right",
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.7))

    ax.set_xlabel("uy")
    ax.set_ylabel("PDF")
    # remove título "Caso X"
    # ax.set_title(...)  # não usado
    ax.legend(fontsize=7)

# remover subplots extras
for j in range(idx+1, len(axes)):
    fig.delaxes(axes[j])

fig.suptitle("Distribuições de uy por caso", fontsize=16)

# --------------------------------------------------
# Heatmap da média
pivot = stats_df.pivot(index="Ly", columns="Lx", values="mean")
fig2, ax2 = plt.subplots(figsize=(6,5))
c = ax2.imshow(pivot.values, origin="lower", cmap="viridis",
               extent=[pivot.columns.min(), pivot.columns.max(),
                       pivot.index.min(), pivot.index.max()],
               aspect="auto")
fig2.colorbar(c, ax=ax2, label="mean(uy)")
ax2.set_xlabel("Lx")
ax2.set_ylabel("Ly")
ax2.set_title("Média de uy em cada (Lx,Ly)")

plt.show()

