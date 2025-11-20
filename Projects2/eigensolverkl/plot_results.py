# plot_results.py
# uso:  python plot_results.py results.csv
# cria os arquivos:
#   err_vs_ndof.png, err_vs_nel.png, time_vs_ndof.png, err_vs_h.png

import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

REQ_COLS = ["ref","nel","ndof","asm_s","solve_s","L2","Energy","eig0"]

def fit_loglog_slope(x, y):
    """retorna o coeficiente angular da regressão log(y) x log(x)."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = (x > 0) & (y > 0)
    if mask.sum() < 2:
        return np.nan
    p = np.polyfit(np.log(x[mask]), np.log(y[mask]), 1)
    return p[0]

def main(csv_path: Path):
    df = pd.read_csv(csv_path)
    # checagem de colunas
    missing = [c for c in REQ_COLS if c not in df.columns]
    if missing:
        raise ValueError(f"CSV faltando colunas: {missing}")

    # ordena por ref e cria passo de malha h (malha inicial 3x3, h0 = 1/3)
    df = df.sort_values("ref").reset_index(drop=True)
    df["h"] = (1.0/3.0) * (0.5 ** df["ref"])

    print("\n=== Dados lidos ===")
    print(df.to_string(index=False))

    # —— Gráfico 1: erros vs ndof (log-log) ——
    plt.figure()
    plt.loglog(df["ndof"], df["L2"], marker="o", linestyle="-", label=r"$\|e\|_{0}$")
    plt.loglog(df["ndof"], df["Energy"], marker="s", linestyle="-", label=r"$\|e\|_{E}$")
    plt.xlabel("NDOF")
    plt.ylabel("Erro")
    plt.title("Erro vs NDOF (log–log)")
    plt.grid(True, which="both", ls=":")
    plt.legend()
    plt.tight_layout()
    plt.savefig("err_vs_ndof.png", dpi=200)

    # —— Gráfico 2: erros vs nel (log-log) ——
    plt.figure()
    plt.loglog(df["nel"], df["L2"], marker="o", linestyle="-", label=r"$\|e\|_{0}$")
    plt.loglog(df["nel"], df["Energy"], marker="s", linestyle="-", label=r"$\|e\|_{E}$")
    plt.xlabel("Número de elementos (nel)")
    plt.ylabel("Erro")
    plt.title("Erro vs nel (log–log)")
    plt.grid(True, which="both", ls=":")
    plt.legend()
    plt.tight_layout()
    plt.savefig("err_vs_nel.png", dpi=200)

    # —— Gráfico 3: tempos vs ndof ——
    plt.figure()
    plt.plot(df["ndof"], df["asm_s"], marker="o", linestyle="-", label="Montagem (s)")
    plt.plot(df["ndof"], df["solve_s"], marker="s", linestyle="-", label="Solve (s)")
    plt.xlabel("NDOF")
    plt.ylabel("Tempo [s]")
    plt.title("Tempos vs NDOF")
    plt.grid(True, ls=":")
    plt.legend()
    plt.tight_layout()
    plt.savefig("time_vs_ndof.png", dpi=200)

    # —— Gráfico 4: erro vs h (log-log) + ordem aproximada ——
    slope_L2     = fit_loglog_slope(df["h"], df["L2"])
    slope_Energy = fit_loglog_slope(df["h"], df["Energy"])
    print(f"\nOrdem (~inclinação log-log)  L2  ≈ {slope_L2:.3f}")
    print(f"Ordem (~inclinação log-log)  E   ≈ {slope_Energy:.3f}")

    plt.figure()
    plt.loglog(df["h"], df["L2"], marker="o", linestyle="-", label=rf"$\|e\|_0$  (slope≈{slope_L2:.2f})")
    plt.loglog(df["h"], df["Energy"], marker="s", linestyle="-", label=rf"$\|e\|_E$  (slope≈{slope_Energy:.2f})")
    plt.xlabel("h (passo de malha ~ 1/3·2^{-ref})")
    plt.ylabel("Erro")
    plt.title("Erro vs h (log–log)")
    plt.grid(True, which="both", ls=":")
    plt.legend()
    plt.tight_layout()
    plt.savefig("err_vs_h.png", dpi=200)

    print("\nArquivos gerados:")
    for fn in ["err_vs_ndof.png","err_vs_nel.png","time_vs_ndof.png","err_vs_h.png"]:
        print(" -", fn)

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("uso: python plot_results.py results.csv")
        sys.exit(1)
    main(Path(sys.argv[1]))

