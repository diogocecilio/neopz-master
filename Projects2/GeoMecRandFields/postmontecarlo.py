# mc_post.py
# uso:  python3 mc_post.py --csv out.csv
import argparse
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

def read_fs_csv(path: Path) -> pd.Series:
    """
    Lê um CSV com 2 colunas (amostra, FS) OU apenas 1 coluna (FS).
    Aceita com/sem cabeçalho. Retorna uma Series de FS (float).
    """
    try:
        df = pd.read_csv(path)
    except Exception:
        # fallback: separador pode ser espaço/; tente leitura genérica
        df = pd.read_csv(path, header=None)

    # tenta detectar a coluna FS
    col = None
    for name in df.columns:
        if str(name).strip().lower() in ("fs", "f_s", "fatorseguranca", "factor_of_safety"):
            col = name
            break

    if col is None:
        # se há exatamente 2 colunas, assuma a 2ª como FS
        if df.shape[1] == 2:
            col = df.columns[1]
        elif df.shape[1] == 1:
            col = df.columns[0]
        else:
            raise ValueError(
                f"Não consegui identificar a coluna de FS em {path}. "
                f"Colunas lidas: {list(df.columns)}"
            )

    fs = pd.to_numeric(df[col], errors="coerce").dropna().reset_index(drop=True)
    if fs.empty:
        raise ValueError("Coluna de FS ficou vazia após coerção numérica.")
    return fs

def running_mean(x: np.ndarray) -> np.ndarray:
    csum = np.cumsum(x, dtype=float)
    n = np.arange(1, len(x)+1, dtype=float)
    return csum / n

def running_std(x: np.ndarray) -> np.ndarray:
    # variância incremental (Welford) para estabilidade
    mean = 0.0
    m2 = 0.0
    out = np.zeros_like(x, dtype=float)
    for i, val in enumerate(x, start=1):
        delta = val - mean
        mean += delta / i
        m2 += delta * (val - mean)
        if i > 1:
            out[i-1] = math.sqrt(m2 / (i - 1))
        else:
            out[i-1] = 0.0
    return out

def ci95(mean: float, std: float, n: int):
    # IC 95% para média: mean ± t_{0.975, n-1} * std/sqrt(n)
    if n < 2 or not math.isfinite(std):
        return (float("nan"), float("nan"))
    try:
        # t crítico; se scipy não estiver disponível, usa z=1.96
        import scipy.stats as st
        tcrit = st.t.ppf(0.975, df=n-1)
    except Exception:
        tcrit = 1.96
    hw = tcrit * std / math.sqrt(n)
    return (mean - hw, mean + hw)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", type=str, default="mc_results_20.000000_2.000000.csv",
                    help="arquivo CSV com amostras (colunas: s,FS) ou apenas FS")
    ap.add_argument("--bins", type=int, default=30, help="número de bins do histograma")
    args = ap.parse_args()

    path = Path(args.csv)
    fs = read_fs_csv(path)
    x = fs.values.astype(float)
    n = x.size

    # estatísticas
    mean = float(np.mean(x))
    std  = float(np.std(x, ddof=1)) if n > 1 else float("nan")
    mn, mx = float(np.min(x)), float(np.max(x))
    se = std / math.sqrt(n) if n > 0 and math.isfinite(std) else float("nan")
    lo, hi = ci95(mean, std, n)

    # convergence (médias e desvios acumulados)
    rm = running_mean(x)
    rs = running_std(x)
    idx = np.arange(1, n+1)

    # IC 95% acumulado para a média: mean_k ± 1.96 * std_k/sqrt(k)  (usa z≈1.96)
    z = 1.96
    with np.errstate(invalid="ignore", divide="ignore"):
        rm_lo = rm - z * rs / np.sqrt(idx)
        rm_hi = rm + z * rs / np.sqrt(idx)

    # ---- Plot 1: Convergência da média ----
    plt.figure()
    plt.plot(idx, rm, marker="o", linestyle="-", label="Média acumulada")
    plt.fill_between(idx, rm_lo, rm_hi, alpha=0.2, label="IC 95% (aprox.)")
    plt.axhline(mean, linestyle="--", label=f"Média final = {mean:.4g}")
    plt.xlabel("Número de amostras")
    plt.ylabel("FS (média acumulada)")
    plt.title("Convergência da média de FS (Monte Carlo)")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    fig1 = path.with_suffix("").as_posix() + "_convergencia.png"
    plt.savefig(fig1, dpi=150)

    # ---- Plot 2: Histograma ----
    plt.figure()
    plt.hist(x, bins=args.bins, density=True, alpha=0.7)
    # curva normal aproximada
    if n > 1 and math.isfinite(std) and std > 0:
        xs = np.linspace(mn, mx, 256)
        pdf = 1.0/(std*np.sqrt(2*np.pi))*np.exp(-0.5*((xs-mean)/std)**2)
        plt.plot(xs, pdf, linewidth=2, label="Normal(mean,std) aprox.")
    plt.xlabel("FS")
    plt.ylabel("Densidade")
    plt.title("Histograma de FS")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    fig2 = path.with_suffix("").as_posix() + "_hist.png"
    plt.savefig(fig2, dpi=150)

    # salva estatísticas em texto
    stats_txt = (
        f"Arquivo: {path}\n"
        f"N        = {n}\n"
        f"mean     = {mean:.10g}\n"
        f"std(ddof=1) = {std:.10g}\n"
        f"stderr   = {se:.10g}\n"
        f"min/max  = {mn:.10g} / {mx:.10g}\n"
        f"95% CI   = [{lo:.10g}, {hi:.10g}]\n"
    )
    out_stats = path.with_suffix("").as_posix() + "_stats.txt"
    with open(out_stats, "w") as f:
        f.write(stats_txt)

    print(stats_txt)
    print(f"Figuras salvas: {fig1}  |  {fig2}")
    print(f"Resumo salvo em: {out_stats}")

    # Mostra na tela (opcional)
    plt.show()

if __name__ == "__main__":
    main()

