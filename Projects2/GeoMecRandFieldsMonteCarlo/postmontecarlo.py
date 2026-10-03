# mc_post.py
# uso:  python3 mc_post.py --csv out.csv
import argparse
import io
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

def read_fs_csv(path: Path, include_all: bool = False) -> pd.Series:
    """
    Lê um CSV com 3 colunas (amostra, FS, status), 2 colunas (amostra, FS) OU apenas 1 coluna (FS).
    Aceita com/sem cabeçalho e arquivos mistos (cabeçalho antigo "sample,FS" seguido de linhas com status).
    O main.cpp grava agora TODAS as amostras com status (ok, lo_fail, hi_ok); os arquivos antigos só tinham
    as amostras com FS <= 10. Se houver coluna de status, por padrão mantém o critério antigo (FS <= 10);
    include_all=True usa todas as amostras. Retorna uma Series de FS (float).
    """
    # última linha sem '\n' = gravação em andamento/interrompida (FS truncado): ignorada
    # (o main.cpp a retira do CSV ao retomar)
    text = path.read_text()
    if text and not text.endswith("\n"):
        cut = text.rfind("\n") + 1
        print(f"ignorada linha final incompleta: {text[cut:]!r}")
        text = text[:cut]
    if not text.strip():
        raise ValueError(f"{path} está vazio.")

    # lê como texto com 3 colunas fixas: linhas com menos campos ficam com NaN (formato antigo/misto)
    df = pd.read_csv(io.StringIO(text), header=None, names=[0, 1, 2], dtype=str, skip_blank_lines=True)
    df = df.dropna(how="all").reset_index(drop=True)
    if df.empty:
        raise ValueError(f"{path} está vazio.")

    # cabeçalho: primeira linha sem nenhum valor numérico
    names = {}
    if pd.to_numeric(df.iloc[0], errors="coerce").isna().all():
        names = {str(v).strip().lower(): k for k, v in df.iloc[0].items() if isinstance(v, str)}
        df = df.iloc[1:].reset_index(drop=True)

    # coluna de FS
    col = None
    for name in ("fs", "f_s", "fatorseguranca", "factor_of_safety"):
        if name in names:
            col = names[name]
            break
    if col is None:
        col = 1 if df[1].notna().any() else 0

    # coluna da amostra: a retomada do main.cpp calcula as lacunas depois (linhas fora de ordem);
    # ordena pelo índice da amostra e mantém só a 1ª ocorrência de cada uma
    scol = names.get("sample", names.get("s", 0 if col == 1 else None))
    if scol is not None and scol != col:
        sid = pd.to_numeric(df[scol], errors="coerce")
        dup = sid.notna() & sid.duplicated(keep="first")
        if dup.any():
            print(f"ignoradas {int(dup.sum())} linhas com amostra repetida (mantida a primeira)")
        sid = sid[~dup]
        df = df.loc[sid.sort_values(kind="stable").index].reset_index(drop=True)

    fs = pd.to_numeric(df[col], errors="coerce")
    stcol = names.get("status", 2)
    has_status = col != stcol and df[stcol].notna().any()
    if has_status:
        status = df[stcol].fillna("ok").str.strip()  # linhas antigas (sem status) = ok
        print("amostras por status:", status[fs.notna()].value_counts().to_dict())
        if not include_all:
            nexcl = int((fs > 10.0).sum())
            if nexcl:
                print(f"excluídas {nexcl} amostras com FS > 10 (critério antigo; use --all para incluir)")
            fs = fs.where(fs <= 10.0)

    fs = fs.dropna().reset_index(drop=True)
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
                    help="arquivo CSV com amostras (colunas: s,FS[,status]) ou apenas FS")
    ap.add_argument("--bins", type=int, default=30, help="número de bins do histograma")
    ap.add_argument("--all", action="store_true",
                    help="usa todas as amostras (por padrão, com coluna de status, só FS <= 10 como antes)")
    args = ap.parse_args()

    path = Path(args.csv)
    fs = read_fs_csv(path, include_all=args.all)
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

