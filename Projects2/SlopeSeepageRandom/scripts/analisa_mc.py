#!/usr/bin/env python3
"""Estatísticas das simulações de Monte Carlo do SlopeSeepageRandom (comando mc) e comparação com o artigo.

Uso: analisa_mc.py arquivo.csv [rotulo=arquivo.csv ...] [--ref "mu,sigma,Pf%"] [--figura saida.png]

Para cada CSV: N, média, desvio, CoV do fator (Γ ou FS), Pf = P(fator < 1) e CoV(Pf) = sqrt((1 - Pf)/(N Pf)).
O fator de cada amostra é o ponto médio do intervalo [último convergido, primeiro sem convergência] quando o
colapso foi encontrado (status ok) e o último valor convergido nos demais casos (status != ok são contados à
parte). Com --figura, desenha a densidade (KDE gaussiana), a convergência de Pf e a da média.
"""
import argparse
import csv
import math
import sys

# Tabelas 3-6 do artigo (Vargas Ceron et al., 2025): (μ, σ, Pf %) da referência e de algumas variações
ARTIGO = {
    "referencia": (1.353, 0.318, 11.50),
    "covk0": (1.329, 0.288, 11.21),
    "covc10": (1.367, 0.180, 0.45),
    "covc50": (1.324, 0.472, 25.78),
    "covc70": (1.296, 0.620, 36.34),
    "covphi5": (1.351, 0.300, 10.20),
    "covphi20": (1.361, 0.384, 15.18),
    "s1.5": (1.356, 0.352, 13.94),
    "s2": (1.361, 0.372, 14.99),
    "s5": (1.350, 0.404, 18.93),
    "s20": (1.355, 0.430, 20.23),
    "covk100": (1.391, 0.360, 11.44),
    "alfa2": (1.562, 0.390, 3.90),
    "alfa3": (1.704, 0.432, 1.83),
    "alfa4": (1.821, 0.474, 1.00),
    "alfa5": (1.910, 0.504, 0.63),
    "hw0.5": (1.670, 0.378, 1.38),
    "hw0.6": (1.502, 0.345, 4.20),
    "hw0.7": (1.393, 0.314, 8.54),
    "hw0.8": (1.338, 0.309, 11.99),
    "hw0.9": (1.325, 0.306, 12.58),
    # seção 5.3 (Cho 2010): só Pf é dado no texto (Cho: 7.9% e 6.37%)
    "cho_coesivo": (None, None, 6.5),
    "cho_cphi": (None, None, 5.50),
}


def le(arquivo):
    vals, outros = [], {}
    with open(arquivo) as f:
        for r in csv.DictReader(f):
            st = r["status"]
            lo, up = float(r["fator"]), float(r["limite_superior"])
            if st == "ok" and up > lo:
                vals.append(0.5 * (lo + up))
            else:
                # sem colapso até o fator máximo (limite_maximo) ou colapso antes do primeiro passo
                vals.append(lo)
                outros[st] = outros.get(st, 0) + 1
    return vals, outros


def estatisticas(v):
    n = len(v)
    mu = sum(v) / n
    sd = math.sqrt(sum((x - mu) ** 2 for x in v) / (n - 1)) if n > 1 else 0.0
    pf = sum(1 for x in v if x < 1.0) / n
    covpf = math.sqrt((1 - pf) / (n * pf)) if pf > 0 else float("inf")
    return n, mu, sd, pf, covpf


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("arquivos", nargs="+", help="arquivo.csv ou rotulo=arquivo.csv (rotulo de ARTIGO compara)")
    ap.add_argument("--figura", default=None)
    args = ap.parse_args()
    series = []
    print(f"{'caso':>14} {'N':>6} {'media':>7} {'desvio':>7} {'CoV%':>6} {'Pf%':>7} {'CoV(Pf)%':>8}   "
          f"{'artigo: media':>13} {'desvio':>7} {'Pf%':>6}  outros status")
    for a in args.arquivos:
        rot, arq = (a.split("=", 1) if "=" in a else (a, a))
        v, outros = le(arq)
        if not v:
            print(f"{rot:>14} sem amostras")
            continue
        n, mu, sd, pf, covpf = estatisticas(v)
        ref = ARTIGO.get(rot)
        fmt = lambda v, w, d: f"{v:{w}.{d}f}" if v is not None else f"{'-':>{w}}"
        refs = (f"{fmt(ref[0], 13, 3)} {fmt(ref[1], 7, 3)} {fmt(ref[2], 6, 2)}" if ref
                else f"{'-':>13} {'-':>7} {'-':>6}")
        print(f"{rot:>14} {n:6d} {mu:7.3f} {sd:7.3f} {100 * sd / mu:6.1f} {100 * pf:7.2f} {100 * covpf:8.2f}   "
              f"{refs}  {outros if outros else ''}")
        series.append((rot, v))
    if args.figura and series:
        import numpy as np
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
        for rot, v in series:
            x = np.asarray(v)
            grid = np.linspace(max(0.0, x.min() - 0.3), x.max() + 0.3, 400)
            h = 1.06 * x.std() * len(x) ** (-0.2)  # regra de Silverman
            dens = np.exp(-0.5 * ((grid[:, None] - x[None, :]) / h) ** 2).sum(1) / (len(x) * h * math.sqrt(2 * math.pi))
            ax[0].plot(grid, dens, label=rot)
            n = np.arange(1, len(x) + 1)
            ax[1].plot(n, 100 * np.cumsum(x < 1.0) / n, label=rot)
            ax[2].plot(n, np.cumsum(x) / n, label=rot)
            ref = ARTIGO.get(rot)
            if ref:
                ax[1].axhline(ref[2], ls="--", lw=0.8, color=ax[1].lines[-1].get_color())
                if ref[0] is not None:
                    ax[2].axhline(ref[0], ls="--", lw=0.8, color=ax[2].lines[-1].get_color())
        ax[0].axvline(1.0, color="k", lw=0.8)
        ax[0].set_xlabel("fator de estabilidade")
        ax[0].set_ylabel("densidade (KDE)")
        ax[1].set_xlabel("N")
        ax[1].set_ylabel("Pf (%)")
        ax[2].set_xlabel("N")
        ax[2].set_ylabel("média")
        for a_ in ax:
            a_.grid(alpha=0.3)
        ax[0].legend()
        ax[1].set_title("tracejado: artigo")
        fig.tight_layout()
        fig.savefig(args.figura, dpi=130)
        print("figura:", args.figura)


if __name__ == "__main__":
    sys.exit(main())
