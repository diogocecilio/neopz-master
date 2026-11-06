#!/usr/bin/env python3
# plot_v6.py <pasta> [--saida <pasta>]: figuras e tabela de comparacao do V6 a partir das execucoes do executavel
# (uma subpasta por execucao, com resumo.csv, convergencia.csv e correlacao.csv; e tabelas.csv numa subpasta
# "tabelas"). Gera v6_convergencia.png, v6_correlacao.png e v6_tabelas.png e imprime a tabela "artigo x nosso".
import argparse
import csv
import math
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


def le_csv(caminho):
    with open(caminho, newline="") as f:
        linhas = [l for l in f if not l.startswith("#")]
    return list(csv.DictReader(linhas))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("pasta")
    ap.add_argument("--saida", default=".")
    ap.add_argument("--correlacao", default="b40_h10", help="execucao usada na figura de correlacao")
    a = ap.parse_args()
    os.makedirs(a.saida, exist_ok=True)
    execs = {}
    for d in sorted(os.listdir(a.pasta)):
        p = os.path.join(a.pasta, d)
        if os.path.isfile(os.path.join(p, "resumo.csv")) and os.path.isfile(os.path.join(p, "convergencia.csv")):
            execs[d] = {"resumo": le_csv(os.path.join(p, "resumo.csv"))[0],
                        "conv": le_csv(os.path.join(p, "convergencia.csv")),
                        "corr": le_csv(os.path.join(p, "correlacao.csv"))}
    # tabela: Monte Carlo contra a Eq. 14
    print("execucao        M        fracao  Gamma(Eq14) Gamma_KL  KL/Eq14   N       dp_MC     ep      MC/Eq14  (MC-KL)/ep")
    for nome, e in execs.items():
        r = e["resumo"]
        G, GKL, sd, ep = (float(r[k]) for k in ("Gamma_eq14", "Gamma_KL", "dp_MC", "ep_dp_MC"))
        N = int(r["N"])
        z = (sd - GKL) / ep if ep > 0 else float("nan")
        print(f"{nome:14s} {int(r['M']):8d} {float(r['fracao']):7.4f} {G:9.5f} {GKL:9.5f} {GKL / G:8.5f} {N:7d} "
              f"{sd:9.5f} {ep:8.5f} {sd / G if N else float('nan'):8.5f} {z:7.2f}")
    # figura 1: convergencia com M
    if execs:
        fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
        for nome, e in execs.items():
            if e["resumo"].get("convencao", "eq13") != "eq13":
                continue
            M = [int(c["M"]) for c in e["conv"]]
            ax[0].semilogx(M, [float(c["Gamma_KL_sobre_Gamma"]) for c in e["conv"]], "o-", ms=3, label=nome)
            ax[1].semilogx(M, [float(c["var_P0"]) for c in e["conv"]], "o-", ms=3, label=nome)
        ax[0].axhline(1., color="k", lw=0.8)
        ax[0].set_xlabel("termos da KL, M")
        ax[0].set_ylabel(r"$\Gamma_{KL}(M) / \Gamma_{Eq.14}$ (desvio da media no arco)")
        ax[1].axhline(1., color="k", lw=0.8)
        ax[1].set_xlabel("termos da KL, M")
        ax[1].set_ylabel("variancia no meio do arco (campo truncado)")
        for x in ax:
            x.grid(True, which="both", alpha=0.3)
        ax[0].legend(fontsize=7, ncol=2)
        fig.tight_layout()
        fig.savefig(os.path.join(a.saida, "v6_convergencia.png"), dpi=130)
    # figura 2: correlacao nas defasagens
    if a.correlacao in execs:
        e = execs[a.correlacao]
        Lx, Ly = float(e["resumo"]["Lx_m"]), float(e["resumo"]["Ly_m"])
        fig, ax = plt.subplots(figsize=(6.4, 4.4))
        for dirc, cor, L in (("dx", "C0", Lx), ("dy", "C1", Ly)):
            pts = [c for c in e["corr"] if c["par"].startswith(dirc)]
            t = [abs(float(c["dx_m" if dirc == "dx" else "dy_m"])) / L for c in pts]
            ax.plot(t, [float(c["rho_KL"]) for c in pts], "s", mfc="none", color=cor, label=f"KL truncada exata ({dirc})")
            if pts and pts[0]["rho_MC"]:
                ax.errorbar(t, [float(c["rho_MC"]) for c in pts], yerr=[2 * float(c["ep_rho_MC"]) for c in pts],
                            fmt="x", color=cor, label=f"Monte Carlo ({dirc}, 2 EP)")
        tt = [i * 0.05 for i in range(0, 71)]
        ax.plot(tt, [math.exp(-x) for x in tt], "k-", lw=1, label=r"nucleo $e^{-|\tau|/L}$")
        ax.set_xlabel(r"defasagem $|\tau| / L$ ($L_x$ = %.2g m, $L_y$ = %.2g m)" % (Lx, Ly))
        ax.set_ylabel("correlacao")
        ax.set_title(f"{a.correlacao}: M = {e['resumo']['M']}, fracao {float(e['resumo']['fracao']):.3f}, "
                     f"N = {e['resumo']['N']}", fontsize=9)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(os.path.join(a.saida, "v6_correlacao.png"), dpi=130)
    # figura 3: as Tabelas 5 e 6
    tab = os.path.join(a.pasta, "tabelas", "tabelas.csv")
    if os.path.isfile(tab):
        T = le_csv(tab)
        print("\ncaso      arco    T(m)   FS_B   FS_art  s_cL   art    dif     s_phiL  art    dif    L=delta: s_cL")
        for r in T:
            sc, sca, sp, spa = (float(r[k]) for k in ("sigma_cL", "sigma_cL_artigo", "sigma_phiL", "sigma_phiL_artigo"))
            print(f"{r['caso']:9s} {r['arco']:6s} {float(r['T_m']):6.1f} {float(r['FS_bishop']):6.3f} "
                  f"{float(r['FS_artigo']):6.3f} {sc:6.3f} {sca:6.3f} {100 * (sc / sca - 1):+6.1f}% {sp:6.3f} {spa:6.3f} "
                  f"{100 * (sp / spa - 1):+6.1f}%  {float(r['Gamma_L_igual_delta']) * 8.68:6.3f}")
        fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
        nomes = [r["caso"] for r in T]
        x = range(len(T))
        for i, (k, ka, s) in enumerate((("sigma_cL", "sigma_cL_artigo", "kPa"), ("sigma_phiL", "sigma_phiL_artigo", "graus"))):
            ax[i].plot(x, [float(r[ka]) for r in T], "ks", label="artigo (Tabelas 5 e 6)")
            ax[i].plot(x, [float(r[k]) for r in T], "o", color="C0", label="Eq. 14, L = delta/2 (Eq. 13)")
            sig = 8.68 if i == 0 else 2.99
            ax[i].plot(x, [float(r["Gamma_L_igual_delta"]) * sig for r in T], "^", color="C1", mfc="none",
                       label="Eq. 14, L = delta")
            ax[i].set_xticks(list(x))
            ax[i].set_xticklabels(nomes, rotation=45, fontsize=8)
            ax[i].set_ylabel(f"{k} ({s})")
            ax[i].grid(True, alpha=0.3)
        ax[0].legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(os.path.join(a.saida, "v6_tabelas.png"), dpi=130)
    return 0


if __name__ == "__main__":
    sys.exit(main())
