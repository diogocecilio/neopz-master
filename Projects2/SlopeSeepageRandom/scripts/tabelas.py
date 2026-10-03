#!/usr/bin/env python3
"""Tabelas (markdown) dos resultados do SlopeSeepageRandom comparados com o artigo.

Uso: tabelas.py <diretório de resultados>   (com det/, mc/<nome>/<nome>.csv, mc_ref/ref.csv e rebaix/)
"""
import csv
import glob
import math
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from analisa_mc import ARTIGO, estatisticas, le  # noqa: E402

# colunas determinísticas das Tabelas 5 e 6 do artigo
ARTIGO_ALFA = {1: 1.336, 2: 1.533, 3: 1.674, 4: 1.783, 5: 1.872}
ARTIGO_HW = {0.5: 1.671, 0.6: 1.494, 0.7: 1.383, 0.8: 1.322, 0.9: 1.307, 1.0: 1.336}


def det(log):
    """(Γ, FS) de um log do comando det"""
    g = f = None
    if not os.path.exists(log):
        return g, f
    for ln in open(log, errors="replace"):
        m = re.search(r"Gamma \(fator de carga\) = ([0-9.eE+-]+)", ln)
        if m:
            g = float(m.group(1))
        m = re.search(r"FS \(reducao de resistencia\) = ([0-9.eE+-]+)", ln)
        if m:
            f = float(m.group(1))
    return g, f


def fmt(v, d=3):
    return "—" if v is None or (isinstance(v, float) and math.isnan(v)) else f"{v:.{d}f}"


def tabela_det(res):
    print("### Varreduras determinísticas (Mohr-Coulomb, `adapt=2`)\n")
    rows = []
    for a in (1, 2, 3, 4, 5):
        g, _ = det(f"{res}/det/alfa{a}.log")
        rows.append((a, g, ARTIGO_ALFA[a]))
    if any(r[1] for r in rows):
        print("| α = k_h/k_v | Γ (FE) | Γ artigo (Tab. 5) | FE/artigo |\n|---|---|---|---|")
        for a, g, ref in rows:
            print(f"| {a} | {fmt(g)} | {ref:.3f} | {fmt(g / ref if g else None)} |")
        print()
    rows = []
    for hw in (2.5, 3, 3.5, 4, 4.5, 5):
        g, _ = det(f"{res}/det/hw{hw}.log")
        rows.append((hw / 5., g, ARTIGO_HW[round(hw / 5., 1)]))
    if any(r[1] for r in rows):
        print("| h_w/H | Γ (FE) | Γ artigo (Tab. 6) | FE/artigo |\n|---|---|---|---|")
        for r, g, ref in rows:
            print(f"| {r:.1f} | {fmt(g)} | {ref:.3f} | {fmt(g / ref if g else None)} |")
        print()
    betas = (15, 30, 45, 60, 75, 90)
    tab = {(a, b): det(f"{res}/det/beta{b}_alfa{a}.log")[0] for a in (1, 5, 10) for b in betas}
    if any(tab.values()):
        print("Fig. 9 (Γ × β):\n")
        print("| β (graus) | α = 1 | α = 5 | α = 10 |\n|---|---|---|---|")
        for b in betas:
            print(f"| {b} | " + " | ".join(fmt(tab[(a, b)]) for a in (1, 5, 10)) + " |")
        print()


def tabela_mc(res):
    print("### Monte Carlo (Γ; `h=1 adapt=2`, KL com M modos e compensação da variância)\n")
    print("| caso | N | μ | σ | CoV % | Pf % | CoV(Pf) % | artigo μ | σ | Pf % |\n|---|---|---|---|---|---|---|---|---|---|")
    casos = [("referencia", f"{res}/mc_ref/ref.csv")]
    for d in sorted(glob.glob(f"{res}/mc/*/")):
        nome = os.path.basename(d.rstrip("/"))
        arq = f"{d}{nome}.csv"
        if not os.path.exists(arq):  # ainda em execução: junta as partes
            partes = sorted(glob.glob(f"{d}{nome}_parte*.csv"))
            if not partes:
                continue
            arq = partes[0] if len(partes) == 1 else None
            if arq is None:
                continue
        casos.append((nome, arq))
    for nome, arq in casos:
        if not os.path.exists(arq):
            continue
        v, outros = le(arq)
        if len(v) < 2:
            continue
        n, mu, sd, pf, covpf = estatisticas(v)
        ref = ARTIGO.get(nome, (None, None, None))
        print(f"| {nome} | {n} | {mu:.3f} | {sd:.3f} | {100 * sd / mu:.1f} | {100 * pf:.2f} | "
              f"{fmt(100 * covpf, 1) if pf > 0 else '—'} | {fmt(ref[0])} | {fmt(ref[1])} | {fmt(ref[2], 2)} |")
    print()


def tabela_rebaixamento(res):
    arqs = sorted(glob.glob(f"{res}/rebaix/*.csv"))
    if not arqs:
        return
    print("### Rebaixamento acoplado: FS e Γ com a poropressão congelada\n")
    for arq in arqs:
        linhas = list(csv.DictReader(open(arq)))
        if not linhas:
            continue
        nome = os.path.basename(arq)[:-4]
        est = [r for r in linhas if r["T"] == "inf"]
        print(f"**{nome}** (T_d = {linhas[0]['Td']}; estacionário desacoplado: FS = {fmt(float(est[0]['FS']))}, "
              f"Γ = {fmt(float(est[0]['Gamma'])) if est[0]['Gamma_status'] else '—'})\n" if est else f"**{nome}**\n")
        print("| T = c_v t/H² | z_w (m) | u_A (kPa) | máx. desloc. (m) | pontos plásticos | FS | Γ |\n|---|---|---|---|---|---|---|")
        for r in linhas:
            if r["T"] == "inf":
                continue
            if "colapso" in r["FS_status"]:
                print(f"| {float(r['T']):.4g} | {float(r['zw']):.3f} | colapso no u-p acoplado | | | | |")
                continue
            print(f"| {float(r['T']):.4g} | {float(r['zw']):.2f} | {float(r['u_A']):.2f} | {float(r['umax_desloc']):.4f} | "
                  f"{r['pontos_plasticos']} | {fmt(float(r['FS']))} | {fmt(float(r['Gamma'])) if r['Gamma_status'] else '—'} |")
        print()


def main():
    res = sys.argv[1] if len(sys.argv) > 1 else "."
    tabela_det(res)
    tabela_mc(res)
    tabela_rebaixamento(res)


if __name__ == "__main__":
    main()
