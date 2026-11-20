#!/usr/bin/env python3
"""resumo.py: a tabela analise limite (LA) x MEF dos casos do Cho (Tabelas 4 e 5 do relatorio_comparacao.pdf).

    python3 -I resumo.py <pasta das saidas>        (no fermi: $DATA/cho_la; no PC: resultados/ ou a copia baixada)

Le o que houver em <pasta>:
  direto/r05.csv, direto/r07.csv                    MC direto LA (campo analitico, semente 2010)
  ismc/<r>_<centro>_s<semente>.csv                   ISMC LA (Nystrom; as sementes de um centro sao somadas)
  ls/<r>_<direcao>/ls.json                           line sampling LA
  form/<r>_<M>/form.json                             FORM LA
e imprime ao lado as referencias do MEF (SRM, tri2 + 4 niveis) das Tabelas 4 e 5 (e, com direto/r00.csv, o r = 0 contra
os 5,50 % da analise limite do Ceron)."""
from __future__ import annotations

import glob
import json
import math
import os
import sys

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)
import mc_la  # noqa: E402

# Tabelas 4 e 5 (MEF): P_f (%), CoV (%)
MEF = {"r05": {"direto": (1.700, 4.9), "ismc_M20": (1.740, 2.2), "ismc_M100": (1.745, 2.2), "ismc_gim20": (1.718, 2.2),
               "ls_M20": (1.739, 5.0), "ls_gim2000": (1.762, 4.2)},
       "r07": {"direto": (0.366, 5.0), "ismc_M20": (0.366, 2.2), "ismc_M100": (0.362, 2.2), "ismc_gim20": (0.352, 2.2),
               "ls_M20": (0.364, 6.8), "ls_gim2000": (0.367, 5.6)},
       # r = 0: sem MEF; a referencia e a analise limite do Ceron (2025, Sec. 5.3.2: 5,50 % com S = 50 000, CoV ~1,9 %;
       # o r nao esta escrito no artigo, r = 0 e inferido pela Fig. 9 do Cho; o campo dele e a KL por EF da secao)
       "r00": {"direto": (5.50, 1.9)}}


def linha(metodo, ref, e, extra=""):
    pf = f"{100 * e['Pf']:.3f} %" if e and e.get("Pf") == e.get("Pf") else "--"
    cv = f"{100 * e['CoV']:.1f} %" if e and math.isfinite(e.get("CoV", float("inf"))) else "--"
    n = str(e.get("N", "")) if e else ""
    mr = f"{ref[0]:.3f} % ({ref[1]:.1f} %)" if ref else "--"
    razao = f"{100 * e['Pf'] / ref[0]:.2f}" if e and ref and e.get("Pf", 0) > 0 else "--"
    print(f"  {metodo:34s} {pf:>11s} {cv:>8s} {n:>8s} {mr:>18s} {razao:>8s}  {extra}")


def main(pasta):
    for r in ("r05", "r07", "r00"):
        if r == "r00":
            if not os.path.exists(os.path.join(pasta, "direto", "r00.csv")):
                continue
            print("\nCho, r = 0: analise limite (LA) x Ceron (2025) LA (a coluna 'MEF' e o Ceron)")
        else:
            print(f"\nCho, r = {'-0,5' if r == 'r05' else '-0,7'}: analise limite (LA) x MEF")
        print(f"  {'metodo':34s} {'P_f LA':>11s} {'CoV':>8s} {'N':>8s} {'P_f MEF (CoV)':>18s} {'LA/MEF':>8s}")
        arq = os.path.join(pasta, "direto", f"{r}.csv")
        if os.path.exists(arq):
            L = mc_la.ler([arq])
            e = mc_la.estat(L)
            G = [float(l["Gamma"]) for l in L if math.isfinite(float(l["Gamma"]))]
            pe = sum(l["modo"] == "pe" for l in L) / max(1, len(L))
            linha("MC direto (analitico, sem. 2010)", MEF[r]["direto"], e,
                  f"Gamma medio {sum(G) / len(G):.3f}; pelo pe {100 * pe:.0f} %")
        for cen in sorted({os.path.basename(f).rsplit("_s", 1)[0] for f in glob.glob(os.path.join(pasta, "ismc", f"{r}_*_s*.csv"))}):
            arqs = sorted(glob.glob(os.path.join(pasta, "ismc", f"{cen}_s*.csv")))
            e = mc_la.estat(mc_la.ler(arqs))
            nome = cen[len(r) + 1:]
            ref = MEF[r].get("ismc_" + nome.replace("mef_", ""))
            linha(f"ISMC centro {nome} ({len(arqs)} sem.)", ref, e)
        for f in sorted(glob.glob(os.path.join(pasta, "ls", f"{r}_*", "ls.json"))):
            d = json.load(open(f))
            nome = os.path.basename(os.path.dirname(f))[len(r) + 1:]
            ref = MEF[r].get("ls_" + nome.replace("mef_", ""))
            linha(f"LS direcao {nome}", ref, {"Pf": d["P_f"], "CoV": d["CoV"], "N": d["linhas"]},
                  f"beta {d['beta']:.3f}; {d['analises_linha']:.1f} analises/linha; perdidas {d.get('perdidas', 0)}")
        for f in sorted(glob.glob(os.path.join(pasta, "form", f"{r}_*", "form.json"))):
            d = json.load(open(f))
            print(f"  FORM LA {os.path.basename(os.path.dirname(f)):10s}: beta {d['beta']:.4f}, P_f FORM "
                  f"{100 * d['P_f']:.3f} %, convergiu {d['convergiu']}, {d['avaliacoes']} avaliacoes")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else os.path.join(AQUI, "resultados"))
