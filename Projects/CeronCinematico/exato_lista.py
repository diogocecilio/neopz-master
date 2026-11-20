#!/usr/bin/env python3
"""exato_lista.py: refaz com a busca EXATA (cinematico.minimizar) as realizacoes de um mc.csv com Gamma numa faixa,
para ver se a busca rapida (rapido.minimizar_rapido) perdeu o minimo perto da falha.

    python3 exato_lista.py mc=<mc.csv> gmin=1.0 gmax=1.10 procs=4 saida=<exato.csv>

Os campos sao os mesmos do mc.py (mesma semente, mesmo h e hD; conferencia._uma). Grava uma linha por realizacao
assim que termina (retomada: as ja gravadas sao puladas). No fim: quantas cairam abaixo de 1 e a P_f corrigida
(Gamma = min(rapido, exato) nas realizacoes refeitas)."""
import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import csv  # noqa: E402
import math  # noqa: E402
import sys  # noqa: E402

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)
import campos as cp  # noqa: E402
import conferencia as cf  # noqa: E402
import mc  # noqa: E402

COL = ["real", "G_rapido", "G_exato", "rel", "cl_rapido", "cl_exato", "xA_rapido", "xA_exato", "xB_rapido",
       "xB_exato", "nav_rapido", "nav_exato", "t_campos", "t_darcy", "t_rapido", "t_exato", "t_malha", "fallbacks",
       "n_max_regra", "viol_rapido", "viol_exato"]


def resumo(mc_csv, saida):
    todas = {int(r["real"]): float(r["Gamma"]) for r in csv.DictReader(open(mc_csv))}
    ex = {int(r["real"]): r for r in csv.DictReader(open(saida))} if os.path.exists(saida) else {}
    S = len(todas)
    f0 = sum(g < 1 for g in todas.values())
    novas = [i for i, r in ex.items() if float(r["G_exato"]) < 1 <= todas[i]]
    perdidas = [i for i, r in ex.items() if float(r["G_exato"]) >= 1 > todas[i]]
    acima = [float(r["rel"]) for r in ex.values()]
    f1 = f0 + len(novas)
    pf = f1 / S
    print(f"refeitas com a busca exata: {len(ex)}")
    if acima:
        n5 = sum(x > 0.005 for x in acima)
        n1 = sum(x > 0.001 for x in acima)
        print(f"rapido/exato - 1: max {max(acima):+.2e}, media {sum(acima) / len(acima):+.2e}; acima de 0,1%: {n1},"
              f" acima de 0,5%: {n5}")
    print(f"falhas: rapido {f0}; novas com o exato (Gamma_exato < 1 <= Gamma_rapido): {len(novas)}"
          f"; exato >= 1 com rapido < 1: {len(perdidas)}")
    print(f"P_f rapido = {f0 / S:.4%}; P_f corrigida = {pf:.4%} (+- {1.96 * math.sqrt(pf * (1 - pf) / S):.2%})"
          f"; Ceron 11,50%")


def main(argv):
    a = dict(mc="", gmin="1.0", gmax="1.10", procs="4", saida="", semente="2025", resumo="0")
    for s in argv:
        k, _, v = s.partition("=")
        if k not in a:
            raise SystemExit(f"opcao desconhecida: {k}")
        a[k] = v
    if a["resumo"] == "1":
        return resumo(a["mc"], a["saida"])
    lin = list(csv.DictReader(open(a["mc"])))
    alvo = sorted(int(r["real"]) for r in lin if float(a["gmin"]) <= float(r["Gamma"]) < float(a["gmax"]))
    feitas = set()
    if os.path.exists(a["saida"]):
        feitas = {int(r["real"]) for r in csv.DictReader(open(a["saida"]))}
    falta = [i for i in alvo if i not in feitas]
    print(f"{len(alvo)} realizacoes com {a['gmin']} <= Gamma < {a['gmax']}; {len(falta)} a fazer", flush=True)
    h, hD = 0.05, 0.25
    arq = mc._arquivo_partida(h, hD)
    P0 = cp.Problema(h=h, hD_janela=hD)
    P0.hD_janela = hD
    mc._confere_partida(P0, arq)
    del P0
    import multiprocessing as mp
    novo = not os.path.exists(a["saida"])
    with open(a["saida"], "a", newline="") as f, \
            mp.get_context("fork").Pool(int(a["procs"]), initializer=mc._inicia, initargs=(h, hD, arq)) as pool:
        w = csv.DictWriter(f, fieldnames=COL)
        if novo:
            w.writeheader()
        for k, r in enumerate(pool.imap_unordered(cf._uma, [(i, int(a["semente"])) for i in falta]), 1):
            w.writerow(r)
            f.flush()
            if k % 20 == 0:
                print(f"  {k}/{len(falta)}", flush=True)
    resumo(a["mc"], a["saida"])


if __name__ == "__main__":
    main(sys.argv[1:])
