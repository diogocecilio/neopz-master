#!/usr/bin/env python3
"""Curvas FS(T) e Γ(T) do rebaixamento acoplado (CSV do comando rebaixamento).

Uso: plota_rebaixamento.py saida.png arquivo1.csv [arquivo2.csv ...]
O tempo adimensional T = c_v t / H² é deslocado de T_d/1000 para o eixo logarítmico; a linha tracejada de cada
curva é o valor do fluxo estacionário desacoplado (artigo).
"""
import csv
import math
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        return 1
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.5))
    cores = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    for i, arq in enumerate(sys.argv[2:]):
        linhas = list(csv.DictReader(open(arq)))
        if not linhas:
            continue
        cor = cores[i % len(cores)]
        td = float(linhas[0]["Td"])
        nome = os.path.basename(arq)[:-4]
        est = [r for r in linhas if r["T"] == "inf"]
        pts = [r for r in linhas if r["T"] != "inf" and "colapso_acoplado" not in r["FS_status"]]
        col = [r for r in linhas if r["T"] != "inf" and "colapso_acoplado" in r["FS_status"]]
        for k, (chave, st) in enumerate((("FS", "FS_status"), ("Gamma", "Gamma_status"))):
            ok = [r for r in pts if r[st] in ("ok", "limite_maximo")]
            abaixo = [r for r in pts if r[st] and r[st] not in ("ok", "limite_maximo")]  # instável já no início
            if ok:
                ax[k].semilogx([float(r["T"]) + td / 1000. for r in ok], [float(r[chave]) for r in ok], "o-", ms=3,
                               color=cor, label=nome)
            for r in abaixo:
                ax[k].plot(float(r["T"]) + td / 1000., 1.0, "v", ms=7, color=cor)
            if est and est[0][st]:
                ax[k].axhline(float(est[0][chave]), ls="--", lw=0.8, color=cor)
            for r in col:
                ax[k].plot(float(r["T"]) + td / 1000., 1.0, "x", ms=9, color=cor)
            ax[k].axvline(td, ls=":", lw=0.6, color=cor)
    for a, t in zip(ax, ("FS (redução de resistência)", "Γ (fator de carga)")):
        a.axhline(1.0, color="k", lw=0.8)
        a.set_xlabel("T = c_v t / H²  (pontilhado: fim do rebaixamento; x: colapso no u-p; v: < 1)")
        a.set_ylabel(t)
        a.grid(alpha=0.3, which="both")
    for a in ax:
        if a.lines:
            a.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(sys.argv[1], dpi=130)
    print("figura:", sys.argv[1])
    return 0


if __name__ == "__main__":
    sys.exit(main())
