#!/usr/bin/env python3
"""amostragem_importancia.py: amostragem por importancia (AI) do Monte Carlo do rebaixamento (SlopeDrawdownMonteCarlo).

A AI sorteia os xi de cada campo de N(mu_c, I) em vez de N(0, I), com mu_c = s_c d: d e a direcao da banda do mecanismo
deterministico (o gradiente da media do campo gaussiano na banda; opcao banda= do programa) e s_c um escalar por campo
(c', phi', k). Cada realizacao leva o peso w = phi(xi) / phi(xi - mu) = exp(-xi . mu + |mu|^2 / 2) (coluna log_peso), e
    P_f = (1/N) sum w_i 1(Gamma_i < 1)
e um estimador SEM VIES de P_f para qualquer mu (a densidade de sorteio e positiva em todo o espaco). O mu so muda a
variancia: o melhor e o centro da regiao de falha. Aqui ele vem de um piloto de Monte Carlo direto, pela entropia
cruzada restrita a direcao d: s_c = media de xi_c . d nas realizacoes do piloto que falharam.

Uso:
  python3 -I amostragem_importancia.py desloc <pasta do piloto> [<projecoes.csv>]
      <pasta do piloto>: as partes mc_<a>_<b>.csv do Monte Carlo direto (FS_min = Gamma), rodado com banda= (as
      colunas proj_*) ou, sem elas, com o <projecoes.csv> do programa com so_projecoes=1 e banda= para as mesmas
      realizacoes. Imprime s_c, s_phi,
      s_k (para desloc=) e quantas falhas o piloto tem.
  python3 -I amostragem_importancia.py estima <pasta da AI> [<pasta do Monte Carlo direto>]
      P_f ponderado, o seu CoV, o tamanho efetivo da amostra (ESS = (sum w)^2 / sum w^2) e, com a segunda pasta, o P_f
      do Monte Carlo direto e o CoV que ele teria com o mesmo numero de realizacoes (o ganho da AI)."""
import csv
import glob
import math
import os
import re
import sys


def linhas(pasta):
    """realizacoes aleatorias (real >= 0) com situacao ok, de todas as partes mc_<a>_<b>.csv da pasta"""
    vistas = {}
    for arq in sorted(glob.glob(os.path.join(pasta, "mc_*_*.csv"))):
        if not re.search(r"mc_-?\d+_-?\d+\.csv$", arq):
            continue  # mapas, projecoes
        for r in csv.DictReader(open(arq, newline="")):
            n = int(r["real"])
            if n >= 0 and r["situacao"] == "ok" and float(r["FS_min"]) > 0:
                vistas[n] = r
    return vistas


def desloc(pasta, arq_proj=None):
    rs = linhas(pasta)
    # as projecoes: do arquivo do so_projecoes ou das proprias linhas do piloto (rodado com banda=)
    proj = {int(r["real"]): r for r in csv.DictReader(open(arq_proj, newline=""))} if arq_proj else rs
    falhas = [n for n, r in rs.items() if float(r["FS_min"]) < 1. and n in proj]
    if not falhas:
        sys.exit("o piloto ainda nao tem falhas com projecao")
    s = [sum(float(proj[n][c]) for n in falhas) / len(falhas) for c in ("proj_c", "proj_phi", "proj_k")]
    usadas = sum(1 for n in rs if n in proj)
    print(f"piloto: {usadas} realizacoes com projecao, {len(falhas)} falhas (P_f = {100. * len(falhas) / usadas:.2f}%)")
    print(f"deslocamento (media de xi . d nas falhas): s_c = {s[0]:+.4f}, s_phi = {s[1]:+.4f}, s_k = {s[2]:+.4f}")
    print(f"opcao do programa: desloc={s[0]:.4f},{s[1]:.4f},{s[2]:.4f}")
    print(f"|mu|^2 = {sum(v * v for v in s):.3f}: a variancia dos pesos e exp(|mu|^2) - 1 = {math.exp(sum(v * v for v in s)) - 1:.2f}")


def estima(pasta, direta=None):
    rs = linhas(pasta)
    if not rs:
        sys.exit(f"nenhuma realizacao em {pasta}")
    w = [math.exp(float(r["log_peso"])) for r in rs.values()]
    f = [1. if float(r["FS_min"]) < 1. else 0. for r in rs.values()]
    n = len(w)
    wf = [a * b for a, b in zip(w, f)]
    pf = sum(wf) / n
    var = sum((x - pf) ** 2 for x in wf) / (n - 1) if n > 1 else math.inf
    cov = math.sqrt(var / n) / pf if pf > 0 else math.inf
    ess = sum(w) ** 2 / sum(x * x for x in w)
    print(f"AI: N = {n}, falhas sorteadas = {int(sum(f))} ({100. * sum(f) / n:.1f}%), media dos pesos = {sum(w) / n:.3f}, "
          f"ESS = {ess:.0f}")
    print(f"    P_f = {100. * pf:.3f}%  CoV(P_f) = {100. * cov:.1f}%  IC 95% ~ [{100. * pf * (1 - 1.96 * cov):.3f}; "
          f"{100. * pf * (1 + 1.96 * cov):.3f}]%")
    if direta:
        d = linhas(direta)
        nd = len(d)
        nf = sum(1 for r in d.values() if float(r["FS_min"]) < 1.)
        p = nf / nd
        covd = math.sqrt((1 - p) / (nd * p)) if p > 0 else math.inf
        covd_n = math.sqrt((1 - p) / (n * p)) if p > 0 else math.inf
        print(f"direto: N = {nd}, falhas = {nf}, P_f = {100. * p:.3f}%  CoV(P_f) = {100. * covd:.1f}%")
        print(f"    com N = {n} o direto teria CoV(P_f) = {100. * covd_n:.1f}%: a AI precisa de "
              f"{(cov / covd_n) ** 2:.2f} vezes as realizacoes do direto para o mesmo CoV")
        dif = (pf - p) / math.sqrt((pf * cov) ** 2 + (p * covd) ** 2) if pf > 0 and p > 0 else math.nan
        print(f"    diferenca AI - direto = {100. * (pf - p):+.3f} pontos percentuais ({dif:+.2f} desvios-padrao)")


if __name__ == "__main__":
    if len(sys.argv) >= 3 and sys.argv[1] == "desloc":
        desloc(sys.argv[2], sys.argv[3] if len(sys.argv) > 3 else None)
    elif len(sys.argv) >= 3 and sys.argv[1] == "estima":
        estima(sys.argv[2], sys.argv[3] if len(sys.argv) > 3 else None)
    else:
        sys.exit(__doc__)
