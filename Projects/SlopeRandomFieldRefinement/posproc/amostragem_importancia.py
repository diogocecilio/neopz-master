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
      do Monte Carlo direto e o CoV que ele teria com o mesmo numero de realizacoes (o ganho da AI). Com o deslocamento
      calculado pelo programa (opcao piloto=, as campanhas ceron e chuva_ai do cluster), confere tambem que todos os
      blocos usaram o mesmo (arquivos mc_desloc_<a>_<b>.txt). Com varios instantes (chuva), a falha e FS < 1 em algum
      instante, e ha tambem a P_f de cada instante, com os mesmos pesos; contra o direto, nos instantes em comum."""
import csv
import glob
import math
import os
import re
import sys


def linhas(pasta):
    """realizacoes aleatorias (real >= 0) de todas as partes mc_<a>_<b>.csv da pasta com todos os instantes ok (situacao
    ok e FS_min > 0): por realizacao, a primeira linha (peso e projecoes), o FS de cada instante em fs_t e o menor deles
    em FS_min (no rebaixamento ha um instante so; na chuva, a falha e FS < 1 em algum instante)"""
    vistas, ruins = {}, set()
    for arq in sorted(glob.glob(os.path.join(pasta, "mc_*_*.csv"))):
        if not re.search(r"mc_-?\d+_-?\d+\.csv$", arq):
            continue  # mapas, projecoes
        for r in csv.DictReader(open(arq, newline="")):
            n = int(r["real"])
            if n < 0:
                continue
            fs = float(r["FS_min"])
            if r["situacao"] != "ok" or not fs > 0:
                ruins.add(n)
                continue
            vistas.setdefault(n, dict(r, fs_t={}))["fs_t"][float(r.get("t_h") or 0.)] = fs
    EXCLUIDAS[pasta] = {n: vistas.pop(n, {}).get("log_peso") for n in ruins}
    for v in vistas.values():
        v["FS_min"] = repr(min(v["fs_t"].values()))
    return vistas


EXCLUIDAS = {}  # pasta -> {realizacao com algum instante sem fator: log_peso (None se nenhum instante ok)}


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


def deslocamentos(pasta):
    """os deslocamentos que os blocos com piloto= gravaram (mc_desloc_<a>_<b>.txt): tem de ser o mesmo em todos"""
    arqs = sorted(glob.glob(os.path.join(pasta, "mc_desloc_*_*.txt")))
    if not arqs:
        return
    textos = {}
    for a in arqs:
        textos.setdefault("".join(l for l in open(a) if not l.startswith("#")), []).append(os.path.basename(a))
    d = dict(l.split("=", 1) for l in next(iter(textos)).splitlines() if "=" in l)
    print(f"deslocamento ({d.get('ai')}, pelo piloto): {len(arqs)} blocos, "
          + ("o mesmo em todos" if len(textos) == 1 else f"ATENCAO: {len(textos)} deslocamentos diferentes") + "; "
          f"piloto com {d.get('validas')} realizacoes validas, {d.get('falhas')} falhas, elite de {d.get('elite')} "
          f"(fator < {float(d.get('nivel', 'nan')):.4f}), |mu|^2 = {float(d.get('mu2', 'nan')):.3f}"
          + (f", s = {d['s']}" if "s" in d else ""))


def evento(rs, ts):
    """por realizacao: o menor FS nos instantes ts (sem as realizacoes que nao tem algum deles)"""
    out = {}
    for n, r in rs.items():
        v = [r["fs_t"].get(t) for t in ts]
        if all(x is not None for x in v):
            out[n] = min(v)
    return out


def pf_ponderada(rs, fs):
    """P_f = soma(peso * falha) / N nas realizacoes de fs, o CoV dela, as falhas sorteadas e N (peso 1 sem log_peso)"""
    q = [math.exp(float(rs[n].get("log_peso") or 0.)) * (x < 1.) for n, x in fs.items()]
    n = len(q)
    pf = sum(q) / n if n else 0.
    var = sum((x - pf) ** 2 for x in q) / (n - 1) if n > 1 else math.inf
    cov = math.sqrt(var / n) / pf if pf > 0 else math.inf
    return pf, cov, sum(1 for x in fs.values() if x < 1.), n


def instantes(rs):
    return sorted(set().union(*(r["fs_t"] for r in rs.values())))


def rotulo(ts):
    return f"FS em {ts[0]:g} h" if len(ts) == 1 else "menor FS de " + ", ".join(f"{t:g}" for t in ts) + " h"


def estima(pasta, direta=None):
    rs = linhas(pasta)
    if not rs:
        sys.exit(f"nenhuma realizacao em {pasta}")
    deslocamentos(pasta)
    ts = instantes(rs)
    w = [math.exp(float(r.get("log_peso") or 0.)) for r in rs.values()]
    ess = sum(w) ** 2 / sum(x * x for x in w)
    pf, cov, nf, n = pf_ponderada(rs, evento(rs, ts))
    print(f"AI: N = {n}, falhas sorteadas = {nf} ({100. * nf / n:.1f}%), media dos pesos = {sum(w) / len(w):.3f}, "
          f"ESS = {ess:.0f}")
    print(f"    P_f ({rotulo(ts)}) = {100. * pf:.3f}%  CoV(P_f) = {100. * cov:.1f}%  IC 95% ~ "
          f"[{100. * pf * (1 - 1.96 * cov):.3f}; {100. * pf * (1 + 1.96 * cov):.3f}]%")
    exc = EXCLUIDAS.get(pasta, {})
    if exc:  # sem fator em algum instante: a P_f acima supoe que elas falham na mesma proporcao das outras
        we = sum(math.exp(float(x or 0.)) for x in exc.values())
        nt = n + len(exc)
        print(f"    ATENCAO: {len(exc)} realizacoes sem fator em algum instante (peso total {we:.3g}) ficaram de fora; "
              f"com elas seguras, P_f = {100. * pf * n / nt:.3f}%; com elas falhas, {100. * (pf * n + we) / nt:.3f}%")
    if len(ts) > 1:
        print("    por instante: " + "  ".join(
            "%g h: %.3f%% (CoV %.0f%%)" % (t, 100. * x[0], 100. * x[1]) for t in ts for x in [pf_ponderada(rs, evento(rs, [t]))]))
    if direta:
        d = linhas(direta)
        comuns = [t for t in ts if t in instantes(d)] if d else []
        if not comuns:
            print(f"direto: sem realizacoes ou sem instantes em comum em {direta}")
            return
        if comuns != ts:  # a AI de novo, so nos instantes do direto
            pf, cov, nf, n = pf_ponderada(rs, evento(rs, comuns))
            print(f"    AI nos instantes do direto ({rotulo(comuns)}): P_f = {100. * pf:.3f}%  CoV(P_f) = {100. * cov:.1f}%")
        p, covd, nfd, nd = pf_ponderada(d, evento(d, comuns))
        print(f"direto: N = {nd}, falhas = {nfd}, P_f ({rotulo(comuns)}) = {100. * p:.3f}%  CoV(P_f) = {100. * covd:.1f}%")
        if not p > 0:
            print("    o direto nao tem falhas: sem comparacao de P_f nem de custo")
            return
        if nfd < 10:
            print(f"    ATENCAO: o direto tem so {nfd} falhas; a comparacao abaixo e muito incerta")
        # o custo do direto com a P_f da AI (bem mais precisa que a de um piloto com poucas falhas)
        covd_n = math.sqrt((1 - pf) / (n * pf)) if pf > 0 else math.inf
        print(f"    com N = {n} e P_f = {100. * pf:.3f}%, o direto teria CoV(P_f) = {100. * covd_n:.1f}%: a AI precisa de "
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
