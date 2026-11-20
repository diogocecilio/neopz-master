#!/usr/bin/env python3
"""pareado_gim.py: a comparacao PAREADA SRM x GIM do Cho2010 (campanha etapa3b_gim do cluster): as mesmas realizacoes
(ou as mesmas linhas do line sampling) com o fator pelo SRM (FS) e pelo GIM (Gamma). A falha e a mesma em teoria
(FS < 1 sse Gamma < 1); as diferencas vem da discretizacao e da numerica. Sem numpy (roda no fermi).

Uso:
  python3 -I pareado_gim.py mc <pasta SRM> <pasta GIM> [saida=<csv>]
      <pasta>: as partes mc_<a>_<b>.csv de um caso (Monte Carlo direto ou ISMC; real >= 0). Confere que o campo e o
      mesmo (log_peso igual; c_med e phi_med, medias na malha final adaptativa, a menos de 2%), da as falhas e a P_f de cada um nas realizacoes em
      comum (com o peso exp(log_peso) na ISMC), os pares discordantes (so um falha), o expoente k de Gamma ~ FS^k
      (minimos quadrados de ln Gamma = k ln FS) e grava o <csv> (padrao: <pasta GIM>/pareado.csv) com real, FS_srm,
      Gamma_gim, c_med, phi_med, log_peso: a tabela para juntar com a analise limite.
  python3 -I pareado_gim.py ls <pasta SRM> <pasta GIM>
      <pasta>: a saida do line_sampling.py (linhas.csv). A linha i usa o mesmo Z nos dois (semente 7): compara o c_i
      (a distancia ate a superficie de falha ao longo da direcao) linha a linha e a P_f nas linhas em comum.
Uma realizacao sem fator (situacao diferente de ok) fica de fora, exceto Gamma = 0 no GIM (falha)."""
import csv
import glob
import math
import os
import re
import sys


def le_mc(pasta):
    out, ruins = {}, set()
    arqs = [a for a in glob.glob(os.path.join(pasta, "mc_*_*.csv")) if re.search(r"/mc_-?\d+_-?\d+\.csv$", a)]
    if not arqs:
        sys.exit(f"nenhum mc_<a>_<b>.csv em {pasta}")
    for arq in arqs:
        for r in csv.DictReader(open(arq, newline="")):
            n = int(r["real"])
            if n < 0:
                continue
            x = float(r["FS_min"])
            gim0 = r.get("metodo") == "gim" and r["situacao"] == "fator" and x == 0.
            if not gim0 and (r["situacao"] != "ok" or not x > 0):
                ruins.add(n)
                continue
            out[n] = dict(f=x, c=float(r["c_med"]), phi=float(r["phi_med"]), lw=float(r.get("log_peso") or 0.),
                          s=float(r.get("s_fator") or "nan"), metodo=r.get("metodo", ""))
    for n in ruins:
        out.pop(n, None)
    return out, len(ruins)


def pf(vals):
    """P_f = media de I w e o CoV"""
    n = len(vals)
    if n == 0:
        return math.nan, math.nan
    m = sum(vals) / n
    if n < 2 or m <= 0:
        return m, math.inf
    v = sum((x - m) ** 2 for x in vals) / (n - 1)
    return m, math.sqrt(v / n) / m


def mc(pa, pb, saida):
    A, ra = le_mc(pa)
    B, rb = le_mc(pb)
    comuns = sorted(set(A) & set(B))
    print(f"SRM: {pa} ({len(A)} realizacoes, {ra} sem fator)\nGIM: {pb} ({len(B)} realizacoes, {rb} sem fator)")
    print(f"em comum: {len(comuns)}" + (f" (real {comuns[0]} a {comuns[-1]})" if comuns else ""))
    if not comuns:
        return
    met = {A[comuns[0]]["metodo"], B[comuns[0]]["metodo"]}
    if met != {"srm", "gim"}:
        print(f"AVISO: metodos {sorted(met)} (esperado srm na 1a pasta e gim na 2a)")
    # o mesmo campo? O log_peso (ISMC) e funcao so dos xi: tem de ser igual. O c_med e o phi_med sao medias nos pontos
    # de integracao da malha FINAL, que o refinamento adaptativo faz diferente no SRM e no GIM: diferem ~0,1 a 0,5%
    rel = lambda x, y: abs(x - y) / max(1e-12, abs(x), abs(y))
    dmax = max(max(rel(A[n]["c"], B[n]["c"]), rel(A[n]["phi"], B[n]["phi"])) for n in comuns)
    dif = [n for n in comuns if abs(A[n]["lw"] - B[n]["lw"]) > 1e-6 * max(1., abs(A[n]["lw"]))
           or rel(A[n]["c"], B[n]["c"]) > 0.02 or rel(A[n]["phi"], B[n]["phi"]) > 0.02]
    print(f"campo: log_peso igual e c_med, phi_med a menos de 2% em {len(comuns) - len(dif)} de {len(comuns)} "
          f"(maior diferenca relativa de c_med/phi_med: {100 * dmax:.2f}%, efeito da malha adaptativa)"
          + (f"; DIFERENTES: {dif[:10]}{' ...' if len(dif) > 10 else ''}" if dif else " -> pareamento OK"))
    fa = {n for n in comuns if A[n]["f"] < 1.}
    fb = {n for n in comuns if B[n]["f"] < 1.}
    pesada = any(A[n]["lw"] != 0. for n in comuns)
    pA, cA = pf([math.exp(A[n]["lw"]) * (n in fa) for n in comuns])
    pB, cB = pf([math.exp(B[n]["lw"]) * (n in fb) for n in comuns])
    print(f"falhas: SRM (FS < 1) {len(fa)}, GIM (Gamma < 1) {len(fb)}, as duas {len(fa & fb)}")
    print(f"P_f{' (ponderada, ISMC)' if pesada else ''}: SRM {100 * pA:.4f}% (CoV {100 * cA:.1f}%), "
          f"GIM {100 * pB:.4f}% (CoV {100 * cB:.1f}%)")
    so_b, so_a = sorted(fb - fa), sorted(fa - fb)
    print(f"discordantes: so o GIM falha {len(so_b)}, so o SRM falha {len(so_a)}")
    for nome, conj in (("so GIM", so_b), ("so SRM", so_a)):
        if conj:
            print(f"  {nome}: " + ", ".join(f"{n} (FS {A[n]['f']:.4f}, G {B[n]['f']:.4f})" for n in conj[:15])
                  + (" ..." if len(conj) > 15 else ""))
    # Gamma ~ FS^k
    xs = [(math.log(A[n]["f"]), math.log(B[n]["f"])) for n in comuns if B[n]["f"] > 0 and abs(A[n]["f"] - 1.) > 0.01]
    if xs:
        k = sum(x * y for x, y in xs) / sum(x * x for x, _ in xs)
        perto = [(math.log(A[n]["f"]), math.log(B[n]["f"])) for n in comuns
                 if B[n]["f"] > 0 and 0.9 < A[n]["f"] < 1.1 and abs(A[n]["f"] - 1.) > 0.005]
        kp = sum(x * y for x, y in perto) / sum(x * x for x, _ in perto) if perto else math.nan
        print(f"Gamma ~ FS^k: k = {k:.2f} (todas com |FS - 1| > 0,01), k = {kp:.2f} perto de 1 (0,9 < FS < 1,1; "
              f"{len(perto)} pares)")
    for nome, faixa in (("FS < 1,05", lambda f: f < 1.05), ("FS < 1,10", lambda f: f < 1.10)):
        sel = [n for n in comuns if faixa(A[n]["f"])]
        print(f"  realizacoes com {nome}: {len(sel)}; Gamma < 1 entre elas: {sum(1 for n in sel if n in fb)}")
    ta = [A[n]["s"] for n in comuns if not math.isnan(A[n]["s"])]
    tb = [B[n]["s"] for n in comuns if not math.isnan(B[n]["s"])]
    if ta and tb:
        print(f"tempo do fator, media: SRM {sum(ta) / len(ta):.1f} s, GIM {sum(tb) / len(tb):.1f} s")
    with open(saida, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["real", "FS_srm", "Gamma_gim", "c_med", "phi_med", "log_peso"])
        for n in comuns:
            w.writerow([n, A[n]["f"], B[n]["f"], A[n]["c"], A[n]["phi"], A[n]["lw"]])
    print(f"tabela pareada: {saida}")


def le_ls(pasta):
    arq = os.path.join(pasta, "linhas.csv")
    if not os.path.exists(arq):
        sys.exit(f"sem {arq}")
    return {int(r["linha"]): (float(r["c"]), int(r["analises"])) for r in csv.DictReader(open(arq, newline=""))}


def ls(pa, pb):
    A, B = le_ls(pa), le_ls(pb)
    comuns = sorted(set(A) & set(B))
    print(f"SRM: {pa} ({len(A)} linhas)\nGIM: {pb} ({len(B)} linhas)\nem comum: {len(comuns)}")
    if not comuns:
        return
    cauda = lambda c: 0.5 * math.erfc(c / math.sqrt(2.))
    for nome, D in (("SRM", A), ("GIM", B)):
        p, cv = pf([cauda(D[i][0]) for i in comuns])
        na = sum(D[i][1] for i in comuns) / len(comuns)
        print(f"  {nome}: P_f {100 * p:.4f}% (CoV {100 * cv:.1f}%), analises por linha {na:.2f}")
    d = [B[i][0] - A[i][0] for i in comuns]
    m = sum(d) / len(d)
    s = math.sqrt(sum((x - m) ** 2 for x in d) / (len(d) - 1)) if len(d) > 1 else 0.
    print(f"  c_GIM - c_SRM: media {m:+.4f}, dp {s:.4f}, max |dif| {max(abs(x) for x in d):.4f}")
    for i in comuns[:10]:
        print(f"    linha {i:4d}: c SRM {A[i][0]:.4f}  c GIM {B[i][0]:.4f}")


def main():
    if len(sys.argv) < 4 or sys.argv[1] not in ("mc", "ls"):
        sys.exit(__doc__)
    modo, pa, pb = sys.argv[1:4]
    op = dict(a.split("=", 1) for a in sys.argv[4:] if "=" in a)
    if modo == "mc":
        mc(pa, pb, op.get("saida", os.path.join(pb, "pareado.csv")))
    else:
        ls(pa, pb)


if __name__ == "__main__":
    main()
