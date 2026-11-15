#!/usr/bin/env python3
"""compara_pares.py: dois Monte Carlo diretos com as MESMAS realizacoes (mesma semente: os mesmos campos), comparados
realizacao a realizacao. Serve para o SRM antigo contra o novo (etapa 1 do plano de trabalho: srm_f0=1 contra 0,8), a
mesma versao em duas maquinas ou duas listas de instantes.

Uso: python3 -I compara_pares.py <pasta A> <pasta B> [t=<h>] [rotulos=<A>,<B>]
  <pasta>: as partes mc_<a>_<b>.csv de um caso (colunas real, t_h, FS_min, situacao).
  t=<h>: so o instante t_h = h (padrao: o menor FS dos instantes, como no status.sh).
Uma realizacao entra se tem fator (situacao ok, FS > 0) nas duas pastas. Imprime:
  - a diferenca FS_B - FS_A: longe de 1 (|FS_A - 1| > 0,03) ela e so o ruido do modelo; perto de 1, o efeito;
  - as falhas (FS < 1) de cada um nas realizacoes em comum e os pares discordantes, com o teste de McNemar exato
    (binomial com p = 1/2 nos discordantes);
  - o histograma de FS de 0,990 a 1,010 em faixas de 0,001: o buraco do SRM antigo e uma faixa vazia logo abaixo de 1
    e um acumulo logo acima;
  - o tempo do fator (coluna s_fator, somada nos instantes), longe e perto de 1: o SRM novo anda o caminho de F desde
    0,8, e perto de 1 o antigo parava logo em F = 1."""
import csv
import glob
import math
import os
import re
import sys


def le(pasta, t):
    """{realizacao: FS} das realizacoes com fator em todos os instantes (ou no instante t), {realizacao: s do fator}"""
    fs, seg, ruins = {}, {}, set()
    arqs = [a for a in glob.glob(os.path.join(pasta, "mc_*_*.csv")) if re.search(r"/mc_-?\d+_-?\d+\.csv$", a)]
    if not arqs:
        sys.exit(f"nenhum mc_<a>_<b>.csv em {pasta}")
    for arq in arqs:
        for r in csv.DictReader(open(arq, newline="")):
            n = int(r["real"])
            if n < 0:
                continue
            th = float(r.get("t_h") or 0.)
            if t is not None and abs(th - t) > 1e-9:
                continue
            x = float(r["FS_min"])
            if r["situacao"] != "ok" or not x > 0:
                ruins.add(n)
                continue
            fs[n] = min(x, fs.get(n, math.inf))
            try:
                seg[n] = seg.get(n, 0.) + float(r.get("s_fator") or "nan")
            except ValueError:
                seg[n] = math.nan
    for n in ruins:
        fs.pop(n, None)
    return fs, len(ruins), seg


def mcnemar(n01, n10):
    """p bilateral exato: binomial(n01 + n10, 1/2)"""
    m = n01 + n10
    if m == 0:
        return 1.
    k = min(n01, n10)
    cauda = sum(math.comb(m, i) for i in range(k + 1)) / 2. ** m
    return min(1., 2. * cauda)


def quantil(v, q):
    v = sorted(v)
    if not v:
        return float("nan")
    i = q * (len(v) - 1)
    a = int(math.floor(i))
    b = min(a + 1, len(v) - 1)
    return v[a] + (i - a) * (v[b] - v[a])


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    pa, pb = sys.argv[1], sys.argv[2]
    t, rot = None, None
    for a in sys.argv[3:]:
        if a.startswith("t="):
            t = float(a[2:])
        elif a.startswith("rotulos="):
            rot = a[8:].split(",")
        else:
            sys.exit(f"opcao desconhecida: {a}")
    ra, rb = rot if rot else (os.path.basename(os.path.normpath(pa)), os.path.basename(os.path.normpath(pb)))
    A, ruinsA, sA = le(pa, t)
    B, ruinsB, sB = le(pb, t)
    comuns = sorted(set(A) & set(B))
    print(f"A = {ra} ({pa})\nB = {rb} ({pb})")
    print(f"FS: {'instante t = %g h' % t if t is not None else 'o menor dos instantes'}")
    print(f"realizacoes com fator: A {len(A)} (sem fator {ruinsA}), B {len(B)} (sem fator {ruinsB}); em comum {len(comuns)}")
    if not comuns:
        return
    d = [B[n] - A[n] for n in comuns]
    longe = [B[n] - A[n] for n in comuns if abs(A[n] - 1.) > 0.03]
    perto = [B[n] - A[n] for n in comuns if abs(A[n] - 1.) <= 0.03]

    def resumo(v, nome):
        if not v:
            print(f"  {nome}: nenhuma")
            return
        m = sum(v) / len(v)
        s = math.sqrt(sum((x - m) ** 2 for x in v) / (len(v) - 1)) if len(v) > 1 else 0.
        print(f"  {nome} ({len(v)}): media {m:+.5f}, dp {s:.5f}, mediana {quantil(v, .5):+.5f}, "
              f"5%-95% [{quantil(v, .05):+.5f}, {quantil(v, .95):+.5f}], max |dif| {max(abs(x) for x in v):.5f}")
    print("\ndiferenca FS_B - FS_A:")
    resumo(d, "todas")
    resumo(longe, "|FS_A - 1| > 0,03")
    resumo(perto, "|FS_A - 1| <= 0,03")
    iguais = sum(1 for x in d if abs(x) < 5e-7)
    print(f"  iguais (|dif| < 5e-7): {iguais} de {len(d)}")

    fa = {n for n in comuns if A[n] < 1.}
    fb = {n for n in comuns if B[n] < 1.}
    n01, n10 = len(fb - fa), len(fa - fb)
    N = len(comuns)
    print(f"\nfalhas (FS < 1) nas {N} em comum: A {len(fa)} ({100 * len(fa) / N:.2f}%), B {len(fb)} ({100 * len(fb) / N:.2f}%)")
    print(f"  discordantes: so B falha {n01}, so A falha {n10}; McNemar exato p = {mcnemar(n01, n10):.3g}")
    for nome, conj in (("so B falha", sorted(fb - fa)), ("so A falha", sorted(fa - fb))):
        if conj:
            print(f"  {nome}: " + ", ".join(f"{n} ({A[n]:.4f} -> {B[n]:.4f})" for n in conj[:12]) +
                  (" ..." if len(conj) > 12 else ""))

    print("\ntempo do fator por realizacao (s_fator, somado nos instantes), media:")
    for nome, sel in (("todas", comuns), ("|FS_A - 1| > 0,03", [n for n in comuns if abs(A[n] - 1.) > 0.03]),
                      ("|FS_A - 1| <= 0,03", [n for n in comuns if abs(A[n] - 1.) <= 0.03])):
        ta = [sA[n] for n in sel if not math.isnan(sA.get(n, math.nan))]
        tb = [sB[n] for n in sel if not math.isnan(sB.get(n, math.nan))]
        if ta and tb:
            ma, mb = sum(ta) / len(ta), sum(tb) / len(tb)
            print(f"  {nome} ({len(sel)}): A {ma:.1f} s, B {mb:.1f} s, B / A = {mb / ma:.2f}")

    print("\nhistograma de FS perto de 1 (faixas de 0,001):")
    print(f"  {'faixa':>17}  {ra[:10]:>10}  {rb[:10]:>10}")
    for k in range(-10, 10):
        a0, a1 = 1. + k / 1000., 1. + (k + 1) / 1000.
        ca = sum(1 for n in comuns if a0 <= A[n] < a1)
        cb = sum(1 for n in comuns if a0 <= B[n] < a1)
        marca = "  <- logo abaixo de 1" if k == -1 else ("  <- logo acima de 1" if k == 0 else "")
        print(f"  [{a0:.3f}, {a1:.3f})  {ca:10d}  {cb:10d}{marca}")


if __name__ == "__main__":
    main()
