#!/usr/bin/env python3
"""resumo_tempos.py <pasta da campanha> <casos.txt>: o Monte Carlo direto de cada caso, instante a instante, das partes
mc_<a>_<b>.csv (uma linha por realizacao e instante):
  - em cada instante: realizacoes, FS medio, CoV, P(FS < 1) com o intervalo de Wilson de 95% e a P_f do ajuste lognormal;
  - no menor FS dos instantes (a falha em algum instante): a P_f e em que instante o minimo ocorre.
So as realizacoes com todos os instantes ok (situacao ok e FS > 0); as outras sao contadas. Sem pesos: para a amostragem
por importancia, use amostragem_importancia.py estima."""
import csv
import glob
import math
import os
import re
import sys


def wilson(k, n, z=1.96):
    p = k / n
    c = (p + z * z / (2 * n)) / (1 + z * z / n)
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return 100. * p, 100. * (c - h), 100. * (c + h)


def lognormal(v):
    if len(v) < 2 or min(v) <= 0:
        return math.nan
    ls = [math.log(x) for x in v]
    m = sum(ls) / len(ls)
    s = math.sqrt(sum((x - m) ** 2 for x in ls) / (len(ls) - 1))
    return 100. * 0.5 * math.erfc(m / s / math.sqrt(2.)) if s > 0 else 100. * (m < 0)


def le(pasta):
    fs, ruins, pesos = {}, set(), False
    for arq in glob.glob(os.path.join(pasta, "mc_*_*.csv")):
        if not re.search(r"mc_-?\d+_-?\d+\.csv$", arq):
            continue
        for r in csv.DictReader(open(arq, newline="")):
            n = int(r["real"])
            if n < 0:
                continue
            x = float(r["FS_min"])
            pesos = pesos or float(r.get("log_peso") or 0.) != 0.
            if r["situacao"] != "ok" or not x > 0:
                ruins.add(n)
            else:
                fs.setdefault(n, {})[float(r["t_h"])] = x
    for n in ruins:
        fs.pop(n, None)
    return fs, len(ruins), pesos


def main(pasta, casos):
    for linha in open(casos):
        if not linha.strip() or linha.startswith("#"):
            continue
        cid = linha.split()[0]
        fs, exc, pesos = le(os.path.join(pasta, cid))
        ts = sorted(set().union(*fs.values())) if fs else []
        fs = {n: v for n, v in fs.items() if len(v) == len(ts)}
        print(f"== {cid}: {len(fs)} realizacoes com todos os instantes ok" + (f" ({exc} de fora)" if exc else ""))
        if pesos:
            print("   ATENCAO: o caso tem pesos de amostragem por importancia; os numeros abaixo nao sao a P_f")
        if not fs:
            continue
        print("   t (h)   FS medio   CoV (%)   P(FS < 1) (%)   IC 95% (%)      P_f lognormal (%)")
        for t in ts:
            v = [x[t] for x in fs.values()]
            m = sum(v) / len(v)
            cov = 100. * math.sqrt(sum((x - m) ** 2 for x in v) / (len(v) - 1)) / m if len(v) > 1 else math.nan
            p, lo, hi = wilson(sum(1 for x in v if x < 1.), len(v))
            print(f"   {t:5g}   {m:8.4f}   {cov:7.1f}   {p:13.2f}   {lo:5.1f} a {hi:5.1f}   {lognormal(v):10.2f}")
        mins = {n: min(x.items(), key=lambda a: (a[1], a[0])) for n, x in fs.items()}
        v = [m[1] for m in mins.values()]
        p, lo, hi = wilson(sum(1 for x in v if x < 1.), len(v))
        print(f"   menor FS: media {sum(v) / len(v):.4f}, P_f = {p:.2f}% (IC 95% {lo:.1f} a {hi:.1f}%), "
              f"lognormal {lognormal(v):.2f}%")
        quando = {}
        for t, _ in mins.values():
            quando[t] = quando.get(t, 0) + 1
        print("   instante do menor FS: " + ", ".join(f"{t:g} h: {quando[t]}" for t in sorted(quando)))


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
