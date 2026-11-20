#!/usr/bin/env python3
"""line_sampling.py: line sampling (amostragem por linhas) na direcao do FORM, com o campo completo.

    python3 line_sampling.py <exe> <form.json> <pasta> [n=100] [semente=7] [paralelos=14] [tol=0.01] [K=2000]

A direcao e = u*/beta vem do FORM (form_kl.py, M termos por campo), completada com zeros ate os K termos de cada campo.
Cada linha: sorteia z ~ N(0, I) com 3K termos (c', phi', k), tira a componente em e (zp = z - (z.e) e) e acha, na reta
zp + c e, o c_i em que Gamma = 1. Como zp e c sao independentes, P_f = E[Phi(-c_i)], e cada linha da p_i = Phi(-c_i):
so a variacao de c_i entre as linhas (a nao linearidade da superficie de falha) entra na variancia. O estimador nao tem
vies para nenhuma direcao; uma direcao ruim so aumenta a dispersao.

A raiz e pela secante em g(c) = ln Gamma, quase linear em c: a 1a analise em c = beta, a 2a no passo de Newton com a
inclinacao s medida no ponto de projeto (2 analises por direcao), e a raiz da secante e aceita sem nova analise quando o
passo for menor que tol (em c). Os mesmos z para qualquer direcao com a mesma semente (numeros aleatorios comuns).

Saidas na pasta: linhas.csv (linha, c, analises, p) e analises.csv (linha, c, Gamma), com retomada; a tabela parcial
(P_f, CoV, linhas para CoV de 5%) sai no terminal a cada 10 linhas.
"""
import concurrent.futures as cf
import csv
import json
import math
import os
import sys
import threading

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from form_kl import comando  # noqa: E402

FIM = (-8., 12.)


def phi_normal(z):
    return 0.5 * math.erfc(-z / math.sqrt(2.))


def main():
    exe, arq_form, pasta = sys.argv[1], sys.argv[2], os.path.abspath(sys.argv[3])
    op = {"n": 100, "semente": 7, "paralelos": 14, "tol": 0.01, "K": 2000}
    for a in sys.argv[4:]:
        k, v = a.split("=", 1)
        op[k] = type(op[k])(v)
    os.makedirs(pasta, exist_ok=True)
    form = json.load(open(arq_form))
    p = form["parametros"]
    M, K = p["M"], op["K"]
    u = np.array(form["u"])
    e = np.zeros(3 * K)
    for f in range(3):  # c', phi', k: os M do FORM nos K do campo
        e[f * K:f * K + M] = u[f * M:(f + 1) * M]
    beta = float(np.linalg.norm(e))
    e /= beta
    trava = threading.Lock()
    reg = open(os.path.join(pasta, "analises.csv"), "a")

    def gamma(x, rotulo, linha, c):
        arq = os.path.join(pasta, rotulo + ".xi")
        np.savetxt(arq, x, fmt="%.17g")
        env = dict(os.environ, PZ_THREADS="0", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
        import subprocess
        subprocess.run(comando(exe, p, arq, rotulo), check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                       env=env, cwd=pasta)
        with open(os.path.join(pasta, rotulo + "_0_0.csv"), newline="") as f:
            r = next(csv.DictReader(f))
        for suf in (".xi", "_0_0.csv", "_mapa_0_0.csv", "_parametros.txt"):
            try:
                os.remove(os.path.join(pasta, rotulo + suf))
            except FileNotFoundError:
                pass
        if r["situacao"] != "ok":
            raise RuntimeError(f"{rotulo}: situacao {r['situacao']}")
        G = float(r["FS_min"])
        with trava:
            reg.write(f"{linha},{c!r},{G!r}\n")
            reg.flush()
        return G

    # a inclinacao de g = ln Gamma ao longo de e, no ponto de projeto (2 analises)
    arq_s = os.path.join(pasta, "inclinacao.json")
    if os.path.exists(arq_s):
        s = json.load(open(arq_s))["s"]
    else:
        with cf.ThreadPoolExecutor(2) as ex:
            gp, gm = ex.map(lambda d: math.log(gamma((beta + d) * e, f"s{'p' if d > 0 else 'm'}", -1, beta + d)), (0.1, -0.1))
        s = (gp - gm) / 0.2
        json.dump({"s": s, "beta": beta}, open(arq_s, "w"))
    print(f"direcao: {arq_form}  beta = {beta:.4f}  P_f FORM = {phi_normal(-beta):.4f}  d lnGamma/dc = {s:.4f}", flush=True)

    rng = np.random.default_rng(op["semente"])
    Z = rng.standard_normal((op["n"], 3 * K))

    def linha(i):
        zp = Z[i] - (Z[i] @ e) * e
        g = lambda c, j: math.log(gamma(zp + c * e, f"l{i}_{j}", i, c))
        ca = beta
        ga = g(ca, 0)
        cb = min(max(ca - ga / s, FIM[0]), FIM[1])
        gb = g(cb, 1)
        n = 2
        cn = min(max(cb - gb * (cb - ca) / (gb - ga), FIM[0]), FIM[1])
        while abs(cn - cb) > op["tol"] and n < 8:
            ca, ga, cb = cb, gb, cn
            gb = g(cb, n)
            n += 1
            cn = min(max(cb - gb * (cb - ca) / (gb - ga), FIM[0]), FIM[1])
        return i, cn, n

    arq_l = os.path.join(pasta, "linhas.csv")
    feitas = {}
    if os.path.exists(arq_l):
        for r in csv.DictReader(open(arq_l)):
            feitas[int(r["linha"])] = (float(r["c"]), int(r["analises"]))
    else:
        open(arq_l, "w").write("linha,c,analises,p\n")
    fl = open(arq_l, "a")

    def tabela():
        c = np.array([v[0] for v in feitas.values()])
        na = np.array([v[1] for v in feitas.values()])
        pi = np.array([phi_normal(-x) for x in c])
        if len(pi) < 2:
            return
        pf, cov = pi.mean(), pi.std(ddof=1) / (math.sqrt(len(pi)) * pi.mean())
        nl5 = len(pi) * (cov / 0.05) ** 2
        print(f"  linhas {len(pi):4d}  analises {na.sum() + 2:5d}  P_f = {pf:.4f}  CoV = {100 * cov:5.2f}%  "
              f"CoV por linha = {100 * cov * math.sqrt(len(pi)):5.1f}%  analises/linha = {na.mean():.2f}  "
              f"linhas p/ CoV 5% = {nl5:5.0f}  analises p/ CoV 5% (sem o FORM) = {2 + na.mean() * nl5:5.0f}", flush=True)

    falta = [i for i in range(op["n"]) if i not in feitas]
    with cf.ThreadPoolExecutor(op["paralelos"]) as ex:
        for fut in cf.as_completed([ex.submit(linha, i) for i in falta]):
            try:
                i, c, n = fut.result()
            except Exception as err:  # uma analise que nao convergiu: a linha fica de fora (e avisa)
                print(f"  linha perdida: {err}", flush=True)
                continue
            feitas[i] = (c, n)
            fl.write(f"{i},{c!r},{n},{phi_normal(-c)!r}\n")
            fl.flush()
            if len(feitas) % 10 == 0:
                tabela()
    tabela()


if __name__ == "__main__":
    main()
