#!/usr/bin/env python3
"""line_sampling.py: line sampling (amostragem por linhas) na direcao do FORM, com o campo completo.

    python3 line_sampling.py <exe> <form.json> <pasta> [n=100] [semente=7] [paralelos=14] [tol=0.01] [K=2000]
                             [chave=valor das opcoes do programa: nystrom=, metodo=, srm_f0=, ...]

A direcao e = u*/beta vem do FORM (form_kl.py, M termos por campo: c', phi' no Cho2010; c', phi', k no rebaixamento),
completada com zeros ate os K termos de cada campo. As analises usam o campo com K termos (M=K) e as opcoes do FORM,
trocadas pelas chave=valor dadas (ex.: nystrom=<a base completa>, metodo=srm).
Cada linha: sorteia z ~ N(0, I) com nc K termos, tira a componente em e (zp = z - (z.e) e) e acha, na reta zp + c e, o
c_i em que o fator (FS ou Gamma) vale 1. Como zp e c sao independentes, P_f = E[Phi(-c_i)], e cada linha da
p_i = Phi(-c_i): so a variacao de c_i entre as linhas (a nao linearidade da superficie de falha) entra na variancia. O
estimador nao tem vies para nenhuma direcao; uma direcao ruim so aumenta a dispersao.

A raiz e pela secante em g(c) = ln(fator), quase linear em c: a 1a analise em c = beta, a 2a no passo de Newton com a
inclinacao s medida no ponto de projeto (2 analises por direcao), e a raiz da secante e aceita sem nova analise quando o
passo for menor que tol (em c). Os mesmos z para qualquer direcao com a mesma semente (numeros aleatorios comuns).
Sem numpy (roda no fermi).

Saidas na pasta: linhas.csv (linha, c, analises, p), analises.csv (linha, c, fator) e ls.json (o resultado), com
retomada; a tabela parcial (P_f, CoV, linhas para CoV de 5%) sai no terminal a cada 10 linhas.
"""
import concurrent.futures as cf
import csv
import json
import math
import os
import random
import subprocess
import sys
import threading

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from form_kl import campos, comando  # noqa: E402

FIM = (-8., 12.)


def cauda(c):
    """Phi(-c)"""
    return 0.5 * math.erfc(c / math.sqrt(2.))


def main():
    if len(sys.argv) < 4:
        sys.exit(__doc__)
    exe, arq_form, pasta = os.path.abspath(sys.argv[1]), sys.argv[2], os.path.abspath(sys.argv[3])
    op = {"n": 100, "semente": 7, "paralelos": 14, "tol": 0.01, "K": 2000}
    form = json.load(open(arq_form))
    p = dict(form["parametros"])
    p.setdefault("prog", "cho")          # form.json de versoes antigas do form_kl.py (antes de prog=)
    p.setdefault("transiente", "tri3")
    for a in sys.argv[4:]:
        k, v = a.split("=", 1)
        if k in op:
            op[k] = type(op[k])(v)
        else:
            p[k] = type(p[k])(v) if k in p and not isinstance(p[k], str) else v
    os.makedirs(pasta, exist_ok=True)
    nc = len(campos(p))
    u = form["u"]
    M, K = len(u) // nc, op["K"]
    p["M"] = K  # as analises com o campo de K termos
    dim = nc * K
    e = [0.] * dim
    for f in range(nc):  # os M do FORM nos K do campo
        e[f * K:f * K + M] = u[f * M:(f + 1) * M]
    beta = math.sqrt(sum(x * x for x in e))
    e = [x / beta for x in e]
    trava = threading.Lock()
    reg = open(os.path.join(pasta, "analises.csv"), "a")

    def fator(x, rotulo, linha, c):
        arq = os.path.join(pasta, rotulo + ".xi")
        with open(arq, "w") as f:
            f.write("\n".join(repr(t) for t in x) + "\n")
        env = dict(os.environ, PZ_THREADS="0", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
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

    def ponto(zp, c):
        return [a + c * b for a, b in zip(zp, e)]

    # a inclinacao de g = ln(fator) ao longo de e, no ponto de projeto (2 analises)
    arq_s = os.path.join(pasta, "inclinacao.json")
    if os.path.exists(arq_s):
        s = json.load(open(arq_s))["s"]
    else:
        zero = [0.] * dim
        with cf.ThreadPoolExecutor(2) as ex:
            gp, gm = ex.map(lambda d: math.log(fator(ponto(zero, beta + d), f"s{'p' if d > 0 else 'm'}", -1, beta + d)),
                            (0.1, -0.1))
        s = (gp - gm) / 0.2
        json.dump({"s": s, "beta": beta}, open(arq_s, "w"))
    print(f"direcao: {arq_form}  beta = {beta:.4f}  P_f FORM = {cauda(beta):.4e}  d ln(fator)/dc = {s:.4f}  "
          f"(campo com {K} termos, metodo {p.get('metodo', '-')})", flush=True)

    rng = random.Random(op["semente"])
    Z = [[rng.gauss(0., 1.) for _ in range(dim)] for _ in range(op["n"])]

    def linha(i):
        ze = sum(a * b for a, b in zip(Z[i], e))
        zp = [a - ze * b for a, b in zip(Z[i], e)]
        g = lambda c, j: math.log(fator(ponto(zp, c), f"l{i}_{j}", i, c))
        # g > 0: seguro; g < 0: falha (a falha fica do lado de c grande). Secante protegida: com a raiz cercada
        # (seg < raiz < fal), o candidato fora do intervalo vira bissecao; sem cerca, o candidato fica em FIM e, se
        # chegar em FIM[1] ainda seguro (ou em FIM[0] ja na falha), a linha toda e segura (p = 0) ou falha (p = 1)
        seg, fal = None, None

        def anota(c, gc):
            nonlocal seg, fal
            if gc > 0. and (seg is None or c > seg):
                seg = c
            if gc <= 0. and (fal is None or c < fal):
                fal = c

        ca = beta
        ga = g(ca, 0)
        anota(ca, ga)
        cb = min(max(ca - ga / s, FIM[0]), FIM[1])
        gb = g(cb, 1)
        anota(cb, gb)
        n = 2
        while True:
            if fal is None and seg is not None and seg >= FIM[1]:
                return i, FIM[1], n            # segura ate FIM[1]
            if seg is None and fal is not None and fal <= FIM[0]:
                return i, FIM[0], n            # falha desde FIM[0]
            if gb == 0.:
                return i, cb, n
            cn = cb - gb * (cb - ca) / (gb - ga) if gb != ga else cb + (1. if gb > 0 else -1.)
            cn = min(max(cn, FIM[0]), FIM[1])
            if seg is not None and fal is not None and not (min(seg, fal) <= cn <= max(seg, fal)):
                cn = 0.5 * (seg + fal)         # bissecao dentro da cerca
            if abs(cn - cb) <= op["tol"] or n >= 10:
                return i, cn, n
            ca, ga, cb = cb, gb, cn
            gb = g(cb, n)
            anota(cb, gb)
            n += 1

    arq_l = os.path.join(pasta, "linhas.csv")
    feitas = {}
    if os.path.exists(arq_l):
        for r in csv.DictReader(open(arq_l)):
            feitas[int(r["linha"])] = (float(r["c"]), int(r["analises"]))
    else:
        open(arq_l, "w").write("linha,c,analises,p\n")
    fl = open(arq_l, "a")
    perdidas = [0]

    def tabela(final=False):
        c = [v[0] for v in feitas.values()]
        na = [v[1] for v in feitas.values()]
        pi = [cauda(x) for x in c]
        if len(pi) < 2:
            return
        m = sum(pi) / len(pi)
        sd = math.sqrt(sum((x - m) ** 2 for x in pi) / (len(pi) - 1))
        cov = sd / (math.sqrt(len(pi)) * m)
        nl5 = len(pi) * (cov / 0.05) ** 2
        npl = sum(na) / len(na)
        print(f"  linhas {len(pi):4d}  analises {sum(na) + 2:5d}  P_f = {m:.4e}  CoV = {100 * cov:5.2f}%  "
              f"CoV por linha = {100 * cov * math.sqrt(len(pi)):5.1f}%  analises/linha = {npl:.2f}  "
              f"linhas p/ CoV 5% = {nl5:5.0f}  analises p/ CoV 5% (sem o FORM) = {2 + npl * nl5:5.0f}", flush=True)
        if final:
            json.dump({"form": os.path.abspath(arq_form), "beta": beta, "s": s, "linhas": len(pi), "perdidas": perdidas[0],
                       "P_f": m, "CoV": cov, "CoV_linha": cov * math.sqrt(len(pi)), "analises": sum(na) + 2,
                       "analises_linha": npl, "analises_form": form.get("avaliacoes"),
                       "analises_cov5": 2 + npl * nl5, "parametros": p},
                      open(os.path.join(pasta, "ls.json"), "w"), indent=1)

    falta = [i for i in range(op["n"]) if i not in feitas]
    with cf.ThreadPoolExecutor(op["paralelos"]) as ex:
        for fut in cf.as_completed([ex.submit(linha, i) for i in falta]):
            try:
                i, c, n = fut.result()
            except Exception as err:  # uma analise que nao convergiu: a linha fica de fora (e avisa)
                perdidas[0] += 1
                print(f"  linha perdida: {err}", flush=True)
                continue
            feitas[i] = (c, n)
            fl.write(f"{i},{c!r},{n},{cauda(c)!r}\n")
            fl.flush()
            if len(feitas) % 10 == 0:
                tabela()
    tabela(final=True)


if __name__ == "__main__":
    main()
