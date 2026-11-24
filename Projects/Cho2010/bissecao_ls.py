#!/usr/bin/env python3
"""bissecao_ls.py: confere por bissecao as raizes de um line sampling pronto (line_sampling.py).

    python3 bissecao_ls.py <exe> <pasta do ls> <pasta de saida> [linhas=10] [lista=] [tol=0.001] [meia=0.25]
                           [extras=0,4] [paralelos=24] [threads=0] [form=<form.json>]

O line_sampling.py acha a raiz c_i de cada linha pela secante em g(c) = ln(fator) e aceita a raiz da secante sem nova
analise quando o passo e menor que tol (0,01); as linhas de raiz baixa sao as que mais pesam na P_f (p_i = Phi(-c_i)).
Este script refaz a raiz de algumas linhas por bissecao cercada, para descartar raizes espurias do passo de Newton:
  - as linhas: as 'linhas' de menor c no linhas.csv (padrao 10), ou as da lista= (indices separados por virgulas);
  - o z de cada linha e regerado como no line_sampling.py (a semente e o K do ls.json; a linha i usa o i-esimo vetor de
    nc K normais), e a direcao e do form.json (o do ls.json, ou form=); as opcoes das analises sao as do ls.json
    (parametros), as mesmas do line sampling;
  - conferencia da reconstrucao: a 1a analise refaz a 1a analise da linha no analises.csv (em c = beta) e compara o
    fator (a diferenca relativa vai para a tabela; se passar de 1e-3, o z ou as opcoes nao sao os mesmos);
  - a cerca: g em c_ls - meia e c_ls + meia (se o sinal nao trocar, a meia largura dobra ate 4 vezes) e nos pontos
    extras= (padrao c = 0 e c = 4: o sinal esperado e seguro em 0 e falha em 4, se a falha for so do lado de c grande);
    essas analises rodam ao mesmo tempo;
  - a bissecao: o meio da cerca ate a largura < tol (padrao 1e-3), uma analise por vez em cada linha (as linhas em
    paralelo); a raiz e o meio da cerca final.
Saidas na pasta de saida: bissecao.csv (linha, c_ls, c_bis, dc = c_bis - c_ls, p_ls, p_bis, a cerca final, as analises
e a conferencia), analises.csv (linha, rotulo, c, fator) e bissecao.json (o resumo: max |dc|, a P_f do LS e a P_f com
as raizes da bissecao no lugar das do LS). Os extras com o sinal inesperado saem na tabela (coluna extras) e no aviso.
Gamma = 0 (sem equilibrio nem com 1% da gravidade) e falha, como no line_sampling.py; na coluna (prog=coluna), o
fator e o do line_sampling.py. Sem numpy (roda no fermi).
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
from form_kl import ambiente, campos, comando, le_fator  # noqa: E402
from line_sampling import GAMMA_ZERO, SUFIXOS_COLUNA, cauda  # noqa: E402


def main():
    if len(sys.argv) < 4:
        sys.exit(__doc__)
    exe, pasta_ls, pasta = os.path.abspath(sys.argv[1]), os.path.abspath(sys.argv[2]), os.path.abspath(sys.argv[3])
    op = {"linhas": 10, "lista": "", "tol": 1e-3, "meia": 0.25, "extras": "0,4", "paralelos": 24, "threads": 0,
          "form": ""}
    for a in sys.argv[4:]:
        k, v = a.split("=", 1)
        if k not in op:
            sys.exit(f"opcao desconhecida: {a}")
        op[k] = type(op[k])(v)
    ls = json.load(open(os.path.join(pasta_ls, "ls.json")))
    p = dict(ls["parametros"])
    p["threads"] = op["threads"]
    arq_form = op["form"] or ls["form"]
    form = json.load(open(arq_form))
    os.makedirs(pasta, exist_ok=True)

    # a direcao, como no line_sampling.py (os M do FORM nos K termos de cada campo; -e se a origem ja falha)
    nc, K = len(campos(p)), ls["K"]
    u = form["u"]
    M, dim = len(u) // nc, nc * K
    e = [0.] * dim
    for f in range(nc):
        e[f * K:f * K + M] = u[f * M:(f + 1) * M]
    beta = math.sqrt(sum(x * x for x in e))
    e = [x / beta for x in e]
    if form.get("beta", 1.) < 0.:
        e = [-x for x in e]
        beta = -beta
    if abs(beta - ls["beta"]) > 1e-9 * abs(beta):
        sys.exit(f"o beta do {arq_form} ({beta}) nao e o do ls.json ({ls['beta']}): outro FORM?")

    # as linhas: as de menor c (ou as da lista), com o c, o p e as analises do line sampling
    feitas = {}
    for r in csv.DictReader(open(os.path.join(pasta_ls, "linhas.csv"))):
        feitas[int(r["linha"])] = float(r["c"])
    if op["lista"]:
        escolha = [int(x) for x in op["lista"].split(",")]
    else:
        escolha = sorted(feitas, key=lambda i: feitas[i])[:op["linhas"]]
    primeira = {}  # a 1a analise de cada linha no analises.csv: (c, fator)
    for r in csv.reader(open(os.path.join(pasta_ls, "analises.csv"))):
        if len(r) == 3 and r[0].lstrip("-").isdigit() and int(r[0]) in escolha and int(r[0]) not in primeira:
            primeira[int(r[0])] = (float(r[1]), float(r[2]))

    # os z das linhas: o i-esimo vetor de dim normais da semente (gerados em ordem, guardados so os escolhidos)
    rng = random.Random(ls["semente"])
    Z = {}
    for i in range(max(escolha) + 1):
        z = [rng.gauss(0., 1.) for _ in range(dim)]
        if i in escolha:
            Z[i] = z
    extras = [float(x) for x in op["extras"].split(",") if x.strip()]
    print(f"bissecao das linhas {escolha} de {pasta_ls}: beta = {beta:.4f}, K = {K}, semente {ls['semente']}, "
          f"tol {op['tol']}, meia {op['meia']}, extras {extras}", flush=True)
    print("  comando: " + " ".join(comando(exe, p, "<x>.xi", "<rotulo>")), flush=True)

    trava = threading.Lock()
    reg = open(os.path.join(pasta, "analises.csv"), "w")
    reg.write("linha,rotulo,c,fator\n")
    analises = cf.ThreadPoolExecutor(op["paralelos"])

    def fator(x, rotulo, linha, c):
        """o fator (Gamma ou FS) da analise nos xi x, como no line_sampling.py"""
        arq = os.path.join(pasta, rotulo + ".xi")
        with open(arq, "w") as f:
            f.write("\n".join(repr(t) for t in x) + "\n")
        subprocess.run(comando(exe, p, arq, rotulo), check=True, stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL, env=ambiente(p), cwd=pasta)
        try:
            if p["prog"] == "chuva":
                G = le_fator(os.path.join(pasta, rotulo + "_0_0.csv"), p)[0]
            elif p["prog"] == "coluna":  # como no line_sampling.py
                G = max(le_fator(os.path.join(pasta, rotulo + "_0_0.csv"), p)[0], GAMMA_ZERO)
            else:
                with open(os.path.join(pasta, rotulo + "_0_0.csv"), newline="") as f:
                    r = next(csv.DictReader(f))
                G = float(r["FS_min"])
                if r["situacao"] == "fator" and G == 0.:
                    G = GAMMA_ZERO
                elif r["situacao"] != "ok":
                    raise RuntimeError(f"{rotulo}: situacao {r['situacao']} (FS_min {r['FS_min']})")
        finally:
            for suf in (".xi", "_0_0.csv", "_mapa_0_0.csv", "_parametros.txt") + \
                    (SUFIXOS_COLUNA if p["prog"] == "coluna" else ()):
                try:
                    os.remove(os.path.join(pasta, rotulo + suf))
                except FileNotFoundError:
                    pass
        with trava:
            reg.write(f"{linha},{rotulo},{c!r},{G!r}\n")
            reg.flush()
        return G

    def linha(i):
        z = Z[i]
        ze = sum(a * b for a, b in zip(z, e))
        zp = [a - ze * b for a, b in zip(z, e)]
        n = [0]

        def g_lote(cs, nome):
            """g = ln(fator) nos c dados, ao mesmo tempo"""
            futs = []
            for j, c in enumerate(cs):
                futs.append(analises.submit(fator, [a + c * b for a, b in zip(zp, e)], f"b{i}_{nome}{n[0] + j}", i, c))
            n[0] += len(cs)
            return [math.log(f.result()) for f in futs]

        c_ls = feitas[i]
        c0, G0 = primeira.get(i, (beta, None))
        meia = op["meia"]
        g = g_lote([c0, c_ls - meia, c_ls + meia] + extras, "c")
        conf = abs(math.exp(g[0]) / G0 - 1.) if G0 else None
        g_ext = g[3:]
        lo, glo, hi, ghi = c_ls - meia, g[1], c_ls + meia, g[2]
        for _ in range(4):  # sem troca de sinal: a cerca dobra
            if glo > 0. >= ghi:
                break
            meia *= 2.
            lo, hi = c_ls - meia, c_ls + meia
            glo, ghi = g_lote([lo, hi], "x")
        cercada = glo > 0. >= ghi
        while cercada and hi - lo >= op["tol"]:
            m = 0.5 * (lo + hi)
            gm = g_lote([m], "m")[0]
            if gm > 0.:
                lo, glo = m, gm
            else:
                hi, ghi = m, gm
        c_bis = 0.5 * (lo + hi) if cercada else None
        sinais = ";".join(f"{c:g}:{'seg' if x > 0 else 'falha'}" for c, x in zip(extras, g_ext))
        esperado = all((x > 0.) == (c < c_ls) for c, x in zip(extras, g_ext))
        print(f"  linha {i:4d}: c_ls {c_ls:.5f}  c_bis {c_bis if c_bis is None else round(c_bis, 5)}  "
              f"dc {'-' if c_bis is None else f'{c_bis - c_ls:+.5f}'}  analises {n[0]}  conferencia "
              f"{'-' if conf is None else f'{conf:.1e}'}  extras {sinais}{'' if esperado else '  AVISO: sinal inesperado'}",
              flush=True)
        return {"linha": i, "c_ls": c_ls, "c_bis": c_bis, "dc": None if c_bis is None else c_bis - c_ls,
                "p_ls": cauda(c_ls), "p_bis": None if c_bis is None else cauda(c_bis), "lo": lo, "hi": hi,
                "g_lo": glo, "g_hi": ghi, "cercada": cercada, "analises": n[0], "c_conf": c0, "fator_ls": G0,
                "fator_conf": math.exp(g[0]), "conferencia": conf, "extras": sinais, "extras_ok": esperado}

    with cf.ThreadPoolExecutor(len(escolha)) as ex:
        res = sorted(ex.map(linha, escolha), key=lambda r: r["c_ls"])
    analises.shutdown()
    reg.close()
    cols = list(res[0])
    with open(os.path.join(pasta, "bissecao.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        w.writerows(res)
    # a P_f do LS e a mesma com as raizes da bissecao no lugar (o efeito das raizes conferidas na P_f)
    pf_ls = sum(cauda(c) for c in feitas.values()) / len(feitas)
    troca = {r["linha"]: r["c_bis"] for r in res if r["c_bis"] is not None}
    pf_bis = sum(cauda(troca.get(i, c)) for i, c in feitas.items()) / len(feitas)
    dcs = [abs(r["dc"]) for r in res if r["dc"] is not None]
    confs = [r["conferencia"] for r in res if r["conferencia"] is not None]
    resumo = {"pasta_ls": pasta_ls, "form": arq_form, "linhas": escolha, "tol": op["tol"], "extras": extras,
              "max_abs_dc": max(dcs) if dcs else None, "nao_cercadas": [r["linha"] for r in res if not r["cercada"]],
              "extras_inesperados": [r["linha"] for r in res if not r["extras_ok"]],
              "max_conferencia": max(confs) if confs else None, "P_f_ls": pf_ls, "P_f_com_bissecao": pf_bis,
              "analises": sum(r["analises"] for r in res), "resultado": res}
    json.dump(resumo, open(os.path.join(pasta, "bissecao.json"), "w"), indent=1)
    print(f"max |c_bis - c_ls| = {resumo['max_abs_dc']}; nao cercadas: {resumo['nao_cercadas']}; extras com sinal "
          f"inesperado: {resumo['extras_inesperados']}; conferencia max {resumo['max_conferencia']}", flush=True)
    print(f"P_f do LS {pf_ls:.5e}; com as raizes da bissecao {pf_bis:.5e} ({100 * (pf_bis / pf_ls - 1):+.2f}%); "
          f"{resumo['analises']} analises", flush=True)


if __name__ == "__main__":
    main()
