#!/usr/bin/env python3
"""teste_ls_reuso.py: testes do reuso do transiente (reuso=1) do form_kl.py e do line_sampling.py na chuva
(prog=chuva), com um executavel falso (roda em segundos).

    python3 teste_ls_reuso.py <pasta> antigo=<pasta com o line_sampling.py e o form_kl.py de antes> [n=30] [M=10]

O executavel falso (chuva_falsa.py, gravado na pasta) tem a interface do SlopeRainMonteCarlo que os dois usam: le os
xi (3K: c', phi', k), grava <saida>_0_0.csv com uma linha por instante de tempos= (FS = exp(0,1 (beta - a.xi)) vezes
1 + 0,05 |t - 24| / 24: o menor em 24 h; o k nao entra), o _grad_0.txt com gradiente= (0 nos termos de k) e, como o
patch do campo de k, grava_pressoes= (a chave com o FNV-1a dos bits dos xi de k) e pressoes= (confere a chave; outra
e erro). Com velho=1, recusa grava_pressoes=/pressoes= (o executavel sem o patch). Cada analise acrescenta
"<saida> <modo> <campo de k>" ao registro= (modo calcula, grava ou le).
(a) O padrao igual bit a bit: o form_kl.py (gradiente=modal, 4 instantes) e o line_sampling.py (busca=secante e
    busca=cerca, paralelos=1) de antes e os de agora, sem reuso: todos os arquivos das pastas e a saida do terminal
    byte a byte (o analises.csv da secante tambem ordenado: as 2 analises da inclinacao rodam em 2 threads).
(b) O reuso: o FORM e o LS (cerca, paralelos=3) com reuso=1 contra os mesmos sem reuso: o mesmo u* e as mesmas raizes,
    analises e cercas por linha; no registro, por linha, a 1a analise grava, as outras leem, o campo de k e um so na
    linha e muda entre linhas; nenhum l<i>_pw.txt sobra; o ls.json tem "reuso": 1; o custos.csv tem os modos.
(c) As recusas, sem nenhuma analise de linha: reuso=1 com a direcao com k (um form.json com u_k != 0), com o
    executavel velho (velho=1: a sonda), com pressoes= nas opcoes; o form_kl.py com reuso=1 e gradiente=dif.
(d) Retomada com reuso: n/2 linhas e depois n (retomando) da o mesmo que n de uma vez.
Codigo de saida 1 se algum criterio falhar.
"""
import csv
import filecmp
import json
import os
import shutil
import subprocess
import sys

AQUI = os.path.dirname(os.path.abspath(__file__))
FALSO = r'''#!/usr/bin/python3 -IS
"""chuva_falsa.py (gravado pelo teste_ls_reuso.py): a interface do SlopeRainMonteCarlo para o form_kl.py e o
line_sampling.py (prog=chuva), com o grava_pressoes=/pressoes= do patch do campo de k"""
import math
import os
import struct
import sys

op = {}
for a in sys.argv[1:]:
    k, v = a.split("=", 1)
    op[k] = v
x = [float(t) for t in open(op["xi"]).read().split()]
K = len(x) // 3
a = [float(t) for t in open(op["alvo"]).read().split()]
h = 1469598103934665603
for v in x[2 * K:]:
    for b in struct.pack("<d", 0. if v == 0. else v):
        h = ((h ^ b) * 1099511628211) & 0xFFFFFFFFFFFFFFFF
chave = f"falso;campo_k={K},{h:x}"
gp, pr = op.get("grava_pressoes"), op.get("pressoes")
if gp and pr:
    sys.exit("pressoes e grava_pressoes: um ou outro")
if op.get("velho") == "1" and (gp or pr):
    sys.exit("pressoes/grava_pressoes so com k fixo (covk=0 e, com valores, k = 1)")
modo = "calcula"
if pr:
    with open(pr) as f:
        f.readline()
        if f.readline().strip() != "chave " + chave:
            sys.exit("pressoes: o transiente de " + pr + " e outro")
    modo = "le"
elif gp:
    with open(gp + ".tmp", "w") as f:
        f.write("# falso\nchave " + chave + "\n")
    os.replace(gp + ".tmp", gp)
    modo = "grava"
if op.get("registro"):
    with open(op["registro"], "a") as f:
        f.write(f"{os.path.basename(op['saida'])} {modo} {h:x}\n")
fs0 = math.exp(float(op.get("escala", 0.1)) * (float(op.get("beta", 2.)) - sum(p * q for p, q in zip(a, x))))
with open(op["saida"] + "_0_0.csv", "w") as f:
    f.write("real,t_h,FS_min,situacao,s_fator,s_transiente\n")
    for t in [float(s) for s in op.get("tempos", "0,24").split(",")]:
        f.write(f"0,{t!r},{fs0 * (1 + 0.05 * abs(t - 24) / 24)!r},ok,0.5,{0. if modo == 'le' else 1.}\n")
if "gradiente" in op:
    G = int(op["gradiente"])
    g = [-float(op.get("escala", 0.1)) * fs0 * ai for ai in a[:2 * G]] + [0.] * G
    with open(op["saida"] + "_grad_0.txt", "w") as f:
        f.write("# falso\n" + repr(fs0) + "\n" + "\n".join(repr(v) for v in g) + "\n")
'''


def roda(cmd, log):
    r = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    with open(log, "w") as f:
        f.write(r.stdout)
    return r.returncode, r.stdout


def arquivos(p):
    return sorted(x for x in os.listdir(p) if os.path.isfile(os.path.join(p, x)))


def linhas_ls(p):
    with open(os.path.join(p, "linhas.csv"), newline="") as f:
        return {int(r["linha"]): (r["c"], r["analises"], r.get("seg"), r.get("fal")) for r in csv.DictReader(f)}


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    base = os.path.abspath(sys.argv[1])
    op = {"antigo": "", "n": 30, "M": 10}
    for a in sys.argv[2:]:
        k, v = a.split("=", 1)
        op[k] = type(op[k])(v)
    if not op["antigo"]:
        sys.exit("antigo=<pasta com o line_sampling.py e o form_kl.py de antes> (ex.: git show HEAD:...)")
    antigo = os.path.abspath(op["antigo"])
    n, M = op["n"], op["M"]
    shutil.rmtree(base, ignore_errors=True)
    os.makedirs(base)
    exe = os.path.join(base, "chuva_falsa.py")
    with open(exe, "w") as f:
        f.write(FALSO)
    os.chmod(exe, 0o755)
    # a direcao verdadeira: em c' e phi', com |a| = 1 (os K de k em 0); beta 2
    alvo = os.path.join(base, "alvo.txt")
    a = [1. / (i + 1) for i in range(M)] + [0.5 / (i + 1) for i in range(M)]
    na = sum(x * x for x in a) ** 0.5
    with open(alvo, "w") as f:
        f.write("\n".join(repr(x / na) for x in a + [0.] * M) + "\n")
    reg = os.path.join(base, "registro.txt")
    comum = [f"alvo={alvo}", "beta=2", "escala=0.1", "tempos=0,18,24,30", f"registro={reg}"]
    ok = True

    def confere(cond, msg):
        nonlocal ok
        print(("  ok    " if cond else "  FALHA ") + msg, flush=True)
        ok = ok and cond

    # (a) o padrao igual bit a bit
    print("(a) o padrao (sem reuso) igual bit a bit ao de antes")
    form = {}
    for qual, pasta_py in (("antigo", antigo), ("novo", AQUI)):
        d = os.path.join(base, "a_" + qual)
        os.makedirs(d)
        rc, out = roda(["python3", "-I", os.path.join(pasta_py, "form_kl.py"), exe, os.path.join(d, "form"),
                        "prog=chuva", "gradiente=modal", f"M={M}", "iteracoes=6", "paralelos=1"] + comum,
                       os.path.join(d, "form.log"))
        confere(rc == 0, f"form_kl.py {qual} rodou")
        form[qual] = os.path.join(d, "form", "form.json")
        for busca in ("secante", "cerca"):
            rc, out = roda(["python3", "-I", os.path.join(pasta_py, "line_sampling.py"), exe, form["antigo"],
                            os.path.join(d, "ls_" + busca), f"n={n}", f"K={M}", "paralelos=1", f"busca={busca}"],
                           os.path.join(d, f"ls_{busca}.log"))
            confere(rc == 0, f"line_sampling.py {qual} busca={busca} rodou")
    da, dn = os.path.join(base, "a_antigo"), os.path.join(base, "a_novo")
    confere(filecmp.cmp(form["antigo"], form["novo"], shallow=False), "form.json byte a byte")
    sa, sn = open(os.path.join(da, "form.log")).read(), open(os.path.join(dn, "form.log")).read()
    confere(sa.replace(da, "") == sn.replace(dn, ""), "saida do form_kl.py igual (menos a pasta)")
    for busca in ("secante", "cerca"):
        pa, pn = os.path.join(da, "ls_" + busca), os.path.join(dn, "ls_" + busca)
        confere(arquivos(pa) == arquivos(pn), f"{busca}: os mesmos arquivos ({', '.join(arquivos(pn))})")
        for x in arquivos(pn):
            if filecmp.cmp(os.path.join(pa, x), os.path.join(pn, x), shallow=False):
                confere(True, f"{busca}: {x} byte a byte")
            elif x == "analises.csv" and busca == "secante":
                confere(sorted(open(os.path.join(pa, x))) == sorted(open(os.path.join(pn, x))),
                        f"{busca}: {x} igual fora da ordem das 2 analises da inclinacao (2 threads)")
            else:
                confere(False, f"{busca}: {x} byte a byte")
        la = open(os.path.join(da, f"ls_{busca}.log")).read().replace(pa, "")
        ln = open(os.path.join(dn, f"ls_{busca}.log")).read().replace(pn, "")
        confere(la == ln, f"{busca}: saida do terminal igual (menos a pasta)")

    # (b) o reuso contra o mesmo sem reuso
    print("(b) reuso=1 contra reuso=0")
    db = os.path.join(base, "b")
    os.makedirs(db)
    open(reg, "w").close()
    rc, _ = roda(["python3", "-I", os.path.join(AQUI, "form_kl.py"), exe, os.path.join(db, "form_reuso"), "prog=chuva",
                  "gradiente=modal", f"M={M}", "iteracoes=6", "paralelos=1", "reuso=1"] + comum,
                 os.path.join(db, "form_reuso.log"))
    confere(rc == 0, "form_kl.py reuso=1 rodou")
    fr = json.load(open(os.path.join(db, "form_reuso", "form.json")))
    f0 = json.load(open(form["novo"]))
    confere(fr["u"] == f0["u"] and fr["beta"] == f0["beta"], f"FORM: o mesmo u* e beta ({fr['beta']:.6f})")
    confere(fr["parametros"].get("reuso") == 1, "form.json com reuso 1 nos parametros")
    rg = [x.split() for x in open(reg)]
    modos = [m for r, m, _ in rg if r.startswith("it")]
    confere([r for r, _, _ in rg][0] == "sonda_reuso", "a sonda roda antes do FORM")
    confere(modos[0] == "grava" and all(m == "le" for m in modos[1:]),
            f"FORM: a 1a analise grava e as {len(modos) - 1} outras leem")
    confere(len({h for r, _, h in rg if r.startswith("it")}) == 1, "FORM: um campo de k so (xi_k = 0)")
    res = {}
    for reuso in (0, 1):
        p = os.path.join(db, f"ls_r{reuso}")
        open(reg, "w").close()
        rc, _ = roda(["python3", "-I", os.path.join(AQUI, "line_sampling.py"), exe, form["novo"], p, f"n={n}", f"K={M}",
                      "paralelos=3", "busca=cerca", f"reuso={reuso}"], p + ".log")
        confere(rc == 0, f"line_sampling.py reuso={reuso} rodou")
        res[reuso] = (linhas_ls(p), [x.split() for x in open(reg)], p)
    confere(res[0][0] == res[1][0], f"LS: as mesmas raizes, analises e cercas nas {len(res[1][0])} linhas")
    lj = json.load(open(os.path.join(res[1][2], "ls.json")))
    l0 = json.load(open(os.path.join(res[0][2], "ls.json")))
    confere(lj.get("reuso") == 1 and "reuso" not in l0, "ls.json: reuso 1 so com reuso=1")
    confere(lj["P_f"] == l0["P_f"] and lj["analises"] == l0["analises"], f"LS: a mesma P_f ({lj['P_f']:.6e})")
    por = {}
    for r, m, h in res[1][1]:
        if r.startswith("l"):
            i, j = (int(v) for v in r[1:].split("_"))
            por.setdefault(i, {})[j] = (m, h)
    bom = all(v[0][0] == "grava" and all(v[j][0] == "le" for j in v if j) and len({x[1] for x in v.values()}) == 1
              for v in por.values())
    confere(bom and len(por) == len(res[1][0]), "LS: em cada linha a 1a analise grava, as outras leem, um campo de k")
    confere(len({v[0][1] for v in por.values()}) == len(por), "LS: um campo de k diferente por linha")
    confere(all(m == "calcula" for r, m, _ in res[0][1] if r.startswith("l")), "LS sem reuso: todas calculam")
    sobra = [x for x in os.listdir(res[1][2]) if x.endswith(("_pw.txt", "_pw.txt.tmp", ".xi"))]
    confere(not sobra, f"LS: nenhum arquivo do transiente sobrou ({sobra})")
    with open(os.path.join(res[1][2], "custos.csv"), newline="") as f:
        cs = list(csv.DictReader(f))
    confere(len(cs) == lj["analises"] and {r["modo"] for r in cs} == {"grava", "le"},
            f"custos.csv: {len(cs)} analises, modos grava e le")
    confere(not os.path.exists(os.path.join(res[0][2], "custos.csv")), "sem reuso: sem custos.csv (o de antes)")

    # (c) as recusas
    print("(c) recusas")
    dc = os.path.join(base, "c")
    os.makedirs(dc)
    fk = json.load(open(form["novo"]))
    fk["u"] = list(fk["u"])
    fk["u"][2 * M] = 0.1  # um termo de k
    arq_fk = os.path.join(dc, "form_k.json")
    json.dump(fk, open(arq_fk, "w"))
    open(reg, "w").close()
    casos = [("direcao com k", [exe, arq_fk, os.path.join(dc, "ls_k"), "reuso=1"], "termos de k"),
             ("executavel velho", [exe, form["novo"], os.path.join(dc, "ls_velho"), "reuso=1", "velho=1"],
              "recusou grava_pressoes"),
             ("pressoes= nas opcoes", [exe, form["novo"], os.path.join(dc, "ls_pr"), "reuso=1", "pressoes=x.txt"],
              "sem pressoes=")]
    for nome, args, frase in casos:
        rc, out = roda(["python3", "-I", os.path.join(AQUI, "line_sampling.py")] + args + [f"n={n}", f"K={M}",
                                                                                          "busca=cerca"],
                       os.path.join(dc, nome.replace(" ", "_") + ".log"))
        confere(rc != 0 and frase in out, f"LS recusa: {nome}")
    rc, out = roda(["python3", "-I", os.path.join(AQUI, "form_kl.py"), exe, os.path.join(dc, "form_dif"), "prog=chuva",
                    "gradiente=dif", f"M={M}", "reuso=1"] + comum, os.path.join(dc, "form_dif.log"))
    confere(rc != 0 and "so com gradiente=modal" in out, "FORM recusa: reuso=1 com gradiente=dif")
    linhas_rodadas = [r for r in open(reg) if r.startswith("l")]
    confere(not linhas_rodadas, "nenhuma analise de linha nas recusas")

    # (d) retomada com reuso
    print("(d) retomada com reuso")
    pd = os.path.join(base, "d")
    rc1, _ = roda(["python3", "-I", os.path.join(AQUI, "line_sampling.py"), exe, form["novo"], pd, f"n={n // 2}",
                   f"K={M}", "paralelos=3", "busca=cerca", "reuso=1"], pd + "_1.log")
    rc2, _ = roda(["python3", "-I", os.path.join(AQUI, "line_sampling.py"), exe, form["novo"], pd, f"n={n}",
                   f"K={M}", "paralelos=3", "busca=cerca", "reuso=1"], pd + "_2.log")
    confere(rc1 == 0 and rc2 == 0 and linhas_ls(pd) == res[1][0], "n/2 e depois n (retomando) = n de uma vez")
    print("TUDO OK" if ok else "HA FALHAS")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
