#!/usr/bin/env python3
"""teste_ls_cerca.py: testes da busca=cerca do line_sampling.py com executaveis falsos (roda em segundos).

    python3 teste_ls_cerca.py <pasta> [antigo=<line_sampling.py de antes, com o form_kl.py ao lado>] [n=60]
                              [tol=0.01] [paralelos=4] [ang=20] [beta=0.33] [semente=7]

(a) O padrao igual bit a bit (so com antigo=; ex.: git show 0614ddc:Projects/Cho2010/line_sampling.py e o form_kl.py
    numa pasta). O line_sampling.py antigo e o atual, na busca padrao (secante), com paralelos=1 (a ordem fixa), em 4
    casos:
    - lin: o coluna_falsa.py (estado limite linear, beta 2,2, direcao a 60 graus), 40 linhas;
    - cov: o mesmo, parando pelo CoV 5% (nmin 10);
    - origem: a origem na falha (beta -0,5 no form.json e no executavel);
    - degrau: o degrau_falso.py abaixo, que passa pelo teto de 10 analises e pela bissecao da secante.
    Compara byte a byte todos os arquivos da pasta e a saida do terminal. As 2 analises da inclinacao rodam em 2
    threads, e a ordem delas no analises.csv pode trocar entre duas rodadas (ja no codigo antigo); essa troca e
    relatada a parte.
(b) Raizes exatas: a busca=secante (padrao) e a busca=cerca nas mesmas n linhas, com dois executaveis falsos:
    - linear: o coluna_falsa.py, FS = exp(0,1 (beta - a.xi));
    - degrau: o degrau_falso.py (gravado na pasta), FS = 1,03 + 0,02 (beta - a.xi) antes de a.xi = beta e
      0,366875 - 0,001 (a.xi - beta) depois: um penhasco ate um patamar, como a coluna em 24 h.
    A direcao do LS fica a ang graus da verdadeira a, e a raiz exata da linha i e c_i = (beta - a.zp_i) / (a.e).
    Compara com as exatas: o erro maximo, as linhas com erro <= tol/2 e <= tol, as cercas (contem a raiz? largura <=
    tol?), as analises por linha e a P_f, contra a das raizes exatas (as mesmas linhas) e a exata Phi(-beta).
(c) Retomada: a cerca com n/2 linhas e depois com n (retomando) da as mesmas raizes e analises da cerca com n; a
    secante recusa a pasta da cerca, e a cerca recusa a da secante.
Saidas: <pasta>/<caso>/ (as rodadas), <pasta>/resumo.json e a tabela no terminal. Codigo de saida 1 se algum criterio
falhar.
"""
import csv
import json
import math
import os
import random
import subprocess
import sys

AQUI = os.path.dirname(os.path.abspath(__file__))
LS, FALSO = os.path.join(AQUI, "line_sampling.py"), os.path.join(AQUI, "coluna_falsa.py")
K = 100
DEGRAU = '''#!/usr/bin/python3 -IS
"""degrau_falso.py (gravado pelo teste_ls_cerca.py): a interface do coluna_falsa.py, com o FS em degrau em a.xi = beta:
1,03 + 0,02 (beta - a.xi) antes (seguro) e 0,366875 - 0,001 (a.xi - beta) depois (o patamar, falha)"""
import sys

op = {}
for a in sys.argv[1:]:
    k, v = a.split("=", 1)
    op[k] = v
x = [float(t) for t in open(op["xi"]).read().split()]
a = [float(t) for t in open(op["alvo"]).read().split()]
b = float(op["beta"])
t = sum(p * q for p, q in zip(a, x))
fs = 1.03 + 0.02 * (b - t) if t < b else max(0.366875 - 0.001 * (t - b), 0.3)
T = float(op["tempos"].split(",")[-1])
with open(op["saida"] + "_0_0.csv", "w") as f:
    f.write("t_h,situacao,t_passagem,FS_passagem,t_parou\\n")
    f.write(f"{T!r},ok,-1,{fs!r},-1\\n")
'''


def cauda(c):
    return 0.5 * math.erfc(c / math.sqrt(2.))


def unit(v):
    n = math.sqrt(sum(x * x for x in v))
    return [x / n for x in v]


def direcoes(ang, semente):
    """a (a verdadeira, nos modos 0-4) e e0 a ang graus dela"""
    rng = random.Random(semente)
    a = unit([0.25, -0.85, 0.10, -0.15, 0.30] + [0.] * (K - 5))
    w = [rng.gauss(0., 1.) for _ in range(K)]
    d = sum(x * y for x, y in zip(w, a))
    w = unit([x - d * y for x, y in zip(w, a)])
    t = math.radians(ang)
    return a, [math.cos(t) * x + math.sin(t) * y for x, y in zip(a, w)]


def grava_form(arq, alvo, beta, e0):
    """form.json como o do teste_ls_estagios.py, com u = beta e0 (beta < 0, a origem na falha: u = -|beta| e0, e o
    line_sampling.py inverte -e0 de volta para e0) e o beta tambem no executavel (falha em a.xi > beta)"""
    u = [beta * x for x in e0]
    json.dump({"parametros": {"prog": "coluna", "falha": "fs", "T": 24., "M": K, "threads": 0,
                              "exe": [f"alvo={alvo}", f"beta={beta!r}", "escala=0.1"]},
               "beta": beta, "u": u, "origem": "teste_ls_cerca.py"}, open(arq, "w"), indent=1)


def raizes_exatas(a, e0, beta, semente, n):
    """c_i = (beta - a.zp_i) / (a.e), com os Z e a direcao e como no line_sampling.py (e = e0 normalizada)"""
    e = unit(e0)
    rng = random.Random(semente)
    ae = sum(x * y for x, y in zip(a, e))
    out = []
    for _ in range(n):
        z = [rng.gauss(0., 1.) for _ in range(K)]
        ze = sum(x * y for x, y in zip(z, e))
        zp = [x - ze * y for x, y in zip(z, e)]
        out.append((beta - sum(x * y for x, y in zip(a, zp))) / ae)
    return out


def roda(script, exe, form, pasta, args, log):
    with open(log, "w") as f:
        return subprocess.run([sys.executable, "-I", script, exe, form, pasta] + args, stdout=f,
                              stderr=subprocess.STDOUT).returncode


def le_linhas(pasta):
    out = {}
    for r in csv.DictReader(open(os.path.join(pasta, "linhas.csv"))):
        out[int(r["linha"])] = (float(r["c"]), int(r["analises"]),
                                tuple(float(r[k]) if r.get(k) else None for k in ("seg", "fal")))
    return out


def compara(p1, p2, log1, log2):
    """'igual', 'igual a menos da ordem da inclinacao' ou a lista dos arquivos diferentes"""
    dif, so_ordem = [], False
    nomes = sorted(set(os.listdir(p1)) | set(os.listdir(p2)))
    pares = [(os.path.join(p1, x), os.path.join(p2, x), x) for x in nomes] + [(log1, log2, "log")]
    for f1, f2, x in pares:
        if not (os.path.exists(f1) and os.path.exists(f2)):
            dif.append(x + " (falta)")
            continue
        b1, b2 = open(f1, "rb").read(), open(f2, "rb").read()
        if b1 == b2:
            continue
        if x == "analises.csv":
            l1, l2 = b1.split(b"\n"), b2.split(b"\n")
            s1 = sorted(y for y in l1 if y.startswith(b"-1,"))
            s2 = sorted(y for y in l2 if y.startswith(b"-1,"))
            r1, r2 = [y for y in l1 if not y.startswith(b"-1,")], [y for y in l2 if not y.startswith(b"-1,")]
            if s1 == s2 and r1 == r2:
                so_ordem = True
                continue
        dif.append(x)
    if dif:
        return "DIFERENTE: " + ", ".join(dif)
    return "igual a menos da ordem das 2 analises da inclinacao" if so_ordem else "igual"


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    pasta = os.path.abspath(sys.argv[1])
    op = {"antigo": "", "n": 60, "tol": 0.01, "paralelos": 4, "ang": 20., "beta": 0.33, "semente": 7}
    for x in sys.argv[2:]:
        k, v = x.split("=", 1)
        op[k] = type(op[k])(v)
    if os.path.exists(pasta) and os.listdir(pasta):
        sys.exit(f"{pasta} ja existe e nao esta vazia (use uma pasta nova)")
    os.makedirs(pasta, exist_ok=True)
    degrau = os.path.join(pasta, "degrau_falso.py")
    open(degrau, "w").write(DEGRAU)
    os.chmod(degrau, 0o755)
    if not os.access(FALSO, os.X_OK):
        sys.exit(f"{FALSO} sem permissao de execucao (chmod +x)")
    alvo = os.path.join(pasta, "alvo.txt")
    a, e0 = direcoes(op["ang"], 1)
    open(alvo, "w").write("\n".join(repr(x) for x in a) + "\n")
    resumo, ok = {"op": op}, True

    # (a) o padrao igual bit a bit
    if op["antigo"]:
        antigo = os.path.abspath(op["antigo"])
        a60, e60 = direcoes(60., 2)
        alvo60 = os.path.join(pasta, "alvo60.txt")
        open(alvo60, "w").write("\n".join(repr(x) for x in a60) + "\n")
        f_lin, f_ori, f_deg = (os.path.join(pasta, f"form_{x}.json") for x in ("lin", "origem", "degrau"))
        grava_form(f_lin, alvo60, 2.2, e60)
        grava_form(f_ori, alvo, -0.5, e0)
        grava_form(f_deg, alvo, op["beta"], e0)
        comum = ["paralelos=1", f"K={K}"]
        casos = [("lin", FALSO, f_lin, ["n=40", "semente=1000", "cov=0"]),
                 ("cov", FALSO, f_lin, ["n=300", "semente=1001", "cov=0.05", "nmin=10"]),
                 ("origem", FALSO, f_ori, ["n=30", "semente=5", "cov=0"]),
                 ("degrau", degrau, f_deg, ["n=30", f"semente={op['semente']}", "cov=0"])]
        resumo["a"] = {}
        print("(a) busca padrao: o line_sampling.py antigo x o atual (paralelos=1)", flush=True)
        for nome, exe, form, args in casos:
            res = {}
            for versao, script in (("antigo", antigo), ("atual", LS)):
                pv = os.path.join(pasta, f"a_{nome}_{versao}")
                res[versao] = (pv, pv + ".log", roda(script, exe, form, pv, comum + args, pv + ".log"))
            cmp = compara(res["antigo"][0], res["atual"][0], res["antigo"][1], res["atual"][1])
            ls = json.load(open(os.path.join(res["atual"][0], "ls.json")))
            print(f"  {nome:7s} codigos {res['antigo'][2]}/{res['atual'][2]}  linhas {ls['linhas']:3d}  "
                  f"analises {ls['analises']:4d}  P_f {ls['P_f']:.6e}: {cmp}", flush=True)
            resumo["a"][nome] = cmp
            ok &= cmp.startswith("igual") and res["antigo"][2] == 0 and res["atual"][2] == 0

    # (b) raizes exatas: secante x cerca, linear e degrau
    print(f"(b) raizes exatas (n = {op['n']}, tol {op['tol']}, direcao a {op['ang']:g} graus, beta {op['beta']})",
          flush=True)
    exatas = raizes_exatas(a, e0, op["beta"], op["semente"], op["n"])
    p_ex = sum(cauda(c) for c in exatas) / len(exatas)
    form_b = os.path.join(pasta, "form_b.json")
    grava_form(form_b, alvo, op["beta"], e0)
    tol, comum = op["tol"], [f"n={op['n']}", f"semente={op['semente']}", f"paralelos={op['paralelos']}", f"K={K}",
                             "cov=0", f"tol={op['tol']}"]
    resumo["b"] = {"P_raizes_exatas": p_ex, "P_exata": cauda(op["beta"])}
    print("| caso | busca | analises/linha | max analises | max \\|erro\\| | erro <= tol/2 | erro <= tol | cerca com a "
          "raiz e <= tol | P_f | P_f das raizes exatas |", flush=True)
    print("|---|---|---|---|---|---|---|---|---|---|", flush=True)
    for caso, exe in (("linear", FALSO), ("degrau", degrau)):
        for busca in ("secante", "cerca"):
            pb = os.path.join(pasta, f"b_{caso}_{busca}")
            cod = roda(LS, exe, form_b, pb, comum + ([] if busca == "secante" else ["busca=cerca"]), pb + ".log")
            lin = le_linhas(pb)
            err = [abs(lin[i][0] - exatas[i]) for i in sorted(lin)]
            na = [lin[i][1] for i in sorted(lin)]
            boas = sum(1 for i in lin if lin[i][2][0] is not None and lin[i][2][1] is not None and
                       lin[i][2][0] <= exatas[i] + 1e-12 and exatas[i] <= lin[i][2][1] + 1e-12 and
                       lin[i][2][1] - lin[i][2][0] <= tol * (1. + 1e-9))
            ls = json.load(open(os.path.join(pb, "ls.json")))
            x = {"codigo": cod, "linhas": len(lin), "analises_linha": sum(na) / len(na), "max_analises": max(na),
                 "max_erro": max(err), "erro_meio_tol": sum(t <= 0.5 * tol + 1e-12 for t in err),
                 "erro_tol": sum(t <= tol + 1e-12 for t in err), "cercas_ok": boas, "P_f": ls["P_f"],
                 "analises_total": ls["analises"]}
            resumo["b"][f"{caso}_{busca}"] = x
            cerca_txt = f"{boas}/{len(lin)}" if busca == "cerca" else "-"
            print(f"| {caso} | {busca} | {x['analises_linha']:.2f} | {x['max_analises']} | {x['max_erro']:.4f} | "
                  f"{x['erro_meio_tol']}/{len(lin)} | {x['erro_tol']}/{len(lin)} | {cerca_txt} | {100 * x['P_f']:.3f}% | "
                  f"{100 * p_ex:.3f}% |", flush=True)
            if busca == "cerca":
                ok &= cod == 0 and x["erro_meio_tol"] == len(lin) == op["n"] and boas == len(lin)
    print(f"  P_f exata Phi(-beta) = {100 * cauda(op['beta']):.3f}%", flush=True)

    # (c) retomada e recusa
    pr = os.path.join(pasta, "c_retomada")
    args = [f"semente={op['semente']}", f"paralelos={op['paralelos']}", f"K={K}", "cov=0", f"tol={op['tol']}",
            "busca=cerca"]
    c1 = roda(LS, degrau, form_b, pr, args + [f"n={op['n'] // 2}"], pr + "_1.log")
    c2 = roda(LS, degrau, form_b, pr, args + [f"n={op['n']}"], pr + "_2.log")
    ref, ret = le_linhas(os.path.join(pasta, "b_degrau_cerca")), le_linhas(pr)
    mesma = ref == ret
    r1 = roda(LS, degrau, form_b, os.path.join(pasta, "b_degrau_cerca"), [f"n={op['n']}", f"K={K}"],
              pr + "_recusa_secante.log")
    r2 = roda(LS, degrau, form_b, os.path.join(pasta, "b_degrau_secante"), [f"n={op['n']}", f"K={K}", "busca=cerca"],
              pr + "_recusa_cerca.log")
    msg1, msg2 = open(pr + "_recusa_secante.log").read().strip(), open(pr + "_recusa_cerca.log").read().strip()
    resumo["c"] = {"codigos": [c1, c2], "retomada_igual": mesma, "recusa_secante": [r1, msg1],
                   "recusa_cerca": [r2, msg2]}
    print(f"(c) retomada da cerca ({op['n'] // 2} e depois {op['n']} linhas) = a cerca direta: {mesma} "
          f"(codigos {c1}/{c2})", flush=True)
    print(f"    secante na pasta da cerca: codigo {r1} ({msg1})", flush=True)
    print(f"    cerca na pasta da secante: codigo {r2} ({msg2})", flush=True)
    ok &= c1 == 0 and c2 == 0 and mesma and r1 != 0 and r2 != 0
    resumo["ok"] = ok
    json.dump(resumo, open(os.path.join(pasta, "resumo.json"), "w"), indent=1)
    print("TUDO CERTO" if ok else "FALHOU", flush=True)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
