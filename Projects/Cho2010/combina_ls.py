#!/usr/bin/env python3
"""combina_ls.py: junta as estimativas dos estagios de um line sampling em estagios (opcao C da nota
viabilidade_ls_adaptativo.md), com pesos fixados antes de rodar.

    python3 combina_ls.py <pasta>[:ini-fim] [<pasta>[:ini-fim] ...] [n=50,100 | w=0.25,0.75] [mc=<P_f>,<EP>]
                          [json=<saida.json>] [boot=2000]

<pasta>: uma rodada do line_sampling.py (linhas.csv, perdidas.csv, ls.json ou inclinacao.json); :ini-fim usa so as
linhas de indice ini..fim (inclusive), ex.: ls_p500_T18:0-49.

Estagio k: P_k = media dos p_i = Phi(-c_i) das N_k linhas, Var_k = s_k^2 / N_k (s_k o desvio de uma linha).
Combinada: P = sum w_k P_k, Var(P) = sum w_k^2 Var_k (estagios independentes), CoV = raiz(Var) / P, IC 95% normal.
Os pesos tem de ser fixados antes de rodar (nota, sec. 8): por padrao proporcionais aos numeros PLANEJADOS de linhas
(n=, na ordem das pastas; sem n=, o tamanho da faixa ini-fim, ou o n maximo da rodada no ls.json), nunca pelos
resultados (pesos pela variancia amostral dao vies para baixo, Papaioannou e Straub 2021, sec. 3.3). w= da os pesos
direto (normalizados para somar 1; w=0,1 descarta o estagio 1, quando a direcao inicial e sabidamente ruim).
Cada estagio usa Z independentes dos que fizeram a sua direcao e dos outros estagios: semente distinta, ou a mesma
semente com linhas disjuntas. Isso e conferido entre as pastas e, quando o form.json do estagio veio do direcao_ls.py,
contra as fontes da direcao (semente e linhas); se falhar, para com erro.
Analises: as das linhas mais as 2 da inclinacao de cada estagio (o FORM ou piloto nao entra; no estagio de direcao
reestimada, as analises que fizeram a direcao sao as do estagio anterior, ja contadas). Analises p/ CoV 5% = analises
x (CoV / 0,05)^2 (mesmas proporcoes entre estagios). mc=: a P_f e o EP do MC de referencia, para o z = (P - P_mc) /
raiz(Var + EP_mc^2) e o ganho (analises do MC bruto p/ CoV 5%, (1 - P) / (P 0,05^2), sobre as do LS).
boot=: reamostragens (bootstrap das linhas, semente fixa) para o IC 90% do CoV por linha de cada estagio (a cauda dos
p_i e pesada: com ~100-200 linhas o CoV por linha tem +-10-15% de erro amostral; decidir pelo IC, nao pelo valor);
boot=0 desliga.
Sem numpy (roda no fermi).
"""
import csv
import json
import math
import os
import random
import re
import sys


def cauda(c):
    """Phi(-c), como no line_sampling.py"""
    return 0.5 * math.erfc(c / math.sqrt(2.))


def le_estagio(arg):
    m = re.match(r"^(.*):(\d+)-(\d+)$", arg)
    pasta, faixa = (m.group(1), (int(m.group(2)), int(m.group(3)))) if m else (arg, None)
    pasta = os.path.abspath(pasta)
    dentro = (lambda i: faixa[0] <= i <= faixa[1]) if faixa else (lambda i: True)
    arq_ls = os.path.join(pasta, "ls.json")
    ls = json.load(open(arq_ls)) if os.path.exists(arq_ls) else {}
    semente = ls.get("semente")
    if semente is None:
        sys.exit(f"{pasta}: sem ls.json (rodada em andamento?): sem a semente")
    linhas = {}
    for r in csv.DictReader(open(os.path.join(pasta, "linhas.csv"))):
        i = int(r["linha"])
        if dentro(i):
            linhas[i] = (float(r["c"]), int(r["analises"]))
    perdidas = []
    if os.path.exists(os.path.join(pasta, "perdidas.csv")):
        perdidas = [int(r["linha"]) for r in csv.DictReader(open(os.path.join(pasta, "perdidas.csv")))
                    if dentro(int(r["linha"]))]
    if len(linhas) < 2:
        sys.exit(f"{pasta}: menos de 2 linhas")
    pi = [cauda(linhas[i][0]) for i in sorted(linhas)]
    N = len(pi)
    P = sum(pi) / N
    sd = math.sqrt(sum((x - P) ** 2 for x in pi) / (N - 1))
    na = sum(v[1] for v in linhas.values()) + 2
    planejado = (faixa[1] - faixa[0] + 1) if faixa else ls.get("criterio", {}).get("n_max")
    form = ls.get("form")
    fontes = []
    if form and os.path.exists(form):
        fontes = json.load(open(form)).get("direcao_ls", {}).get("fontes", [])
    return {"nome": os.path.basename(pasta) + (f":{faixa[0]}-{faixa[1]}" if faixa else ""), "pasta": pasta,
            "semente": semente, "linhas": linhas, "perdidas": perdidas, "N": N, "P": P, "sd": sd,
            "var": sd * sd / N, "pi": pi, "analises": na, "planejado": planejado, "form": form, "fontes": fontes}


def main():
    args, op = [], {}
    for a in sys.argv[1:]:
        m = re.match(r"^(n|w|mc|json|boot)=(.*)$", a)
        if m:
            op[m.group(1)] = m.group(2)
        else:
            args.append(a)
    if not args:
        sys.exit(__doc__)
    est = [le_estagio(a) for a in args]
    # Z independentes: entre estagios e contra as linhas que fizeram a direcao de cada estagio
    for a in range(len(est)):
        for b in range(a + 1, len(est)):
            if est[a]["semente"] == est[b]["semente"] and set(est[a]["linhas"]) & set(est[b]["linhas"]):
                sys.exit(f"{est[a]['nome']} e {est[b]['nome']}: mesma semente e linhas em comum (os mesmos Z)")
        for f in est[a]["fontes"]:
            if f["semente"] == est[a]["semente"] and set(f["linhas"]) & set(est[a]["linhas"]):
                sys.exit(f"{est[a]['nome']}: a direcao (direcao_ls.py) usou os Z de linhas deste estagio "
                         f"({os.path.basename(f['pasta'])}, semente {f['semente']}): vies")
    if "w" in op:
        w = [float(x) for x in op["w"].split(",")]
        regra = "w= dado"
    else:
        n = [float(x) for x in op["n"].split(",")] if "n" in op else [e["planejado"] for e in est]
        if any(x is None for x in n):
            sys.exit("sem o numero planejado de linhas de algum estagio: de n= ou w=")
        w = n
        regra = "proporcionais ao numero planejado de linhas (" + ("n= dado" if "n" in op else "faixa ou n_max") + \
            ": " + ", ".join(f"{x:g}" for x in n) + ")"
    if len(w) != len(est) or sum(w) <= 0. or min(w) < 0.:
        sys.exit("pesos: um por estagio, >= 0, soma > 0")
    w = [x / sum(w) for x in w]
    mc = [float(x) for x in op["mc"].split(",")] if "mc" in op else None
    nboot = int(op.get("boot", 2000))

    def ic_boot(pi):
        """IC 90% (percentis 5 e 95) do CoV por linha, bootstrap das linhas"""
        if nboot <= 0 or len(pi) < 3:
            return None
        rb, v = random.Random(1), []
        for _ in range(nboot):
            x = [pi[rb.randrange(len(pi))] for _ in pi]
            m = sum(x) / len(x)
            if m > 0.:
                v.append(math.sqrt(sum((t - m) ** 2 for t in x) / (len(x) - 1)) / m)
        v.sort()
        return [v[int(0.05 * len(v))], v[min(int(0.95 * len(v)), len(v) - 1)]] if v else None

    def linha_tab(nome, wk, P, var, N, na, sd=None):
        ep = math.sqrt(var)
        cov = ep / P if P > 0 else math.inf
        a5 = na * (cov / 0.05) ** 2
        cl = f"{sd / P:6.2f}" if sd is not None and P > 0 else "     -"
        z = f"{(P - mc[0]) / math.sqrt(var + mc[1] ** 2):+5.2f}" if mc else "    -"
        ganho = f"{(1. - mc[0]) / (mc[0] * 0.05 ** 2) / a5:6.1f}" if mc and a5 > 0 else "     -"
        print(f"| {nome:28s} | {wk:5.3f} | {N:5d} | {100 * P:6.3f} | [{100 * (P - 1.96 * ep):6.3f}; "
              f"{100 * (P + 1.96 * ep):6.3f}] | {100 * cov:6.2f} | {cl} | {na:6d} | {na / N:5.2f} | {a5:8.0f} | {z} | "
              f"{ganho} |", flush=True)
        return {"P_f": P, "EP": ep, "IC95": [P - 1.96 * ep, P + 1.96 * ep], "CoV": cov, "linhas": N, "analises": na,
                "analises_cov5": a5}

    print(f"pesos {regra}: " + ", ".join(f"{x:.4f}" for x in w), flush=True)
    print("| estagio | peso | linhas | P_f (%) | IC 95% (%) | CoV (%) | CoV/linha | analises | an./linha | "
          "an. p/ CoV 5% | z vs MC | ganho vs MC |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|")
    saida = {"pesos": w, "regra": regra, "estagios": []}
    for e, wk in zip(est, w):
        r = linha_tab(e["nome"] + f" (s{e['semente']})", wk, e["P"], e["var"], e["N"], e["analises"], e["sd"])
        r.update({"pasta": e["pasta"], "semente": e["semente"], "peso": wk, "CoV_linha": e["sd"] / e["P"]
                  if e["P"] > 0 else None, "perdidas": len(e["perdidas"]), "form": e["form"]})
        ic = ic_boot(e["pi"])
        if ic:
            r["CoV_linha_IC90_boot"] = ic
            print(f"  {e['nome']}: CoV por linha {e['sd'] / e['P']:.2f}, IC 90% bootstrap [{ic[0]:.2f}; {ic[1]:.2f}] "
                  f"({nboot} reamostragens)", flush=True)
        saida["estagios"].append(r)
        if e["perdidas"]:
            print(f"  {e['nome']}: {len(e['perdidas'])} linhas perdidas na faixa (fora da conta)", flush=True)
    P = sum(wk * e["P"] for e, wk in zip(est, w))
    var = sum(wk * wk * e["var"] for e, wk in zip(est, w))
    saida["combinada"] = linha_tab("combinada", 1., P, var, sum(e["N"] for e in est),
                                   sum(e["analises"] for e in est))
    cl = [e["sd"] / e["P"] for e, wk in zip(est, w) if wk > 0. and e["P"] > 0.]
    if len(cl) > 1 and max(cl) > 2. * min(cl):
        print(f"AVISO: CoV por linha de {min(cl):.2f} a {max(cl):.2f} entre os estagios com peso > 0. O estagio de cauda "
              f"pesada (direcao ruim) infla a variancia da combinada, e o EP dele, com poucas linhas, sai baixo junto com "
              f"a P_f (T0, partida a 60 graus, 200 repeticoes: |P - exata| < 2 EP em 81% com 50 linhas no estagio 1 e 85% "
              f"com 150, todas as falhas para baixo; so o estagio 2: 95% e 97%, com desvio real 1,8x e 3,3x menor). Se a "
              f"direcao do estagio 1 e sabidamente ruim, a regra fixada antes e w=0,1 (so o estagio 2); nunca pesos "
              f"tirados destes CoV (vies).", flush=True)
    if mc:
        saida["mc"] = {"P_f": mc[0], "EP": mc[1]}
        print(f"MC de referencia: P_f = {100 * mc[0]:.3f}% (EP {100 * mc[1]:.3f}%); MC bruto p/ CoV 5%: "
              f"{(1. - mc[0]) / (mc[0] * 0.05 ** 2):.0f} analises", flush=True)
    if "json" in op:
        json.dump(saida, open(op["json"], "w"), indent=1)


if __name__ == "__main__":
    main()
