#!/usr/bin/env python3
"""line_sampling.py: line sampling (amostragem por linhas) na direcao do FORM, com o campo completo.

    python3 line_sampling.py <exe> <form.json> <pasta> [n=100] [semente=7] [paralelos=14] [tol=0.01] [K=2000]
                             [cov=0] [nmin=30]
                             [chave=valor das opcoes do programa: nystrom=, metodo=, srm_f0=, contorno=, tfim=, ...]

A direcao e = u*/beta vem do FORM (form_kl.py, M termos por campo: c', phi' no Cho2010; c', phi', k no rebaixamento e na
chuva),
completada com zeros ate os K termos de cada campo. As analises usam o campo com K termos (M=K) e as opcoes do FORM,
trocadas pelas chave=valor dadas (ex.: nystrom=<a base completa>, metodo=srm). No rebaixamento, as opcoes contorno=,
agua=, trebaix=, tfim=, dt=, hw=, covc=, covphi=, covk=, s= e alfa= vem do form.json (form.json antigos: os padroes do
SlopeDrawdownMonteCarlo). Na chuva (prog=chuva), as opcoes do executavel vem do form.json (a lista "exe") e as
chave=valor dadas aqui trocam ou acrescentam opcoes nela (ex.: nystrom=<outra base>, pressoes=<outro arquivo>); o fator
de cada analise e o menor FS dos instantes (form_kl.le_fator).
Cada linha: sorteia z ~ N(0, I) com nc K termos, tira a componente em e (zp = z - (z.e) e) e acha, na reta zp + c e, o
c_i em que o fator (FS ou Gamma) vale 1. Como zp e c sao independentes, P_f = E[Phi(-c_i)], e cada linha da
p_i = Phi(-c_i): so a variacao de c_i entre as linhas (a nao linearidade da superficie de falha) entra na variancia. O
estimador nao tem vies para nenhuma direcao; uma direcao ruim so aumenta a dispersao.

A raiz e pela secante em g(c) = ln(fator), quase linear em c: a 1a analise em c = beta, a 2a no passo de Newton com a
inclinacao s medida no ponto de projeto (2 analises por direcao), e a raiz da secante e aceita sem nova analise quando o
passo for menor que tol (em c). Os mesmos z para qualquer direcao com a mesma semente (numeros aleatorios comuns).
Sem numpy (roda no fermi).

Parada pelo CoV (cov=0.05, por exemplo): as linhas vao sendo submetidas aos poucos (no maximo 'paralelos' em voo,
sempre a menor linha que falta) e, a cada linha concluida, se ja ha nmin linhas e o CoV(P_f) estimado e <= cov, nao
se submetem mais linhas; as que estao em voo terminam e entram na conta (o CoV final e um pouco menor que o alvo). Sem
cov= (ou cov=0) roda as n linhas. n e sempre o maximo. A linha i usa sempre o mesmo Z[i] (o i-esimo vetor da semente),
entao o resultado nao depende da ordem de conclusao e a retomada continua do mesmo ponto. Um arquivo PARAR na pasta
tambem para de submeter (motivo "PARAR").
Vies da regra de parada: o CoV estimado, sd / (sqrt(N) P_f), e menor quando a P_f da amostra esta acima da media, entao
parar no primeiro N com CoV <= alvo favorece amostras com P_f um pouco alta. O vies relativo e da ordem do quadrado do
CoV (para CoV 5%, uns 0,1 a 0,3%), desprezivel diante do proprio CoV; o nmin evita parar cedo com uma variancia
estimada com poucas linhas (o caso em que o vies seria maior).

Origem na falha (g(0) < 0, beta negativo no form.json): a direcao e -u*/|beta|, o ponto de projeto fica em c = beta < 0
e a falha continua do lado de c grande, com a mesma formula.
Linhas perdidas: uma analise com situacao diferente de ok tira a linha da conta (o transiente que nao convergiu, ou
Gamma NaN). A excecao e Gamma = 0 (situacao "fator" com FS_min = 0: o GIM nao achou equilibrio nem com 1% da gravidade),
que e falha de verdade: entra como Gamma = GAMMA_ZERO (0,005, abaixo de 0,01) e a linha segue. As perdidas vao para o
perdidas.csv (linha, rotulo, motivo), nao sao refeitas na retomada (o mesmo Z daria a mesma falha), e o ls.json soma
as de todas as rodadas (perdidas, fracao_perdida, linhas_perdidas). Se elas ficarem nos campos fracos, a P_f sai baixa:
olhe a fracao (o aviso sai quando passa de 1%).

Saidas na pasta: linhas.csv (linha, c, analises, p), analises.csv (linha, c, fator), perdidas.csv e ls.json (o
resultado, com o criterio de parada, o motivo e as linhas feitas; com menos de 2 linhas, P_f = null), com retomada (as
linhas ja no linhas.csv ou no perdidas.csv nao sao refeitas; analises de uma linha interrompida ficam repetidas no
analises.csv); a tabela parcial (P_f, CoV, linhas para CoV de 5%) sai no terminal a cada 10 linhas.
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
from form_kl import REBAIX_PADRAO, ambiente, campos, comando, le_fator  # noqa: E402

FIM = (-8., 12.)
GAMMA_ZERO = 0.005  # Gamma = 0 (sem equilibrio nem com 1% da gravidade): falha, com ln finito para a secante


def limpa_csv(arq, cabecalho):
    """cria o csv com o cabecalho, ou descarta a ultima linha incompleta (morto no meio de uma escrita)"""
    if not os.path.exists(arq):
        open(arq, "w").write(cabecalho)
        return
    texto = open(arq).read()
    if texto and not texto.endswith("\n"):
        open(arq, "w").write(texto[:texto.rfind("\n") + 1])


def cauda(c):
    """Phi(-c)"""
    return 0.5 * math.erfc(c / math.sqrt(2.))


def main():
    if len(sys.argv) < 4:
        sys.exit(__doc__)
    exe, arq_form, pasta = os.path.abspath(sys.argv[1]), sys.argv[2], os.path.abspath(sys.argv[3])
    op = {"n": 100, "semente": 7, "paralelos": 14, "tol": 0.01, "K": 2000, "cov": 0., "nmin": 30}
    form = json.load(open(arq_form))
    p = dict(form["parametros"])
    p.setdefault("prog", "cho")          # form.json de versoes antigas do form_kl.py (antes de prog=)
    p.setdefault("transiente", "tri3")
    for k, v in REBAIX_PADRAO.items():  # form.json de antes das opcoes contorno=, agua=, ...
        p.setdefault(k, v)
    for a in sys.argv[4:]:
        k, v = a.split("=", 1)
        if k in op:
            op[k] = type(op[k])(v)
        elif p["prog"] == "chuva" and k not in ("threads",):  # uma opcao do executavel: troca ou acrescenta
            p["exe"] = [x for x in p["exe"] if not x.startswith(k + "=")] + [a]
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
    if form.get("beta", 1.) < 0.:  # g(0) < 0 (as medias ja falham): com -e a falha fica do lado de c grande e o ponto
        e = [-x for x in e]        # de projeto em c = beta < 0
        beta = -beta
    trava = threading.Lock()
    limpa_csv(os.path.join(pasta, "analises.csv"), "")
    reg = open(os.path.join(pasta, "analises.csv"), "a")

    def fator(x, rotulo, linha, c):
        arq = os.path.join(pasta, rotulo + ".xi")
        with open(arq, "w") as f:
            f.write("\n".join(repr(t) for t in x) + "\n")
        subprocess.run(comando(exe, p, arq, rotulo), check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                       env=ambiente(p), cwd=pasta)
        try:
            if p["prog"] == "chuva":  # o menor FS dos instantes (um instante sem equilibrio e falha)
                G = le_fator(os.path.join(pasta, rotulo + "_0_0.csv"), p)[0]
            else:
                with open(os.path.join(pasta, rotulo + "_0_0.csv"), newline="") as f:
                    r = next(csv.DictReader(f))
                G = float(r["FS_min"])
                if r["situacao"] == "fator" and G == 0.:  # GIM sem equilibrio nem com lambda < 0,01: falha
                    G = GAMMA_ZERO
                elif r["situacao"] != "ok":
                    raise RuntimeError(f"{rotulo}: situacao {r['situacao']} (FS_min {r['FS_min']})")
        finally:
            for suf in (".xi", "_0_0.csv", "_mapa_0_0.csv", "_parametros.txt"):
                try:
                    os.remove(os.path.join(pasta, rotulo + suf))
                except FileNotFoundError:
                    pass
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

    if p["prog"] in ("rebaixamento", "chuva"):
        print("  comando: " + " ".join(comando(exe, p, "<x>.xi", "<rotulo>")), flush=True)
    rng = random.Random(op["semente"])
    proximo_z = [0]  # Z[i]: o i-esimo vetor da semente, gerado em ordem (so os das linhas em voo ficam na memoria)

    def z_da_linha(i):
        while proximo_z[0] < i:
            for _ in range(dim):
                rng.gauss(0., 1.)
            proximo_z[0] += 1
        assert proximo_z[0] == i
        proximo_z[0] += 1
        return [rng.gauss(0., 1.) for _ in range(dim)]

    def linha(i, z):
        ze = sum(a * b for a, b in zip(z, e))
        zp = [a - ze * b for a, b in zip(z, e)]
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
        cb = min(max(ca - ga / s if s else ca + (1. if ga > 0 else -1.), FIM[0]), FIM[1])
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

    arq_l, arq_p = os.path.join(pasta, "linhas.csv"), os.path.join(pasta, "perdidas.csv")
    limpa_csv(arq_l, "linha,c,analises,p\n")
    limpa_csv(arq_p, "linha,motivo\n")
    feitas, perdidas = {}, {}
    for r in csv.DictReader(open(arq_l)):
        feitas[int(r["linha"])] = (float(r["c"]), int(r["analises"]))
    for r in csv.DictReader(open(arq_p)):
        perdidas[int(r["linha"])] = r["motivo"]
    fl, fp = open(arq_l, "a"), open(arq_p, "a")

    def estimativa():
        """P_f, desvio de uma linha e CoV(P_f) com as linhas feitas"""
        pi = [cauda(v[0]) for v in feitas.values()]
        if len(pi) < 2:
            return None, None, math.inf
        m = sum(pi) / len(pi)
        sd = math.sqrt(sum((x - m) ** 2 for x in pi) / (len(pi) - 1))
        # sd ~ 0 (todas as linhas com o mesmo c, ex.: no limite FIM; 1e-9 m cobre o arredondamento da media): CoV sem
        # sentido, nao para por ele
        return m, sd, (sd / (math.sqrt(len(pi)) * m) if m > 0 and sd > 1e-9 * m else math.inf)

    def tabela(final=False, motivo=""):
        c = [v[0] for v in feitas.values()]
        na = [v[1] for v in feitas.values()]
        pi = [cauda(x) for x in c]
        npe = len(perdidas)
        fr = npe / (npe + len(pi)) if npe + len(pi) else 0.
        comum = {"form": os.path.abspath(arq_form), "beta": beta, "s": s, "linhas": len(pi), "perdidas": npe,
                 "fracao_perdida": fr, "linhas_perdidas": sorted(perdidas), "analises_form": form.get("avaliacoes"),
                 "criterio": {"cov_alvo": op["cov"], "nmin": op["nmin"], "n_max": op["n"],
                              "regra": "para de submeter no 1o N >= nmin com CoV(P_f) <= cov_alvo (cov_alvo = 0: as "
                                       "n_max linhas); as linhas em voo terminam e entram"},
                 "motivo": motivo, "linhas_feitas": sorted(feitas), "semente": op["semente"], "K": K, "parametros": p}
        if final and npe:
            print(f"  linhas perdidas: {npe} ({100 * fr:.1f}%)" + ("  AVISO: mais de 1%, a P_f pode sair baixa"
                                                                 if fr > 0.01 else ""), flush=True)
        if len(pi) < 2:
            if final:
                print(f"  so {len(pi)} linha(s): sem P_f; parada: {motivo}", flush=True)
                json.dump(dict(comum, P_f=None, CoV=None, IC95=None, CoV_linha=None, analises=sum(na) + 2,
                               analises_linha=None, analises_cov5=None), open(os.path.join(pasta, "ls.json"), "w"),
                          indent=1)
            return
        m, sd, cov = estimativa()
        if not math.isfinite(cov):
            print(f"  linhas {len(pi):4d}  P_f = {m:.4e} (sem falhas ainda, ou todas iguais)", flush=True)
            if not final:
                return
            cov = math.inf
        nl5 = len(pi) * (cov / 0.05) ** 2
        npl = sum(na) / len(na)
        print(f"  linhas {len(pi):4d}  analises {sum(na) + 2:5d}  P_f = {m:.4e}  CoV = {100 * cov:5.2f}%  "
              f"CoV por linha = {100 * cov * math.sqrt(len(pi)):5.1f}%  analises/linha = {npl:.2f}  "
              f"linhas p/ CoV 5% = {nl5:5.0f}  analises p/ CoV 5% (sem o FORM) = {2 + npl * nl5:5.0f}", flush=True)
        if final:
            ic = 1.96 * sd / math.sqrt(len(pi))
            print(f"  P_f = {m:.4e}, IC 95% [{m - ic:.4e}, {m + ic:.4e}]; parada: {motivo}", flush=True)
            json.dump(dict(comum, P_f=m, CoV=cov if math.isfinite(cov) else None, IC95=[m - ic, m + ic],
                           CoV_linha=cov * math.sqrt(len(pi)) if math.isfinite(cov) else None, analises=sum(na) + 2,
                           analises_linha=npl, analises_cov5=2 + npl * nl5 if math.isfinite(cov) else None),
                      open(os.path.join(pasta, "ls.json"), "w"), indent=1)

    def parar():
        """o motivo para nao submeter mais linhas, ou "" """
        if os.path.exists(os.path.join(pasta, "PARAR")):
            return "PARAR"
        if op["cov"] > 0 and len(feitas) >= op["nmin"] and estimativa()[2] <= op["cov"]:
            return f"CoV <= {op['cov']} com {len(feitas)} linhas (>= nmin = {op['nmin']})"
        return ""

    falta = iter([i for i in range(op["n"]) if i not in feitas and i not in perdidas])
    motivo = parar()
    emvoo = set()
    with cf.ThreadPoolExecutor(op["paralelos"]) as ex:
        def submete():
            while not motivo and len(emvoo) < op["paralelos"]:
                i = next(falta, None)
                if i is None:
                    return
                fut = ex.submit(linha, i, z_da_linha(i))
                fut.linha = i
                emvoo.add(fut)

        submete()
        while emvoo:
            prontos, _ = cf.wait(emvoo, return_when=cf.FIRST_COMPLETED)
            for fut in prontos:
                emvoo.discard(fut)
                try:
                    i, c, n = fut.result()
                except Exception as err:  # uma analise que nao convergiu: a linha fica de fora (e avisa)
                    perdidas[fut.linha] = str(err).replace(",", ";").replace("\n", " ")
                    fp.write(f"{fut.linha},{perdidas[fut.linha]}\n")
                    fp.flush()
                    print(f"  linha {fut.linha} perdida: {err}", flush=True)
                    continue
                feitas[i] = (c, n)
                fl.write(f"{i},{c!r},{n},{cauda(c)!r}\n")
                fl.flush()
                if len(feitas) % 10 == 0:
                    tabela()
                if not motivo:
                    motivo = parar()
                    if motivo:
                        print(f"  parada: {motivo}; esperando as {len(emvoo)} linhas em voo", flush=True)
            submete()
    if not motivo:
        motivo = f"n maximo ({op['n']} linhas)" if op["cov"] > 0 else f"as n = {op['n']} linhas"
    tabela(final=True, motivo=motivo)


if __name__ == "__main__":
    main()
