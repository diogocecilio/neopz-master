#!/usr/bin/env python3
"""direcao_ls.py: a direcao do proximo estagio do line sampling em estagios (opcao C da nota
viabilidade_ls_adaptativo.md, secoes 5, 7 e 8), pela media condicional das linhas ja feitas.

    python3 direcao_ls.py <pasta>[:faixas] [<pasta>[:faixas] ...] saida=<form.json> [modos=5] [ref=<form.json>[,...]]
                          [form=<form.json>] [semente=<s>] [K=<K>]

<pasta>: uma pasta do line_sampling.py (linhas.csv; o form, a semente e K vem do ls.json, ou de form=, semente= e K=
quando a pasta ainda nao tem ls.json). faixas: as linhas usadas, ex. ls_p500_T18:50-199 ou ls_x:0-49,100-149 (padrao:
todas as do linhas.csv). As perdidas (perdidas.csv) ficam de fora (sem informacao).

Cada linha i (direcao alfa da pasta, semente s, K termos por campo) tem xi_perp_i = Z_i - (Z_i.alfa) alfa, com Z_i o
i-esimo vetor da semente regenerado como no line_sampling.py (random.Random(s), nc*K gauss por linha, em ordem; a mesma
aritmetica: a direcao do form.json com os M termos de cada campo nos K, e -u/|u| quando o beta do form.json e < 0), e a
raiz c_i (fator = 1 em xi_perp_i + c_i alfa, falha do lado de c grande). Ao longo da linha a falha e c > c_i e
    int_{c_i}^inf (xi_perp_i + t alfa) phi(t) dt = xi_perp_i Phi(-c_i) + alfa phi(c_i),
entao (nota, secao 5; Koutsourelakis et al. 2004)
    m = sum_i [xi_perp_i p_i + alfa phi(c_i)] / sum_i p_i,   p_i = Phi(-c_i),
estima E[xi | falha] (o ponto medio da falha) sem analises extras. Vale para qualquer alfa: linhas de pastas com
direcoes diferentes se somam, cada uma com a sua alfa (com a mesma semente as Z sao as mesmas: nao sao independentes;
o aviso sai). Mas somar um estagio de direcao ruim (p_i de cauda pesada: poucas linhas dominam) piora a direcao: no T0
(teste_ls_estagios.py, linear 100-D, partida a 60 graus, 150 + 100 linhas, 200 repeticoes) a direcao das linhas dos 2
estagios ficou a 8,5 graus (mediana; < 10 graus em 65%) da verdadeira, e a das linhas so do estagio 2 a 5,2 (< 10 graus
em 97%): para o 3o estagio, use so o ultimo estagio.
Supoe raiz unica em cada linha (a do line_sampling.py). Linha com raiz em FIM: c = 12 da peso ~0, c = -8
da xi_perp_i. Soma em ordem (pasta, linha): reprodutivel.
Truncamento: m e truncada nos modos 0..modos-1 de cada campo (k= e sinonimo) e normalizada, e = m_t/|m_t|. Padrao 5
(modos 0-4): a regra da nota (secoes 4 e 12) e manter poucos modos quando ha poucas falhas efetivas (ESS = (sum p)^2 /
sum p^2): sem truncar, o ruido dos ~nc*K modos domina a direcao (176 linhas do p5000, ESS 54: 38 graus com 100 modos,
8 com 0-4; com ESS 5, 0-9 ja piora). modos=0: sem truncar. Para escolher, a tabela dos modos (m_j, EP pelo metodo delta
do estimador de razao, z = m_j/EP) sai no terminal.
O form.json de saida: u = beta e, beta = -Phi^-1(P), P = media dos p_i das linhas usadas (a P_f delas). beta so fixa
onde o line_sampling.py comeca a busca e mede a inclinacao (c = beta e beta +- 0,1, 2 analises): o estimador do LS nao
tem vies para nenhuma direcao nem beta. Com P > 0,5 o beta sai < 0 (origem na falha) e o line_sampling.py usa
-u/|beta| = e (a mesma direcao). parametros: os do form.json da 1a pasta, com a lista "exe" que rodou (a do ls.json,
que inclui as chave=valor dadas ao line_sampling.py). avaliacoes: as analises das linhas usadas (o custo da direcao,
ja contado no estagio delas). origem: texto com as pastas, as faixas, as sementes e o truncamento. "direcao_ls": os
detalhes (fontes com as linhas usadas, ESS, m completa e EP, angulos), lidos pelo combina_ls.py para conferir que o
estagio seguinte nao reusa essas Z.

O proximo estagio TEM de usar uma semente fora de 'sementes_usadas' (as Z que fizeram a direcao nao podem voltar nas
linhas que a usam: secao 8 da nota; senao ha vies).
Sem numpy (roda no fermi).
"""
import csv
import json
import math
import os
import random
import sys
from statistics import NormalDist

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from form_kl import campos  # noqa: E402

MODOS_TABELA = 10  # modos de cada campo na tabela do terminal


def cauda(c):
    """Phi(-c)"""
    return 0.5 * math.erfc(c / math.sqrt(2.))


def dens(c):
    """phi(c)"""
    return math.exp(-0.5 * c * c) / math.sqrt(2. * math.pi)


def faixas(texto):
    """'50-199', '0-49,100-149' ou '7' -> o conjunto dos indices de linha"""
    s = set()
    for parte in texto.split(","):
        if "-" in parte:
            a, b = parte.split("-")
            s.update(range(int(a), int(b) + 1))
        elif parte.strip():
            s.add(int(parte))
    return s


def texto_faixas(idx):
    """[0, 1, 2, 5, 6] -> '0-2,5-6'"""
    idx, partes = sorted(idx), []
    for i in idx:
        if partes and i == partes[-1][1] + 1:
            partes[-1][1] = i
        else:
            partes.append([i, i])
    return ",".join(f"{a}-{b}" if b > a else f"{a}" for a, b in partes)


def direcao_form(arq_form, K):
    """a direcao e do line_sampling.py para o form.json (os M termos de cada campo nos K, normalizada; com beta < 0 no
    form.json, -u/|u|), o beta com sinal, nc, o json e os parametros (a mesma aritmetica do line_sampling.py)"""
    form = json.load(open(arq_form))
    p = dict(form["parametros"])
    p.setdefault("prog", "cho")
    nc = len(campos(p))
    u = form["u"]
    M = len(u) // nc
    e = [0.] * (nc * K)
    for f in range(nc):
        e[f * K:f * K + M] = u[f * M:(f + 1) * M]
    beta = math.sqrt(sum(x * x for x in e))
    e = [x / beta for x in e]
    if form.get("beta", 1.) < 0.:
        e = [-x for x in e]
        beta = -beta
    return e, beta, nc, form, p


def le_linhas(arq, sel):
    """{linha: (c, analises)} do linhas.csv (so as da faixa; uma linha incompleta no fim, de uma corrida em andamento,
    fica de fora)"""
    linhas = {}
    with open(arq, newline="") as f:
        for r in csv.DictReader(f):
            try:
                i, c, n = int(r["linha"]), float(r["c"]), int(r["analises"])
            except (TypeError, ValueError):
                continue
            if sel is None or i in sel:
                linhas[i] = (c, n)
    return linhas


def le_pasta(arg, op):
    """uma pasta do LS, com as faixas opcionais: a direcao, a semente, K, as linhas feitas e perdidas, os parametros"""
    pasta, sel = arg, None
    if ":" in arg:
        pasta, txt = arg.rsplit(":", 1)
        sel = faixas(txt)
    pasta = os.path.abspath(pasta)
    arq_ls = os.path.join(pasta, "ls.json")
    lsj = json.load(open(arq_ls)) if os.path.exists(arq_ls) else {}
    try:  # o ls.json manda; form=, semente= e K= so para a pasta sem ls.json
        arq_form = lsj["form"] if "form" in lsj else op["form"]
        semente = int(lsj["semente"] if "semente" in lsj else op["semente"])
        K = int(lsj["K"] if "K" in lsj else op["K"])
    except KeyError:
        sys.exit(f"{pasta}: sem ls.json (corrida em andamento?): de form=, semente= e K=")
    e, beta, nc, form, p = direcao_form(arq_form, K)
    arq_s = os.path.join(pasta, "inclinacao.json")
    if os.path.exists(arq_s):
        b2 = json.load(open(arq_s))["beta"]
        if abs(b2 - beta) > 1e-9 * max(1., abs(beta)):
            sys.exit(f"{pasta}: o beta do inclinacao.json ({b2}) nao e o do {arq_form} ({beta}): o form.json mudou?")
    linhas = le_linhas(os.path.join(pasta, "linhas.csv"), sel)
    perdidas = set()
    arq_p = os.path.join(pasta, "perdidas.csv")
    if os.path.exists(arq_p):
        with open(arq_p, newline="") as f:
            for r in csv.DictReader(f):
                try:
                    i = int(r["linha"])
                except (TypeError, ValueError):
                    continue
                if sel is None or i in sel:
                    perdidas.add(i)
    faltam = sorted(sel - set(linhas) - perdidas) if sel is not None else []
    exe = lsj.get("parametros", {}).get("exe", p.get("exe"))
    return {"pasta": pasta, "sel": sel, "semente": semente, "K": K, "nc": nc, "e": e, "beta": beta, "form": arq_form,
            "json": form, "p": p, "exe": exe, "linhas": linhas, "perdidas": perdidas, "faltam": faltam}


def angulo(a, b):
    """o angulo entre a e b (graus)"""
    na, nb = math.sqrt(sum(x * x for x in a)), math.sqrt(sum(x * x for x in b))
    return math.degrees(math.acos(max(-1., min(1., sum(x * y for x, y in zip(a, b)) / (na * nb)))))


def media_condicional(fontes, dim):
    """os vetores a_i = xi_perp_i p_i + alfa phi(c_i) e os p_i das linhas usadas, em ordem (pasta, linha)"""
    A, B = [], []
    for f in fontes:
        if not f["linhas"]:
            continue
        rng, e = random.Random(f["semente"]), f["e"]
        for i in range(max(f["linhas"]) + 1):
            z = [rng.gauss(0., 1.) for _ in range(dim)]  # Z[i], como o z_da_linha do line_sampling.py
            if i not in f["linhas"]:
                continue
            ze = sum(a * b for a, b in zip(z, e))
            zp = [a - ze * b for a, b in zip(z, e)]
            c = f["linhas"][i][0]
            pi, fi = cauda(c), dens(c)
            A.append([a * pi + b * fi for a, b in zip(zp, e)])
            B.append(pi)
    return A, B


def main():
    op, args = {"modos": "5", "saida": "", "ref": ""}, []
    for a in sys.argv[1:]:
        if "=" in a:
            k, v = a.split("=", 1)
            op["modos" if k == "k" else k] = v
        else:
            args.append(a)
    if not args or not op["saida"]:
        sys.exit(__doc__)
    fontes = [le_pasta(a, op) for a in args]
    nc, K = fontes[0]["nc"], fontes[0]["K"]
    if any((f["nc"], f["K"]) != (nc, K) for f in fontes):
        sys.exit("as pastas tem campos ou K diferentes")
    for f in fontes[1:]:
        if f["exe"] != fontes[0]["exe"] or f["p"].get("T") != fontes[0]["p"].get("T"):
            print(f"AVISO: {f['pasta']} rodou com outras opcoes do executavel (ou outro T) que {fontes[0]['pasta']}",
                  flush=True)
    sementes = [f["semente"] for f in fontes]
    if len(set(sementes)) < len(sementes):
        print("AVISO: pastas com a mesma semente: as mesmas Z (linhas de mesmo indice nao sao independentes)")
    dim = nc * K
    A, B = media_condicional(fontes, dim)
    n, sB = len(B), sum(B)
    if n < 2 or sB <= 0.:
        sys.exit(f"so {n} linha(s), soma dos p_i {sB}: sem direcao")
    m = [sum(col) / sB for col in zip(*A)]
    P = sB / n
    sd = math.sqrt(sum((x - P) ** 2 for x in B) / (n - 1))
    ess = sB * sB / sum(x * x for x in B)
    ep = [math.sqrt(sum((a[j] - m[j] * b) ** 2 for a, b in zip(A, B)) / (n * (n - 1))) / P for j in range(dim)]
    k = int(op["modos"])
    mt = [x if k <= 0 or j % K < k else 0. for j, x in enumerate(m)]
    nm, nmt = math.sqrt(sum(x * x for x in m)), math.sqrt(sum(x * x for x in mt))
    et = [x / nmt for x in mt]
    beta = -NormalDist().inv_cdf(min(max(P, 1e-300), 1. - 1e-15))
    if abs(beta) < 1e-3:  # P ~ 0,5: u = 0 nao serve ao line_sampling.py
        beta = math.copysign(1e-3, beta)
    nan = sum(f["linhas"][i][1] for f in fontes for i in f["linhas"])

    print(f"direcao_ls.py: {n} linhas, P = {P:.4e} (EP {sd / math.sqrt(n):.2e}, CoV por linha {sd / P:.2f}), "
          f"ESS = {ess:.1f}, |m| = {nm:.3f}, |m truncada| = {nmt:.3f}, beta = -Phi^-1(P) = {beta:.4f}")
    for f in fontes:
        print(f"  {f['pasta']}" + (f":{texto_faixas(f['sel'])}" if f["sel"] is not None else "") +
              f"  semente {f['semente']}  linhas {len(f['linhas'])}  perdidas {len(f['perdidas'])}"
              + (f"  AINDA NAO FEITAS {len(f['faltam'])}" if f["faltam"] else "") + f"  form {f['form']}")
    print(f"  modos (m_j, EP, z = m_j/EP; os {MODOS_TABELA} primeiros de cada campo; truncamento: " +
          (f"modos 0-{k - 1} de cada campo)" if k > 0 else "nenhum)"))
    for fc, nome in enumerate(campos(fontes[0]["p"])):
        print(f"    campo {nome}: " + "  ".join(f"{j}:{m[fc * K + j]:+.3f}({ep[fc * K + j]:.3f},{m[fc * K + j] / ep[fc * K + j] if ep[fc * K + j] > 0 else 0.:+.1f})" for j in range(min(MODOS_TABELA, K))))
    angs = {}
    for f in fontes:
        angs[f["pasta"]] = angulo(et, f["e"])
        print(f"  angulo com a direcao de {os.path.basename(f['pasta'])}: {angs[f['pasta']]:.1f} graus "
              f"(m sem truncar: {angulo(m, f['e']):.1f})")
    for r in filter(None, op["ref"].split(",")):
        er = direcao_form(r, K)[0]
        angs[os.path.abspath(r)] = angulo(et, er)
        print(f"  angulo com {r}: {angs[os.path.abspath(r)]:.1f} graus (m sem truncar: {angulo(m, er):.1f})")

    params = dict(fontes[0]["json"]["parametros"])
    if fontes[0]["exe"] is not None:
        params["exe"] = list(fontes[0]["exe"])
    usadas = sorted(set(sementes))
    origem = ("direcao_ls.py: media condicional E[xi | falha] das linhas de " +
              "; ".join(f"{os.path.basename(f['pasta'])}:{texto_faixas(f['linhas'])} (semente {f['semente']})"
                        for f in fontes) +
              f" ({n} linhas, ESS {ess:.1f}), " + (f"truncada nos modos 0-{k - 1} de cada campo" if k > 0 else
                                                  "sem truncar") +
              f"; beta = -Phi^-1(P_f das linhas) = {beta:.4f}; o proximo estagio usa semente fora de {usadas}")
    det = {"fontes": [{"pasta": f["pasta"], "form": os.path.abspath(f["form"]), "semente": f["semente"],
                       "linhas": sorted(f["linhas"]), "perdidas": sorted(f["perdidas"]), "faltam": f["faltam"]}
                      for f in fontes],
           "sementes_usadas": usadas, "modos": k, "linhas": n, "ESS": ess, "P_f": P, "EP_P_f": sd / math.sqrt(n),
           "CoV_linha": sd / P, "norma_m": nm, "norma_m_truncada": nmt, "m": m, "ep_m": ep, "angulos": angs}
    out = {"parametros": params, "beta": beta, "P_f": P, "u": [beta * x for x in et], "avaliacoes": nan,
           "origem": origem, "direcao_ls": det}
    json.dump(out, open(op["saida"], "w"), indent=1)
    print(f"  gravado {op['saida']}; o proximo estagio: line_sampling.py <exe> {op['saida']} <pasta nova> "
          f"semente=<fora de {usadas}> ...")


if __name__ == "__main__":
    main()
