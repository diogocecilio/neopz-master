#!/usr/bin/env python3
"""line_sampling.py: line sampling (amostragem por linhas) na direcao do FORM, com o campo completo.

    python3 line_sampling.py <exe> <form.json> <pasta> [n=100] [semente=7] [paralelos=14] [tol=0.01] [K=2000]
                             [cov=0] [nmin=30] [busca=secante|cerca] [passo_cerca=0.25] [folga_cerca=1]
                             [reuso=0|1] [custos=0|1]
                             [chave=valor das opcoes do programa: nystrom=, metodo=, srm_f0=, contorno=, tfim=, ...]

A direcao e = u*/beta vem do FORM (form_kl.py, M termos por campo: c', phi' no Cho2010; c', phi', k no rebaixamento e na
chuva),
completada com zeros ate os K termos de cada campo. As analises usam o campo com K termos (M=K) e as opcoes do FORM,
trocadas pelas chave=valor dadas (ex.: nystrom=<a base completa>, metodo=srm). No rebaixamento, as opcoes contorno=,
agua=, trebaix=, tfim=, dt=, hw=, covc=, covphi=, covk=, s= e alfa= vem do form.json (form.json antigos: os padroes do
SlopeDrawdownMonteCarlo). Na chuva (prog=chuva), as opcoes do executavel vem do form.json (a lista "exe") e as
chave=valor dadas aqui trocam ou acrescentam opcoes nela (ex.: nystrom=<outra base>, pressoes=<outro arquivo>); o fator
de cada analise e o menor FS dos instantes (form_kl.le_fator).
A falha na chuva (prog=chuva): g = ln(FS_min), com FS_min o menor FS dos instantes de tempos= (as opcoes do
executavel no form.json, ou trocadas aqui): a primeira passagem nesses instantes. Um instante sem equilibrio no SRM
(situacao fator com FS_min = 0) e falha, FS = 0,005 (ln finito); um transiente que parou antes de um instante (situacao
transiente) tira a linha da conta (perdidas.csv). Com busca=cerca, a raiz e a troca de sinal de FS_min - 1 ao longo da
linha. A P_f e a do mesmo evento que o Monte Carlo direto com os mesmos tempos= (o menor FS dos instantes abaixo de 1);
com um instante so (ex.: tempos=24, o fim da chuva no f50), e P(FS(24 h) < 1).
Cada linha: sorteia z ~ N(0, I) com nc K termos, tira a componente em e (zp = z - (z.e) e) e acha, na reta zp + c e, o
c_i em que o fator (FS ou Gamma) vale 1. Como zp e c sao independentes, P_f = E[Phi(-c_i)], e cada linha da
p_i = Phi(-c_i): so a variacao de c_i entre as linhas (a nao linearidade da superficie de falha) entra na variancia. O
estimador nao tem vies para nenhuma direcao; uma direcao ruim so aumenta a dispersao.

A raiz e pela secante em g(c) = ln(fator), quase linear em c: a 1a analise em c = beta, a 2a no passo de Newton com a
inclinacao s medida no ponto de projeto (2 analises por direcao), e a raiz da secante e aceita sem nova analise quando o
passo for menor que tol (em c). Os mesmos z para qualquer direcao com a mesma semente (numeros aleatorios comuns).
Sem numpy (roda no fermi).
Cuidado: esse teste mede so o passo. Com um penhasco no fator (a coluna em 24 h cai do FS ~1,03 ao patamar 0,366875
da coluna saturada), a secante em ln(fator) fica presa numa ponta e aceita raizes que nunca foram cercadas.

busca=cerca (o padrao continua busca=secante, igual bit a bit): a raiz pelo sinal de fator - 1, com a cerca garantida,
sem as 2 analises da inclinacao (nao servem a esta busca, entao nao sao feitas).
  1. Cerca: a 1a analise em c0 = beta; dali para fora (para c maior se c0 e seguro, para c menor se ja falha), com
     passo crescente (passo_cerca, que dobra a cada analise), ate o sinal de fator - 1 trocar. So o sinal conta, entao
     um patamar nao atrapalha. Da 3a analise em diante, a secante dos dois ultimos pontos mira um pouco alem da raiz
     prevista (tol/4), entre 1/4 e 1 vez o passo da vez. Os limites FIM e a origem na falha seguem as regras de sempre:
     segura ate FIM[1] da c = FIM[1]; falha desde FIM[0] da c = FIM[0].
  2. Refino: ITP (Oliveira e Takahashi 2020, ACM TOMS 47:5), com a regula falsi em ln(fator) e os valores do Illinois
     (o valor da ponta mantida duas vezes seguidas cai a metade).
     - Se o ponto cai a menos de tol de uma ponta, vai para 0,98 tol dela, para fechar a cerca numa analise.
     - Depois, e projetado no intervalo em torno do meio que garante o pior caso: n_1/2 + folga_cerca analises no
       refino, com n_1/2 = ceil(log2(largura inicial / tol)), as da bissecao. Num degrau, e a bissecao + folga_cerca.
     - Para quando a CERCA (o ultimo ponto seguro e o ultimo em falha, ja analisados) tem largura <= tol.
  A raiz devolvida e o meio da cerca (erro <= tol/2). Sem teto de analises: o ITP ja limita o refino.
  O linhas.csv ganha as colunas seg e fal (a cerca final; fal fica vazia na linha segura ate FIM[1], e seg na linha
  em falha desde FIM[0]). O ls.json ganha busca, tol, passo_cerca, folga_cerca, largura_max (a maior cerca final) e
  linhas_sem_cerca, e as analises nao somam as 2 da inclinacao (o combina_ls.py ainda soma 2 por estagio).
  Retomada: a busca da pasta fica no busca.json (e no ls.json). Uma pasta da outra busca e recusada (a da secante tem o
  inclinacao.json).

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

reuso=1 (so na chuva; padrao 0, o de antes bit a bit): o transiente calculado uma vez por linha. Com a direcao sem
componente em k (o FORM com gradiente=modal: o fator congela as pressoes e da 0 nos termos de k), o ponto zp + c e tem
os xi de k iguais aos de z (bit a bit) em toda a linha: so c' e phi' mudam, entao o transiente (as pressoes) e o
mesmo em todos os pontos dela. A 1a analise da linha i calcula o transiente e o grava (grava_pressoes=<pasta>/
l<i>_pw.txt); as outras o leem (pressoes=) e so calculam o fator. O arquivo e apagado no fim da linha (e no comeco,
se sobrou de uma rodada interrompida). O k continua aleatorio ENTRE as linhas (na parte ortogonal zp), e o estimador
nao muda (o mesmo FS bit a bit: as pressoes vao no arquivo com 17 algarismos). Garantias de que um arquivo so serve ao
seu campo de k e as suas opcoes do transiente:
  - a direcao tem de ter os K termos de k exatamente 0 (senao recusa); cada analise confere que os xi de k sao os da
    1a analise da linha;
  - o arquivo e da linha (um por linha, numa pasta so com as opcoes do form.json);
  - o executavel confere a chave do transiente na leitura (ChaveTransiente: geometria, nivel, contornos, solo, chuva,
    dt, malha, Newton...), que tem de ter tambem o campo de k (campo_k: CoV, lh, lv, M, base e os bits dos xi de k).
    Isso precisa do patch do campo de k no montecarlo.cpp; sem ele o executavel recusa grava_pressoes= com k
    aleatorio. Uma analise barata confere antes de tudo (form_kl.sonda_reuso, com tempos=0 fator=tri0 nref=0) e sai
    com a explicacao.
  Sem pressoes=, grava_pressoes= e so_transiente= nas opcoes do executavel (com k fixo, o transiente ja vem de um
  arquivo so). O ls.json ganha "reuso": 1.
custos=1 (padrao: 1 com reuso=1, 0 sem): custos.csv com uma linha por analise (linha, analise, c, modo, s, s_transiente,
s_fator): modo calcula (sem reuso), grava (a 1a da linha, com reuso) ou le; s o tempo da analise (s), s_transiente e
s_fator os do executavel (o transiente ate o ultimo instante e a soma dos fatores dos instantes; s_transiente = 0 na
leitura).

Na coluna (prog=coluna, form_kl.py): como na chuva (as opcoes do executavel na lista "exe" do form.json, trocadas
pelas chave=valor dadas aqui), com o fator de form_kl.le_fator (FS_passagem em T, ou t_passagem / T censurado em 2), o
fator limitado por baixo a GAMMA_ZERO (ln finito) e os arquivos _perfil_, _xi_ e _campo_ tambem apagados.

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
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from form_kl import REBAIX_PADRAO, TRANSIENTE_ARQ, ambiente, campos, comando, le_fator, sonda_reuso  # noqa: E402

FIM = (-8., 12.)
GAMMA_ZERO = 0.005  # Gamma = 0 (sem equilibrio nem com 1% da gravidade): falha, com ln finito para a secante
SUFIXOS_COLUNA = ("_perfil_0.csv", "_xi_0.txt", "_campo_0_0.csv")  # prog=coluna: os outros arquivos da analise


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
    op = {"n": 100, "semente": 7, "paralelos": 14, "tol": 0.01, "K": 2000, "cov": 0., "nmin": 30, "busca": "secante",
          "passo_cerca": 0.25, "folga_cerca": 1, "reuso": 0, "custos": -1}
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
        elif p["prog"] in ("chuva", "coluna") and k not in ("threads",):  # uma opcao do executavel: troca ou acrescenta
            p["exe"] = [x for x in p["exe"] if not x.startswith(k + "=")] + [a]
        else:
            p[k] = type(p[k])(v) if k in p and not isinstance(p[k], str) else v
    if op["busca"] not in ("secante", "cerca"):
        sys.exit("busca: secante (padrao) ou cerca")
    reuso = op["reuso"] != 0
    custos = (op["custos"] != 0 if op["custos"] >= 0 else reuso) and p["prog"] == "chuva"
    if reuso and p["prog"] != "chuva":
        sys.exit("reuso=1: so na chuva (prog=chuva)")
    if reuso and any(x.startswith(TRANSIENTE_ARQ) for x in p["exe"]):
        sys.exit("reuso=1: sem pressoes=, grava_pressoes= e so_transiente= nas opcoes do executavel (com k fixo o "
                 "transiente ja vem de um arquivo so)")
    cerca = op["busca"] == "cerca"
    if cerca and (op["tol"] <= 0. or op["passo_cerca"] <= 0. or op["folga_cerca"] < 0):
        sys.exit("busca=cerca: tol > 0, passo_cerca > 0 e folga_cerca >= 0")
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
    if reuso:  # o transiente so e o mesmo ao longo da linha se a direcao nao tem componente em k (o 3o campo)
        nk = sum(1 for x in e[2 * K:] if x != 0.)
        if nk:
            sys.exit(f"reuso=1: a direcao tem {nk} termos de k nao nulos (o FORM nao e o modal?): o transiente mudaria "
                     "ao longo da linha; rode sem reuso")
    # a busca da pasta (retomada): o busca.json (so a da cerca o grava), o ls.json (sem a chave: secante) ou, sem os
    # dois, uma pasta ja usada (inclinacao.json, linhas.csv ou analises.csv) e da secante. Outra busca: recusa
    arq_b, arq_s = os.path.join(pasta, "busca.json"), os.path.join(pasta, "inclinacao.json")

    def busca_da_pasta():
        for arq in (arq_b, os.path.join(pasta, "ls.json")):
            try:
                return json.load(open(arq)).get("busca", "secante")
            except (OSError, ValueError):
                pass
        usada = [x for x in ("linhas.csv", "analises.csv") if os.path.exists(os.path.join(pasta, x))
                 and os.path.getsize(os.path.join(pasta, x)) > 0]
        return "secante" if os.path.exists(arq_s) or usada else None

    anterior = busca_da_pasta()
    if anterior is not None and anterior != op["busca"]:
        sys.exit(f"{pasta}: pasta da busca={anterior}, esta rodada pede busca={op['busca']}: retomada recusada (use "
                 f"outra pasta)")
    if cerca:
        cfg = {"busca": "cerca", "beta": beta, "tol": op["tol"], "passo_cerca": op["passo_cerca"],
               "folga_cerca": op["folga_cerca"]}
        if os.path.exists(arq_b):
            antes = json.load(open(arq_b))
            dif = [k for k in cfg if antes.get(k) != cfg[k]]
            if dif:
                print(f"  AVISO: retomada com {', '.join(f'{k} {antes.get(k)} -> {cfg[k]}' for k in dif)}", flush=True)
        else:
            json.dump(cfg, open(arq_b, "w"), indent=1)
    if reuso:  # o executavel grava o transiente com k aleatorio e o campo de k na chave? (sai daqui se nao)
        sonda_reuso(exe, p, pasta)
    trava = threading.Lock()
    limpa_csv(os.path.join(pasta, "analises.csv"), "")
    reg = open(os.path.join(pasta, "analises.csv"), "a")
    if custos:
        limpa_csv(os.path.join(pasta, "custos.csv"), "linha,analise,c,modo,s,s_transiente,s_fator\n")
        reg_c = open(os.path.join(pasta, "custos.csv"), "a")

    def custo(arq_csv):
        """s_transiente (o do ultimo instante: o transiente ate ali) e s_fator (a soma dos instantes) do executavel"""
        with open(arq_csv, newline="") as f:
            linhas = list(csv.DictReader(f))
        st = max((float(r.get("s_transiente") or 0.) for r in linhas), default=0.)
        return st, sum(float(r.get("s_fator") or 0.) for r in linhas)

    def fator(x, rotulo, linha, c, extra=()):
        arq = os.path.join(pasta, rotulo + ".xi")
        with open(arq, "w") as f:
            f.write("\n".join(repr(t) for t in x) + "\n")
        t0 = time.monotonic()
        subprocess.run(comando(exe, p, arq, rotulo) + list(extra), check=True, stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL, env=ambiente(p), cwd=pasta)
        s = time.monotonic() - t0
        try:
            if custos:
                modo = "calcula" if not extra else "grava" if extra[0].startswith("grava_pressoes=") else "le"
                st, sf = custo(os.path.join(pasta, rotulo + "_0_0.csv"))
                with trava:
                    reg_c.write(f"{linha},{rotulo.rsplit('_', 1)[-1]},{c!r},{modo},{s:.2f},{st:.2f},{sf:.2f}\n")
                    reg_c.flush()
            if p["prog"] == "chuva":  # o menor FS dos instantes (um instante sem equilibrio e falha)
                G = le_fator(os.path.join(pasta, rotulo + "_0_0.csv"), p)[0]
            elif p["prog"] == "coluna":  # FS_passagem em T ou t_passagem / T (ln finito)
                G = max(le_fator(os.path.join(pasta, rotulo + "_0_0.csv"), p)[0], GAMMA_ZERO)
            else:
                with open(os.path.join(pasta, rotulo + "_0_0.csv"), newline="") as f:
                    r = next(csv.DictReader(f))
                G = float(r["FS_min"])
                if r["situacao"] == "fator" and G == 0.:  # GIM sem equilibrio nem com lambda < 0,01: falha
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
            reg.write(f"{linha},{c!r},{G!r}\n")
            reg.flush()
        return G

    def ponto(zp, c):
        return [a + c * b for a, b in zip(zp, e)]

    # reuso=1: o transiente da linha i no l<i>_pw.txt (gravado pela 1a analise, lido pelas outras) e os xi de k da 1a
    # analise de cada linha em voo (as outras conferem)
    kxi = {}

    def arq_pw(i):
        return os.path.join(pasta, f"l{i}_pw.txt")

    def apaga_pw(i):
        for suf in ("", ".tmp"):
            try:
                os.remove(arq_pw(i) + suf)
            except FileNotFoundError:
                pass

    def transiente(i, j, x):
        """as opcoes do transiente da analise j da linha i no ponto x (vazias sem reuso)"""
        if not reuso:
            return ()
        if j == 0:
            kxi[i] = x[2 * K:]
            return (f"grava_pressoes={arq_pw(i)}",)
        if x[2 * K:] != kxi[i]:
            raise RuntimeError(f"linha {i}: os xi de k mudaram ao longo da linha (reuso=1)")
        return (f"pressoes={arq_pw(i)}",)

    # a inclinacao de g = ln(fator) ao longo de e, no ponto de projeto (2 analises; so na busca=secante)
    if cerca:
        s = None
    elif os.path.exists(arq_s):
        s = json.load(open(arq_s))["s"]
    else:
        zero = [0.] * dim
        with cf.ThreadPoolExecutor(2) as ex:
            gp, gm = ex.map(lambda d: math.log(fator(ponto(zero, beta + d), f"s{'p' if d > 0 else 'm'}", -1, beta + d)),
                            (0.1, -0.1))
        s = (gp - gm) / 0.2
        json.dump({"s": s, "beta": beta}, open(arq_s, "w"))
    extra = 0 if cerca else 2  # as analises fora das linhas (as 2 da inclinacao)
    if cerca:
        print(f"direcao: {arq_form}  beta = {beta:.4f}  P_f FORM = {cauda(beta):.4e}  busca = cerca (tol {op['tol']:g},"
              f" passo_cerca {op['passo_cerca']:g}, folga_cerca {op['folga_cerca']}; sem a inclinacao)  "
              f"(campo com {K} termos, metodo {p.get('metodo', '-')})", flush=True)
    else:
        print(f"direcao: {arq_form}  beta = {beta:.4f}  P_f FORM = {cauda(beta):.4e}  d ln(fator)/dc = {s:.4f}  "
              f"(campo com {K} termos, metodo {p.get('metodo', '-')})", flush=True)

    if p["prog"] in ("rebaixamento", "chuva", "coluna"):
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
        def g(c, j):
            x = ponto(zp, c)
            return math.log(fator(x, f"l{i}_{j}", i, c, transiente(i, j, x)))
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

    def linha_cerca(i, z):
        """busca=cerca: (i, c, analises, seg, fal), com seg < fal a cerca final (None na linha sem troca de sinal ate
        FIM). g = ln(fator): g > 0 seguro, g <= 0 falha (a falha do lado de c grande); a cerca usa so o sinal"""
        ze = sum(a * b for a, b in zip(z, e))
        zp = [a - ze * b for a, b in zip(z, e)]
        tol, conta = op["tol"], [0]

        def g(c):
            x = ponto(zp, c)
            v = math.log(fator(x, f"l{i}_{conta[0]}", i, c, transiente(i, conta[0], x)))
            conta[0] += 1
            return v

        # 1. a cerca: de c0 = beta para fora, com passo crescente, ate o sinal trocar
        cu, gu = beta, g(beta)
        if gu == 0.:
            return i, cu, conta[0], cu, cu
        d = 1. if gu > 0. else -1.  # c0 seguro: a falha a frente (c maior); c0 em falha: o seguro atras
        ant, h = None, op["passo_cerca"]
        while True:
            if gu > 0. and cu >= FIM[1]:
                return i, FIM[1], conta[0], cu, None   # segura ate FIM[1]
            if gu <= 0. and cu <= FIM[0]:
                return i, FIM[0], conta[0], None, cu   # falha desde FIM[0]
            dc = h
            if ant is not None and gu != ant[1]:  # a secante dos dois ultimos (do mesmo lado) mira alem da raiz
                dist = d * (cu - gu * (cu - ant[0]) / (gu - ant[1]) - cu)
                if dist > 0.:
                    dc = min(max(dist + 0.25 * tol, 0.25 * h), h)
            cn = min(max(cu + d * dc, FIM[0]), FIM[1])
            gn = g(cn)
            if gn == 0.:
                return i, cn, conta[0], cn, cn
            if (gn > 0.) != (gu > 0.):
                break
            ant, cu, gu, h = (cu, gu), cn, gn, 2. * h
        a, ga, b, gb = (cu, gu, cn, gn) if d > 0. else (cn, gn, cu, gu)  # a seguro < b em falha
        # 2. o refino (ITP): regula falsi com os valores do Illinois (fa, fb), puxada para fechar a cerca e projetada
        # no intervalo [meio - r, meio + r] que garante largura <= tol em nmax analises
        nmax = max(0, math.ceil(math.log2((b - a) / tol))) + op["folga_cerca"]
        fa, fb, ult, j = ga, gb, 0, 0
        while b - a > tol:
            meio = 0.5 * (a + b)
            r = max(0.5 * tol * 2. ** (nmax - j) - 0.5 * (b - a), 0.)
            x = (a * fb - b * fa) / (fb - fa)
            if b - a <= 1.96 * tol:  # qualquer x em [b - 0,98 tol, a + 0,98 tol] fecha a cerca
                x = min(max(x, b - 0.98 * tol), a + 0.98 * tol)
            elif x - a < tol:
                x = a + 0.98 * tol
            elif b - x < tol:
                x = b - 0.98 * tol
            if abs(x - meio) > r:
                x = meio + math.copysign(r, x - meio)
            gx = g(x)
            if gx == 0.:
                return i, x, conta[0], x, x
            if gx > 0.:
                a, fa = x, gx
                if ult < 0:
                    fb *= 0.5
                ult = -1
            else:
                b, fb = x, gx
                if ult > 0:
                    fa *= 0.5
                ult = 1
            j += 1
        return i, 0.5 * (a + b), conta[0], a, b

    def roda_linha(i, z):
        """a busca da linha i; com reuso, o arquivo do transiente dela apagado antes (sobra de uma rodada interrompida) e
        depois"""
        busca = linha_cerca if cerca else linha
        if not reuso:
            return busca(i, z)
        apaga_pw(i)
        try:
            return busca(i, z)
        finally:
            apaga_pw(i)
            kxi.pop(i, None)

    def num(x):
        return "" if x is None else repr(x)

    arq_l, arq_p = os.path.join(pasta, "linhas.csv"), os.path.join(pasta, "perdidas.csv")
    limpa_csv(arq_l, "linha,c,analises,p,seg,fal\n" if cerca else "linha,c,analises,p\n")
    limpa_csv(arq_p, "linha,motivo\n")
    feitas, perdidas, cercas = {}, {}, {}
    for r in csv.DictReader(open(arq_l)):
        feitas[int(r["linha"])] = (float(r["c"]), int(r["analises"]))
        if cerca:
            cercas[int(r["linha"])] = tuple(float(r[k]) if r[k] else None for k in ("seg", "fal"))
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
        if reuso:
            comum["reuso"] = 1
        if cerca:
            larg = [b - a for a, b in cercas.values() if a is not None and b is not None]
            comum.update({"busca": "cerca", "tol": op["tol"], "passo_cerca": op["passo_cerca"],
                          "folga_cerca": op["folga_cerca"], "largura_max": max(larg) if larg else None,
                          "linhas_sem_cerca": sum(a is None or b is None for a, b in cercas.values())})
        if final and npe:
            print(f"  linhas perdidas: {npe} ({100 * fr:.1f}%)" + ("  AVISO: mais de 1%, a P_f pode sair baixa"
                                                                 if fr > 0.01 else ""), flush=True)
        if len(pi) < 2:
            if final:
                print(f"  so {len(pi)} linha(s): sem P_f; parada: {motivo}", flush=True)
                json.dump(dict(comum, P_f=None, CoV=None, IC95=None, CoV_linha=None, analises=sum(na) + extra,
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
        print(f"  linhas {len(pi):4d}  analises {sum(na) + extra:5d}  P_f = {m:.4e}  CoV = {100 * cov:5.2f}%  "
              f"CoV por linha = {100 * cov * math.sqrt(len(pi)):5.1f}%  analises/linha = {npl:.2f}  "
              f"linhas p/ CoV 5% = {nl5:5.0f}  analises p/ CoV 5% (sem o FORM) = {extra + npl * nl5:5.0f}", flush=True)
        if final:
            ic = 1.96 * sd / math.sqrt(len(pi))
            print(f"  P_f = {m:.4e}, IC 95% [{m - ic:.4e}, {m + ic:.4e}]; parada: {motivo}", flush=True)
            json.dump(dict(comum, P_f=m, CoV=cov if math.isfinite(cov) else None, IC95=[m - ic, m + ic],
                           CoV_linha=cov * math.sqrt(len(pi)) if math.isfinite(cov) else None,
                           analises=sum(na) + extra, analises_linha=npl,
                           analises_cov5=extra + npl * nl5 if math.isfinite(cov) else None),
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
                fut = ex.submit(roda_linha, i, z_da_linha(i))
                fut.linha = i
                emvoo.add(fut)

        submete()
        while emvoo:
            prontos, _ = cf.wait(emvoo, return_when=cf.FIRST_COMPLETED)
            for fut in prontos:
                emvoo.discard(fut)
                try:
                    res = fut.result()
                    i, c, n = res[:3]
                except Exception as err:  # uma analise que nao convergiu: a linha fica de fora (e avisa)
                    perdidas[fut.linha] = str(err).replace(",", ";").replace("\n", " ")
                    fp.write(f"{fut.linha},{perdidas[fut.linha]}\n")
                    fp.flush()
                    print(f"  linha {fut.linha} perdida: {err}", flush=True)
                    continue
                feitas[i] = (c, n)
                if cerca:
                    cercas[i] = res[3:]
                    fl.write(f"{i},{c!r},{n},{cauda(c)!r},{num(res[3])},{num(res[4])}\n")
                else:
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
