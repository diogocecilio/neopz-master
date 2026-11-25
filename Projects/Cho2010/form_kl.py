#!/usr/bin/env python3
"""form_kl.py: FORM (HL-RF) do Exemplo 2 de Cho (2010) com os campos de c' e phi' truncados em M termos da KL cada.

As variaveis sao os 2M coeficientes xi da KL: normais padrao independentes (a correlacao cruzada r(c, phi) e a
transformacao lognormal ficam dentro do Cho2010), entao nao ha transformacao de Nataf. A funcao de estado limite e
    g(u) = FS(u) - 1   (metodo=srm, o padrao)   ou   g(u) = Gamma(u) - 1   (metodo=gim),
com o fator do Cho2010 (SRM ou GIM, com refinamento adaptativo) nos xi dados (opcao xi=). Os dois tem a mesma regiao de
falha (FS < 1 se e so se Gamma < 1), entao o mesmo ponto de projeto; o g e o gradiente e que mudam.
  - gradiente por diferencas centrais com passo h em u: 2n + 1 avaliacoes por iteracao (n = 2M), em paralelo. O FS tem
    ruido de uns 0,001 (malha adaptativa e patamar do arc-length): o passo e grande de proposito;
  - passo HL-RF u_{k+1} = (grad . u_k - g) / |grad|^2 grad, com amortecimento (metade do passo) quando |g| mais que dobra;
  - para quando |g| < tol_g (0,001) e o beta do ponto novo muda menos que tol_beta (0,01), ou no limite de iteracoes. O
    passo |u_{k+1} - u_k| nao serve de criterio: o ruido do FS nos termos fracos o mantem em ~0,1 (mais com M grande);
  - o beta corrigido, a distancia da origem ao plano tangente no ultimo ponto (|u| + g / |grad| com u alinhado ao
    gradiente: o ponto ainda nao esta em g = 0), vai junto no form.json.
Saida: <pasta>/form.json (beta, P_f = Phi(-beta), u*, alfa e alfa^2 por modo e por campo, o historico) e o historico na
tela.

Uso: python3 -I form_kl.py <exe Cho2010> <pasta de trabalho> [M=20] [rcphi=-0.5] [h=0.5] [paralelos=2] [iteracoes=10]
     [fator=tri2] [nref=4] [srm_f0=0.8] [tol_g=0.001] [tol_beta=0.01] [gradiente=dif|modal] [nystrom=<base>]
     [metodo=srm|gim] [algoritmo=hlrf|ihlrf] [armijo=0.5] [meias=6]
algoritmo=ihlrf: o HL-RF melhorado (Zhang e Der Kiureghian 1995): a direcao d e a do passo HL-RF, d = u_HLRF - u, e o
tamanho lambda = 1, 1/2, 1/4, ... (ate meias= metades) pela regra de Armijo na funcao de merito
m(u) = |u|^2 / 2 + c |g(u)|, com c = 2 |u| / |grad g| + 10 (c > |u| / |grad g| garante que d desce em m): aceita o
primeiro lambda com m(u + lambda d) - m(u) <= armijo lambda grad m . d. Cada lambda tentado e uma analise (com o
gradiente modal, a do ponto aceito ja da o gradiente da iteracao seguinte); sem o amortecimento do HL-RF.
prog=rebaixamento: o rebaixamento do Ceron (SlopeDrawdownMonteCarlo, Projects/SlopeRandomFieldRefinement), com 3 campos
(c', phi' e k: 3M variaveis) e g = Gamma - 1 (GIM); as opcoes rcphi, srm_f0, metodo e nystrom nao valem, e transiente=
escolhe a malha do transiente. O gradiente modal dele da 0 nos termos de k (a analise limite congela as pressoes).
As opcoes do rebaixamento contorno= (reservatorio|ceron), agua= (unsat|sat), trebaix= (h), tfim= (h; -1: = trebaix),
dt= (h) e hw= (h_w/H) e as dos casos das Secoes 6.2 a 6.4 do artigo, covc=, covphi=, covk= (CoV de c', phi' e k_v),
s= (escala dos comprimentos de autocorrelacao: 20 s e 2 s m) e alfa= (k_h/k_v), tem os padroes do
SlopeDrawdownMonteCarlo (os de antes) e so vao para o executavel quando diferentes do padrao; vao no form.json, e o
line_sampling.py usa as mesmas. Ex.: a agua do artigo (percolacao permanente depois do rebaixamento instantaneo):
contorno=ceron agua=sat trebaix=0.0001 tfim=20 dt=5.
gradiente=modal: uma analise por iteracao, com o gradiente pela projecao modal da sensibilidade de limite superior
(Cho2010 gradiente=M; no colapso do SRM ou do GIM, conforme o metodo), em vez das 4M + 1 das diferencas centrais.
prog=chuva: a chuva (SlopeRainMonteCarlo, Projects/SlopeRandomFieldRefinement), com 3 campos (c', phi' e k: 3M
variaveis) e g = (o menor FS dos instantes tempos=) - 1, a primeira passagem; o gradiente modal e o do instante do menor
FS (o executavel com gradiente=M), e os termos de k saem 0 (o fator congela as pressoes). Aqui so M=, h=, paralelos=,
iteracoes=, tol_g=, tol_beta=, gradiente= e threads= sao do FORM: todas as outras chave=valor vao para o executavel como
estao (ex.: a geometria, o solo, a chuva, tempos=, nystrom=, rcphi=, pressoes=; os padroes sao os do executavel), e
ficam no form.json para o line_sampling.py. Um instante sem fator por SRM sem equilibrio (FS_min = 0) e falha (FS =
0,005); outro problema numa analise (o transiente que parou) e erro.
  reuso=1 (prog=chuva com gradiente=modal; padrao 0, o de antes): o transiente calculado uma vez so. O gradiente modal
  da 0 nos termos de k, entao o u do FORM nunca sai de xi_k = 0 e todas as analises tem o mesmo campo de k, logo o mesmo
  transiente: a 1a analise o grava (grava_pressoes=<pasta>/pressoes_form.txt) e as outras o leem (pressoes=), so com o
  fator. Precisa de um executavel que aceite grava_pressoes= com k aleatorio e os xi dados e ponha o campo de k na chave
  do arquivo (o patch do campo_k no montecarlo.cpp): uma analise barata confere isso antes (sonda_reuso, tempos=0
  fator=tri0 nref=0). Sem pressoes=, grava_pressoes= e so_transiente= nas opcoes (com k fixo o transiente ja vem de um
  arquivo). O reuso entra no form.json ("reuso": 1); o FS lido e o calculado sao iguais bit a bit (as pressoes vao no
  arquivo com 17 algarismos).
prog=coluna: a coluna 1D de Cho (2014) (SlopeRainColumn, Projects/SlopeRandomFieldRefinement/coluna.cpp), com 1
campo (k: M variaveis) e as opcoes do executavel como na chuva (as outras chave=valor vao para ele; as de cada analise,
perfis=, grava_xi=, so_campo=, para_na_falha= e localiza=, sao recusadas), so com gradiente=dif, e a falha:
  falha=fs (padrao) T= (18): g = FS_passagem(T) - 1, o menor FS ate T (primeira passagem; a analise roda com tempos=T);
  falha=tempo T=: g = t_passagem / T - 1, com para_na_falha=1 chuva=2T tempos=2T acrescentados: a analise para na
     travessia; sem travessia ate 2T o fator fica censurado em 2 (g = 1).
Um transiente que parou (situacao transiente) so vale se a travessia veio antes da parada (t_passagem >= 0); senao e
erro (a analise censurada, sem o sinal de g).
threads (0): PZ_THREADS das analises (0: um processo serial por analise; no FORM modal, so uma analise por iteracao,
mais threads aceleram)."""
import concurrent.futures as cf
import csv
import json
import math
import os
import subprocess
import sys

# as opcoes do rebaixamento passadas ao SlopeDrawdownMonteCarlo, com os padroes dele (so as diferentes vao na linha)
REBAIX_PADRAO = {"contorno": "reservatorio", "agua": "unsat", "trebaix": 0.1, "tfim": -1.0, "dt": 0.1, "hw": 1.0,
                 "covc": 0.3, "covphi": 0.1, "covk": 0.6, "s": 1.0, "alfa": 1.0}
# prog=chuva: as opcoes do FORM (as outras vao para o SlopeRainMonteCarlo) e as que o FORM poe na linha de comando
CHUVA_FORM = {"M", "h", "paralelos", "iteracoes", "tol_g", "tol_beta", "gradiente", "prog", "threads", "algoritmo",
              "armijo", "meias"}  # e o reuso=, tratado a parte (so na chuva)
CHUVA_RESERVADAS = {"n", "primeira", "xi", "saida", "gradiente", "det", "figuras"}
# reuso=1 (prog=chuva): as opcoes do executavel que ja tratam do transiente (com elas o reuso nao vale) e o sinal, na
# chave do arquivo do grava_pressoes=, de que o executavel guarda o campo de k
TRANSIENTE_ARQ = ("pressoes=", "grava_pressoes=", "so_transiente=")
CHAVE_CAMPO_K = "campo_k="
# prog=coluna: as do FORM (as da chuva mais a falha e o T) e as reservadas a mais (as de cada analise)
COLUNA_FORM = CHUVA_FORM | {"falha", "T"}
COLUNA_RESERVADAS = {"perfis", "grava_xi", "so_campo", "para_na_falha", "localiza"}
COLUNA_CENSURA = 2.  # falha=tempo: o fator sem travessia ate 2T
SITUACAO_SEM_EQUILIBRIO = 0.005  # FS (ou Gamma) de um instante sem fator: o SRM (ou o GIM) sem equilibrio, falha


def parametros(argv):
    p = {"M": 20, "rcphi": -0.5, "h": 0.5, "paralelos": 2, "iteracoes": 10, "fator": "tri2", "nref": 4, "srm_f0": 0.8,
         "tol_g": 0.001, "tol_beta": 0.01, "gradiente": "dif", "nystrom": "", "metodo": "srm", "prog": "cho",
         "transiente": "tri3", "threads": 0, "algoritmo": "hlrf", "armijo": 0.5, "meias": 6, **REBAIX_PADRAO}
    if "prog=coluna" in argv:  # como a chuva, com a falha e o T
        p.update({"falha": "fs", "T": 18.})
        p["exe"] = []
        for a in argv:
            k, v = a.split("=", 1)
            if k in COLUNA_FORM:
                p[k] = type(p[k])(v)
            elif k in CHUVA_RESERVADAS | COLUNA_RESERVADAS:
                sys.exit(f"prog=coluna: {k}= e do FORM ou da linha de comando que ele monta")
            else:
                p["exe"] = [x for x in p["exe"] if not x.startswith(k + "=")] + [a]
                if k == "nystrom":
                    p["nystrom"] = v
        if p["falha"] not in ("fs", "tempo"):
            sys.exit("falha: fs ou tempo")
        if not p["T"] > 0.:
            sys.exit("T > 0")
        if p["gradiente"] != "dif":
            sys.exit("prog=coluna: so gradiente=dif (o SlopeRainColumn nao tem o gradiente modal)")
        if p["algoritmo"] not in ("hlrf", "ihlrf"):
            sys.exit("algoritmo: hlrf ou ihlrf")
        return p
    if "prog=chuva" in argv:  # as opcoes do FORM e, o resto, as do executavel, na ordem dada
        p["exe"] = []
        for a in argv:
            k, v = a.split("=", 1)
            if k == "reuso":  # so quando dado: o form.json sem ele fica como antes
                p["reuso"] = int(v)
            elif k in CHUVA_FORM:
                p[k] = type(p[k])(v)
            elif k in CHUVA_RESERVADAS:
                sys.exit(f"prog=chuva: {k}= e do FORM ou da linha de comando que ele monta")
            else:
                p["exe"] = [x for x in p["exe"] if not x.startswith(k + "=")] + [a]
                if k == "metodo":
                    p["metodo"] = v
                if k == "nystrom":
                    p["nystrom"] = v
        if p["metodo"] not in ("srm", "gim"):
            sys.exit("metodo: srm ou gim")
        if p["algoritmo"] not in ("hlrf", "ihlrf"):
            sys.exit("algoritmo: hlrf ou ihlrf")
        if p.get("reuso", 0):
            if p["gradiente"] != "modal":
                sys.exit("reuso=1: so com gradiente=modal (nas diferencas centrais os xi de k mudam: outro transiente)")
            if any(x.startswith(TRANSIENTE_ARQ) for x in p["exe"]):
                sys.exit("reuso=1: sem pressoes=, grava_pressoes= e so_transiente= (com k fixo o transiente ja vem de um "
                         "arquivo)")
        return p
    for a in argv:
        k, v = a.split("=", 1)
        if k not in p:
            sys.exit(f"opcao desconhecida: {a}")
        p[k] = type(p[k])(v)
    if p["algoritmo"] not in ("hlrf", "ihlrf"):
        sys.exit("algoritmo: hlrf ou ihlrf")
    if p["prog"] not in ("cho", "rebaixamento"):
        sys.exit("prog=cho, prog=rebaixamento, prog=chuva ou prog=coluna")
    if p["prog"] == "rebaixamento":
        p["metodo"] = "gim"  # o fator do rebaixamento e o Gamma do GIM
    return p


def campos(p):
    """os nomes dos campos (a ordem dos xi: os M de cada um)"""
    if p.get("prog") == "coluna":
        return ["k"]
    return ["c", "phi", "k"] if p["prog"] in ("rebaixamento", "chuva") else ["c", "phi"]


def ambiente(p):
    """o ambiente de uma analise: um processo serial (ou threads= do PZ), sem threads no MKL"""
    return dict(os.environ, PZ_THREADS=str(p.get("threads", 0)), OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")


def le_fator(arq, p):
    """o fator de uma analise (o <rotulo>_0_0.csv): na chuva, o menor FS dos instantes (um instante sem equilibrio, com
    situacao fator e FS_min = 0, e falha: SITUACAO_SEM_EQUILIBRIO), e o instante dele; nos outros, a linha unica"""
    with open(arq, newline="") as f:
        linhas = list(csv.DictReader(f))
    if not linhas:
        raise RuntimeError(f"{arq}: vazio")
    if p["prog"] == "coluna":
        return fator_coluna(arq, linhas, p)
    if p["prog"] != "chuva":
        if linhas[0]["situacao"] != "ok":
            raise RuntimeError(f"{arq}: situacao {linhas[0]['situacao']} (FS_min {linhas[0]['FS_min']})")
        return float(linhas[0]["FS_min"]), float(linhas[0]["t_h"])
    melhor = None
    for r in linhas:
        fs = float(r["FS_min"])
        if r["situacao"] == "fator" and fs == 0.:
            fs = SITUACAO_SEM_EQUILIBRIO
        elif r["situacao"] != "ok":
            raise RuntimeError(f"{arq}: situacao {r['situacao']} em t = {r['t_h']} h (FS_min {r['FS_min']})")
        if melhor is None or fs < melhor[0]:
            melhor = (fs, float(r["t_h"]))
    return melhor


def fator_coluna(arq, linhas, p):
    """prog=coluna: o fator da linha do instante (T na falha=fs, 2T na falha=tempo) e esse instante. Um transiente que
    parou antes so vale com a travessia ja localizada (t_passagem >= 0); senao e erro"""
    T = float(p["T"])
    tf = T if p["falha"] == "fs" else 2. * T
    r = next((x for x in linhas if abs(float(x["t_h"]) - tf) <= 1e-9 * tf), None)
    if r is None:
        raise RuntimeError(f"{arq}: sem a linha de t = {tf} h")
    tp = float(r["t_passagem"])
    if r["situacao"] not in ("ok", "falha") and tp < 0.:
        raise RuntimeError(f"{arq}: situacao {r['situacao']} sem travessia ate t_parou = {r['t_parou']} h (censurada)")
    if p["falha"] == "fs":
        return float(r["FS_passagem"]), tf
    return (tp / T if tp >= 0. else COLUNA_CENSURA), tf


def comando(exe, p, arq, rotulo):
    """a linha de comando de uma analise nos xi do arquivo"""
    if p["prog"] == "coluna":
        T = float(p["T"])
        fim = [f"tempos={T:g}"] if p["falha"] == "fs" else ["para_na_falha=1", f"chuva={2 * T:g}", f"tempos={2 * T:g}"]
        return [exe, "n=1", "primeira=0", f"xi={arq}", f"M={p['M']}"] + list(p["exe"]) + fim + [f"saida={rotulo}"]
    if p["prog"] == "chuva":
        return [exe, "n=1", "primeira=0", f"xi={arq}", f"M={p['M']}"] + list(p["exe"]) + [f"saida={rotulo}"]
    if p["prog"] == "rebaixamento":
        extra = [f"{k}={p[k]}" for k, v in REBAIX_PADRAO.items() if k in p and p[k] != v]
        return [exe, "n=1", "primeira=0", f"xi={arq}", f"transiente={p['transiente']}", f"fator={p['fator']}",
                f"nref={p['nref']}"] + extra + [f"saida={rotulo}"]
    return [exe, "n=1", "primeira=0", f"xi={arq}", f"M={p['M']}", f"rcphi={p['rcphi']}", f"fator={p['fator']}",
            f"nref={p['nref']}", f"srm_f0={p['srm_f0']}", f"metodo={p['metodo']}", f"saida={rotulo}"] + \
        ([f"nystrom={p['nystrom']}"] if p["nystrom"] else [])


def sonda_reuso(exe, p, pasta):
    """reuso=1 (prog=chuva): o executavel aceita grava_pressoes= com k aleatorio e os xi dados (n=1), e a chave do arquivo
    tem o campo de k (CHAVE_CAMPO_K)? Uma analise barata com as opcoes do caso, M = p["M"]: tempos=0 (so o estado
    inicial do transiente), o fator na tri0 sem refinamento, xi de c' e phi' nulos e um xi de k nao nulo. Se nao, sai do
    programa com a explicacao (antes de qualquer analise de verdade)"""
    os.makedirs(pasta, exist_ok=True)
    rot, M = "sonda_reuso", int(p["M"])
    arq_xi, arq_pw = os.path.join(pasta, rot + ".xi"), os.path.join(pasta, rot + "_pw.txt")
    q = dict(p, exe=[x for x in p["exe"] if not x.startswith(("tempos=", "adapta=", "fator=", "nref="))] +
             ["tempos=0", "fator=tri0", "nref=0"])
    with open(arq_xi, "w") as f:
        f.write("\n".join(["0.0"] * (2 * M) + ["0.5"] + ["0.0"] * (M - 1)) + "\n")
    r = subprocess.run(comando(exe, q, arq_xi, rot) + [f"grava_pressoes={arq_pw}"], stdout=subprocess.DEVNULL,
                       stderr=subprocess.PIPE, text=True, env=ambiente(q), cwd=pasta)
    chave = ""
    if r.returncode == 0 and os.path.exists(arq_pw):
        with open(arq_pw) as f:
            f.readline()  # o comentario
            chave = f.readline()
    for suf in (".xi", "_pw.txt", "_pw.txt.tmp", "_0_0.csv", "_mapa_0_0.csv", "_parametros.txt"):
        try:
            os.remove(os.path.join(pasta, rot + suf))
        except FileNotFoundError:
            pass
    if r.returncode != 0:
        sys.exit("reuso=1: o executavel recusou grava_pressoes= com k aleatorio e os xi dados (" +
                 " ".join(r.stderr.split())[-300:] + "): falta o patch do campo de k (campo_k) no montecarlo.cpp; rode "
                 "sem reuso")
    if not chave.startswith("chave ") or CHAVE_CAMPO_K not in chave:
        sys.exit(f"reuso=1: a chave do arquivo do grava_pressoes= nao tem o campo de k ({CHAVE_CAMPO_K}): um arquivo "
                 "poderia servir a outro campo de k; rode sem reuso")


class Problema:
    def __init__(self, exe, pasta, p):
        self.exe, self.pasta, self.p = exe, pasta, p
        self.n = 0
        os.makedirs(pasta, exist_ok=True)
        # reuso=1: o transiente do FORM (xi_k = 0 em todas as analises) num arquivo, gravado pela primeira analise
        self.pw = os.path.join(pasta, "pressoes_form.txt") if p.get("reuso", 0) else None
        self.pw_pronto = False

    def transiente(self, u):
        """reuso=1: as opcoes do transiente desta analise (grava na primeira, le nas outras), conferindo xi_k = 0"""
        if self.pw is None:
            return []
        if any(x != 0. for x in u[2 * self.p["M"]:]):
            raise RuntimeError("reuso=1: xi de k nao nulo no FORM (o gradiente modal da 0 em k): seria outro transiente")
        return [f"pressoes={self.pw}"] if self.pw_pronto else [f"grava_pressoes={self.pw}"]

    def transiente_feito(self):
        if self.pw is not None and not self.pw_pronto and os.path.exists(self.pw):
            self.pw_pronto = True

    def fs(self, u, rotulo):
        """o fator do Cho2010 (FS ou Gamma) com os xi u (os M de c' e os M de phi')"""
        arq = os.path.join(self.pasta, rotulo + ".xi")
        with open(arq, "w") as f:
            f.write("\n".join(repr(x) for x in u) + "\n")
        cmd = comando(self.exe, self.p, arq, rotulo) + self.transiente(u)
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, env=ambiente(self.p),
                       cwd=self.pasta)
        fs = le_fator(os.path.join(self.pasta, rotulo + "_0_0.csv"), self.p)[0]
        self.transiente_feito()
        return fs

    def g_modal(self, u, it):
        """uma analise: o fator e o gradiente pela projecao modal (Cho2010 gradiente=M: limite superior no colapso)"""
        rotulo = f"it{it}_centro"
        arq = os.path.join(self.pasta, rotulo + ".xi")
        with open(arq, "w") as f:
            f.write("\n".join(repr(x) for x in u) + "\n")
        cmd = comando(self.exe, self.p, arq, rotulo) + [f"gradiente={self.p['M']}"] + self.transiente(u)
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, env=ambiente(self.p),
                       cwd=self.pasta)
        self.transiente_feito()
        # sem o _grad_0.txt: a mensagem clara em vez de FileNotFoundError
        fs, t = le_fator(os.path.join(self.pasta, rotulo + "_0_0.csv"), self.p)
        arq_g = os.path.join(self.pasta, rotulo + "_grad_0.txt")
        if not os.path.exists(arq_g):
            raise RuntimeError(f"{rotulo}: sem o gradiente (FS {fs} em t = {t} h)")
        v = [float(x) for x in open(arq_g).read().split("\n")[1:] if x.strip()]
        self.n += 1
        return v[0] - 1., v[1:], [v[0], t]

    def g_e_gradiente(self, u, it):
        if self.p["gradiente"] == "modal":
            return self.g_modal(u, it)
        h = self.p["h"]
        pontos = [("centro", list(u))]
        for i in range(len(u)):
            for sinal in (+1, -1):
                v = list(u)
                v[i] += sinal * h
                pontos.append((f"d{i}{'p' if sinal > 0 else 'm'}", v))
        with cf.ThreadPoolExecutor(max_workers=self.p["paralelos"]) as ex:
            fs = list(ex.map(lambda pt: self.fs(pt[1], f"it{it}_{pt[0]}"), pontos))
        self.n += len(pontos)
        g = fs[0] - 1.
        grad = [(fs[1 + 2 * i] - fs[2 + 2 * i]) / (2. * h) for i in range(len(u))]
        return g, grad, fs


def phi_normal(z):
    return 0.5 * math.erfc(-z / math.sqrt(2.))


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    exe, pasta = os.path.abspath(sys.argv[1]), os.path.abspath(sys.argv[2])
    p = parametros(sys.argv[3:])
    prob = Problema(exe, pasta, p)
    M = p["M"]
    nomes = campos(p)
    nc = len(nomes)
    u = [0.] * (nc * M)
    hist = []
    titulo = ("rebaixamento do Ceron (2025), " + (" ".join(comando("", p, "", "")[4:-1]) or "padrao")) \
        if p["prog"] == "rebaixamento" else \
        ("chuva (SlopeRainMonteCarlo), " + (" ".join(p["exe"]) or "padrao") + ", g = menor FS dos instantes - 1") \
        if p["prog"] == "chuva" else \
        ("SlopeRainColumn (coluna de Cho 2014), " + (" ".join(p["exe"]) or "padrao") +
         (f", g = FS_passagem({p['T']:g} h) - 1" if p["prog"] == "coluna" and p["falha"] == "fs" else
          f", g = t_passagem / {p['T']:g} h - 1 (censura em 2T)")) \
        if p["prog"] == "coluna" else \
        f"Cho (2010), Exemplo 2, r(c, phi) = {p['rcphi']}"
    tipo = "modal (1 avaliacao por iteracao)" if p["gradiente"] == "modal" else \
        f"por diferencas centrais (h = {p['h']}, {2 * nc * M + 1} avaliacoes por iteracao)"
    fator = "FS" if p["prog"] == "coluna" else p["metodo"].upper()
    print(f"FORM do {titulo}, {fator}: KL com M = {M} termos por campo ({nc * M} variaveis), gradiente "
          f"{tipo}, {p['paralelos']} em paralelo", flush=True)
    if p.get("reuso", 0):  # o executavel tem de aceitar o transiente gravado com k aleatorio (sai daqui se nao)
        sonda_reuso(exe, p, pasta)
        print(f"  reuso do transiente: a 1a analise grava {prob.pw}, as outras so calculam o fator", flush=True)
    print(f"  it   g = {'Gamma' if p['metodo'] == 'gim' else 'FS'} - 1    beta     |grad|   maiores alfa^2 (campo, termo)",
          flush=True)
    gAnterior = None
    convergiu = False
    aceito = None  # iHLRF com gradiente modal: (g, grad, fs) do ponto aceito na busca, ja com o gradiente
    for it in range(p["iteracoes"]):
        if aceito is not None:
            g, grad, fs = aceito
        else:
            g, grad, fs = prob.g_e_gradiente(u, it)
        aceito = None
        beta = math.sqrt(sum(x * x for x in u))
        n2 = sum(x * x for x in grad)
        if n2 == 0. and p["prog"] == "coluna":
            sys.exit(f"gradiente nulo na iteracao {it} (g = {g}): na falha=tempo, todas as analises sem travessia ate "
                     "2T? (um T maior)")
        alfa2 = [x * x / n2 for x in grad] if n2 > 0 else [0.] * len(grad)
        top = sorted(range(len(grad)), key=lambda i: -alfa2[i])[:4]
        print(f"  {it:2d}   {g:+.5f}    {beta:.4f}   {math.sqrt(n2):.4f}   " +
              ", ".join(f"{nomes[i // M]}{i % M}: {alfa2[i]:.2f}" for i in top), flush=True)
        hist.append({"it": it, "u": u, "g": g, "grad": grad, "beta": beta, "fs_pontos": fs})
        with open(os.path.join(pasta, "form_parcial.json"), "w") as f:
            json.dump({"parametros": p, "historico": hist, "avaliacoes": prob.n}, f, indent=1)
        gu = sum(gi * ui for gi, ui in zip(grad, u))
        novo = [(gu - g) / n2 * gi for gi in grad]
        beta_novo = math.sqrt(sum(x * x for x in novo))
        if it > 0 and abs(g) < p["tol_g"] and abs(beta_novo - beta) < p["tol_beta"]:
            convergiu = True
            break
        if p["algoritmo"] == "ihlrf":
            d = [ni - ui for ni, ui in zip(novo, u)]
            c = 2. * beta / math.sqrt(n2) + 10.
            sg = 1. if g >= 0. else -1.
            m0 = 0.5 * beta * beta + c * abs(g)
            dm = sum((ui + c * sg * gi) * di for ui, gi, di in zip(u, grad, d))  # grad m . d (< 0)
            lam, busca = 1., []
            for j in range(p["meias"] + 1):
                ut = [ui + lam * di for ui, di in zip(u, d)]
                rot = f"{it}l{j}"
                if p["gradiente"] == "modal":
                    gt, gradt, fst = prob.g_modal(ut, rot)
                else:
                    gt = prob.fs(ut, f"it{rot}_centro") - 1.
                    prob.n += 1
                    gradt = fst = None
                mt = 0.5 * sum(x * x for x in ut) + c * abs(gt)
                ok = mt - m0 <= p["armijo"] * lam * dm
                busca.append({"lambda": lam, "g": gt, "merito": mt, "aceito": ok})
                print(f"       busca: lambda {lam:.4f}  g {gt:+.5f}  merito {mt:.5f} (antes {m0:.5f}, "
                      f"Armijo {m0 + p['armijo'] * lam * dm:.5f}) {'aceito' if ok else ''}", flush=True)
                if ok or j == p["meias"]:
                    if not ok:
                        print("       busca: nenhum lambda passou na regra de Armijo; fica o menor", flush=True)
                    if gradt is not None:
                        aceito = (gt, gradt, fst)
                    break
                lam *= 0.5
            hist[-1].update({"c": c, "merito": m0, "busca": busca})
            u = ut
            continue
        if gAnterior is not None and abs(g) > 2. * abs(gAnterior) and it > 1:
            novo = [ui + 0.5 * (ni - ui) for ui, ni in zip(u, novo)]  # amortece: metade do passo
        gAnterior = g
        u = novo
    beta = math.sqrt(sum(x * x for x in u))
    n = math.sqrt(sum(x * x for x in grad))
    alfa = [-x / n for x in grad]
    a2 = [x * x for x in alfa]
    beta_corr = beta_novo  # distancia ao plano tangente no ultimo ponto avaliado (o proximo passo do HL-RF)
    if hist[0]["g"] < 0.:  # a origem (as medias) ja esta na falha: beta negativo e P_f = Phi(-beta) > 1/2
        beta, beta_corr = -beta, -beta_corr
        print("  g(0) < 0: as medias ja falham; beta negativo (o line_sampling.py inverte a direcao)", flush=True)
    res = {"parametros": p, "beta": beta, "P_f": phi_normal(-beta), "beta_corrigido": beta_corr,
           "P_f_corrigida": phi_normal(-beta_corr), "convergiu": convergiu, "iteracoes": len(hist), "u": u, "alfa": alfa, "alfa2": a2,
           "alfa2_c": sum(a2[:M]), "alfa2_phi": sum(a2[M:2 * M]), "alfa2_k": sum(a2[2 * M:]), "avaliacoes": prob.n,
           "historico": hist}
    with open(os.path.join(pasta, "form.json"), "w") as f:
        json.dump(res, f, indent=1)
    print(f"\n  beta = {beta:.4f}  P_f = Phi(-beta) = {100 * res['P_f']:.3f}%  ({prob.n} avaliacoes do FS)")
    print(f"  beta corrigido = {beta_corr:.4f}  P_f = {100 * res['P_f_corrigida']:.3f}%")
    print("  convergiu" if convergiu else f"  NAO convergiu em {p['iteracoes']} iteracoes")
    if p["prog"] == "coluna":
        print("  importancia nos termos: " + ", ".join(f"1 a {m}: {sum(a2[:m]):.3f}" for m in (10, 20, 50) if m <= M))
    else:
        print(f"  importancia por campo: c' {res['alfa2_c']:.3f}, phi' {res['alfa2_phi']:.3f}"
              + (f", k {res['alfa2_k']:.3f}" if nc == 3 else ""))


if __name__ == "__main__":
    main()
