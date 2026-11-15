#!/usr/bin/env python3
"""figuras_etapa0.py: as figuras, as tabelas e os numeros do relatorio_etapa0.tex, o relatorio da etapa 0 do plano de
trabalho (../README.md): fechar o que ja rodava antes das etapas no fermi.
  A) a chuva com amostragem por importancia (AI) no fermi, completa (CAMPANHA=chuva_ai, ~/projects/runs/fermi_chuva_ai):
     o solo f50 do Rahardjo a 31,05 mm/h, o piloto de Monte Carlo direto e a AI na direcao da banda;
  B) o Cho (2010), Exemplo 2, no PC (~/projects/runs/cho2010), uma previa com o SRM antigo: o deterministico, o Monte
     Carlo direto e a AI com r(c, phi) = -0,5 e -0,7, o FORM com M = 10 e 20 termos por campo e a AI com o ponto de
     projeto do FORM. O PC parou em 15/11 (plano, versao 2); os casos rodam de novo no fermi (campanha etapa1);
  C) o contexto ja validado, citado sem recalcular: as macros do numeros_rebaixamento.tex (relatorio_rebaixamento.pdf) e
     do numeros_solos.tex (relatorio_solos.pdf), em Projects/SlopeRandomFieldRefinement/relatorio.
So le arquivos: nao roda nenhum programa e nao escreve nas pastas das rodadas. As contas da AI sao as do
posproc/amostragem_importancia.py: peso w = exp(log_peso), P_f = soma(w I) / N e CoV = s_q / (raiz(N) P_f).

Escreve (nesta pasta):
  figs/fig_e0_*.pdf     as figuras
  numeros_etapa0.tex    um \\newcommand para cada numero do texto e das tabelas e os \\newif das rodadas que ainda rodam
  resultados/cho2010/   copias dos resumos pequenos do Cho2010 que ja existem: resultado.txt, resultado_form.txt e o
                        form.log e o form.json de cada FORM (o form_parcial.json, se o FORM nao terminou), e o
                        resumo_previa.txt do direto do PC (o resultado.txt da chuva_ai ja esta em resultados/chuva_ai)
e imprime cada numero com o arquivo de onde sai. Rodado de novo (atualiza.sh), o relatorio pega os numeros novos e
troca "pendente" pelos resultados quando as rodadas acabam.

Nao ha numero de resultado digitado aqui: os do artigo do Cho vem do Projects/Cho2010/README.md, os do plano do
../README.md e do cluster/casos_etapa1.txt e os do piloto antigo da chuva do build (mc/gim_I/piloto_srm_I31.05_0_7.csv).
Escolhidos aqui: as faixas de FS em volta de 1 (FAIXAS) e o CoV alvo da conta do direto da chuva (COV_ALVO).

Uso: python3 -I figuras_etapa0.py [chuva_ai baixada] [cho2010] [build do SlopeRandomFieldRefinement]
     (padrao: ~/projects/runs/fermi_chuva_ai, ~/projects/runs/cho2010 e
     ../../../../neopz-master-build-release/Projects/SlopeRandomFieldRefinement)"""

import csv
import datetime as dt
import glob
import heapq
import json
import math
import os
import re
import shutil
import sys

import matplotlib
from scipy.stats import beta as dist_beta

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.ticker  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from scipy.optimize import brentq  # noqa: E402

AQUI = os.path.dirname(os.path.abspath(__file__))  # PlanoDeTrabalho/ETAPA0
PLANO = os.path.dirname(AQUI)
NEOPZ = os.path.dirname(os.path.dirname(PLANO))
SRF = os.path.join(NEOPZ, "Projects", "SlopeRandomFieldRefinement")
CHO_SRC = os.path.join(NEOPZ, "Projects", "Cho2010")
REL = os.path.join(SRF, "relatorio")
RUNS = os.path.expanduser("~/projects/runs")
CHUVA_PADRAO = os.path.join(RUNS, "fermi_chuva_ai")
CHO_PADRAO = os.path.join(RUNS, "cho2010")
BUILD_PADRAO = os.path.join(os.path.dirname(NEOPZ), "neopz-master-build-release", "Projects",
                            "SlopeRandomFieldRefinement")
DATA_RELATORIO = dt.date(2026, 11, 15)
AGORA = dt.datetime.now()

# ---- estilo (o do figuras_rebaixamento.py)
PALETA = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
TINTA, TINTA2, MUDO, GRADE, SOLO = "#0b0b0b", "#52514e", "#898781", "#dddcd8", "#f3f2ee"
COR_DIRETO, COR_AI, COR_TERCEIRA = PALETA[0], PALETA[1], PALETA[2]
plt.rcParams.update({
    "font.size": 8.5, "axes.labelsize": 8.5, "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "legend.fontsize": 7,
    "axes.edgecolor": TINTA2, "axes.labelcolor": TINTA, "xtick.color": TINTA2, "ytick.color": TINTA2,
    "axes.grid": True, "grid.color": GRADE, "grid.linewidth": 0.5, "lines.linewidth": 1.4, "legend.frameon": False,
    "savefig.dpi": 200, "font.family": "DejaVu Sans", "mathtext.fontset": "dejavusans",
})
LARGURA = 6.5  # polegadas: a largura do texto do relatorio (A4 com margens de 2,2 cm)

NUMEROS = []  # (macro, valor impresso, valor no .tex, arquivo)
FLAGS = {}    # nome do \newif -> verdadeiro ou falso
# as faixas de FS em volta de 1 para o buraco do SRM (a mesma ideia do figuras_solos.py, com faixas de 0,01)
FAIXAS = [(0.98, 0.99), (0.99, 1.00), (1.00, 1.01), (1.01, 1.02)]
COV_ALVO = 0.10  # o CoV(P_f) da conta de quantas realizacoes diretas a chuva pede
MIN_WILSON = [20]  # com menos falhas, o IC exato (Clopper-Pearson): o plano (criterios gerais) diz ~20


# ---------------------------------------------------------------- numeros
def v(x, casas=4):
    """numero com virgula decimal (sem o sinal de -0,0)"""
    t = "%.*f" % (casas, x)
    if t.startswith("-") and float(t) == 0.:
        t = t[1:]
    return t.replace(".", ",")


def pct(x):
    """uma probabilidade em %, com as casas pelo tamanho: 6,22; 0,37; 0,014"""
    return v(x, 2) if x >= 0.1 or x == 0. else v(x, 3)


def mil(n):
    """inteiro com o espaco fino do LaTeX a partir de 10 000"""
    n = int(round(n))
    return "{:,}".format(n).replace(",", r"\,") if n >= 10000 else "%d" % n


def num(macro, texto, fonte, latex=None):
    """registra um numero do texto: \\macro no .tex (virgula decimal como {,}, o sinal de menos em modo matematico) e o
    valor impresso com o arquivo de onde sai. latex: o valor no .tex, quando nao e um numero simples"""
    assert re.fullmatch(r"[A-Za-z]+", macro), macro
    assert macro not in {m for m, _, _, _ in NUMEROS}, "macro repetida: " + macro
    if latex is None:
        latex = re.sub(r"(?<!\\),", "{,}", texto)
        if re.match(r"-\d", latex):
            latex = r"\ensuremath{-}" + latex[1:]
    NUMEROS.append((macro, texto.replace(r"\,", " "), latex, fonte))
    return texto


def txt(macro, texto, fonte):
    """um texto que depende dos dados (uma conclusao, um horario, linhas de tabela): vai para o .tex como esta"""
    return num(macro, texto, fonte, latex=texto)


def flag(nome, valor, fonte):
    FLAGS[nome] = bool(valor)
    NUMEROS.append(("if" + nome, "verdadeiro" if valor else "falso", None, fonte))


def tex(x, casas):
    """o numero no .tex, para montar linhas de tabela e frases no script"""
    t = v(x, casas).replace(",", "{,}")
    return r"\ensuremath{-}" + t[1:] if t.startswith("-") else t


def rel(arq, base):
    return os.path.relpath(arq, base)


def hora(t):
    """'21h05', com o dia quando nao e o do relatorio"""
    s = "%dh%02d" % (t.hour, t.minute)
    return s if t.date() == DATA_RELATORIO else s + " de %d/%d" % (t.day, t.month)


def hora_meia(t):
    """arredondada para a meia hora (uma estimativa)"""
    m = int(round((t.hour * 60 + t.minute + t.second / 60.) / 30.) * 30)
    return hora(dt.datetime.combine(t.date(), dt.time()) + dt.timedelta(minutes=m))


def duracao(td):
    m = int(round(td.total_seconds() / 60.))
    return "%d~h %02d~min" % (m // 60, m % 60) if m >= 60 else "%d~min" % m


UNID = ["Zero", "Um", "Dois", "Tres", "Quatro", "Cinco", "Seis", "Sete", "Oito", "Nove", "Dez", "Onze", "Doze", "Treze",
        "Quatorze", "Quinze", "Dezesseis", "Dezessete", "Dezoito", "Dezenove"]
DEZENAS = ["", "", "Vinte", "Trinta", "Quarenta", "Cinquenta", "Sessenta", "Setenta", "Oitenta", "Noventa"]


def extenso(n):
    """o inteiro n (0 a 99) por extenso, para o nome de uma macro"""
    n = int(round(n))
    if n < 20:
        return UNID[n]
    d, u = divmod(n, 10)
    return DEZENAS[d] + (UNID[u] if u else "")


def falhas(n):
    """'nenhuma falha', '1 falha', '13 falhas'"""
    return "nenhuma falha" if n == 0 else "%d falha%s" % (n, "" if n == 1 else "s")


def norm_cdf(x):
    return 0.5 * math.erfc(-x / math.sqrt(2.))


def norm_pdf(x):
    return math.exp(-0.5 * x * x) / math.sqrt(2. * math.pi)


def norm_inv(p):
    """beta com Phi(-beta) = p"""
    return brentq(lambda x: norm_cdf(-x) - p, -10., 40.)


def wilson(k, n, z=1.959964):
    p = k / n
    c = (p + z * z / (2 * n)) / (1 + z * z / n)
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return max(c - h, 0.), c + h


def clopper(k, n, a=0.05):
    """o IC exato (Clopper-Pearson) de uma proporcao k / n"""
    lo = 0. if k == 0 else float(dist_beta.ppf(a / 2., k, n - k + 1))
    hi = 1. if k == n else float(dist_beta.ppf(1. - a / 2., k + 1, n - k))
    return lo, hi


def intervalo(k, n):
    """o IC de 95% do direto: o exato (Clopper-Pearson) com menos de ~20 falhas, como pede o plano, e o de Wilson com
    mais; devolve as pontas e o nome"""
    return clopper(k, n) + ("exato",) if k < MIN_WILSON[0] else wilson(k, n) + ("de Wilson",)


def le_plano(arq, padrao, nome):
    """um trecho do plano (README.md, que muda): o match ou None, com um aviso"""
    m = re.search(padrao, texto(arq), re.S)
    if not m:
        print("AVISO: nao achei %s em %s" % (nome, arq))
    return m


def ganho_ai(p):
    """estado limite linear com P_f = p: o deslocamento otimo da entropia cruzada, m = phi(beta) / p, na direcao da
    falha, e Var(AI) / Var(direto) por realizacao, com E_h[(I w)^2] = exp(m^2) Phi(-beta - m) (relatorio_chuva.tex,
    Secao 7)"""
    b = norm_inv(p)
    m = norm_pdf(b) / p
    return b, m, (math.exp(m * m) * norm_cdf(-b - m) - p * p) / (p * (1. - p))


def ganho_desloc(p, a):
    """o mesmo estado limite linear com um deslocamento de tamanho a na direcao exata: Var(direto) / Var(AI)"""
    b = norm_inv(p)
    return p * (1. - p) / (math.exp(a * a) * norm_cdf(-b - a) - p * p)


def lognormal(x):
    """a P_f do ajuste lognormal, Phi(-media(ln x) / desvio(ln x)), como no posproc/resumo_tempos.py"""
    ln = np.log(np.asarray(x, float))
    return norm_cdf(-ln.mean() / ln.std(ddof=1)), ln.mean(), ln.std(ddof=1)


# ---------------------------------------------------------------- leitura
def texto(arq):
    with open(arq, errors="replace") as f:
        return f.read()


def le_fonte(arq, padrao, nome, flags=re.S):
    m = re.search(padrao, texto(arq), flags)
    if not m:
        sys.exit("nao achei %s em %s" % (nome, arq))
    return m


def le_chaves(arq):
    return dict(l.strip().split("=", 1) for l in open(arq) if "=" in l and not l.startswith("#"))


def le_casos(arq):
    """as linhas de um casos.txt: id -> (N, opcoes)"""
    out = {}
    for linha in texto(arq).splitlines():
        if linha.strip() and not linha.startswith("#"):
            c = linha.split()
            out[c[0]] = (int(c[2]), dict(o.split("=", 1) for o in c[3:] if "=" in o))
    return out


def le_macros(arq):
    """as macros de um numeros_*.tex: nome -> valor no .tex"""
    out = {}
    for l in texto(arq).splitlines():
        m = re.match(r"\\newcommand\{\\(\w+)\}\{(.*)\} %", l)
        if m:
            out[m[1]] = m[2]
    return out


def cita(macros, nome, novo, fonte):
    """um numero de outro relatorio, como esta la (citado, nao recalculado)"""
    lat = macros[nome]
    impresso = lat.replace("{,}", ",").replace(r"\ensuremath{-}", "-").replace("$", "")
    num(novo, impresso, fonte, latex=lat)
    return float(impresso.replace(",", ".")) if re.fullmatch(r"-?[\d,]+", impresso) else None


def le_rodada(pasta, varios=False):
    """as linhas completas dos mc_<a>_<b>.csv da pasta (a rodada pode estar escrevendo), como o linhas() do
    posproc/amostragem_importancia.py: sem as realizacoes com algum instante sem fator (situacao != ok ou FS <= 0).
    varios: a chuva (uma linha por instante; fs_t = {t: FS}); senao uma linha por realizacao (o Cho2010)"""
    rs, ruins = {}, set()
    for arq in sorted(glob.glob(os.path.join(pasta, "mc_*_*.csv"))):
        if not re.search(r"/mc_-?\d+_-?\d+\.csv$", arq):
            continue  # mapas
        with open(arq, newline="") as f:
            for r in csv.DictReader(f):
                try:
                    n, t, x = int(r["real"]), float(r["t_h"] or 0.), float(r["FS_min"])
                    lp = float(r["log_peso"])
                    pj = [float(r[c]) for c in ("proj_c", "proj_phi", "proj_k")]
                except (TypeError, ValueError, KeyError):
                    continue  # a linha ainda nao chegou ao fim
                if r["situacao"] != "ok" or not x > 0:
                    ruins.add(n)
                    continue
                d = rs.setdefault(n, {"fs_t": {}, "log_peso": lp, "proj": pj, "linha": r})
                d["fs_t"][t] = x
    for n in ruins:
        rs.pop(n, None)
    for d in rs.values():
        d["fs"] = min(d["fs_t"].values())
    return rs, ruins


def blocos(N, B):
    return [(a, min(a + B, N) - 1) for a in range(0, N, B)]


def feitos(pasta, N, B, ate=None):
    """os blocos prontos (com o mc_mapa) das realizacoes 0..ate-1 (todas sem ate)"""
    bl = [ab for ab in blocos(N, B) if ate is None or ab[0] < ate]
    ok = [ab for ab in bl if os.path.exists(os.path.join(pasta, "mc_mapa_%d_%d.csv" % ab))]
    return len(ok), len(bl)


RX_TAREFA = re.compile(r"^(\d{4}-\d\d-\d\d \d\d:\d\d:\d\d) (\S+) (\d+)\.\.(\d+): (inicio|fim em (\d+) s)", re.M)


def le_tarefas(pasta_logs):
    """as tarefas dos logs do cluster/tarefas.sh: (caso, a, z) -> {inicio, fim, s}"""
    t = {}
    for arq in glob.glob(os.path.join(pasta_logs, "*.txt")):
        for m in RX_TAREFA.finditer(texto(arq)):
            k = (m[2], int(m[3]), int(m[4]))
            quando = dt.datetime.strptime(m[1], "%Y-%m-%d %H:%M:%S")
            d = t.setdefault(k, {})
            if m[5] == "inicio":
                d["inicio"] = quando
            else:
                d["fim"], d["s"] = quando, int(m[6])
    return t


def le_hosts(pasta_logs):
    """os nos (hostnames) da primeira linha de cada log do tarefas.sh"""
    h = set()
    for arq in glob.glob(os.path.join(pasta_logs, "*.txt")):
        l = texto(arq).split("\n", 1)[0].split()
        if len(l) > 2:
            h.add(l[2])
    return h


RX_REAL_CHO = re.compile(r"realizacao\s+(-?\d+): FS = ([\d.]+) \((\w+)\) \| niveis [\d.;]+ \| ([\d.]+) s")


def le_tempos_cho(pasta):
    """realizacao -> (FS, s) das linhas dos mc_*.log do Cho2010"""
    t = {}
    for arq in glob.glob(os.path.join(pasta, "mc_*_*.log")):
        for m in RX_REAL_CHO.finditer(texto(arq)):
            t[int(m[1])] = (float(m[2]), float(m[4]))
    return t


# ---------------------------------------------------------------- estimadores
def estima(rs, ts=None):
    """P_f = soma(peso * falha) / N, o CoV dela, as falhas sorteadas, N, os pesos e o ESS (posproc/
    amostragem_importancia.py estima). ts: a falha e FS < 1 em algum desses instantes (todos sem ts)"""
    ns = sorted(rs)
    G = np.array([rs[n]["fs"] if ts is None else min(rs[n]["fs_t"][t] for t in ts) for n in ns])
    w = np.exp(np.array([rs[n]["log_peso"] for n in ns]))
    f = G < 1.
    q = w * f
    N = len(q)
    pf = q.mean()
    cov = q.std(ddof=1) / math.sqrt(N) / pf if pf > 0 else math.inf
    return {"N": N, "pf": pf, "cov": cov, "nf": int(f.sum()), "w_media": w.mean(), "w_se": w.std(ddof=1) / math.sqrt(N),
            "ess": w.sum() ** 2 / (w * w).sum(), "wf_media": w[f].mean() if f.any() else math.nan,
            "wf_max": w[f].max() if f.any() else math.nan, "w_max": w.max(), "G": G, "w": w, "q": q, "f": f,
            "parcela": q.max() / q.sum() if q.sum() > 0 else math.nan,
            "proj": np.array([rs[n]["proj"] for n in ns])}


def corrida(q):
    """P_f e o desvio dela com as n primeiras realizacoes, n = 1..N"""
    n = np.arange(1, len(q) + 1, dtype=float)
    m = np.cumsum(q) / n
    s2 = (np.cumsum(q * q) - n * m * m) / np.maximum(n - 1., 1.)
    return n, m, np.sqrt(np.maximum(s2, 0.) / n)


def compara(ai, p, cov_p):
    """a AI contra o direto (P_f = p, CoV cov_p): diferenca em desvios, quantas vezes menos realizacoes para o mesmo
    CoV (com a P_f da AI no CoV do direto, como o amostragem_importancia.py) e esse CoV do direto"""
    dif = (ai["pf"] - p) / math.sqrt((ai["pf"] * ai["cov"]) ** 2 + (p * cov_p) ** 2)
    cov_d = math.sqrt((1 - ai["pf"]) / (ai["N"] * ai["pf"]))
    return dif, (cov_d / ai["cov"]) ** 2, cov_d


def contagens(fs):
    fs = np.asarray(fs)
    return [int(((fs >= a) & (fs < b)).sum()) for a, b in FAIXAS]


def le_desloc(pasta):
    """os deslocamentos que os blocos da AI gravaram (mc_desloc_<a>_<b>.txt): quantos, se sao iguais e as chaves"""
    arqs = sorted(glob.glob(os.path.join(pasta, "mc_desloc_*_*.txt")))
    if not arqs:
        return 0, True, {}
    textos = {"".join(l for l in open(a) if not l.startswith("#")) for a in arqs}
    return len(arqs), len(textos) == 1, le_chaves(arqs[0])


# ---------------------------------------------------------------- desenho
def letra(ax, t):
    ax.text(-0.02, 1.02, "(%s)" % t, transform=ax.transAxes, ha="right", va="bottom", fontsize=9, color=TINTA)


class Virgula(matplotlib.ticker.ScalarFormatter):
    """os rotulos do ScalarFormatter com virgula decimal"""
    def __call__(self, x, pos=None):
        return super().__call__(x, pos).replace(".", ",")


def salva(fig, pasta, nome):
    for ax in fig.axes:  # virgula decimal nos eixos lineares
        for eixo, escala in ((ax.xaxis, ax.get_xscale()), (ax.yaxis, ax.get_yscale())):
            if escala == "linear" and isinstance(eixo.get_major_formatter(), matplotlib.ticker.ScalarFormatter):
                eixo.set_major_formatter(Virgula())
    arq = os.path.join(pasta, nome + ".pdf")
    fig.savefig(arq)
    plt.close(fig)
    print("figura: %s" % arq)


def eixo_log_virgula(ax, eixo="x"):
    fmt = matplotlib.ticker.FuncFormatter(lambda x, p: ("%g" % x).replace(".", ","))
    (ax.xaxis if eixo == "x" else ax.yaxis).set_major_formatter(fmt)
    (ax.xaxis if eixo == "x" else ax.yaxis).set_minor_formatter(matplotlib.ticker.NullFormatter())


def banda_corrida(ax, q, cor, rot, i0, escala=100.):
    n, m, se = corrida(q)
    ax.fill_between(n[i0:], escala * np.maximum(m - 1.96 * se, 0.)[i0:], escala * (m + 1.96 * se)[i0:], color=cor,
                    alpha=0.15, lw=0)
    ax.plot(n[i0:], escala * m[i0:], color=cor, lw=1.2, label=rot)


def fig_chuva(pasta, ts, pil, ai_t, ep, ea, nP):
    """(a) P_f(t) da AI e do direto; (b) P_f com N; (c) pesos"""
    fig, (a, b, c) = plt.subplots(1, 3, figsize=(LARGURA, 2.3))
    a.axvspan(0, ep["chuva"], color=SOLO, lw=0, zorder=0)
    a.text(ep["chuva"] / 2., 0.035, "chuva", transform=a.get_xaxis_transform(), ha="center", va="bottom",
           fontsize=6.5, color=TINTA2)
    t_ai = [t for t in ts if ai_t[t]["pf"] > 0]
    y = np.array([100. * ai_t[t]["pf"] for t in t_ai])
    e = np.array([1.96 * ai_t[t]["cov"] for t in t_ai]) * y
    lo = np.minimum(e, y - 1.2e-3)  # a ponta de baixo cortada no fundo do eixo quando P_f (1 - 1,96 CoV) < 0
    a.errorbar(t_ai, y, yerr=[lo, e], fmt="o", color=COR_AI, ms=3.5, capsize=2, lw=1., label="AI")
    t_d = [t for t in ts if pil[t]["nf"] > 0]
    yd = np.array([100. * pil[t]["pf"] for t in t_d])
    a.errorbar(np.array(t_d) - 1.2, yd, yerr=[yd - 100. * np.array([pil[t]["lo"] for t in t_d]),
                                              100. * np.array([pil[t]["hi"] for t in t_d]) - yd],
               fmt="s", color=COR_DIRETO, ms=3.2, capsize=2, lw=1., label="direto")
    zero_hi = 100. * intervalo(0, nP)[1]
    a.axhline(zero_hi, color=COR_DIRETO, ls=":", lw=0.9)
    a.text(1., zero_hi / 1.3, "direto sem\nfalhas: até %s%%" % v(zero_hi, 1), fontsize=6.2, color=TINTA2,
           ha="left", va="top")
    a.set_yscale("log")
    a.set_ylim(1e-3, 60)
    a.set_xlim(0, ts[-1] + 3)
    a.set_xticks(ts)
    a.set_xlabel("$t$ (h)")
    a.set_ylabel("$P(FS < 1)$ (%)")
    eixo_log_virgula(a, "y")
    a.legend(loc="upper right", fontsize=6.3, handletextpad=0.2, borderaxespad=0.3)
    letra(a, "a")
    i0 = 19
    banda_corrida(b, ep["q"], COR_DIRETO, "direto (95%)", i0)
    banda_corrida(b, ea["q"], COR_AI, "AI (95%)", i0)
    b.set_xlabel("realizações $N$")
    b.set_ylabel("$P_f$ (%)")
    b.set_xlim(0, max(len(ep["q"]), len(ea["q"])))
    b.set_ylim(0, 20)
    b.legend(loc="upper right", fontsize=6.3)
    letra(b, "b")
    f = ea["f"]
    c.scatter(ea["G"][~f], ea["w"][~f], s=5, color=MUDO, alpha=0.6, lw=0)
    c.scatter(ea["G"][f], ea["w"][f], s=6, color=COR_AI, lw=0, alpha=0.9)
    c.axvline(1., color=TINTA2, lw=0.8)
    c.axhline(1., color=MUDO, lw=0.6, ls=":")
    c.set_yscale("log")
    c.set_xlabel("$FS$ mínimo")
    c.set_ylabel("peso $w$")
    c.legend(handles=[Line2D([], [], marker="o", ls="", color=COR_AI, ms=3.5, label="falha ($FS < 1$)"),
                      Line2D([], [], marker="o", ls="", color=MUDO, ms=3.5, label="$FS \\geq 1$")],
             loc="lower right", fontsize=6.3, handletextpad=0.1)
    letra(c, "c")
    fig.tight_layout()
    salva(fig, pasta, "fig_e0_chuva")


def fig_amostra(pasta, ep, ea, nivel, f24, f8v, ln200, ln8, p_ai, ratio_lin):
    """(a) histogramas do FS minimo; (b) distribuicao de FS(24 h): as 200 e as 8; (c) CoV(P_f) com N"""
    fig, (a, b, c) = plt.subplots(1, 3, figsize=(LARGURA, 2.3))
    bins = np.arange(0.80, 1.5001, 0.02)
    h1 = a.hist(ep["G"], bins=bins, color=COR_DIRETO, alpha=0.75, label="piloto (direto)")[0]
    h2 = a.hist(ea["G"], bins=bins, histtype="step", color=COR_AI, lw=1.3, label="AI (sem peso)")[0]
    a.set_ylim(0, 1.5 * max(h1.max(), h2.max()))
    a.axvline(1., color=TINTA2, lw=0.8)
    a.axvline(nivel, color=TINTA, lw=0.9, ls=":")
    a.text(nivel - 0.012, 0.97, "nível do\nelite", transform=a.get_xaxis_transform(), fontsize=6.2, color=TINTA2,
           va="top", ha="right")
    a.set_xlabel("$FS$ mínimo no tempo")
    a.set_ylabel("realizações")
    a.set_xlim(0.8, 1.5)
    a.legend(loc="upper right", fontsize=6.3)
    letra(a, "a")
    s = np.sort(f24)
    b.step(s, np.arange(1, len(s) + 1) / len(s), where="post", color=COR_DIRETO, lw=1.2,
           label="as %d" % len(s))
    s8 = np.sort(f8v)
    b.plot(s8, (np.arange(1, len(s8) + 1) - 0.5) / len(s8), "x", color=COR_AI, ms=4, label="as %d antigas" % len(s8))
    x = np.linspace(0.85, 1.45, 300)
    for (mu, sg), cor, ls in ((ln200, COR_DIRETO, "--"), (ln8, COR_AI, "--")):
        b.plot(x, [norm_cdf((math.log(xx) - mu) / sg) for xx in x], color=cor, ls=ls, lw=0.9)
    b.axvline(1., color=TINTA2, lw=0.8)
    b.set_yscale("log")
    b.set_ylim(1e-3, 1.)
    b.set_xlim(0.85, 1.45)
    b.set_xlabel("$FS$ em 24 h")
    b.set_ylabel("$P(FS \\leq x)$")
    eixo_log_virgula(b, "y")
    hs, ls_ = b.get_legend_handles_labels()
    b.legend(hs + [Line2D([], [], color=TINTA2, ls="--", lw=0.9)], ls_ + ["log-normal"], loc="lower right",
             fontsize=6.3)
    letra(b, "b")
    for q, cor, rot in ((ep["q"], COR_DIRETO, "direto"), (ea["q"], COR_AI, "AI")):
        n, m, se = corrida(q)
        ok = (m > 0) & (n >= 20)
        c.loglog(n[ok], 100. * se[ok] / m[ok], color=cor, lw=1.2, label=rot)
    nn = np.logspace(math.log10(20), math.log10(len(ea["q"])), 50)
    c.loglog(nn, 100. * np.sqrt((1 - p_ai) / (p_ai * nn)), color=COR_DIRETO, ls=":", lw=0.9,
             label="direto, teoria")
    c.loglog(nn, 100. * np.sqrt(ratio_lin * (1 - p_ai) / (p_ai * nn)), color=TINTA2, ls="--", lw=0.9,
             label="AI linear ótima")
    c.set_xlabel("realizações $N$")
    c.set_ylabel("$CoV(P_f)$ (%)")
    c.set_yticks([5, 10, 20, 50, 100])
    eixo_log_virgula(c, "x")
    eixo_log_virgula(c, "y")
    c.legend(loc="lower left", fontsize=6.0)
    letra(c, "c")
    fig.tight_layout()
    salva(fig, pasta, "fig_e0_amostra")


def fig_cho(pasta, niveis, bishop, casos):
    """(a) FS deterministico por nivel; (b, c) P_f com N do direto (e da AI) com r = -0,5 e -0,7"""
    fig, axs = plt.subplots(1, 3, figsize=(LARGURA, 2.3))
    a = axs[0]
    for (rot, vals), cor in zip(niveis, (COR_TERCEIRA, COR_DIRETO)):
        a.plot(range(len(vals)), vals, "o-", color=cor, ms=3.5, label=rot)
    a.axhline(bishop, color=TINTA, ls="--", lw=0.8)
    a.text(0.0, bishop + 0.003, "Cho (Bishop)", fontsize=6.5, color=TINTA2, va="bottom")
    a.set_xticks(range(max(len(vv) for _, vv in niveis)))
    a.set_xlabel("nível de refinamento")
    a.set_ylabel("$FS$ determinístico")
    a.legend(loc="upper right", fontsize=6.3)
    letra(a, "a")
    for ax, (le, d) in zip(axs[1:], zip("bc", casos)):
        i0 = 49
        banda_corrida(ax, d["q_dir"], COR_DIRETO, "direto (95%)", i0)
        for rot, q, cor in d["ai"]:
            banda_corrida(ax, q, cor, rot, min(i0, len(q) - 1))
        ax.axhline(d["cho"], color=TINTA, ls="--", lw=0.8)
        ax.text(0.02, d["cho"], "Cho", transform=ax.get_yaxis_transform(), fontsize=6.5, color=TINTA2, ha="left",
                va="bottom")
        ax.set_xlim(0, max([len(d["q_dir"])] + [len(q) for _, q, _ in d["ai"]]))
        ax.set_ylim(0, d["ymax"])
        ax.set_xlabel("realizações $N$")
        ax.set_ylabel("$P_f$ (%)")
        ax.set_title(d["titulo"], fontsize=7.5, color=TINTA)
        ax.legend(loc="upper right", fontsize=6.3)
        letra(ax, le)
    fig.tight_layout()
    salva(fig, pasta, "fig_e0_cho")


def fig_form(pasta, forms, beta_cho, beta_dir):
    """(a) beta por iteracao; (b, c) alfa^2 por termo na ultima iteracao"""
    fig = plt.figure(figsize=(LARGURA, 2.3))
    gs = fig.add_gridspec(1, 3, width_ratios=[0.8, 1, 1.35], wspace=0.32, left=0.075, right=0.99, bottom=0.2,
                          top=0.88)
    a = fig.add_subplot(gs[0, 0])
    for d, cor in zip(forms, (COR_DIRETO, COR_AI)):
        its = [h["it"] for h in d["hist"]]
        a.plot(its, [h["beta"] for h in d["hist"]], "o-", color=cor, ms=3.5, label="$M = %d$" % d["M"])
        if not d["pronto"]:  # o proximo ponto do HL-RF, ainda sem avaliar
            a.plot([its[-1], its[-1] + 1], [d["hist"][-1]["beta"], d["beta_prox"]], ":", color=cor, lw=1.)
            a.plot([its[-1] + 1], [d["beta_prox"]], "o", color=cor, mfc="white", ms=3.5)
    a.axhline(beta_cho, color=TINTA, ls="--", lw=0.8)
    a.text(0.02, beta_cho - 0.04, "Cho", transform=a.get_yaxis_transform(), fontsize=6.5, color=TINTA2, va="top")
    a.set_xlabel("iteração")
    a.set_ylabel("$\\beta = |u|$")
    a.set_ylim(0, 2.6)
    a.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
    a.legend(loc="lower right", fontsize=6.3)
    letra(a, "a")
    for k, (d, le) in enumerate(zip(forms, "bc")):
        ax = fig.add_subplot(gs[0, 1 + k])
        M, a2 = d["M"], d["alfa2"]
        xc = np.arange(M)
        xp = np.arange(M) + M + 1
        ax.bar(xc, a2[:M], color=COR_DIRETO, width=0.8, label="$c'$")
        ax.bar(xp, a2[M:], color=COR_AI, width=0.8, label="$\\varphi'$")
        top = sorted(range(2 * M), key=lambda i: -a2[i])[:3]
        for i in top:
            xi = i if i < M else i + 1
            ax.text(xi, a2[i] + 0.015, ("c" if i < M else "φ") + "%d" % (i % M), ha="center", va="bottom",
                    fontsize=6.2, color=TINTA)
        ax.set_xticks([(M - 1) / 2., M + 1 + (M - 1) / 2.])
        ax.set_xticklabels(["$c'$", "$\\varphi'$"])
        ax.set_xlim(-1, 2 * M + 1)
        ax.set_ylim(0, 1.)
        ax.set_xlabel("termos 0 a %d de cada campo" % (M - 1))
        if k == 0:
            ax.set_ylabel("$\\alpha^2$")
        ax.set_title("$M = %d$, iteração %d" % (M, d["hist"][-1]["it"]), fontsize=7.5, color=TINTA)
        ax.grid(axis="x", visible=False)
        letra(ax, le)
    salva(fig, pasta, "fig_e0_form")


# ---------------------------------------------------------------- A) chuva com AI no fermi
def chuva(C, B):
    P = os.path.join(C, "chuva_ai", "ref48_piloto")
    A = os.path.join(C, "chuva_ai", "ref48_ai")
    fP, fA = rel(P, C) + "/mc_*.csv", rel(A, C) + "/mc_*.csv"
    fpar = rel(P, C) + "/mc_parametros.txt"
    par = le_chaves(os.path.join(P, "mc_parametros.txt"))
    parA = le_chaves(os.path.join(A, "mc_parametros.txt"))
    txt("ChSolo", par["SOLO"], fpar)
    num("ChI", v(float(par["I"]), 2), fpar)
    chuva_h = float(par["CHUVA"])
    num("ChChuva", "%g" % chuva_h, fpar)
    num("ChHw", "%g" % float(par["HW"]), fpar)
    ts = [float(x) for x in par["TEMPOS"].split(",")]
    num("ChTempoPasso", "%g" % (ts[1] - ts[0]), fpar)
    num("ChTempoFim", "%g" % ts[-1], fpar)
    num("ChNTempos", "%d" % len(ts), fpar)
    for k, mac in (("COV_C", "ChCovC"), ("COV_PHI", "ChCovPhi"), ("COV_K", "ChCovK")):
        num(mac, v(float(par[k]), 1), fpar)
    num("ChLx", "%g" % float(par["LX"]), fpar)
    num("ChLy", "%g" % float(par["LY"]), fpar)
    num("ChKLM", par["KL_M"], fpar)
    num("ChKLvar", v(100. * float(par["KL_VARIANCIA"]), 1), fpar)
    txt("ChMalhaT", par["MALHA_TRANSIENTE"], fpar)
    txt("ChMalhaF", par["MALHA_FATOR"], fpar)
    num("ChNref", par["NREF"], fpar)
    num("ChUt", par["UT"], fpar)
    num("ChSementeP", par["SEMENTE"], fpar)
    fparA = rel(A, C) + "/mc_parametros.txt"
    num("ChSementeA", parA["SEMENTE"], fparA)
    num("ChBandaT", "%g" % float(parA["BANDA_T"]), fparA)
    elite_frac = float(parA["ELITE"])
    num("ChEliteFrac", "%g" % (100. * elite_frac), fparA)
    casos = le_casos(os.path.join(SRF, "cluster", "casos_chuva_ai.txt"))
    num("ChNAlvoP", "%d" % casos["ref48_piloto"][0], "cluster/casos_chuva_ai.txt")
    num("ChNAlvoA", "%d" % casos["ref48_ai"][0], "cluster/casos_chuva_ai.txt")
    # ---- o deterministico (a realizacao -1 do piloto)
    rsP, ruP = le_rodada(P, varios=True)
    rsA, ruA = le_rodada(A, varios=True)
    det = rsP.pop(-1)
    detA = rsA.pop(-1, None)
    assert detA is None or detA["fs_t"] == det["fs_t"], "o deterministico da AI difere do do piloto"
    assert sorted(det["fs_t"]) == ts
    tdet = min(det["fs_t"], key=det["fs_t"].get)
    num("ChDetMin", v(det["fs"]), fP.replace("*", "0_0") + " (realizacao -1)")
    num("ChDetT", "%g" % tdet, fP.replace("*", "0_0") + " (realizacao -1)")
    # so as realizacoes com todos os instantes (o evento() do amostragem_importancia.py)
    inc = {}
    for nome, rs in (("P", rsP), ("A", rsA)):
        inc[nome] = [n for n in rs if sorted(rs[n]["fs_t"]) != ts]
        for n in inc[nome]:
            rs.pop(n)
    num("ChForaP", "%d" % (len(ruP) + len(inc["P"])), fP + " (sem fator ou sem algum instante)")
    num("ChForaA", "%d" % (len(ruA) + len(inc["A"])), fA + " (sem fator ou sem algum instante)")
    # ---- o piloto (Monte Carlo direto)
    ep = estima(rsP)
    assert abs(ep["w_media"] - 1.) < 1e-12, "o piloto tem pesos"
    NP, nfP, pP = ep["N"], ep["nf"], ep["pf"]
    loP, hiP, icP = intervalo(nfP, NP)
    txt("ChPIC", icP, "IC do piloto (%d falhas)" % nfP)
    covP = ep["cov"]  # s_q / (raiz(N) P_f), como o amostragem_importancia.py (~ raiz((1 - P_f) / (N P_f)))
    num("ChNP", "%d" % NP, fP)
    num("ChPNf", "%d" % nfP, fP + " (FS < 1 em algum instante)")
    num("ChPPf", v(100. * pP, 2), fP)
    num("ChPLo", v(100. * loP, 1), fP + " (IC %s)" % icP)
    num("ChPHi", v(100. * hiP, 1), fP + " (IC %s)" % icP)
    num("ChPCov", v(100. * covP, 1), fP + " (s_q / (raiz(N) P_f))")
    tmin = {}
    for n in rsP:
        tm = min(rsP[n]["fs_t"], key=rsP[n]["fs_t"].get)
        tmin[tm] = tmin.get(tm, 0) + 1
    num("ChMinNoDet", "%d" % tmin.get(tdet, 0), fP + " (realizacoes com o minimo no instante do deterministico)")
    # toda falha de outro instante ja falha no instante do deterministico?
    fora = sum(1 for rs in (rsP, rsA) for n in rs
               if any(x < 1. for t, x in rs[n]["fs_t"].items() if t != tdet) and rs[n]["fs_t"][tdet] >= 1.)
    num("ChFalhaFora", "%d" % fora, fP + " e " + fA + " (FS < 1 noutro instante e >= 1 em %g h)" % tdet)
    txt("ChFalhaForaFrase", "nenhuma realização das duas amostras falha" if fora == 0 else
        ("%d realização das duas amostras falha" if fora == 1 else "%d realizações das duas amostras falham") % fora,
        fP + " e " + fA)
    # ---- a AI
    ea = estima(rsA)
    NA = ea["N"]
    num("ChNA", "%d" % NA, fA)
    num("ChANf", "%d" % ea["nf"], fA + " (falhas sorteadas)")
    num("ChAFrac", v(100. * ea["nf"] / NA, 0), fA)
    num("ChAPf", v(100. * ea["pf"], 2), fA + " (soma(w I) / N)")
    num("ChAPfArred", "%d" % round(100. * ea["pf"]), fA + " (arredondada)")
    num("ChACov", v(100. * ea["cov"], 1), fA + " (s_q / (raiz(N) P_f))")
    num("ChALo", v(100. * ea["pf"] * (1 - 1.96 * ea["cov"]), 1), fA + " (P_f (1 - 1,96 CoV))")
    num("ChAHi", v(100. * ea["pf"] * (1 + 1.96 * ea["cov"]), 1), fA + " (P_f (1 + 1,96 CoV))")
    num("ChAPesoMedio", v(ea["w_media"], 3), fA)
    num("ChAPesoSe", v(ea["w_se"], 3), fA + " (desvio da media)")
    zw = abs(ea["w_media"] - 1.) / ea["w_se"]
    txt("ChAPesoCompat", "compatível com 1" if zw < 2. else "a %s desvios de 1" % tex(zw, 1),
        fA + " (|media - 1| / desvio = %.2f)" % zw)
    num("ChAEss", "%d" % round(ea["ess"]), fA + " ((soma w)^2 / soma w^2)")
    num("ChAEssF", "%d" % round(ea["q"].sum() ** 2 / (ea["q"] ** 2).sum()), fA + " (ESS das falhas: (soma q)^2 / soma q^2)")
    num("ChPEssF", "%d" % round(ep["q"].sum() ** 2 / (ep["q"] ** 2).sum()), fP + " (ESS das falhas: no direto, n_f)")
    num("ChPParcela", v(100. * ep["parcela"], 0), fP + " (maior parcela: no direto, 1 / n_f)")
    ma = le_plano(os.path.join(PLANO, "README.md"), r"com alerta acima de (\d+)%", "o alerta da maior parcela")
    num("AlertaParcela", ma[1] if ma else "--", "../README.md (criterios gerais)")
    num("ChAPesoFalha", v(ea["wf_media"], 2), fA + " (peso medio nas falhas)")
    num("ChAPesoFalhaMax", v(ea["wf_max"], 2), fA + " (maior peso nas falhas)")
    num("ChAPesoMax", v(ea["w_max"], 1), fA + " (maior peso)")
    num("ChAParcela", v(100. * ea["parcela"], 0), fA + " (maior w I / soma w I)")
    dif, gan, cov_d = compara(ea, pP, covP)
    num("ChDif", v(dif, 2), fA + " e " + fP + " (desvios)")
    num("ChDifPontos", v(100. * (ea["pf"] - pP), 2), fA + " e " + fP)
    num("ChGanho", v(gan, 1), fA + " ((CoV do direto com N e P_f da AI / CoV da AI)^2)")
    num("ChRazao", v(1. / gan, 2), fA + " (1 / ganho)")
    num("ChCovMesmoN", v(100. * cov_d, 1), "raiz((1 - P_f) / (N P_f)) com N e P_f da AI")
    num("ChNParaCov", mil((1 - ea["pf"]) / (ea["pf"] * ea["cov"] ** 2)), "(1 - P_f) / (P_f CoV^2) com os da AI")
    txt("ChVies", "sem sinal de viés" if abs(dif) < 2 else "acima de dois desvios, a conferir", fA)
    # a conferencia com o resultado.txt do amostragem_importancia.py (feito no fermi)
    res = os.path.join(C, "resultado.txt")
    mr = re.search(r"AI: N = (\d+).*?P_f \(.*?\) = ([\d.]+)%  CoV\(P_f\) = ([\d.]+)%.*?direto: N = (\d+), falhas = (\d+)",
                   texto(res), re.S)
    print("resultado.txt: AI N = %s, P_f = %s%%, CoV = %s%%; direto N = %s, falhas %s (aqui: %d, %.3f%%, %.1f%%; %d, %d)"
          % (mr[1], mr[2], mr[3], mr[4], mr[5], NA, 100 * ea["pf"], 100 * ea["cov"], NP, nfP))
    assert int(mr[1]) == NA and abs(float(mr[2]) - 100 * ea["pf"]) < 2e-3 and int(mr[5]) == nfP
    # ---- por instante: o piloto (contado) e a AI (os mesmos pesos)
    pil, ai_t = {}, {}
    linhas = []
    for t in ts:
        x = np.array([rsP[n]["fs_t"][t] for n in sorted(rsP)])
        k = int((x < 1.).sum())
        lo, hi, _ = intervalo(k, NP)
        pil[t] = {"nf": k, "pf": k / NP, "lo": lo, "hi": hi, "media": x.mean(), "cov": x.std(ddof=1) / x.mean()}
        ai_t[t] = estima(rsA, [t])
        e = ai_t[t]
        nome = extenso(t)
        num("ChT" + nome, "%g" % t, fpar + " (TEMPOS)")
        num("ChPfA" + nome, pct(100. * e["pf"]), fA + " (P_f em %g h)" % t)
        num("ChCovA" + nome, v(100. * e["cov"], 0) if e["pf"] > 0 else "--", fA + " (CoV em %g h)" % t)
        num("ChNfP" + nome, "%d" % k, fP + " (FS < 1 em %g h)" % t)
        txt("ChFalhasP" + nome, falhas(k), fP + " (FS < 1 em %g h)" % t)
        linhas.append(r"%g & %s & %s & %s & %d & %s & %s & %s \\" % (
            t, tex(det["fs_t"][t], 4), tex(pil[t]["media"], 4), tex(100. * pil[t]["cov"], 1), k,
            tex(100. * k / NP, 2) if k else "0", pct(100. * e["pf"]).replace(",", "{,}") if e["pf"] > 0 else "0",
            tex(100. * e["cov"], 0) if e["pf"] > 0 else "--"))
    txt("ChLinhasInstantes", "\n".join(linhas), fP + ", " + fA + " (por instante)")
    zeros = [t for t in ts if ai_t[t]["pf"] == 0]
    txt("ChZerosA", (" e ".join(", ".join("%g" % t for t in zeros).rsplit(", ", 1)) if zeros else "--"),
        fA + " (instantes sem falha sorteada)")
    num("ChZeroHi", v(100. * intervalo(0, NP)[1], 1), "IC exato com 0 falhas em %d" % NP)
    # ---- o deslocamento (mc_desloc_*.txt) e o nivel do elite
    nb, iguais, d = le_desloc(A)
    fd = rel(A, C) + "/mc_desloc_*.txt"
    num("ChDeslocBlocos", "%d" % nb, fd)
    txt("ChDeslocIguais", "o mesmo" if iguais else "\\textbf{deslocamentos diferentes}", fd)
    s = [float(x) for x in d["s"].split(",")]
    mu2 = float(d["mu2"])
    assert abs(sum(x * x for x in s) - mu2) < 1e-6
    num("ChElite", d["elite"], fd)
    num("ChNivel", v(float(d["nivel"]), 4), fd + " (FS do conjunto elite)")
    num("ChPilotoFalhas", d["falhas"], fd)
    for mac, x in zip(("ChSC", "ChSPhi", "ChSK"), s):
        num(mac, v(x, 2), fd)
    num("ChMuDois", v(mu2, 2), fd)
    num("ChVarPeso", v(math.exp(mu2) - 1., 1), "exp(|mu|^2) - 1")
    # o elite: as falhas ou, com menos que a fracao ELITE das validas, as de menor FS (ImportanceSampling.h)
    m_el = int(round(elite_frac * int(d["validas"])))
    G_ord = sorted((rsP[n]["fs"], n) for n in rsP)
    assert int(d["falhas"]) == nfP and int(d["validas"]) == NP
    assert int(d["elite"]) == (m_el if nfP < m_el else nfP)
    assert abs(G_ord[m_el - 1][0] - float(d["nivel"])) < 5e-5, (G_ord[m_el - 1], d["nivel"])
    pm = ea["proj"].mean(axis=0)
    se_p = 1. / math.sqrt(NA)
    num("ChProjMedia", "; ".join(v(x, 2) for x in pm), fA + " (media de proj_* sorteadas)",
        latex=r"\ensuremath{" + ";\\ ".join(tex(x, 2) for x in pm) + "}")
    num("ChProjSe", v(se_p, 2), "1 / raiz(N)")
    dp = [abs(a - b) / se_p for a, b in zip(pm, s)]
    imax = int(np.argmax(dp))
    num("ChProjMaxDif", v(dp[imax], 2), fA + " e " + fd + " (desvios, o maior dos tres campos)")
    txt("ChProjMaxCampo", ("$c'$", r"$\varphi'$", "$k$")[imax], fA + " (o campo da maior diferenca)")
    num("ChProjOutrasDif", " e ".join(v(x, 2) for i, x in enumerate(dp) if i != imax), fA + " (os outros dois campos)")
    # ---- a teoria: o estado limite linear com a P_f da AI
    b, m, r = ganho_ai(ea["pf"])
    num("ChBeta", v(b, 2), "formula, P_f da AI")
    num("ChMdois", v(m * m, 2), "formula, P_f da AI (m = phi(beta) / P_f)")
    num("ChGanhoLinear", v(1. / r, 1), "formula, P_f da AI (exp(m^2) Phi(-beta - m))")
    gmu = ganho_desloc(ea["pf"], math.sqrt(mu2))
    num("ChGanhoLinearMu", v(gmu, 1), "formula, P_f da AI e |mu| do deslocamento usado")
    num("ChFracGanho", v(100. * gan * r, 0), "ganho medido / ganho linear")
    # ---- custo no fermi (logs do tarefas.sh, um bloco por realizacao)
    L = os.path.join(C, "chuva_ai", "logs")
    tar = le_tarefas(L)
    fl = rel(L, C) + "/p*.txt"
    num("ChProcessos", "%d" % len(glob.glob(os.path.join(L, "*.txt"))), fl)
    num("ChNos", "%d" % len(le_hosts(L)), fl + " (nos)")
    ini = min(x["inicio"] for x in tar.values())
    fim = max(x["fim"] for x in tar.values() if "fim" in x)
    txt("ChInicio", hora(ini), fl)
    txt("ChFim", hora(fim), fl)
    txt("ChRelogio", duracao(fim - ini), fl)
    sP = [x["s"] for k, x in tar.items() if k[0] == "ref48_piloto" and "s" in x]
    sA = [x["s"] for k, x in tar.items() if k[0] == "ref48_ai" and "s" in x]
    num("ChCustoP", "%d" % round(np.mean(sP)), fl + " (s por realizacao, piloto)")
    num("ChCustoA", "%d" % round(np.mean(sA)), fl + " (s por realizacao, AI)")
    num("ChCpu", "%d" % round((sum(sP) + sum(sA)) / 3600.), fl + " (soma dos tempos, h)")
    # ---- por que ~6% e nao ~1%: o piloto antigo de 8 realizacoes (relatorio_chuva.tex)
    arq8 = os.path.join(B, "mc", "gim_I", "piloto_srm_I31.05_0_7.csv")
    with open(arq8, newline="") as f:
        velho = {int(r["real"]): float(r["FS_min"]) for r in csv.DictReader(f) if float(r["t_h"]) == 24.}
    f8 = rel(arq8, os.path.dirname(os.path.dirname(B)))
    n8 = len(velho)
    f8v = np.array([velho[n] for n in sorted(velho)])
    num("ChVelhoN", "%d" % n8, f8)
    num("ChVelhoPrimeira", "%d" % min(velho), f8)
    num("ChVelhoUltima", "%d" % max(velho), f8)
    num("ChVelhoMedia", v(f8v.mean()), f8)
    num("ChVelhoCov", v(100. * f8v.std(ddof=1) / f8v.mean(), 1), f8)
    num("ChVelhoMin", v(f8v.min()), f8)
    pl8, mu8, sg8 = lognormal(f8v)
    num("ChVelhoLogn", v(100. * pl8, 2), f8 + " (ajuste log-normal)")
    mlog = re.search(r"\n\s+24\s+[\d.]+\s+[\d.]+\s+[\d.]+\s+[\d.]+%\s+([\d.]+)%", texto(os.path.join(
        os.path.dirname(arq8), "piloto_srm.log")))
    print("piloto antigo: P_f lognormal no log = %s%%, aqui %.2f%%" % (mlog[1], 100 * pl8))
    assert abs(float(mlog[1]) - 100. * pl8) < 0.006
    nf8 = int((f8v < 1.).sum())
    num("ChVelhoNf", "%d" % nf8, f8 + " (FS < 1)")
    lo8, hi8, ic8 = intervalo(nf8, n8)
    num("ChVelhoHi", v(100. * hi8, 0), "IC %s com %d falhas em %d" % (ic8, nf8, n8))
    assert int(par["SEMENTE"]) == int(le_chaves(os.path.join(os.path.dirname(arq8),
                                                             "piloto_srm_I31.05_parametros.txt"))["SEMENTE"])
    f24 = np.array([rsP[n]["fs_t"][24.] for n in sorted(rsP)])
    o8 = np.array([rsP[n]["fs_t"][24.] for n in sorted(velho)])
    num("ChVelhoDifMax", v(np.abs(o8 - f8v).max(), 4), fP + " e " + f8 + " (maior |diferenca| de FS(24 h))")
    num("ChOitoMedia", v(o8.mean()), fP + " (realizacoes 0..%d, FS(24 h))" % max(velho))
    num("ChOitoCov", v(100. * o8.std(ddof=1) / o8.mean(), 1), fP + " (realizacoes 0..%d)" % max(velho))
    num("ChDuzMedia", v(f24.mean()), fP + " (FS(24 h))")
    num("ChDuzCov", v(100. * f24.std(ddof=1) / f24.mean(), 1), fP + " (FS(24 h))")
    pl200, mu200, sg200 = lognormal(f24)
    num("ChDuzLogn", v(100. * pl200, 2), fP + " (ajuste log-normal de FS(24 h))")
    num("ChCruzA", v(100. * norm_cdf(-mu200 / sg8), 1), "log-normal: media de ln FS das %d, desvio das %d" % (NP, n8))
    num("ChCruzB", v(100. * norm_cdf(-mu8 / sg200), 1), "log-normal: media de ln FS das %d, desvio das %d" % (n8, NP))
    rd = os.path.join(SRF, "cluster", "README.md")
    m1 = le_fonte(rd, r"tem P_f ≈ (\d+)% em 24 h, pelos pilotos", "P_f esperada")
    num("ChPfEsperada", m1[1], "cluster/README.md (pilotos do relatorio_chuva.tex)")
    m2 = le_fonte(rd, r"CoV\(P_f\) de ~(\d+)% \(estado limite linear\) a ~(\d+)% com (\d+) realizações", "precisao")
    num("ChCovEspLin", m2[1], "cluster/README.md")
    num("ChCovEspMax", m2[2], "cluster/README.md")
    # ---- quantas realizacoes diretas a chuva pede para o CoV alvo, com a P_f da AI, e o custo no fermi
    n4 = (1 - ea["pf"]) / (ea["pf"] * COV_ALVO ** 2)
    num("ChCovAlvo", "%g" % (100. * COV_ALVO), "COV_ALVO do script")
    n4r = 100 * round(n4 / 100.)  # em centenas (arredondamento normal); o custo usa este mesmo N
    print("chuva: N para CoV %g%% = %.1f, em centenas %d" % (100. * COV_ALVO, n4, n4r))
    num("ChNAlvoCov", "%d" % n4r, "(1 - P_f) / (P_f CoV^2) com a P_f da AI, em centenas")
    num("ChCpuAlvo", "%d" % round(n4r * np.mean(sP) / 3600.), "N (em centenas) x custo do piloto no fermi")
    num("ChRelogioAlvo", v(n4r * np.mean(sP) / 3600. / len(glob.glob(os.path.join(L, "*.txt"))), 1),
        "N (em centenas) x custo do piloto / processos")
    return {"ts": ts, "pil": pil, "ai_t": ai_t, "ep": dict(ep, chuva=chuva_h), "ea": ea, "nP": NP,
            "nivel": float(d["nivel"]), "f24": f24, "f8v": f8v, "ln200": (mu200, sg200), "ln8": (mu8, sg8),
            "ratio_lin": r, "gan": gan, "frac": gan * r, "dif": dif, "pP": pP, "covP": covP, "pA": ea["pf"],
            "sd_dif": math.sqrt((ea["pf"] * ea["cov"]) ** 2 + (pP * covP) ** 2),
            "covA": ea["cov"], "NA": NA, "ganho_lin": 1. / r, "contP": contagens(ep["G"]), "contA": contagens(ea["G"])}


# ---------------------------------------------------------------- B) Cho (2010) no PC
def le_form(pasta):
    """o FORM do form_kl.py: form.json no fim, form_parcial.json enquanto roda. Por iteracao: g, beta, alfa^2 do
    gradiente; o proximo ponto do HL-RF; a iteracao em curso (avaliacoes feitas)"""
    fj, fp = os.path.join(pasta, "form.json"), os.path.join(pasta, "form_parcial.json")
    pronto = os.path.exists(fj)
    if not pronto and not os.path.exists(fp):
        return None
    d = json.load(open(fj if pronto else fp))
    M = d["parametros"]["M"]
    hist = d["historico"]
    for h in hist:
        n2 = sum(x * x for x in h["grad"])
        h["alfa2"] = np.array([x * x / n2 for x in h["grad"]])
        h["passo"] = None
    ult = hist[-1]
    g, grad, u = ult["g"], ult["grad"], ult["u"]
    n2 = sum(x * x for x in grad)
    gu = sum(a * b for a, b in zip(grad, u))
    novo = [(gu - g) / n2 * x for x in grad]
    passo = math.sqrt(sum((a - b) ** 2 for a, b in zip(novo, u)))
    src = texto(os.path.join(pasta, "form_kl.py"))
    mt = re.search(r"abs\(g\) < ([\d.e-]+) and passo < ([\d.e-]+)", src)
    tol_g, tol_p = float(mt[1]), float(mt[2])
    convergiu = len(hist) > 1 and abs(g) < tol_g and passo < tol_p
    em_curso = len(hist)  # a iteracao que roda agora (sem o form.json)
    feitas = len([a for a in os.listdir(pasta) if re.fullmatch(r"it%d_[a-z0-9]+_0_0\.csv" % em_curso, a)])
    # o passo em u entre as duas ultimas iteracoes (o que o form_kl.py compara com a tolerancia)
    passos = [math.sqrt(sum((a - b) ** 2 for a, b in zip(h1["u"], h0["u"]))) for h0, h1 in zip(hist, hist[1:])]
    out = {"M": M, "pronto": pronto, "hist": hist, "passos": passos, "p": d["parametros"], "beta_prox": math.sqrt(sum(x * x for x in novo)),
           "convergiu": convergiu, "tol_g": tol_g, "tol_p": tol_p, "em_curso": em_curso, "feitas": feitas,
           "aval": d["avaliacoes"], "alfa2": ult["alfa2"], "grad": np.array(grad), "u": np.array(u), "g": g}
    if pronto:
        out.update(beta=d["beta"], pf=d["P_f"], ustar=np.array(d["u"]), alfa2=np.array(d["alfa2"]))
    else:
        out.update(beta=ult["beta"], pf=norm_cdf(-ult["beta"]), ustar=None)
    par = le_chaves(os.path.join(pasta, "it0_centro_parametros.txt"))
    out["var"] = float(par["KL_VARIANCIA"])
    c0 = next(csv.DictReader(open(os.path.join(pasta, "it0_centro_0_0.csv"), newline="")))
    out["centro"] = (float(c0["FS_min"]), float(c0["c_med"]), float(c0["phi_med"]))
    return out


def eta_fila(rs_casos, tar, ordem, t_real, t_ref, nproc):
    """o inicio e o fim de cada bloco que falta na fila do tarefas.sh, simulado: os blocos em curso acabam com o que
    falta das realizacoes deles, e os outros comecam na ordem dada quando um processo fica livre. Custo de t_real s por
    realizacao (o do direto ate agora: para a AI, que sorteia mais falhas, e um limite de baixo)"""
    livres, ev = [], {}
    for (caso, a, z), x in tar.items():
        if "inicio" in x and "fim" not in x:
            falta = sum(1 for n in range(a, z + 1) if n not in rs_casos.get(caso, {}))
            t1 = t_ref + dt.timedelta(seconds=falta * t_real)
            livres.append(t1)
            ev[(caso, a, z)] = (x["inicio"], t1)
    livres += [t_ref] * max(nproc - len(livres), 0)
    heapq.heapify(livres)
    for k in ordem:
        if k in tar:
            continue  # feito ou em curso
        t0 = heapq.heappop(livres)
        t1 = t0 + dt.timedelta(seconds=(k[2] - k[1] + 1) * t_real)
        ev[k] = (t0, t1)
        heapq.heappush(livres, t1)
    return ev


def cho(H, rch):
    K = os.path.join(H, "campanha")
    casos = le_casos(os.path.join(H, "casos.txt"))
    rod = texto(os.path.join(H, "rodar.sh"))
    Bl = int(re.search(r"BLOCO=(\d+)", rod)[1])
    num("ChoBloco", "%d" % Bl, "rodar.sh")
    # ---- o problema (mc_parametros.txt do direto) e os numeros do artigo (Projects/Cho2010/README.md)
    D5 = os.path.join(K, "r05_direto")
    par = le_chaves(os.path.join(D5, "mc_parametros.txt"))
    fpar = "campanha/r05_direto/mc_parametros.txt"
    num("ChoC", "%g" % float(par["C"]), fpar)
    num("ChoPhi", "%g" % float(par["PHI"]), fpar)
    num("ChoGamma", "%g" % float(par["GAMMA"]), fpar)
    num("ChoE", "%g" % (float(par["E"]) / 1000.), fpar + " (MPa)")
    num("ChoNu", v(float(par["NU"]), 1), fpar)
    num("ChoCovC", v(float(par["COV_C"]), 1), fpar)
    num("ChoCovPhi", v(float(par["COV_PHI"]), 1), fpar)
    num("ChoLh", "%g" % float(par["LH"]), fpar)
    num("ChoLv", "%g" % float(par["LV"]), fpar)
    num("ChoKLM", par["KL_M"], fpar)
    num("ChoKLvar", v(100. * float(par["KL_VARIANCIA"]), 1), fpar)
    num("ChoSemente", par["SEMENTE"], fpar)
    txt("ChoMalha", par["MALHA_FATOR"], fpar)
    num("ChoNref", par["NREF"], fpar)
    rd = os.path.join(CHO_SRC, "README.md")
    fr = "Projects/Cho2010/README.md (Cho 2010)"
    bishop = float(le_fonte(rd, r"FS determinístico = ([\d,]+)", "FS do Cho")[1].replace(",", "."))
    pf_c5 = float(le_fonte(rd, r"P_f = ([\d,]+) com r = −0,5 \(Tabela 5\)", "P_f do Cho")[1].replace(",", "."))
    pf_c7 = float(le_fonte(rd, r"P_f ≈ ([\d,]+) com r = −0,7 \(Fig\. 9\)", "P_f do Cho")[1].replace(",", "."))
    n_c = int(le_fonte(rd, r"Monte Carlo com (\d[\d ]*\d)\s+realizações", "N do Cho")[1].replace(" ", ""))
    num("ChoBishop", v(bishop, 3), fr)
    num("ChoPfCinco", v(100. * pf_c5, 2), fr + ": Tabela 5")
    num("ChoPfSete", v(100. * pf_c7, 1), fr + ": Fig. 9 (aproximado)")
    num("ChoNmc", mil(n_c), fr)
    beta_cho = norm_inv(pf_c5)
    num("ChoBetaCinco", v(beta_cho, 2), "-Phi^-1(P_f do Cho)")
    rc5, rc7 = float(casos["r05_direto"][1]["rcphi"]), float(casos["r07_direto"][1]["rcphi"])
    num("ChoRCinco", v(rc5, 1), "casos.txt")
    num("ChoRSete", v(rc7, 1), "casos.txt")
    num("ChoRaiz", v(math.sqrt(1. - rc5 * rc5), 3), "raiz(1 - r^2)")
    elite = float(le_fonte(os.path.join(CHO_SRC, "cho2010.cpp"), r"elite \(([\d.]+)\)", "elite")[1])
    num("ChoEliteFrac", "%g" % (100. * elite), "cho2010.cpp (padrao de elite)")
    num("ChoPilotoN", casos["r05_ai"][1]["piloto_n"], "casos.txt (piloto_n)")
    num("ChoEliteN", "%d" % round(elite * int(casos["r05_ai"][1]["piloto_n"])), "elite x piloto_n")
    # ---- o deterministico (det/) e a malha
    DET = os.path.join(H, "det")
    det = {}
    for nome, arq in (("TriUm", "srm_fator_tri1_nref_4_0_-1.csv"), ("TriDoisTres", "srm_fator_tri2_nref_3_0_-1.csv"),
                      ("TriDois", "srm_fator_tri2_nref_4_0_-1.csv"), ("Gim", "gim_0_-1.csv")):
        r = next(csv.DictReader(open(os.path.join(DET, arq), newline="")))
        det[nome] = r
        fd = "det/" + arq
        num("ChoDet" + nome, v(float(r["FS_min"])), fd)
        num("ChoEq" + nome, r["eq_fator"], fd + " (equacoes no ultimo nivel)")
        num("ChoNiv" + nome, "%d" % (len(r["niveis"].split(";")) - 1), fd + " (niveis)")
        if nome != "Gim":
            num("ChoDif" + nome, ("+" if float(r["FS_min"]) >= bishop else "") +
                v(100. * (float(r["FS_min"]) / bishop - 1.), 2), fd + " e " + fr)
    nv = {k: [float(x) for x in det[k]["niveis"].split(";")] for k in ("TriUm", "TriDois")}
    num("ChoUltNivTriDois", v(100. * (nv["TriDois"][-1] / nv["TriDois"][-2] - 1.), 2), "det/srm_fator_tri2_nref_4 (niveis)")
    num("ChoUltNivTriUm", v(100. * (nv["TriUm"][-1] / nv["TriUm"][-2] - 1.), 2), "det/srm_fator_tri1_nref_4 (niveis)")
    assert nv["TriDois"][:4] == [float(x) for x in det["TriDoisTres"]["niveis"].split(";")]
    # ---- o Monte Carlo direto (parcial) e a AI
    tar = le_tarefas(os.path.join(K, "logs"))
    nproc = len(glob.glob(os.path.join(K, "logs", "*.txt")))
    num("ChoProcessos", "%d" % nproc, "campanha/logs/pp*.txt")
    ini = min(x["inicio"] for x in tar.values())
    txt("ChoInicio", hora(ini), "campanha/logs/pp*.txt")
    rs = {}
    out = {"casos": [], "bishop": bishop, "beta_cho": beta_cho}
    for caso in casos:
        pasta = os.path.join(K, caso)
        rs[caso] = le_rodada(pasta)[0] if os.path.isdir(pasta) else {}
        rs[caso].pop(-1, None)
    t_csv = [os.path.getmtime(a) for a in glob.glob(os.path.join(K, "*", "mc_*_*.csv"))]
    t_ref = dt.datetime.fromtimestamp(max(t_csv))
    txt("ChoHoraDados", hora(t_ref), "campanha/*/mc_*.csv (o mais novo)")
    # custo por realizacao (logs do programa)
    tempos = {}
    for caso in ("r05_direto", "r07_direto"):
        tempos.update({(caso, n): x for n, x in le_tempos_cho(os.path.join(K, caso)).items() if n >= 0})
    tf = [s for (fs, s) in tempos.values() if fs < 1.]
    to = [s for (fs, s) in tempos.values() if fs >= 1.]
    t_real = float(np.mean([s for _, s in tempos.values()]))
    num("ChoCustoFalha", "%d" % round(np.mean(tf)), "campanha/r0*_direto/mc_*.log (s, falhas)")
    num("ChoCustoOutras", "%d" % round(np.mean(to)), "campanha/r0*_direto/mc_*.log (s, as outras)")
    num("ChoCustoMedio", "%d" % round(t_real), "campanha/r0*_direto/mc_*.log (s)")
    num("ChoCustoRazao", v(np.mean(tf) / np.mean(to), 1), "campanha/r0*_direto/mc_*.log")
    feitas = sum(len(x) for x in rs.values())
    total = sum(n for n, _ in casos.values())
    num("ChoFeitas", "%d" % feitas, "campanha/*/mc_*.csv")
    num("ChoTotal", "%d" % total, "casos.txt")
    fim_campanha = os.path.exists(os.path.join(H, "resultado.txt")) and " fim" in texto(os.path.join(H, "rodar.log"))
    flag("chofim", fim_campanha, "resultado.txt e rodar.log")
    mp = le_plano(os.path.join(PLANO, "README.md"), r"O PC parou em (\d+)/(\d+) às (\d+):(\d+)", "a parada do PC")
    parado = not fim_campanha and (mp is not None or (AGORA - t_ref).total_seconds() > 7200.)
    flag("pcparou", parado, "../README.md (O PC parou ...) ou sem escrita nas rodadas ha 2 h")
    txt("PcParou", "%sh%s" % (mp[3], mp[4]) if mp else hora(t_ref), "../README.md" if mp else "campanha/*/mc_*.csv")
    # a ordem da fila do tarefas.sh: o bloco 0 de todos os casos, depois o 1, ... Cada processo anda para a frente na
    # fila e so volta ao comeco quando chega ao fim; como os pilotos acabam depois que os processos passam pelos blocos
    # da AI, a AI fica para a segunda volta, depois de todos os blocos do direto
    nb_max = max(len(blocos(n, Bl)) for n, _ in casos.values())
    diretos = [c for c, (_, o) in casos.items() if "piloto" not in o]
    ais = [c for c, (_, o) in casos.items() if "piloto" in o]
    ordem = [(c, a, z) for grupo in (diretos, ais) for b in range(nb_max) for c in grupo
             for (a, z) in blocos(casos[c][0], Bl)[b:b + 1]]
    ev = eta_fila(rs, tar, ordem, t_real, t_ref, nproc) if not parado else {}
    fim_dir = [t1 for k, (t0, t1) in ev.items() if k[0] in diretos]
    ini_ai = [t0 for k, (t0, t1) in ev.items() if k[0] in ais and k not in tar]
    txt("ChoFimDireto", hora_meia(max(fim_dir)) if fim_dir else "--", "fila simulada (estimativa)")
    txt("ChoInicioAI", hora_meia(min(ini_ai)) if ini_ai else "--", "fila simulada (estimativa)")
    txt("ChoFimMin", hora_meia(max(t1 for _, t1 in ev.values())) if ev else "--", "fila simulada (limite de baixo)")
    # as realizacoes em curso: as falhas custam mais, entao ha mais falhas entre elas
    em_curso = sum(1 for x in tar.values() if "inicio" in x and "fim" not in x)
    num("ChoEmCurso", "%d" % em_curso, "campanha/logs/pp*.txt (blocos em curso)")
    for caso, nome_c in (("r05_direto", "Cinco"), ("r07_direto", "Sete")):
        num("Cho%sEmCurso" % nome_c, "%d" % sum(1 for k, x in tar.items() if k[0] == caso and "inicio" in x
                                                 and "fim" not in x),
            "campanha/logs/pp*.txt (blocos de %s em curso, uma realizacao cada)" % caso)
    for chave, nome_c, rcho, pfc, ymax in (("r05", "Cinco", rc5, pf_c5, 6.), ("r07", "Sete", rc7, pf_c7, 2.)):
        cd, ca = chave + "_direto", chave + "_ai"
        ND, _ = casos[cd]
        rd_ = rs[cd]
        fD = "campanha/%s/mc_*.csv" % cd
        ed = estima(rd_)
        G, N, nf, p = ed["G"], ed["N"], ed["nf"], ed["pf"]
        lo, hi, ic = intervalo(nf, N)
        txt("Cho%sIC" % nome_c, ic, fD + " (%d falhas)" % nf)
        covp = ed["cov"]  # s_q / (raiz(N) P_f), como o amostragem_importancia.py
        nfeitos, nbl = feitos(os.path.join(K, cd), ND, Bl)
        pronto = nfeitos == nbl
        flag("cho" + nome_c.lower() + "pronto", pronto, "campanha/%s/mc_mapa_*.csv (%d blocos)" % (cd, nbl))
        npil, nbp = feitos(os.path.join(K, cd), ND, Bl, ate=int(casos[ca][1]["piloto_n"]))
        flag("cho" + nome_c.lower() + "piloto", npil == nbp, "campanha/%s/mc_mapa_*.csv (as %s do piloto)"
             % (cd, casos[ca][1]["piloto_n"]))
        num("Cho%sN" % nome_c, "%d" % N, fD)
        num("Cho%sAlvo" % nome_c, "%d" % ND, "casos.txt")
        num("Cho%sNf" % nome_c, "%d" % nf, fD + " (FS < 1)")
        txt("Cho%sFalhas" % nome_c, falhas(nf), fD + " (FS < 1)")
        num("Cho%sPf" % nome_c, v(100. * p, 2), fD)
        num("Cho%sLo" % nome_c, v(100. * lo, 2), fD + " (IC %s)" % ic)
        num("Cho%sHi" % nome_c, v(100. * hi, 2), fD + " (IC %s)" % ic)
        # a malha: o FS do penultimo nivel (coluna niveis) contra o do ultimo, nas mesmas realizacoes
        nvl = np.array([[float(x) for x in rd_[n]["linha"]["niveis"].split(";")] for n in sorted(rd_)])
        dn = nvl[:, -2] - nvl[:, -1]
        num("Cho%sNivDif" % nome_c, v(dn.mean(), 4), fD + " (FS do nivel %d - FS do nivel %d, media)"
            % (nvl.shape[1] - 2, nvl.shape[1] - 1))
        num("Cho%sNivDp" % nome_c, v(dn.std(ddof=1), 4), fD + " (o mesmo, desvio)")
        num("Cho%sNivNf" % nome_c, "%d" % int((nvl[:, -2] < 1.).sum()), fD + " (FS < 1 no penultimo nivel)")
        num("Cho%sNivTrocas" % nome_c, "%d" % int(((nvl[:, -2] < 1.) != (nvl[:, -1] < 1.)).sum()),
            fD + " (realizacoes que trocam de lado entre os dois niveis)")
        num("Cho%sCov" % nome_c, v(100. * covp, 0) if p > 0 else "--", fD + " (s_q / (raiz(N) P_f))")
        num("Cho%sMedia" % nome_c, v(G.mean()), fD)
        num("Cho%sCovFS" % nome_c, v(100. * G.std(ddof=1) / G.mean(), 1), fD)
        pln = lognormal(G)[0]
        num("Cho%sLogn" % nome_c, v(100. * pln, 2), fD + " (ajuste log-normal)")
        txt("Cho%sDentro" % nome_c, "dentro" if lo <= pfc <= hi else "fora", fD + " e " + fr)
        fc = p * np.mean(tf) / (p * np.mean(tf) + (1 - p) * np.mean(to)) if p > 0 else 0.
        num("Cho%sFracCurso" % nome_c, v(100. * fc, 0) if p > 0 else "--",
            "P_f t_falha / (P_f t_falha + (1 - P_f) t_outras)")
        nec = sum(1 for k, x in tar.items() if k[0] == cd and "inicio" in x and "fim" not in x)
        num("Cho%sFalhasCurso" % nome_c, v(nec * fc, 2), "realizacoes em curso x fracao esperada de falhas")
        d = {"titulo": "$r$ = " + v(rcho, 1).replace("-", "\u2212"), "q_dir": (G < 1.).astype(float), "cho": 100. * pfc, "ai": [],
             "ymax": ymax, "N": N, "p": p, "cov": covp, "cont": contagens(G), "nome": nome_c, "pronto": pronto,
             "caso": cd, "nf": nf, "lo": lo, "hi": hi, "ic": ic, "alvo": ND}
        # a AI na banda deste r
        rsa = rs[ca]
        NA, opa = casos[ca]
        num("ChoAI%sAlvo" % nome_c, "%d" % NA, "casos.txt")
        num("ChoAI%sSemente" % nome_c, opa.get("semente", "--"), "casos.txt")
        fA = "campanha/%s/mc_*.csv" % ca
        ini_ai = len(rsa) >= 2
        fe_a, nb_a = feitos(os.path.join(K, ca), NA, Bl) if os.path.isdir(os.path.join(K, ca)) else (0, 1)
        flag("choai" + nome_c.lower() + "iniciado", ini_ai, fA)
        flag("choai" + nome_c.lower() + "pronto", fe_a == nb_a, "campanha/%s/mc_mapa_*.csv" % ca)
        flag("choai" + nome_c.lower() + "falhas", ini_ai and any(rsa[n]["fs"] < 1. for n in rsa),
             fA + " (alguma falha sorteada)")
        nbd, ig, dd = le_desloc(os.path.join(K, ca)) if os.path.isdir(os.path.join(K, ca)) else (0, True, {})
        if dd:
            fdd = "campanha/%s/mc_desloc_*.txt" % ca
            num("ChoAI%sElite" % nome_c, dd["elite"], fdd)
            num("ChoAI%sNivel" % nome_c, v(float(dd["nivel"]), 4), fdd)
            num("ChoAI%sPilFalhas" % nome_c, dd["falhas"], fdd)
            num("ChoAI%sMuDois" % nome_c, v(float(dd["mu2"]), 2), fdd)
            sd = [float(x) for x in dd["s"].split(",")]
            num("ChoAI%sS" % nome_c, "; ".join(v(x, 2) for x in sd[:2]), fdd + " (s_c; s_phi)",
                latex=r"\ensuremath{" + ";\\ ".join(tex(x, 2) for x in sd[:2]) + "}")
            txt("ChoAI%sIguais" % nome_c, "o mesmo" if ig else "\\textbf{deslocamentos diferentes}", fdd)
        else:
            for mac in ("Elite", "Nivel", "PilFalhas", "MuDois", "S", "Iguais"):
                txt("ChoAI%s%s" % (nome_c, mac), "--", "a AI ainda nao gravou o deslocamento")
        rot_ai = "AI na banda" + ("" if fe_a == nb_a else " (parcial)")
        if ini_ai and p > 0:
            ea = estima(rsa)
            dif, gan, _ = compara(ea, p, covp) if ea["pf"] > 0 else (math.nan, math.nan, math.nan)
            num("ChoAI%sN" % nome_c, "%d" % ea["N"], fA)
            num("ChoAI%sNf" % nome_c, "%d" % ea["nf"], fA)
            txt("ChoAI%sFalhas" % nome_c, falhas(ea["nf"]), fA)
            num("ChoAI%sPf" % nome_c, v(100. * ea["pf"], 2), fA + " (soma(w I) / N)")
            num("ChoAI%sCov" % nome_c, v(100. * ea["cov"], 1) if ea["pf"] > 0 else "--", fA)
            num("ChoAI%sLo" % nome_c, v(max(100. * ea["pf"] * (1 - 1.96 * ea["cov"]), 0.), 2) if ea["pf"] > 0 else "--",
                fA + " (P_f (1 - 1,96 CoV))")
            num("ChoAI%sHi" % nome_c, v(100. * ea["pf"] * (1 + 1.96 * ea["cov"]), 2) if ea["pf"] > 0 else "--",
                fA + " (P_f (1 + 1,96 CoV))")
            num("ChoAI%sDif" % nome_c, v(dif, 2) if ea["pf"] > 0 else "--", fA + " e " + fD + " (desvios)")
            num("ChoAI%sGanho" % nome_c, v(gan, 1) if ea["pf"] > 0 else "--", fA)
            num("ChoAI%sPeso" % nome_c, v(ea["w_media"], 3), fA)
            num("ChoAI%sEss" % nome_c, "%d" % round(ea["ess"]), fA)
            d["ai"].append((rot_ai, ea["q"], COR_AI))
            d["ea"] = ea
            d["dif"], d["gan"] = dif, gan
            d["cont_ai"] = contagens(ea["G"])
        else:
            for mac in ("N", "Nf", "Falhas", "Pf", "Cov", "Lo", "Hi", "Dif", "Ganho", "Peso", "Ess"):
                txt("ChoAI%s%s" % (nome_c, mac), "--", "a AI ainda nao tem realizacoes")
        p_ref = d["ea"]["pf"] if "ea" in d and d["ea"]["pf"] > 0 else p
        d["ganho_lin"] = 1. / ganho_ai(p_ref)[2] if p_ref > 0 else math.nan
        num("ChoGanhoLin%s" % nome_c, v(d["ganho_lin"], 0) if p_ref > 0 else "--",
            "formula, P_f da AI (ou do direto, sem a AI)")
        out["casos"].append(d)
        out[chave] = d
    # as 200 primeiras do r05_direto: o 'antes' do teste pareado da etapa 1
    ce = os.path.join(SRF, "cluster", "casos_etapa1.txt")
    c1 = le_casos(ce) if os.path.exists(ce) else {}
    npar = c1.get("r05_antigo", (1000, {}))[0]
    num("EtapaUmPar", "%d" % npar, "cluster/casos_etapa1.txt (r05_antigo)")
    num("EtapaUmChuva", "%d" % c1["chuva_antigo"][0] if "chuva_antigo" in c1 else "--",
        "cluster/casos_etapa1.txt (chuva_antigo)")
    mf0 = re.search(r"F0 = ([\d,]+)", texto(ce)) if os.path.exists(ce) else None
    txt("EtapaUmFzero", mf0[1].replace(",", "{,}") if mf0 else "--", "cluster/casos_etapa1.txt")
    r5 = rs["r05_direto"]
    antes = [r5[n]["fs"] for n in range(npar) if n in r5]
    num("ChoAntesN", "%d" % len(antes), "campanha/r05_direto/mc_*.csv (realizacoes 0..%d)" % (npar - 1))
    num("ChoAntesNf", "%d" % sum(1 for x in antes if x < 1.), "campanha/r05_direto/mc_*.csv (0..%d, FS < 1)"
        % (npar - 1))
    num("ChoAntesPerto", "%d" % sum(1 for x in antes if FAIXAS[2][0] <= x < FAIXAS[3][1]),
        "campanha/r05_direto/mc_*.csv (FS em [%g; %g))" % (FAIXAS[2][0], FAIXAS[3][1]))
    txt("ChoAntesFaixa", r"$[%s;\ %s)$" % (tex(FAIXAS[2][0], 0), tex(FAIXAS[3][1], 2)), "FAIXAS do script")
    # ---- o FORM (M = 10 e 20)
    forms = []
    mcrit = le_plano(os.path.join(PLANO, "README.md"), r"Parada: \|g\| < ([\d,]+) e mudança de β < ([\d,]+)",
                     "o criterio de parada da etapa 2")
    if mcrit:
        tg2, tb2 = float(mcrit[1].replace(",", ".")), float(mcrit[2].replace(",", "."))
        num("FormTolGNovo", mcrit[1], "../README.md (etapa 2)")
        num("FormTolBNovo", mcrit[2], "../README.md (etapa 2)")
    for nome_m, pasta in (("Dez", "form_M10_r05"), ("Vinte", "form_M20_r05")):
        F = os.path.join(H, pasta)
        fo = le_form(F)
        fF = pasta + ("/form.json" if fo and fo["pronto"] else "/form_parcial.json")
        flag("form" + nome_m.lower() + "pronto", fo is not None and fo["pronto"], pasta + "/form.json")
        if fo is None:
            continue
        forms.append(fo)
        M = fo["M"]
        num("Form%sM" % nome_m, "%d" % M, fF)
        num("Form%sVar" % nome_m, v(100. * fo["var"], 1), pasta + "/it0_centro_parametros.txt")
        num("Form%sAvalIt" % nome_m, "%d" % (4 * M + 1), "4M + 1")
        num("Form%sIt" % nome_m, "%d" % len(fo["hist"]), fF + " (iteracoes completas)")
        num("Form%sUltIt" % nome_m, "%d" % fo["hist"][-1]["it"], fF)
        num("Form%sAval" % nome_m, "%d" % fo["aval"], fF)
        num("Form%sItAtual" % nome_m, "%d" % fo["em_curso"], pasta + "/it*_0_0.csv")
        num("Form%sFeitas" % nome_m, "%d" % fo["feitas"], pasta + "/it%d_*_0_0.csv" % fo["em_curso"])
        num("Form%sBeta" % nome_m, v(fo["beta"], 2), fF)
        num("Form%sG" % nome_m, v(fo["g"], 4), fF + " (g da ultima iteracao)")
        num("Form%sPf" % nome_m, v(100. * fo["pf"], 2), fF + " (Phi(-beta))")
        num("Form%sBetaProx" % nome_m, v(fo["beta_prox"], 2), fF + " (proximo ponto do HL-RF)")
        num("Form%sPfProx" % nome_m, v(100. * norm_cdf(-fo["beta_prox"]), 2), fF + " (Phi(-beta) do proximo ponto)")
        a2 = fo["alfa2"]
        num("Form%sAlfaC" % nome_m, v(float(a2[:M].sum()), 2), fF + " (soma dos alfa^2 de c')")
        num("Form%sAlfaPhi" % nome_m, v(float(a2[M:].sum()), 2), fF + " (soma dos alfa^2 de phi')")
        top = sorted(range(2 * M), key=lambda i: -a2[i])[:3]
        txt("Form%sTop" % nome_m, ", ".join((r"$c'_{%d}$" if i < M else r"$\varphi'_{%d}$") % (i % M) + " ("
                                            + tex(a2[i], 2) + ")" for i in top), fF + " (maiores alfa^2)")
        txt("Form%sConvergiu" % nome_m, ("convergiu em %d iterações" % len(fo["hist"]) if fo["convergiu"] else
                                          "parou no limite de %d iterações" % len(fo["hist"]))
            if fo["pronto"] else "--", fF)
        txt("Form%sSobe" % nome_m, "sobe" if fo["beta_prox"] > fo["beta"] else "desce", fF)
        db = abs(fo["hist"][-1]["beta"] - fo["hist"][-2]["beta"]) if len(fo["hist"]) > 1 else math.inf
        num("Form%sDBeta" % nome_m, v(db, 3), fF + " (mudanca de beta na ultima iteracao)")
        ps = fo["passos"][2:] if len(fo["passos"]) > 2 else fo["passos"][-1:]
        print("FORM M = %d: passos em u por iteracao %s" % (M, ["%.3f" % x for x in fo["passos"]]))
        num("Form%sPasso" % nome_m, v(fo["passos"][-1], 2) if fo["passos"] else "--",
            fF + " (|u_k - u_k-1| na ultima iteracao)")
        num("Form%sPassoMin" % nome_m, v(min(ps), 2) if ps else "--", fF + " (menor |du| depois da iteracao 2)")
        flag("form%spassoacima" % nome_m.lower(), bool(ps) and min(ps) > fo["tol_p"],
             fF + " e form_kl.py (todos os passos depois da iteracao 2 acima da tolerancia)")
        if mcrit:
            txt("Form%sCritNovo" % nome_m, "já passa" if abs(fo["g"]) < tg2 and db < tb2 else "ainda não passa",
                fF + " e ../README.md (|g| = %.4f, mudanca de beta = %.4f)" % (fo["g"], db))
        else:
            txt("Form%sCritNovo" % nome_m, "--", "o plano nao tem mais o criterio")
        num("Form%sParalelos" % nome_m, "%d" % fo["p"]["paralelos"], fF)
        num("Form%sItMax" % nome_m, "%d" % fo["p"]["iteracoes"], fF)
        fo["nome"] = nome_m
    fo10 = forms[0]
    num("FormH", v(fo10["p"]["h"], 1), "form_M10_r05/form_parcial.json (h)")
    num("FormTolG", v(fo10["tol_g"], 3), "form_kl.py")
    num("FormTolP", v(fo10["tol_p"], 2), "form_kl.py")
    fs0, c0, phi0 = fo10["centro"]
    num("FormCentroFS", v(fs0), "form_M10_r05/it0_centro_0_0.csv")
    num("FormCentroC", v(c0, 2), "form_M10_r05/it0_centro_0_0.csv (c_med)")
    num("FormCentroPhi", v(phi0, 1), "form_M10_r05/it0_centro_0_0.csv (phi_med)")
    # c' pesa pouco: o xi de c' tambem mexe no phi' (r < 0). Regra da cadeia no 1o termo, na ultima iteracao de M = 10:
    # d FS / d xi_c = A + r B e d FS / d xi_phi = raiz(1 - r^2) B
    M = fo10["M"]
    gr = fo10["grad"]
    Bt = gr[M] / math.sqrt(1. - rc5 * rc5)
    At = gr[0] - rc5 * Bt
    num("FormGradC", v(gr[0], 3), "form_M10_r05: gradiente em c'_0")
    num("FormGradPhi", v(gr[M], 3), "form_M10_r05: gradiente em phi'_0")
    num("FormViaC", v(At, 3), "regra da cadeia: o efeito do xi de c'_0 pelo c'")
    num("FormViaPhi", v(rc5 * Bt, 3), "regra da cadeia: o mesmo xi pelo phi' (r B)")
    # FORM contra o direto ate agora (r = -0,5)
    p5 = out["r05"]["p"]
    num("ChoBetaDireto", v(norm_inv(p5), 2) if p5 > 0 else "--", "-Phi^-1(P_f do direto ate agora)")
    if len(forms) > 1 and forms[1]["pronto"] and forms[0]["pronto"]:
        txt("FormCompara", "maior" if forms[1]["beta"] > forms[0]["beta"] else "menor", "form_M*_r05/form.json")
    else:
        txt("FormCompara", "--", "os dois FORM ainda nao terminaram")
    # ---- a AI com o ponto de projeto do FORM (cadeia_form.sh)
    cf = texto(os.path.join(H, "cadeia_form.sh"))
    mcf = re.search(r"r05_ai_form\s+\d+\s+(\d+)\s+.*?semente=(\d+)", cf)
    num("AIFormAlvo", mcf[1], "cadeia_form.sh")
    num("AIFormSemente", mcf[2], "cadeia_form.sh")
    seqs = [(int(a), int(b)) for a, b in re.findall(r"for p in \$\(seq (\d+) (\d+)\)", cf)]
    num("AIFormProc", "%d" % (seqs[0][1] - seqs[0][0] + 1), "cadeia_form.sh (processos no comeco)")
    num("AIFormProcMais", "%d" % (seqs[1][1] - seqs[1][0] + 1), "cadeia_form.sh (processos depois do rodar.sh)")
    AF = os.path.join(H, "campanha_form", "r05_ai_form")
    rsf = le_rodada(AF)[0] if os.path.isdir(AF) else {}
    rsf.pop(-1, None)
    naf = int(mcf[1])
    fe_f, nb_f = feitos(AF, naf, int(re.search(r"BLOCO=(\d+)", cf)[1])) if os.path.isdir(AF) else (0, 1)
    flag("aiforminiciado", len(rsf) >= 2, "campanha_form/r05_ai_form/mc_*.csv")
    flag("aiformfalhas", len(rsf) >= 2 and any(rsf[n]["fs"] < 1. for n in rsf),
         "campanha_form/r05_ai_form/mc_*.csv (alguma falha sorteada)")
    flag("aiformpronto", fe_f == nb_f and os.path.exists(os.path.join(H, "resultado_form.txt")),
         "campanha_form/r05_ai_form/mc_mapa_*.csv e resultado_form.txt")
    nbf, igf, dfo = le_desloc(AF) if os.path.isdir(AF) else (0, True, {})
    if dfo.get("mu2"):  # o deslocamento que a AI gravou (ai=dado: o ponto de projeto termo a termo)
        num("AIFormMuDois", v(float(dfo["mu2"]), 2), "campanha_form/r05_ai_form/mc_desloc_*.txt")
        num("AIFormVarPeso", mil(math.exp(float(dfo["mu2"])) - 1.), "exp(|mu|^2) - 1")
    elif fo10["pronto"]:
        num("AIFormMuDois", v(float((fo10["ustar"] ** 2).sum()), 2), "form_M10_r05/form.json (|u*|^2)")
        num("AIFormVarPeso", mil(math.exp(float((fo10["ustar"] ** 2).sum())) - 1.), "exp(|u*|^2) - 1")
    else:
        num("AIFormMuDois", v(fo10["beta_prox"] ** 2, 2), "form_M10_r05/form_parcial.json (beta previsto ao quadrado)")
        num("AIFormVarPeso", mil(math.exp(fo10["beta_prox"] ** 2) - 1.), "exp(beta^2) - 1, previsto")
    r5d = out["r05"]
    if len(rsf) >= 2 and r5d["p"] > 0:
        ef = estima(rsf)
        dif, gan, _ = compara(ef, r5d["p"], r5d["cov"]) if ef["pf"] > 0 else (math.nan, math.nan, math.nan)
        fA = "campanha_form/r05_ai_form/mc_*.csv"
        num("AIFormN", "%d" % ef["N"], fA)
        num("AIFormNf", "%d" % ef["nf"], fA)
        txt("AIFormFalhas", falhas(ef["nf"]), fA)
        num("AIFormPf", v(100. * ef["pf"], 2), fA)
        num("AIFormCov", v(100. * ef["cov"], 1) if ef["pf"] > 0 else "--", fA)
        num("AIFormLo", v(max(100. * ef["pf"] * (1 - 1.96 * ef["cov"]), 0.), 2) if ef["pf"] > 0 else "--",
            fA + " (P_f (1 - 1,96 CoV))")
        num("AIFormHi", v(100. * ef["pf"] * (1 + 1.96 * ef["cov"]), 2) if ef["pf"] > 0 else "--",
            fA + " (P_f (1 + 1,96 CoV))")
        num("AIFormDif", v(dif, 2) if ef["pf"] > 0 else "--", fA + " e campanha/r05_direto (desvios)")
        num("AIFormGanho", v(gan, 1) if ef["pf"] > 0 else "--", fA)
        num("AIFormPeso", v(ef["w_media"], 3), fA)
        num("AIFormEss", "%d" % round(ef["ess"]), fA)
        r5d["ai"].append(("AI no ponto de projeto" + ("" if FLAGS["aiformpronto"] else " (parcial)"), ef["q"],
                          COR_TERCEIRA))
        out["ef"], out["ef_dif"], out["ef_gan"] = ef, dif, gan
    else:
        for mac in ("N", "Nf", "Falhas", "Pf", "Cov", "Lo", "Hi", "Dif", "Ganho", "Peso", "Ess"):
            txt("AIForm" + mac, "--", "a AI com o ponto de projeto ainda nao tem realizacoes")
    out["forms"] = forms
    out["niveis"] = [("tri1 + %d níveis" % (len(nv["TriUm"]) - 1), nv["TriUm"]),
                     ("tri2 + %d níveis" % (len(nv["TriDois"]) - 1), nv["TriDois"])]
    return out


# ---------------------------------------------------------------- C) o contexto citado e o balanco
def contexto():
    mr = le_macros(os.path.join(REL, "numeros_rebaixamento.tex"))
    fr = "numeros_rebaixamento.tex (relatorio_rebaixamento.pdf)"
    rb = {}
    for nome in ("NDireto", "NfDireto", "PfDireto", "CovPfDireto", "NBanda", "PfBanda", "CovPfBanda", "DifBanda",
                 "GanhoBanda", "GanhoLinear", "FracGanhoBanda", "MuDoisBanda", "SC", "SPhi", "SK", "KModos", "NModos",
                 "PfModos", "CovPfModos", "DifModos", "GanhoModos", "KNPar", "KPfFixo", "KPfAleat", "KMcNemar"):
        rb[nome] = cita(mr, nome, "Rb" + nome, fr)
    ms = le_macros(os.path.join(REL, "numeros_solos.tex"))
    fs = "numeros_solos.tex (relatorio_solos.pdf)"
    for nome in ("BuracoAbaixo", "BuracoAcima", "BuracoFaixaAbaixo", "BuracoFaixaAcima", "BuracoPontos", "BuracoFalta",
                 "NCinquenta"):
        cita(ms, nome, "So" + nome, fs)
    return rb


def vmd(pd, cd, pa, ca):
    """o vies minimo detectavel (em % da P_f do direto) pelo teste de dois desvios: o vies que ele pega metade das
    vezes, 1,96 desvio da diferenca (P_f e CoV do direto e da AI)"""
    return 100. * 1.96 * math.sqrt((pd * cd) ** 2 + (pa * ca) ** 2) / pd


def balanco(ch, co, rb):
    """a tabela do balanco das validacoes da AI e as frases que dependem dela"""
    rows, difs, gans, fracs = [], [], [], []
    fr = "numeros_rebaixamento.tex (P_f e CoV citados)"
    num("RbVmdBanda", v(vmd(rb["PfDireto"], rb["CovPfDireto"] / 100., rb["PfBanda"], rb["CovPfBanda"] / 100.), 0), fr)
    num("RbVmdModos", v(vmd(rb["PfDireto"], rb["CovPfDireto"] / 100., rb["PfModos"], rb["CovPfModos"] / 100.), 0), fr)
    num("ChVmd", v(vmd(ch["pP"], ch["covP"], ch["pA"], ch["covA"]), 0), "chuva_ai (P_f e CoV do piloto e da AI)")
    rows.append(r"rebaixamento & banda & \RbPfDireto{} (\RbNDireto) & \RbPfBanda{} (\RbNBanda) & "
                r"\RbDifBanda & \RbVmdBanda & \RbGanhoBanda & \RbGanhoLinear \\")
    rows.append(r"rebaixamento & \RbKModos{} modos & \RbPfDireto{} (\RbNDireto) & \RbPfModos{} (\RbNModos) "
                r"& \RbDifModos & \RbVmdModos & \RbGanhoModos & \RbGanhoLinear \\")
    rows.append(r"chuva (\ChSolo) & banda (elite) & \ChPPf{} (\ChNP) & \ChAPf{} (\ChNA) & \ChDif & \ChVmd & "
                r"\ChGanho & \ChGanhoLinear \\")
    difs += [rb["DifBanda"], rb["DifModos"], ch["dif"]]
    gans += [("banda", rb["GanhoBanda"], rb["FracGanhoBanda"] / 100.), ("banda", ch["gan"], ch["frac"])]
    for d in co["casos"]:
        n = d["nome"]
        dir_ = r"\Cho%sPf{} (\Cho%sN)%s" % (n, n, "" if d["pronto"] else r"$^\dagger$")
        if "ea" in d:
            ai = r"\ChoAI%sPf{} (\ChoAI%sN)%s" % (n, n, "" if FLAGS["choai%spronto" % n.lower()] else r"$^\dagger$")
            num("ChoAI%sVmd" % n, v(vmd(d["p"], d["cov"], d["ea"]["pf"], d["ea"]["cov"]), 0)
                if d["ea"]["pf"] > 0 and d["p"] > 0 else "--", "campanha/%s (P_f e CoV)" % n)
            rows.append(r"Cho, $r = \ChoR%s$ & banda (elite) & %s & %s & \ChoAI%sDif & \ChoAI%sVmd & \ChoAI%sGanho & "
                        r"\ChoGanhoLin%s \\" % (n, dir_, ai, n, n, n, n))
            if FLAGS["choai%spronto" % n.lower()] and d["ea"]["pf"] > 0:
                difs.append(d["dif"])
                gans.append(("banda", d["gan"], d["gan"] / d["ganho_lin"]))
        else:
            rows.append(r"Cho, $r = \ChoR%s$ & banda (elite) & %s & \multicolumn{4}{l}{%s} & \ChoGanhoLin%s \\"
                        % (n, dir_, r"não rodou no PC" if FLAGS["pcparou"] else "pendente", n))
    d = co["r05"]
    dir_ = r"\ChoCincoPf{} (\ChoCincoN)%s" % ("" if d["pronto"] else r"$^\dagger$")
    if "ef" in co:
        num("AIFormVmd", v(vmd(d["p"], d["cov"], co["ef"]["pf"], co["ef"]["cov"]), 0)
            if co["ef"]["pf"] > 0 and d["p"] > 0 else "--", "campanha_form/r05_ai_form (P_f e CoV)")
        rows.append(r"Cho, $r = \ChoRCinco$ & ponto de projeto & %s & \AIFormPf{} (\AIFormN)%s & \AIFormDif & "
                    r"\AIFormVmd & \AIFormGanho & \ChoGanhoLinCinco \\" % (dir_, "" if FLAGS["aiformpronto"] else
                                                                            r"$^\dagger$"))
        if FLAGS["aiformpronto"] and co["ef"]["pf"] > 0:
            difs.append(co["ef_dif"])
    else:
        rows.append(r"Cho, $r = \ChoRCinco$ & ponto de projeto & %s & \multicolumn{4}{l}{%s} & \ChoGanhoLinCinco \\"
                    % (dir_, r"não rodou no PC" if FLAGS["pcparou"] else "pendente"))
    txt("BalLinhas", "\n".join(rows), "as macros das linhas")
    num("BalNCasos", "%d" % len(difs), "casos completos do balanco")
    num("BalMaxDif", v(max(abs(x) for x in difs), 2), "o maior |diferenca| dos casos completos (desvios)")
    txt("BalVies", "todas abaixo de dois desvios" if max(abs(x) for x in difs) < 2 else
        "uma delas acima de dois desvios", "casos completos do balanco")
    gb = [g for t, g, _ in gans if t == "banda"]
    fb = [f for t, _, f in gans if t == "banda"]
    num("BalGanhoMin", v(min(gb), 1), "ganhos medidos na banda")
    num("BalGanhoMax", v(max(gb), 1), "ganhos medidos na banda")
    num("BalFracMin", v(100. * min(fb), 0), "ganho medido / ganho linear, na banda")
    num("BalFracMax", v(100. * max(fb), 0), "ganho medido / ganho linear, na banda")


def buraco(ch, co):
    """as contagens em volta de FS = 1 (faixas de 0,01)"""
    fx = [r"$[%s;\ %s)$" % (tex(a, 0 if a == 1 else 2), tex(b, 0 if b == 1 else 2)) for a, b in FAIXAS]
    txt("BuracoFaixa", fx[1], "FAIXAS do script")
    txt("BuracoCab", " & ".join(fx), "FAIXAS do script")
    rows = [(r"chuva: piloto (direto)", ch["nP"], ch["contP"]), (r"chuva: AI (sem peso)", ch["NA"], ch["contA"])]
    for d in co["casos"]:
        rows.append((r"Cho, $r = \ChoR%s$: direto%s" % (d["nome"], "" if d["pronto"] else " (parcial)"), d["N"],
                     d["cont"]))
        if "cont_ai" in d:
            rows.append((r"Cho, $r = \ChoR%s$: AI (sem peso)" % d["nome"], d["ea"]["N"], d["cont_ai"]))
    c5 = co["r05"]["cont"]
    num("ChoBuracoAbaixo", "%d" % c5[1], "campanha/r05_direto/mc_*.csv (FS na faixa %s)" % (FAIXAS[1],))
    txt("ChoBuracoViz", "%d e %d" % (c5[0], c5[2]), "campanha/r05_direto/mc_*.csv (faixas %s e %s)"
        % (FAIXAS[0], FAIXAS[2]))
    flag("choburaco", c5[1] < min(c5[0], c5[2]), "campanha/r05_direto/mc_*.csv (faixa abaixo de 1 menor que as vizinhas)")
    txt("BuracoLinhas", "\n".join(r"%s & %d & %s \\" % (r, n, " & ".join("%d" % c for c in cs)) for r, n, cs in rows),
        "contagens de FS nas faixas %s" % FAIXAS)
    print("buraco (faixas %s): %s" % (FAIXAS, [(r, n, cs) for r, n, cs in rows]))


def cho_linhas(co):
    """as linhas da tabela do Monte Carlo do Cho: o direto e as AI (pendente, parcial ou pronta)"""
    rows = []
    for d in co["casos"]:
        n = d["nome"]
        rows.append(r"$r = \ChoR%s$ & direto & \Cho%sN{} de \Cho%sAlvo & \Cho%sNf & \Cho%sPf & \Cho%sLo--\Cho%sHi & "
                    r"\Cho%sCov & %s \\" % (n, n, n, n, n, n, n, n, r"\ChoPfCinco" if n == "Cinco" else
                                               r"$\approx$\ChoPfSete"))
        nao = r"não rodou no PC (vai para o fermi)"
        metodos = [("AI na banda", "ChoAI" + n, "choai%s" % n.lower(),
                    nao if FLAGS["pcparou"] else r"pendente (começa por volta das \ChoInicioAI)")]
        if n == "Cinco":
            metodos.append(("AI no ponto de projeto", "AIForm", "aiform", nao if FLAGS["pcparou"] else
                            r"pendente (%s)" % ("começando" if FLAGS["formdezpronto"] else
                                                r"depois do FORM com $M = \FormDezM$")))
        for rot, pre, fl, pend in metodos:
            if FLAGS[fl + "iniciado"]:
                rows.append(r" & %s%s & \%sN{} de \%sAlvo & \%sNf & \%sPf & \%sLo--\%sHi & \%sCov & \\"
                            % (rot, "" if FLAGS[fl + "pronto"] else " (parcial)", pre, pre, pre, pre, pre, pre, pre))
            else:
                rows.append(r" & %s & \multicolumn{6}{l}{%s} \\" % (rot, pend))
        rows.append(r"\addlinespace")
    txt("ChoLinhas", "\n".join(rows[:-1]), "as macros das linhas")


def form_linhas(forms):
    """as linhas da tabela do FORM: g e beta por iteracao, M = 10 e 20 lado a lado"""
    nit = max(len(f["hist"]) for f in forms)
    rows = []
    for it in range(nit):
        cel = []
        for f in forms:
            h = f["hist"][it] if it < len(f["hist"]) else None
            cel += [tex(h["g"], 4), tex(h["beta"], 4)] if h else ["", ""]
        rows.append("%d & %s \\\\" % (it, " & ".join(cel)))
    txt("FormLinhas", "\n".join(rows), "form_M*_r05/form*.json (historico)")


def copia_resultados(H, co):
    """os resumos pequenos do Cho2010 para resultados/cho2010 (so os que existem; o form_parcial.json quando o FORM
    nao terminou) e o resumo_previa.txt com o direto do PC"""
    dst = os.path.join(AQUI, "resultados", "cho2010")
    os.makedirs(dst, exist_ok=True)
    arqs = ["resultado.txt", "resultado_form.txt", "casos_form.txt"]
    for pf in ("form_M10_r05", "form_M20_r05"):
        arqs += [pf + "/form.log", pf + "/form.json"]
        if not os.path.exists(os.path.join(H, pf, "form.json")):
            arqs.append(pf + "/form_parcial.json")
    with open(os.path.join(dst, "resumo_previa.txt"), "w") as f:
        f.write("# Cho2010 no PC, SRM antigo: o Monte Carlo direto ate %s (%s)\n"
                % (AGORA.strftime("%d/%m/%Y %H:%M"), "parcial" if not FLAGS["chofim"] else "completo"))
        f.write("# gerado por figuras_etapa0.py a partir de %s/campanha/r0*_direto/mc_*.csv\n" % H)
        f.write("%-11s %6s %6s %7s %8s %8s %7s %s\n" % ("caso", "N", "alvo", "falhas", "P_f(%)", "CoV(%)", "IC",
                                                         "IC 95% (%)"))
        for d in co["casos"]:
            f.write("%-11s %6d %6d %7d %8.3f %8.1f %7s %.3f a %.3f\n" % (
                d["caso"], d["N"], d["alvo"], d["nf"], 100. * d["p"], 100. * d["cov"] if d["p"] > 0 else float("nan"),
                d["ic"], 100. * d["lo"], 100. * d["hi"]))
    print("resumo: %s" % os.path.join(dst, "resumo_previa.txt"))
    for rel_ in arqs:
        src = os.path.join(H, rel_)
        if os.path.exists(src):
            alvo = os.path.join(dst, rel_)
            os.makedirs(os.path.dirname(alvo), exist_ok=True)
            shutil.copy2(src, alvo)
            print("copia: %s -> %s" % (src, alvo))


def main(argv):
    C = os.path.abspath(argv[1]) if len(argv) > 1 else CHUVA_PADRAO
    H = os.path.abspath(argv[2]) if len(argv) > 2 else CHO_PADRAO
    B = os.path.abspath(argv[3]) if len(argv) > 3 else BUILD_PADRAO
    pasta = os.path.join(AQUI, "figs")
    os.makedirs(pasta, exist_ok=True)
    txt("HoraDados", AGORA.strftime("%d/%m/%Y, ") + hora(AGORA).split(" de ")[0], "relogio do PC ao rodar o script")
    mw = le_plano(os.path.join(PLANO, "README.md"), r"Com menos de ~(\d+) falhas no direto", "o limite do IC exato")
    if mw:
        MIN_WILSON[0] = int(mw[1])
    num("MinWilson", "%d" % MIN_WILSON[0], "../README.md (criterios gerais)" if mw else "padrao do script")
    rb = contexto()
    ch = chuva(C, B)
    co = cho(H, ch)
    balanco(ch, co, rb)
    buraco(ch, co)
    form_linhas(co["forms"])
    cho_linhas(co)
    copia_resultados(H, co)
    # ---- figuras
    fig_chuva(pasta, ch["ts"], ch["pil"], ch["ai_t"], ch["ep"], ch["ea"], ch["nP"])
    fig_amostra(pasta, ch["ep"], ch["ea"], ch["nivel"], ch["f24"], ch["f8v"], ch["ln200"], ch["ln8"], ch["pA"],
                ch["ratio_lin"])
    fig_cho(pasta, co["niveis"], co["bishop"], co["casos"])
    p5 = co["r05"]["p"]
    fig_form(pasta, co["forms"], co["beta_cho"], norm_inv(p5) if p5 > 0 else None)
    # ---- numeros_etapa0.tex
    arq = os.path.join(AQUI, "numeros_etapa0.tex")
    with open(arq, "w") as f:
        f.write("%% gerado por figuras_etapa0.py em %s: nao edite (rode o script de novo)\n"
                % AGORA.strftime("%d/%m/%Y %H:%M"))
        for nome, val in FLAGS.items():
            f.write("\\newif\\if%s \\%s%s\n" % (nome, nome, "true" if val else "false"))
        for macro, valor, latex, fonte in NUMEROS:
            if latex is not None:
                f.write("\\newcommand{\\%s}{%s} %% %s\n" % (macro, latex.replace("\n", "\n  "),
                                                           fonte.replace("%", r"\%")))
    print("numeros: %s" % arq)
    print("\n% ---- numeros do texto (macro: valor  [arquivo])")
    for macro, valor, latex, fonte in NUMEROS:
        print("%-22s %-34s [%s]" % (macro, valor.replace("\n", " | ")[:200], fonte))


if __name__ == "__main__":
    main(sys.argv)
