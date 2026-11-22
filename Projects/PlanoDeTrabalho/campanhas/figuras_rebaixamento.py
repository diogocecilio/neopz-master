#!/usr/bin/env python3
"""figuras_rebaixamento.py: as figuras, as tabelas e os numeros do relatorio_rebaixamento.tex: o Monte Carlo do
rebaixamento rapido de Vargas Ceron et al. (2025), Secao 6.1, com o transiente (SlopeDrawdownMonteCarlo, o
rebaixamento_mc.cpp), rodado no PC em 14/11/2026, a amostragem por importancia (AI) testada nele e o teste do k fixo.
So le arquivos; nao roda nenhum programa e nao escreve nas pastas das rodadas:
  - as rodadas do PC (~/projects/runs/ceron_transiente): direto/ (o Monte Carlo direto e o deterministico), ai/ (AI na
    direcao da banda), ai_modos/ (AI nos primeiros modos), covk0/ (k fixo) e verificacao/ (os pesos sem elementos
    finitos);
  - a campanha antiga do rebaixamento (percolacao permanente, cluster): ~/projects/runs/_arquivo/2026-11-11_ceron_agua_antiga/fermi_ceron/ceron_mc;
  - o custo da chuva no build (mc/serial_I24.log, para a campanha da chuva no fermi);
  - os tres solos baixados do fermi (CAMPANHA=solos, ~/projects/runs/fermi_chuva_solos), se a pasta existe: as
    realizacoes, as falhas e a P_f finais da Secao 6.1 (as regras do figuras_solos.py, que faz o relatorio_solos.tex);
  - as fontes desta pasta (rebaixamento_mc.cpp, DrawdownTransient.h e cluster/).

Escreve:
  figs/fig_rb_*.pdf          as figuras
  numeros_rebaixamento.tex   um \\newcommand para cada numero do texto e das tabelas e os \\newif das rodadas que ainda
                             rodam (\\ifaipronto, \\ifmodospronto, \\ifkfixopronto e os \\if...iniciado)
e imprime cada numero com o arquivo de onde sai. Rodado de novo depois que as rodadas acabam, o relatorio mostra os
resultados no lugar de "pendente".

Os unicos numeros digitados aqui sao: os do artigo do Ceron (Gamma deterministico da Fig. 24, as S = 10 000 simulacoes
da Secao 6.1, a linha s = 1 da Tabela 4 e a linha CoV(k_v) = 0 da Tabela 3), o fator 2,6 entre um nucleo do fermi e o
do PC (medido antes desta campanha) e o status dos tres solos no fermi, colado pelo usuario as 20h24 (saida do
status.sh, que nao esta neste PC); com a campanha dos solos baixada, a Secao 6.1 usa os arquivos dela
(\\ifsolosbaixados).

Uso: python3 -I figuras_rebaixamento.py [rodadas do PC] [campanha antiga] [build do SlopeRandomFieldRefinement]
                                         [tres solos baixados]
     (padrao: ~/projects/runs/ceron_transiente, ~/projects/runs/_arquivo/2026-11-11_ceron_agua_antiga/fermi_ceron/ceron_mc,
     ../../../../neopz-master-build-release/Projects/SlopeRandomFieldRefinement e ~/projects/runs/fermi_chuva_solos)"""

import csv
import datetime as dt
import glob
import math
import os
import re
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.ticker  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from scipy import stats  # noqa: E402
from scipy.optimize import brentq  # noqa: E402

AQUI = os.path.dirname(os.path.abspath(__file__))
NEOPZ = os.path.dirname(os.path.dirname(os.path.dirname(AQUI)))  # AQUI: Projects/PlanoDeTrabalho/campanhas
PROJ = os.path.join(NEOPZ, "Projects", "SlopeRandomFieldRefinement")  # as fontes do programa
RUNS = os.path.expanduser("~/projects/runs")
R_PADRAO = os.path.join(RUNS, "ceron_transiente")
ANTIGA_PADRAO = os.path.join(RUNS, "_arquivo", "2026-11-11_ceron_agua_antiga", "fermi_ceron", "ceron_mc")
BUILD_PADRAO = os.path.join(os.path.dirname(NEOPZ), "neopz-master-build-release", "Projects",
                            "SlopeRandomFieldRefinement")
SOLOS_PADRAO = os.path.join(RUNS, "fermi_chuva_solos")  # a campanha dos tres solos baixada pelo cluster/baixar.sh
DATA_RELATORIO = dt.date(2026, 11, 14)

# ---- o que e digitado
# Vargas Ceron et al. (2025): Gamma deterministico (Fig. 24); Secao 6.1: S = 10 000 simulacoes na referencia; Tabela 4,
# s = 1 (a referencia): media, desvio e CoV de Gamma, P_f e CoV(P_f), e P_f em s = 400; Tabela 3, CoV(k_v) = 0: media de
# Gamma e P_f
CERON = {"gamma_det": 1.336, "n_s1": 10000, "media_s1": 1.353, "dp_s1": 0.318, "cov_s1": 23.5, "pf_s1": 11.50,
         "covpf_s1": 2.77, "pf_s400": 21.93, "media_k0": 1.329, "pf_k0": 11.21}
FERMI_LENTO = 2.6  # um nucleo do fermi e 2,6 vezes mais lento que o do PC (medido antes desta campanha)
# os tres solos do Rahardjo no fermi (CAMPANHA=solos): status.sh colado pelo usuario as 20h24; os minimos do
# deterministico (FS de 6 em 6 h) e a previsao do fim vem da mesma mensagem
SOLOS = {"hora": (20, 24), "inicio": (19, 21), "nos": 4, "previsao_h": (17, 20),
         "casos": {"f10": {"n": 63, "nf": 61, "det": 0.832, "tdet": "12 a 24~h"},
                   "f50": {"n": 81, "nf": 14, "det": 1.035, "tdet": "24~h"},
                   "f100": {"n": 75, "nf": 0, "det": 1.881, "tdet": "36~h"}}}
FONTE_SOLOS = "status.sh no fermi as 20h24 (digitado)"
# as macros da campanha dos solos baixada (com '--' quando a pasta nao existe)
MACROS_SOLOS_BAIXADOS = ["Solos%s%s" % (a, b) for a in ("Validas", "Completas", "Excluidas", "ExcluidasFalham",
                                                        "Faltam", "Nf", "PfFinal", "PfFinalLo", "PfFinalHi", "DetFinal",
                                                        "TdetFinal") for b in ("Dez", "Cinquenta", "Cem")] + \
    ["SolosValidasTotal", "SolosFaltamTotal", "SolosFaltamLista", "SolosInicioLog", "SolosUltimo", "SolosRelogioFinal",
     "SolosCpuFinal", "SolosMinMedido", "SolosRodando", "SolosTarefas", "SolosPrevFimFinal"]

# ---- estilo (o do figuras_relatorio.py)
PALETA = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
TINTA, TINTA2, MUDO, GRADE, SOLO = "#0b0b0b", "#52514e", "#898781", "#dddcd8", "#f3f2ee"
COR_DIRETO, COR_BANDA, COR_MODOS = PALETA[0], PALETA[1], PALETA[2]
plt.rcParams.update({
    "font.size": 8.5, "axes.labelsize": 8.5, "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "legend.fontsize": 7,
    "axes.edgecolor": TINTA2, "axes.labelcolor": TINTA, "xtick.color": TINTA2, "ytick.color": TINTA2,
    "axes.grid": True, "grid.color": GRADE, "grid.linewidth": 0.5, "lines.linewidth": 1.4, "legend.frameon": False,
    "savefig.dpi": 200, "font.family": "DejaVu Sans", "mathtext.fontset": "dejavusans",
})
LARGURA = 6.5  # polegadas: a largura do texto do relatorio (A4 com margens de 2,2 cm)

NUMEROS = []  # (macro, valor impresso, valor no .tex, arquivo)
FLAGS = {}    # nome do \newif -> verdadeiro ou falso


# ---------------------------------------------------------------- numeros
def v(x, casas=4):
    """numero com virgula decimal (sem o sinal de -0,0)"""
    t = "%.*f" % (casas, x)
    if t.startswith("-") and float(t) == 0.:
        t = t[1:]
    return t.replace(".", ",")


def mil(n):
    """inteiro com o espaco fino do LaTeX a partir de 10 000 (como no relatorio_chuva.tex)"""
    n = int(round(n))
    return "{:,}".format(n).replace(",", r"\,") if n >= 10000 else "%d" % n


def num(macro, texto, fonte, latex=None):
    """registra um numero do texto: \\macro no .tex (virgula decimal como {,}, o sinal de menos em modo matematico) e o
    valor impresso com o arquivo de onde sai. latex: o valor no .tex, quando nao e um numero simples"""
    assert re.fullmatch(r"[A-Za-z]+", macro), macro
    assert macro not in {m for m, _, _, _ in NUMEROS}, "macro repetida: " + macro
    if latex is None:
        latex = re.sub(r"(?<!\\),", "{,}", texto)  # a virgula decimal (o \, do milhar fica)
        if re.match(r"-\d", latex):
            latex = r"\ensuremath{-}" + latex[1:]
    NUMEROS.append((macro, texto.replace(r"\,", " "), latex, fonte))
    return texto


def txt(macro, texto, fonte):
    """um texto que depende dos dados (uma conclusao, um horario): vai para o .tex como esta"""
    return num(macro, texto, fonte, latex=texto)


def p_rel(p, casas=2, aprox=False, piso=0.001):
    """um valor-p com a relacao, para o modo matematico ($p \\Macro$): '= 0{,}51', '\\approx 0{,}76' ou '< 0{,}001'.
    Devolve o texto impresso e o do .tex"""
    if p < piso:
        return "< " + v(piso, 3), "< " + v(piso, 3).replace(",", "{,}")
    t = v(p, casas)
    return ("~ " if aprox else "= ") + t, (r"\approx " if aprox else "= ") + t.replace(",", "{,}")


def acima_abaixo(x, casas=1):
    """uma variacao relativa em %: '2{,}4\\% abaixo' ou '1{,}3\\% acima'"""
    return r"%s\%% %s" % (v(abs(x), casas).replace(",", "{,}"), "abaixo" if x < 0 else "acima")


def lilliefors(x, B=2000, semente=1):
    """Kolmogorov-Smirnov da normal com a media e o desvio da propria amostra: D e o p de Lilliefors por bootstrap
    parametrico (B amostras normais do mesmo tamanho, cada uma com a media e o desvio dela; o D nao depende da media nem
    do desvio verdadeiros). p = (1 + numero de D sorteados >= D) / (B + 1)"""
    n = len(x)
    D = stats.kstest(x, "norm", args=(x.mean(), x.std(ddof=1))).statistic
    rng = np.random.default_rng(semente)
    i = np.arange(1, n + 1) / n
    Db = np.empty(B)
    for a in range(0, B, 250):  # em pedacos, para nao ocupar memoria
        z = np.sort(rng.standard_normal((min(250, B - a), n)), axis=1)
        z = (z - z.mean(axis=1, keepdims=True)) / z.std(axis=1, ddof=1, keepdims=True)
        F = stats.norm.cdf(z)
        Db[a:a + len(z)] = np.maximum((i - F).max(axis=1), (F - (i - 1. / n)).max(axis=1))
    return D, (1. + (Db >= D).sum()) / (B + 1.), float(np.quantile(Db, 0.95))


def flag(nome, valor, fonte):
    FLAGS[nome] = bool(valor)
    NUMEROS.append(("if" + nome, "verdadeiro" if valor else "falso", None, fonte))


def rel(arq, base):
    return os.path.relpath(arq, base)


def hora(t, arred=0):
    """'21h05' (arred: arredonda para esse numero de minutos); com o dia quando nao e o do relatorio"""
    if arred:
        m = int(round((t.hour * 60 + t.minute + t.second / 60.) / arred) * arred)
        t = dt.datetime.combine(t.date(), dt.time()) + dt.timedelta(minutes=m)
    s = "%dh%02d" % (t.hour, t.minute)
    return s if t.date() == DATA_RELATORIO else s + " de %d/%d" % (t.day, t.month)


def duracao(td):
    m = int(round(td.total_seconds() / 60.))
    return "%d~h %02d~min" % (m // 60, m % 60) if m >= 60 else "%d~min" % m


def norm_cdf(x):
    return 0.5 * math.erfc(-x / math.sqrt(2.))


def norm_pdf(x):
    return math.exp(-0.5 * x * x) / math.sqrt(2. * math.pi)


def wilson(k, n, z=1.959964):
    p = k / n
    c = (p + z * z / (2 * n)) / (1 + z * z / n)
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return c - h, c + h


def ganho_ai(p):
    """estado limite linear com P_f = p: o deslocamento otimo da entropia cruzada, m = phi(beta) / p, na direcao da
    falha, e Var(AI) / Var(direto) por realizacao, com E_h[(I w)^2] = exp(m^2) Phi(-beta - m) (relatorio_chuva.tex,
    Secao 7)"""
    b = brentq(lambda x: norm_cdf(-x) - p, -10., 40.)
    m = norm_pdf(b) / p
    return b, m, (math.exp(m * m) * norm_cdf(-b - m) - p * p) / (p * (1. - p))


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


def le_mc(pasta):
    """as linhas dos mc_<a>_<b>.csv da pasta por realizacao (a -1 e o deterministico), so as completas (a rodada pode
    estar escrevendo) e com situacao ok e Gamma > 0; e as realizacoes sem fator"""
    rs, ruins = {}, set()
    for arq in sorted(glob.glob(os.path.join(pasta, "mc_*_*.csv"))):
        if not re.search(r"/mc_-?\d+_-?\d+\.csv$", arq):
            continue  # mapas
        with open(arq, newline="") as f:
            for r in csv.DictReader(f):
                try:
                    n, x = int(r["real"]), float(r["FS_min"])
                    float(r["log_peso"]), float(r["proj_k"])  # a linha chegou ao fim
                except (TypeError, ValueError, KeyError):
                    continue
                if r["situacao"] != "ok" or not x > 0:
                    ruins.add(n)
                    continue
                rs[n] = r
    return rs, ruins


RX_REAL = re.compile(r"realizacao\s+(-?\d+): Gamma = ([\d.]+) \((\w+)\) \| ([\d.]+) s \(transiente ([\d.]+) s\)")


def le_tempos(pasta):
    """realizacao -> (tempo total, tempo do transiente), s, das linhas dos logs do programa"""
    t = {}
    for arq in glob.glob(os.path.join(pasta, "mc_*_*.log")):
        for m in RX_REAL.finditer(texto(arq)):
            t[int(m[1])] = (float(m[4]), float(m[5]))
    return t


RX_TAREFA = re.compile(r"^(\d{4}-\d\d-\d\d \d\d:\d\d:\d\d) \S+ (\d+)\.\.(\d+): (inicio|fim)", re.M)


def le_tarefas(pasta):
    """os instantes de inicio e de fim das tarefas (logs do cluster/tarefas.sh, um por processo)"""
    ini, fim = [], []
    for arq in glob.glob(os.path.join(pasta, "*.txt")):
        for m in RX_TAREFA.finditer(texto(arq)):
            (ini if m[4] == "inicio" else fim).append(dt.datetime.strptime(m[1], "%Y-%m-%d %H:%M:%S"))
    return ini, fim


def le_caso(arq):
    """a linha do caso de um casos.txt: id, N e as opcoes"""
    for linha in texto(arq).splitlines():
        if linha.strip() and not linha.startswith("#"):
            c = linha.split()
            return c[0], int(c[2]), dict(o.split("=", 1) for o in c[3:] if "=" in o)
    sys.exit("sem caso em " + arq)


def le_mapa(pasta):
    """celulas de 1 m: (x, y) -> [realizacoes com a banda ali, quantas delas falham] e o numero de realizacoes"""
    cel, N = {}, 0
    for arq in glob.glob(os.path.join(pasta, "mc_mapa_*_*.csv")):
        with open(arq, newline="") as f:
            linhas = list(csv.DictReader(f))
        if not linhas:
            continue
        N += int(linhas[0]["n_real"])
        for r in linhas:
            c = cel.setdefault((float(r["x"]), float(r["y"])), [0, 0])
            c[0] += int(r["n"])
            c[1] += int(r["n_falha"])
    return cel, N


def gamma(rs):
    return np.array([float(rs[n]["FS_min"]) for n in sorted(rs)])


def col(rs, c):
    return np.array([float(rs[n][c]) for n in sorted(rs)])


def aleatorias(rs):
    return {n: r for n, r in rs.items() if n >= 0}


# ---------------------------------------------------------------- estimadores
def estima_ai(rs):
    """P_f = soma(peso * falha) / N, o CoV dela, as falhas sorteadas, N, a media dos pesos e o desvio dela, o ESS e o
    peso medio e o maior nas falhas (posproc/amostragem_importancia.py estima)"""
    G = gamma(rs)
    w = np.exp(col(rs, "log_peso"))
    f = G < 1.
    q = w * f
    N = len(q)
    pf = q.mean()
    cov = q.std(ddof=1) / math.sqrt(N) / pf if pf > 0 else math.inf
    return {"N": N, "pf": pf, "cov": cov, "nf": int(f.sum()), "w_media": w.mean(), "w_dp": w.std(ddof=1) / math.sqrt(N),
            "ess": w.sum() ** 2 / (w * w).sum(), "wf_media": w[f].mean() if f.any() else math.nan,
            "wf_max": w[f].max() if f.any() else math.nan, "G": G, "w": w, "q": q}


def corrida(q):
    """P_f e o desvio dela com as n primeiras realizacoes, n = 1..N"""
    n = np.arange(1, len(q) + 1, dtype=float)
    m = np.cumsum(q) / n
    s2 = (np.cumsum(q * q) - n * m * m) / np.maximum(n - 1., 1.)
    return n, m, np.sqrt(np.maximum(s2, 0.) / n)


def compara(ai, p, cov_p):
    """a AI contra o direto (P_f = p, CoV cov_p): diferenca em desvios e quantas vezes menos realizacoes para o mesmo
    CoV (com a P_f da AI no CoV do direto, como o amostragem_importancia.py)"""
    dif = (ai["pf"] - p) / math.sqrt((ai["pf"] * ai["cov"]) ** 2 + (p * cov_p) ** 2)
    cov_d = math.sqrt((1 - ai["pf"]) / (ai["N"] * ai["pf"]))
    return dif, (cov_d / ai["cov"]) ** 2, cov_d


def pareado(a, b):
    """as mesmas realizacoes em dois casos (a: k fixo, b: k aleatorio): P_f, McNemar e a diferenca de Gamma
    (a logica do covk0/compara_pareado.py)"""
    fa, fb = a < 1., b < 1.
    n01, n10 = int((~fa & fb).sum()), int((fa & ~fb).sum())
    N = len(a)
    pm = stats.binomtest(min(n01, n10), n01 + n10, 0.5).pvalue if n01 + n10 else 1.
    d = b - a
    return {"N": N, "pa": fa.mean(), "pb": fb.mean(), "n01": n01, "n10": n10, "p": pm, "dif": (n01 - n10) / N,
            "se": math.sqrt(max(n01 + n10 - (n01 - n10) ** 2 / N, 0.)) / N, "d_media": d.mean(),
            "d_se": d.std(ddof=1) / math.sqrt(N), "d_mediana": float(np.median(d)), "d_mais": (d > 0).mean(), "d": d}


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


def fig_mc(pasta, G, det, lnmu, lnsig, Gant, n, pfn, sen, p_ant, n_ant):
    fig, (a, b, c) = plt.subplots(1, 3, figsize=(LARGURA, 2.35))
    # (a) histograma de Gamma com a densidade lognormal
    bins = np.arange(0.5, 3.1001, 0.05)
    h, e = np.histogram(G, bins=bins, density=True)
    cor = [COR_BANDA if x1 <= 1.0 + 1e-9 else COR_DIRETO for x1 in e[1:]]
    a.bar(e[:-1], h, width=0.05, align="edge", color=cor, edgecolor="white", linewidth=0.4)
    x = np.linspace(0.45, 3.1, 300)
    a.plot(x, stats.lognorm.pdf(x, lnsig, scale=math.exp(lnmu)), color=TINTA, lw=1.1)
    a.axvline(det, color=MUDO, ls="--", lw=0.9)
    a.text(det + 0.05, 0.97, "determinístico", transform=a.get_xaxis_transform(), fontsize=6.5, color=TINTA2,
           va="top")
    a.text(0.97, 0.85, "log-normal", transform=a.transAxes, fontsize=6.5, color=TINTA, ha="right")
    a.annotate("$\\Gamma < 1$", (0.92, 0.45), xytext=(0.78, 1.05), fontsize=6.5, color=TINTA2, ha="center",
               arrowprops={"arrowstyle": "-", "color": MUDO, "lw": 0.6})
    a.set_xlabel("$\\Gamma$")
    a.set_ylabel("densidade")
    a.set_xlim(0.45, 3.1)
    a.set_ylim(0, 1.2 * h.max())
    letra(a, "a")
    # (b) distribuicao acumulada: transiente e a campanha antiga (permanente)
    for g, cor_, rot in ((G, COR_DIRETO, "transiente"), (Gant, COR_MODOS, "permanente (antiga)")):
        s = np.sort(g)
        b.plot(s, np.arange(1, len(s) + 1) / len(s), color=cor_, lw=1.3, label=rot)
    b.plot([1.], [CERON["pf_s1"] / 100.], "*", color=TINTA, ms=7, label="Ceron (Tabela 4)")
    # rotulos diretos: a esquerda da curva do transiente, a direita da antiga e ao lado da estrela
    b.text(1.17, 0.72, "transiente", color=COR_DIRETO, fontsize=6.8, ha="right", va="center")
    b.text(1.62, 0.42, "permanente\n(antiga)", color=COR_MODOS, fontsize=6.8, ha="left", va="center")
    b.text(1.07, CERON["pf_s1"] / 100., "Ceron", color=TINTA, fontsize=6.8, ha="left", va="center")
    b.axvline(1., color=MUDO, lw=0.8)
    b.set_xlim(0.5, 2.5)
    b.set_ylim(0, 1)
    b.set_xlabel("$\\Gamma$")
    b.set_ylabel("$P(\\Gamma \\leq x)$")
    letra(b, "b")
    # (c) convergencia de P_f com N, com a faixa de 95%
    i0 = 49
    c.fill_between(n[i0:], 100 * (pfn - 1.96 * sen)[i0:], 100 * (pfn + 1.96 * sen)[i0:], color=COR_DIRETO, alpha=0.18,
                   lw=0)
    c.plot(n[i0:], 100 * pfn[i0:], color=COR_DIRETO, label="direto (95%)")
    c.axhline(CERON["pf_s1"], color=TINTA, ls="--", lw=0.9)
    c.text(n[-1], CERON["pf_s1"] + 0.4, "Ceron", fontsize=6.5, color=TINTA2, ha="right", va="bottom")
    c.axhline(100 * p_ant, color=COR_MODOS, ls=":", lw=1.1)
    c.text(n[-1], 100 * p_ant - 0.4, "permanente (antiga, %d)" % n_ant, fontsize=6.5, color=TINTA2, ha="right",
           va="top")
    c.set_xlabel("realizações $N$")
    c.set_ylabel("$P_f$ (%)")
    c.set_xlim(0, n[-1])
    c.set_ylim(5, 25)
    letra(c, "c")
    fig.tight_layout()
    salva(fig, pasta, "fig_rb_mc")


def fig_malha(pasta, med, pfs, med_ant, pfs_ant, det_niv, det_ant, ext_med, ext_pf):
    fig, (a, b) = plt.subplots(1, 2, figsize=(LARGURA, 2.3))
    k = np.arange(len(med))
    a.plot(k, med, "o-", color=COR_DIRETO, ms=4, label="transiente: média")
    a.plot(k, det_niv, "s--", color=COR_DIRETO, ms=3.5, mfc="white", lw=1.0, label="transiente: determinístico")
    a.plot(k, med_ant, "o-", color=COR_MODOS, ms=4, label="permanente: média")
    a.plot(k, det_ant, "s--", color=COR_MODOS, ms=3.5, mfc="white", lw=1.0, label="permanente: determinístico")
    a.plot([k[-1], k[-1] + 1], [med[-1], ext_med], ":", color=COR_DIRETO, lw=1.0)
    a.plot([k[-1] + 1], [ext_med], "o", color=COR_DIRETO, mfc="white", ms=4)
    a.axhline(CERON["gamma_det"], color=TINTA, ls="--", lw=0.8)
    a.text(-0.1, CERON["gamma_det"] - 0.012, "Ceron: determinístico", fontsize=6.5, color=TINTA2, va="top")
    a.set_xticks(range(len(med) + 1))
    a.set_xlabel("nível de refinamento")
    a.set_ylabel("$\\Gamma$")
    a.legend(loc="upper right", fontsize=6.5)
    letra(a, "a")
    b.plot(k, 100 * np.asarray(pfs), "o-", color=COR_DIRETO, ms=4, label="transiente")
    b.plot(k, 100 * np.asarray(pfs_ant), "o-", color=COR_MODOS, ms=4, label="permanente (antiga)")
    b.plot([k[-1], k[-1] + 1], [100 * pfs[-1], 100 * ext_pf], ":", color=COR_DIRETO, lw=1.0)
    b.plot([k[-1] + 1], [100 * ext_pf], "o", color=COR_DIRETO, mfc="white", ms=4, label="nível 4 (extrapolado)")
    b.axhline(CERON["pf_s1"], color=TINTA, ls="--", lw=0.8)
    b.text(-0.1, CERON["pf_s1"] + 0.3, "Ceron", fontsize=6.5, color=TINTA2, va="bottom")
    b.set_xticks(range(len(med) + 1))
    b.set_xlabel("nível de refinamento")
    b.set_ylabel("$P_f$ (%)")
    b.legend(loc="upper left", fontsize=6.5)
    letra(b, "b")
    fig.tight_layout()
    salva(fig, pasta, "fig_rb_malha")


def fig_mecanismo(pasta, pts, cel, N, nf, banda_det, rhos):
    fig = plt.figure(figsize=(LARGURA, 2.55))
    ext = fig.add_gridspec(1, 2, width_ratios=[2.05, 1], wspace=0.42, left=0.065, right=0.985, bottom=0.15, top=0.9)
    gm = ext[0, 0].subgridspec(2, 2, height_ratios=[1, 0.055], hspace=0.62, wspace=0.12)
    axs = [fig.add_subplot(gm[0, 0]), fig.add_subplot(gm[0, 1]), fig.add_subplot(ext[0, 1])]
    cax = fig.add_subplot(gm[1, :])
    xs = sorted({x for x, _ in cel})
    ys = sorted({y for _, y in cel})
    x0, x1, y0, y1 = xs[0] - 0.5, xs[-1] + 0.5, ys[0] - 0.5, ys[-1] + 0.5
    X = np.arange(x0, x1 + 0.5)
    Y = np.arange(y0, y1 + 0.5)
    p = np.asarray(pts + [pts[0]], float)
    im = None
    for ax, i, tot, rot, le in ((axs[0], 0, N, "todas (%d)" % N, "a"), (axs[1], 1, nf, "falhas (%d)" % nf, "b")):
        Z = np.full((len(Y) - 1, len(X) - 1), np.nan)
        for (x, y), c in cel.items():
            Z[int(round(y - 0.5 - y0)), int(round(x - 0.5 - x0))] = c[i] / tot if c[i] else np.nan
        ax.fill(p[:, 0], p[:, 1], color=SOLO, lw=0, zorder=0)
        im = ax.pcolormesh(X, Y, Z, cmap="Blues", vmin=0., vmax=1., shading="flat", zorder=1)
        ax.plot(p[:, 0], p[:, 1], color=TINTA2, lw=1.0, zorder=3)
        if i == 0:
            ax.plot(banda_det[:, 0], banda_det[:, 1], ".", color=COR_BANDA, ms=0.9, alpha=0.7, zorder=4)
        ax.set_xlim(x0 - 0.5, x1 + 0.5)
        ax.set_ylim(y0 - 0.5, max(y1, pts[4][1]) + 0.5)
        ax.set_aspect("equal")
        ax.grid(False)
        ax.set_title(rot, fontsize=7.5, color=TINTA)
        ax.set_xlabel("$x$ (m)")
        if i == 0:
            ax.set_ylabel("$y$ (m)")
        else:
            ax.set_yticklabels([])
        letra(ax, le)
    cb = fig.colorbar(im, cax=cax, orientation="horizontal")
    cb.set_label("fração das realizações com a banda na célula de 1 m", fontsize=7)
    cb.ax.tick_params(labelsize=6.5)
    cb.ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda x, p: v(x, 1)))
    c = axs[2]
    nomes = [r for r, _ in rhos][::-1]
    vals = [x for _, x in rhos][::-1]
    c.barh(range(len(vals)), vals, color=COR_DIRETO, height=0.6)
    for i, x in enumerate(vals):
        c.text(x + (0.03 if x >= 0 else -0.03), i, v(x, 2), va="center", ha="left" if x >= 0 else "right",
               fontsize=6.5, color=TINTA)
    c.set_yticks(range(len(vals)))
    c.set_yticklabels(nomes, fontsize=7)
    c.axvline(0, color=TINTA2, lw=0.8)
    c.set_xlim(-0.8, 1.2)
    c.set_xlabel("Spearman com $\\Gamma$")
    c.grid(axis="y", visible=False)
    letra(c, "c")
    salva(fig, pasta, "fig_rb_mecanismo")


def fig_ai(pasta, q_dir, ais, cov_lin):
    fig, (a, b, c) = plt.subplots(1, 3, figsize=(LARGURA, 2.4))
    i0 = 49
    series = [("direto", COR_DIRETO, q_dir, None)] + [(rot, cor, d["q"], d) for rot, cor, d in ais]
    for rot, cor, q, d in series:
        n, m, se = corrida(q)
        a.fill_between(n[i0:], 100 * (m - 1.96 * se)[i0:], 100 * (m + 1.96 * se)[i0:], color=cor, alpha=0.15, lw=0)
        a.plot(n[i0:], 100 * m[i0:], color=cor, lw=1.2, label=rot)
        ok = m > 0
        b.loglog(n[i0:][ok[i0:]], 100 * (se / np.where(m > 0, m, np.nan))[i0:][ok[i0:]], color=cor, lw=1.2, label=rot)
    a.axhline(CERON["pf_s1"], color=TINTA, ls="--", lw=0.8)
    a.text(len(q_dir), CERON["pf_s1"] - 0.4, "Ceron", fontsize=6.5, color=TINTA2, ha="right", va="top")
    a.set_xlabel("realizações $N$")
    a.set_ylabel("$P_f$ (%) e faixa de 95%")
    a.set_ylim(5, 36)
    a.set_xlim(0, len(q_dir))
    a.legend(loc="upper right", fontsize=6.3, bbox_to_anchor=(1.0, 1.02))
    letra(a, "a")
    nn = np.logspace(math.log10(50), math.log10(len(q_dir)), 50)
    b.loglog(nn, 100 * cov_lin / np.sqrt(nn), color=TINTA2, ls="--", lw=0.9, label="AI linear ótima")
    b.set_xlabel("realizações $N$")
    b.set_ylabel("$CoV(P_f)$ (%)")
    b.set_yticks([2, 3, 5, 10, 20, 30])
    b.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    eixo_log_virgula(b, "x")
    eixo_log_virgula(b, "y")
    b.legend(loc="lower left", fontsize=6.3)
    letra(b, "b")
    for rot, cor, d in ais:
        f = d["G"] < 1.
        c.scatter(d["G"][~f], d["w"][~f], s=3, color=MUDO, alpha=0.5, lw=0)
        c.scatter(d["G"][f], d["w"][f], s=4, color=cor, lw=0, alpha=0.8)
    c.axvline(1., color=TINTA2, lw=0.8)
    c.axhline(1., color=MUDO, lw=0.6, ls=":")
    c.set_yscale("log")
    c.set_xlabel("$\\Gamma$")
    c.set_ylabel("peso $w$")
    c.legend(handles=[Line2D([], [], marker="o", ls="", color=cor, ms=3.5, label="falha, " + rot)
                      for rot, cor, _ in ais] + [Line2D([], [], marker="o", ls="", color=MUDO, ms=3.5,
                                                        label="$\\Gamma \\geq 1$")],
             loc="upper left", fontsize=6.3, handletextpad=0.1)
    letra(c, "c")
    fig.tight_layout()
    salva(fig, pasta, "fig_rb_ai")


def fig_kfixo(pasta, a, b, par):
    fig, (x, y) = plt.subplots(1, 2, figsize=(LARGURA, 2.5))
    fa, fb = a < 1., b < 1.
    dis = fa != fb
    x.scatter(a[~dis], b[~dis], s=3, color=MUDO, alpha=0.5, lw=0, label="mesma situação")
    x.scatter(a[~fa & fb], b[~fa & fb], s=5, color=COR_BANDA, lw=0, label="só k aleatório falha (%d)" % par["n01"])
    x.scatter(a[fa & ~fb], b[fa & ~fb], s=5, color=COR_DIRETO, lw=0, label="só k fixo falha (%d)" % par["n10"])
    lim = [0.5, max(a.max(), b.max()) * 1.02]
    x.plot(lim, lim, color=TINTA2, lw=0.8)
    x.axvline(1., color=MUDO, lw=0.6)
    x.axhline(1., color=MUDO, lw=0.6)
    x.set_xlim(lim)
    x.set_ylim(lim)
    x.set_xlabel("$\\Gamma$ com $k$ fixo")
    x.set_ylabel("$\\Gamma$ com $k$ aleatório")
    x.legend(loc="upper left", fontsize=6.3)
    letra(x, "a")
    y.hist(par["d"], bins=60, color=COR_DIRETO, edgecolor="white", linewidth=0.3)
    y.axvline(0., color=TINTA2, lw=0.8)
    y.axvline(par["d_media"], color=COR_BANDA, lw=1.0, ls="--")
    y.set_xlabel("$\\Gamma$(k aleatório) $-$ $\\Gamma$(k fixo)")
    y.set_ylabel("realizações")
    letra(y, "b")
    fig.tight_layout()
    salva(fig, pasta, "fig_rb_kfixo")


# ---------------------------------------------------------------- secoes
def modelo(R, D, DET):
    """geometria, material, campos, malhas, o deterministico e o custo"""
    src = os.path.join(PROJ, "rebaixamento_mc.cpp")
    m = le_fonte(src, r"kTalude = \{(.*?)\};", "kTalude")
    pts = [[float(a), float(b)] for a, b in re.findall(r"\{([-\d.]+),\s*([-\d.]+)\}", m[1])]
    (xb0, yb0), (xb1, _), (_, yd), (xpe, ype), (xcr, ycr), _ = pts
    H = ycr - ype
    num("DomX", "%g" % (xb1 - xb0), "rebaixamento_mc.cpp (kTalude)")
    num("DomY", "%g" % (ycr - yb0), "rebaixamento_mc.cpp (kTalude)")
    num("AlturaH", "%g" % H, "rebaixamento_mc.cpp (kTalude)")
    num("Inclinacao", "%g" % round(math.degrees(math.atan2(H, xpe - xcr))), "rebaixamento_mc.cpp (kTalude)")
    num("CristaX", "%g" % xcr, "rebaixamento_mc.cpp (kTalude)")
    num("PeX", "%g" % xpe, "rebaixamento_mc.cpp (kTalude)")
    num("CristaY", "%g" % ycr, "rebaixamento_mc.cpp (kTalude)")
    num("PeY", "%g" % ype, "rebaixamento_mc.cpp (kTalude)")
    solo = dict(re.findall(r"(\w+) = ([\d.e+-]+)", le_fonte(src, r"struct Solo \{(.*?)\};", "Solo")[1]))
    num("Coesao", "%g" % float(solo["c"]), "rebaixamento_mc.cpp (Solo)")
    num("Atrito", "%g" % float(solo["phiDeg"]), "rebaixamento_mc.cpp (Solo)")
    num("PesoSolo", "%g" % float(solo["gamma"]), "rebaixamento_mc.cpp (Solo)")
    num("ModuloE", "%g" % (float(solo["E"]) / 1000.), "rebaixamento_mc.cpp (Solo, MPa)")
    num("Poisson", v(float(solo["nu"]), 1), "rebaixamento_mc.cpp (Solo)")
    num("Porosidade", v(float(solo["n"]), 2), "rebaixamento_mc.cpp (Solo)")
    kv = float(solo["kv"])
    num("Kv", "10^{%d}" % round(math.log10(kv)), "rebaixamento_mc.cpp (Solo)", latex="10^{%d}" % round(math.log10(kv)))
    num("ModuloAgua", v(float(solo["Kw"]) / 1e6, 1), "rebaixamento_mc.cpp (Solo, GPa)")
    cm = dict(re.findall(r"(\w+) = ([\d.e+-]+)", le_fonte(src, r"struct Campos \{.*?\n(.*?)\};", "Campos")[1]))
    m3 = le_fonte(src, r"Gamma = (\d,\d+) no talude ceron2 com (\d) niveis", "Gamma do main.cpp")
    gmain = float(m3[1].replace(",", "."))
    num("GammaMain", v(gmain), "rebaixamento_mc.cpp (cabecalho: TransientDrawDown do main.cpp)")
    num("NivelMain", m3[2], "rebaixamento_mc.cpp (cabecalho)")
    par = le_chaves(os.path.join(D, "mc_parametros.txt"))
    fpar = rel(os.path.join(D, "mc_parametros.txt"), R)
    for chave, macro, casas in (("COV_C", "CovC", 1), ("COV_PHI", "CovPhi", 1), ("COV_K", "CovK", 1)):
        num(macro, v(float(par[chave]), casas), fpar)
    assert abs(float(par["COV_C"]) - float(cm["covC"])) < 1e-12
    num("Lx", "%g" % float(par["LX"]), fpar)
    num("Ly", "%g" % float(par["LY"]), fpar)
    num("KLM", par["KL_M"], fpar)
    num("NVarKL", "%d" % (3 * int(par["KL_M"])), fpar + " (3 campos x KL_M)")
    num("KLvar", v(100. * float(par["KL_VARIANCIA"]), 2), fpar)
    num("Semente", par["SEMENTE"], fpar)
    num("Trebaix", v(float(par["TREBAIX"]), 1), fpar)
    num("NivelFinal", "%g" % float(par["NIVEL_FINAL"]), fpar)
    num("HwH", "%g" % float(par["HW_H"]), fpar)
    num("Nref", par["NREF"], fpar)
    txt("MalhaT", par["MALHA_TRANSIENTE"], fpar)
    txt("MalhaF", par["MALHA_FATOR"], fpar)
    lim = float(le_fonte(src, r"ind\[gel->Index\(\)\] >= ([\d.]+) \* vmax", "limiar do refinamento")[1])
    num("LimiarBanda", v(lim, 1), "rebaixamento_mc.cpp (FatorAdaptativo)")
    dth = float(le_fonte(os.path.join(PROJ, "DrawdownTransient.h"), r"REAL dt = ([\d.]+);", "dt")[1])
    num("DtTrans", v(dth, 1), "DrawdownTransient.h (DrawdownData::dt)")
    num("PassosTrans", "%d" % round(float(par["TREBAIX"]) / dth), "TREBAIX / dt")
    # malha do transiente: os triangulos e o lado (a menor aresta) do deterministico
    dm = os.path.join(DET, "mc_real-1")
    V = np.loadtxt(os.path.join(dm, "vertices.csv"), delimiter=",", skiprows=1)
    T = np.loadtxt(os.path.join(dm, "triangles.csv"), delimiter=",", skiprows=1, dtype=int)
    P = V[:, 1:3]
    lados = np.concatenate([np.hypot(*(P[T[:, i]] - P[T[:, j]]).T) for i, j in ((0, 1), (1, 2), (2, 0))])
    lado = lados.min()
    kT = int(par["MALHA_TRANSIENTE"].replace("tri", ""))
    kF = int(par["MALHA_FATOR"].replace("tri", ""))
    passo = lado * 2 ** kT
    ftri = rel(os.path.join(dm, "triangles.csv"), R)
    num("TriTrans", "%d" % len(T), ftri)
    num("LadoTrans", v(lado, 3), ftri + " (menor aresta)")
    num("PassoMalha", "%g" % passo, ftri + " (lado x 2^%d)" % kT)
    num("LadoFator", v(passo / 2 ** kF, 1), "passo / 2^%d" % kF)
    num("LadoBanda", v(passo / 2 ** kF / 2 ** int(par["NREF"]), 2), "passo / 2^%d / 2^NREF" % kF)
    return pts, H, gmain, par


def direto(R, D, DET, pts):
    """o Monte Carlo direto e o deterministico"""
    fonte = rel(D, R) + "/mc_*.csv"
    rs, ruins = le_mc(D)
    assert not ruins, "realizacoes sem fator no direto: %s" % sorted(ruins)[:5]
    detr = le_mc(DET)[0][-1]
    assert detr["FS_min"] == rs[-1]["FS_min"], "o deterministico do direto difere do de det/"
    fdet = rel(DET, R) + "/mc_0_-1.csv"
    gdet = float(detr["FS_min"])
    niv_det = [float(x) for x in detr["niveis"].split(";")]
    num("GammaDet", v(gdet), fdet)
    num("NiveisDet", "; ".join(v(x, 4) for x in niv_det), fdet)
    num("EqFatorDet", detr["eq_fator"], fdet)
    num("DifDetCeron", v(100. * (gdet / CERON["gamma_det"] - 1.), 2), fdet + " e Ceron (Fig. 24)")
    num("GammaCeron", v(CERON["gamma_det"], 3), "Ceron, Fig. 24 (digitado)")
    num("PfCeron", v(CERON["pf_s1"], 2), "Ceron, Tabela 4 (digitado)")
    num("PfCeronS", v(CERON["pf_s400"], 2), "Ceron, Tabela 4 (digitado)")
    num("CovCeron", v(CERON["cov_s1"], 1), "Ceron, Tabela 4 (digitado)")
    # a referencia do Ceron: S = 10 000 (Secao 6.1), media e desvio de Gamma (Tabela 4) e o caso com k_v fixo (Tabela 3)
    nc, pc = CERON["n_s1"], CERON["pf_s1"] / 100.
    assert abs(100. * math.sqrt((1. - pc) / (nc * pc)) - CERON["covpf_s1"]) < 0.01  # CoV(P_f) da Tabela 4 com S
    num("NCeron", mil(nc), "Ceron, Secao 6.1 (S = 10 000 simulacoes; digitado)")
    num("CovPfCeron", v(CERON["covpf_s1"], 2), "Ceron, Secao 6.1 e Tabela 4 (digitado)")
    lo_c, hi_c = wilson(round(pc * nc), nc)
    num("PfCeronLo", v(100. * lo_c, 1), "Ceron, P_f e S (Wilson)")
    num("PfCeronHi", v(100. * hi_c, 1), "Ceron, P_f e S (Wilson)")
    num("GammaMediaCeron", v(CERON["media_s1"], 3), "Ceron, Tabela 4 (digitado)")
    num("GammaDPCeron", v(CERON["dp_s1"], 3), "Ceron, Tabela 4 (digitado)")
    sob_c = CERON["media_s1"] / CERON["gamma_det"] - 1.
    num("MediaSobreDetCeron", v(100. * sob_c, 1), "Ceron, Tabela 4 e Fig. 24")
    txt("PosicaoCeron", acima_abaixo(100. * sob_c), "Ceron, Tabela 4 e Fig. 24")
    num("PfCeronKfixo", v(CERON["pf_k0"], 2), "Ceron, Tabela 3, CoV(k_v) = 0 (digitado)")
    num("GammaMediaCeronKfixo", v(CERON["media_k0"], 3), "Ceron, Tabela 3, CoV(k_v) = 0 (digitado)")
    num("CeronDifK", v(CERON["pf_s1"] - CERON["pf_k0"], 2), "Ceron, Tabelas 3 e 4 (pontos)")
    al = aleatorias(rs)
    N = len(al)
    G = gamma(al)
    f = G < 1.
    nf = int(f.sum())
    p = nf / N
    lo, hi = wilson(nf, N)
    num("NDireto", "%d" % N, fonte)
    num("NfDireto", "%d" % nf, fonte)
    num("PfDireto", v(100. * p, 2), fonte)
    num("PfDiretoLo", v(100. * lo, 1), fonte + " (Wilson)")
    num("PfDiretoHi", v(100. * hi, 1), fonte + " (Wilson)")
    cov_p = math.sqrt((1 - p) / (N * p))
    num("CovPfDireto", v(100. * cov_p, 1), fonte + " (Eq. 68 do Ceron)")
    num("GammaMedia", v(G.mean()), fonte)
    num("GammaMediana", v(float(np.median(G))), fonte)
    num("GammaDP", v(G.std(ddof=1)), fonte)
    num("GammaCoV", v(100. * G.std(ddof=1) / G.mean(), 1), fonte)
    num("GammaMin", v(G.min()), fonte)
    num("GammaMax", v(G.max()), fonte)
    num("MediaSobreDet", v(100. * (G.mean() / gdet - 1.), 1), fonte)
    txt("PosicaoDireto", acima_abaixo(100. * (G.mean() / gdet - 1.)), fonte + " e " + fdet)
    txt("MediaContraCeron", acima_abaixo(100. * (G.mean() / CERON["media_s1"] - 1.)), fonte + " e Ceron (Tabela 4)")
    ln = np.log(G)
    lnmu, lnsig = ln.mean(), ln.std(ddof=1)
    num("LnMu", v(lnmu), fonte)
    num("LnSigma", v(lnsig), fonte)
    num("PfLogn", v(100. * norm_cdf(-lnmu / lnsig), 2), fonte + " (ajuste log-normal)")
    # a log-normal com a media do Ceron e o CoV daqui: a posicao explica a diferenca de P_f?
    s2c = math.log(1. + (G.std(ddof=1) / G.mean()) ** 2)
    num("PfLognCeron", v(100. * norm_cdf(-(math.log(CERON["media_s1"]) - 0.5 * s2c) / math.sqrt(s2c)), 1),
        fonte + " (CoV) e Ceron (media, Tabela 4): log-normal")
    # Kolmogorov-Smirnov com a media e o desvio da propria amostra: o p e o de Lilliefors (bootstrap parametrico)
    B_l = 2000
    ksD, ksP, ksC = lilliefors(ln, B_l)
    knD, knP, _ = lilliefors(G, B_l, semente=2)
    num("LillieB", "%d" % B_l, "bootstrap parametrico (script)")
    num("KSlogD", v(ksD, 4), fonte + " (Kolmogorov-Smirnov, log-normal)")
    t_, l_ = p_rel(ksP, 2, aprox=True)
    num("KSlogP", t_, fonte + " (Lilliefors, %d amostras)" % B_l, latex=l_)
    num("KSCritico", v(ksC, 4), "Lilliefors: D de 5%% (%d amostras de n = %d)" % (B_l, len(ln)))
    num("KSnormD", v(knD, 4), fonte + " (Kolmogorov-Smirnov, normal)")
    t_, l_ = p_rel(knP, 2, aprox=True)
    num("KSnormP", t_, fonte + " (Lilliefors, %d amostras)" % B_l, latex=l_)
    n, pfn, sen = corrida(f.astype(float))
    for k, nome in ((400, "Quatrocentas"), (1000, "Mil")):  # 400: o piloto do relatorio_chuva.tex
        num("NPrimeiras" + nome, "%d" % k, "janela escolhida no script")
        num("PfPrimeiras" + nome, v(100. * pfn[k - 1], 1), fonte + " (realizacoes 0..%d)" % (k - 1))
    # o que leva a ruptura: Spearman simples e parcial de 2a ordem (Eqs. 70-71 do Ceron) e as medias
    nomes = [("banda_c", "C", "$c'$ na banda", 2), ("banda_phi_deg", "Phi", "$\\varphi'$ na banda", 2),
             ("banda_k_mult", "K", "$k_v$ na banda", 3), ("banda_pw", "Pw", "$p_w$ na banda", 2),
             ("banda_prof", "Prof", "profundidade", 2)]
    rhos = []
    for c, mac, rot, casas in nomes:
        x = col(al, c)
        rho = stats.spearmanr(x, G)[0]
        rhos.append((rot, rho))
        num("Rho" + mac, v(rho, 2), fonte + " (Spearman de Gamma com %s)" % c)
        num("Falha" + mac, v(x[f].mean(), casas), fonte + " (media de %s nas falhas)" % c)
        num("Outras" + mac, v(x[~f].mean(), casas), fonte + " (media de %s nas outras)" % c)
        num("Det" + mac, v(float(detr[c]), casas), fdet + " (%s)" % c)
    X = np.column_stack([G] + [col(al, c) for c in ("banda_c", "banda_phi_deg", "banda_k_mult")])
    Pm = np.linalg.inv(stats.spearmanr(X)[0])
    for i, mac in ((1, "C"), (2, "Phi"), (3, "K")):
        num("Parc" + mac, v(-Pm[0, i] / math.sqrt(Pm[0, 0] * Pm[i, i]), 2),
            fonte + " (Spearman parcial de 2a ordem, Eqs. 70-71 do Ceron)")
    num("RhoKPw", v(stats.spearmanr(col(al, "banda_k_mult"), col(al, "banda_pw"))[0], 2),
        fonte + " (Spearman de banda_k_mult com banda_pw)")
    num("RhoPwProf", v(stats.spearmanr(col(al, "banda_pw"), col(al, "banda_prof"))[0], 2),
        fonte + " (Spearman de banda_pw com banda_prof)")
    # niveis de refinamento
    L = np.array([[float(x) for x in al[n_]["niveis"].split(";")] for n_ in sorted(al)])
    med = L.mean(axis=0)
    pfs = [(L[:, k] < 1.).mean() for k in range(L.shape[1])]
    nomes_niv = ["Zero", "Um", "Dois", "Tres"]
    for k in range(L.shape[1]):
        num("MediaNivel" + nomes_niv[k], v(med[k]), fonte + " (niveis)")
        num("PfNivel" + nomes_niv[k], v(100. * pfs[k], 2), fonte + " (niveis)")
        num("DetNivel" + nomes_niv[k], v(niv_det[k]), fdet + " (niveis)")
    ult = L[:, -1] / L[:, -2] - 1.
    num("UltimoNivel", v(100. * ult.mean(), 1), fonte + " (Gamma_3 / Gamma_2 - 1, media)")
    num("UltimoNivelDet", v(100. * (niv_det[-1] / niv_det[-2] - 1.), 1), fdet)
    r = (med[-1] - med[-2]) / (med[-2] - med[-3])
    num("RazaoNiveis", v(r, 2), fonte + " (razao das duas ultimas variacoes da media)")
    ext_med = med[-1] + r * (med[-1] - med[-2])
    ext_lim = med[-1] + (med[-1] - med[-2]) * r / (1. - r)
    d4, dlim = ext_med / med[-1] - 1., ext_lim / med[-1] - 1.
    num("ExtrapQuatro", v(100. * d4, 1), "extrapolacao geometrica das medias por nivel")
    num("ExtrapLimite", v(100. * dlim, 1), "extrapolacao geometrica das medias por nivel")
    pf4, pflim = (G * (1. + d4) < 1.).mean(), (G * (1. + dlim) < 1.).mean()
    num("PfExtrapQuatro", v(100. * pf4, 1), fonte + " (Gamma x (1 + extrapolacao) < 1)")
    num("PfExtrapLimite", v(100. * pflim, 1), fonte + " (Gamma x (1 + extrapolacao) < 1)")
    for dl, mac in ((0.01, "PfMenosUm"), (0.02, "PfMenosDois")):
        num(mac, v(100. * (G * (1. - dl) < 1.).mean(), 1), fonte + " (Gamma x %g < 1)" % (1. - dl))
    eq = col(al, "eq_fator")
    num("EqFatorMin", "%d" % eq.min(), fonte)
    num("EqFatorMax", "%d" % eq.max(), fonte)
    num("EqTrans", al[0]["eq_transiente"], fonte)
    # custo: os tempos de cada realizacao nos logs e o relogio da rodada nos logs das tarefas
    t = le_tempos(D)
    tot = np.array([t[n_][0] for n_ in sorted(al)])
    tra = np.array([t[n_][1] for n_ in sorted(al)])
    num("CustoReal", v(tot.mean(), 1), rel(D, R) + "/mc_*.log (tempo por realizacao)")
    num("CustoMediana", v(float(np.median(tot)), 1), rel(D, R) + "/mc_*.log")
    num("CustoTrans", v(tra.mean(), 1), rel(D, R) + "/mc_*.log (transiente)")
    custo_fator = col(al, "s_fator").mean()
    num("CustoFator", v(custo_fator, 1), fonte + " (s_fator)")
    num("CustoResto", v((tot - tra - col(al, "s_fator")).mean(), 1),
        rel(D, R) + "/mc_*.log e mc_*.csv (total - transiente - s_fator)")
    num("CustoPrimeiras", v(tot[:400].mean(), 1), rel(D, R) + "/mc_*.log (realizacoes 0..399)")
    num("CustoUltimas", v(tot[-400:].mean(), 1), rel(D, R) + "/mc_*.log (realizacoes 1600..1999)")
    a0, a1 = map(int, le_fonte(os.path.join(R, "direto", "rodar.sh"), r"seq (\d+) (\d+)", "seq").groups())
    nproc = a1 - a0 + 1
    num("ProcPC", "%d" % nproc, "direto/rodar.sh")
    ini, fim = le_tarefas(os.path.join(R, "direto", "logs"))
    T2000 = max(fim) - min(ini)
    txt("RelogioDireto", duracao(T2000), "direto/logs/pp*.txt")
    txt("InicioDireto", hora(min(ini)), "direto/logs/pp*.txt")
    txt("FimDireto", hora(max(fim)), "direto/logs/pp*.txt")
    num("RelogioPorReal", v(T2000.total_seconds() / N, 1), "direto/logs/pp*.txt / N")
    # mapa da banda
    cel, Nm = le_mapa(D)
    assert Nm == N, "o mapa tem %d realizacoes" % Nm
    mec = np.loadtxt(os.path.join(DET, "mc_real-1", "mecanismo.csv"), delimiter=",", skiprows=1)
    banda_det = mec[mec[:, 3] >= 0.1 * mec[:, 3].max()][:, 1:3]
    num("CelulasMapa", "%d" % len(cel), rel(D, R) + "/mc_mapa_*.csv")
    um = sum(1 for c in cel.values() if c[0] >= 0.5 * N)
    umf = sum(1 for c in cel.values() if c[1] >= 0.5 * nf)
    num("CelulasMetade", "%d" % um, rel(D, R) + "/mc_mapa_*.csv (celulas com a banda em >= 50% das realizacoes)")
    num("CelulasMetadeFalha", "%d" % umf, rel(D, R) + "/mc_mapa_*.csv (o mesmo nas falhas)")
    return {"rs": rs, "al": al, "G": G, "p": p, "cov_p": cov_p, "N": N, "nf": nf, "gdet": gdet, "niv_det": niv_det,
            "lnmu": lnmu, "lnsig": lnsig, "n": n, "pfn": pfn, "sen": sen, "med": med, "pfs": pfs,
            "ext_med": ext_med, "pf4": pf4, "custo": tot.mean(), "T2000": T2000, "cel": cel,
            "banda_det": banda_det, "rhos": rhos, "nproc": nproc, "f": f, "custo_fator": custo_fator}


def antiga(ANT, md):
    """a campanha antiga do rebaixamento (percolacao permanente, cluster): distribuicao, niveis, o pareamento e o
    teste do k fixo dela"""
    pasta = os.path.join(ANT, "ceron_mc_fase1", "f50_I0_iii")
    arq = os.path.join(pasta, "mc_resumo.csv")
    fonte = rel(arq, os.path.dirname(ANT))
    with open(arq, newline="") as f:
        old = {int(r["real"]): r for r in csv.DictReader(f) if r["situacao"] == "ok"}
    par = le_chaves(os.path.join(pasta, "parametros.txt"))
    fpar = rel(os.path.join(pasta, "parametros.txt"), os.path.dirname(ANT))
    txt("AguaAntiga", par["AGUA"].replace("_", r"\_"), fpar)
    txt("MalhaAntiga", par["MESH"].replace("_", r"\_"), fpar)
    num("NrefAntiga", par["NREF"], fpar)
    num("PmaxAntiga", par["PMAX"], fpar)
    G = np.array([float(old[n]["FS_min"]) for n in sorted(old)])
    N, nf = len(G), int((G < 1).sum())
    lo, hi = wilson(nf, N)
    num("NAntiga", "%d" % N, fonte)
    num("NfAntiga", "%d" % nf, fonte)
    num("PfAntiga", v(100. * nf / N, 2), fonte)
    num("PfAntigaLo", v(100. * lo, 1), fonte + " (Wilson)")
    num("PfAntigaHi", v(100. * hi, 1), fonte + " (Wilson)")
    num("GammaMediaAntiga", v(G.mean()), fonte)
    num("GammaCoVAntiga", v(100. * G.std(ddof=1) / G.mean(), 1), fonte)
    num("GammaMedianaAntiga", v(float(np.median(G))), fonte)
    fdet = os.path.join(ANT, "ceron_mc_fase1", "f50_I0_iii_det", "r00000", "gim_unsat.csv")
    with open(fdet, newline="") as f:
        dl = list(csv.DictReader(f))
    det_ant = [float(r["gamma"]) for r in dl]
    num("GammaDetAntiga", v(det_ant[-1]), rel(fdet, os.path.dirname(ANT)))
    num("MediaSobreDetAntiga", v(100. * (G.mean() / det_ant[-1] - 1.), 1), fonte)
    # o retangulo da KL: a caixa dos nos da malha (Common/SlopeFields.cpp, Box)
    msh = os.path.join(NEOPZ, "Projects", "LearningMohrCoulomb2", "cluster", "malhas", par["MESH"])
    Lm = texto(msh).split("\n")
    i = Lm.index("$Nodes")
    nb = int(Lm[i + 1].split()[0])
    j, xs, ys = i + 2, [], []
    for _ in range(nb):
        k = int(Lm[j].split()[3])
        j += 1 + k
        for _ in range(k):
            x, y = map(float, Lm[j].split()[:2])
            xs.append(x), ys.append(y)
            j += 1
    num("DomXAntiga", "%g" % (max(xs) - min(xs)), rel(msh, NEOPZ) + " (caixa dos nos)")
    num("DomYAntiga", "%g" % (max(ys) - min(ys)), rel(msh, NEOPZ) + " (caixa dos nos)")
    num("KLvarAntiga", v(100. * float(old[0]["KL_variancia"]), 2), fonte + " (KL_variancia)")
    num("AreaAntiga", v(float(old[0]["area_solo_m2"]), 1), fonte + " (area_solo_m2)")
    # a realizacao n tem os mesmos campos nas duas? correlacao das medias para o mesmo n
    al = md["al"]
    com = sorted(n for n in al if n in old)
    num("NComuns", "%d" % len(com), fonte + " e direto")
    for c, mac in (("c_med", "CorrCmed"), ("phi_med", "CorrPhimed"), ("k_mult_med", "CorrKmed"),
                   ("banda_c", "CorrBandaC"), ("FS_min", "CorrGamma")):
        a = np.array([float(al[n][c]) for n in com])
        b = np.array([float(old[n][c]) for n in com])
        num(mac, v(float(np.corrcoef(a, b)[0, 1]), 2), fonte + " e direto (Pearson de %s, mesmo n)" % c)
    Ga = np.array([float(old[n]["FS_min"]) for n in com])
    num("PfAntigaComuns", v(100. * (Ga < 1).mean(), 2), fonte + " (as realizacoes 0..%d)" % (len(com) - 1))
    # niveis
    Lv = np.array([[float(old[n]["gamma_nivel%d" % k]) for k in range(4)] for n in sorted(old)])
    pfs = [(Lv[:, k] < 1.).mean() for k in range(4)]
    for k, nome in enumerate(["Zero", "Um", "Dois", "Tres"]):
        num("PfAntigaNivel" + nome, v(100. * pfs[k], 2), fonte + " (gamma_nivel%d)" % k)
    # o teste do k fixo na campanha antiga (as mesmas realizacoes, mesma malha e mesmo retangulo)
    a0 = os.path.join(ANT, "ceron_mc_fase2", "f50_I0_iii_covk0", "mc_resumo.csv")
    with open(a0, newline="") as f:
        k0 = {int(r["real"]): r for r in csv.DictReader(f) if r["situacao"] == "ok"}
    cm = sorted(set(k0) & set(old))
    par0 = pareado(np.array([float(k0[n]["FS_min"]) for n in cm]), np.array([float(old[n]["FS_min"]) for n in cm]))
    f0 = rel(a0, os.path.dirname(ANT)) + " e " + fonte
    num("AntNPar", "%d" % par0["N"], f0)
    num("AntPfFixo", v(100. * par0["pa"], 2), f0)
    num("AntPfAleat", v(100. * par0["pb"], 2), f0)
    num("AntNzu", "%d" % par0["n01"], f0 + " (so com k aleatorio)")
    num("AntNuz", "%d" % par0["n10"], f0 + " (so com k fixo)")
    t_, l_ = p_rel(par0["p"], 2)
    num("AntMcNemar", t_, f0 + " (McNemar)", latex=l_)
    num("AntDifGamma", v(par0["d_media"], 4), f0 + " (Gamma(aleatorio) - Gamma(fixo))")
    num("AntDifGammaSe", v(par0["d_se"], 4), f0)
    # o retangulo da KL do transiente (RandomField no rebaixamento_mc.cpp: 0, 0, p[1][0], p[4][1])
    return {"G": G, "N": N, "p": nf / N, "pfs": pfs, "med": Lv.mean(axis=0), "det": det_ant}


def le_desloc_banda(R):
    d = texto(os.path.join(R, "ai", "desloc.txt"))
    m = re.search(r"piloto: (\d+) realizacoes com projecao, (\d+) falhas", d)
    s = re.search(r"s_c = ([-\d.]+), s_phi = ([-\d.]+), s_k = ([-\d.]+)", d)
    mu2 = re.search(r"\|mu\|\^2 = ([\d.]+)", d)
    return int(m[1]), int(m[2]), [float(x) for x in s.groups()], float(mu2[1])


def eta(inicio, n_alvo, n_ref, T_ref):
    return inicio + dt.timedelta(seconds=T_ref.total_seconds() * n_alvo / n_ref)


def ai(R, md):
    """as duas amostragens por importancia, a conferencia dos pesos e o ganho linear"""
    p, cov_p, N = md["p"], md["cov_p"], md["N"]
    al = md["al"]
    f = md["f"]
    fonte_d = rel(os.path.join(R, "direto", "ceron_tr_fase1", "f50_I0_iii"), R) + "/mc_*.csv"
    # ganho esperado num estado limite linear com a P_f do direto (e com 1%, a da chuva)
    b, m, r = ganho_ai(p)
    num("BetaLinear", v(b, 2), "formula (relatorio_chuva.tex, Secao 7), P_f do direto")
    num("MdoisLinear", v(m * m, 2), "formula, P_f do direto")
    num("GanhoLinear", v(1. / r, 1), "formula, P_f do direto")
    pfc = float(le_fonte(os.path.join(PROJ, "cluster", "README.md"), r"tem P_f ≈ (\d+)% em 24 h, pelos pilotos",
                         "P_f da chuva")[1]) / 100.
    num("PfChuva", "%g" % (100. * pfc), "cluster/README.md (pilotos do relatorio_chuva.tex)")
    b1, m1, r1 = ganho_ai(pfc)
    num("GanhoLinearUm", v(1. / r1, 1), "formula, P_f da chuva")
    # o deslocamento na direcao da banda (ai/desloc.txt) e a conferencia com as colunas proj_* do direto
    npil, nfp, s, mu2 = le_desloc_banda(R)
    proj = np.column_stack([col(al, c) for c in ("proj_c", "proj_phi", "proj_k")])
    s_cal = proj[f].mean(axis=0)
    assert npil == N and nfp == md["nf"] and np.allclose(s_cal, s, atol=1e-3), (s, s_cal)
    for mac, x in zip(("SC", "SPhi", "SK"), s):
        num(mac, v(x, 2), "ai/desloc.txt")
    num("MuDoisBanda", v(mu2, 2), "ai/desloc.txt")
    num("VarPesoBanda", v(math.exp(mu2) - 1., 2), "exp(|mu|^2) - 1")
    num("RuidoBanda", v(3. / nfp, 3), "3 / n_f")
    # rodada da AI na banda
    A = os.path.join(R, "ai", "ceron_tr_fase1_ai", "f50_I0_iii")
    idb, nalvo, opb = le_caso(os.path.join(R, "ai", "casos.txt"))
    num("NAlvoBanda", "%d" % nalvo, "ai/casos.txt")
    num("SementeBanda", opb["semente"], "ai/casos.txt")
    rsb, ruinsb = le_mc(A)
    detb = rsb.pop(-1, None)
    if detb is not None:
        num("GammaDetBanda", v(float(detb["FS_min"])), rel(A, R) + "/mc_0_49.csv")
    bpronto = len(rsb) + len(ruinsb) >= nalvo and os.path.exists(os.path.join(R, "ai", "resultado.txt"))
    flag("aipronto", bpronto, "ai/resultado.txt e %d realizacoes" % nalvo)
    fb = rel(A, R) + "/mc_*.csv"
    eb = estima_ai(rsb)
    num("NBanda", "%d" % eb["N"], fb)
    num("SemFatorBanda", "%d" % len(ruinsb), fb)
    num("PfBanda", v(100. * eb["pf"], 2), fb + " (soma(peso * falha) / N)")
    num("CovPfBanda", v(100. * eb["cov"], 1), fb)
    num("PfBandaLo", v(100. * eb["pf"] * (1 - 1.96 * eb["cov"]), 1), fb + " (P_f (1 -+ 1,96 CoV))")
    num("PfBandaHi", v(100. * eb["pf"] * (1 + 1.96 * eb["cov"]), 1), fb)
    num("NfBanda", "%d" % eb["nf"], fb + " (falhas sorteadas)")
    num("FracFalhaBanda", v(100. * eb["nf"] / eb["N"], 0), fb)
    num("PesoMedioBanda", v(eb["w_media"], 3), fb)
    num("PesoMedioBandaSe", v(eb["w_dp"], 3), fb)
    num("EssBanda", "%d" % round(eb["ess"]), fb + " ((soma w)^2 / soma w^2)")
    num("PesoFalhaBanda", v(eb["wf_media"], 2), fb + " (peso medio nas falhas)")
    num("PesoMaxBanda", v(eb["wf_max"], 2), fb + " (maior peso nas falhas)")
    pb = np.column_stack([col(rsb, c) for c in ("proj_c", "proj_phi", "proj_k")]).mean(axis=0)
    num("ProjMediaBanda", "; ".join(v(x, 2) for x in pb), fb + " (media de proj_* sorteadas)",
        latex=r"\ensuremath{" + ";\\ ".join(v(x, 2).replace(",", "{,}") for x in pb) + "}")
    dif, gan, cov_d = compara(eb, p, cov_p)
    num("DifBanda", v(dif, 2), fb + " e " + fonte_d + " (desvios)")
    num("DifBandaPontos", v(100. * (eb["pf"] - p), 2), fb + " e " + fonte_d)
    num("GanhoBanda", v(gan, 1), fb + " ((CoV do direto com N e P_f da AI / CoV da AI)^2)")
    num("CovDiretoMesmoNBanda", v(100. * cov_d, 1), "Eq. 68 com N e P_f da AI")
    num("FracGanhoBanda", v(100. * gan * r, 0), "ganho medido / ganho linear")
    num("NDiretoParaCovBanda", mil((1 - eb["pf"]) / (eb["pf"] * eb["cov"] ** 2)), "Eq. 68 com o CoV da AI")
    txt("ViesBanda", "sem sinal de viés, abaixo de dois desvios" if abs(dif) < 2 else
        "acima de dois desvios, a conferir", fb)
    tb = le_tempos(A)
    num("CustoBanda", v(np.mean([tb[n][0] for n in tb if n >= 0]), 1), rel(A, R) + "/mc_*.log")
    res = os.path.join(R, "ai", "resultado.txt")
    if os.path.exists(res):  # a conferencia com o amostragem_importancia.py estima
        mr = re.search(r"P_f \(.*?\) = ([\d.]+)%  CoV\(P_f\) = ([\d.]+)%", texto(res))
        print("ai/resultado.txt: P_f = %s%%, CoV = %s%% (aqui: %.3f%%, %.1f%%)" % (mr[1], mr[2], 100 * eb["pf"],
                                                                               100 * eb["cov"]))
        assert abs(float(mr[1]) - 100 * eb["pf"]) < 2e-3
    ini, fim = le_tarefas(os.path.join(R, "ai", "logs"))
    T_ref = md["T2000"]
    ini_b = min(ini)
    fim_b = max(fim) if bpronto else eta(ini_b, nalvo, md["N"], T_ref)
    T_b = fim_b - ini_b
    txt("InicioBanda", hora(ini_b), "ai/logs/pp*.txt")
    txt("FimBanda", hora(fim_b, 0 if bpronto else 5), "ai/logs/pp*.txt" + ("" if bpronto else " (estimado)"))
    txt("RelogioBanda", duracao(T_b), "ai/logs/pp*.txt")
    # AI nos modos (ai_modos): deslocamento, rodada
    RM = os.path.join(R, "ai_modos")
    cad = texto(os.path.join(RM, "cadeia2.sh"))
    K = int(re.search(r"media_falhas=falhas.txt modos=(\d+)", cad)[1])
    nalvo_m = int(re.search(r"f50_I0_iii\s+\d+\s+(\d+)\s+banda=", cad)[1])
    sem_m = re.search(r"semente=(\d+)", cad)[1]
    num("KModos", "%d" % K, "ai_modos/cadeia2.sh")
    num("NAlvoModos", "%d" % nalvo_m, "ai_modos/cadeia2.sh")
    num("SementeModos", sem_m, "ai_modos/cadeia2.sh")
    num("RuidoModos", v(3. * K / nfp, 2), "3K / n_f")
    num("FatorRuidoModos", v(math.exp(3. * K / nfp), 2), "exp(3K / n_f)")
    V = os.path.join(R, "verificacao")
    dvec = np.array([float(x["d"]) for x in csv.DictReader(open(os.path.join(V, "banda_direcao.csv"), newline=""))])
    num("DirecaoK", v(100. * float((dvec[:K] ** 2).sum()), 0), "verificacao/banda_direcao.csv")
    arqm = os.path.join(RM, "modos_desloc_modos.csv")
    tem_desloc = os.path.exists(arqm)
    flag("modosdesloc", tem_desloc, "ai_modos/modos_desloc_modos.csv")
    mu2m = pu2 = None
    if tem_desloc:
        with open(arqm, newline="") as fm:
            mu = np.array([[float(r_["c"]), float(r_["phi"]), float(r_["k"])] for r_ in csv.DictReader(fm)])
        mu2m = float((mu ** 2).sum())
        pm = mu.T @ dvec[:len(mu)]
        pu2 = float((pm ** 2).sum())
        num("MuDoisModos", v(mu2m, 2), "ai_modos/modos_desloc_modos.csv")
        num("MuDoisModosUtil", v(pu2, 2), "ai_modos/modos_desloc_modos.csv e verificacao/banda_direcao.csv")
        num("MuDoisModosResto", v(mu2m - pu2, 2), "ai_modos/modos_desloc_modos.csv e verificacao/banda_direcao.csv")
        num("FatorRestoModos", v(math.exp(mu2m - pu2), 1), "exp(resto)")
        nfal = len([l for l in texto(os.path.join(RM, "falhas.txt")).split() if l.strip()])
        num("FalhasModos", "%d" % nfal, "ai_modos/falhas.txt")
    else:
        for mac in ("MuDoisModos", "MuDoisModosUtil", "MuDoisModosResto", "FatorRestoModos", "FalhasModos"):
            txt(mac, "--", "ai_modos/modos_desloc_modos.csv ainda nao existe")
    # a mesma decomposicao com o deslocamento por modos da verificacao (piloto de 400 realizacoes)
    vm = glob.glob(os.path.join(V, "modos_desloc_*.txt"))[0]
    lm = texto(vm).split("termo,c,phi,k\n")[1].strip().split("\n")
    muv = np.array([[float(x) for x in l_.split(",")[1:]] for l_ in lm])
    pv = muv.T @ dvec[:len(muv)]
    num("MuDoisVerif", v(float((muv ** 2).sum()), 2), rel(vm, R))
    num("MuDoisVerifUtil", v(float((pv ** 2).sum()), 2), rel(vm, R) + " e banda_direcao.csv")
    num("PilotoVerif", le_chaves(vm)["piloto_n"], rel(vm, R))
    num("FalhasVerif", le_chaves(vm)["falhas"], rel(vm, R))
    # os pesos sem elementos finitos (verificacao/)
    for tipo, mac in (("banda", "Banda"), ("modos", "Modos")):
        dv = le_chaves(glob.glob(os.path.join(V, "%s_desloc_*.txt" % tipo))[0])
        with open(glob.glob(os.path.join(V, "%s_proj_*.csv" % tipo))[0], newline="") as fq:
            w = np.exp(np.array([float(x["log_peso"]) for x in csv.DictReader(fq)]))
        f2 = "verificacao/%s_*" % tipo
        num("Sorteios" + mac, mil(len(w)), f2)
        num("PesoVerif" + mac, v(w.mean(), 3), f2)
        num("PesoVerifSe" + mac, v(w.std(ddof=1) / math.sqrt(len(w)), 3), f2)
        num("VarVerif" + mac, v(w.var(ddof=1), 2), f2)
        num("VarTeor" + mac, v(math.exp(float(dv["mu2"])) - 1., 2), f2 + " (exp(|mu|^2) - 1)")
    # a rodada da AI nos modos
    Mp = os.path.join(RM, "ceron_tr_fase1_ai_modos", "f50_I0_iii")
    rsm, ruinsm = le_mc(Mp) if os.path.isdir(Mp) else ({}, set())
    rsm.pop(-1, None)
    mpronto = len(rsm) + len(ruinsm) >= nalvo_m and os.path.exists(os.path.join(RM, "resultado.txt"))
    flag("modospronto", mpronto, "ai_modos/resultado.txt e %d realizacoes" % nalvo_m)
    flag("modosiniciado", len(rsm) >= 2, rel(Mp, R) + "/mc_*.csv")
    fm_ = rel(Mp, R) + "/mc_*.csv"
    em = None
    if len(rsm) >= 2:
        em = estima_ai(rsm)
        difm, ganm, cov_dm = compara(em, p, cov_p) if em["pf"] > 0 else (math.nan, math.nan, math.nan)
        num("NModos", "%d" % em["N"], fm_)
        num("SemFatorModos", "%d" % len(ruinsm), fm_)
        num("PfModos", v(100. * em["pf"], 2), fm_ + " (soma(peso * falha) / N)")
        num("CovPfModos", v(100. * em["cov"], 1) if em["pf"] > 0 else "--", fm_)
        num("NfModos", "%d" % em["nf"], fm_)
        num("FracFalhaModos", v(100. * em["nf"] / em["N"], 0), fm_)
        num("PesoMedioModos", v(em["w_media"], 3), fm_)
        num("PesoMedioModosSe", v(em["w_dp"], 3), fm_)
        num("EssModos", "%d" % round(em["ess"]), fm_)
        num("PesoFalhaModos", v(em["wf_media"], 2) if em["nf"] else "--", fm_)
        num("PesoMaxModos", v(em["wf_max"], 2) if em["nf"] else "--", fm_)
        num("DifModos", v(difm, 2) if em["pf"] > 0 else "--", fm_ + " e " + fonte_d + " (desvios)")
        num("GanhoModos", v(ganm, 1) if em["pf"] > 0 else "--", fm_)
        num("CovDiretoMesmoNModos", v(100. * cov_dm, 1) if em["pf"] > 0 else "--", "Eq. 68 com N e P_f da AI")
        txt("ViesModos", ("sem sinal de viés, abaixo de dois desvios" if abs(difm) < 2 else
                          "acima de dois desvios, a conferir") if em["pf"] > 0 else "--", fm_)
        txt("ComparaGanho", ("maior que" if ganm > gan else "menor que") if em["pf"] > 0 else "--", fm_ + " e " + fb)
    else:
        for mac in ("NModos", "SemFatorModos", "PfModos", "CovPfModos", "NfModos", "FracFalhaModos", "PesoMedioModos",
                    "PesoMedioModosSe", "EssModos", "PesoFalhaModos", "PesoMaxModos", "DifModos", "GanhoModos",
                    "CovDiretoMesmoNModos", "ViesModos", "ComparaGanho"):
            txt(mac, "--", "a AI nos modos ainda nao tem realizacoes")
    inim, fimm = le_tarefas(os.path.join(RM, "logs"))
    ini_m = min(inim) if inim else fim_b
    fim_m = max(fimm) if mpronto else eta(ini_m, nalvo_m, nalvo, T_b) if bpronto else eta(ini_m, nalvo_m, md["N"],
                                                                                         T_ref)
    txt("InicioModos", hora(ini_m, 0 if inim else 5), "ai_modos/logs/pp*.txt" + ("" if inim else " (estimado)"))
    txt("FimModos", hora(fim_m, 0 if mpronto else 5), "ai_modos/logs/pp*.txt" + ("" if mpronto else " (estimado)"))
    if em is not None:
        print("ai_modos: %d realizacoes ate agora" % em["N"])
    return {"eb": eb, "em": em, "ganho_lin": 1. / r, "fim_m": fim_m, "T_b": T_b, "bpronto": bpronto,
            "mpronto": mpronto}


def kfixo(R, md, fim_m, T_ref):
    """as mesmas realizacoes do direto com k fixo (covk=0): a comparacao pareada"""
    Kd = os.path.join(R, "covk0")
    idk, nalvo, opk = le_caso(os.path.join(Kd, "casos.txt"))
    num("NAlvoKfixo", "%d" % nalvo, "covk0/casos.txt")
    txt("CovKfixo", opk["covk"], "covk0/casos.txt")
    Kp = os.path.join(Kd, "ceron_tr_covk0", "f50_I0_iii")
    rk, ruk = le_mc(Kp) if os.path.isdir(Kp) else ({}, set())
    detk = rk.pop(-1, None)
    pronto = len(rk) + len(ruk) >= nalvo and os.path.exists(os.path.join(Kd, "resultado.txt"))
    flag("kfixopronto", pronto, "covk0/resultado.txt e %d realizacoes" % nalvo)
    al = md["al"]
    com = sorted(n for n in rk if n in al)
    flag("kfixoiniciado", len(com) >= 2, rel(Kp, R) + "/mc_*.csv")
    fk = rel(Kp, R) + "/mc_*.csv e direto"
    par = None
    if len(com) >= 2:
        a = np.array([float(rk[n]["FS_min"]) for n in com])
        b = np.array([float(al[n]["FS_min"]) for n in com])
        par = pareado(a, b)
        num("KNPar", "%d" % par["N"], fk)
        num("KSemFator", "%d" % len(ruk), rel(Kp, R) + "/mc_*.csv")
        for mac, k_ in (("KPfFixo", "pa"), ("KPfAleat", "pb")):
            num(mac, v(100. * par[k_], 2), fk)
            lo, hi = wilson(int(round(par[k_] * par["N"])), par["N"])
            num(mac + "Lo", v(100. * lo, 1), fk + " (Wilson)")
            num(mac + "Hi", v(100. * hi, 1), fk + " (Wilson)")
        num("KNzu", "%d" % par["n01"], fk + " (so com k aleatorio)")
        num("KNuz", "%d" % par["n10"], fk + " (so com k fixo)")
        t_, l_ = p_rel(par["p"], 3)
        num("KMcNemar", t_, fk + " (McNemar, com a relacao)", latex=l_)
        num("KDifPf", v(100. * par["dif"], 2), fk + " (pontos)")
        num("KDifPfSe", v(100. * par["se"], 2), fk + " (pontos)")
        num("KDifGamma", v(par["d_media"], 4), fk)
        num("KDifGammaSe", v(par["d_se"], 4), fk)
        num("KDifGammaMediana", v(par["d_mediana"], 4), fk)
        num("KMaisAleat", v(100. * par["d_mais"], 0), fk + " (% com Gamma maior com k aleatorio)")
        kb = np.array([float(al[n]["banda_k_mult"]) for n in com])
        km = np.array([float(al[n]["k_mult_med"]) for n in com])
        num("KRhoBanda", v(stats.spearmanr(kb, par["d"])[0], 2), fk + " (Spearman da diferenca com banda_k_mult)")
        num("KRhoSolo", v(stats.spearmanr(km, par["d"])[0], 2), fk + " (Spearman da diferenca com k_mult_med)")
        if detk is not None:
            num("KGammaDet", v(float(detk["FS_min"])), rel(Kp, R) + "/mc_0_49.csv")
        else:
            txt("KGammaDet", "--", "sem o deterministico ainda")
        signif = par["p"] < 0.05
        # o teste so mede o efeito do k_v: c' e phi' tem a estatistica do Ceron, e o resto da diferenca para ele fica
        # com a posicao da distribuicao de Gamma, que ainda depende da malha (Secoes 3.3 e 3.4 do relatorio)
        if signif and par["dif"] > 0:
            conc = ("A diferença é significativa: com o $k_v$ aleatório a $P_f$ sobe, e o $k_v$ explica parte da $P_f$ "
                    "maior do transiente.")
        elif signif:
            conc = "A diferença é significativa, mas no outro sentido: com o $k_v$ aleatório a $P_f$ cai."
        else:
            conc = ("A diferença não é significativa: o $k_v$ aleatório não explica a $P_f$ maior do transiente. O "
                    "excesso sobre o Ceron fica com a posição da distribuição de $\\Gamma$, que ainda depende da malha "
                    "(Seções~\\ref{sec:malha} e~\\ref{sec:comparacao}).")
        txt("KConclusao", conc, fk + " (McNemar p < 0,05?)")
    else:
        for mac in ("KNPar", "KSemFator", "KPfFixo", "KPfFixoLo", "KPfFixoHi", "KPfAleat", "KPfAleatLo", "KPfAleatHi",
                    "KNzu", "KNuz", "KMcNemar", "KDifPf", "KDifPfSe", "KDifGamma", "KDifGammaSe", "KDifGammaMediana",
                    "KMaisAleat", "KRhoBanda", "KRhoSolo", "KGammaDet", "KConclusao"):
            txt(mac, "--", "o k fixo ainda nao tem realizacoes")
    log = os.path.join(Kd, "cadeia3.log")
    ini = fim = None
    if os.path.exists(log):
        m = re.search(r"^(\S+ \S+) inicio", texto(log), re.M)
        ini = dt.datetime.strptime(m[1], "%Y-%m-%d %H:%M:%S") if m else None
        m = re.search(r"^(\S+ \S+) fim", texto(log), re.M)
        fim = dt.datetime.strptime(m[1], "%Y-%m-%d %H:%M:%S") if m else None
    ini_k = ini or fim_m
    fim_k = fim if (pronto and fim) else ini_k + T_ref * (nalvo / md["N"])
    txt("InicioKfixo", hora(ini_k, 0 if ini else 5), "covk0/cadeia3.log" + ("" if ini else " (estimado)"))
    txt("FimKfixo", hora(fim_k, 0 if (pronto and fim) else 5),
        "covk0/cadeia3.log" if (pronto and fim) else "inicio + o relogio do direto (estimado)")
    return par, (np.array([float(rk[n]["FS_min"]) for n in com]) if par else None), \
        (np.array([float(al[n]["FS_min"]) for n in com]) if par else None), fim_k


RX_TAREFA_SOLOS = re.compile(r"^(\d{4}-\d\d-\d\d \d\d:\d\d:\d\d) (\S+) (\d+)\.\.(\d+): (inicio|fim em (\d+) s)", re.M)


def solos_baixados(SD):
    """os tres solos baixados do fermi (CAMPANHA=solos), com as regras do figuras_solos.py (relatorio_solos.tex): entram
    as realizacoes com todos os instantes ok; a falha e FS < 1 em algum instante; o deterministico e a realizacao -1, na
    mesma malha. Sem a pasta, as macros ficam com '--' e a Secao 6.1 usa o status digitado"""
    R = os.path.join(SD, "chuva_solos")
    tem = os.path.isdir(R)
    flag("solosbaixados", tem, R if tem else "sem a pasta " + R)
    if not tem:
        flag("solosprontos", False, "sem a pasta " + R)
        for mac in MACROS_SOLOS_BAIXADOS:
            txt(mac, "--", "sem a campanha dos solos baixada")
        return
    casos = {}
    for linha in texto(os.path.join(PROJ, "cluster", "casos_solos.txt")).splitlines():
        if linha.strip() and not linha.startswith("#"):
            c = linha.split()
            casos[c[0]] = int(c[2])
    faltam_t, val_t, falt_n = [], 0, 0
    for sid, mac in (("f10", "Dez"), ("f50", "Cinquenta"), ("f100", "Cem")):
        pasta = os.path.join(R, sid)
        fonte = "chuva_solos/%s/mc_*.csv" % sid
        ts = [float(x) for x in le_chaves(os.path.join(pasta, "mc_parametros.txt"))["TEMPOS"].split(",")]
        per = {}
        for arq in sorted(glob.glob(os.path.join(pasta, "mc_*_*.csv"))):
            if not re.search(r"/mc_-?\d+_-?\d+\.csv$", arq):
                continue  # mapas
            with open(arq, newline="") as f:
                for r in csv.DictReader(f):
                    try:
                        n, t = int(r["real"]), float(r["t_h"])
                        float(r["FS_min"]), float(r["proj_k"])  # a linha chegou ao fim
                    except (TypeError, ValueError, KeyError):
                        continue
                    per.setdefault(n, {})[t] = r
        det = per.pop(-1)
        comp = {n: d for n, d in per.items() if all(t in d for t in ts)}
        ok = lambda r: r["situacao"] == "ok" and float(r["FS_min"]) > 0  # noqa: E731
        ruins = [n for n, d in comp.items() if not all(ok(d[t]) for t in ts)]
        mins = [min(float(d[t]["FS_min"]) for t in ts) for n, d in comp.items() if n not in ruins]
        N, nf = len(mins), sum(x < 1. for x in mins)
        lo, hi = wilson(nf, N)
        falham = sum(min([float(comp[n][t]["FS_min"]) for t in ts if ok(comp[n][t])] or [math.inf]) < 1. for n in ruins)
        faltam = sorted(set(range(casos[sid])) - set(comp))
        num("SolosValidas" + mac, "%d" % N, fonte + " (todos os instantes ok)")
        num("SolosCompletas" + mac, "%d" % len(comp), fonte + " (todos os instantes)")
        num("SolosExcluidas" + mac, "%d" % len(ruins), fonte + " (situacao diferente de ok)")
        num("SolosExcluidasFalham" + mac, "%d" % falham, fonte + " (as excluidas com FS < 1 antes de parar)")
        num("SolosFaltam" + mac, "%d" % len(faltam), fonte + " (sem todas as linhas)")
        num("SolosNf" + mac, "%d" % nf, fonte + " (FS < 1 em algum instante)")
        num("SolosPfFinal" + mac, v(100. * nf / N, 2), fonte)
        num("SolosPfFinalLo" + mac, v(100. * lo, 1), fonte + " (Wilson)")
        num("SolosPfFinalHi" + mac, v(100. * hi, 1), fonte + " (Wilson)")
        dt_ = [float(det[t]["FS_min"]) for t in ts]
        dmin = min(dt_)
        tdet = ["%g" % t for t, x in zip(ts, dt_) if x <= dmin * (1. + 1e-6)]
        num("SolosDetFinal" + mac, v(dmin, 3), "chuva_solos/%s/mc_0_4.csv (real -1, menor FS)" % sid)
        txt("SolosTdetFinal" + mac, (", ".join(tdet[:-1]) + " e " + tdet[-1] if len(tdet) > 1 else tdet[0]) + "~h",
            "chuva_solos/%s/mc_0_4.csv (real -1, instante do menor FS)" % sid)
        val_t += N
        falt_n += len(faltam)
        if faltam:
            faltam_t.append("%s no %s" % (", ".join(map(str, faltam[:-1])) + " e " + str(faltam[-1])
                                          if len(faltam) > 1 else str(faltam[0]), sid))
    num("SolosValidasTotal", "%d" % val_t, "chuva_solos/*/mc_*.csv")
    num("SolosFaltamTotal", "%d" % falt_n, "chuva_solos/*/mc_*.csv")
    txt("SolosFaltamLista", " e ".join(faltam_t) if faltam_t else "nenhuma", "chuva_solos/*/mc_*.csv")
    # as tarefas (logs do cluster/tarefas.sh): inicio, fim, CPU e as que ainda rodavam quando a pasta foi baixada
    tar = {}
    for arq in glob.glob(os.path.join(R, "logs", "*.txt")):
        for m in RX_TAREFA_SOLOS.finditer(texto(arq)):
            x = tar.setdefault((m[2], int(m[3]), int(m[4])), [None, None, None])
            t = dt.datetime.strptime(m[1], "%Y-%m-%d %H:%M:%S")
            if m[5] == "inicio":
                x[0] = t
            else:
                x[1], x[2] = t, int(m[6])
    feitas = {k: x for k, x in tar.items() if x[1]}
    rodando = {k: x for k, x in tar.items() if not x[1]}
    ini = min(x[0] for x in tar.values() if x[0])
    fim = max(x[1] for x in feitas.values())
    bloco = int(texto(os.path.join(R, "campanha.txt")).split()[1])
    ftar = "chuva_solos/logs/p*.txt"
    txt("SolosInicioLog", hora(ini), ftar)
    txt("SolosUltimo", hora(fim), ftar + " (fim da ultima tarefa terminada)")
    txt("SolosRelogioFinal", duracao(fim - ini), ftar)
    cpu = sum(x[2] for x in feitas.values()) / 3600.
    num("SolosCpuFinal", "%d" % round(cpu), ftar + " (soma das tarefas terminadas)")
    num("SolosMinMedido", v(cpu * 60. / sum(k[2] - k[1] + 1 for k in feitas), 1), ftar + " (CPU / realizacoes)")
    num("SolosRodando", "%d" % len(rodando), ftar + " (inicio sem fim)")
    num("SolosTarefas", "%d" % sum(math.ceil(n / bloco) for n in casos.values()), "casos_solos.txt e campanha.txt")
    if rodando:  # inicio + a duracao media das tarefas do mesmo solo
        med = {sid: np.mean([x[2] for k, x in feitas.items() if k[0] == sid]) for sid in casos}
        prev = max(x[0] + dt.timedelta(seconds=med[k[0]]) for k, x in rodando.items())
        txt("SolosPrevFimFinal", hora(prev, 5), ftar + " (inicio + duracao media das tarefas do solo, estimado)")
    else:
        txt("SolosPrevFimFinal", "--", ftar)
    flag("solosprontos", falt_n == 0 and not rodando, "chuva_solos: todas as realizacoes e nenhuma tarefa rodando")
    print("tres solos baixados: %d realizacoes validas, faltam %d, %d tarefas rodando" % (val_t, falt_n, len(rodando)))


def fermi(B, md, SD):
    """os tres solos (os arquivos baixados, se existem; o status digitado sempre, para o texto sem eles) e as campanhas
    prontas (casos e custos)"""
    S = SOLOS
    h0 = dt.datetime.combine(DATA_RELATORIO, dt.time(*S["hora"]))
    i0 = dt.datetime.combine(DATA_RELATORIO, dt.time(*S["inicio"]))
    txt("SolosHora", hora(h0), FONTE_SOLOS)
    txt("SolosInicio", hora(i0), FONTE_SOLOS)
    num("SolosNos", "%d" % S["nos"], FONTE_SOLOS)
    num("SolosPrevMin", "%d" % S["previsao_h"][0], FONTE_SOLOS)
    num("SolosPrevMax", "%d" % S["previsao_h"][1], FONTE_SOLOS)
    f0, f1 = (h0 + dt.timedelta(hours=x) for x in S["previsao_h"])
    a0, a1 = hora(f0, 30), hora(f1, 30)
    if f0.date() == f1.date() and " de " in a0:  # o dia uma vez so: "13h30 e 16h30 de 15/11"
        a0 = a0.split(" de ")[0]
    txt("SolosFimFaixa", "%s e %s" % (a0, a1), FONTE_SOLOS + " + previsao")
    cs = os.path.join(PROJ, "cluster", "casos_solos.txt")
    tot = 0
    for linha in texto(cs).splitlines():
        if linha.strip() and not linha.startswith("#"):
            c = linha.split()
            sid, N = c[0], int(c[2])
            opc = dict(o.split("=", 1) for o in c[3:])
            mac = {"f10": "Dez", "f50": "Cinquenta", "f100": "Cem"}[sid]
            st = S["casos"][sid]
            num("SolosN" + mac, "%d" % N, "cluster/casos_solos.txt")
            txt("SolosTempos" + mac, opc["tempos"].replace(",", ", "), "cluster/casos_solos.txt")
            num("SolosFeitas" + mac, "%d" % st["n"], FONTE_SOLOS)
            num("SolosFalhas" + mac, "%d" % st["nf"], FONTE_SOLOS)
            num("SolosPf" + mac, v(100. * st["nf"] / st["n"], 2), FONTE_SOLOS)
            lo, hi = wilson(st["nf"], st["n"])
            num("SolosLo" + mac, v(100. * lo, 1), FONTE_SOLOS + " (Wilson)")
            num("SolosHi" + mac, v(100. * hi, 1), FONTE_SOLOS + " (Wilson)")
            num("SolosDet" + mac, v(st["det"], 3), FONTE_SOLOS)
            txt("SolosTdet" + mac, st["tdet"], FONTE_SOLOS)
            tot += N
    mk = re.search(r"I = (\d+) k_s do solo por (\d+) h", texto(cs))
    num("SolosIks", mk[1], "cluster/casos_solos.txt")
    num("SolosDuracao", mk[2], "cluster/casos_solos.txt")
    ts = [int(x) for x in opc["tempos"].split(",")]
    num("SolosPasso", "%d" % (ts[1] - ts[0]), "cluster/casos_solos.txt (tempos)")
    num("SolosFim", "%d" % ts[-1], "cluster/casos_solos.txt (tempos)")
    mI = re.search(r"f10 a ([\d,]+) mm/h, f50 a ([\d,]+) mm/h e f100 a ([\d,]+) mm/h", texto(cs))
    for mac, x in zip(("Dez", "Cinquenta", "Cem"), mI.groups()):
        num("SolosI" + mac, x, "cluster/casos_solos.txt")
    num("SolosTotal", "%d" % tot, "cluster/casos_solos.txt")
    num("SolosFeitasTotal", "%d" % sum(c["n"] for c in S["casos"].values()), FONTE_SOLOS)
    readme = texto(os.path.join(PROJ, "cluster", "README.md"))
    num("SolosMinReal", re.search(r"~(\d+) min por realização no fermi; as \d+ realizações", readme)[1],
        "cluster/README.md")
    solos_baixados(SD)
    # submeter.sh: nos x processos
    sub = texto(os.path.join(PROJ, "cluster", "submeter.sh"))
    nos = int(re.search(r'NOS="\$\{1:-(\d+)\}"', sub)[1])
    paralelo = int(re.search(r'PARALELO="\$\{PARALELO:-(\d+)\}"', sub)[1])
    proc = nos * paralelo
    num("Nos", "%d" % nos, "cluster/submeter.sh")
    num("Paralelo", "%d" % paralelo, "cluster/submeter.sh")
    num("FermiLento", v(FERMI_LENTO, 1), "medido antes (digitado)")
    cf = md["custo"] * FERMI_LENTO
    num("CustoFermi", v(cf, 0), "custo do PC x 2,6")
    # CAMPANHA=ceron
    cc = os.path.join(PROJ, "cluster", "casos_ceron.txt")
    Nc = {}
    for linha in texto(cc).splitlines():
        if linha.strip() and not linha.startswith("#"):
            c = linha.split()
            Nc[c[0]] = int(c[2])
    num("CeronDireto", mil(Nc["ref_direto"]), "cluster/casos_ceron.txt")
    num("CeronBanda", mil(Nc["ref_ai_banda"]), "cluster/casos_ceron.txt")
    num("CeronModos", mil(Nc["ref_ai_modos"]), "cluster/casos_ceron.txt")
    Ntc = sum(Nc.values())
    num("CeronTotal", mil(Ntc), "cluster/casos_ceron.txt")
    cpu = Ntc * cf / 3600.
    num("CeronCpu", "%d" % round(cpu), "casos_ceron.txt x custo no fermi")
    num("CeronRelogio", v(cpu / proc, 0), "CPU / (%d x %d)" % (nos, paralelo))
    pilc = re.search(r"piloto_n=(\d+)", texto(cc))[1]
    num("CeronPiloto", pilc, "cluster/casos_ceron.txt")
    num("CeronPilotoFalhas", "%d" % round(int(pilc) * md["p"]), "piloto_n x P_f do direto (esperado)")
    # CAMPANHA=chuva_ai
    ca = os.path.join(PROJ, "cluster", "casos_chuva_ai.txt")
    ta = texto(ca)
    Na = {}
    for linha in ta.splitlines():
        if linha.strip() and not linha.startswith("#"):
            c = linha.split()
            Na[c[0]] = (int(c[2]), dict(o.split("=", 1) for o in c[3:]))
    num("ChuvaAIPiloto", "%d" % Na["ref48_piloto"][0], "cluster/casos_chuva_ai.txt")
    num("ChuvaAIN", "%d" % Na["ref48_ai"][0], "cluster/casos_chuva_ai.txt")
    tempos = Na["ref48_ai"][1]["tempos"].split(",")
    txt("ChuvaAITempos", "%s a %s" % (tempos[0], tempos[-1]), "cluster/casos_chuva_ai.txt")
    num("ChuvaAIPasso", "%d" % (int(tempos[1]) - int(tempos[0])), "cluster/casos_chuva_ai.txt")
    num("ChuvaAIFim", tempos[-1], "cluster/casos_chuva_ai.txt")
    min_ai = int(re.search(r"~(\d+) min por realizacao no fermi", ta)[1])
    num("ChuvaAIMin", "%d" % min_ai, "cluster/casos_chuva_ai.txt")
    # blocos de 1 realizacao e o ref48_ai so comeca quando o piloto acaba: rodadas de 48 realizacoes, uma depois da outra
    rod_ai = math.ceil(Na["ref48_piloto"][0] / proc) + math.ceil(Na["ref48_ai"][0] / proc)
    rel_ai = rod_ai * min_ai / 60.
    num("ChuvaAIRodadas", "%d" % rod_ai, "casos_chuva_ai.txt: ceil(piloto / %d) + ceil(AI / %d)" % (proc, proc))
    num("ChuvaAIRelogio", v(rel_ai, 1), "rodadas x %d min" % min_ai)
    m = re.search(r"CoV\(P_f\) de ~(\d+)% \(estado limite linear\) a ~(\d+)% com (\d+) realizações", readme)
    num("ChuvaAICovLin", m[1], "cluster/README.md")
    num("ChuvaAICovMax", m[2], "cluster/README.md")
    iref = float(re.search(r"^chuva_ai\).*?IREF_PADRAO=([\d.]+)", texto(os.path.join(PROJ, "cluster", "campanha.sh")),
                           re.M | re.S)[1])
    num("Iref", v(iref, 2), "cluster/campanha.sh")
    # CAMPANHA=chuva (casos.txt), com o custo em serie da estacao (mc/serial_I24.log), como no relatorio_chuva.tex
    mc = os.path.join(B, "mc")
    lg = texto(os.path.join(mc, "serial_I24.log"))
    mt = re.search(r"tempo medio por realizacao: transiente ([\d.]+) s, fator ([\d.]+) s", lg)
    with open(os.path.join(mc, "serial_I24_0_0.csv"), newline="") as fs:
        ninst = len(list(csv.DictReader(fs)))
    tT, tF = float(mt[1]), float(mt[2]) / ninst
    fases = {}
    for linha in texto(os.path.join(PROJ, "cluster", "casos.txt")).splitlines():
        if linha.strip() and not linha.startswith("#"):
            c = linha.split()
            opc = dict(o.split("=", 1) for o in c[3:])
            fase, N = int(c[1]), int(c[2])
            if c[0] == "ref":
                tr = opc["tempos"].split(",")
                txt("ChuvaRefTempos", "%s a %s" % (tr[0], tr[-1]), "cluster/casos.txt (ref)")
            custo = tT + len(opc["tempos"].split(",")) * tF
            x = fases.setdefault(fase, [0, 0., 0])
            x[0] += N
            x[1] += N * custo * FERMI_LENTO / 3600.
            x[2] += 1
    for fase, nome in ((1, "Um"), (2, "Dois"), (3, "Tres")):
        num("ChuvaFase%sN" % nome, mil(fases[fase][0]), "cluster/casos.txt")
        num("ChuvaFase%sCasos" % nome, "%d" % fases[fase][2], "cluster/casos.txt")
        num("ChuvaFase%sCpu" % nome, "%d" % round(fases[fase][1]), "cluster/casos.txt + mc/serial_I24.log x 2,6")
        num("ChuvaFase%sRelogio" % nome, v(fases[fase][1] / proc, 1), "CPU / %d processos" % proc)
    T = sum(x[1] for x in fases.values())
    num("ChuvaTotalDias", v(T / proc / 24., 1), "cluster/casos.txt + mc/serial_I24.log")
    num("ChuvaRelogioFaseUmReadme", re.search(r"a fase 1 leva ~(\d+) h", readme)[1], "cluster/README.md")
    return {"cf": cf, "proc": proc, "min_ai": min_ai, "rel_ai": rel_ai, "nos": nos, "paralelo": paralelo,
            "n_ai": (Na["ref48_piloto"][0], Na["ref48_ai"][0])}


def sugestoes(md, fr, ant_dir):
    """os custos das sugestoes (Secao 7)"""
    cf, proc = fr["cf"], fr["proc"]
    # (b) as Secoes 6.2 a 6.5 do Ceron com o transiente: os casos da campanha antiga (pastas f50_I0_iii_*), menos os de
    # CoV(k_v) (o k fixo roda no PC; os outros vao para a familia do item d); 2000 realizacoes por caso
    casos = sorted(os.path.basename(d) for d in glob.glob(os.path.join(ant_dir, "ceron_mc_fase[23]", "f50_I0_iii_*"))
                   if os.path.isdir(d) and not d.endswith("_det"))
    pfs = {}
    for c in casos:
        arq = glob.glob(os.path.join(ant_dir, "ceron_mc_fase[23]", c, "mc_resumo.csv"))[0]
        with open(arq, newline="") as f:
            G = np.array([float(r["FS_min"]) for r in csv.DictReader(f) if r["situacao"] == "ok"])
        pfs[c] = (G < 1).mean()
    # os de CoV(k_v) ficam na familia do item 4 (que so vale a pena se o k fixo mostrar efeito); o covk0 roda no PC
    covk = [c for c in casos if c.startswith("f50_I0_iii_covk") and c != "f50_I0_iii_covk0"]
    novos = [c for c in casos if not c.startswith("f50_I0_iii_covk")]
    Npc = 2000
    num("SugCasos", "%d" % len(novos), "pastas f50_I0_iii_* da campanha antiga, sem os covk*")
    num("SugNporCaso", "%d" % Npc, "sugestao")
    custo4 = md["custo"] + md["custo_fator"]  # nref = 4: o fator umas 2 vezes mais caro (estimado)
    cpu_b = len(novos) * Npc * cf / 3600.
    cpu_b4 = cpu_b * custo4 / md["custo"]
    num("SugBCpu", "%d" % round(cpu_b), "casos x 2000 x custo no fermi (nref = 3)")
    num("SugBRelogio", v(cpu_b / proc, 0), "CPU / %d processos" % proc)
    num("SugBCpuQuatro", "%d" % (100 * round(cpu_b4 / 100.)), "o mesmo com o custo de nref = 4 (estimado)")
    num("SugBRelogioQuatro", v(cpu_b4 / proc, 0), "CPU / %d processos" % proc)
    limiar = 0.03  # abaixo disso, AI (CoV(P_f) do direto com 2000 > 12%)
    num("SugLimiar", "%g" % (100. * limiar), "sugestao")
    pequenos = sorted((c for c in novos if pfs[c] < limiar), key=lambda c: pfs[c])
    num("SugPequenos", "%d" % len(pequenos), "campanha antiga: P_f < 3%")
    nomes = {"covc0.1": "$CoV(c') = 0{,}1$", "alfa3": r"$\alpha = 3$", "alfa4": r"$\alpha = 4$",
             "alfa5": r"$\alpha = 5$", "nivelfinal7.5": "$h_w/H = 0{,}5$", "nivelfinal7": "$h_w/H = 0{,}6$",
             "alfa2": r"$\alpha = 2$"}
    lista = ["%s com %s\\%%" % (nomes.get(c.replace("f50_I0_iii_", ""), c.replace("_", r"\_")),
                                v(100. * pfs[c], 2).replace(",", "{,}")) for c in pequenos]
    txt("SugPequenosLista", ", ".join(lista[:-1]) + " e " + lista[-1] if len(lista) > 1 else "".join(lista),
        "campanha antiga (mc_resumo.csv de cada caso)")
    num("SugDiretos", "%d" % (len(novos) - len(pequenos)), "campanha antiga: P_f >= 3%")
    pf_max = max(pfs[c] for c in novos)
    num("SugPfMaxAntiga", v(100. * pf_max, 1), "campanha antiga")
    pf_min_dir = min(pfs[c] for c in novos if pfs[c] >= limiar)
    num("SugPfMinDireto", v(100. * pf_min_dir, 1), "campanha antiga")
    num("SugCovPfMax", v(100. * math.sqrt((1 - pf_min_dir) / (Npc * pf_min_dir)), 0), "Eq. 68, N = 2000")
    # (c) malha: o deterministico e 50 realizacoes com nref = 4 (custo do fator estimado em 2x o do nivel 3), em blocos
    # de ceil(50 / 8) = 7 (um por processo); o bloco que comeca na realizacao 0 roda tambem o deterministico
    Nm = 50
    bloco_m = math.ceil(Nm / md["nproc"])
    num("SugMalhaN", "%d" % Nm, "sugestao")
    num("SugMalhaBloco", "%d" % bloco_m, "ceil(%d / %d processos)" % (Nm, md["nproc"]))
    num("SugMalhaCusto", v(custo4, 0), "custo do PC com o fator 2 vezes mais caro (estimado)")
    num("SugMalhaRelogio", v((bloco_m + 1) * custo4 / 60., 0), "(BLOCO + o deterministico) x custo")
    # (d) a familia de CoV(k_v): 0,3 e os casos covk* da campanha antiga (a Tabela 3 do Ceron), 2000 cada
    ks = [0.3] + sorted(float(c.replace("f50_I0_iii_covk", "")) for c in covk)
    lk = ["%g" % x for x in ks]
    lk = [x.replace(".", "{,}") for x in lk]
    txt("SugKLista", "; ".join(lk[:-1]) + " e " + lk[-1], "0,3 e as pastas covk* da campanha antiga")
    num("SugKCasos", "%d" % len(ks), "0,3 e as pastas covk* da campanha antiga")
    num("SugKPontos", "%d" % (len(ks) + 2), "com o k_v fixo (0) e o direto (0,6)")
    Nd = len(ks) * Npc
    num("SugKN", mil(Nd), "sugestao: %d casos x %d" % (len(ks), Npc))
    num("SugKRelogioPC", v(Nd * md["custo"] / md["nproc"] / 3600., 1), "custo do PC / %d processos" % md["nproc"])
    num("SugKRelogioFermi", v(Nd * cf / 3600. / proc, 1), "custo no fermi / %d processos" % proc)
    # (e) curva de fragilidade da chuva: 5 intensidades com piloto e AI de 200 (os casos da chuva_ai), todas juntas:
    # as rodadas dos pilotos e depois as das AI
    nI = 5
    npil, nai = fr["n_ai"]
    rod = math.ceil(nI * npil / proc) + math.ceil(nI * nai / proc)
    num("SugFragI", "%d" % nI, "sugestao")
    num("SugFragRelogio", v(rod * fr["min_ai"] / 60., 0), "%d rodadas de %d min" % (rod, fr["min_ai"]))
    return pfs


def main(argv):
    R = os.path.abspath(argv[1]) if len(argv) > 1 else R_PADRAO
    ANT = os.path.abspath(argv[2]) if len(argv) > 2 else ANTIGA_PADRAO
    B = os.path.abspath(argv[3]) if len(argv) > 3 else BUILD_PADRAO
    SD = os.path.abspath(argv[4]) if len(argv) > 4 else SOLOS_PADRAO
    pasta = os.path.join(AQUI, "figs")
    os.makedirs(pasta, exist_ok=True)
    D = os.path.join(R, "direto", "ceron_tr_fase1", "f50_I0_iii")
    DET = os.path.join(R, "direto", "det")
    pts, H, gmain, par = modelo(R, D, DET)
    md = direto(R, D, DET, pts)
    num("DifDetMain", v(100. * (md["gdet"] / gmain - 1.), 2), "direto/det/mc_0_-1.csv e rebaixamento_mc.cpp")
    an = antiga(ANT, md)
    ra = ai(R, md)
    par_k, a_k, b_k, fim_k = kfixo(R, md, ra["fim_m"], md["T2000"])
    fr = fermi(B, md, SD)
    sugestoes(md, fr, ANT)
    # ---- figuras
    fig_mc(pasta, md["G"], md["gdet"], md["lnmu"], md["lnsig"], an["G"], md["n"], md["pfn"], md["sen"], an["p"],
           an["N"])
    fig_malha(pasta, md["med"], md["pfs"], an["med"], an["pfs"], md["niv_det"], an["det"], md["ext_med"], md["pf4"])
    fig_mecanismo(pasta, pts, md["cel"], md["N"], md["nf"], md["banda_det"], md["rhos"])
    ais = [("AI na banda" + ("" if ra["bpronto"] else " (parcial)"), COR_BANDA, ra["eb"])]
    if ra["em"] is not None and ra["em"]["N"] >= 50:
        ais.append(("AI nos modos" + ("" if ra["mpronto"] else " (parcial)"), COR_MODOS, ra["em"]))
    p = md["p"]
    cov_lin = math.sqrt((1. / ra["ganho_lin"]) * (1 - p) / p)
    fig_ai(pasta, md["f"].astype(float), ais, cov_lin)
    if par_k is not None and FLAGS["kfixopronto"]:
        fig_kfixo(pasta, a_k, b_k, par_k)
    # ---- numeros_rebaixamento.tex
    arq = os.path.join(AQUI, "numeros_rebaixamento.tex")
    with open(arq, "w") as f:
        f.write("%% gerado por figuras_rebaixamento.py em %s: nao edite (rode o script de novo)\n"
                % dt.datetime.now().strftime("%d/%m/%Y %H:%M"))
        for nome, val in FLAGS.items():
            f.write("\\newif\\if%s \\%s%s\n" % (nome, nome, "true" if val else "false"))
        for macro, valor, latex, fonte in NUMEROS:
            if latex is not None:
                f.write("\\newcommand{\\%s}{%s} %% %s\n" % (macro, latex, fonte.replace("%", r"\%")))
    print("numeros: %s" % arq)
    print("\n% ---- numeros do texto (macro: valor  [arquivo])")
    for macro, valor, latex, fonte in NUMEROS:
        print("%-24s %-34s [%s]" % (macro, valor, fonte))


if __name__ == "__main__":
    main(sys.argv)
