#!/usr/bin/env python3
"""figuras_solos.py: as figuras, as tabelas e os numeros do relatorio_solos.tex: o Monte Carlo dos tres solos de
Rahardjo et al. (2007) sob chuva (I = 1 k_s por 24 h), com os campos aleatorios de Vargas Ceron et al. (2025), rodado no
fermi (CAMPANHA=solos) em 14 e 15/11/2026. So le arquivos; nao roda nenhum programa e nao escreve nas pastas das
rodadas:
  - a campanha baixada pelo cluster/baixar.sh (~/projects/runs/fermi_chuva_solos): chuva_solos/f10, f50 e f100
    (mc_<a>_<b>.csv, mc_<a>_<b>.log, mc_mapa_<a>_<b>.csv, mc_parametros.txt e mc_real-1/), chuva_solos/logs/ (inicio e
    fim de cada tarefa), chuva_solos/campanha.txt e logs/ (a saida dos jobs do SLURM);
  - no build: o FS deterministico de hora em hora do SlopeRahardjo (rahardjo_f10.log, rahardjo_f50.log e
    rahardjo_f100.log, com o transiente na TriGMesh(3)) e o f50 com I = 31,05 mm/h (o deterministico,
    mc/gim_I/srm_I31.05_0_-1.csv, e o piloto de 8 realizacoes, mc/gim_I/piloto_srm_I31.05_0_7.csv e
    mc/gim_I/piloto_srm.log, com a P_f lognormal);
  - as fontes: montecarlo.cpp (talude, resistencia e limiar da banda), rahardjo.cpp (malha do transiente),
    cluster/casos_solos.txt, cluster/casos_chuva_ai.txt, cluster/campanha.sh, cluster/README.md,
    posproc/dados_chuva.py (tolerancia das categorias do nivel d'agua) e figuras_relatorio.py (as leituras da Fig. 5 do
    Rahardjo, RAHARDJO_FIG5).

Escreve:
  figs/fig_solos_*.pdf   as figuras
  numeros_solos.tex      um \\newcommand para cada numero do texto e das tabelas e a chave \\ifsolosprontos
                         (verdadeira com todas as realizacoes nos arquivos e nenhum bloco rodando)
e imprime cada numero com o arquivo de onde sai. Nenhum resultado e digitado aqui; as unicas leituras de artigo sao as
da Fig. 5 do Rahardjo, tiradas do figuras_relatorio.py.

Regras (as do posproc/resumo_tempos.py): uma realizacao entra se tem todos os instantes com situacao ok e FS > 0; a
falha e FS < 1 em algum instante (o minimo no tempo); P_f = n_f / N, com o intervalo de Wilson de 95%. A realizacao -1
e o deterministico (as medias em todo o solo), na mesma malha do Monte Carlo.

Uso: python3 -I figuras_solos.py [campanha baixada] [build do SlopeRandomFieldRefinement]
     (padrao: ~/projects/runs/fermi_chuva_solos e
     ../../../../neopz-master-build-release/Projects/SlopeRandomFieldRefinement)"""

import ast
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
from matplotlib.legend_handler import HandlerTuple  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
from scipy import stats  # noqa: E402

AQUI = os.path.dirname(os.path.abspath(__file__))
PROJ = os.path.dirname(AQUI)  # Projects/SlopeRandomFieldRefinement (fontes)
NEOPZ = os.path.dirname(os.path.dirname(PROJ))
CAMPANHA_PADRAO = os.path.expanduser("~/projects/runs/fermi_chuva_solos")
BUILD_PADRAO = os.path.join(os.path.dirname(NEOPZ), "neopz-master-build-release", "Projects",
                            "SlopeRandomFieldRefinement")
DATA_RELATORIO = dt.date(2026, 11, 15)

SOLOS = ("f10", "f50", "f100")
NOME = {"f10": "Dez", "f50": "Cinquenta", "f100": "Cem"}  # o sufixo das macros de cada solo
ROTULO = {"f10": "$f_{10,-4}$", "f50": "$f_{50,-5}$", "f100": "$f_{100,-6}$"}
INST = {0: "Zero", 6: "Seis", 12: "Doze", 18: "Dezoito", 24: "VinteQuatro", 30: "Trinta", 36: "TrintaSeis",
        42: "QuarentaDois", 48: "QuarentaOito"}  # o sufixo das macros de cada instante (h)

# ---- estilo (o do figuras_rebaixamento.py e do figuras_relatorio.py)
PALETA = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
TINTA, TINTA2, MUDO, GRADE, SOLO = "#0b0b0b", "#52514e", "#898781", "#dddcd8", "#f3f2ee"
COR = {"f10": PALETA[0], "f50": PALETA[1], "f100": PALETA[2]}
plt.rcParams.update({
    "font.size": 8.5, "axes.labelsize": 8.5, "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "legend.fontsize": 7,
    "axes.edgecolor": TINTA2, "axes.labelcolor": TINTA, "xtick.color": TINTA2, "ytick.color": TINTA2,
    "axes.grid": True, "grid.color": GRADE, "grid.linewidth": 0.5, "lines.linewidth": 1.4, "legend.frameon": False,
    "savefig.dpi": 200, "font.family": "DejaVu Sans", "mathtext.fontset": "dejavusans",
})
LARGURA = 6.5  # polegadas: a largura do texto do relatorio (A4 com margens de 2,2 cm)

NUMEROS = []  # (macro, valor impresso, valor no .tex, arquivo)


# ---------------------------------------------------------------- numeros
def v(x, casas=4):
    """numero com virgula decimal (sem o sinal de -0,0)"""
    t = "%.*f" % (casas, x)
    if t.startswith("-") and float(t) == 0.:
        t = t[1:]
    return t.replace(".", ",")


def vs(x, casas=1):
    """numero com sinal (+1,2 ou -0,8), para as diferencas em %"""
    t = v(x, casas)
    return t if t.startswith("-") or float(t.replace(",", ".")) == 0. else "+" + t


def g(x):
    return ("%g" % x).replace(".", ",")


def acima_abaixo(x, casas=1):
    """uma variacao relativa em %: '2{,}4\\% abaixo' ou '1{,}3\\% acima'"""
    return r"%s\%% %s" % (v(abs(x), casas).replace(",", "{,}"), "abaixo" if x < 0 else "acima")


def num(macro, texto, fonte, latex=None):
    """registra um numero do texto: \\macro no .tex (virgula decimal como {,}, o sinal de menos em modo matematico) e o
    valor impresso com o arquivo de onde sai. latex: o valor no .tex, quando nao e um numero simples"""
    assert re.fullmatch(r"[A-Za-z]+", macro), macro
    assert macro not in {m for m, _, _, _ in NUMEROS}, "macro repetida: " + macro
    if latex is None:
        latex = re.sub(r"(?<!\\),", "{,}", texto)  # a virgula decimal
        if re.match(r"[-+]\d", latex):
            latex = r"\ensuremath{%s}" % ("-" if latex[0] == "-" else "+") + latex[1:]
    NUMEROS.append((macro, texto, latex, fonte))
    return texto


def txt(macro, texto, fonte):
    """um texto que depende dos dados (uma lista, um horario): vai para o .tex como esta"""
    return num(macro, texto, fonte, latex=texto)


def wilson(k, n, z=1.959964):
    p = k / n
    c = (p + z * z / (2 * n)) / (1 + z * z / n)
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return max(c - h, 0.), min(c + h, 1.)


def lista(itens, sep=", "):
    """'a', 'a e b', 'a, b e c' (sep = "; " para numeros com virgula decimal)"""
    itens = [str(x) for x in itens]
    return itens[0] if len(itens) == 1 else sep.join(itens[:-1]) + " e " + itens[-1]


def hora(t, arred=0):
    """'8h40' (arred: arredonda para esse numero de minutos); com o dia quando nao e o do relatorio"""
    if arred:
        m = int(round((t.hour * 60 + t.minute + t.second / 60.) / arred) * arred)
        t = dt.datetime.combine(t.date(), dt.time()) + dt.timedelta(minutes=m)
    s = "%dh%02d" % (t.hour, t.minute)
    return s if t.date() == DATA_RELATORIO else s + " de %d/%d" % (t.day, t.month)


def duracao(td):
    m = int(round(td.total_seconds() / 60.))
    return "%d~h %02d~min" % (m // 60, m % 60) if m >= 60 else "%d~min" % m


def rel(arq, base):
    return os.path.relpath(arq, base)


# ---------------------------------------------------------------- leitura
def texto(arq):
    with open(arq, errors="replace") as f:
        return f.read()


def le_fonte(arq, padrao, nome, flags=re.S | re.M):
    m = re.search(padrao, texto(arq), flags)
    if not m:
        sys.exit("nao achei %s em %s" % (nome, arq))
    return m


def le_chaves(arq):
    return dict(l.strip().split("=", 1) for l in open(arq) if "=" in l and not l.startswith("#"))


def le_casos(arq):
    """as linhas de um casos*.txt: id -> (N, opcoes)"""
    out = {}
    for linha in texto(arq).splitlines():
        if linha.strip() and not linha.startswith("#"):
            c = linha.split()
            out[c[0]] = (int(c[2]), dict(o.split("=", 1) for o in c[3:] if "=" in o))
    return out


def le_linhas(pasta):
    """as linhas dos mc_<a>_<b>.csv por realizacao e instante (a -1 e o deterministico), so as completas"""
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
    return per


RX_PAROU = re.compile(r"realizacao (\d+): o transiente parou em ([\d.]+) h")


def le_solo(pasta, N):
    """um solo: instantes, deterministico, realizacoes completas, as que pararam e as validas (todas ok), as que faltam,
    a matriz F (validas x instantes) e o instante em que o transiente parou (logs)"""
    par = le_chaves(os.path.join(pasta, "mc_parametros.txt"))
    ts = [float(x) for x in par["TEMPOS"].split(",")]
    per = le_linhas(pasta)
    det = per.pop(-1)
    completas = {n: d for n, d in per.items() if all(t in d for t in ts)}
    ruins = sorted(n for n, d in completas.items()
                   if any(d[t]["situacao"] != "ok" or not float(d[t]["FS_min"]) > 0 for t in ts))
    validas = {n: completas[n] for n in sorted(completas) if n not in ruins}
    faltam = sorted(set(range(N)) - set(completas))
    parou = {}
    for arq in glob.glob(os.path.join(pasta, "mc_*_*.log")):
        for m in RX_PAROU.finditer(texto(arq)):
            parou[int(m[1])] = float(m[2])
    F = np.array([[float(validas[n][t]["FS_min"]) for t in ts] for n in validas])
    return {"par": par, "ts": ts, "det": det, "completas": completas, "ruins": ruins, "validas": validas,
            "faltam": faltam, "parou": parou, "F": F, "reais": list(validas)}


RX_LOG = re.compile(r"t =\s*([\d.]+) h: FS = ([\d.]+) \|")  # o RX_LOG do figuras_relatorio.py, so o t e o FS


def le_log(arq):
    """FS(t) de hora em hora de um log do SlopeRahardjo e a intensidade do cabecalho"""
    tx = texto(arq)
    L = {float(m[1]): float(m[2]) for m in RX_LOG.finditer(tx)}
    m = re.search(r"I_r = ([\d.]+) mm/h", tx)
    if not L or not m:
        sys.exit("sem linhas de FS ou sem I_r: %s" % arq)
    return L, float(m[1])


def le_mapa(pasta, t):
    """celulas de 1 m no instante t: (x, y) -> [realizacoes com a banda ali, quantas delas falham] e o numero de
    realizacoes"""
    cel, N = {}, 0
    for arq in glob.glob(os.path.join(pasta, "mc_mapa_*_*.csv")):
        with open(arq, newline="") as f:
            linhas = [r for r in csv.DictReader(f) if abs(float(r["t_h"]) - t) < 1e-9]
        if not linhas:
            continue
        N += int(linhas[0]["n_real"])
        for r in linhas:
            c = cel.setdefault((float(r["x"]), float(r["y"])), [0, 0])
            c[0] += int(r["n"])
            c[1] += int(r["n_falha"])
    return cel, N


RX_TAREFA = re.compile(r"^(\d{4}-\d\d-\d\d \d\d:\d\d:\d\d) (\S+) (\d+)\.\.(\d+): (inicio|fim em (\d+) s)", re.M)


def le_tarefas(pasta):
    """as tarefas dos logs do cluster/tarefas.sh: (solo, a, b) -> [inicio, fim, duracao em s]"""
    tar = {}
    for arq in glob.glob(os.path.join(pasta, "*.txt")):
        for m in RX_TAREFA.finditer(texto(arq)):
            k = (m[2], int(m[3]), int(m[4]))
            t = dt.datetime.strptime(m[1], "%Y-%m-%d %H:%M:%S")
            x = tar.setdefault(k, [None, None, None])
            if m[5] == "inicio":
                x[0] = t
            else:
                x[1], x[2] = t, int(m[6])
    return tar


# ---------------------------------------------------------------- desenho
def letra(ax, t):
    ax.text(-0.02, 1.02, "(%s)" % t, transform=ax.transAxes, ha="right", va="bottom", fontsize=9, color=TINTA)


def faixa_chuva(ax, fim):
    ax.axvspan(0, fim, color=SOLO, lw=0, zorder=0)
    ax.text(fim / 2, 0.985, "chuva", transform=ax.get_xaxis_transform(), ha="center", va="top", fontsize=7,
            color=TINTA2)


def linha_um(ax, eixo="y"):
    (ax.axhline if eixo == "y" else ax.axvline)(1., color=MUDO, lw=0.8, zorder=1)


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


def fig_fs(pasta, S, logs, fig5, chuva, kr):
    """FS(t): media e faixa de 5% a 95% do Monte Carlo, o deterministico na mesma malha, o rahardjo.cpp de hora em hora
    e as leituras da Fig. 5 do Rahardjo"""
    fig, axs = plt.subplots(1, 3, figsize=(LARGURA, 2.85))
    for k, (ax, s) in enumerate(zip(axs, SOLOS)):
        d = S[s]
        t = np.array(d["ts"])
        faixa_chuva(ax, chuva)
        linha_um(ax)
        q05, q95 = np.percentile(d["F"], [5, 95], axis=0)
        ax.fill_between(t, q05, q95, color=COR[s], alpha=0.2, lw=0, zorder=2)
        L, I = logs[s]
        ax.plot(list(L), list(L.values()), color=TINTA2, lw=0.8, zorder=3)
        ax.plot(t, d["F"].mean(axis=0), "o-", color=COR[s], ms=3.2, lw=1.4, zorder=4)
        ax.plot(t, d["Dt"], "s", ms=3.6, mfc="white", mec=TINTA, mew=0.8, zorder=5)
        R = fig5[s]
        ax.plot(list(R), list(R.values()), "x", ms=3.6, color=TINTA, mew=0.9, zorder=6)
        ax.set_title("%s, %s mm/h" % (ROTULO[s], g(I)), fontsize=8, color=TINTA)
        ax.set_xlim(0, 48)
        ax.set_xticks(np.arange(0, 49, 12))
        ax.set_xlabel("$t$ (h)")
        letra(ax, "abc"[k])
    axs[0].set_ylabel("$FS$")
    axs[0].set_ylim(0.6, 2.3)
    axs[1].set_ylim(0.85, 2.35)
    axs[2].set_ylim(1.6, 2.25)
    # a media e a faixa com a cor de cada solo; as outras series sao as mesmas nos tres paineis
    medias = tuple(Line2D([], [], color=COR[s], marker="o", ms=3.2, lw=1.4) for s in SOLOS)
    faixas = tuple(Patch(color=COR[s], alpha=0.2, lw=0) for s in SOLOS)
    leg = [medias, faixas,
           Line2D([], [], ls="", marker="s", ms=3.6, mfc="white", mec=TINTA, mew=0.8),
           Line2D([], [], color=TINTA2, lw=0.8),
           Line2D([], [], ls="", marker="x", ms=3.6, color=TINTA, mew=0.9)]
    rot = ["média do Monte Carlo", "5% a 95% das realizações", "determinístico (mesma malha do Monte Carlo)",
           "rahardjo.cpp, de hora em hora (TriGMesh(%d))" % kr, "Rahardjo et al. (2007), Fig. 5"]
    fig.legend(leg, rot, loc="lower center", ncol=3, fontsize=6.5, handletextpad=0.5, columnspacing=1.6,
               bbox_to_anchor=(0.5, 0.0), handler_map={tuple: HandlerTuple(ndivide=None, pad=0.2)},
               handlelength=3.2)
    fig.tight_layout(w_pad=0.6, rect=(0, 0.12, 1, 1))
    salva(fig, pasta, "fig_solos_fs")


def fig_pf(pasta, S, chuva):
    """(a) P(FS < 1) em cada instante com a faixa de Wilson; (b) o instante do menor FS"""
    fig, (a, b) = plt.subplots(1, 2, figsize=(LARGURA, 2.35), gridspec_kw={"width_ratios": [1.15, 1]})
    faixa_chuva(a, chuva)
    for s in SOLOS:
        d = S[s]
        t = np.array(d["ts"])
        p, lo, hi = (100. * np.array([x[k] for x in d["inst"]]) for k in ("p", "lo", "hi"))
        a.fill_between(t, lo, hi, color=COR[s], alpha=0.2, lw=0, zorder=2)
        a.plot(t, p, "o-", color=COR[s], ms=3.2, zorder=3, label=ROTULO[s])
    a.set_xlim(0, 48)
    a.set_xticks(np.arange(0, 49, 6))
    a.set_ylim(-2, 112)
    a.set_xlabel("$t$ (h)")
    a.set_ylabel("$P(FS < 1)$ (%) e faixa de 95%")
    a.legend(loc="center right", fontsize=6.8)
    letra(a, "a")
    todos = sorted({t for s in SOLOS for t in S[s]["ts"]})
    w = 1.7
    for i, s in enumerate(SOLOS):
        d = S[s]
        c = [100. * d["quando"].get(t, 0) / d["N"] for t in todos]
        b.bar(np.array(todos) + (i - 1) * w, c, width=w, color=COR[s], edgecolor="white", linewidth=0.4,
              label=ROTULO[s])
    b.set_xlim(-3, 51)
    b.set_xticks(todos)
    b.set_ylim(0, 105)
    b.set_xlabel("instante do menor $FS$ (h)")
    b.set_ylabel("realizações (%)")
    b.grid(axis="x", visible=False)
    b.legend(loc="upper right", fontsize=6.8)
    letra(b, "b")
    fig.tight_layout(w_pad=1.2)
    salva(fig, pasta, "fig_solos_pf")


def fig_pareado(pasta, P, S):
    """(a) e (b) o menor FS de cada realizacao em dois solos; (c) o primeiro instante com FS < 1 no f10"""
    fig, (a, b, c) = plt.subplots(1, 3, figsize=(LARGURA, 2.4))
    m10, m50, m100 = P["m"]["f10"], P["m"]["f50"], P["m"]["f100"]
    f10, f50 = m10 < 1., m50 < 1.
    for ax, x, y, gr, rx, ry in ((a, m10, m50, [(~f10, MUDO, "não falha"), (f10 & ~f50, COR["f10"], "falha só no f10"),
                                              (f10 & f50, COR["f50"], "falha nos dois")], "f10", "f50"),
                                 (b, m50, m100, [(~f50, MUDO, "não falha no f50"), (f50, COR["f50"], "falha no f50")],
                                  "f50", "f100")):
        for msk, cor, rot in gr:
            ax.scatter(x[msk], y[msk], s=3.5, color=cor, lw=0, alpha=0.75, label="%s (%d)" % (rot, msk.sum()))
        lim = [min(x.min(), y.min()) - 0.03, max(x.max(), y.max()) + 0.03]
        ax.plot(lim, lim, color=TINTA2, lw=0.8)
        ax.axvline(1., color=MUDO, lw=0.6)
        ax.axhline(1., color=MUDO, lw=0.6)
        ax.set_xlim(lim)
        ax.set_ylim(lim)
        ax.set_xlabel("menor $FS$ no %s" % ROTULO[rx])
        ax.set_ylabel("menor $FS$ no %s" % ROTULO[ry])
        ax.legend(loc="lower right", fontsize=6.2, handletextpad=0.1, markerscale=2.)
    letra(a, "a")
    letra(b, "b")
    d = S["f10"]
    cats = [t for t in d["ts"] if d["primeiro"].get(t, 0)] + [None]
    vals = [100. * d["primeiro"].get(t, 0) / d["N"] for t in cats]
    xs = np.arange(len(cats))
    c.bar(xs, vals, width=0.65, color=COR["f10"])
    for x, y in zip(xs, vals):
        c.text(x, y + 1.5, v(y, 1), ha="center", va="bottom", fontsize=6.5, color=TINTA)
    c.set_xticks(xs)
    c.set_xticklabels(["%g" % t if t is not None else "nunca" for t in cats])
    c.set_ylim(0, 1.18 * max(vals))
    c.set_xlabel("primeiro $t$ com $FS < 1$ (h)")
    c.set_title(ROTULO["f10"], fontsize=8, color=TINTA)
    c.set_ylabel("realizações (%)")
    c.grid(axis="x", visible=False)
    letra(c, "c")
    fig.tight_layout(w_pad=0.8)
    salva(fig, pasta, "fig_solos_pareado")


def fig_mecanismo(pasta, pts, cel, N, nf, banda_det, rhos, t):
    """o mapa da banda no instante t (todas e falhas) e o Spearman do FS com as medias na banda"""
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
            ax.plot(banda_det[:, 0], banda_det[:, 1], ".", color=COR["f50"], ms=0.5, alpha=0.25, zorder=4)
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
    cb.set_label("fração das realizações com a banda na célula de 1 m (%g h)" % t, fontsize=7)
    cb.ax.tick_params(labelsize=6.5)
    cb.ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda x, p_: v(x, 1)))
    c = axs[2]
    nomes = [r for r, _ in rhos][::-1]
    vals = [x for _, x in rhos][::-1]
    c.barh(range(len(vals)), vals, color=COR["f50"], height=0.6)
    for i, x in enumerate(vals):
        c.text(x + (0.03 if x >= 0 else -0.03), i, v(x, 2), va="center", ha="left" if x >= 0 else "right",
               fontsize=6.5, color=TINTA)
    c.set_yticks(range(len(vals)))
    c.set_yticklabels(nomes, fontsize=7)
    c.axvline(0, color=TINTA2, lw=0.8)
    c.set_xlim(-0.7, 1.0)
    c.set_xlabel("Spearman com $FS$(%g h)" % t)
    c.grid(axis="y", visible=False)
    letra(c, "c")
    salva(fig, pasta, "fig_solos_mecanismo")


# ---------------------------------------------------------------- secoes
def modelo(C, S):
    """talude, resistencia, solos, campos, malhas e instantes"""
    src = os.path.join(PROJ, "montecarlo.cpp")
    m = le_fonte(src, r"kTalude = \{(.*?)\};", "kTalude")
    pts = [[float(a), float(b)] for a, b in re.findall(r"\{([-\d.]+),\s*([-\d.]+)\}", m[1])]
    (xb0, yb0), (xb1, _), (_, yd), (xpe, ype), (xcr, ycr), _ = pts
    H = ycr - ype
    num("AlturaH", "%g" % H, "montecarlo.cpp (kTalude)")
    num("Inclinacao", "%g" % round(math.degrees(math.atan2(H, xpe - xcr))), "montecarlo.cpp (kTalude)")
    num("DomX", "%g" % (xb1 - xb0), "montecarlo.cpp (kTalude)")
    num("DomY", "%g" % (ycr - yb0), "montecarlo.cpp (kTalude)")
    corpo = le_fonte(src, r"struct Resistencia \{.*?\n(.*?)\};", "Resistencia")[1]
    res = dict(re.findall(r"(\w+) = ([\d.e+-]+)", corpo))
    num("Coesao", "%g" % float(res["c"]), "montecarlo.cpp (Resistencia)")
    num("Atrito", "%g" % float(res["phiDeg"]), "montecarlo.cpp (Resistencia)")
    num("PesoSolo", "%g" % float(res["gamma"]), "montecarlo.cpp (Resistencia)")
    lim = float(le_fonte(src, r"if \(p\.j >= ([\d.]+) \* jmax\) \{", "limiar da banda")[1])
    num("LimiarBanda", v(lim, 1), "montecarlo.cpp (Mecanismo: banda)")
    casos = le_casos(os.path.join(PROJ, "cluster", "casos_solos.txt"))
    for s in SOLOS:
        par, fp = S[s]["par"], "chuva_solos/%s/mc_parametros.txt" % s
        nm = NOME[s]
        num("A" + nm, "%g" % float(par["A_FX"]), fp)
        e = round(math.log10(float(par["KS"])))
        assert abs(float(par["KS"]) / 10. ** e - 1.) < 1e-9
        num("Ks" + nm, "1e%d" % e, fp, latex="10^{%d}" % e)
        num("I" + nm, g(float(par["I"])), fp)
        num("Razao" + nm, g(float(par["I"]) / (float(par["KS"]) * 3.6e6)), fp + " (I / k_s)")
        num("NPlan" + nm, "%d" % casos[s][0], "cluster/casos_solos.txt")
        num("NInst" + nm, "%d" % len(S[s]["ts"]), fp + " (TEMPOS)")
        txt("TemposFaixa" + nm, "%g a %g" % (S[s]["ts"][0], S[s]["ts"][-1]), fp + " (TEMPOS)")
    p50 = S["f50"]["par"]
    fpar = "chuva_solos/f50/mc_parametros.txt"
    for s in SOLOS[::2]:  # o mesmo em todos os solos, menos solo, a, k_s e I e os instantes do f50
        q = S[s]["par"]
        dif = {k for k in set(q) | set(p50) if q.get(k) != p50.get(k)}
        assert dif <= {"SOLO", "A_FX", "KS", "I", "TEMPOS", "ADAPTA"}, dif
    num("Chuva", "%g" % float(p50["CHUVA"]), fpar)
    num("Hw", "%g" % float(p50["HW"]), fpar)
    for chave, macro in (("COV_C", "CovC"), ("COV_PHI", "CovPhi"), ("COV_K", "CovK")):
        num(macro, v(float(p50[chave]), 1), fpar)
    num("Lx", "%g" % float(p50["LX"]), fpar)
    num("Ly", "%g" % float(p50["LY"]), fpar)
    num("KLM", p50["KL_M"], fpar)
    num("KLvar", v(100. * float(p50["KL_VARIANCIA"]), 2), fpar)
    num("Semente", p50["SEMENTE"], fpar)
    txt("MalhaT", p50["MALHA_TRANSIENTE"], fpar)
    txt("MalhaF", p50["MALHA_FATOR"], fpar)
    num("Nref", p50["NREF"], fpar)
    num("UT", p50["UT"], fpar)
    ts10 = S["f10"]["ts"]
    num("Passo", "%g" % (ts10[1] - ts10[0]), "chuva_solos/f10/mc_parametros.txt (TEMPOS)")
    num("Fim", "%g" % ts10[-1], "chuva_solos/f10/mc_parametros.txt (TEMPOS)")
    # malha do transiente: os triangulos e o lado (a menor aresta) do deterministico (mc_real-1, a malha do transiente)
    dm = os.path.join(C, "chuva_solos", "f50", "mc_real-1")
    V = np.loadtxt(os.path.join(dm, "vertices.csv"), delimiter=",", skiprows=1)
    T = np.loadtxt(os.path.join(dm, "triangles.csv"), delimiter=",", skiprows=1, dtype=int)
    P = V[:, 1:3]
    lados = np.concatenate([np.hypot(*(P[T[:, i]] - P[T[:, j]]).T) for i, j in ((0, 1), (1, 2), (2, 0))])
    lado = lados.min()
    kT = int(p50["MALHA_TRANSIENTE"].replace("tri", ""))
    kF = int(p50["MALHA_FATOR"].replace("tri", ""))
    passo = lado * 2 ** kT
    assert abs(passo - H) < 1e-6, (passo, H)  # a TriGMesh tem passo H
    ftri = "chuva_solos/f50/mc_real-1/triangles.csv"
    num("TriTrans", "%d" % len(T), ftri)
    num("LadoTrans", v(lado, 1), ftri + " (menor aresta)")
    num("LadoFator", v(passo / 2 ** kF, 0), "passo H / 2^%d" % kF)
    num("LadoBanda", v(passo / 2 ** kF / 2 ** int(p50["NREF"]), 2), "passo H / 2^%d / 2^NREF" % kF)
    eq = {S[s]["det"][t]["eq_transiente"] for s in SOLOS for t in S[s]["ts"]}
    assert len(eq) == 1
    num("EqTrans", eq.pop(), "chuva_solos/*/mc_*.csv (eq_transiente)")
    kr = int(le_fonte(os.path.join(PROJ, "rahardjo.cpp"), r"const int refTransiente = (\d+);", "refTransiente")[1])
    num("RefRahardjo", "%d" % kr, "rahardjo.cpp (refTransiente)")
    num("LadoRahardjo", v(passo / 2 ** kr, 2), "passo H / 2^refTransiente")
    num("RefMC", "%d" % kT, fpar + " (MALHA_TRANSIENTE)")
    return pts, H, float(p50["CHUVA"]), kr


def por_solo(S, logs, C):
    """estatisticas por solo e instante, o minimo no tempo e o deterministico"""
    for s in SOLOS:
        d, nm = S[s], NOME[s]
        fonte = "chuva_solos/%s/mc_*.csv" % s
        F, ts = d["F"], d["ts"]
        N = len(F)
        d["N"] = N
        num("N" + nm, "%d" % N, fonte + " (realizacoes com todos os instantes ok)")
        num("Completas" + nm, "%d" % len(d["completas"]), fonte + " (realizacoes com todos os instantes)")
        num("Ruins" + nm, "%d" % len(d["ruins"]), fonte + " (situacao diferente de ok)")
        num("Faltam" + nm, "%d" % len(d["faltam"]), fonte + " (0..N-1 sem todas as linhas)")
        txt("FaltamLista" + nm, lista(d["faltam"]) if d["faltam"] else "nenhuma", fonte)
        Dt = np.array([float(d["det"][t]["FS_min"]) for t in ts])
        d["Dt"] = Dt
        fdet = "chuva_solos/%s/mc_0_4.csv (real -1)" % s
        L, _ = logs[s]
        flog = "rahardjo_%s.log" % s
        d["inst"] = []
        for j, t in enumerate(ts):
            x = F[:, j]
            nf = int((x < 1.).sum())
            lo, hi = wilson(nf, N)
            sf = INST[int(t)]
            m = x.mean()
            d["inst"].append({"t": t, "media": m, "cov": x.std(ddof=1) / m, "p": nf / N, "lo": lo, "hi": hi,
                              "nf": nf})
            num("Media" + nm + sf, v(m), fonte)
            num("Cov" + nm + sf, v(100. * x.std(ddof=1) / m, 1), fonte)
            num("Nf" + nm + sf, "%d" % nf, fonte)
            num("Pf" + nm + sf, v(100. * nf / N, 2), fonte)
            num("Lo" + nm + sf, v(100. * lo, 1), fonte + " (Wilson)")
            num("Hi" + nm + sf, v(100. * hi, 1), fonte + " (Wilson)")
            num("Det" + nm + sf, v(Dt[j]), fdet)
            num("Dif" + nm + sf, vs(100. * (m / Dt[j] - 1.), 1), fonte + " e " + fdet + " (media / det - 1, %)")
            txt("Pos" + nm + sf, acima_abaixo(100. * (m / Dt[j] - 1.)), fonte + " e " + fdet + " (media / det - 1)")
            num("Log" + nm + sf, v(L[t]), flog)
            num("DifLog" + nm + sf, vs(100. * (Dt[j] / L[t] - 1.), 1), fdet + " e " + flog + " (det MC / log - 1, %)")
        difs = [100. * (x["media"] / Dt[j] - 1.) for j, x in enumerate(d["inst"])]
        num("DifFaixaMin" + nm, vs(min(difs), 1), fonte + " e " + fdet + " (menor media / det - 1)")
        num("DifFaixaMax" + nm, vs(max(difs), 1), fonte + " e " + fdet + " (maior media / det - 1)")
        # o minimo no tempo (a falha): o menor FS de cada realizacao e o instante dele (o primeiro, se empata)
        mn = F.min(axis=1)
        im = F.argmin(axis=1)
        d["min"], d["imin"] = mn, im
        nf = int((mn < 1.).sum())
        lo, hi = wilson(nf, N)
        d["nf"], d["pf"], d["lo"], d["hi"] = nf, nf / N, lo, hi
        num("MediaMin" + nm, v(mn.mean()), fonte + " (menor FS de cada realizacao)")
        num("CovMin" + nm, v(100. * mn.std(ddof=1) / mn.mean(), 1), fonte)
        num("NfMin" + nm, "%d" % nf, fonte + " (FS < 1 em algum instante)")
        num("PfMin" + nm, v(100. * nf / N, 2), fonte)
        num("LoMin" + nm, v(100. * lo, 1), fonte + " (Wilson)")
        num("HiMin" + nm, v(100. * hi, 1), fonte + " (Wilson)")
        dmin = Dt.min()
        tdet = [ts[j] for j in range(len(ts)) if Dt[j] <= dmin * (1. + 1e-6)]
        num("DetMin" + nm, v(dmin), fdet + " (menor FS)")
        txt("TDetMin" + nm, (lista(["%g" % t for t in tdet]) + "~h"), fdet + " (instante do menor FS)")
        num("DifMin" + nm, vs(100. * (mn.mean() / dmin - 1.), 1), fonte + " e " + fdet)
        tl = min(L, key=lambda t: (L[t], t))
        num("LogMin" + nm, v(L[tl]), flog + " (menor FS de hora em hora)")
        # todos os instantes empatados no menor FS (o log tem 4 casas): 'de 8 a 24~h', 'em 24~h' ou 'em 32 e 33~h'
        tls = sorted(t for t in L if L[t] <= L[tl] * (1. + 1e-6))
        seguidas = len(tls) > 2 and all(abs(b - a - 1.) < 1e-9 for a, b in zip(tls[:-1], tls[1:]))
        txt("TLogMin" + nm, ("de %g a %g" % (tls[0], tls[-1]) if seguidas else "em " + lista(["%g" % t for t in tls]))
            + "~h", flog + " (instantes do menor FS, com os empates)")
        num("DifLogMin" + nm, vs(100. * (dmin / L[tl] - 1.), 1), fdet + " e " + flog)
        quando = {}
        for j in im:
            quando[ts[j]] = quando.get(ts[j], 0) + 1
        d["quando"] = quando
        tq = sorted(quando)
        txt("QuandoFaixa" + nm, "%g a %g" % (tq[0], tq[-1]) if len(tq) > 1 else "%g" % tq[0],
            fonte + " (instantes em que o menor FS ocorre)")
        for t in ts:
            num("Quando" + nm + INST[int(t)], "%d" % quando.get(t, 0), fonte + " (instante do menor FS)")
            num("QuandoPct" + nm + INST[int(t)], v(100. * quando.get(t, 0) / N, 1), fonte + " (instante do menor FS)")
        # falha em algum instante sem falhar no fim da chuva
        j24 = ts.index(24.)
        fora = int(((mn < 1.) & (F[:, j24] >= 1.)).sum())
        num("FalhaFora" + nm, "%d" % fora, fonte + " (FS < 1 em algum instante e FS(24 h) >= 1)")
        d["fora"] = fora
        txt("FalhaForaFrase" + nm, "nenhuma realização falha só fora de 24~h" if fora == 0 else
            "%d realizaç%s falha%s só fora de 24~h" % (fora, "ão" if fora == 1 else "ões", "" if fora == 1 else "m"),
            fonte + " (FS < 1 em algum instante e FS(24 h) >= 1)")
        num("CovFaixaMin" + nm, v(min(x["cov"] for x in d["inst"]) * 100., 1), fonte + " (menor CoV dos instantes)")
        num("CovFaixaMax" + nm, v(max(x["cov"] for x in d["inst"]) * 100., 1), fonte + " (maior CoV dos instantes)")
        # as que o transiente parou: o instante (logs) e o menor FS antes de parar
        r = d["ruins"]
        txt("RuinsLista" + nm, lista(r) if r else "nenhuma", fonte)
        if r:
            txt("RuinsParou" + nm, lista([v(d["parou"][n], 1) for n in r], "; ") + "~h", "chuva_solos/%s/mc_*.log" % s)
            fsr = [min(float(d["completas"][n][t]["FS_min"]) for t in ts if d["completas"][n][t]["situacao"] == "ok")
                   for n in r]
            txt("RuinsFS" + nm, lista([v(x, 3) for x in fsr], "; "), fonte + " (menor FS antes de parar)")
            falham = sum(x < 1. for x in fsr)
            num("RuinsFalham" + nm, "%d" % falham, fonte + " (FS < 1 antes de parar)")
            nf2 = nf + falham
            pf1, pf2 = v(100. * nf / N, 2), v(100. * nf2 / (N + len(r)), 2)
            num("PfComRuins" + nm, pf2, fonte + " (as que pararam com FS < 1 como falhas)")
            num("NComRuins" + nm, "%d" % (N + len(r)), fonte)
            num("NfComRuins" + nm, "%d" % nf2, fonte)
            # a frase do texto: se elas ja tinham FS < 1 e o que acontece com a P_f quando contam como falhas
            fs_r = lista([v(x, 3) for x in fsr], "; ")
            if falham == len(r):
                frase = "Elas saem da conta, mas já tinham $FS < 1$ antes de parar (%s): contadas como falhas" % fs_r
            elif falham:
                frase = ("Elas saem da conta, mas %d delas já tinham $FS < 1$ antes de parar (%s): contando essas como "
                         "falhas" % (falham, fs_r))
            else:
                frase = "Elas saem da conta; nenhuma tinha $FS < 1$ antes de parar (%s)" % fs_r
            if falham:
                pct = lambda x: x.replace(",", "{,}") + r"\%"  # noqa: E731
                frase += ", a $P_f$ do %s %s (%d de %d)" % (
                    s, "continua " + pct(pf2) if pf1 == pf2 else "passa de %s para %s" % (pct(pf1), pct(pf2)), nf2,
                    N + len(r))
            txt("RuinsFrase" + nm, frase + ".", fonte + " e chuva_solos/%s/mc_*.log (as que pararam)" % s)
        # custo por realizacao: o transiente (s_transiente acumula: o do ultimo instante) e o fator (soma dos instantes)
        tr = np.array([max(float(d["validas"][n][t]["s_transiente"]) for t in ts) for n in d["reais"]])
        fa = np.array([sum(float(d["validas"][n][t]["s_fator"]) for t in ts) for n in d["reais"]])
        d["custo"] = (tr + fa).mean()
        num("CustoTrans" + nm, v(tr.mean(), 0), fonte + " (s_transiente no ultimo instante)")
        num("CustoFator" + nm, v(fa.mean(), 0), fonte + " (soma de s_fator)")
        num("CustoFatorInst" + nm, v(fa.mean() / len(ts), 0), fonte + " (s_fator por instante)")
        num("CustoTotal" + nm, v((tr + fa).mean() / 60., 1), fonte + " (s_transiente + s_fator, min)")
    fd, fc = S["f10"]["fora"], S["f50"]["fora"]
    txt("FalhaForaDois", "nenhuma realização falha só fora de 24~h" if fd == fc == 0 else
        "falham só fora de 24~h %d realizações do f10 e %d do f50" % (fd, fc),
        "chuva_solos/f10 e f50/mc_*.csv (FS < 1 em algum instante e FS(24 h) >= 1)")
    # f10: o patamar de 12 a 24 h (o talude saturado)
    d = S["f10"]
    ts = d["ts"]
    j12, j18, j24 = ts.index(12.), ts.index(18.), ts.index(24.)
    F = d["F"]
    pl = F[:, [j12, j18, j24]]
    num("PatamarIgual", "%d" % int((np.abs(F[:, j18] / F[:, j24] - 1.) < 1e-6).sum()),
        "chuva_solos/f10/mc_*.csv (FS(18 h) = FS(24 h) a 1e-6)")
    num("PatamarMil", "%d" % int((pl.max(axis=1) / pl.min(axis=1) - 1. < 1e-3).sum()),
        "chuva_solos/f10/mc_*.csv (FS de 12, 18 e 24 h a 0,1%)")
    # f100: o minimo depois da chuva e quanto abaixo de FS(24 h)
    d = S["f100"]
    ts = d["ts"]
    j24 = ts.index(24.)
    depois = int((np.array([ts[j] for j in d["imin"]]) > 24.).sum())
    num("DepoisCem", "%d" % depois, "chuva_solos/f100/mc_*.csv (menor FS depois de 24 h)")
    num("DepoisPctCem", v(100. * depois / d["N"], 1), "chuva_solos/f100/mc_*.csv")
    r = 1. - d["min"] / d["F"][:, j24]
    num("AbaixoMediaCem", v(100. * r.mean(), 2), "chuva_solos/f100/mc_*.csv (1 - menor FS / FS(24 h), media, %)")
    num("AbaixoMaxCem", v(100. * r.max(), 1), "chuva_solos/f100/mc_*.csv (o mesmo, maior, %)")


LARG_BURACO = 0.0025  # largura das faixas de FS em volta de 1 (o buraco logo abaixo de 1)


def buraco(S, t=24.):
    """o buraco logo abaixo de FS = 1 no f50 em 24 h. O UnsatAnalysis::SRMArcLength aplica o peso com F = 1 e, se o
    Newton converge, devolve o maior F do caminho, que comeca em 1: um FS um pouco menor que 1 parece sair como FS >= 1.
    A falta logo abaixo de 1 e estimada pela media das duas faixas vizinhas, de mesma largura"""
    d = S["f50"]
    G = d["F"][:, d["ts"].index(t)]
    w = LARG_BURACO
    conta = lambda a, b: int(((G >= a) & (G < b)).sum())  # noqa: E731
    ab, ac = conta(1. - w, 1.), conta(1., 1. + w)
    vab, vac = conta(1. - 2. * w, 1. - w), conta(1. + w, 1. + 2. * w)
    fonte = "chuva_solos/f50/mc_*.csv (FS(%g h) em " % t
    lim = lambda x: v(x, 4).rstrip("0").rstrip(",").replace(",", "{,}")  # noqa: E731
    txt("BuracoFaixaAbaixo", r"$[%s;\, 1)$" % lim(1. - w), "LARG_BURACO do figuras_solos.py")
    txt("BuracoFaixaAcima", r"$[1;\, %s)$" % lim(1. + w), "LARG_BURACO do figuras_solos.py")
    num("BuracoAbaixo", "%d" % ab, fonte + "[%g, 1))" % (1. - w))
    num("BuracoAcima", "%d" % ac, fonte + "[1, %g))" % (1. + w))
    num("BuracoVizAbaixo", "%d" % vab, fonte + "[%g, %g))" % (1. - 2. * w, 1. - w))
    num("BuracoVizAcima", "%d" % vac, fonte + "[%g, %g))" % (1. + w, 1. + 2. * w))
    falta = 0.5 * (vab + vac) - ab
    assert falta > 0., "o texto diz que faltam realizacoes logo abaixo de 1"
    num("BuracoFalta", "%d" % round(falta), "media das faixas vizinhas - [%g, 1)" % (1. - w))
    num("BuracoPontos", v(100. * falta / len(G), 1), "o mesmo / N, em pontos percentuais")
    exatos = [n for n, x in zip(d["reais"], G) if abs(x - 1.) < 1e-6]
    assert exatos, "o texto cita a realizacao com FS = 1 exato"
    txt("BuracoExato", lista(exatos), fonte + "1 a 1e-6: a realizacao)")
    print("buraco (f50, %g h): [%g, 1): %d, [1, %g): %d, vizinhas %d e %d" % (t, 1. - w, ab, 1. + w, ac, vab, vac))


def pareado(S):
    """a mesma realizacao n nos tres solos"""
    com = sorted(set.intersection(*(set(S[s]["reais"]) for s in SOLOS)))
    fonte = "chuva_solos/{f10,f50,f100}/mc_*.csv"
    num("NComuns", "%d" % len(com), fonte + " (validas nos tres solos)")
    m = {}
    for s in SOLOS:
        idx = {n: i for i, n in enumerate(S[s]["reais"])}
        m[s] = S[s]["min"][[idx[n] for n in com]]
    ordem = int(((m["f10"] < m["f50"]) & (m["f50"] < m["f100"])).sum())
    num("Ordem", "%d" % ordem, fonte + " (menor FS: f10 < f50 < f100)")
    # a ordem em cada instante comum aos tres solos
    tc = sorted(set(S["f10"]["ts"]) & set(S["f50"]["ts"]) & set(S["f100"]["ts"]))
    F = {}
    for s in SOLOS:
        idx = {n: i for i, n in enumerate(S[s]["reais"])}
        cols = [S[s]["ts"].index(t) for t in tc]
        F[s] = S[s]["F"][[idx[n] for n in com]][:, cols]
    ok_t = [int(((F["f10"][:, j] < F["f50"][:, j]) & (F["f50"][:, j] < F["f100"][:, j])).sum()) for j in range(len(tc))]
    ate = [tc[j] for j in range(len(tc)) if ok_t[j] == len(com)]
    num("OrdemAte", "%g" % max(t for t in ate if all(tt in ate for tt in tc if tt <= t)),
        fonte + " (ultimo instante ate o qual a ordem vale em todas)")
    j48 = tc.index(48.)
    num("InverteQuarentaOito", "%d" % int((F["f10"][:, j48] > F["f50"][:, j48]).sum()),
        fonte + " (FS(48 h): f10 > f50)")
    j30 = tc.index(30.)
    num("InverteTrinta", "%d" % int((F["f10"][:, j30] > F["f50"][:, j30]).sum()), fonte + " (FS(30 h): f10 > f50)")
    for a, b, mac in (("f10", "f50", "RhoDezCinquenta"), ("f50", "f100", "RhoCinquentaCem"),
                      ("f10", "f100", "RhoDezCem")):
        num(mac, v(stats.spearmanr(m[a], m[b])[0], 2), fonte + " (Spearman do menor FS)")
    j24 = tc.index(24.)
    num("RhoDezCinquentaVinteQuatro", v(stats.spearmanr(F["f10"][:, j24], F["f50"][:, j24])[0], 2),
        fonte + " (Spearman de FS(24 h))")
    f10, f50, f100 = m["f10"] < 1., m["f50"] < 1., m["f100"] < 1.
    num("NfComunsDez", "%d" % f10.sum(), fonte)
    num("NfComunsCinquenta", "%d" % f50.sum(), fonte)
    num("NfComunsCem", "%d" % f100.sum(), fonte)
    num("NfCinquentaEDez", "%d" % (f50 & f10).sum(), fonte + " (falhas do f50 que falham no f10)")
    num("NfSoDez", "%d" % (f10 & ~f50).sum(), fonte + " (so o f10 falha)")
    # o texto diz "em todas" e "tambem falham no f10": se os dados novos mudarem isso, o texto tem que mudar
    assert ordem == len(com), "a ordem f10 < f50 < f100 nao vale em todas: reescreva o texto"
    assert int((F["f10"][:, j48] > F["f50"][:, j48]).sum()) == len(com), "FS(48 h): f10 > f50 nao vale em todas"
    assert (f50 & f10).sum() == f50.sum(), "ha falha do f50 que nao falha no f10: reescreva o texto"
    # as falhas do f50 fora das comuns: o que houve no f10
    fora = [n for i, n in enumerate(S["f50"]["reais"]) if S["f50"]["min"][i] < 1. and n not in set(com)]
    d10 = S["f10"]
    partes = []
    parou = [n for n in fora if n in d10["ruins"]]
    if parou:
        fs = [min(float(d10["completas"][n][t]["FS_min"]) for t in d10["ts"]
                  if d10["completas"][n][t]["situacao"] == "ok") for n in parou]
        partes.append("%s (o transiente do f10 parou depois de $FS < 1$: %s)" % (
            lista(parou), lista([v(x, 3) for x in fs], "; ")) if all(x < 1. for x in fs) else
            "%s (o transiente do f10 parou)" % lista(parou))
    falta = [n for n in fora if n in d10["faltam"]]
    if falta:
        partes.append("%s (ainda rodava no f10)" % lista(falta))
    outras = [n for n in fora if n not in parou and n not in falta]
    if outras:
        partes.append("%s (f100 incompleto)" % lista(outras))
    num("NfCinquentaFora", "%d" % len(fora), "chuva_solos/f50/mc_*.csv (falhas do f50 fora das comuns)")
    txt("CinquentaFora", "; ".join(partes) if partes else "nenhuma", "chuva_solos/f10 e f50/mc_*.csv")
    # f10: o primeiro e o ultimo instante com FS < 1
    d = S["f10"]
    ts = d["ts"]
    pri, ult = {}, {}
    for linha in d["F"]:
        a = [ts[j] for j in range(len(ts)) if linha[j] < 1.]
        k1, k2 = (a[0], a[-1]) if a else (None, None)
        pri[k1] = pri.get(k1, 0) + 1
        ult[k2] = ult.get(k2, 0) + 1
    d["primeiro"] = pri
    # o texto lista so estes instantes
    assert set(pri) <= {6., 12., 18., 24., None}, sorted(pri, key=str)
    assert set(ult) <= {24., 30., None}, sorted(ult, key=str)
    assert set(d["quando"]) <= {12., 18., 24.}, sorted(d["quando"])
    f10f = "chuva_solos/f10/mc_*.csv"
    for t in ts:
        num("Primeiro" + INST[int(t)], "%d" % pri.get(t, 0), f10f + " (primeiro instante com FS < 1)")
        num("PrimeiroPct" + INST[int(t)], v(100. * pri.get(t, 0) / d["N"], 1), f10f)
        num("Ultimo" + INST[int(t)], "%d" % ult.get(t, 0), f10f + " (ultimo instante com FS < 1)")
    num("Nunca", "%d" % pri.get(None, 0), f10f + " (FS >= 1 em todos os instantes)")
    num("NuncaPct", v(100. * pri.get(None, 0) / d["N"], 1), f10f)
    nunca = d["min"][d["min"] >= 1.]
    num("NuncaMin", v(nunca.min(), 3), f10f + " (menor FS das que nunca falham)")
    num("NuncaMax", v(nunca.max(), 3), f10f + " (o maior deles)")
    td = [t for t, x in zip(ts, d["Dt"]) if x < 1.]
    num("DetPrimeiroDez", "%g" % td[0], "chuva_solos/f10/mc_0_4.csv (real -1: primeiro instante com FS < 1)")
    return {"com": com, "m": m}


def cruza_um(L):
    """o instante (interpolacao linear) em que o FS de hora em hora cruza 1 descendo (o cruza_um do
    figuras_relatorio.py)"""
    ts = sorted(L)
    for a, b in zip(ts[:-1], ts[1:]):
        if L[a] >= 1. > L[b]:
            return a + (L[a] - 1.) / (L[a] - L[b]) * (b - a)
    return None


def mecanismo(S, C, pts, t=24.):
    """o que leva a falha no f50 em 24 h: Spearman, medias, categorias do nivel d'agua e o mapa da banda"""
    s = "f50"
    d = S[s]
    j = d["ts"].index(t)
    G = d["F"][:, j]
    f = G < 1.
    fonte = "chuva_solos/f50/mc_*.csv (%g h)" % t
    num("NfMec", "%d" % f.sum(), fonte)
    num("NMec", "%d" % len(G), fonte)
    det = d["det"][t]
    fdet = "chuva_solos/f50/mc_0_4.csv (real -1, %g h)" % t
    col = lambda c: np.array([float(d["validas"][n][t][c]) for n in d["reais"]])  # noqa: E731
    nomes = [("banda_c", "C", "$c'$ na banda", 2), ("banda_phi_deg", "Phi", "$\\varphi'$ na banda", 1),
             ("banda_k_mult", "K", "$k$ na banda", 3), ("massa_k_mult", "KMassa", "$k$ na massa", 3),
             ("banda_prof", "Prof", "profundidade", 2), ("banda_pw", "Pw", "$p_w$ na banda", 1)]
    rhos = []
    for c, mac, rot, casas in nomes:
        x = col(c)
        rho = stats.spearmanr(x, G)[0]
        rhos.append((rot, rho))
        num("Rho" + mac, v(rho, 2), fonte + " (Spearman de FS com %s)" % c)
        num("Falha" + mac, v(x[f].mean(), casas), fonte + " (media de %s nas falhas)" % c)
        num("Outras" + mac, v(x[~f].mean(), casas), fonte + " (media de %s nas outras)" % c)
        num("Det" + mac, v(float(det[c]), casas), fdet + " (%s)" % c)
    X = np.column_stack([G] + [col(c) for c in ("banda_c", "banda_phi_deg", "banda_k_mult")])
    Pm = np.linalg.inv(stats.spearmanr(X)[0])
    for i, mac in ((1, "C"), (2, "Phi"), (3, "K")):
        num("Parc" + mac, v(-Pm[0, i] / math.sqrt(Pm[0, 0] * Pm[i, i]), 2),
            fonte + " (Spearman parcial de 2a ordem, Eqs. 70-71 do Ceron)")
    num("RhoKPw", v(stats.spearmanr(col("banda_k_mult"), col("banda_pw"))[0], 2),
        fonte + " (Spearman de banda_k_mult com banda_pw)")
    # categorias do nivel d'agua (Tabela 3 do Rahardjo; as regras do posproc/dados_chuva.py): I chega a superficie na
    # crista, II no meio da face, III no pe, IV em nenhuma das tres
    tol = float(le_fonte(os.path.join(PROJ, "posproc", "dados_chuva.py"), r"^TOL_NIVEL = ([\d.]+)", "TOL_NIVEL")[1])
    num("TolNivel", v(tol, 2), "posproc/dados_chuva.py (TOL_NIVEL)")
    (_, _), (_, _), (_, _), (xpe, ype), (xcr, ycr), _ = pts
    xm = 0.5 * (xpe + xcr)
    vert = [(xcr, ycr), (xm, 0.5 * (ype + ycr)), (xpe, ype)]
    for x, _ in vert:
        assert ("nivel_x%g" % x) in det, x

    def categoria(r):
        for (x, cota), c in zip(vert, ("I", "II", "III")):
            if float(r["nivel_x%g" % x]) >= cota - tol:
                return c
        return "IV"
    cat = np.array([categoria(d["validas"][n][t]) for n in d["reais"]])
    for c, mac in (("I", "Um"), ("II", "Dois"), ("III", "Tres"), ("IV", "Quatro")):
        k = cat == c
        num("Cat" + mac, "%d" % k.sum(), fonte + " (categoria %s do nivel d'agua)" % c)
        num("CatNf" + mac, "%d" % (k & f).sum(), fonte + " (falhas na categoria %s)" % c)
        num("CatPf" + mac, v(100. * (k & f).sum() / k.sum(), 1) if k.any() else "--", fonte + " (P_f na categoria %s)"
            % c)
    txt("CatDet", categoria(det), fdet + " (categoria do deterministico)")
    xm_col = "nivel_x%g" % xm
    num("NivelMeio", v(col(xm_col).mean(), 1), fonte + " (media de %s)" % xm_col)
    num("NivelMeioFalhas", v(col(xm_col)[f].mean(), 1), fonte + " (%s nas falhas)" % xm_col)
    num("NivelMeioOutras", v(col(xm_col)[~f].mean(), 1), fonte + " (%s nas outras)" % xm_col)
    num("NivelMeioDet", v(float(det[xm_col]), 1), fdet + " (%s)" % xm_col)
    num("CotaMeio", "%g" % (0.5 * (ype + ycr)), "montecarlo.cpp (kTalude: meio da face)")
    num("XMeio", "%g" % xm, "montecarlo.cpp (kTalude: meio da face)")
    num("RhoNivelMeio", v(stats.spearmanr(col(xm_col), G)[0], 2), fonte + " (Spearman de FS com %s)" % xm_col)
    num("RhoNivelK", v(stats.spearmanr(col(xm_col), col("k_mult_med"))[0], 2),
        fonte + " (Spearman de %s com k_mult_med)" % xm_col)
    # o mapa da banda e a banda do deterministico no mesmo instante
    pasta = os.path.join(C, "chuva_solos", s)
    cel, Nm = le_mapa(pasta, t)
    assert Nm == len(G), "o mapa tem %d realizacoes" % Nm
    mec = np.loadtxt(os.path.join(pasta, "mc_real-1", "mecanismo.csv"), delimiter=",", skiprows=1)
    mec = mec[np.abs(mec[:, 0] - t) < 1e-9]
    lim = float(le_fonte(os.path.join(PROJ, "montecarlo.cpp"), r"if \(p\.j >= ([\d.]+) \* jmax\) \{", "limiar")[1])
    banda_det = mec[mec[:, 3] >= lim * mec[:, 3].max()][:, 1:3]
    fm = "chuva_solos/f50/mc_mapa_*.csv (%g h)" % t
    num("CelulasMapa", "%d" % len(cel), fm)
    num("CelulasMetade", "%d" % sum(1 for c in cel.values() if c[0] >= 0.5 * Nm), fm + " (banda em >= 50% de todas)")
    num("CelulasMetadeFalha", "%d" % sum(1 for c in cel.values() if c[1] >= 0.5 * f.sum()),
        fm + " (banda em >= 50% das falhas)")
    num("ProfDetMec", v(float(det["banda_prof"]), 1), fdet + " (banda_prof)")
    return cel, Nm, int(f.sum()), banda_det, rhos


def contra_det(S):
    """o Monte Carlo contra o deterministico: a queda do FS no f50 e o nivel d'agua"""
    d = S["f50"]
    ts = d["ts"]
    j0, j24 = ts.index(0.), ts.index(24.)
    q = d["F"][:, j0] - d["F"][:, j24]
    fonte = "chuva_solos/f50/mc_*.csv"
    fdet = "chuva_solos/f50/mc_0_4.csv (real -1)"
    num("QuedaMedia", v(q.mean(), 3), fonte + " (FS(0 h) - FS(24 h), media)")
    num("QuedaDet", v(d["Dt"][j0] - d["Dt"][j24], 3), fdet + " (FS(0 h) - FS(24 h))")
    num("QuedaMenor", "%d" % int((q < d["Dt"][j0] - d["Dt"][j24]).sum()), fonte + " (queda menor que a do det)")
    for s in SOLOS:  # o maior |det do MC / log - 1| de cada solo e onde
        d = S[s]
        L = d["log"]
        dif = [(abs(d["Dt"][j] / L[t] - 1.), t) for j, t in enumerate(d["ts"])]
        x, t = max(dif)
        num("DifLogMax" + NOME[s], v(100. * x, 1), "chuva_solos/%s/mc_0_4.csv e rahardjo_%s.log" % (s, s))
        num("DifLogMaxT" + NOME[s], "%g" % t, "chuva_solos/%s/mc_0_4.csv e rahardjo_%s.log" % (s, s))
    # o nivel d'agua medio contra o do deterministico no f10: em 6 h no meio da face e em 48 h na crista
    d = S["f10"]
    for t, x, mac in ((6., 35, "NivelMeioDezSeis"), (48., 30, "NivelCristaDezQuarentaOito")):
        c = "nivel_x%d" % x
        num(mac, v(np.mean([float(d["validas"][n][t][c]) for n in d["reais"]]), 1),
            "chuva_solos/f10/mc_*.csv (%s em %g h, media)" % (c, t))
        num(mac.replace("Nivel", "NivelDet", 1), v(float(d["det"][t][c]), 1),
            "chuva_solos/f10/mc_0_4.csv (real -1, %s em %g h)" % (c, t))
    # sem o f10 em 6 h (a frente de saturacao em 360 mm/h): o maior dos outros
    resto = [abs(S[s]["Dt"][j] / S[s]["log"][t] - 1.) for s in SOLOS for j, t in enumerate(S[s]["ts"])
             if not (s == "f10" and t == 6.)]
    num("DifLogResto", v(100. * max(resto), 1), "chuva_solos/*/mc_0_4.csv e rahardjo_*.log (sem o f10 em 6 h)")


def contexto(S, B):
    """a referencia da chuva (I = 31,05 mm/h): a campanha chuva_ai, os pilotos e a sensibilidade a I"""
    readme = texto(os.path.join(PROJ, "cluster", "README.md"))
    pfc = float(re.search(r"tem P_f ≈ (\d+)% em 24 h, pelos pilotos", readme)[1])
    num("PfChuvaPiloto", "%g" % pfc, "cluster/README.md (pilotos do relatorio_chuva.tex)")
    iref = float(re.search(r"^chuva_ai\).*?IREF_PADRAO=([\d.]+)", texto(os.path.join(PROJ, "cluster", "campanha.sh")),
                           re.M | re.S)[1])
    num("Iref", v(iref, 2), "cluster/campanha.sh (chuva_ai)")
    ca = le_casos(os.path.join(PROJ, "cluster", "casos_chuva_ai.txt"))
    tp = ca["ref48_ai"][1]["tempos"].split(",")
    txt("ChuvaAITempos", "%s a %s" % (tp[0], tp[-1]), "cluster/casos_chuva_ai.txt")
    num("ChuvaAIPiloto", "%d" % ca["ref48_piloto"][0], "cluster/casos_chuva_ai.txt")
    num("ChuvaAIN", "%d" % ca["ref48_ai"][0], "cluster/casos_chuva_ai.txt")
    # o f50 com I = 31,05 mm/h: o deterministico e o piloto de 8 realizacoes (as realizacoes 0..7 daqui: mesma semente)
    gi = os.path.join(B, "mc", "gim_I")
    fdet = os.path.join(gi, "srm_I%g_0_-1.csv" % iref)
    with open(fdet, newline="") as f:
        det31 = [float(r["FS_min"]) for r in csv.DictReader(f) if int(r["real"]) == -1 and float(r["t_h"]) == 24.]
    assert len(det31) == 1
    det31 = det31[0]
    d = S["f50"]
    j24 = d["ts"].index(24.)
    det36 = d["Dt"][j24]
    num("DetIref", v(det31), rel(fdet, B) + " (real -1, 24 h)")
    num("RazaoDet", v(det31 / det36, 3), rel(fdet, B) + " e chuva_solos/f50/mc_0_4.csv (det 31,05 / det 36)")
    num("SobeDet", v(100. * (det31 / det36 - 1.), 1), "o mesmo, em %")
    fpil = os.path.join(gi, "piloto_srm_I%g_0_7.csv" % iref)
    ppar = le_chaves(os.path.join(gi, "piloto_srm_I%g_parametros.txt" % iref))
    dif = {k for k in set(ppar) | set(d["par"]) if ppar.get(k) != d["par"].get(k)}
    assert dif <= {"I", "TEMPOS", "ADAPTA", "AI"}, dif  # o mesmo caso, so com outra intensidade
    with open(fpil, newline="") as f:
        pil = {int(r["real"]): float(r["FS_min"]) for r in csv.DictReader(f)
               if int(r["real"]) >= 0 and float(r["t_h"]) == 24. and r.get("situacao", "ok") == "ok"}
    idx = {n: i for i, n in enumerate(d["reais"])}
    ns = sorted(n for n in pil if n in idx)
    fp = rel(fpil, B)
    a36 = np.array([d["F"][idx[n], j24] for n in ns])
    a31 = np.array([pil[n] for n in ns])
    raz = a31 / a36
    num("NPiloto", "%d" % len(ns), fp + " (as mesmas realizacoes daqui)")
    txt("PilotoReais", "%d a %d" % (ns[0], ns[-1]) if ns == list(range(ns[0], ns[-1] + 1)) else lista(ns), fp)
    num("PilotoCov", v(100. * a31.std(ddof=1) / a31.mean(), 1), fp + " (CoV de FS(24 h) em 31,05 mm/h)")
    num("PilotoCovAqui", v(100. * a36.std(ddof=1) / a36.mean(), 1),
        "chuva_solos/f50/mc_*.csv (CoV de FS(24 h) dessas realizacoes em 36 mm/h)")
    num("SobeMin", v(100. * (raz.min() - 1.), 1), fp + " e chuva_solos/f50 (FS31 / FS36 - 1, menor, %)")
    num("SobeMax", v(100. * (raz.max() - 1.), 1), fp + " e chuva_solos/f50 (o maior, %)")
    num("SobeMedia", v(100. * (raz.mean() - 1.), 1), fp + " e chuva_solos/f50 (a media, %)")
    G = d["F"][:, j24]
    fonte = "chuva_solos/f50/mc_*.csv (FS(24 h) x razao < 1)"
    # extrapolacoes pontuais (a mesma subida em todas as realizacoes), nao limites: a subida media e a do
    # deterministico, e a menor e a maior das realizacoes do piloto
    est = []
    for r, mac in ((raz.mean(), "PfEstMedia"), (det31 / det36, "PfEstDet"), (raz.min(), "PfEstSobeMin"),
                   (raz.max(), "PfEstSobeMax")):
        est.append((G * r < 1.).mean())
        num(mac, v(100. * est[-1], 1), fonte)
    num("PfEstMin", v(100. * min(est), 1), fonte + " (a menor das extrapolacoes)")
    num("PfEstMax", v(100. * max(est), 1), fonte + " (a maior das extrapolacoes)")
    assert min(est) > pfc / 100., "o texto diz que todas as extrapolacoes ficam acima da P_f dos pilotos"
    npc = ca["ref48_piloto"][0]
    num("FalhasPilotoMin", "%d" % round(npc * min(est)), "piloto da chuva_ai x a menor P_f extrapolada")
    num("FalhasPilotoMax", "%d" % round(npc * max(est)), "piloto da chuva_ai x a maior P_f extrapolada")
    # a P_f lognormal dos pilotos (o 'P_f ~ 1%' do README.md): 8 realizacoes, sem falha
    tlog = texto(os.path.join(gi, "piloto_srm.log"))
    assert "solo f50, I = %g mm/h" % iref in tlog and os.path.basename(fpil) in tlog, "piloto_srm.log de outro caso"
    m = re.search(r"P\(FS < 1\)\s+P_f lognormal\s*\n\s*24\s+[\d.]+\s+[\d.]+\s+([\d.]+)\s+([\d.]+)%\s+([\d.]+)%", tlog)
    num("PilotoLognormal", v(float(m[3]), 2), rel(os.path.join(gi, "piloto_srm.log"), B) + " (P_f lognormal)")
    assert float(m[2]) == 0., "o texto diz que o piloto nao teve falha"
    assert len(pil) == len(ns), "o piloto_srm.log e o csv tem as mesmas realizacoes"
    # sensibilidade: o FS 1% maior ou menor em todas as realizacoes
    for x, mac in ((1.01, "PfMaisUm"), (0.99, "PfMenosUm")):
        num(mac, v(100. * (G * x < 1.).mean(), 1), "chuva_solos/f50/mc_*.csv (FS(24 h) x %g < 1)" % x)


def cluster(C, S):
    """os jobs, as tarefas, o relogio e a CPU; devolve os blocos rodando, o total de blocos e a previsao do fim"""
    R = os.path.join(C, "chuva_solos")
    outs = sorted(glob.glob(os.path.join(C, "logs", "solosMC_*.out")))
    jobs = {}
    for arq in outs:
        m = re.search(r"job (\d+) \(solos\) em (\S+): (\d+) processos", texto(arq))
        if m:
            jobs[m[1]] = (m[2], int(m[3]))
    nos = len({x[0] for x in jobs.values()})
    proc = {x[1] for x in jobs.values()}
    assert len(proc) == 1
    proc = proc.pop()
    num("Nos", "%d" % nos, "logs/solosMC_*.out (nos dos jobs)")
    num("ProcNo", "%d" % proc, "logs/solosMC_*.out (processos por no)")
    num("Proc", "%d" % (nos * proc), "nos x processos")
    bloco = int(texto(os.path.join(R, "campanha.txt")).split()[1])
    num("Bloco", "%d" % bloco, "chuva_solos/campanha.txt")
    tar = le_tarefas(os.path.join(R, "logs"))
    ini = min(x[0] for x in tar.values() if x[0])
    feitas = {k: x for k, x in tar.items() if x[1]}
    rodando = {k: x for k, x in tar.items() if not x[1]}
    fim = max(x[1] for x in feitas.values())
    ntar = sum(math.ceil(S[s]["NPlan"] / bloco) for s in SOLOS)
    ft = "chuva_solos/logs/p*.txt"
    txt("Inicio", hora(ini), ft)
    txt("UltimaTarefa", hora(fim), ft + " (fim da ultima tarefa terminada)")
    txt("Relogio", duracao(fim - ini), ft)
    num("RelogioHoras", v((fim - ini).total_seconds() / 3600., 1), ft)
    num("Tarefas", "%d" % ntar, "cluster/casos_solos.txt e campanha.txt (N / BLOCO)")
    num("TarefasFeitas", "%d" % len(feitas), ft)
    num("TarefasRodando", "%d" % len(rodando), ft + " (inicio sem fim)")
    cpu = sum(x[2] for x in feitas.values()) / 3600.
    num("CpuHoras", v(cpu, 0), ft + " (soma das tarefas terminadas)")
    nreal = sum(k[2] - k[1] + 1 for k in feitas)
    num("MinPorReal", v(cpu * 60. / nreal, 1), ft + " (CPU / realizacoes das tarefas terminadas)")
    readme = texto(os.path.join(PROJ, "cluster", "README.md"))
    num("MinReadme", re.search(r"~(\d+) min por realização no fermi; as \d+ realizações", readme)[1],
        "cluster/README.md (estimativa)")
    # quando as tarefas que rodavam devem acabar: inicio + a duracao media das tarefas do mesmo solo
    med = {s: np.mean([x[2] for k, x in feitas.items() if k[0] == s]) for s in SOLOS}
    prev = None
    if rodando:
        prev = hora(max(x[0] + dt.timedelta(seconds=med[k[0]]) for k, x in rodando.items()), 5)
        txt("PrevFim", prev, ft + " (inicio + duracao media das tarefas do solo, estimado)")
    else:
        txt("PrevFim", "--", ft)
    return len(rodando), ntar, prev


def status(S, nrod, ntar, prev):
    """a campanha esta pronta (todas as realizacoes e nenhum bloco rodando)? e a frase do que falta"""
    falt = [(s, S[s]["faltam"]) for s in SOLOS if S[s]["faltam"]]
    nfalt = sum(len(x) for _, x in falt)
    if nrod:
        frase = "Quando os dados foram baixados, %d dos %d blocos ainda %s (%s acabar por volta das %s)." % (
            nrod, ntar, "rodava" if nrod == 1 else "rodavam", "deve" if nrod == 1 else "devem", prev)
    else:
        frase = "Quando os dados foram baixados, nenhum bloco rodava."
    if nfalt:
        frase += " %s %d realizaç%s: %s." % ("Falta" if nfalt == 1 else "Faltam", nfalt, "ão" if nfalt == 1 else "ões",
                                           lista(["%s no %s" % (lista(x), s) for s, x in falt]))
    else:
        frase += " Todas as realizações já estavam nos arquivos."
    txt("FaltamFrase", frase, "chuva_solos/logs/p*.txt e chuva_solos/*/mc_*.csv (blocos rodando e realizacoes que "
        "faltam)")
    return nfalt == 0 and nrod == 0


def main(argv):
    C = os.path.abspath(argv[1]) if len(argv) > 1 else CAMPANHA_PADRAO
    B = os.path.abspath(argv[2]) if len(argv) > 2 else BUILD_PADRAO
    pasta = os.path.join(AQUI, "figs")
    os.makedirs(pasta, exist_ok=True)
    casos = le_casos(os.path.join(PROJ, "cluster", "casos_solos.txt"))
    S = {}
    for s in SOLOS:
        S[s] = le_solo(os.path.join(C, "chuva_solos", s), casos[s][0])
        S[s]["NPlan"] = casos[s][0]
    logs = {}
    for s in SOLOS:
        L, I = le_log(os.path.join(B, "rahardjo_%s.log" % s))
        assert abs(I - float(S[s]["par"]["I"])) < 1e-9, (s, I)
        logs[s] = (L, I)
        S[s]["log"] = L
    src = texto(os.path.join(AQUI, "figuras_relatorio.py"))
    fig5 = ast.literal_eval(re.search(r"^RAHARDJO_FIG5 = (\{.*?^\})", src, re.M | re.S)[1])
    txt("LeituraFig", re.search(r"Fig\. 5, .*?leituras a olho \(\+-([\d,]+)\)", src)[1].replace(",", "{,}"),
        "figuras_relatorio.py (precisao das leituras da Fig. 5)")
    pts, H, chuva, kr = modelo(C, S)
    por_solo(S, logs, C)
    buraco(S)
    P = pareado(S)
    c10 = cruza_um(logs["f10"][0])
    num("CruzaLogDez", v(c10, 1), "rahardjo_f10.log (FS = 1 descendo, interpolado)")
    cel, Nm, nf, banda_det, rhos = mecanismo(S, C, pts)
    contra_det(S)
    contexto(S, B)
    prontos = status(S, *cluster(C, S))
    num("NValidasTotal", "%d" % sum(S[s]["N"] for s in SOLOS), "soma das validas")
    num("NPlanTotal", "%d" % sum(S[s]["NPlan"] for s in SOLOS), "cluster/casos_solos.txt")
    num("FaltamTotal", "%d" % sum(len(S[s]["faltam"]) for s in SOLOS), "soma das que faltam")
    # ---- figuras
    fig_fs(pasta, S, logs, fig5, chuva, kr)
    fig_pf(pasta, S, chuva)
    fig_pareado(pasta, P, S)
    fig_mecanismo(pasta, pts, cel, Nm, nf, banda_det, rhos, 24.)
    # ---- numeros_solos.tex
    arq = os.path.join(AQUI, "numeros_solos.tex")
    with open(arq, "w") as f:
        f.write("%% gerado por figuras_solos.py em %s: nao edite (rode o script de novo)\n"
                % dt.datetime.now().strftime("%d/%m/%Y %H:%M"))
        # a campanha pronta: todas as realizacoes nos arquivos e nenhum bloco rodando
        f.write("\\newif\\ifsolosprontos\\solosprontos%s %% chuva_solos/logs/p*.txt e chuva_solos/*/mc_*.csv\n"
                % ("true" if prontos else "false"))
        for macro, valor, latex, fonte in NUMEROS:
            f.write("\\newcommand{\\%s}{%s} %% %s\n" % (macro, latex, fonte.replace("%", r"\%")))
    print("numeros: %s" % arq)
    print("\n% ---- numeros do texto (macro: valor  [arquivo])")
    print("%-28s %-34s [%s]" % ("ifsolosprontos", "true" if prontos else "false",
                                "chuva_solos/logs/p*.txt e chuva_solos/*/mc_*.csv (campanha pronta)"))
    for macro, valor, latex, fonte in NUMEROS:
        print("%-28s %-34s [%s]" % (macro, valor, fonte))


if __name__ == "__main__":
    main(sys.argv)
