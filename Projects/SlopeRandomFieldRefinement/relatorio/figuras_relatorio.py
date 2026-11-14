#!/usr/bin/env python3
"""figuras_relatorio.py: as figuras e as tabelas do relatorio_chuva.tex (Monte Carlo da chuva no talude de Rahardjo et
al. 2007 com o metodo estocastico de Vargas Ceron et al. 2025). So le os arquivos gravados pelos programas desta pasta
(SlopeRahardjo, SlopeRainMonteCarlo e form_variaveis.py) e o casos.txt da campanha; nao roda nenhum deles.

Figuras (PDF em <figs>):
  fig_geometria    talude da Fig. 2 do Rahardjo, contornos, nivel d'agua inicial e a malha do transiente do Monte Carlo
  fig_estado       f50 a 36 mm/h (o deterministico do rahardjo_f50.log): zona saturada e banda de ruptura do SRM em
                   0, 12, 24 e 48 h
  fig_validacao    FS(t) dos tres solos com I = 1 k_s contra as leituras da Fig. 5 do Rahardjo
  fig_gim_srm      FS do SRM e Gamma do GIM no tempo (f50 e f10) e com o nivel de refinamento (f50, 24 h)
  fig_intensidade  FS(24 h) deterministico com I (SRM, GIM e a sigmoide da Fig. 7b do Rahardjo); pilotos; FS(24 h)
                   com k/k_s no solo homogeneo
Na tela: as linhas LaTeX das tabelas do relatorio (copiadas para o .tex) e os numeros do texto, cada um com o arquivo de
onde sai. Os unicos numeros digitados aqui sao os do artigo do Rahardjo (leituras da Fig. 5 a olho, +-0,02; F_ini da
Tabela 7; a sigmoide e a intensidade limiar da Fig. 7b), os do artigo do Ceron (Gamma deterministico e P_f da Tabela 4) e o fator 2,6 entre um
nucleo do fermi e o da estacao (medido antes).

Uso: python3 -I figuras_relatorio.py [pasta do build do SlopeRandomFieldRefinement] [pasta das figuras]
     (padrao: ../../../../neopz-master-build-release/Projects/SlopeRandomFieldRefinement e ./figs)"""

import csv
import glob
import json
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
from matplotlib.tri import Triangulation  # noqa: E402
from scipy.optimize import brentq  # noqa: E402

AQUI = os.path.dirname(os.path.abspath(__file__))
PROJ = os.path.dirname(AQUI)  # Projects/SlopeRandomFieldRefinement (fontes)
NEOPZ = os.path.dirname(os.path.dirname(PROJ))
BUILD_PADRAO = os.path.join(os.path.dirname(NEOPZ), "neopz-master-build-release", "Projects",
                            "SlopeRandomFieldRefinement")

# ---- o que vem dos artigos (digitado)
# Rahardjo et al. (2007), Fig. 5, talude de 45 graus, I = 1 k_s por 24 h: leituras a olho (+-0,02). Depois da chuva os
# marcadores da Fig. 5 estao em 28, 32, 40 e 48 h nos tres paineis
RAHARDJO_FIG5 = {
    "f10": {0: 2.12, 1: 1.58, 2: 1.22, 3: 1.12, 6: 0.97, 12: 0.81, 24: 0.81, 28: 1.04, 32: 1.18, 40: 1.40, 48: 1.57},
    "f50": {0: 2.11, 6: 1.72, 12: 1.40, 24: 1.03, 28: 1.08, 32: 1.12, 40: 1.18, 48: 1.24},
    "f100": {0: 2.11, 6: 2.05, 12: 2.00, 24: 1.90, 48: 1.89},
}
F_INI_TAB7 = 3.153 - 0.023 * 45.  # Tabela 7: F_ini = 3,153 - 0,023 beta, beta = 45 graus
SIGMOIDE = (0.82, 1.83, 10.7, -13.13)  # Fig. 7b (f50): F_min(I) = a + b / (1 + exp(-(I - c) / d))
I_THR_F50 = 51.  # Fig. 7b (f50): intensidade limiar I_r(thr) = 51 mm/h (= 1,42 k_s)
# Vargas Ceron et al. (2025): Gamma deterministico (Fig. 24) e P_f da Tabela 4 (s = 1 e s = 400), CoV(Gamma) em s = 1
CERON = {"gamma_det": 1.336, "pf_s1": 11.50, "pf_s400": 21.93, "cov_s1": 23.5}
FERMI_LENTO = 2.6  # um nucleo do fermi e 2,6 vezes mais lento que o da estacao (medido antes desta campanha)

# ---- estilo (o do posproc/figuras_chuva.py)
PALETA = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
TINTA, TINTA2, MUDO, GRADE, SOLO = "#0b0b0b", "#52514e", "#898781", "#dddcd8", "#f3f2ee"
AZUL_CLARO = "#cde2fb"
plt.rcParams.update({
    "font.size": 8.5, "axes.labelsize": 8.5, "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "legend.fontsize": 7,
    "axes.edgecolor": TINTA2, "axes.labelcolor": TINTA, "xtick.color": TINTA2, "ytick.color": TINTA2,
    "axes.grid": True, "grid.color": GRADE, "grid.linewidth": 0.5, "lines.linewidth": 1.4, "legend.frameon": False,
    "savefig.dpi": 200, "font.family": "DejaVu Sans", "mathtext.fontset": "dejavusans",
})
LARGURA = 6.5  # polegadas: a largura do texto do relatorio (A4 com margens de 2,2 cm)

NUMEROS = []  # (chave, valor no texto, arquivo)


def num(chave, valor, fonte):
    NUMEROS.append((chave, valor, fonte))
    return valor


def v(x, casas=4):
    """numero com virgula decimal (sem o sinal de -0,0)"""
    t = "%.*f" % (casas, x)
    if t.startswith("-") and float(t) == 0.:
        t = t[1:]
    return t.replace(".", ",")


def g(x):
    return ("%g" % x).replace(".", ",")


def rel(arq, base):
    return os.path.relpath(arq, base)


# ---------------------------------------------------------------- leitura
RX_LOG = re.compile(r"t =\s*([\d.]+) h: FS = ([\d.]+) \| p_w de ([-\d.]+) a ([-\d.]+) kPa, (\d+) nos na face \| "
                    r"fator: (\d+) equacoes, ([\d.]+) s \| niveis: (.*)$")
class Log(list):
    """as linhas de FS de um log do SlopeRahardjo (dicts, em ordem de t) e a intensidade I (mm/h) do cabecalho"""
    I = None


def texto(arq):
    with open(arq) as f:
        return f.read()


def le_log(arq):
    """as linhas 't = H h: FS = ...' do SlopeRahardjo e a intensidade do cabecalho"""
    out = Log()
    for linha in texto(arq).splitlines():
        m = RX_LOG.search(linha)
        if m:
            out.append({"t": float(m[1]), "fs": float(m[2]), "pw_min": float(m[3]), "pw_max": float(m[4]),
                        "face": int(m[5]), "eq": int(m[6]), "s": float(m[7]),
                        "niveis": [float(x) for x in m[8].split()]})
        m = re.search(r"I_r = ([\d.]+) mm/h", linha)
        if m and out.I is None:
            out.I = float(m[1])
    if not out or out.I is None:
        sys.exit("sem linhas de FS ou sem I_r no cabecalho: %s" % arq)
    return out


def le_csv(arq):
    """um csv do SlopeRainMonteCarlo; aceita os dois cabecalhos (o antigo: realizacao, FS; o atual: real, FS_min)"""
    with open(arq, newline="") as f:
        linhas = list(csv.DictReader(f))
    out = []
    for r in linhas:
        out.append({"real": int(r["real"] if "real" in r else r["realizacao"]), "t": float(r["t_h"]),
                    "fs": float(r["FS_min"] if "FS_min" in r else r["FS"]), "situacao": r.get("situacao", "ok"),
                    "niveis": [float(x) for x in r["niveis"].split(";")], "s_fator": float(r["s_fator"]),
                    "s_transiente": float(r["s_transiente"]), "eq_transiente": int(r["eq_transiente"]),
                    "pw_min": float(r["pw_min"]), "nivel": [float(r["nivel_x30"]), float(r["nivel_x35"]),
                                                            float(r["nivel_x40"])]})
    if not out:
        sys.exit("csv vazio: %s" % arq)
    return out


def det(arq, t=24.):
    """a linha deterministica (real = -1) no instante t"""
    for r in le_csv(arq):
        if r["real"] == -1 and abs(r["t"] - t) < 1e-9:
            return r
    sys.exit("sem a deterministica em %g h: %s" % (t, arq))


def sem_mapa(arqs):
    """os csv das realizacoes (o glob *_0_-1.csv pega tambem os <saida>_mapa_0_-1.csv)"""
    return sorted(a for a in arqs if "_mapa_" not in os.path.basename(a))


def estat(fs):
    """as estatisticas do resumo do montecarlo.cpp: media, desvio (n - 1), menor, P(FS < 1) e P_f log-normal"""
    x = np.asarray(fs, float)
    ln = np.log(x)
    dl = ln.std(ddof=1)
    return {"n": len(x), "media": x.mean(), "dp": x.std(ddof=1), "cov": x.std(ddof=1) / x.mean(), "min": x.min(),
            "nf": int((x < 1.).sum()), "pf_log": 0.5 * math.erfc(ln.mean() / dl / math.sqrt(2.))}


def le_fonte(arq, padrao, nome):
    m = re.search(padrao, texto(arq), re.S)
    if not m:
        sys.exit("nao achei %s em %s" % (nome, arq))
    return m


def geometria():
    """os cantos do talude (kTalude do montecarlo.cpp) e o estado inicial (RainData do RainTransient.h)"""
    m = le_fonte(os.path.join(PROJ, "montecarlo.cpp"), r"kTalude = \{(.*?)\};", "kTalude")
    pts = [[float(a), float(b)] for a, b in re.findall(r"\{([-\d.]+),\s*([-\d.]+)\}", m[1])]
    h = os.path.join(PROJ, "RainTransient.h")
    hw = float(le_fonte(h, r"REAL Hw = ([\d.]+)", "Hw")[1])
    beta = float(le_fonte(h, r"REAL betaW = ([\d.]+) \* M_PI", "betaW")[1])
    smax = float(le_fonte(h, r"REAL suctionMax = ([\d.]+)", "suctionMax")[1])
    dt = float(le_fonte(h, r"REAL dt = ([\d.]+);", "dt")[1])
    return pts, hw, beta, smax, dt


def kl_fracao(Dx, Dy, Lx, Ly, ms, n1=3000):
    """variancia retida pelos m maiores termos da KL do nucleo exp(-|dx|/Lx - |dy|/Ly) no retangulo Dx x Dy: as formas
    fechadas do KLSeparavel.h (autovalores 1D 2c / (w^2 + c^2), produtos 2D, fracao = soma / (Dx Dy))"""
    def modos(comp, L):
        a, c, eps, lam, k = comp / 2., 1. / L, 1e-9, [], 1
        while len(lam) < n1:
            w = brentq(lambda x: x * math.tan(x * a) - c, ((k - 1) * math.pi + eps) / a, ((k - .5) * math.pi - eps) / a)
            lam.append(2 * c / (w * w + c * c))
            w = brentq(lambda x: x + c * math.tan(x * a), ((k - .5) * math.pi + eps) / a, (k * math.pi - eps) / a)
            lam.append(2 * c / (w * w + c * c))
            k += 1
        return np.sort(lam)[::-1]
    p = np.sort(np.outer(modos(Dx, Lx)[:2500], modos(Dy, Ly)[:2500]).ravel())[::-1]
    return [p[:m].sum() / (Dx * Dy) for m in ms]


def cruza_um(ts, fs, desce=True):
    """primeiro instante (interpolacao linear) em que FS cruza 1 descendo (desce) ou subindo"""
    for i in range(len(ts) - 1):
        a, b = fs[i], fs[i + 1]
        if (desce and a >= 1. > b) or (not desce and a < 1. <= b):
            return ts[i] + (a - 1.) / (a - b) * (ts[i + 1] - ts[i])
    return None


# ---------------------------------------------------------------- desenho
def letra(ax, t):
    ax.text(-0.02, 1.02, "(%s)" % t, transform=ax.transAxes, ha="right", va="bottom", fontsize=9, color=TINTA)


def faixa_chuva(ax, fim=24.):
    ax.axvspan(0, fim, color=SOLO, lw=0, zorder=0)
    ax.text(fim / 2, 0.985, "chuva", transform=ax.get_xaxis_transform(), ha="center", va="top", fontsize=7,
            color=TINTA2)


def linha_um(ax):
    ax.axhline(1., color=MUDO, lw=0.8, zorder=1)


def contorno(ax, pts, preenche=True):
    p = np.asarray(pts + [pts[0]], float)
    if preenche:
        ax.fill(p[:, 0], p[:, 1], color=SOLO, lw=0, zorder=0)
    ax.plot(p[:, 0], p[:, 1], color=TINTA2, lw=1.0, zorder=4)


class Virgula(matplotlib.ticker.ScalarFormatter):
    """os rotulos do ScalarFormatter com virgula decimal"""
    def __call__(self, x, pos=None):
        return super().__call__(x, pos).replace(".", ",")


def salva(fig, pasta, nome):
    for ax in fig.axes:  # virgula decimal nos eixos lineares (os logaritmicos ja tem os rotulos escritos)
        for eixo, escala in ((ax.xaxis, ax.get_xscale()), (ax.yaxis, ax.get_yscale())):
            if escala == "linear" and isinstance(eixo.get_major_formatter(), matplotlib.ticker.ScalarFormatter):
                eixo.set_major_formatter(Virgula())
    arq = os.path.join(pasta, nome + ".pdf")
    fig.savefig(arq)
    plt.close(fig)
    print("figura: %s" % arq)


def le_malha(d):
    V = np.loadtxt(os.path.join(d, "vertices.csv"), delimiter=",", skiprows=1)
    T = np.loadtxt(os.path.join(d, "triangles.csv"), delimiter=",", skiprows=1, dtype=int)
    return Triangulation(V[:, 1], V[:, 2], T)


def fig_geometria(B, pasta, pts, hw, beta):
    d = os.path.join(B, "mc", "teste_casos", "ref", "mc_real-1")  # a malha do transiente do Monte Carlo (tri2)
    tri = le_malha(d)
    num("malha_tri2_triangulos", "%d" % len(tri.triangles), rel(os.path.join(d, "triangles.csv"), B))
    xToe, yToe = pts[3]
    xs = np.array([pts[0][0], pts[1][0]])
    nivel = yToe - hw + (xToe - xs) * math.tan(math.radians(beta))
    fig, ax = plt.subplots(figsize=(LARGURA, 2.75))
    contorno(ax, pts)
    ax.triplot(tri, color="#c3c2b7", lw=0.35, zorder=1)
    ax.plot(xs, nivel, color=PALETA[0], lw=1.6, zorder=3)
    ax.text(2, yToe - hw + (xToe - 2) * math.tan(math.radians(beta)) - 1.3,
            "nível d'água inicial: %s m abaixo do pé, subindo %s° para dentro do maciço" % (g(hw), g(beta)),
            color=TINTA, fontsize=7.2, va="top", rotation=-beta, rotation_mode="anchor")  # eixos iguais
    # contornos da agua
    est = dict(fontsize=7, color=TINTA)
    ax.text(46.5, 27.4, "chuva $I$ na superfície (ou face de\npercolação, $p_w = 0$, onde a água empoçaria)",
            va="center", **est)
    sx, sy = [pts[5][0], pts[4][0], pts[3][0], pts[2][0]], [pts[5][1], pts[4][1], pts[3][1], pts[2][1]]
    for x0 in (8, 16, 24, 33, 37, 48, 56, 64):  # posicoes das setas (desenho)
        y0 = np.interp(x0, sx, sy)
        ax.annotate("", xy=(x0, y0 + 0.2), xytext=(x0, y0 + 2.6),
                    arrowprops=dict(arrowstyle="-|>", color=PALETA[0], lw=0.8, mutation_scale=6))
    ax.text(-1.0, 6, "$p_w$ hidrostática\nabaixo do nível;\nsem fluxo acima", ha="right", va="center", **est)
    ax.text(71.0, 5, "$p_w$ hidrostática\nabaixo do nível;\nsem fluxo acima", ha="left", va="center", **est)
    ax.text(35, -1.2, "base: sem fluxo, $u_x = u_y = 0$;  laterais: $u_x = 0$", ha="center", va="top", **est)
    # cotas
    ax.plot([pts[4][0], 44.], [pts[4][1]] * 2, color=MUDO, lw=0.5, zorder=2)
    ax.annotate("", xy=(43, pts[3][1]), xytext=(43, pts[4][1]),
                arrowprops=dict(arrowstyle="<->", color=TINTA2, lw=0.6, shrinkA=0, shrinkB=0))
    ax.text(43.6, 23.2, "$H = %s$ m" % g(pts[4][1] - pts[3][1]), ha="left", va="center", **est)
    ax.text(38.6, 21.0, "%s°" % g(math.degrees(math.atan2(pts[4][1] - pts[3][1], pts[3][0] - pts[4][0]))),
            ha="right", va="bottom", **est)
    ax.set_aspect("equal")
    ax.set_xlim(-13, 83)
    ax.set_ylim(-3.5, 31)
    ax.spines["bottom"].set_bounds(0, 70)
    ax.spines["left"].set_bounds(0, 30)
    ax.set_xticks(np.arange(0, 71, 10))
    ax.set_yticks(np.arange(0, 31, 10))
    ax.set_xlabel("$x$ (m)")
    ax.set_ylabel("$y$ (m)")
    ax.grid(False)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    salva(fig, pasta, "fig_geometria")


def fig_estado(B, pasta, pts):
    d = os.path.join(B, "f50", "chuva_f50_36")  # a animacao do deterministico do rahardjo_f50.log
    tri = le_malha(d)
    P = np.loadtxt(os.path.join(d, "pontos.csv"), delimiter=",", skiprows=1)
    snaps = {}
    with open(os.path.join(d, "snaps_pw.csv")) as f:
        for linha in f:
            c = linha.strip().split(",")
            try:
                snaps[float(c[0])] = np.array(c[1:], float)
            except ValueError:
                continue
    S = np.loadtxt(os.path.join(d, "srm_pontos.csv"), delimiter=",", skiprows=1, usecols=(0, 1, 2))
    fs = dict(np.loadtxt(os.path.join(d, "fs.csv"), delimiter=",", skiprows=1))
    tempos = (0., 12., 24., 48.)
    fig, axs = plt.subplots(2, 2, figsize=(LARGURA, 3.55), sharex=True, sharey=True)
    for k, (ax, t) in enumerate(zip(axs.ravel(), tempos)):
        contorno(ax, pts)
        pw = snaps[t]
        if pw.max() > 0:
            ax.tricontourf(tri, pw, levels=[0., pw.max() + 1.], colors=[AZUL_CLARO], zorder=1)
            ax.tricontour(tri, pw, levels=[0.], colors=[PALETA[0]], linewidths=1.0, zorder=2)
        sel = S[np.abs(S[:, 0] - t) < 1e-9]
        j, idx = sel[:, 2], sel[:, 1].astype(int)
        banda = j >= 0.1 * j.max()
        ax.scatter(P[idx[banda], 1], P[idx[banda], 2], s=1.6, color=PALETA[1], lw=0, zorder=3)
        ax.set_title("$t = %s$ h: $FS = %s$" % (g(t), v(fs[t], 4).replace(",", "{,}")), fontsize=8, color=TINTA,
                     loc="left")
        num("estado_FS_%gh" % t, v(fs[t], 4), rel(os.path.join(d, "fs.csv"), B))
        ax.set_aspect("equal")
        ax.set_xlim(5, 65)
        ax.set_ylim(0, 31)
        ax.grid(False)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        letra(ax, "abcd"[k])
    for ax in axs[1]:
        ax.set_xlabel("$x$ (m)")
    for ax in axs[:, 0]:
        ax.set_ylabel("$y$ (m)")
    leg = [Line2D([], [], color=AZUL_CLARO, lw=6, label="$p_w \\geq 0$"),
           Line2D([], [], color=PALETA[0], lw=1.0, label="$p_w = 0$"),
           Line2D([], [], color=PALETA[1], marker="o", ms=3, lw=0,
                  label="banda: $\\sqrt{J_2(\\varepsilon_p)} \\geq 0{,}1$ do maior")]
    fig.legend(handles=leg, loc="upper center", ncol=3, bbox_to_anchor=(0.5, 1.0))
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    salva(fig, pasta, "fig_estado")


NOME_SOLO = {"f10": "$f_{10,-4}$", "f50": "$f_{50,-5}$", "f100": "$f_{100,-6}$"}


def titulo_solo(solo, L):
    return "%s, %s mm/h" % (NOME_SOLO[solo], g(L.I))


def fig_validacao(pasta, logs):
    fig, axs = plt.subplots(1, 3, figsize=(LARGURA, 2.35))
    for k, (ax, solo) in enumerate(zip(axs, ("f10", "f50", "f100"))):
        faixa_chuva(ax)
        linha_um(ax)
        L = logs[solo]
        ax.plot([r["t"] for r in L], [r["fs"] for r in L], color=PALETA[0], lw=1.4, zorder=3, label="este trabalho")
        R = RAHARDJO_FIG5[solo]
        ax.plot(list(R), list(R.values()), "o", ms=4.2, color=PALETA[1], mec="white", mew=0.8, zorder=4,
                label="Rahardjo et al. (Fig. 5)")
        ax.set_title(titulo_solo(solo, L), fontsize=8, color=TINTA)
        ax.set_xlim(0, 48)
        ax.set_xticks(np.arange(0, 49, 12))
        ax.set_xlabel("$t$ (h)")
        letra(ax, "abc"[k])
    axs[0].set_ylabel("$FS$")
    axs[0].set_ylim(0.7, 2.25)
    axs[1].set_ylim(0.95, 2.2)
    axs[2].set_ylim(1.84, 2.14)
    axs[1].legend(loc="upper right", bbox_to_anchor=(1.0, 0.93))
    fig.tight_layout(w_pad=0.6)
    salva(fig, pasta, "fig_validacao")


def fig_gim_srm(pasta, srm50, gim50, srm10, gim10):
    fig, axs = plt.subplots(1, 3, figsize=(LARGURA, 2.35))
    for k, (ax, s, gm, fim, nome) in enumerate(((axs[0], srm50, gim50, 48, titulo_solo("f50", gim50)),
                                                (axs[1], srm10, gim10, 30, titulo_solo("f10", gim10)))):
        faixa_chuva(ax)
        linha_um(ax)
        ax.plot([r["t"] for r in s if r["t"] <= fim], [r["fs"] for r in s if r["t"] <= fim], color=PALETA[0],
                label="SRM: $FS$")
        ax.plot([r["t"] for r in gm], [r["fs"] for r in gm], color=PALETA[1], label="GIM: $\\Gamma$")
        ax.set_xlim(0, fim)
        ax.set_xticks(np.arange(0, fim + 1, 12 if fim == 48 else 6))
        ax.set_xlabel("$t$ (h)")
        ax.set_title(nome, fontsize=8, color=TINTA)
        letra(ax, "ab"[k])
    axs[0].set_ylabel("$FS$ ou $\\Gamma$ (escala log)")
    for ax in axs[:2]:
        ax.set_yscale("log")
        ax.set_ylim(0.6, 6.)
        ax.set_yticks([0.6, 0.8, 1, 1.5, 2, 3, 5])
        ax.set_yticklabels(["0,6", "0,8", "1", "1,5", "2", "3", "5"])
        ax.minorticks_off()
    axs[0].legend(loc="upper right", bbox_to_anchor=(1.0, 0.93))
    ax = axs[2]
    for cor, rot, L in ((PALETA[0], "SRM", srm50), (PALETA[1], "GIM", gim50)):
        for t, mk in ((0., "o"), (24., "s")):
            n = next(r["niveis"] for r in L if r["t"] == t)
            y = 100. * (np.array(n) / n[-1] - 1.)
            ax.plot(range(len(n)), y, color=cor, marker=mk, ms=4, mec="white", mew=0.7,
                    label="%s, %s h" % (rot, g(t)))
    ax.axhline(0., color=MUDO, lw=0.8)
    ax.set_xticks(range(5))
    ax.set_xlabel("nível de refinamento")
    ax.set_ylabel("acima do último nível (%)")
    ax.set_title("$f_{50,-5}$: valor por nível", fontsize=8, color=TINTA)
    ax.legend(loc="upper right")
    letra(ax, "c")
    fig.tight_layout(w_pad=0.6)
    salva(fig, pasta, "fig_gim_srm")


def sigmoide(I):
    a, b, c, d = SIGMOIDE
    return a + b / (1. + np.exp(-(np.asarray(I, float) - c) / d))


def fig_intensidade(pasta, srmI, gimI, pilotos, fsk, kmin, Iref, ks_mmh):
    fig, axs = plt.subplots(1, 3, figsize=(LARGURA, 2.45), gridspec_kw={"width_ratios": [1.15, 1.0, 1.0]})
    ax = axs[0]
    Is = np.linspace(15, 40, 200)
    ax.plot(Is, sigmoide(Is), color=TINTA2, lw=1.0, ls=(0, (4, 2)), label="Rahardjo, Fig. 7b")
    ax.plot(sorted(srmI), [srmI[i] for i in sorted(srmI)], color=PALETA[0], marker="o", ms=3.5, mec="white",
            mew=0.6, label="SRM: $FS$")
    ax.plot(sorted(gimI), [gimI[i] for i in sorted(gimI)], color=PALETA[1], marker="s", ms=3.5, mec="white",
            mew=0.6, label="GIM: $\\Gamma$")
    ax.axhline(CERON["gamma_det"], color=MUDO, lw=0.8)
    ax.text(39.6, CERON["gamma_det"] - 0.008, "Ceron:\n$\\Gamma = %s$" % v(CERON["gamma_det"], 3).replace(",", "{,}"),
            fontsize=6.2, color=TINTA2, va="top", ha="right")  # em duas linhas: a direita da linha do I_ref
    ax.vlines(Iref, 0.95, 1.46, color=MUDO, lw=0.8)
    ax.text(Iref, 1.47, "$I_{ref}$", fontsize=6.5, color=TINTA2, ha="center", va="bottom")
    linha_um(ax)
    ax.set_xlim(15, 40)
    ax.set_ylim(0.95, 1.92)
    ax.set_xlabel("$I$ (mm/h)")
    ax.set_ylabel("$FS$ ou $\\Gamma$ em 24 h")
    ax.set_title("determinístico, $f_{50,-5}$", fontsize=8, color=TINTA)
    ax.legend(loc="upper right")
    letra(ax, "a")
    ax = axs[1]
    rng = np.random.default_rng(1)
    rot = []
    for k, (nome, cor, fs, fdet) in enumerate(pilotos):
        x = k + rng.uniform(-0.13, 0.13, len(fs))
        ax.scatter(x, fs, s=13, color=cor, edgecolor="white", linewidth=0.6, zorder=3)
        ax.plot([k - 0.3, k + 0.3], [fdet, fdet], color=TINTA, lw=1.4, zorder=4)
        ax.plot([k - 0.3, k + 0.3], [np.mean(fs)] * 2, color=TINTA, lw=0.8, ls=(0, (2, 1.5)), zorder=4)
        rot.append(nome)
    linha_um(ax)
    ax.set_xticks(range(len(rot)))
    ax.set_xticklabels(rot, fontsize=7)
    ax.set_xlabel("$I$ (mm/h)")
    ax.set_xlim(-0.6, len(rot) - 0.4)
    ax.set_ylabel("$FS$ ou $\\Gamma$ em 24 h")
    ax.set_title("pilotos", fontsize=8, color=TINTA)
    ax.legend(handles=[Line2D([], [], color=TINTA, lw=1.4, label="determinístico"),
                       Line2D([], [], color=TINTA, lw=0.8, ls=(0, (2, 1.5)), label="média")], loc="upper left")
    letra(ax, "b")
    ax = axs[2]
    ks = sorted(fsk)
    ax.plot(ks, [fsk[k] for k in ks], color=PALETA[0], marker="o", ms=3.5, mec="white", mew=0.6)
    ax.axvline(Iref / ks_mmh, color=MUDO, lw=0.8)
    ax.text(Iref / ks_mmh * 1.06, 1.33, "$I/k_s$", fontsize=6.5, color=TINTA2)
    ax.annotate("mínimo: $k/k_s \\approx %s$" % v(kmin[0], 2).replace(",", "{,}"), xy=kmin, xytext=(1.25, 1.135),
                fontsize=6.5, color=TINTA2, arrowprops=dict(arrowstyle="-", color=MUDO, lw=0.5))
    ax.plot([kmin[0]], [kmin[1]], marker="v", ms=5, color=TINTA, lw=0)
    ax.set_xscale("log")
    ax.set_xticks([0.25, 0.5, 1, 2, 4])
    ax.set_xticklabels(["0,25", "0,5", "1", "2", "4"])
    ax.minorticks_off()
    ax.set_xlabel("$k/k_s$ (solo homogêneo)")
    ax.set_ylabel("$FS$ em 24 h")
    ax.set_title("SRM, $I = %s$ mm/h" % g(Iref).replace(",", "{,}"), fontsize=8, color=TINTA)
    letra(ax, "c")
    fig.tight_layout(w_pad=0.5)
    salva(fig, pasta, "fig_intensidade")


# ---------------------------------------------------------------- tabelas e numeros
def tabela(nome, linhas):
    print("\n%% ---- tabela: %s" % nome)
    for linha in linhas:
        print(linha)


def main(argv):
    B = os.path.abspath(argv[1]) if len(argv) > 1 else BUILD_PADRAO
    pasta = os.path.abspath(argv[2]) if len(argv) > 2 else os.path.join(AQUI, "figs")
    os.makedirs(pasta, exist_ok=True)
    mc = os.path.join(B, "mc")
    pts, hw, beta, smax, dt = geometria()
    for nome, val in (("Hw", hw), ("beta_nivel", beta), ("succao_max", smax), ("dt_transiente", dt)):
        num(nome, g(val), "RainTransient.h")

    # ---- 3. deterministico contra o Rahardjo
    logs = {s: le_log(os.path.join(B, "rahardjo_%s.log" % s)) for s in ("f10", "f50", "f100")}
    linhas = []
    for solo in ("f10", "f50", "f100"):
        L, R, fonte = logs[solo], RAHARDJO_FIG5[solo], "rahardjo_%s.log" % solo
        fs = {r["t"]: r["fs"] for r in L}
        dif = [100. * (fs[t] / R[t] - 1.) for t in R]
        menor = min(L, key=lambda r: (r["fs"], r["t"]))
        tmin = [r["t"] for r in L if r["fs"] - menor["fs"] < 1.5e-4]  # ate 0,0001 acima do minimo
        ts, ys = [r["t"] for r in L], [r["fs"] for r in L]
        subida = max(ys[i + 1] - ys[i] for i in range(len(ys) - 1) if ts[i + 1] <= 24.)
        num("%s_FS0" % solo, v(fs[0.]), fonte)
        num("%s_FSmin" % solo, v(menor["fs"]), fonte)
        num("%s_tmin" % solo, "%g-%g h" % (tmin[0], tmin[-1]), fonte)
        num("%s_FS24" % solo, v(fs[24.]), fonte)
        num("%s_FS48" % solo, v(fs[48.]), fonte)
        num("%s_maior_subida_na_chuva" % solo, v(subida), fonte)
        num("%s_dif_media_pct" % solo, v(np.mean(dif), 1), fonte + " e Fig. 5")
        num("%s_dif_max_pct" % solo, v(max(dif, key=abs), 1), fonte + " e Fig. 5")
        rmin = min(R.values())
        linhas.append(r"%s & %s & %s & %s & %s (%s) & %s & %s & %s & %d & %s & %s \\" % (
            NOME_SOLO[solo], g(L.I).replace(",", "{,}"), v(fs[0.]), v(R[0], 2),
            v(menor["fs"]), ("%g" % tmin[0]) + (("--%g" % tmin[-1]) if len(tmin) > 1 else "") + " h", v(rmin, 2),
            v(fs[48.]), v(R[48], 2), len(R), v(np.mean(dif), 1), v(max(dif, key=abs), 1)))
    tabela("validacao (secao 3)", linhas)
    num("F_ini_tab7", v(F_INI_TAB7, 3), "Rahardjo, Tabela 7 (digitado)")
    num("f10_FS0_vs_Fini_pct", v(100. * (logs["f10"][0]["fs"] / F_INI_TAB7 - 1.), 1), "rahardjo_f10.log e Tabela 7")
    L10 = logs["f10"]
    num("f10_cruza_1_h", v(cruza_um([r["t"] for r in L10], [r["fs"] for r in L10]), 2), "rahardjo_f10.log")
    plat = [r for r in L10 if 7 <= r["t"] <= 24]
    num("f10_plateau_pwmax", v(plat[0]["pw_max"], 1), "rahardjo_f10.log")
    num("f10_plateau_face", "%d" % plat[0]["face"], "rahardjo_f10.log")
    L100 = {r["t"]: r["fs"] for r in logs["f100"]}
    num("f100_FSmin_abaixo_FS24_pct", v(100. * (1. - min(L100.values()) / L100[24.]), 2), "rahardjo_f100.log")
    # validacao do Monte Carlo contra o rahardjo.cpp e controle sem chuva
    va = os.path.join(mc, "valida_0_-1.csv")
    v0, v24 = det(va, 0.), det(va, 24.)
    num("valida_mc_tri3_FS0", v(v0["fs"]), rel(va, B))
    num("valida_mc_tri3_FS24", v(v24["fs"]), rel(va, B))
    num("valida_eq_transiente_tri3", "%d" % v24["eq_transiente"], rel(va, B))
    num("valida_dif_pct", v(100. * (v24["fs"] / logs["f50"][24]["fs"] - 1.), 2), rel(va, B) + " e rahardjo_f50.log")
    sc = os.path.join(mc, "semchuva_0_-1.csv")
    s0, s24 = det(sc, 0.), det(sc, 24.)
    num("semchuva_FS0", v(s0["fs"]), rel(sc, B))
    num("semchuva_FS24", v(s24["fs"]), rel(sc, B))
    num("semchuva_dif_pct", v(100. * (s24["fs"] / s0["fs"] - 1.), 3), rel(sc, B))
    num("semchuva_pwmin_24", v(s24["pw_min"], 1), rel(sc, B))
    num("semchuva_niveis_0", "/".join(v(x, 1) for x in s0["nivel"]), rel(sc, B))
    num("semchuva_niveis_24", "/".join(v(x, 1) for x in s24["nivel"]), rel(sc, B))
    num("eq_transiente_tri2", "%d" % s0["eq_transiente"], rel(sc, B))
    num("tri2_vs_tri3_FS0_pct", v(100. * (s0["fs"] / v0["fs"] - 1.), 2), rel(sc, B) + " e " + rel(va, B))
    # malha do transiente do rahardjo.cpp (tri3): triangulos e nos da superficie (o f10 tem todos na face no patamar)
    d3 = os.path.join(B, "f50", "chuva_f50_36")
    V3 = np.loadtxt(os.path.join(d3, "vertices.csv"), delimiter=",", skiprows=1)
    T3 = np.loadtxt(os.path.join(d3, "triangles.csv"), delimiter=",", skiprows=1)
    num("malha_tri3_triangulos", "%d" % len(T3), rel(os.path.join(d3, "triangles.csv"), B))
    sup = [(x, y) for _, x, y in V3 if abs(y - np.interp(x, [pts[5][0], pts[4][0], pts[3][0], pts[2][0]],
                                                         [pts[5][1], pts[4][1], pts[3][1], pts[2][1]])) < 1e-6]
    num("malha_tri3_nos_superficie", "%d" % len(sup), rel(os.path.join(d3, "vertices.csv"), B))
    eqs = [r["eq"] for s in logs for r in logs[s]]
    num("eq_fator_min_max", "%d a %d" % (min(eqs), max(eqs)), "rahardjo_f*.log")

    # ---- 4. GIM x SRM
    gim50 = le_log(os.path.join(B, "gim_f50.log"))
    gim10 = le_log(os.path.join(B, "gim_f10.log"))
    srm6 = le_log(os.path.join(B, "srm_f50_6h_sozinho.log"))
    S50 = {r["t"]: r for r in logs["f50"]}
    S10 = {r["t"]: r for r in logs["f10"]}
    G50 = {r["t"]: r for r in gim50}
    G10 = {r["t"]: r for r in gim10}
    linhas = []
    for t5, t1 in zip((0, 6, 12, 18, 24, 36, 48), (0, 3, 6, 7, 24, 27, 30)):  # f50 e f10 lado a lado
        linhas.append(r"%g & %s & %s & %s & & %g & %s & %s & %s \\" % (
            t5, v(S50[t5]["fs"]), v(G50[t5]["fs"]), v(G50[t5]["fs"] / S50[t5]["fs"], 2),
            t1, v(S10[t1]["fs"]), v(G10[t1]["fs"]), v(G10[t1]["fs"] / S10[t1]["fs"], 2)))
    tabela("gim x srm (secao 4)", linhas)
    for nome, L, fonte in (("srm", logs["f10"], "rahardjo_f10.log"), ("gim", gim10, "gim_f10.log")):
        t_ = [r["t"] for r in L]
        f_ = [r["fs"] for r in L]
        num("f10_%s_cruza_desce_h" % nome, v(cruza_um(t_, f_, True), 2), fonte)
        num("f10_%s_cruza_sobe_h" % nome, v(cruza_um(t_, f_, False), 2), fonte)
        num("f10_%s_plateau" % nome, v(next(r["fs"] for r in L if r["t"] == 12)), fonte)
    for t in (0., 24.):
        for nome, D, fonte in (("srm", S50, "rahardjo_f50.log"), ("gim", G50, "gim_f50.log")):
            n = D[t]["niveis"]
            num("f50_%s_%gh_queda_ultimo_nivel_pct" % (nome, t), v(100. * (n[-2] - n[-1]) / n[-2], 1), fonte)
            num("f50_%s_%gh_nivel0_sobre_final_pct" % (nome, t), v(100. * (n[0] / n[-1] - 1.), 0), fonte)
    num("srm_sozinho_s_por_hora_0_6h", v(np.mean([r["s"] for r in srm6]), 1), "srm_f50_6h_sozinho.log")
    num("gim_s_por_hora_0_6h", v(np.mean([r["s"] for r in gim50 if r["t"] <= 6]), 1), "gim_f50.log")
    num("gim_s_por_hora_0_48h", v(np.mean([r["s"] for r in gim50]), 1), "gim_f50.log")
    num("srm_sozinho_igual_ao_log", str(all(abs(r["fs"] - S50[r["t"]]["fs"]) < 1e-12 for r in srm6)),
        "srm_f50_6h_sozinho.log e rahardjo_f50.log")
    num("gim_f50_FS0", v(G50[0]["fs"]), "gim_f50.log")
    num("gim_f50_FS24", v(G50[24]["fs"]), "gim_f50.log")
    num("gim_f50_tmin", "%g" % min(gim50, key=lambda r: r["fs"])["t"], "gim_f50.log")

    # ---- 5. intensidade
    srmI, gimI, fontesI = {}, {}, {}
    for arq in sem_mapa(glob.glob(os.path.join(mc, "varre_I*_0_-1.csv"))):
        I = float(re.search(r"varre_I([\d.]+)_0_-1", arq)[1])
        srmI[I], fontesI[("srm", I)] = det(arq)["fs"], rel(arq, B)
    for arq in sem_mapa(glob.glob(os.path.join(mc, "gim_I", "srm_I*_0_-1.csv"))):
        I = float(re.search(r"srm_I([\d.]+)_0_-1", arq)[1])
        srmI[I], fontesI[("srm", I)] = det(arq)["fs"], rel(arq, B)
    for arq in sem_mapa(glob.glob(os.path.join(mc, "gim_I", "gim_I*_0_-1.csv"))):
        I = float(re.search(r"gim_I([\d.]+)_0_-1", arq)[1])
        gimI[I], fontesI[("gim", I)] = det(arq)["fs"], rel(arq, B)
    linhas = []
    for I in sorted(set(srmI) | set(gimI)):
        sg = float(sigmoide(I))
        linhas.append(r"%s & %s & %s & %s & %s \\" % (
            g(I), v(srmI[I]) if I in srmI else "--", v(gimI[I]) if I in gimI else "--", v(sg, 3),
            v(100. * (srmI[I] / sg - 1.), 1) if I in srmI else "--"))
    tabela("FS(24 h) com I (secao 5)", linhas)
    difs = [100. * (srmI[I] / float(sigmoide(I)) - 1.) for I in srmI]
    num("srm_vs_sigmoide_pct", "%s a %s" % (v(min(difs), 1), v(max(difs), 1)), "mc/varre_I*, mc/gim_I/srm_I*")
    num("valida_tri2_I36_FS24", v(srmI[36.]), fontesI[("srm", 36.)])
    num("tri2_vs_tri3_I36_pct", v(100. * (srmI[36.] / v24["fs"] - 1.), 2), fontesI[("srm", 36.)] + " e " + rel(va, B))
    for nome, D in (("srm", srmI), ("gim", gimI)):  # ruido: residuo de uma parabola em I
        Is = np.array(sorted(D))
        ys = np.array([D[i] for i in Is])
        res = ys - np.polyval(np.polyfit(Is, ys, 2), Is)
        num("%s_ruido_parabola_rms" % nome, v(math.sqrt(np.mean(res ** 2)), 4), "mc/varre_I*, mc/gim_I/*")
    for I in (30.95, 31.05, 31.1):
        num("srm_I%g" % I, v(srmI[I]), fontesI[("srm", I)])
        num("gim_I%g" % I, v(gimI[I]), fontesI[("gim", I)])
    num("gim_I31.3", v(gimI[31.3]), fontesI[("gim", 31.3)])
    num("gim_I36", v(gimI[36.]), fontesI[("gim", 36.)])
    # referencia da campanha (submeter.sh) e k_s do f50 (mc_parametros.txt de um caso)
    sub = texto(os.path.join(PROJ, "cluster", "submeter.sh"))
    Iref = float(re.search(r'IREF="\$\{IREF:-([\d.]+)\}"', sub)[1])
    nos = int(re.search(r'NOS="\$\{1:-(\d+)\}"', sub)[1])
    paralelo = int(re.search(r'PARALELO="\$\{PARALELO:-(\d+)\}"', sub)[1])
    bloco = int(re.search(r'BLOCO="\$\{BLOCO:-(\d+)\}"', sub)[1])
    num("Iref", g(Iref), "cluster/submeter.sh")
    num("nos_x_processos", "%d x %d" % (nos, paralelo), "cluster/submeter.sh")
    num("bloco", "%d" % bloco, "cluster/submeter.sh")
    jt = re.search(r"#SBATCH --time=(\S+)", texto(os.path.join(PROJ, "cluster", "job.sbatch")))[1]
    num("job_tempo", jt, "cluster/job.sbatch")
    te = re.search(r"FS\(0 h\) ~([\d,]+) e FS\(24 h\) ~([\d,]+)", texto(os.path.join(PROJ, "cluster", "teste.sbatch")))
    num("teste_FS0_FS24_esperados", "%s e %s" % (te[1], te[2]), "cluster/teste.sbatch")
    arqpar = os.path.join(mc, "gim_I", "gim_I31.05_parametros.txt")
    par = dict(lin.split("=", 1) for lin in texto(arqpar).splitlines() if "=" in lin)
    ks_mmh = float(par["KS"]) * 3.6e6
    num("KL_variancia_s1_M2000_pct", v(100. * float(par["KL_VARIANCIA"]), 1), "mc/gim_I/gim_I31.05_parametros.txt")
    num("KL_M", par["KL_M"], "mc/gim_I/gim_I31.05_parametros.txt")
    num("semente", par["SEMENTE"], "mc/gim_I/gim_I31.05_parametros.txt")
    num("sigmoide_Iref", v(float(sigmoide(Iref)), 3), "sigmoide da Fig. 7b (digitada)")
    num("sigmoide_I36", v(float(sigmoide(36.)), 3), "sigmoide da Fig. 7b (digitada)")
    num("gim_Iref_vs_ceron_pct", v(100. * (gimI[Iref] / CERON["gamma_det"] - 1.), 1), fontesI[("gim", Iref)])
    num("gim_Iref_menos_ceron", v(CERON["gamma_det"] - gimI[Iref], 4), fontesI[("gim", Iref)])
    # pilotos
    # (metodo, arquivo das realizacoes, instante, arquivo do deterministico); a intensidade sai do nome do arquivo
    pil = [("SRM", "piloto_I24_0_9.csv", 0., sc),
           ("SRM", "piloto_I24_0_9.csv", 24., os.path.join(mc, "varre_I24_0_-1.csv")),
           ("SRM", "piloto_I27_0_9.csv", 24., os.path.join(mc, "varre_I27_0_-1.csv")),
           ("SRM", os.path.join("gim_I", "piloto_srm_I31.05_0_7.csv"), 24.,
            os.path.join(mc, "gim_I", "srm_I31.05_0_-1.csv")),
           ("GIM", os.path.join("gim_I", "piloto_gim_I31.05_0_7.csv"), 24.,
            os.path.join(mc, "gim_I", "gim_I31.05_0_-1.csv"))]
    linhas, desenho = [], []
    for metodo, arq, t, arqdet in pil:
        A = os.path.join(mc, arq)
        I = float(re.search(r"_I([\d.]+)_\d+_\d+\.csv$", arq)[1])
        fs = [r["fs"] for r in le_csv(A) if r["real"] >= 0 and abs(r["t"] - t) < 1e-9 and r["situacao"] == "ok"]
        e, fd = estat(fs), det(arqdet, t)["fs"]
        rot = "%s_I%g_%gh" % (metodo.lower(), I, t)
        linhas.append(r"%s, $I = %s$, %g~h & %d & %s & %s & %s & %s & %s & %s & %d & %s \\" % (
            metodo, g(I).replace(",", "{,}"), t, e["n"], v(fd), v(e["media"]), v(100. * (e["media"] / fd - 1.), 1),
            v(e["dp"]), v(100. * e["cov"], 1), v(e["min"]), e["nf"], v(100. * e["pf_log"], 2)))
        num("piloto_%s_media" % rot, v(e["media"]), rel(A, B))
        num("piloto_%s_cov_pct" % rot, v(100. * e["cov"], 1), rel(A, B))
        num("piloto_%s_pf_log_pct" % rot, v(100. * e["pf_log"], 2), rel(A, B))
        num("piloto_%s_media_sobre_det_pct" % rot, v(100. * (e["media"] / fd - 1.), 1), rel(A, B))
        if t == 24.:
            desenho.append(("%s\n%s" % (g(I), metodo), PALETA[1] if metodo == "GIM" else PALETA[0], fs, fd))
    tabela("pilotos (secao 5)", linhas)
    sa = {r["real"]: r["fs"] for r in le_csv(os.path.join(mc, "gim_I", "piloto_srm_I31.05_0_7.csv"))}
    ga = {r["real"]: r["fs"] for r in le_csv(os.path.join(mc, "gim_I", "piloto_gim_I31.05_0_7.csv"))}
    comuns = sorted(set(sa) & set(ga))
    num("pilotos_31.05_mesmo_menor", str(min(comuns, key=sa.get) == min(comuns, key=ga.get)) + " (real %d)" %
        min(comuns, key=sa.get), "mc/gim_I/piloto_*_I31.05_0_7.csv")
    # custo em serie (estacao)
    lg = texto(os.path.join(mc, "serial_I24.log"))
    mt = re.search(r"tempo medio por realizacao: transiente ([\d.]+) s, fator ([\d.]+) s", lg)
    mr = re.search(r"realizacao\s+0: .*\| ([\d.]+) s \(transiente", lg)
    sr = le_csv(os.path.join(mc, "serial_I24_0_0.csv"))
    tT, tF = float(mt[1]), float(mt[2]) / len(sr)
    num("serial_total_s", v(float(mr[1]), 1), "mc/serial_I24.log")
    num("serial_transiente_s", v(tT, 1), "mc/serial_I24.log")
    num("serial_srm_por_instante_s", v(tF, 1), "mc/serial_I24.log (fator / %d instantes)" % len(sr))
    for r in sr:
        num("serial_srm_%gh_s" % r["t"], v(r["s_fator"], 1), "mc/serial_I24_0_0.csv")
    # FS com k (solo homogeneo, I = 31,05)
    fsk, fontek, extra = {}, {}, {}
    for arq in sem_mapa(glob.glob(os.path.join(mc, "form", "v_*_0_-1.csv"))):
        c, phi, k = (float(x) for x in re.search(r"v_([\d.]+),([\d.]+),([\d.]+)_0_-1", arq).groups())
        r = det(arq)
        if c == 10. and phi == 26.:
            fsk[k], fontek[k] = r["fs"], rel(arq, B)
        else:
            extra[(c, phi, k)] = (r["fs"], rel(arq, B))
    ks = sorted(fsk)
    i = int(np.argmin([fsk[k] for k in ks]))
    x3, y3 = np.array(ks[i - 1:i + 2]), np.array([fsk[k] for k in ks[i - 1:i + 2]])
    p = np.polyfit(x3, y3, 2)
    kmin = (-p[1] / (2 * p[0]), np.polyval(p, -p[1] / (2 * p[0])))
    num("FSk_pontos", "; ".join("%s: %s" % (g(k), v(fsk[k])) for k in ks), "mc/form/v_10,26,<k>_0_-1.csv")
    num("FSk_min_parabola_k", v(kmin[0], 2), "mc/form/v_10,26,{0.5,0.75,1}_0_-1.csv")
    num("FSk_min_parabola_FS", v(kmin[1]), "mc/form/v_10,26,{0.5,0.75,1}_0_-1.csv")
    num("I_sobre_ks", v(Iref / ks_mmh, 3), "submeter.sh e KS do parametros.txt")
    num("Ithr_f50_mmh_e_sobre_ks", "%s e %s" % (g(I_THR_F50), v(I_THR_F50 / ks_mmh, 2)), "Rahardjo, Fig. 7b (digitado)")
    for (c, phi, k), (fs, fonte) in sorted(extra.items()):
        num("FS(c=%g,phi=%g,k=%g)" % (c, phi, k), v(fs), fonte)
    num("FS(10,26,1) igual a srm_I31.05", v(fsk[1.] - srmI[Iref], 5), fontek[1.])
    c1, c2 = extra[(10.1, 26., 1.)][0], extra[(10.2, 26., 1.)][0]
    num("FS_saltos_c_10_10.1_10.2", "%s e %s" % (v(c1 - fsk[1.]), v(c2 - c1)), "mc/form/v_10{,.1,.2},26,1_0_-1.csv")

    # ---- 6. FORM e KL
    fj = os.path.join(mc, "form", "trabalho_I31.05", "form.json")
    F = json.loads(texto(fj))
    hist = F["historico"]
    linhas = []
    for h_ in hist:
        linhas.append(r"%d & %s & %s & %s & %s & %s & %s \\" % (
            h_["it"], v(h_["beta"], 3), v(h_["x"][0], 2), v(h_["x"][1], 2), v(h_["x"][2], 3), v(h_["g"] + 1., 4),
            " ".join(v(gr, 3) for gr in h_["grad"])))
    tabela("FORM, iteracoes (secao 6)", linhas)
    num("form_beta", v(F["beta"], 3), rel(fj, B))
    num("form_Pf_pct", v(100. * F["P_f"], 1), rel(fj, B))
    num("form_x", "%s kPa, %s graus, %s" % (v(F["x"][0], 2), v(F["x"][1], 2), v(F["x"][2], 3)), rel(fj, B))
    num("form_u", ", ".join(v(x, 3) for x in F["u"]), rel(fj, B))
    num("form_alfa2", ", ".join(v(x, 3) for x in F["alfa2"]), rel(fj, B))
    num("form_avaliacoes", "%d" % F["avaliacoes"], rel(fj, B))
    num("form_iteracoes", "%d" % len(hist), rel(fj, B))
    num("form_h", g(F["parametros"]["h"]), rel(fj, B))
    num("form_FS_medianas", v(hist[0]["g"] + 1., 4), rel(fj, B))
    num("form_x_medianas", ", ".join(v(x, 3) for x in hist[0]["x"]), rel(fj, B))
    sig_c = math.sqrt(math.log(1 + F["parametros"]["covc"] ** 2))
    num("form_dc_por_passo_kPa", v(hist[0]["x"][0] * sig_c * F["parametros"]["h"], 2), rel(fj, B))
    num("form_dFS_c_por_passo", v(hist[0]["grad"][0] * F["parametros"]["h"], 3), rel(fj, B))
    grad0 = np.array(hist[0]["grad"])
    num("form_sigma_FS_linear", v(float(np.linalg.norm(grad0)), 3), rel(fj, B))
    num("form_cov_FS_linear_pct", v(100. * float(np.linalg.norm(grad0)) / (hist[0]["g"] + 1.), 0), rel(fj, B))
    ms = (10, 20, 50, 100, 200, 500, 1000, 2000)
    linhas = []
    Dx, Dy = pts[1][0] - pts[0][0], pts[4][1] - pts[0][1]
    for s in (1., 2., 5.):
        fr = kl_fracao(Dx, Dy, 20. * s, 2. * s, ms)
        linhas.append(r"%s & %s & %s \\" % (g(s), "%s~$\\times$~%s" % (g(20 * s), g(2 * s)),
                                            " & ".join(v(100. * f, 1) for f in fr)))
        if s == 1.:
            num("KL_s1_m50_m200_pct", "%s/%s" % (v(100. * fr[2], 1), v(100. * fr[4], 1)), "calculado (KLSeparavel.h)")
            num("KL_s1_m2000_calc_pct", v(100. * fr[-1], 1), "calculado (KLSeparavel.h)")
    tabela("KL (secao 6)", linhas)
    for m in (50, 200, 2000):  # FORM nos coeficientes da KL: 3 campos x m termos, diferencas centrais
        num("form_campo_avaliacoes_m%d" % m, "%d" % (2 * 3 * m + 1), "formula 2 x 3m + 1")
        num("form_campo_avaliacoes_m%d_%d_iteracoes" % (m, len(hist)), "%d" % (len(hist) * (2 * 3 * m + 1)),
            "formula e " + rel(fj, B))
    num("CoVPf_N2000_Pf0.25_pct", v(100. * math.sqrt(0.75 / (2000 * 0.25)), 1), "Eq. 68 do Ceron")

    # ---- 7. campanha
    casos = []
    for linha in texto(os.path.join(PROJ, "cluster", "casos.txt")).splitlines():
        if not linha.strip() or linha.startswith("#"):
            continue
        c = linha.split()
        opc = dict(o.split("=", 1) for o in c[3:])
        casos.append({"id": c[0], "fase": int(c[1]), "N": int(c[2]), "tempos": opc["tempos"].split(","),
                      "opc": opc})
    familias = []
    for c in casos:
        fam = "ref" if c["id"] == "ref" else c["id"].split("=")[0]
        if not familias or familias[-1]["fam"] != fam or familias[-1]["tempos"] != c["tempos"]:
            familias.append({"fam": fam, "fase": c["fase"], "valores": [], "N": c["N"], "tempos": c["tempos"],
                             "n": 0})
        f = familias[-1]
        f["valores"].append(c["id"].split("=")[1] if "=" in c["id"] else "I$_{ref}$")
        f["n"] += 1
        assert f["N"] == c["N"]
    nome_fam = {"ref": "referência", "covk": r"$CoV(k_s)$", "covc": r"$CoV(c')$", "covphi": r"$CoV(\varphi')$",
                "s": r"$s$", "alfa": r"$\alpha = k_h/k_v$", "I": r"$I$ (mm/h)"}
    linhas, tot = [], {}
    proc = nos * paralelo
    for f in familias:
        custo = tT + len(f["tempos"]) * tF  # s por realizacao em serie na estacao
        cpu = f["n"] * f["N"] * custo * FERMI_LENTO / 3600.
        tot.setdefault(f["fase"], [0, 0.])
        tot[f["fase"]][0] += f["n"] * f["N"]
        tot[f["fase"]][1] += cpu
        vals = ", ".join(x.replace(".", "{,}") for x in f["valores"])
        if f["fam"] == "ref":
            vals = r"$I = %s$~mm/h" % g(Iref).replace(",", "{,}")
        linhas.append(r"%d & %s & %s & %d & %s & %s & %s & %s \\" % (
            f["fase"], nome_fam[f["fam"]], vals, f["N"], ", ".join(f["tempos"]), v(custo, 1), "%d" % round(cpu),
            v(cpu / proc, 1)))
    linhas.append(r"\midrule")
    for fase in sorted(tot):
        linhas.append(r"\multicolumn{3}{l}{fase %d} & %d & & & %d & %s \\" % (fase, tot[fase][0], round(tot[fase][1]),
                                                                            v(tot[fase][1] / proc, 1)))
    T = sum(x[1] for x in tot.values())
    linhas.append(r"\multicolumn{3}{l}{total} & %d & & & %d & %s \\" % (sum(x[0] for x in tot.values()), round(T),
                                                                      v(T / proc, 1)))
    tabela("campanha (secao 7)", linhas)
    for fase in sorted(tot):
        num("campanha_fase%d_cpu_h" % fase, "%d" % round(tot[fase][1]), "cluster/casos.txt + custo")
        num("campanha_fase%d_relogio_h" % fase, v(tot[fase][1] / proc, 1), "cluster/casos.txt + custo")
    num("campanha_total_cpu_h", "%d" % round(T), "cluster/casos.txt + custo")
    num("campanha_total_relogio_dias", v(T / proc / 24., 1), "cluster/casos.txt + custo")
    num("campanha_casos", "%d" % len(casos), "cluster/casos.txt")
    num("campanha_realizacoes", "%d" % sum(x[0] for x in tot.values()), "cluster/casos.txt")
    for n_ in (1, 2, 5):
        num("custo_%d_instantes_s" % n_, v(tT + n_ * tF, 1), "mc/serial_I24.log")
    for pf in (0.01, 0.1):
        num("N_para_CoVPf_5pct_Pf%g" % pf, "%d" % round((1 - pf) / (pf * 0.05 ** 2)), "Eq. 68 do Ceron")
        for N in (2000, 5000):
            num("CoVPf_N%d_Pf%g_pct" % (N, pf), v(100. * math.sqrt((1 - pf) / (N * pf)), 0), "Eq. 68 do Ceron")

    # ---- figuras
    fig_geometria(B, pasta, pts, hw, beta)
    fig_estado(B, pasta, pts)
    fig_validacao(pasta, logs)
    fig_gim_srm(pasta, logs["f50"], gim50, logs["f10"], gim10)
    fig_intensidade(pasta, srmI, gimI, desenho, fsk, kmin, Iref, ks_mmh)

    print("\n% ---- numeros do texto (chave: valor  [arquivo])")
    for chave, valor, fonte in NUMEROS:
        print("%-42s %-28s [%s]" % (chave, valor, fonte))


if __name__ == "__main__":
    main(sys.argv)
