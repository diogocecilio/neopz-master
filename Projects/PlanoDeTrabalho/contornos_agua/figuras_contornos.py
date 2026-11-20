#!/usr/bin/env python3
"""figuras_contornos.py: as figuras e as tabelas do relatorio_contornos.tex (os contornos da agua: Ceron 2025, o PZ
'reservatorio', o PZ 'ceron' e a chuva do Rahardjo 2007).

So LE o Projects/CeronCinematico (figs/, dados/ e as funcoes de estudo_contorno.py e cinematico.py) e as rodadas em
~/projects/runs/ceron_contorno e ~/projects/runs/ceron_instantaneo. Nao escreve nada fora desta pasta.

    python3 -I figuras_contornos.py      -> figs/*.pdf/.png e tabelas.tex

Todas as figuras sao desenhadas com a largura do texto (6,5 pol) e incluidas com width=\\textwidth: a fonte que sai
aqui (>= 10 pt) e a fonte do PDF. As figuras antigas do CeronCinematico sao copiadas para figs/ e, como foram feitas
largas (15 a 18 pol), entram recortadas (um painel por figura) em paginas em paisagem.
"""
import csv
import math
import os
import shutil
import sys

sys.dont_write_bytecode = True  # nao criar __pycache__ no CeronCinematico

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.tri import Triangulation
from PIL import Image

AQUI = os.path.dirname(os.path.abspath(__file__))
NEOPZ = os.path.dirname(os.path.dirname(os.path.dirname(AQUI)))
CC = os.path.join(NEOPZ, "Projects", "CeronCinematico")
RUNS = os.path.expanduser("~/projects/runs")
FIGS = os.path.join(AQUI, "figs")
os.makedirs(FIGS, exist_ok=True)
sys.path.insert(0, CC)
import cinematico as k  # noqa: E402
import estudo_contorno as ec  # noqa: E402

GW = 9.81
LARG = 6.5  # largura do texto em polegadas (A4, margem 2,2 cm)
plt.rcParams.update({"font.size": 11, "axes.titlesize": 11, "axes.labelsize": 11, "xtick.labelsize": 10,
                     "ytick.labelsize": 10, "legend.fontsize": 10, "figure.dpi": 100, "savefig.dpi": 200,
                     "mathtext.fontset": "dejavusans"})

# estilos dos contornos (cor + tipo de linha: a cor nunca e a unica pista; cada trecho tem rotulo escrito)
AZUL, LARANJA, ROSA, VERDE, VERM, CINZA = "#0072B2", "#E69F00", "#CC79A7", "#009E73", "#D55E00", "#7f7f7f"
EST = {
    "u0": dict(color=AZUL, ls="-", lw=4.5),          # u = 0 (hidrostatica com o nivel na crista)
    "ugy": dict(color=LARANJA, ls="-", lw=4.5),      # u = -gw y (p = 0 na face)
    "uhw": dict(color=ROSA, ls="-", lw=4.5),         # u = -gw hw / hidrostatica do nivel final
    "sem": dict(color=CINZA, ls=(0, (4, 2.5)), lw=3),  # sem fluxo
    "perc": dict(color=VERDE, ls=(0, (6, 2, 1.5, 2)), lw=4),  # face de percolacao
    "hid": dict(color=AZUL, ls="-", lw=4.5),         # hidrostatica do NA (Rahardjo)
    "chuva": dict(color=VERM, ls="-", lw=4.5),       # fluxo da chuva
}
NOME = {"u0": r"$u = 0$  ($p = \gamma_w\,y$)", "ugy": r"$u = -\gamma_w\,y$  ($p = 0$)",
        "uhw": r"$u = -\gamma_w\,h_w$", "sem": "sem fluxo", "perc": "face de percolação (só sai água)",
        "hid": r"$p_w = \gamma_w\,(\mathrm{NA} - y)$", "chuva": r"fluxo da chuva $q = I_r$"}
CX = dict(boxstyle="round,pad=0.25", fc="white", ec="0.6", lw=0.6)


def virg(ax, eixo="y"):
    """numeros do eixo com virgula decimal"""
    from matplotlib.ticker import FuncFormatter
    f = FuncFormatter(lambda v, _: f"{v:g}".replace(".", ",").replace("-", "−"))
    (ax.yaxis if eixo == "y" else ax.xaxis).set_major_formatter(f)


def salva(fig, nome):
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(FIGS, nome + "." + ext), bbox_inches="tight", pad_inches=0.03)
    plt.close(fig)
    print("figs/" + nome)


def seg(ax, xs, ys, est, **kw):
    ax.plot(xs, ys, solid_capstyle="butt", zorder=3, **EST[est], **kw)


def legenda(ax, ests, ncol=2, nomes=None, **kw):
    h = [Line2D([], [], **EST[e]) for e in ests]
    nomes = nomes or {}
    ax.legend(h, [nomes.get(e, NOME[e]) for e in ests], ncol=ncol, frameon=False, handlelength=3.2, **kw)


def rot(ax, x, y, s, **kw):
    kw.setdefault("ha", "center")
    kw.setdefault("va", "center")
    ax.text(x, y, s, bbox=CX, zorder=5, **kw)


def solo(ax, pts):
    ax.fill(*zip(*pts), color="#efe6d2", zorder=0)


def nivel(ax, x0, x1, y, txt, xt, cor=AZUL, ls="--", va="bottom", xtxt=None):
    ax.plot([x0, x1], [y, y], color=cor, ls=ls, lw=1.6, zorder=2)
    ax.plot([xt], [y + 0.0], marker="v", color=cor, ms=9, zorder=4)
    xtxt = xt + 0.8 if xtxt is None else xtxt
    ax.text(xtxt, y + (0.4 if va == "bottom" else -0.4), txt, color="k", va=va, ha="left", fontsize=10)


def cota(ax, x0, y0, x1, y1, txt, off=(0, 0), **kw):
    ax.annotate("", (x1, y1), (x0, y0), arrowprops=dict(arrowstyle="<->", color="0.25", lw=1))
    ax.text((x0 + x1) / 2 + off[0], (y0 + y1) / 2 + off[1], txt, ha="center", va="center", fontsize=10, **kw)


# --------------------------------------------------------------------------------------------------------------------
# (a) Ceron: o contorno identificado (Eqs. 1, 19-21 e Fig. 4b); dominio 14H x 2,8H x 8,4H (Gamma 1,337)
def esquema_ceron():
    H = 5.0
    Le, Ld, D = 14 * H, 2.8 * H, 8.4 * H
    xT = H
    xe, xd, yb = -Le, xT + Ld, -H - D
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(LARG, 7.6), height_ratios=[1.0, 1.05])
    for a in (a1, a2):
        solo(a, [(xe, 0), (0, 0), (xT, -H), (xd, -H), (xd, yb), (xe, yb)])
    # dominio inteiro
    a = a1
    seg(a, [xe, 0], [0, 0], "u0")
    seg(a, [0, xT], [0, -H], "ugy")
    seg(a, [xT, xd], [-H, -H], "uhw")
    seg(a, [xe, xe], [0, yb], "u0")
    seg(a, [xe, xd], [yb, yb], "u0")
    seg(a, [xd, xd], [-H, yb], "sem")
    rot(a, -38, -20, "lateral esquerda\n" + r"$u = 0$ (Fig. 4b)")
    rot(a, -35, yb + 4.2, r"base $u = 0$ (Fig. 4b)")
    rot(a, 7, -30, "lateral\ndireita\nsem fluxo\n(Fig. 4b)", fontsize=10)
    rot(a, -38, 5.5, r"crista $u = 0$  ($\partial\Omega_{u1}$)")
    a.annotate("talude\n(zoom abaixo)", (2.5, -2.5), (-14, -14), ha="center", fontsize=10,
               arrowprops=dict(arrowstyle="->", lw=1))
    cota(a, xe, 10.5, 0, 10.5, r"$14H = 70$ m", off=(0, 2.2))
    cota(a, xT, 3.5, xd, 3.5, r"$2{,}8H$", off=(0, 2.6))
    cota(a, xd + 3, -H, xd + 3, yb, r"$8{,}4H$", off=(3.5, 0), rotation=90)
    a.set_xlim(xe - 3, xd + 9)
    a.set_ylim(yb - 2, 16)
    a.set_aspect("equal")
    a.set_title("(a1) contorno que reproduz a Fig. 4b; domínio 14H × 2,8H × 8,4H (Γ = 1,337)", loc="left",
                fontsize=10.5)
    a.set_xlabel("x (m), origem na aresta da crista")
    a.set_ylabel("cota (m)")
    # zoom
    a = a2
    seg(a, [-12, 0], [0, 0], "u0")
    seg(a, [0, xT], [0, -H], "ugy")
    seg(a, [xT, 17], [-H, -H], "uhw")
    a.fill([0, xT, 17, 17, 0], [0, -H, -H, 0, 0], color="#cfe3f3", zorder=0.5)
    nivel(a, -12, 17, 0, "nível antes (na crista)", 8.5, xtxt=6.2)
    nivel(a, -12, 17, -H, r"nível depois ($h_w = H$, no pé)", 8.5, va="top", xtxt=-11.6)
    rot(a, -6.5, 1.9, r"crista: $p = 0$, $u = 0$" + "\n" + r"($\partial\Omega_{u1}$)")
    rot(a, 0.2, -3.6, r"face: $p = 0$" + "\n" + r"$u = -\gamma_w\,y$" + "\n" + r"($\partial\Omega_{u2}$)", ha="right")
    rot(a, 11.5, -7.6, r"pé: $p = \gamma_w(H - h_w) = 0$" + "\n" + r"$u = -\gamma_w h_w$ ($\partial\Omega_{u3}$)")
    a.text(10.2, -2.3, "reservatório\nantes", ha="center", va="center", fontsize=10, color="#335f87")
    a.text(-5.5, -8.3, "regime permanente saturado\n" + r"$\nabla\cdot(-K\nabla u) = 0$, sem tempo",
           ha="center", va="center", fontsize=10)
    a.set_xlim(-12, 17)
    a.set_ylim(-9.5, 3.8)
    a.set_aspect("equal")
    a.set_title(r"(a2) perto do talude ($H = 5$ m, $\beta = 45°$); y do artigo aponta para baixo", loc="left")
    a.set_xlabel("x (m)")
    a.set_ylabel("cota (m)")
    legenda(a2, ["u0", "ugy", "uhw", "sem"], loc="upper center", bbox_to_anchor=(0.5, -0.2))
    fig.tight_layout()
    salva(fig, "esquema_a_ceron")


# --------------------------------------------------------------------------------------------------------------------
# (b) e (c): o PZ (malha kTalude = {(0,0),(70,0),(70,30),(40,30),(35,35),(0,35)}), em coordenadas com origem na crista
def esquema_pz(reserv):
    H = 5.0
    xe, xT, xd, yb = -35.0, 5.0, 35.0, -35.0
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(LARG, 7.8), height_ratios=[1.0, 0.95])
    for a in (a1, a2):
        solo(a, [(xe, 0), (0, 0), (xT, -H), (xd, -H), (xd, yb), (xe, yb)])
    a = a1
    if reserv:
        seg(a, [xe, 0], [0, 0], "perc")
        seg(a, [0, xT], [0, -H], "perc")
        seg(a, [xT, xd], [-H, -H], "uhw")
        seg(a, [xd, xd], [-H, yb], "uhw")
        seg(a, [xe, xd], [yb, yb], "sem")
        rot(a, 18, -20, "lateral direita\n" + r"$p_w = \gamma_w(\mathrm{nível}(t) - y)$")
        rot(a, 0, yb + 3.5, "base: sem fluxo")
        rot(a, -17, 4, "crista: face de percolação")
    else:
        seg(a, [xe, 0], [0, 0], "u0")
        seg(a, [0, xT], [0, -H], "ugy")
        seg(a, [xT, xd], [-H, -H], "uhw")
        seg(a, [xd, xd], [-H, yb], "sem")
        seg(a, [xe, xd], [yb, yb], "u0")
        rot(a, 27, -20, "lateral direita\nsem fluxo")
        rot(a, 0, yb + 3.5, r"base: $u = 0$, $p_w = \gamma_w(35 - y)$")
        rot(a, -17, 4, r"crista: $p_w = 0$ ($u = 0$), fixo")
    seg(a, [xe, xe], [0, yb], "u0")
    rot(a, -24, -18, "lateral esquerda\n" + r"$p_w = \gamma_w(35 - y)$" + "\n" + r"($u = 0$)")
    cota(a, xe, 9.5, 0, 9.5, "7H = 35 m", off=(0, 2.3))
    cota(a, xT, 9.5, xd, 9.5, "6H = 30 m", off=(0, 2.3))
    cota(a, xd + 3, -H, xd + 3, yb, "6H", off=(2.6, 0), rotation=90)
    a.set_xlim(xe - 2, xd + 7)
    a.set_ylim(yb - 2, 14)
    a.set_aspect("equal")
    a.set_xlabel("x (m), origem na crista (no PZ: x + 35)")
    a.set_ylabel("cota (m) (no PZ: y − 35)")
    a.set_title("(%s1) domínio inteiro (malha tri3)" % ("b" if reserv else "c"), loc="left")
    # zoom
    a = a2
    a.fill([0, xT, 17, 17, 0], [0, -H, -H, 0, 0], color="#cfe3f3", zorder=0.5)
    nivel(a, -12, 17, 0, "nível inicial (na crista)", 8.5, xtxt=6.2)
    nivel(a, -12, 17, -H, "nível final (no pé)", 8.5, va="top", xtxt=-11.6)
    a.annotate("", (14.5, -H + 0.2), (14.5, -0.2), arrowprops=dict(arrowstyle="-|>", color=AZUL, lw=2))
    if reserv:
        seg(a, [-12, 0], [0, 0], "perc")
        seg(a, [0, xT], [0, -H], "perc")
        seg(a, [xT, 17], [-H, -H], "uhw")
        a.text(14.0, -2.3, "desce\nem trebaix\n(0,1 h)", fontsize=10, va="center", ha="right")
        rot(a, -6.2, 3.3, "crista e face acima do nível:\n" + r"$p_w = 0$ onde a água sai;" + "\nsem fluxo onde não sai",
            fontsize=10)
        rot(a, -6.3, -3.0, "abaixo do nível (na descida):\n" + r"$p_w = \gamma_w(\mathrm{nível}(t) - y)$", fontsize=10)
        rot(a, 11, -7.4, r"pé: $p_w = \gamma_w(\mathrm{nível}(t) - y)$", fontsize=10)
    else:
        seg(a, [-12, 0], [0, 0], "u0")
        seg(a, [0, xT], [0, -H], "ugy")
        seg(a, [xT, 17], [-H, -H], "uhw")
        a.text(14.0, -2.3, "desce em\n1e-4 h", fontsize=10, va="center", ha="right")
        rot(a, -6.5, 2.6, r"crista: $p_w = 0$, Dirichlet", fontsize=10)
        rot(a, -0.4, -3.6, "face (Dirichlet; a água\npode entrar):\n"
            + r"$p_w = \gamma_w\max(\mathrm{nível}(t) - y,\,0)$", ha="right", fontsize=10)
        rot(a, 10.6, -7.4, r"pé: $p_w = \gamma_w\max(\mathrm{nível}(t) - y,\,0)$", fontsize=10)
    a.set_xlim(-12, 17)
    a.set_ylim(-9.0, 5.6 if reserv else 3.8)
    a.set_aspect("equal")
    a.set_xlabel("x (m)")
    a.set_ylabel("cota (m)")
    a.set_title("(%s2) perto do talude; " % ("b" if reserv else "c")
                + ("transiente, no MC um só passo dt = 0,1 h, $k = 10^{-6}$ m/s" if reserv else
                   "transiente até 10 h e 100 h = o permanente (5e-4 kPa)"), loc="left", fontsize=10.5)
    # legenda na notacao do PZ (y para cima, crista em y = 35)
    ests = ["u0", "uhw", "perc", "sem"] if reserv else ["u0", "ugy", "uhw", "sem"]
    nomes = {"u0": r"$p_w = \gamma_w(35 - y)$, ou seja $u = 0$",
             "uhw": r"$p_w = \gamma_w(\mathrm{nível}(t) - y)$ (no fim, $u = -\gamma_w h_w$)",
             "ugy": r"$p_w = 0$, ou seja $u = -\gamma_w(35 - y)$"}
    legenda(a2, ests, nomes=nomes, loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=1)
    fig.tight_layout()
    salva(fig, "esquema_b_pz_reservatorio" if reserv else "esquema_c_pz_ceron")


# --------------------------------------------------------------------------------------------------------------------
# (d) e (e) Rahardjo: geometria {{0,0},{70,0},{70,20},{40,20},{30,30},{0,30}}, NA a 5 m do pe, 7 graus
def na(x):
    return 20.0 - 5.0 + (40.0 - x) * math.tan(math.radians(7.0))


def base_rahardjo(a, plano):
    P = [(0, 0), (70, 0), (70, 20), (40, 20), (30, 30), (0, 30)]
    solo(a, P)
    xs = np.linspace(0, 70, 2)
    a.plot(xs, na(xs), color=AZUL, ls="--", lw=1.6, zorder=2)
    a.plot([20], [na(20) + 0.05], marker="v", color=AZUL, ms=9, zorder=4)
    a.text(19.8, na(20) - 3.4, "NA inicial (7°)", fontsize=10, ha="center")
    # linha de -75 kPa (7,65 m acima do NA)
    xs = np.linspace(0, 36.8, 50)
    a.plot(xs, na(xs) + 75 / GW, color="0.35", ls=":", lw=1.4)
    a.text(6, na(6) + 75 / GW - 0.5, r"$p_w = -75$ kPa (limite)", fontsize=10, va="top")
    yh, ye = na(0), na(70)
    seg(a, [0, 0], [0, yh], "hid")
    seg(a, [0, 0], [yh, 30], "sem")
    seg(a, [70, 70], [0, ye], "hid")
    seg(a, [70, 70], [ye, 20], "sem")
    seg(a, [0, 70], [0, 0], "sem")
    est_sup = "chuva"
    seg(a, [0, 30], [30, 30], est_sup)
    seg(a, [30, 40], [30, 20], est_sup)
    seg(a, [40, 70], [20, 20], est_sup)
    for p, nm, d in [((0, 30), "a", (-1.6, 0.6)), ((30, 30), "b", (0.2, 1.1)), ((40, 20), "c", (0.3, 0.9)),
                     ((70, 20), "d", (0.8, 0.6)), ((70, ye), "e", (0.9, 0)), ((70, 0), "f", (0.9, -0.6)),
                     ((0, 0), "g", (-1.5, -0.6)), ((0, yh), "h", (-1.5, 0))]:
        a.plot(*p, "ko", ms=4, zorder=6)
        a.text(p[0] + d[0], p[1] + d[1], nm, fontsize=11, fontweight="bold", ha="center", va="center")
    # setas da chuva (verticais)
    for x in np.arange(3, 70, 4.5):
        y = 30 if x <= 30 else (30 - (x - 30) if x <= 40 else 20)
        a.annotate("", (x, y + 0.3), (x, y + 3.3), arrowprops=dict(arrowstyle="-|>", color=VERM, lw=1.4))
    a.set_xlim(-4, 79.5)
    a.set_ylim(-6, 37.5)
    a.set_aspect("equal")
    a.set_xlabel("x (m)")
    a.set_ylabel("y (m)")


def esquema_rahardjo_codigo():
    fig, a = plt.subplots(figsize=(LARG, 6.6))
    base_rahardjo(a, False)
    rot(a, 15, 34.8, r"crista ab: $q = I_r$ por metro")
    a.annotate("face bc: " + r"$q = I_r$" + " por metro\n" + r"inclinado ($\times\sqrt{2}$ da chuva)", (35.6, 24.6),
               (57, 29.2), ha="center", va="center", bbox=CX, zorder=5,
               arrowprops=dict(arrowstyle="->", lw=1.2, color="k", shrinkA=0, shrinkB=2))
    rot(a, 55, 16.5, r"pé cd: $q = I_r$", fontsize=10)
    rot(a, 9, 9.0, "gh: " + r"$p_w = \gamma_w(\mathrm{NA}(x) - y)$" + "\nfixo o tempo todo")
    rot(a, 57, 5.5, "ef: " + r"$p_w = \gamma_w(\mathrm{NA} - y)$")
    rot(a, 35, -3.0, "base fg: sem fluxo")
    rot(a, 7.5, 21.9, "ah: sem fluxo", fontsize=10)
    rot(a, 75.6, 15.8, "de:\nsem\nfluxo", fontsize=10)
    a.set_title("(d) a chuva do Rahardjo como está no código (RainTransient)", loc="left")
    legenda(a, ["chuva", "hid", "sem"], ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.12))
    fig.tight_layout()
    salva(fig, "esquema_d_rahardjo_codigo")


def esquema_rahardjo_plano():
    fig = plt.figure(figsize=(LARG, 8.0))
    g = fig.add_gridspec(2, 2, height_ratios=[1.45, 1.0], width_ratios=[1, 1])
    a = fig.add_subplot(g[0, :])
    base_rahardjo(a, True)
    a.set_title("(e1) o mesmo contorno, com as quatro decisões marcadas", loc="left")
    a.annotate("D1: chuva na face\n(ver e2 e e3)", (35.6, 24.6), (56, 29.2), ha="center", va="center", bbox=CX,
               zorder=5, arrowprops=dict(arrowstyle="->", lw=1.2, color="k", shrinkA=0, shrinkB=2))
    rot(a, 15, 34.8, "D2: face de percolação (mantida)")
    rot(a, 22, 9.0, "D7: NA calculado com o pé\nda malha (yToe = 20)\nem todo driver novo")
    rot(a, 35, -3.4, "D4: FS de 6 em 6 h até 48 h; antes, determinístico até 240 h no f100", fontsize=10)
    legenda(a, ["chuva", "hid", "sem"], ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.17))
    # (e2) a chuva na face: hoje e proposta
    for j, (tit, fator) in enumerate([("(e2) hoje: $q_n = I_r$ por\nmetro inclinado", 1.0),
                                      (r"(e3) proposta: $q_n = I_r\cos\alpha$" + "\n(chuva vertical)", math.cos(math.pi / 4))]):
        b = fig.add_subplot(g[1, j])
        b.fill([0, 4, 0], [4, 0, 0], color="#efe6d2", zorder=0)
        b.plot([0, 4], [4, 0], color=VERM, lw=4)
        if j == 0:   # fluxo normal I por metro inclinado: setas normais a face, todas do mesmo tamanho
            for t in np.linspace(0.5, 3.5, 5):
                b.annotate("", (t, 4 - t), (t + 0.85, 4 - t + 0.85),
                           arrowprops=dict(arrowstyle="-|>", color=VERM, lw=1.8))
            b.text(2.6, 3.9, r"$I_r$ normal à face" + "\nem cada metro\ninclinado: 1,41×\na chuva que cai",
                   fontsize=10, ha="left", va="center")
        else:        # chuva vertical I por metro horizontal
            for t in np.linspace(0.5, 3.5, 5):
                b.annotate("", (t, 4 - t), (t, 4 - t + 1.2), arrowprops=dict(arrowstyle="-|>", color=VERM, lw=1.8))
            b.text(2.4, 3.9, "a chuva que cai\n(BoundaryNodes\ncom horizontal)", fontsize=10, ha="left",
                   va="center")
        b.text(0.25, 0.4, "solo", fontsize=10)
        b.annotate("", (3.3, 0), (4, 0), arrowprops=dict(arrowstyle="-", color="0.3"))
        b.text(3.0, 0.22, r"$\alpha$", fontsize=11)
        b.set_xlim(-0.3, 6.2)
        b.set_ylim(-0.3, 5.3)
        b.set_aspect("equal")
        b.axis("off")
        b.set_title(tit, loc="left", fontsize=11)
    fig.tight_layout()
    salva(fig, "esquema_e_rahardjo_plano")


# --------------------------------------------------------------------------------------------------------------------
# resultados: os campos de u e os Gamma
NIV = np.arange(-50, 0.1, 5.0)


def le_csv(p):
    with open(p) as f:
        return list(csv.DictReader(f))


def campo_pz(pasta):
    tal, campo = k.ler_exportacao(pasta)
    V = np.loadtxt(os.path.join(pasta, "vertices.csv"), delimiter=",", skiprows=1)
    tri = np.loadtxt(os.path.join(pasta, "triangles.csv"), delimiter=",", skiprows=1).astype(int)
    return tal, campo, V, tri


def barra(fig, cs, axs):
    cb = fig.colorbar(cs, ax=axs, orientation="horizontal", shrink=0.8, pad=0.08, aspect=40)
    cb.set_label(r"excesso de poropressão $u = p_w - \gamma_w\,(y_{\mathrm{crista}} - y)$ (kPa)")
    cb.set_ticks(np.arange(-50, 1, 10))


def contornos_permanente():
    """o permanente saturado (Eqs. 19-21) no nosso dominio 7H x 6H x 6H com os contornos A-D do estudo_contorno"""
    lin = le_csv(os.path.join(CC, "dados", "estudo_contorno.csv"))
    G = {r["contorno"][0]: float(r["Gamma"]) for r in lin if r["dominio"].startswith("nosso")}
    tal, xy, tri = ec.malha(35.0, 30.0, 30.0)
    casos = [("A", "A: só crista, face e pé\n(resto sem fluxo)", (None, None, None)),
             ("B", "B: esq. $u = 0$, dir. $-\\gamma_w h_w$,\nbase sem fluxo", ("0", "hw", None)),
             ("C", "C: tudo Dirichlet ($\\partial\\Omega_v = \\emptyset$)\nesq. e base 0, dir. $-\\gamma_w h_w$",
              ("0", "hw", "0")),
             ("D", "D: esq. e base $u = 0$,\ndir. sem fluxo (Fig. 4b)", ("0", None, "0"))]
    fig, ax = plt.subplots(2, 2, figsize=(LARG, 6.4), layout="constrained")
    for a, (c, tit, cc) in zip(ax.ravel(), casos):
        u = ec.permanente(xy, tri, tal, *cc)
        T = Triangulation(xy[:, 0], -xy[:, 1], tri)
        cs = a.tricontourf(T, u, levels=NIV, cmap="viridis", extend="both")
        a.tricontour(T, u, levels=NIV[:-1], colors="k", linewidths=0.35)
        a.plot([-20, 0, 5, 25], [0, 0, -5, -5], color=VERM, lw=1.5)
        a.set_xlim(-20, 25)
        a.set_ylim(-20, 1)
        a.set_aspect("equal")
        a.set_title(tit + f"\nΓ cinemático = {G[c]:.3f}".replace(".", ","), fontsize=10.5)
    for a in ax[1]:
        a.set_xlabel("x (m), origem na crista")
    for a in ax[:, 0]:
        a.set_ylabel("cota (m)")
    barra(fig, cs, ax)
    salva(fig, "permanente_contornos")


def gamma_dominio():
    lin = le_csv(os.path.join(CC, "dados", "estudo_contorno.csv"))
    doms = []
    for r in lin:
        if r["dominio"] not in doms:
            doms.append(r["dominio"])
    rot_d = {doms[0]: "notebook (3,5H; 3H; 3H)", doms[1]: "nosso (7H; 6H; 6H)", doms[2]: "Fig. 4 aprox. (13H; 2H; 8H)"}
    fig, (a, b) = plt.subplots(1, 2, figsize=(LARG, 3.6), layout="constrained", width_ratios=[1.15, 1])
    marc = ["o", "s", "^"]
    cores = [AZUL, LARANJA, ROSA]
    for i, d in enumerate(doms):
        g = [float(r["Gamma"]) for r in lin if r["dominio"] == d]
        a.plot(np.arange(4) + (i - 1) * 0.12, g, ls="none", marker=marc[i], ms=8, color=cores[i], label=rot_d[d])
    a.axhline(1.336, color="k", lw=1, ls="--")
    a.text(3.35, 1.338, "artigo\n1,336", fontsize=10, va="bottom", ha="center")
    a.set_xticks(range(4), ["A", "B", "C", "D"])
    a.set_xlabel("contorno")
    virg(a)
    a.set_ylabel("Γ cinemático (permanente)")
    a.set_xlim(-0.5, 3.8)
    a.legend(loc="lower left", fontsize=9.5, frameon=False, bbox_to_anchor=(0, -0.02))
    a.set_ylim(1.15, 1.45)
    a.grid(axis="y", color="0.9")
    # contorno D x tamanho do dominio (README Etapa 3d, passo 1b)
    Ls = [10, 14, 20, 50]
    Gd = [1.316, 1.337, 1.354, 1.379]
    b.plot(Ls, Gd, "o-", color=AZUL, lw=2, ms=8)
    for L, g in zip(Ls, Gd):
        b.text(L * 1.04, g - 0.004, f"{g:.3f}".replace(".", ","), ha="left", va="top", fontsize=10)
    b.axhline(1.336, color="k", lw=1, ls="--")
    b.set_xscale("log")
    b.set_xticks(Ls, [str(x) for x in Ls])
    b.minorticks_off()
    b.set_xlabel("crista até a lateral esq. (em H)")
    b.set_xlim(8.5, 75)
    virg(b)
    b.set_title("contorno D", fontsize=10.5)
    b.set_ylim(1.30, 1.40)
    b.grid(axis="y", color="0.9")
    salva(fig, "gamma_dominio")


def instantaneo():
    base = os.path.join(RUNS, "ceron_instantaneo")
    lin = le_csv(os.path.join(CC, "dados", "instantaneo_pz.csv"))
    G = {(r["agua"], r["t_h"]): r for r in lin}
    tempos = [("0.0001", "t = 0+"), ("0.1", "t = 0,1 h"), ("1", "t = 1 h"), ("10", "t = 10 h")]
    fig, ax = plt.subplots(4, 2, figsize=(LARG, 8.1), layout="constrained", sharex=True, sharey=True)
    for j, (agua, nome) in enumerate([("sat", "saturada"), ("unsat", "não saturada")]):
        for i, (t, rt) in enumerate(tempos):
            a = ax[i, j]
            tal, campo, V, tri = campo_pz(os.path.join(base, f"{agua}_t{t}_real-1"))
            T = Triangulation(V[:, 1] - 35.0, V[:, 2] - 35.0, tri)
            cs = a.tricontourf(T, campo.un, levels=NIV, cmap="viridis", extend="both")
            a.tricontour(T, campo.un, levels=NIV[:-1], colors="k", linewidths=0.35)
            a.plot([-20, 0, 5, 25], [0, 0, -5, -5], color=VERM, lw=1.5)
            a.set_xlim(-20, 25)
            a.set_ylim(-20, 1)
            a.set_aspect("equal")
            r = G[(agua, t)]
            s = f"{nome}, {rt}\nΓ GIM {float(r['Gamma_GIM_upwpa']):.3f}, cin. {float(r['Gamma_cin']):.3f}"
            a.set_title(s.replace(".", ","), fontsize=10.5)
            if j == 0:
                a.set_ylabel("cota (m)")
            if i == 3:
                a.set_xlabel("x (m), origem na crista")
    barra(fig, cs, ax)
    salva(fig, "instantaneo_campos")


def pz_ceron():
    pasta = os.path.join(RUNS, "ceron_contorno", "sat_ceron_real-1")
    tal, campo, V, tri = campo_pz(pasta)
    xy = np.c_[V[:, 1] - 35.0, 35.0 - V[:, 2]]
    u_perm = ec.permanente(xy, tri, tal, "0", None, "0")
    T = Triangulation(V[:, 1] - 35.0, V[:, 2] - 35.0, tri)
    fig, ax = plt.subplots(2, 1, figsize=(LARG, 6.9), layout="constrained", height_ratios=[1.15, 1])
    for i, a in enumerate(ax):
        cs = a.tricontourf(T, campo.un, levels=NIV, cmap="viridis", extend="both")
        a.tricontour(T, u_perm, levels=NIV[1:-1], colors="w", linewidths=1.3, linestyles="--")
        a.tricontour(T, campo.un, levels=NIV[1:-1], colors="k", linewidths=0.5)
        a.plot([-35, 0, 5, 35], [0, 0, -5, -5], color=VERM, lw=1.5)
        a.set_aspect("equal")
        a.set_ylabel("cota (m)")
        if i == 0:
            a.set_xlim(-35, 35)
            a.set_ylim(-35, 1)
            a.set_title("domínio inteiro: PZ (preto) e Laplace independente (branco tracejado)",
                        fontsize=10.5)
        else:
            a.set_xlim(-20, 25)
            a.set_ylim(-20, 1)
            erro = float(np.max(np.abs(campo.un - u_perm)))
            a.set_title(f"zoom; máx |u_PZ − u_Laplace| = {erro:.0e} kPa".replace("e-0", "e-"), fontsize=10.5)
            a.set_xlabel("x (m), origem na crista")
    barra(fig, cs, ax)
    salva(fig, "pz_ceron_campo")


def richardson(g):
    """extrapolacao com tres niveis (razao constante)"""
    d1, d2 = g[1] - g[0], g[2] - g[1]
    r = d2 / d1
    return g[2] + d2 * r / (1 - r), r


def niveis_ceron():
    """Gamma do GIM por nivel (NREF 0 a 5) do PZ 'ceron' (coluna niveis da rodada com nref=5)"""
    r = le_csv(os.path.join(RUNS, "ceron_contorno", "sat_ceron_n5_0_-1.csv"))[0]
    return [float(v) for v in r["niveis"].split(";")]


def niveis_det():
    """Gamma do GIM por nivel (NREF 0 a 3) do PZ 'reservatorio', exportacao det (rodada padrao)"""
    r = le_csv(os.path.join(RUNS, "ceron_contorno", "padrao_0_-1.csv"))[0]
    return [float(v) for v in r["niveis"].split(";")]


def malha_agua():
    """o ponto aberto 1,342 x 1,322: o cinematico com a agua 'ceron' (contorno D, 7H) em tres combinacoes de malha.
    Devolve (Gamma malha propria fina, Gamma na malha do PZ com o u da malha propria fina, max |u_PZ - u_fino|)."""
    from matplotlib.tri import LinearTriInterpolator
    tal, campo = k.ler_exportacao(os.path.join(RUNS, "ceron_contorno", "sat_ceron_real-1"))
    tal2, xy2, tri2 = ec.malha(35.0, 30.0, 30.0, h=0.0625)
    u2 = ec.permanente(xy2, tri2, tal2, "0", None, "0")
    n2, n = len(xy2), len(campo.xy)
    g_fino = k.minimizar(k.CampoMalha(xy2, tri2, np.full(n2, 10.0), np.full(n2, math.radians(30.0)), u2,
                                      20.0 - GW), tal2).Gamma
    ui = np.asarray(LinearTriInterpolator(Triangulation(xy2[:, 0], xy2[:, 1], tri2), u2)(campo.xy[:, 0],
                                                                                         campo.xy[:, 1]))
    du = float(np.max(np.abs(ui - campo.un)))
    g_pz_ufino = k.minimizar(k.CampoMalha(campo.xy, campo.tri, np.full(n, 10.0), np.full(n, math.radians(30.0)), ui,
                                          20.0 - GW), tal).Gamma
    print(f"malha_agua: fino {g_fino:.4f}, malha PZ com u fino {g_pz_ufino:.4f}, max du {du:.3f} kPa")
    return g_fino, g_pz_ufino, du


def convergencia():
    nr = [0, 1, 2, 3, 4]
    g_res = [1.914, 1.580, 1.396, 1.325, 1.297]
    g_cer = niveis_ceron()
    e_res, _ = richardson(g_res[2:])
    e_cer, _ = richardson(g_cer[3:])
    fig, (a, b) = plt.subplots(1, 2, figsize=(LARG, 3.9), layout="constrained")
    a.plot(nr, g_res, "o-", color=LARANJA, lw=2, ms=7, label="PZ 'reservatorio' (0,1 h)")
    a.plot(range(len(g_cer)), g_cer, "s-", color=AZUL, lw=2, ms=7, label="PZ 'ceron' (permanente)")
    a.plot([6.2], [e_res], "o", color=LARANJA, ms=9, mfc="white", mew=2)
    a.plot([6.2], [e_cer], "s", color=AZUL, ms=9, mfc="white", mew=2)
    a.axhline(1.336, color="k", ls="--", lw=1)
    a.text(0.1, 1.345, "artigo 1,336", fontsize=10, va="bottom")
    a.set_xticks([0, 1, 2, 3, 4, 5, 6.2], ["0", "1", "2", "3", "4", "5", "∞"])
    a.set_xlabel("níveis de refinamento (NREF)")
    a.set_ylabel("Γ determinístico (GIM)")
    a.set_ylim(1.2, 1.7)
    a.text(1.3, 1.675, f"0 níveis: 1,914 e {g_cer[0]:.3f}".replace(".", ","), fontsize=10)
    a.legend(loc="lower left", frameon=False, fontsize=10)
    virg(a)
    a.grid(axis="y", color="0.9")
    pf = [0.65, 6.05, 11.80, 15.70]
    b.plot([0, 1, 2, 3], pf, "o-", color=LARANJA, lw=2, ms=7)
    b.errorbar([3], [15.70], yerr=[[15.70 - 14.1], [17.3 - 15.70]], color=LARANJA, capsize=4, lw=1.5)
    b.axhline(11.5, color="k", ls="--", lw=1)
    b.text(0.05, 11.9, "Ceron 11,5 %", fontsize=10, va="bottom")
    for x, p in zip([0, 1, 2, 3], pf):
        b.text(x + 0.12, p - 1.6, f"{p:.2f}".replace(".", ","), fontsize=10)
    b.set_xticks([0, 1, 2, 3])
    b.set_xlabel("níveis de refinamento (NREF)")
    b.set_ylabel(r"$P_f$ (%), PZ 'reservatorio'")
    b.set_ylim(-1, 19)
    b.set_xlim(-0.3, 3.6)
    virg(b)
    b.grid(axis="y", color="0.9")
    salva(fig, "convergencia")
    return e_res, e_cer


# --------------------------------------------------------------------------------------------------------------------
def copia_e_recorta():
    """as figuras do CeronCinematico: copiadas como estao e, para o texto, um painel recortado (fonte legivel)"""
    for n in ("fig4b_reproducao", "fig4c_reproducao", "pressao_permanente", "instantaneo_pz", "pz_contorno_ceron"):
        for ext in ("png", "pdf"):
            shutil.copy2(os.path.join(CC, "figs", n + "." + ext), os.path.join(FIGS, n + "." + ext))
    im = Image.open(os.path.join(FIGS, "fig4b_reproducao.png"))
    im.crop((0, 556, 900, 1080)).save(os.path.join(FIGS, "fig4b_painel_D.png"), dpi=(120, 120))
    im.crop((900, 28, 1800, 548)).save(os.path.join(FIGS, "fig4b_painel_C.png"), dpi=(120, 120))
    im = Image.open(os.path.join(FIGS, "fig4c_reproducao.png"))
    im.crop((0, 55, 960, 552)).save(os.path.join(FIGS, "fig4c_painel_esq.png"), dpi=(120, 120))
    print("figs/ copias e recortes")


def tabelas(e_res, e_cer, ma):
    f = lambda x, n=3: f"{x:.{n}f}".replace(".", "{,}")
    out = ["% gerado por figuras_contornos.py; nao edite a mao"]
    # dominio x contorno
    lin = le_csv(os.path.join(CC, "dados", "estudo_contorno.csv"))
    doms = []
    for r in lin:
        if r["dominio"] not in doms:
            doms.append(r["dominio"])
    nomes = {doms[0]: "notebook (3,5H; 3H; 3H)", doms[1]: "nosso domínio (7H; 6H; 6H)",
             doms[2]: "Fig.~4 aprox.\\ (13H; 2H; 8H)"}
    out.append(r"\newcommand{\tabdominio}{%")
    out.append(r"\begin{tabular}{lrrrr}\toprule")
    out.append(r"Domínio (atrás da crista; à frente do pé; abaixo do pé) & A & B & C & D \\ \midrule")
    for d in doms:
        g = [float(r["Gamma"]) for r in lin if r["dominio"] == d]
        out.append(nomes[d] + " & " + " & ".join(f(x) for x in g) + r" \\")
    out.append(r"\bottomrule\end{tabular}}")
    # instantaneo
    lin = le_csv(os.path.join(CC, "dados", "instantaneo_pz.csv"))
    rt = {"0.0001": "0+", "0.1": "0,1", "1": "1", "10": "10"}
    out.append(r"\newcommand{\tabinstantaneo}{%")
    out.append(r"\begin{tabular}{llrrr}\toprule")
    out.append(r"Água & $t$ (h) & $\Gamma$ GIM u-pw-pa & $\Gamma$ GIM u-p & $\Gamma$ cinemático \\ \midrule")
    for r in lin:
        out.append(("saturada" if r["agua"] == "sat" else "não saturada") + f" & {rt[r['t_h']]} & "
                   + f(float(r["Gamma_GIM_upwpa"])) + " & " + f(float(r["Gamma_GIM_up"])) + " & "
                   + f(float(r["Gamma_cin"])) + r" \\")
    out.append(r"\bottomrule\end{tabular}}")
    lin = le_csv(os.path.join(CC, "dados", "pz_contorno_ceron.csv"))
    out.append(r"\newcommand{\gcin}{" + f(float(lin[0]["Gamma_cin"])) + "}")
    out.append(r"\newcommand{\gextceron}{" + f(e_cer) + "}")
    out.append(r"\newcommand{\gextres}{" + f(e_res) + "}")
    gc, gd = niveis_ceron(), niveis_det()
    out.append(r"\newcommand{\gcerii}{" + f(gc[2]) + "}")
    out.append(r"\newcommand{\gceriii}{" + f(gc[3]) + "}")
    out.append(r"\newcommand{\gceriv}{" + f(gc[4]) + "}")
    out.append(r"\newcommand{\gcerv}{" + f(gc[5]) + "}")
    out.append(r"\newcommand{\gdetii}{" + f(gd[2]) + "}")
    out.append(r"\newcommand{\gpyfino}{" + f(ma[0]) + "}")
    out.append(r"\newcommand{\gpzufino}{" + f(ma[1]) + "}")
    out.append(r"\newcommand{\dufino}{" + f(ma[2], 2) + "}")
    open(os.path.join(AQUI, "tabelas.tex"), "w", encoding="utf-8").write("\n".join(out) + "\n")
    print("tabelas.tex")


if __name__ == "__main__":
    copia_e_recorta()
    esquema_ceron()
    esquema_pz(True)
    esquema_pz(False)
    esquema_rahardjo_codigo()
    esquema_rahardjo_plano()
    contornos_permanente()
    gamma_dominio()
    instantaneo()
    pz_ceron()
    e_res, e_cer = convergencia()
    tabelas(e_res, e_cer, malha_agua())
