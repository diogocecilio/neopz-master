#!/usr/bin/env python3
"""figuras_sat.py: as figuras da Secao 6 de Vargas Ceron et al. (2025) com as curvas do artigo (extraidas do PDF,
ceron_dig.json) e as deste trabalho (estatisticas de dados_sat.py) nos mesmos eixos: os do artigo, ampliados so quando
os dados nao cabem (eixos_ceron.amplia), e a legenda da figura diz qual era a faixa original.

Convencoes: cada valor do parametro tem uma cor (a paleta categorica de referencia, em ordem fixa); o Ceron e
tracejado e este trabalho e cheio. As densidades sao as log-normais ajustadas, como no artigo."""

import math
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(AQUI), "relatorio_ceron"))
from eixos_ceron import FIG, amplia  # noqa: E402

PALETA = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
TINTA, TINTA2, GRADE = "#0b0b0b", "#52514e", "#dddcd8"
CERON, NOSSO = (0, (5, 2.5)), "-"
plt.rcParams.update({
    "font.size": 8.5, "axes.labelsize": 8.5, "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "legend.fontsize": 7,
    "axes.edgecolor": TINTA2, "axes.labelcolor": TINTA, "xtick.color": TINTA2, "ytick.color": TINTA2,
    "axes.grid": True, "grid.color": GRADE, "grid.linewidth": 0.5, "lines.linewidth": 1.4, "legend.frameon": False,
    "savefig.dpi": 200, "font.family": "DejaVu Sans", "mathtext.fontset": "dejavusans",
})
ROT = {"x_gamma": r"$\Gamma$", "pdf": "Densidade de probabilidade", "pf": r"$P_f$", "mu": r"$\mu(\Gamma)$",
       "cov": r"$CoV(\Gamma)$"}
NOME = {"covk": r"$CoV(k_v)$", "covc": r"$CoV(c)$", "covphi": r"$CoV(\varphi)$", "s": r"$s$", "alfa": r"$\alpha$",
        "hw": r"$h_w/H$"}
# os paineis de cada familia: (figura da densidade, figuras de P_f, CoV e media)
PAINEIS = {"covk": ("27a", "28a", "30a", "29a"), "covc": ("27b", "28b", "30b", "29b"),
           "covphi": ("27c", "28c", "30c", "29c"), "s": ("31a", "31b", "31c", "31d"),
           "alfa": ("32a", "32b", "32c", "32d"), "hw": ("33a", "33b", "33c", "33d")}
VALOR_REF = {"covk": 0.6, "covc": 0.3, "covphi": 0.1, "s": 1.0, "alfa": 1.0, "hw": 1.0}


def lognormal_pdf(x, mu_ln, sigma_ln):
    x = np.asarray(x, float)
    y = np.zeros_like(x)
    m = x > 0
    y[m] = np.exp(-(np.log(x[m]) - mu_ln) ** 2 / (2 * sigma_ln ** 2)) / (x[m] * sigma_ln * math.sqrt(2 * math.pi))
    return y


def _fmt(v):
    return ("%g" % v).replace(".", ",")


class Ampliacoes:
    """registra os eixos ampliados de cada figura, para a legenda"""

    def __init__(self):
        self.notas = {}

    def eixo(self, ax, fig_id, eixo, lim, dados, log=False):
        a, b, p, amp = amplia(lim, dados)
        if amp:
            self.notas.setdefault(fig_id[:2], []).append(
                "(%s) eixo %s ampliado de [%s; %s] para [%s; %s]" % (fig_id[2:] or fig_id, eixo, _fmt(lim[0]), _fmt(lim[1]),
                                                                    _fmt(a), _fmt(b)))
        if eixo == "x":
            ax.set_xlim(a, b)
            if p and not log:
                ax.set_xticks(np.arange(a, b + p / 2, p))
        else:
            ax.set_ylim(a, b)
            if p and not log:
                ax.set_yticks(np.arange(a, b + p / 2, p))


def casos_familia(est, fam):
    """os casos deste trabalho de uma familia, com a referencia no seu valor"""
    out = {}
    for c in est["casos"]:
        if c["familia"] == "ref":
            out[VALOR_REF[fam]] = c
        elif c["familia"] == fam:
            out[float(c["valor"])] = c
    return out


def cores_familia(fam, valores):
    vs = sorted(set(float(v) for v in valores))
    return {v: PALETA[i % len(PALETA)] for i, v in enumerate(vs)}


def _legenda_estilos(ax, extra=()):
    h = [Line2D([], [], color=TINTA2, ls=CERON, lw=1.4), Line2D([], [], color=TINTA2, ls=NOSSO, lw=1.4)]
    r = ["Ceron et al. (2025)", "este trabalho"]
    for hh, rr in extra:
        h.append(hh)
        r.append(rr)
    return h, r


def painel_pdf(ax, fam, fig_id, nossos, ceron_pdfs, amp):
    """densidades: Ceron (curvas extraidas) tracejado, este trabalho (log-normal ajustada) cheio"""
    valores = set(nossos) | set(float(v) for v in ceron_pdfs)
    cores = cores_familia(fam, valores)
    xs_lo, xs_hi, ys_hi = [], [], []
    def visivel(x, y):  # a parte da curva acima de 2% do pico decide a ampliacao do eixo x
        m = y >= 0.02 * y.max()
        xs_lo.append(float(x[m].min()))
        xs_hi.append(float(x[m].max()))
    for v, cur in sorted(ceron_pdfs.items(), key=lambda t: float(t[0])):
        x, y = np.asarray(cur[0]), np.asarray(cur[1])
        ax.plot(x, y, color=cores[float(v)], ls=CERON, lw=1.2)
        ys_hi.append(y.max())
        visivel(x, y)
    for v, c in sorted(nossos.items()):
        e = c.get("estat")
        if not e:
            continue
        lo, hi = math.exp(e["mu_ln"] - 4.5 * e["sigma_ln"]), math.exp(e["mu_ln"] + 4.5 * e["sigma_ln"])
        x = np.linspace(lo, hi, 400)
        y = lognormal_pdf(x, e["mu_ln"], e["sigma_ln"])
        ax.plot(x, y, color=cores[v], ls=NOSSO, lw=1.6)
        visivel(x, y)
        ys_hi.append(y.max())
    lim = FIG[fig_id]
    amp.eixo(ax, fig_id, "x", lim["x"], (min(xs_lo, default=None), max(xs_hi, default=None)))
    amp.eixo(ax, fig_id, "y", lim["y"], (0.0, max(ys_hi, default=None)))
    ax.set_xlabel(ROT["x_gamma"])
    ax.set_ylabel(ROT["pdf"])
    h, r = _legenda_estilos(ax)
    l1 = ax.legend(h, r, loc="upper right", handlelength=2.2, fontsize=6.5)
    ax.add_artist(l1)
    hv = [Line2D([], [], color=cores[v], lw=2.2) for v in sorted(cores)]
    rv = ["%s%s" % (_fmt(v), " (ref.)" if v == VALOR_REF[fam] else "") for v in sorted(cores)]
    ax.legend(hv, rv, title=NOME[fam], title_fontsize=6.8, fontsize=6.5, loc="center right", handlelength=1.4,
              bbox_to_anchor=(1.0, 0.45))


def painel_escalar(ax, fam, fig_id, chave, nossos, ceron_tab, amp, marca_ref=True):
    """P_f, CoV ou media contra o parametro: Ceron (tabela do artigo) e este trabalho; s = 400 como linha horizontal"""
    lim = FIG[fig_id]
    xmin_eixo, xmax_eixo = lim["x"][0], lim["x"][1]

    def fora(x):
        return x < xmin_eixo - 1e-12 or x > xmax_eixo + 1e-12

    def valor(e):
        return {"pf": e["P_f"], "cov": e["CoV"], "mu": e["media"]}[chave]
    # Ceron
    cx, cy, fora_c = [], [], []
    for row in ceron_tab:
        x, y = row["x"], row[chave]
        if fora(x):
            fora_c.append((x, y))
        else:
            cx.append(x)
            cy.append(y)
    o = np.argsort(cx)
    ax.plot(np.array(cx)[o], np.array(cy)[o], color=TINTA2, ls=CERON, lw=1.2, marker="s", ms=4.5, mfc="white",
            mec=TINTA2, mew=1.0)
    # este trabalho
    nx, ny, ne, fora_n = [], [], [], []
    for v, c in sorted(nossos.items()):
        e = c.get("estat")
        if not e:
            continue
        y = valor(e)
        err = 1.96 * math.sqrt(e["P_f"] * (1 - e["P_f"]) / e["S"]) if chave == "pf" else 0.0
        if fora(v):
            fora_n.append((v, y, err))
        else:
            nx.append(v)
            ny.append(y)
            ne.append(err)
    ax.errorbar(nx, ny, yerr=ne if chave == "pf" else None, color=PALETA[0], ls=NOSSO, lw=1.6, marker="o", ms=4.5,
                capsize=2.5, elinewidth=0.9)
    extra = []
    for (x, y) in fora_c:
        ax.axhline(y, color=TINTA2, ls=(0, (1, 1.5)), lw=1.0)
        extra.append((Line2D([], [], color=TINTA2, ls=(0, (1, 1.5))), "Ceron, %s = %s" % (NOME[fam], _fmt(x))))
    for (x, y, err) in fora_n:
        ax.axhline(y, color=PALETA[0], ls=(0, (1, 1.5)), lw=1.2)
        extra.append((Line2D([], [], color=PALETA[0], ls=(0, (1, 1.5))), "este trabalho, %s = %s" % (NOME[fam], _fmt(x))))
    todos_y = cy + ny + [y for _, y in fora_c] + [y for _, y, _ in fora_n] + \
        [y + e for y, e in zip(ny, ne)] + [y - e for y, e in zip(ny, ne)]
    amp.eixo(ax, fig_id, "x", lim["x"], (min(cx + nx, default=None), max(cx + nx, default=None)))
    amp.eixo(ax, fig_id, "y", lim["y"], (min(todos_y, default=None), max(todos_y, default=None)))
    ax.set_xlabel(NOME[fam])
    ax.set_ylabel(ROT[chave])
    h = [Line2D([], [], color=TINTA2, ls=CERON, marker="s", mfc="white", mec=TINTA2),
         Line2D([], [], color=PALETA[0], ls=NOSSO, marker="o")]
    r = ["Ceron et al. (2025)", "este trabalho" + (" (IC 95%)" if chave == "pf" else "")]
    for hh, rr in extra:
        h.append(hh)
        r.append(rr)
    ax.legend(h, r, loc="best")


def letra(ax, t):
    ax.text(-0.02, 1.02, "(%s)" % t, transform=ax.transAxes, ha="right", va="bottom", fontsize=9, color=TINTA)


def figura_familia_pdf_escalares(fam, est, cer, saida, amp):
    """Figs. 31 a 33: (a) densidades, (b) P_f, (c) CoV, (d) media"""
    f_pdf, f_pf, f_cov, f_mu = PAINEIS[fam]
    nossos = casos_familia(est, fam)
    fig, axs = plt.subplots(2, 2, figsize=(7.2, 5.6), constrained_layout=True)
    painel_pdf(axs[0, 0], fam, f_pdf, nossos, cer["pdfs"].get(fam, {}), amp)
    for ax, fid, ch in ((axs[0, 1], f_pf, "pf"), (axs[1, 0], f_cov, "cov"), (axs[1, 1], f_mu, "mu")):
        painel_escalar(ax, fam, fid, ch, nossos, cer["tab"][fam], amp)
    for ax, t in zip(axs.flat, "abcd"):
        letra(ax, t)
    arq = os.path.join(saida, "fig%s.pdf" % f_pdf[:2])
    fig.savefig(arq)
    plt.close(fig)
    return arq


def figura27(est, cer, saida, amp):
    fig, axs = plt.subplots(1, 3, figsize=(10.5, 3.3), constrained_layout=True)
    for ax, fam, t in zip(axs, ("covk", "covc", "covphi"), "abc"):
        painel_pdf(ax, fam, PAINEIS[fam][0], casos_familia(est, fam), cer["pdfs"].get(fam, {}), amp)
        letra(ax, t)
    arq = os.path.join(saida, "fig27.pdf")
    fig.savefig(arq)
    plt.close(fig)
    return arq


def figura_escalar_cov(n, chave, est, cer, saida, amp):
    """Figs. 28 (P_f), 29 (media) e 30 (CoV) contra o CoV de k_v, c e phi"""
    idx = {"pf": 1, "mu": 3, "cov": 2}[chave]
    fig, axs = plt.subplots(1, 3, figsize=(10.5, 3.1), constrained_layout=True)
    for ax, fam, t in zip(axs, ("covk", "covc", "covphi"), "abc"):
        painel_escalar(ax, fam, PAINEIS[fam][idx], chave, casos_familia(est, fam), cer["tab"][fam], amp)
        letra(ax, t)
    arq = os.path.join(saida, "fig%d.pdf" % n)
    fig.savefig(arq)
    plt.close(fig)
    return arq


def figura25(est, cer, saida, amp):
    ref = next(c for c in est["casos"] if c["familia"] == "ref")
    e = ref["estat"]
    g = np.asarray(ref["gamma"])
    fig, axs = plt.subplots(1, 2, figsize=(8.6, 3.3), constrained_layout=True)
    ax = axs[0]
    h25 = cer.get("fig25a", {})
    if h25.get("hist"):
        b, a = np.asarray(h25["hist"]["bordas"]), np.asarray(h25["hist"]["alturas"])
        ax.stairs(a, b, color=TINTA2, lw=0.9, ls=CERON)
    if h25.get("pdf"):
        ax.plot(h25["pdf"][0], h25["pdf"][1], color=TINTA2, ls=CERON, lw=1.4)
    ax.hist(g, bins="auto", density=True, color=PALETA[0], alpha=0.28, edgecolor="white", linewidth=0.4)
    x = np.linspace(max(1e-3, g.min() * 0.7), g.max() * 1.2, 500)
    ax.plot(x, lognormal_pdf(x, e["mu_ln"], e["sigma_ln"]), color=PALETA[0], lw=1.6)
    amp.eixo(ax, "25a", "x", FIG["25a"]["x"], (g.min(), g.max()))
    amp.eixo(ax, "25a", "y", FIG["25a"]["y"], (0, max(lognormal_pdf(x, e["mu_ln"], e["sigma_ln"]).max(),
                                                       np.histogram(g, bins="auto", density=True)[0].max())))
    ax.axvline(1.0, color=TINTA2, lw=0.8)
    ax.set_xlabel(ROT["x_gamma"])
    ax.set_ylabel(ROT["pdf"])
    ax.legend([Patch(facecolor="none", edgecolor=TINTA2, ls="--"), Line2D([], [], color=TINTA2, ls=CERON),
               Patch(color=PALETA[0], alpha=0.28), Line2D([], [], color=PALETA[0])],
              ["Ceron: histograma", "Ceron: log-normal", "este trabalho: histograma", "este trabalho: log-normal"],
              loc="upper right")
    letra(ax, "a")
    ax = axs[1]
    if cer.get("fig25b"):
        ax.plot(cer["fig25b"][0], cer["fig25b"][1], color=TINTA2, ls=CERON, lw=1.2)
    n = np.arange(1, len(ref["conv_pf"]) + 1)
    ax.plot(n, ref["conv_pf"], color=PALETA[0], lw=1.4)
    ax.set_xscale("log")
    ax.set_xlim(*FIG["25b"]["x"][:2])
    amp.eixo(ax, "25b", "y", FIG["25b"]["y"], (0.0, max(ref["conv_pf"])))
    ax.set_xlabel("Número de simulações")
    ax.set_ylabel(ROT["pf"])
    h, r = _legenda_estilos(ax)
    ax.legend(h, r, loc="upper right")
    letra(ax, "b")
    arq = os.path.join(saida, "fig25.pdf")
    fig.savefig(arq)
    plt.close(fig)
    return arq


def figura26(est, cer, saida):
    ref = next(c for c in est["casos"] if c["familia"] == "ref")
    nomes = [("k", r"$k_v$"), ("phi", r"$\varphi$"), ("c", r"$c$")]
    fig, ax = plt.subplots(figsize=(5.4, 2.9), constrained_layout=True)
    y0 = np.arange(len(nomes))[::-1]
    larg = 0.26
    series = [("Ceron et al. (2025)", {k: cer["fig26"].get(k) for k, _ in nomes}, None, TINTA2),
              ("este trabalho: banda e massa", {k: ref["spearman_banda"][k]["rho_parcial"] for k, _ in nomes},
               {k: ref["spearman_banda"][k]["ic95"] for k, _ in nomes}, PALETA[0]),
              ("este trabalho: médias do campo", {k: ref["spearman_campo"][k]["rho_parcial"] for k, _ in nomes},
               {k: ref["spearman_campo"][k]["ic95"] for k, _ in nomes}, PALETA[1])]
    for i, (rot, vals, ics, cor) in enumerate(series):
        ys = y0 + (1 - i) * larg
        v = [vals[k] if vals[k] is not None else np.nan for k, _ in nomes]
        ax.barh(ys, v, height=larg * 0.9, color=cor, label=rot)
        if ics:
            lo = [v[j] - ics[k][0] for j, (k, _) in enumerate(nomes)]
            hi = [ics[k][1] - v[j] for j, (k, _) in enumerate(nomes)]
            ax.errorbar(v, ys, xerr=[lo, hi], fmt="none", ecolor=TINTA, elinewidth=0.8, capsize=2)
    ax.set_yticks(y0)
    ax.set_yticklabels([n for _, n in nomes])
    ax.set_xlim(*FIG["26"]["x"][:2])
    ax.set_xticks(np.arange(-1.0, 1.01, 0.25))
    ax.axvline(0, color=TINTA2, lw=0.8)
    ax.set_xlabel(r"Spearman parcial de 2ª ordem, $\rho_{X\Gamma,YZ}$")
    ax.legend(loc="lower left")
    arq = os.path.join(saida, "fig26.pdf")
    fig.savefig(arq)
    plt.close(fig)
    return arq


def figura24(est, cer, espiral_nossa, saida):
    ref = next(c for c in est["casos"] if c["familia"] == "ref")
    det = ref.get("det") or {}
    fig, ax = plt.subplots(figsize=(7.0, 3.4), constrained_layout=True)
    ax.plot([0, 25, 25, 15, 10, 0, 0], [0, 0, 5, 5, 10, 10, 0], color=TINTA2, lw=1.0)
    if det.get("banda_pontos"):
        p = np.asarray(det["banda_pontos"])
        ax.scatter(p[:, 0], p[:, 1], s=3, color=PALETA[2], alpha=0.6, lw=0,
                   label=r"banda de ruptura do GIM ($\sqrt{J_2} \geq 0{,}1$ do máximo), $\Gamma = %s$" % _fmt(
                       round(det["Gamma"], 3)))
    if cer.get("fig24"):
        ax.plot(cer["fig24"][0], cer["fig24"][1], color=TINTA, ls=CERON, lw=1.4,
                label=r"Ceron et al. (2025), $\Gamma = 1{,}336$")
    if espiral_nossa is not None:
        ax.plot(espiral_nossa["x"], espiral_nossa["y"], color=PALETA[1], lw=1.4,
                label=r"método do Ceron com a pressão deste trabalho, $\Gamma = %s$" % _fmt(round(espiral_nossa["G"], 3)))
    ax.set_aspect("equal")
    ax.set_xlim(0, 25)
    ax.set_ylim(0, 12.5)
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax.legend(loc="upper right", fontsize=6.8)
    arq = os.path.join(saida, "fig24.pdf")
    fig.savefig(arq)
    plt.close(fig)
    return arq


def tem_dados(est, fam):
    """se ha algum caso da familia (alem da referencia) com estatisticas"""
    return any(c["familia"] == fam and c.get("estat") for c in est["casos"])


def todas(est, cer, espiral_nossa, saida):
    """as figuras das familias com dados; as outras ficam de fora (o relatorio diz quais)"""
    os.makedirs(saida, exist_ok=True)
    amp = Ampliacoes()
    arqs = {"24": figura24(est, cer, espiral_nossa, saida), "25": figura25(est, cer, saida, amp),
            "26": figura26(est, cer, saida)}
    if any(tem_dados(est, f) for f in ("covk", "covc", "covphi")):
        arqs["27"] = figura27(est, cer, saida, amp)
        for n, ch in ((28, "pf"), (29, "mu"), (30, "cov")):
            arqs[str(n)] = figura_escalar_cov(n, ch, est, cer, saida, amp)
    for fam in ("s", "alfa", "hw"):
        if tem_dados(est, fam):
            arqs[PAINEIS[fam][0][:2]] = figura_familia_pdf_escalares(fam, est, cer, saida, amp)
    return arqs, amp.notas
