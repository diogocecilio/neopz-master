#!/usr/bin/env python3
"""figuras_chuva.py: as figuras do estudo da chuva (Monte Carlo do talude de Rahardjo et al. 2007 com os campos
aleatorios de Vargas Ceron et al. 2025), analogas as Figs. 24 a 33 da Secao 6 do Ceron que o figuras_sat.py do
rebaixamento desenha (SlopeRandomField/comparacao_ceron), mais as da chuva. Le o JSON do dados_chuva.py.

Figuras (cada uma so quando os dados existem; as outras ficam listadas como puladas, com o motivo):
  fig24_det        talude, banda de ruptura do deterministico antes da chuva e no fim dela, nivel d'agua (p_w = 0)
  fig25_ref        (a) histograma de FS, log-normal ajustada, FS deterministico; (b) convergencia de P_f
  fig26_spearman   Spearman parcial de 2a ordem (banda e massa; medias do campo), antes da chuva e no fim dela
  fig27_cov_pdf    densidades log-normais de FS com o CoV de k_s, c' e phi' (um painel por familia com dados)
  fig28_cov_pf, fig29_cov_mu, fig30_cov_cov   P_f, media e CoV de FS contra o CoV
  fig31_s, fig32_alfa, fig33_I   (a) densidades, (b) P_f, (c) CoV, (d) media de FS; em (d), o FS deterministico (o de
                   cada valor em alfa e I) e, em I, a sigmoide F_min(I) da Fig. 7b do Rahardjo
  fig_duracao      FS e P_f com a duracao da chuva (os instantes da referencia: o FS em t e o do fim de uma chuva de t h)
  fig_fs0_fs24     FS antes da chuva x FS no fim dela, por realizacao, com a categoria do nivel d'agua
  fig_mapa         frequencia da banda do mecanismo em cada celula de 1 m: em todas as realizacoes e so nas rupturas
  fig_nivel        frequencia das categorias do nivel d'agua (Tabela 3 do Rahardjo): com a duracao e em cada caso

Convencoes (as do figuras_sat.py): a paleta categorica de referencia em ordem fixa, um valor do parametro por cor (em
ordem crescente); o valor de referencia leva "(ref.)" na legenda e um anel nos graficos de escalares; as densidades
sao as log-normais ajustadas, como no Ceron. Cada figura sai em PDF (para o LaTeX) e em PNG.

Uso: python3 -I figuras_chuva.py <estatisticas.json> <pasta>"""

import json
import math
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.ticker  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
from matplotlib.tri import Triangulation  # noqa: E402

PALETA = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
TINTA, TINTA2, MUDO, GRADE, SOLO = "#0b0b0b", "#52514e", "#898781", "#dddcd8", "#f3f2ee"
# rampa sequencial azul (100 a 700) para os mapas; a ordinal das categorias do nivel d'agua (mais agua, mais escuro)
RAMPA = ["#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7", "#3987e5", "#2a78d6", "#256abf",
         "#1c5cab", "#184f95", "#104281", "#0d366b"]
CORES_CATEGORIA = {"I": "#0d366b", "II": "#1c5cab", "III": "#3987e5", "IV": "#86b6ef"}
TEXTO_CATEGORIA = {"I": "I: até a crista", "II": "II: até o meio da face", "III": "III: até o pé",
                   "IV": "IV: sem subida ou subida pequena"}
CATEGORIAS = ("I", "II", "III", "IV")
plt.rcParams.update({
    "font.size": 8.5, "axes.labelsize": 8.5, "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "legend.fontsize": 7,
    "axes.edgecolor": TINTA2, "axes.labelcolor": TINTA, "xtick.color": TINTA2, "ytick.color": TINTA2,
    "axes.grid": True, "grid.color": GRADE, "grid.linewidth": 0.5, "lines.linewidth": 1.4, "legend.frameon": False,
    "savefig.dpi": 200, "font.family": "DejaVu Sans", "mathtext.fontset": "dejavusans",
})
ROT = {"x_fs": r"$FS$", "pdf": "Densidade de probabilidade", "pf": r"$P_f$", "mu": r"$\mu(FS)$",
       "cov": r"$CoV(FS)$"}
NOME = {"covk": r"$CoV(k_s)$", "covc": r"$CoV(c')$", "covphi": r"$CoV(\varphi')$", "s": r"$s$",
        "alfa": r"$\alpha = k_h/k_v$", "I": r"$I$ (mm/h)"}
NOME_CURTO = {"covk": r"$CoV(k_s)$", "covc": r"$CoV(c')$", "covphi": r"$CoV(\varphi')$", "s": r"$s$",
              "alfa": r"$\alpha$", "I": r"$I$"}
UNIDADE = {"I": " mm/h"}
FAMILIAS = ("covk", "covc", "covphi", "s", "alfa", "I")
DET_PROPRIO = ("alfa", "I")
S_EIXO_MAX = 20.0  # s: o eixo do Ceron vai ate 20; os valores maiores (s = 400) entram como linha horizontal
MAX_TRAJETORIAS = 300  # fig_duracao: no maximo tantas realizacoes desenhadas


def lognormal_pdf(x, mu_ln, sigma_ln):
    x = np.asarray(x, float)
    y = np.zeros_like(x)
    m = x > 0
    y[m] = np.exp(-(np.log(x[m]) - mu_ln) ** 2 / (2 * sigma_ln ** 2)) / (x[m] * sigma_ln * math.sqrt(2 * math.pi))
    return y


def rahardjo_fmin(I, rah):
    s = rah["sigmoide"]
    return s["a"] + s["b"] / (1. + np.exp(-(np.asarray(I, float) - s["c"]) / s["d"]))


def _fmt(v, casas=None):
    if v is None:
        return "-"
    t = ("%g" % v) if casas is None else ("%.*f" % (casas, v))
    return t.replace(".", ",")


def letra(ax, t):
    ax.text(-0.02, 1.02, "(%s)" % t, transform=ax.transAxes, ha="right", va="bottom", fontsize=9, color=TINTA)


def rotulo_valor(fam, v):
    return "%s = %s%s" % (NOME_CURTO[fam], _fmt(v), UNIDADE.get(fam, ""))


def ref_de(est):
    return next((c for c in est["casos"] if c["familia"] == "ref"), None)


def casos_familia(est, fam):
    """os casos de uma familia, com a referencia no seu valor (um caso explicito com o mesmo valor tem preferencia)"""
    out = {}
    vref = est["valor_ref"][fam]
    for c in est["casos"]:
        if c["familia"] == fam and c.get("estat"):
            out[float(c["valor"])] = c
    ref = ref_de(est)
    if ref and ref.get("estat") and not any(abs(v - vref) < 1e-9 for v in out):
        out[float(vref)] = ref
    return out


def tem_dados(est, fam):
    """se ha algum caso da familia (alem da referencia) com estatisticas"""
    return any(c["familia"] == fam and c.get("estat") for c in est["casos"])


def cores_familia(valores):
    vs = sorted(set(float(v) for v in valores))
    return {v: PALETA[i % len(PALETA)] for i, v in enumerate(vs)}


def instante(caso, t):
    return next((i for i in caso.get("instantes", []) if abs(i["t"] - t) < 1e-6), None)


def salva(fig, pasta, nome):
    pdf = os.path.join(pasta, nome + ".pdf")
    fig.savefig(pdf)
    fig.savefig(os.path.join(pasta, nome + ".png"), dpi=150)
    plt.close(fig)
    return pdf


def contorno_talude(ax, talude, preenche=True):
    p = np.asarray(talude + [talude[0]], float)
    if preenche:
        ax.fill(p[:, 0], p[:, 1], color=SOLO, lw=0, zorder=0)
    ax.plot(p[:, 0], p[:, 1], color=TINTA2, lw=1.0, zorder=3)


# ---------------------------------------------------------------- paineis das familias
def painel_pdf(ax, fam, nossos, vref):
    """densidades log-normais ajustadas, uma cor por valor (crescente)"""
    cores = cores_familia(nossos)
    lo_vis, hi_vis, ymax = [], [], 0.
    for v, c in sorted(nossos.items()):
        e = c["estat"]
        if not e or not e["sigma_ln"] or e["sigma_ln"] <= 0:
            continue
        x = np.linspace(math.exp(e["mu_ln"] - 4.5 * e["sigma_ln"]), math.exp(e["mu_ln"] + 4.5 * e["sigma_ln"]), 500)
        y = lognormal_pdf(x, e["mu_ln"], e["sigma_ln"])
        ax.plot(x, y, color=cores[v], lw=2.0 if abs(v - vref) < 1e-9 else 1.4)
        m = y >= 0.02 * y.max()
        lo_vis.append(x[m].min())
        hi_vis.append(x[m].max())
        ymax = max(ymax, y.max())
    if lo_vis:
        a, b = min(lo_vis), max(hi_vis)
        ax.set_xlim(min(a, 0.95) - 0.02 * (b - a), b + 0.02 * (b - a))
        ax.set_ylim(0, 1.08 * ymax)
    ax.axvline(1.0, color=TINTA2, lw=0.8)
    ax.set_xlabel(ROT["x_fs"])
    ax.set_ylabel(ROT["pdf"])
    hv = [Line2D([], [], color=cores[v], lw=2.2) for v in sorted(cores)]
    rv = ["%s%s" % (_fmt(v), " (ref.)" if abs(v - vref) < 1e-9 else "") for v in sorted(cores)]
    ax.legend(hv, rv, title=NOME[fam], title_fontsize=6.8, fontsize=6.5, loc="upper right", handlelength=1.4)


def valor_escalar(e, chave):
    return {"pf": e["P_f"], "cov": e["CoV"], "mu": e["media"]}[chave]


def painel_escalar(ax, fam, chave, nossos, vref, est):
    """P_f (com o IC de 95% de Wilson e o P_f do ajuste log-normal), CoV ou media de FS contra o parametro; o valor de
    referencia com um anel; em s, os valores alem de S_EIXO_MAX como linhas horizontais"""
    xs, ys, lo, hi, pln, fora = [], [], [], [], [], []
    for v, c in sorted(nossos.items()):
        e = c["estat"]
        y = valor_escalar(e, chave)
        if fam == "s" and v > S_EIXO_MAX + 1e-9:
            fora.append((v, y))
            continue
        xs.append(v)
        ys.append(y)
        if chave == "pf":
            a, b = e["IC95_Pf"]
            lo.append(y - a)
            hi.append(b - y)
            pln.append(e["P_f_lognormal"])
    xs, ys = np.array(xs), np.array(ys)
    h, r = [], []
    if chave == "pf":
        ax.errorbar(xs, ys, yerr=[lo, hi], color=PALETA[0], lw=1.6, marker="o", ms=4.5, capsize=2.5, elinewidth=0.9)
        h.append(Line2D([], [], color=PALETA[0], marker="o"))
        r.append("contagem (IC 95% de Wilson)")
        ax.plot(xs, pln, color=PALETA[0], lw=1.0, ls=(0, (4, 2)), marker="o", ms=3.5, mfc="white", mew=0.9)
        h.append(Line2D([], [], color=PALETA[0], ls=(0, (4, 2)), marker="o", mfc="white", ms=3.5))
        r.append("ajuste log-normal")
    else:
        ax.plot(xs, ys, color=PALETA[0], lw=1.6, marker="o", ms=4.5)
        h.append(Line2D([], [], color=PALETA[0], marker="o"))
        r.append(r"$\mu$ do ajuste log-normal" if chave == "mu" else r"$CoV$ do ajuste log-normal")
    yref = [y for x, y in zip(xs, ys) if abs(x - vref) < 1e-9]
    if yref:
        ax.plot([vref], yref, ls="none", marker="o", ms=9, mfc="none", mec=TINTA, mew=1.1)
        h.append(Line2D([], [], ls="none", marker="o", ms=8, mfc="none", mec=TINTA, mew=1.1))
        r.append("referência")
    extra_y = []
    if chave == "mu":
        # FS deterministico: o de cada valor (alfa, I) ou o da referencia (as mesmas medias)
        if fam in DET_PROPRIO:
            dx = [v for v, c in sorted(nossos.items()) if (c.get("det") or {}).get("FS")]
            dy = [nossos[v]["det"]["FS"] for v in dx]
            if dx:
                ax.plot(dx, dy, color=PALETA[1], lw=1.3, marker="s", ms=4.2)
                h.append(Line2D([], [], color=PALETA[1], marker="s", ms=4.2))
                r.append("FS determinístico")
                extra_y += dy
        else:
            ref = ref_de(est)
            fd = ((ref or {}).get("det") or {}).get("FS")
            if fd:
                ax.axhline(fd, color=PALETA[1], lw=1.0, ls=(0, (1, 1.5)))
                h.append(Line2D([], [], color=PALETA[1], ls=(0, (1, 1.5))))
                r.append("FS determinístico (%s)" % _fmt(round(fd, 3)))
                extra_y.append(fd)
        if fam == "I" and len(xs):
            a, b = xs.min(), xs.max()
            pad = 0.08 * (b - a if b > a else 10.)
            xi = np.linspace(a - pad, b + pad, 200)
            yi = rahardjo_fmin(xi, est["rahardjo"])
            ax.plot(xi, yi, color=TINTA2, lw=1.2, ls=(0, (5, 2.5)))
            h.append(Line2D([], [], color=TINTA2, ls=(0, (5, 2.5))))
            r.append(r"Rahardjo et al. (2007), Fig. 7b: $F_{min}(I)$")
            extra_y += [yi.min(), yi.max()]
    for v, y in fora:
        ax.axhline(y, color=PALETA[0], ls=(0, (1, 1.5)), lw=1.2)
        h.append(Line2D([], [], color=PALETA[0], ls=(0, (1, 1.5))))
        r.append("%s = %s" % (NOME_CURTO[fam], _fmt(v)))
        extra_y.append(y)
    # eixos
    todos = list(ys) + extra_y + ([y + e for y, e in zip(ys, hi)] if chave == "pf" else []) + \
        (list(pln) if chave == "pf" else [])
    todos = [y for y in todos if y is not None and math.isfinite(y)]
    if chave == "pf":
        ax.set_ylim(0, max(0.01, 1.12 * max(todos, default=0.)))
    elif todos:
        a, b = min(todos), max(todos)
        pad = 0.12 * (b - a) if b > a else 0.05 * abs(b) + 1e-3
        ax.set_ylim(a - pad, b + pad)
    if fam == "s":  # como no Ceron, de 0 ate 20 (ou ate o maior s abaixo disso)
        ax.set_xlim(0, 1.05 * min(S_EIXO_MAX, xs.max() if len(xs) else S_EIXO_MAX) + 0.1)
    elif len(xs):
        a, b = xs.min(), xs.max()
        pad = 0.06 * (b - a) if b > a else 0.1 * abs(b) + 0.1
        ax.set_xlim(a - pad, b + pad)
    ax.set_xlabel(NOME[fam])
    ax.set_ylabel(ROT[chave])
    ax.legend(h, r, loc="best", fontsize=6.5)


# ---------------------------------------------------------------- figuras
def figura_det(est, pasta):
    ref = ref_de(est)
    det = (ref or {}).get("det") or {}
    if not det.get("bandas") and not det.get("malha_pw"):
        return None, "sem mc_real-1 na referencia (rodar com figuras=-1)"
    t1 = ref["t"]
    tempos = [0.0, t1] if abs(t1) > 1e-6 else [t1]
    cores = {0.0: PALETA[1], t1: PALETA[0]}
    fs_t = {round(a, 6): b for a, b in det.get("FS_t", [])}
    fig, ax = plt.subplots(figsize=(7.2, 3.6), constrained_layout=True)
    contorno_talude(ax, est["talude"])
    h, r = [], []
    for t in tempos:
        banda = next((b for b in det.get("bandas", []) if abs(b["t"] - t) < 1e-6), None)
        nome_t = "antes da chuva" if abs(t) < 1e-6 else "fim da chuva" if abs(t - t1) < 1e-6 else ""
        if banda and banda["pontos"]:
            p = np.asarray(banda["pontos"])
            ax.scatter(p[:, 0], p[:, 1], s=2.2, color=cores[t], alpha=0.7, lw=0, zorder=4)
            h.append(Line2D([], [], ls="none", marker="o", ms=4, color=cores[t]))
            r.append(r"banda, %s h (%s): $\sqrt{J_2(\varepsilon_p)} \geq 0{,}1$ do máximo, $FS$ = %s" % (
                _fmt(t), nome_t, _fmt(round(fs_t.get(round(t, 6), float("nan")), 3))))
        m = det.get("malha_pw")
        snap = next((s for s in (m or {}).get("instantes", []) if abs(s["t"] - t) < 1e-6), None)
        if snap:
            tri = Triangulation(m["x"], m["y"], m["tri"])
            ax.tricontour(tri, snap["pw"], levels=[0.0], colors=[cores[t]], linewidths=1.3, linestyles="--",
                          zorder=5)
            h.append(Line2D([], [], color=cores[t], ls="--", lw=1.3))
            r.append(r"$p_w = 0$, %s h" % _fmt(t))
    # as verticais do nivel d'agua (categorias da Tabela 3 do Rahardjo) e os niveis do deterministico nelas
    for k, col, cota in est["verticais"]:
        x = float(k[1:])
        ax.plot([x, x], [0, cota], color=MUDO, lw=0.6, ls=(0, (1, 2)), zorder=2)
    for t in tempos:
        nv = next((n for tt, n, _ in det.get("niveis_t", []) if abs(tt - t) < 1e-6), None)
        if nv:
            xs = [float(k[1:]) for k, _, _ in est["verticais"]]
            ys = [nv[k] for k, _, _ in est["verticais"]]
            ax.plot(xs, ys, ls="none", marker="_", ms=9, mew=2.0, color=cores[t], zorder=6)
    h.append(Line2D([], [], ls="none", marker="_", ms=9, mew=2.0, color=TINTA2))
    r.append("nível d'água nas verticais x = 30, 35 e 40 m")
    ax.set_aspect("equal")
    ax.set_xlim(0, 70)
    ax.set_ylim(0, 37)
    ax.set_yticks(np.arange(0, 31, 10))
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax.legend(h, r, loc="upper right", fontsize=6.3, ncol=1)
    return salva(fig, pasta, "fig24_det"), None


def figura_ref(est, pasta):
    ref = ref_de(est)
    if not ref or not ref.get("estat"):
        return None, "sem a referencia"
    e = ref["estat"]
    g = np.asarray(ref["fs"])
    fd = (ref.get("det") or {}).get("FS")
    fig, axs = plt.subplots(1, 2, figsize=(8.6, 3.3), constrained_layout=True)
    ax = axs[0]
    ax.hist(g, bins="auto", density=True, color=PALETA[0], alpha=0.28, edgecolor="white", linewidth=0.4)
    a, b = min(g.min(), fd or g.min(), 1.0), max(g.max(), fd or g.max())
    x = np.linspace(max(1e-3, a - 0.1 * (b - a)), b + 0.15 * (b - a), 500)
    ax.plot(x, lognormal_pdf(x, e["mu_ln"], e["sigma_ln"]), color=PALETA[0], lw=1.6)
    ax.axvline(1.0, color=TINTA2, lw=0.8)
    h = [Patch(color=PALETA[0], alpha=0.28), Line2D([], [], color=PALETA[0])]
    r = ["histograma (S = %d)" % e["S"], "log-normal ajustada"]
    if fd:
        ax.axvline(fd, color=PALETA[1], lw=1.3, ls=(0, (4, 2)))
        h.append(Line2D([], [], color=PALETA[1], ls=(0, (4, 2))))
        r.append("FS determinístico = %s" % _fmt(round(fd, 3)))
    ax.set_xlim(x[0], x[-1])
    ax.set_xlabel(ROT["x_fs"] + " em t = %s h" % _fmt(ref["t"]))
    ax.set_ylabel(ROT["pdf"])
    ax.legend(h, r, loc="upper right")
    letra(ax, "a")
    ax = axs[1]
    n = np.arange(1, len(ref["conv_pf"]) + 1)
    ax.plot(n, ref["conv_pf"], color=PALETA[0], lw=1.4)
    ax.set_xscale("log")
    ax.set_xlim(1, 10 ** max(2, math.ceil(math.log10(max(len(n), 2)))))
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set_ylim(0, max(0.05, 1.15 * max(ref["conv_pf"], default=0.)))
    ax.set_xlabel("Número de simulações")
    ax.set_ylabel(ROT["pf"])
    ax.text(0.97, 0.95, r"$P_f$ = %s%% ($n_f = %d$ de $S = %d$)" % (_fmt(100 * e["P_f"], 2), e["n_f"], e["S"]),
            transform=ax.transAxes, ha="right", va="top", fontsize=7, color=TINTA)
    letra(ax, "b")
    return salva(fig, pasta, "fig25_ref"), None


def figura_spearman(est, pasta):
    ref = ref_de(est)
    if not ref:
        return None, "sem a referencia"
    paineis = []
    i0 = instante(ref, 0.0)
    if i0 and abs(ref["t"]) > 1e-6 and i0.get("spearman_banda") and i0.get("spearman_campo"):
        paineis.append((0.0, i0["spearman_banda"], i0["spearman_campo"]))
    if ref.get("spearman_banda") and ref.get("spearman_campo"):
        paineis.append((ref["t"], ref["spearman_banda"], ref["spearman_campo"]))
    if not paineis:
        n = (ref.get("estat") or {}).get("S", 0)
        return None, "Spearman sem dados (S = %d; precisa de pelo menos 10 realizacoes)" % n
    nomes = [("k", r"$k_s$"), ("phi", r"$\varphi'$"), ("c", r"$c'$")]
    fig, axs = plt.subplots(1, len(paineis), figsize=(5.4 * len(paineis), 2.9), constrained_layout=True, squeeze=False)
    for j, (ax, (t, sb, sc)) in enumerate(zip(axs[0], paineis)):
        y0 = np.arange(len(nomes))[::-1]
        larg = 0.36
        for i, (rot, s, cor) in enumerate((("banda e massa", sb, PALETA[0]), ("médias do campo", sc, PALETA[1]))):
            ys = y0 + (0.5 - i) * larg
            v = [s[k]["rho_parcial"] for k, _ in nomes]
            ax.barh(ys, v, height=larg * 0.9, color=cor, label=rot)
            lo = [max(0., v[q] - s[k]["ic95"][0]) if s[k]["ic95"][0] is not None else 0. for q, (k, _) in
                  enumerate(nomes)]
            hi = [max(0., s[k]["ic95"][1] - v[q]) if s[k]["ic95"][1] is not None else 0. for q, (k, _) in
                  enumerate(nomes)]
            ax.errorbar(v, ys, xerr=[lo, hi], fmt="none", ecolor=TINTA, elinewidth=0.8, capsize=2)
        ax.set_yticks(y0)
        ax.set_yticklabels([n for _, n in nomes])
        ax.set_xlim(-1, 1)
        ax.set_xticks(np.arange(-1.0, 1.01, 0.25))
        ax.axvline(0, color=TINTA2, lw=0.8)
        ax.set_xlabel(r"Spearman parcial de 2ª ordem, $\rho_{X\,FS,YZ}$")
        ax.set_title("t = %s h (S = %d)" % (_fmt(t), sb["n"]), fontsize=8, color=TINTA)
        ax.legend(loc="lower left")
        if len(paineis) > 1:
            letra(ax, "ab"[j])
    return salva(fig, pasta, "fig26_spearman"), None


def familias_com_dados(est, fams):
    return [f for f in fams if tem_dados(est, f)]


def figura_cov_pdf(est, pasta):
    fams = familias_com_dados(est, ("covk", "covc", "covphi"))
    if not fams:
        return None, "sem casos de covk, covc ou covphi"
    fig, axs = plt.subplots(1, len(fams), figsize=(3.5 * len(fams), 3.3), constrained_layout=True, squeeze=False)
    for j, (ax, fam) in enumerate(zip(axs[0], fams)):
        painel_pdf(ax, fam, casos_familia(est, fam), est["valor_ref"][fam])
        letra(ax, "abc"[j])
    return salva(fig, pasta, "fig27_cov_pdf"), None


def figura_cov_escalar(est, pasta, chave, nome):
    fams = familias_com_dados(est, ("covk", "covc", "covphi"))
    if not fams:
        return None, "sem casos de covk, covc ou covphi"
    fig, axs = plt.subplots(1, len(fams), figsize=(3.5 * len(fams), 3.1), constrained_layout=True, squeeze=False)
    for j, (ax, fam) in enumerate(zip(axs[0], fams)):
        painel_escalar(ax, fam, chave, casos_familia(est, fam), est["valor_ref"][fam], est)
        letra(ax, "abc"[j])
    return salva(fig, pasta, nome), None


def figura_familia(est, pasta, fam, nome):
    """(a) densidades, (b) P_f, (c) CoV, (d) media (+ deterministico; em I, a sigmoide do Rahardjo)"""
    if not tem_dados(est, fam):
        return None, "sem casos da familia %s" % fam
    nossos = casos_familia(est, fam)
    vref = est["valor_ref"][fam]
    fig, axs = plt.subplots(2, 2, figsize=(7.2, 5.6), constrained_layout=True)
    painel_pdf(axs[0, 0], fam, nossos, vref)
    for ax, ch in ((axs[0, 1], "pf"), (axs[1, 0], "cov"), (axs[1, 1], "mu")):
        painel_escalar(ax, fam, ch, nossos, vref, est)
    for ax, t in zip(axs.flat, "abcd"):
        letra(ax, t)
    return salva(fig, pasta, nome), None


def figura_duracao(est, pasta):
    ref = ref_de(est)
    ins = [i for i in (ref or {}).get("instantes", []) if i.get("estat")]
    if len(ins) < 2:
        return None, "a referencia tem um instante so (rodar com tempos=0,6,12,18,24)"
    t = np.array([i["t"] for i in ins])
    fig, axs = plt.subplots(1, 2, figsize=(8.6, 3.3), constrained_layout=True)
    ax = axs[0]
    # trajetorias das realizacoes (as que tem todos os instantes)
    por_real = {}
    for k, i in enumerate(ins):
        for real, fs in zip(i["reais"], i["fs"]):
            por_real.setdefault(real, [math.nan] * len(ins))[k] = fs
    reais = [r for r in sorted(por_real) if all(math.isfinite(x) for x in por_real[r])][:MAX_TRAJETORIAS]
    alfa_tr = min(0.35, 12. / max(len(reais), 1))
    for rr in reais:
        ax.plot(t, por_real[rr], color=MUDO, lw=0.5, alpha=alfa_tr, zorder=1)
    q05 = [i["estat"]["q05"] for i in ins]
    q95 = [i["estat"]["q95"] for i in ins]
    ax.fill_between(t, q05, q95, color=PALETA[0], alpha=0.15, lw=0, zorder=2)
    for q in (q05, q95):
        ax.plot(t, q, color=PALETA[0], lw=0.8, zorder=3)
    mu = [i["estat"]["media"] for i in ins]
    ax.plot(t, mu, color=PALETA[0], lw=1.8, marker="o", ms=4.5, zorder=4)
    h = [Line2D([], [], color=MUDO, lw=0.6), Patch(facecolor=PALETA[0], alpha=0.25, edgecolor=PALETA[0]),
         Line2D([], [], color=PALETA[0], marker="o")]
    r = ["realizações (%d)" % len(reais), "quantis de 5% e 95%", r"$\mu(FS)$ do ajuste log-normal"]
    dt = [(i["t"], i["det_FS"]) for i in ins if i.get("det_FS")]
    if dt:
        ax.plot([a for a, _ in dt], [b for _, b in dt], color=PALETA[1], lw=1.3, ls=(0, (4, 2)), marker="s", ms=4,
                zorder=5)
        h.append(Line2D([], [], color=PALETA[1], ls=(0, (4, 2)), marker="s", ms=4))
        r.append("FS determinístico")
    ax.axhline(1.0, color=TINTA2, lw=0.8)
    ax.set_xticks(t)
    ax.set_xlabel("Duração da chuva (h)")
    ax.set_ylabel(ROT["x_fs"])
    ax.legend(h, r, loc="upper right")
    letra(ax, "a")
    ax = axs[1]
    pf = np.array([i["estat"]["P_f"] for i in ins])
    lo = [p - i["estat"]["IC95_Pf"][0] for p, i in zip(pf, ins)]
    hi = [i["estat"]["IC95_Pf"][1] - p for p, i in zip(pf, ins)]
    pln = [i["estat"]["P_f_lognormal"] for i in ins]
    ax.errorbar(t, pf, yerr=[lo, hi], color=PALETA[0], lw=1.6, marker="o", ms=4.5, capsize=2.5, elinewidth=0.9)
    ax.plot(t, pln, color=PALETA[0], lw=1.0, ls=(0, (4, 2)), marker="o", ms=3.5, mfc="white", mew=0.9)
    ax.set_xticks(t)
    ax.set_ylim(0, max(0.01, 1.12 * max(list(pf + np.array(hi)) + pln)))
    ax.set_xlabel("Duração da chuva (h)")
    ax.set_ylabel(ROT["pf"])
    ax.legend([Line2D([], [], color=PALETA[0], marker="o"),
               Line2D([], [], color=PALETA[0], ls=(0, (4, 2)), marker="o", mfc="white", ms=3.5)],
              ["contagem (IC 95% de Wilson)", "ajuste log-normal"], loc="upper left")
    letra(ax, "b")
    return salva(fig, pasta, "fig_duracao"), None


def figura_fs0_fs24(est, pasta):
    ref = ref_de(est)
    p = (ref or {}).get("pares")
    if not p:
        return None, "a referencia nao tem FS em t = 0 e no instante final"
    a0, a1 = np.asarray(p["fs0"]), np.asarray(p["fs1"])
    cats = p["categoria_t1"]
    fig, ax = plt.subplots(figsize=(4.6, 4.2), constrained_layout=True)
    todos = list(a0) + list(a1) + (p["det"] or []) + [1.0]
    lo, hi = min(todos), max(todos)
    pad = 0.06 * (hi - lo)
    lo, hi = lo - pad, hi + pad
    ax.plot([lo, hi], [lo, hi], color=MUDO, lw=0.9, ls=(0, (1, 1.5)))
    ax.axhline(1.0, color=TINTA2, lw=0.8)
    h, r = [Line2D([], [], color=MUDO, ls=(0, (1, 1.5))), Line2D([], [], color=TINTA2, lw=0.8)], \
        ["1:1 (sem mudança)", "FS = 1 no fim da chuva"]
    for c in CATEGORIAS:
        m = np.array([x == c for x in cats])
        if m.any():
            ax.scatter(a0[m], a1[m], s=16, color=CORES_CATEGORIA[c], edgecolor="white", linewidth=0.4, zorder=3)
            h.append(Line2D([], [], ls="none", marker="o", ms=5, color=CORES_CATEGORIA[c]))
            r.append("nível d'água %s (%d)" % (TEXTO_CATEGORIA[c], int(m.sum())))
    if p.get("det"):
        ax.plot([p["det"][0]], [p["det"][1]], ls="none", marker="*", ms=12, color=PALETA[1], mec="white", mew=0.6,
                zorder=4)
        h.append(Line2D([], [], ls="none", marker="*", ms=10, color=PALETA[1]))
        r.append("determinístico")
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_aspect("equal")
    ax.set_xlabel(r"$FS$ antes da chuva (t = 0)")
    ax.set_ylabel(r"$FS$ no fim da chuva (t = %s h)" % _fmt(p["t1"]))
    if p.get("pearson") is not None:
        ax.text(0.97, 0.45, "Pearson %s; Spearman %s\nFS(%s)/FS(0) médio %s" % (
            _fmt(round(p["pearson"], 3)), _fmt(round(p["spearman"], 3)), _fmt(p["t1"]),
            _fmt(round(p["razao_media"], 3))), transform=ax.transAxes, ha="right", va="center", fontsize=6.8,
            color=TINTA)
    ax.legend(h, r, loc="upper left", fontsize=6.3)
    return salva(fig, pasta, "fig_fs0_fs24"), None


def _grade_mapa(celulas, idx, S, nx=70, ny=30):
    """celulas: [x, y, n, n_falha] (centro da celula de 1 m); idx 2: n, idx 3: n_falha; dividido por S"""
    Z = np.zeros((ny, nx))
    for c in celulas:
        ix, iy = int(math.floor(c[0])), int(math.floor(c[1]))
        if 0 <= ix < nx and 0 <= iy < ny:
            Z[iy, ix] += c[idx]
    Z = Z / S if S else Z * np.nan
    return np.ma.masked_where(~(Z > 0), Z)


def figura_mapa(est, pasta):
    ref = ref_de(est)
    mapa = (ref or {}).get("mapa")
    if not mapa or not any(i["celulas"] for i in mapa["instantes"]):
        return None, "sem o mapa da banda na referencia"
    t1 = ref["t"]
    paineis = []
    for i in mapa["instantes"]:
        if (abs(i["t"]) < 1e-6 or abs(i["t"] - t1) < 1e-6) and i["celulas"] and i["S"]:
            paineis.append((i, 2, i["S"], "t = %s h: fração das %d realizações" % (_fmt(i["t"]), i["S"])))
    i1 = next((i for i in mapa["instantes"] if abs(i["t"] - t1) < 1e-6), None)
    if i1 and i1["n_f"]:
        paineis.append((i1, 3, i1["n_f"], "t = %s h: fração das %d rupturas (FS < 1)" % (_fmt(t1), i1["n_f"])))
    if not paineis:
        return None, "o mapa da referencia nao tem realizacoes em 0 h nem no fim da chuva"
    cmap = LinearSegmentedColormap.from_list("azul", RAMPA)
    # recorte: a regiao com celulas, com folga, dentro do dominio
    xs = [c[0] for i, _, _, _ in paineis for c in i["celulas"]]
    ys = [c[1] for i, _, _, _ in paineis for c in i["celulas"]]
    x0, x1 = max(0, min(xs) - 8), min(70, max(xs) + 8)
    y0, y1 = max(0, min(ys) - 5), 31
    larg = 7.0
    alt_painel = larg * (y1 - y0) / (x1 - x0) * 0.82 + 0.45
    fig, axs = plt.subplots(len(paineis), 1, figsize=(larg, alt_painel * len(paineis)), constrained_layout=True,
                            squeeze=False)
    det = ref.get("det") or {}
    for j, (ax, (i, idx, S, titulo)) in enumerate(zip(axs[:, 0], paineis)):
        Z = _grade_mapa(i["celulas"], idx, S)
        contorno_talude(ax, est["talude"])
        pc = ax.pcolormesh(np.arange(71), np.arange(31), Z, cmap=cmap, vmin=0, vmax=1, zorder=1)
        banda = next((b for b in det.get("bandas", []) if abs(b["t"] - i["t"]) < 1e-6), None)
        if banda and banda["pontos"]:
            p = np.asarray(banda["pontos"])
            ax.scatter(p[:, 0], p[:, 1], s=0.8, color=PALETA[1], alpha=0.6, lw=0, zorder=4)
        ax.set_aspect("equal")
        ax.set_xlim(x0, x1)
        ax.set_ylim(y0, y1)
        ax.grid(False)
        ax.set_title(titulo, fontsize=8, color=TINTA, loc="left")
        ax.set_ylabel("y (m)")
        if j == len(paineis) - 1:
            ax.set_xlabel("x (m)")
        cb = fig.colorbar(pc, ax=ax, shrink=0.9, pad=0.01)
        cb.ax.tick_params(labelsize=6.5)
        cb.set_label("fração", fontsize=7)
        if banda:
            ax.legend([Line2D([], [], ls="none", marker="o", ms=3, color=PALETA[1])],
                      ["banda do determinístico"], loc="lower left", fontsize=6.3)
        letra(ax, "abc"[j])
    return salva(fig, pasta, "fig_mapa"), None


def _barra_empilhada(ax, pos, fracoes, horizontal, largura=0.7):
    base = 0.
    for c in CATEGORIAS:
        f = fracoes.get(c) or 0.
        if f <= 0:
            continue
        if horizontal:
            ax.barh(pos, f, left=base, height=largura, color=CORES_CATEGORIA[c], edgecolor="white", linewidth=0.8)
        else:
            ax.bar(pos, f, bottom=base, width=largura, color=CORES_CATEGORIA[c], edgecolor="white", linewidth=0.8)
        if f >= 0.08:
            cor = "white" if c in ("I", "II", "III") else TINTA
            xy = (base + f / 2, pos) if horizontal else (pos, base + f / 2)
            ax.text(*xy, "%d%%" % round(100 * f), ha="center", va="center", fontsize=6.3, color=cor)
        base += f


def figura_nivel(est, pasta):
    """(a) referencia: categorias com a duracao da chuva; (b) cada caso no instante das familias"""
    ref = ref_de(est)
    casos = [c for c in est["casos"] if c.get("nivel") and sum(c["nivel"]["contagem"].values())]
    if not casos:
        return None, "sem os niveis d'agua"
    ins = [i for i in (ref or {}).get("instantes", []) if sum(i["nivel"]["contagem"].values())]
    tem_a = len(ins) > 1
    # posicoes de (b): de cima para baixo, com meio espaco entre familias
    pos, y, fam_ant = [], 0., None
    for c in casos:
        if fam_ant is not None and c["familia"] != fam_ant:
            y += 0.5
        pos.append(y)
        y += 1.
        fam_ant = c["familia"]
    alt_b = 0.2 * (pos[-1] + 1) + 0.9
    alt_a = 2.3 if tem_a else 0.
    fig = plt.figure(figsize=(7.2, alt_a + alt_b + 0.7), constrained_layout=True)
    if tem_a:
        sub = fig.subfigures(2, 1, height_ratios=[alt_a, alt_b + 0.7])
        ax_a, ax_b = sub[0].subplots(), sub[1].subplots()
        sf_leg = sub[1]
    else:
        ax_b, sf_leg = fig.subplots(), fig
    if tem_a:
        ax = ax_a
        for k, i in enumerate(ins):
            _barra_empilhada(ax, k, i["nivel"]["fracao"], horizontal=False, largura=0.6)
            dc = i["nivel"].get("det_categoria")
            if dc:
                ax.text(k, 1.02, "det. %s" % dc, ha="center", va="bottom", fontsize=6.3, color=TINTA2)
        ax.set_xticks(range(len(ins)))
        ax.set_xticklabels([_fmt(i["t"]) for i in ins])
        ax.set_xlabel("Duração da chuva (h)")
        ax.set_ylabel("Fração das realizações")
        ax.set_ylim(0, 1.12)
        ax.set_yticks(np.arange(0, 1.01, 0.25))
        ax.grid(axis="x", visible=False)
        ax.set_title("referência", fontsize=8, color=TINTA)
        letra(ax, "a")
    ax = ax_b
    for c, y in zip(casos, pos):
        _barra_empilhada(ax, y, c["nivel"]["fracao"], horizontal=True, largura=0.72)
        dc = c["nivel"].get("det_categoria")
        ax.text(1.02, y, "det. %s" % dc if dc else "", ha="left", va="center", fontsize=6.3, color=TINTA2)
    ax.set_yticks(pos)
    ax.set_yticklabels(["referência" if c["familia"] == "ref" else rotulo_valor(c["familia"], c["valor"])
                        for c in casos], fontsize=7)
    ax.set_ylim(pos[-1] + 0.6, -0.6)
    ax.set_xlim(0, 1.12)
    ax.set_xticks(np.arange(0, 1.01, 0.25))
    ax.set_xlabel("Fração das realizações (t = %s h)" % _fmt(est["t_familias"]))
    ax.grid(axis="y", visible=False)
    ax.set_title("casos", fontsize=8, color=TINTA)
    if tem_a:
        letra(ax, "b")
    hl = [Patch(color=CORES_CATEGORIA[c]) for c in CATEGORIAS]
    sf_leg.legend(hl, [TEXTO_CATEGORIA[c] for c in CATEGORIAS], loc="outside lower center", ncol=4, fontsize=6.6,
                  title="Categoria do nível d'água (Tabela 3 de Rahardjo et al. 2007)", title_fontsize=6.8)
    return salva(fig, pasta, "fig_nivel"), None


def todas(est, pasta):
    """todas as figuras com dados; devolve ({nome: pdf}, {nome: motivo de ter pulado})"""
    os.makedirs(pasta, exist_ok=True)
    tarefas = [
        ("fig24_det", lambda: figura_det(est, pasta)),
        ("fig25_ref", lambda: figura_ref(est, pasta)),
        ("fig26_spearman", lambda: figura_spearman(est, pasta)),
        ("fig27_cov_pdf", lambda: figura_cov_pdf(est, pasta)),
        ("fig28_cov_pf", lambda: figura_cov_escalar(est, pasta, "pf", "fig28_cov_pf")),
        ("fig29_cov_mu", lambda: figura_cov_escalar(est, pasta, "mu", "fig29_cov_mu")),
        ("fig30_cov_cov", lambda: figura_cov_escalar(est, pasta, "cov", "fig30_cov_cov")),
        ("fig31_s", lambda: figura_familia(est, pasta, "s", "fig31_s")),
        ("fig32_alfa", lambda: figura_familia(est, pasta, "alfa", "fig32_alfa")),
        ("fig33_I", lambda: figura_familia(est, pasta, "I", "fig33_I")),
        ("fig_duracao", lambda: figura_duracao(est, pasta)),
        ("fig_fs0_fs24", lambda: figura_fs0_fs24(est, pasta)),
        ("fig_mapa", lambda: figura_mapa(est, pasta)),
        ("fig_nivel", lambda: figura_nivel(est, pasta)),
    ]
    feitas, puladas = {}, {}
    for nome, f in tarefas:
        arq, motivo = f()
        if arq:
            feitas[nome] = arq
        else:
            puladas[nome] = motivo
    return feitas, puladas


def main(argv):
    if len(argv) != 3:
        sys.exit(__doc__)
    with open(argv[1]) as f:
        est = json.load(f)
    feitas, puladas = todas(est, argv[2])
    for n, a in feitas.items():
        print("figura:", a)
    for n, m in puladas.items():
        print("pulada: %s (%s)" % (n, m))


if __name__ == "__main__":
    main(sys.argv)
