#!/usr/bin/env python3
"""figuras_ls_ceron.py: as figuras da Secao 6 de Vargas Ceron et al. (2025) para a campanha ls_ceron_casos.sbatch
(agua do artigo; FORM modal + line sampling + Monte Carlo direto pareado), nos eixos e com o estilo do figuras_sat.py
(o Ceron tracejado, este trabalho cheio; eixos ampliados so quando os dados nao cabem, com nota):
  - fig25: (a) histograma e log-normal de Gamma da referencia (MC) contra a log-normal da tabela do artigo; (b) a
    convergencia de P_f: o MC pelo numero de simulacoes e o line sampling pelo numero de analises do MEF;
  - fig27, fig29, fig30: densidades log-normais, media e CoV de Gamma contra os CoV de k_v, c e phi (do MC);
  - fig28: P_f contra os CoV: line sampling (IC 95%), FORM (Phi(-beta)) e MC (IC 95%), contra a Tabela 3;
  - fig31 a fig33 (s, alfa, h_w/H): (a) densidades, (b) P_f (LS, FORM, MC), (c) CoV e (d) media de Gamma;
  - fig_det: o Gamma deterministico (as medias) contra alfa e h_w/H (Tabelas 5 e 6).
As densidades, a media e o CoV sao do Monte Carlo (sem MC, esses paineis ficam so com o Ceron).

    python3 -I figuras_ls_ceron.py <resumo.json do resumo_ls_ceron.py> <pasta de saida>"""

import json
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
sys.path.insert(0, AQUI)
import ceron_dig  # noqa: E402
import figuras_sat as F  # noqa: E402

COR_LS, COR_FORM, COR_MC = F.PALETA[0], F.PALETA[1], F.PALETA[2]


def estat_mc(m):
    """o MC no formato do figuras_sat (S, P_f, media, CoV, mu_ln, sigma_ln), ou None"""
    if not m or m.get("S", 0) < 10 or m.get("mu_ln") is None:
        return None
    return {"S": m["S"], "P_f": m["P_f"], "media": m["media"], "CoV": m["CoV"], "mu_ln": m["mu_ln"],
            "sigma_ln": m["sigma_ln"]}


def monta_est(dados):
    casos = []
    for t in dados:
        valor = F.VALOR_REF.get(t["familia"]) if t["familia"] == "ref" else t["valor"]
        casos.append({"id": t["id"], "familia": t["familia"], "valor": valor, "estat": estat_mc(t["mc"]), "t": t})
    return {"casos": casos}


def nossos(est, fam):
    """os casos da familia, com a referencia no seu valor"""
    return F.casos_familia(est, fam)


def painel_pf(ax, fam, fig_id, nos, ceron_tab, amp):
    """P_f contra o parametro: Ceron (Tabelas 3 a 6), line sampling (IC 95%), FORM e MC (IC 95%)"""
    lim = F.FIG[fig_id]
    x0, x1 = lim["x"][0], lim["x"][1]

    def fora(x):
        return x < x0 - 1e-12 or x > x1 + 1e-12
    extra, todos_y, todos_x = [], [], []
    cx = [r["x"] for r in ceron_tab if not fora(r["x"])]
    cy = [r["pf"] for r in ceron_tab if not fora(r["x"])]
    o = np.argsort(cx)
    ax.plot(np.array(cx)[o], np.array(cy)[o], color=F.TINTA2, ls=F.CERON, lw=1.2, marker="s", ms=4.5, mfc="white",
            mec=F.TINTA2, mew=1.0)
    todos_y += cy
    todos_x += cx
    for r in ceron_tab:
        if fora(r["x"]):
            ax.axhline(r["pf"], color=F.TINTA2, ls=(0, (1, 1.5)), lw=1.0)
            extra.append((Line2D([], [], color=F.TINTA2, ls=(0, (1, 1.5))), "Ceron, %s = %s" % (F.NOME[fam], F._fmt(r["x"]))))
            todos_y.append(r["pf"])
    series = []
    for v, c in sorted(nos.items()):
        t = c["t"]
        ls, fo, mc = t["ls"], t["form"], t["mc"]
        series.append((v, ls.get("P_f"), ls.get("IC95"), fo.get("P_f") if fo.get("estado") == "pronto" else None,
                       mc.get("P_f") if mc.get("S", 0) >= 10 else None, mc.get("S", 0)))
    # line sampling
    lx, ly, le = [], [], [[], []]
    fx, fy, mx, my, me = [], [], [], [], [[], []]
    for v, p, ic, pfo, pmc, nmc in series:
        if fora(v):
            if p is not None:
                ax.axhline(p, color=COR_LS, ls=(0, (1, 1.5)), lw=1.2)
                extra.append((Line2D([], [], color=COR_LS, ls=(0, (1, 1.5))),
                              "line sampling, %s = %s" % (F.NOME[fam], F._fmt(v))))
                todos_y.append(p)
            continue
        todos_x.append(v)
        if p is not None:
            lx.append(v)
            ly.append(p)
            lo, hi = (ic if ic else (p, p))
            le[0].append(max(p - lo, 0.))
            le[1].append(max(hi - p, 0.))
            todos_y += [lo, hi]
        if pfo is not None:
            fx.append(v)
            fy.append(pfo)
            todos_y.append(pfo)
        if pmc is not None:
            mx.append(v)
            my.append(pmc)
            d = 1.96 * math.sqrt(max(pmc * (1 - pmc), 1e-12) / nmc)
            me[0].append(min(d, pmc))
            me[1].append(d)
            todos_y += [pmc - d, pmc + d]
    if lx:
        ax.errorbar(lx, ly, yerr=le, color=COR_LS, ls=F.NOSSO, lw=1.6, marker="o", ms=4.5, capsize=2.5, elinewidth=0.9,
                    zorder=4)
    if fx:
        ax.plot(fx, fy, color=COR_FORM, ls=(0, (1, 1.2)), lw=1.1, marker="D", ms=3.8, mfc="white", mec=COR_FORM, zorder=3)
    if mx:
        ax.errorbar(mx, my, yerr=me, color=COR_MC, ls="none", marker="^", ms=3.8, capsize=2, elinewidth=0.8, zorder=3)
    amp.eixo(ax, fig_id, "x", lim["x"], (min(todos_x, default=None), max(todos_x, default=None)))
    amp.eixo(ax, fig_id, "y", lim["y"], (min(todos_y, default=None), max(todos_y, default=None)))
    ax.set_xlabel(F.NOME[fam])
    ax.set_ylabel(F.ROT["pf"])
    h = [Line2D([], [], color=F.TINTA2, ls=F.CERON, marker="s", mfc="white", mec=F.TINTA2),
         Line2D([], [], color=COR_LS, ls=F.NOSSO, marker="o"),
         Line2D([], [], color=COR_FORM, ls=(0, (1, 1.2)), marker="D", mfc="white", mec=COR_FORM)]
    r = ["Ceron et al. (2025)", "line sampling (IC 95%)", r"FORM modal, $\Phi(-\beta)$"]
    if mx:
        h.append(Line2D([], [], color=COR_MC, ls="none", marker="^"))
        r.append("Monte Carlo (IC 95%%), N = %d" % max(s[5] for s in series))
    for hh, rr in extra:
        h.append(hh)
        r.append(rr)
    ax.legend(h, r, loc="best", fontsize=6.3)


def tem_mc(est, fam):
    return any(c["estat"] for c in est["casos"] if c["familia"] in (fam, "ref"))


def painel_escalar_mc(ax, fam, fig_id, chave, nos, ceron_tab, amp):
    """media ou CoV de Gamma: Ceron e o MC deste trabalho (o figuras_sat)"""
    F.painel_escalar(ax, fam, fig_id, chave, nos, ceron_tab, amp)
    leg = ax.get_legend()
    if leg is not None:
        textos = leg.get_texts()
        if len(textos) > 1:
            textos[1].set_text("Monte Carlo (este trabalho)")


def figura28(est, cer, saida, amp):
    fig, axs = plt.subplots(1, 3, figsize=(10.5, 3.3), constrained_layout=True)
    for ax, fam, t in zip(axs, ("covk", "covc", "covphi"), "abc"):
        painel_pf(ax, fam, F.PAINEIS[fam][1], nossos(est, fam), cer["tab"][fam], amp)
        F.letra(ax, t)
    arq = os.path.join(saida, "fig28.pdf")
    fig.savefig(arq)
    plt.close(fig)
    return arq


def figura_cov_mc(n, chave, est, cer, saida, amp):
    idx = {"mu": 3, "cov": 2}[chave]
    fig, axs = plt.subplots(1, 3, figsize=(10.5, 3.1), constrained_layout=True)
    for ax, fam, t in zip(axs, ("covk", "covc", "covphi"), "abc"):
        painel_escalar_mc(ax, fam, F.PAINEIS[fam][idx], chave, nossos(est, fam), cer["tab"][fam], amp)
        F.letra(ax, t)
    arq = os.path.join(saida, "fig%d.pdf" % n)
    fig.savefig(arq)
    plt.close(fig)
    return arq


def figura_familia(fam, est, cer, saida, amp):
    f_pdf, f_pf, f_cov, f_mu = F.PAINEIS[fam]
    nos = nossos(est, fam)
    fig, axs = plt.subplots(2, 2, figsize=(7.4, 5.8), constrained_layout=True)
    F.painel_pdf(axs[0, 0], fam, f_pdf, nos, cer["pdfs"].get(fam, {}), amp)
    painel_pf(axs[0, 1], fam, f_pf, nos, cer["tab"][fam], amp)
    painel_escalar_mc(axs[1, 0], fam, f_cov, "cov", nos, cer["tab"][fam], amp)
    painel_escalar_mc(axs[1, 1], fam, f_mu, "mu", nos, cer["tab"][fam], amp)
    for ax, t in zip(axs.flat, "abcd"):
        F.letra(ax, t)
    arq = os.path.join(saida, "fig%s.pdf" % f_pdf[:2])
    fig.savefig(arq)
    plt.close(fig)
    return arq


def figura27(est, cer, saida, amp):
    fig, axs = plt.subplots(1, 3, figsize=(10.5, 3.3), constrained_layout=True)
    for ax, fam, t in zip(axs, ("covk", "covc", "covphi"), "abc"):
        F.painel_pdf(ax, fam, F.PAINEIS[fam][0], nossos(est, fam), cer["pdfs"].get(fam, {}), amp)
        F.letra(ax, t)
    arq = os.path.join(saida, "fig27.pdf")
    fig.savefig(arq)
    plt.close(fig)
    return arq


def figura25(est, cer, saida, amp):
    ref = next((c for c in est["casos"] if c["familia"] == "ref"), None)
    if ref is None:
        return None
    t = ref["t"]
    fig, axs = plt.subplots(1, 2, figsize=(8.8, 3.4), constrained_layout=True)
    ax = axs[0]
    tab = cer["tab"]["alfa"][0]  # a referencia (alfa = 1): mu e sigma da Tabela 3/5
    xc, yc = ceron_dig.lognormal_da_tabela(tab["mu"], tab["sigma"])
    ax.plot(xc, yc, color=F.TINTA2, ls=F.CERON, lw=1.4)
    g = np.asarray(t["mc"].get("gamma") or [])
    e = ref["estat"]
    ys = [max(yc)]
    if e and len(g):
        ax.hist(g, bins="auto", density=True, color=COR_LS, alpha=0.28, edgecolor="white", linewidth=0.4)
        x = np.linspace(max(1e-3, g.min() * 0.7), g.max() * 1.2, 500)
        y = F.lognormal_pdf(x, e["mu_ln"], e["sigma_ln"])
        ax.plot(x, y, color=COR_LS, lw=1.6)
        ys += [y.max(), np.histogram(g, bins="auto", density=True)[0].max()]
    amp.eixo(ax, "25a", "x", F.FIG["25a"]["x"], (min(xc), max(xc)) if not len(g) else (g.min(), g.max()))
    amp.eixo(ax, "25a", "y", F.FIG["25a"]["y"], (0.0, max(ys)))
    ax.axvline(1.0, color=F.TINTA2, lw=0.8)
    ax.set_xlabel(F.ROT["x_gamma"])
    ax.set_ylabel(F.ROT["pdf"])
    ax.legend([Line2D([], [], color=F.TINTA2, ls=F.CERON), Patch(color=COR_LS, alpha=0.28), Line2D([], [], color=COR_LS)],
              ["Ceron: log-normal (Tabela 3)", "Monte Carlo: histograma (N = %d)" % len(g), "Monte Carlo: log-normal"],
              loc="upper right")
    F.letra(ax, "a")
    ax = axs[1]
    yy = [t["pf_artigo"]]
    if len(g):
        falha = (g < 1.).astype(float)
        n = np.arange(1, len(g) + 1)
        ax.plot(n, np.cumsum(falha) / n, color=COR_MC, lw=1.3, label="Monte Carlo (por simulação)")
        yy.append(float(np.max(np.cumsum(falha) / n)))
    p = t["ls"].get("p") or []
    if p:
        an_linha = (t["ls"].get("analises") or 0) / max(len(p), 1)
        n = np.arange(1, len(p) + 1)
        base = (t["form"].get("analises") or 0) + 2  # o FORM e as 2 analises da inclinacao
        ax.plot(base + n * an_linha, np.cumsum(p) / n, color=COR_LS, lw=1.6, label="line sampling (por análise)")
        yy.append(float(np.max(np.cumsum(p) / n)))
    ax.axhline(t["pf_artigo"], color=F.TINTA2, ls=F.CERON, lw=1.2, label="Ceron: $P_f$ = %s%%" % F._fmt(100 * t["pf_artigo"]))
    if t["form"].get("estado") == "pronto":
        ax.axhline(t["form"]["P_f"], color=COR_FORM, ls=(0, (1, 1.2)), lw=1.1, label=r"FORM modal, $\Phi(-\beta)$")
    ax.set_xscale("log")
    ax.set_xlim(*F.FIG["25b"]["x"][:2])
    amp.eixo(ax, "25b", "y", F.FIG["25b"]["y"], (0.0, max(yy)))
    ax.set_xlabel("Número de simulações (MC) ou de análises do MEF (LS)")
    ax.set_ylabel(F.ROT["pf"])
    ax.legend(loc="upper right")
    F.letra(ax, "b")
    arq = os.path.join(saida, "fig25.pdf")
    fig.savefig(arq)
    plt.close(fig)
    return arq


def figura_det(est, cer, saida):
    fig, axs = plt.subplots(1, 2, figsize=(7.6, 2.9), constrained_layout=True)
    for ax, fam, t in zip(axs, ("alfa", "hw"), "ab"):
        tab = cer["tab"][fam]
        ax.plot([r["x"] for r in tab], [r["det"] for r in tab], color=F.TINTA2, ls=F.CERON, marker="s", ms=4.5,
                mfc="white", mec=F.TINTA2, lw=1.2, label="Ceron et al. (2025)")
        nos = nossos(est, fam)
        xs = sorted(v for v, c in nos.items() if (c["t"].get("det") or {}).get("Gamma"))
        ax.plot(xs, [nos[v]["t"]["det"]["Gamma"] for v in xs], color=COR_LS, ls=F.NOSSO, marker="o", ms=4.5, lw=1.6,
                label="este trabalho")
        ax.set_xlabel(F.NOME[fam])
        ax.set_ylabel(r"$\Gamma$ determinístico")
        ax.legend(loc="best")
        F.letra(ax, t)
    arq = os.path.join(saida, "fig_det.pdf")
    fig.savefig(arq)
    plt.close(fig)
    return arq


def todas(dados, saida):
    os.makedirs(saida, exist_ok=True)
    est = monta_est(dados)
    cer = ceron_dig.carrega(None)
    amp = F.Ampliacoes()
    arqs = {"25": figura25(est, cer, saida, amp), "27": figura27(est, cer, saida, amp),
            "28": figura28(est, cer, saida, amp), "29": figura_cov_mc(29, "mu", est, cer, saida, amp),
            "30": figura_cov_mc(30, "cov", est, cer, saida, amp), "det": figura_det(est, cer, saida)}
    for fam in ("s", "alfa", "hw"):
        arqs[F.PAINEIS[fam][0][:2]] = figura_familia(fam, est, cer, saida, amp)
    with open(os.path.join(saida, "notas_eixos.txt"), "w") as f:
        for k, v in sorted(amp.notas.items()):
            f.write("fig%s: %s\n" % (k, "; ".join(v)))
    return arqs, amp.notas


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    dados = json.load(open(sys.argv[1]))
    arqs, notas = todas(dados, sys.argv[2])
    for k, a in sorted(arqs.items(), key=lambda t: t[0]):
        if a:
            print(a)
    for k, v in sorted(notas.items()):
        print("fig%s: %s" % (k, "; ".join(v)))


if __name__ == "__main__":
    main()
