#!/usr/bin/env python3
"""figura_mc.py: figura e resumo de um Monte Carlo do mc.py.

    python3 -I figura_mc.py csv=/home/diogo/projects/runs/ceron_mc_python/piloto200/mc.csv nome=piloto200 \
        [parede=58] [procs=4] [conf=.../conferencia_30.csv]

Gera figs/<nome>.png/.pdf (histograma de Gamma com Gamma = 1 e a normal/lognormal de media 1,353 e desvio 0,318
do Ceron; P_f acumulada x numero de realizacoes com o IC de 95% de Wilson e a linha 11,5%) e dados/<nome>_resumo.txt.
"""
import csv
import math
import os
import sys

import numpy as np

AQUI = os.path.dirname(os.path.abspath(__file__))
PF_CERON, MED_CERON, DP_CERON, COV_CERON = 0.115, 1.353, 0.318, 0.0277
Z = 1.959963984540054


def wilson(nf, S):
    p = nf / S
    den = 1 + Z * Z / S
    cen = (p + Z * Z / (2 * S)) / den
    meia = Z * np.sqrt(p * (1 - p) / S + Z * Z / (4 * S * S)) / den
    return cen - meia, cen + meia


def br(x, nd=2):
    return f"{x:.{nd}f}".replace(".", ",")


def main(argv):
    a = dict(csv="", nome="piloto200", parede="", procs="4", conf="")
    for s in argv:
        k, _, v = s.partition("=")
        a[k] = v
    linhas = sorted(csv.DictReader(open(a["csv"])), key=lambda r: int(r["real"]))
    G = np.array([float(r["Gamma"]) for r in linhas])
    nnan = int(np.sum(~np.isfinite(G)))
    G = G[np.isfinite(G)]
    S = len(G)
    falha = G < 1.0
    nf = int(falha.sum())
    pf = nf / S
    lo, hi = wilson(nf, S)
    cov = math.sqrt((1 - pf) / (S * pf)) if nf else float("inf")
    med, dp = G.mean(), G.std(ddof=1)
    ep_med = dp / math.sqrt(S)
    t = np.array([float(r["tempo"]) for r in linhas])
    tc = np.array([float(r["t_campos"]) for r in linhas])
    td = np.array([float(r["t_darcy"]) for r in linhas])
    tb = np.array([float(r["t_busca"]) for r in linhas])
    nav = np.array([float(r["avaliacoes"]) for r in linhas])
    sem = sorted({r.get("semente", "") for r in linhas})

    # ---------------- figura ----------------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from scipy import stats
    tinta, tinta2, azul, laranja, verde = "#0b0b0b", "#52514e", "#2a78d6", "#eb6834", "#1baf7a"
    plt.rcParams.update({"font.size": 9, "axes.edgecolor": tinta2, "axes.labelcolor": tinta,
                         "xtick.color": tinta2, "ytick.color": tinta2, "axes.spines.top": False,
                         "axes.spines.right": False})
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 3.8))
    bins = np.arange(0.0, max(3.0, G.max()) + 0.1, 0.1)
    ax1.hist(G, bins=bins, density=True, color=azul, alpha=0.55, edgecolor="#fcfcfb", linewidth=1,
             label=f"este trabalho (S = {S})")
    xx = np.linspace(0.0, bins[-1], 400)
    ax1.plot(xx, stats.norm.pdf(xx, MED_CERON, DP_CERON), color=laranja, lw=2,
             label="normal, Ceron (1,353; 0,318)")
    s_ln = math.sqrt(math.log(1 + (DP_CERON / MED_CERON) ** 2))
    m_ln = math.log(MED_CERON) - 0.5 * s_ln ** 2
    ax1.plot(xx, stats.lognorm.pdf(xx, s_ln, scale=math.exp(m_ln)), color=laranja, lw=2, ls="--",
             label="lognormal, Ceron (1,353; 0,318)")
    ax1.axvline(1.0, color=tinta, lw=1.2)
    ax1.text(0.97, ax1.get_ylim()[1] * 0.95, "Γ = 1", color=tinta, va="top", ha="right")
    ax1.set_xlabel("Γ (fator de seguranca)")
    ax1.set_ylabel("densidade")
    ax1.set_title(f"Histograma de Γ: media {br(med, 3)}, desvio {br(dp, 3)}", fontsize=9, color=tinta)
    ax1.legend(frameon=False, fontsize=8)

    n = np.arange(1, S + 1)
    cum = np.cumsum(falha) / n
    lo_n, hi_n = wilson(np.cumsum(falha), n)
    ax2.fill_between(n, 100 * lo_n, 100 * hi_n, color=azul, alpha=0.18, lw=0, label="IC 95% (Wilson)")
    ax2.plot(n, 100 * cum, color=azul, lw=2, label="P_f acumulada")
    ax2.axhline(100 * PF_CERON, color=laranja, lw=1.5, ls="--", label="Ceron: 11,50%")
    ax2.set_xlabel("numero de realizacoes")
    ax2.set_ylabel("P_f = P(Γ < 1) (%)")
    ax2.set_xlim(1, S)
    ax2.set_ylim(0, min(100, max(40, 100 * hi_n[min(19, S - 1):].max() * 1.1)))
    ax2.set_title(f"P_f = {br(100 * pf)}%  [{br(100 * lo)}%; {br(100 * hi)}%]", fontsize=9, color=tinta)
    ax2.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(AQUI, "figs", f"{a['nome']}.{ext}"), dpi=200)

    # ---------------- resumo ----------------
    L = []
    L.append(f"Monte Carlo do Ceron (2025, Secao 6, configuracao de referencia) com o avaliador rapido (rapido.py)")
    L.append(f"arquivo: {a['csv']}")
    L.append(f"semente: {', '.join(sem)} (realizacao i -> numpy.random.default_rng(semente + i)); "
             f"realizacoes {linhas[0]['real']} a {linhas[-1]['real']}; M = {linhas[0].get('M', '?')} termos KL por campo; "
             f"h = {linhas[0].get('h', '?')} m, hD = {linhas[0].get('hD', '?')} m")
    L.append("")
    L.append(f"{'':28s}{'este trabalho':>18s}{'Ceron (S = 10 000)':>22s}")
    L.append(f"{'S':28s}{S:>18d}{'10 000':>22s}")
    L.append(f"{'falhas (Gamma < 1)':28s}{nf:>18d}{'-':>22s}")
    L.append(f"{'P_f':28s}{br(100 * pf) + '%':>18s}{'11,50%':>22s}")
    L.append(f"{'IC 95% (Wilson)':28s}{'[' + br(100 * lo) + '; ' + br(100 * hi) + ']%':>18s}{'-':>22s}")
    L.append(f"{'CoV(P_f) (Eq. 68)':28s}{br(100 * cov) + '%':>18s}{'2,77%':>22s}")
    L.append(f"{'media de Gamma':28s}{br(med, 3) + ' +- ' + br(ep_med, 3):>18s}{'1,353':>22s}")
    L.append(f"{'desvio de Gamma':28s}{br(dp, 3):>18s}{'0,318':>22s}")
    L.append(f"{'realizacoes sem mecanismo':28s}{nnan:>18d}")
    L.append("")
    L.append(f"(media - 1,353)/erro padrao = {br((med - MED_CERON) / ep_med)}")
    L.append(f"classes criticas: { {str(k): int(v) for k, v in zip(*np.unique([r['classe'] for r in linhas], return_counts=True))} }; "
             f"pelo pe (sB = 0): {sum(float(r['sB']) == 0 for r in linhas)}; "
             f"xA de {br(min(float(r['xA']) for r in linhas))} a {br(max(float(r['xA']) for r in linhas))} m; "
             f"profundidade maxima {br(max(float(r['y_max']) for r in linhas), 3)} m")
    soma = lambda k: int(sum(float(r[k]) for r in linhas))  # noqa: E731
    L.append(f"descartes por normalidade {soma('n_normalidade')}, phi_i sem convergencia {soma('n_phi_nao_conv')}, "
             f"marchas na versao sequencial {soma('fallbacks')}, marchas com a regra 'so sobe' "
             f"{soma('n_max_regra')} (em {sum(int(r['n_max_regra']) > 0 for r in linhas)} realizacoes; "
             f"o Gamma final e de uma marcha fria)")
    L.append("")
    L.append(f"tempos por realizacao (um processo; {a['procs']} processos em paralelo, nice 10, 1 thread BLAS):")
    L.append(f"  total {br(t.mean())} s (mediana {br(np.median(t))}, max {br(t.max())}); campos "
             f"{br(1e3 * tc.mean(), 0)} ms, Darcy {br(1e3 * td.mean(), 0)} ms, busca {br(tb.mean())} s "
             f"({br(nav.mean(), 0)} avaliacoes, {br(1e3 * tb.mean() / nav.mean())} ms cada)")
    if a["parede"]:
        tp = float(a["parede"])
        L.append(f"  tempo de parede total: {br(tp, 0)} s ({br(tp / S)} s por realizacao com {a['procs']} processos)")
        L.append(f"  estimativa para 10 000: {br(1e4 * tp / S / 60, 0)} min aqui ({a['procs']} processos); fermi com 48 "
                 f"processos: ~{br(1e4 * t.mean() / 48 / 60, 0)} min (com o mesmo tempo por realizacao) + ~1 min de partida")
    if a["conf"]:
        c = list(csv.DictReader(open(a["conf"])))
        Gp = {int(r["real"]): float(r["Gamma"]) for r in linhas}
        rel = np.array([Gp[int(r["real"])] / float(r["G_exato"]) - 1 for r in c])
        ge = np.array([float(r["G_exato"]) for r in c])
        gp = np.array([Gp[int(r["real"])] for r in c])
        L.append("")
        L.append(f"conferencia rapido x exato (cinematico.minimizar no CampoMalha dos mesmos campos), realizacoes "
                 f"{c[0]['real']} a {c[-1]['real']}:")
        L.append(f"  Gamma rapido/exato - 1: max {rel.max():+.1e}, min {rel.min():+.1e}, media {rel.mean():+.1e}; "
                 f"acima de 0,5%: {int(np.sum(rel > 0.005))}; acima de 0,1%: {int(np.sum(rel > 0.001))}")
        L.append(f"  falhas: rapido {int(np.sum(gp < 1))}, exato {int(np.sum(ge < 1))}; classificadas diferente: "
                 f"{int(np.sum((gp < 1) != (ge < 1)))}; classe critica igual: "
                 f"{sum(r['cl_rapido'] == r['cl_exato'] for r in c)} de {len(c)}")
        L.append(f"  avaliacoes: rapido {br(np.mean([float(r['nav_rapido']) for r in c]), 0)}, exato "
                 f"{br(np.mean([float(r['nav_exato']) for r in c]), 0)}; busca: rapido "
                 f"{br(np.mean([float(r['t_rapido']) for r in c]))} s, exato {br(np.mean([float(r['t_exato']) for r in c]), 1)} s")
    txt = "\n".join(L) + "\n"
    with open(os.path.join(AQUI, "dados", f"{a['nome']}_resumo.txt"), "w") as f:
        f.write(txt)
    print(txt)


if __name__ == "__main__":
    main(sys.argv[1:])
