#!/usr/bin/env python3
"""plot_v1.py <pasta> [--saida <arquivo.png>]: compara os resultados do V1_LiakopoulosTrifasico (perfis.csv e
historico.csv gravados na <pasta> pelo executavel) com as curvas digitalizadas de Laloui, Klubertanz e Vulliet (2003)
em ../dados. Grava a figura (padrao: ../figuras/v1_comparacao.png) e imprime a tabela "artigo x nosso".
Paineis: (a) p_a em 30 min contra a faixa das cinco equipes do benchmark ALERT (Fig. 6a); (b) p_a no tempo com a lei do
ar do benchmark contra os minimos da Fig. 9a; (c) p_w em 5, 10, 30 e 600 min, duas e tres fases, contra o experimento
(Fig. 5a) e o topo da Fig. 8b; (d) vazao na base (Fig. 5b); (e) recalque do topo (Figs. 6b e 8a); (f) p_w no topo no
tempo, duas e tres fases (Fig. 8b).
Requer numpy e matplotlib."""
import csv
import os
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

AQUI = os.path.dirname(os.path.abspath(__file__))
DADOS = os.path.join(AQUI, "..", "dados")
# paleta categorica (ordem fixa): azul, laranja, aqua, amarelo; referencias em tinta neutra
COR = {"bifasico": "#2a78d6", "trifasico": "#eb6834", "trifasico_alert": "#1baf7a", "trifasico_laloui": "#eda100"}
NOME = {"bifasico": "duas fases (ar passivo)", "trifasico": "tres fases, dados impressos (k_ra >= 1%)",
        "trifasico_alert": "tres fases, ar do benchmark (k_ra >= 1e-4)", "trifasico_laloui": "tres fases, n e k_rw de Laloui"}
TINTA, CINZA = "#0b0b0b", "#8a8984"
# trifasico_laloui quase coincide com trifasico: tracejado fino por cima, para os dois aparecerem
TRACO = {"trifasico_laloui": dict(ls="--", lw=1.2)}


def le_csv(arq):
    with open(arq) as f:
        linhas = [l for l in f if not l.startswith("#")]
    return list(csv.DictReader(linhas))


def curvas(arq):
    out = {}
    for r in le_csv(os.path.join(DADOS, arq)):
        out.setdefault((r["figura"], r["curva"]), []).append((float(r["x"]), float(r["y"])))
    return {k: np.array(v) for k, v in out.items()}


def estilo(a, x, y, titulo):
    a.set_xlabel(x)
    a.set_ylabel(y)
    a.set_title(titulo, fontsize=9, loc="left")
    a.grid(color="0.9", lw=0.6)
    a.set_axisbelow(True)


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        return 1
    pasta = sys.argv[1]
    saida = os.path.join(AQUI, "..", "figuras", "v1_comparacao.png")
    if "--saida" in sys.argv:
        saida = sys.argv[sys.argv.index("--saida") + 1]
    perf, hist = {}, {}
    for r in le_csv(os.path.join(pasta, "perfis.csv")):
        perf.setdefault(r["caso"], {}).setdefault(float(r["t_s"]), []).append(
            [float(r[k]) for k in ("y_m", "pw_kPa", "pa_kPa", "Sw", "uy_m")])
    for r in le_csv(os.path.join(pasta, "historico.csv")):
        hist.setdefault(r["caso"], []).append([float(r[k]) for k in ("t_s", "qb_m_s", "w_topo_m", "pw_topo_kPa")])
    perf = {c: {t: np.array(v) for t, v in d.items()} for c, d in perf.items()}
    hist = {c: np.array(v) for c, v in hist.items()}
    f6 = curvas("laloui2003_fig6a_pa_30min.csv")
    f9 = curvas("laloui2003_fig9a_pa.csv")
    f8 = curvas("laloui2003_fig8b_pw_topo.csv")
    f5 = curvas("laloui2003_fig5.csv")
    fr = curvas("laloui2003_fig6b_8a_recalque.csv")
    tri = [c for c in ("trifasico", "trifasico_alert") if c in perf]

    fig, ax = plt.subplots(2, 3, figsize=(15, 9.5))
    # (a) Fig. 6a
    a = ax[0, 0]
    lo, hi = f6[("6a", "envoltoria min")], f6[("6a", "envoltoria max")]
    a.fill_betweenx(lo[:, 1], lo[:, 0], np.interp(lo[:, 1], hi[:, 1], hi[:, 0]), color="0.85", lw=0,
                    label="cinco equipes ALERT (envoltoria)")
    mins = [(k[1].split()[-1], v[0]) for k, v in f6.items() if k[1].startswith("minimo")]
    a.plot([m[1][0] for m in mins], [m[1][1] for m in mins], "k^", ms=6, label="minimos das equipes")
    for c in tri + [c for c in ("trifasico_laloui",) if c in perf]:
        if 1800.0 in perf[c]:
            p = perf[c][1800.0]
            a.plot(p[:, 2], p[:, 0], color=COR[c], **{"lw": 2, **TRACO.get(c, {})}, label=NOME[c])
    estilo(a, "p_a (kPa)", "altura (m)", "(a) p_a em 30 min: Fig. 6a")
    a.legend(fontsize=7, loc="lower left")
    # (b) Fig. 9a
    a = ax[0, 1]
    c = "trifasico_alert"
    tempos = [300, 600, 1800, 3600, 7200, 36000]
    tons = ["#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
    if c in perf:
        for t, cor in zip(tempos, tons):
            if float(t) in perf[c]:
                p = perf[c][float(t)]
                a.plot(p[:, 2], p[:, 0], color=cor, lw=2, label="%d min" % (t // 60))
    for k, v in f9.items():
        if k[1].startswith("minimo"):
            a.plot(v[0, 0], v[0, 1], "ko", ms=5)
            a.annotate(k[1].split("=")[1], (v[0, 0], v[0, 1]), fontsize=7, xytext=(3, 3), textcoords="offset points")
    for t in ("t=30 min", "t=600 min"):
        v = f9[("9a", t)]
        a.plot(v[:, 0], v[:, 1], "k--", lw=1)
    a.plot([], [], "ko", ms=5, label="Fig. 9a: minimos")
    a.plot([], [], "k--", lw=1, label="Fig. 9a: 30 e 600 min")
    estilo(a, "p_a (kPa)", "altura (m)", "(b) p_a no tempo, ar do benchmark: Fig. 9a")
    a.legend(fontsize=7, loc="lower left", ncol=2)
    # (c) p_w: Figs. 5a e 8b
    a = ax[0, 2]
    for c, ls in (("bifasico", "-"), ("trifasico_alert", "-"), ("trifasico", ":")):
        if c not in perf:
            continue
        for t in (300, 600, 1800, 36000):
            if float(t) in perf[c]:
                p = perf[c][float(t)]
                a.plot(p[:, 1], p[:, 0], color=COR[c], lw=1.6, ls=ls, label=NOME[c] if t == 300 else None)
    marc = {"5": "s", "10": "D", "20": "^", "30": "x", "600": "*"}
    for k, v in f5.items():
        if k[0] == "5a" and k[1].startswith("experimento"):
            tt = k[1].split("=")[1].split()[0]
            a.plot(v[:, 1], v[:, 0], ls="none", marker=marc[tt], mfc="none", mec=TINTA, ms=6,
                   label="experimento %s min" % tt)
    for k, m in ((("8b", "duas fases"), "<"), (("8b", "tres fases"), ">")):
        v = f8[k]
        a.plot(v[:, 1], np.ones(len(v)), ls="none", marker=m, color=CINZA, ms=7, label="Fig. 8b topo, " + k[1])
    estilo(a, "p_w (kPa)", "altura (m)", "(c) p_w em 5, 10, 30 e 600 min: Figs. 5a e 8b")
    a.legend(fontsize=6, loc="lower left", ncol=2)
    # (d) vazao
    a = ax[1, 0]
    for c in hist:
        h = hist[c]
        a.plot(h[1:, 0] / 60, h[1:, 1] * 1e6, color=COR.get(c, "k"), **{"lw": 2, **TRACO.get(c, {})}, label=NOME.get(c, c))
    v = f5[("5b", "experimento")]
    a.plot(v[:, 0], v[:, 1], "ks", ms=5, mfc="none", label="experimento (Fig. 5b)")
    v = f5[("5b", "tres fases")]
    a.plot(v[:, 0], v[:, 1], "k--", lw=1, label="tres fases de Laloui (Fig. 5b)")
    a.set_xlim(0, 200)
    estilo(a, "t (min; o eixo da Fig. 5b diz s)", "q_b (1e-6 m/s)", "(d) vazao de agua na base: Fig. 5b")
    a.legend(fontsize=7)
    # (e) recalque
    a = ax[1, 1]
    for c in hist:
        h = hist[c]
        a.plot(h[:, 0] / 60, h[:, 2] * 1e3, color=COR.get(c, "k"), **{"lw": 2, **TRACO.get(c, {})}, label=NOME.get(c, c))
    lo, hi = fr[("6b", "4 equipes min")], fr[("6b", "4 equipes max")]
    a.fill_between(lo[:, 0], lo[:, 1], hi[:, 1], color="0.8", label="Fig. 6b: EPFL, LCPC, LIEGE, UPD")
    v = fr[("6b", "UPC")]
    a.plot(v[:, 0], v[:, 1], "k^-", ms=5, lw=1, label="Fig. 6b: UPC")
    v = fr[("8a", "topo")]
    a.plot(v[:, 0], v[:, 1], "kD", ms=6, label="Fig. 8a (600 min)")
    estilo(a, "t (min)", "w do topo (mm)", "(e) recalque do topo: Figs. 6b e 8a")
    a.legend(fontsize=7)
    # (f) p_w no topo no tempo
    a = ax[1, 2]
    for c in hist:
        h = hist[c]
        a.plot(h[:, 0] / 60, h[:, 3], color=COR.get(c, "k"), **{"lw": 2, **TRACO.get(c, {})}, label=NOME.get(c, c))
    for k, m in ((("8b", "duas fases"), "o"), (("8b", "tres fases"), "s")):
        v = f8[k]
        a.plot(v[:, 0], v[:, 1], ls="none", marker=m, mfc="none", mec=TINTA, ms=7, label="Fig. 8b, " + k[1])
    a.set_xscale("log")
    estilo(a, "t (min)", "p_w no topo (kPa)", "(f) p_w no topo: Fig. 8b")
    a.legend(fontsize=7)
    fig.tight_layout()
    os.makedirs(os.path.dirname(os.path.abspath(saida)), exist_ok=True)
    fig.savefig(saida, dpi=130)
    print("figura:", saida)

    # tabela artigo x nosso
    print("\n%-34s %-12s %s" % ("grandeza", "artigo", "  ".join("%-16s" % c for c in perf)))

    def linha(nome, art, f):
        vals = []
        for c in perf:
            try:
                vals.append("%-16s" % f(c))
            except (KeyError, IndexError):
                vals.append("%-16s" % "-")
        print("%-34s %-12s %s" % (nome, art, "  ".join(vals)))

    def pamin(c, t):
        p = perf[c][float(t)]
        k = np.argmin(p[:, 2])
        return "%.2f @ %.3f" % (p[k, 2], p[k, 0])

    for t, art in ((300, "-4.65 @ 0.94"), (600, "-6.17 @ 0.88"), (1800, "-6.80 @ 0.82"), (3600, "-6.25 @ 0.75"),
                   (7200, "-5.25 @ 0.65"), (36000, "-1.89 @ 0.42")):
        linha("min p_a %d min (Fig. 9a)" % (t // 60), art, lambda c, t=t: pamin(c, t))
    for t, a2, a3 in ((300, -3.85, -5.70), (600, -5.03, -7.39), (1800, -6.50, -8.56), (7200, -8.72, -9.26),
                      (36000, -9.66, -9.55)):
        linha("p_w topo %d min (8b: 2f / 3f)" % (t // 60), "%.2f / %.2f" % (a2, a3),
              lambda c, t=t: "%.2f" % perf[c][float(t)][-1, 1])
    linha("w topo 600 min (6b / 8a)", "-3.3 / -1.6", lambda c: "%.2f" % (perf[c][36000.0][-1, 4] * 1e3))
    return 0


if __name__ == "__main__":
    sys.exit(main())
