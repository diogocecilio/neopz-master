"""Pós-processamento do LearningLiakopoulos: a coluna de Liakopoulos do NeoPZ contra a referência 1D
(coluna_liakopoulos_1d.py), a solução analítica da coluna saturada (Eq. 48 do resumo) e as curvas calculadas de
Callari e Abati (2009, Figs. 3 e 4), digitalizadas em data/callari2009_fig3_fig4_digitalizado.csv.

Uso (na pasta onde o programa gravou os CSVs; o CMake copia para lá este script, a referência 1D e os dados):
    python3 plot_liakopoulos.py [historico.csv] [perfis.csv] [liakopoulos.png]

Saídas:
  * a figura (como a Fig. 5 do resumo): vazão na base, perfis de p_w, recalque do topo e perfis de recalque;
  * comparacao.csv: a tabela da Seção 10.6, NeoPZ ao lado da referência 1D e da solução analítica;
  * no terminal: (1) a tabela da Seção 10.6; (2) a verificação de cada caso contra a referência 1D com os mesmos
    passos de tempo; (3) a convergência no tempo (metade do passo); (4) a diferença para as curvas do artigo.
Os casos ausentes dos CSVs (o programa pode ter rodado só parte deles) são ignorados.
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)
import coluna_liakopoulos_1d as ref1d  # noqa: E402  (referência 1D: mesmas leis e hipóteses do código)

hist_file = sys.argv[1] if len(sys.argv) > 1 else "historico.csv"
perf_file = sys.argv[2] if len(sys.argv) > 2 else "perfis.csv"
out = sys.argv[3] if len(sys.argv) > 3 else "liakopoulos.png"

hist = pd.read_csv(hist_file)
perf = pd.read_csv(perf_file)
art = pd.read_csv(os.path.join(AQUI, "data", "callari2009_fig3_fig4_digitalizado.csv"), comment="#")

T_TAB = (300, 600, 1200, 1800, 3600, 6400, 7500)          # tabela da Seção 10.6
T_PERF = (300, 600, 1200, 1800, 3600, 6400, 7200, 7500)   # perfis gravados pelo programa
T_PW = (600, 1200, 1800)                                  # Fig. 3b
T_UY = (300, 1200, 7200)                                  # Fig. 4b


def caso(nome):
    return hist[hist.caso == nome].sort_values("t_s")


def perfil(nome, t):
    d = perf[(perf.caso == nome) & (np.abs(perf.t_s - t) < 1e-6)]
    return d.sort_values("y_m")


# ---------------- referência 1D com os mesmos passos de cada caso do NeoPZ ----------------
meio = lambda t: 0.5*ref1d.passo_referencia(t)  # noqa: E731
CASOS = {  # nome no CSV: (rótulo, argumentos da referência 1D)
    "parcial": ("passos da referência", dict(dt_fn=ref1d.passo_referencia)),
    "parcial_meio_passo": ("metade do passo", dict(dt_fn=meio)),
    "parcial_passos_artigo": ("passos do artigo", dict(dt_fn=ref1d.passo_artigo)),
    "saturada": ("saturada", dict(dt_fn=ref1d.passo_referencia, saturado=True)),
}
presentes = [c for c in CASOS if not caso(c).empty]
if not presentes:
    sys.exit(f"nenhum caso conhecido em {hist_file}")
ref = {}
for c in presentes:
    y1, h1, p1 = ref1d.resolve(t_perfis=T_PERF, **CASOS[c][1])
    ref[c] = (y1, h1, p1)


def ref_em(c, col, t):
    """Histórico da referência 1D do caso c no instante t: q_b (m/s) ou w (m)."""
    _, h1, _ = ref[c]
    return np.interp(t, h1["t"], h1[col])


def neo_em(c, col, t):
    d = caso(c)
    return np.interp(t, d.t_s.values, d[col].values)


# ---------------- (1) tabela da Seção 10.6 ----------------
def tabela_secao_10_6():
    linhas = []
    print("\n(1) Tabela da Seção 10.6: NeoPZ (passos da referência) [referência 1D]; saturada: NeoPZ [Eq. 48]")
    print(f"{'t (s)':>6} | {'q_b (1e-6 m/s)':>17} | {'p_w y=0,5 m (kPa)':>17} | {'p_w topo (kPa)':>17} | "
          f"{'w (mm)':>17} || {'q_b sat (1e-6 m/s)':>18} | {'w sat (mm)':>17}")
    for t in T_TAB:
        lin = {"t_s": t}
        txt = f"{t:6d} |"
        if "parcial" in presentes:
            d = perfil("parcial", t)
            y1, h1, p1 = ref["parcial"]
            ymeio = d.y_m.values[np.argmin(np.abs(d.y_m.values - 0.5))]
            lin.update(qb=1e6*neo_em("parcial", "qb_m_s", t), qb_ref=1e6*ref_em("parcial", "q", t),
                       pw_meio=np.interp(ymeio, d.y_m, d.pw_kPa), pw_meio_ref=np.interp(ymeio, y1, p1[t][0]),
                       pw_topo=d.pw_kPa.values[-1], pw_topo_ref=p1[t][0][-1],
                       w_mm=1e3*neo_em("parcial", "w_topo_m", t), w_mm_ref=1e3*ref_em("parcial", "w", t))
            txt += (f" {lin['qb']:7.3f} [{lin['qb_ref']:7.3f}] | {lin['pw_meio']:7.2f} [{lin['pw_meio_ref']:7.2f}] |"
                    f" {lin['pw_topo']:7.2f} [{lin['pw_topo_ref']:7.2f}] |"
                    f" {lin['w_mm']:7.3f} [{lin['w_mm_ref']:7.3f}] ||")
        else:
            txt += f"{'':>79}||"
        if "saturada" in presentes:
            wa, qa = ref1d.saturado_analitico(t)
            lin.update(qb_sat=1e6*neo_em("saturada", "qb_m_s", t), qb_sat_analitica=1e6*qa,
                       w_sat_mm=1e3*neo_em("saturada", "w_topo_m", t), w_sat_mm_analitica=1e3*wa)
            txt += (f" {lin['qb_sat']:8.3f} [{lin['qb_sat_analitica']:7.3f}] |"
                    f" {lin['w_sat_mm']:7.3f} [{lin['w_sat_mm_analitica']:7.3f}]")
        linhas.append(lin)
        print(txt)
    pd.DataFrame(linhas).to_csv("comparacao.csv", index=False, float_format="%.6g")


if {"parcial", "saturada"} & set(presentes):
    tabela_secao_10_6()

# ---------------- (2) verificação: cada caso contra a referência 1D com os mesmos passos ----------------
print("\n(2) NeoPZ contra a referência 1D com os mesmos passos de tempo (máximo em todos os passos e perfis)")
print(f"{'caso':>22} | {'passos':>6} | {'max |dq_b|/k_s':>14} | {'max |dw| (mm)':>13} | {'max |dp_w| (kPa)':>16} |"
      f" {'Newton/passo':>12} | {'pontos plást.':>13}")
for c in presentes:
    d = caso(c)
    y1, h1, p1 = ref[c]
    dq = np.max(np.abs(d.qb_m_s.values - np.interp(d.t_s.values, h1["t"], h1["q"])))/ref1d.ks
    dw = 1e3*np.max(np.abs(d.w_topo_m.values - np.interp(d.t_s.values, h1["t"], h1["w"])))
    dp = 0.
    for t in T_PERF:
        e = perfil(c, t)
        if t in p1 and not e.empty:
            dp = max(dp, np.max(np.abs(e.pw_kPa.values - np.interp(e.y_m.values, y1, p1[t][0]))))
    passos = d[d.passo > 0]
    print(f"{c:>22} | {len(passos):6d} | {100*dq:13.3f}% | {dw:13.4f} | {dp:16.4f} | {passos.newton.mean():12.2f} |"
          f" {int(d.pontos_plasticos.max()):13d}")
if "saturada" in presentes:
    d = caso("saturada")
    an = np.array([ref1d.saturado_analitico(t) for t in d.t_s.values])
    print(f"{'saturada x Eq. (48)':>22} | {'':6} | {100*np.max(np.abs(d.qb_m_s.values - an[:, 1]))/ref1d.ks:13.3f}% |"
          f" {1e3*np.max(np.abs(d.w_topo_m.values - an[:, 0])):13.4f} |")

# ---------------- (3) convergência no tempo ----------------
if {"parcial", "parcial_meio_passo"} <= set(presentes):
    print("\n(3) Convergência no tempo: passos da referência contra metade deles (diferença relativa)")
    print(f"{'t (s)':>6} | {'q_b':>8} | {'w':>8} | {'p_w topo':>8}")
    for t in T_TAB:
        a, b = perfil("parcial", t), perfil("parcial_meio_passo", t)
        rq = neo_em("parcial", "qb_m_s", t)/neo_em("parcial_meio_passo", "qb_m_s", t) - 1
        rw = neo_em("parcial", "w_topo_m", t)/neo_em("parcial_meio_passo", "w_topo_m", t) - 1
        rp = a.pw_kPa.values[-1]/b.pw_kPa.values[-1] - 1
        print(f"{t:6d} | {100*rq:7.3f}% | {100*rw:7.3f}% | {100*rp:7.3f}%")

# ---------------- (4) diferença para as curvas digitalizadas do artigo ----------------
print("\n(4) Diferença para as curvas de Callari e Abati (2009), digitalizadas (incerteza de leitura ~1% do eixo)")
print(f"{'':>34} | {'q_b t<2400 s':>12} | {'q_b t>=2400 s':>13} | {'p_w média':>9} | {'p_w |média|':>11} |"
      f" {'w |média|':>9} | {'w sat |méd.|':>12} | {'u_y |média|':>11}")
q3a = art[art.figura == "3a"]
p3b = art[art.figura == "3b"]
w4a = art[(art.figura == "4a") & art.curva.str.contains("parcialmente")]
w4s = art[(art.figura == "4a") & art.curva.str.contains(r"\(saturada\)")]
u4b = art[art.figura == "4b"]


def tempo_da_curva(curva):
    return float(curva.split("t=")[1].split()[0])


def contra_artigo(fonte, nome):
    """fonte: função (col, t) -> valor no instante t e função t -> (y, p_w, u_y) do perfil (None se ausente)."""
    hist_em, perfil_em = fonte
    rq = 100*(hist_em("qb_m_s", q3a.x.values)*1e6/q3a.y.values - 1)
    cedo, tarde = q3a.x.values < 2400, q3a.x.values >= 2400
    dp, du = [], []
    for curva, g in p3b.groupby("curva"):
        pr = perfil_em(tempo_da_curva(curva))
        if pr is not None:
            dp += list(np.interp(g.y.values, pr[0], pr[1]) - g.x.values)
    for curva, g in u4b.groupby("curva"):
        pr = perfil_em(tempo_da_curva(curva))
        if pr is not None:
            du += list(1e3*np.interp(g.y.values, pr[0], pr[2]) - g.x.values)
    dw = 1e3*hist_em("w_topo_m", w4a.x.values) - w4a.y.values
    print(f"{nome:>34} | {np.mean(np.abs(rq[cedo])):11.1f}% | {np.mean(rq[tarde]):+12.1f}% |"
          f" {np.mean(dp):+9.3f} | {np.mean(np.abs(dp)):11.3f} | {np.mean(np.abs(dw)):9.3f} |", end="")
    return du


def fonte_neo(c):
    def h(col, t):
        return neo_em(c, col, t)

    def p(t):
        e = perfil(c, t)
        return None if e.empty else (e.y_m.values, e.pw_kPa.values, e.uy_m.values)
    return h, p


def fonte_ref(c):
    y1, h1, p1 = ref[c]
    col1 = {"qb_m_s": "q", "w_topo_m": "w"}

    def h(col, t):
        return np.interp(t, h1["t"], h1[col1[col]])

    def p(t):
        return (y1, p1[t][0], p1[t][1]) if t in p1 else None
    return h, p


for c, nome in (("parcial", "NeoPZ, passos da referência"), ("parcial_passos_artigo", "NeoPZ, passos do artigo")):
    if c not in presentes:
        continue
    du = contra_artigo(fonte_neo(c), nome)
    ws = (f"{np.mean(np.abs(1e3*neo_em('saturada', 'w_topo_m', w4s.x.values) - w4s.y.values)):12.3f}"
          if "saturada" in presentes else f"{'':>12}")
    print(f" {ws} | {np.mean(np.abs(du)):11.3f}  (u_y: {len(du)} pontos)")
if "parcial" in presentes:
    du = contra_artigo(fonte_ref("parcial"), "referência 1D, passos da referência")
    wa = np.array([ref1d.saturado_analitico(t)[0] for t in w4s.x.values])
    print(f" {np.mean(np.abs(1e3*wa - w4s.y.values)):12.3f} | {np.mean(np.abs(du)):11.3f}")
print("q_b: diferença relativa (NeoPZ/artigo - 1); p_w, w e u_y: diferença absoluta em kPa e mm "
      "(média com sinal: NeoPZ - artigo)")

# ---------------- figura (como a Fig. 5 do resumo) ----------------
# cores: categóricas na ordem fixa (casos) e rampa ordinal de um só azul (instantes dos perfis)
AZUL, LARANJA, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
RAMPA = ("#86b6ef", "#2a78d6", "#104281")
TINTA, TINTA2 = "#0b0b0b", "#52514e"
GRADE, EIXO, FUNDO = "#e1e0d9", "#c3c2b7", "#fcfcfb"
LW, LW_REF, MS = 1.6, 0.8, 5.5
TRACO = (0, (4, 3))
plt.rcParams.update({
    "figure.facecolor": FUNDO, "axes.facecolor": FUNDO, "savefig.facecolor": FUNDO,
    "axes.edgecolor": EIXO, "axes.linewidth": 0.8, "axes.labelcolor": TINTA2, "axes.titlecolor": TINTA,
    "axes.titlesize": 11, "axes.labelsize": 10, "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.color": GRADE, "grid.linewidth": 0.6, "grid.linestyle": "-",
    "xtick.color": TINTA2, "ytick.color": TINTA2, "xtick.labelsize": 9, "ytick.labelsize": 9,
    "legend.fontsize": 8.5, "legend.frameon": False, "legend.labelcolor": TINTA2,
})
ref_chave = Line2D([], [], color=TINTA, lw=LW_REF, ls=TRACO)
art_chave = Line2D([], [], ls="", marker="o", ms=MS, mfc=FUNDO, mec=TINTA2, mew=1.0)


def pontos(ax, x, y, cor=TINTA2):
    ax.plot(x, y, ls="", marker="o", ms=MS, mfc=FUNDO, mec=cor, mew=1.0, zorder=4)


fig, axs = plt.subplots(2, 2, figsize=(11, 8.6))
SERIES_T = (("parcial", AZUL, "NeoPZ, parcialmente saturada"),
            ("parcial_passos_artigo", LARANJA, "NeoPZ, passos do artigo"),
            ("saturada", AQUA, "NeoPZ, saturada"))

# (a) vazão na base
ax = axs[0, 0]
for c, cor, rot in SERIES_T:
    if c in presentes:
        d = caso(c)
        ax.plot(d.t_s, 1e6*d.qb_m_s, color=cor, lw=LW, label=rot, solid_capstyle="round")
        _, h1, _ = ref[c]
        ax.plot(h1["t"], 1e6*h1["q"], color=TINTA, lw=LW_REF, ls=TRACO, zorder=3)
pontos(ax, q3a.x, q3a.y)
ax.set(title="(a) vazão na base", xlabel="t (s)", ylabel="$q_b$ ($10^{-6}$ m/s)", xlim=(0, 7500), ylim=(0, None))
h, l = ax.get_legend_handles_labels()
ax.legend(h + [ref_chave, art_chave], l + ["referência 1D (mesmos passos)", "Callari e Abati (2009, Fig. 3a)"])

# (b) pressão da água
ax = axs[0, 1]
if "parcial" in presentes:
    y1, _, p1 = ref["parcial"]
    for t, cor in zip(T_PW, RAMPA):
        e = perfil("parcial", t)
        ax.plot(e.pw_kPa, e.y_m, color=cor, lw=LW, label=f"NeoPZ, {t} s")
        ax.plot(p1[t][0], y1, color=TINTA, lw=LW_REF, ls=TRACO, zorder=3)
        g = p3b[p3b.curva == f"p_w em t={t} s"]
        pontos(ax, g.x, g.y, cor)
ax.set(title="(b) pressão da água", xlabel="$p_w$ (kPa)", ylabel="y (m)", xlim=(-8, 0), ylim=(0, 1))
h, l = ax.get_legend_handles_labels()
ax.legend(h + [ref_chave, art_chave], l + ["referência 1D", "Callari e Abati (2009, Fig. 3b)"], loc="lower left")

# (c) recalque do topo
ax = axs[1, 0]
for c, cor, rot in SERIES_T:
    if c in presentes:
        d = caso(c)
        ax.plot(d.t_s, 1e3*d.w_topo_m, color=cor, lw=LW, label=rot, solid_capstyle="round")
        _, h1, _ = ref[c]
        ax.plot(h1["t"], 1e3*h1["w"], color=TINTA, lw=LW_REF, ls=TRACO, zorder=3)
if "saturada" in presentes:
    ta = np.linspace(0, 7500, 301)
    ax.plot(ta, [1e3*ref1d.saturado_analitico(t)[0] for t in ta], color=TINTA, lw=LW_REF, ls=(0, (1, 2)),
            zorder=3, label="saturada, analítica (Eq. 48)")
pontos(ax, w4a.x, w4a.y)
pontos(ax, w4s.x, w4s.y)
ax.set(title="(c) recalque do topo", xlabel="t (s)", ylabel="w (mm)", xlim=(0, 7500), ylim=(-2, 0))
h, l = ax.get_legend_handles_labels()
ax.legend(h + [ref_chave, art_chave], l + ["referência 1D (mesmos passos)", "Callari e Abati (2009, Fig. 4a)"],
          loc="upper right")

# (d) perfis de recalque
ax = axs[1, 1]
if "parcial" in presentes:
    y1, _, p1 = ref["parcial"]
    for t, cor in zip(T_UY, RAMPA):
        e = perfil("parcial", t)
        ax.plot(1e3*e.uy_m, e.y_m, color=cor, lw=LW, label=f"NeoPZ, {t} s")
        ax.plot(1e3*p1[t][1], y1, color=TINTA, lw=LW_REF, ls=TRACO, zorder=3)
        g = u4b[u4b.curva == f"u_y em t={t} s"]
        pontos(ax, g.x, g.y, cor)
ax.set(title="(d) perfis de recalque", xlabel="$u_y$ (mm)", ylabel="y (m)", xlim=(-1.7, 0), ylim=(0, 1))
h, l = ax.get_legend_handles_labels()
ax.legend(h + [ref_chave, art_chave], l + ["referência 1D", "Callari e Abati (2009, Fig. 4b)"], loc="lower left")

fig.suptitle("Dessaturação da coluna de Liakopoulos: NeoPZ (Taylor–Hood P2–P1, 320 triângulos) contra a referência"
             " 1D e o artigo", color=TINTA, fontsize=11.5)
fig.tight_layout()
fig.savefig(out, dpi=130)
print(f"\nfigura: {out}" + ("; tabela: comparacao.csv" if {"parcial", "saturada"} & set(presentes) else ""))
