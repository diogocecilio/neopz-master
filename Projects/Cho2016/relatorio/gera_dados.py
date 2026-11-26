#!/usr/bin/env python3
"""gera_dados.py: as figuras e tabelas do relatorio do Cho2016 (relatorio_cho2016.tex), com os dados que ja existem em
~/projects/runs/cho2016 (as rodadas que faltam ficam de fora, sem erro):
  figs/fs_silte.pdf     o FS contra o tempo no silte, sp e tp: o nosso (nref = 2, malha media; e o refinamento ate
                        convergir, malha grossa) contra a Fig. 12b do artigo (lida);
  figs/fs_areia.pdf     o mesmo na areia (por enquanto so o artigo e o FS inicial);
  figs/perfis_silte.pdf os perfis de p_w na secao C-C' (x = 28 m) a 32, 64, 96, 128 e 160 h, com a frente da Fig. 8;
  figs/malhas.pdf       a malha do transiente e as malhas iniciais do fator;
  tab/fs_silte.tex      o FS nosso / artigo a cada 8 h;
  tab/niveis.tex        o FS em cada nivel do refinamento adaptativo;
  tab/malha_transiente.tex  a profundidade da frente em tres malhas do transiente."""
import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.collections import PolyCollection  # noqa: E402

AQUI = os.path.dirname(os.path.abspath(__file__))
R = os.path.expanduser("~/projects/runs/cho2016")
P2 = os.path.join(R, "passo2")
P1 = os.path.join(R, "passo1")
P3 = os.path.join(R, "passo3")  # a malha fina do transiente (dx 0,5 m, 1a camada 0,5 %), associado e psi ~ 0
DATA = os.path.join(AQUI, "..", "data")
TINTA, TINTA2, GRADE = "#0b0b0b", "#5c5c58", "#e4e3df"
COR = {"sp": "#e34948", "tp": "#2a78d6"}
NOME = {"sp": "sp (só água)", "tp": "tp (água e ar)"}
os.makedirs(os.path.join(AQUI, "figs"), exist_ok=True)
os.makedirs(os.path.join(AQUI, "tab"), exist_ok=True)
v = lambda x: f"{x:.3f}".replace(".", "{,}")
v2 = lambda x: f"{x:.2f}".replace(".", "{,}")


def le(arq):
    return list(csv.DictReader(open(arq))) if os.path.exists(arq) else []


def comentado(arq):
    L = [r for r in csv.reader(open(arq)) if r and not r[0].startswith("#")]
    return [dict(zip(L[0], r)) for r in L[1:]]


lido = comentado(os.path.join(DATA, "fig12_lido.csv"))
frente = comentado(os.path.join(DATA, "fig8_frente_lido.csv"))


def artigo(solo, m):
    L = [r for r in lido if r["solo"] == solo and r["modelo"] == m]
    return [float(r["t_h"]) for r in L], [float(r["FS"]) for r in L]


def eixo(ax):
    ax.grid(True, color=GRADE, lw=0.7)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.tick_params(labelsize=8, colors=TINTA2)


# ---- FS contra o tempo
def fig_fs(solo, tfim, saida):
    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    eixo(ax)
    for m in ("sp", "tp"):
        ta, fa = artigo(solo, m)
        ax.plot(ta, fa, "o" if m == "sp" else "s", ms=4.5, mfc="none", mec=COR[m], mew=1.1,
                label=f"{NOME[m]}: Cho (2016), Fig. 12")
        N = le(os.path.join(P2, f"cho_{solo}_{m}_fs.csv"))
        if N:
            ax.plot([float(r["t_h"]) for r in N], [float(r["FS"]) for r in N], color=COR[m], lw=1.8,
                    label=f"{NOME[m]}: NeoPZ, malha média, nref = 2")
        for pre, ls, rot in (("fina_assoc", (0, (5, 2)), "malha fina do transiente, associado"),
                             ("fina_psi1e-6", (0, (1, 1.3)), "malha fina do transiente, $\\psi \\approx 0$")):
            F = le(os.path.join(P3, f"{pre}_{solo}_{m}_fs.csv"))
            if F:
                ax.plot([float(r["t_h"]) for r in F], [float(r["FS"]) for r in F], color=COR[m], lw=1.4, ls=ls,
                        label=f"{NOME[m]}: NeoPZ, {rot}")
        C = le(os.path.join(P2, f"conv_{solo}_{m}_fs.csv"))
        if C:
            ax.plot([float(r["t_h"]) for r in C], [float(r["FS"]) for r in C], "^" if m == "sp" else "v",
                    color=COR[m], ms=6, ls="none", label=f"{NOME[m]}: NeoPZ, malha grossa, refinado até convergir")
    ax.set_xlim(0, tfim)
    ax.set_xlabel("tempo de chuva (h)", fontsize=9, color=TINTA)
    ax.set_ylabel("fator de segurança (SRM)", fontsize=9, color=TINTA)
    ax.legend(fontsize=7, frameon=False, loc="lower left")
    fig.tight_layout()
    fig.savefig(saida)
    plt.close(fig)


fig_fs("silte", 168, os.path.join(AQUI, "figs", "fs_silte.pdf"))
fig_fs("areia", 30, os.path.join(AQUI, "figs", "fs_areia.pdf"))

# ---- tabela do FS no silte, a cada 8 h
linhas = []
N = {m: {float(r["t_h"]): float(r["FS"]) for r in le(os.path.join(P2, f"cho_silte_{m}_fs.csv"))} for m in ("sp", "tp")}
C = {m: {float(r["t_h"]): float(r["FS"]) for r in le(os.path.join(P2, f"conv_silte_{m}_fs.csv"))} for m in ("sp", "tp")}
A = {m: dict(zip(*artigo("silte", m))) for m in ("sp", "tp")}
with open(os.path.join(AQUI, "tab", "fs_silte.tex"), "w") as f:
    f.write("\\begin{tabular}{r|cccc|cccc}\n\\toprule\n")
    f.write(" & \\multicolumn{4}{c|}{sp (só água)} & \\multicolumn{4}{c}{tp (água e ar)}\\\\\n")
    f.write("$t$ (h) & nref 2 & conv. & Cho & dif. & nref 2 & conv. & Cho & dif.\\\\\n\\midrule\n")
    for t in sorted(set(N["sp"]) | set(N["tp"]) | set(C["sp"]) | set(C["tp"])):
        cel = []
        for m in ("sp", "tp"):
            n, c, a = N[m].get(t), C[m].get(t), A[m].get(t)
            melhor = c if c is not None else n
            cel += [v(n) if n is not None else "--", v(c) if c is not None else "--", v2(a) if a is not None else "--",
                    (f"{melhor - a:+.2f}".replace(".", "{,}") if melhor is not None and a is not None else "--")]
        f.write(f"{t:g} & " + " & ".join(cel) + "\\\\\n")
    f.write("\\bottomrule\n\\end{tabular}\n")

# ---- tabela dos niveis do refinamento
with open(os.path.join(AQUI, "tab", "niveis.tex"), "w") as f:
    f.write("\\begin{tabular}{llr|l}\n\\toprule\n")
    f.write("rodada & modelo & $t$ (h) & FS em cada nível (malha inicial, 1º refinamento, ...)\\\\\n\\midrule\n")
    for pre, nome in (("cho", "média, nref 2"), ("conv", "grossa, até convergir")):
        for m in ("sp", "tp"):
            for r in le(os.path.join(P2, f"{pre}_silte_{m}_fs.csv")):
                niv = [float(x) for x in r["niveis"].split(";")]
                f.write(f"{nome} & {m} & {float(r['t_h']):g} & " + "; ".join(v(x) for x in niv) + "\\\\\n")
    f.write("\\bottomrule\n\\end{tabular}\n")


# ---- perfis de p_w em C-C'
def frente_de(R_, t):
    """profundidade (m, abaixo de y = 22) em que p_w cruza -15 kPa, descendo da superficie"""
    P = sorted(((float(r["y_m"]), float(r["pw_kPa"])) for r in R_ if float(r["t_h"]) == t), reverse=True)
    for (y0, p0), (y1, p1) in zip(P, P[1:]):
        if p0 >= -15. > p1:
            return 22. - (y0 + (y1 - y0) * (-15. - p0) / (p1 - p0))
    return None


HORAS = [32, 64, 96, 128, 160]
CORH = plt.cm.viridis([i / (len(HORAS) - 1) for i in range(len(HORAS))])
fig, axs = plt.subplots(1, 2, figsize=(7.4, 4.6), sharey=True)
for ax, m in zip(axs, ("tp", "sp")):
    eixo(ax)
    Rm = le(os.path.join(P2, f"cho_silte_{m}_cc.csv"))
    for h, cor in zip(HORAS, CORH):
        P = [r for r in Rm if float(r["t_h"]) == h]
        if P:
            ax.plot([float(r["pw_kPa"]) for r in P], [22. - float(r["y_m"]) for r in P], color=cor, lw=1.5, label=f"{h} h")
        F = [float(r["frente_m"]) for r in frente if r["modelo"] == m and float(r["t_h"]) == h]
        if F:
            ax.plot([-15.], F, "D", ms=6, mfc="none", mec=cor, mew=1.3)
    ax.set_xlim(-32, 10)
    ax.set_ylim(10, 0)
    ax.set_title(NOME[m], fontsize=9, loc="left", color=TINTA)
    ax.set_xlabel("$p_w$ (kPa)", fontsize=9, color=TINTA)
axs[0].set_ylabel("profundidade em C–C' (m)", fontsize=9, color=TINTA)
axs[1].plot([], [], "D", ms=6, mfc="none", mec=TINTA2, label="frente na Fig. 8 do artigo")
axs[1].legend(fontsize=7, frameon=False, loc="lower left")
fig.tight_layout()
fig.savefig(os.path.join(AQUI, "figs", "perfis_silte.pdf"))
plt.close(fig)

with open(os.path.join(AQUI, "tab", "frente.tex"), "w") as f:
    f.write("\\begin{tabular}{r|cc|cc}\n\\toprule\n & \\multicolumn{2}{c|}{sp} & \\multicolumn{2}{c}{tp}\\\\\n")
    f.write("$t$ (h) & NeoPZ & Cho & NeoPZ & Cho\\\\\n\\midrule\n")
    for h in HORAS:
        cel = []
        for m in ("sp", "tp"):
            z = frente_de(le(os.path.join(P2, f"cho_silte_{m}_cc.csv")), h)
            a = [float(r["frente_m"]) for r in frente if r["modelo"] == m and float(r["t_h"]) == h]
            cel += [v2(z) if z is not None else "--", v2(a[0]) if a else "--"]
        f.write(f"{h} & " + " & ".join(cel) + "\\\\\n")
    f.write("\\bottomrule\n\\end{tabular}\n")

# ---- frente em tres malhas do transiente (sp, silte)
MT = os.path.join(R, "malha_transiente")
with open(os.path.join(AQUI, "tab", "malha_transiente.tex"), "w") as f:
    f.write("\\begin{tabular}{lcc|cccc}\n\\toprule\n")
    f.write("malha & colunas & 1ª camada & \\multicolumn{4}{c}{frente (m) a 16, 32, 48, 64 h}\\\\\n\\midrule\n")
    for nome, dx, e in (("grossa", "2 m", "2\\,\\%"), ("atual", "1 m", "1\\,\\%"), ("fina", "0,5 m", "0,5\\,\\%")):
        Rm = le(os.path.join(MT, f"{nome}_silte_sp_cc.csv"))
        z = [frente_de(Rm, t) for t in (16, 32, 48, 64)]
        f.write(f"{nome} & {dx} & {e} & " + " & ".join(v2(x) if x is not None else "--" for x in z) + "\\\\\n")
    f.write("\\bottomrule\n\\end{tabular}\n")

# ---- malhas
M = os.path.join(R, "malhas")
fig, axs = plt.subplots(3, 1, figsize=(7.4, 8.6))
for ax, (arq, tit) in zip(axs, (("cho_transiente.csv", "(a) transiente: colunas de 1 m, 1ª camada com 1 % da altura"),
                               ("cho_fator_media.csv", "(b) fator, malha inicial média (nref = 2)"),
                               ("cho_fator_grossa.csv", "(c) fator, malha inicial grossa (refinada até convergir)"))):
    T = le(os.path.join(M, arq))
    pol = [[(float(r["x1"]), float(r["y1"])), (float(r["x2"]), float(r["y2"])), (float(r["x3"]), float(r["y3"]))]
           for r in T]
    ax.add_collection(PolyCollection(pol, facecolors="none", edgecolors=TINTA2, linewidths=0.25))
    ax.plot([0, 60], [6, 6], color="#2a78d6", lw=1.0, ls=(0, (4, 2)))
    for nome, (x, y) in {"A": (0, 28), "B": (20, 26), "C": (28, 22), "D": (40, 12), "E": (60, 10), "F": (60, 0),
                         "G": (0, 0)}.items():
        ax.annotate(nome, (x, y), xytext=(3, 3), textcoords="offset points", fontsize=7, color=TINTA)
    ax.plot([28, 28], [0, 22], color="#e34948", lw=0.9)
    ax.set_xlim(-1, 62)
    ax.set_ylim(-1, 30)
    ax.set_aspect("equal")
    ax.set_title(f"{tit}: {len(T)} triângulos", fontsize=8, loc="left", color=TINTA)
    ax.tick_params(labelsize=7, colors=TINTA2)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
fig.tight_layout()
fig.savefig(os.path.join(AQUI, "figs", "malhas.pdf"))
plt.close(fig)
# ---- pontos H, I e J (Fig. 5b): o NeoPZ da malha fina (associado; o transiente nao depende do psi) contra o artigo
L5 = comentado(os.path.join(DATA, "fig5_lido.csv"))
fig, axs = plt.subplots(1, 3, figsize=(7.6, 3.2), sharey=True)
for ax, pt in zip(axs, "HIJ"):
    eixo(ax)
    for m in ("sp", "tp"):
        A_ = [r for r in L5 if r["modelo"] == m and r["ponto"] == pt]
        ax.plot([float(r["t_h"]) for r in A_], [float(r["pw_kPa"]) for r in A_], "o" if m == "sp" else "s", ms=3.5,
                mfc="none", mec=COR[m], mew=1.0, label=f"{m}: Cho (2016), Fig. 5b")
        for pre in ("fina_assoc", "fina_psi1e-6"):
            N_ = [r for r in le(os.path.join(P3, f"{pre}_silte_{m}_pontos.csv")) if r["ponto"] == pt]
            if N_:
                ax.plot([float(r["t_h"]) for r in N_], [float(r["pw_kPa"]) for r in N_], color=COR[m], lw=1.5,
                        label=f"{m}: NeoPZ, malha fina")
                if m == "tp":
                    ax.plot([float(r["t_h"]) for r in N_], [float(r["pa_kPa"]) for r in N_], color=COR[m], lw=1.0,
                            ls=(0, (3, 1.5)), label="tp: NeoPZ, $p_a$")
                break
    ax.set_title(f"ponto {pt}", fontsize=9, loc="left", color=TINTA)
    ax.set_xlim(0, 170)
    ax.set_ylim(-32, 8)
    ax.set_xlabel("t (h)", fontsize=8, color=TINTA)
axs[0].set_ylabel("$p_w$ (kPa)", fontsize=9, color=TINTA)
h_, l_ = axs[0].get_legend_handles_labels()
fig.legend(h_, l_, loc="upper center", ncol=3, fontsize=7, frameon=False, bbox_to_anchor=(0.5, 1.02))
fig.tight_layout(rect=(0, 0, 1, 0.86))
fig.savefig(os.path.join(AQUI, "figs", "pontos_silte.pdf"))
plt.close(fig)

# ---- dilatancia: o FS inicial com o PV2 (associado) e o PV com psi = phi, 0 e 1e-6 graus
PS = os.path.join(R, "psi")
with open(os.path.join(AQUI, "tab", "psi.tex"), "w") as f:
    f.write("\\begin{tabular}{ll|l|c}\n\\toprule\n")
    f.write("solo & material & FS em cada nível (malha média, 2 níveis) & Cho\\\\\n\\midrule\n")
    for so, cho in (("silte", "1,55"), ("areia", "1,22")):
        lin = [("PV2, $\\psi = \\phi'$", os.path.join(P1, f"cho_{so}_sp_fs.csv")),
               ("PV, $\\psi = \\phi'$", os.path.join(PS, f"assoc_{so}_sp_fs.csv")),
               ("PV, $\\psi = 0$", os.path.join(PS, f"psi0_{so}_sp_fs.csv")),
               ("PV, $\\psi = 10^{-6}$ graus", os.path.join(PS, f"psi1e-6_{so}_sp_fs.csv"))]
        for nome, arq in lin:
            Rr = le(arq)
            niv = "; ".join(v(float(x)) for x in Rr[0]["niveis"].split(";")) if Rr else "(rodando)"
            f.write(f"{so} & {nome} & {niv} & {cho}\\\\\n")
    f.write("\\bottomrule\n\\end{tabular}\n")
print("pronto")
