#!/usr/bin/env python3
"""gera_dados.py: figuras e tabelas do relatorio det_3solos (os 3 solos do Rahardjo, 3 intensidades de chuva de 24 h,
variantes iii, ii e i, FS de hora em hora ate 48 h), a partir das rodadas do PC em ~/projects/runs/det_pc_3solos
(pastas <solo>_I<I> para a iii e <solo>_I<I>_<v> para a ii e a i: det_0_-1.csv, diag.csv, pw.txt, log.txt).
Saidas: figs/*.pdf, tab/*.tex e numeros.tex (macros com os numeros citados no texto).
Uso: python3 gera_dados.py [pasta das rodadas]"""
import csv
import math
import os
import re
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.ticker  # noqa: E402
import matplotlib.tri as mtri  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm  # noqa: E402
from matplotlib.path import Path  # noqa: E402

AQUI = os.path.dirname(os.path.abspath(__file__))
RUNS = sys.argv[1] if len(sys.argv) > 1 else os.path.expanduser("~/projects/runs/det_pc_3solos")
FIGS, TABS = os.path.join(AQUI, "figs"), os.path.join(AQUI, "tab")
os.makedirs(FIGS, exist_ok=True)
os.makedirs(TABS, exist_ok=True)

SOLOS = {"f10": (360, ["36", "90", "180"]), "f50": (36, ["9", "36", "90"]), "f100": (3.6, ["3.6", "9", "36"])}
VARS = [v for v in os.environ.get("VARIANTES", "iii,ii,i").split(",") if v]  # VARIANTES=iii: o relatorio parcial
NOME_VAR = {"iii": "(iii) desacoplado, ar passivo", "ii": "(ii) acoplado, ar passivo", "i": "(i) acoplado, ar ativo"}
# cores: rampa ordinal azul (300/500/700) para as intensidades; categoricas fixas (slots 1-3) + estilo para as variantes
RAMPA = ["#6da7ec", "#256abf", "#0d366b"]
COR_VAR = {"iii": "#2a78d6", "ii": "#eb6834", "i": "#1baf7a"}
ESTILO_VAR = {"iii": "-", "ii": (0, (5, 2)), "i": (0, (6, 2, 1.5, 2))}
TINTA, TINTA2, GRADE, SUPERF, FAIXA = "#0b0b0b", "#5c5c58", "#e4e3df", "#fcfcfb", "#f0efec"
# o artigo: Eq. 6 (FS minimo; os valores da Tabela B do relatorio de outubro) e a Fig. 5b (f50, digitalizada +-0,03)
EQ6 = {"f10": {"36": 1.48, "90": 1.05, "180": 0.85}, "f50": {"9": 1.79, "36": 1.05, "90": 0.82},
       "f100": {"3.6": 1.88, "9": 1.55, "36": 1.45}}
FIG5B = {"9": {6: 2.00, 12: 1.92, 24: 1.78, 48: 1.85}, "36": {6: 1.72, 12: 1.40, 24: 1.04, 48: 1.24},
         "90": {6: 1.35, 12: 1.18, 24: 0.86, 48: 1.15}}
# o relatorio de outubro (docs/nao_saturado, Tabela B: FS minimo; Tabela fig5b: completo / desacoplado)
OUT = {"f10": {"36": (1.58, 1.58, 1.57), "90": (1.15, 1.16, 1.15), "180": (0.86, 0.99, 1.00)},
       "f50": {"9": (1.85, 1.86, 1.74), "36": (1.18, 1.21, 1.01), "90": (0.85, 1.00, 1.00)},
       "f100": {"3.6": (1.94, 1.95, 1.77), "9": (1.69, 1.70, 1.36), "36": (1.48, 1.49, 1.28)}}  # (iii, ii, i)
OUT5B = {"9": {6: (2.03, 2.00), 12: (1.96, 1.90), 24: (1.85, 1.74), 48: (1.91, 1.87)},  # (iii, i)
         "36": {6: (1.82, 1.61), 12: (1.57, 1.18), 24: (1.18, 1.00), 48: (1.33, 1.51)},
         "90": {6: (1.44, 1.24), 12: (1.19, 1.06), 24: (0.85, 1.00), 48: (1.17, 1.51)}}
GEO = [(0, 0), (70, 0), (70, 20), (40, 20), (30, 30), (0, 30)]


def pasta(solo, I, v):
    return os.path.join(RUNS, f"{solo}_I{I}" + ("" if v == "iii" else f"_{v}"))


def le_fs(solo, I, v):
    arq = os.path.join(pasta(solo, I, v), "det_0_-1.csv")
    if not os.path.exists(arq):
        return None
    d = {}
    with open(arq) as f:
        for r in csv.DictReader(f):
            if r["situacao"] == "ok" and float(r["FS_min"]) > 0:
                d[float(r["t_h"])] = {k: r[k] for k in ("FS_min", "banda_prof", "banda_xmin", "banda_xmax", "nos_face")}
    return d or None


def le_diag(solo, I, v):
    arq = os.path.join(pasta(solo, I, v), "diag.csv")
    if not os.path.exists(arq):
        return None
    with open(arq) as f:
        return [{k: float(x) for k, x in r.items()} for r in csv.DictReader(f)]


def le_tempo(solo, I, v):
    arq = os.path.join(pasta(solo, I, v), "log.txt")
    if not os.path.exists(arq):
        return None
    m = re.search(r"realizacao\s+-1:.*\|\s*([\d.]+) s \(transiente ([\d.]+) s\)", open(arq).read())
    return (float(m.group(1)), float(m.group(2))) if m else None


def le_pw(solo, I, v, t):
    """os nos (x, y), p_w e p_a (so na i) no instante t"""
    arq = os.path.join(pasta(solo, I, v), "pw.txt")
    if not os.path.exists(arq):
        return None
    with open(arq) as f:
        linhas = f.read().split("\n")
    i = next(k for k, l in enumerate(linhas) if l.startswith("nos "))
    nn = int(linhas[i].split()[1])
    xy = [tuple(map(float, linhas[i + 1 + k].split())) for k in range(nn)]
    k = i + 1 + nn
    while k < len(linhas):
        if linhas[k].startswith("t "):
            tk = float(linhas[k].split()[1])
            pw = [float(x) for x in linhas[k + 1:k + 1 + nn]]
            pa = None
            k2 = k + 1 + nn
            if v == "i":
                pa = [float(x) for x in linhas[k2:k2 + nn]]
                k2 += nn
            if abs(tk - t) < 1e-9:
                return xy, pw, pa
            k = k2
        else:
            k += 1
    return None


def estilo(ax, titulo, ylabel="FS (SRM)"):
    ax.set_facecolor(SUPERF)
    ax.axvspan(0, 24, color=FAIXA, zorder=0)
    ax.set_xlim(0, 48)
    ax.set_xticks(range(0, 49, 6))
    ax.set_xlabel("t (h)", color=TINTA)
    ax.set_ylabel(ylabel, color=TINTA)
    ax.set_title(titulo, color=TINTA, fontsize=10, loc="left")
    ax.grid(True, color=GRADE, lw=0.8, zorder=0)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(GRADE)
    ax.tick_params(colors=TINTA2, labelsize=8)


def salva(fig, nome):
    fig.savefig(os.path.join(FIGS, nome + ".pdf"))
    fig.savefig(os.path.join(FIGS, nome + ".png"), dpi=130)
    plt.close(fig)


plt.rcParams.update({"font.size": 9, "savefig.facecolor": SUPERF, "figure.facecolor": SUPERF})
nv = len(VARS)
DADOS = {(s, I, v): le_fs(s, I, v) for s in SOLOS for I in SOLOS[s][1] for v in VARS}
DIAG = {(s, I, v): le_diag(s, I, v) for s in SOLOS for I in SOLOS[s][1] for v in VARS}


def serie(d):
    t = sorted(d)
    return t, [float(d[x]["FS_min"]) for x in t]


def cruza(d):
    """a primeira passagem: o instante (h) em que o FS cruza 1 descendo, interpolado linearmente entre os instantes
    calculados (0 se ja comeca abaixo; None se nunca cruza), e o instante em que volta a ficar >= 1 (None se nao volta)"""
    t, fs = serie(d)
    if fs[0] < 1:
        t1 = 0.
        i1 = 0
    else:
        i1 = next((i for i in range(1, len(fs)) if fs[i] < 1), None)
        if i1 is None:
            return None, None
        t1 = t[i1 - 1] + (fs[i1 - 1] - 1) / (fs[i1 - 1] - fs[i1]) * (t[i1] - t[i1 - 1])
    i2 = next((i for i in range(i1 + 1, len(fs)) if fs[i] >= 1), None)
    t2 = None if i2 is None else t[i2 - 1] + (1 - fs[i2 - 1]) / (fs[i2] - fs[i2 - 1]) * (t[i2] - t[i2 - 1])
    return t1, t2


def marca_cruza(ax, d, cor, dy=10):
    """o triangulo na linha FS = 1 no instante da primeira passagem, a linha vertical ate o eixo e o rotulo"""
    t1, _ = cruza(d)
    if t1 is None:
        return
    ax.axvline(t1, color=cor, lw=1, ls=(0, (2, 2)), zorder=2)
    ax.plot(t1, 1.0, "v", ms=10, color=cor, mec=SUPERF, mew=1.5, zorder=5)
    ax.annotate(f"$t_1$ = {t1:.1f} h".replace(".", ","), (t1, 1.0), textcoords="offset points", xytext=(-4, dy),
                ha="right", fontsize=8, color=TINTA)


def minimo(d):
    t, fs = serie(d)
    i = min(range(len(fs)), key=fs.__getitem__)
    return t[i], fs[i]


# ---- F1: FS x t por solo, um painel por variante, as 3 intensidades
for solo, (ks, Is) in SOLOS.items():
    fig, axs = plt.subplots(1, len(VARS), figsize=(5 * len(VARS) + 2 * (len(VARS) == 1), 4.3), sharey=True,
                            squeeze=False)
    for ax, v in zip(axs[0], VARS):
        estilo(ax, NOME_VAR[v], "FS (SRM)" if v == "iii" else "")
        ax.axhline(1.0, color=TINTA2, lw=1, ls=(0, (4, 3)), zorder=1)
        hs = []
        for I, cor in zip(Is, RAMPA):
            d = DADOS[(solo, I, v)]
            if not d:
                continue
            t, fs = serie(d)
            h, = ax.plot(t, fs, color=cor, lw=2, zorder=3, label=f"{I} mm/h ({float(I) / ks:.2g} $k_s$)")
            hs.append(h)
            tm, fm = minimo(d)
            ax.plot(tm, fm, "o", ms=8, color=cor, mec=SUPERF, mew=2, zorder=4)
            ax.annotate(f"{fm:.2f}", (tm, fm), textcoords="offset points", xytext=(5, -11), fontsize=8, color=TINTA)
            ax.plot(24, EQ6[solo][I], "D", ms=7, mfc="none", mec=cor, mew=1.5, zorder=4)
            marca_cruza(ax, d, cor)
            if solo == "f50":
                tp = sorted(FIG5B[I])
                ax.plot(tp, [FIG5B[I][x] for x in tp], "o", ms=7, mfc="none", mec=cor, mew=1.5, zorder=4)
        extra = [plt.Line2D([], [], ls="none", marker="D", ms=6, mfc="none", mec=TINTA2, label="Eq. 6 (FS mín.)"),
                 plt.Line2D([], [], ls="none", marker="v", ms=7, color=TINTA2, label="$t_1$: FS cruza 1")]
        if solo == "f50":
            extra.append(plt.Line2D([], [], ls="none", marker="o", ms=6, mfc="none", mec=TINTA2, label="Fig. 5b"))
        ax.legend(handles=hs + extra, frameon=False, fontsize=7.5, labelcolor=TINTA, loc="upper right")
    fig.suptitle(f"{solo} ($k_s$ = {ks:g} mm/h): FS × t, chuva de 24 h (faixa cinza)", color=TINTA, x=0.01, ha="left")
    fig.tight_layout()
    salva(fig, f"fs_t_{solo}")

# ---- F2: FS x t por solo, um painel por intensidade, as 3 variantes (so com mais de uma variante)
for solo, (ks, Is) in (SOLOS.items() if len(VARS) > 1 else []):
    fig, axs = plt.subplots(1, 3, figsize=(15, 4.3), sharey=True)
    for ax, I in zip(axs, Is):
        estilo(ax, f"{I} mm/h ({float(I) / ks:.2g} $k_s$)", "FS (SRM)" if I == Is[0] else "")
        ax.axhline(1.0, color=TINTA2, lw=1, ls=(0, (4, 3)), zorder=1)
        for v in VARS:
            d = DADOS[(solo, I, v)]
            if not d:
                continue
            t, fs = serie(d)
            ax.plot(t, fs, color=COR_VAR[v], ls=ESTILO_VAR[v], lw=2, zorder=3, label=NOME_VAR[v])
            marca_cruza(ax, d, COR_VAR[v], dy=10 + 12 * VARS.index(v))
        ax.plot(24, EQ6[solo][I], "D", ms=7, mfc="none", mec=TINTA2, mew=1.5, zorder=4, label="Eq. 6 (FS mín.)")
        ax.legend(frameon=False, fontsize=7.5, labelcolor=TINTA, loc="upper right")
    fig.suptitle(f"{solo}: as três variantes do transiente", color=TINTA, x=0.01, ha="left")
    fig.tight_layout()
    salva(fig, f"variantes_{solo}")

# ---- F3: diagnostico do transiente por solo (linhas: infiltracao, nos alagados, p_a max; colunas: intensidades)
for solo, (ks, Is) in SOLOS.items():
    nl = 3 if "i" in VARS else 2  # a linha de p_a so com a variante i
    fig, axs = plt.subplots(nl, 3, figsize=(15, 3 * nl), sharex=True)
    for j, I in enumerate(Is):
        for v in VARS:
            g = DIAG[(solo, I, v)]
            if not g:
                continue
            t = [r["t_h"] for r in g]
            kw = dict(color=COR_VAR[v], ls=ESTILO_VAR[v], lw=2, label=NOME_VAR[v])
            axs[0, j].plot(t, [r["infiltrado_m3m"] for r in g], **kw)
            axs[1, j].plot(t, [r["nos_face"] for r in g], **kw)
            if v == "i":
                axs[2, j].plot(t, [r["pa_max"] for r in g], **kw)
        chuva = [x * float(I) * 1e-3 * 74.142 for x in range(0, 25)] + [24 * float(I) * 1e-3 * 74.142] * 24
        axs[0, j].plot(range(49), chuva, color=TINTA2, lw=1, ls=(0, (1, 2)), label="chuva aplicada")
        for i, yl in enumerate(["água infiltrada líquida (m³/m)", "nós da face alagados", "$p_a$ máx. (kPa), variante (i)"][:nl]):
            estilo(axs[i, j], f"{I} mm/h" if i == 0 else "", yl if j == 0 else "")
            if i < nl - 1:
                axs[i, j].set_xlabel("")
        axs[0, j].legend(frameon=False, fontsize=7, labelcolor=TINTA, loc="upper left")
        axs[1, j].set_ylim(bottom=-0.5, top=max(2, axs[1, j].get_ylim()[1]))
        axs[1, j].yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
    fig.suptitle(f"{solo}: o transiente", color=TINTA, x=0.01, ha="left")
    fig.tight_layout()
    salva(fig, f"diag_{solo}")

# ---- F4: campos de p_w (e p_a na i) no f50 com 36 mm/h em 24 h
cmap = LinearSegmentedColormap.from_list("div", ["#a32d2c", "#e34948", "#f3b8a9", "#f0efec", "#9ec5f4", "#2a78d6",
                                                 "#104281"])
caminho = Path(GEO + [GEO[0]])
for (solo, I, t) in [("f50", "36", 24.0), ("f50", "90", 24.0), ("f100", "9", 24.0)]:
    cs = [(v, le_pw(solo, I, v, t)) for v in VARS]
    cs = [(v, c) for v, c in cs if c]
    if not cs:
        continue
    npan = len(cs) + (1 if any(v == "i" for v, _ in cs) else 0)
    fig, axs = plt.subplots(1, npan, figsize=(4.6 * npan, 3.0))
    axs = list(axs) if npan > 1 else [axs]
    k = 0
    for v, (xy, pw, pa) in cs:
        x, y = zip(*xy)
        tri = mtri.Triangulation(x, y)
        cx = [sum(x[a] for a in tr) / 3 for tr in tri.triangles]
        cy = [sum(y[a] for a in tr) / 3 for tr in tri.triangles]
        tri.set_mask([not caminho.contains_point((a, b), radius=-1e-6) for a, b in zip(cx, cy)])
        campos = [(f"$p_w$, {NOME_VAR[v]}", pw)] + ([("$p_a$, (i) acoplado, ar ativo", pa)] if pa else [])
        for nome, val in campos:
            ax = axs[k]
            k += 1
            if nome.startswith("$p_a"):
                c = ax.tricontourf(tri, val, levels=12, cmap=LinearSegmentedColormap.from_list(
                    "seq", ["#f0efec", "#9ec5f4", "#2a78d6", "#0d366b"]))
            else:
                c = ax.tricontourf(tri, [min(max(p, -75), 225) for p in val], levels=list(range(-75, 226, 15)),
                                   cmap=cmap, norm=TwoSlopeNorm(0, -75, 225))
                ax.tricontour(tri, val, levels=[0], colors=["#e87ba4"], linewidths=1.5)
            ax.plot(*zip(*(GEO + [GEO[0]])), color=TINTA2, lw=0.8)
            ax.set_aspect("equal")
            ax.set_title(nome, fontsize=8.5, color=TINTA, loc="left")
            ax.tick_params(colors=TINTA2, labelsize=7)
            for s in ax.spines.values():
                s.set_color(GRADE)
            fig.colorbar(c, ax=ax, shrink=0.75, label="kPa")
    fig.suptitle(f"{solo}, {I} mm/h, t = {t:g} h (linha rosa: $p_w$ = 0)", color=TINTA, x=0.01, ha="left", fontsize=10)
    fig.tight_layout()
    salva(fig, f"campos_{solo}_I{I}_t{t:g}")


# ---- formatos das tabelas
def f2(x):
    return "--" if x is None else f"{x:.2f}".replace(".", "{,}")


def f3(x):
    return "--" if x is None else f"{x:.3f}".replace(".", "{,}")


# ---- F5: as malhas (rodadas curtas com malhas=, ~/projects/runs/det_pc_3solos/malhas/<solo>/malha_*.csv)
def le_malha(arq):
    els = []
    for l in open(arq):
        if l.startswith("#"):
            continue
        v = [float(x) for x in l.split(",")]
        n = int(v[0])
        els.append([(v[1 + 2 * k], v[2 + 2 * k]) for k in range(n)])
    return els


def lado(e):
    return max(math.dist(e[k], e[(k + 1) % len(e)]) for k in range(len(e)))


MALHAS = os.path.join(RUNS, "malhas")
INFO_MALHA = {}
if os.path.isdir(MALHAS):
    from matplotlib.collections import PolyCollection
    paineis = [("transiente (tri2), todos os instantes", "f50/malha_transiente.csv"),
               ("fator, t = 0 h (tri1 + 4 níveis)", "f50/malha_fator_t0.csv"),
               ("fator, f10 180 mm/h, t = 24 h", "f10/malha_fator_t24.csv"),
               ("fator, f50 36 mm/h, t = 24 h", "f50/malha_fator_t24.csv"),
               ("fator, f100 9 mm/h, t = 24 h", "f100/malha_fator_t24.csv")]
    fig, axs = plt.subplots(2, 3, figsize=(15, 6.2))
    for ax, (tit, arq) in zip(axs.flat, paineis):
        caminho_arq = os.path.join(MALHAS, arq)
        if not os.path.exists(caminho_arq):
            ax.axis("off")
            continue
        els = le_malha(caminho_arq)
        INFO_MALHA[arq] = (len(els), min(map(lado, els)), max(map(lado, els)))
        ax.add_collection(PolyCollection(els, facecolors=SUPERF, edgecolors="#256abf", linewidths=0.35))
        ax.set_xlim(0, 70)
        ax.set_ylim(0, 30)
        ax.set_aspect("equal")
        ax.set_title(f"{tit}: {len(els)} triângulos", fontsize=8.5, color=TINTA, loc="left")
        ax.tick_params(colors=TINTA2, labelsize=7)
        for sp in ax.spines.values():
            sp.set_color(GRADE)
    axs.flat[-1].axis("off")
    fig.tight_layout()
    salva(fig, "malhas")


# ---- custo e convergencia do refinamento (variante iii; com as outras, por variante)
def le_csv(solo, I, v):
    arq = os.path.join(pasta(solo, I, v), "det_0_-1.csv")
    return list(csv.DictReader(open(arq))) if os.path.exists(arq) else None


L, LN = [], []
for solo, (ks, Is) in SOLOS.items():
    for I in Is:
        for v in VARS:
            rows, g, tm = le_csv(solo, I, v), DIAG[(solo, I, v)], le_tempo(solo, I, v)
            if not rows or not g or not tm:
                continue
            eqf = [int(r["eq_fator"]) for r in rows]
            sf = [float(r["s_fator"]) for r in rows]
            L.append(f"{solo} & {I.replace('.', '{,}')} & ({v}) & {rows[0]['eq_transiente']} & {int(g[-1]['passos'])} & "
                     f"{int(g[-1]['newton'])} & {g[-1]['newton'] / g[-1]['passos']:.1f} & {tm[1]:.0f} & "
                     f"{min(eqf)}--{max(eqf)} & {sum(sf) / len(sf):.0f} & {min(sf):.0f}--{max(sf):.0f} & "
                     f"{tm[0]:.0f}\\\\".replace(".", "{,}"))
            r24 = next((r for r in rows if float(r["t_h"]) == 24.0), None)
            if r24:
                nv_ = [float(x) for x in r24["niveis"].split(";")]
                LN.append(f"{solo} & {I.replace('.', '{,}')} & ({v}) & " + " & ".join(f3(x) for x in nv_) +
                          f" & {f3(nv_[-2] - nv_[-1])}\\\\")
    L.append("\\midrule")
    LN.append("\\midrule")
cab = ("\\begin{tabular}{llc ccccc ccc c}\n\\toprule\n & & & \\multicolumn{5}{c}{transiente} & "
       "\\multicolumn{3}{c}{fator (49 instantes)} & \\\\\n\\cmidrule(lr){4-8}\\cmidrule(lr){9-11}\n"
       "solo & $I_r$ & var. & equações & passos & Newton & Newton/passo & s & equações & s/FS & s/FS (mín--máx) & "
       "s total\\\\\n\\midrule\n")
open(os.path.join(TABS, "custo.tex"), "w").write(cab + "\n".join(L[:-1]) + "\n\\bottomrule\n\\end{tabular}\n")
cab = ("\\begin{tabular}{llc ccccc c}\n\\toprule\nsolo & $I_r$ & var. & tri1 & +1 & +2 & +3 & +4 & "
       "(+3) $-$ (+4)\\\\\n\\midrule\n")
open(os.path.join(TABS, "niveis.tex"), "w").write(cab + "\n".join(LN[:-1]) + "\n\\bottomrule\n\\end{tabular}\n")

# ---- tabelas


L = []
for solo, (ks, Is) in SOLOS.items():
    for I in Is:
        cel = []
        for v in VARS:
            d = DADOS[(solo, I, v)]
            if d:
                tm, fm = minimo(d)
                cel.append(f"{f3(fm)} ({tm:.0f})")
            else:
                cel.append("--")
        o = OUT[solo][I]
        L.append(f"{solo} & {I.replace('.', '{,}')} & {f2(EQ6[solo][I])} & " + " & ".join(cel) +
                 f" & {f2(o[0])} & {f2(o[1])} & {f2(o[2])}\\\\")
    L.append("\\midrule")
nv = len(VARS)
cab = ("\\begin{tabular}{llc" + "c" * nv + "ccc}\n\\toprule\n & & & \\multicolumn{" + str(nv) +
       "}{c}{este trabalho} & \\multicolumn{3}{c}{relatório de outubro}\\\\\n\\cmidrule(lr){4-" + str(3 + nv) +
       "}\\cmidrule(lr){" + str(4 + nv) + "-" + str(6 + nv) + "}\nsolo & $I_r$ (mm/h) & Eq.~6 & " +
       " & ".join(f"({v})" for v in VARS) + " & (iii) & (ii) & (i)\\\\\n\\midrule\n")
open(os.path.join(TABS, "fsmin.tex"), "w").write(cab + "\n".join(L[:-1]) + "\n\\bottomrule\n\\end{tabular}\n")

L = []
for I in ["9", "36", "90"]:
    for tt in (6, 12, 24, 48):
        cel = []
        for v in VARS:
            d = DADOS[("f50", I, v)]
            cel.append(f3(float(d[float(tt)]["FS_min"])) if d and float(tt) in d else "--")
        o = OUT5B[I][tt]
        L.append(f"{I} & {tt} & {f2(FIG5B[I][tt])} & " + " & ".join(cel) + f" & {f2(o[0])} & {f2(o[1])}\\\\")
    L.append("\\midrule")
cab = ("\\begin{tabular}{ccc" + "c" * nv + "cc}\n\\toprule\n & & & \\multicolumn{" + str(nv) +
       "}{c}{este trabalho} & \\multicolumn{2}{c}{outubro}\\\\\n\\cmidrule(lr){4-" + str(3 + nv) +
       "}\\cmidrule(lr){" + str(4 + nv) + "-" + str(5 + nv) + "}\n$I_r$ (mm/h) & $t$ (h) & Fig.~5b & " +
       " & ".join(f"({v})" for v in VARS) + " & (iii) & (i)\\\\\n\\midrule\n")
open(os.path.join(TABS, "fig5b.tex"), "w").write(cab + "\n".join(L[:-1]) + "\n\\bottomrule\n\\end{tabular}\n")

L = []
for solo, (ks, Is) in SOLOS.items():
    for I in Is:
        for v in VARS:
            d, g, tm = DADOS[(solo, I, v)], DIAG[(solo, I, v)], le_tempo(solo, I, v)
            if not d or not g:
                continue
            g24 = next((r for r in g if abs(r["t_h"] - 24) < 1e-9), None)
            gmax_pa = max(r["pa_max"] for r in g)
            fs0 = float(d[0.0]["FS_min"]) if 0.0 in d else None
            fs24 = float(d[24.0]["FS_min"]) if 24.0 in d else None
            fs48 = float(d[48.0]["FS_min"]) if 48.0 in d else None
            tmn, _ = minimo(d)
            prof = float(d[tmn]["banda_prof"])
            L.append(f"{solo} & {I.replace('.', '{,}')} & ({v}) & {f3(fs0)} & {f3(fs24)} & {f3(fs48)} & "
                     f"{f2(prof)} & {f2(g24['infiltrado_m3m']) if g24 else '--'} & "
                     f"{int(g24['nos_face']) if g24 else '--'} & {f2(gmax_pa)} & "
                     f"{f3(max(r['u_max_m'] for r in g))} & {int(g[-1]['passos'])} & {int(g[-1]['newton'])} & "
                     f"{tm[1]:.0f} & {tm[0]:.0f}\\\\")
        if nv > 1:
            L.append("\\midrule")
    if nv == 1:
        L.append("\\midrule")
cab = ("\\begin{tabular}{llc ccc c cc c c cc cc}\n\\toprule\nsolo & $I_r$ & var. & FS 0 & FS 24 & FS 48 & banda & "
       "infiltr. & alag. & $p_a$ máx & $|u|$ máx & passos & Newton & s trans. & s total\\\\\n\\midrule\n")
open(os.path.join(TABS, "diag.tex"), "w").write(cab + "\n".join(L[:-1]) + "\n\\bottomrule\n\\end{tabular}\n")

# ---- primeira passagem
L = []
for solo, (ks, Is) in SOLOS.items():
    for I in Is:
        cel = []
        for v in VARS:
            d = DADOS[(solo, I, v)]
            if not d:
                cel += ["--"] * 3
                continue
            t1, t2 = cruza(d)
            if t1 is None:
                cel += ["não cruza", "", ""]
            else:
                cel += [f"{t1:.1f}".replace(".", "{,}"), "--" if t2 is None else f"{t2:.1f}".replace(".", "{,}"),
                        "--" if t2 is None else f"{t2 - t1:.1f}".replace(".", "{,}")]
        L.append(f"{solo} & {I.replace('.', '{,}')} & " + " & ".join(cel) + "\\\\")
    L.append("\\midrule")
cab = ("\\begin{tabular}{ll" + "ccc" * nv + "}\n\\toprule\n & & " +
       " & ".join(f"\\multicolumn{{3}}{{c}}{{({v})}}" for v in VARS) + "\\\\\n" +
       "".join(f"\\cmidrule(lr){{{3 + 3 * k}-{5 + 3 * k}}}" for k in range(nv)) + "\n" +
       "solo & $I_r$ (mm/h) & " + " & ".join("$t_1$ (h) & volta (h) & abaixo (h)" for _ in VARS) +
       "\\\\\n\\midrule\n")
open(os.path.join(TABS, "cruza.tex"), "w").write(cab + "\n".join(L[:-1]) + "\n\\bottomrule\n\\end{tabular}\n")

# ---- numeros citados no texto
N = {}
falt = [f"{s}_I{I}_{v}" for (s, I, v), d in DADOS.items() if not d]
N["nVar"] = str(len(VARS))
N["nCasos"] = str(9 * len(VARS))
N["nFeitos"] = str(sum(1 for d in DADOS.values() if d))
N["faltam"] = ", ".join(falt).replace("_", "\\_") if falt else "nenhum"
for solo, (ks, Is) in SOLOS.items():
    for I in Is:
        for v in VARS:
            d = DADOS[(solo, I, v)]
            if d:
                tm, fm = minimo(d)
                chave = f"{solo}I{I.replace('.', 'p')}{v}"
                chave = re.sub(r"\d", lambda m: "abcdefghij"[int(m.group())], chave)
                N["fs" + chave] = f3(fm)
                N["t" + chave] = f"{tm:.0f}"
N["nTriTrans"] = str(INFO_MALHA.get("f50/malha_transiente.csv", ("--",))[0])
N["nTriFatorZero"] = str(INFO_MALHA.get("f50/malha_fator_t0.csv", ("--",))[0])
with open(os.path.join(AQUI, "numeros.tex"), "w") as f:
    for k, x in N.items():
        f.write(f"\\newcommand{{\\{k}}}{{{x}}}\n")
print(f"figuras em {FIGS}, tabelas em {TABS}; casos com dados: {N['nFeitos']} de 27 (faltam: {N['faltam']})")
for solo, (ks, Is) in SOLOS.items():
    for I in Is:
        s = []
        for v in VARS:
            d = DADOS[(solo, I, v)]
            s.append("   --   " if not d else "{:.3f}@{:>2.0f}h t1={}".format(*minimo(d)[::-1], cruza(d)))
        print(f"{solo:5s} {I:>4s}  " + "  ".join(s) + f"   Eq6 {EQ6[solo][I]:.2f}  out {OUT[solo][I]}")
