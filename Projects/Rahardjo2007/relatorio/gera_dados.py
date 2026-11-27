#!/usr/bin/env python3
"""gera_dados.py: as figuras e tabelas do relatorio do Rahardjo2007 (relatorio_rahardjo2007.tex), com os dados que ja
existem (as rodadas que faltam ficam de fora, sem erro):
  ~/projects/runs/rahardjo2007/fermi/<id>/rah_<solo>_<I>_fs.csv e _nivel.csv  as 28 rodadas do fermi (cluster/casos.txt)
  ~/projects/runs/rahardjo2007/pc_malha/<id>_fina/...                          a conferencia de malha (casos_malha.txt)
  ~/projects/runs/rahardjo2007/malhas/m_transiente.csv, m_fator.csv            as malhas (Rahardjo2007 grava_malha=m)
  ../data/fig5_lido.csv, fig6_lido.csv, fig7_eq6.csv, tabela5.csv               o artigo (lido das figuras e tabelas)
Saidas: figs/malhas.pdf, figs/fs_tempo.pdf (Fig. 5), figs/fs_min.pdf (Figs. 6 e 7), figs/nivel.pdf (o nivel d'agua no
FS minimo), tab/fs_min.tex, tab/eq6.tex, tab/malha.tex, tab/niveis.tex e numeros.tex (as macros do texto)."""
import csv
import math
import os

import numpy as np
np.seterr(over="ignore")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.collections import PolyCollection  # noqa: E402
from scipy.optimize import curve_fit  # noqa: E402

AQUI = os.path.dirname(os.path.abspath(__file__))
R = os.environ.get("RAH_RUNS", os.path.expanduser("~/projects/runs/rahardjo2007"))
FERMI, MALHA, MALHAS = os.path.join(R, "fermi"), os.path.join(R, "pc_malha"), os.path.join(R, "malhas")
DATA = os.path.join(AQUI, "..", "data")
TINTA, TINTA2, GRADE = "#0b0b0b", "#5c5c58", "#e4e3df"
SOLOS = ["f10", "f50", "f100"]
NOMES = {"f10": "$f_{10,-4}$", "f50": "$f_{50,-5}$", "f100": "$f_{100,-6}$"}
KS = {"f10": 360., "f50": 36., "f100": 3.6}  # mm/h
TART = [0, 1, 2, 3, 6, 12, 24, 28, 32, 40, 48]  # os instantes do FS no artigo
os.makedirs(os.path.join(AQUI, "figs"), exist_ok=True)
os.makedirs(os.path.join(AQUI, "tab"), exist_ok=True)
v2 = lambda x: f"{x:.2f}".replace(".", "{,}")
v3 = lambda x: f"{x:.3f}".replace(".", "{,}")
vI = lambda x: f"{x:g}".replace(".", "{,}")
N = {}  # as macros do numeros.tex


def le(arq):
    return list(csv.DictReader(open(arq))) if os.path.exists(arq) else []


def comentado(arq):
    L = [r for r in csv.reader(open(arq)) if r and not r[0].startswith("#")]
    return [dict(zip(L[0], r)) for r in L[1:]]


def eixo(ax):
    ax.grid(True, color=GRADE, lw=0.7)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.tick_params(labelsize=8, colors=TINTA2)


# ---- o artigo
F5 = {}
for r in comentado(os.path.join(DATA, "fig5_lido.csv")):
    F5.setdefault((r["solo"], float(r["I_mmh"])), []).append((float(r["t_h"]), float(r["FS"])))
F6 = {(r["solo"], float(r["I_mmh"])): float(r["FS_min"]) for r in comentado(os.path.join(DATA, "fig6_lido.csv"))}
SB = {(r["solo"], float(r["I_mmh"])): r["serie_b"] == "1" for r in comentado(os.path.join(DATA, "fig6_lido.csv"))}
EQ6 = {r["solo"]: r for r in comentado(os.path.join(DATA, "fig7_eq6.csv"))}
T5 = {(r["solo"], float(r["I_mmh"])): r["categoria"] for r in comentado(os.path.join(DATA, "tabela5.csv"))}
CASOS = sorted(F6, key=lambda k: (SOLOS.index(k[0]), k[1]))  # os 28 (solo, I) das Figs. 6 e 7


# ---- as nossas rodadas: o FS de hora em hora, o nivel d'agua e a categoria
def superficie(x):
    return 30. if x <= 30. else 20. if x >= 40. else 30. - (x - 30.)


def cat_zona(pts, tol=0.05):
    """a categoria da Tabela 3 pelo ponto mais alto da superficie que o nivel d'agua toca (a menos de tol): a crista
    (x <= 30 m) da I, a metade de cima da face (x <= 35 m) II, o resto da face e o pe III; nao tocar da IV. O programa
    (Categoria em main.cpp) testa so tres pontos (x = 30, 35 e 40 m): na quina da crista o nivel desce ~1 m (drena para
    a face) mesmo com a crista ja alagada, e o f10 a 180 mm/h sai II em vez de I"""
    toca = [x for x, y in pts if y >= superficie(x) - tol]
    if not toca:
        return "IV"
    xm = min(toca)
    return "I" if xm <= 30. else "II" if xm <= 35. else "III"


def rodada(pasta, solo, I):
    base = os.path.join(pasta, f"rah_{solo}_{I:g}")
    L = le(base + "_fs.csv")
    if not L:
        return None
    fim = os.path.exists(os.path.join(pasta, "fim")) and open(os.path.join(pasta, "codigo")).read().strip() == "0"
    t = np.array([float(r["t_h"]) for r in L])
    fs = np.array([float(r["FS"]) for r in L])
    niv = [[float(x) for x in r["niveis"].split(";")] if r["niveis"] else [] for r in L]
    nivel = {}
    for r in le(base + "_nivel.csv"):
        nivel.setdefault(float(r["t_h"]), []).append((float(r["x"]), float(r["y_nivel"])))
    cat = [cat_zona(nivel[float(r["t_h"])]) if float(r["t_h"]) in nivel else r["categoria"] for r in L]
    return {"t": t, "fs": fs, "cat": cat, "niveis": niv, "nivel": nivel, "fim": fim,
            "s": float(L[-1]["s_total"]), "pw_max": [float(r["pw_max"]) for r in L],
            "face": [int(r["nos_face"]) for r in L], "infil": [float(r["infiltrado_m3m"]) for r in L],
            "passos": [int(r["passos"]) for r in L], "s_trans": [float(r["s_transiente"]) for r in L],
            "s_fator": [float(r["s_fator"]) for r in L],
            "nc": [float(r["nivel_crista"]) for r in L], "nm": [float(r["nivel_meio"]) for r in L],
            "np": [float(r["nivel_pe"]) for r in L]}


NOS = {}
for (s, I) in CASOS:
    d = rodada(os.path.join(FERMI, f"{s}_I{I:g}"), s, I)
    if d is not None:
        NOS[(s, I)] = d
FINA = {}
for s, I in (("f10", 900.), ("f50", 36.), ("f100", 36.)):
    d = rodada(os.path.join(MALHA, f"{s}_I{I:g}_fina"), s, I)
    if d is not None:
        FINA[(s, I)] = d


def minimo(d, so_artigo=False):
    """o FS minimo e o instante (h); so_artigo: so nos instantes do artigo"""
    ok = [i for i, t in enumerate(d["t"]) if d["fs"][i] > 0 and (not so_artigo or t in TART)]
    if not ok:
        return None, None, None
    i = min(ok, key=lambda j: d["fs"][j])
    return d["fs"][i], d["t"][i], d["cat"][i]


# ---- Eq. 6 e as intensidades de inicio e de limiar
def sig(I, A, B, C, D):
    return A + B / (1. + np.exp(-(I - C) / D))


def ini_thr(A, B, C, D):
    """I_ini e I_thr como na Fig. 6 do artigo: a tangente ao trecho inclinado da curva FS x log10(I) (no ponto de maior
    inclinacao) cruzando os dois patamares (o de cima: o FS da curva em I -> 0; o de baixo: A)"""
    x = np.linspace(-2., 4., 60001)
    f = sig(10. ** x, A, B, C, D)
    df = np.gradient(f, x)
    k = int(np.argmin(df))
    topo = sig(1.e-6, A, B, C, D)
    xi = x[k] + (topo - f[k]) / df[k]
    xt = x[k] + (A - f[k]) / df[k]
    return 10. ** xi, 10. ** xt


def ajusta(Is, F, p0, Cfixo=None):
    """a Eq. 6 ajustada por minimos quadrados. Cfixo: C fica fixo (o f10, onde o C do artigo e negativo: para I >= 0 a
    Eq. 6 vira ~A + B exp(C/|D|) exp(-I/|D|), so o produto B exp(C/|D|) e identificavel, e o ajuste livre vai a B ~ 1e4,
    C ~ -500 com a mesma curva, tanto nos nossos pontos quanto nos lidos do artigo)"""
    Is, F = np.array(Is), np.array(F)
    try:
        if Cfixo is None:
            p, _ = curve_fit(sig, Is, F, p0=p0, maxfev=200000)
        else:
            q, _ = curve_fit(lambda I, A, B, D: sig(I, A, B, Cfixo, D), Is, F, p0=[p0[0], p0[1], p0[3]], maxfev=200000)
            p = np.array([q[0], q[1], Cfixo, q[2]])
        res = F - sig(Is, *p)
        r2 = 1. - (res ** 2).sum() / ((F - F.mean()) ** 2).sum()
        return p, r2
    except Exception:
        return None, None


# ---- figura 1: as malhas e o nivel inicial
def poligonos(arq):
    return [[(float(r["x1"]), float(r["y1"])), (float(r["x2"]), float(r["y2"])), (float(r["x3"]), float(r["y3"]))]
            for r in le(arq)]


PT, PF = poligonos(os.path.join(MALHAS, "m_transiente.csv")), poligonos(os.path.join(MALHAS, "m_fator.csv"))
if PT and PF:
    from matplotlib.patches import Rectangle  # noqa: E402
    import matplotlib.patheffects as pe  # noqa: E402
    fig = plt.figure(figsize=(7.4, 6.4))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.45, 1.], width_ratios=[1., 1.75])
    axA, axB, axC = fig.add_subplot(gs[0, :]), fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])
    xs = np.linspace(0, 70, 141)
    nivel0 = 20. - 5. + (40. - xs) * math.tan(math.radians(7.))
    zoom = (27., 33., 25.5, 30.5)  # o detalhe da crista (x0, x1, y0, y1)
    paineis = (
        (axA, PT, (-1, 71), (-1, 31), f"(a) a malha do transiente, a de todas as rodadas: {len(PT)} triângulos, "
                                      "colunas a cada 1 m, camadas que seguem o terreno"),
        (axB, PT, zoom[:2], zoom[2:], "(b) detalhe da crista: a camada de cima\ntem 9 cm e cresce 15 % até 0,75 m"),
        (axC, PF, (-1, 71), (-1, 31), f"(c) a malha inicial do fator: {len(PF)} triângulos,\nrefinada 4 vezes na banda de "
                                      "ruptura em cada instante (Fig. 6)"))
    for ax, P, xl, yl, tit in paineis:
        ax.add_collection(PolyCollection(P, facecolor="none", edgecolor=TINTA2,
                                         lw=(0.12 if ax is axA else 0.3) if len(P) > 1000 else 0.6))
        ax.plot(xs, nivel0, color="#2a78d6", lw=1.6, ls=(0, (5, 2)), label="nível d'água inicial: 5 m abaixo do pé, subindo 7°",
                path_effects=[pe.Stroke(linewidth=3.6, foreground="white"), pe.Normal()])
        ax.set_xlim(*xl), ax.set_ylim(*yl), ax.set_aspect("equal")
        ax.set_title(tit, fontsize=8, loc="left", color=TINTA)
        ax.tick_params(labelsize=7.5, colors=TINTA2)
    axA.add_patch(Rectangle((zoom[0], zoom[2]), zoom[1] - zoom[0], zoom[3] - zoom[2], fill=False, ec="#e34948", lw=1.2))
    axA.text(zoom[0] - 0.5, zoom[2] - 0.5, "(b)", color="#e34948", fontsize=8, ha="right", va="top",
             bbox=dict(facecolor="white", edgecolor="none", pad=0.6))
    axA.legend(fontsize=8, loc="lower left", frameon=True, framealpha=1., facecolor="white", edgecolor="none")
    for ax in (axB, axC):
        ax.set_xlabel("x (m)", fontsize=8, color=TINTA)
    fig.tight_layout()
    fig.savefig(os.path.join(AQUI, "figs", "malhas.pdf"))
    plt.close(fig)

# ---- figura 2: o FS contra o tempo (Fig. 5)
fig, axs = plt.subplots(3, 1, figsize=(7.4, 10.4), sharex=True)
for ax, s in zip(axs, SOLOS):
    eixo(ax)
    Is = sorted(I for (so, I) in F5 if so == s)
    cores = plt.cm.viridis(np.linspace(0.0, 0.88, len(Is)))
    for I, c in zip(Is, cores):
        a = sorted(F5[(s, I)])
        ax.plot([t for t, _ in a], [f for _, f in a], "o", ms=4.5, mfc="none", mec=c, mew=1.1, zorder=4)
        d = NOS.get((s, I))
        if d is not None:
            ok = d["fs"] > 0
            ax.plot(d["t"][ok], d["fs"][ok], color=c, lw=1.5, zorder=3)
        rot = f"{vI(I).replace('{,}', ',')} mm/h ({I / KS[s]:g} $k_s$)".replace(".", ",")
        ax.plot([], [], color=c, lw=1.5, marker="o", mfc="none", ms=4.5, label=rot)
    ax.axvline(24, color=TINTA2, lw=0.8, ls=(0, (4, 3)))
    ax.set_ylim(0.6, 2.25), ax.set_xlim(0, 48.5)
    ax.set_xticks(range(0, 49, 6))
    ax.set_ylabel("FS", fontsize=9, color=TINTA)
    ax.set_title(f"{NOMES[s]}, $k_s$ = {KS[s]:g} mm/h".replace(".", ","), fontsize=9, loc="left", color=TINTA)
    ax.legend(fontsize=7, frameon=False, loc="lower left", ncol=2)
axs[-1].set_xlabel("t (h)  (chuva até 24 h)", fontsize=9, color=TINTA)
fig.text(0.99, 0.005, "linhas: NeoPZ, modelo (iii), FS de hora em hora; círculos: Rahardjo et al. (2007), Fig. 5",
         fontsize=7, color=TINTA2, ha="right")
fig.tight_layout(rect=(0, 0.01, 1, 1))
fig.savefig(os.path.join(AQUI, "figs", "fs_tempo.pdf"))
plt.close(fig)

# ---- figura 3: o FS minimo contra I (Figs. 6 e 7), com a Eq. 6 do artigo e a nossa
COR = {"f10": "#e34948", "f50": "#2a78d6", "f100": "#3a9a5b"}
MARC = {"f10": "o", "f50": "s", "f100": "^"}
AJ = {}
fig, axs = plt.subplots(1, 3, figsize=(7.6, 3.4), sharey=True)
for ax, s in zip(axs, SOLOS):
    eixo(ax)
    ax.set_xscale("log")
    e = EQ6[s]
    pa = [float(e[k]) for k in "ABCD"]
    x = np.logspace(-1, 4, 400)
    ax.plot(x, sig(x, *pa), color=TINTA2, lw=1.0, ls=(0, (5, 2)), label="Eq. 6 do artigo")
    Ia = sorted(I for (so, I) in F6 if so == s)
    Io = [I for I in Ia if (s, I) in NOS and NOS[(s, I)]["fim"]]
    if Io:
        Fo = [minimo(NOS[(s, I)])[0] for I in Io]
        ax.plot(Io, Fo, MARC[s], ms=4.5, color=COR[s], label="NeoPZ (iii)", zorder=4)
        if len(Io) >= 5:
            p, r2 = ajusta(Io, Fo, pa, pa[2] if pa[2] < 0. else None)
            if p is not None:
                AJ[s] = (p, r2, len(Io), pa[2] < 0.)
                ax.plot(x, sig(x, *p), color=COR[s], lw=1.4, label="Eq. 6 ajustada", zorder=3)
    # o artigo por cima, em anel maior: o nosso ponto aparece dentro do anel quando bate
    ax.plot(Ia, [F6[(s, I)] for I in Ia], MARC[s], ms=8.5, mfc="none", mec=TINTA, mew=0.9, label="artigo (Fig. 7)",
            zorder=5)
    ax.axvline(KS[s], color=GRADE, lw=1.2, zorder=0)
    ax.text(KS[s] * 1.15, 0.15, "$k_s$", fontsize=8, color=TINTA2)
    ax.set_xlim(0.1, 1.e4), ax.set_ylim(0., 2.5)
    ax.set_title(NOMES[s], fontsize=9, loc="left", color=TINTA)
    ax.set_xlabel("$I_r$ (mm/h)", fontsize=9, color=TINTA)
axs[0].set_ylabel("FS mínimo (0 a 48 h)", fontsize=9, color=TINTA)
axs[0].legend(fontsize=7, frameon=False, loc="lower left")
fig.tight_layout()
fig.savefig(os.path.join(AQUI, "figs", "fs_min.pdf"))
plt.close(fig)

# ---- figura 4: o nivel d'agua no instante do FS minimo (as categorias da Tabela 3)
fig, axs = plt.subplots(3, 1, figsize=(7.4, 8.4), sharex=True)
xsup = [0, 30, 40, 70]
ysup = [30, 30, 20, 20]
xs = np.linspace(0, 70, 141)
for ax, s in zip(axs, SOLOS):
    ax.plot(xsup, ysup, color=TINTA, lw=1.2)
    ax.plot(xs, 15. + (40. - xs) * math.tan(math.radians(7.)), color=TINTA2, lw=1.0, ls=(0, (5, 2)),
            label="nível inicial")
    Is = sorted(I for (so, I) in F5 if so == s)
    cores = plt.cm.viridis(np.linspace(0.0, 0.88, len(Is)))
    for I, c in zip(Is, cores):
        d = NOS.get((s, I))
        if d is None or not d["fim"]:
            continue
        fm, tm, cat = minimo(d)
        pts = d["nivel"].get(tm)
        if pts:
            ax.plot([p[0] for p in pts], [p[1] for p in pts], color=c, lw=1.4,
                    label=f"{vI(I).replace('{,}', ',')} mm/h: {tm:g} h, cat. {cat} (artigo {T5.get((s, I), '-')})")
    ax.set_ylim(8, 31), ax.set_xlim(0, 70)
    ax.set_aspect("equal")
    ax.tick_params(labelsize=8, colors=TINTA2)
    ax.set_title(f"{NOMES[s]}: nível d'água no FS mínimo", fontsize=9, loc="left", color=TINTA)
    ax.legend(fontsize=6.5, frameon=False, loc="lower left", ncol=2)
axs[-1].set_xlabel("x (m)", fontsize=9, color=TINTA)
fig.tight_layout()
fig.savefig(os.path.join(AQUI, "figs", "nivel.pdf"))
plt.close(fig)

# ---- tabela: o FS minimo, o instante e a categoria, nosso contra o artigo
dif, cats = [], []
with open(os.path.join(AQUI, "tab", "fs_min.tex"), "w") as f:
    f.write("\\begin{tabular}{lrr|ccc|ccc|cc}\n\\toprule\n")
    f.write("solo & $I_r$ & $I_r/k_s$ & \\multicolumn{3}{c|}{FS inicial e FS mín.} & \\multicolumn{3}{c|}{nosso: FS mín.} "
            "& \\multicolumn{2}{c}{categoria}\\\\\n")
    f.write(" & (mm/h) & & nosso, 0 h & artigo & nosso & instante (h) & nos instantes & dif. & nossa & artigo\\\\\n"
            " & & & & (Fig. 6) & (0--48 h) & & do artigo & & & (Tab. 5)\\\\\n\\midrule\n")
    for k, (s, I) in enumerate(CASOS):
        if k and s != CASOS[k - 1][0]:
            f.write("\\midrule\n")
        d = NOS.get((s, I))
        lin = f"{NOMES[s]} & {vI(I)} & {vI(round(I / KS[s], 4))} & "
        if d is None:
            f.write(lin + "\\multicolumn{8}{c}{(rodando)}\\\\\n")
            continue
        fm, tm, cat = minimo(d)
        fa, ta, _ = minimo(d, True)
        art = F6[(s, I)]
        ok = d["fim"]
        if ok and fm is not None:
            dif.append(fm - art)
            if (s, I) in T5:
                cats.append(cat == T5[(s, I)])
        sd = "0{,}00" if fm and abs(fm - art) < 0.005 else ("$+$" if fm - art > 0 else "$-$") + v2(abs(fm - art))
        f.write(lin + f"{v3(d['fs'][0])} & {v2(art)} & {v3(fm) if fm else '--'}{'' if ok else '$^*$'} & {tm:g} & "
                f"{v3(fa) if fa else '--'} & {sd if fm else '--'} & {cat or '--'} & {T5.get((s, I), '--')}\\\\\n")
    f.write("\\bottomrule\n\\end{tabular}\n")
if dif:
    N["difMedia"] = v2(float(np.mean(dif)))
    N["difMax"] = v2(float(np.max(np.abs(dif))))
    N["nCasos"] = str(len(dif))
    N["catOk"], N["catTot"] = str(sum(cats)), str(len(cats))

# ---- a Fig. 5 ponto a ponto: o nosso FS no instante de cada marcador do artigo
d5 = [dict(zip(NOS[k]["t"], NOS[k]["fs"]))[t] - f for k, a in F5.items() if k in NOS and NOS[k]["fim"] for t, f in a]
if d5:
    d5 = np.array(d5)
    N["cincoN"], N["cincoRms"], N["cincoMax"] = str(len(d5)), v3(float(np.sqrt(np.mean(d5 ** 2)))), v2(float(np.abs(d5).max()))
    N["cincoFrac"] = f"{100. * np.mean(np.abs(d5) <= 0.03):.0f}"

# ---- tabela: os parametros da Eq. 6, do artigo e os nossos
with open(os.path.join(AQUI, "tab", "eq6.tex"), "w") as f:
    f.write("\\begin{tabular}{ll|cccc|c|cc}\n\\toprule\n")
    f.write("solo & & $A$ & $B$ & $C$ & $D$ & $r^2$ & $I_{r(ini)}$ (mm/h) & $I_{r(thr)}$ (mm/h)\\\\\n\\midrule\n")
    for s in SOLOS:
        e = EQ6[s]
        pa = [float(e[k]) for k in "ABCD"]
        ii, it = ini_thr(*pa)
        f.write(f"{NOMES[s]} & artigo & {v2(pa[0])} & {v2(pa[1])} & {vI(pa[2])} & {vI(pa[3])} & {e['r2'].replace('.', '{,}')}"
                f" & {vI(float(e['I_ini']))} ({v2(ii)}) & {vI(float(e['I_thr']))} ({v2(it)})\\\\\n")
        if s in AJ:
            p, r2, n, fixo = AJ[s]
            ii, it = ini_thr(*p)
            f.write(f" & NeoPZ ({n} pontos) & {v2(p[0])} & {v2(p[1])} & {vI(p[2]) + '$^\\dagger$' if fixo else v2(p[2])} & "
                    f"{v2(p[3])} & {v3(r2)} & {v2(ii)} & {v2(it)}\\\\\n")
        else:
            f.write(" & NeoPZ & \\multicolumn{7}{c}{(rodando)}\\\\\n")
        f.write("\\midrule\n" if s != SOLOS[-1] else "")
    f.write("\\bottomrule\n\\end{tabular}\n")

# ---- tabela: a conferencia de malha (o FS nos instantes do artigo, malha base contra malha e passo pela metade)
with open(os.path.join(AQUI, "tab", "malha.tex"), "w") as f:
    f.write("\\begin{tabular}{l|" + "c" * len(TART) + "|c}\n\\toprule\n")
    f.write("caso & " + " & ".join(f"{t} h" for t in TART) + " & máx. $|\\Delta|$\\\\\n\\midrule\n")
    for (s, I), fina in FINA.items():
        base = NOS.get((s, I))
        lin_b = [dict(zip(base["t"], base["fs"])).get(t) if base else None for t in TART]
        lin_f = [dict(zip(fina["t"], fina["fs"])).get(t) for t in TART]
        nome = f"{NOMES[s]}, {vI(I)} mm/h"
        f.write(f"{nome}: base & " + " & ".join(v3(x) if x else "--" for x in lin_b) + " & \\\\\n")
        dm = [abs(a - b) for a, b in zip(lin_b, lin_f) if a and b]
        f.write(f"\\phantom{{{nome}:}} fina & " + " & ".join(v3(x) if x else "--" for x in lin_f) +
                f" & {v3(max(dm)) if dm else '--'}\\\\\n")
    f.write("\\bottomrule\n\\end{tabular}\n")

# ---- tabela: o FS em cada nivel do refinamento adaptativo, no instante do FS minimo
with open(os.path.join(AQUI, "tab", "niveis.tex"), "w") as f:
    f.write("\\begin{tabular}{lr|ccccc|c}\n\\toprule\n")
    f.write("solo & $I_r$ & nível 0 & 1 & 2 & 3 & 4 & 4 $-$ 3\\\\\n\\midrule\n")
    for (s, I) in CASOS:
        d = NOS.get((s, I))
        if d is None or not d["fim"] or I not in [I2 for (s2, I2) in F5 if s2 == s]:
            continue
        fm, tm, _ = minimo(d)
        k = list(d["t"]).index(tm)
        nv = d["niveis"][k]
        f.write(f"{NOMES[s]} & {vI(I)} & " + " & ".join(v3(x) for x in nv) +
                (f" & {v3(nv[-1] - nv[-2])}" if len(nv) > 1 else " & --") + "\\\\\n")
    f.write("\\bottomrule\n\\end{tabular}\n")


# ---- primeira passagem: t_1 (o FS cruza 1 descendo), a volta acima de 1 e o tempo abaixo (interpolacao linear)
def cruza(d):
    t, f = d["t"], d["fs"]
    t1 = volta = None
    for i in range(1, len(t)):
        if t1 is None and f[i - 1] >= 1. > f[i]:
            t1 = t[i - 1] + (t[i] - t[i - 1]) * (f[i - 1] - 1.) / (f[i - 1] - f[i])
        elif t1 is not None and volta is None and f[i - 1] < 1. <= f[i]:
            volta = t[i - 1] + (t[i] - t[i - 1]) * (1. - f[i - 1]) / (f[i] - f[i - 1])
    return t1, volta


LSUP = 30. + 10. * math.sqrt(2.) + 30.  # m de superficie (a chuva entra por metro inclinado)
with open(os.path.join(AQUI, "tab", "cruza.tex"), "w") as f:
    f.write("\\begin{tabular}{lr|ccc}\n\\toprule\nsolo & $I_r$ (mm/h) & $t_1$ (h) & volta (h) & abaixo de 1 (h)\\\\\n\\midrule\n")
    nada = True
    for (s, I) in CASOS:
        d = NOS.get((s, I))
        if d is None or not d["fim"]:
            continue
        t1, vt = cruza(d)
        if t1 is None:
            continue
        nada = False
        f.write(f"{NOMES[s]} & {vI(I)} & {v2(t1)} & {v2(vt) if vt else '$>48$'} & "
                f"{v2((vt or 48.) - t1)}\\\\\n")
    if nada:
        f.write("\\multicolumn{5}{c}{(nenhuma rodada pronta rompe)}\\\\\n")
    f.write("\\bottomrule\n\\end{tabular}\n")

# ---- o transiente: infiltracao, alagamento, nivel d'agua em 24 h e custo
with open(os.path.join(AQUI, "tab", "diag.tex"), "w") as f, open(os.path.join(AQUI, "tab", "custo.tex"), "w") as g:
    f.write("\\begin{tabular}{lr|ccc|cc|ccc}\n\\toprule\n"
            " & & \\multicolumn{3}{c|}{FS} & \\multicolumn{2}{c|}{chuva até 24 h} & \\multicolumn{3}{c}{nível d'água em 24 h (m)}\\\\\n"
            "solo & $I_r$ & 0 h & 24 h & 48 h & infiltrou & alaga em & crista & meio & pé\\\\\n"
            " & (mm/h) & & & & (\\%) & (h) & (30) & (25) & (20)\\\\\n\\midrule\n")
    g.write("\\begin{tabular}{lr|cc|ccc|c}\n\\toprule\n"
            " & & \\multicolumn{2}{c|}{transiente (48 h)} & \\multicolumn{3}{c|}{fator} & \\\\\n"
            "solo & $I_r$ & passos & s & análises & s/FS & s/FS (mín--máx) & s total\\\\\n\\midrule\n")
    for k, (s, I) in enumerate(CASOS):
        d = NOS.get((s, I))
        if d is None or not d["fim"]:
            continue
        T = list(d["t"])
        if 24. not in T:
            continue
        i24 = T.index(24.)
        chuva = I / 1000. * 24. * LSUP
        alaga = next((t for t, n in zip(T, d["face"]) if n > 0), None)
        fs = dict(zip(T, d["fs"]))
        f.write(f"{NOMES[s]} & {vI(I)} & {v3(fs[0.])} & {v3(fs.get(24., float('nan')))} & {v3(fs.get(48., float('nan')))} & "
                f"{100. * d['infil'][i24] / chuva:.0f} & {('%g' % alaga) if alaga is not None else '--'} & "
                f"{v2(d['nc'][i24])} & {v2(d['nm'][i24])} & {v2(d['np'][i24])}\\\\\n")
        sf = [x for x, fsv in zip(d["s_fator"], d["fs"]) if fsv > 0]
        g.write(f"{NOMES[s]} & {vI(I)} & {d['passos'][-1]} & {d['s_trans'][-1]:.0f} & {len(sf)} & {np.mean(sf):.0f} & "
                f"{min(sf):.0f}--{max(sf):.0f} & {d['s']:.0f}\\\\\n")
    f.write("\\bottomrule\n\\end{tabular}\n")
    g.write("\\bottomrule\n\\end{tabular}\n")

# ---- os campos de p_w e as malhas refinadas do fator (PC: pc_campos/<id>/rah_<solo>_<I>_campo_<t>h.csv e _fator_)
from matplotlib.colors import TwoSlopeNorm  # noqa: E402
CAMPOS = os.path.join(R, "pc_campos")
CC = [("f10", 90.), ("f10", 900.), ("f50", 36.), ("f50", 180.), ("f100", 9.), ("f100", 36.)]
TC = [6., 24., 48.]


def le_campo(s, I, t):
    base = os.path.join(CAMPOS, f"{s}_I{I:g}", f"rah_{s}_{I:g}")
    C, F = le(f"{base}_campo_{t:g}h.csv"), le(f"{base}_fator_{t:g}h.csv")
    if not C or not F:
        return None
    tri = np.array([[float(r[k]) for k in ("x1", "y1", "p1", "x2", "y2", "p2", "x3", "y3", "p3")] for r in C])
    fat = np.array([[float(r[k]) for k in ("x1", "y1", "x2", "y2", "x3", "y3", "sqrtJ2_eps_p")] for r in F])
    fs = {float(r["t_h"]): float(r["FS"]) for r in le(f"{base}_fs.csv")}
    return tri, fat, fs.get(t)


from matplotlib.tri import Triangulation  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402


def triangulacao(V):
    """V: (n, 3, 2), os vertices de cada triangulo -> (Triangulation com os nos unicos, T: os nos de cada vertice NA
    ORDEM DE V (para os valores por vertice), area). Na Triangulation os triangulos vao no sentido anti-horario; um no
    pendurado (refinamento nao conforme) so fica sem ligacao com o vizinho"""
    nos, inv = np.unique(np.round(V.reshape(-1, 2), 6), axis=0, return_inverse=True)
    T = np.asarray(inv).reshape(-1, 3)
    ar = (V[:, 1, 0] - V[:, 0, 0]) * (V[:, 2, 1] - V[:, 0, 1]) - (V[:, 2, 0] - V[:, 0, 0]) * (V[:, 1, 1] - V[:, 0, 1])
    Tc = T.copy()
    Tc[ar < 0] = Tc[ar < 0][:, [0, 2, 1]]
    return Triangulation(nos[:, 0], nos[:, 1], Tc), T, 0.5 * np.abs(ar)


def no_max(T, val, nn):
    """o maior valor dos triangulos em volta de cada no (um valor constante por triangulo). A media diluiria a banda,
    que tem um triangulo de espessura: com o maximo, a isolinha de um nivel passa logo fora dos triangulos acima dele"""
    v = np.zeros(nn)
    for k in range(3):
        np.maximum.at(v, T[:, k], val)
    return v


def indicador(fat):
    """o indicador plastico do fator nos nos, relativo ao maior triangulo do instante (o criterio de refinamento:
    os triangulos com >= 0,1 sao divididos)"""
    tr, T, _ = triangulacao(fat[:, :6].reshape(-1, 3, 2))
    return tr, no_max(T, fat[:, 6], len(tr.x)) / fat[:, 6].max()


NPW = [-80., -60., -40., -20.]  # isolinhas de p_w (kPa): succao, o nivel d'agua (0) e positivas
PPW = [50., 100., 150., 200.]
prontos = [(s, I) for (s, I) in CC if all(le_campo(s, I, t) is not None for t in TC)]
if prontos:
    norm = TwoSlopeNorm(vmin=-85., vcenter=0., vmax=225.)
    fig, axs = plt.subplots(len(prontos), len(TC), figsize=(7.6, 1.75 * len(prontos) + 0.7), squeeze=False)
    for i, (s, I) in enumerate(prontos):
        for j, t in enumerate(TC):
            ax = axs[i][j]
            tri, fat, fs = le_campo(s, I, t)
            tr, T, _ = triangulacao(tri[:, [0, 1, 3, 4, 6, 7]].reshape(-1, 3, 2))
            p = np.zeros(len(tr.x))
            p[T.ravel()] = tri[:, [2, 5, 8]].ravel()  # p_w e continua: cada no tem um valor so
            pc = ax.tripcolor(tr, p, cmap="RdBu_r", norm=norm, shading="gouraud", rasterized=True)
            for lv, lw, ls in ((NPW, 0.45, "-"), (PPW, 0.45, "-"), ([0.], 1.1, "-")):
                if p.min() < max(lv) and p.max() > min(lv):
                    cs = ax.tricontour(tr, p, levels=lv, colors=TINTA, linewidths=lw, linestyles=ls)
                    ax.clabel(cs, fmt=lambda v: f"{v:.0f}", fontsize=4.5, inline=True, inline_spacing=1)
            ax.plot([0, 30, 40, 70], [30, 30, 20, 20], color=TINTA, lw=0.6)
            ax.set_xlim(10, 60), ax.set_ylim(8, 31), ax.set_aspect("equal")
            ax.set_xticks([]), ax.set_yticks([])
            for sp in ax.spines.values():
                sp.set_visible(False)
            ax.set_title(f"{NOMES[s]}, {vI(I).replace('{,}', ',')} mm/h, {t:g} h: FS {fs:.3f}".replace(".", ","),
                         fontsize=6.5, loc="left", color=TINTA)
    cb = fig.colorbar(pc, ax=axs, orientation="horizontal", fraction=0.03, pad=0.04, aspect=50,
                      ticks=NPW + [0.] + PPW)
    cb.set_label("$p_w$ (kPa); isolinhas: $-$80 a $-$20 e 50 a 200 kPa, e o nível d'água ($p_w$ = 0, mais grossa)",
                 fontsize=7, color=TINTA)
    cb.ax.tick_params(labelsize=6.5)
    fig.savefig(os.path.join(AQUI, "figs", "campos.pdf"), bbox_inches="tight", dpi=300)
    plt.close(fig)
N["nCampos"] = str(len(prontos))

# ---- as malhas refinadas do fator em 0 e 24 h (a banda de ruptura de cada instante), com o indicador plastico
# sqrt(J2(eps_p)) de cada triangulo relativo ao maior do instante, em faixas de escala log (ele varia ordens de grandeza).
# As faixas de >= 0,1 sao exatamente os triangulos que o criterio do refinamento marca (sem interpolacao)
from matplotlib.colors import BoundaryNorm, ListedColormap  # noqa: E402
CM = [(s, I) for (s, I) in (("f10", 900.), ("f50", 36.), ("f100", 36.)) if le_campo(s, I, 0.) and le_campo(s, I, 24.)]
NJ2 = [0.001, 0.01, 0.1, 0.5, 1.]
cmJ = ListedColormap([plt.cm.Reds(v) for v in (0.12, 0.32, 0.62, 0.88)])  # um tom so, do claro ao escuro
cmJ.set_under("white")
nJ = BoundaryNorm(NJ2, cmJ.N)
if CM:
    fig, axs = plt.subplots(len(CM), 2, figsize=(7.6, 2.1 * len(CM) + 0.5), squeeze=False)
    for i, (s, I) in enumerate(CM):
        for j, t in enumerate((0., 24.)):
            ax = axs[i][j]
            _, fat, fs = le_campo(s, I, t)
            Q = [[(a[0], a[1]), (a[2], a[3]), (a[4], a[5])] for a in fat]
            pj = PolyCollection(Q, array=fat[:, 6] / fat[:, 6].max(), cmap=cmJ, norm=nJ,
                                edgecolor=(0.3, 0.3, 0.29, 0.8), lw=0.22)
            ax.add_collection(pj)
            ax.set_xlim(10, 60), ax.set_ylim(5, 31), ax.set_aspect("equal")
            ax.set_xticks([]), ax.set_yticks([])
            ax.set_title(f"{NOMES[s]}, {vI(I).replace('{,}', ',')} mm/h, {t:g} h: {len(Q)} triângulos, FS {fs:.3f}"
                         .replace(".", ","), fontsize=7, loc="left", color=TINTA)
    fig.tight_layout(rect=(0, 0.11, 1, 1))
    cax = fig.add_axes([0.2, 0.085, 0.6, 0.016])
    cb = fig.colorbar(pj, cax=cax, orientation="horizontal", ticks=NJ2, extend="min", extendfrac=0.04)
    cb.ax.set_xticklabels([f"{v:g}".replace(".", ",") for v in NJ2])
    cb.ax.tick_params(labelsize=7)
    cb.set_label("$\\sqrt{J_2(\\varepsilon^p)}$ de cada triângulo, relativo ao maior do instante\n"
                 "$\\geq$ 0,1 (as duas faixas mais escuras): os triângulos que o critério do refinamento marca para dividir",
                 fontsize=7, color=TINTA)
    fig.savefig(os.path.join(AQUI, "figs", "malhas_fator.pdf"))
    plt.close(fig)
N["nMalhasF"] = str(len(CM))

# ---- o custo: por rodada (min) e a campanha no fermi (do primeiro inicio ao ultimo fim, eventos.log)
from datetime import datetime  # noqa: E402
st = [d["s"] / 60. for d in NOS.values() if d["fim"]]
if st:
    N["custoMedio"], N["custoMin"], N["custoMax"] = f"{np.mean(st):.0f}", f"{min(st):.0f}", f"{max(st):.0f}"
ev = os.path.join(FERMI, "eventos.log")
if os.path.exists(ev):
    hs = [datetime.strptime(lin[:19], "%Y-%m-%d %H:%M:%S") for lin in open(ev) if " inicio em " in lin or " fim em " in lin]
    if hs:
        N["campanha"] = f"{(max(hs) - min(hs)).total_seconds() / 60.:.0f}"

# ---- numeros.tex
N["nFeitos"] = str(sum(1 for d in NOS.values() if d["fim"]))
N["nTotal"] = str(len(CASOS))
with open(os.path.join(AQUI, "numeros.tex"), "w") as f:
    for k, x in N.items():
        f.write(f"\\newcommand{{\\{k}}}{{{x}}}\n")
print(f"pronto: {N['nFeitos']} de {N['nTotal']} rodadas, conferencia de malha {len(FINA)}")
