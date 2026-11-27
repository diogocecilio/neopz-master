#!/usr/bin/env python3
"""gera_dados.py: as figuras e tabelas do relatorio do Rahardjo2007 (relatorio_rahardjo2007.tex), com os dados que ja
existem (as rodadas que faltam ficam de fora, sem erro):
  ~/projects/runs/rahardjo2007/fermi/<id>/rah_<solo>_<I>_fs.csv e _nivel.csv  as 28 rodadas do fermi (cluster/casos.txt)
  ~/projects/runs/rahardjo2007/pc_malha/<id>_fina/...                          a conferencia de malha (casos_malha.txt)
  ~/projects/runs/rahardjo2007/malhas/m_transiente.csv, m_fator.csv            as malhas (Rahardjo2007 grava_malha=m)
  ../data/fig5_lido.csv, fig6_lido.csv, fig7_eq6.csv, tabela5.csv               o artigo (lido das figuras e tabelas)
  ~/projects/runs/rahardjo2007/fermi_i/<id>/rah_i_<solo>_<I>_fs.csv e _campo_   o modelo (i), ar ativo (casos_i.txt)
Saidas: figs/malhas.pdf, figs/fs_tempo.pdf (Fig. 5), figs/fs_min.pdf (Figs. 6 e 7), figs/nivel.pdf (o nivel d'agua no
FS minimo), tab/fs_min.tex, tab/eq6.tex, tab/malha.tex, tab/niveis.tex; o modelo (i): figs/fs_tempo_i.pdf,
figs/fs_min_i.pdf, figs/campos_i.pdf, tab/modelo_i.tex, tab/eq6_i.tex; e numeros.tex (as macros do texto)."""
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
COR_SOLO = {"f10": "#e34948", "f50": "#2a78d6", "f100": "#3a9a5b"}
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


def rodada(pasta, solo, I, mod=""):
    """mod: '' no (iii) (rah_<solo>_<I>), 'i' ou 'ii' nos outros (rah_<mod>_<solo>_<I>)"""
    base = os.path.join(pasta, f"rah_{mod + '_' if mod else ''}{solo}_{I:g}")
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
            "s_fator": [float(r["s_fator"]) for r in L], "pa_max": np.array([float(r.get("pa_max") or 0.) for r in L]),
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

# ---- figura: a curva de retencao (Fredlund-Xing, C = 1, n = m = 1), a permeabilidade a agua (k_s Theta^4) e a
# permeabilidade relativa ao ar (Brooks-Corey, lambda = 3, piso 1e-4) dos tres solos
psi = np.logspace(-2, 6, 600)  # o mesmo eixo da Fig. 3 do artigo
fig, axs = plt.subplots(1, 3, figsize=(7.6, 2.9))
for s, a_, ks_ in (("f10", 10., 1.e-4), ("f50", 50., 1.e-5), ("f100", 100., 1.e-6)):
    Sw = 1. / np.log(np.e + psi / a_)
    axs[0].semilogx(psi, 0.45 * Sw, color=COR_SOLO[s], lw=1.5, label=NOMES[s])
    axs[1].loglog(psi, ks_ * Sw ** 4, color=COR_SOLO[s], lw=1.5)
Se = np.linspace(0., 1., 300)
kra = np.maximum((1. - Se) ** 2 * (1. - Se ** (5. / 3.)), 1.e-4)
axs[2].semilogy(Se, kra, color=TINTA, lw=1.5, label="$k_{ra}$ (ar)")
axs[2].semilogy(Se, np.maximum(Se ** 4, 1.e-12), color=TINTA2, lw=1.2, ls=(0, (4, 2)), label="$k_{rw} = S_w^4$ (água)")
for ax in axs:
    eixo(ax)
for ax in axs[:2]:  # a faixa de succao das rodadas: ate 75 kPa (o limite do artigo); 135 kPa sem o limite
    ax.axvspan(1.e-2, 75., color="#ecebe7", zorder=0)
    ax.axvline(135., color=TINTA2, lw=0.8, ls=(0, (3, 2)), zorder=1)
axs[0].text(0.015, 0.03, "faixa das rodadas\n(até 75 kPa)", fontsize=6.5, color=TINTA2)
axs[0].text(170., 0.465, "135 kPa (sem o limite)", fontsize=6.5, color=TINTA2)
for s_, a_ in (("f10", 10.), ("f50", 50.), ("f100", 100.)):  # o menor S_w das rodadas (s = 75 kPa)
    axs[2].axvline(1. / np.log(np.e + 75. / a_), color=COR_SOLO[s_], lw=0.9, ls=(0, (2, 2)))
axs[0].set_xlabel("sucção $s$ (kPa)", fontsize=8.5, color=TINTA), axs[0].set_ylabel("$\\theta_w$", fontsize=9, color=TINTA)
axs[0].set_ylim(0., 0.5)
for ax in axs[:2]:
    ax.set_xticks([1.e-2, 1., 1.e2, 1.e4, 1.e6])
axs[1].set_xlabel("sucção $s$ (kPa)", fontsize=8.5, color=TINTA), axs[1].set_ylabel("$k_w$ (m/s)", fontsize=9, color=TINTA)
axs[2].set_xlabel("$S_w$", fontsize=8.5, color=TINTA), axs[2].set_ylabel("permeabilidade relativa", fontsize=8.5, color=TINTA)
axs[2].set_ylim(1.e-5, 1.5)
axs[0].set_title("retenção", fontsize=8.5, loc="left", color=TINTA)
axs[1].set_title("permeabilidade à água", fontsize=8.5, loc="left", color=TINTA)
axs[2].set_title("relativas (modelo (i))", fontsize=8.5, loc="left", color=TINTA)
axs[0].legend(fontsize=7, frameon=False, loc="center right")
axs[2].legend(fontsize=6.5, frameon=False, loc="lower center")
fig.tight_layout()
fig.savefig(os.path.join(AQUI, "figs", "solos.pdf"))
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


def isolinhas_pw(ax, tr, p):
    """as isolinhas de p_w (as de NPW e PPW finas, p_w = 0 grossa), so as que o campo atravessa"""
    for lv, lw in ((NPW, 0.45), (PPW, 0.45), ([0.], 1.1)):
        if p.min() < max(lv) and p.max() > min(lv):
            cs = ax.tricontour(tr, p, levels=lv, colors=TINTA, linewidths=lw, linestyles="-")
            ax.clabel(cs, fmt=lambda v: f"{v:.0f}", fontsize=4.5, inline=True, inline_spacing=1)
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
            isolinhas_pw(ax, tr, p)
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

# ======================================================================================================================
# ---- o modelo (i): o ar ativo, as mesmas 28 rodadas (fermi_i/<id>/rah_i_<solo>_<I>_fs.csv), contra o (iii)
FERMI_I = os.path.join(R, "fermi_i")
NOSI = {}
for (s, I) in CASOS:
    d = rodada(os.path.join(FERMI_I, f"{s}_I{I:g}"), s, I, "i")
    if d is not None and d["fim"]:
        NOSI[(s, I)] = d
SOLO_TEX = {"f10": "Dez", "f50": "Cinq", "f100": "Cem"}  # nos nomes das macros (sem digitos)
# o modelo (ii): acoplado (esqueleto elastico) com o ar passivo, as mesmas 28 rodadas (fermi_ii, rah_ii_<solo>_<I>)
FERMI_II = os.path.join(R, "fermi_ii")
NOSII = {}
for (s, I) in CASOS:
    d = rodada(os.path.join(FERMI_II, f"{s}_I{I:g}"), s, I, "ii")
    if d is not None and d["fim"]:
        NOSII[(s, I)] = d


def rot_I(s, I):
    return f"{vI(I).replace('{,}', ',')} mm/h ({I / KS[s]:g} $k_s$)".replace(".", ",")


# ---- figura: o FS e a maior p_a contra o tempo, (i) cheia e (iii) tracejada (as 7 intensidades da Fig. 5)
if NOSI:
    fig, axs = plt.subplots(3, 2, figsize=(7.6, 10.), sharex=True, gridspec_kw=dict(width_ratios=[1.45, 1.]))
    for k, s in enumerate(SOLOS):
        a1, a2 = axs[k]
        eixo(a1), eixo(a2)
        Is = sorted(I for (so, I) in F5 if so == s)
        cores = plt.cm.viridis(np.linspace(0.0, 0.88, len(Is)))
        for I, c in zip(Is, cores):
            d3, di = NOS.get((s, I)), NOSI.get((s, I))
            if d3 is not None:
                ok = d3["fs"] > 0
                a1.plot(d3["t"][ok], d3["fs"][ok], color=c, lw=0.9, ls=(0, (2.5, 1.5)), zorder=3)
            if di is not None:
                ok = di["fs"] > 0
                a1.plot(di["t"][ok], di["fs"][ok], color=c, lw=1.6, zorder=4)
                a2.plot(di["t"], di["pa_max"], color=c, lw=1.6)
            a1.plot([], [], color=c, lw=1.6, label=rot_I(s, I))
        for ax in (a1, a2):
            ax.axvline(24, color=TINTA2, lw=0.8, ls=(0, (4, 3)))
            ax.set_xlim(0, 48.5), ax.set_xticks(range(0, 49, 6))
        a1.axhline(1., color=TINTA2, lw=0.6)
        a1.set_ylim(0.4, 2.25), a2.set_ylim(0., 60.)  # a faixa de 0,4 a 0,8 fica livre para a legenda
        a1.set_ylabel("FS", fontsize=9, color=TINTA)
        a2.set_ylabel("maior $p_a$ no domínio (kPa)", fontsize=9, color=TINTA)
        a1.set_title(f"{NOMES[s]}, $k_s$ = {KS[s]:g} mm/h: o FS".replace(".", ","), fontsize=9, loc="left", color=TINTA)
        a2.set_title("a pressão do ar, modelo (i)", fontsize=9, loc="left", color=TINTA)
        a1.legend(fontsize=6.5, frameon=False, loc="lower left", ncol=2)
    axs[-1][0].set_xlabel("t (h)  (chuva até 24 h)", fontsize=9, color=TINTA)
    axs[-1][1].set_xlabel("t (h)", fontsize=9, color=TINTA)
    fig.text(0.99, 0.005, "linhas cheias: modelo (i), ar ativo; tracejadas: modelo (iii); FS de hora em hora",
             fontsize=7, color=TINTA2, ha="right")
    fig.tight_layout(rect=(0, 0.01, 1, 1))
    fig.savefig(os.path.join(AQUI, "figs", "fs_tempo_i.pdf"))
    plt.close(fig)

# ---- figura: o FS minimo contra I nos dois modelos, cada um com a sua Eq. 6 (no f10, C fixo no do artigo, como no iii)
AJI = {}
if NOSI:
    fig, axs = plt.subplots(1, 3, figsize=(7.6, 3.5), sharey=True)
    for ax, s in zip(axs, SOLOS):
        eixo(ax)
        ax.set_xscale("log")
        x = np.logspace(-1, 4, 400)
        Ia = sorted(I for (so, I) in F6 if so == s)
        pa6 = [float(EQ6[s][k]) for k in "ABCD"]
        Cf = pa6[2] if pa6[2] < 0. else None
        I3 = [I for I in Ia if (s, I) in NOS and NOS[(s, I)]["fim"]]
        ax.plot(I3, [minimo(NOS[(s, I)])[0] for I in I3], MARC[s], ms=5, mfc="none", mec=TINTA2, mew=1.0, zorder=4)
        if s in AJ:
            ax.plot(x, sig(x, *AJ[s][0]), color=TINTA2, lw=1.1, ls=(0, (5, 2)), zorder=3)
        Ii = [I for I in Ia if (s, I) in NOSI]
        Fi = [minimo(NOSI[(s, I)])[0] for I in Ii]
        ax.plot(Ii, Fi, MARC[s], ms=5, color=COR[s], zorder=5)
        if len(Ii) >= 5:
            p, r2 = ajusta(Ii, Fi, AJ[s][0] if s in AJ else pa6, Cf)
            if p is not None:
                AJI[s] = (p, r2, len(Ii), Cf is not None)
                ax.plot(x, sig(x, *p), color=COR[s], lw=1.5, zorder=4)
        ax.axvline(KS[s], color=GRADE, lw=1.2, zorder=0)
        ax.text(KS[s] * 1.15, 0.15, "$k_s$", fontsize=8, color=TINTA2)
        ax.set_xlim(0.1, 1.e4), ax.set_ylim(0., 2.5)
        ax.set_title(NOMES[s], fontsize=9, loc="left", color=TINTA)
        ax.set_xlabel("$I_r$ (mm/h)", fontsize=9, color=TINTA)
        ax.legend(handles=[Line2D([], [], color=COR[s], lw=1.5, marker=MARC[s], ms=5, label="(i), ar ativo"),
                           Line2D([], [], color=TINTA2, lw=1.1, ls=(0, (5, 2)), marker=MARC[s], mfc="none", mec=TINTA2,
                                  ms=5, label="(iii)")], fontsize=7, frameon=False, loc="lower left")
    axs[0].set_ylabel("FS mínimo (0 a 48 h)", fontsize=9, color=TINTA)
    fig.text(0.99, 0.01, "pontos: as rodadas; linhas: a Eq. 6 ajustada a cada modelo", fontsize=7, color=TINTA2,
             ha="right")
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(os.path.join(AQUI, "figs", "fs_min_i.pdf"))
    plt.close(fig)

# ---- a Eq. 6 ajustada ao (ii), como nos outros dois
AJII = {}
for s in SOLOS:
    Ia = sorted(I for (so, I) in F6 if so == s and (so, I) in NOSII)
    pa6 = [float(EQ6[s][k]) for k in "ABCD"]
    if len(Ia) >= 5:
        p, r2 = ajusta(Ia, [minimo(NOSII[(s, I)])[0] for I in Ia], AJ[s][0] if s in AJ else pa6,
                       pa6[2] if pa6[2] < 0. else None)
        if p is not None:
            AJII[s] = (p, r2, len(Ia), pa6[2] < 0.)

# ---- tabela: a Eq. 6 nos tres modelos (e as macros dos parametros, para o texto)
with open(os.path.join(AQUI, "tab", "eq6_i.tex"), "w") as f:
    f.write("\\begin{tabular}{ll|cccc|c|ccc}\n\\toprule\n")
    f.write("solo & modelo & $A$ & $B$ & $C$ & $D$ & $r^2$ & $I_{r(ini)}$ & $I_{r(thr)}$ & $I_{r(thr)}/k_s$\\\\\n"
            " & & & & (mm/h) & (mm/h) & & (mm/h) & (mm/h) & \\\\\n\\midrule\n")
    for s in SOLOS:
        for nome, mt, aj in (("(iii)", "Tres", AJ), ("(ii)", "Dois", AJII), ("(i)", "Um", AJI)):
            cab = f"{NOMES[s] if mt == 'Tres' else ''} & {nome}"
            if s not in aj:
                f.write(cab + " & \\multicolumn{8}{c}{(rodando)}\\\\\n")
                continue
            p, r2, n, fixo = aj[s]
            ii, it = ini_thr(*p)
            f.write(cab + f", {n} pontos & {v2(p[0])} & {v2(p[1])} & {vI(p[2]) + '$^\\dagger$' if fixo else v2(p[2])} & "
                    f"{v2(p[3])} & {v3(r2)} & {v2(ii)} & {v2(it)} & {v2(it / KS[s])}\\\\\n")
            sx = SOLO_TEX[s]
            N[f"A{mt}{sx}"], N[f"D{mt}{sx}"] = v2(p[0]), v2(abs(p[3]))
            N[f"Ini{mt}{sx}"], N[f"Thr{mt}{sx}"], N[f"Thrks{mt}{sx}"] = v2(ii), v2(it), v2(it / KS[s])
        f.write("\\midrule\n" if s != SOLOS[-1] else "")
    f.write("\\bottomrule\n\\end{tabular}\n")

# ---- tabela: cada rodada nos dois modelos: o FS minimo, t_1, a maior p_a, a infiltracao e o alagamento
with open(os.path.join(AQUI, "tab", "modelo_i.tex"), "w") as f:
    f.write("\\begin{tabular}{lrr|cccc|c|cc|c|cc|cc}\n\\toprule\n")
    f.write(" & & & \\multicolumn{4}{c|}{FS mínimo (0--48~h)} & instante & \\multicolumn{2}{c|}{$t_1$ (h)} & "
            "maior $p_a$ & \\multicolumn{2}{c|}{infiltrou até 24~h (\\%)} & \\multicolumn{2}{c}{alaga em (h)}\\\\\n")
    f.write("solo & $I_r$ & $I_r/k_s$ & (iii) & (ii) & (i) & (i) $-$ (iii) & do mín. (i) & (iii) & (i) & (i), kPa & (iii) & "
            "(i) & (iii) & (i)\\\\\n\\midrule\n")
    queda = {}
    for k, (s, I) in enumerate(CASOS):
        if k and s != CASOS[k - 1][0]:
            f.write("\\midrule\n")
        d3, di = NOS.get((s, I)), NOSI.get((s, I))
        lin = f"{NOMES[s]} & {vI(I)} & {vI(round(I / KS[s], 4))} & "
        if d3 is None or not d3["fim"] or di is None:
            f.write(lin + "\\multicolumn{12}{c}{(rodando)}\\\\\n")
            continue
        f3, fi, ti = minimo(d3)[0], *minimo(di)[:2]
        f2 = minimo(NOSII[(s, I)])[0] if (s, I) in NOSII else None
        dlt = fi - f3
        if dlt < queda.get(s, (0., 0.))[0]:
            queda[s] = (dlt, I)
        t13, t1i = cruza(d3)[0], cruza(di)[0]
        chuva = I / 1000. * 24. * LSUP
        inf = [100. * d["infil"][list(d["t"]).index(24.)] / chuva for d in (d3, di)]
        al = [next((t for t, n in zip(d["t"], d["face"]) if n > 0), None) for d in (d3, di)]
        sd = "0{,}00" if abs(dlt) < 0.005 else ("$+$" if dlt > 0 else "$-$") + v2(abs(dlt))
        f.write(lin + f"{v3(f3)} & {v3(f2) if f2 else '--'} & {v3(fi)} & {sd} & {ti:g} & {v2(t13) if t13 else '--'} & {v2(t1i) if t1i else '--'} & "
                f"{di['pa_max'].max():.1f} & {inf[0]:.0f} & {inf[1]:.0f} & "
                f"{'%g' % al[0] if al[0] is not None else '--'} & {'%g' % al[1] if al[1] is not None else '--'}\\\\\n"
                .replace(".", "{,}"))
    f.write("\\bottomrule\n\\end{tabular}\n")
for s, (dq, Iq) in queda.items():
    N[f"queda{SOLO_TEX[s]}"], N[f"queda{SOLO_TEX[s]}I"] = v2(-dq), vI(Iq)
# o (ii) contra o (iii) (a deformacao) e o (i) contra o (ii) (o ar)
DEF = {k: minimo(NOSII[k])[0] - minimo(NOS[k])[0] for k in NOSII if k in NOS and NOS[k]["fim"]}
AR = {k: minimo(NOSI[k])[0] - minimo(NOSII[k])[0] for k in NOSII if k in NOSI}
N["iiFeitos"] = str(len(NOSII))
if DEF:
    kd = max(DEF, key=DEF.get)
    N["defMax"], N["defMin"] = v2(max(DEF.values())), v3(min(DEF.values()))
    N["defMaxSolo"], N["defMaxI"] = NOMES[kd[0]], vI(kd[1])
if AR:
    ka, kb = min(AR, key=AR.get), max(AR, key=AR.get)
    N["arMin"], N["arMinSolo"], N["arMinI"] = v2(-AR[ka]), NOMES[ka[0]], vI(ka[1])
    N["arMax"], N["arMaxSolo"], N["arMaxI"] = v2(AR[kb]), NOMES[kb[0]], vI(kb[1])

# ---- figura: a diferenca de FS minimo separada em deformacao ((ii) - (iii)) e ar ((i) - (ii)), contra I/k_s
if DEF and AR:
    fig, axs = plt.subplots(1, 3, figsize=(7.6, 3.0), sharey=True)
    for ax, s in zip(axs, SOLOS):
        eixo(ax)
        ax.set_xscale("log")
        Is = sorted(I for (so, I) in DEF if so == s)
        xs_ = [I / KS[s] for I in Is]
        ax.axhline(0., color=TINTA2, lw=0.7)
        ax.plot(xs_, [DEF[(s, I)] for I in Is], "s-", color=TINTA2, ms=4, lw=1.2, mfc="white",
                label="deformação: (ii) $-$ (iii)")
        Ia = [I for I in Is if (s, I) in AR]
        ax.plot([I / KS[s] for I in Ia], [AR[(s, I)] for I in Ia], "o-", color=COR[s], ms=4, lw=1.5,
                label="ar: (i) $-$ (ii)")
        ax.set_title(NOMES[s], fontsize=9, loc="left", color=TINTA)
        ax.set_xlabel("$I_r/k_s$", fontsize=9, color=TINTA)
        ax.set_xlim(0.015, 15), ax.set_ylim(-0.42, 0.12)
        ax.legend(fontsize=6.5, frameon=False, loc="lower right" if s == "f50" else "lower left")
    axs[0].set_ylabel("diferença de FS mínimo", fontsize=9, color=TINTA)
    fig.tight_layout()
    fig.savefig(os.path.join(AQUI, "figs", "decomposicao.pdf"))
    plt.close(fig)
N["iFeitos"] = str(len(NOSI))
N["iMinVinte"] = str(sum(1 for d in NOSI.values() if minimo(d)[1] == 24.))  # o FS minimo do (i) no fim da chuva
N["iPamax"] = f"{max(d['pa_max'].max() for d in NOSI.values()):.0f}" if NOSI else "--"

# ---- figura: os campos em 24 h (o instante do FS minimo do (i)): p_w no (iii), p_w no (i) e p_a no (i). No (i) a
# linha grossa (p_w = 0) deixa de ser o limite da zona saturada: o limite e s = p_a - p_w = 0 (tracejada)


def le_campo_i(s, I, t):
    base = os.path.join(FERMI_I, f"{s}_I{I:g}", f"rah_i_{s}_{I:g}")
    C = le(f"{base}_campo_{t:g}h.csv")
    if not C or "a1" not in C[0]:
        return None
    tri = np.array([[float(r[k]) for k in ("x1", "y1", "p1", "x2", "y2", "p2", "x3", "y3", "p3", "a1", "a2", "a3")]
                    for r in C])
    fs = {float(r["t_h"]): float(r["FS"]) for r in le(f"{base}_fs.csv")}
    return tri, fs.get(t)


def nodal(tri, cols):
    """a Triangulation dos nos unicos e o valor de cada no (as colunas cols de cada triangulo, campo continuo)"""
    tr, T, _ = triangulacao(tri[:, [0, 1, 3, 4, 6, 7]].reshape(-1, 3, 2))
    v = np.zeros(len(tr.x))
    v[T.ravel()] = tri[:, cols].ravel()
    return tr, v


CI = [("f10", 90.), ("f10", 900.), ("f50", 36.), ("f100", 9.)]
TI = 24.
prontosI = [(s, I) for (s, I) in CI if le_campo(s, I, TI) is not None and le_campo_i(s, I, TI) is not None]
if prontosI:
    norm = TwoSlopeNorm(vmin=-85., vcenter=0., vmax=225.)
    normA = plt.Normalize(0., 60.)
    PA = [10., 20., 30., 40., 50.]
    fig, axs = plt.subplots(len(prontosI), 3, figsize=(7.6, 1.6 * len(prontosI) + 0.9), squeeze=False)
    for i, (s, I) in enumerate(prontosI):
        tri3, _, fs3 = le_campo(s, I, TI)
        triI, fsI = le_campo_i(s, I, TI)
        tr3, p3 = nodal(tri3, [2, 5, 8])
        trI, pI = nodal(triI, [2, 5, 8])
        _, aI = nodal(triI, [9, 10, 11])
        sI = aI - pI  # a succao no (i)
        amax = aI[sI > 0.].max() if (sI > 0.).any() else 0.  # a maior p_a na zona nao saturada
        for j, ax in enumerate(axs[i]):
            if j < 2:
                tr, p = (tr3, p3) if j == 0 else (trI, pI)
                pc = ax.tripcolor(tr, p, cmap="RdBu_r", norm=norm, shading="gouraud", rasterized=True)
                isolinhas_pw(ax, tr, p)
            else:
                sat = (sI[trI.triangles] <= 0.).any(axis=1)  # o triangulo toca a zona saturada
                V = np.stack([trI.x[trI.triangles], trI.y[trI.triangles]], axis=-1)
                ax.add_collection(PolyCollection(V[sat], facecolor="#dddcd8", edgecolor="none", rasterized=True))
                trA = Triangulation(trI.x, trI.y, trI.triangles, mask=sat)
                pcA = ax.tripcolor(trA, aI, cmap="Purples", norm=normA, shading="gouraud", rasterized=True)
                lv = [v for v in PA if v < amax]
                if lv:
                    cs = ax.tricontour(trA, aI, levels=lv, colors=TINTA, linewidths=0.45)
                    ax.clabel(cs, fmt=lambda v: f"{v:.0f}", fontsize=4.5, inline=True, inline_spacing=1)
            if j >= 1 and sI.min() < 0. < sI.max():
                ax.tricontour(trI, sI, levels=[0.], colors=TINTA, linewidths=1.0, linestyles=[(0, (3, 1.5))])
            ax.plot([0, 30, 40, 70], [30, 30, 20, 20], color=TINTA, lw=0.6)
            ax.set_xlim(10, 60), ax.set_ylim(8, 31), ax.set_aspect("equal")
            ax.set_xticks([]), ax.set_yticks([])
            for sp in ax.spines.values():
                sp.set_visible(False)
        cab = f"{NOMES[s]}, {vI(I).replace('{,}', ',')} mm/h, 24 h. "
        axs[i][0].set_title((cab + f"(iii): $p_w$, FS {fs3:.3f}").replace(".", ","), fontsize=6.5, loc="left",
                            color=TINTA)
        axs[i][1].set_title(f"(i): $p_w$, FS {fsI:.3f}".replace(".", ","), fontsize=6.5, loc="left", color=TINTA)
        axs[i][2].set_title(f"(i): $p_a$, até {amax:.0f} kPa na zona não saturada", fontsize=6.5, loc="left",
                            color=TINTA)
    fig.subplots_adjust(left=0.01, right=0.99, top=0.95, bottom=0.13, wspace=0.04, hspace=0.32)
    c1 = fig.colorbar(pc, cax=fig.add_axes([0.04, 0.075, 0.6, 0.012]), orientation="horizontal", ticks=NPW + [0.] + PPW)
    c1.set_label("$p_w$ (kPa): isolinhas a cada 20 kPa na sucção e a cada 50 kPa acima de zero; grossa: $p_w$ = 0;\n"
                 "tracejada: $s = p_a - p_w$ = 0 (no (i), o limite da zona saturada)", fontsize=6.5, color=TINTA)
    c2 = fig.colorbar(pcA, cax=fig.add_axes([0.7, 0.075, 0.27, 0.012]), orientation="horizontal", ticks=[0] + PA + [60])
    c2.set_label("$p_a$ (kPa), na zona não saturada\n(cinza: saturada, $s \\leq 0$)", fontsize=6.5, color=TINTA)
    for c in (c1, c2):
        c.ax.tick_params(labelsize=6)
    fig.savefig(os.path.join(AQUI, "figs", "campos_i.pdf"), dpi=300)
    plt.close(fig)
N["nCamposI"] = str(len(prontosI))

# ---- figura: linhas de controle, p_w contra a profundidade em duas verticais (sob a crista, x = 20 m, e no meio da
# face, x = 35 m) em 0, 6, 12, 24 e 48 h, (i) cheia e (iii) tracejada; e a taxa media de variacao de p_w entre os
# instantes gravados (kPa/h), na mesma vertical
from matplotlib.tri import LinearTriInterpolator  # noqa: E402
TP = [0., 6., 12., 24., 48.]
VERT = [(20., "sob a crista, $x$ = 20 m"), (35., "meio da face, $x$ = 35 m")]
ZP = np.linspace(0.02, 15., 120)
CP = [(s, I) for (s, I) in CI if all(le_campo(s, I, t) is not None and le_campo_i(s, I, t) is not None for t in TP)]


def perfil(tr, v, x):
    f = LinearTriInterpolator(tr, v)
    return np.array([float(f(x, superficie(x) - z)) for z in ZP])


if CP:
    corT = dict(zip(TP, plt.cm.viridis(np.linspace(0., 0.9, len(TP)))))
    PF, PA = {}, {}
    for (s, I) in CP:
        for t in TP:
            tri3, triI = le_campo(s, I, t)[0], le_campo_i(s, I, t)[0]
            tr3, p3 = nodal(tri3, [2, 5, 8])
            trI, pI = nodal(triI, [2, 5, 8])
            _, aI = nodal(triI, [9, 10, 11])
            for x, _ in VERT:
                PF[(s, I, t, x)] = (perfil(tr3, p3, x), perfil(trI, pI, x))
                PA[(s, I, t, x)] = perfil(trI, aI, x)
    for nome, taxa in (("perfis_i", False), ("taxa_i", True)):
        fig, axs = plt.subplots(len(CP), 2, figsize=(7.4, 2.5 * len(CP) + 0.6), sharey=True, squeeze=False)
        for i, (s, I) in enumerate(CP):
            for j, (x, tx) in enumerate(VERT):
                ax = axs[i][j]
                eixo(ax)
                ts = TP if not taxa else TP[1:]
                for k, t in enumerate(ts):
                    a3, aI = PF[(s, I, t, x)]
                    if taxa:  # a taxa media no intervalo que termina em t
                        t0 = TP[TP.index(t) - 1]
                        b3, bI = PF[(s, I, t0, x)]
                        a3, aI = (a3 - b3) / (t - t0), (aI - bI) / (t - t0)
                        rot = f"{t0:g}–{t:g} h"
                    else:
                        rot = f"{t:g} h"
                    ax.plot(a3, ZP, color=corT[t], lw=1.0, ls=(0, (3, 1.5)))
                    ax.plot(aI, ZP, color=corT[t], lw=1.5, label=rot)
                ax.axvline(0., color=TINTA2, lw=0.7)
                ax.set_ylim(15, 0)
                ax.set_title(f"{NOMES[s]}, {vI(I).replace('{,}', ',')} mm/h, {tx}", fontsize=8, loc="left", color=TINTA)
                if j == 0:
                    ax.set_ylabel("profundidade (m)", fontsize=8, color=TINTA)
                if i == len(CP) - 1:
                    ax.set_xlabel("taxa de variação de $p_w$ (kPa/h)" if taxa else "$p_w$ (kPa)", fontsize=8,
                                  color=TINTA)
            axs[i][0].legend(fontsize=6.5, frameon=False, loc="lower left" if not taxa else "lower right")
        fig.text(0.99, 0.005, "cheias: modelo (i), ar ativo; tracejadas: modelo (iii)" +
                 ("; taxa média entre os instantes gravados" if taxa else ""), fontsize=7, color=TINTA2, ha="right")
        fig.tight_layout(rect=(0, 0.015, 1, 1))
        fig.savefig(os.path.join(AQUI, "figs", nome + ".pdf"))
        plt.close(fig)

    # p_a (so no i) e o grau de saturacao (Fredlund-Xing com C = 1, n = m = 1: S = 1/ln(e + s/a), s = p_a - p_w > 0)
    ASOLO = {"f10": 10., "f50": 50., "f100": 100.}

    def satur(suc, a):
        return np.where(suc > 0., 1. / np.log(np.e + np.maximum(suc, 0.) / a), 1.)

    fig, axs = plt.subplots(len(CP), 4, figsize=(7.6, 2.5 * len(CP) + 0.6), sharey=True, squeeze=False)
    for i, (s, I) in enumerate(CP):
        for j, (x, tx) in enumerate(VERT):
            aa, ab = axs[i][j], axs[i][2 + j]
            eixo(aa), eixo(ab)
            for t in TP:
                p3, pI = PF[(s, I, t, x)]
                aI = PA[(s, I, t, x)]
                aa.plot(aI, ZP, color=corT[t], lw=1.5, label=f"{t:g} h")
                ab.plot(satur(-p3, ASOLO[s]), ZP, color=corT[t], lw=1.0, ls=(0, (3, 1.5)))
                ab.plot(satur(aI - pI, ASOLO[s]), ZP, color=corT[t], lw=1.5)
            aa.set_xlim(-2, 60), ab.set_xlim(0.4, 1.01)
            aa.set_ylim(15, 0)
            lado = "crista" if x < 30 else "face"
            aa.set_title(f"$p_a$ (i), {lado}", fontsize=7.5, loc="left", color=TINTA)
            ab.set_title(f"$S_w$, {lado}", fontsize=7.5, loc="left", color=TINTA)
            if i == len(CP) - 1:
                aa.set_xlabel("$p_a$ (kPa)", fontsize=8, color=TINTA)
                ab.set_xlabel("grau de saturação", fontsize=8, color=TINTA)
        axs[i][0].set_ylabel(f"{NOMES[s]}, {vI(I).replace('{,}', ',')} mm/h\nprofundidade (m)", fontsize=8, color=TINTA)
    axs[0][0].legend(fontsize=6.5, frameon=False, loc="lower right")
    fig.text(0.99, 0.005, "crista: x = 20 m; face: x = 35 m; S_w: cheias (i), tracejadas (iii)", fontsize=7,
             color=TINTA2, ha="right")
    fig.tight_layout(rect=(0, 0.015, 1, 1))
    fig.savefig(os.path.join(AQUI, "figs", "perfis_ar_i.pdf"))
    plt.close(fig)
N["nPerfis"] = str(len(CP))

# ---- a duracao da chuva: o FS no fim de uma chuva de D horas e o FS em t = D da chuva de 24 h (a mesma historia ate
# D). A Eq. 6 ajustada ao FS(D) contra I, para D = 3, 6, 12 e 24 h, nos dois modelos
DUR = [3., 6., 12., 24.]
AJD = {}
for s in SOLOS:
    for mt, NN in (("(iii)", NOS), ("(i)", NOSI)):
        Is = sorted(I for (so, I) in CASOS if so == s and (so, I) in NN and NN[(so, I)]["fim"])
        for D in DUR:
            F = [dict(zip(NN[(s, I)]["t"], NN[(s, I)]["fs"]))[D] for I in Is]
            pa6 = [float(EQ6[s][k]) for k in "ABCD"]
            if len(Is) >= 5:
                p, r2 = ajusta(Is, F, AJ[s][0] if s in AJ else pa6, pa6[2] if pa6[2] < 0. else None)
                if p is not None:
                    AJD[(s, mt, D)] = (p, r2, Is, F)
if AJD:
    fig, axs = plt.subplots(1, 3, figsize=(7.6, 3.6), sharey=True)
    corD = dict(zip(DUR, plt.cm.viridis(np.linspace(0.85, 0.0, len(DUR)))))
    x = np.logspace(-1, 4, 400)
    for ax, s in zip(axs, SOLOS):
        eixo(ax)
        ax.set_xscale("log")
        for D in DUR:
            for mt, ls, lw, mfc in (("(iii)", (0, (3, 1.5)), 1.0, "none"), ("(i)", "-", 1.5, None)):
                if (s, mt, D) in AJD:
                    p, r2, Is, F = AJD[(s, mt, D)]
                    ax.plot(Is, F, MARC[s], ms=3.5, color=corD[D], mfc=mfc or corD[D], mec=corD[D])
                    ax.plot(x, sig(x, *p), color=corD[D], lw=lw, ls=ls, label=f"{D:g} h" if mt == "(i)" else None)
        ax.axvline(KS[s], color=GRADE, lw=1.2, zorder=0)
        ax.set_xlim(0.1, 1.e4), ax.set_ylim(0., 2.5)
        ax.set_title(NOMES[s], fontsize=9, loc="left", color=TINTA)
        ax.set_xlabel("$I_r$ (mm/h)", fontsize=9, color=TINTA)
    axs[0].set_ylabel("FS no fim da chuva", fontsize=9, color=TINTA)
    axs[0].legend(title="duração", title_fontsize=7, fontsize=7, frameon=False, loc="lower left")
    fig.text(0.99, 0.01, "cheias: (i), ar ativo; tracejadas e vazios: (iii); linhas: a Eq. 6 ajustada a cada duração",
             fontsize=7, color=TINTA2, ha="right")
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(os.path.join(AQUI, "figs", "duracao.pdf"))
    plt.close(fig)
    with open(os.path.join(AQUI, "tab", "duracao.tex"), "w") as f:
        f.write("\\begin{tabular}{l|cc|cc|cc|cc}\n\\toprule\n & \\multicolumn{8}{c}{duração da chuva}\\\\\n"
                "solo, grandeza & \\multicolumn{2}{c|}{3 h} & \\multicolumn{2}{c|}{6 h} & \\multicolumn{2}{c|}{12 h} & "
                "\\multicolumn{2}{c}{24 h}\\\\\n & (iii) & (i) & (iii) & (i) & (iii) & (i) & (iii) & (i)\\\\\n\\midrule\n")
        for s in SOLOS:
            for nome, fn in (("$A$", lambda p: v2(p[0])), ("$I_{r(thr)}/k_s$", lambda p: v2(ini_thr(*p)[1] / KS[s])),
                             ("$I_{r(ini)}/k_s$", lambda p: v2(ini_thr(*p)[0] / KS[s]))):
                f.write(f"{NOMES[s]}, {nome} & " + " & ".join(fn(AJD[(s, mt, D)][0]) if (s, mt, D) in AJD else "--"
                                                           for D in DUR for mt in ("(iii)", "(i)")) + "\\\\\n")
            f.write("\\midrule\n" if s != SOLOS[-1] else "")
        f.write("\\bottomrule\n\\end{tabular}\n")
    err = [NN[k]["fs"][list(NN[k]["t"]).index(24.)] - minimo(NN[k])[0] for NN in (NOS, NOSI) for k in NN
           if NN[k]["fim"]]
    N["durErro"] = v3(max(err))

# ---- o custo do (i) e do (ii): so as rodadas da primeira leva (os jobs que comecaram no primeiro minuto, 7 processos de
# 3 threads por no). No (i), o f100 a 5,4 mm/h rodou depois, sozinho num no com 12 threads (a trava dele ficou vazia)


def leva(pasta):
    """os casos rodados pelos jobs da primeira leva e a duracao da campanha (min), do eventos.log"""
    ev = os.path.join(pasta, "eventos.log")
    if not os.path.exists(ev):
        return None, None
    L = [(datetime.strptime(lin[:19], "%Y-%m-%d %H:%M:%S"), lin.split("[")[1].split("]")[0], lin) for lin in open(ev)
         if "[" in lin]
    ini = {j: h for h, j, lin in L if " inicio em " in lin}
    if not ini:
        return None, None
    t0 = min(ini.values())
    jobs = {j for j, h in ini.items() if (h - t0).total_seconds() < 60.}
    casos = {lin.split("] ")[1].split()[0] for h, j, lin in L if j in jobs and " codigo " in lin}
    fins = [h for h, j, lin in L if j in jobs and " fim em " in lin]
    return casos, ((max(fins) - t0).total_seconds() / 60. if fins else None)


casos_i, camp_i = leva(FERMI_I)
DI = [d for (s, I), d in NOSI.items() if casos_i is None or f"{s}_I{I:g}" in casos_i]
if DI:
    sti = [d["s"] / 60. for d in DI]
    tr_i = [d["s_trans"][-1] for d in DI]
    tr_3 = [d["s_trans"][-1] for d in NOS.values() if d["fim"]]
    sf_i = [x for d in DI for x, fv in zip(d["s_fator"], d["fs"]) if fv > 0]
    N["iCustoMin"], N["iCustoMax"], N["iCustoMedio"] = f"{min(sti):.0f}", f"{max(sti):.0f}", f"{np.mean(sti):.0f}"
    N["iTrans"], N["iiiTrans"] = f"{np.mean(tr_i) / 60.:.0f}", f"{np.mean(tr_3) / 60.:.0f}"
    N["iTransRazao"] = v2(np.mean(tr_i) / np.mean(tr_3))
    N["iFator"] = f"{np.mean(sf_i):.0f}"
    N["iLeva"] = str(len(DI))
if camp_i:
    N["iCampanha"] = f"{camp_i:.0f}"
outros = [d for (s, I), d in NOSI.items() if casos_i is not None and f"{s}_I{I:g}" not in casos_i]
if outros:
    N["iRefeita"] = f"{outros[0]['s'] / 60.:.0f}"
casos_ii, camp_ii = leva(FERMI_II)
if NOSII:
    st2 = [d["s"] / 60. for d in NOSII.values()]
    N["iiCustoMin"], N["iiCustoMax"], N["iiCustoMedio"] = f"{min(st2):.0f}", f"{max(st2):.0f}", f"{np.mean(st2):.0f}"
    N["iiTrans"] = f"{np.mean([d['s_trans'][-1] for d in NOSII.values()]) / 60.:.0f}"
ev2 = os.path.join(FERMI_II, "eventos.log")
if os.path.exists(ev2):
    hs = [datetime.strptime(lin[:19], "%Y-%m-%d %H:%M:%S") for lin in open(ev2) if " inicio em " in lin or " fim em " in lin]
    if hs:
        N["iiCampanha"] = f"{(max(hs) - min(hs)).total_seconds() / 60.:.0f}"

# ---- numeros.tex
N["nFeitos"] = str(sum(1 for d in NOS.values() if d["fim"]))
N["nTotal"] = str(len(CASOS))
with open(os.path.join(AQUI, "numeros.tex"), "w") as f:
    for k, x in N.items():
        f.write(f"\\newcommand{{\\{k}}}{{{x}}}\n")
print(f"pronto: {N['nFeitos']} de {N['nTotal']} rodadas, conferencia de malha {len(FINA)}")
