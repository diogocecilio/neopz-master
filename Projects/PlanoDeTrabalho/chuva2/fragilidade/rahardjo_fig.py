#!/usr/bin/env python3
"""rahardjo_fig.py: as Figs. 6 e 7 de Rahardjo et al. (2007) para os tres solos (f10, f50, f100) com as nossas contas,
chamado pelo gera_dados.py (gera(...)):
  - o FS das medias no fim de uma chuva de 24 h contra I (o F_s(min) do artigo, que cai no fim da chuva), das varreduras
    ~/projects/runs/sigmoide_f50 (f50), ~/projects/runs/sigmoide_f100_72h (f100: o FS de 24 h das chuvas de 72 h) e
    ~/projects/runs/rahardjo_3solos (f10 e os patamares que faltavam; <variante>_<solo>_I<I>, e <variante>_<solo>_fs0 com o
    FS inicial), nas variantes (iii) e (i);
  - a Eq. 6 ajustada na forma do artigo, F = A + B / (1 + exp(-(I - C)/D)) em I linear, com D < 0 (A e o patamar de
    baixo), e o r^2;
  - a construcao das tangentes no grafico semilog: a tangente no ponto mais ingreme de F(log10 I) corta o patamar de cima
    (F em I = 0) em I_r(ini) e o de baixo (A) em I_r(thr); V_r = 24 I e I / k_s;
  - as constantes publicadas (Fig. 7 do artigo) para comparar.
Saidas: figs/rahardjo_fig6.pdf|png, figs/rahardjo_fig7.pdf|png, tab/rahardjo_const.tex; devolve as macros."""
import math
import os
import sys

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from scipy.optimize import curve_fit

RUNS = os.path.expanduser("~/projects/runs")
S50 = os.path.join(RUNS, "sigmoide_f50")
S100 = os.path.join(RUNS, "sigmoide_f100_72h")
R3 = os.path.join(RUNS, "rahardjo_3solos")
sys.path.insert(0, os.path.expanduser("~/projects/neopz-master/Projects/SlopeRandomFieldRefinement/cluster"))
import resumo_frag as rf  # noqa: E402

SOLOS = [("f10", 360.), ("f50", 36.), ("f100", 3.6)]  # k_s em mm/h (1e-4, 1e-5 e 1e-6 m/s)
# as cores dos solos: categoricas fora do azul e do verde das variantes (validador da paleta: passa; o magenta tem pouco
# contraste com o fundo, entao todo solo tem rotulo direto e a tabela tem os numeros)
CORS = {"f10": "#eb6834", "f50": "#6250d6", "f100": "#e87ba4"}
# Rahardjo et al. (2007), Fig. 7: as constantes da Eq. 6, o r^2 e as intensidades das tangentes (mm/h; no f100 o texto
# diz 3,8, a figura diz 2,8 = 0,778 k_s, que e a conta certa)
RAH = {"f10": {"A": 0.81, "B": 13.54, "C": -118., "D": -52., "r2": 0.999, "Iini": 11., "Ithr": 140.},
       "f50": {"A": 0.82, "B": 1.83, "C": 10.7, "D": -13.13, "r2": 0.997, "Iini": 8.5, "Ithr": 51.},
       "f100": {"A": 1.45, "B": 0.70, "C": 4.75, "D": -2.35, "r2": 0.998, "Iini": 2.8, "Ithr": 11.}}
NOMEV = {"iii": "(iii)", "i": "(i)"}
XMIN, XMAX = 0.1, 1e4


def num(x, nd=2):
    return f"{x:.{nd}f}".replace(".", "{,}")


def eq6(I, A, B, C, D):
    z = np.clip(-(np.asarray(I, float) - C) / D, -700, 700)
    return A + B / (1 + np.exp(z))


def le_F(solo, var):
    """{I: FS no fim da chuva de 24 h} das varreduras e o FS inicial (ou None)"""
    F = {}
    if solo == "f50":
        F.update(rf.le_det(S50, "" if var == "iii" else "i_", 24.))
    if solo == "f100":
        F.update(rf.le_det(S100, var + "_", 24.))
    if os.path.isdir(R3):
        F.update(rf.le_det(R3, f"{var}_{solo}_", 24.))
    fs0 = None
    a = os.path.join(R3, f"{var}_{solo}_fs0", "det_0_-1.csv")
    if os.path.exists(a):
        import csv
        for r in csv.DictReader(open(a)):
            if abs(float(r["t_h"])) < 1e-9 and r["situacao"] == "ok":
                fs0 = float(r["FS_min"])
    return dict(sorted(F.items())), fs0


def ajusta(F):
    """a Eq. 6 na forma do artigo (D < 0, B > 0, A = o patamar de baixo) e o r^2; None com menos de 6 pontos"""
    if len(F) < 6:
        return None
    I, y = np.array(list(F), float), np.array(list(F.values()), float)
    meio = (y.max() + y.min()) / 2
    melhor = None
    for Dc in (-0.02, -0.05, -0.1, -0.3):  # varias partidas: a escala D em fracoes da faixa de I
        p0 = [y.min(), y.max() - y.min(), float(I[np.argmin(abs(y - meio))]), Dc * (I.max() - I.min())]
        try:
            p, _ = curve_fit(eq6, I, y, p0=p0, maxfev=200000)
        except (RuntimeError, ValueError):
            continue
        ssr = float(np.sum((y - eq6(I, *p)) ** 2))
        if melhor is None or ssr < melhor[1]:
            melhor = (p, ssr)
    if melhor is None:
        return None
    A, B, C, D = melhor[0]
    if D > 0:  # a mesma curva na convencao do artigo
        A, B, D = A + B, -B, -D
    r2 = 1 - melhor[1] / float(np.sum((y - y.mean()) ** 2))
    return {"A": A, "B": B, "C": C, "D": D, "r2": r2}


def tangentes(p):
    """no grafico semilog: a tangente no ponto mais ingreme de F(log10 I) e os cortes com o patamar de cima (F(0)) e o de
    baixo (A): I_r(ini) e I_r(thr)"""
    x = np.linspace(math.log10(XMIN) - 1, math.log10(XMAX) + 1, 12001)
    F = eq6(10 ** x, p["A"], p["B"], p["C"], p["D"])
    dF = np.gradient(F, x)
    k = int(np.argmax(np.abs(dF)))
    Fup, Flo = float(eq6(0., p["A"], p["B"], p["C"], p["D"])), p["A"]
    xi, xt = x[k] + (Fup - F[k]) / dF[k], x[k] + (Flo - F[k]) / dF[k]
    return {"Iini": 10 ** xi, "Ithr": 10 ** xt, "xk": x[k], "Fk": float(F[k]), "dF": float(dF[k]), "Fup": Fup,
            "Flo": Flo}


def I_de_F(p, alvo):
    """a intensidade em que a Eq. 6 vale alvo (None se nao cruza)"""
    I = np.logspace(math.log10(XMIN), math.log10(XMAX * 10), 20001)
    F = eq6(I, p["A"], p["B"], p["C"], p["D"])
    s = np.where(np.diff(np.sign(F - alvo)) != 0)[0]
    if not len(s):
        return None
    j = s[0]
    return float(I[j] + (alvo - F[j]) * (I[j + 1] - I[j]) / (F[j + 1] - F[j]))


def gera(k):
    """k: os utilitarios do gera_dados.py (estilo, salva, TINTA, TINTA2, GRADE, TABS). Devolve as macros."""
    estilo, salva, TINTA, TINTA2, TABS = k["estilo"], k["salva"], k["TINTA"], k["TINTA2"], k["TABS"]
    M = {}
    dados = {}
    for s, ks in SOLOS:
        for v in ("iii", "i"):
            F, fs0 = le_F(s, v)
            p = ajusta(F)
            dados[(s, v)] = {"F": F, "fs0": fs0, "p": p, "t": tangentes(p) if p else None,
                             "Ic": I_de_F(p, 1.) if p else None}
    rah_t = {s: tangentes(RAH[s]) for s, _ in SOLOS}

    # ---- Fig. 6: os tres solos num grafico semilog, uma variante por painel
    fig, axs = plt.subplots(1, 2, figsize=(7.4, 3.5), sharey=True)
    Ix = np.logspace(math.log10(XMIN), math.log10(XMAX), 600)
    for ax, v in zip(axs, ("iii", "i")):
        estilo(ax, f"{NOMEV[v]}: FS no fim de uma chuva de 24 h", "intensidade $I_r$ (mm/h)",
               "$F_{s(min)}$" if v == "iii" else "")
        ax.set_xscale("log")
        ax.axhline(1, color=TINTA2, lw=0.8, ls=(0, (1, 1.5)), zorder=1)
        for s, ks in SOLOS:
            d, cor = dados[(s, v)], CORS[s]
            if v == "iii":
                r = RAH[s]
                ax.plot(Ix, eq6(Ix, r["A"], r["B"], r["C"], r["D"]), color=cor, lw=1.0, ls=(0, (4, 2)), alpha=0.8,
                        zorder=2)
            if d["F"]:
                ax.plot(list(d["F"]), list(d["F"].values()), "o", ms=4, color=cor, mec="white", mew=0.6, zorder=4)
            if d["p"]:
                p = d["p"]
                ax.plot(Ix, eq6(Ix, p["A"], p["B"], p["C"], p["D"]), color=cor, lw=1.4, zorder=3)
                t = d["t"]
                for Ii, mk in ((t["Iini"], "v"), (t["Ithr"], "^")):
                    if XMIN < Ii < XMAX:
                        ax.plot(Ii, eq6(Ii, p["A"], p["B"], p["C"], p["D"]), mk, ms=6, color=cor, mec=TINTA, mew=0.6,
                                zorder=5)
        ax.set_xlim(XMIN, XMAX)
        ax.set_ylim(0, 2.5)
    axs[0].legend(handles=[Line2D([], [], color=TINTA2, marker="o", ls="", ms=4),
                           Line2D([], [], color=TINTA2, lw=1.4), Line2D([], [], color=TINTA2, lw=1.0, ls=(0, (4, 2))),
                           Line2D([], [], color=TINTA2, marker="v", ls="", ms=6),
                           Line2D([], [], color=TINTA2, marker="^", ls="", ms=6)],
                  labels=["nossas contas", "Eq. 6 ajustada", "Rahardjo et al. (2007)", "$I_{r(ini)}$",
                          "$I_{r(thr)}$"], fontsize=6.5, frameon=False, loc="lower left")
    axs[1].legend(handles=[Line2D([], [], color=CORS[s], lw=1.6) for s, _ in SOLOS],
                  labels=[f"{s} ($k_s$ = {ks:g} mm/h)".replace(".", ",") for s, ks in SOLOS], fontsize=6.5,
                  frameon=False, loc="lower left")
    fig.tight_layout()
    salva(fig, "rahardjo_fig6")

    # ---- Fig. 7: um painel por solo e variante, com as constantes e as tangentes
    fig = plt.figure(figsize=(7.4, 8.6))
    gs = GridSpec(3, 2, figure=fig, hspace=0.55, wspace=0.12, left=0.085, right=0.985, top=0.94, bottom=0.06)
    for i, (s, ks) in enumerate(SOLOS):
        for j, v in enumerate(("iii", "i")):
            ax = fig.add_subplot(gs[i, j])
            d, cor = dados[(s, v)], CORS[s]
            estilo(ax, f"{s}, {NOMEV[v]}", "$I_r$ (mm/h)" if i == 2 else "", "$F_{s(min)}$" if j == 0 else "")
            ax.title.set_fontsize(8.5)
            ax.set_xscale("log")
            ax.set_xlim(XMIN, XMAX)
            ax.set_ylim(0, 2.5)
            if j:
                ax.set_yticklabels([])
            ax.axhline(1, color=TINTA2, lw=0.8, ls=(0, (1, 1.5)), zorder=1)
            if v == "iii":
                r = RAH[s]
                ax.plot(Ix, eq6(Ix, r["A"], r["B"], r["C"], r["D"]), color=TINTA2, lw=1.0, ls=(0, (4, 2)), zorder=2)
            if d["F"]:
                ax.plot(list(d["F"]), list(d["F"].values()), "o", ms=4, color=cor, mec="white", mew=0.6, zorder=4)
            if not d["p"]:
                ax.text(0.5, 0.5, "varredura rodando" if not d["F"] else "pontos de menos para o ajuste",
                        transform=ax.transAxes, ha="center", color=TINTA2, fontsize=8)
                continue
            p, t = d["p"], d["t"]
            ax.plot(Ix, eq6(Ix, p["A"], p["B"], p["C"], p["D"]), color=cor, lw=1.4, zorder=3)
            # as tangentes: a inclinada e os dois patamares
            # a inclinada so entre os patamares, com 0,15 de folga em FS
            xa, xb = (t["xk"] + (F_ - t["Fk"]) / t["dF"] for F_ in (t["Fup"] + 0.15, t["Flo"] - 0.15))
            xs = np.linspace(min(xa, xb), max(xa, xb), 50)
            ax.plot(10 ** xs, t["Fk"] + t["dF"] * (xs - t["xk"]), color=TINTA, lw=0.7, zorder=2)
            ax.plot([XMIN, t["Iini"] * 2], [t["Fup"]] * 2, color=TINTA, lw=0.7, zorder=2)
            ax.plot([t["Ithr"] / 2, XMAX], [t["Flo"]] * 2, color=TINTA, lw=0.7, zorder=2)
            txt = []
            for nome, Ii, mk in (("ini", t["Iini"], "v"), ("thr", t["Ithr"], "^")):
                ax.plot(Ii, eq6(Ii, p["A"], p["B"], p["C"], p["D"]), mk, ms=6, color=cor, mec=TINTA, mew=0.6, zorder=5)
                txt.append(f"$I_{{r({nome})}}$ = {Ii:.3g} mm/h = {Ii / ks:.2g} $k_s$; $V_r$ = {24 * Ii:.0f} mm")
            txt.append(f"A = {p['A']:.3g}; B = {p['B']:.3g}; C = {p['C']:.3g}; D = {p['D']:.3g}; $r^2$ = {p['r2']:.3f}")
            ax.text(0.02, 0.03, "\n".join(txt).replace(".", ","), transform=ax.transAxes, fontsize=5.8, color=TINTA,
                    va="bottom", ha="left", linespacing=1.35)
    fig.legend(handles=[Line2D([], [], color=TINTA2, marker="o", ls=""), Line2D([], [], color=TINTA2, lw=1.4),
                        Line2D([], [], color=TINTA, lw=0.7), Line2D([], [], color=TINTA2, lw=1.0, ls=(0, (4, 2)))],
               labels=["nossas contas", "Eq. 6 ajustada", "tangentes", "Rahardjo et al. (2007)"],
               loc="upper center", ncol=4, fontsize=6.5, frameon=False, bbox_to_anchor=(0.5, 0.998))
    salva(fig, "rahardjo_fig7")

    # ---- a tabela das constantes
    L = []
    for s, ks in SOLOS:
        r, rt = RAH[s], rah_t[s]
        Icr = I_de_F(r, 1.)
        L.append(f"{s} & {f'{ks:g}'.replace('.', '{,}')} & Rahardjo et al. & {num(r['A'])} & {num(r['B'])} & "
                 f"{num(r['C'], 2)} & {num(r['D'], 2)} & {num(r['r2'], 3)} & -- & {num(r['Iini'], 1)} "
                 f"({num(r['Iini'] / ks, 2)}) & {num(r['Ithr'], 0)} ({num(r['Ithr'] / ks, 2)}) & "
                 f"{num(Icr, 1) if Icr else '--'}\\\\")
        L.append(f" & & idem, nossas tangentes & & & & & & & {num(rt['Iini'], 1)} ({num(rt['Iini'] / ks, 2)}) & "
                 f"{num(rt['Ithr'], 0)} ({num(rt['Ithr'] / ks, 2)}) & --\\\\")
        for v in ("iii", "i"):
            d = dados[(s, v)]
            fs0 = num(d["fs0"], 3) if d["fs0"] else "--"
            if d["p"]:
                p, t = d["p"], d["t"]
                Ic = num(d["Ic"], 1) if d["Ic"] else "--"
                L.append(f" & & nosso, {NOMEV[v]} & {num(p['A'])} & {num(p['B'])} & {num(p['C'], 2)} & {num(p['D'], 2)} & "
                         f"{num(p['r2'], 3)} & {fs0} & {num(t['Iini'], 1)} ({num(t['Iini'] / ks, 2)}) & "
                         f"{num(t['Ithr'], 0)} ({num(t['Ithr'] / ks, 2)}) & {Ic}\\\\")
            else:
                L.append(f" & & nosso, {NOMEV[v]} & \\multicolumn{{9}}{{l}}{{"
                         f"{'varredura rodando' if not d['F'] else f'{len(d[chr(70)])} pontos: de menos para o ajuste'}"
                         f"}}\\\\")
        L.append("\\addlinespace")
    with open(os.path.join(TABS, "rahardjo_const.tex"), "w") as f:
        f.write("\\begin{tabular}{lrlrrrrrrlll}\n\\toprule\nsolo & $k_s$ (mm/h) & fonte & $A$ & $B$ & $C$ & $D$ & $r^2$ & "
                "FS(0) & $I_{r(ini)}$ ($\\times k_s$) & $I_{r(thr)}$ ($\\times k_s$) & $I_c$\\\\\n\\midrule\n"
                + "\n".join(L[:-1]) + "\n\\bottomrule\n\\end{tabular}\n")
    # as macros do texto: as tangentes recalculadas sobre as constantes do artigo, no f50
    M["rahTangIniCinq"], M["rahTangThrCinq"] = num(rah_t["f50"]["Iini"], 1), num(rah_t["f50"]["Ithr"], 0)
    M["rahSolosProntos"] = str(sum(1 for (s, v), d in dados.items() if d["p"] and v == "iii"))
    return M
