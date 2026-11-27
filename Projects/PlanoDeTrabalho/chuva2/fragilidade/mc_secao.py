#!/usr/bin/env python3
"""mc_secao.py: a secao 'Simulacao de Monte Carlo' do relatorio da fragilidade, chamada pelo gera_dados.py (gera(...)):
  - os casos do PC (~/projects/runs/pc_f50_ls_mc; f50, chuva constante de 24 h, sempre ate CoV(P_f) = 5 %): o MC direto
    pela primeira passagem (mc/saida/<caso>/mc_crpc_*.csv: a realizacao para em t_1) e o line sampling na direcao do FORM
    modal (ls/<caso>/ls/linhas.csv e custos.csv; ls/<caso>/form/form.json);
  - a campanha frag do fermi (~/projects/runs/fermi_frag, parcial), com as funcoes do cluster/resumo_frag.py (o prefixo
    continuo de blocos prontos e a monotonia em I).
Cada caso: a convergencia da P_f (a estimativa corrente contra N, com o IC de 95 %), o CoV(P_f) contra N e o histograma
de t_1 das falhas. As estimativas correntes usam as realizacoes (ou linhas) na ordem do indice, so o prefixo continuo das
ja prontas (sem o vies das que terminam primeiro). Saidas: figs/mc_pc.pdf|png, figs/frag_<v>.pdf|png, tab/mc_pc.tex,
tab/frag_<v>.tex; devolve as macros do numeros.tex."""
import csv
import glob
import json
import math
import os
import sys
import time

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.legend_handler import HandlerTuple
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import MaxNLocator
from scipy.optimize import minimize
from scipy.stats import norm

RUNS = os.path.expanduser("~/projects/runs")
PC = os.path.join(RUNS, "pc_f50_ls_mc")
FRAG = os.path.join(RUNS, "fermi_frag")
REF = os.path.join(RUNS, "fermi_chuva2", "chuva2", "ref")
S50 = os.path.join(RUNS, "sigmoide_f50")
KEF = os.path.join(RUNS, "sigmoide_f50_kef")
sys.path.insert(0, os.path.expanduser("~/projects/neopz-master/Projects/SlopeRandomFieldRefinement/cluster"))
import resumo_frag as rf  # noqa: E402

D = 24.0  # a duracao da chuva (h)
Z95 = 1.959964
# os casos do PC: (rotulo, variante, I, metodo)
CASOS_PC = [("iii_I30", "iii", 30., "ls"), ("iii_I36", "iii", 36., "ls"), ("i_I30", "i", 30., "mc"),
            ("i_I36", "i", 36., "mc")]
NOMEV = {"iii": "(iii)", "i": "(i)"}


def num(x, nd=2):
    return f"{x:.{nd}f}".replace(".", "{,}")


def pcs(p, nd=1):
    """P_f em %, sem o simbolo"""
    if p is None or (isinstance(p, float) and math.isnan(p)):
        return "--"
    return num(100 * p, nd)


def wilson(k, n):
    if n == 0:
        return float("nan"), float("nan")
    p = k / n
    c = (p + Z95 ** 2 / (2 * n)) / (1 + Z95 ** 2 / n)
    h = Z95 * math.sqrt(p * (1 - p) / n + Z95 ** 2 / (4 * n * n)) / (1 + Z95 ** 2 / n)
    return max(0., c - h), min(1., c + h)


# ------------------------------------------------------------------------------------------------- estimativas correntes
def corrida_mc(f):
    """MC: N, P_f, IC de Wilson e CoV = sqrt((1 - P)/(N P)) correntes (as realizacoes na ordem)"""
    f = np.asarray(f, float)
    n = np.arange(1, len(f) + 1, dtype=float)
    k = np.cumsum(f)
    p = k / n
    c = (p + Z95 ** 2 / (2 * n)) / (1 + Z95 ** 2 / n)
    h = Z95 * np.sqrt(p * (1 - p) / n + Z95 ** 2 / (4 * n * n)) / (1 + Z95 ** 2 / n)
    with np.errstate(divide="ignore", invalid="ignore"):
        cv = np.where(k > 0, np.sqrt((1 - p) / np.maximum(k, 1)), np.nan)
    return n, p, np.maximum(c - h, 0), np.minimum(c + h, 1), cv


def corrida_ls(p):
    """LS: N (linhas), P_f = media das p_i, IC normal de 95 % e CoV = dp / (sqrt(N) P_f) correntes"""
    p = np.asarray(p, float)
    n = np.arange(1, len(p) + 1, dtype=float)
    m = np.cumsum(p) / n
    var = np.where(n > 1, (np.cumsum(p * p) - n * m * m) / np.maximum(n - 1, 1), np.nan)
    se = np.sqrt(np.maximum(var, 0)) / np.sqrt(n)
    with np.errstate(divide="ignore", invalid="ignore"):
        cv = np.where(m > 0, se / m, np.nan)
    return n, m, np.maximum(m - Z95 * se, 0), m + Z95 * se, cv


# ----------------------------------------------------------------------------------------------------------- leituras
def le_mc_pc(caso):
    """o MC do PC (uma intensidade): o prefixo continuo das realizacoes prontas, em ordem"""
    d = os.path.join(PC, "mc", "saida", caso)
    r = {}
    for a in glob.glob(os.path.join(d, "mc_crpc_*.csv")):
        for x in csv.DictReader(open(a)):
            n = int(x["real"])
            if n >= 0 and x["falha"] in ("sim", "nao"):
                r[n] = (x["falha"] == "sim", float(x["DC_h"]), (float(x["s_transiente"]) + float(x["s_fator"])) / 60.)
    k = 0
    while k in r:
        k += 1
    est = "rodando"
    at = os.path.join(d, "ATINGIU")
    if os.path.exists(at):
        est = "parado" if open(at).read().startswith("PARADO") else "CoV 5\\,\\%"
    return {"f": np.array([r[n][0] for n in range(k)], bool), "ids": list(range(k)),
            "t1": np.array([r[n][1] for n in range(k) if r[n][0]]), "custo": np.array([r[n][2] for n in range(k)]),
            "est": est, "todas": r}


def le_ls_pc(caso):
    """o line sampling do PC: o prefixo continuo das linhas prontas (a linha i usa sempre o mesmo Z[i]), o custo de cada
    linha (min, a soma das analises do custos.csv) e o FORM"""
    d = os.path.join(PC, "ls", caso)
    p = {}
    arq = os.path.join(d, "ls", "linhas.csv")
    if os.path.exists(arq):
        for x in csv.DictReader(open(arq)):
            p[int(x["linha"])] = (float(x["p"]), int(x["analises"]))
    k = 0
    while k in p:
        k += 1
    custo = {}
    ac = os.path.join(d, "ls", "custos.csv")
    if os.path.exists(ac):
        for x in csv.DictReader(open(ac)):
            custo[int(x["linha"])] = custo.get(int(x["linha"]), 0.) + float(x["s"]) / 60.
    fj = os.path.join(d, "form", "form.json")
    form = json.load(open(fj)) if os.path.exists(fj) else {}
    lg = os.path.join(PC, "ls", caso + ".log")
    txt = open(lg).read() if os.path.exists(lg) else ""
    fim = txt.rfind(" fim") > txt.rfind("inicio LS") >= 0
    est = "rodando"
    if fim:
        lsj = os.path.join(d, "ls", "ls.json")
        mot = str(json.load(open(lsj)).get("motivo", "")) if os.path.exists(lsj) else ""
        est = "CoV 5\\,\\%" if "CoV" in mot else "terminou"
    return {"p": np.array([p[i][0] for i in range(k)]), "an": np.array([p[i][1] for i in range(k)]),
            "custo": np.array([custo[i] for i in range(k) if i in custo]), "beta_form": form.get("beta"),
            "pf_form": form.get("P_f"), "est": est}


def F_em(F, I):
    """F(I): o ponto deterministico ou, entre dois pontos, ln F interpolado em I (fora da faixa: None)"""
    if not F:
        return None
    if I in F:
        return F[I]
    return math.exp(rf.interp_log(F, I)) if min(F) < I < max(F) else None


def calibracao(var):
    """beta e a mediana do FS no fim da chuva do MC de calibracao da variante, e os F das medias (k_s medio e mediano)"""
    if var == "iii":
        fcal, Ical = rf.le_mc(REF, D), 31.05
    else:
        fcal, Ical = rf.le_mc(os.path.join(FRAG, "f50_i_mc20"), D), 20.
    lns = np.log([x for x in fcal.values() if x > 0])
    if len(lns) < 20:
        return None
    beta, med = float(np.std(lns, ddof=1)), float(np.exp(lns.mean()))
    Fm = rf.le_det(S50, "" if var == "iii" else "i_", D)
    Fk = rf.le_det(KEF, ("" if var == "iii" else "i_") + "ks8.57e-6_", D)
    Fmc, Fkc = F_em(Fm, Ical), F_em(Fk, Ical)
    return {"fcal": fcal, "Ical": Ical, "beta": beta, "med": med, "n": len(lns), "Fm": Fm, "Fk": Fk,
            "m": med / Fmc if Fmc else None, "mk": med / Fkc if Fkc else None}


def prev(cal, I, kmed):
    """a formula Phi(-ln(m F)/beta), com o F do k_s medio (kmed=False) ou mediano (True) e o m da calibracao"""
    if cal is None:
        return None
    F, m = (cal["Fk"], cal["mk"]) if kmed else (cal["Fm"], cal["m"])
    FI = F_em(F, I)
    if m is None or FI is None:
        return None
    return float(norm.cdf(-math.log(m * FI) / cal["beta"]))


def le_frag(var):
    """a primeira passagem do frag (f50_<var>): as realizacoes validas do prefixo, em ordem, e por intensidade a falha
    ate D (sem rodar = sobreviveu numa maior) e os t_1"""
    pasta = os.path.join(FRAG, "f50_" + var)
    if not os.path.isdir(pasta):
        return None
    dados = rf.le_crpc(pasta)
    val = sorted(n for n, r in dados.items() if all(f in ("sim", "nao") for f, _ in r.values()))
    if not val:
        return None
    Is = sorted({I for n in val for I in dados[n]}, reverse=True)
    res = {}
    for I in Is:
        fal = [I in dados[n] and dados[n][I][0] == "sim" and dados[n][I][1] <= D + 1e-9 for n in val]
        res[I] = {"f": np.array(fal, bool), "t1": np.array([dados[n][I][1] for n, x in zip(val, fal) if x])}
    return {"val": val, "Is": Is, "res": res, "pasta": pasta, "dados": dados}


def todas_frag(var):
    """todas as linhas prontas do crpc do frag (realizacao x intensidade), com ou sem o bloco terminado: para a
    conferencia realizacao a realizacao, que nao sofre do vies de quem termina primeiro"""
    r = {}
    for a in glob.glob(os.path.join(FRAG, "f50_" + var, "mc_crpc_*.csv")):
        for x in csv.DictReader(open(a)):
            if int(x["real"]) >= 0 and x["falha"] in ("sim", "nao"):
                r.setdefault(int(x["real"]), {})[float(x["I"])] = (x["falha"], float(x["DC_h"]))
    return r


def fs_sobreviventes(fr, I):
    """o FS de 24 h das realizacoes do prefixo que sobreviveram a I (as linhas 'nao' do crpc com t_lo = D)"""
    s, fs = set(fr["val"]), []
    for a in glob.glob(os.path.join(fr["pasta"], "mc_crpc_*.csv")):
        if not rf.bloco_pronto(a, "crpc"):
            continue
        for x in csv.DictReader(open(a)):
            if (int(x["real"]) in s and float(x["I"]) == I and x["falha"] == "nao"
                    and abs(float(x["t_lo"]) - D) < 1e-9):
                fs.append(float(x["FS_lo"]))
    return np.array(fs)


def ajuste_censurado(fs, nf):
    """lognormal do FS por maxima verossimilhanca com censura a esquerda em 1 (das falhas so se sabe FS < 1): mediana,
    beta e os erros padrao (hessiana numerica)"""
    if len(fs) < 10 or nf < 5:
        return None
    x = np.log(fs)

    def nll(t):
        b = math.exp(t[1])
        return -(norm.logpdf(x, t[0], b).sum() + nf * norm.logcdf(-t[0] / b))

    o = minimize(nll, [float(np.median(x)), math.log(0.08)], method="Nelder-Mead",
                 options={"xatol": 1e-8, "fatol": 1e-10, "maxiter": 5000})
    h, H = 1e-4, np.zeros((2, 2))
    for i in range(2):
        for j in range(2):
            ei, ej = np.eye(2)[i] * h, np.eye(2)[j] * h
            H[i, j] = (nll(o.x + ei + ej) - nll(o.x + ei - ej) - nll(o.x - ei + ej) + nll(o.x - ei - ej)) / (4 * h * h)
    cov = np.linalg.inv(H)
    mu, b = float(o.x[0]), math.exp(float(o.x[1]))
    return {"med": math.exp(mu), "med_se": math.exp(mu) * math.sqrt(cov[0, 0]), "beta": b,
            "beta_se": b * math.sqrt(cov[1, 1]), "pf": float(norm.cdf(-mu / b))}


# ----------------------------------------------------------------------------------------------------------- figuras
def faixa_t1(listas):
    """o eixo de t_1 comum a uma figura: de pouco antes da primeira falha (multiplo de 2 h) ate D"""
    tudo = [x for t in listas for x in t]
    return (max(0., 2 * math.floor((min(tudo) - 1) / 2)) if tudo else 0.), D


def hist_t1(ax, t1, cor, xr, TINTA, TINTA2, rotulo=True):
    ax.hist(t1, bins=np.arange(math.floor(xr[0]), D + 1, 1), color=cor, rwidth=0.85, zorder=3)
    md = float(np.median(t1))
    ax.axvline(md, color=TINTA, lw=0.9, ls=(0, (4, 2)), zorder=4)
    if rotulo:
        ax.text(0.03, 0.97, f"mediana\n{md:.1f} h".replace(".", ","), transform=ax.transAxes, fontsize=6.5,
                color=TINTA2, va="top", ha="left")
    ax.set_xlim(*xr)
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))


def rotulos_finais(ax, pontos, cor_txt, dy_min):
    """rotulos no fim das curvas, afastados na vertical para nao se sobreporem"""
    pts = sorted(pontos, key=lambda t: t[1])
    ys = []
    for _, y, _ in pts:
        ys.append(max(y, ys[-1] + dy_min) if ys else y)
    for (x, _, txt), y in zip(pts, ys):
        ax.annotate(txt, (x, y), xytext=(4, 0), textcoords="offset points", va="center", fontsize=6.5, color=cor_txt)


def fig_pc(casos, k):
    """4 linhas (casos) x 3 colunas: a convergencia da P_f, o CoV e o histograma de t_1"""
    estilo, salva, COR, TINTA, TINTA2, FAIXA = k["estilo"], k["salva"], k["COR"], k["TINTA"], k["TINTA2"], k["FAIXA"]
    fig = plt.figure(figsize=(7.4, 2.35 * len(casos) + 0.5))
    gs = GridSpec(len(casos), 3, figure=fig, width_ratios=[1.3, 1., 1.], hspace=0.75, wspace=0.36, left=0.075,
                  right=0.985, top=0.935, bottom=0.05)
    xr = faixa_t1([c["t1"] for c in casos])
    fig.legend(handles=[Line2D([], [], color=COR["iii"], lw=1.4), Line2D([], [], color=COR["i"], lw=1.4),
                        Patch(facecolor=COR["iii"], alpha=0.18, lw=0),
                        Line2D([], [], color=TINTA, lw=0.9, ls=(0, (4, 2))),
                        Line2D([], [], color=TINTA2, lw=0.9, ls=(0, (1, 1.5))),
                        (Patch(facecolor=FAIXA, lw=0), Line2D([], [], color=TINTA2, lw=0.9))],
               labels=["estimativa (iii)", "estimativa (i)", "IC 95%", "fórmula, $k_s$ mediano", "FORM",
                       "MC do fermi (IC 95%)"],
               handler_map={tuple: HandlerTuple(ndivide=None, pad=0)}, loc="upper center", ncol=6, fontsize=6.5,
               frameon=False, bbox_to_anchor=(0.5, 0.995), handlelength=2.4, columnspacing=1.2)
    for i, c in enumerate(casos):
        cor = COR[c["var"]]
        ax0, ax1, ax2 = (fig.add_subplot(gs[i, j]) for j in range(3))
        unid = "linhas" if c["metodo"] == "ls" else "realizações"
        estilo(ax0, f"{NOMEV[c['var']]} a {c['I']:g} mm/h: {c['metodo_curto']}", f"N ({unid})", "$P_f$ (%)")
        estilo(ax1, "CoV($P_f$)", f"N ({unid})", "CoV (%)")
        estilo(ax2, c["t1_titulo"], "$t_1$ (h)", "falhas")
        for ax in (ax0, ax1, ax2):
            ax.title.set_fontsize(8)
        n, p, lo, hi, cv = c["corrida"]
        if len(n) == 0:
            ax0.text(0.5, 0.5, "sem dados ainda", transform=ax0.transAxes, ha="center", color=TINTA2)
        else:
            if "mc_ref" in c:  # o MC das mesmas condicoes no fermi: a faixa do IC
                pr, a_, b_, rot = c["mc_ref"]
                ax0.axhspan(100 * a_, 100 * b_, color=FAIXA, lw=0, zorder=1)
                ax0.axhline(100 * pr, color=TINTA2, lw=0.9, zorder=2)
            ax0.fill_between(n, 100 * lo, 100 * hi, color=cor, alpha=0.18, lw=0, zorder=2)
            ax0.plot(n, 100 * p, color=cor, lw=1.4, zorder=3)
            ax0.plot(n[-1], 100 * p[-1], "o", ms=3.5, color=cor, zorder=4)
            if c.get("formula") is not None:
                ax0.axhline(100 * c["formula"], color=TINTA, lw=0.9, ls=(0, (4, 2)), zorder=2)
            if c.get("pf_form") is not None:
                ax0.axhline(100 * c["pf_form"], color=TINTA2, lw=0.9, ls=(0, (1, 1.5)), zorder=2)
            ax0.set_xscale("log")
            ax0.set_xlim(1, max(n[-1], 10) * 1.05)
            ok = np.isfinite(cv) & (cv > 0) & (n >= (5 if c["metodo"] == "ls" else 1))
            if ok.any():
                ax1.plot(n[ok], 100 * cv[ok], color=cor, lw=1.4, zorder=3)
                ax1.plot(n[ok][-1], 100 * cv[ok][-1], "o", ms=3.5, color=cor, zorder=4)
            ax1.axhline(5, color=TINTA, lw=0.9, ls=(0, (4, 2)), zorder=2)
            ax1.text(1.15, 5 * 1.12, "alvo 5%", fontsize=6.5, color=TINTA2)
            ax1.set_xscale("log")
            ax1.set_yscale("log")
            ax1.set_xlim(1, max(n[-1], 10) * 1.05)
            ax1.set_ylim(2, max(200, float(np.nanmax(100 * cv[ok])) * 1.3) if ok.any() else 200)
        t1 = c["t1"]
        if len(t1):
            hist_t1(ax2, t1, cor, xr, TINTA, TINTA2)
        else:
            ax2.text(0.5, 0.5, "sem falhas ainda", transform=ax2.transAxes, ha="center", color=TINTA2)
    salva(fig, "mc_pc")


def fig_frag(var, fr, k):
    """a campanha frag de uma variante: a convergencia da P_f e o CoV de cada intensidade (rampa da cor da variante:
    clara = I baixa) e os histogramas de t_1 das intensidades com falhas"""
    estilo, salva, TINTA, TINTA2 = k["estilo"], k["salva"], k["TINTA"], k["TINTA2"]
    Is = sorted(fr["Is"])
    mapa = plt.cm.Blues if var == "iii" else plt.cm.Greens
    cores = {I: mapa(x) for I, x in zip(Is, np.linspace(0.42, 0.95, len(Is)))}
    comfal = [I for I in sorted(Is, reverse=True) if len(fr["res"][I]["t1"])][:6]
    fig = plt.figure(figsize=(7.4, 5.6))
    gs = GridSpec(2, 6, figure=fig, height_ratios=[1.4, 1.], hspace=0.62, wspace=0.75, left=0.075, right=0.975,
                  top=0.94, bottom=0.08)
    axA, axB = fig.add_subplot(gs[0, :3]), fig.add_subplot(gs[0, 3:])
    estilo(axA, f"{NOMEV[var]}: $P_f$ corrente (rótulo: I em mm/h)", "N (realizações)", "$P_f$ (%)")
    estilo(axB, "CoV($P_f$)", "N (realizações)", "CoV (%)")
    fins, nmax = [], len(fr["val"])
    for I in Is:
        n, p, lo, hi, cv = corrida_mc(fr["res"][I]["f"])
        axA.plot(n, 100 * p, color=cores[I], lw=1.3)
        fins.append((n[-1], 100 * p[-1], f"{I:g}".replace(".", ",")))
        ok = np.isfinite(cv) & (cv > 0)
        if ok.any():
            axB.plot(n[ok], 100 * cv[ok], color=cores[I], lw=1.3)
    rotulos_finais(axA, fins, TINTA2, (axA.get_ylim()[1] - axA.get_ylim()[0]) * 0.045)
    axA.set_xlim(0, nmax * 1.12)
    axB.axhline(5, color=TINTA, lw=0.9, ls=(0, (4, 2)))
    axB.text(1.15, 5 * 1.12, "alvo 5%", fontsize=6.5, color=TINTA2)
    axB.set_xscale("log")
    axB.set_yscale("log")
    axB.set_xlim(1, nmax * 1.05)
    xr = faixa_t1([fr["res"][I]["t1"] for I in comfal])
    for j, I in enumerate(comfal):
        ax = fig.add_subplot(gs[1, j])
        t1 = fr["res"][I]["t1"]
        estilo(ax, f"{I:g} mm/h: {len(t1)}".replace(".", ","), "$t_1$ (h)", "falhas" if j == 0 else "")
        ax.title.set_fontsize(7.5)
        hist_t1(ax, t1, cores[I], xr, TINTA, TINTA2, rotulo=False)
        ax.set_xticks([x for x in (0, 6, 12, 18, 24) if x >= xr[0]])
        ax.tick_params(labelsize=6.5)
    salva(fig, f"frag_{var}")


# ----------------------------------------------------------------------------------------------------------- tabelas
def linha_mc(rot, metodo, f, t1, custo, formula, est):
    k_, n_ = int(np.sum(f)), len(f)
    if n_ == 0:
        return f"{rot} & {metodo} & 0 & -- & -- & -- & -- & -- & -- & {pcs(formula)} & {est}\\\\"
    p = k_ / n_
    lo, hi = wilson(k_, n_)
    cv = math.sqrt((1 - p) / k_) if k_ else float("nan")
    tq = (f"{num(np.median(t1), 1)} [{num(np.quantile(t1, .1), 1)}, {num(np.quantile(t1, .9), 1)}]"
          if len(t1) else "--")
    return (f"{rot} & {metodo} & {n_} & {k_} & {pcs(p)} & {num(100 * cv, 1) if k_ else '--'} & {pcs(lo)}--{pcs(hi)} & "
            f"{tq} & {num(float(np.mean(custo)), 1) if len(custo) else '--'} & {pcs(formula)} & {est}\\\\")


def linha_ls(rot, metodo, L, formula):
    n_ = len(L["p"])
    if n_ < 2:
        return f"{rot} & {metodo} & {n_} & -- & -- & -- & -- & -- & -- & {pcs(formula)} & {L['est']}\\\\"
    n, m, lo, hi, cv = corrida_ls(L["p"])
    cvt = num(100 * cv[-1], 1) if n_ >= 10 else "--"  # com poucas linhas o desvio estimado nao vale nada
    return (f"{rot} & {metodo} & {n_} & -- & {pcs(m[-1])} & {cvt} & {pcs(lo[-1])}--{pcs(hi[-1])} & "
            f"-- & {num(float(np.mean(L['custo'])), 1) if len(L['custo']) else '--'} & {pcs(formula)} & "
            f"{L['est']}\\\\")


# ----------------------------------------------------------------------------------------------------------- principal
def gera(k):
    """k: os utilitarios do gera_dados.py (estilo, salva, COR, TINTA, TINTA2, FAIXA, TABS). Devolve as macros."""
    TABS = k["TABS"]
    M = {"mcData": time.strftime("%d/%m/%Y, %H:%M")}
    cal = {v: calibracao(v) for v in ("iii", "i")}
    frag = {v: le_frag(v) for v in ("iii", "i")}

    # ---- os casos do PC
    casos, linhas = [], []
    for rot, var, I, met in CASOS_PC:
        form_k = prev(cal[var], I, True)
        fr = frag[var]
        fr_I = fr["res"].get(I) if fr else None
        nome = f"{NOMEV[var]} {I:g}~mm/h"
        if met == "ls":
            L = le_ls_pc(rot)
            c = {"var": var, "I": I, "metodo": "ls", "metodo_curto": "line sampling", "corrida": corrida_ls(L["p"]),
                 "formula": form_k, "pf_form": L["pf_form"]}
            linhas.append(linha_ls(nome, "LS", L, form_k))
            sufixo = "A" if I == 30. else "B"  # as macros do texto: A = 30 mm/h, B = 36 mm/h
            if L["beta_form"] is not None:
                M["mcFormBeta" + sufixo], M["mcFormPf" + sufixo] = num(L["beta_form"], 2), pcs(L["pf_form"])
            # o t_1 e a faixa do MC vem do frag (o LS nao da t_1)
            if fr_I is not None:
                kf, nf_ = int(fr_I["f"].sum()), len(fr_I["f"])
                lo, hi = wilson(kf, nf_)
                c["mc_ref"] = (kf / nf_, lo, hi, f"MC fermi, N = {nf_} ({100 * kf / nf_:.1f}%)".replace(".", ","))
                c["t1"], c["t1_titulo"] = fr_I["t1"], f"$t_1$: MC fermi (n = {len(fr_I['t1'])})"
                M["mcFragPf" + sufixo] = pcs(kf / nf_)
                linhas.append(linha_mc("", "MC (fermi)", fr_I["f"], fr_I["t1"], [], form_k, "parcial"))
            else:
                c["t1"], c["t1_titulo"] = np.array([]), "$t_1$"
            if rot == "iii_I36":  # o MC do PC antes de passar para o LS (as mesmas realizacoes do frag)
                Q = le_mc_pc(rot)
                if len(Q["f"]):
                    linhas.append(linha_mc("", "MC (PC)", Q["f"], Q["t1"], Q["custo"], form_k, Q["est"]))
        else:
            Q = le_mc_pc(rot)
            c = {"var": var, "I": I, "metodo": "mc", "metodo_curto": "MC até $t_1$", "corrida": corrida_mc(Q["f"]),
                 "formula": form_k, "t1": Q["t1"], "t1_titulo": f"$t_1$ das falhas (n = {len(Q['t1'])})"}
            linhas.append(linha_mc(nome, "MC (PC)", Q["f"], Q["t1"], Q["custo"], form_k, Q["est"]))
        casos.append(c)
    fig_pc(casos, k)
    with open(os.path.join(TABS, "mc_pc.tex"), "w") as f:
        f.write("\\begin{tabular}{llrrrrlllrl}\n\\toprule\ncaso & método & $N$ & falhas & $\\Pf$ (\\%) & CoV (\\%) & "
                "IC 95\\,\\% & $t_1$ (h) & min/$N$ & fórmula (\\%) & estado\\\\\n"
                "\\midrule\n" + "\n".join(linhas) + "\n\\bottomrule\n\\end{tabular}\n")

    # ---- conferencia PC x fermi nas mesmas realizacoes (o frag tem a (iii) a 36 e a (i) a 36; a (i) a 30 pela monotonia)
    conf = []
    for rot, var, I, met in CASOS_PC:
        Q = le_mc_pc(rot) if os.path.isdir(os.path.join(PC, "mc", "saida", rot)) else None
        if not Q or not Q["todas"]:
            continue
        dados, comuns, dif = todas_frag(var), 0, 0
        for n, (fp, _, _) in Q["todas"].items():
            if n not in dados:
                continue
            r = dados[n]
            if I in r:
                ff = r[I][0] == "sim" and r[I][1] <= D + 1e-9
            elif any(r[i][0] == "nao" for i in r if i > I):
                ff = False
            elif any(r[i][0] == "sim" for i in r if i < I):
                ff = True
            else:
                continue
            comuns += 1
            dif += fp != ff
        if comuns:
            conf.append(f"{dif} em {comuns} na {NOMEV[var]} a {I:g}~mm/h")
    M["mcConf"] = "; ".join(conf) if conf else "sem realizações em comum ainda"

    # ---- a campanha frag
    for var in ("iii", "i"):
        fr, cv_ = frag[var], cal[var]
        M["fragTem" + var] = "1" if fr else "0"
        if not fr:
            M["fragN" + var] = "0"
            continue
        M["fragN" + var] = str(len(fr["val"]))
        fig_frag(var, fr, k)
        L, zk, zm, par = [], [], [], []
        for I in sorted(fr["Is"], reverse=True):
            f, t1 = fr["res"][I]["f"], fr["res"][I]["t1"]
            n_, k_ = len(f), int(f.sum())
            p = k_ / n_
            lo, hi = wilson(k_, n_)
            pm, pk = prev(cv_, I, False), prev(cv_, I, True)
            z = (lambda q: (p - q) / math.sqrt(q * (1 - q) / n_) if q is not None and 0 < q < 1 else None)
            zm_, zk_ = z(pm), z(pk)
            if zm_ is not None:
                zm.append((I, zm_))
            if zk_ is not None:
                zk.append((I, zk_))
            # o teste pareado: FS_k(I) = FS_k(I_cal) F(I) / F(I_cal) (F com o k_s mediano)
            dv = None
            if cv_ and F_em(cv_["Fk"], I) and F_em(cv_["Fk"], cv_["Ical"]):
                dv, r_ = 0, F_em(cv_["Fk"], I) / F_em(cv_["Fk"], cv_["Ical"])
                for n, x in zip(fr["val"], f):
                    if n in cv_["fcal"]:
                        dv += (cv_["fcal"][n] * r_ < 1) != bool(x)
                par.append((I, dv))
            Fm, Fk = (F_em(cv_["Fm"], I), F_em(cv_["Fk"], I)) if cv_ else (None, None)
            tq = (f"{num(np.median(t1), 1)} [{num(np.quantile(t1, .1), 1)}, {num(np.quantile(t1, .9), 1)}]"
                  if len(t1) else "--")
            fz = (lambda q, zz: f"{pcs(q)} ({'--' if zz is None else num(zz, 1)})" if q is not None else "--")
            L.append(f"{num(I, 2).replace('{,}00', '')} & {num(Fm, 3) if Fm else '--'} & {num(Fk, 3) if Fk else '--'} & "
                     f"{k_}/{n_} & {pcs(p)} & {pcs(lo)}--{pcs(hi)} & {num(100 * math.sqrt((1 - p) / k_), 0) if k_ else '--'}"
                     f" & {fz(pm, zm_)} & {fz(pk, zk_)} & {'--' if dv is None else dv} & {tq}\\\\")
        with open(os.path.join(TABS, f"frag_{var}.tex"), "w") as fo:
            fo.write("\\begin{tabular}{rrrrrlrllrl}\n\\toprule\n$I$ & $F$ médio & $F$ mediano & falhas/$N$ & "
                     "$\\Pf$ (\\%) & IC 95\\,\\% & CoV (\\%) & prev. médio ($z$) & prev. mediano ($z$) & pareado & "
                     "$t_1$ (h)\\\\\n"
                     "\\midrule\n" + "\n".join(L) + "\n\\bottomrule\n\\end{tabular}\n")
        if zk:
            M["fragZk" + var] = num(max(abs(z) for _, z in zk), 1)
            M["fragZkDentro" + var] = f"{sum(abs(z) <= 2 for _, z in zk)} de {len(zk)}"
        if zm:
            M["fragZm" + var] = num(min(z for _, z in zm), 1)
            M["fragZmMax" + var] = num(max(abs(z) for _, z in zm), 1)
            M["fragZmDentro" + var] = f"{sum(abs(z) <= 2 for _, z in zm)} de {len(zm)}"
        if par:
            M["fragPar" + var] = "; ".join(f"{num(I, 2).replace('{,}00', '')}: {d}" for I, d in sorted(par))
        if cv_:
            M["fragBeta" + var] = num(cv_["beta"], 3)
            M["fragMk" + var] = num(cv_["mk"], 3) if cv_["mk"] else "--"
            M["fragIcal" + var] = num(cv_["Ical"], 2).replace("{,}00", "")
        # a lognormal censurada na maior intensidade (todas as realizacoes a rodam)
        Itop = max(fr["Is"])
        aj = ajuste_censurado(fs_sobreviventes(fr, Itop), int(fr["res"][Itop]["f"].sum()))
        if aj and cv_ and cv_["mk"] and F_em(cv_["Fk"], Itop):
            M.update({"fragTopI" + var: num(Itop, 0), "fragCensMed" + var: num(aj["med"], 3),
                      "fragCensMedSe" + var: num(aj["med_se"], 3), "fragCensBeta" + var: num(aj["beta"], 3),
                      "fragCensBetaSe" + var: num(aj["beta_se"], 3),
                      "fragPrevMed" + var: num(cv_["mk"] * F_em(cv_["Fk"], Itop), 3)})
            if cv_["m"] and F_em(cv_["Fm"], Itop):
                M["fragPrevMedMedio" + var] = num(cv_["m"] * F_em(cv_["Fm"], Itop), 3)
    M.update(f100(TABS))
    return M


def f100(TABS):
    """o f100 do frag (parcial): chuva de ate 120 h, falha contada ate 72 h (a duracao das sigmoides do f100); beta e m
    do MC de calibracao (FS em 72 h a 9 mm/h na (iii) e a 3,6 na (i)); a P_f em cada duracao; tab/frag_f100.tex"""
    M = {"fragCemTem": "0"}
    rj = os.path.join(FRAG, "resumo.json")
    if not os.path.exists(rj):
        return M
    R = json.load(open(rj))
    L = []
    for var, cal, Ical in (("iii", "f100_iii_mc", 9.0), ("i", "f100_i_mc", 3.6)):
        r = R.get(f"f100_{var}")
        if not r:
            continue
        det = rf.le_det(os.path.join(RUNS, "sigmoide_f100_72h"), f"{var}_", 72.)
        b = {}
        for t in (72., 120.):
            ln_ = np.log([x for x in rf.le_mc(os.path.join(FRAG, cal), t).values() if x > 0])
            b[t] = (float(np.std(ln_, ddof=1)), float(np.exp(ln_.mean())), len(ln_))
        beta, med, ncal = b[72.]
        m = med / det[Ical]
        M.update({"fragCemN" + var: str(r["N"]), "fragCemBeta" + var: num(beta, 3), "fragCemM" + var: num(m, 3),
                  "fragCemBetaCxx" + var: num(b[120.][0], 3), "fragCemIcal" + var: num(Ical, 1).replace("{,}0", ""),
                  "fragCemNcal" + var: str(ncal)})
        for I, x in sorted(r["intensidades"].items(), key=lambda a: -float(a[0])):
            N, k_, p, F = x["N"], x["falhas"], x["Pf"], x["F_det"]
            q = float(norm.cdf(-math.log(m * F) / beta))
            z = (p - q) / math.sqrt(q * (1 - q) / N)
            lo, hi = x["IC"]
            t1 = np.array(x["t1_falhas"])
            tq = (f"{num(np.median(t1), 1)} [{num(np.quantile(t1, .1), 1)}, {num(np.quantile(t1, .9), 1)}]"
                  if len(t1) else "--")
            du = x["Pf_duracao"]
            L.append(f"{NOMEV[var]} & {num(float(I), 1).replace('{,}0', '')} & {num(F, 3)} & {k_}/{N} & {pcs(p)} & "
                     f"{pcs(lo)}--{pcs(hi)} & {pcs(q)} ({num(z, 1)}) & {pcs(du.get('48'))} & {pcs(du.get('60'))} & "
                     f"{pcs(du.get('96'))} & {tq}\\\\")
        L.append("\\midrule")
        M["fragCemTem"] = "1"
    L = L[:-1]
    with open(os.path.join(TABS, "frag_f100.tex"), "w") as fo:
        fo.write("\\begin{tabular}{llrrrlllll l}\n\\toprule\nvariante & $I$ & $F$ (72 h) & falhas/$N$ & $\\Pf$ (\\%) & "
                 "IC 95\\,\\% & prev. ($z$) & 48 h & 60 h & 96 h & $t_1$ (h)\\\\\n\\midrule\n" + "\n".join(L)
                 + "\n\\bottomrule\n\\end{tabular}\n")
    return M
