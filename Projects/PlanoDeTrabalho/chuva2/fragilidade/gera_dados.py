#!/usr/bin/env python3
"""gera_dados.py: figuras, tabelas e numeros do relatorio da proposta da curva de fragilidade P_f(I) da chuva
(relatorio_fragilidade.tex), so com os dados que ja existem (25/11/2026):
  - as sigmoides deterministicas (FS das medias no fim da chuva): ~/projects/runs/sigmoide_f50 (f50, 24 h, iii em I*/ e
    i em i_I*/) e ~/projects/runs/sigmoide_f100_72h (f100, 72 h, iii_I*/ e i_I*/, FS em 24, 48 e 72 h);
  - o MC direto da referencia (f50, iii, 31,05 mm/h, 5000 realizacoes): ~/projects/runs/fermi_chuva2/chuva2/ref;
  - a variante i em 16 realizacoes da referencia: ~/projects/runs/variante_i_campo;
  - os testes da primeira passagem (crpc=): ~/projects/runs/frag_teste;
  - a primeira passagem deterministica dos 3 solos: ~/projects/runs/det_pc_3solos.
Saidas: figs/*.pdf|png, tab/*.tex, numeros.tex. Uso: python3 gera_dados.py"""
import csv
import glob
import math
import os
import re

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from scipy.optimize import brentq, curve_fit  # noqa: E402
from scipy.stats import kurtosis, norm, skew  # noqa: E402

AQUI = os.path.dirname(os.path.abspath(__file__))
RUNS = os.path.expanduser("~/projects/runs")
FIGS, TABS = os.path.join(AQUI, "figs"), os.path.join(AQUI, "tab")
os.makedirs(FIGS, exist_ok=True)
os.makedirs(TABS, exist_ok=True)
TINTA, TINTA2, GRADE, SUPERF, FAIXA = "#0b0b0b", "#5c5c58", "#e4e3df", "#fcfcfb", "#f0efec"
COR = {"iii": "#2a78d6", "i": "#1baf7a"}
NOME = {"iii": "(iii) desacoplado, ar passivo", "i": "(i) acoplado, ar ativo"}
plt.rcParams.update({"font.size": 9, "savefig.facecolor": SUPERF, "figure.facecolor": SUPERF})
N = {}  # macros do numeros.tex


def v(x, nd=3):
    return f"{x:.{nd}f}".replace(".", "{,}")


def pct(x, nd=2):
    return f"{100 * x:.{nd}f}".replace(".", "{,}")


def estilo(ax, titulo="", xl="", yl=""):
    ax.set_facecolor(SUPERF)
    ax.set_title(titulo, color=TINTA, fontsize=9.5, loc="left")
    ax.set_xlabel(xl, color=TINTA)
    ax.set_ylabel(yl, color=TINTA)
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


def le_det(pasta, prefixo, t):
    d = {}
    for p in glob.glob(os.path.join(pasta, prefixo + "I*")):
        m = re.match(re.escape(prefixo) + r"I([\d.]+)$", os.path.basename(p))
        arq = os.path.join(p, "det_0_-1.csv")
        if not m or not os.path.exists(arq):
            continue
        for row in csv.DictReader(open(arq)):
            if abs(float(row["t_h"]) - t) < 1e-9 and row["situacao"] == "ok":
                d[float(m.group(1))] = float(row["FS_min"])
    return dict(sorted(d.items()))


def sig(x, A, B, C, D):
    """a Eq. 6 do Rahardjo et al. (2007): F = A + B / (1 + exp(-(I - C) / D))"""
    return A + B / (1 + np.exp(-(np.asarray(x, float) - C) / D))


def sig_log(x, A, B, C, D):
    """a mesma sigmoide em ln I (a escala da Fig. 6 do artigo)"""
    return A + B / (1 + np.exp(-(np.log(np.asarray(x, float)) - C) / D))


# FS inicial (sem chuva o FS nao muda): o limite de F(I) com I -> 0
FS0 = 2.0933
_a0 = os.path.join(RUNS, "det_pc_3solos", "f50_I36", "det_0_-1.csv")
if os.path.exists(_a0):
    FS0 = next(float(r["FS_min"]) for r in csv.DictReader(open(_a0)) if float(r["t_h"]) == 0.)
FORMAS = {  # nome -> (funcao, A livre?)
    "Eq. 6, 4 constantes": (sig, True),
    "Eq. 6, A = FS inicial": (sig, False),
    "em ln I, A = FS inicial": (sig_log, False),
}


def ajusta(det):
    """as tres formas: {forma: (parametros A, B, C, D, rms, maior residuo perto de FS = 1 (0,9 <= F <= 1,25))}"""
    I, F = np.array(list(det)), np.array(list(det.values()))
    perto = (F >= 0.9) & (F <= 1.25)
    out = {}
    for nome, (f, livre) in FORMAS.items():
        x = np.log(I) if f is sig_log else I
        try:
            if livre:
                p, _ = curve_fit(sig, I, F, p0=[F.max(), F.min() - F.max(), np.median(I), 8.], maxfev=100000)
            else:
                g = (lambda xx, B, C, D: FS0 + B / (1 + np.exp(-(xx - C) / D)))
                q, _ = curve_fit(g, x, F, p0=[F.min() - FS0, np.median(x), (x.max() - x.min()) / 6], maxfev=100000)
                p = np.array([FS0, *q])
            r = F - f(I, *p)
            out[nome] = (p, float(np.sqrt(np.mean(r ** 2))), float(np.abs(r[perto]).max()) if perto.any() else float("nan"))
        except RuntimeError:
            continue
    return out


def adotada(aj):
    """a forma com o menor residuo perto de FS = 1"""
    return min(aj, key=lambda k: (aj[k][2], aj[k][1]))


def F_sig(chave, x):
    nome = ADOT[chave]
    f = FORMAS[nome][0]
    return f(x, *AJ[chave][nome][0])


def Ic_curva(chave, d):
    xs = np.linspace(min(d), max(d), 4000)
    ys = F_sig(chave, xs) - 1
    for a, b, ya, yb in zip(xs, xs[1:], ys, ys[1:]):
        if ya >= 0 > yb:
            return brentq(lambda z: float(F_sig(chave, z)) - 1, a, b)
    return None


def interp_lnF(det, I):
    xs = list(det)
    I = min(max(I, xs[0]), xs[-1])
    for a, b in zip(xs, xs[1:]):
        if a <= I <= b:
            la, lb = math.log(det[a]), math.log(det[b])
            return la + (lb - la) * (I - a) / (b - a)
    return math.log(det[xs[-1]])


def Ic_pontos(det):
    xs = list(det)
    for a, b in zip(xs, xs[1:]):
        if det[a] >= 1 > det[b]:
            return a + (det[a] - 1) / (det[a] - det[b]) * (b - a)
    return None


# ---------------------------------------------------------------- 1. sigmoides deterministicas
S50 = os.path.join(RUNS, "sigmoide_f50")
S100 = os.path.join(RUNS, "sigmoide_f100_72h")
SIG = {("f50", "iii"): le_det(S50, "", 24.), ("f50", "i"): le_det(S50, "i_", 24.),
       ("f100", "iii"): le_det(S100, "iii_", 72.), ("f100", "i"): le_det(S100, "i_", 72.)}
F100 = {(var, t): le_det(S100, f"{var}_", t) for var in ("iii", "i") for t in (24., 48., 72.)}
AJ = {k: ajusta(d) for k, d in SIG.items()}
ADOT = {k: adotada(a) for k, a in AJ.items()}
N["FSzero"] = v(FS0, 3)
linhas = []
for (solo, var), d in SIG.items():
    D = 24 if solo == "f50" else 72
    icp = Ic_pontos(d)
    for nome, (p, rms, perto) in AJ[(solo, var)].items():
        A, B, C, Dd = p
        marca = "$\\star$" if nome == ADOT[(solo, var)] else ""
        cpar = f"{v(C, 2)}" if FORMAS[nome][0] is sig else f"{v(C, 3)} ($e^C$ = {v(math.exp(C), 1)})"
        linhas.append(f"{solo}, {D} h, ({var}) & {nome}{marca} & {v(A)} & {v(B)} & {cpar} & {v(Dd, 3)} & {v(A + B)} & "
                      f"{v(rms)} & {v(perto)}\\\\")
    ic = Ic_curva((solo, var), d)
    linhas.append(f"\\multicolumn{{9}}{{l}}{{\\quad $I_c$ (FS = 1): {'não cruza' if icp is None else v(icp, 1) + ' mm/h nos pontos'}"
                  f"{'' if ic is None else ', ' + v(ic, 1) + ' mm/h na forma adotada'}; patamar inferior nos pontos: "
                  f"{v(min(d.values()))}}}\\\\")
    linhas.append("\\midrule")
    k = f"{solo}{var}".replace("f50", "cinq").replace("f100", "cem")
    N[f"pat{k}"] = v(min(d.values()))
    N[f"ic{k}"] = "--" if icp is None else v(icp, 1)
open(os.path.join(TABS, "sigmoides.tex"), "w").write("\n".join(linhas[:-1]) + "\n")

fig, axs = plt.subplots(1, 2, figsize=(14, 4.6))
for ax, solo, D in ((axs[0], "f50", 24), (axs[1], "f100", 72)):
    estilo(ax, f"{solo}: FS das médias no fim de uma chuva de {D} h", "intensidade da chuva I (mm/h)", "FS (SRM)")
    ax.axhline(1, color=TINTA2, lw=1, ls=(0, (4, 3)))
    ks = 36 if solo == "f50" else 3.6
    ax.axvline(ks, color=TINTA2, lw=0.8, ls=(0, (1, 2)))
    ax.text(ks, 0.02, " I = k_s", transform=ax.get_xaxis_transform(), fontsize=8, color=TINTA2)
    for var in ("iii", "i"):
        d = SIG[(solo, var)]
        xs = np.geomspace(min(d), max(d), 400)
        ax.plot(xs, F_sig((solo, var), xs), color=COR[var], lw=2, label=f"{NOME[var]}: {ADOT[(solo, var)]}")
        ax.plot(list(d), list(d.values()), "o", ms=6, color=COR[var], mec=SUPERF, mew=1.2)
        icp = Ic_pontos(d)
        if icp:
            ax.plot(icp, 1, "v", ms=10, color=COR[var], mec=SUPERF, mew=1.5)
            ax.annotate(f"I_c = {icp:.1f}".replace(".", ","), (icp, 1), xytext=(4, 8) if var == "iii" else (-58, -14),
                        textcoords="offset points", fontsize=8, color=TINTA)
    ax.set_xscale("log")
    ax.legend(frameon=False, fontsize=7.5, labelcolor=TINTA)
fig.tight_layout()
salva(fig, "sigmoides")


def sigma_I(det, Ic, beta):
    """a largura da fragilidade em I_c: beta / |d ln F / dI| (diferenca nos dois pontos que cercam I_c)"""
    xs = list(det)
    for a, b in zip(xs, xs[1:]):
        if a <= Ic <= b:
            return beta / abs((math.log(det[b]) - math.log(det[a])) / (b - a))
    return None

# ---------------------------------------------------------------- 2. a referencia (MC de 5000)
REF = os.path.join(RUNS, "fermi_chuva2", "chuva2", "ref")
fs24 = {}
for arq in glob.glob(os.path.join(REF, "mc_*.csv")):
    if "_mapa_" in arq:
        continue
    for row in csv.DictReader(open(arq)):
        n = int(row["real"])
        if n >= 0 and row["t_h"] == "24" and row["situacao"] == "ok":
            fs24[n] = float(row["FS_min"])
ln = np.log(np.array([fs24[k] for k in sorted(fs24)]))
nref, kf = len(ln), int(np.sum(ln < 0))
pf = kf / nref
z = 1.96
c = (pf + z * z / (2 * nref)) / (1 + z * z / nref)
h = z * math.sqrt(pf * (1 - pf) / nref + z * z / (4 * nref * nref)) / (1 + z * z / nref)
beta, med, mu = float(np.std(ln, ddof=1)), float(np.exp(np.median(ln))), float(np.mean(ln))
Fdet = SIG[("f50", "iii")][31.05]
prev = {"lognormal ajustada (média e desvio de ln FS)": norm.cdf(-mu / beta),
        "mediana do MC e $\\beta$": norm.cdf(-math.log(med) / beta),
        "FS determinístico e $\\beta$": norm.cdf(-math.log(Fdet) / beta),
        "sigmoide adotada e $\\beta$": norm.cdf(-math.log(float(F_sig(("f50", "iii"), 31.05))) / beta)}
N.update(nRef=str(nref), falhasRef=str(kf), pfRef=pct(pf), icRefLo=pct(c - h), icRefHi=pct(c + h), betaRef=v(beta, 4),
         medRef=v(med, 4), fdetRef=v(Fdet, 4), assRef=v(float(skew(ln)), 2), curtRef=v(float(kurtosis(ln)), 2),
         prevMed=pct(prev["mediana do MC e $\\beta$"]), prevDet=pct(prev["FS determinístico e $\\beta$"]),
         prevEqSeis=pct(prev["sigmoide adotada e $\\beta$"]),
         prevLogn=pct(prev["lognormal ajustada (média e desvio de ln FS)"]))
L = [f"{k} & {pct(x)}\\,\\% & " + f"{(x - pf) / math.sqrt(pf * (1 - pf) / nref):+.2f}".replace(".", "{,}") + "\\\\"
     for k, x in prev.items()]
for var in ("iii", "i"):
    d = SIG[("f50", var)]
    ic = Ic_pontos(d)
    N[f"sigmaI{var}"] = v(sigma_I(d, ic, beta), 1)
    k1 = [x for x in d if x < ic][-1]
    N[f"declive{var}"] = v(abs(math.log(d[[x for x in d if x > ic][0]]) - math.log(d[k1])) /
                           ([x for x in d if x > ic][0] - k1), 4)
N["fatorMediana"] = v(med / Fdet, 4)
open(os.path.join(TABS, "ref.tex"), "w").write("\n".join(L) + "\n")

fig, axs = plt.subplots(1, 2, figsize=(13, 4.2))
estilo(axs[0], f"ln FS em 24 h, {nref} realizações (f50, 31,05 mm/h, iii)", "ln FS", "densidade")
axs[0].hist(ln, bins=60, density=True, color="#9ec5f4", edgecolor=SUPERF, lw=0.4, zorder=2)
xs = np.linspace(ln.min(), ln.max(), 300)
axs[0].plot(xs, norm.pdf(xs, mu, beta), color=COR["iii"], lw=2, zorder=3, label=f"normal (μ = {mu:.4f}, β = {beta:.4f})")
axs[0].axvline(0, color="#e34948", lw=1.5, zorder=3)
axs[0].text(0.002, 0.95, "FS = 1", transform=axs[0].get_xaxis_transform(), color=TINTA, fontsize=8)
axs[0].legend(frameon=False, fontsize=8, labelcolor=TINTA)
estilo(axs[1], "QQ: ln FS contra a normal", "quantil teórico (normal padrão)", "ln FS")
q = np.sort(ln)
tq = norm.ppf((np.arange(1, nref + 1) - 0.5) / nref)
axs[1].plot(tq, q, ".", ms=2.5, color=COR["iii"], zorder=3)
axs[1].plot(tq, mu + beta * tq, color=TINTA2, lw=1, zorder=2)
axs[1].axhline(0, color="#e34948", lw=1)
fig.tight_layout()
salva(fig, "referencia")

# ---------------------------------------------------------------- 3. a variante i em 16 realizacoes
VI = os.path.join(RUNS, "variante_i_campo")
difs = []
if os.path.isdir(VI):
    for d in sorted(glob.glob(os.path.join(VI, "i_*"))):
        m = re.search(r"i_(\d+)$", d)
        arq = glob.glob(os.path.join(d, "mc_*_*.csv"))
        arq = [a for a in arq if "_mapa_" not in a]
        if not m or not arq:
            continue
        for row in csv.DictReader(open(arq[0])):
            if row["t_h"] == "24" and int(row["real"]) >= 0:
                n = int(m.group(1))
                difs.append((n, fs24.get(n), float(row["FS_min"])))
if difs:
    dl = [math.log(b) - math.log(a) for _, a, b in difs if a]
    N.update(nVi=str(len(difs)), dlnMin=v(min(dl), 3), dlnMax=v(max(dl), 3), dlnMed=v(float(np.mean(dl)), 3),
             dlnSd=v(float(np.std(dl, ddof=1)), 3))

# ---------------------------------------------------------------- 4. as curvas de fragilidade previstas (f50) e o f100
fig, axs = plt.subplots(2, 1, figsize=(8.5, 7.6), sharex=True, gridspec_kw={"height_ratios": [1, 1.15]})
ax0, ax1 = axs
estilo(ax0, "o mecanismo: o FS das realizações espalha em torno de F(I) (faixa de 5 % a 95 %)", "", "FS no fim da chuva")
estilo(ax1, "a P_f prevista, P_f(I) = Φ(−ln F(I)/β)", "intensidade da chuva I (mm/h), f50, 24 h", "P_f (%)")
ax0.axhline(1, color=TINTA2, lw=1, ls=(0, (4, 3)))
for var in ("iii", "i"):
    d = SIG[("f50", var)]
    xs = np.linspace(10, 60, 400)
    Fp = np.array([math.exp(interp_lnF(d, x)) for x in xs])
    ax0.plot(xs, Fp, color=COR[var], lw=2, label=NOME[var])
    ax0.fill_between(xs, Fp * math.exp(-1.645 * beta), Fp * math.exp(1.645 * beta), color=COR[var], alpha=0.15, lw=0)
    ax1.plot(xs, 100 * norm.cdf(-np.log(Fp) / beta), color=COR[var], lw=2, label=f"{NOME[var]}: pontos determinísticos")
    ax1.plot(xs, 100 * norm.cdf(-np.log(F_sig(("f50", var), xs)) / beta), color=COR[var], lw=1.2, ls=(0, (4, 2)),
             label=f"({var}): sigmoide adotada")
    for b in (0.06, 0.10):
        ax1.plot(xs, 100 * norm.cdf(-np.log(Fp) / b), color=COR[var], lw=0.7, ls=(0, (1, 2)))
ax0.legend(frameon=False, fontsize=8, labelcolor=TINTA)
ax1.errorbar([31.05], [100 * pf], yerr=[[100 * (pf - (c - h))], [100 * ((c + h) - pf)]], fmt="s", ms=7,
             color=TINTA, capsize=3, zorder=5, label=f"MC direto (5000): {100 * pf:.2f} %".replace(".", ","))
ax1.text(58, 2, "pontilhado: β = 0,06 e 0,10", ha="right", fontsize=7.5, color=TINTA2)
ax1.legend(frameon=False, fontsize=7.5, labelcolor=TINTA, loc="upper left")
ax1.set_xlim(10, 60)
fig.tight_layout()
salva(fig, "fragilidade_prevista_f50")

# previsoes para a campanha (tabela)
L = []
for var, Is in (("iii", [27, 30, 31.05, 33, 36, 40, 45]), ("i", [18, 20, 22, 24, 36])):
    d = SIG[("f50", var)]
    for I in Is:
        F = math.exp(interp_lnF(d, I))
        F6 = float(F_sig(("f50", var), I))
        L.append(f"({var}) & {str(I).replace('.', '{,}')} & {v(F)} & {pct(norm.cdf(-math.log(F) / beta), 1)} & "
                 f"{pct(norm.cdf(-math.log(F6) / beta), 1)} & {pct(norm.cdf(-math.log(F) / 0.06), 1)} & "
                 f"{pct(norm.cdf(-math.log(F) / 0.10), 1)}\\\\")
    L.append("\\midrule")
open(os.path.join(TABS, "previsoes.tex"), "w").write("\n".join(L[:-1]) + "\n")

# f100: P_f prevista em 24, 48 e 72 h (beta da referencia)
fig, axs = plt.subplots(1, 2, figsize=(13, 4.2), sharey=True)
cores_t = {24.: "#86b6ef", 48.: "#2a78d6", 72.: "#0d366b"}
for ax, var in zip(axs, ("iii", "i")):
    estilo(ax, f"f100, {NOME[var]}: P_f prevista com β = {beta:.3f}".replace(".", ","), "intensidade da chuva I (mm/h)",
           "P_f (%)" if var == "iii" else "")
    for t in (24., 48., 72.):
        d = F100[(var, t)]
        if not d:
            continue
        xs = np.linspace(min(d), max(d), 300)
        ax.plot(xs, [100 * norm.cdf(-interp_lnF(d, x) / beta) for x in xs], color=cores_t[t], lw=2,
                label=f"chuva de {t:.0f} h")
    ax.set_xscale("log")
    ax.legend(frameon=False, fontsize=8, labelcolor=TINTA)
fig.tight_layout()
salva(fig, "fragilidade_prevista_f100")
for var in ("iii", "i"):
    d72 = F100[(var, 72.)]
    if d72:
        N[f"pfMaxCem{var}"] = pct(norm.cdf(-math.log(min(d72.values())) / beta), 1)

# ---------------------------------------------------------------- 5. limiares I-D probabilisticos (niveis de F)
for p_alvo, nome in ((0.01, "Um"), (0.05, "Cinco"), (0.5, "Cinquenta")):
    N[f"nivel{nome}"] = v(math.exp(-beta * norm.ppf(p_alvo)), 3)

# ---------------------------------------------------------------- 6. a primeira passagem: testes do crpc e o deterministico
FT = os.path.join(RUNS, "frag_teste")
L = []
for arq, var in ((os.path.join(FT, "iii_1860_crpc_1860_1860.csv"), "iii"), (os.path.join(FT, "i_1860_crpc_1860_1860.csv"), "i"),
                 (os.path.join(FT, "iii_2521_crpc_2521_2521.csv"), "iii")):
    if not os.path.exists(arq):
        continue
    rows = list(csv.DictReader(open(arq)))
    n = rows[0]["real"]
    cel = ", ".join(f"{float(r['I']):g}".replace(".", "{,}") + ": " + (v(float(r["DC_h"]), 1) if r["falha"] == "sim" else
                                               f"não rompe (FS 24 h = {v(float(r['FS_lo']), 3)})") for r in rows)
    L.append(f"{n} & ({var}) & {v(fs24.get(int(n), float('nan')), 3)} & {cel}\\\\")
open(os.path.join(TABS, "crpc_teste.tex"), "w").write("\n".join(L) + "\n")

D3 = os.path.join(RUNS, "det_pc_3solos")
L = []
for solo, Is in (("f10", ["36", "90", "180"]), ("f50", ["9", "36", "90"]), ("f100", ["3.6", "9", "36"])):
    for I in Is:
        cel = []
        for var in ("iii", "ii", "i"):
            arq = os.path.join(D3, f"{solo}_I{I}" + ("" if var == "iii" else f"_{var}"), "det_0_-1.csv")
            if not os.path.exists(arq):
                cel += ["--", "--"]
                continue
            r = sorted((float(x["t_h"]), float(x["FS_min"])) for x in csv.DictReader(open(arq)) if x["situacao"] == "ok")
            fmin = min(f for _, f in r)
            t1 = None
            for (ta, fa), (tb, fb) in zip(r, r[1:]):
                if fa >= 1 > fb:
                    t1 = ta + (fa - 1) / (fa - fb) * (tb - ta)
                    break
            cel += [v(fmin, 3), "--" if t1 is None else v(t1, 1)]
        L.append(f"{solo} & {I.replace('.', '{,}')} & " + " & ".join(cel) + "\\\\")
    L.append("\\midrule")
open(os.path.join(TABS, "t1_det.tex"), "w").write("\n".join(L[:-1]) + "\n")

# ---------------------------------------------------------------- 7. testes com dados que ja existem
import json  # noqa: E402

# 7a. caso A do Hu et al. (2024): 6 cenarios de chuva (I, D) no mesmo talude, MC e LS por cenario; k_s fixo
HU = os.path.join(RUNS, "fermi_hu2024_ks2e-5")
hu = {}
if os.path.exists(os.path.join(HU, "resumo.json")):
    for r in json.load(open(os.path.join(HU, "resumo.json"))):
        cid, D = r["id"], r["D"]
        fsD, fmin = [], {}
        for arq in glob.glob(os.path.join(HU, cid, "mc", "mc_*.csv")):
            if "_mapa_" in arq:
                continue
            for row in csv.DictReader(open(arq)):
                n = int(row["real"])
                if n < 0 or row["situacao"] != "ok":
                    continue
                fv = float(row["FS_min"])
                fmin[n] = min(fmin.get(n, 9.), fv)
                if abs(float(row["t_h"]) - D) < 1e-9:
                    fsD.append(fv)
        if not fsD:
            continue
        lnD = np.log(np.array(fsD))
        pls = float(np.mean(r["ls"]["p"]))
        hu[cid] = dict(I=r["I"], D=D, F=r["det"]["FS_min"], N=len(fmin), nf=sum(1 for x in fmin.values() if x < 1),
                       med=float(np.exp(np.median(lnD))), beta=float(np.std(lnD, ddof=1)), sk=float(skew(lnD)), pls=pls,
                       cvls=float(np.std(r["ls"]["p"], ddof=1) / math.sqrt(len(r["ls"]["p"])) / pls) if pls > 0 else float("nan"))
if "I30_D24" in hu:
    cal = hu["I30_D24"]
    bc, mc = cal["beta"], cal["med"] / cal["F"]
    N.update(huBeta=v(bc, 4), huM=v(mc, 4))
    L = []
    for cid in ("I30_D12", "I20_D24", "I40_D12", "I20_D36", "I30_D24", "I50_D24"):
        x = hu.get(cid)
        if not x:
            continue
        prev = norm.cdf(-math.log(mc * x["F"]) / bc)
        sem_m = norm.cdf(-math.log(x["F"]) / bc)
        pmc = x["nf"] / x["N"]

        def g(p):
            return (f"{100 * p:.3g}".replace(".", "{,}") + "\\,\\%") if p >= 1e-6 else f"${p * 100:.1e}".replace("e-0", "\\cdot10^{-").replace("e-", "\\cdot10^{-") + "}$\\,\\%"
        L.append(f"{x['I']:g} & {x['D']:g} & {v(x['F'])} & {v(x['med'] / x['F'], 4)} & {v(x['beta'], 4)} & {v(x['sk'], 2)} & "
                 f"{x['nf']}/{x['N']} & {g(x['pls'])} ({x['cvls'] * 100:.0f}\\,\\%) & {g(prev)} & {g(sem_m)}"
                 + (" & calibração" if cid == "I30_D24" else " & ") + "\\\\")
    open(os.path.join(TABS, "hu.tex"), "w").write("\n".join(L) + "\n")

# 7b. caso B do Liu et al. (2024): cauda (P_f ~ 1e-5), MC e LS
LIU = os.path.join(RUNS, "fermi_liu2024")
L = []
if os.path.exists(os.path.join(LIU, "resumo.json")):
    RL = json.load(open(os.path.join(LIU, "resumo.json")))
    RL = RL if isinstance(RL, list) else list(RL.values())
    for r in RL:
        cid = r.get("id") or r.get("cen")
        fsl = [float(row["FS_min"]) for arq in glob.glob(os.path.join(LIU, cid, "mc", "*.csv")) if "mapa" not in arq
               for row in csv.DictReader(open(arq)) if int(row["real"]) >= 0 and row["situacao"] == "ok"]
        if not fsl:
            continue
        lnl = np.log(fsl)
        med_l, b_l = float(np.exp(np.median(lnl))), float(np.std(lnl, ddof=1))
        Fd = r["det"]["FS_min"]
        pls = float(np.mean(r["ls"]["p"]))
        cvl = float(np.std(r["ls"]["p"], ddof=1) / math.sqrt(len(r["ls"]["p"])) / pls)
        e = lambda p: f"${p:.2e}$".replace("e-0", "\\cdot10^{-").replace("e-", "\\cdot10^{-").replace("$", "}$", 1)[::-1].replace("$}", "}$", 1)[::-1]
        L.append(f"{cid} & {len(fsl)} & {v(Fd)} & {v(med_l / Fd, 4)} & {v(b_l, 4)} & {v(float(skew(lnl)), 2)} & "
                 f"{pls:.2e} ({cvl * 100:.0f}\\,\\%) & {norm.cdf(-math.log(med_l) / b_l):.2e} & {norm.cdf(-math.log(Fd) / b_l):.2e}\\\\")
open(os.path.join(TABS, "liu.tex"), "w").write("\n".join(L).replace("e-0", "e$-$").replace("e-", "e$-$") + "\n")

# 7c. a coluna do Cho (2014): FS em 6, 12, 18, 24 h e t_1 por realizacao (o contraexemplo)
CHO = os.path.join(RUNS, "fermi_cho2014")
L = []
for caso in ("c05_l0.5", "c10_l0.5", "c05_l4", "c15_l4", "c05_l16"):
    arq_det = os.path.join(CHO, caso, "det", "det_0_-1.csv")
    if not os.path.exists(arq_det):
        continue
    detc = {float(r["t_h"]): float(r["FS_min"]) for r in csv.DictReader(open(arq_det))}
    fsc, t1c = {}, {}
    for arq in glob.glob(os.path.join(CHO, caso, "mc", "*.csv")):
        if "mapa" in arq:
            continue
        for r in csv.DictReader(open(arq)):
            if r["situacao"] != "ok":
                continue
            n = int(r["real"])
            fsc.setdefault(float(r["t_h"]), {})[n] = float(r["FS_min"])
            t1c[n] = float(r["t_passagem"])
    for t in (18., 24.):
        if t not in fsc or t not in detc:
            continue
        vv = np.array([x for x in fsc[t].values() if x > 0])
        lnc = np.log(vv)
        medc, bcc = float(np.exp(np.median(lnc))), float(np.std(lnc, ddof=1))
        p1 = float(np.mean([1. if 0 <= x <= t + 1e-9 else 0. for x in t1c.values()]))
        L.append(f"{caso.replace('_', '\\_')} & {t:g} & {v(detc[t])} & {v(medc)} & {v(medc / detc[t], 2)} & {v(bcc, 3)} & "
                 f"{v(float(skew(lnc)), 2)} & {pct(p1)} & {pct(norm.cdf(-math.log(medc) / bcc))} & "
                 f"{pct(norm.cdf(-math.log(detc[t]) / bcc))}\\\\")
    L.append("\\midrule")
open(os.path.join(TABS, "cho.tex"), "w").write("\n".join(L[:-1]) + "\n")

# ---------------------------------------------------------------- 8. casos calculados em comum com outros artigos
def fp(p, nd=2):
    """P_f em %, com notacao cientifica abaixo de 0,01 %"""
    if p is None or (isinstance(p, float) and math.isnan(p)):
        return "--"
    if p <= 0:
        return "0"
    if 100 * p < 0.01:
        e = int(math.floor(math.log10(100 * p)))
        return f"${100 * p / 10 ** e:.1f}\\cdot10^{{{e}}}$".replace(".", "{,}")
    if 100 * p < 1:
        return f"{100 * p:.2g}".replace(".", "{,}") if 100 * p >= 0.1 else f"{100 * p:.3f}".replace(".", "{,}")
    return f"{100 * p:.{nd}f}".replace(".", "{,}")


def cv_mc(k, n):
    return math.sqrt((1 - k / n) / k) if k > 0 else float("nan")


def ccv(x):
    return "--" if x is None or (isinstance(x, float) and math.isnan(x)) else f"{100 * x:.0f}"


# 8a. Hu et al. (2024), caso A
L = []
if hu:
    cal = hu["I30_D24"]
    bc, mcal = cal["beta"], cal["med"] / cal["F"]
    RH = {r["id"]: r for r in json.load(open(os.path.join(HU, "resumo.json")))}
    for cid in ("I30_D12", "I20_D24", "I40_D12", "I20_D36", "I30_D24", "I50_D24"):
        x = hu[cid]
        prev = norm.cdf(-math.log(mcal * x["F"]) / bc)
        z = (prev - x["pls"]) / (x["pls"] * x["cvls"]) if x["pls"] > 0 and x["cvls"] == x["cvls"] else float("nan")
        L.append(f"{x['I']:g} & {x['D']:g} & {v(x['F'])} & {fp(RH[cid]['pf_artigo'])} & {x['nf']}/{x['N']} & "
                 f"{fp(x['nf'] / x['N'])} & {ccv(cv_mc(x['nf'], x['N']))} & {fp(x['pls'])} & {ccv(x['cvls'])} & "
                 f"{fp(prev)} & " + ("cal." if cid == "I30_D24" else (f"{z:+.1f}".replace('.', '{,}') if z == z and abs(z) < 99 else "--"))
                 + "\\\\")
open(os.path.join(TABS, "comum_hu.tex"), "w").write("\n".join(L) + "\n")

# 8b. Liu et al. (2024), caso B: o artigo (Tabelas 2 e 4), Yang et al. (2026) e nos
L = []
if os.path.exists(os.path.join(LIU, "resumo.json")):
    for r in RL:
        cid = r.get("id") or r.get("cen")
        fsl = [float(row["FS_min"]) for arq in glob.glob(os.path.join(LIU, cid, "mc", "*.csv")) if "mapa" not in arq
               for row in csv.DictReader(open(arq)) if int(row["real"]) >= 0 and row["situacao"] == "ok"]
        if not fsl:
            continue
        lnl = np.log(fsl)
        med_l, b_l = float(np.exp(np.median(lnl))), float(np.std(lnl, ddof=1))
        Fd = r["det"]["FS_min"]
        pls = float(np.mean(r["ls"]["p"]))
        cvl = float(np.std(r["ls"]["p"], ddof=1) / math.sqrt(len(r["ls"]["p"])) / pls)
        n_mc = len(fsl)
        sup = 1 - 0.05 ** (1 / n_mc)  # limite superior de 95 % com 0 falhas
        tab4 = {"S2": 4.08e-5, "S3": 1.19e-4}.get(cid)
        tab2 = {"S3": 4.0e-3}.get(cid)
        yang = {"S3": 3.75e-3}.get(cid)
        L.append(f"{cid} & {v(Fd)} & {fp(tab4)} & {fp(tab2) if tab2 else '--'} & {fp(yang) if yang else '--'} & "
                 f"0/{n_mc} ($<${fp(sup)}) & {fp(pls)} & {ccv(cvl)} & {fp(norm.cdf(-math.log(med_l) / b_l))} & "
                 f"{fp(norm.cdf(-math.log(Fd) / b_l))}\\\\")
open(os.path.join(TABS, "comum_liu.tex"), "w").write("\n".join(L) + "\n")

# 8c. Cho (2014): os 15 casos (P_f em 18 e 24 h: o artigo e o nosso MC do NeoPZ), e a formula em 24 h
L = []
rt = os.path.join(CHO, "resultado.txt")
if os.path.exists(rt):
    for linha in open(rt):
        m_ = re.match(r"^(c\S+)\s*\|\s*([\d.]+)\s*\|\s*(\d+)\s+\d+\s+\d+\s*\|\s*([\d.]+)\s*\[[^\]]*\]\s*([\d.]+)\s*\|"
                      r"\s*([\d.]+)\s*\[[^\]]*\]\s*([\d.]+)", linha)
        if not m_:
            continue
        caso, n = m_.group(1), int(m_.group(3))
        p18, c18, p24, c24 = (float(m_.group(i)) for i in (4, 5, 6, 7))
        # a formula em 24 h com o FS das medias e o beta do proprio MC
        prev = "--"
        arq_det = os.path.join(CHO, caso, "det", "det_0_-1.csv")
        vals = [float(r["FS_min"]) for a_ in glob.glob(os.path.join(CHO, caso, "mc", "*.csv")) if "mapa" not in a_
                for r in csv.DictReader(open(a_)) if r["situacao"] == "ok" and float(r["t_h"]) == 24. and float(r["FS_min"]) > 0]
        if vals and os.path.exists(arq_det):
            Fd = next(float(r["FS_min"]) for r in csv.DictReader(open(arq_det)) if float(r["t_h"]) == 24.)
            lnc = np.log(vals)
            prev = fp(norm.cdf(-math.log(Fd) / float(np.std(lnc, ddof=1))), 1)
        L.append(f"{caso.replace('_', '\\_')} & {n} & {fp(c18, 1)} & {fp(p18, 2)} & {ccv(cv_mc(round(p18 * n), n))} & "
                 f"{fp(c24, 1)} & {fp(p24, 1)} & {ccv(cv_mc(round(p24 * n), n))} & {prev}\\\\")
open(os.path.join(TABS, "comum_cho.tex"), "w").write("\n".join(L) + "\n")

# 8d. X. Liu et al. (2025), caso 1: o artigo, a nossa reproducao (Python) e a formula ajustada de 0 a 18 h
XL_t = [0, 6, 12, 18, 24]
XL_F = [1.253, 1.135, 0.971, 0.889, 0.316]          # rotulos da Fig. 12 (lidos)
XL_P = [0.0724, 0.429, 0.815, 0.938, 0.975]          # o artigo (texto e tabela)
XL_nos = {0: "7,0", 12: "81,5", 24: "96,5"}           # nossa reproducao (0 h: 10^4; 12 e 24 h: MC de 300)
zz = np.array([norm.ppf(p) for p in XL_P[:4]])
xx = np.log(np.array(XL_F[:4]))
bb, aa = np.polyfit(xx, zz, 1)                         # z = aa + bb ln F; bb = -1/beta, aa = -ln m / beta
beta_xl, m_xl = -1 / bb, math.exp(aa / bb)
N.update(betaXL=v(beta_xl, 3), mXL=v(m_xl, 3))
L = []
for t, F, P in zip(XL_t, XL_F, XL_P):
    L.append(f"{t} & {v(F)} & {fp(P, 1)} & {XL_nos.get(t, '--')} & {fp(norm.cdf(-math.log(m_xl * F) / beta_xl), 1)}"
             + (" (ajuste)" if t <= 18 else " (previsão)") + "\\\\")
open(os.path.join(TABS, "comum_xliu.tex"), "w").write("\n".join(L) + "\n")

# 8e. Vargas Ceron et al. (2025): os 26 casos das Secoes 6.1 a 6.5 (rebaixamento, Gamma pelo GIM)
CER = os.path.join(RUNS, "fermi_ceron_ls_casos")
L = []
if os.path.exists(os.path.join(CER, "resultado.csv")):
    rows = list(csv.DictReader(open(os.path.join(CER, "resultado.csv"))))
    est = {}
    for r in rows:
        g = [float(x["FS_min"]) for a_ in glob.glob(os.path.join(CER, r["caso"], "mc", "*.csv")) if "mapa" not in a_
             for x in csv.DictReader(open(a_)) if int(x["real"]) >= 0 and x["situacao"] == "ok" and float(x["FS_min"]) > 0]
        if g:
            lg = np.log(g)
            est[r["caso"]] = (float(np.exp(np.median(lg))), float(np.std(lg, ddof=1)), len(g))
    if "ref" in est:
        med0, b0, _ = est["ref"]
        m0 = med0 / float(next(r for r in rows if r["caso"] == "ref")["Gamma_det"])
        N.update(ceronBeta=v(b0, 3), ceronM=v(m0, 3))
    for r in rows:
        c_ = r["caso"]
        G = float(r["Gamma_det"])
        pls, cvl = float(r["Pf_LS"]), float(r["CoV_LS"])
        pa = float(r["Pf_artigo"]) if r["Pf_artigo"] else float("nan")
        cva = float(r["CoVPf_artigo"]) if r["CoVPf_artigo"] else float("nan")
        pmc, nmc = float(r["Pf_MC"]), int(float(r["MC_n"]))
        if c_ in est:
            med, bta, _ = est[c_]
            pA = norm.cdf(-math.log(med) / bta)  # a lognormal do proprio caso (H1)
            zA = (pA - pls) / (pls * cvl)
        else:
            bta, pA, zA = float("nan"), float("nan"), float("nan")
        fam_carga = c_.startswith("alfa") or c_.startswith("hw") or c_ == "ref"
        pB = norm.cdf(-math.log(m0 * G) / b0) if fam_carga else None  # beta e m da ref (H2): as familias que mudam Gamma
        zB = (pB - pls) / (pls * cvl) if pB is not None else None
        L.append(f"{c_.replace('_', '\\_')} & {r['secao']} & {v(G, 3)} & {fp(pa)} & {ccv(cva)} & {fp(pls)} & {ccv(cvl)} & "
                 f"{fp(pmc, 1)} & {v(bta, 3)} & {fp(pA)} & {f'{zA:+.1f}'.replace('.', '{,}')} & "
                 + (f"{fp(pB)} & {f'{zB:+.1f}'.replace('.', '{,}')}" if pB is not None else "-- & --") + "\\\\")
open(os.path.join(TABS, "comum_ceron.tex"), "w").write("\n".join(L) + "\n")

with open(os.path.join(AQUI, "numeros.tex"), "w") as f:
    f.write("% gerado por gera_dados.py; nao edite\n")
    for k, x in N.items():
        f.write(f"\\newcommand{{\\{k}}}{{{x}}}\n")
print("pronto:", ", ".join(f"{k}={x}" for k, x in N.items()))
