#!/usr/bin/env python3
"""gera_dados.py: as tabelas, os numeros e as figuras do relatorio_hu2024.tex (o caso A, Hu et al. 2024, Comput.
Geotech. 173:106546), a partir de:

  - <campanha>/resumo.json (resumo_hu2024.py, refeito pelo gera.sh): o deterministico, o FORM modal, o line sampling,
    o Monte Carlo direto e a conferencia mc13 de cada cenario;
  - <campanha>/<caso>/ls/linhas.csv e <caso>/mc/mc_*.csv: as linhas do LS e as realizacoes do MC (convergencia,
    raizes e histogramas), lidas pelas funcoes do proprio resumo_hu2024.py;
  - ~/projects/runs/hu2024 (o PC): os deterministicos com k_s 2e-5 e 5,2e-6 nas malhas tri2, tri3 e tri4, a chuva
    vertical, o FS em 24 h por intensidade e o limite homogeneo (homogeneo/roda.sh);
  - as curvas de Hu lidas dos vetores do PDF (referencias_selecionadas_claude/reports/figs/dados): Fig. 9 (media e
    desvio do FS no tempo a 50 mm/h por 24 h), Fig. 11 (FS em 24 h por intensidade) e Figs. 16 e 19 (fragilidade).

Grava tab/*.tex, numeros.tex e figs/*.pdf. Roda em segundos.

    python3 -I gera_dados.py [pasta da campanha = ~/projects/runs/fermi_hu2024_ks2e-5]"""

import csv
import glob
import json
import math
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
from scipy import stats  # noqa: E402

AQUI = os.path.dirname(os.path.abspath(__file__))
RAIZ = os.path.abspath(os.path.join(AQUI, "..", "..", ".."))
sys.path.insert(0, os.path.join(RAIZ, "Projects", "SlopeRandomFieldRefinement", "cluster"))
import resumo_hu2024 as rh  # noqa: E402

CAMPANHA = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/projects/runs/fermi_hu2024_ks2e-5")
PC = os.path.expanduser("~/projects/runs/hu2024")
HU = os.path.expanduser("~/Dropbox/UnsaturatedSoil/referencias/referencias_selecionadas_claude/reports/figs/dados")
FIGS = os.path.join(AQUI, "figs")
TAB = os.path.join(AQUI, "tab")
ORDEM = ["I30_D12", "I20_D24", "I40_D12", "I20_D36", "I50_D24", "I30_D24"]

# o estilo das figuras dos outros relatorios do plano (Projects/SlopeRandomField/comparacao_ceron/figuras_sat.py)
PALETA = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
TINTA, TINTA2, GRADE = "#0b0b0b", "#52514e", "#dddcd8"
plt.rcParams.update({
    "font.size": 8.5, "axes.labelsize": 8.5, "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "legend.fontsize": 7,
    "axes.edgecolor": TINTA2, "axes.labelcolor": TINTA, "xtick.color": TINTA2, "ytick.color": TINTA2,
    "axes.grid": True, "grid.color": GRADE, "grid.linewidth": 0.5, "lines.linewidth": 1.4, "legend.frameon": False,
    "savefig.dpi": 200, "font.family": "DejaVu Sans", "mathtext.fontset": "dejavusans",
})
COR_D = {12: PALETA[0], 24: PALETA[1], 36: PALETA[2], 48: PALETA[6]}  # a cor segue a duracao da chuva
HU_COR = TINTA2


# ---------------------------------------------------------------------------------------------------------- formatos
def br(x, casas=2):
    """numero com virgula decimal para o LaTeX"""
    if x is None or (isinstance(x, float) and not math.isfinite(x)):
        return "--"
    t = f"{x:.{casas}f}".replace(".", "{,}")
    return "$-$" + t[1:] if t.startswith("-") else t


def pct(x, casas=2):
    return br(None if x is None else 100. * x, casas)


def cient(x, casas=1):
    """x em notacao cientifica para o LaTeX: 3{,}0 \\times 10^{-17}"""
    if x is None:
        return "--"
    if x == 0:
        return "0"
    e = int(math.floor(math.log10(abs(x))))
    m = x / 10 ** e
    if round(m, casas) >= 10:
        m, e = m / 10, e + 1
    return "$" + f"{m:.{casas}f}".replace(".", "{,}") + r"\times 10^{" + str(e) + "}$"


def pf_txt(x):
    """P_f em %: fixo de 0,01 % para cima, cientifico (em fracao) abaixo"""
    if x is None:
        return "--"
    if x == 0:
        return "0"
    if x < 1e-4:
        return cient(x)
    return pct(x, 3 if x < 0.01 else 2)


def painel(ax, t):
    ax.text(-0.02, 1.02, "(%s)" % t, transform=ax.transAxes, ha="right", va="bottom", fontsize=9, color=TINTA)


def salva(fig, nome):
    fig.savefig(os.path.join(FIGS, nome), metadata={"CreationDate": None, "ModDate": None})
    plt.close(fig)


def wilson(k, n, z=1.96):
    return rh.wilson(k, n, z)


def nome_caso(i):
    c = i.split("_")
    return "%s mm/h, %s h" % (c[0][1:], c[1][1:])


# ------------------------------------------------------------------------------------------------------------ dados
def carrega_campanha():
    casos = json.load(open(os.path.join(CAMPANHA, "resumo.json")))
    por = {c["id"]: c for c in casos}
    for i, c in por.items():
        d = os.path.join(CAMPANHA, i)
        # as linhas do LS na ordem do indice (o estimador corrente), com a raiz c e as analises de cada uma
        lin = []
        with open(os.path.join(d, "ls", "linhas.csv"), newline="") as f:
            for r in csv.DictReader(f):
                lin.append((int(r["linha"]), float(r["c"]), int(float(r["analises"])), float(r["p"])))
        c["linhas_csv"] = sorted(lin)
        c["ls_json"] = json.load(open(os.path.join(d, "ls", "ls.json")))
        # as realizacoes do MC: o FS no fim da chuva e o menor dos instantes, na ordem da realizacao
        arqs = [a for a in glob.glob(os.path.join(d, "mc", "mc_*_*.csv")) if not os.path.basename(a).startswith("mc_mapa")]
        reais, _ = rh.le_linhas(arqs)
        nts = max(len(v) for v in reais.values())
        reais = {n: v for n, v in reais.items() if len(v) == nts}
        ns = sorted(reais)
        c["mc_fsD"] = np.array([reais[n][c["D"]] for n in ns])
        c["mc_min"] = np.array([min(reais[n].values()) for n in ns])
    return casos, por


def le_det(arq):
    """{t: FS} do csv do deterministico (realizacao -1)"""
    reais, _ = rh.le_linhas([arq])
    return dict(sorted(reais[-1].items()))


def carrega_pc():
    det = os.path.join(PC, "det")
    pc = {"malha": {}, "fs24": {}}
    for ks in ("2e-5", "5.2e-6"):
        for m in ("tri2", "tri3", "tri4"):
            a = os.path.join(det, f"detB_ks{ks}_{m}_0_-1.csv")
            if os.path.exists(a):
                pc["malha"][(ks, m)] = le_det(a)
    pc["I30D12"] = {m: le_det(os.path.join(det, f"detI30D12_{m}_0_-1.csv"))[12.0] for m in ("tri3", "tri4")}
    pc["vertical"] = le_det(os.path.join(det, "vertical_I50_0_-1.csv"))[24.0]
    for I in (10, 30, 50, 70, 90):
        pc["fs24"][("2e-5", "tri2", I)] = le_det(os.path.join(det, "fs24", f"fs24_ks2e-5_I{I}_0_-1.csv"))[24.0]
        pc["fs24"][("5.2e-6", "tri2", I)] = le_det(os.path.join(det, "fs24", f"fs24_ks5.2e-6_I{I}_0_-1.csv"))[24.0]
        pc["fs24"][("2e-5", "tri3", I)] = le_det(os.path.join(det, "fs24_tri3", f"fs24_I{I}_0_-1.csv"))[24.0]
    # o limite homogeneo: FS(c', phi' = 24) em 0 e 24 h (50 mm/h por 24 h)
    hom = []
    for a in glob.glob(os.path.join(PC, "homogeneo", "c*_0_-1.csv")):
        if "mapa" in a:
            continue
        c = float(os.path.basename(a)[1:].split("_")[0])
        fs = le_det(a)
        hom.append((c, fs[0.0], fs[24.0]))
    pc["homogeneo"] = sorted(hom)
    return pc


def carrega_hu():
    def le(nome):
        with open(os.path.join(HU, nome), newline="") as f:
            return list(csv.DictReader(f))
    f9 = le("hu_fig9_media_dp.csv")
    hu = {"fig9": np.array([[float(r["t_h"]), float(r["media_FS"]), float(r["dp_FS"])] for r in f9]),
          "fig11": {int(float(r["I_mm_h"])): {k: float(v) for k, v in r.items()} for r in le("hu_fig11_fs24.csv")}}
    frag = {}
    with open(os.path.join(HU, "hu_fragilidade.csv"), newline="") as f:
        rd = csv.reader(f)
        cab = next(rd)
        xs = np.array([float(x) for x in cab[1:]])
        for r in rd:
            frag[r[0]] = np.array([float(x) for x in r[1:]])
    hu["frag_I"], hu["frag"] = xs, frag
    return hu


# ------------------------------------------------------------------------------------------------- limite homogeneo
def limite_homogeneo(pc, n=400000, semente=2025):
    """c' e phi' variaveis aleatorias (solo homogeneo em cada realizacao). O FS do SRM e homogeneo de grau 1 em
    (c', tan phi'): FS = tan phi' g(c'/tan phi'), com g lida dos deterministicos com phi' = 24 graus (interpolacao
    linear na razao; extrapolacao linear fora da grade)"""
    h = np.array(pc["homogeneo"])
    t24 = math.tan(math.radians(24.))
    r_k = h[:, 0] / t24
    rng = np.random.default_rng(semente)
    z1 = rng.standard_normal(n)
    z2 = -0.5 * z1 + math.sqrt(1 - 0.25) * rng.standard_normal(n)

    def lognormal(z, media, cov):
        s = math.sqrt(math.log(1 + cov * cov))
        return np.exp(math.log(media) - s * s / 2 + s * z)
    c, phi = lognormal(z1, 8., 0.15), lognormal(z2, 24., 0.10)
    tp = np.tan(np.radians(phi))
    r = c / tp
    out = {"rmin": r_k.min(), "rmax": r_k.max(), "fora": float(np.mean((r < r_k.min()) | (r > r_k.max())))}
    for j, t in ((1, 0), (2, 24)):
        g = h[:, j] / t24
        gi = np.interp(r, r_k, g)
        lo, hi = r < r_k[0], r > r_k[-1]
        gi[lo] = g[0] + (r[lo] - r_k[0]) * (g[1] - g[0]) / (r_k[1] - r_k[0])
        gi[hi] = g[-1] + (r[hi] - r_k[-1]) * (g[-1] - g[-2]) / (r_k[-1] - r_k[-2])
        fs = tp * gi
        out[t] = {"media": float(fs.mean()), "dp": float(fs.std(ddof=1)), "pf": float(np.mean(fs < 1.)),
                  "amostra": fs[:20000]}
        out[t]["cov"] = out[t]["dp"] / out[t]["media"]
    # o peso de c' no FS das medias (Euler: w_c + w_tan = 1), pela derivada central em c = 8
    for j, t in ((1, 0), (2, 24)):
        cs, fs = h[:, 0], h[:, j]
        i8 = int(np.argmin(abs(cs - 8.)))
        d = (fs[i8 + 1] - fs[i8 - 1]) / (cs[i8 + 1] - cs[i8 - 1])
        out[t]["w_c"] = 8. * d / fs[i8]
    return out


# ----------------------------------------------------------------------------------------------------------- tabelas
def tabela_pf(casos, por):
    lin = [r"\begin{table}[htbp]\centering\footnotesize\setlength{\tabcolsep}{3pt}",
           r"\caption{$P(FS < 1)$ por cenário (primeira passagem até 72~h, em \%; abaixo de $10^{-4}$, em fração). "
           r"FORM: $\beta$ e $\Phi(-\beta)$ (* = parou nas 10 iterações sem convergir). LS: line sampling na direção do "
           r"FORM, com o $CoV$ e o IC de 95\,\% normal. MC: Monte Carlo direto, IC de Wilson. $z$: LS contra MC. "
           r"Hu: Figs.~16 e 19 (0,7781 impresso).}\label{tab:pf}",
           r"\begin{tabular}{lrrrrrrrrrr}\toprule",
           r" & \multicolumn{2}{c}{FORM} & \multicolumn{3}{c}{LS} & \multicolumn{3}{c}{MC} & & Hu \\"
           r"\cmidrule(lr){2-3}\cmidrule(lr){4-6}\cmidrule(lr){7-9}",
           r"cenário & $\beta$ & $P_f$ & $P_f$ & $CoV$ & IC 95\,\% & $n_f/N$ & $P_f$ & IC 95\,\% & $z$ & $P_f$ \\\midrule"]
    for i in ORDEM:
        c = por[i]
        f, l, m = c["form"], c["ls"], c["mc"]
        b = br(f["beta"], 2) + ("*" if not f["convergiu"] else "")
        if l["P_f"] < 1e-4:
            ic = "--"
        else:
            ic = "[%s; %s]" % (pf_txt(max(l["IC95"][0], 0.)), pf_txt(min(l["IC95"][1], 1.)))
        icm = "[%s; %s]" % (pf_txt(m["IC95"][0]), pf_txt(m["IC95"][1]))
        lin.append("%s & %s & %s & %s & %s & %s & %d/%d & %s & %s & %s & %s \\\\" % (
            nome_caso(i), b, pf_txt(f["P_f"]), pf_txt(l["P_f"]), pct(l["CoV"], 1 if l["CoV"] >= 0.001 else 2), ic, m["n_f"], m["S"],
            pf_txt(m["P_f"]), icm, ("$%+.1f$" % c["z"]).replace(".", "{,}") if c["z"] is not None else "--",
            pct(c["pf_artigo"], 2 if c["pf_artigo"] > 0.1 else 1)))
    lin += [r"\bottomrule\end{tabular}\end{table}"]
    return "\n".join(lin) + "\n"


def tabela_custo(por):
    lin = [r"\begin{table}[htbp]\centering\footnotesize",
           r"\caption{O custo por cenário, em análises do MEF (uma análise: o transiente lido e o SRM no fim da chuva e "
           r"6~h depois). FORM: iterações de uma análise cada (gradiente modal). LS: linhas e análises (a parada no "
           r"primeiro $N \ge 30$ com $CoV \le 5\,\%$, ou o teto). MC: realizações; mc13: as 20 primeiras com o FS de "
           r"6 em 6~h.}\label{tab:custo}",
           r"\begin{tabular}{lrrrrrl}\toprule",
           r"cenário & FORM & linhas & LS & MC & mc13 & parada do LS \\\midrule"]
    for i in ORDEM:
        c = por[i]
        mot = c["ls"]["motivo"]
        mot = ("teto de %d linhas" % c["ls"]["linhas"]) if mot.startswith("n maximo") else \
            (r"$CoV \le 5\,\%%$ com %s" % mot.split(" com ")[1].split(" (")[0])
        lin.append("%s & %d & %d & %d & %d & %d & %s \\\\" % (nome_caso(i), c["form"]["analises"], c["ls"]["linhas"],
                                                             c["ls"]["analises"], c["mc"]["S"], c["mc13"]["S"], mot))
    tot = [sum(por[i]["form"]["analises"] for i in ORDEM), sum(por[i]["ls"]["linhas"] for i in ORDEM),
           sum(por[i]["ls"]["analises"] for i in ORDEM), sum(por[i]["mc"]["S"] for i in ORDEM),
           sum(por[i]["mc13"]["S"] for i in ORDEM)]
    lin.append(r"\midrule total & %d & %d & %d & %d & %d & \\" % tuple(tot))
    lin += [r"\bottomrule\end{tabular}\end{table}"]
    return "\n".join(lin) + "\n"


def tabela_det(pc, por, hu):
    m = pc["malha"]
    f9 = hu["fig9"]

    def hu_t(t):
        return float(np.interp(t, f9[:, 0], f9[:, 1]))
    lin = [r"\begin{table}[htbp]\centering\footnotesize",
           r"\caption{FS determinístico (as médias) a 50~mm/h por 24~h: as malhas do transiente (tri2, tri3, tri4: "
           r"\texttt{TriGMesh} de níveis 2 a 4) e os dois $k_s$, contra a média de Hu (Fig.~9). O fator é sempre "
           r"na \texttt{TriGMesh(1)} com 4 níveis adaptativos.}\label{tab:det}",
           r"\begin{tabular}{llrrrrr}\toprule",
           r"$k_s$ (m/s) & transiente & 0 h & 12 h & 24 h & 48 h & 72 h \\\midrule"]
    for ks, rot in (("2e-5", r"$2 \times 10^{-5}$"), ("5.2e-6", r"$5{,}2 \times 10^{-6}$")):
        for mm in ("tri2", "tri3", "tri4"):
            if (ks, mm) in m:
                v = m[(ks, mm)]
                lin.append("%s & %s & %s \\\\" % (rot, mm, " & ".join(br(v[t], 3) for t in (0., 12., 24., 48., 72.))))
    lin.append(r"\midrule Hu (média das 180) & -- & %s \\" % " & ".join(br(hu_t(t), 3) for t in (0, 12, 24, 48, 72)))
    lin += [r"\bottomrule\end{tabular}\end{table}"]
    return "\n".join(lin) + "\n"


def tabela_normal(por, hu):
    """a aproximacao normal P_f = Phi((1 - mu)/sigma) no fim da chuva: com a media e o desvio de Hu (Fig. 11, 24 h de
    chuva) contra a P_f dele (Fig. 16); com os nossos (MC) contra a nossa (LS e MC)"""
    lin = [r"\begin{table}[htbp]\centering\footnotesize",
           r"\caption{A aproximação normal $P_f \approx \Phi\bigl((1 - \mu)/\sigma\bigr)$ com a média e o desvio do FS "
           r"no fim da chuva. Em cima, Hu: $\mu$ e $\sigma$ da Fig.~11 (24~h de chuva) contra a $P_f$ das "
           r"Figs.~16 e 19. Embaixo, este trabalho: $\mu$ e $\sigma$ do MC direto contra o LS e o MC. A última coluna "
           r"põe o desvio de Hu na nossa média.}\label{tab:normal}",
           r"\begin{tabular}{lrrrrrr}\toprule",
           r"cenário & $\mu$ & $\sigma$ & $CoV$ (\%) & $\Phi\bigl((1 - \mu)/\sigma\bigr)$ & $P_f$ & nossa $\mu$, "
           r"$\sigma$ de Hu \\\midrule",
           r"\multicolumn{7}{l}{\textit{Hu et al.\ (2024), 24~h de chuva}} \\"]
    f11, frag = hu["fig11"], hu["frag"]["D24h"]
    for I in (10, 30, 50, 70, 90):
        mu, sd = f11[I]["media_t24"], f11[I]["dp_t24"]
        pf = 0.7781 if I == 50 else float(np.interp(I, hu["frag_I"], frag))
        lin.append("%d mm/h, 24 h & %s & %s & %s & %s & %s & -- \\\\" % (
            I, br(mu, 3), br(sd, 3), br(100 * sd / mu, 1), pct(stats.norm.cdf((1 - mu) / sd), 1), pct(pf, 1)))
    lin.append(r"\midrule\multicolumn{7}{l}{\textit{este trabalho (MC; $P_f$ do LS e, entre parênteses, do MC)}} \\")
    for i in ORDEM:
        c = por[i]
        mu, sd = c["mc"]["media_D"], c["mc"]["dp_D"]
        mist = "--"
        if c["D"] == 24.:
            sdh = float(np.interp(c["I"], sorted(f11), [f11[k]["dp_t24"] for k in sorted(f11)]))
            mist = pct(stats.norm.cdf((1 - mu) / sdh), 1)
        lin.append("%s & %s & %s & %s & %s & %s (%s) & %s \\\\" % (
            nome_caso(i), br(mu, 3), br(sd, 3), br(100 * sd / mu, 1), pf_txt(stats.norm.cdf((1 - mu) / sd)),
            pf_txt(c["ls"]["P_f"]), pf_txt(c["mc"]["P_f"]), mist))
    lin += [r"\bottomrule\end{tabular}\end{table}"]
    return "\n".join(lin) + "\n"


def tabela_desvio(por, hu, hom):
    """o desvio do FS: o nosso (campos), o de Hu e o limite homogeneo, em 0 e 24 h a 50 mm/h por 24 h"""
    f9 = hu["fig9"]
    c = por["I50_D24"]
    sd0 = c["mc13"]["por_instante"]["0.0"]
    n13 = c["mc13"]["S"]
    q = stats.chi2.ppf([0.975, 0.025], n13 - 1)
    ic0 = [sd0[1] * math.sqrt((n13 - 1) / x) for x in q]
    lin = [r"\begin{table}[htbp]\centering\footnotesize",
           r"\caption{Média, desvio e $CoV$ do FS a 50~mm/h por 24~h, em 0~h e no fim da chuva: os nossos campos "
           r"(0~h: as %d realizações do mc13, com o IC de 95\,\%% do desvio; 24~h: as %d do MC), o limite homogêneo "
           r"($c'$ e $\varphi'$ como variáveis aleatórias, com os mesmos parâmetros) e Hu (Fig.~9).}\label{tab:desvio}" %
           (n13, c["mc"]["S"]),
           r"\begin{tabular}{lrrrrrr}\toprule",
           r" & \multicolumn{3}{c}{0 h} & \multicolumn{3}{c}{24 h} \\\cmidrule(lr){2-4}\cmidrule(lr){5-7}",
           r" & $\mu$ & $\sigma$ & $CoV$ (\%) & $\mu$ & $\sigma$ & $CoV$ (\%) \\\midrule"]
    lin.append("campos de Nyström ($M = 150$) & %s & %s [%s; %s] & %s & %s & %s & %s \\\\" % (
        br(sd0[0], 3), br(sd0[1], 3), br(ic0[0], 3), br(ic0[1], 3), br(100 * sd0[1] / sd0[0], 1),
        br(c["mc"]["media_D"], 3), br(c["mc"]["dp_D"], 3), br(100 * c["mc"]["dp_D"] / c["mc"]["media_D"], 1)))
    lin.append("limite homogêneo & %s & %s & %s & %s & %s & %s \\\\" % (
        br(hom[0]["media"], 3), br(hom[0]["dp"], 3), br(100 * hom[0]["cov"], 1), br(hom[24]["media"], 3),
        br(hom[24]["dp"], 3), br(100 * hom[24]["cov"], 1)))
    h0, h24 = f9[0], f9[np.argmin(abs(f9[:, 0] - 24))]
    lin.append("Hu (Fig.~9) & %s & %s & %s & %s & %s & %s \\\\" % (
        br(h0[1], 3), br(h0[2], 3), br(100 * h0[2] / h0[1], 1), br(h24[1], 3), br(h24[2], 3),
        br(100 * h24[2] / h24[1], 1)))
    lin += [r"\bottomrule\end{tabular}\end{table}"]
    return "\n".join(lin) + "\n"


# ------------------------------------------------------------------------------------------------------------ figuras
def figura_fs_t(por):
    fig, axs = plt.subplots(2, 3, figsize=(7.0, 4.4), sharex=True, sharey=True, constrained_layout=True)
    for k, (ax, i) in enumerate(zip(axs.flat, ORDEM)):
        c = por[i]
        cor = COR_D[int(c["D"])]
        pi = c["mc13"]["por_instante"]
        ts = sorted(float(t) for t in pi)
        mu = np.array([pi[str(t)][0] for t in ts])
        sd = np.array([pi[str(t)][1] for t in ts])
        ax.axvspan(0, c["D"], color=GRADE, alpha=0.45, lw=0)
        ax.fill_between(ts, mu - sd, mu + sd, color=cor, alpha=0.22, lw=0)
        ax.plot(ts, mu, color=cor, lw=1.4, ls="--")
        td = sorted(float(t) for t in c["det"]["FS"])
        ax.plot(td, [c["det"]["FS"][str(t)] for t in td], color=TINTA,
                lw=1.2, marker="o", ms=2.5)
        ax.axhline(1, color=TINTA2, lw=0.8)
        ax.set_title(nome_caso(i), fontsize=8.5)
        ax.set_xticks(range(0, 73, 12))
        if k >= 3:
            ax.set_xlabel("t (h)")
        if k % 3 == 0:
            ax.set_ylabel("FS")
    axs.flat[0].legend([Line2D([], [], color=TINTA, marker="o", ms=2.5), Line2D([], [], color=TINTA2, ls="--"),
                        Patch(color=TINTA2, alpha=0.3), Patch(color=GRADE, alpha=0.6)],
                       ["determinístico (médias)", "média do mc13", r"média $\pm$ desvio (mc13)", "chuva"],
                       loc="lower left", fontsize=6.5)
    salva(fig, "fig_fs_t.pdf")


def figura_fig9(por, pc, hu):
    f9 = hu["fig9"]
    c = por["I50_D24"]
    pi = c["mc13"]["por_instante"]
    ts = sorted(float(t) for t in pi)
    mu = np.array([pi[str(t)][0] for t in ts])
    sd = np.array([pi[str(t)][1] for t in ts])
    fig, ax = plt.subplots(figsize=(6.4, 4.0), constrained_layout=True)
    ax.fill_between(f9[:, 0], f9[:, 1] - f9[:, 2], f9[:, 1] + f9[:, 2], color=HU_COR, alpha=0.18, lw=0)
    ax.plot(f9[:, 0], f9[:, 1], color=HU_COR, lw=1.4, ls=(0, (5, 2)))
    ax.fill_between(ts, mu - sd, mu + sd, color=PALETA[1], alpha=0.3, lw=0)
    ax.plot(ts, mu, color=PALETA[1], lw=1.2, ls="--")
    m4 = pc["malha"][("2e-5", "tri4")]
    ax.plot(list(m4), list(m4.values()), color=TINTA, lw=1.3, marker="o", ms=3)
    m5 = pc["malha"][("5.2e-6", "tri3")]
    ax.plot(list(m5), list(m5.values()), color=PALETA[0], lw=1.2, marker="s", ms=3, ls=":")
    ax.errorbar([24.6], [c["mc"]["media_D"]], yerr=[c["mc"]["dp_D"]], color=PALETA[1], marker="D", ms=4, capsize=2,
                lw=1.2)
    ax.plot([24.6], [pc["vertical"]], marker="^", color=PALETA[2], ms=6, ls="none")
    ax.axvline(24, color=TINTA2, lw=0.8, ls=":")
    ax.text(23.5, 1.88, "fim da chuva", ha="right", va="top", fontsize=7, color=TINTA2)
    ax.axhline(1, color=TINTA2, lw=0.8)
    ax.set_xlabel("t (h)")
    ax.set_ylabel("FS")
    ax.set_xticks(range(0, 73, 6))
    ax.set_xlim(0, 72)
    ax.legend([Line2D([], [], color=HU_COR, ls=(0, (5, 2))), Patch(color=HU_COR, alpha=0.25),
               Line2D([], [], color=TINTA, marker="o", ms=3),
               Line2D([], [], color=PALETA[1], ls="--"), Patch(color=PALETA[1], alpha=0.3),
               Line2D([], [], color=PALETA[1], marker="D", ms=4),
               Line2D([], [], color=PALETA[0], marker="s", ms=3, ls=":"),
               Line2D([], [], color=PALETA[2], marker="^", ls="none", ms=6)],
              ["Hu, média (Fig. 9)", r"Hu, média $\pm$ desvio", r"nosso determinístico, $k_s = 2\times10^{-5}$",
               "nossa média (mc13, 20 realizações)", r"nossa média $\pm$ desvio (mc13)",
               r"nosso MC em 24 h (%d realizações)" % c["mc"]["S"], r"determinístico, $k_s = 5{,}2\times10^{-6}$ (tri3)",
               "determinístico, chuva vertical (tri3)"],
              loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=2, fontsize=6.5)
    ax.set_ylim(0.55, 1.95)
    salva(fig, "fig_fig9.pdf")


def figura_malha(pc):
    m = pc["malha"]
    fig, ax = plt.subplots(figsize=(6.0, 3.0), constrained_layout=True)
    est = {"tri2": (":", "s"), "tri3": ("--", "^"), "tri4": ("-", "o")}
    hs = []
    for (ks, mm), v in sorted(m.items()):
        cor = TINTA if ks == "2e-5" else PALETA[0]
        ls, mk = est[mm]
        ax.plot(list(v), list(v.values()), color=cor, ls=ls, marker=mk, ms=3, lw=1.1)
        hs.append((Line2D([], [], color=cor, ls=ls, marker=mk, ms=3),
                   (r"$k_s = 2\times10^{-5}$" if ks == "2e-5" else r"$k_s = 5{,}2\times10^{-6}$") + ", " + mm))
    ax.axvspan(0, 24, color=GRADE, alpha=0.35, lw=0)
    ax.axhline(1, color=TINTA2, lw=0.8)
    ax.set_xlabel("t (h)")
    ax.set_ylabel("FS (médias)")
    ax.set_xticks(range(0, 73, 6))
    ax.set_xlim(0, 72)
    ax.legend(*zip(*hs), loc="lower right", ncol=2, fontsize=6.5)
    salva(fig, "fig_malha.pdf")


def figura_fig11(por, pc, hu):
    f11 = hu["fig11"]
    Is = sorted(f11)
    fig, ax = plt.subplots(figsize=(5.6, 3.3), constrained_layout=True)
    mu = np.array([f11[I]["media_t24"] for I in Is])
    sd = np.array([f11[I]["dp_t24"] for I in Is])
    ax.fill_between(Is, mu - sd, mu + sd, color=HU_COR, alpha=0.18, lw=0)
    ax.plot(Is, mu, color=HU_COR, ls=(0, (5, 2)), marker="s", ms=4, mfc="white")
    for (ks, mm), cor, ls, mk in ((("2e-5", "tri3"), TINTA, "-", "o"), (("2e-5", "tri2"), TINTA, ":", "o"),
                                  (("5.2e-6", "tri2"), PALETA[0], ":", "s")):
        ax.plot(Is, [pc["fs24"][(ks, mm, I)] for I in Is], color=cor, ls=ls, marker=mk, ms=3.5, lw=1.1,
                mfc="white" if mm == "tri2" else cor)
    ids = ["I30_D24", "I50_D24"]
    ax.plot([por[i]["I"] for i in ids], [por[i]["det"]["FS"]["24.0"] for i in ids], ls="none", marker="*", ms=9,
            color=PALETA[3], mec=TINTA2, mew=0.5)
    ax.errorbar([por[i]["I"] + 1.2 for i in ids], [por[i]["mc"]["media_D"] for i in ids],
                yerr=[por[i]["mc"]["dp_D"] for i in ids], ls="none", marker="D", ms=4, color=PALETA[1], capsize=2)
    ax.axhline(1, color=TINTA2, lw=0.8)
    ax.set_xlabel("intensidade da chuva (mm/h), 24 h de chuva")
    ax.set_ylabel("FS em t = 24 h")
    ax.set_xticks(Is)
    ax.legend([Line2D([], [], color=HU_COR, ls=(0, (5, 2)), marker="s", mfc="white"), Patch(color=HU_COR, alpha=0.25),
               Line2D([], [], color=TINTA, marker="o", ms=3.5), Line2D([], [], color=TINTA, ls=":", marker="o",
                                                                          mfc="white", ms=3.5),
               Line2D([], [], color=PALETA[3], marker="*", ls="none", ms=9, mec=TINTA2, mew=0.5),
               Line2D([], [], color=PALETA[1], marker="D", ms=4),
               Line2D([], [], color=PALETA[0], ls=":", marker="s", mfc="white", ms=3.5)],
              ["Hu, média (Fig. 11)", r"Hu, média $\pm$ desvio", "nosso determinístico, tri3",
               "nosso determinístico, tri2", "nosso determinístico, tri4 (campanha)",
               r"nosso MC, média $\pm$ desvio", r"determinístico, $k_s = 5{,}2\times10^{-6}$, tri2"],
              loc="upper right", fontsize=6.3)
    salva(fig, "fig_fig11.pdf")


def figura_fragilidade(por, hu):
    xs, fr = hu["frag_I"], hu["frag"]
    fig, axs = plt.subplots(1, 2, figsize=(7.0, 3.4), constrained_layout=True)
    for k, ax in enumerate(axs):
        for D in (12, 24, 36, 48):
            y = fr["D%dh" % D]
            if k == 1:
                y = np.where(y > 0, y, np.nan)
            ax.plot(xs, y, color=COR_D[D], ls=(0, (5, 2)), lw=1.1, alpha=0.9 if D < 48 else 0.5)
        for i in ORDEM:
            c = por[i]
            cor = COR_D[int(c["D"])]
            I = c["I"]
            l, m, f = c["ls"], c["mc"], c["form"]
            lo, hi = max(l["IC95"][0], 0.), min(l["IC95"][1], 1.)
            if k == 1 and l["P_f"] > 0:
                lo = max(lo, l["P_f"] / 10) if lo <= 0 else lo
            ax.errorbar([I - 1.2], [l["P_f"]], yerr=[[l["P_f"] - lo], [hi - l["P_f"]]], color=cor, marker="o", ms=4.5,
                        capsize=2, lw=1.0, ls="none")
            if m["n_f"] > 0:
                ax.errorbar([I + 1.2], [m["P_f"]], yerr=[[m["P_f"] - m["IC95"][0]], [m["IC95"][1] - m["P_f"]]],
                            color=cor, marker="s", ms=4.5, mfc="white", capsize=2, lw=1.0, ls="none")
            else:  # nenhuma falha: so o limite superior de Wilson
                ax.plot([I + 1.2], [m["IC95"][1]], marker="v", color=cor, ms=5, mfc="white", ls="none")
            ax.plot([I], [f["P_f"]], marker="x", color=cor, ms=5, mew=1.3, ls="none")
            ax.plot([I], [c["pf_artigo"]], marker="o", ms=7, mfc="none", mec=cor, mew=0.8, ls="none")
        ax.set_xlabel("intensidade da chuva (mm/h)")
        ax.set_xticks(range(0, 101, 10))
        ax.set_xlim(0, 100)
        if k == 0:
            ax.set_ylabel(r"$P(FS < 1)$")
            ax.set_ylim(-0.02, 1.02)
        else:
            ax.set_yscale("log")
            ax.set_ylim(1e-30, 5)
            ax.set_yticks([1e-30, 1e-25, 1e-20, 1e-15, 1e-10, 1e-5, 1])
        painel(ax, "ab"[k])
    hs = [Line2D([], [], color=COR_D[D], ls=(0, (5, 2))) for D in (12, 24, 36, 48)]
    rot = ["Hu, chuva de %d h" % D for D in (12, 24, 36, 48)]
    hs += [Line2D([], [], color=TINTA2, marker="o", ls="none", ms=4.5), Line2D([], [], color=TINTA2, marker="s",
                                                                              mfc="white", ls="none", ms=4.5),
           Line2D([], [], color=TINTA2, marker="v", mfc="white", ls="none", ms=5),
           Line2D([], [], color=TINTA2, marker="x", ls="none", ms=5, mew=1.3),
           Line2D([], [], color=TINTA2, marker="o", mfc="none", ls="none", ms=7, mew=0.8)]
    rot += ["LS (IC 95 %)", "MC (IC de Wilson)", "MC sem falha: limite de Wilson", "FORM", "Hu no cenário"]
    axs[0].legend(hs, rot, loc="upper left", fontsize=6.0)
    salva(fig, "fig_fragilidade.pdf")


def figura_convergencia(por):
    ids = ["I30_D24", "I20_D36"]
    fig, axs = plt.subplots(1, 2, figsize=(7.0, 3.1), constrained_layout=True)
    for k, (ax, i) in enumerate(zip(axs, ids)):
        c = por[i]
        lin = c["linhas_csv"]
        p = np.array([x[3] for x in lin])
        an = c["form"]["analises"] + np.cumsum([x[2] for x in lin])
        n = np.arange(1, len(p) + 1)
        mu = np.cumsum(p) / n
        var = np.array([np.var(p[:j], ddof=1) if j > 1 else np.nan for j in n])
        se = np.sqrt(var / n)
        ax.fill_between(an, np.maximum(mu - 1.96 * se, 0), mu + 1.96 * se, color=PALETA[0], alpha=0.22, lw=0)
        ax.plot(an, mu, color=PALETA[0], lw=1.3)
        fl = (c["mc_min"] < 1.).astype(float)
        nm = np.arange(1, len(fl) + 1)
        km = np.cumsum(fl)
        pm = km / nm
        w = np.array([wilson(int(a), int(b)) for a, b in zip(km, nm)])
        ax.fill_between(nm, w[:, 0], w[:, 1], color=PALETA[1], alpha=0.18, lw=0)
        ax.plot(nm, pm, color=PALETA[1], lw=1.1)
        ax.axhline(c["ls"]["P_f"], color=PALETA[0], lw=0.7, ls=":")
        ax.axhline(c["form"]["P_f"], color=TINTA2, lw=0.8, ls="-.")
        ax.set_xscale("log")
        ax.set_xlim(8, max(nm[-1], an[-1]) * 1.1)
        top = {"I30_D24": 0.2, "I20_D36": 0.012}[i]
        ax.set_ylim(0, top)
        ax.set_xlabel("análises do MEF")
        ax.set_ylabel(r"$P_f$")
        ax.set_title(nome_caso(i), fontsize=8.5)
        painel(ax, "ab"[k])
    axs[0].legend([Line2D([], [], color=PALETA[0]), Patch(color=PALETA[0], alpha=0.3), Line2D([], [], color=PALETA[1]),
                   Patch(color=PALETA[1], alpha=0.25), Line2D([], [], color=TINTA2, ls="-.")],
                  ["LS (FORM + linhas)", "LS, IC 95 %", "MC direto", "MC, IC de Wilson", r"FORM, $\Phi(-\beta)$"],
                  loc="upper right", fontsize=6.3)
    salva(fig, "fig_convergencia.pdf")


def figura_raizes(por):
    dados = [[x[1] for x in por[i]["linhas_csv"]] for i in ORDEM]
    betas = [por[i]["form"]["beta"] for i in ORDEM]
    bls = [-stats.norm.ppf(por[i]["ls"]["P_f"]) for i in ORDEM]
    fig, ax = plt.subplots(figsize=(6.4, 3.1), constrained_layout=True)
    pos = np.arange(len(ORDEM))
    ax.boxplot(dados, positions=pos, widths=0.5, showfliers=True, flierprops={"markersize": 2},
               medianprops={"color": TINTA}, boxprops={"color": TINTA2}, whiskerprops={"color": TINTA2},
               capprops={"color": TINTA2})
    ax.plot(pos, betas, "D", color=PALETA[1], ms=5, label=r"$\beta$ do FORM")
    ax.plot(pos, bls, "x", color=PALETA[0], ms=7, mew=1.5, label=r"$-\Phi^{-1}(P_f^{LS})$")
    ax.axhline(0, color=TINTA2, lw=0.8)
    ax.set_xticks(pos)
    ax.set_xticklabels([nome_caso(i).replace(", ", "\n") for i in ORDEM], fontsize=7)
    ax.set_ylabel(r"raiz $c$ de cada linha")
    ax.legend(loc="upper right")
    salva(fig, "fig_raizes.pdf")


def figura_hist(por, hu, hom):
    f11 = hu["fig11"]
    fig, axs = plt.subplots(1, 2, figsize=(7.0, 3.0), constrained_layout=True)
    for k, (ax, i) in enumerate(zip(axs, ["I30_D24", "I50_D24"])):
        c = por[i]
        x = c["mc_fsD"]
        I = int(c["I"])
        muh, sdh = f11[I]["media_t24"], f11[I]["dp_t24"]
        xx = np.linspace(0.55, 1.45, 400)
        ax.axvspan(0.5, 1., color=GRADE, alpha=0.45, lw=0)
        ax.hist(x, bins=30, density=True, color=PALETA[1], alpha=0.35, edgecolor="white", linewidth=0.4)
        ax.plot(xx, stats.norm.pdf(xx, x.mean(), x.std(ddof=1)), color=PALETA[1], lw=1.4)
        ax.plot(xx, stats.norm.pdf(xx, muh, sdh), color=HU_COR, lw=1.4, ls=(0, (5, 2)))
        if i == "I50_D24":
            ax.plot(xx, stats.gaussian_kde(hom[24]["amostra"])(xx), color=PALETA[2], lw=1.2, ls=":")
        ax.set_xlim(0.55, 1.45)
        ax.set_xlabel("FS em t = 24 h (fim da chuva)")
        ax.set_ylabel("densidade")
        ax.set_title("%s (%d realizações)" % (nome_caso(i), len(x)), fontsize=8.5)
        painel(ax, "ab"[k])
    axs[1].legend([Patch(color=PALETA[1], alpha=0.4), Line2D([], [], color=PALETA[1]),
                   Line2D([], [], color=HU_COR, ls=(0, (5, 2))), Line2D([], [], color=PALETA[2], ls=":"),
                   Patch(color=GRADE, alpha=0.6)],
                  ["nosso MC", "normal com a nossa média e desvio", "normal com a média e o desvio de Hu (Fig. 11)",
                   "limite homogêneo", "FS < 1"], loc="upper right", fontsize=6.0)
    salva(fig, "fig_hist.pdf")


# ----------------------------------------------------------------------------------------------------------- numeros
def numeros(casos, por, pc, hu, hom):
    d = {}
    f9 = hu["fig9"]
    c50 = por["I50_D24"]
    d["nCasos"] = str(len(casos))
    d["nTot"] = str(sum(c["form"]["analises"] + c["ls"]["analises"] + c["mc"]["S"] + c["mc13"]["S"] for c in casos))
    d["nFORMtot"] = str(sum(c["form"]["analises"] for c in casos))
    d["nLStot"] = str(sum(c["ls"]["analises"] for c in casos))
    d["nMCtot"] = str(sum(c["mc"]["S"] for c in casos))
    d["nTreze"] = str(sum(c["mc13"]["S"] for c in casos))
    d["nFormConv"] = str(sum(1 for c in casos if c["form"]["convergiu"]))
    d["anPorLinha"] = br(sum(c["ls"]["analises"] for c in casos) / sum(c["ls"]["linhas"] for c in casos), 2)
    d["fsDetZero"] = br(c50["det"]["FS"]["0.0"], 3)
    d["fsDetVinte"] = br(c50["det"]["FS"]["24.0"], 3)
    d["fsDetSetenta"] = br(c50["det"]["FS"]["72.0"], 3)
    hu0, hu24 = f9[0], f9[np.argmin(abs(f9[:, 0] - 24))]
    hu72 = f9[np.argmin(abs(f9[:, 0] - 72))]
    d["huZero"], d["huVinte"], d["huSetenta"] = br(hu0[1], 3), br(hu24[1], 3), br(hu72[1], 3)
    d["huSdZero"], d["huSdVinte"], d["huSdSetenta"] = br(hu0[2], 3), br(hu24[2], 3), br(hu72[2], 3)
    d["difZero"] = br(100 * (c50["det"]["FS"]["0.0"] / hu0[1] - 1), 1)
    d["difVinte"] = br(100 * (c50["det"]["FS"]["24.0"] / hu24[1] - 1), 1)
    d["difSetenta"] = br(100 * (c50["det"]["FS"]["72.0"] / hu72[1] - 1), 1)
    m = pc["malha"]
    d["ksBaixoVinte"] = br(m[("5.2e-6", "tri3")][24.0], 3)
    d["triDoisVinte"], d["triTresVinte"], d["triQuatroVinte"] = (br(m[("2e-5", x)][24.0], 3)
                                                                 for x in ("tri2", "tri3", "tri4"))
    d["triTresIXII"], d["triQuatroIXII"] = br(pc["I30D12"]["tri3"], 3), br(pc["I30D12"]["tri4"], 3)
    d["vertical"] = br(pc["vertical"], 3)
    d["inclinadaTri"] = br(m[("2e-5", "tri3")][24.0], 3)
    d["fsNoventaDois"], d["fsNoventaTres"] = br(pc["fs24"][("2e-5", "tri2", 90)], 3), br(pc["fs24"][("2e-5", "tri3", 90)], 3)
    pi = c50["mc13"]["por_instante"]
    d["sdZero"] = br(pi["0.0"][1], 3)
    d["muZero"] = br(pi["0.0"][0], 3)
    d["razaoSdZero"] = br(pi["0.0"][1] / hu0[2], 2)
    d["sdVinteMC"] = br(c50["mc"]["dp_D"], 3)
    d["razaoSdVinte"] = br(c50["mc"]["dp_D"] / hu24[2], 2)
    d["cvHuZero"] = br(100 * hu0[2] / hu0[1], 1)
    d["cvNossoZero"] = br(100 * pi["0.0"][1] / pi["0.0"][0], 1)
    d["cvHomZero"] = br(100 * hom[0]["cov"], 1)
    d["cvHomVinte"] = br(100 * hom[24]["cov"], 1)
    d["cvHuVinte"] = br(100 * hu24[2] / hu24[1], 1)
    d["wcZero"], d["wcVinte"] = br(hom[0]["w_c"], 2), br(hom[24]["w_c"], 2)
    d["pfHomVinte"] = pct(hom[24]["pf"], 1)
    d["homFora"] = pct(hom["fora"], 2)
    # a aproximacao normal: Hu a 50 e 30 mm/h e nos
    f11 = hu["fig11"]
    d["huNormalCinquenta"] = pct(stats.norm.cdf((1 - f11[50]["media_t24"]) / f11[50]["dp_t24"]), 1)
    d["huNormalTrinta"] = pct(stats.norm.cdf((1 - f11[30]["media_t24"]) / f11[30]["dp_t24"]), 1)
    c30 = por["I30_D24"]
    d["nossoNormalTrinta"] = pct(stats.norm.cdf((1 - c30["mc"]["media_D"]) / c30["mc"]["dp_D"]), 2)
    d["pfLSTrinta"], d["cvLSTrinta"] = pct(c30["ls"]["P_f"], 2), pct(c30["ls"]["CoV"], 1)
    d["pfMCTrinta"] = pct(c30["mc"]["P_f"], 2)
    d["nfMCTrinta"], d["nMCTrinta"] = str(c30["mc"]["n_f"]), str(c30["mc"]["S"])
    d["linhasTrinta"], d["anLSTrinta"] = str(c30["ls"]["linhas"]), str(c30["ls"]["analises"])
    d["zTrinta"] = ("%+.1f" % c30["z"]).replace(".", "{,}")
    c36 = por["I20_D36"]
    d["pfLSTrintaSeis"], d["cvLSTrintaSeis"] = pct(c36["ls"]["P_f"], 3), pct(c36["ls"]["CoV"], 1)
    d["pfMCTrintaSeis"] = pct(c36["mc"]["P_f"], 3)
    d["nfMCTrintaSeis"], d["nMCTrintaSeis"] = str(c36["mc"]["n_f"]), str(c36["mc"]["S"])
    d["linhasTrintaSeis"], d["anLSTrintaSeis"] = str(c36["ls"]["linhas"]), str(c36["ls"]["analises"])
    d["zTrintaSeis"] = ("%+.1f" % c36["z"]).replace(".", "{,}")
    # o MC direto que daria o CoV do LS em I20_D36: N = (1 - p) / (p CoV^2)
    p, cv = c36["ls"]["P_f"], c36["ls"]["CoV"]
    d["nMCEquivTrintaSeis"] = "{:,}".format(int(round((1 - p) / (p * cv * cv), -3))).replace(",", ".")
    d["razaoTrintaSeis"] = str(int(round((1 - p) / (p * cv * cv) / (c36["ls"]["analises"] + c36["form"]["analises"]))))
    p, cv = c30["ls"]["P_f"], c30["ls"]["CoV"]
    nmc = (1 - p) / (p * cv * cv)
    d["nMCEquivTrinta"] = "{:,}".format(int(round(nmc, -2))).replace(",", ".")
    d["razaoTrinta"] = str(int(round(nmc / (c30["ls"]["analises"] + c30["form"]["analises"]))))
    c40 = por["I40_D12"]
    d["pfLSQuarenta"], d["cvLSQuarenta"] = pct(c40["ls"]["P_f"], 3), pct(c40["ls"]["CoV"], 1)
    d["pfMCQuarenta"], d["nfMCQuarenta"] = pct(c40["mc"]["P_f"], 3), str(c40["mc"]["n_f"])
    d["pfFORMQuarenta"] = pct(c40["form"]["P_f"], 3)
    d["pfLSCinquenta"], d["pfMCCinquenta"] = pct(c50["ls"]["P_f"], 2), pct(c50["mc"]["P_f"], 2)
    zs = [c["z"] for c in casos if c["z"] is not None]
    d["zMax"] = br(max(abs(z) for z in zs), 1)
    d["nZ"] = str(len(zs))
    d["mcTrezeOk"] = str(sum(1 for c in casos if c["mc13"]["argmin_D"] == c["mc13"]["S"]))
    d["mcTrezeN"] = str(casos[0]["mc13"]["S"])
    d["mcTrezeDif"] = br(max(c["mc13"]["dif_max"] for c in casos), 4)
    # o nosso CoV do FS no fim da chuva e o de Hu
    d["cvNossoMin"] = br(min(100 * c["mc"]["dp_D"] / c["mc"]["media_D"] for c in casos), 1)
    d["cvNossoMax"] = br(max(100 * c["mc"]["dp_D"] / c["mc"]["media_D"] for c in casos), 1)
    d["linhaTeto"] = br(max(x[1] for x in por["I30_D12"]["linhas_csv"]), 0)
    linhas = ["% gerado por gera_dados.py; nao edite"]
    linhas += [r"\newcommand{\%s}{%s}" % (k, v) for k, v in d.items()]
    return "\n".join(linhas) + "\n", d


def main():
    os.makedirs(FIGS, exist_ok=True)
    os.makedirs(TAB, exist_ok=True)
    casos, por = carrega_campanha()
    pc, hu = carrega_pc(), carrega_hu()
    hom = limite_homogeneo(pc)
    blocos = {"pf": tabela_pf(casos, por), "custo": tabela_custo(por), "det": tabela_det(pc, por, hu),
              "normal": tabela_normal(por, hu), "desvio": tabela_desvio(por, hu, hom)}
    for nome, texto in blocos.items():
        with open(os.path.join(TAB, nome + ".tex"), "w") as f:
            f.write("% gerado por gera_dados.py; nao edite\n" + texto)
    texto, d = numeros(casos, por, pc, hu, hom)
    with open(os.path.join(AQUI, "numeros.tex"), "w") as f:
        f.write(texto)
    figura_fs_t(por)
    figura_fig9(por, pc, hu)
    figura_malha(pc)
    figura_fig11(por, pc, hu)
    figura_fragilidade(por, hu)
    figura_convergencia(por)
    figura_raizes(por)
    figura_hist(por, hu, hom)
    print(f"{len(casos)} cenarios; {d['nTot']} analises do MEF; desvio do FS em 0 h {d['sdZero']} contra "
          f"{d['huSdZero']} de Hu (razao {d['razaoSdZero']}); CoV homogeneo {d['cvHomZero']} %")


if __name__ == "__main__":
    main()
