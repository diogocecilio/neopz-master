#!/usr/bin/env python3
"""gera_dados.py: as tabelas, os numeros e as figuras do relatorio_liu2024.tex (o caso B, Liu, Han, Huang, Huang e Cheng
2024, Comput. Geotech. 165:105938), a partir de:

  - <campanha>/resumo.json (resumo_liu2024.py, refeito pelo gera.sh): o deterministico, o FORM modal, o line sampling e
    o Monte Carlo direto dos cenarios 2 (c', phi') e 3 (c', phi', k_s);
  - <campanha>/<S>/ls/linhas.csv, <S>/mc/mc_*.csv, <S>/form/ e S2/pressoes.txt: as linhas do LS (o estimador, a
    convergencia e as raizes; refeitos com as linhas que existirem, entao uma campanha estendida entra sem mudar nada),
    as realizacoes do MC (o FS, o nivel d'agua na crista e na meia altura), o FS no centro do FORM e o p_w permanente;
  - ~/projects/runs/liu2024/det (o PC): o FS deterministico com as malhas, a chuva vertical, o psi_aev e sem succao;
  - dados/liu_fig3_pwp.csv e dados/liu_fig11.csv: a Fig. 3 (p_w nas secoes X-X e Y-Y) e a Fig. 11 (o desvio do FS)
    de Liu, lidas dos pixels das figuras raster do PDF (ver o cabecalho de cada arquivo).

Grava tab/*.tex, numeros.tex e figs/*.pdf. Roda em segundos.

    python3 -I gera_dados.py [pasta da campanha = ~/projects/runs/fermi_liu2024]"""

import csv
import glob
import json
import math
import os
import re
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
from matplotlib.path import Path  # noqa: E402
from matplotlib.tri import Triangulation, LinearTriInterpolator  # noqa: E402
from scipy import optimize, stats  # noqa: E402

AQUI = os.path.dirname(os.path.abspath(__file__))
RAIZ = os.path.abspath(os.path.join(AQUI, "..", "..", ".."))
sys.path.insert(0, os.path.join(RAIZ, "Projects", "SlopeRandomFieldRefinement", "cluster"))
import resumo_hu2024 as rh  # noqa: E402

CAMPANHA = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/projects/runs/fermi_liu2024")
PC = os.path.expanduser("~/projects/runs/liu2024")
DADOS = os.path.join(AQUI, "dados")
FIGS = os.path.join(AQUI, "figs")
TAB = os.path.join(AQUI, "tab")
CEN = ["S2", "S3"]
NOME = {"S2": "cenário 2", "S3": "cenário 3"}
TALUDE = [(0, 0), (30, 0), (30, 5), (20, 5), (10, 10), (0, 10)]

# o artigo
FS_BISHOP = 1.525                                    # Fig. 2, Bishop com "Entry and Exit"
TAB4 = {"S2": (4.08e-5, 3.95), "S3": (1.19e-4, 3.68)}  # alfa = 0, q/k_s = 0,01: (P_f, beta) das lognormais ajustadas
TAB2_LHS, TAB2_CDF = 4.0e-3, 4.58e-3                 # a verificacao (o modelo do cenario 3): LHS com 2500 e a CDF ajustada
N_LHS = 2500

# o estilo das figuras dos outros relatorios do plano (hu2024/gera_dados.py)
PALETA = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
TINTA, TINTA2, GRADE = "#0b0b0b", "#52514e", "#dddcd8"
plt.rcParams.update({
    "font.size": 8.5, "axes.labelsize": 8.5, "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "legend.fontsize": 7,
    "axes.edgecolor": TINTA2, "axes.labelcolor": TINTA, "xtick.color": TINTA2, "ytick.color": TINTA2,
    "axes.grid": True, "grid.color": GRADE, "grid.linewidth": 0.5, "lines.linewidth": 1.4, "legend.frameon": False,
    "savefig.dpi": 200, "font.family": "DejaVu Sans", "mathtext.fontset": "dejavusans",
})
COR = {"S2": PALETA[0], "S3": PALETA[1]}
LIU_COR = TINTA2


# ---------------------------------------------------------------------------------------------------------- formatos
def br(x, casas=2):
    """numero com virgula decimal para o LaTeX"""
    if x is None or (isinstance(x, float) and not math.isfinite(x)):
        return "--"
    t = f"{x:.{casas}f}".replace(".", "{,}")
    return r"\ensuremath{-}" + t[1:] if t.startswith("-") else t


def cient(x, casas=1, dolar=True):
    """x em notacao cientifica para o LaTeX: 3{,}0 \\times 10^{-5}"""
    if x is None:
        return "--"
    if x == 0:
        return "0"
    e = int(math.floor(math.log10(abs(x))))
    m = x / 10 ** e
    if round(m, casas) >= 10:
        m, e = m / 10, e + 1
    t = f"{m:.{casas}f}".replace(".", "{,}") + r"\times 10^{" + str(e) + "}"
    return r"\ensuremath{" + t + "}" if dolar else t


def inteiro(n):
    """inteiro com o espaco fino dos milhares: 10\\,000"""
    s = f"{int(n):,}".replace(",", r"\,")
    return s


def painel(ax, t):
    ax.text(-0.02, 1.02, "(%s)" % t, transform=ax.transAxes, ha="right", va="bottom", fontsize=9, color=TINTA)


def salva(fig, nome):
    fig.savefig(os.path.join(FIGS, nome), metadata={"CreationDate": None, "ModDate": None})
    plt.close(fig)


# ------------------------------------------------------------------------------------------------------- estatistica
def zeta(cv):
    return math.sqrt(math.log(1. + cv * cv))


def beta_ln(mu, sd):
    """o beta da Eq. 13 de Liu: FS lognormal com media mu e desvio sd, falha FS < 1"""
    z = zeta(sd / mu)
    return (math.log(mu) - 0.5 * z * z) / z


def mu_de_beta(beta, sd):
    """a media do FS que, com o desvio sd, da o beta (Eq. 13)"""
    return optimize.brentq(lambda m: beta_ln(m, sd) - beta, 1.0001, 5.)


def phi(x):
    return float(stats.norm.cdf(x))


def ic_log(pf, cv, z=1.96):
    """IC de 95 % na escala log (delta), para os graficos log: pf exp(+-z s), s = sqrt(ln(1 + CoV^2))"""
    s = zeta(cv)
    return pf * math.exp(-z * s), pf * math.exp(z * s)


def ls_estatistica(p):
    """o estimador do LS (media das probabilidades das linhas), o CoV e o IC normal"""
    p = np.asarray(p, float)
    n = len(p)
    m = p.mean()
    se = p.std(ddof=1) / math.sqrt(n)
    return {"P_f": m, "CoV": se / m, "IC95": (m - 1.96 * se, m + 1.96 * se), "linhas": n,
            "beta_eq": -stats.norm.ppf(m)}


# ------------------------------------------------------------------------------------------------------------ dados
def le_mc(pasta):
    """{real: linha do csv} de todas as realizacoes do MC (t = 0: o permanente)"""
    out = {}
    for a in glob.glob(os.path.join(pasta, "mc", "mc_*_*.csv")):
        if os.path.basename(a).startswith("mc_mapa"):
            continue
        with open(a, newline="") as f:
            for r in csv.DictReader(f):
                if r.get("situacao", "ok") == "ok" and r.get("FS_min") not in (None, ""):
                    out[int(r["real"])] = r
    return out


def le_pressoes(arq):
    """os nos (x, y) e o p_w (kPa) do arquivo do grava_pressoes (um instante)"""
    with open(arq) as f:
        f.readline()
        f.readline()
        nn = int(f.readline().split()[1])
        xy = np.array([[float(v) for v in f.readline().split()] for _ in range(nn)])
        f.readline()
        pw = np.array([float(f.readline()) for _ in range(nn)])
    return xy, pw


def carrega_campanha():
    casos = json.load(open(os.path.join(CAMPANHA, "resumo.json")))
    por = {c["id"]: c for c in casos}
    for i in CEN:
        c, d = por[i], os.path.join(CAMPANHA, i)
        lin = []
        with open(os.path.join(d, "ls", "linhas.csv"), newline="") as f:
            for r in csv.DictReader(f):
                lin.append((int(r["linha"]), float(r["c"]), int(float(r["analises"])), float(r["p"])))
        c["linhas_csv"] = sorted(lin)
        c["lsx"] = ls_estatistica([x[3] for x in c["linhas_csv"]])
        c["form_json"] = json.load(open(os.path.join(d, "form", "form.json")))
        reais = le_mc(d)
        ns = sorted(reais)
        c["mc_n"] = np.array(ns)
        c["mc_fs"] = np.array([float(reais[n]["FS_min"]) for n in ns])
        c["mc_nivel_xx"] = np.array([float(reais[n]["nivel_x30"]) for n in ns])   # a crista, x = 10 m (X-X)
        c["mc_nivel_yy"] = np.array([float(reais[n]["nivel_x35"]) for n in ns])   # x = 15 m (Y-Y)
        c["mc_ymin"] = np.array([float(reais[n]["banda_ymin"]) for n in ns])     # o fundo da banda de plastificacao
        c["mc_s"] = np.array([float(reais[n]["s_fator"]) + float(reais[n]["s_transiente"]) for n in ns])
        # o FS no centro do FORM (xi = 0: as medianas das lognormais)
        with open(os.path.join(d, "form", "it0_centro_0_0.csv"), newline="") as f:
            c["fs_xi0"] = float(next(csv.DictReader(f))["FS_min"])
    xy, pw = le_pressoes(os.path.join(CAMPANHA, "S2", "pressoes.txt"))
    return casos, por, (xy, pw)


def fs_log(nome):
    """o FS da realizacao -1 (as medias) no log do PC"""
    t = open(os.path.join(PC, "det", nome + ".log")).read()
    return float(re.findall(r"realizacao\s+-1: FS 0 h: ([0-9.]+)", t)[-1])


def carrega_pc():
    return {"final": fs_log("detB_final"), "t4": fs_log("detB_t4"), "f2": fs_log("detB_f2"),
            "vert": fs_log("detB_aev2.6_vert"), "aev0": fs_log("detB_aev0.001_tri3"), "phib0": fs_log("detB_phib0")}


def carrega_liu():
    def le(nome):
        with open(os.path.join(DADOS, nome), newline="") as f:
            return list(csv.DictReader(l for l in f if not l.startswith("#")))
    f3 = {"XX": [], "YY": []}
    for r in le("liu_fig3_pwp.csv"):
        f3[r["secao"]].append((float(r["y_m"]), float(r["pw_kPa"])))
    f3 = {k: np.array(v) for k, v in f3.items()}
    f11 = {r["cenario"]: (float(r["dp_FS"]), float(r["FS_min_2000"])) for r in le("liu_fig11.csv")}
    return {"fig3": f3, "fig11": f11}


# ---------------------------------------------------------------------------------------------------------- perfis
def interpolador(xy, pw):
    tri = Triangulation(xy[:, 0], xy[:, 1])
    cx = xy[tri.triangles].mean(axis=1)
    tri.set_mask(~Path(TALUDE).contains_points(cx, radius=1e-6))
    return tri, LinearTriInterpolator(tri, pw)


def perfil(it, x, topo, n=401):
    ys = np.linspace(0., topo, n)
    v = np.ma.filled(it(np.full_like(ys, x), ys), np.nan)
    return ys, v


def nivel(ys, v):
    """a cota em que p_w passa de positivo a negativo (interpolacao linear)"""
    for k in range(len(v) - 1):
        if v[k] >= 0 > v[k + 1]:
            return ys[k] + (ys[k + 1] - ys[k]) * v[k] / (v[k] - v[k + 1])
    return None


def nivel_liu(a):
    return nivel(a[:, 0], a[:, 1])


# ----------------------------------------------------------------------------------------------------------- contas
def contas(por, pc, liu, agua):
    """tudo o que as tabelas, as figuras e o texto usam"""
    r = {}
    fsdet = por["S2"]["det"]["FS_min"]
    r["fsdet"] = fsdet
    r["razao_det"] = FS_BISHOP / fsdet
    for i in CEN:
        c = por[i]
        x = c["mc_fs"]
        mu, sd = float(x.mean()), float(x.std(ddof=1))
        e = {"mu": mu, "sd": sd, "cv": sd / mu, "N": len(x), "nf": int((x < 1).sum()), "min": float(x.min()),
             "assim": float(stats.skew(x))}
        e["wilson"] = rh.wilson(e["nf"], e["N"])
        e["b_ln"] = beta_ln(mu, sd)
        e["pf_ln"] = phi(-e["b_ln"])
        e["b_n"] = (mu - 1.) / sd
        e["pf_n"] = phi(-e["b_n"])
        # o nosso FS reduzido pela razao Bishop/SRM do deterministico (mesmo CoV)
        k = r["razao_det"]
        e["b_ln_red"] = beta_ln(k * mu, k * sd)
        e["pf_ln_red"] = phi(-e["b_ln_red"])
        e["pf_n_red"] = phi(-(k * mu - 1.) / (k * sd))
        # a reducao da media que leva o nosso beta (mesmo CoV) ao beta do artigo
        bl = TAB4[i][1]
        e["mu_req"] = math.exp(bl * zeta(e["cv"]) + 0.5 * zeta(e["cv"]) ** 2)
        e["red_req"] = 1. - e["mu_req"] / mu
        e["dbeta_dlnmu"] = 1. / zeta(e["cv"])
        # o artigo: desvio da Fig. 11 e a media que a Eq. 13 exige para o beta da Tabela 4
        sdl, minl = liu["fig11"][i]
        e["liu_sd"], e["liu_min"] = sdl, minl
        e["liu_mu"] = mu_de_beta(bl, sdl)
        e["liu_cv"] = sdl / e["liu_mu"]
        e["liu_pf_n"] = phi(-(e["liu_mu"] - 1.) / sdl)
        # LS
        ls = c["lsx"]
        e["ls"] = ls
        e["razao_ln_ls"] = e["pf_ln"] / ls["P_f"]
        e["p0_tab2"] = (1. - TAB2_LHS) ** e["N"]
        e["p0_tab4"] = (1. - TAB4[i][0]) ** e["N"]
        e["n_tab2"] = TAB2_LHS * e["N"]
        # as linhas: a raiz c e a contribuicao das maiores
        p = np.array([y[3] for y in c["linhas_csv"]])
        cs = np.array([y[1] for y in c["linhas_csv"]])
        ps = np.sort(p)[::-1]
        e["frac_maior"] = float(ps[0] / p.sum())
        e["frac_10"] = float(ps[:10].sum() / p.sum())
        e["n_c4"] = int((cs < 4).sum())
        e["frac_c4"] = float(p[cs < 4].sum() / p.sum())
        e["frac_base"] = float((c["mc_ymin"] < 1.).mean())
        e["c_med"], e["c_min"], e["c_max"] = float(np.median(cs)), float(cs.min()), float(cs.max())
        e["c_p05"], e["c_p95"] = float(np.percentile(cs, 5)), float(np.percentile(cs, 95))
        e["an_linha"] = float(np.mean([y[2] for y in c["linhas_csv"]]))
        e["an_ls"] = int(sum(y[2] for y in c["linhas_csv"]))
        e["form"] = c["form"]
        e["fs_xi0"] = c["fs_xi0"]
        e["s_med"] = float(np.median(c["mc_s"]))
        # o nivel d'agua nas secoes
        for s, v in (("xx", c["mc_nivel_xx"]), ("yy", c["mc_nivel_yy"])):
            e["niv_" + s] = (float(v.mean()), float(v.std(ddof=1)), float(np.percentile(v, 2.5)),
                             float(np.percentile(v, 97.5)), float(v.min()), float(v.max()))
        r[i] = e
    # o efeito do campo de k: a razao das P_f e o MC pareado (a realizacao n tem os mesmos c' e phi' nos dois)
    a, b = r["S2"]["ls"], r["S3"]["ls"]
    r["razao_ls"] = b["P_f"] / a["P_f"]
    r["razao_ls_cv"] = math.sqrt(a["CoV"] ** 2 + b["CoV"] ** 2)
    r["razao_tab4"] = TAB4["S3"][0] / TAB4["S2"][0]
    n2 = dict(zip(por["S2"]["mc_n"], por["S2"]["mc_fs"]))
    n3 = dict(zip(por["S3"]["mc_n"], por["S3"]["mc_fs"]))
    comum = sorted(set(n2) & set(n3))
    f2, f3 = np.array([n2[n] for n in comum]), np.array([n3[n] for n in comum])
    d = f3 - f2
    r["par"] = {"N": len(comum), "f2": f2, "f3": f3, "d_mu": float(d.mean()), "d_sd": float(d.std(ddof=1)),
                "corr": float(np.corrcoef(f2, f3)[0, 1]), "sd2": float(f2.std(ddof=1)), "sd3": float(f3.std(ddof=1)),
                "mu2": float(f2.mean()), "mu3": float(f3.mean()),
                "d_p": (float(np.percentile(d, 2.5)), float(np.percentile(d, 97.5)))}
    # o deterministico: os perfis e os niveis contra a Fig. 3
    tri, it = interpolador(*agua)
    r["tri"], r["it"] = tri, it
    for s, x, topo in (("XX", 10., 10.), ("YY", 15., 7.5)):
        ys, v = perfil(it, x, topo)
        r["perfil_" + s] = (ys, v)
        r["nivel_" + s] = nivel(ys, v)
        r["nivel_liu_" + s] = nivel_liu(liu["fig3"][s])
        lf = liu["fig3"][s]
        nosso = np.interp(lf[:, 0], ys, v)
        r["dif_max_" + s] = float(np.max(np.abs(nosso - lf[:, 1])))
        r["succ_" + s] = float(v[-1])
        r["succ_liu_" + s] = float(lf[-1, 1])
        sat = lf[:, 1] > 0
        r["dif_sat_" + s] = float(np.max(np.abs(nosso[sat] - lf[sat, 1])))
    r["pc"] = pc
    return r


# ----------------------------------------------------------------------------------------------------------- tabelas
def tabela_det(r, por):
    pc = r["pc"]
    lin = [r"\begin{table}[htbp]\centering\footnotesize",
           r"\caption{FS determinístico do caso B (as médias: $c' = 10$~kPa, $\varphi' = 15^\circ$, $k_s = 10^{-5}$~m/s), "
           r"contra o Bishop de Liu (Fig.~2). Permanente: a malha da percolação; fator: a malha do SRM e os níveis "
           r"adaptativos. A linha da campanha é a do fermi; as outras, do PC.}\label{tab:det}",
           r"\begin{tabular}{lllr}\toprule",
           r"variante & permanente & fator & FS \\\midrule",
           r"\textbf{campanha} (a base) & tri3 & tri1 + 3 & \textbf{%s} \\" % br(r["fsdet"], 4),
           r"mais um nível no fator & tri3 & tri1 + 4 & %s \\" % br(pc["final"], 4),
           r"permanente mais fino & tri4 & tri1 + 4 & %s \\" % br(pc["t4"], 4),
           r"fator mais fino & tri3 & tri2 + 4 & %s \\" % br(pc["f2"], 4),
           r"chuva vertical (projeção horizontal) & tri3 & tri1 + 4 & %s \\" % br(pc["vert"], 4),
           r"$\psi_{aev} = 0{,}001$~kPa na integral de Fredlund & tri3 & tri1 + 4 & %s \\" % br(pc["aev0"], 4),
           r"sem a parcela da sucção ($\varphi^b = 0$) & tri3 & tri1 + 4 & %s \\" % br(pc["phib0"], 4),
           r"\midrule",
           r"no centro do FORM ($\xi = 0$: medianas), cenário 2 & tri3 & tri1 + 3 & %s \\" % br(por["S2"]["fs_xi0"], 4),
           r"no centro do FORM ($\xi = 0$: medianas), cenário 3 & tri3 & tri1 + 3 & %s \\" % br(por["S3"]["fs_xi0"], 4),
           r"\midrule Liu et al.\ (2024): Bishop, ``Entry and Exit'' & 0,5~m & -- & %s \\" % br(FS_BISHOP, 3),
           r"\bottomrule\end{tabular}\end{table}"]
    return "\n".join(lin) + "\n"


def tabela_pf(r):
    s2, s3 = r["S2"], r["S3"]

    def ls_txt(e):
        ls = e["ls"]
        return r"%s ($CoV$ %s\,\%%)" % (cient(ls["P_f"]), br(100 * ls["CoV"], 0))

    def ic_txt(e):
        lo, hi = ls_ic(e)
        return r"[%s; %s]" % (cient(lo), cient(hi))
    lin = [r"\begin{table}[htbp]\centering\footnotesize\setlength{\tabcolsep}{4pt}",
           r"\caption{$P_f = P(FS < 1)$ do caso B ($\alpha = 0^\circ$, $q/k_s = 0{,}01$) por método, contra as Tabelas 2 e "
           r"4 de Liu. LS: IC de 95\,\% na escala log ($P_f\, e^{\pm 1{,}96 s}$, $s^2 = \ln(1 + CoV^2)$; o normal "
           r"tem limite inferior negativo com $CoV > 51\,\%$). MC: IC de Wilson. Ajustes: a Eq.~13 de Liu (lognormal) e "
           r"a normal com a média e o desvio do nosso MC. A Tabela 2 é a verificação no modelo do cenário 3.}"
           r"\label{tab:pf}",
           r"\begin{tabular}{llll}\toprule",
           r"método & detalhe & cenário 2 ($c'$, $\varphi'$) & cenário 3 ($c'$, $\varphi'$, $k_s$) \\\midrule",
           r"FORM modal & $\beta$; $\Phi(-\beta)$ & %s; %s & %s; %s \\" % (
               br(s2["form"]["beta"], 2), cient(s2["form"]["P_f"]), br(s3["form"]["beta"], 2), cient(s3["form"]["P_f"])),
           r"line sampling & %s e %s linhas & %s & %s \\" % (inteiro(s2["ls"]["linhas"]), inteiro(s3["ls"]["linhas"]),
                                                          ls_txt(s2), ls_txt(s3)),
           r" & IC 95\,\%% (log) & %s & %s \\" % (ic_txt(s2), ic_txt(s3)),
           r" & $-\Phi^{-1}(P_f)$ & %s & %s \\" % (br(s2["ls"]["beta_eq"], 2), br(s3["ls"]["beta_eq"], 2)),
           r"MC direto & falhas & %d em %s & %d em %s \\" % (s2["nf"], inteiro(s2["N"]), s3["nf"], inteiro(s3["N"])),
           r" & IC 95\,\%% (Wilson) & [0; %s] & [0; %s] \\" % (cient(s2["wilson"][1]), cient(s3["wilson"][1])),
           r" & média $\pm$ desvio do FS & %s $\pm$ %s & %s $\pm$ %s \\" % (br(s2["mu"], 3), br(s2["sd"], 3),
                                                                           br(s3["mu"], 3), br(s3["sd"], 3)),
           r"ajuste ao nosso MC & lognormal: $\beta$; $P_f$ & %s; %s & %s; %s \\" % (
               br(s2["b_ln"], 2), cient(s2["pf_ln"]), br(s3["b_ln"], 2), cient(s3["pf_ln"])),
           r" & normal: $P_f$ & %s & %s \\" % (cient(s2["pf_n"]), cient(s3["pf_n"])),
           r"\midrule Liu, Tabela 4 & lognormal (2000): $\beta$; $P_f$ & %s; %s & %s; %s \\" % (
               br(TAB4["S2"][1], 2), cient(TAB4["S2"][0], 2), br(TAB4["S3"][1], 2), cient(TAB4["S3"][0], 2)),
           r"Liu, Tabela 2 & LHS, %s amostras & -- & %s \\" % (inteiro(N_LHS), cient(TAB2_LHS)),
           r" & CDF ajustada & -- & %s \\" % cient(TAB2_CDF, 2),
           r"\bottomrule\end{tabular}\end{table}"]
    return "\n".join(lin) + "\n"


def ls_ic(e):
    return ic_log(e["ls"]["P_f"], e["ls"]["CoV"])


def tabela_aprox(r):
    lin = [r"\begin{table}[htbp]\centering\footnotesize\setlength{\tabcolsep}{4pt}",
           r"\caption{As aproximações lognormal (Eq.~13 de Liu) e normal a partir da média e do desvio do FS. "
           r"``Reduzido'': o nosso FS multiplicado por $1{,}525/%s$ (Bishop/SRM no determinístico), com o mesmo $CoV$. "
           r"Liu: o desvio da Fig.~11 e a média que a Eq.~13 exige para o $\beta$ da Tabela 4. Última coluna: a "
           r"$P_f$ de referência de cada linha (o nosso LS; a Tabela 4).}\label{tab:aprox}" % br(r["fsdet"], 4),
           r"\begin{tabular}{llrrrrrrr}\toprule",
           r"cenário & FS & média & desvio & $CoV$ (\%) & $\beta_{LN}$ & $P_f$ lognormal & $P_f$ normal & referência \\"
           r"\midrule"]
    for i in CEN:
        e = r[i]
        k = r["razao_det"]
        lin.append(r"%s & nosso MC & %s & %s & %s & %s & %s & %s & LS %s \\" % (
            NOME[i], br(e["mu"], 3), br(e["sd"], 3), br(100 * e["cv"], 1), br(e["b_ln"], 2), cient(e["pf_ln"]),
            cient(e["pf_n"]), cient(e["ls"]["P_f"])))
        lin.append(r" & reduzido & %s & %s & %s & %s & %s & %s & \\" % (
            br(k * e["mu"], 3), br(k * e["sd"], 3), br(100 * e["cv"], 1), br(e["b_ln_red"], 2), cient(e["pf_ln_red"]),
            cient(e["pf_n_red"])))
        lin.append(r" & Liu & %s & %s & %s & %s & %s & %s & Tab.~4: %s \\" % (
            br(e["liu_mu"], 3), br(e["liu_sd"], 3), br(100 * e["liu_cv"], 1), br(TAB4[i][1], 2),
            cient(phi(-TAB4[i][1])), cient(e["liu_pf_n"]), cient(TAB4[i][0], 2)))
        if i == "S2":
            lin.append(r"\midrule")
    lin += [r"\bottomrule\end{tabular}\end{table}"]
    return "\n".join(lin) + "\n"


def tabela_custo(r, por):
    lin = [r"\begin{table}[htbp]\centering\footnotesize",
           r"\caption{O custo da campanha \texttt{liu2024}, em análises do MEF (uma análise: o permanente, lido no "
           r"cenário 2 e resolvido no 3, e o SRM). Tempo: a mediana por análise do MC num núcleo do fermi.}"
           r"\label{tab:custo}",
           r"\begin{tabular}{lrrrrrr}\toprule",
           r"cenário & FORM & linhas & LS & análises/linha & MC & tempo (s) \\\midrule"]
    for i in CEN:
        e = r[i]
        lin.append(r"%s & %d & %s & %s & %s & %s & %s \\" % (
            NOME[i], e["form"]["analises"], inteiro(e["ls"]["linhas"]), inteiro(e["an_ls"]), br(e["an_linha"], 2),
            inteiro(e["N"]), br(e["s_med"], 0)))
    tot = [sum(r[i]["form"]["analises"] for i in CEN), sum(r[i]["ls"]["linhas"] for i in CEN),
           sum(r[i]["an_ls"] for i in CEN), sum(r[i]["N"] for i in CEN)]
    lin.append(r"\midrule total & %d & %s & %s & & %s & \\" % (tot[0], inteiro(tot[1]), inteiro(tot[2]), inteiro(tot[3])))
    lin += [r"\bottomrule\end{tabular}\end{table}"]
    return "\n".join(lin) + "\n"


# ---------------------------------------------------------------------------------------------------------- figuras
def figura_pressoes(r, liu):
    tri, it = r["tri"], r["it"]
    fig = plt.figure(figsize=(7.0, 4.6), constrained_layout=True)
    g = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.25])
    ax = fig.add_subplot(g[0, :])
    z = np.ma.filled(it(tri.x, tri.y), np.nan)
    niveis = np.arange(-15, 61, 5)
    cs = ax.tricontourf(tri, z, levels=niveis, cmap="Blues")
    ax.tricontour(tri, z, levels=niveis, colors="white", linewidths=0.3)
    ax.tricontour(tri, z, levels=[0.], colors=PALETA[1], linewidths=1.6)
    ax.plot(*zip(*(TALUDE + [TALUDE[0]])), color=TINTA2, lw=0.9)
    for x, topo, nome in ((10., 10., "X–X"), (15., 7.5, "Y–Y")):
        ax.plot([x, x], [0, topo], color=TINTA, lw=0.9, ls="--")
        ax.text(x + 0.25, 0.4, nome, fontsize=7, color=TINTA)
    ax.set_aspect("equal")
    ax.set_xlim(-0.5, 30.5)
    ax.set_ylim(-0.5, 10.5)
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax.grid(False)
    cb = fig.colorbar(cs, ax=ax, shrink=0.9, pad=0.01)
    cb.set_label(r"$p_w$ (kPa)")
    ax.legend([Line2D([], [], color=PALETA[1], lw=1.6)], [r"nível d'água ($p_w = 0$)"], loc="upper right", fontsize=7)
    painel(ax, "a")
    for k, (s, x, topo) in enumerate((("XX", 10., 10.), ("YY", 15., 7.5))):
        ax = fig.add_subplot(g[1, k])
        ys, v = r["perfil_" + s]
        lf = liu["fig3"][s]
        ax.axvline(0, color=TINTA2, lw=0.6, ls="--")
        ax.plot(v, ys, color=PALETA[0], lw=1.5, label="este trabalho (tri3)")
        ax.plot(lf[:, 1], lf[:, 0], "o", ms=3.2, mfc="none", mec=LIU_COR, mew=0.9, label="Liu et al., Fig. 3")
        ax.set_xlim(-20, 52)
        ax.set_ylim(0, topo + 0.2)
        ax.set_xlabel(r"$p_w$ (kPa)")
        ax.set_ylabel("cota y (m)")
        ax.set_title("seção %s (x = %g m): nível %s m (Liu: %s m)" % (
            s[0] + "–" + s[1], x, ("%.2f" % r["nivel_" + s]).replace(".", ","),
            ("%.2f" % r["nivel_liu_" + s]).replace(".", ",")), fontsize=7.5)
        if k == 0:
            ax.legend(loc="upper right", fontsize=6.8)
        painel(ax, "bc"[k])
    salva(fig, "fig_pressoes.pdf")


def figura_hist(r, por):
    fig, axs = plt.subplots(2, 2, figsize=(7.0, 5.9), layout="constrained")
    for k, i in enumerate(CEN):
        e, x = r[i], por[i]["mc_fs"]
        mu, sd = e["mu"], e["sd"]
        zl = zeta(e["cv"])
        ml = math.log(mu) - 0.5 * zl * zl
        zL = zeta(e["liu_cv"])
        mL = math.log(e["liu_mu"]) - 0.5 * zL * zL
        xx = np.linspace(0.85, 2.3, 600)
        # (a, b) densidade
        ax = axs[0, k]
        ax.axvspan(0.8, 1., color=GRADE, alpha=0.5, lw=0)
        ax.hist(x, bins=45, density=True, color=COR[i], alpha=0.3, edgecolor="white", linewidth=0.3)
        ax.plot(xx, stats.lognorm.pdf(xx, zl, scale=math.exp(ml)), color=COR[i], lw=1.4)
        ax.plot(xx, stats.norm.pdf(xx, mu, sd), color=COR[i], lw=1.0, ls=":")
        ax.plot(xx, stats.lognorm.pdf(xx, zL, scale=math.exp(mL)), color=LIU_COR, lw=1.3, ls=(0, (5, 2)))
        ax.axvline(e["liu_min"], color=LIU_COR, lw=0.7, ls="-.")
        ax.set_xlim(0.85, 2.3)
        ax.set_xlabel("FS")
        ax.set_ylabel("densidade")
        ax.set_title("%s: %s realizações, sem falha" % (NOME[i], inteiro(len(x)).replace(r"\,", " ")), fontsize=8.5)
        painel(ax, "ab"[k])
        # (c, d) a cauda: P(FS < f) em escala log
        ax = axs[1, k]
        xs = np.sort(x)
        ff = np.linspace(0.9, 1.6, 400)
        emp = np.arange(1, len(xs) + 1) / len(xs)
        ax.step(xs, emp, where="post", color=COR[i], alpha=0.5, lw=1.0)
        ax.plot(ff, stats.lognorm.cdf(ff, zl, scale=math.exp(ml)), color=COR[i], lw=1.4)
        ax.plot(ff, stats.norm.cdf(ff, mu, sd), color=COR[i], lw=1.0, ls=":")
        ax.plot(ff, stats.lognorm.cdf(ff, zL, scale=math.exp(mL)), color=LIU_COR, lw=1.3, ls=(0, (5, 2)))
        lo, hi = ls_ic(e)
        ax.errorbar([1.], [e["ls"]["P_f"]], yerr=[[e["ls"]["P_f"] - lo], [hi - e["ls"]["P_f"]]], fmt="o",
                    color=PALETA[2], ms=4.5, capsize=2.5, lw=1.1, zorder=5)
        ax.plot([1.], [TAB4[i][0]], "s", ms=5, mfc="white", mec=LIU_COR, mew=1.1, zorder=5)
        ax.plot([1.], [e["wilson"][1]], "v", ms=5.5, color=COR[i], zorder=5)
        if i == "S3":
            ax.plot([1.], [TAB2_LHS], "D", ms=5, mfc="white", mec=PALETA[7], mew=1.1, zorder=5)
        ax.axvline(1., color=TINTA2, lw=0.6)
        ax.set_yscale("log")
        ax.set_ylim(1e-7, 1)
        ax.set_xlim(0.9, 1.6)
        ax.set_xlabel("valor f do FS")
        ax.set_ylabel(r"$P(FS < f)$")
        painel(ax, "cd"[k])
    fig.legend([Patch(color=TINTA2, alpha=0.3), Line2D([], [], color=TINTA2, lw=1.4),
                      Line2D([], [], color=TINTA2, lw=1.0, ls=":"), Line2D([], [], color=LIU_COR, ls=(0, (5, 2))),
                      Line2D([], [], color=LIU_COR, lw=0.7, ls="-.")],
                     ["nosso MC", "lognormal (nossa média e desvio)", "normal (nossa média e desvio)",
                      "lognormal de Liu (desvio da Fig. 11, média da Eq. 13)",
                      "menor FS das 2000 de Liu (Fig. 11)"], loc="outside lower center", ncol=3, fontsize=6.6)
    axs[1, 1].legend([Line2D([], [], color=PALETA[2], marker="o", ls="none"),
                      Line2D([], [], color=LIU_COR, marker="s", mfc="white", ls="none"),
                      Line2D([], [], color=TINTA2, marker="v", ls="none"),
                      Line2D([], [], color=PALETA[7], marker="D", mfc="white", ls="none")],
                     ["LS (IC 95 %)", "Liu, Tabela 4", "MC: limite de Wilson (0 falhas)", "Liu, Tabela 2 (LHS)"],
                     loc="lower right", fontsize=6.0)
    salva(fig, "fig_hist.pdf")


def figura_pf(r):
    fig, axs = plt.subplots(1, 2, figsize=(7.0, 3.2), sharey=True, constrained_layout=True)
    rot = ["FORM modal", "line sampling", "MC direto (0 falhas)", "lognormal do nosso MC", "normal do nosso MC",
           "lognormal reduzida (Bishop/SRM)", "Liu, Tabela 4", "Liu, Tabela 2 (LHS)", "Liu, Tabela 2 (CDF)"]
    ys = np.arange(len(rot))[::-1]
    for k, i in enumerate(CEN):
        ax, e = axs[k], r[i]
        ax.axvspan(*ls_ic(e), color=PALETA[2], alpha=0.12, lw=0)
        y = dict(zip(rot, ys))
        ax.plot(e["form"]["P_f"], y["FORM modal"], "D", color=COR[i], ms=5)
        lo, hi = ls_ic(e)
        ax.errorbar(e["ls"]["P_f"], y["line sampling"], xerr=[[e["ls"]["P_f"] - lo], [hi - e["ls"]["P_f"]]], fmt="o",
                    color=PALETA[2], ms=5, capsize=2.5)
        w = e["wilson"][1]
        ax.annotate("", xy=(1e-8, y["MC direto (0 falhas)"]), xytext=(w, y["MC direto (0 falhas)"]),
                    arrowprops={"arrowstyle": "->", "color": COR[i], "lw": 1.1})
        ax.plot(w, y["MC direto (0 falhas)"], "|", color=COR[i], ms=9, mew=1.6)
        ax.plot(e["pf_ln"], y["lognormal do nosso MC"], "o", color=COR[i], ms=5)
        ax.plot(e["pf_n"], y["normal do nosso MC"], "o", mfc="white", mec=COR[i], ms=5)
        ax.plot(e["pf_ln_red"], y["lognormal reduzida (Bishop/SRM)"], "o", color=COR[i], ms=5, alpha=0.55)
        ax.plot(TAB4[i][0], y["Liu, Tabela 4"], "s", mfc="white", mec=LIU_COR, ms=5.5, mew=1.1)
        if i == "S3":
            ax.plot(TAB2_LHS, y["Liu, Tabela 2 (LHS)"], "s", mfc="white", mec=PALETA[7], ms=5.5, mew=1.1)
            ax.plot(TAB2_CDF, y["Liu, Tabela 2 (CDF)"], "s", mfc="white", mec=PALETA[7], ms=5.5, mew=1.1)
        ax.set_xscale("log")
        ax.set_xlim(1e-8, 3e-2)
        ax.set_yticks(ys)
        ax.set_yticklabels(rot, fontsize=7)
        ax.set_xlabel(r"$P_f$")
        ax.set_title(NOME[i] + (r" ($c'$, $\varphi'$)" if i == "S2" else r" ($c'$, $\varphi'$, $k_s$)"), fontsize=8.5)
        ax.grid(axis="y", visible=False)
        painel(ax, "ab"[k])
    salva(fig, "fig_pf.pdf")


def figura_convergencia(r, por):
    fig, axs = plt.subplots(1, 2, figsize=(7.0, 3.0), constrained_layout=True)
    for k, i in enumerate(CEN):
        ax, e = axs[k], r[i]
        p = np.array([x[3] for x in por[i]["linhas_csv"]])
        n = np.arange(1, len(p) + 1)
        mu = np.cumsum(p) / n
        cv = np.array([np.std(p[:j], ddof=1) / math.sqrt(j) / mu[j - 1] if j > 1 else np.nan for j in n])
        s = np.sqrt(np.log1p(cv ** 2))
        ax.fill_between(n, mu * np.exp(-1.96 * s), mu * np.exp(1.96 * s), color=COR[i], alpha=0.2, lw=0)
        ax.plot(n, mu, color=COR[i], lw=1.3)
        ax.axhline(TAB4[i][0], color=LIU_COR, lw=1.0, ls=(0, (5, 2)))
        ax.axhline(e["pf_ln"], color=COR[i], lw=0.8, ls=":")
        ax.axhline(e["form"]["P_f"], color=TINTA2, lw=0.8, ls="-.")
        ax.set_yscale("log")
        ax.set_ylim(1e-7, 1e-3)
        ax.set_xlim(1, len(p))
        ax.set_xlabel("linhas")
        ax.set_ylabel(r"$P_f$ do LS")
        ax.set_title("%s: %d linhas, CoV final %s %%" % (NOME[i], len(p), ("%.0f" % (100 * e["ls"]["CoV"]))),
                     fontsize=8.5)
        painel(ax, "ab"[k])
    axs[0].legend([Line2D([], [], color=TINTA2), Patch(color=TINTA2, alpha=0.25),
                   Line2D([], [], color=LIU_COR, ls=(0, (5, 2))), Line2D([], [], color=TINTA2, ls=":"),
                   Line2D([], [], color=TINTA2, ls="-.")],
                  ["LS acumulado", "IC 95 % (log)", "Liu, Tabela 4", "lognormal do nosso MC", r"FORM, $\Phi(-\beta)$"],
                  loc="upper right", fontsize=6.3, ncol=2)
    salva(fig, "fig_convergencia.pdf")


def figura_raizes(r, por):
    fig, axs = plt.subplots(1, 2, figsize=(7.0, 2.9), constrained_layout=True)
    ax = axs[0]
    bins = np.arange(2.5, 8.01, 0.2)
    for i in CEN:
        cs = [x[1] for x in por[i]["linhas_csv"]]
        ax.hist(cs, bins=bins, histtype="step", color=COR[i], lw=1.3, label=NOME[i])
        ax.axvline(r[i]["form"]["beta"], color=COR[i], lw=1.0, ls="--")
        ax.axvline(r[i]["ls"]["beta_eq"], color=COR[i], lw=1.0, ls=":")
        ax.axvline(TAB4[i][1], color=COR[i], lw=0.8, ls="-.", alpha=0.7)
    ax.set_xlabel(r"raiz $c$ da linha (distância ao hiperplano, em $\sigma$)")
    ax.set_ylabel("linhas")
    h, l = ax.get_legend_handles_labels()
    h += [Line2D([], [], color=TINTA2, ls="--"), Line2D([], [], color=TINTA2, ls=":"),
          Line2D([], [], color=TINTA2, ls="-.")]
    l += [r"$\beta$ do FORM", r"$-\Phi^{-1}(P_f^{LS})$", r"$\beta$ de Liu (Tab. 4)"]
    ax.legend(h, l, loc="upper right", fontsize=6.3)
    painel(ax, "a")
    ax = axs[1]
    for i in CEN:
        p = np.sort([x[3] for x in por[i]["linhas_csv"]])[::-1]
        ax.plot(np.arange(1, len(p) + 1), np.cumsum(p) / p.sum(), color=COR[i], lw=1.3, label=NOME[i])
    ax.set_xscale("log")
    ax.set_xlabel("as k linhas de maior probabilidade")
    ax.set_ylabel(r"fração da $P_f$ do LS")
    ax.set_ylim(0, 1.02)
    ax.legend(loc="lower right", fontsize=6.5)
    painel(ax, "b")
    salva(fig, "fig_raizes.pdf")


def figura_k(r, por):
    par = r["par"]
    fig, axs = plt.subplots(1, 2, figsize=(7.0, 3.0), constrained_layout=True)
    ax = axs[0]
    ax.plot(par["f2"], par["f3"], ".", ms=1.6, color=PALETA[6], alpha=0.5)
    lim = [min(par["f2"].min(), par["f3"].min()) - 0.02, max(par["f2"].max(), par["f3"].max()) + 0.02]
    ax.plot(lim, lim, color=TINTA2, lw=0.8)
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_aspect("equal")
    ax.set_xlabel(r"FS no cenário 2 ($k_s$ fixo)")
    ax.set_ylabel(r"FS no cenário 3 (campo de $k_s$)")
    ax.set_title("as %s realizações pareadas" % inteiro(par["N"]).replace(r"\,", " "), fontsize=8.5)
    painel(ax, "a")
    ax = axs[1]
    v = por["S3"]["mc_nivel_xx"]
    w = por["S3"]["mc_nivel_yy"]
    bins = np.arange(min(v.min(), w.min()) - 0.05, max(v.max(), w.max()) + 0.15, 0.1)
    ax.hist(v, bins=bins, color=PALETA[0], alpha=0.45, label="X–X (x = 10 m)")
    ax.hist(w, bins=bins, color=PALETA[2], alpha=0.45, label="Y–Y (x = 15 m)")
    ax.axvline(float(np.median(por["S2"]["mc_nivel_xx"])), color=PALETA[0], lw=1.2, ls="--")
    ax.axvline(float(np.median(por["S2"]["mc_nivel_yy"])), color=PALETA[2], lw=1.2, ls="--")
    ax.set_xlabel("cota do nível d'água (m)")
    ax.set_ylabel("realizações")
    ax.set_title(r"cenário 3; tracejado: $k_s$ fixo (cenário 2)", fontsize=8.5)
    ax.legend(loc="upper left", fontsize=6.5)
    painel(ax, "b")
    salva(fig, "fig_k.pdf")


# ----------------------------------------------------------------------------------------------------------- numeros
def numeros(r, por):
    d = {}
    s2, s3 = r["S2"], r["S3"]
    d["fsDet"] = br(r["fsdet"], 4)
    d["fsDetTres"] = br(r["fsdet"], 3)
    d["difDet"] = br(100 * (r["fsdet"] / FS_BISHOP - 1), 1)
    d["difDetRel"] = br(100 * (1 - r["razao_det"]), 1)
    d["fsSemSuccao"] = br(r["pc"]["phib0"], 3)
    d["fsVert"] = br(r["pc"]["vert"], 4)
    d["fsFinal"] = br(r["pc"]["final"], 4)
    d["fsTquatro"] = br(r["pc"]["t4"], 4)
    d["fsAevZero"] = br(r["pc"]["aev0"], 4)
    d["fsXiZeroDois"] = br(por["S2"]["fs_xi0"], 3)
    d["fsXiZeroTres"] = br(por["S3"]["fs_xi0"], 3)
    for s, nm in (("XX", "XX"), ("YY", "YY")):
        d["nivel" + nm] = br(r["nivel_" + s], 2)
        d["nivelLiu" + nm] = br(r["nivel_liu_" + s], 2)
        d["succ" + nm] = br(-r["succ_" + s], 1)
        d["succLiu" + nm] = br(-r["succ_liu_" + s], 1)
        d["difSat" + nm] = br(r["dif_sat_" + s], 1)
        d["difMax" + nm] = br(r["dif_max_" + s], 1)
    for i, nm in (("S2", "Dois"), ("S3", "Tres")):
        e = r[i]
        d["beta" + nm] = br(e["form"]["beta"], 2)
        d["pfFORM" + nm] = cient(e["form"]["P_f"])
        d["alfaPhi" + nm] = br(e["form"]["alfa2_phi"], 2)
        d["alfaC" + nm] = br(e["form"]["alfa2_c"], 2)
        d["pfLS" + nm] = cient(e["ls"]["P_f"])
        d["cvLS" + nm] = br(100 * e["ls"]["CoV"], 0)
        d["linhas" + nm] = inteiro(e["ls"]["linhas"])
        d["anLS" + nm] = inteiro(e["an_ls"])
        d["anLinha" + nm] = br(e["an_linha"], 2)
        d["betaLS" + nm] = br(e["ls"]["beta_eq"], 2)
        lo, hi = ls_ic(e)
        d["icLSlo" + nm], d["icLShi" + nm] = cient(lo), cient(hi)
        d["nMC" + nm] = inteiro(e["N"])
        d["nfMC" + nm] = str(e["nf"])
        d["wilson" + nm] = cient(e["wilson"][1])
        d["mu" + nm] = br(e["mu"], 3)
        d["sd" + nm] = br(e["sd"], 3)
        d["cv" + nm] = br(100 * e["cv"], 1)
        d["minMC" + nm] = br(e["min"], 3)
        d["assim" + nm] = br(e["assim"], 2)
        d["betaLN" + nm] = br(e["b_ln"], 2)
        d["pfLN" + nm] = cient(e["pf_ln"])
        d["pfN" + nm] = cient(e["pf_n"])
        d["razaoLNLS" + nm] = br(e["razao_ln_ls"], 2)
        d["betaLNred" + nm] = br(e["b_ln_red"], 2)
        d["pfLNred" + nm] = cient(e["pf_ln_red"])
        d["pfNred" + nm] = cient(e["pf_n_red"])
        d["redReq" + nm] = br(100 * e["red_req"], 1)
        d["dbeta" + nm] = br(e["dbeta_dlnmu"], 1)
        d["liuMu" + nm] = br(e["liu_mu"], 3)
        d["liuSd" + nm] = br(e["liu_sd"], 3)
        d["liuCv" + nm] = br(100 * e["liu_cv"], 1)
        d["liuMin" + nm] = br(e["liu_min"], 2)
        d["liuPfN" + nm] = cient(e["liu_pf_n"])
        d["razaoTabQuatro" + nm] = br(TAB4[i][0] / e["ls"]["P_f"], 1)
        d["razaoTabQuatroLN" + nm] = br(TAB4[i][0] / e["pf_ln"], 1)
        d["razaoTabQuatroRed" + nm] = br(TAB4[i][0] / e["pf_ln_red"], 1)
        d["pZeroTabDois" + nm] = cient(e["p0_tab2"])
        d["pZeroTabQuatro" + nm] = br(e["p0_tab4"], 2)
        d["nTabDois" + nm] = br(e["n_tab2"], 0)
        d["razaoFORM" + nm] = br(e["ls"]["P_f"] / e["form"]["P_f"], 0)
        d["nCquatro" + nm] = str(e["n_c4"])
        d["fracCquatro" + nm] = br(100 * e["frac_c4"], 0)
        d["fracBase" + nm] = br(100 * e["frac_base"], 0)
        d["razaoNLS" + nm] = br(e["pf_n"] / e["ls"]["P_f"], 0)
        d["razaoTabDoisRed" + nm] = br(TAB2_LHS / e["pf_ln_red"], 0)
        d["fracMaior" + nm] = br(100 * e["frac_maior"], 0)
        d["fracDez" + nm] = br(100 * e["frac_10"], 0)
        d["cMin" + nm], d["cMax" + nm] = br(e["c_min"], 2), br(e["c_max"], 2)
        d["cMed" + nm] = br(e["c_med"], 2)
        d["cPcinco" + nm], d["cPnoventa" + nm] = br(e["c_p05"], 2), br(e["c_p95"], 2)
        d["tempo" + nm] = br(e["s_med"], 0)
        m = e["niv_xx"]
        d["nivXXmu" + nm], d["nivXXsd" + nm] = br(m[0], 2), br(m[1], 2)
        d["nivXXlo" + nm], d["nivXXhi" + nm] = br(m[2], 1), br(m[3], 1)
        m = e["niv_yy"]
        d["nivYYmu" + nm], d["nivYYsd" + nm] = br(m[0], 2), br(m[1], 2)
        d["nivYYlo" + nm], d["nivYYhi" + nm] = br(m[2], 1), br(m[3], 1)
    d["assimLN"] = br(3 * s2["cv"] + s2["cv"] ** 3, 2)
    d["pZeroLiu"] = cient((1. - TAB2_LHS) ** 2000)
    d["redLiuMu"] = br(100 * (1 - s2["liu_mu"] / s2["mu"]), 1)
    d["razaoDet"] = br(r["razao_det"], 3)
    d["razaoLS"] = br(r["razao_ls"], 1)
    d["razaoLScv"] = br(100 * r["razao_ls_cv"], 0)
    d["razaoTabQuatro"] = br(r["razao_tab4"], 1)
    d["razaoTabDois"] = br(TAB2_LHS / TAB4["S3"][0], 0)
    d["pfTabDoisCDF"] = cient(TAB2_CDF, 2)
    par = r["par"]
    d["nPar"] = inteiro(par["N"])
    d["parDmu"] = br(par["d_mu"], 4)
    d["parDsd"] = br(par["d_sd"], 3)
    d["parDlo"], d["parDhi"] = br(par["d_p"][0], 3), br(par["d_p"][1], 3)
    d["parCorr"] = br(par["corr"], 3)
    d["parSdDois"], d["parSdTres"] = br(par["sd2"], 3), br(par["sd3"], 3)
    d["parMuDois"], d["parMuTres"] = br(par["mu2"], 3), br(par["mu3"], 3)
    d["nTot"] = inteiro(sum(r[i]["form"]["analises"] + r[i]["an_ls"] + r[i]["N"] for i in CEN))
    linhas = ["% gerado por gera_dados.py; nao edite"]
    linhas += [r"\newcommand{\%s}{%s}" % (k, v) for k, v in d.items()]
    return "\n".join(linhas) + "\n", d


def main():
    os.makedirs(FIGS, exist_ok=True)
    os.makedirs(TAB, exist_ok=True)
    casos, por, agua = carrega_campanha()
    pc, liu = carrega_pc(), carrega_liu()
    r = contas(por, pc, liu, agua)
    blocos = {"det": tabela_det(r, por), "pf": tabela_pf(r), "aprox": tabela_aprox(r), "custo": tabela_custo(r, por)}
    for nome, texto in blocos.items():
        with open(os.path.join(TAB, nome + ".tex"), "w") as f:
            f.write("% gerado por gera_dados.py; nao edite\n" + texto)
    texto, d = numeros(r, por)
    with open(os.path.join(AQUI, "numeros.tex"), "w") as f:
        f.write(texto)
    figura_pressoes(r, liu)
    figura_hist(r, por)
    figura_pf(r)
    figura_convergencia(r, por)
    figura_raizes(r, por)
    figura_k(r, por)
    for i in CEN:
        e = r[i]
        print("%s: LS %.3e (CoV %.0f%%, %d linhas); MC 0/%d? nf=%d, FS %.3f +- %.3f; LN %.2e; LN reduzida %.2e; Tab.4 %.2e"
              % (i, e["ls"]["P_f"], 100 * e["ls"]["CoV"], e["ls"]["linhas"], e["N"], e["nf"], e["mu"], e["sd"],
                 e["pf_ln"], e["pf_ln_red"], TAB4[i][0]))
    print("FS det %.4f; nivel X-X %.2f (Liu %.2f), Y-Y %.2f (Liu %.2f); razao S3/S2 do LS %.2f (Liu %.2f)"
          % (r["fsdet"], r["nivel_XX"], r["nivel_liu_XX"], r["nivel_YY"], r["nivel_liu_YY"], r["razao_ls"],
             r["razao_tab4"]))


if __name__ == "__main__":
    main()
