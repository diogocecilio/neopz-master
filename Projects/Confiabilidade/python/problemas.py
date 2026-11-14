"""problemas.py: os resultados do documento doc/confiabilidade.tex, com a biblioteca confiab.py.

Problema 1 (o do codigo original): g = X1 X2 - X3, normais independentes N(40; 5), N(50; 2,5), N(1000; 200).
Problema 2 (proposto): talude infinito sob chuva, com os parametros do talude de Vargas Ceron et al. (2025): c' e phi'
lognormais correlacionados e a carga de pressao h no plano de ruptura normal, com a media mu_h dada pelo cenario de chuva.

Para cada um: probabilidade de falha de referencia (integracao numerica), FORM (HL-RF com e sem busca de Armijo), Monte
Carlo e amostragem por importancia no ponto de projeto; no problema 2, tambem o efeito dos erros do codigo original e a
curva de fragilidade P_f(mu_h). Grava doc/tabelas/*.tex, doc/figuras/*.pdf e dados/*.csv.
Uso: python3 problemas.py   (cerca de 1 min)
"""
import os
import sys

import numpy as np
from scipy import integrate, stats

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)
import confiab as cf  # noqa: E402

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.ticker import FuncFormatter  # noqa: E402

DOC = os.path.join(AQUI, "..", "doc")
FIG, TAB, DADOS = os.path.join(DOC, "figuras"), os.path.join(DOC, "tabelas"), os.path.join(AQUI, "..", "dados")
for d in (FIG, TAB, DADOS):
    os.makedirs(d, exist_ok=True)
SEMENTE = 20261114
rng = np.random.default_rng(SEMENTE)

PALETA = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#6250d6", "#e34948"]
TINTA, TINTA2, GRADE, NEUTRO = "#0b0b0b", "#52514e", "#dddcd8", "#8a8984"
plt.rcParams.update({
    "font.size": 8.5, "axes.labelsize": 8.5, "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "legend.fontsize": 7,
    "axes.titlesize": 8.5, "axes.edgecolor": TINTA2, "axes.labelcolor": TINTA, "xtick.color": TINTA2,
    "ytick.color": TINTA2, "axes.grid": True, "grid.color": GRADE, "grid.linewidth": 0.5, "lines.linewidth": 1.4,
    "legend.frameon": False, "savefig.dpi": 200, "font.family": "DejaVu Sans", "mathtext.fontset": "dejavusans",
    "axes.unicode_minus": True})
LARGURA = 6.3
virg = FuncFormatter(lambda v, _: ("%g" % v).replace(".", ",").replace("-", "−"))


# ---------------------------------------------------------------- numeros para o LaTeX (virgula decimal)
def tf(v, casas=3):
    s = "%.*f" % (casas, v)
    return "$" + s.replace("-", "-").replace(".", "{,}") + "$"


def ts(v, sig=3):
    """notacao cientifica: 1,18 x 10^-3"""
    if v == 0:
        return "$0$"
    if not np.isfinite(v):
        return "---"
    m, e = ("%.*e" % (sig - 1, v)).split("e")
    e = int(e)
    if -2 <= e <= 0:
        return tf(v, sig - 1 - e)
    return "$%s\\times 10^{%d}$" % (m.replace(".", "{,}"), e)


def tpc(v, casas=1):
    return tf(100 * v, casas)[:-1] + "\\,\\%$"


def grava(nome, linhas):
    with open(os.path.join(TAB, nome), "w") as f:
        f.write("\n".join(linhas) + "\n")


macros = []


def macro(nome, valor):
    macros.append("\\newcommand{\\%s}{%s}" % (nome, valor))


# ================================================================ problema 1
m1 = cf.Modelo([cf.Normal(40, 5), cf.Normal(50, 2.5), cf.Normal(1000, 200)])


def g1(x):
    return x[..., 0] * x[..., 1] - x[..., 2]


def G1(u):
    return g1(m1.x_de_u(u))


def pf_exato_1():
    def f(u2, u1):
        return stats.norm.pdf(u1) * stats.norm.pdf(u2) * stats.norm.cdf((1000 - (40 + 5 * u1) * (50 + 2.5 * u2)) / 200)
    return integrate.dblquad(f, -9, 9, -9, 9, epsabs=0, epsrel=1e-10)[0]


pe1 = pf_exato_1()
f1 = cf.form(G1, 3)
f1s = cf.form(G1, 3, armijo=False)
x1s = m1.x_de_u(f1["u"])
print("P1: exato %.6e (beta %.5f); FORM beta %.5f pf %.5e (%d it, %d aval.)" % (pe1, cf.beta_de_pf(pe1), f1["beta"],
                                                                               f1["pf"], f1["iteracoes"],
                                                                               f1["avaliacoes"]))
mc1 = {N: cf.monte_carlo(G1, 3, N, rng) for N in (10 ** 3, 10 ** 4, 10 ** 5, 10 ** 6)}
is1 = {N: cf.importancia(G1, f1["u"], N, rng) for N in (10 ** 3, 10 ** 4)}
# repeticoes com N = 1000: espalhamento das duas estimativas
R = 2000
rep_mc = np.array([cf.monte_carlo(G1, 3, 1000, rng)["pf"] for _ in range(R)])
rep_is = np.array([cf.importancia(G1, f1["u"], 1000, rng)["pf"] for _ in range(R)])
zero_mc = (rep_mc == 0).mean()
print("P1: MC N=1000: P(pf=0) = %.3f (teorico %.3f); IS N=1000: CoV empirico %.4f" % (
    zero_mc, (1 - pe1) ** 1000, rep_is.std() / rep_is.mean()))

lin = ["\\begin{tabular}{@{}l r r r r r@{}}", "\\toprule",
       "Método & $N$ (avaliações de $g$) & $P_f$ & $\\beta$ & CoV estimado & erro de $P_f$ \\\\", "\\midrule",
       "Referência (integração numérica) & --- & %s & %s & --- & --- \\\\" % (ts(pe1, 4), tf(cf.beta_de_pf(pe1), 4)),
       "FORM, HL-RF (diferenças centrais) & %d & %s & %s & --- & %s \\\\" % (
           f1["avaliacoes"], ts(f1["pf"], 4), tf(f1["beta"], 4), tpc(f1["pf"] / pe1 - 1))]
for N, r in mc1.items():
    lin.append("Monte Carlo & $10^{%d}$ & %s & %s & %s & %s \\\\" % (
        round(np.log10(N)), ts(r["pf"], 3), tf(cf.beta_de_pf(r["pf"]), 3), tpc(r["cov"]), tpc(r["pf"] / pe1 - 1)))
for N, r in is1.items():
    lin.append("Amostragem por importância & $10^{%d}$ & %s & %s & %s & %s \\\\" % (
        round(np.log10(N)), ts(r["pf"], 4), tf(cf.beta_de_pf(r["pf"]), 4), tpc(r["cov"]), tpc(r["pf"] / pe1 - 1)))
lin += ["\\bottomrule", "\\end{tabular}"]
grava("p1_resultados.tex", lin)
macro("pUmExato", ts(pe1, 4))
macro("pUmBetaExato", tf(cf.beta_de_pf(pe1), 4))
macro("pUmBeta", tf(f1["beta"], 5))
macro("pUmPfForm", ts(f1["pf"], 4))
macro("pUmErroForm", tpc(f1["pf"] / pe1 - 1))
macro("pUmIter", str(f1["iteracoes"]))
macro("pUmAval", str(f1["avaliacoes"]))
macro("pUmUestrela", "(%s)" % "; ".join(tf(v, 4).strip("$") for v in f1["u"]))
macro("pUmXestrela", "(%s)" % "; ".join(tf(v, 2).strip("$") for v in x1s))
macro("pUmAlfaDois", "(%s)" % "; ".join(tf(v, 3).strip("$") for v in f1["alfa"] ** 2))
macro("pUmZeroMC", tpc(zero_mc, 0))
macro("pUmZeroMCteo", tpc((1 - pe1) ** 1000, 0))
macro("pUmCovMCmil", tpc(np.sqrt((1 - pe1) / (1000 * pe1)), 0))
macro("pUmCovISmil", tpc(rep_is.std() / rep_is.mean(), 1))
macro("pUmNmcDez", ts((1 - pe1) / (pe1 * 0.01), 2))
macro("pUmRep", str(R))
q5, q95 = np.percentile(rep_mc, [5, 95])
macro("pUmMCcinco", ts(q5, 2))
macro("pUmMCnoventaecinco", ts(q95, 2))
q5, q95 = np.percentile(rep_is, [5, 95])
macro("pUmIScinco", ts(q5, 3))
macro("pUmISnoventaecinco", ts(q95, 3))

# figura: estimativa acumulada (MC e IS) e CoV contra N
fig, axs = plt.subplots(1, 2, figsize=(LARGURA, 2.5))
a = axs[0]
qmc = cf.monte_carlo(G1, 3, 10 ** 5, rng)["q"]
qis = cf.importancia(G1, f1["u"], 10 ** 5, rng)["q"]
for q, cor, nome in ((qmc, PALETA[1], "Monte Carlo"), (qis, PALETA[0], "importância (ponto de projeto)")):
    mm, ss = cf.media_acumulada(q)
    k = np.unique(np.round(np.logspace(0, np.log10(len(q)), 1500)).astype(int))   # pontos espacados em log
    a.fill_between(k, np.maximum(mm[k - 1] - 2 * ss[k - 1], 1e-6), mm[k - 1] + 2 * ss[k - 1], color=cor, alpha=0.18,
                   lw=0)
    a.plot(k, mm[k - 1], color=cor, lw=1.1, label=nome)
a.axhline(pe1, color=TINTA, lw=1.0, label="referência")
a.axhline(f1["pf"], color=NEUTRO, lw=1.0, ls="--", label="FORM")
a.set_xscale("log")
a.set_ylim(0, 3e-3)
a.set_xlim(30, 1e5)
a.set_xlabel("amostras N")
a.set_ylabel("estimativa acumulada de $P_f$")
a.yaxis.set_major_formatter(FuncFormatter(lambda v, _: ("%g" % (v * 1e3)).replace(".", ",")))
a.set_title("(a) $P_f$ ($\\times 10^{-3}$) e faixa de ±2 desvios", loc="left")
a.legend(loc="upper right")
a = axs[1]
Ns = np.array([100, 300, 1000, 3000, 10000, 30000, 100000])
a.loglog(Ns, np.sqrt((1 - pe1) / (Ns * pe1)), color=PALETA[1], label="Monte Carlo")
cis = [np.mean([cf.importancia(G1, f1["u"], int(N), rng)["cov"] for _ in range(20)]) for N in Ns]
a.loglog(Ns, cis, "o-", color=PALETA[0], ms=3, label="importância")
a.axhline(0.1, color=NEUTRO, lw=0.8, ls=":")
a.set_xlabel("amostras N")
a.set_ylabel("CoV de $P_f$")
a.yaxis.set_major_formatter(virg)
a.set_title("(b) coeficiente de variação da estimativa", loc="left")
a.legend()
fig.tight_layout()
fig.savefig(os.path.join(FIG, "p1_convergencia.pdf"), bbox_inches="tight")
plt.close(fig)

# ================================================================ problema 2: talude infinito sob chuva
GAMA, GAMA_W, BETA_S, ZF, SIG_H = 20., 9.81, np.radians(45.), 2.0, 0.5
SIG_N = GAMA * ZF * np.cos(BETA_S) ** 2           # tensao normal total no plano (kPa)
TAU = GAMA * ZF * np.sin(BETA_S) * np.cos(BETA_S)  # tensao cisalhante (kPa)
C_MED, C_COV, PHI_MED, PHI_COV, RHO = 10., 0.3, 30., 0.1, -0.5


def fs2(x):
    c, phi, h = x[..., 0], x[..., 1], x[..., 2]
    return (c + (SIG_N - GAMA_W * h) * np.tan(np.radians(phi))) / TAU   # Bishop com chi = 1 (como no SlopeRandomField)


def g2(x):
    return fs2(x) - 1.


def modelo2(muh, rho=RHO):
    return cf.Modelo([cf.LogNormal(C_MED, C_COV), cf.LogNormal(PHI_MED, PHI_COV), cf.Normal(muh, SIG_H)],
                     [[1, rho, 0], [rho, 1, 0], [0, 0, 1]])


def pf_exato_2(muh, rho=RHO):
    """P(c' < c_crit(phi', h)): c' | phi' e lognormal; integra em phi' (z2) e h"""
    a, b = cf.LogNormal(C_MED, C_COV), cf.LogNormal(PHI_MED, PHI_COV)
    rz = cf.rho_nataf(a, b, rho)

    def f(h, z2):
        cc = TAU - (SIG_N - GAMA_W * h) * np.tan(np.radians(b.x_de_z(z2)))
        if cc <= 0:
            return 0.
        return (stats.norm.pdf(z2) * stats.norm.pdf(h, muh, SIG_H) *
                stats.norm.cdf((np.log(cc) - a.lam - a.zeta * rz * z2) / (a.zeta * np.sqrt(1 - rz * rz))))
    return integrate.dblquad(f, -8, 8, muh - 8 * SIG_H, muh + 8 * SIG_H, epsabs=0, epsrel=1e-9)[0]


class ModeloComErro(cf.Modelo):
    """o erro do FromXtoY original: y = L z no lugar de y = L^-1 z, isto e, z = L^-1 y na volta"""

    def __init__(self, dists, RX):
        super().__init__(dists, RX)
        self.L = np.linalg.inv(self.L)


class ModeloSemNataf(cf.Modelo):
    """a correlacao das variaveis fisicas usada direto no espaco normal (sem a correcao de Nataf)"""

    def __init__(self, dists, RX):
        super().__init__(dists, RX)
        self.RZ = np.asarray(RX, float)
        self.L = np.linalg.cholesky(self.RZ)


def dists2(muh):
    return [cf.LogNormal(C_MED, C_COV), cf.LogNormal(PHI_MED, PHI_COV), cf.Normal(muh, SIG_H)]


RX2 = [[1, RHO, 0], [RHO, 1, 0], [0, 0, 1]]
cen = {"A": -2.0, "B": -0.5}
res2 = {}
for k, muh in cen.items():
    m = modelo2(muh)

    def G(u, m=m):
        return g2(m.x_de_u(u))
    r = dict(muh=muh, pe=pf_exato_2(muh), form=cf.form(G, 3), form_sem=cf.form(G, 3, armijo=False))
    r["x"] = m.x_de_u(r["form"]["u"])
    r["mc"] = {N: cf.monte_carlo(G, 3, N, rng) for N in (10 ** 3, 10 ** 6)}
    r["is"] = {N: cf.importancia(G, r["form"]["u"], N, rng) for N in (10 ** 3, 10 ** 4)}
    r["fs_med"] = fs2(np.array([C_MED, PHI_MED, muh]))
    # os erros do codigo original
    for nome, M in (("erro_L", ModeloComErro(dists2(muh), RX2)), ("sem_nataf", ModeloSemNataf(dists2(muh), RX2)),
                    ("indep", modelo2(muh, 0.))):
        r[nome] = cf.form(lambda u, M=M: g2(M.x_de_u(u)), 3)
    r["pe_indep"] = pf_exato_2(muh, 0.)
    res2[k] = r
    print("P2 %s (mu_h=%g): FS med %.3f exato %.4e FORM beta %.4f pf %.4e it %d; IS1e4 %.4e CoV %.3f; MC1e6 %.4e" % (
        k, muh, r["fs_med"], r["pe"], r["form"]["beta"], r["form"]["pf"], r["form"]["iteracoes"],
        r["is"][10 ** 4]["pf"], r["is"][10 ** 4]["cov"], r["mc"][10 ** 6]["pf"]))
    print("    erro L: beta %.4f; sem Nataf %.4f; independentes %.4f (exato %.4e)" % (
        r["erro_L"]["beta"], r["sem_nataf"]["beta"], r["indep"]["beta"], r["pe_indep"]))

Mbug = ModeloComErro(dists2(-1.), RX2)
cov_bug = np.linalg.inv(cf.Modelo(dists2(-1.), RX2).L) @ np.linalg.inv(cf.Modelo(dists2(-1.), RX2).L).T
macro("pDoisRhoZ", tf(cf.rho_nataf(cf.LogNormal(C_MED, C_COV), cf.LogNormal(PHI_MED, PHI_COV), RHO), 4))
macro("pDoisRhoErro", tf(cov_bug[0, 1] / np.sqrt(cov_bug[0, 0] * cov_bug[1, 1]), 3))
macro("pDoisDpErro", tf(np.sqrt(cov_bug[1, 1]), 3))
macro("pDoisSigN", tf(SIG_N, 1))
macro("pDoisTau", tf(TAU, 1))

lin = ["\\begin{tabular}{@{}l r r r r r r@{}}", "\\toprule",
       " & \\multicolumn{3}{c}{Cenário A: $\\mu_h = %s$\\,m} & \\multicolumn{3}{c@{}}{Cenário B: $\\mu_h = %s$\\,m} \\\\" % (
           tf(cen["A"], 1).strip("$"), tf(cen["B"], 1).strip("$")),
       "\\cmidrule(lr){2-4}\\cmidrule(l){5-7}",
       "Método & $N$ & $P_f$ & CoV ou erro & $N$ & $P_f$ & CoV ou erro \\\\", "\\midrule"]
A, B = res2["A"], res2["B"]
lin.append("Referência (integração) & --- & %s & --- & --- & %s & --- \\\\" % (ts(A["pe"], 4), ts(B["pe"], 4)))
lin.append("FORM (HL-RF com Armijo) & %d & %s & %s & %d & %s & %s \\\\" % (
    A["form"]["avaliacoes"], ts(A["form"]["pf"], 4), tpc(A["form"]["pf"] / A["pe"] - 1),
    B["form"]["avaliacoes"], ts(B["form"]["pf"], 4), tpc(B["form"]["pf"] / B["pe"] - 1)))
for N in (10 ** 3, 10 ** 6):
    lin.append("Monte Carlo & $10^{%d}$ & %s & %s & $10^{%d}$ & %s & %s \\\\" % (
        round(np.log10(N)), ts(A["mc"][N]["pf"], 3), tpc(A["mc"][N]["cov"]) if np.isfinite(A["mc"][N]["cov"]) else "---",
        round(np.log10(N)), ts(B["mc"][N]["pf"], 3), tpc(B["mc"][N]["cov"])))
for N in (10 ** 3, 10 ** 4):
    lin.append("Importância (ponto de projeto) & $10^{%d}$ & %s & %s & $10^{%d}$ & %s & %s \\\\" % (
        round(np.log10(N)), ts(A["is"][N]["pf"], 4), tpc(A["is"][N]["cov"]),
        round(np.log10(N)), ts(B["is"][N]["pf"], 4), tpc(B["is"][N]["cov"])))
lin += ["\\midrule", "$\\beta$ do FORM & & %s & & & %s & \\\\" % (tf(A["form"]["beta"], 4), tf(B["form"]["beta"], 4)),
        "FS com as médias & & %s & & & %s & \\\\" % (tf(A["fs_med"], 3), tf(B["fs_med"], 3)),
        "$\\alpha^2$ de $c'$, $\\phi'$, $h$ & & \\multicolumn{2}{l}{$%s$} & & \\multicolumn{2}{l@{}}{$%s$} \\\\" % (
            ";\\ ".join(tf(v, 2).strip("$") for v in A["form"]["alfa"] ** 2),
            ";\\ ".join(tf(v, 2).strip("$") for v in B["form"]["alfa"] ** 2)),
        "Ponto de projeto $x^*$ & & \\multicolumn{2}{l}{$%s$} & & \\multicolumn{2}{l@{}}{$%s$} \\\\" % (
            "%s\\,\\mathrm{kPa};\\ %s^\\circ;\\ %s\\,\\mathrm{m}" % tuple(tf(v, c).strip("$") for v, c in zip(A["x"], (2, 1, 2))),
            "%s\\,\\mathrm{kPa};\\ %s^\\circ;\\ %s\\,\\mathrm{m}" % tuple(tf(v, c).strip("$") for v, c in zip(B["x"], (2, 1, 2)))),
        "\\bottomrule", "\\end{tabular}"]
grava("p2_cenarios.tex", lin)

lin = ["\\begin{tabular}{@{}l r r r r@{}}", "\\toprule",
       " & \\multicolumn{2}{c}{Cenário A} & \\multicolumn{2}{c@{}}{Cenário B} \\\\",
       "\\cmidrule(lr){2-3}\\cmidrule(l){4-5}",
       "Transformação usada no FORM & $\\beta$ & $P_f$ do FORM & $\\beta$ & $P_f$ do FORM \\\\", "\\midrule"]
for nome, chave in (("correta (Nataf, $u = L^{-1} z$)", "form"), ("com o erro do \\texttt{FromXtoY} ($L$ no lugar de $L^{-1}$)", "erro_L"),
                    ("sem a correção de Nataf ($R_Z = R_X$)", "sem_nataf")):
    lin.append("%s & %s & %s & %s & %s \\\\" % (nome, tf(A[chave]["beta"], 4), ts(A[chave]["pf"], 3),
                                                tf(B[chave]["beta"], 4), ts(B[chave]["pf"], 3)))
lin += ["\\midrule",
        "$c'$ e $\\phi'$ independentes ($\\rho = 0$), FORM & %s & %s & %s & %s \\\\" % (
            tf(A["indep"]["beta"], 4), ts(A["indep"]["pf"], 3), tf(B["indep"]["beta"], 4), ts(B["indep"]["pf"], 3)),
        "$c'$ e $\\phi'$ independentes, referência & %s & %s & %s & %s \\\\" % (
            tf(cf.beta_de_pf(A["pe_indep"]), 4), ts(A["pe_indep"], 3), tf(cf.beta_de_pf(B["pe_indep"]), 4),
            ts(B["pe_indep"], 3)),
        "\\bottomrule", "\\end{tabular}"]
grava("p2_erros.tex", lin)
for k, r in res2.items():
    s = "A" if k == "A" else "B"
    macro("pDois%sExato" % s, ts(r["pe"], 3))
    macro("pDois%sBeta" % s, tf(r["form"]["beta"], 3))
    macro("pDois%sErroForm" % s, tpc(r["form"]["pf"] / r["pe"] - 1))
    macro("pDois%sCovIS" % s, tpc(r["is"][10 ** 3]["cov"]))
    macro("pDois%sFS" % s, tf(r["fs_med"], 3))
    macro("pDois%sErroL" % s, tf(r["erro_L"]["pf"] / r["form"]["pf"], 1))
    macro("pDois%sMCmil" % s, ts(r["mc"][10 ** 3]["pf"], 3))
    macro("pDois%sNmcDez" % s, ts((1 - r["pe"]) / (r["pe"] * 0.01), 2))
    macro("pDois%sNisDez" % s, ts(1000 * r["is"][10 ** 3]["cov"] ** 2 / 0.01, 2))

# ---- HL-RF puro contra HL-RF com Armijo, e a curva de fragilidade
mus = np.round(np.arange(-4.0, 0.51, 0.25), 2)
fr = []
for muh in mus:
    m = modelo2(muh)

    def G(u, m=m):
        return g2(m.x_de_u(u))
    fa, fp = cf.form(G, 3), cf.form(G, 3, armijo=False)
    isr = cf.importancia(G, fa["u"], 1000, rng)
    mcr = cf.monte_carlo(G, 3, 1000, rng)
    fr.append(dict(muh=muh, pe=pf_exato_2(muh), pe0=pf_exato_2(muh, 0.), beta=fa["beta"], pf_form=fa["pf"],
                   it=fa["iteracoes"], conv=fa["convergiu"], beta_puro=fp["beta"], it_puro=fp["iteracoes"],
                   conv_puro=fp["convergiu"], pf_is=isr["pf"], cov_is=isr["cov"], pf_mc=mcr["pf"],
                   a2c=fa["alfa"][0] ** 2, a2p=fa["alfa"][1] ** 2, a2h=fa["alfa"][2] ** 2,
                   fs=fs2(np.array([C_MED, PHI_MED, muh]))))
cols = list(fr[0].keys())
with open(os.path.join(DADOS, "p2_fragilidade.csv"), "w") as f:
    f.write(",".join(cols) + "\n")
    for r in fr:
        f.write(",".join("%.10g" % float(r[c]) for c in cols) + "\n")

lin = ["\\begin{tabular}{@{}r r r r r r r@{}}", "\\toprule",
       "$\\mu_h$ (m) & FS médio & $P_f$ (referência) & $\\beta$ (HL-RF com Armijo) & iterações & $\\beta$ (HL-RF puro) & iterações \\\\",
       "\\midrule"]
for r in fr:
    if r["muh"] in (-4., -3.5, -3., -2.5, -2., -1., 0., 0.5):
        lin.append("%s & %s & %s & %s & %d & %s & %s \\\\" % (
            tf(r["muh"], 2), tf(r["fs"], 3), ts(r["pe"], 3), tf(r["beta"], 4), r["it"], tf(r["beta_puro"], 4),
            ("%d" % r["it_puro"]) if r["conv_puro"] else "não converge (%d)" % r["it_puro"]))
lin += ["\\bottomrule", "\\end{tabular}"]
grava("p2_hlrf.tex", lin)
nao = [r["muh"] for r in fr if not r["conv_puro"]]
macro("pDoisNaoConvMax", tf(max(nao), 2) if nao else "---")

# figura: fragilidade e fatores de importancia
fig, axs = plt.subplots(1, 2, figsize=(LARGURA, 2.6))
a = axs[0]
mu = np.array([r["muh"] for r in fr])
a.semilogy(mu, [r["pe"] for r in fr], color=TINTA, lw=1.6, label="referência, $\\rho = -0{,}5$")
a.semilogy(mu, [r["pe0"] for r in fr], color=NEUTRO, lw=1.2, ls="--", label="referência, $\\rho = 0$")
a.semilogy(mu, [r["pf_form"] for r in fr], color=PALETA[3], lw=1.0, ls="-.", label="FORM")
pis = np.array([r["pf_is"] for r in fr])
cis = np.array([r["cov_is"] for r in fr])
a.errorbar(mu, pis, yerr=2 * cis * pis, fmt="o", ms=3, color=PALETA[0], elinewidth=0.8, capsize=1.5,
           label="importância, N = 1000 (±2 desvios)")
pmc = np.array([r["pf_mc"] for r in fr])
ok = pmc > 0
a.semilogy(mu[ok], pmc[ok], "s", ms=3, color=PALETA[1], label="Monte Carlo, N = 1000")
a.set_xlabel("média da carga de pressão no plano de ruptura, $\\mu_h$ (m)")
a.set_ylabel("$P_f$")
a.set_ylim(1e-13, 1)
a.set_title("(a) curva de fragilidade", loc="left")
a.legend(loc="upper left", fontsize=6.3)
a.xaxis.set_major_formatter(virg)
a = axs[1]
a.stackplot(mu, [r["a2c"] for r in fr], [r["a2p"] for r in fr], [r["a2h"] for r in fr],
            colors=[PALETA[0], PALETA[2], PALETA[1]], labels=["$c'$", "$\\phi'$", "$h$ (chuva)"], alpha=0.85)
a.set_xlim(mu[0], mu[-1])
a.set_ylim(0, 1)
a.set_xlabel("$\\mu_h$ (m)")
a.set_ylabel("$\\alpha_i^2$")
a.set_title("(b) fatores de importância do FORM", loc="left")
a.legend(loc="center left", fontsize=7, frameon=True, facecolor="white", edgecolor=GRADE)
a.xaxis.set_major_formatter(virg)
a.yaxis.set_major_formatter(virg)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "p2_fragilidade.pdf"), bbox_inches="tight")
plt.close(fig)

# figura: superficie de estado limite no plano (u_c, u_h) com u_phi no ponto de projeto, cenario A, e amostras
r = res2["A"]
m = modelo2(r["muh"])
us = r["form"]["u"]
ua, ub = np.meshgrid(np.linspace(-5, 5, 301), np.linspace(-3, 6, 301))
U = np.stack([ua, np.full_like(ua, us[1]), ub], axis=-1)
Gv = g2(m.x_de_u(U))
fig, ax = plt.subplots(figsize=(LARGURA * 0.75, 2.7))
ax.contourf(ua, ub, Gv, levels=[-50, 0], colors=["#f6d3cf"])
ax.contour(ua, ub, Gv, levels=[0], colors=[PALETA[6]], linewidths=1.2)
amc = rng.standard_normal((400, 3))
ais = us + rng.standard_normal((400, 3))
ax.plot(amc[:, 0], amc[:, 2], ".", ms=2, color=PALETA[1], label="amostras de Monte Carlo")
ax.plot(ais[:, 0], ais[:, 2], ".", ms=2, color=PALETA[0], label="amostras de importância")
ax.plot([0, us[0]], [0, us[2]], color=TINTA, lw=1)
ax.plot(us[0], us[2], "k*", ms=8, label="ponto de projeto $u^*$")
for rr in (1, 2, 3):
    t = np.linspace(0, 2 * np.pi, 200)
    ax.plot(rr * np.cos(t), rr * np.sin(t), color=NEUTRO, lw=0.5, ls=":")
ax.set_aspect("equal")
ax.set_xlim(-4, 4)
ax.set_ylim(-3, 6)
ax.set_xlabel("$u_1$ (coesão)")
ax.set_ylabel("$u_3$ (carga de pressão)")
ax.text(-3.8, 5.3, "falha, $G(u) \\leq 0$", fontsize=7, color=PALETA[6], bbox=dict(fc="white", ec="none", pad=1))
from matplotlib.lines import Line2D  # noqa: E402
ax.legend(handles=[Line2D([], [], ls="", marker=".", ms=6, color=PALETA[1], label="amostras de Monte Carlo"),
                   Line2D([], [], ls="", marker=".", ms=6, color=PALETA[0], label="amostras de importância"),
                   Line2D([], [], ls="", marker="*", ms=8, color="k", label="ponto de projeto $u^*$")],
          loc="center left", bbox_to_anchor=(1.02, 0.5), fontsize=6.5)
ax.xaxis.set_major_formatter(virg)
ax.yaxis.set_major_formatter(virg)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "p2_espaco_u.pdf"), bbox_inches="tight")
plt.close(fig)

with open(os.path.join(TAB, "valores.tex"), "w") as f:
    f.write("% gerado por python/problemas.py\n" + "\n".join(macros) + "\n")
print("ok:", len(macros), "macros")
