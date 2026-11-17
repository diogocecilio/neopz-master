#!/usr/bin/env python3
"""relatorio_etapas.py: todos os resultados das etapas 0 a 3 do plano de trabalho e a referencia da chuva da etapa 4
(numeros, figuras e o texto LaTeX), lidos das campanhas baixadas do fermi (~/projects/runs/fermi_etapa1 e fermi_etapa3)
e das rodadas do PC (~/projects/runs/cho2010: grossa, nystrom, gim_gradiente e diagnostico; o teste linear em
~/projects/runs/etapa2_linear). O relatorio de situacao (../situacao/gera_situacao.py) roda este arquivo e usa os mesmos
numeros no resumo. Para atualizar os tres PDFs (a formulacao, este e a situacao):
    Projects/PlanoDeTrabalho/situacao/atualiza.sh"""
import csv
import glob
import json
import math
import os
import re
import runpy
import subprocess

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

AQUI = os.path.dirname(os.path.abspath(__file__))  # Projects/PlanoDeTrabalho/ETAPAS1a3
PLANO = os.path.dirname(AQUI)
NEOPZ = os.path.dirname(os.path.dirname(PLANO))
RUNS = os.path.expanduser("~/projects/runs")
E1 = f"{RUNS}/fermi_etapa1/etapa1"
E3 = f"{RUNS}/fermi_etapa3/etapa3"
LIN = f"{RUNS}/etapa2_linear"
CHO2010 = os.path.expanduser("~/projects/neopz-master-build-release/Projects/Cho2010/Cho2010")
os.makedirs(f"{AQUI}/figs", exist_ok=True)

def le(pasta, t=None):
    """{real: linha} com fator (situacao ok, FS > 0); com t, so o instante t"""
    out = {}
    for a in glob.glob(pasta + "/mc_*_*.csv"):
        if not re.search(r"/mc_-?\d+_-?\d+\.csv$", a):
            continue
        for r in csv.DictReader(open(a, newline="")):
            n = int(r["real"])
            if n < 0 or r["situacao"] != "ok" or not float(r["FS_min"]) > 0:
                continue
            if t is not None and abs(float(r.get("t_h") or 0) - t) > 1e-9:
                continue
            out[n] = r
    return out


def phi(z):
    return 0.5 * math.erfc(-z / math.sqrt(2))


def clopper(k, n, a=0.05):
    """IC exato de uma proporcao (bissecao na binomial)"""
    def cdf(p, kk):  # em logaritmos: math.comb(4000, i) nao cabe num float
        if p <= 0.:
            return 1.
        if p >= 1.:
            return 0. if kk < n else 1.
        lg = lambda i: (math.lgamma(n + 1) - math.lgamma(i + 1) - math.lgamma(n - i + 1) + i * math.log(p)
                        + (n - i) * math.log1p(-p))
        return sum(math.exp(lg(i)) for i in range(kk + 1))
    lo = 0.
    if k > 0:
        a0, b0 = 0., 1.
        for _ in range(60):
            m = (a0 + b0) / 2
            (a0, b0) = (m, b0) if 1 - cdf(m, k - 1) < a / 2 else (a0, m)
        lo = a0
    hi = 1.
    if k < n:
        a0, b0 = 0., 1.
        for _ in range(60):
            m = (a0 + b0) / 2
            (a0, b0) = (m, b0) if cdf(m, k) > a / 2 else (a0, m)
        hi = b0
    return lo, hi


def direto(pasta, t=None):
    rs = le(pasta, t)
    n = len(rs)
    k = sum(float(r["FS_min"]) < 1 for r in rs.values())
    p = k / n
    lo, hi = clopper(k, n)
    return dict(n=n, k=k, p=p, cov=math.sqrt((1 - p) / (n * p)) if k else float("inf"), lo=lo, hi=hi,
                s=sum(float(r["s_fator"]) for r in rs.values()) / n)


def ai(pasta, t=None):
    rs = le(pasta, t)
    n = len(rs)
    wi = [math.exp(float(r["log_peso"])) * (float(r["FS_min"]) < 1) for r in rs.values()]
    p = sum(wi) / n
    s2 = sum(x * x for x in wi) / n
    sd = math.sqrt(max(s2 - p * p, 0) / (n - 1))
    f = [x for x in wi if x > 0]
    return dict(n=n, k=len(f), p=p, cov=sd / p, lo=p - 1.96 * sd, hi=p + 1.96 * sd, sd=sd,
                essF=sum(f) ** 2 / sum(x * x for x in f), maxp=max(f) / sum(f),
                mw=sum(math.exp(float(r["log_peso"])) for r in rs.values()) / n)


def ganho(a, d):
    """quantas realizacoes do direto valem uma da AI (mesmo CoV), com a P_f da AI"""
    covd1 = math.sqrt((1 - a["p"]) / a["p"])  # CoV do direto com N = 1
    return covd1 ** 2 / (a["cov"] ** 2 * a["n"])


def desvios(a, d):
    sdd = math.sqrt(d["p"] * (1 - d["p"]) / d["n"])
    return (a["p"] - d["p"]) / math.sqrt(a["sd"] ** 2 + sdd ** 2)


def form(pasta):
    import json
    j = json.load(open(f"{pasta}/form.json"))
    return dict(beta=j["beta"], bc=j.get("beta_corrigido", j["beta"]), pf=phi(-j.get("beta_corrigido", j["beta"])),
                n=j["avaliacoes"], it=j.get("iteracoes", len(j["historico"])), a2c=j["alfa2_c"], a2p=j["alfa2_phi"])


def pct(x, d=2):
    return f"{100 * x:.{d}f}".replace(".", "{,}")


def num(x, d=2):
    return f"{x:.{d}f}".replace(".", "{,}")


def barras(arq, d, casos, titulo, ref=None):
    fig, ax = plt.subplots(figsize=(6.4, 2.6))
    nomes = ["direto"] + [c[0] for c in casos]
    vals = [d] + [c[1] for c in casos]
    for i, v in enumerate(vals):
        ax.errorbar(100 * v["p"], i, xerr=[[100 * (v["p"] - v["lo"])], [100 * (v["hi"] - v["p"])]], fmt="o",
                    color="k" if i == 0 else "tab:blue", capsize=3)
    ax.axvspan(100 * d["lo"], 100 * d["hi"], color="0.85", zorder=0)
    if ref:
        ax.axvline(100 * ref, ls="--", color="tab:red", lw=1)
        titulo += "  (tracejado vermelho: Cho, 2010)"
    ax.set_yticks(range(len(vals)))
    ax.set_yticklabels(nomes, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("$P_f$ (%) com IC 95%")
    ax.set_title(titulo, fontsize=9)
    fig.tight_layout()
    fig.savefig(f"{AQUI}/figs/{arq}")
    plt.close(fig)


# ===================== etapa 1: SRM antigo contra o novo =====================
def pares(pa, pb, t=None):
    A, B = le(pa, t), le(pb, t)
    c = sorted(set(A) & set(B))
    fa = {n for n in c if float(A[n]["FS_min"]) < 1}
    fb = {n for n in c if float(B[n]["FS_min"]) < 1}
    dif = [float(B[n]["FS_min"]) - float(A[n]["FS_min"]) for n in c]
    m = sum(dif) / len(dif)
    sd = math.sqrt(sum((x - m) ** 2 for x in dif) / (len(dif) - 1))
    perto = [n for n in c if abs(float(A[n]["FS_min"]) - 1) <= 0.03]
    tA = sum(float(A[n]["s_fator"]) for n in c) / len(c)
    tB = sum(float(B[n]["s_fator"]) for n in c) / len(c)
    tAp = sum(float(A[n]["s_fator"]) for n in perto) / len(perto)
    tBp = sum(float(B[n]["s_fator"]) for n in perto) / len(perto)
    hist = [(sum(1 for n in c if 1 + k / 1000 <= float(A[n]["FS_min"]) < 1 + (k + 1) / 1000),
             sum(1 for n in c if 1 + k / 1000 <= float(B[n]["FS_min"]) < 1 + (k + 1) / 1000)) for k in range(-10, 10)]
    return dict(n=len(c), fa=len(fa), fb=len(fb), n01=len(fb - fa), n10=len(fa - fb), m=m, sd=sd,
                mx=max(abs(x) for x in dif), tA=tA, tB=tB, tAp=tAp, tBp=tBp, np=len(perto), hist=hist)


def mcnemar(n01, n10):
    m = n01 + n10
    if m == 0:
        return 1.
    k = min(n01, n10)
    return min(1., 2 * sum(math.comb(m, i) for i in range(k + 1)) / 2 ** m)


cho = pares(f"{E1}/r05_direto", f"{E1}/r05_antigo")  # A = novo, B = antigo
chuva = pares(f"{E1}/chuva_antigo", f"{E1}/chuva_novo")  # A = antigo, B = novo
# o buraco nos solos (campanha antiga, SRM antigo): f50, todos os instantes
solos = []
for a in glob.glob(f"{RUNS}/fermi_chuva_solos/chuva_solos/f50/mc_*_*.csv"):
    if re.search(r"/mc_\d+_\d+\.csv$", a):
        solos += [float(r["FS_min"]) for r in csv.DictReader(open(a)) if int(r["real"]) >= 0 and r["situacao"] == "ok"]
s_abaixo = sum(0.995 <= x < 1 for x in solos)
s_acima = sum(1 <= x < 1.005 for x in solos)

fig, axs = plt.subplots(1, 2, figsize=(7.2, 2.6), sharey=True)
ks = [1 + (k + 0.5) / 1000 for k in range(-10, 10)]
for ax, (rot, h) in zip(axs, (("chuva f50, 24 h (600 pares)", chuva["hist"]), ("Cho2010, r = -0,5 (1000 pares)", cho["hist"]))):
    antigo = [x[0] for x in h] if rot.startswith("chuva") else [x[1] for x in h]
    novo = [x[1] for x in h] if rot.startswith("chuva") else [x[0] for x in h]
    ax.bar([k - 0.00022 for k in ks], antigo, width=0.00042, label="SRM antigo ($F_0 = 1$)", color="0.6")
    ax.bar([k + 0.00022 for k in ks], novo, width=0.00042, label="SRM novo ($F_0 = 0{,}8$)", color="tab:blue")
    ax.axvline(1, color="k", lw=0.8)
    ax.set_title(rot, fontsize=9)
    ax.set_xlabel("FS (faixas de 0,001)")
axs[0].set_ylabel("realizações")
axs[0].legend(fontsize=7)
fig.tight_layout()
fig.savefig(f"{AQUI}/figs/srm_hist.pdf")
plt.close(fig)

# ===================== etapa 2: r = -0,5 =====================
d05 = direto(f"{E1}/r05_direto")
a05 = [("AI pela banda (entropia cruzada)", ai(f"{E1}/r05_ai")),
       ("AI pelo FORM $M = 10$", ai(f"{E1}/r05_ai_form")),
       ("AI pelo FORM $M = 20$", ai(f"{E1}/r05_ai_form20"))]
for s in ("s4", "s5"):  # repeticoes independentes (etapa3), se ja baixadas
    p = f"{E3}/r05_ai_form20_{s}"
    if glob.glob(p + "/mc_*_*.csv"):
        a05.append((f"AI pelo FORM $M = 20$, semente {s[1]}", ai(p)))
f05 = (form(f"{E1}/form_M10"), form(f"{E1}/form_M20"))
barras("pf_r05.pdf", d05, [(re.sub(r"[\$\\]", "", c[0]), c[1]) for c in a05], "Cho2010, r(c, phi) = -0,5", ref=0.0171)

# ===================== etapa 3: r = -0,7 =====================
d07 = direto(f"{E3}/r07_direto")
a07 = [("AI pela banda (entropia cruzada)", ai(f"{E1}/r07_ai")),
       ("AI pelo FORM $M = 10$", ai(f"{E3}/r07_ai_form10")),
       ("AI pelo FORM $M = 20$", ai(f"{E3}/r07_ai_form20"))]
for c, rot in (("r07_ai_form20_s4", "AI pelo FORM $M = 20$, semente 4"), ("r07_ai_form20_s5", "AI pelo FORM $M = 20$, semente 5"),
               ("r07_ai_banda_s4", "AI pela banda, semente 4")):
    if glob.glob(f"{E3}/{c}/mc_*_*.csv"):
        a07.append((rot, ai(f"{E3}/{c}")))
f07 = (form(f"{E3}/form_M10"), form(f"{E3}/form_M20"))
barras("pf_r07.pdf", d07, [(re.sub(r"[\$\\]", "", c[0]), c[1]) for c in a07], "Cho2010, r(c, phi) = -0,7", ref=0.004)

# o r07 da etapa1 (1000) contra as 1000 primeiras do r07 da etapa3: a mesma realizacao, versoes do programa diferentes
r07a, r07b = le(f"{E1}/r07_direto"), le(f"{E3}/r07_direto")
cm = sorted(set(r07a) & set(r07b))
rep_max = max(abs(float(r07a[n]["FS_min"]) - float(r07b[n]["FS_min"])) for n in cm)

# ===================== diagnostico de outro modo de ruptura =====================
tmp = f"{RUNS}/cho2010/diagnostico/modo2_rel"  # o modo linear do Cho2010 (sem elementos finitos)
diag = {}
for nome, pasta_d, ustar, rc, n in (("r05", f"{E1}/r05_direto", f"{E1}/form_M20/ustar.txt", "-0.5", d05["n"]),
                                   ("r07", f"{E3}/r07_direto", f"{E3}/form_M20/ustar.txt", "-0.7", d07["n"])):
    os.makedirs(f"{tmp}/{nome}", exist_ok=True)
    n = max(le(pasta_d)) + 1  # o maior indice (os blocos em andamento deixam buracos)
    lin = le(f"{tmp}/{nome}")
    if len(lin) < n:  # o direto cresceu: refaz (sem elementos finitos, segundos)
        for a_ in glob.glob(f"{tmp}/{nome}/mc_*"):
            os.remove(a_)
        subprocess.run([CHO2010, f"n={n}", "linear=1", f"alfa={ustar}", f"rcphi={rc}", f"saida={tmp}/{nome}/mc"],
                       check=True, stdout=subprocess.DEVNULL)
        lin = le(f"{tmp}/{nome}")
    t = {k: 1 - math.log(float(r["FS_min"])) for k, r in lin.items()}
    dd = le(pasta_d)
    fal = [k for k, r in dd.items() if float(r["FS_min"]) < 1]
    seg = [k for k, r in dd.items() if float(r["FS_min"]) >= 1]
    diag[nome] = dict(tf=[t[k] for k in fal], ts=[t[k] for k in seg], tmin=min(t[k] for k in fal), nf=len(fal),
                      prof=sum(float(dd[k]["banda_prof"]) for k in fal) / len(fal),
                      xc=sum(float(dd[k]["banda_xc"]) for k in fal) / len(fal))
fig, axs = plt.subplots(1, 2, figsize=(7.2, 2.5))
for ax, (nome, rot, beta) in zip(axs, (("r05", "r = -0,5", f05[1]["bc"]), ("r07", "r = -0,7", f07[1]["bc"]))):
    g = diag[nome]
    ax.hist(g["ts"], bins=50, density=True, color="0.75", label="seguras")
    ax.hist(g["tf"], bins=15, density=True, color="tab:red", alpha=0.7, label="falhas")
    ax.axvline(beta, color="k", ls="--", lw=1)
    ax.text(beta, ax.get_ylim()[1] * 0.9, r" $\beta$ (FORM)", fontsize=8)
    ax.set_title(f"{rot}: t = alfa . xi nas realizações do direto", fontsize=8)
    ax.set_xlabel("t")
axs[0].legend(fontsize=7)
fig.tight_layout()
fig.savefig(f"{AQUI}/figs/modo.pdf")
plt.close(fig)

# ===================== etapa 4: referencia da chuva =====================
dch = direto(f"{E3}/chuva31_24h", 24)
ach = [("AI pela banda, semente 1", ai(f"{E3}/chuva31_ai", 24))]
if glob.glob(f"{E3}/chuva31_ai_s2/mc_*_*.csv"):
    ach.append(("AI pela banda, semente 2", ai(f"{E3}/chuva31_ai_s2", 24)))

# ===================== teste linear =====================
lin_txt = None
if os.path.exists(f"{LIN}/resultado.txt"):
    lin_txt = open(f"{LIN}/resultado.txt").read()

# ===================== teste linear (tabela) =====================
def linear_tab():
    if not os.path.isdir(LIN):
        return None
    lin = []
    for b in (2, 3, 4):
        ex = phi(-b)
        lin_d = None
        try:
            txt = open(f"{LIN}/b{b}/direto.log").read()
            m = re.search(r"P\(FS < 1\) contada = ([\d.]+)%", txt)
            lin_d = float(m.group(1)) / 100 if m else None
        except OSError:
            pass
        linha = {"b": b, "ex": ex, "d": lin_d}
        for a in ("ai20", "ai10"):
            w = []
            for f in glob.glob(f"{LIN}/b{b}/{a}/mc_*_*.csv"):
                if "mapa" in f:
                    continue
                w += [math.exp(float(r["log_peso"])) * (float(r["FS_min"]) < 1) for r in csv.DictReader(open(f))]
            n = len(w)
            pp = sum(w) / n
            sd = math.sqrt((sum(x * x for x in w) / n - pp * pp) / (n - 1))
            linha[a] = (pp, sd / pp, (pp - ex) / sd)
        lin.append(linha)
    a10 = [float(x) for x in open(f"{LIN}/alfa_M10.txt").read().split()]
    a20 = [float(x) for x in open(f"{LIN}/alfa_M20.txt").read().split()]
    cosang = sum(a10[k] * a20[k] for k in range(10)) + sum(a10[10 + k] * a20[20 + k] for k in range(10))
    return lin, cosang


LINTAB = linear_tab()


# ===================== tamanho do piloto: sub-pilotos disjuntos do direto =====================
def projecoes(nome, n, ustar, banda, rc):
    """proj_c e proj_phi das realizacoes 0..n-1 (o modo linear do Cho2010, sem elementos finitos, com banda=)"""
    pasta = f"{tmp}/proj_{nome}"
    os.makedirs(pasta, exist_ok=True)
    if not glob.glob(f"{pasta}/mc_*_*.csv"):
        subprocess.run([CHO2010, f"n={n}", "linear=1", f"alfa={ustar}", f"banda={banda}", f"rcphi={rc}",
                        f"saida={pasta}/mc"], check=True, stdout=subprocess.DEVNULL)
    return {k: (float(r["proj_c"]), float(r["proj_phi"])) for k, r in le(pasta).items()}


def subpilotos(dd, proj, nmax):
    out = []
    for m in (100, 250, 500, 1000):
        ss = []
        for c0 in range(0, nmax - m + 1, m):
            grupo = [k for k in range(c0, c0 + m) if k in dd]
            fs = sorted(grupo, key=lambda k: float(dd[k]["FS_min"]))
            nf = sum(float(dd[k]["FS_min"]) < 1 for k in grupo)
            nel = max(nf, round(0.1 * len(grupo)))
            el = fs[:nel]
            sc = sum(proj[k][0] for k in el) / nel
            sp = sum(proj[k][1] for k in el) / nel
            ss.append((sc, sp, sc * sc + sp * sp))
        mu = [sum(x[i] for x in ss) / len(ss) for i in range(3)]
        sd = [math.sqrt(sum((x[i] - mu[i]) ** 2 for x in ss) / (len(ss) - 1)) if len(ss) > 1 else 0. for i in range(3)]
        out.append((m, len(ss), mu, sd, ss[0]))
    return out


SUB = {
    "r05": subpilotos(le(f"{E1}/r05_direto"), projecoes("r05", 4000, f"{E1}/form_M20/ustar.txt",
                                                         f"{E1}/r05_direto/mc_real-1/mecanismo.csv", "-0.5"), 4000),
    "r07": subpilotos(le(f"{E3}/r07_direto"), projecoes("r07", 6000, f"{E3}/form_M20/ustar.txt",
                                                         f"{E3}/r07_direto/mc_real-1/mecanismo.csv", "-0.7"), 6000),
}


# ===================== custo total para CoV(P_f) = 10% =====================
def custo_tab(d, casos, extra):
    """[(metodo, N para CoV 10%, custo extra, total)]"""
    nd = (1 - d["p"]) / (d["p"] * 0.01)
    linhas = [("Monte Carlo direto", nd, 0, nd)]
    for (rot, a), ex in zip(casos, extra):
        na = a["n"] * (a["cov"] / 0.10) ** 2
        linhas.append((rot, na, ex, na + ex))
    return linhas


CUSTO05 = custo_tab(d05, a05[:3], [1000, f05[0]["n"], f05[1]["n"]])
CUSTO07 = custo_tab(d07, a07[:3], [1000, f07[0]["n"], f07[1]["n"]])

# ===================== previsao: Var(M = 10) / Var(M = 20) ~ exp(beta10^2 - beta20^2) =====================
PREV = []
for nome, f, a10_, a20_ in (("$r = -0{,}5$", f05, a05[1][1], a05[2][1]), ("$r = -0{,}7$", f07, a07[1][1], a07[2][1])):
    PREV.append((nome, math.exp(f[0]["bc"] ** 2 - f[1]["bc"] ** 2), (a10_["cov"] / a20_["cov"]) ** 2))

# ===================== LaTeX =====================
def tab_simples(cab, linhas, fmt):
    return "\\resizebox{\\textwidth}{!}{%\n" + "\n".join([r"\begin{tabular}{" + fmt + "}", r"\toprule", cab + r" \\", r"\midrule"] +
                     [l + r" \\" for l in linhas] + [r"\bottomrule", r"\end{tabular}"]) + "}"


# rodadas do PC: AI com o piloto e o FORM na malha Tri(0) e com o gradiente modal
G = f"{RUNS}/cho2010/grossa/campanha"
pc = []
for c, rot in (("ai_banda_tri0", "AI pela banda, piloto na Tri(0) [PC]"), ("ai_banda_tri1", "AI pela banda, piloto na Tri(1) [PC]"),
               ("ai_banda_tri2", "AI pela banda, piloto na Tri(2) [PC]"),
               ("ai_form100_tri0", "AI pelo FORM $M = 100$, direção da Tri(0) [PC]"),
               ("ai_form100_modal", "AI pelo FORM $M = 100$, gradiente modal [PC]"),
               ("ai_form200_modal", "AI pelo FORM $M = 200$, gradiente modal [PC]")):
    if glob.glob(f"{G}/{c}/mc_*_*.csv"):
        pc.append((rot, ai(f"{G}/{c}")))
    else:
        pc.append((rot, None))
# as novas do fermi (sementes proprias)
E3 = f"{RUNS}/fermi_etapa3/etapa3"
for c, rot in (("ai_form100_modal_s6", "AI pelo FORM $M = 100$ modal, semente 6 [fermi]"),
               ("ai_form200_modal_s7", "AI pelo FORM $M = 200$ modal, semente 7 [fermi]"),
               ("ai_banda_tri1_s10", "AI pela banda, piloto Tri(1), semente 10 [fermi]")):
    if glob.glob(f"{E3}/{c}/mc_*_*.csv"):
        pc.append((rot, ai(f"{E3}/{c}")))
# campos com a mesma variancia (≈0,93): Nystrom 200 e a analitica truncada em 200
iguais = []
for c, rot in (("ai_nys_form200", "Nyström, 200 termos (ponto do FORM modal Nyström $M = 200$)"),
               ("ai_nys_form20", "Nyström, 200 termos (ponto do FORM modal Nyström $M = 20$)"),
               ("ai_form200_trunc", "analítica truncada em 200 termos (ponto do FORM modal $M = 200$)")):
    iguais.append((rot, ai(f"{G}/{c}") if glob.glob(f"{G}/{c}/mc_*_*.csv") else None))

# a fase 4 (cancelada em 17/11 de manha, com ~65 realizacoes por semente): as sementes parciais com o mesmo
# deslocamento, juntas. Sementes diferentes dao realizacoes diferentes, e o estimador e a media de I w de todas
def ai_juntas(pastas):
    wi, w = [], []
    for p in pastas:
        for r in le(p).values():
            w.append(math.exp(float(r["log_peso"])))
            wi.append(w[-1] * (float(r["FS_min"]) < 1))
    n = len(wi)
    f = [x for x in wi if x > 0]
    if n < 2 or not f:
        return None
    p = sum(wi) / n
    sd = math.sqrt(max(sum(x * x for x in wi) / n - p * p, 0) / (n - 1))
    return dict(n=n, k=len(f), p=p, cov=sd / p, lo=p - 1.96 * sd, hi=p + 1.96 * sd, sd=sd,
                essF=sum(f) ** 2 / sum(x * x for x in f), maxp=max(f) / sum(f), mw=sum(w) / n, pastas=len(pastas))


def existentes(fins):
    return [f"{E3}/{c}" for c in fins if glob.glob(f"{E3}/{c}/mc_*_*.csv")]


F4_KL = existentes([f"ai_form100_modal_s{s}" for s in range(11, 19)])
F4_NY = existentes([f"ai_nys_form200_s{s}" for s in [9] + list(range(11, 19))])
F4_NY20 = existentes(["ai_nys_form20_s8"])
F4 = [(f"AI pelo FORM modal $M = 100$, fase 4: {len(F4_KL)} sementes parciais [fermi]", ai_juntas(F4_KL)),
      (f"Nyström, FORM modal $M = 200$, fase 4: {len(F4_NY)} sementes parciais [fermi]", ai_juntas(F4_NY)),
      ("Nyström, FORM modal $M = 20$, semente 8 parcial [fermi]", ai_juntas(F4_NY20))]
F4 = [(r, a) for r, a in F4 if a]



# media das sementes do FORM M = 20 (independentes)
def media_sementes(casos):
    v = [a for a in casos if a]
    if not v:
        return None
    p = sum(a["p"] for a in v) / len(v)
    sd = math.sqrt(sum(a["sd"] ** 2 for a in v)) / len(v)
    return dict(p=p, sd=sd, k=len(v))


SEM05 = [ai(f"{E1}/r05_ai_form20")] + [ai(f"{E3}/r05_ai_form20_s{k}") for k in (4, 5)
                                                  if glob.glob(f"{E3}/r05_ai_form20_s{k}/mc_*_*.csv")]
SEM07 = [ai(f"{E3}/r07_ai_form20")] + [ai(f"{E3}/r07_ai_form20_s{k}") for k in (4, 5)
                                             if glob.glob(f"{E3}/r07_ai_form20_s{k}/mc_*_*.csv")]
M05, M07 = media_sementes(SEM05), media_sementes(SEM07)
MODAL = [ai(f"{E3}/{c}") for c in ("ai_form100_modal_s6", "ai_form200_modal_s7")
         if glob.glob(f"{E3}/{c}/mc_*_*.csv")]
SEM05T = SEM05 + MODAL
M05T = media_sementes(SEM05T)


def dif_media(m, d):
    return (m["p"] - d["p"]) / math.hypot(m["sd"], d["p"] * d["cov"])


# HIP05 e HIP07: depois do diagnostico (dependem dele)


# FORM com o gradiente modal: as duas bases
def form_json(d):
    if not os.path.exists(f"{d}/form.json"):
        return None
    j = json.load(open(f"{d}/form.json"))
    bc = j.get("beta_corrigido", j["beta"])
    return dict(b=bc, pf=0.5 * math.erfc(bc / math.sqrt(2)), n=j["avaliacoes"])


NY = f"{RUNS}/cho2010/nystrom"
FM = {}
for M in (10, 20, 50, 100, 200):
    an = f"{NY}/form_analitica_M{M}" if M <= 50 else f"{RUNS}/cho2010/grossa/form_M{M}_modal"
    FM[M] = (form_json(an), form_json(f"{NY}/form_nystrom_M{M}"))
# variancia retida da base de Nystrom
VAR_NY = {10: 0.655, 20: 0.769, 50: 0.876, 100: 0.933, 200: 0.970}

# grafico: beta(M) das duas bases e as diferencas centrais
fig, ax = plt.subplots(figsize=(6, 3.2))
Ms = sorted(FM)
for i, (rot, cor) in enumerate((("base analítica (retângulo)", "tab:blue"), ("base de Nyström (talude)", "tab:red"))):
    xs = [M for M in Ms if FM[M][i]]
    ax.plot(xs, [FM[M][i]["b"] for M in xs], "o-", color=cor, label=rot + ", gradiente modal")
ax.plot([10, 20], [f05[0]["bc"], f05[1]["bc"]], "s--", color="0.4", label="analítica, diferenças centrais")
ax.set_xscale("log")
ax.set_xlabel("M (termos por campo)")
ax.set_ylabel(r"$\beta$ (FORM)")
ax.grid(alpha=0.3)
ax.legend(fontsize=7)
fig.tight_layout()
fig.savefig(f"{AQUI}/figs/beta_M.pdf")
plt.close(fig)

# ===================== o diagnostico da AI pelo FORM com o direto (Projects/Cho2010/diagnostico_ai.py) =====================
DGM = runpy.run_path(f"{NEOPZ}/Projects/Cho2010/diagnostico_ai.py", run_name="diagnostico")
DG = DGM["roda"]()
PP = DGM["ponto_de_projeto"]()
DD5, DD7 = DG["diretos"]["r05"], DG["diretos"]["r07"]
TTd = DG["TT"]


def acha(fim):
    """a AI do diagnostico cuja pasta termina em fim"""
    v = [a for a in DG["ai"] if a["pasta"].endswith("/" + fim)]
    return v[0] if v else None


# os campos no ponto de projeto, nas duas bases
fig, axs = plt.subplots(2, 2, figsize=(7.4, 4.3), sharex=True, sharey=True)
for i, (b, rot) in enumerate((("kl", "KL analítica"), ("nys", "Nyström"))):
    for j, (q, un, cmap) in enumerate((("c", "kPa", "viridis"), ("phi", "graus", "magma"))):
        ax = axs[i, j]
        vmin = min(np.nanmin(PP[bb][q]) for bb in ("kl", "nys"))
        vmax = max(np.nanmax(PP[bb][q]) for bb in ("kl", "nys"))
        im = ax.pcolormesh(PP["x"], PP["y"], PP[b][q], shading="nearest", vmin=vmin, vmax=vmax, cmap=cmap,
                           rasterized=True, linewidth=0)
        ax.plot([0, 30, 30, 20, 10, 0, 0], [0, 0, 5, 5, 15, 15, 0], "k-", lw=0.6)
        ax.set_aspect("equal")
        ax.set_title(f"{rot}: " + ("$c'$" if q == "c" else r"$\varphi'$") + f" ({un}); " + r"$\beta$ = " +
                     f"{PP[b]['beta']:.3f}".replace(".", ","), fontsize=8)
        fig.colorbar(im, ax=ax, shrink=0.8)
for ax in axs[1]:
    ax.set_xlabel("x (m)")
for ax in axs[:, 0]:
    ax.set_ylabel("y (m)")
fig.tight_layout()
fig.savefig(f"{AQUI}/figs/ponto_projeto.pdf")
plt.close(fig)

# t = alfa . xi: a chance de falha do direto em cada faixa de t, e as densidades do direto e da AI
fig, axs = plt.subplots(1, 2, figsize=(7.4, 2.7), sharey=True)
for ax, (fim, rot) in zip(axs, (("ai_form100_modal", "KL analítica, FORM modal $M = 100$"),
                                ("ai_nys_form200", "Nyström, FORM modal $M = 200$"))):
    a = acha(fim)
    h = a["hist"]
    bd = np.array(h["bordas"]); cen = 0.5 * (bd[1:] + bd[:-1])
    tod, fal = np.array(h["todas"]), np.array(h["falhas"])
    ok = tod > 0
    ax.bar(cen[ok], fal[ok] / tod[ok], width=0.25, color="tab:red", alpha=0.6, label="P(falha | t), direto")
    tt = np.linspace(-4, 5, 300)
    ax.plot(tt, np.exp(-tt ** 2 / 2), color="0.4", lw=1.2, label="densidade do direto (escala)")
    ax.plot(tt, np.exp(-(tt - a["mu"]) ** 2 / 2), color="tab:blue", lw=1.2, label="densidade da AI (escala)")
    ax.set_title(rot + f"; t mínimo das falhas: {a['t_min']:.2f}".replace(".", ","), fontsize=8)
    ax.set_xlabel(r"$t = \xi\cdot\xi^*/\beta$")
    ax.grid(alpha=0.3)
axs[0].legend(fontsize=6.5, loc="upper left")
fig.tight_layout()
fig.savefig(f"{AQUI}/figs/vies_t.pdf")
plt.close(fig)


# as sementes independentes: uma AI por semente
def sementes(fins):
    v = [acha(f) for f in fins]
    return [a for a in v if a and a["n"] >= 1000]


KL05 = sementes(["ai_form100_modal", "r05_ai_form20_s4", "r05_ai_form20_s5", "ai_form100_modal_s6", "ai_form200_modal_s7"] +
                [f"ai_form100_modal_s{s}" for s in range(11, 19)])
NY05 = sementes(["ai_nys_form200", "ai_nys_form20_s8", "ai_nys_form200_s9"] + [f"ai_nys_form200_s{s}" for s in range(11, 19)])
KL07 = sementes(["r07_ai_form20", "r07_ai_form20_s4", "r07_ai_form20_s5"])


def resume(v, d):
    k = len(v)
    p, pc = sum(a["p"] for a in v) / k, sum(a["p_cv"] for a in v) / k
    sd, sdc = math.sqrt(sum(a["sd"] ** 2 for a in v)) / k, math.sqrt(sum(a["sd_cv"] ** 2 for a in v)) / k
    emp = (math.sqrt(sum((a["p"] - p) ** 2 for a in v) / (k - 1)) if k > 1 else float("nan"))
    rms = math.sqrt(sum(a["sd"] ** 2 for a in v) / k)
    return dict(k=k, p=p, sd=sd, pc=pc, sdc=sdc, emp=emp, rms=rms, dif=(p - d["p"]) / math.hypot(sd, d["sd"]),
                difc=(pc - d["cv"]["p"]) / math.hypot(sdc, d["cv"]["sd"]))


RS = [("$r = -0{,}5$, KL (modal; $M = 20$ em s4 e s5)", KL05, DD5),
      ("$r = -0{,}5$, Nyström (modal)", NY05, DD5), ("$r = -0{,}7$, KL ($M = 20$)", KL07, DD7)]
RS = [(rot, v, d, resume(v, d)) for rot, v, d in RS if v]

fig, axs = plt.subplots(1, 2, figsize=(7.4, 3.0), gridspec_kw=dict(width_ratios=[3, 1]))
for ax, grupos, d in ((axs[0], [("KL", KL05, "tab:blue"), ("Nyström", NY05, "tab:red")], DD5), (axs[1], [("KL", KL07, "tab:blue")], DD7)):
    x0 = 0
    for nome, v, cor in grupos:
        xs = np.arange(x0, x0 + len(v))
        ax.errorbar(xs - 0.12, [100 * a["p"] for a in v], yerr=[200 * a["sd"] for a in v], fmt="o", mfc="none", color=cor,
                    ms=4, lw=0.8, capsize=2, label=f"{nome}, AI")
        ax.errorbar(xs + 0.12, [100 * a["p_cv"] for a in v], yerr=[200 * a["sd_cv"] for a in v], fmt="o", color=cor, ms=4,
                    lw=0.8, capsize=2, label=f"{nome}, AI corrigida")
        x0 += len(v) + 1
    ax.axhspan(100 * (d["p"] - 2 * d["sd"]), 100 * (d["p"] + 2 * d["sd"]), color="0.85", zorder=0, label="direto, IC 95%")
    ax.axhspan(100 * (d["cv"]["p"] - 2 * d["cv"]["sd"]), 100 * (d["cv"]["p"] + 2 * d["cv"]["sd"]), color="0.65", zorder=0,
               label="direto corrigido, IC 95%")
    ax.axhline(100 * d["p"], color="k", lw=0.6)
    lab, pos, x0 = [], [], 0
    for nome, v, cor in grupos:
        pos += list(range(x0, x0 + len(v)))
        lab += [f"s{a['semente']}" for a in v]
        x0 += len(v) + 1
    ax.set_xticks(pos)
    ax.set_xticklabels(lab, fontsize=6, rotation=90)
    ax.set_ylabel("$P_f$ (%)")
    ax.grid(axis="y", alpha=0.3)
axs[0].set_title("r = -0,5: uma AI por semente (KL à esquerda, Nyström à direita)", fontsize=8)
axs[1].set_title("r = -0,7", fontsize=8)
axs[0].legend(fontsize=6, ncol=2, loc="upper left")
fig.tight_layout()
fig.savefig(f"{AQUI}/figs/vies_sementes.pdf")
plt.close(fig)



def linha_ai(rot, a, d):
    if a is None:
        return f"{rot} & --- & na fila & --- & --- & --- & --- & --- \\\\"
    return (f"{rot} & {a['n']} & {pct(a['p'])} & {pct(a['lo'])}--{pct(a['hi'])} & {num(100 * a['cov'], 1)} & "
            f"{num(desvios(a, d), 1)} & {a['essF']:.0f} & {num(ganho(a, d), 1)} \\\\")


def valor(cel):
    """o primeiro numero de uma celula ('1{,}70', '-1{,}3', '0{,}8\\,\\% (0{,}5)'); None se a celula nao comeca por um"""
    m = re.match(r"\s*\$?(-?\d+(?:(?:\{,\}|,)\d+)?)", cel)
    return float(m.group(1).replace("{,}", ".").replace(",", ".")) if m else None


def pinta(tab, grupo, regra):
    """entre as celulas (linha, coluna) do grupo, o melhor valor fica azul e o pior vermelho. regra 'min': o menor e o
    melhor; 'max': o maior; 'abs': o menor em modulo. A comparacao usa o valor impresso, e os empates ganham a mesma cor."""
    vals = {}
    for i, c in grupo:
        v = valor(tab[i][c])
        if v is not None:
            vals[(i, c)] = abs(v) if regra == "abs" else v
    if len(set(vals.values())) < 2:
        return
    lo, hi = min(vals.values()), max(vals.values())
    melhor, pior = (hi, lo) if regra == "max" else (lo, hi)
    for (i, c), v in vals.items():
        if v in (melhor, pior):
            tab[i][c] = "\\textcolor{" + ("azul" if v == melhor else "verm") + "}{\\textbf{" + tab[i][c].strip() + "}}"


def separa(linhas):
    tab = []
    for l in linhas:
        l = l.rstrip()
        assert l.endswith("\\\\"), l
        tab.append(l[:-2].rstrip().split(" & "))
    return tab


def junta(tab):
    return [" & ".join(c) + " \\\\" for c in tab]


def destaca_colunas(linhas, regras, entra):
    """regras = {coluna: 'min' | 'max' | 'abs'}; entra(i, celulas) diz se a linha entra na comparacao"""
    tab = separa(linhas)
    for c, regra in regras.items():
        pinta(tab, [(i, c) for i, cels in enumerate(tab) if entra(i, cels)], regra)
    return junta(tab)


def destaca_linhas(linhas, colunas, regra):
    """a comparacao e entre as colunas, linha a linha (os metodos estao nas colunas)"""
    tab = separa(linhas)
    for i in range(len(tab)):
        pinta(tab, [(i, c) for c in colunas], regra)
    return junta(tab)


def completa(i, cels):
    """as AI completas: o direto (linha 0) e a referencia, e as parciais e as da fila ficam de fora"""
    return i > 0 and cels[1].strip().isdigit() and int(cels[1]) >= 1000


def destaca_ai(linhas):
    """tabelas de AI: CoV (o menor), dif. (o menor em modulo) e ganho (o maior)"""
    return destaca_colunas(linhas, {4: "min", 5: "abs", 7: "max"}, completa)


def tabela(cab, linhas, fmt):
    return ("\\resizebox{\\ifdim\\width>\\textwidth\\textwidth\\else\\width\\fi}{!}{%\n\\begin{tabular}{" + fmt + "}\n\\toprule\n" + cab + " \\\\\n\\midrule\n" +
            "\n".join(linhas) + "\n\\bottomrule\n\\end{tabular}}")


CAB_AI = "Método & $N$ & $P_f$ (\\%) & IC 95\\,\\% (\\%) & CoV (\\%) & dif. (desvios) & ESS$_F$ & ganho"


def linha_direto(d):
    return (f"Monte Carlo direto (fermi) & {d['n']} & {pct(d['p'])} & {pct(d['lo'])}--{pct(d['hi'])} & "
            f"{num(100 * d['cov'], 1)} & --- & --- & 1 \\\\")


def linhas_form(fs):
    return [f"FORM, $M = {M}$ ($\\beta = {num(f['bc'], 3)}$) & {f['n']} análises & {pct(f['pf'])} & --- & --- & --- & "
            f"--- & --- \\\\" for M, f in zip((10, 20), fs) if f]


T_R05 = tabela(CAB_AI, destaca_ai([linha_direto(d05)] + [linha_ai(r, a, d05) for r, a in a05 + pc + F4]) +
               linhas_form(f05), "lrrrrrrr")
T_R07 = tabela(CAB_AI, destaca_ai([linha_direto(d07)] + [linha_ai(r, a, d07) for r, a in a07]) + linhas_form(f07),
               "lrrrrrrr")
T_CH = tabela(CAB_AI, destaca_ai([linha_direto(dch)] + [linha_ai(r, a, dch) for r, a in ach]), "lrrrrrrr")
T_FORM = tabela("$M$ & variância retida (Nyström) & $\\beta$ analítica & $\\beta$ Nyström & análises por FORM",
                [f"{M} & {num(VAR_NY[M], 3)} & " + (num(FM[M][0]['b'], 3) if FM[M][0] else "---") + " & " +
                 (num(FM[M][1]['b'], 3) if FM[M][1] else "---") + " & " +
                 f"{(FM[M][1] or FM[M][0] or {'n': '---'})['n']} \\\\" for M in Ms], "rrrrr")


def tabela_vies(r, d, dec):
    """as AI pelo ponto de projeto de um r: o diagnostico pelo direto e a estimativa com a variavel de controle"""
    lin = [f"Monte Carlo direto & --- & 2010 & {d['n']} & {pct(d['p'], dec)} & {num(100 * d['sd'] / d['p'], 1)} & --- & --- & "
           f"{pct(d['cv']['p'], dec)} & {num(100 * d['cv']['sd'] / d['cv']['p'], 1)} \\\\"]
    for a in DG["ai"]:
        if a["r"] != r or a["n"] < 1000:  # as completas
            continue
        lin.append(f"{a['rotulo']} & {'KL' if a['base'] == 'kl' else 'Nyström'} & {a['semente']} & {a['n']} & "
                   f"{pct(a['p'], dec)} & {num(100 * a['sd'] / a['p'], 1)} & {num(100 * a['cov_verdadeiro'], 1)} & "
                   f"{num(a['t_min'], 2)} & {pct(a['p_cv'], dec)} & {num(100 * a['sd_cv'] / a['p_cv'], 1)} \\\\")
    return tabela("AI & base & semente & $N$ & $P_f$ (\\%) & CoV (\\%) & CoV verd. (\\%) & $t_{\\min}$ & $P_f$ corr. (\\%) & "
                  "CoV corr. (\\%)",
                  destaca_colunas(lin, {5: "min", 9: "min"}, lambda i, c: i > 0 and c[3].strip().isdigit() and int(c[3]) >= 1000),
                  "llrrrrrrrr")


T_VIES05, T_VIES07 = tabela_vies("r05", DD5, 2), tabela_vies("r07", DD7, 3)
T_SEM = tabela("Grupo & sementes & $P_f$ média (\\%) & dif. & $P_f$ corr. média (\\%) & dif. & desvio entre sementes / estimado",
               [f"{rot} & {s['k']} & {pct(s['p'], 3)} $\\pm$ {pct(s['sd'], 3)} & ${num(s['dif'], 1)}$ & "
                f"{pct(s['pc'], 3)} $\\pm$ {pct(s['sdc'], 3)} & ${num(s['difc'], 1)}$ & " +
                (f"{pct(s['emp'], 3)} / {pct(s['rms'], 3)}" if s["k"] > 2 else "---") + " \\\\" for rot, v, d, s in RS],
               "lrrrrrr")

# o mapa dos documentos: a mesma lista do indice do README (../documentos.py)


def ponto(arq):
    """beta, a fracao de alfa^2 em phi', os 3 maiores termos e os 10 primeiros modos de cada campo"""
    v = np.loadtxt(arq)
    k = len(v) // 2
    a2 = v ** 2 / (v ** 2).sum()
    top = np.argsort(-a2)[:3]
    txt = "os maiores termos são " + ", ".join(("$c'_{" if i < k else "$\\varphi'_{") + f"{i % k + 1}" + "}$ (" + num(a2[i], 2) + ")"
                                             for i in top)
    return dict(phi=a2[k:].sum(), txt=txt, dez=a2[:10].sum() + a2[k:k + 10].sum())


UKL, UNY = ponto(f"{RUNS}/cho2010/grossa/ustar_M100_modal.txt"), ponto(f"{NY}/ustar_nys_M200.txt")


def cient(x):
    m, e = f"{x:.0e}".split("e")
    return f"{m}\\times10^{{{int(e)}}}"


AM, A20, AT0 = acha("ai_form100_modal"), acha("r05_ai_form20"), acha("ai_form100_tri0")
AN20, AN200 = acha("ai_nys_form20"), acha("ai_nys_form200")
A07 = [acha(f) for f in ("r07_ai_form20", "r07_ai_form20_s4", "r07_ai_form20_s5")]
A0710 = acha("r07_ai_form10")
GANHOS = [(a["sd"] / a["sd_cv"]) ** 2 for a in DG["ai"] if a["r"] == "r05" and a["n"] >= 1000]


# o gradiente modal do GIM (PC, 17/11): Projects/Cho2010/valida_gradiente.py e form_kl.py metodo=gim
def le_json(a):
    with open(a) as f:
        return json.load(f)


def angulo(a, b):
    a, b = np.array(a), np.array(b)
    return float(np.degrees(np.arccos(np.clip(a @ b / np.linalg.norm(a) / np.linalg.norm(b), -1, 1))))


GIMD = f"{RUNS}/cho2010/gim_gradiente"
GV0, GVU = le_json(f"{GIMD}/valida_zero/valida.json"), le_json(f"{GIMD}/valida_ustar_srm/valida.json")
SRM_CD = le_json(f"{RUNS}/cho2010/form_M10_r05/form_parcial.json")["historico"][-1]["grad"]
SRM_MOD = le_json(f"{GIMD}/srm_modal_ustar/srm_modal.json")["grad_modal"]
FGIM = {M: le_json(f"{GIMD}/form_gim_M{M}_modal/form.json") for M in (10, 20, 100)}
DOM0, DOMU = GV0["dominante"]["phi0"], GVU["dominante"]["phi0"]  # [modal, dif.]
GIM_TXT = (r"""  \item \textbf{O gradiente modal do GIM} (\path{Projects/Cho2010/valida_gradiente.py}), no PC: fica a """ +
           num(GV0["angulo_graus"], 1) + r"""$^\circ$ das diferenças centrais em $\bm\xi = \bm 0$ ($\Gamma = """ +
           num(GV0["fator_modal"], 3) + r"""$) e a """ + num(GVU["angulo_graus"], 1) + r"""$^\circ$ no ponto de projeto
    do SRM com $M = 10$, onde o modal do SRM fica a """ + num(angulo(SRM_MOD, SRM_CD), 1) + r"""$^\circ$. O termo
    dominante ($\varphi'_0$) sai """ + num(100 * (DOM0[0] / DOM0[1] - 1), 1) + r""" e """ +
           num(100 * (DOMU[0] / DOMU[1] - 1), 1) + r"""\,\% acima. Nesse ponto, $\Gamma = """ + num(GVU["fator_modal"], 3) +
           r"""$, e as diferenças centrais do GIM e do SRM ficam a """ + num(angulo(GVU["grad_dif"], SRM_CD), 1) +
           r"""$^\circ$ uma da outra: as duas são normais à mesma superfície de falha ($\Gamma < 1$ se e só se
    $\mathrm{FS} < 1$). O FORM do GIM com o gradiente modal converge em """ + str(FGIM[10]["avaliacoes"]) +
           r""" análises: $\beta = """ + num(FGIM[10]["beta_corrigido"], 3) + r"""$ ($M = 10$), $""" + num(FGIM[20]["beta_corrigido"], 3) +
           r"""$ ($M = 20$) e $""" + num(FGIM[100]["beta_corrigido"], 3) + r"""$ ($M = 100$), mais perto das diferenças centrais
    do SRM (2,115 e 2,012) que o modal do SRM (Tabela~\ref{tab:gim}).""")

# a variavel de controle com o substituto ajustado na propria amostra (PC, 17/11): Projects/Cho2010/teste_controle.py
TC = le_json(f"{RUNS}/cho2010/diagnostico/teste_controle.json")
LOTES = TC["lotes"]
T_CV_LOTES = tabela("Estimador & média (\\%) & desvio verdadeiro (\\%) & desvio estimado (\\%) & ganho & $\\rho$",
                    [f"{k.replace('realista', 'próprio lote,').replace('ideal', 'fora do lote,').replace('K=', '$K$ = ')} & "
                     f"{num(100 * v['media'], 3)} & {num(100 * v['desvio'], 3)} & {num(100 * v['desvio_estimado'], 3)} & "
                     f"{num(v['ganho'], 2)} & {num(v['rho'], 2) if k != 'P' else '---'} \\\\" for k, v in LOTES.items()],
                    "lrrrrr")
TC_AI = [a for a in TC["ai"] if "Tri(0)" not in a["rotulo"]]
GANHO_K50 = [(a["sd"] / a["realista"]["50"][1]) ** 2 for a in TC_AI]
GANHO_ID = [(a["sd"] / a["ideal"][1]) ** 2 for a in TC_AI]
T_CV_AI = tabela("AI ($r = -0{,}5$) & semente & $P$ (\\%) & própria amostra, $K = 50$ (\\%) & substituto do direto (\\%)",
                 [f"{a['rotulo']}{' (Nyström)' if a['base'] == 'nys' else ''} & {a['semente']} & "
                  f"{num(100 * a['p'], 3)} $\\pm$ {num(100 * a['sd'], 3)} & "
                  f"{num(100 * a['realista']['50'][0], 3)} $\\pm$ {num(100 * a['realista']['50'][1], 3)} & "
                  f"{num(100 * a['ideal'][0], 3)} $\\pm$ {num(100 * a['ideal'][1], 3)} \\\\" for a in TC["ai"]], "lrrrr")
S07 = [s for rot, v, d, s in RS if v is KL07]
S07 = S07[0] if S07 else None

S05 = [s for rot, v, d, s in RS if v is KL05]
S05 = S05[0] if S05 else None
DIF_T0 = (AT0["p"] - DD5["p"]) / math.hypot(AT0["cov_verdadeiro"] * DD5["p"], DD5["sd"])


# o FORM com cada gradiente: o GIM (modal) e o SRM (modal e diferencas centrais), com o beta corrigido
def b_srm_dif(M):
    f = dict(zip((10, 20), f05)).get(M)
    return f"{num(f['bc'], 3)} & {f['n']}" if f else "--- & ---"


T_GIM = tabela("$M$ & GIM, modal & análises & SRM, modal & análises & SRM, diferenças centrais & análises",
               [f"{M} & {num(FGIM[M]['beta_corrigido'], 3)} & {FGIM[M]['avaliacoes']} & " +
                (f"{num(FM[M][0]['b'], 3)} & {FM[M][0]['n']}" if M in FM and FM[M][0] else "--- & ---") + " & " +
                b_srm_dif(M) + " \\\\" for M in sorted(FGIM)], "rrrrrrr")
T_F4 = tabela("AI (fase 4, parcial) & sementes & $N$ & $P_f$ (\\%) & CoV (\\%) & dif. (desvios)",
              [f"{r.split(', fase 4')[0].split(', semente')[0]} & {a['pastas']} & {a['n']} & {pct(a['p'])} $\\pm$ "
               f"{pct(a['sd'])} & {num(100 * a['cov'], 1)} & {num(desvios(a, d05), 1)} \\\\" for r, a in F4], "lrrrrr")
DIFS05 = [abs(desvios(a, d05)) for r, a in a05 + pc if a and a["n"] >= 1000 and "Tri(0)" not in r]

tex = (r"""\documentclass[11pt,a4paper]{article}
\usepackage[utf8]{inputenc}
\DeclareUnicodeCharacter{2248}{\ensuremath{\approx}}
\usepackage[T1]{fontenc}
\usepackage[brazil]{babel}
\usepackage{lmodern}
\usepackage[margin=2cm]{geometry}
\usepackage{amsmath,amssymb,bm}
\usepackage{graphicx,booktabs,float,url}
\usepackage{xcolor}
\definecolor{azul}{RGB}{0,70,200}
\definecolor{verm}{RGB}{200,0,0}
\usepackage{xr-hyper}
\usepackage[hidelinks]{hyperref}
% as secoes da formulacao (docs/confiabilidade_campos/formulacao), citadas como F-<rotulo>: compile a formulacao antes
\externaldocument[F-]{../../../docs/confiabilidade_campos/formulacao/formulacao}
% a data e o commit vem do atualiza.sh (Projects/PlanoDeTrabalho/situacao); compilado a mao, fica so a data
\providecommand{\versao}{17/11/2026}
\title{Plano de trabalho: os resultados das etapas 0 a 3\\\large SRM corrigido, FORM e gradiente modal, as duas bases,
amostragem por importância e variável de controle no Cho (2010), e a referência da chuva}
\date{\versao}
\begin{document}
\maketitle
\tableofcontents

\section*{Resumo}
\begin{itemize}
  \item \textbf{SRM corrigido (etapa 1):} começar o caminho de $F$ em 0,8, e não em 1, quase não muda a $P_f$. No
    Cho2010, as falhas são as mesmas nas """ + f"{cho['n']}" + r""" realizações pareadas. Na chuva, """ +
    f"{chuva['n01']}" + r""" realizações passam a falhar e """ + f"{chuva['n10']}" + r""" deixa de falhar (McNemar
    $p = """ + num(mcnemar(chuva['n01'], chuva['n10']), 2) + r"""$).
  \item \textbf{AI contra o direto (etapas 0, 2 e 3):} as AI completas concordam com o Monte Carlo direto do mesmo modelo,
    a menos de """ + num(max(DIFS05), 1) + r""" desvio com $r = -0{,}5$; só a direção errada da malha Tri(0) fica mais
    longe. O ganho cresce quando a $P_f$ cai: até """ + num(max(ganho(a, d05) for _, a in a05), 0) + r"""$\times$ com
    $P_f \approx 1{,}7\,\%$ e até """ + num(max(ganho(a, d07) for _, a in a07), 0) + r"""$\times$ com
    $P_f \approx 0{,}4\,\%$.
  \item \textbf{FORM sozinho:} com mais termos, $\beta \to 2{,}00$, e a $P_f$ do FORM ($\approx 2{,}3\,\%$) superestima
    a do direto em ≈30\,\%. Ele serve para escolher o deslocamento da AI, e não como estimativa final.
  \item \textbf{As duas bases (etapa 2b):} a de Nyström chega ao $\beta$ convergido com menos termos, e o ponto de projeto
    é o mesmo nas duas (a """ + num(PP["angulo"], 1) + r"""$^\circ$ um do outro).
  \item \textbf{O gradiente modal (etapa 2c):} uma análise por iteração, no lugar de $4M + 1$. Fica a 9 a 11$^\circ$ das
    diferenças centrais no SRM e a 3 a 6$^\circ$ no GIM, e o FORM converge em 3 a 5 análises.
  \item \textbf{Viés:} não há. Nenhuma falha do direto fica fora da região que a AI amostra, a variância verdadeira da AI
    é a que ela estima, e as sementes independentes acertam o direto. A semente 3 é uma tiragem baixa.
  \item \textbf{Variável de controle:} de """ + num(min(GANHOS), 1) + " a " + num(max(GANHOS), 1) + r""" vezes menos
    realizações com o substituto ajustado no direto, e de """ + num(min(GANHO_K50), 1) + " a " +
    num(max(GANHO_K50), 1) + r""" com o ajustado na própria amostra, sem viés.
  \item \textbf{A referência da chuva (etapa 4):} o direto de """ + f"{dch['n']:,}".replace(",", "\\,") + r"""
    realizações dá $P_f = """ + pct(dch['p']) + r"""\,\%$, e a AI pela banda concorda com ele.
\end{itemize}

\section{O que rodou}\label{sec:rodou}
O problema é o talude seco do Cho (2010), Exemplo 2: $c'$ com média 10\,kPa e CoV 0,3; $\varphi'$ com média $30^\circ$ e
CoV 0,2; campos lognormais com $l_h = 20$\,m e $l_v = 2$\,m; correlação cruzada $r = -0{,}5$ ou $-0{,}7$. O FS
determinístico é 1,2044, contra 1,204 do Cho. A formulação está em
\path{docs/confiabilidade_campos/formulacao/formulacao.pdf}, citada aqui como ``formulação, Seção $n$''.
\begin{itemize}
  \item \textbf{No fermi:} as campanhas \texttt{etapa1} (15/11) e \texttt{etapa3} (de 15/11 em diante), com o SRM novo,
    a malha tri2 + 4 níveis e a KL analítica de 2000 termos por campo (ou a de Nyström, de 200). O direto de
    $r = -0{,}7$ continua rodando até 110\,000 realizações. A fase 4 (as sementes 11 a 18 do FORM modal, nas duas
    bases) foi cancelada em 17/11, com ≈65 realizações por semente; as parciais entram juntas
    (Seção~\ref{sec:fase4}).
  \item \textbf{No PC} (\path{~/projects/runs/cho2010}): os pilotos e as AI com piloto em malha grossa, o FORM modal nas
    duas bases, as AI de Nyström, a validação do gradiente do GIM e os testes sem elementos finitos (o diagnóstico de viés
    e a variável de controle).
\end{itemize}
\textbf{Cores nas tabelas:} em cada coluna de desempenho, \textcolor{azul}{\textbf{azul}} é o melhor valor e
\textcolor{verm}{\textbf{vermelho}} o pior. CoV e custo: o menor é o melhor; ganho: o maior; dif.: o menor em módulo. Só
entram as AI completas (1000 realizações): o direto é a referência, e as parciais e as da fila ficam de fora. Atenção: um $|$dif.$|$ abaixo de 2 é compatível com o
acaso; o azul na dif. não quer dizer que o método é mais exato, e o vermelho só é um sinal quando passa de 2.

\section{Etapa 1: o SRM corrigido}\label{sec:srm}
\begin{itemize}
  \item \textbf{O defeito:} com $F_0 = 1$, o peso às vezes convergia num talude com FS logo abaixo de 1, e o FS saía $\ge 1$.
    Na campanha antiga dos solos (f50, todos os instantes), há """ + f"{s_abaixo}" + r""" FS em $[0{,}995;\,1)$ e """ + \
    f"{s_acima}" + r""" em $[1;\,1{,}005)$, e nenhum em $[0{,}999;\,1)$.
  \item \textbf{O teste:} as mesmas realizações com as duas versões, na mesma máquina e com as mesmas opções.
\end{itemize}
\begin{table}[H]\centering\small
""" + tab_simples("Caso & pares & falhas antigo & falhas novo & passam a falhar & deixam & McNemar $p$ & dif. FS: média (dp)",
                  [f"Cho2010, $r = -0{{,}}5$ & {cho['n']} & {cho['fb']} & {cho['fa']} & {cho['n10']} & {cho['n01']} & "
                   f"{num(mcnemar(cho['n01'], cho['n10']), 2)} & {num(-cho['m'], 5)} ({num(cho['sd'], 4)})",
                   f"chuva f50, 36 mm/h, 24 h & {chuva['n']} & {chuva['fa']} & {chuva['fb']} & {chuva['n01']} & {chuva['n10']} & "
                   f"{num(mcnemar(chuva['n01'], chuva['n10']), 2)} & {num(chuva['m'], 5)} ({num(chuva['sd'], 4)})"],
                  "lrrrrrrr") + r"""
\caption{SRM antigo contra o novo, realização a realização.}
\end{table}
\begin{figure}[H]\centering\includegraphics[width=\textwidth]{figs/srm_hist.pdf}
\caption{FS perto de 1 nas duas versões, nas mesmas realizações.}\end{figure}
\begin{itemize}
  \item \textbf{Efeito na $P_f$:} nulo no Cho2010 e de cerca de $+0{,}5$ ponto na chuva, onde a densidade de FS perto de 1 é
    maior.
  \item \textbf{Custo:} no Cho2010, perto de FS = 1, o fator novo leva """ + num(cho['tAp'] / cho['tBp'], 1) + r""" vezes o tempo do
    antigo (""" + f"{cho['np']}" + r""" realizações); na média, """ + num(cho['tA'] / cho['tB'], 2) + r""" vezes. Na chuva, o novo foi
    mais rápido (""" + num(chuva['tB'] / chuva['tA'], 2) + r""" vezes o tempo do antigo).
  \item \textbf{Decisão:} o SRM novo fica, porque é o correto.
  \item \textbf{Conferência de versão:} as 1000 primeiras realizações do direto com $r = -0{,}7$ rodaram duas vezes, com versões
    diferentes do programa (campanhas \texttt{etapa1} e \texttt{etapa3}). A maior diferença de FS foi """ + num(rep_max, 4) + r""".
\end{itemize}

\subsection{O piso do SRM em $F_0 = 0{,}8$}
Nenhuma realização das AI deu FS abaixo de 0,8 (o mínimo foi 0,802). É o mesmo mecanismo do buraco em $F = 1$: um talude
com FS verdadeiro logo abaixo de $F_0$ às vezes converge o peso com $F_0$ e devolve $\approx F_0$. Para a $P_f$ não muda
nada (o FS continua $< 1$), mas a cauda inferior do FS fica distorcida (≈1 realização em 1000). Correção proposta:
recomeçar de $F_0/2$ quando o caminho de $F$ colapsa logo no início.

\section{Etapa 2: FORM e AI com $r = -0{,}5$ ($P_f \approx 1{,}7\,\%$)}\label{sec:r05}
\subsection{AI contra o direto (etapas 0, 2 e 2c)}\label{sec:aidireto}
\begin{table}[H]\centering\small
""" + T_R05 + r"""
\caption{$r = -0{,}5$ (Cho: $P_f = 1{,}71\,\%$, Monte Carlo de 50\,000 com o equilíbrio limite). O direto tem
""" + f"{d05['n']:,}".replace(",", "\\,") + r""" realizações. ``dif.'' é a diferença AI $-$ direto em desvios;
ESS$_F = (\sum wI)^2/\sum (wI)^2$; o ganho é quantas realizações do direto valem uma da AI, para o mesmo CoV. Sem
indicação, o FORM usa a semente 3 e a banda a semente 1; [PC] e [fermi] dizem onde rodou. As de Nyström usam o campo de
200 termos (variância ≈0,93), e a comparação com o direto é aproximada (Seção~\ref{sec:mesmavar}). As da fase 4 juntam
as sementes parciais (Seção~\ref{sec:fase4}).}\label{tab:r05}
\end{table}
\begin{figure}[H]\centering\includegraphics[width=0.85\textwidth]{figs/pf_r05.pdf}
\caption{$P_f$ e IC 95\,\% das AI do fermi. A faixa cinza é o IC do direto; a linha tracejada, o valor do Cho.}\end{figure}

\subsection{O teste linear (sem elementos finitos)}\label{sec:linear}
""" + (r"""O Cho2010 sorteia os campos como sempre, com a mesma semente, o mesmo deslocamento e o mesmo peso, mas troca o FS por
$\exp(\beta - \bm\alpha\cdot\bm\xi)$. A falha é $\bm\alpha\cdot\bm\xi > \beta$, com $P_f = \Phi(-\beta)$ exata. A direção
$\bm\alpha$ é a do ponto de projeto do FORM com $M = 20$. A AI ``M = 20'' desloca exatamente nessa direção. A AI ``M = 10'' usa
a direção do FORM com $M = 10$, a """ + num(math.degrees(math.acos(LINTAB[1])), 1) + r"""$^\circ$ da correta. Cada AI tem
$10^4$ realizações, e o direto, $10^6$.
\begin{table}[H]\centering\small
""" + tab_simples(r"$\beta$ & $P_f$ exata & direto ($10^6$) & AI $M = 20$: erro (desvios) & AI $M = 10$: erro (desvios)",
                  [f"{x['b']} & {x['ex']:.3e} & " + (f"{pct(x['d'])}\\,\\%" if x['d'] else "---") +
                   f" & {num(100 * (x['ai20'][0] / x['ex'] - 1), 1)}\\,\\% ({num(x['ai20'][2], 1)}) & "
                   f"{num(100 * (x['ai10'][0] / x['ex'] - 1), 1)}\\,\\% ({num(x['ai10'][2], 1)})" for x in LINTAB[0]],
                  "rrrrr") + r"""
\caption{O teste linear.}
\end{table}
\begin{itemize}
  \item \textbf{O código da AI está certo:} com o deslocamento na direção exata, o erro fica abaixo de 1\,\% até
    $P_f = 3\times10^{-5}$.
  \item \textbf{Com o deslocamento um pouco fora da direção, a AI sai baixa,} e mais baixa quando $P_f$ cai. O estimador
    não tem viés, mas a região que mais pesa fica na cauda do sorteio: poucas amostras, com pesos enormes. Com $N$ finito,
    a estimativa típica fica abaixo do valor certo, e o CoV calculado subestima o erro.
  \item Isso explica a AI pelo FORM um pouco abaixo do direto nos elementos finitos (Tabelas~\ref{tab:r05} e~\ref{tab:r07}): o ponto de projeto
    nunca está exatamente na direção ótima.
\end{itemize}
""" if LINTAB else r"""O resultado corrigido ainda não foi baixado.
""") + r"""

\subsection{A previsão teórica: $M = 10$ contra $M = 20$}
No modelo linear, deslocar para o ponto de projeto truncado multiplica a variância da AI por cerca de
$\exp(\beta_{10}^2 - \beta_{20}^2)$. As duas AI de cada $r$ usam a mesma semente (números aleatórios comuns).
\begin{table}[H]\centering\small
""" + tab_simples(r"caso & previsto: $\exp(\beta_{10}^2-\beta_{20}^2)$ & medido: $(\mathrm{CoV}_{10}/\mathrm{CoV}_{20})^2$",
                  [f"{n} & {num(pv, 2)} & {num(md, 2)}" for n, pv, md in PREV], "lrr") + r"""
\caption{A razão das variâncias das AI com o ponto de projeto de $M = 10$ e de $M = 20$.}
\end{table}
A ordem de grandeza confere: o $M = 20$ reduz a variância à metade. Com 1000 amostras, o próprio CoV tem cerca de 10\,\% de
erro, e as razões medidas são compatíveis com as previstas.

\section{Etapa 2b: as duas bases de Karhunen--Loève}\label{sec:bases}
A base analítica é a do retângulo que contém o talude, e a de Nyström é a do próprio talude (formulação,
Seção~\ref{F-sec:bases}).
\subsection{O FORM com o gradiente modal nas duas bases}\label{sec:nystrom}
\begin{table}[H]\centering\small
""" + T_FORM + r"""
\caption{$\beta$ do FORM do SRM com o gradiente modal na malha usual. Com diferenças centrais: 2,115 ($M = 10$, 123
análises) e 2,012 ($M = 20$, 243).}\label{tab:form}
\end{table}
\begin{figure}[H]\centering\includegraphics[width=0.6\textwidth]{figs/beta_M.pdf}\caption{$\beta(M)$: a base de
Nyström converge com menos termos.}\end{figure}
A base de Nyström chega mais perto do $\beta$ convergido com poucos termos ($M = 10$: 2,062 contra 2,151), e as duas
coincidem a partir de $M = 50$. É o mesmo efeito da placa curva do estudo do Ghanem \& Spanos (\path{docs/confiabilidade_campos/kl_form_placas}). O $\beta$ converge para
$\approx 2{,}00$, mas a $P_f$ do FORM sozinho ($\approx 2{,}3\,\%$) superestima a do direto em ≈30\,\%: é o erro de
linearizar uma superfície de falha curva. Por isso o FORM entra só como deslocamento da AI.

\subsection{Comparação das bases com a mesma variância}\label{sec:mesmavar}
A variância do campo nos pontos do talude depende da base e de $M$ (\path{runs/cho2010/nystrom/variancia_bases.py}):
\begin{table}[H]\centering\small
""" + tabela("$M$ & analítica (retângulo) & Nyström (talude)",
             ["10 & 0,60 & 0,65 \\\\", "20 & 0,72 & 0,76 \\\\", "50 & 0,84 & 0,86 \\\\", "100 & 0,90 & 0,90 \\\\",
              "200 & 0,935 & 0,931 \\\\", "400 & 0,96 & 0,94 \\\\", "2000 & 0,99 & 0,95 (teto) \\\\"], "rrr") + r"""
\caption{Variância média do campo nos pontos do talude (alvo 1). O teto de Nyström vem da quadratura de 2400 nós: a
interpolação entre os nós não reproduz o bico do kernel exponencial.}
\end{table}
Com poucos termos, Nyström retém mais variância; em $M \approx 200$ as duas se igualam (≈0,93). A AI analítica usa 2000
termos (0,99), e a de Nyström 200 (0,93): são modelos de campo diferentes. Mas a variância que falta está nos modos
altos, que quase não mudam o FS (a média ao longo da banda de ruptura os filtra): o $\beta$ do FORM vai de 2,001
($M = 100$) para 1,999 ($M = 200$). Por isso a $P_f$ do campo com 0,93 deve ficar perto da do campo completo; a
estimativa inicial, de 1,3 a 1,4\,\% pela variância, estava errada. A comparação justa das bases usa os dois campos
com a mesma variância:
\begin{table}[H]\centering\small
""" + tabela("Campo (variância ≈0,93) & $N$ & $P_f$ (\\%) & CoV (\\%)",
             destaca_colunas([f"{r} & " + (f"{a['n']} & {pct(a['p'])} & {num(100 * a['cov'], 1)}" if a else "--- & na fila & ---") +
                             " \\\\" for r, a in iguais], {3: "min"}, lambda i, c: c[1].strip().isdigit() and int(c[1]) >= 1000),
             "lrrr") + r"""
\caption{AI no PC, semente 3, a mesma das AI analíticas da Tabela 1 (os sorteios são os mesmos). A referência de
24\,000 é do campo completo (2000 termos) e só vale aproximadamente aqui.}
\end{table}

\subsection{O ponto de projeto nas duas bases}\label{sec:pp}
Com o gradiente modal (formulação, Seção~\ref{F-sec:form}), o FORM converge em 3 análises nas duas bases
(Figura~\ref{fig:pp}).
\begin{itemize}
  \item \emph{Analítica} ($M = 100$): $\beta = """ + num(FM[100][0]["b"], 3) + r"""$, com """ + num(UKL["phi"], 2) + r"""
    de $\alpha^2$ em $\varphi'$; """ + UKL["txt"] + r""". Os 10 primeiros modos de cada campo somam
    """ + num(UKL["dez"], 2) + r""".
  \item \emph{Nyström} ($M = 200$): $\beta = """ + num(FM[200][1]["b"], 3) + r"""$, com """ + num(UNY["phi"], 2) + r"""
    em $\varphi'$; """ + UNY["txt"] + r""". Os 10 primeiros somam """ + num(UNY["dez"], 2) + r""".
\end{itemize}
Os coeficientes são outros, mas o ponto é o mesmo. Levado para as coordenadas de Nyström por $T$, o ponto analítico
tem $|\bm\xi| = """ + num(PP["beta_kl_em_nys"], 3) + r"""$ e fica a """ + num(PP["angulo"], 1) + r"""$^\circ$ do de
Nyström. O campo mais provável de ruptura é uma camada fraca horizontal sob a face do talude:
\begin{itemize}
  \item $\varphi'$ mínimo de 23,4$^\circ$ em $y \approx 6{,}6$\,m, 1,6\,m acima do pé;
  \item $c'$ quase uniforme, de 8,6 a 9,6\,kPa;
  \item os perfis de $\varphi'$ das duas bases diferem em menos de 0,1$^\circ$.
\end{itemize}
\begin{figure}[H]\centering\includegraphics[width=\textwidth]{figs/ponto_projeto.pdf}
\caption{$c'$ e $\varphi'$ no ponto de projeto do FORM modal: KL analítica ($M = 100$, em cima) e Nyström ($M = 200$,
embaixo), com as mesmas escalas de cor.}\label{fig:pp}\end{figure}

\section{Etapa 2c: o gradiente modal}\label{sec:modal}
O FORM com diferenças centrais custa $4M + 1$ análises por iteração: 123 com $M = 10$ e 243 com $M = 20$. O gradiente
modal sai da própria análise de colapso, pela análise limite e pelo teorema do envelope (formulação,
Seção~\ref{F-sec:gradmodal}): uma análise por iteração, qualquer que seja $M$.
\subsection{A validação do gradiente: SRM e GIM}\label{sec:validacao}
\begin{itemize}
  \item \textbf{Validação do gradiente modal} nos pontos de projeto do fermi: ele fica a 9,6$^\circ$ ($M = 10$) e
    8,9$^\circ$ ($M = 20$) do gradiente das diferenças centrais, e o termo dominante ($\varphi'_0$) bate a menos de 4\,\%.
""" + GIM_TXT + r"""
\end{itemize}
\begin{table}[H]\centering\small
""" + T_GIM + r"""
\caption{$\beta$ do FORM (corrigido: a distância da origem ao plano tangente no último ponto) com cada gradiente, na malha
usual. GIM: $g = \Gamma - 1$; SRM: $g = \mathrm{FS} - 1$. A região de falha é a mesma nos dois, então o $\beta$ também
deveria ser.}\label{tab:gim}
\end{table}

\subsection{A AI pelo ponto de projeto modal}\label{sec:aimodal}
O deslocamento é o ponto de projeto do FORM modal (formulação, Seção~\ref{F-sec:aiform}). As AI estão na
Tabela~\ref{tab:r05}: $M = 100$ e 200 com a semente 3 [PC], e as sementes 6 e 7 [fermi]; as de Nyström, na
Seção~\ref{sec:mesmavar}. O ponto modal de $M = 100$ dá o mesmo resultado que o das diferenças centrais com muito menos
análises, e o viés está descartado na Seção~\ref{sec:vies}.

\subsection{A fase 4, parcial}\label{sec:fase4}
A fase 4 roda 10 sementes do FORM modal analítico ($M = 100$) e 10 do de Nyström ($M = 200$), as mesmas nas duas bases,
para ter a distribuição do estimador. Ela foi cancelada em 17/11 de manhã, para dar o fermi ao direto de $r = -0{,}7$, com
≈65 realizações por semente. As sementes com o mesmo deslocamento são juntadas numa AI só: sementes diferentes dão
realizações diferentes, e o estimador é a média de $I w$ de todas.
\begin{table}[H]\centering\small
""" + T_F4 + r"""
\caption{A fase 4 até o cancelamento. A dif. é contra o direto de $r = -0{,}5$; nas de Nyström, ela é aproximada
(Seção~\ref{sec:mesmavar}).}\label{tab:fase4}
\end{table}


\section{A AI pelo FORM tem viés?}\label{sec:vies}
\subsection{A semente 3 e as repetições independentes}
As primeiras AI pelo FORM usavam todas a semente 3 (números aleatórios comuns): com pontos de projeto quase iguais, elas dão
quase a mesma estimativa e não são independentes. As repetições com as sementes 4 e 5 (FORM $M = 20$) e 6 e 7 (FORM
modal, $M = 100$ e 200) separam o acaso do viés:
\begin{table}[H]\centering\small
""" + tabela("Caso & sementes & $P_f$ de cada (\\%) & média (\\%) & direto (\\%) & dif. (desvios)",
             [f"$r = -0{{,}}5$, FORM $M = 20$ & {M05['k']} & " + ", ".join(pct(a['p']) for a in SEM05 if a) +
              f" & {pct(M05['p'])} & {pct(d05['p'])} & {num((M05['p'] - d05['p']) / math.hypot(M05['sd'], d05['p'] * d05['cov']), 1)} \\\\",
              f"$r = -0{{,}}5$, FORM $M = 20$, 100 e 200 modal & {M05T['k']} & " + ", ".join(pct(a['p']) for a in SEM05T if a) +
              f" & {pct(M05T['p'])} & {pct(d05['p'])} & {num(dif_media(M05T, d05), 1)} \\\\",
              f"$r = -0{{,}}7$, FORM $M = 20$ & {M07['k']} & " + ", ".join(pct(a['p']) for a in SEM07 if a) +
              f" & {pct(M07['p'])} & {pct(d07['p'])} & {num((M07['p'] - d07['p']) / math.hypot(M07['sd'], d07['p'] * d07['cov']), 1)} \\\\"],
             "lrrrrr") + r"""
\caption{AI pelo FORM com sementes independentes, 1000 realizações cada. A dif. usa o desvio da média,
$\sqrt{\sum s_i^2}/k$, e o do direto.}
\end{table}

\subsection{Um modo de ruptura só?}\label{sec:modo}
Se houvesse um segundo mecanismo, a AI pelo FORM quase não o sortearia, e a $P_f$ sairia baixa, com um CoV pequeno enganoso.
Para verificar, os $\bm\xi$ das realizações do direto foram sorteados de novo, sem elementos finitos, com o modo
\texttt{linear} do Cho2010. Em cada realização calculou-se $t = \bm\alpha\cdot\bm\xi$, a coordenada na direção do ponto de
projeto do FORM ($M = 20$).
\begin{itemize}
  \item \textbf{Direção:} nenhuma falha tem $t < 1$. O menor $t$ entre as falhas é """ + num(diag['r05']['tmin'], 2) + \
    r""" com $r = -0{,}5$ (""" + f"{diag['r05']['nf']}" + r""" falhas) e """ + num(diag['r07']['tmin'], 2) + \
    r""" com $r = -0{,}7$ (""" + f"{diag['r07']['nf']}" + r""" falhas).
  \item \textbf{Geometria:} a banda de ruptura das falhas é a mesma no direto e nas AI ponderadas: profundidade de
    """ + num(diag['r05']['prof'], 1) + r"""\,m e centro em $x \approx """ + num(diag['r05']['xc'], 1) + r"""$\,m.
  \item \textbf{Conclusão:} há um modo só. Os desvios negativos da AI pelo FORM (de $-0{,}6$ a $-1{,}3$) são compatíveis com o
    acaso. As duas AI de cada $r$ usam a mesma semente e não são independentes. As repetições com sementes novas
    (\texttt{\_s4}, \texttt{\_s5}) vão confirmar.
\end{itemize}
\begin{figure}[H]\centering\includegraphics[width=\textwidth]{figs/modo.pdf}
\caption{$t = \bm\alpha\cdot\bm\xi$ nas realizações do direto. As falhas (vermelho) ficam todas à direita, perto ou além de
$\beta$.}\end{figure}

\subsection{O teste com o próprio direto (KL e Nyström)}\label{sec:diag}
Com a semente 3, todas as AI pelo FORM ficaram abaixo do direto (1,53 a 1,57\,\% contra 1,70\,\%), e a AI modal, a mais
eficiente, parecia convergir para o valor errado. Uma AI com deslocamento fixo não tem viés, mas pode se comportar como se
tivesse. Se a densidade de amostragem quase não visita uma parte da região de falha (outro modo de ruptura, longe do
ponto de projeto), a estimativa sai quase sempre baixa, e o CoV estimado fica menor que o verdadeiro. O direto permite
testar isso sem nenhuma análise nova (\path{Projects/Cho2010/diagnostico_ai.py}):
\begin{enumerate}
  \item \textbf{Os $\bm\xi$ de cada realização} saem de uma cópia, em Python, do gerador do \texttt{RandomField}:
    xoshiro256** e Box--Muller, com uma semente por realização e campo. A cópia confere com o peso de cada amostra de cada
    AI (erro $\le """ + cient(max(a["erro_peso"] for a in DG["ai"])) + r"""$). Junto com a KL analítica em Python, ela
    também confere com o campo que o programa grava (diferença $\le 5\times10^{-6}$, a precisão do CSV).
  \item \textbf{Modo de ruptura fora da AI.} As """ + str(DD5["falhas"]) + r""" falhas do direto, projetadas na direção
    do ponto de projeto ($t = \bm\xi\cdot\bm\xi^*/\beta$), ficam todas onde a AI amostra (Figura~\ref{fig:vies_t}). Na direção
    modal, o menor $t$ é """ + num(AM["t_min"], 2) + r""", e só """ + str(AM["abaixo_15"]) + r""" falhas ficam abaixo de
    1,5. A variância verdadeira de cada AI sai do direto pela identidade $E_h[I w^2] = E_\varphi[I w]$, com as falhas do
    direto, e é a que a AI estima (``CoV verd.'' na Tabela~\ref{tab:vies05}): """ + num(100 * AM["cov_verdadeiro"], 1) +
    r"""\,\% contra """ + num(100 * AM["sd"] / AM["p"], 1) + r"""\,\% no modal $M = 100$. Só a direção errada da Tri(0)
    subestima o próprio CoV (""" + num(100 * AT0["sd"] / AT0["p"], 1) + r"""\,\% estimado contra """ +
    num(100 * AT0["cov_verdadeiro"], 1) + r"""\,\% verdadeiro): o $-2{,}2$ desvio da Tabela~\ref{tab:r05} vinha daí.
  \item \textbf{Base de Nyström.} Com $T$ (formulação, Seção~\ref{F-sec:mudanca}), as realizações do direto vão para as coordenadas de
    Nyström, e o mesmo teste vale para as AI de Nyström: $t \ge """ + num(min(AN20["t_min"], AN200["t_min"]), 2) + r"""$, e
    CoV verdadeiro de """ + num(100 * min(AN20["cov_verdadeiro"], AN200["cov_verdadeiro"]), 1) + r""" a """ +
    num(100 * max(AN20["cov_verdadeiro"], AN200["cov_verdadeiro"]), 1) + r"""\,\%, contra """ +
    num(100 * min(AN20["sd"] / AN20["p"], AN200["sd"] / AN200["p"]), 1) + r""" a """ +
    num(100 * max(AN20["sd"] / AN20["p"], AN200["sd"] / AN200["p"]), 1) + r"""\,\% estimados. Os FS são os do campo
    analítico, uma aproximação: o campo de Nyström é a projeção dele. A regressão de $\ln\mathrm{FS}$ explica o mesmo nas
    duas bases ($R^2 = """ + num(DD5["r2_kl"], 3) + r"""$ e $""" + num(DD5["r2_nys"], 3) + r"""$).
  \item \textbf{A sorte de cada amostra.} A variável de controle mede a sorte de cada amostra no que a regressão explica.
    A amostra da semente 3 é azarada: o controle fica em """ + num(100 * AM["sorte"], 0) + r"""\,\% da sua média exata na
    AI modal, """ + num(100 * A20["sorte"], 0) + r"""\,\% no FORM $M = 20$ e """ + num(100 * AN200["sorte"], 0) + r"""\,\%
    no Nyström $M = 200$. Corrigida, a AI modal da semente 3 vai de """ + pct(AM["p"]) + r""" para """ + pct(AM["p_cv"]) +
    r"""\,\%. O direto corrigido fica em """ + pct(DD5["cv"]["p"], 3) + r""" $\pm$ """ + pct(DD5["cv"]["sd"], 3) + r"""\,\%
    (CoV de """ + num(100 * DD5["cv"]["sd"] / DD5["cv"]["p"], 1) + r"""\,\%), o equivalente a """ +
    f"{DD5['n_equiv']:,.0f}".replace(",", "\\,") + r""" realizações.
\end{enumerate}
\begin{figure}[H]\centering\includegraphics[width=\textwidth]{figs/vies_t.pdf}
\caption{A chance de falha do direto em cada faixa de $t = \bm\xi\cdot\bm\xi^*/\beta$ (barras) e as densidades de $t$ no
direto e na AI (curvas, em escala). A AI amostra toda a região onde o direto falha.}\label{fig:vies_t}\end{figure}
\begin{table}[H]\centering\small
""" + T_VIES05 + r"""
\caption{$r = -0{,}5$: as AI pelo ponto de projeto. ``CoV verd.'' é o CoV verdadeiro de cada AI, calculado pelas falhas
do direto; $t_{\min}$ é o menor $t$ entre as falhas do direto; ``corr.'' é a estimativa com a variável de controle.
Em Nyström, o direto vai para a base de Nyström pela mudança de base.}\label{tab:vies05}
\end{table}
\begin{table}[H]\centering\small
""" + T_SEM + r"""
\caption{As sementes independentes (uma AI por semente; a 3 é a do PC). A dif. é a diferença para o direto, em desvios;
a corrigida usa o direto corrigido. A última coluna compara a dispersão real entre as sementes com o desvio que cada AI
estima.}\label{tab:sementes}
\end{table}
\begin{figure}[H]\centering\includegraphics[width=\textwidth]{figs/vies_sementes.pdf}
\caption{Uma AI por semente, com IC 95\,\%: vazado, a estimativa da AI; cheio, a corrigida. As faixas são o IC 95\,\% do
direto (claro) e do direto corrigido (escuro).}\label{fig:vies_sementes}\end{figure}
\textbf{$r = -0{,}7$} (direto de """ + f"{DD7['n']:,}".replace(",", "\\,") + r""", """ + str(DD7["falhas"]) + r"""
falhas): $t \ge """ + num(A07[0]["t_min"], 2) + r"""$ na direção do FORM $M = 20$ (""" + num(A0710["t_min"], 2) + r"""
no $M = 10$). O CoV verdadeiro é """ + num(100 * A07[0]["cov_verdadeiro"], 1) + r"""\,\%, e as AI estimam de """ +
num(100 * min(a["sd"] / a["p"] for a in A07), 1) + r""" a """ + num(100 * max(a["sd"] / a["p"] for a in A07), 1) + r"""\,\%.
As sementes 3, 4 e 5, corrigidas, dão """ + ", ".join(pct(a["p_cv"], 3) for a in A07) + r"""\,\%, contra """ +
pct(DD7["cv"]["p"], 3) + r""" $\pm$ """ + pct(DD7["cv"]["sd"], 3) + r"""\,\% do direto corrigido ($""" +
(num(S07["difc"], 1) if S07 else "---") + r"""$ desvio).
\begin{table}[H]\centering\small
""" + T_VIES07 + r"""
\caption{$r = -0{,}7$: as AI pelo ponto de projeto, como na Tabela~\ref{tab:vies05}.}\label{tab:vies07}
\end{table}
\textbf{Conclusão.}
\begin{itemize}
  \item Nenhum modo de ruptura fica fora da AI, e a variância verdadeira é a estimada.
  \item As sementes independentes acertam o direto, nas duas bases. A semente 3 é uma tiragem baixa, comum a todas as
    AI que a usam.
  \item As repetições da fase 4 no fermi fecham a questão: 10 sementes do FORM modal analítico e 10 do de Nyström, com as
    mesmas sementes nas duas bases. Elas dão a distribuição do estimador.
  \item A variável de controle é um ganho à parte, sem custo: de """ + num(min(GANHOS), 1) + r""" a """ +
    num(max(GANHOS), 1) + r""" vezes menos realizações nas AI e """ + num(DD5["n_equiv"] / DD5["n"], 1) + r""" no direto.
\end{itemize}

\section{A variável de controle}\label{sec:cv}
O método está na formulação (Seção~\ref{F-sec:controle}). Com o substituto ajustado no direto inteiro, o ganho vai de
""" + num(min(GANHOS), 1) + " a " + num(max(GANHOS), 1) + r""" vezes menos realizações nas AI de $r = -0{,}5$, e é de
""" + num(DD5["n_equiv"] / DD5["n"], 1) + r""" no direto (Seção~\ref{sec:diag} e Tabela~\ref{tab:vies05}).
\subsection{A variável de controle com o que se tem na prática}\label{sec:cvlocal}
O diagnóstico da Seção~\ref{sec:diag} ajusta o substituto (formulação, Seção~\ref{F-sec:controle}) no direto inteiro, que
não existe num problema real: se existisse, a AI não faria falta. O teste com o que se tem na prática
(\path{Projects/Cho2010/teste_controle.py}, no PC, sem análises novas):
\begin{enumerate}
  \item \textbf{O direto em lotes.} As """ + f"{TC['n_direto']:,}".replace(",", "\\,") + r""" realizações do direto de
    $r = -0{,}5$ em """ + str(len(LOTES["P"]["por_lote"])) + r""" lotes de """ + str(TC["lote"]) + r""" (o tamanho de uma
    AI). Em cada lote, a $P$ do lote e a $P_{cv}$ com o substituto ajustado só no próprio lote (5 dobras, $K$ termos por
    campo) ou fora dele (o ideal, $K = 200$). O desvio entre os lotes é o verdadeiro (Tabela~\ref{tab:cvlotes}).
  \item \textbf{As AI}, com o substituto ajustado na própria amostra de 1000 (5 dobras, $K = 50$) e, para comparar, no
    direto inteiro (Tabela~\ref{tab:cvai}).
\end{enumerate}
\begin{table}[H]\centering\small
""" + T_CV_LOTES + r"""
\caption{O direto de $r = -0{,}5$ em lotes de """ + str(TC["lote"]) + r""". O ganho é a razão das variâncias entre os
lotes. A $P_f$ do direto inteiro é """ + pct(TC["P_direto"], 3) + r"""\,\%.}\label{tab:cvlotes}
\end{table}
\begin{table}[H]\centering\small
""" + T_CV_AI + r"""
\caption{As AI de $r = -0{,}5$ com a variável de controle: o substituto ajustado na própria amostra (5 dobras, $K = 50$
por campo) e no direto inteiro (o ideal). A semente 3 é a mesma em várias linhas.}\label{tab:cvai}
\end{table}
\begin{itemize}
  \item \textbf{Sem viés:} as médias dos lotes ficam todas perto da $P_f$ do direto inteiro.
  \item \textbf{O ganho realista depende de $K$:} """ + ", ".join(num(LOTES[f"realista K={K}"]["ganho"], 1) + f" com $K = {K}$" for K in (10, 20, 50)) + r""",
    contra """ + num(LOTES["ideal K=200"]["ganho"], 1) + r""" do ideal. Com 800 realizações por ajuste, os $2K + 1 = 101$
    coeficientes de $K = 50$ ainda cabem.
  \item \textbf{O desvio estimado da $P_{cv}$ sai menor que o verdadeiro} (de """ +
    num(100 * (1 - max(LOTES[k]["desvio_estimado"] / LOTES[k]["desvio"] for k in LOTES if k != "P")), 0) + r""" a """ +
    num(100 * (1 - min(LOTES[k]["desvio_estimado"] / LOTES[k]["desvio"] for k in LOTES if k != "P")), 0) + r"""\,\%), enquanto o da
    $P$ bate. Com """ + str(len(LOTES["P"]["por_lote"])) + r""" lotes, o desvio verdadeiro tem uma incerteza de ≈15\,\%:
    a confirmar com as 10 sementes da fase 4. Até lá, o desvio da $P_{cv}$ é otimista.
  \item \textbf{Nas AI,} o substituto da própria amostra dá ganho de """ + num(min(GANHO_K50), 1) + r""" a """ +
    num(max(GANHO_K50), 1) + r""", contra """ + num(min(GANHO_ID), 1) + r""" a """ + num(max(GANHO_ID), 1) + r""" com o
    do direto. A direção errada da Tri(0) não se corrige com o substituto da própria amostra: ele só vê a região que a AI
    amostra.
\end{itemize}

\section{Etapa 3: $P_f$ pequena, $r = -0{,}7$ ($P_f \approx 0{,}4\,\%$)}\label{sec:r07}
\begin{table}[H]\centering\small
""" + T_R07 + r"""
\caption{$r = -0{,}7$. Cho (2010): $P_f \approx 0{,}4\,\%$ (Fig.~9). O direto tem """ +
f"{d07['n']:,}".replace(",", "\\,") + r""" realizações e """ + f"{d07['k']}" + r""" falhas (o IC dele é o exato, de
Clopper--Pearson) e continua rodando até 110\,000. O FORM convergiu em """ + f"{f07[0]['it']}" + r""" ($M = 10$) e """ +
f"{f07[1]['it']}" + r""" ($M = 20$) iterações.}\label{tab:r07}
\end{table}
\begin{figure}[H]\centering\includegraphics[width=0.85\textwidth]{figs/pf_r07.pdf}
\caption{$P_f$ e IC 95\,\% com $r = -0{,}7$.}\end{figure}
O teste de viés com $r = -0{,}7$ está na Seção~\ref{sec:diag} (Tabela~\ref{tab:vies07}).


\section{O piloto e o custo}\label{sec:custo}
\subsection{Tamanho do piloto (AI pela banda)}
O deslocamento da AI pela banda é $\bm\mu = s_c\,\bm d_c + s_\varphi\,\bm d_\varphi$, com $s$ = a média de $\bm\xi\cdot\bm d$ no
conjunto elite do piloto (os 10\,\% de menor FS, ou as falhas se forem mais). Para ver quanto $s$ depende do tamanho do
piloto, o direto foi dividido em sub-pilotos disjuntos de 100, 250, 500 e 1000 realizações. As projeções foram refeitas
sem elementos finitos, pelo modo \texttt{linear} com \texttt{banda=}. O sub-piloto 0--999 reproduz o deslocamento que o
programa usou no fermi: $s = (""" + num(SUB['r05'][3][4][0], 3) + "; " + num(SUB['r05'][3][4][1], 3) + r""")$.
\begin{table}[H]\centering\small
""" + tab_simples(r"caso & piloto & sub-pilotos & $s_c$: média (dp) & $s_\varphi$: média (dp) & $|\bm\mu|^2$: média (dp)",
                  [f"{nome} & {m} & {k} & {num(mu[0], 2)} ({num(sd[0], 2)}) & {num(mu[1], 2)} ({num(sd[1], 2)}) & "
                   f"{num(mu[2], 2)} ({num(sd[2], 2)})" for nome, rot in (("r05", "$r = -0{,}5$"), ("r07", "$r = -0{,}7$"))
                   for (m, k, mu, sd, _) in SUB[nome] for nome in [rot]], "lrrrrr") + r"""
\caption{O deslocamento da banda com pilotos de tamanhos diferentes.}
\end{table}
A média de $s$ não depende do tamanho do piloto; só a dispersão cai. Um erro $\bm\delta$ no deslocamento multiplica a
variância da AI por cerca de $e^{|\bm\delta|^2}$. Com piloto de 100, $|\bm\delta|^2 \approx """ + \
    num(SUB['r05'][0][3][0] ** 2 + SUB['r05'][0][3][1] ** 2, 2) + r"""$ (variância $\times""" + \
    num(math.exp(SUB['r05'][0][3][0] ** 2 + SUB['r05'][0][3][1] ** 2), 2) + r"""$); com 250, $""" + \
    num(SUB['r05'][1][3][0] ** 2 + SUB['r05'][1][3][1] ** 2, 2) + r"""$ ($\times""" + \
    num(math.exp(SUB['r05'][1][3][0] ** 2 + SUB['r05'][1][3][1] ** 2), 2) + r"""$). Valores com $r = -0{,}5$. Um piloto de 250 basta.

\subsection{Custo total para CoV$(P_f) = 10\,\%$}
O número de análises para chegar a CoV de 10\,\% é $N\,(\mathrm{CoV}/0{,}10)^2$. A ele soma-se o custo fixo de cada método:
o piloto de 1000 na banda e as análises do FORM.
\begin{table}[H]\centering\small
""" + tab_simples(r"Método & $r = -0{,}5$: amostragem & + fixo & total & $r = -0{,}7$: amostragem & + fixo & total",
                  [f"{a[0]} & {a[1]:.0f} & {a[2]} & {a[3]:.0f} & {b[1]:.0f} & {b[2]} & {b[3]:.0f}" for a, b in zip(CUSTO05, CUSTO07)],
                  "lrrrrrr") + r"""
\caption{Análises de elementos finitos para CoV$(P_f) = 10\,\%$.}
\end{table}
Mesmo pagando o FORM, a AI pelo ponto de projeto custa de """ + num(CUSTO05[0][3] / CUSTO05[2][3], 0) + " a " + \
    num(CUSTO07[0][3] / CUSTO07[3][3], 0) + r""" vezes menos que o direto. A vantagem cresce quando a $P_f$ cai. O piloto de
1000 da banda pesa muito. Com um piloto de 250, que basta (tabela anterior), a banda iria a """ + \
    f"{CUSTO05[1][1] + 250:.0f}" + " e " + f"{CUSTO07[1][1] + 250:.0f}" + r""" análises, ainda 1,5 a 3 vezes o FORM.

\subsection{O piloto em malha grossa}
O elite do piloto é escolhido por ordem de FS, e a ordem das realizações quase não muda com a malha. Com o piloto na
Tri(0), 300 vezes mais barata, a AI pela banda acertou o direto, com o mesmo CoV do piloto na malha usual; Tri(1) e
Tri(2) também estão na Tabela~\ref{tab:r05}. O FORM, ao contrário, não pode rodar na malha grossa: na Tri(0) (76
equações, FS 1,47 contra 1,20), $\beta = 4{,}08$ e a direção fica a 31,7$^\circ$ da certa. O erro da malha desloca o
ponto de projeto, e não só o FS. Com a direção da Tri(0) e a distância da malha usual ($\beta = 2{,}32$), a AI funciona,
mas pior que o FORM na malha usual.


\section{Etapa 4: a referência da chuva}\label{sec:chuva}
O caso é o f50 de Rahardjo a 31,05\,mm/h, com o FS só em 24\,h.
\begin{table}[H]\centering\small
""" + T_CH + r"""
\caption{Chuva: o Monte Carlo direto de """ + f"{dch['n']:,}".replace(",", "\\,") + r""" é a referência para validar o
FORM e a AI na chuva (etapa 4).}\label{tab:chuva}
\end{table}
A AI pela banda concorda com o direto (""" + num(desvios(ach[0][1], dch), 1) + r""" desvio), com ganho de """ +
num(ganho(ach[0][1], dch), 1) + r"""$\times$. Falta o FORM na chuva, que precisa de código novo (as opções \texttt{M},
\texttt{xi} e \texttt{desloc\_xi} no \texttt{montecarlo.cpp}).

\section{O que falta}
Os próximos passos estão no relatório de situação (\path{Projects/PlanoDeTrabalho/situacao/situacao.pdf}).
\end{document}
""")
open(f"{AQUI}/relatorio_etapas.tex", "w").write(tex)
print("ok: relatorio_etapas.tex e figs/")
print(f"r05: direto {pct(d05['p'])}%;", "; ".join(f"{re.sub(r'[$\\\\]', '', c)} {pct(a['p'])}% ganho {ganho(a, d05):.1f}" for c, a in a05))
print(f"r07: direto {pct(d07['p'])}%;", "; ".join(f"{re.sub(r'[$\\\\]', '', c)} {pct(a['p'])}% ganho {ganho(a, d07):.1f}" for c, a in a07))
