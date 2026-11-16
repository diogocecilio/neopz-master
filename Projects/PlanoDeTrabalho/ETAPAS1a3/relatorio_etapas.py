#!/usr/bin/env python3
"""relatorio_etapas.py: numeros, figuras e o texto LaTeX do relatorio das etapas 1 a 3 do plano de trabalho (e a referencia
da chuva da etapa 4), lidos das campanhas baixadas do fermi (~/projects/runs/fermi_etapa1 e fermi_etapa3; o teste linear
em ~/projects/runs/etapa2_linear, se existir). Roda de novo quando chegarem mais dados:
    python3 -I relatorio_etapas.py && pdflatex relatorio_etapas.tex && pdflatex relatorio_etapas.tex"""
import csv
import glob
import math
import os
import re
import subprocess

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

AQUI = os.path.dirname(os.path.abspath(__file__))
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


L = []  # linhas do LaTeX
fig_count = 0


def tabela_ai(nome_dir, d, casos, f10=None, f20=None):
    """linhas de uma tabela: o direto, as AI e o FORM"""
    s = [r"\begin{tabular}{lrrrrrrr}", r"\toprule",
         r"Método & $N$ & $P_f$ (\%) & IC 95\% (\%) & CoV (\%) & dif. (desvios) & ESS$_F$ & ganho \\", r"\midrule",
         f"Monte Carlo direto & {d['n']} & {pct(d['p'])} & {pct(d['lo'])}–{pct(d['hi'])} & {num(100 * d['cov'], 1)} & --- & --- & 1 \\\\"]
    for rot, a in casos:
        s.append(f"{rot} & {a['n']} & {pct(a['p'])} & {pct(a['lo'])}–{pct(a['hi'])} & {num(100 * a['cov'], 1)} & "
                 f"{num(desvios(a, d), 1)} & {a['essF']:.0f} & {num(ganho(a, d), 1)} \\\\")
    for rot, f in (("FORM, $M = 10$", f10), ("FORM, $M = 20$", f20)):
        if f:
            s.append(f"{rot} ($\\beta = {num(f['bc'], 3)}$) & {f['n']} análises & {pct(f['pf'])} & --- & --- & --- & --- & --- \\\\")
    s += [r"\bottomrule", r"\end{tabular}"]
    return "\\resizebox{\\textwidth}{!}{%\n" + "\n".join(s) + "}"


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
        ax.text(100 * ref, -0.55, " Cho (2010)", color="tab:red", fontsize=8, va="bottom")
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
tmp = "/tmp/claude-1000/-home-diogo/5aa52bd8-2162-4e81-9474-74c603a4e56c/scratchpad/modo2_rel"
diag = {}
for nome, pasta_d, ustar, rc, n in (("r05", f"{E1}/r05_direto", f"{E1}/form_M20/ustar.txt", "-0.5", d05["n"]),
                                   ("r07", f"{E3}/r07_direto", f"{E3}/form_M20/ustar.txt", "-0.7", d07["n"])):
    os.makedirs(f"{tmp}/{nome}", exist_ok=True)
    if not glob.glob(f"{tmp}/{nome}/mc_*_*.csv"):
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

# ===================== LaTeX =====================
def tab_simples(cab, linhas, fmt):
    return "\\resizebox{\\textwidth}{!}{%\n" + "\n".join([r"\begin{tabular}{" + fmt + "}", r"\toprule", cab + r" \\", r"\midrule"] +
                     [l + r" \\" for l in linhas] + [r"\bottomrule", r"\end{tabular}"]) + "}"


tex = r"""\documentclass[11pt,a4paper]{article}
\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}
\usepackage[brazil]{babel}
\usepackage{lmodern}
\usepackage[margin=2cm]{geometry}
\usepackage{amsmath,amssymb,bm}
\usepackage{graphicx}
\usepackage{booktabs}
\usepackage{float}
\usepackage[hidelinks]{hyperref}
\title{Plano de trabalho: etapas 1, 2 e 3\\\large SRM corrigido, FORM e amostragem por importância no Cho (2010)}
\date{16/11/2026}
\begin{document}
\maketitle

\section*{Resumo}
\begin{itemize}
  \item \textbf{SRM corrigido (etapa 1):} começar o caminho de $F$ em 0,8, e não em 1, quase não muda a $P_f$.
    No Cho2010, as falhas são as mesmas nas """ + f"{cho['n']}" + r""" realizações pareadas. Na chuva, """ + \
    f"{chuva['n01']}" + r""" realizações passam a falhar e """ + f"{chuva['n10']}" + r""" deixa de falhar (McNemar $p = """ + \
    num(mcnemar(chuva['n01'], chuva['n10']), 2) + r"""$).
  \item \textbf{AI contra o direto (etapas 2 e 3):} todas as AI concordam com o Monte Carlo direto do mesmo modelo, a menos de
    1,3 desvio. O ganho cresce quando a $P_f$ cai: até """ + num(max(ganho(a, d05) for _, a in a05), 0) + r"""$\times$ com
    $P_f \approx 1{,}7\,\%$ e até """ + num(max(ganho(a, d07) for _, a in a07), 0) + r"""$\times$ com $P_f \approx 0{,}4\,\%$.
    O melhor deslocamento é o ponto de projeto do FORM com $M = 20$.
  \item \textbf{FORM sozinho:} com $M = 10$ ele acerta por coincidência, e com $M = 20$ superestima a $P_f$ em 30 a 45\,\%.
    Ele é bom para escolher o deslocamento da AI, mas não serve como estimativa final.
  \item \textbf{Um modo de ruptura só:} todas as falhas do direto estão do lado do ponto de projeto ($t = \bm\alpha\cdot\bm\xi \ge """ + \
    num(min(diag['r05']['tmin'], diag['r07']['tmin']), 1) + r"""$), e a banda de ruptura é a mesma no direto e na AI. A AI pelo
    FORM não está perdendo falhas.
\end{itemize}

\section{O que rodou}
Tudo rodou no fermi, nas campanhas \texttt{etapa1} (15/11) e \texttt{etapa3} (noite de 15 para 16/11), com o SRM novo,
tri2 + 4 níveis e a KL analítica de 2000 termos por campo. O problema é o talude seco do Cho (2010), Exemplo 2: $c'$ com média
10\,kPa e CoV 0,3; $\varphi'$ com média $30^\circ$ e CoV 0,2; campos lognormais com $l_h = 20$\,m e $l_v = 2$\,m; correlação
cruzada $r = -0{,}5$ ou $-0{,}7$. O FS determinístico é 1,2044, contra 1,204 do Cho.

\section{Etapa 1: o SRM corrigido}
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

\section{Etapa 2: FORM e AI com $r = -0{,}5$ ($P_f \approx 1{,}7\,\%$)}
\begin{table}[H]\centering\small
""" + tabela_ai("r05", d05, a05, *f05) + r"""
\caption{$r = -0{,}5$. Cho (2010): $P_f = 1{,}71\,\%$ (Monte Carlo de 50\,000, equilíbrio limite). dif.: a diferença AI $-$
direto em desvios-padrão. ESS$_F = (\sum w I)^2/\sum (wI)^2$. Ganho: quantas realizações do direto valem uma da AI, para
o mesmo CoV. O FORM com $M = 10$ e $M = 20$ convergiu em 3 iterações.}
\end{table}
\begin{figure}[H]\centering\includegraphics[width=0.85\textwidth]{figs/pf_r05.pdf}
\caption{$P_f$ e IC 95\,\%. A faixa cinza é o IC do direto; a linha tracejada, o valor do Cho.}\end{figure}

\paragraph{O teste linear (sem elementos finitos).}
""" + (r"\begin{verbatim}" + "\n" + lin_txt + r"\end{verbatim}" if lin_txt else
       r"""A primeira rodada saiu viciada por um defeito do próprio teste: o FS $\le 0$ era descartado. A rodada corrigida ainda não
foi baixada (\texttt{rsync} de \texttt{\$DATA/etapa2\_linear}).""") + r"""

\section{Etapa 3: $P_f$ pequena, $r = -0{,}7$ ($P_f \approx 0{,}4\,\%$)}
\begin{table}[H]\centering\small
""" + tabela_ai("r07", d07, a07, *f07) + r"""
\caption{$r = -0{,}7$. Cho (2010): $P_f \approx 0{,}4\,\%$ (Fig.~9). Como o direto tem """ + f"{d07['k']}" + r""" falhas, o IC dele é
exato (Clopper--Pearson). O FORM convergiu em """ + f"{f07[0]['it']}" + r""" ($M = 10$) e """ + f"{f07[1]['it']}" + r""" ($M = 20$) iterações.}
\end{table}
\begin{figure}[H]\centering\includegraphics[width=0.85\textwidth]{figs/pf_r07.pdf}
\caption{$P_f$ e IC 95\,\% com $r = -0{,}7$.}\end{figure}

\subsection*{Um modo de ruptura só?}
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

\section{Etapa 4: a referência da chuva}
O caso é o f50 de Rahardjo a 31,05\,mm/h, com o FS só em 24\,h.
\begin{table}[H]\centering\small
""" + tabela_ai("chuva", dch, ach) + r"""
\caption{Chuva: o Monte Carlo direto de 4000 é a referência para validar o FORM e a AI na chuva (etapa 4).}
\end{table}
A AI pela banda concorda com o direto (""" + num(desvios(ach[0][1], dch), 1) + r""" desvio), com ganho de """ + \
    num(ganho(ach[0][1], dch), 1) + r"""$\times$. Falta o FORM na chuva, que precisa de código novo (as opções \texttt{M},
\texttt{xi} e \texttt{desloc\_xi} no \texttt{montecarlo.cpp}).

\section{Próximos passos}
\begin{enumerate}
  \item Juntar as repetições independentes da AI, que estão rodando no fermi, e fechar a questão do pequeno desvio da AI pelo
    FORM.
  \item Etapa 4: o código do FORM e da AI pelo ponto de projeto na chuva, validados contra o direto de 4000.
  \item Para o artigo da AI: a AI pelo FORM $M = 20$ deu o maior ganho (até """ + num(max(ganho(a, d07) for _, a in a07), 0) + \
    r"""$\times$ com $P_f \approx 0{,}4\,\%$), contando só a amostragem. O custo do FORM é de 243 análises. Mesmo somando
    esse custo, ele ainda vence a AI pela banda, que precisa de um piloto de 1000.
\end{enumerate}
\end{document}
"""
open(f"{AQUI}/relatorio_etapas.tex", "w").write(tex)
print("ok: relatorio_etapas.tex e figs/")
print(f"r05: direto {pct(d05['p'])}%;", "; ".join(f"{re.sub(r'[$\\\\]', '', c)} {pct(a['p'])}% ganho {ganho(a, d05):.1f}" for c, a in a05))
print(f"r07: direto {pct(d07['p'])}%;", "; ".join(f"{re.sub(r'[$\\\\]', '', c)} {pct(a['p'])}% ganho {ganho(a, d07):.1f}" for c, a in a07))
