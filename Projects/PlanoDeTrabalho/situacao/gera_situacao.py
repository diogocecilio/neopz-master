#!/usr/bin/env python3
"""gera_situacao.py: o relatorio de situacao (o plano de trabalho, os testes e os resultados, as hipoteses, os proximos
passos e o mapa dos documentos), com os numeros lidos das campanhas baixadas do fermi e das rodadas do PC. A formulacao
fica em docs/confiabilidade_campos/formulacao: este relatorio cita as secoes dela pelo .aux (pacote xr). O mapa sai da
lista unica ../documentos.py, a mesma do indice do README. Para atualizar os dois PDFs (e o das etapas 1 a 3):
    Projects/PlanoDeTrabalho/situacao/atualiza.sh
"""
import glob
import json
import math
import os
import re
import runpy
import shutil

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

AQUI = os.path.dirname(os.path.abspath(__file__))  # Projects/PlanoDeTrabalho/situacao
PLANO = os.path.dirname(AQUI)
NEOPZ = os.path.dirname(os.path.dirname(PLANO))
ETAPAS = f"{PLANO}/ETAPAS1a3"
RUNS = os.path.expanduser("~/projects/runs")
os.makedirs(f"{AQUI}/figs", exist_ok=True)

# os numeros das etapas 1 a 3 (o mesmo gerador do relatorio das etapas: le os dados baixados do fermi)
cwd = os.getcwd()
os.chdir(ETAPAS)
E = runpy.run_path(f"{ETAPAS}/relatorio_etapas.py")
os.chdir(cwd)
for f in ("srm_hist.pdf", "pf_r05.pdf", "pf_r07.pdf", "modo.pdf"):
    shutil.copy(f"{ETAPAS}/figs/{f}", f"{AQUI}/figs/{f}")
pct, num, ai, direto, desvios, ganho = E["pct"], E["num"], E["ai"], E["direto"], E["desvios"], E["ganho"]
d05, a05, d07, a07, dch, ach = E["d05"], E["a05"], E["d07"], E["a07"], E["dch"], E["ach"]

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


# media das sementes do FORM M = 20 (independentes)
def media_sementes(casos):
    v = [a for a in casos if a]
    if not v:
        return None
    p = sum(a["p"] for a in v) / len(v)
    sd = math.sqrt(sum(a["sd"] ** 2 for a in v)) / len(v)
    return dict(p=p, sd=sd, k=len(v))


SEM05 = [E["ai"](f"{E['E1']}/r05_ai_form20")] + [E["ai"](f"{E3}/r05_ai_form20_s{k}") for k in (4, 5)
                                                  if glob.glob(f"{E3}/r05_ai_form20_s{k}/mc_*_*.csv")]
SEM07 = [E["ai"](f"{E3}/r07_ai_form20")] + [E["ai"](f"{E3}/r07_ai_form20_s{k}") for k in (4, 5)
                                             if glob.glob(f"{E3}/r07_ai_form20_s{k}/mc_*_*.csv")]
M05, M07 = media_sementes(SEM05), media_sementes(SEM07)
MODAL = [E["ai"](f"{E3}/{c}") for c in ("ai_form100_modal_s6", "ai_form200_modal_s7")
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
ax.plot([10, 20], [E["f05"][0]["bc"], E["f05"][1]["bc"]], "s--", color="0.4", label="analítica, diferenças centrais")
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
    ax.set_xlabel(r"$t = \mathbf{y}\cdot\mathbf{y}^*/\beta$")
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
        return f"{rot} & --- & na fila & --- & --- & --- \\\\"
    return (f"{rot} & {a['n']} & {pct(a['p'])} & {num(100 * a['cov'], 1)} & {num(desvios(a, d), 1)} & "
            f"{num(ganho(a, d), 1)} \\\\")


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
    return destaca_colunas(linhas, {3: "min", 4: "abs", 5: "max"}, completa)


def tabela(cab, linhas, fmt):
    return ("\\resizebox{\\ifdim\\width>\\textwidth\\textwidth\\else\\width\\fi}{!}{%\n\\begin{tabular}{" + fmt + "}\n\\toprule\n" + cab + " \\\\\n\\midrule\n" +
            "\n".join(linhas) + "\n\\bottomrule\n\\end{tabular}}")


T_R05 = tabela("Método & $N$ & $P_f$ (\\%) & CoV (\\%) & dif. (desvios) & ganho",
               destaca_ai([f"Monte Carlo direto (fermi) & {d05['n']} & {pct(d05['p'])} & {num(100 * d05['cov'], 1)} & --- & 1 \\\\"] +
               [linha_ai(r, a, d05) for r, a in a05] + [linha_ai(r, a, d05) for r, a in pc]), "lrrrrr")
T_R07 = tabela("Método & $N$ & $P_f$ (\\%) & CoV (\\%) & dif. (desvios) & ganho",
               destaca_ai([f"Monte Carlo direto (fermi) & {d07['n']} & {pct(d07['p'])} & {num(100 * d07['cov'], 1)} & --- & 1 \\\\"] +
               [linha_ai(r, a, d07) for r, a in a07]), "lrrrrr")
T_CH = tabela("Método & $N$ & $P_f$ (\\%) & CoV (\\%) & dif. (desvios) & ganho",
              destaca_ai([f"Monte Carlo direto (fermi) & {dch['n']} & {pct(dch['p'])} & {num(100 * dch['cov'], 1)} & --- & 1 \\\\"] +
              [linha_ai(r, a, dch) for r, a in ach]), "lrrrrr")
T_FORM = tabela("$M$ & variância retida (Nyström) & $\\beta$ analítica & $\\beta$ Nyström & análises por FORM",
                [f"{M} & {num(VAR_NY[M], 3)} & " + (num(FM[M][0]['b'], 3) if FM[M][0] else "---") + " & " +
                 (num(FM[M][1]['b'], 3) if FM[M][1] else "---") + " & " +
                 f"{(FM[M][1] or FM[M][0] or {'n': '---'})['n']} \\\\" for M in Ms], "rrrrr")
L = E["LINTAB"]
T_LIN = tabela("$\\beta$ & $P_f$ exata & direto ($10^6$) & AI na direção exata & AI a " +
               num(math.degrees(math.acos(L[1])), 1) + "$^\\circ$ da direção",
               destaca_linhas([f"{x['b']} & {x['ex']:.3e} & " + (f"{pct(x['d'])}\\,\\%" if x['d'] else "---") +
                f" & {num(100 * (x['ai20'][0] / x['ex'] - 1), 1)}\\,\\% ({num(x['ai20'][2], 1)}) & "
                f"{num(100 * (x['ai10'][0] / x['ex'] - 1), 1)}\\,\\% ({num(x['ai10'][2], 1)})" + " \\\\" for x in L[0]],
                [3, 4], "abs"),
               "rrrrr") if L else "(teste linear não baixado)"
C05, C07 = E["CUSTO05"], E["CUSTO07"]
T_CUSTO = tabela("Método & $r = -0{,}5$ & $r = -0{,}7$",
                 destaca_colunas([f"{a[0]} & {a[3]:.0f} & {b[3]:.0f} \\\\" for a, b in zip(C05, C07)],
                                 {1: "min", 2: "min"}, lambda i, c: i > 0), "lrr")
cho, chuva, mc = E["cho"], E["chuva"], E["mcnemar"]


def tabela_vies(r, d, dec):
    """as AI pelo ponto de projeto de um r: o diagnostico pelo direto e a estimativa com a variavel de controle"""
    lin = [f"Monte Carlo direto & --- & 2010 & {d['n']} & {pct(d['p'], dec)} & {num(100 * d['sd'] / d['p'], 1)} & --- & --- & "
           f"{pct(d['cv']['p'], dec)} & {num(100 * d['cv']['sd'] / d['cv']['p'], 1)} \\\\"]
    for a in DG["ai"]:
        if a["r"] != r:
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
MAPA = runpy.run_path(f"{PLANO}/documentos.py", run_name="documentos")["mapa_latex"]()


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
S07 = [s for rot, v, d, s in RS if v is KL07]
S07 = S07[0] if S07 else None

S05 = [s for rot, v, d, s in RS if v is KL05]
S05 = S05[0] if S05 else None
DIF_T0 = (AT0["p"] - DD5["p"]) / math.hypot(AT0["cov_verdadeiro"] * DD5["p"], DD5["sd"])
HIP05 = (r"  \item \textbf{A AI pelo FORM não tem viés} ($r = -0{,}5$), nas duas bases (Seção~\ref{sec:vies}). Com a semente 3, "
         r"as AI saem todas em ≈1,55\,\% (os sorteios são os mesmos). As " + str(S05["k"]) + " sementes independentes "
         r"da KL dão, em média, " + pct(S05["p"], 3) + r"\,\% contra " + pct(DD5["p"], 3) + r"\,\% do direto e, com a "
         r"variável de controle, " + pct(S05["pc"], 3) + r" $\pm$ " + pct(S05["sdc"], 3) + r"\,\% contra " +
         pct(DD5["cv"]["p"], 3) + r" $\pm$ " + pct(DD5["cv"]["sd"], 3) + r"\,\% ($" + num(S05["difc"], 1) + "$ desvio). "
         r"A direção errada da Tri(0) é a única que subestima o próprio CoV (" + num(100 * AT0["sd"] / AT0["p"], 1) +
         r" contra " + num(100 * AT0["cov_verdadeiro"], 1) + r"\,\%); com o CoV verdadeiro, o seu desvio cai de $-2{,}2$ "
         r"para $" + num(DIF_T0, 1) + "$.")
HIP07 = (r"  \item \textbf{A AI pelo FORM com $r = -0{,}7$ também não tem viés detectável:} as sementes 3, 4 e 5, corrigidas, "
         r"dão " + ", ".join(pct(a["p_cv"], 3) for a in A07) + r"\,\%, contra " + pct(DD7["cv"]["p"], 3) + r" $\pm$ " +
         pct(DD7["cv"]["sd"], 3) + r"\,\% do direto corrigido (" + f"{DD7['n']:,}".replace(",", "\\,") + " realizações, " +
         str(DD7["falhas"]) + " falhas): $" + (num(S07["difc"], 1) if S07 else "---") + "$ desvio.")

tex = (r"""\documentclass[11pt,a4paper]{article}
\usepackage[utf8]{inputenc}
\DeclareUnicodeCharacter{2248}{\ensuremath{\approx}}
\usepackage[T1]{fontenc}
\usepackage[brazil]{babel}
\usepackage{lmodern}
\usepackage[margin=2cm]{geometry}
\usepackage{amsmath,amssymb,bm}
\usepackage{graphicx,booktabs,float,longtable,url}
\usepackage{xcolor}
\definecolor{azul}{RGB}{0,70,200}
\definecolor{verm}{RGB}{200,0,0}
\usepackage{xr-hyper}
\usepackage[hidelinks]{hyperref}
% as secoes da formulacao (docs/confiabilidade_campos/formulacao), citadas como F-<rotulo>: compile a formulacao antes
\externaldocument[F-]{../../../docs/confiabilidade_campos/formulacao/formulacao}
% a data e o commit vem do atualiza.sh; compilado a mao, fica so a data
\providecommand{\versao}{17/11/2026}
\title{Confiabilidade de taludes com campos aleatórios no NeoPZ\\\large Situação: plano de trabalho, testes e resultados,
hipóteses e próximos passos}
\date{\versao}
\begin{document}
\maketitle
\tableofcontents

\section{Objetivo e onde estamos}
O trabalho serve a dois artigos:
\begin{itemize}
  \item \textbf{Amostragem por importância (AI)} para taludes com campos aleatórios no MEF. A validação é no Ceron (2025),
    no rebaixamento, e no Cho (2010), um talude seco com $c'$ e $\varphi'$ como campos correlacionados e $P_f$ de 0,4\,\% a
    2\,\%.
  \item \textbf{Chuva em talude não saturado} (Rahardjo et al., 2007), com as mesmas ferramentas: Monte Carlo, FORM e AI.
\end{itemize}
O plano de trabalho (\path{Projects/PlanoDeTrabalho/README.md}) valida as ferramentas uma de cada vez, antes de usá-las
na chuva. As etapas 0 a 3 estão rodadas, e a 2c (gradiente modal) foi testada no PC. A etapa 4 (FORM e AI na chuva) tem
a referência pronta, e a etapa 5 (a investigação da chuva) ainda não começou. Tudo roda no cluster fermi (4 nós, 48
processos); o PC serve para testes curtos.

A formulação (o modelo, as duas bases de Karhunen--Loève, o Monte Carlo, a AI, a variável de controle e o FORM) está em
\path{docs/confiabilidade_campos/formulacao/formulacao.pdf} e é citada aqui como ``formulação, Seção $n$''. Este
relatório traz só os resultados. A Seção~\ref{sec:mapa} lista todos os documentos.

\section{O plano de trabalho e a situação}
\begin{table}[H]\centering\small
""" + tabela("Etapa & O quê & Situação",
             ["0 & Cho2010 no fermi: direto, AI pela banda e pelo FORM & pronta \\\\",
              "1 & SRM corrigido ($F_0 = 0{,}8$): o antigo contra o novo nas mesmas realizações & pronta \\\\",
              "2 & FORM ($M = 10$, 20) e AI pelo ponto de projeto; teste linear sem MEF & pronta; sem viés (Seção~\\ref{sec:vies}) \\\\",
              "2b & base de Nyström contra a analítica & o mesmo ponto de projeto (Seção~\\ref{sec:pp}); AI de Nyström no fermi (fase 4, pausada) \\\\",
              "2c & gradiente do FORM pela projeção modal & validado; a AI modal não tem viés; 10 sementes no fermi (fase 4, pausada) \\\\",
              "3 & $P_f$ pequena ($r = -0{,}7$), tamanho do piloto, custo & o direto vai de 25\\,000 (93 falhas) a 110\\,000 (CoV 5\\,\\%), com os 4 nós do fermi \\\\",
              "4 & FORM e AI na chuva & o direto de 10\\,000 está terminando; falta o código do FORM na chuva \\\\",
              "5 & investigação da chuva & não começou \\\\"], "lp{0.55\\textwidth}p{0.3\\textwidth}") + r"""
\end{table}

\section{Testes e resultados}
\textbf{Cores nas tabelas:} em cada coluna de desempenho, \textcolor{azul}{\textbf{azul}} é o melhor valor e
\textcolor{verm}{\textbf{vermelho}} o pior. CoV e custo: o menor é o melhor; ganho: o maior; dif.: o menor em módulo. Só
entram as AI completas (1000 realizações): o direto é a referência, e as parciais e as da fila ficam de fora. No teste
linear, a comparação é entre as duas direções, linha a linha. Atenção: um $|$dif.$|$ abaixo de 2 é compatível com o
acaso; o azul na dif. não quer dizer que o método é mais exato, e o vermelho só é um sinal quando passa de 2.

\subsection{SRM corrigido (etapa 1)}
Com as mesmas realizações nas duas versões:
\begin{itemize}
  \item \textbf{Cho2010:} as falhas são idênticas (""" + f"{cho['fa']} e {cho['fb']}" + r""" em """ + f"{cho['n']}" + r""").
  \item \textbf{Chuva:} """ + f"{chuva['n01']}" + r""" realizações passam a falhar e """ + f"{chuva['n10']}" + \
    r""" deixa de falhar (McNemar $p = """ + num(mc(chuva['n01'], chuva['n10']), 2) + r"""$), cerca de $+0{,}5$ ponto.
\end{itemize}
O efeito é pequeno, e o SRM novo é o correto. As duas versões do programa reproduzem o mesmo FS (diferença de 0,0000), e o
PC reproduz o fermi.
\begin{figure}[H]\centering\includegraphics[width=0.9\textwidth]{figs/srm_hist.pdf}\caption{FS perto de 1 nas duas
versões do SRM.}\end{figure}

\subsection{Cho2010, $r = -0{,}5$: AI contra o direto (etapas 0, 2 e 2c)}
\begin{table}[H]\centering\small
""" + T_R05 + r"""
\caption{$r = -0{,}5$ (Cho: $P_f = 1{,}71\,\%$). O direto terminou: 24\,000 realizações. ``dif.'' é a diferença AI
$-$ direto em desvios. Sem indicação, o FORM usa a semente 3 e a banda a semente 1; [PC] e [fermi] dizem onde rodou.}
\end{table}
\begin{figure}[H]\centering\includegraphics[width=0.75\textwidth]{figs/pf_r05.pdf}\caption{$P_f$ e IC 95\,\%, $r = -0{,}5$.}\end{figure}

\subsection{$P_f$ pequena, $r = -0{,}7$ (etapa 3)}
\begin{table}[H]\centering\small
""" + T_R07 + r"""
\caption{$r = -0{,}7$ (Cho: $P_f \approx 0{,}4\,\%$).}
\end{table}

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

\subsection{A AI pelo FORM tem viés? O teste com o próprio direto (KL e Nyström)}\label{sec:vies}
Com a semente 3, todas as AI pelo FORM ficaram abaixo do direto (1,53 a 1,57\,\% contra 1,70\,\%), e a AI modal, a mais
eficiente, parecia convergir para o valor errado. Uma AI com deslocamento fixo não tem viés, mas pode se comportar como se
tivesse. Se a densidade de amostragem quase não visita uma parte da região de falha (outro modo de ruptura, longe do
ponto de projeto), a estimativa sai quase sempre baixa, e o CoV estimado fica menor que o verdadeiro. O direto permite
testar isso sem nenhuma análise nova (\path{Projects/Cho2010/diagnostico_ai.py}):
\begin{enumerate}
  \item \textbf{Os $\bm y$ de cada realização} saem de uma cópia do gerador do \texttt{RandomField} em Python
    (xoshiro256** e Box--Muller, uma semente por realização e campo). A cópia confere com o peso de cada amostra de cada
    AI (erro $\le """ + cient(max(a["erro_peso"] for a in DG["ai"])) + r"""$). Junto com a KL analítica em Python, ela
    também confere com o campo que o programa grava (diferença $\le 5\times10^{-6}$, a precisão do CSV).
  \item \textbf{Modo de ruptura fora da AI.} As """ + str(DD5["falhas"]) + r""" falhas do direto, projetadas na direção
    do ponto de projeto ($t = \bm y\cdot\bm y^*/\beta$), ficam todas onde a AI amostra (Figura~\ref{fig:vies_t}). Na direção
    modal, o menor $t$ é """ + num(AM["t_min"], 2) + r""", e só """ + str(AM["abaixo_15"]) + r""" falhas ficam abaixo de
    1,5. A variância verdadeira de cada AI sai do direto pela identidade $E_h[I w^2] = E_\varphi[I w]$, com as falhas do
    direto, e é a que a AI estima (``CoV verd.'' na Tabela~\ref{tab:vies05}): """ + num(100 * AM["cov_verdadeiro"], 1) +
    r"""\,\% contra """ + num(100 * AM["sd"] / AM["p"], 1) + r"""\,\% no modal $M = 100$. Só a direção errada da Tri(0)
    subestima o próprio CoV (""" + num(100 * AT0["sd"] / AT0["p"], 1) + r"""\,\% estimado contra """ +
    num(100 * AT0["cov_verdadeiro"], 1) + r"""\,\% verdadeiro): o $-2{,}2$ desvio da Tabela~1 vinha daí.
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
\caption{A chance de falha do direto em cada faixa de $t = \bm y\cdot\bm y^*/\beta$ (barras) e as densidades de $t$ no
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

\subsection{Comparação das bases com a mesma variância}
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

\subsection{O piso do SRM em $F_0 = 0{,}8$}
Nenhuma realização das AI deu FS abaixo de 0,8 (o mínimo foi 0,802). É o mesmo mecanismo do buraco em $F = 1$: um talude
com FS verdadeiro logo abaixo de $F_0$ às vezes converge o peso com $F_0$ e devolve $\approx F_0$. Para a $P_f$ não muda
nada (o FS continua $< 1$), mas a cauda inferior do FS fica distorcida (≈1 realização em 1000). Correção proposta:
recomeçar de $F_0/2$ quando o caminho de $F$ colapsa logo no início.

\subsection{Custo total para CoV$(P_f) = 10\,\%$}
\begin{table}[H]\centering\small
""" + T_CUSTO + r"""
\caption{Análises de elementos finitos, com o piloto (1000) e o FORM incluídos. Com o piloto na malha grossa (custo
1/20, Seção~\ref{sec:hip}), a AI pela banda cai para ≈1300 e ≈2700.}
\end{table}

\subsection{O teste linear (sem elementos finitos)}
\begin{table}[H]\centering\small
""" + T_LIN + r"""
\caption{Falha quando $\bm e\cdot\bm y > \beta$, com $\bm e$ unitário e $P_f = \Phi(-\beta)$ exata. Erro relativo, com o erro em desvios entre parênteses.}
\end{table}
Com o deslocamento na direção exata, a AI acerta até $3\times10^{-5}$. Com a direção um pouco errada, ela sai baixa, e
mais baixa quando $\beta$ cresce.

\subsection{FORM: gradiente modal e as duas bases (etapas 2b e 2c)}\label{sec:nystrom}
\begin{itemize}
  \item \textbf{Validação do gradiente modal} nos pontos de projeto do fermi: ele fica a 9,6$^\circ$ ($M = 10$) e
    8,9$^\circ$ ($M = 20$) do gradiente das diferenças centrais, e o termo dominante ($\varphi'_0$) bate a menos de 4\,\%.
\end{itemize}
\begin{table}[H]\centering\small
""" + T_FORM + r"""
\caption{$\beta$ do FORM com o gradiente modal na malha usual. Com diferenças centrais: 2,115 ($M = 10$, 123 análises) e
2,012 ($M = 20$, 243).}
\end{table}
\begin{figure}[H]\centering\includegraphics[width=0.6\textwidth]{figs/beta_M.pdf}\caption{$\beta(M)$: a base de
Nyström converge com menos termos.}\end{figure}
A base de Nyström chega mais perto do $\beta$ convergido com poucos termos ($M = 10$: 2,062 contra 2,151), e as duas
coincidem a partir de $M = 50$. É o mesmo efeito da placa curva do estudo do Ghanem \& Spanos (\path{docs/confiabilidade_campos/kl_form_placas}). O $\beta$ converge para
$\approx 2{,}00$, mas a $P_f$ do FORM sozinho ($\approx 2{,}3\,\%$) superestima a do direto em ≈30\,\%: é o erro de
linearizar uma superfície de falha curva. Por isso o FORM entra só como deslocamento da AI.

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
tem $|\bm u| = """ + num(PP["beta_kl_em_nys"], 3) + r"""$ e fica a """ + num(PP["angulo"], 1) + r"""$^\circ$ do de
Nyström. O campo mais provável de ruptura é uma camada fraca horizontal sob a face do talude:
\begin{itemize}
  \item $\varphi'$ mínimo de 23,4$^\circ$ em $y \approx 6{,}6$\,m, 1,6\,m acima do pé;
  \item $c'$ quase uniforme, de 8,6 a 9,6\,kPa;
  \item os perfis de $\varphi'$ das duas bases diferem em menos de 0,1$^\circ$.
\end{itemize}
\begin{figure}[H]\centering\includegraphics[width=\textwidth]{figs/ponto_projeto.pdf}
\caption{$c'$ e $\varphi'$ no ponto de projeto do FORM modal: KL analítica ($M = 100$, em cima) e Nyström ($M = 200$,
embaixo), com as mesmas escalas de cor.}\label{fig:pp}\end{figure}

\subsection{A referência da chuva (etapa 4)}
\begin{table}[H]\centering\small
""" + T_CH + r"""
\caption{f50 a 31,05\,mm/h, FS em 24\,h.}
\end{table}

\section{Hipóteses}\label{sec:hip}
\subsection*{Confirmadas}
\begin{itemize}
  \item A AI não tem viés com qualquer deslocamento fixado antes; o deslocamento só muda a variância (teste linear e Ceron).
  \item O SRM antigo tinha um buraco logo abaixo de FS = 1; o efeito na $P_f$ é pequeno.
  \item \textbf{Há um só modo de ruptura} no Cho2010, nas duas bases: todas as falhas do direto ficam onde a AI pelo
    FORM amostra ($t = \bm y\cdot\bm y^*/\beta \ge """ + num(AM["t_min"], 2) + r"""$ na direção modal, com só
    """ + str(AM["abaixo_15"]) + r""" das """ + str(DD5["falhas"]) + r""" abaixo de 1,5;
    $t \ge """ + num(A07[0]["t_min"], 2) + r"""$ com $r = -0{,}7$), e a variância verdadeira da AI,
    calculada pelas falhas do direto, é a que a AI estima (Seção~\ref{sec:vies}). A banda de ruptura das falhas é a mesma
    no direto e na AI.
  \item \textbf{O ponto de projeto é o mesmo nas duas bases:} levado para as coordenadas de Nyström, o ponto analítico
    fica a """ + num(PP["angulo"], 1) + r"""$^\circ$ do de Nyström, e os campos de $\varphi'$ diferem em
    menos de 0,1$^\circ$ (Figura~\ref{fig:pp}).
  \item \textbf{A variável de controle dá um ganho extra sem custo:} de """ + num(min(GANHOS), 1) + r""" a
    """ + num(max(GANHOS), 1) + r""" nas AI de $r = -0{,}5$ e """ + num(DD5["n_equiv"] / DD5["n"], 1) + r"""
    no direto, cuja referência vai a """ + pct(DD5["cv"]["p"], 3) + r""" $\pm$ """ + pct(DD5["cv"]["sd"], 3) + r"""\,\%.
  \item \textbf{O piloto pode ser pequeno:} com 250 realizações, o deslocamento da banda custa só ≈5\,\% de variância.
  \item \textbf{O piloto pode rodar em malha grossa:} o elite é escolhido por ordem, e a ordem das realizações quase não
    muda. Com o piloto na Tri(0) (300 vezes mais barata), a AI pela banda acertou o direto ($-0{,}2$ desvio), com o mesmo
    CoV do piloto na malha usual. Tri(1) e Tri(2) também estão nas tabelas.
  \item \textbf{O gradiente modal dá a direção do FORM com 1 análise,} a menos de ≈9$^\circ$ das diferenças centrais, e o
    FORM converge em 3 análises.
  \item \textbf{A base de Nyström precisa de menos termos} que a analítica para o mesmo $\beta$.
""" + HIP05 + r"""
""" + HIP07 + r"""
\end{itemize}
\subsection*{Descartadas}
\begin{itemize}
  \item \textbf{O FORM na malha Tri(0)} (76 equações, FS 1,47 contra 1,20): $\beta = 4{,}08$ e direção a 31,7$^\circ$ da
    certa. O erro da malha desloca o ponto de projeto, e não só o FS. Com a direção da Tri(0) e a distância da malha usual
    ($\beta = 2{,}32$), a AI funciona, mas pior que o FORM na malha usual.
\end{itemize}
\subsection*{Sob investigação}
\begin{itemize}
  \item \textbf{O campo de $k$ na chuva:} o gradiente modal pega $c'$ e $\varphi'$; o $k$ age pelo transiente e
    precisaria do adjunto. Hoje o $k$ pesa pouco na chuva f50 ($\alpha^2_k \approx 0$), e a AI continua sem viés sem
    deslocá-lo.
\end{itemize}

\section{Próximos passos}
\begin{enumerate}
  \item \textbf{Direto de $r = -0{,}7$ até 110\,000} (CoV de 5\,\%; com os 4 nós do fermi, ≈33\,h). A fase 4 fica
    pausada até ele acabar.
  \item \textbf{Fase 4 no fermi} (≈8\,h): 10 sementes do FORM modal analítico e 10 do de Nyström, as mesmas nas duas
    bases. Dão a distribuição do estimador e fecham a questão do viés (Seção~\ref{sec:vies}).
  \item \textbf{Etapa 4:} o FORM com gradiente modal e a AI pelo ponto de projeto na chuva, contra o direto de 10\,000
    (terminando no fermi).
  \item \textbf{Variável de controle} nas AI pela banda e na chuva. Para o artigo da AI: o ganho contra o custo total
    (FORM, piloto e controle), de $P_f \approx 2\,\%$ a $0{,}4\,\%$.
\end{enumerate}

\section{Mapa dos documentos}\label{sec:mapa}
É a mesma lista do índice do \path{Projects/PlanoDeTrabalho/README.md}: as duas saem de
\path{Projects/PlanoDeTrabalho/documentos.py}. Os caminhos relativos são a partir de \path{~/projects/neopz-master}.
""" + MAPA + r"""
\end{document}
""")
open(f"{AQUI}/situacao.tex", "w").write(tex)
print("ok: situacao.tex e figs/")
