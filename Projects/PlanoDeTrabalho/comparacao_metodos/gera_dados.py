#!/usr/bin/env python3
"""gera_dados.py: le os resultados parciais (Monte Carlo direto, ISMC e line sampling) dos tres problemas de talude e da
viga, e grava as figuras (figs/*.pdf) e as tabelas e numeros do relatorio (tabelas.tex, numeros.tex).

Fontes (no PC, em ~/projects/runs; a viga vem dos relatorios em Mathematica):
  Cho r = -0,5 e -0,7: direto r05 (fermi_etapa1/etapa1/r05_direto, 24 000) e r07 (fermi_etapa3/etapa3/r07_direto,
      110 000); ISMC da etapa 3b (fermi_etapa3b/etapa3b/r0x_nys_<centro>_s21..25, base de Nystrom de 2000 modos); FORM
      da etapa 3b (form_r0x_<M20|M100|gim20>/form.json); line sampling no PC (etapa3b_ls/ls/r0x_<dir>/linhas.csv) com as
      direcoes etapa3b_ls/form_r0x_<gim2000|M20>/form.json.
  Ceron (rebaixamento, Pf ~16%): direto de 2000 (ceron_transiente/direto/ceron_tr_fase1/f50_I0_iii); ISMC
      (ceron_form/ismc_dif, ismc_modal); line sampling (ceron_form/ls_dif, ls_modal, ls_modal2000) e os FORM
      (ceron_form/form_*/form.json).
  Viga: MC direto e ISMC do relatorio da viga (viga_mathematica/relatorio_viga.tex, "O ponto otimo do M do FORM");
      line sampling e FORM adjunto das secoes 16 e 17 do codigo viga_mathematica/codigo/viga_FORM5_AI200.txt (saida do
      Mathematica de 19/11/2026, ainda fora do relatorio da viga).
Uso: python3 -I gera_dados.py   (roda em segundos; refaz tudo com o que houver de dados)."""
import csv
import json
import math
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

AQUI = os.path.dirname(os.path.abspath(__file__))
RAIZ = os.path.abspath(os.path.join(AQUI, "..", "..", ".."))
sys.path.insert(0, os.path.join(RAIZ, "Projects", "Cho2010"))
import ismc_parada as ip  # noqa: E402

R = os.path.expanduser("~/projects/runs")
FIGS = os.path.join(AQUI, "figs")
os.makedirs(FIGS, exist_ok=True)
ALVO = 0.05
SEM = [21, 22, 23, 24, 25]


def cauda(c):
    return 0.5 * math.erfc(c / math.sqrt(2.))


def form(arq):
    d = json.load(open(arq))
    return {"beta": d["beta"], "n": d["avaliacoes"], "M": d["parametros"]["M"],
            "grad": d["parametros"].get("gradiente", "dif")}


def ls(pasta):
    """o line sampling: as linhas na ordem em que terminaram"""
    arq = os.path.join(pasta, "linhas.csv")
    if not os.path.exists(arq):
        return None
    rs = list(csv.DictReader(open(arq)))
    if len(rs) < 2:
        return None
    c = [float(r["c"]) for r in rs]
    na = [int(r["analises"]) for r in rs]
    p = [cauda(x) for x in c]
    n = len(p)
    m = sum(p) / n
    sd = math.sqrt(sum((x - m) ** 2 for x in p) / (n - 1))
    s = json.load(open(os.path.join(pasta, "inclinacao.json")))
    with open(os.path.join(pasta, "analises.csv")) as f:
        feitas = sum(1 for x in f if x.strip())   # todas as analises, inclusive as das linhas ainda em andamento
    return {"nu": c, "na": na, "p": p, "N": n, "Pf": m, "CoV": sd / (math.sqrt(n) * m), "CoVlinha": sd / m,
            "porlinha": sum(na) / n, "analises": max(feitas, sum(na) + 2), "beta": s["beta"], "s": s["s"]}


def ismc(pastas):
    """as sementes: a P_f de cada uma e todas juntas; o custo para CoV 5% por semente (media)"""
    sem = []
    for pa in pastas:
        if not os.path.isdir(pa):
            continue
        iw, _ = ip.serie(pa)
        if len(iw) < 10:
            continue
        p, cv, f, ess, maior = ip.estatisticas(iw)
        sem.append({"iw": iw, "N": len(iw), "Pf": p, "CoV": cv})
    if not sem:
        return None
    todos = [x for s in sem for x in s["iw"]]
    p, cv, f, ess, maior = ip.estatisticas(todos)
    n5 = sum(s["N"] * (s["CoV"] / ALVO) ** 2 for s in sem) / len(sem)
    return {"sem": sem, "N": len(todos), "Pf": p, "CoV": cv, "n5": n5, "nsem": len(sem)}


def direto(pasta):
    d = ip.direto(pasta)
    d["n5"] = (1 - d["Pf"]) / (d["Pf"] * ALVO ** 2)
    return d


def z(p, cv, ref):
    return (p - ref["Pf"]) / math.sqrt((p * cv) ** 2 + (ref["Pf"] * ref["CoV"]) ** 2)


# ---------------- os dados ----------------
E3B = os.path.join(R, "fermi_etapa3b", "etapa3b")
LSC = os.path.join(R, "etapa3b_ls")
CER = os.path.join(R, "ceron_form")
P = {}
for r, rotulo, ref in (("r05", "Cho, $r = -0{,}5$", os.path.join(R, "fermi_etapa1", "etapa1", "r05_direto")),
                       ("r07", "Cho, $r = -0{,}7$", os.path.join(R, "fermi_etapa3", "etapa3", "r07_direto"))):
    pr = {"rotulo": rotulo, "ref": direto(ref), "ismc": [], "ls": []}
    for cen, nome in (("M20", "FORM $M = 20$, numérico"), ("M100", "FORM $M = 100$, numérico"),
                      ("gim20", "FORM $M = 20$, modal")):
        fo = form(os.path.join(E3B, f"form_{r}_{cen}", "form.json"))
        d = ismc([os.path.join(E3B, f"{r}_nys_{cen}_s{s}") for s in SEM])
        if d:
            pr["ismc"].append({"nome": nome, "form": fo, **d})
    for dirc, nome in (("M20", "FORM $M = 20$, numérico"), ("gim2000", "FORM $M = 2000$, modal")):
        fo = form(os.path.join(LSC, f"form_{r}_{dirc}", "form.json")) if os.path.exists(
            os.path.join(LSC, f"form_{r}_{dirc}", "form.json")) else None
        d = ls(os.path.join(LSC, "ls", f"{r}_{dirc}"))
        if d and fo:
            pr["ls"].append({"nome": nome, "form": fo, **d})
    P[r] = pr
pc = {"rotulo": "Ceron, rebaixamento", "ref": direto(os.path.join(R, "ceron_transiente", "direto", "ceron_tr_fase1",
                                                                   "f50_I0_iii")), "ismc": [], "ls": []}
for cen, nome in (("dif", "FORM $M = 20$, numérico"), ("modal", "FORM $M = 20$, modal")):
    fo = form(os.path.join(CER, f"form_{cen}", "form.json"))
    d = ismc([os.path.join(CER, f"ismc_{cen}")])
    pc["ismc"].append({"nome": nome, "form": fo, **d})
for dirc, fdir, nome in (("dif", "form_dif", "FORM $M = 20$, numérico"), ("modal", "form_modal", "FORM $M = 20$, modal"),
                         ("modal2000", "form_modal2000", "FORM $M = 2000$, modal")):
    fo = form(os.path.join(CER, fdir, "form.json"))
    d = ls(os.path.join(CER, f"ls_{dirc}"))
    if d:
        pc["ls"].append({"nome": nome, "form": fo, **d})
P["ceron"] = pc
ORDEM = ["ceron", "r05", "r07"]
# a viga (Mathematica; relatorio_viga e secoes 16-17 do codigo): Pf, CoV da rodada, total p/ CoV 5% (com o FORM)
CV_VIGA = math.sqrt((1 - 0.02098) / (50000 * 0.02098))
VIGA = {"ref": {"N": 50000, "falhas": 1049, "Pf": 0.02098, "CoV": CV_VIGA, "n5": 18666},
        "linhas": [("MC direto", "--", 50000, 0.02098, CV_VIGA, 18666),
                   ("ISMC, centro do FORM $M = 8$ (numérico)", 85, 1085, 0.0229, None, 1198),
                   ("ISMC, centro do FORM $M = 20$ (numérico)", 205, 1205, 0.0220, None, 1224),
                   ("LS, FORM $M = 8$ (numérico)", 85, 1099, 0.02155, 0.0219, 281),
                   ("LS, FORM $M = 20$ (numérico)", 205, 1209, 0.02104, 0.0135, 280),
                   ("LS, FORM $M = 200$ (numérico)", 2005, 3008, 0.02045, 0.0059, 2021),
                   ("LS, FORM $M = 200$ (adjunto)", 15, 1021, 0.02107, 0.0137, 92)]}

# ---------------- as linhas de comparacao ----------------
def linhas_problema(k):
    pr = P[k]
    ref = pr["ref"]
    out = [{"metodo": "MC direto", "dir": "--", "form": None, "gasto": ref["N"], "Pf": ref["Pf"], "CoV": ref["CoV"],
            "z": None, "n5": ref["n5"], "tipo": "mc"}]
    for d in pr["ismc"]:
        out.append({"metodo": "ISMC" + (f" ({d['nsem']} sementes)" if d["nsem"] > 1 else ""), "dir": d["nome"],
                    "form": d["form"]["n"], "gasto": d["N"] + d["form"]["n"], "Pf": d["Pf"], "CoV": d["CoV"],
                    "z": z(d["Pf"], d["CoV"], ref), "n5": d["n5"] + d["form"]["n"], "tipo": "ismc"})
    for d in pr["ls"]:
        n5 = 2 + d["porlinha"] * d["N"] * (d["CoV"] / ALVO) ** 2
        out.append({"metodo": f"LS ({d['N']} linhas{', parcial' if d['N'] < 100 else ''})", "dir": d["nome"], "form": d["form"]["n"],
                    "gasto": d["analises"] + d["form"]["n"], "Pf": d["Pf"], "CoV": d["CoV"],
                    "z": z(d["Pf"], d["CoV"], ref), "n5": n5 + d["form"]["n"], "tipo": "ls", "parcial": d["N"] < 100,
                    "nlin": d["N"]})
    melhor_ismc = min([x["n5"] for x in out if x["tipo"] == "ismc"] or [float("nan")])
    for x in out:
        x["ganho"] = ref["n5"] / x["n5"]
        x["ganho_ismc"] = melhor_ismc / x["n5"] if x["tipo"] == "ls" else None
    return out


def num(x, d=0):
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return "--"
    s = f"{x:,.{d}f}".replace(",", "X").replace(".", ",").replace("X", "\\,")
    return s


def pct(x, d=2):
    return "--" if x is None else num(100 * x, d) + "\\,\\%"


T = []
T.append("% gerado por gera_dados.py: nao editar a mao\n")
# tabela de referencias
T.append("\\begin{table}[H]\\centering\\small\n\\caption{As referências: o Monte Carlo direto de cada problema. "
         "O custo para CoV de 5\\,\\% é $(1 - P_f)/(0{,}05^2 P_f)$.}\\label{tab:ref}\n"
         "\\begin{tabular}{lrrrrr}\\toprule\nproblema & $N$ & falhas & $P_f$ & CoV & $N$ p/ CoV 5\\,\\% \\\\\\midrule\n")
for k in ORDEM:
    d = P[k]["ref"]
    T.append(f"{P[k]['rotulo']} & {num(d['N'])} & {num(d['falhas'])} & {pct(d['Pf'], 3)} & {pct(d['CoV'], 1)} & "
             f"{num(d['n5'])} \\\\\n")
T.append(f"viga (Mathematica) & {num(VIGA['ref']['N'])} & {VIGA['ref']['falhas']} & {pct(VIGA['ref']['Pf'], 3)} & "
         f"{pct(VIGA['ref']['CoV'], 1)} & {num(VIGA['ref']['n5'])} \\\\\n")
T.append("\\bottomrule\\end{tabular}\\end{table}\n\n")
# tabela dos FORM
T.append("\\begin{table}[H]\\centering\\small\n\\caption{Os FORM que dão os centros da ISMC e as direções do line "
         "sampling. ``Numérico'': diferenças centrais ($4M + 1$ análises por iteração no Cho, $6M + 1$ no Ceron); "
         "``modal'': o gradiente pela projeção modal do GIM (1 análise por iteração).}\\label{tab:form}\n"
         "\\begin{tabular}{llrrrr}\\toprule\nproblema & FORM & $\\beta$ & $\\Phi(-\\beta)$ & análises & usado em \\\\\\midrule\n")
for k in ORDEM:
    vistos = {}
    for d in P[k]["ismc"]:
        vistos.setdefault(d["nome"], [d["form"], set()])[1].add("ISMC")
    for d in P[k]["ls"]:
        vistos.setdefault(d["nome"], [d["form"], set()])[1].add("LS")
    for nome, (fo, uso) in vistos.items():
        T.append(f"{P[k]['rotulo']} & {nome} & {num(fo['beta'], 3)} & {pct(cauda(fo['beta']), 2)} & {fo['n']} & "
                 f"{', '.join(sorted(uso))} \\\\\n")
T.append("\\bottomrule\\end{tabular}\\end{table}\n\n")
# a comparacao, um problema por tabela
for k in ORDEM:
    T.append(f"\\begin{{table}}[H]\\centering\\small\n\\caption{{{P[k]['rotulo']}: todos os métodos. ``Gasto'': as "
             f"análises do MEF já feitas (com o FORM). ``Total p/ CoV 5\\,\\%'': o FORM mais as amostras ou linhas para "
             f"CoV de 5\\,\\%, pelo CoV estimado (na ISMC, a média das sementes). $z$: a diferença contra o direto, em "
             f"desvios padrão.}}\\label{{tab:{k}}}\n\\resizebox{{\\textwidth}}{{!}}{{%\n"
             "\\begin{tabular}{llrrrrrrrr}\\toprule\nmétodo & centro ou direção & FORM & gasto & $P_f$ & CoV & $z$ & "
             "total p/ CoV 5\\,\\% & ganho vs direto & ganho vs ISMC \\\\\\midrule\n")
    for x in linhas_problema(k):
        neg = x["tipo"] == "ls" and x["n5"] == min(y["n5"] for y in linhas_problema(k))
        b = (lambda s: "\\textbf{" + s + "}") if neg else (lambda s: s)
        T.append(f"{x['metodo']} & {x['dir']} & {x['form'] if x['form'] is not None else '--'} & {num(x['gasto'])} & "
                 f"{pct(x['Pf'], 3)} & {pct(x['CoV'], 1)} & {num(x['z'], 2) if x['z'] is not None else '--'} & "
                 f"{b(num(x['n5']))} & {b(num(x['ganho'], 1) + chr(215))} & "
                 f"{(num(x['ganho_ismc'], 1) + chr(215)) if x['ganho_ismc'] else '--'} \\\\\n")
    T.append("\\bottomrule\\end{tabular}}\\end{table}\n\n")
# viga
T.append("\\begin{table}[H]\\centering\\small\n\\caption{A viga com campo aleatório (Mathematica). MC direto e ISMC: "
         "relatório da viga, ``O ponto ótimo do M do FORM''; na ISMC, o total usa o CoV verdadeiro, calculado com as falhas "
         "do direto de 50\\,000. Line sampling e FORM adjunto: seções 16 e 17 do código "
         "\\texttt{viga\\_mathematica/codigo/viga\\_FORM5\\_AI200.txt}, ainda fora do relatório da viga. FORM: só as análises do "
         "FORM; ``gasto'' e ``total'' do line sampling incluem as 2 da inclinação.}"
         "\\label{tab:viga}\n\\resizebox{\\textwidth}{!}{%\n\\begin{tabular}{lrrrrrr}\\toprule\nmétodo & FORM & gasto & $P_f$ & CoV & total p/ CoV 5\\,\\% & "
         "ganho vs direto \\\\\\midrule\n")
for nome, fo, g, p, cv, n5 in VIGA["linhas"]:
    T.append(f"{nome} & {fo} & {num(g)} & {pct(p, 3)} & {pct(cv, 1)} & {num(n5)} & {num(18666 / n5, 1)}{chr(215)} \\\\\n")
T.append("\\bottomrule\\end{tabular}}\\end{table}\n\n")
# as linhas do LS
T.append("\\begin{table}[H]\\centering\\small\n\\caption{O line sampling, direção a direção: o $\\nu_\\ell$ de cada linha "
         "(onde ela entra na falha) e o CoV de uma linha, $\\mathrm{sd}(p_\\ell)/P_f$. Quanto menor a dispersão de "
         "$\\nu_\\ell$, melhor a direção.}\\label{tab:ls}\n\\resizebox{\\textwidth}{!}{%\n\\begin{tabular}{llrrrrrr}\\toprule\nproblema & direção & "
         "linhas & $\\beta$ & média de $\\nu_\\ell$ & dp de $\\nu_\\ell$ & CoV por linha & análises por linha \\\\\\midrule\n")
for k in ORDEM:
    for d in P[k]["ls"]:
        m = sum(d["nu"]) / d["N"]
        sd = math.sqrt(sum((x - m) ** 2 for x in d["nu"]) / (d["N"] - 1))
        T.append(f"{P[k]['rotulo']} & {d['nome']} & {d['N']} & {num(d['beta'], 3)} & {num(m, 2)} & {num(sd, 2)} & "
                 f"{pct(d['CoVlinha'], 0)} & {num(d['porlinha'], 2)} \\\\\n")
T.append("\\bottomrule\\end{tabular}}\\end{table}\n\n")
# ISMC por semente (Cho)
T.append("\\begin{table}[H]\\centering\\small\n\\caption{A ISMC da etapa 3b, semente a semente (Cho, base de Nyström de "
         "2000 modos, sementes 21 a 25, as mesmas nos três centros). Nenhuma semente chegou ainda ao CoV de 5\\,\\%.}"
         "\\label{tab:ismcsem}\n\\resizebox{\\textwidth}{!}{%\n\\begin{tabular}{ll" + "r" * len(SEM) + "rr}\\toprule\n"
         "caso & centro & " + " & ".join(f"$P_f$ s{s} (CoV)" for s in SEM) + " & média & dp entre sementes \\\\\\midrule\n")
for k in ("r05", "r07"):
    for d in P[k]["ismc"]:
        ps = [s["Pf"] for s in d["sem"]]
        m = sum(ps) / len(ps)
        sd = math.sqrt(sum((x - m) ** 2 for x in ps) / (len(ps) - 1))
        cel = " & ".join(f"{num(100 * s['Pf'], 3)} ({num(100 * s['CoV'], 1)})" for s in d["sem"])
        T.append(f"{P[k]['rotulo']} & {d['nome']} & {cel} & {num(100 * m, 3)}\\,\\% & {num(100 * sd / m, 1)}\\,\\% \\\\\n")
T.append("\\bottomrule\\end{tabular}}\\end{table}\n")
open(os.path.join(AQUI, "tabelas.tex"), "w").write("".join(T))

# numeros para o texto
Nn = []


def cmd(nome, val):
    Nn.append(f"\\newcommand{{\\{nome}}}{{{val}}}\n")


for k, sig in (("ceron", "Ce"), ("r05", "Ra"), ("r07", "Rb")):
    L = linhas_problema(k)
    ls_ = [x for x in L if x["tipo"] == "ls"]
    is_ = [x for x in L if x["tipo"] == "ismc"]
    bl = min(ls_, key=lambda x: x["n5"])
    bi = min(is_, key=lambda x: x["n5"])
    cmd(f"refPf{sig}", pct(P[k]["ref"]["Pf"], 3))
    cmd(f"refN{sig}", num(P[k]["ref"]["n5"]))
    cmd(f"lsN{sig}", num(bl["n5"]))
    cmd(f"lsPf{sig}", pct(bl["Pf"], 3))
    cmd(f"lsGanho{sig}", num(bl["ganho"], 0))
    cmd(f"lsGanhoIsmc{sig}", num(bl["ganho_ismc"], 1))
    cmd(f"ismcN{sig}", num(bi["n5"]))
open(os.path.join(AQUI, "numeros.tex"), "w").write("".join(Nn))

# ---------------- as figuras ----------------
COR = {"ismc": ["#e08214", "#b35806", "#fdb863"], "ls": ["#2166ac", "#4393c3", "#92c5de"]}


def curvas(k):
    pr = P[k]
    out = []
    for i, d in enumerate(pr["ismc"]):
        iw = [x for s in d["sem"] for x in s["iw"]]             # as sementes em sequencia
        pk, cvk = ip.acumulada(iw)
        x = [d["form"]["n"] + j + 1 for j in range(len(iw))]
        out.append(("ISMC, " + d["nome"].replace("$", "").replace("\\,", " "), x, pk, cvk, COR["ismc"][i % 3], "--"))
    for i, d in enumerate(pr["ls"]):
        s1 = s2 = 0.
        x, pk, cvk = [], [], []
        acc = d["form"]["n"] + 2
        for j, (p, na) in enumerate(zip(d["p"], d["na"]), 1):
            s1 += p
            s2 += p * p
            acc += na
            m = s1 / j
            x.append(acc)
            pk.append(m)
            cvk.append(math.sqrt(max((s2 - j * m * m) / (j - 1), 0.) / j) / m if j > 1 and m > 0 else float("inf"))
        out.append(("LS, " + d["nome"].replace("$", "").replace("\\,", " "), x, pk, cvk, COR["ls"][i % 3], "-"))
    return out


fig, ax = plt.subplots(2, 3, figsize=(15, 8.5))
for j, k in enumerate(ORDEM):
    ref = P[k]["ref"]
    a1, a2 = ax[0, j], ax[1, j]
    for nome, x, pk, cvk, c, ls_ in curvas(k):
        a1.semilogx(x, [100 * v for v in pk], ls_, color=c, lw=1.6, label=nome)
        xs = [a for a, b in zip(x, cvk) if math.isfinite(b)]
        a2.loglog(xs, [b for b in cvk if math.isfinite(b)], ls_, color=c, lw=1.6, label=nome)
    lo, hi = ref["Pf"] * (1 - 1.96 * ref["CoV"]), ref["Pf"] * (1 + 1.96 * ref["CoV"])
    a1.axhspan(100 * lo, 100 * hi, color="0.85", zorder=0)
    a1.axhline(100 * ref["Pf"], color="k", lw=1)
    a1.set_ylim(0, 100 * ref["Pf"] * 2.2)
    a1.set_title(P[k]["rotulo"].replace("$", "").replace("{,}", ","))
    a1.set_ylabel("$P_f$ acumulada (%)")
    a2.axhline(0.05, color="k", ls=":", lw=1)
    a2.axvline(ref["n5"], color="0.5", ls=":", lw=1)
    a2.text(ref["n5"], 0.9, " MC direto\n p/ CoV 5%", fontsize=8, color="0.4", va="top")
    a2.set_ylim(0.02, 1.2)
    a2.set_xlabel("análises do MEF (com o FORM)")
    a2.set_ylabel("CoV estimado")
    a1.legend(fontsize=7, loc="upper right")
fig.suptitle("A estimativa (em cima; faixa cinza: IC de 95% do MC direto) e o CoV (embaixo; pontilhado: 5%) contra as "
             "análises do MEF", fontsize=11)
fig.tight_layout()
fig.savefig(os.path.join(FIGS, "convergencia.pdf"))
plt.close(fig)

# o custo para CoV 5%
fig, ax = plt.subplots(1, 4, figsize=(16, 4.6), sharey=True)
grupos = [(P[k]["rotulo"].replace("$", "").replace("{,}", ","), linhas_problema(k)) for k in ORDEM]
grupos.append(("viga (Mathematica)", [{"metodo": n.split(",")[0], "dir": n, "n5": n5,
                                        "tipo": "mc" if n.startswith("MC") else ("ls" if n.startswith("LS") else "ismc")}
                                       for n, fo, g, p, cv, n5 in VIGA["linhas"]]))
for a, (tit, L) in zip(ax, grupos):
    nomes = [(x["metodo"].split(" (")[0] + "\n" + x["dir"].replace("$", "").replace("FORM ", "").replace(
        "\\,", " ").replace("ISMC, centro do ", "").replace("LS, ", "")) if x["tipo"] != "mc" else "MC direto"
             for x in L]
    cores = ["0.4" if x["tipo"] == "mc" else ("#e08214" if x["tipo"] == "ismc" else "#2166ac") for x in L]
    for i, x in enumerate(L):
        a.bar([i], [x["n5"]], color=cores[i], hatch="//" if x.get("parcial") else None,
              edgecolor="k" if x.get("parcial") else cores[i])
        if x.get("parcial"):
            nomes[i] += f"\n(parcial, {x['nlin']} linhas)"
    for i, x in enumerate(L):
        a.text(i, x["n5"] * 1.1, f"{x['n5']:.0f}", ha="center", fontsize=7)
    a.set_yscale("log")
    a.set_xticks(range(len(L)))
    a.set_xticklabels(nomes, rotation=90, fontsize=7)
    a.set_title(tit, fontsize=10)
ax[0].set_ylabel("análises do MEF para CoV de 5% (com o FORM)")
fig.tight_layout()
fig.savefig(os.path.join(FIGS, "custo.pdf"))
plt.close(fig)

# as linhas: o nu de cada linha, por direcao
fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
for a, k in zip(ax, ORDEM):
    for i, d in enumerate(P[k]["ls"]):
        a.hist(d["nu"], bins=20, histtype="step", lw=1.8, color=COR["ls"][i % 3],
               label=d["nome"].replace("$", "") + f" (dp {math.sqrt(sum((x - sum(d['nu']) / d['N']) ** 2 for x in d['nu']) / (d['N'] - 1)):.2f})")
        a.axvline(d["beta"], color=COR["ls"][i % 3], ls=":", lw=1)
    a.set_title(P[k]["rotulo"].replace("$", "").replace("{,}", ","))
    a.set_xlabel(r"$\nu_\ell$ (onde a linha entra na falha); pontilhado: $\beta$ da direção")
    a.set_ylabel("linhas")
    a.legend(fontsize=7)
fig.tight_layout()
fig.savefig(os.path.join(FIGS, "linhas_nu.pdf"))
plt.close(fig)

# ISMC: a P_f de cada semente
fig, ax = plt.subplots(1, 2, figsize=(12, 4))
for a, k in zip(ax, ("r05", "r07")):
    ref = P[k]["ref"]
    for i, d in enumerate(P[k]["ismc"]):
        xs = [i + 0.12 * (j - 2) for j in range(len(d["sem"]))]
        a.errorbar(xs, [100 * s["Pf"] for s in d["sem"]], yerr=[196 * s["Pf"] * s["CoV"] for s in d["sem"]], fmt="o",
                   color=COR["ismc"][i % 3], ms=4, capsize=2)
    for i, d in enumerate(P[k]["ls"]):
        x = len(P[k]["ismc"]) + i
        a.errorbar([x], [100 * d["Pf"]], yerr=[196 * d["Pf"] * d["CoV"]], fmt="s", color=COR["ls"][i % 3], ms=6, capsize=3)
    a.axhspan(100 * ref["Pf"] * (1 - 1.96 * ref["CoV"]), 100 * ref["Pf"] * (1 + 1.96 * ref["CoV"]), color="0.85", zorder=0)
    a.axhline(100 * ref["Pf"], color="k", lw=1)
    lab = [("ISMC\n" + d["nome"]) for d in P[k]["ismc"]] + [("LS\n" + d["nome"]) for d in P[k]["ls"]]
    a.set_xticks(range(len(lab)))
    a.set_xticklabels([s.replace("$", "").replace("FORM ", "") for s in lab], fontsize=8)
    a.set_ylabel("$P_f$ (%) e IC de 95%")
    a.set_title(P[k]["rotulo"].replace("$", "").replace("{,}", ",") + ": ISMC (5 sementes, ~1040 cada) e LS (1 semente)")
fig.tight_layout()
fig.savefig(os.path.join(FIGS, "sementes.pdf"))
plt.close(fig)
print("pronto:", ", ".join(sorted(os.listdir(FIGS))))
for k in ORDEM:
    for x in linhas_problema(k):
        print(f"{k:6s} {x['metodo'][:22]:22s} {x['dir'][:28]:28s} Pf {100 * x['Pf']:.3f}% CoV {100 * x['CoV']:.1f}% "
              f"z {x['z'] if x['z'] is None else round(x['z'], 2)} n5 {x['n5']:.0f} ganho {x['ganho']:.1f}")
