#!/usr/bin/env python3
"""gera_dados.py: le os resultados parciais (Monte Carlo direto, ISMC e line sampling) dos tres problemas de talude e da
viga, e grava as figuras (figs/*.pdf) e as tabelas e numeros do relatorio (tabelas.tex, numeros.tex).

Fontes (no PC, em ~/projects/runs; a viga vem dos relatorios em Mathematica):
  Cho r = -0,5 e -0,7: direto r05 (fermi_etapa1/etapa1/r05_direto, 24 000) e r07 (fermi_etapa3/etapa3/r07_direto,
      110 000); ISMC da etapa 3b (fermi_etapa3b/etapa3b/r0x_nys_<centro>_s21..25, base de Nystrom de 2000 modos); FORM
      da etapa 3b (form_r0x_<M20|M100|gim20>/form.json); line sampling no PC (etapa3b_ls/ls/r0x_<dir>/linhas.csv) com as
      direcoes etapa3b_ls/form_r0x_<gim2000|M20>/form.json.
  Ceron (rebaixamento com a agua antiga, Pf ~16%): direto de 2000 (ceron_transiente/direto/ceron_tr_fase1/f50_I0_iii); ISMC
      (ceron_form/ismc_dif, ismc_modal); line sampling (ceron_form/ls_dif, ls_modal, ls_modal2000) e os FORM
      (ceron_form/form_*/form.json).
  Viga: MC direto e ISMC do relatorio da viga (viga_mathematica/relatorio_viga.tex, "O ponto otimo do M do FORM");
      line sampling e FORM adjunto das secoes 16 e 17 do codigo viga_mathematica/codigo/viga_FORM5_AI200.txt (saida do
      Mathematica de 19/11/2026, ainda fora do relatorio da viga).
  SRM x GIM pareado (srm_gim, no fim): a campanha etapa3b_gim do fermi (fermi_etapa3b_gim/etapa3b_gim: r0x_nys_<centro>_s21
      e ls/r0x_<dir>, PARCIAL enquanto roda) contra as mesmas realizacoes e linhas do SRM; o piloto da analise limite
      (Projects/CeronCinematico/cho/resultados/piloto_h0.025).
  Ceron com a agua do artigo (ceron_agua, no fim): line sampling do MEF (ceron_ls), MC de 10 000 pela analise limite em
      Python (ceron_mc_python/mc10000) e os deterministicos do GIM (ceron_contorno).
Uso: python3 -I gera_dados.py   (roda em segundos; refaz tudo com o que houver de dados)."""
import csv
import json
import math
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.ticker  # noqa: E402,F401

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
    pr = {"rotulo": rotulo, "ref": direto(ref), "ismc": [], "ls": [], "refpasta": ref}
    for cen, nome in (("M20", "FORM $M = 20$, numérico"), ("M100", "FORM $M = 100$, numérico"),
                      ("gim20", "FORM $M = 20$, modal")):
        fo = form(os.path.join(E3B, f"form_{r}_{cen}", "form.json"))
        d = ismc([os.path.join(E3B, f"{r}_nys_{cen}_s{s}") for s in SEM])
        if d:
            pr["ismc"].append({"nome": nome, "form": fo, "cen": cen, "fpasta": os.path.join(E3B, f"form_{r}_{cen}"),
                               "pastas": [os.path.join(E3B, f"{r}_nys_{cen}_s{s}") for s in SEM], **d})
    for dirc, nome in (("M20", "FORM $M = 20$, numérico"), ("gim2000", "FORM $M = 2000$, modal")):
        fo = form(os.path.join(LSC, f"form_{r}_{dirc}", "form.json")) if os.path.exists(
            os.path.join(LSC, f"form_{r}_{dirc}", "form.json")) else None
        d = ls(os.path.join(LSC, "ls", f"{r}_{dirc}"))
        if d and fo:
            # o FORM M = 20 do LS e copia do da etapa 3b (sem os CSV das analises): o custo vem da pasta do fermi
            fp = os.path.join(E3B, f"form_{r}_M20") if dirc == "M20" else os.path.join(LSC, f"form_{r}_{dirc}")
            pr["ls"].append({"nome": nome, "form": fo, "cen": dirc, "fpasta": fp,
                             "pasta": os.path.join(LSC, "ls", f"{r}_{dirc}"), **d})
    P[r] = pr
CER_DIR = os.path.join(R, "ceron_transiente", "direto", "ceron_tr_fase1", "f50_I0_iii")
# o rotulo diz a agua: tudo do Ceron nas tabelas e figuras de antes e do modelo antigo (transiente de 0,1 h, nao saturado,
# contorno do reservatorio, 3 niveis); com a agua do artigo, ceron_agua() no fim
pc = {"rotulo": "Ceron (água antiga)", "ref": direto(CER_DIR), "ismc": [], "ls": [], "refpasta": CER_DIR}
for cen, nome in (("dif", "FORM $M = 20$, numérico"), ("modal", "FORM $M = 20$, modal")):
    fo = form(os.path.join(CER, f"form_{cen}", "form.json"))
    d = ismc([os.path.join(CER, f"ismc_{cen}")])
    pc["ismc"].append({"nome": nome, "form": fo, "cen": cen, "fpasta": os.path.join(CER, f"form_{cen}"),
                       "pastas": [os.path.join(CER, f"ismc_{cen}")], **d})
for dirc, fdir, nome in (("dif", "form_dif", "FORM $M = 20$, numérico"), ("modal", "form_modal", "FORM $M = 20$, modal"),
                         ("modal2000", "form_modal2000", "FORM $M = 2000$, modal")):
    fo = form(os.path.join(CER, fdir, "form.json"))
    d = ls(os.path.join(CER, f"ls_{dirc}"))
    if d:
        pc["ls"].append({"nome": nome, "form": fo, "cen": dirc, "fpasta": os.path.join(CER, fdir),
                         "pasta": os.path.join(CER, f"ls_{dirc}"), **d})
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
# a comparacao, um problema por tabela (numa pagina em paisagem, sem encolher)
T.append("\\begin{landscape}\n")
for k in ORDEM:
    T.append(f"\\begin{{table}}[H]\\centering\\small\n\\caption{{{P[k]['rotulo']}: todos os métodos. ``Gasto'': as "
             f"análises do MEF já feitas (com o FORM). ``Total p/ CoV 5\\,\\%'': o FORM mais as amostras ou linhas para "
             f"CoV de 5\\,\\%, pelo CoV estimado (na ISMC, a média das sementes). $z$: a diferença contra o direto, em "
             f"desvios padrão.}}\\label{{tab:{k}}}\n"
             "\\begin{tabular}{llrrrrrrrr}\\toprule\nmétodo & centro ou direção & FORM & gasto & $P_f$ & CoV & $z$ & "
             "total p/ CoV 5\\,\\% & ganho vs direto & ganho vs ISMC \\\\\\midrule\n")
    for x in linhas_problema(k):
        neg = x["tipo"] == "ls" and x["n5"] == min(y["n5"] for y in linhas_problema(k))
        b = (lambda s: "\\textbf{" + s + "}") if neg else (lambda s: s)
        T.append(f"{x['metodo']} & {x['dir']} & {x['form'] if x['form'] is not None else '--'} & {num(x['gasto'])} & "
                 f"{pct(x['Pf'], 3)} & {pct(x['CoV'], 1)} & {num(x['z'], 2) if x['z'] is not None else '--'} & "
                 f"{b(num(x['n5']))} & {b(num(x['ganho'], 1) + chr(215))} & "
                 f"{(num(x['ganho_ismc'], 1) + chr(215)) if x['ganho_ismc'] else '--'} \\\\\n")
    T.append("\\bottomrule\\end{tabular}\\end{table}\n\n")
T.append("\\end{landscape}\n\n")
# viga
T.append("\\begin{table}[H]\\centering\\small\n\\caption{A viga com campo aleatório (Mathematica). MC direto e ISMC: "
         "relatório da viga, ``O ponto ótimo do M do FORM''; na ISMC, o total usa o CoV verdadeiro, calculado com as falhas "
         "do direto de 50\\,000. Line sampling e FORM adjunto: seções 16 e 17 do código "
         "\\texttt{viga\\_mathematica/codigo/viga\\_FORM5\\_AI200.txt}, ainda fora do relatório da viga. FORM: só as análises do "
         "FORM; ``gasto'' e ``total'' do line sampling incluem as 2 da inclinação.}"
         "\\label{tab:viga}\n\\begin{tabular}{lrrrrrr}\\toprule\nmétodo & FORM & gasto & $P_f$ & CoV & p/ CoV 5\\,\\% & "
         "ganho \\\\\\midrule\n")
for nome, fo, g, p, cv, n5 in VIGA["linhas"]:
    T.append(f"{nome} & {fo} & {num(g)} & {pct(p, 3)} & {pct(cv, 1)} & {num(n5)} & {num(18666 / n5, 1)}{chr(215)} \\\\\n")
T.append("\\bottomrule\\end{tabular}\\end{table}\n\n")
T.append("\\begin{landscape}\n")
# as linhas do LS
T.append("\\begin{table}[H]\\centering\\small\n\\caption{O line sampling, direção a direção: o $\\nu_\\ell$ de cada linha "
         "(onde ela entra na falha) e o CoV de uma linha, $\\mathrm{sd}(p_\\ell)/P_f$. Quanto menor a dispersão de "
         "$\\nu_\\ell$, melhor a direção.}\\label{tab:ls}\n\\begin{tabular}{llrrrrrr}\\toprule\nproblema & direção & "
         "linhas & $\\beta$ & média de $\\nu_\\ell$ & dp de $\\nu_\\ell$ & CoV por linha & análises por linha \\\\\\midrule\n")
for k in ORDEM:
    for d in P[k]["ls"]:
        m = sum(d["nu"]) / d["N"]
        sd = math.sqrt(sum((x - m) ** 2 for x in d["nu"]) / (d["N"] - 1))
        T.append(f"{P[k]['rotulo']} & {d['nome']} & {d['N']} & {num(d['beta'], 3)} & {num(m, 2)} & {num(sd, 2)} & "
                 f"{pct(d['CoVlinha'], 0)} & {num(d['porlinha'], 2)} \\\\\n")
T.append("\\bottomrule\\end{tabular}\\end{table}\n\n")
# ISMC por semente (Cho)
T.append("\\begin{table}[H]\\centering\\small\n\\caption{A ISMC da etapa 3b, semente a semente (Cho, base de Nyström de "
         "2000 modos, sementes 21 a 25, as mesmas nos três centros). Todas pararam com CoV de 5\\,\\%; aqui, com todas as "
         "realizações de cada semente, inclusive as que terminaram depois da parada (na parada: Tabela~\\ref{tab:covobs}).}"
         "\\label{tab:ismcsem}\n\\begin{tabular}{ll" + "r" * len(SEM) + "rr}\\toprule\n"
         "caso & centro & " + " & ".join(f"$P_f$ s{s} (CoV)" for s in SEM) + " & média & dp entre sementes \\\\\\midrule\n")
for k in ("r05", "r07"):
    for d in P[k]["ismc"]:
        ps = [s["Pf"] for s in d["sem"]]
        m = sum(ps) / len(ps)
        sd = math.sqrt(sum((x - m) ** 2 for x in ps) / (len(ps) - 1))
        cel = " & ".join(f"{num(100 * s['Pf'], 3)} ({num(100 * s['CoV'], 1)})" for s in d["sem"])
        T.append(f"{P[k]['rotulo']} & {d['nome']} & {cel} & {num(100 * m, 3)}\\,\\% & {num(100 * sd / m, 1)}\\,\\% \\\\\n")
T.append("\\bottomrule\\end{tabular}\\end{table}\n\\end{landscape}\n")
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
COR = {"ismc": ["#e08214", "#b35806", "#fdb863"], "ls": ["#2166ac", "#4393c3", "#92c5de", "#762a83"]}


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


plt.rcParams.update({"font.size": 10, "axes.titlesize": 11, "axes.labelsize": 10, "legend.fontsize": 8.5,
                     "xtick.labelsize": 9, "ytick.labelsize": 9})
LARG = 16.6 / 2.54   # polegadas: a largura do texto (A4, margens de 2,2 cm), entao as fontes saem no tamanho real


def tit(k):
    return P[k]["rotulo"].replace("$", "").replace("{,}", ",")


# a convergencia: uma figura por problema, a P_f e o CoV lado a lado, a legenda embaixo
for k in ORDEM:
    ref = P[k]["ref"]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(LARG, 3.9))
    for nome, x, pk, cvk, c, ls_ in curvas(k):
        a1.semilogx(x, [100 * v for v in pk], ls_, color=c, lw=1.6, label=nome)
        xs = [a for a, b in zip(x, cvk) if math.isfinite(b)]
        a2.loglog(xs, [b for b in cvk if math.isfinite(b)], ls_, color=c, lw=1.6)
    lo, hi = ref["Pf"] * (1 - 1.96 * ref["CoV"]), ref["Pf"] * (1 + 1.96 * ref["CoV"])
    a1.axhspan(100 * lo, 100 * hi, color="0.85", zorder=0)
    a1.axhline(100 * ref["Pf"], color="k", lw=1)
    a1.set_ylim(0, 100 * ref["Pf"] * 2.2)
    a1.set_ylabel("$P_f$ acumulada (%)")
    a1.set_xlabel("análises do MEF")
    a1.set_title("a estimativa (cinza: IC 95% do direto)", fontsize=10)
    a2.axhline(0.05, color="k", ls=":", lw=1)
    a2.axvline(ref["n5"], color="0.5", ls=":", lw=1)
    a2.text(ref["n5"] * 0.95, 0.9, "MC direto\np/ CoV 5%", fontsize=8, color="0.35", va="top", ha="right")
    a2.set_ylim(0.02, 1.2)
    a2.set_xlabel("análises do MEF")
    a2.set_ylabel("CoV estimado")
    a2.set_title("o CoV (pontilhado: 5%)", fontsize=10)
    h, l = a1.get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=2, frameon=False)
    fig.suptitle(tit(k), fontsize=12)
    fig.tight_layout(rect=(0, 0.03 + 0.045 * math.ceil(len(l) / 2), 1, 1))
    fig.savefig(os.path.join(FIGS, f"convergencia_{k}.pdf"))
    plt.close(fig)

# o custo para CoV 5%: grade 2 x 2
fig, ax = plt.subplots(4, 1, figsize=(LARG, 8.7))
grupos = [(tit(k), linhas_problema(k)) for k in ORDEM]
grupos.append(("viga (Mathematica)", [{"metodo": n.split(",")[0], "dir": n, "n5": n5,
                                        "tipo": "mc" if n.startswith("MC") else ("ls" if n.startswith("LS") else "ismc")}
                                       for n, fo, g, p, cv, n5 in VIGA["linhas"]]))
for a, (t, L) in zip(ax.flat, grupos):
    nomes = [(x["metodo"].split(" (")[0] + ", " + x["dir"].replace("$", "").replace("FORM ", "").replace(
        "\\,", " ").replace("ISMC, centro do ", "").replace("LS, ", "").replace("numérico", "num.")) if x["tipo"] != "mc"
             else "MC direto" for x in L]
    cores = ["0.4" if x["tipo"] == "mc" else ("#e08214" if x["tipo"] == "ismc" else "#2166ac") for x in L]
    for i, x in enumerate(L):
        a.barh([i], [x["n5"]], color=cores[i], hatch="//" if x.get("parcial") else None,
               edgecolor="k" if x.get("parcial") else cores[i])
        if x.get("parcial"):
            nomes[i] += f" (parcial, {x['nlin']} lin.)"
        a.text(x["n5"] * 1.12, i, f"{x['n5']:.0f}", va="center", fontsize=8.5)
    a.set_xscale("log")
    a.set_yticks(range(len(L)))
    a.set_yticklabels(nomes, fontsize=8.5)
    a.invert_yaxis()
    a.set_xlim(50, 6e5)
    a.set_title(t, fontsize=11)
ax[-1].set_xlabel("análises do MEF para CoV de 5% (com o FORM)")
fig.tight_layout()
fig.savefig(os.path.join(FIGS, "custo.pdf"))
plt.close(fig)

# as linhas: o nu de cada linha, por direcao (um problema por linha da figura)
fig, ax = plt.subplots(3, 1, figsize=(LARG, 8.2))
for a, k in zip(ax, ORDEM):
    for i, d in enumerate(P[k]["ls"]):
        m = sum(d["nu"]) / d["N"]
        sd = math.sqrt(sum((x - m) ** 2 for x in d["nu"]) / (d["N"] - 1))
        a.hist(d["nu"], bins=20, histtype="step", lw=1.8, color=COR["ls"][i % 3],
               label=d["nome"].replace("$", "") + f": dp de nu = {sd:.2f} ({d['N']} linhas)")
        a.axvline(d["beta"], color=COR["ls"][i % 3], ls=":", lw=1.2)
    a.set_title(tit(k), fontsize=11)
    a.set_ylabel("linhas")
    a.legend(loc="upper right")
ax[-1].set_xlabel(r"$\nu_\ell$ (onde a linha entra na falha); pontilhado: o $\beta$ da direção")
fig.tight_layout()
fig.savefig(os.path.join(FIGS, "linhas_nu.pdf"))
plt.close(fig)

# ISMC: a P_f de cada semente, e o LS
fig, ax = plt.subplots(2, 1, figsize=(LARG, 6.9))
for a, k in zip(ax, ("r05", "r07")):
    ref = P[k]["ref"]
    for i, d in enumerate(P[k]["ismc"]):
        xs = [i + 0.12 * (j - 2) for j in range(len(d["sem"]))]
        a.errorbar(xs, [100 * s["Pf"] for s in d["sem"]], yerr=[196 * s["Pf"] * s["CoV"] for s in d["sem"]], fmt="o",
                   color=COR["ismc"][i % 3], ms=5, capsize=2)
    for i, d in enumerate(P[k]["ls"]):
        x = len(P[k]["ismc"]) + i
        a.errorbar([x], [100 * d["Pf"]], yerr=[196 * d["Pf"] * d["CoV"]], fmt="s", color=COR["ls"][i % 3], ms=7, capsize=3)
    a.axhspan(100 * ref["Pf"] * (1 - 1.96 * ref["CoV"]), 100 * ref["Pf"] * (1 + 1.96 * ref["CoV"]), color="0.85", zorder=0)
    a.axhline(100 * ref["Pf"], color="k", lw=1)
    lab = [("ISMC\n" + d["nome"]) for d in P[k]["ismc"]] + [("LS" + (" (parcial)" if d["N"] < 100 else "") + "\n" + d["nome"])
                                                            for d in P[k]["ls"]]
    a.set_xticks(range(len(lab)))
    a.set_xticklabels([s.replace("$", "").replace("FORM ", "").replace("numérico", "num.") for s in lab], fontsize=9.5)
    a.set_ylabel("$P_f$ (%) e IC 95%")
    ns = [s["N"] for d in P[k]["ismc"] for s in d["sem"]]
    a.set_title(tit(k) + f": ISMC (5 sementes, {min(ns)} a {max(ns)} cada) e LS (1 semente)", fontsize=10.5)
fig.tight_layout()
fig.savefig(os.path.join(FIGS, "sementes.pdf"))
plt.close(fig)
# =====================================================================================================================
# Diagnosticos: tempo de CPU, iteracoes do FORM, convergencia e eficiencia, ISMC (pesos), line sampling (linhas),
# CoV estimado x observado e a tabela-sintese. Tabelas em tab/*.tex, numeros em numeros.tex (anexados), figuras em figs/.
# Fontes de tempo: s_fator (+ s_transiente no Ceron) de cada realizacao (1 thread por processo: e o tempo de CPU da
# analise); as linhas "fim em <s> s" dos logs dos blocos (fermi); os mtimes dos arquivos (PC; parede).
# =====================================================================================================================
import glob  # noqa: E402
import re  # noqa: E402

TAB = os.path.join(AQUI, "tab")
os.makedirs(TAB, exist_ok=True)
FIM_LS = (-8., 12.)          # os limites da busca do line_sampling.py
GAMMA_ZERO = 0.005
# processos em paralelo do line sampling (rodar.log da etapa3b_ls e rodar_ls.sh do ceron_form)
PROC_LS = {"r05_gim2000": 14, "r05_M20": 14, "r07_gim2000": 12, "r07_M20": 14, "ceron_dif": 14, "ceron_modal": 14,
           "ceron_modal2000": 14}
# quantis da qui-quadrado (95 %, bilateral) para a razao CoV observado / CoV estimado com k - 1 graus de liberdade
QUI2 = {4: (0.4844, 11.143), 9: (2.700, 19.023)}


def le_mc(pasta):
    """as realizacoes real >= 0 (situacao ok) de mc_<a>_<b>.csv, na ordem, ate o primeiro buraco"""
    vis = {}
    for arq in glob.glob(os.path.join(pasta, "mc_*_*.csv")):
        if not re.search(r"mc_-?\d+_-?\d+\.csv$", arq):
            continue
        for r in csv.DictReader(open(arq, newline="")):
            n = int(r["real"])
            if n < 0:
                continue
            fs = float(r["FS_min"])
            vis[n] = {"fs": fs, "lw": float(r.get("log_peso") or 0.), "ok": r["situacao"] == "ok" and fs > 0,
                      "s": float(r.get("s_fator") or 0.) + float(r.get("s_transiente") or 0.)}
    seq, n = [], 0
    while n in vis:
        if vis[n]["ok"]:
            seq.append(vis[n])
        n += 1
    return seq


def est(v):
    n = len(v)
    p = sum(v) / n
    var = sum((x - p) ** 2 for x in v) / (n - 1)
    return p, (math.sqrt(var / n) / p if p > 0 else float("inf"))


def pesos(iw, w):
    t1, t2 = sum(iw), sum(x * x for x in iw)
    s1, s2 = sum(w), sum(x * x for x in w)
    return {"mw": s1 / len(w), "ess_w": s1 * s1 / s2, "ess_iw": t1 * t1 / t2 if t2 > 0 else 0.,
            "maior": max(iw) / t1 if t1 > 0 else float("nan"), "ffal": sum(1 for x in iw if x > 0) / len(iw)}


def semente(pasta):
    seq = le_mc(pasta)
    iw = [math.exp(x["lw"]) * (x["fs"] < 1.) for x in seq]
    w = [math.exp(x["lw"]) for x in seq]
    arq = os.path.join(pasta, "ATINGIU")
    alvo = int(open(arq).read().split()[0]) if os.path.exists(arq) else None
    p, cv = est(iw)
    pa, cva = est(iw[:alvo]) if alvo else (p, cv)
    return {"N": len(iw), "alvo": alvo, "Pf": p, "CoV": cv, "Pa": pa, "CoVa": cva, "iw": iw, "lw": [x["lw"] for x in seq],
            "lw_f": [x["lw"] for x in seq if x["fs"] < 1.], "s": sum(x["s"] for x in seq) / len(seq),
            "h": sum(x["s"] for x in seq) / 3600., **pesos(iw, w)}


def custo_form(pasta):
    """as analises do FORM (it*_0_0.csv, sem os mapas): CPU somada e a parede (do 1o parametros.txt ao form.json)"""
    s, n = 0., 0
    for a in glob.glob(os.path.join(pasta, "it*_0_0.csv")):
        if "_mapa_" in a:
            continue
        for r in csv.DictReader(open(a, newline="")):
            s += float(r.get("s_fator") or 0.) + float(r.get("s_transiente") or 0.)
            n += 1
    par = glob.glob(os.path.join(pasta, "it*_parametros.txt"))
    fj = os.path.join(pasta, "form.json")
    parede = (os.path.getmtime(fj) - min(os.path.getmtime(a) for a in par)) / 3600. if par else None
    d = json.load(open(fj))
    hist = [{"it": h["it"], "g": h["g"], "beta": h["beta"], "grad": math.sqrt(sum(x * x for x in h["grad"]))}
            for h in d["historico"]]
    return {"n": n, "h": s / 3600., "s": s / n if n else float("nan"), "parede": parede, "hist": hist,
            "it": d["iteracoes"], "conv": d["convergiu"], "aval": d["avaliacoes"], "metodo": d["parametros"].get("metodo"),
            "alfa2": {k[6:]: d[k] for k in d if k.startswith("alfa2_")}}


def blocos(pastas_logs):
    """parede somada dos blocos (linhas '<caso> a..b: fim em <s> s' dos logs dos processos), por caso"""
    out = {}
    for pl in pastas_logs:
        for a in glob.glob(os.path.join(pl, "p*.txt")):
            for ln in open(a, errors="replace"):
                m = re.search(r"\d\d:\d\d:\d\d (\S+) -?\d+\.\.-?\d+: fim em (\d+) s", ln)
                if m:
                    out[m.group(1)] = out.get(m.group(1), 0.) + int(m.group(2)) / 3600.
    return out


def ls_diag(pasta, beta):
    """reconstroi a busca de cada linha pelo analises.csv (linha, c, fator): passos de secante e de bissecao, linhas
    nos limites FIM, no maximo de 10 analises, fatores = 0; e as curvas ln(fator) x c"""
    por = {}
    for ln in open(os.path.join(pasta, "analises.csv")):
        if ln.strip():
            li, c, G = ln.split(",")
            por.setdefault(int(li), []).append((float(c), float(G)))
    por.pop(-1, None)
    sec = bis = desvio = nmax = nzero = 0
    for li, seq in por.items():
        g = [(c, math.log(max(G, 1e-300))) for c, G in seq]
        nzero += sum(1 for c, G in seq if G <= GAMMA_ZERO * 1.01)
        if len(seq) >= 10:
            nmax += 1
        for k in range(2, len(g)):
            (ca, ga), (cb, gb) = g[k - 2], g[k - 1]
            cn = cb - gb * (cb - ca) / (gb - ga) if gb != ga else cb + (1. if gb > 0 else -1.)
            cn = min(max(cn, FIM_LS[0]), FIM_LS[1])
            segs = [c for c, x in g[:k] if x > 0]
            fals = [c for c, x in g[:k] if x <= 0]
            if segs and fals and not (min(max(segs), min(fals)) <= cn <= max(max(segs), min(fals))):
                cn = 0.5 * (max(segs) + min(fals))
                bis += 1
            else:
                sec += 1
            if abs(cn - g[k][0]) > 1e-6:
                desvio += 1
    nu = [float(r["c"]) for r in csv.DictReader(open(os.path.join(pasta, "linhas.csv")))]
    nfim = sum(1 for x in nu if x <= FIM_LS[0] or x >= FIM_LS[1])
    ini, fim = os.path.getmtime(os.path.join(pasta, "inclinacao.json")), os.path.getmtime(os.path.join(pasta, "linhas.csv"))
    return {"sec": sec, "bis": bis, "desvio": desvio, "nmax": nmax, "nzero": nzero, "nfim": nfim, "curvas": por,
            "parede": (fim - ini) / 3600.}


def hist_na(na):
    h = {}
    for x in na:
        h[x] = h.get(x, 0) + 1
    return h


# ---------------- leitura ----------------
LOGS = {"r05": [os.path.join(R, "fermi_etapa1", "etapa1", "logs"), os.path.join(E3B, "logs")],
        "r07": [os.path.join(R, "fermi_etapa3", "etapa3", "logs"), os.path.join(E3B, "logs")]}
BLK = {k: blocos(v) for k, v in LOGS.items()}
for k in ORDEM:
    pr = P[k]
    seq = le_mc(pr["refpasta"])
    pr["ref"]["seq_i"] = [1. * (x["fs"] < 1.) for x in seq]
    pr["ref"]["s"] = sum(x["s"] for x in seq) / len(seq)
    pr["ref"]["h"] = sum(x["s"] for x in seq) / 3600.
    pr["ref"]["parede_blocos"] = BLK.get(k, {}).get(os.path.basename(pr["refpasta"]))
    for d in pr["ismc"]:
        d["sd"] = [semente(pa) for pa in d["pastas"]]
        d["fc"] = custo_form(d["fpasta"])
        d["s"] = sum(x["s"] * x["N"] for x in d["sd"]) / sum(x["N"] for x in d["sd"])
        d["h"] = sum(x["h"] for x in d["sd"])
        d["parede_blocos"] = sum(BLK.get(k, {}).get(os.path.basename(pa), 0.) for pa in d["pastas"]) or None
        d["iw_all"] = [x for s in d["sd"] for x in s["iw"]]
    for d in pr["ls"]:
        d["fc"] = custo_form(d["fpasta"])
        d["dg"] = ls_diag(d["pasta"], d["beta"])
        d["proc"] = PROC_LS.get(f"{k}_{d['cen']}")
# o custo unitario: no Cho, o s_fator medio das ISMC do mesmo r (fermi, SRM, base de 2000 modos) vale tambem para o LS
# (que rodou no PC, sem tempo por analise); no Ceron, o direto (8 processos, sem disputa por nucleo) vale para todos
for k in ORDEM:
    pr = P[k]
    if k == "ceron":
        pr["u_amostra"] = pr["ref"]["s"]
        for d in pr["ismc"] + pr["ls"]:
            d["u"] = pr["ref"]["s"]
    else:
        tot = [x for d in pr["ismc"] for x in d["sd"]]
        pr["u_amostra"] = sum(x["s"] * x["N"] for x in tot) / sum(x["N"] for x in tot)
        for d in pr["ismc"]:
            d["u"] = d["s"]
        for d in pr["ls"]:
            d["u"] = pr["u_amostra"]
    pr["ref"]["u"] = pr["ref"]["s"]


def h_form(k, d):
    """o custo do FORM em h-nucleo: o numerico pela CPU somada (no Ceron, pelas analises x o custo unitario comum); o
    modal pela parede (uma analise por iteracao, em serie, e a parede inclui o campo e o gradiente modal)"""
    fc = d["fc"]
    if d["form"]["grad"] == "modal":
        return fc["parede"]
    if k == "ceron":
        return fc["n"] * P[k]["u_amostra"] / 3600.
    return fc["h"]


# a viga (Mathematica): os .mx exportados por viga/extrai_viga.wls
VJ = os.path.join(AQUI, "viga", "viga_dados.json")
VD = json.load(open(VJ)) if os.path.exists(VJ) else None
DLIM = 0.00075
if VD:
    for kk, a in VD["ai"].items():
        iw = [w * (f > DLIM) for f, w in zip(a["flechas"], a["pesos"])]
        a["iw"] = iw
        a["Pf"], a["CoV"] = est(iw)
        a["lw"] = [math.log(w) for w in a["pesos"]]
        a["lw_f"] = [math.log(w) for f, w in zip(a["flechas"], a["pesos"]) if f > DLIM]
        a.update(pesos(iw, a["pesos"]))
    for kk, a in VD["ls"].items():
        a["p"] = [cauda(x) for x in a["c"]]
        a["Pf"], a["CoV"] = est(a["p"])
        a["N"] = len(a["c"])
        a["porlinha"] = sum(a["nLinha"]) / a["N"]
        a["analises"] = sum(a["nLinha"]) + 2
        a["n5"] = 2 + a["porlinha"] * a["N"] * (a["CoV"] / ALVO) ** 2
    VD["ref_i"] = [1. * (f > DLIM) for f in VD["mc_flechas"]]
VFORM = {"f8": "M8", "f20": "M20", "f200": "M200", "adj200": "adj200", "f200_s2": "M200"}
VNOME = {"f8": "FORM $M = 8$, numérico", "f20": "FORM $M = 20$, numérico", "f200": "FORM $M = 200$, numérico",
         "adj200": "FORM $M = 200$, adjunto", "f200_s2": "FORM $M = 200$, numérico"}


def fmt_h(x, d=1):
    return "--" if x is None or (isinstance(x, float) and math.isnan(x)) else num(x, d)


# ---------------- as metricas de eficiencia, por metodo ----------------
def metricas(k):
    """uma linha por metodo: analises do FORM e da amostragem, P_f, CoV, delta_1 = CoV sqrt(n) (o CoV de uma
    analise de amostragem), n5 (com o FORM), h5 (h-nucleo para CoV 5 %), eficiencias relativas ao direto"""
    pr = P[k]
    ref = pr["ref"]
    out = []
    d1d = math.sqrt((1 - ref["Pf"]) / ref["Pf"])
    out.append({"tipo": "mc", "metodo": "MC direto", "dir": "--", "nf": 0, "itf": 0, "n": ref["N"], "Pf": ref["Pf"],
                "CoV": ref["CoV"], "z": None, "d1": d1d, "n5": ref["n5"], "hf": 0., "h5": ref["n5"] * ref["u"] / 3600.,
                "diag": "--", "u": ref["u"]})
    for d in pr["ismc"]:
        n = len(d["iw_all"])
        d1 = d["CoV"] * math.sqrt(n)
        hf = h_form(k, d)
        n5s = d["n5"]       # media das sementes: N (CoV/0,05)^2
        ess = sum(s["ess_iw"] for s in d["sd"]) / sum(s["N"] for s in d["sd"])
        out.append({"tipo": "ismc", "metodo": "ISMC" + (f" ({d['nsem']} sem.)" if d["nsem"] > 1 else ""),
                    "dir": d["nome"], "nf": d["form"]["n"], "itf": d["fc"]["it"], "n": n, "Pf": d["Pf"], "CoV": d["CoV"],
                    "z": z(d["Pf"], d["CoV"], ref), "d1": d1, "n5": n5s + d["form"]["n"], "hf": hf,
                    "h5": hf + n5s * d["u"] / 3600., "diag": f"ESS(Iw)/N = {num(100 * ess, 0)}\\,\\%", "u": d["u"]})
    for d in pr["ls"]:
        n = d["analises"]
        d1 = d["CoV"] * math.sqrt(n)
        hf = h_form(k, d)
        n5s = 2 + d["porlinha"] * d["N"] * (d["CoV"] / ALVO) ** 2
        out.append({"tipo": "ls", "metodo": f"LS ({d['N']} lin.)", "dir": d["nome"], "nf": d["form"]["n"],
                    "itf": d["fc"]["it"], "n": n, "Pf": d["Pf"], "CoV": d["CoV"], "z": z(d["Pf"], d["CoV"], ref),
                    "d1": d1, "n5": n5s + d["form"]["n"], "hf": hf, "h5": hf + n5s * d["u"] / 3600.,
                    "diag": f"{num(d['porlinha'], 2)} an./linha", "u": d["u"]})
    for x in out:
        x["ef_amostra"] = (d1d / x["d1"]) ** 2
        x["ef_n"] = out[0]["n5"] / x["n5"]
        x["ef_h"] = out[0]["h5"] / x["h5"]
    return out


def metricas_viga():
    if not VD:
        return []
    ref = VIGA["ref"]
    d1d = math.sqrt((1 - ref["Pf"]) / ref["Pf"])
    out = [{"tipo": "mc", "metodo": "MC direto", "dir": "--", "nf": 0, "itf": 0, "n": ref["N"], "Pf": ref["Pf"],
            "CoV": ref["CoV"], "z": None, "d1": d1d, "n5": ref["n5"], "diag": "--"}]
    for kk in ("f8", "f20", "f200_s2"):
        a = VD["ai"][kk]
        fo = VD["form"][VFORM[kk]]
        n5s = a["analises"] * (a["CoV"] / ALVO) ** 2
        out.append({"tipo": "ismc", "metodo": "ISMC", "dir": VNOME[kk], "nf": fo["analises"], "itf": 4,
                    "n": a["analises"], "Pf": a["Pf"], "CoV": a["CoV"], "z": z(a["Pf"], a["CoV"], ref),
                    "d1": a["CoV"] * math.sqrt(a["analises"]), "n5": n5s + fo["analises"],
                    "diag": f"ESS(Iw)/N = {num(100 * a['ess_iw'] / a['analises'], 0)}\\,\\%"})
    for kk in ("f8", "f20", "f200", "adj200"):
        a = VD["ls"][kk]
        fo = VD["form"][VFORM[kk]]
        out.append({"tipo": "ls", "metodo": f"LS ({a['N']} lin.)", "dir": VNOME[kk], "nf": fo["analises"], "itf": 4,
                    "n": a["analises"], "Pf": a["Pf"], "CoV": a["CoV"], "z": z(a["Pf"], a["CoV"], ref),
                    "d1": a["CoV"] * math.sqrt(a["analises"]), "n5": a["n5"] + fo["analises"],
                    "diag": f"{num(a['porlinha'], 2)} an./linha"})
    for x in out:
        x["ef_amostra"] = (d1d / x["d1"]) ** 2
        x["ef_n"] = out[0]["n5"] / x["n5"]
        x["ef_h"] = None
        x["h5"] = None
        x["hf"] = None
    return out


MET = {k: metricas(k) for k in ORDEM}
MET["viga"] = metricas_viga()
ROT = {k: P[k]["rotulo"] for k in ORDEM}
ROT["viga"] = "viga (Mathematica)"


def curto(s):
    return s.replace("FORM ", "").replace("numérico", "num.")


# ---------------- tabela: tempo de CPU ----------------
T = ["% gerado por gera_dados.py: nao editar a mao\n\\begin{landscape}\n\\begin{table}[H]\\centering\\footnotesize\n"
     "\\caption{O tempo de CPU (h$\\cdot$núcleo; cada processo roda com 1 \\emph{thread}, então a CPU de uma análise é o "
     "seu tempo de parede). ``s/análise'': a média medida (\\texttt{s\\_fator}, mais \\texttt{s\\_transiente} no Ceron). "
     "FORM: numérico, CPU somada das análises; modal, a parede (em série; inclui o campo e o gradiente modal). "
     "No Ceron (PC), o FORM numérico e as colunas ``p/ CoV 5\\,\\%'' usam as análises $\\times$ o custo unitário do "
     "direto, \\uDirCe{}\\,s (a CPU medida do FORM numérico, \\formCeronCpu{}\\,h, inflada pela disputa de núcleo, está na "
     "Tabela~\\ref{tab:formit}); ``s/análise'' e ``amostragem'' são as medidas. "
     "Amostragem: o que já foi gasto. ``Parede dos blocos'': a soma dos tempos dos blocos nos logs do fermi "
     "(inclui ler a malha e gravar). LS: sem tempo por análise; a estimativa usa o custo unitário da última coluna, e a "
     "parede do PC $\\times$ processos é um teto. ``p/ CoV 5\\,\\%'': o FORM mais as análises para CoV de 5\\,\\% "
     "vezes o custo unitário.}\\label{tab:tempo}\n"
     "\\begin{tabular}{llrrrrrrr}\\toprule\nproblema & método e centro/direção & \\makecell{s/análise\\\\medido} & "
     "\\makecell{FORM\\\\(h)} & \\makecell{amostragem\\\\(h)} & \\makecell{parede dos\\\\blocos (h)} & "
     "\\makecell{PC: parede\\\\$\\times$ proc. (h)} & \\makecell{p/ CoV\\\\5\\,\\% (h)} & \\makecell{custo\\\\unitário (s)} "
     "\\\\\\midrule\n"]
for k in ORDEM:
    pr = P[k]
    ref = pr["ref"]
    T.append(f"{ROT[k]} & MC direto ({num(ref['N'])}) & {num(ref['s'], 1)} & -- & {fmt_h(ref['h'])} & "
             f"{fmt_h(ref.get('parede_blocos'))} & -- & {fmt_h(MET[k][0]['h5'])} & {num(ref['u'], 1)} \\\\\n")
    for d, x in zip(pr["ismc"], [m for m in MET[k] if m["tipo"] == "ismc"]):
        T.append(f" & ISMC, {curto(d['nome'])} & {num(d['s'], 1)} & {fmt_h(x['hf'], 2)} & {fmt_h(d['h'])} & "
                 f"{fmt_h(d.get('parede_blocos'))} & -- & {fmt_h(x['h5'])} & {num(d['u'], 1)} \\\\\n")
    for d, x in zip(pr["ls"], [m for m in MET[k] if m["tipo"] == "ls"]):
        teto = d["dg"]["parede"] * d["proc"] if d.get("proc") else None
        T.append(f" & LS, {curto(d['nome'])} & -- & {fmt_h(x['hf'], 2)} & $\\approx${fmt_h(d['analises'] * d['u'] / 3600.)} & "
                 f"-- & {fmt_h(teto)} ({d['proc']}) & {fmt_h(x['h5'])} & {num(d['u'], 1)} \\\\\n")
    T.append("\\midrule\n" if k != ORDEM[-1] else "")
T.append("\\bottomrule\\end{tabular}\\end{table}\n\\end{landscape}\n")
open(os.path.join(TAB, "tempo.tex"), "w").write("".join(T))

# ---------------- tabela: o FORM ----------------
T = ["% gerado por gera_dados.py\n\\begin{table}[H]\\centering\\footnotesize\\setlength{\\tabcolsep}{3pt}\n\\caption{Os FORM, iteração a iteração "
     "resumidos: iterações, análises, $g$ inicial e final ($\\mathrm{FS} - 1$ no SRM, $\\Gamma - 1$ no GIM), $\\beta$, "
     "a CPU das análises e a parede (do primeiro arquivo de parâmetros ao \\texttt{form.json}). $\\alpha^2$: a parte da "
     "variância linearizada em cada campo.}\\label{tab:formit}\n"
     "\\begin{tabular}{llrrrrrrrl}\\toprule\nproblema & FORM & it. & análises & $g_0$ & $g$ final & $\\beta$ & "
     "\\makecell{CPU\\\\(h)} & \\makecell{parede\\\\(h)} & \\makecell{$\\alpha^2$\\\\($c'$, $\\varphi'$, $k$)} \\\\\\midrule\n"]
FORMS = {}
for k in ORDEM:
    vis = {}
    for d in P[k]["ismc"] + P[k]["ls"]:
        vis.setdefault(d["nome"], d)
    FORMS[k] = vis
    for j, (nome, d) in enumerate(vis.items()):
        fc = d["fc"]
        h = fc["hist"]
        a2 = ", ".join(num(fc["alfa2"][c], 2) for c in ("c", "phi", "k") if c in fc["alfa2"])
        T.append(f"{ROT[k] if j == 0 else ''} & {curto(nome)} ({fc['metodo'].upper()}) & {fc['it']} & {fc['n']} & {num(h[0]['g'], 3)} & "
                 f"{num(h[-1]['g'], 4)} & {num(h[-1]['beta'], 3)} & {fmt_h(fc['h'], 2)} & {fmt_h(fc['parede'], 2)} & "
                 f"{a2} \\\\\n")
if VD:
    for kk in ("M8", "M20", "M200", "adj200"):
        fo = VD["form"][kk]
        T.append(f"{'viga' if kk == 'M8' else ''} & {curto(VNOME['adj200' if kk == 'adj200' else 'f' + kk[1:]])} & 4 & {fo['analises']} & -- & -- & "
                 f"{num(fo['beta'], 3)} & \\multicolumn{{2}}{{c}}{{{num(fo['tempo'], 1)} s (parede)}} & -- \\\\\n")
T.append("\\bottomrule\\end{tabular}\\end{table}\n")
open(os.path.join(TAB, "form_it.tex"), "w").write("".join(T))

# ---------------- tabela: eficiencia ----------------
T = ["% gerado por gera_dados.py\n\\begin{landscape}\n\\begin{table}[H]\\centering\\footnotesize\n\\caption{A eficiência. "
     "$\\delta_1 = \\mathrm{CoV}\\sqrt{n}$: o CoV de uma análise de amostragem ($n$: as análises da amostragem, sem o "
     "FORM); no direto, $\\sqrt{(1 - P_f)/P_f}$. Com $\\mathrm{CoV} = \\delta_1/\\sqrt{n}$, o $n$ para CoV de 5\\,\\% é "
     "$(\\delta_1/0{,}05)^2$. Figura de mérito: $1/(\\mathrm{CoV}^2 \\times \\text{custo})$; a eficiência relativa é a "
     "figura de mérito do método dividida pela do direto: por análise de amostragem, $(\\delta_{1,\\text{dir}}/\\delta_1)^2$; "
     "com o FORM, a razão das análises para CoV de 5\\,\\% (o ``ganho'' das tabelas anteriores); em tempo, a razão das "
     "h$\\cdot$núcleo para CoV de 5\\,\\% (Tabela~\\ref{tab:tempo}).}\\label{tab:efic}\n"
     "\\begin{tabular}{lllrrrrrr}\\toprule\nproblema & método & centro ou direção & $n$ & CoV & $\\delta_1$ & "
     "efic. por análise & efic. com o FORM & efic. em tempo \\\\\\midrule\n"]
for k in ORDEM + ["viga"]:
    for i, x in enumerate(MET[k]):
        T.append(f"{ROT[k] if i == 0 else ''} & {x['metodo']} & {curto(x['dir'])} & {num(x['n'])} & {pct(x['CoV'], 1)} & "
                 f"{num(x['d1'], 2)} & {num(x['ef_amostra'], 1)}{chr(215)} & {num(x['ef_n'], 1)}{chr(215)} & "
                 f"{(num(x['ef_h'], 1) + chr(215)) if x['ef_h'] else '--'} \\\\\n")
    T.append("\\midrule\n" if k != "viga" else "")
T.append("\\bottomrule\\end{tabular}\\end{table}\n\\end{landscape}\n")
open(os.path.join(TAB, "eficiencia.tex"), "w").write("".join(T))

# ---------------- tabela: diagnosticos da ISMC ----------------
T = ["% gerado por gera_dados.py\n\\begin{landscape}\n\\begin{table}[H]\\centering\\footnotesize\n\\caption{Diagnósticos da ISMC, por "
     "centro (no Cho, a média das 5 sementes e, entre parênteses, a menor e a maior). $\\bar w$: a média dos pesos (deve "
     "ficar perto de 1); ESS$(w) = (\\sum w)^2/\\sum w^2$ com todos os pesos e ESS$(Iw)$ só com os que entram na $P_f$; "
     "``maior'': o maior $I w$ sobre a soma; ``falhas'': a fração da amostra deslocada que falha.}\\label{tab:ismcdiag}\n"
     "\\begin{tabular}{llrrrrrr}\\toprule\nproblema & centro & $N$ & $\\bar w$ & ESS$(w)$ & ESS$(Iw)$ & maior & falhas "
     "\\\\\\midrule\n"]


def faixa(v, d=0, p=False):
    m = sum(v) / len(v)
    f = (lambda x: num(100 * x, d)) if p else (lambda x: num(x, d))
    return f(m) + (f" ({f(min(v))}--{f(max(v))})" if len(v) > 1 else "")


for k in ORDEM:
    for d in P[k]["ismc"]:
        sd = d["sd"]
        T.append(f"{ROT[k]} & {curto(d['nome'])} & {faixa([s['N'] for s in sd])} & {faixa([s['mw'] for s in sd], 2)} & "
                 f"{faixa([s['ess_w'] for s in sd])} & {faixa([s['ess_iw'] for s in sd])} & "
                 f"{faixa([s['maior'] for s in sd], 1, True)}\\,\\% & {faixa([s['ffal'] for s in sd], 0, True)}\\,\\% \\\\\n")
if VD:
    for kk in ("f8", "f20", "f200_s2"):
        a = VD["ai"][kk]
        T.append(f"viga & {curto(VNOME[kk])} & {len(a['iw'])} & {num(a['mw'], 2)} & {num(a['ess_w'])} & "
                 f"{num(a['ess_iw'])} & {num(100 * a['maior'], 1)}\\,\\% & {num(100 * a['ffal'])}\\,\\% \\\\\n")
    v = [VD["ai"][f"f200_s{s}"] for s in range(101, 111)]
    T.append(f"viga & {curto(VNOME['f200'])}, 10 sem. & {faixa([len(a['iw']) for a in v])} & {faixa([a['mw'] for a in v], 2)} & "
             f"{faixa([a['ess_w'] for a in v])} & {faixa([a['ess_iw'] for a in v])} & "
             f"{faixa([a['maior'] for a in v], 1, True)}\\,\\% & {faixa([a['ffal'] for a in v], 0, True)}\\,\\% \\\\\n")
T.append("\\bottomrule\\end{tabular}\\end{table}\n\\end{landscape}\n")
open(os.path.join(TAB, "ismc_diag.tex"), "w").write("".join(T))

# ---------------- tabela: diagnosticos do LS ----------------
T = ["% gerado por gera_dados.py\n\\begin{landscape}\n\\begin{table}[H]\\centering\\footnotesize\n\\caption{Diagnósticos "
     "do line sampling. Análises por linha: quantas linhas pararam com 2, 3, 4 e 5 ou mais análises. A busca de cada "
     "linha foi reconstruída do \\texttt{analises.csv} com a regra do \\texttt{line\\_sampling.py}: passos depois do 2.º "
     "(secante ou bissecção), linhas no máximo de 10 análises, nos limites $(-8, 12)$ e análises com $\\Gamma = 0$ "
     "(não houve nenhuma). ``Conferência'': passos em que a reconstrução não bate com o $c$ gravado (0 = a regra foi "
     "reproduzida). $s$: a inclinação $d\\ln(\\text{fator})/dc$ no ponto de projeto. A viga só guardou $\\nu_\\ell$ e "
     "as análises por linha.}\\label{tab:lsdiag}\n"
     "\\begin{tabular}{llrrrrrrrrrrr}\\toprule\nproblema & direção & $s$ & 2 & 3 & 4 & $\\ge 5$ & secante & bissecção & "
     "máx./limites & conferência & $\\nu$ mín.--máx. & CoV por linha \\\\\\midrule\n"]
for k in ORDEM + (["viga"] if VD else []):
    if k == "viga":
        itens = [(VNOME[kk], VD["ls"][kk]["nLinha"], VD["ls"][kk]["c"], VD["ls"][kk]["s"], None,
                  est(VD["ls"][kk]["p"])[1] * math.sqrt(VD["ls"][kk]["N"])) for kk in ("f8", "f20", "f200", "adj200")]
    else:
        itens = [(d["nome"], d["na"], d["nu"], d["s"], d["dg"], d["CoVlinha"]) for d in P[k]["ls"]]
    for nome, na, nu, s, dg, cvl in itens:
        h = hist_na(na)
        cel = " & ".join(str(h.get(j, 0)) for j in (2, 3, 4)) + f" & {sum(v for j, v in h.items() if j >= 5)}"
        dgs = (f"{dg['sec']} & {dg['bis']} & {dg['nmax']}/{dg['nfim']} & {dg['desvio']}" if dg else "-- & -- & -- & --")
        T.append(f"{ROT[k]} & {curto(nome)} & {num(s, 3)} & {cel} & {dgs} & {num(min(nu), 2)}--{num(max(nu), 2)} & "
                 f"{pct(cvl, 0)} \\\\\n")
T.append("\\bottomrule\\end{tabular}\\end{table}\n\\end{landscape}\n")
open(os.path.join(TAB, "ls_diag.tex"), "w").write("".join(T))

# ---------------- CoV estimado x observado e cobertura ----------------
COBS = []   # (rotulo, k sementes, N medio, P media, CoV estimado medio, CoV observado, razao, faixa da razao)
for k in ("r05", "r07"):
    for d in P[k]["ismc"]:
        ps = [s["Pa"] for s in d["sd"]]
        m, _ = est(ps) if len(ps) > 1 else (ps[0], 0)
        sdp = math.sqrt(sum((x - m) ** 2 for x in ps) / (len(ps) - 1))
        cve = sum(s["CoVa"] for s in d["sd"]) / len(ps)
        COBS.append((f"{ROT[k]} & {curto(d['nome'])}", len(ps), sum(s["alvo"] for s in d["sd"]) / len(ps), m, cve,
                     sdp / m, k))
if VD:
    v = [VD["ai"][f"f200_s{s}"] for s in range(101, 111)]
    ps = [a["Pf"] for a in v]
    m = sum(ps) / len(ps)
    sdp = math.sqrt(sum((x - m) ** 2 for x in ps) / (len(ps) - 1))
    COBS.append(("viga & $M = 200$, num.", 10, 1000, m, sum(a["CoV"] for a in v) / 10, sdp / m, "viga"))
T = ["% gerado por gera_dados.py\n\\begin{table}[H]\\centering\\footnotesize\n\\caption{O CoV estimado contra o "
     "observado entre sementes. Cho: as 5 sementes de cada centro, na parada (o $N$ do arquivo \\texttt{ATINGIU}); viga: "
     "10 sementes de 1000 amostras com o centro do FORM $M = 200$. Se o CoV estimado está certo, a razão "
     "observado/estimado cai, com 95\\,\\% de chance, na faixa da última coluna ($\\sqrt{\\chi^2_{k-1}/(k-1)}$).}"
     "\\label{tab:covobs}\n\\begin{tabular}{llrrrrrrr}\\toprule\ncaso & centro da ISMC & sem. & $N$ médio & $P_f$ média & "
     "\\makecell{CoV\\\\estimado} & \\makecell{CoV\\\\observado} & razão & \\makecell{faixa\\\\de 95\\,\\%} \\\\\\midrule\n"]
for rot, ks, nm, m, cve, cvo, _ in COBS:
    lo, hi = QUI2[ks - 1]
    T.append(f"{rot} & {ks} & {num(nm)} & {pct(m, 3)} & {pct(cve, 1)} & {pct(cvo, 1)} & {num(cvo / cve, 2)} & "
             f"{num(math.sqrt(lo / (ks - 1)), 2)}--{num(math.sqrt(hi / (ks - 1)), 2)} \\\\\n")
T.append("\\bottomrule\\end{tabular}\\end{table}\n")
# o vies da parada e a cobertura do IC de 95 %
ZS = []   # (grupo, rotulo, z)
for k in ("r05", "r07"):
    ref = P[k]["ref"]
    for d in P[k]["ismc"]:
        for s in d["sd"]:
            ZS.append((k, "ISMC sem.", z(s["Pa"], s["CoVa"], ref), d["cen"]))
for k in ORDEM:
    for x in MET[k]:
        if x["tipo"] == "ls" or (k == "ceron" and x["tipo"] == "ismc"):
            ZS.append((k, x["tipo"], x["z"], x["dir"]))
if VD:
    for kk, a in VD["ai"].items():
        ZS.append(("viga", "ISMC sem.", z(a["Pf"], a["CoV"], VIGA["ref"]), kk))
    for x in MET["viga"]:
        if x["tipo"] == "ls":
            ZS.append(("viga", "ls", x["z"], x["dir"]))
dentro = sum(1 for g in ZS if abs(g[2]) < 1.96)
vies = [(s["Pa"] - s["Pf"]) / s["Pf"] for k in ("r05", "r07") for d in P[k]["ismc"] for s in d["sd"]]
T.append("\\begin{table}[H]\\centering\\footnotesize\n\\caption{A parada pelo CoV na ISMC do Cho: a $P_f$ na parada "
         "contra a $P_f$ com todas as realizações da mesma semente (os blocos que já rodavam terminaram depois da "
         "parada).}\\label{tab:parada}\n\\begin{tabular}{llrrrr}\\toprule\ncaso & centro & $N$ na parada & $N$ com tudo & "
         "$P_f$ na parada & $P_f$ com tudo \\\\\\midrule\n")
for k in ("r05", "r07"):
    for d in P[k]["ismc"]:
        sd = d["sd"]
        T.append(f"{ROT[k]} & {curto(d['nome'])} & {num(sum(s['alvo'] for s in sd))} & {num(sum(s['N'] for s in sd))} & "
                 f"{pct(sum(s['Pa'] for s in sd) / 5, 3)} & {pct(sum(s['Pf'] for s in sd) / 5, 3)} \\\\\n")
T.append("\\bottomrule\\end{tabular}\\end{table}\n")
open(os.path.join(TAB, "cov_obs.tex"), "w").write("".join(T))

# ---------------- tabela-sintese ----------------
T = ["% gerado por gera_dados.py\n\\begin{landscape}\n\\begin{table}[p]\\centering\\footnotesize\\setlength{\\tabcolsep}{4pt}\n\\caption{Tabela-síntese: "
     "todas as métricas por método e problema. FORM: iterações/análises. $n$: análises da amostragem já feitas. "
     "$\\delta_1$: o CoV de uma análise. $N_5$: análises para CoV de 5\\,\\%, com o FORM; $H_5$: as h$\\cdot$núcleo "
     "correspondentes (Tabela~\\ref{tab:tempo}; a viga só por análises). Efic.: a eficiência relativa ao direto, por "
     "análise de amostragem / com o FORM / em tempo. Diagnóstico: ESS$(Iw)/N$ na ISMC, análises por linha no LS. "
     "Na ISMC da viga, $N_5$ usa o CoV estimado da rodada (na Tabela~\\ref{tab:viga}, o CoV calculado com o direto).}"
     "\\label{tab:sintese}\n\\begin{tabular}{lllrrrrrrrrrl}\\toprule\nproblema & método & centro/direção & FORM & $n$ & "
     "$P_f$ & CoV & $z$ & $\\delta_1$ & $N_5$ & $H_5$ (h) & efic. & diagnóstico \\\\\\midrule\n"]
for k in ORDEM + ["viga"]:
    for i, x in enumerate(MET[k]):
        fo = f"{x['itf']}/{x['nf']}" if x["nf"] else "--"
        ef = f"{num(x['ef_amostra'], 1)}/{num(x['ef_n'], 1)}/{num(x['ef_h'], 1) if x['ef_h'] else '--'}"
        T.append(f"{ROT[k] if i == 0 else ''} & {x['metodo']} & {curto(x['dir'])} & {fo} & {num(x['n'])} & "
                 f"{pct(x['Pf'], 3)} & {pct(x['CoV'], 1)} & {num(x['z'], 2) if x['z'] is not None else '--'} & "
                 f"{num(x['d1'], 2)} & {num(x['n5'])} & {fmt_h(x['h5'])} & {ef} & {x['diag']} \\\\\n")
    T.append("\\midrule\n" if k != "viga" else "")
T.append("\\bottomrule\\end{tabular}\\end{table}\n\\end{landscape}\n")
open(os.path.join(TAB, "sintese.tex"), "w").write("".join(T))

# ---------------- numeros para o texto ----------------
Nn = []
for k, sig in (("ceron", "Ce"), ("r05", "Ra"), ("r07", "Rb")):
    M_ = MET[k]
    ls_ = min([x for x in M_ if x["tipo"] == "ls"], key=lambda x: x["h5"])
    is_ = min([x for x in M_ if x["tipo"] == "ismc"], key=lambda x: x["h5"])
    cmd(f"hDir{sig}", fmt_h(M_[0]["h5"], 0))
    cmd(f"hMedDir{sig}", fmt_h(P[k]["ref"]["h"], 1))
    cmd(f"hLs{sig}", fmt_h(ls_["h5"], 1))
    cmd(f"hIsmc{sig}", fmt_h(is_["h5"], 1))
    cmd(f"efHLs{sig}", num(ls_["ef_h"], 0))
    cmd(f"efALs{sig}", num(max(x["ef_amostra"] for x in M_ if x["tipo"] == "ls"), 0))
    cmd(f"efAIsmc{sig}", num(max(x["ef_amostra"] for x in M_ if x["tipo"] == "ismc"), 0))
    cmd(f"uDir{sig}", num(P[k]["ref"]["s"], 1))
    cmd(f"uAm{sig}", num(P[k]["u_amostra"], 1))
cmd("hIsmcTotal", fmt_h(sum(d["h"] for k in ("r05", "r07") for d in P[k]["ismc"]), 0))
cmd("hIsmcBlocos", fmt_h(sum(d["parede_blocos"] or 0 for k in ("r05", "r07") for d in P[k]["ismc"]), 0))
cmd("hDirTotalCho", fmt_h(P["r05"]["ref"]["h"] + P["r07"]["ref"]["h"], 0))
cmd("zDentro", str(dentro))
cmd("zSd", num(math.sqrt(sum(g[2] ** 2 for g in ZS) / len(ZS)), 2))
# o IC de 95 % da propria semente (sem o erro do direto) contem a P_f do direto?
ICS = [(s["Pa"], s["CoVa"], P[k]["ref"]["Pf"]) for k in ("r05", "r07") for d in P[k]["ismc"] for s in d["sd"]]
if VD:
    ICS += [(a["Pf"], a["CoV"], VIGA["ref"]["Pf"]) for a in VD["ai"].values()]
cmd("icDentro", str(sum(1 for p, cv, r in ICS if abs(p - r) <= 1.96 * p * cv)))
cmd("icTotal", str(len(ICS)))
cmd("zTotal", str(len(ZS)))
cmd("viesParada", num(100 * sum(vies) / len(vies), 2))
cmd("viesParadaMin", num(100 * min(vies), 1))
cmd("viesParadaMax", num(100 * max(vies), 1))
gim = [d["fc"] for k in ("r05", "r07") for d in P[k]["ls"] if d["cen"] == "gim2000"]
cmd("gimParedeIt", num(60 * sum(f["parede"] for f in gim) / sum(f["it"] for f in gim), 1))
cmd("gimCpuIt", num(sum(f["h"] * 3600 for f in gim) / sum(f["n"] for f in gim), 0))
ess_r = [s["ess_iw"] / s["N"] for k in ("r05", "r07") for d in P[k]["ismc"] for s in d["sd"]]
cmd("essIwMin", num(100 * min(ess_r), 0))
cmd("essIwMax", num(100 * max(ess_r), 0))
# o GIM no fermi (FORM modal de 20 termos da etapa 3b): CPU por analise, para comparar com os 26 s do PC
gimf = [d["fc"] for k in ("r05", "r07") for d in P[k]["ismc"] if d["cen"] == "gim20"]
cmd("gimCpuItFermi", num(sum(f["h"] * 3600 for f in gimf) / sum(f["n"] for f in gimf), 0))
# por analise de amostragem, o LS contra a ISMC: (delta_1 ISMC / delta_1 LS)^2, melhor contra melhor e todos os pares (Cho)
for k, sig in (("r05", "Ra"), ("r07", "Rb")):
    d_is = [x["d1"] for x in MET[k] if x["tipo"] == "ismc"]
    d_ls = [x["d1"] for x in MET[k] if x["tipo"] == "ls"]
    cmd(f"lsIsmcMelhor{sig}", num((min(d_is) / min(d_ls)) ** 2, 1))
pares = [(a / b) ** 2 for k in ("r05", "r07") for a in [x["d1"] for x in MET[k] if x["tipo"] == "ismc"]
         for b in [x["d1"] for x in MET[k] if x["tipo"] == "ls"]]
cmd("lsIsmcMin", num(min(pares), 1))
cmd("lsIsmcMax", num(max(pares), 1))
cmd("zMaxTaludes", num(max(abs(x["z"]) for k in ORDEM for x in MET[k] if x["z"] is not None), 2))
cmd("formCeronCpu", fmt_h(P["ceron"]["ismc"][0]["fc"]["h"], 1))
open(os.path.join(AQUI, "numeros.tex"), "a").write("".join(Nn))

# ---------------- figuras ----------------
plt.rcParams.update({"font.size": 10, "axes.titlesize": 11, "axes.labelsize": 10, "legend.fontsize": 9.5,
                     "xtick.labelsize": 10, "ytick.labelsize": 10})
CF = {"dif": "#1b9e77", "M20": "#1b9e77", "M100": "#7570b3", "gim20": "#d95f02", "modal": "#d95f02",
      "modal2000": "#e7298a", "gim2000": "#e7298a"}


# (1) o FORM: beta e |g/g0| por iteracao
fig, ax = plt.subplots(3, 2, figsize=(LARG, 7.8))
for (a1, a2), k in zip(ax, ORDEM):
    for nome, d in FORMS[k].items():
        h = d["fc"]["hist"]
        c = CF.get(d["cen"], "k")
        rot = curto(nome).replace("$", "") + f" ({d['fc']['n']} an.)"
        a1.plot([x["it"] for x in h], [x["beta"] for x in h], "o-", color=c, lw=1.6, ms=4, label=rot)
        a2.semilogy([x["it"] for x in h], [max(abs(x["g"] / h[0]["g"]), 1e-5) for x in h], "o-", color=c, lw=1.6, ms=4)
    a1.set_ylabel(r"$\beta$")
    a2.set_ylabel(r"$|g/g_0|$")
    a1.set_title(tit(k), fontsize=11, loc="left")
    a1.legend(loc="lower right", fontsize=9)
    for a in (a1, a2):
        a.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
for a in ax[-1]:
    a.set_xlabel("iteração do FORM")
fig.tight_layout()
fig.savefig(os.path.join(FIGS, "form_iteracoes.pdf"))
plt.close(fig)


# (2) a eficiencia: CoV x n (log-log, com as retas delta_1/sqrt(n)) e delta_1 corrente, por problema
def correntes(k):
    """(rotulo, cor, estilo, n, CoV_n) para o direto, as ISMC (sementes em sequencia) e o LS, n sem o FORM"""
    out = []
    if k == "viga":
        if not VD:
            return out
        pk, cvk = ip.acumulada(VD["ref_i"])
        out.append(("MC direto", "0.3", "-", list(range(1, len(pk) + 1)), cvk))
        for i, kk in enumerate(("f8", "f20", "f200_s2")):
            pk, cvk = ip.acumulada(VD["ai"][kk]["iw"])
            out.append(("ISMC, " + curto(VNOME[kk]).replace("$", ""), COR["ismc"][i], "--", list(range(1, len(pk) + 1)), cvk))
        seqs = [(kk, VD["ls"][kk]) for kk in ("f8", "f20", "f200", "adj200")]
    else:
        pk, cvk = ip.acumulada(P[k]["ref"]["seq_i"])
        out.append(("MC direto", "0.3", "-", list(range(1, len(pk) + 1)), cvk))
        for i, d in enumerate(P[k]["ismc"]):
            pk, cvk = ip.acumulada(d["iw_all"])
            out.append(("ISMC, " + curto(d["nome"]).replace("$", ""), COR["ismc"][i % 3], "--", list(range(1, len(pk) + 1)),
                        cvk))
        seqs = [(d["cen"], {"c": d["nu"], "nLinha": d["na"], "nome": d["nome"]}) for d in P[k]["ls"]]
    for i, (kk, a) in enumerate(seqs):
        s1 = s2 = 0.
        acc = 2
        xs, cv = [], []
        for j, (c, na) in enumerate(zip(a["c"], a["nLinha"]), 1):
            p = cauda(c)
            s1 += p
            s2 += p * p
            acc += na
            m = s1 / j
            xs.append(acc)
            cv.append(math.sqrt(max((s2 - j * m * m) / (j - 1), 0.) / j) / m if j > 1 and m > 0 else float("inf"))
        nome = a.get("nome") or VNOME[kk]
        out.append(("LS, " + curto(nome).replace("$", ""), COR["ls"][i % 4], "-", xs, cv))
    return out


for k in ORDEM + ["viga"]:
    cs = correntes(k)
    if not cs:
        continue
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(LARG, 3.9))
    for nome, c, st, xs, cv in cs:
        par = [(x, y) for x, y in zip(xs, cv) if math.isfinite(y) and y > 0 and x >= 5]
        if not par:
            continue
        x_, y_ = zip(*par)
        a1.loglog(x_, y_, st, color=c, lw=1.5, label=nome)
        d1 = y_[-1] * math.sqrt(x_[-1])
        a1.loglog([5, 2e5], [d1 / math.sqrt(5), d1 / math.sqrt(2e5)], ":", color=c, lw=0.9)
        a2.loglog(x_, [y * math.sqrt(x) for x, y in par], st, color=c, lw=1.5)
    a1.axhline(0.05, color="k", lw=0.8)
    a1.set_ylim(0.01, 1.5)
    a1.set_xlim(5, 2e5 if k != "ceron" else 5e3)
    a2.set_xlim(5, 2e5 if k != "ceron" else 5e3)
    a1.set_xlabel("análises da amostragem $n$")
    a1.set_ylabel("CoV estimado")
    a1.set_title(r"CoV; pontilhado: $\delta_1/\sqrt{n}$", fontsize=10)
    a2.set_xlabel("análises da amostragem $n$")
    a2.set_ylabel(r"$\delta_1 = \mathrm{CoV}\,\sqrt{n}$")
    fmt_ = matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}".replace(".", ","))
    a2.yaxis.set_major_formatter(fmt_)
    a2.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    a2.yaxis.set_major_locator(matplotlib.ticker.LogLocator(subs=(1., 2., 5.)))
    a2.set_title(r"$\delta_1$ constante: CoV $\propto 1/\sqrt{n}$", fontsize=10)
    h, l = a1.get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=2, frameon=False, fontsize=9.5)
    fig.suptitle(ROT[k].replace("$", "").replace("{,}", ","), fontsize=12)
    fig.tight_layout(rect=(0, 0.01 + 0.05 * math.ceil(len(l) / 2), 1, 1))
    fig.savefig(os.path.join(FIGS, f"eficiencia_{k}.pdf"))
    plt.close(fig)

# (3) eficiencia relativa: barras (por analise de amostragem, com o FORM, em tempo)
fig, ax = plt.subplots(4, 1, figsize=(LARG, 8.8))
for a, k in zip(ax, ORDEM + ["viga"]):
    L = [x for x in MET[k] if x["tipo"] != "mc"]
    nomes = [x["metodo"].split(" (")[0] + ", " + curto(x["dir"]).replace("$", "") for x in L]
    yy = list(range(len(L)))
    a.barh([y - 0.27 for y in yy], [x["ef_amostra"] for x in L], height=0.26, color="#bdbdbd", label="por análise")
    a.barh(yy, [x["ef_n"] for x in L], height=0.26, color="#4393c3", label="com o FORM (análises)")
    if L and L[0]["ef_h"]:
        a.barh([y + 0.27 for y in yy], [x["ef_h"] for x in L], height=0.26, color="#d6604d", label="com o FORM (h·núcleo)")
    a.set_xscale("log")
    a.axvline(1, color="k", lw=0.8)
    a.set_yticks(yy)
    a.set_yticklabels(nomes, fontsize=9.5)
    a.invert_yaxis()
    a.set_title(ROT[k].replace("$", "").replace("{,}", ","), fontsize=11, loc="left")
    a.set_xlim(0.5, 3e3)
ax[0].legend(loc="lower right", fontsize=9)
ax[-1].set_xlabel("eficiência relativa ao MC direto (figura de mérito / a do direto)")
fig.tight_layout()
fig.savefig(os.path.join(FIGS, "eficiencia_barras.pdf"))
plt.close(fig)

# (4) ISMC: log w (todas as amostras e as que falham)
fig, ax = plt.subplots(2, 2, figsize=(LARG, 6.2))
for a, k in zip(ax.flat, ORDEM + ["viga"]):
    if k == "viga":
        if not VD:
            continue
        itens = [(curto(VNOME[kk]).replace("$", ""), VD["ai"][kk]["lw"], VD["ai"][kk]["lw_f"]) for kk in ("f8", "f20", "f200_s2")]
    else:
        itens = [(curto(d["nome"]).replace("$", ""), [x for s in d["sd"] for x in s["lw"]],
                  [x for s in d["sd"] for x in s["lw_f"]]) for d in P[k]["ismc"]]
    for i, (nome, lw, lwf) in enumerate(itens):
        c = COR["ismc"][i % 3]
        a.hist(lw, bins=40, histtype="step", color=c, lw=1.2, ls=":", density=True)
        a.hist(lwf, bins=40, histtype="step", color=c, lw=1.8, density=True, label=nome)
    a.axvline(0, color="k", lw=0.8)
    a.set_title(ROT[k].replace("$", "").replace("{,}", ","), fontsize=11)
    a.set_xlabel(r"$\ln w$")
    a.legend(fontsize=9, loc="best")
ax[0][0].set_ylabel("densidade")
ax[1][0].set_ylabel("densidade")
fig.tight_layout()
fig.savefig(os.path.join(FIGS, "ismc_pesos.pdf"))
plt.close(fig)

# (5) LS: ln(fator) ao longo de cada linha, na melhor direcao de cada talude
fig, ax = plt.subplots(1, 3, figsize=(LARG, 3.3), sharey=False)
for a, k in zip(ax, ORDEM):
    d = min(P[k]["ls"], key=lambda x: x["CoV"])
    for li, seq in d["dg"]["curvas"].items():
        cs = sorted(seq)
        a.plot([c for c, G in cs], [math.log(G) for c, G in cs], "-", color="#4393c3", lw=0.6, alpha=0.6)
        a.plot([c for c, G in cs], [math.log(G) for c, G in cs], ".", color="#2166ac", ms=2.5)
    a.axhline(0, color="k", lw=0.8)
    a.axvline(d["beta"], color="#d6604d", ls=":", lw=1.2)
    a.set_title(tit(k) + "\n" + curto(d["nome"]).replace("$", ""), fontsize=10)
    a.set_xlabel(r"$c$ ao longo da linha")
ax[0].set_ylabel(r"$\ln(\mathrm{fator})$")
fig.tight_layout()
fig.savefig(os.path.join(FIGS, "ls_curvas.pdf"))
plt.close(fig)

# (6) z de todas as estimativas contra o direto
fig, a = plt.subplots(figsize=(LARG, 3.6))
grupos = [("r05", "ISMC sem."), ("r07", "ISMC sem."), ("viga", "ISMC sem."), ("ceron", "ismc"), ("ceron", "ls"),
          ("r05", "ls"), ("r07", "ls"), ("viga", "ls")]
rot_g = []
for i, (k, t) in enumerate(grupos):
    zz = [g[2] for g in ZS if g[0] == k and g[1] == t]
    if not zz:
        continue
    rot_g.append(f"{ {'r05': 'Cho r=-0,5', 'r07': 'Cho r=-0,7', 'ceron': 'Ceron', 'viga': 'viga'}[k]}\n"
                 f"{'ISMC' if t.startswith('ISMC') or t == 'ismc' else 'LS'} ({len(zz)})")
    y = len(rot_g) - 1
    a.plot(zz, [y + 0.15 * ((j % 5) - 2) / 2 for j in range(len(zz))], "o",
           color="#e08214" if "ISMC" in t or t == "ismc" else "#2166ac", ms=5, alpha=0.8)
a.axvspan(-1.96, 1.96, color="0.9", zorder=0)
a.axvline(0, color="k", lw=0.8)
a.set_yticks(range(len(rot_g)))
a.set_yticklabels(rot_g, fontsize=9)
a.invert_yaxis()
a.set_xlabel("$z$ contra o MC direto (cinza: $|z| < 1{,}96$)")
fig.tight_layout()
fig.savefig(os.path.join(FIGS, "z_todos.pdf"))
plt.close(fig)
for k in ORDEM + ["viga"]:
    for x in MET[k]:
        print(f"{k:6s} {x['metodo'][:14]:14s} {x['dir'][:26]:26s} n {x['n']:7d} d1 {x['d1']:.2f} n5 {x['n5']:.0f} "
              f"h5 {x['h5'] if x['h5'] is None else round(x['h5'], 1)} ef {x['ef_amostra']:.1f}/{x['ef_n']:.1f}/"
              f"{x['ef_h'] if x['ef_h'] is None else round(x['ef_h'], 1)}")
print("z dentro:", dentro, "/", len(ZS), "vies parada %:", [round(100 * v, 2) for v in vies])
for c in COBS:
    print(c[0], c[1], round(c[4], 4), round(c[5], 4), round(c[5] / c[4], 2))
print("pronto:", ", ".join(sorted(os.listdir(FIGS))))
for k in ORDEM:
    for x in linhas_problema(k):
        print(f"{k:6s} {x['metodo'][:22]:22s} {x['dir'][:28]:28s} Pf {100 * x['Pf']:.3f}% CoV {100 * x['CoV']:.1f}% "
              f"z {x['z'] if x['z'] is None else round(x['z'], 2)} n5 {x['n5']:.0f} ganho {x['ganho']:.1f}")


# =====================================================================================================================
# Uma semente, pareada: SRM x GIM x analise limite (decisao de 20/11/2026). A ISMC da semente 21 (os 3 centros, r = -0,5
# e -0,7) e o line sampling da semente 7 (as mesmas 100 linhas e direcoes) com o fator pelo GIM (campanha etapa3b_gim do
# fermi, ainda rodando: PARCIAL) contra as mesmas realizacoes e linhas do SRM (etapa 3b e etapa3b_ls); o MC direto do SRM
# e as 5 sementes juntas ficam como referencia. A analise limite (Ceron 2025, Projects/CeronCinematico/cho) roda no
# fermi: aqui so o piloto local (FS_LA x FS_MEF em 200 realizacoes do direto de cada r). E o Ceron com a agua do artigo
# (contorno=ceron agua=sat trebaix=0.0001): o line sampling do MEF e o MC de 10 000 pela analise limite em Python.
# Saidas: tab/srm_gim_ismc.tex, tab/srm_gim_ls.tex, tab/ceron_agua.tex, figs/srm_gim.pdf e os numeros (numeros.tex).
# =====================================================================================================================
import time  # noqa: E402

import matplotlib.patches  # noqa: E402
import pareado_gim as pg  # noqa: E402  (Projects/Cho2010: a leitura pareada dos mc_<a>_<b>.csv e dos linhas.csv)

GIMC = os.path.join(R, "fermi_etapa3b_gim", "etapa3b_gim")
LA_PILOTO = os.path.join(RAIZ, "Projects", "CeronCinematico", "cho", "resultados", "piloto_h0.025")
LA_LOG = os.path.join(R, "fermi_etapa3b_gim", "logs", "laCho_14991.out")   # o MC direto da LA no fermi (job 14991)
RTEX = {"r05": "$-0{,}5$", "r07": "$-0{,}7$", "r00": "$0$"}


def pf_txt(x):
    """P_f em %, com 4 algarismos significativos (1,617 ou 0,3380)"""
    return num(100 * x, 3 if x >= 0.01 else 4)


def sinal(x, d=1):
    """o numero com o sinal ($-$1,7 ou $+$1,5; o zero sem sinal)"""
    s = num(abs(x), d)
    return s if round(abs(x), d) == 0 else ("$-$" if x < 0 else "$+$") + s


def teto(x, d):
    """x arredondado PARA CIMA com d casas: o numero de um limite ("todos tem |FS - 1| <= ...", "no maximo ...")"""
    return num(math.ceil(x * 10 ** d - 1e-9) / 10 ** d, d)


def mdp(v):
    """media e desvio padrao amostral"""
    n = len(v)
    m = sum(v) / n
    return m, (math.sqrt(sum((x - m) ** 2 for x in v) / (n - 1)) if n > 1 else 0.)


def k_pot(pares, perto):
    """o k de Gamma ~ FS^k (minimos quadrados de ln Gamma = k ln FS, como o pareado_gim.py): todas com |FS - 1| > 0,01
    ou so perto de 1 (0,9 < FS < 1,1 e |FS - 1| > 0,005)"""
    xs = [(math.log(f), math.log(g)) for f, g in pares if g > 0 and
          ((0.9 < f < 1.1 and abs(f - 1.) > 0.005) if perto else abs(f - 1.) > 0.01)]
    return sum(x * y for x, y in xs) / sum(x * x for x, _ in xs) if xs else float("nan")


def srm_gim():
    """a comparacao pareada SRM x GIM com uma semente (ISMC s21 e line sampling s7): as duas tabelas, a figura e os
    numeros; a coluna da analise limite (LA) fica "rodando" ate o fermi terminar"""
    ism, lin, datas = [], [], []
    for k in ("r05", "r07"):
        for d5 in P[k]["ismc"]:
            pa = os.path.join(E3B, f"{k}_nys_{d5['cen']}_s21")
            pb = os.path.join(GIMC, f"{k}_nys_{d5['cen']}_s21")
            if not (os.path.isdir(pa) and os.path.isdir(pb)):
                continue
            A, _ = pg.le_mc(pa)
            B, _ = pg.le_mc(pb)
            com = sorted(set(A) & set(B))
            if len(com) < 10:
                continue
            datas += [os.path.getmtime(a) for a in glob.glob(os.path.join(pb, "mc_*_*.csv"))]
            w = [math.exp(A[n]["lw"]) for n in com]
            iS = [x * (A[n]["f"] < 1.) for x, n in zip(w, com)]
            iG = [x * (B[n]["f"] < 1.) for x, n in zip(w, com)]
            pS, cS = est(iS)
            pG, cG = est(iG)
            _, sdd = mdp([g - s for g, s in zip(iG, iS)])      # o erro pareado: so os pares discordantes contam
            pares = [(A[n]["f"], B[n]["f"]) for n in com]
            tS = [A[n]["s"] for n in com if math.isfinite(A[n]["s"])]
            tG = [B[n]["s"] for n in com if math.isfinite(B[n]["s"])]
            ism.append({"k": k, "cen": d5["cen"], "nome": d5["nome"], "N": len(com), "nG": len(B), "pS": pS, "cS": cS,
                        "pG": pG, "cG": cG, "dif": (pG - pS) / pS, "edif": sdd / math.sqrt(len(com)) / pS,
                        "soS": [p for p in pares if p[0] < 1. <= p[1]], "soG": [p for p in pares if p[1] < 1. <= p[0]],
                        # os numeros das realizacoes discordantes: com a mesma semente, a realizacao n dos 3 centros
                        # tem o mesmo z (xi = xi* + z), entao a mesma realizacao pode discordar em mais de um centro
                        "disc_n": [n for n in com if (A[n]["f"] < 1.) != (B[n]["f"] < 1.)],
                        "dlw": max(abs(A[n]["lw"] - B[n]["lw"]) for n in com),
                        "dcampo": max(abs(A[n][c] - B[n][c]) / max(abs(A[n][c]), abs(B[n][c]), 1e-12)
                                      for n in com for c in ("c", "phi")),
                        "k_longe": k_pot(pares, False), "k_perto": k_pot(pares, True), "pares": pares,
                        "tS": sum(tS) / len(tS), "tG": sum(tG) / len(tG), "p5": d5["Pf"], "ref": P[k]["ref"]})
    for k in ("r05", "r07"):
        for dl in P[k]["ls"]:
            pa = os.path.join(LSC, "ls", f"{k}_{dl['cen']}")
            pb = os.path.join(GIMC, "ls", f"{k}_{dl['cen']}")
            if not os.path.exists(os.path.join(pb, "linhas.csv")):
                continue
            A, B = pg.le_ls(pa), pg.le_ls(pb)
            com = sorted(set(A) & set(B))
            if len(com) < 2:
                continue
            datas.append(os.path.getmtime(os.path.join(pb, "linhas.csv")))
            pS, cS = est([cauda(A[i][0]) for i in com])
            pG, cG = est([cauda(B[i][0]) for i in com])
            _, sdd = mdp([cauda(B[i][0]) - cauda(A[i][0]) for i in com])
            dc = [B[i][0] - A[i][0] for i in com]
            mdc, sdc = mdp(dc)
            # a diferenca de P_f que o proprio Delta c preve (linearizada): -sum phi(c) Delta c / sum Phi(-c)
            prev = -sum(math.exp(-0.5 * A[i][0] ** 2) / math.sqrt(2 * math.pi) * (B[i][0] - A[i][0]) for i in com) / sum(
                cauda(A[i][0]) for i in com)
            lin.append({"k": k, "cen": dl["cen"], "nome": dl["nome"], "N": len(com), "nG": len(B), "pS": pS, "cS": cS,
                        "pG": pG, "cG": cG, "dif": (pG - pS) / pS, "edif": sdd / math.sqrt(len(com)) / pS, "prev": prev,
                        "nS": len(A),
                        "dc": dc, "c": [A[i][0] for i in com], "mdc": mdc, "sdc": sdc,
                        "naS": sum(A[i][1] for i in com) / len(com), "naG": sum(B[i][1] for i in com) / len(com),
                        "sS": dl["s"]})
    if not ism or not lin:
        sys.exit("srm_gim: faltam os dados do GIM em " + GIMC)
    data = time.strftime("%d/%m/%Y %H:%M", time.localtime(max(datas)))
    nmin, nmax = min(x["N"] for x in ism), max(x["N"] for x in ism)

    def rot_r(L, i):
        return RTEX[L[i]["k"]] if i == 0 or L[i]["k"] != L[i - 1]["k"] else ""

    # ---- tabela: ISMC
    T = ["% gerado por gera_dados.py (srm_gim): nao editar a mao\n\\begin{table}[H]\\centering\\footnotesize"
         "\\setlength{\\tabcolsep}{3pt}\n"
         f"\\caption{{A ISMC com uma semente (a 21), pareada: o SRM e o GIM nas mesmas realizações (Cho, base de Nyström "
         f"de 2000 modos). \\textbf{{GIM parcial}}: dados de {data}, com $N$ = {nmin} a {nmax} realizações em comum por "
         "caso (a campanha vai a pelo menos 2000). $P_f$ em \\%: a média de $I w$ nas realizações em comum (entre "
         "parênteses, o CoV, "
         "em \\%). Dif.: $(P_\\text{GIM} - P_\\text{SRM})/P_\\text{SRM}$ $\\pm$ o erro pareado, "
         "$\\mathrm{dp}[w\\,(I_\\text{GIM} - I_\\text{SRM})]/(\\sqrt N\\,P_\\text{SRM})$. Discord.: os pares em que só o "
         "SRM falha / só o GIM falha. Referências (SRM): as 5 sementes juntas e o MC direto. LA: a análise limite na "
         "semente 21 (job 14992 do fermi, na fila).}\\label{tab:sgismc}\n"
         "\\begin{tabular}{llrrrrcrrc}\\toprule\n"
         "$r$ & centro & $N$ & \\makecell{SRM s21\\\\$P_f$ (CoV)} & \\makecell{GIM s21\\\\$P_f$ (CoV)} & "
         "\\makecell{dif.\\\\(\\%)} & \\makecell{discord.\\\\SRM/GIM} & \\makecell{SRM\\\\5 sementes} & "
         "\\makecell{MC\\\\direto} & \\makecell{LA\\\\s21} \\\\\\midrule\n"]
    for i, x in enumerate(ism):
        if i and x["k"] != ism[i - 1]["k"]:
            T.append("\\midrule\n")
        T.append(f"{rot_r(ism, i)} & {curto(x['nome'])} & {x['N']} & {pf_txt(x['pS'])} ({num(100 * x['cS'], 1)}) & "
                 f"{pf_txt(x['pG'])} ({num(100 * x['cG'], 1)}) & {sinal(100 * x['dif'])} $\\pm$ {num(100 * x['edif'], 1)} & "
                 f"{len(x['soS'])}/{len(x['soG'])} & {pf_txt(x['p5'])} & {pf_txt(x['ref']['Pf'])} & rodando \\\\\n")
    T.append("\\bottomrule\\end{tabular}\\end{table}\n")
    open(os.path.join(TAB, "srm_gim_ismc.tex"), "w").write("".join(T))
    # ---- tabela: line sampling
    T = ["% gerado por gera_dados.py (srm_gim): nao editar a mao\n\\begin{table}[H]\\centering\\footnotesize"
         "\\setlength{\\tabcolsep}{3pt}\n"
         f"\\caption{{O line sampling com uma semente (a 7), pareado: as mesmas linhas (o mesmo $\\bm z$ de cada uma) e as "
         f"mesmas direções, com o SRM e com o GIM. \\textbf{{GIM parcial}}: dados de {data}; o GIM já tem "
         f"{min(x['nG'] for x in lin)} a {max(x['nG'] for x in lin)} linhas por direção, e entram só as "
         f"{min(x['N'] for x in lin)} em comum com o SRM"
         + (" (todas as linhas do SRM já têm o par no GIM)" if all(x["N"] == x["nS"] for x in lin) else "")
         + ". $P_f$ em \\% (CoV em \\%). "
         "Dif.: como na "
         "Tabela~\\ref{tab:sgismc}, com o erro pareado das $\\Phi(-c_\\ell)$. "
         "$\\Delta c = c_\\text{GIM} - c_\\text{SRM}$: onde a linha entra na falha, GIM menos SRM (média e dp das linhas). "
         "LA: o line sampling pela análise limite nas mesmas linhas (job 14992 do fermi, na fila).}\\label{tab:sgls}\n"
         "\\begin{tabular}{llrrrrrcc}\\toprule\n"
         "$r$ & direção & linhas & \\makecell{SRM\\\\$P_f$ (CoV)} & \\makecell{GIM\\\\$P_f$ (CoV)} & "
         "\\makecell{dif.\\\\(\\%)} & \\makecell{$\\Delta c$\\\\média (dp)} & "
         "\\makecell{análises por linha\\\\SRM / GIM} & LA \\\\\\midrule\n"]
    for i, x in enumerate(lin):
        if i and x["k"] != lin[i - 1]["k"]:
            T.append("\\midrule\n")
        T.append(f"{rot_r(lin, i)} & {curto(x['nome'])} & {x['N']} & {pf_txt(x['pS'])} ({num(100 * x['cS'], 1)}) & "
                 f"{pf_txt(x['pG'])} ({num(100 * x['cG'], 1)}) & {sinal(100 * x['dif'])} $\\pm$ {num(100 * x['edif'], 1)} & "
                 f"{sinal(x['mdc'], 3)} ({num(x['sdc'], 3)}) & {num(x['naS'], 2)} / {num(x['naG'], 2)} & rodando \\\\\n")
    T.append("\\bottomrule\\end{tabular}\\end{table}\n")
    open(os.path.join(TAB, "srm_gim_ls.tex"), "w").write("".join(T))

    # ---- numeros para o texto
    Nn.clear()
    # a conferencia de versao da campanha: as 10 primeiras realizacoes de cada direto, refeitas com o SRM (o mesmo FS?)
    iguais = total = 0
    for k in ("r05", "r07"):
        a, b = os.path.join(GIMC, f"conf_srm_{k}", "mc_0_9.csv"), os.path.join(P[k]["refpasta"], "mc_0_9.csv")
        if os.path.exists(a) and os.path.exists(b):
            fa = {int(r["real"]): float(r["FS_min"]) for r in csv.DictReader(open(a, newline="")) if int(r["real"]) >= 0}
            fb = {int(r["real"]): float(r["FS_min"]) for r in csv.DictReader(open(b, newline="")) if int(r["real"]) >= 0}
            total += len(set(fa) & set(fb))
            iguais += sum(1 for n in set(fa) & set(fb) if abs(fa[n] - fb[n]) < 5e-7)
    cmd("sgConfIguais", str(iguais))
    cmd("sgConfTotal", str(total))
    dd = [p for x in ism for p in x["soS"] + x["soG"]]
    todos = [p for x in ism for p in x["pares"]]
    kp = k_pot(todos, True)
    ks = [x[c] for x in ism for c in ("k_longe", "k_perto")]
    tS, tG = [x["tS"] for x in ism], [x["tG"] for x in ism]
    cmd("sgData", data)
    cmd("sgNmin", str(nmin))
    cmd("sgNmax", str(nmax))
    cmd("sgPares", num(sum(x["N"] for x in ism)))
    cmd("sgDisc", str(len(dd)))
    cmd("sgDiscSrm", str(sum(len(x["soS"]) for x in ism)))
    cmd("sgDiscGim", str(sum(len(x["soG"]) for x in ism)))
    # as realizacoes distintas (a mesma n do mesmo r e a mesma amostra z nos 3 centros) e os discordantes por caso
    cmd("sgDiscDist", str(len({(x["k"], n) for x in ism for n in x["disc_n"]})))
    cmd("sgDiscCasoMin", str(min(len(x["disc_n"]) for x in ism)))
    cmd("sgDiscCasoMax", str(max(len(x["disc_n"]) for x in ism)))
    # limites ("todos tem |FS - 1| <= ..."): arredondados para cima
    cmd("sgDiscFS", teto(max(abs(f - 1.) for f, g in dd), 4) if dd else "--")
    cmd("sgDiscG", teto(max(abs(g - 1.) for f, g in dd), 3) if dd else "--")
    # o |Gamma - 1| dos discordantes em FS equivalente, Gamma^(1/k) - 1
    cmd("sgDiscGfs", teto(100 * max(abs(g ** (1. / kp) - 1.) for f, g in dd if g > 0), 1) if dd else "--")
    cmd("sgKmin", num(min(ks), 1))
    cmd("sgKmax", num(max(ks), 1))
    cmd("sgKperto", num(kp, 2))
    cmd("sgDifMin", sinal(100 * min(x["dif"] for x in ism)))
    cmd("sgDifMax", sinal(100 * max(x["dif"] for x in ism)))
    cmd("sgZmax", num(max([abs(x["dif"]) / x["edif"] for x in ism if x["edif"] > 0] or [0.]), 1))
    cmd("sgEdifMin", num(100 * min([x["edif"] for x in ism if x["edif"] > 0] or [0.]), 1))
    cmd("sgEdifMax", num(100 * max(x["edif"] for x in ism), 1))
    cmd("sgCovMin", num(100 * min(x["cS"] for x in ism), 1))
    cmd("sgCovMax", num(100 * max(x["cS"] for x in ism), 1))
    cmd("sgAbaixo", str(sum(1 for x in ism if x["dif"] < 0)))
    cmd("sgCasos", str(len(ism)))
    cmd("sgCampo", teto(100 * max(x["dcampo"] for x in ism), 2))
    # a maior diferenca relativa das P_f, ISMC e line sampling
    cmd("sgDifAbsMax", num(100 * max(abs(x["dif"]) for x in ism + lin), 1))
    # a P_f da semente 21 nas realizacoes em comum (as primeiras da semente), SRM e GIM, por r
    for k, sig in (("r05", "Ra"), ("r07", "Rb")):
        sub = [p for x in ism if x["k"] == k for p in (x["pS"], x["pG"])]
        cmd(f"sgSub{sig}Min", pf_txt(min(sub)) if sub else "--")
        cmd(f"sgSub{sig}Max", pf_txt(max(sub)) if sub else "--")
    cmd("sgTgimMin", num(min(tG), 0))
    cmd("sgTgimMax", num(max(tG), 0))
    cmd("sgTsrmMin", num(min(tS), 0))
    cmd("sgTsrmMax", num(max(tS), 0))
    cmd("sgEcon", num(100 * (1. - sum(tG) / sum(tS)), 0))
    cmd("sglN", str(min(x["N"] for x in lin)))
    cmd("sglNGmin", str(min(x["nG"] for x in lin)))
    cmd("sglNGmax", str(max(x["nG"] for x in lin)))
    cmd("sglDifMin", sinal(100 * min(x["dif"] for x in lin)))
    cmd("sglDifMax", sinal(100 * max(x["dif"] for x in lin)))
    cmd("sglPrevMin", sinal(100 * min(x["prev"] for x in lin)))
    cmd("sglPrevMax", sinal(100 * max(x["prev"] for x in lin)))
    cmd("sglEdifMin", num(100 * min(x["edif"] for x in lin), 1))
    cmd("sglEdifMax", num(100 * max(x["edif"] for x in lin), 1))
    cmd("sglZmin", num(min(abs(x["dif"]) / x["edif"] for x in lin), 0))
    cmd("sglZmax", num(max(abs(x["dif"]) / x["edif"] for x in lin), 0))
    cmd("sglCovMin", num(100 * min(x["cS"] for x in lin), 1))
    cmd("sglCovMax", num(100 * max(x["cS"] for x in lin), 1))
    cmd("sglDcMin", sinal(min(x["mdc"] for x in lin), 3))
    cmd("sglDcMax", sinal(max(x["mdc"] for x in lin), 3))
    cmd("sglDcSdMin", num(min(x["sdc"] for x in lin), 3))
    cmd("sglDcSdMax", num(max(x["sdc"] for x in lin), 3))
    cmd("sglSmin", num(min(abs(x["sS"]) for x in lin), 2))
    cmd("sglSmax", num(max(abs(x["sS"]) for x in lin), 2))
    cmd("sglDfs", num(100 * max(abs(x["sS"] * x["mdc"]) for x in lin), 2))   # o Delta c em ln FS: |s| Delta c
    cmd("sglNaSrmMin", num(min(x["naS"] for x in lin), 1))
    cmd("sglNaSrmMax", num(max(x["naS"] for x in lin), 1))
    cmd("sglNaGimMin", num(min(x["naG"] for x in lin), 1))
    cmd("sglNaGimMax", num(max(x["naG"] for x in lin), 1))
    # o custo por linha (analises por linha x o custo por analise da ISMC do mesmo r, nas mesmas maquinas do fermi)
    tr = {k: (sum(x["tS"] for x in ism if x["k"] == k), sum(x["tG"] for x in ism if x["k"] == k)) for k in ("r05", "r07")}
    cpl = [x["naG"] * tr[x["k"]][1] / (x["naS"] * tr[x["k"]][0]) - 1. for x in lin if tr[x["k"]][0] > 0]
    cmd("sglCustoMin", sinal(100 * min(cpl), 0))
    cmd("sglCustoMax", sinal(100 * max(cpl), 0))
    # o FORM modal (GIM, 2000 modos por campo) e o numerico M = 20 (SRM) que dao as direcoes do line sampling
    for k, sig in (("r05", "Ra"), ("r07", "Rb")):
        fg = next((d["form"] for d in P[k]["ls"] if d["cen"] == "gim2000"), None)
        fn = next((d["form"] for d in P[k]["ls"] if d["cen"] == "M20"), None)
        cmd(f"fgN{sig}", str(fg["n"]) if fg else "--")
        cmd(f"fgBeta{sig}", num(fg["beta"], 3) if fg else "--")
        cmd(f"fnN{sig}", str(fn["n"]) if fn else "--")
        cmd(f"fnBeta{sig}", num(fn["beta"], 3) if fn else "--")
        # a semente 21 do SRM entre as 5 (todas as realizacoes de cada semente, como na Tabela tab:ismcsem)
        ds = P[k]["ismc"]
        p21 = [d["sem"][0]["Pf"] for d in ds]
        outras = [s["Pf"] for d in ds for s in d["sem"][1:]]
        cmd(f"sUm{sig}Min", pf_txt(min(p21)))
        cmd(f"sUm{sig}Max", pf_txt(max(p21)))
        cmd(f"sOut{sig}Min", pf_txt(min(outras)))
        cmd(f"sOut{sig}Max", pf_txt(max(outras)))
        cmd(f"sJun{sig}Min", pf_txt(min(d["Pf"] for d in ds)))
        cmd(f"sJun{sig}Max", pf_txt(max(d["Pf"] for d in ds)))
        # a posicao da semente 21 entre as 5 de cada centro (1 = a mais baixa)
        cmd(f"sPos{sig}", ", ".join(str(sorted(s["Pf"] for s in d["sem"]).index(d["sem"][0]["Pf"]) + 1) for d in ds))
    # o piloto da analise limite: FS_LA x FS_MEF nas mesmas 200 realizacoes do direto de cada r
    la = []
    for r in ("r05", "r07"):
        arq = os.path.join(LA_PILOTO, f"direto_{r}_pares.csv")
        if not os.path.exists(arq):
            continue
        rs = [(float(x["FS_LA"]), float(x["FS_MEF"])) for x in csv.DictReader(open(arq, newline=""))]
        ma, mb = mdp([a for a, b in rs])[0], mdp([b for a, b in rs])[0]
        corr = sum((a - ma) * (b - mb) for a, b in rs) / math.sqrt(sum((a - ma) ** 2 for a, b in rs) *
                                                                     sum((b - mb) ** 2 for a, b in rs))
        la.append({"n": len(rs), "dif": sum(a - b for a, b in rs) / len(rs), "rel": sum(a / b - 1. for a, b in rs) / len(rs),
                   "corr": corr, "acima": sum(1 for a, b in rs if a > b),
                   "rmin": min(a / b - 1. for a, b in rs), "rmax": max(a / b - 1. for a, b in rs),
                   "fla": sum(1 for a, b in rs if a < 1.), "fmef": sum(1 for a, b in rs if b < 1.)})
    cmd("laN", str(min(x["n"] for x in la)) if la else "--")
    cmd("laTotal", str(sum(x["n"] for x in la)) if la else "--")
    cmd("laAcima", str(sum(x["acima"] for x in la)) if la else "--")
    cmd("laDif", sinal(sum(x["dif"] * x["n"] for x in la) / sum(x["n"] for x in la), 3) if la else "--")
    cmd("laRelMin", num(100 * min(x["rel"] for x in la), 1) if la else "--")
    cmd("laRelMax", num(100 * max(x["rel"] for x in la), 1) if la else "--")
    cmd("laCorr", num(min(x["corr"] for x in la), 3) if la else "--")
    # a faixa por realizacao (FS_LA / FS_MEF - 1) e as falhas (FS < 1) de cada um no piloto
    cmd("laRelLo", num(100 * min(x["rmin"] for x in la), 1) if la else "--")
    cmd("laRelHi", num(100 * max(x["rmax"] for x in la), 1) if la else "--")
    cmd("laFalLa", str(sum(x["fla"] for x in la)) if la else "--")
    cmd("laFalMef", str(sum(x["fmef"] for x in la)) if la else "--")
    # o andamento do MC direto da LA no fermi (o ultimo "N = " do log baixado)
    hora, nla, rla = "--", "--", "--"
    if os.path.exists(LA_LOG):
        for ln in open(LA_LOG, errors="replace"):
            m = re.search(r"===== direto (r\d\d)", ln)
            if m:
                rla = RTEX.get(m.group(1), m.group(1))
            m = re.search(r"(\d\d:\d\d):\d\d N = (\d+),", ln)
            if m:
                hora, nla = m.group(1), num(int(m.group(2)))
    cmd("laLogHora", hora)
    cmd("laLogN", nla)
    cmd("laLogR", rla)
    open(os.path.join(AQUI, "numeros.tex"), "a").write("".join(Nn))

    # ---- figura: Gamma x FS (ISMC) e Delta c x c (line sampling)
    fig, ax = plt.subplots(2, 2, figsize=(LARG, 6.2))
    (a1, a2), (a3, a4) = ax
    n_s, n_g = sum(len(x["soS"]) for x in ism), sum(len(x["soG"]) for x in ism)
    for a, (lo, hi), ms in ((a1, (0.9, 1.1), 2.0), (a2, (0.99, 1.01), 3.2)):
        xs = [p for p in todos if lo < p[0] < hi]
        a.plot([f for f, g in xs], [g for f, g in xs], "o", color="#b35806", ms=ms, alpha=0.45, mew=0,
               label=f"pares ({len(xs)})")
        fx = [lo + (hi - lo) * j / 200. for j in range(201)]
        a.plot(fx, [f ** kp for f in fx], "-", color="k", lw=1.0,
               label=r"$\Gamma = \mathrm{FS}^{" + f"{kp:.2f}".replace(".", "{,}") + "}$")
        a.plot([f for f, g in dd], [g for f, g in dd], "o", mfc="none", mec="k", ms=7, mew=1.1,
               label=f"discordantes ({len(dd)})")
        a.axvline(1, color="0.5", ls=":", lw=0.9)
        a.axhline(1, color="0.5", ls=":", lw=0.9)
        a.set_xlim(lo, hi)
        a.set_xlabel("FS do SRM")
        a.set_ylabel(r"$\Gamma$ do GIM")
    a2.set_ylim(0.97, 1.03)
    for xy in ((0.99, 1.0), (1.0, 0.97)):        # os quadrantes discordantes
        a2.add_patch(matplotlib.patches.Rectangle(xy, 0.01, 0.03, color="0.92", zorder=0, lw=0))
    a2.text(0.9905, 1.0285, f"só o SRM falha ({n_s})", fontsize=8.5, va="top", ha="left", color="0.25")
    a2.text(1.0095, 0.9715, f"só o GIM falha ({n_g})", fontsize=8.5, va="bottom", ha="right", color="0.25")
    a1.legend(loc="upper left", fontsize=8.5)
    a1.set_title(f"(a) ISMC, semente 21: os {len(ism)} casos", fontsize=10, loc="left")
    a2.set_title("(b) o mesmo, perto de (1, 1)", fontsize=10, loc="left")
    for a in (a1, a2):
        a.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}".replace(".", ",")))
        a.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}".replace(".", ",")))
    a2.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.005))
    cor_d = {"M20": "#2166ac", "gim2000": "#d6604d"}
    mar_d = {"M20": "o", "gim2000": "^"}
    for a, k, t in ((a3, "r05", "(c)"), (a4, "r07", "(d)")):
        for x in [y for y in lin if y["k"] == k]:
            a.plot(x["c"], x["dc"], mar_d.get(x["cen"], "s"), color=cor_d.get(x["cen"], "0.3"), ms=3.6, alpha=0.8, mew=0,
                   label=curto(x["nome"]).replace("$", "") + " (média " + f"{x['mdc']:+.4f}".replace(".", ",") + ")")
            a.axhline(x["mdc"], color=cor_d.get(x["cen"], "0.3"), ls="--", lw=1.0)
        a.axhline(0, color="k", lw=0.8)
        a.set_xlabel(r"$c_\mathrm{SRM}$")
        a.set_ylabel(r"$c_\mathrm{GIM} - c_\mathrm{SRM}$")
        a.set_title(f"{t} LS, semente 7: " + tit(k), fontsize=10, loc="left")
        a.legend(loc="lower left", fontsize=8)
        a.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}".replace(".", ",")))
        a.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}".replace(".", ",")))
    yl = max(max(abs(v) for v in x["dc"]) for x in lin) * 1.6
    for a in (a3, a4):
        a.set_ylim(-yl, yl)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGS, "srm_gim.pdf"))
    plt.close(fig)
    for x in ism:
        print(f"SRM x GIM ISMC {x['k']} {x['cen']:6s} N {x['N']} SRM {100 * x['pS']:.4f}% GIM {100 * x['pG']:.4f}% dif "
              f"{100 * x['dif']:+.2f} +- {100 * x['edif']:.2f}% so SRM {len(x['soS'])} so GIM {len(x['soG'])} k "
              f"{x['k_longe']:.2f}/{x['k_perto']:.2f} t {x['tS']:.1f}/{x['tG']:.1f} s dlw {x['dlw']:.1e} campo "
              f"{100 * x['dcampo']:.2f}%")
    for x in lin:
        print(f"SRM x GIM LS {x['k']} {x['cen']:7s} {x['N']} linhas SRM {100 * x['pS']:.4f}% GIM {100 * x['pG']:.4f}% dif "
              f"{100 * x['dif']:+.2f} +- {100 * x['edif']:.2f}% (pelo dc {100 * x['prev']:+.2f}%) dc {x['mdc']:+.4f} "
              f"({x['sdc']:.4f}) an/linha {x['naS']:.2f}/{x['naG']:.2f}")


# ---------------- o Ceron com a agua do artigo ----------------
CERON_ARTIGO = {"Pf": 0.1150, "CoV": 0.0277, "N": 10000, "Gdet": 1.336, "Gmed": 1.353}   # Ceron (2025), Secao 6
CA_LS = os.path.join(R, "ceron_ls")                    # LS do MEF (GIM) com a agua do artigo, FORM modal M = 2000
CA_MC = os.path.join(R, "ceron_mc_python", "mc10000")  # MC de 10 000 pela analise limite em Python
CA_DET = os.path.join(R, "ceron_contorno")             # deterministicos do GIM: sat_ceron_n5 (artigo), padrao (antiga)


def niveis(arq):
    """o Gamma de cada nivel de refinamento (a coluna niveis da realizacao -1, a deterministica)"""
    return [float(x) for x in next(csv.DictReader(open(arq, newline="")))["niveis"].split(";")]


def aitken(a, b, c):
    """a extrapolacao de tres niveis com convergencia geometrica"""
    q = (c - b) / (b - a)
    return c + (c - b) * q / (1. - q)


def ceron_agua():
    """o rebaixamento do Ceron com a agua do artigo: tab/ceron_agua.tex e os numeros"""
    art = CERON_ARTIGO
    d = ls(os.path.join(CA_LS, "ls"))
    fo = json.load(open(os.path.join(CA_LS, "form", "form.json")))
    G = {int(r["real"]): float(r["Gamma"]) for r in csv.DictReader(open(os.path.join(CA_MC, "mc.csv"), newline=""))}
    ex = {int(r["real"]): (float(r["G_rapido"]), float(r["G_exato"]))
          for r in csv.DictReader(open(os.path.join(CA_MC, "exato.csv"), newline=""))}
    Ge = dict(G)
    Ge.update({n: b for n, (a, b) in ex.items()})          # a busca exata no lugar da rapida onde foi refeita
    nmc = len(G)
    p_r = sum(1 for g in G.values() if g < 1.) / nmc
    p_e = sum(1 for g in Ge.values() if g < 1.) / nmc
    novas = sum(1 for a, b in ex.values() if b < 1. <= a)
    sel = [n for n, g in G.items() if 1. <= g < 1.1]
    if any(n not in ex for n in sel):
        print("AVISO ceron_agua: o exato.csv nao cobre todas as realizacoes com 1 <= Gamma < 1,1")
    gla = float(re.search(r"deterministico: Gamma = ([\d.]+)", open(os.path.join(CA_MC, "rodar.log")).read()).group(1))
    nv = niveis(os.path.join(CA_DET, "sat_ceron_n5_0_-1.csv"))
    gpz = aitken(*nv[-3:])
    gvelho = niveis(os.path.join(CA_DET, "padrao_0_-1.csv"))[-1]
    ref = P["ceron"]["ref"]
    fs_velho = [x["fs"] for x in le_mc(P["ceron"]["refpasta"])]

    def cv(p, n):
        return math.sqrt((1. - p) / (n * p))

    def zart(p, c):
        return (p - art["Pf"]) / math.sqrt((p * c) ** 2 + (art["Pf"] * art["CoV"]) ** 2)

    linhas = [("Ceron (2025): análise limite, MC de " + num(art["N"]), art["Pf"], art["CoV"], None,
               num(art["Gdet"], 3), num(art["Gmed"], 3)),
              (f"MEF (GIM, {len(nv) - 1} níveis): LS, {d['N']} linhas e {fo['avaliacoes']} + {d['analises']} análises",
               d["Pf"], d["CoV"], zart(d["Pf"], d["CoV"]), num(gpz, 3) + "$^*$", "--"),
              (f"análise limite (Python): MC de {num(nmc)}, busca rápida", p_r, cv(p_r, nmc), zart(p_r, cv(p_r, nmc)),
               num(gla, 3), num(sum(G.values()) / nmc, 3)),
              (f"\\quad e a busca exata nas {len(ex)} com $1 \\le \\Gamma < 1{{,}}1$", p_e, cv(p_e, nmc),
               zart(p_e, cv(p_e, nmc)), "", num(sum(Ge.values()) / nmc, 3))]
    T = ["% gerado por gera_dados.py (ceron_agua): nao editar a mao\n\\begin{table}[H]\\centering\\footnotesize"
         "\\setlength{\\tabcolsep}{4pt}\n"
         "\\caption{O rebaixamento do Ceron com a água do artigo (\\texttt{agua=sat trebaix=0.0001 contorno=ceron}). "
         "IC de 95\\,\\%: $P_f(1 \\pm 1{,}96\\,\\mathrm{CoV})$ (no artigo, com o CoV de 2,77\\,\\% que ele dá). $z$: a "
         "diferença contra o artigo, em desvios padrão. $\\Gamma$ det.: o determinístico (as médias de $c'$, $\\varphi'$ e "
         f"$k$); $^*$: extrapolado dos níveis {len(nv) - 3} a {len(nv) - 1} do GIM ("
         + "; ".join(num(x, 3) for x in nv[-3:]) + "). Na última linha, a referência das outras seções, com a água "
         "antiga.}\\label{tab:ceronagua}\n"
         "\\begin{tabular}{lrrrrr}\\toprule\nmétodo e amostra & $P_f$ & IC de 95\\,\\% & $z$ & $\\Gamma$ det. & "
         "média de $\\Gamma$ \\\\\\midrule\n"]
    for nome, p, c, zz, gd, gm in linhas:
        T.append(f"{nome} & {pct(p, 2)} & {num(100 * p * (1 - 1.96 * c), 2)}--{num(100 * p * (1 + 1.96 * c), 2)}\\,\\% & "
                 f"{sinal(zz, 1) if zz is not None else '--'} & {gd} & {gm} \\\\\n")
    T.append(f"\\midrule\nMEF (GIM, 3 níveis): MC direto de {num(ref['N'])}, água antiga & {pct(ref['Pf'], 2)} & "
             f"{num(100 * ref['Pf'] * (1 - 1.96 * ref['CoV']), 2)}--{num(100 * ref['Pf'] * (1 + 1.96 * ref['CoV']), 2)}"
             f"\\,\\% & {sinal(zart(ref['Pf'], ref['CoV']), 1)} & {num(gvelho, 3)} & "
             f"{num(sum(fs_velho) / len(fs_velho), 3)} \\\\\n")
    T.append("\\bottomrule\\end{tabular}\\end{table}\n")
    open(os.path.join(TAB, "ceron_agua.tex"), "w").write("".join(T))
    Nn.clear()
    cmd("ceaLsPf", pct(d["Pf"], 2))
    cmd("ceaLsLo", pct(d["Pf"] * (1 - 1.96 * d["CoV"]), 2))
    cmd("ceaLsHi", pct(d["Pf"] * (1 + 1.96 * d["CoV"]), 2))
    cmd("ceaLsCov", pct(d["CoV"], 1))
    cmd("ceaLsN", str(d["N"]))
    cmd("ceaLsAn", str(d["analises"]))
    cmd("ceaLsZ", sinal(zart(d["Pf"], d["CoV"]), 1))
    cmd("ceaForm", str(fo["avaliacoes"]))
    cmd("ceaBeta", num(fo["beta"], 3))
    cmd("ceaFormPf", pct(cauda(fo["beta"]), 1))
    cmd("ceaMcN", num(nmc))
    cmd("ceaMcPf", pct(p_r, 2))
    cmd("ceaMcPfEx", pct(p_e, 2))
    cmd("ceaMcIc", pct(1.96 * math.sqrt(p_e * (1 - p_e) / nmc), 2))
    cmd("ceaMcRef", str(len(ex)))
    cmd("ceaMcNovas", str(novas))
    # quanto a busca exata mudou: (Gamma rapido - Gamma exato) / Gamma rapido; e a P_f, em pontos
    red = [(a - b) / a for a, b in ex.values()]
    cmd("ceaExMuda", str(sum(1 for x in red if x > 0.01)))
    cmd("ceaExMax", num(100 * max(red), 1))
    cmd("ceaMcSobe", num(100 * (p_e - p_r), 2))
    cmd("ceaMcZ", sinal(zart(p_e, cv(p_e, nmc)), 1))
    cmd("ceaGdetPz", num(gpz, 3))
    cmd("ceaGniveis", "; ".join(num(x, 3) for x in nv[-3:]))
    cmd("ceaGnivelFim", num(nv[-1], 3))
    cmd("ceaGdetLa", num(gla, 3))
    cmd("ceaGdetArt", num(art["Gdet"], 3))
    cmd("ceaGmedLa", num(sum(G.values()) / nmc, 3))
    cmd("ceaGmedArt", num(art["Gmed"], 3))
    cmd("ceaPfArt", pct(art["Pf"], 2))
    cmd("ceaGvelho", num(gvelho, 3))
    cmd("ceaZvelho", sinal(zart(ref["Pf"], ref["CoV"]), 1))
    cmd("ceaLsDifPt", num(100 * abs(d["Pf"] - art["Pf"]), 1))     # em pontos percentuais
    cmd("ceaMcDifPt", num(100 * abs(p_e - art["Pf"]), 1))
    open(os.path.join(AQUI, "numeros.tex"), "a").write("".join(Nn))
    print(f"Ceron, agua do artigo: LS {d['N']} linhas {100 * d['Pf']:.3f}% (CoV {100 * d['CoV']:.2f}%), {d['analises']} "
          f"analises + FORM {fo['avaliacoes']}; MC LA {nmc}: {100 * p_r:.2f}% -> exato {100 * p_e:.2f}% ({novas} novas); "
          f"Gamma det PZ {gpz:.4f} LA {gla:.4f} artigo {art['Gdet']}")


srm_gim()
ceron_agua()
