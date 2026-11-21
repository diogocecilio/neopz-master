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
  SRM x GIM x analise limite, pareado (srm_gim e la_mef, no fim): a campanha etapa3b_gim do fermi, terminada
      (fermi_etapa3b_gim/etapa3b_gim: r0x_nys_<centro>_s21 e ls/r0x_<dir>), contra as mesmas realizacoes e linhas do SRM;
      a analise limite do fermi (fermi_cho_la/h0.025: direto, ismc, ls, form; o direto r07 ainda roda), os pilotos com
      h = 0,1 e 0,025 (Projects/CeronCinematico/cho/resultados/piloto e piloto_h0.025) e o MC simples da analise limite
      com o campo de Nystrom, sem deslocamento, rodado no PC (cho_la_pc/h0.025/direto_nys_r05_s4242.csv); o MEF com
      r = 0 e o da verificacao de 21/11 V12 (constantes MEF_R0).
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
# fermi, terminada) e pela analise limite (Ceron 2025, Projects/CeronCinematico/cho; fermi_cho_la) contra as mesmas
# realizacoes e linhas do SRM (etapa 3b e etapa3b_ls); o MC direto do SRM e as 5 sementes juntas ficam como referencia.
# O resto da analise limite (MC direto, FORM, o deslocamento FS_LA/FS_MEF e o teste da grade) esta em la_mef, no fim. E o
# Ceron com a agua do artigo (contorno=ceron agua=sat trebaix=0.0001): o line sampling do MEF e o MC de 10 000 pela
# analise limite em Python.
# Saidas: tab/srm_gim_ismc.tex, tab/srm_gim_ls.tex, tab/ceron_agua.tex, figs/srm_gim.pdf e os numeros (numeros.tex).
# =====================================================================================================================
import time  # noqa: E402

import matplotlib.patches  # noqa: E402
import pareado_gim as pg  # noqa: E402  (Projects/Cho2010: a leitura pareada dos mc_<a>_<b>.csv e dos linhas.csv)

GIMC = os.path.join(R, "fermi_etapa3b_gim", "etapa3b_gim")
LA_PILOTO = os.path.join(RAIZ, "Projects", "CeronCinematico", "cho", "resultados", "piloto_h0.025")
LA_H01 = os.path.join(RAIZ, "Projects", "CeronCinematico", "cho", "resultados", "piloto")   # o piloto antigo, grade h = 0,1
# a analise limite no fermi (jobs 14991 e 14992, baixados em 21/11/2026): direto/<r>.csv (KL analitica, semente 2010; o
# r07 ainda roda), ismc/<r>_<centro>_s21.csv (Nystrom), ls/<r>_<direcao>/ (semente 7), form/<r>_M<M>/
LA_RES = os.path.join(R, "fermi_cho_la", "h0.025")
RTEX = {"r05": "$-0{,}5$", "r07": "$-0{,}7$", "r00": "$0$"}


def le_la(arq):
    """real -> {g: Gamma, fs: FS (nan sem FS), lw: log_peso, modo} de um CSV do mc_la.py ({} se nao existe; linhas cortadas
    ou incompletas ficam de fora, como no mc_la.ler)"""
    out = {}
    if not os.path.exists(arq):
        return out
    for x in csv.DictReader(open(arq, newline="")):
        try:
            out[int(x["real"])] = {"g": float(x["Gamma"]), "fs": float(x["FS"]) if x.get("FS") else float("nan"),
                                   "lw": float(x["log_peso"] or 0.), "modo": x["modo"]}
        except (ValueError, KeyError, TypeError):
            pass
    return out


def falha_la(x):
    """falha na LA: Gamma < 1; com Gamma = inf (nenhum mecanismo admissivel) o FS decide (como o mc_la.estat)"""
    return x["g"] < 1. or (not math.isfinite(x["g"]) and math.isfinite(x["fs"]) and x["fs"] < 1.)


def razao(y, x):
    """r = media(y)/media(x) nas mesmas amostras e o seu erro padrao (estimador de razao, pareado):
    sqrt(var(y - r x)/N)/media(x)"""
    n = len(x)
    mx, my = sum(x) / n, sum(y) / n
    r = my / mx
    d = [a - r * b for a, b in zip(y, x)]
    md = sum(d) / n
    return r, math.sqrt(sum((t - md) ** 2 for t in d) / (n - 1) / n) / mx


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
    """a comparacao pareada SRM x GIM x analise limite com uma semente (ISMC s21 e line sampling s7): as duas tabelas, a
    figura e os numeros; devolve os casos (ism, lin) para o la_mef"""
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
            # o GIM com todas as suas realizacoes (a campanha foi alem das do SRM)
            pGa, cGa = est([math.exp(B[n]["lw"]) * (B[n]["f"] < 1.) for n in sorted(B)])
            ism[-1].update(nGa=len(B), pGa=pGa, cGa=cGa, soS_n=[n for n in com if A[n]["f"] < 1. <= B[n]["f"]],
                           soG_n=[n for n in com if B[n]["f"] < 1. <= A[n]["f"]])
            # a analise limite (LA) nas mesmas realizacoes: o mesmo xi e o mesmo peso (o log_peso da LA confere com o do
            # SRM); FS da LA so onde Gamma < 1,3 (todas as perto da falha)
            L = le_la(os.path.join(LA_RES, "ismc", f"{k}_mef_{d5['cen']}_s21.csv"))
            cl = [n for n in com if n in L]
            if len(cl) >= 10:
                wl = [math.exp(A[n]["lw"]) for n in cl]
                iSl = [x * (A[n]["f"] < 1.) for x, n in zip(wl, cl)]
                iL = [x * falha_la(L[n]) for x, n in zip(wl, cl)]
                pL, cL = est(iL)
                rL, eL = razao(iL, iSl)
                pLa, cLa = est([math.exp(x["lw"]) * falha_la(x) for _, x in sorted(L.items())])
                ism[-1].update(la=True, NL=len(cl), pL=pL, cL=cL, rL=rL, eL=eL, NLa=len(L), pLa=pLa, cLa=cLa,
                               soM=sum(1 for n in cl if A[n]["f"] < 1. and not falha_la(L[n])),
                               soL=sum(1 for n in cl if falha_la(L[n]) and not A[n]["f"] < 1.),
                               dlwL=max(abs(A[n]["lw"] - L[n]["lw"]) for n in cl),
                               # os pares (FS_MEF, FS_LA, Gamma_LA, w) com FS da LA, para o deslocamento
                               parLA=[(A[n]["f"], L[n]["fs"], L[n]["g"], math.exp(A[n]["lw"])) for n in cl
                                      if math.isfinite(L[n]["fs"])],
                               # (FS_MEF, Gamma_LA, w) de todas, para o ajuste de ln Gamma_LA x ln FS_MEF
                               parG=[(A[n]["f"], L[n]["g"], math.exp(A[n]["lw"])) for n in cl
                                     if math.isfinite(L[n]["g"]) and L[n]["g"] > 0.])
            else:
                ism[-1].update(la=False)
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
            # o GIM com todas as linhas (parada pelo CoV de 5 % com nmin = 100; as linhas em voo terminam e entram) e o
            # total para CoV de 5 % como nas tabelas dos metodos: FORM + 2 da inclinacao + linhas x analises por linha
            pGa, cGa = est([cauda(B[i][0]) for i in B])
            naGa = sum(B[i][1] for i in B) / len(B)
            lin[-1].update(nGa=len(B), pGa=pGa, cGa=cGa, naGa=naGa, cTodasG=[B[i][0] for i in B],
                           n5Ga=dl["form"]["n"] + 2 + naGa * len(B) * (cGa / ALVO) ** 2,
                           n5S=dl["form"]["n"] + 2 + dl["porlinha"] * dl["N"] * (dl["CoV"] / ALVO) ** 2)
            # a analise limite nas mesmas linhas (o mesmo z e a mesma direcao do MEF; ls_la.py, semente 7)
            pl = os.path.join(LA_RES, "ls", f"{k}_mef_{dl['cen']}")
            if os.path.exists(os.path.join(pl, "linhas.csv")):
                Lc = pg.le_ls(pl)
                cl = [i for i in com if i in Lc]
                yL, xS = [cauda(Lc[i][0]) for i in cl], [cauda(A[i][0]) for i in cl]
                pL, cL = est(yL)
                rL, eL = razao(yL, xS)
                mdl, sdl = mdp([Lc[i][0] - A[i][0] for i in cl])
                inc = json.load(open(os.path.join(pl, "inclinacao.json")))
                lin[-1].update(la=True, NL=len(cl), pL=pL, cL=cL, rL=rL, eL=eL, mdcL=mdl, sdcL=sdl,
                               naL=sum(Lc[i][1] for i in cl) / len(cl), incL=inc)
            else:
                lin[-1].update(la=False)
    if not ism or not lin:
        sys.exit("srm_gim: faltam os dados do GIM em " + GIMC)
    data = time.strftime("%d/%m/%Y %H:%M", time.localtime(max(datas)))
    nmin, nmax = min(x["N"] for x in ism), max(x["N"] for x in ism)

    def rot_r(L, i):
        return RTEX[L[i]["k"]] if i == 0 or L[i]["k"] != L[i - 1]["k"] else ""

    # ---- tabela: ISMC (SRM x GIM x LA; numa pagina em paisagem com a do line sampling: o .tex poe as duas no
    # mesmo landscape)
    def cel(p, c):
        return f"{pf_txt(p)} ({num(100 * c, 1)})"

    T = ["% gerado por gera_dados.py (srm_gim): nao editar a mao\n\\begin{table}[H]\\centering\\footnotesize"
         "\\setlength{\\tabcolsep}{3pt}\n"
         f"\\caption{{A ISMC com uma semente (a 21), pareada: o SRM, o GIM e a análise limite (LA) nas mesmas realizações "
         f"(Cho, base de Nyström de 2000 modos; dados de {data}). $N$: todas as realizações do SRM na semente 21 ({nmin} a "
         f"{nmax}), que têm todas o par no GIM e na LA. $P_f$ em \\%: a média de $I w$ (entre parênteses, o CoV, em \\%). "
         "Dif.\\ GIM: $(P_\\text{GIM} - P_\\text{SRM})/P_\\text{SRM}$ $\\pm$ o erro pareado, "
         "$\\mathrm{dp}[w\\,(I_\\text{GIM} - I_\\text{SRM})]/(\\sqrt N\\,P_\\text{SRM})$. Discord.: os pares em que só o "
         "SRM falha / só o outro falha. LA/SRM: a razão das $P_f$ $\\pm$ o seu erro pareado (estimador de razão). "
         "``GIM, todas'': a ISMC do GIM com todas as suas realizações (a campanha foi além das do SRM). Referências (SRM): "
         "as 5 sementes juntas e o MC direto.}\\label{tab:sgismc}\n"
         "\\begin{tabular}{llrrrrcrrcrrr}\\toprule\n"
         "$r$ & centro & $N$ & \\makecell{SRM s21\\\\$P_f$ (CoV)} & \\makecell{GIM s21\\\\$P_f$ (CoV)} & "
         "\\makecell{dif.\\ GIM\\\\(\\%)} & \\makecell{discord.\\\\SRM/GIM} & \\makecell{LA s21\\\\$P_f$ (CoV)} & "
         "LA/SRM & \\makecell{discord.\\\\SRM/LA} & \\makecell{GIM, todas\\\\$N$: $P_f$ (CoV)} & "
         "\\makecell{SRM\\\\5 sementes} & \\makecell{MC direto\\\\(SRM)} \\\\\\midrule\n"]
    for i, x in enumerate(ism):
        if i and x["k"] != ism[i - 1]["k"]:
            T.append("\\midrule\n")
        la_ = (f"{cel(x['pL'], x['cL'])} & {num(x['rL'], 3)} $\\pm$ {num(x['eL'], 3)} & {x['soM']}/{x['soL']}"
               if x["la"] else "-- & -- & --")
        T.append(f"{rot_r(ism, i)} & {curto(x['nome'])} & {x['N']} & {cel(x['pS'], x['cS'])} & {cel(x['pG'], x['cG'])} & "
                 f"{sinal(100 * x['dif'])} $\\pm$ {num(100 * x['edif'], 1)} & {len(x['soS'])}/{len(x['soG'])} & {la_} & "
                 f"{x['nGa']}: {cel(x['pGa'], x['cGa'])} & {pf_txt(x['p5'])} & {pf_txt(x['ref']['Pf'])} \\\\\n")
    T.append("\\bottomrule\\end{tabular}\\end{table}\n")
    open(os.path.join(TAB, "srm_gim_ismc.tex"), "w").write("".join(T))
    # ---- tabela: line sampling
    T = ["% gerado por gera_dados.py (srm_gim): nao editar a mao\n\\begin{table}[H]\\centering\\footnotesize"
         "\\setlength{\\tabcolsep}{2pt}\n"
         f"\\caption{{O line sampling com uma semente (a 7), pareado: as mesmas {min(x['N'] for x in lin)} linhas (o mesmo "
         "$\\bm z$ de cada uma) e as mesmas direções, com o SRM, o GIM e a análise limite (LA). $P_f$ em \\% (CoV em \\%). "
         "Dif.\\ GIM e LA/SRM: como na Tabela~\\ref{tab:sgismc}, com as $\\Phi(-c_\\ell)$ de cada linha. "
         "$\\Delta c$: onde a linha entra na falha, menos o $c_\\text{SRM}$ (média e dp das linhas). ``GIM, todas as "
         "linhas'': o line sampling do GIM até o fim (parada no CoV de 5\\,\\% com pelo menos 100 linhas; as linhas "
         "em voo terminam e entram) e o total de análises para CoV de 5\\,\\% como na Tabela~\\ref{tab:r05} (FORM, "
         "inclinação e linhas; o do SRM, com as 100 linhas, entre parênteses).}\\label{tab:sgls}\n"
         "\\begin{tabular}{llrrrrrrrrrr}\\toprule\n"
         "$r$ & direção & linhas & \\makecell{SRM\\\\$P_f$ (CoV)} & \\makecell{GIM\\\\$P_f$ (CoV)} & "
         "\\makecell{dif.\\ GIM\\\\(\\%)} & \\makecell{$\\Delta c$ GIM\\\\média (dp)} & \\makecell{LA\\\\$P_f$ (CoV)} & "
         "LA/SRM & \\makecell{$\\Delta c$ LA\\\\média (dp)} & \\makecell{análises/linha\\\\SRM/GIM/LA} & "
         "\\makecell{GIM, todas as linhas: linhas,\\\\$P_f$ (CoV), p/ CoV 5\\,\\% (SRM)} \\\\\\midrule\n"]
    for i, x in enumerate(lin):
        if i and x["k"] != lin[i - 1]["k"]:
            T.append("\\midrule\n")
        la_ = (f"{cel(x['pL'], x['cL'])} & {num(x['rL'], 3)} $\\pm$ {num(x['eL'], 3)} & {sinal(x['mdcL'], 3)} "
               f"({num(x['sdcL'], 3)})" if x["la"] else "-- & -- & --")
        naL = f"/{num(x['naL'], 2)}" if x["la"] else ""
        T.append(f"{rot_r(lin, i)} & {curto(x['nome'])} & {x['N']} & {cel(x['pS'], x['cS'])} & {cel(x['pG'], x['cG'])} & "
                 f"{sinal(100 * x['dif'])} $\\pm$ {num(100 * x['edif'], 1)} & {sinal(x['mdc'], 3)} ({num(x['sdc'], 3)}) & "
                 f"{la_} & {num(x['naS'], 2)}/{num(x['naG'], 2)}{naL} & {x['nGa']}, {cel(x['pGa'], x['cGa'])}, "
                 f"{num(x['n5Ga'])} ({num(x['n5S'])}) \\\\\n")
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
    # o GIM com todas as suas realizacoes (ISMC) e todas as linhas (LS): a campanha foi alem das do SRM
    cmd("sgNRaMin", str(min(x["N"] for x in ism if x["k"] == "r05")))     # as realizacoes do SRM s21 no r = -0,5
    cmd("sgNRaMax", str(max(x["N"] for x in ism if x["k"] == "r05")))
    cmd("sgGaNmin", str(min(x["nGa"] for x in ism)))
    cmd("sgGaNmax", str(max(x["nGa"] for x in ism)))
    for k, sig in (("r05", "Ra"), ("r07", "Rb")):
        sub = [x for x in ism if x["k"] == k]
        cmd(f"sgGa{sig}Min", pf_txt(min(x["pGa"] for x in sub)))
        cmd(f"sgGa{sig}Max", pf_txt(max(x["pGa"] for x in sub)))
        cmd(f"sgGaCov{sig}Min", num(100 * min(x["cGa"] for x in sub), 1))
        cmd(f"sgGaCov{sig}Max", num(100 * max(x["cGa"] for x in sub), 1))
        # o z do GIM sozinho contra o MC direto do SRM
        cmd(f"sgGaZ{sig}Max", num(max(abs(z(x["pGa"], x["cGa"], P[k]["ref"])) for x in sub), 1))
        for x in [y for y in lin if y["k"] == k]:
            D = {"M20": "Mv", "gim2000": "Gm"}[x["cen"]]
            cmd(f"sglGaPf{sig}{D}", pf_txt(x["pGa"]))
            cmd(f"sglGaN{sig}{D}", str(x["nGa"]))
            cmd(f"sglGaCov{sig}{D}", num(100 * x["cGa"], 1))
            cmd(f"sglGaNc{sig}{D}", num(x["n5Ga"]))
            cmd(f"sglSNc{sig}{D}", num(x["n5S"]))
    cmd("sgGaZMax", num(max(abs(z(x["pGa"], x["cGa"], P[x["k"]]["ref"])) for x in ism), 1))

    # o custo para CoV de 5 % com todas as linhas do GIM contra o estimado com as 100 do SRM: a razao, as analises por
    # linha (GIM com todas / SRM), o CoV de uma linha (GIM nas 100 em comum -> todas) e o erro relativo de CoV^2 que a
    # curtose kappa das proprias 100 linhas do SRM preve, sqrt((kappa - 1)/N); e a maior curtose com todas as linhas do GIM
    def curtose(v):
        n = len(v)
        m = sum(v) / n
        m2 = sum((a - m) ** 2 for a in v) / n
        return sum((a - m) ** 4 for a in v) / n / m2 ** 2

    cmd("sglGaCustoMin", num(100 * min(x["n5Ga"] / x["n5S"] - 1. for x in lin), 0))
    cmd("sglGaCustoMax", num(100 * max(x["n5Ga"] / x["n5S"] - 1. for x in lin), 0))
    cmd("sglNaRazMin", num(100 * min(x["naGa"] / x["naS"] - 1. for x in lin), 0))
    cmd("sglNaRazMax", num(100 * max(x["naGa"] / x["naS"] - 1. for x in lin), 0))
    # o CoV^2 de uma linha: o do GIM com todas as linhas contra o do SRM com as 100 (a parte do custo que nao e analises
    # por linha)
    c2 = [(x["cGa"] ** 2 * x["nGa"]) / (x["cS"] ** 2 * x["N"]) - 1. for x in lin]
    cmd("sglCovDoisMin", num(100 * min(c2), 0))
    cmd("sglCovDoisMax", num(100 * max(c2), 0))
    c2g = [v for v, x in zip(c2, lin) if x["cen"] == "gim2000"]          # so a direcao modal 2000 (a das tabelas)
    cmd("sglCovDoisGmMin", num(100 * min(c2g), 0))
    cmd("sglCovDoisGmMax", num(100 * max(c2g), 0))
    ek = [math.sqrt((curtose([cauda(c) for c in x["c"]]) - 1.) / len(x["c"])) for x in lin]
    cmd("sglCurtErrMin", num(100 * min(ek), 0))
    cmd("sglCurtErrMax", num(100 * max(ek), 0))
    for x in lin:
        if x["k"] == "r07":
            D = {"M20": "Mv", "gim2000": "Gm"}[x["cen"]]
            cmd(f"sglClRb{D}Cem", num(100 * x["cG"] * math.sqrt(x["N"]), 0))
            cmd(f"sglClRb{D}Todas", num(100 * x["cGa"] * math.sqrt(x["nGa"]), 0))
            cmd(f"sglCurtRb{D}Todas", num(curtose([cauda(B) for B in x["cTodasG"]]), 0))
    # os discordantes de todos os casos juntos: o sinal (so o SRM falha x so o GIM falha), em realizacoes distintas (a
    # mesma realizacao n de um r pode discordar em mais de um centro), e a chance de um sinal tao forte com 50 % para
    # cada lado (cauda binomial)
    cmd("sgDiscFSmed", teto(sum(abs(f - 1.) for f, g in dd) / len(dd), 4) if dd else "--")
    ds = {(x["k"], n) for x in ism for n in x["soS_n"]}
    dg = {(x["k"], n) for x in ism for n in x["soG_n"]}
    na_, nb_ = len(ds), len(dg)
    cmd("sgDistSrm", str(na_))
    cmd("sgDistGim", str(nb_))
    cauda_b = sum(math.comb(na_ + nb_, j) for j in range(max(na_, nb_), na_ + nb_ + 1)) / 2 ** (na_ + nb_) if na_ + nb_ else 1.
    cmd("sgDistP", f"{num(1e4 * cauda_b, 0)} em 10\\,000" if cauda_b < 0.01 else f"{num(100 * cauda_b, 0)}\\,\\%")
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
    return ism, lin


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


# =====================================================================================================================
# A analise limite (LA: Ceron 2025, Projects/CeronCinematico/cho, grade h = 0,025 m) contra o MEF (SRM, tri2 + 4 niveis):
# todos os metodos (MC direto com a KL analitica; ISMC s21, LS s7 e FORM com a base de Nystrom), o deslocamento
# FS_LA/FS_MEF - 1 nos pares (por faixa de FS_MEF, nos dois campos), a linha do ponto de projeto, o teste da grade
# (h = 0,1 x 0,025 nas mesmas realizacoes, pelos pilotos), a P_f da LA prevista pelo MEF deslocado e os modos de ruptura.
# Saidas: tab/la_mef.tex, tab/la_teste.tex, tab/la_modos.tex, figs/la_mef.pdf e os numeros lm* (numeros.tex).
# =====================================================================================================================
FS_DET = {"LA": 1.2032, "MEF": 1.2044}      # campo medio: FS da LA (forma fechada; testes.py det) e do SRM (tri2 + 4)
CERON_R0 = {"Pf": 0.0550, "CoV": 0.019, "N": 50000, "modos": (80.80, 19.20, 0.00)}   # Ceron (2025), Sec. 5.3.2 e Fig. 20
CHO_R0 = 0.0637                             # o que o Ceron (2025) da como a P_f do Cho com r = 0
# os modos do MEF pelas bandas de cisalhamento (verificacoes de 21/11, V07): (todas, falhas), cada um (pe, face, base) %
MEF_MODOS = {"r05": ((97.99, 0.06, 1.95), (99.51, 0.00, 0.49)), "r07": ((99.26, 0.01, 0.73), (99.75, 0.00, 0.25))}
FAIXAS = [(0., 0.95), (0.95, 1.0), (1.0, 1.05), (1.05, 1.10), (1.10, 1.20), (1.20, 9.)]
# a LA com o campo de Nystrom SEM deslocamento: MC simples, so Gamma, r = -0,5, semente 4242 (que nenhuma outra rodada
# usa), 1500 realizacoes, no PC em 21/11/2026 (revisao do relatorio; mc_la.py campo=nystrom sem centro=, 2 processos,
# ~38 min; o comando esta no LEIAME.txt da pasta). A distribuicao de Gamma sem FORM nem ISMC, contra o direto analitico.
LA_NYS_PC = os.path.join(R, "cho_la_pc", "h0.025", "direto_nys_r05_s4242.csv")
NIV_G = (1.0, 1.05, 1.1, 1.2, 1.3, 1.4, 1.5)   # niveis de Gamma da LA
NIV_FS = (1.0, 1.05, 1.1)                  # niveis de FS do MEF
QS = (0.01, 0.05, 0.10, 0.25, 0.50)        # quantis de Gamma
# o MEF com r = 0 (verificacoes de 21/11, V12/resumo_v12.txt): line sampling do SRM, base de Nystrom, semente 7, 84
# linhas (P_f, IC 95 %, CoV); a P_f da LA prevista pela V12 com o deslocamento do piloto analitico (regressao dos 400
# pares, ln FS_LA - ln FS_MEF = 0,0243 em FS_MEF = 1) aplicado linha a linha (e o IC das linhas); o Cho, Fig. 9
# digitalizada na V12
MEF_R0 = {"Pf": 0.074045, "lo": 0.06807, "hi": 0.08002, "CoV": 0.041159, "prev": 0.04783, "prev_lo": 0.04371,
          "prev_hi": 0.05195, "cho_fig9": 0.0662}


def modos(xs):
    """(pe, face, base, nenhum) em % das realizacoes dadas (sem peso)"""
    n = len(xs)
    return tuple(100. * sum(1 for x in xs if x["modo"] == m) / n for m in ("pe", "face", "base", "nenhum")) if n else None


def faixa_rel(pares, lo, hi):
    """FS_LA/FS_MEF - 1 nos pares (FS_MEF, FS_LA, Gamma_LA, w) com lo <= FS_MEF < hi: n, media com peso, dp"""
    s = [(fl / fm - 1., w) for fm, fl, g, w in pares if lo <= fm < hi and math.isfinite(fl)]
    if len(s) < 2:
        return None
    sw = sum(w for _, w in s)
    return {"n": len(s), "m": sum(v * w for v, w in s) / sw, "mu": sum(v for v, _ in s) / len(s),
            "dp": mdp([v for v, _ in s])[1], "ess": sw * sw / sum(w * w for _, w in s)}


def gamma_em_1(pares, lo=0.9, hi=1.1):
    """o ajuste (minimos quadrados com o peso w) de ln Gamma_LA = a + k ln FS_MEF com lo < FS_MEF < hi: Gamma_LA quando
    FS_MEF = 1 (exp a), k e os pares. Usa todas as realizacoes (o FS da LA so existe onde Gamma < 1,3)."""
    s = [(math.log(f), math.log(g), w) for f, g, w in pares if lo < f < hi]
    if len(s) < 10:
        return None
    sw = sum(w for _, _, w in s)
    mx, my = sum(x * w for x, _, w in s) / sw, sum(y * w for _, y, w in s) / sw
    k = sum(w * (x - mx) * (y - my) for x, y, w in s) / sum(w * (x - mx) ** 2 for x, _, w in s)
    return {"g1": math.exp(my - k * mx), "k": k, "n": len(s)}


def le_pares_la(arq):
    """real -> (Gamma_LA, FS_LA, FS_MEF) de um *_pares.csv do pareado.py juntar"""
    return {int(x["real"]): (float(x["Gamma_LA"]), float(x["FS_LA"]) if x["FS_LA"] else float("nan"),
                             float(x["FS_MEF"])) for x in csv.DictReader(open(arq, newline=""))}


def quantil(v, p):
    """o quantil p da lista ordenada v (o valor de ordem ceil(p N))"""
    return v[min(len(v) - 1, max(0, int(math.ceil(p * len(v))) - 1))]


def la_mef(ism, lin):
    """a analise limite contra o MEF: tabelas, figura e numeros (o que faltar fica de fora)"""
    if not os.path.isdir(LA_RES):
        print("la_mef: sem", LA_RES)
        return
    # ---- (1) MC direto da LA (KL analitica, semente 2010: as mesmas realizacoes do direto do MEF) e os pares
    DIR = {}
    for k in ("r05", "r07", "r00"):
        L = le_la(os.path.join(LA_RES, "direto", f"{k}.csv"))
        if not L:
            continue
        xs = list(L.values())
        I = [1. * falha_la(x) for x in xs]
        p, cv = est(I)
        G = sorted(x["g"] for x in xs)                       # com os inf (nenhum mecanismo) no fim
        Gf = [g for g in G if math.isfinite(g)]
        d = {"N": len(xs), "fal": int(sum(I)), "Pf": p, "CoV": cv, "Gmed": sum(Gf) / len(Gf), "G": G,
             "mod": modos(xs), "modf": modos([x for x in xs if falha_la(x)]), "ninf": len(G) - len(Gf)}
        if k in P:
            Mf, _ = pg.le_mc(P[k]["refpasta"])
            com = sorted(set(L) & set(Mf))
            iL = [1. * falha_la(L[n]) for n in com]
            iM = [1. * (Mf[n]["f"] < 1.) for n in com]
            r_, e_ = razao(iL, iM)
            d.update(Npar=len(com), fLpar=int(sum(iL)), fMpar=int(sum(iM)), rpar=r_, epar=e_,
                     soL=sum(1 for a, b in zip(iL, iM) if a and not b), pMpar=sum(iM) / len(com),
                     pares=[(Mf[n]["f"], L[n]["fs"], L[n]["g"], 1.) for n in com if math.isfinite(L[n]["fs"])],
                     parG=[(Mf[n]["f"], L[n]["g"], 1.) for n in com if math.isfinite(L[n]["g"]) and L[n]["g"] > 0.],
                     fsM=sorted(Mf[n]["f"] for n in Mf), acima=all(L[n]["fs"] > Mf[n]["f"] for n in com
                                                                  if math.isfinite(L[n]["fs"])),
                     # (FS_MEF, real, tem FS da LA?, FS_LA) de todas as em comum: o FS da LA so existe nas 2400
                     # primeiras e onde Gamma < 1,3 (r05) ou 1,1 (r07)
                     parN=[(Mf[n]["f"], n, math.isfinite(L[n]["fs"]), L[n]["fs"]) for n in com])
        DIR[k] = d
    # ---- (2) ISMC da LA (semente 21, 2500 realizacoes, campo de Nystrom) nos 4 centros
    ISM = {}
    for k in ("r05", "r07"):
        for cen in ("mef_M20", "mef_M100", "mef_gim20", "la_M20"):
            L = le_la(os.path.join(LA_RES, "ismc", f"{k}_{cen}_s21.csv"))
            if not L:
                continue
            xs = [L[n] for n in sorted(L)]
            p, cv = est([math.exp(x["lw"]) * falha_la(x) for x in xs])
            ISM[(k, cen)] = {"N": len(xs), "Pf": p, "CoV": cv, "modf": modos([x for x in xs if falha_la(x)]),
                             "ref": next((d for d in P[k]["ismc"] if "mef_" + d["cen"] == cen), None), "xs": xs}
    # ---- (2b) a distribuicao de Gamma da LA nas duas bases, r = -0,5: o direto analitico (fermi), o MC simples de
    # Nystrom (PC, sem FORM nem ISMC) e, na cauda, as 4 ISMC de Nystrom com o peso (a media das 4, que tem a mesma
    # semente: o CoV medio, sem dividir por raiz de 4); e a do FS do MEF nas duas bases (o direto analitico x as 5
    # sementes da ISMC de Nystrom, por centro, e a media dos 3 centros)
    DG = {}
    Lny = le_la(LA_NYS_PC)
    if Lny and "r05" in DIR:
        def qic(G, p):
            """o quantil p da lista ordenada G e o seu IC de 95 % pelas estatisticas de ordem"""
            n = len(G)
            h = 1.96 * math.sqrt(n * p * (1. - p))
            return (quantil(G, p), G[max(0, int(math.floor(n * p - h)) - 1)], G[min(n - 1, int(math.ceil(n * p + h)) - 1)])

        for f, G in (("an", DIR["r05"]["G"]), ("ny", sorted(x["g"] for x in Lny.values()))):
            n = len(G)
            lg = [math.log(g) for g in G if math.isfinite(g)]
            m_, s_ = mdp(lg)
            DG[f] = {"N": n, "F": {g: sum(1 for x in G if x < g) / n for g in NIV_G}, "q": {p: qic(G, p) for p in QS},
                     "mlg": m_, "elg": s_ / math.sqrt(len(lg))}
        I = [1. * falha_la(x) for x in Lny.values()]
        p_, c_ = est(I)
        DG["ny"].update(fal=int(sum(I)), Pf=p_, CoV=c_)
        xs4 = [ISM[("r05", c)]["xs"] for c in ("mef_M20", "mef_M100", "mef_gim20", "la_M20") if ("r05", c) in ISM]
        if xs4:
            DG["is"] = {}
            for g0 in NIV_G:
                v = [est([math.exp(x["lw"]) * (x["g"] < g0) for x in xs]) for xs in xs4]
                DG["is"][g0] = (sum(p for p, _ in v) / len(v), sum(c for _, c in v) / len(v))
        fsM = DIR["r05"].get("fsM")
        if fsM:
            vals = {s: [] for s in NIV_FS}
            for cen in ("M20", "M100", "gim20"):
                pool = []
                for sm in SEM:
                    pa = os.path.join(E3B, f"r05_nys_{cen}_s{sm}")
                    if os.path.isdir(pa):
                        pool += list(pg.le_mc(pa)[0].values())
                for s in NIV_FS:
                    vals[s].append(est([math.exp(x["lw"]) * (x["f"] < s) for x in pool]))
            DG["mef"] = {"an": {s: sum(1 for v in fsM if v < s) / len(fsM) for s in NIV_FS}, "Nan": len(fsM),
                         "ny": {s: (sum(p for p, _ in v) / len(v), sum(c for _, c in v) / len(v)) for s, v in vals.items()}}
    # ---- (3) line sampling da LA (100 linhas, semente 7) nas direcoes do MEF e da LA
    LSD = {}
    for k in ("r05", "r07"):
        for dr in ("mef_M20", "mef_gim2000", "la_M20", "la_M100"):
            f = os.path.join(LA_RES, "ls", f"{k}_{dr}", "ls.json")
            if not os.path.exists(f):
                continue
            j = json.load(open(f))
            LSD[(k, dr)] = {"N": j["linhas"], "Pf": j["P_f"], "CoV": j["CoV"], "na": j["analises_linha"],
                            "beta": j["beta"], "perd": j.get("perdidas", 0),
                            "ref": next((d for d in P[k]["ls"] if "mef_" + d["cen"] == dr), None),
                            "par": next((x for x in lin if x["k"] == k and "mef_" + x["cen"] == dr and x.get("la")), None)}
    # ---- (4) FORM da LA (diferencas centrais, M termos por campo) contra o FORM numerico do MEF com o mesmo M
    FO = {}
    for k in ("r05", "r07"):
        for M in (20, 100):
            fl, fm = os.path.join(LA_RES, "form", f"{k}_M{M}", "form.json"), os.path.join(E3B, f"form_{k}_M{M}", "form.json")
            if not (os.path.exists(fl) and os.path.exists(fm)):
                continue
            a, b = json.load(open(fl)), json.load(open(fm))
            cs = sum(x * y for x, y in zip(a["u"], b["u"])) / math.sqrt(sum(x * x for x in a["u"]) * sum(y * y for y in b["u"]))
            FO[(k, M)] = {"bL": a["beta"], "nL": a["avaliacoes"], "bM": b["beta"], "nM": b["avaliacoes"], "cos": cs,
                          "itL": a["iteracoes"]}
    # ---- (5) o deslocamento FS_LA/FS_MEF - 1 por faixa de FS_MEF: KL analitica (direto) x Nystrom (ISMC s21, os 3
    # centros juntos, com o peso w: estima a media sob a densidade original, como no direto)
    PP = {("an", k): DIR[k]["pares"] for k in ("r05", "r07") if k in DIR and "pares" in DIR[k]}
    PP.update({("ny", k): [p for x in ism if x["k"] == k and x.get("la") for p in x["parLA"]] for k in ("r05", "r07")})
    BAND = {(f, k, i): faixa_rel(v, lo, hi) for (f, k), v in PP.items() for i, (lo, hi) in enumerate(FAIXAS)}
    PERTO = {(f, k): faixa_rel(v, 0.95, 1.05) for (f, k), v in PP.items()}
    PG = {("an", k): DIR[k]["parG"] for k in ("r05", "r07") if k in DIR and "parG" in DIR[k]}
    PG.update({("ny", k): [p for x in ism if x["k"] == k and x.get("la") for p in x["parG"]] for k in ("r05", "r07")})
    G1 = {fk: gamma_em_1(v) for fk, v in PG.items()}
    # ---- (6) a linha do ponto de projeto do MEF (z_perp = 0: o campo so com os modos do FORM): onde a LA entra na falha,
    # pela inclinacao medida em c = beta +- 0,1 (ln Gamma linear), contra as 100 linhas sorteadas; em FS, |s_SRM| Delta c
    PTO = {}
    for x in lin:
        if x.get("la"):
            gm, gp = math.log(x["incL"]["Gamma_m"]), math.log(x["incL"]["Gamma_p"])
            c0 = -0.1 + 0.2 * gm / (gm - gp)
            PTO[(x["k"], x["cen"])] = {"dc": c0, "dfs": abs(x["sS"]) * c0, "dcl": x["mdcL"], "dfsl": abs(x["sS"]) * x["mdcL"]}
    # ---- (7) a grade: FS_LA com h = 0,1 (piloto antigo) e h = 0,025 nas mesmas realizacoes
    GR = {"an": [], "ny": [], "ny01": [], "ny025": [], "an_fm": [], "ny_fm": []}   # *_fm: o FS do MEF das mesmas
    GRr = {}                                  # por r (KL analitica): a razao das grades e o deslocamento em cada grade
    for k in ("r05", "r07"):
        a, b = os.path.join(LA_H01, f"direto_{k}_pares.csv"), os.path.join(LA_PILOTO, f"direto_{k}_pares.csv")
        if os.path.exists(a) and os.path.exists(b):
            A1, A2 = le_pares_la(a), le_pares_la(b)
            g = {"an": [], "an01": [], "an025": []}
            for n in sorted(set(A1) & set(A2)):
                f1, f2, fm = A1[n][1], A2[n][1], A2[n][2]
                if math.isfinite(f1) and math.isfinite(f2):
                    g["an"].append(f1 / f2 - 1.)
                    g["an01"].append(f1 / fm - 1.)
                    g["an025"].append(f2 / fm - 1.)
                    GR["an_fm"].append(fm)
            GRr[k] = g
            GR["an"] += g["an"]
    a = os.path.join(LA_H01, "ismc_r05_mef_M20_s21_pares.csv")
    L2 = le_la(os.path.join(LA_RES, "ismc", "r05_mef_M20_s21.csv"))
    if os.path.exists(a) and L2:
        A1 = le_pares_la(a)
        for n in sorted(A1):
            if n in L2 and math.isfinite(A1[n][1]) and math.isfinite(L2[n]["fs"]):
                GR["ny"].append(A1[n][1] / L2[n]["fs"] - 1.)
                GR["ny01"].append(A1[n][1] / A1[n][2] - 1.)
                GR["ny025"].append(L2[n]["fs"] / A1[n][2] - 1.)
                GR["ny_fm"].append(A1[n][2])
    # ---- (8) a P_f da LA prevista pelo MEF: P(FS_MEF (1 + delta) < 1) no MC direto do MEF, com o delta perto da falha
    # de cada campo, contra a P_f medida (direto analitico; a media das 3 ISMC da LA nos centros do MEF)
    PREV = {}
    for k in ("r05", "r07"):
        if k not in DIR or "fsM" not in DIR[k]:
            continue
        fsM = DIR[k]["fsM"]
        ny = [ISM[(k, c)] for c in ("mef_M20", "mef_M100", "mef_gim20") if (k, c) in ISM]
        for f in ("an", "ny"):
            if PERTO.get((f, k)):
                dl = PERTO[(f, k)]["m"]
                PREV[(k, f)] = {"d": dl, "p": sum(1 for v in fsM if v * (1. + dl) < 1.) / len(fsM)}
        PREV[(k, "obs_an")] = (DIR[k]["Pf"], DIR[k]["CoV"])
        if ny:
            # as 3 ISMC usam a mesma semente (erros correlacionados): o CoV medio, sem dividir por raiz de 3
            PREV[(k, "obs_ny")] = (sum(x["Pf"] for x in ny) / len(ny), sum(x["CoV"] for x in ny) / len(ny))
    # ---- (9) r = 0 (sem MEF): o mesmo fator em Gamma que leva o direto analitico do r = -0,5 a ISMC de Nystrom
    R0 = None
    if "r00" in DIR and "r05" in DIR and (("r05", "obs_ny") in PREV):
        fny = quantil(DIR["r05"]["G"], PREV[("r05", "obs_ny")][0])
        G0 = DIR["r00"]["G"]
        R0 = {"f": fny, "p": sum(1 for g in G0 if g < fny) / len(G0), "fc": quantil(G0, CERON_R0["Pf"])}

    # ---------------- tabela: o resumo LA x MEF
    def pc(p, c):
        return f"{pf_txt(p)} ({num(100 * c, 1)})"

    NOMEC = {"M20": "$M = 20$", "M100": "$M = 100$", "gim20": "modal 20", "gim2000": "modal 2000"}
    T = ["% gerado por gera_dados.py (la_mef): nao editar a mao\n\\begin{table}[H]\\centering\\footnotesize"
         "\\setlength{\\tabcolsep}{3pt}\n"
         "\\caption{A análise limite (LA, Ceron 2025; grade de $h = 0{,}025$\\,m) contra o MEF (SRM, \\texttt{tri2} e 4 "
         "níveis), método a método. $P_f$ em \\% (CoV em \\%). O MC direto usa a KL analítica de 2000 termos (as "
         "realizações do direto do MEF, semente 2010); o resto, a base de Nyström de 2000 modos (a ISMC da semente 21 e o "
         "line sampling da semente 7 do MEF, as mesmas realizações e linhas). MEF: no direto, o MC direto do SRM; na "
         "ISMC, as 5 sementes do SRM juntas; no line sampling, o do SRM nas mesmas linhas. LA/MEF: a razão das $P_f$; "
         "entre parênteses, a razão pareada nas mesmas amostras (Tabelas~\\ref{tab:sgismc} e~\\ref{tab:sgls}). FORM: "
         "$\\beta$ e, entre parênteses, $\\Phi(-\\beta)$ em \\%; na última coluna, a razão dos $\\Phi(-\\beta)$ e o "
         "cosseno entre os pontos de projeto da LA e do MEF. Com $r = 0$ não há MEF pareado: o MEF é o line sampling do "
         "SRM das verificações de 21/11 (V12: base de Nyström, semente 7, 84 linhas), e o Ceron e o Cho são os "
         "publicados. O direto da LA com $r = -0{,}7$ ainda roda: é parcial.}\\label{tab:lamef}\n"
         "\\begin{tabular}{llrrrl}\\toprule\n"
         "$r$ & método & amostra & \\makecell{LA\\\\$P_f$ (CoV)} & \\makecell{MEF\\\\$P_f$ (CoV)} & LA/MEF \\\\\\midrule\n"]
    for k in ("r05", "r07"):
        if k not in DIR:
            continue
        d = DIR[k]
        ref = P[k]["ref"]
        par = f" ({num(d['rpar'], 2)})" if "rpar" in d else ""
        T.append(f"{RTEX[k]} & MC direto (KL analítica){', parcial' if k == 'r07' else ''} & {num(d['N'])} & "
                 f"{pc(d['Pf'], d['CoV'])} & {pc(ref['Pf'], ref['CoV'])} & {num(d['Pf'] / ref['Pf'], 2)}{par} \\\\\n")
        for cen in ("mef_M20", "mef_M100", "mef_gim20", "la_M20"):
            x = ISM.get((k, cen))
            if not x:
                continue
            pi = next((y for y in ism if y["k"] == k and "mef_" + y["cen"] == cen and y.get("la")), None)
            nome = ("centro MEF " if cen.startswith("mef") else "centro LA ") + NOMEC[cen.split("_")[1]]
            mef = pc(x["ref"]["Pf"], x["ref"]["CoV"]) if x["ref"] else "--"
            rz = (num(x["Pf"] / x["ref"]["Pf"], 2) + (f" ({num(pi['rL'], 2)})" if pi else "")) if x["ref"] else "--"
            T.append(f" & ISMC s21, {nome} & {num(x['N'])} & {pc(x['Pf'], x['CoV'])} & {mef} & {rz} \\\\\n")
        for dr in ("mef_M20", "mef_gim2000", "la_M20", "la_M100"):
            x = LSD.get((k, dr))
            if not x:
                continue
            nome = ("direção MEF " if dr.startswith("mef") else "direção LA ") + NOMEC[dr.split("_")[1]]
            mef = pc(x["ref"]["Pf"], x["ref"]["CoV"]) if x["ref"] else "--"
            rz = (num(x["Pf"] / x["ref"]["Pf"], 2) + (f" ({num(x['par']['rL'], 2)})" if x["par"] else "")) if x["ref"] else "--"
            T.append(f" & LS s7, {nome} & {x['N']} linhas & {pc(x['Pf'], x['CoV'])} & {mef} & {rz} \\\\\n")
        for M in (20, 100):
            x = FO.get((k, M))
            if x:
                T.append(f" & FORM, $M = {M}$: $\\beta$ ($\\Phi(-\\beta)$) & {x['nL']}/{x['nM']} análises & "
                         f"{num(x['bL'], 3)} ({num(100 * cauda(x['bL']), 2)}) & {num(x['bM'], 3)} "
                         f"({num(100 * cauda(x['bM']), 2)}) & {num(cauda(x['bL']) / cauda(x['bM']), 2)}; cos "
                         f"{num(x['cos'], 4)} \\\\\n")
        T.append("\\midrule\n")
    if "r00" in DIR:
        d = DIR["r00"]
        T.append(f"$0$ & MC direto (KL analítica) & {num(d['N'])} & {pc(d['Pf'], d['CoV'])} & "
                 f"\\makecell[r]{{MEF (V12): {num(100 * MEF_R0['Pf'], 2)}\\\\"
                 f"Ceron (LA): {num(100 * CERON_R0['Pf'], 2)}\\\\Cho: {num(100 * CHO_R0, 2)}}} & "
                 f"\\makecell[l]{{{num(d['Pf'] / MEF_R0['Pf'], 2)} (MEF)\\\\{num(d['Pf'] / CERON_R0['Pf'], 2)} (Ceron)\\\\"
                 f"{num(d['Pf'] / CHO_R0, 2)} (Cho)}} \\\\\n")
    T.append("\\bottomrule\\end{tabular}\\end{table}\n")
    open(os.path.join(TAB, "la_mef.tex"), "w").write("".join(T))

    # ---------------- tabela: o teste (o deslocamento por faixa, a linha do ponto de projeto, a grade) e a previsao
    def cel_b(b):
        return f"{sinal(100 * b['m'], 2)} ({b['n']})" if b else "--"

    def cobertura(k, lo, hi):
        """as realizacoes do direto com lo <= FS_MEF < hi: quantas, quantas sem o FS da LA, e o deslocamento medio so
        nas 2400 primeiras (que tem FS em todas) e quantas sao"""
        v = [x for x in DIR[k].get("parN", []) if lo <= x[0] < hi]
        s_ = [fl / fm - 1. for fm, n, ok, fl in v if ok and n < 2400]
        return len(v), sum(1 for x in v if not x[2]), (sum(s_) / len(s_) if s_ else float("nan")), len(s_)

    c05, c07 = cobertura("r05", 0.95, 1.05), cobertura("r07", 1.0, 1.05)
    c05b = cobertura("r05", 1.05, 1.10)
    b05b = BAND.get(("an", "r05", 3))
    b07 = BAND.get(("an", "r07", 2))
    rot_f = ["$< 0{,}95$", "0,95--1,00", "1,00--1,05", "1,05--1,10", "1,10--1,20", "$\\ge 1{,}20$"]
    T = ["% gerado por gera_dados.py (la_mef): nao editar a mao\n\\begin{table}[H]\\centering\\footnotesize"
         "\\setlength{\\tabcolsep}{4pt}\n"
         "\\caption{O teste com os dados: $\\mathrm{FS}_\\text{LA}/\\mathrm{FS}_\\text{MEF} - 1$ (em \\%) nas mesmas "
         "realizações, com a KL analítica (os pares do MC direto; o FS da LA foi calculado nas 2400 primeiras e onde "
         "$\\Gamma < 1{,}3$ no $r = -0{,}5$ e $\\Gamma < 1{,}1$ no $r = -0{,}7$) e com a base de Nyström (os pares da ISMC "
         "da semente 21, os 3 centros juntos, com o peso $w$; FS da LA onde $\\Gamma < 1{,}3$). Entre parênteses, os "
         f"pares. No $r = -0{{,}}5$, até 1,05 quase todas as realizações têm o FS da LA (faltam {c05[1]} das {num(c05[0])} "
         "entre 0,95 e 1,05); acima, faltam as de $\\Gamma$ alto, e a média da faixa fica um pouco baixa (entre 1,05 e "
         f"1,10: {sinal(100 * b05b['m'], 2) if b05b else '--'}\\,\\% com todas, {sinal(100 * c05b[2], 2)}\\,\\% só nas "
         f"2400 primeiras, que têm FS em todas). No $r = -0{{,}}7$, o FS da LA falta já entre 1,00 e 1,05 ({c07[1]} das "
         f"{c07[0]}), e a média dessa faixa ({sinal(100 * b07['m'], 2) if b07 else '--'}\\,\\%) também fica um pouco "
         f"baixa: só nas 2400 primeiras, {sinal(100 * c07[2], 2)}\\,\\% ({c07[3]} pares). "
         "$\\Gamma_\\text{LA}$ em $\\mathrm{FS}_\\text{MEF} = 1$: o ajuste de "
         "$\\ln\\Gamma_\\text{LA} = a + k\\ln\\mathrm{FS}_\\text{MEF}$ com $0{,}9 < \\mathrm{FS}_\\text{MEF} < 1{,}1$, "
         "com todas as realizações (e o peso $w$ na ISMC). Embaixo: o campo médio; a linha que passa pelo ponto de projeto do MEF ($\\bm z_\\perp = 0$: um campo "
         "só com os modos do FORM), pela inclinação medida em $c = \\beta \\pm 0{,}1$, contra a média das 100 linhas "
         "sorteadas, as duas em FS ($|s|\\,\\Delta c$, com o $s$ do SRM; direções $M = 20$ / modal 2000); e a grade "
         "da LA ($h = 0{,}1$ e $0{,}025$\\,m nas mesmas realizações: 200 do direto de cada $r$ e 37 da ISMC do centro "
         "$M = 20$, pelos pilotos).}\\label{tab:lateste}\n"
         "\\begin{tabular}{lrrrr}\\toprule\n"
         " & \\multicolumn{2}{c}{KL analítica (MC direto)} & \\multicolumn{2}{c}{Nyström (ISMC s21)} \\\\"
         "\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}\n"
         "$\\mathrm{FS}_\\text{MEF}$ & $r = -0{,}5$ & $r = -0{,}7$ (parcial) & $r = -0{,}5$ & $r = -0{,}7$ \\\\\\midrule\n"]
    for i, rf in enumerate(rot_f):
        T.append(f"{rf} & {cel_b(BAND.get(('an', 'r05', i)))} & {cel_b(BAND.get(('an', 'r07', i)))} & "
                 f"{cel_b(BAND.get(('ny', 'r05', i)))} & {cel_b(BAND.get(('ny', 'r07', i)))} \\\\\n")
    T.append(f"\\textbf{{0,95--1,05}} & \\textbf{{{cel_b(PERTO.get(('an', 'r05')))}}} & "
             f"\\textbf{{{cel_b(PERTO.get(('an', 'r07')))}}} & \\textbf{{{cel_b(PERTO.get(('ny', 'r05')))}}} & "
             f"\\textbf{{{cel_b(PERTO.get(('ny', 'r07')))}}} \\\\\n")
    def cel_g(b):
        return f"{num(b['g1'], 3)} ({b['n']})" if b else "--"

    T.append(f"$\\Gamma_\\text{{LA}}$ em $\\mathrm{{FS}}_\\text{{MEF}} = 1$ (ajuste) & {cel_g(G1.get(('an', 'r05')))} & "
             f"{cel_g(G1.get(('an', 'r07')))} & {cel_g(G1.get(('ny', 'r05')))} & {cel_g(G1.get(('ny', 'r07')))} "
             "\\\\\\midrule\n")
    T.append(f"campo médio (FS {num(FS_DET['LA'], 4)} $\\times$ {num(FS_DET['MEF'], 4)}) & "
             f"\\multicolumn{{4}}{{c}}{{{sinal(100 * (FS_DET['LA'] / FS_DET['MEF'] - 1.), 2)}}} \\\\\n")

    def pto(k, campo):
        a, b = PTO.get((k, "M20")), PTO.get((k, "gim2000"))
        return f"{sinal(100 * a[campo], 2)} / {sinal(100 * b[campo], 2)}" if a and b else "--"

    T.append(f"linha do ponto de projeto ($\\bm z_\\perp = 0$) & -- & -- & {pto('r05', 'dfs')} & {pto('r07', 'dfs')} \\\\\n")
    T.append(f"100 linhas sorteadas (LS s7) & -- & -- & {pto('r05', 'dfsl')} & {pto('r07', 'dfsl')} \\\\\\midrule\n")

    def m_(v):
        return f"{sinal(100 * sum(v) / len(v), 2)} ({len(v)})" if v else "--"

    def g_(k, c):
        return m_(GRr[k][c]) if k in GRr else "--"

    T.append(f"grade $h = 0{{,}}1$ & {g_('r05', 'an01')} & {g_('r07', 'an01')} & {m_(GR['ny01'])} & -- \\\\\n")
    T.append(f"grade $h = 0{{,}}025$ (as mesmas) & {g_('r05', 'an025')} & {g_('r07', 'an025')} & {m_(GR['ny025'])} & "
             "-- \\\\\n")
    T.append(f"$\\mathrm{{FS}}_\\text{{LA}}(0{{,}}1)/\\mathrm{{FS}}_\\text{{LA}}(0{{,}}025) - 1$ & {g_('r05', 'an')} & "
             f"{g_('r07', 'an')} & {m_(GR['ny'])} & -- \\\\\n")
    T.append("\\bottomrule\\end{tabular}\\end{table}\n")
    # a previsao
    T.append("\\begin{table}[H]\\centering\\footnotesize\\setlength{\\tabcolsep}{4pt}\n"
             "\\caption{A $P_f$ da análise limite prevista pelo MEF: $P(\\mathrm{FS}_\\text{MEF}\\,(1 + \\delta) < 1)$ no MC "
             "direto do MEF (24\\,000 e 110\\,000 realizações), com o $\\delta$ perto da falha de cada campo "
             "(Tabela~\\ref{tab:lateste}, linha em negrito), contra a $P_f$ medida da LA: o MC direto com a KL analítica "
             "e a média das 3 ISMC da LA nos centros do MEF (Nyström). $P_f$ em \\%; IC de 95\\,\\% entre colchetes.}"
             "\\label{tab:laprev}\n\\begin{tabular}{lrrrrrrr}\\toprule\n"
             " & & \\multicolumn{3}{c}{com o $\\delta$ da KL analítica} & \\multicolumn{3}{c}{com o $\\delta$ de Nyström} \\\\"
             "\\cmidrule(lr){3-5}\\cmidrule(lr){6-8}\n"
             "$r$ & MEF & $\\delta$ & prevista & medida (direto) & $\\delta$ & prevista & medida (ISMC) \\\\\\midrule\n")
    for k in ("r05", "r07"):
        if (k, "an") not in PREV or (k, "ny") not in PREV:
            continue
        ref = P[k]["ref"]
        oa, on = PREV[(k, "obs_an")], PREV.get((k, "obs_ny"), (float("nan"), float("nan")))
        T.append(f"{RTEX[k]} & {pf_txt(ref['Pf'])} & {sinal(100 * PREV[(k, 'an')]['d'], 2)}\\,\\% & {pf_txt(PREV[(k, 'an')]['p'])} & "
                 f"{pf_txt(oa[0])} [{pf_txt(max(0., oa[0] * (1 - 1.96 * oa[1])))}; {pf_txt(oa[0] * (1 + 1.96 * oa[1]))}]"
                 f"{' (parcial)' if k == 'r07' else ''} & "
                 f"{sinal(100 * PREV[(k, 'ny')]['d'], 2)}\\,\\% & {pf_txt(PREV[(k, 'ny')]['p'])} & {pf_txt(on[0])} "
                 f"[{pf_txt(on[0] * (1 - 1.96 * on[1]))}; {pf_txt(on[0] * (1 + 1.96 * on[1]))}] \\\\\n")
    T.append("\\bottomrule\\end{tabular}\\end{table}\n")
    open(os.path.join(TAB, "la_teste.tex"), "w").write("".join(T))

    # ---------------- tabela: a distribuicao de Gamma da LA nas duas bases, sem FORM nem ISMC (r = -0,5)
    if DG.get("an") and DG.get("ny"):
        def pb(F, n):
            """P em % e o erro padrao binomial"""
            return f"{num(100 * F, 2)} ({num(100 * math.sqrt(F * (1. - F) / n), 2)})"

        def gnum(g):
            return num(g, 3) if math.isfinite(g) else "--"

        a_, b_ = DG["an"], DG["ny"]
        T = ["% gerado por gera_dados.py (la_mef): nao editar a mao\n\\begin{table}[H]\\centering\\footnotesize"
             "\\setlength{\\tabcolsep}{4pt}\n"
             f"\\caption{{A distribuição de $\\Gamma$ da LA nas duas bases, com $r = -0{{,}}5$, sem FORM nem ISMC. KL "
             f"analítica: o MC direto do fermi ({num(a_['N'])} realizações). Nyström, MC simples: {num(b_['N'])} "
             "realizações sem deslocamento (semente 4242, só $\\Gamma$), rodadas no PC em 21/11/2026 para esta revisão. "
             "Nyström, ISMC: a média das 4 ISMC da LA da semente 21, com o peso $w$ (só na cauda, onde estão as "
             "amostras). $P$ em \\%; entre parênteses, o erro padrão binomial (no MC) ou o CoV em \\% (nas ISMC, o CoV "
             "médio, porque as 4 têm a mesma semente); nos quantis, o IC de 95\\,\\% pelas estatísticas de ordem. "
             "Embaixo, o mesmo para o FS do MEF: o MC direto (KL analítica) e as 5 sementes da ISMC de Nyström (a média "
             "dos 3 centros).}\\label{tab:lanys}\n"
             "\\begin{tabular}{lrrrr}\\toprule\n"
             " & \\makecell{KL analítica\\\\MC direto} & \\makecell{Nyström\\\\MC simples} & \\makecell{Nyström\\\\4 ISMC "
             "(com $w$)} & \\makecell{Nyström/analítica\\\\(MC simples; ISMC)} \\\\\\midrule\n"]
        for g0 in NIV_G:
            fa, fb = a_["F"][g0], b_["F"][g0]
            fi = DG.get("is", {}).get(g0)
            rz = (num(fb / fa, 2) if fa > 0 else "--") + ("; " + num(fi[0] / fa, 2) if fi and fa > 0 else "")
            T.append(f"$P(\\Gamma < {num(g0, 2 if round(g0, 1) != g0 else 1)})$ & {pb(fa, a_['N'])} & {pb(fb, b_['N'])} & "
                     f"{(num(100 * fi[0], 2) + ' (' + num(100 * fi[1], 1) + ')') if fi else '--'} & {rz} \\\\\n")
        T.append("\\midrule\n")
        for p in QS:
            qa, qb = a_["q"][p], b_["q"][p]
            rot = "mediana de $\\Gamma$" if p == 0.5 else f"quantil de {num(100 * p, 0)}\\,\\% de $\\Gamma$"
            T.append(f"{rot} & {gnum(qa[0])} [{gnum(qa[1])}; {gnum(qa[2])}] & {gnum(qb[0])} [{gnum(qb[1])}; "
                     f"{gnum(qb[2])}] & -- & {num(qb[0] / qa[0], 3)} \\\\\n")
        T.append(f"média de $\\ln\\Gamma$ (erro padrão) & {num(a_['mlg'], 4)} ({num(a_['elg'], 4)}) & "
                 f"{num(b_['mlg'], 4)} ({num(b_['elg'], 4)}) & -- & $\\times${num(math.exp(b_['mlg'] - a_['mlg']), 3)} "
                 "\\\\\n")
        if DG.get("mef"):
            m_ = DG["mef"]
            T.append("\\midrule\n")
            for s in NIV_FS:
                fa, fi = m_["an"][s], m_["ny"][s]
                T.append(f"MEF: $P(\\mathrm{{FS}} < {num(s, 2)})$ & {pb(fa, m_['Nan'])} & -- & {num(100 * fi[0], 2)} "
                         f"({num(100 * fi[1], 1)}) & {num(fi[0] / fa, 2)} \\\\\n")
        T.append("\\bottomrule\\end{tabular}\\end{table}\n")
        open(os.path.join(TAB, "la_nys.tex"), "w").write("".join(T))

    # ---------------- tabela: os modos de ruptura
    def md(t, n=3):
        return " & ".join(num(v, 1) for v in t[:n]) if t else " & ".join(["--"] * n)

    T = ["% gerado por gera_dados.py (la_mef): nao editar a mao\n\\begin{table}[H]\\centering\\footnotesize"
         "\\setlength{\\tabcolsep}{4pt}\n"
         "\\caption{Os modos de ruptura, em \\% (como na Fig.~20 do Ceron): pelo pé, pela face (acima do pé) e pela base "
         "(abaixo do pé). LA: o mecanismo crítico de cada realização (sem peso na ISMC; no direto, sobram as com "
         "nenhum mecanismo admissível, até 0,2\\,\\%). MEF: pela banda de cisalhamento (verificações de 21/11, V07; a "
         "banda sempre pega o pé, então uma saída a poucos centímetros acima conta como pé).}\\label{tab:lamodos}\n"
         "\\begin{tabular}{llrrrrrrr}\\toprule\n"
         " & & & \\multicolumn{3}{c}{todas as realizações} & \\multicolumn{3}{c}{as que falham} \\\\"
         "\\cmidrule(lr){4-6}\\cmidrule(lr){7-9}\n"
         "$r$ & caso & $N$ (falhas) & pé & face & base & pé & face & base \\\\\\midrule\n"]
    for k in ("r00", "r05", "r07"):
        if k in DIR:
            d = DIR[k]
            T.append(f"{RTEX[k]} & LA, MC direto{' (parcial)' if k == 'r07' else ''} & {num(d['N'])} ({d['fal']}) & "
                     f"{md(d['mod'])} & {md(d['modf'])} \\\\\n")
        if k == "r00":
            T.append(f" & Ceron (2025), Fig.~20 (LA) & {num(CERON_R0['N'])} & {md(CERON_R0['modos'])} & -- & -- & -- \\\\\n")
            continue
        xs = [ISM[(k, c)] for c in ("mef_M20", "mef_M100", "mef_gim20", "la_M20") if (k, c) in ISM]
        if xs:
            mf = tuple(sum(x["modf"][j] for x in xs) / len(xs) for j in range(4))
            T.append(f" & LA, ISMC s21 (Nyström, 4 centros) & {num(sum(x['N'] for x in xs))} & -- & -- & -- & {md(mf)} \\\\\n")
        if k in MEF_MODOS:
            T.append(f" & MEF, MC direto (V07) & {num(P[k]['ref']['N'])} ({P[k]['ref']['falhas']}) & {md(MEF_MODOS[k][0])} & "
                     f"{md(MEF_MODOS[k][1])} \\\\\n")
    T.append("\\bottomrule\\end{tabular}\\end{table}\n")
    open(os.path.join(TAB, "la_modos.tex"), "w").write("".join(T))

    # ---------------- numeros para o texto
    Nn.clear()
    for k, sig in (("r05", "Ra"), ("r07", "Rb"), ("r00", "Z")):
        if k not in DIR:
            continue
        d = DIR[k]
        cmd(f"lmDir{sig}Pf", pct(d["Pf"], 3 if d["Pf"] < 0.01 else 2))
        cmd(f"lmDir{sig}Cov", pct(d["CoV"], 1))
        cmd(f"lmDir{sig}N", num(d["N"]))
        cmd(f"lmDir{sig}Fal", str(d["fal"]))
        cmd(f"lmDir{sig}Lo", pct(max(0., d["Pf"] * (1 - 1.96 * d["CoV"])), 2 if d["Pf"] > 0.01 else 3))
        cmd(f"lmDir{sig}Hi", pct(d["Pf"] * (1 + 1.96 * d["CoV"]), 2 if d["Pf"] > 0.01 else 3))
        cmd(f"lmDir{sig}Gmed", num(d["Gmed"], 3))
        if k in P:
            cmd(f"lmDir{sig}Raz", num(d["Pf"] / P[k]["ref"]["Pf"], 2))
            cmd(f"lmDir{sig}ParN", num(d["Npar"]))
            cmd(f"lmDir{sig}ParFl", str(d["fLpar"]))
            cmd(f"lmDir{sig}ParFm", str(d["fMpar"]))
            cmd(f"lmDir{sig}ParPm", pct(d["pMpar"], 2))
            cmd(f"lmDir{sig}ParRaz", num(d["rpar"], 2))
            cmd(f"lmDir{sig}ParErr", num(d["epar"], 2))
            cmd(f"lmDir{sig}SoL", str(d["soL"]))
        for j, m in enumerate(("Pe", "Face", "Base")):
            cmd(f"lmModo{sig}{m}", num(d["mod"][j], 1))
            if d["modf"]:
                cmd(f"lmModoF{sig}{m}", num(d["modf"][j], 1))
    cmd("lmZRazCeron", num(DIR["r00"]["Pf"] / CERON_R0["Pf"], 2) if "r00" in DIR else "--")
    cmd("lmZz", sinal((DIR["r00"]["Pf"] - CERON_R0["Pf"]) / math.sqrt((DIR["r00"]["Pf"] * DIR["r00"]["CoV"]) ** 2 +
                                                                      (CERON_R0["Pf"] * CERON_R0["CoV"]) ** 2), 1)
        if "r00" in DIR else "--")
    cmd("lmCeronPe", num(CERON_R0["modos"][0], 1))
    cmd("lmCeronFace", num(CERON_R0["modos"][1], 1))
    for k, sig in (("r05", "Ra"), ("r07", "Rb")):
        xs = [ISM[(k, c)] for c in ("mef_M20", "mef_M100", "mef_gim20", "la_M20") if (k, c) in ISM]
        if xs:
            cmd(f"lmIs{sig}Min", pf_txt(min(x["Pf"] for x in xs)))
            cmd(f"lmIs{sig}Max", pf_txt(max(x["Pf"] for x in xs)))
            cmd(f"lmIsCov{sig}Min", num(100 * min(x["CoV"] for x in xs), 1))
            cmd(f"lmIsCov{sig}Max", num(100 * max(x["CoV"] for x in xs), 1))
            r5 = [x["Pf"] / x["ref"]["Pf"] for x in xs if x["ref"]]
            cmd(f"lmIsRazCinco{sig}Min", num(min(r5), 2))
            cmd(f"lmIsRazCinco{sig}Max", num(max(r5), 2))
        pi = [y for y in ism if y["k"] == k and y.get("la")]
        if pi:
            cmd(f"lmIsRaz{sig}Min", num(min(y["rL"] for y in pi), 3))
            cmd(f"lmIsRaz{sig}Max", num(max(y["rL"] for y in pi), 3))
            cmd(f"lmIsErr{sig}Max", num(max(y["eL"] for y in pi), 3))
            cmd(f"lmIsSoM{sig}Min", str(min(y["soM"] for y in pi)))
            cmd(f"lmIsSoM{sig}Max", str(max(y["soM"] for y in pi)))
        ls_ = [LSD[(k, dr)] for dr in ("mef_M20", "mef_gim2000", "la_M20", "la_M100") if (k, dr) in LSD]
        if xs or ls_:
            cmd(f"lmTodos{sig}Min", pf_txt(min(x["Pf"] for x in xs + ls_)))
            cmd(f"lmTodos{sig}Max", pf_txt(max(x["Pf"] for x in xs + ls_)))
        if ls_:
            cmd(f"lmLs{sig}Min", pf_txt(min(x["Pf"] for x in ls_)))
            cmd(f"lmLs{sig}Max", pf_txt(max(x["Pf"] for x in ls_)))
            cmd(f"lmLsNa{sig}Min", num(min(x["na"] for x in ls_), 1))
            cmd(f"lmLsNa{sig}Max", num(max(x["na"] for x in ls_), 1))
            cmd(f"lmLsPerd{sig}", str(sum(x["perd"] for x in ls_)))
        pl = [x for x in lin if x["k"] == k and x.get("la")]
        if pl:
            cmd(f"lmLsRaz{sig}Min", num(min(x["rL"] for x in pl), 3))
            cmd(f"lmLsRaz{sig}Max", num(max(x["rL"] for x in pl), 3))
            cmd(f"lmLsDc{sig}Min", sinal(min(x["mdcL"] for x in pl), 2))
            cmd(f"lmLsDc{sig}Max", sinal(max(x["mdcL"] for x in pl), 2))
            cmd(f"lmLsDfs{sig}Min", num(100 * min(abs(x["sS"]) * x["mdcL"] for x in pl), 2))
            cmd(f"lmLsDfs{sig}Max", num(100 * max(abs(x["sS"]) * x["mdcL"] for x in pl), 2))
            cmd(f"lmSSrm{sig}", num(abs(pl[0]["sS"]), 3))
        for M in (20, 100):
            x = FO.get((k, M))
            if x:
                cmd(f"lmBetaLa{sig}M{'v' if M == 20 else 'c'}", num(x["bL"], 3))
                cmd(f"lmBetaMef{sig}M{'v' if M == 20 else 'c'}", num(x["bM"], 3))
                cmd(f"lmCos{sig}M{'v' if M == 20 else 'c'}", num(x["cos"], 4))
                cmd(f"lmPhiRaz{sig}M{'v' if M == 20 else 'c'}", num(cauda(x["bL"]) / cauda(x["bM"]), 2))
        for f, nome in (("an", "An"), ("ny", "Ny")):
            b = PERTO.get((f, k))
            if b:
                cmd(f"lmDes{nome}{sig}", sinal(100 * b["m"], 2))
                cmd(f"lmDes{nome}{sig}N", num(b["n"]))
                cmd(f"lmDes{nome}{sig}Dp", num(100 * b["dp"], 2))
                cmd(f"lmDes{nome}{sig}Sp", sinal(100 * b["mu"], 2))       # sem peso
            g1 = G1.get((f, k))
            if g1:
                cmd(f"lmGum{nome}{sig}", num(g1["g1"], 3))               # Gamma_LA em FS_MEF = 1 (ajuste)
            if (k, f) in PREV:
                cmd(f"lmPrev{nome}{sig}", pf_txt(PREV[(k, f)]["p"]))
                cmd(f"lmPrevRaz{nome}{sig}", num(PREV[(k, f)]["p"] / P[k]["ref"]["Pf"], 2))
        if (k, "obs_ny") in PREV:
            cmd(f"lmObsNy{sig}", pf_txt(PREV[(k, "obs_ny")][0]))
        a, b = PTO.get((k, "M20")), PTO.get((k, "gim2000"))
        if a and b:
            cmd(f"lmPtoDc{sig}Min", sinal(min(a["dc"], b["dc"]), 3))
            cmd(f"lmPtoDc{sig}Max", sinal(max(a["dc"], b["dc"]), 3))
            cmd(f"lmPtoFs{sig}Max", num(100 * max(a["dfs"], b["dfs"]), 2))
    if PTO:
        cmd("lmPtoFsMax", num(100 * max(v["dfs"] for v in PTO.values()), 2))
        cmd("lmPtoDcMin", sinal(min(v["dc"] for v in PTO.values()), 3))
        cmd("lmPtoDcMax", sinal(max(v["dc"] for v in PTO.values()), 3))
    # as bandas: a faixa do deslocamento em todas as faixas de FS_MEF, por campo
    for f, nome in (("an", "An"), ("ny", "Ny")):
        v = [b["m"] for (ff, k, i), b in BAND.items() if ff == f and b and b["n"] >= 30]
        if v:
            cmd(f"lmBanda{nome}Min", sinal(100 * min(v), 1))
            cmd(f"lmBanda{nome}Max", sinal(100 * max(v), 1))
    cmd("lmDet", sinal(100 * (FS_DET["LA"] / FS_DET["MEF"] - 1.), 2))
    cmd("lmFsDetLa", num(FS_DET["LA"], 4))
    cmd("lmFsDetMef", num(FS_DET["MEF"], 4))
    # FS_LA > FS_MEF em todos os pares com o FS da LA (o direto e as ISMC)?
    dd_ = [x for x in DIR.values() if "acima" in x]
    ok_is = all(fl > fm for x in ism if x.get("la") for fm, fl, g, w in x["parLA"])
    cmd("lmAcima", "todos" if dd_ and all(x["acima"] for x in dd_) and ok_is else "quase todos")
    if GR["an"]:
        cmd("lmGradeAn", sinal(100 * sum(GR["an"]) / len(GR["an"]), 2))
        cmd("lmGradeAnDp", num(100 * mdp(GR["an"])[1], 2))
        cmd("lmGradeAnN", str(len(GR["an"])))
        um = [v for g in GRr.values() for v in g["an01"]]
        dois = [v for g in GRr.values() for v in g["an025"]]
        cmd("lmGradeAnUm", sinal(100 * sum(um) / len(um), 2))
        cmd("lmGradeAnDois", sinal(100 * sum(dois) / len(dois), 2))
    if GR["ny"]:
        cmd("lmGradeNy", sinal(100 * sum(GR["ny"]) / len(GR["ny"]), 2))
        cmd("lmGradeNyDp", num(100 * mdp(GR["ny"])[1], 2))
        cmd("lmGradeNyN", str(len(GR["ny"])))
        cmd("lmGradeNyUm", sinal(100 * sum(GR["ny01"]) / len(GR["ny01"]), 2))
        cmd("lmGradeNyDois", sinal(100 * sum(GR["ny025"]) / len(GR["ny025"]), 2))
        # as faixas de FS_MEF dos pilotos (nao sao as mesmas: o analitico sao as primeiras realizacoes do direto; o de
        # Nystrom, as primeiras da ISMC, perto da falha)
        cmd("lmGradeNyFsMin", num(min(GR["ny_fm"]), 2))
        cmd("lmGradeNyFsMax", num(max(GR["ny_fm"]), 2))
    if GR["an_fm"]:
        cmd("lmGradeAnFsMin", num(min(GR["an_fm"]), 2))
        cmd("lmGradeAnFsMax", num(max(GR["an_fm"]), 2))
    if R0:
        cmd("lmZfator", num(R0["f"], 3))
        cmd("lmZprev", pct(R0["p"], 2))
        cmd("lmZfatorCeron", num(R0["fc"], 3))
    # r = 0 contra o nosso MEF (V12) e a previsao da V12; o Cho (Fig. 9) digitalizado
    if "r00" in DIR:
        d = DIR["r00"]
        rz = d["Pf"] / MEF_R0["Pf"]
        cmd("lmZmef", pct(MEF_R0["Pf"], 2))
        cmd("lmZmefLo", pct(MEF_R0["lo"], 2))
        cmd("lmZmefHi", pct(MEF_R0["hi"], 2))
        cmd("lmZrazMef", num(rz, 2))
        cmd("lmZrazMefErr", num(rz * math.sqrt(d["CoV"] ** 2 + MEF_R0["CoV"] ** 2), 2))
        cmd("lmZprevV", pct(MEF_R0["prev"], 2))
        cmd("lmZprevVLo", pct(MEF_R0["prev_lo"], 2))
        cmd("lmZprevVHi", pct(MEF_R0["prev_hi"], 2))
        cmd("lmZprevVraz", num(MEF_R0["prev"] / MEF_R0["Pf"], 2))
        cmd("lmZchoFig", pct(MEF_R0["cho_fig9"], 2))
    # a distribuicao de Gamma da LA nas duas bases (r = -0,5) e a do FS do MEF
    if DG.get("an") and DG.get("ny"):
        a_, b_ = DG["an"], DG["ny"]
        cmd("lnyN", num(b_["N"]))
        cmd("lnyFal", str(b_["fal"]))
        cmd("lnyPf", pct(b_["Pf"], 2))
        cmd("lnyPfLo", pct(max(0., b_["Pf"] * (1. - 1.96 * b_["CoV"])), 2))
        cmd("lnyPfHi", pct(b_["Pf"] * (1. + 1.96 * b_["CoV"]), 2))
        if DG.get("is"):
            # o MC simples de Nystrom contra a P_f das ISMC de Nystrom: z com o erro binomial do MC simples sob a P_f das
            # ISMC e o erro das ISMC
            fi, ci = DG["is"][1.0]
            cmd("lnyPfIs", pct(fi, 2))
            cmd("lnyZis", sinal((b_["Pf"] - fi) / math.sqrt(fi * (1. - fi) / b_["N"] + (fi * ci) ** 2), 1))
            cmd("lnyEsp", num(fi * b_["N"], 0))      # as falhas que se esperariam com a P_f das ISMC
            # Gamma < 1,05: o MC simples, as ISMC e o direto analitico
            fb, cb = DG["is"][1.05]
            cmd("lnyFalB", str(round(b_["F"][1.05] * b_["N"])))
            cmd("lnyPfB", pct(b_["F"][1.05], 2))
            cmd("lnyPfIsB", pct(fb, 2))
            cmd("lanPfB", pct(a_["F"][1.05], 2))
            cmd("lnyZisB", sinal((b_["F"][1.05] - fb) / math.sqrt(fb * (1. - fb) / b_["N"] + (fb * cb) ** 2), 1))
        cauda_ = [g for g in NIV_G if g >= 1.1]          # a cauda com bastante realizacoes no MC simples
        rr = [b_["F"][g] / a_["F"][g] for g in cauda_]
        cmd("lnyRazMin", num(min(rr), 2))
        cmd("lnyRazMax", num(max(rr), 2))
        zz = [(b_["F"][g] - a_["F"][g]) / math.sqrt(a_["F"][g] * (1 - a_["F"][g]) / a_["N"] + b_["F"][g] * (1 - b_["F"][g])
                                                    / b_["N"]) for g in cauda_]
        cmd("lnyZmin", num(min(zz), 1))
        cmd("lnyZmax", num(max(zz), 1))
        if DG.get("is"):
            ri = [DG["is"][g][0] / a_["F"][g] for g in cauda_]
            cmd("lisRazMin", num(min(ri), 2))
            cmd("lisRazMax", num(max(ri), 2))
            # o MC simples de Nystrom contra as ISMC de Nystrom (a mesma base): a razao das P(Gamma < g)
            rb = [b_["F"][g] / DG["is"][g][0] for g in cauda_]
            cmd("lnyIsMin", num(min(rb), 2))
            cmd("lnyIsMax", num(max(rb), 2))
        cmd("lnyMed", sinal(100 * (b_["q"][0.5][0] / a_["q"][0.5][0] - 1.), 1))
        cmd("lnyQvc", sinal(100 * (b_["q"][0.25][0] / a_["q"][0.25][0] - 1.), 1))      # o quantil de 25 %
        cmd("lnyQcinco", sinal(100 * (b_["q"][0.05][0] / a_["q"][0.05][0] - 1.), 1))  # o de 5 %
        cmd("lnyQum", sinal(100 * (b_["q"][0.01][0] / a_["q"][0.01][0] - 1.), 1))     # o de 1 %
        cmd("lnyLg", sinal(100 * (math.exp(b_["mlg"] - a_["mlg"]) - 1.), 1))
        cmd("lnyLgErr", num(100 * math.sqrt(a_["elg"] ** 2 + b_["elg"] ** 2), 1))
        cmd("lnyLgZ", num(abs(b_["mlg"] - a_["mlg"]) / math.sqrt(a_["elg"] ** 2 + b_["elg"] ** 2), 1))
        if DG.get("mef"):
            rm = [DG["mef"]["ny"][s][0] / DG["mef"]["an"][s] for s in NIV_FS]
            cmd("lmefRazMin", num(min(rm), 2))
            cmd("lmefRazMax", num(max(rm), 2))
    open(os.path.join(AQUI, "numeros.tex"), "a").write("".join(Nn))

    # ---------------- figura: (a, b) o deslocamento nos pares, (c) a P_f prevista contra o delta, (d) a grade
    fig, ax = plt.subplots(2, 2, figsize=(LARG, 6.6))
    cor = {"an": "#7f7f7f", "ny": "#d95f02"}
    for a, k, t in ((ax[0][0], "r05", "(a)"), (ax[0][1], "r07", "(b)")):
        for f, rot in (("an", "KL analítica (direto)"), ("ny", "Nyström (ISMC s21)")):
            pp = [(fm, 100 * (fl / fm - 1.)) for fm, fl, g, w in PP.get((f, k), []) if 0.85 < fm < 1.45]
            if not pp:
                continue
            a.plot([x for x, _ in pp], [y for _, y in pp], ".", color=cor[f], ms=1.6, alpha=0.35, mew=0)
            xs, ys = [], []
            for i, (lo, hi) in enumerate(FAIXAS):
                b = BAND.get((f, k, i))
                if b and b["n"] >= 30:
                    xs.append(0.5 * (max(lo, 0.85) + min(hi, 1.45)))
                    ys.append(100 * b["m"])
            a.plot(xs, ys, "o-", color=cor[f], lw=1.6, ms=4.5,
                   label=rot + (f": {100 * PERTO[(f, k)]['m']:+.2f} % perto de 1".replace(".", ",") if PERTO.get((f, k)) else ""))
        a.plot([FS_DET["MEF"]], [100 * (FS_DET["LA"] / FS_DET["MEF"] - 1.)], "*", color="k", ms=10,
               label="campo médio (determinístico)")
        p0 = PTO.get((k, "M20"))
        if p0:
            a.plot([1.0], [100 * p0["dfs"]], "D", color="#1b9e77", ms=6, label="linha do ponto de projeto (LS, $M = 20$)")
        a.axhline(0, color="k", lw=0.7)
        a.axvline(1, color="0.6", ls=":", lw=0.8)
        a.set_xlim(0.85, 1.45)
        a.set_ylim(-1, 6)
        a.set_xlabel(r"$\mathrm{FS}_\mathrm{MEF}$")
        a.set_ylabel(r"$\mathrm{FS}_\mathrm{LA}/\mathrm{FS}_\mathrm{MEF} - 1$ (%)")
        a.set_title(f"{t} " + tit(k) + (" (direto parcial)" if k == "r07" else ""), fontsize=10, loc="left")
        a.legend(loc="upper left", fontsize=7.5, framealpha=0.9)
    a = ax[1][0]
    for k, c in (("r05", "#2166ac"), ("r07", "#b2182b")):
        if k not in DIR or "fsM" not in DIR[k]:
            continue
        fsM = DIR[k]["fsM"]
        dls = [j * 0.0005 for j in range(0, 81)]
        a.semilogy([100 * x for x in dls], [100 * max(sum(1 for v in fsM if v * (1 + x) < 1.) / len(fsM), 1e-5) for x in dls],
                   "-", color=c, lw=1.5, label=tit(k) + ": MEF deslocado")
        for f, mk in (("an", "s"), ("ny", "o")):
            ob = PREV.get((k, "obs_an" if f == "an" else "obs_ny"))
            if ob and PERTO.get((f, k)):
                x = 100 * PERTO[(f, k)]["m"]
                a.errorbar([x], [100 * ob[0]], yerr=[[100 * ob[0] * min(0.999, 1.96 * ob[1])], [100 * ob[0] * 1.96 * ob[1]]],
                           fmt=mk, color=c, ms=6, capsize=3, mfc="white" if f == "an" else c)
    a.plot([], [], "s", color="0.3", mfc="white", label="LA medida, KL analítica (direto)")
    a.plot([], [], "o", color="0.3", label="LA medida, Nyström (ISMC)")
    a.set_xlabel(r"$\delta$ (%): $\mathrm{FS}_\mathrm{LA} = \mathrm{FS}_\mathrm{MEF}\,(1 + \delta)$")
    a.set_ylabel(r"$P(\mathrm{FS}_\mathrm{MEF}\,(1 + \delta) < 1)$ (%)")
    a.set_title("(c) a $P_f$ da LA prevista pelo MEF", fontsize=10, loc="left")
    a.legend(loc="lower left", fontsize=7.5)
    a.set_xlim(0, 4)
    a.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}".replace(".", ",")))
    a = ax[1][1]
    if GR["an"] and GR["ny"]:
        bins = [-1.6 + 0.05 * j for j in range(33)]
        a.hist([100 * v for v in GR["an"]], bins=bins, histtype="step", lw=1.6, color=cor["an"], density=True,
               label=f"KL analítica ({len(GR['an'])}): média {100 * sum(GR['an']) / len(GR['an']):+.2f} %".replace(".", ","))
        a.hist([100 * v for v in GR["ny"]], bins=bins, histtype="step", lw=1.6, color=cor["ny"], density=True,
               label=f"Nyström ({len(GR['ny'])}): média {100 * sum(GR['ny']) / len(GR['ny']):+.2f} %".replace(".", ","))
        a.axvline(0, color="k", lw=0.7)
        a.set_xlabel(r"$\mathrm{FS}_\mathrm{LA}(h = 0{,}1)/\mathrm{FS}_\mathrm{LA}(h = 0{,}025) - 1$ (%)")
        a.set_ylabel("densidade")
        a.set_title("(d) a grade da LA, $h$: 0,1 x 0,025 m", fontsize=10, loc="left")
        a.legend(loc="upper left", fontsize=7.5)
    for a in (ax[0][0], ax[0][1], ax[1][1]):
        a.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}".replace(".", ",")))
    fig.tight_layout()
    fig.savefig(os.path.join(FIGS, "la_mef.pdf"))
    plt.close(fig)
    for k in DIR:
        d = DIR[k]
        print(f"LA direto {k}: N {d['N']} falhas {d['fal']} P_f {100 * d['Pf']:.4f}% (CoV {100 * d['CoV']:.1f}%)"
              + (f"; pareado {d['Npar']}: LA {d['fLpar']} MEF {d['fMpar']} razao {d['rpar']:.3f} +- {d['epar']:.3f}" if "Npar" in d else ""))
    for (f, k), b in PERTO.items():
        if b:
            print(f"LA desloc perto da falha {f} {k}: n {b['n']} com peso {100 * b['m']:+.2f}% sem peso {100 * b['mu']:+.2f}% "
                  f"dp {100 * b['dp']:.2f}% ESS {b['ess']:.0f}")
    for kk, v in PREV.items():
        print("LA prev", kk, v)
    print("LA ponto de projeto", {kk: {a: round(b, 4) for a, b in v.items()} for kk, v in PTO.items()})
    print("LA grade", {kk: (len(v), round(100 * sum(v) / len(v), 3) if v else None) for kk, v in GR.items()})
    print("LA r = 0", R0)


ISM_, LIN_ = srm_gim()
ceron_agua()
la_mef(ISM_, LIN_)
