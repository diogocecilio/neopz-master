#!/usr/bin/env python3
"""gera_dados.py: as tabelas, os numeros e as figuras do relatorio_ceron_secao6.tex (Secoes 6.1 a 6.5 de Vargas Ceron
et al. 2025 com a agua do artigo), a partir da campanha lsCeronC baixada pelo baixar_ls_ceron.sh:

  - <campanha>/resumo.json (resumo_ls_ceron.py): o deterministico, o FORM modal, o line sampling e o MC pareado de cada
    caso, com os numeros das Tabelas 3 a 6 do artigo;
  - <campanha>/figs/*.pdf (figuras_ls_ceron.py): as Figs. 25 e 27 a 33 do artigo refeitas, copiadas para figs/;
  - <campanha>/ref/mc/mc_*.csv: as realizacoes do MC da referencia, para o Spearman parcial (Fig. 26);
  - <campanha>/<caso>/ls/linhas.csv e form/form.json: as raizes das linhas (fig_raizes);
  - o relatorio antigo de 13/11 (agua permanente antiga, Gamma 1,403), para a coluna "agua antiga":
    ~/projects/runs/_arquivo/2026-11-11_ceron_agua_antiga/fermi_ceron/relatorio_parcial/estatisticas.json.

Grava tab/*.tex, numeros.tex e figs/ (fig26.pdf, fig_z.pdf, fig_raizes.pdf e as copias). Roda em segundos.

    python3 -I gera_dados.py [pasta da campanha = ~/projects/runs/fermi_ceron_ls_casos]"""

import csv
import glob
import json
import math
import os
import shutil
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from scipy import stats  # noqa: E402

AQUI = os.path.dirname(os.path.abspath(__file__))
RAIZ = os.path.abspath(os.path.join(AQUI, "..", "..", ".."))
sys.path.insert(0, os.path.join(RAIZ, "Projects", "SlopeRandomField", "comparacao_ceron"))
import ceron_dig  # noqa: E402
import dados_sat  # noqa: E402
import figuras_sat as fs  # noqa: E402  (o estilo das figuras do artigo refeitas)

CAMPANHA = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/projects/runs/fermi_ceron_ls_casos")
ANTIGA = os.path.expanduser("~/projects/runs/_arquivo/2026-11-11_ceron_agua_antiga/fermi_ceron/relatorio_parcial/"
                            "estatisticas.json")
FIGS = os.path.join(AQUI, "figs")
# Fig. 26 do artigo (Spearman parcial de 2a ordem), lida nos pixels da figura do PDF (pagina 2481): as barras vao de 0
# a +0,879 (c), +0,509 (phi) e -0,036 (k_v), na escala de -1 a 1 da moldura
FIG26_CERON = {"c": 0.879, "phi": 0.509, "k": -0.036}
FAMILIAS = [("covk", "6.2", r"$CoV(k_v)$"), ("covc", "6.2", r"$CoV(c)$"), ("covphi", "6.2", r"$CoV(\varphi)$"),
            ("s", "6.3", r"$s$"), ("alfa", "6.4", r"$\alpha$"), ("hw", "6.5", r"$h_w/H$")]
VALOR_REF = fs.VALOR_REF
NOMES_CASO = {"alfa=2": "AlfaII", "alfa=3": "AlfaIII", "alfa=4": "AlfaIV", "alfa=5": "AlfaV", "s=400": "SCD",
              "covk=0": "CovkZero", "covc=0.1": "CovcUm", "hw=0.9": "HwIX"}


def br(x, casas=2):
    """numero com virgula decimal para o LaTeX"""
    if x is None or (isinstance(x, float) and not math.isfinite(x)):
        return "--"
    t = f"{x:.{casas}f}".replace(".", "{,}")
    return "$-$" + t[1:] if t.startswith("-") else t  # o menos de verdade, em modo texto


def pct(x, casas=2):
    return br(None if x is None else 100. * x, casas)


def valor_txt(fam, v):
    if fam in ("covk", "covc", "covphi"):
        return br(v, 2 if fam == "covphi" else 1 if v * 10 == int(v * 10) else 2)
    if fam == "s":
        return br(v, 1) if v != int(v) else str(int(v))
    if fam == "alfa":
        return str(int(v))
    return br(v, 1)


def clopper_pearson(k, n):
    if not n:
        return (math.nan, math.nan)
    lo = 0. if k == 0 else stats.beta.ppf(0.025, k, n - k + 1)
    hi = 1. if k == n else stats.beta.ppf(0.975, k + 1, n - k)
    return lo, hi


def carrega():
    with open(os.path.join(CAMPANHA, "resumo.json")) as f:
        casos = json.load(f)
    por_id = {c["id"]: c for c in casos}
    antiga = {}
    if os.path.isfile(ANTIGA):
        with open(ANTIGA) as f:
            for c in json.load(f)["casos"]:
                antiga[c["id"]] = c
    return casos, por_id, antiga


def caso_de(por_id, fam, x):
    if abs(x - VALOR_REF[fam]) < 1e-9:
        return por_id.get("ref")
    return por_id.get(f"{fam}={x:g}")


def mc_txt(mc):
    if not mc or not mc.get("S"):
        return "--"
    return f"{pct(mc.get('P_f'))} ({mc.get('n_f', 0)}/{mc['S']})"


def tabela_pf(fam, sec, nome, tab, por_id, antiga):
    """P_f: o artigo, o line sampling, o FORM, o MC pareado e a agua antiga (13/11)"""
    lin = [r"\begin{table}[htbp]\centering\scriptsize\setlength{\tabcolsep}{3pt}",
           rf"\caption{{Seção {sec}, {nome}: $P_f$ (\%). Ceron: Tabelas 3 a 6 do artigo, com o $CoV(P_f)$. LS: line "
           r"sampling com o IC de 95\,\% e o $CoV$; $z = (P_f^{LS} - P_f^{C})/\sqrt{\sigma_{LS}^2 + \sigma_C^2}$. "
           r"FORM: $\beta$ e $\Phi(-\beta)$. MC: Monte Carlo direto pareado, falhas/realizações (parcial). Água "
           r"antiga: o MC do relatório de 13/11.}" + f"\\label{{tab:pf_{fam}}}",
           r"\begin{tabular}{rrrrrrrrr}\toprule",
           rf"{nome} & Ceron & LS [IC 95\,\%] & $CoV$ & linhas/análises & $z$ & FORM: $\beta$; $P_f$ & MC & "
           r"água antiga \\\midrule"]
    for r in tab:
        c = caso_de(por_id, fam, r["x"])
        if c is None:
            continue
        ls, fo, mc = c["ls"], c["form"], c["mc"]
        ic = ls.get("IC95") or [None, None]
        ant = antiga.get(c["id"], {}).get("estat", {})
        ref = " (ref.)" if c["id"] == "ref" else ""
        lin.append(
            f"{valor_txt(fam, r['x'])}{ref} & {pct(r['pf'])} ({br(100 * c['cov_pf_artigo'], 1)}) & "
            f"{pct(ls.get('P_f'))} [{pct(ic[0])}; {pct(ic[1])}] & {pct(ls.get('CoV'), 1)} & "
            f"{ls.get('linhas', '--')}/{ls.get('analises', '--')} & "
            f"{'--' if c['z'] is None else br(c['z'], 1)} & "
            f"{br(fo.get('beta'), 3)}; {pct(fo.get('P_f'))} & {mc_txt(mc)} & {pct(ant.get('P_f'))} \\\\")
    lin += [r"\bottomrule\end{tabular}\end{table}", ""]
    return "\n".join(lin)


def tabela_gamma(fam, sec, nome, tab, por_id, antiga):
    """a distribuicao de Gamma: o deterministico, a media e o CoV (artigo x MC pareado x agua antiga)"""
    tem_det = fam in ("alfa", "hw")
    cab = [nome]
    if tem_det:
        cab += [r"$\Gamma_{det}$ C", r"$\Gamma_{det}$ T", r"$\Gamma_{det}$ antiga"]
    cab += [r"$\mu$ C", r"$\mu$ MC", r"$\mu$ antiga", r"$CoV$ C", r"$CoV$ MC", r"$CoV$ antiga"]
    lin = [r"\begin{table}[htbp]\centering\scriptsize\setlength{\tabcolsep}{3pt}",
           rf"\caption{{Seção {sec}, {nome}: $\Gamma$. C: Ceron (Tabelas 3 a 6). T: este trabalho (o determinístico "
           r"com as médias). MC: média e $CoV$ amostrais do Monte Carlo pareado (parcial). Antiga: o relatório de "
           r"13/11 (amostrais). $CoV$ em \%.}" + f"\\label{{tab:g_{fam}}}",
           r"\begin{tabular}{" + "r" * len(cab) + r"}\toprule", " & ".join(cab) + r" \\\midrule"]
    for r in tab:
        c = caso_de(por_id, fam, r["x"])
        if c is None:
            continue
        mc = c["mc"]
        ant_c = antiga.get(c["id"], {})
        ant = ant_c.get("estat", {})
        cel = [valor_txt(fam, r["x"]) + (" (ref.)" if c["id"] == "ref" else "")]
        if tem_det:
            cel += [br(r.get("det"), 3), br((c.get("det") or {}).get("Gamma"), 3),
                    br((ant_c.get("det") or {}).get("Gamma"), 3)]
        cel += [br(r["mu"], 3), br(mc.get("media"), 3), br(ant.get("media_amostral"), 3), br(100 * r["cov"], 1),
                br(None if mc.get("CoV") is None else 100 * mc["CoV"], 1),
                br(None if ant.get("CoV_amostral") is None else 100 * ant["CoV_amostral"], 1)]
        lin.append(" & ".join(cel) + r" \\")
    lin += [r"\bottomrule\end{tabular}\end{table}", ""]
    return "\n".join(lin)


def tabela_ref(c, ant):
    """a referencia (Secao 6.1): os metodos lado a lado"""
    ls, fo, mc = c["ls"], c["form"], c["mc"]
    lo, hi = clopper_pearson(mc.get("n_f", 0), mc.get("S", 0))
    a = ant.get("estat", {})
    ic = ls["IC95"]
    lin = [r"\begin{table}[htbp]\centering\footnotesize",
           r"\caption{Referência (Seção 6.1). Ceron: $S = 10\,000$ simulações pela análise limite cinemática. Este "
           r"trabalho: o MEF (GIM) com a água do artigo, pelo FORM modal, pelo line sampling (LS) e pelo Monte Carlo "
           r"direto pareado (MC, parcial; IC exato de Clopper--Pearson~\cite{clopper1934}). Água antiga: o MC do "
           r"relatório de 13/11.}"
           r"\label{tab:ref}",
           r"\begin{tabular}{lrrrrr}\toprule",
           r"& Ceron & FORM modal & LS & MC & água antiga \\\midrule",
           rf"$\Gamma_{{det}}$ & 1{{,}}336 & \multicolumn{{3}}{{c}}{{{br(c['det']['Gamma'], 4)}}} & "
           rf"{br((ant.get('det') or {}).get('Gamma'), 3)} \\",
           rf"$P_f$ (\%) & 11{{,}}50 & {pct(fo['P_f'])} & {pct(ls['P_f'])} & {pct(mc.get('P_f'))} & {pct(a.get('P_f'))} \\",
           rf"IC 95\,\% (\%) & -- & -- & [{pct(ic[0])}; {pct(ic[1])}] & [{pct(lo)}; {pct(hi)}] & -- \\",
           rf"$CoV(P_f)$ (\%) & 2{{,}}77 & -- & {pct(ls['CoV'], 1)} & {pct(mc.get('CoV_Pf'), 1)} & "
           rf"{pct(a.get('CoV_Pf'), 1)} \\",
           rf"$\beta$ & -- & {br(fo['beta'], 3)} & {br(-stats.norm.ppf(ls['P_f']), 3)} & -- & -- \\",
           rf"análises do MEF & -- & {fo['analises']} & {ls['analises']} ({ls['linhas']} linhas) & {mc.get('S')} & "
           rf"{a.get('S', '--')} \\",
           rf"$\mu(\Gamma)$ & 1{{,}}353 & -- & -- & {br(mc.get('media'), 3)} & {br(a.get('media_amostral'), 3)} \\",
           rf"$CoV(\Gamma)$ (\%) & 23{{,}}5 & -- & -- & {pct(mc.get('CoV'), 1)} & {pct(a.get('CoV_amostral'), 1)} \\",
           r"\bottomrule\end{tabular}\end{table}", ""]
    return "\n".join(lin)


def spearman_ref():
    """Spearman parcial de 2a ordem na referencia (dados_sat.spearman, as mesmas colunas do relatorio antigo)"""
    linhas = []
    for arq in sorted(glob.glob(os.path.join(CAMPANHA, "ref", "mc", "mc_*_*.csv"))):
        if os.path.basename(arq).startswith("mc_mapa_"):
            continue
        with open(arq, newline="") as f:
            for r in csv.DictReader(f):
                try:
                    g = float(r["FS_min"])
                except (TypeError, ValueError):
                    continue
                if r.get("situacao") == "ok" and math.isfinite(g) and g > 0:
                    linhas.append(r)
    g = np.array([float(r["FS_min"]) for r in linhas])
    return dados_sat.spearman(linhas, g, "banda"), dados_sat.spearman(linhas, g, "campo"), len(linhas)


def tabela_spearman(sb, sc, n):
    nomes = [("k", r"$k_v$"), ("phi", r"$\varphi$"), ("c", r"$c$")]
    lin = [r"\begin{table}[htbp]\centering\footnotesize",
           r"\caption{Spearman parcial de 2ª ordem entre $\Gamma$ e cada variável, controlando as outras duas (Eqs.~70 "
           r"e 71 do artigo), com o IC de 95\,\% (transformação de Fisher), nas " + str(n) + r" realizações do MC da "
           r"referência. Ceron: lido nas barras da Fig.~26. Este trabalho, duas definições das variáveis de cada "
           r"realização: $c'$ e $\varphi'$ médios na banda de ruptura e $k$ médio na massa que desliza; e as médias de "
           r"cada campo no solo todo.}\label{tab:spearman}",
           r"\begin{tabular}{lrrr}\toprule",
           r"& Ceron (Fig.~26) & banda e massa & médias do campo \\\midrule"]
    for k, nome in nomes:
        def cel(s):
            if not s:
                return "--"
            v = s[k]
            return f"{br(v['rho_parcial'], 3)} [{br(v['ic95'][0], 2)}; {br(v['ic95'][1], 2)}]"
        lin.append(f"{nome} & {br(FIG26_CERON[k], 2)} & {cel(sb)} & {cel(sc)} \\\\")
    lin += [r"\bottomrule\end{tabular}\end{table}", ""]
    return "\n".join(lin)


def figura_z(casos):
    ordem = [c for c in casos if c["z"] is not None]
    fig, ax = plt.subplots(figsize=(6.4, 4.6), constrained_layout=True)
    cores = {"6.1": fs.PALETA[0], "6.2": fs.PALETA[1], "6.3": fs.PALETA[2], "6.4": fs.PALETA[3], "6.5": fs.PALETA[4]}
    y = np.arange(len(ordem))[::-1]
    ax.axvspan(-2, 2, color=fs.GRADE, alpha=0.6, lw=0)
    ax.barh(y, [c["z"] for c in ordem], color=[cores.get(c["secao"], fs.TINTA2) for c in ordem], height=0.7)
    ax.axvline(0, color=fs.TINTA2, lw=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels([c["id"] for c in ordem], fontsize=7)
    ax.set_xlabel(r"$z = (P_f^{LS} - P_f^{Ceron}) / \sqrt{\sigma_{LS}^2 + \sigma_{Ceron}^2}$")
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color=v, label=f"Seção {k}") for k, v in cores.items()] +
              [Patch(color=fs.GRADE, label=r"$|z| < 2$")], loc="lower right")
    arq = os.path.join(FIGS, "fig_z.pdf")
    fig.savefig(arq)
    plt.close(fig)


def figura_raizes(por_id):
    ids = ["ref", "alfa=2", "alfa=3", "alfa=4", "alfa=5", "s=400"]
    dados, betas, bls = [], [], []
    for i in ids:
        with open(os.path.join(CAMPANHA, i, "ls", "linhas.csv"), newline="") as f:
            cs = [float(r["c"]) for r in csv.DictReader(f) if r["c"] not in ("", "inf", "nan")]
        dados.append(cs)
        betas.append(por_id[i]["form"]["beta"])
        bls.append(-stats.norm.ppf(por_id[i]["ls"]["P_f"]))
    fig, ax = plt.subplots(figsize=(6.4, 3.0), constrained_layout=True)
    pos = np.arange(len(ids))
    ax.boxplot(dados, positions=pos, widths=0.5, showfliers=True, flierprops={"markersize": 2},
               medianprops={"color": fs.TINTA}, boxprops={"color": fs.TINTA2}, whiskerprops={"color": fs.TINTA2},
               capprops={"color": fs.TINTA2})
    ax.plot(pos, betas, "D", color=fs.PALETA[1], ms=5, label=r"$\beta$ do FORM")
    ax.plot(pos, bls, "x", color=fs.PALETA[0], ms=7, mew=1.5, label=r"$-\Phi^{-1}(P_f^{LS})$")
    ax.set_xticks(pos)
    ax.set_xticklabels(ids)
    ax.set_ylabel(r"raiz $c$ de cada linha")
    ax.legend(loc="upper left")
    arq = os.path.join(FIGS, "fig_raizes.pdf")
    fig.savefig(arq)
    plt.close(fig)


def numeros(casos, por_id, antiga, sb, n_sp, tab):
    ref = por_id["ref"]
    zs = [c["z"] for c in casos if c["z"] is not None]
    n_mc = [c["mc"].get("S", 0) for c in casos]
    tot_form = sum(c["form"].get("analises") or 0 for c in casos)
    tot_ls = sum(c["ls"].get("analises") or 0 for c in casos)
    tot_mc = sum(n_mc)
    fora = [c for c in casos if c["z"] is not None and abs(c["z"]) >= 2]
    d = {
        "nCasos": len(casos), "nZdois": sum(1 for z in zs if abs(z) < 2), "nZum": sum(1 for z in zs if abs(z) < 1),
        "casosFora": ", ".join(f"{c['id']} ({br(c['z'], 1)})" for c in fora),
        "nMCmin": min(n_mc), "nMCmax": max(n_mc), "nMCtot": tot_mc, "nFORMtot": tot_form, "nLStot": tot_ls,
        "nTot": tot_form + tot_ls + tot_mc,
        "linhasMin": min(c["ls"]["linhas"] for c in casos), "linhasMax": max(c["ls"]["linhas"] for c in casos),
        "anLSmin": min(c["ls"]["analises"] for c in casos), "anLSmax": max(c["ls"]["analises"] for c in casos),
        "anPorLinha": br(tot_ls / sum(c["ls"]["linhas"] for c in casos), 2),
        "nLinhasMax": sum(1 for c in casos if c["ls"]["linhas"] >= 600),
        "formIt": br(np.mean([c["form"]["iteracoes"] for c in casos]), 1),
        "formAnMin": min(c["form"]["analises"] for c in casos), "formAnMax": max(c["form"]["analises"] for c in casos),
        "gdetRef": br(ref["det"]["Gamma"], 4), "pfLSref": pct(ref["ls"]["P_f"]), "cvLSref": pct(ref["ls"]["CoV"], 1),
        "pfFORMref": pct(ref["form"]["P_f"]), "betaFORMref": br(ref["form"]["beta"], 3),
        "pfMCref": pct(ref["mc"].get("P_f")), "nfMCref": ref["mc"].get("n_f", 0), "nMCref": ref["mc"].get("S", 0),
        "muMCref": br(ref["mc"].get("media"), 3), "cvMCref": pct(ref["mc"].get("CoV"), 1),
        "pfAntigaRef": pct(antiga.get("ref", {}).get("estat", {}).get("P_f")),
        "nSpearman": n_sp,
    }
    for i, suf in NOMES_CASO.items():
        c = por_id[i]
        d[f"pfLS{suf}"] = pct(c["ls"]["P_f"])
        d[f"pfFORM{suf}"] = pct(c["form"]["P_f"])
        d[f"pfC{suf}"] = pct(c["pf_artigo"])
        d[f"z{suf}"] = br(c["z"], 1)
        d[f"cvLS{suf}"] = pct(c["ls"]["CoV"], 1)
        d[f"cvMC{suf}"] = pct(c["mc"].get("CoV"), 1)
        d[f"gdet{suf}"] = br(c["det"]["Gamma"], 4)
    if sb:
        for k in ("c", "phi", "k"):
            d[f"sp{k}"] = br(sb[k]["rho_parcial"], 2)
    # a media do MC pareado contra a do artigo (%), por familia e no total, e contra o deterministico na referencia
    def sinal(x):
        return ("+" if x > 0 else "") + br(x, 1)
    difs = {}
    for fam, _, _ in FAMILIAS:
        for r in tab[fam]:
            c = caso_de(por_id, fam, r["x"])
            if c is not None and c["mc"].get("media"):
                difs.setdefault(fam, []).append(100. * (c["mc"]["media"] / r["mu"] - 1.))
    nomes_fam = {"covk": "Covk", "covc": "Covc", "covphi": "Covphi", "s": "S", "alfa": "Alfa", "hw": "Hw"}
    for fam, v in difs.items():
        d[f"muDifMin{nomes_fam[fam]}"] = sinal(min(v))
        d[f"muDifMax{nomes_fam[fam]}"] = sinal(max(v))
    covs = difs.get("covk", []) + difs.get("covc", []) + difs.get("covphi", [])
    todas = [x for v in difs.values() for x in v]
    d["muDifMinCovs"], d["muDifMaxCovs"] = sinal(min(covs)), sinal(max(covs))
    d["muDifMinAll"], d["muDifMaxAll"] = sinal(min(todas)), sinal(max(todas))
    unicos = {}
    for fam, _, _ in FAMILIAS:
        for r in tab[fam]:
            c = caso_de(por_id, fam, r["x"])
            if c is not None and c["mc"].get("media"):
                unicos[c["id"]] = 100. * (c["mc"]["media"] / r["mu"] - 1.)
    d["nMuAbaixo"], d["nMuTotal"] = sum(1 for x in unicos.values() if x < 0), len(unicos)
    d["muDifRef"] = sinal(100. * (ref["mc"]["media"] / 1.353 - 1.))
    d["muVsDetRef"] = sinal(100. * (ref["mc"]["media"] / ref["det"]["Gamma"] - 1.))
    d["cvMCCovphiII"] = pct(por_id["covphi=0.2"]["mc"].get("CoV"), 1)
    # erro padrao do CoV amostral (em pontos), com CoV ~0,27 e o menor N do MC
    d["cvErroPts"] = br(100. * 0.27 / math.sqrt(2. * min(n_mc)), 1)
    lin = ["% gerado por gera_dados.py; nao edite"]
    for k, v in d.items():
        assert k.isalpha(), k  # macros do LaTeX: so letras
        lin.append(rf"\newcommand{{\{k}}}{{{v}}}")
    return "\n".join(lin) + "\n", d


def main():
    os.makedirs(FIGS, exist_ok=True)
    for arq in glob.glob(os.path.join(CAMPANHA, "figs", "*.pdf")):
        shutil.copy2(arq, FIGS)
    casos, por_id, antiga = carrega()
    tab = ceron_dig.tabelas()
    sb, sc, n_sp = spearman_ref()
    os.makedirs(os.path.join(AQUI, "tab"), exist_ok=True)
    blocos = {"ref": tabela_ref(por_id["ref"], antiga.get("ref", {})), "spearman": tabela_spearman(sb, sc, n_sp)}
    for fam, sec, nome in FAMILIAS:
        blocos[f"pf_{fam}"] = tabela_pf(fam, sec, nome, tab[fam], por_id, antiga)
        blocos[f"g_{fam}"] = tabela_gamma(fam, sec, nome, tab[fam], por_id, antiga)
    for nome, texto in blocos.items():
        with open(os.path.join(AQUI, "tab", nome + ".tex"), "w") as f:
            f.write("% gerado por gera_dados.py; nao edite\n" + texto)
    texto, d = numeros(casos, por_id, antiga, sb, n_sp, tab)
    with open(os.path.join(AQUI, "numeros.tex"), "w") as f:
        f.write(texto)
    est = {"casos": [{"familia": "ref", "spearman_banda": sb, "spearman_campo": sc}]}
    cer = ceron_dig.carrega(None)
    cer["fig26"] = FIG26_CERON
    fs.figura26(est, cer, FIGS)
    figura_z(casos)
    figura_raizes(por_id)
    print(f"{len(casos)} casos; |z| < 2 em {d['nZdois']}; MC de {d['nMCmin']} a {d['nMCmax']} por caso; "
          f"{d['nTot']} analises do MEF; Spearman com {n_sp} realizacoes")


if __name__ == "__main__":
    main()
