#!/usr/bin/env python3
"""Relatório comparativo em LaTeX a partir de analise_artigo.py.

Uso: relatorio_tex.py <saída de analise_artigo.py> <diretório do relatório>
Escreve <dir>/relatorio.tex e copia as figuras (PDF) para <dir>/figuras/; compilar com pdflatex três vezes. Os
números do texto são calculados de resultados.json, de modo que o relatório pode ser refeito com mais amostras.
Tabelas e figuras do relatório são numeradas R1, R2, ...; "Tabela 3 do artigo" etc. referem-se ao artigo.
"""
import json
import math
import os
import shutil
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from analise_artigo import ARTIGO_DET, GRUPOS  # noqa: E402

ROTULO = {
    "ref": "referência", "covk0": r"$\mathrm{CoV}(k_v) = 0$", "covk75": r"$\mathrm{CoV}(k_v) = 75\,\%$",
    "covk90": r"$\mathrm{CoV}(k_v) = 90\,\%$", "covk100": r"$\mathrm{CoV}(k_v) = 100\,\%$",
    "covc10": r"$\mathrm{CoV}(c) = 10\,\%$", "covc50": r"$\mathrm{CoV}(c) = 50\,\%$",
    "covc70": r"$\mathrm{CoV}(c) = 70\,\%$", "covphi5": r"$\mathrm{CoV}(\varphi) = 5\,\%$",
    "covphi15": r"$\mathrm{CoV}(\varphi) = 15\,\%$", "covphi20": r"$\mathrm{CoV}(\varphi) = 20\,\%$",
    "s1.5": "$s = 1{,}5$", "s2": "$s = 2$", "s5": "$s = 5$", "s10": "$s = 10$", "s20": "$s = 20$",
    "s400": "$s = 400$", "alfa2": r"$\alpha = 2$", "alfa3": r"$\alpha = 3$", "alfa4": r"$\alpha = 4$",
    "alfa5": r"$\alpha = 5$", "hw0.5": "$h_w/H = 0{,}5$", "hw0.6": "$h_w/H = 0{,}6$", "hw0.7": "$h_w/H = 0{,}7$",
    "hw0.8": "$h_w/H = 0{,}8$", "hw0.9": "$h_w/H = 0{,}9$", "cho_coesivo": "Cho, coesivo",
    "cho_cphi": r"Cho, $c$-$\varphi$",
}
TITULO_GRUPO = {
    "Tabela 3": r"Monte Carlo, casos da Tabela 3 do artigo: coeficientes de variação",
    "Tabela 4": r"Monte Carlo, casos da Tabela 4 do artigo: distâncias de autocorrelação "
                r"$(L_x, L_y) = s\,(20\ \mathrm{m}, 2\ \mathrm{m})$",
    "Tabela 5": r"Monte Carlo, casos da Tabela 5 do artigo: anisotropia $\alpha = k_h/k_v$",
    "Tabela 6": r"Monte Carlo, casos da Tabela 6 do artigo: rebaixamento $h_w/H$",
    "Seção 5.3": r"Monte Carlo, exemplos de Cho (2010) da seção 5.3 do artigo, sem percolação",
}


def f(v, d=3):
    """número com vírgula decimal"""
    if v is None or (isinstance(v, float) and (math.isnan(v) or math.isinf(v))):
        return "---"
    s = f"{abs(v):.{d}f}".replace(".", "{,}")
    return ("$-$" + s) if v < 0 and float(f"{abs(v):.{d}f}") != 0 else s


def p(v, d=1):
    """fração -> porcentagem com vírgula"""
    return "---" if v is None or (isinstance(v, float) and math.isnan(v)) else f(100 * v, d)


def dif(a, b, d=1):
    if a is None or not b:
        return "---"
    x = (a / b - 1) * 100
    return ("$+$" if x >= 0 else "$-$") + f(abs(x), d) + r"\,\%"


def milhar(v):
    """separador de milhar (ponto) a partir de 10 000; números de 4 algarismos sem separador"""
    n = int(round(v))
    return f"{n:,}".replace(",", ".") if abs(n) >= 10000 else str(n)


def dur(h):
    return f"{f(h, 1)}\\,h" if h < 48 else f"{f(h, 0)}\\,h ({f(h / 24, 1)} d)"


def main(dirres, dirtex):
    R = json.load(open(os.path.join(dirres, "resultados.json")))
    mc, det = R["mc"], R["det"]
    conv = det.get("conv", {})
    os.makedirs(os.path.join(dirtex, "figuras"), exist_ok=True)
    for nome in ("fig5_funcional", "fig8_hcrit", "fig9_gamma_beta", "densidades", "tendencias", "convergencia_pf",
                 "mecanismo_A35"):
        src = os.path.join(dirres, nome + ".pdf")
        if os.path.exists(src):
            shutil.copy(src, os.path.join(dirtex, "figuras", nome + ".pdf"))

    # ================================================================== números do texto
    Ns = [st["N"] for st in mc.values()]
    ntot, nmin, nmax = sum(Ns), min(Ns), max(Ns)
    perc = [c for c in mc if c not in ("cho_coesivo", "cho_cphi")]
    mu_e = np.array([mc[c]["escalado"]["mu"] / mc[c]["artigo"]["mu"] - 1 for c in perc])
    mu_b = np.array([mc[c]["mu"] / mc[c]["artigo"]["mu"] - 1 for c in perc])
    sd_e = np.array([mc[c]["escalado"]["sd"] / mc[c]["artigo"]["sd"] - 1 for c in perc])
    com_pf = [c for c in mc if "pf" in mc[c].get("artigo", {})]
    dentro_b = [c for c in com_pf if mc[c]["pf_lo"] <= mc[c]["artigo"]["pf"] <= mc[c]["pf_hi"]]
    entre_pt = [c for c in perc if mc[c]["pf"] <= mc[c]["artigo"]["pf"] <= mc[c]["escalado"]["pf"]]
    fora_pt = [c for c in perc if c not in entre_pt]
    uniao_ic = [c for c in perc if mc[c]["pf_lo"] <= mc[c]["artigo"]["pf"] <= mc[c]["pf_hi"] or
                mc[c]["escalado"]["pf_lo"] <= mc[c]["artigo"]["pf"] <= mc[c]["escalado"]["pf_hi"]]
    ref = mc["ref"]
    rf = {a: conv.get(f"ref_h1_a{a}", {}) for a in range(6)}
    cc = {a: conv.get(f"chocoes_h1_a{a}", {}) for a in range(6)}
    cp = {a: conv.get(f"chocphi_h2_a{a}", {}) for a in range(6)}
    g_conv = rf[5].get("gamma") or rf[4].get("gamma")
    g_sup = [x for x in (rf[4].get("gamma_sup"), rf[5].get("gamma_sup")) if x]
    sp = ref.get("spearman", {})
    modos = R.get("modos", {})
    cmp8 = det.get("fig8_cmp", [])
    bons8 = [r for r in cmp8 if not (r["solo"] == "A" and r["beta"] == 35)]
    d8b = [abs(r["fe"] / r["grad_u_fe"] - 1) for r in bons8]
    a35 = sorted((r for r in cmp8 if r["solo"] == "A" and r["beta"] == 35), key=lambda r: r["hwH"])
    d8r = [r["fe"] / r["grad_u_fe"] - 1 for r in a35 if r["hwH"] >= 0.2 - 1e-9]
    r01 = next((r for r in a35 if abs(r["hwH"] - 0.1) < 1e-9), None)
    r05 = next((r for r in a35 if abs(r["hwH"] - 0.5) < 1e-9), None)
    c8 = det.get("conv8", {})
    a35_l = [c8.get("A_b35_hw2.5_a2", {}).get("gamma"), r05["fe"] / 5 if r05 else None,
             c8.get("A_b35_hw2.5_a4", {}).get("gamma")]
    faixa8 = f"{f(-100 * max(d8r), 0)} a {f(-100 * min(d8r), 0)}"
    hw0 = {(r["solo"], r["beta"]): r for r in cmp8 if r["hwH"] == 0}
    cmp9 = [r for r in det.get("fig9_cmp", []) if r["curva"] == "grad"]
    d9 = {a: [r["fe"] / r["artigo"] - 1 for r in cmp9 if r["alfa"] == a] for a in (1, 5, 10)}
    t56 = det.get("tab56", {})
    dom = {k: det.get("dom", {}).get(k, {}) for k in ("alfa3_dom50", "alfa5_dom25", "alfa5_dom50")}
    f5 = det.get("fig5", {})
    ch, cf = mc["cho_coesivo"], mc["cho_cphi"]
    acf = cf["artigo"]

    # previsão (tempo medido por caso, com os processos simultâneos da máquina do autor)
    def falta(alvo):
        cpu, tot = 0., 0
        for c, st in mc.items():
            S = st["S_artigo"]
            if alvo == "artigo":
                a = S
            elif alvo == "cov5":
                a = S if st["pf"] <= 0 else max(1000, math.ceil((1 - st["pf"]) / (st["pf"] * 0.05 ** 2)))
                a = -(-a // 50) * 50
            else:
                a = min(alvo, S)
            tot += a
            cpu += max(a - st["N"], 0) * st["tempo_medio_s"]
        return tot, cpu / 3600.

    prev = [(alvo, *falta(alvo)) for alvo in (1000, 2000, "cov5", "artigo")]
    t_perc = [mc[c]["tempo_medio_s"] for c in perc]
    t_coes, t_cphi = ch["tempo_medio_s"], cf["tempo_medio_s"]
    cpu_coes = (ch["S_artigo"] - ch["N"]) * t_coes / 3600.
    a5 = mc["alfa5"]
    nf_a5 = round(a5["pf"] * a5["N"])
    sem_falha = [c for c in mc if mc[c]["pf"] <= 0]
    nota_sem_falha = (" Casos ainda sem falhas usam o $S$ do artigo: " + "; ".join(
        f"{ROTULO[c]}, {milhar(mc[c]['S_artigo'])} amostras (cerca de "
        f"{milhar((mc[c]['S_artigo'] - mc[c]['N']) * mc[c]['tempo_medio_s'] / 3600)}\\,h de CPU)"
        for c in sem_falha) + ".") if sem_falha else ""
    rot_fora = [ROTULO[c] for c in fora_pt]
    excl = rot_fora[0] if len(rot_fora) == 1 else ", ".join(rot_fora[:-1]) + " e " + rot_fora[-1] if rot_fora else ""

    L = []
    w = L.append
    w(r"""\documentclass[11pt,a4paper]{article}
\usepackage[T1]{fontenc}
\usepackage[utf8]{inputenc}
\usepackage[brazil]{babel}
\usepackage{lmodern}
\usepackage[margin=2.2cm]{geometry}
\usepackage{amsmath,amssymb}
\usepackage{graphicx}
\usepackage{booktabs}
\usepackage{array}
\usepackage{caption}
\usepackage{float}
\usepackage{xcolor}
\usepackage[hidelinks]{hyperref}
\captionsetup{font=small,labelfont=bf}
\renewcommand{\thetable}{R\arabic{table}}
\renewcommand{\thefigure}{R\arabic{figure}}
\newcommand{\Gam}{\Gamma}
\newcommand{\gradu}{-\mathrm{grad}\,u'_{\mathrm{FE}}}
\newcommand{\vopt}{\underline{\underline{K}}\cdot\underline{v}'_{\mathrm{opt}}}
\newcommand{\Gart}{\Gamma_{\mathrm{artigo}}}
\newcommand{\Gfe}{\Gamma_{\mathrm{FE}}}
\newcommand{\Hcrit}{H_{\mathrm{crit}}}
\definecolor{art}{RGB}{184,57,43}

\title{Talude sob rebaixamento com variabilidade espacial:\\ reprodução de Vargas Ceron et al.\ (2025) com o NeoPZ}
\author{Projects2/SlopeSeepageRandom --- relatório comparativo}
\date{\today}

\begin{document}
\maketitle
""")
    faixaN = milhar(nmin) if nmin == nmax else f"{milhar(nmin)} a {milhar(nmax)}"
    w(r"""\begin{abstract}
Os exemplos de M.~Vargas Ceron, D.~L.~Cecílio, R.~V.~Linn e S.~Maghous, \emph{Stability Analysis of Slope Subjected to
Seepage Forces Considering Spatial Variability of Soil Properties} (IJNAMG 49, 2025, 2459--2491) foram recalculados
com o NeoPZ: elastoplasticidade Mohr-Coulomb incremental, fluxo de Darcy por elementos finitos e campos aleatórios de
Karhunen-Loève. O artigo usa a análise limite cinemática com mecanismos log-espirais. Este relatório compara as
análises determinísticas (Figs.~5, 8 e 9, Tabelas~5 e 6 e exemplos de Cho do artigo) e o Monte Carlo de todos os casos
das Tabelas~3 a 6 e da seção~5.3 do artigo, com """ + milhar(ntot) + r""" amostras (""" + faixaN + r""" por caso), e
estima o tempo para chegar ao número de amostras do artigo. Tabelas e figuras deste relatório são numeradas R1, R2,
\ldots; as do artigo são citadas como ``Tabela~3 do artigo'', ``Fig.~8 do artigo'' etc.
\end{abstract}

\tableofcontents
""")

    # ================================================================== resumo
    w(r"\section{Resumo dos resultados}")
    w(r"\begin{itemize}")
    w(r"\item \textbf{Reprodutibilidade.} Nas amostras em comum entre duas máquinas (o container de desenvolvimento e "
      r"a máquina do autor, 631 amostras) e entre os dois pacotes de resultados recebidos (2800 amostras), os "
      r"resultados são idênticos em todos os dígitos gravados ($\Gam$ e médias dos campos; modo de ruptura nas 431 "
      r"amostras do container que o registram); cada "
      r"amostra depende só de (semente, índice da amostra, campo).")
    w(r"\item \textbf{Exemplos secos (Cho).} No talude $c$-$\varphi$ o FE converge para logo abaixo do limite superior "
      rf"do artigo, como esperado: $\Gam = {f(cp[4].get('gamma'))}$ (nível 4) contra 1,777 e "
      rf"FS $= {f(cp[4].get('fs'))}$ contra 1,203 (artigo) e 1,204 (Cho). No coesivo ainda decresce: "
      rf"$\Gam = {f(cc[3].get('gamma'))}$ e FS $= {f(cc[3].get('fs'))}$ no nível 3 (cerca de 0,6\,\% por nível, "
      rf"{dif(cc[3].get('gamma'), 1.354)} em relação a 1,354; Cho: 1,356).")
    w(rf"\item \textbf{{Caso de referência com percolação.}} O último fator convergido no FE refinado é "
      rf"$\Gam = {f(g_conv)}$ e o colapso ocorre em {f(min(g_sup))}--{f(max(g_sup))}, isto é, "
      rf"{dif(g_conv, 1.336)} a {dif(max(g_sup), 1.336)} acima do $\Gam = 1{{,}}336$ do artigo, que é um limite "
      r"superior. Com $\gamma_w = 9{,}81$ (inferido da Fig.~8 do artigo), a diferença não vem da malha mecânica nem "
      rf"do tamanho do domínio para $\alpha = 1$; $\gamma_w = 10$ reduziria $\Gam$ em "
      rf"{f(-100 * (conv.get('ref_gw10_a3', {}).get('gamma', 1) / rf[3].get('gamma', 1) - 1), 1)}\,\%, o que não "
      r"fecha a diferença. O restante é atribuído, por exclusão, às forças de percolação. Com anisotropia soma-se o "
      r"efeito do domínio truncado (laterais impermeáveis a 10\,m do talude): com 50\,m de crista, pé e base, a "
      rf"diferença de $\alpha = 3$ e $5$ cai para {dif(dom['alfa3_dom50'].get('gamma'), 1.674)} e "
      rf"{dif(dom['alfa5_dom50'].get('gamma'), 1.872)} (nível 3, $\beta = 45^\circ$; $\alpha = 1$ no mesmo nível e "
      rf"domínio: {dif(conv.get('ref_dom_a3', {}).get('gamma'), 1.336)}).")
    w(r"\item \textbf{Fig.~8 do artigo.} Com os solos trocados em relação à Tabela~1 do artigo (seção~\ref{sec:obs}), "
      rf"o FE reproduz as curvas ${{\gradu}}$ com diferença de no máximo {f(100 * max(d8b), 1)}\,\% em três das "
      rf"quatro curvas (nível 3). Na quarta ($\varphi = 32^\circ$, $\beta = 35^\circ$) o FE fica {faixa8}\,\% abaixo "
      rf"para $h_w/H \ge 0{{,}}2$ (nível 3) e continua diminuindo com o refinamento; em $h_w/H = 0{{,}}1$ fica "
      rf"{dif(r01['fe'], r01['grad_u_fe'], 0) if r01 else '---'} e sem rebaixamento não converge com a malha. Como as "
      r"outras três curvas concordam, o limite superior do artigo provavelmente não é justo nesse caso ($\varphi$ "
      r"próximo de $\beta$, mecanismo raso junto à face).")
    w(rf"\item \textbf{{Monte Carlo, média e desvio.}} Nos {len(perc)} casos com percolação (Tabelas~3 a 6 do artigo) "
      rf"a média de $\Gam$ do FE é em média {p(mu_b.mean())}\,\% maior que a do artigo, o mesmo viés do problema médio "
      r"na malha do Monte Carlo. Reescalada amostra a amostra por $\Gart/\Gfe$ do problema médio, coincide com a do "
      rf"artigo (diferença média {p(mu_e.mean())}\,\%, máxima {p(abs(mu_e).max())}\,\%); o desvio padrão reescalado "
      rf"é em média {p(sd_e.mean(), 0)}\,\% maior.")
    w(rf"\item \textbf{{Monte Carlo, probabilidade de falha.}} O Pf do artigo está no intervalo de confiança de 95\,\% "
      rf"do FE em {len(dentro_b)} dos {len(com_pf)} casos. Nos casos com percolação ele fica entre o Pf do FE (menor, "
      rf"pelo viés da média) e o Pf reescalado (maior, pelo desvio maior) em {len(entre_pt)} dos {len(perc)} casos"
      + (rf" (exceções: {excl}, em que o Pf do FE já é maior que o do artigo)" if fora_pt else "") +
      rf", e dentro de pelo menos um dos dois intervalos de confiança em {len(uniao_ic)} dos {len(perc)}. As "
      r"tendências principais das Tabelas~3 a 6 do artigo se repetem (Pf com CoV$(c)$, CoV$(\varphi)$, $s$, $\alpha$ e "
      r"$h_w/H$; $\sigma$ e CoV com $s$; $\mu$ com $\alpha$ e $h_w/H$); variações pequenas ($\mu$ com CoV$(\varphi)$ e "
      r"com $s$, oscilações de $\sigma$ entre valores vizinhos de $s$) ficam dentro do ruído amostral.")
    w(rf"\item \textbf{{Cho.}} Coesivo: $\mu = {f(ch['mu'])}$, $\sigma = {f(ch['sd'])}$ contra "
      rf"{f(ch['artigo']['mu'])} e {f(ch['artigo']['sd'])} da Fig.~17 do artigo; Pf $= {p(ch['pf'])}\,\%$ contra "
      r"6,5\,\% (artigo) e 7,9\,\% (Cho). $c$-$\varphi$: a mediana coincide "
      rf"({f(cf['mediana'])} contra {f(acf.get('mediana'))}); $\sigma = {f(cf['sd'])}$ é dominado por uma amostra "
      rf"extrema ($\Gam = {f(cf['max'], 1)}$; sem ela, $\sigma = {f(cf['sd_sem_max'])}$), contra "
      rf"{f(acf['sd'], 2)}--{f(acf.get('sd_cdf'), 2)} no artigo (Figs.~21 e 22); Pf $= {p(cf['pf'])}\,\%$ contra "
      r"5,5\,\% e 6,37\,\%.")
    if sp.get("mecanismo"):
        m = sp["mecanismo"]
        w(r"\item \textbf{Sensibilidade e modos de ruptura.} Correlações parciais de Spearman de 2ª ordem com médias no "
          rf"mecanismo: $c$ {f(m['c'][0], 2)}, $\varphi$ {f(m['phi'][0], 2)}, $k_v$ {f(m['kv'][0], 2)} (artigo: 0,88, "
          r"0,51, $-$0,04); a correlação com $k_v$ é fraca nos dois, mais negativa no FE.")
    if modos:
        mcs, mcp = modos.get("cho_coesivo", {}), modos.get("cho_cphi", {})
        tc_, tp_ = sum(mcs.values()), sum(mcp.values())
        w(rf"Cho coesivo: {p(mcs.get('abaixo', 0) / tc_)}\,\% das rupturas abaixo do pé (artigo 92,7\,\%); Cho "
          rf"$c$-$\varphi$: {p(mcp.get('pe', 0) / tp_)}\,\% pelo pé e {p(mcp.get('acima', 0) / tp_)}\,\% acima "
          r"(artigo 80,8\,\% e 19,2\,\%).")
    w(rf"\item \textbf{{Previsão.}} Na máquina do autor (8 processos; {f(min(t_perc), 1)}--{f(max(t_perc), 1)}\,s por "
      rf"amostra nos casos com percolação, {f(t_cphi, 1)}\,s no Cho $c$-$\varphi$ e {f(t_coes, 0)}\,s no Cho coesivo) "
      rf"faltam {f(prev[0][2] / 8, 1)}\,h para $N = 1000$ por caso, {f(prev[1][2] / 8, 1)}\,h para $N = 2000$ e "
      rf"{f(prev[3][2] / 8 / 24, 1)} dias para o número de amostras do artigo, {p(cpu_coes / prev[3][2], 0)}\,\% "
      r"disso no Cho coesivo.")
    w(r"\end{itemize}")

    # ================================================================== método
    w(r"""\section{O que foi calculado}
\subsection{Formulação}
O artigo avalia a estabilidade pelo fator $\Gam$ (eq.~54 do artigo), multiplicador das cargas $\gamma'\,\underline{e}_y
- \mathrm{grad}\,u$ no colapso, obtido pela análise limite cinemática com mecanismos log-espirais rotacionais (limite
superior). As forças de percolação vêm de uma solução semianalítica em velocidade de filtração ($\vopt$) ou de
elementos finitos em poropressão (${\gradu}$). O $\Gam = 1{,}336$ do caso de referência (seção~6 e Tabelas~5 e 6 do
artigo) corresponde à curva ${\gradu}$: a Fig.~9 do artigo dá 1,31 em $\beta = 45^\circ$, $\alpha = 1$ para
${\gradu}$, contra 1,89 para $\vopt$ (a diferença de 1,7\,\% entre o valor digitalizado, 1,314, e 1,336 excede um
pouco a precisão de leitura, $\pm 0{,}01$).

O código (\texttt{Projects2/SlopeSeepageRandom} do NeoPZ) calcula o mesmo $\Gam$ como o multiplicador das mesmas
cargas no colapso de uma análise elastoplástica incremental: Mohr-Coulomb associado
(\texttt{TPZPlasticStepVoigt<TPZYCMohrCoulombPV2>}), quadriláteros de ordem 2, colapso definido pela não
convergência do Newton com bissecção do passo até 0,5\,\%. O excesso de poropressão $u$ é a solução de Darcy
estacionário (H1, ordem 2) com as condições de contorno da eq.~21 do artigo e laterais e base impermeáveis. Para a
plasticidade perfeita associada a carga de colapso é única: o FE em deslocamentos converge para ela de cima com o
refinamento, e o mecanismo log-espiral é um limite superior. A comparação direta é, portanto, com as curvas ${\gradu}$
do artigo, e um FE convergido abaixo do valor do artigo indica um limite superior pouco justo.

\subsection{Parâmetros numéricos}
\begin{itemize}
\item $\gamma_w = 9{,}81$\,kN/m$^3$. O artigo não informa $\gamma_w$; os valores da Fig.~8 do artigo sem
rebaixamento só são reproduzidos com 9,81 (com 10 a diferença é de 2,3 a 2,5\,\%).
\item Domínio: crista de 10\,m, pé de 10\,m e base 5\,m abaixo do pé ($25 \times 10$\,m no talude de referência, com
$H = 5$\,m e $\beta = 45^\circ$). Malha estruturada com $h = 1$\,m e refinamento adaptativo guiado pelo mecanismo de
colapso do problema médio. No Monte Carlo, dois níveis: 4914 equações no caso de referência, 18.014 no Cho coesivo
e 3194 no Cho $c$-$\varphi$ ($h = 2$\,m, $H = 10$\,m).
\item Campos aleatórios: KL de Galerkin (quadriláteros de 9 nós, $h_{\mathrm{KL}} = 1$\,m) com todos os modos da KL
discreta ($M$ = número de equações da malha KL: 871, 1081 e 1491); lognormais, covariância exponencial
$\exp(-|\Delta x|/L_x - |\Delta y|/L_y)$, $c$, $\varphi$ e $k_v$ independentes. O erro de discretização da variância
é $\varepsilon_M \approx 3{,}6\,\%$ para $(L_x, L_y) = (20, 2)$\,m (0,007 a 2,3\,\% nos casos da Tabela~4 do artigo),
compensado ponto a ponto. O artigo usa $M = 2000$ termos com erro abaixo de 6\,\%, sem compensação.
\item Monte Carlo direto, $P_f = P(\Gam < 1)$, intervalo de confiança de 95\,\% de Wilson. Colunas com~$^*$:
cada amostra multiplicada por $\Gart/\Gfe$ do problema médio na mesma malha e domínio, o que retira a diferença
determinística e deixa só o efeito da variabilidade.
\item Dados do artigo: tabelas e texto transcritos; Figs.~5, 8, 17, 21, 22 e 26 extraídas dos vetores do PDF e
Fig.~9 digitalizada da imagem, cada uma conferida por um segundo método (diferença máxima de 1,7\,\% nos pontos
lidos). Os momentos $\mu$ e $\sigma$ dos exemplos de Cho foram calculados dessas curvas e têm incerteza maior
(seção~4.1).
\end{itemize}
""")

    # ================================================================== determinísticos
    w(r"\section{Análises determinísticas}")
    w(r"\subsection{Convergência com a malha}\label{sec:conv}")
    w(r"A Tabela~\ref{tab:conv} mostra $\Gam$ e FS (redução de resistência) do problema médio para níveis crescentes "
      r"de refinamento adaptativo.")
    w(r"\begin{table}[H]\centering\small\caption{Convergência com a malha (nível = número de refinamentos "
      r"adaptativos; $h = 1$\,m, ou 2\,m no exemplo $c$-$\varphi$). $\Gam$ é o último fator de carga convergido.}"
      r"\label{tab:conv}")
    w(r"\begin{tabular}{llrrrrrrr}\toprule exemplo & & nível 0 & 1 & 2 & 3 & 4 & 5 & artigo\\\midrule")
    for nome, dd, gk, art in ((r"Cho coesivo", cc, "gamma", "1,354"), ("", cc, "fs", "1,354$^a$"),
                              (r"Cho $c$-$\varphi$ ($H = 10$\,m)", cp, "gamma", "1,777"), ("", cp, "fs", "1,203$^b$"),
                              (r"referência (Tab.~2 do artigo)", rf, "gamma", "1,336"), ("", rf, "fs", "---")):
        rot = r"$\Gam$" if gk == "gamma" else "FS"
        w(f"{nome} & {rot} & " + " & ".join(f(dd[a].get(gk)) if dd[a].get(gk) else "---" for a in range(6))
          + f" & {art}\\\\")
    w(r"\bottomrule\end{tabular}\\[2pt]{\footnotesize $^a$ $\Gam = F_s$ para $\varphi = 0$; Cho (2010), equilíbrio "
      r"limite: 1,356. $^b$ Cho (2010): 1,204.}\end{table}")
    w(r"Nos exemplos secos o FE decresce com o refinamento e fica abaixo do limite superior do artigo (logo abaixo no "
      r"$c$-$\varphi$; no coesivo ainda decresce cerca de 0,6\,\% por nível); para o "
      r"$c$-$\varphi$, a solução log-espiral de Chen (1975), reimplementada aqui, dá $\Gam = 1{,}7770$ e "
      r"$F_s = 1{,}203$ com $H = 10$\,m, os valores do artigo. No caso com percolação o último fator convergido é "
      rf"{f(g_conv)} nos níveis 4 e 5 e com $h = 0{{,}}5$\,m e 3 níveis, e o colapso ocorre em "
      rf"{f(min(g_sup))}--{f(max(g_sup))}; a diferença em relação ao artigo é de {dif(g_conv, 1.336)} a "
      rf"{dif(max(g_sup), 1.336)}.")
    w(rf"Variações do caso de referência (nível 3): $\gamma_w = 10$ dá $\Gam = {f(conv.get('ref_gw10_a3', {}).get('gamma'))}$ "
      rf"(com 9,81: {f(rf[3].get('gamma'))}); o domínio de $105 \times 55$\,m dá "
      rf"{f(conv.get('ref_dom_a3', {}).get('gamma'))}. Nenhuma das duas fecha a diferença em relação ao artigo; o "
      r"restante é atribuído, por exclusão, às forças de percolação (malha e domínio do FE hidráulico do artigo não "
      r"são informados).")
    w(r"Nos níveis mais finos aparecem quedas isoladas fora da sequência: FS = "
      rf"{f(rf[4].get('fs'))} (referência, nível 4) e {f(conv.get('ref_h05_a3', {}).get('fs'))} ($h = 0{{,}}5$\,m, "
      rf"nível 3), e $\Gam = {f(cc[4].get('gamma'))}$ no Cho coesivo, nível 4 (20 cortes de passo; a execução foi "
      r"interrompida antes do FS). As quedas de FS são colapsos prematuros do "
      r"critério de não convergência do Newton em malhas muito finas; o valor do Cho coesivo no nível 4 não foi "
      r"confirmado e não é usado. Com $h = 1$\,m, os níveis 2 e 3 usados nas comparações não apresentam o problema. "
      r"Na seção~\ref{sec:f8} dois valores de nível 4 (9 e 10 cortes de passo) aparecem só como indicação da "
      r"tendência; as conclusões usam o nível 3.")

    w(r"\subsection{Anisotropia e tamanho do domínio}\label{sec:dom}")
    w(r"\begin{table}[H]\centering\small\caption{$\Gam$ (nível 3, $\beta = 45^\circ$) para três domínios "
      r"(crista/pé/base). Laterais e base impermeáveis.}\label{tab:dom}")
    w(r"\begin{tabular}{lrrrr}\toprule $\alpha$ & 10/10/5\,m & 25/25/15\,m & 50/50/50\,m & artigo\\\midrule")
    w(f"1 & {f(rf[3].get('gamma'))} & --- & {f(conv.get('ref_dom_a3', {}).get('gamma'))} & 1,336\\\\")
    w(f"3 & {f(t56.get('alfa3', {}).get('gamma'))} & --- & {f(dom['alfa3_dom50'].get('gamma'))} & 1,674\\\\")
    w(f"5 & {f(t56.get('alfa5', {}).get('gamma'))} & {f(dom['alfa5_dom25'].get('gamma'))} & "
      f"{f(dom['alfa5_dom50'].get('gamma'))} & 1,872\\\\")
    w(r"\bottomrule\end{tabular}\end{table}")
    w(r"Com $k_h > k_v$ as laterais impermeáveis próximas bloqueiam o fluxo horizontal e reduzem as forças de "
      r"percolação; com o domínio grande a diferença em relação ao artigo volta ao patamar de $\alpha = 1$. O efeito "
      r"do domínio foi verificado só em $\beta = 45^\circ$ e $\alpha \le 5$. No Monte Carlo ele é compensado pelas "
      r"colunas reescaladas, que usam o $\Gam$ determinístico no mesmo domínio.")

    w(r"\subsection{Fig.~5 do artigo: funcional hidráulico}")
    w(r"\begin{figure}[H]\centering\includegraphics[width=\textwidth]{figuras/fig5_funcional}"
      r"\caption{$J(u)/(k_h H^2\gamma_w^2)$ do FE em poropressão (malha convergida) para dois domínios (crista/pé/base "
      r"10/10/5\,m e 50/50/50\,m), com as estimativas do artigo $J(u'_{\mathrm{FE}})$ (tracejada) e "
      r"$-J^*(\underline{v}'_{\mathrm{opt}})$ (contínua); só o domínio grande fica entre elas em todos os pontos.}"
      r"\label{fig:f5}"
      r"\end{figure}")
    n_ab = sum(1 for a in (1, 2, 4, 10) for b in (15, 30, 45, 60, 75, 90) if f5.get(f"a{a}_b{b}", {}).get("J"))
    w(r"O funcional depende do tamanho do domínio, que o artigo não informa para o FE. Com crista, pé e base de 50\,m "
      r"(o $L_m = 10H$ do modelo semianalítico), o valor do FE fica entre as duas estimativas do artigo para todos "
      rf"os $\beta$ e $\alpha$, como exige a eq.~28 do artigo; com 10/10/5\,m ele fica abaixo da estimativa inferior "
      rf"em {n_ab - 1} dos {n_ab} pontos (a exceção é $\alpha = 10$, $\beta = 15^\circ$).")

    w(r"\subsection{Fig.~8 do artigo: altura crítica com rebaixamento}\label{sec:f8}")
    w(r"Como o problema é autossemelhante, $\Hcrit = \Gam \cdot H$ com domínio e malha escalados com $H$ "
      r"($\gamma = 18$\,kN/m$^3$, $\alpha = 1$, 20\,m de crista e de pé, 10\,m de base). O solo A tem $c = 6$\,kPa e "
      r"$\varphi = 32^\circ$; o B, $c = 11{,}7$\,kPa e $\varphi = 24{,}7^\circ$.")
    w(r"\begin{figure}[H]\centering\includegraphics[width=\textwidth]{figuras/fig8_hcrit}"
      r"\caption{Fig.~8 do artigo (cinza: contínua $r_p = 0{,}25$ de Wu et al.; tracejada $\vopt$; traço-ponto "
      r"${\gradu}$) e FE (azul: solo que reproduz o artigo em $h_w = 0$; laranja: solo da Tabela~1 do artigo para o "
      r"rótulo do painel), nível 3.}\label{fig:f8}\end{figure}")
    xs = (0., 0.2, 0.4, 0.6, 0.8, 1.0)
    w(r"\begin{table}[H]\centering\small\caption{$\Hcrit$ (m): FE (nível 3, que ainda diminui cerca de 1 a 3\,\% por "
      r"nível de refinamento) contra a curva ${\gradu}$ do artigo. n.c.: sem convergência com a malha.}\label{tab:f8}")
    w(r"\begin{tabular}{llr" + "r" * len(xs) + r"}\toprule painel (solo) & $\beta$ & & \multicolumn{"
      + str(len(xs)) + r"}{c}{$h_w/H$}\\\cmidrule(l){4-" + str(3 + len(xs)) + "} & & & "
      + " & ".join(f(x, 1) for x in xs) + r"\\\midrule")
    for panel, solo, betas in (("London", "B", (30, 60)), ("Israeli", "A", (35, 60))):
        for b in betas:
            pts = {round(r["hwH"], 2): r for r in cmp8 if r["painel"].startswith(panel) and r["beta"] == b}
            for k, (rot, fn) in enumerate((("FE", lambda r: "n.c." if r["nao_conv"] else f(r["fe"], 1)),
                                           ("artigo", lambda r: f(r["grad_u_fe"], 1) + ("$^c$" if r["grad_u_fe"] > 200
                                                                                        else "")),
                                           ("dif.", lambda r: "---" if r["nao_conv"] else dif(r["fe"], r["grad_u_fe"]))
                                           )):
                cab = f"{panel} ({solo}) & ${b}^\\circ$" if k == 0 else " & "
                w(f"{cab} & {rot} & " + " & ".join(fn(pts[x]) if x in pts else "---" for x in xs) + r"\\")
            w(r"\addlinespace")
    w(r"\bottomrule\end{tabular}\\[2pt]{\footnotesize $^c$ fora da escala impressa (até 200\,m), lido do vetor do "
      r"PDF.}\end{table}")
    w(r"Sem rebaixamento, fora o solo A com $\beta = 35^\circ$, o FE reproduz a solução log-espiral (e o artigo) nas "
      r"três curvas (nível 3, que ainda diminui 1,7 a 3,1\,\% por nível): "
      + ", ".join(f"{f(hw0[k]['fe'], 2)}\\,m (solo {k[0]}, ${k[1]}^\\circ$; {dif(hw0[k]['fe'], hw0[k]['grad_u_fe'])})"
                  for k in (("A", 60), ("B", 60), ("B", 30)) if k in hw0) + ".")
    w(r"No solo A com $\beta = 35^\circ$ o FE decresce com o refinamento (em $h_w/H = 0{,}5$: "
      + ", ".join(f"{f(5 * g, 2)}\\,m" for g in a35_l if g) + r" nos níveis 2, 3 e 4, com 7, 8 e 9 cortes de "
      r"passo) e fica abaixo da curva do artigo (" + (f(r05['grad_u_fe'], 2) if r05 else "---") +
      r"\,m) já na malha mais grossa. O mecanismo do FE (Fig.~\ref{fig:mec}) é uma banda curva rasa junto à face, da "
      r"aresta da crista ao pé; sem rebaixamento, ainda mais rasa e quase paralela à face. Como as outras três curvas da "
      r"Fig.~8 do artigo concordam com o FE dentro de 4\,\%, a diferença é específica desse caso ($\varphi$ próximo de "
      r"$\beta$) e indica que o limite superior do artigo não é justo aqui; não foi possível separar o efeito da forma "
      r"do mecanismo do das forças de percolação junto à face, onde são mais intensas. Em "
      rf"$h_w/H = 0{{,}}1$ o FE fica {dif(r01['fe'], r01['grad_u_fe'], 0) if r01 else '---'} e sem rebaixamento esse "
      r"caso não converge com a malha. No solo B com $\beta = 60^\circ$ e $h_w/H = 1$ o nível 4 dá "
      rf"{f(5 * c8.get('B_b60_hw5_a4', {}).get('gamma', float('nan')), 2)}\,m, contra 5,41\,m do artigo.")
    if os.path.exists(os.path.join(dirtex, "figuras", "mecanismo_A35.pdf")):
        w(r"\begin{figure}[H]\centering\includegraphics[width=\textwidth]{figuras/mecanismo_A35}"
          r"\caption{Norma da deformação plástica no colapso, solo A ($c = 6$\,kPa, $\varphi = 32^\circ$), "
          r"$\beta = 35^\circ$, nível 3: sem rebaixamento (esquerda; caso que não converge com a malha) e "
          r"$h_w/H = 0{,}5$ (direita). Modelo com $H = 5$\,m e $\Hcrit = \Gam H$.}\label{fig:mec}\end{figure}")

    w(r"\subsection{Fig.~9 e Tabelas~5 e 6 do artigo: inclinação, anisotropia e rebaixamento}")
    w(r"\begin{figure}[H]\centering\includegraphics[width=\textwidth]{figuras/fig9_gamma_beta}"
      r"\caption{$\Gam \times \beta$ para $\alpha = 1, 5, 10$ ($H = h_w = 5$\,m, $c = 10$\,kPa, $\varphi = 30^\circ$, "
      r"$\gamma = 20$\,kN/m$^3$). Cinza: Fig.~9 do artigo (contínua ${\gradu}$, tracejada $\vopt$). Azul: FE, nível 3, "
      r"crista 10\,m, pé 10\,m, base 5\,m.}\label{fig:f9}\end{figure}")
    bs = (30, 45, 60, 75, 90)
    fig9 = det.get("fig9", {})
    w(r"\begin{table}[H]\centering\small\caption{Fig.~9 do artigo: FE (nível 3) contra a curva ${\gradu}$ do artigo "
      r"(digitalizada da imagem, $\pm 0{,}01$; dif.\ calculada com os valores digitalizados sem arredondamento).}\label{tab:f9}\begin{tabular}{llrrrrr}\toprule $\alpha$ & & " +
      " & ".join(f"$\\beta = {b}^\\circ$" for b in bs) + r"\\\midrule")
    for a9 in (1, 5, 10):
        pts = {r["beta"]: r for r in cmp9 if r["alfa"] == a9}

        def fe9(b, a9=a9):
            g = fig9.get(f"beta{b}_alfa{a9}", {}).get("gamma")
            return f(g) if g else "---"
        linhas9 = (("FE", fe9),
                   ("artigo", lambda b, pts=pts: f(pts[b]["artigo"], 2) if b in pts else "$> 5$"),
                   ("dif.", lambda b, pts=pts: dif(pts[b]["fe"], pts[b]["artigo"]) if b in pts else "---"))
        for k, (rot, fn) in enumerate(linhas9):
            w(f"{a9 if k == 0 else ''} & {rot} & " + " & ".join(fn(b) for b in bs) + r"\\")
        w(r"\addlinespace")
    w(r"\bottomrule\end{tabular}\end{table}")
    w(rf"Para $\alpha = 1$ o FE fica {f(100 * min(d9[1]), 0)} a {f(100 * max(d9[1]), 0)}\,\% acima da curva "
      rf"${{\gradu}}$, o mesmo resíduo do caso de referência; para $\alpha = 5$ e 10, {f(100 * min(d9[5] + d9[10]), 0)} "
      rf"a {f(100 * max(d9[5] + d9[10]), 0)}\,\% acima. Na Tabela~\ref{{tab:dom}} ($\beta = 45^\circ$, $\alpha = 3$ e 5) "
      r"o domínio grande reduz a diferença ao patamar de $\alpha = 1$; para $\alpha = 10$ e outros $\beta$ o efeito do "
      r"domínio não foi verificado.")
    w(r"\begin{table}[H]\centering\small\caption{Coluna determinística das Tabelas~5 e 6 do artigo "
      r"(domínio 10/10/5\,m).}\label{tab:t56}"
      r"\begin{tabular}{lrrrr}\toprule caso & FE nível 2 (malha do MC) & FE nível 3 & artigo & dif. (nível 3)"
      r"\\\midrule")
    linhas = [(r"$\alpha = 1$", rf[2], rf[3], 1.336)]
    for a in (2, 3, 4, 5):
        linhas.append((rf"$\alpha = {a}$", t56.get(f"alfa{a}_a2", {}), t56.get(f"alfa{a}", {}), ARTIGO_DET[f"alfa{a}"]))
    for h, hw in ((0.5, "2.5"), (0.6, "3"), (0.7, "3.5"), (0.8, "4"), (0.9, "4.5")):
        linhas.append((f"$h_w/H = {f(h, 1)}$", t56.get(f"hw{hw}_a2", {}), t56.get(f"hw{hw}", {}),
                       ARTIGO_DET[f"hw{h}"]))
    for nome, d2, d3, a in linhas:
        w(f"{nome} & {f(d2.get('gamma'))} & {f(d3.get('gamma'))} & {f(a)} & {dif(d3.get('gamma'), a)}\\\\")
    w(r"\bottomrule\end{tabular}\end{table}")

    # ================================================================== Monte Carlo
    w(r"\section{Monte Carlo}")
    w(r"Nas tabelas, $\mu$, $\sigma$ e CoV são de $\Gam$; Pf com o intervalo de confiança de 95\,\%; as colunas "
      r"com~$^*$ são as reescaladas pelo $\Gam$ determinístico (seção~2.2); as três últimas colunas, em vermelho, são "
      r"os valores do artigo.")
    for titulo, casos in GRUPOS:
        chave = next(k for k in TITULO_GRUPO if titulo.startswith(k))
        rotulo_tab = {"Tabela 3": "tab:mc3", "Tabela 4": "tab:mc4", "Tabela 5": "tab:mc5", "Tabela 6": "tab:mc6",
                      "Seção 5.3": "tab:mccho"}[chave]
        w(r"\begin{table}[H]\centering\footnotesize\setlength{\tabcolsep}{3.2pt}\caption{" + TITULO_GRUPO[chave] +
          r"}\label{" + rotulo_tab + "}")
        w(r"\begin{tabular}{lrrrrlrrr>{\color{art}}r>{\color{art}}r>{\color{art}}r}\toprule & "
          r"\multicolumn{5}{c}{FE} & \multicolumn{3}{c}{FE reescalado} & \multicolumn{3}{c}{\color{art}artigo}\\"
          r"\cmidrule(lr){2-6}\cmidrule(lr){7-9}\cmidrule(l){10-12}"
          r"caso & $N$ & $\mu$ & $\sigma$ & CoV\,\% & Pf\,\% (IC 95\,\%) & $\mu^*$ & $\sigma^*$ & Pf$^*$\,\% & "
          r"$\mu$ & $\sigma$ & Pf\,\%\\\midrule")
        for c in casos:
            st = mc.get(c)
            if not st:
                continue
            a = st.get("artigo", {})
            e = st.get("escalado") or {}
            pfa = p(a.get("pf"), 2) if "pf" in a else "---"
            if "fonte_momentos" in a:
                mua = f(a["mu"]) + "$^\\dagger$"
                sda = f(a["sd"]) if c == "cho_coesivo" else f"{f(a['sd'], 2)}--{f(a.get('sd_cdf'), 2)}"
            else:
                mua, sda = f(a.get("mu")), f(a.get("sd"))
            w(f"{ROTULO.get(c, c)} & {milhar(st['N'])} & {f(st['mu'])} & {f(st['sd'])} & {p(st['cov'])} & "
              f"{p(st['pf'], 2)} ({p(st['pf_lo'])}--{p(st['pf_hi'])}) & {f(e.get('mu'))} & {f(e.get('sd'))} & "
              f"{p(e.get('pf'), 2) if e else '---'} & {mua} & {sda} & {pfa}\\\\")
        w(r"\bottomrule\end{tabular}")
        if chave == "Seção 5.3":
            w(r"\\[2pt]\parbox{\linewidth}{\footnotesize $^\dagger$ $\mu$ e $\sigma$ do artigo calculados das curvas das Figs.~17, 21 e 22 "
              r"do artigo (não tabelados); no $c$-$\varphi$ a densidade da Fig.~21 é cortada na cauda e a CDF da Fig.~22 "
              rf"dá $\mu = {f(acf.get('mu_cdf'), 2)}$ e $\sigma = {f(acf.get('sd_cdf'), 2)}$. Pf de Cho (2010): 7,9\,\% "
              r"(coesivo) e 6,37\,\% ($c$-$\varphi$).}")
        w(r"\end{table}")
    w(r"\subsection{Exemplos de Cho}")
    w(rf"No talude coesivo a distribuição de $\Gam$ do FE é quase a do artigo ($\mu = {f(ch['mu'])}$ e "
      rf"$\sigma = {f(ch['sd'])}$ contra {f(ch['artigo']['mu'])} e {f(ch['artigo']['sd'])}); a média um pouco menor "
      rf"acompanha o $\Gam$ determinístico do FE na malha do Monte Carlo ({f(ch['gamma_det_fe'])} contra 1,354) e "
      rf"leva a um Pf maior ({p(ch['pf'])}\,\%, ou {p(ch['escalado']['pf'])}\,\% reescalado, contra 6,5\,\% no artigo "
      r"e 7,9\,\% em Cho). No talude $c$-$\varphi$ os quantis do FE e do artigo (CDF da Fig.~22) são próximos: "
      + ", ".join(f"{int(q[1:])}\\,\\%: {f(cf[q], 2)} e {f(acf['quantis'][q], 2)}" for q in ("q05", "q25", "q75", "q95"))
      + rf"; mediana {f(cf['mediana'], 2)} e {f(acf['mediana'], 2)}. O desvio padrão maior do FE "
      rf"({f(cf['sd'])}) vem da cauda superior: uma amostra tem $\Gam = {f(cf['max'], 1)}$ e, sem ela, "
      rf"$\sigma = {f(cf['sd_sem_max'])}$.")
    w(r"\begin{figure}[H]\centering\includegraphics[width=\textwidth]{figuras/densidades}"
      r"\caption{Densidade de $\Gam$: referência (lognormal do artigo, $\mu = 1{,}353$, $\sigma = 0{,}318$) e exemplos "
      r"de Cho (curvas do artigo extraídas das Figs.~17 e 21, com a densidade de $F_s$ de Cho).}\label{fig:dens}"
      r"\end{figure}")
    w(r"\subsection{Casos com percolação (Tabelas~3 a 6 do artigo)}")
    w(r"\begin{figure}[H]\centering\includegraphics[width=\textwidth]{figuras/tendencias}"
      r"\caption{Tendências das Tabelas~3 a 6 do artigo (Tabelas~\ref{tab:mc3} a \ref{tab:mc6}): Pf, $\mu(\Gam)$ e "
      r"CoV$(\Gam)$. Vermelho: artigo; azul: FE (barras: IC de 95\,\% de Pf); verde: FE reescalado. Eixo de $s$ em "
      r"escala logarítmica.}\label{fig:tend}\end{figure}")
    w(rf"Nos {len(perc)} casos com percolação: (i) a média do FE é {p(mu_b.mean())}\,\% maior que a do artigo, o "
      r"mesmo viés do problema médio na malha do Monte Carlo (Tabela~\ref{tab:t56}); reescalada, coincide com a do "
      rf"artigo (diferença média {p(mu_e.mean())}\,\%); (ii) o desvio padrão reescalado é {p(sd_e.mean(), 0)}\,\% "
      r"maior em média. Hipóteses para (ii), não testadas aqui: o FE forma mecanismos que seguem as zonas fracas, "
      r"enquanto a família de superfícies log-espirais média a resistência ao longo de curvas suaves (efeito "
      r"conhecido do RFEM); o artigo trunca a KL sem compensar a variância perdida; e usa em cada arco o "
      r"$\varphi$ máximo do subdomínio, o que suaviza os valores baixos. (iii) Com isso o Pf do artigo fica entre o Pf "
      rf"do FE e o reescalado em {len(entre_pt)} dos {len(perc)} casos"
      + (rf" (exceções: {excl}, com Pf do FE já maior que o do artigo)" if fora_pt else "") +
      r". As tendências principais se repetem: Pf cresce muito com CoV$(c)$, pouco com CoV$(\varphi)$ e quase nada "
      r"com CoV$(k_v)$; cresce com a escala $s$ das distâncias de autocorrelação; cai com $\alpha$; e o mínimo de "
      r"estabilidade ocorre perto de $h_w/H = 0{,}8$--$0{,}9$; $\sigma$ e CoV crescem com $s$. Variações pequenas "
      r"($\mu$ com CoV$(\varphi)$ e com $s$, oscilações de $\sigma$ entre valores vizinhos de $s$) ficam dentro do "
      r"ruído amostral.")
    w(rf"Nos casos com Pf pequeno ainda há poucas falhas (por exemplo, $\alpha = 5$: {nf_a5} em {a5['N']}), e o Pf "
      r"deles só fica bem determinado com muito mais amostras (seção~\ref{sec:prev}).")
    w(r"\begin{figure}[H]\centering\includegraphics[width=0.62\textwidth]{figuras/convergencia_pf}"
      r"\caption{Convergência de Pf com o número de amostras (tracejadas: valores do artigo).}\label{fig:conv}"
      r"\end{figure}")
    if sp.get("mecanismo"):
        m, d = sp["mecanismo"], sp["dominio"]
        w(r"\subsection{Sensibilidade (seção 6.1 e Fig.~26 do artigo)}")
        w(r"\begin{table}[H]\centering\small\caption{Correlação parcial de Spearman de 2ª ordem entre $\Gam$ e cada "
          rf"variável (caso de referência, $N = {sp['N']}$).}}\begin{{tabular}}{{lrrr}}\toprule variável & FE, médias "
          r"no mecanismo & FE, médias no domínio & artigo\\\midrule")
        for k, lab, art in (("c", "$c$", "0,879"), ("phi", r"$\varphi$", "0,508"), ("kv", "$k_v$", "$-$0,036")):
            w(f"{lab} & {f(m[k][0])} & {f(d[k][0])} & {art}\\\\")
        w(r"\bottomrule\end{tabular}\end{table}")
        w(r"Médias no mecanismo: $c$ e $\varphi$ médios na banda de cisalhamento ($\|\Delta\varepsilon^p\| \ge 10\,\%$ "
          r"do máximo no colapso, ponderados por $\|\Delta\varepsilon^p\|$) e $k_v$ médio na massa que se move "
          r"($\|\Delta u\| \ge 30\,\%$ do máximo), análogas às médias ao longo da superfície de ruptura e no volume "
          r"deslizante usadas no artigo (Hu et al.).")
    if modos:
        w(r"\subsection{Modos de ruptura (Figs.~16 e 20 do artigo)}")
        w(r"\begin{table}[H]\centering\small\caption{Modos de ruptura: banda $\|\Delta\varepsilon^p\| \ge 20\,\%$ do "
          r"máximo no colapso; abaixo do pé se ela desce mais de $0{,}1H$ abaixo do nível do pé, pelo pé se passa a "
          r"menos de $0{,}1H$ do pé.}\begin{tabular}{lrrrrl}\toprule exemplo & $N$ & acima & pelo pé & abaixo & "
          r"artigo (acima/pelo/abaixo)\\\midrule")
        art = {"cho_coesivo": r"0,55 / 6,76 / 92,68\,\%", "cho_cphi": r"19,20 / 80,80 / 0\,\%", "ref": "---"}
        for c in ("cho_coesivo", "cho_cphi", "ref"):
            mm = modos.get(c)
            if not mm:
                continue
            t = sum(mm.values())
            w(f"{ROTULO[c]} & {t} & " + " & ".join(f"{p(mm.get(k, 0) / t)}\\,\\%" for k in ("acima", "pe", "abaixo"))
              + f" & {art[c]}\\\\")
        w(r"\bottomrule\end{tabular}\end{table}")

    # ================================================================== observações
    hw0_txt = ", ".join(f"{f(hw0[k]['fe'], 2)}\\,m" for k in (("A", 60), ("B", 60), ("B", 30)) if k in hw0)
    w(r"""\section{Observações sobre o artigo}\label{sec:obs}
\begin{enumerate}
\item \textbf{Fig.~8 com os solos trocados.} Sem rebaixamento ($h_w/H = 0$), a solução log-espiral de Chen (1975)
com $\gamma' = 18 - 9{,}81$ dá, para $c = 6$\,kPa e $\varphi = 32^\circ$ (``London clay'' na Tabela~1 do artigo),
$\Hcrit = 13{,}00$\,m ($\beta = 60^\circ$) e $229{,}1$\,m ($\beta = 35^\circ$), que são as curvas do painel ``Israeli
Clay'' (12,98\,m e 229,3\,m, este fora da escala impressa, lido do vetor do PDF); para $c = 11{,}7$\,kPa e
$\varphi = 24{,}7^\circ$ dá 17,97\,m ($\beta = 60^\circ$) e 156,7\,m ($\beta = 30^\circ$), as do painel ``London
Clay'' (17,95 e 156,6\,m). Para $\varphi = 32^\circ$, $\beta = 35^\circ$ o FE não converge com a malha
(seção~\ref{sec:f8}); nas outras três curvas ele dá os mesmos valores (""" + hw0_txt + r""", diferença de 0,7 a
2,1\,\%, nível 3).
\item \textbf{$\gamma_w$ não informado.} Os valores da Fig.~8 do artigo sem rebaixamento só são reproduzidos com
$\gamma_w = 9{,}81$\,kN/m$^3$.
\item \textbf{Altura do exemplo $c$-$\varphi$ de Cho.} O texto da seção~5.3.2 e o rótulo da cota da Fig.~19 do
artigo dão $H = 5$\,m (o desenho, em escala, tem $H = 10$\,m e $30 \times 15$\,m), mas $\Gam = 1{,}777$ e
$F_s = 1{,}203/1{,}204$ correspondem a $H = 10$\,m: Chen dá $\Gam = 1{,}7770$ e $F_s = 1{,}203$ com $H = 10$\,m, e
3,554 e 1,60 com $H = 5$\,m.
\item \textbf{Seção 6.1.} O texto fala em correlação ``perfeita'' de $\Gam$ com a coesão, mas a barra da Fig.~26 do
artigo vale 0,88.
\item \textbf{Forças de percolação do caso de referência.} Mesmo usando a mesma aproximação (${\gradu}$), o FE
convergido fica acima do $\Gam = 1{,}336$ do artigo (seção~\ref{sec:conv}), o que aponta para uma diferença na solução
hidráulica (malha e domínio do FE hidráulico do artigo não são informados).
\item \textbf{Limite superior com $\varphi$ próximo de $\beta$.} No solo A com $\beta = 35^\circ$ o FE fica
""" + faixa8 + r"""\,\% abaixo do limite superior do artigo para $h_w/H \ge 0{,}2$ e continua diminuindo com o
refinamento, com um mecanismo raso junto à face (seção~\ref{sec:f8}).
\end{enumerate}
""")

    # ================================================================== previsão
    w(r"\section{Previsão de tempo}\label{sec:prev}")
    w(rf"Tempo por amostra medido na máquina do autor com 8 processos simultâneos: {f(min(t_perc), 1)} a "
      rf"{f(max(t_perc), 1)}\,s nos casos com percolação, {f(t_cphi, 1)}\,s no Cho $c$-$\varphi$ e {f(t_coes, 1)}\,s "
      rf"no Cho coesivo (malha de 18 mil equações). A Tabela~\ref{{tab:prev}} dá o que falta a partir das "
      rf"{milhar(ntot)} amostras já calculadas.")
    w(r"\begin{table}[H]\centering\small\caption{Tempo que falta, a partir das amostras atuais. A coluna de 16 "
      r"processos supõe 16 núcleos físicos e o mesmo tempo por amostra.}\label{tab:prev}"
      r"\begin{tabular}{lrrrr}\toprule meta por caso & amostras no total & CPU (h) & 8 processos & 16 processos"
      r"\\\midrule")
    nomes = {1000: "$N = 1000$", 2000: "$N = 2000$", "cov5": r"CoV(Pf) $< 5\,\%$ com o Pf do FE$^d$",
             "artigo": r"$S$ do artigo (10 mil a 100 mil)"}
    for alvo, tot, cpu in prev:
        w(f"{nomes[alvo]} & {milhar(tot)} & {milhar(cpu)} & {dur(cpu / 8)} & {dur(cpu / 16)}\\\\")
    w(r"\bottomrule\end{tabular}\\[2pt]\parbox{\linewidth}{\footnotesize $^d$ fórmula aplicada ao Pf atual de cada caso (no mínimo 1000, "
      r"em blocos de 50); muito incerta onde há poucas falhas: $\alpha = 5$ tem " + str(nf_a5) + " falha"
      + ("s" if nf_a5 != 1 else "") + " em " + str(a5["N"]) + r" amostras, e o alvo dele vai de cerca de 30 mil a 1 "
      r"milhão de amostras no intervalo de confiança de Pf." + nota_sem_falha + r"}\end{table}")
    w(rf"No alvo do artigo o Cho coesivo responde por {milhar(cpu_coes)} das {milhar(prev[3][2])}\,h de CPU que "
      r"faltam. Com $h = 2$\,m nesse caso (cerca de 7 mil equações) o custo por amostra cai para cerca de um quinto "
      r"(medição no container de desenvolvimento, 4 amostras simultâneas: 6,4\,s com $h = 2$\,m contra 29,7\,s "
      r"com $h = 1$\,m), com $\Gam$ cerca de 1\,\% maior. Na ferramenta, "
      r"\texttt{-{}-alvo cov5} leva primeiro cada caso a 1000 amostras (\texttt{-{}-min}) e reavalia o alvo cada vez "
      r"que a fila é gerada de novo.")

    # ================================================================== reprodução
    w(r"""\section{Reprodução}
{\sloppy Branch \texttt{claude/great-clarke-xist30} do neopz-master; dados, logs, figuras e este relatório em
\texttt{Projects2/SlopeSeepageRandom/resultados/artigo2025/} (\texttt{campanha/}: amostras juntadas;
\texttt{det/}: logs determinísticos; \texttt{ref/}: dados extraídos do artigo; \texttt{relatorio\_tex/}: este
documento).\par}
\par\noindent\begin{minipage}{\textwidth}\footnotesize
\begin{verbatim}
sudo apt install cmake g++ make liblapack-dev libblas-dev liblapacke-dev \
     python3-numpy python3-scipy python3-matplotlib texlive-latex-extra texlive-lang-portuguese
# Monte Carlo: compila, gera a fila (blocos de 50 amostras, casos intercalados) e roda
bash Projects2/SlopeSeepageRandom/scripts/rodar_campanha.sh 1000   # [N|cov5|artigo] [proc.] [dir]
S=$PWD/Projects2/SlopeSeepageRandom/scripts
A=$PWD/Projects2/SlopeSeepageRandom/resultados/artigo2025
PROCESSOS=8 python3 $S/campanha_artigo.py status ~/campanha_artigo
# análise (campanha: ~/campanha_artigo ou $A/campanha) e relatórios
python3 $S/analise_artigo.py ~/campanha_artigo $A/det $A/ref saida
python3 $S/relatorio_tex.py saida relatorio_tex
cd relatorio_tex && for i in 1 2 3; do pdflatex relatorio.tex; done
\end{verbatim}
\end{minipage}
\end{document}
""")
    os.makedirs(dirtex, exist_ok=True)
    with open(os.path.join(dirtex, "relatorio.tex"), "w") as fo:
        fo.write("\n".join(L))
    print(os.path.join(dirtex, "relatorio.tex"))


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(1)
    main(*sys.argv[1:])
