#!/usr/bin/env python3
"""gera_relatorio.py: o relatorio parcial consolidado (plano de trabalho, formulacoes, testes, hipoteses e o mapa dos
relatorios), com os numeros lidos das campanhas baixadas do fermi e das rodadas do PC. Para atualizar:
    cd docs/confiabilidade_campos && python3 -I gera_relatorio.py && pdflatex relatorio_parcial.tex && pdflatex relatorio_parcial.tex
"""
import glob
import json
import math
import os
import runpy
import shutil

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

AQUI = os.path.dirname(os.path.abspath(__file__))
NEOPZ = os.path.abspath(f"{AQUI}/../..")
ETAPAS = f"{NEOPZ}/Projects/PlanoDeTrabalho/ETAPAS1a3"
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
for c, rot in (("ai_banda_tri0", "AI pela banda, piloto na Tri(0)"), ("ai_form100_tri0", "AI pelo FORM $M = 100$, direção da Tri(0)"),
               ("ai_form100_modal", "AI pelo FORM $M = 100$, gradiente modal"),
               ("ai_form200_modal", "AI pelo FORM $M = 200$, gradiente modal")):
    if glob.glob(f"{G}/{c}/mc_*_*.csv"):
        pc.append((rot, ai(f"{G}/{c}")))
    else:
        pc.append((rot, None))


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


def linha_ai(rot, a, d):
    if a is None:
        return f"{rot} & --- & na fila & --- & --- & --- \\\\"
    return (f"{rot} & {a['n']} & {pct(a['p'])} & {num(100 * a['cov'], 1)} & {num(desvios(a, d), 1)} & "
            f"{num(ganho(a, d), 1)} \\\\")


def tabela(cab, linhas, fmt):
    return ("\\resizebox{\\textwidth}{!}{%\n\\begin{tabular}{" + fmt + "}\n\\toprule\n" + cab + " \\\\\n\\midrule\n" +
            "\n".join(linhas) + "\n\\bottomrule\n\\end{tabular}}")


T_R05 = tabela("Método & $N$ & $P_f$ (\\%) & CoV (\\%) & dif. (desvios) & ganho",
               [f"Monte Carlo direto (fermi) & {d05['n']} & {pct(d05['p'])} & {num(100 * d05['cov'], 1)} & --- & 1 \\\\"] +
               [linha_ai(r, a, d05) for r, a in a05] + [linha_ai(r, a, d05) for r, a in pc], "lrrrrr")
T_R07 = tabela("Método & $N$ & $P_f$ (\\%) & CoV (\\%) & dif. (desvios) & ganho",
               [f"Monte Carlo direto (fermi) & {d07['n']} & {pct(d07['p'])} & {num(100 * d07['cov'], 1)} & --- & 1 \\\\"] +
               [linha_ai(r, a, d07) for r, a in a07], "lrrrrr")
T_CH = tabela("Método & $N$ & $P_f$ (\\%) & CoV (\\%) & dif. (desvios) & ganho",
              [f"Monte Carlo direto (fermi) & {dch['n']} & {pct(dch['p'])} & {num(100 * dch['cov'], 1)} & --- & 1 \\\\"] +
              [linha_ai(r, a, dch) for r, a in ach], "lrrrrr")
T_FORM = tabela("$M$ & variância retida (Nyström) & $\\beta$ analítica & $\\beta$ Nyström & análises por FORM",
                [f"{M} & {num(VAR_NY[M], 3)} & " + (num(FM[M][0]['b'], 3) if FM[M][0] else "---") + " & " +
                 (num(FM[M][1]['b'], 3) if FM[M][1] else "---") + " & " +
                 f"{(FM[M][1] or FM[M][0] or {'n': '---'})['n']} \\\\" for M in Ms], "rrrrr")
L = E["LINTAB"]
T_LIN = tabela("$\\beta$ & $P_f$ exata & direto ($10^6$) & AI na direção exata & AI a " +
               num(math.degrees(math.acos(L[1])), 1) + "$^\\circ$ da direção",
               [f"{x['b']} & {x['ex']:.3e} & " + (f"{pct(x['d'])}\\,\\%" if x['d'] else "---") +
                f" & {num(100 * (x['ai20'][0] / x['ex'] - 1), 1)}\\,\\% ({num(x['ai20'][2], 1)}) & "
                f"{num(100 * (x['ai10'][0] / x['ex'] - 1), 1)}\\,\\% ({num(x['ai10'][2], 1)})" + " \\\\" for x in L[0]],
               "rrrrr") if L else "(teste linear não baixado)"
C05, C07 = E["CUSTO05"], E["CUSTO07"]
T_CUSTO = tabela("Método & $r = -0{,}5$ & $r = -0{,}7$",
                 [f"{a[0]} & {a[3]:.0f} & {b[3]:.0f} \\\\" for a, b in zip(C05, C07)], "lrr")
cho, chuva, mc = E["cho"], E["chuva"], E["mcnemar"]

MAPA = r"""
\begin{longtable}{p{0.40\textwidth}p{0.55\textwidth}}
\toprule
Documento & O que tem \\
\midrule
\endhead
\multicolumn{2}{l}{\textbf{Formulação e modelo}} \\
\path{docs/biot_nao_saturado/biot_nao_saturado_agua_ar5.pdf} & formulação poroelástica não saturada (esqueleto, água e ar), a implementação no NeoPZ, Rahardjo et al. (2007), benchmarks (49 p.) \\
\path{~/Dropbox/UnsaturatedSoil/formulacao/formulacao/formulacao.pdf} & formulação acoplada $u$--$p_w$--$p_a$ com Mohr--Coulomb associado: discretização, resíduos e jacobiano (25 p.) \\
\path{~/Dropbox/UnsaturatedSoil/formulacao/formulacao/resumo_nao_saturado.pdf} & resumo da formulação $u$--$p_w$--$p_a$ (13 p.) \\
\path{docs/slope_mohr_coulomb/relatorio.pdf} & talude Mohr--Coulomb no NeoPZ: correções da plasticidade, FS e rebaixamento (18 p.) \\
\path{docs/slope_mohr_coulomb/nota_uqp.pdf} & nota: formulação de três campos $u$--$q$--$p$ (6 p.) \\
\path{docs/fs_arc_length/relatorio_fs.pdf} & FS de talude não saturado sob chuva: SRM e GIM com e sem arc-length, sensibilidade e parâmetros (34 p.) \\
\multicolumn{2}{l}{\textbf{Campos aleatórios, FORM e amostragem por importância}} \\
\path{~/relatorio_kl_form/relatorio.pdf} & estudo do Ghanem \& Spanos: KL analítica, Nyström e Galerkin, MC e FORM nas placas (quadrada e curva) \\
\path{Projects/SlopeRandomFieldRefinement/relatorio/implementacao_form_ai.pdf} & como o FORM e a AI foram implementados no NeoPZ (5 p.) \\
\path{Projects/Cho2010/nystrom_ai/ai_nystrom.pdf} & AI com a KL de Nyström, um roteiro testado em Mathematica (5 p.) \\
\path{Projects/Confiabilidade/doc/confiabilidade.pdf} & FORM, Monte Carlo e AI em variáveis; talude infinito \\
\path{Projects/GuiaEstudo/guia_estudo.pdf} & guia de estudo que liga tudo (inclui o FORM com a base da KL, Eq.~14) \\
\multicolumn{2}{l}{\textbf{Resultados das campanhas}} \\
\path{Projects/SlopeRandomFieldRefinement/relatorio/relatorio_chuva.pdf} & Monte Carlo da chuva (Rahardjo), Seção 7: AI \\
\path{Projects/SlopeRandomFieldRefinement/relatorio/relatorio_rebaixamento.pdf} & rebaixamento do Ceron, Monte Carlo direto e AI (validação 15,70\,\% contra 15,87\,\%) \\
\path{Projects/SlopeRandomFieldRefinement/relatorio/relatorio_solos.pdf} & os três solos do Rahardjo sob chuva \\
\path{Projects/PlanoDeTrabalho/README.md} & o plano de trabalho (etapas 0 a 5, 2b e 2c) \\
\path{Projects/PlanoDeTrabalho/ETAPA0/relatorio_etapa0.pdf} & etapa 0: a AI da chuva no fermi e a prévia do Cho2010 (SRM antigo; superado pelo das etapas 1 a 3) \\
\path{Projects/PlanoDeTrabalho/ETAPAS1a3/relatorio_etapas.pdf} & etapas 1 a 3: SRM corrigido, FORM e AI, $P_f$ pequena, teste linear, piloto, custo \\
\path{docs/confiabilidade_campos/relatorio_parcial.pdf} & \textbf{este}: o consolidado \\
\bottomrule
\end{longtable}
"""

tex = r"""\documentclass[11pt,a4paper]{article}
\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}
\usepackage[brazil]{babel}
\usepackage{lmodern}
\usepackage[margin=2cm]{geometry}
\usepackage{amsmath,amssymb,bm}
\usepackage{graphicx,booktabs,float,longtable,url}
\usepackage[hidelinks]{hyperref}
\title{Confiabilidade de taludes com campos aleatórios no NeoPZ\\\large Relatório parcial: plano de trabalho, formulações,
testes, hipóteses e o mapa dos relatórios}
\date{16/11/2026}
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
processos); o PC serve para testes curtos. Este relatório reúne o que está espalhado em vários documentos (Seção~\ref{sec:mapa}).

\section{Formulações}
\subsection{O modelo mecânico e o fator de segurança}
\begin{itemize}
  \item \textbf{Material:} Mohr--Coulomb elastoplástico com fluxo associado, dentro do material não saturado
    \texttt{TPZMatUnsaturatedUPwPa} ($u$--$p_w$--$p_a$; \path{docs/biot_nao_saturado}, \path{formulacao.pdf}). No talude
    seco, as pressões ficam congeladas em zero.
  \item \textbf{FS pelo SRM:} $c'/F$ e $\tan\varphi'/F$ até o colapso, com arc-length no parâmetro $F$ e refinamento
    adaptativo em $h$ (tri2 + 4 níveis; \path{docs/fs_arc_length}). Esse FS é o do equilíbrio limite (Ceron 2025, Eq.~58).
  \item \textbf{Correção do SRM (etapa 1):} o peso é aplicado com $F_0 = 0{,}8$, e não com 1. Com 1, o peso às vezes
    convergia num talude com FS logo abaixo de 1, e o FS saía $\ge 1$, o que deixava um buraco no histograma.
\end{itemize}

\subsection{Campos aleatórios: Karhunen--Loève}
O campo gaussiano é $Z(\bm x) = \sum_{k=1}^M \sqrt{\lambda_k}\, f_k(\bm x)\, \xi_k$, com $\xi_k \sim N(0,1)$
independentes, e $c' = \exp(\mu + \sigma Z_c)$, $\varphi' = \exp(\mu + \sigma(\rho Z_c + \sqrt{1-\rho^2} Z_\varphi))$.
Sortear um campo é sortear o vetor $\bm\xi$. Há duas bases, as duas com o kernel exponencial separável ($l_h = 20$\,m,
$l_v = 2$\,m no Cho):
\begin{itemize}
  \item \textbf{Analítica do retângulo} que contém o talude (\texttt{RandomField}, 2000 termos): restrita ao talude, tem a
    covariância certa, mas precisa de mais termos que a base do próprio domínio (\path{relatorio_kl_form}, placa curva).
  \item \textbf{Nyström no próprio talude} (novo, \path{Projects/Cho2010/nystrom_base.py}): quadratura de 2400 nós e
    $f_k(\bm x) = \lambda_k^{-1}\sum_j w_j C(\bm x, \bm x_j) f_k(\bm x_j)$. No \texttt{RandomField}, o campo sai de
    $Z(\bm x) = \sum_j C(\bm x, \bm x_j)\, a_j$, com $a_j = w_j \sum_k \xi_k f_k(\bm x_j)/\sqrt{\lambda_k}$ calculado uma vez
    por realização.
\end{itemize}

\subsection{Monte Carlo direto e amostragem por importância}
\begin{itemize}
  \item \textbf{Direto:} $\hat P_f = \frac1N\sum I(\mathrm{FS}_i < 1)$, com CoV $= \sqrt{(1-P_f)/(N P_f)}$.
  \item \textbf{AI:} sorteia $\bm\xi \sim N(\bm\mu, I)$ e usa o peso exato
    $w = \exp(-\bm\xi\cdot\bm\mu + \tfrac12|\bm\mu|^2)$, com $\hat P_f = \frac1N \sum I_i w_i$. A estimativa não tem viés
    para qualquer $\bm\mu$ fixado antes; $\bm\mu$ só muda a variância.
  \item \textbf{Deslocamentos:}
    \begin{itemize}
      \item \emph{pela banda} (entropia cruzada): $\bm\mu = s\,\bm d$, com $\bm d$ a direção da banda de ruptura do
        determinístico e $s$ a média de $\bm\xi\cdot\bm d$ no conjunto elite de um piloto direto (as falhas, ou os 10\,\%
        de menor FS);
      \item \emph{pelo FORM}: $\bm\mu = \bm u^*$, o ponto de projeto.
    \end{itemize}
  \item \textbf{Diagnósticos:} o CoV, o ESS das falhas e a maior parcela de peso.
\end{itemize}

\subsection{FORM e os dois gradientes}
O FORM usa o algoritmo HL-RF nos $2M$ coeficientes $\bm\xi$ (os $M$ primeiros de cada campo), com $g = \mathrm{FS} - 1$.
\begin{itemize}
  \item \textbf{Diferenças centrais} (\texttt{form\_kl.py}): $4M + 1$ análises do SRM por iteração, com $h = 0{,}5$. Os
    termos fracos ficam dominados pelo ruído do FS (~0,001).
  \item \textbf{Projeção modal} (etapa 2c, \texttt{gradiente=} no Cho2010): o limite superior com Mohr--Coulomb associado
    e o teorema do envelope dão, no colapso do SRM ($c^* = c'/F$, $\tan\varphi^* = \tan\varphi'/F$),
    \begin{equation*}
      \frac{\partial D}{\partial c} = \dot\varepsilon_v\cot\varphi^*, \qquad
      \frac{\partial D}{\partial \varphi} = -\frac{\dot\varepsilon_v}{2\sin\varphi^*}\big[(\sigma_1+\sigma_3)\cos\varphi^* +
      2c^*\sin\varphi^*\big].
    \end{equation*}
    A primeira é a Eq.~14 do guia de estudo: a versão contínua da potência resistente da Eq.~48 do Ceron (2025). A
    segunda a estende ao $\varphi'$ com as tensões de colapso. Com $F$ definido por $\Gamma(c/F, \varphi^*) = 1$ (Ceron,
    Eq.~58), $\partial F/\partial c'(\bm x) = D_c/\mathrm{Den}$ e $\partial F/\partial\tan\varphi'(\bm x) =
    D_\varphi\cos^2\varphi^*/\mathrm{Den}$, com $\mathrm{Den} = \int (D_c c^* + D_\varphi\sin\varphi^*\cos\varphi^*)$. A
    projeção nos modos usa $\partial Z/\partial\xi_k = \sqrt{\lambda_k} f_k(\bm x)$ exata. Custa \textbf{uma} análise
    por gradiente, contra $4M + 1$. A taxa $\dot\varepsilon_v$ do mecanismo é aproximada pela deformação plástica
    volumétrica no colapso, e a escala se cancela.
\end{itemize}

\section{O plano de trabalho e a situação}
\begin{table}[H]\centering\small
""" + tabela("Etapa & O quê & Situação",
             ["0 & Cho2010 no fermi: direto, AI pela banda e pelo FORM & pronta \\\\",
              "1 & SRM corrigido ($F_0 = 0{,}8$): o antigo contra o novo nas mesmas realizações & pronta \\\\",
              "2 & FORM ($M = 10$, 20) e AI pelo ponto de projeto; teste linear sem MEF & pronta; repetições com sementes novas na fila \\\\",
              "2b & base de Nyström contra a analítica & FORM feito no PC (Seção~\\ref{sec:nystrom}) \\\\",
              "2c & gradiente do FORM pela projeção modal & validado no PC; AI com o ponto modal rodando \\\\",
              "3 & $P_f$ pequena ($r = -0{,}7$), tamanho do piloto, custo & pronta; repetições na fila \\\\",
              "4 & FORM e AI na chuva & referência pronta (direto 4000); falta o código do FORM na chuva \\\\",
              "5 & investigação da chuva & não começou \\\\"], "lp{0.55\\textwidth}p{0.3\\textwidth}") + r"""
\end{table}

\section{Testes e resultados}
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
\caption{$r = -0{,}5$ (Cho: $P_f = 1{,}71\,\%$). O direto ainda está rodando (meta: 24\,000). ``dif.'' é a diferença AI
$-$ direto em desvios. As linhas com ``Tri(0)'' e ``modal'' são do PC.}
\end{table}
\begin{figure}[H]\centering\includegraphics[width=0.75\textwidth]{figs/pf_r05.pdf}\caption{$P_f$ e IC 95\,\%, $r = -0{,}5$.}\end{figure}

\subsection{$P_f$ pequena, $r = -0{,}7$ (etapa 3)}
\begin{table}[H]\centering\small
""" + T_R07 + r"""
\caption{$r = -0{,}7$ (Cho: $P_f \approx 0{,}4\,\%$).}
\end{table}

\subsection{Custo total para CoV$(P_f) = 10\,\%$}
\begin{table}[H]\centering\small
""" + T_CUSTO + r"""
\caption{Análises de elementos finitos, com o piloto (1000) e o FORM incluídos. Com o piloto na malha grossa (custo
1/20, Seção~\ref{sec:hip}), a AI pela banda cai para ~1300 e ~2700.}
\end{table}

\subsection{O teste linear (sem elementos finitos)}
\begin{table}[H]\centering\small
""" + T_LIN + r"""
\caption{Falha $\bm\alpha\cdot\bm\xi > \beta$, com $P_f$ exata. Erro relativo, com o erro em desvios entre parênteses.}
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
coincidem a partir de $M = 50$. É o mesmo efeito da placa curva do \path{relatorio_kl_form}. O $\beta$ converge para
$\approx 2{,}00$, mas a $P_f$ do FORM sozinho ($\approx 2{,}3\,\%$) superestima a do direto em ~30\,\%: é o erro de
linearizar uma superfície de falha curva. Por isso o FORM entra só como deslocamento da AI.

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
  \item \textbf{Há um só modo de ruptura} no Cho2010: todas as falhas do direto têm $t = \bm\alpha\cdot\bm\xi \ge 1{,}6$, e
    a banda de ruptura das falhas é a mesma no direto e na AI.
  \item \textbf{O piloto pode ser pequeno:} com 250 realizações, o deslocamento da banda custa só ~5\,\% de variância.
  \item \textbf{O piloto pode rodar em malha grossa:} o elite é escolhido por ordem, e a ordem das realizações quase não
    muda. No nível 0 da tri2 (1/20 do custo), a perda é de 1,5 a 4\,\% de variância.
  \item \textbf{O gradiente modal dá a direção do FORM com 1 análise,} a menos de ~9$^\circ$ das diferenças centrais, e o
    FORM converge em 3 análises.
  \item \textbf{A base de Nyström precisa de menos termos} que a analítica para o mesmo $\beta$.
\end{itemize}
\subsection*{Descartadas}
\begin{itemize}
  \item \textbf{O FORM na malha Tri(0)} (76 equações, FS 1,47 contra 1,20): $\beta = 4{,}08$ e direção a 31,7$^\circ$ da
    certa. O erro da malha desloca o ponto de projeto, e não só o FS. Com a direção da Tri(0) e a distância da malha usual
    ($\beta = 2{,}32$), a AI funciona, mas pior que o FORM na malha usual.
\end{itemize}
\subsection*{Sob investigação}
\begin{itemize}
  \item \textbf{A AI pelo FORM sai um pouco abaixo do direto} (de $-0{,}9$ a $-1{,}5$ desvio nos quatro casos). Explicação
    provável, pelo teste linear: o ponto de projeto fica um pouco fora da direção ótima, e com $N$ finito a estimativa
    típica fica baixa. Testes em andamento:
    \begin{itemize}
      \item repetições com sementes novas, independentes (fermi, na fila);
      \item AI com o ponto do gradiente modal, de direção mais certa (PC, rodando).
    \end{itemize}
  \item \textbf{A AI pela banda com o piloto na Tri(0)} (PC, rodando): a parcial está abaixo do direto.
  \item \textbf{O campo de $k$ na chuva:} o gradiente modal pega $c'$ e $\varphi'$; o $k$ age pelo transiente e
    precisaria do adjunto. Hoje o $k$ pesa pouco na chuva f50 ($\alpha^2_k \approx 0$), e a AI continua sem viés sem
    deslocá-lo.
\end{itemize}

\section{Próximos passos}
\begin{enumerate}
  \item Fechar a questão do pequeno desvio da AI pelo FORM com as repetições e com a AI pelo ponto modal.
  \item Etapa 4: o FORM com gradiente modal e a AI pelo ponto de projeto na chuva, contra o direto de 4000.
  \item AI com a base de Nyström e o seu próprio direto de referência (a mesma semente dá outras realizações).
  \item Para o artigo da AI: o ganho contra o custo total (FORM e piloto incluídos), de $P_f \approx 2\,\%$ a $0{,}4\,\%$.
\end{enumerate}

\section{Mapa dos relatórios}\label{sec:mapa}
Os documentos estão espalhados; esta lista diz onde está cada um e o que tem. Os caminhos relativos são a partir de
\path{~/projects/neopz-master}.
""" + MAPA + r"""
\end{document}
"""
open(f"{AQUI}/relatorio_parcial.tex", "w").write(tex)
print("ok")
