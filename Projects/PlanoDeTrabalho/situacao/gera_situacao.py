#!/usr/bin/env python3
"""gera_situacao.py: o relatorio de situacao: o plano de trabalho e a situacao de cada etapa, os resultados principais, as
hipoteses, os proximos passos e o mapa dos documentos. Todos os resultados, com as tabelas e as figuras, ficam no
relatorio das etapas (../ETAPAS1a3/relatorio_etapas.pdf): este arquivo roda o gerador dele (relatorio_etapas.py), usa os
mesmos numeros e cita as secoes dele (E-) e as da formulacao (F-) pelo .aux (pacote xr). O mapa sai da lista unica
../documentos.py, a mesma do indice do README. Para atualizar os tres PDFs:
    Projects/PlanoDeTrabalho/situacao/atualiza.sh
"""
import os
import runpy

AQUI = os.path.dirname(os.path.abspath(__file__))  # Projects/PlanoDeTrabalho/situacao
PLANO = os.path.dirname(AQUI)
ETAPAS = f"{PLANO}/ETAPAS1a3"

# todos os numeros: o gerador do relatorio das etapas (grava relatorio_etapas.tex e as figuras dele)
cwd = os.getcwd()
os.chdir(ETAPAS)
E = runpy.run_path(f"{ETAPAS}/relatorio_etapas.py")
os.chdir(cwd)
pct, num, ganho, desvios, mcnemar, tabela = E["pct"], E["num"], E["ganho"], E["desvios"], E["mcnemar"], E["tabela"]
d05, a05, d07, a07, dch, ach, cho, chuva = (E[k] for k in ("d05", "a05", "d07", "a07", "dch", "ach", "cho", "chuva"))
DD5, DD7, AM, AT0, A07, PP, FGIM, GV0, GVU = (E[k] for k in ("DD5", "DD7", "AM", "AT0", "A07", "PP", "FGIM", "GV0", "GVU"))
S05, S07, F4, GANHOS, GANHO_K50, LOTES = (E[k] for k in ("S05", "S07", "F4", "GANHOS", "GANHO_K50", "LOTES"))
DIF_T0, DIFS05, CUSTO05, CUSTO07, f05 = (E[k] for k in ("DIF_T0", "DIFS05", "CUSTO05", "CUSTO07", "f05"))
SRM_MOD, SRM_CD, angulo = E["SRM_MOD"], E["SRM_CD"], E["angulo"]
F4KL = [a for r, a in F4 if "M = 100" in r]
F4NY = [a for r, a in F4 if "M = 200" in r]
OTIMISMO = [100 * (1 - LOTES[k]["desvio_estimado"] / LOTES[k]["desvio"]) for k in LOTES if k != "P"]


def mil(n):
    return f"{n:,}".replace(",", "\\,")


# o mapa dos documentos: a mesma lista do indice do README (../documentos.py)
MAPA = runpy.run_path(f"{PLANO}/documentos.py", run_name="documentos")["mapa_latex"]()

HIP05 = (r"  \item \textbf{A AI pelo FORM não tem viés} ($r = -0{,}5$), nas duas bases (etapas, Seção~\ref{E-sec:vies}). "
         r"Com a semente 3, as AI saem todas em ≈1,55\,\% (os sorteios são os mesmos). As " + str(S05["k"]) +
         r" sementes independentes da KL dão, em média, " + pct(S05["p"], 3) + r"\,\% contra " + pct(DD5["p"], 3) +
         r"\,\% do direto e, com a variável de controle, " + pct(S05["pc"], 3) + r" $\pm$ " + pct(S05["sdc"], 3) +
         r"\,\% contra " + pct(DD5["cv"]["p"], 3) + r" $\pm$ " + pct(DD5["cv"]["sd"], 3) + r"\,\% ($" + num(S05["difc"], 1) +
         r"$ desvio). A direção errada da Tri(0) é a única que subestima o próprio CoV (" + num(100 * AT0["sd"] / AT0["p"], 1) +
         r" contra " + num(100 * AT0["cov_verdadeiro"], 1) + r"\,\%); com o CoV verdadeiro, o seu desvio cai de $-2{,}2$ "
         r"para $" + num(DIF_T0, 1) + "$.")
HIP07 = (r"  \item \textbf{A AI pelo FORM com $r = -0{,}7$ também não tem viés detectável:} as sementes 3, 4 e 5, corrigidas, "
         r"dão " + ", ".join(pct(a["p_cv"], 3) for a in A07) + r"\,\%, contra " + pct(DD7["cv"]["p"], 3) + r" $\pm$ " +
         pct(DD7["cv"]["sd"], 3) + r"\,\% do direto corrigido (" + mil(DD7["n"]) + " realizações, " + str(DD7["falhas"]) +
         " falhas): $" + (num(S07["difc"], 1) if S07 else "---") + "$ desvio.")
FASE4 = ("; ".join(f"{pct(a['p'])} $\\pm$ {pct(a['sd'])}\\,\\% ({a['pastas']} sementes, {a['n']} realizações) "
                   f"{'na KL' if 'M = 100' in r else 'em Nyström'}" for r, a in F4 if a["pastas"] > 1))

tex = (r"""\documentclass[11pt,a4paper]{article}
\usepackage[utf8]{inputenc}
\DeclareUnicodeCharacter{2248}{\ensuremath{\approx}}
\usepackage[T1]{fontenc}
\usepackage[brazil]{babel}
\usepackage{lmodern}
\usepackage[margin=2cm]{geometry}
\usepackage{amsmath,amssymb,bm}
\usepackage{graphicx,booktabs,float,longtable,url}
\usepackage{xr-hyper}
\usepackage[hidelinks]{hyperref}
% as secoes da formulacao (F-) e do relatorio das etapas (E-): compile os dois antes (o atualiza.sh faz isso)
\externaldocument[F-]{../../../docs/confiabilidade_campos/formulacao/formulacao}
\externaldocument[E-]{../ETAPAS1a3/relatorio_etapas}
% a data e o commit vem do atualiza.sh; compilado a mao, fica so a data
\providecommand{\versao}{17/11/2026}
\title{Confiabilidade de taludes com campos aleatórios no NeoPZ\\\large Situação: plano de trabalho, resultados
principais, hipóteses e próximos passos}
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
na chuva. As etapas 0 a 3 estão rodadas, e a 2c (gradiente modal) foi testada no PC, no SRM e no GIM. A etapa 4 (FORM e
AI na chuva) tem a referência pronta, e a etapa 5 (a investigação da chuva) ainda não começou. Tudo roda no cluster
fermi (4 nós, 48 processos); o PC serve para testes curtos.

Três documentos:
\begin{itemize}
  \item a \textbf{formulação}, citada como ``formulação, Seção $n$'':\newline
    \path{docs/confiabilidade_campos/formulacao/formulacao.pdf};
  \item o \textbf{relatório das etapas}, com todos os resultados, as tabelas e as figuras, citado como ``etapas, Seção
    $n$'':\newline \path{Projects/PlanoDeTrabalho/ETAPAS1a3/relatorio_etapas.pdf};
  \item este, com a situação, o resumo dos resultados, as hipóteses e os próximos passos. A Seção~\ref{sec:mapa} lista
    todos os documentos.
\end{itemize}

\section{O plano de trabalho e a situação}
\begin{table}[H]\centering\small
""" + tabela("Etapa & O quê & Situação",
             ["0 & Cho2010 no fermi: direto, AI pela banda e pelo FORM & pronta \\\\",
              "1 & SRM corrigido ($F_0 = 0{,}8$): o antigo contra o novo nas mesmas realizações & pronta (etapas, "
              "Seção~\\ref{E-sec:srm}) \\\\",
              "2 & FORM ($M = 10$, 20) e AI pelo ponto de projeto; teste linear sem MEF & pronta; sem viés (etapas, "
              "Seção~\\ref{E-sec:vies}) \\\\",
              "2b & base de Nyström contra a analítica & o mesmo ponto de projeto (etapas, Seção~\\ref{E-sec:pp}); fase 4 "
              "de Nyström parcial \\\\",
              "2c & gradiente do FORM pela projeção modal & validado no SRM e no GIM (etapas, "
              "Seção~\\ref{E-sec:validacao}); fase 4 parcial \\\\",
              f"3 & $P_f$ pequena ($r = -0{{,}}7$), tamanho do piloto, custo & o direto vai "
              f"de {mil(d07['n'])} ({d07['k']} falhas) a 110\\,000 (CoV 5\\,\\%), com os 4 nós do fermi \\\\",
              f"4 & FORM e AI na chuva & o direto de {mil(dch['n'])} terminou ($P_f = {pct(dch['p'])}\\,\\%$); falta o "
              "código do FORM na chuva \\\\",
              "5 & investigação da chuva & não começou \\\\"], "lp{0.5\\textwidth}p{0.35\\textwidth}") + r"""
\end{table}

\section{Resultados principais}\label{sec:resultados}
Cada item cita a seção do relatório das etapas com o resultado completo.
\begin{itemize}
  \item \textbf{SRM corrigido} (etapas, Seção~\ref{E-sec:srm}). No Cho2010, as falhas são as mesmas nas duas versões
    (""" + f"{cho['fa']} e {cho['fb']} em {cho['n']}" + r"""); na chuva, """ + f"{chuva['n01']}" + r""" realizações passam a
    falhar e """ + f"{chuva['n10']}" + r""" deixa de falhar (McNemar $p = """ + num(mcnemar(chuva['n01'], chuva['n10']), 2) +
    r"""$). O efeito é pequeno, e o SRM novo, o correto, fica.
  \item \textbf{AI contra o direto, $r = -0{,}5$} (etapas, Seção~\ref{E-sec:aidireto}). O direto de """ + mil(d05['n']) +
    r""" dá $P_f = """ + pct(d05['p']) + r"""\,\%$ (Cho: 1,71\,\%). As AI completas ficam a menos de """ +
    num(max(DIFS05), 1) + r""" desvio dele, com ganho de até """ + num(max(ganho(a, d05) for _, a in a05), 0) +
    r"""$\times$; só a direção errada da malha Tri(0) fica mais longe.
  \item \textbf{$P_f$ pequena, $r = -0{,}7$} (etapas, Seção~\ref{E-sec:r07}). O direto de """ + mil(d07['n']) + r""" dá
    $P_f = """ + pct(d07['p']) + r"""\,\%$ (""" + f"{d07['k']}" + r""" falhas; Cho: ≈0,4\,\%). As AI pelo FORM acertam,
    com ganho de até """ + num(max(ganho(a, d07) for _, a in a07), 0) + r"""$\times$.
  \item \textbf{FORM e as duas bases} (etapas, Seções~\ref{E-sec:nystrom} e~\ref{E-sec:pp}). O $\beta$ converge para
    $\approx 2{,}00$ nas duas bases, e a de Nyström chega lá com menos termos. O FORM sozinho superestima a $P_f$ em
    ≈30\,\%: ele serve para o deslocamento da AI. O ponto de projeto é o mesmo nas duas bases (a """ +
    num(PP["angulo"], 1) + r"""$^\circ$).
  \item \textbf{O gradiente modal} (etapas, Seção~\ref{E-sec:validacao}). Uma análise por iteração, no lugar de $4M + 1$.
    No SRM, fica a 9 a """ + num(angulo(SRM_MOD, SRM_CD), 0) + r"""$^\circ$ das diferenças centrais; no GIM, a """ +
    num(GV0["angulo_graus"], 1) + " e " + num(GVU["angulo_graus"], 1) + r"""$^\circ$. O FORM do GIM converge em 5 análises
    para $\beta = """ + num(FGIM[10]["beta_corrigido"], 3) + r"""$ ($M = 10$) e $""" + num(FGIM[20]["beta_corrigido"], 3) +
    r"""$ ($M = 20$), contra """ + num(f05[0]["bc"], 3) + " e " + num(f05[1]["bc"], 3) + r""" das diferenças centrais no SRM.
  \item \textbf{Sem viés} (etapas, Seção~\ref{E-sec:vies}). As falhas do direto ficam todas onde a AI amostra, a
    variância verdadeira é a estimada, e as sementes independentes acertam o direto. A fase 4 parcial (etapas,
    Seção~\ref{E-sec:fase4}) dá """ + FASE4 + r""", contra """ + pct(d05['p']) + r"""\,\% do direto.
  \item \textbf{Variável de controle} (etapas, Seção~\ref{E-sec:cv}). Com o substituto ajustado no direto inteiro, de
    """ + num(min(GANHOS), 1) + " a " + num(max(GANHOS), 1) + r""" vezes menos realizações; com o ajustado na própria
    amostra (o que se tem na prática), de """ + num(min(GANHO_K50), 1) + " a " + num(max(GANHO_K50), 1) + r""". Sem viés,
    mas o desvio que ela estima saiu de """ + num(min(OTIMISMO), 0) + " a " + num(max(OTIMISMO), 0) + r"""\,\% menor que o
    verdadeiro, nos lotes do direto.
  \item \textbf{Piloto e custo} (etapas, Seção~\ref{E-sec:custo}). Um piloto de 250 basta, e ele pode rodar em malha
    grossa; o FORM, não. Para CoV de 10\,\%, a AI pelo ponto de projeto custa de """ +
    num(CUSTO05[0][3] / CUSTO05[2][3], 0) + " a " + num(CUSTO07[0][3] / CUSTO07[3][3], 0) + r""" vezes menos análises que
    o direto, já contando o FORM.
  \item \textbf{A referência da chuva} (etapas, Seção~\ref{E-sec:chuva}). O direto de """ + mil(dch['n']) + r""" dá
    $P_f = """ + pct(dch['p']) + r"""\,\%$, e a AI pela banda concorda com ele (""" + num(desvios(ach[0][1], dch), 1) +
    r""" desvio).
\end{itemize}

\section{Hipóteses}\label{sec:hip}
\subsection*{Confirmadas}
\begin{itemize}
  \item A AI não tem viés com qualquer deslocamento fixado antes; o deslocamento só muda a variância (teste linear e Ceron).
  \item O SRM antigo tinha um buraco logo abaixo de FS = 1; o efeito na $P_f$ é pequeno.
  \item \textbf{Há um só modo de ruptura} no Cho2010, nas duas bases: todas as falhas do direto ficam onde a AI pelo
    FORM amostra ($t = \bm\xi\cdot\bm\xi^*/\beta \ge """ + num(AM["t_min"], 2) + r"""$ na direção modal, com só
    """ + str(AM["abaixo_15"]) + r""" das """ + str(DD5["falhas"]) + r""" abaixo de 1,5;
    $t \ge """ + num(A07[0]["t_min"], 2) + r"""$ com $r = -0{,}7$), e a variância verdadeira da AI,
    calculada pelas falhas do direto, é a que a AI estima (etapas, Seção~\ref{E-sec:diag}). A banda de ruptura das falhas
    é a mesma no direto e na AI.
  \item \textbf{O ponto de projeto é o mesmo nas duas bases:} levado para as coordenadas de Nyström, o ponto analítico
    fica a """ + num(PP["angulo"], 1) + r"""$^\circ$ do de Nyström, e os campos de $\varphi'$ diferem em
    menos de 0,1$^\circ$ (etapas, Figura~\ref{E-fig:pp}).
  \item \textbf{A variável de controle dá um ganho extra sem custo:} de """ + num(min(GANHOS), 1) + r""" a
    """ + num(max(GANHOS), 1) + r""" nas AI de $r = -0{,}5$ com o substituto do direto, e de """ +
    num(min(GANHO_K50), 1) + " a " + num(max(GANHO_K50), 1) + r""" com o da própria amostra. No direto, a referência vai
    a """ + pct(DD5["cv"]["p"], 3) + r""" $\pm$ """ + pct(DD5["cv"]["sd"], 3) + r"""\,\%.
  \item \textbf{O piloto pode ser pequeno:} com 250 realizações, o deslocamento da banda custa só ≈5\,\% de variância.
  \item \textbf{O piloto pode rodar em malha grossa:} o elite é escolhido por ordem, e a ordem das realizações quase não
    muda. Com o piloto na Tri(0) (300 vezes mais barata), a AI pela banda acertou o direto ($-0{,}2$ desvio), com o mesmo
    CoV do piloto na malha usual.
  \item \textbf{O gradiente modal dá a direção do FORM com 1 análise,} a menos de ≈9 a 11$^\circ$ das diferenças centrais
    no SRM e de 3 a 6$^\circ$ no GIM, e o FORM converge em 3 a 5 análises.
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
  \item \textbf{O desvio estimado pela variável de controle} saiu otimista nos lotes do direto (etapas,
    Seção~\ref{E-sec:cvlocal}); as sementes da fase 4 vão dizer se é acaso.
\end{itemize}

\section{Próximos passos}
\begin{enumerate}
  \item \textbf{Direto de $r = -0{,}7$ até 110\,000} (CoV de 5\,\%; com os 4 nós do fermi). Está em """ + mil(d07['n']) +
    r""" realizações. A fase 4 fica pausada até ele acabar.
  \item \textbf{Fase 4 no fermi} (≈8\,h): 10 sementes do FORM modal analítico e 10 do de Nyström, as mesmas nas duas
    bases. Dão a distribuição do estimador e fecham a questão do viés (etapas, Seção~\ref{E-sec:vies}). As partes que
    rodaram antes do cancelamento já estão juntas (etapas, Seção~\ref{E-sec:fase4}).
  \item \textbf{Etapa 4:} o FORM com gradiente modal e a AI pelo ponto de projeto na chuva, contra o direto de
    """ + mil(dch['n']) + r""".
  \item \textbf{Variável de controle realista} (etapas, Seção~\ref{E-sec:cvlocal}): o substituto ajustado na própria
    amostra ($K = 50$, 5 dobras) em todas as AI, inclusive pela banda e na chuva; conferir o desvio estimado com as 10
    sementes da fase 4. Para o artigo da AI: o ganho contra o custo total (FORM, piloto e controle), de
    $P_f \approx 2\,\%$ a $0{,}4\,\%$.
  \item \textbf{O FORM da chuva pelo GIM?} O gradiente modal do GIM é mais preciso que o do SRM (""" +
    num(GVU["angulo_graus"], 1) + " contra " + num(angulo(SRM_MOD, SRM_CD), 1) + r"""$^\circ$ no mesmo ponto), e a região de
    falha é a mesma. Testar se vale na chuva, com as pressões do transiente.
\end{enumerate}

\section{Mapa dos documentos}\label{sec:mapa}
É a mesma lista do índice do \path{Projects/PlanoDeTrabalho/README.md}: as duas saem de
\path{Projects/PlanoDeTrabalho/documentos.py}. Os caminhos relativos são a partir de \path{~/projects/neopz-master}.
""" + MAPA + r"""
\end{document}
""")
open(f"{AQUI}/situacao.tex", "w").write(tex)
print("ok: situacao.tex")
