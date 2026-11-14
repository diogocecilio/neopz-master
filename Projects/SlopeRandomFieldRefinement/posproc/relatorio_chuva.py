#!/usr/bin/env python3
"""relatorio_chuva.py: o relatorio (LaTeX e PDF) do estudo da chuva: o Monte Carlo do talude de Rahardjo et al. (2007)
com os campos aleatorios de Vargas Ceron et al. (2025), nas Secoes 6.1 a 6.5 do Ceron (a 6.5 com a intensidade da
chuva no lugar do rebaixamento), como o relatorio_sat.py do rebaixamento (SlopeRandomField/comparacao_ceron), mais a
duracao da chuva, FS antes x depois da chuva, as categorias do nivel d'agua (Tabela 3 do Rahardjo) e o mapa da banda.
So os fatos: valores, condicoes e o que falta.

Le o JSON do dados_chuva.py, faz as figuras (figuras_chuva.py) em <pasta>/figs e escreve <pasta>/relatorio.tex; compila
com pdflatex se ele existir.

Uso: python3 -I relatorio_chuva.py <estatisticas.json> <pasta>"""

import json
import math
import os
import shutil
import subprocess
import sys

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)
import figuras_chuva  # noqa: E402

CATEGORIAS = figuras_chuva.CATEGORIAS
NOME_TEX = {"covk": r"$CoV(k_s)$", "covc": r"$CoV(c')$", "covphi": r"$CoV(\varphi')$", "s": r"$s$",
            "alfa": r"$\alpha$", "I": r"$I$ (mm/h)"}
SECAO = {"covk": "6.2", "covc": "6.2", "covphi": "6.2", "s": "6.3", "alfa": "6.4", "I": "6.5"}
TEXTO_CATEGORIA = {"I": "o nível chega à superfície na crista ($x = 30$~m)",
                   "II": "chega à superfície no meio da face ($x = 35$~m), não na crista",
                   "III": "chega à superfície só no pé ($x = 40$~m)",
                   "IV": "não chega à superfície em nenhuma das três verticais (sem subida ou subida pequena)"}


def n(v, casas=3):
    if v is None or (isinstance(v, float) and not math.isfinite(v)):
        return "--"
    return ("%.*f" % (casas, v)).replace(".", "{,}")


def g(v):
    """numero curto (%g) com virgula"""
    if v is None:
        return "--"
    return ("%g" % v).replace(".", "{,}")


def num_tex(texto):
    """numero do parametros.txt para dentro de $...$: 1e-05 -> 10^{-5}, 2.5e-06 -> 2{,}5 \\times 10^{-6}, 0.3 -> 0{,}3"""
    try:
        v = float(texto)
    except (TypeError, ValueError):
        return r"\text{%s}" % tex_seguro(texto if texto is not None else "--")
    if v != 0 and (abs(v) < 1e-3 or abs(v) >= 1e5):
        e = int(math.floor(math.log10(abs(v))))
        m = v / 10 ** e
        if abs(m - 1) < 1e-9:
            return r"10^{%d}" % e
        return r"%s \times 10^{%d}" % (("%g" % m).replace(".", "{,}"), e)
    return g(v)


def lista_tex(texto):
    """'0.0000,24.000' -> '0, 24'"""
    try:
        return ", ".join(g(float(x)) for x in str(texto).split(",") if x.strip())
    except ValueError:
        return tex_seguro(texto)


def caminho_tex(c):
    """caminho em \\texttt que pode quebrar depois das barras"""
    return tex_seguro(c).replace("/", "/\\allowbreak{}")


def pct(v, casas=2):
    return "--" if v is None or not math.isfinite(v) else n(100. * v, casas)


def ic(e, casas=1):
    a, b = e["IC95_Pf"]
    return "[%s; %s]" % (pct(a, casas), pct(b, casas))


def tex_seguro(t):
    """texto qualquer como texto do LaTeX (um passe so, para nao reescapar o que ja foi escapado)"""
    mapa = {"\\": r"\textbackslash{}", "_": r"\_", "%": r"\%", "&": r"\&", "#": r"\#", "$": r"\$", "{": r"\{",
            "}": r"\}", "~": r"\textasciitilde{}", "^": r"\textasciicircum{}"}
    return "".join(mapa.get(ch, ch) for ch in str(t))


def origem_tex(o):
    """a origem dos niveis iniciais (texto do dados_chuva.py) no texto do relatorio"""
    if not o:
        return "--"
    o = o.replace("deterministico", "o determinístico").replace("realizacao", "a realização")
    o = o.replace("referencia (mesmo nivel inicial)", "os da referência, com o mesmo nível inicial")
    return tex_seguro(o).replace("t = 0", "$t = 0$")


def figura(arq, legenda, rotulo, largura=r"\linewidth"):
    return "\n".join([r"\begin{figure}[htbp]\centering",
                      r"\includegraphics[width=%s,height=0.8\textheight,keepaspectratio]{%s}" % (largura, arq),
                      r"\caption{%s}\label{%s}" % (legenda, rotulo), r"\end{figure}"])


def fig_rel(figs, nome):
    return "figs/%s.pdf" % nome if nome in figs else None


def rotulo_caso(c):
    if c["familia"] == "ref":
        return "ref"
    return "%s = %s" % (NOME_TEX[c["familia"]].replace(" (mm/h)", ""), g(c["valor"]))


def cats(nv):
    if not nv:
        return "--"
    return "/".join("%d" % round(100. * (nv["fracao"].get(c) or 0.)) for c in CATEGORIAS)


# ---------------------------------------------------------------- tabelas
def tabela_casos(est):
    L = [r"\begin{table}[htbp]\centering\scriptsize\setlength{\tabcolsep}{3pt}",
         r"\caption{Casos lidos: $S$ realizações usadas no instante das famílias e excluídas nele (situação diferente "
         r"de \texttt{ok} ou $FS \le 0$), instantes calculados (h) e o determinístico (próprio, da referência ou "
         r"ausente).}",
         r"\begin{tabular}{llrrrll}\toprule",
         r"caso & família & $S$ & excluídas & realizações & instantes (h) & determinístico \\\midrule"]
    for c in est["casos"]:
        e = c.get("estat") or {}
        det = c.get("det") or {}
        L.append(" & ".join([tex_seguro(c["id"]), tex_seguro(c["familia"]), str(e.get("S", 0)),
                             str(c.get("S_excluidas_t", 0)), str(c["realizacoes"]),
                             ", ".join(g(t) for t in c["instantes_h"]),
                             {"proprio": "próprio", "ref": "da referência"}.get(det.get("origem"), "--")]) + r" \\")
    L += [r"\bottomrule\end{tabular}\end{table}"]
    return "\n".join(L)


def tabela_det(est, ref):
    det = ref.get("det") or {}
    rah = est["rahardjo"]
    fs_t = {round(t, 6): fs for t, fs in det.get("FS_t", [])}
    nv_t = {round(t, 6): (nv, cat) for t, nv, cat in det.get("niveis_t", [])}
    L = [r"\begin{table}[htbp]\centering\footnotesize\setlength{\tabcolsep}{4pt}",
         r"\caption{Determinístico (médias, sem campos) da referência: $FS$ e cota do nível d'água (m) nas verticais "
         r"da crista, do meio da face e do pé (superfície a 30, 25 e 20~m), com a categoria da Tabela~3 do Rahardjo. "
         r"Rahardjo et al.\ (2007), solo f50 a $45^\circ$: $F_{ini} = %s$ (Tabela~7) e, com $I = %s$~mm/h, "
         r"$F_{min} = %s$ (sigmoide da Fig.~7b).}" % (n(rah["F_ini"]), g(est["valor_ref"]["I"]), n(rah["F_min_ref"])),
         r"\begin{tabular}{rrrrrc}\toprule",
         r"$t$ (h) & $FS$ & nível em $x = 30$ & $x = 35$ & $x = 40$ & categoria \\\midrule"]
    for t in sorted(set(fs_t) | set(nv_t)):
        nv, cat = nv_t.get(t, ({}, None))
        L.append(" & ".join([g(t), n(fs_t.get(t), 4), n(nv.get("x30"), 1), n(nv.get("x35"), 1), n(nv.get("x40"), 1),
                             cat or "--"]) + r" \\")
    L += [r"\bottomrule\end{tabular}\end{table}"]
    return "\n".join(L)


def tabela_ref(ref):
    e = ref["estat"]
    q2n, q2l = e.get("qui2_normal"), e.get("qui2_lognormal")
    linhas = [
        ("$S$", str(e["S"])), (r"$n_f = \#(FS < 1)$", str(e["n_f"])), (r"$P_f$ (\%)", pct(e["P_f"])),
        (r"IC de 95\% de $P_f$ (Wilson, \%)", ic(e, 2)), (r"$CoV(P_f)$ (\%)", pct(e["CoV_Pf"], 1)),
        (r"$P_f$ do ajuste log-normal (\%)", pct(e["P_f_lognormal"])),
        (r"$\mu$", n(e["media"]) + " (" + n(e["media_amostral"]) + ")"), ("mediana", n(e["mediana"])),
        (r"$\sigma$", n(e["desvio"]) + " (" + n(e["desvio_amostral"]) + ")"),
        (r"$CoV$ (\%)", pct(e["CoV"], 1) + " (" + pct(e["CoV_amostral"], 1) + ")"),
        (r"mínimo; máximo", n(e["minimo"]) + "; " + n(e["maximo"])),
        (r"quantis de 5\% e 95\%", n(e["q05"]) + "; " + n(e["q95"])),
        ("qui-quadrado, $p$ da normal", n(q2n["p"]) if q2n else "--"),
        ("qui-quadrado, $p$ da log-normal", n(q2l["p"]) if q2l else "--"),
    ]
    det = ref.get("det") or {}
    if det.get("FS"):
        linhas.insert(0, ("$FS$ determinístico", n(det["FS"], 4)))
    L = [r"\begin{table}[htbp]\centering\footnotesize",
         r"\caption{Referência em $t = %s$~h. $\mu$, mediana, $\sigma$ e $CoV$ da log-normal ajustada (máxima "
         r"verossimilhança), como no Ceron; entre parênteses, os amostrais. Qui-quadrado com $k = \lceil 2 S^{0,4} "
         r"\rceil$ classes equiprováveis.}" % g(ref["t"]),
         r"\begin{tabular}{lr}\toprule"]
    L += ["%s & %s \\\\" % (a, b) for a, b in linhas]
    L += [r"\bottomrule\end{tabular}\end{table}"]
    return "\n".join(L)


def tabela_spearman(ref):
    blocos = []
    i0 = next((i for i in ref.get("instantes", []) if abs(i["t"]) < 1e-6), None)
    if i0 and abs(ref["t"]) > 1e-6 and i0.get("spearman_banda"):
        blocos.append((0.0, i0["spearman_banda"], i0["spearman_campo"]))
    if ref.get("spearman_banda"):
        blocos.append((ref["t"], ref["spearman_banda"], ref["spearman_campo"]))
    if not blocos:
        return None
    L = [r"\begin{table}[htbp]\centering\footnotesize",
         r"\caption{Spearman parcial de 2ª ordem entre $FS$ e cada variável, controlando as outras duas (Eqs.~70 e 71 "
         r"do Ceron), com o intervalo de 95\% (transformação de Fisher). Duas definições das variáveis de cada "
         r"realização: $c'$ e $\varphi'$ médios na banda de ruptura e $k$ médio na massa que desliza; e as médias de "
         r"cada campo no solo todo.}",
         r"\begin{tabular}{llrr}\toprule", r"$t$ (h) & & banda e massa & médias do campo \\\midrule"]
    for t, sb, sc in blocos:
        for k, nome in (("k", r"$k_s$"), ("phi", r"$\varphi'$"), ("c", r"$c'$")):
            cel = []
            for s in (sb, sc):
                if s and s.get(k):
                    a, b = s[k]["ic95"]
                    cel.append("%s [%s; %s]" % (n(s[k]["rho_parcial"]), n(a, 2), n(b, 2)))
                else:
                    cel.append("--")
            L.append(" & ".join([g(t) if k == "k" else "", nome] + cel) + r" \\")
        L.append(r"\midrule" if t != blocos[-1][0] else "")
    L += [r"\bottomrule\end{tabular}\end{table}"]
    return "\n".join(L)


def tabela_duracao(ref):
    L = [r"\begin{table}[htbp]\centering\scriptsize\setlength{\tabcolsep}{3pt}",
         r"\caption{Referência em cada instante: o $FS$ em $t$ é o do fim de uma chuva de $t$ horas com a mesma "
         r"intensidade. $\mu$, $\sigma$ e $CoV$ da log-normal ajustada; categorias do nível d'água (I/II/III/IV, \% "
         r"das realizações).}",
         r"\begin{tabular}{rrrrrrrrrl}\toprule",
         r"$t$ (h) & $S$ & $FS_{det}$ & $\mu$ & $\sigma$ & $CoV$ (\%) & $P_f$ (\%) & IC 95\% (\%) & $P_f$ LN (\%) & "
         r"categorias (\%) \\\midrule"]
    for i in ref["instantes"]:
        e = i.get("estat")
        if not e:
            continue
        L.append(" & ".join([g(i["t"]), str(e["S"]), n(i.get("det_FS"), 4), n(e["media"]), n(e["desvio"]),
                             pct(e["CoV"], 1), pct(e["P_f"]), ic(e), pct(e["P_f_lognormal"]), cats(i["nivel"])])
                 + r" \\")
    L += [r"\bottomrule\end{tabular}\end{table}"]
    return "\n".join(L)


def tabela_nivel(est):
    L = [r"\begin{table}[htbp]\centering\scriptsize\setlength{\tabcolsep}{3pt}",
         r"\caption{Categorias do nível d'água em $t = %s$~h (\%% das realizações), a do determinístico, a subida "
         r"média do nível em cada vertical desde o nível inicial (m) e $P_f$ em cada categoria "
         r"($n_f$/realizações da categoria).}" % g(est["t_familias"]),
         r"\begin{tabular}{lrrrrcrrrl}\toprule",
         r"caso & I & II & III & IV & det. & subida $x = 30$ & $x = 35$ & $x = 40$ & $P_f$ por categoria "
         r"(I; II; III; IV) \\\midrule"]
    for c in est["casos"]:
        nv = c.get("nivel")
        if not nv or not sum(nv["contagem"].values()):
            continue
        sub = nv.get("subida_media") or {}
        pfc = "; ".join("%d/%d" % (nv["n_f"][k], nv["contagem"][k]) if nv["contagem"][k] else "--" for k in CATEGORIAS)
        L.append(" & ".join([rotulo_caso(c)] + [pct(nv["fracao"][k], 0) for k in CATEGORIAS] +
                            [nv.get("det_categoria") or "--", n(sub.get("x30"), 1), n(sub.get("x35"), 1),
                             n(sub.get("x40"), 1), pfc]) + r" \\")
    L += [r"\bottomrule\end{tabular}\end{table}"]
    return "\n".join(L)


def tabela_familia(fam, est):
    nossos = figuras_chuva.casos_familia(est, fam)
    vref = est["valor_ref"][fam]
    tem_det = fam in figuras_chuva.DET_PROPRIO
    cab = [NOME_TEX[fam]]
    if tem_det:
        cab += [r"$FS_{det}$"]
    if fam == "I":
        cab += [r"$F_{min}$ Rahardjo"]
    cab += [r"$\mu$", r"$\sigma$", r"$CoV$ (\%)", r"$P_f$ (\%)", r"IC 95\% (\%)", r"$P_f$ LN (\%)",
            r"$CoV(P_f)$ (\%)", "$S$", r"nível (\%)"]
    L = [r"\begin{table}[htbp]\centering\scriptsize\setlength{\tabcolsep}{2.6pt}",
         r"\caption{Análogo da Seção %s do Ceron, %s, em $t = %s$~h. $\mu$, $\sigma$ e $CoV$ da log-normal ajustada; "
         r"$P_f$ contado, com o intervalo de Wilson, e o do ajuste log-normal (LN); nível: categorias I/II/III/IV "
         r"da Tabela~3 do Rahardjo.}" % (SECAO[fam], NOME_TEX[fam], g(est["t_familias"])),
         r"\begin{tabular}{" + "r" * (len(cab) - 1) + "l}", r"\toprule", " & ".join(cab) + r" \\", r"\midrule"]
    for v, c in sorted(nossos.items()):
        e = c["estat"]
        row = [g(v) + (" (ref.)" if abs(v - vref) < 1e-9 else "")]
        if tem_det:
            row += [n((c.get("det") or {}).get("FS"), 4)]
        if fam == "I":
            row += [n(float(figuras_chuva.rahardjo_fmin(v, est["rahardjo"])))]
        row += [n(e["media"]), n(e["desvio"]), pct(e["CoV"], 1), pct(e["P_f"]), ic(e), pct(e["P_f_lognormal"]),
                pct(e["CoV_Pf"], 1), str(e["S"]), cats(c.get("nivel"))]
        L.append(" & ".join(row) + r" \\")
    L += [r"\bottomrule", r"\end{tabular}", r"\end{table}"]
    return "\n".join(L)


# ---------------------------------------------------------------- documento
def documento(est, figs, puladas):
    ref = figuras_chuva.ref_de(est)
    T = []
    A = T.append
    A(r"""\documentclass[11pt,a4paper]{article}
\usepackage[T1]{fontenc}\usepackage[utf8]{inputenc}\usepackage[brazil]{babel}
\usepackage[margin=2.2cm]{geometry}\usepackage{booktabs}\usepackage{graphicx}\usepackage{caption}
\usepackage{float}\usepackage[hidelinks]{hyperref}\usepackage{amsmath}\usepackage[section]{placeins}
\captionsetup{font=small}
\author{SlopeRandomFieldRefinement, \texttt{posproc}}
\title{Talude de Rahardjo et al.\ (2007) sob chuva com campos aleatórios de $c'$, $\varphi'$ e $k_s$: Monte Carlo no
formato da Seção~6 de Vargas Ceron et al.\ (2025)}
\date{}
\begin{document}\maketitle
""")
    if not ref:
        A(r"Não há a pasta \texttt{ref} nos dados lidos: o relatório precisa da configuração de referência.")
    p = (ref or {}).get("parametros", {})
    det = (ref or {}).get("det") or {}

    def par(k, padrao="--"):
        return tex_seguro(p.get(k, padrao))
    if ref:
        cm, phim = det.get("c_med"), det.get("phi_med")
        # instantes adaptativos (ADAPTA) e a malha de cada outro instante: a do adaptativo mais proximo (no empate, o
        # primeiro, como no montecarlo.cpp)
        try:
            adapta = [float(x) for x in p.get("ADAPTA", "").split(",") if x.strip()]
        except ValueError:
            adapta = []
        outros = []
        for t in ref["instantes_h"]:
            if adapta and not any(abs(t - a) < 1e-6 for a in adapta):
                perto = adapta[0]
                for a in adapta:
                    if abs(a - t) < abs(perto - t):
                        perto = a
                outros.append("%s~h na de %s~h" % (g(t), g(perto)))
        txt_outros = (" (" + "; ".join(outros) + ")") if outros else ""
        A(r"""\section{O que foi rodado}
O talude é o da Fig.~2 de Rahardjo et al.\ (2007): $H = 10$~m a $45^\circ$, domínio de 70~m $\times$ 30~m, solo
\texttt{%s} ($k_s = %s$~m/s, parâmetro $a = %s$~kPa de Fredlund--Xing), nível d'água inicial %s~m abaixo do pé e chuva de
$I = %s$~mm/h durante %s~h. Os campos de $c'$, $\varphi'$ e $k_s$ são log-normais (Karhunen--Loève com $M = %s$ termos,
%s\%% da variância retida), com médias %s~kPa e %s$^\circ$ (e $k_s$ do solo), $CoV(c') = %s$, $CoV(\varphi') = %s$,
$CoV(k_s) = %s$ e distâncias de autocorrelação $L_x = %s$~m e $L_y = %s$~m ($s = %s$; Tabela~2 do Ceron com $s = 1$);
$\alpha = k_h/k_v = %s$. A realização $n$ tem os mesmos números aleatórios em todos os casos.

Em cada realização: o transiente da chuva na malha \texttt{%s} ($u$ de ordem %s) até o último instante, e em cada
instante (%s~h) o fator de segurança por redução da resistência (\texttt{%s}) com as pressões do transiente, na malha
\texttt{%s} com %s níveis de refinamento adaptativo em $h$ feitos do zero nos instantes %s~h; os outros instantes usam a
malha refinada do instante adaptativo mais próximo%s. Nas realizações gravadas para figuras (\texttt{%s\_real<n>}),
todos os instantes são adaptativos. Falha: $FS < 1$ no instante. Custo médio por realização na referência: transiente
%s~s, fator %s~s.
""" % (par("SOLO"), num_tex(p.get("KS")), num_tex(p.get("A_FX")), num_tex(p.get("HW")), num_tex(p.get("I")),
       num_tex(p.get("CHUVA")), par("KL_M"), pct(ref.get("KL_variancia_media"), 1), g(cm), g(phim),
       num_tex(p.get("COV_C")), num_tex(p.get("COV_PHI")), num_tex(p.get("COV_K")), num_tex(p.get("LX")),
       num_tex(p.get("LY")), g(est["valor_ref"]["s"]), num_tex(p.get("ALFA")), par("MALHA_TRANSIENTE"), par("UT"),
       ", ".join(g(t) for t in ref["instantes_h"]), par("METODO").upper(), par("MALHA_FATOR"), par("NREF"),
       lista_tex(p.get("ADAPTA", "--")), txt_outros, tex_seguro(est["prefixo"]),
       n(ref["custo_medio_s"]["transiente"], 1), n(ref["custo_medio_s"]["fator"], 1)))
    A(tabela_casos(est))

    if ref:
        A(r"\section{Determinístico (análogo da Fig.~24)}")
        if det:
            A(tabela_det(est, ref))
        if fig_rel(figs, "fig24_det"):
            A(figura(fig_rel(figs, "fig24_det"), r"Talude da Fig.~2 do Rahardjo, banda de ruptura do determinístico "
                     r"($\sqrt{J_2(\varepsilon_p)} \geq 0{,}1$ do máximo, nos pontos de integração do fator) antes da "
                     r"chuva e no fim dela, linha $p_w = 0$ do transiente e o nível d'água nas três verticais das "
                     r"categorias.", "fig:det"))
        else:
            A(r"\par Sem a figura: %s." % tex_seguro(puladas.get("fig24_det", "")))

    if ref and ref.get("estat"):
        A(r"\section{Configuração de referência (Seção~6.1; análogos das Figs.~25 e 26)}")
        e = ref["estat"]
        A(r"No fim da chuva ($t = %s$~h), %d realizações: $n_f = %d$, $P_f = %s\%%$ (IC de 95\%% de Wilson %s\%%); "
          r"$\mu(FS) = %s$ e $CoV(FS) = %s\%%$ (log-normal ajustada)%s." % (
              g(ref["t"]), e["S"], e["n_f"], pct(e["P_f"]), ic(e, 2), n(e["media"]), pct(e["CoV"], 1),
              ("; $FS$ determinístico %s" % n(det["FS"], 4)) if det.get("FS") else ""))
        A(tabela_ref(ref))
        if fig_rel(figs, "fig25_ref"):
            A(figura(fig_rel(figs, "fig25_ref"), r"(a) Histograma e densidade log-normal ajustada de $FS$ no fim da "
                     r"chuva, com o $FS$ determinístico; (b) convergência de $P_f$ com o número de simulações.",
                     "fig:ref"))
        ts = tabela_spearman(ref)
        if ts:
            A(ts)
        if fig_rel(figs, "fig26_spearman"):
            A(figura(fig_rel(figs, "fig26_spearman"), r"Spearman parcial de 2ª ordem entre $FS$ e $k_s$, $\varphi'$ e "
                     r"$c'$, antes da chuva e no fim dela.", "fig:spearman"))
        else:
            A(r"\par Spearman: %s." % tex_seguro(puladas.get("fig26_spearman", "")))

    if ref and len([i for i in ref.get("instantes", []) if i.get("estat")]) > 1:
        A(r"""\section{Duração da chuva}
A chuva é a mesma até cada instante, então o $FS$ em $t$ é o do fim de uma chuva de $t$ horas com a mesma intensidade.
O mínimo de uma chuva curta pode vir depois do fim dela (a água continua descendo), o que estes instantes não medem.""")
        A(tabela_duracao(ref))
        if fig_rel(figs, "fig_duracao"):
            A(figura(fig_rel(figs, "fig_duracao"), r"(a) $FS$ de cada realização, quantis de 5\% e 95\%, média da "
                     r"log-normal ajustada e $FS$ determinístico contra a duração da chuva; (b) $P_f$ contado (IC de "
                     r"95\% de Wilson) e do ajuste log-normal.", "fig:duracao"))
        pr = ref.get("pares")
        if pr:
            A(r"\section{$FS$ antes e no fim da chuva}")
            A(r"Nas %d realizações com os dois instantes, a correlação entre $FS(0)$ e $FS(%s)$ é %s (Pearson) e %s "
              r"(Spearman); $FS(%s)/FS(0)$ tem média %s e desvio %s." % (
                  len(pr["fs0"]), g(pr["t1"]), n(pr.get("pearson")), n(pr.get("spearman")), g(pr["t1"]),
                  n(pr["razao_media"]), n(pr.get("razao_desvio"))))
            if fig_rel(figs, "fig_fs0_fs24"):
                A(figura(fig_rel(figs, "fig_fs0_fs24"), r"$FS$ antes da chuva contra $FS$ no fim dela, por "
                         r"realização, com a categoria do nível d'água no fim da chuva; linha 1:1 e $FS = 1$.",
                         "fig:fs0", largura=r"0.62\linewidth"))

    if any(c.get("nivel") and sum(c["nivel"]["contagem"].values()) for c in est["casos"]):
        A(r"\section{Nível d'água (Tabela~3 do Rahardjo)}")
        A(r"""As categorias da Tabela~3 de Rahardjo et al.\ (2007) são lidas nas cotas do nível d'água (topo da zona
saturada ligada à base, procurado de 0{,}1 em 0{,}1~m) nas verticais $x = 30$, 35 e 40~m, onde a superfície está a 30,
25 e 20~m; o nível \emph{chega à superfície} se fica a menos de %s~m dela:
\begin{itemize}\small
\item[I:] %s;
\item[II:] %s;
\item[III:] %s;
\item[IV:] %s.
\end{itemize}
A subida é medida desde o nível inicial (%s).""" % (
            n(est["tol_nivel"], 2), TEXTO_CATEGORIA["I"], TEXTO_CATEGORIA["II"], TEXTO_CATEGORIA["III"],
            TEXTO_CATEGORIA["IV"], origem_tex((ref or {}).get("niveis_iniciais_origem"))))
        A(tabela_nivel(est))
        if fig_rel(figs, "fig_nivel"):
            A(figura(fig_rel(figs, "fig_nivel"), r"Frequência das categorias do nível d'água: (a) na referência, com "
                     r"a duração da chuva; (b) em cada caso, no fim da chuva. ``det.'': a categoria do "
                     r"determinístico.", "fig:nivel"))

    if fig_rel(figs, "fig_mapa"):
        A(r"\section{Onde passa a banda de ruptura}")
        A(r"Cada célula de 1~m conta as realizações cuja banda de ruptura passa por ela (\texttt{%s\_mapa}, somado "
          r"nas partes), dividido pelo número de realizações (ou de rupturas) no instante." % tex_seguro(est["prefixo"]))
        A(figura(fig_rel(figs, "fig_mapa"), r"Fração das realizações da referência cuja banda de ruptura "
                 r"($\sqrt{J_2(\varepsilon_p)} \geq 0{,}1$ do máximo) passa por cada célula de 1~m, antes da chuva e "
                 r"no fim dela, e entre as rupturas ($FS < 1$) no fim da chuva, com a banda do determinístico.",
                 "fig:mapa"))

    A(r"\section{Coeficientes de variação (Seção~6.2; análogos das Figs.~27 a 30)}")
    fams_cov = [f for f in ("covk", "covc", "covphi") if figuras_chuva.tem_dados(est, f)]
    if not fams_cov:
        A(r"Não há casos destas famílias nos dados lidos.")
    for fam in fams_cov:
        A(tabela_familia(fam, est))
    letras = "abc"
    paineis = ", ".join("%s (%s)" % (NOME_TEX[f], letras[i]) for i, f in enumerate(fams_cov))
    for nome, txt in (("fig27_cov_pdf", r"Densidades log-normais de $FS$"), ("fig28_cov_pf", r"$P_f$"),
                      ("fig29_cov_mu", r"$\mu(FS)$"), ("fig30_cov_cov", r"$CoV(FS)$")):
        if fig_rel(figs, nome):
            A(figura(fig_rel(figs, nome), txt + " com " + paineis + "." +
                     (r" Anel: a referência; linha pontilhada: o $FS$ determinístico." if nome == "fig29_cov_mu"
                      else ""), "fig:" + nome))
    for fam, nome, titulo in (("s", "fig31_s", "Distâncias de autocorrelação (Seção~6.3; análogo da Fig.~31)"),
                              ("alfa", "fig32_alfa", "Anisotropia da permeabilidade (Seção~6.4; análogo da Fig.~32)"),
                              ("I", "fig33_I", "Intensidade da chuva (no lugar da Seção~6.5; análogo da Fig.~33)")):
        A(r"\section{%s}" % titulo)
        if not figuras_chuva.tem_dados(est, fam):
            A(r"Não há casos desta família nos dados lidos.")
            continue
        if fam == "s":
            A(r"$L_x = 20\,s$~m e $L_y = 2\,s$~m; o determinístico é o da referência.")
        elif fam == "alfa":
            A(r"$k_h = \alpha\,k_v$, com $k_v$ o do campo; cada $\alpha$ tem o seu determinístico.")
        else:
            A(r"Cada intensidade tem o seu determinístico, comparado à sigmoide $F_{min}(I) = 0{,}82 + 1{,}83/(1 + "
              r"e^{-(I - 10{,}7)/(-13{,}13)})$ ajustada por Rahardjo et al.\ (2007, Fig.~7b) aos seus resultados do "
              r"solo f50 a $45^\circ$ (equilíbrio limite, Bishop).")
        A(tabela_familia(fam, est))
        if fig_rel(figs, nome):
            extra = {"s": r" Os valores de $s$ além de 20 aparecem em (b) a (d) como linha horizontal.",
                     "alfa": r" Em (d), o $FS$ determinístico de cada $\alpha$.",
                     "I": r" Em (d), o $FS$ determinístico de cada $I$ e a sigmoide do Rahardjo."}[fam]
            A(figura(fig_rel(figs, nome), r"(a) Densidades log-normais de $FS$; (b) $P_f$; (c) $CoV(FS)$; "
                     r"(d) $\mu(FS)$, com %s. Anel: a referência.%s" % (NOME_TEX[fam], extra), "fig:" + nome))

    A(r"\appendix\section{Realizações e dados}")
    A(r"{\sloppy Dados lidos: \texttt{%s} (prefixo \texttt{%s}). Uma realização entra num instante se terminou "
      r"com \texttt{situacao = ok} e $FS > 0$.\par}" % (caminho_tex(est["raiz"]), tex_seguro(est["prefixo"])))
    excl = [(c["id"], len(c["excluidas"])) for c in est["casos"] if c["excluidas"]]
    if excl:
        A(" Fora (todas as linhas e instantes): " + "; ".join("%s: %d" % (tex_seguro(i), k) for i, k in excl) + ".")
    else:
        A(" Todas as realizações entraram.")
    if est.get("avisos"):
        A(r"\par Avisos da leitura:\begin{itemize}\small\raggedright" +
          "".join(r"\item %s" % tex_seguro(a) for a in est["avisos"][:40]) + r"\end{itemize}")
    if puladas:
        A(r"\par Figuras sem dados:\begin{itemize}\small\raggedright" +
          "".join(r"\item \texttt{%s}: %s" % (tex_seguro(k), tex_seguro(v)) for k, v in puladas.items()) +
          r"\end{itemize}")
    A(r"\end{document}")
    return "\n".join(T)


def main(argv):
    if len(argv) != 3:
        sys.exit(__doc__)
    with open(argv[1]) as f:
        est = json.load(f)
    saida = argv[2]
    os.makedirs(saida, exist_ok=True)
    figs, puladas = figuras_chuva.todas(est, os.path.join(saida, "figs"))
    tex = documento(est, figs, puladas)
    arq = os.path.join(saida, "relatorio.tex")
    with open(arq, "w") as f:
        f.write(tex)
    print(arq)
    for k, v in puladas.items():
        print("figura pulada: %s (%s)" % (k, v))
    if not shutil.which("pdflatex"):
        print("pdflatex nao encontrado: so o .tex")
        return
    for _ in range(2):
        r = subprocess.run(["pdflatex", "-interaction=nonstopmode", "-halt-on-error", "relatorio.tex"], cwd=saida,
                           stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, errors="replace")
        if r.returncode != 0:
            print(r.stdout[-3000:])
            sys.exit("pdflatex falhou")
    print(os.path.join(saida, "relatorio.pdf"))


if __name__ == "__main__":
    main(sys.argv)
