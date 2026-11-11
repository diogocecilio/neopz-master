#!/usr/bin/env python3
"""relatorio_sat.py: o relatorio comparativo (LaTeX e PDF) do exemplo da Secao 6 de Vargas Ceron et al. (2025) com o
material TPZMatUnsaturatedUPwPa no limite saturado: as figuras e tabelas do artigo (extraidas do PDF) e os resultados
deste trabalho nos mesmos graficos, so com os fatos (valores, diferencas e as condicoes de cada lado).

Uso: python3 -I relatorio_sat.py <raiz dos dados> <prefixo (ceron_sat)> <pasta com ceron_dig.json ou -> <saida>
     [rotulo da malha, por exemplo "focada 0,375 m (ceron_h5_a45_foco375)"]"""

import json
import math
import os
import subprocess
import sys

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)
sys.path.insert(0, os.path.join(os.path.dirname(AQUI), "verificacao_ceron"))
import ceron_dig  # noqa: E402
import dados_sat  # noqa: E402
import figuras_sat  # noqa: E402

NOME_TEX = {"covk": r"$CoV(k_v)$", "covc": r"$CoV(c)$", "covphi": r"$CoV(\varphi)$", "s": r"$s$",
            "alfa": r"$\alpha$", "hw": r"$h_w/H$"}
SECAO = {"covk": "6.2", "covc": "6.2", "covphi": "6.2", "s": "6.3", "alfa": "6.4", "hw": "6.5"}


def n(v, casas=3):
    if v is None or (isinstance(v, float) and not math.isfinite(v)):
        return "--"
    return ("%.*f" % (casas, v)).replace(".", "{,}")


def pct(v, casas=2):
    return "--" if v is None or not math.isfinite(v) else n(100. * v, casas)


def dif(a, b):
    """diferenca relativa de b (este trabalho) a a (Ceron), em %"""
    if a is None or b is None or not math.isfinite(a) or not math.isfinite(b) or a == 0:
        return "--"
    return ("%+.1f" % (100. * (b - a) / a)).replace(".", "{,}") + r"\%"


def espiral_da_referencia(est, raiz):
    """o metodo do Ceron (espiral logaritmica) com a pressao do deterministico de referencia deste trabalho"""
    ref = next(c for c in est["casos"] if c["familia"] == "ref")
    det = ref.get("det")
    if not det:
        return None
    pasta = os.path.join(raiz, det["pasta"])
    if not (os.path.isfile(os.path.join(pasta, "pw.csv")) and os.path.isfile(os.path.join(pasta, "triangulos.csv"))):
        return None
    import numpy as np
    import espiral_ceron as E
    pres = E.Pressao(pasta, 0.01)
    (G, v), _ = E.busca(10., math.radians(30.), 20., pres)
    Bx, By = E.ponto_B(min(max(v[2], -1.), 0.99))
    cur = E.curva_mecanismo(v[0], v[1], Bx, By, math.radians(30.))
    if cur is None:
        return None
    return {"G": float(G), "O": [float(v[0]), float(v[1])], "B": [float(Bx), float(By)],
            "x": np.asarray(cur[0]).tolist(), "y": np.asarray(cur[1]).tolist()}


def tabela_familia(fam, est, cer):
    nossos = figuras_sat.casos_familia(est, fam)
    linhas = {r["x"]: r for r in cer["tab"][fam]}
    valores = sorted(set(linhas) | set(nossos))
    tem_det = fam in ("alfa", "hw")
    cab = [NOME_TEX[fam]]
    if tem_det:
        cab += [r"$\Gamma_{det}$ C", r"$\Gamma_{det}$ T"]
    cab += [r"$\mu$ C", r"$\mu$ T", r"$\sigma$ C", r"$\sigma$ T", r"$CoV$ C (\%)", r"$CoV$ T (\%)",
            r"$P_f$ C (\%)", r"$P_f$ T (\%)", r"$CoV(P_f)$ T (\%)", "$S$ T"]
    L = [r"\begin{table}[htbp]\centering\scriptsize\setlength{\tabcolsep}{2.6pt}",
         r"\caption{Seção %s, %s: Ceron et al.\ (C) e este trabalho (T).}" % (SECAO[fam], NOME_TEX[fam]),
         r"\begin{tabular}{" + "r" * len(cab) + "}", r"\toprule", " & ".join(cab) + r" \\", r"\midrule"]
    for v in valores:
        c = linhas.get(v, {})
        t = nossos.get(v, {})
        e = t.get("estat") or {}
        ref = " (ref.)" if v == figuras_sat.VALOR_REF[fam] else ""
        row = [n(v, 2 if fam != "s" else 1).replace("{,}00", "") + ref]
        if tem_det:
            row += [n(c.get("det")), n((t.get("det") or {}).get("Gamma"))]
        row += [n(c.get("mu")), n(e.get("media")), n(c.get("sigma")), n(e.get("desvio")),
                pct(c.get("cov"), 1), pct(e.get("CoV"), 1), pct(c.get("pf")), pct(e.get("P_f")),
                pct(e.get("CoV_Pf"), 1), str(e.get("S", "--"))]
        L.append(" & ".join(row) + r" \\")
    L += [r"\bottomrule", r"\end{tabular}", r"\end{table}"]
    return "\n".join(L)


def figura(arq, legenda, rotulo, notas=None, largura=r"\linewidth"):
    extra = ""
    if notas:
        extra = " Eixos ampliados em relação aos do artigo, para conter os dados: " + "; ".join(notas) + "."
    return "\n".join([r"\begin{figure}[htbp]\centering", r"\includegraphics[width=%s]{%s}" % (largura, arq),
                      r"\caption{%s%s}\label{%s}" % (legenda, extra, rotulo), r"\end{figure}"])


SEM_DADOS = (r"Não há resultados deste trabalho para estes casos no conjunto de dados lido (ver os avisos no "
             r"apêndice).")


def tex_seguro(t):
    """texto qualquer como texto do LaTeX (um passe so, para nao reescapar o que ja foi escapado)"""
    mapa = {"\\": r"\textbackslash{}", "_": r"\_", "%": r"\%", "&": r"\&", "#": r"\#", "$": r"\$", "{": r"\{",
            "}": r"\}", "~": r"\textasciitilde{}", "^": r"\textasciicircum{}"}
    return "".join(mapa.get(ch, ch) for ch in t)


def documento(est, cer, esp, figs, notas, rotulo_malha):
    ref = next(c for c in est["casos"] if c["familia"] == "ref")
    e = ref["estat"]
    det = ref.get("det") or {}
    versoes = ", ".join("%s (%d)" % (k.replace("_", r"\_"), v) for k, v in sorted(est["versoes"].items()))
    excl = [(c["id"], len(c["excluidas"])) for c in est["casos"] if c["excluidas"]]
    T = []
    A = T.append
    A(r"""\documentclass[11pt,a4paper]{article}
\usepackage[T1]{fontenc}\usepackage[utf8]{inputenc}\usepackage[brazil]{babel}
\usepackage[margin=2.2cm]{geometry}\usepackage{booktabs}\usepackage{graphicx}\usepackage{caption}
\usepackage{float}\usepackage[hidelinks]{hyperref}\usepackage{amsmath}
\captionsetup{font=small}
\author{SlopeRandomField, \texttt{comparacao\_ceron}}
\title{Exemplo da Seção 6 de Vargas Ceron et al.\ (2025) com o material \texttt{TPZMatUnsaturatedUPwPa} no limite
saturado: comparação}
\date{}
\begin{document}\maketitle
""")
    A(r"""\section{O que foi rodado}
O exemplo é o da Seção 6 do artigo (talude de $H = 5$~m e $45^\circ$; rebaixamento instantâneo do reservatório da crista
até o pé, $h_w = H$; $c'$, $\varphi'$ e $k_v$ log-normais com as médias, os coeficientes de variação e as distâncias de
autocorrelação da Tabela 2 do artigo). Os casos das Seções 6.1 a 6.5 rodaram em três fases no fermi, com 1000
realizações por caso.

\begin{table}[htbp]\centering\footnotesize
\caption{As duas análises.}
\begin{tabular}{p{0.2\linewidth}p{0.36\linewidth}p{0.36\linewidth}}\toprule
& Ceron et al.\ (2025) & este trabalho \\\midrule
água & percolação permanente saturada, acoplamento de uma via (Eqs.~19 a 21) & \texttt{TPZMatUnsaturatedUPwPa}, variante
iii, com $S_w = 1$, $dS_w/dp_w = 0$, $k_{rw} = 1$ (\texttt{SATURADO=1}), água incompressível (\texttt{KW=1e30}) e
$p_w = 0$ na superfície acima do nível (\texttt{SUPERFICIE\_P0=1}): a mesma percolação permanente \\
domínio da água & 20~m $\times$ 10~m nas realizações (Figs.~12 e 13), bordas impermeáveis & 25~m $\times$ 10~m:
esquerda com a carga inicial ($u = 0$), direita com o reservatório, base impermeável \\
fator & $\Gamma$ pela análise limite cinemática, mecanismos em espiral logarítmica por partes & $\Gamma$ pelo aumento da
gravidade (GIM) com \emph{arc-length}, elementos finitos (deslocamentos quadráticos), em $t = 0{,}01$~h \\
malha & -- & """ + rotulo_malha + r""" \\
simulações & $S = 10\,000$ na referência (até 30\,000 em outros casos) & 1000 por caso \\
campos & Karhunen--Loève, $M = 2000$ & Karhunen--Loève, $M = 2000$ (variância retida média """ +
      pct(ref.get("KL_variancia_media"), 1) + r"""\%) \\
\bottomrule\end{tabular}\end{table}

Versões do executável nas realizações: """ + versoes + r""".
""")
    A(r"""\section{Determinístico (Fig.~24)}
Com as médias, o GIM deste trabalho dá $\Gamma = """ + n(det.get("Gamma")) + r"""$; o artigo dá $\Gamma = 1{,}336$.""")
    if esp:
        A(r""" O método do próprio Ceron (\texttt{verificacao\_ceron/espiral\_ceron.py}, mecanismos rotacionais em espiral
logarítmica, Eqs.~42 a 56), com a pressão de água deste trabalho, dá $\Gamma = """ + n(esp["G"]) + r"""$. A mesma
implementação dá $\Gamma = 1{,}7770$ no exemplo seco da Seção 5.3.2 do artigo, que publica $1{,}777$.""")
    A("\n" + figura("figs/fig24.pdf", r"Talude, superfície crítica do artigo (extraída da Fig.~24), espiral crítica do "
                    r"método do Ceron com a pressão deste trabalho e a banda de ruptura do GIM no determinístico.", "fig:24"))
    A(r"""
\section{Configuração de referência (Seção 6.1, Figs.~25 e 26)}
\begin{table}[htbp]\centering\footnotesize
\caption{Referência: Ceron et al.\ (C) e este trabalho (T). $\mu$, mediana, $\sigma$ e $CoV$ são os da log-normal
ajustada, como no artigo; entre parênteses, os amostrais.}
\begin{tabular}{lrr}\toprule & C & T \\\midrule
$S$ & 10\,000 & """ + str(e["S"]) + r""" \\
$n_f = \#(\Gamma < 1)$ & 1150 & """ + str(e["n_f"]) + r""" \\
$P_f$ (\%) & 11{,}50 & """ + pct(e["P_f"]) + r""" \\
$CoV(P_f)$ (\%) & 2{,}77 & """ + pct(e["CoV_Pf"]) + r""" \\
$\mu$ & 1{,}353 & """ + n(e["media"]) + " (" + n(e["media_amostral"]) + r""") \\
mediana & 1{,}313 & """ + n(e["mediana"]) + r""" \\
$\sigma$ & 0{,}318 & """ + n(e["desvio"]) + " (" + n(e["desvio_amostral"]) + r""") \\
$CoV$ (\%) & 23{,}5 & """ + pct(e["CoV"], 1) + " (" + pct(e["CoV_amostral"], 1) + r""") \\
qui-quadrado, $p$ da normal & -- & """ + n(e["qui2_normal"]["p"], 3) + r""" \\
qui-quadrado, $p$ da log-normal & -- & """ + n(e["qui2_lognormal"]["p"], 3) + r""" \\
\bottomrule\end{tabular}\end{table}
""")
    A(figura("figs/fig25.pdf", r"(a) Histograma e densidade log-normal ajustada de $\Gamma$; (b) convergência de $P_f$ "
             r"com o número de simulações.", "fig:25", notas.get("25")))
    sb, sc = ref.get("spearman_banda"), ref.get("spearman_campo")
    if sb and sc:
        linhas = []
        for k, nome in (("k", r"$k_v$"), ("phi", r"$\varphi$"), ("c", r"$c$")):
            cv = cer["fig26"].get(k) if cer.get("fig26") else None
            linhas.append(" & ".join([nome, n(cv, 2), n(sb[k]["rho_parcial"], 3) + " [" + n(sb[k]["ic95"][0], 2) + "; " +
                                      n(sb[k]["ic95"][1], 2) + "]", n(sc[k]["rho_parcial"], 3) + " [" +
                                      n(sc[k]["ic95"][0], 2) + "; " + n(sc[k]["ic95"][1], 2) + "]"]) + r" \\")
        A(r"""\begin{table}[htbp]\centering\footnotesize
\caption{Spearman parcial de 2ª ordem entre $\Gamma$ e cada variável, controlando as outras duas (Eqs.~70 e 71), com
o intervalo de 95\% (transformação de Fisher). Este trabalho, duas definições das variáveis de cada realização: $c'$ e
$\varphi'$ médios na banda de ruptura e $k$ médio na massa que desliza; e as médias de cada campo no solo todo.}
\begin{tabular}{lrrr}\toprule
& Ceron (Fig.~26) & banda e massa & médias do campo \\\midrule
""" + "\n".join(linhas) + r"""
\bottomrule\end{tabular}\end{table}
""")
    A(figura("figs/fig26.pdf", r"Spearman parcial de 2ª ordem (Fig.~26 do artigo e este trabalho).", "fig:26",
             largura=r"0.7\linewidth"))
    A(r"\section{Coeficientes de variação (Seção 6.2, Figs.~27 a 30)}")
    if "27" not in figs:
        A(SEM_DADOS)
    for fam in ("covk", "covc", "covphi"):
        if figuras_sat.tem_dados(est, fam):
            A(tabela_familia(fam, est, cer))
    if "27" in figs:
        A(figura("figs/fig27.pdf", r"Densidades log-normais de $\Gamma$ com o $CoV$ de $k_v$ (a), $c$ (b) e $\varphi$ (c). "
             r"Tracejado: Ceron et al.; cheio: este trabalho.", "fig:27", notas.get("27")))
    for nf, txt in ((28, r"$P_f$"), (29, r"$\mu(\Gamma)$"), (30, r"$CoV(\Gamma)$")):
        if str(nf) not in figs:
            continue
        A(figura("figs/fig%d.pdf" % nf, txt + r" contra o $CoV$ de $k_v$ (a), $c$ (b) e $\varphi$ (c). O valor fora da "
                 r"faixa do eixo x do artigo ($CoV(k_v) = 0$) aparece como linha horizontal, como no artigo.",
                 "fig:%d" % nf, notas.get(str(nf))))
    for fam, nf, sec, txt in (("s", 31, "6.3", "distâncias de autocorrelação"),
                              ("alfa", 32, "6.4", "anisotropia da permeabilidade"),
                              ("hw", 33, "6.5", "altura do rebaixamento")):
        A(r"\section{%s%s (Seção %s, Fig.~%d)}" % (txt[0].upper(), txt[1:], sec, nf))
        if str(nf) not in figs:
            A(SEM_DADOS)
            continue
        A(tabela_familia(fam, est, cer))
        A(figura("figs/fig%d.pdf" % nf, r"(a) Densidades log-normais de $\Gamma$; (b) $P_f$; (c) $CoV(\Gamma)$; "
                 r"(d) $\mu(\Gamma)$, com %s." % NOME_TEX[fam] +
                 (r" O caso $s = 400$ aparece em (b) a (d) como linha horizontal, como no artigo." if fam == "s" else ""),
                 "fig:%d" % nf, notas.get(str(nf))))
    A(r"\appendix\section{Realizações e dados}")
    A(r"Uma realização entra se terminou com \texttt{situacao = ok} e $\Gamma > 0$. ")
    if excl:
        A("Fora: " + "; ".join(r"%s: %d" % (i.replace("_", r"\_"), k) for i, k in excl) + ".")
    else:
        A("Todas as realizações entraram.")
    sem = list(dict.fromkeys(a.split(":")[0] for a in est.get("avisos", []) if a.split(":", 1)[-1].startswith(" sem ")))
    if sem:
        A(r"\par Sem resultados no conjunto lido: " + ", ".join(tex_seguro(x) for x in sem) + ".")
    outros = [a for a in est.get("avisos", []) if a.split(":")[0] not in sem]
    if outros:
        A(r"\par Outros avisos da leitura: " + "; ".join(tex_seguro(a) for a in outros) + ".")
    for fase in (1, 2, 3):
        arq = os.path.join(est["raiz"], "%s_fase%d" % (est["prefixo"], fase), "AVISOS")
        if os.path.isfile(arq):
            linhas = [l.strip() for l in open(arq, errors="replace") if l.strip()]
            A(r"\par Avisos da submissão da fase %d (\texttt{AVISOS}):" % fase)
            A(r"\begin{itemize}\small\raggedright" + "".join(r"\item \texttt{%s}" % tex_seguro(l) for l in linhas[:10]) +
              r"\end{itemize}")
    A(r"""\par As curvas do artigo foram extraídas dos objetos vetoriais do PDF (ver \texttt{ceron\_dig.json} e os
arquivos de extração), e os valores escalares vêm das Tabelas 3 a 6. Onde uma densidade do artigo não foi extraída, a
figura usa a log-normal com o $\mu$ e o $\sigma$ da tabela (o artigo desenha as log-normais ajustadas).
\end{document}""")
    return "\n".join(T)


def main():
    if len(sys.argv) < 5:
        sys.exit(__doc__)
    raiz, prefixo, digit, saida = sys.argv[1:5]
    rot_malha = sys.argv[5] if len(sys.argv) > 5 else "--"
    os.makedirs(saida, exist_ok=True)
    est = dados_sat.calcula(raiz, prefixo)
    with open(os.path.join(saida, "estatisticas.json"), "w") as f:
        json.dump(dados_sat._limpa(est), f)
    cer = ceron_dig.carrega(None if digit == "-" else digit)
    esp = espiral_da_referencia(est, raiz)
    figs, notas = figuras_sat.todas(est, cer, esp, os.path.join(saida, "figs"))
    tex = documento(est, cer, esp, figs, notas, rot_malha)
    with open(os.path.join(saida, "relatorio.tex"), "w") as f:
        f.write(tex)
    for _ in range(2):
        r = subprocess.run(["pdflatex", "-interaction=nonstopmode", "-halt-on-error", "relatorio.tex"], cwd=saida,
                           stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, errors="replace")
        if r.returncode != 0:
            print(r.stdout[-3000:])
            sys.exit("pdflatex falhou")
    print(os.path.join(saida, "relatorio.pdf"))


if __name__ == "__main__":
    main()
