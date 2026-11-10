#!/usr/bin/env python3
"""tabelas.py: as tabelas do relatorio em LaTeX (booktabs), com as colunas das Tabelas 2 a 6 de Vargas Ceron et al.
(2025) e uma ultima coluna S (numero de realizacoes), mais as tabelas do apendice e a de proveniencia.

gera(est, saida) le o dicionario de dados.calcula (ou de <saida>/estatisticas.json) e grava em <saida>/tabelas/:
<solo>_tab2.tex (dados de referencia, Tabela 2), <solo>_tab3.tex (CoV, Tabela 3), <solo>_tab4.tex (distancias de
autocorrelacao, Tabela 4), <solo>_tab5.tex (anisotropia, Tabela 5), <solo>_tab6.tex (intensidade da chuva, no lugar do
rebaixamento da Tabela 6), <solo>_tabA1.tex (estatisticas amostrais), <solo>_tabA2.tex (qui-quadrado),
<solo>_tabA3.tex (Spearman) e tab_proveniencia.tex. Cada arquivo e um ambiente table completo (rotulo tab:<solo>_tabN,
tab:proveniencia). Formatos do Ceron: mu e sigma com 3 casas, CoV com 1, P_f e CoV(P_f) com 2, o FS deterministico com
os seus digitos; virgula decimal."""

import json
import math
import os
import sys

_AQUI = os.path.dirname(os.path.abspath(__file__))
_PAI = os.path.dirname(_AQUI)
for _p in (_PAI, _AQUI):
    if _p not in sys.path:
        sys.path.insert(0, _p)

SUFIXOS = ["covk0", "covk0.75", "covk0.9", "covk1", "covc0.1", "covc0.5", "covc0.7", "covphi0.05", "covphi0.15",
           "covphi0.2", "lx40_ly4", "lx100_ly10", "lx8000_ly800", "alfa2", "alfa5", "alfa10"]


# ------------------------------------------------------------------------------------------------------- formatos
def tex(s):
    """texto -> LaTeX (escapa os caracteres especiais dos nomes de pastas)"""
    s = str(s)
    for a, b in (("\\", r"\textbackslash{}"), ("_", r"\_"), ("%", r"\%"), ("&", r"\&"), ("#", r"\#"),
                 ("$", r"\$"), ("{", r"\{"), ("}", r"\}")):
        s = s.replace(a, b)
    return s


def _ok(v):
    return v is not None and isinstance(v, (int, float)) and math.isfinite(v)


def num(v, casas):
    """numero com virgula decimal e casas fixas; None -> '--'"""
    if not _ok(v):
        return "--"
    s = f"{v:.{casas}f}"
    if s.startswith("-"):
        s = s[1:]
        if s.strip("0.") == "":
            return s.replace(".", ",")
        return "$-$" + s.replace(".", ",")
    return s.replace(".", ",")


def num_g(v):
    """numero com os seus digitos (ate 7 casas, sem zeros a direita), virgula decimal"""
    if not _ok(v):
        return "--"
    s = f"{v:.7f}".rstrip("0").rstrip(".")
    if s.startswith("-"):
        return "$-$" + s[1:].replace(".", ",")
    return s.replace(".", ",")


def nm(v):
    """numero para o modo matematico: os seus digitos (ate 7 casas), virgula decimal sem espaco ({,})"""
    if not _ok(v):
        return r"\textrm{--}"
    return f"{v:.7f}".rstrip("0").rstrip(".").replace(".", "{,}")


def txt_pt(s):
    """o texto de um numero como esta no arquivo (os digitos exatos), com virgula decimal"""
    s = str(s).strip()
    if s.startswith("-"):
        return "$-$" + s[1:].replace(".", ",")
    return s.replace(".", ",")


def pct(v, casas):
    return num(100.0 * v, casas) if _ok(v) else "--"


def num_p(p):
    """valor p: 3 casas a partir de 0,001; abaixo, notacao cientifica com 2 algarismos"""
    if not _ok(p):
        return "--"
    if p >= 0.001:
        return num(p, 3)
    if p == 0:
        return "0"
    e = int(math.floor(math.log10(p)))
    m = p / 10 ** e
    if round(m, 1) >= 10:
        m, e = m / 10, e + 1
    return f"${m:.1f}".replace(".", "{,}") + rf"\times 10^{{{e}}}$"


def negrito(s, sim):
    return rf"\textbf{{{s}}}" if sim and s not in ("", "--") else s


def _versao(v):
    """'a62653a_2026-11-06_08:56:16' -> 'a62653a (2026-11-06)'"""
    partes = str(v).split("_")
    return f"{partes[0]} ({partes[1]})" if len(partes) >= 2 else str(v)


def _lista_pt(itens):
    itens = list(itens)
    if not itens:
        return ""
    return itens[0] if len(itens) == 1 else ", ".join(itens[:-1]) + " e " + itens[-1]


NOME_FAMILIA = {"ref": "de referência", "covk": "de $CoV$", "covc": "de $CoV$", "covphi": "de $CoV$",
                "s": "de $s$", "alfa": r"de $\alpha$", "chuva": "de intensidade da chuva $I$"}


def rotulo(c):
    """rotulo LaTeX de um caso"""
    f, v = c.get("familia"), c.get("valor")
    if f == "ref":
        return "Referência"
    if f == "covk":
        return rf"$CoV(k_v) = {nm(100 * v)}\%$"
    if f == "covc":
        return rf"$CoV(c') = {nm(100 * v)}\%$"
    if f == "covphi":
        return rf"$CoV(\varphi') = {nm(100 * v)}\%$"
    if f == "s":
        return f"$s = {nm(v)}$"
    if f == "alfa":
        return rf"$\alpha = {nm(v)}$"
    return f"$I = {nm(v)}$ mm/h"


def _ambiente(rotulo_tab, legenda, corpo, notas=(), fonte=r"\small", colsep=None, fonte_notas=r"\footnotesize"):
    """um ambiente table completo: legenda acima (como no Ceron), tabular e notas abaixo. O tabular vai numa caixa;
    se ela for mais larga que a linha (e graphicx estiver carregado), e reduzida a largura da linha"""
    l = [r"\begin{table}[htbp]", r"\centering", rf"\caption{{{legenda}}}", rf"\label{{{rotulo_tab}}}", fonte]
    if colsep:
        l.append(rf"\setlength{{\tabcolsep}}{{{colsep}}}")
    l += [r"\ifdefined\caixatabelaceron\else\newsavebox{\caixatabelaceron}\fi",
          r"\sbox{\caixatabelaceron}{%"]
    l += corpo
    l[-1] += "}%"
    l += [r"\ifdefined\resizebox\ifdim\wd\caixatabelaceron>\linewidth"
          r"\resizebox{\linewidth}{!}{\usebox{\caixatabelaceron}}\else\usebox{\caixatabelaceron}\fi"
          r"\else\usebox{\caixatabelaceron}\fi"]
    notas = [n for n in notas if n]
    if notas:
        l += [r"\par\smallskip", r"\begin{minipage}{\linewidth}", fonte_notas + r"\raggedright"]
        l += [n + (r"\par" if i < len(notas) - 1 else "") for i, n in enumerate(notas)]
        l += [r"\end{minipage}"]
    l += [r"\end{table}", ""]
    return "\n".join(l)


def _grava(saida, nome, texto, gravados):
    d = os.path.join(saida, "tabelas")
    os.makedirs(d, exist_ok=True)
    arq = os.path.join(d, nome)
    with open(arq, "w", encoding="utf-8") as f:
        f.write(texto)
    gravados.append(arq)


# ----------------------------------------------------------------------------------------------------- auxiliares
def _casos_ordem(S):
    """referencia, CoV (k_v, c', phi'), s, alfa e chuva, sem repetir o caso de referencia nem dados iguais aos dele"""
    casos, fam = S["casos"], S["familias"]
    ids = [S["ref"]]
    for f in ("covk", "covc", "covphi", "s", "alfa", "chuva"):
        for i in fam.get(f, []):
            if i in ids or i not in casos:
                continue
            if i != S["ref"] and casos[i].get("pastas") == casos[S["ref"]].get("pastas"):
                continue  # os mesmos dados da referencia (f100, I = 36 mm/h)
            ids.append(i)
    return ids


def _duplicados_ref(S):
    casos = S["casos"]
    return [i for i, c in casos.items() if i != S["ref"] and c.get("pastas") == casos[S["ref"]].get("pastas")]


def _nota_estat(c):
    return (r"$\mu$, $\sigma$ e $CoV$: média, desvio-padrão e coeficiente de variação da distribuição log-normal "
            r"ajustada ao $FS$ por máxima verossimilhança. $P_f = n_f/S$, com $n_f$ o número de realizações com "
            r"$FS < 1$; $CoV(P_f) = \sqrt{(1 - P_f)/(S\,P_f)}$; $S$: número de realizações. " + _nota_fs(c))


def _nota_fs(c):
    if c.get("fs_def") == "FS_24h":
        return r"$FS$: fator de segurança do SRM em $t = 24$ h."
    inst = _lista_pt([num_g(t) for t in c.get("instantes_h", [])])
    return rf"$FS$: o menor fator de segurança do SRM entre os instantes {inst} h de cada realização."


def _nota_pf0(casos_ids, casos):
    if any(casos[i].get("Pf") == 0 for i in casos_ids):
        return r"--: $P_f = 0$, $CoV(P_f)$ indefinido."
    return ""


def _nota_negativo(casos_ids, casos):
    itens = [f"{rotulo(casos[i])}: {casos[i]['n_fs_negativo']}" for i in casos_ids
             if casos[i].get("n_fs_negativo", 0) > 0]
    if not itens:
        return ""
    return ("Realizações com $FS = -1$ (sem equilíbrio nem com $c'/0{,}6$ e $\\tan\\varphi'/0{,}6$), contadas em "
            "$n_f$ e fora do ajuste: " + "; ".join(itens) + ".")


def _nota_det_versao(casos_ids, casos):
    itens = []
    for i in casos_ids:
        c = casos[i]
        d = c.get("det")
        if not d:
            continue
        vmc = sorted(c.get("versoes", {}))
        if d.get("versao") not in vmc:
            itens.append(f"{rotulo(c)}: rodada determinística {tex(d['pasta'])}, versão {_versao(d['versao'])}; "
                         f"Monte Carlo, versão " + ", ".join(_versao(v) for v in vmc))
    return ("Versões do programa: " + "; ".join(itens) + ".") if itens else ""


def _linha_estat(c):
    """mu, sigma, CoV (%) do ajuste log-normal, P_f (%), CoV(P_f) (%), S"""
    ln = c.get("lognormal") or {}
    return [num(ln.get("media"), 3), num(ln.get("dp"), 3), pct(ln.get("cov"), 1), pct(c.get("Pf"), 2),
            pct(c.get("cov_Pf"), 2), str(c.get("N", "--"))]


def _det_txt(c):
    d = c.get("det")
    if not d:
        return "--"
    return txt_pt(d.get("FS_txt", d.get("FS")))


# ---------------------------------------------------------------------------------------------------------- tabelas
def tab2(solo, S, prov):
    c = S["casos"][S["ref"]]
    k_s = S.get("k_s_m_s")
    k_c = S.get("k_s_ceron", k_s / 9.81 if _ok(k_s) else None)
    const = prov.get("constantes", {})
    gw = const.get("gamma_w_kN_m3", 9.81)
    cmed, phimed = const.get("c_medio_kPa", 10.0), const.get("phi_medio_graus", 26.0)
    gamma = const.get("gamma_kN_m3", 20.0)

    def pot10(v):
        e = int(math.floor(math.log10(v)))
        m = v / 10 ** e
        if abs(m - 1) < 1e-12:
            return rf"10^{{{e}}}"
        return (f"{m:.2f}".replace(".", "{,}") + rf"\times 10^{{{e}}}")

    L = [num_g(c["Lx"]), num_g(c["Ly"])]
    corpo = [
        r"\begin{tabular}{lp{3.9cm}cccc}",
        r"\toprule",
        r" & & & Coeficiente & \multicolumn{2}{c}{Distância de} \\",
        r" & & & de & \multicolumn{2}{c}{autocorrelação} \\",
        r"\cmidrule(lr){5-6}",
        r"\multicolumn{2}{l}{Parâmetros} & Valor médio & variação & $L_x$ (m) & $L_y$ (m) \\",
        r"\midrule",
        rf"$c'$ & Coesão efetiva & {num_g(cmed)} kPa & {num_g(100 * c['cov_c'])}\% & {L[0]} & {L[1]} \\",
        rf"$\varphi'$ & Ângulo de atrito efetivo & {num_g(phimed)}$^\circ$ & {num_g(100 * c['cov_phi'])}\% & "
        rf"{L[0]} & {L[1]} \\",
        rf"$k_v$ & Permeabilidade vertical & ${pot10(k_s)}$ m/s & {num_g(100 * c['cov_k'])}\% & "
        rf"{L[0]} & {L[1]} \\",
        rf" & \quad na unidade do Ceron & ${pot10(k_c)}$ m$^4$/kN/s & & & \\",
        rf"$\alpha$ & Razão de anisotropia da permeabilidade & {num_g(c['alfa'])} & -- & -- & -- \\",
        rf"$\gamma$ & Peso específico do solo & {num_g(gamma)} kN/m$^3$ & -- & -- & -- \\",
        r"\bottomrule",
        r"\end{tabular}",
    ]
    notas = [
        rf"$k_v$: permeabilidade saturada $k_s$ do solo, a média do campo aleatório de $k_v$ ($k_v/k_s$ "
        rf"log-normal de média 1). Na unidade do Ceron (m$^2$ kPa$^{{-1}}$ "
        rf"s$^{{-1}}$ = m$^4$ kN$^{{-1}}$ s$^{{-1}}$): $k$ (m/s)$/\gamma_w$, "
        rf"$\gamma_w = {nm(gw)}$ kN/m$^3$. Demais dados, sem variabilidade espacial: porosidade "
        rf"{num_g(const.get('porosidade', 0.45))}, $E = {nm(const.get('E_MPa', 20))}$ MPa, "
        rf"$\nu = {nm(const.get('nu', 0.3))}$.",
        rf"Campos aleatórios: expansão de Karhunen-Loève com $M = {c['KL_M']}$ termos (variância retida "
        rf"{pct(c['KL_variancia'], 2)}\%), covariância exponencial, distribuição log-normal, os três campos "
        rf"independentes. Chuva de referência: $I = {nm(c['I_mm_h'])}$ mm/h durante "
        rf"{num_g(const.get('chuva_duracao_h', 24))} h. {_nota_fs(c)} $S = {c['N']}$ realizações.",
    ]
    return _ambiente(f"tab:{solo}_tab2", f"Dados estatísticos do modelo de referência do solo {solo}.", corpo, notas)


def tab3(solo, S):
    casos, ref = S["casos"], S["ref"]
    fam = S["familias"]
    ids = [ref] + [i for f in ("covk", "covc", "covphi") for i in fam.get(f, []) if i != ref and i in casos]
    linhas = []
    for i in ids:
        c = casos[i]
        r = i == ref
        cv = [num_g(100 * c["cov_k"]), num_g(100 * c["cov_c"]), num_g(100 * c["cov_phi"])]
        var = {"covk": 0, "covc": 1, "covphi": 2}.get(c["familia"])
        cv = [negrito(v, r or j == var) for j, v in enumerate(cv)]
        est = [negrito(v, r) for v in _linha_estat(c)]
        linhas.append(" & ".join(cv + est) + r" \\")
    corpo = [
        r"\begin{tabular}{ccccccccc}",
        r"\toprule",
        r"\multicolumn{3}{c}{Coeficiente de variação (\%)} & \multicolumn{3}{c}{Fator de segurança $FS$} & & & \\",
        r"\cmidrule(lr){1-3}\cmidrule(lr){4-6}",
        r"$k_v$ & $c'$ & $\varphi'$ & $\mu$ & $\sigma$ & $CoV$ (\%) & $P_f$ (\%) & $CoV(P_f)$ (\%) & $S$ \\",
        r"\midrule",
    ] + linhas + [r"\bottomrule", r"\end{tabular}"]
    notas = [_nota_estat(casos[ref]), _nota_pf0(ids, casos), _nota_negativo(ids, casos)]
    return _ambiente(f"tab:{solo}_tab3", f"Resumo das simulações de Monte Carlo do solo {solo}: efeito do "
                                         "coeficiente de variação.", corpo, notas)


def tab4(solo, S):
    casos, ref = S["casos"], S["ref"]
    ids = [i for i in S["familias"].get("s", []) if i in casos]
    linhas = []
    for i in ids:
        c = casos[i]
        r = i == ref
        linhas.append(" & ".join([negrito(num(c["s"], 1), r)] + [negrito(v, r) for v in _linha_estat(c)]) + r" \\")
    c0 = casos[ref]
    corpo = [
        r"\begin{tabular}{ccccccc}",
        r"\toprule",
        r"Fator de escala & \multicolumn{3}{c}{Fator de segurança $FS$} & & & \\",
        r"\cmidrule(lr){1-1}\cmidrule(lr){2-4}",
        r"$s$ & $\mu$ & $\sigma$ & $CoV$ (\%) & $P_f$ (\%) & $CoV(P_f)$ (\%) & $S$ \\",
        r"\midrule",
    ] + linhas + [r"\bottomrule", r"\end{tabular}"]
    notas = [rf"$(L_x, L_y) = s\,({nm(c0['Lx'])}\ \mathrm{{m}}, {nm(c0['Ly'])}\ \mathrm{{m}})$.",
             _nota_estat(c0), _nota_pf0(ids, casos), _nota_negativo(ids, casos)]
    return _ambiente(f"tab:{solo}_tab4", f"Resumo das simulações de Monte Carlo do solo {solo}: efeito das "
                                         "distâncias de autocorrelação.", corpo, notas)


def tab5(solo, S):
    casos, ref = S["casos"], S["ref"]
    ids = [i for i in S["familias"].get("alfa", []) if i in casos]
    linhas = []
    for i in ids:
        c = casos[i]
        r = i == ref
        linhas.append(" & ".join([negrito(num(c["alfa"], 1), r), negrito(_det_txt(c), r)]
                                 + [negrito(v, r) for v in _linha_estat(c)]) + r" \\")
    corpo = [
        r"\begin{tabular}{cccccccc}",
        r"\toprule",
        r" & \multicolumn{4}{c}{Fator de segurança $FS$} & & & \\",
        r"\cmidrule(lr){2-5}",
        r" & Análise & & & & & & \\",
        r"$\alpha$ & determinística & $\mu$ & $\sigma$ & $CoV$ (\%) & $P_f$ (\%) & $CoV(P_f)$ (\%) & $S$ \\",
        r"\midrule",
    ] + linhas + [r"\bottomrule", r"\end{tabular}"]
    sem_det = [rotulo(casos[i]) for i in ids if not casos[i].get("det")]
    notas = [r"$\alpha = k_h/k_v$. Análise determinística: o $FS$ da rodada com os valores médios, sem campos "
             r"aleatórios, pela mesma definição de $FS$.",
             (f"--: sem rodada determinística ({_lista_pt(sem_det)})." if sem_det else ""),
             _nota_det_versao(ids, casos),
             _nota_estat(casos[ref]), _nota_pf0(ids, casos), _nota_negativo(ids, casos)]
    return _ambiente(f"tab:{solo}_tab5", f"Resumo das simulações de Monte Carlo do solo {solo}: efeito da razão "
                                         "de anisotropia da permeabilidade.", corpo, notas)


def _nota_tmin(c):
    cont = c.get("t_min_h_contagem") or []
    return ", ".join(f"{num_g(t)} h ({n})" for t, n in cont)


def tab6(solo, S, prov):
    casos, ref, refc = S["casos"], S["ref"], S.get("ref_chuva")
    ids = [i for i in S["familias"].get("chuva", []) if i in casos]
    linhas = []
    for i in ids:
        c = casos[i]
        r = i == refc
        linhas.append(" & ".join([negrito(num_g(c["I_mm_h"]), r), negrito(_det_txt(c), r)]
                                 + [negrito(v, r) for v in _linha_estat(c)]) + r" \\")
    corpo = [
        r"\begin{tabular}{cccccccc}",
        r"\toprule",
        r" & \multicolumn{4}{c}{Fator de segurança $FS$} & & & \\",
        r"\cmidrule(lr){2-5}",
        r" & Análise & & & & & & \\",
        r"$I$ (mm/h) & determinística & $\mu$ & $\sigma$ & $CoV$ (\%) & $P_f$ (\%) & $CoV(P_f)$ (\%) & $S$ \\",
        r"\midrule",
    ] + linhas + [r"\bottomrule", r"\end{tabular}"]
    notas = [rf"Chuva de intensidade $I$ durante 24 h. Em negrito, a intensidade de referência ($I = "
             rf"{nm(casos[refc]['I_mm_h'])}$ mm/h)." if refc in casos else ""]
    if refc in casos and ref in casos:
        cr, c0 = casos[refc], casos[ref]
        if S.get("mesmos_dados_ref_chuva") or cr.get("pastas") == c0.get("pastas"):
            notas.append(rf"$I = {nm(cr['I_mm_h'])}$ mm/h: os mesmos dados do caso de referência das Tabelas "
                         rf"\ref{{tab:{solo}_tab2}} a \ref{{tab:{solo}_tab5}} ($S = {cr['N']}$).")
        else:
            camp_r = ", ".join(tex(p.split("/")[0]) for p in cr["pastas"])
            camp_0 = ", ".join(tex(p.split("/")[0]) for p in c0["pastas"])
            txt = (rf"A linha $I = {nm(cr['I_mm_h'])}$ mm/h (campanha {camp_r}, $S = {cr['N']}$, "
                   rf"{_nota_fs(cr)[:-1].replace('$FS$: o menor', 'o menor')}) não é o mesmo conjunto de dados do "
                   rf"caso de referência das Tabelas \ref{{tab:{solo}_tab2}} a \ref{{tab:{solo}_tab5}} (campanha "
                   rf"{camp_0}, $S = {c0['N']}$, {_nota_fs(c0)[:-1].replace('$FS$: fator', 'fator')}).")
            comp = (prov or {}).get("comparacao_f50_ref")
            if comp and set(comp.get("pastas", [])) == set(c0["pastas"] + cr["pastas"]):
                txt += (rf" Nas realizações {comp['faixa_comum'][0]} a {comp['faixa_comum'][1]}, comuns às duas "
                        rf"campanhas, o $FS$ em $t = 24$ h é idêntico em {comp['iguais']['FS_24h']} de "
                        rf"{comp['reais_comuns']}.")
            notas.append(txt)
    notas.append(_nota_fs(casos[ids[0]]) if ids else "")
    tmins = [f"{rotulo(casos[i])}: {_nota_tmin(casos[i])}" for i in ids if casos[i].get("fs_def") == "FS_min"]
    if tmins:
        notas.append("Instante do $FS$ mínimo (número de realizações): " + "; ".join(tmins) + ".")
    dets = [i for i in ids if casos[i].get("det")]
    if dets:
        tds = sorted({casos[i]["det"]["t_h"] for i in dets})
        quando = (f"o mínimo de todas as rodadas determinísticas em $t = {nm(tds[0])}$ h" if len(tds) == 1 else
                  "; ".join(f"{rotulo(casos[i])}: em $t = {nm(casos[i]['det']['t_h'])}$ h" for i in dets))
        notas.append("Análise determinística: o $FS$ da rodada com os valores médios, sem campos aleatórios, pela "
                     f"mesma definição de $FS$ ({quando}).")
    notas += [_nota_det_versao(ids, casos),
              r"$\mu$, $\sigma$ e $CoV$: média, desvio-padrão e coeficiente de variação da distribuição log-normal "
              r"ajustada ao $FS$ por máxima verossimilhança. $P_f = n_f/S$; $CoV(P_f) = \sqrt{(1 - P_f)/(S\,P_f)}$; "
              r"$S$: número de realizações.",
              _nota_pf0(ids, casos), _nota_negativo(ids, casos)]
    return _ambiente(f"tab:{solo}_tab6", f"Resumo das simulações de Monte Carlo do solo {solo}: intensidade da "
                                         "chuva.", corpo, notas)


def _nota_tmin_casos(ids, casos):
    """instante do FS minimo nos casos com FS minimo: '24 h em todas as realizacoes de ...; caso: 24 h (n), 30 h (m)'"""
    ids = [i for i in ids if casos[i].get("fs_def") == "FS_min"]
    if not ids:
        return ""
    todos = [i for i in ids if len(casos[i].get("t_min_h_contagem", [])) == 1]
    partes = []
    for t in sorted({casos[i]["t_min_h_contagem"][0][0] for i in todos}):
        quem = [rotulo(casos[i]) for i in todos if casos[i]["t_min_h_contagem"][0][0] == t]
        partes.append(f"{num_g(t)} h em todas as realizações de {_lista_pt(quem)}")
    for i in ids:
        if i not in todos:
            partes.append(f"{rotulo(casos[i])}: {_nota_tmin(casos[i])}")
    return "Instante do $FS$ mínimo (número de realizações): " + "; ".join(partes) + "."


def tabA1(solo, S):
    casos = S["casos"]
    ids = _casos_ordem(S)
    linhas = []
    for i in ids:
        c = casos[i]
        a = c.get("amostra") or {}
        linhas.append(" & ".join([rotulo(c), str(c["N"]), str(c["n_f"]), num(a.get("media"), 3), num(a.get("dp"), 3),
                                  pct(a.get("cov"), 1), num_g(a.get("mediana")), num_g(a.get("min")),
                                  num_g(a.get("max")), num(a.get("assimetria"), 2), str(c.get("n_fs_igual_1", "--"))])
                      + r" \\")
    corpo = [
        r"\begin{tabular}{lrrrrrrrrrr}",
        r"\toprule",
        r"Caso & $S$ & $n_f$ & Média & Desvio & $CoV$ (\%) & Mediana & Mín. & Máx. & Assim. & "
        r"$n_{FS=1}$ \\",
        r"\midrule",
    ] + linhas + [r"\bottomrule", r"\end{tabular}"]
    defs = sorted({(casos[i]["fs_def"], tuple(casos[i]["instantes_h"])) for i in ids})
    nfs = []
    for d, inst in defs:
        nota = _nota_fs({"fs_def": d, "instantes_h": list(inst)})
        if len(defs) > 1:
            fams = []
            for i in ids:
                if (casos[i]["fs_def"], tuple(casos[i]["instantes_h"])) == (d, inst):
                    f = NOME_FAMILIA.get(casos[i]["familia"], casos[i]["familia"])
                    if f not in fams:
                        fams.append(f)
            nota = nota[:-1] + " (casos " + _lista_pt(fams) + ")."
        nfs.append(nota)
    dup = _duplicados_ref(S)
    tmin = _nota_tmin_casos(ids, casos)
    notas = [r"Estatísticas amostrais das realizações com $FS > 0$: média, desvio-padrão (divisor $n - 1$), "
             r"$CoV$ = desvio/média, mediana, mínimo, máximo e coeficiente de assimetria amostral ajustado de "
             r"Fisher-Pearson ($G_1$, Assim.). $n_f$: realizações com $FS < 1$; $n_{FS=1}$: realizações com "
             r"$FS = 1$ exatamente.",
             " ".join(nfs),
             tmin,
             _nota_negativo(ids, casos),
             (_lista_pt([rotulo(casos[i]) for i in dup]) + ": os mesmos dados da referência.") if dup else ""]
    return _ambiente(f"tab:{solo}_tabA1", f"Estatísticas amostrais do $FS$ por caso, solo {solo}.", corpo, notas,
                     fonte=r"\footnotesize", colsep="4pt")


def tabA2(solo, S):
    casos = S["casos"]
    ids = _casos_ordem(S)
    mesmo_gl = all((casos[i].get("qui2_grade") or {}).get("normal", {}).get("gl") == casos[i]["qui2"]["normal"]["gl"]
                   and (casos[i].get("qui2_grade") or {}).get("lognormal", {}).get("gl")
                   == casos[i]["qui2"]["lognormal"]["gl"] for i in ids)
    tem_grade = all(casos[i].get("qui2_grade") for i in ids)
    linhas = []
    for i in ids:
        c = casos[i]
        q = c["qui2"]
        cel = [rotulo(c), str(q.get("n", c["N"])), str(q["k"]), str(q["normal"]["gl"])]
        for qq in ([q] + ([c["qui2_grade"]] if tem_grade else [])):
            if qq is not q and not mesmo_gl:
                cel.append(str(qq["normal"]["gl"]))
            xn, xl = qq["normal"]["chi2"], qq["lognormal"]["chi2"]
            cel += [negrito(num(xn, 1), xn < xl), num_p(qq["normal"]["p"]),
                    negrito(num(xl, 1), xl < xn), num_p(qq["lognormal"]["p"])]
        linhas.append(" & ".join(cel) + r" \\")
    ncol_b = (4 if mesmo_gl else 5) if tem_grade else 0
    cab1 = r" & & & & \multicolumn{4}{c}{(a) Classes pelo valor de $F(FS)$}"
    cmid = r"\cmidrule(lr){5-8}"
    cab2 = r" & & & & \multicolumn{2}{c}{Normal} & \multicolumn{2}{c}{Log-normal}"
    cmid2 = r"\cmidrule(lr){5-6}\cmidrule(lr){7-8}"
    cab3 = r"Caso & $n$ & $k$ & gl & $\chi^2$ & $p$ & $\chi^2$ & $p$"
    if tem_grade:
        a = 9
        cab1 += rf" & \multicolumn{{{ncol_b}}}{{c}}{{(b) Bordas na grade do $FS$}}"
        cmid += rf"\cmidrule(lr){{{a}-{a + ncol_b - 1}}}"
        if not mesmo_gl:
            cab2 += " &"
            cab3 += " & gl"
            a += 1
        cab2 += r" & \multicolumn{2}{c}{Normal} & \multicolumn{2}{c}{Log-normal}"
        cmid2 += rf"\cmidrule(lr){{{a}-{a + 1}}}\cmidrule(lr){{{a + 2}-{a + 3}}}"
        cab3 += r" & $\chi^2$ & $p$ & $\chi^2$ & $p$"
    ncols = 8 + ncol_b
    corpo = [
        r"\begin{tabular}{l" + "r" * (ncols - 1) + "}",
        r"\toprule",
        cab1 + r" \\", cmid,
        cab2 + r" \\", cmid2,
        cab3 + r" \\",
        r"\midrule",
    ] + linhas + [r"\bottomrule", r"\end{tabular}"]
    nf = sum(casos[i].get("n_fs_fora_grade", 0) for i in ids)
    nt = sum((casos[i].get("qui2") or {}).get("n", 0) for i in ids)
    h = (casos[ids[0]].get("qui2_grade") or {}).get("h", 0.00625) if ids else 0.00625
    notas = [r"$n$: realizações com $FS > 0$. Distribuições normal ($\mu$ e $\sigma$ de máxima verossimilhança) e "
             r"log-normal (Tabelas 3 a 6) ajustadas ao $FS$; $k = \lceil 2\,n^{0{,}4} \rceil$ classes; "
             r"gl $= k - 3$ graus de liberdade. Em negrito, o menor $\chi^2$ de cada par.",
             r"(a) $k$ classes de igual probabilidade sob a distribuição ajustada $F$: cada $FS$ vai para a classe "
             r"$\lfloor k\,F(FS) \rfloor$; frequência esperada $n/k$.",
             (rf"(b) As bordas das classes são os pontos médios da grade do $FS$ (múltiplos ímpares de "
              rf"{num_g(h / 2)}) mais próximos dos quantis $j/k$ da distribuição ajustada; a frequência esperada de "
              rf"cada classe é $n$ vezes a probabilidade ajustada dela"
              + ("; o número de classes é $k$ em todos os casos." if mesmo_gl else
                 "; bordas repetidas se fundem (gl = classes $- 3$).")) if tem_grade else "",
             (rf"O $FS$ do SRM (bissecção com resolução {num_g(h)}) é múltiplo de {num_g(h)} em {nt - nf} "
              rf"dos {nt} valores desta tabela." if tem_grade else ""),
             (_lista_pt([rotulo(casos[i]) for i in _duplicados_ref(S)]) + ": os mesmos dados da referência.")
             if _duplicados_ref(S) else ""]
    return _ambiente(f"tab:{solo}_tabA2", "Teste qui-quadrado de aderência das distribuições normal e log-normal "
                                          f"ajustadas ao $FS$, solo {solo}.", corpo, notas,
                     fonte=r"\footnotesize", colsep="3pt")


def tabA3(solo, S):
    c = S["casos"][S["ref"]]
    sp = c.get("spearman")
    nomes = (("c", "$c'$", "média na banda de ruptura"),
             ("phi", r"$\varphi'$", "média na banda de ruptura"),
             ("k", "$k_v$", "média na massa que desliza"))
    linhas = []
    for chave, simb, desc in nomes:
        if not sp or sp["segunda_ordem"].get(chave) is None:
            linhas.append(f"{simb} & {desc} & -- & -- & -- & -- \\\\")
            continue
        a, b = sp["ic95"][chave]
        a1, b1 = sp["ic95_simples"][chave]
        linhas.append(f"{simb} & {desc} & {num(sp['segunda_ordem'][chave], 3)} & [{num(a, 3)}; {num(b, 3)}] & "
                      f"{num(sp['simples'][chave], 3)} & [{num(a1, 3)}; {num(b1, 3)}] \\\\")
    corpo = [
        r"\begin{tabular}{llcccc}",
        r"\toprule",
        r" & & \multicolumn{2}{c}{Parcial de 2ª ordem} & \multicolumn{2}{c}{Simples} \\",
        r"\cmidrule(lr){3-4}\cmidrule(lr){5-6}",
        r"$X$ & Variável equivalente & $\rho_{X\,FS,YZ}$ & IC 95\% & $\rho_{X\,FS}$ & IC 95\% \\",
        r"\midrule",
    ] + linhas + [r"\bottomrule", r"\end{tabular}"]
    notas = [(rf"$n = {sp['n']}$ realizações (as com $FS > 0$). " if sp else "")
             + r"Correlação de postos de Spearman entre o $FS$ e as variáveis de entrada equivalentes de cada "
               r"realização, no instante do $FS$: $c'$ e $\varphi'$ médios na banda de ruptura (pontos com "
               r"$\sqrt{J_2}$ da deformação plástica $\geq 0{,}1\times$ o máximo, no último estado convergido do "
               r"SRM) e $k_v$ médio na massa que desliza (pontos com $|\Delta u| \geq 0{,}5\times$ o máximo).",
             r"$\rho_{X\,FS,YZ}$: correlação parcial de 2ª ordem entre $X$ e o $FS$, controlando as outras duas "
             r"variáveis (eqs. 70 e 71 do Ceron); $\rho_{X\,FS}$: correlação simples. IC 95\%: intervalo de 95\% "
             r"pela transformação de Fisher, com erro-padrão $1/\sqrt{n - 3 - q}$, $q = 2$ na parcial e $q = 0$ "
             r"na simples. " + _nota_fs(c)]
    return _ambiente(f"tab:{solo}_tabA3", "Correlação de postos de Spearman entre o $FS$ e as variáveis de entrada "
                                          f"equivalentes, solo {solo}, configuração de referência.", corpo, notas)


# ---------------------------------------------------------------------------------------------------- proveniencia
def _linhas_prov():
    """{campanha: [(pastas, uso)]} na ordem da tabela; uso None = nao usada (o motivo vem dos dados). O que estiver nos
    dados e nao aqui vai no fim da sua campanha, uma linha por pasta"""
    s50 = [f"f50_I36_i_{s}" for s in SUFIXOS]
    s100 = [f"f100_I36_i_{s}" for s in SUFIXOS]
    ch100 = ["f100_I3.6_i", "f100_I9_i", "f100_I80_i"]
    return {
        "sensibilidade": [
            (["f50_I36_i"], r"f50: referência"),
            (s50, r"f50: $CoV$, $s$ e $\alpha$")],
        "monte_carlo": [
            (["f50_I9_i_ftrans0.5"], r"f50: $I = 9$ mm/h"),
            (["f50_I36_i_ftrans0.5"], r"f50: $I = 36$ mm/h"),
            (["f50_I80_i_ftrans0.5"], r"f50: $I = 80$ mm/h"),
            (["f50_I36_i_ftrans0.5_det"], r"f50: determinístico (referência e $I = 36$ mm/h)"),
            (["f50_I9_i_ftrans0.5_det", "f50_I80_i_ftrans0.5_det"], r"f50: determinísticos ($I = 9$ e 80 mm/h)"),
            (["f50_I36_i", "f50_I36_i_det"], None),
            (["f50_I36_i_ftrans0.5_kramin0.01", "f50_I36_i_ftrans0.5_kramin0.01_det"], None)],
        "f100": [
            (["f100_I36_i"], r"f100: referência e $I = 36$ mm/h (com f100\_fase1)"),
            (s100, r"f100: $CoV$, $s$ e $\alpha$"),
            (ch100, r"f100: $I = 3{,}6$, 9 e 80 mm/h (com f100\_fase3)"),
            (["f100_I36_i_det", "f100_I36_i_alfa2_det", "f100_I36_i_alfa5_det", "f100_I36_i_alfa10_det"],
             r"f100: determinísticos (referência e $\alpha = 2$, 5 e 10)"),
            (["f100_I3.6_i_det", "f100_I9_i_det", "f100_I80_i_det"],
             r"f100: determinísticos ($I = 3{,}6$, 9 e 80 mm/h)"),
            (["f100_I36_i_kramin0.01", "f100_I36_i_kramin0.01_det"], None)],
        "f100_fase1": [
            (["f100_I36_i"], r"f100: referência e $I = 36$ mm/h (com f100)")],
        "f100_fase3": [
            (ch100, r"f100: $I = 3{,}6$, 9 e 80 mm/h (com f100)"),
            (["f100_I36_i_kramin0.01"], None)],
    }


def _faixa(fx):
    if not fx:
        return "--"
    return str(fx[0]) if fx[0] == fx[1] else f"{fx[0]}--{fx[1]}"


def _instantes_txt(inst):
    """instantes do FS: '0 a 48, de 6 em 6' quando sao mais de 4 e igualmente espacados; senao a lista"""
    inst = list(inst)
    if len(inst) > 4:
        passos = {round(b - a, 9) for a, b in zip(inst, inst[1:])}
        if len(passos) == 1:
            p = passos.pop()
            return f"{num_g(inst[0])} a {num_g(inst[-1])}, de {num_g(p)} em {num_g(p)}"
    return _lista_pt([num_g(t) for t in inst])


def _motivo_curto(m):
    """'F_TRANS = 1 (...): sem analogo no Ceron' -> 'F_TRANS = 1'"""
    return m.split(":")[0].split("(")[0].strip()


def _linha_prov(combs, rs, uso):
    if len(combs) == len(SUFIXOS) and len(rs) > 1:
        pasta = tex(os.path.commonprefix(combs)) + r"\textit{sufixo}"
    else:
        nomes = [c for c in combs if any(r["combinacao"] == c for r in rs)]
        pasta = ", ".join(tex(c) for c in nomes)
    faixas = {tuple(r["faixa_reais"]) if r.get("faixa_reais") else None for r in rs}
    faixa = _faixa(next(iter(faixas))) if len(faixas) == 1 else "; ".join(_faixa(r.get("faixa_reais")) for r in rs)
    concl = sorted({r["concluidas"] for r in rs})
    concl_txt = (str(concl[0]) + (" cada" if len(rs) > 1 else "")) if len(concl) == 1 else \
        "; ".join(str(r["concluidas"]) for r in rs)
    fal = [f"{len(r['falhou'])} (" + ", ".join(str(n) for n in r["falhou"]) + ")" for r in rs if r.get("falhou")]
    fal_txt = "; ".join(fal) if fal else "0"
    vers = sorted({v for r in rs for v in r.get("versoes", {})})
    ver_txt = ", ".join(str(v).split("_")[0] for v in vers) if vers else "--"
    if uso is None:
        mot = sorted({_motivo_curto(r["excluida"]) for r in rs if r.get("excluida")})
        uso = "não usada (" + "; ".join(tex(m) for m in mot) + ")" if mot else "não usada"
    return (rf"\raggedright {pasta} & {faixa} & {concl_txt} & \raggedright {fal_txt} & {ver_txt} & "
            rf"\raggedright {uso} \tabularnewline")


def tab_proveniencia(est):
    prov = est.get("proveniencia", {})
    reg = {p["pasta"]: p for p in prov.get("pastas", [])}
    camps = prov.get("campanhas", {})
    especificadas = _linhas_prov()
    ordem = list(especificadas) + sorted({r["campanha"] for r in reg.values()} - set(especificadas))
    usados = set()
    linhas = []
    sufixos_usados = False
    for camp in ordem:
        rs_camp = [r for r in reg.values() if r["campanha"] == camp]
        if not rs_camp:
            continue
        ct = camps.get(camp, {})
        inst = sorted({tuple(r.get("instantes_h", [])) for r in rs_camp})
        cab = rf"\textbf{{{tex(camp)}}}"
        info = []
        if len(inst) == 1:
            info.append(f"$FS$ em {_instantes_txt(inst[0])} h")
        if ct.get("VARIANTE"):
            info.append(f"variante {tex(ct['VARIANTE'])}")
        if info:
            cab += " (" + "; ".join(info) + ")"
        linhas += [(r"\midrule" if linhas else ""), rf"\multicolumn{{6}}{{l}}{{{cab}}} \tabularnewline"]
        for combs, uso in especificadas.get(camp, []):
            rs = [reg[f"{camp}/{c}"] for c in combs if f"{camp}/{c}" in reg]
            if not rs:
                continue
            usados.update(r["pasta"] for r in rs)
            linhas.append(_linha_prov(combs, rs, uso))
            sufixos_usados |= len(combs) == len(SUFIXOS) and len(rs) > 1
        for r in sorted(rs_camp, key=lambda r: r["combinacao"]):  # o que nao estava previsto
            if r["pasta"] not in usados:
                usados.add(r["pasta"])
                linhas.append(_linha_prov([r["combinacao"]], [r], None if r.get("excluida") else
                                          "uso: " + ", ".join(tex(u) for u in r.get("uso", []))))
    linhas = [l for l in linhas if l]
    corpo = [
        r"\begin{tabular}{p{5.1cm}cc p{2.0cm}l p{3.6cm}}",
        r"\toprule",
        r"Pasta & Realizações & Concluídas & \texttt{falhou.txt} & Versão & Uso \tabularnewline",
    ] + linhas + [r"\bottomrule", r"\end{tabular}"]
    notas = []
    if sufixos_usados:
        notas.append(r"\textit{sufixo}: " + ", ".join(tex(s) for s in SUFIXOS) + r" (16 pastas; covk, covc e "
                     r"covphi: $CoV$ de $k_v$, $c'$ e $\varphi'$; lx$a$\_ly$b$: $(L_x, L_y) = (a, b)$ m; alfa: "
                     r"$\alpha$).")
    notas.append(r"Realizações: os números das pastas \texttt{r}$nnnnn$. Concluída: \texttt{resultado.txt} não "
                 r"vazio; todas as concluídas das pastas usadas têm situação ok. \texttt{falhou.txt}: realização que "
                 r"terminou com erro três vezes, fora de $S$ (números entre parênteses). Versão: o commit do programa "
                 r"(\texttt{versao} em \texttt{resumo.csv}). $F_{TRANS} = 0{,}5$ em todas as pastas usadas.")
    vers = sorted({v for r in reg.values() for v in r.get("versoes", {})})
    if vers:
        notas.append("Versões: " + "; ".join(_versao(v).replace(" (", ", compilada em ").rstrip(")")
                                             for v in vers) + ".")
    fus = prov.get("fusoes", [])
    if fus:
        grupos = {}  # (campanhas, S) -> combinacoes
        for f in fus:
            camps_f = tuple(p.split("/")[0] for p in f["pastas"])
            grupos.setdefault((camps_f, sum(f["N_por_pasta"])), []).append(f["pastas"][0].split("/", 1)[1])
        itens = [" + ".join(tex(c) for c in cs) + ": " + _lista_pt([tex(x) for x in combs])
                 + f" ($S = {n}$" + (" cada)" if len(combs) > 1 else ")") for (cs, n), combs in grupos.items()]
        notas.append(r"Pastas juntadas num caso, com \texttt{parametros.txt} idêntico e números de realização "
                     "disjuntos: " + "; ".join(itens) + ".")
    corr = {}
    for S in est.get("solos", {}).values():
        for c in S["casos"].values():
            for x in c.get("fs_csv_corrompido", []):
                corr[(x["pasta"], x["real"])] = x
    for (p, n), x in sorted(corr.items()):
        notas.append(rf"{tex(p)}/r{n:05d}: \texttt{{fs.csv}} com {x['linhas_corrompidas']} linhas corrompidas "
                     rf"(caracteres não imprimíveis); o $FS$ de {_lista_pt([num_g(t) for t in x['instantes_de_resumo']])}"
                     rf" h vem de \texttt{{resumo.csv}}, íntegro, que confere com \texttt{{resultado.txt}} no $FS$ "
                     rf"mínimo e no seu instante.")
    for p, r in reg.items():
        extra = []
        if r.get("sem_resultado"):
            extra.append(r"realizações sem \texttt{resultado.txt} nem \texttt{falhou.txt}: "
                         + ", ".join(str(n) for n in r["sem_resultado"]))
        if r.get("concluidas_sem_FS"):
            extra.append(r"\texttt{resultado.txt} sem $FS$: " + ", ".join(str(n) for n in r["concluidas_sem_FS"]))
        if extra:
            notas.append(f"{tex(p)}: " + "; ".join(extra) + ".")
    comp = prov.get("comparacao_f50_ref")
    if comp:
        notas.append(rf"{tex(comp['pastas'][0])} (versão {', '.join(str(v).split('_')[0] for v in comp['versoes'][0])})"
                     rf" e {tex(comp['pastas'][1])} (versão "
                     rf"{', '.join(str(v).split('_')[0] for v in comp['versoes'][1])}), realizações "
                     rf"{comp['faixa_comum'][0]} a {comp['faixa_comum'][1]}: $FS$ em $t = 24$ h idêntico em "
                     rf"{comp['iguais']['FS_24h']} de {comp['reais_comuns']}; médias de $c'$, $\varphi'$ e $k$ no "
                     rf"solo idênticas em {min(comp['iguais'][k] for k in ('c_med', 'phi_med', 'k_mult_med'))} de "
                     rf"{comp['reais_comuns']}.")
    excl = sorted({r["excluida"] for r in reg.values() if r.get("excluida")})
    if excl:
        notas.append("Não usadas: " + "; ".join(tex(m) for m in excl) + ".")
    solos_min = [solo for solo, S in est.get("solos", {}).items()
                 if any(c.get("fs_def") == "FS_min" for c in S["casos"].values())]
    if solos_min:
        refs = [rf"\ref{{tab:{solo}_tabA1}}" for solo in solos_min]
        notas.append("Instante do $FS$ mínimo em cada caso: notas " + ("da Tabela " if len(refs) == 1 else
                                                                         "das Tabelas ") + _lista_pt(refs) + ".")
    return _ambiente("tab:proveniencia", "Proveniência dos dados: campanhas, pastas, realizações, versões do "
                                         "programa e combinações não usadas.", corpo, notas,
                     fonte=r"\scriptsize", colsep="3pt", fonte_notas=r"\scriptsize")


# ------------------------------------------------------------------------------------------------------------ gera
def gera(est, saida):
    if isinstance(est, str):  # caminho do estatisticas.json
        with open(est, encoding="utf-8") as f:
            est = json.load(f)
    gravados = []
    prov = est.get("proveniencia", {})
    for solo, S in est["solos"].items():
        _grava(saida, f"{solo}_tab2.tex", tab2(solo, S, prov), gravados)
        _grava(saida, f"{solo}_tab3.tex", tab3(solo, S), gravados)
        _grava(saida, f"{solo}_tab4.tex", tab4(solo, S), gravados)
        _grava(saida, f"{solo}_tab5.tex", tab5(solo, S), gravados)
        _grava(saida, f"{solo}_tab6.tex", tab6(solo, S, prov), gravados)
        _grava(saida, f"{solo}_tabA1.tex", tabA1(solo, S), gravados)
        _grava(saida, f"{solo}_tabA2.tex", tabA2(solo, S), gravados)
        _grava(saida, f"{solo}_tabA3.tex", tabA3(solo, S), gravados)
    _grava(saida, "tab_proveniencia.tex", tab_proveniencia(est), gravados)
    return gravados


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("uso: python3 -I tabelas.py <saida> (le <saida>/estatisticas.json)")
        sys.exit(2)
    for a in gera(os.path.join(sys.argv[1], "estatisticas.json"), sys.argv[1]):
        print(a)
