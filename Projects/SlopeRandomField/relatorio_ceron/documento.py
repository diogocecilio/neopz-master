#!/usr/bin/env python3
"""documento.py: o relatorio em LaTeX com os resultados das campanhas de Monte Carlo do SlopeRandomField no formato da
secao 6 (e da secao 5.2) de Vargas Ceron et al. (2025), "Stability Analysis of Slope Subjected to Seepage Forces
Considering Spatial Variability of Soil Properties", Int J Numer Anal Methods Geomech 49:2459-2491 (doi
10.1002/nag.3993).

gera(est, saida) escreve <saida>/relatorio.tex e o compila com latexmk -pdf em <saida>; devolve o caminho do PDF.
est e o dicionario de dados.calcula (ou o caminho de estatisticas.json). O documento usa, por estes nomes:
- <saida>/tabelas/<solo>_tab{2,3,4,5,6,A1,A2,A3}.tex e tab_proveniencia.tex (tabelas.py);
- <saida>/figuras/<solo>_fig{12,13,14,16,24,25,...,33}.pdf (figuras_campos.py e figuras_estat.py);
- as notas de legenda de <saida>/figuras/legendas_estat.json e legendas_campos.json e as porcentagens dos modos de
  ruptura de <saida>/figuras/modos_ruptura.json.
Cada \\includegraphics e cada \\input fica num \\IfFileExists com um aviso visivel no lugar do arquivo ausente: o
documento compila mesmo sem algum deles.

Regras do texto: so fatos (valores, contagens, variacoes que aparecem nos numeros, definicoes e metodos), na ordem do
texto do Ceron, sem conjecturas, projecoes ou interpretacoes. Todo numero de resultado vem de est ou das notas JSON e
e formatado aqui, com virgula decimal."""

import json
import math
import os
import re
import shutil
import subprocess
import sys
import time

_AQUI = os.path.dirname(os.path.abspath(__file__))
_PAI = os.path.dirname(_AQUI)  # Projects/SlopeRandomField
for _p in (_PAI, _AQUI):
    if _p not in sys.path:
        sys.path.insert(0, _p)

try:
    import eixos_ceron  # noqa: E402  (valores das familias nas figuras do Ceron)
except Exception:  # pragma: no cover - o documento sai sem as comparacoes de valores com o Ceron
    eixos_ceron = None

NOME_TEX = "relatorio"
SOLOS = ("f50", "f100")

# dados fixos do modelo (README e main.cpp); usados so se a nota da Fig. 24 (legendas_campos.json) nao os trouxer
GEOMETRIA_PADRAO = {"H": 10.0, "beta": 45.0, "crista": (30.0, 30.0), "pe": (40.0, 20.0), "dom": (70.0, 30.0),
                    "area": 1750.0}
# parametros do metodo (main.cpp: BANDA_LIMIAR, MASSA_LIMIAR, reducao da carga de gravidade); usados so se o texto de
# est["proveniencia"]["metodos"] nao os trouxer
METODO_PADRAO = {"banda": 0.1, "massa": 0.5, "c_grav": 0.6, "h": 0.00625}

# o que cada figura mostra, nas palavras do Ceron (legendas traduzidas)
LEGENDAS = {
    "12": ("Realização típica dos campos aleatórios com $L_x = <<Lx>>$~m e $L_y = <<Ly>>$~m: (a) campo de "
           "permeabilidade vertical; (b) campo de coesão; (c) campo de ângulo de atrito",
           "Realização típica dos campos aleatórios"),
    "13": ("Distribuição da poropressão obtida numericamente com a realização típica do campo aleatório de "
           "permeabilidade", "Poropressão da realização típica"),
    "14": ("Superfície crítica associada a um ciclo típico de simulação de Monte Carlo",
           "Superfície crítica de uma realização típica"),
    "16": ("Superfícies críticas de ruptura associadas às realizações dos campos aleatórios de coesão, ângulo de "
           "atrito e permeabilidade", "Superfícies críticas das realizações"),
    "24": ("Geometria do talude e superfície crítica de ruptura obtida da análise determinística",
           "Geometria e superfície crítica determinística"),
    "25": ("Resultados das simulações de Monte Carlo na configuração de referência: (a) função densidade de "
           "probabilidade do fator de segurança; (b) convergência da probabilidade de falha",
           "Configuração de referência: densidade e convergência"),
    "26": ("Coeficiente de correlação de postos de Spearman parcial de segunda ordem para o fator de segurança",
           "Spearman parcial de segunda ordem"),
    "27": ("Efeito do coeficiente de variação das propriedades do solo na função densidade de probabilidade do fator "
           "de segurança: (a) permeabilidade vertical; (b) coesão; (c) ângulo de atrito",
           "Coeficiente de variação: densidade do fator de segurança"),
    "28": ("Probabilidade de falha em função do coeficiente de variação da permeabilidade vertical (a), da coesão (b) "
           "e do ângulo de atrito (c)", "Coeficiente de variação: probabilidade de falha"),
    "29": ("Valor médio do fator de segurança em função do coeficiente de variação da permeabilidade vertical (a), da "
           "coesão (b) e do ângulo de atrito (c)", "Coeficiente de variação: média do fator de segurança"),
    "30": ("Coeficiente de variação do fator de segurança em função do coeficiente de variação da permeabilidade "
           "vertical (a), da coesão (b) e do ângulo de atrito (c)",
           "Coeficiente de variação: coeficiente de variação do fator de segurança"),
    "31": ("Efeito da distância de autocorrelação na função densidade de probabilidade do fator de segurança (a), na "
           "probabilidade de falha global (b), no coeficiente de variação do fator de segurança (c) e no valor médio "
           "do fator de segurança (d)", "Distâncias de autocorrelação"),
    "32": ("Efeito da anisotropia da permeabilidade na função densidade de probabilidade do fator de segurança (a), "
           "na probabilidade de falha global (b), no coeficiente de variação do fator de segurança (c) e no valor "
           "médio do fator de segurança (d)", "Anisotropia da permeabilidade"),
    "33": ("Efeito da intensidade da chuva na função densidade de probabilidade do fator de segurança (a), na "
           "probabilidade de falha global (b), no coeficiente de variação do fator de segurança (c) e no valor médio "
           "do fator de segurança (d)", "Intensidade da chuva"),
}
# a figura do Ceron de que cada uma e analoga (texto depois de "analoga a")
ANALOGA = {"16": ("às Figs.~16 e 20", ""),
           "33": ("à Fig.~33", ", em que a variável é a razão de rebaixamento $h_w/H$")}

FAMILIAS = ("covk", "covc", "covphi", "s", "alfa", "chuva")
VAR = {"covk": r"CoV(k_v)", "covc": r"CoV(c')", "covphi": r"CoV(\varphi')", "s": "s", "alfa": r"\alpha",
       "chuva": "I"}
UNID = {"covk": r"\%", "covc": r"\%", "covphi": r"\%", "s": "", "alfa": "", "chuva": "~mm/h"}
ESCALA = {"covk": 100.0, "covc": 100.0, "covphi": 100.0, "s": 1.0, "alfa": 1.0, "chuva": 1.0}
CAMPO_PADRAO = {"covk": "cov_k", "covc": "cov_c", "covphi": "cov_phi", "s": "s", "alfa": "alfa", "chuva": "I_mm_h"}
NOME_PROP = {"covk": "permeabilidade vertical $k_v$", "covc": "coesão $c'$", "covphi": r"ângulo de atrito $\varphi'$"}

AVISOS = []


# ================================================================================================= formatos numericos
def _ok(v):
    return isinstance(v, (int, float)) and not isinstance(v, bool) and math.isfinite(v)


def num(v, casas):
    """numero em modo texto, com virgula decimal e casas fixas; None -> '--'"""
    if not _ok(v):
        return "--"
    s = f"{v:.{casas}f}"
    neg = s.startswith("-")
    if neg:
        s = s[1:]
    s = s.replace(".", ",")
    if neg and s.strip("0,") != "":
        return "$-$" + s
    return s


def numg(v, casas=7):
    """numero em modo texto com os seus digitos (ate 'casas' casas, sem zeros a direita)"""
    if not _ok(v):
        return "--"
    s = f"{v:.{casas}f}".rstrip("0").rstrip(".")
    if s in ("-0", ""):
        s = "0"
    if s.startswith("-"):
        return "$-$" + s[1:].replace(".", ",")
    return s.replace(".", ",")


def mat(s):
    """converte um numero em modo texto (de num/numg) para o modo matematico"""
    return s.replace("$-$", "-").replace(",", "{,}")


def pct(v, casas=2):
    return (num(100.0 * v, casas) + r"\%") if _ok(v) else "--"


def txt_digitos(s):
    """o texto de um numero como esta no arquivo (os digitos exatos), com virgula decimal"""
    s = str(s).strip()
    if s.startswith("-"):
        return "$-$" + s[1:].replace(".", ",")
    return s.replace(".", ",")


def fs_det_txt(det):
    if not det:
        return "--"
    if det.get("FS_txt"):
        return txt_digitos(det["FS_txt"])
    return numg(det.get("FS"))


def lista_pt(itens, sep=", "):
    itens = [str(i) for i in itens]
    if not itens:
        return ""
    if len(itens) == 1:
        return itens[0]
    return sep.join(itens[:-1]) + " e " + itens[-1]


def lista_num(itens):
    """lista de numeros ja formatados: separados por ';' quando algum tem virgula decimal"""
    itens = [str(i) for i in itens]
    return lista_pt(itens, "; " if any("," in i for i in itens) else ", ")


def num_p(p):
    """valor p para o modo matematico: 3 casas a partir de 0,001; abaixo, notacao cientifica"""
    if not _ok(p):
        return r"\textrm{--}"
    if p >= 0.001:
        return mat(num(p, 3))
    if p == 0:
        return "0"
    e = int(math.floor(math.log10(p)))
    m = p / 10 ** e
    if round(m, 1) >= 10:
        m, e = m / 10, e + 1
    return f"{m:.1f}".replace(".", "{,}") + r"\times 10^{" + str(e) + "}"


def sci_m(v, alg=2):
    """numero positivo em notacao cientifica para o modo matematico: 1{,}02\\times 10^{-6}"""
    if not _ok(v) or v <= 0:
        return r"\textrm{--}"
    e = int(math.floor(math.log10(v)))
    m = v / 10 ** e
    if round(m, alg) >= 10:
        m, e = m / 10, e + 1
    ms = f"{m:.{alg}f}".rstrip("0").rstrip(".")
    if ms == "1":
        return r"10^{" + str(e) + "}"
    return ms.replace(".", "{,}") + r"\times 10^{" + str(e) + "}"


def horas(t):
    return numg(t, 3)


# ======================================================================================================= texto/LaTeX
_ESPECIAIS = {"\\": r"\textbackslash{}", "{": r"\{", "}": r"\}", "%": r"\%", "&": r"\&", "#": r"\#", "_": r"\_",
              "^": r"\^{}", "~": r"\textasciitilde{}", "$": r"\$", "<": r"\textless{}", ">": r"\textgreater{}"}
_UNICODE = {
    "−": "$-$", "–": "--", "—": "---", "×": r"$\times$", "±": r"$\pm$", "≥": r"$\geq$", "≤": r"$\leq$",
    "≈": r"$\approx$", "≠": r"$\neq$", "°": r"$^\circ$", "²": r"$^2$", "³": r"$^3$", "¹": r"$^1$", "⁻": r"$^-$",
    "ª": r"\textsuperscript{a}", "º": r"\textsuperscript{o}", "µ": r"$\mu$", "μ": r"$\mu$", "σ": r"$\sigma$",
    "φ": r"$\varphi$", "ϕ": r"$\phi$", "α": r"$\alpha$", "β": r"$\beta$", "γ": r"$\gamma$", "Γ": r"$\Gamma$",
    "ρ": r"$\rho$", "Δ": r"$\Delta$", "δ": r"$\delta$", "χ": r"$\chi$", "θ": r"$\theta$", "λ": r"$\lambda$",
    "ε": r"$\varepsilon$", "η": r"$\eta$", "π": r"$\pi$", "τ": r"$\tau$", "ψ": r"$\psi$", "ω": r"$\omega$",
    "Ω": r"$\Omega$", "ν": r"$\nu$", "κ": r"$\kappa$", "ξ": r"$\xi$", "√": r"$\surd$", "→": r"$\rightarrow$",
    "←": r"$\leftarrow$", "⇔": r"$\Leftrightarrow$", "⇒": r"$\Rightarrow$", "∞": r"$\infty$", "·": r"$\cdot$",
    "⌈": r"$\lceil$", "⌉": r"$\rceil$", "⌊": r"$\lfloor$", "⌋": r"$\rfloor$", "…": r"\dots{}", "“": "``",
    "”": "''", "‘": "`", "’": "'", "\u00a0": "~", "∈": r"$\in$", "∑": r"$\sum$", "∫": r"$\int$",
}


def tex(s):
    """texto simples -> LaTeX: escapa os caracteres especiais e troca os simbolos Unicode que o pdflatex nao tem"""
    out = []
    for ch in str(s):
        if ch in _ESPECIAIS:
            out.append(_ESPECIAIS[ch])
        elif ch in _UNICODE:
            out.append(_UNICODE[ch])
        elif ord(ch) < 128 or (0xC0 <= ord(ch) <= 0xFF and ch not in "×÷"):
            out.append(ch)
        else:
            AVISOS.append(f"caractere sem tradução para LaTeX trocado por '?': {ch!r} (U+{ord(ch):04X})")
            out.append("?")
    return "".join(out).replace("$$", "${}$")


def nota_tex(s):
    """nota de legenda -> LaTeX. As notas sao texto simples; se tiverem trechos entre $...$ (em numero par), esses
    trechos sao mantidos como modo matematico e o resto e escapado. Sem $, os simbolos escritos em texto simples
    (FS, Pf, nf, N, kv, phi, CoV(c), 1e-5, -50, ...) viram a mesma notacao matematica do texto (_nota_plana_tex)"""
    s = str(s).strip()
    if s.count("$") >= 2 and s.count("$") % 2 == 0:
        partes = s.split("$")
        ok = all(p.count("{") == p.count("}") for p in partes[1::2])
        if ok:
            return "".join(tex(p) if i % 2 == 0 else "$" + p + "$" for i, p in enumerate(partes))
    return _nota_plana_tex(s)


# simbolos das notas em texto simples (figuras_estat.py e figuras_campos.py) -> modo matematico, com a notacao do
# texto: S (numero de realizacoes, eq. 68 do Ceron) no lugar de N, \FS, P_f, n_f, k_v, \varphi, ...
_NUM_NOTA = r"-?\d+(?:,\d+)?(?:e[-+]?\d+)?"
_SIMBOLOS_NOTA = [  # (texto, LaTeX), os mais longos primeiro
    ("alfa = kh/kv", r"\alpha = k_h/k_v"), ("Pf = nf/N", r"P_f = n_f/S"), ("CoV(kv)", r"CoV(k_v)"),
    ("CoV(phi)", r"CoV(\varphi)"), ("CoV(c)", r"CoV(c)"), ("raiz de J2", r"\sqrt{J_2}"), ("hw/H", r"h_w/H"),
    ("y min", r"y_{\min}"), ("CoV", "CoV"), ("FS", r"\FS"), ("Pf", "P_f"), ("nf", "n_f"), ("kv", "k_v"),
    ("kh", "k_h"), ("ks", "k_s"), ("pw", "p_w"), ("phi", r"\varphi"), ("alfa", r"\alpha"), ("Lx", "L_x"),
    ("Ly", "L_y"), ("N", "S"), ("n", "n"),
]
_RE_NOTA = re.compile(
    # (L_x, L_y) = s (20 m, 2 m)
    r"(?P<ls>\(Lx, Ly\) = s \((?P<lx>\d+(?:,\d+)?) m, (?P<ly>\d+(?:,\d+)?) m\))"
    # simbolo, talvez seguido de '= numero', '< numero' ou '> numero'
    r"|(?<![\w'])(?P<sym>" + "|".join(re.escape(a) for a, _ in _SIMBOLOS_NOTA) + r")(?![\w'])"
    r"(?P<rel>\s*[=<>]\s*(?P<num>" + _NUM_NOTA + r")(?![\w]|,\d))?"
    # letras isoladas so com '= numero' (I = 36, t = 24, H = 10, s = 400) ou c e k nos contextos das notas
    r"|(?<![\w'])(?P<let>[IHts])\s*=\s*(?P<num1>" + _NUM_NOTA + r")(?![\w]|,\d)"
    r"|(?<![\w'])(?P<let2>I(?= \(mm/h\))|k(?= \(m/s\))|c(?=\s*=|:))"
    # numero em notacao cientifica (2,01e-7) ou negativo (-3,4) solto no texto
    r"|(?<![\w,.])(?P<sci>-?\d+(?:,\d+)?e[-+]?\d+)(?![\w])"
    r"|(?<![\w,.-])(?P<neg>-\d+(?:,\d+)?)(?![\w])"
    # unidades com expoente negativo (kPa-1 s-1)
    r"|(?<![\w])(?P<uni>kPa|s|m)-1(?![\w])"
    # nomes de arquivos
    r"|(?<![\w/.])(?P<arq>[\w.-]+\.(?:csv|txt))(?![\w])")


def _num_nota_mat(t):
    """numero de uma nota (virgula decimal, talvez com expoente 'e-7') -> modo matematico"""
    m = re.fullmatch(r"(-?)(\d+(?:,\d+)?)(?:e([-+]?\d+))?", t)
    if not m:
        return tex(t)
    sinal, mant, exp = m.groups()
    if exp is None:
        return sinal + mant.replace(",", "{,}")
    e = str(int(exp))
    if mant == "1":
        return sinal + r"10^{" + e + "}"
    return sinal + mant.replace(",", "{,}") + r"\times 10^{" + e + "}"


def _nota_plana_tex(s):
    """nota em texto simples -> LaTeX: escapa o texto e poe os simbolos e numeros especiais em modo matematico"""
    simb = dict(_SIMBOLOS_NOTA)
    out, i = [], 0
    for m in _RE_NOTA.finditer(s):
        out.append(tex(s[i:m.start()]))
        i = m.end()
        g = m.groupdict()
        if g["ls"]:
            out.append(rf"$(L_x, L_y) = s\,({_num_nota_mat(g['lx'])}~\mathrm{{m}}, "
                       rf"{_num_nota_mat(g['ly'])}~\mathrm{{m}})$")
        elif g["sym"]:
            mm = simb[g["sym"]]
            if g["rel"]:
                op = g["rel"].strip()[0]
                mm += f" {op} " + _num_nota_mat(g["num"])
            out.append("$" + mm + "$")
        elif g["let"]:
            out.append("$" + g["let"] + " = " + _num_nota_mat(g["num1"]) + "$")
        elif g["let2"]:
            out.append("$" + g["let2"] + "$")
        elif g["sci"]:
            out.append("$" + _num_nota_mat(g["sci"]) + "$")
        elif g["neg"]:
            out.append("$-$" + g["neg"][1:])
        elif g["uni"]:
            out.append(g["uni"] + "$^{-1}$")
        elif g["arq"]:
            out.append(r"\texttt{" + tex(g["arq"]) + "}")
    out.append(tex(s[i:]))
    # dois trechos matematicos vizinhos viram um so
    return "".join(out).replace("$$", "")


def pasta_tt(p):
    """nome de pasta em \\texttt, com quebras de linha permitidas depois de / e _"""
    s = tex(p).replace(r"\_", r"\_\allowbreak{}").replace("/", r"/\allowbreak{}")
    return r"\texttt{" + s + "}"


def T(modelo, **kw):
    """modelo com campos <<nome>>"""
    s = modelo
    for k, v in kw.items():
        s = s.replace("<<" + k + ">>", str(v))
    m = re.search(r"<<(\w+)>>", s)
    if m:
        raise KeyError(f"campo sem valor no modelo: {m.group(1)}")
    return s


def frase(s):
    """termina a frase com ponto"""
    s = s.strip()
    return s if s.endswith((".", ":", "?", "!")) else s + "."


def maiuscula(s):
    """primeira letra maiuscula (so se o texto comeca com uma letra; o resto fica como esta)"""
    return s[0].upper() + s[1:] if s and s[0].isalpha() else s


# ================================================================================================== leitura de dados
def _carrega(est):
    if isinstance(est, (str, os.PathLike)):
        with open(est, encoding="utf-8") as f:
            return json.load(f)
    return est


def _le_json(arq):
    try:
        with open(arq, encoding="utf-8") as f:
            return json.load(f)
    except FileNotFoundError:
        return None
    except Exception as e:  # arquivo ilegivel: o documento sai sem ele
        AVISOS.append(f"{arq}: {e}")
        return None


def _num_pt(s):
    return float(s.replace(",", "."))


def _geometria(notas):
    """geometria do talude lida da nota da Fig. 24 (legendas_campos.json); o que faltar vem de GEOMETRIA_PADRAO"""
    g = dict(GEOMETRIA_PADRAO)
    txt = " ".join(" ".join(notas.get(f"{s}_fig24", [])) for s in SOLOS)
    n = r"(\d+(?:,\d+)?)"
    m = re.search(r"H = " + n + " m", txt)
    if m:
        g["H"] = _num_pt(m.group(1))
    m = re.search(r"inclinação de " + n + "°", txt)
    if m:
        g["beta"] = _num_pt(m.group(1))
    m = re.search(r"crista em \(" + n + " m, " + n + r" m\)", txt)
    if m:
        g["crista"] = (_num_pt(m.group(1)), _num_pt(m.group(2)))
    m = re.search(r"pé em \(" + n + " m, " + n + r" m\)", txt)
    if m:
        g["pe"] = (_num_pt(m.group(1)), _num_pt(m.group(2)))
    m = re.search(r"domínio de " + n + " m por " + n + r" m \(" + n + " m²\\)", txt)
    if m:
        g["dom"] = (_num_pt(m.group(1)), _num_pt(m.group(2)))
        g["area"] = _num_pt(m.group(3))
    return g


def _metodo(est):
    """limiares da banda e da massa, reducao da carga de gravidade e resolucao do SRM, lidos de est"""
    md = dict(METODO_PADRAO)
    met = est.get("proveniencia", {}).get("metodos", {})
    n = r"(\d+(?:,\d+)?)"
    m = re.search(r"≥\s*" + n + r"\s*×", met.get("banda", ""))
    if m:
        md["banda"] = _num_pt(m.group(1))
    m = re.search(r"≥\s*" + n + r"\s*×", met.get("massa", ""))
    if m:
        md["massa"] = _num_pt(m.group(1))
    m = re.search(r"c'/" + n, met.get("FS", ""))
    if m:
        md["c_grav"] = _num_pt(m.group(1))
    for S in est.get("solos", {}).values():
        for c in S.get("casos", {}).values():
            h = (c.get("qui2_grade") or {}).get("h")
            if _ok(h):
                md["h"] = h
                return md
    m = re.search(r"resolução " + n, met.get("FS", ""))
    if m:
        md["h"] = _num_pt(m.group(1))
    return md


# ============================================================================================= grandezas de um caso
def caso_valor(S, c, fam):
    campo = (S.get("campo_familia") or CAMPO_PADRAO).get(fam, CAMPO_PADRAO[fam])
    return c.get(campo)


def x_txt(S, c, fam):
    return numg(ESCALA[fam] * caso_valor(S, c, fam), 4)


def ln(c, chave):
    return (c.get("lognormal") or {}).get(chave)


def moda_ln(c):
    mu, s = ln(c, "mu_ln"), ln(c, "sigma_ln")
    if not (_ok(mu) and _ok(s)):
        return None
    return math.exp(mu - s * s)


# metricas: (valor, texto, valor arredondado como aparece no texto)
METRICAS = {
    "Pf": (lambda c: c.get("Pf"), lambda v: pct(v, 2), lambda v: round(100 * v, 2)),
    "media": (lambda c: ln(c, "media"), lambda v: num(v, 3), lambda v: round(v, 3)),
    "cov": (lambda c: ln(c, "cov"), lambda v: pct(v, 1), lambda v: round(100 * v, 1)),
    "moda": (moda_ln, lambda v: num(v, 3), lambda v: round(v, 3)),
    "det": (lambda c: (c.get("det") or {}).get("FS"), None, lambda v: round(v, 7)),
}
NOME_METRICA = {"Pf": "$P_f$", "media": r"$\mu(\FS)$", "cov": r"$CoV(\FS)$", "moda": "moda",
                "det": r"$\FS$ determinístico"}


def serie(S, fam, metrica, ids=None, marca_ref=True):
    """'<metrica> = v1; v2 e v3 para <var> = x1; x2 e x3<unid>, respectivamente: <tendencia>' para os casos da
    familia (ou os ids dados), na ordem crescente da variavel"""
    casos = S["casos"]
    ids = list(ids if ids is not None else S["familias"].get(fam, []))
    f_val, f_txt, f_red = METRICAS[metrica]
    pares = []
    for i in ids:
        c = casos.get(i)
        if c is None:
            continue
        v = f_val(c)
        if not _ok(v):
            continue
        yt = fs_det_txt(c.get("det")) if metrica == "det" else f_txt(v)
        xt = x_txt(S, c, fam)
        if marca_ref and i == S.get("ref") and fam != "chuva":
            xt += " (referência)"
        pares.append((xt, yt, f_red(v)))
    if not pares:
        return None
    nome = NOME_METRICA[metrica]
    if len(pares) == 1:
        return f"{nome} = {pares[0][1]} para ${VAR[fam]} =$ {pares[0][0]}{UNID[fam]}"
    ys = [p[2] for p in pares]
    txt = (f"{nome} = {lista_num([p[1] for p in pares])} para ${VAR[fam]} =$ "
           f"{lista_num([p[0] for p in pares])}{UNID[fam]}, respectivamente")
    return txt + tendencia(ys, VAR[fam])


def tendencia(ys, var):
    """a variacao de uma sequencia de valores (ja arredondados como no texto) na ordem crescente da variavel"""
    pares = list(zip(ys, ys[1:]))
    if all(b == a for a, b in pares):
        return " (o mesmo valor em todos)"
    if all(b > a for a, b in pares):
        return f": aumenta a cada aumento de ${var}$"
    if all(b < a for a, b in pares):
        return f": diminui a cada aumento de ${var}$"
    if all(b >= a for a, b in pares):
        return f": não diminui com o aumento de ${var}$ (valores iguais em casos consecutivos)"
    if all(b <= a for a, b in pares):
        return f": não aumenta com o aumento de ${var}$ (valores iguais em casos consecutivos)"
    return " (variação não monotônica)"


def rotulo_caso(S, c):
    f = c.get("familia")
    if f == "ref":
        return "referência"
    if f in VAR:
        return f"${VAR[f]} = {mat(x_txt(S, c, f))}${UNID[f]}"
    return tex(c.get("id", "?"))


def fs_def_tex(c):
    inst = [horas(t) for t in c.get("instantes_h", [])]
    if c.get("fs_def") == "FS_24h" or len(inst) == 1:
        return r"o $\FS$ do SRM em $t = " + (mat(inst[0]) if inst else "--") + "$~h"
    return r"o menor $\FS$ do SRM entre os instantes " + lista_num(inst) + r"~h ($\FS$ mínimo)"


def tmin_tex(c):
    """'instante do FS minimo: 24 h em todas as N realizacoes' ou '... (numero de realizacoes): 24 h (995) e ...'"""
    cont = c.get("t_min_h_contagem") or []
    if not cont:
        return None
    if len(cont) == 1:
        return f"instante do $\\FS$ mínimo: {horas(cont[0][0])}~h em todas as {c.get('N')} realizações"
    return ("instante do $\\FS$ mínimo (número de realizações): "
            + lista_pt([f"{horas(t)}~h ({n})" for t, n in cont]))


def falhou_tex(c):
    reais = []
    for f in c.get("falhou") or []:
        reais += f.get("reais", [])
    if not reais:
        return None
    return f"{len(reais)} realizações marcadas em \\texttt{{falhou.txt}} (números {lista_num(sorted(reais))}), fora de $S$"


def faixa_tex(c):
    fx = c.get("faixa_reais")
    if not fx:
        return ""
    return f"{fx[0]} a {fx[1]}"


def versao_curta(v):
    return str(v).split("_")[0]


def cov_pf_resumo(S, ids):
    """frases sobre CoV(P_f) e o criterio de 5% do Ceron nos casos dados"""
    casos = S["casos"]
    abaixo, acima, zero = [], [], []
    for i in ids:
        c = casos[i]
        if not _ok(c.get("cov_Pf")):
            zero.append(rotulo_caso(S, c))
        elif c["cov_Pf"] < 0.05:
            abaixo.append(rotulo_caso(S, c))
        else:
            acima.append(rotulo_caso(S, c))
    n = len(ids)
    partes = []
    if len(abaixo) == n:
        partes.append(f"$CoV(P_f) < 5\\%$ em todos os {n} casos")
    elif len(acima) == n:
        partes.append(f"$CoV(P_f) \\geq 5\\%$ em todos os {n} casos")
    else:
        if abaixo:
            partes.append(f"$CoV(P_f) < 5\\%$ em {len(abaixo)} dos {n} casos")
        if acima:
            partes.append(f"$CoV(P_f) \\geq 5\\%$ em {len(acima)}" + ("" if abaixo else f" dos {n} casos")
                          + " (" + lista_pt(acima, "; ") + ")")
        if zero:
            partes.append(f"$P_f = 0$ e $CoV(P_f)$ indefinido em {len(zero)} (" + lista_pt(zero, "; ") + ")")
    return "; ".join(partes)


def s_resumo(S, ids):
    """'S = 300 em todos os casos' ou 'S = 1000 na referencia e 300 nos demais casos'"""
    casos = S["casos"]
    grupos = {}
    for i in ids:
        grupos.setdefault(casos[i].get("N"), []).append(i)
    if len(grupos) == 1:
        return f"$S = {next(iter(grupos))}$ em todos os casos"
    maior = max(grupos, key=lambda n: len(grupos[n]))
    outros = [f"$S = {n}$ " + ("na referência" if g == [S.get("ref")] else
                               "em " + lista_pt([rotulo_caso(S, casos[i]) for i in g], "; "))
              for n, g in grupos.items() if n != maior]
    return lista_pt(outros, "; ") + f" e {maior} nos demais casos"


def negativos_tex(S, ids):
    casos = S["casos"]
    itens = [f"{rotulo_caso(S, casos[i])} ({casos[i]['n_fs_negativo']})" for i in ids
             if (casos[i].get("n_fs_negativo") or 0) > 0]
    if not itens:
        return None
    return ("Realizações com $FS = -1$ (sem equilíbrio), contadas em $n_f$ e fora do ajuste log-normal: "
            + lista_pt(itens) + ".")


# ================================================================================================= blocos de LaTeX
def figura(solo, nn, notas_fig, subst=None, largura=r"\linewidth", extras=()):
    """um ambiente figure com o PDF de <saida>/figuras/<solo>_fig<nn>.pdf (ou um aviso, se ele nao existir), a
    legenda do Ceron traduzida, a figura analoga do Ceron e as notas do JSON"""
    arq = f"figuras/{solo}_fig{nn}.pdf"
    leg, curta = LEGENDAS[nn]
    if subst:
        leg = T(leg, **subst)
    notas = notas_fig.get(f"{solo}_fig{nn}") or []
    notas_txt = "".join(" " + frase(e) for e in extras if e)
    notas_txt += "".join(" " + frase(maiuscula(nota_tex(n))) for n in notas if str(n).strip())
    analoga, sufixo = ANALOGA.get(nn, (f"à Fig.~{nn}", ""))
    return T(r"""\numerofig{<<nn>>}
\begin{figure}[htbp]
\centering
\IfFileExists{<<arq>>}{\includegraphics[max width=<<larg>>,max height=0.8\textheight]{<<arq>>}}{\semarquivo{<<arqtt>>}}
\caption[<<curta>> (<<solo>>)]{<<leg>> (análoga <<analoga>> de \citealp{ceron2025}<<sufixo>>).<<notas>>}
\label{fig:<<solo>>_fig<<nn>>}
\end{figure}
""", nn=int(nn), arq=arq, larg=largura, arqtt=tex(arq), curta=curta, solo=solo, leg=leg, analoga=analoga,
             sufixo=sufixo, notas=notas_txt)


def tabela(nome, numero=None, legenda_falta="Tabela"):
    """\\input de <saida>/tabelas/<nome>.tex dentro de \\IfFileExists; sem o arquivo, um ambiente table com aviso,
    legenda e o mesmo rotulo"""
    arq = f"tabelas/{nome}.tex"
    rot = "tab:proveniencia" if nome == "tab_proveniencia" else f"tab:{nome}"
    pre = rf"\numerotab{{{numero}}}" + "\n" if numero is not None else ""
    return pre + T(r"""\IfFileExists{<<arq>>}{\input{<<arq>>}}{\begin{table}[htbp]
\centering
\caption{<<leg>>}\label{<<rot>>}
\semarquivo{<<arqtt>>}
\end{table}}
""", arq=arq, arqtt=tex(arq), leg=legenda_falta, rot=rot)


PREAMBULO = r"""\documentclass[a4paper,11pt]{article}
\usepackage[T1]{fontenc}
\usepackage{ae,aecompl}
\usepackage[brazil]{babel}
\usepackage[a4paper,top=2.5cm,bottom=2.5cm,left=2.5cm,right=2.5cm]{geometry}
\usepackage{amsmath,amssymb}
\usepackage{array}
\usepackage{booktabs}
\usepackage{graphicx}
\usepackage[export]{adjustbox}
\usepackage{caption}
\usepackage{subcaption}
\usepackage{siunitx}
\usepackage{xcolor}
\usepackage{enumitem}
\usepackage{placeins}
\usepackage[round]{natbib}
\usepackage[hyphens]{url}
\usepackage{hyperref}
\hypersetup{colorlinks=true,linkcolor=blue!50!black,citecolor=blue!50!black,urlcolor=blue!50!black,
  pdftitle={Análise probabilística de estabilidade do talude de Rahardjo sob chuva},
  pdfauthor={<<autor_pdf>>}}
\sisetup{output-decimal-marker={,},group-digits=false,per-mode=symbol}
\captionsetup{font=small,labelfont=bf,labelsep=colon}
\setlist{nosep,leftmargin=*}
\renewcommand{\labelitemi}{\ensuremath{\bullet}}% o \textbullet da fonte T1 ae viria em bitmap (TS1)
\setlength{\parskip}{0.4\baselineskip}
\setlength{\parindent}{0pt}
\renewcommand{\topfraction}{0.9}
\renewcommand{\bottomfraction}{0.8}
\renewcommand{\textfraction}{0.07}
\renewcommand{\floatpagefraction}{0.75}
\setcounter{topnumber}{3}
\setcounter{totalnumber}{4}
\urlstyle{same}
% sem hifenizacao nos nomes de pastas e arquivos (\texttt)
\AtBeginDocument{{\ttfamily\hyphenchar\font=-1}{\small\ttfamily\hyphenchar\font=-1}%
  {\footnotesize\ttfamily\hyphenchar\font=-1}{\scriptsize\ttfamily\hyphenchar\font=-1}}
\emergencystretch=2em
% caracteres do Unicode que podem vir nos arquivos de tabelas e notas e que a fonte T1 nao tem
\DeclareUnicodeCharacter{00AA}{\textsuperscript{a}}
\DeclareUnicodeCharacter{00BA}{\textsuperscript{o}}
\DeclareUnicodeCharacter{00B0}{\ensuremath{^\circ}}
\DeclareUnicodeCharacter{00B2}{\ensuremath{^2}}
\DeclareUnicodeCharacter{00B3}{\ensuremath{^3}}
\DeclareUnicodeCharacter{00B9}{\ensuremath{^1}}
\DeclareUnicodeCharacter{00B5}{\ensuremath{\mu}}
\DeclareUnicodeCharacter{00D7}{\ensuremath{\times}}
\DeclareUnicodeCharacter{00B1}{\ensuremath{\pm}}
\DeclareUnicodeCharacter{00B7}{\ensuremath{\cdot}}
\DeclareUnicodeCharacter{2212}{\ensuremath{-}}
\DeclareUnicodeCharacter{2264}{\ensuremath{\leq}}
\DeclareUnicodeCharacter{2265}{\ensuremath{\geq}}
\DeclareUnicodeCharacter{2248}{\ensuremath{\approx}}
\DeclareUnicodeCharacter{221A}{\ensuremath{\surd}}
\DeclareUnicodeCharacter{0394}{\ensuremath{\Delta}}
\DeclareUnicodeCharacter{0393}{\ensuremath{\Gamma}}
\DeclareUnicodeCharacter{03B1}{\ensuremath{\alpha}}
\DeclareUnicodeCharacter{03B2}{\ensuremath{\beta}}
\DeclareUnicodeCharacter{03B3}{\ensuremath{\gamma}}
\DeclareUnicodeCharacter{03BC}{\ensuremath{\mu}}
\DeclareUnicodeCharacter{03C1}{\ensuremath{\rho}}
\DeclareUnicodeCharacter{03C3}{\ensuremath{\sigma}}
\DeclareUnicodeCharacter{03C6}{\ensuremath{\varphi}}
\DeclareUnicodeCharacter{03C7}{\ensuremath{\chi}}
% arquivo de figura ou tabela ausente na compilacao: aviso visivel no lugar dele
\newcommand{\semarquivo}[1]{\fbox{\parbox[c][3.2cm][c]{0.8\linewidth}{\centering\small
  \textbf{Arquivo ausente na compilação:}\\[3pt]\texttt{#1}}}}
% numeracao das tabelas e figuras: o nome do solo e o numero da tabela ou figura analoga do Ceron (f50.25 = Fig. 25)
\newcommand{\numeracao}[1]{%
  \renewcommand{\thetable}{#1\arabic{table}}\renewcommand{\thefigure}{#1\arabic{figure}}%
  \renewcommand{\theHtable}{#1\arabic{table}}\renewcommand{\theHfigure}{#1\arabic{figure}}}
\newcommand{\numeracaoromana}{%
  \renewcommand{\thetable}{\Roman{table}}\renewcommand{\thefigure}{\Roman{figure}}%
  \renewcommand{\theHtable}{R\Roman{table}}\renewcommand{\theHfigure}{R\Roman{figure}}}
\newcommand{\numerofig}[1]{\setcounter{figure}{\numexpr#1-1\relax}}
\newcommand{\numerotab}[1]{\setcounter{table}{\numexpr#1-1\relax}}
\newcommand{\FS}{\mathit{FS}}
"""


# ================================================================================================ secoes do texto
def titulo(autor):
    linha_autor = rf"\author{{{autor}}}" if autor else r"\author{}"
    return T(r"""\title{Análise probabilística de estabilidade do talude de Rahardjo sob chuva: resultados no formato de
Vargas Ceron et al.\ (2025)}
<<autor>>
\date{\today}
\begin{document}
\maketitle
\tableofcontents
\numeracaoromana
""", autor=linha_autor)


def _campanhas(est):
    """resumo das campanhas usadas: {campanha: {"solos", "comb", "reais", "det", "inst"}}"""
    camp = {}
    for p in est.get("proveniencia", {}).get("pastas", []):
        if p.get("excluida"):
            continue
        if not p.get("uso"):
            continue
        d = camp.setdefault(p["campanha"], {"solos": set(), "comb": [], "reais": 0, "det": 0, "inst": set()})
        for u in p["uso"]:
            d["solos"].add(u.split("_")[0])
        if p.get("deterministica"):
            d["det"] += 1
        else:
            d["comb"].append(p)
            d["reais"] += int(p.get("concluidas", 0))
        d["inst"].add(tuple(p.get("instantes_h", [])))
    return camp


def _excluidas_tex(excl):
    """frase com as pastas nao usadas, agrupadas pelo motivo (de est)"""
    if not excl:
        return ""
    sufixo = ": sem análogo no Ceron"
    mot = {}
    for e in excl:
        m = str(e.get("motivo", ""))
        mot.setdefault(m[:-len(sufixo)] if m.endswith(sufixo) else m, []).append(e.get("pasta", ""))
    todos = all(str(e.get("motivo", "")).endswith(sufixo) for e in excl)
    partes = [f"{len(ps)} com {tex(m)}" for m, ps in mot.items()]
    return (f"Ficaram fora {len(excl)} pastas" + (", sem análogo no Ceron" if todos else "") + ": "
            + lista_pt(partes, "; ") + ". ")


def _rotulo_pasta(est, p):
    """o rotulo do caso que usa a pasta (o primeiro que nao e a referencia), para identificar contagens"""
    for u in p.get("uso", []):
        for S in est["solos"].values():
            c = S.get("casos", {}).get(u)
            if c is not None and u != S.get("ref"):
                return rotulo_caso(S, c)
    return tex(p.get("combinacao", ""))


def abertura(est, notas_fig):
    camp = _campanhas(est)
    solos = [s for s in SOLOS if s in est["solos"]]
    ks = {s: est["solos"][s].get("k_s_m_s") for s in solos}
    ks_txt = lista_pt([f"{s} ($k_s = {sci_m(ks[s])}$~m/s, Seção~\\ref{{sec:{s}}})" for s in solos])
    out = [r"\section{Conteúdo do documento}\label{sec:conteudo}"]
    out.append(T(r"""Este documento apresenta os resultados das campanhas de Monte Carlo do talude de \citet{rahardjo2007}
sob chuva, com campos aleatórios de coesão efetiva $c'$, ângulo de atrito efetivo $\varphi'$ e permeabilidade $k$, na
forma da Seção~6 de \citet{ceron2025} e, para a realização típica, da sua Seção~5.2: as mesmas tabelas, as mesmas
figuras e os mesmos eixos das figuras do artigo, ampliados só quando os dados não cabem neles
(Seção~\ref{sec:eixos}). O fator de estabilidade $\Gamma$ da análise limite do Ceron é
substituído pelo fator de segurança $\FS$ do método de redução de resistência (SRM). Há um capítulo para cada solo:
<<ks>>. O texto registra valores, contagens e variações que aparecem nos números, na ordem do texto do Ceron, sem
interpretação dos resultados; as definições estão na Seção~\ref{sec:definicoes}.""", ks=ks_txt))
    # campanhas e realizacoes
    itens = []
    tot_r, tot_d = 0, 0
    def _ordem_camp(k):
        d = camp[k]
        solo = min((SOLOS.index(s) if s in SOLOS else 9) for s in d["solos"])
        tem_ref = any(u.endswith("_ref") for p in d["comb"] for u in p.get("uso", []))
        return (solo, 0 if tem_ref else 1, k)

    def _ordem_comb(p):
        Is = [S["casos"][u].get("I_mm_h") for u in p.get("uso", []) for S in est["solos"].values()
              if u in S.get("casos", {})]
        return (min([v for v in Is if _ok(v)], default=0.0), p.get("pasta", ""))

    for nome in sorted(camp, key=_ordem_camp):
        d = camp[nome]
        d["comb"].sort(key=_ordem_comb)
        tot_r += d["reais"]
        tot_d += d["det"]
        ncomb = len(d["comb"])
        ns = sorted({int(p.get("concluidas", 0)) for p in d["comb"]})
        if len(ns) == 1:
            reais_txt = (f"{ns[0]} realizações cada" if ncomb > 1 else f"{ns[0]} realizações")
        else:
            reais_txt = (lista_pt([f"{int(p.get('concluidas', 0))} ({_rotulo_pasta(est, p)})" for p in d["comb"]],
                                  "; ") + " realizações")
        inst = sorted(d["inst"], key=len)
        inst_txt = lista_num([horas(t) for t in inst[-1]]) + "~h" if inst else "--"
        det_txt = (f"{d['det']} rodada{'s' if d['det'] > 1 else ''} determinística{'s' if d['det'] > 1 else ''}"
                   if d["det"] else "sem rodadas determinísticas")
        itens.append(T(r"\item \texttt{<<nome>>} (solo <<solos>>): <<ncomb>> combinaç<<plural>> usada<<s>>, <<reais>> "
                       r"(<<tot>> no total); <<det>>; $\FS$ calculado em <<inst>>.",
                       nome=tex(nome), solos=lista_pt(sorted(d["solos"], key=lambda s: SOLOS.index(s)
                                                             if s in SOLOS else 9)),
                       ncomb=ncomb, plural="ões" if ncomb > 1 else "ão", s="s" if ncomb > 1 else "",
                       reais=reais_txt, tot=d["reais"], det=det_txt, inst=inst_txt))
    excl = est.get("proveniencia", {}).get("excluidas", [])
    out.append(T(r"""Campanhas usadas (realizações concluídas, sem as marcadas em \texttt{falhou.txt}):
\begin{itemize}
<<itens>>
\end{itemize}
No total, <<tot_r>> realizações de Monte Carlo e <<tot_d>> rodadas determinísticas (com os valores médios das
propriedades, sem campos aleatórios). <<excl>>As pastas de cada caso, as versões do programa e as realizações que
falharam estão na Tabela~\ref{tab:proveniencia}.""", itens="\n".join(itens), tot_r=tot_r, tot_d=tot_d,
                 excl=_excluidas_tex(excl)))
    # numeracao e correspondencia
    ex = "f50" if "f50" in solos else (solos[0] if solos else "f50")
    out.append(T(r"""As tabelas e figuras de cada solo têm o número da tabela ou da figura análoga do Ceron, precedido do
nome do solo: a Figura~\ref{fig:<<ex>>_fig25}, por exemplo, é a análoga da Fig.~25 do Ceron para o solo <<ex>>. A
Tabela~\ref{tab:correspondencia} dá a correspondência completa; as tabelas sem análogo no Ceron têm números romanos
(Tabelas~\ref{tab:correspondencia} a \ref{tab:proveniencia}) ou a letra A (apêndice).""", ex=ex))
    linhas = []

    def cel(solo, tabs, figs):
        partes = []
        if tabs:
            partes.append(("Tab.~" if len(tabs) == 1 else "Tabs.~") +
                          lista_pt([rf"\ref{{tab:{solo}_{t}}}" for t in tabs]))
        if figs:
            if len(figs) > 2:
                partes.append(rf"Figs.~\ref{{fig:{solo}_fig{figs[0]}}} a \ref{{fig:{solo}_fig{figs[-1]}}}")
            else:
                partes.append(("Fig.~" if len(figs) == 1 else "Figs.~") +
                              lista_pt([rf"\ref{{fig:{solo}_fig{f}}}" for f in figs]))
        return "; ".join(partes)

    corr = [
        ("Tabela~2", "Dados estatísticos do modelo de referência", ["tab2"], []),
        ("Fig.~24", "Geometria e superfície crítica da análise determinística", [], ["24"]),
        ("Fig.~12", r"Realização típica dos campos de $k_v$, $c'$ e $\varphi'$", [], ["12"]),
        ("Fig.~13", "Excesso de poropressão da realização típica (aqui: poropressão $p_w$)", [], ["13"]),
        ("Fig.~14", "Superfície crítica da realização típica", [], ["14"]),
        ("Figs.~16 e 20", "Superfícies críticas de todas as realizações e modos de ruptura", [], ["16"]),
        ("Fig.~25", "Configuração de referência: densidade do fator e convergência de $P_f$", [], ["25"]),
        ("Fig.~26", "Correlação de Spearman parcial de segunda ordem (Seção~6.1)", [], ["26"]),
        ("Tabela~3, Figs.~27 a 30", "Coeficiente de variação das propriedades (Seção~6.2)", ["tab3"],
         ["27", "28", "29", "30"]),
        ("Tabela~4, Fig.~31", "Distâncias de autocorrelação (Seção~6.3)", ["tab4"], ["31"]),
        ("Tabela~5, Fig.~32", "Anisotropia da permeabilidade (Seção~6.4)", ["tab5"], ["32"]),
        ("Tabela~6, Fig.~33", "Razão de rebaixamento $h_w/H$ (Seção~6.5; aqui: intensidade da chuva $I$)", ["tab6"],
         ["33"]),
        ("--", r"Estatísticas amostrais, teste qui-quadrado e correlações de Spearman", ["tabA1", "tabA2", "tabA3"],
         []),
    ]
    for ceron, conteudo, tabs, figs in corr:
        celulas = [cel(s, tabs, figs) for s in solos]
        linhas.append(" & ".join([ceron, conteudo] + celulas) + r" \\")
    linhas.append(" & ".join(["--", "Proveniência dos dados"] +
                             [r"\multicolumn{" + str(len(solos)) + r"}{l}{Tab.~\ref{tab:proveniencia}}"]) + r" \\")
    cab = " & ".join(["Vargas Ceron et al.", "Conteúdo"] + [f"Solo {s}" for s in solos]) + r" \\"
    out.append(T(r"""\numerotab{1}
\begin{table}[htbp]
\centering
\caption{Correspondência entre as tabelas e figuras de \citet{ceron2025} e as deste documento.}
\label{tab:correspondencia}
\small
\begin{tabular}{>{\raggedright\arraybackslash}p{2.3cm}>{\raggedright\arraybackslash}p{5.3cm}<<colsolos>>}
\toprule
<<cab>>
\midrule
<<linhas>>
\bottomrule
\end{tabular}
\end{table}
""", cab=cab, linhas="\n".join(linhas),
                 colsolos="".join([r">{\raggedright\arraybackslash}p{" + ("3.25" if len(solos) > 1 else "6.5")
                                   + "cm}" for _ in solos])))
    return "\n\n".join(out)


def definicoes(est, notas_fig, md, geo):
    solos = [s for s in SOLOS if s in est["solos"]]
    cst = est.get("proveniencia", {}).get("constantes", {})
    ref = est["solos"][solos[0]]["casos"][est["solos"][solos[0]]["ref"]] if solos else {}
    out = [r"\section{Definições}\label{sec:definicoes}"]
    # modelo
    ks_txt = lista_pt([f"$k_s = {sci_m(est['solos'][s]['k_s_m_s'])}$~m/s ({s})" for s in solos])
    out.append(r"\subsection{Modelo e campos aleatórios}")
    out.append(T(r"""Talude de \citet{rahardjo2007} com altura $H = <<H>>$~m e inclinação $\beta = <<beta>>^\circ$ (crista
em $(<<cx>>; <<cy>>)$~m e pé em $(<<px>>; <<py>>)$~m), no domínio de <<dx>>~m $\times$ <<dy>>~m (área de solo de
<<area>>~m$^2$), sob chuva de intensidade $I$ durante <<tch>>~h. O fluxo não saturado de água e ar é acoplado à
mecânica (variante i); no transiente, $c'$ e $\tan\varphi'$ são divididos por $F_{TRANS} = <<ftr>>$, e o SRM usa a
resistência real. Valores médios: $c' = <<c>>$~kPa, $\varphi' = <<phi>>^\circ$, $\gamma = <<gam>>$~kN/m$^3$,
porosidade <<por>>, $E = <<E>>$~MPa, $\nu = <<nu>>$ e permeabilidade saturada <<ks>>.""",
                 H=mat(numg(geo["H"])), beta=mat(numg(geo["beta"])), cx=mat(numg(geo["crista"][0])),
                 cy=mat(numg(geo["crista"][1])), px=mat(numg(geo["pe"][0])), py=mat(numg(geo["pe"][1])),
                 dx=numg(geo["dom"][0]), dy=numg(geo["dom"][1]), area=numg(geo["area"]),
                 tch=numg(cst.get("chuva_duracao_h")), ftr=mat(numg(ref.get("f_trans"))),
                 c=mat(numg(cst.get("c_medio_kPa"))), phi=mat(numg(cst.get("phi_medio_graus"))),
                 gam=mat(numg(cst.get("gamma_kN_m3"))), por=numg(cst.get("porosidade")),
                 E=mat(numg(cst.get("E_MPa"))), nu=mat(numg(cst.get("nu"))), ks=ks_txt))
    Lr = cst.get("L_ref_m") or [ref.get("Lx"), ref.get("Ly")]
    out.append(T(r"""Os campos de $c'$, $\varphi'$ e $k_v$ (o multiplicador $k/k_s$ do tensor de permeabilidade, com
$k_h = \alpha\,k_v$) vêm da expansão de Karhunen-Loève com $M = <<M>>$ termos e covariância exponencial (eq.~64 do
Ceron), com distribuição log-normal (eqs.~65 e 66), os três campos independentes. A estatística de referência é a da
Tabela~2 do Ceron, com as médias de Rahardjo: $CoV(c') = <<cc>>\%$, $CoV(\varphi') = <<cp>>\%$,
$CoV(k_v) = <<ck>>\%$, $L_x = <<Lx>>$~m, $L_y = <<Ly>>$~m e $\alpha = k_h/k_v = <<alfa>>$. A realização $n$ usa os
mesmos números aleatórios em todos os casos e em todos os solos (semente <<sem>>).""",
                 M=ref.get("KL_M"), cc=mat(numg(100 * ref.get("cov_c", float("nan")), 3)),
                 cp=mat(numg(100 * ref.get("cov_phi", float("nan")), 3)),
                 ck=mat(numg(100 * ref.get("cov_k", float("nan")), 3)), Lx=mat(numg(Lr[0])), Ly=mat(numg(Lr[1])),
                 alfa=mat(numg(ref.get("alfa"))), sem=ref.get("semente")))
    # fator de seguranca
    fora, tot, vistos = 0, 0, set()
    for s in solos:
        for c in est["solos"][s]["casos"].values():
            k = tuple(c.get("pastas", []))
            if k in vistos:
                continue
            vistos.add(k)
            fora += c.get("n_fs_fora_grade") or 0
            tot += (c.get("amostra") or {}).get("n") or 0
    out.append(r"\subsection{Fator de segurança}")
    out.append(T(r"""No Ceron, a estabilidade é medida pelo fator de estabilidade $\Gamma$ da análise limite (abordagem
cinemática). Aqui, o fator de segurança $\FS$ vem do método de redução de resistência (SRM): $c'$ e $\tan\varphi'$ de
cada ponto de integração são divididos por $F$, e $\FS$ é o maior $F$ com equilíbrio, obtido por bissecção com
resolução <<h>>. $\FS = -1$ é o valor gravado quando não há equilíbrio nem com $c'/<<cg>>$ e
$\tan\varphi'/<<cg>>$. Pela eq.~(57)
do Ceron, $\Gamma \geq 1 \Leftrightarrow \FS \geq 1$; nas tabelas e figuras, $\FS$ ocupa o lugar de $\Gamma$. Dos
<<tot>> valores de $\FS$ positivos usados neste documento, <<fora>> não são múltiplos de <<h>>.""",
                 h=num(md["h"], 5), cg=mat(numg(md["c_grav"])), tot=tot, fora=fora))
    # conjuntos (instantes e definicao do FS)
    grupos = {}
    for s in solos:
        S = est["solos"][s]
        for cid, c in S["casos"].items():
            camps = sorted({p.split("/")[0] for p in c.get("pastas", [])})
            chave = (s, c.get("fs_def"), tuple(c.get("instantes_h", [])), tuple(camps))
            grupos.setdefault(chave, []).append(cid)
    linhas = []
    for (s, fsd, inst, camps), ids in grupos.items():
        S = est["solos"][s]
        fams = []
        for f in ("ref", "covk", "covc", "covphi", "s", "alfa", "chuva"):
            if any(S["casos"][i].get("familia") == f for i in ids):
                fams.append({"ref": "referência", "covk": "$CoV(k_v)$", "covc": "$CoV(c')$",
                             "covphi": r"$CoV(\varphi')$", "s": "$s$", "alfa": r"$\alpha$",
                             "chuva": "$I$"}[f])
        Ns = []
        for i in ids:
            if S["casos"][i].get("N") not in [n for n, _ in Ns]:
                Ns.append((S["casos"][i].get("N"), i))
        if len(Ns) == 1:
            S_txt = str(Ns[0][0])
        else:
            S_txt = lista_pt([f"{n} ({rotulo_caso(S, S['casos'][i])})" for n, i in Ns], "; ")
        if fsd == "FS_24h" or len(inst) == 1:
            fs_txt = r"$\FS$ em $t = " + mat(horas(inst[0] if inst else 24)) + "$~h"
        else:
            fs_txt = r"$\FS$ mínimo"
        linhas.append(" & ".join([s, lista_pt([r"\texttt{" + tex(c) + "}" for c in camps]), lista_pt(fams),
                                  lista_num([horas(t) for t in inst]), fs_txt, S_txt])
                      + r" \\")
    out.append(T(r"""O $\FS$ é calculado nos instantes da Tabela~\ref{tab:conjuntos}. O $\FS$ de cada realização é o
$\FS$ desse instante, quando há um só, ou o menor $\FS$ entre os instantes calculados ($\FS$ mínimo; em caso de
empate, o primeiro instante).
\numerotab{2}
\begin{table}[htbp]
\centering
\caption{Conjuntos de dados: instantes de cálculo do $\FS$ e definição do $\FS$ de cada realização.}
\label{tab:conjuntos}
\small
\begin{tabular}{l>{\raggedright\arraybackslash}p{2.4cm}>{\raggedright\arraybackslash}p{3.2cm}%
>{\raggedright\arraybackslash}p{2.4cm}>{\raggedright\arraybackslash}p{1.7cm}>{\raggedright\arraybackslash}p{2.8cm}}
\toprule
Solo & Campanhas & Casos & Instantes (h) & $\FS$ usado & $S$ \\
\midrule
<<linhas>>
\bottomrule
\end{tabular}
\end{table}""", linhas="\n".join(linhas)))
    # probabilidade de falha
    out.append(r"\subsection{Probabilidade de falha}")
    out.append(r"""Como nas eqs.~(67) e (69) do Ceron,
\begin{equation*}
P_f = P(\FS < 1) = P(\Gamma < 1) \approx \frac{n_f}{S},
\end{equation*}
com $n_f$ o número de realizações com $\FS < 1$ ($\FS = -1$ incluído) e $S$ o número de realizações usadas;
$\FS = 1$ exatamente não conta como falha. O coeficiente de variação de $P_f$ é o da eq.~(68):
\begin{equation*}
CoV(P_f) \approx \sqrt{\frac{1 - P_f}{S\,P_f}};
\end{equation*}
com $P_f = 0$, $CoV(P_f)$ é indefinido (``--'' nas tabelas). O Ceron aumenta $S$ até atingir o critério de
convergência $CoV(P_f) < 5\%$; aqui, $S$ é o número de realizações de cada campanha. A convergência de $P_f$ (Figura~25b de cada solo) é $P_f(n)$, a fração de $\FS < 1$ nas $n$ primeiras
realizações, em ordem crescente do número da realização.""")
    # ajuste log-normal e qui-quadrado
    out.append(r"\subsection{Distribuição ajustada e teste qui-quadrado}")
    tabsA = {n: lista_pt([rf"Tabela~\ref{{tab:{s}_tab{n}}}" for s in solos]) for n in ("A1", "A2")}
    out.append(T(r"""Como no Ceron, a distribuição do $\FS$ de cada caso é ajustada por uma log-normal de dois parâmetros,
por máxima verossimilhança, sobre as realizações com $\FS > 0$: $\mu_{\ln}$ é a média de $\ln \FS$ e $\sigma_{\ln}$ o
desvio-padrão populacional de $\ln \FS$. Da distribuição ajustada vêm a média $\mu = \exp(\mu_{\ln} +
\sigma_{\ln}^2/2)$, a mediana $\exp(\mu_{\ln})$, o coeficiente de variação $CoV = \sqrt{\exp(\sigma_{\ln}^2) - 1}$, o
desvio-padrão $\sigma = \mu\,CoV$ e a moda $\exp(\mu_{\ln} - \sigma_{\ln}^2)$. Os valores $\mu$, $\sigma$ e $CoV$ das
tabelas análogas às Tabelas~3 a 6 do Ceron são os da distribuição ajustada, e as funções densidade de probabilidade
das figuras são as dessa distribuição. As estatísticas amostrais estão no apêndice (<<tA1>>).

O Ceron escolhe a log-normal pelo teste qui-quadrado de aderência contra a normal. Aqui, as duas distribuições
ajustadas (a normal com $\mu$ e $\sigma$ de máxima verossimilhança) são comparadas pelo teste qui-quadrado com
$k = \lceil 2\,n^{0{,}4} \rceil$ classes e $k - 3$ graus de liberdade, $n$ o número de realizações com $\FS > 0$, em
duas variantes (<<tA2>>):
\begin{enumerate}[label=(\alph*)]
\item $k$ classes de igual probabilidade sob a distribuição ajustada $F$: cada $\FS$ vai para a classe
$\lfloor k\,F(\FS) \rfloor$, com frequência esperada $n/k$;
\item as bordas das classes são os pontos médios da grade do $\FS$ mais próximos dos quantis $j/k$ da distribuição
ajustada, e a frequência esperada de cada classe é $n$ vezes a probabilidade ajustada dela; nenhum valor da grade cai
sobre uma borda.
\end{enumerate}
O texto informa, para a configuração de referência, qual distribuição tem o menor $\chi^2$ em cada variante.""",
                 tA1=tabsA["A1"], tA2=tabsA["A2"]))
    # Spearman
    out.append(r"\subsection{Correlação de postos de Spearman}")
    out.append(T(r"""Como em \citet{hu2022} e na Seção~6.1 do Ceron, cada campo aleatório de uma realização é
transformado num valor escalar por média espacial, no instante do $\FS$: $c$ é a coesão $c'$ média na banda de
ruptura, $\varphi$ o ângulo de atrito $\varphi'$ médio na banda e $k_v$ a permeabilidade média ($k/k_s$) na massa que
desliza. A banda de ruptura são os pontos de integração com $\sqrt{J_2}$ da deformação plástica $\geq <<b>>$ vez o
máximo, no último estado convergido do SRM; a massa que desliza são os pontos com $|\Delta u| \geq <<m>>$ vez o
máximo, $\Delta u$ o incremento de deslocamento entre os dois últimos estados convergidos do SRM. O coeficiente de
correlação de postos de Spearman \citep{conover1999} é o coeficiente de correlação dos postos (postos médios nos
empates). O coeficiente parcial de segunda ordem $\rho_{X\,\FS,\,YZ}$ entre $X$ e o $\FS$, controlando as outras duas
variáveis $Y$ e $Z$ (eqs.~70 e 71 do Ceron), é calculado da inversa $P$ da matriz $4\times4$ de correlação de postos,
$\rho_{X\,\FS,\,YZ} = -P_{X\FS}/\sqrt{P_{XX}P_{\FS\FS}}$, o que equivale a aplicar as eqs.~70 e 71. O intervalo de 95\%
vem da transformação de Fisher, com erro-padrão $1/\sqrt{n - 3 - q}$ e $q = 2$ variáveis controladas.\footnote{No
texto de \citet{ceron2025}, Hu et al.\ são citados com o número [81] e Conover com o número [82]. Na lista de
referências do artigo, o item 81 é Wang, Meng e Qin (2024), o item 82 é Phoon e Kulhawy (1999), Conover é o item 84 e
Hu, Takahashi e Kasama (2022) é o item 85, o único de Hu et al.\ na lista.}""",
                 b=mat(numg(md["banda"])), m=mat(numg(md["massa"]))))
    # eixos e estilo
    out.append(r"\subsection{Eixos e estilo das figuras}\label{sec:eixos}")
    out.append(r"""Cada painel tem os limites, o passo das marcações, a escala (linear ou logarítmica) e os rótulos do
painel análogo do Ceron, lidos nas figuras do artigo; os rótulos trazem $\FS$ no lugar de $\Gamma$. Quando os dados
não cabem nesses limites, o eixo é estendido até a união da faixa do Ceron com a dos dados, arredondada para fora no
passo das marcações (ou num múltiplo dele, para limitar o número de marcações), e a legenda da figura diz qual eixo
foi ampliado e qual é a faixa do Ceron; o Apêndice~\ref{ap:eixos} lista todas as ampliações. O mesmo vale para as
barras de cores dos campos (análogas às Figs.~12 a 14). Nas famílias de curvas de densidade, cada valor tem o estilo de
linha que o Ceron usa para o mesmo valor; os valores que o Ceron não tem recebem os estilos seguintes. O caso de
referência é marcado por uma estrela. As linhas tracejadas horizontais são as do Ceron: $CoV(k_v) = 0$ nas figuras
análogas às Figs.~28a a 30a e $s = 400$ nas análogas às Figs.~31b a 31d. As marcações dos eixos usam vírgula decimal.

Formatos dos números, como nas tabelas do Ceron: $\mu$ e $\sigma$ com 3 casas decimais, $CoV$ com 1, $P_f$ e
$CoV(P_f)$ com 2; o $\FS$ determinístico com todos os seus dígitos.""")
    return "\n\n".join(out)


# --------------------------------------------------------------------------------------------- capitulo de um solo
def capitulo(solo, est, notas_fig, modos, md, geo):
    S = est["solos"][solo]
    casos = S["casos"]
    ref = casos[S["ref"]]
    ks = S.get("k_s_m_s")
    ks_ceron = S.get("k_s_ceron")
    out = []
    out.append(T(r"""\clearpage
\section{Solo <<solo>>}\label{sec:<<solo>>}
\numeracao{<<solo>>.}""", solo=solo))
    ks_mmh = ks * 3.6e6 if _ok(ks) else None
    out.append(T(r"""Solo <<solo>> de \citet{rahardjo2007}: permeabilidade saturada $k_s = <<ks>>$~m/s (<<ksmmh>>~mm/h),
$<<ksc>>$~m$^4$~kN$^{-1}$~s$^{-1}$ na unidade do Ceron ($k/\gamma_w$). $\FS$ da configuração de referência:
<<fsdef>>.""", solo=solo, ks=sci_m(ks), ksmmh=numg(ks_mmh, 3), ksc=sci_m(ks_ceron, 2), fsdef=fs_def_tex(ref)))

    # ---------------------------------------------------------------- 6 (inicio): Tabela 2 e Fig. 24
    out.append(T(r"""\subsection{Dados de referência e análise determinística (Ceron, Seção~6)}\label{sec:<<solo>>_ref_det}
A geometria do talude é fixa: altura $H = <<H>>$~m e inclinação $\beta = <<beta>>^\circ$ (Figura~\ref{fig:<<solo>>_fig24}).
Os dados estatísticos do modelo de referência estão na Tabela~\ref{tab:<<solo>>_tab2} (análoga à Tabela~2 do Ceron).
Na configuração de referência, o talude está sob chuva de $I = <<I>>$~mm/h durante <<tch>>~h, com permeabilidade
isotrópica ($\alpha = k_h/k_v = <<alfa>>$).""", solo=solo, H=mat(numg(geo["H"])), beta=mat(numg(geo["beta"])),
                 I=mat(numg(ref.get("I_mm_h"))),
                 tch=numg(est.get("proveniencia", {}).get("constantes", {}).get("chuva_duracao_h")),
                 alfa=mat(numg(ref.get("alfa")))))
    out.append(tabela(f"{solo}_tab2", 2, "Dados estatísticos do modelo de referência (arquivo ausente)."))
    det = ref.get("det")
    if det:
        serie_det = det.get("serie") or []
        serie_txt = lista_num([f"{txt_digitos(s[2]) if len(s) > 2 else numg(s[1])} ({horas(s[0])}~h)"
                               for s in serie_det])
        txt = T(r"""A Figura~\ref{fig:<<solo>>_fig24} mostra a superfície crítica (banda de ruptura do SRM) e o $\FS$ da
análise determinística, com os valores médios das propriedades: $\FS = <<fs>>$, em $t = <<t>>$~h (pasta <<pasta>>).
$\FS$ determinístico em cada instante calculado: <<serie>>.""", solo=solo, fs=mat(fs_det_txt(det)),
                t=mat(horas(det.get("t_h"))), pasta=pasta_tt(det.get("pasta", "")), serie=serie_txt)
        vdet = versao_curta(det.get("versao", ""))
        vmc = sorted({versao_curta(v) for v in (ref.get("versoes") or {})})
        if vmc and (len(vmc) > 1 or vmc[0] != vdet):
            txt += T(r""" A rodada determinística é da versão <<vd>> do programa; as <<N>> realizações do caso de referência
(<<pastas>>), da versão <<vm>>.""", vd=tex(vdet), N=ref.get("N"),
                     pastas=lista_pt([pasta_tt(p) for p in ref.get("pastas", [])]), vm=tex(lista_pt(vmc)))
        out.append(txt)
    else:
        out.append(T(r"""Não há rodada determinística da configuração de referência do solo <<solo>>.""", solo=solo))
    out.append(figura(solo, "24", notas_fig))
    out.append(r"\FloatBarrier")

    # ---------------------------------------------------------------- 5.2: realizacao tipica (Figs. 12 a 14)
    am = ref.get("amostras") or {}
    reais = am.get("real") or []
    r0 = reais[0] if reais else None
    fs0 = (am.get("FS") or [None])[0]
    t0 = (am.get("t_min_h") or [None])[0]
    kl = ref.get("KL_variancia")
    out.append(T(r"""\subsection{Uma realização típica do Monte Carlo (Ceron, Seção~5.2)}\label{sec:<<solo>>_tipica}
Como na Seção~5.2 do Ceron, a Figura~\ref{fig:<<solo>>_fig12} mostra uma realização dos campos aleatórios, a de número
<<r0>> do caso de referência ($L_x = <<Lx>>$~m e $L_y = <<Ly>>$~m): (a) permeabilidade vertical, (b) coesão e (c)
ângulo de atrito. A expansão de Karhunen-Loève tem $M = <<M>>$ termos e retém <<kl>> da variância.""", solo=solo,
                 r0=r0, Lx=mat(numg(ref.get("Lx"))), Ly=mat(numg(ref.get("Ly"))), M=ref.get("KL_M"),
                 kl=(num(100 * kl, 2) + r"\%") if _ok(kl) else "--"))
    out.append(figura(solo, "12", notas_fig, subst={"Lx": mat(numg(ref.get("Lx"))), "Ly": mat(numg(ref.get("Ly")))}))
    out.append(T(r"""Com essa realização da permeabilidade vertical e $\alpha = <<alfa>>$, o problema de fluxo de água e ar
acoplado à mecânica é resolvido pelo método dos elementos finitos. A Figura~\ref{fig:<<solo>>_fig13} mostra a
poropressão $p_w$ em $t = <<t0>>$~h, o instante do $\FS$ dessa realização. A Fig.~13 do Ceron mostra o excesso de
poropressão do fluxo permanente: as duas grandezas são diferentes.""", solo=solo, alfa=mat(numg(ref.get("alfa"))),
                 t0=mat(horas(t0))))
    out.append(figura(solo, "13", notas_fig))
    out.append(T(r"""A Figura~\ref{fig:<<solo>>_fig14} mostra a superfície crítica dessa realização (banda de ruptura do SRM
em $t = <<t0>>$~h), com as distribuições de $p_w$, $c'$ e $\varphi'$: $\FS = <<fs0>>$.""", solo=solo,
                 t0=mat(horas(t0)), fs0=mat(numg(fs0))))
    out.append(figura(solo, "14", notas_fig))
    out.append(r"\FloatBarrier")

    # ---------------------------------------------------------------- 6 (referencia): Fig. 25
    lnr = ref.get("lognormal") or {}
    q = ref.get("qui2") or {}
    qg = ref.get("qui2_grade") or {}
    out.append(T(r"""\subsection{Configuração de referência (Ceron, Seção~6)}\label{sec:<<solo>>_config_ref}
O Monte Carlo da configuração de referência tem $S = <<N>>$ realizações (números <<faixa>><<fusao>>), com
$n_f = <<nf>>$, $P_f = <<pf>>$ e $CoV(P_f) = <<cpf>>$ (<<crit>>). <<igual1>>A Figura~\ref{fig:<<solo>>_fig25}a mostra
o histograma do $\FS$ e a função densidade de probabilidade da distribuição log-normal ajustada, com o $\FS$
determinístico (linha tracejada vertical); a Figura~\ref{fig:<<solo>>_fig25}b mostra a convergência de $P_f$ com o
número de simulações.<<fnota>>""", solo=solo, N=ref.get("N"), faixa=faixa_tex(ref), fusao=_fusao_tex(ref),
                 fnota=_fusao_nota(ref), nf=ref.get("n_f"),
                 pf=mat(pct(ref.get("Pf"), 2)), cpf=mat(pct(ref.get("cov_Pf"), 2)),
                 crit=("abaixo de 5\\%" if _ok(ref.get("cov_Pf")) and ref["cov_Pf"] < 0.05 else
                       ("$P_f = 0$" if not _ok(ref.get("cov_Pf")) else "acima de 5\\%")),
                 igual1=(f"Em {ref['n_fs_igual_1']} realizações, $\\FS = 1$ exatamente. "
                         if (ref.get("n_fs_igual_1") or 0) > 0 else "")))
    out.append(figura(solo, "25", notas_fig))
    tm = tmin_tex(ref)
    if tm and ref.get("fs_def") == "FS_min":
        out.append(frase(maiuscula(tm)))
    out.append(_frase_qui2(solo, q, qg))
    out.append(T(r"""Como no Ceron, as funções densidade das figuras seguintes são as da distribuição log-normal ajustada.
Pelo ajuste, a distribuição do $\FS$ tem média <<media>>, mediana <<mediana>> e coeficiente de variação <<cov>>.""",
                 media=num(lnr.get("media"), 3), mediana=num(lnr.get("mediana"), 3), cov=pct(lnr.get("cov"), 2)))
    out.append(r"\FloatBarrier")

    # ---------------------------------------------------------------- Figs. 16 e 20: superficies criticas
    out.append(_secao_superficies(solo, S, notas_fig, modos, md))
    out.append(r"\FloatBarrier")

    # ---------------------------------------------------------------- 6.1: Spearman
    out.append(_secao_spearman(solo, S, notas_fig))
    out.append(r"\FloatBarrier")

    # ---------------------------------------------------------------- 6.2: CoV
    out.append(_secao_cov(solo, S, notas_fig))
    out.append(r"\FloatBarrier")

    # ---------------------------------------------------------------- 6.3: distancias de autocorrelacao
    out.append(_secao_s(solo, S, notas_fig))
    out.append(r"\FloatBarrier")

    # ---------------------------------------------------------------- 6.4: anisotropia
    out.append(_secao_alfa(solo, S, notas_fig))
    out.append(r"\FloatBarrier")

    # ---------------------------------------------------------------- 6.5 (analogo): intensidade da chuva
    out.append(_secao_chuva(solo, S, est, notas_fig))
    out.append(r"\FloatBarrier")
    return "\n\n".join(out)


def _fusao_tex(c):
    fu = c.get("fusao")
    if not fu:
        return ""
    partes = [f"{fx[0]} a {fx[1]} de {pasta_tt(p.split('/')[0])}" for p, fx in zip(fu.get("pastas", []),
                                                                                   fu.get("faixas", []))]
    return ": " + lista_pt(partes)


def _fusao_nota(c):
    fu = c.get("fusao")
    if fu and fu.get("parametros_identicos") and fu.get("reais_disjuntos"):
        return r" As pastas juntadas têm \texttt{parametros.txt} idêntico e números de realização disjuntos."
    return ""


def _frase_qui2(solo, q, qg):
    def par(d):
        no, lo = d.get("normal") or {}, d.get("lognormal") or {}
        return no, lo

    def menor(no, lo):
        if not (_ok(no.get("chi2")) and _ok(lo.get("chi2"))):
            return None
        if lo["chi2"] < no["chi2"]:
            return "o da log-normal"
        if no["chi2"] < lo["chi2"]:
            return "o da normal"
        return "o mesmo nas duas"

    no, lo = par(q)
    m = menor(no, lo)
    if m is None:
        return ""
    txt = T(r"""Teste qui-quadrado de aderência (Tabela~\ref{tab:<<solo>>_tabA2}), variante (a): $\chi^2 = <<cn>>$ para a
normal e $\chi^2 = <<cl>>$ para a log-normal ($k = <<k>>$ classes, <<gl>> graus de liberdade; $p = <<pn>>$ e
$p = <<pl>>$): o menor $\chi^2$ é <<m>>.""", solo=solo, cn=mat(num(no.get("chi2"), 1)), cl=mat(num(lo.get("chi2"), 1)),
            k=q.get("k"), gl=no.get("gl"), pn=num_p(no.get("p")), pl=num_p(lo.get("p")), m=m)
    gno, glo = par(qg)
    mg = menor(gno, glo)
    if mg is not None:
        txt += T(r""" Variante (b), com as bordas das classes na grade do $\FS$: $\chi^2 = <<cn>>$ (normal) e
$\chi^2 = <<cl>>$ (log-normal); o menor é <<m>>.""", cn=mat(num(gno.get("chi2"), 1)), cl=mat(num(glo.get("chi2"), 1)),
                 m=mg)
    return txt


def _secao_superficies(solo, S, notas_fig, modos, md):
    ref = S["casos"][S["ref"]]
    out = [T(r"""\subsection{Superfícies críticas das realizações (Ceron, Figs.~16 e 20)}\label{sec:<<solo>>_superficies}""", solo=solo)]
    m = (modos or {}).get(solo)
    if m:
        reais = m.get("realizacoes") or []
        faixa = f" (números {reais[0]} a {reais[-1]})" if len(reais) == 2 else ""
        Nd = m.get("N")
        Nc = m.get("N_caso") or ref.get("N")
        base = (f"das {Nd} realizações{faixa} da configuração de referência"
                if Nd == Nc else f"de {Nd} das {Nc} realizações{faixa} da configuração de referência")
        det = m.get("det") or {}
        nome_modo = {"acima": "acima do pé", "pe": "pelo pé", "abaixo": "abaixo do pé"}
        txt = T(r"""Como nas Figs.~16 e 20 do Ceron, a Figura~\ref{fig:<<solo>>_fig16} mostra as superfícies críticas
(bandas de ruptura) <<base>>, com a banda da análise determinística, e a porcentagem de cada modo de ruptura: acima do
pé do talude <<a>>, pelo pé <<p>> e abaixo do pé <<b>> (<<na>>, <<np>> e <<nb>> realizações). O modo é dado pelo ponto
mais baixo da banda, $y_{\min}$, em relação à cota do pé ($y = <<ype>>$~m), com tolerância $\delta = <<d>>$~m: acima
do pé se $y_{\min} > <<ype>> + \delta$, abaixo do pé se $y_{\min} < <<ype>> - \delta$ e pelo pé nos demais casos.""",
                solo=solo, base=base, a=pct(m.get("acima", float("nan")) / 100, 2),
                p=pct(m.get("pe", float("nan")) / 100, 2), b=pct(m.get("abaixo", float("nan")) / 100, 2),
                na=m.get("n_acima", "--"), np=m.get("n_pe", "--"), nb=m.get("n_abaixo", "--"),
                ype=mat(numg(m.get("y_pe_m"))), d=mat(numg(m.get("delta_m"))))
        if det:
            txt += T(r""" Banda determinística: $y_{\min} = <<y>>$~m, <<modo>> ($\FS = <<fs>>$, $t = <<t>>$~h).""",
                     y=mat(num(det.get("banda_ymin"), 2)), modo=nome_modo.get(det.get("modo"), tex(det.get("modo"))),
                     fs=mat(numg(det.get("FS"))), t=mat(horas(det.get("t_h"))))
        cr = ((modos or {}).get(f"{solo}_casos") or {}).get(S["ref"])
        if cr and cr.get("N") != Nd:
            txt += T(r""" Com o $y_{\min}$ da banda de \texttt{resumo.csv}, nas <<N>> realizações da configuração de
referência: acima do pé <<a>>, pelo pé <<p>> e abaixo do pé <<b>> (<<na>>, <<np>> e <<nb>> realizações).""",
                     N=cr.get("N"), a=pct(cr.get("acima", float("nan")) / 100, 2),
                     p=pct(cr.get("pe", float("nan")) / 100, 2), b=pct(cr.get("abaixo", float("nan")) / 100, 2),
                     na=cr.get("n_acima"), np=cr.get("n_pe"), nb=cr.get("n_abaixo"))
        out.append(txt)
    else:
        out.append(T(r"""A Figura~\ref{fig:<<solo>>_fig16} mostra as superfícies críticas (bandas de ruptura) das realizações
da configuração de referência, com a banda da análise determinística (as porcentagens dos modos de ruptura não estavam
disponíveis na geração deste documento: \texttt{figuras/modos\_ruptura.json} ausente).""", solo=solo))
    out.append(figura(solo, "16", notas_fig))
    return "\n\n".join(out)


def _secao_spearman(solo, S, notas_fig):
    ref = S["casos"][S["ref"]]
    sp = ref.get("spearman") or {}
    out = [T(r"""\subsection{Análise de sensibilidade (Ceron, Seção~6.1)}\label{sec:<<solo>>_sensibilidade}""", solo=solo)]
    so, ic = sp.get("segunda_ordem") or {}, sp.get("ic95") or {}
    if not so:
        out.append("Não há correlações de Spearman para a configuração de referência.")
        out.append(figura(solo, "26", notas_fig))
        return "\n\n".join(out)
    simb = {"c": (r"\rho_{c\,\FS,\,\varphi k_v}", "$c'$"), "phi": (r"\rho_{\varphi\,\FS,\,c k_v}", r"$\varphi'$"),
            "k": (r"\rho_{k_v\FS,\,c\varphi}", "$k_v$")}
    itens = []
    for v in ("c", "phi", "k"):
        r = so.get(v)
        if not _ok(r):
            continue
        i95 = ic.get(v) or [None, None]
        itens.append(T(r"""$<<s>> = <<r>>$ (intervalo de 95\%: $[<<lo>>; <<hi>>]$)""", s=simb[v][0],
                       r=mat(num(r, 3)), lo=mat(num(i95[0], 3)), hi=mat(num(i95[1], 3))))
    pos = [simb[v][1] for v in ("c", "phi", "k") if _ok(so.get(v)) and so[v] > 0]
    neg = [simb[v][1] for v in ("c", "phi", "k") if _ok(so.get(v)) and so[v] < 0]
    ordem = sorted([v for v in ("c", "phi", "k") if _ok(so.get(v))], key=lambda v: -abs(so[v]))
    sinais = []
    if pos:
        sinais.append(("o coeficiente de " if len(pos) == 1 else "os coeficientes de ") + lista_pt(pos) +
                      (" é positivo" if len(pos) == 1 else " são positivos"))
    if neg:
        sinais.append(("o" if len(neg) == 1 else "os") + ("" if pos else (" coeficiente" if len(neg) == 1 else
                                                                        " coeficientes")) + " de " +
                      lista_pt(neg) + (" é negativo" if len(neg) == 1 else " são negativos"))
    out.append(T(r"""Como no Ceron, a análise de sensibilidade usa a configuração de referência, com permeabilidade
isotrópica ($\alpha = <<alfa>>$), e as variáveis equivalentes de cada realização da Seção~\ref{sec:definicoes}:
$c'$ e $\varphi'$ médios na banda de ruptura e $k_v$ médio na massa que desliza. A Figura~\ref{fig:<<solo>>_fig26}
mostra os coeficientes de correlação de postos de Spearman parciais de segunda ordem entre essas variáveis e o $\FS$
($n = <<n>>$ realizações):
\begin{itemize}
<<itens>>
\end{itemize}
<<sinais>>. Em valor absoluto, a ordem é <<ordem>>. As correlações simples e os seus intervalos estão na
Tabela~\ref{tab:<<solo>>_tabA3}.""", alfa=mat(numg(ref.get("alfa"))), solo=solo, n=sp.get("n"),
                 itens="\n".join(r"\item " + i + ";" for i in itens[:-1]) + "\n" + r"\item " + itens[-1] + ".",
                 sinais=maiuscula(lista_pt(sinais)) if sinais else "",
                 ordem=lista_pt([f"{simb[v][1]} (${mat(num(abs(so[v]), 3))}$)" for v in ordem])))
    out.append(figura(solo, "26", notas_fig))
    return "\n\n".join(out)


def _valores_ceron(fam):
    if eixos_ceron is None:
        return None
    return list(getattr(eixos_ceron, "VALORES_CERON", {}).get(fam, [])) or None


def _comparacao_valores(S, fam, nome_var, escala=1.0, unid="", so_aqui=True):
    """frase com os valores do Ceron (eixos_ceron.VALORES_CERON) que nao foram rodados aqui e, se so_aqui, os daqui que
    nao estao entre eles"""
    vc = _valores_ceron(fam)
    if not vc:
        return ""
    aqui = sorted({round(caso_valor(S, S["casos"][i], fam), 9) for i in S["familias"].get(fam, [])})
    ceron = sorted({round(v, 9) for v in vc})
    if aqui == ceron:
        return " Os valores são os mesmos do Ceron."
    so_ceron = [v for v in ceron if v not in aqui]
    extra = [v for v in aqui if v not in ceron] if so_aqui else []
    partes = []
    if so_ceron:
        partes.append(f"o Ceron tem também ${nome_var} =$ {lista_num([numg(escala * v, 4) for v in so_ceron])}{unid}, "
                      "não rodado" + ("s" if len(so_ceron) > 1 else "") + " aqui")
    if extra:
        partes.append(f"${nome_var} =$ {lista_num([numg(escala * v, 4) for v in extra])}{unid} não "
                      + ("estão" if len(extra) > 1 else "está") + " no Ceron")
    return (" " + frase(maiuscula(lista_pt(partes, "; ")))) if partes else ""


def _lista_series(S, fam_metricas):
    itens = []
    for fam, metrica in fam_metricas:
        s = serie(S, fam, metrica)
        if s:
            itens.append(r"\item " + maiuscula(s))
    if not itens:
        return ""
    return "\\begin{itemize}\n" + ";\n".join(itens) + ".\n\\end{itemize}"


def _ja_nas_notas(notas_fig, chave, padrao):
    return any(re.search(padrao, str(n)) for n in notas_fig.get(chave, []))


def _secao_cov(solo, S, notas_fig):
    casos = S["casos"]
    fams = S["familias"]
    ref = casos[S["ref"]]
    ids = []
    for f in ("covk", "covc", "covphi"):
        for i in fams.get(f, []):
            if i not in ids:
                ids.append(i)
    if S["ref"] in ids:  # referencia primeiro, como na Tabela 3 do Ceron
        ids.remove(S["ref"])
        ids.insert(0, S["ref"])

    def faixa(f):
        vals = [ESCALA[f] * caso_valor(S, casos[i], f) for i in fams.get(f, [])]
        if not vals:
            return "--"
        return (f"de {numg(min(vals), 4)} a {numg(max(vals), 4)}\\% "
                f"({numg(ESCALA[f] * caso_valor(S, ref, f), 4)}\\% na referência)")

    out = [T(r"""\subsection{Coeficiente de variação das propriedades do solo (Ceron, Seção~6.2)}\label{sec:<<solo>>_cov}""",
             solo=solo)]
    comps = [_comparacao_valores(S, f, VAR[f], 100.0, r"\%") for f in ("covk", "covc", "covphi")]
    comp = (" Os valores são os mesmos do Ceron." if all(c == " Os valores são os mesmos do Ceron." for c in comps)
            else "".join(c for c in comps if c != " Os valores são os mesmos do Ceron."))
    out.append(T(r"""Como na Seção~6.2 do Ceron, o coeficiente de variação de cada propriedade varia numa faixa:
$CoV(k_v)$ <<fk>>, $CoV(c')$ <<fc>> e $CoV(\varphi')$ <<fp>>; os demais parâmetros ficam nos valores de
referência.<<comp>> A Tabela~\ref{tab:<<solo>>_tab3} (análoga à Tabela~3 do Ceron) resume os resultados: <<sres>>;
<<covpf>>.""", fk=faixa("covk"), fc=faixa("covc"), fp=faixa("covphi"), comp=comp, solo=solo, sres=s_resumo(S, ids),
                 covpf=cov_pf_resumo(S, ids)))
    neg = negativos_tex(S, ids)
    if neg:
        out.append(neg)
    out.append(tabela(f"{solo}_tab3", 3, "Efeito do coeficiente de variação (arquivo ausente)."))
    out.append(T(r"""A Figura~\ref{fig:<<solo>>_fig27} mostra a função densidade de probabilidade do $\FS$ (log-normal
ajustada) para os valores de $CoV$ da permeabilidade vertical (a), da coesão (b) e do ângulo de atrito (c). Moda das
distribuições ajustadas (o coeficiente de variação $CoV(\FS)$ de cada uma está na Tabela~\ref{tab:<<solo>>_tab3} e na
Figura~\ref{fig:<<solo>>_fig30}):""", solo=solo))
    out.append(figura(solo, "27", notas_fig))
    out.append(_lista_series(S, [("covk", "moda"), ("covc", "moda"), ("covphi", "moda")]))
    out.append(T(r"""A Figura~\ref{fig:<<solo>>_fig28} mostra a probabilidade de falha $P_f$ em função do coeficiente de
variação de cada propriedade, com o caso de referência:""", solo=solo))
    out.append(figura(solo, "28", notas_fig))
    out.append(_lista_series(S, [("covk", "Pf"), ("covc", "Pf"), ("covphi", "Pf")]))
    k0 = [i for i in fams.get("covk", []) if _ok(caso_valor(S, casos[i], "covk"))
          and caso_valor(S, casos[i], "covk") == 0.0]
    if k0:
        c0 = casos[k0[0]]
        out.append(T(r"""A linha tracejada horizontal das Figuras~\ref{fig:<<solo>>_fig28}a, \ref{fig:<<solo>>_fig29}a e
\ref{fig:<<solo>>_fig30}a é o caso $CoV(k_v) = 0$, como no Ceron: $P_f = <<pf>>$, $\mu(\FS) = <<mu>>$ e
$CoV(\FS) = <<cv>>$.""", solo=solo, pf=mat(pct(c0.get("Pf"), 2)), mu=mat(num(ln(c0, "media"), 3)),
                     cv=mat(pct(ln(c0, "cov"), 1))))
    out.append(T(r"""As Figuras~\ref{fig:<<solo>>_fig29} e \ref{fig:<<solo>>_fig30} mostram a média $\mu(\FS)$ e o
coeficiente de variação $CoV(\FS)$ da distribuição ajustada em função do coeficiente de variação de cada
propriedade:""", solo=solo))
    out.append(figura(solo, "29", notas_fig))
    out.append(figura(solo, "30", notas_fig))
    out.append(_lista_series(S, [("covk", "media"), ("covc", "media"), ("covphi", "media"), ("covk", "cov"),
                                 ("covc", "cov"), ("covphi", "cov")]))
    return "\n\n".join(out)


def _secao_s(solo, S, notas_fig):
    casos = S["casos"]
    ids = S["familias"].get("s", [])
    ref = casos[S["ref"]]
    out = [T(r"""\subsection{Distâncias de autocorrelação (Ceron, Seção~6.3)}\label{sec:<<solo>>_s}""", solo=solo)]
    sval = lista_num([numg(caso_valor(S, casos[i], "s"), 4) for i in ids])
    # o maior s e o da linha tracejada do Ceron (Figs. 31b a 31d), fora da legenda da Fig. 31a (VALORES_CERON)
    comp = _comparacao_valores(S, "s", "s", so_aqui=False)
    klv = [(numg(caso_valor(S, casos[i], "s"), 4), casos[i].get("KL_variancia")) for i in ids]
    kl_txt = (" Variância retida pela expansão de Karhunen-Loève com $M = " + str(ref.get("KL_M")) + "$ termos: " +
              lista_num([num(100 * v, 3) + r"\%" for _, v in klv if _ok(v)]) + " para $s =$ " +
              lista_num([s for s, v in klv if _ok(v)]) + ", respectivamente.") if klv else ""
    out.append(T(r"""Como na Seção~6.3 do Ceron, as distâncias de autocorrelação são proporcionais às da configuração
de referência, $(L_x, L_y) = s\,(<<Lx>>~\mathrm{m}, <<Ly>>~\mathrm{m})$ (eq.~72 do Ceron), com o fator de escala
$s = <<sref>>$ na referência. Valores de $s$: <<sval>>.<<comp>><<kl>> A Tabela~\ref{tab:<<solo>>_tab4} (análoga à
Tabela~4 do Ceron) resume os resultados: <<sres>>; <<covpf>>.""", Lx=mat(numg(ref.get("Lx"))),
                 Ly=mat(numg(ref.get("Ly"))), sref=mat(numg(caso_valor(S, ref, "s"))), sval=sval, comp=comp,
                 kl=kl_txt, solo=solo, sres=s_resumo(S, ids), covpf=cov_pf_resumo(S, ids)))
    out.append(tabela(f"{solo}_tab4", 4, "Efeito das distâncias de autocorrelação (arquivo ausente)."))
    smax = max((caso_valor(S, casos[i], "s") for i in ids), default=None)
    out.append(T(r"""A Figura~\ref{fig:<<solo>>_fig31} mostra a função densidade de probabilidade do $\FS$ (a), a
probabilidade de falha (b), o coeficiente de variação do $\FS$ (c) e a média do $\FS$ (d) em função de $s$, com o caso
de referência; como no Ceron, o caso $s = <<smax>>$ é a linha tracejada horizontal de (b) a (d):""", solo=solo,
                 smax=mat(numg(smax, 4))))
    out.append(figura(solo, "31", notas_fig))
    out.append(_lista_series(S, [("s", "Pf"), ("s", "cov"), ("s", "media")]))
    return "\n\n".join(out)


def _secao_alfa(solo, S, notas_fig):
    casos = S["casos"]
    ids = S["familias"].get("alfa", [])
    out = [T(r"""\subsection{Anisotropia da permeabilidade (Ceron, Seção~6.4)}\label{sec:<<solo>>_alfa}""", solo=solo)]
    aval = lista_num([numg(caso_valor(S, casos[i], "alfa"), 4) for i in ids])
    comp = _comparacao_valores(S, "alfa", r"\alpha")
    com_det = [i for i in ids if casos[i].get("det")]
    sem_det = [i for i in ids if not casos[i].get("det")]
    if sem_det and com_det:
        det_txt = (" Há rodada determinística para " + lista_pt([rotulo_caso(S, casos[i]) for i in com_det], "; ") +
                   "; para " + lista_pt([rotulo_caso(S, casos[i]) for i in sem_det], "; ") +
                   ", não (``--'' na tabela).")
    elif sem_det:
        det_txt = " Não há rodadas determinísticas desta família (``--'' na tabela)."
    else:
        det_txt = ""
    out.append(T(r"""Como na Seção~6.4 do Ceron, o tensor de permeabilidade tem a componente vertical $k_v$ e a razão de
anisotropia $\alpha = k_h/k_v$; o campo de $k_v$ e os demais parâmetros ficam nos dados de referência. Valores de
$\alpha$: <<aval>>.<<comp>> A Tabela~\ref{tab:<<solo>>_tab5} (análoga à Tabela~5 do Ceron) resume os resultados, com
os da análise determinística (o $\FS$ da rodada com os valores médios, pela mesma definição).<<det>> <<sres>>;
<<covpf>>.""", aval=aval, comp=comp, solo=solo, det=det_txt, sres=maiuscula(s_resumo(S, ids)),
                 covpf=cov_pf_resumo(S, ids)))
    out.append(tabela(f"{solo}_tab5", 5, "Efeito da razão de anisotropia da permeabilidade (arquivo ausente)."))
    out.append(T(r"""A Figura~\ref{fig:<<solo>>_fig32} mostra a função densidade de probabilidade do $\FS$ (a), a
probabilidade de falha (b), o coeficiente de variação do $\FS$ (c) e a média do $\FS$ (d) em função de $\alpha$, com o
caso de referência:""", solo=solo))
    out.append(figura(solo, "32", notas_fig))
    series = [("alfa", "Pf"), ("alfa", "cov"), ("alfa", "media")]
    if len(com_det) >= 2:
        series.append(("alfa", "det"))
    out.append(_lista_series(S, series))
    return "\n\n".join(out)


def _secao_chuva(solo, S, est, notas_fig):
    casos = S["casos"]
    ids = S["familias"].get("chuva", [])
    ref = casos[S["ref"]]
    out = [T(r"""\subsection{Intensidade da chuva (no lugar da Seção~6.5 do Ceron)}\label{sec:<<solo>>_chuva}""",
             solo=solo)]
    if not ids:
        out.append("Não há casos de intensidade da chuva neste solo.")
        return "\n\n".join(out)
    hw = _valores_ceron("hw")
    hw_txt = (f"de {numg(min(hw), 3)} a {numg(max(hw), 3)} " if hw else "")
    ks = S.get("k_s_m_s")
    ks_mmh = ks * 3.6e6 if _ok(ks) else None
    ivals = [caso_valor(S, casos[i], "chuva") for i in ids]
    rel = (" (" + lista_num([numg(v / ks_mmh, 3) for v in ivals]) + " vezes $k_s$, respectivamente)"
           if _ok(ks_mmh) else "")
    c0 = casos[ids[0]]
    out.append(T(r"""Na Seção~6.5 do Ceron, a componente de carregamento é o rebaixamento: a razão $h_w/H$ varia
<<hw>>(Tabela~6 e Fig.~33 do Ceron). Aqui, a intensidade $I$ da chuva ocupa o lugar de $h_w/H$:
$I =$ <<ivals>>~mm/h<<rel>>, durante <<tch>>~h, com os demais dados de referência. $\FS$ de cada realização:
<<fsdef>>. Realizações de cada caso:""", hw=hw_txt, ivals=lista_num([numg(v, 4) for v in ivals]), rel=rel,
                 tch=numg(est.get("proveniencia", {}).get("constantes", {}).get("chuva_duracao_h")),
                 fsdef=fs_def_tex(c0)))
    itens = []
    fundidos = False
    for i in ids:
        c = casos[i]
        fundidos = fundidos or bool(c.get("fusao"))
        partes = [f"$S = {c.get('N')}$ (números {faixa_tex(c)}{_fusao_tex(c)})"]
        f = falhou_tex(c)
        if f:
            partes.append(f)
        for cc in c.get("fs_csv_corrompido") or []:
            partes.append(T(r"""na realização <<r>>, \texttt{fs.csv} tem <<n>> linhas corrompidas, e o $\FS$ de
<<inst>>~h vem de \texttt{resumo.csv}""", r=cc.get("real"), n=cc.get("linhas_corrompidas", "--"),
                            inst=lista_num([horas(t) for t in cc.get("instantes_de_resumo", [])])))
        tm = tmin_tex(c)
        if tm and c.get("fs_def") == "FS_min":
            partes.append(tm)
        itens.append(r"\item " + rotulo_caso(S, c) + ": " + "; ".join(partes))
    out.append("\\begin{itemize}\n" + ";\n".join(itens) + ".\n\\end{itemize}")
    if fundidos:
        out.append(r"As pastas juntadas num caso têm \texttt{parametros.txt} idêntico e números de realização "
                   r"disjuntos.")
    # a linha de I = 36 mm/h e a referencia
    rc = S.get("ref_chuva")
    extras = []
    if rc in casos:
        c36 = casos[rc]
        if S.get("mesmos_dados_ref_chuva"):
            out.append(T(r"""A linha $I = <<I>>$~mm/h da Tabela~\ref{tab:<<solo>>_tab6} é a configuração de referência das
seções anteriores (os mesmos dados, $S = <<N>>$).""", I=mat(numg(c36.get("I_mm_h"))), solo=solo, N=c36.get("N")))
            if not _ja_nas_notas(notas_fig, f"{solo}_fig33", r"mesmos dados|configuração de referência"):
                extras.append(T(r"""O caso $I = <<I>>$~mm/h é a configuração de referência (os mesmos dados)""",
                                I=mat(numg(c36.get("I_mm_h")))))
        else:
            cmp_ = est.get("proveniencia", {}).get("comparacao_f50_ref") if solo == "f50" else None
            txt = T(r"""A linha $I = <<I>>$~mm/h da Tabela~\ref{tab:<<solo>>_tab6} não é o mesmo conjunto de dados da
configuração de referência das seções anteriores: aqui, <<p36>>, $S = <<N36>>$, com <<fs36>>; na referência,
<<pref>>, $S = <<Nref>>$, com <<fsref>>.""", I=mat(numg(c36.get("I_mm_h"))), solo=solo,
                    p36=lista_pt([pasta_tt(p) for p in c36.get("pastas", [])]), N36=c36.get("N"),
                    fs36=fs_def_tex(c36), pref=lista_pt([pasta_tt(p) for p in ref.get("pastas", [])]),
                    Nref=ref.get("N"), fsref=fs_def_tex(ref))
            chaves_fs = [k for k in (cmp_ or {}).get("iguais", {}) if re.fullmatch(r"FS_[\d.]+h", k)]
            if chaves_fs:
                k = chaves_fs[0]
                fx = cmp_.get("faixa_comum") or []
                txt += T(r""" Nas realizações <<a>> a <<b>>, comuns aos dois conjuntos, o $\FS$ em $t = <<t>>$~h é
idêntico em <<ig>> de <<n>>.""", a=fx[0] if fx else "--", b=fx[1] if fx else "--", ig=cmp_["iguais"][k],
                         n=cmp_.get("reais_comuns"), t=mat(horas(float(k[3:-1]))))
            out.append(txt)
            if not _ja_nas_notas(notas_fig, f"{solo}_fig33", r"não é o mesmo conjunto"):
                extras.append(T(r"""O caso $I = <<I>>$~mm/h desta figura (<<p36>>, $S = <<N36>>$, <<fs36>>) não é o
mesmo conjunto de dados da configuração de referência das Figuras~\ref{fig:<<solo>>_fig25} a
\ref{fig:<<solo>>_fig32} (<<pref>>, $S = <<Nref>>$, <<fsref>>)""", I=mat(numg(c36.get("I_mm_h"))), solo=solo,
                                p36=lista_pt([pasta_tt(p) for p in c36.get("pastas", [])]), N36=c36.get("N"),
                                fs36=fs_def_tex(c36), pref=lista_pt([pasta_tt(p) for p in ref.get("pastas", [])]),
                                Nref=ref.get("N"), fsref=fs_def_tex(ref)))
    com_det = [i for i in ids if casos[i].get("det")]
    out.append(T(r"""A Tabela~\ref{tab:<<solo>>_tab6} (análoga à Tabela~6 do Ceron) resume os resultados, com os da
análise determinística (o $\FS$ da rodada com os valores médios, pela mesma definição); <<covpf>>.""",
                 solo=solo, covpf=cov_pf_resumo(S, ids)))
    out.append(tabela(f"{solo}_tab6", 6, "Intensidade da chuva (arquivo ausente)."))
    out.append(T(r"""A Figura~\ref{fig:<<solo>>_fig33} mostra a função densidade de probabilidade do $\FS$ (a; legenda em
$I$, mm/h), a probabilidade de falha (b), o coeficiente de variação do $\FS$ (c) e a média do $\FS$ (d) em função de
$I$, num eixo linear próprio, em mm/h:""", solo=solo))
    out.append(figura(solo, "33", notas_fig, extras=extras))
    series = [("chuva", "Pf"), ("chuva", "cov"), ("chuva", "media")]
    if len(com_det) >= 2:
        series.append(("chuva", "det"))
    out.append(_lista_series(S, series))
    return "\n\n".join(out)


# --------------------------------------------------------------------------------------------------------- apendices
def apendices(est, notas_fig, tem_estat, tem_campos):
    solos = [s for s in SOLOS if s in est["solos"]]
    out = [r"""\clearpage
\appendix
\section{Tabelas complementares}\label{ap:tabelas}
Estatísticas amostrais do $\FS$ por caso (Tabelas A1), teste qui-quadrado de aderência das distribuições normal e
log-normal ajustadas (Tabelas A2), correlações de postos de Spearman da configuração de referência (Tabelas A3) e
proveniência dos dados (Tabela~\ref{tab:proveniencia})."""]
    for s in solos:
        out.append(rf"\numeracao{{{s}.A}}")
        out.append(tabela(f"{s}_tabA1", 1, f"Estatísticas amostrais do $\\FS$ por caso, solo {s} (arquivo ausente)."))
        out.append(tabela(f"{s}_tabA2", 2, f"Teste qui-quadrado, solo {s} (arquivo ausente)."))
        out.append(tabela(f"{s}_tabA3", 3, f"Correlações de Spearman, solo {s} (arquivo ausente)."))
        out.append(r"\FloatBarrier")
    out.append(r"\numeracaoromana")
    out.append(tabela("tab_proveniencia", 3, "Proveniência dos dados (arquivo ausente)."))
    out.append(r"\FloatBarrier")
    # eixos ampliados
    out.append(r"""\section{Eixos ampliados}\label{ap:eixos}""")
    ordem = ["12", "13", "14", "16", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33"]
    itens = []
    for s in solos:
        for nn in ordem:
            ns = [nota_tex(n).rstrip(".") for n in (notas_fig.get(f"{s}_fig{nn}") or [])
                  if re.search(r"ampliad[oa]", str(n))]
            if ns:
                itens.append(rf"\item Figura~\ref{{fig:{s}_fig{nn}}} (análoga à Fig.~{int(nn)} do Ceron): "
                             + frase("; ".join(ns)))
    faltam = []
    if not tem_estat:
        faltam.append(r"\texttt{figuras/legendas\_estat.json} (Figuras 25 a 33)")
    if not tem_campos:
        faltam.append(r"\texttt{figuras/legendas\_campos.json} (Figuras 12 a 14, 16 e 24)")
    txt = (r"Eixos e barras de cores ampliados em relação aos do painel análogo do Ceron (regra na "
           r"Seção~\ref{sec:eixos}), como nas notas das legendas das figuras:")
    if itens:
        out.append(txt + "\n\\begin{itemize}\n" + "\n".join(itens) + "\n\\end{itemize}")
    else:
        out.append(r"Nenhuma nota de legenda registra eixo ou barra de cores ampliados.")
    if faltam:
        out.append("Notas de legenda ausentes na geração deste documento: " + lista_pt(faltam) + ".")
    return "\n\n".join(out)


def referencias():
    return r"""\begin{thebibliography}{4}
\bibitem[Conover(1999)]{conover1999}
CONOVER, W.~J. \emph{Practical Nonparametric Statistics}. 3.~ed. Wiley, 1999. (Wiley Series in Probability and
Statistics, Applied Probability and Statistics Section.)

\bibitem[Hu et~al.(2022)]{hu2022}
HU, L.; TAKAHASHI, A.; KASAMA, K. Effect of spatial variability on stability and failure mechanisms of 3D slope using
random limit equilibrium method. \emph{Soils and Foundations}, v.~62, 101225, 2022.
DOI: \href{https://doi.org/10.1016/j.sandf.2022.101225}{10.1016/j.sandf.2022.101225}.

\bibitem[Rahardjo et~al.(2007)]{rahardjo2007}
RAHARDJO, H.; ONG, T.~H.; REZAUR, R.~B.; LEONG, E.~C. Factors controlling instability of homogeneous soil slopes under
rainfall. \emph{Journal of Geotechnical and Geoenvironmental Engineering}, v.~133, n.~12, p.~1532--1543, 2007.
DOI: \href{https://doi.org/10.1061/(ASCE)1090-0241(2007)133:12(1532)}{10.1061/(ASCE)1090-0241(2007)133:12(1532)}.

\bibitem[Vargas Ceron et~al.(2025)]{ceron2025}
VARGAS CERON, M.; CECÍLIO, D.~L.; LINN, R.~V.; MAGHOUS, S. Stability analysis of slope subjected to seepage forces
considering spatial variability of soil properties. \emph{International Journal for Numerical and Analytical Methods in
Geomechanics}, v.~49, n.~11, p.~2459--2491, 2025.
DOI: \href{https://doi.org/10.1002/nag.3993}{10.1002/nag.3993}.
\end{thebibliography}
"""


# ============================================================================================================ autor
def _autor():
    """o nome do usuario, se aparecer na documentacao do repositorio (autor do artigo em docs/); senao, vazio"""
    raiz = os.path.dirname(os.path.dirname(_PAI))  # neopz-master
    arq = os.path.join(raiz, "docs", "biot_nao_saturado", "artigo", "artigo.tex")
    try:
        with open(arq, encoding="utf-8") as f:
            txt = f.read()
    except OSError:
        return ""
    m = re.search(r"\\author(?:\[[^\]]*\])?\{(Diogo[^\\}]*)", txt)
    if not m:
        return ""
    nome = m.group(1).strip()
    # a grafia do sobrenome no README do projeto (Cecílio, com acento)
    try:
        with open(os.path.join(_PAI, "README.md"), encoding="utf-8") as f:
            if "Cecílio" in f.read():
                nome = nome.replace("Cecilio", "Cecílio")
    except OSError:
        pass
    return nome


# ======================================================================================================== compilacao
def _erros_log(arq_log):
    erros, overfull = [], []
    try:
        with open(arq_log, encoding="latin-1") as f:
            linhas = f.read().splitlines()
    except OSError:
        return ["sem log"], []
    for i, l in enumerate(linhas):
        if l.startswith("!") or re.match(r"^\S+:\d+: ", l):
            erros.append(" | ".join(linhas[i:i + 3]))
        elif l.startswith("Overfull \\hbox"):
            m = re.search(r"\(([\d.]+)pt too wide\)", l)
            if m and float(m.group(1)) > 1.0:
                overfull.append(l.strip())
    return erros, overfull


def compila(saida, nome=NOME_TEX):
    """latexmk -pdf -interaction=nonstopmode em <saida>; sem latexmk, pdflatex tres vezes"""
    pdf = os.path.join(saida, nome + ".pdf")
    if shutil.which("latexmk"):
        cmd = ["latexmk", "-pdf", "-g", "-interaction=nonstopmode", "-file-line-error", nome + ".tex"]
        r = subprocess.run(cmd, cwd=saida, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
                           encoding="utf-8", errors="replace", timeout=1800)
        cod = r.returncode
        saida_txt = r.stdout
    else:
        cod, saida_txt = 0, ""
        for _ in range(3):
            r = subprocess.run(["pdflatex", "-interaction=nonstopmode", "-file-line-error", nome + ".tex"],
                               cwd=saida, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
                               encoding="utf-8", errors="replace", timeout=900)
            cod = cod or r.returncode
            saida_txt += r.stdout
    erros, overfull = _erros_log(os.path.join(saida, nome + ".log"))
    return {"pdf": pdf if os.path.isfile(pdf) else None, "codigo": cod, "erros": erros, "overfull": overfull,
            "saida": saida_txt}


# ============================================================================================================== gera
def monta(est, saida):
    """o texto LaTeX do relatorio"""
    AVISOS.clear()
    est = _carrega(est)
    fig_dir = os.path.join(saida, "figuras")
    ne = _le_json(os.path.join(fig_dir, "legendas_estat.json"))
    nc = _le_json(os.path.join(fig_dir, "legendas_campos.json"))
    notas_fig = {}
    for d in (nc, ne):
        if isinstance(d, dict):
            for k, v in d.items():
                notas_fig.setdefault(k, [])
                notas_fig[k] += [str(x) for x in (v if isinstance(v, list) else [v])]
    modos = _le_json(os.path.join(fig_dir, "modos_ruptura.json")) or {}
    md = _metodo(est)
    geo = _geometria(notas_fig)
    autor = _autor()
    partes = [T(PREAMBULO, autor_pdf=autor), titulo(autor), abertura(est, notas_fig), definicoes(est, notas_fig, md, geo)]
    for s in SOLOS:
        if s in est.get("solos", {}):
            partes.append(capitulo(s, est, notas_fig, modos, md, geo))
    partes.append(apendices(est, notas_fig, ne is not None, nc is not None))
    partes.append(referencias())
    partes.append(r"\end{document}")
    return "\n\n".join(partes) + "\n"


def gera(est, saida):
    """escreve <saida>/relatorio.tex, compila com latexmk e devolve o caminho do PDF"""
    saida = os.path.abspath(saida)
    os.makedirs(saida, exist_ok=True)
    texto = monta(est, saida)
    arq = os.path.join(saida, NOME_TEX + ".tex")
    tmp = arq + ".tmp"
    with open(tmp, "w", encoding="utf-8") as f:
        f.write(texto)
    os.replace(tmp, arq)
    for a in AVISOS:
        print("AVISO (documento):", a, file=sys.stderr)
    r = compila(saida)
    if r["pdf"] is None or r["codigo"] != 0:
        # uma segunda tentativa: um arquivo de figura gravado por outro processo durante a compilacao (lido pela
        # metade) da erro fatal no pdflatex; a segunda leitura encontra o arquivo inteiro
        time.sleep(2)
        r = compila(saida)
    if r["erros"]:
        print("documento: erros do LaTeX em " + NOME_TEX + ".log:", file=sys.stderr)
        for e in r["erros"][:30]:
            print("  " + e, file=sys.stderr)
    if r["overfull"]:
        print(f"documento: {len(r['overfull'])} linhas mais largas que a página (overfull > 1 pt):", file=sys.stderr)
        for e in r["overfull"][:20]:
            print("  " + e, file=sys.stderr)
    if r["pdf"] is None:
        raise RuntimeError("o PDF não foi gerado; fim da saída do latexmk:\n" + r["saida"][-4000:])
    if r["codigo"] != 0:
        print(f"documento: latexmk terminou com código {r['codigo']}", file=sys.stderr)
    return r["pdf"]


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("uso: python3 -I documento.py <saida>  (le <saida>/estatisticas.json, as tabelas e as figuras)")
        sys.exit(2)
    print(gera(os.path.join(sys.argv[1], "estatisticas.json"), sys.argv[1]))
