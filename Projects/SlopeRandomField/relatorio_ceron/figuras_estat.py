"""figuras_estat.py: as Figs. 25 a 33 de Vargas Ceron et al. (2025), "Stability Analysis of Slope Subjected to Seepage
Forces Considering Spatial Variability of Soil Properties", Int J Numer Anal Methods Geomech 49:2459-2491
(doi 10.1002/nag.3993), Secao 6, para cada solo (f50 e f100), a partir de estatisticas.json (dados.calcula):

  <solo>_fig25  caso de referencia: (a) histograma do FS (densidade), PDF log-normal ajustada e FS da solucao
                deterministica; (b) convergencia de P_f com o numero de simulacoes
  <solo>_fig26  correlacao de Spearman parcial de segunda ordem entre FS e (c, phi, k_v) do caso de referencia
  <solo>_fig27  PDFs log-normais ajustadas com o CoV de k_v (a), c (b) e phi (c)
  <solo>_fig28  P_f contra o CoV de k_v (a), c (b) e phi (c); em (a), CoV(k_v) = 0 e uma linha tracejada horizontal
  <solo>_fig29  media mu(FS) da log-normal ajustada contra o CoV (como a Fig. 28)
  <solo>_fig30  CoV(FS) da log-normal ajustada contra o CoV (como a Fig. 28)
  <solo>_fig31  distancias de autocorrelacao (L_x, L_y) = s (20 m, 2 m): PDFs (a), P_f (b), CoV(FS) (c) e mu(FS) (d);
                em (b) a (d), s = 400 e uma linha tracejada horizontal (o eixo x vai de 0 a 20, como no Ceron)
  <solo>_fig32  anisotropia alfa = k_h/k_v: PDFs (a), P_f (b), CoV(FS) (c) e mu(FS) (d)
  <solo>_fig33  intensidade da chuva I (no lugar da razao de rebaixamento h_w/H do Ceron): PDFs (a), P_f (b),
                CoV(FS) (c) e mu(FS) (d); o eixo x de (b) a (d) e I (mm/h)

mais as copias PNG (150 dpi) e legendas_estat.json = {"<solo>_figNN": [nota, ...]} (notas das legendas: texto simples,
sem caracteres especiais do LaTeX e so com caracteres Latin-1, para entrar na legenda com ou sem escape; as notas de
eixo ampliado vem por ultimo).

Eixos: os do Ceron (eixos_ceron.FIG); quando os dados desenhados nao cabem, eixos_ceron.amplia() estende o eixo a
uniao da faixa do Ceron com a dos dados e a nota "(painel) eixo <x|y> ampliado; no Ceron: <a> a <b>" vai para a
legenda. Dados de cada eixo: os pontos e linhas desenhados; nas PDFs, o maximo das curvas e a faixa de x em que a PDF
passa de 0,1% do seu pico. Se o passo devolvido por amplia() nao e redondo (0,25 x 2,5 = 0,625; 0,025 x 5 = 0,125;
0,02 x 2 = 0,04), usa-se o menor multiplo redondo (1, 2, 2,5 ou 5 x 10^n) do passo do Ceron com ate 11 marcacoes.

Estatisticas (as do Ceron, calculadas em dados.py): mu e CoV = media e coeficiente de variacao da log-normal ajustada
por maxima verossimilhanca (as das Tabelas 3 a 6); P_f = n_f/N, n_f = numero de realizacoes com FS < 1.

Desenho como o do Ceron: paineis (a) e (b) lado a lado, (c) centrado abaixo (Figs. 27 a 30) ou (c) e (d) abaixo
(Figs. 31 a 33), com a letra do painel abaixo do rotulo do eixo x; pontos pretos ligados por segmentos de reta, o caso
de referencia como estrela (asterisco); estilos das PDFs pelo valor (eixos_ceron.ESTILOS e VALORES_CERON); legendas
nas posicoes do Ceron (acima a direita) ou, se cobrirem dados, na primeira posicao livre. Todas as figuras tem 16 cm
de largura (a largura do texto do relatorio): incluidas com width=\\linewidth, as letras ficam no tamanho desenhado
(marcacoes 7,5 pt, rotulos 10 pt).

Uso: gera(est, saida) grava em <saida>/figuras; est e o dicionario de dados.calcula() ou o caminho de
estatisticas.json. Linha de comando:
    python3 -I figuras_estat.py <saida>            (le <saida>/estatisticas.json)"""

import json
import math
import os
import sys

_AQUI = os.path.dirname(os.path.abspath(__file__))
for _p in (os.path.dirname(_AQUI), _AQUI):  # o pacote e Projects/SlopeRandomField (o pipeline roda com python3 -I)
    if _p not in sys.path:
        sys.path.insert(0, _p)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import FixedLocator, FuncFormatter, LogLocator, NullFormatter
from matplotlib.transforms import blended_transform_factory

import eixos_ceron as EC

# ---------------------------------------------------------------------------------------------------------------------
# aparencia (a das figuras do Ceron: fontes com serifa, mathtext cm, eixos pretos, sem grade, marcacoes para fora)

CM = 1.0 / 2.54
LARG = 16.0  # largura de todas as figuras (cm)
AW = 5.9  # largura dos eixos de um painel (cm)
ESQ = 1.45  # a esquerda dos eixos: marcacoes e rotulo do eixo y (cm)
VAO = 2.0  # do eixo direito do painel da esquerda ao eixo esquerdo do painel da direita (cm)
TOPO = 0.25  # acima da primeira linha de paineis (cm)
BASE = 1.75  # abaixo dos eixos: marcacoes, rotulo do eixo x e "(a)" (cm)
BASE_SEM_LETRA = 1.3
ENTRE = 0.3  # entre linhas de paineis (cm)
ASPECTO = 1.82  # largura/altura dos eixos nas Figs. 27 a 33 do Ceron
ASPECTO_25 = 1.52  # Fig. 25
ASPECTO_26 = 1.51  # Fig. 26
AW_26 = 7.8  # largura dos eixos da Fig. 26 (painel unico)

LW_PDF = 1.15  # curvas das PDFs (Figs. 27, 31a, 32a, 33a)
LW_25 = 1.8  # curva e linha tracejada da Fig. 25a
LW_LINHA = 1.0  # linhas que ligam os pontos e linhas tracejadas
MS_PONTO = 6.0  # circulos pretos (no Ceron, cerca de 3,5% da largura dos eixos)
MS_REF = 6.5  # estrela (asterisco de 6 bracos) do caso de referencia
MEW_REF = 1.1
MARCA_REF = (6, 2, 0)
CINZA_HIST = "0.83"  # histograma (o cinza claro da Fig. 25a)
CINZA_SETA = "0.5"  # seta da Fig. 25b
MARGEM_X = 0.025  # margem do eixo x dos graficos de pontos alem da marcacao extrema, quando ha ponto nela
MARGEM_X_31 = 0.05  # Fig. 31 (b a d): o eixo do Ceron vai de 0 a 21
PICO_FRACAO = 1.0e-3  # faixa de x de uma PDF: onde ela passa de 0,1% do pico
MAX_MARCAS = 11
TAM_ROTULO_LINHA = 9.5  # rotulos das linhas tracejadas
TAM_LETRA = 10.5  # (a), (b), ...

ESTILO = {
    "font.family": "serif",
    "font.serif": ["LM Roman 10", "Latin Modern Roman", "DejaVu Serif"],  # o Computer Modern de 10 pt do Ceron
    "mathtext.fontset": "cm",
    "axes.unicode_minus": True,
    "font.size": 9.5,
    "axes.linewidth": 0.9,
    "axes.edgecolor": "black",
    "axes.labelcolor": "black",
    "axes.labelsize": 10.0,
    "axes.labelpad": 3.0,
    "axes.grid": False,
    "text.color": "black",
    "xtick.color": "black",
    "ytick.color": "black",
    "xtick.labelsize": 7.5,
    "ytick.labelsize": 7.5,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "xtick.major.size": 3.5,
    "ytick.major.size": 3.5,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "xtick.minor.size": 2.0,
    "xtick.minor.width": 0.6,
    "xtick.top": False,
    "ytick.right": False,
    "legend.fontsize": 9.0,
    "legend.title_fontsize": 9.0,
    "legend.handlelength": 2.0,  # o padrao do matplotlib, o das legendas do Ceron
    "legend.fancybox": True,
    "legend.framealpha": 0.8,
    "legend.edgecolor": "0.8",
    "lines.dash_capstyle": "butt",
}

REF_ROTULO = "Caso de referência"

# eixo x das Figs. 33b a 33d: a intensidade da chuva (sem correspondente no Ceron, cujo eixo e h_w/H de 0,5 a 1,0)
LIM_CHUVA = (0.0, 80.0, 10.0)
ROT_CHUVA = r"$I$ (mm/h)"

NOME_FAMILIA = {"covk": "CoV(kv)", "covc": "CoV(c)", "covphi": "CoV(phi)", "s": "s", "alfa": "alfa", "chuva": "I"}
UNID_FAMILIA = {"chuva": " mm/h"}
NOME_Y = {"Pf": "Pf", "media": "média", "cov": "CoV"}


# ---------------------------------------------------------------------------------------------------------------------
# numeros com virgula decimal


def _casas(v, maximo=6):
    """o numero de casas decimais de v (0,025 -> 3; 0,5 -> 1; 5 -> 0)"""
    for d in range(maximo + 1):
        w = v * 10 ** d
        if abs(w - round(w)) < 1e-7 * max(1.0, abs(w)):
            return d
    return maximo


def _virg(v, casas=None):
    """numero com virgula decimal; casas = None: o minimo de casas (ate 6)"""
    s = f"{v:.6f}".rstrip("0").rstrip(".") if casas is None else f"{v:.{casas}f}"
    if s.startswith("-") and float(s) == 0:
        s = s[1:]
    return s.replace(".", ",")


def _virg_marca(v, casas):
    return _virg(v, casas).replace("-", "\N{MINUS SIGN}")


def _pct(v, casas=None):
    """fracao -> porcentagem com virgula (0,6 -> '60'; 0,075 -> '7,5')"""
    return _virg(100.0 * v, casas)


def _y_txt(chave, v):
    """valor de P_f (4 casas, as 2 casas da porcentagem das tabelas), da media (3 casas) ou do CoV (3 casas)"""
    return _virg(v, 4 if chave == "Pf" else 3)


def _lista(itens):
    itens = list(itens)
    return itens[0] if len(itens) == 1 else ", ".join(itens[:-1]) + " e " + itens[-1]


def _lista_num(valores, unid=""):
    """'1, 2, 5 e 400'; com algum decimal, separado por ponto e virgula: '3,6; 9; 36 e 80'"""
    txt = [_virg(v) for v in valores]
    sep = "; " if any("," in t for t in txt) else ", "
    s = txt[0] if len(txt) == 1 else sep.join(txt[:-1]) + " e " + txt[-1]
    return s + unid


# ---------------------------------------------------------------------------------------------------------------------
# eixos do Ceron, ampliados quando os dados nao cabem


def _redondo(q):
    e = math.floor(math.log10(q) + 1e-12)
    m = q / 10 ** e
    return any(abs(m - r) < 1e-6 for r in (1.0, 2.0, 2.5, 5.0, 10.0))


def _amplia(lim, dados):
    """eixos_ceron.amplia(); se o passo devolvido nao e redondo (0,625, 0,125, 0,04, ...), o menor multiplo redondo
    do passo do Ceron com ate MAX_MARCAS marcacoes. Devolve (a, b, passo, ampliado, passo_corrigido)."""
    a, b, p, ampl = EC.amplia(lim, dados)
    if not ampl or p is None or _redondo(p):
        return a, b, p, ampl, False
    a0, b0, p0 = lim
    lo, hi = min(a0, dados[0]), max(b0, dados[1])
    q = None
    for m in (1, 2, 2.5, 4, 5, 10, 20, 25, 40, 50, 100, 200, 250, 400, 500, 1000):
        q = p0 * m
        if _redondo(q) and (hi - lo) / q + 1 <= MAX_MARCAS + 1e-9:
            break
    na = math.floor(lo / q + 1e-9) * q
    nb = math.ceil(hi / q - 1e-9) * q
    return min(na, a0), max(nb, b0), q, True, True


def _faixa_txt(lim):
    a, b, p = lim
    if p is None:
        return f"{_virg(a)} a {_virg(b)}"
    c = max(_casas(p), _casas(a), _casas(b))
    return f"{_virg(a, c)} a {_virg(b, c)}"


class Eixos:
    """aplica a um painel os eixos do Ceron (ou os ampliados) e guarda as notas de eixo ampliado"""

    def __init__(self, registro, fig_id):
        self.notas = []
        self.registro = registro  # [(fig_id, painel, eixo, ceron, novo, passo_corrigido)]
        self.fig_id = fig_id

    def aplica(self, ax, painel, eixo, lim, dados, margem=(0.0, 0.0), log=False):
        """eixo 'x' ou 'y'; lim = (min, max, passo) do Ceron; dados = (min, max) do que o painel desenha; margem =
        fracao da faixa somada alem da marcacao extrema (esquerda/baixo, direita/cima)"""
        a, b, p, ampl, corr = _amplia(lim, dados)
        if ampl:
            pref = f"({painel}) " if painel else ""
            self.notas.append(f"{pref}eixo {eixo} ampliado; no Ceron: {_faixa_txt(lim)}")
            self.registro.append((self.fig_id, painel, eixo, _faixa_txt(lim), _faixa_txt((a, b, p)), corr))
        eix = ax.xaxis if eixo == "x" else ax.yaxis
        if log:
            (ax.set_xscale if eixo == "x" else ax.set_yscale)("log")
            (ax.set_xlim if eixo == "x" else ax.set_ylim)(a, b)
            eix.set_major_locator(LogLocator(base=10.0, numticks=20))
            eix.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10), numticks=100))
            eix.set_minor_formatter(NullFormatter())
            return a, b, p
        # marcacoes nos multiplos do passo dentro da faixa (a ponta de baixo pode nao ter marcacao, como na Fig. 32d)
        k0, k1 = math.ceil(a / p - 1e-9), math.floor(b / p + 1e-9)
        marcas = [round(k * p, 10) for k in range(k0, k1 + 1)]
        c = _casas(p)
        fx = b - a
        (ax.set_xlim if eixo == "x" else ax.set_ylim)(a - margem[0] * fx, b + margem[1] * fx)
        eix.set_major_locator(FixedLocator(marcas))
        eix.set_major_formatter(FuncFormatter(lambda v, pos, c=c: _virg_marca(v, c)))
        return a, b, p


# ---------------------------------------------------------------------------------------------------------------------
# log-normal


def _pdf_lognormal(x, mu, sg):
    x = np.asarray(x, dtype=float)
    f = np.zeros_like(x)
    m = x > 0
    lx = np.log(x[m])
    f[m] = np.exp(-0.5 * ((lx - mu) / sg) ** 2) / (x[m] * sg * math.sqrt(2.0 * math.pi))
    return f


def _suporte_pdf(mu, sg):
    """(pico, x_min, x_max): o maximo da PDF e a faixa de x em que ela passa de PICO_FRACAO do pico"""
    moda = math.exp(mu - sg * sg)
    pico = float(_pdf_lognormal([moda], mu, sg)[0])
    t = np.linspace(mu - 12.0 * sg, mu + 12.0 * sg, 40001)
    x = np.exp(t)
    f = _pdf_lognormal(x, mu, sg)
    ok = np.nonzero(f >= PICO_FRACAO * pico)[0]
    return pico, float(x[ok[0]]), float(x[ok[-1]])


# ---------------------------------------------------------------------------------------------------------------------
# disposicao dos paineis (cm)


def _figura(linhas, aspecto, base=BASE, larg_eixos=AW):
    ah = larg_eixos / aspecto
    alt = TOPO + linhas * (ah + base) + (linhas - 1) * ENTRE
    fig = plt.figure(figsize=(LARG * CM, alt * CM))
    return fig, ah, alt


def _x_colunas(larg_eixos=AW):
    total = ESQ + larg_eixos + VAO + larg_eixos
    x0 = (LARG - total) / 2.0 + ESQ
    return [x0, x0 + larg_eixos + VAO]


def _x_centro(larg_eixos=AW):
    return (LARG - (ESQ + larg_eixos)) / 2.0 + ESQ


def _eixos(fig, alt, x0, linha, ah, base=BASE, larg_eixos=AW):
    """eixos na linha 'linha' (0 = a de cima), a partir de x0 (cm)"""
    y0 = alt - TOPO - ah - linha * (ah + base + ENTRE)
    return fig.add_axes([x0 / LARG, y0 / alt, larg_eixos / LARG, ah / alt])


def _rotulo_y(ax, texto, **kw):
    """rotulo do eixo y; se for mais comprido que 92% da altura dos eixos, em duas linhas (quebra no espaco que
    equilibra as linhas)"""
    lab = ax.set_ylabel(texto, **kw)
    r = ax.figure.canvas.get_renderer()
    if "$" in texto or " " not in texto:
        return lab
    if lab.get_window_extent(r).height <= 0.92 * ax.get_window_extent(r).height:
        return lab
    pal = texto.split(" ")
    i = min(range(1, len(pal)), key=lambda k: abs(len(" ".join(pal[:k])) - len(" ".join(pal[k:]))))
    return ax.set_ylabel(" ".join(pal[:i]) + "\n" + " ".join(pal[i:]), linespacing=1.05, **kw)


def _letras(fig, paineis):
    """'(a)', '(b)', ... abaixo do rotulo do eixo x de cada painel, como no Ceron"""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    larg_px, alt_px = fig.canvas.get_width_height()
    for ax, letra in paineis:
        bb_ax = ax.get_window_extent(r)
        bb_rot = ax.xaxis.label.get_window_extent(r)
        y = bb_rot.y0 - r.points_to_pixels(5.0)
        x = (bb_ax.x0 + bb_ax.x1) / 2.0 - 0.35 * CM * fig.dpi
        fig.text(x / larg_px, y / alt_px, f"({letra})", ha="center", va="top", fontsize=TAM_LETRA)


# ---------------------------------------------------------------------------------------------------------------------
# legendas e rotulos sem sobrepor os dados


def _obstaculos(ax, r, ignorar=()):
    """pontos (x, y, raio) em pixels das linhas e marcadores dos eixos e retangulos das barras"""
    pts = []
    for ln in ax.get_lines():
        if ln in ignorar or not ln.get_visible():
            continue
        xy = np.asarray(ln.get_xydata(), dtype=float)
        if xy.size == 0:
            continue
        d = ln.get_transform().transform(xy)
        d = d[np.all(np.isfinite(d), axis=1)]
        if len(d) == 0:
            continue
        if ln.get_linestyle() not in ("None", "none", "", " "):
            meia = r.points_to_pixels(ln.get_linewidth()) / 2.0
            for (x0, y0), (x1, y1) in zip(d[:-1], d[1:]):
                n = int(max(abs(x1 - x0), abs(y1 - y0)) / 2.0) + 1
                t = np.linspace(0.0, 1.0, n + 1)
                pts.extend(zip(x0 + t * (x1 - x0), y0 + t * (y1 - y0), np.full(t.size, meia)))
        if ln.get_marker() not in (None, "None", "none", "", " "):
            raio = r.points_to_pixels(ln.get_markersize()) / 2.0
            pts.extend((xx, yy, raio) for xx, yy in d)
    caixas = [p.get_window_extent(r) for p in ax.patches if p.get_visible()]
    return np.array(pts, dtype=float).reshape(-1, 3), caixas


def _conta_colisoes(bb, pts, caixas, folga):
    n = 0
    if len(pts):
        x, y, rr = pts[:, 0], pts[:, 1], pts[:, 2] + folga
        n += int(np.count_nonzero((x + rr > bb.x0) & (x - rr < bb.x1) & (y + rr > bb.y0) & (y - rr < bb.y1)))
    for c in caixas:
        if c.width <= 0 or c.height <= 0:
            continue
        if c.x1 > bb.x0 - folga and c.x0 < bb.x1 + folga and c.y1 > bb.y0 - folga and c.y0 < bb.y1 + folga:
            n += 100
    return n


def _legenda(ax, handles, rotulos, opcoes, **kw):
    """opcoes = [loc ou (loc, {kw})]: a primeira em que a legenda cabe nos eixos e nao cobre linhas, marcadores ou
    barras; senao, a que cobre menos"""
    fig = ax.figure
    r = fig.canvas.get_renderer()
    pts, caixas = _obstaculos(ax, r)
    folga = r.points_to_pixels(1.5)
    bb_ax = ax.get_window_extent(r)
    melhor = None
    for op in opcoes:
        loc, extra = (op, {}) if isinstance(op, str) else op
        leg = ax.legend(handles, rotulos, loc=loc, **{**kw, **extra})
        bb = leg.get_window_extent(r)
        fora = bb.x0 < bb_ax.x0 - 1 or bb.x1 > bb_ax.x1 + 1 or bb.y0 < bb_ax.y0 - 1 or bb.y1 > bb_ax.y1 + 1
        n = _conta_colisoes(bb, pts, caixas, folga) + (10 ** 9 if fora else 0)  # fora dos eixos: so em ultimo caso
        if n == 0:
            return leg
        if melhor is None or n < melhor[0]:
            melhor = (n, loc, extra)
        leg.remove()
    return ax.legend(handles, rotulos, loc=melhor[1], **{**kw, **melhor[2]})


def _rotulo_linha(ax, y, texto, preferencias, linha):
    """rotulo de uma linha tracejada horizontal (y em dados), junto da linha (3 pt acima ou abaixo), a esquerda ou a
    direita, na primeira posicao que nao cobre dados nem a legenda e fica dentro dos eixos; None se nao ha posicao
    livre (dados sobre a linha ou linha rente ao eixo)"""
    fig = ax.figure
    r = fig.canvas.get_renderer()
    pts, caixas = _obstaculos(ax, r, ignorar=(linha,))
    if ax.get_legend() is not None:
        caixas = caixas + [ax.get_legend().get_window_extent(r)]
    bb_ax = ax.get_window_extent(r)
    folga = r.points_to_pixels(1.0)
    tr = blended_transform_factory(ax.transAxes, ax.transData)
    # posicoes: perto da ponta do eixo (como no Ceron) e, se os dados estao la, deslocadas em passos de 2% da largura
    xs = {"esq": np.arange(0.04, 0.42, 0.02), "dir": np.arange(0.96, 0.58, -0.02)}

    def poe(lado_x, x, lado_y, dy):
        return ax.annotate(texto, xy=(x, y), xycoords=tr, xytext=(0, dy if lado_y == "acima" else -dy),
                           textcoords="offset points", ha="left" if lado_x == "esq" else "right",
                           va="bottom" if lado_y == "acima" else "top", fontsize=TAM_ROTULO_LINHA)

    for lado_x, lado_y in preferencias:
        for x in xs[lado_x]:
            t = poe(lado_x, x, lado_y, 3.0)
            bb = t.get_window_extent(r)
            fora = bb.x0 < bb_ax.x0 or bb.x1 > bb_ax.x1 or bb.y0 < bb_ax.y0 or bb.y1 > bb_ax.y1
            if not fora and _conta_colisoes(bb, pts, caixas, folga) == 0:
                return t
            t.remove()
    return None


def _verifica_textos(fig, nome):
    """avisos: textos fora da figura ou sobrepostos (marcacoes, rotulos, legendas, rotulos das linhas e letras)"""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    larg_px, alt_px = fig.canvas.get_width_height()
    caixas = []
    for k, ax in enumerate(fig.axes):
        for eixo in (ax.xaxis, ax.yaxis):
            v0, v1 = sorted(eixo.get_view_interval())
            tol = 1e-9 * max(1.0, abs(v1 - v0))
            for loc, t in zip(eixo.get_majorticklocs(), eixo.get_ticklabels()):
                if v0 - tol <= loc <= v1 + tol and t.get_visible() and t.get_text():
                    caixas.append((f"eixos{k} marcação {t.get_text()}", t.get_window_extent(r)))
            if eixo.label.get_text():
                caixas.append((f"eixos{k} rótulo", eixo.label.get_window_extent(r)))
        if ax.get_legend() is not None:
            caixas.append((f"eixos{k} legenda", ax.get_legend().get_window_extent(r)))
        for t in ax.texts:
            caixas.append((f"eixos{k} texto {t.get_text()}", t.get_window_extent(r)))
    for t in fig.texts:
        caixas.append((f"letra {t.get_text()}", t.get_window_extent(r)))
    avisos = []
    for nome_c, bb in caixas:
        if bb.x0 < -0.5 or bb.y0 < -0.5 or bb.x1 > larg_px + 0.5 or bb.y1 > alt_px + 0.5:
            avisos.append(f"{nome}: {nome_c} sai da figura")
    for i in range(len(caixas)):
        for j in range(i + 1, len(caixas)):
            a, b = caixas[i][1], caixas[j][1]
            if a.x1 > b.x0 + 0.5 and b.x1 > a.x0 + 0.5 and a.y1 > b.y0 + 0.5 and b.y1 > a.y0 + 0.5:
                avisos.append(f"{nome}: {caixas[i][0]} sobre {caixas[j][0]}")
    return avisos


# ---------------------------------------------------------------------------------------------------------------------
# dados de um solo


def _valor(S, fam, C):
    return float(C[S["campo_familia"][fam]])


def _familia(S, fam):
    """[(valor, caso, id)] da familia em ordem crescente do valor"""
    lst = [(_valor(S, fam, S["casos"][i]), S["casos"][i], i) for i in S["familias"][fam]]
    return sorted(lst, key=lambda t: t[0])


def _id_ref(S, fam):
    return S["ref_chuva"] if fam == "chuva" else S["ref"]


def _novos_valores(est):
    """valores de cada familia que o Ceron nao tem (uniao dos solos, em ordem crescente): recebem os proximos
    estilos livres, os mesmos nos dois solos"""
    novos = {}
    for fam in ("covk", "covc", "covphi", "s", "alfa", "chuva"):
        vals = set()
        for S in est["solos"].values():
            if fam in S.get("familias", {}):
                vals.update(_valor(S, fam, S["casos"][i]) for i in S["familias"][fam])
        ceron = EC.VALORES_CERON.get(fam, [])
        novos[fam] = sorted(v for v in vals if not any(abs(v - c) < 1e-9 for c in ceron))
    return novos


def _estilo(fam, valor, novos):
    """(cor, tracejado) de um valor: o que o Ceron usa para o mesmo valor; um valor novo recebe o proximo livre"""
    ceron = EC.VALORES_CERON.get(fam, [])
    for i, c in enumerate(ceron):
        if abs(valor - c) < 1e-9:
            return EC.ESTILOS[i]
    i = len(ceron) + min(range(len(novos[fam])), key=lambda k: abs(novos[fam][k] - valor))
    return EC.ESTILOS[i % len(EC.ESTILOS)]


def _rotulo_valor(fam, v):
    if fam in ("covk", "covc", "covphi"):
        return f"{_pct(v)}%"
    return _virg(v)


def _txt_valor(fam, v):
    """'CoV(c) = 0,7', 's = 400', 'I = 9 mm/h' (texto simples das notas)"""
    return f"{NOME_FAMILIA[fam]} = {_virg(v)}{UNID_FAMILIA.get(fam, '')}"


def _txt_N(fam, serie, id_ref):
    """'N = 300 realizações por caso', 'N = 1000 no caso de referência e 300 nos demais' ou a lista por valor"""
    ns = [int(C["N"]) for v, C, i in serie]
    if len(set(ns)) == 1:
        return f"N = {ns[0]} realizações por caso"
    outros = {int(C["N"]) for v, C, i in serie if i != id_ref}
    nref = [int(C["N"]) for v, C, i in serie if i == id_ref]
    if len(outros) == 1 and nref:
        return f"N = {nref[0]} realizações no caso de referência e {outros.pop()} nos demais"
    return "N = " + _lista([f"{int(C['N'])} ({_txt_valor(fam, v)})" for v, C, i in serie])


def _txt_fs(casos):
    return "FS: " + "; ".join(sorted({C["fs_def_texto"] for C in casos}))


def _notas_negativos(fam, serie, painel, efeito):
    """realizacoes com FS = -1 (sem equilibrio nem com c'/0,6): contadas em n_f, fora do ajuste"""
    out = []
    for v, C, i in serie:
        neg = int(C.get("n_fs_negativo") or 0)
        if neg:
            q = "1 realização" if neg == 1 else f"{neg} realizações"
            out.append(f"({painel}) {_txt_valor(fam, v)}: {q} com FS = -1 (sem equilíbrio), {efeito}.")
    return out


# ---------------------------------------------------------------------------------------------------------------------
# paineis


def _handle_ref():
    return Line2D([], [], linestyle="none", marker=MARCA_REF, markersize=MS_REF, markeredgewidth=MEW_REF, color="k")


def _painel_pdf(ax, chave, painel, serie, fam, novos, eixos, titulo=None):
    """PDFs log-normais ajustadas: serie = [(valor, caso)]"""
    e = EC.FIG[chave]
    sup = [_suporte_pdf(C["lognormal"]["mu_ln"], C["lognormal"]["sigma_ln"]) for v, C in serie]
    xa, xb, _ = eixos.aplica(ax, painel, "x", e["x"], (min(s[1] for s in sup), max(s[2] for s in sup)))
    eixos.aplica(ax, painel, "y", e["y"], (0.0, max(s[0] for s in sup)))
    x = np.linspace(xa, xb, 6001)
    handles, rotulos = [], []
    for v, C in serie:
        cor, ls = _estilo(fam, v, novos)
        ln = C["lognormal"]
        h, = ax.plot(x, _pdf_lognormal(x, ln["mu_ln"], ln["sigma_ln"]), color=cor, linestyle=ls, linewidth=LW_PDF,
                     zorder=3)
        handles.append(h)
        rotulos.append(_rotulo_valor(fam, v))
    ax.set_xlabel(e["xrot"])
    _rotulo_y(ax, e["yrot"])
    _legenda(ax, handles, rotulos, ["upper right", "upper left", "center right", "center left", "lower right",
                                    ("upper right", {"ncol": 2}), ("upper left", {"ncol": 2})],
             title=titulo if titulo is not None else e["legenda"])


def _painel_pontos(ax, chave, painel, pontos, eixos, xlim=None, xmargem=None, tracejada=None, rot_x=None,
                   pref_rotulo=None):
    """pontos = [(x, y, e_referencia)] ligados por segmentos de reta em ordem crescente de x; tracejada = (y, texto)
    de uma linha horizontal (caso fora do eixo x, como no Ceron)"""
    e = EC.FIG[chave]
    lim_x = xlim if xlim is not None else e["x"]
    xs = np.array([p[0] for p in pontos], dtype=float)
    ys = np.array([p[1] for p in pontos], dtype=float)
    ordem = np.argsort(xs, kind="stable")
    xs, ys = xs[ordem], ys[ordem]
    refs = np.array([pontos[i][2] for i in ordem], dtype=bool)
    ydados = list(ys) + ([tracejada[0]] if tracejada is not None else [])
    # margem do eixo x alem da marcacao extrema so do lado em que ha ponto nela (o Ceron nao corta os marcadores)
    a, b, _, _, _ = _amplia(lim_x, (float(xs.min()), float(xs.max())))
    tol = 1e-9 * max(1.0, abs(b - a))
    if xmargem is None:
        xmargem = (MARGEM_X if np.any(np.abs(xs - a) < tol) else 0.0,
                   MARGEM_X if np.any(np.abs(xs - b) < tol) else 0.0)
    eixos.aplica(ax, painel, "x", lim_x, (float(xs.min()), float(xs.max())), margem=xmargem)
    eixos.aplica(ax, painel, "y", e["y"], (float(min(ydados)), float(max(ydados))))
    ax.plot(xs, ys, color="k", linewidth=LW_LINHA, zorder=3)
    ax.plot(xs[~refs], ys[~refs], linestyle="none", marker="o", markersize=MS_PONTO, color="k", zorder=4,
            clip_on=False)
    ax.plot(xs[refs], ys[refs], linestyle="none", marker=MARCA_REF, markersize=MS_REF, markeredgewidth=MEW_REF,
            color="k", zorder=5, clip_on=False)
    linha = None
    if tracejada is not None:
        linha = ax.axhline(tracejada[0], color="k", linestyle="--", linewidth=LW_LINHA, zorder=2)
    ax.set_xlabel(rot_x if rot_x is not None else e["xrot"])
    _rotulo_y(ax, e["yrot"])
    locais = ["upper right", "upper left", "lower right", "lower left", "center right", "center left", "upper center",
              "lower center"]
    _legenda(ax, [_handle_ref()], [REF_ROTULO], locais)
    if tracejada is not None and _rotulo_linha(ax, tracejada[0], tracejada[1], pref_rotulo, linha) is None:
        # sem lugar livre junto da linha (dados sobre ela ou linha rente ao eixo): o rotulo vai para a legenda
        h = Line2D([], [], color="k", linestyle="--", linewidth=LW_LINHA)
        _legenda(ax, [_handle_ref(), h], [REF_ROTULO, tracejada[1]], locais)


# ---------------------------------------------------------------------------------------------------------------------
# figuras


def _salva(fig, saida_fig, nome, feitas, avisos):
    avisos += _verifica_textos(fig, nome)
    fig.savefig(os.path.join(saida_fig, nome + ".pdf"))
    fig.savefig(os.path.join(saida_fig, nome + ".png"), dpi=150)
    plt.close(fig)
    feitas.append(os.path.join(saida_fig, nome + ".pdf"))


def _fig25(solo, S, saida_fig, reg, feitas, avisos):
    C = S["casos"][S["ref"]]
    eixos = Eixos(reg, f"{solo}_fig25")
    fig, ah, alt = _figura(1, ASPECTO_25)
    xc = _x_colunas()
    ax_a = _eixos(fig, alt, xc[0], 0, ah)
    ax_b = _eixos(fig, alt, xc[1], 0, ah)
    fs, pos = _fs_validos(C)
    N = len(fs)
    # (a) histograma em densidade, classes de 0,025 com bordas em k x 0,025 - 0,003125 (4 valores da grade por classe)
    h, d0 = 0.025, 0.003125
    k0 = int(math.floor((pos.min() + d0) / h + 1e-9))
    k1 = int(math.floor((pos.max() + d0) / h + 1e-9))
    bordas = np.arange(k0, k1 + 2) * h - d0
    cont, _ = np.histogram(pos, bins=bordas)
    dens = cont / (len(pos) * h)
    ln = C["lognormal"]
    pico, xlo, xhi = _suporte_pdf(ln["mu_ln"], ln["sigma_ln"])
    det = C.get("det")
    nz = np.nonzero(cont)[0]
    xs_d = [bordas[nz[0]], bordas[nz[-1] + 1], xlo, xhi] + ([float(det["FS"])] if det else [])
    e = EC.FIG["25a"]
    xa, xb, _ = eixos.aplica(ax_a, "a", "x", e["x"], (min(xs_d), max(xs_d)))
    eixos.aplica(ax_a, "a", "y", e["y"], (0.0, max(float(dens.max()), pico)))
    ax_a.bar(bordas[:-1], dens, width=h, align="edge", color=CINZA_HIST, linewidth=0, zorder=1)
    x = np.linspace(xa, xb, 6001)
    h_ln, = ax_a.plot(x, _pdf_lognormal(x, ln["mu_ln"], ln["sigma_ln"]), color="k", linewidth=LW_25, zorder=3)
    handles = [Patch(facecolor=CINZA_HIST, edgecolor="none"), h_ln]
    rotulos = ["Histograma", "Log-normal"]
    if det:
        h_det = ax_a.axvline(float(det["FS"]), color="k", linestyle="--", linewidth=LW_25, zorder=4)
        handles.append(h_det)
        rotulos.append("Solução\ndeterminística")
    ax_a.set_xlabel(e["xrot"])
    _rotulo_y(ax_a, e["yrot"])
    estreita = {"handlelength": 1.6, "borderpad": 0.3, "handletextpad": 0.6}
    _legenda(ax_a, handles, rotulos, ["upper right", ("upper right", estreita),
                                      ("upper right", {**estreita, "fontsize": 8.5}), "upper left",
                                      ("upper left", estreita), "center right", "center left"])
    # (b) convergencia de P_f: fracao de FS < 1 entre as n primeiras realizacoes (ordem crescente do numero)
    n = np.arange(1, N + 1)
    pf = np.cumsum(fs < 1.0) / n
    e = EC.FIG["25b"]
    xa, xb, _ = eixos.aplica(ax_b, "b", "x", e["x"], (1.0, float(N)), log=True)
    eixos.aplica(ax_b, "b", "y", e["y"], (float(pf.min()), float(pf.max())))
    ax_b.plot(n, pf, color="k", linewidth=LW_LINHA, zorder=3)
    ax_b.set_xlabel(e["xrot"])
    _rotulo_y(ax_b, e["yrot"])
    yb0, yb1 = ax_b.get_ylim()
    pf_fim = float(pf[-1])
    # valor final com seta cinza ate o fim da curva (no Ceron, acima e a esquerda do fim, em N = 10^4; aqui, com
    # N < 10^4, a direita do fim da curva, onde nao ha curva)
    lx, lxb = math.log10(N), math.log10(xb)
    xt = 10 ** min((lx + lxb) / 2.0, lxb - 0.35) if lxb - lx > 0.7 else 10 ** (lx - 1.3)
    dy = 0.25 * (yb1 - yb0)
    yt = pf_fim + dy if pf_fim + dy <= yb1 - 0.12 * (yb1 - yb0) else pf_fim - dy
    ax_b.annotate(f"{_pct(pf_fim, 2)}%", xy=(N, pf_fim), xytext=(xt, yt), ha="center", va="center", fontsize=9.5,
                  arrowprops=dict(arrowstyle="->", color=CINZA_SETA, linewidth=1.1, shrinkA=4, shrinkB=1))
    _letras(fig, [(ax_a, "a"), (ax_b, "b")])
    _salva(fig, saida_fig, f"{solo}_fig25", feitas, avisos)
    # notas
    reais = C["amostras"]["real"]
    notas = [f"Caso de referência do solo {solo}: {C['fs_def_texto']} (fator de segurança do SRM); N = {N} "
             f"realizações (números {min(reais)} a {max(reais)})."]
    txt = ("(a) Histograma em densidade, com classes de largura 0,025 e bordas em múltiplos de 0,025 menos 0,003125 "
           "(4 valores da grade de 0,00625 do FS por classe)")
    fora = int(C.get("n_fs_fora_grade") or 0)
    if fora == 1:
        txt += f"; 1 dos {N} valores de FS não é múltiplo de 0,00625"
    elif fora > 1:
        txt += f"; {fora} dos {N} valores de FS não são múltiplos de 0,00625"
    neg = len(fs) - len(pos)
    if neg:
        txt += ("; 1 realização" if neg == 1 else f"; {neg} realizações") + " com FS = -1 fora do histograma e do ajuste"
    txt += (f". Curva: distribuição log-normal ajustada por máxima verossimilhança (média {_virg(ln['media'], 3)}, "
            f"mediana {_virg(ln['mediana'], 3)}, CoV {_virg(ln['cov'], 4)}).")
    if det:
        txt += (f" Linha tracejada: FS da solução determinística (valores médios das propriedades, sem campos "
                f"aleatórios), {_virg(float(det['FS']))}.")
    notas.append(txt)
    nf = int(C["n_f"])
    igual1 = int(C.get("n_fs_igual_1") or 0)
    txt = (f"(b) Pf após n simulações: número de realizações com FS < 1 entre as n primeiras (em ordem crescente do "
           f"número da realização) dividido por n; valor final {_virg(pf_fim, 4)} ({nf} de {N})")
    if igual1:
        txt += f"; {igual1} realizações com FS igual a 1 não entram na contagem"
    notas.append(txt + ".")
    if det:
        vers = [v.split("_")[0] for v in sorted(C.get("versoes", {}))]
        vd = str(det.get("versao", "")).split("_")[0]
        if vd and vers:
            notas.append(f"Versão do programa: realizações {_lista(vers)}; rodada determinística {vd}.")
    return notas + eixos.notas


def _fs_validos(C):
    fs = np.asarray(C["amostras"]["FS"], dtype=float)
    return fs, fs[fs > 0]


def _fig26(solo, S, saida_fig, reg, feitas, avisos):
    C = S["casos"][S["ref"]]
    sp = C.get("spearman") or {}
    eixos = Eixos(reg, f"{solo}_fig26")
    fig, ah, alt = _figura(1, ASPECTO_26, base=BASE_SEM_LETRA, larg_eixos=AW_26)
    ax = _eixos(fig, alt, _x_centro(AW_26), 0, ah, base=BASE_SEM_LETRA, larg_eixos=AW_26)
    e = EC.FIG["26"]
    vals = sp.get("segunda_ordem") or {}
    ordem = [("c", 0.0), ("phi", 1.0), ("k", 2.0)]  # de baixo para cima: c, phi, k_v (como no Ceron)
    presentes = [vals[k] for k, _ in ordem if vals.get(k) is not None]
    eixos.aplica(ax, "", "x", e["x"], (min(presentes + [0.0]), max(presentes + [0.0])))
    for k, y in ordem:
        if vals.get(k) is not None:
            ax.barh(y, vals[k], height=0.5, color=EC.CORES_SPEARMAN[k], linewidth=0, zorder=2)
    ax.axvline(0.0, color="k", linewidth=ESTILO["axes.linewidth"], zorder=3)
    ax.set_ylim(-0.375, 2.375)
    ax.set_yticks([])
    ax.set_xlabel(e["xrot"], fontsize=12.0)
    _rotulo_y(ax, e["yrot"], fontsize=11.0)
    rot = {"k": r"$k_v$", "phi": r"$\varphi$", "c": r"$c$"}
    hs = [Patch(facecolor=EC.CORES_SPEARMAN[k], edgecolor="none") for k in ("k", "phi", "c")]
    _legenda(ax, hs, [rot[k] for k in ("k", "phi", "c")],
             ["upper left", "lower left", "upper right", "lower right", "center left", "center right"],
             fontsize=11.0)
    _salva(fig, saida_fig, f"{solo}_fig26", feitas, avisos)
    notas = []
    if sp:
        partes = [f"{nome}: {_virg(vals[k], 3)}" for k, nome in (("c", "c"), ("phi", "phi"), ("k", "kv"))
                  if vals.get(k) is not None]
        notas.append("Coeficiente de correlação de postos de Spearman parcial de segunda ordem (Eqs. 70 e 71 do "
                     "Ceron) entre o FS e as variáveis equivalentes de cada realização, no instante do FS: c = coesão "
                     "média ao longo da banda de ruptura, phi = ângulo de atrito médio ao longo da banda, kv = "
                     "permeabilidade média da massa deslizante. "
                     f"Caso de referência do solo {solo}, {sp['n']} realizações. Valores: " + "; ".join(partes) + ".")
        falta = [nome for k, nome in (("c", "c"), ("phi", "phi"), ("k", "kv")) if vals.get(k) is None]
        if falta:
            notas.append("Sem barra (variável constante em todas as realizações): " + _lista(falta) + ".")
    return notas + eixos.notas


def _fig27_30(solo, S, novos, saida_fig, legs, reg, feitas, avisos):
    fams = [("covk", "a"), ("covc", "b"), ("covphi", "c")]
    series = {fam: _familia(S, fam) for fam, _ in fams}
    ref = S["ref"]
    todos = [(v, C, i) for fam, _ in fams for v, C, i in series[fam]]
    casos = [C for v, C, i in todos]
    xc = _x_colunas()
    pos = {"a": (xc[0], 0), "b": (xc[1], 0), "c": (_x_centro(), 1)}
    ns = {int(C["N"]) for C in casos}
    nref = int(S["casos"][ref]["N"])
    outros = sorted({int(C["N"]) for v, C, i in todos if i != ref})
    txt_n = (f"N = {nref} realizações por caso" if len(ns) == 1 else
             f"N = {nref} realizações no caso de referência e " + _lista([str(n) for n in outros]) + " nos demais")
    # Fig. 27: PDFs
    eixos = Eixos(reg, f"{solo}_fig27")
    fig, ah, alt = _figura(2, ASPECTO)
    paineis = []
    for fam, p in fams:
        ax = _eixos(fig, alt, pos[p][0], pos[p][1], ah)
        _painel_pdf(ax, f"27{p}", p, [(v, C) for v, C, _ in series[fam]], fam, novos, eixos)
        paineis.append((ax, p))
    _letras(fig, paineis)
    _salva(fig, saida_fig, f"{solo}_fig27", feitas, avisos)
    notas = [f"Curvas: distribuições log-normais ajustadas ao FS de cada caso por máxima verossimilhança. {txt_n}; "
             f"{_txt_fs(casos)}."]
    for fam, p in fams:
        notas += _notas_negativos(fam, series[fam], p, "fora do ajuste")
    legs[f"{solo}_fig27"] = notas + eixos.notas
    # Figs. 28 a 30: P_f, mu e CoV contra o CoV
    for num, chave_y in ((28, "Pf"), (29, "media"), (30, "cov")):
        eixos = Eixos(reg, f"{solo}_fig{num}")
        fig, ah, alt = _figura(2, ASPECTO)
        paineis, tracos = [], []
        for fam, p in fams:
            ax = _eixos(fig, alt, pos[p][0], pos[p][1], ah)
            e = EC.FIG[f"{num}{p}"]
            pontos, tracejada = [], None
            for v, C, i in series[fam]:
                y = _y_estat(C, chave_y)
                if v < e["x"][0] - 1e-12:  # CoV(k_v) = 0: fora do eixo x (de 0,6 a 1,0), linha tracejada
                    tracejada = (y, r"$CoV(k_v) = 0$")
                    tracos.append(f"({p}) Linha tracejada: {_txt_valor(fam, v)} ({NOME_Y[chave_y]} = "
                                  f"{_y_txt(chave_y, y)}), fora da faixa do eixo x, como no Ceron.")
                    continue
                pontos.append((v, y, i == ref))
            _painel_pontos(ax, f"{num}{p}", p, pontos, eixos, tracejada=tracejada,
                           pref_rotulo=[("esq", "abaixo"), ("esq", "acima"), ("dir", "abaixo"), ("dir", "acima")])
            paineis.append((ax, p))
        _letras(fig, paineis)
        _salva(fig, saida_fig, f"{solo}_fig{num}", feitas, avisos)
        if chave_y == "Pf":
            txt = "Pf = nf/N, com nf o número de realizações com FS < 1"
        elif chave_y == "media":
            txt = "Média do FS: a da distribuição log-normal ajustada por máxima verossimilhança"
        else:
            txt = "CoV do FS: o da distribuição log-normal ajustada por máxima verossimilhança"
        notas = [f"{txt}; {txt_n}; {_txt_fs(casos)}."] + tracos
        if chave_y == "Pf":
            zeros = [_txt_valor(fam, v) for fam, _ in fams for v, C, i in series[fam] if int(C["n_f"]) == 0]
            if zeros:
                notas.append("Pf = 0 (nenhuma realização com FS < 1) em " + _lista(zeros) + ".")
        for fam, p in fams:
            notas += _notas_negativos(fam, series[fam], p, "contada em nf" if chave_y == "Pf" else "fora do ajuste")
        legs[f"{solo}_fig{num}"] = notas + eixos.notas


def _y_estat(C, chave):
    if chave == "Pf":
        return float(C["Pf"])
    return float(C["lognormal"][chave])


def _fig31_33(num, fam, solo, S, novos, saida_fig, legs, reg, feitas, avisos):
    serie = _familia(S, fam)
    ref = _id_ref(S, fam)
    eixos = Eixos(reg, f"{solo}_fig{num}")
    fig, ah, alt = _figura(2, ASPECTO)
    xc = _x_colunas()
    axs = {"a": _eixos(fig, alt, xc[0], 0, ah), "b": _eixos(fig, alt, xc[1], 0, ah),
           "c": _eixos(fig, alt, xc[0], 1, ah), "d": _eixos(fig, alt, xc[1], 1, ah)}
    titulo = ROT_CHUVA if fam == "chuva" else None
    _painel_pdf(axs["a"], f"{num}a", "a", [(v, C) for v, C, _ in serie], fam, novos, eixos, titulo=titulo)
    tracos = []
    for p, chave_y in (("b", "Pf"), ("c", "cov"), ("d", "media")):
        e = EC.FIG[f"{num}{p}"]
        pontos, tracejada = [], None
        for v, C, i in serie:
            y = _y_estat(C, chave_y)
            if fam == "s" and v > e["x"][1] + 1e-12:  # s = 400: fora do eixo x (de 0 a 20), linha tracejada
                tracejada = (y, rf"$s = {_virg(v)}$")
                tracos.append((p, v, chave_y, y))
                continue
            pontos.append((v, y, i == ref))
        if fam == "s":
            pref = [("esq", "acima"), ("dir", "abaixo"), ("esq", "abaixo"), ("dir", "acima")]
            _painel_pontos(axs[p], f"{num}{p}", p, pontos, eixos, xmargem=(0.0, MARGEM_X_31), tracejada=tracejada,
                           pref_rotulo=pref)
        elif fam == "chuva":
            _painel_pontos(axs[p], f"{num}{p}", p, pontos, eixos, xlim=LIM_CHUVA, rot_x=ROT_CHUVA)
        else:
            _painel_pontos(axs[p], f"{num}{p}", p, pontos, eixos)
    _letras(fig, [(axs[p], p) for p in "abcd"])
    _salva(fig, saida_fig, f"{solo}_fig{num}", feitas, avisos)
    casos = [C for v, C, i in serie]
    vals = [v for v, C, i in serie]
    base = (f"(a) PDFs log-normais ajustadas por máxima verossimilhança; (b) Pf = nf/N, com nf o número de "
            f"realizações com FS < 1; (c) e (d) CoV e média da log-normal ajustada. {_txt_N(fam, serie, ref)}; "
            f"{_txt_fs(casos)}.")
    if fam == "s":
        notas = [f"Distâncias de autocorrelação (Lx, Ly) = s (20 m, 2 m), com s = {_lista_num(vals)} (no Ceron: 1; "
                 f"1,5; 2; 5; 10; 20 e 400). " + base]
        if tracos:
            dentro = [v for v in vals if v <= EC.FIG[f"{num}b"]["x"][1]]
            partes = [f"{NOME_Y[ch]} = {_y_txt(ch, y)} em ({p})" for p, v, ch, y in tracos]
            v400 = tracos[0][1]
            notas.append(f"(b) a (d) Linha tracejada: s = {_virg(v400)}, fora da faixa do eixo x, como no Ceron; "
                         + _lista(partes) + f". Os pontos s = {_lista_num(dentro)} estão ligados por segmentos de "
                         "reta.")
    elif fam == "alfa":
        notas = [f"Anisotropia de permeabilidade alfa = kh/kv = {_lista_num(vals)} (no Ceron: 1 a 5). " + base]
    else:
        notas = [f"Intensidade da chuva I = {_lista_num(vals, ' mm/h')} (chuva de 24 h), no lugar da razão de "
                 f"rebaixamento hw/H do Ceron; eixo x de (b) a (d): I (mm/h), sem correspondente no Ceron (no Ceron, "
                 f"hw/H de 0,5 a 1,0). " + base]
        if not S.get("mesmos_dados_ref_chuva", False):
            notas.append(_nota_ref_chuva(solo, S["casos"][ref], S["casos"][S["ref"]]))
        else:
            notas.append(f"Solo {solo}: o caso I = {_virg(float(S['casos'][ref]['I_mm_h']))} mm/h (estrela) é o caso "
                         f"de referência das Figs. 25 a 32 (os mesmos dados, {int(S['casos'][ref]['N'])} "
                         f"realizações).")
        zeros = [_txt_valor(fam, v) for v, C, i in serie if int(C["n_f"]) == 0]
        if zeros:
            notas.append("(b) Pf = 0 (nenhuma realização com FS < 1) em " + _lista(zeros) + ".")
    for p, efeito in (("b", "contada em nf"), ("c", "fora do ajuste"), ("d", "fora do ajuste")):
        notas += _notas_negativos(fam, serie, p, efeito)
    legs[f"{solo}_fig{num}"] = notas + eixos.notas


def _nota_ref_chuva(solo, Cc, C0):
    """o caso I = 36 mm/h da Fig. 33 e o caso de referencia das Figs. 25 a 32, quando nao sao os mesmos dados"""
    txt = (f"Solo {solo}: o caso I = {_virg(float(Cc['I_mm_h']))} mm/h (estrela) não é o mesmo conjunto de dados do "
           f"caso de referência das Figs. 25 a 32: aqui, {int(Cc['N'])} realizações e {Cc['fs_def_texto']}; lá, "
           f"{int(C0['N'])} realizações e {C0['fs_def_texto']}.")
    da = dict(zip(Cc["amostras"]["real"], Cc["amostras"]["FS"]))
    pares = [(da[r], f) for r, f in zip(C0["amostras"]["real"], C0["amostras"]["FS"]) if r in da]
    if pares:
        comuns = [r for r in C0["amostras"]["real"] if r in da]
        iguais = sum(1 for a, b in pares if abs(a - b) < 1e-9)
        txt += (f" Nas {len(pares)} realizações comuns aos dois casos (números {min(comuns)} a {max(comuns)}), o FS "
                f"daqui é igual ao de lá em {iguais}.")
    cont = Cc.get("t_min_h_contagem")
    if cont:
        partes = [f"{_virg(float(t))} h em {int(n)}" for t, n in cont]
        txt += " Instante do FS mínimo, aqui: " + _lista(partes) + " realizações."
    return txt


# ---------------------------------------------------------------------------------------------------------------------


def gera(est, saida):
    """grava em <saida>/figuras as Figs. 25 a 33 de cada solo (PDF e PNG) e legendas_estat.json; devolve
    {"figuras": [...], "ampliacoes": [...], "avisos": [...]}"""
    if isinstance(est, (str, os.PathLike)):
        with open(est, encoding="utf-8") as f:
            est = json.load(f)
    saida_fig = os.path.join(saida, "figuras")
    os.makedirs(saida_fig, exist_ok=True)
    novos = _novos_valores(est)
    legs, reg, feitas, avisos = {}, [], [], []
    with plt.rc_context(ESTILO):
        for solo in ("f50", "f100"):
            S = est["solos"].get(solo)
            if S is None:
                continue
            fams = S.get("familias", {})
            if S.get("ref") in S.get("casos", {}):
                legs[f"{solo}_fig25"] = _fig25(solo, S, saida_fig, reg, feitas, avisos)
                legs[f"{solo}_fig26"] = _fig26(solo, S, saida_fig, reg, feitas, avisos)
            else:
                avisos.append(f"{solo}: sem caso de referência; Figs. 25 e 26 não feitas")
            if all(fams.get(f) for f in ("covk", "covc", "covphi")):
                _fig27_30(solo, S, novos, saida_fig, legs, reg, feitas, avisos)
            else:
                avisos.append(f"{solo}: falta uma das famílias covk, covc, covphi; Figs. 27 a 30 não feitas")
            for num, fam in ((31, "s"), (32, "alfa"), (33, "chuva")):
                if fams.get(fam):
                    _fig31_33(num, fam, solo, S, novos, saida_fig, legs, reg, feitas, avisos)
                else:
                    avisos.append(f"{solo}: sem a família {fam}; Fig. {num} não feita")
    legs = {k: legs[k] for k in sorted(legs, key=lambda k: (k.split("_")[0] != "f50", k))}
    with open(os.path.join(saida_fig, "legendas_estat.json"), "w", encoding="utf-8") as f:
        json.dump(legs, f, ensure_ascii=False, indent=1)
    for a in avisos:
        print("AVISO:", a, file=sys.stderr)
    amp = [{"figura": f, "painel": p, "eixo": e, "ceron": c, "novo": n, "passo_corrigido": pc}
           for f, p, e, c, n, pc in reg]
    return {"figuras": feitas, "ampliacoes": amp, "avisos": avisos}


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(__doc__)
        sys.exit(1)
    saida = sys.argv[1]
    r = gera(os.path.join(saida, "estatisticas.json"), saida)
    print("\n".join(r["figuras"]))
    for a in r["ampliacoes"]:
        print(f"{a['figura']} ({a['painel']}) eixo {a['eixo']}: {a['ceron']} -> {a['novo']}"
              + (" (passo corrigido)" if a["passo_corrigido"] else ""))
