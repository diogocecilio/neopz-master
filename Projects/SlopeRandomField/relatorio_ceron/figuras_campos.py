"""figuras_campos.py: as figuras de campos e de mecanismos de Vargas Ceron et al. (2025), Int J Numer Anal Methods
Geomech 49:2459-2491 (doi 10.1002/nag.3993), para cada solo (f50 e f100), a partir dos arquivos de cada realizacao:

  <solo>_fig12  realizacao tipica (r00000 do caso de referencia): k_v (na unidade do Ceron, k / 9,81), c' e phi'
  <solo>_fig13  poropressao p_w no instante do FS dessa realizacao (o Ceron mostra o excesso de poropressao)
  <solo>_fig14  p_w, c' e phi' com a banda de ruptura critica dessa realizacao e o seu FS
  <solo>_fig16  as bandas de ruptura de todas as realizacoes da pasta de referencia com mecanismo.csv (cinza
                translucido), a da solucao deterministica (vermelho) e a porcentagem dos modos de ruptura (acima, pelo
                e abaixo do pe do talude); Figs. 16 e 20 do Ceron
  <solo>_fig24  geometria do talude (H, 45 graus, chuva) e a banda de ruptura da solucao deterministica com o seu FS

mais legendas_campos.json (notas das legendas: texto simples, sem caracteres especiais do LaTeX e so com caracteres
Latin-1, para entrar na legenda com ou sem escape) e modos_ruptura.json (modos de ruptura do caso de referencia e de
todos os casos).

Uso: gera(est, raiz_dados, saida) grava em <saida>/figuras. est e o dicionario de dados.calcula() (ou o caminho de
estatisticas.json); com est = None, o modulo le o que precisa (t_min_h, FS deterministico, banda_ymin) dos resumo.csv,
com o mesmo mapa de casos do contrato. Linha de comando:
    python3 -I figuras_campos.py <raiz_dados> <saida> [--sem-est]
(usa <saida>/estatisticas.json se existir, a menos de --sem-est).

Desenho dos campos de campos.csv (valores nos pontos de integracao, 3 por triangulo): triangulacao de Delaunay dos
pontos de integracao e dos nos do contorno da malha (o valor de um no do contorno e o do ponto de integracao mais
proximo; sem eles, sobraria uma faixa em branco entre os pontos e o contorno), com mascara nos triangulos cujo centroide
fica fora do contorno. Verificado em cada desenho: todos os pontos entram na triangulacao, os pontos medios das arestas
dos triangulos restantes estao dentro do dominio ou no contorno, e a area dos triangulos restantes e a area da malha
(1750 m2). No interior, a interpolacao e linear entre os pontos de integracao, sem valores inventados. p_w esta nos nos
(pw.csv) e e desenhado sobre a propria malha (triangulos.csv).

Banda de ruptura de uma realizacao num instante: as linhas de mecanismo.csv desse instante com
epsp_raizJ2 >= 0,1 x maximo (BANDA_LIMIAR = 0,1 de main.cpp; mecanismo.csv guarda os pontos com >= 0,01 x maximo, e o
maximo esta entre eles)."""

import csv
import glob
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
import matplotlib.tri as mtri
import numpy as np
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LinearSegmentedColormap, Normalize, to_rgb
from matplotlib.lines import Line2D
from matplotlib.path import Path
from scipy.spatial import cKDTree

import eixos_ceron as EC

# ---------------------------------------------------------------------------------------------------------------------
# constantes do problema (README e main.cpp de Projects/SlopeRandomField)
BANDA_LIMIAR = 0.1  # banda de ruptura: sqrt(J2) da deformacao plastica >= 0,1 x maximo (main.cpp, BANDA_LIMIAR)
TOL_FS = 5e-6  # mecanismo.csv grava o FS com 7 algarismos significativos (1,3453125 -> 1,345312)
X_PE, Y_PE = 40.0, 20.0  # pe do talude (m)
X_CRISTA, Y_CRISTA = 30.0, 30.0  # crista (m)
H_TALUDE, BETA_GRAUS = 10.0, 45.0
RAIO_PE = 1.0  # tolerancia do modo de ruptura: arestas dos triangulos com vertice a ate 1 m do pe
MEDIA_C, MEDIA_PHI = 10.0, 26.0  # medias de c' (kPa) e phi' (graus) dos campos (README: c' = 10 kPa, phi' = 26 graus)
K_S = {"f50": 1.0e-5, "f100": 1.0e-6}  # k_s (m/s) de cada solo, se est nao o trouxer
ALFA_SOBREPOSICAO = 0.02  # opacidade de um ponto de banda de uma realizacao na Fig. 16
DIAM_PONTO_16 = 0.7  # diametro (m) do ponto de banda das realizacoes na Fig. 16
DIAM_DET_16 = 0.35  # diametro (m) do ponto da banda deterministica na Fig. 16
CINZA_BANDAS = "0.30"
VERMELHO = "#d7191c"
LINHA_TERRENO = "#8b2323"  # o vermelho escuro da superficie do terreno nas Figs. 16 e 20 do Ceron
AZUL_CHUVA = "#1f5fbf"
ROTULO_PW = r"Poropressão $p_w$ (kPa)"  # a grandeza aqui e p_w, nao o excesso de poropressao do Ceron
NIVEIS = 49  # niveis do preenchimento dos campos (o mapa de cores e continuo na barra)

ESTILO = {
    "font.family": "serif",
    "font.serif": ["LM Roman 10", "Latin Modern Roman", "DejaVu Serif"],  # a mesma fonte de figuras_estat.py
    "mathtext.fontset": "cm",
    "axes.unicode_minus": True,
    "font.size": 8,
    "axes.linewidth": 0.6,
    "axes.edgecolor": "black",
    "axes.labelcolor": "black",
    "text.color": "black",
    "xtick.color": "black",
    "ytick.color": "black",
    "xtick.major.width": 0.5,
    "xtick.major.size": 2.0,
    "xtick.labelsize": 6.5,
    "legend.fontsize": 8,
    "axes.grid": False,
}

MODOS = ("acima", "pe", "abaixo")
REGRA_MODO = ("modo de ruptura pelo ponto mais baixo da banda de ruptura (y min) em relacao a cota do pe do talude "
              "(y pe = 20 m), com tolerancia delta = mediana do comprimento das arestas (distintas) dos triangulos da "
              "malha com vertice a ate 1 m do pe (40 m, 20 m): acima do pe se y min maior que y pe + delta; abaixo do "
              "pe se y min menor que y pe - delta; pelo pe nos demais casos")
NOME_MODO = {"acima": "acima do pé", "pe": "pelo pé", "abaixo": "abaixo do pé", None: "sem banda"}

# ---------------------------------------------------------------------------------------------------------------------
# mapa de casos para o modo autonomo (est = None): o mesmo do contrato. Ids: <solo>_ref, <solo>_<familia><valor>, ...
SUFIXOS = [("covk0", "covk", 0.0), ("covk0.75", "covk", 0.75), ("covk0.9", "covk", 0.9), ("covk1", "covk", 1.0),
           ("covc0.1", "covc", 0.1), ("covc0.5", "covc", 0.5), ("covc0.7", "covc", 0.7),
           ("covphi0.05", "covphi", 0.05), ("covphi0.15", "covphi", 0.15), ("covphi0.2", "covphi", 0.2),
           ("lx40_ly4", "s", 2.0), ("lx100_ly10", "s", 5.0), ("lx8000_ly800", "s", 400.0),
           ("alfa2", "alfa", 2.0), ("alfa5", "alfa", 5.0), ("alfa10", "alfa", 10.0)]


def _id_sufixo(solo, suf, fam, val):
    return f"{solo}_s{val:g}" if fam == "s" else f"{solo}_{suf}"


def _mapa_padrao():
    """{solo: {"ref": id, "casos": {id: {"pastas": [...], "fs_def": ..., "det": pasta ou None}}}}, com as pastas
    relativas a raiz_dados"""
    m = {}
    f50 = {"f50_ref": {"pastas": ["sensibilidade/f50_I36_i"], "fs_def": "FS_24h",
                       "det": "monte_carlo/f50_I36_i_ftrans0.5_det"}}
    for suf, fam, val in SUFIXOS:
        f50[_id_sufixo("f50", suf, fam, val)] = {"pastas": [f"sensibilidade/f50_I36_i_{suf}"], "fs_def": "FS_24h",
                                                 "det": None}
    for i in ("9", "36", "80"):
        f50[f"f50_I{i}"] = {"pastas": [f"monte_carlo/f50_I{i}_i_ftrans0.5"], "fs_def": "FS_min",
                            "det": f"monte_carlo/f50_I{i}_i_ftrans0.5_det"}
    m["f50"] = {"ref": "f50_ref", "casos": f50}
    f100 = {"f100_ref": {"pastas": ["f100/f100_I36_i", "f100_fase1/f100_I36_i"], "fs_def": "FS_min",
                         "det": "f100/f100_I36_i_det"}}
    for suf, fam, val in SUFIXOS:
        det = f"f100/f100_I36_i_{suf}_det" if fam == "alfa" else None
        f100[_id_sufixo("f100", suf, fam, val)] = {"pastas": [f"f100/f100_I36_i_{suf}"], "fs_def": "FS_min",
                                                   "det": det}
    for i in ("3.6", "9", "80"):
        f100[f"f100_I{i}"] = {"pastas": [f"f100/f100_I{i}_i", f"f100_fase3/f100_I{i}_i"], "fs_def": "FS_min",
                              "det": f"f100/f100_I{i}_i_det"}
    f100["f100_I36"] = dict(f100["f100_ref"])
    m["f100"] = {"ref": "f100_ref", "casos": f100}
    return m


# ---------------------------------------------------------------------------------------------------------------------
# leitura


def _caminho(raiz, p):
    return p if os.path.isabs(p) else os.path.join(raiz, p)


def _concluida(pasta_real):
    r = os.path.join(pasta_real, "resultado.txt")
    return os.path.isfile(r) and os.path.getsize(r) > 0


def _realizacoes(pasta_comb):
    """[(numero, pasta)] das realizacoes concluidas (com resultado.txt nao vazio), em ordem de numero"""
    out = []
    for p in glob.glob(os.path.join(pasta_comb, "r[0-9]*")):
        b = os.path.basename(p)
        if os.path.isdir(p) and b[1:].isdigit() and _concluida(p):
            out.append((int(b[1:]), p))
    return sorted(out)


def _pasta_det(raiz, p):
    """a pasta da realizacao deterministica: p pode ser a da combinacao (..._det) ou a da realizacao (r00000)"""
    p = _caminho(raiz, p)
    if _concluida(p):
        return p
    r = _realizacoes(p)
    if len(r) != 1:
        raise RuntimeError(f"pasta deterministica {p}: esperava 1 realizacao concluida, achei {len(r)}")
    return r[0][1]


def _resumo(pasta_real):
    with open(os.path.join(pasta_real, "resumo.csv"), newline="") as f:
        linhas = list(csv.DictReader(f))
    if len(linhas) != 1:
        raise RuntimeError(f"{pasta_real}/resumo.csv: esperava 1 linha, achei {len(linhas)}")
    return linhas[0]


def _float(s):
    try:
        return float(s)
    except (TypeError, ValueError):
        return float("nan")


def _le_num(arq):
    """CSV numerico com cabecalho -> {coluna: array}"""
    with open(arq) as f:
        cab = f.readline().strip().split(",")
    d = np.loadtxt(arq, delimiter=",", skiprows=1, ndmin=2)
    if d.size == 0:
        d = np.zeros((0, len(cab)))
    return {c: d[:, i] for i, c in enumerate(cab)}


def _instante(fs_def, res):
    """(t_h, FS) que definem o FS da realizacao: FS_min -> (t_min_h, FS_min); FS_<t>h -> (t, FS_<t>h)"""
    if fs_def == "FS_min":
        return _float(res["t_min_h"]), _float(res["FS_min"])
    if fs_def.startswith("FS_") and fs_def.endswith("h"):
        return float(fs_def[3:-1]), _float(res[fs_def])
    raise ValueError(f"fs_def desconhecida: {fs_def}")


def _nas(t, alvo):
    return np.abs(t - alvo) < 1e-6


def _banda(pasta_real, t):
    """banda de ruptura no instante t: ip, x, y, epsp e FS das linhas de mecanismo.csv com epsp >= 0,1 x maximo"""
    m = _le_num(os.path.join(pasta_real, "mecanismo.csv"))
    s = _nas(m["t_h"], t)
    if not s.any():
        raise RuntimeError(f"{pasta_real}/mecanismo.csv sem linhas em t = {t:g} h")
    e = m["epsp_raizJ2"][s]
    b = e >= BANDA_LIMIAR * e.max()
    fs = np.unique(m["FS"][s])
    return {"ip": m["ip"][s][b].astype(int), "x": m["x"][s][b], "y": m["y"][s][b], "epsp": e[b],
            "FS": float(fs[0]) if fs.size == 1 else float("nan")}


class Malha:
    """malha geometrica (nos de pw.csv num instante e triangulos.csv), o seu contorno e a triangulacao dos campos"""

    def __init__(self, pasta_real):
        pw = _le_num(os.path.join(pasta_real, "pw.csv"))
        self.pw = pw
        s = _nas(pw["t_h"], pw["t_h"][0])
        no = pw["no"][s].astype(int)
        n = no.max() + 1
        if no.size != n or np.unique(no).size != n:
            raise RuntimeError(f"{pasta_real}/pw.csv: numeracao de nos nao contigua")
        self.X = np.empty(n)
        self.Y = np.empty(n)
        self.X[no] = pw["x"][s]
        self.Y[no] = pw["y"][s]
        sup = pw["superficie"][s] > 0
        o = np.argsort(pw["x"][s][sup])
        self.superficie = (pw["x"][s][sup][o], pw["y"][s][sup][o])  # a superficie do terreno (chuva), de x = 0 a 70
        tr = _le_num(os.path.join(pasta_real, "triangulos.csv"))
        T = np.column_stack([tr["n0"], tr["n1"], tr["n2"]]).astype(int)
        T = self._anti_horario(self.X, self.Y, T)
        self.T = T
        self.tri = mtri.Triangulation(self.X, self.Y, T)
        self.area = float(np.sum(self._areas(self.X, self.Y, T)))
        # arestas de contorno: as que pertencem a um triangulo so
        a = np.sort(np.vstack([T[:, [0, 1]], T[:, [1, 2]], T[:, [2, 0]]]), axis=1)
        u, c = np.unique(a, axis=0, return_counts=True)
        self.contorno = u[c == 1]
        self.lacos = self._encadeia(self.contorno)
        if len(self.lacos) != 1:
            raise RuntimeError(f"contorno da malha com {len(self.lacos)} lacos (esperava 1)")
        laco = self.lacos[0]
        self.poligono = Path(np.column_stack([np.append(self.X[laco], self.X[laco[0]]),
                                              np.append(self.Y[laco], self.Y[laco[0]])]), closed=True)
        self.trip = None

    @staticmethod
    def _areas(X, Y, T):
        return 0.5 * ((X[T[:, 1]] - X[T[:, 0]]) * (Y[T[:, 2]] - Y[T[:, 0]])
                      - (X[T[:, 2]] - X[T[:, 0]]) * (Y[T[:, 1]] - Y[T[:, 0]]))

    @classmethod
    def _anti_horario(cls, X, Y, T):
        T = T.copy()
        neg = cls._areas(X, Y, T) < 0
        T[neg] = T[neg][:, [0, 2, 1]]
        return T

    @staticmethod
    def _encadeia(arestas):
        """as arestas de contorno em lacos fechados de nos"""
        viz = {}
        for a, b in arestas:
            viz.setdefault(int(a), []).append(int(b))
            viz.setdefault(int(b), []).append(int(a))
        if any(len(v) != 2 for v in viz.values()):
            raise RuntimeError("contorno da malha com no de grau diferente de 2")
        livres = set(viz)
        lacos = []
        while livres:
            ini = min(livres)
            laco = [ini]
            ant, at = None, ini
            while True:
                prox = viz[at][0] if viz[at][0] != ant else viz[at][1]
                if prox == ini:
                    break
                laco.append(prox)
                ant, at = at, prox
            livres -= set(laco)
            lacos.append(np.array(laco))
        return lacos

    def dentro(self, x, y):
        """pontos dentro do dominio ou no contorno (tolerancia de 1e-6 m)"""
        p = np.column_stack([x, y])
        return self.poligono.contains_points(p, radius=1e-6) | self.poligono.contains_points(p, radius=-1e-6)

    def tolerancia_pe(self):
        """mediana do comprimento das arestas (distintas) dos triangulos com vertice a ate RAIO_PE do pe"""
        d = np.hypot(self.X - X_PE, self.Y - Y_PE)
        sel = np.any(d[self.T] <= RAIO_PE, axis=1)
        Ts = self.T[sel]
        a = np.unique(np.sort(np.vstack([Ts[:, [0, 1]], Ts[:, [1, 2]], Ts[:, [2, 0]]]), axis=1), axis=0)
        L = np.hypot(self.X[a[:, 0]] - self.X[a[:, 1]], self.Y[a[:, 0]] - self.Y[a[:, 1]])
        return float(np.median(L)), int(sel.sum()), int(L.size)

    def prepara_campos(self, xip, yip):
        """triangulacao dos campos (ver o cabecalho do modulo), com as verificacoes"""
        if np.any(~self.dentro(xip, yip)):
            raise RuntimeError("pontos de integracao fora do dominio da malha")
        nb = np.unique(self.contorno.ravel())
        _, viz = cKDTree(np.column_stack([xip, yip])).query(np.column_stack([self.X[nb], self.Y[nb]]))
        xs = np.concatenate([xip, self.X[nb]])
        ys = np.concatenate([yip, self.Y[nb]])
        tri = mtri.Triangulation(xs, ys)
        T = tri.triangles
        if np.unique(T).size != xs.size:
            raise RuntimeError("pontos que nao entraram na triangulacao de Delaunay")
        fora = ~self.poligono.contains_points(np.column_stack([xs[T].mean(axis=1), ys[T].mean(axis=1)]))
        tri.set_mask(fora)
        Tv = T[~fora]
        for a, b in ((0, 1), (1, 2), (2, 0)):
            if np.any(~self.dentro(0.5 * (xs[Tv[:, a]] + xs[Tv[:, b]]), 0.5 * (ys[Tv[:, a]] + ys[Tv[:, b]]))):
                raise RuntimeError("triangulo dos campos atravessa o contorno do dominio")
        area = float(np.sum(np.abs(self._areas(xs, ys, Tv))))
        if abs(area - self.area) > 1e-9 * self.area:
            raise RuntimeError(f"triangulacao dos campos com area {area} diferente da malha ({self.area})")
        self.trip = tri
        self._viz = viz
        self.verif_campos = {"pontos": int(xs.size), "triangulos": int(T.shape[0]), "mascarados": int(fora.sum()),
                             "area_m2": area, "area_malha_m2": self.area}

    def valores_campo(self, v_ip):
        """valores nos pontos da triangulacao dos campos: os pontos de integracao e os nos do contorno"""
        return np.concatenate([v_ip, v_ip[self._viz]])


# ---------------------------------------------------------------------------------------------------------------------
# formatacao (virgula decimal)


def _virg(v, casas=None):
    """numero com virgula decimal; casas = None: o minimo de casas (ate 6)"""
    s = f"{v:.6f}".rstrip("0").rstrip(".") if casas is None else f"{v:.{casas}f}"
    if s.startswith("-") and float(s) == 0:
        s = s[1:]
    return s.replace(".", ",")


def _virg_mat(v, casas=None):
    """o mesmo, para mathtext: {,} evita o espaco depois da virgula"""
    return _virg(v, casas).replace(",", "{,}")


def _menos(s):
    return s.replace("-", "\N{MINUS SIGN}")


def _mant_exp(v, alg=3):
    e = int(math.floor(math.log10(abs(v)) + 1e-9))
    m = round(v / 10 ** e, alg - 1)
    if abs(m) >= 10:
        m, e = m / 10, e + 1
    return m, e


def _sci_mat(v, alg=3):
    """1,0194e-6 -> $1{,}02\\times10^{-6}$ (alg algarismos significativos, sem zeros a direita)"""
    if v == 0:
        return "$0$"
    m, e = _mant_exp(v, alg)
    return rf"${_virg_mat(m)}\times10^{{{e}}}$"


def _sci_txt(v, alg=3):
    """1,0194e-6 -> 1,02e-6 (texto simples para as notas)"""
    m, e = _mant_exp(v, alg)
    return f"{_virg(m)}e{e}"


def _fs_txt(v):
    return _virg(v)  # o minimo de casas (ate 6) da o valor exato do FS do SRM


# ---------------------------------------------------------------------------------------------------------------------
# barras de cores do Ceron (eixos_ceron.CAMPOS), ampliadas por eixos_ceron.amplia() quando os dados nao cabem


def _cmap():
    return LinearSegmentedColormap.from_list("ceron", EC.CORES_CAMPOS, N=512)


def _barra(chave, dmin, dmax):
    """(lo, hi, marcas, nota): a faixa do Ceron, ou a ampliada por amplia() para conter [dmin, dmax]"""
    cfg = EC.CAMPOS[chave]
    a, b = cfg["lim"]
    marcas = list(cfg["marcas"])
    passo = marcas[-1] - marcas[-2]
    lo, hi, q, ampliou = EC.amplia((a, b, passo), (dmin, dmax))
    if not ampliou:
        return a, b, marcas, None
    if dmin >= a:  # so o topo passou: o inicio da barra continua o do Ceron
        lo = a
    m = [lo] if abs(lo / q - round(lo / q)) > 1e-9 else []
    k = math.ceil(lo / q - 1e-9)
    while k * q <= hi + 1e-9 * abs(q):
        if not m or k * q - m[-1] > 0.3 * q:
            m.append(k * q)
        k += 1
    fmt = _sci_txt if chave == "k" else _virg
    return lo, hi, m, f"barra de cores ampliada; no Ceron: {fmt(a)} a {fmt(b)}"


def _rotulos_marcas(chave, marcas):
    if chave == "k":
        return [_sci_mat(v) for v in marcas]
    return [_menos(_virg(v)) for v in marcas]


# ---------------------------------------------------------------------------------------------------------------------
# paineis de campo: o dominio (70 m x 30 m) em escala 1:1, sem eixos, com a barra de cores horizontal logo abaixo,
# como nas Figs. 12 a 14 do Ceron. Posicoes em polegadas a partir do canto superior esquerdo da figura.
LARG_DOM, ALT_DOM = 70.0, 30.0
ALT_BARRA, ESP_BARRA = 0.085, 0.025  # altura da barra e espaco entre o campo e a barra (pol)


def _eixos_painel(fig, esq, topo, larg):
    W, H = fig.get_size_inches()
    alt = larg * ALT_DOM / LARG_DOM
    ax = fig.add_axes([esq / W, 1 - (topo + alt) / H, larg / W, alt / H])
    yb = topo + alt + ESP_BARRA
    cax = fig.add_axes([esq / W, 1 - (yb + ALT_BARRA) / H, larg / W, ALT_BARRA / H])
    return ax, cax


def _desenha_campo(fig, ax, cax, tri, z, chave, lo, hi, marcas, rotulo, letra=None):
    """campo preenchido com o mapa de cores do Ceron e a barra de cores horizontal abaixo"""
    cmap = _cmap()
    cs = ax.tricontourf(tri, np.clip(z, lo, hi), levels=np.linspace(lo, hi, NIVEIS), cmap=cmap, vmin=lo, vmax=hi)
    # contorno de cada poligono com a sua propria cor: sem isso, os leitores de PDF mostram frestas brancas entre os
    # niveis (suavizacao das bordas)
    cs.set_edgecolor("face")
    cs.set_linewidth(0.25)
    ax.set_xlim(0, LARG_DOM)
    ax.set_ylim(0, ALT_DOM)
    ax.set_aspect("equal")
    ax.set_axis_off()
    cb = fig.colorbar(ScalarMappable(norm=Normalize(lo, hi), cmap=cmap), cax=cax, orientation="horizontal",
                      ticks=marcas)
    cb.ax.set_xticklabels(_rotulos_marcas(chave, marcas))
    cb.outline.set_linewidth(0.5)
    cb.ax.tick_params(labelsize=6.5, length=2, width=0.5, pad=1.5)
    cb.set_label(rotulo, fontsize=7.5, labelpad=2)
    if letra:
        cax.annotate(letra, xy=(0.5, 0.0), xycoords=cax.xaxis.label, xytext=(0, -6), textcoords="offset points",
                     ha="center", va="top", fontsize=8)


def _texto_estat(ax, txt, tam=6.3):
    """media e CoV dentro do campo, junto da base, como na Fig. 12 do Ceron (fundo claro para ler sobre o azul)"""
    ax.text(35.0, 1.0, txt, ha="center", va="bottom", fontsize=tam, color="black",
            bbox=dict(boxstyle="square,pad=0.12", facecolor="white", edgecolor="none", alpha=0.6))


def _s_marcador(diam_m, larg_pol, faixa_m):
    """area do marcador (pt2) de um ponto com diametro diam_m num eixo de larg_pol polegadas para faixa_m metros"""
    return (diam_m * larg_pol * 72.0 / faixa_m) ** 2


def _txt_k(ks, cov):
    return (rf"$\mu = {_sci_mat(ks / EC.GAMMA_W)[1:-1]}\ \mathrm{{m^2\,kPa^{{-1}}\,s^{{-1}}}}\quad "
            rf"CoV = {_virg_mat(100 * cov, 0)}\%$")


def _txt_c(cov):
    return rf"$\mu = {_virg_mat(MEDIA_C)}\ \mathrm{{kPa}}\quad CoV = {_virg_mat(100 * cov, 0)}\%$"


def _txt_phi(cov):
    return rf"$\mu = {_virg_mat(MEDIA_PHI)}^\circ\quad CoV = {_virg_mat(100 * cov, 0)}\%$"


def _salva(fig, saida_fig, nome):
    fig.savefig(os.path.join(saida_fig, nome + ".pdf"), bbox_inches="tight", pad_inches=0.02)
    fig.savefig(os.path.join(saida_fig, nome + ".png"), dpi=150, bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)


# ---------------------------------------------------------------------------------------------------------------------
# dados de um solo


def _info_solo(solo, est, raiz):
    """caso de referencia, pastas, deterministico e banda_ymin de todos os casos de um solo"""
    if est is not None:
        se = est["solos"][solo]
        casos = {}
        for cid, c in se["casos"].items():
            det = c.get("det")
            a = c.get("amostras") or {}
            casos[cid] = {"pastas": list(c["pastas"]), "fs_def": c["fs_def"],
                          "det": det.get("pasta") if det else None, "det_est": det,
                          "real": list(a.get("real") or []),
                          "banda_ymin": [_float(v) for v in (a.get("banda_ymin") or [])]}
        return {"ref": se["ref"], "casos": casos, "k_s": float(se.get("k_s_m_s", K_S[solo])),
                "fonte": "estatisticas.json (amostras.banda_ymin)"}
    mapa = _mapa_padrao()[solo]
    casos = {}
    for cid, c in mapa["casos"].items():
        reais, ymin = _amostras_resumo(raiz, c["pastas"])
        casos[cid] = {"pastas": c["pastas"], "fs_def": c["fs_def"], "det": c["det"], "det_est": None,
                      "real": reais, "banda_ymin": ymin}
    return {"ref": mapa["ref"], "casos": casos, "k_s": K_S[solo], "fonte": "resumo.csv (banda_ymin)"}


def _amostras_resumo(raiz, pastas):
    """modo autonomo: numero e banda_ymin (de resumo.csv) das realizacoes concluidas das pastas, unidas so depois de
    conferir parametros.txt identicos e numeros de realizacao disjuntos"""
    pastas = [_caminho(raiz, p) for p in pastas]
    par = []
    for p in pastas:
        with open(os.path.join(p, "parametros.txt"), "rb") as f:
            par.append(f.read())
    if any(x != par[0] for x in par[1:]):
        raise RuntimeError(f"parametros.txt diferentes nas pastas a unir: {pastas}")
    vistos = {}
    for p in pastas:
        for n, pr in _realizacoes(p):
            if n in vistos:
                raise RuntimeError(f"realizacao {n} em {vistos[n]} e em {pr}: numeros nao disjuntos")
            vistos[n] = pr
    reais = sorted(vistos)
    return reais, [_float(_resumo(vistos[n]).get("banda_ymin")) for n in reais]


def _modo(ymin, delta):
    if not math.isfinite(ymin):
        return None
    if ymin > Y_PE + delta:
        return "acima"
    if ymin < Y_PE - delta:
        return "abaixo"
    return "pe"


def _porcentagens(ymins, delta):
    m = [_modo(y, delta) for y in ymins]
    n = sum(1 for x in m if x is not None)
    cont = {k: sum(1 for x in m if x == k) for k in MODOS}
    out = {k: (round(100.0 * cont[k] / n, 6) if n else None) for k in MODOS}
    out.update({f"n_{k}": cont[k] for k in MODOS})
    out["N"] = n
    out["n_sem_banda"] = len(m) - n
    return out


def _versao(res):
    return (res.get("versao") or "?").split("_")[0]


def _carrega_solo(solo, est, raiz):
    """le os dados das figuras de um solo"""
    info = _info_solo(solo, est, raiz)
    ref = info["casos"][info["ref"]]
    fs_def = ref["fs_def"]
    avisos = []
    # realizacao tipica: r00000 do caso de referencia (a pasta que tem campos.csv)
    tipica = None
    for p in ref["pastas"]:
        r0 = os.path.join(_caminho(raiz, p), "r00000")
        if _concluida(r0) and os.path.isfile(os.path.join(r0, "campos.csv")):
            tipica = r0
            break
    if tipica is None:
        raise RuntimeError(f"{solo}: r00000 com campos.csv nao encontrada nas pastas de referencia {ref['pastas']}")
    res0 = _resumo(tipica)
    t0, fs0 = _instante(fs_def, res0)
    malha = Malha(tipica)
    cam = _le_num(os.path.join(tipica, "campos.csv"))
    malha.prepara_campos(cam["x"], cam["y"])
    pw = malha.pw
    s = _nas(pw["t_h"], t0)
    if s.sum() != malha.X.size:
        raise RuntimeError(f"{tipica}/pw.csv sem os {malha.X.size} nos em t = {t0:g} h")
    pw0 = np.empty(malha.X.size)
    pw0[pw["no"][s].astype(int)] = pw["pw_kPa"][s]
    banda0 = _banda(tipica, t0)
    if abs(banda0["FS"] - fs0) > TOL_FS:
        raise RuntimeError(f"{tipica}: FS de mecanismo.csv ({banda0['FS']}) diferente do resumo ({fs0})")
    # solucao deterministica do caso de referencia, no instante que define o FS do caso
    if not ref.get("det"):
        raise RuntimeError(f"{solo}: caso de referencia sem pasta deterministica")
    det_dir = _pasta_det(raiz, ref["det"])
    rdet = _resumo(det_dir)
    tdet, fsdet = _instante(fs_def, rdet)
    de = ref.get("det_est")
    if de and (abs(_float(de.get("FS")) - fsdet) > 1e-9 or abs(_float(de.get("t_h")) - tdet) > 1e-9):
        avisos.append(f"{solo}: FS e t deterministicos do JSON ({de.get('FS')}, {de.get('t_h')}) diferentes dos de "
                      f"resumo.csv ({fsdet}, {tdet}); as figuras usam resumo.csv")
    banda_det = _banda(det_dir, tdet)
    if abs(banda_det["FS"] - fsdet) > TOL_FS:
        raise RuntimeError(f"{det_dir}: FS de mecanismo.csv ({banda_det['FS']}) diferente do resumo ({fsdet})")
    # bandas de todas as realizacoes da pasta de referencia que tem mecanismo.csv em todas
    desenho = None
    for p in ref["pastas"]:
        reais = _realizacoes(_caminho(raiz, p))
        if reais and all(os.path.isfile(os.path.join(pr, "mecanismo.csv")) for _, pr in reais):
            desenho = (p, reais)
            break
    if desenho is None:
        raise RuntimeError(f"{solo}: nenhuma pasta de referencia com mecanismo.csv em todas as realizacoes")
    p_des, reais = desenho
    contagem, coord = {}, {}
    ymins, dif, versoes, insts = [], [], {}, {}
    for n, pr in reais:
        r = _resumo(pr)
        t, fs = _instante(fs_def, r)
        b = _banda(pr, t)
        if abs(b["FS"] - fs) > TOL_FS:
            raise RuntimeError(f"{pr}: FS de mecanismo.csv ({b['FS']}) diferente do resumo ({fs})")
        for ip, x, y in zip(b["ip"].tolist(), b["x"].tolist(), b["y"].tolist()):
            contagem[ip] = contagem.get(ip, 0) + 1
            if ip not in coord:
                coord[ip] = (x, y)
            elif abs(coord[ip][0] - x) > 1e-3 or abs(coord[ip][1] - y) > 1e-3:
                raise RuntimeError(f"ponto de integracao {ip} com coordenadas diferentes entre realizacoes")
        ym = float(b["y"].min())
        ymins.append(ym)
        # banda_ymin de resumo.csv e a do instante do FS minimo; noutro instante, a comparacao e com fs.csv
        if abs(t - _float(r["t_min_h"])) < 1e-9:
            yref = _float(r["banda_ymin"])
        else:
            f = _le_num(os.path.join(pr, "fs.csv"))
            yref = float(f["banda_ymin"][_nas(f["t_h"], t)][0])
        dif.append(abs(ym - yref))
        v = _versao(r)
        versoes[v] = versoes.get(v, 0) + 1
        insts[t] = insts.get(t, 0) + 1
    ips = sorted(contagem)
    sobre = {"x": np.array([coord[i][0] for i in ips]), "y": np.array([coord[i][1] for i in ips]),
             "n": np.array([contagem[i] for i in ips])}
    dif = np.array(dif)
    return {
        "solo": solo, "info": info, "fs_def": fs_def, "tipica": tipica, "res0": res0, "t0": t0, "fs0": fs0,
        "malha": malha, "campos": cam, "pw0": pw0, "banda0": banda0,
        "det_dir": det_dir, "rdet": rdet, "tdet": tdet, "fsdet": fsdet, "banda_det": banda_det,
        "p_desenho": p_des, "reais_desenho": [n for n, _ in reais], "sobre": sobre, "ymins": ymins,
        "verif": {"n": int(dif.size), "iguais": int(np.sum(dif <= 1e-4)), "dif_max_m": float(dif.max()),
                  "tolerancia_m": 1e-4},
        "versoes": versoes, "instantes": insts, "avisos": avisos,
        "N_ref": len(ref["real"]),
    }


# ---------------------------------------------------------------------------------------------------------------------
# figuras


def _desc_ref(d):
    t = _virg(d["t0"])
    if d["fs_def"] == "FS_min":
        return f"solo {d['solo']}, I = 36 mm/h; instante do FS mínimo da realização, t = {t} h"
    return f"solo {d['solo']}, I = 36 mm/h; FS em t = {t} h"


def _fig12(d, saida_fig, notas):
    """Fig. 12 do Ceron: (a) k_v, (b) c', (c) phi' da realizacao tipica, (a) e (b) em cima, (c) embaixo a esquerda"""
    solo, m, cam, res0 = d["solo"], d["malha"], d["campos"], d["res0"]
    ks = d["info"]["k_s"]
    kv = cam["k_mult"] * ks / EC.GAMMA_W
    larg, esq2, alt_linha = 3.1, 3.5, 1.95
    fig = plt.figure(figsize=(esq2 + larg, 2 * alt_linha))
    paineis = [("k", kv, "(a)", 0.0, 0.0, _txt_k(ks, _float(res0["cov_k"]))),
               ("c", cam["c"], "(b)", esq2, 0.0, _txt_c(_float(res0["cov_c"]))),
               ("phi", cam["phi_deg"], "(c)", 0.0, alt_linha, _txt_phi(_float(res0["cov_phi"])))]
    lim_k = None
    for chave, z, letra, esq, topo, txt in paineis:
        ax, cax = _eixos_painel(fig, esq, topo, larg)
        lo, hi, marcas, nota = _barra(chave, float(z.min()), float(z.max()))
        _desenha_campo(fig, ax, cax, m.trip, m.valores_campo(z), chave, lo, hi, marcas, EC.CAMPOS[chave]["rot"],
                       letra)
        _texto_estat(ax, txt)
        if nota:
            notas.append(f"{letra} {nota}.")
        if chave == "k":
            lim_k = (lo, hi)
    _salva(fig, saida_fig, f"{solo}_fig12")
    v = m.verif_campos
    notas[:0] = [
        f"Realização 0 do caso de referência ({_desc_ref(d)}). Os campos estão nos pontos de integração (3 por "
        f"triângulo da malha); desenho por interpolação linear numa triangulação de Delaunay desses pontos e dos nós "
        f"do contorno (valor do ponto de integração mais próximo), restrita ao domínio "
        f"({_virg(v['area_m2'], 0)} m²).",
        f"(a) Permeabilidade vertical na unidade do Ceron (m² kPa-1 s-1): k (m/s) dividida por 9,81 kN/m³. Nesta "
        f"realização, de "
        f"{_sci_txt(float(kv.min()))} a {_sci_txt(float(kv.max()))}; barra de cores de {_sci_txt(lim_k[0])} a "
        f"{_sci_txt(lim_k[1])}.",
        f"Média e CoV escritos nos painéis: os dos campos lognormais do caso de referência (coesão "
        f"{_virg(MEDIA_C)} kPa, ângulo de atrito {_virg(MEDIA_PHI)}°, permeabilidade igual a ks = "
        f"{_sci_txt(ks)} m/s).",
        f"Coesão desta realização de {_virg(float(cam['c'].min()), 2)} a {_virg(float(cam['c'].max()), 2)} kPa; "
        f"ângulo de atrito de {_virg(float(cam['phi_deg'].min()), 2)}° a {_virg(float(cam['phi_deg'].max()), 2)}°."]


def _fig13(d, saida_fig, notas):
    """Fig. 13 do Ceron: a poropressao (aqui p_w) da realizacao tipica"""
    solo, m, pw0 = d["solo"], d["malha"], d["pw0"]
    larg = 3.5
    fig = plt.figure(figsize=(larg, larg * ALT_DOM / LARG_DOM + 0.5))
    ax, cax = _eixos_painel(fig, 0.0, 0.0, larg)
    lo, hi, marcas, nota = _barra("u", float(pw0.min()), float(pw0.max()))
    _desenha_campo(fig, ax, cax, m.tri, pw0, "u", lo, hi, marcas, ROTULO_PW)
    _salva(fig, saida_fig, f"{solo}_fig13")
    notas.append(f"Poropressão da água pw nos nós da malha, na realização 0 do caso de referência ({_desc_ref(d)}); "
                 f"de {_virg(float(pw0.min()), 1)} a {_virg(float(pw0.max()), 1)} kPa. O Ceron mostra o excesso de "
                 f"poropressão do fluxo permanente (de -50 a 0 kPa na barra): as grandezas são diferentes.")
    if nota:
        notas.append(nota + ".")


def _fig14(d, saida_fig, notas):
    """Fig. 14 do Ceron: p_w, c' e phi' com a superficie (banda) critica da realizacao tipica e o seu FS"""
    solo, m, cam, b, res0 = d["solo"], d["malha"], d["campos"], d["banda0"], d["res0"]
    larg, esp = 2.15, 0.26  # o espaco evita que as marcacoes de barras vizinhas se encostem
    fig = plt.figure(figsize=(3 * larg + 2 * esp, larg * ALT_DOM / LARG_DOM + 0.45))
    nomes = {"u": "poropressão", "c": "coesão", "phi": "ângulo de atrito"}
    paineis = (("u", d["pw0"], m.tri, None), ("c", cam["c"], m.trip, _txt_c(_float(res0["cov_c"]))),
               ("phi", cam["phi_deg"], m.trip, _txt_phi(_float(res0["cov_phi"]))))
    for i, (chave, z, tri, txt) in enumerate(paineis):
        ax, cax = _eixos_painel(fig, i * (larg + esp), 0.0, larg)
        zz = z if tri is m.tri else m.valores_campo(z)
        lo, hi, marcas, nota = _barra(chave, float(z.min()), float(z.max()))
        rot = ROTULO_PW if chave == "u" else EC.CAMPOS[chave]["rot"]
        _desenha_campo(fig, ax, cax, tri, zz, chave, lo, hi, marcas, rot)
        ax.scatter(b["x"], b["y"], s=_s_marcador(0.4, larg, LARG_DOM), c="black", linewidths=0, zorder=4)
        if txt:
            _texto_estat(ax, txt, tam=5.5)
        if chave == "phi":  # rotulo da superficie critica no vazio acima da berma, como na Fig. 14 do Ceron
            ax.annotate(f"Superfície crítica\n$FS = {_virg_mat(d['fs0'])}$", xy=_ponto_rotulo(b),
                        xytext=(56.0, 26.0), ha="center", va="center", fontsize=6.5,
                        arrowprops=dict(arrowstyle="-", color="black", lw=0.5, shrinkA=1, shrinkB=1))
        if nota:
            notas.append(f"({nomes[chave]}) {nota}.")
    _salva(fig, saida_fig, f"{solo}_fig14")
    notas[:0] = [f"Realização 0 do caso de referência ({_desc_ref(d)}); FS = {_fs_txt(d['fs0'])}. Superfície "
                 f"crítica: banda de ruptura, os {b['x'].size} pontos de integração com raiz de J2 da deformação "
                 f"plástica maior ou igual a 0,1 vez o máximo no último estado convergido do SRM nesse instante. "
                 f"Poropressão pw no lugar do excesso de poropressão do Ceron."]


def _ponto_rotulo(b):
    """o ponto da banda mais a direita e acima entre 24 m e 28,5 m de altura, para a linha do rotulo"""
    x, y = b["x"], b["y"]
    sel = (y > 24.0) & (y < 28.5)
    if not sel.any():
        sel = np.ones_like(x, dtype=bool)
    i = int(np.argmax(np.where(sel, x + y, -np.inf)))
    return float(x[i]), float(y[i])


def _limites_fig16(ds):
    """o mesmo recorte para os dois solos: as bandas, o pe e espaco para a legenda e as porcentagens"""
    xs = np.concatenate([np.concatenate([d["sobre"]["x"], d["banda_det"]["x"]]) for d in ds])
    ys = np.concatenate([np.concatenate([d["sobre"]["y"], d["banda_det"]["y"]]) for d in ds])
    x0 = max(min(math.floor(xs.min()) - 3.0, 15.0), 0.0)
    x1 = min(max(math.ceil(xs.max()) + 3.0, 66.0), LARG_DOM)
    y0 = max(min(math.floor(ys.min()) - 1.5, 12.0), 0.0)
    return x0, x1, y0, 36.5


def _fig16(d, saida_fig, notas, lims, delta, modos_json):
    """Figs. 16 e 20 do Ceron: bandas de todas as realizacoes, a deterministica e os modos de ruptura"""
    solo, m = d["solo"], d["malha"]
    x0, x1, y0, y1 = lims
    w = 6.6
    fig = plt.figure(figsize=(w, w * (y1 - y0) / (x1 - x0)))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)
    ax.set_aspect("equal")
    ax.set_axis_off()
    s = _s_marcador(DIAM_PONTO_16, w, x1 - x0)
    sb = d["sobre"]
    cor = np.zeros((sb["n"].size, 4))
    cor[:, :3] = to_rgb(CINZA_BANDAS)
    cor[:, 3] = 1.0 - (1.0 - ALFA_SOBREPOSICAO) ** sb["n"]  # n pontos sobrepostos de opacidade ALFA_SOBREPOSICAO
    o = np.argsort(cor[:, 3])
    ax.scatter(sb["x"][o], sb["y"][o], s=s, c=cor[o], linewidths=0, zorder=2)
    bd = d["banda_det"]
    # a deterministica com pontos menores, para o cinza das realizacoes continuar visivel em volta
    ax.scatter(bd["x"], bd["y"], s=_s_marcador(DIAM_DET_16, w, x1 - x0), c=VERMELHO, linewidths=0, zorder=3)
    xs, ys = m.superficie
    ax.plot(xs, ys, color=LINHA_TERRENO, lw=1.3, zorder=4, solid_capstyle="butt")
    leg = [Line2D([0], [0], marker="o", ls="none", color=VERMELHO, markersize=4.5, label="Solução determinística")]
    ax.legend(handles=leg, loc="upper left", bbox_to_anchor=(0.01, 0.985), frameon=True, fancybox=True,
              edgecolor="0.55", framealpha=1.0, handletextpad=0.4, borderpad=0.6)
    pc = modos_json[solo]
    xa, xb, ytop = 0.66, 0.97, 0.86
    ax.text((xa + xb) / 2, ytop, "Porcentagem do modo de ruptura", transform=ax.transAxes, ha="center",
            va="bottom", fontsize=8.5)
    ax.text(xa, ytop - 0.015, "Acima do pé do talude:\nPelo pé do talude:\nAbaixo do pé do talude:",
            transform=ax.transAxes, ha="left", va="top", fontsize=8.5, linespacing=1.3)
    ax.text(xb, ytop - 0.015, "\n".join(f"{_virg(pc[k], 2)}%" for k in MODOS), transform=ax.transAxes,
            ha="right", va="top", fontsize=8.5, linespacing=1.3)
    ax.text((xa + xb) / 2, ytop - 0.235, f"({pc['N']} realizações)", transform=ax.transAxes, ha="center",
            va="top", fontsize=7.5)
    _salva(fig, saida_fig, f"{solo}_fig16")
    N_ref, nd = d["N_ref"], len(d["reais_desenho"])
    r0, r1 = min(d["reais_desenho"]), max(d["reais_desenho"])
    if nd < N_ref:
        quem = (f"das {nd} realizações de número {r0} a {r1} do caso de referência, as únicas das {N_ref} com "
                f"mecanismo.csv")
    else:
        quem = f"das {nd} realizações do caso de referência (números {r0} a {r1})"
    t = "; ".join(f"{v} {'realizações' if v > 1 else 'realização'} em t = {_virg(k)} h"
                  for k, v in sorted(d["instantes"].items()))
    notas.append(f"Bandas de ruptura {quem}, no instante do FS de cada uma ({t}), em cinza translúcido (opacidade "
                 f"{_virg(ALFA_SOBREPOSICAO)} por realização em cada ponto de integração, acumulada), e banda da "
                 f"solução determinística (médias das propriedades; t = {_virg(d['tdet'])} h, "
                 f"FS = {_fs_txt(d['fsdet'])}) em vermelho. Banda de ruptura: pontos de integração com raiz de J2 da "
                 f"deformação plástica maior ou igual a 0,1 vez o máximo no último estado convergido do SRM.")
    notas.append(f"Modo de ruptura pelo ponto mais baixo da banda (y min) em relação à cota do pé (20 m), com "
                 f"tolerância de {_virg(delta)} m (mediana das arestas dos triângulos da malha com vértice a até 1 m "
                 f"do pé): acima do pé se y min maior que {_virg(Y_PE + delta)} m, abaixo do pé se y min menor que "
                 f"{_virg(Y_PE - delta)} m, pelo pé nos demais casos. Contagens: acima do pé {pc['n_acima']}, pelo "
                 f"pé {pc['n_pe']}, abaixo do pé {pc['n_abaixo']}. Solução determinística: y min = "
                 f"{_virg(float(bd['y'].min()), 2)} m ({NOME_MODO[pc['det']['modo']]}).")
    v = d["verif"]
    notas.append(f"Em {v['iguais']} de {v['n']} realizações, o y min da banda reconstruída de mecanismo.csv coincide "
                 f"com o de resumo.csv dentro de {_virg(v['tolerancia_m'])} m (diferença máxima: "
                 f"{_virg(v['dif_max_m'], 6)} m).")
    cr = modos_json[f"{solo}_casos"].get(d["info"]["ref"])
    if cr and cr["N"] and cr["N"] != pc["N"]:
        notas.append(f"Com o y min da banda de resumo.csv, nas {cr['N']} realizações do caso de referência: acima do "
                     f"pé {cr['n_acima']}, pelo pé {cr['n_pe']}, abaixo do pé {cr['n_abaixo']}.")
    vd = _versao(d["rdet"])
    if len(d["versoes"]) != 1 or vd not in d["versoes"]:
        vs = ", ".join(f"{k} ({c})" for k, c in sorted(d["versoes"].items()))
        notas.append(f"Versão do programa: rodada determinística {vd}; realizações desenhadas {vs}.")
    notas.append("O SRM não define centro de rotação: a nuvem de centros das Figs. 16 e 20 do Ceron não tem "
                 "equivalente aqui.")


def _fig24(d, saida_fig, notas):
    """Fig. 24 do Ceron: geometria do talude e superficie critica da solucao deterministica"""
    solo, m, bd = d["solo"], d["malha"], d["banda_det"]
    x0, x1, y0, y1 = -1.0, 71.0, -1.0, 36.0
    w = 5.4
    fig = plt.figure(figsize=(w, w * (y1 - y0) / (x1 - x0)))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)
    ax.set_aspect("equal")
    ax.set_axis_off()
    laco = m.lacos[0]
    ax.fill(m.X[laco], m.Y[laco], facecolor="0.93", edgecolor="none", zorder=0)
    ax.plot(np.append(m.X[laco], m.X[laco[0]]), np.append(m.Y[laco], m.Y[laco[0]]), color="black", lw=0.8,
            zorder=3)
    # chuva: setas verticais sobre a superficie, fora das cotas de H e de 45 graus
    xs, ys = m.superficie
    for xa in (2.5, 7.5, 12.5, 17.5, 22.5, 27.5, 36.0, 51.0, 56.0, 61.0, 66.0):
        ya = float(np.interp(xa, xs, ys))
        ax.annotate("", xy=(xa, ya + 0.4), xytext=(xa, ya + 3.4),
                    arrowprops=dict(arrowstyle="-|>", color=AZUL_CHUVA, lw=0.7, mutation_scale=5, shrinkA=0,
                                    shrinkB=0))
    I = _float(d["res0"]["I_mm_h"])
    ax.text(15.0, 34.2, f"Chuva: $I = {_virg_mat(I)}$ mm/h durante 24 h", ha="center", va="bottom", fontsize=7.5,
            color=AZUL_CHUVA)
    # H: cotas da crista e do pe, com seta dupla a direita do pe
    xh = 45.0
    ax.plot([X_CRISTA, xh + 1.0], [Y_CRISTA, Y_CRISTA], ls=(0, (3, 2)), color="0.35", lw=0.6, zorder=2)
    ax.annotate("", xy=(xh, Y_PE), xytext=(xh, Y_CRISTA),
                arrowprops=dict(arrowstyle="<|-|>", color="black", lw=0.7, mutation_scale=6, shrinkA=0, shrinkB=0))
    ax.text(xh + 0.8, (Y_PE + Y_CRISTA) / 2, rf"$H = {_virg_mat(H_TALUDE)}$ m", ha="left", va="center", fontsize=8)
    # 45 graus na crista, entre a horizontal e a face (como na Fig. 24 do Ceron)
    r = 3.0
    ang = np.radians(np.linspace(-BETA_GRAUS, 0.0, 30))
    ax.plot(X_CRISTA + r * np.cos(ang), Y_CRISTA + r * np.sin(ang), color="black", lw=0.6, zorder=3)
    ax.text(X_CRISTA + r + 0.5, Y_CRISTA - 1.2, rf"${BETA_GRAUS:.0f}^\circ$", ha="left", va="center", fontsize=8)
    ax.scatter(bd["x"], bd["y"], s=_s_marcador(0.45, w, x1 - x0), c="black", linewidths=0, zorder=4)
    ax.annotate(f"Superfície crítica, $FS = {_virg_mat(d['fsdet'])}$", xy=_ponto_meio(bd), xytext=(14.0, 21.5),
                ha="center", va="center", fontsize=7.5,
                arrowprops=dict(arrowstyle="-", color="black", lw=0.5, shrinkA=1, shrinkB=1))
    _salva(fig, saida_fig, f"{solo}_fig24")
    notas.append(f"Talude de altura H = {_virg(H_TALUDE)} m e inclinação de {BETA_GRAUS:.0f}° (crista em "
                 f"(30 m, 30 m), pé em (40 m, 20 m)), no domínio de {_virg(float(np.ptp(m.X)))} m por "
                 f"{_virg(float(np.ptp(m.Y)))} m ({_virg(m.area, 0)} m²); chuva de {_virg(I)} mm/h durante 24 h na "
                 f"superfície do terreno.")
    notas.append(f"Superfície crítica: banda de ruptura da solução determinística (médias das propriedades, "
                 f"ks = {_sci_txt(d['info']['k_s'])} m/s) em t = {_virg(d['tdet'])} h, FS = {_fs_txt(d['fsdet'])}: "
                 f"os {bd['x'].size} pontos de integração com raiz de J2 da deformação plástica maior ou igual a 0,1 "
                 f"vez o máximo no último estado convergido do SRM.")
    vd = _versao(d["rdet"])
    if len(d["versoes"]) != 1 or vd not in d["versoes"]:
        notas.append(f"Versão do programa da rodada determinística: {vd}; das realizações do caso de referência: "
                     f"{', '.join(sorted(d['versoes']))}.")


def _ponto_meio(b):
    """o ponto da banda mais a esquerda perto da metade da sua altura, para a linha do rotulo"""
    x, y = b["x"], b["y"]
    sel = np.abs(y - 0.5 * (y.min() + y.max())) < 1.0
    if not sel.any():
        sel = np.ones_like(x, dtype=bool)
    i = int(np.argmin(np.where(sel, x, np.inf)))
    return float(x[i]), float(y[i])


# ---------------------------------------------------------------------------------------------------------------------


def gera(est, raiz_dados, saida):
    """grava em <saida>/figuras as Figs. 12, 13, 14, 16 e 24 de cada solo (PDF e PNG), legendas_campos.json e
    modos_ruptura.json; devolve {"figuras": [...], "avisos": [...]}"""
    if isinstance(est, (str, os.PathLike)):
        with open(est) as f:
            est = json.load(f)
    saida_fig = os.path.join(saida, "figuras")
    os.makedirs(saida_fig, exist_ok=True)
    solos = [s for s in ("f50", "f100") if est is None or s in est.get("solos", {})]
    ds = [_carrega_solo(s, est, raiz_dados) for s in solos]
    delta, ntri, nar = ds[0]["malha"].tolerancia_pe()
    for d in ds[1:]:
        if abs(d["malha"].tolerancia_pe()[0] - delta) > 1e-12:
            raise RuntimeError("malhas diferentes entre os solos")
    modos = {}
    for d in ds:
        solo = d["solo"]
        pc = _porcentagens(d["ymins"], delta)
        ydet = float(d["banda_det"]["y"].min())
        pc.update({"delta_m": delta, "regra": REGRA_MODO, "y_pe_m": Y_PE, "triangulos_perto_do_pe": ntri,
                   "arestas_perto_do_pe": nar, "pasta": d["p_desenho"],
                   "realizacoes": [min(d["reais_desenho"]), max(d["reais_desenho"])], "N_caso": d["N_ref"],
                   "fonte": "y min da banda reconstruida de mecanismo.csv no instante do FS de cada realizacao",
                   "det": {"modo": _modo(ydet, delta), "banda_ymin": ydet, "FS": d["fsdet"], "t_h": d["tdet"],
                           "pasta": os.path.relpath(d["det_dir"], raiz_dados)},
                   "verificacao_banda_ymin": d["verif"]})
        modos[solo] = pc
        casos = {}
        for cid, c in d["info"]["casos"].items():
            casos[cid] = _porcentagens(c["banda_ymin"], delta)
            casos[cid]["fonte"] = d["info"]["fonte"]
        modos[f"{solo}_casos"] = casos
    lims = _limites_fig16(ds)
    legendas, avisos = {}, []
    with plt.rc_context(ESTILO):
        for d in ds:
            solo = d["solo"]
            avisos += d["avisos"]
            for nome, f in (("fig12", _fig12), ("fig13", _fig13), ("fig14", _fig14), ("fig24", _fig24)):
                notas = []
                f(d, saida_fig, notas)
                legendas[f"{solo}_{nome}"] = notas
            notas = []
            _fig16(d, saida_fig, notas, lims, delta, modos)
            legendas[f"{solo}_fig16"] = notas
    legendas = {k: legendas[k] for k in sorted(legendas, key=lambda k: (k.split("_")[0] != "f50", k))}
    with open(os.path.join(saida_fig, "legendas_campos.json"), "w", encoding="utf-8") as f:
        json.dump(legendas, f, ensure_ascii=False, indent=1)
    with open(os.path.join(saida_fig, "modos_ruptura.json"), "w", encoding="utf-8") as f:
        json.dump(modos, f, ensure_ascii=False, indent=1)
    for a in avisos:
        print("AVISO:", a, file=sys.stderr)
    figs = [os.path.join(saida_fig, f"{d['solo']}_{n}.pdf") for d in ds
            for n in ("fig12", "fig13", "fig14", "fig16", "fig24")]
    return {"figuras": figs, "avisos": avisos}


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    if len(args) != 2:
        print(__doc__)
        sys.exit(1)
    raiz, saida = args
    j = os.path.join(saida, "estatisticas.json")
    est = j if "--sem-est" not in sys.argv and os.path.isfile(j) else None
    print("estatisticas:", est or "nenhuma (le os resumo.csv)")
    r = gera(est, raiz, saida)
    print("\n".join(r["figuras"]))
