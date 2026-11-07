#!/usr/bin/env python3
"""Animacao de um caso de talude sob chuva (Rahardjo et al. 2007): campos hidraulicos, estado do SRM e FS(t).

Uso: python3 animacao.py <pasta do caso> [--passo N] [--fps F] [--campos a,b,...] [--largura-gif PX] [--u-total]
                         [--manter-quadros]
     --passo N        usa um a cada N passos gravados (padrao: o que der cerca de 100 quadros); todo FS vira um quadro
     --fps F          quadros por segundo (padrao 4)
     --campos         paineis, na ordem pedida (padrao: todos): pw, succao, saturacao, ar, u, fs, j2, hardening, falha
     --largura-gif PX largura do GIF em pixels (padrao 1400; 0: nao grava o GIF)
     --u-total        desenha o deslocamento total, com o recalque de ligar a gravidade (padrao: desde o estado geostatico)
     --manter-quadros guarda os PNG dos quadros em <pasta do caso>/quadros

Le da pasta do caso, gravados pelo main.cpp:
  vertices.csv    no,x,y
  triangles.csv   n0,n1,n2 (numeros de no)
  pontos.csv      tri,x,y,peso: cada ponto de integracao, com a linha do seu triangulo em triangles.csv
  snaps_pw.csv, snaps_pa.csv, snaps_sw.csv
                  em cada linha: t (h) e o valor nos vertices, na ordem de vertices.csv
  fs.csv          t_h,FS
  srm_nos.csv     t_h,no,ux,uy,ux_geo,uy_geo: ultimo estado convergido do SRM e estado geostatico, nos vertices
  srm_pontos.csv  t_h,ponto,raizJ2_epsp,hardening,tipo: ultimo estado convergido do SRM, nos pontos de integracao
  summary.txt     "I <mm/h>" e "train <h>", um por linha

Grava na mesma pasta: fs_tempo.png, chuva.mp4, chuva.gif.

O transiente tem esqueleto rigido (u = 0, sem plasticidade): deslocamento, deformacao plastica, hardening e tipo de
falha so existem nos estados do SRM, entao esses paineis mudam apenas nos instantes em que o FS foi calculado. O ultimo
estado convergido do SRM esta num F proximo do FS, a uma distancia que varia de um instante para outro, e o pico junto
ao pe depende da malha: por isso u e normalizado pelo seu maximo em cada estado, e sqrt(J2) e hardening usam escala
logaritmica fixa. O painel de u mostra, por padrao, o deslocamento desde o estado geostatico (o mecanismo): o SRM parte
de uma malha sem tensoes, e o recalque de ligar a gravidade esconderia o mecanismo nos instantes com |u| pequeno.
"""
import argparse
import math
import os
import shutil
import subprocess
import tempfile
from types import SimpleNamespace

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np
from matplotlib.cm import ScalarMappable
from matplotlib.collections import PolyCollection
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm, from_levels_and_colors
from matplotlib.patches import Patch
from matplotlib.ticker import FuncFormatter

# ---------------------------------- aparencia ----------------------------------
FUNDO = "#fcfcfb"
TINTA = "#0b0b0b"
TINTA2 = "#52514e"
APAGADO = "#898781"
GRADE = "#e1e0d9"
EIXO = "#c3c2b7"
NEUTRO = "#f0efec"  # meio da escala divergente; tambem "nada aqui": saturado, sem plastificacao, elastico
MAGENTA = "#d6246e"  # p_w = 0: o lencol e a face de percolacao
AZUL_LINHA = "#2a78d6"
CHUVA = "#104281"
PASSO_ISOLINHA = 25.0  # kPa
PASSO_COR = 5.0  # kPa: a cor varia mais fino que as isolinhas
LOG_MIN = 1.0e-3  # sqrt(J2) e hardening abaixo disso ficam sem cor

# Divergente com o meio cinza em p = 0: vermelho para p < 0 (succao) e azul para p > 0 (agua sob pressao)
DIVERGENTE = LinearSegmentedColormap.from_list(
    "pressao",
    [(0.00, "#7f1d1d"), (0.25, "#c8393b"), (0.42, "#eda18f"), (0.50, NEUTRO),
     (0.60, "#9ec5f4"), (0.78, "#3987e5"), (1.00, "#104281")],
)

# Tipo de resposta do Mohr-Coulomb PV2 (m_m_type): cada tipo tem sempre a mesma cor. Os tres tons plasticos passam na
# validacao de daltonismo para qualquer par (cores 1 a 3 da paleta); os dois cantos sao o mesmo laranja em dois tons.
TIPOS = {
    0: ("elástico", NEUTRO),
    1: ("plano principal", "#2a78d6"),
    3: ("canto esquerdo (β = π/3)", "#eb6834"),
    2: ("canto direito (β = 0)", "#a8431c"),
    -1: ("ápice (tração)", "#1baf7a"),
}


def virgula(v, casas=None):
    """Numero em portugues: virgula decimal e sinal de menos tipografico."""
    texto = f"{v:g}" if casas is None else f"{v:.{casas}f}"
    return texto.replace(".", ",").replace("-", "−")


FORMATO = FuncFormatter(lambda v, _: virgula(v))


def rotulo_tempo(t):
    return virgula(t, 2 if t < 1 else 1) + " h"


# ---------------------------------- dados ----------------------------------
def contorno_do_dominio(tri):
    """Indices dos vertices do contorno em ordem: as arestas que pertencem a um so triangulo, encadeadas."""
    contagem = {}
    for a, b, c in tri:
        for u, v in ((a, b), (b, c), (c, a)):
            chave = (min(u, v), max(u, v))
            contagem[chave] = contagem.get(chave, 0) + 1
    vizinhos = {}
    for (u, v), n in contagem.items():
        if n == 1:
            vizinhos.setdefault(u, []).append(v)
            vizinhos.setdefault(v, []).append(u)
    inicio = next(iter(vizinhos))
    caminho, anterior, atual = [inicio], None, inicio
    while True:
        proximo = next(n for n in vizinhos[atual] if n != anterior)
        if proximo == inicio:
            return caminho + [inicio]
        caminho.append(proximo)
        anterior, atual = atual, proximo


def le_caso(pasta):
    arq = lambda nome: os.path.join(pasta, nome)
    le = lambda nome, pula=1: np.loadtxt(arq(nome), delimiter=",", skiprows=pula, ndmin=2)

    necessarios = ("vertices.csv", "triangles.csv", "pontos.csv", "snaps_pw.csv", "snaps_pa.csv", "snaps_sw.csv", "fs.csv",
                   "srm_nos.csv", "srm_pontos.csv", "summary.txt")
    faltando = [n for n in necessarios if not os.path.exists(arq(n))]
    if faltando:
        raise SystemExit(f"faltam em {pasta}: {', '.join(faltando)}. Rode o main.cpp atual para gravar esses arquivos.")

    c = SimpleNamespace()
    vert = le("vertices.csv")
    c.nos, c.xy = vert[:, 0].astype(int), vert[:, 1:3]
    c.tri = np.searchsorted(c.nos, le("triangles.csv").astype(int))
    c.malha = mtri.Triangulation(c.xy[:, 0], c.xy[:, 1], c.tri)
    c.localiza = c.malha.get_trifinder()
    c.borda = contorno_do_dominio(c.tri)
    c.contorno = c.xy[c.borda]
    c.L, c.H = c.xy[:, 0].max(), c.xy[:, 1].max()
    c.x0, c.y0 = c.xy[:, 0].min(), c.xy[:, 1].min()

    # transiente: o numero de instantes e o do arquivo mais curto (uma rodada interrompida pode cortar a ultima linha)
    pw, pa, sw = (le(n, 0) for n in ("snaps_pw.csv", "snaps_pa.csv", "snaps_sw.csv"))
    n = min(len(pw), len(pa), len(sw))
    c.t, c.pw, c.pa, c.sw = pw[:n, 0], pw[:n, 1:], pa[:n, 1:], sw[:n, 1:]

    fs = le("fs.csv")
    c.fs_t, c.fs = np.round(fs[fs[:, 1] > 0, 0], 6), fs[fs[:, 1] > 0, 1]  # FS <= 0: nao calculado
    if len(c.fs) == 0:
        raise SystemExit("fs.csv nao tem nenhum FS calculado")
    resumo = dict(linha.split() for linha in open(arq("summary.txt")) if linha.strip())
    c.intensidade, c.t_chuva = float(resumo["I"]), float(resumo["train"])

    # pontos de integracao: funcoes de forma P1 nos pontos, para projetar nos vertices, e a celula de cada ponto
    pts = le("pontos.csv")
    c.pt_tri, c.pt_xy, c.pt_peso = pts[:, 0].astype(int), pts[:, 1:3], pts[:, 3]
    c.pt_N = coordenadas_baricentricas(c.xy[c.tri[c.pt_tri]], c.pt_xy)
    c.pt_celulas = celulas_dos_pontos(c.xy[c.tri[c.pt_tri]], c.pt_N)
    if np.any(np.bincount(c.pt_tri * 3 + c.pt_N.argmax(axis=1)) > 1):
        raise SystemExit("cada triangulo deve ter um ponto de integracao perto de cada vertice (regra de 3 pontos)")
    c.pt_soma = np.bincount(c.tri[c.pt_tri].ravel(), (c.pt_N * c.pt_peso[:, None]).ravel(), minlength=len(c.xy))

    # estados do SRM, por instante
    nos = le("srm_nos.csv")
    pontos = le("srm_pontos.csv")
    nos[:, 0], pontos[:, 0] = np.round(nos[:, 0], 6), np.round(pontos[:, 0], 6)
    c.srm_t = np.array(sorted(set(nos[:, 0]) & set(pontos[:, 0]) & set(c.fs_t)))
    if len(c.srm_t) == 0:
        raise SystemExit("srm_nos.csv e srm_pontos.csv nao tem nenhum instante com FS calculado")
    c.srm_u, c.srm_j2, c.srm_alfa, c.srm_tipo = [], [], [], []
    for t in c.srm_t:
        a = nos[nos[:, 0] == t]
        c.srm_u.append(a[np.argsort(np.searchsorted(c.nos, a[:, 1].astype(int)))][:, 2:6])
        b = pontos[pontos[:, 0] == t]
        b = b[np.argsort(b[:, 1])]
        c.srm_j2.append(b[:, 2])
        c.srm_alfa.append(b[:, 3])
        c.srm_tipo.append(b[:, 4].astype(int))
    return c


def coordenadas_baricentricas(triangulos, p):
    """Coordenadas baricentricas de p (n, 2) nos triangulos (n, 3, 2): sao as funcoes de forma P1."""
    a, b, c = triangulos[:, 0], triangulos[:, 1], triangulos[:, 2]
    v0, v1, v2 = b - a, c - a, p - a
    det = v0[:, 0] * v1[:, 1] - v1[:, 0] * v0[:, 1]
    l1 = (v2[:, 0] * v1[:, 1] - v1[:, 0] * v2[:, 1]) / det
    l2 = (v0[:, 0] * v2[:, 1] - v2[:, 0] * v0[:, 1]) / det
    return np.column_stack([1.0 - l1 - l2, l1, l2])


def celulas_dos_pontos(triangulos, N):
    """Divide cada triangulo em tres quadrilateros, um por ponto de integracao: o do vertice mais proximo do ponto."""
    k = N.argmax(axis=1)
    i = np.arange(len(N))
    v = triangulos[i, k]
    v1 = triangulos[i, (k + 1) % 3]
    v2 = triangulos[i, (k + 2) % 3]
    centro = triangulos.mean(axis=1)
    return np.stack([v, 0.5 * (v + v1), centro, 0.5 * (v + v2)], axis=1)


def projeta_nos_vertices(c, valor):
    """Projecao concentrada dos pontos de integracao nos vertices: s_i = soma(N_i w s) / soma(N_i w). Campo continuo."""
    soma = np.bincount(c.tri[c.pt_tri].ravel(), (c.pt_N * (c.pt_peso * valor)[:, None]).ravel(), minlength=len(c.xy))
    return soma / c.pt_soma


# ---------------------------------- escalas ----------------------------------
def escala_discreta(mapa, niveis, de=0.12, ate=0.92, abaixo=None):
    """Escala por faixas de um mapa de um so tom (claro para escuro); abaixo: cor dos valores sob o primeiro nivel."""
    cores = list(plt.get_cmap(mapa)(np.linspace(de, ate, len(niveis) - 1)))
    if abaixo:
        return from_levels_and_colors(niveis, [abaixo] + cores, extend="min")
    return from_levels_and_colors(niveis, cores)


def escala_pressao(niveis, p_min, p_max):
    """Escala divergente por faixas, com o cinza em p = 0 e cada braco usando toda a sua metade do mapa."""
    centro = TwoSlopeNorm(vcenter=0.0, vmin=min(p_min, -1.0), vmax=max(p_max, 1.0))
    meios = 0.5 * (niveis[:-1] + niveis[1:])
    return from_levels_and_colors(niveis, list(DIVERGENTE(centro(meios))))


def niveis_por_passo(minimo, maximo, passo):
    return np.arange(math.floor(minimo / passo) * passo, math.ceil(maximo / passo) * passo + 0.5 * passo, passo)


# ---------------------------------- paineis ----------------------------------
def formata_mapa(c, ax, esquerda, embaixo):
    ax.set_xlim(c.x0, c.L)
    ax.set_ylim(c.y0, c.y0 + 1.13 * (c.H - c.y0))
    ax.set_aspect("equal", adjustable="box")
    ax.set_facecolor(FUNDO)
    for lado in ("top", "right", "left", "bottom"):
        ax.spines[lado].set_visible(False)
    ax.tick_params(length=2.5, labelleft=esquerda, labelbottom=embaixo, left=esquerda, bottom=embaixo)
    if esquerda:
        ax.set_ylabel("y (m)")
    if embaixo:
        ax.set_xlabel("x (m)")
    ax.xaxis.set_major_formatter(FORMATO)
    ax.yaxis.set_major_formatter(FORMATO)


def contorno_talude(c, ax):
    ax.plot(c.contorno[:, 0], c.contorno[:, 1], color=TINTA, lw=0.8, zorder=6, solid_joinstyle="miter")


def cota_superficie(c, x):
    """A maior cota em que a vertical por x cruza o contorno."""
    cotas = []
    for (x0, y0), (x1, y1) in zip(c.contorno[:-1], c.contorno[1:]):
        if x0 != x1 and min(x0, x1) <= x <= max(x0, x1):
            cotas.append(y0 + (y1 - y0) * (x - x0) / (x1 - x0))
    return max(cotas)


def nos_da_superficie(c):
    """Vertices do contorno que nao estao so na base nem nas laterais."""
    i = np.array(c.borda[:-1])
    a, b = c.xy[i], c.xy[np.roll(i, -1)]
    lateral = (a[:, 0] == b[:, 0]) & ((a[:, 0] == c.x0) | (a[:, 0] == c.L))
    base = (a[:, 1] == c.y0) & (b[:, 1] == c.y0)
    sup = ~(lateral | base)
    return np.unique(np.concatenate([i[sup], np.roll(i, -1)[sup]]))


def rotula_isolinhas(c, ax, campo, niveis, colunas):
    """Valor das isolinhas em colunas fixas: os rotulos deslizam com a linha em vez de saltar de um quadro para outro."""
    interp = mtri.LinearTriInterpolator(c.malha, campo, trifinder=c.localiza)
    for x in colunas:
        ys = np.linspace(c.y0, cota_superficie(c, x), 300)
        zs = interp(np.full_like(ys, x), ys)
        ocupado = []
        for nivel in niveis:
            if nivel == 0:
                continue
            for k in np.nonzero(np.ma.filled((zs[:-1] - nivel) * (zs[1:] - nivel), 1.0) < 0)[0]:
                y = ys[k] + (ys[k + 1] - ys[k]) * (nivel - zs[k]) / (zs[k + 1] - zs[k])
                if all(abs(y - o) > 1.5 for o in ocupado) and c.y0 + 0.7 < y < ys[-1] - 0.7:  # m
                    ocupado.append(y)
                    ax.text(x, y, virgula(nivel), fontsize=7, color=TINTA, ha="center", va="center", zorder=5,
                            bbox=dict(boxstyle="round,pad=0.12", fc=FUNDO, ec="none", alpha=0.8))


def painel_pw(c):
    p_min, p_max = c.pw.min(), c.pw.max()
    faixas = niveis_por_passo(p_min, p_max, PASSO_COR)
    linhas = niveis_por_passo(p_min, p_max, PASSO_ISOLINHA)
    linhas = linhas[linhas != 0.0]
    mapa, norma = escala_pressao(faixas, p_min, p_max)
    sup = nos_da_superficie(c)
    x_setas = np.linspace(c.x0 + 0.02 * (c.L - c.x0), c.x0 + 0.98 * (c.L - c.x0), 28)
    y_setas = np.array([cota_superficie(c, x) for x in x_setas])
    altura = 0.08 * (c.H - c.y0)

    def desenha(ax, q):
        pw = c.pw[q.i]
        ax.tricontourf(c.malha, pw, levels=faixas, cmap=mapa, norm=norma)
        visiveis = linhas[(linhas > pw.min()) & (linhas < pw.max())]
        if len(visiveis):
            ax.tricontour(c.malha, pw, levels=visiveis, colors=TINTA, linewidths=0.45, alpha=0.8)
        if pw.min() < 0.0 < pw.max():
            ax.tricontour(c.malha, pw, levels=[0.0], colors=MAGENTA, linewidths=2.0)
        rotula_isolinhas(c, ax, pw, visiveis, (c.x0 + 0.26 * (c.L - c.x0), c.x0 + 0.84 * (c.L - c.x0)))
        # a face de percolacao tambem e p_w = 0, mas fica sobre o contorno, onde a isolinha nao aparece
        face = sup[np.abs(pw[sup]) < 1.0e-6]
        ax.plot(c.xy[face, 0], c.xy[face, 1], ls="none", marker="o", ms=3.0, mfc=MAGENTA, mec="none", zorder=7,
                clip_on=False)
        contorno_talude(c, ax)
        if q.chove:
            ax.quiver(x_setas, y_setas + altura + 0.4, 0 * x_setas, np.full_like(x_setas, -altura), angles="xy",
                      scale_units="xy", scale=1.0, color=CHUVA, width=0.0022, headwidth=4.5, headlength=5, zorder=7)

    titulo = (r"Pressão de água $p_w$ (kPa): isolinhas a cada 25 kPa; magenta: $p_w$ = 0")
    return SimpleNamespace(chave="pw", titulo=titulo, desenha=desenha, mapa=mapa, norma=norma, srm=False,
                           marcas=niveis_por_passo(p_min, p_max, PASSO_ISOLINHA), rotulo="kPa", mapa_barra=True)


def painel_succao(c):
    s = c.pa - c.pw
    s_max = max(s.max(), 1.0)
    faixas = niveis_por_passo(0.0, s_max, PASSO_COR)
    linhas = niveis_por_passo(0.0, s_max, PASSO_ISOLINHA)[1:]
    mapa, norma = escala_discreta("Oranges", faixas, abaixo=NEUTRO)

    def desenha(ax, q):
        campo = s[q.i]
        # s no limite superior (o patamar de -75 kPa em t = 0): um pouco abaixo, para a faixa de cima ficar fechada
        ax.tricontourf(c.malha, np.minimum(campo, faixas[-1] - 1.0e-9), levels=faixas, cmap=mapa, norm=norma, extend="min")
        visiveis = linhas[(linhas > campo.min()) & (linhas < campo.max())]
        if len(visiveis):
            ax.tricontour(c.malha, campo, levels=visiveis, colors=TINTA, linewidths=0.45, alpha=0.8)
            rotula_isolinhas(c, ax, campo, visiveis, (c.x0 + 0.26 * (c.L - c.x0), c.x0 + 0.84 * (c.L - c.x0)))
        contorno_talude(c, ax)

    titulo = r"Sucção $s = p_a - p_w$ (kPa); cinza: saturado ($s \leq 0$)"
    return SimpleNamespace(chave="succao", titulo=titulo, desenha=desenha, mapa=mapa, norma=norma, srm=False,
                           marcas=niveis_por_passo(0.0, s_max, PASSO_ISOLINHA), rotulo="kPa", mapa_barra=True)


def painel_saturacao(c):
    faixas = niveis_por_passo(c.sw.min(), 1.0, 0.05)
    mapa, norma = escala_discreta("Blues", faixas)

    def desenha(ax, q):
        # S_w = 1 no limite superior: um pouco abaixo, para a faixa de cima ficar fechada
        ax.tricontourf(c.malha, np.minimum(c.sw[q.i], 1.0 - 1.0e-9), levels=faixas, cmap=mapa, norm=norma)
        contorno_talude(c, ax)

    titulo = r"Saturação $S_w$"
    return SimpleNamespace(chave="saturacao", titulo=titulo, desenha=desenha, mapa=mapa, norma=norma, srm=False,
                           marcas=faixas[::2], rotulo="", mapa_barra=True)


def painel_ar(c):
    p_abs = np.abs(c.pa).max()
    imposto = p_abs < 1.0e-6
    if imposto:
        mapa, norma = None, None
        marcas = None
    else:
        passo = max(5.0, math.ceil(p_abs / 10.0 / 5.0) * 5.0)
        faixas = niveis_por_passo(-p_abs, p_abs, passo)
        mapa, norma = escala_pressao(faixas, -p_abs, p_abs)
        marcas = faixas[::2]

    def desenha(ax, q):
        if imposto:
            ax.fill(c.contorno[:, 0], c.contorno[:, 1], color=NEUTRO, zorder=1)
            ax.text(0.5 * (c.x0 + c.L), c.y0 + 0.4 * (c.H - c.y0), r"$p_a$ = 0 em todo o domínio" "\n(imposto nesta variante)",
                    ha="center", va="center", color=TINTA2, fontsize=10)
        else:
            ax.tricontourf(c.malha, c.pa[q.i], levels=faixas, cmap=mapa, norm=norma)
        contorno_talude(c, ax)

    titulo = r"Pressão de ar $p_a$ (kPa)"
    return SimpleNamespace(chave="ar", titulo=titulo, desenha=desenha, mapa=mapa, norma=norma, srm=False,
                           marcas=marcas, rotulo="kPa", mapa_barra=not imposto, estatico=imposto)


def titulo_srm(c, q):
    return f"t = {rotulo_tempo(c.srm_t[q.j])}, F = {virgula(c.fs[np.argmin(np.abs(c.fs_t - c.srm_t[q.j]))], 2)}"


def painel_u(c, total):
    faixas = np.linspace(0.0, 1.0, 11)
    mapa, norma = escala_discreta("Purples", faixas, de=0.05, ate=0.65)
    passo = 1.8  # m entre as setas
    gx, gy = np.meshgrid(np.arange(c.x0 + 0.5 * passo, c.L, passo), np.arange(c.y0 + 0.5 * passo, c.H, passo))
    dentro = c.localiza(gx, gy) >= 0
    gx, gy = gx[dentro], gy[dentro]

    def desenha(ax, q):
        u = c.srm_u[q.j][:, 0:2] if total else c.srm_u[q.j][:, 0:2] - c.srm_u[q.j][:, 2:4]
        modulo = np.hypot(u[:, 0], u[:, 1])
        maximo = max(modulo.max(), 1.0e-12)
        ax.tricontourf(c.malha, modulo / maximo, levels=faixas, cmap=mapa, norm=norma)
        ux = mtri.LinearTriInterpolator(c.malha, u[:, 0], trifinder=c.localiza)(gx, gy).filled(0.0)
        uy = mtri.LinearTriInterpolator(c.malha, u[:, 1], trifinder=c.localiza)(gx, gy).filled(0.0)
        visivel = np.hypot(ux, uy) > 0.02 * maximo
        ax.quiver(gx[visivel], gy[visivel], ux[visivel], uy[visivel], angles="xy", scale_units="xy",
                  scale=maximo / (2.2 * passo), color=TINTA, width=0.0022, headwidth=3.5, headlength=4,
                  headaxislength=3.5, pivot="tail", zorder=5)
        contorno_talude(c, ax)
        nome = "Deslocamento total" if total else "Deslocamento desde o geostático"
        ax.set_title(f"{nome} ({titulo_srm(c, q)}); máx. = {virgula(maximo, 2)} m", loc="left")

    return SimpleNamespace(chave="u", titulo=None, desenha=desenha, mapa=mapa, norma=norma, srm=True,
                           marcas=faixas[::2], rotulo=r"$|\mathbf{u}|\ /\ |\mathbf{u}|_{máx}$", mapa_barra=True)


def painel_log(c, chave, valores, nome, mapa_base):
    campos = [projeta_nos_vertices(c, v) for v in valores]
    maior = max(max(v.max() for v in valores), 10.0 * LOG_MIN)
    faixas = np.arange(math.log10(LOG_MIN), math.ceil(2.0 * math.log10(maior)) / 2.0 + 0.25, 0.5)
    mapa, norma = escala_discreta(mapa_base, faixas, de=0.15, ate=0.95, abaixo=NEUTRO)

    def desenha(ax, q):
        campo = np.log10(np.maximum(campos[q.j], 1.0e-30))
        ax.tricontourf(c.malha, np.minimum(campo, faixas[-1] - 1.0e-9), levels=faixas, cmap=mapa, norm=norma, extend="min")
        contorno_talude(c, ax)
        ax.set_title(f"{nome} no SRM ({titulo_srm(c, q)}); cinza: < " + rf"$10^{{{round(math.log10(LOG_MIN))}}}$", loc="left")

    marcas = faixas[np.isclose(faixas, np.round(faixas))]
    return SimpleNamespace(chave=chave, titulo=None, desenha=desenha, mapa=mapa, norma=norma, srm=True, marcas=marcas,
                           rotulo="", mapa_barra=True,
                           formato=FuncFormatter(lambda v, _: rf"$10^{{{int(round(v))}}}$"))


def painel_falha(c):
    presentes = sorted({int(t) for estado in c.srm_tipo for t in np.unique(estado)}, key=list(TIPOS).index)

    def desenha(ax, q):
        cores = [TIPOS[int(t)][1] for t in c.srm_tipo[q.j]]
        ax.add_collection(PolyCollection(c.pt_celulas, facecolors=cores, edgecolors=cores, linewidths=0.15, zorder=2))
        contorno_talude(c, ax)
        ax.set_title(f"Tipo de falha no SRM ({titulo_srm(c, q)})", loc="left")
        area = np.array([c.pt_peso[c.srm_tipo[q.j] == t].sum() for t in presentes]) / c.pt_peso.sum()
        marcas = [Patch(fc=TIPOS[t][1], ec=EIXO if t == 0 else "none", lw=0.5,
                        label=f"{TIPOS[t][0]}: {virgula(100.0 * a, 0)}%") for t, a in zip(presentes, area)]
        ax.legend(handles=marcas, loc="upper right", fontsize=8, frameon=False, handlelength=1.0, handleheight=1.0,
                  borderaxespad=0.1, labelspacing=0.25)

    return SimpleNamespace(chave="falha", titulo=None, desenha=desenha, mapa=None, norma=None, srm=True, marcas=None,
                           rotulo="", mapa_barra=False)


def eixo_fs(c, ax, t_fim):
    ax.axvspan(0.0, min(c.t_chuva, t_fim), color=AZUL_LINHA, alpha=0.08, lw=0, zorder=0)
    ax.axhline(1.0, color=TINTA2, lw=0.9, zorder=1)
    ax.plot(c.fs_t, c.fs, "-", color=AZUL_LINHA, lw=2.0, zorder=3, solid_capstyle="round", solid_joinstyle="round")
    ax.plot(c.fs_t, c.fs, ls="none", marker="o", ms=5.5, mfc=AZUL_LINHA, mec=FUNDO, mew=1.0, zorder=4)
    ax.set_xlim(0.0, t_fim)
    ax.set_ylim(min(0.9, c.fs.min() - 0.05), 1.05 * c.fs.max())
    passo = 6.0 if t_fim > 12.0 else 1.0
    ax.set_xticks(np.arange(0.0, t_fim + 0.5 * passo, passo))
    ax.grid(axis="y", color=GRADE, lw=0.8)
    ax.set_axisbelow(True)
    for lado in ("top", "right"):
        ax.spines[lado].set_visible(False)
    ax.tick_params(length=2.5)
    ax.set_xlabel("t (h)")
    ax.set_ylabel("FS")
    ax.xaxis.set_major_formatter(FORMATO)
    ax.yaxis.set_major_formatter(FORMATO)
    topo = ax.get_ylim()[1]
    ax.text(0.5 * min(c.t_chuva, t_fim), topo, f"chuva de {virgula(c.intensidade)} mm/h", color=TINTA2, fontsize=8,
            ha="center", va="top")
    ax.text(t_fim, 1.0, "FS = 1 ", color=TINTA2, fontsize=8, ha="right", va="bottom")


def painel_fs(c):
    t_fim = max(c.t.max(), c.fs_t.max(), 1.0e-3)

    def desenha(ax, q):
        eixo_fs(c, ax, t_fim)
        ax.axvline(q.t, color=TINTA, lw=1.0, zorder=2)
        ax.plot([q.t], [q.fs], ls="none", marker="o", ms=9, mfc=TINTA, mec=FUNDO, mew=2.0, zorder=6, clip_on=False)

    k = int(np.argmin(c.fs))
    titulo = f"FS por redução de resistência; mínimo: {virgula(c.fs[k], 2)} em t = {virgula(c.fs_t[k])} h"
    return SimpleNamespace(chave="fs", titulo=titulo, desenha=desenha, mapa=None, norma=None, srm=False, marcas=None,
                           rotulo="", mapa_barra=False, eixo_livre=True)


# ---------------------------------- montagem ----------------------------------
def monta_paineis(c, campos, u_total):
    construtores = {
        "pw": lambda: painel_pw(c),
        "succao": lambda: painel_succao(c),
        "saturacao": lambda: painel_saturacao(c),
        "ar": lambda: painel_ar(c),
        "u": lambda: painel_u(c, u_total),
        "fs": lambda: painel_fs(c),
        "j2": lambda: painel_log(c, "j2", c.srm_j2, r"$\sqrt{J_2(\varepsilon^p)}$", "Reds"),
        "hardening": lambda: painel_log(c, "hardening", c.srm_alfa, r"Hardening $\alpha$", "Reds"),
        "falha": lambda: painel_falha(c),
    }
    desconhecidos = [n for n in campos if n not in construtores]
    if desconhecidos:
        raise SystemExit(f"campos desconhecidos: {', '.join(desconhecidos)}; opcoes: {', '.join(construtores)}")
    return [construtores[n]() for n in campos]


def monta_figura(c, paineis):
    n = len(paineis)
    colunas = min(3, n)
    linhas = math.ceil(n / colunas)
    celula_l, celula_a = 6.4, 3.25
    topo, base = 0.72, 0.33
    largura, altura = colunas * celula_l, linhas * celula_a + topo + base
    # o video pede dimensoes pares em pixels
    largura = 2 * math.ceil(largura * 100 / 2) / 100
    altura = 2 * math.ceil(altura * 100 / 2) / 100
    fig = plt.figure(figsize=(largura, altura), dpi=100, facecolor=FUNDO)
    razao = 1.13 * (c.H - c.y0) / (c.L - c.x0)
    mapa_l = 4.86
    mapa_a = mapa_l * razao
    for k, p in enumerate(paineis):
        col, lin = k % colunas, k // colunas
        x = (col * celula_l + 0.56) / largura
        y = (altura - topo - (lin + 1) * celula_a + 0.46) / altura
        p.ax = fig.add_axes([x, y, mapa_l / largura, mapa_a / altura])
        p.esquerda = col == 0 or getattr(p, "eixo_livre", False)
        p.embaixo = lin == linhas - 1 or (lin + 1) * colunas + col >= n or getattr(p, "eixo_livre", False)
        if p.mapa_barra:
            cax = fig.add_axes([x + (mapa_l + 0.09) / largura, y, 0.13 / largura, mapa_a / altura])
            barra = fig.colorbar(ScalarMappable(p.norma, p.mapa), cax=cax, ticks=p.marcas,
                                 format=getattr(p, "formato", FORMATO), spacing="proportional")
            barra.outline.set_visible(False)
            barra.ax.tick_params(length=2, labelsize=8, pad=2)
            barra.minorticks_off()
            if p.rotulo and p.rotulo != "kPa":
                barra.set_label(p.rotulo, fontsize=8, labelpad=2)
    return fig, largura, altura


def fs_no_tempo(c, t):
    return float(np.interp(t, c.fs_t, c.fs))


def escolhe_quadros(c, passo, alvo=100):
    n = len(c.t)
    if passo is None:
        passo = max(1, round(n / alvo))
    idx = set(range(0, n, passo)) | {n - 1}
    idx |= {int(np.argmin(np.abs(c.t - t))) for t in c.fs_t}  # todo FS calculado aparece num quadro
    return sorted(idx)


def principal():
    ap = argparse.ArgumentParser(description="Animacao de um caso de talude sob chuva")
    ap.add_argument("pasta")
    ap.add_argument("--passo", type=int, default=None)
    ap.add_argument("--fps", type=float, default=4.0)
    ap.add_argument("--campos", default="pw,succao,saturacao,ar,fs,u,j2,hardening,falha")
    ap.add_argument("--largura-gif", type=int, default=1400)
    ap.add_argument("--u-total", action="store_true")
    ap.add_argument("--manter-quadros", action="store_true")
    args = ap.parse_args()

    plt.rcParams.update({
        "font.size": 10, "text.color": TINTA, "axes.labelcolor": TINTA2, "axes.edgecolor": EIXO, "axes.linewidth": 0.8,
        "axes.titlesize": 10.5, "axes.titlecolor": TINTA, "axes.titlelocation": "left", "axes.titlepad": 5,
        "xtick.color": APAGADO, "ytick.color": APAGADO, "xtick.labelcolor": TINTA2, "ytick.labelcolor": TINTA2,
        "xtick.labelsize": 8.5, "ytick.labelsize": 8.5, "axes.labelsize": 9, "mathtext.default": "regular",
    })

    arq = lambda nome: os.path.join(args.pasta, nome)
    c = le_caso(args.pasta)
    paineis = monta_paineis(c, [n.strip() for n in args.campos.split(",") if n.strip()], args.u_total)

    # FS(t) sozinho
    fig_fs, ax_fs = plt.subplots(figsize=(7.0, 4.0), layout="constrained", facecolor=FUNDO)
    eixo_fs(c, ax_fs, max(c.t.max(), c.fs_t.max()))
    k = int(np.argmin(c.fs))
    ax_fs.set_title(f"FS por redução de resistência; mínimo: {virgula(c.fs[k], 2)} em t = {virgula(c.fs_t[k])} h")
    fig_fs.savefig(arq("fs_tempo.png"), dpi=150)
    plt.close(fig_fs)

    fig, largura, altura = monta_figura(c, paineis)
    texto_tempo = fig.text(0.56 / largura, 1.0 - 0.34 / altura, "", fontsize=16, fontweight="bold", color=TINTA, va="center")
    if any(p.srm for p in paineis):
        fig.text(0.56 / largura, 0.12 / altura,
                 "Deslocamento, deformação plástica, hardening e tipo de falha: último F convergido do SRM; mudam só nos "
                 "instantes em que o FS é calculado (transiente com esqueleto rígido).",
                 fontsize=8.5, color=TINTA2, va="center", wrap=True)

    quadros = escolhe_quadros(c, args.passo)
    pasta_png = arq("quadros") if args.manter_quadros else tempfile.mkdtemp(prefix="animacao_")
    os.makedirs(pasta_png, exist_ok=True)
    feitos = {}
    for numero, i in enumerate(quadros):
        t = c.t[i]
        j = max(np.searchsorted(c.srm_t, t + 1.0e-6) - 1, 0)  # o ultimo estado do SRM ate este instante
        q = SimpleNamespace(t=t, i=i, j=int(j), fs=fs_no_tempo(c, t), chove=0.0 < t <= c.t_chuva + 1.0e-9)
        exato = np.any(np.abs(c.fs_t - t) < 1.0e-6)
        estado = f"chuva de {virgula(c.intensidade)} mm/h" if q.chove else "sem chuva"
        texto_tempo.set_text(f"t = {rotulo_tempo(t)}   ·   {estado}   ·   FS {'=' if exato else '≈'} {virgula(q.fs, 2)}")
        for p in paineis:
            # paineis do SRM e paineis estaticos so mudam quando ha um estado novo
            chave = q.j if p.srm else (0 if getattr(p, "estatico", False) else i)
            if feitos.get(p.chave) == chave:
                continue
            feitos[p.chave] = chave
            p.ax.clear()
            if not getattr(p, "eixo_livre", False):
                formata_mapa(c, p.ax, p.esquerda, p.embaixo)
            p.desenha(p.ax, q)
            if p.titulo:
                p.ax.set_title(p.titulo, loc="left")
        fig.savefig(os.path.join(pasta_png, f"q{numero:05d}.png"), dpi=100, facecolor=FUNDO,
                    pil_kwargs={"compress_level": 1})
        if (numero + 1) % 20 == 0:
            print(f"  {numero + 1}/{len(quadros)} quadros", flush=True)
    plt.close(fig)

    # video e GIF a partir dos PNG; o GIF usa uma so paleta para todos os quadros, sem pontilhado
    entrada = ["-framerate", f"{args.fps:g}", "-i", os.path.join(pasta_png, "q%05d.png")]
    pausa = "tpad=stop_mode=clone:stop_duration=2"  # segura o ultimo quadro antes de recomecar
    ffmpeg = ["ffmpeg", "-y", "-loglevel", "error"]
    subprocess.run(ffmpeg + entrada + ["-vf", pausa, "-c:v", "libx264", "-crf", "18", "-pix_fmt", "yuv420p",
                                       arq("chuva.mp4")], check=True)
    gravados = ["fs_tempo.png", "chuva.mp4"]
    if args.largura_gif > 0:
        filtro = (f"{pausa},scale={args.largura_gif}:-1:flags=lanczos,split[a][b];"
                  "[a]palettegen=max_colors=256:stats_mode=full[p];[b][p]paletteuse=dither=none:diff_mode=rectangle")
        subprocess.run(ffmpeg + entrada + ["-filter_complex", filtro, "-loop", "0", arq("chuva.gif")], check=True)
        gravados.append("chuva.gif")
    if not args.manter_quadros:
        shutil.rmtree(pasta_png)
    print(f"{len(quadros)} quadros: {', '.join(gravados)} em {args.pasta}")


if __name__ == "__main__":
    principal()
