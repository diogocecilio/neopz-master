#!/usr/bin/env python3
"""Figuras do GIM com refinamento adaptativo (saidas da SlopeAnalysis da Common e do SlopeRandomFieldRefinement).

Uso: python3 figuras_gim.py [pasta]   (padrao: pasta atual, onde o programa rodou)
Para cada rodada <r> = <modo>_<caso> (modo seco ou rebaix; caso ceron ou talude10) le gim_<r>_<mat>.csv,
_nivel<k>.vtk e _nivel<k>_pontos.csv (mat = up, unsat) e grava:
  gim_<r>_<mat>.png          um painel por nivel: malha, pontos da banda e Gamma
  gim_<r>_convergencia.png   Gamma x equacoes livres, com o alvo (coluna alvo do csv)
  gim_<r>_percolacao.png     so no rebaixamento: excesso de poropressao u e a banda do ultimo nivel
"""
import csv
import glob
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.collections import LineCollection, PolyCollection
from matplotlib.patches import Patch
import numpy as np

LIMIAR = 0.1  # banda: sqrt(J2(eps_p)) >= LIMIAR x maximo, como na marcacao do refinamento
NOMES = {"up": "material u-p (TPZMatPoroElastoPlasticUP)", "unsat": "material nao saturado (TPZMatUnsaturatedUPwPa)"}
MODOS = {"seco": "GIM seco", "rebaix": "GIM com rebaixamento ($h_w = H$)",
         "rebaixforca": "GIM com rebaixamento, forca $-$grad u do MEF",
         "semianalitico": "GIM com rebaixamento, campo semi-analitico"}
CORES_P = {3: "#c6dbef", 4: "#6baed6", 5: "#2171b5"}  # triangulos com u de ordem maior que 2


def br(v, n=3):
    """numero com virgula decimal"""
    return f"{v:.{n}f}".replace(".", ",")


def le_escalares(linhas, i0, n):
    """SCALARS a partir da linha i0 (CELL_DATA ou POINT_DATA), n valores cada, ate a proxima secao"""
    campos = {}
    k = i0 + 1
    while k < len(linhas) and not linhas[k].startswith(("POINT_DATA", "CELL_DATA")):
        if linhas[k].startswith("SCALARS"):
            campos[linhas[k].split()[1]] = np.array([float(linhas[k + 2 + j]) for j in range(n)])
            k += 2 + n
        else:
            k += 1
    return campos


def le_vtk(arq):
    """nos (n x 2), triangulos (m x 3) e os campos de POINT_DATA e CELL_DATA do VTK gravado por WriteLevelVTK"""
    linhas = open(arq).read().split("\n")
    i = next(k for k, l in enumerate(linhas) if l.startswith("POINTS"))
    nn = int(linhas[i].split()[1])
    nos = np.array([[float(v) for v in linhas[i + 1 + j].split()[:2]] for j in range(nn)])
    i = next(k for k, l in enumerate(linhas) if l.startswith("CELLS"))
    nt = int(linhas[i].split()[1])
    tri = np.array([[int(v) for v in linhas[i + 1 + j].split()[1:]] for j in range(nt)])
    campos, celulas = {}, {}
    for k, l in enumerate(linhas):
        if l.startswith("POINT_DATA"):
            campos = le_escalares(linhas, k, nn)
        elif l.startswith("CELL_DATA"):
            celulas = le_escalares(linhas, k, nt)
    return nos, tri, campos, celulas


def arestas(tri):
    """todas as arestas e as do contorno (as que so tem um triangulo)"""
    cont = {}
    for t in tri:
        for a, b in ((t[0], t[1]), (t[1], t[2]), (t[2], t[0])):
            k = (min(a, b), max(a, b))
            cont[k] = cont.get(k, 0) + 1
    return np.array(list(cont.keys())), np.array([k for k, n in cont.items() if n == 1])


def le_resumo(arq):
    with open(arq) as f:
        linhas = [{k: float(v) for k, v in r.items()} for r in csv.DictReader(f)]
    return linhas if linhas and "alvo" in linhas[0] else None  # csv antigo (sem alvo): ignora


def banda(base):
    """pontos de integracao com sqrt(J2) >= LIMIAR x maximo"""
    p = np.loadtxt(base + "_pontos.csv", delimiter=",", skiprows=1)
    return p[p[:, 2] >= LIMIAR * p[:, 2].max()]


def janela(prefixo, n):
    """ampliacao em volta das bandas de todos os niveis (no ceron da 4 a 21 m por 2 a 10,5 m)"""
    b = np.vstack([banda(f"{prefixo}_nivel{k}") for k in range(n)])
    x0, x1, y0, y1 = b[:, 0].min(), b[:, 0].max(), b[:, 1].min(), b[:, 1].max()
    w = max(x1 - x0, y1 - y0)
    return x0 - 0.9 * w, x1 + 1.0 * w, y0 - 0.55 * w, y1 + 0.1 * w


def malha(ax, nos, tri):
    todas, borda = arestas(tri)
    ax.add_collection(LineCollection(nos[todas], colors="0.55", linewidths=0.35))
    ax.add_collection(LineCollection(nos[borda], colors="k", linewidths=1.2))


def figura_material(pasta, rodada, mat, resumo):
    n = len(resumo)
    prefixo = os.path.join(pasta, f"gim_{rodada}_{mat}")
    jan = janela(prefixo, n)
    fig, eixos = plt.subplots(n, 1, figsize=(9, 3.6 * n), squeeze=False)
    for k, (ax, r) in enumerate(zip(eixos[:, 0], resumo)):
        base = f"{prefixo}_nivel{k}"
        nos, tri, _, cel = le_vtk(base + ".vtk")
        ordem = cel.get("ordem_u", np.full(len(tri), 2.)).astype(int)
        legenda = []
        for o, cor in CORES_P.items():
            sel = ordem == o
            if sel.any():
                ax.add_collection(PolyCollection(nos[tri[sel]], facecolors=cor, edgecolors="none", zorder=0))
                legenda.append(Patch(facecolor=cor, label=f"u P{o} ({sel.sum()} triângulos)"))
        malha(ax, nos, tri)
        b = banda(base)
        ax.scatter(b[:, 0], b[:, 1], s=7, c="tab:red", marker="o", linewidths=0, zorder=3,
                   label=r"pontos da banda ($\sqrt{J_2}$ de $\varepsilon_p \geq 0{,}1$ do máximo)")
        ax.set_xlim(jan[0], jan[1])
        ax.set_ylim(jan[2], jan[3])
        ax.set_aspect("equal")
        ax.set_ylabel("y (m)")
        ax.set_title(f"({chr(97 + k)}) nível {k}: {int(r['triangulos'])} triângulos, {int(r['livres'])} equações livres, "
                     f"$\\Gamma$ = {br(r['gamma'])}; {len(b)} pontos da banda", fontsize=10)
        h, _ = ax.get_legend_handles_labels()
        ax.legend(handles=h + legenda, loc="lower left", fontsize=8, framealpha=0.9)
    eixos[-1, 0].set_xlabel("x (m)")
    modo, caso = rodada.split("_", 1)
    fig.suptitle(f"{MODOS[modo]}, caso {caso}, {NOMES[mat]}\n$\\Gamma$ final = {br(resumo[-1]['gamma'])} "
                 f"(alvo {br(resumo[0]['alvo'])})", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.975))
    saida = prefixo + ".png"
    fig.savefig(saida, dpi=130)
    plt.close(fig)
    return saida


def figura_convergencia(pasta, rodada, resumos):
    fig, ax = plt.subplots(figsize=(7, 4.2))
    estilos = {"up": ("tab:blue", "o", "u-p"), "unsat": ("tab:orange", "s", "não saturado")}
    for mat, res in resumos.items():
        cor, marca, rot = estilos[mat]
        neq = [r["livres"] for r in res]
        g = [r["gamma"] for r in res]
        ax.plot(neq, g, marker=marca, color=cor, label=rot, mfc="none" if mat == "unsat" else cor)
        for x, y in zip(neq, g):
            ax.annotate(br(y), (x, y), textcoords="offset points", xytext=(5, 6), fontsize=8, color=cor)
    alvo = next(iter(resumos.values()))[0]["alvo"]
    ax.axhline(alvo, color="k", ls="--", lw=1, label=f"alvo {br(alvo)}")
    ax.set_xscale("log")
    ax.margins(x=0.12, y=0.1)
    ax.set_xlabel("equações livres")
    ax.set_ylabel(r"$\Gamma$ (GIM)")
    modo, caso = rodada.split("_", 1)
    ax.set_title(f"{MODOS[modo]} ({caso}): convergência com o refinamento adaptativo", fontsize=10)
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    saida = os.path.join(pasta, f"gim_{rodada}_convergencia.png")
    fig.savefig(saida, dpi=130)
    plt.close(fig)
    return saida


def figura_percolacao(pasta, rodada, mat, resumo):
    """excesso de poropressao u = p - gamma_w y_artigo (Fig. 4b do artigo) no ultimo nivel, com a banda"""
    k = len(resumo) - 1
    base = os.path.join(pasta, f"gim_{rodada}_{mat}_nivel{k}")
    nos, tri, campos, _ = le_vtk(base + ".vtk")
    if "u_excesso" not in campos:
        return None
    fig, ax = plt.subplots(figsize=(9, 4.2))
    tr = mtri.Triangulation(nos[:, 0], nos[:, 1], tri)
    c = ax.tricontourf(tr, campos["u_excesso"], levels=20, cmap="viridis")
    ax.tricontour(tr, campos["u_excesso"], levels=10, colors="w", linewidths=0.4)
    fig.colorbar(c, ax=ax, label="excesso de poropressão u (kPa)")
    _, borda = arestas(tri)
    ax.add_collection(LineCollection(nos[borda], colors="k", linewidths=1.2))
    b = banda(base)
    ax.scatter(b[:, 0], b[:, 1], s=3, c="tab:red", linewidths=0, label=f"banda do nível {k}")
    ax.set_aspect("equal")
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax.set_title(f"Percolação após o rebaixamento ({rodada.split('_', 1)[1]}), {mat}: "
                 f"$\\Gamma$ = {br(resumo[-1]['gamma'])}", fontsize=10)
    ax.legend(loc="lower left", fontsize=8)
    fig.tight_layout()
    saida = os.path.join(pasta, f"gim_{rodada}_percolacao.png")
    fig.savefig(saida, dpi=130)
    plt.close(fig)
    return saida


def main():
    pasta = sys.argv[1] if len(sys.argv) > 1 else "."
    rodadas = set()
    for mat in ("up", "unsat"):
        for a in glob.glob(os.path.join(pasta, f"gim_*_*_{mat}.csv")):
            r = os.path.basename(a)[len("gim_"):-len(f"_{mat}.csv")]
            if r.split("_", 1)[0] in MODOS:
                rodadas.add(r)
    for rodada in sorted(rodadas):
        resumos = {}
        for mat in ("up", "unsat"):
            arq = os.path.join(pasta, f"gim_{rodada}_{mat}.csv")
            res = le_resumo(arq) if os.path.isfile(arq) else None
            if res:
                resumos[mat] = res
                print("gravado", figura_material(pasta, rodada, mat, res))
        if resumos:
            print("gravado", figura_convergencia(pasta, rodada, resumos))
            mat, res = next(iter(resumos.items()))
            perc = figura_percolacao(pasta, rodada, mat, res)
            if perc:
                print("gravado", perc)
    if not rodadas:
        sys.exit("figuras_gim: nenhum gim_<modo>_<caso>_<mat>.csv em " + pasta)


if __name__ == "__main__":
    main()
