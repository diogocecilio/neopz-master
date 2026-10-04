#!/usr/bin/env python3
"""Figura do mecanismo de colapso do solo A com beta = 35 graus (Fig. 8 do artigo), sem rebaixamento e com hw/H = 0,5.

Uso: mecanismo_artigo.py <vtk hw = 0> <vtk hw/H = 0,5> <saída sem extensão>
Os VTK são o percolacao_mc_gamma.scal_vec.0.vtk de duas execuções do comando det com vtk=1:
  SlopeSeepageRandom det caso=percolacao h=1 adapt=3 fs=0 fatormax=50 gw=9.81 gam=18 c=6 phi=32 beta=35 hw=0 \
      Lc=20 Lt=20 Hb=10 vtk=1
e o mesmo com hw=2.5. Escreve <saída>.pdf e <saída>.png; analise_artigo.py copia det/mech/mecanismo_A35.* para a
sua saída, de onde relatorio_tex.py a inclui.
"""
import sys

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.tri as mtri  # noqa: E402


def le_vtk(fn):
    """pontos, células e campos de um VTK legado ASCII (POINTS, CELLS, SCALARS, VECTORS)"""
    L = open(fn).read().split("\n")
    i, d, P, C = 0, {}, None, []
    while i < len(L):
        t = L[i].split()
        if t and t[0] == "POINTS":
            n = int(t[1])
            P = np.array([list(map(float, L[i + 1 + k].split())) for k in range(n)])
            i += n + 1
            continue
        if t and t[0] == "CELLS":
            n = int(t[1])
            C = [list(map(int, L[i + 1 + k].split()))[1:] for k in range(n)]
            i += n + 1
            continue
        if t and t[0] == "SCALARS":
            nome = t[1]
            i += 2 if L[i + 1].startswith("LOOKUP") else 1
            d[nome] = np.array([float(L[i + k]) for k in range(len(P))])
            i += len(P)
            continue
        if t and t[0] == "VECTORS":
            nome = t[1]
            i += 1
            d[nome] = np.array([list(map(float, L[i + k].split())) for k in range(len(P))])
            i += len(P)
            continue
        i += 1
    return P, C, d


def main(vtk0, vtk25, saida):
    casos = [(vtk0, r"$h_w = 0$ (sem rebaixamento; não converge com a malha)"), (vtk25, r"$h_w/H = 0{,}5$")]
    fig, axs = plt.subplots(1, 2, figsize=(10, 3.2))
    for ax, (fn, tit) in zip(axs, casos):
        P, C, d = le_vtk(fn)
        tris = []
        for c in C:
            if len(c) >= 4:
                tris += [[c[0], c[1], c[2]], [c[0], c[2], c[3]]]
        T = mtri.Triangulation(P[:, 0], P[:, 1], tris)
        v = d["PlasticStrainNorm"]
        ax.tripcolor(T, v, shading="gouraud", cmap="magma_r", vmax=np.percentile(v, 99.5))
        # contorno do talude: crista de 20 m, H = 5 m, base 10 m abaixo do pé
        ax.plot([0, 20, 20 + 5 / np.tan(np.radians(35)), 47.1], [15, 15, 10, 10], "k-", lw=0.8)
        ax.set_xlim(14, 32)
        ax.set_ylim(6, 15.5)
        ax.set_aspect("equal")
        ax.set_title(tit, fontsize=9)
        ax.set_xlabel("x (m)")
        ax.set_ylabel("y (m)")
    fig.suptitle(r"Solo A ($c$ = 6 kPa, $\varphi$ = 32°), $\beta$ = 35°, nível 3: "
                 r"$\|\varepsilon^p\|$ no colapso", fontsize=10)
    fig.tight_layout()
    fig.savefig(saida + ".pdf")
    fig.savefig(saida + ".png", dpi=130)


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print(__doc__)
        sys.exit(1)
    main(*sys.argv[1:])
