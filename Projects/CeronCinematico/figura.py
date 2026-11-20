"""Figuras do mecanismo critico: campo de poropressao (excesso u) ou de phi,
mecanismo cinematico critico e banda plastica do GIM (mecanismo.csv).

Uso:  python3 -I figura.py            (gera figs/mecanismo_det.pdf e
                                       figs/mecanismo_real0.pdf)"""
import math
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)
import cinematico as k  # noqa: E402

plt.rcParams.update({"font.size": 9, "axes.titlesize": 9, "axes.labelsize": 9,
                     "legend.fontsize": 8.5, "xtick.labelsize": 8.5, "ytick.labelsize": 8.5})

COR_KIN = "#c0392b"     # mecanismo cinematico
COR_GIM = "#333333"     # banda do GIM


def banda_gim(pasta, origem=(35.0, 35.0), frac=0.1):
    """Pontos de integracao do GIM com raiz(J2(eps_p)) >= frac do maior."""
    arq = os.path.join(pasta, "mecanismo.csv")
    if not os.path.exists(arq):
        return None
    M = np.loadtxt(arq, delimiter=",", skiprows=1, ndmin=2)
    M = M[M[:, 3] >= frac * M[:, 3].max()]
    return M[:, 1] - origem[0], origem[1] - M[:, 2]


def desenhar(ax, talude, campo, res, gim=None, campo_plot="u", titulo=""):
    T = campo.T
    # coordenadas de desenho: x, elevacao z = -y (crista em z = 0)
    from matplotlib.tri import Triangulation
    Td = Triangulation(T.x, -T.y, campo.tri)
    if campo_plot == "u":
        v, rot = campo.un, "excesso de poropressao u = p - $\\gamma_w y$ (kPa)"
        cmap, niveis = "Blues_r", np.linspace(-50, 0, 11)
    else:
        v, rot = np.degrees(campo.phin), "$\\phi'$ (graus)"
        cmap, niveis = "Greens", np.arange(2 * np.floor(v.min() / 2), np.ceil(v.max()) + 2, 2.0)
    cf = ax.tricontourf(Td, v, levels=niveis, cmap=cmap, extend="both")
    if campo_plot == "u":
        cs = ax.tricontour(Td, v, levels=niveis[1:-1:2], colors="k", linewidths=0.4, alpha=0.5)
        ax.clabel(cs, fmt="%g", fontsize=7.5)
    cb = plt.colorbar(cf, ax=ax, orientation="horizontal", pad=0.13, fraction=0.06, aspect=40)
    cb.set_label(rot)
    # superficie do terreno
    xs = np.array([-10, 0, talude.xT, 15])
    ax.plot(xs, -talude.y_sup(xs), color="k", lw=1.0)
    if gim is not None:
        ax.plot(gim[0], -gim[1], ".", ms=2.5, color=COR_GIM, alpha=0.55,
                label="GIM: pontos com $\\sqrt{J_2(\\varepsilon_p)}\\geq 0{,}1$ max")
    m = res.mecanismo
    X, Y = m.pontos()
    ax.plot(X, -Y, color=COR_KIN, lw=2.0,
            label=f"cinematico critico: $\\Gamma$ = {res.Gamma:.3f}")
    ax.plot([m.Cx], [m.Cy], "+", color=COR_KIN, ms=8, mew=1.5)
    ax.annotate("C", (m.Cx, m.Cy), xytext=(4, 2), textcoords="offset points", color=COR_KIN)
    ax.set_xlim(-6, 11)
    ax.set_ylim(-9, 4.5)
    ax.set_aspect("equal")
    ax.set_xlabel("x (m), origem na aresta da crista")
    ax.set_ylabel("z = -y (m)")
    ax.set_title(titulo, loc="left")
    ax.legend(loc="upper left", frameon=True, framealpha=0.9)


def figura(pasta, saida, campo_plot, titulo, res=None):
    if res is None:
        res = k.avaliar_realizacao(pasta)
    fig, ax = plt.subplots(figsize=(7.0, 5.4))
    desenhar(ax, res.talude, res.campo, res, banda_gim(pasta), campo_plot, titulo)
    fig.tight_layout()
    fig.savefig(saida)
    plt.close(fig)
    return res


if __name__ == "__main__":
    os.makedirs(os.path.join(AQUI, "figs"), exist_ok=True)
    d = os.path.join(AQUI, "dados")
    r = figura(os.path.join(d, "det"), os.path.join(AQUI, "figs", "mecanismo_det.pdf"), "u",
               "Rebaixamento total (h_w = H), campos medios, p_w do transiente NeoPZ")
    print("det:", r.Gamma)
    r = figura(os.path.join(d, "real0"), os.path.join(AQUI, "figs", "mecanismo_real0.pdf"),
               "phi", "Realizacao 0 (campos lognormais de c, phi, k)")
    print("real0:", r.Gamma)
