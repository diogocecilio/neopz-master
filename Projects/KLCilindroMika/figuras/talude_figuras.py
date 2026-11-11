#!/usr/bin/env python3
"""talude_figuras.py: figuras da comparacao no talude entre a KL em forma fechada do retangulo [0,25] x [0,10] (restrita
ao talude, a do SlopeRandomField), a KL de Galerkin do NeoPZ no talude (KLTaludeCeron) e a KL exata do talude (Nystrom).

  kl_truncamento.pdf  epsilon_M = 1 - media de sigma_M^2 no talude contra M (log-log): forma fechada restrita ao talude
                      (klfechada.KL2D.erro_talude, integrais exatas, M = 1 a 20000), Nystrom (M <= 210, extrapolado
                      de n = 80, 160, 320) e Galerkin, 1 - sum_{k<=M} lambda_k / 187,5 com os autovalores de cada malha/regra
  kl_espectros.pdf    lambda_k contra k: retangulo (forma fechada), talude (Nystrom n = 160) e Galerkin; e a diferenca
                      relativa dos autovalores de Galerkin contra os de Nystrom, k <= 200
  kl_modos_talude.png os 9 primeiros modos da KL exata do talude (Nystrom, dados/kl_nystrom_modos.csv)
  dados/kl_truncamento.csv e dados/kl_espectros.csv: os valores das duas primeiras figuras
  kl_variancia.png    sigma_M^2 com M = 200 e 2000: forma fechada restrita ao talude contra Galerkin P2 ceron_h5_a45_h025
                      (regra de ordem 10), nos pontos da rede de 15 pontos de cada triangulo (pt_*.csv do KLTaludeCeron),
                      escala divergente centrada em 1

Dados (so leitura): $KL_SCRATCH/kl_rel_tal/galerkin/av_<malha>_p<P>_o<ordem>.txt (todos os autovalores de Galerkin, um por
linha), pt_ceron_h5_a45_h025_p2_o10.csv (sigma_M^2 nos pontos), $KL_SCRATCH/verif_kl/mf_talude_n<n>_k210.npy (210
autovalores de Nystrom com n = 40, 80, 160, 320, da verificacao independente; sem eles, kl_rel_tal/nystrom_n160_lam.npy)
e ../../SlopeRandomFieldRefina/relatorio/dados/kl_nystrom_modos.csv. Os autovalores ficam guardados em dados/kl_espectros.csv.

Uso: talude_figuras.py [truncamento] [espectros] [variancia] [modos] [dados]   (padrao: todos)
"""
import os
import sys

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.tri import Triangulation

import estilo as E
from klfechada import KL2D, AREA_TALUDE

BASE = os.path.join(E.SCRATCH, "kl_rel_tal")
GAL = os.path.join(BASE, "galerkin")
# (arquivo, rotulo, cor, estilo de linha): regra do elemento e regra de ordem 10 (25 pontos)
CASOS_ELEM = [("ceron_h5_a45_p1_o2", "P1, a45 (609 eq.)", E.PALETA[3], "-"),
              ("ceron_h5_a45_h025_p1_o2", "P1, h025 (1314 eq.)", E.PALETA[4], "-"),
              ("ceron_h5_a45_fina_p1_o2", "P1, fina (2301 eq.)", E.PALETA[1], "-"),
              ("ceron_h5_a45_p2_o4", "P2, a45 (2320 eq.)", E.PALETA[2], "-"),
              ("ceron_h5_a45_h025_p2_o4", "P2, h025 (5085 eq.)", E.PALETA[0], "-")]
CASOS_O10 = [("ceron_h5_a45_fina_p1_o10", "P1, fina (2301 eq.)", E.PALETA[1], "-"),
             ("ceron_h5_a45_p2_o10", "P2, a45 (2320 eq.)", E.PALETA[2], "-"),
             ("ceron_h5_a45_h025_p2_o10", "P2, h025 (5085 eq.)", E.PALETA[0], "-")]


def nystrom():
    """autovalores de Nystrom do talude: extrapolacao de Richardson em h^2 e h^4 com n = 80, 160, 320 (rodadas da
    verificacao independente, $KL_SCRATCH/verif_kl/mf_talude_n<n>_k210.npy); a diferenca para a extrapolacao com
    n = 40, 80, 160 e a incerteza da referencia. Sem esses arquivos, usa a grade n = 160 sem extrapolar."""
    arqs = [os.path.join(E.SCRATCH, "verif_kl", "mf_talude_n%d_k210.npy" % n) for n in (40, 80, 160, 320)]
    if not all(os.path.exists(a) for a in arqs):
        lam = np.sort(np.load(os.path.join(BASE, "nystrom_n160_lam.npy")))[::-1]
        return lam, None, "Nyström, n = 160"
    l40, l80, l160, l320 = [np.sort(np.load(a))[::-1] for a in arqs]

    def extrap(a, b, c):  # a: h, b: 2h, c: 4h
        r1, r2 = (4. * a - b) / 3., (4. * b - c) / 3.
        return (16. * r1 - r2) / 15.
    L = extrap(l320, l160, l80)
    Lb = extrap(l160, l80, l40)
    return L, np.abs(Lb / L - 1.), "Nyström extrapolado"


def autovalores(nome):
    lam = np.loadtxt(os.path.join(GAL, "av_%s.txt" % nome))
    return np.sort(lam)[::-1]


def truncamento(kl, lam_ny):
    M = np.arange(1, kl.M + 1)
    eps_ret = kl.erro_talude()
    eps_ny = 1. - np.cumsum(lam_ny) / AREA_TALUDE
    fig, axs = plt.subplots(1, 2, figsize=(E.LARGURA, 3.0), sharey=True)
    for ax, casos, titulo in ((axs[0], CASOS_ELEM, "(a) Galerkin com a regra do elemento (ordem 2P)"),
                              (axs[1], CASOS_O10, "(b) Galerkin com a regra de ordem 10 (25 pontos)")):
        for nome, rot, cor, ls in casos:
            lam = autovalores(nome)
            e = 1. - np.cumsum(lam) / AREA_TALUDE
            ax.plot(np.arange(1, len(e) + 1), e, ls, color=cor, lw=1.1, label=rot)
        ax.plot(M, eps_ret, color=E.TINTA, lw=1.5, label="forma fechada de R, restrita a Ω")
        ax.plot(np.arange(1, len(eps_ny) + 1), eps_ny, "--", color=E.NEUTRO, lw=1.5, dashes=(3, 1.5),
                label="exata do talude (Nyström)", zorder=5)
        ax.axvline(2000, color=E.TINTA2, lw=0.6, ls=":")
        ax.text(2000 * 1.1, 0.25, "M = 2000", fontsize=6.5, color=E.TINTA2, rotation=90, va="center")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(1, 2.0e4)
        ax.set_ylim(1e-3, 1.0)
        E.eixo_log(ax)
        ax.set_xlabel("número de termos M")
        ax.set_title(titulo, loc="left", fontsize=7.5)
        ax.legend(loc="lower left", fontsize=6.2, handlelength=1.8)
    axs[0].set_ylabel(r"$\varepsilon_M = 1 - $ média de $\sigma_M^2$ no talude")
    fig.tight_layout(w_pad=0.6)
    E.salva(fig, "kl_truncamento")
    # alguns valores para o texto
    for m in (50, 200, 1000, 2000, 20000):
        print("  forma fechada: eps_%d = %.4f" % (m, eps_ret[m - 1]))
    for m in (50, 200):
        print("  Nystrom: eps_%d = %.4f" % (m, eps_ny[m - 1]))


def espectros(kl, lam_ny, inc):
    fig, axs = plt.subplots(1, 2, figsize=(E.LARGURA, 2.9))
    ax = axs[0]
    k = np.arange(1, kl.M + 1)
    ax.plot(k, kl.lam, color=E.TINTA, lw=1.5, label="retângulo R (forma fechada)")
    for nome, rot, cor, ls in (CASOS_ELEM[2], CASOS_O10[2]):
        lam = autovalores(nome)
        ax.plot(np.arange(1, len(lam) + 1), lam, color=cor, lw=1.1,
                label="Galerkin %s, %s" % (rot.split(" (")[0], "ordem 2P" if nome.endswith("o2") else "ordem 10"))
    ax.plot(np.arange(1, len(lam_ny) + 1), lam_ny, "--", color=E.NEUTRO, lw=1.5, dashes=(3, 1.5),
            label="talude (Nyström)")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(1, 2e4)
    ax.set_ylim(1e-4, 100)
    E.eixo_log(ax)
    ax.set_xlabel("índice k")
    ax.set_ylabel(r"$\lambda_k$ (m$^2$)")
    ax.set_title("(a) autovalores", loc="left", fontsize=7.5)
    ax.legend(loc="lower left", fontsize=6.2)
    ax = axs[1]
    KMAX = 100
    kk = np.arange(1, KMAX + 1)
    for casos in (CASOS_ELEM, CASOS_O10):
        for nome, rot, cor, ls in casos:
            if nome in ("ceron_h5_a45_p1_o2", "ceron_h5_a45_h025_p1_o2"):
                continue
            lam = autovalores(nome)[:KMAX]
            d = np.abs(lam / lam_ny[:KMAX] - 1.)
            o10 = nome.endswith("o10")
            ax.plot(kk, d, color=cor, lw=1.0, ls="--" if o10 else "-",
                    label="%s, %s" % (rot.split(" (")[0], "ordem 10" if o10 else "ordem 2P"))
    if inc is not None:
        ax.fill_between(kk, 1e-7, np.maximum(inc[:KMAX], 1e-7), color=E.NEUTRO, alpha=0.25, lw=0,
                        label="incerteza da referência")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(1, KMAX)
    ax.set_ylim(1e-6, 30)
    E.eixo_log(ax)
    ax.set_xlabel("índice k")
    ax.set_ylabel(r"$|\lambda_k^h / \lambda_k - 1|$ (Nyström como referência)")
    ax.set_title("(b) Galerkin no talude contra a KL exata do talude", loc="left", fontsize=7.5)
    ax.legend(loc="upper left", fontsize=5.6, ncol=3, handlelength=1.6, columnspacing=0.8)
    fig.tight_layout(w_pad=0.8)
    E.salva(fig, "kl_espectros")


def variancia():
    arq = os.path.join(GAL, "pt_ceron_h5_a45_h025_p2_o10.csv")
    cols = ["x", "y", "s2_num_M200", "s2_num_M2000", "s2_ret_M200", "s2_ret_M2000"]
    cab = open(arq).readline().strip().split(",")
    d = np.loadtxt(arq, delimiter=",", skiprows=1, usecols=[cab.index(c) for c in cols])
    x, y = d[:, 0], d[:, 1]
    # pontos repetidos (arestas e vertices comuns): fica um de cada
    _, u = np.unique(np.round(np.c_[x, y], 6), axis=0, return_index=True)
    d = d[np.sort(u)]
    x, y = d[:, 0], d[:, 1]
    tri = Triangulation(x, y)
    xc, yc = x[tri.triangles].mean(axis=1), y[tri.triangles].mean(axis=1)
    fora = (xc + yc > 20. + 1e-9) & (yc > 5. + 1e-9)
    # arestas longas da triangulacao de Delaunay no canto reentrante
    tri.set_mask(fora)
    norm = TwoSlopeNorm(vcenter=1.0, vmin=0.85, vmax=1.15)
    fig, axs = plt.subplots(2, 2, figsize=(E.LARGURA, 3.55))
    paineis = [("s2_ret_M200", "(a) forma fechada de R, M = 200"), ("s2_num_M200", "(b) Galerkin P2 h025, M = 200"),
               ("s2_ret_M2000", "(c) forma fechada de R, M = 2000"), ("s2_num_M2000", "(d) Galerkin P2 h025, M = 2000")]
    for ax, (col, tit) in zip(axs.ravel(), paineis):
        v = d[:, cols.index(col)]
        im = ax.tripcolor(tri, v, cmap=E.DIVERGENTE, norm=norm, shading="gouraud", rasterized=True)
        E.contorno_talude(ax, lw=0.6)
        ax.set_aspect("equal")
        ax.grid(False)
        ax.set_xlim(-0.2, 25.2)
        ax.set_ylim(-0.2, 10.2)
        ax.set_xticks([0, 5, 10, 15, 20, 25])
        ax.set_yticks([0, 5, 10])
        ax.tick_params(labelsize=6.5, length=2)
        ax.set_title("%s: %s a %s" % (tit, E.num(v.min(), 3), E.num(v.max(), 3)), loc="left", fontsize=7, pad=2)
    fig.subplots_adjust(right=0.88, hspace=0.35, wspace=0.12)
    cax = fig.add_axes([0.9, 0.15, 0.015, 0.7])
    cb = fig.colorbar(im, cax=cax, extend="both")
    cb.set_label(r"$\sigma_M^2$", fontsize=8)
    cb.ax.tick_params(labelsize=6.5)
    cb.ax.yaxis.set_major_formatter(E.fmt_virgula(2))
    E.salva(fig, "kl_variancia", png=True, dpi=200)


def modos_talude(lam_ny, nmodos=9):
    """os 9 primeiros modos da KL exata do talude (Nystrom, dados/kl_nystrom_modos.csv, grade de 0,25 m), no mesmo
    formato de kl_modos_retangulo.png (analitica.py)"""
    d = np.genfromtxt(os.path.join(E.DADOS, "kl_nystrom_modos.csv"), delimiter=",", names=True)
    xs, ys = np.unique(d["x"]), np.unique(d["y"])
    ix = np.searchsorted(xs, d["x"])
    iy = np.searchsorted(ys, d["y"])
    fig, axs = plt.subplots(3, 3, figsize=(E.LARGURA, 3.55), sharex=True, sharey=True)
    for k, ax in enumerate(axs.flat[:nmodos]):
        Z = np.full((len(ys), len(xs)), np.nan)
        Z[iy, ix] = d["phi%d" % (k + 1)]
        imx = np.unravel_index(np.nanargmax(np.abs(Z)), Z.shape)
        if Z[imx] < 0:
            Z = -Z
        v = np.nanmax(np.abs(Z))
        ax.pcolormesh(xs, ys, np.ma.masked_invalid(Z), cmap=E.DIVERGENTE, vmin=-v, vmax=v, shading="gouraud",
                      rasterized=True)
        E.contorno_talude(ax, lw=0.7)
        ax.set_aspect("equal")
        ax.grid(False)
        ax.set_xlim(0, 25)
        ax.set_ylim(0, 10)
        ax.set_title(r"$k$ = %d, $\lambda_k$ = %s" % (k + 1, E.num(lam_ny[k], casas=2)), fontsize=7, loc="left", pad=2)
        ax.set_xticks([0, 10, 15, 25])
        ax.set_yticks([0, 5, 10])
        ax.tick_params(labelsize=6.5, length=2)
    for ax in axs[-1]:
        ax.set_xlabel("x (m)", fontsize=7)
    for ax in axs[:, 0]:
        ax.set_ylabel("y (m)", fontsize=7)
    fig.tight_layout(w_pad=0.3, h_pad=0.25)
    E.salva(fig, "kl_modos_talude", png=True)


def grava_dados(kl, lam_ny, inc):
    """dados/kl_truncamento.csv (epsilon_M por M) e dados/kl_espectros.csv (lambda_k por k) das figuras"""
    casos = CASOS_ELEM + CASOS_O10
    lams = {nome: autovalores(nome) for nome, _, _, _ in casos}
    eps_ret = kl.erro_talude()
    Ms = sorted(set(list(range(1, 101)) + list(range(110, 1001, 10)) + list(range(1100, 20001, 100))))
    with open(os.path.join(E.DADOS, "kl_truncamento.csv"), "w") as f:
        f.write("M,forma_fechada_R_no_talude,nystrom_talude," + ",".join("galerkin_" + n for n, _, _, _ in casos) + "\n")
        for m in Ms:
            lin = [str(m), "%.6g" % eps_ret[m - 1], "%.6g" % (1. - lam_ny[:m].sum() / AREA_TALUDE) if m <= len(lam_ny) else ""]
            for nome, _, _, _ in casos:
                lam = lams[nome]
                lin.append("%.6g" % (1. - lam[:m].sum() / AREA_TALUDE) if m <= len(lam) else "")
            f.write(",".join(lin) + "\n")
    K = 6000
    with open(os.path.join(E.DADOS, "kl_espectros.csv"), "w") as f:
        f.write("k,retangulo_forma_fechada,nystrom_talude,incerteza_nystrom," +
                ",".join("galerkin_" + n for n, _, _, _ in casos) + "\n")
        for k in range(1, K + 1):
            lin = [str(k), "%.8g" % kl.lam[k - 1]]
            lin += ["%.8g" % lam_ny[k - 1], "%.2g" % inc[k - 1]] if (k <= len(lam_ny) and inc is not None) else ["", ""]
            for nome, _, _, _ in casos:
                lam = lams[nome]
                lin.append("%.8g" % lam[k - 1] if k <= len(lam) else "")
            f.write(",".join(lin) + "\n")
    print("  dados/kl_truncamento.csv, dados/kl_espectros.csv")
    if inc is not None:
        print("  incerteza da referencia, k = 50 a 100: max %.1e" % inc[49:100].max())


if __name__ == "__main__":
    etapas = sys.argv[1:] or ["truncamento", "espectros", "variancia", "modos", "dados"]
    lam_ny, inc, rot_ny = nystrom()
    print("  referencia:", rot_ny, "lambda_1..3 =", lam_ny[:3])
    kl = KL2D(20000) if ("truncamento" in etapas or "espectros" in etapas) else None
    if "truncamento" in etapas:
        truncamento(kl, lam_ny)
    if "espectros" in etapas:
        espectros(kl, lam_ny, inc)
    if "variancia" in etapas:
        variancia()
    if "modos" in etapas:
        modos_talude(lam_ny)
    if "dados" in etapas:
        grava_dados(kl if kl is not None else KL2D(20000), lam_ny, inc)
