#!/usr/bin/env python3
"""analitica.py: como funciona a KL em forma fechada (sem dados externos; so klfechada.py).

  kl_analitica_1d.pdf      (A1) as equacoes transcendentes par e impar com as raizes, as primeiras autofuncoes 1D e
                           lambda = 2c/(omega^2 + c^2), em x (L = 20 m, a = 12,5 m) e em y (L = 2 m, a = 5 m)
  kl_modos_retangulo.png   (A2) os 9 primeiros modos 2D f_i(x) g_j(y) no retangulo [0,25] x [0,10], com o talude
  dados/kl_modos_1d.csv    as primeiras raizes e autovalores 1D
"""
import os

import numpy as np

import estilo as E
from estilo import plt
from klfechada import KL2D, Modos1D, TALUDE

CASOS = [("x", 25., 20.), ("y", 10., 2.)]  # direcao, comprimento do retangulo, comprimento de correlacao


def figura_1d():
    fig, axs = plt.subplots(2, 3, figsize=(E.LARGURA, 4.3), gridspec_kw=dict(width_ratios=[1.15, 1.1, 0.95]))
    linhas_csv = []
    for lin, (d, comp, L) in enumerate(CASOS):
        m = Modos1D(comp, L, 40)
        a, c = m.a, m.c
        ca = c * a
        # (1) equacoes transcendentes em theta = omega a: par theta tan(theta) = c a; impar tan(theta) = -theta / (c a)
        ax = axs[lin, 0]
        tmax, ylim = (2.25, 9.) if d == "x" else (3., 4.)
        for k in range(int(np.ceil(tmax)) + 1):
            t0, t1 = (k - 0.5) * np.pi, min((k + 0.5) * np.pi, tmax * np.pi)
            t = np.linspace(max(t0, 0.) + 1e-4, t1 - 1e-4, 400)
            ax.plot(t / np.pi, np.tan(t), color=E.NEUTRO, lw=1.0, label=r"$\tan\theta$" if k == 0 else None)
        t = np.linspace(1e-3, tmax * np.pi, 600)
        ax.plot(t / np.pi, ca / t, color=E.PALETA[0], lw=1.4, label=r"$ca/\theta$: pares")
        ax.plot(t / np.pi, -t / ca, color=E.PALETA[1], lw=1.4, label=r"$-\theta/(ca)$: ímpares")
        th = m.omega * a
        sel = th < tmax * np.pi
        for par, cor, mk in ((True, E.PALETA[0], "o"), (False, E.PALETA[1], "s")):
            s = sel & (m.par == par) & (np.abs(np.tan(th)) < ylim)
            ax.plot(th[s] / np.pi, np.tan(th[s]), mk, color=cor, ms=4.5, mec="white", mew=0.6, zorder=5)
        ax.set_xlim(0, tmax)
        ax.set_ylim(-ylim, ylim)
        ax.set_xticks(np.arange(0, tmax + 0.01, 0.5))
        ax.xaxis.set_major_formatter(E.FuncFormatter(lambda v, _p: "0" if v == 0 else
                                                     (r"$\pi$" if v == 1 else (r"$%s\pi$" % E.num(v)))))
        ax.yaxis.set_major_formatter(E.fmt_virgula())
        ax.set_xlabel(r"$\theta = \omega a$")
        ax.set_title("%s: $L_%s$ = %s m, $a$ = %s m, $ca$ = %s" % (d, d, E.num(L), E.num(a), E.num(ca)), loc="left")
        if lin == 0:
            ax.legend(loc="upper right", fontsize=6.3, handlelength=1.5, borderaxespad=0.2, labelspacing=0.25,
                      facecolor="white", framealpha=0.9, frameon=True, edgecolor="none")
        # (2) as primeiras autofuncoes 1D no intervalo [0, D]
        ax = axs[lin, 1]
        s = np.linspace(-a, a, 400)
        for i in range(4):
            f = m.avalia(i, s)[0]
            k = np.argmax(np.abs(f))
            if f[k] < 0:
                f = -f
            ax.plot(s + a, f, color=E.PALETA[i], lw=1.3, ls="-" if m.par[i] else (0, (4, 1.6)),
                    label=r"$%s_%d$, $\lambda$ = %s" % ("f" if d == "x" else "g", i + 1, E.num(m.lam[i], sig=4)))
        ax.axhline(0, color=E.TINTA2, lw=0.6)
        ax.set_xlim(0, comp)
        ax.set_xlabel("%s (m)" % d)
        ax.set_ylabel(r"$%s_i(%s)$ (m$^{-1/2}$)" % ("f" if d == "x" else "g", d))
        E.virgula(ax)
        ylo, yhi = ax.get_ylim()
        ax.set_ylim(ylo - 0.75 * (yhi - ylo), yhi)
        ax.set_yticks([t for t in ax.get_yticks() if ylo - 1e-9 <= t <= yhi])
        ax.legend(loc="lower center", fontsize=6.2, handlelength=1.8, ncol=2, borderaxespad=0.2, labelspacing=0.25,
                  columnspacing=0.8)
        # (3) lambda = 2c / (omega^2 + c^2) e os autovalores nas raizes
        ax = axs[lin, 2]
        w = np.linspace(0, m.omega[11] * 1.05, 400)
        ax.plot(w, 2 * c / (w * w + c * c), color=E.NEUTRO, lw=1.0, label=r"$2c/(\omega^2+c^2)$")
        for par, cor, mk, nome in ((True, E.PALETA[0], "o", "pares"), (False, E.PALETA[1], "s", "ímpares")):
            s_ = (m.par == par) & (np.arange(len(m.lam)) < 12)
            ax.plot(m.omega[s_], m.lam[s_], mk, color=cor, ms=4.5, mec="white", mew=0.6, label=nome, zorder=5)
        ax.set_xlabel(r"$\omega$ (rad/m)")
        ax.set_ylabel(r"$\lambda$ (m)")
        ax.set_xlim(0, w[-1])
        ax.set_ylim(0, 2 / c * 1.04)
        E.virgula(ax)
        if lin == 0:
            ax.legend(loc="upper right", fontsize=6.5, handlelength=1.2, borderaxespad=0.2)
        for i in range(10):
            linhas_csv.append((d, comp, L, i + 1, "par" if m.par[i] else "impar", m.omega[i], m.omega[i] * a, m.lam[i],
                               m.norma[i]))
    fig.tight_layout(w_pad=0.6, h_pad=0.9)
    E.salva(fig, "kl_analitica_1d")
    with open(os.path.join(E.DADOS, "kl_modos_1d.csv"), "w") as f:
        f.write("direcao,comprimento_m,L_m,i,paridade,omega_rad_m,theta,lambda_m,norma\n")
        for r in linhas_csv:
            f.write("%s,%g,%g,%d,%s,%.12g,%.12g,%.12g,%.12g\n" % r)


def figura_modos_2d(nmodos=9):
    kl = KL2D(2000)
    xg = np.linspace(0, 25, 501)
    yg = np.linspace(0, 10, 201)
    fig, axs = plt.subplots(3, 3, figsize=(E.LARGURA, 3.55), sharex=True, sharey=True)
    for k, ax in enumerate(axs.flat[:nmodos]):
        Z = (kl.mx.avalia(kl.ix[k], xg - 12.5)[0][None, :] * kl.my.avalia(kl.iy[k], yg - 5.)[0][:, None])
        imx = np.unravel_index(np.argmax(np.abs(Z)), Z.shape)
        if Z[imx] < 0:
            Z = -Z
        v = np.abs(Z).max()
        ax.pcolormesh(xg, yg, Z, cmap=E.DIVERGENTE, vmin=-v, vmax=v, shading="gouraud", rasterized=True)
        fora = np.array([[25., 5.], [15., 5.], [10., 10.], [25., 10.]])
        ax.fill(fora[:, 0], fora[:, 1], color="white", alpha=0.55, lw=0)
        E.contorno_talude(ax, lw=0.7)
        ax.set_aspect("equal")
        ax.grid(False)
        ax.set_title(r"$k$ = %d: $(i,j)$ = (%d, %d), $\lambda_k$ = %s" % (k + 1, kl.ix[k] + 1, kl.iy[k] + 1,
                                                                        E.num(kl.lam[k], casas=2)),
                     fontsize=7, loc="left", pad=2)
        ax.set_xticks([0, 10, 15, 25])
        ax.set_yticks([0, 5, 10])
        ax.tick_params(labelsize=6.5, length=2)
    for ax in axs[-1]:
        ax.set_xlabel("x (m)", fontsize=7)
    for ax in axs[:, 0]:
        ax.set_ylabel("y (m)", fontsize=7)
    fig.tight_layout(w_pad=0.3, h_pad=0.25)
    E.salva(fig, "kl_modos_retangulo", png=True)
    return kl


if __name__ == "__main__":
    figura_1d()
    figura_modos_2d()
