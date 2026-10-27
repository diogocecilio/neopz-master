"""FS ao longo do tempo e animação das isolinhas de p_w (com o FS) de um caso do SlopeRainfall.

Equivalente ao pós-processamento em Mathematica (mma/animacao.wl):
  * fs_tempo.png: FS por redução de resistência e por aumento da gravidade em escala log(t), t = 0 desenhado em
    0,01 h, linha vertical vermelha no instante de ruptura do transiente (se houver);
  * chuva.mp4 / chuva.gif: em cada instante gravado (snaps_pw.csv), preenchimento de p_w a cada 10 kPa, isolinhas
    rotuladas a cada 50 kPa, p_w = 0 em magenta, e ao lado o gráfico do FS com um marcador no instante. O FS do
    marcador é interpolado linearmente em log(t) (constante se houver um só ponto).

Uso: python animacao.py <pasta do caso> [--pmin -80] [--pmax 0 (automático)] [--fps 4] [--sem-video]
"""
import argparse
import os
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt                     # noqa: E402
import matplotlib.tri as mtri                       # noqa: E402
from matplotlib import animation                    # noqa: E402
from matplotlib.colors import Normalize             # noqa: E402

T0 = 0.01   # t = 0 desenhado em 0,01 h


def ler_caso(pasta):
    v = np.loadtxt(os.path.join(pasta, "vertices.csv"), delimiter=",", skiprows=1)
    t = np.loadtxt(os.path.join(pasta, "triangles.csv"), delimiter=",", skiprows=1, dtype=int)
    o = np.loadtxt(os.path.join(pasta, "outline.csv"), delimiter=",", skiprows=1)
    sw = np.atleast_2d(np.loadtxt(os.path.join(pasta, "snaps_pw.csv"), delimiter=","))
    fs = np.atleast_2d(np.loadtxt(os.path.join(pasta, "fs.csv"), delimiter=",", skiprows=1))
    summ = {}
    p = os.path.join(pasta, "summary.txt")
    if os.path.exists(p):
        for linha in open(p):
            k, val = linha.split(None, 1)
            summ[k] = val.strip()
    return dict(xy=v, tri=t, outline=o, ts=sw[:, 0], pw=sw[:, 1:], fs=fs, summary=summ)


def interp_fs(t, f):
    """InterpFS: FS interpolado em log(t) (t em h, t = 0 em 0,01 h, fora dos dados: extremos)"""
    t = np.maximum(np.asarray(t, float), T0)
    ok = f > 0
    t, f = t[ok], f[ok]
    if len(t) < 2:
        c = f[0] if len(f) else np.nan
        return lambda x: c
    lt = np.log(t)
    return lambda x: float(np.interp(np.log(np.clip(max(x, T0), t.min(), t.max())), lt, f))


def grafico_fs(ax, c, titulo=True):
    fs = c["fs"]
    t = np.maximum(fs[:, 0], T0)
    linhas = []
    for col, cor, rot in ((1, "#1f3fa8", "redução de resistência"), (2, "#1d7a1d", "aumento da gravidade")):
        ok = fs[:, col] > 0
        if ok.any():
            linhas += ax.semilogx(t[ok], fs[ok, col], "-o", ms=4, color=cor, label=rot)
    vals = fs[:, 1:3][fs[:, 1:3] > 0]
    top = 1.05 * vals.max() if len(vals) else 3.
    ax.set_ylim(min(0.9, vals.min() - 0.05) if len(vals) else 0.9, top)
    ax.axhline(1.0, color="gray", ls="--", lw=0.8)
    ax.axvline(1.0, color="gray", ls="--", lw=0.8)
    tfail = float(c["summary"].get("t_fail", -1))
    if tfail >= 0:
        # ruptura só se o FS no instante for ~1; senão o passo acoplado apenas não convergiu (dificuldade numérica)
        f_t = interp_fs(fs[:, 0], fs[:, 1])(tfail)
        rot = " ruptura" if f_t < 1.05 else " transiente sem convergência"
        ax.axvline(max(tfail, T0), color="red" if f_t < 1.05 else "orange", lw=2)
        ax.text(max(tfail, T0), 1.0, rot, color="red" if f_t < 1.05 else "orange", va="bottom", fontsize=8)
    ax.set_xlabel("t (h)  [t = 0 desenhado em 0,01 h]")
    ax.set_ylabel("FS")
    ax.grid(False)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=2, frameon=False, fontsize=8)
    if titulo:
        s = c["summary"]
        ax.set_title(f"Chuva de {float(s.get('I', 0)):g} mm/h por {float(s.get('train', 24)):g} h: "
                     f"FS ao longo do tempo", fontsize=10)
    return linhas


def rotulo_t(t):
    return f"{t:.2f} h" if t < 1 else f"{t:.0f} h"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("pasta")
    ap.add_argument("--pmin", type=float, default=-80.)
    ap.add_argument("--pmax", type=float, default=0., help="0: máximo dos instantes, arredondado a 50 kPa")
    ap.add_argument("--fps", type=float, default=4.)
    ap.add_argument("--sem-video", action="store_true")
    a = ap.parse_args()
    c = ler_caso(a.pasta)
    s = c["summary"]
    I, train = float(s.get("I", 0)), float(s.get("train", 24))

    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    grafico_fs(ax, c)
    fig.tight_layout()
    fig.savefig(os.path.join(a.pasta, "fs_tempo.png"), dpi=130)
    plt.close(fig)
    if a.sem_video:
        return

    fsAt = interp_fs(c["fs"][:, 0], c["fs"][:, 1])
    fsGAt = interp_fs(c["fs"][:, 0], c["fs"][:, 2])
    tri = mtri.Triangulation(c["xy"][:, 0], c["xy"][:, 1], c["tri"])
    pmax = a.pmax if a.pmax > 0 else 50. * np.ceil(c["pw"].max() / 50.)
    pmin = a.pmin
    norm = Normalize(pmin, pmax)
    cmap = plt.get_cmap("RdYlBu_r")
    niveis = np.arange(pmin, pmax + 1e-9, 10.)
    linhas50 = np.arange(0., pmax + 1e-9, 50.)
    o = c["outline"]
    L, H = o[:, 0].max(), o[:, 1].max()

    fig = plt.figure(figsize=(13, 4.8))
    axp = fig.add_axes([0.05, 0.13, 0.50, 0.78])
    axc = fig.add_axes([0.565, 0.13, 0.012, 0.78])
    axf = fig.add_axes([0.685, 0.25, 0.30, 0.62])
    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    fig.colorbar(sm, cax=axc, label="p_w (kPa)")
    grafico_fs(axf, c, titulo=False)
    marc = axf.axvline(T0, color="red", ls="--")
    pts, = axf.plot([], [], "o", color="red", ms=7)
    titf = axf.set_title("", fontsize=10)

    def quadro(i):
        t, p = c["ts"][i], c["pw"][i]
        axp.clear()
        axp.tricontourf(tri, np.clip(p, pmin, pmax), levels=niveis, cmap=cmap, norm=norm)
        cs = axp.tricontour(tri, p, levels=linhas50, colors="k", linewidths=0.5)
        axp.clabel(cs, fmt="%.0f", fontsize=7)
        if p.min() < 0 < p.max():
            axp.tricontour(tri, p, levels=[0.], colors="magenta", linewidths=2)
        axp.plot(np.r_[o[:, 0], o[0, 0]], np.r_[o[:, 1], o[0, 1]], "k-", lw=0.8)
        chovendo = 0 < t <= train + 1e-9 or (t == 0 and I > 0)
        if chovendo and t > 0:
            xs = np.linspace(0.02 * L, 0.98 * L, 28)
            ys = np.interp(xs, np.r_[0, o[4, 0], o[3, 0], L], np.r_[o[5, 1], o[4, 1], o[3, 1], o[3, 1]])
            for x, y in zip(xs, ys):
                axp.annotate("", xy=(x, y + 0.4), xytext=(x, y + 0.08 * H),
                             arrowprops=dict(arrowstyle="->", color="tab:blue", lw=0.8))
        axp.set_xlim(0, L)
        axp.set_ylim(0, 1.13 * H)
        axp.set_aspect("equal")
        axp.set_xlabel("x (m)")
        axp.set_ylabel("y (m)")
        estado = f"chuva {I:g} mm/h" if chovendo and t > 0 else "sem chuva"
        axp.set_title(f"t = {rotulo_t(t)}    {estado}    (magenta: p_w = 0)", fontsize=10)
        th = max(t, T0)
        marc.set_xdata([th, th])
        f1, f2 = fsAt(t), fsGAt(t)
        pts.set_data([th, th] if np.isfinite(f2) else [th], [f1, f2] if np.isfinite(f2) else [f1])
        titf.set_text(f"FS = {f1:.2f}" + (f"    FSg = {f2:.2f}" if np.isfinite(f2) else ""))
        return []

    anim = animation.FuncAnimation(fig, quadro, frames=len(c["ts"]), blit=False)
    base = os.path.join(a.pasta, "chuva")
    try:
        anim.save(base + ".mp4", writer=animation.FFMpegWriter(fps=a.fps, bitrate=2400), dpi=96)
    except Exception as e:                                  # noqa: BLE001
        print("mp4 não gerado:", e, file=sys.stderr)
    anim.save(base + ".gif", writer=animation.PillowWriter(fps=a.fps), dpi=72)
    plt.close(fig)
    print("gerados:", base + ".mp4", base + ".gif", os.path.join(a.pasta, "fs_tempo.png"))


if __name__ == "__main__":
    main()
