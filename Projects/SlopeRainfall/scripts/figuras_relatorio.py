"""Figuras adicionais do relatório (docs/biot_nao_saturado): comparação das três variantes do modelo na Série B e campos
de p_w e p_a em instantes escolhidos.

Uso: python figuras_relatorio.py [--resultados resultados] [--figuras figuras]
"""
import argparse
import glob
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt     # noqa: E402
import matplotlib.tri as mtri       # noqa: E402
import numpy as np                  # noqa: E402

import ref_artigo as ref            # noqa: E402

plt.rcParams.update({"font.size": 9, "axes.linewidth": 0.6, "legend.fontsize": 7.5, "axes.titlesize": 9.5})
ROT_SOLO = {"f10": "f$_{10,-4}$", "f50": "f$_{50,-5}$", "f100": "f$_{100,-6}$"}


def resumo(pasta):
    s = {}
    for linha in open(os.path.join(pasta, "summary.txt")):
        k, v = linha.split(None, 1)
        try:
            s[k] = float(v)
        except ValueError:
            s[k] = v.strip()
    return s


def variantes_B(raiz, saida):
    series = [("B", "acoplado, ar ativo", "s", "#1b6ca8"), ("B_passivo", "acoplado, $p_a = 0$", "^", "#e9a23b"),
              ("B_passivo_desacoplado", "desacoplado, $p_a = 0$ (hipóteses do artigo)", "o", "#c8553d")]
    fig, axs = plt.subplots(1, 3, figsize=(14, 4.2), sharey=True)
    for ax, solo in zip(axs, ("f10", "f50", "f100")):
        Ir = np.logspace(-1, 3.2, 200)
        ax.semilogx(Ir, ref.fs_min_eq6(solo, Ir), "k-", lw=1.2, label="Rahardjo et al. (2007), Eq. 6")
        for s, rot, m, c in series:
            pts = []
            for f in glob.glob(os.path.join(raiz, s, f"{solo}_*", "summary.txt")):
                d = resumo(os.path.dirname(f))
                pts.append((d["I"], d["fs_min_srm"], d["t_fail"] >= 0 and d["fs_min_srm"] > 1.02))
            pts.sort()
            if not pts:
                continue
            x, y, inc = map(np.array, zip(*pts))
            ax.semilogx(x, y, m + "--", color=c, ms=5, lw=0.7, label=rot)
            if inc.any() and s.endswith("desacoplado"):
                ax.semilogx(x[inc], y[inc], m, mfc="w", mec=c, ms=9, label="percolação interrompida")
        ax.axhline(1, color="#999", lw=0.6, ls=":")
        ax.set(xlabel="I$_r$ (mm/h)", title=f"Série B, solo {ROT_SOLO[solo]}", ylim=(0.6, 2.3))
    axs[0].set_ylabel("FS mínimo")
    axs[0].legend(frameon=False, loc="lower left")
    fig.tight_layout()
    fig.savefig(os.path.join(saida, "B_variantes_fsmin.png"), dpi=160)
    plt.close(fig)


def campos(pasta, arquivo, instantes, titulo, saida, vmin=-80., vmax=None, rotulo="p$_w$ (kPa)", cmap="RdYlBu_r"):
    v = np.loadtxt(os.path.join(pasta, "vertices.csv"), delimiter=",", skiprows=1)
    t = np.loadtxt(os.path.join(pasta, "triangles.csv"), delimiter=",", skiprows=1, dtype=int)
    o = np.loadtxt(os.path.join(pasta, "outline.csv"), delimiter=",", skiprows=1)
    S = np.loadtxt(os.path.join(pasta, arquivo), delimiter=",")
    tri = mtri.Triangulation(v[:, 0], v[:, 1], t)
    vmax = vmax or 50. * np.ceil(S[:, 1:].max() / 50.)
    niveis = np.arange(vmin, vmax + 1e-9, 5. if vmax - vmin < 120 else 10.)
    fig, axs = plt.subplots(2, (len(instantes) + 1) // 2, figsize=(13, 5.2))
    for ax, ti in zip(axs.ravel(), instantes):
        k = int(np.argmin(np.abs(S[:, 0] - ti)))
        p = S[k, 1:]
        cf = ax.tricontourf(tri, np.clip(p, vmin, vmax), levels=niveis, cmap=cmap)
        if p.min() < 0 < p.max() and arquivo == "snaps_pw.csv":
            ax.tricontour(tri, p, levels=[0.], colors="magenta", linewidths=1.5)
        ax.plot(np.r_[o[:, 0], o[0, 0]], np.r_[o[:, 1], o[0, 1]], "k-", lw=0.7)
        ax.set_aspect("equal")
        ax.set_title(f"t = {S[k, 0]:g} h", fontsize=9)
        ax.set_xticks([])
        ax.set_yticks([])
    fig.colorbar(cf, ax=axs.ravel().tolist(), shrink=0.85, label=rotulo)
    fig.suptitle(titulo)
    fig.savefig(os.path.join(saida, arquivo.replace("snaps", "campos").replace(".csv", "") + "_" +
                             os.path.basename(pasta.rstrip("/")) + "_" + os.path.basename(os.path.dirname(pasta.rstrip("/")))
                             + ".png"), dpi=150)
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--resultados", default="resultados")
    p.add_argument("--figuras", default="figuras")
    a = p.parse_args()
    os.makedirs(a.figuras, exist_ok=True)
    variantes_B(a.resultados, a.figuras)
    inst = (0, 6, 12, 24, 36, 48)
    c1 = os.path.join(a.resultados, "B", "f50_a45_H10_Hw5_I36")
    c2 = os.path.join(a.resultados, "B_passivo_desacoplado", "f50_a45_H10_Hw5_I90")
    if os.path.exists(os.path.join(c1, "snaps_pw.csv")):
        campos(c1, "snaps_pw.csv", inst, "Modelo completo, f$_{50,-5}$, 36 mm/h: p$_w$ (magenta: p$_w$ = 0)", a.figuras)
        campos(c1, "snaps_pa.csv", inst, "Modelo completo, f$_{50,-5}$, 36 mm/h: pressão do ar p$_a$", a.figuras,
               vmin=0., vmax=55., rotulo="p$_a$ (kPa)", cmap="viridis")
    if os.path.exists(os.path.join(c2, "snaps_pw.csv")):
        campos(c2, "snaps_pw.csv", inst, "Desacoplado, p$_a$ = 0, f$_{50,-5}$, 90 mm/h: p$_w$ (magenta: p$_w$ = 0)",
               a.figuras)


if __name__ == "__main__":
    main()
