"""Figuras de comparação com Rahardjo et al. (2007) a partir de resultados/<série>/<caso>/ (SlopeRainfall).

O FS do artigo é o de Bishop simplificado (SEEP/W + SLOPE/W); o deste modelo é o da redução de resistência (SRM)
com as pressões do instante congeladas. O FS por aumento da gravidade (GI) aparece tracejado onde foi calculado.

Uso: python figuras.py --serie B          (A, B, C, D, B_passivo ou todas; grava figuras/<série>_*.png)
"""
import argparse
import glob
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt   # noqa: E402
import numpy as np                # noqa: E402

import ref_artigo as ref          # noqa: E402

plt.rcParams.update({"font.size": 9, "axes.linewidth": 0.6, "legend.fontsize": 7.5, "axes.titlesize": 9.5})
CORES = ["#1b6ca8", "#2a9d8f", "#e9a23b", "#c8553d", "#6c3483", "#4c9be8", "#7a7a7a", "#2e7d32"]
SOLOS = ("f10", "f50", "f100")
ROT_SOLO = {"f10": "f$_{10,-4}$", "f50": "f$_{50,-5}$", "f100": "f$_{100,-6}$"}


def carrega(raiz, serie):
    out = []
    for d in sorted(glob.glob(os.path.join(raiz, serie, "*", "summary.txt"))):
        pasta = os.path.dirname(d)
        s = {}
        for linha in open(d):
            k, v = linha.split(None, 1)
            v = v.strip()
            try:
                s[k] = float(v)
            except ValueError:
                s[k] = v
        fs = np.atleast_2d(np.loadtxt(os.path.join(pasta, "fs.csv"), delimiter=",", skiprows=1))
        s["t"], s["fs"], s["fsg"] = fs[:, 0], fs[:, 1], fs[:, 2]
        s["pasta"] = pasta
        out.append(s)
    return out


def fs_min(d):
    return d["fs"][d["fs"] > 0].min()


def serie_A(R, pasta):
    # Fig. 4: FS(t) por solo e ângulo
    fig, axs = plt.subplots(3, 4, figsize=(14, 8.5), sharex=True, sharey=True)
    for i, solo in enumerate(SOLOS):
        for j, a in enumerate(ref_alphas()):
            ax = axs[i, j]
            rs = sorted([d for d in R if d["solo"] == solo and abs(d["alpha"] - a) < 1e-6], key=lambda d: d["I"])
            for k, d in enumerate(rs):
                ax.plot(d["t"], d["fs"], "-o", ms=2.5, color=CORES[k], label=f"{d['I']:g} mm/h")
            ax.axhline(1, color="#999", lw=0.6, ls=":")
            ax.set(xlim=(0, 240), ylim=(0.5, 3.0), title=f"{ROT_SOLO[solo]}, α = {a:g}°")
            if rs:
                ax.legend(frameon=False, fontsize=6.5)
            if i == 2:
                ax.set_xlabel("t (h)")
            if j == 0:
                ax.set_ylabel("FS (SRM)")
    fig.suptitle("Série A (compare com a Fig. 4 do artigo)")
    fig.tight_layout()
    fig.savefig(os.path.join(pasta, "A_fs_tempo.png"), dpi=150)
    plt.close(fig)
    # Fig. 9 e Tabela 7
    fig, axs = plt.subplots(1, 3, figsize=(13, 4), sharey=True)
    al = np.linspace(20, 70, 50)
    for ax, solo in zip(axs, SOLOS):
        ax.plot(al, ref.FS_INI[0] + ref.FS_INI[1] * al, "k-", lw=1, label="FS inicial, Tabela 7")
        rs = [d for d in R if d["solo"] == solo]
        if rs:
            a0 = sorted({d["alpha"] for d in rs})
            ax.plot(a0, [next(d for d in rs if d["alpha"] == a)["fs"][0] for a in a0], "ko", mfc="w",
                    label="FS inicial, NeoPZ")
        for k, I in enumerate(sorted({d["I"] for d in rs})):
            c = CORES[k]
            if (solo, I) in ref.TABELA7:
                E, F = ref.TABELA7[(solo, I)]
                ax.plot(al, E + F * al, "-", color=c, lw=1, label=f"{I:g} mm/h, Tabela 7")
            if (solo, I) in ref.FIG9:
                ax.plot(ref.FIG9_ALPHA, ref.FIG9[(solo, I)], "x--", color=c, lw=0.8, ms=6,
                        label=f"{I:g} mm/h, Fig. 9")
            pts = sorted((d["alpha"], fs_min(d)) for d in rs if d["I"] == I)
            ax.plot(*zip(*pts), "s", color=c, label=f"{I:g} mm/h, NeoPZ")
        ax.axhline(1, color="#999", lw=0.6, ls=":")
        ax.set(xlabel="α (graus)", title=f"Série A, solo {ROT_SOLO[solo]}: FS mínimo (Fig. 9)", ylim=(0.2, 3.0))
        ax.legend(frameon=False, fontsize=5.5, ncol=2, loc="lower left")
    axs[0].set_ylabel("FS")
    fig.tight_layout()
    fig.savefig(os.path.join(pasta, "A_fs_min_alpha.png"), dpi=150)
    plt.close(fig)


def ref_alphas():
    return (26.6, 33.7, 45.0, 63.4)


def serie_B(R, pasta, nome="B"):
    fig, axs = plt.subplots(2, 3, figsize=(14, 7.5))
    f5 = ref.fig5b()
    for j, solo in enumerate(SOLOS):
        rs = sorted([d for d in R if d["solo"] == solo], key=lambda d: d["I"])
        ax = axs[0, j]
        ax.axvspan(0, 24, color="#e3effb", lw=0)
        for k, d in enumerate(rs):
            I = d["I"]
            ax.plot(d["t"], d["fs"], "-", color=CORES[k % 8], lw=1.3, label=f"{I:g} mm/h")
            if (d["fsg"] > 0).any():
                ax.plot(d["t"], d["fsg"], "--", color=CORES[k % 8], lw=0.7)
            if solo == "f50" and I in f5:
                ax.plot(*f5[I], "o", mfc="w", mec=CORES[k % 8], ms=4)
        ax.axhline(1, color="#999", lw=0.6, ls=":")
        ax.set(xlabel="t (h)", ylabel="FS", xlim=(0, 48), ylim=(0.5, 2.4),
               title=f"Série {nome}, {ROT_SOLO[solo]} (Fig. 5)" + (" — círculos: Fig. 5b" if solo == "f50" else ""))
        ax.legend(frameon=False, ncol=2, fontsize=6.5)
        # Fig. 7: FS mínimo x I_r com a Eq. 6
        ax = axs[1, j]
        Ir = np.logspace(-1, 3.2, 200)
        ax.semilogx(Ir, ref.fs_min_eq6(solo, Ir), "k-", lw=1, label="Eq. 6 (Fig. 7)")
        if rs:
            ax.semilogx([d["I"] for d in rs], [fs_min(d) for d in rs], "s", color=CORES[0], label="NeoPZ (SRM)")
            ax.semilogx([d["I"] for d in rs], [d["fs"][0] for d in rs], "o", mfc="w", color="k", ms=3,
                        label="NeoPZ, FS inicial")
        ax.axhline(1, color="#999", lw=0.6, ls=":")
        ax.set(xlabel="I$_r$ (mm/h)", ylabel="FS mínimo", ylim=(0.5, 2.4), title=f"{ROT_SOLO[solo]}: FS mínimo × I$_r$")
        ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(pasta, f"{nome}_fs.png"), dpi=150)
    plt.close(fig)


def serie_C(R, pasta):
    fig, ax = plt.subplots(figsize=(6, 4.2))
    ax.plot(ref.FIG8_HW, ref.FIG8["ini"], "k-", lw=1, label="FS inicial, Fig. 8")
    hw = sorted({d["hw"] for d in R})
    ini = [next(d for d in R if d["hw"] == h)["fs"][0] for h in hw]
    ax.plot(hw, ini, "ko", mfc="w", label="FS inicial, NeoPZ")
    for k, I in enumerate((9.0, 36.0, 80.0)):
        ax.plot(ref.FIG8_HW, ref.FIG8[I], "-", color=CORES[k], lw=1, label=f"{I:g} mm/h, Fig. 8")
        pts = sorted((d["hw"], fs_min(d)) for d in R if d["I"] == I)
        if pts:
            ax.plot(*zip(*pts), "s", color=CORES[k], label=f"{I:g} mm/h, NeoPZ")
    ax.axhline(1, color="#999", lw=0.6, ls=":")
    ax.set(xlabel="H$_w$ (m)", ylabel="FS", title="Série C: FS inicial e mínimo × H$_w$ (Fig. 8)", ylim=(0.5, 3.0))
    ax.legend(frameon=False, fontsize=6.5, ncol=2)
    fig.tight_layout()
    fig.savefig(os.path.join(pasta, "C_fs_hw.png"), dpi=150)
    plt.close(fig)


def serie_D(R, pasta):
    fig, ax = plt.subplots(figsize=(6, 4.2))
    ax.plot(ref.FIG10_HS, ref.FIG10["ini"], "k-", lw=1, label="FS inicial, Fig. 10")
    hs = sorted({d["hs"] for d in R})
    ax.plot(hs, [next(d for d in R if d["hs"] == h)["fs"][0] for h in hs], "ko", mfc="w", label="FS inicial, NeoPZ")
    for k, I in enumerate((9.0, 36.0, 80.0)):
        ax.plot(ref.FIG10_HS, ref.FIG10[I], "-", color=CORES[k], lw=1, label=f"{I:g} mm/h, Fig. 10")
        pts = sorted((d["hs"], fs_min(d)) for d in R if d["I"] == I)
        if pts:
            ax.plot(*zip(*pts), "s", color=CORES[k], label=f"{I:g} mm/h, NeoPZ")
    ax.axhline(1, color="#999", lw=0.6, ls=":")
    ax.set(xlabel="H$_s$ (m)", ylabel="FS", title="Série D: FS inicial e mínimo × H$_s$ (Fig. 10)", ylim=(0.5, 3.5))
    ax.legend(frameon=False, fontsize=6.5, ncol=2)
    fig.tight_layout()
    fig.savefig(os.path.join(pasta, "D_fs_hs.png"), dpi=150)
    plt.close(fig)


def tabela(R, serie, pasta):
    with open(os.path.join(pasta, f"{serie}_tabela.csv"), "w") as f:
        f.write("caso,solo,alpha,hs,hw,I,fs_ini_srm,fs_min_srm,t_min_srm,fs_ini_gi,fs_min_gi,t_fail,minutos\n")
        for d in R:
            f.write(f"{os.path.basename(d['pasta'])},{d['solo']},{d['alpha']:g},{d['hs']:g},{d['hw']:g},{d['I']:g},"
                    f"{d['fs'][0]:.3f},{fs_min(d):.3f},{d['t_min_srm']:g},{d['fs_ini_gi']:.3f},{d['fs_min_gi']:.3f},"
                    f"{d['t_fail']:g},{d['minutes']:.1f}\n")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--serie", default="B")
    p.add_argument("--resultados", default="resultados")
    p.add_argument("--figuras", default="figuras")
    a = p.parse_args()
    os.makedirs(a.figuras, exist_ok=True)
    for s in (["A", "B", "C", "D", "B_passivo", "B_passivo_desacoplado"] if a.serie == "todas" else [a.serie]):
        R = carrega(a.resultados, s)
        if not R:
            print(f"série {s}: sem resultados")
            continue
        if s.startswith("B"):
            serie_B(R, a.figuras, s)
        else:
            {"A": serie_A, "C": serie_C, "D": serie_D}[s[0]](R, a.figuras)
        tabela(R, s, a.figuras)
        print(f"série {s}: {len(R)} casos")


if __name__ == "__main__":
    main()
