"""Figuras e tabelas do patamar inferior do FS e dos benchmarks com a literatura (Zhang et al. 2014, Tang et al.
2018, Pradhan et al. 2022).

Uso: python benchmarks.py --runs <pasta com patamar/, zhang/, tang/, pradhan/> --figuras ../figuras
"""
import argparse
import glob
import os

import numpy as np

# leituras dos artigos (as de figuras são aproximadas)
ZHANG = {"A": {"chi0": 0.898, "phib15": 1.414, "chi1": 1.920, "sw_swcc1": 1.032},
         "B": {"chi0": 1.086, "phib15": 1.133, "chi1": 1.159}}
# Tang et al. 2018: FS inicial 1,20 (Fig. 9); FSmin = 1 em Ia ~ 2,70 (uniforme, 25,6 h) e 9,70 mm/h (6,4 h; a
# probabilidade condicional de 43,2 % do texto cai em 9,65 mm/h na mesma figura); inclinação de FSmin x volume na
# Fig. 10: -0,0039/mm (25,6 h: 1,06 em 2,0 mm/h -> 0,91 em 3,5) e -0,0047/mm (6,4 h)
TANG = {"fs_ini": 1.20, "Ic": {25.6: 2.70, 6.4: 9.70}, "Ic_pad": {"adiantado": 2.77, "uniforme": 2.70, "atrasado": 2.55},
        "dfsdv": {25.6: -0.0039, 6.4: -0.0047}}
# Pradhan et al. 2022, Fig. 10 (leitura): FS com o corte e sem o corte
PRADHAN_CORTE = [(1, 1.79), (34, 1.75), (43, 1.46), (96, 1.71), (101, 1.42), (167, 1.82), (173, 1.46), (183, 1.62),
                 (190, 1.32), (192, 1.45), (200, 1.12), (215, 1.42), (222, 1.13), (260, 1.29), (264, 1.10),
                 (267, 0.97), (271, 1.00)]
PRADHAN_SEM = [(1, 2.75), (35, 2.72), (98, 2.71), (168, 2.72), (190, 2.70), (218, 2.68), (262, 2.65), (271, 2.62)]


def ler_fs(pasta):
    f = os.path.join(pasta, "fs.csv")
    if not os.path.exists(f):
        return None
    a = np.atleast_2d(np.loadtxt(f, delimiter=",", skiprows=1))
    if a.size == 0:  # case still running: header only
        return None
    return a[:, 0], a[:, 1]


def resumo(pasta):
    s = {}
    f = os.path.join(pasta, "summary.txt")
    if os.path.exists(f):
        for l in open(f):
            k, v = l.split(None, 1)
            try:
                s[k] = float(v)
            except ValueError:
                s[k] = v.strip()
    return s


def patamar(a, plt):
    casos = [("f50_iso", "isotrópico $f_{50}$", 1e-5, 0.84), ("f50_th135", r"$\theta=135^\circ$ ($k_1$ paralelo à face)", 1e-5, 0.63),
             ("f50_th45", r"$\theta=45^\circ$", 1e-5, 0.87), ("f50_th0", r"$\theta=0^\circ$ (horizontal)", 1e-5, 1.05),
             ("f10_iso", "isotrópico $f_{10}$", 1e-4, 0.84), ("f100_iso", "isotrópico $f_{100}$", 1e-6, None)]
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.3))
    cores = ["#1b6ca8", "#c8553d", "#e9a23b", "#2a9d8f", "#6c3483", "#555555"]
    linhas = []
    for (nome, rot, ks, plat), cor in zip(casos, cores):
        r = ler_fs(os.path.join(a.runs, "patamar", nome))
        if r is None:
            continue
        t, fs = r
        m = t > 0
        axs[0].semilogx(t[m], fs[m], "o-", color=cor, label=rot)
        axs[1].semilogx(t[m] * ks * 3600 / 10.0, fs[m], "o-", color=cor, label=rot)
        if plat:
            axs[0].axhline(plat, color=cor, ls=":", lw=0.8)
        linhas.append((nome, fs[0], fs[m].min(), t[m][np.argmin(fs[m])], fs[m][-1], t[m][-1], plat))
    axs[0].set(xlabel="tempo com a superfície alagada (h)", ylabel="FS (SRM)",
               title="Superfície alagada ($I\\to\\infty$); pontilhado: patamar com chuva (A3, série B)")
    axs[1].set(xlabel="$k_s t / H_s$ (k_s: condutividade de referência)", ylabel="FS (SRM)", title="Tempo adimensional")
    for ax in axs:
        ax.axhline(1, color="#999", lw=0.6, ls="--")
        ax.legend(frameon=False, fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(a.figuras, "patamar_alagado.png"), dpi=150)
    plt.close(fig)
    with open(os.path.join(a.figuras, "patamar_tabela.csv"), "w") as f:
        f.write("caso,fs_ini,fs_min,t_min_h,fs_final,t_final_h,patamar_chuva\n")
        for l in linhas:
            f.write(",".join(f"{v}" if isinstance(v, str) else (f"{v:.4g}" if v is not None else "") for v in l) + "\n")
    return linhas


def zhang(a):
    out = []
    for caso, ref in ZHANG.items():
        for k, v in ref.items():
            s = resumo(os.path.join(a.runs, "zhang", f"{caso}_{k}"))
            out.append((caso, k, v, s.get("fs_ini_srm")))
    s = resumo(os.path.join(a.runs, "zhang", "A_chi1_smax75"))
    out.append(("A", "chi1_smax75", None, s.get("fs_ini_srm")))
    with open(os.path.join(a.figuras, "zhang_tabela.csv"), "w") as f:
        f.write("caso,lei,fs_artigo,fs_neopz\n")
        for c in out:
            f.write(",".join("" if v is None else (f"{v:.4g}" if not isinstance(v, str) else v) for v in c) + "\n")
    return out


def tang(a, plt):
    pasta = os.path.join(a.runs, "tang")
    fig, axs = plt.subplots(1, 3, figsize=(15, 4.2))
    res = {}
    for Dr, ax, Is in ((25.6, axs[0], (2.0, 2.4, 2.7, 3.0, 3.5)), (6.4, axs[1], (8.0, 10.0, 12.0))):
        pts = []
        for I in Is:
            s = resumo(os.path.join(pasta, f"u{Dr:g}_I{I:g}"))
            if "fs_min_srm" in s:
                pts.append((I, s["fs_min_srm"], s["fs_ini_srm"]))
        if pts:
            I_, fm, fi = map(np.array, zip(*pts))
            ax.plot(I_, fm, "o-", color="#1b6ca8", label="NeoPZ (SRM)")
            res[Dr] = pts
        Ig = np.linspace(min(Is) - 0.5, max(Is) + 0.5, 50)
        ax.plot(Ig, 1 + TANG["dfsdv"][Dr] * Dr * (Ig - TANG["Ic"][Dr]), "k--", lw=0.8,
                label="Tang, Fig. 10 (FS = 1 em %.2f mm/h, dFS/dV = %.4f/mm)" % (TANG["Ic"][Dr], TANG["dfsdv"][Dr]))
        ax.axhline(1, color="#999", lw=0.6, ls=":")
        ax.set(xlabel="$I_a$ (mm/h)", ylabel="FS mínimo", title=f"Chuva uniforme de {Dr:g} h")
        ax.legend(frameon=False, fontsize=7)
    for pad, cor in (("adiantado", "#2a9d8f"), ("uniforme", "#1b6ca8"), ("atrasado", "#c8553d")):
        r = ler_fs(os.path.join(pasta, f"{pad}25.6_I2.64"))
        if r is not None:
            axs[2].plot(r[0], r[1], "o-", ms=3, color=cor, label=pad)
    r = ler_fs(os.path.join(pasta, "controle"))
    if r is not None:
        axs[2].plot(r[0], r[1], "k:", label="sem chuva")
    axs[2].axhline(1, color="#999", lw=0.6, ls=":")
    axs[2].set(xlabel="t (h)", ylabel="FS (SRM)", title="$I_a$ = 2,64 mm/h, 25,6 h: padrões")
    axs[2].legend(frameon=False, fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(a.figuras, "tang_fs.png"), dpi=150)
    plt.close(fig)
    return res


def pradhan(a, plt):
    pasta = os.path.join(a.runs, "pradhan")
    chuva = np.loadtxt(os.path.join(os.path.dirname(__file__), "..", "pradhan", "chuva_julho2018.csv"),
                       delimiter=",", comments="#")
    fig, ax = plt.subplots(figsize=(11, 4.3))
    ax2 = ax.twinx()
    ax2.bar(0.5 * (chuva[:, 0] + chuva[:, 1]), chuva[:, 2], width=1.0, color="#bbbbbb", edgecolor="#888888", lw=0.3)
    ax2.set_ylim(30, 0)
    ax2.set_ylabel("chuva (mm/h)")
    for nome, rot, cor, ref in (("corte", "com o corte (NeoPZ)", "#c8553d", PRADHAN_CORTE),
                                ("semcorte", "sem o corte (NeoPZ)", "#1b6ca8", PRADHAN_SEM)):
        r = ler_fs(os.path.join(pasta, nome))
        if r is not None:
            ax.plot(r[0], r[1], "o-", ms=3, color=cor, label=rot)
        t, f = zip(*ref)
        ax.plot(t, f, "x--", color=cor, lw=0.8, label=rot.replace("NeoPZ", "Pradhan, Fig. 10"))
    ax.axhline(1, color="#999", lw=0.6, ls=":")
    ax.set(xlabel="horas a partir de 01/07/2018", ylabel="FS", title="Pradhan et al. (2022): chuva de 1 a 12/07/2018")
    ax.set_xlim(0, 272)
    ax.legend(frameon=False, fontsize=7, loc="center left", bbox_to_anchor=(0.2, 0.64), ncol=2)
    fig.tight_layout()
    fig.savefig(os.path.join(a.figuras, "pradhan_fs.png"), dpi=150)
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--runs", default="runs")
    p.add_argument("--figuras", default="figuras")
    a = p.parse_args()
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 9, "axes.linewidth": 0.6, "legend.fontsize": 7.5, "axes.titlesize": 9.5})
    os.makedirs(a.figuras, exist_ok=True)
    for l in patamar(a, plt):
        print("patamar", l)
    for z in zhang(a):
        print("zhang", z)
    print("tang", tang(a, plt))
    pradhan(a, plt)


if __name__ == "__main__":
    main()
