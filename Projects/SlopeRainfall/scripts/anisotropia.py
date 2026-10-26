"""Estudo da permeabilidade anisotrópica K = diag(k_x, k_y): impacto no FS do talude sob chuva.

Caso base: solo f50 (k_s = 1e-5 m/s), alpha = 45 graus, Hs = 10 m, Hw = 5 m, chuva de 36 mm/h (= 1 k_s) por 24 h,
48 h simuladas; FS por redução de resistência e por aumento da gravidade. Varia-se k_x e k_y (a curva de
retenção e k_rw(S) não mudam; k_x = k_y = k_s é o caso isotrópico do artigo).

Uso:
  python anisotropia.py rodar --jobs 4 --exe <build>/Projects/SlopeRainfall/SlopeRainfall
  python anisotropia.py figura                    # figuras/anisotropia_fs.png e anisotropia_tabela.csv
  python anisotropia.py rodar --casos 1e-5,1e-5 1e-4,1e-5 --chuva 80
"""
import argparse
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np

CASOS = [(1e-5, 1e-5), (3e-5, 1e-5), (1e-4, 1e-5), (1e-5, 1e-6), (1e-6, 1e-5), (1e-5, 1e-4)]
TFS = (3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 36, 48)


def rotulo(kx, ky, I):
    return f"f50_a45_I{I:g}_kx{kx:.0e}_ky{ky:.0e}".replace("e-0", "e-")


def rodar(a):
    casos = [tuple(map(float, c.split(","))) for c in a.casos] if a.casos else CASOS
    tarefas = []
    for kx, ky in casos:
        out = os.path.join(a.saida, rotulo(kx, ky, a.chuva))
        if os.path.exists(os.path.join(out, "summary.txt")) and not a.refazer:
            continue
        os.makedirs(out, exist_ok=True)
        args = [a.exe, "solo=f50", "alpha=45", "hs=10", "hw=5", f"I={a.chuva:g}", f"kx={kx:g}", f"ky={ky:g}",
                "tfs=" + ",".join(f"{t:g}" for t in TFS), f"out={out}", f"gi={0 if a.sem_gi else 1}"]
        tarefas.append((out, args))
    print(f"{len(tarefas)} casos")

    def roda(t):
        out, args = t
        with open(os.path.join(out, "log.txt"), "w") as log:
            return out, subprocess.run(args, stdout=log, stderr=subprocess.STDOUT).returncode

    with ThreadPoolExecutor(max_workers=a.jobs) as ex:
        for f in as_completed([ex.submit(roda, t) for t in tarefas]):
            print(*f.result(), flush=True)


def ler(pasta):
    s = {}
    for linha in open(os.path.join(pasta, "summary.txt")):
        k, v = linha.split(None, 1)
        try:
            s[k] = float(v)
        except ValueError:
            s[k] = v.strip()
    fs = np.atleast_2d(np.loadtxt(os.path.join(pasta, "fs.csv"), delimiter=",", skiprows=1))
    h = np.atleast_2d(np.loadtxt(os.path.join(pasta, "hist.csv"), delimiter=",", skiprows=1))
    s.update(t=fs[:, 0], fs=fs[:, 1], fsg=fs[:, 2], hist=h, pasta=pasta)
    return s


def figura(a):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    R = []
    for d in sorted(os.listdir(a.saida)):
        p = os.path.join(a.saida, d)
        if os.path.exists(os.path.join(p, "summary.txt")):
            R.append(ler(p))
    R.sort(key=lambda r: (r["ky"], r["kx"]))
    os.makedirs(a.figuras, exist_ok=True)
    fig, axs = plt.subplots(1, 3, figsize=(15, 4.3))
    cores = plt.cm.viridis(np.linspace(0, 0.9, len(R)))
    for c, r in zip(cores, R):
        rot = f"k$_x$ = {r['kx']:.0e}, k$_y$ = {r['ky']:.0e}  (k$_x$/k$_y$ = {r['kx'] / r['ky']:g})"
        axs[0].plot(r["t"], r["fs"], "-o", ms=3, color=c, label=rot)
        if (r["fsg"] > 0).any():
            axs[1].plot(r["t"], r["fsg"], "-o", ms=3, color=c, label=rot)
        h = r["hist"]
        axs[2].plot(h[:, 0], h[:, 4] / np.maximum(h[:, 3], 1e-30) * (h[:, 2] > 0), "-", color=c, label=rot)
    for ax, tit in zip(axs[:2], ("FS por redução de resistência", "FS por aumento da gravidade")):
        ax.axvspan(0, 24, color="#e3effb", lw=0)
        ax.axhline(1, color="#999", lw=0.6, ls=":")
        ax.set(xlabel="t (h)", ylabel="FS", title=tit, xlim=(0, 48))
    axs[2].set(xlabel="t (h)", ylabel="infiltração / chuva", title="Fração da chuva que infiltra", xlim=(0, 24),
               ylim=(0, 1.05))
    axs[0].legend(frameon=False, fontsize=7)
    I = R[0]["I"] if R else 0
    fig.suptitle(f"Solo f50 (k$_s$ = 1e-5 m/s), α = 45°, H$_s$ = 10 m, H$_w$ = 5 m, chuva {I:g} mm/h por 24 h")
    fig.tight_layout()
    fig.savefig(os.path.join(a.figuras, "anisotropia_fs.png"), dpi=150)
    with open(os.path.join(a.figuras, "anisotropia_tabela.csv"), "w") as f:
        f.write("kx,ky,kx_ky,fs_ini_srm,fs_min_srm,t_min_srm,fs_48h_srm,fs_ini_gi,fs_min_gi,t_min_gi,t_fail,"
                "infil_24h_m3m,pa_max_kPa,minutos\n")
        for r in R:
            h = r["hist"]
            infil = np.sum(h[:, 4] * h[:, 1] * 3600 * (h[:, 2] > 0))
            f.write(f"{r['kx']:g},{r['ky']:g},{r['kx'] / r['ky']:g},{r['fs'][0]:.3f},{r['fs_min_srm']:.3f},"
                    f"{r['t_min_srm']:g},{r['fs'][-1]:.3f},{r['fs_ini_gi']:.3f},{r['fs_min_gi']:.3f},"
                    f"{r['t_min_gi']:g},{r['t_fail']:g},{infil:.2f},{h[:, 7].max():.1f},{r['minutes']:.1f}\n")
    print(open(os.path.join(a.figuras, "anisotropia_tabela.csv")).read())


def main():
    aqui = os.path.dirname(os.path.abspath(__file__))
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("acao", choices=["rodar", "figura"])
    p.add_argument("--exe", default=os.path.join(aqui, "..", "..", "..", "..", "build", "Projects", "SlopeRainfall",
                                                  "SlopeRainfall"))
    p.add_argument("--saida", default=os.path.join("resultados", "anisotropia"))
    p.add_argument("--figuras", default="figuras")
    p.add_argument("--casos", nargs="*", help="pares kx,ky (m/s)")
    p.add_argument("--chuva", type=float, default=36.0)
    p.add_argument("--jobs", type=int, default=4)
    p.add_argument("--sem-gi", action="store_true")
    p.add_argument("--refazer", action="store_true")
    a = p.parse_args()
    rodar(a) if a.acao == "rodar" else figura(a)


if __name__ == "__main__":
    main()
