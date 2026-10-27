"""Estudo da permeabilidade anisotrópica: impacto no FS do talude sob chuva.

K = R(theta) diag(k_x, k_y) R(theta)^T: k_x na direção a theta graus (anti-horário) do eixo x, k_y normal a ela.

Caso base: solo f50 (k_s = 1e-5 m/s), alpha = 45 graus, Hs = 10 m, Hw = 5 m, chuva de 36 mm/h (= 1 k_s) por 24 h,
48 h simuladas; FS por redução de resistência e por aumento da gravidade. Varia-se k_x e k_y (a curva de
retenção e k_rw(S) não mudam; k_x = k_y = k_s é o caso isotrópico do artigo).

Uso:
  python anisotropia.py rodar --jobs 4 --exe <build>/Projects/SlopeRainfall/SlopeRainfall
  python anisotropia.py figura                    # figuras/anisotropia_fs.png e anisotropia_tabela.csv
  python anisotropia.py rodar --casos 1e-5,1e-5 1e-4,1e-5 --chuva 80
  python anisotropia.py rotacao --jobs 4          # k_x = 1e-4, k_y = 1e-5 girados de 0 a 165 graus (15 em 15)
  python anisotropia.py figura-rotacao            # figuras/rotacao_fs.png e rotacao_tabela.csv
A face do talude desce a -45 graus (estratos paralelos à face: theta = 135 graus).
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


THETAS = tuple(range(0, 180, 15))


def rodar(a):
    if a.acao == "rotacao":
        k1, k2 = (tuple(map(float, a.casos[0].split(","))) if a.casos else (1e-4, 1e-5))
        casos = [(k1, k2, th) for th in THETAS]
        saida = a.saida.replace("anisotropia", "rotacao")
    else:
        casos = [tuple(map(float, c.split(","))) + (0.,) for c in a.casos] if a.casos else [c + (0.,) for c in CASOS]
        saida = a.saida
    tarefas = []
    for kx, ky, th in casos:
        nome = rotulo(kx, ky, a.chuva) + (f"_th{th:g}" if a.acao == "rotacao" else "")
        out = os.path.join(saida, nome)
        if os.path.exists(os.path.join(out, "summary.txt")) and not a.refazer:
            continue
        os.makedirs(out, exist_ok=True)
        args = [a.exe, "solo=f50", "alpha=45", "hs=10", "hw=5", f"I={a.chuva:g}", f"kx={kx:g}", f"ky={ky:g}",
                f"theta={th:g}", "tfs=" + ",".join(f"{t:g}" for t in TFS), f"out={out}",
                f"gi={0 if a.sem_gi else 1}"]
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


def figura_rotacao(a):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    pasta = a.saida.replace("anisotropia", "rotacao")
    R = [ler(os.path.join(pasta, d)) for d in sorted(os.listdir(pasta))
         if os.path.exists(os.path.join(pasta, d, "summary.txt"))]
    R.sort(key=lambda r: r["theta"])
    iso = [ler(os.path.join(a.saida, d)) for d in os.listdir(a.saida) if "kx1e-5_ky1e-5" in d and
           os.path.exists(os.path.join(a.saida, d, "summary.txt"))] if os.path.isdir(a.saida) else []
    os.makedirs(a.figuras, exist_ok=True)
    fig, axs = plt.subplots(1, 3, figsize=(15, 4.3))
    cores = plt.cm.twilight(np.linspace(0, 1, len(R) + 1))
    for c, r in zip(cores, R):
        axs[0].plot(r["t"], r["fs"], "-", lw=1.2, color=c, label=f"θ = {r['theta']:g}°")
    if iso:
        axs[0].plot(iso[0]["t"], iso[0]["fs"], "k--", lw=1, label="isotrópico (k$_s$)")
    axs[0].axvspan(0, 24, color="#e3effb", lw=0)
    axs[0].axhline(1, color="#999", lw=0.6, ls=":")
    axs[0].set(xlabel="t (h)", ylabel="FS (SRM)", title="FS por redução de resistência", xlim=(0, 48))
    axs[0].legend(frameon=False, fontsize=6.5, ncol=2)
    th = np.array([r["theta"] for r in R])
    axs[1].plot(th, [r["fs_min_srm"] for r in R], "o-", color="#1f3fa8", label="FS mínimo, SRM")
    if any(r["fs_min_gi"] > 0 for r in R):
        axs[1].plot(th, [r["fs_min_gi"] for r in R], "s-", color="#1d7a1d", label="FS mínimo, aumento da gravidade")
    if iso:
        axs[1].axhline(iso[0]["fs_min_srm"], color="#1f3fa8", ls="--", lw=0.8, label="isotrópico, SRM")
        axs[1].axhline(iso[0]["fs_min_gi"], color="#1d7a1d", ls="--", lw=0.8, label="isotrópico, GI")
    axs[1].axvline(135, color="gray", ls=":", lw=0.8)
    axs[1].text(135, axs[1].get_ylim()[1], " k$_x$ paralelo à face", va="top", fontsize=7, color="gray")
    axs[1].set(xlabel="θ (graus): direção de k$_x$ a partir do eixo x", ylabel="FS mínimo", xticks=THETAS,
               title="FS mínimo × orientação do tensor")
    axs[1].legend(frameon=False, fontsize=7)
    infil = [np.sum(r["hist"][:, 4] * r["hist"][:, 1] * 3600 * (r["hist"][:, 2] > 0)) for r in R]
    axs[2].plot(th, infil, "o-", color="#c8553d")
    axs[2].set(xlabel="θ (graus)", ylabel="volume infiltrado em 24 h (m³/m)", xticks=THETAS,
               title="Infiltração durante a chuva")
    if R:
        r = R[0]
        fig.suptitle(f"Solo f50, α = 45°, H$_s$ = 10 m, H$_w$ = 5 m, chuva {r['I']:g} mm/h por 24 h; "
                     f"k$_x$ = {r['kx']:.0e}, k$_y$ = {r['ky']:.0e} m/s girados de θ")
    fig.tight_layout()
    fig.savefig(os.path.join(a.figuras, "rotacao_fs.png"), dpi=150)
    with open(os.path.join(a.figuras, "rotacao_tabela.csv"), "w") as f:
        f.write("theta,kx,ky,fs_ini_srm,fs_min_srm,t_min_srm,fs_48h_srm,fs_min_gi,t_fail,infil_24h_m3m,pa_max_kPa\n")
        for r, v in zip(R, infil):
            f.write(f"{r['theta']:g},{r['kx']:g},{r['ky']:g},{r['fs'][0]:.3f},{r['fs_min_srm']:.3f},"
                    f"{r['t_min_srm']:g},{r['fs'][-1]:.3f},{r['fs_min_gi']:.3f},{r['t_fail']:g},{v:.2f},"
                    f"{r['hist'][:, 7].max():.1f}\n")
    print(open(os.path.join(a.figuras, "rotacao_tabela.csv")).read())


def main():
    aqui = os.path.dirname(os.path.abspath(__file__))
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("acao", choices=["rodar", "figura", "rotacao", "figura-rotacao"])
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
    {"rodar": rodar, "rotacao": rodar, "figura": figura, "figura-rotacao": figura_rotacao}[a.acao](a)


if __name__ == "__main__":
    main()
