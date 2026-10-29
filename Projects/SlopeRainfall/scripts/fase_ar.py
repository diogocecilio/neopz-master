"""Estudo A2: a fase de ar contradiz "taludes pouco permeáveis são seguros em chuvas curtas" (Rahardjo et al. 2007)?

Solo f100 (k_s = 1e-6 m/s), chuva de 24 h, 48 h simuladas, FS por redução de resistência. Três modelos:
  ativo    acoplado, p_a calculada (p_a = 0 só na superfície)
  passivo  acoplado, p_a = 0 em todo o domínio
  artigo   desacoplado e p_a = 0 (hipóteses de Rahardjo et al., SEEP/W)
Grupos de casos:
  G1  alpha = 63,4 graus, Hs = 10 m, Hw = 5 m; I_r = 9, 18, 36, 80 mm/h; os três modelos
  G2  alpha = 45 graus, Hs = 20 m, Hw = 5 m; I_r = 9, 36, 80 mm/h; ativo e artigo
  G3  G1 com I_r = 36 mm/h, ar ativo, k_ra,min = 1e-6, 1e-5, 1e-4, 1e-3, 1e-2 (saída do ar)
  G4  G1 com 36 mm/h e ar ativo na malha refinada: elementos de h_sup = 0,75 (a malha padrão), 0,5 e 0,35 m até 2 m
      de profundidade e h_max = 0,3 Hs h_sup/0,75 (malha_talude.py, gerada em meshes/ se faltar); tolerância do SRM
      0,002

Uso:
  python fase_ar.py rodar --jobs 4 --exe <build>/Projects/SlopeRainfall/SlopeRainfall
  python fase_ar.py figuras                      # figuras/fase_ar_*.png e fase_ar_tabela.csv
"""
import argparse
import os
import re
import subprocess
import sys

# um caso por núcleo: o MKL (Pardiso) e o OpenMP/TBB não devem abrir threads em todos os núcleos em cada processo
ENV1 = dict(os.environ, MKL_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_DYNAMIC="FALSE", TBB_NUM_THREADS="1")
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np

TFS = (3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 36, 48)
MODOS = {"ativo": dict(passivo=0, acoplado=1), "passivo": dict(passivo=1, acoplado=1),
         "artigo": dict(passivo=1, acoplado=0)}
ROT_MODO = {"ativo": "acoplado, ar ativo", "passivo": "acoplado, ar passivo",
            "artigo": "desacoplado, ar passivo (artigo)"}
COR_MODO = {"ativo": "#1b6ca8", "passivo": "#e9a23b", "artigo": "#c8553d"}


def casos():
    out = []
    for I in (9.0, 18.0, 36.0, 80.0):
        for m in MODOS:
            out.append(dict(grupo="G1", alpha=63.4, hs=10.0, I=I, modo=m, kramin=1e-4))
    for I in (9.0, 36.0, 80.0):
        for m in ("ativo", "artigo"):
            out.append(dict(grupo="G2", alpha=45.0, hs=20.0, I=I, modo=m, kramin=1e-4))
    for kr in (1e-6, 1e-5, 1e-3, 1e-2):
        out.append(dict(grupo="G3", alpha=63.4, hs=10.0, I=36.0, modo="ativo", kramin=kr))
    for h in MALHAS:
        c = dict(grupo="G4", alpha=63.4, hs=10.0, I=36.0, modo="ativo", kramin=1e-4, malha=h, tolfs=0.002)
        if h < 0.4:
            c["threads"] = 2        # a malha mais fina é o caminho crítico do lote
        out.append(c)
    return out


MALHAS = (0.75, 0.5, 0.35)     # h_sup (m) do G4; 0,75 é a malha padrão (h_max = 0,3 Hs)


def malha(c, pasta_malhas):
    """arquivo da malha do G4 (None: a malha padrão do programa); gera com malha_talude.py se faltar"""
    h = c.get("malha", 0.75)
    if h == 0.75:
        return None
    arq = os.path.join(pasta_malhas, f"talude_H{c['hs']:g}_a{c['alpha']:g}_hsup{h:g}.msh")
    if not os.path.exists(arq):
        aqui = os.path.dirname(os.path.abspath(__file__))
        subprocess.run([sys.executable, os.path.join(aqui, "malha_talude.py"), f"{c['hs']:g}", f"{c['alpha']:g}", f"{h:g}",
                        f"{0.3 * h / 0.75:g}", arq], check=True)
    return arq


def rotulo(c):
    r = f"f100_a{c['alpha']:g}_H{c['hs']:g}_I{c['I']:g}_{c['modo']}_kra{c['kramin']:.0e}".replace("e-0", "e-")
    return r + (f"_m{c['malha']:g}" if "malha" in c else "")


def rodar(a):
    exes = [a.exe] + ([a.exe + "PV"] if "A1" in getattr(a, "itens", []) else [])
    for e in exes:
        if not (os.path.isfile(e) and os.access(e, os.X_OK)):
            raise SystemExit(f"executável não encontrado: {e}\ncompile antes (ninja SlopeRainfall"
                             f"{' SlopeRainfallPV' if e.endswith('PV') else ''}) e passe o caminho em --exe")
    tarefas = []
    for c in casos():
        out = os.path.join(a.saida, rotulo(c))
        rodando = subprocess.run(["pgrep", "-f", f"out={re.escape(out)}( |$)"], capture_output=True).returncode == 0
        if (os.path.exists(os.path.join(out, "summary.txt")) and not a.refazer) or rodando:
            continue
        os.makedirs(out, exist_ok=True)
        m = MODOS[c["modo"]]
        args = [a.exe, "solo=f100", f"alpha={c['alpha']:g}", f"hs={c['hs']:g}", "hw=5", f"I={c['I']:g}",
                f"passivo={m['passivo']}", f"acoplado={m['acoplado']}", f"kramin={c['kramin']:g}", "gi=0",
                "tfs=" + ",".join(f"{t:g}" for t in TFS), f"out={out}"]
        arq = malha(c, a.malhas)
        if arq:
            args.append(f"mesh={arq}")
        if "tolfs" in c:
            args.append(f"tolfs={c['tolfs']:g}")
        env = ENV1
        if "threads" in c:
            n = str(c["threads"])
            args.append(f"threads={n}")
            env = dict(ENV1, MKL_NUM_THREADS=n, OMP_NUM_THREADS=n, TBB_NUM_THREADS=n)
        tarefas.append(((c["hs"] / c.get("malha", 0.75)) ** 2, out, args, env))
    tarefas.sort(key=lambda t: -t[0])            # os casos mais pesados primeiro
    print(f"{len(tarefas)} casos")

    def roda(t):
        _, out, args, env = t
        with open(os.path.join(out, "log.txt"), "w") as log:
            return out, subprocess.run(args, stdout=log, stderr=subprocess.STDOUT, env=env).returncode

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
    s.update(t=fs[:, 0], fs=fs[:, 1], hist=h, pasta=pasta)
    s["infil"] = float(np.sum(h[:, 4] * h[:, 1] * 3600 * (h[:, 2] > 0)))
    s["pa_max"] = float(h[:, 7].max())
    return s


def figuras(a):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.tri as mtri
    plt.rcParams.update({"font.size": 9, "axes.linewidth": 0.6, "legend.fontsize": 7.5, "axes.titlesize": 9.5})
    R = {}
    for c in casos():
        p = os.path.join(a.saida, rotulo(c))
        if os.path.exists(os.path.join(p, "summary.txt")):
            R[rotulo(c)] = dict(c, **ler(p))
    os.makedirs(a.figuras, exist_ok=True)
    get = lambda **k: [r for r in R.values() if all(r[x] == v for x, v in k.items())]

    # G1: FS(t) por intensidade
    fig, axs = plt.subplots(1, 4, figsize=(15, 3.8), sharey=True)
    for ax, I in zip(axs, (9.0, 18.0, 36.0, 80.0)):
        ax.axvspan(0, 24, color="#e3effb", lw=0)
        for m in MODOS:
            for r in get(grupo="G1", I=I, modo=m):
                ax.plot(r["t"], r["fs"], "-o", ms=2.5, color=COR_MODO[m], label=ROT_MODO[m])
                if r["t_fail"] >= 0:
                    ax.axvline(r["t_fail"], color=COR_MODO[m], ls=":", lw=1)
        ax.axhline(1, color="#999", lw=0.6, ls=":")
        ax.set(xlabel="t (h)", title=f"I$_r$ = {I:g} mm/h", xlim=(0, 48))
    axs[0].set_ylabel("FS (SRM)")
    axs[0].legend(frameon=False)
    fig.suptitle("G1: solo f$_{100,-6}$, α = 63,4°, H$_s$ = 10 m, H$_w$ = 5 m")
    fig.tight_layout()
    fig.savefig(os.path.join(a.figuras, "fase_ar_G1_fs_tempo.png"), dpi=150)
    plt.close(fig)

    # G1 e G2: FS mínimo x I_r, p_a máxima e infiltração
    fig, axs = plt.subplots(1, 3, figsize=(15, 4))
    for g, mk, tit in (("G1", "o", "α = 63,4°, H$_s$ = 10 m"), ("G2", "s", "α = 45°, H$_s$ = 20 m")):
        for m in MODOS:
            rs = sorted(get(grupo=g, modo=m), key=lambda r: r["I"])
            if not rs:
                continue
            I = [r["I"] for r in rs]
            ls = "-" if g == "G1" else "--"
            axs[0].semilogx(I, [r["fs_min_srm"] for r in rs], mk + ls, color=COR_MODO[m], label=f"{tit}: {ROT_MODO[m]}")
            axs[1].semilogx(I, [r["pa_max"] for r in rs], mk + ls, color=COR_MODO[m])
            axs[2].semilogx(I, [r["infil"] for r in rs], mk + ls, color=COR_MODO[m])
    import ref_artigo as ref
    axs[0].plot([3.6, 9.0, 80.0], [ref.FIG9[("f100", 3.6)][3], ref.FIG9[("f100", 9.0)][3], ref.FIG9[("f100", 80.0)][3]],
                "kx", ms=8, label="Fig. 9 do artigo (α = 63,4°)")
    axs[0].plot([9.0, 36.0, 80.0], [1.37, 0.97, 0.88], "k+", ms=9, label="Fig. 10 do artigo (f$_{50}$, H$_s$ = 20 m)")
    axs[0].axhline(1, color="#999", lw=0.6, ls=":")
    axs[0].set(xlabel="I$_r$ (mm/h)", ylabel="FS mínimo", title="FS mínimo (24 h de chuva)")
    axs[0].legend(frameon=False, fontsize=6.5)
    axs[1].set(xlabel="I$_r$ (mm/h)", ylabel="p$_a$ máxima (kPa)", title="Pressão máxima do ar")
    axs[2].set(xlabel="I$_r$ (mm/h)", ylabel="volume infiltrado em 24 h (m³/m)", title="Infiltração")
    fig.tight_layout()
    fig.savefig(os.path.join(a.figuras, "fase_ar_fsmin.png"), dpi=150)
    plt.close(fig)

    # G3: sensibilidade a k_ra,min
    rs = sorted(get(grupo="G3") + get(grupo="G1", I=36.0, modo="ativo"), key=lambda r: r["kramin"])
    if rs:
        fig, axs = plt.subplots(1, 3, figsize=(15, 4))
        kr = [r["kramin"] for r in rs]
        axs[0].semilogx(kr, [r["fs_min_srm"] for r in rs], "o-", color=COR_MODO["ativo"], label="ar ativo")
        for m in ("passivo", "artigo"):
            for r in get(grupo="G1", I=36.0, modo=m):
                axs[0].axhline(r["fs_min_srm"], color=COR_MODO[m], ls="--", label=ROT_MODO[m])
        axs[0].axhline(1, color="#999", lw=0.6, ls=":")
        axs[0].set(xlabel="k$_{ra,min}$", ylabel="FS mínimo", title="FS mínimo × saída do ar")
        axs[0].legend(frameon=False)
        axs[1].semilogx(kr, [r["pa_max"] for r in rs], "o-", color=COR_MODO["ativo"])
        axs[1].set(xlabel="k$_{ra,min}$", ylabel="p$_a$ máxima (kPa)", title="Pressão máxima do ar")
        cores = plt.cm.viridis(np.linspace(0, 0.9, len(rs)))
        for c, r in zip(cores, rs):
            axs[2].plot(r["t"], r["fs"], "-", color=c, label=f"k$_{{ra,min}}$ = {r['kramin']:.0e}")
        axs[2].axvspan(0, 24, color="#e3effb", lw=0)
        axs[2].axhline(1, color="#999", lw=0.6, ls=":")
        axs[2].set(xlabel="t (h)", ylabel="FS", title="FS(t)", xlim=(0, 48))
        axs[2].legend(frameon=False, fontsize=7)
        fig.suptitle("G3: f$_{100,-6}$, α = 63,4°, H$_s$ = 10 m, 36 mm/h, ar ativo")
        fig.tight_layout()
        fig.savefig(os.path.join(a.figuras, "fase_ar_G3_kramin.png"), dpi=150)
        plt.close(fig)

    # campos em t = 24 h para I = 36 mm/h (G1): p_w nos três modelos e p_a no ativo
    rs = [get(grupo="G1", I=36.0, modo=m) for m in MODOS]
    if all(rs):
        fig, axs = plt.subplots(2, 2, figsize=(12, 6))
        paineis = [(rs[0][0], "snaps_pw.csv", "p$_w$, " + ROT_MODO["ativo"]),
                   (rs[1][0], "snaps_pw.csv", "p$_w$, " + ROT_MODO["passivo"]),
                   (rs[2][0], "snaps_pw.csv", "p$_w$, " + ROT_MODO["artigo"]),
                   (rs[0][0], "snaps_pa.csv", "p$_a$, " + ROT_MODO["ativo"])]
        for ax, (r, arq, tit) in zip(axs.ravel(), paineis):
            p = r["pasta"]
            v = np.loadtxt(os.path.join(p, "vertices.csv"), delimiter=",", skiprows=1)
            t = np.loadtxt(os.path.join(p, "triangles.csv"), delimiter=",", skiprows=1, dtype=int)
            o = np.loadtxt(os.path.join(p, "outline.csv"), delimiter=",", skiprows=1)
            S = np.loadtxt(os.path.join(p, arq), delimiter=",")
            k = int(np.argmin(np.abs(S[:, 0] - 24.)))
            tri = mtri.Triangulation(v[:, 0], v[:, 1], t)
            if arq == "snaps_pw.csv":
                cf = ax.tricontourf(tri, np.clip(S[k, 1:], -80, 200), levels=np.arange(-80, 201, 10), cmap="RdYlBu_r")
                if S[k, 1:].min() < 0 < S[k, 1:].max():
                    ax.tricontour(tri, S[k, 1:], levels=[0.], colors="magenta", linewidths=1.5)
            else:
                cf = ax.tricontourf(tri, S[k, 1:], levels=12, cmap="viridis")
            fig.colorbar(cf, ax=ax, shrink=0.8, label="kPa")
            ax.plot(np.r_[o[:, 0], o[0, 0]], np.r_[o[:, 1], o[0, 1]], "k-", lw=0.7)
            ax.set_aspect("equal")
            ax.set_title(f"{tit} (t = 24 h)", fontsize=8.5)
            ax.set_xticks([])
            ax.set_yticks([])
        fig.tight_layout()
        fig.savefig(os.path.join(a.figuras, "fase_ar_campos.png"), dpi=150)
        plt.close(fig)

    # G4: refinamento da malha
    rs = get(grupo="G4")
    if rs:
        fig, axs = plt.subplots(1, 2, figsize=(12, 4.2))
        for m in MODOS:
            for I, ls in ((36.0, "-"), (80.0, "--")):
                q = sorted((r["malha"], r["fs_min_srm"]) for r in rs if r["modo"] == m and r["I"] == I)
                if q:
                    axs[0].plot(*zip(*q), "o" + ls, color=COR_MODO[m], mfc="w" if I == 80 else COR_MODO[m],
                                label=f"{ROT_MODO[m]}, {I:g} mm/h")
        axs[0].invert_xaxis()
        axs[0].axhline(1, color="#999", lw=0.6, ls=":")
        axs[0].set(xlabel="h$_{sup}$ (m), tamanho dos elementos até 2 m de profundidade", ylabel="FS mínimo",
                   title="G4: FS mínimo × malha (tolerância do SRM 0,002)")
        axs[0].legend(frameon=False, fontsize=7)
        cores = plt.cm.viridis(np.linspace(0, 0.85, len(MALHAS)))
        for cor, h in zip(cores, MALHAS):
            for r in rs:
                if r["malha"] == h and r["I"] == 36.0 and r["modo"] == "ativo":
                    axs[1].plot(r["t"], r["fs"], "-o", ms=2.5, color=cor, label=f"h$_{{sup}}$ = {h:g} m")
        axs[1].axvspan(0, 24, color="#e3effb", lw=0)
        axs[1].axhline(1, color="#999", lw=0.6, ls=":")
        axs[1].set(xlabel="t (h)", ylabel="FS", title="FS(t), ar ativo, 36 mm/h", xlim=(0, 48))
        axs[1].legend(frameon=False, fontsize=7)
        fig.suptitle("G4: f$_{100,-6}$, α = 63,4°, H$_s$ = 10 m")
        fig.tight_layout()
        fig.savefig(os.path.join(a.figuras, "fase_ar_G4_malha.png"), dpi=150)
        plt.close(fig)

    with open(os.path.join(a.figuras, "fase_ar_tabela.csv"), "w") as f:
        f.write("grupo,alpha,hs,I,modo,kramin,malha,elementos,fs_ini,fs_min,t_min,t_fail,pa_max_kPa,infil_24h_m3m,"
                "minutos\n")
        for r in sorted(R.values(), key=lambda r: (r["grupo"], r["hs"], r["I"], r["modo"], r["kramin"],
                                                   -r.get("malha", 0.75))):
            f.write(f"{r['grupo']},{r['alpha']:g},{r['hs']:g},{r['I']:g},{r['modo']},{r['kramin']:g},"
                    f"{r.get('malha', 0.75):g},{r['elements']:g},{r['fs'][0]:.3f},"
                    f"{r['fs_min_srm']:.3f},{r['t_min_srm']:g},{r['t_fail']:g},{r['pa_max']:.1f},{r['infil']:.2f},"
                    f"{r['minutes']:.1f}\n")
    print(open(os.path.join(a.figuras, "fase_ar_tabela.csv")).read())


def main():
    aqui = os.path.dirname(os.path.abspath(__file__))
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("acao", choices=["rodar", "figuras"])
    p.add_argument("--exe", default=os.path.join(aqui, "..", "..", "..", "..", "build", "Projects", "SlopeRainfall",
                                                  "SlopeRainfall"))
    p.add_argument("--saida", default=os.path.join("resultados", "fase_ar"))
    p.add_argument("--figuras", default="figuras")
    p.add_argument("--malhas", default=os.path.join(aqui, "..", "meshes"), help="pasta das malhas do G4")
    p.add_argument("--jobs", type=int, default=4)
    p.add_argument("--refazer", action="store_true")
    a = p.parse_args()
    rodar(a) if a.acao == "rodar" else figuras(a)


if __name__ == "__main__":
    main()
