"""Contraexemplos às conclusões de Rahardjo et al. (2007): itens A1, A3, A4, A5 e A6 (o A2 está em fase_ar.py).

Modelo de base: o do artigo (percolação desacoplada e p_a = 0; acoplado=0 passivo=1), para que só a variável estudada
mude; o A1 usa o modelo acoplado com ar ativo (o patamar FS ~ 1 só existe com acoplamento). FS por redução de
resistência (SRM). Altura Hs = 10 m e Hw = 5 m salvo indicação.

  A1  dilatância: f50, alpha = 45 graus, I_r = 18, 36, 54, 90, 180 mm/h; psi = 0, 13 e 26 graus
      (SlopeRainfallPV, Mohr-Coulomb iterativo não associativo)
  A3  anisotropia na sigmoide de I_r: f50, alpha = 45 graus, k_1 = 1e-4 e k_2 = 1e-5 m/s, theta = 0, 45, 135 graus;
      I_r = 9, 18, 36, 54, 90, 180, 360 mm/h
  A4  limite da sucção inicial: f50, alpha = 45 graus, Hw = 2,5 / 5 / 7,5 / 10 / 15 m, I_r = 9 e 36 mm/h,
      smax = 75, 150 e 300 kPa
  A5  ângulo: alpha = 20 ... 70 graus a cada 5; (f10, 80 mm/h), (f50, 36 mm/h), (f100, 9 mm/h)
  A6  mesmo volume (864 mm = 36 mm/h x 24 h) em 3, 6, 12, 24 e 72 h, solos f10, f50 e f100

Uso:
  python contraexemplos.py rodar [--itens A1 A3 A4 A5 A6] --jobs 4 --exe <build>/Projects/SlopeRainfall/SlopeRainfall
  python contraexemplos.py figuras               # figuras/contra_*.png e contra_tabela_*.csv
"""
import argparse
import glob
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np

ARTIGO = dict(passivo=1, acoplado=0)
COMPLETO = dict(passivo=0, acoplado=1)


def tfs_para(train, ttotal):
    t = [x for x in (1, 2, 3, 6, 9, 12, 15, 18, 21, 24, 36, 48, 60, 72) if x <= train]
    if train not in t:
        t.append(train)
    for d in (3, 6, 12, 24):
        if train + d <= ttotal:
            t.append(train + d)
    if ttotal not in t:
        t.append(ttotal)
    return sorted(set(t))


def casos(itens):
    out = []
    if "A1" in itens:
        for psi in (0.0, 13.0, 26.0):
            for I in (18.0, 36.0, 54.0, 90.0, 180.0):
                out.append(dict(item="A1", solo="f50", alpha=45.0, hw=5.0, I=I, psi=psi, pv=True, **COMPLETO))
    if "A3" in itens:
        for th in (0.0, 45.0, 135.0):
            for I in (9.0, 18.0, 36.0, 54.0, 90.0, 180.0, 360.0):
                out.append(dict(item="A3", solo="f50", alpha=45.0, hw=5.0, I=I, kx=1e-4, ky=1e-5, theta=th, **ARTIGO))
    if "A4" in itens:
        for smax in (75.0, 150.0, 300.0):
            for hw in (2.5, 5.0, 7.5, 10.0, 15.0):
                for I in (9.0, 36.0):
                    out.append(dict(item="A4", solo="f50", alpha=45.0, hw=hw, I=I, smax=smax, **ARTIGO))
    if "A5" in itens:
        for solo, I in (("f10", 80.0), ("f50", 36.0), ("f100", 9.0)):
            for a in range(20, 71, 5):
                out.append(dict(item="A5", solo=solo, alpha=float(a), hw=5.0, I=I, **ARTIGO))
    if "A6" in itens:
        for solo in ("f10", "f50", "f100"):
            for d in (3.0, 6.0, 12.0, 24.0, 72.0):
                out.append(dict(item="A6", solo=solo, alpha=45.0, hw=5.0, I=864.0 / d, train=d, **ARTIGO))
    for c in out:
        c.setdefault("train", 24.0)
        c["ttotal"] = max(48.0, c["train"] + 24.0)
    return out


def rotulo(c):
    r = f"{c['item']}_{c['solo']}_a{c['alpha']:g}_Hw{c['hw']:g}_I{c['I']:.4g}_tr{c['train']:g}"
    if "psi" in c:
        r += f"_psi{c['psi']:g}"
    if "theta" in c:
        r += f"_th{c['theta']:g}"
    if "smax" in c:
        r += f"_smax{c['smax']:g}"
    return r


def rodar(a):
    exes = [a.exe] + ([a.exe + "PV"] if "A1" in getattr(a, "itens", []) else [])
    for e in exes:
        if not (os.path.isfile(e) and os.access(e, os.X_OK)):
            raise SystemExit(f"executável não encontrado: {e}\ncompile antes (ninja SlopeRainfall"
                             f"{' SlopeRainfallPV' if e.endswith('PV') else ''}) e passe o caminho em --exe")
    tarefas = []
    for c in casos(a.itens):
        out = os.path.join(a.saida, rotulo(c))
        if os.path.exists(os.path.join(out, "summary.txt")) and not a.refazer:
            continue
        os.makedirs(out, exist_ok=True)
        exe = a.exe + ("PV" if c.get("pv") else "")
        args = [exe, f"solo={c['solo']}", f"alpha={c['alpha']:g}", "hs=10", f"hw={c['hw']:g}", f"I={c['I']:.6g}",
                f"train={c['train']:g}", f"ttotal={c['ttotal']:g}", f"passivo={c['passivo']}",
                f"acoplado={c['acoplado']}", "gi=0", "tfs=" + ",".join(f"{t:g}" for t in tfs_para(c["train"], c["ttotal"])),
                f"out={out}"]
        for k in ("psi", "theta", "smax", "kx", "ky"):
            if k in c:
                args.append(f"{k}={c[k]:g}")
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
    s.update(t=fs[:, 0], fs=fs[:, 1], hist=h, pasta=pasta)
    s["infil"] = float(np.sum(h[:, 4] * h[:, 1] * 3600 * (h[:, 2] > 0)))
    s["completo"] = bool(len(fs) and fs[-1, 0] >= s["ttotal"] - 1e-6)
    return s


def figuras(a):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import ref_artigo as ref
    plt.rcParams.update({"font.size": 9, "axes.linewidth": 0.6, "legend.fontsize": 7.5, "axes.titlesize": 9.5})
    R = []
    for c in casos(["A1", "A3", "A4", "A5", "A6"]):
        p = os.path.join(a.saida, rotulo(c))
        if os.path.exists(os.path.join(p, "summary.txt")):
            R.append(dict(c, **ler(p)))
    os.makedirs(a.figuras, exist_ok=True)
    cores = ["#1b6ca8", "#e9a23b", "#c8553d", "#2a9d8f", "#6c3483"]

    def base(serie, solo, I=None):
        """resultados das séries do artigo (resultados/<serie>/<caso>)"""
        out = []
        for f in glob.glob(os.path.join(a.resultados, serie, f"{solo}_a45_H10_Hw5_I*", "summary.txt")):
            d = ler(os.path.dirname(f))
            out.append((d["I"], d["fs_min_srm"]))
        return sorted(out)

    def tabela(nome, rs, chaves):
        with open(os.path.join(a.figuras, f"contra_tabela_{nome}.csv"), "w") as f:
            f.write(",".join(chaves) + ",fs_ini,fs_min,t_min,t_fail,infil_m3m,completo,minutos\n")
            for r in rs:
                f.write(",".join(f"{r.get(k, '')}" for k in chaves) +
                        f",{r['fs'][0]:.3f},{r['fs_min_srm']:.3f},{r['t_min_srm']:g},{r['t_fail']:g},{r['infil']:.2f},"
                        f"{int(r['completo'])},{r['minutes']:.1f}\n")

    # A1
    rs = [r for r in R if r["item"] == "A1"]
    if rs:
        fig, ax = plt.subplots(figsize=(6.5, 4.3))
        Ir = np.logspace(0.5, 2.6, 100)
        ax.semilogx(Ir, ref.fs_min_eq6("f50", Ir), "k-", lw=1, label="Eq. 6 do artigo")
        for k, psi in enumerate((0.0, 13.0, 26.0)):
            p = sorted((r["I"], r["fs_min_srm"]) for r in rs if r["psi"] == psi)
            if p:
                ax.semilogx(*zip(*p), "o-", color=cores[k], label=f"acoplado, ar ativo, ψ = {psi:g}°")
        for s, m, rot in (("B", "s", "acoplado, ψ = φ' (modelo analítico)"), ("B_passivo_desacoplado", "^",
                                                                              "desacoplado, p$_a$ = 0 (artigo)")):
            p = [x for x in base(s, "f50") if x[0] >= 18]
            if p:
                ax.semilogx(*zip(*p), m + "--", color="gray", mfc="w", label=rot)
        ax.axhline(1, color="#999", lw=0.6, ls=":")
        ax.set(xlabel="I$_r$ (mm/h)", ylabel="FS mínimo", title="A1: dilatância (f$_{50}$, α = 45°)")
        ax.legend(frameon=False, fontsize=7)
        fig.tight_layout()
        fig.savefig(os.path.join(a.figuras, "contra_A1_dilatancia.png"), dpi=150)
        plt.close(fig)
        tabela("A1", sorted(rs, key=lambda r: (r["psi"], r["I"])), ["psi", "I"])

    # A3
    rs = [r for r in R if r["item"] == "A3"]
    if rs:
        fig, ax = plt.subplots(figsize=(6.5, 4.3))
        Ir = np.logspace(0.5, 2.7, 100)
        ax.semilogx(Ir, ref.fs_min_eq6("f50", Ir), "k-", lw=1, label="Eq. 6 (f$_{50}$, isotrópico)")
        p = base("B_passivo_desacoplado", "f50")
        if p:
            ax.semilogx(*zip(*p), "^--", color="gray", mfc="w", label="isotrópico k$_s$ = 1e-5 (NeoPZ)")
        for k, th in enumerate((0.0, 45.0, 135.0)):
            q = sorted((r["I"], r["fs_min_srm"]) for r in rs if r["theta"] == th)
            if q:
                ax.semilogx(*zip(*q), "o-", color=cores[k], label=f"k$_1$ = 1e-4, k$_2$ = 1e-5, θ = {th:g}°")
        ax.axhline(1, color="#999", lw=0.6, ls=":")
        ax.set(xlabel="I$_r$ (mm/h)", ylabel="FS mínimo", title="A3: anisotropia e a intensidade-limite")
        ax.legend(frameon=False, fontsize=7)
        fig.tight_layout()
        fig.savefig(os.path.join(a.figuras, "contra_A3_anisotropia.png"), dpi=150)
        plt.close(fig)
        tabela("A3", sorted(rs, key=lambda r: (r["theta"], r["I"])), ["theta", "I"])

    # A4
    rs = [r for r in R if r["item"] == "A4"]
    if rs:
        fig, axs = plt.subplots(1, 2, figsize=(12, 4.2))
        for k, smax in enumerate((75.0, 150.0, 300.0)):
            for I, ls in ((9.0, "-"), (36.0, "--")):
                q = sorted((r["hw"], r["fs"][0], r["fs_min_srm"]) for r in rs if r["smax"] == smax and r["I"] == I)
                if not q:
                    continue
                hw, ini, mn = map(np.array, zip(*q))
                if I == 9.0:
                    axs[0].plot(hw, ini, "o-", color=cores[k], label=f"FS inicial, s$_{{max}}$ = {smax:g} kPa")
                axs[0].plot(hw, mn, "s" + ls, color=cores[k], mfc="w" if I == 36 else cores[k],
                            label=f"FS mínimo, {I:g} mm/h, s$_{{max}}$ = {smax:g}")
                axs[1].plot(hw, ini - mn, "s" + ls, color=cores[k], label=f"{I:g} mm/h, s$_{{max}}$ = {smax:g} kPa")
        axs[0].plot(ref.FIG8_HW, ref.FIG8["ini"], "k-", lw=0.8, label="FS inicial, Fig. 8")
        axs[0].set(xlabel="H$_w$ (m)", ylabel="FS", title="A4: FS inicial e mínimo")
        axs[1].set(xlabel="H$_w$ (m)", ylabel="FS inicial − FS mínimo", title="Queda do FS pela chuva")
        axs[0].legend(frameon=False, fontsize=6.5, ncol=2)
        axs[1].legend(frameon=False, fontsize=7)
        fig.tight_layout()
        fig.savefig(os.path.join(a.figuras, "contra_A4_succao.png"), dpi=150)
        plt.close(fig)
        tabela("A4", sorted(rs, key=lambda r: (r["smax"], r["I"], r["hw"])), ["smax", "I", "hw"])

    # A5
    rs = [r for r in R if r["item"] == "A5"]
    if rs:
        fig, axs = plt.subplots(1, 2, figsize=(12, 4.2))
        al = np.linspace(20, 70, 50)
        axs[0].plot(al, ref.FS_INI[0] + ref.FS_INI[1] * al, "k-", lw=0.8, label="FS inicial, Tabela 7")
        for k, (solo, I) in enumerate((("f10", 80.0), ("f50", 36.0), ("f100", 9.0))):
            q = sorted((r["alpha"], r["fs"][0], r["fs_min_srm"]) for r in rs if r["solo"] == solo)
            if not q:
                continue
            a_, ini, mn = map(np.array, zip(*q))
            if k == 0:
                axs[0].plot(a_, ini, "ko", mfc="w", label="FS inicial, NeoPZ")
            axs[0].plot(a_, mn, "s-", color=cores[k], label=f"{solo}, {I:g} mm/h, NeoPZ")
            if (solo, I) in ref.FIG9:
                axs[0].plot(ref.FIG9_ALPHA, ref.FIG9[(solo, I)], "x--", color=cores[k], lw=0.8, label=f"{solo}, {I:g} mm/h, Fig. 9")
            if (solo, I) in ref.TABELA7:
                E, F = ref.TABELA7[(solo, I)]
                axs[0].plot(al, E + F * al, ":", color=cores[k], lw=0.8)
            # inclinação local dFS/dalpha
            axs[1].plot(0.5 * (a_[1:] + a_[:-1]), np.diff(mn) / np.diff(a_), "o-", color=cores[k], label=f"{solo}, {I:g} mm/h")
        axs[1].plot(0.5 * (a_[1:] + a_[:-1]), np.diff(ini) / np.diff(a_), "ko--", mfc="w", label="FS inicial")
        axs[1].axhline(-0.025, color="gray", ls=":", label="F = −0,025 (Tabela 7)")
        axs[0].axhline(1, color="#999", lw=0.6, ls=":")
        axs[0].set(xlabel="α (graus)", ylabel="FS", title="A5: FS mínimo × α (pontilhado: Tabela 7)")
        axs[1].set(xlabel="α (graus)", ylabel="dFS/dα (1/grau)", title="Inclinação local")
        axs[0].legend(frameon=False, fontsize=6.5)
        axs[1].legend(frameon=False, fontsize=7)
        fig.tight_layout()
        fig.savefig(os.path.join(a.figuras, "contra_A5_angulo.png"), dpi=150)
        plt.close(fig)
        tabela("A5", sorted(rs, key=lambda r: (r["solo"], r["alpha"])), ["solo", "alpha", "I"])

    # A6
    rs = [r for r in R if r["item"] == "A6"]
    if rs:
        fig, axs = plt.subplots(1, 2, figsize=(12, 4.2))
        for k, solo in enumerate(("f10", "f50", "f100")):
            q = sorted((r["train"], r["I"], r["fs_min_srm"], r["infil"]) for r in rs if r["solo"] == solo)
            if not q:
                continue
            d, I, mn, inf = map(np.array, zip(*q))
            axs[0].semilogx(d, mn, "o-", color=cores[k], label=solo)
            axs[1].semilogx(d, inf, "o-", color=cores[k], label=solo)
        axs[0].axhline(1, color="#999", lw=0.6, ls=":")
        axs[0].set(xlabel="duração da chuva (h)  [I$_r$ = 864 mm / duração]", ylabel="FS mínimo",
                   title="A6: mesmo volume (864 mm), durações diferentes")
        axs[1].set(xlabel="duração da chuva (h)", ylabel="volume infiltrado (m³/m)", title="Infiltração")
        axs[0].legend(frameon=False)
        fig.tight_layout()
        fig.savefig(os.path.join(a.figuras, "contra_A6_volume.png"), dpi=150)
        plt.close(fig)
        tabela("A6", sorted(rs, key=lambda r: (r["solo"], r["train"])), ["solo", "train", "I"])
    print(f"{len(R)} casos lidos")


def main():
    aqui = os.path.dirname(os.path.abspath(__file__))
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("acao", choices=["rodar", "figuras"])
    p.add_argument("--itens", nargs="*", default=["A1", "A3", "A4", "A5", "A6"])
    p.add_argument("--exe", default=os.path.join(aqui, "..", "..", "..", "..", "build", "Projects", "SlopeRainfall",
                                                  "SlopeRainfall"))
    p.add_argument("--saida", default=os.path.join("resultados", "contraexemplos"))
    p.add_argument("--resultados", default="resultados")
    p.add_argument("--figuras", default="figuras")
    p.add_argument("--jobs", type=int, default=4)
    p.add_argument("--refazer", action="store_true")
    a = p.parse_args()
    rodar(a) if a.acao == "rodar" else figuras(a)


if __name__ == "__main__":
    main()
