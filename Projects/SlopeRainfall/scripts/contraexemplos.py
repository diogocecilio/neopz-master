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
  A6  mesmo volume (864 mm = 36 mm/h x 24 h) em 3, 6, 12, 24 e 72 h, solos f10, f50 e f100; durações acrescentadas em
      torno de t* = V/k_s: 0,75, 1,5 e 2,4 h (f10, passo de chuva = duração/12), 18, 36 e 48 h (f50), 144, 240 e 480 h
      (f100); f10 com 3 h e passo de 0,0625 h (verificação do passo de tempo)
Itens da discussão dos resultados (critério de ruptura FS <= 1; tolfs = 0,002 nos casos perto de FS = 1):
  A1c intensidade crítica: f50, alpha = 45 graus; modelo completo com psi = 0, 13 e 26 graus e I_r = 18, 22, 26, 30,
      33, 36 e 54 mm/h; modelo do artigo com I_r = 36, 40, 44, 48, 52 e 54 mm/h; FS de hora em hora até 30 h
  A3c anisotropia no modelo completo (acoplado, ar ativo, Mohr-Coulomb associativo): isotrópico (k_s = 1e-5) e
      k_1 = 1e-4, k_2 = 1e-5 m/s com theta = 0, 45 e 135 graus; I_r = 9, 18, 36, 54, 90 e 180 mm/h
  A4c o A4 com chi = S_e (para Fredlund-Xing, S_e = S_w) em vez de chi = 1 na sucção

Uso:
  python contraexemplos.py rodar [--itens A1 A1c A3 A3c A4 A4c A5 A6] --jobs 4 --exe <build>/Projects/SlopeRainfall/SlopeRainfall
  python contraexemplos.py figuras               # figuras/contra_*.png e contra_tabela_*.csv
"""
import argparse
import glob
import os
import subprocess

# um caso por núcleo: o MKL (Pardiso) e o OpenMP/TBB não devem abrir threads em todos os núcleos em cada processo
ENV1 = dict(os.environ, MKL_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_DYNAMIC="FALSE", TBB_NUM_THREADS="1")
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np

ARTIGO = dict(passivo=1, acoplado=0)
COMPLETO = dict(passivo=0, acoplado=1)


ITENS = ["A1", "A1c", "A3", "A3c", "A4", "A4c", "A5", "A6"]
TOL_FINA = 0.002         # tolerância do SRM nos casos perto de FS = 1 (padrão do programa: 0,01)


def tfs_para(train, ttotal):
    """instantes do FS; os de depois da chuva caem na grade de 1 h que começa no fim da chuva"""
    t = [x for x in (1, 2, 3, 6, 9, 12, 15, 18, 21, 24, 36, 48, 60, 72) if x <= train]
    t += [x for x in range(96, int(train), 24)]      # chuvas longas (f100 com 144 a 480 h)
    if train not in t:
        t.append(train)
    for d in ((1, 2, 3, 6, 12, 24) if train < 3 else (3, 6, 12, 24)):
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
        novas = {"f10": (0.75, 1.5, 2.4), "f50": (18.0, 36.0, 48.0), "f100": (144.0, 240.0, 480.0)}
        for solo in ("f10", "f50", "f100"):
            for d in sorted((3.0, 6.0, 12.0, 24.0, 72.0) + novas[solo]):
                c = dict(item="A6", solo=solo, alpha=45.0, hw=5.0, I=864.0 / d, train=d, **ARTIGO)
                if d < 3:
                    c["dtchuva"] = d / 12          # 12 passos na chuva, como na de 3 h
                out.append(c)
        out.append(dict(item="A6", solo="f10", alpha=45.0, hw=5.0, I=288.0, train=3.0, dtchuva=0.0625,
                        verif_dt=True, **ARTIGO))
    if "A1c" in itens:
        for psi in (0.0, 13.0, 26.0):
            for I in (18.0, 22.0, 26.0, 30.0, 33.0, 36.0, 54.0):
                out.append(dict(item="A1c", solo="f50", alpha=45.0, hw=5.0, I=I, psi=psi, pv=True, tolfs=TOL_FINA,
                                horario=True, **COMPLETO))
        for I in (36.0, 40.0, 44.0, 48.0, 52.0, 54.0):
            out.append(dict(item="A1c", solo="f50", alpha=45.0, hw=5.0, I=I, tolfs=TOL_FINA, horario=True, **ARTIGO))
    if "A3c" in itens:
        for I in (9.0, 18.0, 36.0, 54.0, 90.0, 180.0):
            out.append(dict(item="A3c", solo="f50", alpha=45.0, hw=5.0, I=I, tolfs=TOL_FINA, **COMPLETO))
            for th in (0.0, 45.0, 135.0):
                out.append(dict(item="A3c", solo="f50", alpha=45.0, hw=5.0, I=I, kx=1e-4, ky=1e-5, theta=th,
                                tolfs=TOL_FINA, **COMPLETO))
    if "A4c" in itens:
        for smax in (75.0, 150.0, 300.0):
            for hw in (2.5, 5.0, 7.5, 10.0, 15.0):
                for I in (9.0, 36.0):
                    out.append(dict(item="A4c", solo="f50", alpha=45.0, hw=hw, I=I, smax=smax, chi="se", **ARTIGO))
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
    if "chi" in c:
        r += f"_chi{c['chi']}"
    if "dtchuva" in c:
        r += f"_dt{c['dtchuva']:g}"
    return r


def tfs_caso(c):
    if c.get("horario"):     # FS de hora em hora até 30 h: instante em que o FS cruza 1
        return list(range(1, 31)) + [33, 36, 42, 48]
    return tfs_para(c["train"], c["ttotal"])


def rodar(a):
    exes = [a.exe] + ([a.exe + "PV"] if any(c.get("pv") for c in casos(a.itens)) else [])
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
                f"acoplado={c['acoplado']}", "gi=0", "tfs=" + ",".join(f"{t:g}" for t in tfs_caso(c)), f"out={out}"]
        for k in ("psi", "theta", "smax", "kx", "ky", "tolfs", "dtchuva"):
            if k in c:
                args.append(f"{k}={c[k]:g}")
        if "chi" in c:
            args.append(f"chi={c['chi']}")
        # os mais longos primeiro: chuvas longas, FS de hora em hora, modelo acoplado
        custo = c["ttotal"] * (3 if c.get("horario") else 1) * (2 if c["acoplado"] else 1)
        tarefas.append((custo, out, args))
    tarefas.sort(key=lambda t: -t[0])
    print(f"{len(tarefas)} casos")

    def roda(t):
        _, out, args = t
        with open(os.path.join(out, "log.txt"), "w") as log:
            return out, subprocess.run(args, stdout=log, stderr=subprocess.STDOUT, env=ENV1).returncode

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
    for c in casos(ITENS):
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
        fig, axs = plt.subplots(1, 3, figsize=(16, 4.2))
        for k, solo in enumerate(("f10", "f50", "f100")):
            q = sorted((r["train"], r["I"], r["fs_min_srm"], r["infil"]) for r in rs
                       if r["solo"] == solo and not r.get("verif_dt"))
            if not q:
                continue
            d, I, mn, inf = map(np.array, zip(*q))
            axs[0].semilogx(d, mn, "o-", color=cores[k], label=f"{solo} (t* = V/k$_s$ = {864 / KS_MMH[solo]:g} h)")
            axs[0].axvline(864 / KS_MMH[solo], color=cores[k], lw=0.6, ls=":")
            axs[1].semilogx(d, inf, "o-", color=cores[k], label=solo)
            axs[2].semilogx(I / KS_MMH[solo], mn, "o-", color=cores[k], label=solo)
        for r in rs:
            if r.get("verif_dt"):
                axs[0].semilogx(r["train"], r["fs_min_srm"], "*", ms=10, color=cores[0], mfc="w",
                                label=f"f10, 3 h, passo de {r['dtchuva']:g} h na chuva")
        axs[0].axhline(1, color="#999", lw=0.6, ls=":")
        axs[2].axvline(1, color="#999", lw=0.6, ls=":")
        axs[0].set(xlabel="duração da chuva (h)  [I$_r$ = 864 mm / duração]", ylabel="FS mínimo",
                   title="A6: mesmo volume (864 mm), durações diferentes")
        axs[1].set(xlabel="duração da chuva (h)", ylabel="volume infiltrado (m³/m)", title="Infiltração")
        axs[2].set(xlabel="I$_r$/k$_s$", ylabel="FS mínimo", title="FS mínimo × I$_r$/k$_s$ (mesmo volume)")
        axs[0].legend(frameon=False, fontsize=7)
        axs[2].legend(frameon=False)
        fig.tight_layout()
        fig.savefig(os.path.join(a.figuras, "contra_A6_volume.png"), dpi=150)
        plt.close(fig)
        tabela("A6", sorted(rs, key=lambda r: (r["solo"], r["train"])), ["solo", "train", "I"])
    # A1c: critério de ruptura FS <= 1 (intensidade crítica e instante em que o FS chega a 1)
    rs = [r for r in R if r["item"] == "A1c"]
    if rs:
        fig, axs = plt.subplots(1, 2, figsize=(12, 4.3))
        linhas = []
        for k, (rot, sel) in enumerate(MODELOS_A1C):
            q = sorted((r for r in rs if sel(r)), key=lambda r: r["I"])
            if not q:
                continue
            I = np.array([r["I"] for r in q])
            mn = np.array([r["fs_min_srm"] for r in q])
            axs[0].plot(I, mn, "o-", color=cores[k], label=rot)
            t1 = np.array([t_cruza(r, 1.0) for r in q])
            t101 = np.array([t_cruza(r, 1.01) for r in q])
            axs[1].plot(I, t1, "o-", color=cores[k], label=rot)
            axs[1].plot(I, t101, "s:", color=cores[k], mfc="w")
            for r, a1, a2 in zip(q, t1, t101):
                linhas.append((rot, r["I"], r["fs_min_srm"], r["t_min_srm"], a1, a2, r["t_fail"]))
            for nivel in (1.0, 1.01):
                linhas.append((rot, f"I_c(FS={nivel:g})", intensidade_critica(I, mn, nivel), "", "", "", ""))
        axs[0].axhline(1, color="#999", lw=0.6, ls=":")
        axs[0].set(xlabel="I$_r$ (mm/h)", ylabel="FS mínimo (tolerância do SRM 0,002)",
                   title="A1c: FS mínimo × I$_r$ (f$_{50}$, α = 45°, 24 h)")
        axs[1].set(xlabel="I$_r$ (mm/h)", ylabel="t (h)", title="instante em que o FS chega a 1 (cheio) e a 1,01 (vazio)",
                   ylim=(0, 30))
        axs[0].legend(frameon=False, fontsize=7)
        fig.tight_layout()
        fig.savefig(os.path.join(a.figuras, "contra_A1_criterio.png"), dpi=150)
        plt.close(fig)
        with open(os.path.join(a.figuras, "contra_tabela_A1c.csv"), "w") as f:
            f.write("modelo,I,fs_min,t_min,t_fs1,t_fs1.01,t_fail\n")
            for l in linhas:
                f.write(",".join(f"{x:.4g}" if isinstance(x, float) else str(x) for x in l) + "\n")

    # A3c: anisotropia no modelo completo
    rs = [r for r in R if r["item"] == "A3c"]
    if rs:
        fig, ax = plt.subplots(figsize=(6.5, 4.3))
        A3 = [r for r in R if r["item"] == "A3"]
        for k, th in enumerate((None, 0.0, 45.0, 135.0)):
            def sel(r):
                return r["kx"] == r["ky"] if th is None else (r["kx"] != r["ky"] and r["theta"] == th)
            q = sorted((r["I"], r["fs_min_srm"]) for r in rs if sel(r))
            rot = "isotrópico 1e-5" if th is None else f"θ = {th:g}°"
            if q:
                ax.semilogx(*zip(*q), "o-", color=cores[k + 1 if th is not None else 4], label=f"acoplado, ar ativo, {rot}")
            q = sorted((r["I"], r["fs_min_srm"]) for r in A3 if th is not None and r["theta"] == th)
            if q:
                ax.semilogx(*zip(*q), "^--", color=cores[k + 1], mfc="w", lw=0.8, label=f"desacoplado, p$_a$ = 0, {rot}")
        ax.axhline(1, color="#999", lw=0.6, ls=":")
        ax.set(xlabel="I$_r$ (mm/h)", ylabel="FS mínimo", title="A3c: anisotropia no modelo acoplado (f$_{50}$, α = 45°)")
        ax.legend(frameon=False, fontsize=6.5)
        fig.tight_layout()
        fig.savefig(os.path.join(a.figuras, "contra_A3_acoplado.png"), dpi=150)
        plt.close(fig)
        for r in rs:
            r["orient"] = "iso" if r["kx"] == r["ky"] else f"th{r['theta']:g}"
        tabela("A3c", sorted(rs, key=lambda r: (r["orient"], r["I"])), ["orient", "I"])

    # A4c: lei de chi
    rs = [r for r in R if r["item"] == "A4c"]
    if rs:
        A4 = [r for r in R if r["item"] == "A4"]
        fig, axs = plt.subplots(1, 2, figsize=(12, 4.2))
        for k, smax in enumerate((75.0, 150.0, 300.0)):
            for grupo, ls, mk in ((rs, "-", "o"), (A4, "--", "s")):
                lei = "χ = S$_e$" if grupo is rs else "χ = 1"
                for I in (9.0, 36.0):
                    q = sorted((r["hw"], r["fs"][0], r["fs_min_srm"]) for r in grupo if r["smax"] == smax and r["I"] == I)
                    if not q:
                        continue
                    hw, ini, mn = map(np.array, zip(*q))
                    if I == 36.0:
                        axs[0].plot(hw, ini, mk + ls, color=cores[k], mfc="w", label=f"inicial, {lei}, s$_{{max}}$ = {smax:g}")
                        axs[0].plot(hw, mn, mk + ls, color=cores[k], label=f"mínimo 36 mm/h, {lei}, s$_{{max}}$ = {smax:g}")
                    axs[1].plot(hw, ini - mn, mk + ls, color=cores[k], mfc="w" if I == 9.0 else cores[k],
                                label=f"{I:g} mm/h, {lei}, s$_{{max}}$ = {smax:g}")
        axs[0].set(xlabel="H$_w$ (m)", ylabel="FS", title="A4c: FS inicial (vazio) e mínimo com 36 mm/h (cheio)")
        axs[1].set(xlabel="H$_w$ (m)", ylabel="FS inicial − FS mínimo", title="Queda do FS pela chuva")
        axs[0].legend(frameon=False, fontsize=5.5, ncol=2)
        axs[1].legend(frameon=False, fontsize=5.5, ncol=2)
        fig.tight_layout()
        fig.savefig(os.path.join(a.figuras, "contra_A4_chi.png"), dpi=150)
        plt.close(fig)
        tabela("A4c", sorted(rs, key=lambda r: (r["smax"], r["I"], r["hw"])), ["smax", "I", "hw"])
    print(f"{len(R)} casos lidos")


KS_MMH = {"f10": 360.0, "f50": 36.0, "f100": 3.6}
MODELOS_A1C = [("completo, ψ = 0", lambda r: r["acoplado"] == 1 and r["psi"] == 0),
               ("completo, ψ = φ'/2", lambda r: r["acoplado"] == 1 and r["psi"] == 13),
               ("completo, ψ = φ'", lambda r: r["acoplado"] == 1 and r["psi"] == 26),
               ("artigo (desacoplado, p$_a$ = 0)", lambda r: r["acoplado"] == 0)]


def t_cruza(r, nivel):
    """primeiro instante (interpolação linear entre as saídas do FS) em que o FS chega a nivel; nan se não chega"""
    t, fs = r["t"], r["fs"]
    for k in range(1, len(t)):
        if fs[k] <= nivel < fs[k - 1]:
            return float(t[k - 1] + (fs[k - 1] - nivel) / (fs[k - 1] - fs[k]) * (t[k] - t[k - 1]))
    return float("nan")


def intensidade_critica(I, fsmin, nivel):
    """menor I_r com FS mínimo <= nivel, por interpolação linear entre as intensidades calculadas"""
    for k in range(len(I)):
        if fsmin[k] <= nivel:
            if k == 0:
                return float(I[0])
            return float(I[k - 1] + (fsmin[k - 1] - nivel) / (fsmin[k - 1] - fsmin[k]) * (I[k] - I[k - 1]))
    return float("nan")


def main():
    aqui = os.path.dirname(os.path.abspath(__file__))
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("acao", choices=["rodar", "figuras"])
    p.add_argument("--itens", nargs="*", default=ITENS, choices=ITENS)
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
