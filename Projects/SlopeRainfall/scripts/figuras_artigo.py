"""Figuras e tabelas do artigo sobre a intensidade crítica (docs/biot_nao_saturado/artigo).

Lê os casos de artigo.py (resultados/artigo), de contraexemplos.py (A1, A1c), de fase_ar.py (G1, G3, G4), da Série B
(variante iii, resultados/B_passivo_desacoplado) e da revisão (resultados/artigo_rev: FS mínimo sem parada antecipada,
Liakopoulos com a curva de retenção de Callari e Abati 2009, Eq. 58). A intensidade crítica I_c,delta é calculada para
os limiares FS <= 1 + delta de DELTAS.

Uso: python3 figuras_artigo.py [--runs resultados] [--figuras <pasta>]
"""
import argparse
import csv
import glob
import os

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import matplotlib.ticker

AQUI = os.path.dirname(os.path.abspath(__file__))
REPO_RES = os.path.normpath(os.path.join(AQUI, "..", "resultados"))
LIMIAR = 1.01           # FS que define a ruptura (delta = 0,01, com a resolução do SRM e da malha)
DELTAS = (0.02, 0.01, 0.005, 0.0)
plt.rcParams.update({"font.size": 8.5, "axes.linewidth": 0.6, "legend.fontsize": 7, "axes.titlesize": 8.5,
                     "lines.linewidth": 1.1, "savefig.bbox": "tight"})
COR = {"i": "#b2182b", "i0": "#e08214", "i13": "#f4a582", "ii": "#2166ac", "ii0": "#92c5de", "iii": "#1a1a1a",
       "iv": "#5aae61"}
NOME = {"i": r"(i) coupled, active air, $\psi=\phi'$", "i0": r"(i) coupled, active air, $\psi=0$",
        "i13": r"(i) coupled, active air, $\psi=\phi'/2$",
        "ii": r"(ii) coupled, passive air, $\psi=\phi'$", "ii0": r"(ii) coupled, passive air, $\psi=0$",
        "iii": "(iii) uncoupled, passive air", "iv": "(iv) uncoupled, active air"}


def ler(pasta):
    p = os.path.join(pasta, "summary.txt")
    if not os.path.exists(p):
        return None
    s = {}
    for linha in open(p):
        k, v = linha.split(None, 1)
        try:
            s[k] = float(v)
        except ValueError:
            s[k] = v.strip()
    fs = np.atleast_2d(np.loadtxt(os.path.join(pasta, "fs.csv"), delimiter=",", skiprows=1))
    s["t"], s["fs"] = fs[:, 0], fs[:, 1]
    s["pasta"] = pasta
    s["t_cruza"] = cruzamento(s["t"], s["fs"])
    # rodada encerrada no primeiro FS <= parar: o FS mínimo dela só vale para limiares >= parar
    s["parado"] = bool(s.get("parar", 0)) and s["t"][-1] < s.get("ttotal", 0) - 1e-6
    return s


def cruzamento(t, fs, lim=LIMIAR):
    """primeiro instante com FS <= lim (interpolação linear entre as saídas); -1 se não cruza"""
    for k in range(len(fs)):
        if fs[k] <= lim:
            if k == 0:
                return t[0]
            return t[k - 1] + (fs[k - 1] - lim) / (fs[k - 1] - fs[k]) * (t[k] - t[k - 1])
    return -1.


def intensidade_critica(pares, lim=LIMIAR):
    """I_c por interpolação linear do FS mínimo entre a maior intensidade acima do limiar e a menor abaixo dele"""
    pares = sorted(pares)
    for (i0, f0), (i1, f1) in zip(pares[:-1], pares[1:]):
        if f0 > lim >= f1:
            return i0 + (f0 - lim) / (f0 - f1) * (i1 - i0), (i0, i1)
    if pares and pares[0][1] <= lim:
        return -1., (None, pares[0][0])                   # abaixo da menor intensidade calculada
    return np.inf, (pares[-1][0] if pares else None, None)  # acima da maior


def ic_delta(L, d):
    """I_c,delta com L = [(I, FS mínimo, parado)], interpolando o FS mínimo entre as intensidades que cercam 1 + delta.
    Uma rodada parada no primeiro FS <= 1,01 tem FS mínimo verdadeiro <= o registrado: se o registrado já é <= 1 + delta,
    a interpolação dá um limite superior de I_c (devolvido com cota=True); se não, I_c não é determinável (None).
    Devolve (valor, (I abaixo, I acima), cota)."""
    lim = 1.0 + d
    acima = None
    for I, F, parado in sorted(L):
        if F > lim:
            if parado:
                return None, (acima[0] if acima else None, I), False
            acima = (I, F)
            continue
        if acima is None:
            return -1., (None, I), parado
        return acima[0] + (acima[1] - lim) / (acima[1] - F) * (I - acima[0]), (acima[0], I), parado
    return np.inf, (acima[0] if acima else None, None), False


def texto_ic(v, br, cota=False):
    if v is None:
        return "n.c."
    if np.isinf(v):
        return f"$>{br[0]:g}$"
    if v < 0:
        return f"$<{br[1]:g}$"
    return (r"$\le$" if cota else "") + f"{v:.0f}"


# --------------------------------------------------------------------------------------------------- casos
def series(R):
    """{modelo: [(I, FS mínimo, t_cruza, t_f, caso)]} do f50 (alpha = 45 graus, chi = 1, malha padrão)"""
    art, con, far = R["artigo"], R["contra"], R["fase_ar"]
    S = {k: [] for k in ("i", "i0", "i13", "ii", "ii0", "iii", "iv")}

    def add(k, I, s):
        if s is not None:
            S[k].append((I, s["fs_min_srm"], s["t_cruza"], s.get("t_fail", -1.), s))

    for p in glob.glob(os.path.join(con, "A1c_f50_a45_Hw5_I*_tr24*")):
        b = os.path.basename(p)
        I = float(b.split("_I")[1].split("_")[0])
        k = "i0" if b.endswith("psi0") else ("i" if b.endswith("psi26") else "iii")
        add(k, I, ler(p))
    for p in glob.glob(os.path.join(con, "A1_f50_a45_Hw5_I*_tr24_psi*")):
        b = os.path.basename(p)
        I = float(b.split("_I")[1].split("_")[0])
        if b.endswith("psi26") and not any(abs(x[0] - I) < 1e-9 for x in S["i"]):
            add("i", I, ler(p))
        elif b.endswith("psi13"):
            add("i13", I, ler(p))
        elif b.endswith("psi0") and I > 54:
            add("i0", I, ler(p))
    for p in glob.glob(os.path.join(REPO_RES, "B_passivo_desacoplado", "f50_a45_H10_Hw5_I*")):
        I = float(os.path.basename(p).split("_I")[1])
        if not any(abs(x[0] - I) < 1e-9 for x in S["iii"]):
            add("iii", I, ler(p))
    rev = R.get("rev")
    novos = {}
    if rev:
        for p in glob.glob(os.path.join(rev, "D_*_f50_a45_I*_semparar")):
            b = os.path.basename(p)
            s = ler(p)
            if s is not None:
                novos[(b.split("_")[1], float(b.split("_I")[1].split("_")[0]))] = s
    for p in glob.glob(os.path.join(art, "D_*_f50_a45_I*")):
        b = os.path.basename(p)
        I = float(b.split("_I")[1].split("_")[0])
        m = b.split("_")[1]
        k = m + ("0" if b.endswith("psi0") else "")
        if k == m and (m, I) in novos:   # caso refeito sem parar no primeiro FS <= 1,01
            continue
        add(k, I, ler(p))
    for (m, I), s in novos.items():
        add(m, I, s)
    for k in S:
        S[k].sort(key=lambda x: x[0])
    return S


# --------------------------------------------------------------------------------------------------- Liakopoulos
def perfil(pasta, arq, t_h):
    V = np.loadtxt(os.path.join(pasta, "vertices.csv"), delimiter=",", skiprows=1)
    A = np.atleast_2d(np.loadtxt(os.path.join(pasta, arq), delimiter=","))
    k = int(np.argmin(np.abs(A[:, 0] - t_h)))
    m = np.abs(V[:, 0]) < 1e-9
    o = np.argsort(V[m, 1])
    return V[m, 1][o], A[k, 1:][m][o]


def fig_liakopoulos(R, saida):
    """Coluna de Liakopoulos: modelo bifásico (i) e ar passivo (ii) com a retenção e a permeabilidade de Callari e Abati
    (2009, Eq. 58), contra a vazão e as pressões medidas (digitalizadas das Figs. 3a e 3b de Callari e Abati 2009)."""
    base = R.get("rev") or R["artigo"]
    pa_, pp_ = os.path.join(base, "L_liakopoulos_i_callari"), os.path.join(base, "L_liakopoulos_ii_callari")
    if not (os.path.exists(os.path.join(pa_, "summary.txt")) and os.path.exists(os.path.join(pp_, "summary.txt"))):
        print("Liakopoulos: casos ainda não calculados")
        return
    ref = os.path.join(AQUI, "ref")
    exp_q = np.loadtxt(os.path.join(ref, "liakopoulos_vazao_exp.csv"), delimiter=",", comments="#")
    exp_p = np.loadtxt(os.path.join(ref, "liakopoulos_pw_exp.csv"), delimiter=",", comments="#")
    fig, ax = plt.subplots(1, 4, figsize=(7.4, 2.35))
    cores = plt.cm.viridis(np.linspace(0, 0.9, 6))
    tmin = [5, 10, 20, 30, 60, 120]
    for j, tm in enumerate(tmin):
        z, p = perfil(pa_, "snaps_pw.csv", tm / 60)
        ax[0].plot(p, z, color=cores[j], label=f"{tm} min")
        z, a = perfil(pa_, "snaps_pa.csv", tm / 60)
        ax[1].plot(a, z, color=cores[j])
    ax[0].set(xlabel="$p_w$ (kPa)", ylabel="height (m)", title="(a) water pressure, two-phase")
    ax[0].legend(frameon=False, fontsize=6, loc="lower left")
    ax[1].set(xlabel="$p_a$ (kPa)", title="(b) air pressure, two-phase")
    erros = []
    for j, (ts, tm) in enumerate(((600, 10), (1200, 20), (1800, 30))):
        c = plt.cm.plasma(0.15 + 0.3 * j)
        z, p = perfil(pp_, "snaps_pw.csv", tm / 60)
        ax[2].plot(p, z, color=c, label=f"{ts} s, passive air")
        e = exp_p[np.abs(exp_p[:, 0] - ts) < 1]
        d = np.interp(e[:, 1], z, p) - e[:, 2]
        z2, p2 = perfil(pa_, "snaps_pw.csv", tm / 60)
        ax[2].plot(p2, z2, color=c, ls="--", lw=0.9)
        d2 = np.interp(e[:, 1], z2, p2) - e[:, 2]
        erros.append((ts, np.sqrt(np.mean(d ** 2)), np.abs(d).max(), np.sqrt(np.mean(d2 ** 2)), np.abs(d2).max()))
        ax[2].plot(e[:, 2], e[:, 1], "s", mfc="none", color=c, ms=3.5)
    ax[2].set(xlabel="$p_w$ (kPa)", title="(c) water pressure vs. test")
    ax[2].legend(frameon=False, fontsize=6, loc="lower left")
    Q = {}
    for pasta, ls, lab in ((pp_, "-", "passive air"), (pa_, "--", "two-phase")):
        h = np.loadtxt(os.path.join(pasta, "hist.csv"), delimiter=",", skiprows=1)
        Q[lab] = (h[:, 0] * 3600, -h[:, -1] / 0.1 * 1e6)
        ax[3].plot(*Q[lab], color="k", ls=ls, label=lab)
    ax[3].plot(exp_q[:, 0], exp_q[:, 1], "s", mfc="none", color="#b2182b", ms=3.5, label="test")
    ax[3].set(xlabel="t (s)", ylabel=r"outflow ($10^{-6}$ m/s)", title="(d) outflow at the base", xlim=(0, 7500))
    ax[3].legend(frameon=False, fontsize=6)
    for a in ax[1:3]:
        a.set_yticklabels([])
    for a in ax[:3]:
        a.set_ylim(0, 1)
    fig.tight_layout(w_pad=0.6)
    fig.savefig(os.path.join(saida, "liakopoulos.pdf"))
    plt.close(fig)
    linhas = []
    for lab, (t, q) in Q.items():
        e = np.interp(exp_q[:, 0], t, q) - exp_q[:, 1]
        linhas.append((lab, "vazão RMS", np.sqrt(np.mean(e ** 2)), "máx", np.abs(e).max()))
    tp, qp = Q["passive air"]
    tb, qb = Q["two-phase"]
    for ts in (600, 1200, 1800):
        linhas.append(("razão vazão passivo/bifásico", ts, np.interp(ts, tp, qp) / np.interp(ts, tb, qb)))
    for e in erros:
        linhas.append(("p_w", e[0], "RMS passivo", e[1], "máx", e[2], "RMS bifásico", e[3], "máx", e[4]))
    return linhas


# --------------------------------------------------------------------------------------------------- decomposição
def fig_decomposicao(R, saida):
    S = series(R)
    fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.9))
    linhas = []
    for k in ("iii", "iv", "ii", "ii0", "i", "i13", "i0"):
        L = S[k]
        if not L:
            continue
        I = np.array([x[0] for x in L])
        F = np.array([x[1] for x in L])
        par = np.array([x[4]["parado"] for x in L])
        ax[0].plot(I, F, "-", color=COR[k], label=NOME[k])
        ax[0].plot(I[~par], F[~par], "o", color=COR[k], ms=3)
        ax[0].plot(I[par], F[par], "o", mfc="none", color=COR[k], ms=3)
        Ic, br = intensidade_critica([(x[0], x[1]) for x in L])
        if np.isfinite(Ic) and Ic > 0:
            ax[0].plot(Ic, LIMIAR, "v", color=COR[k], ms=6)
        tc = [(x[0], x[3] if x[3] >= 0 and (x[2] < 0 or x[3] < x[2]) else x[2]) for x in L if x[2] >= 0 or x[3] >= 0]
        if tc:
            ax[1].plot([t[0] for t in tc], [t[1] for t in tc], "o-", color=COR[k], ms=3)
        linhas.append((k, Ic, br, L))
    ax[0].axhline(LIMIAR, color="0.5", lw=0.6, ls=":")
    ax[0].set(xscale="log", xlabel="$I_r$ (mm/h)", ylabel="minimum FS (24-h rainfall)", ylim=(0.8, 2.1),
              title="(a) minimum FS ($\\blacktriangledown$: $I_c$)")
    ax[0].legend(frameon=False, fontsize=6.0, loc="upper right")
    ax[1].set(xscale="log", xlabel="$I_r$ (mm/h)", ylabel="time to failure (h)", ylim=(0, 30),
              title="(b) first FS $\\leq$ 1.01 or loss of equilibrium")
    for a in ax:
        a.set_xticks([5, 10, 20, 50, 100, 200])
        a.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
        a.get_xaxis().set_minor_formatter(matplotlib.ticker.NullFormatter())
    fig.tight_layout()
    fig.savefig(os.path.join(saida, "decomposicao.pdf"))
    plt.close(fig)
    return linhas


def fig_fs_tempo(R, saida):
    S = series(R)
    fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.4), sharey=True)
    for a, k in zip(ax, ("iii", "i", "i0")):
        cores = plt.cm.viridis(np.linspace(0, 0.9, 8))
        Ls = [x for x in S[k] if x[0] >= 18]
        for j, x in enumerate(Ls[:8]):
            s = x[4]
            a.plot(s["t"], s["fs"], color=cores[j], label=f"{x[0]:g} mm/h")
            if k == "i0" and s.get("t_fail", -1) > 0:
                a.plot(s["t_fail"], np.interp(s["t_fail"], s["t"], s["fs"]), "x", color=cores[j], ms=5)
        a.axhline(1, color="0.5", lw=0.6, ls=":")
        a.axvline(24, color="0.8", lw=0.6)
        a.set(xlabel="t (h)", title=NOME[k], xlim=(0, 48), ylim=(0.8, 2.2))
        a.legend(frameon=False, fontsize=5.8, loc="upper right")
    ax[0].set_ylabel("FS")
    fig.tight_layout(w_pad=0.4)
    fig.savefig(os.path.join(saida, "fs_tempo.pdf"))
    plt.close(fig)


def tabela_decomposicao(R, saida):
    S = series(R)
    ks = 36.0
    rot = {"iii": r"(iii) uncoupled, passive air", "iv": r"(iv) uncoupled, active air",
           "ii": r"(ii) coupled, passive, $\psi=\phi'$", "ii0": r"(ii) coupled, passive, $\psi=0$",
           "i": r"(i) coupled, active, $\psi=\phi'$", "i13": r"(i) coupled, active, $\psi=\phi'/2$",
           "i0": r"(i) coupled, active, $\psi=0$"}
    linhas, dados = [], {}
    for k in ("iii", "iv", "ii", "ii0", "i", "i13", "i0"):
        L = S[k]
        if not L:
            continue
        L3 = [(x[0], x[1], x[4]["parado"]) for x in L]
        ics = {d: ic_delta(L3, d) for d in DELTAS}
        dados[k] = {"ic": {d: v[0] for d, v in ics.items()}, "L": [(x[0], x[1], x[4]["parado"], x[2], x[3]) for x in L]}
        Ic, (a, b), _ = ics[0.01]
        br = f"{a:g}--{b:g}" if a is not None and b is not None else ""
        rat = f"{Ic / ks:.2f}" if Ic is not None and np.isfinite(Ic) and Ic > 0 else ""
        sem = [x[1] for x in L if not x[4]["parado"]]
        fmin = f"{min(sem):.3f}" if sem else ""

        def tempo(I):
            for x in L:
                if abs(x[0] - I) < 1e-9:
                    tc, tf = x[2], x[3]
                    if tf > 0 and (tc < 0 or tf < tc):
                        return f"{tf:.1f}$^\\dagger$"
                    return f"{tc:.1f}" if tc >= 0 else "--"
            return ""
        linhas.append(f"{rot[k]} & " + " & ".join(texto_ic(*ics[d]) for d in DELTAS) + f" & {br} & {rat} & {fmin} & "
                      + " & ".join(tempo(I) for I in (36.0, 54.0, 90.0)) + "\\\\")
    with open(os.path.join(saida, "tabela_decomposicao.tex"), "w") as f:
        f.write("\\begin{tabular}{lccccccccc c}\n\\toprule\n")
        f.write("model & \\multicolumn{4}{c}{$I_{c,\\delta}$ (mm/h) for $\\delta=$} & bracket & $I_c/k_s$ & lowest & "
                "\\multicolumn{3}{c}{time to failure (h), $I_r=$}\\\\\n")
        f.write(" & " + " & ".join(f"{d:g}" for d in DELTAS) + " & ($\\delta=0.01$) & & min.\\ FS & 36 & 54 & 90\\\\\n"
                "\\midrule\n")
        f.write("\n".join(linhas) + "\n\\bottomrule\n\\end{tabular}\n")
    return linhas, dados


def generalizacao(R):
    """{(caso, modelo): [(I, FS mínimo, parado)]} do grupo G, da Série B (iii, f10) e das rodadas sem parada da revisão"""
    art, rev = R["artigo"], R.get("rev")
    out = {}

    def junta(chave, padrao, filtro=lambda b: True, extra=()):
        d = {}
        for p in glob.glob(os.path.join(art, padrao)):
            b = os.path.basename(p)
            if not filtro(b):
                continue
            s = ler(p)
            if s:
                d[float(b.split("_I")[1].split("_")[0])] = (s["fs_min_srm"], s["parado"])
        for I, p in extra:
            s = ler(p)
            if s:
                d[I] = (s["fs_min_srm"], s["parado"])
        out[chave] = sorted((I, v[0], v[1]) for I, v in d.items())

    revf = (lambda n: os.path.join(rev, n)) if rev else (lambda n: "")
    junta(("f10", "i"), "G_i_f10_a45_I*", lambda b: not b.endswith("chise"),
          [(140.0, revf("G_i_f10_a45_I140_semparar"))])
    junta(("f10", "iii"), "G_iii_f10_a45_I*", lambda b: not b.endswith("chise"),
          [(I, os.path.join(REPO_RES, "B_passivo_desacoplado", f"f10_a45_H10_Hw5_I{I:g}")) for I in (90.0, 180.0, 360.0)]
          + [(140.0, revf("G_iii_f10_a45_I140_semparar"))])
    junta(("f50_se", "i"), "G_i_f50_a45_I*_chise", extra=[(30.0, revf("G_i_f50_a45_I30_chise_semparar"))])
    junta(("f50_se", "iii"), "G_iii_f50_a45_I*_chise", extra=[(54.0, revf("G_iii_f50_a45_I54_chise_semparar"))])
    far = R["fase_ar"]
    for modo, m in (("artigo", "iii"), ("ativo", "i")):
        L = []
        for I in (9, 18, 36, 80):
            x = ler(os.path.join(far, f"f100_a63.4_H10_I{I}_{modo}_kra1e-4"))
            if x:
                L.append((float(I), x["fs_min_srm"], x["parado"]))
        out[("f100_63", m)] = L
        L = []
        for I in (9, 36, 80):
            x = ler(os.path.join(far, f"f100_a45_H20_I{I}_{modo}_kra1e-4"))
            if x:
                L.append((float(I), x["fs_min_srm"], x["parado"]))
        out[("f100_H20", m)] = L
    return out


def tabela_generalizacao(R, saida):
    S = series(R)
    G = generalizacao(R)
    G[("f50", "iii")] = [(x[0], x[1], x[4]["parado"]) for x in S["iii"]]
    G[("f50", "i")] = [(x[0], x[1], x[4]["parado"]) for x in S["i"]]
    casos = [("f50", r"$f_{50}$, $45^\circ$, $\chi=1$ (reference)", 36),
             ("f10", r"$f_{10}$, $45^\circ$, $\chi=1$", 360),
             ("f50_se", r"$f_{50}$, $45^\circ$, $\chi=S_e$", 36),
             ("f100_63", r"$f_{100}$, $63.4^\circ$, $\chi=1$", 3.6),
             ("f100_H20", r"$f_{100}$, $45^\circ$, $H_s=20$~m, $\chi=1$", 3.6)]
    linhas = []
    for c, rot, ks in casos:
        cel = []
        for m in ("iii", "i"):
            L = G.get((c, m), [])
            cel += [texto_ic(*ic_delta(L, d)) for d in DELTAS]
            sem = [x[1] for x in L if not x[2]]
            cel.append(f"{min(sem):.3f}" if sem else "")
        linhas.append((rot, ks, cel, {m: G.get((c, m), []) for m in ("iii", "i")}))
    with open(os.path.join(saida, "tabela_generalizacao.tex"), "w") as f:
        f.write("\\begin{tabular}{lc ccccc ccccc}\n\\toprule\n")
        f.write(" & & \\multicolumn{5}{c}{model (iii)} & \\multicolumn{5}{c}{model (i), $\\psi=\\phi'$}\\\\\n")
        f.write("\\cmidrule(lr){3-7}\\cmidrule(lr){8-12}\n")
        f.write("slope and soil & $k_s$ & " + " & ".join(f"{d:g}" for d in DELTAS) + " & min.\\ FS & "
                + " & ".join(f"{d:g}" for d in DELTAS) + " & min.\\ FS\\\\\n\\midrule\n")
        for rot, ks, cel, _ in linhas:
            f.write(f"{rot} & {ks:g} & " + " & ".join(cel) + "\\\\\n")
        f.write("\\bottomrule\n\\end{tabular}\n")
    return linhas


# --------------------------------------------------------------------------------------------------- ar aprisionado
def coluna(pasta, t_h, x0, y0, prof=10.0, n=60):
    V = np.loadtxt(os.path.join(pasta, "vertices.csv"), delimiter=",", skiprows=1)
    T = np.loadtxt(os.path.join(pasta, "triangles.csv"), delimiter=",", skiprows=1, dtype=int)
    tri = mtri.Triangulation(V[:, 0], V[:, 1], T)
    out = []
    for arq in ("snaps_pw.csv", "snaps_pa.csv"):
        A = np.atleast_2d(np.loadtxt(os.path.join(pasta, arq), delimiter=","))
        k = int(np.argmin(np.abs(A[:, 0] - t_h)))
        f = mtri.LinearTriInterpolator(tri, A[k, 1:])
        d = np.linspace(0.01, prof, n)
        out.append(np.asarray(f(np.full(n, x0), y0 - d).filled(np.nan)))
    return d, out[0], out[1]


def fig_ar(R, saida):
    far = R["fase_ar"]
    nomes = {"artigo": ("iii", "(iii) uncoupled, passive air"), "passivo": ("ii", "(ii) coupled, passive air"),
             "ativo": ("i", "(i) coupled, active air")}
    casos = {m: os.path.join(far, f"f100_a63.4_H10_I36_{m}_kra1e-4") for m in nomes}
    if not all(os.path.exists(os.path.join(p, "vertices.csv")) for p in casos.values()):
        print("ar: instantâneos do G1 ausentes (figura mantida)")
        return
    # vertical sob o meio da face (alpha = 63,4 graus, H = 10 m): x = 30 + 5/tan(alpha), y = 25
    xm = 30 + 5 / np.tan(np.radians(63.4))
    fig, ax = plt.subplots(1, 4, figsize=(7.4, 2.4))
    for m, p in casos.items():
        k, lab = nomes[m]
        s = ler(p)
        ax[0].plot(s["t"], s["fs"], color=COR[k], label=lab)
        h = np.loadtxt(os.path.join(p, "hist.csv"), delimiter=",", skiprows=1)
        ax[1].plot(h[:, 0], np.cumsum(h[:, 4] * h[:, 1] * 3600), color=COR[k])
        d, pw, pa = coluna(p, 24.0, xm, 25.0)
        ax[2].plot(pw, -d, color=COR[k])
        if m == "ativo":
            ax[2].plot(pa, -d, color=COR[k], ls="--")
        ax[3].plot(pa - pw, -d, color=COR[k])
    ax[0].axhline(1, color="0.5", lw=0.6, ls=":")
    ax[0].set(xlabel="t (h)", ylabel="FS", title="(a) FS, 36 mm/h", xlim=(0, 48))
    ax[0].legend(frameon=False, fontsize=5.8)
    ax[1].set(xlabel="t (h)", ylabel="infiltrated volume (m$^3$/m)", title="(b) infiltration", xlim=(0, 48))
    ax[2].set(xlabel="pressure (kPa)", ylabel="depth below the face (m)", title="(c) $p_w$ (solid), $p_a$ (dashed)")
    ax[3].set(xlabel="$p_a-p_w$ (kPa)", title="(d) matric suction")
    ax[3].set_yticklabels([])
    fig.tight_layout(w_pad=0.5)
    fig.savefig(os.path.join(saida, "ar_aprisionado.pdf"))
    plt.close(fig)


# --------------------------------------------------------------------------------------------------- malha
def fig_malha(R, saida):
    far, art = R["fase_ar"], R["artigo"]
    pts = {"G1": [], "f50_i": [], "f50_iii": []}
    for h, suf in ((0.75, "_m0.75"), (0.5, "_m0.5"), (0.35, "_m0.35")):
        s = ler(os.path.join(far, f"f100_a63.4_H10_I36_ativo_kra1e-4{suf}"))
        if s:
            pts["G1"].append((h, s["fs_min_srm"], s["fs_ini_srm"]))
    base = {"f50_i": os.path.join(R["contra"], "A1c_f50_a45_Hw5_I30_tr24_psi26"),
            "f50_iii": os.path.join(R["contra"], "A1c_f50_a45_Hw5_I48_tr24")}
    fina = {"f50_i": os.path.join(art, "M_i_f50_a45_I30_m0.35"), "f50_iii": os.path.join(art, "M_iii_f50_a45_I48_m0.35")}
    for k in base:
        for h, p in ((0.75, base[k]), (0.35, fina[k])):
            s = ler(p)
            if s:
                pts[k].append((h, s["fs_min_srm"], s["fs_ini_srm"]))
            elif os.path.exists(os.path.join(p, "fs.csv")):    # caso encerrado depois do último FS (PARADO.txt)
                fs = np.atleast_2d(np.loadtxt(os.path.join(p, "fs.csv"), delimiter=",", skiprows=1))
                pts[k].append((h, fs[:, 1].min(), fs[0, 1]))
    fig, ax = plt.subplots(figsize=(3.4, 2.5))
    rot = {"G1": r"$f_{100}$, 63.4$^\circ$, 36 mm/h, (i)", "f50_i": r"$f_{50}$, 45$^\circ$, 30 mm/h, (i)",
           "f50_iii": r"$f_{50}$, 45$^\circ$, 48 mm/h, (iii)"}
    for k, L in pts.items():
        if L:
            ax.plot([x[0] for x in L], [x[1] for x in L], "o-", ms=3, label=rot[k])
    ax.axhline(1, color="0.5", lw=0.6, ls=":")
    ax.set(xlabel="element size near the surface (m)", ylabel="minimum FS", xlim=(0.8, 0.3))
    ax.legend(frameon=False, fontsize=6.5)
    fig.tight_layout()
    fig.savefig(os.path.join(saida, "malha.pdf"))
    plt.close(fig)
    return pts


# --------------------------------------------------------------------------------------------------- problema
def fig_problema(R, saida):
    p = os.path.join(R["contra"], "A1c_f50_a45_Hw5_I30_tr24_psi26")
    if not os.path.exists(os.path.join(p, "vertices.csv")):
        print("problema: malha da rodada ausente (figura mantida)")
        return
    V = np.loadtxt(os.path.join(p, "vertices.csv"), delimiter=",", skiprows=1)
    T = np.loadtxt(os.path.join(p, "triangles.csv"), delimiter=",", skiprows=1, dtype=int)
    fig, ax = plt.subplots(figsize=(7.0, 2.6))
    ax.triplot(V[:, 0], V[:, 1], T, color="0.6", lw=0.25)
    x = np.array([0, 40, 70.0])
    xc = 40.0
    ax.plot([0, 70], [(20 - 5) - np.tan(np.radians(7)) * (0 - xc), (20 - 5) - np.tan(np.radians(7)) * (70 - xc)],
            "b--", lw=1, label="initial water table")
    for xx in np.linspace(2, 68, 23):
        yy = 30 if xx <= 30 else (20 if xx >= 40 else 30 - (xx - 30))
        ax.annotate("", (xx, yy + 0.3), (xx, yy + 3.0), arrowprops=dict(arrowstyle="->", lw=0.5, color="#2166ac"))
    ax.text(35, 33.8, "rainfall $I_r$ on the horizontal projection; $p_a=0$; $p_w\\leq0$, inflow $\\leq I_r$",
            ha="center", fontsize=6.8)
    ax.text(1, 1.5, "base fixed, closed", fontsize=6.5)
    ax.text(0.5, 8, "rollers;\nhydrostatic\nbelow the WT", fontsize=6.0)
    ax.set_aspect("equal")
    ax.set(xlabel="x (m)", ylabel="y (m)", xlim=(-1, 71), ylim=(-1, 35))
    ax.legend(frameon=False, fontsize=6.5, loc="lower right")
    fig.tight_layout()
    fig.savefig(os.path.join(saida, "problema.pdf"))
    plt.close(fig)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--runs", default=REPO_RES)
    p.add_argument("--figuras", default=os.path.normpath(os.path.join(AQUI, "..", "..", "..", "docs",
                                                                     "biot_nao_saturado", "artigo", "figuras")))
    a = p.parse_args()
    os.makedirs(a.figuras, exist_ok=True)
    R = {"artigo": os.path.join(a.runs, "artigo"), "contra": os.path.join(a.runs, "contraexemplos"),
         "fase_ar": os.path.join(a.runs, "fase_ar"), "rev": os.path.join(a.runs, "artigo_rev")}
    for l in fig_liakopoulos(R, a.figuras) or []:
        print("Liakopoulos", *[f"{v:.3f}" if isinstance(v, float) else v for v in l])
    fig_problema(R, a.figuras)
    fig_fs_tempo(R, a.figuras)
    for k, Ic, br, L in fig_decomposicao(R, a.figuras):
        print(f"{k:4s} I_c = {Ic:.1f} entre {br}: " + "  ".join(f"{x[0]:g}:{x[1]:.4f}{'*' if x[4]['parado'] else ''}"
                                                              for x in L))
    linhas, dados = tabela_decomposicao(R, a.figuras)
    for l in linhas:
        print(l)
    for k, d in dados.items():
        print("I_c,delta", k, "  ".join(f"{dl:g}: {v:.2f}" if v is not None else f"{dl:g}: n.c." for dl, v in d["ic"].items()))
    fig_ar(R, a.figuras)
    print("malha", fig_malha(R, a.figuras))
    for rot, ks, cel, L in tabela_generalizacao(R, a.figuras):
        print("generalização", rot, cel, L)


if __name__ == "__main__":
    main()
