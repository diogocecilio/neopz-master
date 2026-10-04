#!/usr/bin/env python3
"""Comparação com o artigo (Vargas Ceron et al., IJNAMG 2025): tabelas (JSON + markdown) e figuras.

Uso: analise_artigo.py <campanha> <det> <ref> <saída>
  <campanha>: diretório de campanha_artigo.py (um subdiretório por caso, blocos b*.csv e b*.csv.mec)
  <det>:      logs do comando det (fig5/, fig8/, fig9/, conv/, tab56/)
  <ref>:      dados extraídos do artigo (fig5.json, fig8.json, fig9.json, dist.json, fig26_values.json)
  <saída>:    resultados.json, tabelas.md e figuras PNG
"""
import glob
import json
import math
import os
import re
import shutil
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from campanha_artigo import CASOS, ler  # noqa: E402

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

GW = 9.81
# pontos da Fig. 8 sem convergência com a malha (mecanismo translacional raso, φ próximo de β): fora das figuras
NAO_CONV = {"A_b35_hw0"}
# ---------------------------------------------------------------- valores do artigo (Tabelas 3-6 e texto)
# (μ, σ, CoV %, Pf %, CoV(Pf) %)
T_REF = (1.353, 0.318, 23.5, 11.50, 2.77)
ARTIGO_MC = {
    "ref": T_REF,
    "covk0": (1.329, 0.288, 21.7, 11.21, 2.81), "covk75": (1.375, 0.341, 24.8, 11.41, 2.79),
    "covk90": (1.378, 0.354, 25.7, 11.63, 2.76), "covk100": (1.391, 0.360, 25.9, 11.44, 2.78),
    "covc10": (1.367, 0.180, 13.2, 0.45, 4.97), "covc50": (1.324, 0.472, 35.7, 25.78, 1.70),
    "covc70": (1.296, 0.620, 47.9, 36.34, 1.32),
    "covphi5": (1.351, 0.300, 22.2, 10.20, 2.97), "covphi15": (1.357, 0.346, 25.5, 12.71, 2.62),
    "covphi20": (1.361, 0.384, 28.2, 15.18, 2.36),
    "s1.5": (1.356, 0.352, 25.9, 13.94, 2.48), "s2": (1.361, 0.372, 27.3, 14.99, 2.38),
    "s5": (1.350, 0.404, 29.9, 18.93, 2.07), "s10": (1.348, 0.419, 31.1, 19.85, 2.01),
    "s20": (1.355, 0.430, 31.7, 20.23, 1.99), "s400": (1.348, 0.432, 32.1, 21.93, 1.89),
    "alfa2": (1.562, 0.390, 24.9, 3.90, 3.51), "alfa3": (1.704, 0.432, 25.3, 1.83, 4.23),
    "alfa4": (1.821, 0.474, 26.0, 1.00, 4.98), "alfa5": (1.910, 0.504, 26.4, 0.63, 4.74),
    "hw0.5": (1.670, 0.378, 22.6, 1.38, 4.88), "hw0.6": (1.502, 0.345, 22.9, 4.20, 4.78),
    "hw0.7": (1.393, 0.314, 22.5, 8.54, 3.27), "hw0.8": (1.338, 0.309, 23.1, 11.99, 2.71),
    "hw0.9": (1.325, 0.306, 23.1, 12.58, 2.64),
}
ARTIGO_PF_CHO = {"cho_coesivo": (6.5, 7.9), "cho_cphi": (5.50, 6.37)}  # (artigo, Cho 2010)
ARTIGO_DET = {"ref": 1.336, "alfa2": 1.533, "alfa3": 1.674, "alfa4": 1.783, "alfa5": 1.872, "hw0.5": 1.671,
              "hw0.6": 1.494, "hw0.7": 1.383, "hw0.8": 1.322, "hw0.9": 1.307, "cho_coesivo": 1.354,
              "cho_cphi": 1.777}
# caso de MC -> log determinístico na malha do MC (adapt=2, γw = 9.81)
DET_MC = {"ref": "conv/ref_h1_a2", "cho_coesivo": "conv/chocoes_h1_a2", "cho_cphi": "conv/chocphi_h2_a2"}
for _a in (2, 3, 4, 5):
    DET_MC[f"alfa{_a}"] = f"tab56/alfa{_a}_a2"
for _h, _hw in ((0.5, 2.5), (0.6, 3), (0.7, 3.5), (0.8, 4), (0.9, 4.5)):
    DET_MC[f"hw{_h}"] = f"tab56/hw{_hw}_a2"

# observações sobre o artigo (verificadas: ver relatório)
ACHADOS = [
    "Fig. 8: os painéis estão com os solos trocados em relação à Tabela 1. Em h<sub>w</sub>/H = 0 (sem percolação) a "
    "solução log-espiral de Chen (1975) com γ' = 18 − 9,81 dá, para c = 6 kPa e φ = 32° (London clay na Tabela 1), "
    "H<sub>crit</sub> = 13,00 m (β = 60°) e 229,1 m (β = 35°), que são as curvas do painel “Israeli Clay” (12,98 m e "
    "229,3 m); para c = 11,7 kPa e φ = 24,7° dá 17,97 m (β = 60°) e 156,7 m (β = 30°), as do painel “London Clay” "
    "(17,95 m e 156,6 m). O FE reproduz os mesmos valores (diferenças de 0,5–2 %).",
    "γ<sub>w</sub> não é informado; os valores da Fig. 8 em h<sub>w</sub>/H = 0 só são reproduzidos com γ<sub>w</sub> = "
    "9,81 kN/m³ (com 10 kN/m³ a diferença é de 2,5 %). As contas aqui usam 9,81.",
    "Exemplo c-φ de Cho (seção 5.3.2): o texto e a Fig. 19 indicam H = 5 m, mas Γ = 1,777 e F<sub>s</sub> = 1,203/1,204 "
    "correspondem a H = 10 m (Chen dá Γ = 1,7770 com H = 10 m e 3,554 com H = 5 m).",
    "Seção 6.1: o texto fala em correlação “perfeita” de Γ com a coesão, mas a barra da Fig. 26 vale 0,88.",
    "As curvas −grad u'<sub>FE</sub> (Fig. 9) e não K·v'<sub>opt</sub> são as do Γ = 1,336 das Tabelas 5 e 6 "
    "(a Fig. 9 dá 1,31 em β = 45°, α = 1, contra 1,89 para K·v'<sub>opt</sub>).",
]

GRUPOS = [
    ("Tabela 3 — coeficientes de variação", ["ref", "covk0", "covk75", "covk90", "covk100", "covc10", "covc50",
                                             "covc70", "covphi5", "covphi15", "covphi20"]),
    ("Tabela 4 — distâncias de autocorrelação (L_x, L_y) = s (20 m, 2 m)",
     ["ref", "s1.5", "s2", "s5", "s10", "s20", "s400"]),
    ("Tabela 5 — anisotropia α = k_h/k_v", ["ref", "alfa2", "alfa3", "alfa4", "alfa5"]),
    ("Tabela 6 — rebaixamento h_w/H", ["hw0.5", "hw0.6", "hw0.7", "hw0.8", "hw0.9", "ref"]),
    ("Seção 5.3 — Cho (2010), sem percolação", ["cho_coesivo", "cho_cphi"]),
]


def det_log(path):
    """Γ, FS, funcional hidráulico, equações e tempo de um log do comando det"""
    out = {}
    if not os.path.exists(path):
        return out
    txt = open(path, errors="replace").read()
    for key, pat in (("gamma", r"Gamma \(fator de carga\) = ([0-9.eE+-]+)"),
                     ("gamma_sup", r"Gamma \(fator de carga\) = [0-9.eE+-]+ \(colapso em ([0-9.eE+-]+)"),
                     ("fs", r"FS \(reducao de resistencia\) = ([0-9.eE+-]+)"),
                     ("J", r"J\(u_FE\)/\(kh H\^2 gw\^2\) = ([0-9.eE+-]+)"),
                     ("neq", r"\n\s+(\d+) equacoes"),
                     ("tempo", r"tempo total ([0-9.eE+-]+)")):
        m = re.search(pat, txt)
        if m:
            out[key] = float(m.group(1))
    out["limite_maximo"] = "limite_maximo" in txt
    return out


def estat(g):
    g = np.asarray(g, float)
    n = len(g)
    if n == 0:
        return None
    pf = float(np.mean(g < 1.))
    mu, sd = float(np.mean(g)), float(np.std(g, ddof=1)) if n > 1 else 0.
    covpf = math.sqrt((1 - pf) / (n * pf)) if pf > 0 else float("nan")
    # intervalo de confiança de 95 % de Pf (Wilson)
    z = 1.96
    den = 1 + z * z / n
    cen = (pf + z * z / (2 * n)) / den
    hw = z * math.sqrt(pf * (1 - pf) / n + z * z / (4 * n * n)) / den
    lg = np.log(g[g > 0])
    gs = np.sort(g)
    rob = {"q05": float(np.quantile(g, .05)), "q25": float(np.quantile(g, .25)), "q75": float(np.quantile(g, .75)),
           "q95": float(np.quantile(g, .95)), "max": float(gs[-1]),
           "mu_sem_max": float(np.mean(gs[:-1])) if n > 2 else float("nan"),
           "sd_sem_max": float(np.std(gs[:-1], ddof=1)) if n > 2 else float("nan")}
    return {**rob, "N": n, "mu": mu, "sd": sd, "cov": sd / mu, "pf": pf, "covpf": covpf, "pf_lo": cen - hw,
            "pf_hi": cen + hw, "mediana": float(np.median(g)), "mu_ln": float(np.mean(lg)),
            "sd_ln": float(np.std(lg, ddof=1)) if n > 1 else 0.,
            "pf_lognormal": float(0.5 * math.erfc(np.mean(lg) / (np.std(lg, ddof=1) * math.sqrt(2))))
            if n > 1 else float("nan")}


def spearman_parcial(x, y, z, w):
    """ρ_xy,zw de segunda ordem (eqs. 70-71) a partir dos postos"""
    from scipy.stats import spearmanr
    M = np.column_stack([x, y, z, w])
    r = spearmanr(M).correlation

    def p1(a, b, c):  # ρ_ab,c
        return (r[a, b] - r[a, c] * r[b, c]) / math.sqrt((1 - r[a, c] ** 2) * (1 - r[b, c] ** 2))

    def p2(a, b, c, d):  # ρ_ab,cd
        ab, ad, bd = p1(a, b, c), p1(a, d, c), p1(b, d, c)
        return (ab - ad * bd) / math.sqrt((1 - ad ** 2) * (1 - bd ** 2))

    return p2(0, 1, 2, 3), r[0, 1]


def comparacoes(R, ref):
    """Fig. 8 e Fig. 9: FE contra as curvas do artigo (interpoladas nos pontos do FE); Cho: momentos das densidades"""
    D = R["det"]
    try:
        f8 = json.load(open(os.path.join(ref, "fig8.json")))

        def curva(panel, superior, estilo):
            ss = [x for x in f8["series"] if x["panel"] == panel and x["style"].startswith(estilo)]
            ss.sort(key=lambda x: -x["points"][-1][1])
            p = np.array(ss[0 if superior else 1]["points"], float)
            return lambda x: float(np.exp(np.interp(x, p[:, 0], np.log(p[:, 1]))))

        cmp8 = []
        for panel, solo, betas in (("London Clay", "B", (30, 60)), ("Israeli Clay", "A", (35, 60))):
            for i, b in enumerate(betas):
                dd, ds, so = (curva(panel, i == 0, e) for e in ("dash-dot", "dashed", "solid"))
                for hw in (0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5):
                    k = f"{solo}_b{b}_hw{hw:g}"
                    d = D["fig8"].get(k, {})
                    if not d.get("gamma") or d.get("limite_maximo"):
                        continue
                    x = hw / 5.
                    cmp8.append({"painel": panel, "solo": solo, "beta": b, "hwH": x, "fe": 5 * d["gamma"],
                                 "grad_u_fe": dd(x), "v_opt": ds(x), "wu_rp": so(x), "nao_conv": k in NAO_CONV})
        D["fig8_cmp"] = cmp8
    except Exception as e:  # noqa: BLE001
        print("fig8_cmp:", e)
    try:
        f9 = json.load(open(os.path.join(ref, "fig9.json")))
        cmp9 = []
        for s9 in f9["series"]:
            for b, y in s9["points"]:
                if y is None or b % 15:
                    continue
                g = D["fig9"].get(f"beta{int(b)}_alfa{s9['alpha']}", {}).get("gamma")
                if g:
                    cmp9.append({"alfa": s9["alpha"], "beta": int(b), "fe": g, "curva": "grad" if "grad" in s9["name"]
                                 else "vopt", "artigo": y})
        D["fig9_cmp"] = cmp9
    except Exception as e:  # noqa: BLE001
        print("fig9_cmp:", e)
    try:
        dist = json.load(open(os.path.join(ref, "dist.json")))
        for s5 in dist["summary"]:
            c = "cho_coesivo" if "5.3.1" in s5["example"] else "cho_cphi"
            if s5["variable"] == "Gamma" and c in R["mc"]:
                a = R["mc"][c].setdefault("artigo", {})
                a["mu"], a["sd"] = s5["from_pdf"]["mean"], s5["from_pdf"]["std"]
                fc = s5.get("from_cdf", {})
                a["mu_cdf"], a["sd_cdf"] = fc.get("mean_crosscheck"), fc.get("std_crosscheck")
                a["mediana"] = fc.get("median")
                a["quantis"] = {k: fc.get(k) for k in ("q05", "q25", "q75", "q95")}
                a["cov"] = a["sd"] / a["mu"]
                a["fonte_momentos"] = "Fig. 17" if c == "cho_coesivo" else "Fig. 21"
    except Exception as e:  # noqa: BLE001
        print("dist:", e)


def main(camp, det, ref, out):
    os.makedirs(out, exist_ok=True)
    R = {"mc": {}, "det": {}, "notas": []}
    md = []
    # ------------------------------------------------------------ Monte Carlo
    dados = {}
    for c in CASOS:
        d = os.path.join(camp, c)
        if not os.path.isdir(d) and not os.path.exists(d + ".csv"):
            continue
        rows, mec = ler(d)
        if not rows:
            continue
        ks = sorted(rows)
        g = np.array([float(rows[k]["fator"]) for k in ks])
        t = np.array([float(rows[k]["tempo_s"]) for k in ks])
        st = estat(g)
        st["tempo_medio_s"] = float(np.mean(t))
        st["status"] = {s: sum(1 for k in ks if rows[k]["status"] == s) for s in set(r["status"] for r in rows.values())}
        st["S_artigo"] = CASOS[c][1]
        dl = det_log(os.path.join(det, DET_MC.get(c, "conv/ref_h1_a2") + ".log"))
        st["gamma_det_fe"] = dl.get("gamma")
        gd_art = ARTIGO_DET.get(c, ARTIGO_DET["ref"] if c.startswith(("cov", "s")) else None)
        if st["gamma_det_fe"] and gd_art:
            k = gd_art / st["gamma_det_fe"]
            st["escala"] = k
            st["escalado"] = estat(g * k)
        if c in ARTIGO_MC:
            a = ARTIGO_MC[c]
            st["artigo"] = {"mu": a[0], "sd": a[1], "cov": a[2] / 100, "pf": a[3] / 100, "covpf": a[4] / 100}
        if c in ARTIGO_PF_CHO:
            st["artigo"] = {"pf": ARTIGO_PF_CHO[c][0] / 100, "pf_cho": ARTIGO_PF_CHO[c][1] / 100}
        # convergência de Pf com N (ordem das amostras)
        st["conv_pf"] = [[int(n), float(np.mean(g[:n] < 1))] for n in
                         sorted(set(np.unique(np.logspace(0, np.log10(len(g)), 60).astype(int))))]
        if mec:
            km = [k for k in ks if k in mec]
            if len(km) > 20:
                G = np.array([float(rows[k]["fator"]) for k in km])
                C = np.array([float(mec[k]["c_banda"]) for k in km])
                P = np.array([float(mec[k]["phi_banda"]) for k in km])
                K = np.array([float(mec[k]["kv_massa"]) for k in km])
                Cd = np.array([float(rows[k]["c_medio"]) for k in km])
                Pd = np.array([float(rows[k]["phi_medio"]) for k in km])
                Kd = np.array([float(rows[k]["kv_medio"]) for k in km])
                sp = {"N": len(km)}
                if np.std(K) > 0:
                    sp["mecanismo"] = {"c": spearman_parcial(G, C, P, K), "phi": spearman_parcial(G, P, C, K),
                                       "kv": spearman_parcial(G, K, C, P)}
                    sp["dominio"] = {"c": spearman_parcial(G, Cd, Pd, Kd), "phi": spearman_parcial(G, Pd, Cd, Kd),
                                     "kv": spearman_parcial(G, Kd, Cd, Pd)}
                st["spearman"] = sp
        R["mc"][c] = st
        dados[c] = g
    # modos de ruptura (<bloco>.csv.modo)
    R["modos"] = {}
    for c in CASOS:
        cont = {}
        blocos = glob.glob(os.path.join(camp, c, "b*.csv.modo"))  # blocos; senão o arquivo juntado <caso>.modo
        for f in blocos or glob.glob(os.path.join(camp, c + ".modo")):
            for ln in open(f):
                p = ln.strip().split(",")
                if len(p) == 4 and p[0].isdigit():
                    cont[p[1]] = cont.get(p[1], 0) + 1
        if cont and c in ("cho_coesivo", "cho_cphi", "ref"):
            R["modos"][c] = cont
    R["achados"] = ACHADOS
    # ------------------------------------------------------------ determinísticos
    D = R["det"]
    D["conv"] = {os.path.basename(f)[:-4]: det_log(f) for f in sorted(glob.glob(os.path.join(det, "conv/*.log")))}
    D["tab56"] = {os.path.basename(f)[:-4]: det_log(f) for f in sorted(glob.glob(os.path.join(det, "tab56/*.log")))}
    D["fig9"] = {os.path.basename(f)[:-4]: det_log(f) for f in sorted(glob.glob(os.path.join(det, "fig9/*.log")))}
    D["fig8"] = {os.path.basename(f)[:-4]: det_log(f) for f in sorted(glob.glob(os.path.join(det, "fig8/*.log")))}
    D["fig5"] = {os.path.basename(f)[:-4]: det_log(f) for f in sorted(glob.glob(os.path.join(det, "fig5/*.log")))}
    for sub in ("dom", "conv8"):  # domínio × anisotropia; convergência de pontos da Fig. 8
        D[sub] = {os.path.basename(f)[:-4]: det_log(f) for f in sorted(glob.glob(os.path.join(det, sub, "*.log")))}
    comparacoes(R, ref)
    json.dump(R, open(os.path.join(out, "resultados.json"), "w"), indent=1, default=float)

    # ------------------------------------------------------------ markdown
    f3 = lambda v, d=3: "—" if v is None or (isinstance(v, float) and math.isnan(v)) else f"{v:.{d}f}"  # noqa: E731
    for titulo, casos in GRUPOS:
        md.append(f"### {titulo}\n")
        md.append("| caso | N | μ | σ | CoV % | Pf % (IC 95 %) | CoV(Pf) % | μ* | σ* | Pf* % | artigo μ | σ | CoV % | Pf % |")
        md.append("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
        for c in casos:
            st = R["mc"].get(c)
            if not st:
                continue
            a = st.get("artigo", {})
            e = st.get("escalado") or {}
            md.append(
                f"| {c} | {st['N']} | {f3(st['mu'])} | {f3(st['sd'])} | {f3(st['cov'] * 100, 1)} | "
                f"{f3(st['pf'] * 100, 2)} ({f3(st['pf_lo'] * 100, 1)}–{f3(st['pf_hi'] * 100, 1)}) | "
                f"{f3(st['covpf'] * 100, 1)} | {f3(e.get('mu'))} | {f3(e.get('sd'))} | "
                f"{f3(e['pf'] * 100 if e else None, 2)} | {f3(a.get('mu'))} | {f3(a.get('sd'))} | "
                f"{f3(a['cov'] * 100 if 'cov' in a else None, 1)} | "
                f"{f3(a['pf'] * 100 if 'pf' in a else None, 2)}"
                f"{' (Cho ' + f3(a['pf_cho'] * 100, 2) + ')' if 'pf_cho' in a else ''} |")
        md.append("")
    open(os.path.join(out, "tabelas.md"), "w").write("\n".join(md) + "\n")

    # ------------------------------------------------------------ figuras
    figuras(R, dados, det, ref, out)
    print("\n".join(md))


def salva(fig, out, nome):
    """PNG (relatório HTML) e PDF vetorial (relatório LaTeX)"""
    fig.savefig(os.path.join(out, nome + ".png"))
    fig.savefig(os.path.join(out, nome + ".pdf"))


def figuras(R, dados, det, ref, out):
    # figura do mecanismo (feita por mecanismo_artigo.py a partir dos VTK, que não ficam no repositório)
    for ext in ("pdf", "png"):
        src = os.path.join(det, "mech", "mecanismo_A35." + ext)
        if os.path.exists(src):
            shutil.copy(src, out)
    plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.alpha": 0.3, "figure.dpi": 130})
    C_FE, C_ART, C_ART2 = "#1f6fb4", "#c0392b", "#7f7f7f"

    # Fig. 5 — funcional hidráulico
    try:
        f5 = json.load(open(os.path.join(ref, "fig5.json")))
        fig, axs = plt.subplots(1, 4, figsize=(12, 3), sharex=True)
        for ax, a in zip(axs, (1, 2, 4, 10)):
            for s in f5["series"]:
                if s["alpha"] != a:
                    continue
                p = np.array(s["resampled_beta_5deg"], float)
                ls = "--" if "FE" in s["name"] else "-"
                nm = r"$J(u'_{FE})$" if "FE" in s["name"] else r"$-J^*(v'_{opt})$"
                ax.plot(p[:, 0], p[:, 1], ls, color=C_ART2, label=f"artigo {nm}")
            for pref, lab, mk in (("g_", "FE (crista/pé/base 50 m)", "o"), ("", "FE (crista/pé/base 10/10/5 m)", "s")):
                xs, ys = [], []
                for b in (15, 30, 45, 60, 75, 90):
                    v = R["det"]["fig5"].get(f"{pref}a{a}_b{b}", {}).get("J")
                    if v is not None:
                        xs.append(b)
                        ys.append(v)
                ax.plot(xs, ys, mk + "-", color=C_FE if pref else "#76b7e5", ms=4, label=lab)
            ax.set_title(rf"$\alpha = {a}$")
            ax.set_xlabel(r"$\beta$ (°)")
        axs[0].set_ylabel(r"$J/(k_h H^2 \gamma_w^2)$")
        # legenda abaixo dos painéis, para não cobrir os pontos
        h, lb = axs[0].get_legend_handles_labels()
        fig.legend(h, lb, loc="lower center", ncol=4, fontsize=8, frameon=False)
        fig.tight_layout(rect=(0, 0.08, 1, 1))
        salva(fig, out, "fig5_funcional")
        plt.close(fig)
    except Exception as e:  # noqa: BLE001
        print("fig5:", e)

    # Fig. 8 — altura crítica
    try:
        f8 = json.load(open(os.path.join(ref, "fig8.json")))
        fig, axs = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
        # painel do artigo -> (solo usado no cálculo do artigo, betas)
        paineis = (("London Clay", "B", (30, 60)), ("Israeli Clay", "A", (35, 60)))
        for ax, (pn, solo, betas) in zip(axs, paineis):
            for s in f8["series"]:
                if s["panel"] != pn:
                    continue
                p = np.array(s["points"], float)
                ls = {"solid": "-", "dashed": "--"}.get(s["style"].split()[0], "-.")
                ax.plot(p[:, 0], p[:, 1], ls, color=C_ART2, lw=1)
            for b in betas:
                xs, ys = [], []
                for hw in (0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5):
                    d = R["det"]["fig8"].get(f"{solo}_b{b}_hw{hw:g}", {})
                    if d.get("gamma") and not d.get("limite_maximo") and f"{solo}_b{b}_hw{hw:g}" not in NAO_CONV:
                        xs.append(hw / 5)
                        ys.append(5 * d["gamma"])
                ax.plot(xs, ys, "o-", color=C_FE, ms=4, label=rf"FE, $\beta = {b}°$")
                xs, ys = [], []
                other = "A" if solo == "B" else "B"
                for hw in (0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5):
                    d = R["det"]["fig8"].get(f"{other}_b{b}_hw{hw:g}", {})
                    if d.get("gamma") and not d.get("limite_maximo") and f"{other}_b{b}_hw{hw:g}" not in NAO_CONV:
                        xs.append(hw / 5)
                        ys.append(5 * d["gamma"])
                if xs:
                    ax.plot(xs, ys, "x:", color="#e67e22", ms=4, label=rf"FE, solo da Tab. 1 do artigo, $\beta = {b}°$")
            ax.set_yscale("log")
            ax.set_ylim(3, 300)
            ax.set_title(f"{pn} (painel do artigo)")
            ax.set_xlabel(r"$h_w/H$")
            ax.legend(fontsize=6)
        axs[0].set_ylabel(r"$H_{crit}$ (m)")
        fig.tight_layout()
        salva(fig, out, "fig8_hcrit")
        plt.close(fig)
    except Exception as e:  # noqa: BLE001
        print("fig8:", e)

    # Fig. 9 — Γ × β
    try:
        f9 = json.load(open(os.path.join(ref, "fig9.json")))
        fig, axs = plt.subplots(1, 3, figsize=(11, 3.3), sharey=True)
        for ax, a in zip(axs, (1, 5, 10)):
            for s in f9["series"]:
                if s["alpha"] != a:
                    continue
                p = np.array([(x, y) for x, y in s["points"] if y is not None], float)
                ls = "-" if "grad" in s["name"] else "--"
                ax.plot(p[:, 0], p[:, 1], ls, color=C_ART2,
                        label="artigo " + (r"$-\mathrm{grad}\,u'_{FE}$" if ls == "-" else r"$K\cdot v'_{opt}$"))
            xs, ys = [], []
            for b in (15, 30, 45, 60, 75, 90):
                d = R["det"]["fig9"].get(f"beta{b}_alfa{a}", {})
                if d.get("gamma"):
                    xs.append(b)
                    ys.append(d["gamma"])
            ax.plot(xs, ys, "o-", color=C_FE, ms=4, label="FE (nível 3)")
            ax.set_ylim(0, 5)
            ax.set_title(rf"$\alpha = {a}$")
            ax.set_xlabel(r"$\beta$ (°)")
        axs[0].set_ylabel(r"$\Gamma$")
        axs[0].legend(fontsize=7)
        fig.tight_layout()
        salva(fig, out, "fig9_gamma_beta")
        plt.close(fig)
    except Exception as e:  # noqa: BLE001
        print("fig9:", e)

    # densidades: referência (lognormal do artigo), Cho (curvas extraídas)
    try:
        from scipy.stats import gaussian_kde
        dist = json.load(open(os.path.join(ref, "dist.json")))
        fig, axs = plt.subplots(1, 3, figsize=(12, 3.3))
        x = np.linspace(0.01, 4, 400)
        if "ref" in dados:
            g = dados["ref"]
            ax = axs[0]
            ax.hist(g, bins=60, range=(0, 4), density=True, color="#cfe0f1", label=f"FE, N = {len(g)}")
            ax.plot(x, gaussian_kde(g)(x), color=C_FE, label="FE (KDE)")
            m, s = T_REF[0], T_REF[1]
            sl = math.sqrt(math.log(1 + (s / m) ** 2))
            ml = math.log(m) - sl * sl / 2
            ax.plot(x, np.exp(-(np.log(x) - ml) ** 2 / (2 * sl * sl)) / (x * sl * math.sqrt(2 * math.pi)), "--",
                    color=C_ART, label=r"artigo: lognormal ($\mu$ = 1,353, $\sigma$ = 0,318)")
            ax.axvline(1, color="k", lw=0.6)
            ax.set_title("Referência (Tabela 2 do artigo)")
            ax.set_xlabel(r"$\Gamma$")
            ax.legend(fontsize=6)
        for ax, c, fig_ in ((axs[1], "cho_coesivo", "Figure 17"), (axs[2], "cho_cphi", "Figure 21")):
            for s in dist["series"]:
                if s["figure"] == fig_ and s["curve_type"] == "PDF":
                    p = np.array(s["points"], float)
                    ax.plot(p[:, 0], p[:, 1], "--" if "Gamma" in s["name"] else ":",
                            color=C_ART if "Gamma" in s["name"] else C_ART2,
                            label=r"artigo, $\Gamma$ (análise limite)" if "Gamma" in s["name"] else r"Cho (2010), $F_s$")
            if c in dados:
                g = dados[c]
                xx = np.linspace(0.3, 6 if c == "cho_cphi" else 2.5, 400)
                ax.plot(xx, gaussian_kde(g)(xx), color=C_FE, label=rf"FE, $\Gamma$, N = {len(g)}")
            ax.set_title({"cho_coesivo": "Cho, coesivo", "cho_cphi": r"Cho, $c$-$\varphi$"}[c])
            ax.set_xlabel(r"$\Gamma$ ou $F_s$")
            ax.legend(fontsize=6)
        fig.tight_layout()
        salva(fig, out, "densidades")
        plt.close(fig)
    except Exception as e:  # noqa: BLE001
        print("densidades:", e)

    # tendências (Tabelas 3-6): Pf e μ
    try:
        series = [(r"CoV($c$)", ["covc10", "ref", "covc50", "covc70"], [0.1, 0.3, 0.5, 0.7]),
                  (r"CoV($\varphi$)", ["covphi5", "ref", "covphi15", "covphi20"], [0.05, 0.1, 0.15, 0.2]),
                  (r"CoV($k_v$)", ["covk0", "ref", "covk75", "covk90", "covk100"], [0, 0.6, 0.75, 0.9, 1.0]),
                  (r"escala $s$ (log)", ["ref", "s1.5", "s2", "s5", "s10", "s20", "s400"], [1, 1.5, 2, 5, 10, 20, 400]),
                  (r"$\alpha$", ["ref", "alfa2", "alfa3", "alfa4", "alfa5"], [1, 2, 3, 4, 5]),
                  (r"$h_w/H$", ["hw0.5", "hw0.6", "hw0.7", "hw0.8", "hw0.9", "ref"], [0.5, 0.6, 0.7, 0.8, 0.9, 1.0])]
        fig, axs = plt.subplots(3, 6, figsize=(16, 8))
        for j, (nome, casos, xv) in enumerate(series):
            for i, (key, lab) in enumerate((("pf", "Pf"), ("mu", r"$\mu(\Gamma)$"), ("cov", r"CoV($\Gamma$)"))):
                ax = axs[i, j]
                xa, ya, xf, yf, lo, hi, xs_, ys_ = [], [], [], [], [], [], [], []
                for c, xx in zip(casos, xv):
                    st = R["mc"].get(c)
                    if c in ARTIGO_MC:
                        a = ARTIGO_MC[c]
                        xa.append(xx)
                        ya.append({"pf": a[3] / 100, "mu": a[0], "cov": a[2] / 100}[key])
                    if st:
                        xf.append(xx)
                        yf.append(st[key])
                        if key == "pf":
                            lo.append(max(st["pf"] - st["pf_lo"], 0.))
                            hi.append(max(st["pf_hi"] - st["pf"], 0.))
                        if st.get("escalado"):
                            xs_.append(xx)
                            ys_.append(st["escalado"][key])
                ax.plot(xa, ya, "s--", color=C_ART, ms=4, label="artigo")
                if key == "pf" and xf:
                    ax.errorbar(xf, yf, yerr=[lo, hi], fmt="o-", color=C_FE, ms=4, capsize=2, label="FE (IC 95 %)")
                else:
                    ax.plot(xf, yf, "o-", color=C_FE, ms=4, label="FE")
                if xs_ and key != "cov":
                    ax.plot(xs_, ys_, "^:", color="#2e8b57", ms=4, label=r"FE $\times\,\Gamma_{art}/\Gamma_{FE}$")
                if "escala" in nome:
                    ax.set_xscale("log")
                if i == 2:
                    ax.set_xlabel(nome)
                if j == 0:
                    ax.set_ylabel(lab)
                if i == 0 and j == 0:
                    ax.legend(fontsize=6)
        fig.tight_layout()
        salva(fig, out, "tendencias")
        plt.close(fig)
    except Exception as e:  # noqa: BLE001
        print("tendencias:", e)

    # convergência de Pf (referência e Cho)
    try:
        fig, ax = plt.subplots(figsize=(6, 3.3))
        for c, col in (("ref", C_FE), ("cho_coesivo", "#e67e22"), ("cho_cphi", "#2e8b57")):
            st = R["mc"].get(c)
            if st:
                p = np.array(st["conv_pf"], float)
                ax.semilogx(p[:, 0], p[:, 1], color=col,
                            label={"ref": "referência", "cho_coesivo": "Cho, coesivo", "cho_cphi": r"Cho, $c$-$\varphi$"}[c])
        for v, col in ((0.115, C_FE), (0.065, "#e67e22"), (0.055, "#2e8b57")):
            ax.axhline(v, color=col, ls="--", lw=0.7)
        ax.set_xlabel("número de amostras")
        ax.set_ylabel("Pf")
        ax.set_ylim(0, 0.3)
        ax.legend(fontsize=7)
        fig.tight_layout()
        salva(fig, out, "convergencia_pf")
        plt.close(fig)
    except Exception as e:  # noqa: BLE001
        print("convergencia:", e)


if __name__ == "__main__":
    if len(sys.argv) != 5:
        print(__doc__)
        sys.exit(1)
    main(*sys.argv[1:])
