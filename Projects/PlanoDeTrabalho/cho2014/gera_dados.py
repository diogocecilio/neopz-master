#!/usr/bin/env python3
"""gera_dados.py: as tabelas, os numeros e as figuras do relatorio_cho2014.tex (a coluna 1D de Cho 2014, Eng. Geol.
171: 11-20), a partir de:

  - <campanha> (~/projects/runs/fermi_cho2014): o MC do NeoPZ (SlopeRainColumn, quad ny 80) dos 15 casos da Fig. 10
    (5000 cada) e dos estudos ny200, tri e tricons (2000) e eps1e-4 (5000, 30 % parados: descartado);
  - ~/projects/runs/cho2014_python: o MC de referencia em Python (coluna_ref.py, FV N = 200), 5000 por caso e 55 000 no
    caso base;
  - ~/projects/runs/coluna_estudos (V3.2: k* do limite de variavel aleatoria, quad ny 80) e
    ~/projects/runs/coluna_ls_piloto (line sampling: as linhas, a busca=cerca de 24 h e as raizes C++ x Python);
  - ~/projects/runs/cho2014_python/det: o deterministico em Python (N = 200);
  - ~/projects/runs/coluna_verificacao_independente/det/q80: o deterministico quad ny 80 (o padrao de producao) de 0 a
    24 h, com os perfis; o 16 h vem do deterministico da campanha;
  - ~/projects/runs/coluna_estudos/v32 e v33: FS_min(T) contra k_s/q na coluna homogenea (bissecao e grade);
  - dados/: os recortes das figuras do PDF de Cho (pdftoppm, 200 dpi).
Nenhuma simulacao nova: so le o que ja existe.

Os valores do artigo sao os de cho2014.md (research_notes/Reproducao dos artigos novos): a Fig. 10 digitalizada
(+-0,002) e as Figs. 3 e 5 a 9 lidas a olho.

Grava tab/*.tex, numeros.tex e figs/*.pdf.

    python3 -I gera_dados.py [campanha = ~/projects/runs/fermi_cho2014]"""

import csv
import glob
import json
import math
import os
import subprocess
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from scipy import stats  # noqa: E402

AQUI = os.path.dirname(os.path.abspath(__file__))
RAIZ = os.path.abspath(os.path.join(AQUI, "..", "..", ".."))
sys.path.insert(0, os.path.join(RAIZ, "Projects", "SlopeRandomFieldRefinement", "cluster"))
from resumo_hu2024 import wilson  # noqa: E402

CAMPANHA = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/projects/runs/fermi_cho2014")
RUNS = os.path.expanduser("~/projects/runs")
PY = os.path.join(RUNS, "cho2014_python")
EST = os.path.join(RUNS, "coluna_estudos")
LSP = os.path.join(RUNS, "coluna_ls_piloto")
DET = os.path.join(RUNS, "coluna_verificacao_independente", "det")  # q80: quad ny 80, o padrao de producao
PDF = os.path.expanduser("~/Dropbox/UnsaturatedSoil/referencias/referencias_selecionadas_claude/novas/Cho_2014.pdf")
DADOS = os.path.join(AQUI, "dados")
FIGS = os.path.join(AQUI, "figs")
TAB = os.path.join(AQUI, "tab")
TEMPOS = [6., 12., 16., 18., 20., 22., 24.]
Q = 3e-6  # a chuva e o k_s medio (m/s)

# ------------------------------------------------------------------------------------------------------- o artigo
COVS = [0.5, 1.0, 1.5]
LVS = [0.5, 2., 4., 8., 16.]
# Fig. 10 digitalizada (cho2014.md, Secao 4.3): (CoV, lv) -> (P_f 18 h, P_f 24 h); o texto da 0,434 no caso base
FIG10 = {(0.5, 0.5): (0.058, 0.435), (1.0, 0.5): (0.091, 0.275), (1.5, 0.5): (0.103, 0.218),
         (0.5, 2.): (0.081, 0.501), (1.0, 2.): (0.062, 0.351), (1.5, 2.): (0.066, 0.282),
         (0.5, 4.): (0.091, 0.520), (1.0, 4.): (0.053, 0.371), (1.5, 4.): (0.049, 0.302),
         (0.5, 8.): (0.103, 0.533), (1.0, 8.): (0.061, 0.391), (1.5, 8.): (0.049, 0.321),
         (0.5, 16.): (0.113, 0.543), (1.0, 16.): (0.070, 0.402), (1.5, 16.): (0.051, 0.331)}
# Fig. 9: a fracao de superficies criticas na base em 12/18/24 h
FIG9 = {"c05_l0.5": (0.29, 0.10, 0.44), "c10_l0.5": (0.555, 0.36, 0.43), "c15_l0.5": (0.665, 0.51, 0.52),
        "c05_l2": (0.255, 0.13, 0.51), "c05_l8": (0.235, 0.15, 0.54)}
FIG7_MED = {12: 1.08, 18: 1.38, 24: 1.77}   # medianas de z_crit, caso base
FIG5_MODA = {12: (1.05, 1.10, 0.10), 18: (1.65, 1.70, 0.063)}  # classe modal acima da base e a sua frequencia
FIG6 = {"media": 0.8376, "dp": 0.3706, "min": 0.3626, "max": 1.4806, "pf": 0.434}
FIG6_PICOS = [(0.36, 0.39, 0.29), (1.05, 1.08, 0.10), (1.45, 1.48, 0.016)]
# Fig. 3, lida: FS minimo (faixa) e a frente (faixa) por instante
FIG3_FS = {0: (1.49, 1.51), 6: (1.45, 1.49), 12: (1.25, 1.27), 18: (1.05, 1.07), 24: (0.36, 0.37)}
FIG3_FRENTE = {6: (0.60, 0.65), 12: (1.15, 1.25), 18: (1.75, 1.85)}
N_CHO = 5000

# --------------------------------------------------------------------------------------------------------- estilo
PALETA = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
TINTA, TINTA2, GRADE = "#0b0b0b", "#52514e", "#dddcd8"
plt.rcParams.update({
    "font.size": 8.5, "axes.labelsize": 8.5, "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "legend.fontsize": 7,
    "axes.edgecolor": TINTA2, "axes.labelcolor": TINTA, "xtick.color": TINTA2, "ytick.color": TINTA2,
    "axes.grid": True, "grid.color": GRADE, "grid.linewidth": 0.5, "lines.linewidth": 1.4, "legend.frameon": False,
    "savefig.dpi": 200, "font.family": "DejaVu Sans", "mathtext.fontset": "dejavusans",
})
COR_COV = {0.5: PALETA[7], 1.0: PALETA[0], 1.5: PALETA[2]}   # as cores de Cho: vermelho, azul, verde
MARC_COV = {0.5: "o", 1.0: "D", 1.5: "^"}
COR_T = {0: TINTA2, 6: PALETA[7], 12: PALETA[0], 16: PALETA[6], 18: PALETA[2], 20: PALETA[3], 22: PALETA[4], 24: TINTA}
CHO = TINTA2


def br(x, casas=2):
    if x is None or (isinstance(x, float) and not math.isfinite(x)):
        return "--"
    t = f"{x:.{casas}f}".replace(".", "{,}")
    return r"\ensuremath{-}" + t[1:] if t.startswith("-") else t


def pct(x, casas=1):
    return br(100 * x, casas)


def sinal(x, casas=1):
    t = br(abs(x), casas)
    return ("+" if x >= 0 else r"\ensuremath{-}") + t


def inteiro(n):
    return f"{int(n):,}".replace(",", r"\,")


def painel(ax, t):
    ax.text(-0.02, 1.02, "(%s)" % t, transform=ax.transAxes, ha="right", va="bottom", fontsize=9, color=TINTA)


def salva(fig, nome):
    fig.savefig(os.path.join(FIGS, nome), metadata={"CreationDate": None, "ModDate": None})
    plt.close(fig)


def cid(cov, lv):
    return "c%02d_l%g" % (round(10 * cov), lv)


def pf_ln(kstar, cov):
    """P(k_s > k*) com k_s lognormal de media q e CoV cov (k* em unidades de q)"""
    s = math.sqrt(math.log(1 + cov * cov))
    return float(stats.norm.sf((math.log(kstar) + 0.5 * s * s) / s))


# ------------------------------------------------------------------------------------------------------------ dados
SIT = {"ok": 0, "falha": 1, "transiente": 2}


def le_mc(arqs):
    """as realizacoes de varios csv (SlopeRainColumn ou coluna_ref.py) em matrizes (n, len(TEMPOS))"""
    reais, linhas = set(), []
    for a in arqs:
        with open(a, newline="") as f:
            r = csv.reader(f)
            cab = next(r, None)
            if not cab:
                continue
            i = {k: cab.index(k) for k in ("real", "t_h", "FS_min", "situacao", "z_crit", "base", "FS_passagem",
                                           "t_passagem", "KL_variancia")}
            for v in r:
                if len(v) < len(cab) or not v[i["KL_variancia"]]:
                    continue
                t = float(v[i["t_h"]])
                if t not in TEMPOS:
                    continue
                n = int(v[i["real"]])
                reais.add(n)
                linhas.append((n, TEMPOS.index(t), float(v[i["FS_min"]]), SIT.get(v[i["situacao"]], 2),
                               float(v[i["z_crit"]]), float(v[i["base"]]), float(v[i["FS_passagem"]]),
                               float(v[i["t_passagem"]])))
    reais = np.array(sorted(reais))
    pos = {n: k for k, n in enumerate(reais)}
    nt = len(TEMPOS)
    d = {"real": reais}
    for k in ("fs", "z", "base", "fsp", "tp"):
        d[k] = np.full((len(reais), nt), np.nan)
    d["sit"] = np.full((len(reais), nt), -1, int)
    for n, j, fs, s, z, b, fsp, tp in linhas:
        k = pos[n]
        d["fs"][k, j], d["sit"][k, j], d["z"][k, j], d["base"][k, j], d["fsp"][k, j], d["tp"][k, j] = fs, s, z, b, fsp, tp
    return d


def falha(d, t):
    """1 falha, 0 segura, nan censurada (o transiente parou antes, sem travessia) no instante t"""
    j = TEMPOS.index(t)
    s, fs, tp = d["sit"][:, j], d["fs"][:, j], d["tp"][:, j]
    out = np.where(s == 1, 1., np.where(fs < 1., 1., 0.))
    cens = (s == 2) | (s == -1)
    out[cens] = np.where((tp[cens] >= 0) & (tp[cens] <= t), 1., np.nan)
    return out


def pf(d, t, sel=None):
    f = falha(d, t)
    if sel is not None:
        f = f[sel]
    ok = np.isfinite(f)
    n, k = int(ok.sum()), int(np.nansum(f))
    p = k / n
    return {"P": p, "k": k, "n": n, "SE": math.sqrt(p * (1 - p) / n), "IC": wilson(k, n), "cens": int((~ok).sum())}


def le_caso_neopz(c):
    return le_mc(sorted(glob.glob(os.path.join(CAMPANHA, c, "mc", "mc_*.csv"))))


def le_caso_py(c):
    return le_mc([a for a in sorted(glob.glob(os.path.join(PY, c, "*.csv")))
                  if "_parametros" not in a and "_perfil_" not in a])


def le_perfil(arq):
    out = {}
    for r in csv.DictReader(open(arq)):
        t = float(r["t_h"])
        out.setdefault(t, []).append((float(r["z"]), float(r["pw"]), float(r["FS"])))
    return {t: np.array(v) for t, v in out.items()}


def le_det_csv(arq):
    return {float(r["t_h"]): r for r in csv.DictReader(open(arq))}


def recortes():
    """as figuras do artigo, recortadas do PDF a 200 dpi (pagina, x, y, largura, altura em pixels)"""
    caixas = {"cho_fig3": (5, 851, 137, 714, 1229), "cho_fig5": (6, 851, 137, 714, 1449),
              "cho_fig6": (7, 137, 137, 686, 506), "cho_fig7": (7, 137, 1457, 686, 548),
              "cho_fig8": (7, 857, 137, 709, 450), "cho_fig9": (8, 114, 137, 714, 1795),
              "cho_fig10": (8, 857, 137, 709, 1062)}
    for nome, (p, x, y, w, h) in caixas.items():
        if not os.path.exists(os.path.join(DADOS, nome + ".png")) and os.path.exists(PDF):
            subprocess.run(["pdftoppm", "-r", "200", "-f", str(p), "-l", str(p), "-x", str(x), "-y", str(y), "-W", str(w),
                            "-H", str(h), "-png", "-singlefile", PDF, os.path.join(DADOS, nome)], check=True)


def kgrade():
    """FS_min(T) contra k_s/q na coluna homogenea, dos pontos ja rodados: a bissecao do V3.2 (quad ny 80, T = 18, 20 e
    24 h) e a grade do V3.3 (quad ny 200, dt 0,1, T = 18 h)"""
    out = {}
    for T in (18, 20, 24):
        v = [(float(x["kmult"]), float(x["FS_T"])) for x in
             csv.DictReader(open(os.path.join(EST, "v32", "quad_ny80_dt0.1", f"T{T}_kstar.csv")))]
        if T == 18:
            for a in glob.glob(os.path.join(EST, "v33", "quad_ny200_dt0.1", "k_k*_0_-1.csv")):
                k = float(os.path.basename(a)[3:].split("_0_-1")[0])
                v.append((k, float(le_det_csv(a)[18.]["FS_min"])))
        out[T] = sorted(v)
    return out


def kstar_v32():
    """k*/q do V3.2 (quad ny 80) em 18, 20 e 24 h"""
    out = {}
    for T in (18, 20, 24):
        txt = open(os.path.join(EST, "v32", "quad_ny80_dt0.1", f"T{T}_kstar.txt")).read()
        out[T] = float(txt.split("k* = ")[1].split(" q")[0])
    return out


def ls_linhas(pasta, ini=None, fim=None):
    p = [(int(r["linha"]), float(r["p"]), int(r["analises"])) for r in csv.DictReader(open(os.path.join(pasta, "linhas.csv")))
         if r.get("p") not in (None, "", "nan")]
    if ini is not None:
        p = [x for x in p if ini <= x[0] <= fim]
    v = np.array([x[1] for x in p])
    n, m = len(v), v.mean()
    sd = v.std(ddof=1)
    return {"P": m, "SE": sd / math.sqrt(n), "CoV": sd / m / math.sqrt(n), "CoVl": sd / m, "n": n,
            "an": float(np.mean([x[2] for x in p]))}


def raizes_cpp_py():
    """FS_passagem(18 h) do C++ (quad ny 80) e do Python (N 200) nos pontos de raiz das 16 linhas (k 0) e a +-1 em c"""
    sr = os.path.join(LSP, "segunda_raiz")
    out = []
    for r in csv.DictReader(open(os.path.join(sr, "varredura.csv"))):
        if r["tipo"] != "grade" or r["k"] not in ("0", "4", "-4"):
            continue
        rot = f"{r['corrida']}_l{r['linha']}_k{r['k']}"
        fs = sorted(a for a in glob.glob(os.path.join(sr, "py_raizes", rot + "_*.csv")) if "perfil" not in a)
        if not fs:
            continue
        L = [x for x in csv.DictReader(open(fs[0])) if abs(float(x["t_h"]) - 18) < 1e-9][0]
        out.append((r["corrida"], int(r["linha"]), int(r["k"]), float(r["c"]) - float(r["c_i"]),
                    float(r["FS_passagem"]), float(L["FS_passagem"])))
    return out


# ------------------------------------------------------------------------------------------------------------ contas
def carrega():
    recortes()
    r = {"casos": {}, "py": {}, "est": {}}
    for cov in COVS:
        for lv in LVS:
            c = cid(cov, lv)
            r["casos"][c] = le_caso_neopz(c)
            r["py"][c] = le_caso_py(c)
    for c in ("c05_l0.5_ny200", "c05_l16_ny200", "c05_l0.5_tri", "c05_l0.5_tricons", "c05_l0.5_eps1e-4"):
        r["est"][c] = le_caso_neopz(c)
    r["det_cpp"] = le_perfil(os.path.join(DET, "q80_perfil_-1.csv"))
    r["det_csv"] = le_det_csv(os.path.join(DET, "q80_0_-1.csv"))
    # o 16 h do deterministico da campanha (o mesmo modelo) completa a serie FS(t)
    for t, v in le_det_csv(os.path.join(CAMPANHA, "c05_l0.5", "det", "det_0_-1.csv")).items():
        r["det_csv"].setdefault(t, v)
    r["det_py"] = le_perfil(os.path.join(PY, "det", "det_N200_dt0.1_perfil_-1.csv"))
    r["det_py_csv"] = le_det_csv(os.path.join(PY, "det", "det_N200_dt0.1_0_-1.csv"))
    r["kgrade"] = kgrade()
    r["kstar"] = kstar_v32()
    r["resumo"] = json.load(open(os.path.join(DADOS, "resumo.json")))
    return r


def contas(r):
    a = {}
    # P_f por caso e instante
    a["pf"], a["pfpy"] = {}, {}
    for c, d in r["casos"].items():
        a["pf"][c] = {t: pf(d, t) for t in TEMPOS}
        a["pfpy"][c] = {t: pf(r["py"][c], t) for t in TEMPOS}
    # conferencia contra o resumo_cho2014.py
    for e in r["resumo"]["casos"]:
        if e["id"] in a["pf"]:
            for t in ("18", "24"):
                assert abs(e["instante"][t]["P_f"] - a["pf"][e["id"]][float(t)]["P"]) < 1e-12, (e["id"], t)
    # o instante em que a nossa P_f iguala a de Cho em 18 h (interpolacao linear entre os instantes gravados)
    a["t_equiv"] = {}
    for (cov, lv), (p18, _) in FIG10.items():
        c = cid(cov, lv)
        ps = [a["pf"][c][t]["P"] for t in TEMPOS]
        a["t_equiv"][c] = float(np.interp(p18, ps, TEMPOS)) if max(ps) > p18 else None
    # z_crit e base
    a["alvo"] = {}
    for c, d in r["casos"].items():
        x = {}
        for t in (12., 18., 24.):
            j = TEMPOS.index(t)
            ok = d["sit"][:, j] == 0
            z = d["z"][ok, j]
            x[t] = {"base": float(np.mean(d["base"][ok, j])), "med": float(np.median(z)), "z": z,
                    "acima075": float(np.mean(z < 0.75))}
        j = TEMPOS.index(24.)
        fs = d["fs"][d["sit"][:, j] == 0, j]
        x["fs24"] = {"media": fs.mean(), "dp": fs.std(ddof=1), "min": fs.min(), "max": fs.max(), "fs": fs,
                     "pico": float(np.mean((fs >= 0.36) & (fs < 0.39)))}
        a["alvo"][c] = x
    # os estudos, pareados com as realizacoes de mesmo numero do caso de producao
    a["est"] = {}
    for c, d in r["est"].items():
        prod = r["casos"]["c05_l16" if "l16" in c else "c05_l0.5"]
        comuns = np.intersect1d(d["real"], prod["real"])
        sa, sb = np.isin(d["real"], comuns), np.isin(prod["real"], comuns)
        e = {"n": len(d["real"]), "pararam": int(np.any(d["sit"] == 2, axis=1).sum())}
        for t in (18., 24.):
            fa, fb = falha(d, t)[sa], falha(prod, t)[sb]
            ok = np.isfinite(fa) & np.isfinite(fb)
            e[t] = {"P": float(np.nanmean(fa)), "Pprod": float(np.mean(fb[ok])), "Pa_ok": float(np.mean(fa[ok])),
                    "b": int(((fa == 1) & (fb == 0) & ok).sum()), "c": int(((fa == 0) & (fb == 1) & ok).sum()),
                    "n_ok": int(ok.sum()), "cens": int((~np.isfinite(fa)).sum()),
                    "lim": (float(np.nansum(fa) / len(fa)), float((np.nansum(fa) + np.isnan(fa).sum()) / len(fa)))}
        j = TEMPOS.index(12.)
        e["base12"] = float(np.mean(d["base"][d["sit"][:, j] == 0, j]))
        a["est"][c] = e
    # o limite de variavel aleatoria
    a["lim"] = {T: {cov: pf_ln(k, cov) for cov in COVS} for T, k in r["kstar"].items()}
    # o line sampling
    a["ls"] = {
        "rec": ls_linhas(os.path.join(LSP, "ls_p5000_T18_rec")),
        "rec100": ls_linhas(os.path.join(LSP, "ls_p5000_T18_rec"), 0, 99),
        "p500": ls_linhas(os.path.join(LSP, "ls_p500_T18")),
        "tr10": ls_linhas(os.path.join(LSP, "ls_p5000tr10_T18")),
        "tr5": ls_linhas(os.path.join(LSP, "ls_p5000tr5_T18")),
        "est2": ls_linhas(os.path.join(LSP, "ls_est2_50a199_T18")),
        "T24": ls_linhas(os.path.join(LSP, "ls_p5000_T24")),
        "T24c": ls_linhas(os.path.join(LSP, "ls_p5000_T24_cerca")),
    }
    rz = json.load(open(os.path.join(LSP, "diagnostico", "raizes_T24.json")))
    a["ls"]["T24corr_P"] = float(np.mean([x.get("p_corr", x.get("dp", 0) + 0) for x in rz])) if rz and "p_corr" in rz[0] \
        else None
    # T1 (tr10) e T2 (est2) tem sementes independentes: pelo inverso da variancia
    t1, t2 = a["ls"]["tr10"], a["ls"]["est2"]
    w1, w2 = 1 / t1["SE"] ** 2, 1 / t2["SE"] ** 2
    a["ls"]["pool"] = {"P": (w1 * t1["P"] + w2 * t2["P"]) / (w1 + w2), "SE": 1 / math.sqrt(w1 + w2)}
    # as raizes C++ x Python e o deslocamento da P_f do MC em Python (55 000) por um FS deslocado
    a["raizes"] = raizes_cpp_py()
    d0 = np.array([x[4] - x[5] for x in a["raizes"] if x[2] == 0])
    dall = np.array([x[4] - x[5] for x in a["raizes"]])
    a["dfs0"], a["dfs0_se"], a["dfsall"] = float(d0.mean()), float(d0.std(ddof=1) / math.sqrt(len(d0))), float(dall.mean())
    fsp = r["py"]["c05_l0.5"]["fs"][:, TEMPOS.index(18.)]
    a["py_fs18"] = fsp
    a["desloc"] = lambda dd: float(np.mean(fsp < 1 + dd) - np.mean(fsp < 1))
    a["mc18"] = {"py": a["pfpy"]["c05_l0.5"][18.], "py5000": pf(r["py"]["c05_l0.5"], 18., r["py"]["c05_l0.5"]["real"] < 5000),
                 "cpp": a["pf"]["c05_l0.5"][18.]}
    a["mc24"] = {"py": a["pfpy"]["c05_l0.5"][24.], "py5000": pf(r["py"]["c05_l0.5"], 24., r["py"]["c05_l0.5"]["real"] < 5000),
                 "cpp": a["pf"]["c05_l0.5"][24.]}
    return a


# ----------------------------------------------------------------------------------------------------------- figuras
def fig3(r):
    P, Ppy = r["det_cpp"], r["det_py"]
    fig = plt.figure(figsize=(6.6, 6.4))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.25, 1], hspace=0.33, wspace=0.28)
    ax1, ax2, ax3 = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[1, :])
    for t in (0, 6, 12, 18, 24):
        v = P[float(t)]
        ax1.plot(v[:, 1] / 9.81, v[:, 0], color=COR_T[t], label=f"{t} h")
        fs = np.where(np.isfinite(v[:, 2]), v[:, 2], np.nan)
        ax2.plot(fs, v[:, 0], color=COR_T[t], label=f"{t} h")
        if float(t) in Ppy:
            w = Ppy[float(t)]
            ax1.plot(w[:, 1] / 9.81, w[:, 0], color=COR_T[t], ls=":", lw=1.0)
            ax2.plot(w[:, 2], w[:, 0], color=COR_T[t], ls=":", lw=1.0)
    for t, (z0, z1) in FIG3_FRENTE.items():
        ax1.plot([0.35, 0.35], [z0, z1], color=COR_T[t], lw=4, alpha=0.45, solid_capstyle="butt")
    ax1.plot([], [], color=TINTA2, lw=4, alpha=0.45, label="frente (Cho, lida)")
    ax1.plot([], [], color=TINTA2, ls=":", lw=1, label="Python N 200")
    ax1.set(xlim=(-3, 3), ylim=(2, 0), xlabel="carga de pressão (m)", ylabel="profundidade (m)")
    ax1.legend(loc="lower left", fontsize=6.3)
    ax2.set(xlim=(0, 10), ylim=(2, 0), xlabel="FS", ylabel="profundidade (m)")
    # o FS minimo de Cho (faixa lida) na profundidade critica
    zc = {0: 2., 6: 2., 12: 1.15, 18: 1.72, 24: 2.}
    for t, (f0, f1) in FIG3_FS.items():
        ax2.errorbar(0.5 * (f0 + f1), zc[t], xerr=0.5 * (f1 - f0) + 0.03, fmt="s", ms=4, mfc="none", color=COR_T[t],
                     capsize=2, lw=0.9)
    ax2.plot([], [], "s", mfc="none", color=TINTA2, label="FS mín. (Cho, lido)")
    ax2.axvline(1, color=TINTA2, lw=0.6, ls="--")
    ax2.legend(loc="lower right", fontsize=6.3)
    painel(ax1, "a")
    painel(ax2, "b")
    # FS_min(t)
    D = r["det_csv"]
    ts = sorted(D)
    ax3.plot(ts, [float(D[t]["FS_min"]) for t in ts], color=PALETA[0], label="NeoPZ (quad ny 80)")
    Dp = r["det_py_csv"]
    tp = sorted(Dp)
    ax3.plot(tp, [float(Dp[t]["FS_min"]) for t in tp], "o", ms=3.5, mfc="none", color=PALETA[1], label="Python N 200")
    for t, (f0, f1) in FIG3_FS.items():
        ax3.errorbar(t, 0.5 * (f0 + f1), yerr=0.5 * (f1 - f0) + 0.005, fmt="s", ms=4, color=TINTA, capsize=2, lw=0.9,
                     label="Cho (Fig. 3b, lido)" if t == 0 else None)
    tpas = float(D[24.]["t_passagem"])
    ax3.axvline(tpas, color=TINTA2, lw=0.6, ls="--")
    ax3.text(tpas - 0.2, 1.3, f"a frente chega à base:\n$t$ = {tpas:.2f} h".replace(".", ","), ha="right", fontsize=7,
             color=TINTA2)
    ax3.axhline(1, color=TINTA2, lw=0.6)
    ax3.set(xlim=(0, 24), ylim=(0.3, 1.6), xticks=range(0, 25, 3), xlabel="tempo (h)", ylabel="FS mínimo na coluna")
    ax3.legend(loc="lower left")
    painel(ax3, "c")
    salva(fig, "fig3_det.pdf")


def hist_z(ax, z, cor, lw=0, preenche=True, rot=None):
    """o histograma horizontal de z_crit como o da Fig. 5: classes de 0,05 m acima da base e a barra da base"""
    n = len(z)
    base = z >= 1.9999
    bordas = np.arange(0., 2.0001, 0.05)
    h, _ = np.histogram(z[~base], bins=bordas)
    h = h / n
    if preenche:
        ax.barh(bordas[:-1] + 0.025, h, height=0.05, color=cor, alpha=0.55, label=rot)
        ax.barh([2.0], [base.mean()], height=0.035, color=cor, alpha=0.9)
    else:
        ax.step(np.r_[h, h[-1]], bordas, where="post", color=cor, lw=0.9, label=rot)
        ax.plot([0, base.mean()], [2.0, 2.0], color=cor, lw=2.2)


def fig5(r, a):
    fig, axs = plt.subplots(1, 3, figsize=(7.0, 3.0), sharey=True)
    for ax, t, L in zip(axs, (12., 18., 24.), "abc"):
        hist_z(ax, a["alvo"]["c05_l0.5"][t]["z"], PALETA[0], rot="NeoPZ")
        zp = r["py"]["c05_l0.5"]["z"][:, TEMPOS.index(t)]
        hist_z(ax, zp[np.isfinite(zp)], PALETA[1], preenche=False, rot="Python (55 000)")
        b = FIG9["c05_l0.5"][(12., 18., 24.).index(t)]
        ax.plot(b, 2.0, "s", ms=6, mfc="none", mec=TINTA, mew=1.2, label="Cho: base")
        if int(t) in FIG5_MODA:
            z0, z1, f = FIG5_MODA[int(t)]
            ax.plot(f, 0.5 * (z0 + z1), "D", ms=5, mfc="none", mec=TINTA, mew=1.2, label="Cho: classe modal")
        ax.set(xlim=(0, 0.5), ylim=(2.05, 0), xlabel="frequência relativa", title=f"{t:g} h")
        painel(ax, L)
    axs[0].set_ylabel("profundidade crítica (m)")
    axs[0].legend(loc="upper right", fontsize=6.3)
    fig.tight_layout()
    salva(fig, "fig5_zcrit.pdf")


def fig6(r, a):
    x = a["alvo"]["c05_l0.5"]["fs24"]
    fs = x["fs"]
    fsp = r["py"]["c05_l0.5"]["fs"][:, TEMPOS.index(24.)]
    fsp = fsp[np.isfinite(fsp)]
    bordas = np.arange(0.0, 1.62, 0.03)
    fig, ax = plt.subplots(figsize=(6.4, 3.0))
    h, _ = np.histogram(fs, bins=bordas)
    ax.bar(bordas[:-1] + 0.015, h / len(fs), width=0.03, color=PALETA[0], alpha=0.6, label="NeoPZ (5000)")
    hp, _ = np.histogram(fsp, bins=bordas)
    ax.step(bordas, np.r_[hp / len(fsp), 0], where="post", color=PALETA[1], lw=0.9, label="Python (55 000)")
    for i, (f0, f1, v) in enumerate(FIG6_PICOS):
        ax.plot(0.5 * (f0 + f1), v, "s", ms=6, mfc="none", mec=TINTA, mew=1.2, label="Cho: picos lidos" if i == 0 else None)
    ax.axvline(1, color=TINTA2, lw=0.6, ls="--")
    txt = ("              NeoPZ    Cho\n"
           f"P(FS<1)   {a['pf']['c05_l0.5'][24.]['P']:.3f}    {FIG6['pf']:.3f}\n"
           f"média      {x['media']:.3f}    {FIG6['media']:.3f}\n"
           f"desvio      {x['dp']:.3f}    {FIG6['dp']:.3f}\n"
           f"mín.        {x['min']:.3f}    {FIG6['min']:.3f}\n"
           f"máx.        {x['max']:.3f}    {FIG6['max']:.3f}").replace(".", ",")
    ax.text(0.47, 0.95, txt, transform=ax.transAxes, va="top", fontsize=6.8, family="DejaVu Sans Mono")
    ax.set(xlim=(0.2, 1.6), ylim=(0, 0.4), xlabel="FS em 24 h", ylabel="frequência relativa")
    ax.legend(loc="upper left")
    fig.tight_layout()
    salva(fig, "fig6_fs24.pdf")


def cdf_prof(z, d):
    """P(z_crit >= d): a curva das Figs. 7 e 9 (vale a fracao na base em d = 2 m e sobe para as rasas)"""
    z = np.sort(z)
    return 1. - np.searchsorted(z, d - 1e-9, side="left") / len(z)


def fig7(r, a):
    d = np.linspace(0, 2.0, 801)
    fig, ax = plt.subplots(figsize=(4.4, 3.2))
    for t, ls in ((12., "--"), (18., "-."), (24., "-")):
        z = a["alvo"]["c05_l0.5"][t]["z"]
        ax.plot(d, cdf_prof(z, d), color=COR_T[int(t)], ls=ls, drawstyle="steps-post", label=f"NeoPZ {t:g} h")
        ax.plot(2.0, FIG9["c05_l0.5"][(12., 18., 24.).index(t)], "s", ms=6, mfc="none", mec=COR_T[int(t)], mew=1.2)
        ax.plot(FIG7_MED[int(t)], 0.5, "D", ms=5, mfc="none", mec=COR_T[int(t)], mew=1.2)
    ax.plot([], [], "s", mfc="none", mec=TINTA, label="Cho: na base")
    ax.plot([], [], "D", mfc="none", mec=TINTA, label="Cho: mediana")
    ax.axvline(0.75, color=TINTA2, lw=0.6, ls=":")
    ax.set(xlim=(2.0, 0.0), ylim=(0, 1), xlabel="profundidade da superfície crítica (m)", ylabel="probabilidade acumulada")
    ax.legend(loc="center right", fontsize=6.5)
    fig.tight_layout()
    salva(fig, "fig7_cdf.pdf")


def fig8(r, a):
    fig, axs = plt.subplots(1, 2, figsize=(7.0, 2.8))
    d = r["casos"]["c05_l0.5"]
    N = np.arange(1, len(d["real"]) + 1)
    for t in (6., 12., 18., 24.):
        f = falha(d, t)
        axs[0].plot(N, np.cumsum(f) / N, color=COR_T[int(t)], lw=1.0, label=f"{t:g} h")
    axs[0].plot([5000], [0.058], "^", ms=6, mfc="none", mec=COR_T[18], mew=1.2, label="Cho 18 h (lido)")
    axs[0].plot([5000], [0.434], "x", ms=6, color=TINTA, mew=1.2, label="Cho 24 h (texto)")
    axs[0].set(xscale="log", xlim=(1, 1e4), ylim=(0, 1), xlabel="número de simulações",
               ylabel="probabilidade de falha", title="NeoPZ, caso base (5000)")
    axs[0].legend(loc="upper right", fontsize=6.3, ncol=2)
    painel(axs[0], "a")
    dp = r["py"]["c05_l0.5"]
    Np = np.arange(1, len(dp["real"]) + 1)
    for t in (18., 24.):
        f = falha(dp, t)
        m = np.cumsum(f) / Np
        se = np.sqrt(m * (1 - m) / Np)
        axs[1].plot(Np, m, color=COR_T[int(t)], lw=1.0, label=f"Python {t:g} h")
        axs[1].fill_between(Np, m - 1.96 * se, m + 1.96 * se, color=COR_T[int(t)], alpha=0.15, lw=0)
    axs[1].axhline(0.058, color=COR_T[18], ls="--", lw=0.8, label="Cho 18 h")
    axs[1].axhline(0.434, color=TINTA, ls="--", lw=0.8, label="Cho 24 h")
    axs[1].set(xscale="log", xlim=(10, 6e4), yscale="log", ylim=(3e-3, 1), xlabel="número de simulações",
               title="Python, caso base (55 000), IC 95 %")
    axs[1].legend(loc="lower right", fontsize=6.3, ncol=2)
    painel(axs[1], "b")
    fig.tight_layout()
    salva(fig, "fig8_conv.pdf")


def fig9(r, a):
    casos = [("c05_l0.5", "--", PALETA[7]), ("c10_l0.5", "-.", TINTA), ("c15_l0.5", "-", PALETA[7]),
             ("c05_l2", "--", PALETA[0]), ("c05_l8", "-", TINTA)]
    nome = {"c05_l0.5": "$l_v$ 0,5; CoV 0,5", "c10_l0.5": "$l_v$ 0,5; CoV 1,0", "c15_l0.5": "$l_v$ 0,5; CoV 1,5",
            "c05_l2": "$l_v$ 2; CoV 0,5", "c05_l8": "$l_v$ 8; CoV 0,5"}
    d = np.linspace(0, 2.0, 801)
    fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.9), sharey=True)
    for ax, t, L in zip(axs, (12., 18., 24.), "abc"):
        for c, ls, cor in casos:
            z = a["alvo"][c][t]["z"]
            ax.plot(d, cdf_prof(z, d), ls=ls, color=cor, lw=1.0, drawstyle="steps-post", label=nome[c])
            ax.plot(2.0, FIG9[c][(12., 18., 24.).index(t)], "s", ms=5, mfc="none", mec=cor, mew=1.1, clip_on=False)
        ax.set(xlim=(2.0, 0.0), ylim=(0, 1), xlabel="profundidade crítica (m)", title=f"{t:g} h")
        painel(ax, L)
    axs[0].set_ylabel("probabilidade acumulada")
    axs[0].plot([], [], "s", mfc="none", mec=TINTA, label="Cho: na base")
    axs[0].legend(loc="upper left", fontsize=5.8)
    fig.tight_layout()
    salva(fig, "fig9_cdf.pdf")
    # a fracao na base nos 15 casos
    fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.6), sharey=True)
    for ax, t, L in zip(axs, (12., 18., 24.), "abc"):
        for cov in COVS:
            v = [a["alvo"][cid(cov, lv)][t]["base"] for lv in LVS]
            ax.plot(LVS, v, marker=MARC_COV[cov], ms=4, color=COR_COV[cov], label=f"NeoPZ CoV {cov:g}".replace(".", ","))
            for lv in LVS:
                c = cid(cov, lv)
                if c in FIG9:
                    ax.plot(lv, FIG9[c][(12., 18., 24.).index(t)], marker=MARC_COV[cov], ms=7, mfc="none", ls="none",
                            mec=TINTA, mew=1.1)
        ax.set(xlim=(0, 17), ylim=(0, 0.8), xlabel="$l_v$ (m)", title=f"{t:g} h")
        painel(ax, L)
    axs[0].plot([], [], "o", mfc="none", mec=TINTA, ls="none", label="Cho (5 casos)")
    axs[0].set_ylabel("fração na base")
    axs[0].legend(loc="lower right", fontsize=6)
    fig.tight_layout()
    salva(fig, "fig9_base.pdf")


def fig10(r, a):
    fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.1))
    for ax, t, L in zip(axs, (18., 24.), "ab"):
        k = (18., 24.).index(t)
        for cov in COVS:
            cor = COR_COV[cov]
            p = [a["pf"][cid(cov, lv)][t] for lv in LVS]
            ax.errorbar(LVS, [x["P"] for x in p], yerr=[[x["P"] - x["IC"][0] for x in p], [x["IC"][1] - x["P"] for x in p]],
                        color=cor, marker=MARC_COV[cov], ms=4, capsize=2, lw=1.3,
                        label=f"NeoPZ CoV {cov:g}".replace(".", ","))
            ax.plot(LVS, [a["pfpy"][cid(cov, lv)][t]["P"] for lv in LVS], ls="none", marker="+", ms=7, color=cor)
            ax.plot(LVS, [FIG10[(cov, lv)][k] for lv in LVS], ls="--", lw=0.9, marker=MARC_COV[cov], ms=5, mfc="none",
                    color=cor)
            if int(t) in a["lim"]:
                ax.plot(18.5, a["lim"][int(t)][cov], marker="*", ms=8, color=cor, clip_on=False)
        ax.plot([], [], "+", color=TINTA2, ms=7, label="Python (+)")
        ax.plot([], [], "--o", mfc="none", color=TINTA2, lw=0.9, label="Cho (tracejado)")
        ax.plot([], [], "*", color=TINTA2, ms=8, label="limite v.a. ($l_v\\to\\infty$)")
        ax.set(xlim=(0, 19.5), xticks=[0, 4, 8, 12, 16], xlabel="$l_v$ (m)", title=f"{t:g} h")
        painel(ax, L)
    axs[0].set(ylim=(0, 0.16), ylabel="probabilidade de falha")
    axs[1].set(ylim=(0, 0.6))
    axs[1].legend(loc="lower right", fontsize=6.2, ncol=2)
    fig.tight_layout()
    salva(fig, "fig10_pf.pdf")


def fig_pf_t(r, a):
    """P_f(t) nos 15 casos: os 18 h de Cho caem entre os nossos 19 e 21 h"""
    fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.7), sharey=True)
    for ax, cov, L in zip(axs, COVS, "abc"):
        for i, lv in enumerate(LVS):
            c = cid(cov, lv)
            ax.plot(TEMPOS, [a["pf"][c][t]["P"] for t in TEMPOS], marker=".", color=PALETA[[0, 1, 2, 3, 6][i]], lw=1.0,
                    label=f"$l_v$ {lv:g}".replace(".", ","))
            for k, t in enumerate((18., 24.)):
                ax.plot(t, FIG10[(cov, lv)][k], marker="o", mfc="none", ms=5, color=PALETA[[0, 1, 2, 3, 6][i]], ls="none")
        ax.set(xlim=(11, 25), xticks=[12, 16, 18, 20, 22, 24], ylim=(0, 0.6), xlabel="tempo (h)",
               title=f"CoV {cov:g}".replace(".", ","))
        painel(ax, L)
    axs[0].set_ylabel("$P_f$ no instante")
    axs[0].plot([], [], "o", mfc="none", color=TINTA2, label="Cho (18 e 24 h)")
    axs[0].legend(loc="upper left", fontsize=6)
    fig.tight_layout()
    salva(fig, "fig_pf_t.pdf")


def fig_limite(r, a):
    fig, axs = plt.subplots(1, 2, figsize=(7.2, 2.9), gridspec_kw={"width_ratios": [1, 0.8]})
    for T in (18., 20., 24.):
        v = r["kgrade"][int(T)]
        axs[0].plot([x[0] for x in v], [x[1] for x in v], marker="o", ms=2.5, lw=0.8, color=COR_T[int(T)],
                    label=f"{T:g} h")
        axs[0].axvline(r["kstar"][int(T)], color=COR_T[int(T)], lw=0.6, ls=":")
    axs[0].axhline(1, color=TINTA2, lw=0.6)
    axs[0].set(xscale="log", xlim=(0.3, 12), ylim=(0.3, 1.6), xlabel="$k_s/q$ (coluna homogênea)",
               ylabel="FS mínimo na coluna")
    axs[0].legend(loc="lower left")
    painel(axs[0], "a")
    cv = np.linspace(0.2, 1.7, 61)
    for T in (18, 20, 24):
        axs[1].plot(cv, [pf_ln(r["kstar"][T], x) for x in cv], color=COR_T[T], lw=1.2,
                    label=f"$P(k_s>k^*)$, {T} h, $k^*$ = {r['kstar'][T]:.2f}$q$".replace(".", ","))
    for k, T in enumerate((18., 24.)):
        axs[1].plot(COVS, [a["pf"][cid(cov, 16.)][T]["P"] for cov in COVS], "o", ms=4, color=COR_T[int(T)],
                    label=f"NeoPZ $l_v$ 16, {T:g} h")
        axs[1].plot(COVS, [FIG10[(cov, 16.)][k] for cov in COVS], "s", ms=6, mfc="none", mec=COR_T[int(T)], mew=1.2,
                    label=f"Cho $l_v$ 16, {T:g} h")
    axs[1].set(xlim=(0.2, 1.7), ylim=(-0.01, 0.65), xlabel="CoV de $k_s$", ylabel="probabilidade de falha")
    axs[1].legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), fontsize=6)
    painel(axs[1], "b")
    fig.tight_layout()
    salva(fig, "fig_limite.pdf")


def fig_ls(r, a):
    fig, axs = plt.subplots(2, 2, figsize=(7.2, 5.6))
    L = a["ls"]
    ax = axs[0, 0]
    m18 = a["mc18"]
    itens = [("MC Python 5000", m18["py5000"]["P"], m18["py5000"]["SE"], PALETA[1]),
             ("MC Python 55 000", m18["py"]["P"], m18["py"]["SE"], PALETA[1]),
             ("MC NeoPZ 5000", m18["cpp"]["P"], m18["cpp"]["SE"], PALETA[0]),
             ("LS piloto p5000 (200)", L["rec"]["P"], L["rec"]["SE"], TINTA2),
             ("LS p5000 truncado 0-9 (T1)", L["tr10"]["P"], L["tr10"]["SE"], TINTA),
             ("LS p5000 truncado 0-4", L["tr5"]["P"], L["tr5"]["SE"], TINTA2),
             ("LS em estágios (T2)", L["est2"]["P"], L["est2"]["SE"], TINTA),
             ("T1 + T2 (inverso da variância)", L["pool"]["P"], L["pool"]["SE"], PALETA[7])]
    y = np.arange(len(itens))[::-1]
    for yy, (n, p, se, c) in zip(y, itens):
        ax.errorbar(100 * p, yy, xerr=196 * se, fmt="o", ms=4, color=c, capsize=2)
    ax.axvspan(100 * m18["py"]["P"] + 100 * a["desloc"](0.002), 100 * m18["py"]["P"] + 100 * a["desloc"](0.004),
               color=PALETA[3], alpha=0.25, lw=0, label="Python + desvio do FS C++")
    ax.set(yticks=y, yticklabels=[x[0] for x in itens], xlabel="$P_f$(18 h) (%)", xlim=(0.8, 3.0))
    ax.tick_params(axis="y", labelsize=6.5)
    ax.legend(loc="upper right", fontsize=6)
    painel(ax, "a")
    ax = axs[0, 1]
    m24 = a["mc24"]
    itens = [("MC Python 5000", m24["py5000"]["P"], m24["py5000"]["SE"], PALETA[1]),
             ("MC NeoPZ 5000", m24["cpp"]["P"], m24["cpp"]["SE"], PALETA[0]),
             ("LS secante (48)", L["T24"]["P"], L["T24"]["SE"], TINTA2),
             ("LS busca=cerca (48)", L["T24c"]["P"], L["T24c"]["SE"], TINTA)]
    y = np.arange(len(itens))[::-1]
    for yy, (n, p, se, c) in zip(y, itens):
        ax.errorbar(100 * p, yy, xerr=196 * se, fmt="o", ms=4, color=c, capsize=2)
    ax.set(yticks=y, yticklabels=[x[0] for x in itens], xlabel="$P_f$(24 h) (%)", xlim=(30, 41))
    ax.tick_params(axis="y", labelsize=6.5)
    painel(ax, "b")
    ax = axs[1, 0]
    rz = a["raizes"]
    for k, mk, nome in ((-4, "v", "$c_i - 1$ (seguro)"), (0, "o", "$c_i$ (raiz)"), (4, "^", "$c_i + 1$ (falha)")):
        v = [(x[5], x[4] - x[5]) for x in rz if x[2] == k]
        ax.plot([x[0] for x in v], [1000 * x[1] for x in v], mk, ms=4, mfc="none", label=nome)
    ax.axhline(1000 * a["dfs0"], color=PALETA[7], lw=0.8, ls="--", label=f"média na raiz: {1000 * a['dfs0']:.1f}".replace(".", ","))
    ax.axhline(0, color=TINTA2, lw=0.5)
    ax.set(xlabel="FS(18 h), Python", ylabel="FS C++ − FS Python (×10$^{-3}$)")
    ax.legend(loc="lower left", fontsize=6)
    painel(ax, "c")
    ax = axs[1, 1]
    dd = np.linspace(0, 0.006, 61)
    ax.plot(1000 * dd, [100 * a["desloc"](x) for x in dd], color=PALETA[0])
    ax.axvspan(1000 * (-a["dfs0"] - 2 * a["dfs0_se"]), 1000 * (-a["dfs0"] + 2 * a["dfs0_se"]), color=PALETA[3], alpha=0.25,
               lw=0, label="desvio medido nas raízes (±2 EP)")
    ax.set(xlabel="FS deslocado para baixo (×10$^{-3}$)", ylabel="aumento de $P_f$(18 h) (pontos %)")
    ax.legend(loc="upper left", fontsize=6)
    painel(ax, "d")
    fig.tight_layout()
    salva(fig, "fig_ls.pdf")


# ----------------------------------------------------------------------------------------------------------- tabelas
def escreve(nome, txt):
    open(os.path.join(TAB, nome), "w").write("% gerado por gera_dados.py; nao edite\n" + txt)


def tab_fig10(a):
    L = [r"\begin{table}[htbp]\centering\footnotesize\setlength{\tabcolsep}{3.5pt}",
         r"\caption{A Fig.~10: $P_f$ no instante (FS mínimo $< 1$) em 18 e 24~h, em \%. NeoPZ: 5000 realizações, IC~95\,\% "
         r"de Wilson; Python: as mesmas realizações (55\,000 no caso base); Cho: a figura digitalizada ($\pm 0{,}2$ ponto). "
         r"$\Delta$: NeoPZ $-$ Cho, em pontos; $t_{eq}$: o instante em que a nossa $P_f$ iguala a de Cho em 18~h.}"
         r"\label{tab:fig10}",
         r"\begin{tabular}{llrrrrrrrrr}\toprule",
         r"& & \multicolumn{4}{c}{18 h} & & \multicolumn{4}{c}{24 h} \\\cmidrule{3-6}\cmidrule{8-11}",
         r"CoV & $l_v$ (m) & NeoPZ & Python & Cho & $\Delta$ & $t_{eq}$ (h) & NeoPZ [IC 95\,\%] & Python & Cho & $\Delta$ \\\midrule"]
    for cov in COVS:
        for lv in LVS:
            c = cid(cov, lv)
            p18, p24 = a["pf"][c][18.], a["pf"][c][24.]
            y18, y24 = a["pfpy"][c][18.], a["pfpy"][c][24.]
            c18, c24 = FIG10[(cov, lv)]
            te = a["t_equiv"][c]
            L.append(f"{br(cov, 1)} & {lv:g} & {pct(p18['P'], 2)} & {pct(y18['P'], 2)} & {pct(c18)} & "
                     f"{sinal(100 * (p18['P'] - c18))} & {br(te, 1) if te else '--'} & "
                     f"{pct(p24['P'])} [{pct(p24['IC'][0])}; {pct(p24['IC'][1])}] & {pct(y24['P'])} & {pct(c24)} & "
                     f"{sinal(100 * (p24['P'] - c24))} \\\\".replace("lv", ""))
        if cov != COVS[-1]:
            L.append(r"\addlinespace")
    L.append(r"\bottomrule\end{tabular}\end{table}")
    escreve("fig10.tex", "\n".join(L) + "\n")


def tab_alvos(a):
    L = [r"\begin{table}[htbp]\centering\footnotesize\setlength{\tabcolsep}{4pt}",
         r"\caption{Os alvos das Figs.~5 a 9: a fração de superfícies críticas na base e a mediana de $z_{crit}$ em "
         r"12/18/24~h (NeoPZ, 5000 realizações; Cho entre parênteses, lido).}\label{tab:alvos}",
         r"\begin{tabular}{lccc}\toprule",
         r"caso & base, 12 h & base, 18 h & base, 24 h \\\midrule"]
    for c, nome in (("c05_l0.5", r"$l_v$ 0,5; CoV 0,5"), ("c10_l0.5", r"$l_v$ 0,5; CoV 1,0"),
                    ("c15_l0.5", r"$l_v$ 0,5; CoV 1,5"), ("c05_l2", r"$l_v$ 2; CoV 0,5"), ("c05_l8", r"$l_v$ 8; CoV 0,5")):
        L.append(nome + " & " + " & ".join(f"{br(a['alvo'][c][t]['base'], 3)} ({br(FIG9[c][i], 3)})"
                                          for i, t in enumerate((12., 18., 24.))) + r" \\")
    L.append(r"\midrule mediana de $z_{crit}$, caso base (m) & " + " & ".join(
        f"{br(a['alvo']['c05_l0.5'][t]['med'], 3)} ({br(FIG7_MED[int(t)])})" for t in (12., 18., 24.)) + r" \\")
    L.append(r"$z_{crit} < 0{,}75$~m, caso base & " + " & ".join(
        f"{pct(a['alvo']['c05_l0.5'][t]['acima075'])}\\,\\% (0)" for t in (12., 18., 24.)) + r" \\")
    L.append(r"\bottomrule\end{tabular}\end{table}")
    escreve("alvos.tex", "\n".join(L) + "\n")


def tab_estudos(a):
    E = a["est"]
    nome = {"c05_l0.5_ny200": r"ny 200 ($l_v$ 0,5)", "c05_l16_ny200": r"ny 200 ($l_v$ 16)",
            "c05_l0.5_tri": r"tri, massa concentrada", "c05_l0.5_tricons": r"tri, massa consistente",
            "c05_l0.5_eps1e-4": r"$\varepsilon = 10^{-4}$ (descartado)"}
    L = [r"\begin{table}[htbp]\centering\footnotesize\setlength{\tabcolsep}{3.5pt}",
         r"\caption{Os estudos do método (CoV 0,5), pareados com as realizações de mesmo número da produção (quad, ny 80, "
         r"$\varepsilon = 10^{-3}$). $b/c$: pares discordantes (só o estudo falha / só a produção falha). $P_f$ em \%.}"
         r"\label{tab:estudos}",
         r"\begin{tabular}{lrrrrrrrr}\toprule",
         r"estudo & $N$ & paradas & $P_f$ 18 h & produção & $b/c$ & $P_f$ 24 h & produção & $b/c$ \\\midrule"]
    for c in ("c05_l0.5_ny200", "c05_l16_ny200", "c05_l0.5_tri", "c05_l0.5_tricons", "c05_l0.5_eps1e-4"):
        e = E[c]
        L.append(f"{nome[c]} & {inteiro(e['n'])} & {inteiro(e['pararam'])} & {pct(e[18.]['Pa_ok'], 2)} & "
                 f"{pct(e[18.]['Pprod'], 2)} & {e[18.]['b']}/{e[18.]['c']} & {pct(e[24.]['Pa_ok'])} & {pct(e[24.]['Pprod'])} "
                 f"& {e[24.]['b']}/{e[24.]['c']} \\\\")
    L.append(r"\bottomrule\end{tabular}\end{table}")
    escreve("estudos.tex", "\n".join(L) + "\n")


def tab_ls(a):
    L = a["ls"]
    m = a["mc18"]["py"]
    linhas = [("piloto p5000, sem truncar", "rec", "7"), ("o mesmo, linhas 0--99", "rec100", "7"),
              ("piloto p500 (6 falhas)", "p500", "7"), ("T1: p5000 truncado 0--9", "tr10", "7"),
              ("T1: p5000 truncado 0--4", "tr5", "7"), ("T2: estágio 2 (direção das linhas)", "est2", "9")]
    T = [r"\begin{table}[htbp]\centering\footnotesize\setlength{\tabcolsep}{4pt}",
         r"\caption{O line sampling no caso base (CoV 0,5; $l_v$ 0,5~m), $M = 100$, NeoPZ quad ny 80. $z$: contra o MC "
         r"em Python de 55\,000 ($P_f(18\,\mathrm{h}) = " + pct(m["P"], 2) + r"\,\%$). CoV por linha: o desvio de $p$ "
         r"entre as linhas sobre a média.}\label{tab:ls}",
         r"\begin{tabular}{llrrrrrr}\toprule",
         r"direção & semente & linhas & $P_f$ (\%) & CoV (\%) & CoV por linha & análises/linha & $z$ \\\midrule",
         r"\multicolumn{8}{l}{\emph{18 h}} \\"]
    for nome, k, s in linhas:
        e = L[k]
        z = (e["P"] - m["P"]) / math.hypot(e["SE"], m["SE"])
        T.append(f"{nome} & {s} & {e['n']} & {pct(e['P'], 2)} & {pct(e['CoV'])} & {br(e['CoVl'])} & {br(e['an'], 2)} & "
                 f"{sinal(z, 2)} \\\\")
    p = L["pool"]
    z = (p["P"] - m["P"]) / math.hypot(p["SE"], m["SE"])
    T.append(f"T1 0--9 + T2 (inverso da variância) & & 200 & {pct(p['P'], 2)} & {pct(p['SE'] / p['P'])} & & & {sinal(z, 2)} \\\\")
    m24 = a["mc24"]["py5000"]
    T.append(r"\multicolumn{8}{l}{\emph{24 h} ($z$ contra o MC em Python de 5000, $" + pct(m24["P"], 2) + r"\,\%$)} \\")
    for nome, k in (("piloto p5000, secante", "T24"), ("piloto p5000, \\texttt{busca=cerca}", "T24c")):
        e = L[k]
        z = (e["P"] - m24["P"]) / math.hypot(e["SE"], m24["SE"])
        T.append(f"{nome} & 7 & {e['n']} & {pct(e['P'], 2)} & {pct(e['CoV'])} & {br(e['CoVl'])} & {br(e['an'], 2)} & "
                 f"{sinal(z, 2)} \\\\")
    T.append(r"\bottomrule\end{tabular}\end{table}")
    escreve("ls.tex", "\n".join(T) + "\n")


# ----------------------------------------------------------------------------------------------------------- numeros
def numeros(r, a):
    N = {}
    D = r["det_csv"]
    for t, nm in ((0., "Zero"), (6., "Seis"), (12., "Doze"), (18., "Dezoito"), (24., "Vquatro")):
        N["fs" + nm] = br(float(D[t]["FS_min"]), 3)
    for t, nm in ((6., "Seis"), (12., "Doze"), (18., "Dezoito")):
        N["frente" + nm] = br(float(D[t]["frente_m"]), 2)
    N["tPassagem"] = br(float(D[24.]["t_passagem"]), 2)
    N["zcDoze"], N["zcDezoito"] = br(float(D[12.]["z_crit"]), 3), br(float(D[18.]["z_crit"]), 3)
    # o limite de variavel aleatoria
    for T, nm in ((18, "Dezoito"), (20, "Vinte"), (24, "Vquatro")):
        N["kstar" + nm] = br(r["kstar"][T], 2 if T < 24 else 3)
        for cov, cn in zip(COVS, ("A", "B", "C")):
            N["lim" + nm + cn] = br(a["lim"][T][cov], 3)
    for cov, cn in zip(COVS, ("A", "B", "C")):
        N["mcDezesseis" + cn] = br(a["pf"][cid(cov, 16.)][24.]["P"], 3)
        N["mcDezesseisDezoito" + cn] = br(a["pf"][cid(cov, 16.)][18.]["P"], 4)
    # 24 h contra Cho
    d24 = [100 * (a["pf"][cid(cov, lv)][24.]["P"] - FIG10[(cov, lv)][1]) for cov in COVS for lv in LVS]
    d18 = [100 * (a["pf"][cid(cov, lv)][18.]["P"] - FIG10[(cov, lv)][0]) for cov in COVS for lv in LVS]
    N["difVqMin"], N["difVqMax"] = br(-max(d24), 1), br(-min(d24), 1)
    N["difDzMin"], N["difDzMax"] = br(-max(d18), 1), br(-min(d18), 1)
    zpy = [abs(a["pf"][c][t]["P"] - a["pfpy"][c][t]["P"]) / math.hypot(a["pf"][c][t]["SE"], a["pfpy"][c][t]["SE"])
           for c in a["pf"] for t in (18., 24.) if a["pf"][c][t]["SE"] > 0 or a["pfpy"][c][t]["SE"] > 0]
    N["zPyMax"] = br(max(zpy), 2)
    dpy = [100 * abs(a["pf"][c][t]["P"] - a["pfpy"][c][t]["P"]) for c in a["pf"] for t in (18., 24.)]
    N["difPyMax"] = br(max(dpy), 2)
    te = [v for v in a["t_equiv"].values() if v]
    N["tEqMin"], N["tEqMax"] = br(min(te), 1), br(max(te), 1)
    N["nTeqSem"] = str(sum(1 for v in a["t_equiv"].values() if not v))
    b = a["pf"]["c05_l0.5"]
    N["pfBaseDz"], N["pfBaseVq"] = pct(b[18.]["P"], 2), pct(b[24.]["P"], 1)
    N["pfBaseVqLo"], N["pfBaseVqHi"] = pct(b[24.]["IC"][0], 1), pct(b[24.]["IC"][1], 1)
    N["pfBaseDzLo"], N["pfBaseDzHi"] = pct(b[18.]["IC"][0], 2), pct(b[18.]["IC"][1], 2)
    N["pfPyBaseDz"] = pct(a["mc18"]["py"]["P"], 2)
    N["pfPyBaseDzSE"] = pct(a["mc18"]["py"]["SE"], 2)
    N["pfPyBaseVq"] = pct(a["mc24"]["py"]["P"], 1)
    N["pfPyCincoDz"] = pct(a["mc18"]["py5000"]["P"], 2)
    N["pfPyCincoVq"] = pct(a["mc24"]["py5000"]["P"], 2)
    N["pfBaseDoze"] = pct(b[12.]["P"], 2)
    N["pfDozeMax"] = pct(max(a["pf"][c][12.]["P"] for c in a["pf"]), 2)
    x = a["alvo"]["c05_l0.5"]
    N["baseDoze"], N["baseDezoito"], N["baseVq"] = br(x[12.]["base"], 3), br(x[18.]["base"], 3), br(x[24.]["base"], 3)
    N["medDoze"], N["medDezoito"], N["medVq"] = br(x[12.]["med"], 3), br(x[18.]["med"], 3), br(x[24.]["med"], 3)
    f = x["fs24"]
    N["fsMedia"], N["fsDp"], N["fsMin"], N["fsMax"], N["fsPico"] = (br(f["media"], 3), br(f["dp"], 3), br(f["min"], 4),
                                                                   br(f["max"], 3), br(f["pico"], 3))
    # estudos
    E = a["est"]
    e = E["c05_l0.5_eps1e-4"]
    N["epsN"], N["epsParou"] = inteiro(e["n"]), inteiro(e["pararam"])
    N["epsParouPct"] = br(100 * e["pararam"] / e["n"], 0)
    N["epsCens"] = inteiro(e[24.]["cens"])
    N["epsLimLo"], N["epsLimHi"] = pct(e[24.]["lim"][0], 1), pct(e[24.]["lim"][1], 1)
    # LS
    L = a["ls"]
    for k, nm in (("rec", "Rec"), ("rec100", "RecCem"), ("tr10", "Tdez"), ("tr5", "Tcinco"), ("est2", "Est"),
                  ("p500", "Pqui"), ("T24", "Vq"), ("T24c", "VqC")):
        N["ls" + nm] = pct(L[k]["P"], 2)
        N["lsCv" + nm] = br(L[k]["CoVl"], 2)
        N["lsAn" + nm] = br(L[k]["an"], 1)
    N["lsPool"], N["lsPoolSE"] = pct(L["pool"]["P"], 2), pct(L["pool"]["SE"], 2)
    m = a["mc18"]["py"]
    N["zPool"] = br((L["pool"]["P"] - m["P"]) / math.hypot(L["pool"]["SE"], m["SE"]), 1)
    N["dfsRaiz"] = br(-1000 * a["dfs0"], 1)
    N["dfsRaizSE"] = br(1000 * a["dfs0_se"], 1)
    N["deslocRaiz"] = br(100 * a["desloc"](-a["dfs0"]), 2)
    N["deslocTres"] = br(100 * a["desloc"](0.003), 2)
    N["deslocQuatro"] = br(100 * a["desloc"](0.004), 2)
    corr = m["P"] + a["desloc"](-a["dfs0"])
    N["pyCorr"] = pct(corr, 2)
    N["zPoolCorr"] = br((L["pool"]["P"] - corr) / math.hypot(L["pool"]["SE"], m["SE"]), 1)
    m24 = a["mc24"]["py5000"]
    N["zVq"] = br((L["T24"]["P"] - m24["P"]) / math.hypot(L["T24"]["SE"], m24["SE"]), 2)
    N["zVqC"] = br((L["T24c"]["P"] - m24["P"]) / math.hypot(L["T24c"]["SE"], m24["SE"]), 2)
    N["nRaizes"] = str(len([x for x in a["raizes"] if x[2] == 0]))
    txt = "% gerado por gera_dados.py; nao edite\n" + "".join(f"\\newcommand{{\\{k}}}{{{v}}}\n" for k, v in N.items())
    open(os.path.join(AQUI, "numeros.tex"), "w").write(txt)
    return N


def main():
    for p in (DADOS, FIGS, TAB):
        os.makedirs(p, exist_ok=True)
    r = carrega()
    a = contas(r)
    fig3(r)
    fig5(r, a)
    fig6(r, a)
    fig7(r, a)
    fig8(r, a)
    fig9(r, a)
    fig10(r, a)
    fig_pf_t(r, a)
    fig_limite(r, a)
    fig_ls(r, a)
    tab_fig10(a)
    tab_alvos(a)
    tab_estudos(a)
    tab_ls(a)
    N = numeros(r, a)
    for k in sorted(N):
        print(f"{k:22s} {N[k]}")


if __name__ == "__main__":
    main()
