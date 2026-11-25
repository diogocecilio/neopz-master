#!/usr/bin/env python3
"""resumo_frag.py: a verificacao da curva de fragilidade P_f(I) = Phi(-ln F(I) / beta) (CAMPANHA=frag, casos_frag.txt).

Para cada caso crpc= (primeira passagem): a P_f de cada intensidade e de cada duracao D <= crpc_dmax, P(D_C(I) <= D),
com IC de Wilson de 95 %, contra a prevista Phi(-ln F(I) / beta), com F(I) do deterministico (os pontos e a Eq. 6 do
Rahardjo ajustada a eles) e beta do MC de calibracao da variante (o desvio de ln FS no fim da chuva). Quem nao rompe numa
intensidade nao roda as menores: conta como sobrevivente nelas (monotonia em I). Imprime a tabela, grava o JSON e as
figuras (P_f x I e a distribuicao de t_1).

Uso: python3 resumo_frag.py <pasta da campanha (com f50_iii/, f50_i/, ...)> [--ref <pasta da referencia chuva2/ref>]
     [--det50 <runs/sigmoide_f50>] [--det100 <runs/sigmoide_f100_72h>] [--json <arquivo>] [--figs <pasta>]
"""
import argparse
import csv
import glob
import json
import math
import os
import re

from scipy.optimize import brentq, curve_fit
from scipy.stats import norm

CASOS = {  # caso crpc -> (solo, variante, D da P_f prevista, caso de calibracao de beta (MC), instante do FS no MC,
    #           a intensidade da conferencia contra a calibracao)
    "f50_iii": ("f50", "iii", 24, "REF", 24.0, 31.05),
    "f50_i": ("f50", "i", 24, "f50_i_mc20", 24.0, 20.0),
    "f100_iii": ("f100", "iii", 72, "f100_iii_mc", 72.0, 9.0),
    "f100_i": ("f100", "i", 72, "f100_i_mc", 72.0, 3.6),
}


def wilson(k, n, z=1.96):
    if n == 0:
        return float("nan"), float("nan")
    p = k / n
    c = (p + z * z / (2 * n)) / (1 + z * z / n)
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return c - h, c + h


def le_crpc(pasta):
    """{real: {I: (falha, DC)}}: falha 'sim' (DC = t_1), 'nao', ou erro; as linhas de todos os blocos"""
    r = {}
    for arq in glob.glob(os.path.join(pasta, "*_crpc_*.csv")):
        with open(arq) as f:
            for row in csv.DictReader(f):
                n, I = int(row["real"]), float(row["I"])
                if n < 0:  # a deterministica (det=1 no primeiro bloco): fica de fora da P_f
                    continue
                r.setdefault(n, {})[I] = (row["falha"], float(row["DC_h"]))
    return r


def le_mc(pasta, t):
    """{real: FS no instante t} de um MC (mc_*.csv, sem os mapas)"""
    fs = {}
    for arq in glob.glob(os.path.join(pasta, "mc_*.csv")):
        if "_mapa_" in arq or "_crpc_" in arq or "_chuva_" in arq:
            continue
        with open(arq) as f:
            for row in csv.DictReader(f):
                n = int(row["real"])
                if n >= 0 and abs(float(row["t_h"]) - t) < 1e-9 and row["situacao"] == "ok":
                    fs[n] = float(row["FS_min"])
    return fs


def le_det(pasta, prefixo, t):
    """{I: FS das medias no instante t} das pastas <prefixo>I<I>/det_0_-1.csv"""
    d = {}
    for p in glob.glob(os.path.join(pasta, prefixo + "I*")):
        m = re.match(re.escape(prefixo) + r"I([\d.]+)$", os.path.basename(p))
        arq = os.path.join(p, "det_0_-1.csv")
        if not m or not os.path.exists(arq):
            continue
        with open(arq) as f:
            for row in csv.DictReader(f):
                if abs(float(row["t_h"]) - t) < 1e-9 and row["situacao"] == "ok":
                    d[float(m.group(1))] = float(row["FS_min"])
    return dict(sorted(d.items()))


def sig(x, A, B, C, D):
    return A + B / (1 + math.exp(-(x - C) / D)) if abs((x - C) / D) < 700 else (A if x < C else A + B)


def ajusta_eq6(det):
    import numpy as np
    I, F = np.array(list(det)), np.array(list(det.values()))
    f = np.vectorize(sig)
    p, _ = curve_fit(lambda x, A, B, C, D: f(x, A, B, C, D), I, F, p0=[F.max(), F.min() - F.max(), np.median(I), 10],
                     maxfev=50000)
    return tuple(p)


def interp_log(det, I):
    """ln F interpolado linearmente em I entre os pontos deterministicos"""
    xs = list(det)
    if I <= xs[0]:
        return math.log(det[xs[0]])
    for a, b in zip(xs, xs[1:]):
        if a <= I <= b:
            la, lb = math.log(det[a]), math.log(det[b])
            return la + (lb - la) * (I - a) / (b - a)
    return math.log(det[xs[-1]])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("pasta")
    ap.add_argument("--ref", default=os.path.expanduser("~/projects/runs/fermi_chuva2/chuva2/ref"))
    ap.add_argument("--det50", default=os.path.expanduser("~/projects/runs/sigmoide_f50"))
    ap.add_argument("--det100", default=os.path.expanduser("~/projects/runs/sigmoide_f100_72h"))
    ap.add_argument("--json")
    ap.add_argument("--figs")
    a = ap.parse_args()
    saida = {}
    for caso, (solo, var, D, cal, tcal, Iconf) in CASOS.items():
        pasta = os.path.join(a.pasta, caso)
        if not os.path.isdir(pasta):
            continue
        dados = le_crpc(pasta)
        if not dados:
            continue
        # o deterministico (F no fim da chuva) e a Eq. 6
        if solo == "f50":
            det = le_det(a.det50, "" if var == "iii" else "i_", float(D))
        else:
            det = le_det(a.det100, f"{var}_", float(D))
        eq6 = ajusta_eq6(det) if len(det) >= 5 else None
        # beta e a mediana da calibracao
        fcal = le_mc(a.ref if cal == "REF" else os.path.join(a.pasta, cal), tcal)
        lns = [math.log(v) for v in fcal.values() if v > 0]
        beta = (sum((x - sum(lns) / len(lns)) ** 2 for x in lns) / (len(lns) - 1)) ** 0.5 if len(lns) > 10 else None
        Is = sorted({I for r in dados.values() for I in r}, reverse=True)
        # realizacoes validas: sem erro em nenhuma intensidade rodada
        validas = [n for n, r in dados.items() if all(f in ("sim", "nao") for f, _ in r.values())]
        erros = len(dados) - len(validas)
        res = {"solo": solo, "variante": var, "D": D, "N": len(validas), "erros": erros, "beta": beta,
               "n_calibracao": len(lns), "eq6": eq6, "intensidades": {}}
        print(f"\n== {caso}: {len(validas)} realizacoes validas ({erros} com erro), chuva de {D} h; beta = "
              + (f"{beta:.4f} ({len(lns)} do MC de calibracao)" if beta else "sem calibracao")
              + ("; Eq. 6: A={:.4f} B={:.4f} C={:.3f} D={:.3f}".format(*eq6) if eq6 else ""))
        print(f"  {'I':>6} {'F det':>7} {'P_f MC':>8} {'IC 95%':>17} {'prev. pontos':>13} {'z':>6} {'prev. Eq.6':>11}"
              f" {'z':>6}   t_1 (h): mediana [10%, 90%] das falhas")
        for I in Is:
            dc = []
            for n in validas:
                r = dados[n]
                if I in r:
                    dc.append(r[I][1] if r[I][0] == "sim" else math.inf)
                else:  # nao rodada: sobreviveu numa intensidade maior (monotonia em I)
                    dc.append(math.inf)
            N = len(dc)
            k = sum(1 for x in dc if x <= D + 1e-9)
            p = k / N
            lo, hi = wilson(k, N)
            se = math.sqrt(max(p * (1 - p), 1e-12) / N)
            Fd = math.exp(interp_log(det, I)) if det else None
            prev = norm.cdf(-math.log(Fd) / beta) if (beta and Fd) else None
            prev6 = norm.cdf(-math.log(sig(I, *eq6)) / beta) if (beta and eq6) else None
            fal = sorted(x for x in dc if x <= D)
            q = (lambda f: fal[min(len(fal) - 1, int(f * len(fal)))]) if fal else None
            zp = (p - prev) / se if prev is not None else None
            z6 = (p - prev6) / se if prev6 is not None else None
            durs = {str(d): sum(1 for x in dc if x <= d) / N for d in (6, 12, 18, 24, 36, 48, 60, 72, 96, 120)
                    if d <= max(D, max((x for x in dc if x < math.inf), default=0))}
            print(f"  {I:6g} {Fd if Fd else float('nan'):7.4f} {100 * p:7.2f}% [{100 * lo:6.2f}, {100 * hi:6.2f}]%"
                  f" {100 * prev if prev is not None else float('nan'):12.2f}% {zp if zp is not None else float('nan'):6.2f}"
                  f" {100 * prev6 if prev6 is not None else float('nan'):10.2f}% {z6 if z6 is not None else float('nan'):6.2f}"
                  + (f"   {q(0.5):.1f} [{q(0.1):.1f}, {q(0.9):.1f}] ({k} falhas)" if fal else ""))
            res["intensidades"][str(I)] = {"N": N, "falhas": k, "Pf": p, "IC": [lo, hi], "F_det": Fd, "Pf_prev": prev,
                                           "z": zp, "Pf_prev_eq6": prev6, "z_eq6": z6, "Pf_duracao": durs,
                                           "t1_falhas": fal}
        # conferencia embutida: a referencia (iii, 31,05) ou o MC de calibracao, nas mesmas realizacoes
        if Iconf in Is and fcal:
            comuns = [n for n in validas if n in fcal]
            difer = 0
            for n in comuns:
                r = dados[n]
                f_crpc = Iconf in r and r[Iconf][0] == "sim" and r[Iconf][1] <= tcal + 1e-9
                difer += f_crpc != (fcal[n] < 1)
            print(f"  conferencia em {Iconf:g} mm/h contra o MC ({len(comuns)} realizacoes em comum): "
                  f"{difer} com falha diferente")
            res["conferencia"] = {"I": Iconf, "comuns": len(comuns), "diferentes": difer}
        saida[caso] = res
    if a.json:
        with open(a.json, "w") as f:
            json.dump(saida, f, indent=1)
    if a.figs and saida:
        figuras(saida, a.figs)


def figuras(saida, pasta):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    os.makedirs(pasta, exist_ok=True)
    cor = {"iii": "#2a78d6", "i": "#1baf7a"}
    for solo in ("f50", "f100"):
        casos = [c for c, r in saida.items() if r["solo"] == solo]
        if not casos:
            continue
        fig, ax = plt.subplots(figsize=(7.5, 4.6))
        for c in casos:
            r = saida[c]
            Is = sorted(float(x) for x in r["intensidades"])
            P = [r["intensidades"][str(I)] for I in Is]
            ax.errorbar(Is, [100 * x["Pf"] for x in P], yerr=[[100 * (x["Pf"] - x["IC"][0]) for x in P],
                        [100 * (x["IC"][1] - x["Pf"]) for x in P]], fmt="o", ms=7, color=cor[r["variante"]], capsize=3,
                        label=f"MC, variante ({r['variante']})")
            if r["beta"] and r["eq6"]:
                import numpy as np
                xs = np.linspace(min(Is) * 0.8, max(Is) * 1.1, 200)
                ax.plot(xs, [100 * norm.cdf(-math.log(sig(x, *r["eq6"])) / r["beta"]) for x in xs], color=cor[r["variante"]],
                        lw=2, label=f"Φ(−ln F/β), ({r['variante']}): Eq. 6 + β = {r['beta']:.3f}")
        ax.set_xlabel("intensidade da chuva (mm/h)")
        ax.set_ylabel("P_f (%)")
        ax.set_title(f"{solo}: curva de fragilidade (chuva de {saida[casos[0]]['D']} h)", loc="left")
        ax.grid(True, color="#e4e3df")
        ax.legend(frameon=False, fontsize=8)
        fig.tight_layout()
        fig.savefig(os.path.join(pasta, f"fragilidade_{solo}.png"), dpi=140)
        fig.savefig(os.path.join(pasta, f"fragilidade_{solo}.pdf"))
        plt.close(fig)


if __name__ == "__main__":
    main()
