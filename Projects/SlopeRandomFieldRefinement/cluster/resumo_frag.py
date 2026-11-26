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


def bloco_pronto(arq, tipo):
    """o bloco do arquivo terminou? (o mapa mc_mapa_<a>_<z>.csv e gravado no fim do bloco; antes disso o csv tem as
    realizacoes e as intensidades ja feitas, e uma realizacao em curso pareceria ter sobrevivido nas intensidades que
    ainda nao rodou)"""
    m = re.search(r"_" + tipo + r"_(\d+)_(\d+)\.csv$", arq) if tipo else re.search(r"mc_(\d+)_(\d+)\.csv$", arq)
    return bool(m) and os.path.exists(os.path.join(os.path.dirname(arq), f"mc_mapa_{m.group(1)}_{m.group(2)}.csv"))


def prefixo(pasta):
    """o fim do prefixo continuo de blocos terminados (realizacoes 0..k-1 sem buraco): os blocos que terminam primeiro
    sao os mais baratos (realizacoes que nao rompem rodam uma intensidade so), e usar todos os prontos enviesaria a P_f
    para baixo; o prefixo e uma amostra sem vies (a fila pega os blocos em ordem)"""
    feitos = {}
    for m_ in (re.search(r"mc_mapa_(\d+)_(\d+)\.csv$", a) for a in glob.glob(os.path.join(pasta, "mc_mapa_*.csv"))):
        if m_:
            feitos[int(m_.group(1))] = int(m_.group(2))
    k = 0
    while k in feitos:
        k = feitos[k] + 1
    return k


def le_crpc(pasta):
    """{real: {I: (falha, DC)}}: falha 'sim' (DC = t_1), 'nao', ou erro; so o prefixo continuo de blocos terminados"""
    r = {}
    kmax = prefixo(pasta)
    for arq in glob.glob(os.path.join(pasta, "*_crpc_*.csv")):
        if not bloco_pronto(arq, "crpc"):
            continue
        with open(arq) as f:
            for row in csv.DictReader(f):
                n, I = int(row["real"]), float(row["I"])
                if n < 0 or n >= kmax:  # a deterministica (det=1 no primeiro bloco) e o que passa do prefixo
                    continue
                r.setdefault(n, {})[I] = (row["falha"], float(row["DC_h"]))
    return r


def le_mc(pasta, t):
    """{real: FS no instante t} de um MC (mc_*.csv, sem os mapas; so o prefixo continuo de blocos terminados)"""
    fs = {}
    kmax = prefixo(pasta)
    for arq in glob.glob(os.path.join(pasta, "mc_*.csv")):
        if "_mapa_" in arq or "_crpc_" in arq or "_chuva_" in arq or not bloco_pronto(arq, ""):
            continue
        with open(arq) as f:
            for row in csv.DictReader(f):
                n = int(row["real"])
                if 0 <= n < kmax and abs(float(row["t_h"]) - t) < 1e-9 and row["situacao"] == "ok":
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
    ap.add_argument("--det50k", default=os.path.expanduser("~/projects/runs/sigmoide_f50_kef"),
                    help="o F das medias com o k_s mediano (geometrico, 8,57e-6 m/s): pastas [i_]ks8.57e-6_I<I>")
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
        # o F com o k_s mediano (a condutividade efetiva do campo lognormal; so o f50 por enquanto)
        detk = le_det(a.det50k, ("" if var == "iii" else "i_") + "ks8.57e-6_", float(D)) if solo == "f50" else {}
        # beta e a mediana da calibracao
        fcal = le_mc(a.ref if cal == "REF" else os.path.join(a.pasta, cal), tcal)
        lns = [math.log(v) for v in fcal.values() if v > 0]
        beta = (sum((x - sum(lns) / len(lns)) ** 2 for x in lns) / (len(lns) - 1)) ** 0.5 if len(lns) > 10 else None
        # m' = mediana da calibracao / F(k_s mediano) na intensidade da calibracao (a da conferencia)
        mk = (math.exp(sum(lns) / len(lns)) / detk[Iconf]) if (beta and Iconf in detk) else None
        Is = sorted({I for r in dados.values() for I in r}, reverse=True)
        # realizacoes validas: sem erro em nenhuma intensidade rodada
        validas = [n for n, r in dados.items() if all(f in ("sim", "nao") for f, _ in r.values())]
        erros = len(dados) - len(validas)
        res = {"solo": solo, "variante": var, "D": D, "N": len(validas), "erros": erros, "beta": beta, "m_kmed": None,
               "n_calibracao": len(lns), "eq6": eq6, "intensidades": {}}
        print(f"\n== {caso}: {len(validas)} realizacoes validas ({erros} com erro), chuva de {D} h; beta = "
              + (f"{beta:.4f} ({len(lns)} do MC de calibracao)" if beta else "sem calibracao")
              + ("; Eq. 6: A={:.4f} B={:.4f} C={:.3f} D={:.3f}".format(*eq6) if eq6 else "")
              + (f"; k_s mediano: m' = {mk:.4f}" if mk else ""))
        print("  (z = (P_f MC - prev.) / dp binomial da prev.; prev. pontos: F das medias com o k_s medio, m = 1; "
              "prev. k_s med.: F com o k_s mediano e o m' da calibracao, Phi(-ln(m' F)/beta))")
        print(f"  {'I':>6} {'F det':>7} {'P_f MC':>8} {'IC 95%':>17} {'prev. pontos':>13} {'z':>6} {'prev. Eq.6':>11}"
              f" {'z':>6} {'prev. k_s med.':>14} {'z':>6}   t_1 (h): mediana [10%, 90%] das falhas")
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
            # z com o desvio padrao binomial da previsao (o teste da hipotese "a previsao esta certa"); com o da P_f
            # observada, P_f = 0 dava dp 0 e z absurdo
            se = (lambda q: math.sqrt(max(q * (1 - q), 1e-12) / N))
            Fd = math.exp(interp_log(det, I)) if det else None
            prev = norm.cdf(-math.log(Fd) / beta) if (beta and Fd) else None
            prev6 = norm.cdf(-math.log(sig(I, *eq6)) / beta) if (beta and eq6) else None
            Fk = detk.get(I) if detk else None
            prevk = norm.cdf(-math.log(mk * Fk) / beta) if (mk and Fk) else None
            fal = sorted(x for x in dc if x <= D)
            q = (lambda f: fal[min(len(fal) - 1, int(f * len(fal)))]) if fal else None
            zp = (p - prev) / se(prev) if prev is not None else None
            z6 = (p - prev6) / se(prev6) if prev6 is not None else None
            zk = (p - prevk) / se(prevk) if prevk is not None else None
            durs = {str(d): sum(1 for x in dc if x <= d) / N for d in (6, 12, 18, 24, 36, 48, 60, 72, 96, 120)
                    if d <= max(D, max((x for x in dc if x < math.inf), default=0))}
            print(f"  {I:6g} {Fd if Fd else float('nan'):7.4f} {100 * p:7.2f}% [{100 * lo:6.2f}, {100 * hi:6.2f}]%"
                  f" {100 * prev if prev is not None else float('nan'):12.2f}% {zp if zp is not None else float('nan'):6.2f}"
                  f" {100 * prev6 if prev6 is not None else float('nan'):10.2f}% {z6 if z6 is not None else float('nan'):6.2f}"
                  f" {100 * prevk if prevk is not None else float('nan'):13.2f}% {zk if zk is not None else float('nan'):6.2f}"
                  + (f"   {q(0.5):.1f} [{q(0.1):.1f}, {q(0.9):.1f}] ({k} falhas)" if fal else ""))
            res["intensidades"][str(I)] = {"N": N, "falhas": k, "Pf": p, "IC": [lo, hi], "F_det": Fd, "Pf_prev": prev,
                                           "z": zp, "Pf_prev_eq6": prev6, "z_eq6": z6, "F_kmed": Fk, "Pf_prev_kmed": prevk,
                                           "z_kmed": zk, "Pf_duracao": durs,
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
        res["m_kmed"] = mk
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
