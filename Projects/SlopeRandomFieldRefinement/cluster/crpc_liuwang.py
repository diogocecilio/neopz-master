#!/usr/bin/env python3
"""crpc_liuwang.py: a P_F1 por evento de chuva do caso C (Liu e Wang 2023, Comput. Geotech. 162:105706) a partir das
curvas criticas (CRPC) do SlopeRainMonteCarlo (opcao crpc=, arquivos <saida>_crpc_<a>_<b>.csv), pelo metodo
semi-analitico do artigo (Secao 3, Eqs. 9 a 14), e o Monte Carlo direto (chuva_aleatoria=1, <saida>_<a>_<b>.csv com
<saida>_chuva_<a>_<b>.csv). Sem numpy (roda no fermi).

A chuva (Tabela 2 do artigo, Jun et al. 2017): I e D Pareto generalizadas, F(x) = 1 - (1 + xi (x - mu)/sigma)^(-1/xi),
ligadas pela copula de Frank, C(u, v) = -ln(1 + (e^-wu - 1)(e^-wv - 1)/(e^-w - 1)) / w. O artigo imprime sigma_D = -0,30 h,
impossivel; o padrao usa |sigma_D| (gpd= troca).

A CRPC de uma realizacao: os pares (I, D_C) com falha ate crpc_dmax, a maior I sem falha (I_nf; as menores nao sao
feitas) e o FS0. Ajuste do artigo (Eq. 9): D_C = a I^b + c, minimos quadrados (para cada b numa grade, a e c lineares),
com 3 pontos ou mais; I_lb = a I em que o ajuste da D_C = dmax, limitada a [I_nf, a menor I com falha]. Com menos de 3
pontos (ou ajuste=interp), a interpolacao linear de ln D_C em ln I entre os pontos, constante acima da maior I, e I_lb
pela interpolacao linear de D_C entre (I_nf, dmax) e a menor I com falha. Acima de dmax a curva nao e extrapolada: a
falha exige D >= D_C(I) e I >= I_lb, sem teto em D.
P_1 = P(I >= I_lb, D >= D_C(I)) = int_{F_I(I_lb)}^1 [1 - dC/du(u, F_D(D_C(F_I^-1(u))))] du (Simpson com 20 000
intervalos; dC/du e a condicional da copula). P_F1 = media de P_1 nas realizacoes com FS0 >= 1 (Eq. 13), CoV pela Eq. 14
e o IC de 95% normal (Eq. 8).

Uso: python3 -I crpc_liuwang.py <pasta ou arquivos crpc/mc> [gpd=xiI,muI,sI,xiD,muD,sD,omega] [dmax=200]
                                 [ajuste=potencia|interp] [lambda=3.88] [--json arq]
"""
import csv
import glob
import json
import math
import os
import sys

GPD_PADRAO = (-0.04, 3.44, 5.13, 0.43, 5.15, 0.30, -2.38)


def cdf_gpd(x, xi, mu, s):
    if x <= mu:
        return 0.0
    z = (x - mu) / s
    if abs(xi) < 1e-12:
        return 1.0 - math.exp(-z)
    t = 1.0 + xi * z
    if t <= 0.0:  # acima do limite superior (xi < 0)
        return 1.0
    return 1.0 - t ** (-1.0 / xi)


def inv_gpd(u, xi, mu, s):
    return mu - s * math.log1p(-u) if abs(xi) < 1e-12 else mu + s * ((1.0 - u) ** (-xi) - 1.0) / xi


def cond_frank(v, u, w):
    """dC/du (u, v) da copula de Frank: P(V <= v | U = u)"""
    a, b, c = math.expm1(-w * u), math.expm1(-w * v), math.expm1(-w)
    return math.exp(-w * u) * b / (c + a * b)


def ajusta_potencia(pts):
    """D_C = a I^b + c por minimos quadrados: b numa grade de -6 a -0,01 (refinada), a e c lineares; (a, b, c, sse)"""
    melhor = None

    def lin(b):
        xs = [i ** b for i, _ in pts]
        ys = [d for _, d in pts]
        n = len(xs)
        mx, my = sum(xs) / n, sum(ys) / n
        sxx = sum((x - mx) ** 2 for x in xs)
        if sxx <= 0:
            return None
        a = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sxx
        c = my - a * mx
        sse = sum((a * x + c - y) ** 2 for x, y in zip(xs, ys))
        return a, b, c, sse

    grade = [-6.0 + k * 0.01 for k in range(600)]
    for b in grade:
        r = lin(b)
        if r and (melhor is None or r[3] < melhor[3]):
            melhor = r
    if melhor is None:
        return None
    b0 = melhor[1]
    for k in range(-100, 101):
        r = lin(b0 + k * 1e-4)
        if r and r[3] < melhor[3]:
            melhor = r
    return melhor


class Curva:
    """a CRPC de uma realizacao: D_C(I) e I_lb"""

    def __init__(self, pontos, inf, dmax, ajuste):
        self.p = sorted(pontos)  # (I, D_C) com falha
        self.inf = inf  # maior I sem falha (None: todas as I feitas falharam)
        self.dmax = dmax
        self.fit = None
        if ajuste == "potencia" and len(self.p) >= 3:
            self.fit = ajusta_potencia(self.p)
        imin = self.p[0][0]
        if self.fit:
            a, b, c, _ = self.fit
            # I_lb: D_C(I) = dmax, entre I_nf e a menor I com falha
            lo = self.inf if self.inf is not None else 1e-6
            f = lambda i: a * i ** b + c - dmax
            if f(imin) >= 0:
                self.ilb = imin
            elif f(lo) <= 0:
                self.ilb = lo
            else:
                x0, x1 = lo, imin
                for _ in range(200):
                    xm = 0.5 * (x0 + x1)
                    if f(xm) > 0:
                        x0 = xm
                    else:
                        x1 = xm
                self.ilb = 0.5 * (x0 + x1)
        else:
            if self.inf is None:
                self.ilb = imin
            else:
                self.ilb = self.inf  # D_C linear em I de (I_nf, dmax) ate a menor I com falha

    def dc(self, i):
        if i < self.ilb:
            return math.inf
        if self.fit:
            a, b, c, _ = self.fit
            return max(a * i ** b + c, 0.0)
        p = self.p
        if i >= p[-1][0]:
            return p[-1][1]
        if i <= p[0][0]:
            if self.inf is None:
                return p[0][1]
            # entre I_nf (D_C = dmax) e a menor I com falha
            t = (i - self.inf) / (p[0][0] - self.inf)
            return self.dmax + t * (p[0][1] - self.dmax)
        for (i0, d0), (i1, d1) in zip(p, p[1:]):
            if i0 <= i <= i1:
                t = (math.log(i) - math.log(i0)) / (math.log(i1) - math.log(i0))
                return math.exp(math.log(d0) + t * (math.log(d1) - math.log(d0)))
        return p[-1][1]


def p1(curva, gpd, n=20000):
    xiI, muI, sI, xiD, muD, sD, w = gpd
    ulb = cdf_gpd(curva.ilb, xiI, muI, sI)
    if ulb >= 1.0:
        return 0.0
    h = (1.0 - ulb) / n

    def f(u):
        if u >= 1.0:
            u = 1.0 - 1e-15
        i = inv_gpd(u, xiI, muI, sI)
        d = curva.dc(i)
        if math.isinf(d):
            return 0.0
        v = cdf_gpd(d, xiD, muD, sD)
        return 1.0 - cond_frank(v, u, w)

    s = f(ulb) + f(1.0 - 1e-12)
    for k in range(1, n):
        s += (4 if k % 2 else 2) * f(ulb + k * h)
    return s * h / 3.0


def le_crpc(arqs):
    """realizacao -> {'fs0', 'pontos', 'inf', 'q0', 'erro'}"""
    reais = {}
    for a in arqs:
        for l in csv.DictReader(open(a)):
            r = reais.setdefault(int(l["real"]), {"fs0": None, "pontos": [], "inf": None, "q0": float(l["q0_mmh"]),
                                                   "erro": None, "srm": 0, "s": 0.0})
            r["srm"] += int(l["n_srm"])
            r["s"] += float(l["s_transiente"]) + float(l["s_fator"])
            if l["falha"] == "descartada":
                r["fs0"] = float(l["FS0"])
            elif l["falha"].startswith("erro"):
                r["erro"] = l["falha"]
            else:
                r["fs0"] = float(l["FS0"])
                if l["falha"] == "sim":
                    r["pontos"].append((float(l["I"]), float(l["DC_h"])))
                else:  # a maior I sem falha (o caso deterministico roda uma I por arquivo)
                    r["inf"] = max(r["inf"] or 0.0, float(l["I"]))
    return reais


def resumo_crpc(arqs, gpd, dmax, ajuste):
    reais = le_crpc(arqs)
    ef, desc, err, p1s, linhas, srm = [], 0, 0, [], [], 0
    for n in sorted(reais):
        r = reais[n]
        srm += r["srm"]
        if r["erro"]:
            err += 1
            continue
        if r["fs0"] is None:
            continue
        if r["fs0"] < 1.0:
            desc += 1
            continue
        if r["pontos"]:
            c = Curva(r["pontos"], r["inf"], dmax, ajuste)
            v = p1(c, gpd)
            fit = c.fit
        else:
            c, v, fit = None, 0.0, None
        p1s.append(v)
        linhas.append({"real": n, "FS0": r["fs0"], "q0": r["q0"], "P1": v, "I_lb": c.ilb if c else None,
                       "ajuste": fit[:3] if fit else None, "pontos": r["pontos"], "I_nf": r["inf"]})
    N = len(p1s)
    out = {"N": N, "descartadas": desc, "erros": err, "srm": srm, "realizacoes": linhas}
    if N:
        m = sum(p1s) / N
        out["P_F1"] = m
        if N > 1 and m > 0:
            sd = math.sqrt(sum((x - m) ** 2 for x in p1s) / (N - 1))
            out["CoV"] = sd / (math.sqrt(N) * m)
            out["IC95"] = [m - 1.96 * sd / math.sqrt(N), m + 1.96 * sd / math.sqrt(N)]
    return out


def resumo_mc(pasta, gpd=None):
    """MC direto: FS em t = 0 e t = D por realizacao (<pasta>/mc_*_*.csv); descarta FS0 < 1 (Eqs. 6 a 8)"""
    fs = {}
    for a in glob.glob(os.path.join(pasta, "*_[0-9]*_[0-9]*.csv")):
        b = os.path.basename(a)
        if "_mapa_" in b or "_chuva_" in b or "_crpc_" in b:
            continue
        for l in csv.DictReader(open(a)):
            n = int(l["real"])
            if n < 0:
                continue
            fs.setdefault(n, {})[float(l["t_h"])] = (float(l["FS_min"]), l["situacao"])
    N = NF = desc = perd = 0
    for n, d in fs.items():
        if len(d) < 2 or any(s != "ok" and not (s == "fator" and v == 0.0) for v, s in d.values()):
            perd += 1
            continue
        t0 = min(d)
        if d[t0][0] < 1.0:
            desc += 1
            continue
        N += 1
        NF += min(v for v, _ in d.values()) < 1.0
    out = {"NT": len(fs), "N": N, "falhas": NF, "descartadas": desc, "perdidas": perd}
    if N:
        p = NF / N
        out["P_F1"] = p
        if NF:
            cov = math.sqrt((1 - p) / (N * p))
            out["CoV"] = cov
            out["IC95"] = [p - 1.96 * p * cov, p + 1.96 * p * cov]
    return out


def main(argv):
    opts = {"gpd": ",".join(str(x) for x in GPD_PADRAO), "dmax": "200", "ajuste": "potencia", "lambda": "3.88"}
    arqs, js = [], None
    it = iter(argv)
    for a in it:
        if a == "--json":
            js = next(it)
        elif "=" in a and a.split("=", 1)[0] in opts:
            k, v = a.split("=", 1)
            opts[k] = v
        else:
            arqs.append(a)
    gpd = tuple(float(x) for x in opts["gpd"].split(","))
    dmax, lam = float(opts["dmax"]), float(opts["lambda"])
    crpcs, pastas = [], []
    for a in arqs:
        if os.path.isdir(a):
            crpcs += sorted(glob.glob(os.path.join(a, "*_crpc_*.csv")))
            if glob.glob(os.path.join(a, "*_chuva_*.csv")):
                pastas.append(a)
        else:
            crpcs.append(a)
    res = {}
    if crpcs:
        r = resumo_crpc(crpcs, gpd, dmax, opts["ajuste"])
        res["crpc"] = r
        print(f"CRPC ({opts['ajuste']}): N = {r['N']} efetivas, {r['descartadas']} descartadas (FS0 < 1), "
              f"{r['erros']} com erro, {r['srm']} SRM")
        for l in r["realizacoes"][:40]:
            fit = l["ajuste"]
            print(f"  real {l['real']:4d}: FS0 {l['FS0']:.4f}, q0 {l['q0']:.4f} mm/h, pontos "
                  + " ".join(f"({i:g}; {d:.1f})" for i, d in sorted(l["pontos"]))
                  + (f", I_nf {l['I_nf']:g}" if l["I_nf"] else "")
                  + (f", D_C = {fit[0]:.4g} I^{fit[1]:.3f} + {fit[2]:.2f}" if fit else "")
                  + (f", I_lb {l['I_lb']:.2f}" if l["I_lb"] else "") + f": P_1 = {l['P1']:.3e}")
        if "P_F1" in r:
            pf = r["P_F1"]
            print(f"  P_F1 = {pf:.3e}" + (f", CoV {r['CoV']:.3f}, IC95 [{r['IC95'][0]:.3e}, {r['IC95'][1]:.3e}]"
                                          if "CoV" in r else "") + f"; P_FA = 1 - exp(-lambda P_F1) = "
                  f"{1 - math.exp(-lam * pf):.3e}")
    for p in pastas:
        r = resumo_mc(p)
        res.setdefault("mc", {})[p] = r
        print(f"MC direto {p}: {r['NT']} realizacoes, {r['N']} efetivas, {r['descartadas']} descartadas, "
              f"{r['perdidas']} perdidas, {r['falhas']} falhas"
              + (f": P_F1 = {r['P_F1']:.3e}" if "P_F1" in r else "")
              + (f", CoV {r['CoV']:.3f}, IC95 [{r['IC95'][0]:.3e}, {r['IC95'][1]:.3e}]" if "CoV" in r else ""))
    if js:
        json.dump(res, open(js, "w"), indent=1)


if __name__ == "__main__":
    main(sys.argv[1:])
