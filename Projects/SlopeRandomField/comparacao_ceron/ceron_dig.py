#!/usr/bin/env python3
"""ceron_dig.py: os dados de Vargas Ceron et al. (2025), Secao 6, para as figuras comparativas, num dicionario so:

  tab[familia]  -> linhas das Tabelas 3 a 6: x (o valor do parametro), mu, sigma, cov, pf (fracoes) e det (Gamma do
                   deterministico, quando a tabela o tem);
  pdfs[familia] -> {valor: (x, y)}: as densidades das Figs. 27, 31a, 32a e 33a extraidas do PDF; sem elas, as log-normais
                   com o mu e o sigma da tabela (o artigo desenha as log-normais ajustadas, Secao 6);
  fig24, fig25a, fig25b, fig26 -> a superficie critica (coordenadas do modelo), o histograma e a densidade, a
                   convergencia de P_f e os coeficientes de Spearman, quando extraidos.

As tabelas vem do texto do artigo (TAB_TEXTO, conferidas com a extracao em digit/out/tabelas quando ela existe)."""

import csv
import glob
import json
import math
import os

# Tabelas 3 a 6 do artigo: (valor, mu, sigma, CoV %, P_f %) e, nas 5 e 6, o deterministico antes de mu
TAB3 = [  # CoV de k_v, c, phi (%)
    ((60, 30, 10), 1.353, 0.318, 23.5, 11.50), ((0, 30, 10), 1.329, 0.288, 21.7, 11.21),
    ((75, 30, 10), 1.375, 0.341, 24.8, 11.41), ((90, 30, 10), 1.378, 0.354, 25.7, 11.63),
    ((100, 30, 10), 1.391, 0.360, 25.9, 11.44), ((60, 10, 10), 1.367, 0.180, 13.2, 0.45),
    ((60, 50, 10), 1.324, 0.472, 35.7, 25.78), ((60, 70, 10), 1.296, 0.620, 47.9, 36.34),
    ((60, 30, 5), 1.351, 0.300, 22.2, 10.20), ((60, 30, 15), 1.357, 0.346, 25.5, 12.71),
    ((60, 30, 20), 1.361, 0.384, 28.2, 15.18)]
TAB4 = [(1.0, 1.353, 0.318, 23.5, 11.50), (1.5, 1.356, 0.352, 25.9, 13.94), (2.0, 1.361, 0.372, 27.3, 14.99),
        (5.0, 1.350, 0.404, 29.9, 18.93), (10.0, 1.348, 0.419, 31.1, 19.85), (20.0, 1.355, 0.430, 31.7, 20.23),
        (400.0, 1.348, 0.432, 32.1, 21.93)]
TAB5 = [(1.0, 1.336, 1.353, 0.318, 23.5, 11.50), (2.0, 1.533, 1.562, 0.390, 24.9, 3.90),
        (3.0, 1.674, 1.704, 0.432, 25.3, 1.83), (4.0, 1.783, 1.821, 0.474, 26.0, 1.00),
        (5.0, 1.872, 1.910, 0.504, 26.4, 0.63)]
TAB6 = [(0.5, 1.671, 1.670, 0.378, 22.6, 1.38), (0.6, 1.494, 1.502, 0.345, 22.9, 4.20),
        (0.7, 1.383, 1.393, 0.314, 22.5, 8.54), (0.8, 1.322, 1.338, 0.309, 23.1, 11.99),
        (0.9, 1.307, 1.325, 0.306, 23.1, 12.58), (1.0, 1.336, 1.353, 0.318, 23.5, 11.50)]
REF = {"covk": 60, "covc": 30, "covphi": 10}


def tabelas():
    tab = {"covk": [], "covc": [], "covphi": []}
    for (k, c, p), mu, sg, cov, pf in TAB3:
        for fam, v, idx in (("covk", k, 0), ("covc", c, 1), ("covphi", p, 2)):
            outros = [x for j, x in enumerate((k, c, p)) if j != idx]
            refs = [REF[f] for j, f in enumerate(("covk", "covc", "covphi")) if j != idx]
            if outros == refs:
                tab[fam].append({"x": v / 100., "mu": mu, "sigma": sg, "cov": cov / 100., "pf": pf / 100.})
    tab["s"] = [{"x": s, "mu": mu, "sigma": sg, "cov": cov / 100., "pf": pf / 100.} for s, mu, sg, cov, pf in TAB4]
    tab["alfa"] = [{"x": a, "det": d, "mu": mu, "sigma": sg, "cov": cov / 100., "pf": pf / 100.}
                   for a, d, mu, sg, cov, pf in TAB5]
    tab["hw"] = [{"x": h, "det": d, "mu": mu, "sigma": sg, "cov": cov / 100., "pf": pf / 100.}
                 for h, d, mu, sg, cov, pf in TAB6]
    for fam in tab:
        tab[fam].sort(key=lambda r: r["x"])
    return tab


def lognormal_da_tabela(mu, sigma, n=400):
    s2 = math.log(1. + (sigma / mu) ** 2)
    m = math.log(mu) - s2 / 2.
    s = math.sqrt(s2)
    lo, hi = math.exp(m - 4.5 * s), math.exp(m + 4.5 * s)
    xs = [lo + (hi - lo) * i / (n - 1) for i in range(n)]
    ys = [math.exp(-(math.log(x) - m) ** 2 / (2 * s2)) / (x * s * math.sqrt(2 * math.pi)) for x in xs]
    return xs, ys


def carrega(pasta_digit=None):
    tab = tabelas()
    cer = {"tab": tab, "pdfs": {}, "origem_pdfs": {}, "fig24": None, "fig25a": {}, "fig25b": None, "fig26": {}}
    for fam, linhas in tab.items():
        cer["pdfs"][fam] = {}
        for r in linhas:
            if fam == "s" and r["x"] > 20:
                continue  # s = 400 nao esta na Fig. 31a
            cer["pdfs"][fam][r["x"]] = lognormal_da_tabela(r["mu"], r["sigma"])
            cer["origem_pdfs"][(fam, r["x"])] = "tabela"
    if pasta_digit and os.path.isfile(os.path.join(pasta_digit, "ceron_dig.json")):
        with open(os.path.join(pasta_digit, "ceron_dig.json")) as f:
            d = json.load(f)
        for fam, curvas in d.get("pdfs", {}).items():
            for v, xy in curvas.items():
                cer["pdfs"].setdefault(fam, {})[float(v)] = (xy[0], xy[1])
                cer["origem_pdfs"][(fam, float(v))] = "extraida"
        for k in ("fig24", "fig25a", "fig25b", "fig26"):
            if d.get(k):
                cer[k] = d[k]
    return cer
