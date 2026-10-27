"""Resultados de Rahardjo et al. (2007) usados na comparação (FS por Bishop simplificado, SEEP/W + SLOPE/W).

Tabela 7 e Eq. 6 (com os parâmetros das Figs. 7a-c) são transcritas do artigo; as Figs. 5b, 8 e 10 foram
digitalizadas à mão a partir do PDF (incerteza da ordem de +-0,03 no FS)."""
import numpy as np

# Tabela 7: FS = E + F alpha (graus)
FS_INI = (3.153, -0.023)
TABELA7 = {("f10", 9.0): (3.006, -0.025), ("f10", 80.0): (2.830, -0.027), ("f10", 360.0): (1.495, -0.014),
           ("f50", 9.0): (3.000, -0.025), ("f50", 36.0): (1.863, -0.018), ("f50", 80.0): (1.908, -0.017),
           ("f100", 3.6): (3.007, -0.025), ("f100", 9.0): (2.910, -0.025), ("f100", 80.0): (1.455, -0.013)}

# Eq. 6: FS_min = A + B / (1 + exp(-(I_r - C)/D)), Figs. 7a, 7b e 7c (alpha = 45 graus, Hs = 10 m, Hw = 5 m)
EQ6 = {"f10": (0.81, 13.54, -118.0, -52.0), "f50": (0.82, 1.83, 10.7, -13.13), "f100": (1.45, 0.70, 4.75, -2.35)}


def fs_min_eq6(solo, I):
    A, B, C, D = EQ6[solo]
    return A + B / (1 + np.exp(-(np.asarray(I, float) - C) / D))


def fig5b():
    """{I_r: (t, FS)} para o solo f50,-5"""
    d = np.loadtxt(__file__.replace("ref_artigo.py", "ref/fig5b_rahardjo2007.csv"), delimiter=",", comments="#")
    return {I: (d[d[:, 0] == I, 1], d[d[:, 0] == I, 2]) for I in np.unique(d[:, 0])}


# Fig. 8 (Série C): f50,-5, alpha = 45 graus, Hs = 10 m
FIG8_HW = np.array([2.5, 5.0, 7.5, 10.0, 15.0])
FIG8 = {"ini": [1.80, 2.11, 2.33, 2.35, 2.35], 9.0: [1.43, 1.79, 2.00, 2.13, 2.20],
        36.0: [0.85, 1.04, 1.19, 1.26, 1.29], 80.0: [0.81, 0.82, 1.05, 1.15, 1.16]}

# Fig. 10 (Série D): f50,-5, alpha = 45 graus, Hw = 5 m
FIG10_HS = np.array([5.0, 10.0, 20.0, 40.0])
FIG10 = {"ini": [3.22, 2.11, 1.50, 1.12], 9.0: [2.48, 1.79, 1.37, 1.08],
         36.0: [1.20, 1.04, 0.97, 0.90], 80.0: [1.19, 0.82, 0.88, 0.83]}

# Fig. 9 (Série A, Hs = 10 m, Hw = 5 m): FS mínimo x alpha, lido da figura (incerteza ~ +-0,05). A Tabela 7 não
# reproduz a Fig. 9 em algumas linhas (p. ex. f10 com 80 mm/h: 2,11 na tabela e ~1,45 na figura em 26,6 graus).
FIG9_ALPHA = np.array([26.6, 33.7, 45.0, 63.4])
FIG9 = {("f10", 9.0): [2.35, 2.16, 1.88, 1.45], ("f10", 80.0): [1.45, 1.32, 1.12, 0.83],
        ("f10", 360.0): [1.15, 1.00, 0.82, 0.62],
        ("f50", 9.0): [2.27, 2.07, 1.79, 1.36], ("f50", 36.0): [1.43, 1.27, 1.04, 0.76],
        ("f50", 80.0): [1.15, 1.00, 0.82, 0.67],
        ("f100", 3.6): [2.37, 2.19, 1.89, 1.45], ("f100", 9.0): [2.12, 1.92, 1.55, 1.11],
        ("f100", 80.0): [2.05, 1.83, 1.45, 1.05]}
