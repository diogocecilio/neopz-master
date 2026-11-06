#!/usr/bin/env python3
# Equilibrio limite em superficies circulares nos dois taludes de Cho (2010) com a geometria reconstruida das Figs. 2 e
# 6 de Liao et al. (2025): referencia independente para o FS do SRM. Bishop simplificado (para phi = 0 e o metodo sueco,
# exato para o circulo: FS = c L R / sum(W x)). Busca em grade de centros e raios, so circulos dentro do dominio
# (0 <= x <= 30, fundo z >= 0, a base rigida), depois refinamento local. Uso: python3 equilibrio_limite.py [cho1|cho2]
import sys
import numpy as np

CASOS = {
    # L, xc, zc, xp, zp, c (kPa), phi (graus), gamma (kN/m3)
    "cho1": (30., 10., 10., 20., 5., 23., 0., 20.),
    "cho2": (30., 10., 15., 20., 5., 10., 30., 20.),
}


def superficie(x, L, xc, zc, xp, zp):
    return np.where(x <= xc, zc, np.where(x >= xp, zp, zc + (zp - zc) * (x - xc) / (xp - xc)))


def fs_circulo(x0, y0, R, caso, n=200):
    L, xc, zc, xp, zp, c, phig, g = caso
    if y0 - R < -1e-9:  # abaixo da base rigida
        return np.inf
    a, b = max(0., x0 - R), min(L, x0 + R)
    if b - a < 1e-6:
        return np.inf
    xs = np.linspace(a, b, 2001)
    arco = y0 - np.sqrt(np.maximum(R * R - (xs - x0) ** 2, 0.))
    d = superficie(xs, L, xc, zc, xp, zp) - arco
    dentro = d > 0
    if not dentro.any():
        return np.inf
    i0, i1 = np.argmax(dentro), len(dentro) - 1 - np.argmax(dentro[::-1])
    if dentro[0] or dentro[-1]:  # a massa toca uma lateral, ou o circulo nao sai pela superficie
        return np.inf
    if np.count_nonzero(dentro[i0:i1 + 1]) != i1 - i0 + 1:  # duas massas: descarta
        return np.inf
    # pontos de entrada e saida por interpolacao linear de d
    xa = xs[i0 - 1] + (xs[i0] - xs[i0 - 1]) * (-d[i0 - 1]) / (d[i0] - d[i0 - 1])
    xb = xs[i1] + (xs[i1 + 1] - xs[i1]) * d[i1] / (d[i1] - d[i1 + 1])
    if xb - xa < 0.5:
        return np.inf
    e = np.linspace(xa, xb, n + 1)
    xm, bw = 0.5 * (e[1:] + e[:-1]), np.diff(e)
    raiz = np.sqrt(np.maximum(R * R - (xm - x0) ** 2, 1e-12))
    h = superficie(xm, L, xc, zc, xp, zp) - (y0 - raiz)
    W = g * bw * np.maximum(h, 0.)
    sena, cosa = (x0 - xm) / R, raiz / R  # alfa > 0 onde a base desce no sentido do movimento (para a direita)
    motor = np.sum(W * sena)
    if motor <= 0.:
        return np.inf
    tphi = np.tan(np.radians(phig))
    F = 1.
    for _ in range(100):
        m = cosa + sena * tphi / F
        if np.any(m <= 0.05):
            return np.inf
        Fn = np.sum((c * bw + W * tphi) / m) / motor
        if abs(Fn - F) < 1e-10:
            break
        F = Fn
    return Fn


def busca(nome):
    caso = CASOS[nome]
    L, xc, zc = caso[0], caso[1], caso[2]
    melhor = (np.inf, None)
    for x0 in np.arange(5., 31., 1.):
        for y0 in np.arange(zc - 4., zc + 31., 1.):
            for R in np.arange(2., y0 + 1e-9, 0.25):
                F = fs_circulo(x0, y0, R, caso, 100)
                if F < melhor[0]:
                    melhor = (F, (x0, y0, R))
    # refinamento local (coordenadas alternadas com passo decrescente)
    F, (x0, y0, R) = melhor
    passo = 0.5
    while passo > 1e-3:
        mudou = False
        for dx, dy, dr in [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1), (1, 1, 1), (-1, -1, -1)]:
            t = (x0 + dx * passo, y0 + dy * passo, R + dr * passo)
            Ft = fs_circulo(*t, caso, 400)
            if Ft < F - 1e-12:
                F, (x0, y0, R), mudou = Ft, t, True
        if not mudou:
            passo /= 2.
    return F, (x0, y0, R)


if __name__ == "__main__":
    for nome in (sys.argv[1:] or ["cho1", "cho2"]):
        F, (x0, y0, R) = busca(nome)
        print(f"{nome}: FS = {F:.4f}, centro = ({x0:.3f}, {y0:.3f}), raio = {R:.3f}, ponto mais baixo z = {y0 - R:.3f}")
