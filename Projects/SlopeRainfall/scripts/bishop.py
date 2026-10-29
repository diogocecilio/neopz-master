"""FS de Bishop simplificado com busca de círculos no talude de Rahardjo et al. (2007), para conferir o FS por redução
de resistência: poro-pressão inicial (NA hw abaixo do pé, inclinado a 7 graus, sucção limitada a smax) ou o campo p_w
de uma rodada (snaps_pw.csv). Resistência: tau = c' + (sigma_n - u+) tan phi' + s tan phi_b, com u+ = max(u, 0) e
s a sucção.

Uso: python3 bishop.py --c 0 --phi 26 --phib 17.3 --smax 50                      (campo inicial; Tang: 1,298)
     python3 bishop.py --c 0 --phi 26 --phib 17.3 --rodada <pasta> --t 0 25.6   (campos da rodada)
"""
import argparse
import os

import numpy as np


class Talude:
    def __init__(self, hs, alpha, hw):
        self.hs, self.hw = hs, hw
        self.ta = np.tan(np.radians(alpha))
        self.xb = 3 * hs
        self.xc = self.xb + hs / self.ta
        self.L = self.xc + 3 * hs

    def terreno(self, x):
        return np.where(x <= self.xb, 3 * self.hs,
                        np.where(x >= self.xc, 2 * self.hs, 3 * self.hs - (x - self.xb) * self.ta))

    def na(self, x):
        return (2 * self.hs - self.hw) - np.tan(np.radians(7.0)) * (x - self.xc)


def bishop(tal, xo, yo, R, pw, c, phi, phib, gamma, nfatias=50):
    """FS do círculo (xo, yo, R) ou None se o círculo não corta o talude em um único arco"""
    xs = np.linspace(max(0.0, xo - R), min(tal.L, xo + R), 400)
    ya = yo - np.sqrt(np.maximum(R * R - (xs - xo) ** 2, 0.0))
    abaixo = (ya < tal.terreno(xs) - 1e-6) & (np.abs(xs - xo) < R)
    idx = np.where(abaixo)[0]
    if len(idx) < 5 or idx[-1] - idx[0] + 1 != len(idx):
        return None
    xl, xr = xs[idx[0]], xs[idx[-1]]
    if xl <= 0.5 or xr >= tal.L - 0.5:
        return None
    xe = np.linspace(xl, xr, nfatias + 1)
    xm = 0.5 * (xe[1:] + xe[:-1])
    b = np.diff(xe)
    yb = yo - np.sqrt(np.maximum(R * R - (xm - xo) ** 2, 0.0))
    h = tal.terreno(xm) - yb
    if yb.min() < 0.0 or h.max() < 0.3:
        return None
    u = pw(xm, yb)
    if np.isnan(u).any():
        return None
    W = gamma * np.maximum(h, 0.0) * b
    sa = (xo - xm) / R  # seno da inclinação da base (positivo onde a base desce no sentido do escorregamento, +x)
    ca = np.sqrt(np.maximum(1 - sa * sa, 1e-6))
    tf, tb = np.tan(np.radians(phi)), np.tan(np.radians(phib))
    instab = np.sum(W * sa)
    if instab <= 0:
        return None
    resist = c * b + np.maximum(-u, 0.0) * tb * b + (W - np.maximum(u, 0.0) * b) * tf
    F = 1.5
    for _ in range(100):
        Fn = np.sum(resist / np.maximum(ca + sa * tf / F, 0.2)) / instab
        if abs(Fn - F) < 1e-5:
            break
        F = Fn
    return F, xl, xr


def busca(tal, pw, a):
    melhor = (np.inf, None)
    for xo in np.arange(tal.xb - 4.0, tal.xc + 17.0, a.passo):
        for yo in np.arange(3 * tal.hs, 3 * tal.hs + 3.6 * tal.hs, a.passo):
            for R in np.arange(0.4 * tal.hs, 4.2 * tal.hs, a.passo / 2):
                r = bishop(tal, xo, yo, R, pw, a.c, a.phi, a.phib, a.gamma)
                if r and r[0] < melhor[0]:
                    melhor = (r[0], (xo, yo, R, round(r[1], 2), round(r[2], 2)))
    return melhor


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--c", type=float, default=10.0)
    p.add_argument("--phi", type=float, default=26.0)
    p.add_argument("--phib", type=float, default=26.0)
    p.add_argument("--gamma", type=float, default=20.0)
    p.add_argument("--hs", type=float, default=10.0)
    p.add_argument("--alpha", type=float, default=45.0)
    p.add_argument("--hw", type=float, default=5.0)
    p.add_argument("--smax", type=float, default=75.0, help="limite da sucção inicial (kPa)")
    p.add_argument("--rodada", help="pasta de uma rodada do SlopeRainfall (usa snaps_pw.csv)")
    p.add_argument("--t", type=float, nargs="*", default=[0.0], help="instantes (h) dos campos da rodada")
    p.add_argument("--passo", type=float, default=1.0, help="passo da grade de centros (m)")
    a = p.parse_args()
    tal = Talude(a.hs, a.alpha, a.hw)
    if not a.rodada:
        pw = lambda x, y: np.maximum(9.81 * (tal.na(x) - y), -a.smax)
        F, circ = busca(tal, pw, a)
        print(f"campo inicial (smax = {a.smax:g} kPa): FS Bishop = {F:.4f}, círculo (xo, yo, R, x1, x2) = {circ}")
        return
    import matplotlib.tri as mtri
    V = np.loadtxt(os.path.join(a.rodada, "vertices.csv"), delimiter=",", skiprows=1)
    T = np.loadtxt(os.path.join(a.rodada, "triangles.csv"), delimiter=",", skiprows=1, dtype=int)
    S = np.atleast_2d(np.loadtxt(os.path.join(a.rodada, "snaps_pw.csv"), delimiter=","))
    tri = mtri.Triangulation(V[:, 0], V[:, 1], T)
    for t in a.t:
        k = int(np.argmin(np.abs(S[:, 0] - t)))
        interp = mtri.LinearTriInterpolator(tri, S[k, 1:])
        pw = lambda x, y: np.asarray(interp(x, y).filled(np.nan))
        F, circ = busca(tal, pw, a)
        print(f"{a.rodada}, t = {S[k, 0]:g} h: FS Bishop = {F:.4f}, círculo (xo, yo, R, x1, x2) = {circ}", flush=True)


if __name__ == "__main__":
    main()
