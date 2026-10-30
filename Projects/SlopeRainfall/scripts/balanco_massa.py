"""Balanço global de massa de água e de ar de uma rodada do SlopeRainfall com balanco=1, nos instantes dos snapshots.

Armazenamento nos vértices com as massas condensadas dos triângulos lineares (área/3 por vértice), como o termo de
armazenamento do material: água n S_w (V_w) e ar n (1 - S_w) (p_a + p_atm)/p_atm (V_a, volume à pressão atmosférica).
Fluxos de hist.csv: infiltração pela superfície (infil_m3s), água pelas laterais e pela base (water_bnd_m3s) e ar pelos
contornos com p_a prescrita (air_in_m3s). Só para o modelo desacoplado (rígido): no acoplado o armazenamento tem também a
deformação volumétrica, que não está nos snapshots.

Uso: python3 balanco_massa.py <pasta da rodada> [--n 0.45] [--a 50] [--nfx 1] [--mfx 1]
"""
import argparse
import os

import numpy as np

PATM = 101.325


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("pasta")
    p.add_argument("--n", type=float, default=0.45, help="porosidade")
    p.add_argument("--a", type=float, default=50.0, help="Fredlund-Xing a (kPa)")
    p.add_argument("--nfx", type=float, default=1.0)
    p.add_argument("--mfx", type=float, default=1.0)
    a = p.parse_args()
    V = np.loadtxt(os.path.join(a.pasta, "vertices.csv"), delimiter=",", skiprows=1)
    T = np.loadtxt(os.path.join(a.pasta, "triangles.csv"), delimiter=",", skiprows=1, dtype=int)
    area = 0.5 * np.abs((V[T[:, 1], 0] - V[T[:, 0], 0]) * (V[T[:, 2], 1] - V[T[:, 0], 1]) -
                        (V[T[:, 2], 0] - V[T[:, 0], 0]) * (V[T[:, 1], 1] - V[T[:, 0], 1]))
    M = np.zeros(len(V))
    for k in range(3):
        np.add.at(M, T[:, k], area / 3.0)
    PW = np.atleast_2d(np.loadtxt(os.path.join(a.pasta, "snaps_pw.csv"), delimiter=","))
    PA = np.atleast_2d(np.loadtxt(os.path.join(a.pasta, "snaps_pa.csv"), delimiter=","))
    H = np.genfromtxt(os.path.join(a.pasta, "hist.csv"), delimiter=",", names=True)

    def S(pc):
        pc = np.maximum(pc, 0.0)
        return np.log(np.e + (pc / a.a) ** a.nfx) ** (-a.mfx)

    def armazenado(i):
        pw, pa = PW[i, 1:], PA[i, 1:]
        s = S(pa - pw)
        return np.sum(M * a.n * s), np.sum(M * a.n * (1.0 - s) * (pa + PATM) / PATM)

    tw0, ta0 = armazenado(0)
    t_h = H["t_h"]
    dt_s = H["dt_h"] * 3600.0
    cw = np.cumsum((H["infil_m3s"] + H["water_bnd_m3s"]) * dt_s)
    ca = np.cumsum(H["air_in_m3s"] * dt_s)
    print("t_h, dV_agua (m3/m), entrada de agua (m3/m), erro agua (%), dV_ar (m3/m), entrada de ar (m3/m), "
          "erro ar (% de V_ar inicial)")
    for i in range(1, len(PW)):
        t = PW[i, 0]
        k = np.argmin(np.abs(t_h - t))
        if abs(t_h[k] - t) > 1e-6:
            continue
        w, ar = armazenado(i)
        dw, dar = w - tw0, ar - ta0
        ew = 100.0 * (dw - cw[k]) / max(abs(cw[k]), 1e-12)
        ea = 100.0 * (dar - ca[k]) / ta0
        print(f"{t:g}, {dw:.4f}, {cw[k]:.4f}, {ew:.2f}, {dar:.4f}, {ca[k]:.4f}, {ea:.3f}")


if __name__ == "__main__":
    main()
