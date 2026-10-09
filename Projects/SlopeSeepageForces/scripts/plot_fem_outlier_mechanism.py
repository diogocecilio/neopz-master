#!/usr/bin/env python3
"""Picture of the FEM collapse mechanism of an outlier run (scripts/run_fem_outliers.sh with vtk=): the plastic strain
sqrt(J2(eps_p)) of the last accepted state of a refinement cycle (log colour scale, from the <prefix>_GI_ref<k>.scal_vec.0.vtk
written by SlopeAnalysis::PostPlasticity) and the displacement field of that state (arrows, direction only), with the
slope outline, the stability box and the log-spiral of the limit analysis of the same case (from an 'la' log of
results/fem/outliers/) overlaid, with the velocity directions of its rigid rotation about C at the same points (green
arrows). Prints the mean direction (degrees above the horizontal) of the FEM displacement and of the rigid rotation in
the toe region of the mechanism, and how much of the area of the limit-analysis block (between the spiral and the
ground) is plastic in the FEM (cells with sqrt(J2) above 10 %, 1 % of the max): a rigid-block mechanism shows as a
thin band along the spiral, a deforming one fills the block. NeoPZ coordinates (y up, crest edge O at the origin).

    python3 scripts/plot_fem_outlier_mechanism.py <scal_vec vtk> <la log> <png> [title]
"""
import math
import re
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PolyCollection
from matplotlib.colors import LogNorm


def read_vtk(path):
    L = open(path).read().split('\n')
    i = next(k for k, l in enumerate(L) if l.startswith('POINTS'))
    n = int(L[i].split()[1])
    P = np.array([list(map(float, L[i + 1 + k].split()))[:2] for k in range(n)])
    i = next(k for k, l in enumerate(L) if l.startswith('CELLS'))
    nc = int(L[i].split()[1])
    cells = [list(map(int, L[i + 1 + k].split()))[1:] for k in range(nc)]

    def scal(name):
        j = next(k for k, l in enumerate(L) if l.startswith('SCALARS ' + name))
        return np.array([float(L[j + 2 + k]) for k in range(n)])

    j = next(k for k, l in enumerate(L) if l.startswith('VECTORS Displacement'))
    U = np.array([list(map(float, L[j + 1 + k].split()))[:2] for k in range(n)])
    return P, cells, scal('StrainPlasticJ2'), U


def read_la(path):
    """H, beta, phi, A, C, theta1, theta2 (paper coordinates) of an 'la' log"""
    txt = open(path).read()
    g = {}
    m = re.search(r'theta1 ([-\d.e+]+), theta2 ([-\d.e+]+)', txt)
    g['theta1'], g['theta2'] = float(m.group(1)), float(m.group(2))
    m = re.search(r'A = \(([-\d.e+]+), ([-\d.e+]+)\), B = \(([-\d.e+]+), ([-\d.e+]+)\), C = \(([-\d.e+]+), ([-\d.e+]+)\)', txt)
    g['A'] = (float(m.group(1)), float(m.group(2)))
    g['B'] = (float(m.group(3)), float(m.group(4)))
    g['C'] = (float(m.group(5)), float(m.group(6)))
    m = re.search(r'Gamma = ([-\d.e+]+)', txt)
    g['Gamma'] = float(m.group(1))
    return g


def spiral(g, phi_deg):
    """points of the log-spiral r = r0 exp((theta - theta1) tan phi) about C: x = C_x - r cos theta, y_paper = C_y +
    r sin theta (checked on A and B of the la logs); returned in NeoPZ coordinates"""
    Cx, Cy = g['C']
    r0 = math.hypot(g['A'][0] - Cx, g['A'][1] - Cy)
    th = np.linspace(g['theta1'], g['theta2'], 300)
    r = r0 * np.exp((th - g['theta1']) * math.tan(math.radians(phi_deg)))
    return Cx - r * np.cos(th), -(Cy + r * np.sin(th))


def main():
    vtk, lalog, png = sys.argv[1:4]
    title = sys.argv[4] if len(sys.argv) > 4 else ''
    P, cells, J2, U = read_vtk(vtk)
    g = read_la(lalog)
    txt = open(lalog).read()
    H = float(re.search(r'H ([\d.]+) m', txt).group(1)) if re.search(r'H ([\d.]+) m', txt) else 5.
    beta = float(re.search(r'beta ([\d.]+) deg', txt).group(1))
    phi = float(re.search(r'phi ([\d.]+)', txt).group(1)) if re.search(r'phi ([\d.]+)', txt) else 30.
    xT = H / math.tan(math.radians(beta))
    fig, ax = plt.subplots(figsize=(11, 6.5))
    vmax = J2.max()
    cval = np.array([max(J2[c].max(), 1e-4 * vmax) for c in cells])
    polys = [P[c] for c in cells]
    pc = PolyCollection(polys, array=cval / vmax, cmap='magma_r', norm=LogNorm(1e-3, 1.), edgecolor='none')
    ax.add_collection(pc)
    cb = fig.colorbar(pc, ax=ax, shrink=0.8)
    cb.set_label('sqrt(J2(eps_p)) / max (last accepted state)')
    xs, ys = spiral(g, phi)
    Cx, Cy = g['C'][0], -g['C'][1]  # NeoPZ coordinates
    # displacement directions at points of the plastic zone, and the rigid rotation about C of the limit analysis
    # (counterclockwise in y-up coordinates: the crest moves down, the toe out) at the points inside its block
    sel = np.where(J2 > 1e-3 * vmax)[0]
    if len(sel):
        rng = np.random.default_rng(0)
        sel = rng.choice(sel, size=min(400, len(sel)), replace=False)
        nrm = np.hypot(U[sel, 0], U[sel, 1])
        ok = nrm > 0
        ax.quiver(P[sel[ok], 0], P[sel[ok], 1], U[sel[ok], 0] / nrm[ok], U[sel[ok], 1] / nrm[ok], color='tab:blue',
                  scale=35, width=0.002, alpha=0.8, label='FEM displacement direction')
        inside = (P[sel, 0] > xs.min()) & (P[sel, 0] < xs.max()) & (P[sel, 1] > np.interp(P[sel, 0], xs, ys)) & ok
        q = sel[inside]
        vx, vy = -(P[q, 1] - Cy), P[q, 0] - Cx
        vn = np.hypot(vx, vy)
        ax.quiver(P[q, 0], P[q, 1], vx / vn, vy / vn, color='green', scale=35, width=0.002, alpha=0.6,
                  label='rigid rotation of the limit analysis')
        # mean directions in the toe region of the mechanism
        for x0 in (xT - H, xT - 0.5 * H, xT):
            toe = sel[ok & (P[sel, 0] > x0)]
            if len(toe):
                ang = np.degrees(np.arctan2(U[toe, 1], U[toe, 0]))
                angla = np.degrees(np.arctan2(P[toe, 0] - Cx, -(P[toe, 1] - Cy)))
                print(f'x > {x0:.2f} m ({len(toe)} points of the plastic zone): mean direction of the FEM displacement '
                      f'{ang.mean():.1f} deg (std {ang.std():.1f}), of the rigid rotation {angla.mean():.1f} deg '
                      f'(std {angla.std():.1f}), above the horizontal')
    # plastic fraction of the area of the limit-analysis block (cells by their centroid)
    cen = np.array([P[c].mean(axis=0) for c in cells])
    area = np.array([0.5 * abs(np.dot(P[c][:, 0], np.roll(P[c][:, 1], -1)) - np.dot(P[c][:, 1], np.roll(P[c][:, 0], -1)))
                     for c in cells])
    inb = (cen[:, 0] > xs.min()) & (cen[:, 0] < xs.max()) & (cen[:, 1] > np.interp(cen[:, 0], xs, ys))
    ablock = area[inb].sum()
    if ablock > 0:
        print(f'limit-analysis block: area {ablock:.2f} m2 ({ablock / H**2:.3f} H^2); FEM plastic area inside it: '
              + ', '.join(f'{100 * area[inb & (cval >= fr * vmax)].sum() / ablock:.0f} % above {100 * fr:g} % of max'
                          for fr in (0.1, 0.01)))
        # thickness of the band normal to the spiral: the plastic area divided by the length of the spiral
        length = np.hypot(np.diff(xs), np.diff(ys)).sum()
        print(f'spiral length {length:.2f} m; equivalent thickness of the plastic zone (> 10 % / 1 % of max, whole mesh): '
              + ', '.join(f'{area[cval >= fr * vmax].sum() / length:.2f} m' for fr in (0.1, 0.01)))
    # slope outline and box
    ax.plot(xs, ys, 'g-', lw=2, label=f"limit analysis log-spiral (Gamma_LA {g['Gamma']:.3f})")
    ax.plot([g['C'][0]], [-g['C'][1]], 'g+', ms=12, mew=2)
    xmin, xmax, ymin = P[:, 0].min(), P[:, 0].max(), P[:, 1].min()
    ax.plot([xmin, 0, xT, xmax], [0, 0, -H, -H], 'k-', lw=1.5)
    ax.plot([xmin, xmin, xmax, xmax], [0, ymin, ymin, -H], 'k--', lw=0.8)
    ax.set_xlim(min(-1.5 * H, g['A'][0] - 0.3 * H), max(xT + 1.6 * H, g['B'][0] + 0.3 * H))
    ax.set_ylim(-H - 1.2 * H, 0.4 * H)
    ax.set_aspect('equal')
    ax.set_xlabel('x (m), crest edge O at 0')
    ax.set_ylabel('y (m), up')
    ax.set_title(title)
    ax.legend(loc='lower left', fontsize=9)
    fig.tight_layout()
    fig.savefig(png, dpi=130)
    print('written', png)


if __name__ == '__main__':
    main()
