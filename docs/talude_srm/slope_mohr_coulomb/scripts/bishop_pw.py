# Bishop simplified with pore pressure (effective stress) for the slope of TriGMesh, independent check of
# SlopeDrawdown. Pore pressure p+ = max(p, 0) interpolated from the u-p VTK (PorePressure point data).
# Usage: python3 bishop_pw.py <file.vtk | dry | crest | drained>
import sys
import numpy as np

c0, phi0, gsat, gw = 10.0, np.radians(30.0), 20.0, 10.0


def ground(x):
    return np.where(x <= 30, 40.0, np.where(x >= 40, 30.0, 40.0 - (x - 30.0)))


def read_vtk(fname):
    lines = open(fname).read().split('\n')
    i = next(k for k, l in enumerate(lines) if l.startswith('POINTS'))
    n = int(lines[i].split()[1])
    vals = []
    j = i + 1
    while len(vals) < 3 * n:
        vals += [float(v) for v in lines[j].split()]
        j += 1
    pts = np.array(vals).reshape(n, 3)[:, :2]
    i = next(k for k, l in enumerate(lines) if l.startswith('CELLS'))
    ncell = int(lines[i].split()[1])
    cells = [list(map(int, lines[i + 1 + k].split()))[1:] for k in range(ncell)]
    i = next(k for k, l in enumerate(lines) if l.startswith('SCALARS PorePressure'))
    vals = []
    j = i + 2
    while len(vals) < n:
        vals += [float(v) for v in lines[j].split()]
        j += 1
    tris = []
    for c in cells:  # graph mesh of NeoPZ: quadrilaterals (possibly degenerate) -> two triangles
        for t in ([c] if len(c) == 3 else [[c[0], c[1], c[2]], [c[0], c[2], c[3]]]):
            a, b, d = pts[t[0]], pts[t[1]], pts[t[2]]
            if abs((b[0] - a[0]) * (d[1] - a[1]) - (b[1] - a[1]) * (d[0] - a[0])) > 1e-12:
                tris.append(t)
    return pts, np.array(tris), np.array(vals)


class Field:
    """p+ on a 0.25 m grid (exact linear interpolation of the VTK triangles at the grid nodes), bilinear between"""
    def __init__(self, fname, h=0.25):
        pts, tris, p = read_vtk(fname)
        P = pts[tris]
        x0 = P[:, 0]
        inv = np.linalg.inv(np.stack([P[:, 1] - P[:, 0], P[:, 2] - P[:, 0]], axis=2))
        self.h, self.gx, self.gy = h, np.arange(0, 70 + h / 2, h), np.arange(0, 40 + h / 2, h)
        self.g = np.zeros((len(self.gx), len(self.gy)))
        for i, xx in enumerate(self.gx):
            for j, yy in enumerate(self.gy):
                if yy > ground(np.array([xx]))[0] + 1e-9:
                    continue
                xi = np.einsum('kij,kj->ki', inv, np.array([xx, yy]) - x0)
                ok = (xi[:, 0] >= -1e-9) & (xi[:, 1] >= -1e-9) & (xi.sum(1) <= 1 + 1e-9)
                t = np.argmax(ok)
                if ok[t]:
                    pv = p[tris[t]]
                    self.g[i, j] = pv[0] + (pv[1] - pv[0]) * xi[t, 0] + (pv[2] - pv[0]) * xi[t, 1]
        self.g = np.maximum(self.g, 0.0)

    def __call__(self, x, y):
        fx, fy = np.clip(x / self.h, 0, len(self.gx) - 1.001), np.clip(y / self.h, 0, len(self.gy) - 1.001)
        i, j = fx.astype(int), fy.astype(int)
        a, b = fx - i, fy - j
        g = self.g
        return (1 - a) * (1 - b) * g[i, j] + a * (1 - b) * g[i + 1, j] + (1 - a) * b * g[i, j + 1] + a * b * g[i + 1, j + 1]


def fs_circle(xc, yc, R, c, tanphi, pore, gam, n):
    xs = np.linspace(xc - R + 1e-9, xc + R - 1e-9, 2000)
    yb = yc - np.sqrt(np.maximum(R ** 2 - (xs - xc) ** 2, 0))
    if yc < 30.0:
        return None
    idx = np.where(yb < ground(xs))[0]
    if len(idx) < 10:
        return None
    blk = max(np.split(idx, np.where(np.diff(idx) > 1)[0] + 1), key=len)
    xa, xb = xs[blk[0]], xs[blk[-1]]
    if xa < 0 or xb > 70 or np.any(yb[blk] < 0.0) or blk[0] == 0 or blk[-1] == len(xs) - 1:
        return None
    edges = np.linspace(xa, xb, n + 1)
    xm, b = 0.5 * (edges[1:] + edges[:-1]), np.diff(edges)
    ybm = yc - np.sqrt(np.maximum(R ** 2 - (xm - xc) ** 2, 0))
    h = np.maximum(ground(xm) - ybm, 0)
    W = gam * b * h
    alpha = np.arcsin(np.clip((xc - xm) / R, -1, 1))
    u = pore(xm, ybm) if pore else np.zeros_like(xm)
    drive = np.sum(W * np.sin(alpha))
    if drive <= 1e-6:
        return None
    F = 1.5
    for _ in range(200):
        ma = np.cos(alpha) * (1 + np.tan(alpha) * tanphi / F)
        if np.any(ma <= 0.2):
            return None
        Fn = np.sum((c * b + np.maximum(W - u * b, 0) * tanphi) / ma) / drive
        if abs(Fn - F) < 1e-9:
            break
        F = Fn
    return F


def search(c, tanphi, pore, gam, coarse=1):
    """grid of circles, then a pattern search started from the 8 best distinct grid circles"""
    cand = []
    for xc in np.linspace(20, 60, 41 // coarse):
        for yc in np.linspace(31, 91, 61 // coarse):
            for R in np.linspace(2, 60, 59 // coarse):
                F = fs_circle(xc, yc, R, c, tanphi, pore, gam, 40)
                if F is not None:
                    cand.append((F, xc, yc, R))
    cand.sort()
    starts = []
    for F, xc, yc, R in cand:
        if all(abs(xc - s[1]) + abs(yc - s[2]) + abs(R - s[3]) > 6 for s in starts):
            starts.append((F, xc, yc, R))
        if len(starts) == 8:
            break
    best = (1e9, None)
    for F0, xc, yc, R in starts:
        F0 = fs_circle(xc, yc, R, c, tanphi, pore, gam, 150) or F0
        for step in [2.0, 1.0, 0.5, 0.2, 0.1, 0.05]:
            improved = True
            while improved:
                improved = False
                for d in [(step, 0, 0), (-step, 0, 0), (0, step, 0), (0, -step, 0), (0, 0, step), (0, 0, -step)]:
                    F = fs_circle(xc + d[0], yc + d[1], R + d[2], c, tanphi, pore, gam, 150)
                    if F is not None and F < F0 - 1e-10:
                        F0, xc, yc, R = F, xc + d[0], yc + d[1], R + d[2]
                        improved = True
        if F0 < best[0]:
            best = (F0, (xc, yc, R))
    return best


arg = sys.argv[1]
pore, gam = None, gsat
if arg == 'crest':
    gam = gsat - gw  # submerged, hydrostatic: buoyant weight, u = 0
elif arg == 'drained':
    pore = lambda x, y: gw * np.maximum(30.0 - y, 0.0)  # water table at the toe (y = 30), hydrostatic
elif arg != 'dry':
    pore = Field(arg)
F, circ = search(c0, np.tan(phi0), pore, gam)
lo, hi = 0.2, 6.0  # gravity increase: c / lambda (weight and pore pressure scaled together)
while hi - lo > 1e-3:
    lam = 0.5 * (lo + hi)
    Fl, _ = search(c0 / lam, np.tan(phi0), pore, gam, coarse=2)
    if Fl < 1.0:
        hi = lam
    else:
        lo = lam
print("%s: Bishop FS(SRM) = %.4f circle (xc, yc, R) = (%.2f, %.2f, %.2f); FS(GI) = %.3f" % (arg, F, *circ, lo))
