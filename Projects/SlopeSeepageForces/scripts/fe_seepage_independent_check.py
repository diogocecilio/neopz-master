#!/usr/bin/env python3
"""Independent cross-check of fe_seepage.py (adversarial review): a second P2 FE solver of Eqs. (20)-(22) of
Ceron et al. (IJNAMG 2025) written from scratch, compared with fe_seepage.FESeepage.

Independent ingredients (nothing is imported from fe_seepage for the solution):
  - mesh: graded point cloud made of hexagonal lattices (one per size level s = hmax / 2^k, kept where
    s <= h < 2 s), boundary nodes by equidistribution of int ds / h, scipy Delaunay, triangles outside the soil
    removed, conformity of every boundary segment checked (midpoint insertion otherwise), exact area checked;
    size h = min(hmax, max(hmin, c r), hface + gface d_face), r = distance to O, W, T (pure geometric grading);
  - P2 assembly on the reference triangle with the 3-point Gauss rule (exact for the P2 stiffness), edge dofs
    from sorted vertex pairs;
  - Dirichlet dofs classified GEOMETRICALLY (crest y = 0, x <= 0; face segment O-T; toe ground y = H, x >= x_T;
    far sides x = x_l, y = y_b), data of Eq. (21) coded separately;
  - point location with scipy.spatial.Delaunay.find_simplex.

Checks (log in results/fe_seepage/independent_check.txt):
  A) J/(k_h H^2 gamma_w^2) and the force -grad u at grid points within 2 H of the face, 5 cases
     (beta, alpha, h_w/H) = (30, 1, 1) zero_lb and impermeable, (60, 5, 0.4), (90, 5, 1), (90, 1, 0.4);
  B) Fig. 5 dashed curve spot values (vector data) at 7 (alpha, beta) points, zero_lb and impermeable;
  C) field API of fe_seepage: timing of 1e5 points, points on the ground surface, shapes, outside points.

usage: python3 Projects/SlopeSeepageForces/scripts/fe_seepage_independent_check.py   (about 1-2 min, 1 CPU)
Coordinates: x right, y DOWN, origin at the crest edge O.  H = 1.
"""
import os
import sys
import time

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.spatial import Delaunay, cKDTree


class Geo:
    def __init__(self, beta, hw, H=1.0, L=50.0, R=10.0, D=30.0):
        self.beta, self.H, self.hw = np.radians(beta), H, hw
        self.xT = H * np.cos(self.beta) / np.sin(self.beta) if beta < 90 else 0.0
        self.xl, self.xr, self.yb = -L * H, self.xT + R * H, H + D * H
        self.O, self.T = np.array([0.0, 0.0]), np.array([self.xT, H])
        self.W = np.array([hw / H * self.xT, hw])
        self.sing = [self.O, self.T] + ([self.W] if 0 < hw < H else [])

    def ground(self, x):
        x = np.asarray(x, float)
        if self.xT == 0.0:
            return np.where(x <= 0, 0.0, self.H)
        return np.clip(x / self.xT, 0.0, 1.0) * self.H

    def corners(self):
        """polygon vertices (closed loop): crest-left, O, [W], T, toe-right, bottom-right, bottom-left"""
        P = [np.array([self.xl, 0.0]), self.O] + ([self.W] if 0 < self.hw < self.H else []) + \
            [self.T, np.array([self.xr, self.H]), np.array([self.xr, self.yb]), np.array([self.xl, self.yb])]
        return np.array(P)


def size_fun(geo, c, hmin, hmax, hface, gface):
    """h = min(hmax, max(hmin, c r_sing), hface + gface d_face)"""
    def h(x, y):
        x, y = np.asarray(x, float), np.asarray(y, float)
        r = np.full(x.shape, np.inf)
        for p in geo.sing:
            r = np.minimum(r, np.hypot(x - p[0], y - p[1]))
        d = geo.T - geo.O
        t = np.clip(((x - geo.O[0]) * d[0] + (y - geo.O[1]) * d[1]) / (d @ d), 0, 1)
        df = np.hypot(x - geo.O[0] - t * d[0], y - geo.O[1] - t * d[1])
        return np.minimum(np.minimum(hmax, np.maximum(hmin, c * r)), hface + gface * df)
    h.c, h.hface, h.gface = c, hface, gface
    return h


def march(P, Q, h):
    """points on segment P->Q with local spacing ~ h (equidistribution of int ds / h, integrated numerically
    on a fine grid refined geometrically towards both ends)"""
    L = np.hypot(*(Q - P))
    s = np.unique(np.concatenate([np.linspace(0, L, 20001), L * np.geomspace(1e-9, 1, 3000),
                                  L - L * np.geomspace(1e-9, 1, 3000)]))
    s = np.clip(s, 0, L)
    s = np.unique(s)
    xy = P + np.outer(s / L, Q - P)
    f = 1.0 / h(xy[:, 0], xy[:, 1])
    N = np.concatenate([[0], np.cumsum(0.5 * (f[1:] + f[:-1]) * np.diff(s))])
    n = max(1, int(round(N[-1])))
    sk = np.interp(np.linspace(0, N[-1], n + 1), N, s)
    return P + np.outer(sk / L, Q - P)


def hex_points(x0, x1, y0, y1, s, X0, Y0):
    """points of the global hexagonal lattice (spacing s, anchored near (X0, Y0)) inside the window"""
    dy = s * np.sqrt(3) / 2
    oy = Y0 + 0.071 * s
    ox = X0 + 0.137 * s
    j0, j1 = int(np.floor((y0 - oy) / dy)), int(np.ceil((y1 - oy) / dy))
    i0, i1 = int(np.floor((x0 - ox) / s)) - 1, int(np.ceil((x1 - ox) / s)) + 1
    if (i1 - i0) * (j1 - j0) > 6e6:
        raise RuntimeError("too many lattice points")
    j = np.arange(j0, j1 + 1)
    i = np.arange(i0, i1 + 1)
    X = ox + s * (i[None, :] + 0.5 * (j[:, None] % 2))
    Y = oy + dy * j[:, None] + 0 * i[None, :]
    P = np.stack([X.ravel(), Y.ravel()], 1)
    return P[(P[:, 0] >= x0) & (P[:, 0] <= x1) & (P[:, 1] >= y0) & (P[:, 1] <= y1)]


def make_mesh(geo, h, hmin, hmax):
    V = geo.corners()
    B = []
    for k in range(len(V)):
        B.append(march(V[k], V[(k + 1) % len(V)], h)[:-1])
    B = np.vstack(B)
    # interior: hex lattice for each level s = hmax / 2^k, kept where s <= h < 2 s
    pts = []
    s = hmax
    X0, Y0 = geo.xl, 0.0
    while s > 0.5 * hmin:
        if s > 0.6:
            P = hex_points(geo.xl, geo.xr, 0.0, geo.yb, s, X0, Y0)
        else:
            wins = []
            if 2 * s > h.hface:
                wf = (2 * s - h.hface) / h.gface + 2 * s
                wins.append((min(0, geo.xT) - wf, max(0, geo.xT) + wf, -wf, geo.H + wf))
            ws = 2 * s / h.c + 2 * s
            for p in geo.sing:
                wins.append((p[0] - ws, p[0] + ws, p[1] - ws, p[1] + ws))
            P = np.vstack([hex_points(max(a, geo.xl), min(b, geo.xr), max(c, 0.0), min(d, geo.yb), s, X0, Y0)
                           for a, b, c, d in wins])
            P = np.unique(np.round(P / (1e-6 * s)).astype(np.int64), axis=0) * (1e-6 * s)
        hv = h(P[:, 0], P[:, 1])
        if s <= hmin:
            keep = hv < 2 * s
        elif s >= hmax:
            keep = hv >= s
        else:
            keep = (hv >= s) & (hv < 2 * s)
        pts.append(P[keep])
        s *= 0.5
    P = np.vstack(pts)
    inside = (P[:, 0] > geo.xl) & (P[:, 0] < geo.xr) & (P[:, 1] < geo.yb) & (P[:, 1] > geo.ground(P[:, 0]))
    P = P[inside]
    # distance to the boundary polygon
    dmin = np.full(len(P), np.inf)
    for k in range(len(V)):
        a, b = V[k], V[(k + 1) % len(V)]
        e = b - a
        t = np.clip(((P[:, 0] - a[0]) * e[0] + (P[:, 1] - a[1]) * e[1]) / (e @ e), 0, 1)
        dmin = np.minimum(dmin, np.hypot(P[:, 0] - a[0] - t * e[0], P[:, 1] - a[1] - t * e[1]))
    P = P[dmin > 0.6 * h(P[:, 0], P[:, 1])]
    for it in range(12):
        nb = len(B)
        # remove interior points inside the diametral circles of boundary segments
        mid = 0.5 * (B + np.roll(B, -1, 0))
        rad = 0.5 * np.hypot(*(np.roll(B, -1, 0) - B).T) * 1.02
        tree = cKDTree(P)
        bad = set()
        for lst in tree.query_ball_point(mid, rad):
            bad.update(lst)
        Pin = np.delete(P, sorted(bad), 0)
        X = np.vstack([B, Pin])
        tri = Delaunay(X)
        S = tri.simplices
        c = X[S].mean(1)
        ok = (c[:, 1] > geo.ground(c[:, 0])) & (c[:, 0] > geo.xl) & (c[:, 0] < geo.xr) & (c[:, 1] < geo.yb)
        T = S[ok]
        edges = set(map(tuple, np.sort(np.vstack([T[:, [0, 1]], T[:, [1, 2]], T[:, [0, 2]]]), 1).tolist()))
        miss = [k for k in range(nb) if tuple(sorted((k, (k + 1) % nb))) not in edges]
        if not miss:
            break
        newB = []
        for k in range(nb):
            newB.append(B[k])
            if k in miss:
                newB.append(0.5 * (B[k] + B[(k + 1) % nb]))
        B = np.array(newB)
    else:
        raise RuntimeError("boundary not conforming")
    # orientation
    x = X[T]
    det = (x[:, 1, 0] - x[:, 0, 0]) * (x[:, 2, 1] - x[:, 0, 1]) - (x[:, 1, 1] - x[:, 0, 1]) * (x[:, 2, 0] - x[:, 0, 0])
    area = 0.5 * np.abs(det).sum()
    exact = (geo.xr - geo.xl) * geo.yb - (geo.xr - geo.xT) * geo.H - 0.5 * geo.xT * geo.H
    assert abs(area / exact - 1) < 1e-9, (area, exact)
    simplex_map = -np.ones(len(S), int)
    simplex_map[np.nonzero(ok)[0]] = np.arange(len(T))
    return X, T, tri, simplex_map


# P2 on the reference triangle: lam0 = 1 - xi - eta, lam1 = xi, lam2 = eta; local dofs v0 v1 v2 m01 m12 m20
def ref_grad(xi, eta):
    l0, l1, l2 = 1 - xi - eta, xi, eta
    dl = np.array([[-1.0, -1.0], [1.0, 0.0], [0.0, 1.0]])
    l = [l0, l1, l2]
    g = [(4 * l[i] - 1) * dl[i] for i in range(3)]
    for i, j in ((0, 1), (1, 2), (2, 0)):
        g.append(4 * (l[i] * dl[j] + l[j] * dl[i]))
    return np.array(g)          # (6, 2) d/dxi, d/deta


def ref_val(xi, eta):
    l = [1 - xi - eta, xi, eta]
    return np.array([l[0] * (2 * l[0] - 1), l[1] * (2 * l[1] - 1), l[2] * (2 * l[2] - 1),
                     4 * l[0] * l[1], 4 * l[1] * l[2], 4 * l[2] * l[0]])


GP = [((1 / 6, 1 / 6), 1 / 6), ((2 / 3, 1 / 6), 1 / 6), ((1 / 6, 2 / 3), 1 / 6)]


class P2Solver:
    def __init__(self, geo, mesh, kx, ky, gw=9.81, bc="zero_lb"):
        X, T, self.tri, self.smap = mesh
        self.geo, self.X, self.T, self.kx, self.ky, self.gw = geo, X, T, kx, ky, gw
        nv, ne = len(X), len(T)
        # edge dofs
        E = np.sort(np.vstack([T[:, [0, 1]], T[:, [1, 2]], T[:, [2, 0]]]), 1)
        key = E[:, 0].astype(np.int64) * nv + E[:, 1]
        uk, inv = np.unique(key, return_inverse=True)
        self.dofs = np.hstack([T, nv + inv.reshape(3, ne).T])
        ea, eb = uk // nv, uk % nv
        self.Xd = np.vstack([X, 0.5 * (X[ea] + X[eb])])
        nd = len(self.Xd)
        # Jacobians: x = x0 + J (xi, eta)
        x0, x1, x2 = X[T[:, 0]], X[T[:, 1]], X[T[:, 2]]
        Jm = np.stack([x1 - x0, x2 - x0], 2)                 # (ne, 2, 2): columns dx/dxi, dx/deta
        detJ = Jm[:, 0, 0] * Jm[:, 1, 1] - Jm[:, 0, 1] * Jm[:, 1, 0]
        Jinv = np.empty_like(Jm)
        Jinv[:, 0, 0], Jinv[:, 1, 1] = Jm[:, 1, 1] / detJ, Jm[:, 0, 0] / detJ
        Jinv[:, 0, 1], Jinv[:, 1, 0] = -Jm[:, 0, 1] / detJ, -Jm[:, 1, 0] / detJ
        self.Jinv, self.detJ, self.x0 = Jinv, detJ, x0
        K = np.diag([kx, ky])
        Ke = np.zeros((ne, 6, 6))
        for (xi, eta), w in GP:
            gr = ref_grad(xi, eta)                             # (6, 2)
            # physical gradient: grad N = Jinv^T grad_ref N
            gp = np.einsum("eba,ib->eia", Jinv, gr)            # (ne, 6, 2)
            Ke += w * np.abs(detJ)[:, None, None] * np.einsum("eia,ab,ejb->eij", gp, K, gp)
        r = np.repeat(self.dofs, 6, 1).ravel()
        c = np.tile(self.dofs, (1, 6)).ravel()
        A = sp.coo_matrix((Ke.ravel(), (r, c)), shape=(nd, nd)).tocsr()
        self.A = A
        # geometric classification of the boundary dofs
        x, y = self.Xd[:, 0], self.Xd[:, 1]
        tol = 1e-9
        g = geo
        crest = (np.abs(y) < tol) & (x <= tol)
        d = g.T - g.O
        t = np.clip((x * d[0] + y * d[1]) / (d @ d), 0, 1)
        face = np.hypot(x - t * d[0], y - t * d[1]) < tol
        toe = (np.abs(y - g.H) < tol) & (x >= g.xT - tol)
        surf = crest | face | toe
        val = np.zeros(nd)
        val[face] = -gw * np.minimum(y[face], g.hw)
        val[toe] = -gw * g.hw
        val[crest] = 0.0
        fixed = surf.copy()
        if bc == "zero_lb":
            far = ((np.abs(x - g.xl) < tol) | (np.abs(y - g.yb) < tol)) & ~surf
            fixed |= far
        elif bc != "impermeable":
            raise ValueError(bc)
        u = np.where(fixed, val, 0.0)
        free = ~fixed
        rhs = -(A @ u)
        Aff = A[free][:, free].tocsc()
        u[free] = spla.splu(Aff, permc_spec="COLAMD").solve(rhs[free])
        self.u = u
        self.J = 0.5 * u @ (A @ u)
        self.nfix = fixed.sum()

    def grad(self, px, py):
        px, py = np.asarray(px, float).ravel(), np.asarray(py, float).ravel()
        s = self.tri.find_simplex(np.stack([px, py], 1))
        e = np.where(s >= 0, self.smap[np.maximum(s, 0)], -1)
        gx, gy = np.zeros(len(px)), np.zeros(len(px))
        ok = e >= 0
        eo = e[ok]
        dx = np.stack([px[ok], py[ok]], 1) - self.x0[eo]
        ref = np.einsum("eab,eb->ea", self.Jinv[eo], dx)        # (xi, eta)
        xi, eta = ref[:, 0], ref[:, 1]
        l = [1 - xi - eta, xi, eta]
        dl = np.array([[-1.0, -1.0], [1.0, 0.0], [0.0, 1.0]])
        gr = np.zeros((len(eo), 6, 2))
        for i in range(3):
            gr[:, i] = (4 * l[i] - 1)[:, None] * dl[i]
        for k, (i, j) in enumerate(((0, 1), (1, 2), (2, 0))):
            gr[:, 3 + k] = 4 * (l[i][:, None] * dl[j] + l[j][:, None] * dl[i])
        gp = np.einsum("eba,eib->eia", self.Jinv[eo], gr)
        ue = self.u[self.dofs[eo]]
        g = np.einsum("ei,eia->ea", ue, gp)
        gx[ok], gy[ok] = g[:, 0], g[:, 1]
        return gx, gy

    def force(self, px, py):
        gx, gy = self.grad(px, py)
        return -gx, -gy


# ======================================================================================================
# comparison with fe_seepage.py
# ======================================================================================================
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.dirname(SCRIPT_DIR)
GW = 9.81


def indep_solve(beta, alpha, hw, bc, lev):
    """independent solution at level lev (sizes halved per level, hmin quartered)"""
    geo = Geo(beta, hw)
    f = 0.5 ** lev
    h = size_fun(geo, c=0.5 * f, hmin=0.016 * f * f, hmax=6.0 * f, hface=0.1 * f, gface=0.4 * f)
    mesh = make_mesh(geo, h, 0.016 * f * f, 6.0 * f)
    return geo, P2Solver(geo, mesh, 1.0, 1.0 / alpha, GW, bc)


def samples(geo, rmax=2.0, step=0.05):
    """grid points of the soil within rmax of the face segment O-T"""
    xs = np.arange(-rmax - 0.5, geo.xT + rmax + 0.5, step) + 0.0123
    ys = np.arange(0, geo.H + rmax, step) + 0.0171
    X, Y = np.meshgrid(xs, ys)
    X, Y = X.ravel(), Y.ravel()
    d = geo.T - geo.O
    t = np.clip((X * d[0] + Y * d[1]) / (d @ d), 0, 1)
    df = np.hypot(X - t * d[0], Y - t * d[1])
    keep = (df <= rmax) & (Y > geo.ground(X) + 1e-3)
    return X[keep], Y[keep]


def main():
    sys.path.insert(0, SCRIPT_DIR)
    import fe_seepage as fs
    lines = []

    def log(s=""):
        print(s, flush=True)
        lines.append(s)

    t0 = time.time()
    log("=" * 110)
    log("A) independent P2 solver (levels 2, 3) vs fe_seepage.FESeepage (P2 ref 1 = default, ref 2), box 50/10/30 H")
    log("   J = J/(k_h H^2 gw^2); force compared on a 0.05 H grid of soil points within 2 H of the face;")
    log("   rel = RMS|f - f_indep| / RMS|f_indep|, also for the points farther than 0.05 H from O, W, T")
    log("=" * 110)
    cases = [(30, 1, 1.0, "zero_lb"), (30, 1, 1.0, "impermeable"), (60, 5, 0.4, "zero_lb"), (90, 5, 1.0, "zero_lb"),
             (90, 1, 0.4, "zero_lb")]
    for beta, a, hw, bc in cases:
        geo, s2 = indep_solve(beta, a, hw, bc, 2)
        _, s3 = indep_solve(beta, a, hw, bc, 3)
        fe1 = fs.FESeepage(beta, 1.0, hw, a, gamma_w=GW, bc=bc, order=2, ref=1)
        fe2 = fs.FESeepage(beta, 1.0, hw, a, gamma_w=GW, bc=bc, order=2, ref=2)
        J2, J3 = s2.J / GW ** 2, s3.J / GW ** 2
        log(f"  beta={beta:2d} alpha={a:2d} hw/H={hw:3.1f} {bc:11s}: J indep lev2 {J2:.7f} lev3 {J3:.7f} (ndof {len(s3.u)}); "
            f"fe_seepage ref1 {fe1.J_normalized():.7f} ref2 {fe2.J_normalized():.7f}; ref2 / lev3 - 1 = "
            f"{fe2.J_normalized() / J3 - 1:+.1e}")
        X, Y = samples(geo)
        f0 = np.array(s3.force(X, Y))
        r = np.min([np.hypot(X - p[0], Y - p[1]) for p in geo.sing], 0)
        far = r > 0.05
        nb, nbf = np.sqrt(np.mean(np.sum(f0 ** 2, 0))), np.sqrt(np.mean(np.sum(f0[:, far] ** 2, 0)))
        for lab, fe in (("ref1", fe1), ("ref2", fe2)):
            dn = np.hypot(*(np.array(fe.force(X, Y)) - f0))
            log(f"      force fe_seepage {lab} vs indep lev3 ({len(X)} pts): rel {100 * np.sqrt(np.mean(dn ** 2)) / nb:.3f}% "
                f"(beyond 0.05 H of O/W/T {100 * np.sqrt(np.mean(dn[far] ** 2)) / nbf:.3f}%), max |df| {dn.max():.4f} "
                f"kN/m^3 (max |f| {np.hypot(*f0).max():.2f})")
        dn = np.hypot(*(np.array(s2.force(X, Y)) - f0))
        log(f"      own discretisation: indep lev2 vs lev3 rel {100 * np.sqrt(np.mean(dn ** 2)) / nb:.3f}%")
    log("=" * 110)
    log("B) Fig. 5 dashed curve (data/fig5_vector_fill_polygons.csv), h_w = H: independent solver (level 2)")
    log("=" * 110)
    paper = np.loadtxt(os.path.join(PROJECT_DIR, "data", "fig5_vector_fill_polygons.csv"), delimiter=",", comments="#")
    for a, b in ((1, 15), (1, 45), (1, 90), (2, 60), (4, 30), (10, 15), (10, 90)):
        p = paper[(paper[:, 0] == a) & (paper[:, 1] == b)][0, 3]
        out = []
        for bc in ("zero_lb", "impermeable"):
            out.append(indep_solve(b, a, 1.0, bc, 2)[1].J / GW ** 2)
        log(f"  alpha={a:2d} beta={b:2d}: paper {p:.5f}; indep zero_lb {out[0]:.6f} ({100 * (out[0] / p - 1):+.3f}%), "
            f"impermeable {out[1]:.6f} ({100 * (out[1] / p - 1):+.2f}%)")
    log("=" * 110)
    log("C) fe_seepage field API (default P2 ref 1, H = 5, beta = 60, h_w = H, alpha = 5)")
    log("=" * 110)
    rng = np.random.default_rng(0)
    fe = fs.FESeepage(60.0, 5.0, 5.0, 5.0, gamma_w=GW)
    x, y = rng.uniform(-15, 15, 100000), rng.uniform(-2, 15, 100000)
    t = time.time()
    fx, fy = fe.force(x, y)
    t1 = time.time() - t
    t = time.time()
    fe.force(x, y)
    log(f"  1e5 random points: first call {t1:.2f} s (incl. locator build), second {time.time() - t:.2f} s; "
        f"zero outside the box/soil: {bool(np.all((fx[~fe.dom.in_box(x, y)] == 0) & (fy[~fe.dom.in_box(x, y)] == 0)))}")
    d = fe.dom
    s = rng.uniform(0, 1, 2000)
    for name, (xs, ys) in (("crest", (rng.uniform(-10, 0, 2000), np.zeros(2000))), ("face", (s * d.xT, s * d.H)),
                           ("toe ground", (d.xT + rng.uniform(0, 10, 2000), np.full(2000, d.H)))):
        fx, fy = fe.force(xs, ys)
        uu = fe.u(xs, ys)
        log(f"  2000 points on the {name:10s}: f = 0 at {int(np.sum((fx == 0) & (fy == 0)))}, u = nan at "
            f"{int(np.isnan(uu).sum())}, max |u - Eq. 21| = {np.nanmax(np.abs(uu + GW * np.clip(ys, 0, fe.hw))):.1e}")
    shapes = [np.shape(fe.force(*arg)[0]) for arg in ((0.0, 3.0), ([0.0, -1.0], [3.0, 3.0]),
                                                      (np.zeros((3, 4)), np.full((3, 4), 3.0)),
                                                      (np.zeros((2, 1)), np.full((1, 5), 3.0)))]
    log(f"  output shapes for inputs (), (2,), (3, 4), (2, 1) x (1, 5): {shapes}; non-finite input -> "
        f"{fe.force(np.array([np.nan, 1.0]), np.array([1.0, np.nan]))}")
    log(f"total time {time.time() - t0:.1f} s")
    out = os.path.join(PROJECT_DIR, "results", "fe_seepage", "independent_check.txt")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, "w") as fh:
        fh.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
