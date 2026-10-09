#!/usr/bin/env python3
"""Finite element solution of the uncoupled hydraulic problem of Ceron et al. (IJNAMG 2025), Eqs. (20)-(22).

Problem (paper coordinates of the shared spec: origin O at the crest edge, x to the right, y DOWNWARD):

    div(-K grad u) = 0 in Omega,      K = k_v e_y (x) e_y + k_h (1 - e_y (x) e_y) = diag(k_h, k_h / alpha)
    u = 0              on the crest (y = 0, x <= 0)                         dOmega_u1
    u = -gamma_w y     on the face for 0 <= y <= h_w                         dOmega_u2
    u = -gamma_w h_w   on the face for y >= h_w and on the toe ground y = H  dOmega_u2 / dOmega_u3

u is the excess pore pressure (kPa).  The FE solution minimises J(u') = 1/2 int grad u' . K . grad u' dOmega
(Eq. 22 with dOmega_v = empty on the ground surface) over continuous piecewise P1 or P2 fields that take
the Dirichlet data; the data are piecewise linear along the ground surface with kinks at O, W (water line
on the face) and T (toe), which are mesh nodes, so the Lagrange interpolant of the data is exact and the
discrete minimum is an UPPER bound of the exact J (it converges from above under refinement).

The infinite half-space of the paper is truncated to a box (FE domain of Fig. 4a): `left` to the left of O,
`right` to the right of the toe T and `depth` below T (in units of H, as SlopeGeometry.h of the C++ code).
Each far side (left, bottom, right) is either impermeable (zero flux, natural BC) or carries Dirichlet data
u = 0 ("zero") or u = -gamma_w h_w ("toe"); named combinations in BC_PRESETS.

Paper configuration (DEFAULT, decided by the checks 5 and 6 below): box left 50 H / right 10 H / depth 30 H below
the toe (measured on the native 1000 x 1728 px image of Fig. 4 embedded in the PDF: box 938 x 477 px = 61 H x 31 H
for beta = 45 deg) and preset "zero_lb": u = 0 on the left side and on the bottom, zero flux on the right side.
It reproduces the dashed curves J(u'_FE) of Fig. 5 (h_w = H) within 0.25 % for alpha = 1, 2, 4, 10 and
beta = 15..90 deg (always slightly BELOW the paper, as expected from a converged lower value of an upper-bound
functional against the paper's coarse 6-node mesh), and the iso-lines of Fig. 4b (they end on the right side; the
"0.00" at the bottom-right corner is u there).  The fully impermeable box ("impermeable") gives J 3-11 % lower and
iso-lines that end on the bottom: it does NOT reproduce the paper.

Mesh (default method "quadtree"): conforming Delaunay triangulation (scipy.spatial.Delaunay) of
  - boundary nodes distributed on the polygon edges by equidistribution of int ds / h (O, W, T are nodes), and
  - the vertices of an axis-aligned quadtree whose leaves have side <= h(centre) (graded, deterministic),
    without the vertices that are closer than 0.5 h to the boundary or inside the diametral circle of a boundary
    segment (Gabriel condition: every boundary segment is a Delaunay edge; checked, with segment splitting as a
    fallback), triangles outside the soil removed.  Min angle about 18-25 deg, max angle < 135 deg for any beta.
The target size h is graded towards O, W and T (T and W are singular points of grad u) and the face:
      h = min(hmax, h0 + grade * dist(O, W, T), hs + grade * dist(face));
refinement level `ref` divides hs, grade and hmax by 2**ref and h0 by 4**ref (so that the corner layer does
not limit the convergence rate of the singular solution at the re-entrant toe T).
Alternative method "blocks": three mapped (transfinite) blocks of quadrilaterals split into triangles, below the
crest, the face and the toe ground (the previous generator; strongly stretched elements away from the slope, min
angle < 1 deg, max angle up to 163 deg for beta = 90; kept for the independent-mesh cross-check of J).

Field interface of the shared spec:  field.force(x, y) = -grad u_FE (kN/m^3), field.u(x, y),
field.velocity(x, y) = -K grad u_FE; point location with matplotlib.tri.TrapezoidMapTriFinder, plus a retry
for points within 1e-9 H of the soil (points ON the crest, face and toe ground are part of the soil);
0 force (nan u) outside the FE domain and for non-finite input.  The mesh is generated for H = 1 and scaled,
so u / (gamma_w H) and f / gamma_w depend on (x / H, y / H) only (exact similarity in H).

Run as a script to execute the verification checks and the comparison with Figs. 4b and 5 of the paper:
    python3 Projects/SlopeSeepageForces/scripts/fe_seepage.py [--quick] [--no-plot] [--only fig5,fig4,...]
(output in results/fe_seepage/checks_output.txt, figures in results/fe_seepage and data/fe_seepage_fig4_check.png)
"""
from __future__ import annotations

import argparse
import importlib.util
import os
import sys
import time

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.dirname(SCRIPT_DIR)
DATA_DIR = os.path.join(PROJECT_DIR, "data")
RESULTS_DIR = os.path.join(PROJECT_DIR, "results", "fe_seepage")
REPO_DIR = os.path.dirname(os.path.dirname(PROJECT_DIR))

# boundary ids (same numbering as SlopeGeometry.h / TriGMesh)
BASE, RIGHT, TOE, CREST, LEFT, FACE = -1, -2, -3, -4, -5, -6
SURFACE_IDS = (CREST, FACE, TOE)
FAR_IDS = {"left": LEFT, "bottom": BASE, "right": RIGHT}

# named far-side boundary conditions ("paper" = "zero_lb", the default: it reproduces Figs. 4b and 5)
BC_PRESETS = {
    "impermeable": {"left": "neumann", "bottom": "neumann", "right": "neumann"},
    "zero_lb": {"left": "zero", "bottom": "zero", "right": "neumann"},      # u = 0 left + bottom, right no-flow
    "paper": {"left": "zero", "bottom": "zero", "right": "neumann"},        # alias of zero_lb
    "zero_b": {"left": "neumann", "bottom": "zero", "right": "neumann"},    # u = 0 on the bottom only
    "zero_l": {"left": "zero", "bottom": "neumann", "right": "neumann"},    # u = 0 on the left side only
    "zero_lbr": {"left": "zero", "bottom": "zero", "right": "zero"},        # data jump at the right top corner
    "toe_r": {"left": "zero", "bottom": "neumann", "right": "toe"},         # u = -gw hw on the right side
}


def _bc_dict(bc):
    if isinstance(bc, str):
        return dict(BC_PRESETS[bc])
    d = dict(BC_PRESETS["impermeable"])
    d.update(bc)
    return d


# ======================================================================================================
# geometry
# ======================================================================================================
class SlopeDomain:
    """Slope of height H and angle beta (paper coordinates, y down) truncated to a box.

    left, right, depth : extents in units of H (left of O, right of the toe T, below T).
    """

    def __init__(self, beta_deg, H=1.0, hw=None, left=50.0, right=10.0, depth=30.0):
        self.beta_deg = float(beta_deg)
        if not (0.0 < self.beta_deg <= 90.0):
            raise ValueError("beta must be in (0, 90] deg")
        self.beta = np.radians(self.beta_deg)
        self.H = float(H)
        self.hw = self.H if hw is None else float(hw)
        if not (0.0 <= self.hw <= self.H * (1 + 1e-12)):
            raise ValueError("hw must be in [0, H]")
        self.hw = min(self.hw, self.H)
        self.sb, self.cb = np.sin(self.beta), np.cos(self.beta)
        if self.beta_deg >= 90.0:
            self.cb = 0.0
        self.cot = self.cb / self.sb
        self.xT = self.H * self.cot
        self.left, self.right, self.depth = float(left), float(right), float(depth)
        self.xl = -self.left * self.H
        self.xr = self.xT + self.right * self.H
        self.yb = self.H + self.depth * self.H
        self.O = np.array([0.0, 0.0])
        self.T = np.array([self.xT, self.H])
        self.W = np.array([self.hw * self.cot, self.hw])
        self.has_W = 1e-9 * self.H < self.hw < self.H * (1 - 1e-9)

    def ground(self, x):
        """depth of the ground surface below the crest level at abscissa x"""
        x = np.asarray(x, float)
        if self.xT <= 0.0:
            return np.where(x <= 0.0, 0.0, self.H)
        return np.where(x <= 0.0, 0.0, np.where(x >= self.xT, self.H, x * self.sb / self.cb))

    def in_box(self, x, y, tol=0.0):
        x, y = np.asarray(x, float), np.asarray(y, float)
        return (x >= self.xl - tol) & (x <= self.xr + tol) & (y <= self.yb + tol) & (y >= self.ground(x) - tol)

    def surface_u(self, x, y, gamma_w):
        """Dirichlet data of Eq. (21) at points of the ground surface (crest, face, toe ground)."""
        y = np.asarray(y, float)
        return -gamma_w * np.clip(y, 0.0, self.hw)

    def polygon(self):
        """boundary polygon (positive orientation in (x, y)) and the boundary id of each edge V[k] -> V[k+1]"""
        H = self.H
        V = [np.array([self.xl, 0.0]), self.O] + ([self.W] if self.has_W else []) + \
            [self.T, np.array([self.xr, H]), np.array([self.xr, self.yb]), np.array([self.xl, self.yb])]
        ids = [CREST] + [FACE] * (2 if self.has_W else 1) + [TOE, RIGHT, BASE, LEFT]
        return np.array(V), ids

    def area(self):
        return (self.xr - self.xl) * self.yb - (self.xr - self.xT) * self.H - 0.5 * self.xT * self.H

    def dist_boundary(self, x, y):
        """distance of the points (x, y) to the boundary polygon"""
        V, _ = self.polygon()
        d = np.full(np.shape(x), np.inf)
        for k in range(len(V)):
            P, Q = V[k], V[(k + 1) % len(V)]
            e = Q - P
            t = np.clip(((x - P[0]) * e[0] + (y - P[1]) * e[1]) / (e @ e), 0.0, 1.0)
            d = np.minimum(d, np.hypot(x - P[0] - t * e[0], y - P[1] - t * e[1]))
        return d


# box of the paper's FE model (Fig. 4a/b), in units of H: left of O, right of T, below T
PAPER_BOX = dict(left=50.0, right=10.0, depth=30.0)


# ======================================================================================================
# mesh
# ======================================================================================================
class MeshSize:
    """Target element size (lengths in units of H): h0 at O, W, T; hs along the face; growth `grade` per unit
    distance; at most hmax.  Level `ref` divides hs, grade, hmax by 2**ref and h0 by 4**ref."""

    def __init__(self, h0=0.02, hs=0.1, grade=0.25, hmax=4.0, ref=0, uniform=None):
        self.h0, self.hs, self.grade, self.hmax, self.ref = h0, hs, grade, hmax, int(ref)
        self.uniform = uniform      # if set: constant size uniform * H / 2**ref (no grading)

    def __call__(self, dom, x, y):
        H = dom.H
        f = 2.0 ** (-self.ref)
        if self.uniform is not None:
            return np.full(np.shape(x), self.uniform * H * f)
        h = np.full(np.shape(x), self.hmax * H * f)
        pts = [dom.O, dom.T] + ([dom.W] if dom.has_W else [])
        h0 = self.h0 * H * 4.0 ** (-self.ref)
        g = self.grade * f
        for p in pts:
            h = np.minimum(h, h0 + g * np.hypot(x - p[0], y - p[1]))
        # distance to the face segment O-T
        d = dom.T - dom.O
        t = np.clip(((x - dom.O[0]) * d[0] + (y - dom.O[1]) * d[1]) / (d @ d), 0.0, 1.0)
        dist = np.hypot(x - dom.O[0] - t * d[0], y - dom.O[1] - t * d[1])
        return np.minimum(h, self.hs * H * f + g * dist)


def _edge_nodes(P, Q, sizefun, n=None, foci=()):
    """Nodes on the segment P -> Q distributed with the size function (equidistribution of int ds / h).
    n = number of segments (None: from the size function).  Returns ((n+1, 2) nodes, natural count)."""
    P, Q = np.asarray(P, float), np.asarray(Q, float)
    ell = float(np.hypot(*(Q - P)))
    e = (Q - P) / ell
    fs = [0.0, ell]
    for F in foci:
        s = float(np.dot(np.asarray(F, float) - P, e))
        if 0.0 < s < ell:
            fs.append(s)
    g = np.logspace(-11, 0, 700) * ell
    s = [np.linspace(0.0, ell, 4001)] + [f + g for f in fs] + [f - g for f in fs]
    s = np.unique(np.clip(np.concatenate(s), 0.0, ell))
    xy = P[None, :] + s[:, None] * e[None, :]
    inv = 1.0 / sizefun(xy[:, 0], xy[:, 1])
    N = np.concatenate([[0.0], np.cumsum(0.5 * (inv[1:] + inv[:-1]) * np.diff(s))])
    nat = max(1, int(np.ceil(N[-1] - 1e-6)))
    if n is None:
        n = nat
    sk = np.interp(np.linspace(0.0, N[-1], n + 1), N, s)
    sk[0], sk[-1] = 0.0, ell
    nodes = P[None, :] + sk[:, None] * e[None, :]
    nodes[0], nodes[-1] = P, Q
    return nodes, nat


def _block_nodes(top, bot, left, right):
    """Structured block: node (i, j) = intersection of the straight lines top[i]-bot[i] and left[j]-right[j].
    Corners: top[0] = left[0], top[-1] = right[0], bot[0] = left[-1], bot[-1] = right[-1]."""
    P1, P2 = top[:, None, :], bot[:, None, :]
    Q1, Q2 = left[None, :, :], right[None, :, :]
    d1, d2, r = P2 - P1, Q2 - Q1, Q1 - P1
    cross = lambda a, b: a[..., 0] * b[..., 1] - a[..., 1] * b[..., 0]
    a = cross(r, d2) / cross(d1, d2)
    X = P1 + a[..., None] * d1
    X[:, 0], X[:, -1], X[0, :], X[-1, :] = top, bot, left, right
    return X                                                  # (ni, nj, 2)


class Mesh:
    """Triangle mesh: X (n, 2) vertices, T (m, 3) triangles (positive orientation in (x, y)),
    edges (k, 2) boundary edges with ids eid (k,)."""

    def __init__(self, X, T, edges, eid):
        self.X, self.T, self.edges, self.eid = X, T, edges, eid

    @property
    def nv(self):
        return len(self.X)

    def quality(self):
        x = self.X[self.T]
        ang = []
        for k in range(3):
            a = x[:, (k + 1) % 3] - x[:, k]
            b = x[:, (k + 2) % 3] - x[:, k]
            c = np.sum(a * b, 1) / (np.linalg.norm(a, axis=1) * np.linalg.norm(b, axis=1))
            ang.append(np.degrees(np.arccos(np.clip(c, -1, 1))))
        ang = np.array(ang)
        L = np.stack([np.linalg.norm(x[:, (k + 1) % 3] - x[:, k], axis=1) for k in range(3)])
        return dict(min_angle=ang.min(), max_angle=ang.max(), hmin=L.min(), hmax=L.max())


def _boundary_loop(dom, sf):
    """boundary nodes (closed loop, no repetition) on the polygon edges and the id of each segment k -> k+1"""
    V, ids = dom.polygon()
    foc = [dom.O, dom.T] + ([dom.W] if dom.has_W else [])
    nodes, eid = [], []
    for k in range(len(V)):
        e, _ = _edge_nodes(V[k], V[(k + 1) % len(V)], sf, foci=foc)
        nodes.append(e[:-1])
        eid.append(np.full(len(e) - 1, ids[k]))
    return np.vstack(nodes), np.concatenate(eid)


def _quadtree_points(dom, sf):
    """vertices of an axis-aligned quadtree over the box whose leaves (cells meeting the soil) have side
    <= h(centre); returns (n, 2) points"""
    width = dom.xr - dom.xl
    corners = np.array([[dom.xl, 0.0], [dom.xr, dom.H], [dom.xl, dom.yb], [dom.xr, dom.yb]])
    hmax = float(np.max(sf(corners[:, 0], corners[:, 1])))
    nx = int(np.ceil(width / hmax - 1e-9))
    s0 = width / nx
    ny = int(np.ceil(dom.yb / s0 - 1e-9))
    i, j = (a.ravel() for a in np.meshgrid(np.arange(nx), np.arange(ny), indexing="ij"))
    leaves = []
    lev = 0
    while True:
        s = s0 / 2 ** lev
        x0, y0 = dom.xl + i * s, j * s
        # the ground depth is non-decreasing in x: the cell meets the soil iff its bottom is below ground(x0)
        meets = (y0 + s > dom.ground(x0) + 1e-12 * dom.H) & (x0 < dom.xr) & (y0 < dom.yb)
        split = meets & (s > sf(x0 + 0.5 * s, y0 + 0.5 * s))
        leaves.append((lev, i[meets & ~split], j[meets & ~split]))
        if not split.any() or lev > 40:
            break
        i2, j2 = i[split], j[split]
        i = np.concatenate([2 * i2, 2 * i2 + 1, 2 * i2, 2 * i2 + 1])
        j = np.concatenate([2 * j2, 2 * j2, 2 * j2 + 1, 2 * j2 + 1])
        lev += 1
    # leaf corners in integer coordinates of the finest level -> unique vertices
    C = []
    for lv, il, jl in leaves:
        f = 2 ** (lev - lv)
        for di in (0, 1):
            for dj in (0, 1):
                C.append(np.stack([(il + di) * f, (jl + dj) * f], 1))
    C = np.unique(np.vstack(C), axis=0)
    sL = s0 / 2 ** lev
    return np.stack([dom.xl + C[:, 0] * sL, C[:, 1] * sL], 1)


def _quadtree_mesh(dom: SlopeDomain, size: MeshSize, excl=0.5, max_fix=10):
    """graded, deterministic conforming Delaunay mesh of the slope box (see module docstring)"""
    from scipy.spatial import Delaunay, cKDTree
    sf = lambda x, y: size(dom, x, y)
    H = dom.H
    B, beid = _boundary_loop(dom, sf)
    P = _quadtree_points(dom, sf)
    x, y = P[:, 0], P[:, 1]
    inside = (x > dom.xl) & (x < dom.xr) & (y < dom.yb) & (y > dom.ground(x))
    P = P[inside]
    P = P[dom.dist_boundary(P[:, 0], P[:, 1]) > excl * sf(P[:, 0], P[:, 1])]
    tree = cKDTree(P)
    for it in range(max_fix + 1):
        nb = len(B)
        bedges = np.stack([np.arange(nb), (np.arange(nb) + 1) % nb], 1)
        # Gabriel condition: no interior point inside (1.05 x) the diametral circle of a boundary segment
        mid = 0.5 * (B[bedges[:, 0]] + B[bedges[:, 1]])
        rad = 0.525 * np.linalg.norm(B[bedges[:, 1]] - B[bedges[:, 0]], axis=1)
        hit = tree.query_ball_point(mid, rad)
        bad = np.unique(np.concatenate([np.asarray(h, dtype=int) for h in hit] + [np.zeros(0, int)]))
        Pin = np.delete(P, bad, axis=0)
        X = np.vstack([B, Pin])
        T = Delaunay(X).simplices.astype(np.int64)
        c = X[T].mean(1)
        T = T[(c[:, 1] > dom.ground(c[:, 0])) & (c[:, 0] > dom.xl) & (c[:, 0] < dom.xr) & (c[:, 1] < dom.yb)]
        # every boundary segment must be an edge of the triangulation (then no triangle crosses the boundary)
        E = np.sort(np.concatenate([T[:, [0, 1]], T[:, [1, 2]], T[:, [2, 0]]]), 1)
        key = E[:, 0] * len(X) + E[:, 1]
        be = np.sort(bedges, 1)
        missing = ~np.isin(be[:, 0] * len(X) + be[:, 1], key)
        if not missing.any():
            break
        if it == max_fix:
            raise RuntimeError(f"quadtree mesh: {missing.sum()} boundary segments not conforming")
        # split the missing segments at their midpoints and retry
        k = np.nonzero(missing)[0]
        newB = 0.5 * (B[bedges[k, 0]] + B[bedges[k, 1]])
        order = np.argsort(np.concatenate([np.arange(nb), k + 0.5]), kind="stable")
        B = np.vstack([B, newB])[order]
        beid = np.concatenate([beid, beid[k]])[order]
    x = X[T]
    det = (x[:, 1, 0] - x[:, 0, 0]) * (x[:, 2, 1] - x[:, 0, 1]) - (x[:, 1, 1] - x[:, 0, 1]) * (x[:, 2, 0] - x[:, 0, 0])
    T[det < 0] = T[det < 0][:, [0, 2, 1]]
    if np.any(np.abs(det) < 1e-12 * H * H * 4.0 ** (-size.ref)) or abs(0.5 * np.abs(det).sum() / dom.area() - 1) > 1e-10:
        raise RuntimeError("quadtree mesh: degenerate triangles or wrong area")
    # drop unused points (cannot happen with a conforming triangulation, kept for safety)
    used = np.zeros(len(X), bool)
    used[T.ravel()] = True
    if not used.all():
        new = np.cumsum(used) - 1
        if not used[:len(B)].all():
            raise RuntimeError("quadtree mesh: unused boundary node")
        X, T = X[used], new[T]
    nb = len(B)
    return Mesh(X, T, np.stack([np.arange(nb), (np.arange(nb) + 1) % nb], 1), beid)


_MESH_CACHE = {}


def build_mesh(dom: SlopeDomain, size: MeshSize | None = None, method="quadtree"):
    """triangle mesh of the slope box, method 'quadtree' (default) or 'blocks'; cached (alpha, bc and the
    polynomial order do not change the mesh, which is shared read-only by the solutions)"""
    size = size or MeshSize()
    # The mesh is generated for H = 1 and scaled by H: the Delaunay triangulation of the (cocircular) quadtree
    # points depends on round-off, so generating it at the actual H gave a different choice of diagonals for
    # each H (e.g. 3518 of 8469 triangles differ between H = 1 and H = 5 at beta = 60, ref 1) and a field that
    # was not exactly the scaled one (0.1 % RMS in the force near the slope).  Now u / (gamma_w H) and f / gamma_w
    # are exactly functions of (x / H, y / H), as the similarity argument H_crit = Gamma(H) H requires.
    hwr = dom.hw / dom.H
    key = (method, dom.beta_deg, hwr, dom.left, dom.right, dom.depth,
           size.h0, size.hs, size.grade, size.hmax, size.ref, size.uniform)
    if key not in _MESH_CACHE:
        if len(_MESH_CACHE) > 24:
            _MESH_CACHE.clear()
        unit = SlopeDomain(dom.beta_deg, 1.0, hwr, dom.left, dom.right, dom.depth)
        if method == "quadtree":
            _MESH_CACHE[key] = _quadtree_mesh(unit, size)
        elif method == "blocks":
            _MESH_CACHE[key] = _block_mesh(unit, size)
        else:
            raise ValueError(f"unknown mesh method {method!r}")
    m = _MESH_CACHE[key]
    if dom.H == 1.0:
        return m
    return Mesh(m.X * dom.H, m.T, m.edges, m.eid)


def _block_mesh(dom: SlopeDomain, size: MeshSize | None = None):
    """Three-block structured triangle mesh of the slope box (method 'blocks', see module docstring)."""
    size = size or MeshSize()
    sf = lambda x, y: size(dom, x, y)
    H, yb = dom.H, dom.yb
    # interface lines from O and T in the bisector direction (dx/dy = -tan(beta/2)), feet on the bottom
    tb2 = np.tan(0.5 * dom.beta)
    xOf = max(-yb * tb2, dom.xl + 0.25 * (0.0 - dom.xl))
    xTf = xOf + H / dom.sb
    xTf = min(xTf, dom.xr - 0.25 * (dom.xr - dom.xT))
    Of, Tf = np.array([xOf, yb]), np.array([xTf, yb])
    TL, BL, TR, BR = np.array([dom.xl, 0.0]), np.array([dom.xl, yb]), np.array([dom.xr, H]), np.array([dom.xr, yb])
    O, T = dom.O, dom.T
    foc = [O, T] + ([dom.W] if dom.has_W else [])
    # natural counts
    _, nA1 = _edge_nodes(TL, O, sf, foci=foc)
    _, nA2 = _edge_nodes(BL, Of, sf, foci=foc)
    _, nB1 = _edge_nodes(O, T, sf, foci=foc)
    _, nB2 = _edge_nodes(Of, Tf, sf, foci=foc)
    _, nC1 = _edge_nodes(T, TR, sf, foci=foc)
    _, nC2 = _edge_nodes(Tf, BR, sf, foci=foc)
    nvs = [_edge_nodes(P, Q, sf, foci=foc)[1] for P, Q in ((TL, BL), (O, Of), (T, Tf), (TR, BR))]
    nA, nB, nC, nV = max(nA1, nA2), max(nB1, nB2, 2 if dom.has_W else 1), max(nC1, nC2), max(nvs)
    crest = _edge_nodes(TL, O, sf, nA, foc)[0]
    botA = _edge_nodes(BL, Of, sf, nA, foc)[0]
    if dom.has_W:   # the water-line point must be a node: split the face at W
        sW = dom.hw / dom.H
        nB1a = max(1, int(round(nB * sW)))
        nB1a = min(nB1a, nB - 1)
        fa = _edge_nodes(O, dom.W, sf, nB1a, foc)[0]
        fb = _edge_nodes(dom.W, T, sf, nB - nB1a, foc)[0]
        face = np.vstack([fa, fb[1:]])
    else:
        face = _edge_nodes(O, T, sf, nB, foc)[0]
    botB = _edge_nodes(Of, Tf, sf, nB, foc)[0]
    toe = _edge_nodes(T, TR, sf, nC, foc)[0]
    botC = _edge_nodes(Tf, BR, sf, nC, foc)[0]
    lside = _edge_nodes(TL, BL, sf, nV, foc)[0]
    oline = _edge_nodes(O, Of, sf, nV, foc)[0]
    tline = _edge_nodes(T, Tf, sf, nV, foc)[0]
    rside = _edge_nodes(TR, BR, sf, nV, foc)[0]
    blocks = [
        _block_nodes(crest, botA, lside, oline),    # A
        _block_nodes(face, botB, oline, tline),     # B
        _block_nodes(toe, botC, tline, rside),      # C
    ]
    # global numbering (shared edges are bitwise identical arrays)
    allX = np.concatenate([b.reshape(-1, 2) for b in blocks])
    key = np.round(allX / (H * 1e-11)).astype(np.int64)
    _, first, inv = np.unique(key, axis=0, return_index=True, return_inverse=True)
    inv = inv.ravel()
    order = np.argsort(first)                       # keep a deterministic, block-wise ordering
    remap = np.empty_like(order)
    remap[order] = np.arange(len(order))
    X = allX[first[order]]
    gid = remap[inv]
    tris, edges, eids = [], [], []
    off = 0
    for bi, b in enumerate(blocks):
        ni, nj = b.shape[:2]
        G = gid[off:off + ni * nj].reshape(ni, nj)
        off += ni * nj
        a, bq, c, d = G[:-1, :-1], G[1:, :-1], G[1:, 1:], G[:-1, 1:]
        xa, xb, xc, xd = X[a], X[bq], X[c], X[d]
        diag_ac = np.sum((xa - xc) ** 2, -1) <= np.sum((xb - xd) ** 2, -1)
        t1 = np.where(diag_ac[..., None], np.stack([a, bq, c], -1), np.stack([a, bq, d], -1))
        t2 = np.where(diag_ac[..., None], np.stack([a, c, d], -1), np.stack([bq, c, d], -1))
        tris += [t1.reshape(-1, 3), t2.reshape(-1, 3)]
        top_id = (CREST, FACE, TOE)[bi]
        edges += [np.stack([G[:-1, 0], G[1:, 0]], 1), np.stack([G[:-1, -1], G[1:, -1]], 1)]
        eids += [np.full(ni - 1, top_id), np.full(ni - 1, BASE)]
        if bi == 0:
            edges.append(np.stack([G[0, :-1], G[0, 1:]], 1))
            eids.append(np.full(nj - 1, LEFT))
        if bi == 2:
            edges.append(np.stack([G[-1, :-1], G[-1, 1:]], 1))
            eids.append(np.full(nj - 1, RIGHT))
    T = np.concatenate(tris)
    x = X[T]
    det = (x[:, 1, 0] - x[:, 0, 0]) * (x[:, 2, 1] - x[:, 0, 1]) - (x[:, 1, 1] - x[:, 0, 1]) * (x[:, 2, 0] - x[:, 0, 0])
    if np.any(np.abs(det) < 1e-14 * H * H):
        raise RuntimeError("degenerate triangle in the block mesh")
    neg = det < 0
    T[neg] = T[neg][:, [0, 2, 1]]
    return Mesh(X, T, np.concatenate(edges), np.concatenate(eids))


# ======================================================================================================
# finite elements (P1 / P2 Lagrange on straight triangles)
# ======================================================================================================
# 7-point degree-5 rule (barycentric coordinates, weights sum to 1)
_A1, _B1, _A2, _B2 = 0.0597158717, 0.4701420641, 0.7974269853, 0.1012865073
Q7_L = np.array([[1 / 3, 1 / 3, 1 / 3], [_A1, _B1, _B1], [_B1, _A1, _B1], [_B1, _B1, _A1],
                 [_A2, _B2, _B2], [_B2, _A2, _B2], [_B2, _B2, _A2]])
Q7_W = np.array([0.225] + [0.1323941527] * 3 + [0.1259391805] * 3)
GL3_T = 0.5 * (1 + np.array([-np.sqrt(0.6), 0.0, np.sqrt(0.6)]))
GL3_W = np.array([5 / 18, 8 / 18, 5 / 18])


def _lambda_grads(X, T):
    x = X[T]
    det = (x[:, 1, 0] - x[:, 0, 0]) * (x[:, 2, 1] - x[:, 0, 1]) - (x[:, 1, 1] - x[:, 0, 1]) * (x[:, 2, 0] - x[:, 0, 0])
    e = np.stack([x[:, 2] - x[:, 1], x[:, 0] - x[:, 2], x[:, 1] - x[:, 0]], axis=1)    # edge opposite i
    G = np.stack([-e[:, :, 1], e[:, :, 0]], axis=2) / det[:, None, None]                # grad lambda_i
    return G, 0.5 * det


def _p2_grads(L, G):
    """gradients of the 6 P2 shape functions at barycentric point L (3,) for all elements; (ne, 6, 2).
    local order: vertices 0, 1, 2, edges (0,1), (1,2), (2,0)."""
    out = np.empty((G.shape[0], 6, 2))
    for i in range(3):
        out[:, i] = (4 * L[i] - 1) * G[:, i]
    for k, (i, j) in enumerate(((0, 1), (1, 2), (2, 0))):
        out[:, 3 + k] = 4 * (L[i] * G[:, j] + L[j] * G[:, i])
    return out


def _p2_vals(L):
    v = [L[0] * (2 * L[0] - 1), L[1] * (2 * L[1] - 1), L[2] * (2 * L[2] - 1),
         4 * L[0] * L[1], 4 * L[1] * L[2], 4 * L[2] * L[0]]
    return np.array(v)


class FESpace:
    """P1 or P2 Lagrange space on a Mesh; dof coordinates Xd; element dofs Ed (ne, 3 or 6)."""

    def __init__(self, mesh: Mesh, order=1):
        if order not in (1, 2):
            raise ValueError("order must be 1 or 2")
        self.mesh, self.order = mesh, order
        X, T = mesh.X, mesh.T
        self.G, self.area = _lambda_grads(X, T)
        if order == 1:
            self.Xd, self.Ed = X, T
            self.bdofs = lambda e: e                   # boundary edge -> dofs (end, end)
        else:
            E = np.concatenate([T[:, [0, 1]], T[:, [1, 2]], T[:, [2, 0]]])
            Es = np.sort(E, axis=1)
            uniq, inv = np.unique(Es, axis=0, return_inverse=True)
            inv = inv.ravel()
            ne = len(T)
            mid = len(X) + inv.reshape(3, ne).T
            self.Ed = np.concatenate([T, mid], axis=1)
            self.Xd = np.concatenate([X, 0.5 * (X[uniq[:, 0]] + X[uniq[:, 1]])])
            self._uniq = uniq
        self.ndof = len(self.Xd)

    def edge_midpoints(self, edges):
        """global midpoint dof of each boundary edge (P2)"""
        es = np.sort(edges, axis=1)
        # vectorised lookup through searchsorted on the lexicographically sorted unique edges
        key_u = self._uniq[:, 0].astype(np.int64) * (len(self.mesh.X) + 1) + self._uniq[:, 1]
        key_e = es[:, 0].astype(np.int64) * (len(self.mesh.X) + 1) + es[:, 1]
        pos = np.searchsorted(key_u, key_e)
        assert np.all(key_u[pos] == key_e)
        return len(self.mesh.X) + pos

    def stiffness(self, kx, ky):
        G, A = self.G, self.area
        Kd = np.array([kx, ky])
        if self.order == 1:
            Ke = A[:, None, None] * np.einsum("eik,k,ejk->eij", G, Kd, G)
        else:
            Ke = 0.0
            for L in ([0.5, 0.5, 0.0], [0.0, 0.5, 0.5], [0.5, 0.0, 0.5]):
                B = _p2_grads(np.array(L), G)
                Ke = Ke + (A / 3.0)[:, None, None] * np.einsum("eik,k,ejk->eij", B, Kd, B)
        nl = self.Ed.shape[1]
        rows = np.repeat(self.Ed, nl, axis=1).ravel()
        cols = np.tile(self.Ed, (1, nl)).ravel()
        return sp.csr_matrix((Ke.ravel(), (rows, cols)), shape=(self.ndof, self.ndof))

    def load(self, f):
        """int f phi_i (7-point rule); f(x, y) vectorised"""
        X, T = self.mesh.X, self.mesh.T
        x = X[T]
        b = np.zeros(self.ndof)
        for L, w in zip(Q7_L, Q7_W):
            xq = np.einsum("k,ekd->ed", L, x)
            fq = f(xq[:, 0], xq[:, 1]) * w * self.area
            phi = L if self.order == 1 else _p2_vals(L)
            np.add.at(b, self.Ed, fq[:, None] * phi[None, :])
        return b

    def neumann(self, edges, g):
        """int_edges g phi_i ds (3-point Gauss); g(x, y, nx, ny) with the outward normal computed by the
        caller-independent rule: n = rotation of the edge tangent pointing away from the interior."""
        X = self.mesh.X
        b = np.zeros(self.ndof)
        P, Q = X[edges[:, 0]], X[edges[:, 1]]
        ell = np.linalg.norm(Q - P, axis=1)
        if self.order == 2:
            mids = self.edge_midpoints(edges)
        for t, w in zip(GL3_T, GL3_W):
            xq = P + t * (Q - P)
            gq = g(xq[:, 0], xq[:, 1]) * w * ell
            if self.order == 1:
                np.add.at(b, edges[:, 0], gq * (1 - t))
                np.add.at(b, edges[:, 1], gq * t)
            else:
                np.add.at(b, edges[:, 0], gq * (1 - t) * (1 - 2 * t))
                np.add.at(b, edges[:, 1], gq * t * (2 * t - 1))
                np.add.at(b, mids, gq * 4 * t * (1 - t))
        return b

    def boundary_dofs(self, ids):
        m = np.isin(self.mesh.eid, ids)
        e = self.mesh.edges[m]
        d = np.unique(e)
        if self.order == 2:
            d = np.concatenate([d, self.edge_midpoints(e)])
        return np.unique(d)


def solve_dirichlet(space: FESpace, A, dofs, values, rhs=None):
    """solve A u = rhs with u[dofs] = values; returns u"""
    u = np.zeros(space.ndof)
    u[dofs] = values
    b = (np.zeros(space.ndof) if rhs is None else rhs.copy()) - A @ u
    free = np.ones(space.ndof, bool)
    free[dofs] = False
    Aff = A[free][:, free].tocsc()
    u[free] = spla.spsolve(Aff, b[free], permc_spec="COLAMD")     # (MMD_AT_PLUS_A is erratic on P1 meshes)
    return u


# ======================================================================================================
# the hydraulic problem of the paper and the seepage-force field object
# ======================================================================================================
class FESeepage:
    """FE solution u'_FE of Eqs. (20)-(21) and seepage force field f = -grad u'_FE.

    Parameters
    ----------
    beta_deg, H, hw, alpha : slope angle (deg), height, water level after drawdown (0 <= hw <= H), k_h / k_v
    kh, gamma_w            : horizontal permeability (only scales J and v) and unit weight of water
    left, right, depth     : FE box in units of H (left of O, right of T, below T); default: the paper's box
    bc                     : far-side conditions, preset name of BC_PRESETS or dict side -> neumann|zero|toe;
                             default "zero_lb" (= "paper": u = 0 on the left side and the bottom, right no-flow)
    order, ref, size       : P1/P2 (default P2), refinement level, MeshSize (default MeshSize(ref=ref))
    mesh, mesh_method      : prebuilt Mesh, or the generator "quadtree" (default) / "blocks"
    """

    def __init__(self, beta_deg, H=1.0, hw=None, alpha=1.0, kh=1.0, gamma_w=9.81, left=PAPER_BOX["left"],
                 right=PAPER_BOX["right"], depth=PAPER_BOX["depth"], bc="zero_lb", order=2, ref=1, size=None,
                 mesh=None, mesh_method="quadtree"):
        t0 = time.time()
        self.dom = SlopeDomain(beta_deg, H, hw, left, right, depth)
        self.H, self.hw, self.alpha = self.dom.H, self.dom.hw, float(alpha)
        if self.alpha < 1.0:
            raise ValueError("alpha = k_h / k_v must be >= 1")
        self.kh, self.kv, self.gamma_w = float(kh), float(kh) / self.alpha, float(gamma_w)
        self.bc = _bc_dict(bc)
        self.order = order
        self.mesh = mesh if mesh is not None else build_mesh(self.dom, size or MeshSize(ref=ref), mesh_method)
        self.space = FESpace(self.mesh, order)
        self.A = self.space.stiffness(self.kh, self.kv)
        sp_ = self.space
        # Dirichlet dofs: far sides first, the ground surface last (it wins at shared corners)
        dofs, vals = [], []
        for side, kind in self.bc.items():
            if kind == "neumann":
                continue
            d = sp_.boundary_dofs([FAR_IDS[side]])
            v = 0.0 if kind == "zero" else -self.gamma_w * self.hw
            if kind not in ("zero", "toe"):
                raise ValueError(f"unknown boundary condition {kind!r}")
            dofs.append(d)
            vals.append(np.full(len(d), v))
        ds = sp_.boundary_dofs(list(SURFACE_IDS))
        dofs.append(ds)
        vals.append(self.dom.surface_u(sp_.Xd[ds, 0], sp_.Xd[ds, 1], self.gamma_w))
        alld = np.concatenate(dofs)
        allv = np.concatenate(vals)
        # keep the LAST value of each dof (surface data override far-side data at corners)
        _, idx = np.unique(alld[::-1], return_index=True)
        idx = len(alld) - 1 - idx
        self.ddofs, self.dvals = alld[idx], allv[idx]
        self.uh = solve_dirichlet(sp_, self.A, self.ddofs, self.dvals)
        self.J = 0.5 * float(self.uh @ (self.A @ self.uh))
        self.solve_time = time.time() - t0
        self._tri = None
        self._trifinder = None
        if order == 1:
            self._grad_tri = np.einsum("eik,ei->ek", self.space.G, self.uh[self.mesh.T])

    # ------------------------------------------------------------------ scalar outputs
    def J_normalized(self):
        """J(u'_FE) / (k_h H^2 gamma_w^2), the normalisation of Fig. 5"""
        return self.J / (self.kh * self.H ** 2 * self.gamma_w ** 2)

    def pore_pressure_dofs(self):
        """total pore pressure p = u + gamma_w y at the dofs"""
        return self.uh + self.gamma_w * self.space.Xd[:, 1]

    def min_pore_pressure(self):
        p = self.pore_pressure_dofs()
        i = int(np.argmin(p))
        return float(p[i]), self.space.Xd[i].copy()

    def boundary_flux(self, ids):
        """outflow int (-K grad u) . n ds through the boundary parts `ids` from the residual (reaction)
        of the discrete equations at their Dirichlet dofs: R = A u (= inflow K grad u . n integrated)"""
        R = self.A @ self.uh
        d = self.space.boundary_dofs(ids)
        return -float(R[d].sum())

    # ------------------------------------------------------------------ point evaluation
    # points of the closed soil domain within LOCATE_TOL * H of the boundary are always located (see _locate)
    LOCATE_TOL = 1e-9

    def _locate(self, x, y):
        """element index of each point (-1 outside the FE box).  The trapezoid-map finder misses most points
        that lie ON an inclined boundary edge (round-off puts about half of them a few ulps outside, and the
        finder rejects some of the others), e.g. 1890 of 2000 random points of the face returned f = 0 and
        u = nan before this fix.  Points within LOCATE_TOL * H of the soil that are not found are therefore
        searched again after a shift of 4 LOCATE_TOL * H into the soil (inward normals of the face, of the
        horizontal and vertical sides and of the box corners)."""
        if self._trifinder is None:
            import matplotlib.tri as mtri
            self._tri = mtri.Triangulation(self.mesh.X[:, 0], self.mesh.X[:, 1], self.mesh.T)
            self._trifinder = self._tri.get_trifinder()
        e = np.asarray(self._trifinder(x, y))
        bad = ~(np.isfinite(x) & np.isfinite(y))          # the finder returns an element for y = nan
        if np.any(bad):
            e = np.where(bad, -1, e)
        miss = np.nonzero((e < 0) & ~bad)[0]
        if len(miss):
            tol = self.LOCATE_TOL * self.H
            xm, ym = x[miss], y[miss]
            with np.errstate(invalid="ignore"):
                near = self.dom.in_box(xm, ym, tol) | (self.dom.dist_boundary(xm, ym) <= tol)
            miss = miss[near]
            if len(miss):
                e = e.copy()
                s = 4.0 * tol
                r2 = np.sqrt(0.5)
                dirs = ((-self.dom.sb, self.dom.cb), (0.0, 1.0), (1.0, 0.0), (-1.0, 0.0), (0.0, -1.0),
                        (r2, r2), (-r2, r2), (r2, -r2), (-r2, -r2))
                for dx, dy in dirs:
                    ee = np.asarray(self._trifinder(x[miss] + s * dx, y[miss] + s * dy))
                    ok = ee >= 0
                    e[miss[ok]] = ee[ok]
                    miss = miss[~ok]
                    if not len(miss):
                        break
        return e

    def _bary(self, e, x, y):
        X = self.mesh.X
        x0 = X[self.mesh.T[e, 0]]
        G = self.space.G[e]
        l1 = G[:, 1, 0] * (x - x0[:, 0]) + G[:, 1, 1] * (y - x0[:, 1])
        l2 = G[:, 2, 0] * (x - x0[:, 0]) + G[:, 2, 1] * (y - x0[:, 1])
        return np.stack([1.0 - l1 - l2, l1, l2])

    def grad(self, x, y):
        """grad u'_FE (kPa/m) at arbitrary points; 0 outside the mesh"""
        x, y = np.broadcast_arrays(np.asarray(x, float), np.asarray(y, float))
        shape = x.shape
        xf, yf = x.ravel(), y.ravel()
        e = self._locate(xf, yf)
        gx, gy = np.zeros(xf.shape), np.zeros(xf.shape)
        ok = e >= 0
        if np.any(ok):
            eo = e[ok]
            if self.order == 1:
                g = self._grad_tri[eo]
            else:
                L = self._bary(eo, xf[ok], yf[ok])
                G = self.space.G[eo]
                ue = self.uh[self.space.Ed[eo]]
                g = np.zeros((len(eo), 2))
                for i in range(3):
                    g += ((4 * L[i] - 1) * ue[:, i])[:, None] * G[:, i]
                for k, (i, j) in enumerate(((0, 1), (1, 2), (2, 0))):
                    g += (4 * ue[:, 3 + k])[:, None] * (L[i][:, None] * G[:, j] + L[j][:, None] * G[:, i])
            gx[ok], gy[ok] = g[:, 0], g[:, 1]
        return gx.reshape(shape), gy.reshape(shape)

    def force(self, x, y):
        """seepage force f = -grad u'_FE (kN/m^3), paper coordinates; 0 outside the soil/FE box"""
        gx, gy = self.grad(x, y)
        return -gx, -gy

    def velocity(self, x, y):
        gx, gy = self.grad(x, y)
        return -self.kh * gx, -self.kv * gy

    def u(self, x, y):
        """u'_FE at arbitrary points (nan outside the FE box)"""
        x, y = np.broadcast_arrays(np.asarray(x, float), np.asarray(y, float))
        shape = x.shape
        xf, yf = x.ravel(), y.ravel()
        e = self._locate(xf, yf)
        out = np.full(xf.shape, np.nan)
        ok = e >= 0
        if np.any(ok):
            eo = e[ok]
            L = self._bary(eo, xf[ok], yf[ok])
            ue = self.uh[self.space.Ed[eo]]
            if self.order == 1:
                out[ok] = np.sum(L.T * ue, 1)
            else:
                out[ok] = np.sum(_p2_vals(L).T * ue, 1)
        return out.reshape(shape)

    def summary(self):
        q = self.mesh.quality()
        pmin, xp = self.min_pore_pressure()
        return (f"beta={self.dom.beta_deg:5.1f} alpha={self.alpha:5.2f} hw/H={self.hw / self.H:5.3f} "
                f"box(L,R,D)/H=({self.dom.left:g},{self.dom.right:g},{self.dom.depth:g}) bc={self.bc} P{self.order} "
                f"nv={self.mesh.nv} ntri={len(self.mesh.T)} ndof={self.space.ndof} angles=[{q['min_angle']:.1f},"
                f"{q['max_angle']:.1f}] J/(kh H^2 gw^2)={self.J_normalized():.6f} min p={pmin:.3e} at {xp}")


class NoSeepage:
    """f = 0 everywhere (h_w = 0)."""

    def force(self, x, y):
        x, y = np.broadcast_arrays(np.asarray(x, float), np.asarray(y, float))
        return np.zeros(x.shape), np.zeros(x.shape)

    def u(self, x, y):
        x, y = np.broadcast_arrays(np.asarray(x, float), np.asarray(y, float))
        return np.zeros(x.shape)


# ======================================================================================================
# verification helpers
# ======================================================================================================
def fe_errors(space: FESpace, uh, uex, gex, kx, ky):
    """L2 error and energy-norm error ||grad(u_h - u)||_K with the 7-point rule"""
    X, T = space.mesh.X, space.mesh.T
    x = X[T]
    ue = uh[space.Ed]
    eL2 = eE = nE = 0.0
    for L, w in zip(Q7_L, Q7_W):
        xq = np.einsum("k,ekd->ed", L, x)
        if space.order == 1:
            uq = ue @ L
            gq = np.einsum("ei,eik->ek", ue, space.G)
        else:
            uq = ue @ _p2_vals(L)
            gq = np.einsum("ei,eik->ek", ue, _p2_grads(L, space.G))
        gx, gy = gex(xq[:, 0], xq[:, 1])
        wa = w * space.area
        eL2 += np.sum(wa * (uq - uex(xq[:, 0], xq[:, 1])) ** 2)
        eE += np.sum(wa * (kx * (gq[:, 0] - gx) ** 2 + ky * (gq[:, 1] - gy) ** 2))
        nE += np.sum(wa * (kx * gx ** 2 + ky * gy ** 2))
    return np.sqrt(eL2), np.sqrt(eE), np.sqrt(nE)


def _rate(e, h):
    e, h = np.asarray(e), np.asarray(h)
    r = np.full(len(e), np.nan)
    r[1:] = np.log(e[1:] / e[:-1]) / np.log(h[1:] / h[:-1])
    return r


def _load_laplace_check():
    path = os.path.join(REPO_DIR, "docs", "slope_mohr_coulomb", "scripts", "laplace_check.py")
    if not os.path.exists(path):
        return None
    spec = importlib.util.spec_from_file_location("laplace_check", path)
    mod = importlib.util.module_from_spec(spec)
    argv = sys.argv
    sys.argv = [path]            # the module reads sys.argv at import time
    try:
        spec.loader.exec_module(mod)
    finally:
        sys.argv = argv
    return mod


def _load_fig5_paper():
    path = os.path.join(DATA_DIR, "fig5_vector_fill_polygons.csv")
    if not os.path.exists(path):
        return None
    arr = np.loadtxt(path, delimiter=",", comments="#")
    return {(int(a), int(round(b))): (s, d) for a, b, s, d in arr}


def _soil_samples(dom, rmax, step):
    """grid points of the soil within distance rmax of O (excluding a thin layer at the surface)"""
    g = np.arange(-rmax, rmax + 1e-9, step)
    x, y = np.meshgrid(g + 0.5 * step, np.arange(0.0, rmax + 1e-9, step) + 0.5 * step)
    x, y = x.ravel(), y.ravel()
    keep = (np.hypot(x, y) < rmax) & (y > dom.ground(x) + 1e-6 * dom.H)
    return x[keep], y[keep]


# ---- Fig. 4 of the paper (measured on the 1034 x 1630 crop fig4.png of PDF page 8, panel b) -------------
# box edges and slope in pixels: left 42.5, right 824.5, crest 604.5, toe ground 617.5, bottom 1002,
# crest edge O at x = 684.5, toe at x ~ 696 (beta ~ 45-50 deg, H ~ 13 px);
# rows where the 9 iso-lines -0.98 k (k = 9 .. 1) reach the right side (column 821):
FIG4B_PX = dict(left=42.5, right=824.5, crest=604.5, toe=617.5, bottom=1002.0, xO=684.5)
FIG4B_RIGHT_ROWS = [(-8.83, 637.5), (-7.85, 658.5), (-6.87, 681.0), (-5.89, 707.5), (-4.91, 736.5),
                    (-3.92, 772.0), (-2.94, 815.5), (-1.96, 868.5), (-0.98, 931.5)]
COARSE_MESHES = ((0.25, 0.35, 8.0), (0.15, 0.3, 6.0), (0.1, 0.3, 6.0))    # (hc, grade, hmax) like Fig. 4a
FIG4_IMAGE_DEFAULT = "/tmp/claude-0/-home-user-neopz-master/339d491a-b062-561e-a71b-775db071f418/scratchpad/pdf/fig4.png"


PAPER_PDF_DEFAULT = os.path.join(os.path.dirname(FIG4_IMAGE_DEFAULT), "paper.pdf")


def measure_fig4_native(pdf):
    """Box of Fig. 4a/4b on the native image embedded in the PDF (page 8, 1000 x 1728 px grey JPEG at 300 ppi,
    extracted with pdfimages): box columns, crest / toe-ground / bottom rows of panel (b) and the extents in H for
    a box width of 61 H (left 50 H + face 1 H (beta = 45) + right 10 H) and for H = crest-to-toe distance."""
    import shutil
    import subprocess
    import tempfile
    from PIL import Image
    if not (os.path.exists(pdf) and shutil.which("pdfimages")):
        return None
    with tempfile.TemporaryDirectory() as tmp:
        subprocess.run(["pdfimages", "-j", "-f", "8", "-l", "8", pdf, os.path.join(tmp, "f4")], check=True)
        files = sorted(f for f in os.listdir(tmp) if f.startswith("f4"))
        if not files:
            return None
        im = np.array(Image.open(os.path.join(tmp, files[0])).convert("L")).astype(float)
    if im.shape != (1728, 1000):
        return None
    dark = im < 150
    b = dark[540:1100]                                         # panel (b)
    cols = np.nonzero(b.mean(0) > 0.6)[0]
    left, right = cols[cols < 500].mean(), cols[cols > 500].mean()
    prof = im[540:1100, 20:700].mean(1)
    rows = np.nonzero(prof < 200)[0] + 540                     # crest line (top) and bottom line
    crest, bottom = rows[rows < 800].mean(), rows[rows > 800].mean()
    prof_r = im[560:640, 860:940].mean(1)
    rr = np.nonzero(prof_r < 200)[0] + 560                     # toe-ground line = first dark run (then iso-lines)
    toe = rr[rr <= rr[0] + 2].mean()
    # crest edge O: right end of the crest line
    xO = np.nonzero(dark[int(round(crest)), :int(right) - 50])[0].max()
    out = dict(left_col=left, right_col=right, crest_row=crest, toe_row=toe, bottom_row=bottom, xO_col=float(xO))
    for name, Hpx in (("width61", (right - left) / 61.0), ("crest_to_toe", toe - crest)):
        out[name] = dict(H_px=Hpx, left=(xO - left) / Hpx, right_of_O=(right - xO) / Hpx,
                         depth_below_crest=(bottom - crest) / Hpx, depth_below_toe=(bottom - toe) / Hpx)
    return out


def measure_fig4b(path):
    """re-measure the box of Fig. 4b and the right-side ends of the iso-lines from the page crop"""
    from PIL import Image
    im = np.array(Image.open(path).convert("L")).astype(float)
    dark = im < 160
    rows = np.arange(590, 1010)
    cols = np.arange(30, 840)
    sub = dark[np.ix_(rows, cols)]
    vcol = cols[sub.mean(0) > 0.9]
    hrow = rows[dark[np.ix_(rows, np.arange(100, 600))].mean(1) > 0.9]
    hrow_r = rows[dark[np.ix_(rows, np.arange(720, 820))].mean(1) > 0.9]
    col = dark[618:1002, 821]
    idx = np.nonzero(col)[0] + 618
    runs = np.split(idx, np.nonzero(np.diff(idx) > 1)[0] + 1)
    ends = [0.5 * (r[0] + r[-1]) for r in runs if len(r)]
    return dict(vertical_lines=vcol.tolist(), crest_rows=hrow.tolist(), right_top_rows=hrow_r.tolist(),
                right_side_crossings=ends)


# ======================================================================================================
# checks
# ======================================================================================================
def check_mms(log, quick):
    log("=" * 110)
    log("1) Manufactured solution, anisotropic K = diag(1, 1/4) (alpha = 4), beta = 60, h_w = 0.6 H, box 2/1.5/1.5 H:")
    log("   u = cos(1.1x+0.3) sinh(0.6y+0.2) + 0.3xy, f = -div K grad u, Dirichlet on crest/face/toe, NON-zero")
    log("   Neumann flux g = K grad u . n on the left, bottom and right sides; uniform-size quadtree meshes (not nested:")
    log("   rates with respect to h_eff = sqrt(|Omega| / n_vertices); expected P1: 2 (L2), 1 (energy); P2: 3, 2)")
    log("=" * 110)
    kx, ky = 1.0, 0.25
    uex = lambda x, y: np.cos(1.1 * x + 0.3) * np.sinh(0.6 * y + 0.2) + 0.3 * x * y
    gex = lambda x, y: (-1.1 * np.sin(1.1 * x + 0.3) * np.sinh(0.6 * y + 0.2) + 0.3 * y,
                        0.6 * np.cos(1.1 * x + 0.3) * np.cosh(0.6 * y + 0.2) + 0.3 * x)
    fsrc = lambda x, y: -(kx * (-1.21) + ky * 0.36) * np.cos(1.1 * x + 0.3) * np.sinh(0.6 * y + 0.2)
    dom = SlopeDomain(60.0, 1.0, 0.6, left=2.0, right=1.5, depth=1.5)
    normals = {LEFT: (-1.0, 0.0), BASE: (0.0, 1.0), RIGHT: (1.0, 0.0)}
    res = {}
    for order in (1, 2):
        es, hs, ns = [], [], []
        for r in range(4 if quick else 5):
            mesh = build_mesh(dom, MeshSize(uniform=0.25, ref=r))
            spc = FESpace(mesh, order)
            A = spc.stiffness(kx, ky)
            b = spc.load(fsrc)
            for i, n in normals.items():
                e = mesh.edges[mesh.eid == i]
                b += spc.neumann(e, lambda x, y, n=n: kx * gex(x, y)[0] * n[0] + ky * gex(x, y)[1] * n[1])
            d = spc.boundary_dofs(list(SURFACE_IDS))
            uh = solve_dirichlet(spc, A, d, uex(spc.Xd[d, 0], spc.Xd[d, 1]), rhs=b)
            eL2, eE, nE = fe_errors(spc, uh, uex, gex, kx, ky)
            es.append((eL2, eE))
            hs.append(np.sqrt(dom.area() / mesh.nv))
            ns.append(spc.ndof)
        es = np.array(es)
        rL, rE = _rate(es[:, 0], hs), _rate(es[:, 1], hs)
        log(f"  P{order}:  {'ref':>3} {'ndof':>8} {'h_eff':>8} {'L2 error':>11} {'rate':>5} {'energy err':>11} {'rate':>5}")
        for r in range(len(hs)):
            log(f"        {r:3d} {ns[r]:8d} {hs[r]:8.4f} {es[r, 0]:11.3e} {rL[r]:5.2f} {es[r, 1]:11.3e} {rE[r]:5.2f}")
        res[order] = (rL[-1], rE[-1])
    # patch test: linear u, mixed BCs, P1 and P2 must be exact (both mesh generators)
    lin = lambda x, y: 2.0 - 0.7 * x + 1.3 * y
    for order, method in ((1, "quadtree"), (2, "quadtree"), (1, "blocks"), (2, "blocks")):
        mesh = build_mesh(dom, MeshSize(ref=0, h0=0.05, hs=0.2, hmax=1.0), method)
        spc = FESpace(mesh, order)
        A = spc.stiffness(kx, ky)
        b = np.zeros(spc.ndof)
        for i, n in normals.items():
            e = mesh.edges[mesh.eid == i]
            b += spc.neumann(e, lambda x, y, n=n: np.full(np.shape(x), kx * (-0.7) * n[0] + ky * 1.3 * n[1]))
        d = spc.boundary_dofs(list(SURFACE_IDS))
        uh = solve_dirichlet(spc, A, d, lin(spc.Xd[d, 0], spc.Xd[d, 1]), rhs=b)
        log(f"  patch test P{order} (u linear, Dirichlet surface + Neumann far sides, graded {method} mesh): "
            f"max |u_h - u| = {np.abs(uh - lin(spc.Xd[:, 0], spc.Xd[:, 1])).max():.2e}")
    return res


def check_mesh(log, quick, plot):
    log("=" * 110)
    log("0) Meshes: quadtree conforming Delaunay generator (default) vs the block generator, box 50/10/30 H")
    log("=" * 110)
    log(f"  {'beta':>4} {'hw/H':>5} {'method':>8} {'ref':>3} {'nv':>7} {'ntri':>7} {'min ang':>7} {'max ang':>7} "
        f"{'h_min/H':>9} {'h_max/H':>7} {'#bnd':>5} {'time':>6}")
    for b in ((15, 45, 90) if quick else (15, 30, 45, 60, 75, 90)):
        for hwr in (1.0, 0.5):
            dom = SlopeDomain(b, 1.0, hwr)
            for method, refs in (("quadtree", (0, 1, 2)), ("blocks", (1,))):
                for r in refs:
                    t = time.time()
                    m = build_mesh(dom, MeshSize(ref=r), method)
                    dt = time.time() - t
                    q = m.quality()
                    d1, d2 = m.X[m.T[:, 1]] - m.X[m.T[:, 0]], m.X[m.T[:, 2]] - m.X[m.T[:, 0]]
                    area = 0.5 * np.abs(d1[:, 0] * d2[:, 1] - d1[:, 1] * d2[:, 0]).sum()
                    assert abs(area / dom.area() - 1) < 1e-10
                    log(f"  {b:4.0f} {hwr:5.2f} {method:>8} {r:3d} {m.nv:7d} {len(m.T):7d} {q['min_angle']:7.2f} "
                        f"{q['max_angle']:7.1f} {q['hmin']:9.2e} {q['hmax']:7.2f} {len(m.edges):5d} {dt:5.2f}s")
    if plot:
        _plot_meshes()


def _plot_meshes():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(3, 3, figsize=(15, 10.5), constrained_layout=True)
    cols = {CREST: "#2a78d6", FACE: "#eb6834", TOE: "#1baf7a", LEFT: "0.3", BASE: "0.3", RIGHT: "0.3"}
    for row, (b, hwr) in zip(axs, ((15, 0.5), (45, 1.0), (90, 0.5))):
        dom = SlopeDomain(b, 1.0, hwr)
        m = build_mesh(dom, MeshSize(ref=0))
        for ax, lim in zip(row, (None, (-6.0, 8.0, 6.0, -0.4), (-0.3, dom.xT + 0.3, 1.3, -0.1))):
            ax.triplot(m.X[:, 0], m.X[:, 1], m.T, lw=0.25, color="0.45")
            for i, c in cols.items():
                e = m.edges[m.eid == i]
                ax.plot(m.X[e].transpose(2, 1, 0)[0], m.X[e].transpose(2, 1, 0)[1], color=c, lw=1.4)
            ax.set_aspect("equal")
            if lim:
                ax.set_xlim(lim[0], lim[1])
                ax.set_ylim(lim[2], lim[3])
            else:
                ax.invert_yaxis()
            ax.set_title(f"beta = {b}, h_w = {hwr:g} H, ref 0 ({m.nv} nodes)", fontsize=9)
    fig.suptitle("quadtree conforming Delaunay meshes (paper coordinates, y down; x, y in units of H)", fontsize=10)
    os.makedirs(RESULTS_DIR, exist_ok=True)
    fig.savefig(os.path.join(RESULTS_DIR, "fe_meshes.png"), dpi=90)
    plt.close(fig)


def check_laplace_cross(log):
    log("=" * 110)
    log("2) Cross-check with the independent P1 solver docs/slope_mohr_coulomb/scripts/laplace_check.py")
    log("=" * 110)
    lc = _load_laplace_check()
    if lc is None:
        log("  laplace_check.py not found - skipped")
        return None
    out = {}
    # (a) its own problem: TriGMesh(2) of SlopeMohrCoulomb, head h = y on -3,-4,-6, no flow on -1,-2,-5
    X, T, Lb = lc.tri_gmesh(2)
    h_ref, dm, R, Aref = lc.solve_dirichlet(X, T, Lb, {-3, -4, -6}, lambda Y: Y[:, 1].copy())
    edges = np.array([e for e, i in Lb])
    eid = np.array([i for e, i in Lb])
    mesh = Mesh(X, np.asarray(T), edges, eid)
    spc = FESpace(mesh, 1)
    A = spc.stiffness(1.0, 1.0)
    d = spc.boundary_dofs([-3, -4, -6])
    h = solve_dirichlet(spc, A, d, X[d, 1])
    out["a"] = np.abs(h - h_ref).max()
    log(f"  (a) SlopeDrawdown steady head on TriGMesh(2) ({len(X)} nodes): max |h_this - h_laplace_check| = "
        f"{out['a']:.2e} m (h in [{h.min():.3f}, {h.max():.3f}])")
    # (b) this module's problem (alpha = 1, zero_lb and impermeable) solved by laplace_check.solve_dirichlet
    gw = 9.81
    for bc in ("zero_lb", "impermeable"):
        fe = FESeepage(45.0, 1.0, 1.0, 1.0, gamma_w=gw, bc=bc, order=1, ref=0)
        m = fe.mesh
        lines = [((int(a), int(b)), int(i)) for (a, b), i in zip(m.edges, m.eid)]
        ids = set(SURFACE_IDS) | {FAR_IDS[s] for s, k in fe.bc.items() if k == "zero"}
        surf = lambda Y: np.abs(Y[:, 1] - fe.dom.ground(Y[:, 0])) < 1e-9
        g = lambda Y: np.where(surf(Y), fe.dom.surface_u(Y[:, 0], Y[:, 1], gw), 0.0)
        u_lc, dm, R, Alc = lc.solve_dirichlet(m.X, m.T, lines, ids, g)
        J_lc = 0.5 * float(u_lc @ Alc.mv(u_lc))
        du = np.abs(u_lc - fe.uh).max()
        out["b_" + bc] = (du, J_lc / (gw ** 2), fe.J_normalized())
        log(f"  (b) beta=45, alpha=1, h_w=H, bc={bc:11s} mesh ref 0 ({m.nv} nodes): max |u_this - u_lc| = {du:.2e} kPa, "
            f"J/(kh H^2 gw^2): laplace_check {J_lc / gw ** 2:.8f}, this {fe.J_normalized():.8f}")
    # (c) anisotropy by the affine map y' = sqrt(alpha) y: the anisotropic problem (beta, H, alpha) is the ISOTROPIC
    # problem on the stretched slope beta' = atan(sqrt(alpha) tan beta), H' = sqrt(alpha) H, h_w' = sqrt(alpha) h_w,
    # gamma_w' = gamma_w / sqrt(alpha), box left 50 / sqrt(alpha), right 10 / sqrt(alpha), depth 30 (units of H'),
    # and J_aniso / (k_h H^2 gw^2) = J_iso' / (k_h H'^2 gw'^2) / sqrt(alpha)  (different, non-affine meshes)
    log("  (c) anisotropy vs the isotropic problem on the stretched slope (y' = sqrt(alpha) y), P2 ref 1, zero_lb:")
    for b, a, hwr in ((30.0, 4.0, 1.0), (60.0, 10.0, 1.0), (45.0, 4.0, 0.5)):
        fa = FESeepage(b, 1.0, hwr, a, bc="zero_lb", order=2, ref=1)
        sa = np.sqrt(a)
        bs = np.degrees(np.arctan(sa * np.tan(np.radians(b))))
        fi = FESeepage(bs, sa, sa * hwr, 1.0, gamma_w=9.81 / sa, bc="zero_lb", order=2, ref=1,
                       left=PAPER_BOX["left"] / sa, right=PAPER_BOX["right"] / sa, depth=PAPER_BOX["depth"])
        Ji = fi.J_normalized() / sa
        # the force: f_x = f'_x, f_y = sqrt(alpha) f'_y at (x, y' = sqrt(alpha) y)
        xs, ys = _soil_samples(fa.dom, 3.0, 0.1)
        fx, fy = fa.force(xs, ys)
        gx, gy = fi.force(xs, sa * ys)
        rel = np.sqrt(np.mean((fx - gx) ** 2 + (fy - sa * gy) ** 2) / np.mean(fx ** 2 + fy ** 2))
        out[f"c_{b}_{a}"] = (fa.J_normalized(), Ji, rel)
        log(f"      beta={b:4.0f} alpha={a:4.0f} hw/H={hwr}: J anisotropic {fa.J_normalized():.7f}, isotropic stretched "
            f"(beta' = {bs:.2f}) {Ji:.7f}, rel. diff {Ji / fa.J_normalized() - 1:+.1e}; force within 3H of O: RMS rel. "
            f"diff {100 * rel:.3f}%")
    return out


def check_convergence(log, quick):
    log("=" * 110)
    log("3) Mesh convergence of J(u'_FE)/(k_h H^2 gamma_w^2), box 50/10/30 H, bc zero_lb (u = 0 left + bottom)")
    log("   (J converges from above: the discrete space with the exact piecewise-linear Dirichlet data is conforming)")
    log("   quadtree meshes; last line: P2 on the independent 'blocks' mesh generator (ref 2)")
    log("=" * 110)
    cases = [(45, 1, 1.0), (90, 10, 1.0)] if quick else [(45, 1, 1.0), (90, 1, 1.0), (15, 1, 1.0), (30, 10, 1.0),
                                                         (60, 4, 0.5)]
    rows = []
    for b, a, hwr in cases:
        log(f"  beta={b} alpha={a} hw/H={hwr}")
        for order, refs in ((1, range(3 if quick else 4)), (2, range(2 if quick else 3))):
            Js, nd = [], []
            for r in refs:
                fe = FESeepage(b, 1.0, hwr, a, bc="zero_lb", order=order, ref=r)
                Js.append(fe.J_normalized())
                nd.append(fe.space.ndof)
                pmin, xp = fe.min_pore_pressure()
                rows.append((b, a, hwr, order, r, fe.space.ndof, Js[-1]))
                rate = np.nan
                if len(Js) >= 3:
                    rate = np.log((Js[-3] - Js[-2]) / (Js[-2] - Js[-1])) / np.log(nd[-1] / nd[-2])
                log(f"    P{order} ref {r}: ndof {fe.space.ndof:8d}  J = {Js[-1]:.7f}  dJ = "
                    f"{(Js[-1] - Js[-2]) if len(Js) > 1 else np.nan:+.2e}  rate(N) {rate:5.2f}  min p = {pmin:.2e} kPa"
                    f"  ({fe.solve_time:.1f} s)")
            if len(Js) >= 3:
                q = (Js[-2] - Js[-1]) / (Js[-3] - Js[-2])
                Jx = Js[-1] - (Js[-2] - Js[-1]) * q / (1 - q)
                log(f"    P{order} Aitken extrapolation J_inf = {Jx:.7f}")
        rb = 1 if quick else 2
        fb = FESeepage(b, 1.0, hwr, a, bc="zero_lb", order=2, ref=rb, mesh_method="blocks")
        q = fb.mesh.quality()
        log(f"    P2 blocks ref {rb}: ndof {fb.space.ndof:8d}  J = {fb.J_normalized():.7f}  (minus the last quadtree P2 value: "
            f"{fb.J_normalized() - Js[-1]:+.2e}; block-mesh angles [{q['min_angle']:.2f}, {q['max_angle']:.1f}] deg)")
        rows.append((b, a, hwr, -2, rb, fb.space.ndof, fb.J_normalized()))
    return rows


def check_sensitivity(log, quick):
    log("=" * 110)
    log("4) Sensitivity to the box extents and the far-side BC (P2 ref 1): J and the seepage force -grad u'_FE")
    log("   in the soil within 3H of O (grid 0.05H): rel = RMS|f - f_base| / RMS|f_base|, max = max|f - f_base| / max|f_base|")
    log("=" * 110)
    configs = [("base 50/10/30 zero_lb", dict()),
               ("depth 29 (30H below crest)", dict(depth=29.0)),
               ("left 49.4, depth 29.6 (raw px)", dict(left=49.4, depth=29.6, right=9.9)),
               ("left 25", dict(left=25.0)), ("left 100", dict(left=100.0)),
               ("right 5", dict(right=5.0)), ("right 20", dict(right=20.0)),
               ("depth 15", dict(depth=15.0)), ("depth 60", dict(depth=60.0)),
               ("all x2 (100/20/60)", dict(left=100.0, right=20.0, depth=60.0)),
               ("bc zero_b", dict(bc="zero_b")), ("bc zero_l", dict(bc="zero_l")),
               ("bc impermeable", dict(bc="impermeable")),
               ("bc impermeable 100/20/60", dict(bc="impermeable", left=100.0, right=20.0, depth=60.0)),
               ("bc toe_r (u=-gw hw right)", dict(bc="toe_r"))]
    if quick:
        configs = configs[:2] + configs[9:13]
    out = []
    for b, a in ((45, 1), (90, 10)):
        dom = SlopeDomain(b, 1.0, 1.0)
        xs, ys = _soil_samples(dom, 3.0, 0.05)
        base = None
        log(f"  beta={b} alpha={a} h_w=H  ({len(xs)} sample points)")
        for name, kw in configs:
            kw = dict(kw)
            bc = kw.pop("bc", "zero_lb")
            fe = FESeepage(b, 1.0, 1.0, a, bc=bc, order=2, ref=1, **kw)
            fx, fy = fe.force(xs, ys)
            if base is None:
                base = (fe.J_normalized(), fx, fy)
            dn = np.hypot(fx - base[1], fy - base[2])
            nb = np.hypot(base[1], base[2])
            rel = np.sqrt(np.mean(dn ** 2) / np.mean(nb ** 2))
            mx = dn.max() / nb.max()
            out.append((b, a, name, fe.J_normalized(), rel, mx))
            log(f"    {name:32s} J = {fe.J_normalized():.5f} ({100 * (fe.J_normalized() / base[0] - 1):+6.2f}%)  "
                f"force rel {100 * rel:6.2f}%  max {100 * mx:6.2f}%")
    return out


FIG5_VARIANTS = (   # (column label, FESeepage keyword arguments besides beta, H, hw, alpha)
    ("zero_lb", dict(bc="zero_lb")),
    ("zero_b", dict(bc="zero_b")),
    ("zero_l", dict(bc="zero_l")),
    ("impermeable", dict(bc="impermeable")),
    ("zero_lb_d29", dict(bc="zero_lb", depth=29.0)),            # 30 H below the CREST
    ("zero_lb_rO", dict(bc="zero_lb", right="from_O")),         # right side 10 H beyond O (not beyond the toe)
)


def check_fig5(log, quick, plot):
    log("=" * 110)
    log("5) J(u'_FE)/(k_h H^2 gamma_w^2) at h_w = H vs beta, Fig. 5 dashed curves (vector data), box 50/10/30 H, P2 ref 1")
    log("   far sides: zero_lb (u = 0 left + bottom, right no-flow; the default), zero_b, zero_l, impermeable;")
    log("   box variants with zero_lb: _d29 = depth 30 H below the crest (29 H below the toe), _rO = right side 10 H")
    log("   beyond O instead of beyond the toe.  rel = J_FE / J_paper - 1.  J_FE (conforming, exact piecewise-linear")
    log("   data) is an upper bound of the exact J of its box: the paper's values must lie ABOVE the converged J of")
    log("   the paper's box (P2 ref 1 here: converged to ~1e-5 relative, see check 3)")
    log("=" * 110)
    paper = _load_fig5_paper()
    betas = np.arange(15, 91, 15 if quick else 5)
    names = [n for n, _ in FIG5_VARIANTS]
    rows = []
    pmins = []
    for a in (1, 2, 4, 10):
        log(f"  alpha = {a}:  {'beta':>4} {'paper':>7} " + " ".join(f"{n:>17s}" for n in names))
        for b in betas:
            ref = paper[(a, int(b))][1] if paper else np.nan
            vals = []
            for _, kw in FIG5_VARIANTS:
                kw = dict(kw)
                if kw.get("right") == "from_O":
                    kw["right"] = PAPER_BOX["right"] - 1.0 / np.tan(np.radians(b))
                fe = FESeepage(b, 1.0, 1.0, a, order=2, ref=0 if quick else 1, **kw)
                vals.append(fe.J_normalized())
                pmin, _ = fe.min_pore_pressure()
                pmins.append(pmin)
            rows.append([a, b, ref] + vals)
            log(f"              {b:4.0f} {ref:7.4f} " + " ".join(f"{v:9.5f} {100 * (v / ref - 1):+6.2f}%" for v in vals))
    rows = np.array(rows)
    log(f"  min total pore pressure p = u + gw y over all dofs of all runs: {min(pmins):.3e} kPa (must be >= 0)")
    for k, n in enumerate(names):
        rel = rows[:, 3 + k] / rows[:, 2] - 1
        log(f"  {n:12s}: rel. diff to the paper dashed curves: mean {100 * rel.mean():+.3f}%, min {100 * rel.min():+.3f}%, "
            f"max {100 * rel.max():+.3f}%, RMS {100 * np.sqrt(np.mean(rel ** 2)):.3f}%; points with J_FE > J_paper "
            f"(impossible for the paper's box): {int(np.sum(rel > 0))} of {len(rel)}")
    out = os.path.join(DATA_DIR, "fe_seepage_fig5_J.csv")
    np.savetxt(out, rows, delimiter=",", fmt="%.6f",
               header="J(u'_FE)/(k_h H^2 gamma_w^2), h_w = H, box left 50H / right 10H beyond the toe / depth 30H below "
                      "the toe unless stated, P2 ref 1 quadtree mesh (fe_seepage.py); _d29: depth 29 H below the toe; "
                      "_rO: right side 10 H beyond O\nalpha,beta_deg,paper_fig5_dashed," + ",".join(f"J_{n}" for n in names))
    log(f"  -> written {out}")
    # the paper's own mesh (Fig. 4a: a few hundred 6-node triangles, about H/5 at the slope, several H far away)
    # gives an upper bound above the converged J: emulate it with coarse P2 meshes of the default generator
    log("  Effect of a coarse P2 mesh like Fig. 4a (h = hc at O, W, T and the face, grade 0.3, hmax 6-8 H) on J:")
    log(f"  {'alpha':>5} {'beta':>4} {'J fine':>9} {'paper':>8}  " +
        "  ".join(f"hc={hc:g} (nv)" for hc, _, _ in COARSE_MESHES))
    for a in (1, 10):
        for b in ((15, 90) if quick else (15, 45, 90)):
            fine = FESeepage(b, 1.0, 1.0, a, bc="zero_lb", order=2, ref=1).J_normalized()
            cs = []
            for hc, g, hm in COARSE_MESHES:
                co = FESeepage(b, 1.0, 1.0, a, bc="zero_lb", order=2, size=MeshSize(h0=hc, hs=hc, grade=g, hmax=hm))
                cs.append(f"{100 * (co.J_normalized() / fine - 1):+.3f}% ({co.mesh.nv})")
            ref = paper[(a, int(b))][1] if paper else np.nan
            log(f"  {a:5d} {b:4d} {fine:9.5f} {100 * (ref / fine - 1):+.3f}%  " + "  ".join(f"{c:>15s}" for c in cs))
    if plot:
        _plot_fig5(rows, names)
    return rows


def _plot_fig5(rows, names):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    paper = np.loadtxt(os.path.join(DATA_DIR, "fig5_vector_fill_polygons.csv"), delimiter=",", comments="#")
    ana = os.path.join(PROJECT_DIR, "results", "analytical_seepage", "fig5_analytical_Jstar.csv")
    ana = np.loadtxt(ana, delimiter=",", comments="#") if os.path.exists(ana) else None
    fig, axs = plt.subplots(2, 4, figsize=(17, 7.6), constrained_layout=True, gridspec_kw=dict(height_ratios=[1.5, 1]))
    col = {"zero_lb": "#2a78d6", "impermeable": "#eb6834", "zero_b": "#1baf7a", "zero_l": "#a35bd6",
           "zero_lb_d29": "#c9a227", "zero_lb_rO": "0.5"}
    lab = {"zero_lb": "FE, u = 0 left + bottom (default)", "impermeable": "FE, impermeable sides",
           "zero_b": "FE, u = 0 bottom only", "zero_l": "FE, u = 0 left only",
           "zero_lb_d29": "FE default BC, depth 30 H below the crest", "zero_lb_rO": "FE default BC, right side 10 H beyond O"}
    for k, (a, ymax) in enumerate(zip((1, 2, 4, 10), (1.0, 0.7, 0.5, 0.3))):
        ax, axr = axs[0, k], axs[1, k]
        p = paper[paper[:, 0] == a]
        ax.fill_between(p[:, 1], p[:, 2], p[:, 3], color="0.88", lw=0)
        ax.plot(p[:, 1], p[:, 3], "k--", lw=1.6, label="paper J(u'_FE) (dashed)")
        ax.plot(p[:, 1], p[:, 2], "k-", lw=1.2, label="paper -J*(v'_opt) (solid)")
        r = rows[rows[:, 0] == a]
        for i, n in enumerate(names):
            if n in ("zero_lb", "zero_l", "impermeable"):
                ax.plot(r[:, 1], r[:, 3 + i], "o", ms=4.5 if i == 0 else 3.5, mfc="none", mew=1.4 if i == 0 else 1.0,
                        color=col[n], label=lab[n])
            if n in ("zero_lb", "zero_b", "zero_lb_d29", "zero_l"):
                axr.plot(r[:, 1], 100 * (r[:, 3 + i] / r[:, 2] - 1), "-o", ms=3, lw=1.2 if i == 0 else 0.9,
                         color=col[n], label=lab[n])
        if ana is not None:
            q = ana[ana[:, 0] == a]
            ax.plot(q[:, 1], q[:, 2], "x", ms=4, color="0.35", label="analytical -J* (analytical_seepage.py)")
        ax.set_title(f"alpha = {a}", fontsize=10)
        ax.set_ylim(0, ymax)
        ax.set_ylabel("J / (k_h H^2 gamma_w^2)")
        axr.axhline(0.0, color="k", lw=0.8)
        axr.set_ylim(-1.3, 0.6)
        axr.set_ylabel("J_FE / J_paper - 1 (%)")
        axr.set_xlabel("slope inclination beta (deg)")
        for x in (ax, axr):
            x.set_xlim(15, 90)
            x.set_xticks(range(15, 91, 15))
            x.grid(color="0.92", lw=0.6)
    axs[0, 0].legend(fontsize=7.5, loc="lower right", frameon=False)
    axs[1, 0].legend(fontsize=7, loc="lower left", frameon=False)
    fig.suptitle("Fig. 5 check: FE hydraulic functional at h_w = H (P2, converged), box 50 H left of O / 10 H right of the "
                 "toe / 30 H below the toe; bottom: relative difference to the dashed curves (u = 0 left only: off scale "
                 "for small alpha; impermeable: -3.4 to -11 %, off scale)", fontsize=9.5)
    os.makedirs(RESULTS_DIR, exist_ok=True)
    fig.savefig(os.path.join(RESULTS_DIR, "fe_fig5_vs_paper.png"), dpi=120)
    plt.close(fig)


def check_fig4(log, plot, image=None):
    log("=" * 110)
    log("6) Fig. 4b: iso-lines of u'_FE (levels -0.98 k, gamma_w = 9.81, H = 1, h_w = H, beta = 45, alpha = 1);")
    log("   depth below the toe ground where each iso-line reaches the right side, s = depth / (30 H)")
    log("=" * 110)
    image = image or FIG4_IMAGE_DEFAULT
    nat = measure_fig4_native(PAPER_PDF_DEFAULT)
    if nat is not None:
        log(f"  native Fig. 4 image of the PDF (1000 x 1728 px): panel b box columns {nat['left_col']:.1f} - "
            f"{nat['right_col']:.1f}, crest row {nat['crest_row']:.1f}, toe-ground row {nat['toe_row']:.1f}, bottom row "
            f"{nat['bottom_row']:.1f}, crest edge O at column {nat['xO_col']:.0f}")
        for name in ("width61", "crest_to_toe"):
            e = nat[name]
            log(f"    H = {e['H_px']:.2f} px ({name}): left of O {e['left']:.2f} H, right of O {e['right_of_O']:.2f} H "
                f"(right of the toe {e['right_of_O'] - 1:.2f} H for beta = 45), depth below crest "
                f"{e['depth_below_crest']:.2f} H = below toe {e['depth_below_toe']:.2f} H")
    if os.path.exists(image):
        m = measure_fig4b(image)
        log(f"  re-measured on {image}: vertical box lines at columns {m['vertical_lines']}, crest rows "
            f"{m['crest_rows']}, toe-ground rows {m['right_top_rows']}")
        log(f"  right-side crossings (rows): {m['right_side_crossings']}")
    px = FIG4B_PX
    s_paper = np.array([(r - px["toe"]) / (px["bottom"] - px["toe"]) for _, r in FIG4B_RIGHT_ROWS])
    lev = np.array([v for v, _ in FIG4B_RIGHT_ROWS])
    sc = (px["right"] - px["left"]) / (50.0 + 10.0 + 1.0)
    log(f"  box measured in Fig. 4b (px): width {px['right'] - px['left']:.1f}, left of O {px['xO'] - px['left']:.1f}, "
        f"depth below toe {px['bottom'] - px['toe']:.1f}, crest-to-toe {px['toe'] - px['crest']:.1f}; with H = {sc:.2f} px "
        f"(box width = 61 H for beta = 45): left {(px['xO'] - px['left']) / sc:.1f} H, right of toe "
        f"{(px['right'] - px['xO']) / sc - 1:.1f} H, depth below toe {(px['bottom'] - px['toe']) / sc:.1f} H")
    sols = {bc: FESeepage(45.0, 1.0, 1.0, 1.0, gamma_w=9.81, bc=bc, order=2, ref=1) for bc in ("zero_lb", "impermeable")}
    yy = np.linspace(1.0, 31.0, 6001)
    log(f"  {'level':>7} {'paper s':>8} " + " ".join(f"{bc:>12s}" for bc in sols))
    s_fe = {}
    for bc, fe in sols.items():
        ur = fe.u(np.full_like(yy, fe.dom.xr - 1e-9), yy)
        s_fe[bc] = np.array([np.interp(v, ur, (yy - 1.0) / 30.0) if ur.min() <= v <= ur.max() else np.nan
                             for v in lev])
        s_fe[bc + "_curve"] = (yy, ur)
    for k, v in enumerate(lev):
        log(f"  {v:7.2f} {s_paper[k]:8.3f} " + " ".join(f"{s_fe[bc][k]:12.3f}" for bc in sols))
    for bc, fe in sols.items():
        ur = s_fe[bc + "_curve"][1]
        log(f"  {bc:11s}: u at the bottom-right corner = {ur[-1]:.3f} kPa (paper legend: 0.00); "
            f"RMS(s_FE - s_paper) = {np.sqrt(np.nanmean((s_fe[bc] - s_paper) ** 2)):.3f}")
    if plot:
        _plot_fig4(sols, s_paper, lev, s_fe, image)
    return s_paper, s_fe


def _plot_fig4(sols, s_paper, lev, s_fe, image):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.tri as mtri
    px = FIG4B_PX
    have_img = os.path.exists(image)
    fig = plt.figure(figsize=(13, 9.6), constrained_layout=True)
    gsp = fig.add_gridspec(2, 2)
    levels = np.sort(np.append(lev, -9.81 * 0.999))
    for k, (bc, fe) in enumerate(sols.items()):
        ax = fig.add_subplot(gsp[0, k])
        # map paper coordinates to the pixels of the page crop: box width = 61 H
        sc = (px["right"] - px["left"]) / (fe.dom.xr - fe.dom.xl)
        X = px["xO"] + fe.mesh.X[:, 0] * sc
        Y = px["crest"] + fe.mesh.X[:, 1] * sc
        if have_img:
            from PIL import Image
            im = np.array(Image.open(image).convert("L"))
            ax.imshow(im[560:1040, 0:900], cmap="gray", extent=(0, 900, 1040, 560), alpha=0.75)
        tri = mtri.Triangulation(X, Y, fe.mesh.T)
        ax.tricontour(tri, fe.uh[:fe.mesh.nv], levels=lev, colors="#2a78d6" if bc == "zero_lb" else "#eb6834",
                      linewidths=1.0)
        ax.set_xlim(0, 900)
        ax.set_ylim(1040, 560)
        ax.set_aspect("equal")
        ax.set_title(f"Fig. 4b (grey) + FE iso-lines -0.98k ({bc})", fontsize=10)
        ax.set_xticks([])
        ax.set_yticks([])
    ax = fig.add_subplot(gsp[1, 0])
    ax.plot(-lev, s_paper, "ks", ms=6, mfc="none", label="paper Fig. 4b (iso-line ends, measured)")
    for bc, c in (("zero_lb", "#2a78d6"), ("impermeable", "#eb6834")):
        yy, ur = s_fe[bc + "_curve"]
        ax.plot(-ur, (yy - 1.0) / 30.0, color=c, lw=2, label=f"FE right side, {bc}")
    ax.invert_yaxis()
    ax.set_xlabel("-u on the right side (kPa)")
    ax.set_ylabel("depth below toe ground / 30 H")
    ax.grid(color="0.92", lw=0.6)
    ax.legend(fontsize=8, frameon=False)
    ax.set_title("right side of the box: where the iso-lines end", fontsize=10)
    ax = fig.add_subplot(gsp[1, 1])
    fe = sols["zero_lb"]
    gx, gy = np.meshgrid(np.linspace(-1.0, 2.2, 23), np.linspace(0.04, 2.0, 15))
    keep = gy > fe.dom.ground(gx) + 0.02
    fx, fy = fe.force(gx[keep], gy[keep])
    ax.plot([-1.2, 0, fe.dom.xT, 2.4], [0, 0, 1, 1], "k-", lw=1.2)
    # angles="xy": arrows drawn in data coordinates, so with the inverted (downward) y axis fy > 0 points down
    ax.quiver(gx[keep], gy[keep], fx, fy, angles="xy", scale_units="xy", scale=60.0, width=0.0025, color="#2a78d6")
    ax.invert_yaxis()
    ax.set_aspect("equal")
    ax.set_title("seepage force -grad u'_FE near the slope (cf. Fig. 4c), zero_lb, P2 ref 1", fontsize=10)
    ax.set_xlabel("x / H")
    ax.set_ylabel("y / H (down)")
    fig.savefig(os.path.join(DATA_DIR, "fe_seepage_fig4_check.png"), dpi=110)
    plt.close(fig)


def check_field(log):
    log("=" * 110)
    log("7) Field interface: force(x, y) = -grad u'_FE, vectorised point location, P1 vs P2, outside points")
    log("=" * 110)
    f1 = FESeepage(60.0, 5.0, 5.0, 5.0, gamma_w=9.81, bc="zero_lb", order=1, ref=2)
    f2 = FESeepage(60.0, 5.0, 5.0, 5.0, gamma_w=9.81, bc="zero_lb", order=2, ref=2)
    rng = np.random.default_rng(3)
    x = rng.uniform(-15, 15, 200000)
    y = rng.uniform(-2, 15, 200000)
    t = time.time()
    fx, fy = f2.force(x, y)
    dt = time.time() - t
    ins = f2.dom.in_box(x, y)
    log(f"  H = 5 m, beta = 60, alpha = 5: P2 ref 2 force at 2e5 random points in {dt:.2f} s; "
        f"points outside the soil returning 0: {np.all((fx[~ins] == 0) & (fy[~ins] == 0))}")
    xs, ys = _soil_samples(f2.dom, 15.0, 0.25)
    a1 = np.array(f1.force(xs, ys))
    a2 = np.array(f2.force(xs, ys))
    rel = np.sqrt(np.mean(np.sum((a1 - a2) ** 2, 0)) / np.mean(np.sum(a2 ** 2, 0)))
    log(f"  P1 ref 2 vs P2 ref 2 within 3H of O: RMS rel. diff of f = {100 * rel:.2f}%;  J/(kh H^2 gw^2): "
        f"P1 {f1.J_normalized():.6f}, P2 {f2.J_normalized():.6f}")
    f2c = FESeepage(60.0, 5.0, 5.0, 5.0, gamma_w=9.81, bc="zero_lb", order=2, ref=1)
    a2c = np.array(f2c.force(xs, ys))
    rel_c = np.sqrt(np.mean(np.sum((a2c - a2) ** 2, 0)) / np.mean(np.sum(a2 ** 2, 0)))
    far = np.hypot(xs - f2.dom.xT, ys - f2.dom.H) > 0.25 * f2.H        # away from the singular toe
    rel_cf = np.sqrt(np.mean(np.sum((a2c - a2)[:, far] ** 2, 0)) / np.mean(np.sum(a2[:, far] ** 2, 0)))
    log(f"  P2 ref 1 vs P2 ref 2 within 3H of O: RMS rel. diff of f = {100 * rel_c:.3f}% "
        f"({100 * rel_cf:.3f}% beyond 0.25 H from the toe)")
    shape_ok = f2.force(np.zeros((3, 4)), np.ones((3, 4)))[0].shape == (3, 4)
    log(f"  arbitrary array shapes: {shape_ok};  u(-H/2, H) = {float(f2.u(-2.5, 5.0)):.4f} kPa, f = "
        f"{tuple(float(v) for v in f2.force(-2.5, 5.0))}")
    # first call on a fresh default field (P2 ref 1) includes the construction of the point locator
    f3 = FESeepage(45.0, 5.0, 2.0, 4.0, gamma_w=9.81)
    x3, y3 = rng.uniform(-15, 15, 100000), rng.uniform(-2, 15, 100000)
    t = time.time()
    f3.force(x3, y3)
    log(f"  default field (P2 ref 1), 1e5 random points: first call (incl. point-locator build) {time.time() - t:.2f} s")
    # points ON the ground surface (crest, face, toe ground) belong to the closed soil domain: f must not be
    # 0 there and u must equal the Dirichlet data (before the _locate fix, 1890 of 2000 face points gave 0 / nan)
    for fe in (f2, f3):
        d = fe.dom
        s = rng.uniform(0.0, 1.0, 2000)
        xs_ = np.concatenate([rng.uniform(-2.0, 0.0, 1000) * fe.H, s * d.xT, d.xT + rng.uniform(0.0, 2.0, 1000) * fe.H])
        ys_ = np.concatenate([np.zeros(1000), s * d.H, np.full(1000, d.H)])
        fx_, fy_ = fe.force(xs_, ys_)
        uu = fe.u(xs_, ys_)
        err = np.nanmax(np.abs(uu - d.surface_u(xs_, ys_, fe.gamma_w)))
        log(f"  points on the ground surface (beta = {d.beta_deg:g}, h_w/H = {fe.hw / fe.H:g}): {len(xs_)}; with f = 0: "
            f"{int(np.sum((fx_ == 0) & (fy_ == 0)))}; u = nan: {int(np.isnan(uu).sum())}; max |u - Eq. 21 data| = {err:.1e} kPa")
    # exact similarity in H: the mesh is generated for H = 1 and scaled
    f4 = FESeepage(45.0, 1.0, 0.4, 4.0, gamma_w=9.81)
    xs4, ys4 = _soil_samples(f4.dom, 3.0, 0.05)
    d4 = np.array(f3.force(5.0 * xs4, 5.0 * ys4)) - np.array(f4.force(xs4, ys4))
    log(f"  similarity: |f(H = 5)(5x, 5y) - f(H = 1)(x, y)| <= {np.abs(d4).max():.1e} kN/m^3 at {len(xs4)} points; "
        f"J/(kh H^2 gw^2) {f3.J_normalized():.12f} vs {f4.J_normalized():.12f}")
    return rel


def write_reference_points(log):
    """u and f at probe points for the C++ (NeoPZ) regression, P2 ref 2"""
    rows = []
    pts = [(-2.0, 0.5), (-0.5, 0.25), (-0.5, 1.0), (0.2, 0.5), (0.5, 1.2), (1.0, 2.0), (2.0, 1.5), (-1.0, 3.0),
           (0.0, 5.0), (4.0, 3.0)]
    cases = [(45, 1, 1.0), (90, 1, 1.0), (30, 10, 1.0), (60, 5, 1.0), (60, 4, 0.5)]
    for ib, bc in enumerate(("zero_lb", "impermeable")):
        for b, a, hwr in cases:
            fe = FESeepage(b, 1.0, hwr, a, gamma_w=9.81, bc=bc, order=2, ref=2)
            for x, y in pts:
                if not fe.dom.in_box(x, y) or y <= fe.dom.ground(x) + 0.05:
                    continue
                u = float(fe.u(x, y))
                fx, fy = (float(v) for v in fe.force(x, y))
                rows.append((ib, b, a, hwr, x, y, u, fx, fy, fe.J_normalized()))
    out = os.path.join(DATA_DIR, "fe_seepage_reference_points.csv")
    np.savetxt(out, np.array(rows), delimiter=",", fmt="%.8g",
               header="FE reference (fe_seepage.py, P2 ref 2 quadtree mesh, box 50/10/30 H; bc 0 = zero_lb: u = 0 on "
                      "left + bottom, right no-flow; bc 1 = impermeable), H = 1, gamma_w = 9.81, k_h = 1, paper "
                      "coordinates (y down)\nbc,beta_deg,alpha,hw_over_H,x,y,u_kPa,fx,fy,J_over_kh_H2_gw2")
    log(f"  -> written {out} ({len(rows)} rows)")


ALL_CHECKS = "mesh,mms,cross,conv,sens,fig5,fig4,field,ref"


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--no-plot", action="store_true")
    ap.add_argument("--only", default=ALL_CHECKS, help="comma list of " + ALL_CHECKS)
    ap.add_argument("--fig4-image", default=None, help="page crop of Fig. 4 (PDF page 8) for the overlay")
    args = ap.parse_args(argv)
    os.makedirs(RESULTS_DIR, exist_ok=True)
    only = set(args.only.split(","))
    lines = []

    def log(s=""):
        print(s, flush=True)
        lines.append(s)

    t0 = time.time()
    if "mesh" in only:
        check_mesh(log, args.quick, not args.no_plot)
    if "mms" in only:
        check_mms(log, args.quick)
    if "cross" in only:
        check_laplace_cross(log)
    if "conv" in only:
        check_convergence(log, args.quick)
    if "sens" in only:
        check_sensitivity(log, args.quick)
    if "fig5" in only:
        check_fig5(log, args.quick, not args.no_plot)
    if "fig4" in only:
        check_fig4(log, not args.no_plot, args.fig4_image)
    if "field" in only:
        check_field(log)
    if "ref" in only:
        write_reference_points(log)
    log(f"total time {time.time() - t0:.1f} s")
    tag = "_quick" if args.quick else ""
    if only == set(ALL_CHECKS.split(",")):
        with open(os.path.join(RESULTS_DIR, f"checks_output{tag}.txt"), "w") as f:
            f.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
