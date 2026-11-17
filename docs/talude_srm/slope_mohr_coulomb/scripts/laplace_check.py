"""Independent check (numpy only, no NeoPZ) of the hydraulic part of Projects/SlopeDrawdown.

Mesh: TriGMesh(ref) of Projects/SlopeMohrCoulomb/SlopeModel.h re-implemented (36 nodes, 49 triangles,
21 boundary lines, uniform refinement by edge midpoints, boundary ids kept).
Problem: total head h = y + p/gw, div grad h = 0, h = y on -3, -4, -6, no flow on -1, -2, -5 (steady,
reservoir drawn down to y = 30, surface drained); crest state h = 40 (p = gw (40 - y)).
Linear (P1) finite elements, own assembly; sparse CG (Jacobi) for the convergence study.
"""
import sys, os, csv
import numpy as np

GW = 10.0
HERE = os.path.dirname(os.path.abspath(__file__))
# usage: python laplace_check.py <directory of a SlopeDrawdown run> [<output directory>]
DD3 = sys.argv[1] if len(sys.argv) > 1 else "."
OUT = sys.argv[2] if len(sys.argv) > 2 else "."

CO = [(0, 0), (10, 0), (20, 0), (30, 0), (40, 0), (50, 0), (60, 0), (70, 0),
      (0, 10), (10, 10), (20, 10), (30, 10), (40, 10), (50, 10), (60, 10), (70, 10),
      (0, 20), (10, 20), (20, 20), (30, 20), (40, 20), (50, 20), (60, 20), (70, 20),
      (0, 30), (10, 30), (20, 30), (30, 30), (40, 30), (50, 30), (60, 30), (70, 30),
      (0, 40), (10, 40), (20, 40), (30, 40)]
TRI = [(0, 1, 8), (1, 9, 8), (1, 2, 9), (2, 10, 9), (2, 3, 10), (3, 11, 10), (3, 4, 11),
       (4, 12, 11), (4, 5, 12), (5, 13, 12), (5, 6, 13), (6, 14, 13), (6, 7, 14), (7, 15, 14),
       (8, 9, 16), (9, 17, 16), (9, 10, 17), (10, 18, 17), (10, 11, 18), (11, 19, 18), (11, 12, 19),
       (12, 20, 19), (12, 13, 20), (13, 21, 20), (13, 14, 21), (14, 22, 21), (14, 15, 22), (15, 23, 22),
       (16, 17, 24), (17, 25, 24), (17, 18, 25), (18, 26, 25), (18, 19, 26), (19, 27, 26), (19, 20, 27),
       (20, 28, 27), (20, 21, 28), (21, 29, 28), (21, 22, 29), (22, 30, 29), (22, 23, 30), (23, 31, 30),
       (24, 25, 32), (25, 33, 32), (25, 26, 33), (26, 34, 33), (26, 27, 34), (27, 35, 34), (27, 28, 35)]
LINES = [((0, 1), -1), ((1, 2), -1), ((2, 3), -1), ((3, 4), -1), ((4, 5), -1), ((5, 6), -1), ((6, 7), -1),
         ((7, 15), -2), ((15, 23), -2), ((23, 31), -2),
         ((31, 30), -3), ((30, 29), -3), ((29, 28), -3),
         ((35, 34), -4), ((34, 33), -4), ((33, 32), -4),
         ((32, 24), -5), ((24, 16), -5), ((16, 8), -5), ((8, 0), -5),
         ((28, 35), -6)]


def tri_gmesh(ref):
    X = [list(map(float, c)) for c in CO]
    tri = [tuple(t) for t in TRI]
    lines = [(l, i) for l, i in LINES]
    for _ in range(ref):
        mid = {}

        def m(a, b):
            k = (min(a, b), max(a, b))
            if k not in mid:
                X.append([(X[a][0] + X[b][0]) / 2, (X[a][1] + X[b][1]) / 2])
                mid[k] = len(X) - 1
            return mid[k]

        nt = []
        for a, b, c in tri:
            ab, bc, ca = m(a, b), m(b, c), m(c, a)
            nt += [(a, ab, ca), (ab, b, bc), (ca, bc, c), (bc, ca, ab)]
        nl = []
        for (a, b), i in lines:
            mm = m(a, b)
            nl += [((a, mm), i), ((mm, b), i)]
        tri, lines = nt, nl
    return np.array(X), np.array(tri, dtype=int), lines


def p1_stiffness(X, T):
    """element gradients and areas; returns COO (rows, cols, vals) of the P1 Laplace matrix"""
    x = X[T]  # (ne, 3, 2)
    d1 = x[:, 1] - x[:, 0]
    d2 = x[:, 2] - x[:, 0]
    det = d1[:, 0] * d2[:, 1] - d1[:, 1] * d2[:, 0]
    assert np.all(det > 0), "triangles must be counter-clockwise"
    area = det / 2
    # gradients of barycentric functions: grad L_i = rot90(edge opposite i) / (2A)
    e = np.stack([x[:, 2] - x[:, 1], x[:, 0] - x[:, 2], x[:, 1] - x[:, 0]], axis=1)  # (ne,3,2)
    G = np.stack([-e[:, :, 1], e[:, :, 0]], axis=2) / det[:, None, None]  # (ne,3,2)
    K = area[:, None, None] * np.einsum("eik,ejk->eij", G, G)
    rows = np.repeat(T, 3, axis=1).ravel()
    cols = np.tile(T, (1, 3)).ravel()
    return rows, cols, K.ravel(), G, area


class Sparse:
    def __init__(self, n, r, c, v):
        key = r.astype(np.int64) * n + c
        u, inv = np.unique(key, return_inverse=True)
        self.v = np.bincount(inv, weights=v)
        self.r = (u // n).astype(int)
        self.c = (u % n).astype(int)
        self.n = n

    def mv(self, x):
        return np.bincount(self.r, weights=self.v * x[self.c], minlength=self.n)

    def diag(self):
        d = np.zeros(self.n)
        m = self.r == self.c
        d[self.r[m]] = self.v[m]
        return d

    def dense(self):
        A = np.zeros((self.n, self.n))
        A[self.r, self.c] = self.v
        return A


def solve_dirichlet(X, T, lines, dir_ids, gfun, dense_max=3000):
    n = len(X)
    r, c, v, G, area = p1_stiffness(X, T)
    A = Sparse(n, r, c, v)
    dnodes = sorted({k for (a, b), i in lines if i in dir_ids for k in (a, b)})
    dmask = np.zeros(n, bool)
    dmask[dnodes] = True
    u = np.zeros(n)
    u[dmask] = gfun(X[dmask])
    rhs = -A.mv(u)
    free = ~dmask
    if n <= dense_max:
        Ad = A.dense()
        u[free] = np.linalg.solve(Ad[np.ix_(free, free)], rhs[free])
    else:  # Jacobi-preconditioned CG on the free unknowns
        dg = A.diag()
        x = np.zeros(n)
        b = np.where(free, rhs, 0.)

        def Af(y):
            y = np.where(free, y, 0.)
            return np.where(free, A.mv(y), 0.)

        rr = b - Af(x)
        z = np.where(free, rr / dg, 0.)
        p = z.copy()
        rz = rr @ z
        nb = np.linalg.norm(b)
        for it in range(20 * n):
            Ap = Af(p)
            al = rz / (p @ Ap)
            x += al * p
            rr -= al * Ap
            if np.linalg.norm(rr) < 1e-13 * nb:
                break
            z = np.where(free, rr / dg, 0.)
            rz2 = rr @ z
            p = z + rz2 / rz * p
            rz = rz2
        u[free] = x[free]
    # residual = reaction (boundary flux) at Dirichlet nodes: R = A u  (outflow of k grad h integrated, k = 1)
    R = A.mv(u)
    return u, dmask, R, A


def read_vtk(path):
    with open(path) as f:
        L = f.read().split("\n")
    i = next(k for k, s in enumerate(L) if s.startswith("POINTS"))
    npnt = int(L[i].split()[1])
    P = np.array([list(map(float, L[i + 1 + k].split())) for k in range(npnt)])[:, :2]
    j = next(k for k, s in enumerate(L) if s.startswith("CELLS"))
    ncell = int(L[j].split()[1])
    C = np.array([list(map(int, L[j + 1 + k].split()))[1:] for k in range(ncell)])
    s = next(k for k, s in enumerate(L) if s.startswith("SCALARS PorePressure"))
    vals = []
    k = s + 2
    while len(vals) < npnt:
        vals += [float(t) for t in L[k].split()]
        k += 1
    return P, C, np.array(vals[:npnt])


def match_vertices(X, P, pv, tol=1e-6):
    """for every mesh vertex: all VTK points within tol -> mean value, spread, count"""
    key = lambda a: (np.round(a[:, 0] / tol).astype(np.int64), np.round(a[:, 1] / tol).astype(np.int64))
    kx, ky = key(P)
    d = {}
    for idx, (a, b) in enumerate(zip(kx, ky)):
        d.setdefault((a, b), []).append(idx)
    vx, vy = key(X)
    val = np.full(len(X), np.nan)
    spread = np.zeros(len(X))
    cnt = np.zeros(len(X), int)
    for i, (a, b) in enumerate(zip(vx, vy)):
        ids = d.get((a, b))
        if ids is None:  # fallback: brute force within tol
            dist = np.hypot(P[:, 0] - X[i, 0], P[:, 1] - X[i, 1])
            ids = list(np.nonzero(dist < tol)[0])
        if ids:
            vals = pv[ids]
            val[i] = vals.mean()
            spread[i] = vals.max() - vals.min()
            cnt[i] = len(ids)
    return val, spread, cnt


def interp_p1(X, T, u, Q):
    """evaluate the P1 field u at the points Q (brute force barycentric search)"""
    x = X[T]
    out = np.full(len(Q), np.nan)
    a = x[:, 0]
    d1 = x[:, 1] - a
    d2 = x[:, 2] - a
    det = d1[:, 0] * d2[:, 1] - d1[:, 1] * d2[:, 0]
    for k, q in enumerate(Q):
        dq = q - a
        l1 = (dq[:, 0] * d2[:, 1] - dq[:, 1] * d2[:, 0]) / det
        l2 = (d1[:, 0] * dq[:, 1] - d1[:, 1] * dq[:, 0]) / det
        l0 = 1 - l1 - l2
        ok = np.nonzero((l0 > -1e-9) & (l1 > -1e-9) & (l2 > -1e-9))[0]
        if len(ok):
            e = ok[0]
            out[k] = l0[e] * u[T[e, 0]] + l1[e] * u[T[e, 1]] + l2[e] * u[T[e, 2]]
    return out


def write_csv(path, header, cols):
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        for row in zip(*cols):
            w.writerow([f"{v:.10g}" if isinstance(v, (float, np.floating)) else v for v in row])


def main():
    out = OUT
    res = {}
    # ---------------- mesh ----------------
    X0, T0, L0 = tri_gmesh(0)
    X, T, Lb = tri_gmesh(2)
    ids = sorted({i for _, i in Lb})
    nlin = {i: sum(1 for _, j in Lb if j == i) for i in ids}
    print(f"TriGMesh(0): {len(X0)} nodes, {len(T0)} triangles, {len(L0)} lines")
    print(f"TriGMesh(2): {len(X)} nodes, {len(T)} triangles, {len(Lb)} lines {nlin}")
    # sanity: total area 70*30 + 30*10 + 0.5*10*10 = 2400 + ...
    _, _, _, _, area = p1_stiffness(X, T)
    print(f"area = {area.sum():.6f} (expected 70*30 + 30*10 + 10*10/2 = {70*30 + 30*10 + 50})")
    # boundary length check
    blen = {i: sum(np.hypot(*(X[b] - X[a])) for (a, b), j in Lb if j == i) for i in ids}
    print("boundary lengths:", {i: round(v, 6) for i, v in blen.items()})

    # ---------------- patch test: h = 3 + 2x - 5y with Dirichlet on all boundaries ----------------
    allids = {-1, -2, -3, -4, -5, -6}
    lin = lambda Y: 3 + 2 * Y[:, 0] - 5 * Y[:, 1]
    up, _, _, _ = solve_dirichlet(X, T, Lb, allids, lin)
    print(f"patch test (linear h, all Dirichlet): max err = {np.abs(up - lin(X)).max():.3e}")
    # patch test with natural BC: h = y is exact for Dirichlet on -3,-4,-6 only if flux zero on -1,-2,-5 -> not;
    # use h = const 40 (crest state) which is exact
    hc, dm, Rc, _ = solve_dirichlet(X, T, Lb, {-3, -4, -6}, lambda Y: np.full(len(Y), 40.))
    pc = GW * (hc - X[:, 1])
    pc_exact = GW * (40. - X[:, 1])
    print(f"crest: max |p_py - gw(40-y)| = {np.abs(pc - pc_exact).max():.3e}")

    # ---------------- steady drawn-down state ----------------
    h, dm, R, A = solve_dirichlet(X, T, Lb, {-3, -4, -6}, lambda Y: Y[:, 1].copy())
    p = GW * (h - X[:, 1])
    # boundary discharge: reactions at Dirichlet nodes (k = 1, units m^2/unit k per m)
    Q_in = R[dm][R[dm] > 0].sum()
    Q_out = -R[dm][R[dm] < 0].sum()
    print(f"steady: h in [{h.min():.4f}, {h.max():.4f}], p in [{p.min():.4f}, {p.max():.4f}] kPa, "
          f"sum of reactions = {R[dm].sum():.3e}, |interior residual| = {np.abs(R[~dm]).max():.3e}")
    # reactions by boundary group
    def node_group(i):
        g = set(j for (a, b), j in Lb if i in (a, b) and j in (-3, -4, -6))
        return g
    Rg = {-3: 0., -4: 0., -6: 0.}
    for i in np.nonzero(dm)[0]:
        g = sorted(node_group(i))
        for j in g:
            Rg[j] += R[i] / len(g)
    print(f"steady discharge (k=1): R>0 sum {Q_in:.6f}, R<0 sum {-Q_out:.6f}; by group (corner nodes split) {Rg}")
    res["Q"] = Q_in

    # ---------------- comparison with NeoPZ ----------------
    for state, fn, pyv in [("steady", "drawdown_up.scal_vec.5.vtk", p), ("crest", "drawdown_up.scal_vec.0.vtk", pc),
                           ("T10", "drawdown_up.scal_vec.4.vtk", p)]:
        P, C, pv = read_vtk(os.path.join(DD3, fn))
        val, spread, cnt = match_vertices(X, P, pv)
        nm = np.isfinite(val).sum()
        dif = val - pyv
        print(f"[{state}] {fn}: {len(P)} VTK points, {len(C)} cells; matched {nm}/{len(X)} vertices "
              f"(VTK copies per vertex {cnt.min()}..{cnt.max()}, max spread among copies {spread.max():.3e})")
        print(f"[{state}] vertices: max|p_neopz - p_py| = {np.nanmax(np.abs(dif)):.6e}, "
              f"RMS = {np.sqrt(np.nanmean(dif ** 2)):.6e}, max|p_py| = {np.abs(pyv).max():.6f}, "
              f"max|p_neopz| = {np.nanmax(np.abs(val)):.6f}, argmax at {X[np.nanargmax(np.abs(dif))]}")
        # all VTK points (also edge midpoints / interior points of the graph mesh): P1 interpolation of p_py
        pi = interp_p1(X, T, pyv, P)
        d2 = pv - pi
        print(f"[{state}] all {len(P)} VTK points (P1 interp of p_py): max|diff| = {np.nanmax(np.abs(d2)):.6e}, "
              f"RMS = {np.sqrt(np.nanmean(d2 ** 2)):.6e}, nan = {np.isnan(pi).sum()}")
        if state == "crest":
            de = val - pc_exact
            print(f"[crest] NeoPZ vs exact gw(40-y): max = {np.nanmax(np.abs(de)):.6e}, RMS = {np.sqrt(np.nanmean(de**2)):.6e}")
        res[state] = (val, dif)
        if state != "T10":
            name = os.path.join(out, f"vertices_{state}.csv")
            hcol = h if state == "steady" else hc
            cols = [X[:, 0], X[:, 1], hcol, pyv, val, dif]
            hdr = ["x", "y", "h_python", "p_python", "p_neopz", "p_neopz_minus_p_python"]
            if state == "crest":
                cols.append(pc_exact)
                hdr.append("p_exact")
            write_csv(name, hdr, cols)
            name2 = os.path.join(out, f"vtkpoints_{state}.csv")
            write_csv(name2, ["x", "y", "p_neopz", "p_python_interp"], [P[:, 0], P[:, 1], pv, pi])
    # T10 vs steady from NeoPZ
    d = res["T10"][0] - res["steady"][0]
    print(f"NeoPZ T10 - steady at vertices: max|.| = {np.nanmax(np.abs(d)):.6e}")

    # head field + mesh
    write_csv(os.path.join(out, "head_steady.csv"), ["x", "y", "h_python", "p_python", "dirichlet"],
              [X[:, 0], X[:, 1], h, p, dm.astype(int)])
    write_csv(os.path.join(out, "triangles_ref2.csv"), ["n0", "n1", "n2"], [T[:, 0], T[:, 1], T[:, 2]])
    write_csv(os.path.join(out, "boundary_ref2.csv"), ["n0", "n1", "bc_id"],
              [[a for (a, b), i in Lb], [b for (a, b), i in Lb], [i for _, i in Lb]])

    # ---------------- mesh convergence of the P1 solution (coarse vertices) ----------------
    conv = []
    hs = {}
    for r in range(0, 6):
        Xr, Tr, Lr = tri_gmesh(r)
        hr, dmr, Rr, _ = solve_dirichlet(Xr, Tr, Lr, {-3, -4, -6}, lambda Y: Y[:, 1].copy())
        Qr = Rr[dmr][Rr[dmr] > 0].sum()
        hs[r] = (Xr, Tr, hr)
        conv.append((r, len(Xr), len(Tr), hr[:36].copy(), Qr))
        print(f"ref {r}: {len(Xr)} nodes, {len(Tr)} tri, max p = {GW * (hr - Xr[:, 1]).max():.6f}, "
              f"p(0,0) = {GW * (hr[0] - 0):.6f}, p(70,0) = {GW * hr[7]:.6f}, Q(k=1) = {Qr:.6f}")
    href = conv[-1][3]
    rows = []
    for r, nn, nt, h36, Qr in conv:
        e = GW * np.abs(h36 - href).max()
        rows.append((r, nn, nt, e, GW * h36[0], GW * h36[7], Qr))
        print(f"ref {r}: max |p_r - p_ref5| at the 36 coarse nodes = {e:.6e}")
    # also: ref2 P1 field vs ref5 at all ref2 vertices
    X5, T5, h5 = hs[5]
    X2 = hs[2][0]
    # vertices of ref2 are the first len(X2) vertices of ref5 (nested numbering)
    assert np.allclose(X5[:len(X2)], X2)
    e2 = GW * np.abs(hs[2][2] - h5[:len(X2)])
    print(f"ref2 vs ref5 at all {len(X2)} ref2 vertices: max |dp| = {e2.max():.6e}, RMS = {np.sqrt((e2**2).mean()):.6e}")
    write_csv(os.path.join(out, "convergence.csv"),
              ["ref", "nodes", "triangles", "max_abs_dp_coarse_nodes_vs_ref5", "p_00", "p_700", "Q_k1"],
              list(zip(*rows)))
    np.savez(os.path.join(out, "fields.npz"), X=X, T=T, h=h, p=p, pc=pc,
             p_neopz_steady=res["steady"][0], p_neopz_crest=res["crest"][0])


if __name__ == "__main__":
    main()
