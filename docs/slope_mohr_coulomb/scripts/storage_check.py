"""Explain p_NeoPZ - p_Laplace by the storage term of the last backward-Euler step of the u-p analysis.

Mass balance of TPZMatPoroElastoPlasticUP (Rp): int Np alpha (tr eps - tr eps_n) + dt k int grad Np . (grad p - rho_w g) = 0
=> p_{n+1} = p_Laplace + dp,  k K dp = -(1/dt) int Np alpha (div u_{n+1} - div u_n)   (dp = 0 on -3,-4,-6)
u is P2 (Taylor-Hood); its vertex and edge-midpoint values are read from the NeoPZ VTK (resolution 1).
"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from laplace_check import tri_gmesh, p1_stiffness, Sparse, read_vtk, match_vertices, write_csv, DD3, GW

E, NU, KH = 20000., 0.3, 1.e-2
MOB = KH / GW
M = E * (1 - NU) / ((1 + NU) * (1 - 2 * NU))
TC = 100. / (MOB * M)


def read_disp(path):
    with open(path) as f:
        L = f.read().split("\n")
    i = next(k for k, s in enumerate(L) if s.startswith("POINTS"))
    n = int(L[i].split()[1])
    P = np.array([list(map(float, L[i + 1 + k].split())) for k in range(n)])[:, :2]
    j = next(k for k, s in enumerate(L) if s.startswith("VECTORS Displacement"))
    U = np.array([list(map(float, L[j + 1 + k].split())) for k in range(n)])[:, :2]
    return P, U


def nodal_lookup(P, V, Q, tol=1e-6):
    d = {}
    for idx, (a, b) in enumerate(np.round(P / tol).astype(np.int64)):
        d.setdefault((a, b), []).append(idx)
    out = np.zeros((len(Q), V.shape[1]))
    spread = 0.
    for k, (a, b) in enumerate(np.round(Q / tol).astype(np.int64)):
        ids = d[(a, b)]
        out[k] = V[ids].mean(0)
        spread = max(spread, np.abs(V[ids] - V[ids].mean(0)).max())
    return out, spread


def storage_load(X, T, Uv, Um):
    """b_i = int div(u) phi_i (P2 u, P1 phi), exact (3-point edge-midpoint rule); Uv: vertex values (n,2),
    Um: per element midpoint values (ne,3,2) in order m01, m12, m20"""
    n = len(X)
    b = np.zeros(n)
    x = X[T]
    d1 = x[:, 1] - x[:, 0]
    d2 = x[:, 2] - x[:, 0]
    det = d1[:, 0] * d2[:, 1] - d1[:, 1] * d2[:, 0]
    area = det / 2
    e = np.stack([x[:, 2] - x[:, 1], x[:, 0] - x[:, 2], x[:, 1] - x[:, 0]], axis=1)
    G = np.stack([-e[:, :, 1], e[:, :, 0]], axis=2) / det[:, None, None]  # grad L_i (ne,3,2)
    uv = Uv[T]  # (ne,3,2)
    pairs = [(0, 1), (1, 2), (2, 0)]
    divs = []
    for q, (a, c) in enumerate(pairs):  # quadrature points: edge midpoints, L_a = L_c = 1/2
        L = np.zeros(3)
        L[a] = L[c] = 0.5
        div = np.zeros(len(T))
        for i in range(3):  # vertex functions: grad N_i = (4 L_i - 1) grad L_i
            div += (4 * L[i] - 1) * np.einsum("ek,ek->e", uv[:, i], G[:, i])
        for m, (i, j) in enumerate(pairs):  # midpoint functions: grad N_ij = 4 (L_i grad L_j + L_j grad L_i)
            g = 4 * (L[i] * G[:, j] + L[j] * G[:, i])
            div += np.einsum("ek,ek->e", Um[:, m], g)
        divs.append(div)
        for i in range(3):
            np.add.at(b, T[:, i], area / 3 * div * L[i])
    return b, np.array(divs).T


def main():
    X, T, Lb = tri_gmesh(2)
    mids = np.stack([(X[T[:, 0]] + X[T[:, 1]]) / 2, (X[T[:, 1]] + X[T[:, 2]]) / 2, (X[T[:, 2]] + X[T[:, 0]]) / 2], 1)
    r, c, v, _, _ = p1_stiffness(X, T)
    A = Sparse(len(X), r, c, v).dense()
    dnodes = sorted({k for (a, b), i in Lb if i in (-3, -4, -6) for k in (a, b)})
    free = np.ones(len(X), bool)
    free[dnodes] = False
    d = np.load(os.path.join(os.path.dirname(os.path.abspath(__file__)), "fields.npz"))
    print(f"M = {M:.4f} kPa, mobility = {MOB}, Tc = {TC:.6f} day")
    cases = [("crest", None, "drawdown_up.scal_vec.0.vtk", 1.e4 * TC, d["pc"]),
             ("steady", "drawdown_up.scal_vec.4.vtk", "drawdown_up.scal_vec.5.vtk", 1.e4 * TC - 10. * TC, d["p"])]
    for name, f0, f1, dt, plap in cases:
        P1, U1 = read_disp(os.path.join(DD3, f1))
        if f0 is None:
            U0 = np.zeros_like(U1)
        else:
            P0, U0 = read_disp(os.path.join(DD3, f0))
            assert np.allclose(P0, P1)
        dU = U1 - U0
        Uv, s1 = nodal_lookup(P1, dU, X)
        Um, s2 = nodal_lookup(P1, dU, mids.reshape(-1, 2))
        Um = Um.reshape(-1, 3, 2)
        b, divs = storage_load(X, T, Uv, Um)
        dp = np.zeros(len(X))
        rhs = -b / (dt * MOB)
        dp[free] = np.linalg.solve(A[np.ix_(free, free)], rhs[free])
        _, _, pz = read_vtk(os.path.join(DD3, f1))
        pv, _, _ = match_vertices(X, P1, pz)
        pred = plap + dp
        e0 = pv - plap
        e1 = pv - pred
        print(f"[{name}] dt = {dt:.4f} day; delta div u in [{divs.min():.4e}, {divs.max():.4e}] "
              f"(continuity spread of u copies {max(s1, s2):.1e})")
        print(f"[{name}] predicted storage correction dp in [{dp.min():.4f}, {dp.max():.4f}] kPa")
        print(f"[{name}] |p_neopz - p_laplace|: max {np.abs(e0).max():.4e}, RMS {np.sqrt((e0**2).mean()):.4e}")
        print(f"[{name}] |p_neopz - (p_laplace + dp)|: max {np.abs(e1).max():.4e}, RMS {np.sqrt((e1**2).mean()):.4e}")
        write_csv(os.path.join(os.path.dirname(os.path.abspath(__file__)), f"storage_{name}.csv"),
                  ["x", "y", "p_laplace", "dp_storage", "p_laplace_plus_dp", "p_neopz"],
                  [X[:, 0], X[:, 1], plap, dp, pred, pv])
        # where is the volumetric change concentrated
        de = divs.mean(1)
        o = np.argsort(-np.abs(de))[:3]
        print(f"[{name}] largest |delta div u| elements (centroid, mean delta div u):",
              [(tuple(np.round(X[T[k]].mean(0), 2)), f"{de[k]:.3e}") for k in o])


if __name__ == "__main__":
    main()
