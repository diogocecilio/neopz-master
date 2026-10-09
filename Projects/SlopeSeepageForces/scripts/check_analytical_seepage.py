#!/usr/bin/env python3
"""Independent (adversarial) check of analytical_seepage.py (Ceron et al., IJNAMG 2025, Eqs. 29-40).

Nothing below reuses the closed forms of analytical_seepage.py (A, B, C, D, m, h1..h4, Eq. 40); the module is
only queried through its public API (velocity, polar_velocity, force, Jstar, F, m, h2/dh2 for the printed-Eq. 31
test).  Checks:

  A) h2 problem by a Chebyshev-Ritz method (no ODE, no shooting):
     F = sup_h2 h2(Theta)^2 / (sqrt C + sqrt D)^2 = 1 / min_{m>=0} Phi(m),
     Phi(m) = (1 + 1/m) min_{h2(Theta)=1} [int c h2'^2 + m int d h2^2],  Phi(0+) = int_0^Theta d = A.
     Small-m expansion Phi(m) = A + m (A - P) + O(m^2), P = int_0^Theta (int_0^t d)^2 / c dt: the m -> 0 limit
     is a local optimum iff A > P (threshold angles for the degenerate case).
  B) J* of the IMPLEMENTED field evaluated by an independent 2-D quadrature of Eq. 23 (Cartesian v, Cartesian
     K^-1, Cartesian outward normals, u^d of Eq. 21 on the face and on the toe ground).
  C) Direct minimisation of J* (Eq. 23) over a discretisation of the class (29) with scipy.optimize (BFGS):
     g(r) = r h1(r): C1 cubic Hermite on a geometrically graded mesh in s = r/R_w (g(0) = g(R_w) = 0);
     h2(theta): Chebyshev series normalised by h2(Theta) = 1; h3(r): Chebyshev in r;
     h4(r) = (1/r) * Chebyshev in w = ln(x_toe(r) + H) (x_toe(r) = sqrt(r^2 - H^2)).
     Compare J*_direct with Jstar() and the direct-optimum velocity with velocity() at sample points.
  D) The printed Eq. 31 (prefactor 1/(C - D), e_theta exponent and coefficient C/D) vs the corrected form.
  E) API contract of force(): shapes, zero outside the soil and for r >= R_e, h_w = 0, f = K^-1 v, timing.

Run:  python3 Projects/SlopeSeepageForces/scripts/check_analytical_seepage.py [--quick] [--sweep]
"""
from __future__ import annotations

import os

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import argparse
import sys
import time

import numpy as np
from numpy.polynomial import chebyshev as Ch
from scipy.integrate import quad
from scipy.optimize import brentq, minimize, minimize_scalar

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPT_DIR)
from analytical_seepage import AnalyticalSeepage, NoSeepage  # noqa: E402

RESULTS_DIR = os.path.join(os.path.dirname(SCRIPT_DIR), "results", "analytical_seepage")


def gl(a, b, n):
    x, w = np.polynomial.legendre.leggauss(n)
    return 0.5 * (b - a) * (x + 1.0) + a, 0.5 * (b - a) * w


def cheb_vander_d(t, N):
    """T_k(t) and dT_k/dt, k = 0..N."""
    V = Ch.chebvander(t, N)
    dV = np.zeros_like(V)
    for k in range(1, N + 1):
        e = np.zeros(N + 1)
        e[k] = 1.0
        dV[:, k] = Ch.chebval(t, Ch.chebder(e))
    return V, dV


def polar_to_xy(r, th):
    """Fig. 3 polar frame mapped to paper coordinates (x right, y down)."""
    return -r * np.cos(th), r * np.sin(th)


def er_et(th):
    """e_r = d(x,y)/dr, e_theta = (1/r) d(x,y)/dtheta."""
    return (-np.cos(th), np.sin(th)), (np.sin(th), np.cos(th))


def cart(vr, vt, th):
    (erx, ery), (etx, ety) = er_et(th)
    return vr * erx + vt * etx, vr * ery + vt * ety


# ================================================================================================
# A) h2 problem by Chebyshev-Ritz
# ================================================================================================
class RitzH2:
    def __init__(self, beta_deg, alpha, N=28, nq=120):
        self.Theta = np.pi - np.radians(beta_deg)
        self.alpha = alpha
        th, w = gl(0.0, self.Theta, nq)
        t = 2.0 * th / self.Theta - 1.0
        V, dV = cheb_vander_d(t, N)
        self.P = V[:, 1:] - 1.0                 # h2 = 1 + P c  ->  h2(Theta) = 1
        self.Pd = dV[:, 1:] * 2.0 / self.Theta
        s2 = np.sin(th) ** 2
        self.c = (1 - s2) + alpha * s2           # weight of v_r^2  (cos^2 + alpha sin^2)
        self.d = s2 + alpha * (1 - s2)           # weight of v_t^2  (sin^2 + alpha cos^2)
        self.w = w
        self.th = th
        self.Kc = self.Pd.T @ (self.Pd * (w * self.c)[:, None])
        self.Md = self.P.T @ (self.P * (w * self.d)[:, None])
        self.fd = self.P.T @ (w * self.d)
        self.A = np.sum(w * self.d)

    def V(self, m):
        cc = np.linalg.solve(self.Kc + m * self.Md, -m * self.fd)
        h = 1.0 + self.P @ cc
        hp = self.Pd @ cc
        C = np.sum(self.w * self.c * hp ** 2)
        D = np.sum(self.w * self.d * h ** 2)
        return C + m * D, C, D

    def Phi(self, m):
        return self.A if m <= 0 else (1.0 + 1.0 / m) * self.V(m)[0]

    def solve(self):
        lg = np.linspace(np.log(1e-6), np.log(30.0), 70)
        ph = np.array([self.Phi(np.exp(t)) for t in lg])
        i = int(np.argmin(ph))
        if ph[i] >= self.A:
            return 0.0, self.A, True
        res = minimize_scalar(lambda t: self.Phi(np.exp(t)), bounds=(lg[max(i - 1, 0)], lg[min(i + 1, len(lg) - 1)]),
                              method="bounded", options={"xatol": 1e-10})
        m = float(np.exp(res.x))
        return m, self.Phi(m), False


def small_m_slope(beta_deg, alpha):
    """A - P (Phi(m) = A + m (A - P) + O(m^2)); >0 -> Phi increases from m = 0 (degenerate limit is a local min)."""
    Th = np.pi - np.radians(beta_deg)

    def d(t):
        return np.sin(t) ** 2 + alpha * np.cos(t) ** 2

    def c(t):
        return np.cos(t) ** 2 + alpha * np.sin(t) ** 2

    def Dcum(t):
        return quad(d, 0.0, t, epsabs=1e-14, epsrel=1e-13)[0]

    A = Dcum(Th)
    P = quad(lambda t: Dcum(t) ** 2 / c(t), 0.0, Th, epsabs=1e-13, epsrel=1e-12)[0]
    return A - P


# ================================================================================================
# B/C) J* of Eq. 23 for fields of the class (29) by brute-force quadrature
# ================================================================================================
class ClassFunctional:
    """J*(v') = 1/2 int v.K^-1.v + int_{dOmega_u} u^d v.n for v of the class (29), discretised as described in the
    module docstring, with analytic gradient.  Also evaluates J* of any velocity callable (polar or Cartesian)
    on the same quadrature (check B)."""

    def __init__(self, beta_deg, H, hw, alpha, kh=1.0, gw=9.81, Lm_over_H=10.0,
                 n_el=100, sigma=0.5, ng=8, nth=40, N2=12, N3=10, N4=16, g0_free=False):
        b = np.radians(beta_deg)
        self.sb, self.cb = np.sin(b), np.cos(b)
        self.H, self.hw, self.al, self.kh, self.gw = H, hw, alpha, kh, gw
        self.Th = np.pi - b
        self.Rw, self.R = hw / self.sb, H / self.sb
        self.xtoe = H * self.cb / self.sb
        self.Lm = Lm_over_H * H
        self.n_face = np.array([self.sb, -self.cb])          # outward normal of the face
        # ---------------- zone 1: graded Hermite mesh in s
        nodes = np.concatenate([[0.0], sigma ** np.arange(n_el - 1, -1, -1.0)])   # s_0 = 0 < ... < s_n = 1
        n = len(nodes) - 1
        self.nodes = nodes
        hs = np.diff(nodes)
        dscale = np.concatenate([[hs[0]], hs])                # derivative dof = g_s(s_j) * local size
        val_idx = {j: j - 1 for j in range(1, n)}             # value dofs (g(0) = g(1) = 0)
        if g0_free:
            val_idx[0] = n - 1                                # extra dof: g(0) (outside the class, check only)
        nval = len(val_idx)
        der_idx = {j: nval + j for j in range(n + 1)}
        self.nq_g = nval + n + 1
        xg, wg = np.polynomial.legendre.leggauss(ng)
        S, WS, PG, PGs = [], [], [], []
        for e in range(n):
            a_, b_ = nodes[e], nodes[e + 1]
            h = b_ - a_
            xi = 0.5 * (xg + 1.0)
            S.append(a_ + h * xi)
            WS.append(0.5 * h * wg)
            H00, H10, H01, H11 = 2 * xi ** 3 - 3 * xi ** 2 + 1, xi ** 3 - 2 * xi ** 2 + xi, -2 * xi ** 3 + 3 * xi ** 2, xi ** 3 - xi ** 2
            D00, D10, D01, D11 = (6 * xi ** 2 - 6 * xi) / h, 3 * xi ** 2 - 4 * xi + 1, (-6 * xi ** 2 + 6 * xi) / h, 3 * xi ** 2 - 2 * xi
            pg = np.zeros((ng, self.nq_g))
            pgs = np.zeros((ng, self.nq_g))
            for j, (Hv, Dv, Hd, Dd) in ((e, (H00, D00, H10, D10)), (e + 1, (H01, D01, H11, D11))):
                if j in val_idx:
                    pg[:, val_idx[j]] += Hv
                    pgs[:, val_idx[j]] += Dv
                k = der_idx[j]
                pg[:, k] += Hd * h / dscale[j]
                pgs[:, k] += Dd / dscale[j]
            PG.append(pg)
            PGs.append(pgs)
        self.s = np.concatenate(S)
        self.ws = np.concatenate(WS)
        self.PG = np.vstack(PG)
        self.PGs = np.vstack(PGs)
        # theta quadrature on [0, Theta] and h2 basis
        self.th1, self.wth1 = gl(0.0, self.Th, nth)
        V, dV = cheb_vander_d(2 * self.th1 / self.Th - 1.0, N2)
        self.Ph = V[:, 1:] - 1.0
        self.Phd = dV[:, 1:] * 2.0 / self.Th
        _, dVe = cheb_vander_d(np.array([1.0]), N2)
        self.Phd_e = dVe[0, 1:] * 2.0 / self.Th
        self.nc = N2
        # ---------------- zone 2
        self.has2 = self.R > self.Rw * (1 + 1e-12)
        if self.has2:
            self.r2, self.wr2 = gl(self.Rw, self.R, 24)
            self.P3 = Ch.chebvander(2 * (self.r2 - self.Rw) / (self.R - self.Rw) - 1.0, N3)
        self.nb = N3 + 1 if self.has2 else 0
        # ---------------- zone 3 (r from x_toe(r) = e^w - H)
        w0, w1 = np.log(self.xtoe + H), np.log(self.xtoe + self.Lm + H)
        W, WW = [], []
        for k in range(4):
            ww, wq = gl(w0 + (w1 - w0) * k / 4, w0 + (w1 - w0) * (k + 1) / 4, 20)
            W.append(ww)
            WW.append(wq)
        w = np.concatenate(W)
        ww = np.concatenate(WW)
        x = np.exp(w) - H
        self.x3 = x
        self.r3 = np.hypot(x, H)
        self.wx3 = ww * np.exp(w)                         # dx
        self.wr3 = self.wx3 * x / self.r3                 # dr = (x / r) dx
        self.thm = np.arctan2(H, -x)                      # toe-ground polar angle (independent of Eq. 35)
        t3 = 2 * (w - w0) / (w1 - w0) - 1.0
        self.P4 = Ch.chebvander(t3, N4) / self.r3[:, None]
        xg3, wg3 = np.polynomial.legendre.leggauss(nth)
        self.th3 = 0.5 * self.thm[:, None] * (xg3[None, :] + 1.0)
        self.wth3 = 0.5 * self.thm[:, None] * wg3[None, :]
        self.nd = N4 + 1
        self.Re = np.hypot(H, self.xtoe + self.Lm)

    # ---------------------------------------------------------------- parameter vector
    def split(self, p):
        i = 0
        q = p[i:i + self.nq_g]; i += self.nq_g
        c = p[i:i + self.nc]; i += self.nc
        bb = p[i:i + self.nb]; i += self.nb
        d = p[i:i + self.nd]
        return q, c, bb, d

    def initial(self, rng=None):
        nodes = self.nodes
        q = np.zeros(self.nq_g)
        # g(s) ~ -2 R_w s (1 - s): values at interior nodes, scaled derivatives
        n = len(nodes) - 1
        hs = np.diff(nodes)
        dscale = np.concatenate([[hs[0]], hs])
        q[:n - 1] = -2 * self.Rw * nodes[1:n] * (1 - nodes[1:n])
        q[n - 1 + (self.nq_g - (n - 1) - (n + 1)):][: n + 1] = -2 * self.Rw * (1 - 2 * nodes) * dscale
        c = np.zeros(self.nc)
        bb = np.zeros(self.nb)
        if self.nb:
            bb[0] = 3.0
        d = np.zeros(self.nd)
        d[0] = 1.5
        p = np.concatenate([q, c, bb, d])
        if rng is not None:
            p = p + 0.3 * rng.standard_normal(p.size) * (np.abs(p) + 0.1)
        return p

    # ---------------------------------------------------------------- energy helper
    def _energy(self, vx, vy, W):
        Ex = W * vx / self.kh
        Ey = self.al * W * vy / self.kh
        return 0.5 * np.sum(vx * Ex + vy * Ey), Ex, Ey

    def J(self, p, grad=True):
        q, c, bb, d = self.split(p)
        Rw, al, kh, gw = self.Rw, self.al, self.kh, self.gw
        gq = np.zeros_like(q); gc = np.zeros_like(c); gb = np.zeros_like(bb); gd = np.zeros_like(d)
        J = 0.0
        # ---------------- zone 1
        if Rw > 0:
            G = self.PG @ q
            Gs = self.PGs @ q
            h = 1.0 + self.Ph @ c
            hp = self.Phd @ c
            s = self.s
            fr = -(G / (Rw * s))                 # v_r = fr * h2'
            ft = Gs / Rw                         # v_t = ft * h2
            th = self.th1
            vr = fr[:, None] * hp[None, :]
            vt = ft[:, None] * h[None, :]
            (erx, ery), (etx, ety) = er_et(th)
            vx = vr * erx + vt * etx
            vy = vr * ery + vt * ety
            W = (self.ws * Rw ** 2 * s)[:, None] * self.wth1[None, :]
            E, Ex, Ey = self._energy(vx, vy, W)
            J += E
            Er = Ex * erx + Ey * ery
            Et = Ex * etx + Ey * ety
            dfr = Er @ hp
            dft = Et @ h
            dhp = fr @ Er
            dh = ft @ Et
            # face 0 < r < R_w: u^d = -gw y = -gw r sin(beta), n = (sin b, -cos b)
            hpe = self.Phd_e @ c
            (erx, ery), (etx, ety) = er_et(self.Th)
            vr_f = fr * hpe
            vt_f = ft * 1.0
            vn = (vr_f * erx + vt_f * etx) * self.n_face[0] + (vr_f * ery + vt_f * ety) * self.n_face[1]
            bn = self.ws * Rw * (-gw * Rw * s * self.sb)
            J += np.sum(bn * vn)
            dvr_f = bn * (erx * self.n_face[0] + ery * self.n_face[1])
            dvt_f = bn * (etx * self.n_face[0] + ety * self.n_face[1])
            dfr = dfr + dvr_f * hpe
            dft = dft + dvt_f
            dhpe = np.sum(dvr_f * fr)
            if grad:
                dG = dfr * (-1.0 / (Rw * s))
                dGs = dft / Rw
                gq = self.PG.T @ dG + self.PGs.T @ dGs
                gc = self.Ph.T @ dh + self.Phd.T @ dhp + self.Phd_e * dhpe
        # ---------------- zone 2: v = h3 e_theta
        if self.has2:
            h3 = self.P3 @ bb
            th = self.th1
            (_, _), (etx, ety) = er_et(th)
            vx = h3[:, None] * etx[None, :]
            vy = h3[:, None] * ety[None, :]
            W = (self.wr2 * self.r2)[:, None] * self.wth1[None, :]
            E, Ex, Ey = self._energy(vx, vy, W)
            dh3 = (Ex * etx + Ey * ety).sum(axis=1)
            (_, _), (etx, ety) = er_et(self.Th)
            vn_unit = etx * self.n_face[0] + ety * self.n_face[1]
            bn = self.wr2 * (-gw * self.hw)
            J += E + np.sum(bn * h3 * vn_unit)
            dh3 = dh3 + bn * vn_unit
            gb = self.P3.T @ dh3
        # ---------------- zone 3: v = h4 e_theta on 0 < theta < theta_m(r); toe ground y = H, n = (0, -1)
        if self.hw > 0:
            h4 = self.P4 @ d
            th = self.th3
            (_, _), (etx, ety) = er_et(th)
            vx = h4[:, None] * etx
            vy = h4[:, None] * ety
            W = (self.wr3 * self.r3)[:, None] * self.wth3
            E, Ex, Ey = self._energy(vx, vy, W)
            dh4 = (Ex * etx + Ey * ety).sum(axis=1)
            (_, _), (etx, ety) = er_et(self.thm)
            vn_unit = etx * 0.0 + ety * (-1.0)
            bn = self.wx3 * (-gw * self.hw)
            J += E + np.sum(bn * h4 * vn_unit)
            dh4 = dh4 + bn * vn_unit
            gd = self.P4.T @ dh4
        if not grad:
            return J
        return J, np.concatenate([gq, gc, gb, gd])

    def velocity_polar(self, p, r, th):
        """(v_r, v_theta) of the discretised field at polar points (inside the soil)."""
        q, c, bb, d = self.split(p)
        r = np.asarray(r, float)
        th = np.asarray(th, float)
        vr = np.zeros(r.shape)
        vt = np.zeros(r.shape)
        z1 = r < self.Rw
        if np.any(z1):
            s = r[z1] / self.Rw
            # evaluate the Hermite field at s
            G, Gs = self._hermite(q, s)
            tt = 2 * np.clip(th[z1], 0, self.Th) / self.Th - 1
            V, dV = cheb_vander_d(tt, self.nc)
            h = 1.0 + (V[:, 1:] - 1.0) @ c
            hp = (dV[:, 1:] * 2.0 / self.Th) @ c
            vr[z1] = -(G / (self.Rw * s)) * hp
            vt[z1] = Gs / self.Rw * h
        z2 = (r >= self.Rw) & (r < self.R)
        if np.any(z2):
            vt[z2] = Ch.chebvander(2 * (r[z2] - self.Rw) / (self.R - self.Rw) - 1.0, self.nb - 1) @ bb
        z3 = (r >= self.R) & (r < self.Re)
        if np.any(z3):
            x = np.sqrt(r[z3] ** 2 - self.H ** 2)
            w = np.log(x + self.H)
            w0, w1 = np.log(self.xtoe + self.H), np.log(self.xtoe + self.Lm + self.H)
            vt[z3] = Ch.chebvander(2 * (w - w0) / (w1 - w0) - 1.0, self.nd - 1) @ d / r[z3]
        return vr, vt

    def _hermite(self, q, s):
        nodes = self.nodes
        n = len(nodes) - 1
        e = np.clip(np.searchsorted(nodes, s, side="right") - 1, 0, n - 1)
        # rebuild by sampling the precomputed element matrices is awkward: re-evaluate the Hermite basis
        hs = np.diff(nodes)
        dscale = np.concatenate([[hs[0]], hs])
        nval = self.nq_g - (n + 1)
        g0_free = nval == n
        vals = np.zeros(n + 1)
        vals[1:n] = q[:n - 1]
        if g0_free:
            vals[0] = q[n - 1]
        ders = q[nval:] / dscale
        a_, h = nodes[e], hs[e]
        xi = (s - a_) / h
        H00, H10, H01, H11 = 2 * xi ** 3 - 3 * xi ** 2 + 1, xi ** 3 - 2 * xi ** 2 + xi, -2 * xi ** 3 + 3 * xi ** 2, xi ** 3 - xi ** 2
        D00, D10, D01, D11 = (6 * xi ** 2 - 6 * xi) / h, 3 * xi ** 2 - 4 * xi + 1, (-6 * xi ** 2 + 6 * xi) / h, 3 * xi ** 2 - 2 * xi
        G = vals[e] * H00 + ders[e] * h * H10 + vals[e + 1] * H01 + ders[e + 1] * h * H11
        Gs = vals[e] * D00 + ders[e] * D10 + vals[e + 1] * D01 + ders[e + 1] * D11
        return G, Gs

    # ---------------------------------------------------------------- J* of an arbitrary field (check B)
    def J_of_field(self, vel_xy, vel_polar_boundary):
        """vel_xy(x, y) -> (vx, vy) at interior points (Cartesian API of the module);
        vel_polar_boundary(r, th) -> (v_r, v_theta) for boundary points (avoids in/out rounding on the boundary)."""
        Rw, gw = self.Rw, self.gw
        J = 0.0
        parts = []
        if Rw > 0:
            r = Rw * self.s[:, None] * np.ones_like(self.th1)[None, :]
            th = np.ones_like(self.s)[:, None] * self.th1[None, :]
            vx, vy = vel_xy(*polar_to_xy(r, th))
            W = (self.ws * Rw ** 2 * self.s)[:, None] * self.wth1[None, :]
            E1 = self._energy(vx, vy, W)[0]
            rf = Rw * self.s
            vr, vt = vel_polar_boundary(rf, np.full_like(rf, self.Th))
            vx, vy = cart(vr, vt, self.Th)
            B1 = np.sum(self.ws * Rw * (-gw * rf * self.sb) * (vx * self.n_face[0] + vy * self.n_face[1]))
            parts += [E1, B1]
        if self.has2:
            r = self.r2[:, None] * np.ones_like(self.th1)[None, :]
            th = np.ones_like(self.r2)[:, None] * self.th1[None, :]
            vx, vy = vel_xy(*polar_to_xy(r, th))
            W = (self.wr2 * self.r2)[:, None] * self.wth1[None, :]
            E2 = self._energy(vx, vy, W)[0]
            vr, vt = vel_polar_boundary(self.r2, np.full_like(self.r2, self.Th))
            vx, vy = cart(vr, vt, self.Th)
            B2 = np.sum(self.wr2 * (-gw * self.hw) * (vx * self.n_face[0] + vy * self.n_face[1]))
            parts += [E2, B2]
        r = self.r3[:, None] * np.ones_like(self.th3)
        vx, vy = vel_xy(*polar_to_xy(r, self.th3))
        W = (self.wr3 * self.r3)[:, None] * self.wth3
        E3 = self._energy(vx, vy, W)[0]
        vr, vt = vel_polar_boundary(self.r3, self.thm)
        vx, vy = cart(vr, vt, self.thm)
        B3 = np.sum(self.wx3 * (-gw * self.hw) * (vx * 0.0 + vy * (-1.0)))
        parts += [E3, B3]
        return float(sum(parts)), parts


def direct_minimisation(cf, p0, maxiter=20000):
    res = minimize(cf.J, p0, jac=True, method="BFGS", options={"gtol": 1e-10, "maxiter": maxiter})
    return res


def compare_fields(md, cf, p, H, zones=(1, 2, 3)):
    """Velocity of the direct optimum (parameters p of cf) vs md.velocity() at sample points of each zone."""
    rows = []
    for zone, rr, fr in ((1, md.Rw * np.array([0.1, 0.3, 0.6, 0.9]), (0.05, 0.5, 0.95)),
                         (2, md.Rw + (md.R - md.Rw) * np.array([0.2, 0.5, 0.8]), (0.05, 0.5, 0.95)),
                         (3, np.array([1.3 * md.R, 3.0 * md.R, 0.8 * md.Re]), (0.05, 0.5, 0.95))):
        if zone not in zones or (zone == 2 and md.R <= md.Rw):
            continue
        for r in rr:
            thmax = md.Theta if zone < 3 else np.pi - np.arcsin(H / r)
            for f in fr:
                th = f * thmax
                x, y = polar_to_xy(r, th)
                vx_m, vy_m = md.velocity(x, y)
                vr_d, vt_d = cf.velocity_polar(p, np.array([r]), np.array([th]))
                vx_d, vy_d = cart(vr_d[0], vt_d[0], th)
                rows.append((zone, r, th, float(vx_m), float(vy_m), vx_d, vy_d))
    rows = np.array(rows)
    for zone in zones:
        z = rows[rows[:, 0] == zone] if len(rows) else rows
        if len(z) == 0:
            continue
        dv = np.hypot(z[:, 3] - z[:, 5], z[:, 4] - z[:, 6])
        vm = np.hypot(z[:, 3], z[:, 4])
        print(f"     zone {zone}: max |v_direct - v_module| / max|v_module| = {dv.max() / vm.max():.2e}"
              f"  (max pointwise rel {np.max(dv / np.maximum(vm, 1e-12 * vm.max())):.2e}, {len(z)} points)")


# ================================================================================================
# driver
# ================================================================================================
def main(quick=False, sweep=False):
    os.makedirs(RESULTS_DIR, exist_ok=True)
    H, gw = 1.0, 9.81
    t_start = time.time()

    print("=" * 100)
    print("A) h2 problem by Chebyshev-Ritz (no ODE): m, F = 1/Phi_min vs module; small-m slope A - P")
    print("=" * 100)
    cases = [(15, 1), (30, 1), (60, 1), (75, 1), (80, 1), (85, 1), (90, 1), (30, 2), (85, 2), (60, 4), (85, 4),
             (90, 4), (15, 10), (60, 10), (90, 10)]
    print(f"{'beta':>5} {'alpha':>5} {'m Ritz':>10} {'m module':>10} {'F Ritz':>11} {'F module':>11} {'rel':>9} {'A-P':>9}")
    for b, a in cases:
        rz = RitzH2(b, a)
        m, phi, deg = rz.solve()
        md = AnalyticalSeepage(b, H, H, a, gamma_w=gw)
        F = 1.0 / phi
        print(f"{b:5d} {a:5d} {m:10.6f} {md.m:10.6f} {F:11.8f} {md.F:11.8f} {(F - md.F) / md.F:9.1e} {small_m_slope(b, a):+9.4f}"
              + ("  degenerate (m->0)" if deg else ""))
    # thresholds of the local degeneracy A = P
    print("  threshold beta (A = P, Phi'(0+) = 0): ", end="")
    for a in (1, 2, 4, 5, 10):
        f = lambda bb: small_m_slope(bb, a)
        if f(89.999) < 0:
            print(f"alpha={a}: none <90;  ", end="")
        else:
            print(f"alpha={a}: {brentq(f, 60.0, 89.999, xtol=1e-6):.3f} deg;  ", end="")
    print(f"\n  analytic alpha=1: pi - sqrt(3) = {np.degrees(np.pi - np.sqrt(3.0)):.3f} deg")
    if sweep:
        # Fig. 9 range: beta = 15..90 deg every 1 deg, alpha = 1, 5, 10 (h2 problem only; ~2 min)
        for a in (1, 5, 10):
            worst, ms, mism, degs = 0.0, [], 0, []
            for b in range(15, 91):
                md = AnalyticalSeepage(b, H, H, a, gamma_w=gw)
                m, phi, deg = RitzH2(b, a).solve()
                worst = max(worst, abs(1.0 / phi - md.F) / md.F)
                mism += int(deg != md.degenerate)
                ms.append(md.m)
                if md.degenerate:
                    degs.append(b)
            print(f"  sweep alpha={a:2d}, beta=15..90 by 1 deg: max rel |F_Ritz - F_module| = {worst:.1e}, "
                  f"degenerate-flag mismatches = {mism}, max |m(b+1)-m(b)| = {np.max(np.abs(np.diff(ms))):.4f}, "
                  f"degenerate betas = {degs if degs else 'none'}")

    print()
    print("=" * 100)
    print("B) J* of the IMPLEMENTED field by independent 2-D quadrature (Cartesian v, K^-1, normals) vs Jstar()")
    print("=" * 100)
    bcases = [(30, 1, 0.5), (30, 1, 1.0), (60, 4, 1.0), (20, 10, 0.8), (75, 2, 0.3), (85, 1, 1.0), (90, 10, 0.6)]
    for (b, a, hwr) in bcases:
        md = AnalyticalSeepage(b, H, hwr * H, a, gamma_w=gw)
        cf = ClassFunctional(b, H, hwr * H, a, gw=gw)
        Jq, parts = cf.J_of_field(md.velocity, md.polar_velocity)
        print(f"  beta={b:3d} alpha={a:3d} hw/H={hwr:3.1f} m={md.m:.5f}: J*(quadrature) = {Jq:+.10f}  Jstar() = {md.Jstar():+.10f}"
              f"  rel = {(Jq - md.Jstar()) / abs(md.Jstar()):+.1e}")
    # Fig. 5 (h_w = H, L_m = 10 H) with the independent quadrature vs the paper's solid curves (vector data)
    path = os.path.join(os.path.dirname(SCRIPT_DIR), "data", "fig5_vector_fill_polygons.csv")
    if os.path.exists(path):
        ref = {(int(r[0]), int(round(r[1]))): r[2] for r in np.loadtxt(path, delimiter=",", comments="#")}
        print("  Fig. 5 solid curve, -J*/(k_h H^2 gw^2) by this quadrature (paper value):")
        for a in (1, 2, 4, 10):
            line = []
            for b in (15, 45, 80, 85, 90):
                md = AnalyticalSeepage(b, H, H, a, gamma_w=gw)
                Jq, _ = ClassFunctional(b, H, H, a, gw=gw).J_of_field(md.velocity, md.polar_velocity)
                line.append(f"b={b}: {-Jq / gw ** 2:.4f} ({ref.get((a, b), np.nan):.4f})")
            print(f"    alpha={a:2d}  " + "  ".join(line))

    print()
    print("=" * 100)
    print("C) direct minimisation of J* over the discretised class (BFGS), vs Jstar() and velocity()")
    print("=" * 100)
    ccases = [(30, 1, 0.5), (60, 4, 1.0), (20, 10, 0.8), (85, 1, 1.0)] if not quick else [(30, 1, 0.5), (85, 1, 1.0)]
    rng = np.random.default_rng(7)
    for (b, a, hwr) in ccases:
        md = AnalyticalSeepage(b, H, hwr * H, a, gamma_w=gw)
        cf = ClassFunctional(b, H, hwr * H, a, gw=gw)
        # gradient check (central differences, a few components)
        p0 = cf.initial(rng)
        _, g0 = cf.J(p0)
        idx = rng.choice(p0.size, 6, replace=False)
        fd = []
        for i in idx:
            e = np.zeros_like(p0)
            e[i] = 1e-6 * max(1.0, abs(p0[i]))
            fd.append((cf.J(p0 + e, grad=False) - cf.J(p0 - e, grad=False)) / (2 * e[i]))
        gerr = np.max(np.abs(np.array(fd) - g0[idx])) / (np.max(np.abs(g0[idx])) + 1e-30)
        t0 = time.time()
        best = None
        for start in (cf.initial(), cf.initial(rng)):
            res = direct_minimisation(cf, start)
            if best is None or res.fun < best.fun:
                best = res
            print(f"    start: J = {res.fun:+.10f}  ({res.nit} it, |grad| = {np.max(np.abs(res.jac)):.1e}, {res.message})")
        Jd = best.fun
        print(f"  beta={b} alpha={a} hw/H={hwr}: grad-check {gerr:.1e};  J*_direct = {Jd:+.10f}  Jstar() = {md.Jstar():+.10f}"
              f"  rel = {(Jd - md.Jstar()) / abs(md.Jstar()):+.2e}  ({time.time() - t0:.1f} s)"
              + ("  [module: degenerate m->0]" if md.degenerate else ""))
        compare_fields(md, cf, best.x, H)
        if (b, a, hwr) == (30, 1, 0.5) and not quick:
            # refinement: the zone-1 discrepancy is discretisation error of the direct solution (it shrinks)
            cfr = ClassFunctional(b, H, hwr * H, a, gw=gw, sigma=0.65, n_el=170)
            rr_ = direct_minimisation(cfr, cfr.initial())
            print(f"     refined mesh (sigma=0.65, 170 el.): J*_direct = {rr_.fun:+.10f}"
                  f"  rel = {(rr_.fun - md.Jstar()) / abs(md.Jstar()):+.2e}")
            compare_fields(md, cfr, rr_.x, H, zones=(1,))
        if md.degenerate:
            # (i) the m -> 0 limit field itself: h2 = const (C = 0), g(0) free, i.e. v = a(r) e_theta in zone 1.
            #     (g(0) free with C > 0 is NOT tested: its energy diverges like C g(0)^2 ln(1/r_min) and a
            #     truncated quadrature would report spuriously low values.)
            cf0 = ClassFunctional(b, H, hwr * H, a, gw=gw, N2=0, g0_free=True)
            res0 = direct_minimisation(cf0, cf0.initial())
            print(f"     h2 = const, g(0) free (m -> 0 limit field): J*_direct = {res0.fun:+.10f}"
                  f"  rel = {(res0.fun - md.Jstar()) / abs(md.Jstar()):+.2e}")
            # (ii) inside the class (g(0) = 0): the infimum is approached from above, only logarithmically in
            #      the smallest radius the discretisation can resolve (the infimum is not attained in the class)
            for n_el in ((40, 100) if quick else (40, 100, 250)):
                cfn = ClassFunctional(b, H, hwr * H, a, gw=gw, n_el=n_el)
                rn = direct_minimisation(cfn, cfn.initial())
                rel = (rn.fun - md.Jstar()) / abs(md.Jstar())
                lns = np.log(1.0 / cfn.nodes[1])
                print(f"     class with g(0)=0, mesh down to s_1 = {cfn.nodes[1]:.1e}: J*_direct = {rn.fun:+.10f}"
                      f"  rel = {rel:+.2e}  (rel * ln(1/s_1) = {rel * lns:.4f})")

    print()
    print("=" * 100)
    print("D) printed Eq. 31 vs corrected form (beta=30, alpha=1, hw=H), same quadrature as B)")
    print("=" * 100)
    b, a = 30, 1
    md = AnalyticalSeepage(b, H, H, a, gamma_w=gw)
    cf = ClassFunctional(b, H, H, a, gw=gw)
    K0 = md.kh * gw * md.sb * md.h2e
    r_CD = md.C / md.D
    m = np.sqrt(r_CD)

    def make(pref, ex_t, co_t, ex_r):
        def vpol(r, th):
            r = np.asarray(r, float)
            th = np.asarray(th, float)
            vr, vt = md.polar_velocity(r, th)          # zones 2-3 unchanged
            z1 = r < md.Rw
            s = r[z1] / md.Rw
            t1 = np.clip(th[z1], 0, md.Theta)
            vt[z1] = pref * (1 - co_t * s ** (ex_t - 1)) * md.h2(t1)
            vr[z1] = -pref * (1 - s ** (ex_r - 1)) * md.dh2(t1)
            return vr, vt

        def vxy(x, y):
            r = np.hypot(x, y)
            th = np.arctan2(y, -x)
            vr, vt = vpol(r, th)
            return cart(vr, vt, th)
        return vxy, vpol

    for label, args in (("printed (1/(C-D), C/D in e_theta)", (K0 / (md.C - md.D), r_CD, r_CD, m)),
                        ("sign fixed only (1/(D-C), C/D in e_theta)", (K0 / (md.D - md.C), r_CD, r_CD, m)),
                        ("exponent fixed only (1/(C-D), sqrt)", (K0 / (md.C - md.D), m, m, m)),
                        ("corrected (1/(D-C), sqrt(C/D))", (K0 / (md.D - md.C), m, m, m))):
        vxy, vpol = make(*args)
        Jq, _ = cf.J_of_field(vxy, vpol)
        # divergence residual in zone 1 (polar formula, central differences)
        rr = md.Rw * np.linspace(0.1, 0.9, 9)[:, None] * np.ones(7)[None, :]
        tt = md.Theta * np.linspace(0.1, 0.9, 7)[None, :] * np.ones(9)[:, None]
        hr, ht = 1e-6 * rr, 1e-6
        vr_p, _ = vpol(rr + hr, tt)
        vr_m, _ = vpol(rr - hr, tt)
        _, vt_p = vpol(rr, tt + ht)
        _, vt_m = vpol(rr, tt - ht)
        div = ((rr + hr) * vr_p - (rr - hr) * vr_m) / (2 * hr) / rr + (vt_p - vt_m) / (2 * ht) / rr
        vr0, vt0 = vpol(rr, tt)
        sc = np.max(np.hypot(vr0, vt0) / rr)
        print(f"  {label:45s} J* = {Jq:+.6f}  (optimum {md.Jstar():+.6f});  max|div v|/max(|v|/r) = {np.max(np.abs(div)) / sc:.1e}")

    print()
    print("=" * 100)
    print("E) API contract of force()")
    print("=" * 100)
    md = AnalyticalSeepage(45, 5.0, 3.0, 4.0, kh=2.5, gamma_w=9.81)
    xs = np.array([[-2.0, 1.0], [3.0, 7.0]])
    ys = np.array([[1.0, 3.0], [1.0, 6.0]])
    fx, fy = md.force(xs, ys)
    vx, vy = md.velocity(xs, ys)
    print(f"  shapes: in {xs.shape} -> out {fx.shape}, {fy.shape}; scalar -> {np.shape(md.force(1.0, 2.0)[0])};"
          f" broadcast (2,2)x() -> {md.force(xs, 2.0)[0].shape}")
    print(f"  f = K^-1 v: max |fx - vx/kh| = {np.max(np.abs(fx - vx / 2.5)):.1e}, max |fy - alpha vy/kh| = {np.max(np.abs(fy - 4 * vy / 2.5)):.1e}")
    out_pts = (np.array([-1.0, 1.0, 6.0, 0.5, -60.0, 1.0]), np.array([-0.01, 0.5, 4.9, -1.0, 2.0, 60.0]))
    print(f"  outside soil (above crest/face/toe ground) and r >= R_e: f = {np.array(md.force(*out_pts)).ravel()}")
    f_kh1 = AnalyticalSeepage(45, 5.0, 3.0, 4.0, kh=1.0).force(xs, ys)
    print(f"  k_h independence of f: {np.max(np.abs(np.array(f_kh1) - np.array([fx, fy]))):.1e}")
    print(f"  hw = 0: max|f| = {np.max(np.abs(AnalyticalSeepage(45, 5.0, 0.0, 4.0).force(xs, ys))):.1e};"
          f"  NoSeepage: {NoSeepage().force(xs, ys)[0].shape}")
    # sign / units: far under the crest f ~ +y (downwards), magnitude O(gamma_w h_w / r)
    fx1, fy1 = md.force(-3.0, 0.3)
    print(f"  under the crest at (-3, 0.3): f = ({float(fx1):+.3f}, {float(fy1):+.3f}) kN/m^3 (y down -> downward)")
    rng = np.random.default_rng(0)
    X = rng.uniform(-10, 25, 100000)
    Y = rng.uniform(-1, 15, 100000)
    t0 = time.time()
    for _ in range(5):
        md.force(X, Y)
    print(f"  timing: force() on 1e5 points: {(time.time() - t0) / 5 * 1e3:.1f} ms per call;"
          f" construction: ", end="")
    t0 = time.time()
    AnalyticalSeepage(45, 5.0, 5.0, 1.0)
    print(f"{time.time() - t0:.2f} s")
    print(f"\ntotal time {time.time() - t_start:.0f} s")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--quick", action="store_true", help="fewer cases")
    ap.add_argument("--sweep", action="store_true", help="also sweep beta = 15..90 by 1 deg (alpha = 1, 5, 10)")
    args = ap.parse_args()
    main(args.quick, args.sweep)
