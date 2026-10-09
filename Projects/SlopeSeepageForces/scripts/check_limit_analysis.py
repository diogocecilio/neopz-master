#!/usr/bin/env python3
"""Independent (adversarial) checks of limit_analysis.py.

Written separately from limit_analysis.py: it shares only the definition of the mechanism parameters
(theta1, theta2, s) and calls limit_analysis.py only to obtain the values being checked.

Independent work integrals: Cartesian strip quadrature.  For each vertical line x = x0 the crossings
with the boundary of the rotating region (ground surface A-O-(T)-B, log spiral B-A) are found directly
(surface: piecewise linear; spiral: bisection on its two x-monotone branches, split at theta = phi), sorted
and paired into inside intervals; each interval is split where the vertical line crosses known field
discontinuities (circles, straight lines) and integrated with graded composite Gauss rules; the x integral
uses composite Gauss panels with breakpoints at all vertices, at the x-extremes of the spiral and of the
discontinuity circles, at the intersections of the circles with the boundary, and graded towards O.
Velocity of the rigid rotation about C (paper coordinates, y down) from its definition: perpendicular to
X - C, |U| = omega |X - C|, pointing to +x directly below C.

Paper coordinates of the shared spec: O at the crest edge, x right, y DOWN.

    python3 Projects/SlopeSeepageForces/scripts/check_limit_analysis.py [--only kin,work,plot,lit,opt,ext]
"""
from __future__ import annotations

import argparse
import os
import sys
import time

import numpy as np
from numpy.polynomial.legendre import leggauss

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.dirname(SCRIPT_DIR)
DATA_DIR = os.path.join(PROJECT_DIR, "data")
RESULTS_DIR = os.path.join(PROJECT_DIR, "results", "limit_analysis")
if SCRIPT_DIR not in sys.path:
    sys.path.insert(0, SCRIPT_DIR)


# ======================================================================================================
# independent geometry
# ======================================================================================================
class Mech:
    """One log-spiral mechanism (theta1, theta2, s): s <= 1 B on the face at depth s H, s > 1 B on the toe
    ground (s - 1) H beyond the toe.  Angles at C from the -x direction turning downwards."""

    def __init__(self, t1, t2, s, beta_deg, H, phi_deg):
        self.t1, self.t2, self.s = float(t1), float(t2), float(s)
        self.H, self.beta_deg, self.phi_deg = float(H), float(beta_deg), float(phi_deg)
        b = np.radians(beta_deg)
        self.tb = np.inf if beta_deg == 90 else np.tan(b)
        self.xt = 0.0 if beta_deg == 90 else H / np.tan(b)
        self.k = np.tan(np.radians(phi_deg))
        if s <= 1.0:
            self.B = np.array([s * self.xt, s * H])
        else:
            self.B = np.array([self.xt + (s - 1.0) * H, H])
        E = np.exp((t2 - t1) * self.k)
        # A = C + r0 (-cos t1, sin t1) has y = 0; B = C + r0 E (-cos t2, sin t2)
        self.r0 = self.B[1] / (E * np.sin(t2) - np.sin(t1))
        self.C = np.array([self.B[0] + self.r0 * E * np.cos(t2), -self.r0 * np.sin(t1)])
        self.A = self.C + self.r0 * np.array([-np.cos(t1), np.sin(t1)])
        self.II = s > 1.0
        surf = [self.A, np.zeros(2)]
        if self.II:
            surf.append(np.array([self.xt, H]))
        surf.append(self.B)
        self.surf = np.array(surf)

    def spiral(self, th):
        th = np.asarray(th, float)
        r = self.r0 * np.exp((th - self.t1) * self.k)
        return self.C[0] - r * np.cos(th), self.C[1] + r * np.sin(th)

    def velocity(self, x, y):
        """rigid rotation about C, omega = 1: U perpendicular to X - C, +x below C"""
        return (y - self.C[1]), -(x - self.C[0])

    def ground_y(self, x):
        x = np.asarray(x, float)
        if self.beta_deg == 90:
            return np.where(x <= 0.0, 0.0, self.H)
        return np.where(x <= 0.0, 0.0, np.where(x < self.xt, x * self.tb, self.H))

    def admissible(self, n=4001, tol=1e-9):
        """independent admissibility: r0 > 0, A on the crest left of O, every interior spiral point strictly
        inside the soil (below the ground surface), theta1 < theta2."""
        if not (self.t2 > self.t1 and np.isfinite(self.r0) and self.r0 > 0.0 and self.A[0] < 0.0):
            return False
        if self.C[1] >= 0.0:
            return False
        th = np.linspace(self.t1, self.t2, n)[1:-1]
        x, y = self.spiral(th)
        if self.beta_deg == 90:
            margin = np.where(x < 0.0, y, np.where(x > 0.0, y - self.H, -1.0))
        else:
            margin = y - self.ground_y(x)
        return bool(np.all(margin > -tol * self.H))

    def polygon(self, n=4000):
        th = np.linspace(self.t2, self.t1, n)
        x, y = self.spiral(th)
        return np.vstack([self.surf[:-1], np.stack([x, y], 1)])

    def Pmr(self, c):
        """direct quadrature of c cos(phi) |[U]| along the spiral, dS = r dtheta / cos(phi)"""
        th, w = leggauss(200)
        th = self.t1 + 0.5 * (self.t2 - self.t1) * (th + 1)
        r = self.r0 * np.exp((th - self.t1) * self.k)
        return c * np.sum(0.5 * (self.t2 - self.t1) * w * r * r)

    # --------------------------------------------------------------------- crossings of x = x0
    def _spiral_crossings(self, x0):
        """y of the spiral at x = x0 (vector x0), up to 2 per x0 (NaN when absent)."""
        x0 = np.asarray(x0, float)
        phi = np.arctan(self.k)
        branches = [(self.t1, self.t2)]
        if self.t1 < phi < self.t2:
            branches = [(self.t1, phi), (phi, self.t2)]
        out = []
        for a, b in branches:
            xa, _ = self.spiral(a)
            xb, _ = self.spiral(b)
            lo = np.full(x0.shape, a)
            hi = np.full(x0.shape, b)
            has = (x0 - xa) * (x0 - xb) <= 0.0
            inc = xb > xa
            for _ in range(70):
                mid = 0.5 * (lo + hi)
                xm, _ = self.spiral(mid)
                go_right = (xm < x0) == inc
                lo = np.where(go_right, mid, lo)
                hi = np.where(go_right, hi, mid)
            _, y = self.spiral(0.5 * (lo + hi))
            out.append(np.where(has, y, np.nan))
        while len(out) < 2:
            out.append(np.full(x0.shape, np.nan))
        return np.stack(out, 1)

    def _surface_crossings(self, x0):
        x0 = np.asarray(x0, float)
        out = np.full(x0.shape, np.nan)
        P = self.surf
        for (xa, ya), (xb, yb) in zip(P[:-1], P[1:]):
            if xb == xa:
                continue
            t = (x0 - xa) / (xb - xa)
            sel = (t >= 0.0) & (t < 1.0)
            out = np.where(sel & np.isnan(out), ya + t * (yb - ya), out)
        return out

    def inside_intervals(self, x0):
        """(n, 2) pairs of y intervals (up to 2 per vertical line)."""
        Y = np.concatenate([self._surface_crossings(x0)[:, None], self._spiral_crossings(x0)], 1)
        Y = np.sort(Y, 1)                       # NaN last
        cnt = np.sum(np.isfinite(Y), 1)
        I1 = np.where((cnt >= 2)[:, None], Y[:, 0:2], np.nan)
        I2 = np.where((cnt >= 4)[:, None], Y[:, 2:4], np.nan)
        odd = cnt % 2 == 1
        return I1, I2, odd

    def x_extremes(self):
        th = np.linspace(self.t1, self.t2, 20001)
        x, _ = self.spiral(th)
        return x.min(), x.max()


def _graded_panels(a, b, n_uniform, n_grade, ratio=0.5, grade_a=True, grade_b=True):
    br = set(np.linspace(0.0, 1.0, n_uniform + 1).tolist())
    h = 1.0 / n_uniform
    for j in range(1, n_grade + 1):
        if grade_a:
            br.add(h * ratio ** j)
        if grade_b:
            br.add(1.0 - h * ratio ** j)
    br = np.array(sorted(br))
    return a + (b - a) * br


def strip_integral(mech: Mech, integrand, splits=None, x_extra=(), nx=60, ngx=14, qx=8, ny=4, ngy=14, qy=8,
                   grade_O=True):
    """int over the rotating region of integrand(x, y) by Cartesian strips.
    splits(x) -> list of arrays of y values where the integrand jumps along x = const;
    x_extra: additional x breakpoints."""
    xmin, xmax = mech.x_extremes()
    xmin = min(xmin, mech.A[0])
    bx = {xmin, mech.A[0], 0.0, mech.B[0], xmax}
    if mech.II:
        bx.add(mech.xt)
    bx |= {float(v) for v in x_extra if xmin < v < xmax}
    bx = np.array(sorted(v for v in bx if xmin <= v <= xmax))
    gx, wx = leggauss(qx)
    gx, wx = 0.5 * (gx + 1), 0.5 * wx
    gy, wy = leggauss(qy)
    gy, wy = 0.5 * (gy + 1), 0.5 * wy
    X, WX = [], []
    for a, b in zip(bx[:-1], bx[1:]):
        if b - a <= 1e-14 * mech.H:
            continue
        near_O_a = grade_O and abs(a) < 1e-12
        near_O_b = grade_O and abs(b) < 1e-12
        pa = _graded_panels(a, b, max(2, int(nx * (b - a) / (xmax - xmin)) + 2),
                            ngx + (12 if near_O_a or near_O_b else 0), 0.5, True, True)
        L = np.diff(pa)
        X.append((pa[:-1, None] + L[:, None] * gx).ravel())
        WX.append((L[:, None] * wx).ravel())
    X = np.concatenate(X)
    WX = np.concatenate(WX)
    I1, I2, odd = mech.inside_intervals(X)
    total = 0.0
    for I in (I1, I2):
        ok = np.isfinite(I[:, 0])
        if not np.any(ok):
            continue
        xs, ws, ya, yb = X[ok], WX[ok], I[ok, 0], I[ok, 1]
        cuts = [ya, yb]
        if splits is not None:
            for yy in splits(xs):
                cuts.append(np.clip(np.where(np.isfinite(yy), yy, ya), ya, yb))
        Cc = np.sort(np.stack(cuts, 1), 1)
        for j in range(Cc.shape[1] - 1):
            a, b = Cc[:, j], Cc[:, j + 1]
            br = _graded_panels(0.0, 1.0, ny, ngy, 0.5)
            for p0, p1 in zip(br[:-1], br[1:]):
                yy = a[:, None] + (b - a)[:, None] * (p0 + (p1 - p0) * gy)
                ww = (b - a)[:, None] * (p1 - p0) * wy
                val = integrand(np.broadcast_to(xs[:, None], yy.shape), yy)
                total += np.sum(ws[:, None] * ww * val)
    return total, int(np.sum(odd))


def work(mech, force, **kw):
    def integ(x, y):
        fx, fy = force(x, y)
        ux, uy = mech.velocity(x, y)
        return fx * ux + fy * uy
    return strip_integral(mech, integ, **kw)


# ======================================================================================================
# test fields
# ======================================================================================================
class SmoothField:
    """non-uniform, non-potential smooth field"""

    def __init__(self, H):
        self.H = H

    def force(self, x, y):
        H = self.H
        return 3.0 * np.sin(x / H) + 2.0 * y / H, 5.0 * np.cos(y / H) + x * y / H ** 2 - 1.0


class CircleJumpField:
    """piecewise smooth: different smooth fields inside / outside a circle centred at Q (not at O)"""

    def __init__(self, H, Q=(0.3, 0.4), R=0.55):
        self.H = H
        self.qx, self.qy, self.R = Q[0] * H, Q[1] * H, R * H

    def force(self, x, y):
        H = self.H
        ins = (x - self.qx) ** 2 + (y - self.qy) ** 2 < self.R ** 2
        fx = np.where(ins, 4.0 + x / H, -1.0 + 0.5 * y / H)
        fy = np.where(ins, 2.0 - y / H, 6.0 + 0.3 * x / H)
        return fx, fy

    def splits(self, x):
        d = self.R ** 2 - (x - self.qx) ** 2
        s = np.sqrt(np.maximum(d, 0.0))
        return [np.where(d > 0, self.qy - s, np.nan), np.where(d > 0, self.qy + s, np.nan)]

    def circle(self):
        return (self.qx, self.qy, self.R)

    def x_extra(self):
        return [self.qx - self.R, self.qx + self.R]


class LineJumpField:
    """discontinuous across the straight line y = y0 + a x (not a circle: limit_analysis cannot split it)"""

    def __init__(self, H, y0=0.35, a=0.4):
        self.H, self.y0, self.a = H, y0 * H, a

    def force(self, x, y):
        above = y < self.y0 + self.a * x
        return np.where(above, 3.0, -2.0), np.where(above, 1.0, 7.0)

    def splits(self, x):
        return [self.y0 + self.a * x]


def analytical_splits(field):
    """discontinuity circles of AnalyticalSeepage about O: R_w, R, R_e"""
    radii = sorted({r for r in (field.Rw, field.R, field.Re) if r > 0})

    def splits(x):
        out = []
        for R in radii:
            d = R * R - x * x
            s = np.sqrt(np.maximum(d, 0.0))
            out.append(np.where(d > 0, s, np.nan))
        return out
    return splits, radii


def circle_boundary_x(mech, qx, qy, R):
    """x of the intersections of a circle with the boundary of the region (surface + spiral)"""
    xs = []
    th = np.linspace(mech.t1, mech.t2, 20001)
    x, y = mech.spiral(th)
    g = (x - qx) ** 2 + (y - qy) ** 2 - R * R
    for i in np.nonzero(np.signbit(g[:-1]) != np.signbit(g[1:]))[0]:
        a, b = th[i], th[i + 1]
        for _ in range(60):
            m = 0.5 * (a + b)
            xm, ym = mech.spiral(m)
            gm = (xm - qx) ** 2 + (ym - qy) ** 2 - R * R
            if np.signbit(gm) == np.signbit(g[i]):
                a = m
            else:
                b = m
        xs.append(mech.spiral(0.5 * (a + b))[0])
    P = mech.surf
    for (xa, ya), (xb, yb) in zip(P[:-1], P[1:]):
        dx, dy = xb - xa, yb - ya
        A_ = dx * dx + dy * dy
        B_ = 2 * (dx * (xa - qx) + dy * (ya - qy))
        C_ = (xa - qx) ** 2 + (ya - qy) ** 2 - R * R
        disc = B_ * B_ - 4 * A_ * C_
        if disc > 0:
            for t in ((-B_ - np.sqrt(disc)) / (2 * A_), (-B_ + np.sqrt(disc)) / (2 * A_)):
                if 0 <= t <= 1:
                    xs.append(xa + t * dx)
    return xs


# ======================================================================================================
# checks
# ======================================================================================================
def sample_admissible(beta, phi, H, kind, n, rng, la):
    """random mechanisms admissible for BOTH my test and limit_analysis.Mechanisms.ok"""
    out = []
    tries = 0
    while len(out) < n and tries < 200000:
        tries += 1
        t1 = rng.uniform(0.02, 2.5)
        t2 = rng.uniform(t1 + 0.05, np.pi - 0.02)
        s = rng.uniform(0.2, 1.0) if kind == "I" else rng.uniform(1.0 + 1e-3, 3.0)
        if kind == "I" and rng.random() < 0.5:
            s = 1.0
        m = Mech(t1, t2, s, beta, H, phi)
        if not np.isfinite(m.r0) or m.r0 <= 0 or m.r0 * np.exp((t2 - t1) * m.k) > 20 * H:
            continue
        mine = m.admissible()
        theirs = bool(la.Mechanisms(t1, t2, s, beta, H, phi).ok[0])
        if mine and theirs:
            out.append(m)
    return out


def dilatancy_ok(m, n=41, tol=1e-5):
    """normality with dilatancy on the spiral: U . n_out = -|U| sin(phi), n_out pointing out of the region"""
    from matplotlib.path import Path
    th = np.linspace(m.t1, m.t2, n)[1:-1]
    x, y = m.spiral(th)
    h = 1e-6
    xp, yp = m.spiral(th + h)
    xm, ym = m.spiral(th - h)
    tx, ty = xp - xm, yp - ym
    tn = np.hypot(tx, ty)
    nx_, ny_ = ty / tn, -tx / tn
    out = ~Path(m.polygon(3000)).contains_points(np.stack([x + 1e-6 * m.H * nx_, y + 1e-6 * m.H * ny_], 1))
    nx_, ny_ = np.where(out, nx_, -nx_), np.where(out, ny_, -ny_)
    ux, uy = m.velocity(x, y)
    return bool(np.all(np.abs((ux * nx_ + uy * ny_) / np.hypot(ux, uy) + np.sin(np.radians(m.phi_deg))) < tol))


def check_kinematics(log):
    import limit_analysis as la
    log("=" * 100)
    log("K) kinematics: rotation sense, normality on the spiral, region orientation, admissibility agreement")
    log("=" * 100)
    rng = np.random.default_rng(2024)
    for beta, phi in ((30.0, 24.7), (60.0, 32.0), (90.0, 30.0), (15.0, 30.0)):
        worst_norm, worst_dir, n_tested = 0.0, 0.0, 0
        for kind in ("I", "II"):
            for m in sample_admissible(beta, phi, 1.0, kind, 20, rng, la):
                th = np.linspace(m.t1, m.t2, 401)[1:-1]
                x, y = m.spiral(th)
                h = 1e-6
                xp, yp = m.spiral(th + h)
                xm, ym = m.spiral(th - h)
                tx, ty = (xp - xm) / (2 * h), (yp - ym) / (2 * h)
                tn = np.hypot(tx, ty)
                tx, ty = tx / tn, ty / tn
                ux, uy = m.velocity(x, y)
                un = np.hypot(ux, uy)
                # normal pointing out of the rotating region = away from C side? decide by point test
                nx_, ny_ = ty, -tx
                from matplotlib.path import Path
                poly = Path(m.polygon(8000))
                eps = 1e-4 * m.H
                outside = ~poly.contains_points(np.stack([x + eps * nx_, y + eps * ny_], 1))
                nx_ = np.where(outside, nx_, -nx_)
                ny_ = np.where(outside, ny_, -ny_)
                # U . n_out = -|U| sin(phi): rotating block moves away from the rigid soil (dilatancy)
                cosang = (ux * nx_ + uy * ny_) / un
                worst_norm = max(worst_norm, np.max(np.abs(cosang + np.sin(np.radians(phi)))))
                # information: soil at A moves down (theta1 < 90 deg)
                worst_dir += m.velocity(*m.A)[1] > 0
                n_tested += 1
        log(f"   beta={beta:4.1f} phi={phi:4.1f}: {n_tested} mechanisms; max |U.n_out/|U| + sin(phi)| = "
            f"{worst_norm:.1e}; A moving down in {int(worst_dir)}")
    # admissibility agreement on random draws (also mechanisms with C on the soil side of the face line)
    log("   admissibility: mine (spiral strictly inside the soil, A on the crest) vs limit_analysis ok")
    for beta, phi in ((30.0, 24.7), (60.0, 30.0), (90.0, 30.0), (45.0, 0.0)):
        n_both = n_mine_only = n_theirs_only = 0
        mine_only = []
        for _ in range(6000):
            t1 = rng.uniform(0.0, np.pi)
            t2 = rng.uniform(t1, np.pi)
            s = rng.uniform(0.0, 1.0) if rng.random() < 0.5 else rng.uniform(1.0, 4.0)
            m = Mech(t1, t2, s, beta, 1.0, phi)
            if not np.isfinite(m.r0) or m.r0 <= 0 or m.r0 * np.exp((t2 - t1) * m.k) > 50:
                continue
            a = m.admissible(n=2001)
            b = bool(la.Mechanisms(t1, t2, s, beta, 1.0, phi).ok[0])
            n_both += a and b
            n_mine_only += a and not b
            n_theirs_only += b and not a
            if a and not b:
                mine_only.append((t1, t2, s))
        log(f"   beta={beta:4.1f} phi={phi:4.1f}: both {n_both}, mine only {n_mine_only}, limit_analysis only "
            f"{n_theirs_only}" + (f"; e.g. {np.round(mine_only[0], 4)}" if mine_only else ""))


def check_work(log):
    import limit_analysis as la
    from analytical_seepage import AnalyticalSeepage
    log("=" * 100)
    log("W) work integrals: independent Cartesian strip quadrature vs limit_analysis (3 random admissible")
    log("   mechanisms per class).  rel = limit_analysis / independent - 1")
    log("=" * 100)
    rng = np.random.default_rng(77)
    H = 5.0
    for beta, phi in ((30.0, 24.7), (60.0, 30.0), (90.0, 32.0)):
        for kind in ("I", "II"):
            mechs = sample_admissible(beta, phi, H, kind, 3, rng, la)
            for m in mechs:
                M = la.Mechanisms(m.t1, m.t2, m.s, beta, H, phi)
                line = f"   beta={beta:4.0f} phi={phi:4.1f} {kind:>2} x=({m.t1:.3f},{m.t2:.3f},{m.s:.3f}) "
                # geometry
                dgeo = max(abs(M.Cx[0] - m.C[0]), abs(-M.Cy[0] - m.C[1]), abs(M.Ax[0] - m.A[0])) / H
                pmr_rel = M.Pmr(1.0)[0] / m.Pmr(1.0) - 1
                # gamma work (gamma' = 1)
                gi, odd = work(m, lambda x, y: (np.zeros_like(x), np.ones_like(y)))
                gl = M.Pgamma(1.0)[0]
                line += f"geo {dgeo:.0e} Pmr {pmr_rel:+.0e} Pgam {gl / gi - 1:+.1e}"
                # smooth field
                sf = SmoothField(H)
                si, _ = work(m, sf.force)
                sl = la.domain_power(M, sf.force, "fine")[0]
                line += f" smooth {sl / si - 1:+.1e}"
                # circle-jump field (circle passed to limit_analysis)
                cf = CircleJumpField(H)
                ci, _ = work(m, cf.force, splits=cf.splits,
                             x_extra=cf.x_extra() + circle_boundary_x(m, *cf.circle()))
                cl = la.domain_power(M, cf.force, "fine", [cf.circle()])[0]
                cl0 = la.domain_power(M, cf.force, "fine", [])[0]
                line += f" circ {cl / ci - 1:+.1e} (no split {cl0 / ci - 1:+.1e})"
                # line-jump field (limit_analysis cannot split it)
                lf = LineJumpField(H)
                li, _ = work(m, lf.force, splits=lf.splits)
                ll = la.domain_power(M, lf.force, "fine")[0]
                llx = la.domain_power(M, lf.force, "ref")[0]
                line += f" line fine {ll / li - 1:+.1e} ref {llx / li - 1:+.1e}"
                if odd:
                    line += f" [odd crossings {odd}]"
                log(line)
    log("   analytical v'_opt field (alpha = 1), P_u: limit_analysis 'fine'/'ref' vs independent strips")
    for beta, phi, hwr in ((30.0, 24.7, 0.2), (30.0, 24.7, 1.0), (60.0, 30.0, 0.5), (60.0, 30.0, 1.0),
                           (90.0, 30.0, 1.0), (45.0, 32.0, 0.7)):
        f = AnalyticalSeepage(beta, H, hwr * H, 1.0, gamma_w=9.81)
        spl, radii = analytical_splits(f)
        for kind in ("I", "II"):
            for m in sample_admissible(beta, phi, H, kind, 3, rng, la):
                M = la.Mechanisms(m.t1, m.t2, m.s, beta, H, phi)
                xe = []
                for R in radii:
                    xe += [-R, R] + circle_boundary_x(m, 0.0, 0.0, R)
                t0 = time.time()
                ia, _ = work(m, f.force, splits=spl, x_extra=xe)
                ib, _ = work(m, f.force, splits=spl, x_extra=xe, nx=120, ngx=20, ny=6, ngy=20, qy=10)
                ti = time.time() - t0
                res = {}
                for lev in ("fine", "ref"):
                    qd = la.get_quad(lev)
                    res[lev] = la.domain_power(M, f.force, qd, la.field_circles(f, qd.n_circle_grade))[0]
                log(f"   beta={beta:4.0f} hw/H={hwr:.1f} {kind:>2} x=({m.t1:.3f},{m.t2:.3f},{m.s:.3f}): "
                    f"indep {ib:.8e} (self-conv {ia / ib - 1:+.1e}) LA fine {res['fine'] / ib - 1:+.1e} "
                    f"ref {res['ref'] / ib - 1:+.1e}  ({ti:.1f} s)")


def plot_mechanisms(log, path=os.path.join(DATA_DIR, "la_mechanisms_check.png")):
    """Draw mechanisms: random admissible ones (I, II), optima of limit_analysis, rejected ones."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import limit_analysis as la
    from analytical_seepage import AnalyticalSeepage
    rng = np.random.default_rng(3)
    panels = []
    for beta, phi, kind in ((30.0, 24.7, "I"), (60.0, 30.0, "II"), (90.0, 30.0, "I"), (45.0, 20.0, "II")):
        for m in sample_admissible(beta, phi, 1.0, kind, 1, rng, la):
            panels.append((f"random admissible {kind}\nbeta={beta:.0f} phi={phi:.1f}", m, True))
    for beta in (30.0, 60.0, 90.0):
        f = AnalyticalSeepage(beta, 5.0, 5.0, 1.0, gamma_w=9.81)
        r = la.stability_factor(beta, 5.0, 10.0, 30.0, 20.0, 9.81, f, seeds=(0,))
        m = Mech(*r["x"], beta, 1.0 * 5.0, 30.0)
        panels.append((f"optimum Fig. 9 a=1 beta={beta:.0f}\nGamma={r['Gamma']:.4f} ({r['mechanism']})", m, True))
    r = la.stability_factor(60.0, 10.0, 6.0, 32.0, 18.0, 9.81, None, seeds=(0,))
    panels.append((f"optimum dry c=6 phi=32 beta=60\nHcrit={r['Hcrit']:.3f} m", Mech(*r["x"], 60.0, 10.0, 32.0), True))
    # rejected by limit_analysis: spiral above O, A right of O, C on the soil side
    for (x, beta, phi, lab) in (((0.3, 1.2, 1.0), 60.0, 30.0, "spiral crosses the face"),
                                ((0.1868, 1.3119, 0.6008), 60.0, 30.0, "A right of O (in the air)"),
                                ((0.6504, 2.6984, 0.7204), 60.0, 30.0, "C on the soil side (theta2 > pi-beta)")):
        m = Mech(*x, beta, 1.0, phi)
        ok = bool(la.Mechanisms(*x, beta, 1.0, phi).ok[0])
        panels.append((f"{lab}\nlimit_analysis ok={ok}, mine={m.admissible()}", m, False))
    n = len(panels)
    nc = 4
    nr = (n + nc - 1) // nc
    fig, axs = plt.subplots(nr, nc, figsize=(3.6 * nc, 3.2 * nr))
    for ax, (lab, m, good) in zip(axs.ravel(), panels):
        H = m.H
        P = m.polygon(600)
        xr = max(m.B[0], m.xt) + 0.6 * H
        xl = min(P[:, 0].min(), m.C[0], 0.0) - 0.3 * H
        ax.fill(P[:, 0], P[:, 1], color="#2a78d6" if good else "#eb6834", alpha=0.25, lw=0)
        ax.plot([xl, 0, m.xt, xr], [0, 0, H, H], color="#0b0b0b", lw=1.5)
        th = np.linspace(m.t1, m.t2, 400)
        xs, ys = m.spiral(th)
        ax.plot(xs, ys, color="#2a78d6" if good else "#eb6834", lw=1.8)
        ax.plot([m.A[0], m.C[0], m.B[0]], [m.A[1], m.C[1], m.B[1]], ls=":", color="#555", lw=0.8)
        ax.plot(*m.C, "o", color="#0b0b0b", ms=3)
        # velocity arrows at a few interior points
        Pin = P[::60]
        cx = Pin[:, 0] * 0.6 + np.mean(P[:, 0]) * 0.4
        cy = Pin[:, 1] * 0.6 + np.mean(P[:, 1]) * 0.4
        ux, uy = m.velocity(cx, cy)
        sc = 0.15 * H / max(np.max(np.hypot(ux, uy)), 1e-12)
        ax.quiver(cx, cy, ux * sc, uy * sc, angles="xy", scale_units="xy", scale=1, width=0.006, color="#333")
        ax.set_aspect("equal")
        ax.invert_yaxis()
        ax.set_title(lab, fontsize=7)
        ax.tick_params(labelsize=6)
    for ax in axs.ravel()[n:]:
        ax.axis("off")
    fig.suptitle("Mechanisms (paper coordinates, y down; arrows = rigid-rotation velocity)", fontsize=9)
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)
    log(f"   mechanism drawings written to {path}")


# ======================================================================================================
# independent dry stability numbers (vectorised, no limit_analysis code)
# ======================================================================================================
def dry_N_batch(t1, t2, s, beta, phi, n=1200, smax=11.0, allow_C_below=False):
    """N = gamma H / c * Gamma (H = c = gamma = 1) for a batch of mechanisms: own geometry, admissibility by
    sampling the spiral, shoelace moments of the polygon A-O-(T)-B + spiral; inf if inadmissible / P <= 0."""
    t1, t2, s = (np.atleast_1d(np.asarray(v, float)) for v in (t1, t2, s))
    b = np.radians(beta)
    k = np.tan(np.radians(phi))
    xt = 0.0 if beta == 90 else 1.0 / np.tan(b)
    II = s > 1.0
    Bx = np.where(II, xt + (s - 1.0), s * xt)
    By = np.where(II, 1.0, s)
    E = np.exp((t2 - t1) * k)
    with np.errstate(all="ignore"):
        r0 = By / (E * np.sin(t2) - np.sin(t1))
        Cx = Bx + r0 * E * np.cos(t2)
        Cy = -r0 * np.sin(t1)                                    # y of C (y down)
        Ax = Cx - r0 * np.cos(t1)
        def spiral(u):
            th = t1[:, None] + (t2 - t1)[:, None] * u
            r = r0[:, None] * np.exp((th - t1[:, None]) * k)
            return Cx[:, None] - r * np.cos(th), Cy[:, None] + r * np.sin(th)
        # admissibility samples: uniform + clustered at both ends (small excursions into the air near A, B, T)
        ue = 10.0 ** np.linspace(-10, -1, 50)
        xa, ya = spiral(np.concatenate([np.linspace(0.0, 1.0, n)[1:-1], ue, 1.0 - ue]))
        if beta == 90:
            mg = np.where(xa < 0, ya, np.where(xa > 0, ya - 1.0, -1.0))
        else:
            mg = ya - np.where(xa <= 0, 0.0, np.where(xa < xt, xa * np.tan(b), 1.0))
        ok = (t2 > t1) & np.isfinite(r0) & (r0 > 0) & (Ax < 0) & ((Cy < 0) | allow_C_below) & (s > 0) \
            & (s <= smax) & np.all(mg > -1e-12, 1) & (r0 * E < 300)
        xs, ys = spiral(np.linspace(0.0, 1.0, n))
        M = len(t1)
        Tx, Ty = np.where(II, xt, 0.0), np.where(II, 1.0, 0.0)
        X = np.concatenate([Ax[:, None], np.zeros((M, 1)), Tx[:, None], Bx[:, None], xs[:, -2:0:-1]], 1)
        Y = np.concatenate([np.zeros((M, 2)), Ty[:, None], By[:, None], ys[:, -2:0:-1]], 1)
        X1, Y1 = np.roll(X, -1, 1), np.roll(Y, -1, 1)
        cr = X * Y1 - X1 * Y
        A = 0.5 * cr.sum(1)
        Mx = ((X + X1) * cr).sum(1) / 6.0
        sg = np.sign(A)
        Pg = Cx * A * sg - Mx * sg                                # int (Cx - x) dA
        fac = (t2 - t1) if k == 0 else np.expm1(2 * (t2 - t1) * k) / (2 * k)
        return np.where(ok & (Pg > 0), r0 ** 2 * fac / Pg, np.inf)


def dry_N_solve(beta, phi, nrand=300000, seed=0):
    from scipy.optimize import minimize
    rng = np.random.default_rng(seed)
    cand = []
    for _ in range(nrand // 20000):
        t1, t2 = rng.uniform(0, np.pi, 20000), rng.uniform(0, np.pi, 20000)
        s = np.where(rng.random(20000) < 0.3, 1.0,
                     np.where(rng.random(20000) < 0.5, rng.uniform(0, 1, 20000), rng.uniform(1, 11, 20000)))
        N = dry_N_batch(t1, t2, s, beta, phi, n=300)
        cand += [(N[i], (t1[i], t2[i], s[i])) for i in np.argsort(N)[:3] if np.isfinite(N[i])]
    cand.sort(key=lambda v: v[0])
    best = (np.inf, None)
    for _, z0 in cand[:4]:
        def f(z):
            return float(dry_N_batch(*z, beta, phi, n=4000)[0])
        r = minimize(f, z0, method="Nelder-Mead", options=dict(xatol=1e-10, fatol=1e-12, maxiter=3000))
        r = minimize(f, r.x, method="Nelder-Mead", options=dict(xatol=1e-11, fatol=1e-13, maxiter=3000))
        if r.fun < best[0]:
            best = (r.fun, r.x)
    if best[1] is None:
        return np.inf, None
    N4, N8 = dry_N_batch(*best[1], beta, phi, n=4000)[0], dry_N_batch(*best[1], beta, phi, n=8000)[0]
    return N8 + (N8 - N4) / 3.0, best[1]          # Richardson in the number of spiral chords


def check_literature(log):
    import limit_analysis as la
    log("=" * 100)
    log("L) dry stability numbers N = gamma H_c / c: independent (own geometry + shoelace + random search + NM)")
    log("   vs limit_analysis; Fig. 8 h_w = 0 ends from N: H_crit = N c / (gamma - gamma_w)")
    log("=" * 100)
    for beta, phi in ((90, 0), (90, 20), (90, 40), (75, 0), (75, 30), (60, 0), (60, 15), (60, 30), (45, 20),
                      (45, 40), (30, 25), (45, 0)):
        t0 = time.time()
        N, x = dry_N_solve(beta, phi)
        r = la.stability_factor(beta, 1.0, 1.0, phi, 1.0, 0.0, None, seeds=(0,))
        log(f"   beta={beta:3d} phi={phi:3d}: N indep {N:10.5f} at {np.round(x, 4)}; limit_analysis {r['Gamma']:10.5f} "
            f"rel {r['Gamma'] / N - 1:+.1e} ({time.time() - t0:.0f} s)")
    paper = {("London", 60): 17.949, ("London", 30): 156.555, ("Israeli", 60): 12.986, ("Israeli", 35): 229.319}
    soils = {"London (Table 1)": (6.0, 32.0), "Israeli (Table 1)": (11.7, 24.7)}
    for (panel, beta), pap in paper.items():
        row = []
        for name, (c, phi) in soils.items():
            N, _ = dry_N_solve(float(beta), phi, nrand=200000)
            row.append(f"{name}: {N * c / (18.0 - 9.8):9.4f} (gw 9.8) {N * c / (18.0 - 9.81):9.4f} (gw 9.81)")
        log(f"   Fig. 8 panel {panel:7s} beta={beta}: paper {pap:8.3f} | " + " | ".join(row))


def check_optimizer(log):
    """differential evolution (scipy, other seeds, init from admissible random points) + the same NM polish,
    on the limit_analysis objective, vs stability_factor (PSO + NM)"""
    import limit_analysis as la
    from analytical_seepage import AnalyticalSeepage
    from scipy.optimize import differential_evolution
    log("=" * 100)
    log("O) optimiser: limit_analysis.stability_factor vs scipy differential_evolution (+ same NM polish)")
    log("=" * 100)
    cases = [("Fig9 a=1 beta=30", 30.0, 5.0, 10.0, 30.0, 20.0, 9.81, 5.0),
             ("Fig9 a=1 beta=90", 90.0, 5.0, 10.0, 30.0, 20.0, 9.81, 5.0),
             ("Israeli beta=35 hw=0.5H", 35.0, 10.0, 11.7, 24.7, 18.0, 9.81, 5.0),
             ("London(Table 1) beta=30 hw=0.2H", 30.0, 10.0, 6.0, 32.0, 18.0, 9.81, 2.0),
             ("London(Table 1) beta=60 hw=0", 60.0, 10.0, 6.0, 32.0, 18.0, 9.81, 0.0)]
    for name, beta, H, c, phi, g, gw, hw in cases:
        f = AnalyticalSeepage(beta, H, hw, 1.0, gamma_w=gw) if hw > 0 else None
        t0 = time.time()
        r = la.stability_factor(beta, H, c, phi, g, gw, f)
        prob = la.Problem(beta, H, c, phi, g, gw, f)
        best = np.inf
        out = []
        for kind in ("I", "II"):
            lb, ub = prob.bounds(kind)
            init = la._random_admissible(prob, kind, 60, np.random.default_rng(21))
            if len(init) < 5:
                continue
            res = differential_evolution(lambda X, kind=kind: prob.gamma_factor(np.asarray(X).T, kind, "coarse"),
                                         list(zip(lb, ub)), init=init, maxiter=300, tol=1e-10, seed=21,
                                         polish=False, vectorized=True, updating="deferred", mutation=(0.5, 1.0),
                                         recombination=0.9)
            xp, fp = la.nelder_mead(lambda z, kind=kind: prob.gamma_factor(z, kind, "fine")[0], res.x, lb, ub)
            out.append(f"{kind}: DE {res.fun:.6g} -> NM {fp:.8f}")
            best = min(best, fp)
        log(f"   {name}: stability_factor {r['Gamma']:.8f} ({r['mechanism']}, x={np.round(r['x'], 5)}); "
            f"DE+NM best {best:.8f} rel {best / r['Gamma'] - 1:+.1e}  [" + "; ".join(out) + f"] ({time.time() - t0:.0f} s)")


def check_extended_class(log):
    """mechanisms rejected by limit_analysis but geometrically valid (C on the soil side of the face line,
    C below the crest): random search with my dry evaluator, then NM in the union of classes"""
    import limit_analysis as la
    from scipy.optimize import minimize
    log("=" * 100)
    log("E) mechanisms outside the limit_analysis class (C on the soil side of the face, C below the crest),")
    log("   dry, H = c = gamma = 1: random search, NM from the best (may walk back into the class)")
    log("=" * 100)
    rng = np.random.default_rng(8)
    for beta, phi in ((60.0, 32.0), (90.0, 0.0), (30.0, 20.0), (60.0, 24.7)):
        n = 200000
        t1 = rng.uniform(-0.8, np.pi, n)
        t2 = rng.uniform(0.0, np.pi, n)
        s = np.where(rng.random(n) < 0.5, rng.uniform(0, 1, n), rng.uniform(1, 6, n))
        N = np.concatenate([dry_N_batch(t1[i:i + 10000], t2[i:i + 10000], s[i:i + 10000], beta, phi, n=400,
                                        allow_C_below=True) for i in range(0, n, 10000)])
        okLA = la.Mechanisms(t1, t2, s, beta, 1.0, phi).ok
        ext = np.isfinite(N) & ~okLA
        r = la.stability_factor(beta, 1.0, 1.0, phi, 1.0, 0.0, None, seeds=(0,))
        if not np.any(ext):
            log(f"   beta={beta:4.1f} phi={phi:4.1f}: limit_analysis {r['Gamma']:.6f}; no extended mechanism with P > 0")
            continue
        j = np.argmin(np.where(ext, N, np.inf))
        z = (t1[j], t2[j], s[j])
        res = minimize(lambda q: float(dry_N_batch(*q, beta, phi, n=4000, allow_C_below=True)[0]), z,
                       method="Nelder-Mead", options=dict(xatol=1e-9, fatol=1e-11, maxiter=3000))
        inLA = bool(la.Mechanisms(*res.x, beta, 1.0, phi).ok[0])
        dil = dilatancy_ok(Mech(*z, beta, 1.0, phi))
        log(f"   beta={beta:4.1f} phi={phi:4.1f}: limit_analysis {r['Gamma']:.6f}; {int(ext.sum())} extended mechanisms "
            f"with P > 0, best {N[j]:.5f} at {np.round(z, 4)} -> NM {res.fun:.6f} at {np.round(res.x, 4)} "
            f"(inside the limit_analysis class: {inLA}; random best dilatant: {dil})")


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default="kin,work,plot,lit,opt,ext")
    a = ap.parse_args(argv)
    sel = a.only.split(",")
    os.makedirs(RESULTS_DIR, exist_ok=True)
    fh = open(os.path.join(RESULTS_DIR, "independent_checks_output.txt"), "w" if a.only == ap.get_default("only") else "a")

    def log(s=""):
        print(s, flush=True)
        fh.write(s + "\n")
        fh.flush()
    log(f"# check_limit_analysis.py --only {a.only}  ({time.strftime('%Y-%m-%d %H:%M:%S')})")
    if "kin" in sel:
        check_kinematics(log)
    if "work" in sel:
        check_work(log)
    if "plot" in sel:
        plot_mechanisms(log)
    if "lit" in sel:
        check_literature(log)
    if "opt" in sel:
        check_optimizer(log)
    if "ext" in sel:
        check_extended_class(log)
    fh.close()


if __name__ == "__main__":
    main()
