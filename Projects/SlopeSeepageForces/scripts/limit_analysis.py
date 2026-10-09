#!/usr/bin/env python3
"""Kinematic limit analysis of Ceron et al. (IJNAMG 2025), Sect. 4.2 (Figs. 6-7, Eqs. 42-58).

Rotational log-spiral failure mechanisms about a centre C (paper coordinates of the shared spec:
origin O at the crest edge, x to the right, y DOWNWARD, toe T = (H / tan(beta), H)).

Angles theta are measured at C from the direction pointing to -x, turning downwards (Fig. 7):
    point(rho, theta) = C + rho * e_r,   e_r = (-cos theta, sin theta),   e_theta = (sin theta, cos theta)
    C = (Cx, -Cy) in (x, y-down) coordinates, i.e. Cy (Eq. 45) is the HEIGHT of C above the crest.
The log spiral (Eq. 42)  r(theta) = r0 exp((theta - theta1) tan phi),  theta1 <= theta <= theta2,
goes from A = C + r0 e_r(theta1) on the crest to B = C + r_h e_r(theta2) on the ground surface.
Velocity field: U = omega r e_theta = omega (y + Cy, Cx - x)  (rigid rotation, omega = 1 below).

Mechanism I : B on the face,  B = (eta H / tan beta, eta H), 0 < eta <= 1  (Fig. 7, parameters theta1, theta2, eta)
Mechanism II: B on the toe ground beyond the toe, B = (H / tan beta + d, H), d >= 0  (Fig. 6 right)
Unified parameter s: s = eta in (0, 1] (I),  s = 1 + d / H > 1 (II).  Given (theta1, theta2, s), B is known and
    r0 = B_y / (E sin theta2 - sin theta1),  E = exp((theta2 - theta1) tan phi)      (Eq. 43 for I)
    Cx = B_x + r0 E cos theta2,  Cy = r0 sin theta1,  L = OA = r0 cos theta1 - Cx      (Eqs. 44-46)

Admissibility (Eq. 56 plus the geometric conditions it implies):
    0 < theta1 < theta2 < pi, r0 > 0, L > 0 (A on the crest left of O),
    C on the air side of the face line (Cx sin beta + Cy cos beta > 0; for I this is theta2 < pi - beta),
    so that the ground surface A-O-(T)-B is seen from C with a monotonically increasing angle,
    spiral below the ground surface between A and B, size r_h <= 1000 H (guards the closed forms against
    round-off for degenerate, nearly flat huge arcs).  Along every straight piece of the surface the
    function ln(r(theta) / rho_surface(theta)) is concave, so it suffices to check the spiral below O
    (and below T for mechanism II).   P_gamma + P_u > 0.

Rates of work (per unit omega):
    P_mr   = c r0^2 (exp(2 (theta2 - theta1) tan phi) - 1) / (2 tan phi)   (Eq. 48; -> c r0^2 (theta2 - theta1), phi -> 0)
    P_gamma = gamma' (f1 - f2 - f3 [- f4 for II])                           (Eqs. 50-53, f4: triangle C-T-B)
    P_u    = int_Omega f . U dOmega,  f = field.force(x, y)                  (Eq. 49)
P_gamma and P_u are ALSO computed by quadrature in polar coordinates about C: for each ray theta, rho from the
ground surface to the spiral; the angular range is split at the rays through O and T, with composite
Gauss-Legendre panels graded towards O and T, and the rays are split where they cross the discontinuity
circles of the field (analytical field: r = R_w, R, R_e about O, plus circles graded towards the r^(m-1)
singularity at O).  For a potential field f = -grad u (FE field) P_u can also be computed exactly by the
boundary formula  P_u = -int_dOmega u U.n dS  (div U = 0): the u values on the ground surface and along the
spiral only (pu_method="boundary").

Stability factor Gamma = min P_mr / (P_gamma + P_u) (Eq. 55) by Particle Swarm Optimisation (fixed seeds) for
each mechanism class, followed by a Nelder-Mead polish with a finer quadrature.  H_crit = Gamma * H (spec).

Run as a script for the verification checks and the comparisons with Figs. 8 (h_w = 0) and 9:
    python3 Projects/SlopeSeepageForces/scripts/limit_analysis.py [--quick] [--only a,b,...]
"""
from __future__ import annotations

import argparse
import os
import sys
import time

import numpy as np
from scipy.optimize import minimize

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.dirname(SCRIPT_DIR)
DATA_DIR = os.path.join(PROJECT_DIR, "data")
RESULTS_DIR = os.path.join(PROJECT_DIR, "results", "limit_analysis")
if SCRIPT_DIR not in sys.path:
    sys.path.insert(0, SCRIPT_DIR)

GAMMA_W = 9.81
BIG = 1e30          # objective value of inadmissible mechanisms (finite: keeps Nelder-Mead arithmetic clean)
_EPS_ANG = 1e-6


# ======================================================================================================
# fields used by the checks (same interface as the seepage fields of the spec)
# ======================================================================================================
class NoSeepage:
    """f = 0 everywhere (h_w = 0)."""

    def force(self, x, y):
        x, y = np.broadcast_arrays(np.asarray(x, float), np.asarray(y, float))
        return np.zeros(x.shape), np.zeros(x.shape)

    def u(self, x, y):
        x, y = np.broadcast_arrays(np.asarray(x, float), np.asarray(y, float))
        return np.zeros(x.shape)


class UniformForce:
    """f = (gx, gy) everywhere (e.g. an additional unit weight g: f = (0, g))."""

    force_is_minus_grad_u = True

    def __init__(self, gx=0.0, gy=0.0):
        self.gx, self.gy = float(gx), float(gy)

    def force(self, x, y):
        x, y = np.broadcast_arrays(np.asarray(x, float), np.asarray(y, float))
        return np.full(x.shape, self.gx), np.full(x.shape, self.gy)

    def u(self, x, y):
        return -(self.gx * np.asarray(x, float) + self.gy * np.asarray(y, float))


class PolynomialPotential:
    """Smooth test potential u = a x^2 y - b y^3 + c x y + d x, f = -grad u (domain vs boundary formula)."""

    force_is_minus_grad_u = True

    def __init__(self, a=0.7, b=0.3, c=1.1, d=-0.4):
        self.a, self.b, self.c, self.d = a, b, c, d

    def u(self, x, y):
        x, y = np.asarray(x, float), np.asarray(y, float)
        return self.a * x * x * y - self.b * y ** 3 + self.c * x * y + self.d * x

    def force(self, x, y):
        x, y = np.broadcast_arrays(np.asarray(x, float), np.asarray(y, float))
        return -(2 * self.a * x * y + self.c * y + self.d), -(self.a * x * x - 3 * self.b * y * y + self.c * x)


def _is_zero_field(field):
    return field is None or type(field).__name__ == "NoSeepage"


def _is_potential_field(field):
    """True when field.force == -grad field.u exactly (FE P1/P2 fields, test potentials)."""
    return callable(getattr(field, "u", None)) and (
        getattr(field, "force_is_minus_grad_u", False) or type(field).__name__ == "FESeepage")


def field_circles(field, n_grade=0, grade_ratio=0.3):
    """Circles (cx, cy, R, is_discontinuity) of a field: every circle splits the radial quadrature; the
    discontinuity circles also add angular panel breaks (rays tangent to the circle, rays through its
    intersections with the ground surface and with the spiral).  field.quadrature_circles() if provided;
    AnalyticalSeepage: r = R_w, R, R_e about O (discontinuities) plus n_grade circles R_w * grade_ratio^k
    towards the singular point O (r^(m-1) field)."""
    if field is None:
        return []
    if callable(getattr(field, "quadrature_circles", None)):
        return list(field.quadrature_circles())
    if all(hasattr(field, a) for a in ("Rw", "R", "Re")):
        out = [(0.0, 0.0, r, True) for r in sorted(set((field.Rw, field.R, field.Re))) if r > 0.0]
        if field.Rw > 0.0:
            out += [(0.0, 0.0, field.Rw * grade_ratio ** k, False) for k in range(1, n_grade + 1)]
        return out
    return []


# ======================================================================================================
# geometry of a batch of mechanisms
# ======================================================================================================
class Mechanisms:
    """Geometry of a batch of log-spiral mechanisms given (theta1, theta2, s) arrays (see module doc)."""

    R_MAX = 1000.0   # admissible mechanisms: r_h <= R_MAX * H (beyond, the closed forms lose all digits)

    def __init__(self, theta1, theta2, s, beta_deg, H, phi_deg):
        t1, t2, s = (np.atleast_1d(np.asarray(v, float)) for v in (theta1, theta2, s))
        t1, t2, s = np.broadcast_arrays(t1, t2, s)
        self.t1, self.t2, self.s = t1.copy(), t2.copy(), s.copy()
        self.beta_deg, self.H, self.phi_deg = float(beta_deg), float(H), float(phi_deg)
        b = np.radians(self.beta_deg)
        self.beta = b
        self.sb, self.cb = np.sin(b), (0.0 if abs(self.beta_deg - 90.0) < 1e-12 else np.cos(b))
        self.xtoe = H * self.cb / self.sb
        self.k = np.tan(np.radians(self.phi_deg))
        k, xt = self.k, self.xtoe
        self.II = s > 1.0
        self.Bx = np.where(self.II, xt + (s - 1.0) * H, s * xt)
        self.By = np.where(self.II, H, s * H)
        with np.errstate(all="ignore"):
            self.E = np.exp((t2 - t1) * k)
            self.den = self.E * np.sin(t2) - np.sin(t1)
            self.r0 = self.By / self.den
            self.rh = self.r0 * self.E
            self.Cx = self.Bx + self.rh * np.cos(t2)
            self.Cy = self.r0 * np.sin(t1)
            self.Ax = self.Cx - self.r0 * np.cos(t1)
            self.L = -self.Ax
            self.thO = np.arctan2(self.Cy, self.Cx)                     # Eq. (47)
            self.dO = np.hypot(self.Cx, self.Cy)
            self.thT = np.arctan2(self.Cy + H, self.Cx - xt)
            self.dT = np.hypot(self.Cx - xt, self.Cy + H)
            self.D = self.Cx * self.sb + self.Cy * self.cb             # distance from C to the face line
            self.thF = np.where(self.II, self.thT, t2)                  # end of the face piece
            # admissibility
            ok = (t1 > 0.0) & (t2 > t1) & (t2 < np.pi) & (s > 0.0)
            ok &= self.den > 1e-12 * np.maximum(self.E * np.abs(np.sin(t2)), 1e-300)     # no round-off r0
            ok &= np.isfinite(self.r0) & (self.r0 > 0.0) & (self.rh <= self.R_MAX * H)
            ok &= (self.L > 0.0) & (self.D > 0.0)
            ok &= (self.thO >= t1) & (self.thO <= self.thF)
            ok &= np.where(self.II, True, t2 < np.pi - b + 1e-12)                        # Eq. (56)
            lr0 = np.log(self.r0)
            ok &= (self.thO - t1) * k + lr0 >= np.log(self.dO) - 1e-12               # spiral below O
            ok &= np.where(self.II, ((self.thT - t1) * k + lr0 >= np.log(self.dT) - 1e-12)
                           & (self.thT <= t2), True)                                  # spiral below T
        self.ok = ok

    def __len__(self):
        return len(self.t1)

    def subset(self, idx):
        m = Mechanisms(self.t1[idx], self.t2[idx], self.s[idx], self.beta_deg, self.H, self.phi_deg)
        return m

    def spiral_r(self, theta):
        return self.r0 * np.exp((theta - self.t1) * self.k)

    # ------------------------------------------------------------------ closed forms (per unit omega)
    def Pmr(self, c):
        """Eq. (48), divided by omega; phi -> 0 limit c r0^2 (theta2 - theta1)."""
        dth = self.t2 - self.t1
        k = self.k
        fac = dth if k < 1e-14 else np.expm1(2.0 * dth * k) / (2.0 * k)
        return c * self.r0 ** 2 * fac

    def f_terms(self):
        """f1, f2, f3 (Eqs. 51-53, f3 in a form valid for beta = 90 deg) and f4 (triangle C-T-B, mechanism II)."""
        k, t1, t2, r0, rh = self.k, self.t1, self.t2, self.r0, self.rh
        f1 = (rh ** 3 * (3 * k * np.cos(t2) + np.sin(t2)) - r0 ** 3 * (3 * k * np.cos(t1) + np.sin(t1))) \
            / (3.0 * (9.0 * k * k + 1.0))
        f2 = r0 ** 3 * np.sin(t1) / 6.0 * (1.0 - np.sin(t1) ** 2 / np.sin(self.thO) ** 2)
        cb, sb = self.cb, self.sb

        def G(psi):
            return -cb / (2.0 * np.sin(psi) ** 2) - sb * np.cos(psi) / np.sin(psi)

        f3 = self.D ** 3 / 3.0 * (G(self.thF + self.beta) - G(self.thO + self.beta))
        with np.errstate(all="ignore"):
            f4 = np.where(self.II, (self.H + self.Cy) ** 3 / 6.0
                          * (1.0 / np.sin(self.thT) ** 2 - 1.0 / np.sin(t2) ** 2), 0.0)
        return f1, f2, f3, f4

    def f3_paper(self):
        """Eq. (53) literally (mechanism I, beta < 90 deg)."""
        tb = np.tan(self.beta)
        return (self.Cy + tb * self.Cx) ** 3 / 6.0 * (1.0 / (np.tan(self.thO) + tb) ** 2
                                                      - 1.0 / (np.tan(self.t2) + tb) ** 2)

    def Pgamma(self, gamma_p):
        """Eq. (50) (+ f4 for mechanism II), divided by omega."""
        f1, f2, f3, f4 = self.f_terms()
        return gamma_p * (f1 - f2 - f3 - f4)

    # ------------------------------------------------------------------ points
    def A(self):
        return np.stack([self.Ax, np.zeros_like(self.Ax)], -1)

    def B(self):
        return np.stack([self.Bx, self.By], -1)

    def C(self):
        """centre in paper coordinates (y down): (Cx, -Cy)"""
        return np.stack([self.Cx, -self.Cy], -1)

    def kind(self):
        return np.where(self.II, "II", "I")

    def pieces(self):
        """Angular pieces (M, 3): crest [t1, thO], face [thO, thF], toe ground [thT, t2] (empty for I)."""
        ta = np.stack([self.t1, self.thO, np.where(self.II, self.thT, self.t2)], -1)
        tb = np.stack([self.thO, self.thF, self.t2], -1)
        return ta, np.maximum(tb, ta)

    def rho_surface(self, piece, theta):
        """distance from C to the ground surface along the ray theta, for piece 0 (crest), 1 (face), 2 (toe)"""
        sl = (slice(None),) + (None,) * (theta.ndim - 1)
        if piece == 0:
            return self.Cy[sl] / np.sin(theta)
        if piece == 1:
            return self.D[sl] / np.sin(theta + self.beta)
        return (self.H + self.Cy[sl]) / np.sin(theta)


# ======================================================================================================
# quadrature
# ======================================================================================================
def _gauss01(q):
    x, w = np.polynomial.legendre.leggauss(q)
    return 0.5 * (x + 1.0), 0.5 * w


class Quadrature:
    """Composite Gauss rule in polar coordinates about C.

    theta: n_panels uniform panels per piece + n_grade geometric panels (ratio grade_ratio) at both ends
           of each piece (rays through A, O, T, B), q_theta points per panel;
    rho  : n_rho uniform panels between the ground surface and the spiral, further split at the
           intersections with the field circles, q_rho points per sub-segment."""

    def __init__(self, n_panels=6, q_theta=6, n_grade=2, grade_ratio=0.25, n_rho=2, q_rho=5,
                 n_circle_grade=0, name=""):
        self.name = name
        self.n_panels, self.q_theta, self.n_grade, self.grade_ratio = n_panels, q_theta, n_grade, grade_ratio
        self.n_rho, self.q_rho, self.n_circle_grade = n_rho, q_rho, n_circle_grade
        h = 1.0 / n_panels
        br = set(np.linspace(0.0, 1.0, n_panels + 1).tolist())
        for j in range(1, n_grade + 1):
            br.add(h * grade_ratio ** j)
            br.add(1.0 - h * grade_ratio ** j)
        br = np.array(sorted(br))
        self.br = br
        self.x, self.w = _gauss01(q_theta)
        self.t = (br[:-1, None] + np.diff(br)[:, None] * self.x).ravel()
        self.wt = (np.diff(br)[:, None] * self.w).ravel()
        self.u, self.wu = _gauss01(q_rho)

    def __repr__(self):
        return (f"Quadrature({self.name}: panels={self.n_panels}x{self.q_theta} grade={self.n_grade} "
                f"rho={self.n_rho}x{self.q_rho} circle_grade={self.n_circle_grade}, {len(self.t)} rays/piece)")


QUAD_LEVELS = {
    "coarse": dict(n_panels=4, q_theta=5, n_grade=1, n_rho=1, q_rho=4, n_circle_grade=0),
    "medium": dict(n_panels=8, q_theta=6, n_grade=3, n_rho=2, q_rho=5, n_circle_grade=2),
    "fine": dict(n_panels=12, q_theta=8, n_grade=5, n_rho=2, q_rho=6, n_circle_grade=4),
    "xfine": dict(n_panels=24, q_theta=10, n_grade=8, n_rho=4, q_rho=8, n_circle_grade=7),
    "ref": dict(n_panels=48, q_theta=12, n_grade=12, n_rho=6, q_rho=10, n_circle_grade=10),
}


def get_quad(q):
    if isinstance(q, Quadrature):
        return q
    return Quadrature(name=q, **QUAD_LEVELS[q])


def domain_power(mech: Mechanisms, force, quad="fine", circles=(), chunk=64):
    """int_Omega f . U dOmega / omega for every mechanism of the batch (assumed admissible).
    force(x, y) -> (fx, fy) (vectorised); circles: (cx, cy, R[, is_discontinuity]) where f may jump
    (3-tuples are discontinuities).  Jumps across curves that are NOT declared circles (e.g. element edges
    of a FE gradient) are not split: the error then decreases only algebraically (check_limit_analysis.py:
    a straight-line jump gives 1e-3..3e-2 at 'fine', 1e-5..2e-3 at 'ref'); use the boundary formula for
    potential fields (pu_method='auto' does so for FESeepage)."""
    quad = get_quad(quad)
    circles = [tuple(c) + (True,) if len(c) == 3 else tuple(c) for c in circles]
    M = len(mech)
    out = np.zeros(M)
    for i0 in range(0, M, chunk):
        idx = np.arange(i0, min(M, i0 + chunk))
        out[idx] = _domain_power_chunk(mech, idx, force, quad, circles)
    return out


def _circle_angles(mech, idx, circle, n_sample=48, n_bisect=48):
    """Angles (seen from C) where the integrand of the angular integral has kinks because of the circle
    (qx, qy, R): rays tangent to the circle, rays through its intersections with the ground surface and
    with the spiral (up to 4).  Returns (m, 12) with NaN where absent."""
    qx, qy, R = circle[:3]
    Cx, Cy = mech.Cx[idx], mech.Cy[idx]
    m = len(idx)
    H, xt, sb, cb = mech.H, mech.xtoe, mech.sb, mech.cb
    out = []
    # tangent rays
    dx, dy = Cx - qx, Cy + qy
    dist = np.hypot(dx, dy)
    thc = np.arctan2(dy, dx)
    with np.errstate(invalid="ignore"):
        a = np.arcsin(np.where(dist > R, R / dist, np.nan))
    out += [thc - a, thc + a]
    # intersections with the ground surface (fixed points, independent of the mechanism)
    pts = []
    if R > abs(qy):
        h = np.sqrt(R * R - qy * qy)
        pts += [(x0, 0.0) for x0 in (qx - h, qx + h) if x0 <= 0.0]
    dq = qx * cb + qy * sb
    disc = dq * dq - (qx * qx + qy * qy - R * R)
    if disc > 0:
        for t in (dq - np.sqrt(disc), dq + np.sqrt(disc)):
            if 0.0 <= t <= H / sb:
                pts.append((t * cb, t * sb))
    if R > abs(H - qy):
        h = np.sqrt(R * R - (H - qy) ** 2)
        pts += [(x0, H) for x0 in (qx - h, qx + h) if x0 >= xt]
    for (px, py) in pts[:6]:
        out.append(np.arctan2(py + Cy, Cx - px))
    out += [np.full(m, np.nan)] * (6 - min(len(pts), 6))
    # intersections with the spiral: sign changes of g(theta) = |S(theta) - Q|^2 - R^2, then bisection
    t1, t2, r0, k = mech.t1[idx], mech.t2[idx], mech.r0[idx], mech.k

    def g(th):
        r = r0[:, None] * np.exp((th - t1[:, None]) * k)
        return (Cx[:, None] - r * np.cos(th) - qx) ** 2 + (-Cy[:, None] + r * np.sin(th) - qy) ** 2 - R * R

    ts = t1[:, None] + (t2 - t1)[:, None] * np.linspace(0.0, 1.0, n_sample)
    gs = g(ts)
    ch = np.signbit(gs[:, :-1]) != np.signbit(gs[:, 1:])
    roots = np.full((m, 4), np.nan)
    order = np.argsort(~ch, axis=1, kind="stable")[:, :4]       # first (up to) 4 sign changes
    for j in range(4):
        jj = order[:, j]
        has = ch[np.arange(m), jj]
        if not np.any(has):
            continue
        lo = ts[np.arange(m), jj].copy()
        hi = ts[np.arange(m), jj + 1].copy()
        glo = gs[np.arange(m), jj]
        for _ in range(n_bisect):
            mid = 0.5 * (lo + hi)
            gm = g(mid[:, None])[:, 0]
            left = np.signbit(gm) == np.signbit(glo)
            lo = np.where(left, mid, lo)
            glo = np.where(left, gm, glo)
            hi = np.where(left, hi, mid)
        roots[:, j] = np.where(has, 0.5 * (lo + hi), np.nan)
    out += [roots[:, j] for j in range(4)]
    return np.stack(out, 1)


def _domain_power_chunk(mech, idx, force, quad, circles):
    m = len(idx)
    ta, tb = mech.pieces()
    ta, tb = ta[idx], tb[idx]
    t1, t2 = mech.t1[idx], mech.t2[idx]
    # angular breakpoints: per-piece uniform + graded panels, plus the kinks due to discontinuity circles
    bl = [(ta[:, p, None] + (tb - ta)[:, p, None] * quad.br) for p in range(3)]
    for c in circles:
        if c[3]:
            sa = _circle_angles(mech, idx, c)
            bl.append(np.where(np.isfinite(sa), np.clip(sa, t1[:, None], t2[:, None]), t1[:, None]))
    BRK = np.sort(np.concatenate(bl, 1), 1)                               # (m, Nb)
    PA, PL = BRK[:, :-1], np.diff(BRK, axis=1)
    TH = (PA[:, :, None] + PL[:, :, None] * quad.x).reshape(m, -1)        # (m, Nt)
    WT = (PL[:, :, None] * quad.w).reshape(m, -1)
    Cx, Cy, r0 = mech.Cx[idx], mech.Cy[idx], mech.r0[idx]
    thO, thF = mech.thO[idx][:, None], mech.thF[idx][:, None]
    with np.errstate(all="ignore"):
        RS = np.where(TH < thO, Cy[:, None] / np.sin(TH),
                      np.where(TH < thF, mech.D[idx][:, None] / np.sin(TH + mech.beta),
                               (mech.H + Cy[:, None]) / np.sin(TH)))
    RS = np.where(WT > 0, RS, 0.0)
    RP = r0[:, None] * np.exp((TH - t1[:, None]) * mech.k)
    RP = np.where(WT > 0, np.maximum(RP, RS), 0.0)
    cT, sT = np.cos(TH), np.sin(TH)
    brk = [RS + (RP - RS) * (j / quad.n_rho) for j in range(quad.n_rho + 1)]
    for c in circles:
        qx, qy, R = c[:3]
        dx = (Cx - qx)[:, None]
        dy = (Cy + qy)[:, None]
        bb = dx * cT + dy * sT
        disc = bb * bb - (dx * dx + dy * dy - R * R)
        sq = np.sqrt(np.maximum(disc, 0.0))
        for root in (bb - sq, bb + sq):
            brk.append(np.where(disc > 0, np.clip(root, RS, RP), RS))
    BR = np.sort(np.stack(brk, -1), -1)                                  # (m, Nt, S+1)
    SA = BR[..., :-1]
    SL = np.diff(BR, axis=-1)
    RHO = SA[..., None] + SL[..., None] * quad.u                         # (m, Nt, S, Q)
    W = (WT[..., None] * SL)[..., None] * quad.wu
    mask = W > 0.0
    mid = np.broadcast_to(np.arange(m)[:, None, None, None], W.shape)[mask]
    rho = RHO[mask]
    cc = np.broadcast_to(cT[..., None, None], W.shape)[mask]
    ss = np.broadcast_to(sT[..., None, None], W.shape)[mask]
    x = Cx[mid] - rho * cc
    y = -Cy[mid] + rho * ss
    fx, fy = force(x, y)
    val = (np.asarray(fx) * ss + np.asarray(fy) * cc) * rho * rho * W[mask]
    return np.bincount(mid, weights=val, minlength=m)


def boundary_power(mech: Mechanisms, u, hw=None, n_panels=24, q=8):
    """P_u / omega = -int_dOmega u U.n dS for a potential field f = -grad u (div U = 0).
    Spiral part: tan(phi) int u r^2 dtheta;  ground surface part: crest A-O, face O-F, toe ground T-B.
    hw (optional): water line depth on the face, used as a panel break (kink of the surface data)."""
    M = len(mech)
    x01, w01 = _gauss01(q)
    br = np.linspace(0.0, 1.0, n_panels + 1)
    tq = (br[:-1, None] + np.diff(br)[:, None] * x01).ravel()
    wq = (np.diff(br)[:, None] * w01).ravel()
    eps = 1e-10 * mech.H

    def u_surf(x, y, nx, ny):
        """u on the ground surface; points where u is undefined (e.g. outside a FE mesh by round-off) are
        re-evaluated after a 1e-10 H shift into the soil (-n)"""
        uu = np.asarray(u(x, y), float)
        bad = ~np.isfinite(uu)
        if np.any(bad):
            uu = np.where(bad, u(x - eps * nx, y - eps * ny), uu)
        return uu

    total = np.zeros(M)
    # spiral (points strictly inside the soil)
    if mech.k > 0.0:
        th = mech.t1[:, None] + (mech.t2 - mech.t1)[:, None] * tq
        r = mech.r0[:, None] * np.exp((th - mech.t1[:, None]) * mech.k)
        xs = mech.Cx[:, None] - r * np.cos(th)
        ys = -mech.Cy[:, None] + r * np.sin(th)
        uu = u(xs, ys)
        total += mech.k * np.sum(uu * r * r * wq, 1) * (mech.t2 - mech.t1)
    # crest A -> O (y = 0, outward normal (0, -1), U.n = -(Cx - x))
    xs = mech.Ax[:, None] * (1.0 - tq)
    uu = u_surf(xs, np.zeros(xs.shape), 0.0, -1.0)
    total -= np.sum(uu * (-(mech.Cx[:, None] - xs)) * wq, 1) * mech.L
    # face O -> F (F = B for I, T for II), outward normal (sin b, -cos b)
    sb, cb = mech.sb, mech.cb
    lF = np.where(mech.II, mech.H, mech.By) / sb                       # face length used
    breaks = [np.zeros(M), lF]
    if hw is not None and hw > 0.0:
        breaks.insert(1, np.clip(hw / sb, 0.0, lF))
    for a, b in zip(breaks[:-1], breaks[1:]):
        ll = b - a
        sq = a[:, None] + ll[:, None] * tq
        xs, ys = sq * cb, sq * sb
        uu = u_surf(xs, ys, sb, -cb)
        un = (ys + mech.Cy[:, None]) * sb - (mech.Cx[:, None] - xs) * cb
        total -= np.sum(uu * un * wq, 1) * ll
    # toe ground T -> B (y = H)
    d = np.where(mech.II, mech.Bx - mech.xtoe, 0.0)
    xs = mech.xtoe + d[:, None] * tq
    uu = u_surf(xs, np.full(xs.shape, mech.H), 0.0, -1.0)
    total -= np.sum(uu * (-(mech.Cx[:, None] - xs)) * wq, 1) * d
    return total


# ======================================================================================================
# objective
# ======================================================================================================
class Problem:
    """Slope + soil + seepage field; evaluates Gamma(theta1, theta2, s) for batches of mechanisms."""

    def __init__(self, beta_deg, H, c, phi_deg, gamma, gamma_w=GAMMA_W, field=None, pu_method="auto",
                 d_max=10.0):
        self.beta_deg, self.H, self.c, self.phi_deg = float(beta_deg), float(H), float(c), float(phi_deg)
        self.gamma, self.gamma_w = float(gamma), float(gamma_w)
        self.gamma_p = self.gamma - self.gamma_w
        self.field = field
        self.zero_field = _is_zero_field(field)
        if pu_method == "auto":
            pu_method = "boundary" if (not self.zero_field and _is_potential_field(field)) else "domain"
        self.pu_method = pu_method
        self.d_max = float(d_max)
        self.n_eval = 0

    def bounds(self, kind):
        b = np.radians(self.beta_deg)
        if kind == "I":
            return np.array([_EPS_ANG, _EPS_ANG, 1e-4]), np.array([np.pi - _EPS_ANG, np.pi - b, 1.0])
        return np.array([_EPS_ANG, _EPS_ANG, 1.0]), np.array([np.pi - _EPS_ANG, np.pi - _EPS_ANG, 1.0 + self.d_max])

    def mechanisms(self, X):
        X = np.atleast_2d(X)
        return Mechanisms(X[:, 0], X[:, 1], X[:, 2], self.beta_deg, self.H, self.phi_deg)

    def powers(self, mech, quad="fine", pgamma="closed"):
        """(P_mr, P_gamma, P_u) per unit omega for admissible mechanisms."""
        Pmr = mech.Pmr(self.c)
        if pgamma == "closed":
            Pg = mech.Pgamma(self.gamma_p)
        else:
            Pg = domain_power(mech, UniformForce(0.0, self.gamma_p).force, quad)
        if self.zero_field or len(mech) == 0:
            Pu = np.zeros(len(mech))
        elif self.pu_method == "boundary":
            Pu = boundary_power(mech, self.field.u, hw=getattr(self.field, "hw", None),
                                **({} if quad in ("coarse", "medium") else dict(n_panels=48, q=8)))
        else:
            qd = get_quad(quad)
            Pu = domain_power(mech, self.field.force, qd, field_circles(self.field, qd.n_circle_grade))
        return Pmr, Pg, Pu

    def gamma_factor(self, X, kind=None, quad="fine"):
        """Gamma = P_mr / (P_gamma + P_u) for rows X = (theta1, theta2, s); BIG if inadmissible."""
        X = np.atleast_2d(np.asarray(X, float))
        mech = self.mechanisms(X)
        ok = mech.ok.copy()
        if kind == "I":
            ok &= X[:, 2] <= 1.0
        elif kind == "II":
            ok &= X[:, 2] >= 1.0
        out = np.full(len(X), BIG)
        if np.any(ok):
            sub = mech.subset(np.nonzero(ok)[0])
            Pmr, Pg, Pu = self.powers(sub, quad)
            Pext = Pg + Pu
            with np.errstate(all="ignore"):
                G = np.where(Pext > 0.0, Pmr / Pext, BIG)
            out[ok] = np.where(np.isfinite(G), G, BIG)
        self.n_eval += len(X)
        return out


# ======================================================================================================
# optimisation
# ======================================================================================================
def pso(fun, lb, ub, n_particles=40, n_iter=150, seed=0, feasible=None, w=0.7298, c1=1.49618, c2=1.49618,
        stall=40, vmax_frac=0.2):
    """Global-best Particle Swarm Optimisation (constriction coefficients of Clerc & Kennedy) on a box.
    fun(X) -> values (n,), BIG for inadmissible points; feasible(X) -> bool (cheap geometric test) is used
    to draw admissible initial particles.  Deterministic for a given seed."""
    rng = np.random.default_rng(seed)
    lb, ub = np.asarray(lb, float), np.asarray(ub, float)
    dim = len(lb)
    span = ub - lb
    # initial swarm: admissible random points (rejection sampling), half of them the best of a larger draw
    cand = np.empty((0, dim))
    for _ in range(50):                      # up to 50 x 2000 * n_particles candidates (cheap geometric test)
        draw = lb + rng.random((2000 * n_particles if feasible is not None else 4 * n_particles, dim)) * span
        if feasible is not None:
            draw = draw[feasible(draw)]
        cand = np.vstack([cand, draw])
        if len(cand) >= 4 * n_particles:
            break
    if len(cand) == 0:
        return None, BIG, 0
    pool = cand[: 4 * n_particles]
    fp = fun(pool)
    good = pool[fp < BIG]
    fgood = fp[fp < BIG]
    if len(good) == 0:
        return pool[0], BIG, len(pool)
    order = np.argsort(fgood)
    nb = min(n_particles // 2, len(good))
    pick = list(order[:nb])
    rest = order[nb:]
    if len(rest):
        pick += list(rng.choice(rest, size=min(len(rest), n_particles - nb), replace=False))
    X = good[pick]
    F = fgood[pick]
    if len(X) < n_particles:  # fill with random points of the box
        extra = lb + rng.random((n_particles - len(X), dim)) * span
        X = np.vstack([X, extra])
        F = np.concatenate([F, fun(extra)])
    V = (rng.random(X.shape) - 0.5) * 0.2 * span
    vmax = vmax_frac * span
    P, PF = X.copy(), F.copy()
    g = int(np.argmin(PF))
    best_hist = [PF[g]]
    n_eval = len(pool) + len(X)
    it_last = 0
    for it in range(n_iter):
        r1, r2 = rng.random(X.shape), rng.random(X.shape)
        V = w * V + c1 * r1 * (P - X) + c2 * r2 * (P[g] - X)
        V = np.clip(V, -vmax, vmax)
        X = np.clip(X + V, lb, ub)
        F = fun(X)
        n_eval += len(X)
        better = F < PF
        P[better], PF[better] = X[better], F[better]
        g = int(np.argmin(PF))
        if PF[g] < best_hist[-1] * (1.0 - 1e-9):
            it_last = it
        best_hist.append(PF[g])
        if it - it_last > stall:
            break
    return P[g].copy(), float(PF[g]), n_eval


def nelder_mead(fun1, x0, lb, ub, step=(0.02, 0.02, 0.01), restarts=2):
    """Bounded Nelder-Mead polish (restarted) of a scalar objective."""
    x = np.asarray(x0, float)
    f = fun1(x)
    for _ in range(restarts):
        sim = [x]
        for i in range(len(x)):
            e = x.copy()
            e[i] = e[i] + step[i] if e[i] + step[i] <= ub[i] else e[i] - step[i]
            sim.append(e)
        res = minimize(fun1, x, method="Nelder-Mead", bounds=list(zip(lb, ub)),
                       options=dict(initial_simplex=np.array(sim), xatol=1e-10, fatol=1e-13, maxiter=4000,
                                    maxfev=8000))
        if res.fun <= f:
            x, f = res.x, res.fun
        step = tuple(0.2 * s for s in step)
    return x, f


def stability_factor(beta_deg, H, c, phi_deg, gamma, gamma_w=GAMMA_W, field=None, mechanisms=("I", "II"),
                     seeds=(0, 1, 2), n_particles=40, n_iter=150, quad_search="coarse", quad_final="fine",
                     pu_method="auto", d_max=10.0, polish=True, verbose=False):
    """Upper-bound stability factor Gamma = min P_mr / (P_gamma + P_u) (Eq. 55), gamma' = gamma - gamma_w.

    field: seepage-force field (field.force(x, y) = -grad u, paper coordinates), None / NoSeepage for h_w = 0.
    Returns dict(Gamma, Hcrit = Gamma * H, mechanism, params, A, B, C, ...); Gamma = inf when no admissible
    mechanism with positive external power exists (e.g. beta < phi without seepage)."""
    prob = Problem(beta_deg, H, c, phi_deg, gamma, gamma_w, field, pu_method, d_max)
    t0 = time.time()
    runs = []
    for kind in mechanisms:
        lb, ub = prob.bounds(kind)

        def fsearch(X, kind=kind):
            return prob.gamma_factor(X, kind, quad_search)

        def feas(X, kind=kind):
            m = prob.mechanisms(X)
            ok = m.ok
            return ok & ((X[:, 2] <= 1.0) if kind == "I" else (X[:, 2] >= 1.0))

        for seed in seeds:
            x, f, ne = pso(fsearch, lb, ub, n_particles, n_iter, seed, feas)
            if x is None or f >= BIG:
                runs.append(dict(kind=kind, seed=seed, x=None, Gamma_pso=np.inf, Gamma=np.inf, n_eval=ne))
                continue
            if polish:
                xp, fp = nelder_mead(lambda z, kind=kind: prob.gamma_factor(z, kind, quad_final)[0], x, lb, ub)
            else:
                xp, fp = x, prob.gamma_factor(x, kind, quad_final)[0]
            runs.append(dict(kind=kind, seed=seed, x=xp, x_pso=x, Gamma_pso=f, Gamma=fp, n_eval=ne))
            if verbose:
                print(f"   {kind} seed {seed}: PSO {f:.6f} ({ne} evals) -> NM {fp:.6f}  x={np.round(xp, 6)}",
                      flush=True)
    best = min(runs, key=lambda r: r["Gamma"]) if runs else None
    res = dict(beta_deg=beta_deg, H=H, c=c, phi_deg=phi_deg, gamma=gamma, gamma_w=gamma_w, runs=runs,
               pu_method=prob.pu_method if not prob.zero_field else "none", time=time.time() - t0)
    if best is None or not np.isfinite(best["Gamma"]) or best["Gamma"] >= BIG:
        res.update(Gamma=np.inf, Hcrit=np.inf, mechanism=None, params=None, A=None, B=None, C=None)
        return res
    x = best["x"]
    m = prob.mechanisms(x)
    Pmr, Pg, Pu = prob.powers(m, quad_final)
    G = float(Pmr[0] / (Pg[0] + Pu[0]))
    kind = "II" if x[2] > 1.0 else "I"
    params = dict(theta1=float(x[0]), theta2=float(x[1]))
    if kind == "I":
        params["eta"] = float(x[2])
    else:
        params["d_over_H"] = float(x[2] - 1.0)
    lb, ub = prob.bounds(best["kind"])
    names = ("theta1", "theta2", "s")
    at_bound = [f"{names[i]}={'lower' if abs(x[i] - lb[i]) < 1e-6 * (ub[i] - lb[i]) else 'upper'}"
                for i in range(3) if min(abs(x[i] - lb[i]), abs(x[i] - ub[i])) < 1e-6 * (ub[i] - lb[i])]
    gam = [r["Gamma"] for r in runs if np.isfinite(r["Gamma"]) and r["kind"] == best["kind"]]
    res.update(Gamma=G, Hcrit=G * H, mechanism=kind, params=params, x=x, A=m.A()[0], B=m.B()[0], C=m.C()[0],
               r0=float(m.r0[0]), rh=float(m.rh[0]), L=float(m.L[0]), theta_O=float(m.thO[0]),
               P_mr=float(Pmr[0]), P_gamma=float(Pg[0]), P_u=float(Pu[0]), at_search_bound=at_bound,
               d_max_reached=bool(kind == "II" and x[2] > 1.0 + 0.999 * d_max),
               seed_spread=float(max(gam) / min(gam) - 1.0) if gam else np.nan,
               best_by_class={k: min([r["Gamma"] for r in runs if r["kind"] == k] or [np.inf]) for k in mechanisms},
               n_eval=prob.n_eval)
    return res


def critical_height(beta_deg, c, phi_deg, gamma, gamma_w=GAMMA_W, field_factory=None, H_ref=10.0, **kw):
    """H_crit = Gamma(H_ref) * H_ref (exact by similarity at fixed h_w/H, see the spec).
    field_factory(H) -> field for the reference height (None: no seepage)."""
    field = field_factory(H_ref) if field_factory is not None else None
    r = stability_factor(beta_deg, H_ref, c, phi_deg, gamma, gamma_w, field, **kw)
    return r["Hcrit"], r


# ======================================================================================================
# verification checks
# ======================================================================================================
def _random_admissible(prob, kind, n, rng, max_draw=400000):
    lb, ub = prob.bounds(kind)
    out = np.empty((0, 3))
    for _ in range(20):
        X = lb + rng.random((max_draw // 20, 3)) * (ub - lb)
        m = prob.mechanisms(X)
        ok = m.ok & ((X[:, 2] <= 1.0) if kind == "I" else (X[:, 2] >= 1.0))
        out = np.vstack([out, X[ok]])
        if len(out) >= n:
            break
    return out[:n]


def polygon_moment(mech, i, n_spiral=20001):
    """Independent check of the region: polygon A-O-(T)-B + discretised spiral B->A (shoelace formulas),
    Richardson-extrapolated in the number of spiral chords.  Returns int_Omega (Cx - x) dA."""
    a, b = _polygon_moment(mech, i, n_spiral), _polygon_moment(mech, i, 2 * n_spiral - 1)
    return b + (b - a) / 3.0


def _polygon_moment(mech, i, n_spiral):
    pts = [(mech.Ax[i], 0.0), (0.0, 0.0)]
    if mech.II[i]:
        pts.append((mech.xtoe, mech.H))
    pts.append((mech.Bx[i], mech.By[i]))
    th = np.linspace(mech.t2[i], mech.t1[i], n_spiral)[1:-1]
    r = mech.r0[i] * np.exp((th - mech.t1[i]) * mech.k)
    P = np.vstack([np.array(pts), np.stack([mech.Cx[i] - r * np.cos(th), -mech.Cy[i] + r * np.sin(th)], 1)])
    x, y = P[:, 0], P[:, 1]
    x1, y1 = np.roll(x, -1), np.roll(y, -1)
    cr = x * y1 - x1 * y
    A = 0.5 * cr.sum()
    Mx = (x + x1) @ cr / 6.0
    if A < 0:
        A, Mx = -A, -Mx
    return mech.Cx[i] * A - Mx


def check_closed_form(log, quick):
    log("=" * 100)
    log("1) P_gamma: closed form f1 - f2 - f3 (- f4 for II) vs polar quadrature ('xfine') vs polygon (shoelace)")
    log("   gamma' = 1, H = 1; random admissible mechanisms; rel = |a - b| / |closed form|; 'max rel (r_h<=20H)'")
    log("   = max over mechanisms of moderate size; 'max scaled' = max |a - b| / (|f1| + |f2| + |f3| + |f4|) over all")
    log("   (huge near-degenerate random mechanisms, r0 ~ 1e3 H, have P_gamma ~ 1e-4 of the f terms: cancellation)")
    log("=" * 100)
    rng = np.random.default_rng(12345)
    n = 60 if quick else 300
    log(f"{'beta':>5} {'phi':>5} {'mech':>4} {'n':>4} {'max rel(r_h<=20H)':>18} {'median rel':>11} {'max scaled':>11} "
        f"{'Eq53 vs sinform':>16} {'polygon max rel':>16}")
    worst = 0.0
    for beta in (30.0, 45.0, 60.0, 90.0):
        for phi in (0.0, 20.0, 32.0):
            for kind in ("I", "II"):
                prob = Problem(beta, 1.0, 1.0, phi, 1.0, 0.0)
                X = _random_admissible(prob, kind, n, rng)
                m = prob.mechanisms(X)
                cf = m.Pgamma(1.0)
                qd = domain_power(m, UniformForce(0.0, 1.0).force, "xfine")
                rel = np.abs(qd - cf) / np.abs(cf)
                ft = m.f_terms()
                scale = sum(np.abs(t) for t in ft)
                small = m.rh <= 20.0
                worst = max(worst, rel[small].max())
                f3s = ft[2]
                e53 = (f"{np.max(np.abs(m.f3_paper() - f3s) / np.abs(f3s)):16.2e}"
                       if kind == "I" and beta < 90 else f"{'-':>16}")
                ip = np.nonzero(small)[0][:8]
                pol = np.array([polygon_moment(m, i) for i in ip])
                prel = np.max(np.abs(pol - cf[ip]) / np.abs(cf[ip]))
                log(f"{beta:5.0f} {phi:5.1f} {kind:>4} {len(X):4d} {rel[small].max():18.2e} {np.median(rel):11.2e} "
                    f"{np.max(np.abs(qd - cf) / scale):11.2e} {e53} {prel:16.2e}")
    log(f"   worst closed-form vs quadrature relative difference (r_h <= 20 H): {worst:.2e}")
    # P_mr phi -> 0
    m0 = Mechanisms(0.5, 1.6, 1.0, 60.0, 1.0, 0.0)
    m1 = Mechanisms(0.5, 1.6, 1.0, 60.0, 1.0, 1e-7)
    log(f"   P_mr continuity at phi -> 0: P_mr(phi=0) = {m0.Pmr(1.0)[0]:.12f}, P_mr(phi=1e-7 deg) = "
        f"{m1.Pmr(1.0)[0]:.12f}")
    # P_mr by quadrature of c cos(phi) |U| dS along the spiral
    for phi in (0.0, 30.0):
        m = Mechanisms(0.5, 1.6, 1.0, 60.0, 1.0, phi)
        th, w = np.polynomial.legendre.leggauss(40)
        th = m.t1[0] + (m.t2[0] - m.t1[0]) * 0.5 * (th + 1)
        r = m.r0[0] * np.exp((th - m.t1[0]) * m.k)
        q = np.sum(0.5 * (m.t2[0] - m.t1[0]) * w * np.cos(np.radians(phi)) * r * r / np.cos(np.radians(phi)))
        log(f"   P_mr Eq. 48 vs int c cos(phi)|U| dS (phi = {phi:4.1f}): {m.Pmr(1.0)[0]:.12f} {q:.12f}")


def check_potential(log, quick):
    log("=" * 100)
    log("2) P_u: polar domain quadrature ('xfine') vs boundary formula -int u U.n dS for smooth potentials")
    log("=" * 100)
    rng = np.random.default_rng(7)
    n = 40 if quick else 120
    fields = [("uniform (0,1)", UniformForce(0.0, 1.0)), ("uniform (0.7,-0.3)", UniformForce(0.7, -0.3)),
              ("x^2y-y^3+xy+x", PolynomialPotential())]
    for beta, phi in ((45.0, 30.0), (90.0, 20.0), (30.0, 10.0)):
        for kind in ("I", "II"):
            prob = Problem(beta, 1.0, 1.0, phi, 1.0, 0.0)
            m = prob.mechanisms(_random_admissible(prob, kind, n, rng))
            row = []
            for name, f in fields:
                d = domain_power(m, f.force, "xfine")
                b = boundary_power(m, f.u, n_panels=48, q=10)
                rel = np.abs(d - b) / np.abs(d)
                row.append(f"{name}: med {np.median(rel):.1e} max(r_h<=20H) {np.max(rel[m.rh <= 20.0]):.1e}")
            log(f"   beta={beta:4.0f} phi={phi:4.1f} {kind:>2}: " + " | ".join(row))


def check_vertical_cut_and_literature(log, quick):
    log("=" * 100)
    log("3) Dry homogeneous slopes: N = gamma H_c / c (gamma_w = 0, no seepage) vs literature")
    log("=" * 100)
    # phi = 0: Taylor (1948) stability numbers 1/0.261, 1/0.219, 1/0.191, 1/0.181 (beta <= 53 deg, deep base
    # failure, unlimited depth); Chen (1975) log-spiral (= circle for phi = 0) 3.83, 4.57, 5.25, 5.53.
    # beta = 90 row: Chen (1975) vertical-cut values as remembered (not re-checked against the book);
    # beta = 60 row: Chen (1975) table values as remembered, low confidence.  (Checker note: 16.18 may be the
    # (beta=45, phi=20) entry of that table; the independent re-computation in check_limit_analysis.py gives
    # N(60, 30) = 16.035 and N(45, 20) = 16.161, identical to this module.  Not verified against the book.)
    lit_high = {(90, 0): 3.83, (75, 0): 4.57, (60, 0): 5.24, (45, 0): 5.52, (30, 0): 5.52}
    lit_low = {(60, 5): 6.17, (60, 10): 7.26, (60, 15): 8.65, (60, 20): 10.39, (60, 25): 12.75, (60, 30): 16.18}
    lit_med = {(90, 5): 4.19, (90, 10): 4.59, (90, 15): 5.02, (90, 20): 5.51, (90, 25): 6.06, (90, 30): 6.69,
               (90, 35): 7.43, (90, 40): 8.30}
    betas = (90, 75, 60, 45, 30) if not quick else (90, 60, 45)
    phis = (0, 5, 10, 15, 20, 25, 30, 35, 40) if not quick else (0, 20, 30)
    seeds = (0, 1) if quick else (0, 1, 2)
    rows = []
    log(f"{'beta':>5} {'phi':>4} {'N ours':>10} {'mech':>4} {'params (theta1, theta2, eta | d/H)':>40} "
        f"{'literature':>10} {'conf':>6} {'3.83tan(45+phi/2)':>18}")
    for beta in betas:
        for phi in phis:
            r = stability_factor(beta, 1.0, 1.0, phi, 1.0, 0.0, None, seeds=seeds)
            if (beta, phi) in lit_high:
                lit, conf = lit_high[(beta, phi)], "high"
            elif (beta, phi) in lit_med:
                lit, conf = lit_med[(beta, phi)], "medium"
            elif (beta, phi) in lit_low:
                lit, conf = lit_low[(beta, phi)], "low"
            else:
                lit, conf = np.nan, ""
            p = r["params"]
            ps = (f"({p['theta1']:.4f}, {p['theta2']:.4f}, " + (f"eta={p['eta']:.4f})" if "eta" in p else
                                                                 f"d/H={p['d_over_H']:.3f})")) if p else "-"
            extra = f"{3.83 * np.tan(np.radians(45 + phi / 2)):18.3f}" if beta == 90 else ""
            flag = "  [d = d_max: deeper mechanisms lower N further]" if r.get("d_max_reached") else ""
            log(f"{beta:5d} {phi:4d} {r['Gamma']:10.4f} {str(r['mechanism']):>4} {ps:>40} "
                f"{lit:10.2f} {conf:>6} {extra}{flag}")
            rows.append((beta, phi, r["Gamma"], r["mechanism"], lit, conf))
    return rows


def check_scale(log, quick):
    log("=" * 100)
    log("4) Scale invariance: Gamma(H) / Gamma(2H) = 2 (P_mr ~ c H^2, P_gamma, P_u ~ gamma H^3 at fixed h_w/H)")
    log("=" * 100)
    from analytical_seepage import AnalyticalSeepage
    cases = [("dry London soil beta=60", 60.0, 6.0, 32.0, 18.0, 9.81, None),
             ("analytical seepage beta=60 alpha=1 hw=H", 60.0, 10.0, 30.0, 20.0, 9.81, "an")]
    for name, beta, c, phi, g, gw, fld in cases:
        out = []
        for H in (5.0, 10.0):
            f = AnalyticalSeepage(beta, H, H, 1.0, gamma_w=gw) if fld == "an" else None
            out.append(stability_factor(beta, H, c, phi, g, gw, f, seeds=(0,)))
        r1, r2 = out
        log(f"   {name}: Gamma(H=5) = {r1['Gamma']:.8f}, Gamma(H=10) = {r2['Gamma']:.8f}, ratio = "
            f"{r1['Gamma'] / r2['Gamma']:.8f}; Hcrit = {r1['Hcrit']:.6f} / {r2['Hcrit']:.6f}; angles "
            f"{np.round(r1['x'][:2], 5)} / {np.round(r2['x'][:2], 5)}")


def check_uniform_field(log, quick):
    log("=" * 100)
    log("5) Uniform body-force field f = (0, g): Gamma(gamma', f) must equal Gamma(gamma' + g, no field)")
    log("=" * 100)
    for beta, c, phi, gam, gw, g in ((60.0, 6.0, 32.0, 18.0, 9.81, 4.0), (45.0, 10.0, 30.0, 20.0, 9.81, 9.81),
                                     (90.0, 10.0, 0.0, 18.0, 0.0, 2.0)):
        ra = stability_factor(beta, 5.0, c, phi, gam, gw, UniformForce(0.0, g), seeds=(0,), pu_method="domain")
        rb = stability_factor(beta, 5.0, c, phi, gam, gw, UniformForce(0.0, g), seeds=(0,), pu_method="boundary")
        rc = stability_factor(beta, 5.0, c, phi, gam + g, gw, None, seeds=(0,))
        log(f"   beta={beta:4.0f} c={c:5.1f} phi={phi:4.1f} gamma'={gam - gw:6.2f} g={g:5.2f}: "
            f"domain {ra['Gamma']:.10f}  boundary {rb['Gamma']:.10f}  gamma'+g {rc['Gamma']:.10f}  "
            f"rel {ra['Gamma'] / rc['Gamma'] - 1:+.1e} / {rb['Gamma'] / rc['Gamma'] - 1:+.1e}")


def _grid_min(prob, kind, n1, n2, n3, quad, chunk=200000):
    lb, ub = prob.bounds(kind)
    g1 = np.linspace(lb[0], ub[0], n1)
    g2 = np.linspace(lb[1], ub[1], n2)
    g3 = np.linspace(lb[2], ub[2], n3)
    G = np.stack(np.meshgrid(g1, g2, g3, indexing="ij"), -1).reshape(-1, 3)
    best, xb = BIG, None
    for i0 in range(0, len(G), chunk):
        X = G[i0:i0 + chunk]
        f = prob.gamma_factor(X, kind, quad)
        j = int(np.argmin(f))
        if f[j] < best:
            best, xb = f[j], X[j]
    return best, xb, len(G)


def check_brute_force(log, quick):
    log("=" * 100)
    log("6) PSO (+ Nelder-Mead) vs brute-force grid search over (theta1, theta2, s) for each mechanism class")
    log("   (grid best then polished by Nelder-Mead from the grid point)")
    log("=" * 100)
    from analytical_seepage import AnalyticalSeepage
    cases = [("dry London soil (c=6, phi=32), gamma'=8.19, beta=60", 60.0, 10.0, 6.0, 32.0, 18.0, 9.81, None,
              (160, 160, 60), "fine"),
             ("dry phi=0, beta=45 (mechanism II, d<=10H)", 45.0, 1.0, 1.0, 0.0, 1.0, 0.0, None,
              (160, 160, 60), "fine"),
             ("Fig. 9 alpha=1 beta=60, analytical v'_opt field", 60.0, 5.0, 10.0, 30.0, 20.0, 9.81, "an",
              (50, 50, 20) if not quick else (30, 30, 12), "coarse"),
             ("beta=30, c=11.7, phi=24.7, gamma=18, h_w/H=0.2, analytical field (R_w inside the mechanism)",
              30.0, 10.0, 11.7, 24.7, 18.0, 9.81, "an0.2", (50, 50, 20) if not quick else (30, 30, 12), "coarse")]
    for name, beta, H, c, phi, g, gw, fld, ng, gq in cases:
        f = None
        if fld is not None:
            hw = H * (float(fld[2:]) if len(fld) > 2 else 1.0)
            f = AnalyticalSeepage(beta, H, hw, 1.0, gamma_w=gw)
        if f is None:
            gq_label = "closed-form P_gamma"
        else:
            gq_label = f"{gq} quadrature"
        prob = Problem(beta, H, c, phi, g, gw, f)
        t0 = time.time()
        r = stability_factor(beta, H, c, phi, g, gw, f, seeds=(0, 1, 2))
        tp = time.time() - t0
        log(f"   {name}: PSO+NM Gamma = {r['Gamma']:.8f} ({r['mechanism']}, x = {np.round(r['x'], 6)}), "
            f"seed spread {r['seed_spread']:.1e}, {r['n_eval']} evaluations, {tp:.1f} s")
        for kind in ("I", "II"):
            t0 = time.time()
            gb, xb, ntot = _grid_min(prob, kind, *ng, quad=gq)
            if xb is None or gb >= BIG:
                log(f"      grid {kind}: no admissible mechanism with P_ext > 0")
                continue
            lb, ub = prob.bounds(kind)
            xp, fp = nelder_mead(lambda z: prob.gamma_factor(z, kind, "fine")[0], xb, lb, ub)
            log(f"      grid {kind} {ng[0]}x{ng[1]}x{ng[2]} = {ntot} points ({gq_label}): min {gb:.6f} at "
                f"{np.round(xb, 4)} -> NM {fp:.8f} at {np.round(xp, 6)}; PSO class best "
                f"{r['best_by_class'][kind]:.8f}; diff {fp / r['best_by_class'][kind] - 1:+.1e}  ({time.time() - t0:.1f} s)")


def check_pu_convergence(log, quick):
    log("=" * 100)
    log("7) Convergence of the P_u quadrature at the optimal mechanisms (analytical field, Fig. 9 data alpha=1,")
    log("   h_w = H: discontinuities at r = R_w = R, R_e about O, r^(m-1) singularity at O)")
    log("=" * 100)
    log("   (last case: beta=30, H=10, h_w/H=0.2, c=11.7, phi=24.7, gamma=18: circle R_w inside the mechanism)")
    from analytical_seepage import AnalyticalSeepage
    for beta, H, hw, c, phi, g in ((45.0, 5.0, 5.0, 10.0, 30.0, 20.0), (60.0, 5.0, 5.0, 10.0, 30.0, 20.0),
                                   (90.0, 5.0, 5.0, 10.0, 30.0, 20.0), (30.0, 10.0, 2.0, 11.7, 24.7, 18.0)):
        f = AnalyticalSeepage(beta, H, hw, 1.0, gamma_w=9.81)
        r = stability_factor(beta, H, c, phi, g, 9.81, f, seeds=(0,))
        prob = Problem(beta, H, c, phi, g, 9.81, f)
        m = prob.mechanisms(r["x"])
        Pg = m.Pgamma(prob.gamma_p)[0]
        Pmr = m.Pmr(prob.c)[0]
        vals = {}
        for lev in ("coarse", "medium", "fine", "xfine", "ref"):
            qd = get_quad(lev)
            t0 = time.time()
            vals[lev] = (domain_power(m, f.force, qd, field_circles(f, qd.n_circle_grade))[0], time.time() - t0)
        ref = vals["ref"][0]
        log(f"   beta={beta:4.0f} h_w/H={hw / H:.1f} (m = {f.m:.4f}) P_gamma = {Pg:.6f}, P_u(ref) = {ref:.10f}")
        for lev, (v, t) in vals.items():
            log(f"      {lev:>6}: P_u = {v:.10f}  rel {v / ref - 1:+.2e}  Gamma {Pmr / (Pg + v):.8f}  ({t * 1e3:.1f} ms)")
        # without the circle splits (generic field treatment)
        for lev in ("fine", "xfine"):
            v = domain_power(m, f.force, lev, [])[0]
            log(f"      {lev:>6} without circle splits: rel {v / ref - 1:+.2e}")


def check_fe_field(log, quick):
    log("=" * 100)
    log("8) FE field (-grad u'_FE, fe_seepage.py defaults, ref=1, H=5, h_w=H, alpha=1): domain quadrature vs the")
    log("   boundary formula (exact for f = -grad u_h); the optimiser uses the boundary formula for FE fields")
    log("=" * 100)
    try:
        from fe_seepage import FESeepage
    except Exception as e:  # pragma: no cover
        log(f"   fe_seepage not importable: {e}")
        return
    dense = Quadrature(n_panels=120, q_theta=6, n_grade=6, grade_ratio=0.25, n_rho=60, q_rho=4, name="dense")
    for beta in ((45.0,) if quick else (45.0, 90.0)):
        fe = FESeepage(beta, H=5.0, hw=5.0, alpha=1.0, gamma_w=9.81, ref=1)
        log("   FE field (fe_seepage.py defaults at run time): " + fe.summary())
        r = stability_factor(beta, 5.0, 10.0, 30.0, 20.0, 9.81, fe, seeds=(0,))
        prob = Problem(beta, 5.0, 10.0, 30.0, 20.0, 9.81, fe)
        m = prob.mechanisms(r["x"])
        b = [boundary_power(m, fe.u, hw=fe.hw, n_panels=n, q=8)[0] for n in (24, 48, 192)]
        log(f"   beta={beta:4.0f}: Gamma_FE = {r['Gamma']:.6f} ({r['pu_method']}); boundary P_u with 24/48/192 panels: "
            f"{b[0]:.8f} {b[1]:.8f} {b[2]:.8f}")
        for lev in ("coarse", "fine", "xfine", "ref", dense):
            t0 = time.time()
            d = domain_power(m, fe.force, lev)[0]
            nm = lev if isinstance(lev, str) else lev.name
            log(f"      domain {nm:>6}: {d:.8f}  rel to boundary {d / b[2] - 1:+.2e}  ({time.time() - t0:.2f} s)")


def _paper_fig8():
    p = os.path.join(DATA_DIR, "paper_fig8.csv")
    out = {}
    if os.path.exists(p):
        import csv
        with open(p) as fh:
            for row in csv.DictReader(fh):
                if float(row["hw_over_H"]) == 0.0 and row["curve"] in ("vopt", "FE"):
                    out[(row["soil"], int(float(row["beta_deg"])), row["curve"])] = float(row["Hcrit_m"])
    pv = os.path.join(DATA_DIR, "paper_fig8_vertices.csv")
    if os.path.exists(pv):  # clipped end points (above the axis range) are only in the vertex file
        import csv
        with open(pv) as fh:
            for row in csv.DictReader(fh):
                key = (row["soil"], int(float(row["beta_deg"])), row["curve"])
                if float(row["hw_over_H"]) == 0.0 and row["curve"] in ("vopt", "FE") and key not in out:
                    out[key] = float(row["Hcrit_m"])
    return out


TABLE1 = {"London": dict(gamma=18.0, c=6.0, phi=32.0), "Israeli": dict(gamma=18.0, c=11.7, phi=24.7)}
# Parameters that reproduce the h_w/H = 0 ends of the Fig. 8 panels (check 9): the (c, phi) pairs of Table 1
# EXCHANGED between the two panels, gamma' = 18 - 9.8 = 8.2 kN/m^3 (60 deg points within 3e-5).
FIG8_PANEL_PARAMS_FITTED = {"London": dict(gamma=18.0, c=11.7, phi=24.7, gamma_w=9.8),
                            "Israeli": dict(gamma=18.0, c=6.0, phi=32.0, gamma_w=9.8)}
FIG8_PANELS = (("London", 30), ("London", 60), ("Israeli", 35), ("Israeli", 60))


def check_fig8_hw0(log, quick):
    log("=" * 100)
    log("9) Fig. 8 at h_w/H = 0: buoyant slope (gamma' = gamma - gamma_w, NoSeepage), H_crit = Gamma * H")
    log("   paper values: dashed (v'_opt) / dash-dot (FE) ends of paper_fig8*.csv (identical at h_w = 0)")
    log("=" * 100)
    paper = _paper_fig8()
    rows = []
    log(f"{'panel':>8} {'beta':>4} | {'soil params':>22} {'gamma_w':>7} | {'Hcrit ours':>11} {'paper':>9} "
        f"{'ours/paper':>10} | mechanism")
    for panel, beta in FIG8_PANELS:
        pap = paper.get((panel, beta, "vopt"), np.nan)
        for soil in ("London", "Israeli"):
            s = TABLE1[soil]
            for gw in (9.81, 9.8):
                r = stability_factor(beta, 10.0, s["c"], s["phi"], s["gamma"], gw, None, seeds=(0, 1, 2))
                tag = "Table 1" if soil == panel else "SWAPPED"
                desc = f"{soil} c={s['c']:.1f} phi={s['phi']:.1f}"
                if np.isfinite(r["Gamma"]):
                    p = r["params"]
                    mech = f"{r['mechanism']} theta1={p['theta1']:.4f} theta2={p['theta2']:.4f} " + \
                        (f"eta={p['eta']:.4f}" if "eta" in p else f"d/H={p['d_over_H']:.4f}") + \
                        f" L/H={r['L'] / 10.0:.4f}"
                else:
                    mech = "no admissible mechanism with P_gamma > 0 (Gamma = inf)"
                log(f"{panel:>8} {beta:4d} | {desc:>22} {gw:7.2f} | {r['Hcrit']:11.4f} {pap:9.3f} "
                    f"{r['Hcrit'] / pap:10.5f} | {tag} {mech}")
                rows.append(dict(panel=panel, beta_deg=beta, soil_params=soil, gamma_w=gw, Hcrit=r["Hcrit"],
                                 paper=pap, mechanism=r["mechanism"]))
    return rows


def check_beta_lt_phi(log, quick):
    log("=" * 100)
    log("10) beta < phi without seepage (London clay of Table 1, beta = 30 < phi = 32)")
    log("=" * 100)
    rng = np.random.default_rng(99)
    for beta, phi in ((30.0, 32.0), (31.9, 32.0), (32.1, 32.0)):
        prob = Problem(beta, 1.0, 1.0, phi, 1.0, 0.0, d_max=100.0)
        tot, npos, worst = 0, 0, -np.inf
        for kind in ("I", "II"):
            X = _random_admissible(prob, kind, 200000 if not quick else 40000, rng, max_draw=4000000)
            m = prob.mechanisms(X)
            pg = m.Pgamma(1.0) / m.rh ** 3          # P_gamma normalised by the mechanism size
            tot += len(X)
            npos += int(np.sum(pg > 0))
            worst = max(worst, float(pg.max()))
        log(f"   beta={beta:5.1f} phi={phi:4.1f}: {tot} random admissible mechanisms (I and II, d <= 100 H): "
            f"{npos} with P_gamma > 0, max P_gamma / (gamma' r_h^3) = {worst:+.3e}")
    for beta in (30.0, 31.0, 32.5, 33.0, 35.0):
        r = stability_factor(beta, 10.0, 6.0, 32.0, 18.0, 9.81, None, seeds=(0, 1), d_max=100.0)
        p = r["params"]
        log(f"   London soil beta={beta:4.1f}: Hcrit = {r['Hcrit']:.4f}" +
            (f" ({r['mechanism']}, theta1={p['theta1']:.4f} theta2={p['theta2']:.4f} L/H={r['L'] / 10:.4f})"
             if p else ""))


def _paper_fig9():
    p = os.path.join(DATA_DIR, "paper_fig9.csv")
    out = {}
    if os.path.exists(p):
        import csv
        with open(p) as fh:
            for row in csv.DictReader(fh):
                out[(int(row["alpha"]), int(float(row["beta_deg"])), row["curve"])] = float(row["Gamma"])
    return out


def check_fig9_smoke(log, quick):
    log("=" * 100)
    log("11) Fig. 9 smoke test: alpha = 1, H = 5, c = 10, phi = 30, gamma = 20, h_w = H, analytical v'_opt field")
    log("    (dashed curve of Fig. 9, paper_fig9.csv; digitisation uncertainty ~0.002-0.02)")
    log("=" * 100)
    from analytical_seepage import AnalyticalSeepage
    paper = _paper_fig9()
    betas = (45, 60, 90) if quick else (30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90)
    rows = []
    log(f"{'beta':>5} {'gamma_w':>7} {'Gamma ours':>11} {'paper':>8} {'ours/paper':>10} {'mech':>4} "
        f"{'theta1':>8} {'theta2':>8} {'eta|d/H':>8} {'P_gamma':>10} {'P_u':>10} {'spread':>8} {'time':>6}")
    for beta in betas:
        for gw in ((9.81, 9.8) if beta in (45, 60, 90) else (9.81,)):
            f = AnalyticalSeepage(beta, 5.0, 5.0, 1.0, gamma_w=gw)
            t0 = time.time()
            r = stability_factor(beta, 5.0, 10.0, 30.0, 20.0, gw, f, seeds=(0, 1, 2) if not quick else (0,))
            pap = paper.get((1, beta, "vopt"), np.nan)
            p = r["params"]
            log(f"{beta:5d} {gw:7.2f} {r['Gamma']:11.5f} {pap:8.4f} {r['Gamma'] / pap:10.4f} {r['mechanism']:>4} "
                f"{p['theta1']:8.4f} {p['theta2']:8.4f} {p.get('eta', p.get('d_over_H', np.nan)):8.4f} "
                f"{r['P_gamma']:10.4f} {r['P_u']:10.4f} {r['seed_spread']:8.1e} {time.time() - t0:6.1f}")
            rows.append(dict(beta_deg=beta, gamma_w=gw, Gamma=r["Gamma"], paper_vopt=pap, mechanism=r["mechanism"],
                             theta1=p["theta1"], theta2=p["theta2"], s=r["x"][2], P_gamma=r["P_gamma"], P_u=r["P_u"],
                             mech=r))
    return rows


def _plot_mechanisms(results, path, title):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    n = len(results)
    fig, axs = plt.subplots(1, n, figsize=(3.4 * n, 3.4))
    axs = np.atleast_1d(axs)
    for ax, (lab, r) in zip(axs, results):
        H, beta = r["H"], np.radians(r["beta_deg"])
        xt = H / np.tan(beta) if r["beta_deg"] < 90 else 0.0
        if not np.isfinite(r["Gamma"]):
            ax.set_title(lab + "\nGamma = inf", fontsize=8)
            continue
        m = Mechanisms(*r["x"], r["beta_deg"], H, r["phi_deg"])
        xl = min(m.Ax[0], 0) - 0.3 * H
        xr = max(m.Bx[0], xt) + 0.6 * H
        ax.plot([xl, 0, xt, xr], [0, 0, H, H], color="#0b0b0b", lw=1.5)
        th = np.linspace(m.t1[0], m.t2[0], 300)
        rr = m.r0[0] * np.exp((th - m.t1[0]) * m.k)
        ax.plot(m.Cx[0] - rr * np.cos(th), -m.Cy[0] + rr * np.sin(th), color="#2a78d6", lw=2)
        ax.plot([m.Ax[0], m.Cx[0], m.Bx[0]], [0, -m.Cy[0], m.By[0]], ls=":", color="#eb6834", lw=1)
        ax.plot([m.Cx[0]], [-m.Cy[0]], "o", color="#eb6834", ms=4)
        ax.set_aspect("equal")
        ax.invert_yaxis()
        ax.set_title(f"{lab}\nGamma={r['Gamma']:.4f} ({r['mechanism']})", fontsize=8)
        ax.tick_params(labelsize=7)
    fig.suptitle(title, fontsize=9)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def _plot_fig9(rows, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    paper = _paper_fig9()
    fig, ax = plt.subplots(figsize=(5.2, 3.6))
    pb = sorted(b for (a, b, cv) in paper if a == 1 and cv == "vopt")
    ax.plot(pb, [paper[(1, b, "vopt")] for b in pb], "s--", color="#eb6834", ms=5, lw=1.5,
            label="paper Fig. 9, dashed (K^-1 v'_opt)")
    rr = [r for r in rows if r["gamma_w"] == 9.81]
    ax.plot([r["beta_deg"] for r in rr], [r["Gamma"] for r in rr], "o-", color="#2a78d6", ms=5, lw=2,
            label="limit_analysis.py + analytical_seepage.py")
    ax.set_xlabel("slope inclination beta (deg)")
    ax.set_ylabel("stability factor Gamma")
    ax.set_title("alpha = 1, H = 5 m, c = 10 kPa, phi = 30, gamma = 20, h_w = H", fontsize=9)
    ax.grid(color="#e4e3df", lw=0.6)
    ax.legend(fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


CHECKS = ("closed", "potential", "literature", "scale", "uniform", "brute", "pu_conv", "fe", "fig8", "beta_lt_phi",
          "fig9")


def run_checks(quick=False, only=None, plot=True):
    os.makedirs(RESULTS_DIR, exist_ok=True)
    out_path = os.path.join(RESULTS_DIR, "checks_output.txt" if not quick else "checks_output_quick.txt")
    fh = open(out_path, "w")

    def log(s=""):
        print(s, flush=True)
        fh.write(s + "\n")
        fh.flush()

    t0 = time.time()
    log(f"limit_analysis.py checks ({'quick' if quick else 'full'}), quadrature levels: "
        + "; ".join(f"{k}: {get_quad(k)}" for k in QUAD_LEVELS))
    sel = CHECKS if not only else tuple(only)
    import csv
    if "closed" in sel:
        check_closed_form(log, quick)
    if "potential" in sel:
        check_potential(log, quick)
    if "literature" in sel:
        rows = check_vertical_cut_and_literature(log, quick)
        with open(os.path.join(RESULTS_DIR, "dry_stability_numbers.csv"), "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["beta_deg", "phi_deg", "N_gammaHc_over_c", "mechanism", "literature", "confidence"])
            for r in rows:
                w.writerow([r[0], r[1], f"{r[2]:.6f}", r[3], "" if not np.isfinite(r[4]) else r[4], r[5]])
    if "scale" in sel:
        check_scale(log, quick)
    if "uniform" in sel:
        check_uniform_field(log, quick)
    if "brute" in sel:
        check_brute_force(log, quick)
    if "pu_conv" in sel:
        check_pu_convergence(log, quick)
    if "fe" in sel:
        check_fe_field(log, quick)
    if "fig8" in sel:
        rows = check_fig8_hw0(log, quick)
        with open(os.path.join(RESULTS_DIR, "fig8_hw0_buoyant.csv"), "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["panel", "beta_deg", "soil_params", "gamma_w", "Hcrit_m", "paper_Hcrit_m", "mechanism"])
            for r in rows:
                w.writerow([r["panel"], r["beta_deg"], r["soil_params"], r["gamma_w"], f"{r['Hcrit']:.5f}",
                            f"{r['paper']:.3f}", r["mechanism"]])
    if "beta_lt_phi" in sel:
        check_beta_lt_phi(log, quick)
    if "fig9" in sel:
        rows = check_fig9_smoke(log, quick)
        with open(os.path.join(RESULTS_DIR, "fig9_alpha1_vopt.csv"), "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["alpha", "beta_deg", "gamma_w", "Gamma", "paper_vopt", "mechanism", "theta1", "theta2", "s",
                        "P_gamma", "P_u"])
            for r in rows:
                w.writerow([1, r["beta_deg"], r["gamma_w"], f"{r['Gamma']:.6f}", r["paper_vopt"], r["mechanism"],
                            f"{r['theta1']:.6f}", f"{r['theta2']:.6f}", f"{r['s']:.6f}", f"{r['P_gamma']:.6f}",
                            f"{r['P_u']:.6f}"])
        if plot:
            _plot_fig9(rows, os.path.join(RESULTS_DIR, "fig9_alpha1_vopt.png"))
            sel_rows = [r for r in rows if r["gamma_w"] == 9.81 and r["beta_deg"] in (45, 60, 90)]
            _plot_mechanisms([(f"beta={r['beta_deg']}", r["mech"]) for r in sel_rows],
                             os.path.join(RESULTS_DIR, "fig9_alpha1_mechanisms.png"),
                             "Critical mechanisms, Fig. 9 data, alpha=1, analytical field")
    log(f"total time {time.time() - t0:.1f} s")
    fh.close()
    return out_path


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--no-plot", action="store_true")
    ap.add_argument("--only", default="", help="comma separated subset of " + ",".join(CHECKS))
    a = ap.parse_args(argv)
    only = [s for s in a.only.split(",") if s] or None
    run_checks(a.quick, only, not a.no_plot)


if __name__ == "__main__":
    main()
