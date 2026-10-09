#!/usr/bin/env python3
"""Semi-analytical optimal filtration-velocity field of Ceron et al. (IJNAMG 2025), Sect. 4.1.

The class of admissible velocity fields is Eq. (29) of the paper (polar coordinates of Fig. 3):

    zone 1, 0   <= r < R_w : v = -h1(r) h2'(theta) e_r + (r h1(r))' h2(theta) e_theta,  h1(R_w) = 0
    zone 2, R_w <= r < R   : v = h3(r) e_theta
    zone 3, R   <= r < R_e : v = h4(r) e_theta
    zone 4, r >= R_e       : v = 0

with R_w = h_w / sin(beta), R = H / sin(beta), R_e = sqrt(H^2 + (L_m + H / tan(beta))^2).
The optimum of J*(v') (Eq. 23) over this class is derived in analytical_seepage_derivation.md
(it does NOT rely on the printed Eqs. 31 and 39, which contain typos / a normalisation artefact).

Polar coordinates (origin O = crest edge) mapped to the shared paper coordinates (x right, y DOWN):
    x = -r cos(theta),  y = r sin(theta),  theta = atan2(y, -x)
    e_r = (-cos(theta), sin(theta)),  e_theta = (sin(theta), cos(theta))
theta = 0 is the crest (negative x axis), theta = pi - beta is the slope face, on which
e_theta = (sin(beta), -cos(beta)) is the OUTWARD normal (n = e_theta, Fig. 3).

Permeability K = k_v e_y (x) e_y + k_h (1 - e_y (x) e_y), alpha = k_h / k_v, so
K^-1 = diag(1/k_h, alpha/k_h) in (x, y).  Seepage force f = K^-1 . v'_opt  (approximates -grad u).

Run as a script to execute all the verification checks:
    python3 Projects/SlopeSeepageForces/scripts/analytical_seepage.py [--quick] [--no-plot]
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.interpolate import CubicHermiteSpline
from scipy.optimize import brentq, minimize, minimize_scalar

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.dirname(SCRIPT_DIR)
DATA_DIR = os.path.join(PROJECT_DIR, "data")
RESULTS_DIR = os.path.join(PROJECT_DIR, "results", "analytical_seepage")


# ----------------------------------------------------------------------------------------------
# helpers
# ----------------------------------------------------------------------------------------------
def _E(s, m):
    """E(s, m) = (1 - s^(m-1)) / (1 - m), continuous at m = 1 (-> ln s). s > 0."""
    ls = np.log(s)
    if abs(m - 1.0) < 1e-12:
        return ls
    return np.expm1((m - 1.0) * ls) / (m - 1.0)


def _cd(theta, alpha):
    """Angular weights of the anisotropic energy:
    c(theta) = cos^2 + alpha sin^2  (multiplies v_r^2, i.e. h2'^2)
    d(theta) = sin^2 + alpha cos^2  (multiplies v_theta^2, i.e. h2^2)."""
    s2 = np.sin(theta) ** 2
    c2 = 1.0 - s2
    return c2 + alpha * s2, s2 + alpha * c2


class NoSeepage:
    """f = 0 everywhere (h_w = 0)."""

    def force(self, x, y):
        x, y = np.broadcast_arrays(np.asarray(x, float), np.asarray(y, float))
        return np.zeros(x.shape), np.zeros(x.shape)


# ----------------------------------------------------------------------------------------------
# main class
# ----------------------------------------------------------------------------------------------
class AnalyticalSeepage:
    """Optimal velocity field v'_opt of the class (29) and seepage force f = K^-1 . v'_opt.

    Parameters
    ----------
    beta_deg : slope angle in degrees, 0 < beta <= 90
    H        : slope height (m)
    hw       : drawdown depth below the crest, 0 <= hw <= H (m)
    alpha    : anisotropy ratio k_h / k_v >= 1
    kh       : horizontal permeability (only scales v and J*, not f)
    gamma_w  : unit weight of water (kN/m^3)
    Lm_over_H: L_m / H, defines R_e (paper: 10)
    ode_method: integrator for the h2 shooting (paper: LSODA)
    """

    def __init__(self, beta_deg, H, hw, alpha, kh=1.0, gamma_w=9.81, Lm_over_H=10.0,
                 ode_method="LSODA", n_theta=2001):
        beta = np.radians(float(beta_deg))
        if not (0.0 < beta <= 0.5 * np.pi + 1e-12):
            raise ValueError("beta must be in (0, 90] degrees")
        if not (0.0 <= hw <= H * (1 + 1e-12)):
            raise ValueError("hw must be in [0, H]")
        if alpha < 1.0:
            raise ValueError("alpha = k_h/k_v must be >= 1")
        self.beta_deg = float(beta_deg)
        self.beta = beta
        self.H = float(H)
        self.hw = float(min(hw, H))
        self.alpha = float(alpha)
        self.kh = float(kh)
        self.kv = self.kh / self.alpha
        self.gamma_w = float(gamma_w)
        self.Lm = float(Lm_over_H) * self.H
        self.ode_method = ode_method
        self.n_theta = int(n_theta)

        self.sb, self.cb = np.sin(beta), np.cos(beta)
        self.tb = np.inf if abs(self.cb) < 1e-15 else self.sb / self.cb
        self.Theta = np.pi - beta                       # angular opening of zones 1 and 2
        self.xtoe = self.H * self.cb / self.sb          # toe abscissa H / tan(beta)
        self.R = self.H / self.sb
        self.Rw = self.hw / self.sb
        self.Re = np.hypot(self.H, self.Lm + self.xtoe)
        a = self.alpha
        # Eq. (34): A = int_0^Theta d(theta) dtheta
        self.A = (a + 1.0) * self.Theta / 2.0 - (a - 1.0) * np.sin(2.0 * beta) / 4.0
        # Eq. (40), zone-3 integral int_R^Re 4 / (B r) dr
        #    computed with r = sqrt(H^2 + X^2) (X = horizontal distance along the toe ground), which removes
        #    the sqrt(r^2 - H^2) endpoint singularity of B at r = R when beta -> 90 deg:  dr / r = X dX / r^2
        self.I3 = quad(lambda X: 4.0 * X / (self.B(np.hypot(self.H, X)) * (self.H ** 2 + X ** 2)),
                       self.xtoe, self.Lm + self.xtoe, limit=400, epsabs=1e-14, epsrel=1e-12)[0] \
            if self.Re > self.R else 0.0
        if self.hw > 0.0:
            self._solve_h2()
        else:
            self.m, self.C, self.D, self.h2e, self.Phi, self.F = 0.0, 0.0, self.A, 1.0, self.A, 1.0 / self.A
            self.degenerate = True
            self._h2 = lambda t: np.ones_like(np.asarray(t, float))
            self._dh2 = lambda t: np.zeros_like(np.asarray(t, float))
            self.ap = 0.0

    # ------------------------------------------------------------------ geometry
    def B(self, r):
        """Eq. (35): B(r, alpha) = 2 int_0^{pi - arcsin(H/r)} d(theta) dtheta."""
        a, H = self.alpha, self.H
        r = np.asarray(r, float)
        q = np.clip(H / r, -1.0, 1.0)
        return (a + 1.0) * (np.pi - np.arcsin(q)) - (a - 1.0) * H * np.sqrt(np.maximum(r * r - H * H, 0.0)) / (r * r)

    def in_soil(self, x, y):
        """Soil region in paper coordinates (y down): crest y = 0 (x <= 0), face, toe ground y = H."""
        x = np.asarray(x, float)
        y = np.asarray(y, float)
        face = np.where(x < self.xtoe, y >= x * (self.tb if np.isfinite(self.tb) else 0.0), y >= self.H)
        return np.where(x <= 0.0, y >= 0.0, face)

    # ------------------------------------------------------------------ h2 problem
    def _shoot(self, m, dense=False, t_eval=None, method=None):
        """Integrate (c h2')' = m d h2, h2(0) = 1, h2'(0) = 0 (Eq. 37-38) together with
        C = int c h2'^2 and D = int d h2^2 (Eq. 36).  State y = [h2, q = c h2', C(theta), D(theta)]."""
        a = self.alpha

        def rhs(t, y):
            c, d = _cd(t, a)
            h, q = y[0], y[1]
            return [q / c, m * d * h, q * q / c, h * h * d]

        method = method or self.ode_method
        sol = solve_ivp(rhs, (0.0, self.Theta), [1.0, 0.0, 0.0, 0.0], method=method,
                        rtol=1e-12, atol=1e-14, dense_output=dense, t_eval=t_eval)
        if not sol.success:
            raise RuntimeError("h2 shooting failed: " + sol.message)
        return sol

    def Phi_of_m(self, m):
        """Reduced objective (see derivation): Phi(m) = (1 + 1/m) * c_e h2'(Theta) / h2(Theta)
        = (1 + 1/m) * min_{h2(Theta)=1} (C + m D);  min_m Phi = min_h2 (sqrt C + sqrt D)^2 / h2(Theta)^2."""
        if m <= 0.0:
            return self.A
        y = self._shoot(m).y[:, -1]
        return (1.0 + 1.0 / m) * y[1] / y[0]

    def fixed_point_residual(self, m):
        """m - sqrt(C(m)/D(m)) for the ODE solution with parameter m (zero at the optimum)."""
        y = self._shoot(m).y[:, -1]
        return m - np.sqrt(y[2] / y[3])

    def _solve_h2(self):
        # 1) scan Phi(m) on a log grid, 2) refine with a bounded Brent search, 3) cross-check with
        #    the fixed point m = sqrt(C/D) (the stationarity condition of Phi).
        tgrid = np.linspace(np.log(1e-7), np.log(30.0), 81)
        ph = np.array([self.Phi_of_m(np.exp(t)) for t in tgrid])
        i = int(np.argmin(ph))
        self.degenerate = False
        if i == 0:
            # Phi increases from its limit Phi(0+) = A: the infimum is the m -> 0 limit of the class,
            # i.e. the purely tangential zone-1 field v = (k_h gw sin(beta) / A) e_theta (m = 0)
            self.degenerate = True
            m = 0.0
        else:
            lo, hi = tgrid[max(i - 1, 0)], tgrid[min(i + 1, len(tgrid) - 1)]
            res = minimize_scalar(lambda t: self.Phi_of_m(np.exp(t)), bounds=(lo, hi), method="bounded",
                                  options={"xatol": 1e-11})
            m = float(np.exp(res.x))
            # sharpen with the fixed-point equation (smooth root, Phi is flat at its minimum)
            try:
                f_lo, f_hi = self.fixed_point_residual(np.exp(lo)), self.fixed_point_residual(np.exp(hi))
                if f_lo * f_hi < 0.0:
                    m = brentq(self.fixed_point_residual, np.exp(lo), np.exp(hi), xtol=1e-14, rtol=1e-14)
            except (RuntimeError, ValueError):
                pass
        self.m = m
        if self.degenerate:
            self.C, self.D, self.h2e = 0.0, self.A, 1.0
            self._h2 = lambda t: np.ones_like(np.asarray(t, float))
            self._dh2 = lambda t: np.zeros_like(np.asarray(t, float))
            self.Phi = self.A
        else:
            th = np.linspace(0.0, self.Theta, self.n_theta)
            sol = self._shoot(m, t_eval=th, method="DOP853")
            h, q, Cc, Dc = sol.y
            c, d = _cd(th, self.alpha)
            dh = q / c
            cp = (self.alpha - 1.0) * np.sin(2.0 * th)       # c'(theta)
            d2h = (m * d * h - cp * dh) / c                  # h2''
            self._h2 = CubicHermiteSpline(th, h, dh)
            self._dh2 = CubicHermiteSpline(th, dh, d2h)
            self.C, self.D, self.h2e = Cc[-1], Dc[-1], h[-1]
            self.q_e = q[-1]
            self.Phi = (1.0 + 1.0 / m) * q[-1] / h[-1]
        # F = h2e^2 / (sqrt C + sqrt D)^2  (maximised);  J1* = -k_h gw^2 hw^2 F / 4
        self.F = self.h2e ** 2 / (np.sqrt(self.C) + np.sqrt(self.D)) ** 2
        # amplitude a' = a (1 - m) of zone 1:  v_r = -a' E(s,m) h2', v_th = a' (E + s^(m-1)) h2
        self.ap = self.kh * self.gamma_w * self.sb * self.h2e / (self.D * (1.0 + self.m))

    def h2(self, theta):
        return self._h2(theta)

    def dh2(self, theta):
        return self._dh2(theta)

    def h1(self, r):
        """h1(r) = a' E(r/R_w, m), with h1(R_w) = 0 (pairs with h2 normalised by h2(0) = 1)."""
        r = np.asarray(r, float)
        return self.ap * _E(r / self.Rw, self.m)

    def g(self, r):
        """g = r h1(r) (stream-function amplitude) and g' = (r h1)'."""
        r = np.asarray(r, float)
        s = r / self.Rw
        E = _E(s, self.m)
        return self.ap * r * E, self.ap * (E + s ** (self.m - 1.0))

    def h3(self, r):
        """Eq. (32)."""
        return self.kh * self.gamma_w * self.hw / (self.A * np.asarray(r, float))

    def h4(self, r):
        """Eq. (33)."""
        r = np.asarray(r, float)
        return 2.0 * self.kh * self.gamma_w * self.hw / (self.B(r) * r)

    # ------------------------------------------------------------------ functional
    def Jstar_parts(self):
        """Contributions of zones 1, 2, 3 to J*(v'_opt) (Eq. 40), absolute units."""
        if self.hw <= 0.0:
            return 0.0, 0.0, 0.0
        pre = -self.kh * self.gamma_w ** 2 * self.hw ** 2 / 4.0
        return pre * self.F, pre * 2.0 / self.A * np.log(self.H / self.hw), pre * self.I3

    def Jstar(self):
        return float(sum(self.Jstar_parts()))

    def Jstar_normalized(self):
        """J*(v'_opt) / (k_h H^2 gamma_w^2)  (negative; Fig. 5 plots minus this)."""
        return self.Jstar() / (self.kh * self.H ** 2 * self.gamma_w ** 2)

    # ------------------------------------------------------------------ fields
    def polar_velocity(self, r, th):
        """(v_r, v_theta) of v'_opt at polar points assumed to be inside the soil."""
        r, th = np.broadcast_arrays(np.asarray(r, float), np.asarray(th, float))
        vr = np.zeros(r.shape)
        vt = np.zeros(r.shape)
        if self.hw <= 0.0:
            return vr, vt
        ok = (r > 0.0) & (r < self.Re)
        z1 = ok & (r < self.Rw)
        z2 = ok & (r >= self.Rw) & (r < self.R)
        z3 = ok & (r >= self.R)
        if np.any(z1):
            s = r[z1] / self.Rw
            t1 = np.clip(th[z1], 0.0, self.Theta)
            E = _E(s, self.m)
            vr[z1] = -self.ap * E * self.dh2(t1)
            vt[z1] = self.ap * (E + s ** (self.m - 1.0)) * self.h2(t1)
        if np.any(z2):
            vt[z2] = self.h3(r[z2])
        if np.any(z3):
            vt[z3] = self.h4(r[z3])
        return vr, vt

    def velocity(self, x, y):
        """Darcy velocity v'_opt in paper coordinates (x right, y down); 0 outside soil / r >= R_e."""
        x, y = np.broadcast_arrays(np.asarray(x, float), np.asarray(y, float))
        vx = np.zeros(x.shape)
        vy = np.zeros(x.shape)
        if self.hw <= 0.0:
            return vx, vy
        r = np.hypot(x, y)
        th = np.arctan2(y, -x)
        act = self.in_soil(x, y) & (r > 0.0) & (r < self.Re)
        if np.any(act):
            vr, vt = self.polar_velocity(r[act], th[act])
            ct, st = np.cos(th[act]), np.sin(th[act])
            vx[act] = -vr * ct + vt * st
            vy[act] = vr * st + vt * ct
        return vx, vy

    def force(self, x, y):
        """Seepage force f = K^-1 . v'_opt (kN/m^3), paper coordinates (x right, y down)."""
        vx, vy = self.velocity(x, y)
        return vx / self.kh, self.alpha * vy / self.kh

    def summary(self):
        return (f"beta={self.beta_deg:5.1f} alpha={self.alpha:5.2f} hw/H={self.hw / self.H:5.3f} "
                f"m=sqrt(C/D)={self.m:.6f} C/h2e^2={self.C / self.h2e ** 2:.6f} D/h2e^2={self.D / self.h2e ** 2:.6f} "
                f"F={self.F:.6f} 1/A={1 / self.A:.6f} I3={self.I3:.6f} -J*n={-self.Jstar_normalized():.6f}"
                + ("  [degenerate m->0: purely tangential zone-1 field]" if self.degenerate and self.hw > 0 else ""))


# ----------------------------------------------------------------------------------------------
# independent checks
# ----------------------------------------------------------------------------------------------
def jstar_by_quadrature(md, g=None, h2=None, h3=None, h4=None, ngl=96):
    """Evaluate J*(v') = 1/2 int v.K^-1.v dOmega + int_{dOmega_u} u^d v.n dS by brute-force quadrature
    for a field of the class (29) given by callables
        g(r) -> (g, g')  with g = r h1(r), g(R_w) = 0;   h2(theta) -> (h2, h2');   h3(r);   h4(r).
    The energy uses the Cartesian K^-1 = diag(1/k_h, alpha/k_h) (cross terms included) and the boundary
    terms use the Cartesian outward normals of the face and of the toe ground (independent of the
    simplifications made in the derivation)."""
    g = g or md.g
    h2 = h2 or (lambda t: (md.h2(t), md.dh2(t)))
    h3 = h3 or md.h3
    h4 = h4 or md.h4
    kh, al, gw, hw = md.kh, md.alpha, md.gamma_w, md.hw
    xg, wg = np.polynomial.legendre.leggauss(ngl)

    def energy_density(vr, vt, th):
        vx = -vr * np.cos(th) + vt * np.sin(th)
        vy = vr * np.sin(th) + vt * np.cos(th)
        return 0.5 * (vx * vx + al * vy * vy) / kh

    def ang_int(r, th_max, vfun):
        th = 0.5 * th_max * (xg + 1.0)
        vr, vt = vfun(r, th)
        return 0.5 * th_max * np.sum(wg * energy_density(vr, vt, th)) * r

    n_face = np.array([md.sb, -md.cb])                  # outward normal of the face = e_theta(Theta)
    Th = md.Theta
    J = 0.0
    # zone 1: substitute r = R_w t^p to tame the r^(m-1) singularity at the corner
    if md.Rw > 0.0:
        mm = md.m
        p = 1.0 if (mm <= 0.0 or mm >= 1.0) else min(1.0 / mm, 25.0)

        def v1(r, th):
            gg, dg = g(r)
            hh, dhh = h2(th)
            return -gg / r * dhh, dg * hh

        def z1(t):
            r = md.Rw * t ** p
            jac = md.Rw * p * t ** (p - 1.0)
            vr, vt = v1(r, np.array([Th]))
            vn = (-vr[0] * np.cos(Th) + vt[0] * np.sin(Th)) * n_face[0] + (vr[0] * np.sin(Th) + vt[0] * np.cos(Th)) * n_face[1]
            ud = -gw * r * md.sb                          # u^d = -gamma_w y on the face above water
            return (ang_int(r, Th, v1) + ud * vn) * jac

        J += quad(z1, 0.0, 1.0, limit=500, epsabs=0.0, epsrel=1e-11)[0]
        # zone 2
        if md.R > md.Rw:
            def v2(r, th):
                return np.zeros_like(th), h3(r) * np.ones_like(th)

            def z2(r):
                vn = h3(r) * 1.0                          # v = h3 e_theta, n = e_theta on the face
                return ang_int(r, Th, v2) + (-gw * hw) * vn

            J += quad(z2, md.Rw, md.R, limit=200, epsabs=0.0, epsrel=1e-11)[0]
    # zone 3 (energy) and the toe-ground boundary term parametrised by x
    if md.Re > md.R and hw > 0.0:
        def v3(r, th):
            return np.zeros_like(th), h4(r) * np.ones_like(th)

        def z3(r):
            thm = np.pi - np.arcsin(min(md.H / r, 1.0))
            return ang_int(r, thm, v3)

        J += quad(z3, md.R, md.Re, limit=400, epsabs=0.0, epsrel=1e-11)[0]

        def bt(x):
            r = np.hypot(x, md.H)
            th = np.arctan2(md.H, -x)
            vy = h4(r) * np.cos(th)                       # v = h4 e_theta, e_theta.e_y = cos(theta)
            return (-gw * hw) * (-vy)                     # outward normal of the toe ground = -e_y
        J += quad(bt, md.xtoe, np.sqrt(md.Re ** 2 - md.H ** 2), limit=400, epsabs=0.0, epsrel=1e-11)[0]
    return J


def direct_h2_minimisation(md, n_el=300):
    """Direct minimisation of (sqrt C + sqrt D)^2 / h2(Theta)^2 over P1 finite-element h2 (h2(Theta)=1),
    without using the ODE: an independent check of F = 1 / Phi_min.  Returns (min value, sqrt(C/D))."""
    th = np.linspace(0.0, md.Theta, n_el + 1)
    hl = np.diff(th)
    xg, wg = np.polynomial.legendre.leggauss(4)
    N = n_el + 1
    ce = np.zeros(n_el)                     # C = sum_e ce_e (h_{e+1} - h_e)^2   (exactly >= 0)
    Md = np.zeros((N, N))
    for e in range(n_el):
        t = th[e] + 0.5 * hl[e] * (xg + 1.0)
        c, d = _cd(t, md.alpha)
        w = 0.5 * hl[e] * wg
        ce[e] = np.sum(w * c) / hl[e] ** 2
        phi = np.vstack([1.0 - (t - th[e]) / hl[e], (t - th[e]) / hl[e]])
        Md[e:e + 2, e:e + 2] += (phi * (w * d)) @ phi.T

    def CD(h):
        dh = np.diff(h)
        C = np.sum(ce * dh * dh)
        gC = np.zeros(N)
        gC[:-1] -= 2 * ce * dh
        gC[1:] += 2 * ce * dh
        Mh = Md @ h
        return C, gC, h @ Mh, 2 * Mh

    eps = 1e-14 * md.A

    def fun(z):
        h = np.append(z, 1.0)
        C, gC, D, gD = CD(h)
        sC, sD = np.sqrt(C + eps), np.sqrt(D)
        val = (sC + sD) ** 2
        grad = (sC + sD) * (gC / sC + gD / sD)
        return val, grad[:-1]

    best = None
    starts = [np.cosh(th) / np.cosh(th[-1]), 1.0 + 0.3 * ((th / th[-1]) ** 2 - 1.0), np.ones(N)]
    for start in starts:
        res = minimize(fun, start[:-1], jac=True, method="L-BFGS-B",
                       options={"maxiter": 50000, "maxcor": 50, "ftol": 1e-15, "gtol": 1e-12})
        if best is None or res.fun < best.fun:
            best = res
    C, _, D, _ = CD(np.append(best.x, 1.0))
    return best.fun, np.sqrt(C / D)


def _numerical_div(md, x, y, rel=1e-5):
    h = rel * np.hypot(x, y)
    vxp, _ = md.velocity(x + h, y)
    vxm, _ = md.velocity(x - h, y)
    _, vyp = md.velocity(x, y + h)
    _, vym = md.velocity(x, y - h)
    vx, vy = md.velocity(x, y)
    div = (vxp - vxm) / (2 * h) + (vyp - vym) / (2 * h)
    scale = np.hypot(vx, vy) / np.hypot(x, y)
    return div, scale


def _load_fig5_paper():
    path = os.path.join(DATA_DIR, "fig5_vector_fill_polygons.csv")
    if not os.path.exists(path):
        return None
    arr = np.loadtxt(path, delimiter=",", comments="#")
    return {(int(a), int(round(b))): (s, d) for a, b, s, d in arr}


# ----------------------------------------------------------------------------------------------
# driver
# ----------------------------------------------------------------------------------------------
def run_checks(quick=False, plot=True):
    np.set_printoptions(linewidth=160)
    os.makedirs(RESULTS_DIR, exist_ok=True)
    H, gw = 1.0, 9.81

    print("=" * 100)
    print("1) h2 problem: reduced 1-D minimisation over m (ODE, %s) vs fixed point m = sqrt(C/D) vs direct FE"
          % "LSODA")
    print("=" * 100)
    cases = [(15, 1), (30, 1), (60, 1), (90, 1), (30, 4), (60, 10), (90, 10)] if quick else \
        [(b, a) for a in (1, 2, 4, 10) for b in (15, 30, 45, 60, 75, 90)]
    print(f"{'beta':>5} {'alpha':>5} {'m*':>12} {'m-sqrt(C/D)':>12} {'Phi_min':>12} {'A':>10} "
          f"{'FE-direct (sC+sD)^2':>20} {'FE sqrt(C/D)':>12} {'identity':>10}")
    for b, a in cases:
        md = AnalyticalSeepage(b, H, H, a, gamma_w=gw)
        fe_val, fe_m = direct_h2_minimisation(md, n_el=200 if quick else 400)
        if md.degenerate:
            res, ident = 0.0, 0.0
        else:
            res = md.m - np.sqrt(md.C / md.D)
            ident = (md.q_e * md.h2e - (md.C + md.m * md.D)) / md.C
        print(f"{b:5.0f} {a:5.0f} {md.m:12.8f} {res:12.2e} {md.Phi:12.8f} {md.A:10.6f} {fe_val:20.8f} {fe_m:12.6f} {ident:10.1e}"
              + ("  (m->0 limit)" if md.degenerate else ""))

    print()
    print("=" * 100)
    print("2) -J*(v'_opt)/(k_h H^2 gamma_w^2) at h_w = H, L_m = 10 H   vs   Fig. 5 solid curves (vector data)")
    print("=" * 100)
    paper = _load_fig5_paper()
    betas = np.arange(15, 91, 5)
    rows = []
    for a in (1, 2, 4, 10):
        print(f"alpha = {a}")
        print(f"  {'beta':>5} {'-J*n ours':>10} {'paper':>8} {'rel.diff':>9} {'m':>9} {'F':>8} {'1/A':>8} {'I3':>8}")
        for b in betas:
            md = AnalyticalSeepage(b, H, H, a, gamma_w=gw)
            val = -md.Jstar_normalized()
            ref = paper.get((a, int(b)), (np.nan, np.nan))[0] if paper else np.nan
            rows.append((a, b, val, ref, md.m, md.F, 1 / md.A, md.I3))
            print(f"  {b:5.0f} {val:10.4f} {ref:8.4f} {100 * (val - ref) / ref:8.2f}% {md.m:9.5f} {md.F:8.5f} {1 / md.A:8.5f} {md.I3:8.5f}")
    out = os.path.join(RESULTS_DIR, "fig5_analytical_Jstar.csv")
    np.savetxt(out, np.array(rows), delimiter=",", fmt="%.6f",
               header="alpha,beta_deg,minusJstar_norm_ours,minusJstar_norm_paper_solid,m,F,inv_A,I3")
    print(f"  -> written {out}")

    print()
    print("=" * 100)
    print("3) dependence on L_m (h_w = H): -J*n; field inside r < R_e does not depend on L_m")
    print("=" * 100)
    for a in (1, 10):
        for b in (30, 60, 90):
            vals = [-AnalyticalSeepage(b, H, H, a, gamma_w=gw, Lm_over_H=L).Jstar_normalized() for L in (2, 5, 10, 20, 50, 100)]
            print(f"  alpha={a:3d} beta={b:3d}  Lm/H = 2,5,10,20,50,100 : " + " ".join(f"{v:.4f}" for v in vals))
    md10 = AnalyticalSeepage(45, H, 0.8, 4, Lm_over_H=10)
    md50 = AnalyticalSeepage(45, H, 0.8, 4, Lm_over_H=50)
    xs = np.array([-0.5, 0.2, 0.8, 1.5, -3.0])
    ys = np.array([0.3, 0.5, 1.1, 1.5, 2.0])
    f10 = np.array(md10.force(xs, ys))
    f50 = np.array(md50.force(xs, ys))
    print(f"  max |f(Lm=10H) - f(Lm=50H)| at 5 points inside R_e(10H): {np.max(np.abs(f10 - f50)):.2e}")

    print()
    print("=" * 100)
    print("4) admissibility: div v = 0 (central differences), continuity of v.e_r across r = R_w, R, R_e")
    print("=" * 100)
    rng = np.random.default_rng(1)
    for (b, a, hwr) in [(30, 1, 1.0), (60, 4, 0.6), (45, 10, 0.3), (90, 2, 1.0), (20, 2, 0.5)]:
        md = AnalyticalSeepage(b, H, hwr * H, a, gamma_w=gw)
        worst = []
        for (r0, r1) in [(0.05 * md.Rw, 0.95 * md.Rw), (1.02 * md.Rw, 0.98 * md.R), (1.02 * md.R, 0.9 * md.Re)]:
            if r1 <= r0:
                continue
            r = rng.uniform(r0, r1, 200)
            th = rng.uniform(0.02, 0.98, 200) * md.Theta
            x, y = -r * np.cos(th), r * np.sin(th)
            keep = md.in_soil(x, y)
            div, sc = _numerical_div(md, x[keep], y[keep])
            worst.append(np.max(np.abs(div) / sc))
        jumps = []
        for rr in (md.Rw, md.R):
            th = np.linspace(0.01, 0.99, 50) * md.Theta
            vin = md.polar_velocity(rr * (1 - 1e-10), th)
            vout = md.polar_velocity(rr * (1 + 1e-10), th)
            vscale = np.max(np.abs(vout[1])) + np.max(np.abs(vin[1]))
            jumps.append(np.max(np.abs(vin[0] - vout[0])) / vscale)
        th = np.linspace(0.01, 0.99, 50) * md.Theta
        vin = md.polar_velocity(md.Re * (1 - 1e-12), th)
        jumps.append(np.max(np.abs(vin[0])))
        print(f"  beta={b:3d} alpha={a:3d} hw/H={hwr:4.2f} m={md.m:.4f}: max|div v|/(|v|/r) per zone = "
              + " ".join(f"{w:.1e}" for w in worst)
              + f";  [v.e_r] jump at R_w, R (rel) = {jumps[0]:.1e}, {jumps[1]:.1e}; v.e_r at R_e- = {jumps[2]:.1e}")

    print()
    print("=" * 100)
    print("5) stationarity: J* by brute-force quadrature at v'_opt and under perturbations of h1, h2, h3, h4")
    print("=" * 100)
    pcases = [(30, 1, 0.5), (60, 4, 1.0), (20, 10, 0.8), (90, 1, 1.0)] if not quick else [(30, 1, 0.5), (90, 1, 1.0)]
    for (b, a, hwr) in pcases:
        md = AnalyticalSeepage(b, H, hwr * H, a, gamma_w=gw)
        J0 = jstar_by_quadrature(md)
        print(f"  beta={b} alpha={a} hw/H={hwr}: m={md.m:.6f}  J*(closed form)={md.Jstar():.10f}  "
              f"J*(quadrature)={J0:.10f}  rel.diff={(J0 - md.Jstar()) / abs(md.Jstar()):.1e}")
        Rw, Re, R = md.Rw, md.Re, md.R
        g0, h20 = md.g, (lambda t: (md.h2(t), md.dh2(t)))
        ap = md.ap

        def gp(eps, k):
            def f(r):
                G, dG = g0(r)
                s = r / Rw
                # delta h1 = ap * (1 - s) * s^k  (vanishes at R_w)  ->  delta g = r delta h1
                return G + eps * ap * r * (1 - s) * s ** k, dG + eps * ap * ((1 - s) * s ** k + s * (k * s ** (k - 1) * (1 - s) - s ** k) if k > 0 else (1 - 2 * s))
            return f

        def h2p(eps, kind):
            def f(t):
                h, dh = h20(t)
                if kind == 0:
                    return h + eps * np.cos(t), dh - eps * np.sin(t)
                return h + eps * (t / md.Theta) ** 2, dh + eps * 2 * t / md.Theta ** 2
            return f

        def h3p(eps):
            return lambda r: md.h3(r) * (1 + eps * (r / R))

        def h4p(eps):
            return lambda r: md.h4(r) * (1 + eps * np.sin(np.pi * (r - R) / (Re - R)))

        tests = [("h1 += (1-s)", lambda e: dict(g=gp(e, 0))),
                 ("h1 += (1-s)s", lambda e: dict(g=gp(e, 1)))]
        if not md.degenerate:
            tests += [("h2 += cos(th)", lambda e: dict(h2=h2p(e, 0))),
                      ("h2 += (th/Th)^2", lambda e: dict(h2=h2p(e, 1)))]
        tests += [("h3 *= 1+e r/R", lambda e: dict(h3=h3p(e))),
                  ("h4 *= 1+e sin", lambda e: dict(h4=h4p(e)))]
        for name, kw in tests:
            if "h3" in name and md.R <= md.Rw:
                continue
            out = []
            for e in (0.1, -0.1, 0.02, -0.02):
                out.append(jstar_by_quadrature(md, **kw(e)) - J0)
            d1 = (out[2] - out[3]) / 0.04            # first-order coefficient (central difference)
            d2 = (out[2] + out[3]) / 0.02 ** 2 / 2  # second-order coefficient
            print(f"     {name:18s} dJ(+-0.1)={out[0]:+.3e},{out[1]:+.3e}  dJ(+-0.02)={out[2]:+.3e},{out[3]:+.3e}"
                  f"  dJ/de={d1:+.1e}  (1/2)d2J/de2={d2:+.3e}")
        if md.degenerate:
            # perturb (h1, h2) jointly inside the class: the in-class optimum for a given m > 0 has
            # J1* = -k gw^2 hw^2 / (4 Phi(m)); Phi(m) > Phi(0) = A shows the m -> 0 limit is optimal.
            print("     joint (h1,h2) family m>0: Phi(m) - A = " +
                  " ".join(f"m={mm:g}:{md.Phi_of_m(mm) - md.A:+.2e}" for mm in (1e-3, 1e-2, 0.1, 0.5, 1.0)))

    # printed Eq. (31): exponent C/D (not sqrt(C/D)) and coefficient C/D in the e_theta term, prefactor 1/(C-D)
    md = AnalyticalSeepage(30, H, H, 1, gamma_w=gw)
    mm, rCD = md.m, md.C / md.D
    K0 = md.kh * md.gamma_w * md.sb * md.h2e
    sgrid = np.linspace(1e-3, 0.999, 999)
    for label, pref in (("as printed, K0/(C-D)", K0 / (md.C - md.D)), ("sign fixed, K0/(D-C)", K0 / (md.D - md.C))):
        def g_pr(r, pref=pref):
            ss = r / md.Rw
            return pref * r * (1 - ss ** (mm - 1)), pref * (1 - rCD * ss ** (rCD - 1))
        G_true = pref * (1 - mm * sgrid ** (mm - 1))          # (r h1)' implied by the printed e_r term
        G_pr = pref * (1 - rCD * sgrid ** (rCD - 1))          # printed e_theta coefficient
        divres = np.max(np.abs(G_pr - G_true) / np.abs(G_pr).max())
        Jpr = jstar_by_quadrature(md, g=g_pr)
        print(f"  printed Eq.31 ({label}) beta=30 alpha=1 hw=H: J* = {Jpr:+.6f} vs J*(opt) = {md.Jstar():+.6f}; "
              f"max |div v| r/|h2'| / max|v_th| = {divres:.2e} (non-zero -> not admissible)")
    print(f"  face inflow (v.n < 0) of v'_opt for r < s* R_w with s* = m^(1/(1-m)) = {mm ** (1 / (1 - mm)):.4f} (beta=30, alpha=1)")

    print()
    print("=" * 100)
    print("6) h_w = 0 gives zero force; force pattern of K^-1 v'_opt (|f|/gamma_w and direction)")
    print("=" * 100)
    md0 = AnalyticalSeepage(45, H, 0.0, 2)
    fx, fy = md0.force(np.array([-0.3, 0.2, 1.5]), np.array([0.2, 0.5, 1.2]))
    print(f"  hw=0: max|f| = {max(np.max(np.abs(fx)), np.max(np.abs(fy))):.1e}")
    for (b, a) in ((30, 1), (60, 1), (60, 10)):
        md = AnalyticalSeepage(b, H, H, a, gamma_w=gw)
        print(f"  beta={b}, alpha={a}, hw=H (angle measured from +x towards +y=down; 90 = vertical downward)")
        print(f"   {'x/H':>6} {'y/H':>6} {'zone':>5} {'|f|/gw':>8} {'angle':>7} {'fx/gw':>8} {'fy/gw':>8}")
        pts = [(-1.0, 0.25), (-0.5, 0.5), (-0.2, 0.2), (-0.05, 0.6), (0.1, 0.12 + 0.1 * md.tb if np.isfinite(md.tb) else 0.5),
               (0.5 * md.xtoe, 0.5 * H + 0.05), (0.9 * md.xtoe, 0.95 * H + 0.05), (md.xtoe + 0.2, 1.1 * H),
               (md.xtoe + 1.0, 1.05 * H), (0.0, 1.5 * H), (-2.0, 2.0)]
        for (x, y) in pts:
            if not md.in_soil(x, y):
                continue
            fx, fy = md.force(x, y)
            r = np.hypot(x, y)
            z = 1 if r < md.Rw else (2 if r < md.R else (3 if r < md.Re else 4))
            print(f"   {x:6.2f} {y:6.2f} {z:5d} {np.hypot(fx, fy) / gw:8.4f} {np.degrees(np.arctan2(fy, fx)):7.1f} "
                  f"{float(fx) / gw:8.4f} {float(fy) / gw:8.4f}")
    if plot:
        _plot_force_field(gw)


def _plot_force_field(gw):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return
    fig, axs = plt.subplots(1, 3, figsize=(16, 4.6))
    for ax, (b, a, hwr) in zip(axs, ((45, 1, 1.0), (45, 10, 1.0), (60, 1, 0.5))):
        md = AnalyticalSeepage(b, 1.0, hwr, a, gamma_w=gw)
        X, Y = np.meshgrid(np.linspace(-2.0, md.xtoe + 2.0, 46), np.linspace(0.02, 2.2, 30))
        fx, fy = md.force(X, Y)
        mag = np.hypot(fx, fy)
        sel = mag > 0
        ax.quiver(X[sel], Y[sel], fx[sel] / mag[sel], fy[sel] / mag[sel], np.log10(mag[sel] / gw),
                  angles="xy", scale_units="xy", scale=12, cmap="viridis", width=0.0025)
        ax.plot([-2.0, 0.0, md.xtoe, md.xtoe + 2.0], [0.0, 0.0, 1.0, 1.0], "k-", lw=1)
        th = np.linspace(0, md.Theta, 50)
        for rr in (md.Rw, md.R):
            ax.plot(-rr * np.cos(th), rr * np.sin(th), "k--", lw=0.6)
        ax.axhline(hwr, color="tab:blue", lw=0.6, ls=":")
        ax.set_aspect("equal")
        ax.invert_yaxis()
        ax.set_title(f"K^-1 v'_opt, beta={b}, alpha={a}, hw/H={hwr} (colour log10|f|/gw)", fontsize=9)
        ax.set_xlabel("x/H")
        ax.set_ylabel("y/H (down)")
    fig.tight_layout()
    path = os.path.join(RESULTS_DIR, "analytical_force_field.png")
    fig.savefig(path, dpi=130)
    print(f"  -> written {path}")

    # Fig. 5 comparison plot
    data = os.path.join(RESULTS_DIR, "fig5_analytical_Jstar.csv")
    if os.path.exists(data):
        arr = np.loadtxt(data, delimiter=",", comments="#")
        fig, axs = plt.subplots(2, 2, figsize=(9, 6.5))
        for ax, a in zip(axs.ravel(), (1, 2, 4, 10)):
            s = arr[arr[:, 0] == a]
            ax.plot(s[:, 1], s[:, 3], "k-", lw=2.5, alpha=0.35, label="paper Fig. 5 solid")
            ax.plot(s[:, 1], s[:, 2], "r--", lw=1.2, label="this script")
            ax.set_title(f"alpha = {a}")
            ax.set_xlabel("beta (deg)")
            ax.set_ylabel("-J*(v'_opt)/(k_h H^2 gw^2)")
            ax.grid(alpha=0.3)
        axs[0, 0].legend(fontsize=8)
        fig.tight_layout()
        path = os.path.join(RESULTS_DIR, "fig5_analytical_vs_paper.png")
        fig.savefig(path, dpi=130)
        print(f"  -> written {path}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--quick", action="store_true", help="fewer cases")
    ap.add_argument("--no-plot", action="store_true", help="skip PNG output")
    args = ap.parse_args()
    run_checks(quick=args.quick, plot=not args.no_plot)
    sys.exit(0)
