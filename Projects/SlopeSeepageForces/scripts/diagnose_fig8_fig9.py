#!/usr/bin/env python3
"""Diagnosis of the remaining disagreements of the Python reproduction (reproduce_python.py) with the paper
(Ceron, Cecilio, Linn & Maghous, IJNAMG 2025, Figs. 8 and 9):

  (1) Fig. 8, Israeli panel, beta = 35 deg (c = 6, phi = 32 after the swap), -grad u'_FE curve: ours 19-27 % LOWER
      than the paper for h_w/H >= 0.2;
  (2) Israeli 35 / London 30 K^-1 v'_opt curves at small h_w/H: ours HIGHER by 3-57 %;
  (3) Fig. 9, K^-1 v'_opt, beta = 30-50 deg: ours higher by 1-3 % (all alpha).

Gamma is an upper bound, so "ours lower" means either a better mechanism than the paper's (the paper's optimiser or
class missed it) or a different field/setting, and "ours higher" means the paper found a lower value of ITS
objective (better mechanism, different field, or numerical error of its objective exploited by its optimiser).
Each experiment below tests one hypothesis quantitatively (all in paper coordinates, see the shared spec):

  mechanisms   plots of our optimal mechanisms on the slope with the seepage-force field, the circles R_w, R and
               the recirculation radius s* R_w of the analytical field (results/diagnose/mechanisms_*.png)
  selfsimilar  local (self-similar) mechanisms of the analytical field: H_crit * h_w/H along the curves, ours vs paper
  zones        split of P_u of our optimal K^-1 v'_opt mechanisms by zone (r < 0.05 R_w, < s* R_w, < R_w, < R, < R_e)
  quadconv     convergence of our P_u quadrature at the optimal mechanisms (coarse ... ref)
  grid         brute-force grids over (theta1, theta2, s): (a) for the "ours higher" points, is there any mechanism
               below the paper value in our objective? (b) for the "ours lower" Israeli 35 FE points, which fraction of
               the admissible mechanisms beats the paper value (what the paper's optimiser would have had to miss)
  restrict     Israeli 35 FE vs the London 30 / Israeli 60 FE controls: minimum of Gamma over restricted mechanism
               sets (toe only, L >= L_min, mechanism II with d >= d_min, theta2 <= bound, r_h <= bound)
  eq56         mechanism class of Eq. 56 only (no "spiral below O" check, signed integrals) re-optimised
  fe_settings  Israeli 35 FE: alpha, far-side conditions, box size, box fixed in metres with H_ref = 10 m, coarse P1
               meshes, gamma_w, FE field of another slope angle, other slope angle (with controls)
  soil         (c, phi) family that keeps the h_w = 0 end (229.3 m) of the Israeli 35 panel: both curves
  loading      gamma' -> gamma (no buoyancy), saturated weight above the lowered water level
  fe_scale     uniform factor k on the FE seepage force of the Israeli 35 panel that reproduces the paper
  vopt_field   analytical-field variants: m of one fixed-point step from a perturbed m, degenerate m -> 0 branch for
               every beta, L_m (R_e), zone 3 removed, gamma_w = 10
  gamma_w      gamma_w = 9.8 / 9.81 / 10 for Fig. 9 (both curves) and in the Fig. 8 seepage field
  quadrature   emulation of a fixed-resolution numerical P_u (cell-midpoint grid of spacing h = H/25, H/50 in absolute
               coordinates, P_mr and P_gamma closed form) minimised by the same PSO: reported (objective) minimum
  digitization Fig. 9 residual in absolute Gamma vs Gamma (a constant stroke-centre offset would give a constant
               absolute residual); Fig. 8 vector accuracy from data/digitize_meta.json
  summary      ranked list of explanations with the numbers (also written to results/diagnose/diagnose_output.txt)

Runs are cached in results/diagnose/cache.jsonl (one JSON line per finished task, keyed by the task spec), so an
interrupted run resumes where it stopped; tables are rebuilt from the cache.

    python3 Projects/SlopeSeepageForces/scripts/diagnose_fig8_fig9.py [--only a,b,...] [--workers 2] [--plot-only]
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")          # shared CPUs: one thread per worker process

import argparse
import csv
import json
import multiprocessing as mp
import sys
import time

import numpy as np

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.dirname(SCRIPT_DIR)
DATA_DIR = os.path.join(PROJECT_DIR, "data")
OUT_DIR = os.path.join(PROJECT_DIR, "results", "diagnose")
REPRO_DIR = os.path.join(PROJECT_DIR, "results", "reproduce_python")
if SCRIPT_DIR not in sys.path:
    sys.path.insert(0, SCRIPT_DIR)

from analytical_seepage import AnalyticalSeepage, _cd  # noqa: E402
from fe_seepage import FESeepage, MeshSize  # noqa: E402
from limit_analysis import (BIG, Mechanisms, Problem, _gauss01, domain_power, field_circles, get_quad,  # noqa: E402
                            nelder_mead, pso, stability_factor)
from reproduce_python import interp_polyline, paper_fig8_polylines, paper_fig9_polylines  # noqa: E402

# ======================================================================================================
# cases (the settings of reproduce_python.py: Fig. 8 'fitted' soils, gamma_w = 9.8, H_ref = 10 m, FE box 50/10/30 H;
# Fig. 9 H = 5 m, gamma_w = 9.81, FE box fixed in metres 50/10/30 m = 10/2/6 H)
# ======================================================================================================
SOIL8 = {"London": dict(c=11.7, phi=24.7), "Israeli": dict(c=6.0, phi=32.0)}
PANELS8 = (("London", 30.0), ("London", 60.0), ("Israeli", 35.0), ("Israeli", 60.0))
H8 = 10.0
H9 = 5.0
BOX9 = dict(left=10.0, right=2.0, depth=6.0)
EXPERIMENTS = ("mechanisms", "selfsimilar", "zones", "quadconv", "grid", "restrict", "eq56", "fe_settings", "soil",
               "loading", "fe_scale", "vopt_field", "gamma_w", "quadrature", "digitization", "summary")


def case8(soil, beta_, hw_, approach_, **kw):
    """task spec of a Fig. 8 point (approach: vopt | FE | none); kw overrides any entry (e.g. beta, c, phi)"""
    d = dict(fig=8, soil=soil, beta=float(beta_), H=H8, hw=float(hw_), alpha=1.0, gamma=18.0, gamma_w=9.8,
             approach=approach_, **SOIL8[soil])
    d.update(kw)
    return d


def case9(alpha_, beta_, approach_, **kw):
    d = dict(fig=9, beta=float(beta_), H=H9, hw=1.0, alpha=float(alpha_), gamma=20.0, gamma_w=9.81, c=10.0, phi=30.0,
             approach=approach_)
    if approach_ == "FE":
        d["fe"] = dict(BOX9)
    d.update(kw)
    return d


# ======================================================================================================
# field variants
# ======================================================================================================
class MPerturbedSeepage(AnalyticalSeepage):
    """K^-1 v' of the class (29) for h2 = solution of Eq. 37 with a PERTURBED parameter m1 = m_opt (1 + dm), and the
    optimal zone-1 radial function for that h2 (exponent m2 = sqrt(C(m1)/D(m1)), amplitude k_h gw sin(beta) h2e /
    (D (1 + m2))): what one fixed-point step m1 -> sqrt(C/D) from a wrong m1 produces.  Admissible (div-free,
    same boundary fluxes), J* slightly above the optimum.  dm = 0 reproduces the optimal field."""

    def __init__(self, *args, dm=0.0, m1=None, **kw):
        self.dm = float(dm)
        self.m1_abs = m1          # absolute starting value m1 (overrides dm), e.g. 1 = one step from m = 1
        super().__init__(*args, **kw)

    def _solve_h2(self):
        super()._solve_h2()
        if self.degenerate or (self.dm == 0.0 and self.m1_abs is None):
            return
        from scipy.interpolate import CubicHermiteSpline
        m1 = self.m * (1.0 + self.dm) if self.m1_abs is None else float(self.m1_abs)
        th = np.linspace(0.0, self.Theta, self.n_theta)
        sol = self._shoot(m1, t_eval=th, method="DOP853")
        h, q, Cc, Dc = sol.y
        c, d = _cd(th, self.alpha)
        dh = q / c
        cp = (self.alpha - 1.0) * np.sin(2.0 * th)
        d2h = (m1 * d * h - cp * dh) / c
        self._h2 = CubicHermiteSpline(th, h, dh)
        self._dh2 = CubicHermiteSpline(th, dh, d2h)
        self.C, self.D, self.h2e = Cc[-1], Dc[-1], h[-1]
        self.m1 = m1
        self.m = float(np.sqrt(self.C / self.D))
        self.F = self.h2e ** 2 / (np.sqrt(self.C) + np.sqrt(self.D)) ** 2
        self.ap = self.kh * self.gamma_w * self.sb * self.h2e / (self.D * (1.0 + self.m))


class DegenerateSeepage(AnalyticalSeepage):
    """the degenerate m -> 0 optimum (zone 1 = k_h gw sin(beta) / A e_theta) used for EVERY beta (hypothesis: the
    paper's h2 solver always ended on that branch)"""

    def _solve_h2(self):
        self.degenerate = True
        self.m = 0.0
        self.C, self.D, self.h2e = 0.0, self.A, 1.0
        self._h2 = lambda t: np.ones_like(np.asarray(t, float))
        self._dh2 = lambda t: np.zeros_like(np.asarray(t, float))
        self.Phi = self.A
        self.F = 1.0 / self.A
        self.ap = self.kh * self.gamma_w * self.sb * self.h2e / (self.D * (1.0 + self.m))


class NoZone3Seepage(AnalyticalSeepage):
    """analytical field with f = 0 in zone 3 (r >= R)"""

    def polar_velocity(self, r, th):
        vr, vt = super().polar_velocity(r, th)
        z3 = np.asarray(r) >= self.R
        return np.where(z3, 0.0, vr), np.where(z3, 0.0, vt)

    def quadrature_circles(self):
        return [(0.0, 0.0, rr, True) for rr in sorted({self.Rw, self.R}) if rr > 0]

    def Jstar_parts(self):
        j1, j2, _ = super().Jstar_parts()
        return j1, j2, 0.0


class ScaledField:
    """k * f of a potential field (u scaled too, so the boundary formula of P_u stays valid)"""

    force_is_minus_grad_u = True

    def __init__(self, base, k):
        self.base, self.k = base, float(k)
        self.hw = getattr(base, "hw", None)

    def force(self, x, y):
        fx, fy = self.base.force(x, y)
        return self.k * fx, self.k * fy

    def u(self, x, y):
        return self.k * self.base.u(x, y)


class SaturatedAboveWater:
    """FE seepage force plus gw e_y above the lowered water level (y < h_w): the weight above the water level taken
    as the saturated gamma instead of gamma' (domain quadrature only)"""

    def __init__(self, base, hw, gw):
        self.base, self.hw, self.gw = base, float(hw), float(gw)

    def force(self, x, y):
        fx, fy = self.base.force(x, y)
        y = np.asarray(y, float)
        inside = (np.asarray(fx) != 0.0) | (np.asarray(fy) != 0.0)
        return fx, fy + np.where(inside & (y < self.hw), self.gw, 0.0)


def make_field(spec):
    """seepage field of a task spec (None: no seepage)"""
    ap = spec["approach"]
    H, hw = spec["H"], spec["hw"] * spec["H"]
    if ap == "none" or hw <= 0.0:
        return None
    gwf = spec.get("gamma_w_field", spec["gamma_w"])
    beta_f = spec.get("beta_field", spec["beta"])
    if ap == "vopt":
        v = spec.get("vopt", {})
        kw = dict(gamma_w=gwf, Lm_over_H=v.get("Lm_over_H", 10.0))
        var = v.get("variant", "opt")
        if var == "opt":
            return AnalyticalSeepage(beta_f, H, hw, spec["alpha"], **kw)
        if var == "dm":
            return MPerturbedSeepage(beta_f, H, hw, spec["alpha"], dm=v.get("dm", 0.0), m1=v.get("m1"), **kw)
        if var == "degenerate":
            return DegenerateSeepage(beta_f, H, hw, spec["alpha"], **kw)
        if var == "nozone3":
            return NoZone3Seepage(beta_f, H, hw, spec["alpha"], **kw)
        raise ValueError(var)
    if ap == "FE":
        fe = dict(spec.get("fe", {}))
        size = fe.pop("size", None)
        if size is not None:
            fe["size"] = MeshSize(**size)
        f = FESeepage(beta_f, H=H, hw=hw, alpha=spec.get("alpha_field", spec["alpha"]), gamma_w=gwf, **fe)
        if "fe_scale" in spec:
            f = ScaledField(f, spec["fe_scale"])
        if spec.get("loading") == "sat_above_water":
            f = SaturatedAboveWater(f, hw, spec["gamma_w"])
        return f
    raise ValueError(ap)


# ======================================================================================================
# objectives: exact (limit_analysis) with optional restrictions, emulated grid quadrature, Eq. 56-only class
# ======================================================================================================
def grid_pu(mech, field, h):
    """P_u / omega by a cell-midpoint rule on the fixed grid x, y in (k + 1/2) h (absolute coordinates, origin O):
    cells whose centre lies in the soil and inside the mechanism (between the ground surface and the spiral)"""
    out = np.zeros(len(mech))
    tanb = np.tan(mech.beta) if mech.beta_deg < 90.0 else np.inf
    for i in range(len(mech)):
        th = np.linspace(mech.t1[i], mech.t2[i], 200)
        r = mech.r0[i] * np.exp((th - mech.t1[i]) * mech.k)
        xs, ys = mech.Cx[i] - r * np.cos(th), -mech.Cy[i] + r * np.sin(th)
        x0, x1, y1 = min(xs.min(), mech.Ax[i]), max(xs.max(), mech.Bx[i]), max(ys.max(), mech.By[i])
        ix = np.arange(np.floor(x0 / h) - 1, np.ceil(x1 / h) + 1)
        iy = np.arange(-1, np.ceil(y1 / h) + 1)
        X, Y = np.meshgrid((ix + 0.5) * h, (iy + 0.5) * h)
        X, Y = X.ravel(), Y.ravel()
        dx, dy = mech.Cx[i] - X, Y + mech.Cy[i]
        ang = np.arctan2(dy, dx)
        rsp = mech.r0[i] * np.exp((ang - mech.t1[i]) * mech.k)
        with np.errstate(invalid="ignore"):
            soil = np.where(X <= 0.0, Y >= 0.0, np.where(X < mech.xtoe, Y >= X * tanb, Y >= mech.H))
        ins = soil & (ang >= mech.t1[i]) & (ang <= mech.t2[i]) & (np.hypot(dx, dy) < rsp)
        if not np.any(ins):
            continue
        X, Y = X[ins], Y[ins]
        fx, fy = field.force(X, Y)
        out[i] = np.sum(fx * (Y + mech.Cy[i]) + fy * (mech.Cx[i] - X)) * h * h
    return out


def eq56_ok(m):
    """admissibility of Eq. 56 only (mechanism I): 0 < eta <= 1, 0 < theta1 < theta2 < pi - beta, L > 0"""
    with np.errstate(all="ignore"):
        ok = (m.t1 > 0) & (m.t2 > m.t1) & (m.t2 < np.pi - m.beta) & (m.s > 0) & (m.s <= 1.0)
        ok &= np.isfinite(m.r0) & (m.r0 > 0) & (m.rh <= m.R_MAX * m.H) & (m.L > 0)
    return ok


def signed_power(m, force, circles, npan=16, q=8, qr=8, nr=2):
    """signed polar integral of f.U between the ground surface (A-O-B) and the spiral (valid also when the spiral
    passes above the ground surface; the region between them then counts negatively, as in f1 - f2 - f3)"""
    out = np.zeros(len(m))
    x01, w01 = _gauss01(q)
    u01, wu01 = _gauss01(qr)
    for i in range(len(m)):
        t1, t2, thO = m.t1[i], m.t2[i], m.thO[i]
        br = set(np.linspace(t1, thO, npan // 2 + 1)) | set(np.linspace(thO, t2, npan + 1))
        for a, b in ((t1, thO), (thO, t2)):
            for k in range(1, 6):
                br.add(a + (b - a) * 0.25 ** k)
                br.add(b - (b - a) * 0.25 ** k)
        ts = np.linspace(t1, t2, 2001)
        g = m.r0[i] * np.exp((ts - t1) * m.k) - np.where(ts < thO, m.Cy[i] / np.sin(ts), m.D[i] / np.sin(ts + m.beta))
        for j in np.nonzero(np.signbit(g[:-1]) != np.signbit(g[1:]))[0]:
            br.add(0.5 * (ts[j] + ts[j + 1]))
        br = np.array(sorted(br))
        TH = (br[:-1, None] + np.diff(br)[:, None] * x01).ravel()
        WT = (np.diff(br)[:, None] * w01).ravel()
        RS = np.where(TH < thO, m.Cy[i] / np.sin(TH), m.D[i] / np.sin(TH + m.beta))
        RP = m.r0[i] * np.exp((TH - t1) * m.k)
        ub = [np.zeros_like(TH), np.ones_like(TH)] + [np.full_like(TH, j / nr) for j in range(1, nr)]
        cT, sT = np.cos(TH), np.sin(TH)
        for c in circles:
            qx, qy, R = c[:3]
            ddx, ddy = m.Cx[i] - qx, m.Cy[i] + qy
            bb = ddx * cT + ddy * sT
            disc = bb * bb - (ddx * ddx + ddy * ddy - R * R)
            sq = np.sqrt(np.maximum(disc, 0))
            for root in (bb - sq, bb + sq):
                with np.errstate(all="ignore"):
                    uu = (root - RS) / (RP - RS)
                ub.append(np.where((disc > 0) & np.isfinite(uu), np.clip(uu, 0, 1), 0.0))
        UB = np.sort(np.stack(ub, -1), -1)
        UA, UL = UB[:, :-1], np.diff(UB, axis=-1)
        RHO = RS[:, None, None] + (RP - RS)[:, None, None] * (UA[..., None] + UL[..., None] * u01)
        W = WT[:, None, None] * (RP - RS)[:, None, None] * UL[..., None] * wu01
        fx, fy = force(m.Cx[i] - RHO * cT[:, None, None], -m.Cy[i] + RHO * sT[:, None, None])
        out[i] = np.sum((fx * sT[:, None, None] + fy * cT[:, None, None]) * RHO * RHO * W)
    return out


class Objective:
    """Gamma(X) for rows X = (theta1, theta2, s) of one mechanism class with optional restrictions and P_u methods.
    la options: quad = 'exact' (limit_analysis polar quadrature / boundary formula) | 'cellgrid' (grid_pu with
    h = H / n_cell) | 'eq56' (Eq. 56-only class, signed integrals); restrict = dict(L_min, d_min, t2_max, rh_max)
    in units of H / rad."""

    def __init__(self, prob, kind, la):
        self.p, self.kind, self.la = prob, kind, la
        self.quad = la.get("quad", "exact")
        self.rs = la.get("restrict", {})
        self.circles = field_circles(prob.field, 4) if prob.field is not None else []

    def ok(self, X):
        m = self.p.mechanisms(X)
        ok = eq56_ok(m) if self.quad == "eq56" else m.ok
        ok = ok & ((X[:, 2] <= 1.0) if self.kind == "I" else (X[:, 2] >= 1.0))
        H = self.p.H
        if "L_min" in self.rs:
            ok &= m.L >= self.rs["L_min"] * H
        if "d_min" in self.rs:
            ok &= (X[:, 2] - 1.0) >= self.rs["d_min"]
        if "t2_max" in self.rs:
            ok &= X[:, 1] <= self.rs["t2_max"]
        if "rh_max" in self.rs:
            ok &= m.rh <= self.rs["rh_max"] * H
        return ok, m

    def __call__(self, X, level="coarse"):
        X = np.atleast_2d(np.asarray(X, float))
        ok, m = self.ok(X)
        out = np.full(len(X), BIG)
        idx = np.nonzero(ok)[0]
        if len(idx) == 0:
            return out
        sub = m.subset(idx)
        p = self.p
        if self.quad == "exact":
            Pmr, Pg, Pu = p.powers(sub, level)
        else:
            Pmr, Pg = sub.Pmr(p.c), sub.Pgamma(p.gamma_p)
            if p.field is None:
                Pu = np.zeros(len(idx))
            elif self.quad == "cellgrid":
                Pu = grid_pu(sub, p.field, p.H / self.la["n_cell"])
            else:
                Pu = signed_power(sub, p.field.force, self.circles)
        Pe = Pg + Pu
        with np.errstate(all="ignore"):
            G = np.where(Pe > 0.0, Pmr / Pe, BIG)
        out[idx] = np.where(np.isfinite(G), G, BIG)
        return out


def optimise(prob, la):
    """min Gamma over the classes of la['classes'] (default I and II) by PSO (search level) + Nelder-Mead polish
    (final level), as limit_analysis.stability_factor; with a non-exact quadrature both use the same objective.
    Returns dict(Gamma, x, kind, Gamma_exact (exact objective at the found mechanism), seed_spread)."""
    seeds = tuple(la.get("seeds", (0, 1)))
    npart, nit = la.get("n_particles", 40), la.get("n_iter", 150)
    best = None
    vals = []
    for kind in la.get("classes", ("I", "II")):
        obj = Objective(prob, kind, la)
        lb, ub = prob.bounds(kind)
        if kind == "II" and "d_max" in la:
            ub = ub.copy()
            ub[2] = 1.0 + la["d_max"]
        if la.get("toe_only"):
            lb, ub = lb.copy(), ub.copy()
            lb[2] = ub[2] = 1.0
        fixed = lb[2] == ub[2]
        for sd in seeds:
            if fixed:
                l2, u2 = lb[:2], ub[:2]

                def f2(Z, obj=obj):
                    Z = np.atleast_2d(Z)
                    return obj(np.column_stack([Z, np.full(len(Z), lb[2])]))

                def feas2(Z, obj=obj):
                    return obj.ok(np.column_stack([Z, np.full(len(Z), lb[2])]))[0]

                z, fv, _ = pso(f2, l2, u2, npart, nit, sd, feas2)
                if z is None or fv >= BIG:
                    continue
                fin = "fine" if obj.quad == "exact" else None
                zz, fv = nelder_mead(lambda w: obj(np.r_[w, lb[2]], fin or "coarse")[0], z, l2, u2, step=(0.02, 0.02))
                x = np.r_[zz, lb[2]]
            else:
                x, fv, _ = pso(obj, lb, ub, npart, nit, sd, lambda X, obj=obj: obj.ok(X)[0])
                if x is None or fv >= BIG:
                    continue
                fin = "fine" if obj.quad == "exact" else "coarse"
                x, fv = nelder_mead(lambda w: obj(w, fin)[0], x, lb, ub)
            vals.append(fv)
            if best is None or fv < best[0]:
                best = (fv, x, kind)
    if best is None:
        return dict(Gamma=np.inf, x=None, kind=None, Gamma_exact=np.inf, seed_spread=np.nan)
    G, x, kind = best
    gex = float(Objective(prob, kind, dict(quad="exact"))(x, "fine")[0])
    m = prob.mechanisms(x)
    return dict(Gamma=float(G), x=[float(v) for v in x], kind=kind, Gamma_exact=gex,
                L_over_H=float(m.L[0] / prob.H), seed_spread=float(max(vals) / min(vals) - 1.0))


# ======================================================================================================
# tasks
# ======================================================================================================
def _problem(spec, field):
    gw_gamma = spec.get("gamma_w_weight", spec["gamma_w"])
    return Problem(spec["beta"], spec["H"], spec["c"], spec["phi"], spec["gamma"], gw_gamma, field,
                   spec.get("la", {}).get("pu_method", "auto"))


def run_task(spec):
    t0 = time.time()
    field = make_field(spec)
    t_field = time.time() - t0
    task = spec.get("task", "opt")
    out = dict(spec=spec)
    if task == "opt":
        la = spec.get("la", {})
        if not la or set(la) <= {"seeds", "classes", "d_max", "n_particles", "n_iter", "pu_method"}:
            r = stability_factor(spec["beta"], spec["H"], spec["c"], spec["phi"], spec["gamma"],
                                 spec.get("gamma_w_weight", spec["gamma_w"]), field,
                                 mechanisms=tuple(la.get("classes", ("I", "II"))), seeds=tuple(la.get("seeds", (0, 1))),
                                 n_particles=la.get("n_particles", 40), n_iter=la.get("n_iter", 150),
                                 d_max=la.get("d_max", 10.0), pu_method=la.get("pu_method", "auto"))
            out.update(Gamma=r["Gamma"], Gamma_exact=r["Gamma"], x=r.get("x"), kind=r["mechanism"],
                       L_over_H=(r["L"] / spec["H"]) if r["mechanism"] else np.nan, seed_spread=r.get("seed_spread"))
        else:
            out.update(optimise(_problem(spec, field), la))
        out["Hcrit"] = out["Gamma"] * spec["H"]
    elif task == "eval":
        prob = _problem(spec, field)
        X = np.array(spec["X"], float)
        m = prob.mechanisms(X)
        Pmr, Pg, Pu = prob.powers(m, "fine")
        out.update(Gamma=[float(v) for v in Pmr / (Pg + Pu)], P_mr=Pmr.tolist(), P_gamma=Pg.tolist(), P_u=Pu.tolist())
    elif task == "grid":
        out.update(grid_task(spec, field))
    elif task == "field":
        out.update(field_info(field))
    else:
        raise ValueError(task)
    out.update(t_field=t_field, time=time.time() - t0)
    if isinstance(field, AnalyticalSeepage):
        out["m"] = float(field.m)
        out["Jn"] = float(-field.Jstar_normalized())
    return _jsonable(out)


def field_info(field):
    if field is None:
        return {}
    if isinstance(field, AnalyticalSeepage):
        return dict(m=float(field.m), Jn=float(-field.Jstar_normalized()), F=float(field.F))
    return dict(Jn=float(field.J_normalized())) if hasattr(field, "J_normalized") else {}


def grid_task(spec, field):
    """brute-force grid of class I (log-spaced eta) and II: min Gamma, count and fraction below a threshold, and the
    minimum over restricted subsets (restrict experiment)"""
    prob = _problem(spec, field)
    g = spec["grid"]
    thr = g.get("threshold", np.nan)
    H = spec["H"]
    res = dict(n_adm=0, n_below=0, min=np.inf, x_min=None, subsets={})
    sub_defs = g.get("subsets", {})
    sub_min = {k: np.inf for k in sub_defs}
    for kind in ("I", "II"):
        lb, ub = prob.bounds(kind)
        g1 = np.linspace(lb[0], ub[0], g.get("n1", 100))
        g2 = np.linspace(lb[1], ub[1], g.get("n2", 100))
        if kind == "I":
            g3 = np.exp(np.linspace(np.log(g.get("eta_min", 0.01)), 0.0, g.get("n3", 40)))
        else:
            g3 = 1.0 + np.linspace(0.0, g.get("d_max", 3.0), g.get("n3ii", 30))
        G = np.stack(np.meshgrid(g1, g2, g3, indexing="ij"), -1).reshape(-1, 3)
        m = prob.mechanisms(G)
        G = G[m.ok]
        if len(G) == 0:
            continue
        f = np.concatenate([prob.gamma_factor(G[i:i + 20000], kind, g.get("quad", "coarse"))
                            for i in range(0, len(G), 20000)])
        fin = f < BIG / 10
        res["n_adm"] += int(fin.sum())
        res["n_below"] += int((f < thr).sum()) if np.isfinite(thr) else 0
        j = int(np.argmin(f))
        if f[j] < res["min"]:
            res["min"], res["x_min"] = float(f[j]), G[j].tolist()
        mm = prob.mechanisms(G)
        for name, sd in sub_defs.items():
            sel = fin.copy()
            if sd.get("toe_only"):
                sel &= np.abs(G[:, 2] - 1.0) < 1e-12
            if "L_min" in sd:
                sel &= mm.L >= sd["L_min"] * H
            if "d_min" in sd:
                sel &= (G[:, 2] - 1.0) >= sd["d_min"]
            if "t2_max" in sd:
                sel &= G[:, 1] <= sd["t2_max"]
            if "rh_max" in sd:
                sel &= mm.rh <= sd["rh_max"] * H
            if "class" in sd:
                sel &= kind == sd["class"]
            if np.any(sel):
                sub_min[name] = min(sub_min[name], float(f[sel].min()))
    # polish the grid minimum with the exact objective (Nelder-Mead, class of the grid minimum)
    if res["x_min"] is not None:
        kind = "II" if res["x_min"][2] > 1.0 else "I"
        lb, ub = prob.bounds(kind)
        x, fv = nelder_mead(lambda w: prob.gamma_factor(w, kind, "fine")[0], np.array(res["x_min"]), lb, ub)
        res["min_polished"], res["x_polished"] = float(fv), x.tolist()
    res["subsets"] = sub_min
    return res


def _jsonable(v):
    if isinstance(v, dict):
        return {str(k): _jsonable(x) for k, x in v.items()}
    if isinstance(v, (list, tuple)):
        return [_jsonable(x) for x in v]
    if isinstance(v, (np.floating, float)):
        return float(v)
    if isinstance(v, np.integer):
        return int(v)
    if isinstance(v, np.bool_):
        return bool(v)
    if isinstance(v, np.ndarray):
        return [_jsonable(x) for x in v.tolist()]
    return v


def key(spec):
    return json.dumps(spec, sort_keys=True)


class Cache:
    def __init__(self, path):
        self.path = path
        self.data = {}
        self.refresh()

    def refresh(self):
        if os.path.exists(self.path):
            with open(self.path) as fh:
                for line in fh:
                    try:
                        d = json.loads(line)
                    except json.JSONDecodeError:
                        continue
                    self.data[key(d["spec"])] = d

    def get(self, spec):
        return self.data.get(key(spec))

    def put(self, res):
        self.data[key(res["spec"])] = res
        with open(self.path, "a") as fh:
            fh.write(json.dumps(res) + "\n")


def run_all(specs, args, cache, log, label):
    """run the specs missing from the cache (in parallel, results appended as they finish) and return all results"""
    todo = []
    seen = set()
    for s in specs:
        k = key(s)
        if cache.get(s) is None and k not in seen:
            todo.append(s)
            seen.add(k)
    if todo and not args.plot_only:
        log(f"[{label}] {len(specs)} tasks, {len(todo)} to run on {args.workers} worker(s)")
        t0 = time.time()
        if args.workers <= 1:
            for i, r in enumerate(map(run_task, todo)):
                cache.put(r)
                _progress(log, label, i + 1, len(todo), r, t0)
        else:
            with mp.get_context("fork").Pool(args.workers, maxtasksperchild=20) as pool:
                for i, r in enumerate(pool.imap_unordered(run_task, todo, chunksize=1)):
                    cache.put(r)
                    _progress(log, label, i + 1, len(todo), r, t0)
    return [cache.get(s) for s in specs]


def _progress(log, label, i, n, r, t0):
    s = r["spec"]
    v = r.get("Hcrit", r.get("min", r.get("Gamma", "")))
    v = f"{v:.4f}" if isinstance(v, float) else str(v)[:40]
    log(f"   [{label} {i}/{n} {time.time() - t0:6.0f} s] fig{s['fig']} {s.get('soil', '')} beta={s['beta']:g} "
        f"hw/H={s['hw']:g} alpha={s['alpha']:g} {s['approach']} {s.get('task', 'opt')}: {v} ({r['time']:.1f} s)")


# ======================================================================================================
# paper values and our reproduce_python baseline
# ======================================================================================================
P8 = paper_fig8_polylines()
P9 = paper_fig9_polylines()


def paper8(soil, beta, hw, curve):
    return interp_polyline(P8[(soil, float(beta), curve)], hw, log=True)[0]


def paper9(alpha, beta, curve):
    return interp_polyline(P9[(int(alpha), curve)], beta)[0]


def baseline8():
    """{(soil, beta, curve, hw): row} of results/reproduce_python/comparison_fig8_fitted.csv"""
    out = {}
    with open(os.path.join(REPRO_DIR, "comparison_fig8_fitted.csv")) as fh:
        for r in csv.DictReader(fh):
            out[(r["soil"], float(r["beta_deg"]), r["curve"], float(r["hw_over_H"]))] = r
    return out


def baseline9():
    out = {}
    with open(os.path.join(REPRO_DIR, "comparison_fig9.csv")) as fh:
        for r in csv.DictReader(fh):
            if r["curve"] == "vopt":
                out[(int(r["alpha"]), float(r["beta_deg"]), "vopt")] = r
    with open(os.path.join(REPRO_DIR, "diag_fig9_fe_box_metres.csv")) as fh:
        for r in csv.DictReader(fh):
            out[(int(r["alpha"]), float(r["beta_deg"]), "FE")] = dict(ours=r["Gamma_box_metres"], paper=r["paper"])
    return out


def x_of(r):
    s = float(r["eta_or_dH"]) if r["mechanism"] == "I" else 1.0 + float(r["eta_or_dH"])
    return [float(r["theta1"]), float(r["theta2"]), s]


def _rel(a, b):
    return a / b - 1.0 if (np.isfinite(a) and np.isfinite(b) and b != 0) else np.nan


def write_csv(name, rows):
    if not rows:
        return
    path = os.path.join(OUT_DIR, name)
    fields = list(rows[0].keys())
    for r in rows[1:]:
        for k in r:
            if k not in fields:
                fields.append(k)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.6g}" if isinstance(v, float) else v) for k, v in r.items()})


# ======================================================================================================
# experiments
# ======================================================================================================
def exp_mechanisms(args, cache, log):
    """plots of the optimal mechanisms (from the reproduce_python baseline) with the seepage-force fields"""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    b8 = baseline8()
    col = {"vopt": "#2a78d6", "FE": "#eb6834", "slope": "#0b0b0b", "field": "#a9a8a2", "circ": "#52514e"}

    def draw(ax, beta, H, hw, field, prob, x, title, zoom=None):
        m = prob.mechanisms(np.array(x))
        xt = H / np.tan(np.radians(beta))
        ext = zoom or (-0.8 * H, xt + 0.6 * H, -0.5 * H, 1.35 * H)
        ax.plot([ext[0], 0, xt, ext[1]], [0, 0, H, H], "-", color=col["slope"], lw=1.0)
        n = 30
        X, Y = np.meshgrid(np.linspace(ext[0], ext[1], n), np.linspace(max(ext[2], 0.0), ext[3], int(n * 0.8)))
        fx, fy = field.force(X, Y) if field is not None else (np.zeros_like(X), np.zeros_like(Y))
        mag = np.hypot(fx, fy)
        sc = np.where(mag > 0, np.minimum(1.0, 1.2 * prob.gamma_w / np.maximum(mag, 1e-30)), 0.0)   # clip the arrows
        # arrow length H/6 for |f| = gamma_w (clipped at 1.2 gamma_w); same relative length in the zoomed plots
        ax.quiver(X, Y, fx * sc, fy * sc, angles="xy", color=col["field"], width=0.0022, scale_units="xy",
                  scale=prob.gamma_w * 6.0 / H * (1.0 if zoom is None else 2.4 * H / (ext[1] - ext[0])))
        if isinstance(field, AnalyticalSeepage) and field.hw > 0:
            t = np.linspace(0.0, np.pi - np.radians(beta), 120)
            radii = [(field.Rw, "--", "R_w"), (field.R, ":", "R")]
            if 0 < field.m < 1:
                radii.append((field.Rw * field.m ** (1 / (1 - field.m)), "-.", "s* R_w"))
            for R, ls, lab in radii:
                ax.plot(-R * np.cos(t), R * np.sin(t), ls, color=col["circ"], lw=0.8, label=lab)
        th = np.linspace(m.t1[0], m.t2[0], 400)
        r = m.r0[0] * np.exp((th - m.t1[0]) * m.k)
        ax.plot(m.Cx[0] - r * np.cos(th), -m.Cy[0] + r * np.sin(th), "-", color="#d03b3b", lw=1.8)
        if hw > 0:
            ax.plot([hw / np.tan(np.radians(beta))], [hw], marker="v", color="#2a78d6", ms=6, ls="none")
        ax.set_aspect("equal")
        ax.set_xlim(ext[0], ext[1])
        ax.set_ylim(ext[3], ext[2])
        ax.set_title(title, fontsize=7.5, loc="left")
        ax.tick_params(labelsize=7)

    sets = [("Israeli", 35.0, (0.1, 0.2, 0.5, 1.0)), ("London", 30.0, (0.05, 0.1, 0.5, 1.0)),
            ("London", 60.0, (0.1, 0.5, 1.0)), ("Israeli", 60.0, (0.1, 0.5, 1.0))]
    for soil, beta, hws in sets:
        fig, axs = plt.subplots(2, len(hws), figsize=(3.6 * len(hws), 6.0), squeeze=False)
        for j, hw in enumerate(hws):
            for i, curve in enumerate(("vopt", "FE")):
                r = b8[(soil, beta, curve, hw)]
                spec = case8(soil, beta, hw, curve)
                field = make_field(spec)
                prob = _problem(spec, field)
                draw(axs[i][j], beta, H8, hw * H8, field, prob, x_of(r),
                     f"{curve}  h_w/H={hw}: H_crit {float(r['ours']):.1f} m (paper {float(r['paper']):.1f})\n"
                     f"{r['mechanism']} s={float(r['eta_or_dH']):.3f}  L/H={float(r['L_over_H']):.3f}")
        fig.suptitle(f"{soil} panel, beta = {beta:.0f} deg (c = {SOIL8[soil]['c']}, phi = {SOIL8[soil]['phi']}): our optimal "
                     "mechanisms (red), seepage force (grey, clipped), water level (triangle); circles R_w, R, s*R_w",
                     fontsize=8.5)
        fig.tight_layout(rect=(0, 0, 1, 0.96))
        path = os.path.join(OUT_DIR, f"mechanisms_{soil.lower()}{beta:.0f}.png")
        fig.savefig(path, dpi=120)
        plt.close(fig)
        log(f"   wrote {path}")
    # Fig. 9: vopt and FE (box in metres) at beta = 35, 45, 60, alpha = 1, 5
    rows9 = {}
    with open(os.path.join(REPRO_DIR, "comparison_fig9.csv")) as fh:
        for r in csv.DictReader(fh):
            rows9[(int(r["alpha"]), float(r["beta_deg"]), r["curve"])] = r
    fig, axs = plt.subplots(2, 3, figsize=(11.0, 6.4), squeeze=False)
    for j, beta in enumerate((35.0, 45.0, 60.0)):
        for i, a in enumerate((1, 5)):
            r = rows9[(a, beta, "vopt")]
            spec = case9(a, beta, "vopt")
            field = make_field(spec)
            prob = _problem(spec, field)
            draw(axs[i][j], beta, H9, H9, field, prob, x_of(r),
                 f"Fig. 9 vopt alpha={a} beta={beta:.0f}: Gamma {float(r['ours']):.4f} (paper {float(r['paper']):.4f})\n"
                 f"L/H = {prob.mechanisms(np.array(x_of(r))).L[0] / H9:.4f}")
    fig.tight_layout()
    path = os.path.join(OUT_DIR, "mechanisms_fig9_vopt.png")
    fig.savefig(path, dpi=120)
    plt.close(fig)
    log(f"   wrote {path}")
    # zoom on O: vopt optimal mechanisms that start next to O
    fig, axs = plt.subplots(1, 3, figsize=(11.0, 3.9), squeeze=False)
    for j, (lab, spec, r) in enumerate((("Israeli 35, h_w/H = 0.5", case8("Israeli", 35.0, 0.5, "vopt"),
                                         b8[("Israeli", 35.0, "vopt", 0.5)]),
                                        ("Israeli 35, h_w/H = 1", case8("Israeli", 35.0, 1.0, "vopt"),
                                         b8[("Israeli", 35.0, "vopt", 1.0)]),
                                        ("Fig. 9 alpha = 5, beta = 35", case9(5, 35.0, "vopt"), rows9[(5, 35.0, "vopt")]))):
        field = make_field(spec)
        prob = _problem(spec, field)
        H = spec["H"]
        R = field.Rw
        draw(axs[0][j], spec["beta"], H, spec["hw"] * H, field, prob, x_of(r),
             f"zoom on O: {lab}\nL/H = {prob.mechanisms(np.array(x_of(r))).L[0] / H:.4f}",
             zoom=(-0.45 * R, 0.6 * R, -0.25 * R, 0.6 * R))
    fig.tight_layout()
    path = os.path.join(OUT_DIR, "mechanisms_zoom_O.png")
    fig.savefig(path, dpi=130)
    plt.close(fig)
    log(f"   wrote {path}")


def exp_selfsimilar(args, cache, log):
    """H_crit * h_w / H along the Fig. 8 curves (constant for a mechanism family that scales with R_w)"""
    b8 = baseline8()
    rows = []
    log("   H_crit * (h_w/H) [m] (constant along a family of mechanisms that scales with R_w = h_w / sin beta);")
    log("   'B/R_w' = distance of the exit point B from O over R_w")
    for soil, beta in PANELS8:
        for curve in ("vopt", "FE"):
            line = []
            for hw in (0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0):
                r = b8.get((soil, beta, curve, hw))
                if r is None:
                    continue
                ours, pap = float(r["ours"]), float(r["paper"])
                s = float(r["eta_or_dH"]) if r["mechanism"] == "I" else 1.0
                b_over_rw = s / hw
                rows.append(dict(soil=soil, beta=beta, curve=curve, hw_over_H=hw, ours_x_hw=ours * hw,
                                 paper_x_hw=pap * hw, rel=_rel(ours, pap), mechanism=r["mechanism"], s=s,
                                 B_over_Rw=b_over_rw, L_over_H=float(r["L_over_H"])))
                line.append(f"{hw:.1f}:{ours * hw:6.2f}/{pap * hw:6.2f}" + ("*" if s < 0.999 else " "))
            log(f"   {soil:7s} {beta:3.0f} {curve:4s} ours/paper  " + " ".join(line))
    log("   (* = mechanism exits on the face above the toe, eta < 1)")
    write_csv("selfsimilar.csv", rows)
    return rows


def exp_zones(args, cache, log):
    """split of P_u of our optimal vopt mechanisms by zone (evaluated here, no optimisation)"""
    b8 = baseline8()
    rows9 = {}
    with open(os.path.join(REPRO_DIR, "comparison_fig9.csv")) as fh:
        for r in csv.DictReader(fh):
            rows9[(int(r["alpha"]), float(r["beta_deg"]), r["curve"])] = r
    cases = [(f"{s} {b:.0f} hw/H={hw}", case8(s, b, hw, "vopt"), b8[(s, b, "vopt", hw)])
             for s, b in PANELS8 for hw in (0.1, 0.2, 0.5, 1.0)]
    cases += [(f"Fig9 a={a} b={b:.0f}", case9(a, b, "vopt"), rows9[(a, b, "vopt")])
              for a in (1, 5, 10) for b in (30.0, 35.0, 40.0, 45.0, 50.0, 60.0, 75.0) if (a, b, "vopt") in rows9]
    qd = get_quad("fine")
    out = []
    log(f"   {'case':22s} {'m':>6s} {'L/H':>7s} | fraction of P_u from r < 0.05R_w, [0.05, s*]R_w, [s*, 1]R_w, "
        f"zone 2, zone 3 | Gamma if f = 0 for r < s*R_w | paper/ours - 1")
    for name, spec, r in cases:
        if r.get("mechanism") in (None, ""):
            continue
        fld = make_field(spec)
        prob = _problem(spec, fld)
        m = prob.mechanisms(np.array(x_of(r)))
        sst = fld.m ** (1.0 / (1.0 - fld.m)) if 0.0 < fld.m < 1.0 else 0.0
        edges = [0.0, 0.05 * fld.Rw, max(sst, 0.05) * fld.Rw, fld.Rw, fld.R, fld.Re]
        circ = field_circles(fld, qd.n_circle_grade) + [(0.0, 0.0, e, True) for e in edges[1:3]]
        parts = []
        for lo, hi in zip(edges[:-1], edges[1:]):
            def fz(x, y, lo=lo, hi=hi):
                fx, fy = fld.force(x, y)
                rr = np.hypot(x, y)
                k = (rr >= lo) & (rr < hi)
                return np.where(k, fx, 0.0), np.where(k, fy, 0.0)
            parts.append(float(domain_power(m, fz, qd, circ)[0]))
        Pmr, Pg, _ = prob.powers(m, "fine")
        tot = sum(parts)
        G = Pmr[0] / (Pg[0] + tot)
        Gc = Pmr[0] / (Pg[0] + tot - parts[0] - parts[1])
        pap = float(r["paper"])
        ours = float(r["ours"])
        row = dict(case=name, m=fld.m, L_over_H=m.L[0] / spec["H"], P_u=tot, P_gamma=float(Pg[0]), P_mr=float(Pmr[0]),
                   frac_core005=parts[0] / tot, frac_core_sstar=parts[1] / tot, frac_zone1_rest=parts[2] / tot,
                   frac_zone2=parts[3] / tot, frac_zone3=parts[4] / tot, Gamma=G,
                   rel_Gamma_without_sstar_core=Gc / G - 1.0, paper_over_ours=pap / ours - 1.0)
        out.append(row)
        log(f"   {name:22s} {fld.m:6.3f} {row['L_over_H']:7.4f} | {parts[0] / tot:+.3f} {parts[1] / tot:+.3f} "
            f"{parts[2] / tot:+.3f} {parts[3] / tot:+.3f} {parts[4] / tot:+.3f} | {Gc / G - 1:+.3f} | {pap / ours - 1:+.4f}")
    write_csv("zones_pu_split.csv", out)


def exp_quadconv(args, cache, log):
    b8 = baseline8()
    rows9 = {}
    with open(os.path.join(REPRO_DIR, "comparison_fig9.csv")) as fh:
        for r in csv.DictReader(fh):
            rows9[(int(r["alpha"]), float(r["beta_deg"]), r["curve"])] = r
    cases = [(f"{s} {b:.0f} hw/H={hw}", case8(s, b, hw, "vopt"), b8[(s, b, "vopt", hw)])
             for s, b in (("Israeli", 35.0), ("London", 30.0)) for hw in (0.1, 0.2, 0.5, 1.0)]
    cases += [(f"Fig9 a={a} b={b:.0f}", case9(a, b, "vopt"), rows9[(a, b, "vopt")]) for a in (1, 5) for b in (35.0, 45.0)]
    out = []
    levels = ("coarse", "medium", "fine", "xfine", "ref")
    log("   Gamma(level) / Gamma(ref) - 1 at our optimal K^-1 v'_opt mechanisms (P_u by polar quadrature)")
    for name, spec, r in cases:
        fld = make_field(spec)
        prob = _problem(spec, fld)
        m = prob.mechanisms(np.array(x_of(r)))
        Pmr, Pg, _ = prob.powers(m, "fine")
        G = []
        for q in levels:
            qd = get_quad(q)
            G.append(float(Pmr[0] / (Pg[0] + domain_power(m, fld.force, qd, field_circles(fld, qd.n_circle_grade))[0])))
        out.append(dict(case=name, **{f"rel_{q}": g / G[-1] - 1 for q, g in zip(levels, G)}, Gamma_ref=G[-1]))
        log(f"   {name:22s} " + " ".join(f"{q}:{g / G[-1] - 1:+.1e}" for q, g in zip(levels, G)))
    write_csv("quadrature_convergence.csv", out)


def exp_grid(args, cache, log):
    """(a) ours-higher points: any mechanism below the paper value? (b) Israeli 35 FE: fraction below the paper"""
    specs = []
    pts = [("London", 30.0, hw, "vopt") for hw in (0.05, 0.1, 0.2, 0.5)] + \
          [("Israeli", 35.0, hw, "vopt") for hw in (0.1, 0.2, 0.5, 1.0)] + [("London", 30.0, 0.05, "FE")] + \
          [("Israeli", 35.0, hw, "FE") for hw in (0.1, 0.2, 0.5, 1.0)]
    for soil, beta, hw, ap in pts:
        thr = paper8(soil, beta, hw, ap) / H8
        specs.append(case8(soil, beta, hw, ap, task="grid", grid=dict(threshold=thr)))
    for a, b in ((1, 35.0), (5, 35.0), (10, 35.0), (1, 45.0)):
        specs.append(case9(a, b, "vopt", task="grid", grid=dict(threshold=paper9(a, b, "vopt"))))
    res = run_all(specs, args, cache, log, "grid")
    out = []
    log(f"   {'case':34s} {'paper':>9s} {'ours(PSO)':>9s} {'grid min':>9s} {'polished':>9s} {'#adm':>7s} "
        f"{'#<paper':>8s} {'frac<paper':>10s}")
    b8 = baseline8()
    b9 = baseline9()
    for s, r in zip(specs, res):
        if r is None:
            continue
        sc = s["H"] if s["fig"] == 8 else 1.0
        if s["fig"] == 8:
            name = f"Fig8 {s['soil']} {s['beta']:.0f} {s['approach']} hw/H={s['hw']}"
            ours = float(b8[(s["soil"], s["beta"], s["approach"], s["hw"])]["ours"])
        else:
            name = f"Fig9 a={s['alpha']:g} b={s['beta']:.0f} vopt"
            ours = float(b9[(int(s["alpha"]), s["beta"], "vopt")]["ours"])
        thr = s["grid"]["threshold"] * sc
        row = dict(case=name, paper=thr, ours_pso=ours, grid_min=r["min"] * sc, grid_min_polished=r.get("min_polished",
                   np.nan) * sc, n_admissible=r["n_adm"], n_below_paper=r["n_below"],
                   frac_below_paper=r["n_below"] / max(r["n_adm"], 1))
        out.append(row)
        log(f"   {name:34s} {thr:9.3f} {ours:9.3f} {row['grid_min']:9.3f} {row['grid_min_polished']:9.3f} "
            f"{r['n_adm']:7d} {r['n_below']:8d} {row['frac_below_paper']:10.4f}")
    write_csv("grid_search.csv", out)
    return out


def exp_restrict(args, cache, log):
    """minimum of Gamma over restricted mechanism sets (grid): which restriction would raise Israeli 35 FE to the
    paper while keeping the matching FE curves (London 30, Israeli 60) unchanged?"""
    subsets = {"all": {}, "class I": {"class": "I"}, "toe only": {"toe_only": True},
               "L>=0.25H": {"L_min": 0.25}, "L>=0.5H": {"L_min": 0.5}, "L>=1H": {"L_min": 1.0},
               "II d>=0.25H": {"class": "II", "d_min": 0.25}, "II d>=0.5H": {"class": "II", "d_min": 0.5},
               "theta2<=pi/2": {"t2_max": np.pi / 2}, "theta2<=1.9": {"t2_max": 1.9},
               "r_h<=1.5H": {"rh_max": 1.5}, "r_h<=1H": {"rh_max": 1.0}}
    pts = [("Israeli", 35.0, 0.5), ("Israeli", 35.0, 1.0), ("London", 30.0, 0.5), ("London", 30.0, 1.0),
           ("Israeli", 60.0, 0.5), ("London", 60.0, 0.5)]
    specs = [case8(s, b, hw, "FE", task="grid", grid=dict(threshold=paper8(s, b, hw, "FE") / H8, subsets=subsets,
                                                          n1=90, n2=90, n3=30, eta_min=0.05))
             for s, b, hw in pts]
    res = run_all(specs, args, cache, log, "restrict")
    out = []
    log("   min Gamma * H [m] over restricted sets (grid, then rel. to the paper value):")
    log(f"   {'case':28s} {'paper':>8s} " + " ".join(f"{k:>12s}" for k in subsets))
    for (s, b, hw), sp, r in zip(pts, specs, res):
        if r is None:
            continue
        pv = paper8(s, b, hw, "FE")
        row = dict(case=f"{s} {b:.0f} FE hw/H={hw}", paper=pv)
        cells = []
        for k in subsets:
            v = r["subsets"][k] * H8
            row[k] = v
            cells.append(f"{v:6.2f}({_rel(v, pv):+.2f})" if np.isfinite(v) else f"{'-':>12s}")
        out.append(row)
        log(f"   {row['case']:28s} {pv:8.3f} " + " ".join(f"{c:>12s}" for c in cells))
    write_csv("restricted_classes_fe.csv", out)


def exp_eq56(args, cache, log):
    la = dict(quad="eq56", classes=("I",), seeds=(0,), n_particles=40, n_iter=120)
    pts = [case8("London", 30.0, 0.1, "vopt"), case8("Israeli", 35.0, 0.2, "vopt"), case8("Israeli", 35.0, 0.5, "vopt"),
           case8("Israeli", 35.0, 0.5, "FE"), case8("London", 60.0, 0.5, "vopt"), case9(1, 35.0, "vopt"),
           case9(5, 35.0, "vopt")]
    specs = [dict(p, la=la) for p in pts]
    res = run_all(specs, args, cache, log, "eq56")
    ref = run_all([dict(p, la=dict(classes=("I",), seeds=(0,))) for p in pts], args, cache, log, "eq56-ref")
    out = []
    for p, r, r0 in zip(pts, res, ref):
        if r is None or r0 is None:
            continue
        pv = paper8(p["soil"], p["beta"], p["hw"], p["approach"]) / H8 if p["fig"] == 8 else paper9(p["alpha"], p["beta"], "vopt")
        name = (f"Fig8 {p['soil']} {p['beta']:.0f} {p['approach']} hw/H={p['hw']}" if p["fig"] == 8
                else f"Fig9 a={p['alpha']:g} b={p['beta']:.0f}")
        out.append(dict(case=name, Gamma_strict=r0["Gamma"], Gamma_eq56_only=r["Gamma"], paper=pv,
                        rel_eq56_vs_strict=_rel(r["Gamma"], r0["Gamma"]), x_eq56=json.dumps(np.round(r["x"], 5).tolist())))
        log(f"   {name:34s} strict class I {r0['Gamma']:.5f} | Eq.56-only (signed) {r['Gamma']:.5f} "
            f"({_rel(r['Gamma'], r0['Gamma']):+.2e}) | paper {pv:.5f}")
    write_csv("eq56_only_class.csv", out)


def exp_fe_settings(args, cache, log):
    variants = {
        "baseline": {},
        "alpha=2": dict(alpha_field=2.0), "alpha=3": dict(alpha_field=3.0), "alpha=5": dict(alpha_field=5.0),
        "impermeable box": dict(fe=dict(bc="impermeable")), "toe_r": dict(fe=dict(bc="toe_r")),
        "box x2 (100/20/60 H)": dict(fe=dict(left=100.0, right=20.0, depth=60.0)),
        "box 50/10/30 m at H_ref=10 m": dict(fe=dict(left=5.0, right=1.0, depth=3.0)),
        "P1 coarse h0=0.25H hs=0.5H": dict(fe=dict(order=1, size=dict(h0=0.25, hs=0.5, grade=0.5, hmax=8.0))),
        "P1 very coarse h0=0.5H hs=H": dict(fe=dict(order=1, size=dict(h0=0.5, hs=1.0, grade=0.5, hmax=8.0))),
        "gamma_w field 9.81": dict(gamma_w_field=9.81), "gamma_w field 10": dict(gamma_w_field=10.0),
        "FE field of beta=30": dict(beta_field=30.0), "FE field of beta=32.5": dict(beta_field=32.5),
        "slope beta=30 (Israeli soil)": dict(beta=30.0), "slope beta=32.5": dict(beta=32.5),
    }
    controls = [("London", 30.0, 0.5), ("London", 60.0, 0.5), ("Israeli", 60.0, 0.5)]
    hws = (0.2, 0.5, 1.0)
    specs, idx = [], []
    for name, kw in variants.items():
        la = dict(seeds=(0,))
        for hw in hws:
            specs.append(case8("Israeli", 35.0, hw, "FE", la=la, **kw))
            idx.append((name, "Israeli", 35.0, hw))
        if name.startswith(("box", "P1", "impermeable", "gamma_w")):
            for s, b, hw in controls:
                kw2 = {k: v for k, v in kw.items() if k not in ("beta",)}
                specs.append(case8(s, b, hw, "FE", la=la, **kw2))
                idx.append((name, s, b, hw))
    res = run_all(specs, args, cache, log, "fe_settings")
    out = []
    for (name, s, b, hw), r in zip(idx, res):
        if r is None:
            continue
        pv = paper8(s, b, hw, "FE")
        out.append(dict(variant=name, soil=s, beta=b, hw_over_H=hw, Hcrit=r["Hcrit"], paper=pv, rel=_rel(r["Hcrit"], pv)))
    log(f"   {'variant':30s} | Israeli 35 FE rel to paper at h_w/H = 0.2, 0.5, 1.0 | controls (FE, h_w/H = 0.5): "
        "London 30, London 60, Israeli 60")
    for name in variants:
        a = [o for o in out if o["variant"] == name and o["soil"] == "Israeli" and o["hw_over_H"] in hws
             and abs(o["beta"] - 35.0) < 1e-9]
        c = [o for o in out if o["variant"] == name and not (o["soil"] == "Israeli" and abs(o["beta"] - 35.0) < 1e-9)]
        log(f"   {name:30s} | " + " ".join(f"{o['rel']:+.3f}" for o in a) + " | " + " ".join(f"{o['rel']:+.3f}" for o in c))
    write_csv("fe_settings_israeli35.csv", out)


def exp_soil(args, cache, log):
    phis = (30.0, 31.0, 31.5, 32.0, 32.5)
    hws = (0.2, 0.5, 1.0)
    la = dict(seeds=(0,))
    specs = []
    for phi in phis:
        specs.append(case8("Israeli", 35.0, 0.0, "none", c=1.0, phi=phi, la=la))
        for hw in hws:
            for ap in ("vopt", "FE"):
                specs.append(case8("Israeli", 35.0, hw, ap, c=1.0, phi=phi, la=la))
    res = {key(s): r for s, r in zip(specs, run_all(specs, args, cache, log, "soil"))}
    dry = paper8("Israeli", 35.0, 0.0, "FE")
    out = []
    log(f"   c chosen so that H_crit(h_w = 0) = {dry:.2f} m (the paper's clipped end); rel to the paper:")
    for phi in phis:
        r0 = res.get(key(case8("Israeli", 35.0, 0.0, "none", c=1.0, phi=phi, la=la)))
        if r0 is None:
            continue
        c = dry / r0["Hcrit"]
        cells = []
        for ap in ("vopt", "FE"):
            for hw in hws:
                r = res.get(key(case8("Israeli", 35.0, hw, ap, c=1.0, phi=phi, la=la)))
                v = c * r["Hcrit"] if r else np.nan
                pv = paper8("Israeli", 35.0, hw, ap)
                out.append(dict(phi=phi, c=c, curve=ap, hw_over_H=hw, Hcrit=v, paper=pv, rel=_rel(v, pv)))
                cells.append(f"{ap} {hw}: {_rel(v, pv):+.3f}")
        log(f"   phi = {phi:4.1f}, c = {c:6.3f}: " + " | ".join(cells))
    write_csv("soil_family_israeli35.csv", out)


def exp_loading(args, cache, log):
    pts = [("Israeli", 35.0, 0.2), ("Israeli", 35.0, 0.5), ("Israeli", 35.0, 1.0), ("London", 30.0, 0.5),
           ("London", 30.0, 1.0), ("London", 60.0, 0.5)]
    specs, idx = [], []
    for s, b, hw in pts:
        specs.append(case8(s, b, hw, "FE", gamma_w_weight=0.0, la=dict(seeds=(0,))))
        idx.append(("gamma' -> gamma = 18 (seepage field unchanged)", s, b, hw))
        specs.append(case8(s, b, hw, "FE", loading="sat_above_water", la=dict(seeds=(0,), pu_method="domain")))
        idx.append(("saturated gamma above the lowered water level", s, b, hw))
    res = run_all(specs, args, cache, log, "loading")
    out = []
    for (name, s, b, hw), r in zip(idx, res):
        if r is None:
            continue
        pv = paper8(s, b, hw, "FE")
        out.append(dict(variant=name, soil=s, beta=b, hw_over_H=hw, Hcrit=r["Hcrit"], paper=pv, rel=_rel(r["Hcrit"], pv)))
        log(f"   {name:48s} {s:7s} {b:3.0f} FE hw/H={hw}: {r['Hcrit']:8.3f} m vs paper {pv:8.3f} ({_rel(r['Hcrit'], pv):+.3f})")
    write_csv("loading_variants.csv", out)


def exp_fe_scale(args, cache, log):
    ks = (0.75, 0.8, 0.85)
    hws = (0.1, 0.2, 0.3, 0.5, 0.7, 1.0)
    specs = [case8("Israeli", 35.0, hw, "FE", fe_scale=k, la=dict(seeds=(0,))) for k in ks for hw in hws]
    res = run_all(specs, args, cache, log, "fe_scale")
    out = []
    for s, r in zip(specs, res):
        if r is None:
            continue
        pv = paper8("Israeli", 35.0, s["hw"], "FE")
        out.append(dict(k=s["fe_scale"], hw_over_H=s["hw"], Hcrit=r["Hcrit"], paper=pv, rel=_rel(r["Hcrit"], pv)))
    for k in ks:
        rr = [o["rel"] for o in out if o["k"] == k]
        log(f"   FE seepage force x {k}: rel to paper at h_w/H = {hws}: " + " ".join(f"{v:+.3f}" for v in rr)
            + (f"  (h_w >= 0.2: mean {np.mean(rr[1:]):+.3f}, spread {np.ptp(rr[1:]):.3f})" if len(rr) == len(hws) else ""))
    write_csv("fe_scale_israeli35.csv", out)


def exp_vopt_field(args, cache, log):
    la = dict(seeds=(0,))
    pts = [("Fig9 a=1 b=35", case9(1, 35.0, "vopt")), ("Fig9 a=1 b=45", case9(1, 45.0, "vopt")),
           ("Fig9 a=1 b=60", case9(1, 60.0, "vopt")), ("Fig9 a=5 b=35", case9(5, 35.0, "vopt")),
           ("Isr35 hw=0.2", case8("Israeli", 35.0, 0.2, "vopt")), ("Isr35 hw=0.5", case8("Israeli", 35.0, 0.5, "vopt")),
           ("Isr35 hw=1", case8("Israeli", 35.0, 1.0, "vopt")), ("Lon30 hw=0.1", case8("London", 30.0, 0.1, "vopt")),
           ("Lon30 hw=0.5", case8("London", 30.0, 0.5, "vopt")), ("Lon60 hw=0.5", case8("London", 60.0, 0.5, "vopt"))]
    variants = {"optimal (baseline)": {}, "dm=-0.2": dict(vopt=dict(variant="dm", dm=-0.2)),
                "dm=-0.1": dict(vopt=dict(variant="dm", dm=-0.1)), "dm=+0.1": dict(vopt=dict(variant="dm", dm=0.1)),
                "dm=+0.2": dict(vopt=dict(variant="dm", dm=0.2)), "one step from m=1": dict(vopt=dict(variant="dm", m1=1.0)),
                "degenerate m->0": dict(vopt=dict(variant="degenerate")),
                "L_m=2H": dict(vopt=dict(Lm_over_H=2.0)), "L_m=0.5H": dict(vopt=dict(Lm_over_H=0.5)),
                "zone 3 removed": dict(vopt=dict(variant="nozone3")), "gamma_w field 10": dict(gamma_w_field=10.0)}
    specs, idx = [], []
    for vn, kw in variants.items():
        for pn, p in pts:
            specs.append(dict(p, la=la, **kw))
            idx.append((vn, pn, p))
    res = run_all(specs, args, cache, log, "vopt_field")
    # J* of the variants at h_w = H (Fig. 5 check)
    jspecs = [dict(case9(1, b, "vopt"), task="field", rev=2, **kw) for kw in variants.values() for b in (30.0, 45.0, 60.0)]
    jres = run_all(jspecs, args, cache, log, "vopt_field J*")
    out = []
    base = {pn: r for (vn, pn, p), r in zip(idx, res) if vn == "optimal (baseline)" and r is not None}
    for (vn, pn, p), r in zip(idx, res):
        if r is None:
            continue
        pv = paper8(p["soil"], p["beta"], p["hw"], "vopt") / H8 if p["fig"] == 8 else paper9(p["alpha"], p["beta"], "vopt")
        out.append(dict(variant=vn, case=pn, Gamma=r["Gamma"], m=r.get("m"), paper=pv, rel_paper=_rel(r["Gamma"], pv),
                        rel_baseline=_rel(r["Gamma"], base[pn]["Gamma"]) if pn in base else np.nan))
    log("   rel. to the paper (rel. to our optimal field) per case:")
    log(f"   {'variant':20s} " + " ".join(f"{pn:>17s}" for pn, _ in pts))
    for vn in variants:
        cells = [o for o in out if o["variant"] == vn]
        log(f"   {vn:20s} " + " ".join(f"{o['rel_paper']:+.3f}({o['rel_baseline']:+.3f})" for o in cells))
    jb = {}
    for s, r in zip(jspecs, jres):
        if r is not None:
            jb.setdefault(s["beta"], {})[json.dumps({k: v for k, v in s.items() if k in ("vopt", "gamma_w_field")})] = r["Jn"]
    log("   -J*/(k_h H^2 gw^2) at h_w = H, alpha = 1 (Fig. 5 solid curve) of each variant, rel. to the optimal field:")
    for (vn, kw) in variants.items():
        kk = json.dumps({k: v for k, v in kw.items() if k in ("vopt", "gamma_w_field")})
        cells = []
        for b in (30.0, 45.0, 60.0):
            j0 = jb.get(b, {}).get(json.dumps({}))
            j1 = jb.get(b, {}).get(kk)
            cells.append(f"beta {b:.0f}: {_rel(j1, j0) if (j0 and j1) else np.nan:+.2e}")
        log(f"   {vn:20s} " + "  ".join(cells))
    write_csv("vopt_field_variants.csv", out)


def exp_gamma_w(args, cache, log):
    la = dict(seeds=(0,))
    specs, idx = [], []
    for a in (1, 5, 10):
        for b in (30.0, 35.0, 40.0, 45.0, 50.0, 60.0, 75.0, 90.0):
            for ap in ("vopt", "FE"):
                for name, kw in (("9.81 (baseline)", {}), ("9.8", dict(gamma_w=9.8)), ("10", dict(gamma_w=10.0)),
                                 ("10 in the field only", dict(gamma_w_field=10.0))):
                    specs.append(case9(a, b, ap, la=la, **kw))
                    idx.append((name, a, b, ap))
    for s, b in PANELS8:
        for hw in (0.1, 0.5, 1.0):
            for ap in ("vopt", "FE"):
                for name, kw in (("9.8 (baseline)", {}), ("9.81 in the field", dict(gamma_w_field=9.81)),
                                 ("10 in the field", dict(gamma_w_field=10.0))):
                    specs.append(case8(s, b, hw, ap, la=la, **kw))
                    idx.append((name, s, b, hw, ap))
    res = run_all(specs, args, cache, log, "gamma_w")
    out9, out8 = [], []
    for k, r in zip(idx, res):
        if r is None:
            continue
        if len(k) == 4:
            name, a, b, ap = k
            pv, st = interp_polyline(P9[(int(a), ap)], b)
            out9.append(dict(gamma_w=name, alpha=a, beta=b, curve=ap, Gamma=r["Gamma"], paper=pv, status=st,
                             rel=_rel(r["Gamma"], pv) if st == "vertex" else np.nan))
        else:
            name, s, b, hw, ap = k
            pv = paper8(s, b, hw, ap)
            out8.append(dict(gamma_w=name, soil=s, beta=b, hw_over_H=hw, curve=ap, Hcrit=r["Hcrit"], paper=pv,
                             rel=_rel(r["Hcrit"], pv)))
    log("   Fig. 9 (nodes with a visible paper value): mean / rms of ours/paper - 1")
    for name in ("9.81 (baseline)", "9.8", "10", "10 in the field only"):
        for ap in ("vopt", "FE"):
            e = np.array([o["rel"] for o in out9 if o["gamma_w"] == name and o["curve"] == ap and np.isfinite(o["rel"])])
            e3 = np.array([o["rel"] for o in out9 if o["gamma_w"] == name and o["curve"] == ap and np.isfinite(o["rel"])
                           and 30.0 <= o["beta"] <= 50.0])
            if len(e):
                log(f"   gamma_w {name:22s} {ap:4s}: all beta n={len(e):2d} mean {e.mean():+.4f} rms {np.sqrt(np.mean(e ** 2)):.4f}"
                    f" | beta 30-50: mean {e3.mean():+.4f} rms {np.sqrt(np.mean(e3 ** 2)):.4f}")
    log("   Fig. 8 (h_w/H = 0.1, 0.5, 1): mean / rms of ours/paper - 1, panels London 60 / Israeli 60 / London 30 + Israeli 35")
    for name in ("9.8 (baseline)", "9.81 in the field", "10 in the field"):
        for grp, sel in (("beta 60 panels", lambda o: o["beta"] == 60.0), ("London 30", lambda o: o["soil"] == "London" and o["beta"] == 30.0),
                         ("Israeli 35", lambda o: o["soil"] == "Israeli" and o["beta"] == 35.0)):
            e = np.array([o["rel"] for o in out8 if o["gamma_w"] == name and sel(o)])
            if len(e):
                log(f"   gamma_w {name:20s} {grp:15s}: mean {e.mean():+.4f} rms {np.sqrt(np.mean(e ** 2)):.4f}")
    write_csv("gamma_w_fig9.csv", out9)
    write_csv("gamma_w_fig8.csv", out8)


def quadrature_specs():
    la0 = dict(classes=("I",), seeds=(0, 1), n_particles=30, n_iter=100)
    pts = []
    for s, b in PANELS8:
        for hw in (0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0):
            for ap in ("vopt", "FE"):
                pts.append(case8(s, b, hw, ap))
    for a in (1, 5, 10):
        for b in (30.0, 35.0, 40.0, 45.0, 50.0, 60.0, 75.0, 90.0):
            for ap in ("vopt", "FE"):
                pts.append(case9(a, b, ap))
    specs = []
    for p in pts:
        for n in (25, 50):
            specs.append(dict(p, la=dict(la0, quad="cellgrid", n_cell=n)))
    return pts, specs


def exp_quadrature(args, cache, log):
    pts, specs = quadrature_specs()
    # cheap first: Fig. 9 and the vopt curves, the FE curves of Fig. 8 last
    order = sorted(range(len(specs)), key=lambda i: (specs[i]["approach"] == "FE", specs[i]["fig"] == 8,
                                                     specs[i]["la"]["n_cell"]))
    res_o = run_all([specs[i] for i in order], args, cache, log, "quadrature")
    res = [None] * len(specs)
    for j, i in enumerate(order):
        res[i] = res_o[j]
    b8 = baseline8()
    b9 = baseline9()
    out = []
    for s, r in zip(specs, res):
        if r is None:
            continue
        if s["fig"] == 8:
            bl = b8.get((s["soil"], s["beta"], s["approach"], s["hw"]))
            ours = float(bl["ours"]) / H8
            pv = paper8(s["soil"], s["beta"], s["hw"], s["approach"]) / H8
            name = f"Fig8 {s['soil']} {s['beta']:.0f}"
            st = bl["status"]
        else:
            bl = b9.get((int(s["alpha"]), s["beta"], s["approach"]))
            ours = float(bl["ours"])
            pv = paper9(s["alpha"], s["beta"], s["approach"])
            name = f"Fig9 a={s['alpha']:g}"
            st = "vertex" if np.isfinite(pv) else "none"
        out.append(dict(case=name, beta=s["beta"], hw_over_H=s["hw"], curve=s["approach"], n_cell=s["la"]["n_cell"],
                        ours_exact=ours, Gamma_emulated=r["Gamma"], Gamma_exact_at_emulated=r["Gamma_exact"],
                        bias=_rel(r["Gamma"], ours), paper=pv, observed=_rel(pv, ours), status=st))
    write_csv("quadrature_emulation.csv", out)
    log("   reported minimum of the emulated objective rel. to our exact minimum ('bias') vs the observed paper/ours - 1")
    log("   (h = H/25 and H/50; PSO 2 seeds x 30 particles x 100 iterations + Nelder-Mead on the same objective)")
    groups = [("Fig8 London 30", "vopt"), ("Fig8 Israeli 35", "vopt"), ("Fig8 London 60", "vopt"), ("Fig8 Israeli 60", "vopt"),
              ("Fig8 London 30", "FE"), ("Fig8 Israeli 35", "FE"), ("Fig8 London 60", "FE"), ("Fig8 Israeli 60", "FE")]
    for g, ap in groups:
        rows = sorted([o for o in out if o["case"] == g and o["curve"] == ap], key=lambda o: (o["hw_over_H"], o["n_cell"]))
        if not rows:
            continue
        cells = {}
        for o in rows:
            cells.setdefault(o["hw_over_H"], {})[o["n_cell"]] = o
        txt = []
        for hw, d in sorted(cells.items()):
            ob = next(iter(d.values()))["observed"]
            txt.append(f"{hw:g}: obs {ob:+.3f} H/25 {d.get(25, {}).get('bias', np.nan):+.3f} "
                       f"H/50 {d.get(50, {}).get('bias', np.nan):+.3f}")
        log(f"   {g:16s} {ap:4s} " + " | ".join(txt))
    for a in (1, 5, 10):
        for ap in ("vopt", "FE"):
            rows = sorted([o for o in out if o["case"] == f"Fig9 a={a}" and o["curve"] == ap], key=lambda o: (o["beta"], o["n_cell"]))
            cells = {}
            for o in rows:
                cells.setdefault(o["beta"], {})[o["n_cell"]] = o
            txt = []
            for b, d in sorted(cells.items()):
                ob = next(iter(d.values()))["observed"]
                txt.append(f"{b:g}: obs {ob:+.3f} /25 {d.get(25, {}).get('bias', np.nan):+.3f} "
                           f"/50 {d.get(50, {}).get('bias', np.nan):+.3f}")
            if txt:
                log(f"   Fig9 a={a:<2d} {ap:4s} " + " | ".join(txt))
    # correlation between emulated bias and observed difference
    for n in (25, 50):
        e = np.array([(o["bias"], o["observed"]) for o in out if o["n_cell"] == n and np.isfinite(o["observed"])
                      and np.isfinite(o["bias"]) and o["status"] in ("vertex", "segment")
                      and not (o["case"] == "Fig8 Israeli 35" and o["curve"] == "FE")])
        if len(e) > 3:
            rms0 = np.sqrt(np.mean(e[:, 1] ** 2))
            rms1 = np.sqrt(np.mean((e[:, 1] - e[:, 0]) ** 2))
            log(f"   h = H/{n}: {len(e)} points (Israeli 35 FE excluded): corr(bias, observed) = "
                f"{np.corrcoef(e[:, 0], e[:, 1])[0, 1]:+.3f}; rms of observed {rms0:.4f} -> rms of observed - bias {rms1:.4f}")
    return out


def exp_digitization(args, cache, log):
    b9 = baseline9()
    out = []
    for a in (1, 5, 10):
        for ap in ("vopt", "FE"):
            rows = []
            for (aa, b, cc), r in sorted(b9.items()):
                if aa != a or cc != ap:
                    continue
                pv = paper9(a, b, ap)
                st = interp_polyline(P9[(a, ap)], b)[1]
                if st != "vertex":
                    continue
                ours = float(r["ours"])
                rows.append((b, ours, pv))
            if len(rows) < 3:
                continue
            arr = np.array(rows)
            d = arr[:, 1] - arr[:, 2]
            # constant offset model (stroke-centre bias in Gamma units) vs proportional model
            off = d.mean()
            res_off = np.sqrt(np.mean((d - off) ** 2))
            k = np.sum(d * arr[:, 2]) / np.sum(arr[:, 2] ** 2)
            res_prop = np.sqrt(np.mean((d - k * arr[:, 2]) ** 2))
            out.append(dict(alpha=a, curve=ap, n=len(rows), mean_abs_diff=off, rms_about_constant_offset=res_off,
                            proportional_factor=k, rms_about_proportional=res_prop,
                            one_pixel_in_Gamma=1.0 / 87.608))
            log(f"   Fig. 9 alpha={a:2d} {ap:4s}: ours - paper (Gamma units) " +
                " ".join(f"{b:g}:{x:+.3f}" for (b, _, _), x in zip(rows, d)))
            log(f"        constant offset {off:+.4f} leaves rms {res_off:.4f}; proportional {100 * k:+.2f} % leaves rms "
                f"{res_prop:.4f}  (1 px of the 300 ppi bitmap = {1 / 87.608:.4f})")
    with open(os.path.join(DATA_DIR, "digitize_meta.json")) as fh:
        meta = json.load(fh)
    v = meta["fig8"]["validation"]["vector_polyline_fit"]
    log("   Fig. 8 (vector strokes, centre lines of the PDF paths): max |log10 misfit| of the vertex polylines: " +
        ", ".join(f"{k}: {val['max_log10']:.1e}" for k, val in v.items() if "Wu" not in k))
    write_csv("digitization_fig9_offsets.csv", out)


def _read(name):
    path = os.path.join(OUT_DIR, name)
    if not os.path.exists(path):
        return []
    with open(path) as fh:
        rows = list(csv.DictReader(fh))
    for r in rows:
        for k, v in r.items():
            try:
                r[k] = float(v)
            except (TypeError, ValueError):
                pass
    return rows


def exp_summary(args, cache, log):
    """ranked list of the tested explanations with their numbers (from the CSV files of the other experiments)"""
    L = log
    g = {r["case"]: r for r in _read("grid_search.csv")}
    zs = _read("zones_pu_split.csv")
    qe = _read("quadrature_emulation.csv")
    vf = _read("vopt_field_variants.csv")
    gw9 = _read("gamma_w_fig9.csv")
    gw8 = _read("gamma_w_fig8.csv")
    fs = _read("fe_settings_israeli35.csv")
    so = _read("soil_family_israeli35.csv")
    sc = _read("fe_scale_israeli35.csv")
    e56 = _read("eq56_only_class.csv")
    dg = _read("digitization_fig9_offsets.csv")
    ss = _read("selfsimilar.csv")
    rc = {r["case"]: r for r in _read("restricted_classes_fe.csv")}
    lo = _read("loading_variants.csv")

    def gl(case, k):
        return g.get(case, {}).get(k, np.nan)

    def qrow(case, curve, x, n):
        key_ = "beta" if case.startswith("Fig9") else "hw_over_H"
        r = [q for q in qe if q["case"] == case and q["curve"] == curve and abs(q[key_] - x) < 1e-9 and q["n_cell"] == n]
        return r[0] if r else None

    def vrow(v, c, k="rel_baseline"):
        r = [q for q in vf if q["variant"] == v and q["case"] == c]
        return r[0][k] if r else np.nan

    L("Ranked list (most to least supported).  'obs' = paper / ours - 1 at the paper's plotted vertices.")
    L("")
    L("1. SUPPORTED (explains (2), (3) and the +0.7 % of the Fig. 9 FE curve): the paper's minima are minima of a")
    L("   NOISY objective.  Its P_u was evaluated with a fixed (absolute) resolution, and the PSO converges to")
    L("   mechanisms where the integration error favours failure, so every reported Gamma is biased LOW; the bias")
    L("   is large for small or thin mechanisms next to O (local K^-1 v'_opt mechanisms at small h_w, optima that")
    L("   start 0.004-0.03 H behind O) and negligible for the large mechanisms of the beta = 60 panels.")
    L("   a) Our objective is exact.  P_u at the optimal mechanisms converges to <= 8e-6 ('fine' vs 'ref',")
    L("      quadrature_convergence.csv).  Brute-force grids with log-spaced eta >= 0.01 contain no mechanism below")
    L("      the paper value (grid_search.csv): " + "; ".join(
        f"{c.replace('Fig8 ', '').replace(' vopt', '').replace('hw/H=', '')} {gl(c, 'grid_min'):.4g} vs paper {gl(c, 'paper'):.4g}"
        for c in ("Fig8 Israeli 35 vopt hw/H=0.1", "Fig8 Israeli 35 vopt hw/H=0.2", "Fig8 London 30 vopt hw/H=0.1",
                  "Fig9 a=5 b=35 vopt") if c in g) + ".")
    L("      The Eq. 56-only class (no 'spiral below O' check, signed integrals) changes Gamma by <= "
      f"{max(abs(r['rel_eq56_vs_strict']) for r in e56):.1e} (eq56_only_class.csv)." if e56 else "")
    if ss:
        for soil, beta in (("Israeli", 35.0), ("London", 30.0)):
            rr = [r for r in ss if r["soil"] == soil and r["beta"] == beta and r["curve"] == "vopt"]
            loc = [r["hw_over_H"] for r in rr if r["s"] < 0.999]
            L(f"   b) {soil} {beta:.0f} vopt: for h_w/H = {min(loc):g}-{max(loc):g} our optimum is ONE mechanism scaled with "
              f"R_w (exit B on the face at {np.median([r['B_over_Rw'] for r in rr if r['s'] < 0.999]):.3f} R_w from O): H_crit h_w/H = "
              + "/".join(f"{r['ours_x_hw']:.2f}" for r in rr if 0.1 <= r["hw_over_H"] <= 0.7)
              + " m (h_w/H = 0.1..0.7); the paper's product falls as h_w decreases: "
              + "/".join(f"{r['paper_x_hw']:.2f}" for r in rr if 0.1 <= r["hw_over_H"] <= 0.7)
              + " m.  An exact objective is scale invariant; the paper's is not (absolute resolution).")
    if qe:
        L("   c) Emulation: P_u by a cell-midpoint grid of FIXED spacing h (absolute coordinates) minimised by the same")
        L("      PSO (2 seeds x 30 particles x 100 iterations, class I).  Reported-minimum bias vs observed (obs):")
        for case, curve, xs in (("Fig8 Israeli 35", "vopt", (0.1, 0.2, 0.3, 0.5, 0.7, 1.0)),
                                ("Fig8 London 30", "vopt", (0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0)),
                                ("Fig8 London 60", "vopt", (0.1, 0.5, 1.0)), ("Fig8 Israeli 60", "vopt", (0.1, 0.5, 1.0)),
                                ("Fig8 London 30", "FE", (0.05, 0.1, 0.5, 1.0)),
                                ("Fig9 a=1", "vopt", (30.0, 35.0, 40.0, 50.0, 60.0, 90.0)),
                                ("Fig9 a=5", "vopt", (35.0, 40.0, 50.0, 60.0, 90.0)),
                                ("Fig9 a=10", "vopt", (35.0, 40.0, 50.0, 60.0, 90.0)),
                                ("Fig9 a=1", "FE", (30.0, 35.0, 40.0, 50.0, 60.0, 90.0))):
            cells = []
            for x in xs:
                a, b = qrow(case, curve, x, 25), qrow(case, curve, x, 50)
                if a is None or b is None:
                    continue
                cells.append(f"{x:g}: obs {100 * a['observed']:+.1f} | {100 * a['bias']:+.1f} / {100 * b['bias']:+.1f}")
            L(f"      {case:15s} {curve:4s} [%] obs | h=H/25 / H/50:  " + ";  ".join(cells))
        for n in (25, 50):
            e = np.array([(r["bias"], r["observed"]) for r in qe if r["n_cell"] == n and np.isfinite(r["observed"])
                          and np.isfinite(r["bias"]) and r["status"] in ("vertex", "segment")
                          and not (r["case"] == "Fig8 Israeli 35" and r["curve"] == "FE")])
            if len(e) > 3:
                L(f"      h = H/{n}: {len(e)} points (Israeli 35 FE excluded): corr(bias, obs) = "
                  f"{np.corrcoef(e[:, 0], e[:, 1])[0, 1]:+.2f}; rms(obs) {100 * np.sqrt(np.mean(e[:, 1] ** 2)):.2f} % -> "
                  f"rms(obs - bias) {100 * np.sqrt(np.mean((e[:, 1] - e[:, 0]) ** 2)):.2f} %")
        L("      The paper's resolution lies between H/25 and H/50 (about H/35).  beta = 45 in Fig. 9 is an artefact of")
        L("      the emulation (face along the cell diagonals: -5..-13 %), which itself shows how irregular such errors are.")
    if zs:
        e9 = np.array([(r["frac_core005"] + r["frac_core_sstar"], r["paper_over_ours"]) for r in zs
                       if r["case"].startswith("Fig9") and np.isfinite(r["paper_over_ours"])])
        if len(e9) > 3:
            L(f"   d) Fig. 9: obs is more negative where a larger share of P_u comes from the singular fan r < s* R_w around O"
              f" (corr {np.corrcoef(e9[:, 0], e9[:, 1])[0, 1]:+.2f}, {len(e9)} optima; share 25-55 % for beta <= 50, "
              "~0 for beta = 75) -- zones_pu_split.csv, mechanisms_zoom_O.png, mechanisms_fig9_vopt.png.")
    L("")
    L("2. POSSIBLE COMPLEMENT for (3), not decisive: a zone-1 exponent m ~10 % above the optimum (e.g. a fixed point")
    L("   m = sqrt(C/D) not converged).  J* is stationary in m, so Fig. 5 cannot see it (-J* changes by "
      "1.5e-4 / 0.8e-4 at beta 30 / 60 for dm = +0.1), but the A-at-O optima are sensitive:")
    if vf:
        cases = list(dict.fromkeys(r["case"] for r in vf))
        L("   ours/paper - 1 [%], optimal m -> m (1 + 0.1): " + "; ".join(
            f"{c}: {100 * vrow('optimal (baseline)', c, 'rel_paper'):+.1f} -> {100 * vrow('dm=+0.1', c, 'rel_paper'):+.1f}"
            for c in cases))
        L("   It removes the Fig. 9 residual at beta 35-60 for alpha = 1 but not the small-h_w residuals (Lon30 0.1,")
        L("   Isr35 0.2) and worsens London 60.  'One step from m = 1' changes -J* by -2.8e-3 at beta 60 (Fig. 5 match is")
        L("   2e-4): excluded.")
    L("")
    if gw9:
        def gm(name, curve, lo_=0.0, hi_=90.0):
            e = np.array([r["rel"] for r in gw9 if str(r["gamma_w"]).rstrip("0").rstrip(".") == name and r["curve"] == curve
                          and np.isfinite(r["rel"]) and lo_ <= r["beta"] <= hi_])
            return (100 * e.mean(), 100 * np.sqrt(np.mean(e ** 2))) if len(e) else (np.nan, np.nan)
        L("3. UNLIKELY: gamma_w = 10 in Fig. 9.  vopt beta 30-50 mean/rms %.2f/%.2f %% -> %.2f/%.2f %%, but FE (all beta)"
          " %.2f/%.2f %% -> %.2f/%.2f %%, and the paper's own Fig. 4b legend (-9.81 for a unit drawdown) and the Fig. 8"
          " dry ends (9.8) say otherwise; gamma_w = 9.8 changes Fig. 9 by +0.05 %%." % (
              *gm("9.81 (baseline)", "vopt", 30, 50), *gm("10", "vopt", 30, 50), *gm("9.81 (baseline)", "FE"), *gm("10", "FE")))
    L("")
    L("4. REFUTED")
    L("   - optimiser of the paper missed OUR better mechanisms (for the 'ours lower' Israeli 35 FE points): "
      f"{100 * gl('Fig8 Israeli 35 FE hw/H=0.2', 'frac_below_paper'):.1f} / {100 * gl('Fig8 Israeli 35 FE hw/H=0.5', 'frac_below_paper'):.1f}"
      f" / {100 * gl('Fig8 Israeli 35 FE hw/H=1.0', 'frac_below_paper'):.1f} % of the admissible grid beats the paper at h_w/H = 0.2 / 0.5 / 1;"
      " of the whole PSO box 0.8 % (0.5) and 0.1 % (1.0), whereas the paper reached the 1 %-of-optimum regions of the"
      " London 30 / Israeli 60 FE problems, which are 2e-5 of the box (exploratory count, 4e5 random points).")
    if rc:
        r1, r2 = rc.get("Israeli 35 FE hw/H=0.5"), rc.get("Israeli 35 FE hw/H=1.0")
        L("   - restricted mechanism sets (toe only, L >= 0.25/0.5/1 H, II with d >= 0.25/0.5 H, theta2 <= pi/2 or 1.9,"
          " r_h <= 1 or 1.5 H): none raises Israeli 35 FE to the paper at both h_w/H = 0.5 and 1 without moving the"
          " matching London 30 / London 60 / Israeli 60 FE points by +10 to +86 % (restricted_classes_fe.csv).")
    L("   - mechanism class of Eq. 56 only / mechanism II: <= 6e-5 change; class II never lower at these points.")
    if fs:
        def frow(v):
            rr = [r for r in fs if r["variant"] == v and r["soil"] == "Israeli" and r["beta"] == 35.0]
            cc = [r for r in fs if r["variant"] == v and not (r["soil"] == "Israeli" and r["beta"] == 35.0)]
            return " ".join(f"{100 * r['rel']:+.1f}" for r in rr) + (" | controls " + " ".join(f"{100 * r['rel']:+.1f}" for r in cc) if cc else "")
        L("   - FE settings for Israeli 35 FE (ours/paper - 1 [%] at h_w/H = 0.2 0.5 1 | controls London 30, London 60,"
          " Israeli 60 FE at 0.5):")
        for v in ("baseline", "impermeable box", "toe_r", "box x2 (100/20/60 H)", "box 50/10/30 m at H_ref=10 m",
                  "P1 coarse h0=0.25H hs=0.5H", "P1 very coarse h0=0.5H hs=H", "gamma_w field 9.81", "gamma_w field 10",
                  "alpha=2", "alpha=3", "alpha=5", "FE field of beta=30", "FE field of beta=32.5",
                  "slope beta=30 (Israeli soil)", "slope beta=32.5"):
            L(f"       {v:30s} {frow(v)}")
    if so:
        L("   - soil swap / (c, phi) with the dry end 229.3 m kept: " + "; ".join(
            f"phi {phi:g} c {c:.2f}: vopt " + " ".join(f"{100 * r['rel']:+.0f}" for r in so if r["phi"] == phi and r["curve"] == "vopt")
            + " FE " + " ".join(f"{100 * r['rel']:+.0f}" for r in so if r["phi"] == phi and r["curve"] == "FE")
            for phi, c in sorted({(r["phi"], r["c"]) for r in so})) + " (ours/paper - 1 [%] at h_w/H = 0.2 0.5 1): no single"
          " soil fits both curves.")
    if lo:
        L("   - loading (ours/paper - 1): " + "; ".join(f"{r['variant'][:20]} {r['soil']} {r['beta']:.0f} {r['hw_over_H']:g}: {100 * r['rel']:+.0f} %"
                                    for r in lo))
    if vf:
        L(f"   - analytical-field variants: degenerate m -> 0 for every beta (-J* -7.2 % at beta 30: excluded by Fig. 5; Gamma "
          f"{100 * vrow('degenerate m->0', 'Fig9 a=1 b=35'):+.0f} % at Fig. 9 beta 35); L_m = 2 H / 0.5 H and zone 3 removed: Gamma "
          f"unchanged (|change| <= {max(abs(vrow(v, c)) for v in ('L_m=2H', 'L_m=0.5H', 'zone 3 removed') for c in dict.fromkeys(r['case'] for r in vf)):.0e}"
          ", the optimal mechanisms never reach R; J* would change by -43..-75 %); printed Eq. 31: -28..+1500 %"
          " (results/reproduce_python/diag_eq31_variants.csv).")
    if dg:
        L("   - digitization / stroke-centre offsets: Fig. 8 = PDF path centre lines, polyline misfit <= 1.5e-4 in log10 "
          "(0.03 %); Fig. 9 residuals grow with Gamma (a constant offset leaves rms " + ", ".join(
              f"{r['rms_about_constant_offset']:.3f}" for r in dg if r["curve"] == "vopt") + " vs 1 px = 0.011).")
    L("   - FE box / H_ref for Fig. 8: box 50/10/30 m at H_ref = 10 m lowers Gamma (wrong sign for Israeli 35 FE) and moves"
      " the matching controls by -3..-10 %; box x2 or impermeable: <= 0.4 %.")
    L("")
    L("5. NOT EXPLAINED")
    L("   - (1) Israeli 35 FE, ours 19-27 % lower at h_w/H >= 0.2: every tested FE setting gives the same Gamma within 2 %,")
    L("     the paper's optimiser cannot have missed 3-9 % of the admissible space, and numerical noise biases LOW (emulated")
    L("     bias -0.5..-6 %, opposite sign).  The paper's curve corresponds to FE seepage forces scaled by ~0.8 ("
      + (" ".join(f"{100 * r['rel']:+.1f}" for r in sc if r["k"] == 0.8) if sc else "") + " % at h_w/H = 0.1 0.2 0.3 0.5 0.7 1)")
    L("     or to (phi 31.5, c 7.59) for this curve only.  The paper's own FE/vopt ratio for this panel (0.76-0.91)")
    L("     differs from that of every other panel and of Fig. 9 at beta = 35 (0.69), where ours equals the paper's:")
    L("     most likely an inconsistent setting of that single run in the paper.  Our curve should stay as it is.")
    L("   - isolated paper points London 30 vopt h_w/H = 0.1 (obs -9.6 %, emulated -2.9 %) and London 30 FE 0.05 (obs")
    L("     -11.8 %, emulated -1.6 %): kinks of the paper's curves, consistent with lucky outliers of a noisy objective but")
    L("     not reproduced by the emulation.  Israeli 35 h_w/H = 0.05 has NO paper vertex (vertex step 0.1): the 'paper'")
    L("     value there is the plotted chord, not a computed point.")
    plot_summary(zs, qe, vf, log)


def plot_summary(zs, qe, vf, log):
    """observed residual paper/ours - 1 of the K^-1 v'_opt points vs (a) share of P_u from r < s* R_w around O,
    (b) the emulated quadrature bias (h = H/25, H/50), (c) Gamma change of the m-perturbed fields"""
    if not (zs and qe):
        return
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(1, 3, figsize=(12.6, 4.0))
    ax = axs[0]
    for r in zs:
        col = "#2a78d6" if r["case"].startswith("Fig9") else "#eb6834"
        ax.plot(100 * (r["frac_core005"] + r["frac_core_sstar"]), 100 * r["paper_over_ours"], "o", color=col, ms=4)
    ax.plot([], [], "o", color="#2a78d6", label="Fig. 9 vopt optima")
    ax.plot([], [], "o", color="#eb6834", label="Fig. 8 vopt optima")
    ax.set_xlabel("share of P_u from r < s* R_w around O (%)")
    ax.set_ylabel("paper / ours - 1 (%)")
    ax.legend(fontsize=7)
    ax = axs[1]
    for n, mk in ((25, "s"), (50, "o")):
        for curve, col in (("vopt", "#2a78d6"), ("FE", "#eb6834")):
            e = np.array([(r["bias"], r["observed"]) for r in qe if r["n_cell"] == n and r["curve"] == curve
                          and np.isfinite(r["observed"]) and np.isfinite(r["bias"]) and r["status"] in ("vertex", "segment")
                          and not (r["case"] == "Fig8 Israeli 35" and r["curve"] == "FE")])
            if len(e):
                ax.plot(100 * e[:, 0], 100 * e[:, 1], mk, color=col, ms=3.5, mfc="none" if n == 25 else col,
                        label=f"{curve}, h = H/{n}")
    lim = (-65, 5)
    ax.plot(lim, lim, "-", color="#52514e", lw=0.8)
    ax.set_xlim(*lim)
    ax.set_ylim(*lim)
    ax.set_xlabel("emulated bias of the reported minimum (%)")
    ax.set_ylabel("paper / ours - 1 (%)")
    ax.legend(fontsize=7)
    ax = axs[2]
    if vf:
        cases = [c for c in dict.fromkeys(r["case"] for r in vf)]
        for v, col in (("dm=-0.1", "#a9a8a2"), ("dm=+0.1", "#2a78d6"), ("dm=+0.2", "#1baf7a"), ("one step from m=1", "#eb6834")):
            rr = {r["case"]: r for r in vf if r["variant"] == v}
            ax.plot(range(len(cases)), [100 * rr[c]["rel_paper"] for c in cases], "o-", color=col, ms=3.5, lw=1, label=v)
        rr = {r["case"]: r for r in vf if r["variant"] == "optimal (baseline)"}
        ax.plot(range(len(cases)), [100 * rr[c]["rel_paper"] for c in cases], "k^-", ms=4, lw=1.2, label="optimal m (ours)")
        ax.axhline(0.0, color="#52514e", lw=0.8)
        ax.set_xticks(range(len(cases)))
        ax.set_xticklabels(cases, rotation=60, fontsize=6.5, ha="right")
        ax.set_ylabel("ours / paper - 1 (%)")
        ax.legend(fontsize=6.5)
    for a in axs:
        a.grid(True, color="#e6e5df", lw=0.6)
    fig.suptitle("K^-1 v'_opt residuals vs the near-O share of P_u, the emulated quadrature bias and m perturbations",
                 fontsize=9)
    fig.tight_layout()
    path = os.path.join(OUT_DIR, "summary_vopt_residuals.png")
    fig.savefig(path, dpi=130)
    plt.close(fig)
    log(f"   wrote {path}")


# ======================================================================================================
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--only", default=",".join(EXPERIMENTS), help="comma list of experiments (default: all)")
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--plot-only", action="store_true", help="no new computations: tables from the cache")
    args = ap.parse_args(argv)
    os.makedirs(OUT_DIR, exist_ok=True)
    cache = Cache(os.path.join(OUT_DIR, "cache.jsonl"))
    names = [n.strip() for n in args.only.split(",") if n.strip()]
    logpath = os.path.join(OUT_DIR, "diagnose_output.txt" if set(names) >= set(EXPERIMENTS) - {"summary"}
                           else "diagnose_" + "_".join(names) + ".txt")
    fh = open(logpath, "w")

    def log(s=""):
        print(s, flush=True)
        fh.write(s + "\n")
        fh.flush()

    t0 = time.time()
    log(f"diagnose_fig8_fig9.py  {time.strftime('%Y-%m-%d %H:%M:%S')}  workers={args.workers}")
    for n in names:
        fun = globals().get("exp_" + n)
        if fun is None:
            log(f"unknown experiment {n!r}; available: {', '.join(EXPERIMENTS)}")
            continue
        log("")
        log("=" * 110)
        log(f"{n}: {fun.__doc__.strip().splitlines()[0] if fun.__doc__ else ''}")
        log("=" * 110)
        t1 = time.time()
        fun(args, cache, log)
        log(f"   [{n}: {time.time() - t1:.0f} s]")
    log(f"total {time.time() - t0:.0f} s")
    fh.close()


if __name__ == "__main__":
    main()
