#!/usr/bin/env python3
"""Full Python reproduction of the deterministic results of Ceron et al. (IJNAMG 2025), Sect. 4.3.

Integrates the verified reference modules of this directory:
    analytical_seepage.py  -> f = K^-1 . v'_opt   (semi-analytical optimal velocity field, Sect. 4.1)
    fe_seepage.py          -> f = -grad u'_FE     (P2 FE solution of Eqs. 20-22, paper box, "zero_lb")
    limit_analysis.py      -> Gamma = min P_mr / (P_gamma + P_u)  (log-spiral mechanisms I and II, PSO + NM)
and compares with the digitised paper data (data/paper_fig{5,8,9}.csv and the *_vertices.csv polylines).

Computed:
    Fig. 5 : -J*(v'_opt) / (k_h H^2 gw^2) and J(u'_FE) / (k_h H^2 gw^2), h_w = H, alpha = 1, 2, 4, 10,
             beta = 15..90 deg in steps of 7.5 deg.
    Fig. 8 : H_crit = Gamma(H_ref) H_ref vs h_w/H in {0, 0.05, 0.1, ..., 1.0}, alpha = 1, London (beta 30, 60)
             and Israeli (beta 35, 60) panels, both seepage approaches.  Two soil-parameter sets:
               "fitted": the (c, phi) pairs of Table 1 exchanged between the panels and gamma_w = 9.8
                         (limit_analysis.FIG8_PANEL_PARAMS_FITTED; reproduces the h_w = 0 ends of Fig. 8)
               "table1": Table 1 as printed (London c = 6, phi = 32; Israeli c = 11.7, phi = 24.7), gamma_w = 9.8
             data/python_fig8.csv holds the "fitted" set (the one that reproduces the paper),
             results/reproduce_python/python_fig8_table1.csv the Table 1 set.
    Fig. 9 : Gamma vs beta, H = 5 m, c = 10 kPa, phi = 30 deg, gamma = 20 kN/m^3, h_w = H, gamma_w = 9.81,
             alpha = 1, 5, 10, both approaches; beta in {15, 20, 25, 30, 37.5, 45, 52.5, 60, 67.5, 75, 82.5, 90}
             plus the remaining 5-deg nodes of the paper polylines (35, 40, 50, 55, 65, 70, 80, 85) so that every
             comparison with the paper is made at a vertex of its plotted polyline.

Outputs (same columns as the paper_*.csv files):
    data/python_fig5.csv, data/python_fig8.csv, data/python_fig9.csv
    data/python_vs_paper_fig{5,8,9}.png
    results/reproduce_python/: comparison_fig{5,8,9}.csv (ours, paper, relative difference, mechanism ...),
        python_fig8_table1.csv, python_vs_paper_fig8_table1.png, reproduce_output.txt (log), cache.jsonl
The limit-analysis runs are cached in results/reproduce_python/cache.jsonl (keyed by the full task spec);
--fresh ignores the cache.  Diagnostic experiments: --diagnose (see DIAGNOSTICS below).

    python3 Projects/SlopeSeepageForces/scripts/reproduce_python.py [--only fig5,fig8,fig9] [--workers 2]
            [--seeds 0,1] [--quick] [--fresh] [--plot-only] [--diagnose name,...]
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
OUT_DIR = os.path.join(PROJECT_DIR, "results", "reproduce_python")
if SCRIPT_DIR not in sys.path:
    sys.path.insert(0, SCRIPT_DIR)

from analytical_seepage import AnalyticalSeepage  # noqa: E402
from fe_seepage import FESeepage  # noqa: E402
from limit_analysis import FIG8_PANEL_PARAMS_FITTED, TABLE1, stability_factor  # noqa: E402

# ======================================================================================================
# configuration
# ======================================================================================================
FIG5_ALPHAS = (1, 2, 4, 10)
FIG5_BETAS = tuple(15.0 + 7.5 * k for k in range(11))
FIG8_HW = (0.0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
FIG8_PANELS = (("London", 30.0), ("London", 60.0), ("Israeli", 35.0), ("Israeli", 60.0))
FIG8_H_REF = 10.0
FIG8_PARAM_SETS = {
    "fitted": {s: dict(FIG8_PANEL_PARAMS_FITTED[s]) for s in ("London", "Israeli")},
    "table1": {s: dict(TABLE1[s], gamma_w=9.8) for s in ("London", "Israeli")},
}
FIG9_BETAS_REQ = (15.0, 20.0, 25.0, 30.0, 37.5, 45.0, 52.5, 60.0, 67.5, 75.0, 82.5, 90.0)
FIG9_BETAS = tuple(sorted(set(FIG9_BETAS_REQ) | {35.0, 40.0, 50.0, 55.0, 65.0, 70.0, 80.0, 85.0}))
FIG9_ALPHAS = (1, 5, 10)
FIG9_DATA = dict(H=5.0, c=10.0, phi=30.0, gamma=20.0, gamma_w=9.81)
APPROACHES = ("vopt", "FE")

# plotting tokens (categorical slots 1-2 of the validated reference palette; text and paper in ink tokens)
COL = {"vopt": "#2a78d6", "FE": "#eb6834", "FE_box_m": "#1baf7a", "paper": "#52514e", "wu": "#a9a8a2", "grid": "#e6e5df",
       "text": "#0b0b0b", "text2": "#52514e", "surface": "#fcfcfb"}
LS8 = {"vopt": (0, (5, 2.2)), "FE": (0, (7, 2, 1.2, 2)), "Wu_rp025": "-"}       # Fig. 8 paper styles
LS9 = {"vopt": (0, (5, 2.2)), "FE": "-", "FE_box_m": "-"}                                          # Fig. 9 paper styles


# ======================================================================================================
# diagnostic field variants
# ======================================================================================================
class Eq31VariantSeepage(AnalyticalSeepage):
    """Zone-1 velocity of Eq. (31) with the typos of the PRINTED equation re-introduced (diagnostic only).

    Corrected (analytical_seepage.py):  v1 = a [ (1 - m s^(m-1)) h2 e_theta - (1 - s^(m-1)) h2' e_r ],
        a = k_h h2(pi - beta) gw sin(beta) / (D - C),  m = sqrt(C/D),  s = r / R_w.
    Printed Eq. (31): prefactor 1/(C - D) (sign) and C/D instead of sqrt(C/D) in the exponent AND the
    coefficient of the e_theta term.  variant = "printed" (both), "sign" (sign only), "expo" (C/D only),
    "none" (= corrected; consistency check).  Zones 2-3 are unchanged (Eqs. 32-33)."""

    def __init__(self, *args, variant="printed", **kw):
        super().__init__(*args, **kw)
        if variant not in ("printed", "sign", "expo", "none"):
            raise ValueError(variant)
        self.variant = variant

    def polar_velocity(self, r, th):
        vr, vt = super().polar_velocity(r, th)
        if self.hw <= 0.0:
            return vr, vt
        r, th = np.broadcast_arrays(np.asarray(r, float), np.asarray(th, float))
        z1 = (r > 0.0) & (r < self.Rw) & (r < self.Re)
        if np.any(z1):
            s = r[z1] / self.Rw
            t1 = np.clip(th[z1], 0.0, self.Theta)
            m = self.m
            a = self.kh * self.gamma_w * self.sb * self.h2e / (self.D - self.C)
            me = m * m if self.variant in ("printed", "expo") else m
            sign = -1.0 if self.variant in ("printed", "sign") else 1.0
            pw = me * s ** (me - 1.0) if me > 0.0 else np.zeros_like(s)
            # degenerate case m = 0: h2' = 0, so the e_r term vanishes (s > 0 keeps s^(m-1) finite)
            vr[z1] = -sign * a * (1.0 - s ** (m - 1.0)) * self.dh2(t1)
            vt[z1] = sign * a * (1.0 - pw) * self.h2(t1)
        return vr, vt


def make_field(approach, beta, H, hw, alpha, gamma_w, field_kw=None):
    """Seepage-force field for one case (None = no seepage, h_w = 0)."""
    field_kw = dict(field_kw or {})
    if hw <= 0.0 or approach == "none":
        return None
    if approach == "vopt":
        return AnalyticalSeepage(beta, H, hw, alpha, gamma_w=gamma_w, **field_kw)
    if approach.startswith("vopt_eq31_"):
        return Eq31VariantSeepage(beta, H, hw, alpha, gamma_w=gamma_w, variant=approach[len("vopt_eq31_"):],
                                  **field_kw)
    if approach == "FE":
        return FESeepage(beta, H=H, hw=hw, alpha=alpha, gamma_w=gamma_w, **field_kw)
    raise ValueError(f"unknown approach {approach!r}")


def _field_info(field):
    if field is None:
        return {}
    if isinstance(field, AnalyticalSeepage):
        return dict(m=float(field.m), degenerate=bool(field.degenerate), Jn=float(-field.Jstar_normalized()))
    if isinstance(field, FESeepage):
        return dict(Jn=float(field.J_normalized()), ndof=int(field.space.ndof), bc=str(field.bc),
                    order=int(field.order))
    return {}


# ======================================================================================================
# tasks (top level: picklable for multiprocessing)
# ======================================================================================================
def _jsonable(v):
    if isinstance(v, dict):
        return {str(k): _jsonable(x) for k, x in v.items()}
    if isinstance(v, (list, tuple)):
        return [_jsonable(x) for x in v]
    if isinstance(v, (np.floating, float)):
        return float(v)
    if isinstance(v, (np.integer,)):
        return int(v)
    if isinstance(v, np.bool_):
        return bool(v)
    if isinstance(v, np.ndarray):
        return [_jsonable(x) for x in v.tolist()]
    return v


def run_task(spec):
    """One task: Fig. 5 functionals, or one limit analysis (Fig. 8 / Fig. 9 / diagnostics)."""
    t0 = time.time()
    if spec["kind"] == "fig5":
        an = AnalyticalSeepage(spec["beta"], 1.0, 1.0, spec["alpha"])
        t1 = time.time()
        fe = FESeepage(spec["beta"], H=1.0, hw=1.0, alpha=spec["alpha"], **spec.get("field_kw", {}))
        return _jsonable(dict(spec=spec, minus_Jstar_opt=-an.Jstar_normalized(), J_FE=fe.J_normalized(),
                              m=an.m, degenerate=an.degenerate, ndof=fe.space.ndof, t_an=t1 - t0,
                              t_fe=time.time() - t1, time=time.time() - t0))
    H = spec["H"]
    hw = spec["hw_over_H"] * H
    field = make_field(spec["approach"], spec["beta"], H, hw, spec["alpha"], spec["gamma_w"], spec.get("field_kw"))
    t_field = time.time() - t0
    r = stability_factor(spec["beta"], H, spec["c"], spec["phi"], spec["gamma"], spec["gamma_w"], field,
                         seeds=tuple(spec["seeds"]), **spec.get("la_kw", {}))
    out = dict(spec=spec, Gamma=r["Gamma"], Hcrit=r["Hcrit"], mechanism=r["mechanism"], params=r["params"],
               pu_method=r["pu_method"], t_field=t_field, t_la=r["time"], time=time.time() - t0,
               field_info=_field_info(field))
    if r["mechanism"] is not None:
        out.update(x=r["x"], P_mr=r["P_mr"], P_gamma=r["P_gamma"], P_u=r["P_u"], L=r["L"], r0=r["r0"],
                   seed_spread=r["seed_spread"], at_search_bound=r["at_search_bound"],
                   d_max_reached=r["d_max_reached"], best_by_class=r["best_by_class"], A=r["A"], B=r["B"],
                   C=r["C"])
    return _jsonable(out)


def spec_key(spec):
    return json.dumps(spec, sort_keys=True)


class Cache:
    def __init__(self, path, fresh=False):
        self.path = path
        self.data = {}
        if not fresh:
            self.refresh()

    def refresh(self):
        """(re-)read the cache file: picks up results appended by another process running concurrently"""
        if os.path.exists(self.path):
            with open(self.path) as fh:
                for line in fh:
                    line = line.strip()
                    if line:
                        try:
                            d = json.loads(line)
                        except json.JSONDecodeError:
                            continue
                        self.data[spec_key(d["spec"])] = d

    def get(self, spec):
        return self.data.get(spec_key(spec))

    def put(self, res):
        self.data[spec_key(res["spec"])] = res
        with open(self.path, "a") as fh:
            fh.write(json.dumps(res) + "\n")


def _cost(spec):
    """rough relative cost, used to start the slow (flat-slope) cases first"""
    if spec["kind"] == "fig5":
        return 1.0
    return 10.0 * len(spec["seeds"]) * (1.0 + 400.0 / spec["beta"] ** 2) * (1.4 if spec["approach"] != "FE" else 1.0)


def run_specs(specs, workers, cache, log, label, compute=True):
    """Run all specs not in the cache (in parallel), return {key: result}."""
    todo = [s for s in specs if cache.get(s) is None]
    t0 = time.time()
    if todo and compute:
        log(f"[{label}] {len(specs)} tasks, {len(specs) - len(todo)} cached, running {len(todo)} on {workers} worker(s)")
        todo.sort(key=_cost, reverse=True)
        done = 0
        if workers <= 1:
            it = map(run_task, todo)
            for res in it:
                cache.put(res)
                done += 1
                _progress(log, label, done, len(todo), res, t0)
        else:
            ctx = mp.get_context("fork")
            with ctx.Pool(workers, maxtasksperchild=25) as pool:
                for res in pool.imap_unordered(run_task, todo, chunksize=1):
                    cache.put(res)
                    done += 1
                    _progress(log, label, done, len(todo), res, t0)
    elif todo:
        log(f"[{label}] {len(todo)} of {len(specs)} tasks missing from the cache (--plot-only: skipped)")
    wall = time.time() - t0
    out = {spec_key(s): cache.get(s) for s in specs if cache.get(s) is not None}
    cpu = sum(r["time"] for r in out.values())
    log(f"[{label}] wall {wall:.1f} s for the new tasks; summed task time of all {len(out)} results {cpu:.1f} s")
    return out, wall


def _progress(log, label, done, n, res, t0):
    s = res["spec"]
    if s["kind"] == "fig5":
        msg = f"alpha={s['alpha']} beta={s['beta']}: -J*={res['minus_Jstar_opt']:.5f} J_FE={res['J_FE']:.5f}"
    else:
        tag = s.get("soil", f"alpha={s['alpha']}")
        msg = (f"{tag} beta={s['beta']} hw/H={s['hw_over_H']} {s['approach']}: Gamma={res['Gamma']:.5f} "
               f"Hcrit={res['Hcrit']:.4f} {res['mechanism']}")
    log(f"   [{label} {done}/{n} {time.time() - t0:7.1f} s] {msg} ({res['time']:.1f} s)")


# ======================================================================================================
# paper data
# ======================================================================================================
def _read_csv(name):
    with open(os.path.join(DATA_DIR, name)) as fh:
        return list(csv.DictReader(fh))


def paper_fig5_polylines():
    out = {}
    for r in _read_csv("paper_fig5_vertices.csv"):
        out.setdefault((int(r["alpha"]), r["curve"]), []).append((float(r["beta_deg"]), float(r["value"])))
    return {k: np.array(sorted(v)) for k, v in out.items()}


def paper_fig8_polylines():
    out = {}
    for r in _read_csv("paper_fig8_vertices.csv"):
        out.setdefault((r["soil"], float(r["beta_deg"]), r["curve"]), []).append(
            (float(r["hw_over_H"]), float(r["Hcrit_m"]), int(r["visible"])))
    return {k: np.array(sorted(v)) for k, v in out.items()}


def paper_fig9_polylines():
    out = {}
    for r in _read_csv("paper_fig9_vertices.csv"):
        out.setdefault((int(r["alpha"]), r["curve"]), []).append(
            (float(r["beta_deg"]), float(r["Gamma"]), int(r["visible"]), float(r["u_Gamma"])))
    return {k: np.array(sorted(v)) for k, v in out.items()}


def interp_polyline(P, x, log=False):
    """value of the paper polyline P (rows x, y, visible[, u]) at x, and a status string:
    'vertex' / 'segment' (both ends visible), 'extrapolated' (a hidden end, i.e. a clipped part read by
    extrapolation of the visible segment), 'none' (outside the polyline)."""
    xs, ys = P[:, 0], P[:, 1]
    vis = P[:, 2] > 0.5 if P.shape[1] > 2 else np.ones(len(xs), bool)
    if x < xs[0] - 1e-9 or x > xs[-1] + 1e-9:
        return np.nan, "none"
    j = int(np.argmin(np.abs(xs - x)))
    if abs(xs[j] - x) < 1e-9:
        return float(ys[j]), ("vertex" if vis[j] else "extrapolated")
    i = int(np.searchsorted(xs, x)) - 1
    w = (x - xs[i]) / (xs[i + 1] - xs[i])
    if log:
        v = 10.0 ** ((1 - w) * np.log10(ys[i]) + w * np.log10(ys[i + 1]))
    else:
        v = (1 - w) * ys[i] + w * ys[i + 1]
    return float(v), ("segment" if vis[i] and vis[i + 1] else "extrapolated")


def _rel(ours, paper):
    if not (np.isfinite(ours) and np.isfinite(paper)) or paper == 0:
        return np.nan
    return ours / paper - 1.0


# ======================================================================================================
# figure drivers
# ======================================================================================================
def fig5(args, cache, log):
    alphas = (1, 10) if args.quick else FIG5_ALPHAS
    betas = (15.0, 45.0, 90.0) if args.quick else FIG5_BETAS
    specs = [dict(kind="fig5", alpha=a, beta=b, field_kw={}) for a in alphas for b in betas]
    res, wall = run_specs(specs, args.workers, cache, log, "fig5", compute=not args.plot_only)
    paper = paper_fig5_polylines()
    rows, comp = [], []
    for s in specs:
        r = res.get(spec_key(s))
        if r is None:
            continue
        for curve in ("J_FE", "minus_Jstar_opt"):
            v = r[curve]
            rows.append(dict(alpha=s["alpha"], beta_deg=_fmt_beta(s["beta"]), curve=curve, value=f"{v:.5f}"))
            pv, st = interp_polyline(paper[(s["alpha"], curve)], s["beta"])
            comp.append(dict(alpha=s["alpha"], beta_deg=s["beta"], curve=curve, ours=v, paper=pv, status=st,
                             rel_diff=_rel(v, pv)))
    _write_csv(os.path.join(DATA_DIR, "python_fig5.csv"), ["alpha", "beta_deg", "curve", "value"], rows)
    _write_csv(os.path.join(OUT_DIR, "comparison_fig5.csv"),
               ["alpha", "beta_deg", "curve", "ours", "paper", "status", "rel_diff"], comp)
    log("")
    log("Fig. 5: normalised functionals, h_w = H (ours / paper vertex polyline, rel = ours/paper - 1)")
    log(f"{'alpha':>5} {'beta':>6} | {'-J*(vopt)':>9} {'paper':>8} {'rel':>9} | {'J(uFE)':>8} {'paper':>8} {'rel':>9}")
    for a in alphas:
        for b in betas:
            c1 = [c for c in comp if c["alpha"] == a and c["beta_deg"] == b and c["curve"] == "minus_Jstar_opt"]
            c2 = [c for c in comp if c["alpha"] == a and c["beta_deg"] == b and c["curve"] == "J_FE"]
            if c1 and c2:
                c1, c2 = c1[0], c2[0]
                log(f"{a:5d} {b:6.1f} | {c1['ours']:9.5f} {c1['paper']:8.5f} {c1['rel_diff']:+9.2e} | "
                    f"{c2['ours']:8.5f} {c2['paper']:8.5f} {c2['rel_diff']:+9.2e}")
    _curve_stats(log, comp, ("alpha", "curve"))
    if not args.no_plot:
        plot_fig5(comp, alphas, paper, os.path.join(DATA_DIR, "python_vs_paper_fig5.png"))
    return comp, wall


def fig8_specs(pset, args):
    hws = (0.0, 0.5, 1.0) if args.quick else FIG8_HW
    specs = []
    for soil, beta in FIG8_PANELS:
        p = FIG8_PARAM_SETS[pset][soil]
        base = dict(kind="la", fig="fig8", pset=pset, soil=soil, beta=beta, alpha=1, H=FIG8_H_REF, c=p["c"],
                    phi=p["phi"], gamma=p["gamma"], gamma_w=p["gamma_w"], seeds=list(args.seeds), field_kw={},
                    la_kw={})
        for hw in hws:
            if hw == 0.0:
                specs.append(dict(base, hw_over_H=0.0, approach="none"))
            else:
                for ap in APPROACHES:
                    specs.append(dict(base, hw_over_H=hw, approach=ap))
    return specs


def fig8(args, cache, log, pset="fitted"):
    specs = fig8_specs(pset, args)
    res, wall = run_specs(specs, args.workers, cache, log, f"fig8-{pset}", compute=not args.plot_only)
    paper = paper_fig8_polylines()
    rows, comp = [], []
    for s in specs:
        r = res.get(spec_key(s))
        if r is None:
            continue
        for curve in (APPROACHES if s["approach"] == "none" else (s["approach"],)):
            v = r["Hcrit"]
            rows.append(dict(soil=s["soil"], beta_deg=int(s["beta"]), curve=curve, hw_over_H=s["hw_over_H"],
                             Hcrit_m=f"{v:.3f}" if np.isfinite(v) else "inf"))
            pv, st = interp_polyline(paper[(s["soil"], s["beta"], curve)], s["hw_over_H"], log=True)
            comp.append(dict(pset=pset, soil=s["soil"], beta_deg=s["beta"], curve=curve, hw_over_H=s["hw_over_H"],
                             ours=v, paper=pv, status=st, rel_diff=_rel(v, pv), mechanism=r["mechanism"],
                             eta_or_dH=_s_param(r), theta1=(r["params"] or {}).get("theta1", np.nan),
                             theta2=(r["params"] or {}).get("theta2", np.nan), L_over_H=r.get("L", np.nan) / s["H"],
                             at_bound=";".join(r.get("at_search_bound", []) or []),
                             seed_spread=r.get("seed_spread", np.nan), time_s=r["time"]))
    if pset == "fitted":
        _write_csv(os.path.join(DATA_DIR, "python_fig8.csv"), ["soil", "beta_deg", "curve", "hw_over_H", "Hcrit_m"],
                   rows)
    else:
        _write_csv(os.path.join(OUT_DIR, f"python_fig8_{pset}.csv"),
                   ["soil", "beta_deg", "curve", "hw_over_H", "Hcrit_m"], rows)
    fields = ["pset", "soil", "beta_deg", "curve", "hw_over_H", "ours", "paper", "status", "rel_diff", "mechanism",
              "eta_or_dH", "theta1", "theta2", "L_over_H", "at_bound", "seed_spread", "time_s"]
    _write_csv(os.path.join(OUT_DIR, f"comparison_fig8_{pset}.csv"), fields, comp)
    log("")
    ps = FIG8_PARAM_SETS[pset]
    log(f"Fig. 8, parameter set '{pset}': " + "; ".join(f"{k} panel c={v['c']} phi={v['phi']} gamma={v['gamma']} "
                                                       f"gamma_w={v['gamma_w']}" for k, v in ps.items()))
    log("   H_crit (m) ours vs paper polyline (log-interpolated), rel = ours/paper - 1; mechanism (eta | d/H)")
    hws = sorted({c["hw_over_H"] for c in comp})
    for soil, beta in FIG8_PANELS:
        for curve in APPROACHES:
            log(f"   {soil} beta={beta:.0f} {curve}:")
            for hw in hws:
                c = [q for q in comp if q["soil"] == soil and q["beta_deg"] == beta and q["curve"] == curve
                     and q["hw_over_H"] == hw]
                if c:
                    c = c[0]
                    log(f"      hw/H={hw:4.2f}  ours={c['ours']:9.3f}  paper={c['paper']:9.3f} ({c['status']:>12})  "
                        f"rel={c['rel_diff']:+8.4f}  {c['mechanism']} {c['eta_or_dH']:.4f}")
    _curve_stats(log, comp, ("soil", "beta_deg", "curve"))
    if not args.no_plot:
        path = os.path.join(DATA_DIR, "python_vs_paper_fig8.png") if pset == "fitted" else \
            os.path.join(OUT_DIR, f"python_vs_paper_fig8_{pset}.png")
        plot_fig8(comp, paper, path, pset)
    return comp, wall


def fig9_specs(args):
    betas = (30.0, 60.0, 90.0) if args.quick else FIG9_BETAS
    alphas = (1,) if args.quick else FIG9_ALPHAS
    specs = []
    for a in alphas:
        for b in betas:
            for ap in APPROACHES:
                specs.append(dict(kind="la", fig="fig9", alpha=a, beta=b, approach=ap, hw_over_H=1.0,
                                  seeds=list(args.seeds), field_kw={}, la_kw={}, **FIG9_DATA))
            # FE variant with the FE box fixed in metres (50/10/30 m: the Fig. 4 box of the H = 1 m slope), see D1
            specs.append(dict(kind="la", fig="fig9", alpha=a, beta=b, approach="FE", hw_over_H=1.0,
                              seeds=list(args.seeds), field_kw=box_in_metres(FIG9_DATA["H"]), la_kw={}, **FIG9_DATA))
    return specs


def _curve9(spec):
    """curve label of a Fig. 9 spec: vopt, FE (box scaled with H) or FE_box_m (box fixed in metres)"""
    if spec["approach"] == "FE" and spec.get("field_kw"):
        return "FE_box_m"
    return spec["approach"]


def fig9(args, cache, log):
    specs = fig9_specs(args)
    res, wall = run_specs(specs, args.workers, cache, log, "fig9", compute=not args.plot_only)
    paper = paper_fig9_polylines()
    rows, comp = [], []
    for s in specs:
        r = res.get(spec_key(s))
        if r is None:
            continue
        v = r["Gamma"]
        curve = _curve9(s)
        rows.append(dict(alpha=s["alpha"], beta_deg=_fmt_beta(s["beta"]), curve=curve,
                         Gamma=f"{v:.4f}" if np.isfinite(v) else "inf"))
        P = paper[(s["alpha"], "vopt" if curve == "vopt" else "FE")]
        pv, st = interp_polyline(P, s["beta"])
        if st == "none":
            st = "offplot(>5)"      # left of the first (hidden) node: the paper curve is above the plot range
        comp.append(dict(alpha=s["alpha"], beta_deg=s["beta"], curve=curve, ours=v, paper=pv, status=st,
                         rel_diff=_rel(v, pv), node=bool(abs(s["beta"] / 5.0 - round(s["beta"] / 5.0)) < 1e-9),
                         mechanism=r["mechanism"], eta_or_dH=_s_param(r),
                         theta1=(r["params"] or {}).get("theta1", np.nan),
                         theta2=(r["params"] or {}).get("theta2", np.nan), P_gamma=r.get("P_gamma", np.nan),
                         P_u=r.get("P_u", np.nan), P_mr=r.get("P_mr", np.nan),
                         seed_spread=r.get("seed_spread", np.nan), time_s=r["time"]))
    _write_csv(os.path.join(DATA_DIR, "python_fig9.csv"), ["alpha", "beta_deg", "curve", "Gamma"], rows)
    fields = ["alpha", "beta_deg", "curve", "ours", "paper", "status", "rel_diff", "node", "mechanism", "eta_or_dH",
              "theta1", "theta2", "P_gamma", "P_u", "P_mr", "seed_spread", "time_s"]
    _write_csv(os.path.join(OUT_DIR, "comparison_fig9.csv"), fields, comp)
    log("")
    log("Fig. 9: Gamma, H = 5, c = 10, phi = 30, gamma = 20, h_w = H, gamma_w = 9.81; ours vs paper polyline")
    log("   ('vertex' = 5-deg node of the paper polyline, 'segment' = between nodes (chord), 'extrapolated' =")
    log("    hidden node above the plot, 'offplot(>5)' = no paper curve visible there)")
    for a in sorted({c["alpha"] for c in comp}):
        log(f"   alpha = {a}")
        log(f"   {'beta':>6} | {'FE ours':>8} {'paper':>8} {'rel':>8} {'mech':>9} | {'vopt ours':>9} {'paper':>8} "
            f"{'rel':>8} {'mech':>9} | {'FE_box_m':>8} {'rel':>8} {'mech':>9}")
        for b in sorted({c["beta_deg"] for c in comp}):
            c1 = [c for c in comp if c["alpha"] == a and c["beta_deg"] == b and c["curve"] == "FE"]
            c2 = [c for c in comp if c["alpha"] == a and c["beta_deg"] == b and c["curve"] == "vopt"]
            c3 = [c for c in comp if c["alpha"] == a and c["beta_deg"] == b and c["curve"] == "FE_box_m"]
            if not (c1 and c2):
                continue
            c1, c2 = c1[0], c2[0]
            c3 = c3[0] if c3 else dict(ours=np.nan, rel_diff=np.nan, mechanism=None)
            log(f"   {b:6.1f} | {c1['ours']:8.4f} {c1['paper']:8.4f} {c1['rel_diff']:+8.4f} "
                f"{_mtag(c1):>9} | {c2['ours']:9.4f} {c2['paper']:8.4f} {c2['rel_diff']:+8.4f} {_mtag(c2):>9} | "
                f"{c3['ours']:8.4f} {c3['rel_diff']:+8.4f} {_mtag(c3):>9}   [{c1['status']}/{c2['status']}]")
    _curve_stats(log, comp, ("alpha", "curve"))
    if not args.no_plot:
        plot_fig9(comp, paper, os.path.join(DATA_DIR, "python_vs_paper_fig9.png"))
    return comp, wall


def _mtag(c):
    if c["mechanism"] is None:
        return "-"
    return f"{c['mechanism']}{c['eta_or_dH']:.3f}"


def _s_param(r):
    p = r.get("params") or {}
    return p.get("eta", p.get("d_over_H", np.nan))


def _fmt_beta(b):
    return f"{b:g}"


def _curve_stats(log, comp, keys):
    """per-curve statistics of rel_diff on points with a visible paper value"""
    log("   per-curve relative differences (points with a visible paper value: vertex/segment):")
    groups = {}
    for c in comp:
        groups.setdefault(tuple(c[k] for k in keys), []).append(c)
    for g, cs in groups.items():
        rel = np.array([c["rel_diff"] for c in cs if c["status"] in ("vertex", "segment") and np.isfinite(c["rel_diff"])])
        if len(rel) == 0:
            continue
        i = int(np.argmax(np.abs(rel)))
        log(f"      {' '.join(str(x) for x in g):>24}: n={len(rel):2d}  mean={rel.mean():+.4f}  "
            f"rms={np.sqrt(np.mean(rel ** 2)):.4f}  max|rel|={abs(rel[i]):.4f} ({rel[i]:+.4f})")


def _write_csv(path, fields, rows):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.6g}" if isinstance(v, float) else v) for k, v in r.items()})


# ======================================================================================================
# plots
# ======================================================================================================
def _mpl():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 9, "axes.edgecolor": COL["text2"], "axes.labelcolor": COL["text"],
                         "xtick.color": COL["text2"], "ytick.color": COL["text2"], "axes.titlesize": 10,
                         "figure.facecolor": "white", "axes.facecolor": "white", "legend.frameon": False})
    return plt


def _style_ax(ax):
    ax.grid(True, color=COL["grid"], lw=0.7)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def plot_fig5(comp, alphas, paper, path):
    plt = _mpl()
    n = len(alphas)
    nc = 2 if n > 1 else 1
    nr = int(np.ceil(n / nc))
    fig, axs = plt.subplots(nr, nc, figsize=(7.2, 2.7 * nr + 0.6), squeeze=False)
    for k, a in enumerate(alphas):
        ax = axs[k // nc][k % nc]
        _style_ax(ax)
        for curve, ls, lab in (("minus_Jstar_opt", "-", r"$-J^*(v'_{opt})$"), ("J_FE", (0, (5, 2.2)), r"$J(u'_{FE})$")):
            P = paper[(a, curve)]
            ax.plot(P[:, 0], P[:, 1], ls=ls, color=COL["paper"], lw=0.8, label=f"paper {lab}", zorder=4)
            cs = sorted([c for c in comp if c["alpha"] == a and c["curve"] == curve], key=lambda c: c["beta_deg"])
            col = COL["vopt"] if curve == "minus_Jstar_opt" else COL["FE"]
            ax.plot([c["beta_deg"] for c in cs], [c["ours"] for c in cs], ls=ls, color=col, lw=2.0, alpha=0.85,
                    marker="o", ms=3.5, label=f"ours {lab}")
            mx = max(abs(c["rel_diff"]) for c in cs)
            ax.text(0.98, 0.04 + (0.0 if curve == "J_FE" else 0.09), f"{lab}: max |rel| = {100 * mx:.2f} %",
                    transform=ax.transAxes, ha="right", va="bottom", fontsize=7.5, color=COL["text2"])
        ax.set_title(rf"$\alpha$ = {a}", loc="left", color=COL["text"])
        ax.set_xlim(15, 90)
        ax.set_xticks([15, 30, 45, 60, 75, 90])
        ax.set_ylim(0, None)
        if k // nc == nr - 1:
            ax.set_xlabel(r"slope inclination $\beta$ (deg)")
        if k % nc == 0:
            ax.set_ylabel(r"functional / ($k_h H^2 \gamma_w^2$)")
    for k in range(n, nr * nc):
        axs[k // nc][k % nc].axis("off")
    h, l = axs[0][0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=4, fontsize=8)
    fig.suptitle(r"Fig. 5, $h_w = H$: Python reproduction (thick) vs paper (thin grey)", color=COL["text"], fontsize=10)
    fig.tight_layout(rect=(0, 0.06, 1, 0.97))
    fig.savefig(path, dpi=150)
    plt.close(fig)


def plot_fig8(comp, paper, path, pset):
    plt = _mpl()
    fig, axs = plt.subplots(2, 2, figsize=(9.0, 7.0), gridspec_kw=dict(height_ratios=(3, 1.25)), sharex=True)
    for j, soil in enumerate(("London", "Israeli")):
        ax, axr = axs[0][j], axs[1][j]
        _style_ax(ax)
        _style_ax(axr)
        betas = [b for s, b in FIG8_PANELS if s == soil]
        for beta in betas:
            P = paper[(soil, beta, "Wu_rp025")]
            ax.plot(P[:, 0], P[:, 1], ls="-", color=COL["wu"], lw=0.9, label="paper: Wu et al. $r_p$ = 0.25"
                    if beta == betas[0] else None)
            for curve in APPROACHES:
                P = paper[(soil, beta, curve)]
                vis = P[:, 2] > 0.5
                ax.plot(P[:, 0], P[:, 1], ls=LS8[curve], color=COL["paper"], lw=0.9, zorder=4,
                        label=f"paper: {_lab(curve)}" if beta == betas[0] else None)
                if np.any(~vis):
                    ax.plot(P[~vis, 0], P[~vis, 1], ls="none", marker="x", ms=4, color=COL["paper"],
                            label="paper: clipped end (extrapolated)" if beta == betas[0] and curve == "FE" else None)
                cs = sorted([c for c in comp if c["soil"] == soil and c["beta_deg"] == beta and c["curve"] == curve],
                            key=lambda c: c["hw_over_H"])
                x = np.array([c["hw_over_H"] for c in cs])
                y = np.array([c["ours"] for c in cs])
                fin = np.isfinite(y)
                ax.plot(x[fin], y[fin], ls=LS8[curve], color=COL[curve], lw=2.0, marker="o", ms=3.5, alpha=0.9,
                        label=f"ours: {_lab(curve)}" if beta == betas[0] else None)
                rel = np.array([c["rel_diff"] for c in cs])
                ok = np.isfinite(rel)
                axr.plot(x[ok], 100 * rel[ok], ls=LS8[curve], color=COL[curve], lw=1.4,
                         marker="o" if beta == betas[0] else "s", ms=3.5, mfc="white" if beta != betas[0] else COL[curve],
                         label=rf"{_lab(curve)}, $\beta$ = {beta:.0f}")
            ax.annotate(rf"$\beta$ = {beta:.0f}$\degree$", xy=(0.62, _interp_y(comp, soil, beta, "vopt", 0.6)),
                        xytext=(6, 8), textcoords="offset points", fontsize=8.5, color=COL["text"])
        ax.set_yscale("log")
        ax.set_ylim(3.0, 260.0)
        ax.set_xlim(0, 1)
        ax.set_title(f"{soil} clay panel", loc="left", color=COL["text"])
        if j == 0:
            ax.set_ylabel(r"$H_{crit}$ (m)")
            axr.set_ylabel("ours / paper - 1 (%)")
        axr.axhline(0.0, color=COL["text2"], lw=0.8)
        axr.set_xlabel(r"$h_w / H$")
        axr.legend(fontsize=7, ncol=2, loc="best")
    h, l = axs[0][0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=3, fontsize=8)
    ps = FIG8_PARAM_SETS[pset]
    sub = "; ".join(f"{k} panel: c={v['c']}, phi={v['phi']}" for k, v in ps.items())
    fig.suptitle(f"Fig. 8, parameter set '{pset}' ({sub}; gamma=18, gamma_w={ps['London']['gamma_w']}, alpha=1)",
                 color=COL["text"], fontsize=9)
    fig.tight_layout(rect=(0, 0.07, 1, 0.97))
    fig.savefig(path, dpi=150)
    plt.close(fig)


def _interp_y(comp, soil, beta, curve, x):
    cs = sorted([c for c in comp if c["soil"] == soil and c["beta_deg"] == beta and c["curve"] == curve
                 and np.isfinite(c["ours"])], key=lambda c: c["hw_over_H"])
    if not cs:
        return 10.0
    return float(np.interp(x, [c["hw_over_H"] for c in cs], [c["ours"] for c in cs]))


def _lab(curve):
    return {"vopt": r"$K^{-1} v'_{opt}$", "FE": r"$-\nabla u'_{FE}$",
            "FE_box_m": r"$-\nabla u'_{FE}$, box 50/10/30 m"}[curve]


def plot_fig9(comp, paper, path):
    plt = _mpl()
    alphas = sorted({c["alpha"] for c in comp})
    fig, axs = plt.subplots(2, len(alphas), figsize=(3.3 * len(alphas) + 0.6, 6.2),
                            gridspec_kw=dict(height_ratios=(3, 1.3)), sharex=True, squeeze=False)
    for j, a in enumerate(alphas):
        ax, axr = axs[0][j], axs[1][j]
        _style_ax(ax)
        _style_ax(axr)
        for curve in ("FE", "vopt"):
            P = paper[(a, curve)]
            vis = P[:, 2] > 0.5
            ax.plot(P[:, 0], P[:, 1], ls=LS9[curve], color=COL["paper"], lw=0.9, label=f"paper: {_lab(curve)}", zorder=4)
            ax.plot(P[vis, 0], P[vis, 1], ls="none", marker="o", ms=3, mfc="white", mec=COL["paper"], mew=0.8, zorder=5)
        for curve in ("FE", "vopt", "FE_box_m"):
            cs = sorted([c for c in comp if c["alpha"] == a and c["curve"] == curve], key=lambda c: c["beta_deg"])
            if not cs:
                continue
            x = np.array([c["beta_deg"] for c in cs])
            y = np.array([c["ours"] for c in cs])
            ax.plot(x, y, ls=LS9[curve], color=COL[curve], lw=2.0 if curve != "FE_box_m" else 1.4, alpha=0.9,
                    label=f"ours: {_lab(curve)}")
            ax.plot(x, y, ls="none", marker="o", ms=3.0, color=COL[curve])
            rel = np.array([c["rel_diff"] for c in cs])
            ok = np.isfinite(rel) & np.array([c["status"] == "vertex" for c in cs])   # paper nodes only (no chords)
            axr.plot(x[ok], 100 * rel[ok], ls=LS9[curve], color=COL[curve], lw=1.4 if curve != "FE_box_m" else 1.1,
                     marker="o", ms=3, label=_lab(curve))
        ax.set_ylim(0, 5)
        ax.set_xlim(15, 90)
        ax.set_xticks([15, 30, 45, 60, 75, 90])
        ax.set_title(rf"$\alpha$ = {a}", loc="left", color=COL["text"])
        axr.axhline(0.0, color=COL["text2"], lw=0.8)
        axr.set_xlabel(r"slope inclination $\beta$ (deg)")
        if j == 0:
            ax.set_ylabel(r"stability factor $\Gamma$")
            axr.set_ylabel("ours / paper - 1 (%)\nat the paper's 5-deg nodes")
            axr.legend(fontsize=6.5)
    h, l = axs[0][0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=3, fontsize=8)
    fig.suptitle(r"Fig. 9: H = 5 m, c = 10 kPa, $\varphi$ = 30$\degree$, $\gamma$ = 20, $h_w = H$, $\gamma_w$ = 9.81"
                 " (ours thick, paper thin grey)", color=COL["text"], fontsize=9.5)
    fig.tight_layout(rect=(0, 0.06, 1, 0.96))
    fig.savefig(path, dpi=150)
    plt.close(fig)


# ======================================================================================================
# diagnostics (targeted experiments, see the report in the log)
# ======================================================================================================
PAPER_BOX_M = dict(left=50.0, right=10.0, depth=30.0)   # FE box of Fig. 4 in metres (for H = 1 m: 50/10/30 H)


def box_in_metres(H):
    """FE box kept at 50 m / 10 m / 30 m (left of O / right of the toe / below the toe) for a slope of height H,
    i.e. NOT scaled with H (hypothesis for Fig. 9, H = 5 m: 10 H / 2 H / 6 H)."""
    return {k: v / H for k, v in PAPER_BOX_M.items()}


def _fig9_base(args, alpha, beta, approach, **kw):
    d = dict(kind="la", fig="fig9", alpha=alpha, beta=beta, approach=approach, hw_over_H=1.0,
             seeds=list(args.seeds), field_kw={}, la_kw={}, **FIG9_DATA)
    d.update(kw)
    return d


def _fig8_base(args, pset, soil, beta, hw, approach, **kw):
    p = FIG8_PARAM_SETS[pset][soil]
    d = dict(kind="la", fig="fig8", pset=pset, soil=soil, beta=beta, alpha=1, H=FIG8_H_REF, c=p["c"], phi=p["phi"],
             gamma=p["gamma"], gamma_w=p["gamma_w"], seeds=list(args.seeds), field_kw={}, la_kw={},
             hw_over_H=hw, approach=approach)
    d.update(kw)
    return d


def _get(res, spec):
    r = res.get(spec_key(spec))
    return r if r is not None else {"Gamma": np.nan, "Hcrit": np.nan, "mechanism": None, "params": None}


DIAGNOSTICS = ("fe_box_m", "fe_box_m_fig8", "eq31", "gamma_w", "fe_bc", "classes", "beta15", "robust")
ROBUST_POINTS = (("Israeli", 35.0, 0.1, "vopt"), ("Israeli", 35.0, 0.2, "vopt"), ("Israeli", 35.0, 0.2, "FE"),
                 ("Israeli", 35.0, 0.5, "FE"), ("London", 30.0, 0.1, "vopt"), ("London", 30.0, 0.05, "FE"))


def diagnostics(names, args, cache, log):
    """Targeted experiments for the curves that disagree with the paper (see the module docstring)."""
    if "all" in names:
        names = list(DIAGNOSTICS)
    p9 = paper_fig9_polylines()
    p8 = paper_fig8_polylines()
    compute = not args.plot_only
    for name in names:
        log("")
        log("=" * 110)
        if name == "fe_box_m":
            log("D1) Fig. 9, -grad u'_FE with the FE box fixed in METRES (50 m left of O, 10 m right of the toe, 30 m")
            log("    below the toe, i.e. the Fig. 4 box for H = 1 m) instead of scaled with H (H = 5 m: 10 H / 2 H / 6 H)")
            log("=" * 110)
            betas = [b for b in FIG9_BETAS if abs(b / 5.0 - round(b / 5.0)) < 1e-9 and b >= 20.0]
            specs, rows = [], []
            for a in FIG9_ALPHAS:
                for b in betas:
                    s0 = _fig9_base(args, a, b, "FE")
                    s1 = _fig9_base(args, a, b, "FE", field_kw=box_in_metres(FIG9_DATA["H"]))
                    specs += [s0, s1]
                    rows.append((a, b, s0, s1))
            run_specs([r[3] for r in rows], args.workers, cache, log, "D1 (box in m)", compute)
            cache.refresh()
            res, _ = run_specs(specs, args.workers, cache, log, "D1", compute)
            log(f"{'alpha':>5} {'beta':>5} | {'paper':>7} {'status':>12} | {'box ~ H':>8} {'rel':>8} | {'box in m':>8} "
                f"{'rel':>8} {'mech':>8}")
            st = {"H": [], "m": []}
            out = []
            for a, b, s0, s1 in rows:
                pv, sts = interp_polyline(p9[(a, "FE")], b)
                r0, r1 = _get(res, s0), _get(res, s1)
                e0, e1 = _rel(r0["Gamma"], pv), _rel(r1["Gamma"], pv)
                if sts == "vertex":
                    st["H"].append(e0)
                    st["m"].append(e1)
                log(f"{a:5d} {b:5.0f} | {pv:7.4f} {sts:>12} | {r0['Gamma']:8.4f} {e0:+8.4f} | {r1['Gamma']:8.4f} "
                    f"{e1:+8.4f} {_mtag(dict(mechanism=r1['mechanism'], eta_or_dH=_s_param(r1))):>8}")
                out.append(dict(alpha=a, beta_deg=b, paper=pv, status=sts, Gamma_box_scaled=r0["Gamma"],
                                rel_box_scaled=e0, Gamma_box_metres=r1["Gamma"], rel_box_metres=e1,
                                mechanism_box_metres=r1["mechanism"]))
            for k, lab in (("H", "box scaled with H"), ("m", "box fixed in metres")):
                e = np.array([x for x in st[k] if np.isfinite(x)])
                if len(e):
                    log(f"   {lab:>22}: n={len(e)} mean rel {e.mean():+.4f}, rms {np.sqrt(np.mean(e ** 2)):.4f}, "
                        f"max |rel| {np.abs(e).max():.4f}")
            _write_csv(os.path.join(OUT_DIR, "diag_fig9_fe_box_metres.csv"), list(out[0].keys()), out)
            if not args.no_plot:
                _plot_diag_box(out, p9, os.path.join(OUT_DIR, "diag_fig9_fe_box_metres.png"))
        elif name == "fe_box_m_fig8":
            log("D1b) Counter-check on Fig. 8: had the FE box been fixed in metres AND H_crit been found at the actual")
            log("     height, Gamma(H = H_crit,paper) with the 50/10/30 m box would be 1.  Gamma at H = H_crit,paper:")
            log("=" * 110)
            pts = (("London", 30.0, 1.0), ("London", 60.0, 1.0), ("Israeli", 60.0, 0.5))
            specs, rows = [], []
            for soil, beta, hw in pts:
                pv, _ = interp_polyline(p8[(soil, beta, "FE")], hw, log=True)
                s0 = _fig8_base(args, "fitted", soil, beta, hw, "FE", H=round(pv, 3))
                s1 = _fig8_base(args, "fitted", soil, beta, hw, "FE", H=round(pv, 3), field_kw=box_in_metres(round(pv, 3)))
                specs += [s0, s1]
                rows.append((soil, beta, hw, pv, s0, s1))
            res, _ = run_specs(specs, args.workers, cache, log, "D1b", compute)
            for soil, beta, hw, pv, s0, s1 in rows:
                r0, r1 = _get(res, s0), _get(res, s1)
                log(f"   {soil} beta={beta:.0f} hw/H={hw}: H = H_crit,paper = {pv:.3f} m: Gamma(box ~ H) = {r0['Gamma']:.4f}, "
                    f"Gamma(box 50/10/30 m) = {r1['Gamma']:.4f}   (consistent hypothesis gives 1)")
        elif name == "eq31":
            log("D2) Eq. (31) as PRINTED vs corrected: v'_opt zone-1 field with the sign typo (prefactor 1/(C-D)), the")
            log("    exponent/coefficient typo (C/D instead of sqrt(C/D)) or both; Gamma after full optimisation")
            log("=" * 110)
            variants = ("vopt", "vopt_eq31_printed", "vopt_eq31_sign", "vopt_eq31_expo")
            specs9 = {(a, b, v): _fig9_base(args, a, b, v) for a in (1, 10) for b in (30.0, 45.0, 60.0, 75.0, 90.0)
                      for v in variants}
            specs8 = {(soil, beta, hw, v): _fig8_base(args, "fitted", soil, beta, hw, v)
                      for soil, beta in (("London", 30.0), ("London", 60.0), ("Israeli", 60.0)) for hw in (0.5, 1.0)
                      for v in variants}
            allspecs = list(specs9.values()) + list(specs8.values())
            run_specs([q for q in allspecs if q["approach"] != "vopt"], args.workers, cache, log, "D2 (variants)", compute)
            cache.refresh()
            res, _ = run_specs(allspecs, args.workers, cache, log, "D2", compute)
            log("   Fig. 9 (Gamma):")
            log(f"   {'alpha':>5} {'beta':>5} {'paper':>8} | " + " | ".join(f"{v:>18}" for v in variants))
            out = []
            for (a, b) in sorted({(k[0], k[1]) for k in specs9}):
                pv, _ = interp_polyline(p9[(a, "vopt")], b)
                cells = []
                for v in variants:
                    g = _get(res, specs9[(a, b, v)])["Gamma"]
                    cells.append(f"{g:8.4f} ({_rel(g, pv):+7.3f})")
                    out.append(dict(fig="fig9", case=f"alpha={a} beta={b:g}", variant=v, value=g, paper=pv,
                                    rel=_rel(g, pv)))
                log(f"   {a:5d} {b:5.0f} {pv:8.4f} | " + " | ".join(f"{c:>18}" for c in cells))
            log("   Fig. 8 ('fitted' parameters, H_crit in m):")
            for (soil, beta, hw) in sorted({k[:3] for k in specs8}):
                pv, _ = interp_polyline(p8[(soil, beta, "vopt")], hw, log=True)
                cells = []
                for v in variants:
                    g = _get(res, specs8[(soil, beta, hw, v)])["Hcrit"]
                    cells.append(f"{g:8.3f} ({_rel(g, pv):+7.3f})")
                    out.append(dict(fig="fig8", case=f"{soil} beta={beta:g} hw/H={hw}", variant=v, value=g, paper=pv,
                                    rel=_rel(g, pv)))
                log(f"   {soil:>7} {beta:3.0f} hw/H={hw:3.1f} {pv:8.3f} | " + " | ".join(f"{c:>18}" for c in cells))
            _write_csv(os.path.join(OUT_DIR, "diag_eq31_variants.csv"), list(out[0].keys()), out)
        elif name == "gamma_w":
            log("D3) gamma_w = 9.8 instead of 9.81 in Fig. 9 (both the buoyant weight and the seepage field)")
            log("=" * 110)
            specs = {(b, ap, gw): _fig9_base(args, 1, b, ap, gamma_w=gw) for b in (30.0, 60.0, 90.0)
                     for ap in APPROACHES for gw in (9.81, 9.8)}
            res, _ = run_specs(list(specs.values()), args.workers, cache, log, "D3", compute)
            for b in (30.0, 60.0, 90.0):
                for ap in APPROACHES:
                    g1, g2 = _get(res, specs[(b, ap, 9.81)])["Gamma"], _get(res, specs[(b, ap, 9.8)])["Gamma"]
                    pv, _ = interp_polyline(p9[(1, ap)], b)
                    log(f"   alpha=1 beta={b:4.0f} {ap:>4}: gw 9.81 -> {g1:.5f}, gw 9.8 -> {g2:.5f} "
                        f"(change {g2 / g1 - 1:+.5f}); paper {pv:.4f}")
        elif name == "fe_bc":
            log("D4) Fig. 9 FE, alpha = 1: far-side condition / mesh / P_u evaluation (box scaled with H)")
            log("=" * 110)
            variants = (("default (P2 ref1 zero_lb, boundary P_u)", {}, {}),
                        ("impermeable sides", dict(bc="impermeable"), {}),
                        ("P1 ref0 (coarse)", dict(order=1, ref=0), {}),
                        ("domain quadrature P_u", {}, dict(pu_method="domain")))
            specs = {(b, i): _fig9_base(args, 1, b, "FE", field_kw=fk, la_kw=lk) for b in (30.0, 45.0, 60.0)
                     for i, (_, fk, lk) in enumerate(variants)}
            res, _ = run_specs(list(specs.values()), args.workers, cache, log, "D4", compute)
            for b in (30.0, 45.0, 60.0):
                pv, _ = interp_polyline(p9[(1, "FE")], b)
                for i, (lab, _, _) in enumerate(variants):
                    g = _get(res, specs[(b, i)])["Gamma"]
                    log(f"   beta={b:4.0f} {lab:>42}: Gamma = {g:.4f} (rel to paper {_rel(g, pv):+.4f})")
        elif name == "classes":
            log("D5) Mechanism classes: best Gamma of class I (B on the face) and II (B below the toe) for every Fig. 8/9 run")
            log("=" * 110)
            n_ii, n_tot, close = 0, 0, []
            for r in cache.data.values():
                s = r["spec"]
                if s.get("kind") != "la" or not r.get("best_by_class") or s.get("field_kw") or s.get("la_kw") or \
                        s["approach"] not in ("vopt", "FE", "none"):
                    continue
                n_tot += 1
                bb = r["best_by_class"]
                gi, gii = bb.get("I", np.inf), bb.get("II", np.inf)
                if gii < gi * (1 - 1e-9):
                    n_ii += 1
                    log(f"      class II lower: {s.get('soil', '')} alpha={s['alpha']} beta={s['beta']} "
                        f"hw/H={s['hw_over_H']} {s['approach']}: Gamma_I={gi:.4f} Gamma_II={gii:.4f} "
                        f"(d/H={r['params'].get('d_over_H', np.nan):.3f})")
                if np.isfinite(gi) and np.isfinite(gii):
                    close.append((gii / gi - 1.0, s))
            log(f"   {n_tot} runs; class II strictly lower in {n_ii}.  Optimum of class II = toe mechanism (d = 0) when "
                f"Gamma_II = Gamma_I:")
            gaps = np.array([c[0] for c in close])
            if len(gaps):
                log(f"   Gamma_II / Gamma_I - 1: min {gaps.min():+.2e}, median {np.median(gaps):+.2e}, max {gaps.max():+.2e}")
            etas = [r["params"].get("eta", np.nan) for r in cache.data.values() if r["spec"].get("kind") == "la"
                    and r.get("params") and "eta" in r["params"] and not r["spec"].get("field_kw")
                    and not r["spec"].get("la_kw") and r["spec"]["approach"] in ("vopt", "FE", "none")]
            log(f"   class-I optima with eta < 0.999 (B above the toe): "
                f"{sum(1 for e in etas if e < 0.999)} of {len(etas)}")
        elif name == "robust":
            log("D7) Optimiser robustness at the Fig. 8 outliers: 6 seeds, 80 particles, 300 iterations (default: 2 seeds,")
            log("    40 particles, 150 iterations), both mechanism classes, d_max = 30 H")
            log("=" * 110)
            specs = {}
            for soil, beta, hw, ap in ROBUST_POINTS:
                specs[(soil, beta, hw, ap, "default")] = _fig8_base(args, "fitted", soil, beta, hw, ap)
                specs[(soil, beta, hw, ap, "robust")] = _fig8_base(
                    args, "fitted", soil, beta, hw, ap, seeds=[0, 1, 2, 3, 4, 5],
                    la_kw=dict(n_particles=80, n_iter=300, d_max=30.0))
            res, _ = run_specs(list(specs.values()), args.workers, cache, log, "D7", compute)
            for soil, beta, hw, ap in ROBUST_POINTS:
                pv, _ = interp_polyline(p8[(soil, beta, ap)], hw, log=True)
                r0, r1 = _get(res, specs[(soil, beta, hw, ap, "default")]), _get(res, specs[(soil, beta, hw, ap, "robust")])
                log(f"   {soil} beta={beta:.0f} hw/H={hw:4.2f} {ap:>4}: paper {pv:8.3f} | default {r0['Hcrit']:8.3f} "
                    f"(spread {r0.get('seed_spread', np.nan):.1e}) | robust {r1['Hcrit']:8.3f} "
                    f"(spread {r1.get('seed_spread', np.nan):.1e}, {r1['mechanism']} {_s_param(r1):.4f}, "
                    f"theta1={(r1['params'] or {}).get('theta1', np.nan):.4f} theta2={(r1['params'] or {}).get('theta2', np.nan):.4f}"
                    f" L/H={r1.get('L', np.nan) / FIG8_H_REF:.4f})")
        elif name == "beta15":
            log("D6) beta = 15 deg (< phi = 30): admissible mechanisms with P_ext > 0 by a brute-force grid (vopt and FE)")
            log("=" * 110)
            _beta15_grid(args, log)
        else:
            log(f"unknown diagnostic {name!r}; available: {', '.join(DIAGNOSTICS)}")


def _beta15_grid(args, log):
    from limit_analysis import Problem, _grid_min
    for a in FIG9_ALPHAS:
        for ap in APPROACHES:
            t0 = time.time()
            f = make_field(ap, 15.0, FIG9_DATA["H"], FIG9_DATA["H"], a, FIG9_DATA["gamma_w"])
            prob = Problem(15.0, FIG9_DATA["H"], FIG9_DATA["c"], FIG9_DATA["phi"], FIG9_DATA["gamma"],
                           FIG9_DATA["gamma_w"], f)
            best = []
            for kind in ("I", "II"):
                g, x, _ = _grid_min(prob, kind, 40, 40, 16, "coarse")
                best.append((g, kind, x))
            g, kind, x = min(best, key=lambda t: t[0])
            txt = "none with P_ext > 0 (Gamma = inf)" if g >= 1e29 else f"{g:.4g} ({kind}, x = {np.round(x, 4)})"
            log(f"   alpha={a:2d} {ap:>4}: grid (40x40x16 per class) min Gamma = {txt}  [{time.time() - t0:.0f} s]")


def _plot_diag_box(out, p9, path):
    plt = _mpl()
    fig, axs = plt.subplots(1, 3, figsize=(10.5, 3.6), sharey=True)
    for j, a in enumerate(FIG9_ALPHAS):
        ax = axs[j]
        _style_ax(ax)
        P = p9[(a, "FE")]
        ax.plot(P[:, 0], P[:, 1], color=COL["paper"], lw=0.9, label=r"paper $-\nabla u'_{FE}$")
        cs = sorted([c for c in out if c["alpha"] == a], key=lambda c: c["beta_deg"])
        b = [c["beta_deg"] for c in cs]
        ax.plot(b, [c["Gamma_box_scaled"] for c in cs], color=COL["FE"], lw=1.8, ls=(0, (2, 1.5)), marker="o", ms=3,
                label="ours, box 50/10/30 H (scaled)")
        ax.plot(b, [c["Gamma_box_metres"] for c in cs], color=COL["vopt"], lw=1.8, marker="s", ms=3,
                label="ours, box 50/10/30 m (H = 5 m)")
        ax.set_ylim(0, 5)
        ax.set_xlim(15, 90)
        ax.set_xticks([15, 30, 45, 60, 75, 90])
        ax.set_title(rf"$\alpha$ = {a}", loc="left", color=COL["text"])
        ax.set_xlabel(r"$\beta$ (deg)")
    axs[0].set_ylabel(r"$\Gamma$")
    h, l = axs[0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=3, fontsize=8)
    fig.suptitle("Fig. 9 FE curve: FE box scaled with H vs fixed in metres", fontsize=9.5, color=COL["text"])
    fig.tight_layout(rect=(0, 0.1, 1, 0.95))
    fig.savefig(path, dpi=150)
    plt.close(fig)


# ======================================================================================================
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--only", default="fig5,fig8,fig9", help="comma list of fig5, fig8, fig8table1, fig9")
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--seeds", default="0,1", help="PSO seeds per mechanism class")
    ap.add_argument("--quick", action="store_true", help="reduced grids (smoke test)")
    ap.add_argument("--fresh", action="store_true", help="ignore the cache")
    ap.add_argument("--plot-only", action="store_true", help="no computation: tables/plots from the cache")
    ap.add_argument("--no-plot", action="store_true")
    ap.add_argument("--diagnose", default="", help="comma list of diagnostic experiments")
    args = ap.parse_args(argv)
    args.seeds = tuple(int(s) for s in args.seeds.split(",") if s.strip())
    os.makedirs(OUT_DIR, exist_ok=True)
    cache = Cache(os.path.join(OUT_DIR, "cache.jsonl" if not args.quick else "cache_quick.jsonl"), args.fresh)
    logpath = os.path.join(OUT_DIR, "reproduce_output.txt" if not args.diagnose else
                           "diagnostics_" + "_".join(n.strip() for n in args.diagnose.split(",") if n.strip()) + ".txt")
    if args.quick:
        logpath = logpath.replace(".txt", "_quick.txt")
    fh = open(logpath, "w")

    def log(s=""):
        print(s, flush=True)
        fh.write(s + "\n")
        fh.flush()

    t0 = time.time()
    log(f"reproduce_python.py  {time.strftime('%Y-%m-%d %H:%M:%S')}  workers={args.workers} seeds={args.seeds} "
        f"quick={args.quick}")
    walls = {}
    if args.diagnose:
        diagnostics([n.strip() for n in args.diagnose.split(",") if n.strip()], args, cache, log)
    else:
        only = {s.strip() for s in args.only.split(",")}
        if "fig5" in only:
            walls["fig5"] = fig5(args, cache, log)[1]
        if "fig9" in only:
            walls["fig9"] = fig9(args, cache, log)[1]
        if "fig8" in only:
            walls["fig8-fitted"] = fig8(args, cache, log, "fitted")[1]
        if "fig8" in only or "fig8table1" in only:
            walls["fig8-table1"] = fig8(args, cache, log, "table1")[1]
    log("")
    log("wall time of new computations: " + ", ".join(f"{k} {v:.1f} s" for k, v in walls.items()))
    log(f"total wall time {time.time() - t0:.1f} s")
    fh.close()


if __name__ == "__main__":
    main()
