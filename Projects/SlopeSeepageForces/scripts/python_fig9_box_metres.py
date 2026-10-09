#!/usr/bin/env python3
"""Python reference of the Fig. 9 FE curve with the hydraulic box FIXED IN METRES (50 m / 10 m / 30 m left of O /
right of the toe / below the toe, u = 0 on the left side and the base), at every beta of the Fig. 9 grid of
reproduce_python.py.

data/python_fig9.csv holds the FE curve with the box scaled with H (50 H / 10 H / 30 H); reproduce_python.py's
diagnostic D1 computed the box-in-metres variant (field_kw = box_in_metres(5) = 10 H / 2 H / 6 H) at the 5-deg nodes
20..90 only.  This script runs the missing betas (15, 37.5, 52.5, 67.5, 82.5) through the same task function and cache
(results/reproduce_python/cache.jsonl, so the existing results are reused and new ones are appended: resumable) and
writes results/cpp/python_fig9_FE_box_m.csv (columns of data/paper_fig9.csv plus the mechanism).

    python3 Projects/SlopeSeepageForces/scripts/python_fig9_box_metres.py [--workers 1]
"""
import argparse
import csv
import os
import sys
import time

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPT_DIR)

import reproduce_python as rp  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--seeds", default="0,1", help="as reproduce_python.py (default 0,1: reuses its cache)")
    args = ap.parse_args()
    args.seeds = [int(s) for s in args.seeds.split(",")]
    specs = [dict(kind="la", fig="fig9", alpha=a, beta=b, approach="FE", hw_over_H=1.0, seeds=list(args.seeds),
                  field_kw=rp.box_in_metres(rp.FIG9_DATA["H"]), la_kw={}, **rp.FIG9_DATA)
             for a in rp.FIG9_ALPHAS for b in rp.FIG9_BETAS]
    out_dir = os.path.join(rp.PROJECT_DIR, "results", "cpp")
    os.makedirs(out_dir, exist_ok=True)
    log_path = os.path.join(out_dir, "python_fig9_FE_box_m.log")

    def log(s=""):
        print(s, flush=True)
        with open(log_path, "a") as fh:
            fh.write(time.strftime("%Y-%m-%d %H:%M:%S  ", time.gmtime()) + s + "\n")

    cache = rp.Cache(os.path.join(rp.OUT_DIR, "cache.jsonl"))
    res, wall = rp.run_specs(specs, args.workers, cache, log, "fig9-FE-box-m")
    rows = []
    for s in specs:
        r = res.get(rp.spec_key(s))
        if r is None:
            continue
        g = r["Gamma"]
        p = r.get("params") or {}
        rows.append(dict(alpha=s["alpha"], beta_deg=rp._fmt_beta(s["beta"]), curve="FE",
                         Gamma=f"{g:.6f}" if g < 1e29 else "inf", mechanism=r["mechanism"] or "",
                         eta_or_dH=p.get("eta", p.get("d_over_H", "")), seed_spread=r.get("seed_spread", ""),
                         time_s=f"{r['time']:.1f}"))
    path = os.path.join(out_dir, "python_fig9_FE_box_m.csv")
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["alpha", "beta_deg", "curve", "Gamma", "mechanism", "eta_or_dH",
                                           "seed_spread", "time_s"])
        w.writeheader()
        w.writerows(rows)
    log(f"wrote {path}: {len(rows)} of {len(specs)} cases (wall {wall:.1f} s for the new ones)")


if __name__ == "__main__":
    main()
