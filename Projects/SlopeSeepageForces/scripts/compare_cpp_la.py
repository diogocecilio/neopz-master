#!/usr/bin/env python3
"""C++ limit analysis (LimitAnalysis.h, command 'la') against the Python reference (limit_analysis.py).

Cases: limit analyses already computed by reproduce_python.py (results/reproduce_python/cache.jsonl), with
    FE field  (-grad u'_FE, P_u by the boundary formula; C++ hydraulic mesh href = 0 and 1): the FE box of the paper
              fixed in METRES, 50 / 10 / 30 m
              (left of O / right of the toe / below the toe), u = 0 on the left side and the base ("zero_lb"):
              Fig. 9 cases (H = 5 m, field_kw = box 10 H / 2 H / 6 H) run at H = 5 m; Fig. 8 cases (Python at
              H = 10 m with the box 50 H / 10 H / 30 H) run in C++ at H = 1 m with the box 50 / 10 / 30 m, i.e. the
              same problem by similarity: H_crit = Gamma(1 m) * 1 m = Gamma(10 m) * 10 m;
    analytical field (K^-1 v'_opt, P_u by the domain quadrature with the circle splits): same H as the Python.
The C++ FE field (SeepageFE.h, Delaunay mesh, href) and the Python one (fe_seepage.py, ref = 1) are different
discretisations of the same problem, so the FE comparison mixes the optimiser tolerance with the FE discretisation
error; the analytical field is the same function in both (AnalyticalSeepage.h matches analytical_seepage.py to ~1e-8).

    python3 Projects/SlopeSeepageForces/scripts/compare_cpp_la.py cases     # writes the case file (once)
    <build>/Projects/SlopeSeepageForces/SlopeSeepageForces labatch \
        cases=Projects/SlopeSeepageForces/results/limit_analysis_cpp/cases_vs_python.txt \
        out=Projects/SlopeSeepageForces/results/limit_analysis_cpp/cpp_vs_python.csv threads=2   # resumable
    python3 Projects/SlopeSeepageForces/scripts/compare_cpp_la.py compare   # table + comparison csv
"""
from __future__ import annotations

import csv
import json
import os
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.dirname(SCRIPT_DIR)
CACHE = os.path.join(PROJECT_DIR, "results", "reproduce_python", "cache.jsonl")
OUT_DIR = os.path.join(PROJECT_DIR, "results", "limit_analysis_cpp")
CASES = os.path.join(OUT_DIR, "cases_vs_python.txt")
MAP = os.path.join(OUT_DIR, "cases_vs_python_map.json")
CPP_CSV = os.path.join(OUT_DIR, "cpp_vs_python.csv")
CMP_CSV = os.path.join(OUT_DIR, "comparison_cpp_vs_python.csv")

BOX_M_FIG9 = {"left": 10.0, "right": 2.0, "depth": 6.0}       # 50 / 10 / 30 m at H = 5 m

# (fig, pset, soil, beta, alpha, hw_over_H, approach, field_kw) of the selected Python results
SELECT = [
    # Fig. 9, FE field with the box in metres (FE_box_m), H = 5 m
    ("fig9", None, None, 20.0, 1, 1.0, "FE", BOX_M_FIG9), ("fig9", None, None, 30.0, 1, 1.0, "FE", BOX_M_FIG9),
    ("fig9", None, None, 45.0, 1, 1.0, "FE", BOX_M_FIG9), ("fig9", None, None, 60.0, 1, 1.0, "FE", BOX_M_FIG9),
    ("fig9", None, None, 75.0, 1, 1.0, "FE", BOX_M_FIG9), ("fig9", None, None, 90.0, 1, 1.0, "FE", BOX_M_FIG9),
    ("fig9", None, None, 25.0, 5, 1.0, "FE", BOX_M_FIG9), ("fig9", None, None, 40.0, 5, 1.0, "FE", BOX_M_FIG9),
    ("fig9", None, None, 60.0, 5, 1.0, "FE", BOX_M_FIG9), ("fig9", None, None, 90.0, 5, 1.0, "FE", BOX_M_FIG9),
    ("fig9", None, None, 20.0, 10, 1.0, "FE", BOX_M_FIG9), ("fig9", None, None, 30.0, 10, 1.0, "FE", BOX_M_FIG9),
    ("fig9", None, None, 50.0, 10, 1.0, "FE", BOX_M_FIG9), ("fig9", None, None, 90.0, 10, 1.0, "FE", BOX_M_FIG9),
    # Fig. 8 ('fitted' = swapped Table 1 parameters, gamma_w = 9.8), FE field, box 50 H / 10 H / 30 H
    ("fig8", "fitted", "London", 30.0, 1, 0.1, "FE", {}), ("fig8", "fitted", "London", 30.0, 1, 0.5, "FE", {}),
    ("fig8", "fitted", "London", 60.0, 1, 0.2, "FE", {}), ("fig8", "fitted", "London", 60.0, 1, 1.0, "FE", {}),
    ("fig8", "fitted", "Israeli", 35.0, 1, 0.1, "FE", {}), ("fig8", "fitted", "Israeli", 35.0, 1, 0.6, "FE", {}),
    ("fig8", "fitted", "Israeli", 60.0, 1, 0.5, "FE", {}), ("fig8", "fitted", "Israeli", 60.0, 1, 0.9, "FE", {}),
    # analytical field K^-1 v'_opt
    ("fig9", None, None, 15.0, 1, 1.0, "vopt", {}), ("fig9", None, None, 20.0, 1, 1.0, "vopt", {}),
    ("fig9", None, None, 45.0, 1, 1.0, "vopt", {}), ("fig9", None, None, 25.0, 5, 1.0, "vopt", {}),
    ("fig9", None, None, 50.0, 5, 1.0, "vopt", {}), ("fig9", None, None, 90.0, 5, 1.0, "vopt", {}),
    ("fig9", None, None, 35.0, 10, 1.0, "vopt", {}), ("fig9", None, None, 90.0, 10, 1.0, "vopt", {}),
    ("fig8", "fitted", "London", 30.0, 1, 0.05, "vopt", {}), ("fig8", "fitted", "London", 30.0, 1, 0.5, "vopt", {}),
    ("fig8", "fitted", "London", 60.0, 1, 0.1, "vopt", {}), ("fig8", "fitted", "Israeli", 35.0, 1, 0.2, "vopt", {}),
    ("fig8", "fitted", "Israeli", 60.0, 1, 0.7, "vopt", {}),
]


def _load_cache():
    out = []
    with open(CACHE) as fh:
        for line in fh:
            line = line.strip()
            if line:
                try:
                    out.append(json.loads(line))
                except json.JSONDecodeError:
                    pass
    return out


def _find(rows, sel):
    fig, pset, soil, beta, alpha, hw, ap, fkw = sel
    for r in rows:
        s = r["spec"]
        if (s.get("kind") == "la" and s.get("fig") == fig and s.get("pset") == pset and s.get("soil") == soil
                and s["beta"] == beta and s["alpha"] == alpha and s["hw_over_H"] == hw and s["approach"] == ap
                and (s.get("field_kw") or {}) == fkw and not s.get("la_kw")
                and s["H"] == (10.0 if fig == "fig8" else 5.0)):     # (the cache also holds diagnostic runs at other H)
            return r
    return None


def _case_line(spec):
    """options of the C++ 'la' command for a Python spec (see the module doc)"""
    s = spec
    if s["approach"] == "FE":
        H = 1.0 if s["fig"] == "fig8" else s["H"]
        water = "water=fe hboxm=50,10,30"
    else:
        H = s["H"]
        water = "water=analytical"
    return (f"H={H:g} beta={s['beta']:g} hw={s['hw_over_H']:g} alpha={s['alpha']:g} c={s['c']:g} phi={s['phi']:g} "
            f"gamma={s['gamma']:g} gammaw={s['gamma_w']:g} {water} seeds={','.join(str(x) for x in s['seeds'])}")


def write_cases():
    os.makedirs(OUT_DIR, exist_ok=True)
    rows = _load_cache()
    lines, mapping = [], {}
    for sel in SELECT:
        r = _find(rows, sel)
        if r is None:
            print("missing in the Python cache:", sel)
            continue
        base = _case_line(r["spec"])
        # FE cases twice: C++ hydraulic mesh href = 0 (default) and href = 1 (one uniform refinement)
        for line in ([base, base + " href=1"] if r["spec"]["approach"] == "FE" else [base]):
            lines.append(line)
            mapping[line] = dict(spec=r["spec"], Gamma=r["Gamma"], Hcrit=r["Hcrit"], mechanism=r["mechanism"],
                                 x=r.get("x"), P_mr=r.get("P_mr"), P_gamma=r.get("P_gamma"), P_u=r.get("P_u"),
                                 time=r.get("time"))
    with open(CASES, "w") as fh:
        fh.write("# C++ 'la' cases mirroring Python results of reproduce_python.py (see scripts/compare_cpp_la.py)\n")
        fh.write("\n".join(lines) + "\n")
    with open(MAP, "w") as fh:
        json.dump(mapping, fh, indent=1)
    print(f"{len(lines)} cases written to {CASES}")


def compare():
    mapping = json.load(open(MAP))
    cpp = {}
    with open(CPP_CSV) as fh:
        for row in csv.DictReader(fh):
            cpp[row["id"]] = row
    out = []
    print(f"{'case':<96} {'Python':>12} {'C++':>12} {'rel':>9} {'mech py/c++':>11} {'t_py':>6} {'t_c++':>6}")
    for line, ref in mapping.items():
        row = cpp.get(line)
        if row is None:
            print(f"{line:<96} (not run yet)")
            continue
        s = ref["spec"]
        # compare H_crit (= Gamma H): identical to Gamma for equal H, and the similarity-invariant for the Fig. 8 FE
        # cases run at H = 1 m (Python at H = 10 m)
        py, cc = float(ref["Hcrit"]), float(row["Hcrit"])
        rel = cc / py - 1.0 if py not in (0.0, float("inf")) and cc != float("inf") else float("nan")
        mpy = (ref["mechanism"] or "-") + (f"{ref['x'][2]:.3f}" if ref.get("x") else "")
        mc = row["mechanism"] + (f"{float(row['s']):.3f}" if row["s"] else "")
        tpy, tc = float(ref.get("time") or float("nan")), float(row["t_la"]) + float(row["t_field"])
        print(f"{line:<96} {py:12.6f} {cc:12.6f} {rel:+9.2e} {mpy:>5}/{mc:<6} {tpy:6.1f} {tc:6.1f}")
        out.append(dict(case=line, fig=s["fig"], approach=s["approach"], beta=s["beta"], alpha=s["alpha"],
                        hw_over_H=s["hw_over_H"], H_python=s["H"], H_cpp=float(row["H"]), Hcrit_python=py, Hcrit_cpp=cc,
                        Gamma_python=ref["Gamma"], Gamma_cpp=float(row["Gamma"]), rel_diff=rel,
                        mechanism_python=mpy, mechanism_cpp=mc, time_python=tpy, time_cpp=tc))
    if out:
        with open(CMP_CSV, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(out[0].keys()))
            w.writeheader()
            for r in out:
                w.writerow({k: (f"{v:.8g}" if isinstance(v, float) else v) for k, v in r.items()})
        for ap, tag in (("FE", ""), ("FE", "href=1"), ("vopt", "")):
            rel = [abs(r["rel_diff"]) for r in out if r["approach"] == ap and r["rel_diff"] == r["rel_diff"]
                   and (tag in r["case"] if tag else "href=" not in r["case"])]
            if rel:
                print(f"{ap} {tag or '(href=0)' if ap == 'FE' else ''}: {len(rel)} cases, max |rel| {max(rel):.2e}, "
                      f"mean |rel| {sum(rel) / len(rel):.2e}")
        print("written", CMP_CSV)


if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 else "compare"
    if what == "cases":
        write_cases()
    else:
        compare()
