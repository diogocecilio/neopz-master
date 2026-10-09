#!/usr/bin/env python3
"""Figures and comparison tables of the C++ reproduction of Figs. 5, 8 and 9 (Ceron et al., IJNAMG 2025, Sect. 4.3).

Reads the resumable CSVs of the C++ commands (results/cpp/, see scripts/run_cpp_figures.sh):
    fig5.csv                  SlopeSeepageForces fig5 out=results/cpp/fig5.csv
    fig8.csv, fig8_table1.csv SlopeSeepageForces fig8 [soil=table1]
    fig9.csv                  SlopeSeepageForces fig9
    fig8_hw0_gammaw.csv       h_w = 0 ends with gamma_w = 9.8 and 9.81 (both soil sets)
    fig9_gammaw9.8.csv        Fig. 9 with gamma_w = 9.8
and compares them with
    the paper  : data/paper_fig{5,8,9}_vertices.csv (digitized polylines; Fig. 8 interpolated in log scale)
    the Python : results/reproduce_python/comparison_fig{5,8_fitted,8_table1,9}.csv ("ours" columns, the verified
                 Python reproduction) and results/cpp/python_fig9_FE_box_m.csv (Python FE curve of Fig. 9 with the
                 paper's box fixed in metres, scripts/python_fig9_box_metres.py).
Writes
    results/fig5.png, results/fig8.png, results/fig9.png  (paper layout: paper curves thin grey, ours bold;
        Fig. 5 solid -J*(v'_opt) / dashed J(u'_FE); Fig. 8 dashed vopt / dash-dot FE / solid Wu et al. (paper only),
        log H_crit; Fig. 9 solid FE / dashed vopt; Python results as small open markers)
    results/cpp/comparison_fig{5,8,8_table1,9}.csv and results/cpp/comparison_summary.txt

    python3 Projects/SlopeSeepageForces/scripts/plot_results.py
"""
import csv
import math
import os
import sys

import numpy as np

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.dirname(SCRIPT_DIR)
DATA = os.path.join(PROJECT_DIR, "data")
CPP = os.path.join(PROJECT_DIR, "results", "cpp")
PYRES = os.path.join(PROJECT_DIR, "results", "reproduce_python")
OUT = os.path.join(PROJECT_DIR, "results")

# colours: categorical slots 1-2 of the validated reference palette (as reproduce_python.py), ink tokens for text and
# the paper's curves
COL = {"vopt": "#2a78d6", "FE": "#eb6834", "paper": "#8d8c86", "band": "#ececea", "text": "#0b0b0b",
       "text2": "#52514e", "grid": "#e6e5df"}
LS = {"vopt": (0, (5, 2.2)), "FE8": (0, (7, 2, 1.2, 2)), "FE9": "-", "Wu_rp025": "-", "J_FE": (0, (5, 2.2)),
      "minus_Jstar_opt": "-"}
LAB = {"vopt": r"seepage forces $\underline{K}^{-1}\cdot\underline{v}'_{opt}$",
       "FE": r"seepage forces $-\mathrm{grad}\,u'_{FE}$", "Wu_rp025": r"$r_p$ = 0.25 (Wu et al.)",
       "J_FE": r"$J(u'_{FE})$", "minus_Jstar_opt": r"$-J^*(\underline{v}'_{opt})$"}
FIG8_PANELS = (("London", 30.0), ("London", 60.0), ("Israeli", 35.0), ("Israeli", 60.0))


# ======================================================================================================
# input
# ======================================================================================================
def fnum(s):
    try:
        return float(s)
    except (TypeError, ValueError):
        return math.nan


def read_rows(path):
    """rows of a CSV as dicts; rows with a wrong number of columns (cut by a kill) are skipped"""
    if not os.path.exists(path):
        return []
    with open(path, newline="") as fh:
        r = csv.reader(fh)
        header = next(r, None)
        if header is None:
            return []
        return [dict(zip(header, row)) for row in r if len(row) == len(header)]


def cpp_table(path, key_fields):
    """C++ rows by key; when a key appears with several settings the last row wins (reported)"""
    out, settings = {}, {}
    for row in read_rows(path):
        k = tuple(row[f] for f in key_fields)
        if k in out and row.get("settings") != settings[k]:
            print(f"   note: {os.path.basename(path)}: {k} present with several settings, the last row is used")
        out[k], settings[k] = row, row.get("settings")
    return out


def polylines(name, key_fields, x_field, y_field):
    out = {}
    for r in read_rows(os.path.join(DATA, name)):
        out.setdefault(tuple(r[f] if f in ("curve", "soil") else fnum(r[f]) for f in key_fields), []).append(
            (fnum(r[x_field]), fnum(r[y_field]), int(r.get("visible", 1))))
    return {k: np.array(sorted(v)) for k, v in out.items()}


def interp_polyline(P, x, log=False):
    """value of the paper polyline P (rows x, y, visible) at x and its status: 'vertex' / 'segment' (visible),
    'extrapolated' (a clipped end), 'none' (outside); the same rule as reproduce_python.interp_polyline"""
    xs, ys, vis = P[:, 0], P[:, 1], P[:, 2] > 0.5
    if x < xs[0] - 1e-9 or x > xs[-1] + 1e-9:
        return math.nan, "none"
    j = int(np.argmin(np.abs(xs - x)))
    if abs(xs[j] - x) < 1e-9:
        return float(ys[j]), ("vertex" if vis[j] else "extrapolated")
    i = int(np.searchsorted(xs, x)) - 1
    w = (x - xs[i]) / (xs[i + 1] - xs[i])
    v = 10.0 ** ((1 - w) * np.log10(ys[i]) + w * np.log10(ys[i + 1])) if log else (1 - w) * ys[i] + w * ys[i + 1]
    return float(v), ("segment" if vis[i] and vis[i + 1] else "extrapolated")


def rel(a, b):
    return a / b - 1.0 if (math.isfinite(a) and math.isfinite(b) and b != 0.0) else math.nan


def write_csv(path, fields, rows):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.7g}" if isinstance(v, float) else v) for k, v in r.items()})


def stats(rows, field, keep=lambda r: True):
    v = np.array([r[field] for r in rows if keep(r) and math.isfinite(r[field])])
    if len(v) == 0:
        return "n=0"
    i = int(np.argmax(np.abs(v)))
    return f"n={len(v):2d} mean={100 * v.mean():+7.3f}% rms={100 * np.sqrt(np.mean(v ** 2)):6.3f}% max|.|={100 * abs(v[i]):6.3f}% ({100 * v[i]:+.3f}%)"


# ======================================================================================================
# comparisons
# ======================================================================================================
def compare_fig5():
    cpp = cpp_table(os.path.join(CPP, "fig5.csv"), ("alpha", "beta_deg", "curve"))
    py = {(fnum(r["alpha"]), fnum(r["beta_deg"]), r["curve"]): fnum(r["ours"])
          for r in read_rows(os.path.join(PYRES, "comparison_fig5.csv"))}
    paper = polylines("paper_fig5_vertices.csv", ("alpha", "curve"), "beta_deg", "value")
    rows = []
    for (a, b, c), r in cpp.items():
        a, b = fnum(a), fnum(b)
        v = fnum(r["value"])
        pv, st = interp_polyline(paper[(a, c)], b)
        p = py.get((a, b, c), math.nan)
        rows.append(dict(alpha=int(a), beta_deg=b, curve=c, cpp=v, python=p, rel_cpp_python=rel(v, p), paper=pv,
                         paper_status=st, rel_cpp_paper=rel(v, pv), m=r["m"], degenerate=r["degenerate"], neq=r["neq"]))
    rows.sort(key=lambda r: (r["alpha"], r["curve"], r["beta_deg"]))
    write_csv(os.path.join(CPP, "comparison_fig5.csv"), list(rows[0].keys()) if rows else [], rows)
    return rows


def compare_fig8(name, py_name):
    cpp = cpp_table(os.path.join(CPP, name), ("soil", "beta_deg", "curve", "hw_over_H"))
    py = {(r["soil"], fnum(r["beta_deg"]), r["curve"], fnum(r["hw_over_H"])): (fnum(r["ours"]), r["mechanism"],
                                                                            fnum(r["eta_or_dH"]))
          for r in read_rows(os.path.join(PYRES, py_name))}
    paper = polylines("paper_fig8_vertices.csv", ("soil", "beta_deg", "curve"), "hw_over_H", "Hcrit_m")
    rows = []
    for (s, b, c, hw), r in cpp.items():
        b, hw = fnum(b), fnum(hw)
        v = fnum(r["Hcrit_m"])
        pv, st = interp_polyline(paper[(s, b, c)], hw, log=True)
        p, pm, pe = py.get((s, b, c, hw), (math.nan, "", math.nan))
        rows.append(dict(soil=s, beta_deg=b, curve=c, hw_over_H=hw, cpp=v, python=p, rel_cpp_python=rel(v, p),
                         paper=pv, paper_status=st, rel_cpp_paper=rel(v, pv), mechanism=r["mechanism"],
                         eta_or_dH=fnum(r["eta_or_dH"]), python_mechanism=pm, python_eta_or_dH=pe,
                         seed_spread=fnum(r["seed_spread"]), at_bound=r["at_bound"], gamma_w=fnum(r["gamma_w"]),
                         c=fnum(r["c"]), phi_deg=fnum(r["phi_deg"]), t_s=fnum(r["t_field_s"]) + fnum(r["t_la_s"])))
    order = {p: i for i, p in enumerate(FIG8_PANELS)}
    rows.sort(key=lambda r: (order.get((r["soil"], r["beta_deg"]), 9), r["curve"], r["hw_over_H"]))
    out = "comparison_" + name
    write_csv(os.path.join(CPP, out), list(rows[0].keys()) if rows else [], rows)
    return rows


def compare_fig9(name="fig9.csv"):
    cpp = cpp_table(os.path.join(CPP, name), ("alpha", "beta_deg", "curve"))
    py = {}
    for r in read_rows(os.path.join(PYRES, "comparison_fig9.csv")):
        if r["curve"] == "vopt":
            py[(fnum(r["alpha"]), fnum(r["beta_deg"]), "vopt")] = (fnum(r["ours"]), r["mechanism"], fnum(r["eta_or_dH"]))
        elif r["curve"] == "FE":  # box scaled with H: not the paper's setting (kept for information)
            py[(fnum(r["alpha"]), fnum(r["beta_deg"]), "FE_box_H")] = (fnum(r["ours"]), r["mechanism"],
                                                                       fnum(r["eta_or_dH"]))
    for r in read_rows(os.path.join(CPP, "python_fig9_FE_box_m.csv")):
        py[(fnum(r["alpha"]), fnum(r["beta_deg"]), "FE")] = (fnum(r["Gamma"]), r["mechanism"], fnum(r["eta_or_dH"]))
    paper = polylines("paper_fig9_vertices.csv", ("alpha", "curve"), "beta_deg", "Gamma")
    rows = []
    for (a, b, c), r in cpp.items():
        a, b = fnum(a), fnum(b)
        v = fnum(r["Gamma"])
        pv, st = interp_polyline(paper[(a, c)], b)
        if st == "none":
            st = "offplot(>5)"
        p, pm, pe = py.get((a, b, c), (math.nan, "", math.nan))
        ph = py.get((a, b, "FE_box_H"), (math.nan,))[0] if c == "FE" else math.nan
        rows.append(dict(alpha=int(a), beta_deg=b, curve=c, cpp=v, python=p, rel_cpp_python=rel(v, p), paper=pv,
                         paper_status=st, rel_cpp_paper=rel(v, pv), mechanism=r["mechanism"],
                         eta_or_dH=fnum(r["eta_or_dH"]), python_mechanism=pm, python_eta_or_dH=pe,
                         python_FE_box_scaled=ph, seed_spread=fnum(r["seed_spread"]), at_bound=r["at_bound"],
                         gamma_w=fnum(r["gamma_w"]), t_s=fnum(r["t_field_s"]) + fnum(r["t_la_s"])))
    rows.sort(key=lambda r: (r["alpha"], r["curve"], r["beta_deg"]))
    write_csv(os.path.join(CPP, "comparison_" + name), list(rows[0].keys()) if rows else [], rows)
    return rows


# ======================================================================================================
# plots (paper layout)
# ======================================================================================================
def _mpl():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 9, "axes.edgecolor": COL["text2"], "axes.labelcolor": COL["text"],
                         "xtick.color": COL["text2"], "ytick.color": COL["text2"], "figure.facecolor": "white",
                         "axes.facecolor": "white", "legend.frameon": False, "mathtext.fontset": "cm"})
    return plt


def _box_label(ax, text, x=0.04, y=0.95, va="top", ha="left"):
    ax.text(x, y, text, transform=ax.transAxes, ha=ha, va=va, fontsize=10, color=COL["text"],
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec=COL["text2"], lw=0.9))


def _deg_axis(ax):
    ax.set_xlim(15, 90)
    ax.set_xticks([15, 30, 45, 60, 75, 90])
    ax.set_xticklabels([f"{t}°" for t in (15, 30, 45, 60, 75, 90)])
    ax.grid(True, color=COL["grid"], lw=0.7)
    ax.set_axisbelow(True)


def plot_fig5(rows, path):
    plt = _mpl()
    from matplotlib.lines import Line2D
    paper = polylines("paper_fig5_vertices.csv", ("alpha", "curve"), "beta_deg", "value")
    fig, axs = plt.subplots(2, 2, figsize=(7.6, 6.0))
    ylim = {1: 1.0, 2: 0.70, 4: 0.5, 10: 0.30}
    for k, a in enumerate((1, 2, 4, 10)):
        ax = axs[k // 2][k % 2]
        _deg_axis(ax)
        curves = {}
        for c in ("minus_Jstar_opt", "J_FE"):
            cs = sorted([r for r in rows if r["alpha"] == a and r["curve"] == c], key=lambda r: r["beta_deg"])
            curves[c] = (np.array([r["beta_deg"] for r in cs]), np.array([r["cpp"] for r in cs]),
                         np.array([r["python"] for r in cs]))
        x1, y1, _ = curves["minus_Jstar_opt"]
        x2, y2, _ = curves["J_FE"]
        if len(x1) and len(x1) == len(x2) and np.allclose(x1, x2):
            ax.fill_between(x1, y1, y2, color=COL["band"], lw=0, zorder=1)
        for c in ("minus_Jstar_opt", "J_FE"):
            col = COL["vopt"] if c == "minus_Jstar_opt" else COL["FE"]
            x, y, p = curves[c]
            ax.plot(x, y, ls=LS[c], color=col, lw=2.2, zorder=3)
            ax.plot(x, p, ls="none", marker="o", ms=4.5, mfc="none", mec=COL["text2"], mew=0.8, zorder=5)
            P = paper[(float(a), c)]
            ax.plot(P[:, 0], P[:, 1], ls=LS[c], color=COL["paper"], lw=0.9, zorder=4)
            ok = [r for r in rows if r["alpha"] == a and r["curve"] == c and math.isfinite(r["rel_cpp_paper"])]
            mx = max(abs(r["rel_cpp_paper"]) for r in ok) if ok else math.nan
            ax.text(0.97, 0.05 + (0.0 if c == "J_FE" else 0.08), f"{LAB[c]}: max |C++/paper − 1| = {100 * mx:.2f} %",
                    transform=ax.transAxes, ha="right", va="bottom", fontsize=7, color=COL["text2"])
        ax.set_ylim(0, ylim[a])
        ax.set_yticks(np.linspace(0, ylim[a], 6 if a != 2 and a != 10 else 6))
        _box_label(ax, rf"$\alpha$ = {a}")
        if k % 2 == 0:
            ax.set_ylabel("Normalized Functional")
        if k // 2 == 1:
            ax.set_xlabel("Slope Inclination")
    handles = [Line2D([], [], ls=LS["minus_Jstar_opt"], color=COL["vopt"], lw=2.2, label="C++ " + LAB["minus_Jstar_opt"]),
               Line2D([], [], ls=LS["J_FE"], color=COL["FE"], lw=2.2, label="C++ " + LAB["J_FE"]),
               Line2D([], [], ls="-", color=COL["paper"], lw=0.9, label="paper (digitized)"),
               Line2D([], [], ls="none", marker="o", ms=4.5, mfc="none", mec=COL["text2"], label="Python reference")]
    fig.legend(handles=handles, loc="lower center", ncol=4, fontsize=8)
    fig.suptitle(r"Fig. 5: hydraulic functionals / ($k_h H^2 \gamma_w^2$), $h_w = H$ — C++ (bold) vs paper (thin grey)",
                 fontsize=9.5, color=COL["text"])
    fig.tight_layout(rect=(0, 0.05, 1, 0.97))
    fig.savefig(path, dpi=160)
    plt.close(fig)


def plot_fig8(rows, path, title):
    plt = _mpl()
    from matplotlib.lines import Line2D
    paper = polylines("paper_fig8_vertices.csv", ("soil", "beta_deg", "curve"), "hw_over_H", "Hcrit_m")
    fig, axs = plt.subplots(1, 2, figsize=(8.6, 4.6))
    for j, soil in enumerate(("London", "Israeli")):
        ax = axs[j]
        ax.grid(True, which="both", color=COL["grid"], lw=0.6)
        ax.set_axisbelow(True)
        for s, beta in FIG8_PANELS:
            if s != soil:
                continue
            for c in ("Wu_rp025", "vopt", "FE"):
                P = paper[(soil, beta, c)]
                ls = LS["FE8"] if c == "FE" else LS[c]
                ax.plot(P[:, 0], P[:, 1], ls=ls, color=COL["paper"], lw=0.9, zorder=4)
            for c in ("vopt", "FE"):
                cs = sorted([r for r in rows if r["soil"] == soil and r["beta_deg"] == beta and r["curve"] == c],
                            key=lambda r: r["hw_over_H"])
                x = np.array([r["hw_over_H"] for r in cs])
                y = np.array([r["cpp"] for r in cs])
                p = np.array([r["python"] for r in cs])
                fin = np.isfinite(y)
                ax.plot(x[fin], y[fin], ls=LS["FE8"] if c == "FE" else LS["vopt"], color=COL[c], lw=2.2, zorder=3)
                ax.plot(x, p, ls="none", marker="o", ms=4, mfc="none", mec=COL["text2"], mew=0.8, zorder=5)
            cs = [r for r in rows if r["soil"] == soil and r["beta_deg"] == beta and r["curve"] == "vopt"
                  and abs(r["hw_over_H"] - 1.0) < 1e-9]
            if cs and math.isfinite(cs[0]["cpp"]):
                ax.annotate(rf"$\beta$ = {beta:.0f}°", xy=(1.0, cs[0]["cpp"]), xytext=(4, 0), textcoords="offset points",
                            fontsize=9, color=COL["text"], va="center", annotation_clip=False)
        ax.set_yscale("log")
        ax.set_ylim(3.0, 220.0)
        ax.set_xlim(0, 1)
        ax.set_xlabel(r"$h_w/H$")
        if j == 0:
            ax.set_ylabel(r"$H_{crit}$ (m)")
        _box_label(ax, f"{soil} Clay", x=0.96, ha="right")
    handles = [Line2D([], [], ls=LS["vopt"], color=COL["vopt"], lw=2.2, label="C++ " + LAB["vopt"]),
               Line2D([], [], ls=LS["FE8"], color=COL["FE"], lw=2.2, label="C++ " + LAB["FE"]),
               Line2D([], [], ls="-", color=COL["paper"], lw=0.9, label="paper: " + LAB["Wu_rp025"]),
               Line2D([], [], ls=LS["vopt"], color=COL["paper"], lw=0.9, label="paper: vopt (dashed), FE (dash-dot)"),
               Line2D([], [], ls="none", marker="o", ms=4, mfc="none", mec=COL["text2"], label="Python reference")]
    fig.legend(handles=handles, loc="lower center", ncol=3, fontsize=8)
    fig.suptitle(title, fontsize=9, color=COL["text"])
    fig.tight_layout(rect=(0, 0.11, 0.97, 0.95))
    fig.savefig(path, dpi=160)
    plt.close(fig)


def plot_fig9(rows, path):
    plt = _mpl()
    from matplotlib.lines import Line2D
    paper = polylines("paper_fig9_vertices.csv", ("alpha", "curve"), "beta_deg", "Gamma")
    fig, axs = plt.subplots(1, 3, figsize=(9.4, 3.9), sharey=True)
    for j, a in enumerate((1, 5, 10)):
        ax = axs[j]
        _deg_axis(ax)
        for c in ("FE", "vopt"):
            P = paper[(float(a), c)]
            ax.plot(P[:, 0], P[:, 1], ls=LS["FE9"] if c == "FE" else LS["vopt"], color=COL["paper"], lw=0.9, zorder=4)
            cs = sorted([r for r in rows if r["alpha"] == a and r["curve"] == c], key=lambda r: r["beta_deg"])
            x = np.array([r["beta_deg"] for r in cs])
            y = np.array([r["cpp"] for r in cs])
            p = np.array([r["python"] for r in cs])
            fin = np.isfinite(y)
            ax.plot(x[fin], y[fin], ls=LS["FE9"] if c == "FE" else LS["vopt"], color=COL[c], lw=2.2, zorder=3)
            ax.plot(x, p, ls="none", marker="o", ms=4, mfc="none", mec=COL["text2"], mew=0.8, zorder=5)
        ax.set_ylim(0, 5)
        _box_label(ax, rf"$\alpha$ = {a}", x=0.07, y=0.06, va="bottom")
        if j == 0:
            ax.set_ylabel(r"Stability Factor $\Gamma$")
        if j == 1:
            ax.set_xlabel("Slope inclination")
    handles = [Line2D([], [], ls="-", color=COL["FE"], lw=2.2, label="C++ " + LAB["FE"]),
               Line2D([], [], ls=LS["vopt"], color=COL["vopt"], lw=2.2, label="C++ " + LAB["vopt"]),
               Line2D([], [], ls="-", color=COL["paper"], lw=0.9, label="paper (digitized; same line styles)"),
               Line2D([], [], ls="none", marker="o", ms=4, mfc="none", mec=COL["text2"], label="Python reference")]
    fig.legend(handles=handles, loc="lower center", ncol=4, fontsize=8)
    gw = sorted({r["gamma_w"] for r in rows})
    fig.suptitle(r"Fig. 9: H = 5 m, c = 10 kPa, $\varphi$ = 30°, $\gamma$ = 20 kN/m³, $h_w = H$, $\gamma_w$ = "
                 + "/".join(f"{g:g}" for g in gw) + "; FE box 50/10/30 m — C++ (bold) vs paper (thin grey)",
                 fontsize=9, color=COL["text"])
    fig.tight_layout(rect=(0, 0.08, 1, 0.94))
    fig.savefig(path, dpi=160)
    plt.close(fig)


# ======================================================================================================
# summary
# ======================================================================================================
def main():
    lines = []

    def log(s=""):
        print(s)
        lines.append(s)

    vis = lambda r: r["paper_status"] in ("vertex", "segment")   # noqa: E731  points with a visible paper value
    r5 = compare_fig5()
    if r5:
        plot_fig5(r5, os.path.join(OUT, "fig5.png"))
        log("Fig. 5 (h_w = H): C++ vs paper (digitized polylines) and vs Python, per alpha and curve")
        for a in (1, 2, 4, 10):
            for c in ("minus_Jstar_opt", "J_FE"):
                g = [r for r in r5 if r["alpha"] == a and r["curve"] == c]
                log(f"   alpha {a:2d} {c:16s} paper: {stats(g, 'rel_cpp_paper', vis)}")
                log(f"   {'':8s}{'':16s} Python: {stats(g, 'rel_cpp_python')}")
    for name, py_name, title in (
            ("fig8.csv", "comparison_fig8_fitted.csv",
             r"Fig. 8 ($\alpha$ = 1, $\gamma$ = 18, $\gamma_w$ = 9.8; (c, $\varphi$) of Table 1 swapped: London 11.7 kPa / 24.7°,"
             " Israeli 6 kPa / 32°) — C++ (bold) vs paper (thin grey)"),
            ("fig8_table1.csv", "comparison_fig8_table1.csv",
             r"Fig. 8 with Table 1 as printed (London 6 kPa / 32°, Israeli 11.7 kPa / 24.7°, $\gamma_w$ = 9.8) — C++ (bold)"
             " vs paper (thin grey)")):
        r8 = compare_fig8(name, py_name)
        if not r8:
            continue
        plot_fig8(r8, os.path.join(OUT, "fig8.png" if name == "fig8.csv" else "fig8_table1.png"), title)
        log("")
        log(f"Fig. 8 ({name}): H_crit, C++ vs paper (visible points, log-interpolated) and vs Python")
        for s, b in FIG8_PANELS:
            for c in ("vopt", "FE"):
                g = [r for r in r8 if r["soil"] == s and r["beta_deg"] == b and r["curve"] == c]
                log(f"   {s:7s} {b:2.0f} {c:4s} paper: {stats(g, 'rel_cpp_paper', vis)}")
                log(f"   {'':16s}Python: {stats(g, 'rel_cpp_python')}")
    r9 = compare_fig9()
    if r9:
        plot_fig9(r9, os.path.join(OUT, "fig9.png"))
        log("")
        log("Fig. 9: Gamma, C++ vs paper (visible nodes) and vs Python (FE: box fixed in metres)")
        for a in (1, 5, 10):
            for c in ("FE", "vopt"):
                g = [r for r in r9 if r["alpha"] == a and r["curve"] == c]
                log(f"   alpha {a:2d} {c:4s} paper: {stats(g, 'rel_cpp_paper', vis)}")
                log(f"   {'':13s}Python: {stats(g, 'rel_cpp_python')}")
    # gamma_w evidence: the h_w = 0 ends
    g0 = read_rows(os.path.join(CPP, "fig8_hw0_gammaw.csv"))
    if g0:
        paper = polylines("paper_fig8_vertices.csv", ("soil", "beta_deg", "curve"), "hw_over_H", "Hcrit_m")
        log("")
        log("Fig. 8, h_w = 0 ends (f = 0, gamma' = 18 - gamma_w): H_crit (m), C++ / paper - 1")
        log(f"   {'panel':12s} {'soil set':8s} {'paper':>9s} | {'gamma_w = 9.8':>22s} | {'gamma_w = 9.81':>22s}")
        tab = {}
        for r in g0:
            if r["curve"] == "vopt":
                tab[(r["soil"], fnum(r["beta_deg"]), r["settings"].split("soil=")[-1], fnum(r["gamma_w"]))] = fnum(r["Hcrit_m"])
        for s, b in FIG8_PANELS:
            pv, st = interp_polyline(paper[(s, b, "vopt")], 0.0, log=True)
            for pset in ("swapped", "table1"):
                cells = []
                for gw in (9.8, 9.81):
                    v = tab.get((s, b, pset, gw), math.nan)
                    cells.append(f"{v:10.4f} ({100 * rel(v, pv):+7.3f}%)" if math.isfinite(v) else f"{'inf':>22s}")
                log(f"   {s:7s} {b:2.0f}   {pset:8s} {pv:9.3f} | {cells[0]:>22s} | {cells[1]:>22s}"
                    + ("" if st == "vertex" else f"  [paper {st}]"))
    g9 = compare_fig9("fig9_gammaw9.8.csv") if os.path.exists(os.path.join(CPP, "fig9_gammaw9.8.csv")) else []
    if g9 and r9:
        log("")
        log("Fig. 9 with gamma_w = 9.8 instead of 9.81: C++ vs paper (visible nodes); Gamma(9.8) / Gamma(9.81) - 1")
        base = {(r["alpha"], r["beta_deg"], r["curve"]): r["cpp"] for r in r9}
        for a in (1, 5, 10):
            for c in ("FE", "vopt"):
                g = [r for r in g9 if r["alpha"] == a and r["curve"] == c]
                d = [rel(r["cpp"], base.get((a, r["beta_deg"], c), math.nan)) for r in g]
                d = np.array([x for x in d if math.isfinite(x)])
                log(f"   alpha {a:2d} {c:4s} paper: {stats(g, 'rel_cpp_paper', vis)}; change "
                    f"{100 * d.min():+.3f}..{100 * d.max():+.3f}%" if len(d) else "")
    with open(os.path.join(CPP, "comparison_summary.txt"), "w") as fh:
        fh.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    sys.exit(main())
