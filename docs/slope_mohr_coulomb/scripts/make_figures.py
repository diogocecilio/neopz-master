"""Figures of the report (docs/slope_mohr_coulomb/relatorio.tex) from the outputs of the programs.

Usage: python make_figures.py <slope_dir> <drawdown_dir> <taylor_csv> <laplace_dir> <out_dir>
  slope_dir     run of SlopeMohrCoulomb (SLOPE_VERBOSE=1 ... nref=4 > out.txt) with its VTK files
  drawdown_dir  run of SlopeDrawdown (nref=3 > out.txt) with its VTK files and bishop_*.txt
  taylor_csv    Taylor remainders (model, case, m_type, alpha, err2, err1)
  laplace_dir   independent P1 seepage solution (vertices_steady.csv, vertices_crest.csv)
"""
import csv
import glob
import os
import re
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np

plt.rcParams.update({"font.size": 9, "axes.titlesize": 9, "legend.fontsize": 8, "figure.dpi": 150,
                     "savefig.bbox": "tight", "font.family": "serif"})
slope_dir, dd_dir, taylor_csv, laplace_dir, out = sys.argv[1:6]
os.makedirs(out, exist_ok=True)
BISHOP_DRY = (1.206, 1.806)  # Bishop simplified: SRM, gravity increase (scripts/bishop_pw.py)
GROUND = np.array([[0, 40], [30, 40], [40, 30], [70, 30], [70, 0], [0, 0], [0, 40]])


def save(fig, name):
    fig.savefig(os.path.join(out, name))
    plt.close(fig)


def read_vtk(fname):
    """NeoPZ graph-mesh VTK: points, triangles (quads split, degenerate removed), point scalars and vectors"""
    lines = open(fname).read().split("\n")
    def block(start, count):
        vals, j = [], start
        while len(vals) < count:
            vals += [float(v) for v in lines[j].split()]
            j += 1
        return np.array(vals)
    i = next(k for k, l in enumerate(lines) if l.startswith("POINTS"))
    n = int(lines[i].split()[1])
    pts = block(i + 1, 3 * n).reshape(n, 3)[:, :2]
    i = next(k for k, l in enumerate(lines) if l.startswith("CELLS"))
    tris = []
    for k in range(int(lines[i].split()[1])):
        c = [int(v) for v in lines[i + 1 + k].split()][1:]
        for t in ([c] if len(c) == 3 else [[c[0], c[1], c[2]], [c[0], c[2], c[3]]]):
            a, b, d = pts[t]
            if abs((b[0] - a[0]) * (d[1] - a[1]) - (b[1] - a[1]) * (d[0] - a[0])) > 1e-12:
                tris.append(t)
    fields = {}
    for k, l in enumerate(lines):
        if l.startswith("SCALARS"):
            fields[l.split()[1]] = block(k + 2, n)
        elif l.startswith("VECTORS"):
            fields[l.split()[1]] = block(k + 1, 3 * n).reshape(n, 3)
    # the graph mesh repeats the nodes of each element: merge coincident points (averaging the fields)
    key = np.round(pts, 6)
    uniq, inv = np.unique(key, axis=0, return_inverse=True)
    inv = inv.ravel()
    cnt = np.bincount(inv, minlength=len(uniq)).astype(float)
    for name, v in fields.items():
        if v.ndim == 1:
            fields[name] = np.bincount(inv, weights=v, minlength=len(uniq)) / cnt
        else:
            fields[name] = np.stack([np.bincount(inv, weights=v[:, k], minlength=len(uniq)) / cnt
                                     for k in range(v.shape[1])], axis=1)
    return uniq, inv[np.array(tris)], fields


def outline(ax, zoom=False):
    ax.plot(GROUND[:, 0], GROUND[:, 1], "k-", lw=0.8)
    ax.set_aspect("equal")
    ax.set_xlim(*((12, 62) if zoom else (-1, 71)))
    ax.set_ylim(*((16, 41) if zoom else (-1, 41)))
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")


def circle(ax, xc, yc, R, **kw):
    th = np.linspace(0, 2 * np.pi, 400)
    x, y = xc + R * np.cos(th), yc + R * np.sin(th)
    ground = np.where(x <= 30, 40.0, np.where(x >= 40, 30.0, 70.0 - x))
    m = (y <= ground + 1e-9) & (x >= 0) & (x <= 70)
    x, y = np.where(m, x, np.nan), np.where(m, y, np.nan)
    ax.plot(x, y, **kw)


def bishop(fname):
    """(FS SRM, (xc, yc, R), FS GI) from the output of bishop_pw.py"""
    s = open(fname).read()
    m = re.search(r"FS\(SRM\) = ([\d.]+) circle \(xc, yc, R\) = \(([\d.]+), ([\d.]+), ([\d.]+)\); FS\(GI\) = ([\d.]+)", s)
    return float(m.group(1)), tuple(float(m.group(k)) for k in (2, 3, 4)), float(m.group(5))


# ---------------------------------------------------------------- mesh and boundary ids (TriGMesh(1) of SlopeModel.h)
co = [(x, y) for y in (0, 10, 20, 30) for x in range(0, 80, 10)] + [(x, 40) for x in (0, 10, 20, 30)]
tri = [(0, 1, 8), (1, 9, 8), (1, 2, 9), (2, 10, 9), (2, 3, 10), (3, 11, 10), (3, 4, 11), (4, 12, 11), (4, 5, 12),
       (5, 13, 12), (5, 6, 13), (6, 14, 13), (6, 7, 14), (7, 15, 14), (8, 9, 16), (9, 17, 16), (9, 10, 17), (10, 18, 17),
       (10, 11, 18), (11, 19, 18), (11, 12, 19), (12, 20, 19), (12, 13, 20), (13, 21, 20), (13, 14, 21), (14, 22, 21),
       (14, 15, 22), (15, 23, 22), (16, 17, 24), (17, 25, 24), (17, 18, 25), (18, 26, 25), (18, 19, 26), (19, 27, 26),
       (19, 20, 27), (20, 28, 27), (20, 21, 28), (21, 29, 28), (21, 22, 29), (22, 30, 29), (22, 23, 30), (23, 31, 30),
       (24, 25, 32), (25, 33, 32), (25, 26, 33), (26, 34, 33), (26, 27, 34), (27, 35, 34), (27, 28, 35)]
co = np.array(co, float)
X, T = co, np.array(tri)
for _ in range(1):  # uniform refinement: TriGMesh(1)
    nodes, mid, new = list(map(tuple, X)), {}, []
    def m(a, b):
        key = (min(a, b), max(a, b))
        if key not in mid:
            mid[key] = len(nodes)
            nodes.append(tuple((np.array(nodes[a]) + np.array(nodes[b])) / 2))
        return mid[key]
    for a, b, c in T:
        ab, bc, ca = m(a, b), m(b, c), m(c, a)
        new += [(a, ab, ca), (ab, b, bc), (ca, bc, c), (ab, bc, ca)]
    X, T = np.array(nodes), np.array(new)
fig, ax = plt.subplots(figsize=(6.2, 3.6))
ax.triplot(X[:, 0], X[:, 1], T, lw=0.4, color="0.55")
outline(ax)
box = dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.9)
for txt, x, y, kw in [("$-1$ base (fixa)", 35, 1.2, dict(va="bottom", ha="center")),
                      ("$-2$ rolete", 68.8, 15, dict(rotation=90, va="center", ha="right")),
                      ("$-3$ pé", 56, 28.8, dict(ha="center", va="top")),
                      ("$-4$ topo", 15, 38.8, dict(ha="center", va="top")),
                      ("$-5$ rolete", 1.2, 20, dict(rotation=90, va="center", ha="left")),
                      ("$-6$ face", 33.6, 34.2, dict(rotation=-45, ha="center", va="center"))]:
    ax.text(x, y, txt, fontsize=8, bbox=box, **kw)
ax.set_title("TriGMesh(1): 196 triângulos (P2), talude de 10 m a 45°")
save(fig, "fig_mesh.pdf")

# ---------------------------------------------------------------- Taylor test
if os.path.exists(taylor_csv):
    rows = list(csv.DictReader(open(taylor_csv)))
    labels = {}
    info = os.path.join(os.path.dirname(taylor_csv), "case_info.csv")
    if os.path.exists(info):
        labels = {(r["model"], r["case"]): r["label"] for r in csv.DictReader(open(info))}
    floor = os.path.join(os.path.dirname(taylor_csv), "diag_floor.csv")
    panels = [("PV", "PV (antigo, corrigido)", "err2"), ("Voigt", "RHW (artigo)", "err2")]
    if os.path.exists(floor):
        panels.append(("Jacobi", "RHW, autovalores por Jacobi", "err2_jacobi"))
    fig, axs = plt.subplots(1, len(panels), figsize=(6.9, 2.7), sharey=True)
    for ax, (model, title, col) in zip(axs, panels):
        if model == "Jacobi":
            rr_all = [dict(r, model="Voigt") for r in csv.DictReader(open(floor))]
        else:
            rr_all = [r for r in rows if r["model"] == model]
        for case in dict.fromkeys(r["case"] for r in rr_all):
            rr = [r for r in rr_all if r["case"] == case]
            seen, uniq = set(), []
            for r in rr:
                if r["alpha"] not in seen:
                    seen.add(r["alpha"])
                    uniq.append(r)
            a_ = np.array([float(r["alpha"]) for r in uniq])
            e_ = np.array([float(r[col]) for r in uniq])
            lab = {"1": "geral 1", "2": r"aresta $\varepsilon_2=\varepsilon_3$", "3": r"aresta $\varepsilon_1=\varepsilon_2$",
                   "4": "ápice", "5": "geral 2", "6": "elástico"}.get(case, case)
            ax.loglog(a_, np.maximum(e_, 1e-17), "o-", ms=2.2, lw=0.8, label=lab)
        a_ = np.array([1e-6, 1e-2])
        ax.loglog(a_, 1e2 * a_ ** 2, "k--", lw=0.7, label="ordem 2")
        ax.set_xlabel(r"$\alpha$")
        ax.set_title(title, fontsize=8)
        ax.grid(alpha=0.3, which="both")
    axs[0].set_ylabel(r"$\|\sigma(\varepsilon+\alpha\Delta\varepsilon)-\sigma(\varepsilon)-\alpha\,\mathbb{D}\Delta\varepsilon\|$")
    axs[-1].legend(loc="lower right", fontsize=5.5)
    save(fig, "fig_taylor.pdf")

# ---------------------------------------------------------------- Newton convergence (SLOPE_VERBOSE)
hist, cur = [], []
for l in open(os.path.join(slope_dir, "out.txt")):
    m = re.match(r"\s+it (\d+) \|R\|/\|F\| (\S+)", l)
    if m:
        cur.append(float(m.group(2)))
    elif l.startswith("[") and cur:
        hist.append((l.strip(), cur))
        cur = []
conv = [(h, c) for h, c in hist if "converged" in h and len(c) >= 3]
fig, ax = plt.subplots(figsize=(4.0, 2.9))
for h, c in conv[:: max(1, len(conv) // 8)][:8]:
    m_ = re.match(r"\[(\S+)\] p = (\S+)", h)
    lab_ = "%s, p = %s" % (m_.group(1), m_.group(2).replace(".", ",")[:6]) if m_ else h[:20]
    ax.semilogy(range(1, len(c) + 1), c, "o-", ms=2.5, lw=0.8, label=lab_)
ax.axhline(1e-8, color="k", lw=0.6, ls=":")
ax.set_xlabel("iteração de Newton")
ax.set_ylabel(r"$\|R\|/\|F_{ext}\|$")
ax.grid(alpha=0.3)
ax.legend(fontsize=5.5)
save(fig, "fig_newton.pdf")

# ---------------------------------------------------------------- FS x refinement (dry slope)
tab = []
for l in open(os.path.join(slope_dir, "out.txt")):
    m = re.match(r"\s+(\d+)\s+(\d+)\s+([\d.]+)\s+([\d.]+)\s*$", l)
    if m:
        tab.append([float(v) for v in m.groups()])
tab = np.array(tab)
pv = np.array([[870, 3.043, 1.401], [918, 2.516, 1.312], [1140, 2.094, 1.258], [1814, 1.918, 1.229],
               [3408, 1.840, 1.212], [7142, 1.797, 1.203]])  # SlopeMohrCoulomb pv (README, nref=5)
rhw = [list(r[1:]) for r in tab]
if len(rhw) == 5:
    rhw.append([7684, 1.793, 1.203])  # cycle 5 of the README (nref=5)
rhw = np.array(rhw)
fig, axs = plt.subplots(1, 2, figsize=(6.6, 2.7))
for ax, col, lab, b in [(axs[0], 1, "aumento de gravidade", BISHOP_DRY[1]), (axs[1], 2, "redução de resistência", BISHOP_DRY[0])]:
    ax.plot(range(len(rhw)), rhw[:, col], "o-", ms=3, label="RHW (artigo)")
    ax.plot(range(len(pv)), pv[:, col], "s--", ms=3, mfc="none", label="PV (antigo, corrigido)")
    ax.axhline(b, color="k", lw=0.7, ls=":", label="Bishop simplificado")
    ax.set_xticks(range(len(rhw)))
    ax.set_xticklabels(["%d\n%d" % (k, n) for k, n in enumerate(rhw[:, 0])], fontsize=7)
    ax.set_xlabel("ciclo / equações")
    ax.set_title(lab)
    ax.grid(alpha=0.3)
axs[0].set_ylabel("FS")
axs[1].legend()
save(fig, "fig_fs_refinement.pdf")

# ---------------------------------------------------------------- mechanisms of the dry slope
nref = int(tab[-1, 0])
fig, axs = plt.subplots(1, 2, figsize=(6.8, 2.4))
bdry = os.path.join(dd_dir, "bms", "bishop_dry.txt")
circ_dry = bishop(bdry)[1] if os.path.exists(bdry) else None
for ax, kind, circ in [(axs[0], "GI", None), (axs[1], "SRM", circ_dry)]:
    pts, tris, f = read_vtk(os.path.join(slope_dir, "slope_rhw_%s_ref%d.scal_vec.0.vtk" % (kind, nref)))
    v = np.maximum(f["StrainPlasticJ2"], 0)
    tp = ax.tripcolor(pts[:, 0], pts[:, 1], tris, v / v.max(), shading="gouraud", cmap="inferno_r", vmin=0, vmax=0.2)
    ax.triplot(pts[:, 0], pts[:, 1], tris, lw=0.1, color="0.75")
    outline(ax, zoom=True)
    if circ:
        circle(ax, *circ, color="c", lw=1.0, ls="--")
    ax.set_title(("%s: " % kind) + r"$\sqrt{J_2(\varepsilon^p)}/\max$ no colapso")
fig.colorbar(tp, ax=axs, shrink=0.8)
save(fig, "fig_mechanism_dry.pdf")

# ---------------------------------------------------------------- pore pressure of the u-p analysis
states = ["crest", "T0", "T0.1", "T1", "T10", "steady"]
titles = {"crest": "reservatório no topo", "T0": "logo após o rebaixamento (T = 0)", "T0.1": "T = 0,1", "T1": "T = 1",
          "T10": "T = 10", "steady": "regime permanente"}
up = {s: read_vtk(os.path.join(dd_dir, "drawdown_up.scal_vec.%d.vtk" % k)) for k, s in enumerate(states)}
fig, axs = plt.subplots(2, 2, figsize=(6.8, 4.6), layout="constrained")
lev = np.arange(0, 401, 25)
for ax, s in zip(axs.flat, ["crest", "T0", "T1", "steady"]):
    pts, tris, f = up[s]
    tr = mtri.Triangulation(pts[:, 0], pts[:, 1], tris)
    cs = ax.tricontourf(tr, f["PorePressure"], levels=lev, cmap="Blues", extend="both")
    ax.tricontour(tr, f["PorePressure"], levels=lev[::2], colors="k", linewidths=0.3)
    outline(ax)
    ax.set_title(titles[s])
fig.colorbar(cs, ax=axs, shrink=0.9, label="p (kPa)")
save(fig, "fig_pressure_fields.pdf")


def profile(pts, tris, val, x0, ys):
    tr = mtri.Triangulation(pts[:, 0], pts[:, 1], tris)
    it = mtri.LinearTriInterpolator(tr, val)
    return np.array(it(np.full_like(ys, x0), ys))


fig, axs = plt.subplots(1, 3, figsize=(6.8, 2.8), sharey=True)
for ax, x0, top in zip(axs, (15.0, 35.0, 55.0), (40.0, 35.0, 30.0)):
    ys = np.linspace(0.01, top - 0.01, 120)
    for s in states:
        pts, tris, f = up[s]
        ax.plot(profile(pts, tris, f["PorePressure"], x0, ys), ys, lw=0.9, label=titles[s])
    ax.plot(10 * (40 - ys), ys, "k:", lw=0.7, label=r"$\gamma_w(40-y)$")
    ax.plot(np.maximum(10 * (30 - ys), 0), ys, "k--", lw=0.7, label=r"$\gamma_w(30-y)^+$")
    ax.set_title("x = %g m" % x0)
    ax.set_xlabel("p (kPa)")
    ax.grid(alpha=0.3)
axs[0].set_ylabel("y (m)")
h, l = axs[0].get_legend_handles_labels()
fig.legend(h, l, loc="upper center", ncol=4, fontsize=6.5, bbox_to_anchor=(0.5, -0.02))
save(fig, "fig_pressure_profiles.pdf")

# ---------------------------------------------------------------- FS of the drawdown states
fs = {}
for l in open(os.path.join(dd_dir, "out.txt")):
    m = re.match(r"\[(\S+)\] refinement (\d+): (\d+) equations, FS GI (\S+), FS SRM (\S+)", l)
    if m:
        fs.setdefault(m.group(1), []).append([int(m.group(2)), int(m.group(3)), float(m.group(4)), float(m.group(5))])
fs = {k: np.array(v) for k, v in fs.items()}
bish = {}
bdir = os.path.join(dd_dir, "bms") if os.path.isdir(os.path.join(dd_dir, "bms")) else dd_dir
for s, fname in [("dry", "bishop_dry.txt"), ("crest", "bishop_crest.txt"), ("drained", "bishop_drained.txt")] + \
        [(st, "bishop_%d.txt" % k) for k, st in enumerate(states) if k > 0]:
    p = os.path.join(bdir, fname)
    if os.path.exists(p) and "FS(SRM)" in open(p).read():
        bish[s] = bishop(p)
# LaTeX table of the FS of every state (FE cycles and Bishop)
num = lambda v: ("%.3f" % v).replace(".", ",")
lab = {"dry": "seco", "crest": "reservatório no topo", "T0": "logo após o rebaixamento", "T0.1": "adensamento",
       "T1": "adensamento", "T10": "adensamento", "steady": "regime permanente", "drained": "drenado (freática no pé)"}
Tlab = {"T0": "0", "T0.1": "0,1", "T1": "1", "T10": "10", "steady": r"$\infty$"}
with open(os.path.join(out, "tab_fs_drawdown.tex"), "w") as ftab:
    for st in ["dry", "crest", "T0", "T0.1", "T1", "T10", "steady", "drained"]:
        if st not in fs:
            continue
        b = bish.get(st)
        ftab.write("%s & %s & %s & %s & %s & %s \\\\\n" % (lab[st], Tlab.get(st, "---"),
                   " / ".join(num(v) for v in fs[st][:, 2]), " / ".join(num(v) for v in fs[st][:, 3]),
                   num(b[2]) if b else "", num(b[0]) if b else ""))
Tval = {"T0": 1e-2, "T0.1": 0.1, "T1": 1.0, "T10": 10.0, "steady": 1e3}
order = [s for s in ["T0", "T0.1", "T1", "T10", "steady"] if s in fs]
fig, axs = plt.subplots(1, 2, figsize=(6.8, 2.8))
for ax, col, bcol, lab in [(axs[0], 3, 0, "redução de resistência"), (axs[1], 2, 2, "aumento de gravidade")]:
    ncyc = min(fs[s].shape[0] for s in order)
    for cyc in range(ncyc):
        ax.semilogx([Tval[s] for s in order], [fs[s][cyc, col] for s in order], "o-", ms=3, lw=0.8,
                    alpha=0.35 + 0.65 * cyc / max(1, ncyc - 1), color="C0", label="EF, ciclo %d" % cyc)
    bs = [s for s in order if s in bish]
    ax.semilogx([Tval[s] for s in bs], [bish[s][bcol] for s in bs], "k^", ms=5, mfc="none", label="Bishop")
    for s, ls, c, name in [("dry", ":", "C3", "seco"), ("crest", "--", "C2", "reservatório no topo"),
                           ("drained", "-.", "C1", "drenado (freática no pé)")]:
        if s in fs:
            ax.axhline(fs[s][-1, col], color=c, lw=0.8, ls=ls, label=name + " (EF)")
    ax.axhline(1.0, color="k", lw=0.5)
    ax.set_xticks([1e-2, 1e-1, 1, 10, 1e3])
    ax.set_xticklabels(["0", "0,1", "1", "10", r"$\infty$"])
    ax.set_xlabel(r"$T = c_v t/H^2$ após o rebaixamento")
    ax.set_title(lab)
    ax.grid(alpha=0.3)
axs[0].set_ylabel("FS")
axs[0].legend(fontsize=5.5)
save(fig, "fig_fs_time.pdf")

# ---------------------------------------------------------------- mechanisms after the drawdown
fig, axs = plt.subplots(1, 2, figsize=(6.8, 2.4))
for ax, s in zip(axs, ["T0", "steady"]):
    pts, tris, f = read_vtk(os.path.join(dd_dir, "drawdown_%s.scal_vec.0.vtk" % s))
    v = np.maximum(f["StrainPlasticJ2"], 0)
    tp = ax.tripcolor(pts[:, 0], pts[:, 1], tris, v / v.max(), shading="gouraud", cmap="inferno_r", vmin=0, vmax=0.2)
    ax.triplot(pts[:, 0], pts[:, 1], tris, lw=0.1, color="0.75")
    outline(ax, zoom=True)
    if s in bish:
        circle(ax, *bish[s][1], color="c", lw=1.0, ls="--")
    ax.set_title("SRM, %s" % titles[s])
fig.colorbar(tp, ax=axs, shrink=0.8)
save(fig, "fig_mechanism_drawdown.pdf")

# ---------------------------------------------------------------- independent seepage solution
for name, fname in [("steady", "vertices_steady.csv"), ("crest", "vertices_crest.csv")]:
    p = os.path.join(laplace_dir, fname)
    if not os.path.exists(p):
        continue
    d = np.genfromtxt(p, delimiter=",", names=True)
    print(name, "max |p_neopz - p_python| =", np.max(np.abs(d["p_neopz"] - d["p_python"])))
if os.path.exists(os.path.join(laplace_dir, "vertices_steady.csv")):
    d = np.genfromtxt(os.path.join(laplace_dir, "vertices_steady.csv"), delimiter=",", names=True)
    fig, axs = plt.subplots(1, 2, figsize=(6.8, 2.6))
    axs[0].plot(d["p_python"], d["p_neopz"], ".", ms=2)
    axs[0].plot([0, d["p_python"].max()], [0, d["p_python"].max()], "k-", lw=0.6)
    axs[0].set_xlabel("p, P1 independente (kPa)")
    axs[0].set_ylabel("p, NeoPZ u-p (kPa)")
    axs[0].set_aspect("equal")
    axs[0].grid(alpha=0.3)
    tp = axs[1].tricontourf(d["x"], d["y"], d["p_neopz"] - d["p_python"], 21, cmap="RdBu_r")
    outline(axs[1])
    fig.colorbar(tp, ax=axs[1], shrink=0.8, label="diferença (kPa)")
    axs[1].set_title("NeoPZ $-$ independente")
    save(fig, "fig_laplace_check.pdf")
print("figures in", out)
