"""Figures 8, 9 and 10 of the article (Abaqus benchmark 1.15.2, Sect. 6.5) and two supplementary figures (softening
study and global iterations of Sect. 6.7) from the files written by AbaqusTriaxialConsolidation.

Usage (run the executable first; it writes its CSV files to the current directory):

    python3 <neopz>/Projects/AbaqusTriaxialConsolidation/plot_figures.py [run directory] [-o output directory]

The run directory (default: the current directory) holds the CSV files; the figures are written as PDF and PNG
to <run directory>/figures, or to the output directory. A figure whose CSV files are missing (the executable was
run with only some of its parts) is skipped with a message.

Figures produced:

- fig08_abaqus_model: Fig. 8, the model of the upper half of the specimen: (a) quarter of the specimen with
  2 x 2 x 4 divisions (48 Hex20-Hex8 elements, curved lateral faces), the boundary conditions (platen, cell
  pressure on the lateral face, symmetry planes x = 0 and y = 0, mid-plane), the vertex (u and p_w) and mid-edge
  (u only) nodes of the visible faces and the point A; (b) the refined mesh 4 x 4 x 8 (384 elements) of the
  softening study. The meshes are read from the CSV files of mcc::WriteMeshCSV (abaqus_mesh_2x2x4_*.csv and
  abaqus_mesh_4x4x8_*.csv, part "mesh"); the curved edges are the parabolas through their mid-edge nodes. The
  drawing is an orthographic projection with back-face culling (the quarter of cylinder is convex).
- fig09_abaqus_states: Fig. 9, drained tests at a material point from p'0 = 100 kPa (R = 1.17, subcritical) and
  p'0 = 20 kPa (R = 5.83, supercritical) with p'c0 = 116.6 kPa: (a) stress paths p'-q, (b) q-eps_1 and
  (c) eps_v-eps_1 in the four-quadrant layout; 600 increments (abaqus_material_point_p0_<p'0>.csv) against the
  closed form of Appendix B.1 (abaqus_closed_form_p0_<p'0>.csv). Part "mp".
- fig10_abaqus_results: Fig. 10, q at A against delta/H: (a) smooth platen, Hex20-Hex8 model with 150 increments
  (abaqus_smooth_2x2x2.csv) and material point with 30 increments (abaqus_material_point_30.csv); (b) rough
  platen with 2 x 2 x 2 and 3 x 3 x 3 points (abaqus_rough_2x2x2.csv, abaqus_rough_3x3x3.csv), with the range of q
  at the 27 integration points of the element that contains A for the full integration (columns q_elA_min and
  q_elA_max: the oscillation of the locked stresses; A lies on an edge of that element); (c) stress paths at A. Markers: the Abaqus curves of Figs. 1.15.2-2 and -3 of the Abaqus Benchmarks Manual, digitized
  (reference/abaqus_1_15_2_digitalizado.json). Parts "mp" and "fe".
- supplementary_abaqus_softening: smooth platen from p'0 = 20 kPa (600 increments): (a) global deviatoric stress
  sigma_a - p'0 from the platen force and q at A with the meshes 2 x 2 x 4 and 4 x 4 x 8, against the material
  point; (b) deformed outer generatrix (y = 0, r = R) at delta/H = 0.6 (abaqus_smooth_600_p0_20*.csv,
  abaqus_smooth_600_p0_20*_profile.csv). Parts "states" and "softening" (refined mesh, drawn if present).
- supplementary_abaqus_tangents: rough platen, cumulative number of evaluations of the residual against delta/H
  with the five tangent operators of Table 10: (a) tolerance 1e-8 (abaqus_tangents_evaluations.csv,
  abaqus_table10.csv, part "tangents"); (b) tolerance 1e-9, drawn if the part "tolerance" was run
  (abaqus_table10_tolerance.csv and the histories abaqus_rough_2x2x2_tolerance_tangent_<operator>.csv).

The figures use only the CSV files; the VTK series of the executable are not needed (the argument "novtk" of the
executable disables them).
"""
import csv
import os
import sys

import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch
from matplotlib.patches import Polygon as MPoly

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'Common'))
from mcc_figstyle import C1, C2, C3, C4, INK, INK2, MUTED, TEXTW, arguments, plt, read_csv, reference, save  # noqa

# AbaqusTriaxialConsolidation.h: upper half of the specimen (mm), point A and the material of the benchmark
R_MM, H_MM = 20.0, 60.0
A_MM = (5.0, 0.0, 7.5)
M, PC0 = 1.0, 116.6
# boundary ids of AbaqusTriaxialConsolidation.h
EX0, ELATERAL, EY0, EZ0, EZH, EPZH = -21, -22, -23, -25, -26, -36


def read_table(path):
    """CSV file with text columns (summary and numbers files): list of dicts."""
    if not os.path.exists(path):
        sys.exit(f'file not found: {path}')
    with open(path, newline='') as f:
        return list(csv.DictReader(f))


# =========================================================================================== Fig. 8 model
def read_mesh(run, prefix):
    """Mesh written by mcc::WriteMeshCSV: node coordinates (mm), boundary faces (node lists) and their ids."""
    N = read_csv(os.path.join(run, prefix + '_nodes.csv'))
    X = 1e3 * np.column_stack([N['x'], N['y'], N['z']])
    F = read_csv(os.path.join(run, prefix + '_faces.csv'))
    E = read_csv(os.path.join(run, prefix + '_elements.csv'))
    faces = np.column_stack([F[f'n{k}'] for k in range(8)]).astype(int)
    return X, faces, F['matid'].astype(int), len(E['element'])


class View:
    """Orthographic projection with the camera of matplotlib's view_init(elev, azim)."""

    def __init__(self, elev, azim):
        e, a = np.radians(elev), np.radians(azim)
        self.d = np.array([np.cos(e) * np.cos(a), np.cos(e) * np.sin(a), np.sin(e)])   # towards the camera
        self.e1 = np.array([-np.sin(a), np.cos(a), 0.0])
        self.e2 = np.cross(self.d, self.e1)

    def __call__(self, P):
        P = np.atleast_2d(P)
        return np.column_stack([P @ self.e1, P @ self.e2])

    def visible(self, n):
        return float(np.dot(n, self.d)) > 1e-9


def face_normal(X, f, mid):
    """Outward normal of a boundary face of the quarter of cylinder (from its boundary id)."""
    if mid == EX0:
        return np.array([-1.0, 0, 0])
    if mid == EY0:
        return np.array([0, -1.0, 0])
    if mid == EZ0:
        return np.array([0, 0, -1.0])
    if mid in (EZH, EPZH):
        return np.array([0, 0, 1.0])
    c = X[f[:4]].mean(axis=0)
    return np.array([c[0], c[1], 0.0]) / np.hypot(c[0], c[1])


def face_ring(X, f, n=6):
    """Boundary of a Hex20 face: the four edges as parabolas through their mid-edge nodes (nodes 4-7)."""
    t = np.linspace(0, 1, n + 1)
    Nq = np.c_[(1 - t) * (1 - 2 * t), 4 * t * (1 - t), t * (2 * t - 1)]
    ring = []
    for i, j, m in ((0, 1, 4), (1, 2, 5), (2, 3, 6), (3, 0, 7)):
        ring += list((Nq @ X[[f[i], f[m], f[j]]])[:-1])
    return np.array(ring)


FACECOL = {EX0: '#e9e8e3', EY0: '#e9e8e3', ELATERAL: '#c4e5d6', EZ0: '#d9d8d2', EZH: '#b7d3f6'}


def draw_mesh(ax, X, faces, ids, view, lw=0.45, nodes=False, shift=(0.0, 0.0)):
    """Draws the visible boundary faces (the platen once: its two coincident boundary elements have the ids EZH and
    EPZH), shifted on the screen by shift; with nodes, also the vertex and mid-edge nodes of the visible faces.
    Returns the screen bounding box (xmin, xmax, ymin, ymax) of the visible faces."""
    seen = set()
    vis_nodes = set()
    corner = set()
    box = [np.inf, -np.inf, np.inf, -np.inf]
    for f, mid in zip(faces, ids):
        corner.update(int(v) for v in f[:4])
        if mid == EPZH:
            continue
        key = tuple(sorted(f[:4]))
        if key in seen:
            continue
        seen.add(key)
        if not view.visible(face_normal(X, f, mid)):
            continue
        P = view(face_ring(X, f)) + shift
        ax.add_patch(MPoly(P, closed=True, fc=FACECOL[mid], ec=INK2, lw=lw, joinstyle='round'))
        box = [min(box[0], P[:, 0].min()), max(box[1], P[:, 0].max()), min(box[2], P[:, 1].min()),
               max(box[3], P[:, 1].max())]
        vis_nodes.update(int(v) for v in f)
    if nodes:
        V = sorted(vis_nodes)
        mm = [v for v in V if v not in corner]
        vv = [v for v in V if v in corner]
        P = view(X[mm]) + shift
        ax.plot(P[:, 0], P[:, 1], 'o', ms=1.3, color=INK2, zorder=4)
        P = view(X[vv]) + shift
        ax.plot(P[:, 0], P[:, 1], 's', ms=2.0, color=C1, zorder=5)
    return box


def fig_abaqus_model(run, out):
    """Fig. 8: (a) the 3D model with 48 Hex20-Hex8 elements and its boundary conditions; (b) refined mesh."""
    X, faces, ids, nel = read_mesh(run, 'abaqus_mesh_2x2x4')
    fine = os.path.exists(os.path.join(run, 'abaqus_mesh_4x4x8_faces.csv'))
    view = View(22, -28)
    fig = plt.figure(figsize=(TEXTW * 0.80, 3.45))
    ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
    box = draw_mesh(ax, X, faces, ids, view, nodes=True)
    P = lambda *x: view(np.array(x, dtype=float))[0]
    sil = np.radians(90.0 - 28.0)          # silhouette of the lateral face: radial direction parallel to the screen
    # axis of the specimen (edge x = y = 0), extended above the platen
    a0, a1 = P(0, 0, -3), P(0, 0, 72)
    ax.plot([a0[0], a1[0]], [a0[1], a1[1]], color=INK2, lw=0.6, ls=(0, (6, 2, 1, 2)), zorder=3)
    ax.text(*P(0, 0, 73), 'axis', fontsize=6.6, color=INK2, ha='center', va='bottom')
    # platen: prescribed vertical displacement (arrows above the top face)
    for (x, y) in ((4, 1.5), (12, 1.5), (18.5, 2.5), (9, 9), (2.5, 12), (2.5, 18.5), (13.5, 13.5)):
        p0, p1 = P(x, y, 69.0), P(x, y, 61.0)
        ax.add_patch(FancyArrowPatch(tuple(p0), tuple(p1), arrowstyle='-|>', mutation_scale=5.5, color=INK, lw=0.7,
                                     zorder=6))
    # cell pressure: radial arrows pointing to the silhouette of the lateral face
    for z in np.linspace(4, 56, 7):
        c, s = np.cos(sil), np.sin(sil)
        p0, p1 = P(27.0 * c, 27.0 * s, z), P(20.8 * c, 20.8 * s, z)
        ax.add_patch(FancyArrowPatch(tuple(p0), tuple(p1), arrowstyle='-|>', mutation_scale=5.5, color=C2, lw=0.7,
                                     zorder=6))
    # point A on the symmetry plane y = 0
    pa = P(*A_MM)
    ax.plot(*pa, '*', ms=8.5, color=C3, mec=INK, mew=0.4, zorder=8)
    ax.annotate('A', xy=tuple(pa), xytext=(-8, 3), textcoords='offset points', fontsize=8, fontweight='bold',
                ha='center', va='center', zorder=9)
    # dimensions
    dim = dict(arrowstyle='<->', color=INK2, lw=0.6, mutation_scale=6, shrinkA=0, shrinkB=0)
    d0, d1 = P(0, -6, 0), P(20, -6, 0)
    ax.annotate('', xy=tuple(d0), xytext=tuple(d1), arrowprops=dim)
    for x in (0, 20):
        e0, e1 = P(x, -0.8, 0), P(x, -7, 0)
        ax.plot([e0[0], e1[0]], [e0[1], e1[1]], color=INK2, lw=0.4)
    ang = np.degrees(np.arctan2(d1[1] - d0[1], d1[0] - d0[0]))
    ax.annotate('$R$ = 20 mm', xy=tuple((d0 + d1) / 2), xytext=(0, -2), textcoords='offset points', fontsize=6.8,
                ha='center', va='top', rotation=ang, rotation_mode='anchor')
    h0, h1 = P(0, -14, 0), P(0, -14, 60)
    ax.annotate('', xy=tuple(h0), xytext=tuple(h1), arrowprops=dim)
    ax.text(*((h0 + h1) / 2 + [-1.0, 0]), '$H$ = 60 mm (upper half)', rotation=90, fontsize=6.8, ha='right',
            va='center')
    for z in (0, 60):   # extension lines of the dimension H
        e0, e1 = P(0, -1.0, z), P(0, -15, z)
        ax.plot([e0[0], e1[0]], [e0[1], e1[1]], color=INK2, lw=0.4)
    # labels of the boundary conditions, in a column to the right of the mesh
    xt = box[1] + 9.0
    lab = [(P(10, 6, 60), P(0, 0, 70)[1], 'platen: $u_z$ prescribed, $p_w = 0$\n(rough platen: also $u_x = u_y = 0$)',
            C1),
           (P(20 * np.cos(sil - 0.35), 20 * np.sin(sil - 0.35), 36), P(0, 0, 44)[1],
            "lateral face: cell pressure\n$P = p'_0$ = 100 kPa", C2),
           (P(14, 0, 22), P(0, 0, 21)[1], 'symmetry planes $x = 0$ (hidden)\nand $y = 0$: $u_x = 0$, $u_y = 0$',
            INK2)]
    for xy, y, s, c in lab:
        ax.annotate(s, xy=tuple(xy), xytext=(xt, y), textcoords='data', fontsize=6.8, color=INK, ha='left',
                    va='center', arrowprops=dict(arrowstyle='-', color=c, lw=0.7, shrinkA=2, shrinkB=0), zorder=7)
    y = P(0, 0, 4)[1]
    ax.text(xt, y, 'mid-plane $z = 0$ (hidden):\n$u_z = 0$, impermeable', fontsize=6.8, color=INK2, ha='left',
            va='center')
    ax.text(xt, y - 8.5, 'A: $x$ = 5 mm, $y$ = 0, $z$ = 7.5 mm', fontsize=6.8, color=INK, ha='left', va='center')
    yl = y - 15.0
    ax.plot(xt + 0.8, yl, 's', ms=2.0, color=C1)
    ax.text(xt + 2.6, yl, 'vertex node: $u$ and $p_w$', fontsize=6.6, va='center')
    ax.plot(xt + 0.8, yl - 3.6, 'o', ms=1.3, color=INK2)
    ax.text(xt + 2.6, yl - 3.6, 'mid-edge node: $u$ only', fontsize=6.6, va='center')
    ytitle = P(0, 0, 70)[1] + 9.0
    ax.text(box[0] - 14, ytitle, '(a) %d elements, 2 × 2 × 4 divisions' % nel, fontsize=8, va='center')
    nfine, right = 0, xt + 44.0
    if fine:
        Xf, ff, idf, nfine = read_mesh(run, 'abaqus_mesh_4x4x8')
        sh = (right - box[0], 0.0)
        boxf = draw_mesh(ax, Xf, ff, idf, view, lw=0.3, shift=sh)
        ax.text(boxf[1], ytitle, '(b) %d elements, 4 × 4 × 8 divisions' % nfine, fontsize=8, va='center',
                ha='right')
        right = boxf[1]
    ax.set_aspect('equal')
    ax.set_xlim(box[0] - 15, right + 1.5)
    ax.set_ylim(yl - 6, ytitle + 4)
    ax.axis('off')
    print('Fig. 8: coarse mesh %d elements, %d nodes; refined mesh %d elements' % (nel, len(X), nfine))
    save(fig, out, 'fig08_abaqus_model')


# =========================================================================================== Fig. 9 two initial states
def material_point_states(run):
    """Material point (600 increments) and closed form of the two initial states: p'0 -> (point, closed), arrays
    with the columns eps_1, p', q, eps_v (compression positive)."""
    Z = {}
    for p0 in (100.0, 20.0):
        cur = []
        for kind in ('material_point', 'closed_form'):
            name = f'abaqus_{kind}_p0_{p0:.0f}.csv'
            T = read_csv(os.path.join(run, name))
            A = np.column_stack([T['eps_1'], T['p'], T['q'], T['eps_v']])
            if abs(A[0, 1] - p0) > 1e-9 * p0:
                sys.exit(f"{name} starts at p' = {A[0, 1]:g} kPa, expected {p0:g} kPa")
            cur.append(A)
        Z[p0] = dict(point=cur[0], closed=cur[1])
    return Z


def fig_abaqus_states(run, out):
    """Fig. 9: four quadrants, p'-q (right), q-eps_1 (eps_1 to the left) and eps_v-eps_1 (eps_v downwards)."""
    Z = material_point_states(run)
    cases = [(100.0, C1), (20.0, C2)]
    fig = plt.figure(figsize=(TEXTW * 0.92, 4.6))
    gs = fig.add_gridspec(2, 2, width_ratios=(1, 1), height_ratios=(1, 0.78), wspace=0.035, hspace=0.035,
                          left=0.09, right=0.985, top=0.985, bottom=0.10)
    axQ = fig.add_subplot(gs[0, 1])
    axE = fig.add_subplot(gs[0, 0], sharey=axQ)
    axV = fig.add_subplot(gs[1, 0], sharex=axE)
    axT = fig.add_subplot(gs[1, 1])
    axT.axis('off')
    dash = dict(color=INK, lw=0.8, ls=(0, (3, 2)))
    # (a) stress paths in the p'-q plane
    pp = np.linspace(0, PC0, 400)
    axQ.plot(pp, M * np.sqrt(np.maximum(pp * (PC0 - pp), 0)), color=MUTED, lw=0.9, ls=(0, (4, 2)))
    axQ.plot([0, 165], [0, 165 * M], color=MUTED, lw=0.8)
    for p0, col in cases:
        mp, cf = Z[p0]['point'], Z[p0]['closed']
        axQ.plot(mp[:, 1], mp[:, 2], color=col, lw=1.6)
        axQ.plot(cf[:, 1], cf[:, 2], **dash)
        axQ.plot(p0, 0, 'o', ms=4, color=col, clip_on=False, zorder=5)
        axQ.plot(1.5 * p0, 1.5 * p0, 'o', ms=4.5, mfc='white', mec=col, mew=1.0, zorder=5)
    k = int(np.argmax(Z[20.0]['closed'][:, 2]))
    pk = Z[20.0]['closed'][k]
    axQ.annotate(f'peak, $q$ = {pk[2]:.1f} kPa', xy=(pk[1], pk[2]), xytext=(6, 98), fontsize=6.8, color=INK,
                 ha='left', va='bottom', arrowprops=dict(arrowstyle='-', color=MUTED, lw=0.6))
    axQ.annotate('critical state', xy=(150, 150), xytext=(100, 157), fontsize=6.8, color=INK, ha='left',
                 va='center', arrowprops=dict(arrowstyle='-', color=MUTED, lw=0.6))
    axQ.text(95, 99, 'CSL', fontsize=6.8, color=INK2, ha='right', va='bottom')
    axQ.text(75, 22, 'initial yield surface', fontsize=6.8, color=INK2, ha='center')
    axQ.text(97, 3, "$p'_0$ = 100", fontsize=6.6, color=INK, ha='right', va='bottom')
    axQ.text(23, 3, "$p'_0$ = 20", fontsize=6.6, color=INK, ha='left', va='bottom')
    axQ.set(xlim=(0, 165), ylim=(0, 165), xlabel="$p'$ (kPa)")
    axQ.text(0.025, 0.975, '(a)', transform=axQ.transAxes, ha='left', va='top', fontsize=8.5)
    # (b) q against eps_1 (eps_1 grows to the left)
    for p0, col in cases:
        mp, cf = Z[p0]['point'], Z[p0]['closed']
        axE.plot(100 * mp[:, 0], mp[:, 2], color=col, lw=1.6)
        axE.plot(100 * cf[:, 0], cf[:, 2], **dash)
    axE.text(40, 132, 'hardening', fontsize=7, color=INK, ha='center')
    axE.text(40, 17, 'softening', fontsize=7, color=INK, ha='center')
    axE.set_xlim(60, 0)
    axE.set_ylabel('$q$ (kPa)')
    axE.text(0.975, 0.975, '(b)', transform=axE.transAxes, ha='right', va='top', fontsize=8.5)
    # (c) eps_v against eps_1 (eps_v, compression positive, grows downwards)
    axV.axhline(0, color=MUTED, lw=0.6, ls=(0, (1, 1.5)))
    for p0, col in cases:
        mp, cf = Z[p0]['point'], Z[p0]['closed']
        axV.plot(100 * mp[:, 0], 100 * mp[:, 3], color=col, lw=1.6)
        axV.plot(100 * cf[:, 0], 100 * cf[:, 3], **dash)
    axV.text(30, -4.4, 'dilation', fontsize=7, color=INK, ha='center', va='bottom')
    axV.set_ylim(8, -5.5)
    axV.set(xlabel=r'$\varepsilon_1$ (%)', ylabel=r'$\varepsilon_v$ (%)')
    axV.text(0.975, 0.04, '(c)', transform=axV.transAxes, ha='right', va='bottom', fontsize=8.5)
    # shared axes: tick labels on the outer sides only; the inner axes drawn as in the classical scheme
    plt.setp(axQ.get_yticklabels(), visible=False)
    plt.setp(axE.get_xticklabels(), visible=False)
    for ax, sides in ((axQ, ('left', 'bottom')), (axE, ('right', 'bottom')), (axV, ('right', 'top'))):
        for sd in sides:
            ax.spines[sd].set_color(INK)
            ax.spines[sd].set_linewidth(1.0)
    # (d) legend and data
    hh = [Line2D([], [], color=C1, lw=1.6, label="$p'_0$ = 100 kPa, $R$ = %.2f (subcritical)" % (PC0 / 100.0)),
          Line2D([], [], color=C2, lw=1.6, label="$p'_0$ = 20 kPa, $R$ = %.2f (supercritical)" % (PC0 / 20.0)),
          Line2D([], [], **dash, label='closed form'),
          Line2D([], [], color=MUTED, lw=0.9, ls=(0, (4, 2)),
                 label="initial yield surface, $p'_{c0}$ = %g kPa" % PC0),
          Line2D([], [], color=MUTED, lw=0.8, label='critical state line, $q = Mp\'$')]
    axT.legend(handles=hh, loc='upper left', bbox_to_anchor=(0.04, 0.80), fontsize=6.9, handlelength=2.0,
               labelspacing=0.45, borderaxespad=0)
    axT.text(0.06, 0.06, "drained, $\\sigma_3 = p'_0$ constant, material point\n$M$ = 1, $\\lambda$ = 0.174, "
             "$\\kappa$ = 0.026, $v_0$ = 2.08, $\\nu$ = 0.3", transform=axT.transAxes, fontsize=6.9, color=INK2,
             va='bottom', linespacing=1.4)
    for ax in (axQ, axE, axV):
        ax.tick_params(labelsize=7)
    for p0, _ in cases:
        mp, cf = Z[p0]['point'], Z[p0]['closed']
        kq = int(np.argmax(mp[:, 2]))
        print("Fig. 9: p'0 = %g kPa, peak q = %.3f kPa at eps_1 = %.3f (closed form %.3f); end q = %.3f, eps_v = %.4f;"
              " max |q - closed| = %.3f kPa" % (p0, mp[kq, 2], mp[kq, 0], cf[:, 2].max(), mp[-1, 2], mp[-1, 3],
                                               np.abs(mp[:, 2] - np.interp(mp[:, 0], cf[:, 0], cf[:, 2])).max()))
    save(fig, out, 'fig09_abaqus_states')


# =========================================================================================== Fig. 10 results
COLS = ('delta_H', 'p_A', 'q_A', 'sigma_a_platen', 'q_platen', 'max_pw', 'eps_v', 'evaluations', 'q_elA_min',
        'q_elA_mean', 'q_elA_max')


def history(run, name):
    """History of a finite element run (abaqus_<name>.csv): columns delta/H, p' at A, q at A, axial stress on the
    platen, sigma_a - p'0, largest |p_w|, eps_v, evaluations of the increment and the smallest, mean and largest q
    at the integration points of the element that contains A."""
    T = read_csv(os.path.join(run, f'abaqus_{name}.csv'))
    missing = [k for k in COLS if k not in T]
    if missing:
        sys.exit(f'abaqus_{name}.csv has no column {", ".join(missing)}: rerun the current AbaqusTriaxialConsolidation')
    return np.column_stack([T[k] for k in COLS])


def fig_abaqus(run, out):
    """Fig. 10: q at A against delta/H, smooth (a) and rough (b) platens, and the stress paths at A (c)."""
    ab = reference(HERE, 'abaqus_1_15_2_digitalizado.json')
    fig, axs = plt.subplots(1, 3, figsize=(TEXTW, 2.35))
    hs, hr, hr4 = history(run, 'smooth_2x2x2'), history(run, 'rough_2x2x2'), history(run, 'rough_3x3x3')
    T = read_csv(os.path.join(run, 'abaqus_material_point_30.csv'))     # material point, 30 increments (0.02)
    mp = np.column_stack([T['delta_H'], T['p'], T['q'], T['eps_v']])
    ax = axs[0]
    ax.axhline(150.0, color=MUTED, lw=0.7, ls=(0, (1, 1.5)))
    ax.text(0.01, 153, 'critical state: $q$ = 150 kPa', fontsize=6.8, color=INK2, va='bottom')
    ax.plot(hs[:, 0], hs[:, 2], color=C1, lw=1.5, label='Hex20–Hex8, %d incr.' % (len(hs) - 1))
    ax.plot(mp[:, 0], mp[:, 2], color=INK, lw=0.8, ls=(0, (3, 2)), label='material point, %d incr.' % (len(mp) - 1))
    d = np.array(ab['lisa']['qd'])
    ax.plot(d[:, 0], d[:, 1], 'o', ms=3.8, mfc='white', mec=C2, mew=0.9, label='Abaqus')
    ax.set(xlabel='$\\delta/H$', ylabel='$q$ at A (kPa)', xlim=(0, 0.6), ylim=(0, 170))
    ax.set_title('(a) smooth platen (homogeneous)', fontsize=8)
    ax.legend(loc='lower right', fontsize=6.8)
    ax = axs[1]
    # range of q at the 27 points of the element that contains A (A lies on an edge of this element)
    ax.fill_between(hr4[:, 0], hr4[:, 8], hr4[:, 10], color=MUTED, alpha=0.3, lw=0, zorder=1,
                    label='3×3×3: range in the\nelement of A')
    ax.plot(hr4[:, 0], hr4[:, 2], color=MUTED, lw=1.1, label='3×3×3 points (full)')
    ax.plot(hr[:, 0], hr[:, 2], color=C3, lw=1.5, label='2×2×2 points (reduced)')
    d = np.array(ab['rugosa']['qd'])
    ax.plot(d[:, 0], d[:, 1], 's', ms=3.6, mfc='white', mec=C2, mew=0.9, label='Abaqus')
    k = int(np.argmin(abs(hr4[:, 0] - 0.55)))
    ax.annotate('3×3×3: volumetric\nlocking', xy=(hr4[k, 0], hr4[k, 2]), xytext=(0.40, 112), fontsize=6.8,
                color=INK, ha='left', va='top', arrowprops=dict(arrowstyle='-', color=MUTED, lw=0.7))
    ax.set(xlabel='$\\delta/H$', ylabel='$q$ at A (kPa)', xlim=(0, 0.6), ylim=(0, 170))
    ax.set_title('(b) rough platen', fontsize=8)
    ax.legend(loc='lower right', fontsize=6.8)
    ax = axs[2]
    ax.plot(hs[:, 1], hs[:, 2], color=C1, lw=1.5, label='smooth')
    ax.plot(hr[:, 1], hr[:, 2], color=C3, lw=1.5, label='rough')
    d = np.array(ab['lisa']['pq'])
    ax.plot(d[:, 0], d[:, 1], 'o', ms=3.8, mfc='white', mec=C2, mew=0.9, label='smooth, Abaqus')
    d = np.array(ab['rugosa']['pq'])
    ax.plot(d[:, 0], d[:, 1], 's', ms=3.6, mfc='white', mec=C2, mew=0.9, label='rough, Abaqus')
    pm = np.linspace(95, 175, 2)
    ax.plot(pm, M * pm, color=MUTED, lw=0.8)
    ax.text(127, 133, "CSL: $q = Mp'$", fontsize=6.8, color=INK2, ha='right', va='bottom')
    ax.set(xlabel="$p'$ at A (kPa)", ylabel='$q$ at A (kPa)', xlim=(95, 175), ylim=(0, 170))
    ax.set_title("(c) stress paths at A", fontsize=8)
    ax.legend(loc='lower right', fontsize=6.8, handlelength=1.4, borderaxespad=0.3, labelspacing=0.35)
    fig.tight_layout(w_pad=0.9)
    # label along the drained path (angle measured on the screen, after the layout)
    fig.canvas.draw()
    a1, a2 = ax.transData.transform((100, 0)), ax.transData.transform((150, 150))
    ang = np.degrees(np.arctan2(a2[1] - a1[1], a2[0] - a1[0]))
    off = 4.0 * np.array([-np.sin(np.radians(ang)), np.cos(np.radians(ang))])
    ax.annotate('drained: slope 3', xy=(114.5, 43.5), xytext=tuple(off), textcoords='offset points',
                rotation=ang, rotation_mode='anchor', ha='center', va='bottom', fontsize=6.6, color=INK)
    print('Fig. 10: q_A at delta/H = 0.6: smooth %.3f, material point (30 incr.) %.3f, rough 2x2x2 %.3f, rough 3x3x3 '
          '%.3f (max %.3f at %.3f) kPa' % (hs[-1, 2], mp[-1, 2], hr[-1, 2], hr4[-1, 2], hr4[:, 2].max(),
                                          hr4[np.argmax(hr4[:, 2]), 0]))
    save(fig, out, 'fig10_abaqus_results')


# =========================================================================================== supplementary figures
def fig_softening(run, out):
    """Softening state p'0 = 20 kPa (smooth platen, 600 increments): global and local deviatoric stresses with the
    two meshes against the material point, and the deformed outer generatrix at the end."""
    T = read_csv(os.path.join(run, 'abaqus_material_point_p0_20.csv'))
    meshes = [('smooth_600_p0_20', '2×2×4 (48 elements)', C1), ('smooth_600_p0_20_mesh4x4x8', '4×4×8 (384 elements)',
                                                                  C3)]
    meshes = [m for m in meshes if os.path.exists(os.path.join(run, f'abaqus_{m[0]}.csv'))]
    fig, axs = plt.subplots(1, 2, figsize=(TEXTW * 0.9, 2.6), gridspec_kw=dict(width_ratios=(1.5, 1)))
    ax = axs[0]
    ax.plot(T['eps_1'], T['q'], color=INK, lw=0.9, ls=(0, (3, 2)), label='material point')
    for name, lab, col in meshes:
        h = history(run, name)
        ax.plot(h[:, 0], h[:, 4], color=col, lw=1.4, label=f"$\\sigma_a - p'_0$, {lab}")
        ax.plot(h[:, 0], h[:, 2], color=col, lw=0.8, ls=(0, (1, 1.2)), label=f'$q$ at A, {lab}')
        print("softening, %s: global q peak %.3f at %.3f, end %.3f; q_A end %.3f (material point %.3f)"
              % (name, h[:, 4].max(), h[np.argmax(h[:, 4]), 0], h[-1, 4], h[-1, 2], T['q'][-1]))
    ax.set(xlabel='$\\delta/H$', ylabel='deviatoric stress (kPa)', xlim=(0, 0.6), ylim=(0, 60))
    ax.set_title("(a) $p'_0$ = 20 kPa, smooth platen", fontsize=8)
    ax.legend(loc='lower right', fontsize=6.4, handlelength=1.8, labelspacing=0.3)
    ax = axs[1]
    ax.plot([R_MM, R_MM], [0, H_MM], color=MUTED, lw=0.8, ls=(0, (3, 2)), label='initial')
    for name, lab, col in meshes:
        P = read_csv(os.path.join(run, f'abaqus_{name}_profile.csv'))
        ax.plot(1e3 * (0.02 + P['u_r']), 1e3 * (P['z'] + P['u_z']), color=col, lw=1.4, label=lab.split(' ')[0])
    ax.set(xlabel='$r$ (mm)', ylabel='$z$ (mm)', xlim=(15, 45), ylim=(0, 62))
    ax.set_aspect('equal', adjustable='box')
    ax.set_title('(b) outer generatrix at $\\delta/H$ = 0.6', fontsize=8)
    ax.legend(loc='upper right', fontsize=6.4, handlelength=1.4)
    fig.tight_layout(w_pad=1.0)
    save(fig, out, 'supplementary_abaqus_softening')


def fig_tangents(run, out):
    """Rough platen: cumulative evaluations of the residual with the five operators of Table 10, with the tolerance
    1e-8 of the article (a) and, if the part "tolerance" was run, with 1e-9 (b)."""
    names = {'D': ('consistent $D$', C1, '-'), 'fd': ('finite differences', INK, (0, (2, 2))),
             'sym': ('$(D + D^T)/2$', C2, '-'), 'cont': ('continuum tangent', C4, '-'), 'DT': ('$D^T$', C3, '-')}
    T = read_csv(os.path.join(run, 'abaqus_tangents_evaluations.csv'))
    rows = {r['tangent']: r for r in read_table(os.path.join(run, 'abaqus_table10.csv'))}
    curves = [(rows, {k: T[k] for k in names if k in T})]
    tolfile = os.path.join(run, 'abaqus_table10_tolerance.csv')
    if os.path.exists(tolfile):
        rt = {r['tangent']: r for r in read_table(tolfile)}
        curves.append((rt, {k: history(run, r['run'])[1:, 7] for k, r in rt.items()}))
    fig, axs = plt.subplots(1, len(curves), figsize=(TEXTW * (0.55 if len(curves) == 1 else 1.0), 2.4), squeeze=False)
    top = max(np.sum(ev) for _, E in curves for ev in E.values())
    for i, (ax, (rr, E)) in enumerate(zip(axs[0], curves)):
        tol = float(next(iter(rr.values()))['tolerance']) if 'tolerance' in next(iter(rr.values())) else 1e-8
        base = float(rr['D']['mean_evaluations'])
        for key, (lab, col, ls) in names.items():
            if key not in E:
                continue
            r = rr[key]
            m = float(r['mean_evaluations'])
            ax.plot(np.r_[0, T['delta_H']], np.r_[0, np.cumsum(E[key])], color=col, lw=1.3 if key != 'fd' else 1.0,
                    ls=ls, label='%s: %.2f (%s), %.2f' % (lab, m, r['max_evaluations'], m / base))
        ax.set(xlabel='$\\delta/H$', ylabel='cumulative evaluations', xlim=(0, 0.6), ylim=(0, 1.05 * top))
        ax.set_title('(%s) tolerance $10^{%d}$' % ('ab'[i], round(np.log10(tol))), fontsize=8)
        ax.legend(loc='upper left', fontsize=6.4, title='mean per increment (largest), ratio to $D$',
                  title_fontsize=6.4)
        print('Table 10, tolerance %g: ' % tol + ', '.join('%s %.3f (ratio %.3f)' % (k, float(rr[k]['mean_evaluations']),
                                                                                      float(rr[k]['mean_evaluations']) / base)
                                                          for k in names if k in rr))
    fig.tight_layout()
    save(fig, out, 'supplementary_abaqus_tangents')


# parts of the executable (command line arguments, in the order of the executable)
PARTS = ('mesh', 'mp', 'fe', 'tangents', 'tolerance', 'states', 'softening')
# (function, {CSV file: part of the executable that writes it})
FIGURES = (
    (fig_abaqus_model, {'abaqus_mesh_2x2x4_nodes.csv': 'mesh', 'abaqus_mesh_2x2x4_faces.csv': 'mesh',
                        'abaqus_mesh_2x2x4_elements.csv': 'mesh'}),
    (fig_abaqus_states, {'abaqus_material_point_p0_100.csv': 'mp', 'abaqus_closed_form_p0_100.csv': 'mp',
                         'abaqus_material_point_p0_20.csv': 'mp', 'abaqus_closed_form_p0_20.csv': 'mp'}),
    (fig_abaqus, {'abaqus_smooth_2x2x2.csv': 'fe', 'abaqus_rough_2x2x2.csv': 'fe', 'abaqus_rough_3x3x3.csv': 'fe',
                  'abaqus_material_point_30.csv': 'mp'}),
    (fig_softening, {'abaqus_smooth_600_p0_20.csv': 'states', 'abaqus_smooth_600_p0_20_profile.csv': 'states',
                     'abaqus_material_point_p0_20.csv': 'mp'}),
    (fig_tangents, {'abaqus_tangents_evaluations.csv': 'tangents', 'abaqus_table10.csv': 'tangents'}),
)

if __name__ == '__main__':
    args = arguments('Figs. 8 to 10 of the article (Abaqus benchmark 1.15.2) and the supplementary figures from the '
                     'CSV files of AbaqusTriaxialConsolidation.')
    for function, files in FIGURES:
        missing = [f for f in files if not os.path.exists(os.path.join(args.rundir, f))]
        if missing:
            parts = ' '.join(p for p in PARTS if any(files[f] == p for f in missing))
            print(f'skipped {function.__name__}: missing {", ".join(missing)} '
                  f'(run AbaqusTriaxialConsolidation {parts} in {args.rundir})')
            continue
        function(args.rundir, args.outdir)
