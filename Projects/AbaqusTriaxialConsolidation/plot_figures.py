"""Figures 7, 8 and 9 of the article (Abaqus benchmark 1.15.2) from the files written by AbaqusTriaxialConsolidation.

Usage (run the executable first; it writes its CSV files to the current directory):

    python3 <neopz>/Projects/AbaqusTriaxialConsolidation/plot_figures.py [run directory] [-o output directory]

The run directory (default: the current directory) holds the CSV files; the figures are written as PDF and PNG
to <run directory>/figures, or to the output directory. A figure whose CSV files are missing (the executable was
run with only some of its parts) is skipped with a message.

Figures produced (same panels, curves, markers and annotations as fig_abaqus_model, fig_abaqus_states and
fig_abaqus of figs.py of the Python code):

- fig07_abaqus_model: Fig. 7, the models of the upper half of the specimen: (a) axisymmetric 2 x 4 Q8-Q4 mesh
  with the boundary conditions, the vertex (u and p_w) and mid-side (u only) nodes and the point A; (b) quarter
  of the specimen with 48 Hex20-Hex8 elements (curved lateral faces) and the boundary faces. The executable
  writes no file with the quadratic meshes (the VTK files hold the linear cells only), so the meshes are built
  here with the generators of the C++ code (mcc::CreateRectangleMesh and mcc::CreateQuarterCylinderMesh of
  Projects/Common/MCCPaperTools.h, with the dimensions of AbaqusTriaxialConsolidation::CreateGeoMesh).
- fig08_abaqus_states: Fig. 8, drained tests at a material point from p'0 = 100 kPa (R = 1.17, subcritical) and
  p'0 = 20 kPa (R = 5.83, supercritical) with p'c0 = 116.6 kPa: (a) stress paths p'-q, (b) q-eps_1 and
  (c) eps_v-eps_1 in the four-quadrant layout; 600 increments (abaqus_material_point_p0_<p'0>.csv) against the
  closed form of Appendix B.1 (abaqus_closed_form_p0_<p'0>.csv). Part "mp" of the executable.
- fig09_abaqus_results: Fig. 9, q at A against delta/H: (a) smooth platen, axisymmetric model with 150
  increments (abaqus_smooth_2x2.csv) and material point with 30 increments (abaqus_material_point_30.csv);
  (b) rough platen, axisymmetric model with 2 x 2 and 3 x 3 points (abaqus_rough_2x2.csv, abaqus_rough_3x3.csv)
  and 3D Hex20 model (abaqus_3d_rough.csv); (c) stress paths at A. Markers: the Abaqus curves of Figs. 1.15.2-2
  and -3 of the Abaqus Benchmarks Manual, digitized (reference/abaqus_1_15_2_digitalizado.json, copied unchanged
  from dados/ of the Python code). Parts "mp", "axi" and "3d" of the executable.

The figures use only the CSV files; the VTK series of the executable are not needed (the argument "novtk" of the
executable disables them).
"""
import os
import sys

import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch
from matplotlib.patches import Polygon as MPoly

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'Common'))
from mcc_figstyle import C1, C2, C3, INK, INK2, MUTED, TEXTW, arguments, plt, read_csv, reference, save  # noqa: E402

# AbaqusTriaxialConsolidation.h: upper half of the specimen (mm) and the material of the benchmark
R_MM, H_MM = 20.0, 60.0
M, PC0 = 1.0, 116.6


# =========================================================================================== meshes of Fig. 7
def rectangle_mesh(R, H, nr, nz):
    """Axisymmetric Q8 mesh of nr x nz elements on [0, R] x [0, H] (mcc::CreateRectangleMesh, nodes row by row).

    Returns the vertex coordinates, the coordinates of the mid-side nodes (one per edge) and the quadrilaterals
    (four vertex indices, counter-clockwise)."""
    V = np.array([[R * i / nr, H * j / nz] for j in range(nz + 1) for i in range(nr + 1)])
    idx = lambda i, j: j * (nr + 1) + i
    quads = [[idx(i, j), idx(i + 1, j), idx(i + 1, j + 1), idx(i, j + 1)] for j in range(nz) for i in range(nr)]
    edges = sorted({tuple(sorted((q[k], q[(k + 1) % 4]))) for q in quads for k in range(4)})
    mid = np.array([(V[a] + V[b]) / 2 for a, b in edges])
    return V, mid, quads


def quarter_disk(R, nc, nr, a_ratio=0.5):
    """Quarter of the disk of radius R (as mcc::CreateQuarterCylinderMesh): central square [0, a R]^2 with
    nc x nc cells and two outer blocks with nc cells along the arc and nr radially.

    Returns the points (x, y) and the quadrilaterals (counter-clockwise)."""
    a = a_ratio * R
    keys, pts = {}, []

    def node(x, y):
        key = (int(round(x / R * 1e9)), int(round(y / R * 1e9)))
        if key not in keys:
            keys[key] = len(pts)
            pts.append((x, y))
        return keys[key]
    S = [[node(a * i / nc, a * j / nc) for j in range(nc + 1)] for i in range(nc + 1)]
    quads = [[S[i][j], S[i + 1][j], S[i + 1][j + 1], S[i][j + 1]] for i in range(nc) for j in range(nc)]
    for block in (1, 2):
        G = []
        for j in range(nc + 1):
            if block == 1:
                ix, iy, th = a, a * j / nc, np.pi / 4 * j / nc
            else:
                ix, iy, th = a * j / nc, a, np.pi / 2 - np.pi / 4 * j / nc
            ox, oy = R * np.cos(th), R * np.sin(th)
            G.append([node(ix + k / nr * (ox - ix), iy + k / nr * (oy - iy)) for k in range(nr + 1)])
        quads += [[G[j][k], G[j][k + 1], G[j + 1][k + 1], G[j + 1][k]] for j in range(nc) for k in range(nr)]
    P = np.array(pts)
    out = []
    for q in quads:
        x, y = P[q, 0], P[q, 1]
        out.append(q if np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y) > 0 else [q[0], q[3], q[2], q[1]])
    return P, out


# local vertices of the faces of a hexahedron and boundary markers of the 3D model (numbers of the Python code,
# the last digit of the ids EX0, ELateral, EY0, EZ0 and EZH of AbaqusTriaxialConsolidation.h): 1 (x = 0) and
# 3 (y = 0) symmetry planes, 2 lateral face (cell pressure), 5 mid-plane, 6 platen
HEXFACES = [(0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]


def quarter_cylinder_faces(R, H, nc, nr, nz, a_ratio=0.5):
    """Boundary faces of the Hex20 mesh of the quarter of cylinder (mcc::CreateQuarterCylinderMesh with the marker
    of AbaqusTriaxialConsolidation::CreateGeoMesh): quarter disk extruded in nz layers, the mid-edge nodes of the
    lateral faces moved radially to r = R.

    Returns, for each boundary face, the eight nodes (four vertices, then the mid-edge nodes of the edges 0-1, 1-2,
    2-3 and 3-0) and its marker."""
    P, quads = quarter_disk(R, nc, nr, a_ratio)
    n2 = len(P)
    X = np.array([[p[0], p[1], H * k / nz] for k in range(nz + 1) for p in P])
    cells = [[k * n2 + q[i] for i in range(4)] + [(k + 1) * n2 + q[i] for i in range(4)]
             for k in range(nz) for q in quads]
    count = {}
    for c in cells:
        for f in HEXFACES:
            count.setdefault(tuple(sorted(c[i] for i in f)), []).append([c[i] for i in f])
    tol = 1e-9 * R

    def marker(F):
        if np.all(abs(F[:, 0]) < tol): return 1
        if np.all(abs(F[:, 1]) < tol): return 3
        if np.all(abs(F[:, 2]) < tol): return 5
        if np.all(abs(F[:, 2] - H) < tol): return 6
        return 2
    faces, marks = [], []
    for lst in count.values():
        if len(lst) != 1:
            continue
        f = lst[0]
        mid = []
        for i, j in ((0, 1), (1, 2), (2, 3), (3, 0)):
            m = (X[f[i]] + X[f[j]]) / 2
            if all(abs(np.hypot(*X[v, :2]) - R) < tol for v in (f[i], f[j])):   # edge on the lateral face
                m[:2] *= R / np.hypot(m[0], m[1])
            mid.append(m)
        faces.append(np.vstack([X[f], mid]))
        marks.append(marker(X[f]))
    return faces, marks


# =========================================================================================== Fig. 7 models
def fig_abaqus_model(run, out):
    """Fig. 7: (a) axisymmetric mesh with the boundary conditions and (b) quarter of the specimen in 3D (the run
    directory is not used: the meshes are built from the definition of the models)."""
    V, mid, quads = rectangle_mesh(R_MM, H_MM, 2, 4)
    fig = plt.figure(figsize=(TEXTW * 0.72, 2.9))
    ax = fig.add_axes([0.0, 0.05, 0.42, 0.9])
    for q in quads:
        ax.add_patch(MPoly(V[q], closed=True, fc='#f4f3ef', ec=INK2, lw=0.7))
    C = np.vstack([V, mid])
    ax.plot(C[:, 0], C[:, 1], 'o', ms=1.6, color=INK2)
    ax.plot(V[:, 0], V[:, 1], 's', ms=2.6, color=C1)
    # boundary conditions
    for y in np.linspace(3, 57, 7):
        ax.add_patch(FancyArrowPatch((27, y), (21.5, y), arrowstyle='-|>', mutation_scale=6, color=C2, lw=0.7))
    ax.text(27.5, 30, '$P$ = 100 kPa', rotation=90, va='center', fontsize=7, color=INK)
    for x in np.linspace(0, 20, 5):
        ax.add_patch(FancyArrowPatch((x, 66.5), (x, 61.5), arrowstyle='-|>', mutation_scale=6, color=INK, lw=0.7))
    ax.text(10, 68, 'platen: $u_z$ prescribed, $p_w = 0$', ha='center', fontsize=7)
    ax.plot([0, 20], [-1.2, -1.2], color=INK, lw=1.2)
    ax.text(10, -4.5, 'mid-plane: $u_z = 0$, impermeable', ha='center', fontsize=7)
    ax.plot([-1.2, -1.2], [0, 60], color=INK, lw=1.2)
    ax.text(-2.4, 30, 'axis: $u_r = 0$', rotation=90, va='center', ha='right', fontsize=7)
    ax.plot(5, 7.5, '*', ms=8, color=C3, mec=INK, mew=0.4)
    ax.text(6, 9, 'A', fontsize=8, fontweight='bold')
    # dimensions (upper half of the specimen)
    dim = dict(arrowstyle='<->', color=INK2, lw=0.6, mutation_scale=6, shrinkA=0, shrinkB=0)
    ax.annotate('', xy=(0, -9.5), xytext=(20, -9.5), arrowprops=dim)
    ax.text(10, -10.5, '$R$ = 20 mm', ha='center', va='top', fontsize=7)
    ax.annotate('', xy=(31.5, 0), xytext=(31.5, 60), arrowprops=dim)
    ax.text(32.3, 30, '$H$ = 60 mm (upper half)', rotation=90, ha='left', va='center', fontsize=7)
    # legend of the nodes
    ax.plot(-0.5, -17.0, 's', ms=2.6, color=C1, clip_on=False)
    ax.text(1.2, -17.0, 'vertex node: $u$ and $p_w$', fontsize=6.8, va='center')
    ax.plot(-0.5, -21.5, 'o', ms=1.6, color=INK2, clip_on=False)
    ax.text(1.2, -21.5, 'mid-side node: $u$ only', fontsize=6.8, va='center')
    ax.set_xlim(-7, 36)
    ax.set_ylim(-23, 76)
    ax.set_aspect('equal')
    ax.axis('off')
    ax.text(-7, 76, '(a)', fontsize=8.5, va='top')
    # 3D: quarter of cylinder (boundary faces with the quadratic edges)
    from mpl_toolkits.mplot3d import proj3d
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    faces, marks = quarter_cylinder_faces(R_MM, H_MM, 2, 2, 4)
    ax3 = fig.add_axes([0.42, 0.0, 0.58, 1.0], projection='3d')
    facecol = {1: '#e9e8e3', 3: '#e9e8e3', 2: '#c4e5d6', 5: '#d9d8d2', 6: '#b7d3f6'}
    t = np.linspace(0, 1, 7)
    Nq = np.c_[(1 - t) * (1 - 2 * t), 4 * t * (1 - t), t * (2 * t - 1)]
    polys, cols = [], []
    for F, mk in zip(faces, marks):
        ring = []
        for (i, j, mm) in ((0, 1, 4), (1, 2, 5), (2, 3, 6), (3, 0, 7)):
            ring += list((Nq @ F[[i, mm, j]])[:-1])
        polys.append(np.array(ring))
        cols.append(facecol[mk])
    ax3.add_collection3d(Poly3DCollection(polys, facecolors=cols, edgecolors=INK2, linewidths=0.35, alpha=1.0))
    ax3.set_xlim(0, 20)
    ax3.set_ylim(0, 20)
    ax3.set_zlim(-2, 62)
    ax3.view_init(elev=22, azim=-28)
    ax3.set_box_aspect((20, 20, 60))
    ax3.set_axis_off()
    ax3.text2D(0.05, 0.97, '(b)', transform=ax3.transAxes, fontsize=8.5, va='top')
    # labels of the faces (projected 3D point -> 2D annotation)
    fig.canvas.draw()
    P2 = lambda x, y, z: proj3d.proj_transform(x, y, z, ax3.get_proj())[:2]
    c45 = 20 * np.cos(np.pi / 4)
    lab = [((10.0, 6.0, 60.0), 'platen: $u_z$ prescribed,\n$p_w = 0$', (0.70, 0.93), C1),
           ((c45 + 1.5, c45 - 1.5, 38.0), 'lateral face:\ncell pressure\n$P$ = 100 kPa', (0.70, 0.62), C2),
           ((12.0, 0.0, 14.0), 'symmetry planes\n$x = 0$ and $y = 0$', (0.70, 0.27), INK2)]
    for (x, y, z), s, xyt, c in lab:
        ax3.annotate(s, xy=P2(x, y, z), xycoords='data', xytext=xyt, textcoords='axes fraction', fontsize=7,
                     color=INK, ha='left', va='top', arrowprops=dict(arrowstyle='-', color=c, lw=0.7))
    ax3.text2D(0.70, 0.08, 'mid-plane (hidden):\n$u_z = 0$, impermeable', transform=ax3.transAxes, fontsize=7,
               color=INK2, ha='left', va='top')
    print('Fig. 7: axisymmetric mesh %d elements, %d vertex and %d mid-side nodes; 3D mesh %d boundary faces'
          % (len(quads), len(V), len(mid), len(faces)))
    save(fig, out, 'fig07_abaqus_model')


# =========================================================================================== Fig. 8 two initial states
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
    """Fig. 8: four quadrants, p'-q (right), q-eps_1 (eps_1 to the left) and eps_v-eps_1 (eps_v downwards)."""
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
        print("Fig. 8: p'0 = %g kPa, peak q = %.3f kPa at eps_1 = %.3f (closed form %.3f); end q = %.3f, eps_v = %.4f;"
              " max |q - closed| = %.3f kPa" % (p0, mp[kq, 2], mp[kq, 0], cf[:, 2].max(), mp[-1, 2], mp[-1, 3],
                                               np.abs(mp[:, 2] - np.interp(mp[:, 0], cf[:, 0], cf[:, 2])).max()))
    save(fig, out, 'fig08_abaqus_states')


# =========================================================================================== Fig. 9 results
def history(run, name):
    """History of a finite element run (abaqus_<name>.csv): columns delta/H, p' at A, q at A, axial stress on the
    platen, largest |p_w| and eps_v (the columns of the hist arrays of gen_data.py)."""
    T = read_csv(os.path.join(run, f'abaqus_{name}.csv'))
    return np.column_stack([T[k] for k in ('delta_H', 'p_A', 'q_A', 'sigma_a_platen', 'max_pw', 'eps_v')])


def fig_abaqus(run, out):
    """Fig. 9: q at A against delta/H, smooth (a) and rough (b) platens, and the stress paths at A (c)."""
    ab = reference(HERE, 'abaqus_1_15_2_digitalizado.json')
    fig, axs = plt.subplots(1, 3, figsize=(TEXTW, 2.35))
    hs, hr, hr4 = history(run, 'smooth_2x2'), history(run, 'rough_2x2'), history(run, 'rough_3x3')
    t3h = history(run, '3d_rough')
    t3 = dict(dh=t3h[:, 0], q=t3h[:, 2])
    T = read_csv(os.path.join(run, 'abaqus_material_point_30.csv'))     # material point, 30 increments (0.02)
    mp = np.column_stack([T['delta_H'], T['p'], T['q'], T['eps_v']])
    ax = axs[0]
    ax.axhline(150.0, color=MUTED, lw=0.7, ls=(0, (1, 1.5)))
    ax.text(0.01, 153, 'critical state: $q$ = 150 kPa', fontsize=6.8, color=INK2, va='bottom')
    ax.plot(hs[:, 0], hs[:, 2], color=C1, lw=1.5, label='this work, %d increments' % (len(hs) - 1))
    ax.plot(mp[:, 0], mp[:, 2], color=INK, lw=0.8, ls=(0, (3, 2)), label='material point, %d incr.' % (len(mp) - 1))
    d = np.array(ab['lisa']['qd'])
    ax.plot(d[:, 0], d[:, 1], 'o', ms=3.8, mfc='white', mec=C2, mew=0.9, label='Abaqus')
    ax.set(xlabel='$\\delta/H$', ylabel='$q$ at A (kPa)', xlim=(0, 0.6), ylim=(0, 170))
    ax.set_title('(a) smooth platen (homogeneous)', fontsize=8)
    ax.legend(loc='lower right', fontsize=6.8)
    ax = axs[1]
    ax.plot(hr4[:, 0], hr4[:, 2], color=MUTED, lw=1.1, label='axisym., 3×3 points')
    ax.plot(hr[:, 0], hr[:, 2], color=C3, lw=1.5, label='axisym., 2×2 points')
    ax.plot(t3['dh'], t3['q'], color=INK, lw=0.8, ls=(0, (3, 2)), label='3D Hex20, 2×2×2 points')
    d = np.array(ab['rugosa']['qd'])
    ax.plot(d[:, 0], d[:, 1], 's', ms=3.6, mfc='white', mec=C2, mew=0.9, label='Abaqus')
    k = int(np.argmin(abs(hr4[:, 0] - 0.56)))
    ax.annotate('3×3: volumetric\nlocking', xy=(hr4[k, 0], hr4[k, 2]), xytext=(0.42, 112), fontsize=6.8,
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
    pm = np.linspace(95, 165, 2)
    ax.plot(pm, M * pm, color=MUTED, lw=0.8)
    ax.text(127, 133, "CSL: $q = Mp'$", fontsize=6.8, color=INK2, ha='right', va='bottom')
    ax.set(xlabel="$p'$ at A (kPa)", ylabel='$q$ at A (kPa)', xlim=(95, 165), ylim=(0, 170))
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
    print('Fig. 9: q_A at delta/H = 0.6: smooth %.3f, material point (30 incr.) %.3f, rough 2x2 %.3f, rough 3x3 %.3f'
          ' (max %.3f at %.3f), 3D rough %.3f kPa' % (hs[-1, 2], mp[-1, 2], hr[-1, 2], hr4[-1, 2], hr4[:, 2].max(),
                                                     hr4[np.argmax(hr4[:, 2]), 0], t3h[-1, 2]))
    save(fig, out, 'fig09_abaqus_results')


# parts of the executable (command line arguments, in the order of the executable)
PARTS = ('mp', 'axi', 'states', '3d')
# (function, {CSV file: part of the executable that writes it})
FIGURES = (
    (fig_abaqus_model, {}),
    (fig_abaqus_states, {'abaqus_material_point_p0_100.csv': 'mp', 'abaqus_closed_form_p0_100.csv': 'mp',
                         'abaqus_material_point_p0_20.csv': 'mp', 'abaqus_closed_form_p0_20.csv': 'mp'}),
    (fig_abaqus, {'abaqus_smooth_2x2.csv': 'axi', 'abaqus_rough_2x2.csv': 'axi', 'abaqus_rough_3x3.csv': 'axi',
                  'abaqus_3d_rough.csv': '3d', 'abaqus_material_point_30.csv': 'mp'}),
)

if __name__ == '__main__':
    args = arguments('Figs. 7 to 9 of the article (Abaqus benchmark 1.15.2) from the CSV files of '
                     'AbaqusTriaxialConsolidation.')
    for function, files in FIGURES:
        missing = [f for f in files if not os.path.exists(os.path.join(args.rundir, f))]
        if missing:
            parts = ' '.join(p for p in PARTS if any(files[f] == p for f in missing))
            print(f'skipped {function.__name__}: missing {", ".join(missing)} '
                  f'(run AbaqusTriaxialConsolidation {parts} in {args.rundir})')
            continue
        function(args.rundir, args.outdir)
