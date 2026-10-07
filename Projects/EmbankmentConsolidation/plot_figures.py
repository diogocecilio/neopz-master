"""Figures 11, 12 and 13 of the article (embankment on a Cam-Clay foundation, Sect. 6.6) from the files written by
EmbankmentConsolidation (3D model: 20 x 10 x 1 Hex20-Hex8 elements, the 1 m slice of the FLAC3D model).

Usage (run the executable first; it writes its CSV files to the current directory):

    python3 <neopz>/Projects/EmbankmentConsolidation/plot_figures.py [run directory] [-o output directory]

The run directory (default: the current directory) holds the CSV files; the figures are written as PDF and PNG
to <run directory>/figures, or to the output directory. A figure whose CSV files are missing is skipped with a
message. This script draws the three states of the article; the fields at every converged state (36 states) are
in the VTK series written by the executable in vtk/embankment and vtk/embankment_elastic, to be opened in
ParaView (README.md, section "Viewing the solution in ParaView").

Figures produced:

- fig12_embankment_model: Fig. 12, the model in an oblique (cavalier) projection with x to the right, y up and z
  towards the reader: the mesh of the slab (embankment_mesh_*.csv, written by mcc::WriteMeshCSV: the faces z = 1 m,
  y = 10 m and x = 20 m are visible), the strip load, the drained top, the boundary conditions (u_x = 0 at x = 0
  and x = 20 m, u_z = 0 on the faces z = 0 and z = 1 m, fixed and impermeable base), the settlement points at
  x = 0, 2, 4, 6 m (vertices of the face z = 0) and the elements of the zones pp1 and pp2
  (embankment_monitor.csv); on the right, one Hex20-Hex8 element with its 20 displacement nodes and its 8 pore
  pressure nodes.
- fig13_embankment_history: Fig. 13, histories from the end of the undrained loading (t = 0) to t = 1e8 s
  (embankment_history.csv): (a) settlements of the top at x = 0, 2, 4, 6 m, (b) pore pressures in the zones pp1
  and pp2, with pp2 against log t in the inset (Mandel-Cryer effect). Markers: the FLAC3D histories, digitized
  (reference/flac_historicos_digitalizados.json, the file dados/ of the Python code), with the FLAC3D values at
  the end of the loading (Table 8) at t = 0. The numbers of the discussion read from the FLAC3D histories (share of
  the excess pore pressure of pp2 dissipated and of the consolidation settlement at x = 0 developed at t = 2.5e5 s
  and 1e6 s), with those of this work, are printed and written to fig13_flac3d_numbers.csv in the output directory.
- fig14_embankment_fields: Fig. 14, fields on the face z = 0 (embankment_nodal_<state>.csv, vertices of the face,
  contours on the two triangles of each quadrilateral, as in the Python code; the fields do not depend on z):
  (a) excess pore pressure at the end of the undrained loading, (b) at t = 1e6 s, (c) plastic integration points
  at t = 1e8 s in the layer of points nearest to the face z = 0 (embankment_gauss_t1e8.csv: type 1 subcritical,
  2 supercritical; the three layers of points are identical), (d) settlement at t = 1e8 s.
"""
import csv
import os
import sys

import matplotlib.tri as mtri
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch
from matplotlib.patches import Polygon as MPoly

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'Common'))
from mcc_figstyle import (C1, C2, C3, C4, ENGGREEN, GRID, INK, INK2, MUTED, SEQ, TEXTW, arguments, plt,  # noqa: E402
                          read_csv, reference, save)

# FLAC3D values at the end of the undrained loading (t = 0; Table 8 of the article), the first point of the
# FLAC3D histories (the digitized curves start at t = 2.5e5 s)
FLAC_T0 = dict(uz_x0=0.140, uz_x2=0.135, uz_x4=0.055, uz_x6=-0.0416, pp1=18.1, pp2=62.4)
X_SETTLEMENT = (0, 2, 4, 6)  # settlement points of the top (m), columns s_x0 ... s_x6 of embankment_history.csv
# material ids of the boundary faces (EmbankmentConsolidation.h)
EBOTTOM, ERIGHT, ELEFT, ELOAD, EZ0, EZ1, EPTOP = -1, -2, -4, -5, -6, -7, -13


# =========================================================================================== mesh and data
def read_table(path):
    """CSV file with a header and text or numeric columns (embankment_monitor.csv): list of dicts."""
    if not os.path.exists(path):
        sys.exit(f'file not found: {path} (run the executable in this directory first)')
    with open(path, newline='') as f:
        return list(csv.DictReader(f))


def face_z0(nodal):
    """Rows of a nodal CSV file (vertices of the slab) on the face z = 0, as a dict of arrays."""
    sel = np.abs(nodal['z']) < 1e-9
    return {k: v[sel] for k, v in nodal.items()}


def vertex_mesh(nodal):
    """Quadrilaterals of the face from its vertices (structured grid of EmbankmentConsolidation::CreateGeoMesh).

    Returns the quadrilaterals as four vertex indices (rows of the arrays), counter-clockwise from the lower
    left corner, as the first four nodes of the elements of the Python code."""
    x, y = np.round(nodal['x'], 9), np.round(nodal['y'], 9)
    xs, ys = np.unique(x), np.unique(y)
    index = {(a, b): k for k, (a, b) in enumerate(zip(x, y))}
    if len(index) != len(x) or len(x) != len(xs) * len(ys):
        sys.exit('the vertices of the face z = 0 do not form a structured grid')
    return np.array([[index[xs[i], ys[j]], index[xs[i + 1], ys[j]], index[xs[i + 1], ys[j + 1]],
                      index[xs[i], ys[j + 1]]] for j in range(len(ys) - 1) for i in range(len(xs) - 1)])


def triangulation(nodal, quads):
    """Two triangles per quadrilateral (diagonal from the first to the third vertex), as quad_tri of figs.py."""
    tris = np.stack([quads[:, [0, 1, 2]], quads[:, [0, 2, 3]]], axis=1).reshape(-1, 3)
    return mtri.Triangulation(nodal['x'], nodal['y'], tris)


def consolidation_history(run):
    """Rows of embankment_history.csv from the end of the undrained loading (the last state with t = 0) on, as an
    array with the columns t, s(x = 0), s(x = 2), s(x = 4), s(x = 6), pp1, pp2 (hist[nU - 1:] of figs.py)."""
    h = read_csv(os.path.join(run, 'embankment_history.csv'))
    iu = int(np.nonzero(h['t'] <= 0.0)[0][-1])
    cols = ['t'] + [f's_x{x}' for x in X_SETTLEMENT] + ['pp1', 'pp2']
    return np.column_stack([h[c] for c in cols])[iu:]


def flac_history(flac, key, t0=True):
    """Digitized FLAC3D history sorted by t; with t0, the value at the end of the loading is prepended at t = 0."""
    t, v = np.array(flac[key][0]), np.array(flac[key][1])
    o = np.argsort(t)
    t, v = t[o], v[o]
    return (np.r_[0.0, t], np.r_[FLAC_T0[key], v]) if t0 else (t, v)


def flac_numbers(flac, hc, out):
    """Numbers of the discussion of Sect. 6.6 read from the digitized FLAC3D histories (linear interpolation in t),
    with the values of this work at the same times from the history hc: share of the excess pore pressure of pp2
    dissipated and share of the consolidation settlement at x = 0 developed at the first FLAC3D point
    (t = 2.5e5 s) and at t = 1e6 s. Printed and written to <out>/fig13_flac3d_numbers.csv."""
    hyd2 = 25.0  # hydrostatic pore pressure at the centre of pp2 (y = 7.5 m)
    tp, vp = flac_history(flac, 'pp2', t0=False)
    ts, vs = flac_history(flac, 'uz_x0', t0=False)
    sfin = vs[-1]
    rows = []
    for t in (tp[0], 1.0e6):
        p = float(np.interp(t, tp, vp))
        sx0 = float(np.interp(t, ts, vs))
        pw = float(np.interp(np.log(t), np.log(hc[1:, 0]), hc[1:, 6]))
        sw = float(np.interp(np.log(t), np.log(hc[1:, 0]), hc[1:, 1]))
        rows.append((t, p, 100 * (FLAC_T0['pp2'] - p) / (FLAC_T0['pp2'] - hyd2), sx0,
                     100 * (sx0 - FLAC_T0['uz_x0']) / (sfin - FLAC_T0['uz_x0']), pw,
                     100 * (hc[0, 6] - pw) / (hc[0, 6] - hyd2), sw, 100 * (sw - hc[0, 1]) / (hc[-1, 1] - hc[0, 1])))
        print('FLAC3D at t = %.3g s: pp2 %.2f kPa (%.0f%% of the excess dissipated), settlement x = 0 %.4f m (%.0f%% of '
              'the consolidation settlement); this work: pp2 %.2f kPa (%.1f%%), settlement %.4f m (%.1f%%)' % rows[-1])
    with open(os.path.join(out, 'fig13_flac3d_numbers.csv'), 'w') as f:
        f.write('t,flac_pp2,flac_pp2_dissipated_pct,flac_s_x0,flac_s_x0_consolidation_pct,'
                'this_pp2,this_pp2_dissipated_pct,this_s_x0,this_s_x0_consolidation_pct\n')
        for r in rows:
            f.write(','.join('%.6g' % v for v in r) + '\n')


# =========================================================================================== Fig. 12
# oblique (cavalier) projection: x to the right, y up, z towards the reader drawn down and to the left at 45 degrees
# with its true length, so that the visible faces are z = 1 m (front), y = 10 m (top) and x = 20 m (right)
OBLIQUE = np.array([-np.cos(np.pi / 4), -np.sin(np.pi / 4)])


def project(X):
    """Oblique projection of points (..., 3) to the plane of the drawing (..., 2)."""
    X = np.asarray(X, dtype=float)
    return X[..., :2] + X[..., 2:3] * OBLIQUE


def fig_embankment_model(run, out):
    """Fig. 12: the 3D model, the boundary conditions, the monitoring points and the Hex20-Hex8 element."""
    nodes = read_csv(os.path.join(run, 'embankment_mesh_nodes.csv'))
    faces = read_csv(os.path.join(run, 'embankment_mesh_faces.csv'))
    elements = read_csv(os.path.join(run, 'embankment_mesh_elements.csv'))
    edges = read_csv(os.path.join(run, 'embankment_mesh_edges.csv'))
    monitor = read_table(os.path.join(run, 'embankment_monitor.csv'))
    X = np.column_stack([nodes['x'], nodes['y'], nodes['z']])
    W, H, T = X[:, 0].max(), X[:, 1].max(), X[:, 2].max()
    fmat = faces['matid'].astype(int)
    fnodes = np.column_stack([faces[f'n{k}'] for k in range(4)]).astype(int)

    fig = plt.figure(figsize=(TEXTW, 3.05))
    ax = fig.add_axes([0.0, 0.0, 0.80, 1.0])
    axe = fig.add_axes([0.80, 0.20, 0.20, 0.62])

    def faces_of(matid):
        return fnodes[fmat == matid]

    # visible faces: z = T (front), y = H (top), x = W (right); the loaded strip of the top in green
    loaded = {tuple(sorted(f)) for f in faces_of(ELOAD)}
    for f in faces_of(EZ1):
        ax.add_patch(MPoly(project(X[f]), closed=True, fc='#f4f3ef', ec=INK2, lw=0.4, zorder=1))
    for f in faces_of(EPTOP):
        fc = '#cfe8dc' if tuple(sorted(f)) in loaded else '#e3eefa'
        ax.add_patch(MPoly(project(X[f]), closed=True, fc=fc, ec=INK2, lw=0.4, zorder=1))
    for f in faces_of(ERIGHT):
        ax.add_patch(MPoly(project(X[f]), closed=True, fc='#e8e7e1', ec=INK2, lw=0.4, zorder=1))
    # outline of the slab
    corners = np.array([[0, 0, T], [W, 0, T], [W, 0, 0], [W, H, 0], [0, H, 0], [0, H, T]])
    ax.add_patch(MPoly(project(corners), closed=True, fc='none', ec=INK, lw=0.9, zorder=3))
    ax.plot(*project(np.array([[0, H, T], [W, H, T], [W, 0, T]])).T, color=INK, lw=0.6, zorder=3)
    ax.plot(*project(np.array([[W, H, T], [W, H, 0]])).T, color=INK, lw=0.6, zorder=3)
    # hidden edges of the face z = 0 and of the face x = 0
    ax.plot(*project(np.array([[0, H, 0], [0, 0, 0], [W, 0, 0]])).T, color=MUTED, lw=0.5, ls=(0, (2, 2)), zorder=2)
    ax.plot(*project(np.array([[0, 0, 0], [0, 0, T]])).T, color=MUTED, lw=0.5, ls=(0, (2, 2)), zorder=2)

    # load on the strip 0 <= x <= 4 m (arrows onto the middle of the top face)
    xl = max(X[f, 0].max() for f in faces_of(ELOAD))
    for x in np.linspace(0.25, xl - 0.25, 9):
        p = project([x, H, T / 2])
        ax.add_patch(FancyArrowPatch(p + [0, 1.25], p + [0, 0.05], arrowstyle='-|>', mutation_scale=5, color=C2,
                                     lw=0.7, zorder=6))
    pa, pb = project([0.25, H, T / 2]) + [0, 1.25], project([xl - 0.25, H, T / 2]) + [0, 1.25]
    ax.plot([pa[0], pb[0]], [pa[1], pb[1]], color=C2, lw=0.9, zorder=6)
    ax.text(pb[0] + 0.5, pb[1], '$q$ = 50 kPa', ha='left', va='center', fontsize=7.3)
    # drained top: water table symbol and label
    pw = project([W - 3.0, H, T / 2])
    ax.add_patch(MPoly([pw + [-0.3, 0.85], pw + [0.3, 0.85], pw + [0, 0.15]], closed=True, fc='white', ec=C1,
                       lw=0.8, zorder=6))
    ax.text(pw[0] - 0.6, pw[1] + 0.8, 'drained top, water table ($p_w = 0$)', ha='right', va='center', fontsize=7.3,
            color=INK)

    # settlement points (vertices of the face z = 0) and zones pp1, pp2 (front face of their elements)
    for m in monitor:
        if m['name'].startswith('s_x'):
            p = project([float(m['x']), float(m['y']), float(m['z'])])
            ax.plot(*p, 'v', ms=4.6, color=C3, mec=INK, mew=0.4, zorder=7)
    p6 = project([6.0, H, 0.0])
    ax.annotate('settlement points ($z$ = 0)', xy=p6 + [0.1, -0.1], xytext=(p6[0] + 1.3, p6[1] - 1.75), fontsize=7.0,
                color=INK, ha='left', va='center', arrowprops=dict(arrowstyle='-', color=INK2, lw=0.5), zorder=7,
                bbox=dict(boxstyle='round,pad=0.2', fc='white', ec='none', alpha=0.9))
    for m in monitor:
        if m['name'] in ('pp1', 'pp2'):
            x0, x1, y0, y1 = (float(m[k]) for k in ('xmin', 'xmax', 'ymin', 'ymax'))
            quad = project(np.array([[x0, y0, T], [x1, y0, T], [x1, y1, T], [x0, y1, T]]))
            ax.add_patch(MPoly(quad, closed=True, fc=C4, ec=INK, lw=0.5, alpha=0.85, zorder=4))
            ax.text(quad[1, 0] + 0.15, 0.5 * (quad[1, 1] + quad[2, 1]), m['name'], fontsize=7.3, va='center', zorder=6)

    # material
    ax.text(12.0, 4.4, "clay layer (MCC), saturated\n"
                       "$p'_{c0}$ = 160 kPa (uniform)\n"
                       "$\\sigma'_v = 13d$, $\\sigma'_h = 6.1d$ kPa\n"
                       "($d$: depth in m)",
            ha='center', va='center', fontsize=7.0, color=INK, linespacing=1.35,
            bbox=dict(boxstyle='round,pad=0.35', fc='white', ec=GRID, lw=0.6, alpha=0.95), zorder=7)
    ax.text(12.0, 1.35, '$u_z = 0$ on the faces $z$ = 0 and $z$ = 1 m (plane strain)', ha='center', va='center',
            fontsize=7.0, color=INK, bbox=dict(boxstyle='round,pad=0.25', fc='white', ec='none', alpha=0.9), zorder=7)

    # supports: u_x = 0 at x = 0 (rollers left of the front face) and x = W (right of the right face)
    # (rollers between the slab and the support line on both sides; the right face x = W ends at x = W in the drawing)
    for y in np.linspace(0.6, H - 0.6, 8):
        pl = project([0, y, T])
        ax.plot(pl[0] - 0.3, pl[1], 'o', ms=3.2, mfc='white', mec=INK, mew=0.7, zorder=5)
        pr = np.array([W + 0.3, project([W, y, T / 2])[1]])
        ax.plot(pr[0], pr[1], 'o', ms=3.2, mfc='white', mec=INK, mew=0.7, zorder=5)
    pl0, pl1 = project([0, 0, T]), project([0, H, T])
    ax.plot([pl0[0] - 0.55] * 2, [pl0[1], pl1[1]], color=INK, lw=0.9)
    pr0, pr1 = project([W, 0, 0]), project([W, H, 0])
    ax.plot([pr0[0] + 0.55] * 2, [pr0[1] - 0.35, pr1[1] - 0.35], color=INK, lw=0.9)
    ax.text(pl0[0] - 0.95, 0.5 * (pl0[1] + pl1[1]), 'symmetry: $u_x = 0$', rotation=90, ha='right', va='center',
            fontsize=7.3)
    ax.text(pr0[0] + 1.2, 0.5 * (pr0[1] + pr1[1]), '$u_x = 0$', rotation=90, ha='left', va='center', fontsize=7.3)
    # fixed, impermeable base (hatch under the front edge)
    b0, b1 = project([0, 0, T]), project([W, 0, T])
    ax.plot([b0[0], b1[0]], [b0[1] - 0.22, b1[1] - 0.22], color=INK, lw=1.4)
    for x in np.linspace(b0[0], b1[0], 41):
        ax.plot([x, x - 0.3], [b0[1] - 0.22, b0[1] - 0.57], color=INK, lw=0.55)
    ax.text(0.5 * (b0[0] + b1[0]), b0[1] - 0.75, 'fixed, impermeable base', ha='center', va='top', fontsize=7.3)

    # dimensions
    def dim(p, q, text, offset, **kw):
        p, q = np.asarray(p, float) + offset, np.asarray(q, float) + offset
        ax.annotate('', xy=p, xytext=q, arrowprops=dict(arrowstyle='<->', color=INK2, lw=0.6, mutation_scale=6))
        ax.text(*(0.5 * (p + q)), text, fontsize=7.3, **kw)

    dim(project([0, 0, T]), project([W, 0, T]), '20 m', [0, -2.15], ha='center', va='top')
    dim(project([W, 0, 0]), project([W, H, 0]), '10 m', [2.75, 0], ha='left', va='center', rotation=90)
    dim(project([0, H, 0]), project([xl, H, 0]), '4 m', [0, 2.15], ha='center', va='bottom')
    # thickness of the slab: the edge from (W, H, 0) to (W, H, T), dimensioned above the top face (offset normal to the
    # edge in the drawing, with extension lines)
    pz0, pz1 = project([W, H, 0]), project([W, H, T])
    nz = np.array([-1.0, 1.0]) / np.sqrt(2.0)
    for p in (pz0, pz1):
        ax.plot(*np.column_stack([p + 0.15 * nz, p + 1.25 * nz]), color=INK2, lw=0.4)
    dim(pz0, pz1, '', 1.1 * nz)
    pm = 0.5 * (pz0 + pz1) + 1.45 * nz
    ax.text(pm[0], pm[1], '1 m', fontsize=7.3, ha='right', va='bottom')
    # axes (outside the slab, upper left)
    o = np.array([-2.7, H + 0.2])
    for v, lab in (([1.3, 0, 0], '$x$'), ([0, 1.3, 0], '$y$'), ([0, 0, 1.3], '$z$')):
        d = project(np.array(v, dtype=float))
        ax.add_patch(FancyArrowPatch(o, o + d, arrowstyle='-|>', mutation_scale=5, color=INK, lw=0.6))
        ax.text(*(o + d * (1 + 0.28 / np.linalg.norm(d))), lab, ha='center', va='center', fontsize=7.3)
    ax.set_xlim(-4.0, W + 3.6)
    ax.set_ylim(-3.2, H + 2.9)
    ax.set_aspect('equal')
    ax.axis('off')

    # one Hex20-Hex8 element (unit cube; nodes from the edges file: vertices and mid-edge nodes)
    e0 = elements['element'].astype(int)
    first = np.column_stack([elements[f'n{k}'] for k in range(8)]).astype(int)[np.argsort(e0)[0]]
    V = X[first] - X[first].min(axis=0)
    cube = [(0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4), (0, 4), (1, 5), (2, 6), (3, 7)]
    # the three edges at the hidden vertex, the corner (xmin, ymin, zmin) of the element, are dashed
    h0 = int(np.argmin(V.sum(axis=1)))
    PV = project(V)
    for a, b in cube:
        hid = h0 in (a, b)
        axe.plot(*PV[[a, b]].T, color=MUTED if hid else INK2, lw=0.7, ls=(0, (2, 2)) if hid else '-')
    mids = np.array([0.5 * (PV[a] + PV[b]) for a, b in cube])
    axe.plot(mids[:, 0], mids[:, 1], 'o', ms=3.4, mfc='white', mec=C1, mew=0.9, zorder=4)
    axe.plot(PV[:, 0], PV[:, 1], 'o', ms=4.4, mfc=C1, mec=INK, mew=0.5, zorder=5)
    axe.plot(PV[:, 0], PV[:, 1], 's', ms=8.0, mfc='none', mec=C2, mew=0.9, zorder=5)
    axe.set_aspect('equal')
    axe.set_xlim(-1.05, 1.25)
    axe.set_ylim(-2.95, 1.25)
    axe.axis('off')
    axe.text(0.1, 1.2, 'Hex20–Hex8 element', ha='center', va='bottom', fontsize=7.5)
    hh = [Line2D([], [], ls='', marker='o', ms=4.4, mfc=C1, mec=INK, mew=0.5, label='vertex node: $u$'),
          Line2D([], [], ls='', marker='o', ms=3.4, mfc='white', mec=C1, mew=0.9, label='mid-edge node: $u$'),
          Line2D([], [], ls='', marker='s', ms=8.0, mfc='none', mec=C2, mew=0.9, label='pore pressure node: $p_w$')]
    axe.legend(handles=hh, loc='upper center', bbox_to_anchor=(0.5, 0.40), fontsize=6.8, handletextpad=0.4,
               borderaxespad=0.0, labelspacing=0.7)
    axe.text(0.1, -2.98, '3 × 3 × 3 Gauss points', ha='center', va='bottom', fontsize=6.8, color=INK2)
    nel = len(elements['element'])
    nmid = int(np.sum(edges['nmid'] >= 0)) if 'nmid' in edges else 0
    print('Fig. 12: %d hexahedra, %d vertices, %d edges (%d with a geometric mid-edge node), %d boundary faces; '
          'slab %g x %g x %g m' % (nel, len(X), len(edges['n0']), nmid, len(fmat), W, H, T))
    save(fig, out, 'fig12_embankment_model')


# =========================================================================================== Fig. 13
def fig_embankment_history(run, out):
    """Fig. 13: settlements and pore pressures against t, with the FLAC3D histories (fig_aterro_hist of figs.py)."""
    hc = consolidation_history(run)
    flac = reference(HERE, 'flac_historicos_digitalizados.json')
    fig, axs = plt.subplots(1, 2, figsize=(TEXTW, 2.7))
    ax = axs[0]
    cols = [C1, C2, C3, C4]
    hh = []
    for k, x in enumerate(X_SETTLEMENT):
        ax.plot(hc[:, 0] / 1e6, hc[:, 1 + k], color=cols[k], lw=1.4)
        t, v = flac_history(flac, f'uz_x{x}')
        sel = np.linspace(0, len(t) - 1, 14).astype(int)
        ax.plot(t[sel] / 1e6, v[sel], 'o', ms=3.3, mfc='white', mec=cols[k], mew=0.8)
        hh.append(Line2D([], [], color=cols[k], lw=1.4, label=f'$x$ = {x} m'))
    hh += [Line2D([], [], color=INK2, lw=1.4, label='this work'),
           Line2D([], [], ls='', marker='o', ms=3.3, mfc='white', mec=INK2, label='FLAC3D')]
    ax.set(xlabel='$t$ (10$^6$ s)', ylabel='settlement (m)', xlim=(0, 100), ylim=(-0.06, 0.40))
    ax.set_yticks([0, 0.1, 0.2, 0.3])
    ax.set_title('(a) surface settlement', fontsize=8)
    ax.legend(handles=hh, loc='upper center', ncol=3, fontsize=6.9, columnspacing=1.0, handlelength=1.6)
    # final values at x = 0 (the difference discussed in the text); the FLAC3D label is placed above the value it
    # prints (last digitized point, rounded to mm), as in figs.py
    s_flac = round(float(flac_history(flac, 'uz_x0', t0=False)[1][-1]), 3)
    ax.text(99, hc[-1, 1] + 0.012, f'this work, $x$ = 0: {hc[-1, 1]:.3f} m', ha='right', va='bottom',
            fontsize=6.8, color=INK)
    ax.text(99, s_flac + 0.010, f'FLAC3D, $x$ = 0: {s_flac:.3f} m', ha='right', va='bottom', fontsize=6.8,
            color=INK)
    ax = axs[1]
    for k, key in enumerate(('pp1', 'pp2')):
        ax.plot(hc[:, 0] / 1e6, hc[:, 5 + k], color=cols[k], lw=1.4)
        t, v = flac_history(flac, key)
        sel = np.linspace(0, len(t) - 1, 14).astype(int)
        ax.plot(t[sel] / 1e6, v[sel], 'o', ms=3.3, mfc='white', mec=cols[k], mew=0.8)
        ax.text(101, hc[-1, 5 + k] + (1.4 if k == 1 else -1.4), key, fontsize=7.5, va='center', color=INK)
    ax.set(xlabel='$t$ (10$^6$ s)', ylabel='pore pressure (kPa)', xlim=(0, 112), ylim=(0, 70))
    ax.set_xticks([0, 20, 40, 60, 80, 100])
    ax.set_title('(b) pore pressure in zones pp1 and pp2', fontsize=8)
    ax.legend(handles=hh[-2:], loc='center', bbox_to_anchor=(0.62, 0.22), fontsize=7)
    # detail in log scale: Mandel-Cryer effect in pp2 and early dissipation of FLAC3D
    axin = ax.inset_axes([0.40, 0.53, 0.57, 0.37])
    axin.semilogx(hc[1:, 0], hc[1:, 6], color=C2, lw=1.3)
    t, v = flac_history(flac, 'pp2', t0=False)
    tk = np.unique(np.searchsorted(t, np.logspace(np.log10(t[0]), 8, 9)).clip(0, len(t) - 1))
    axin.semilogx(t[tk], v[tk], 'o', ms=3.0, mfc='white', mec=C2, mew=0.8)
    km = int(np.argmax(hc[:, 6]))
    axin.annotate(f'{hc[0, 6]:.1f} → {hc[km, 6]:.1f} kPa\n(Mandel–Cryer)', xy=(hc[km, 0], hc[km, 6]),
                  xytext=(2e2, 46), fontsize=6.5, color=INK, ha='left', va='top',
                  arrowprops=dict(arrowstyle='-', color=MUTED, lw=0.6))
    axin.set_xlim(1e2, 1e8)
    axin.set_ylim(20, 65)
    axin.set_xticks([1e2, 1e4, 1e6, 1e8])
    axin.set_yticks([20, 40, 60])
    axin.tick_params(labelsize=6.3, length=2, pad=1.5)
    axin.set_title('pp2 against log $t$ (s)', fontsize=6.8, pad=2)
    axin.grid(True, which='major', lw=0.4)
    fig.tight_layout(w_pad=1.2)
    flac_numbers(flac, hc, out)
    print('Fig. 13: t = 1e8 s: settlements %s m, pp1 %.4f, pp2 %.4f kPa; pp2 %.3f -> %.3f kPa at t = %.3g s'
          % (np.array2string(hc[-1, 1:5], precision=6), hc[-1, 5], hc[-1, 6], hc[0, 6], hc[km, 6], hc[km, 0]))
    save(fig, out, 'fig13_embankment_history')


# =========================================================================================== Fig. 14
def fig_embankment_fields(run, out):
    """Fig. 14: excess pore pressure at the end of the loading and at t = 1e6 s, plastic integration points and
    settlement at t = 1e8 s, on the face z = 0 (fig_aterro_fields of figs.py)."""
    nodal = {tag: face_z0(read_csv(os.path.join(run, f'embankment_nodal_{tag}.csv')))
             for tag in ('undrained', 't1e6', 't1e8')}
    gauss = read_csv(os.path.join(run, 'embankment_gauss_t1e8.csv'))
    quads = vertex_mesh(nodal['t1e8'])
    cmap = LinearSegmentedColormap.from_list('seq', ['#ffffff'] + SEQ)
    fig = plt.figure(figsize=(TEXTW, 3.6))
    gs = fig.add_gridspec(2, 3, width_ratios=(1, 1, 0.03), wspace=0.16, hspace=0.42, left=0.06, right=0.94,
                          bottom=0.09, top=0.94)
    axa, axb, cax1 = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[0, 2])
    axd, axc, cax2 = fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1]), fig.add_subplot(gs[1, 2])
    lev = np.linspace(0, 60, 13)
    for ax, tag, title in ((axa, 'undrained', '(a) excess pore pressure, end of undrained loading'),
                           (axb, 't1e6', '(b) excess pore pressure, $t$ = 10$^6$ s')):
        n = nodal[tag]
        tri = triangulation(n, vertex_mesh(n))
        pex = n['p_excess']
        cs = ax.tricontourf(tri, pex, levels=lev, cmap=cmap, extend='max')
        ax.tricontour(tri, pex, levels=lev[1::2], colors=INK2, linewidths=0.35)
        ax.set_title(title, fontsize=8)
        ax.text(19.6, 0.6, f'max {pex.max():.1f} kPa', ha='right', va='bottom', fontsize=6.8, color=INK)
        k = int(np.argmax(pex))
        print('Fig. 14: %s: largest excess pore pressure on the face z = 0: %.2f kPa at (%g, %g)'
              % (tag, pex[k], n['x'][k], n['y'][k]))
    cb = fig.colorbar(cs, cax=cax1)
    cb.set_label('kPa', fontsize=7.5)
    cb.ax.tick_params(labelsize=7)
    # settlement at t = 1e8 s; diverging: blue (heave) | white at zero | enggreen ramp (settlement)
    n = nodal['t1e8']
    tri = triangulation(n, quads)
    s = -n['uy']
    top = np.abs(n['y'] - n['y'].max()) < 1e-9
    heave = max(0.0, -s[top].min())
    levu = np.r_[-0.01, np.linspace(0.0, 0.28, 15)]
    cmapu = LinearSegmentedColormap.from_list('div', [(0.0, '#256abf'), (0.5, '#ffffff'), (0.62, '#d2ebe0'),
                                                      (0.75, '#93cdb3'), (0.88, '#3f9e7a'), (1.0, ENGGREEN)])
    cs2 = axc.tricontourf(tri, s, levels=levu, cmap=cmapu, norm=TwoSlopeNorm(vmin=-0.02, vcenter=0.0, vmax=0.28))
    axc.tricontour(tri, s, levels=levu[2::2], colors=INK2, linewidths=0.35)
    axc.annotate(f'slight heave\n(up to {1e3 * heave:.0f} mm)', xy=(15.0, 8.5), xytext=(12.2, 5.6), fontsize=6.8,
                 color=INK, ha='left', va='top', arrowprops=dict(arrowstyle='-', color=INK2, lw=0.6))
    axc.set_title('(d) settlement, $t$ = 10$^8$ s', fontsize=8)
    cb2 = fig.colorbar(cs2, cax=cax2)
    cb2.set_label('settlement (m)', fontsize=7.5)
    cb2.ax.tick_params(labelsize=7)
    cb2.set_ticks([0, 0.1, 0.2])
    # plastic integration points at t = 1e8 s: the layer of points nearest to the face z = 0
    C = np.column_stack([n['x'], n['y']])
    for q in quads:
        axd.add_patch(MPoly(C[q], closed=True, fc='white', ec=GRID, lw=0.4))
    zl = np.round(gauss['z'], 9)
    layers = np.unique(zl)
    first = zl == layers[0]
    typ = np.rint(gauss['type']).astype(int)
    sel1, sel2 = (typ == 1) & first, (typ == 2) & first
    axd.plot(gauss['x'][sel1], gauss['y'][sel1], 'o', ms=2.6, mfc=C1, mec='white', mew=0.3,
             label=f'subcritical ({sel1.sum()})')
    axd.plot(gauss['x'][sel2], gauss['y'][sel2], '^', ms=3.0, mfc=C2, mec='white', mew=0.3,
             label=f'supercritical ({sel2.sum()})')
    axd.set_title('(c) plastic Gauss points, $t$ = 10$^8$ s', fontsize=8)
    axd.legend(loc='upper right', fontsize=6.8, frameon=True, framealpha=0.9, edgecolor='none', handletextpad=0.3,
               title=f'layer $z$ = {layers[0]:.3f} m', title_fontsize=6.6)
    for ax in (axa, axb, axc, axd):
        ax.plot([0, 4], [10, 10], color=C2, lw=3.2, solid_capstyle='butt', clip_on=False, zorder=6)  # loaded strip
        ax.set_aspect('equal')
        ax.set_xlim(0, 20)
        ax.set_ylim(0, 10)
        ax.set_xticks([0, 5, 10, 15, 20])
        ax.set_yticks([0, 5, 10])
        ax.grid(False)
        ax.tick_params(labelsize=7)
    for ax in (axc, axd):
        ax.set_xlabel('$x$ (m)', fontsize=8)
    for ax in (axa, axd):
        ax.set_ylabel('$y$ (m)', fontsize=8)
    # colour bars with the height of the axes (equal aspect)
    fig.canvas.draw()
    for ax, cax in ((axb, cax1), (axc, cax2)):
        pa, pc = ax.get_position(), cax.get_position()
        cax.set_position([pc.x0, pa.y0, pc.width, pa.height])
    per_layer = ', '.join('z = %.4f: %d sub, %d super' % (z, np.sum((typ == 1) & (zl == z)), np.sum((typ == 2) & (zl == z)))
                          for z in layers)
    print('Fig. 14: t = 1e8 s: plastic points (%s; total %d, %d of %d); settlement from %.4f to %.4f m, heave of the '
          'top up to %.4f mm' % (per_layer, np.sum(typ == 1), np.sum(typ == 2), len(typ), s.min(), s.max(), 1e3 * heave))
    save(fig, out, 'fig14_embankment_fields')


# (function, CSV files)
FIGURES = (
    (fig_embankment_model, ('embankment_mesh_nodes.csv', 'embankment_mesh_elements.csv', 'embankment_mesh_faces.csv',
                            'embankment_mesh_edges.csv', 'embankment_monitor.csv')),
    (fig_embankment_history, ('embankment_history.csv',)),
    (fig_embankment_fields, ('embankment_nodal_undrained.csv', 'embankment_nodal_t1e6.csv',
                             'embankment_nodal_t1e8.csv', 'embankment_gauss_t1e8.csv')),
)

if __name__ == '__main__':
    args = arguments('Figs. 12 to 14 of the article (embankment on a Cam-Clay foundation) from the CSV files of '
                     'EmbankmentConsolidation.')
    for function, files in FIGURES:
        missing = [f for f in files if not os.path.exists(os.path.join(args.rundir, f))]
        if missing:
            print(f'skipped {function.__name__}: missing {", ".join(missing)} '
                  f'(run EmbankmentConsolidation in {args.rundir})')
            continue
        function(args.rundir, args.outdir)
