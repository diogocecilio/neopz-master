"""Figures of the FLAC3D triaxial tests (Sect. 6.3 of the article) from the CSV files of FLAC3DTriaxial.

Usage (run the executable first; it writes its CSV files to the current directory):

    python3 <neopz>/Projects/FLAC3DTriaxial/plot_figures.py [run directory] [-o output directory]

The run directory (default: the current directory) holds the CSV files; the figures are written as PDF and PNG
to <run directory>/figures, or to the output directory.

Figures produced:

- fig_single_element_model: the single Hex20-Hex8 element (unit cube, 2 x 2 x 2 Gauss points) of the triaxial
  tests with its boundary conditions: (a) drained tests (RS2, Sect. 6.1, and FLAC3D, Sect. 6.3): pore pressure
  prescribed as zero at the eight vertices; (b) undrained tests (FLAC3D): no drained face, k = 0, Dt = 0. Symmetry
  planes x = 0 and y = 0 and base z = 0 (hidden, grey), cell pressure on the faces x = 1 and y = 1 (green),
  controlled vertical displacement of the top (blue); vertex nodes (u and p_w) and mid-edge nodes (u only). The
  meshes are read from flac3d_mesh_<drained|undrained>_*.csv (mcc::WriteMeshCSV), the faces are coloured by their
  boundary ids (FLAC3DTriaxial.h) and the cell pressure of the FLAC3D tests is read from flac3d_summary.csv.
- fig06_flac3d_triaxial: Fig. 6, drained (panels a, b) and undrained (panels c, d) triaxial tests with R = 1.6
  and R = 8: q against eps_a and stress paths in the p'-q plane, with the critical state line and the initial
  yield surface of R = 1.6. Solid lines: this work (flac3d_<drained|undrained>_R<1.6|8>.csv, columns eps_a, p_eff,
  q, v, u); dashed lines: closed-form solutions written by the executable (flac3d_<test>_closed.csv); squares:
  final states of FLAC3D (column flac3d of flac3d_table5.csv); M and p'0 from flac3d_summary.csv.

The script also prints the final states, the peaks and the closed-form values at the same eps_a (Table 5).
"""
import csv
import os
import sys

import numpy as np
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'Common'))
from mcc_figstyle import C1, C2, INK, INK2, MUTED, TEXTW, arguments, plt, read_csv, save  # noqa: E402
import mcc_hexmodel as hexmodel  # noqa: E402

# boundary ids of FLAC3DTriaxial.h
EX0, EX1, EY0, EY1, EZ0, EZ1, EPDRAINED = -21, -22, -23, -24, -25, -26, -30
TESTS = [(t, R) for t in ('drained', 'undrained') for R in (1.6, 8.0)]


def read_rows(path):
    """CSV file with text columns: list of dicts."""
    if not os.path.exists(path):
        sys.exit(f'file not found: {path} (run the executable in this directory first)')
    with open(path, newline='') as f:
        return list(csv.DictReader(f))


def summary(run):
    """flac3d_summary.csv: dict test name -> row (numbers as float)."""
    return {r['test']: {k: float(v) for k, v in r.items() if k != 'test'}
            for r in read_rows(os.path.join(run, 'flac3d_summary.csv'))}


def name(test, R):
    return f'{test}_R{R:g}'


# =========================================================================================== model figure
def cube_facecolor(ids):
    """Colour of a face of the cube from its boundary ids."""
    if EZ1 in ids:
        return hexmodel.FACE_TOP
    if EX1 in ids or EY1 in ids:
        return hexmodel.FACE_LATERAL
    return hexmodel.FACE_RESTRAINED


def draw_element(ax, mesh, view, labels, title, drained=True, legend=True):
    """One panel of the model figure: the cube with the arrows of the loading and the labels of the conditions.

    labels: texts of the top (controlled displacement), of the lateral faces (cell pressure), of the hidden faces and
    of the drainage. The vertex nodes are filled squares where p_w = 0 is prescribed (drained) and open squares where
    p_w is unknown (undrained). Returns the limits (xmin, xmax, ymin, ymax) of the drawing and labels."""
    P = view.point
    box = hexmodel.draw_model(ax, mesh, view, cube_facecolor, open_vertices=() if drained else mesh.vertices)
    # controlled vertical displacement of the top: arrows above the face z = 1
    for x, y in ((0.15, 0.15), (0.85, 0.15), (0.5, 0.5), (0.15, 0.85), (0.85, 0.85)):
        hexmodel.arrow(ax, P(x, y, 1.38), P(x, y, 1.03), color=INK)
    # cell pressure on the faces x = 1 and y = 1: arrows normal to the faces at their outer edges (silhouette)
    for z in (0.2, 0.5, 0.8):
        hexmodel.arrow(ax, P(1.5, 0.0, z), P(1.02, 0.0, z), color=C2)
        hexmodel.arrow(ax, P(0.0, 1.5, z), P(0.0, 1.02, z), color=C2)
    # axes, below the element on the left
    o = np.array([box[0] - 0.30, box[2] - 0.06])
    for k, lab in enumerate('xyz'):
        d = np.zeros(3)
        d[k] = 0.30
        v = view(d)[0] - view(np.zeros(3))[0]
        hexmodel.arrow(ax, o, o + v, color=INK2, lw=0.6, ms=5)
        ax.text(*(o + 1.35 * v), f'${lab}$', fontsize=7, color=INK2, ha='center', va='center')
    # labels in a column to the right of the drawing
    xt = P(0.0, 1.5, 0.5)[0] + 0.12
    ax.annotate(labels['top'], xy=tuple(P(0.3, 0.75, 1.0)), xytext=(xt, P(0.5, 0.5, 1.38)[1] - 0.04),
                textcoords='data', fontsize=6.6, color=INK, ha='left', va='center',
                arrowprops=dict(arrowstyle='-', color=C1, lw=0.7, shrinkA=2, shrinkB=0), zorder=8)
    ax.annotate(labels['lateral'], xy=tuple(P(0.3, 1.0, 0.62)), xytext=(xt, P(0, 1, 0.62)[1]),
                textcoords='data', fontsize=6.6, color=INK, ha='left', va='center',
                arrowprops=dict(arrowstyle='-', color=C2, lw=0.7, shrinkA=2, shrinkB=0), zorder=8)
    y = P(0, 1, 0.0)[1] - 0.06
    ax.text(xt, y, labels['hidden'], fontsize=6.6, color=INK2, ha='left', va='center')
    ax.text(xt, y - 0.34, labels['drainage'], fontsize=6.6, color=INK, ha='left', va='center')
    yl = y - 0.66
    if legend:
        hexmodel.node_legend(ax, xt + 0.03, yl, 0.13, open_vertex=not drained,
                             vertex='vertex node: $u$, $p_w$ = 0' if drained else 'vertex node: $u$ and $p_w$')
        ax.text(xt - 0.0, yl - 0.31, '2 × 2 × 2 Gauss points', fontsize=6.6, color=INK2, ha='left', va='center')
    ytop = P(0.5, 0.5, 1.38)[1] + 0.22
    ax.text(box[0] - 0.42, ytop, title, fontsize=8, va='bottom', ha='left')
    return [box[0] - 0.45, xt + 1.62, yl - 0.40, ytop + 0.18]


def fig_model(run, out):
    """Single element of the triaxial tests with the boundary conditions of the drained and undrained tests."""
    S = summary(run)
    p0 = S['drained_R1.6']['p0']
    view = hexmodel.View(24, 33)
    fig, axs = plt.subplots(1, 2, figsize=(TEXTW, 2.45))
    panels = (
        ('drained', True, '(a) drained tests: RS2 (Sect. 6.1) and FLAC3D (Sect. 6.3)',
         {'top': 'top $z$ = 1 m: $u_z$ prescribed\n($\\varepsilon_a = -u_z$)',
          'lateral': "faces $x$ = 1, $y$ = 1: total cell\npressure $\\sigma_c = p'_0$ (RS2: 200,\n"
                     "100 kPa; FLAC3D: %g kPa)" % p0,
          'hidden': 'hidden faces: $u_x$ = 0 on $x$ = 0,\n$u_y$ = 0 on $y$ = 0, $u_z$ = 0 on $z$ = 0',
          'drainage': 'drained: $p_w$ = 0 at the eight\nvertices (all faces)'}),
        ('undrained', False, '(b) undrained tests: FLAC3D (Sect. 6.3)',
         {'top': 'top $z$ = 1 m: $u_z$ prescribed\n($\\varepsilon_a = -u_z$)',
          'lateral': "faces $x$ = 1, $y$ = 1: total cell\npressure $\\sigma_c = p'_0$ = %g kPa" % p0,
          'hidden': 'hidden faces: $u_x$ = 0 on $x$ = 0,\n$u_y$ = 0 on $y$ = 0, $u_z$ = 0 on $z$ = 0',
          'drainage': 'undrained: no drained face,\n$k$ = 0, $\\Delta t$ = 0, $1/M_B = n/K_w$'}))
    for ax, (tag, drained, title, labels) in zip(axs, panels):
        mesh = hexmodel.Mesh(run, f'flac3d_mesh_{tag}')
        drained_faces = len(mesh.faces_with(EPDRAINED))
        if drained != (drained_faces == 6):
            sys.exit(f'flac3d_mesh_{tag}_faces.csv: {drained_faces} drained faces, expected {6 if drained else 0}')
        lim = draw_element(ax, mesh, view, labels, title, drained)
        ax.set_aspect('equal')
        ax.set_xlim(lim[0], lim[1])
        ax.set_ylim(lim[2], lim[3])
        ax.axis('off')
        print('model %-9s: %d element, %d vertices, %d edges, %d boundary faces, %d drained faces'
              % (tag, mesh.nelements, len(mesh.vertices), len(mesh.edges), len(mesh.faces), drained_faces))
    fig.tight_layout(w_pad=0.2)
    save(fig, out, 'fig_single_element_model')


# =========================================================================================== Fig. 6
def fig_flac3d(run, out):
    """Fig. 6: q-eps_a and p'-q of the drained and undrained tests with R = 1.6 and R = 8."""
    S = summary(run)
    M = S['drained_R1.6']['M']
    hist = {k: read_csv(os.path.join(run, f'flac3d_{name(*k)}.csv')) for k in TESTS}
    closed = {k: read_csv(os.path.join(run, f'flac3d_{name(*k)}_closed.csv')) for k in TESTS}
    t5 = {}
    for r in read_rows(os.path.join(run, 'flac3d_table5.csv')):
        t5[(r['test'], r['quantity'])] = {k: float(v) for k, v in r.items() if k not in ('test', 'quantity')}
    p0 = hist[('drained', 1.6)]['p_eff'][0]
    pc16 = S['drained_R1.6']['pc0']                                # p'c0 of R = 1.6 (initial surface)
    fig, axs = plt.subplots(1, 4, figsize=(TEXTW, 1.95))
    cols = {1.6: C1, 8.0: C2}
    for R in (8.0, 1.6):                       # R = 1.6 on top (the paths coincide at the beginning)
        z = 3 if R == 1.6 else 2
        for j, t in enumerate(('drained', 'undrained')):
            h, a = hist[(t, R)], closed[(t, R)]
            n = name(t, R)
            flac = (t5[(n, 'p_eff')]['flac3d'], t5[(n, 'q')]['flac3d'])
            axs[2 * j].plot(h['eps_a'], h['q'], color=cols[R], lw=1.5, zorder=z)
            axs[2 * j].plot(a['eps_a'], a['q'], color=INK, lw=0.7, ls=(0, (3, 2)), zorder=5)
            axs[2 * j + 1].plot(h['p_eff'], h['q'], color=cols[R], lw=1.5, zorder=z)
            axs[2 * j + 1].plot(a['p_eff'], a['q'], color=INK, lw=0.7, ls=(0, (3, 2)), zorder=5)
            axs[2 * j + 1].plot(*flac, 's', ms=4.8, mfc='none', mec=cols[R], mew=1.0, zorder=6)
            s = S[n]
            k = int(np.argmax(h['q']))
            peak = '' if np.isnan(s['q_peak_closed']) else ' (closed form %.4f at %.4f)' % (s['q_peak_closed'],
                                                                                            s['eps_a_peak_closed'])
            print("%-9s R = %-3g eps_a = %g: p' = %.5f q = %.5f kPa (closed form %.5f, %.5f; FLAC3D %g, %g), "
                  "peak q = %.4f kPa at eps_a = %g%s; %.4f evaluations per increment"
                  % (t, R, h['eps_a'][-1], h['p_eff'][-1], h['q'][-1], t5[(n, 'p_eff')]['closed_form'],
                     t5[(n, 'q')]['closed_form'], *flac, h['q'][k], h['eps_a'][k], peak, s['mean_evaluations']))
    for ax in (axs[1], axs[3]):
        pm = np.linspace(0, 16, 2)
        ax.plot(pm, M * pm, color=MUTED, lw=0.8, zorder=1)
        ax.text(12.6, M * 12.6 - 2.6, 'CSL', fontsize=7, color=INK2, ha='left')
        pe = np.linspace(0, pc16, 200)                                # initial surface with R = 1.6 (p'c0 = 8)
        ax.plot(pe, M * np.sqrt(np.clip(pe * (pc16 - pe), 0, None)), color=MUTED, lw=0.7, ls=(0, (3, 2)),
                zorder=1)
    axs[1].annotate('initial surface\n($R$ = 1.6)', xy=(2.2, 3.6), xytext=(0.3, 13.6), fontsize=6.3, color=INK2,
                    ha='left', va='top', arrowprops=dict(arrowstyle='-', color=MUTED, lw=0.6))
    axs[1].text(9.2, 4.0, "drained path\n$q = 3(p' - p'_0)$", fontsize=6.3, color=INK2)
    axs[3].annotate("undrained: $p'$ constant\nuntil yield", xy=(p0, 2.2), xytext=(8.3, 4.6), fontsize=6.3,
                    color=INK2, ha='left', va='top', arrowprops=dict(arrowstyle='-', color=MUTED, lw=0.6))
    for ax, key in ((axs[0], 'drained'), (axs[2], 'undrained')):
        h8 = hist[(key, 8.0)]
        k8 = int(np.argmax(h8['q']))
        ax.annotate('peak %.1f kPa' % h8['q'][k8], (h8['eps_a'][k8], h8['q'][k8]),
                    xytext=(0.30, 0.93), textcoords='axes fraction', fontsize=6.5, color=INK2,
                    arrowprops=dict(arrowstyle='-', color=INK2, lw=0.6))
        h16 = hist[(key, 1.6)]
        ax.text(h16['eps_a'][-1] * 0.97, h16['q'][-1] - 2.4, '$R$ = 1.6', fontsize=6.8, color=C1, ha='right')
        ax.text(h8['eps_a'][-1] * 0.97, h8['q'][-1] + 1.0, '$R$ = 8', fontsize=6.8, color=C2, ha='right')
    axs[0].set(xlabel='$\\varepsilon_a$', ylabel='$q$ (kPa)', xlim=(0, 0.5), ylim=(0, 19))
    axs[1].set(xlabel="$p'$ (kPa)", ylabel='$q$ (kPa)', xlim=(0, 16), ylim=(0, 19))
    axs[2].set(xlabel='$\\varepsilon_a$', ylabel='$q$ (kPa)', xlim=(0, 0.1), ylim=(0, 19))
    axs[3].set(xlabel="$p'$ (kPa)", ylabel='$q$ (kPa)', xlim=(0, 16), ylim=(0, 19))
    for k, (ax, t) in enumerate(zip(axs, ('drained', 'drained', 'undrained', 'undrained'))):
        ax.set_title(f'({"abcd"[k]}) {t}', fontsize=8)
    handles = [Line2D([], [], color=C1, lw=1.3, label='$R = 1.6$'), Line2D([], [], color=C2, lw=1.3, label='$R = 8$'),
               Line2D([], [], color=INK, lw=0.8, ls=(0, (4, 2)), label='closed form'),
               Line2D([], [], ls='', marker='s', ms=4.5, mfc='white', mec=INK2, label='FLAC3D (final state)')]
    fig.legend(handles=handles, loc='lower center', ncol=4, bbox_to_anchor=(0.5, -0.07), fontsize=7.5)
    fig.tight_layout(rect=(0, 0.05, 1, 1), w_pad=0.6)
    save(fig, out, 'fig06_flac3d_triaxial')


if __name__ == '__main__':
    args = arguments('Figures of the FLAC3D triaxial tests (model and Fig. 6) from the CSV files of FLAC3DTriaxial.')
    fig_model(args.rundir, args.outdir)
    fig_flac3d(args.rundir, args.outdir)
