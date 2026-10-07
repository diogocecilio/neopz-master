"""Figure 7 of the article (Terzaghi consolidation) from the CSV files written by TerzaghiConsolidation.

Usage (run the executable first; it writes its CSV files to the current directory):

    python3 <neopz>/Projects/TerzaghiConsolidation/plot_figures.py [run directory] [-o output directory]

The run directory (default: the current directory) holds the CSV files; the figures are written as PDF and PNG
to <run directory>/figures, or to the output directory.

Figure produced:

- fig07_terzaghi_consolidation: Fig. 7, (a) the model: column of 1 x 1 x 10 Hex20-Hex8 elements (3 x 3 x 3 Gauss
  points) with the loaded and drained top, the lateral faces with zero normal displacement and impermeable and the
  fixed impermeable base, the vertices of the edge x = y = 0 where the pore pressure is monitored and the top
  vertex of the settlement (terzaghi_mesh_*.csv, mcc::WriteMeshCSV, drawn with Common/mcc_hexmodel.py);
  (b) excess pore pressure p_w/q along the edge x = y = 0 for T = 0.001, 0.01, 0.1 and 0.5 at the vertices
  (terzaghi_fig7a.csv) against the exact series solution (terzaghi_exact_isochrones.csv); (c) degree of
  consolidation U = w/w_inf of the top settlement against T = c_v t/H^2 (terzaghi_fig7b.csv) against the exact
  series solution (terzaghi_exact_degree.csv). c_v, H and w_inf of the annotation are read from
  terzaghi_summary.csv.

The same panels (b) and (c), without the model, were Fig. 7 of v0.6 of the article (plane strain Q8-Q4 column);
the nodal values of the 3D column are identical (see README.md).
"""
import os
import sys

import numpy as np
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'Common'))
from mcc_figstyle import C1, C2, C3, C4, INK, INK2, TEXTW, arguments, plt, read_csv, save  # noqa: E402
import mcc_hexmodel as hexmodel  # noqa: E402  (Common/mcc_hexmodel.py: drawing of the hexahedral models)

# boundary ids of TerzaghiConsolidation.h
EX0, EX1, EY0, EY1, EZ0, EZ1, EPZ1 = -21, -22, -23, -24, -25, -26, -36


def column_facecolor(ids):
    """Colour of a boundary face of the column from its boundary ids."""
    return hexmodel.FACE_TOP if EZ1 in ids else hexmodel.FACE_RESTRAINED


def draw_column(ax, run, S):
    """Panel (a): the column with its boundary conditions and the monitored nodes."""
    mesh = hexmodel.Mesh(run, 'terzaghi_mesh')
    H, q = S['H'], S['q']
    if len(mesh.faces_with(EPZ1)) != 1 or len(mesh.faces_with(EZ1)) != 1:
        sys.exit('terzaghi_mesh_faces.csv: the top must be one loaded and drained face')
    view = hexmodel.View(16, 215)                     # the monitored edge x = y = 0 in front
    P = view.point
    box = hexmodel.draw_model(ax, mesh, view, column_facecolor, lw=0.5, nodes=False, hidden=False)
    # load on the top
    for x, y in ((0.0, 0.0), (1.0, 0.0), (0.0, 1.0)):
        hexmodel.arrow(ax, P(x, y, H + 1.7), P(x, y, H + 0.08), color=INK, ms=4.5, lw=0.6)
    # monitored nodes: vertices of the edge x = y = 0 (p_w, panel b) and its top vertex (settlement, panel c)
    edge = [k for k in mesh.vertices if abs(mesh.X[k, 0]) < 1e-9 and abs(mesh.X[k, 1]) < 1e-9]
    Pe = view(mesh.X[edge])
    ax.plot(Pe[:, 0], Pe[:, 1], 'o', ms=3.3, mfc='none', mec=C3, mew=0.8, zorder=8)
    ktop = max(edge, key=lambda k: mesh.X[k, 2])
    ax.plot(*view(mesh.X[ktop])[0], '*', ms=7.5, color=C3, mec=INK, mew=0.4, zorder=9)
    # base: fixed and impermeable (ground hatch)
    b0, b1 = P(1, 0, 0), P(0, 1, 0)
    xb = np.linspace(min(b0[0], b1[0]) - 0.15, max(b0[0], b1[0]) + 0.15, 2)
    yb = box[2] - 0.15
    ax.plot(xb, [yb, yb], color=INK2, lw=0.8)
    for x in np.linspace(xb[0], xb[1], 7)[:-1]:
        ax.plot([x, x + 0.18], [yb, yb - 0.3], color=INK2, lw=0.5)
    # dimension H
    xd = box[0] - 0.35
    ax.annotate('', xy=(xd, P(0, 0, 0)[1]), xytext=(xd, P(0, 0, H)[1]),
                arrowprops=dict(arrowstyle='<->', color=INK2, lw=0.6, mutation_scale=6, shrinkA=0, shrinkB=0))
    ax.text(xd - 0.12, 0.5 * (P(0, 0, 0)[1] + P(0, 0, H)[1]), '$H$ = %g m' % H, rotation=90, fontsize=6.6,
            ha='right', va='center', color=INK2)
    # labels to the right of the column
    xt = box[1] + 0.35
    fs = 6.4
    lab = [(P(0.5, 0.5, H), P(0, 0, H + 1.0)[1], 'top: load $q$ = %g kPa,\ndrained ($p_w$ = 0)' % q, C1),
           (P(0.0, 0.5, 0.62 * H), P(0, 0, 0.66 * H)[1], 'lateral faces: $u_n$ = 0,\nimpermeable', INK2),
           (tuple(view(mesh.X[edge])[3]), P(0, 0, 0.36 * H)[1], 'vertices of $x = y = 0$:\n$p_w$ in (b)', C3),
           (tuple(view(mesh.X[ktop])[0]), P(0, 0, 0.86 * H)[1], 'settlement $w$ (c)', C3)]
    for xy, y, s, c in lab:
        ax.annotate(s, xy=tuple(xy), xytext=(xt, y), textcoords='data', fontsize=fs, color=INK, ha='left',
                    va='center', arrowprops=dict(arrowstyle='-', color=c, lw=0.6, shrinkA=2, shrinkB=2), zorder=10)
    ax.text(xt, P(0, 0, 0.12 * H)[1], 'base: $u_z$ = 0,\nimpermeable', fontsize=fs, color=INK2, ha='left',
            va='center')
    ax.text(box[0] - 0.9, yb - 0.75, '1 × 1 × %d Hex20–Hex8 elements,\n3 × 3 × 3 Gauss points' % mesh.nelements,
            fontsize=fs, color=INK2, ha='left', va='top')
    ax.set_aspect('equal')
    ax.set_xlim(box[0] - 0.95, xt + 3.6)
    ax.set_ylim(yb - 2.3, P(0, 0, H + 1.8)[1] + 0.2)
    ax.axis('off')
    ax.set_title('(a) model', fontsize=8, loc='left')
    print('model: %d elements, %d vertices, %d edges, %d boundary faces' % (mesh.nelements, len(mesh.vertices),
                                                                           len(mesh.edges), len(mesh.faces)))


def fig_terzaghi(run, out):
    """Fig. 7: model, pore pressure isochrones at x = y = 0 and degree of consolidation."""
    a7 = read_csv(os.path.join(run, 'terzaghi_fig7a.csv'))
    b7 = read_csv(os.path.join(run, 'terzaghi_fig7b.csv'))
    iso = read_csv(os.path.join(run, 'terzaghi_exact_isochrones.csv'))
    deg = read_csv(os.path.join(run, 'terzaghi_exact_degree.csv'))
    S = {k: v[0] for k, v in read_csv(os.path.join(run, 'terzaghi_summary.csv')).items()}
    H = S['H']
    fig = plt.figure(figsize=(TEXTW, 2.75))
    gs = fig.add_gridspec(1, 3, width_ratios=(0.92, 1.0, 1.25), left=0.0, right=0.995, bottom=0.15, top=0.9,
                          wspace=0.3)
    draw_column(fig.add_subplot(gs[0]), run, S)
    ax = fig.add_subplot(gs[1])
    cols = [C1, C2, C3, C4]
    rows = []
    for k, T in enumerate(np.unique(a7['T'])):
        m = a7['T'] == T
        o = np.argsort(a7['z'][m])
        zs, p = a7['z'][m][o], a7['pw_over_q'][m][o]
        rows.append((T, zs, p))
        e = iso['T'] == T
        ax.plot(iso['pw_over_q'][e], iso['z'][e], color=cols[k], lw=1.1)
        ax.plot(p, zs, 'o', ms=3.6, mfc='white', mec=cols[k], mew=0.9)
        print('T = %-5g max |p_w/q - exact| at the vertices = %.5f' % (T, np.abs(p - a7['exact'][m][o]).max()))
    ax.set(xlabel='$p_w/q$', ylabel='$z$ (m)', xlim=(0, 1.62), ylim=(0, H))
    _, zs1, p1 = rows[0]
    ax.annotate('oscillation next to\nthe drained face', (p1[-2], zs1[-2]), xytext=(1.14, 8.5), fontsize=6.3,
                color=INK2, arrowprops=dict(arrowstyle='-', color=INK2, lw=0.6))
    ax.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_title('(b) excess pore pressure at $x = y = 0$', fontsize=8)
    handles = [Line2D([], [], color=cols[k], lw=1.1, marker='o', ms=3.6, mfc='white', mec=cols[k],
                      label=f'$T$ = {r[0]:g}') for k, r in enumerate(rows)]
    handles += [Line2D([], [], color=INK2, lw=1.1, label='exact'),
                Line2D([], [], ls='', marker='o', ms=3.6, mfc='white', mec=INK2, label='vertices')]
    ax.legend(handles=handles, loc='center left', bbox_to_anchor=(0.665, 0.4), fontsize=6.6, handlelength=1.5,
              labelspacing=0.3, borderaxespad=0.0)
    ax = fig.add_subplot(gs[2])
    Tt, Un = b7['T'][1:], b7['U_numerical'][1:]       # the first row is the initial state (before the load)
    ax.semilogx(deg['T'], deg['U'], color=INK, lw=0.9, ls=(0, (4, 2)), label='exact')
    ax.semilogx(Tt, Un, color=C1, lw=1.4, label='Hex20–Hex8')
    ax.set(xlabel='$T = c_v t / H^2$', ylabel='$U = w/w_\\infty$', xlim=(1e-5, 1), ylim=(0, 1))
    ax.text(0.03, 0.55, '$c_v = kE_{oed}$ = %.3g m$^2$/s, $H$ = %g m\n$w_\\infty = qH/E_{oed}$ = %.2f mm'
            % (S['cv'], H, 1e3 * S['w_inf']), fontsize=6.5, color=INK2, transform=ax.transAxes)
    Ti = 1e-4
    ax.annotate('offset at small $T$: diffusion layer\nthinner than one element', (Ti, np.interp(Ti, Tt, Un)),
                xytext=(2e-5, 0.22), fontsize=6.3, color=INK2, arrowprops=dict(arrowstyle='-', color=INK2, lw=0.6))
    ax.set_title('(c) degree of consolidation (top settlement)', fontsize=8)
    ax.legend(loc='upper left', fontsize=7)
    print('c_v = %.4g m2/s, w_inf = %.4f mm, U(T = 1) = %.4f (exact %.4f); %d increments, %.0f evaluations per '
          'increment' % (S['cv'], 1e3 * S['w_inf'], Un[-1], b7['U_exact'][-1], S['increments'], S['mean_evaluations']))
    save(fig, out, 'fig07_terzaghi_consolidation')


if __name__ == '__main__':
    args = arguments('Fig. 7 of the article (Terzaghi consolidation) from the CSV files of TerzaghiConsolidation.')
    fig_terzaghi(args.rundir, args.outdir)
