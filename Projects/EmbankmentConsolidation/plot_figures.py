"""Figures 10, 11 and 12 of the article (embankment on a Cam-Clay foundation) from the files written by
EmbankmentConsolidation.

Usage (run the executable first; it writes its CSV files to the current directory):

    python3 <neopz>/Projects/EmbankmentConsolidation/plot_figures.py [run directory] [-o output directory]

The run directory (default: the current directory) holds the CSV files; the figures are written as PDF and PNG
to <run directory>/figures, or to the output directory. A figure whose CSV files are missing is skipped with a
message. This script draws the three states of the article; the fields at every converged state (36 states) are
in the VTK series written by the executable in vtk/embankment and vtk/embankment_elastic, to be opened in
ParaView (README.md, section "Viewing the solution in ParaView").

Figures produced (same panels, curves, markers and annotations as fig_aterro_model, fig_aterro_hist and
fig_aterro_fields of figs.py of the Python code):

- fig10_embankment_model: Fig. 10, the model: the 20 x 10 element mesh (vertices of the Q8-Q4 elements, read
  from embankment_nodal_undrained.csv), the strip load, the drained top, the boundary conditions, the
  settlement points at x = 0, 2, 4, 6 m and the zones pp1 and pp2.
- fig11_embankment_history: Fig. 11, histories from the end of the undrained loading (t = 0) to t = 1e8 s
  (embankment_history.csv): (a) settlements of the top at x = 0, 2, 4, 6 m, (b) pore pressures in the zones pp1
  and pp2, with pp2 against log t in the inset (Mandel-Cryer effect). Markers: the FLAC3D histories, digitized
  (reference/flac_historicos_digitalizados.json, the file dados/ of the Python code), with the FLAC3D values at
  the end of the loading (Table 7) at t = 0.
- fig12_embankment_fields: Fig. 12, fields on the vertices of the mesh (embankment_nodal_<state>.csv,
  contours on the two triangles of each quadrilateral, as in the Python code): (a) excess pore pressure at
  the end of the undrained loading, (b) at t = 1e6 s, (c) plastic integration points at t = 1e8 s
  (embankment_gauss_t1e8.csv: type 1 subcritical, 2 supercritical), (d) settlement at t = 1e8 s.
"""
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

# FLAC3D values at the end of the undrained loading (t = 0; Table 7 of the article), the first point of the
# FLAC3D histories (the digitized curves start at t = 2.5e5 s)
FLAC_T0 = dict(uz_x0=0.140, uz_x2=0.135, uz_x4=0.055, uz_x6=-0.0416, pp1=18.1, pp2=62.4)
X_SETTLEMENT = (0, 2, 4, 6)  # settlement points of the top (m), columns s_x0 ... s_x6 of embankment_history.csv


# =========================================================================================== mesh and data
def vertex_mesh(nodal):
    """Quadrilaterals of the mesh from the vertices of embankment_nodal_<state>.csv (structured grid of
    EmbankmentConsolidation::CreateGeoMesh, elements row by row).

    Returns the quadrilaterals as four vertex indices (rows of the CSV file), counter-clockwise from the lower
    left corner, as the first four nodes of the elements of the Python code."""
    x, y = np.round(nodal['x'], 9), np.round(nodal['y'], 9)
    xs, ys = np.unique(x), np.unique(y)
    index = {(a, b): k for k, (a, b) in enumerate(zip(x, y))}
    if len(index) != len(x) or len(x) != len(xs) * len(ys):
        sys.exit('the vertices of the nodal CSV file do not form a structured grid')
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


# =========================================================================================== Fig. 10
def fig_embankment_model(run, out):
    """Fig. 10: the model, the boundary conditions and the monitoring points (fig_aterro_model of figs.py)."""
    nodal = read_csv(os.path.join(run, 'embankment_nodal_undrained.csv'))
    C = np.column_stack([nodal['x'], nodal['y']])
    quads = vertex_mesh(nodal)
    fig, ax = plt.subplots(figsize=(TEXTW, 2.75))
    for q in quads:
        ax.add_patch(MPoly(C[q], closed=True, fc='#f4f3ef', ec=INK2, lw=0.45))
    # load
    for x in np.linspace(0.2, 3.8, 10):
        ax.add_patch(FancyArrowPatch((x, 11.3), (x, 10.15), arrowstyle='-|>', mutation_scale=5, color=C2, lw=0.7))
    ax.plot([0, 4], [11.3, 11.3], color=C2, lw=0.9)
    ax.text(5.0, 11.3, '$q$ = 50 kPa', ha='left', va='center', fontsize=7.5)
    ax.plot([4.05, 20], [10.07, 10.07], color=C1, lw=1.6)
    ax.plot([0, 4.0], [10.07, 10.07], color=C1, lw=1.6)
    ax.text(12.6, 10.35, 'drained surface, water table ($p_w = 0$)', fontsize=7.5, color=INK, ha='center')
    # water table symbol
    ax.add_patch(MPoly([[18.95, 10.95], [19.55, 10.95], [19.25, 10.2]], closed=True, fc='white', ec=C1, lw=0.8))
    # material (white box over the mesh)
    ax.text(13.0, 4.6, "clay layer (MCC), saturated\n"
                       "$p'_{c0}$ = 160 kPa (uniform)\n"
                       "$\\sigma'_v = 13d$, $\\sigma'_h = 6.1d$ kPa\n"
                       "($d$: depth in m)",
            ha='center', va='center', fontsize=7.2, color=INK, linespacing=1.35,
            bbox=dict(boxstyle='round,pad=0.35', fc='white', ec=GRID, lw=0.6, alpha=0.95), zorder=7)
    # supports
    for y in np.linspace(0.5, 9.5, 8):
        ax.plot(-0.35, y, 'o', ms=3.5, mfc='white', mec=INK, mew=0.7)
        ax.plot(20.35, y, 'o', ms=3.5, mfc='white', mec=INK, mew=0.7)
    ax.plot([-0.6, -0.6], [0, 10], color=INK, lw=0.9)
    ax.plot([20.6, 20.6], [0, 10], color=INK, lw=0.9)
    ax.plot([0, 20], [-0.25, -0.25], color=INK, lw=1.6)
    for x in np.linspace(0, 20, 41):
        ax.plot([x, x - 0.3], [-0.25, -0.6], color=INK, lw=0.6)
    ax.text(10, -0.75, 'fixed, impermeable base', ha='center', va='top', fontsize=7.5)
    ax.text(-1.0, 5, 'symmetry: $u_x = 0$', rotation=90, ha='right', va='center', fontsize=7.5)
    ax.text(21.0, 5, '$u_x = 0$', rotation=90, ha='left', va='center', fontsize=7.5)
    # monitoring
    for x in X_SETTLEMENT:
        ax.plot(x, 10, 'v', ms=5, color=C3, mec=INK, mew=0.4, zorder=5)
    ax.text(6.3, 9.55, 'settlement points', fontsize=7, color=INK, va='center')
    for (x, y), lab in (((0.5, 9.5), 'pp1'), ((1.5, 7.5), 'pp2')):
        x0, y0 = np.floor(x), np.floor(y)
        ax.add_patch(MPoly([[x0, y0], [x0 + 1, y0], [x0 + 1, y0 + 1], [x0, y0 + 1]], closed=True, fc=C4, ec=INK,
                           lw=0.5, alpha=0.8, zorder=4))
        ax.text(x0 + 1.15, y0 + 0.5, lab, fontsize=7.5, va='center', zorder=6)
    # dimensions
    ax.annotate('', xy=(0, -2.0), xytext=(20, -2.0), arrowprops=dict(arrowstyle='<->', color=INK2, lw=0.6,
                                                                   mutation_scale=6))
    ax.text(10, -2.15, '20 m', ha='center', va='top', fontsize=7.5)
    ax.annotate('', xy=(22.2, 0), xytext=(22.2, 10), arrowprops=dict(arrowstyle='<->', color=INK2, lw=0.6,
                                                                    mutation_scale=6))
    ax.text(22.4, 5, '10 m', ha='left', va='center', fontsize=7.5, rotation=90)
    ax.annotate('', xy=(0, 12.4), xytext=(4, 12.4), arrowprops=dict(arrowstyle='<->', color=INK2, lw=0.6,
                                                                   mutation_scale=6))
    ax.text(4.2, 12.4, '4 m', ha='left', va='center', fontsize=7.5)
    ax.set_xlim(-2.2, 23.2)
    ax.set_ylim(-2.9, 12.9)
    ax.set_aspect('equal')
    ax.axis('off')
    print('Fig. 10: %d elements, %d vertices' % (len(quads), len(C)))
    save(fig, out, 'fig10_embankment_model')


# =========================================================================================== Fig. 11
def fig_embankment_history(run, out):
    """Fig. 11: settlements and pore pressures against t, with the FLAC3D histories (fig_aterro_hist of figs.py)."""
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
    print('Fig. 11: t = 1e8 s: settlements %s m, pp1 %.4f, pp2 %.4f kPa; pp2 %.3f -> %.3f kPa at t = %.3g s'
          % (np.array2string(hc[-1, 1:5], precision=6), hc[-1, 5], hc[-1, 6], hc[0, 6], hc[km, 6], hc[km, 0]))
    save(fig, out, 'fig11_embankment_history')


# =========================================================================================== Fig. 12
def fig_embankment_fields(run, out):
    """Fig. 12: excess pore pressure at the end of the loading and at t = 1e6 s, plastic integration points and
    settlement at t = 1e8 s (fig_aterro_fields of figs.py)."""
    nodal = {tag: read_csv(os.path.join(run, f'embankment_nodal_{tag}.csv')) for tag in ('undrained', 't1e6', 't1e8')}
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
        print('Fig. 12: %s: largest excess pore pressure %.2f kPa at (%g, %g)' % (tag, pex[k], n['x'][k], n['y'][k]))
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
    # plastic integration points at t = 1e8 s
    C = np.column_stack([n['x'], n['y']])
    for q in quads:
        axd.add_patch(MPoly(C[q], closed=True, fc='white', ec=GRID, lw=0.4))
    typ = np.rint(gauss['type']).astype(int)
    sel1, sel2 = typ == 1, typ == 2
    axd.plot(gauss['x'][sel1], gauss['y'][sel1], 'o', ms=2.6, mfc=C1, mec='white', mew=0.3,
             label=f'subcritical ({sel1.sum()})')
    axd.plot(gauss['x'][sel2], gauss['y'][sel2], '^', ms=3.0, mfc=C2, mec='white', mew=0.3,
             label=f'supercritical ({sel2.sum()})')
    axd.set_title('(c) plastic Gauss points, $t$ = 10$^8$ s', fontsize=8)
    axd.legend(loc='upper right', fontsize=6.8, frameon=True, framealpha=0.9, edgecolor='none', handletextpad=0.3)
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
    print('Fig. 12: t = 1e8 s: plastic points %d subcritical, %d supercritical; settlement from %.4f to %.4f m, '
          'heave of the top up to %.4f mm' % (sel1.sum(), sel2.sum(), s.min(), s.max(), 1e3 * heave))
    save(fig, out, 'fig12_embankment_fields')


# (function, CSV files)
FIGURES = (
    (fig_embankment_model, ('embankment_nodal_undrained.csv',)),
    (fig_embankment_history, ('embankment_history.csv',)),
    (fig_embankment_fields, ('embankment_nodal_undrained.csv', 'embankment_nodal_t1e6.csv',
                             'embankment_nodal_t1e8.csv', 'embankment_gauss_t1e8.csv')),
)

if __name__ == '__main__':
    args = arguments('Figs. 10 to 12 of the article (embankment on a Cam-Clay foundation) from the CSV files of '
                     'EmbankmentConsolidation.')
    for function, files in FIGURES:
        missing = [f for f in files if not os.path.exists(os.path.join(args.rundir, f))]
        if missing:
            print(f'skipped {function.__name__}: missing {", ".join(missing)} '
                  f'(run EmbankmentConsolidation in {args.rundir})')
            continue
        function(args.rundir, args.outdir)
