"""Figures 1 and 2 of the article from the CSV files written by the YieldSurfaceProjection executable.

Usage (run the executable first; it writes its CSV files to the current directory):

    python3 <neopz>/Projects/YieldSurfaceProjection/plot_figures.py [run directory] [-o output directory]

The run directory (default: the current directory) holds the CSV files; the figures are written as PDF and PNG
to <run directory>/figures, or to the output directory.

Figures produced (same panels, curves and annotations as fig_surface.py of the Python code of the article):

- fig01_mcc_surface: Fig. 1, MCC yield surface (M = 1, p'c = 100 kPa, pt = 0, omega = 1) with the critical
  state circle, (a) in principal stresses (compression axes) and (b) in rotated Haigh-Westergaard space.
  Files: fig1_surface.csv (121 x 97 grid, beta outer loop and xi inner loop) and fig1_critical_state_circle.csv.
- fig02_meridian_projection: Fig. 2, closest-point projection in the meridian plane, (a) subcritical trial state
  (compaction and hardening) and (b) supercritical trial state (dilation and softening): ellipses at the start
  and at the end of the step, critical state line, energy-norm contour through the projected state, trial and
  projected states and flow direction. Files: fig2a_subcritical_*.csv, fig2b_supercritical_*.csv and
  fig2_critical_state_line.csv.
"""
import os
import sys

import matplotlib
import numpy as np
from matplotlib.colors import LightSource
from mpl_toolkits.mplot3d import proj3d

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'Common'))
from mcc_figstyle import (C1, C2, GREEN_EDGE, GREEN_FACE, INK, INK2, MUTED, TEXTW, arguments,  # noqa: E402
                          plt, read_csv, save)

SQ3 = np.sqrt(3.0)


def grid(S, name):
    """Column of the Fig. 1 grid as a (n_beta, n_xi) array, the numpy meshgrid(xis, th) layout of fig_surface.py."""
    beta = S['beta']
    nxi = int(np.argmax(beta != beta[0])) or beta.size          # length of the inner (xi) loop
    return S[name].reshape(-1, nxi)


# ------------------------------------------------------------------------------------------- Fig. 1
def fig_surface(run, out):
    """Fig. 1: MCC yield surface in principal stresses (compression axes) and in RHW space (fig_surface.py)."""
    S = read_csv(os.path.join(run, 'fig1_surface.csv'))
    C = read_csv(os.path.join(run, 'fig1_critical_state_circle.csv'))
    pc = -S['xi'].min() / SQ3                                     # xi runs from sqrt3 pt = 0 to -sqrt3 p'c
    nb, nxi = grid(S, 'xi').shape                                 # 97 Lode angles x 121 values of xi
    xi_cs = C['sstar1'][0]                                        # section of the critical state circle
    rcs = np.hypot(C['sstar2'][0], C['sstar3'][0])                # rho_max

    fig = plt.figure(figsize=(TEXTW, 3.25))
    face = np.array(matplotlib.colors.to_rgb(GREEN_FACE))
    ls = LightSource(azdeg=300, altdeg=35)
    # (a) principal stresses, compression axes
    ax = fig.add_axes([0.0, 0.0, 0.5, 1.0], projection='3d')
    X, Y, Z = -grid(S, 'sigma1'), -grid(S, 'sigma2'), -grid(S, 'sigma3')
    rgb = ls.shade_rgb(np.ones(X.shape + (3,)) * face, Z + 0.4 * X - 0.4 * Y, blend_mode='soft', fraction=0.6)
    ax.plot_surface(X, Y, Z, facecolors=rgb, rstride=4, cstride=6, linewidth=0.25, edgecolor=GREEN_EDGE,
                    antialiased=True, shade=False, alpha=0.82)
    ax.plot(-C['sigma1'], -C['sigma2'], -C['sigma3'], color=C1, lw=1.6, zorder=10)
    L = 128
    for v, lab, f in ((np.array([1, 0, 0]), r'$-\sigma_1$', 1.1), (np.array([0, 1, 0]), r'$-\sigma_2$', 1.3),
                      (np.array([0, 0, 1]), r'$-\sigma_3$', 1.1)):
        ax.plot([0, f * L * v[0] / 1.1], [0, f * L * v[1] / 1.1], [0, f * L * v[2] / 1.1], color=INK2, lw=0.7)
        ax.text(*(f * L * v), lab, color=INK, fontsize=9, ha='center', va='center')
    hyd = np.array([0, 1.18 * pc])
    ax.plot(hyd, hyd, hyd, color=INK, lw=0.9)
    ax.text(1.24 * pc, 1.24 * pc, 1.24 * pc, r'$-\xi$', fontsize=9, color=INK)
    ax.scatter([0, pc], [0, pc], [0, pc], s=9, color=INK, depthshade=False, zorder=11)
    ax.text(-6, -6, -14, "$p'=0$", fontsize=7.5, color=INK, ha='right')
    ax.text(pc + 6, pc + 2, pc - 16, "$p'=p'_c$", fontsize=7.5, color=INK)
    lim = (-12, 132)
    ax.set_xlim(lim); ax.set_ylim(lim); ax.set_zlim(lim)
    ax.set_box_aspect((1, 1, 1))
    ax.view_init(elev=22, azim=-25)                               # about 60 degrees from the hydrostatic axis
    ax.set_axis_off()
    # label of the critical state section: point of the circle projected to 2D, label in a free corner
    pts = np.column_stack((-C['sigma1'], -C['sigma2'], -C['sigma3']))
    P2 = np.array([proj3d.proj_transform(*p, ax.get_proj())[:2] for p in pts])
    k = int(np.argmax(P2[:, 0] - P2[:, 1]))                      # rightmost and lowest point of the drawing
    ax.annotate("critical state section\n" + r"$\bar p = 0$, $p' = p'_c/2$", xy=tuple(P2[k]), xycoords='data',
                xytext=(0.66, 0.16), textcoords='axes fraction', fontsize=7.2, color=INK, ha='left', va='top',
                arrowprops=dict(arrowstyle='-', color=C1, lw=0.7))
    ax.text2D(0.02, 0.95, '(a) principal stresses (compression axes)', transform=ax.transAxes, fontsize=8)
    # (b) RHW space: sigma* = {xi, rho cos(beta), rho sin(beta)}
    ax = fig.add_axes([0.5, 0.0, 0.5, 1.0], projection='3d')
    Xr, Yr, Zr = grid(S, 'sstar1'), grid(S, 'sstar2'), grid(S, 'sstar3')
    rgb = ls.shade_rgb(np.ones(Xr.shape + (3,)) * face, Zr - 0.3 * Yr, blend_mode='soft', fraction=0.6)
    ax.plot_surface(Xr, Yr, Zr, facecolors=rgb, rstride=4, cstride=6, linewidth=0.25, edgecolor=GREEN_EDGE,
                    antialiased=True, shade=False, alpha=0.82)
    ax.plot(C['sstar1'], C['sstar2'], C['sstar3'], color=C1, lw=1.6, zorder=10)
    ax.plot([-1.12 * SQ3 * pc, 0.16 * SQ3 * pc], [0, 0], [0, 0], color=INK, lw=0.9)
    ax.text(0.2 * SQ3 * pc, 0, 0, r'$\sigma^*_1=\xi$', fontsize=9, color=INK, va='center')
    ax.plot([0, 0], [0, 70], [0, 0], color=INK2, lw=0.7)
    ax.text(0, 78, 0, r'$\sigma^*_2$', fontsize=9, color=INK, ha='center')
    ax.plot([0, 0], [0, 0], [0, 70], color=INK2, lw=0.7)
    ax.text(0, 0, 78, r'$\sigma^*_3$', fontsize=9, color=INK, ha='center')
    ax.scatter([0, -SQ3 * pc], [0, 0], [0, 0], s=9, color=INK, depthshade=False, zorder=11)
    ax.text(0, 0, -16, r'$\xi = 0$', fontsize=7.5, color=INK, ha='center', va='top')
    ax.text(-1.08 * SQ3 * pc, 0, -30, r"$\xi = -\sqrt{3}\,p'_c$", fontsize=7.5, color=INK, ha='center', va='top')
    ax.plot([xi_cs, xi_cs], [0, 0], [rcs, rcs + 26], color=C1, lw=0.7)
    ax.text(xi_cs, 0, rcs + 30, r"critical state: $\rho_{max} = %.1f$ kPa" % rcs, fontsize=7.2,
            color=INK, ha='center', va='bottom')
    ax.set_xlim(-1.15 * SQ3 * pc, 0.3 * SQ3 * pc); ax.set_ylim(-75, 75); ax.set_zlim(-75, 75)
    ax.set_box_aspect((1.45 * SQ3 * pc, 150, 150))
    ax.view_init(elev=22, azim=-62)
    ax.set_axis_off()
    ax.text2D(0.02, 0.95, '(b) rotated Haigh–Westergaard space', transform=ax.transAxes, fontsize=8)
    print('Fig. 1: grid %d x %d, p\'c = %.6g kPa, rho_max = %.6g kPa, max |Phi|/a^2 = %.2g'
          % (nxi, nb, pc, rcs, np.abs(S['phi_over_a2']).max()))
    save(fig, out, 'fig01_mcc_surface')


# ------------------------------------------------------------------------------------------- Fig. 2
def fig_meridian(run, out):
    """Fig. 2: closest-point projection in the meridian plane, linear elasticity (fig_meridian of fig_surface.py)."""
    cases = [('(a) subcritical region: compaction and hardening', 'fig2a_subcritical'),
             ('(b) supercritical region: dilation and softening', 'fig2b_supercritical')]
    csl = read_csv(os.path.join(run, 'fig2_critical_state_line.csv'))
    fig, axs = plt.subplots(1, 2, figsize=(TEXTW, 2.75))
    for ax, (title, stem) in zip(axs, cases):
        cv = read_csv(os.path.join(run, stem + '_curves.csv'))
        pts = read_csv(os.path.join(run, stem + '_points.csv'))
        (ptr, pp, ptip), (qtr, qq, qtip) = pts['p'][:3], pts['q'][:3]       # trial, projected, arrow tip
        an, a = cv['p_start'].max() / 2, cv['p_end'].max() / 2            # the ellipses run from p' = 0 to 2a
        ax.plot(cv['p_start'], cv['q_start'], color=MUTED, lw=1.0, label='$\\Phi_n = 0$ (start of step)')
        ax.plot(cv['p_end'], cv['q_end'], color=C1, lw=1.5, label='$\\Phi = 0$ (end of step)')
        ax.plot(csl['p'], csl['q'], color=INK2, lw=0.8, ls=(0, (5, 2)), label='critical state line $q = Mp\'$')
        ax.plot(cv['p_energy'], cv['q_energy'], color=C2, lw=0.9, ls=(0, (2, 1.5)),
                label='energy-norm contour through $\\sigma_{proj}$')
        ax.plot([ptr, pp], [qtr, qq], color=INK, lw=0.7)
        ax.plot(ptr, qtr, 'o', ms=5, mfc='white', mec=INK, mew=0.9, zorder=5)
        ax.plot(pp, qq, 'o', ms=5, mfc=C1, mec='white', mew=0.8, zorder=5)
        ax.annotate('trial', (ptr, qtr), xytext=(6, 4), textcoords='offset points', fontsize=7.5, color=INK)
        ax.annotate('projected', (pp, qq), xytext=(-8, -11), textcoords='offset points', fontsize=7.5, color=INK,
                    ha='right')
        # normal (flow direction) at the projected state, 28 kPa long
        ax.annotate('', xy=(ptip, qtip), xytext=(pp, qq),
                    arrowprops=dict(arrowstyle='-|>', color=C1, lw=0.9, mutation_scale=7))
        ax.annotate(r'$\partial\Phi/\partial\sigma$', (ptip, qtip), xytext=(4, -2),
                    textcoords='offset points', fontsize=7, color=C1)
        # what happens to the surface in the step
        grow = a > an
        txt = (('surface expands:' + '\n' + r'$\Delta\alpha > 0$, $a = %.1f > a_n = %.0f$') % (a, an) if grow else
               ('surface contracts:' + '\n' + r'$\Delta\alpha < 0$, $a = %.1f < a_n = %.0f$') % (a, an))
        ax.text(185, 8, txt, fontsize=6.8, color=INK2, ha='right', va='bottom')
        ax.set_xlim(0, 240); ax.set_ylim(0, 175)
        ax.set_xlabel("$p' = -\\xi/\\sqrt{3}$ (kPa)")
        ax.set_ylabel('$q = \\sqrt{3/2}\\,\\rho$ (kPa)')
        ax.set_aspect('equal')
        ax.set_title(title, fontsize=8.5)
        print('Fig. 2 %s: trial (%.6g, %.6g), projected (%.12g, %.12g), a_n = %.6g, a = %.12g'
              % (title[:3], ptr, qtr, pp, qq, an, a))
    axs[0].legend(loc='upper left', fontsize=6.8, handlelength=2.2)
    fig.tight_layout(w_pad=1.5)
    save(fig, out, 'fig02_meridian_projection')


if __name__ == '__main__':
    args = arguments('Figs. 1 and 2 of the article from the CSV files of YieldSurfaceProjection.')
    fig_surface(args.rundir, args.outdir)
    fig_meridian(args.rundir, args.outdir)
