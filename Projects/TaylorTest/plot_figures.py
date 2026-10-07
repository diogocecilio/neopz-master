"""Figure 3 of the article (Taylor test of the consistent tangent) from the CSV files written by TaylorTest.

Usage (run the executable first; it writes its CSV files to the current directory):

    python3 <neopz>/Projects/TaylorTest/plot_figures.py [run directory] [-o output directory]

The run directory (default: the current directory) holds the CSV files; the figures are written as PDF and PNG
to <run directory>/figures, or to the output directory.

Figures produced (same panels, fits, slope triangles and annotations as fig_taylor of figs.py of the Python code):

- fig03_taylor_test: Fig. 3, log E(alpha) against log alpha for 300 random perturbations at the subcritical
  (compaction) and supercritical (dilation) states, with the consistent tangent D (panels a, b: second order)
  and with its transpose D^T (panels c, d: first order). Files: taylor_pcg64_<kind>_<op>.csv and
  taylor_pcg64_summary.csv (numpy-compatible PCG64 stream, the draws of gen_data.py, i.e. the data of Fig. 3).
- fig03_taylor_test_mt19937: the same figure for the independent std::mt19937_64 sample
  (taylor_mt19937_*.csv; not in the article, other states, same orders). Written only when these files exist.
- supplementary_taylor_operators: Taylor test of the tangent operators of Sect. 6.7 at the two plastic states of
  Fig. 3: consistent tangent D, transpose D^T, symmetric part (D + D^T)/2 and continuum operator (the slopes are
  the last column of Table 10; not a figure of the article). Files: taylor_operators_<kind>_<op>.csv and
  taylor_operators_summary.csv (draws of numpy's default_rng(7), as tangentes() of gen_data.py).
"""
import os
import sys

import numpy as np
from matplotlib.patches import Polygon as MPoly

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'Common'))
from mcc_figstyle import C1, C2, INK, INK2, TEXTW, arguments, plt, read_csv, save  # noqa: E402

KINDS = {1: 'subcritical', 2: 'supercritical'}                    # kind index of the summary -> file name
NAMES = {1: 'subcritical (compaction)', 2: 'supercritical (dilation)'}


def summary_row(S, kind, transp):
    """Row of taylor_<stream>_summary.csv of the state kind (0 elastic, 1 subcritical, 2 supercritical)."""
    i = np.flatnonzero((S['kind'] == kind) & (S['transposed'] == int(transp)))
    if i.size != 1:
        sys.exit(f'summary file: no single row for kind {kind}, transposed {int(transp)}')
    return {k: v[i[0]] for k, v in S.items()}


def fig_taylor(run, out, stream='pcg64', name='fig03_taylor_test'):
    """Fig. 3: Taylor test of D (second order) and of D^T (first order) at the two plastic states."""
    S = read_csv(os.path.join(run, f'taylor_{stream}_summary.csv'))
    fig, axs = plt.subplots(2, 2, figsize=(TEXTW * 0.86, 4.6))
    for row, transp in enumerate((False, True)):
        for col, kind in enumerate((1, 2)):
            ax = axs[row, col]
            r = summary_row(S, kind, transp)
            P = read_csv(os.path.join(run, f'taylor_{stream}_{KINDS[kind]}_{"DT" if transp else "D"}.csv'))
            x, y = P['log_alpha1'], P['log_E1']
            b1, b0 = r['fit_slope'], r['fit_intercept']
            ax.plot(x, y, 'o', ms=1.9, mfc=C1, mec='none', alpha=0.55)
            xx = np.array([x.min() - 0.1, x.max() + 0.1])
            ax.plot(xx, b0 + b1 * xx, color=C2, lw=1.3)
            x0 = x.min() + 0.5
            p0, p1, p2 = (x0, b0 + b1 * x0), (x0 + 1, b0 + b1 * (x0 + 1)), (x0 + 1, b0 + b1 * x0)
            ax.add_patch(MPoly([p0, p1, p2], closed=True, facecolor='#d9d8d2', edgecolor=INK, lw=0.7, alpha=0.9))
            ax.text(p2[0] + 0.15, p2[1] + 0.3 * b1, f'slope = {b1:.3f}', fontsize=7.5, color=INK, va='center',
                    bbox=dict(boxstyle='square,pad=0.15', fc='white', ec='none'))
            tg = '$\\mathbb{D}$' if not transp else '$\\mathbb{D}^{\\mathsf{T}}$'
            ax.set_title('(' + 'abcd'[2 * row + col] + ') ' + NAMES[kind] + ', ' + tg, fontsize=8)
            info = "$p' = %.1f$ kPa, $q = %.1f$ kPa" % (r['p_eff'], r['q'])
            info += ('\n' + 'consistent operator: second order' if not transp else
                     '\n' + r'$\|\mathbb{D}-\mathbb{D}^{\mathsf{T}}\|/\|\mathbb{D}\| = %.1f\%%$: first order'
                     % (100 * r['asym']))
            ax.text(0.03, 0.96, info, transform=ax.transAxes, fontsize=6.8, color=INK2, va='top', ha='left',
                    bbox=dict(boxstyle='square,pad=0.2', fc='white', ec='none', alpha=0.85))
            ax.set_xlabel('log $\\alpha$')
            if col == 0:
                ax.set_ylabel('log $E(\\alpha)$')
            print('Taylor %s %-13s %-3s fit %.3f, median pair slope %.4f, asym %.4f, points %d'
                  % (stream, KINDS[kind], 'D^T' if transp else 'D', b1, np.median(P['pair_slope']), r['asym'],
                     x.size))
    fig.tight_layout(h_pad=1.0, w_pad=1.0)
    save(fig, out, name)


OPERATORS = (('D', 0, '$\\mathbb{D}$'), ('DT', 1, '$\\mathbb{D}^{\\mathsf{T}}$'),
             ('sym', 2, '$(\\mathbb{D}+\\mathbb{D}^{\\mathsf{T}})/2$'), ('cont', 3, 'continuum'))


def fig_operators(run, out):
    """Supplementary figure: Taylor test of the four tangent operators of Table 10 at the states of Fig. 3."""
    path = os.path.join(run, 'taylor_operators_summary.csv')
    if not os.path.exists(path):
        return
    S = read_csv(path)
    fig, axs = plt.subplots(2, 4, figsize=(TEXTW, 3.5), sharex=True)
    for row, kind in enumerate((1, 2)):
        for col, (op, code, label) in enumerate(OPERATORS):
            ax = axs[row, col]
            i = np.flatnonzero((S['kind'] == kind) & (S['operator'] == code))
            if i.size != 1:
                sys.exit(f'operator summary: no single row for kind {kind}, operator {op}')
            b1, b0, med = S['fit_slope'][i[0]], S['fit_intercept'][i[0]], S['median_slope'][i[0]]
            P = read_csv(os.path.join(run, f'taylor_operators_{KINDS[kind]}_{op}.csv'))
            x, y = P['log_alpha1'], P['log_E1']
            ax.plot(x, y, 'o', ms=1.4, mfc=C1, mec='none', alpha=0.55)
            xx = np.array([x.min() - 0.1, x.max() + 0.1])
            ax.plot(xx, b0 + b1 * xx, color=C2, lw=1.1)
            ax.text(0.04, 0.95, f'fit {b1:.3f}\nmedian {med:.3f}', transform=ax.transAxes, fontsize=6.4, color=INK,
                    va='top', ha='left', bbox=dict(boxstyle='square,pad=0.15', fc='white', ec='none', alpha=0.85))
            ax.set_title('(' + 'abcdefgh'[4 * row + col] + ') ' + label, fontsize=7.4)
            if row == 1:
                ax.set_xlabel('log $\\alpha$')
            if col == 0:
                ax.set_ylabel(NAMES[kind].split(' ')[0] + '\nlog $E(\\alpha)$')
            ax.tick_params(labelsize=6.6)
            print('Taylor operator %-13s %-4s fit %.4f, median %.4f' % (KINDS[kind], op, b1, med))
    fig.tight_layout(h_pad=0.6, w_pad=0.5)
    save(fig, out, 'supplementary_taylor_operators')


if __name__ == '__main__':
    args = arguments('Fig. 3 of the article (Taylor test) from the CSV files of TaylorTest.')
    fig_taylor(args.rundir, args.outdir)
    if os.path.exists(os.path.join(args.rundir, 'taylor_mt19937_summary.csv')):
        fig_taylor(args.rundir, args.outdir, 'mt19937', 'fig03_taylor_test_mt19937')
    fig_operators(args.rundir, args.outdir)
