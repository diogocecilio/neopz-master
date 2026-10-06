"""Supplementary figure of Table 3 of the article (exact versus frozen integration of the porous law) from the CSV
files written by FrozenBulkModulus.

Usage (run the executable first; it writes its CSV files to the current directory):

    python3 <neopz>/Projects/FrozenBulkModulus/plot_figures.py [run directory] [-o output directory]

The run directory (default: the current directory) holds the CSV files; the figures are written as PDF and PNG
to <run directory>/figures, or to the output directory.

The article has no figure for this example: Sect. 6.1 reports its results in Table 3 and in the text. The figure
below is supplementary (NOT in the article); it is labelled as such and draws the numbers of Table 3 with the
style of the figures of the article (Common/mcc_figstyle.py):

- supplementary_table03_frozen_bulk_modulus: triaxial tests at a material point with the RS2 clay (G = 20 MPa),
  eps_a up to 20% in n = 10 ... 1600 increments, NC (p'0 = p'c0 = 200 kPa) and OCR = 5 (p'0 = 100 kPa,
  p'c0 = 500 kPa); exact integration of the porous law (this work) against the bulk modulus frozen at its trial
  value during the plastic correction (Sanei et al. 2020). Panels: (a) drained, largest error in q along the
  path against the closed form of Appendix B.1; (b) drained, mean number of local Newton iterations per plastic
  projection; (c) undrained, largest error in p' against the closed-form path (B.7); (d) undrained, mean local
  iterations. Table 3 holds the drained NC curves of (a) and (b) and the frozen curves of (c) (rows n = 15 and 25
  come from the Python code, not from the article); the other curves are quoted in the text of Sect. 6.1 or
  printed by the executable. The drained NC tests of the frozen form with
  n = 10 and 15 have no solution (dashes of Table 3): they are marked with crosses at the bottom of (a) and (b).
  File: frozen_runs.csv (one row per test: test, state, integration, n, converged, errors, local iterations).
"""
import os
import sys

import numpy as np
from matplotlib.lines import Line2D
from matplotlib.transforms import blended_transform_factory

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'Common'))
from mcc_figstyle import C1, C3, INK2, MUTED, TEXTW, arguments, plt, read_csv, save  # noqa: E402

# column names of frozen_runs.csv (the codes of the first three columns are given in the header)
TEST, STATE, INTEG = 'test(0=drained;1=undrained)', 'state(0=NC;1=OCR5)', 'integration(0=exact;1=frozen)'
INTEGRATIONS = ((0, 'exact', C1, 'exact integration (this work)'), (1, 'frozen', C3, 'frozen bulk modulus'))
STATES = ((0, 'NC', 'o', '-'), (1, 'OCR = 5', 's', (0, (4, 2))))


def series(R, test, state, integ, column):
    """Values of a column of frozen_runs.csv for one test, state and integration: n, value, converged."""
    i = np.flatnonzero((R[TEST] == test) & (R[STATE] == state) & (R[INTEG] == integ))
    i = i[np.argsort(R['n'][i])]
    return R['n'][i].astype(int), R[column][i], R['converged'][i] > 0.5


def style(integ, state):
    """Line and marker of a curve: colour by integration, marker and dashes by initial state. The exact curves
    (filled markers) are drawn on top of the frozen ones (larger open markers), so that both stay visible where
    they coincide."""
    c = INTEGRATIONS[integ][2]
    _, _, mk, ls = STATES[state]
    if integ == 0:
        return dict(color=c, ls=ls, lw=1.0, marker=mk, ms=2.8, mec=c, mew=0.6, mfc=c, zorder=3)
    return dict(color=c, ls=ls, lw=1.0, marker=mk, ms=4.8, mec=c, mew=0.8, mfc='white', zorder=2)


def fig_table3(run, out):
    """Supplementary figure (not in the article): the data of Table 3 against the number of increments."""
    R = read_csv(os.path.join(run, 'frozen_runs.csv'))
    fig, axs = plt.subplots(2, 2, figsize=(TEXTW, 4.7))
    panels = (('drained', 0, 'max_error', "largest error in $q$ (kPa)"),
              ('drained', 0, 'mean_local_its', 'mean local iterations'),
              ('undrained', 1, 'max_error', "largest error in $p'$ (kPa)"),
              ('undrained', 1, 'mean_local_its', 'mean local iterations'))
    for k, (ax, (name, test, column, ylabel)) in enumerate(zip(axs.flat, panels)):
        failed = []
        for integ, _, _, _ in INTEGRATIONS:
            for state, _, _, _ in STATES:
                n, v, ok = series(R, test, state, integ, column)
                ax.plot(n[ok], v[ok], **style(integ, state))
                failed += [(m, integ, state) for m in n[~ok]]
        if failed:                                    # tests without solution (dashes of Table 3)
            tr = blended_transform_factory(ax.transData, ax.transAxes)
            for m, integ, state in failed:
                ax.plot(m, 0.05, 'x', ms=4.5, mew=1.0, color=INTEGRATIONS[integ][2], transform=tr)
            m = max(f[0] for f in failed)
            who = sorted({(INTEGRATIONS[i][1], STATES[s][1]) for _, i, s in failed})
            ax.text(m * 1.25, 0.05, 'no solution: ' + ', '.join(f'{a}, {b}' for a, b in who) + f', $n \\leq {m}$',
                    transform=tr, fontsize=6.5, color=INK2, va='center', ha='left')
        ax.set_xscale('log')
        ax.set_xticks([10, 20, 50, 100, 200, 400, 800, 1600])
        ax.set_xticklabels(['10', '20', '50', '100', '200', '400', '800', '1600'])
        ax.minorticks_off()
        ax.set_xlim(8, 2000)
        if column == 'max_error':
            ax.set_yscale('log')
        ax.set_title(f'({"abcd"[k]}) {name}', fontsize=8)
        ax.set_ylabel(ylabel)
        if k >= 2:
            ax.set_xlabel('number of increments $n$')
    # slope -1 (first order) next to the drained curves
    ax = axs[0, 0]
    x0, y0 = 320., 30.
    ax.plot([x0, 2 * x0, x0, x0], [y0, y0 / 2, y0 / 2, y0], color=MUTED, lw=0.7)
    ax.text(2.25 * x0, y0 / 1.6, 'slope $-1$', fontsize=6.5, color=INK2, ha='left', va='center')
    # undrained: the exact integration follows the closed-form path (B.7) to round-off
    big = max(v[ok].max() for _, v, ok in (series(R, 1, s, 0, 'max_error') for s, _, _, _ in STATES))
    axs[1, 0].text(0.97, 0.42, 'exact integration: round-off\n(largest %.1e kPa)' % big,
                   transform=axs[1, 0].transAxes, fontsize=6.5, color=INK2, ha='right', va='center')
    axs[1, 1].text(0.97, 0.94, 'exact and frozen: almost the same iterations', transform=axs[1, 1].transAxes,
                   fontsize=6.5, color=INK2, ha='right', va='top')
    handles = [Line2D([], [], **style(i, s)) for i, _, _, _ in INTEGRATIONS for s, _, _, _ in STATES]
    labels = [f'{lab}, {st}' for _, _, _, lab in INTEGRATIONS for _, st, _, _ in STATES]
    fig.legend(handles, labels, loc='lower center', ncol=2, bbox_to_anchor=(0.5, -0.01), fontsize=7.3,
               handlelength=2.8, columnspacing=3.0)
    fig.text(0.0, 1.0, 'Supplementary figure, NOT in the article: Table 3 (Sect. 6.1) and the other tests of the '
             'example FrozenBulkModulus.\nRS2 clay, triaxial tests at a material point to $\\varepsilon_a$ = 20% in '
             '$n$ increments; exact integration of the porous law (this work)\nand bulk modulus frozen at its '
             'trial value during the plastic correction (Sanei et al. 2020).',
             fontsize=7.3, color=INK2, ha='left', va='bottom', linespacing=1.3)
    fig.tight_layout(rect=(0, 0.075, 1, 0.97), h_pad=1.2, w_pad=1.5)
    save(fig, out, 'supplementary_table03_frozen_bulk_modulus')
    for state, st, _, _ in STATES:                  # ratios quoted in the text of Sect. 6.1
        e, f = series(R, 0, state, 0, 'max_error'), series(R, 0, state, 1, 'max_error')
        r = f[1][e[2] & f[2]] / e[1][e[2] & f[2]]
        print('drained %-7s largest error frozen / exact: %.2f to %.2f' % (st, r.min(), r.max()))
    print('undrained, exact integration: largest error in p\' %.1e kPa' % big)


if __name__ == '__main__':
    args = arguments('Supplementary figure of Table 3 (not in the article) from the CSV files of '
                     'FrozenBulkModulus.')
    fig_table3(args.rundir, args.outdir)
