"""Figure 6 of the article (Terzaghi consolidation) from the CSV files written by TerzaghiConsolidation.

Usage (run the executable first; it writes its CSV files to the current directory):

    python3 <neopz>/Projects/TerzaghiConsolidation/plot_figures.py [run directory] [-o output directory]

The run directory (default: the current directory) holds the CSV files; the figures are written as PDF and PNG
to <run directory>/figures, or to the output directory.

Figure produced (same panels, curves, markers and annotations as fig_terzaghi of figs.py of the Python code):

- fig06_terzaghi_consolidation: Fig. 6, (a) excess pore pressure p_w/q along the column at x = 0 for
  T = 0.001, 0.01, 0.1 and 0.5, Q8-Q4 at the vertices (terzaghi_fig6a.csv) against the exact series solution;
  (b) degree of consolidation U = w/w_inf of the top settlement against T = c_v t/H^2 (terzaghi_fig6b.csv),
  Q8-Q4 against the exact series solution. The time and settlement of terzaghi_q8q4_history.csv give c_v and
  w_inf of the annotation; terzaghi_hex20hex8_history.csv is checked to be identical to Q8-Q4 (a warning is
  printed when it differs or is missing).
"""
import os
import sys

import numpy as np
from matplotlib.lines import Line2D

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'Common'))
from mcc_figstyle import C1, C2, C3, C4, INK, INK2, TEXTW, arguments, plt, read_csv, save  # noqa: E402

NTERMS = 300                                       # terms of the series of the exact solutions (figs.py)
A = (2 * np.arange(NTERMS) + 1) * np.pi / 2


def pw_exact(y, H, T):
    """Exact excess pore pressure p_w/q at the height y (drained top y = H, impermeable base y = 0)."""
    return np.array([np.sum(2 / A * np.sin(A * (H - yi) / H) * np.exp(-A * A * T)) for yi in y])


def u_exact(T):
    """Exact degree of consolidation U(T)."""
    return np.array([1 - np.sum(2 / A ** 2 * np.exp(-A * A * t)) for t in T])


def consolidation_data(run, H):
    """c_v (m^2/s) and w_inf = qH/E_oed (m) from the history (t, settlement) and terzaghi_fig6b.csv (T, U).

    Also compares the settlement history of the Hex20-Hex8 model with that of the Q8-Q4 model.
    """
    h = read_csv(os.path.join(run, 'terzaghi_q8q4_history.csv'))
    b = read_csv(os.path.join(run, 'terzaghi_fig6b.csv'))
    if h['t'].size != b['T'].size:
        sys.exit('terzaghi_q8q4_history.csv and terzaghi_fig6b.csv have different numbers of rows')
    k = int(np.argmax(h['t']))
    cv = b['T'][k] * H ** 2 / h['t'][k]
    winf = h['settlement'][-1] / b['U_numerical'][-1]
    hex_file = os.path.join(run, 'terzaghi_hex20hex8_history.csv')
    if not os.path.exists(hex_file):
        print('warning: terzaghi_hex20hex8_history.csv not found; the legend "Hex20-Hex8 identical" is not checked')
    else:
        w3 = read_csv(hex_file)['settlement']
        if w3.size != h['settlement'].size:
            print('warning: the Hex20-Hex8 and Q8-Q4 histories have %d and %d rows (legend "Hex20-Hex8 identical")'
                  % (w3.size, h['settlement'].size))
        else:
            dw = np.abs(w3 - h['settlement']).max()
            print('Hex20-Hex8 against Q8-Q4: max |settlement difference| = %.2e m' % dw)
            if dw > 1e-9 * np.abs(h['settlement']).max():
                print('warning: the Hex20-Hex8 history differs from Q8-Q4 (legend "Hex20-Hex8 identical")')
    return cv, winf


def fig_terzaghi(run, out):
    """Fig. 6: pore pressure isochrones at x = 0 and degree of consolidation."""
    a6 = read_csv(os.path.join(run, 'terzaghi_fig6a.csv'))
    b6 = read_csv(os.path.join(run, 'terzaghi_fig6b.csv'))
    H = a6['y'].max()
    cv, winf = consolidation_data(run, H)
    fig, axs = plt.subplots(1, 2, figsize=(TEXTW, 2.55), gridspec_kw=dict(width_ratios=(1, 1.25)))
    ax = axs[0]
    yy = np.linspace(0, H, 200)
    cols = [C1, C2, C3, C4]
    rows = []
    for k, T in enumerate(np.unique(a6['T'])):
        m = a6['T'] == T
        o = np.argsort(a6['y'][m])
        ys, p = a6['y'][m][o], a6['pw_over_q'][m][o]
        rows.append((T, ys, p))
        ax.plot(pw_exact(yy, H, T), yy, color=cols[k], lw=1.1)
        ax.plot(p, ys, 'o', ms=3.6, mfc='white', mec=cols[k], mew=0.9)
        print('T = %-5g max |p_w/q - exact| at the vertices = %.5f' % (T, np.abs(p - a6['exact'][m][o]).max()))
    ax.set(xlabel='$p_w/q$', ylabel='$y$ (m)', xlim=(0, 1.62), ylim=(0, H))
    _, ys1, p1 = rows[0]
    ax.annotate('oscillation next to\nthe drained face', (p1[-2], ys1[-2]), xytext=(1.12, 8.6), fontsize=6.3,
                color=INK2, arrowprops=dict(arrowstyle='-', color=INK2, lw=0.6))
    ax.text(1.06, 0.35, 'drained top: $p_w = 0$\nimpermeable base', fontsize=6.3, color=INK2, va='bottom',
            linespacing=1.4)
    ax.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_title('(a) excess pore pressure at $x = 0$', fontsize=8)
    handles = [Line2D([], [], color=cols[k], lw=1.1, marker='o', ms=3.6, mfc='white', mec=cols[k],
                      label=f'$T$ = {r[0]:g}') for k, r in enumerate(rows)]
    handles += [Line2D([], [], color=INK2, lw=1.1, label='exact'),
                Line2D([], [], ls='', marker='o', ms=3.6, mfc='white', mec=INK2, label='Q8–Q4 (vertices)')]
    ax.legend(handles=handles, loc='center left', bbox_to_anchor=(0.665, 0.42), fontsize=6.8, handlelength=1.6,
              labelspacing=0.3, borderaxespad=0.0)
    ax = axs[1]
    Tt, Un = b6['T'][1:], b6['U_numerical'][1:]       # the first row is the initial state (before the load)
    Tf = np.logspace(-5, 0, 300)
    ax.semilogx(Tf, u_exact(Tf), color=INK, lw=0.9, ls=(0, (4, 2)), label='exact')
    ax.semilogx(Tt, Un, color=C1, lw=1.4, label='Q8–Q4 (Hex20–Hex8 identical)')
    ax.set(xlabel='$T = c_v t / H^2$', ylabel='$U = w/w_\\infty$', xlim=(1e-5, 1), ylim=(0, 1))
    ax.text(0.03, 0.55, '$c_v = kE_{oed}$ = %.3g m$^2$/s, $H$ = %g m\n$w_\\infty = qH/E_{oed}$ = %.2f mm'
            % (cv, H, 1e3 * winf), fontsize=6.5, color=INK2, transform=ax.transAxes)
    Ti = 1e-4
    ax.annotate('offset at small $T$: diffusion layer\nthinner than one element', (Ti, np.interp(Ti, Tt, Un)),
                xytext=(2e-5, 0.22), fontsize=6.3, color=INK2, arrowprops=dict(arrowstyle='-', color=INK2, lw=0.6))
    ax.set_title('(b) degree of consolidation (top settlement)', fontsize=8)
    ax.legend(loc='upper left', fontsize=7)
    fig.tight_layout(w_pad=1.2)
    save(fig, out, 'fig06_terzaghi_consolidation')


if __name__ == '__main__':
    args = arguments('Fig. 6 of the article (Terzaghi consolidation) from the CSV files of TerzaghiConsolidation.')
    fig_terzaghi(args.rundir, args.outdir)
