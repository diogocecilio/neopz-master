"""Figure 5 of the article (error against work of the integration schemes) from the CSV files written by
IntegrationSchemes.

Usage (run the executable first; it writes its CSV files to the current directory):

    python3 <neopz>/Projects/IntegrationSchemes/plot_figures.py [run directory] [-o output directory]

The run directory (default: the current directory) holds the CSV files; the figures are written as PDF and PNG
to <run directory>/figures, or to the output directory.

Figures produced:

- fig06_integration_schemes: Fig. 6 (fig_rivais of figs.py of the Python code), error against work (local Newton
  iterations or evaluations of the elastoplastic operator) in the tests of Table 4: (a) test B of Xie et al.
  (relative error of the stress; BE with 1 to 1024 increments, RK in one increment with tolerances 1e-1 to 1e-8);
  (b) undrained test with OCR = 10 of Krabbenhoft and Lyamin (|p' - p'_exact| at eps_a = 8%; BE with 5 to 1000
  increments, RK in 10 increments with tolerances 1e-2 to 1e-7); (c) drained test on the normally consolidated
  clay (|q - q_exact| at eps_a = 25%; BE with 10 to 500 increments, RK with 10 to 100 increments, tolerance 1e-4).
  File: schemes_fig05.csv.
- supplementary_secant_first_increment (not in the article): first increment of the drained test with OCR = 10 and
  10 increments, a single backward-Euler step for lateral strains eps_r in [0.010, 0.030]: sigma_r + p'_0 and the
  local iterations with this work and with the secant shear modulus (failures and spurious states p' = 0 marked).
  File: schemes_secant_first_increment_sweep.csv.
"""
import csv
import os
import sys

import numpy as np
from matplotlib.transforms import blended_transform_factory

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'Common'))
from mcc_figstyle import C1, C2, C3, C4, INK2, MUTED, TEXTW, arguments, plt, save  # noqa: E402

# style of each scheme (fig_rivais of figs.py)
STY = {'exact_n': dict(color=C1, marker='o', ls='-', mfc=C1, ms=3.0, label='This work: BE, exact $K$, $G(p_n)$',
                       zorder=4),
       'exact_secant': dict(color=C2, marker='s', ls='--', mfc='white', ms=4.6, label='BE, exact $K$, secant $G$',
                            zorder=3),
       'frozen_n': dict(color=C3, marker='^', ls=':', mfc='white', ms=3.8, label='BE, frozen $K$', zorder=3),
       'ME2(1)': dict(color=C4, marker='D', ls='-.', mfc='white', ms=3.2, label='Explicit RK, ME2(1)', zorder=3),
       'RKDP5(4)': dict(color=INK2, marker='p', ls='-', mfc='white', ms=3.8, label='Explicit RK, RKDP5(4)',
                        zorder=3)}


def read_rows(path):
    """Rows of a CSV file with text columns (list of dicts)."""
    if not os.path.exists(path):
        sys.exit(f'file not found: {path} (run the executable in this directory first)')
    with open(path, newline='') as f:
        return list(csv.DictReader(f))


def curve(ax, pts, key):
    """pts = [(parameter, work, error)], joined in the order of the parameter (increments or -tolerance)."""
    pts = sorted((p for p in pts if np.isfinite(p[2])), key=lambda t: t[0])
    if pts:
        ax.plot([p[1] for p in pts], [p[2] for p in pts], lw=1.0, mew=0.9, **STY[key])


def fig_schemes(run, out):
    """Fig. 6: error against work in the tests of Table 4."""
    rows = read_rows(os.path.join(run, 'schemes_fig05.csv'))
    fig, axs = plt.subplots(1, 3, figsize=(TEXTW, 2.45))
    labels = {'a': ('relative error of $\\sigma$', '(a) undrained, NC (Xie et al. test B)'),
              'b': ("$|p' - p'_{\\mathrm{exact}}|$ at $\\varepsilon_a = 8\\%$ (kPa)",
                    '(b) undrained, OCR = 10 (K&L example 1)'),
              'c': ('$|q - q_{\\mathrm{exact}}|$ at $\\varepsilon_a = 25\\%$ (kPa)', '(c) drained, NC (K&L example 1)')}
    for ax, panel in zip(axs, 'abc'):
        sel = [r for r in rows if r['panel'] == panel and r['converged'] == '1']
        for key in ('exact_secant', 'exact_n', 'frozen_n'):
            curve(ax, [(int(r['n']), int(r['work']), float(r['error'])) for r in sel if r['scheme'] == key], key)
        for key in ('ME2(1)', 'RKDP5(4)'):
            # (a), (b): the curve runs over the tolerances (one or 10 increments); (c): over the increments
            par = (lambda r: int(r['n'])) if panel == 'c' else (lambda r: -float(r['tol']))
            curve(ax, [(par(r), int(r['work']), float(r['error'])) for r in sel if r['scheme'] == key], key)
        ax.set_ylabel(labels[panel][0])
        ax.set_title(labels[panel][1], fontsize=7.6)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel('local evaluations')
        failed = [r for r in rows if r['panel'] == panel and r['converged'] == '0']
        for r in failed:
            print(f"panel ({panel}): {r['scheme']} with n = {r['n']} has no solution (not drawn)")
    h, lab = axs[0].get_legend_handles_labels()
    order = [1, 0, 2, 3, 4]
    fig.legend([h[i] for i in order], [lab[i] for i in order], loc='upper center', ncol=5, fontsize=6.6,
               bbox_to_anchor=(0.5, 1.05), handlelength=2.4, columnspacing=1.0)
    fig.tight_layout(w_pad=0.6)
    save(fig, out, 'fig06_integration_schemes')


def fig_secant_sweep(run, out):
    """Supplementary figure: single step in the first increment of the drained test with OCR = 10."""
    path = os.path.join(run, 'schemes_secant_first_increment_sweep.csv')
    if not os.path.exists(path):
        return
    rows = read_rows(path)
    fig, axs = plt.subplots(1, 2, figsize=(TEXTW, 2.5))
    rug = blended_transform_factory(axs[0].transData, axs[0].transAxes)
    ref = {}
    for key, color, lab in (('exact_n', C1, 'This work: BE, exact $K$, $G(p_n)$'),
                            ('exact_secant_single_step', C2, 'BE, exact $K$, secant $G$, single step')):
        sel = [r for r in rows if r['scheme'] == key]
        x = np.array([float(r['eps_r']) for r in sel])
        conv = np.array([r['converged'] == '1' for r in sel])
        spur = np.array([r['spurious_p_zero'] == '1' for r in sel])
        ok = conv & ~spur
        sr = np.array([float(r['sigma_r_plus_p0']) for r in sel])
        its = np.array([int(r['iterations']) for r in sel])
        axs[0].plot(100 * x[ok], sr[ok], 'o', ms=1.6, mec='none', color=color, label=lab)
        axs[1].plot(100 * x[conv], its[conv], 'o', ms=1.6, mec='none', color=color, label=lab)
        if key == 'exact_n':
            ref = dict(zip(x, sr))
            print(f'{key}: converged {ok.sum()} of {len(sel)}, {its.min()} to {its.max()} local iterations')
            continue
        off = np.array([o and abs(v - ref.get(xx, v)) > 50.0 for o, v, xx in zip(ok, sr, x)])
        axs[0].plot(100 * x[spur], sr[spur], 'x', color=C3, ms=3.5, mew=0.8, label="secant $G$: spurious state $p'=0$")
        axs[0].plot(100 * x[~conv], np.full((~conv).sum(), 0.03), '|', color=MUTED, ms=5, mew=0.6, transform=rug,
                    label='secant $G$: no convergence (50 iterations)')
        print(f'{key}: converged {ok.sum()} (of which {off.sum()} on another branch, more than 50 kPa from this '
              f'work), spurious p = 0 {spur.sum()}, failed {(~conv).sum()} of {len(sel)}')
    axs[0].axhline(0.0, color=INK2, lw=0.6)
    for ax in axs:
        ax.set_xlabel('lateral strain $\\varepsilon_r$ of the increment (%)')
    axs[0].set_ylabel("$\\sigma_r + p'_0$ (kPa)")
    axs[0].set_title('(a) radial stress after one step', fontsize=7.6)
    axs[1].set_ylabel('local Newton iterations')
    axs[1].set_title('(b) local iterations of the converged steps', fontsize=7.6)
    h, lab = axs[0].get_legend_handles_labels()
    fig.legend(h, lab, loc='upper center', ncol=2, fontsize=6.4, bbox_to_anchor=(0.5, 1.12), markerscale=1.6)
    fig.tight_layout(w_pad=1.0)
    save(fig, out, 'supplementary_secant_first_increment')


if __name__ == '__main__':
    args = arguments('Fig. 6 of the article (integration schemes) from the CSV files of IntegrationSchemes.')
    fig_schemes(args.rundir, args.outdir)
    fig_secant_sweep(args.rundir, args.outdir)
