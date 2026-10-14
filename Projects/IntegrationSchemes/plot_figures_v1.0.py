"""Figures of the comparison of integration schemes (Sect. 5.2 of the article, v1.0) from the CSV files written by
IntegrationSchemes, and the legacy figure of v0.7 (error against work of the schemes of Table 4).

Usage (run the executable first; it writes its CSV files to the current directory):

    python3 <neopz>/Projects/IntegrationSchemes/plot_figures.py [run directory] [-o output directory]

The run directory (default: the current directory) holds the CSV files; the figures are written as PDF and PNG
to <run directory>/figures, or to the output directory.

Figures produced:

- fig_error_work (article, v1.0): error against work (local Newton iterations of the implicit schemes, evaluations
  of the elastoplastic operator of the explicit ones) of every scheme in four tests: (a) test B of Xie et al.
  (undrained, NC, one increment of 10%); (b) undrained, OCR = 10, 10 increments; (c) drained, NC; (d) drained,
  OCR = 10 (softening). Implicit schemes: the number of increments varies along the curve; adaptive RK: the
  tolerance (a, b) or the increments (c, d); fixed-step explicit schemes: the sub-steps of the increments of the
  test. File schemes_error_work.csv.
- fig_order (article, v1.0): error against the number of increments in the strain-controlled undrained tests
  (test B and OCR = 10), with the slopes 1, 2 and 4. File schemes_error_work.csv.
- fig_local_convergence (article, v1.0): norm of the residual against the local Newton iteration in four trial
  states for the return mapping with two unknowns (this work), with four unknowns, with the secant shear modulus
  and for the trapezoidal rule. File schemes_local_convergence.csv.
- fig06_integration_schemes (v0.7, Fig. 6): error against work of the schemes of Table 4 of v0.7. File
  schemes_fig05.csv.
- supplementary_secant_first_increment (not in the article): first increment of the drained test with OCR = 10 and
  10 increments, a single step for lateral strains eps_r in [0.010, 0.030]: sigma_r + p'_0 and the local iterations
  with this work and with the secant shear modulus (failures and spurious states p' = 0 marked). File
  schemes_secant_first_increment_sweep.csv.
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

# styles of the figure of v1.0 (scheme, curve) -> style
ORANGE = '#d97a1f'
STY2 = {('library_exact', 'n'): dict(color=C1, marker='o', ls='-', mfc=C1, ms=3.0, zorder=6,
                                     label='BE in RHW, two unknowns (this work)'),
        ('library_full', 'n'): dict(color=C1, marker='o', ls='--', mfc='white', ms=3.4, zorder=5,
                                    label='BE in RHW, four unknowns'),
        ('exact_secant', 'n'): dict(color=C2, marker='s', ls='--', mfc='white', ms=3.6, zorder=4,
                                    label='BE, secant $G$ (sub-stepped)'),
        ('frozen_n', 'n'): dict(color=C3, marker='^', ls=':', mfc='white', ms=3.6, zorder=4, label='BE, frozen $K$'),
        ('trapezoidal', 'n'): dict(color=ORANGE, marker='v', ls='-', mfc='white', ms=3.6, zorder=4,
                                   label='trapezoidal rule, secant $G$'),
        ('FE', 'nsub'): dict(color=C4, marker='x', ls='-.', mfc=C4, ms=3.6, zorder=3,
                             label='explicit Euler, fixed sub-steps'),
        ('RK4', 'nsub'): dict(color=C4, marker='D', ls='-', mfc='white', ms=3.0, zorder=3,
                              label='RK4, fixed sub-steps'),
        ('ME2(1)', 'tol'): dict(color=INK2, marker='+', ls='-.', mfc=INK2, ms=4.2, zorder=3,
                                label='adaptive ME2(1)'),
        ('ME2(1)', 'n'): dict(color=INK2, marker='+', ls='-.', mfc=INK2, ms=4.2, zorder=3, label='adaptive ME2(1)'),
        ('RKDP5(4)', 'tol'): dict(color=INK2, marker='p', ls='-', mfc='white', ms=3.6, zorder=3,
                                  label='adaptive RKDP5(4)'),
        ('RKDP5(4)', 'n'): dict(color=INK2, marker='p', ls='-', mfc='white', ms=3.6, zorder=3,
                                label='adaptive RKDP5(4)')}
ORDER2 = [('library_exact', 'n'), ('library_full', 'n'), ('exact_secant', 'n'), ('frozen_n', 'n'),
          ('trapezoidal', 'n'), ('FE', 'nsub'), ('RK4', 'nsub'), ('ME2(1)', 'tol'), ('RKDP5(4)', 'tol')]


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


def curve2(ax, pts, sty):
    pts = sorted((p for p in pts if np.isfinite(p[2]) and p[2] > 0), key=lambda t: t[0])
    if pts:
        ax.plot([p[1] for p in pts], [p[2] for p in pts], lw=1.0, mew=0.9, **sty)


def fig_error_work(run, out):
    """Error against work of every scheme in four tests (article, v1.0)."""
    rows = read_rows(os.path.join(run, 'schemes_error_work.csv'))
    panels = [('xieB', '(a) undrained, NC, one increment of 10% (test B of Xie et al.)',
               'relative error of $\\sigma$'),
              ('kl_undrained_ocr10', '(b) undrained, OCR = 10, 10 increments to 8%',
               "$|p' - p'_{\\mathrm{exact}}|$ (kPa)"),
              ('kl_drained_ocr1', '(c) drained, NC, 10 increments to 25%', '$|q - q_{\\mathrm{exact}}|$ (kPa)'),
              ('kl_drained_ocr10', '(d) drained, OCR = 10, 10 increments to 25%', '$|q - q_{\\mathrm{exact}}|$ (kPa)')]
    fig, axs = plt.subplots(2, 2, figsize=(TEXTW, 5.0))
    for ax, (test, title, ylabel) in zip(axs.ravel(), panels):
        sel = [r for r in rows if r['test'] == test and r['converged'] == '1']
        drained = test.startswith('kl_drained')
        for key in ORDER2:
            scheme, kind = key
            if scheme in ('ME2(1)', 'RKDP5(4)') and drained:
                kind = 'n'
            pts = []
            for r in sel:
                if r['scheme'] != scheme or r['curve'] != kind:
                    continue
                par = int(r['n']) if kind == 'n' else (int(r['nsub']) if kind == 'nsub' else -float(r['tol']))
                pts.append((par, int(r['work']), float(r['error'])))
            # the fixed-step curves start at one sub-step (curve "n" at the increments of the test)
            if kind == 'nsub':
                narticle = 1 if test == 'xieB' else 10
                for r in sel:
                    if r['scheme'] == scheme and r['curve'] == 'n' and int(r['n']) == narticle:
                        pts.append((1, int(r['work']), float(r['error'])))
            curve2(ax, pts, STY2[key])
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel('local evaluations (Newton iterations or operator evaluations)')
        ax.set_ylabel(ylabel)
        ax.set_title(title, fontsize=7.6)
        for r in [r for r in rows if r['test'] == test and r['converged'] == '0']:
            print(f"{test}: {r['scheme']} (n = {r['n']}, nsub = {r['nsub']}, tol = {r['tol']}) has no solution")
    h, lab = axs[0, 0].get_legend_handles_labels()
    fig.legend(h, lab, loc='upper center', ncol=3, fontsize=6.4, bbox_to_anchor=(0.5, 1.0), handlelength=2.6,
               columnspacing=1.0)
    fig.tight_layout(w_pad=0.8, h_pad=0.9, rect=(0, 0, 1, 0.925))
    save(fig, out, 'fig_error_work')


def fig_order(run, out):
    """Error against the number of increments in the strain-controlled tests, with the reference slopes."""
    rows = read_rows(os.path.join(run, 'schemes_error_work.csv'))
    panels = [('xieB', '(a) undrained, NC (test B of Xie et al.), $\\varepsilon_a = 10\\%$',
               'relative error of $\\sigma$'),
              ('kl_undrained_ocr10', '(b) undrained, OCR = 10, $\\varepsilon_a = 8\\%$',
               "$|p' - p'_{\\mathrm{exact}}|$ (kPa)")]
    keys = [('library_exact', 'n'), ('exact_secant', 'n'), ('frozen_n', 'n'), ('trapezoidal', 'n'), ('FE', 'n'),
            ('RK4', 'n')]
    fig, axs = plt.subplots(1, 2, figsize=(TEXTW, 2.6))
    for ax, (test, title, ylabel) in zip(axs, panels):
        sel = [r for r in rows if r['test'] == test and r['converged'] == '1' and r['curve'] == 'n']
        for scheme, kind in keys:
            pts = [(int(r['n']), int(r['n']), float(r['error'])) for r in sel
                   if r['scheme'] == scheme and (r['nsub'] in ('0', '1'))]
            sty = dict(STY2.get((scheme, 'nsub' if scheme in ('FE', 'RK4') else 'n')))
            if scheme in ('FE', 'RK4'):
                sty['label'] = 'explicit Euler, one step per increment' if scheme == 'FE' else 'RK4, one step per increment'
            curve2(ax, pts, sty)
        # reference slopes
        ns = np.array([64.0, 1024.0]) if test == 'xieB' else np.array([50.0, 1000.0])
        y0 = {'xieB': 3e-3, 'kl_undrained_ocr10': 1e-1}[test]
        for order, dy in ((1, 1.0), (2, 1.0), (4, 1.0)):
            y = y0 * dy * (ns / ns[0]) ** (-order)
            ax.plot(ns, y, color=MUTED, lw=0.7)
            ax.annotate(f'{order}', (ns[1], y[1]), fontsize=6.4, color=MUTED, ha='left', va='center',
                        xytext=(2, 0), textcoords='offset points')
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel('number of increments $n$')
        ax.set_ylabel(ylabel)
        ax.set_title(title, fontsize=7.6)
    h, lab = axs[0].get_legend_handles_labels()
    fig.legend(h, lab, loc='upper center', ncol=3, fontsize=6.4, bbox_to_anchor=(0.5, 1.12), handlelength=2.6)
    fig.tight_layout(w_pad=1.0)
    save(fig, out, 'fig_order')


def fig_local_convergence(run, out):
    """Residual norm against the local Newton iteration in four trial states."""
    rows = read_rows(os.path.join(run, 'schemes_local_convergence.csv'))
    states = [('subcritical', '(a) NC, undrained 10%\n(subcritical, test B)'),
              ('supercritical', '(b) OCR = 10, undrained 2%\n(supercritical, softening)'),
              ('apex', '(c) NC, isotropic $\\varepsilon_v = -1\\%$\n(apex, $\\rho_{tr} = 0$)'),
              ('near_apex', '(d) as (c), $\\rho_{tr}/a_n = 5\\times10^{-5}$\n(near the apex)')]
    schemes = [('library_exact', dict(color=C1, marker='o', mfc=C1, ms=3.4, ls='-', label='two unknowns (this work)')),
               ('library_full', dict(color=C1, marker='o', mfc='white', ms=3.8, ls='--', label='four unknowns')),
               ('exact_secant', dict(color=C2, marker='s', mfc='white', ms=3.4, ls='--', label='four unknowns, secant $G$')),
               ('trapezoidal', dict(color=ORANGE, marker='v', mfc='white', ms=3.4, ls='-', label='trapezoidal rule, secant $G$'))]
    fig, axs = plt.subplots(1, 4, figsize=(TEXTW, 2.3))
    for ax, (state, title) in zip(axs, states):
        for scheme, sty in schemes:
            sel = [r for r in rows if r['state'] == state and r['scheme'] == scheme]
            if not sel:
                continue
            it = np.array([int(r['iteration']) for r in sel])
            res = np.array([float(r['residual']) for r in sel])
            res = np.where(res > 0, res, np.nan)
            ax.plot(it, res, lw=1.0, mew=0.9, **sty)
        ax.set_yscale('log')
        ax.set_xlabel('Newton iteration')
        ax.set_title(title, fontsize=6.8)
        ax.set_ylim(1e-17, 1e3)
    axs[0].set_ylabel('norm of the residual')
    h, lab = axs[0].get_legend_handles_labels()
    fig.legend(h, lab, loc='upper center', ncol=4, fontsize=6.4, bbox_to_anchor=(0.5, 1.0), handlelength=2.4)
    fig.tight_layout(w_pad=0.6, rect=(0, 0, 1, 0.9))
    save(fig, out, 'fig_local_convergence')


def fig_schemes(run, out):
    """Fig. 6 of v0.7: error against work in the tests of Table 4."""
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
    args = arguments('Figures of the comparison of integration schemes from the CSV files of IntegrationSchemes.')
    fig_error_work(args.rundir, args.outdir)
    fig_order(args.rundir, args.outdir)
    fig_local_convergence(args.rundir, args.outdir)
    fig_schemes(args.rundir, args.outdir)
    fig_secant_sweep(args.rundir, args.outdir)
