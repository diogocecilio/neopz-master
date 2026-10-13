"""Figure of Sect. 7.1 of the article (v1.0): triaxial tests of example 1 of Krabbenhoft and Lyamin (2012) with one
Hex20-Hex8 element against the closed forms of Appendix B, from the CSV files written by KLTriaxial.

Usage (run the executable first; it writes its CSV files to the current directory):

    python3 <neopz>/Projects/KLTriaxial/plot_figures.py [run directory] [-o output directory]

Figures produced (PDF and PNG in <run directory>/figures or in the output directory):

- fig_kl_triaxial: (a) q against eps_a of the drained tests (NC and OCR = 10), (b) eps_v against eps_a of the
  drained tests, (c) q and the excess pore pressure u against eps_a of the undrained tests, (d) effective stress
  paths of the four tests with the initial yield surfaces and the critical state line. Lines: closed forms; markers:
  one Hex20-Hex8 element with 100 increments (every fifth increment drawn). Files kl_<case>_closed.csv and
  kl_<case>_fe_n100.csv.
"""
import csv
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'Common'))
from mcc_figstyle import C1, C2, C3, C4, INK2, MUTED, TEXTW, arguments, plt, save  # noqa: E402

M = 3 * np.sin(np.radians(24)) / np.sqrt(3 + np.sin(np.radians(24)) ** 2)
P0 = 200.0
CASES = {'drained_nc': ('drained, NC', C1, 'o'), 'drained_ocr10': ('drained, OCR = 10', C2, 's'),
         'undrained_nc': ('undrained, NC', C3, '^'), 'undrained_ocr10': ('undrained, OCR = 10', C4, 'D')}


def read(path):
    if not os.path.exists(path):
        sys.exit(f'file not found: {path} (run the executable in this directory first)')
    with open(path, newline='') as f:
        rows = list(csv.DictReader(f))
    return {k: np.array([float(r[k]) for r in rows]) for k in rows[0]}


def fig_kl(run, out):
    fig, axs = plt.subplots(2, 2, figsize=(TEXTW, 5.0))
    (ax_q, ax_v), (ax_u, ax_p) = axs
    every = 5
    for tag, (label, color, marker) in CASES.items():
        cl = read(os.path.join(run, f'kl_{tag}_closed.csv'))
        fe = read(os.path.join(run, f'kl_{tag}_fe_n100.csv'))
        drained = tag.startswith('drained')
        sel = slice(0, None, every)
        if drained:
            ax_q.plot(100 * cl['eps_a'], cl['q'], '-', color=color, lw=1.0, label=f'{label}, closed form')
            ax_q.plot(100 * fe['eps_a'][sel], fe['q'][sel], marker, color=color, mfc='white', ms=3.2, mew=0.8,
                      label=f'{label}, one element')
            ax_v.plot(100 * cl['eps_a'], 100 * cl['eps_v'], '-', color=color, lw=1.0)
            ax_v.plot(100 * fe['eps_a'][sel], 100 * fe['eps_v'][sel], marker, color=color, mfc='white', ms=3.2, mew=0.8)
        else:
            ax_u.plot(100 * cl['eps_a'], cl['q'], '-', color=color, lw=1.0, label=f'{label}: $q$')
            ax_u.plot(100 * fe['eps_a'][sel], fe['q'][sel], marker, color=color, mfc='white', ms=3.2, mew=0.8)
            ax_u.plot(100 * cl['eps_a'], cl['u'], '--', color=color, lw=1.0, label=f'{label}: $u$')
            ax_u.plot(100 * fe['eps_a'][sel], fe['u'][sel], marker, color=color, mfc=color, ms=2.6, mew=0.8)
        ax_p.plot(cl['p_eff'], cl['q'], '-', color=color, lw=1.0, label=label)
        ax_p.plot(fe['p_eff'][sel], fe['q'][sel], marker, color=color, mfc='white', ms=3.2, mew=0.8)
    # initial yield surfaces and the critical state line
    for pc0, ls in ((P0, ':'), (10 * P0, ':')):
        p = np.linspace(0, pc0, 300)
        q = M * np.sqrt(np.maximum(p * (pc0 - p), 0.0))
        ax_p.plot(p, q, ls, color=MUTED, lw=0.8)
    p = np.linspace(0, 2 * P0 * 10, 2)
    ax_p.plot(p, M * p, '-', color=INK2, lw=0.7, label='CSL')
    ax_p.set_xlim(0, 10 * P0)
    ax_p.set_ylim(0, None)
    ax_q.set_xlabel('axial strain $\\varepsilon_a$ (%)')
    ax_q.set_ylabel('$q$ (kPa)')
    ax_q.set_title('(a) drained: deviatoric stress', fontsize=7.6)
    ax_v.set_xlabel('axial strain $\\varepsilon_a$ (%)')
    ax_v.set_ylabel('volumetric strain $\\varepsilon_v$ (%)')
    ax_v.set_title('(b) drained: volumetric strain (compression positive)', fontsize=7.6)
    ax_u.set_xlabel('axial strain $\\varepsilon_a$ (%)')
    ax_u.set_ylabel('$q$, $u$ (kPa)')
    ax_u.set_title('(c) undrained: deviatoric stress and excess pore pressure', fontsize=7.6)
    ax_p.set_xlabel("$p'$ (kPa)")
    ax_p.set_ylabel('$q$ (kPa)')
    ax_p.set_title('(d) effective stress paths, initial surfaces (dotted)', fontsize=7.6)
    ax_q.legend(fontsize=6.0, loc='lower right')
    ax_u.legend(fontsize=6.0, loc='center right', ncol=1)
    ax_p.legend(fontsize=6.0, loc='upper left')
    fig.tight_layout(w_pad=1.0, h_pad=1.0)
    save(fig, out, 'fig_kl_triaxial')


if __name__ == '__main__':
    args = arguments('Figure of the Krabbenhoft-Lyamin triaxial tests with one element (Sect. 7.1 of the article).')
    fig_kl(args.rundir, args.outdir)
