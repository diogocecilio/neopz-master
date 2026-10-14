"""Figure of Sect. 7.1 of the article (v1.0): triaxial tests of example 1 of Krabbenhoft and Lyamin (2012) with one
Hex20-Hex8 element against the closed forms of Appendix B, from the CSV files written by KLTriaxial.

Usage (run the executable first; it writes its CSV files to the current directory):

    python3 <neopz>/Projects/KLTriaxial/plot_figures.py [run directory] [-o output directory]

Figures produced (PDF and PNG in <run directory>/figures or in the output directory):

- fig_kl_triaxial (layout of v1.1): (a) q and (b) eps_v against eps_a of the drained tests (NC and OCR = 10),
  (c) effective stress paths of the four tests normalized by p'_c0, with the initial yield surface (the same in
  normalized coordinates) and the critical state line, (d) and (e) q and the excess pore pressure against eps_a of
  the undrained tests, NC and OCR = 10, each with its own scale, (f) the legend. Lines: closed forms; markers: one
  Hex20-Hex8 element with 100 increments (every fifth increment drawn). Files kl_<case>_closed.csv and
  kl_<case>_fe_n100.csv. The layout of v1.0 is kept in plot_figures_v1.0.py.
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
    """v1.1 layout: (a) q and (b) eps_v of the drained tests; (c) effective stress paths normalized by p'_c0, where
    the initial yield surfaces of the NC and OCR = 10 clays coincide; (d) and (e) q and excess pore pressure of the
    undrained tests, NC and OCR = 10, each with its own scale; (f) the legend. Lines: closed forms; markers: one
    Hex20-Hex8 element with 100 increments (every fifth increment)."""
    fig = plt.figure(figsize=(TEXTW, 4.35))
    gs = fig.add_gridspec(2, 3, wspace=0.42, hspace=0.52, left=0.075, right=0.985, top=0.94, bottom=0.10)
    ax_q, ax_v, ax_p = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[0, 2])
    ax_unc, ax_uoc, ax_leg = fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1]), fig.add_subplot(gs[1, 2])
    every = 5
    ms, mew = 3.0, 0.8
    for tag, (label, color, marker) in CASES.items():
        cl = read(os.path.join(run, f'kl_{tag}_closed.csv'))
        fe = read(os.path.join(run, f'kl_{tag}_fe_n100.csv'))
        pc0 = (10 if tag.endswith('ocr10') else 1) * P0
        sel = slice(0, None, every)
        if tag.startswith('drained'):
            ax_q.plot(100 * cl['eps_a'], cl['q'], '-', color=color, lw=1.0)
            ax_q.plot(100 * fe['eps_a'][sel], fe['q'][sel], marker, color=color, mfc='white', ms=ms, mew=mew)
            ax_v.plot(100 * cl['eps_a'], 100 * cl['eps_v'], '-', color=color, lw=1.0)
            ax_v.plot(100 * fe['eps_a'][sel], 100 * fe['eps_v'][sel], marker, color=color, mfc='white', ms=ms, mew=mew)
        else:
            ax = ax_unc if tag == 'undrained_nc' else ax_uoc
            ax.plot(100 * cl['eps_a'], cl['q'], '-', color=color, lw=1.0)
            ax.plot(100 * fe['eps_a'][sel], fe['q'][sel], marker, color=color, mfc='white', ms=ms, mew=mew)
            ax.plot(100 * cl['eps_a'], cl['u'], '--', color=color, lw=1.0)
            ax.plot(100 * fe['eps_a'][sel], fe['u'][sel], marker, color=color, mfc=color, ms=ms - 0.4, mew=mew)
            ax.axhline(0.0, color=MUTED, lw=0.5, zorder=0)
        ax_p.plot(cl['p_eff'] / pc0, cl['q'] / pc0, '-', color=color, lw=1.0)
        sp = slice(0, None, 2 * every)  # fewer markers where the paths of the OCR = 10 tests overlap
        ax_p.plot(fe['p_eff'][sp] / pc0, fe['q'][sp] / pc0, marker, color=color, mfc='white', ms=ms - 0.4, mew=mew)
    # initial yield surface (the same ellipse for both clays in normalized coordinates) and critical state line
    x = np.linspace(0.0, 1.0, 400)
    ax_p.plot(x, M * np.sqrt(np.maximum(x * (1.0 - x), 0.0)), ':', color=MUTED, lw=0.9)
    xc = np.array([0.0, 1.4])
    ax_p.plot(xc, M * xc, '-', color=INK2, lw=0.7)
    ax_p.annotate('CSL', (1.30, M * 1.30), xytext=(-14, 2), textcoords='offset points', fontsize=6.4, color=INK2,
                  ha='right')
    ax_p.set_xlim(0.0, 1.4)
    ax_p.set_ylim(0.0, 0.95)
    eps = '$\\varepsilon_a$ (%)'
    ax_q.set_xlabel('axial strain ' + eps)
    ax_q.set_ylabel('$q$ (kPa)')
    ax_q.set_title('(a) drained: deviatoric stress', fontsize=7.6)
    ax_v.set_xlabel('axial strain ' + eps)
    ax_v.set_ylabel('$\\varepsilon_v$ (%), compression positive')
    ax_v.set_title('(b) drained: volumetric strain', fontsize=7.6)
    ax_p.set_xlabel("$p\\,'/p\\,'_{c0}$")
    ax_p.set_ylabel("$q/p\\,'_{c0}$")
    ax_p.set_title("(c) effective stress paths, normalized", fontsize=7.6)
    for ax, title in ((ax_unc, '(d) undrained, NC'), (ax_uoc, '(e) undrained, OCR = 10')):
        ax.set_xlabel('axial strain ' + eps)
        ax.set_ylabel('$q$, $\\Delta p_w$ (kPa)')
        ax.set_title(title, fontsize=7.6)
    # legend panel: colours = tests, line/marker = closed form/element, solid/dashed = q/excess pore pressure
    ax_leg.axis('off')
    from matplotlib.lines import Line2D
    handles = [Line2D([], [], color=c, marker=mk, mfc='white', ms=ms, mew=mew, lw=1.0, label=lab)
               for lab, c, mk in CASES.values()]
    handles += [Line2D([], [], color=INK2, lw=1.0, label='closed form (Appendix B)'),
                Line2D([], [], color=INK2, lw=0, marker='o', mfc='white', ms=ms, mew=mew,
                       label='NeoPZ, one element, 100 increments'),
                Line2D([], [], color=INK2, lw=1.0, ls='--', label='excess pore pressure $\\Delta p_w$ (dashed)'),
                Line2D([], [], color=INK2, lw=0, marker='o', mfc=INK2, ms=ms - 0.4, mew=mew,
                       label='$\\Delta p_w$ of the element (filled)'),
                Line2D([], [], color=MUTED, lw=0.9, ls=':', label='initial yield surface')]
    ax_leg.legend(handles=handles, loc='center left', fontsize=6.4, frameon=False, handlelength=2.4,
                  borderaxespad=0.0)
    save(fig, out, 'fig_kl_triaxial')


if __name__ == '__main__':
    args = arguments('Figure of the Krabbenhoft-Lyamin triaxial tests with one element (Sect. 7.1 of the article).')
    fig_kl(args.rundir, args.outdir)
