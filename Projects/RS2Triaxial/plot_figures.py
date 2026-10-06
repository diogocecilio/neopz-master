"""Figure 4 of the article (drained triaxial tests of the RS2 manual) from the CSV files written by RS2Triaxial.

Usage (run the executable first; it writes its CSV files to the current directory):

    python3 <neopz>/Projects/RS2Triaxial/plot_figures.py [run directory] [-o output directory]

The run directory (default: the current directory) holds the CSV files; the figures are written as PDF and PNG
to <run directory>/figures, or to the output directory.

Figure produced (same panels, curves, markers and annotations as fig_rs2 of figs.py of the Python code):

- fig04_rs2_triaxial: Fig. 4, drained triaxial tests at a material point for the four cases of the RS2 manual
  (NC with constant nu, NC with constant G, OCR = 2 and OCR = 5): q against eps_q (panels a-d) and eps_v against
  eps_a (panels e-h). Solid lines: this work with 400 increments (rs2_<case>_n400.csv); dashed lines: closed form
  of Appendix B.1 up to eps_a = 20% (rs2_<case>_closed.csv); markers: curves "Analytical" and FE of Figs. 8.5-8.8
  of the RS2 manual, digitized (reference/rs2_fig85_88_digitized.json, an unchanged copy of
  dados/rs2_fig85_88_digitized.json of the Python code).
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'Common'))
from mcc_figstyle import C1, C2, C3, INK, INK2, TEXTW, arguments, plt, read_csv, reference, save  # noqa: E402

# (tag of the CSV files, figure of the RS2 manual, panel title, p'0 and p'c0 in kPa); see RS2Triaxial::Cases()
CASES = (('nc_nu', 'Fig8.5', 'NC, constant $\\nu$', 200., 200.),
         ('nc_g', 'Fig8.6', 'NC, constant $G$', 200., 200.),
         ('ocr2', 'Fig8.7', 'OCR = 2', 100., 200.),
         ('ocr5', 'Fig8.8', 'OCR = 5', 100., 500.))


def columns(T):
    """Columns eps_a, p', q, eps_v, eps_q of a CSV table, in the order of the arrays of gen_data.py."""
    return np.column_stack([T['eps_a'], T['p_eff'], T['q'], T['eps_v'], T['eps_q']])


def fig_rs2(run, out):
    """Fig. 4: q-eps_q and eps_v-eps_a of the four drained tests, this work, closed form and RS2."""
    rs = reference(HERE, 'rs2_fig85_88_digitized.json')
    fig, axs = plt.subplots(2, 4, figsize=(TEXTW, 3.75))
    for j, (tag, name, title, p0, pc0) in enumerate(CASES):
        num = columns(read_csv(os.path.join(run, f'rs2_{tag}_n400.csv')))
        ana = columns(read_csv(os.path.join(run, f'rs2_{tag}_closed.csv')))
        ana = ana[ana[:, 0] <= 0.2 + 1e-9]
        if abs(num[0, 1] - p0) > 1e-9 * p0:
            sys.exit(f'rs2_{tag}_n400.csv starts at p\' = {num[0, 1]:g} kPa, expected {p0:g} kPa (case table)')
        d = rs[name]
        ax = axs[0, j]
        ax.plot(num[:, 4], num[:, 2], color=C1, lw=1.6, label='this work (400 steps)', zorder=2)
        ax.plot(ana[:, 4], ana[:, 2], color=INK, lw=0.8, ls=(0, (3, 2)), label='closed form', zorder=3)
        an_pts = np.array(d['q-epsq']['RS2Analytical'])
        ax.plot(an_pts[:, 0], an_pts[:, 1], '^', ms=3.3, mfc='white', mec=C3, mew=0.8, label='RS2 (analytical)')
        fe_pts = np.array(d['q-epsq']['RS2FE'])
        ax.plot(fe_pts[:, 0], fe_pts[:, 1], 'o', ms=3.3, mfc='white', mec=C2, mew=0.8, label='RS2 (FE)')
        ax.set_xlim(0, 0.2)
        ax.set_ylim(0, None)
        ax.set_title(f'({"abcd"[j]}) {title}', fontsize=8)
        ax.text(0.97, 0.06, "$p'_0$ = %g, $p'_{c0}$ = %g kPa" % (p0, pc0), transform=ax.transAxes,
                fontsize=6.3, color=INK2, ha='right', va='bottom')
        if name == 'Fig8.8':
            kq = int(np.argmax(num[:, 2]))
            ax.annotate('peak %.0f kPa' % num[kq, 2], (num[kq, 4], num[kq, 2]), xytext=(0.05, 255), fontsize=6.5,
                        color=INK2, arrowprops=dict(arrowstyle='-', color=INK2, lw=0.6))
        ax.set_xlabel('$\\varepsilon_q$')
        if j == 0:
            ax.set_ylabel('$q$ (kPa)')
        ax = axs[1, j]
        ax.plot(num[:, 0], num[:, 3], color=C1, lw=1.6, zorder=2)
        ax.plot(ana[:, 0], ana[:, 3], color=INK, lw=0.8, ls=(0, (3, 2)), zorder=3)
        an_pts = np.array(d['epsv-epsa']['RS2Analytical'])
        ax.plot(an_pts[:, 0], an_pts[:, 1], '^', ms=3.3, mfc='white', mec=C3, mew=0.8)
        fe_pts = np.array(d['epsv-epsa']['RS2FE'])
        ax.plot(fe_pts[:, 0], fe_pts[:, 1], 'o', ms=3.3, mfc='white', mec=C2, mew=0.8)
        ax.set_xlim(0, 0.2)
        ax.set_title(f'({"efgh"[j]})', fontsize=8)
        if name == 'Fig8.8':
            kmax = int(np.argmax(num[:, 3]))
            ax.annotate('max. compression\n%.2f%%' % (100 * num[kmax, 3]), (num[kmax, 0], num[kmax, 3]),
                        xytext=(0.06, 0.0005), fontsize=6.5, color=INK2,
                        arrowprops=dict(arrowstyle='-', color=INK2, lw=0.6))
            ax.text(0.19, -0.0115, 'dilation', fontsize=6.5, color=INK2, ha='right')
        if j == 0:
            ax.text(0.19, 0.012, 'compaction', fontsize=6.5, color=INK2, ha='right')
        ax.set_xlabel('$\\varepsilon_a$')
        if j == 0:
            ax.set_ylabel('$\\varepsilon_v$')
        print('RS2 %-6s q(20%%) = %.6f kPa, eps_v(20%%) = %.8f, max|q - q_closed| = %.4f kPa'
              % (tag, num[-1, 2], num[-1, 3], np.abs(num[:, 2] - np.interp(num[:, 0], ana[:, 0], ana[:, 2])).max()))
    h, lab = axs[0, 0].get_legend_handles_labels()
    fig.legend(h, lab, loc='lower center', ncol=4, bbox_to_anchor=(0.5, -0.035), fontsize=7.5)
    fig.tight_layout(rect=(0, 0.04, 1, 1), h_pad=0.6, w_pad=0.5)
    save(fig, out, 'fig04_rs2_triaxial')


if __name__ == '__main__':
    args = arguments('Fig. 4 of the article (RS2 drained triaxial tests) from the CSV files of RS2Triaxial.')
    fig_rs2(args.rundir, args.outdir)
