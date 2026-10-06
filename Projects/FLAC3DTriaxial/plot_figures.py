"""Figure 5 of the article (triaxial tests of the FLAC3D verification problem) from the CSV files of FLAC3DTriaxial.

Usage (run the executable first; it writes its CSV files to the current directory):

    python3 <neopz>/Projects/FLAC3DTriaxial/plot_figures.py [run directory] [-o output directory]

The run directory (default: the current directory) holds the CSV files; the figures are written as PDF and PNG
to <run directory>/figures, or to the output directory.

Figure produced (same panels, curves, markers and annotations as fig_itasca of figs.py of the Python code):

- fig05_flac3d_triaxial: Fig. 5, drained (panels a, b) and undrained (panels c, d) triaxial tests with R = 1.6
  and R = 8: q against eps_a and stress paths in the p'-q plane, with the critical state line and the initial
  yield surface of R = 1.6. Solid lines: this work (flac3d_<drained|undrained>_R<1.6|8>.csv, columns eps_a, p_eff,
  q, v, u); dashed lines: closed-form solutions (drained: Appendix B.1 with constant G, the TriaxialDrainedClosedCC
  of camclay_hw.py; undrained: constant volume, the undrained_closed of figs.py), evaluated here; squares: final
  states of FLAC3D (Table 4 of the article, as in figs.py).
"""
import os
import sys

import numpy as np
from matplotlib.lines import Line2D

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'Common'))
from mcc_figstyle import C1, C2, INK, INK2, MUTED, TEXTW, arguments, plt, read_csv, save  # noqa: E402

# Material of Table 1 (FLAC3DTriaxial.h): M, lambda, kappa, shear modulus G (kPa); initial p'0 = 5 kPa
M, LAM, KAP, G = 1.02, 0.2, 0.05, 250.0
# final states (p', q) in kPa of FLAC3D (Itasca verification problem, Table 4 of the article)
FLAC = {('drained', 1.6): (7.573, 7.718), ('drained', 8.0): (7.583, 7.747),
        ('undrained', 1.6): (4.234, 4.312), ('undrained', 8.0): (14.05, 14.42)}


def history(run, test, R):
    """Columns eps_a, p', q, v, u of flac3d_<test>_R<R>.csv, in the order of the arrays of gen_data.py."""
    T = read_csv(os.path.join(run, f'flac3d_{test}_R{R:g}.csv'))
    return np.column_stack([T['eps_a'], T['p_eff'], T['q'], T['v'], T['u']])


def drained_closed(p0, pc0, v0, npts=600):
    """Closed form of the drained test (q = 3(p' - p'0), porous law, constant G): rows {eps_a, p', q, eps_v, eps_q}.

    Same expressions as triaxial_closed of camclay_hw.py (TriaxialDrainedClosedCC, Appendix B.1) with pt = 0,
    beta = 1 and the specific volume v0 frozen in the plastic strains.
    """
    def F(x):
        """Primitive of the plastic distortion against the stress ratio eta = q/p' on the drained path (App. B.1)."""
        return ((1 / M) * np.log(abs((M + x) / (M - x))) - (2 / M) * np.arctan(x / M)
                - np.log(abs(M - x)) / (3 - M) - np.log(M + x) / (3 + M) + 6 * np.log(3 - x) / (9 - M ** 2))

    def el(pp, q):
        """Elastic volumetric and distortional strains (porous law, constant G) at (p', q)."""
        return KAP / v0 * np.log(pp / p0), q / (3 * G)

    A = 9 + M ** 2
    B = -(18 * p0 + M ** 2 * pc0)
    C = 9 * p0 ** 2
    disc = np.sqrt(B * B - 4 * A * C)
    py = min(x for x in ((-B + disc) / (2 * A), (-B - disc) / (2 * A)) if x >= p0 - 1e-9)
    etay = 3 * (py - p0) / py
    rows = [(pp, 3 * (pp - p0)) + el(pp, 3 * (pp - p0)) for pp in np.linspace(p0, py, 41)]
    etas = np.unique(np.r_[etay + (M - etay) * np.arange(1, npts) / npts,
                           M - (M - etay) * np.exp(-12.0 * np.arange(1, npts + 1) / npts)])
    if etay > M:
        etas = etas[::-1]
    for eta in etas:
        pp = 3 * p0 / (3 - eta)
        q = eta * pp
        pc = pp * (1 + eta ** 2 / M ** 2)
        ev, eq = el(pp, q)
        ev += (LAM - KAP) / v0 * np.log(pc / pc0)
        eq += (LAM - KAP) / v0 * (F(eta) - F(etay))
        rows.append((pp, q, ev, eq))
    R = np.array(rows)
    return np.c_[R[:, 3] + R[:, 2] / 3, R[:, 0], R[:, 1], R[:, 2], R[:, 3]]


def undrained_closed(p0, R, v0, npts=600):
    """Closed form of the undrained test (constant volume, constant G): rows {eps_a, p', q, u}.

    Same expressions as undrained_closed of fig_itasca of figs.py: elastic up to the yield at p' = p'0, then
    p' = p'0 ((M^2 + eta^2) / (M^2 R))^(-Lambda) with Lambda = (lambda - kappa)/lambda.
    """
    Lam = (LAM - KAP) / LAM

    def g(x):
        """Primitive of the plastic axial strain against the stress ratio eta = q/p' (undrained path)."""
        return 0.5 * np.log(abs((M + x) / (M - x))) - np.arctan(x / M)

    etay = M * np.sqrt(R - 1.0)
    qy = etay * p0
    rows = [(q / (3 * G), p0, q, q / 3) for q in np.linspace(0, qy, 40)]
    etas = np.unique(np.r_[etay + (M - etay) * np.arange(1, npts) / npts,
                           M - (M - etay) * np.exp(-14.0 * np.arange(1, npts + 1) / npts)])
    if etay > M:
        etas = etas[::-1]
    for eta in etas:
        p = p0 * ((M ** 2 + eta ** 2) / (M ** 2 * R)) ** (-Lam)
        q = eta * p
        rows.append((q / (3 * G) + 2 * Lam * KAP / (M * v0) * (g(eta) - g(etay)), p, q, q / 3 + p0 - p))
    return np.array(rows)


def fig_flac3d(run, out):
    """Fig. 5: q-eps_a and p'-q of the drained and undrained tests with R = 1.6 and R = 8."""
    hist = {(t, R): history(run, t, R) for t in ('drained', 'undrained') for R in (1.6, 8.0)}
    p0 = hist[('drained', 1.6)][0, 1]
    pc16 = 1.6 * p0                                                    # p'c0 of R = 1.6 (initial surface)
    fig, axs = plt.subplots(1, 4, figsize=(TEXTW, 1.95))
    cols = {1.6: C1, 8.0: C2}
    summary = {}
    for R in (8.0, 1.6):                       # R = 1.6 on top (the paths coincide at the beginning)
        z = 3 if R == 1.6 else 2
        h = hist[('drained', R)]
        v0 = h[0, 3]
        an_all = drained_closed(p0, R * p0, v0)                     # whole curve (summary at the final eps_a)
        an = an_all[an_all[:, 0] <= 0.5]
        axs[0].plot(h[:, 0], h[:, 2], color=cols[R], lw=1.5, zorder=z)
        axs[0].plot(an[:, 0], an[:, 2], color=INK, lw=0.7, ls=(0, (3, 2)), zorder=5)
        axs[1].plot(h[:, 1], h[:, 2], color=cols[R], lw=1.5, zorder=z)
        axs[1].plot(an[:, 1], an[:, 2], color=INK, lw=0.7, ls=(0, (3, 2)), zorder=5)
        axs[1].plot(*FLAC[('drained', R)], 's', ms=4.8, mfc='none', mec=cols[R], mew=1.0, zorder=6)
        hu = hist[('undrained', R)]
        v0u = hu[0, 3]
        au_all = undrained_closed(p0, R, v0u)
        au = au_all[au_all[:, 0] <= 0.1]
        axs[2].plot(hu[:, 0], hu[:, 2], color=cols[R], lw=1.5, zorder=z)
        axs[2].plot(au[:, 0], au[:, 2], color=INK, lw=0.7, ls=(0, (3, 2)), zorder=5)
        axs[3].plot(hu[:, 1], hu[:, 2], color=cols[R], lw=1.5, zorder=z)
        axs[3].plot(au[:, 1], au[:, 2], color=INK, lw=0.7, ls=(0, (3, 2)), zorder=5)
        axs[3].plot(*FLAC[('undrained', R)], 's', ms=4.8, mfc='none', mec=cols[R], mew=1.0, zorder=6)
        for t, num, ana in (('drained', h, an_all), ('undrained', hu, au_all)):
            k = int(np.argmax(num[:, 2]))
            ea = num[-1, 0]
            summary[(t, R)] = ("%-9s R = %-3g eps_a = %g: p' = %.5f q = %.5f kPa (closed form %.5f, %.5f; FLAC3D "
                               "%g, %g), peak q = %.4f kPa at eps_a = %g"
                               % (t, R, ea, num[-1, 1], num[-1, 2], np.interp(ea, ana[:, 0], ana[:, 1]),
                                  np.interp(ea, ana[:, 0], ana[:, 2]), *FLAC[(t, R)], num[k, 2], num[k, 0]))
    for key in sorted(summary):
        print(summary[key])
    for ax in (axs[1], axs[3]):
        pm = np.linspace(0, 16, 2)
        ax.plot(pm, M * pm, color=MUTED, lw=0.8, zorder=1)
        ax.text(12.6, M * 12.6 - 2.6, 'CSL', fontsize=7, color=INK2, ha='left')
        pe = np.linspace(0, pc16, 200)                                # initial surface with R = 1.6 (p'c0 = 8)
        ax.plot(pe, M * np.sqrt(np.clip(pe * (pc16 - pe), 0, None)), color=MUTED, lw=0.7, ls=(0, (3, 2)),
                zorder=1)
    axs[1].annotate('initial surface\n($R$ = 1.6)', xy=(2.2, 3.6), xytext=(0.3, 13.6), fontsize=6.3, color=INK2,
                    ha='left', va='top', arrowprops=dict(arrowstyle='-', color=MUTED, lw=0.6))
    axs[1].text(9.2, 4.0, "drained path\n$q = 3(p' - p'_0)$", fontsize=6.3, color=INK2)
    axs[3].annotate("undrained: $p'$ constant\nuntil yield", xy=(5.0, 2.2), xytext=(8.3, 4.6), fontsize=6.3,
                    color=INK2, ha='left', va='top', arrowprops=dict(arrowstyle='-', color=MUTED, lw=0.6))
    for ax, key in ((axs[0], 'drained'), (axs[2], 'undrained')):
        h8 = hist[(key, 8.0)]
        k8 = int(np.argmax(h8[:, 2]))
        ax.annotate('peak %.1f kPa' % h8[k8, 2], (h8[k8, 0], h8[k8, 2]),
                    xytext=(0.30, 0.93), textcoords='axes fraction', fontsize=6.5, color=INK2,
                    arrowprops=dict(arrowstyle='-', color=INK2, lw=0.6))
        h16 = hist[(key, 1.6)]
        ax.text(h16[-1, 0] * 0.97, h16[-1, 2] - 2.4, '$R$ = 1.6', fontsize=6.8, color=C1, ha='right')
        ax.text(h8[-1, 0] * 0.97, h8[-1, 2] + 1.0, '$R$ = 8', fontsize=6.8, color=C2, ha='right')
    axs[0].set(xlabel='$\\varepsilon_a$', ylabel='$q$ (kPa)', xlim=(0, 0.5), ylim=(0, 19))
    axs[1].set(xlabel="$p'$ (kPa)", ylabel='$q$ (kPa)', xlim=(0, 16), ylim=(0, 19))
    axs[2].set(xlabel='$\\varepsilon_a$', ylabel='$q$ (kPa)', xlim=(0, 0.1), ylim=(0, 19))
    axs[3].set(xlabel="$p'$ (kPa)", ylabel='$q$ (kPa)', xlim=(0, 16), ylim=(0, 19))
    for k, (ax, t) in enumerate(zip(axs, ('drained', 'drained', 'undrained', 'undrained'))):
        ax.set_title(f'({"abcd"[k]}) {t}', fontsize=8)
    handles = [Line2D([], [], color=C1, lw=1.3, label='$R = 1.6$'), Line2D([], [], color=C2, lw=1.3, label='$R = 8$'),
               Line2D([], [], color=INK, lw=0.8, ls=(0, (4, 2)), label='closed form'),
               Line2D([], [], ls='', marker='s', ms=4.5, mfc='white', mec=INK2, label='FLAC3D (final state)')]
    fig.legend(handles=handles, loc='lower center', ncol=4, bbox_to_anchor=(0.5, -0.07), fontsize=7.5)
    fig.tight_layout(rect=(0, 0.05, 1, 1), w_pad=0.6)
    save(fig, out, 'fig05_flac3d_triaxial')


if __name__ == '__main__':
    args = arguments('Fig. 5 of the article (FLAC3D triaxial tests) from the CSV files of FLAC3DTriaxial.')
    fig_flac3d(args.rundir, args.outdir)
