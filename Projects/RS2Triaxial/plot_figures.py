"""Figure 4 of the article (drained triaxial tests of the RS2 manual) from the CSV files written by RS2Triaxial.

Usage (run the executable first; it writes its CSV files to the current directory):

    python3 <neopz>/Projects/RS2Triaxial/plot_figures.py [run directory] [-o output directory]

The run directory (default: the current directory) holds the CSV files; the figures are written as PDF and PNG
to <run directory>/figures, or to the output directory.

Figures produced:

- fig04_rs2_triaxial (same panels, curves, markers and annotations as fig_rs2 of figs.py of the Python code): Fig. 4, drained triaxial tests at a material point for the four cases of the RS2 manual
  (NC with constant nu, NC with constant G, OCR = 2 and OCR = 5): q against eps_q (panels a-d) and eps_v against
  eps_a (panels e-h). Solid lines: this work with 400 increments (rs2_<case>_n400.csv); dashed lines: closed form
  of Appendix B.1 up to eps_a = 20% (rs2_<case>_closed.csv); markers: curves "Analytical" and FE of Figs. 8.5-8.8
  of the RS2 manual, digitized (reference/rs2_fig85_88_digitized.json, an unchanged copy of
  dados/rs2_fig85_88_digitized.json of the Python code).
- supplementary_rs2_element_check: finite element check of the material point solution: (a) the single Hex20-Hex8
  element (unit cube, 2 x 2 x 2 Gauss points) with the boundary conditions of the drained RS2 tests, drawn from
  rs2_mesh_*.csv with the module Common/mcc_hexmodel.py (the same element is the model of the FLAC3D tests,
  figure fig_single_element_model of FLAC3DTriaxial); (b) q against eps_a of the four cases, element (markers) and
  material point (lines), 400 increments (rs2_<case>_fe.csv); (c) |q_FE - q_point| along the path.

The script also compares this work with the last points of the digitized RS2 curves (numbers of the text of
Sect. 6.1: the end values of q differ from the RS2 analytical curves by less than 0.3 %; the RS2 finite element
curves of the normally consolidated cases are softer in q and have larger volumetric strains) and writes the
comparison to rs2_digitized_comparison.csv in the output directory: for each case and RS2 curve, the last digitized
point of q-eps_q and of eps_v-eps_a, this work interpolated at the same abscissa and the relative difference
(RS2 - this work)/this work.
"""
import csv
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'Common'))
from mcc_figstyle import C1, C2, C3, C4, INK, INK2, TEXTW, arguments, plt, read_csv, reference, save  # noqa: E402
import mcc_hexmodel as hexmodel  # noqa: E402  (Common/mcc_hexmodel.py: drawing of the hexahedral models)

# boundary ids of RS2Triaxial.h (the same as FLAC3DTriaxial.h)
EX0, EX1, EY0, EY1, EZ0, EZ1, EPDRAINED = -21, -22, -23, -24, -25, -26, -30

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


def compare_digitized(run, out):
    """Relative difference between the last points of the digitized RS2 curves and this work (400 increments)."""
    rs = reference(HERE, 'rs2_fig85_88_digitized.json')
    rows = []
    for tag, name, title, p0, pc0 in CASES:
        num = read_csv(os.path.join(run, f'rs2_{tag}_n400.csv'))
        for curve in ('RS2Analytical', 'RS2FE'):
            a = np.array(rs[name]['q-epsq'][curve])[-1]
            b = np.array(rs[name]['epsv-epsa'][curve])[-1]
            q = np.interp(a[0], num['eps_q'], num['q'])
            ev = np.interp(b[0], num['eps_a'], num['eps_v'])
            rows.append([tag, curve, a[0], a[1], q, 100 * (a[1] - q) / q, b[0], b[1], ev, 100 * (b[1] - ev) / ev])
            print('RS2 %-6s %-13s last point: q = %.2f kPa at eps_q = %.4f (this work %.2f, %+.2f%%), eps_v = %.5f at '
                  'eps_a = %.4f (this work %.5f, %+.2f%%)' % (tag, curve, a[1], a[0], q, rows[-1][5], b[1], b[0], ev,
                                                              rows[-1][9]))
    with open(os.path.join(out, 'rs2_digitized_comparison.csv'), 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['case', 'curve', 'eps_q', 'q_rs2', 'q_this_work', 'q_diff_percent', 'eps_a', 'eps_v_rs2',
                    'eps_v_this_work', 'eps_v_diff_percent'])
        for r in rows:
            w.writerow(r[:2] + ['%.10g' % v for v in r[2:]])


def fig_element_check(run, out):
    """Finite element check: the Hex20-Hex8 element and its boundary conditions, q of the element against the
    material point and the difference along the path."""
    mesh = hexmodel.Mesh(run, 'rs2_mesh')
    if len(mesh.faces_with(EPDRAINED)) != 6:
        sys.exit('rs2_mesh_faces.csv: the six faces must be drained')
    check = read_csv(os.path.join(run, 'rs2_fe_check.csv'))
    fig = plt.figure(figsize=(TEXTW, 2.25))
    ax = fig.add_axes([0.0, 0.02, 0.36, 0.92])
    view = hexmodel.View(24, 33)
    P = view.point

    def color(ids):
        if EZ1 in ids:
            return hexmodel.FACE_TOP
        if EX1 in ids or EY1 in ids:
            return hexmodel.FACE_LATERAL
        return hexmodel.FACE_RESTRAINED
    box = hexmodel.draw_model(ax, mesh, view, color)
    for x, y in ((0.15, 0.15), (0.85, 0.15), (0.5, 0.5), (0.15, 0.85), (0.85, 0.85)):
        hexmodel.arrow(ax, P(x, y, 1.38), P(x, y, 1.03), color=INK)
    for z in (0.2, 0.5, 0.8):
        hexmodel.arrow(ax, P(1.5, 0.0, z), P(1.02, 0.0, z), color=C2)
        hexmodel.arrow(ax, P(0.0, 1.5, z), P(0.0, 1.02, z), color=C2)
    xt = P(0.0, 1.5, 0.5)[0] + 0.1
    ax.text(xt, P(0.5, 0.5, 1.3)[1], '$u_z$ prescribed\nto $\\varepsilon_a$ = 20%', fontsize=6.4, va='center')
    ax.text(xt, P(0, 1, 0.55)[1], "cell pressure\n$\\sigma_c = p'_0$", fontsize=6.4, va='center', color=INK)
    ax.text(xt, P(0, 1, 0.0)[1] - 0.12, '$p_w$ = 0 at the\nvertices (drained)', fontsize=6.4, va='center')
    ax.text(box[0] - 0.4, box[2] - 0.32, 'hidden: $u_x$ = 0 on $x$ = 0, $u_y$ = 0 on $y$ = 0,\n$u_z$ = 0 on $z$ = 0; '
            '2 × 2 × 2 Gauss points', fontsize=6.2, color=INK2, va='center')
    ax.text(box[0] - 0.4, P(0.5, 0.5, 1.38)[1] + 0.12, '(a) one Hex20–Hex8 element', fontsize=8, va='bottom')
    ax.set_aspect('equal')
    ax.set_xlim(box[0] - 0.45, xt + 1.05)
    ax.set_ylim(box[2] - 0.5, P(0.5, 0.5, 1.38)[1] + 0.3)
    ax.axis('off')
    axb = fig.add_axes([0.435, 0.2, 0.235, 0.66])
    axc = fig.add_axes([0.785, 0.2, 0.205, 0.66])
    for j, ((tag, name, title, p0, pc0), col) in enumerate(zip(CASES, (C1, C2, C3, C4))):
        fe = read_csv(os.path.join(run, f'rs2_{tag}_fe.csv'))
        lab = title.replace('constant ', 'const. ')
        axb.plot(fe['eps_a'], fe['q_point'], color=col, lw=1.1, label=lab)
        axb.plot(fe['eps_a'][::16], fe['q'][::16], 'o', ms=2.6, mfc='white', mec=col, mew=0.7)
        d = np.abs(fe['q_minus_q_point'])
        axc.semilogy(fe['eps_a'][1:], np.maximum(d[1:], 1e-12), color=col, lw=0.9)
        print('RS2 element check %-6s q(20%%) = %.9f (element) | %.9f (material point) kPa, max |q_FE - q_point| = '
              '%.2e kPa, %.4f evaluations per increment'
              % (tag, check['q_fe'][j], check['q_point'][j], check['max_diff_q'][j], check['mean_evaluations'][j]))
    axb.set(xlabel='$\\varepsilon_a$', ylabel='$q$ (kPa)', xlim=(0, 0.2), ylim=(0, None))
    axb.set_title('(b) element and material point', fontsize=8)
    axb.plot([], [], 'o', ms=2.6, mfc='white', mec=INK2, mew=0.7, label='element')
    axb.plot([], [], color=INK2, lw=1.1, label='material point')
    axb.legend(loc='lower right', fontsize=5.9, handlelength=1.3, labelspacing=0.2, ncol=2, columnspacing=0.8,
               borderaxespad=0.3)
    axc.set(xlabel='$\\varepsilon_a$', ylabel='$|q_{FE} - q_{point}|$ (kPa)', xlim=(0, 0.2), ylim=(1e-10, 1e-4))
    axc.set_title('(c) difference', fontsize=8)
    save(fig, out, 'supplementary_rs2_element_check')


if __name__ == '__main__':
    args = arguments('Figures of the RS2 drained triaxial tests (Fig. 4 and the finite element check) from the CSV '
                     'files of RS2Triaxial.')
    fig_rs2(args.rundir, args.outdir)
    compare_digitized(args.rundir, args.outdir)
    fig_element_check(args.rundir, args.outdir)
