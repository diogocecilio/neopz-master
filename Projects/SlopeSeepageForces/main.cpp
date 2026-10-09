// Slope stability under the seepage forces of a rapid drawdown (Ceron, Cecilio, Linn & Maghous, IJNAMG 2025,
// doi:10.1002/nag.3993, Section 4.3), built on Projects/SlopeDrawdown and Projects/SlopeMohrCoulomb:
//  1. parametric slope meshes (SlopeGeometry.h, DelaunayMesher.h): large hydraulic domain and stability domain;
//  2. steady anisotropic seepage of the excess pore pressure u after the drawdown (SeepageFE.h,
//     AnisotropicDarcy.h): J(u_FE) = 1/2 int grad u . K grad u;
//  3. seepage force field f = -grad u_FE at any point (SeepageForceField.h);
//  4. semi-analytical optimal field K^-1 v'_opt of the paper (AnalyticalSeepage.h, port of
//     scripts/analytical_seepage.py: Eqs. 29-40 with the corrected Eq. 31 and the degenerate m -> 0 optimum of steep
//     slopes), a slope::ForceField like the FE field;
//  5. kinematic limit analysis (LimitAnalysis.h, port of scripts/limit_analysis.py): upper bound Gamma = min P_mr /
//     (P_gamma + P_u) over rotational log-spiral mechanisms with the seepage forces of any slope::ForceField (P_u by the
//     boundary formula for the FE field f = -grad u_h); H_crit = Gamma H;
//  6. FEM stability: gravity increase of SlopeAnalysis.h with b = lambda (gamma' g + f) (FEMStability.h);
//     Gamma_FEM = lambda_crit, H_crit = lambda_crit H.
// Commands: SeepageCommands.cpp (mesh, seepage, probe, fig5, analytical), LimitAnalysisCommands.cpp (la, labatch,
// fig8, fig9), FEMCommands.cpp (fs, fembatch), SelfTests.cpp (check, verify); shared: Commands.h (options, problem,
// helpers), FigureCases.h (limit-analysis cases of the figures), ResumableCSV.h (production output).
//
// Usage: SlopeSeepageForces <command> [key=value ...]
//   mesh     meshes, quality and consistency (vtk=1 writes them; sweep=1 checks beta = 15..90, h_w / H = 0..1)
//   verify   manufactured solutions of the anisotropic seepage solver (exact for P2) and point location (alpha=5)
//   seepage  drawdown seepage: J, J / (k_h H^2 gamma_w^2), min p, for href=<list> uniform refinements (with 3 or more:
//            observed order and Richardson extrapolation)
//   probe    u and f = -grad u at the points of pts=<file> (lines x_paper,y_paper), paper coordinates
//   fig5     J / (k_h H^2 gamma_w^2) for alphas=<list> betas=<list> against the dashed curves of the paper's Fig. 5,
//            and -J*(v'_opt) / (k_h H^2 gamma_w^2) of the analytical field against the solid curves; out=<csv>: the
//            production grid into a resumable csv (see fig8 / fig9)
//   analytical  analytical field K^-1 v'_opt (AnalyticalSeepage.h): m, C, D, F, J* (Eq. 40) and the force at sample
//            points; ref=<file> compares with the Python reference of scripts/analytical_seepage_reference.py
//   fs       gravity-increase factor of the stability domain with the seepage forces (FE or analytical field, none or
//            dry), refinement cycles of the plastic zone and the extrapolations of lambda to h -> 0
//   la       kinematic limit analysis: Gamma = min P_mr / (P_gamma + P_u) over log-spiral mechanisms I (B on the face)
//            and II (B on the toe ground), H_crit = Gamma H; prints Gamma, the mechanism (theta1, theta2, eta or d/H;
//            A, B, C), P_mr, P_gamma, P_u
//   labatch  la for every line of cases=<file> (key=value options per line), rows appended to out=<csv>; cases
//            already in out are skipped (resumable)
//   fig8     paper Fig. 8: H_crit = Gamma(H) H versus h_w / H (alpha = 1, H = 1 m), London (beta = 30, 60) and Israeli
//            (35, 60) clay panels, curves vopt (K^-1 v'_opt) and FE (-grad u'_FE, box 50 / 10 / 30 m), limit analysis
//   fig9     paper Fig. 9: Gamma versus beta (H = 5 m, h_w = H) for alphas=1,5,10, curves vopt and FE (box in metres)
//            fig5 out=, fig8 and fig9 append one row per case to a resumable csv (default results/cpp/fig*.csv;
//            cases already there with the same settings are skipped) and log the progress to the .log next to it
//   fembatch FEM gravity-increase factor (FEMStability.h) of cases of Fig. 9 (fig=9: alphas=1,5,10 x betas=30,45,60,
//            75,90 x curves FE, vopt) or Fig. 8 (fig=8: London30/60, Israeli35/60 x hws=0,0.2,0.5,1, h_w = 0 once),
//            an independent check of the limit-analysis curves: each row holds lambda of every refinement cycle, the
//            extrapolations to h -> 0 (Gamma_FEM = order-1 extrapolation 2 lambda_last - lambda_previous) and the
//            limit analysis of the same case; resumable csv results/cpp/fem_fig9.csv / fem_fig8.csv
//            (scripts/run_fem_batch.sh)
//   check    self-tests (about 20-30 s, exit code 1 on failure): manufactured solutions and orientation of K; Eq. 21
//            and the far-side data on the solved field, boundary ids, p >= 0, f = -grad u, point location on
//            vertices / edges / just outside, continuity of the P2 field, evaluator from 4 threads; load vector of
//            the stability problem: resultant = -int u n ds - gamma' A e_y (divergence theorem), lambda scales
//            gamma' and f, forms u and p identical, threaded assembly = serial; limit analysis: closed forms vs
//            quadrature and polygon, domain vs boundary P_u, limit_analysis.py values at fixed mechanisms, dry
//            stability numbers, Fig. 8 h_w = 0 ends, scale invariance, uniform field, threads, FE field vs Python;
//            analytical field: Python reference (data/analytical_seepage_reference.csv), Fig. 5 solid curves,
//            degenerate optimum vs A > P, J* by quadrature vs Eq. 40, div v = 0, f = 0 outside, threads; figure
//            drivers: resumable csv (keys, cut rows), hydraulic box fixed in metres, Fig. 8 soil sets, tiny fig8 / fig9
//            runs repeated (the second must skip every case); fembatch: extrapolations, stability box inside the
//            hydraulic box, similarity (load vectors at H and 4 H, lambda H of tiny runs), tiny runs repeated
// Options (default; numeric options must be numbers, unknown options are an error):
//   slope/soil: H=5 beta=45 hw=1 (h_w / H) gamma=20 gammaw=9.81 c=10 phi=30 E=20000 nu=0.3
//   seepage:    alpha=1 (k_h / k_v) kv=1 horder=2 hbc=zero_lb|impermeable|zero_b|zero_l|zero_lbr|toe_r (far sides,
//               presets of scripts/fe_seepage.py; zero_lb: u = 0 left and base, see SeepageFE.h) hbcleft= hbcbottom=
//               hbcright= (noflow|zero|toe, override of the preset)
//   hydraulic mesh: hmesh=gen|trig hleft=50 hright=10 hdepth=30 (units of H, from O, T, T) hh0=0.025 hhs=0.0625
//               hgrade=0.15 hhmax=2 (sizes in units of H) href=0 (list for seepage, e.g. href=0,1,2)
//   stability mesh: smesh=gen|trig sa=2 (extents sa H + H / tan(beta)) sleft= sright= sdepth= (override, units of
//               H) sh0=0.25 shs=0.25 sgrade=0.25 shmax=1 sref=0
//               mesh=trig: TriGMesh(1 + ref) of SlopeMohrCoulomb (H = 10 m, beta = 45 deg, 70 x 40 m), the mesh of
//               SlopeDrawdown
//   output:     seepage: csv=<file> (u and f on a grid around the slope, paper coordinates), vtk=<file> (u, -grad u
//               and the Darcy velocity on the hydraulic mesh, NeoPZ coordinates), both for the last href
//   fig5:       alphas=1,2,4,10 betas=15,30,45,60,75,90 paper=<csv> (default data/fig5_vector_fill_polygons.csv)
//               fe=1 (0: analytical curve only) lm=10 (L_m / H of the analytical field); out=<csv>: production grid
//               alphas=1,2,4,10 betas=15,22.5,..,90 href=1
//   analytical: lm=10 (L_m / H, R_e) pts=<file> (lines x_paper,y_paper: zone and f, paper coordinates) csv=<file>
//               (grid) ref=<file> (Python reference, e.g. data/analytical_seepage_reference.csv) bench=<n> (timing)
//   fs:         water=seepage (= fe) | analytical (K^-1 v'_opt, lm=10) | none (h_w = 0: f = 0 with gamma' = gamma -
//               gamma_w) | dry; form=u|p|p+ (water=seepage: b = lambda (gamma' g - grad u) | lambda (gamma_sat g -
//               grad p) | lambda (gamma_sat g - grad max(p, 0)) as SlopeDrawdown) nref=3 srm=0 vtk=<prefix>
//               mark=0.1 (refinement of the elements with sqrt(J2(eps_p)) >= mark * max at collapse)
//               maxnewton=100 tolfs=0.002 (driver: Newton iterations, relative step of the continuation;
//               maxnewton=30 is the setting of SlopeMohrCoulomb / SlopeDrawdown, see FEMStability.h)
//               checkforms=1 (compares the load vectors of the forms u and p) hboxm=<left>,<right>,<depth>
//               (hydraulic box in metres, as la; overrides hleft, hright, hdepth)
//   la:         water=fe|analytical|none|dry (default fe, none if hw=0: FE field -grad u_h | K^-1 v'_opt (lm=10) |
//               f = 0 with gamma' = gamma - gamma_w (h_w = 0) | f = 0 and gamma_w = 0) hboxm=<left>,<right>,<depth>
//               (hydraulic box in metres, e.g. 50,10,30; overrides hleft..) href=0 pu=auto|domain|boundary (P_u: auto
//               = boundary formula for the FE field) mech=I,II seeds=0,1,2 np=40 niter=150 pools=25 (PSO; pools:
//               further initial pools of 4 np admissible mechanisms when fewer than np / 2 have P_ext > 0; 0 =
//               limit_analysis.py) qsearch=coarse qfinal=fine (coarse|medium|fine|xfine|ref|dense) polish=1 dmax=10
//               (d / H of mechanism II) threads=4 (at most the CPUs; runs class x seed in parallel) verbose=0
//               out=<csv> (append a row; A, B, C in paper coordinates) x=theta1,theta2,s (no optimisation: rates of
//               work of this mechanism by every rule)
//   labatch:    cases=<file> out=<csv> threads=4
//   fig8, fig9: limit analysis as la (seeds=0,1,2 np=40 niter=150 pools=25 qsearch=coarse qfinal=fine polish=1
//               dmax=10 mech=I,II lm=10, but threads=2) with hboxm=50,10,30 (hydraulic box in metres; hleft, hright,
//               hdepth are not accepted) href=1 hbc=zero_lb; curves=vopt,FE out=<csv>
//   fig8:       soil=swapped|table1 (Table 1 (c, phi) pairs exchanged between the panels | as printed) gammaw=9.8 H=1
//               panels=London30,London60,Israeli35,Israeli60 hws=0,0.05,0.1,0.2,..,1 (default results/cpp/fig8.csv,
//               fig8_table1.csv for soil=table1; beta, hw, c, phi, gamma come from the panels)
//   fig9:       H=5 c=10 phi=30 gamma=20 gammaw=9.81 hw=1 alphas=1,5,10 betas=15,20,..,90 (with 37.5, 52.5, 67.5, 82.5)
//               (default results/cpp/fig9.csv)
//   fembatch:   fig=9|8; FEM nref=3 mark=0.05 sa=2 (stability box sa H + H / tan(beta), capped at the hydraulic box)
//               sh0=0.25 shs=0.25 sgrade=0.25 shmax=1 maxnewton=100 tolfs=0.002 (E, nu: written to the settings when
//               not 20000, 0.3); limit analysis and seepage field as fig8 / fig9 (href=1 also for the FE field of the
//               FEM); cases: fig=9 H=5 c=10 phi=30 gamma=20 gammaw=9.81 hw=1 alphas= betas= curves=; fig=8
//               soil=swapped gammaw=9.8 panels= hws= curves= (limit analysis at H = 1 m, FEM at H_ref = H_crit of the
//               limit analysis, 3 digits, box 50 / 10 / 30 H_ref); out=<csv>
#include "Commands.h"

#include <iostream>
#include <string>

int main(int argc, char *argv[]) {
    using namespace slope;
    if (argc < 2) {
        std::cout << "usage: SlopeSeepageForces mesh|verify|seepage|fig5|fig8|fig9|analytical|probe|fs|check|la|labatch|fembatch "
                     "[key=value ...] (see main.cpp)\n";
        return 1;
    }
    Options o(argc, argv, 2);
    const std::string cmd = argv[1];
    if (cmd == "mesh") return CmdMesh(o);
    if (cmd == "verify") return CmdVerify(o);
    if (cmd == "seepage") return CmdSeepage(o);
    if (cmd == "probe") return CmdProbe(o);
    if (cmd == "fig5") return CmdFig5(o);
    if (cmd == "analytical") return CmdAnalytical(o);
    if (cmd == "fs") return CmdFS(o);
    if (cmd == "la") return CmdLA(o, argc, argv);
    if (cmd == "labatch") return CmdLABatch(o);
    if (cmd == "fig8") return CmdFig8(o);
    if (cmd == "fig9") return CmdFig9(o);
    if (cmd == "fembatch") return CmdFEMBatch(o);
    if (cmd == "check") return CmdCheck(o);
    std::cerr << "unknown command " << cmd << "\n";
    return 1;
}
