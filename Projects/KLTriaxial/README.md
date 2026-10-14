# KLTriaxial: triaxial tests of Krabbenhoft and Lyamin (2012) with one Hex20-Hex8 element (Sect. 7.1, v1.0)

This example runs the triaxial tests of example 1 of Krabbenhoft and Lyamin (2012) — the clay of the
material-point comparison of the integration schemes (`IntegrationSchemes`) — with the return mapping of the
article in one Hex20-Hex8 u-p element of NeoPZ, and compares the element with the closed forms of Appendix B and
with the material point solution. Four tests: drained and undrained, normally consolidated (NC) and with OCR = 10.
The OCR = 10 clay yields in the supercritical region: in the drained test q peaks and the clay softens and dilates,
in the undrained test the excess pore pressure becomes negative.

## Material and tests

M = 3 sin(24°)/sqrt(3 + sin²(24°)) = 0.6858, lambda = 0.07, kappa = 0.008, v0 = 1.5, porous elasticity with
nu = 0.3, p_t = 0, omega = 1, p'0 = 200 kPa (cell pressure), p'c0 = 200 (NC) or 2000 kPa (OCR = 10).

| test (tag of the files) | control | final eps_a | closed form |
|---|---|---|---|
| `drained_nc` | sigma_r = -p'0 constant (total traction on the faces x = 1 and y = 1), eps_zz prescribed, pore pressure zero at the vertices | 25% | Appendix B.1 (`mcc::TriaxialDrainedClosed`, 20000 points) |
| `drained_ocr10` | the same | 25% | the same |
| `undrained_nc` | the same tractions, no drainage: the pore pressure is the Lagrange multiplier of the incompressibility (Biot alpha = 1, 1/Q = 0, permeability 0, dt = 0) | 2.5% | Appendix B.2 (`IntegrationSchemes::UndrainedClosed`, eq. (B.7) with the elastic shear strain (B.9)) |
| `undrained_ocr10` | the same | 8% | the same |

The element is the unit cube with 2 x 2 x 2 Gauss points, u_x = 0 on x = 0, u_y = 0 on y = 0, u_z = 0 on z = 0 and the
vertical displacement of the top controlled by the analysis (`TPZPoroElastoPlasticUPAnalysis`, skyline LU). Each
test is run with 10, 100 and 1000 increments; the material point solution (`mcc::TriaxialDrained`,
`mcc::TriaxialUndrained`) is run with the same increments. The excess pore pressure of the undrained tests is
computed from the stresses, u = p'0 + q/3 - p' (the nodal value of the pressure space is printed next to it).

## How to run

    ninja -C build KLTriaxial
    mkdir run_kl && cd run_kl && ../build/Projects/KLTriaxial/KLTriaxial
    python3 <neopz>/Projects/KLTriaxial/plot_figures.py

The program writes its files to the current directory and runs in about 20 s. With the argument `tangents`
(`../build/Projects/KLTriaxial/KLTriaxial tangents`, about 6 s) it runs only the drained tests with 100 increments
and the five tangent operators of `TPZPlasticStepModifiedCamClay::SetTangentMode` (D, fd, sym, cont, DT) and writes
`kl_tangents.csv` (column "Drained, OCR = 10" of Table 7 of the article v1.1): with OCR = 10, 3.05 evaluations per
increment (max 6) with D, fd, sym and DT and 9.52 (max 13) with the continuum tangent; on the NC clay 3.08 to 3.10
with D, fd, sym and DT and 1297 evaluations in all (four bisections) with the continuum tangent.

## Output files

| file | content |
|---|---|
| `kl_<case>_closed.csv` | closed form: `eps_a, p_eff, q, eps_v` (drained) or `eps_a, p_eff, q, u` (undrained) |
| `kl_<case>_fe_n<n>.csv` | element history with n increments: `eps_a, p_eff, q, eps_v, u, evaluations` (global residual evaluations of the increment), the material point `p_point, q_point, eps_v_point` or `u_point` and `q_minus_q_point` |
| `kl_table.csv` | one row per test and number of increments: errors of q, p' and eps_v at the end (element and material point), largest \|q_FE - q_point\| along the path, spread of q over the 8 points, global evaluations per increment (mean, largest), global iterations, bisections, wall time, local Newton iterations (mean, largest), closed-form values and the peak of q |
| `kl_<case>_fe.vtk`, `kl_<case>_fe_gauss.vtk` | last state of the 100-increment runs |

Figure (`figures/fig_kl_triaxial`, layout of v1.1): (a) q and (b) eps_v against eps_a of the drained tests, (c)
effective stress paths normalized by p'_c0 (the initial yield surfaces of both clays coincide) with the critical state
line, (d) and (e) q and the excess pore pressure against eps_a of the undrained tests, NC and OCR = 10, each with its
own scale, (f) the legend; lines are the closed forms and markers the element with 100 increments. The layout of v1.0
is kept in `plot_figures_v1.0.py`.

## Results

Errors at the end of the tests (element; the material point gives the same values to 1e-6 kPa):

| test | 10 increments | 100 | 1000 | closed form at the end |
|---|---|---|---|---|
| drained, NC: q - q_exact (kPa) | -3.587 | -0.370 | -0.0383 | q = 172.96 kPa, eps_v = 3.96% |
| drained, OCR = 10: q - q_exact (kPa) | 6.284 | 0.691 | 0.0697 | q = 190.03 kPa, eps_v = -5.15%; peak 537.7 kPa at eps_a = 0.85% |
| undrained, NC: p' - p'_exact (kPa) | 1.225 | 0.115 | 0.0115 | p' = 109.60, q = 74.11, u = 115.10 kPa |
| undrained, OCR = 10: p' - p'_exact (kPa) | -1.698 | -0.0472 | -0.00396 | p' = 831.92, q = 570.61, u = -441.72 kPa |

The errors fall by a factor of 9 to 12 per decade of increments (first order; 36 between 10 and 100 increments in the undrained OCR = 10 test, whose first increment crosses the surface), the largest difference between the
element and the material point is 1.4e-6 kPa, the spread of q over the eight integration points is at most 1e-7 kPa
(homogeneous deformation), and no increment is bisected: 2.0 to 4.7 global evaluations per increment (the
predictor included; 2.0 means that the predictor converges in one correction), 1.8 to 3.3 local Newton iterations
per projection on average (at most 4). The peak of q of the drained OCR = 10 test is 460.8 kPa with 10 increments
(at the end of the first increment, 2.5%), 531.9 kPa with 100 and 536.6 kPa with 1000, against 537.7 kPa at
0.853% in the closed form.

## Files

* `KLTriaxial.h`: the class `KLTriaxial` (cases, closed forms, material point, element, table).
* `main.cpp`: creates the example and calls `RunAll()`.
* `plot_figures.py`: the figure.
