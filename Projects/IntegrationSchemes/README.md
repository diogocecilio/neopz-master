# IntegrationSchemes: the return mapping against other integration schemes (Sect. 5.2 of the article, v1.0)

This example compares the backward-Euler return mapping of the article *"Backward-Euler integration of Modified
Cam-Clay in rotated Haigh-Westergaard space"* (D. Lira Cecilio; Sect. 5.2 of v1.0, Sect. 6.2 of v0.7) with other
integration schemes of the same model, at a material point, in tests with closed-form references: implicit
backward-Euler variants (this work with two and with four unknowns, the secant shear modulus of the recent implicit
schemes, the frozen bulk modulus of Sanei et al. 2020), the generalized trapezoidal rule, explicit Euler and RK4 with
fixed sub-steps, and the adaptive Runge-Kutta pairs ME2(1) and RKDP5(4). It measures accuracy against the number of
increments (order), accuracy against work, computing time, robustness against the size of the increment and the
convergence of the local Newton iterations. Up to v0.6 of the article the schemes of Table 4 existed only in the
Python code (`rivais.py`); the class `IntegrationSchemes` is their C++ transcription, extended in v1.0.

This is a material-point example: there is no finite element mesh and no VTK output.

## What is computed

All the schemes integrate the same Modified Cam-Clay model (omega = 1, p_t = 0) with the same laws: porous
volumetric elasticity p = p_n exp(-v0 Delta eps_v^e / kappa), shear modulus G = r K with
r = 3(1 - 2 nu) / (2 (1 + nu)), associated flow and the hardening law p_c = p_c,n exp(v0 Delta alpha / (lambda - kappa))
with a constant specific volume v0.

### Schemes

| scheme (name in the CSV files) | method | in the article |
|---|---|---|
| `library_exact` | the return mapping of this work: backward Euler (BE) in rotated Haigh-Westergaard space with the reduced local system of two unknowns (theta, Delta alpha), `TPZPlasticStepModifiedCamClay` (`LibraryUpdate`); the work is the number of local Newton iterations; the time includes the consistent tangent | this work |
| `library_full` | the same stress update with the four-unknown local system (xi, rho, Delta alpha, Delta gamma), option `TPZYCModifiedCamClayRHW::EFull` | four unknowns |
| `exact_n` | BE in tensors with the radial return (`BackwardEuler`, `be_tensor` of the Python code): the four-unknown system without the tangent; exact integral of the porous law, G = r K(p_n) kept in the step; same iterations as `library_full` | (check) |
| `exact_secant` | BE with the secant shear modulus G = r K_secant, K_secant = p_n (exp(-v0 x / kappa) - 1) / x, x = Delta eps_v + Delta alpha, which is the implicit system of Krabbenhoft and Lyamin 2012, Zhou et al. 2022, Lu et al. 2023 and Bui et al. 2026 at convergence; the increment is halved recursively only when the single step fails (`BackwardEulerSubstep`, strategy of Bui et al.) | secant G |
| `frozen_n` | BE with the bulk modulus frozen at its trial value in the plastic correction (Sanei et al. 2020); `library_frozen` is the same with the option `EFrozen` of the library | frozen K |
| `trapezoidal` | generalized trapezoidal rule of Ortiz and Popov (1985), theta = 1/2, for the flow and hardening rules (`GeneralizedTrapezoidal`): Delta eps^p = Delta gamma [(1 - theta) N_n + theta N_n+1], Delta alpha from the same weights, yield condition at the end; exact porous law, secant shear modulus of the step, so that the elastic law is integrated to second order; N_n is evaluated at the intersection of the elastic path with the yield surface when the step starts inside it; four unknowns, Newton with partial pivoting; theta = 1 reproduces `exact_secant` (single step) | trapezoidal rule |
| `trapezoidal_gn` | the same with G = r K(p_n) kept in the step (as this work); theta = 1 reproduces `exact_n` to round-off; first order (the deviatoric elastic law is then integrated only to first order within the step) | trapezoidal rule, G(p_n) |
| `FE` | explicit (forward) Euler with a fixed number of equal sub-steps per increment (Nayak and Zienkiewicz 1972; Sloan 1987 for the adaptive version), elastic part up to the yield surface, drift correction back to the surface after every sub-step (`RungeKuttaFixed` with the tableau `ForwardEuler`); the work is the number of evaluations of the elastoplastic operator (one per sub-step); the state must stay admissible (p' > 0, finite) and every drift correction must converge, otherwise the increment fails | explicit Euler |
| `RK4` | classical fourth-order Runge-Kutta with fixed sub-steps, same driver (four evaluations per sub-step) | RK4 |
| `ME2(1)`, `RKDP5(4)` | explicit adaptive Runge-Kutta pairs with local error control (`RungeKutta`; Sloan et al. 2001, Xie et al. 2026): exact integration of the elastic part, Newton search of the intersection with the yield surface, one drift correction at the end of the increment, step control of Sloan et al. (the sub-step grows by up to 1.1 after an accepted one) | adaptive RK |

The step control and model of Xie et al. 2026 (the next sub-step is the rest of the increment; void ratio updated)
are run in test B to validate the transcription against Fig. 3 of Xie et al. (control `xie` in the CSV files).

### Tests

| test (name in the CSV files) | material | path | error measured | runs |
|---|---|---|---|---|
| `xieB`: test B of Xie et al. 2026 | M = 0.896, lambda = 0.240, kappa = 0.045, v0 = 2.27 (e0 = 1.27), nu = 0.2, p'0 = p'c0 = 689.02 kPa | undrained (eps_v = 0) to eps_a = 10% | relative error of the stress, ‖sigma - sigma_exact‖ / ‖sigma_exact‖ | implicit schemes with 1, 2, 4, ..., 1024 increments; adaptive RK in one increment with tolerances 1e-1 to 1e-8 (and 1e-1 to 1e-5 with the control of Xie et al.); FE and RK4 with one sub-step over the same increments and with 2 to 1000 (FE) or 100 (RK4) sub-steps of the single increment |
| `kl_undrained_ocr1`, `kl_undrained_ocr10`: example 1 of Krabbenhoft and Lyamin 2012 | M = 3 sin(24°) / sqrt(3 + sin²(24°)) = 0.6858, lambda = 0.07, kappa = 0.008, v0 = 1.5, nu = 0.3, p'0 = 200 kPa | undrained, OCR = 1 to eps_a = 2.5% and OCR = 10 to 8% | p' - p'_exact and q - q_exact at the end | implicit with 1, 2, 4, 5, 10, ..., 1000 increments; adaptive RK in 10 increments with tolerances 1e-2 to 1e-7 (and with 1, 2, 4 increments, tolerance 1e-4); FE and RK4 with one sub-step over the increments and with 2 to 1000 (100) sub-steps of 10 increments |
| `kl_drained_ocr1`, `kl_drained_ocr10` | the same clay | drained, OCR = 1 and 10, to eps_a = 25% | q - q_exact and eps_v - eps_v,exact at the end, largest q | implicit with 1, 2, 4, 10, ..., 2000 increments; adaptive RK with 10 to 100 increments, tolerance 1e-4 (and 1e-5), and with 1, 2, 4; FE and RK4 as above |

The closed forms are those of Appendix B: the undrained path (B.7) with the elastic shear strain (B.9) of a
constant Poisson ratio (`UndrainedClosed`), solved for the stress ratio at which eps_q = eps_a reaches the final
strain by Brent's method (`Brent`, transcription of `scipy.optimize.brentq`), and the drained solution
(`mcc::TriaxialDrainedClosed` with 20000 points, interpolated at eps_a = 25%).

### Strain control and mixed control

* **Strain control** (undrained tests, `UndrainedPath`): all the components of the strain increment are prescribed,
  Delta eps_xx = Delta eps_yy = -Delta eps_zz / 2, so the material point follows a given strain path, as each
  integration point does within one global iteration of a finite element analysis. The observed orders of accuracy
  are those of the schemes.
* **Mixed control** (drained tests, `DrainedPath`): the lateral stress sigma_r = -p'0 is prescribed together with
  eps_zz, so the lateral strain of each increment is unknown and is found by equilibrium iterations (Newton with a
  finite-difference derivative and, if it fails, bisection in an interval where sigma_r + p'0 changes sign). The
  strain path inside the increment is then a straight line, a load-stepping error of second order that is common to
  all the schemes: with 10 increments every explicit scheme converges to the same error (q - q_exact = -2.01 kPa in
  the NC test and 1.59 kPa with OCR = 10) however many sub-steps it takes, and the second- and fourth-order schemes
  show their order only with many increments (`trapezoidal` and RK4 both reach order 2 beyond about 250 increments).

### Work and time

The work of a run is the sum, over its increments, of the local Newton iterations of the converged implicit updates
or of the evaluations of the elastoplastic operator in the stages of the explicit ones (the drift corrections are
counted apart, column `drift_iterations`); in the drained tests only the converged update of each increment is
counted (as in the Python code). `RunTimingStudy` measures the wall time per increment of every scheme in test B
(one increment), in the undrained test with OCR = 10 (10 increments) and in the drained test with OCR = 10 (10
increments, driver included): the path is repeated until it lasts about 0.5 s, five times, and the median is kept.
The library times include the spectral decomposition and the consistent tangent; `exact_n` is the bare
four-unknown Newton solve.

### Robustness, local convergence and the single step

* `RunRobustness` integrates the four tests of Krabbenhoft and Lyamin in 1, 2 and 4 increments with every scheme
  (one step per increment for the implicit ones; 1 and 10 sub-steps or the adaptive control for the explicit ones);
  the failures are listed by scheme.
* `RunSecantSweep` applies a single step of each implicit scheme (this work in tensors and in the library with two
  and with four unknowns, the secant shear modulus, the frozen bulk modulus, the two trapezoidal rules) to the first
  increment of the drained test with OCR = 10 and 10 increments (Delta eps_a = 2.5%) for 2001 lateral strains in
  [0.010, 0.030], and counts the failures and the spurious converged states (p' = 0).
* `RunLocalConvergence` records the norm of the residual at every local Newton iteration in four trial states:
  (a) NC clay of test B, undrained step of 10% (subcritical); (b) OCR = 10 clay at p' = 200, q = 400 kPa (just below
  the yield surface), undrained step of 2% (crosses the surface and softens, supercritical); (c) NC clay at the
  compression apex, isotropic step eps_v = -1% (rho_tr = 0, return to the apex); (d) the same with a shear strain of
  2e-7 (rho_tr/a_n = 5e-5); for this work with two unknowns (library), with four unknowns (tensors and library),
  the secant shear modulus and the trapezoidal rule.
* `RunToleranceStudy` repeats test B with the scheme of this work, 1, 4, 16 and 64 increments and tolerances of the
  local Newton iterations from 1e-2 to 1e-14 (default 1e-12): the error of a backward-Euler step is its truncation
  error, which only smaller increments reduce.

## How to run

The example is built with the other examples of the article (`-DBUILD_PLASTICITY_MATERIALS=ON -DBUILD_PROJECTS=ON`):

    cmake -G Ninja -B build -DBUILD_PLASTICITY_MATERIALS=ON -DBUILD_PROJECTS=ON -DCMAKE_BUILD_TYPE=Release
    ninja -C build IntegrationSchemes
    mkdir run_schemes && cd run_schemes && ../build/Projects/IntegrationSchemes/IntegrationSchemes

The program takes no arguments, writes its files to the current directory and runs in about 2 min (about 750 runs
in a few seconds; the timing study takes the rest). It prints the checks (theta = 1 of the trapezoidal rules
against backward Euler, the library against `be_tensor`, the comparison with the Python results of v0.6), the
observed orders, the robustness, the local convergence histories and the timings.

## Output files

| file | content |
|---|---|
| `schemes_results.csv` | One row per run: `test, scheme, family` (BE, library, explicit, RK), `control` (`-`, `fixed`, `sloan`, `xie`), `model, n, tol, converged`, the values `v0, v1, v2` (xieB: relative error; undrained: p' - p'_exact, q - q_exact; drained: q - q_exact, eps_v - eps_v,exact, largest q), `work, substeps, substepped_increments, rk_attempts, rk_rejections, p_end, q_end, eps_v_end, time_s`, for the drained tests `bisection_increments, unbalanced_increments, max_abs_sigma_r_plus_p0`, and `nsub` (sub-steps per increment of FE and RK4), `drift_iterations, message` (reason of a failure) |
| `schemes_error_work.csv` | The points of the figure error against work: `test, scheme, family, curve` (`n`: the increments vary, one step or sub-step per increment; `tol`: the tolerances of the adaptive RK at the increments of the test; `nsub`: the sub-steps at the increments of the test), `n, nsub, tol, converged, work, error, time_s, error2` |
| `schemes_order.csv` | Observed orders: least-squares slope of log(error) against log(n) over the four largest n with a solution (one step per increment), per test and scheme |
| `schemes_robustness.csv` | Per test, scheme and (fixed) sub-steps: increments run, with a solution, the smallest n that converged and the failed n |
| `schemes_local_convergence.csv` | `RunLocalConvergence`: `state, scheme, converged, iteration, residual, p_eff, q, pc` |
| `schemes_timing.csv` | `RunTimingStudy`: `test, scheme, n, nsub, tol, converged, error, work, microseconds_per_increment, spread` ((max - min)/median of the five repetitions) |
| `schemes_table_v1.csv` | The table of the article: one row per scheme with the errors, the work and the convergence in the four tests at the increments of the article (test B in one increment; the three tests of Krabbenhoft and Lyamin in 10) |
| `schemes_table4.csv`, `schemes_fig05.csv` | Table 4 and Fig. 6 of v0.7 (the schemes of v0.7 only) |
| `schemes_references.csv` | Closed-form references: `test, eps_a, eta, p_eff, q, eps_v, q_peak` |
| `schemes_library_check.csv` | `TPZPlasticStepModifiedCamClay` (two and four unknowns, frozen) against `be_tensor` for every BE run: final p', q, the largest difference along the path and the work of both |
| `schemes_python_comparison.csv` | Each result of `data_rivais.pkl` next to the C++ value and the difference, with the work of both |
| `schemes_secant_first_increment_sweep.csv` | `RunSecantSweep`: `eps_r, scheme, converged, iterations, spurious_p_zero, p_eff, q, sigma_r_plus_p0` |
| `schemes_newton_tolerance.csv` | `RunToleranceStudy`: `n, newton_tol, converged, error, work` |
| `schemes_paths_<test>.csv` | Paths of all the converged runs (`scheme, n, nsub, tol, increment, eps_a, p_eff, q, eps_v`) |

## Figures

```
python3 <neopz>/Projects/IntegrationSchemes/plot_figures.py [run directory] [-o output directory]
```

| File | Article | Data |
|---|---|---|
| `fig_error_work` | error against work of every scheme: (a) test B in one increment, (b) undrained OCR = 10, (c) drained NC, (d) drained OCR = 10 | `schemes_error_work.csv` |
| `fig_order` | error against the number of increments in the two strain-controlled tests with the slopes 1, 2 and 4 | `schemes_error_work.csv` |
| `fig_local_convergence` | residual norm against the local iteration in the four states of `RunLocalConvergence` | `schemes_local_convergence.csv` |
| `fig06_integration_schemes` | Fig. 6 of v0.7 (the schemes of Table 4) | `schemes_fig05.csv` |
| `supplementary_secant_first_increment` | a single step of this work and of the secant scheme in the first increment of the drained test with OCR = 10 | `schemes_secant_first_increment_sweep.csv` |

## Results (v1.0)

### Errors and work at the increments of the article

(a) test B in one increment of 10%: relative error of the stress; (b) undrained, OCR = 10, 8% in 10 increments:
p' - p'_exact (kPa); (c) drained, NC, 25% in 10 increments: q - q_exact (kPa), eps_v - eps_v,exact (%); (d) drained,
OCR = 10, 25% in 10 increments: the same. Work: local Newton iterations or operator evaluations.

| scheme | (a) error | work | (b) Delta p' | work | (c) Delta q | Delta eps_v | work | (d) Delta q | Delta eps_v | work |
|---|---|---|---|---|---|---|---|---|---|---|
| this work: BE in RHW, two unknowns | 8.04e-2 | 3 | -1.70 | 26 | -3.59 | -0.087 | 30 | 6.28 | 0.145 | 33 |
| BE in RHW, four unknowns | 8.04e-2 | 9 | -1.70 | 56 | -3.59 | -0.087 | 84 | 6.28 | 0.145 | 78 |
| BE, secant G (sub-stepped) | 8.12e-2 | 9 | -1.61 | 56 | -3.58 | -0.087 | 85 | 5.28 | 0.122 | 79 |
| BE, frozen K | 4.90e-2 | 9 | -28.7 | 56 | no solution | | | 7.41 | 1.92 | 78 |
| trapezoidal rule, secant G | 4.20e-1 | 9 | 0.048 | 56 | 0.261 | 0.0063 | 85 | no solution | | |
| trapezoidal rule, G(p_n) | 4.24e-1 | 9 | 0.045 | 56 | 0.252 | 0.0061 | 82 | -2.64 | -0.062 | 79 |
| explicit Euler, 1 sub-step | 2.27e-1 | 1 | -46.1 | 10 | 2.29 | 0.62 | 10 | -8.90 | 1.24 | 10 |
| explicit Euler, 10 sub-steps | 1.38e-2 | 10 | -4.25 | 100 | -2.02 | -0.010 | 100 | 0.91 | 0.113 | 100 |
| explicit Euler, 100 sub-steps | 1.39e-3 | 100 | -0.42 | 1000 | -2.01 | -0.045 | 1000 | 1.52 | 0.044 | 1000 |
| RK4, 1 sub-step | 1.17e-1 | 4 | -0.45 | 40 | 50.8 | 5.6 | 40 | 60.7 | 4.4 | 40 |
| RK4, 10 sub-steps | 7.39e-4 | 40 | -4.9e-5 | 400 | -2.01 | -0.048 | 400 | 1.59 | 0.037 | 400 |
| adaptive ME2(1), tol 1e-4 | 3.25e-5 | 288 | -0.044 | 252 | -2.01 | -0.049 | 734 | 1.60 | 0.037 | 612 |
| adaptive RKDP5(4), tol 1e-4 | 2.35e-7 | 150 | -0.010 | 72 | -2.01 | -0.049 | 432 | 1.59 | 0.037 | 306 |

The two- and the four-unknown systems give the same stresses (largest difference along all the paths 2.8e-6 kPa,
relative 1.5e-8, accumulated over up to 2000 increments) with 2.6 to 3 times fewer local iterations
(17151 against 26887 in the 51 runs of the library check; the four-unknown library run has exactly the iterations of
`be_tensor` in every run). The values of the v0.7 schemes are unchanged (Table 4 of v0.7 reproduced to the last
digit, `schemes_table4.csv`).

### Orders of accuracy (strain control)

Least-squares slopes over the four largest numbers of increments (`schemes_order.csv`):

| test | BE (this work) | BE secant G | BE frozen K | trapezoidal, secant G | trapezoidal, G(p_n) | explicit Euler | RK4 |
|---|---|---|---|---|---|---|---|
| test B, undrained NC | 1.01 | 1.01 | 0.98 | 2.00 | 1.24 | 0.99 | 3.96 |
| undrained, NC (K&L) | 1.00 | 1.00 | 1.01 | 2.00 | 1.19 | 0.98 | 3.02 |
| undrained, OCR = 10 | 1.07 | 1.07 | 0.98 | 2.00 | (0.12: sign change) | 1.00 | 3.99 |

Backward Euler is first order whatever the elastic law of the correction; the trapezoidal rule is second order only
when the deviatoric elastic law is integrated with the secant modulus (with G(p_n) its error curve bends to first
order, 1.24 in test B and 0.12 in the OCR = 10 test, where the error changes sign); RK4 is fourth order under strain
control. In the drained tests the mixed control limits every scheme to second order: with 250 to 4000 increments
the errors in q of the drained NC test are 1.51e-1 to 9.6e-3 kPa for backward Euler (order 1), 1.76e-3 to 7.7e-6 kPa
for the trapezoidal rule and 1.05e-2 to 3.3e-5 kPa for RK4 (both order 2).

### Robustness

* This work (two and four unknowns, single step per increment): every run converged, in all the tests from 1 to 2000
  increments (25% of axial strain in a single increment in the drained tests, 8% in the undrained one with OCR = 10),
  and all the 2001 single steps of the sweep (two unknowns: 3 to 4 iterations; four unknowns: 5 to 8).
* Secant G: the single step fails at 1255 of the 2001 lateral strains of the sweep and converges to a spurious
  state (p' = 0) at 14; with the sub-stepping of Bui et al. the scheme solves every test (two sub-steps in the first
  increment of the drained OCR = 10 test with 10 increments).
* Frozen K: no solution in the drained NC test with 1, 2, 4 and 10 increments (the bound Delta alpha < kappa/v0 of
  the frozen law); converges at the 2001 single steps of the sweep (5 to 11 iterations).
* Trapezoidal rule with the secant G: no solution in the drained OCR = 10 test with 1, 2, 4 and 10 increments and in
  the undrained OCR = 10 test with 1 and 2 increments; 277 of the 2001 single steps converge (33 to p' = 0,
  1691 fail). With G(p_n): no solution in the drained OCR = 10 test with 1, 2 and 4 increments; all the 2001
  single steps converge (6 to 8 iterations).
* Explicit Euler with one sub-step: no solution in the undrained OCR = 10 test with 1, 2 and 4 increments and in the
  drained OCR = 10 test with 4 (the sub-step leaves the admissible domain or the drift correction fails); with 10
  sub-steps, no solution with 1 increment in both drained tests. RK4 with one sub-step: no solution in the drained
  OCR = 10 test with 1 and 4 increments. The adaptive pairs solve every test (their step control sub-steps as needed).

### Local convergence (norm of the residual at each Newton iteration)

| state | two unknowns (this work) | four unknowns | four unknowns, secant G | trapezoidal rule |
|---|---|---|---|---|
| (a) subcritical, test B, 10% | 1.6e-2, 8.6e-5, 6.6e-9, 3.7e-17 (3 corrections) | 6.4e2 ... 8.6e-8, 2.6e-15 (9) | 9 | 9 |
| (b) supercritical, OCR = 10, 2% | 6.6e-3, 2.5e-5, 6.7e-9, 3.6e-16 (3) | 4.0, 0.92, 0.19, 1.7e-2, 5.6e-5, 6.9e-9, 1.6e-16 (6) | 6 | 7 |
| (c) apex, isotropic -1% | 1.0e-2, 1.7e-18 (1) | 1.4e2 ... 1.9e-11, 1.3e-15 (8) | 8 | 8 |
| (d) near the apex, rho_tr/a_n = 5e-5 | 1.0e-2, 4.2e-10, 4.0e-23 (2) | 8 | 8 | 8 |

Both systems converge quadratically at the end; the four-unknown iterations spend 4 to 6 corrections reducing the
residual from 1e2 to 1e0 before the quadratic regime, because they start from the trial state, far from the
surface, whereas the two-unknown iterations start on the surface at the angle of the trial state. The residuals
have different scales (strains for the reduced system; stresses scaled by a_n for the four-unknown one).

### Computing time per increment (microseconds, median of five; Release build, one core)

| scheme | test B, 1 increment | undrained OCR = 10, 10 increments | drained OCR = 10, 10 increments (driver included) |
|---|---|---|---|
| this work, two unknowns (library, tangent included) | 3.6 | 3.6 | 30 |
| four unknowns (library, tangent included) | 4.9 | 4.1 | 40 |
| four unknowns, tensors, no tangent (`exact_n`) | 1.5 | 1.1 | 12 |
| secant G (sub-stepped) | 1.8 | 1.3 | 24 |
| frozen K | 1.5 | 0.9 | 9 |
| trapezoidal, secant G / G(p_n) | 2.3 / 1.8 | 1.7 / 1.6 | - / 17 |
| explicit Euler, 1 / 10 / 100 sub-steps | 4.6 / 7.6 / 33 | 0.6 / 3.4 / 28 | 39 / 27 / 204 |
| RK4, 1 / 10 sub-steps | 5.0 / 9.8 | 0.8 / 6.3 | 9 / 51 |
| ME2(1) / RKDP5(4), tolerance 1e-4 | 49 / 25 | 4.7 / 1.4 | 91 / 39 |

The spread of the five repetitions is 5 to 50% on the shared machine of the runs: only ratios larger than about 1.5
are meaningful. The reduced system costs 15 to 25% less than the four-unknown one in the library (where the spectral
decomposition and the tangent are a large part of the time) and 2.3 to 2.9 times less in the bare local solve
(`Projects/YieldSurfaceProjection`, `two_vs_four_summary.csv`: 0.84 against 1.47 us with the linear law, 0.68
against 2.08 us with the porous law).

### Checks

* theta = 1 of the trapezoidal rule with G(p_n) against `BackwardEuler` (exact K, G(p_n)) in 205 steps (the four
  states above and 201 steps of the sweep): same convergence, same iteration counts, largest difference of the
  stresses 2.3e-13 kPa. With the secant G the comparison is against the single secant step, whose Newton
  iterations are chaotic where the step fails: 203 of the 205 steps agree, two differ in convergence.
* Library (two unknowns) against `be_tensor`: 51 runs, same convergence, largest difference of p', q along the paths
  2.8e-6 kPa (relative 1.5e-8); library with four unknowns against `be_tensor`: same iterations in every run,
  largest difference 2.3e-10 kPa.
* Python results of v0.6 (`data_rivais.pkl`, 188 results): the work counts and the convergence are identical in all
  of them and the values agree to round-off, except the drained OCR = 10 test with the secant G in 10 increments
  (q - q_exact = 5.278 kPa here against 5.316 in Python: the Python driver accepted an unbalanced first increment at
  a jump of sigma_r(eps_r) between branches of the sub-stepped secant step; see v0.7 of this README in the history).

## Files

* `IntegrationSchemes.h`: the class `IntegrationSchemes`: materials, backward-Euler schemes (`BackwardEuler`,
  `BackwardEulerSubstep`, `GeneralizedTrapezoidal`, `ElasticIntersection`), explicit schemes (`TTableau`,
  `TRKModel`, `RungeKutta`, `RungeKuttaFixed`), drivers (`UndrainedPath`, `DrainedPath`), closed forms
  (`UndrainedClosed`, `Brent`), update functions (`BEUpdate`, `RKUpdate`, `RKFixedUpdate`, `TrapezoidalUpdate`,
  `LibraryUpdate`), tests (`RunXieTestB`, `RunKLUndrained`, `RunKLDrained`, `RunRobustness`, `RunSecantSweep`,
  `RunToleranceStudy`, `RunLocalConvergence`, `RunTimingStudy`) and post-processing.
* `IntegrationSchemesReference.h`: the results of `data_rivais.pkl` (generated by `make_reference.py`; used only
  for the comparison).
* `make_reference.py`: writes `IntegrationSchemesReference.h` from `data_rivais.pkl` of the Python code.
* `main.cpp`: creates the example and calls `RunAll()`.
* `plot_figures.py`: the figures.
