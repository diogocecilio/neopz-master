# IntegrationSchemes: the return mapping against other integration schemes (Sect. 6.2, Table 4, Fig. 5)

This example compares the return mapping of the article *"Return mapping for Modified Cam-Clay plasticity in
rotated Haigh-Westergaard space with consistent tangent operator and coupled u-p consolidation"* (D. Lira Cecilio)
with the integration schemes of the recent literature for the same model, in material-point tests that have
closed-form references. It reproduces Table 4, Fig. 5 and the numbers of the text of Sect. 6.2. Up to version
v0.6 of the article these schemes existed only in the Python code (`rivais.py` and the function `rivais()` of
`gen_data.py`); the class `IntegrationSchemes` is their C++ transcription, and the return mapping of this work is
also run with the library class `TPZPlasticStepModifiedCamClay`.

This is a material-point example: there is no finite element mesh and no VTK output.

## What is computed

All the schemes integrate the same Modified Cam-Clay model (omega = 1, p_t = 0) with the same laws: porous
volumetric elasticity p = p_n exp(-v0 Delta eps_v^e / kappa), shear modulus G = r K with
r = 3(1 - 2 nu) / (2 (1 + nu)), associated flow and the hardening law p_c = p_c,n exp(v0 Delta alpha / (lambda - kappa))
with a constant specific volume v0.

### Schemes

| scheme (name in the CSV files) | method | in the article |
|---|---|---|
| `exact_n` | backward Euler (BE) in tensors with the radial return (`BackwardEuler`, `be_tensor`): exact integral of the porous law, G = r K(p_n) kept in the step | this work |
| `library_exact` | the same return mapping in rotated Haigh-Westergaard space, spectral form: `TPZPlasticStepModifiedCamClay` (`LibraryUpdate`) | this work (check) |
| `exact_secant` | BE with the secant shear modulus G = r K_secant, K_secant = p_n (exp(-v0 x / kappa) - 1) / x, x = Delta eps_v + Delta alpha, which is the implicit system of Krabbenhoft and Lyamin 2012, Zhou et al. 2022, Lu et al. 2023 and Bui et al. 2026 at convergence; the increment is halved recursively only when the single step fails (`BackwardEulerSubstep`, `be_substep`, strategy of Bui et al.) | secant G [10-13] |
| `frozen_n` | BE with the bulk modulus frozen at its trial value in the plastic correction (Sanei et al. 2020) | frozen K [26] |
| `library_frozen` | the same with the option `TPZYCModifiedCamClayRHW::EFrozen` of the library | frozen K (check) |
| `ME2(1)`, `RKDP5(4)` | explicit adaptive Runge-Kutta pairs with local error control (`RungeKutta`, `rk_update`; class `TRKModel` = `RKModel`): exact integration of the elastic part, Newton search of the intersection with the yield surface (with a bisection safeguard), one drift correction at the end of the increment, step control of Sloan et al. 2001 (the sub-step grows by up to 1.1 after an accepted one) | RK [14, 16] |

The step control and model of Xie et al. 2026 (the next sub-step is the rest of the increment; void ratio updated)
are run in test B to validate the transcription against Fig. 3 of Xie et al. (control `xie` in the CSV files).

### Tests (`rivais()` of `gen_data.py`)

| test (name in the CSV files) | material | path | error measured | runs |
|---|---|---|---|---|
| `xieB`: test B of Xie et al. 2026 | M = 0.896, lambda = 0.240, kappa = 0.045, v0 = 2.27 (e0 = 1.27), nu = 0.2, p'0 = p'c0 = 689.02 kPa | undrained (eps_v = 0) to eps_a = 10% | relative error of the stress, ‖sigma - sigma_exact‖ / ‖sigma_exact‖ | BE with 1, 2, 4, ..., 1024 increments; RK in one increment with tolerances 1e-1 to 1e-8 (and 1e-1 to 1e-5 with the control of Xie et al.) |
| `kl_undrained_ocr1`, `kl_undrained_ocr10`: example 1 of Krabbenhoft and Lyamin 2012 | M = 3 sin(24°) / sqrt(3 + sin²(24°)) = 0.6858, lambda = 0.07, kappa = 0.008, v0 = 1.5, nu = 0.3, p'0 = 200 kPa | undrained, OCR = 1 to eps_a = 2.5% and OCR = 10 to 8% | p' - p'_exact and q - q_exact at the end | BE with 5 to 1000 increments; RK in 10 increments with tolerances 1e-2 to 1e-7 |
| `kl_drained_ocr1`, `kl_drained_ocr10` | the same clay | drained, OCR = 1 and 10, to eps_a = 25% | q - q_exact and eps_v - eps_v,exact at the end, largest q | BE with 10 to 500 increments; RK with 10 to 100 increments, tolerance 1e-4 (and 1e-5, not in the Python results) |

The closed forms are those of Appendix B: the undrained path (B.7) with the elastic shear strain (B.9) of a
constant Poisson ratio (`UndrainedClosed`), solved for the stress ratio at which eps_q = eps_a reaches the final
strain by Brent's method (`Brent`, transcription of `scipy.optimize.brentq`), and the drained solution
(`mcc::TriaxialDrainedClosed` with 20000 points, interpolated at eps_a = 25%).

### Strain control and mixed control

* **Strain control** (undrained tests, `UndrainedPath`): all the components of the strain increment are prescribed,
  Delta eps_xx = Delta eps_yy = -Delta eps_zz / 2, so the material point follows a given strain path, as each
  integration point does within one global iteration of a finite element analysis.
* **Mixed control** (drained tests, `DrainedPath`): the lateral stress sigma_r = -p'0 is prescribed together with
  eps_zz, so the lateral strain of each increment is unknown and is found by equilibrium iterations (Newton with a
  finite-difference derivative and, if it fails, bisection in an interval where sigma_r + p'0 changes sign). The
  strain path inside the increment is then a straight line, an error of the load stepping that is common to all
  the schemes.

### Work

The work of a run is the sum, over its increments, of the local Newton iterations of the converged backward-Euler
updates or of the evaluations of the elastoplastic operator in the Runge-Kutta stages; in the drained tests only
the converged update of each increment is counted (as in the Python code).

### Tolerance of the local Newton iterations

`RunToleranceStudy` repeats test B with the scheme of this work, 1, 4, 16 and 64 increments and tolerances of the
local Newton iterations from 1e-2 to 1e-14 (default 1e-12): the error of a backward-Euler step is its truncation
error, which only smaller increments reduce.

### Single step per increment

The return mapping of this work uses a single backward-Euler step per increment, without sub-stepping. The program
counts the sub-steps of the secant scheme and, in addition (`RunSecantSweep`), applies a single step of this work
and of the secant scheme to the first increment of the drained test with OCR = 10 and 10 increments
(Delta eps_a = 2.5%) for 2001 lateral strains in [0.010, 0.030].

## How to run

The example is built with the other examples of the article (`-DBUILD_PLASTICITY_MATERIALS=ON -DBUILD_PROJECTS=ON`):

    cmake -G Ninja -B build -DBUILD_PLASTICITY_MATERIALS=ON -DBUILD_PROJECTS=ON -DCMAKE_BUILD_TYPE=Release
    ninja -C build IntegrationSchemes
    mkdir run_schemes && cd run_schemes && ../build/Projects/IntegrationSchemes/IntegrationSchemes

The program takes no arguments, writes its files to the current directory and runs in about 0.3 s (277 runs, the
sweep of the first increment and the tolerance study; Release build). It prints Table 4 and the numbers of the text of Sect. 6.2 next to the values of the
article (v0.6), the check of the library against `be_tensor` and the comparison with the Python results
(`IntegrationSchemesReference.h`, generated from `data_rivais.pkl`).

## Output files

| file | content |
|---|---|
| `schemes_results.csv` | One row per run: `test, scheme, family` (BE, library, RK), `control, model` (RK), `n, tol, converged`, the values `v0, v1, v2` (xieB: relative error; undrained: p' - p'_exact, q - q_exact; drained: q - q_exact, eps_v - eps_v,exact, largest q), `work, substeps, substepped_increments, rk_attempts, rk_rejections, p_end, q_end, eps_v_end, time_s`, and for the drained tests `bisection_increments` (increments solved by the bisection safeguard of the driver), `unbalanced_increments` (bisections stopped without equilibrium) and `max_abs_sigma_r_plus_p0` (largest \|sigma_r + p'0\| of the accepted increments, kPa; nan for the other tests) |
| `schemes_table4.csv` | Table 4: one row per scheme (and the two library checks): (a) error and work of test B in one increment; (b) p' and q errors and work of the undrained test with OCR = 10 in 10 increments; (c) q error, eps_v error (%) and work of the drained NC test in 10 increments, `c_converged` |
| `schemes_fig05.csv` | The points of Fig. 5: `panel` (a, b, c), `test, scheme, family, n, tol, converged, work, error` (the absolute value of the error drawn) |
| `schemes_references.csv` | Closed-form references: `test, eps_a, eta, p_eff, q, eps_v, q_peak` |
| `schemes_library_check.csv` | `TPZPlasticStepModifiedCamClay` (exact and frozen) against `be_tensor` for every BE run: final p', q, the largest difference along the path and the work of both |
| `schemes_python_comparison.csv` | Each result of `data_rivais.pkl` next to the C++ value and the difference, with the work of both |
| `schemes_secant_first_increment_sweep.csv` | `RunSecantSweep`: `eps_r, scheme, converged, iterations, spurious_p_zero, p_eff, q, sigma_r_plus_p0` (`iterations` 0 and nan values when the step did not converge in 50 iterations) |
| `schemes_newton_tolerance.csv` | `RunToleranceStudy`: `n, newton_tol, converged, error, work` |
| `schemes_paths_<test>.csv` | Paths of all the converged runs (`scheme, n, tol, increment, eps_a, p_eff, q, eps_v`) |

## Figures

```
python3 <neopz>/Projects/IntegrationSchemes/plot_figures.py [run directory] [-o output directory]
```

| File | Article | Data |
|---|---|---|
| `fig05_integration_schemes` | Fig. 5: error against work in the tests of Table 4, (a) test B of Xie et al., (b) undrained, OCR = 10, (c) drained, NC; port of `fig_rivais()` of `figs.py` | `schemes_fig05.csv` |
| `supplementary_secant_first_increment` | not in the article: a single step of this work and of the secant scheme in the first increment of the drained test with OCR = 10, (a) sigma_r + p'0 against the lateral strain, with the failures and the spurious states, (b) local iterations | `schemes_secant_first_increment_sweep.csv` |

## Results

### Table 4 (the values of the article v0.6 in brackets)

| Scheme | (a) error | work | (b) Delta p' (kPa) | work | (c) Delta q (kPa) | Delta eps_v (%) | work |
|---|---|---|---|---|---|---|---|
| This work: BE, exact K, G(p_n) | 8.04e-2 [8.04e-2] | 9 [9] | -1.70 [-1.70] | 56 [56] | -3.59 [-3.59] | -0.09 [-0.09] | 84 [84] |
| BE, exact K, secant G | 8.12e-2 [8.12e-2] | 9 [9] | -1.61 [-1.61] | 56 [56] | -3.58 [-3.58] | -0.09 [-0.09] | 85 [85] |
| BE, frozen K | 4.90e-2 [4.90e-2] | 9 [9] | -28.74 [-28.74] | 56 [56] | no convergence | | |
| RK ME2(1) | 3.25e-5 [3.25e-5] | 288 [288] | -0.044 [-0.044] | 252 [252] | -2.01 [-2.01] | -0.05 [-0.05] | 734 [734] |
| RK RKDP5(4) | 2.35e-7 [2.35e-7] | 150 [150] | -0.010 [-0.010] | 72 [72] | -2.01 [-2.01] | -0.05 [-0.05] | 432 [432] |

Every value of Table 4 and of Fig. 5 is the same as in the article: the tests are at a material point, so the
three-dimensional finite element models of v0.7 do not change them.

### Numbers of the text of Sect. 6.2

| quantity | this work | article v0.6 |
|---|---|---|
| test B, closed form at eps_q = 10% | p' = 393.0023 kPa, q = 351.3815 kPa | |
| control of Xie et al., tolerance 1e-4: ME2(1), RKDP5(4) | 4.30e-5, 2.87e-6 | 4.30e-5, 2.87e-6 |
| one increment: this work, secant G, frozen K | 8.04%, 8.12%, 4.90%, 9 iterations each | 8.0%, 8.1%, 4.9%, 9 |
| 128 increments, this work | 2.30e-4 in 511 iterations | 2.3e-4, 511 |
| RKDP5(4), tolerance 1e-3 / 1e-4 | 2.23e-6 with 126 / 2.35e-7 with 150 evaluations | 2.2e-6, 126 / 2.3e-7, 150 |
| ME2(1), tolerance 1e-4 | 3.25e-5 with 288 evaluations | 3.2e-5, 288 |
| drained NC, 10 increments, error in q: both RK, tolerance 1e-4 and 1e-5 | -2.01 kPa (all four) | 2.0 kPa |
| drained NC, 10 increments, BE | -3.59 kPa | 3.6 kPa |
| drained NC, 100 increments: RKDP5(4) / BE | -0.05 kPa with 666 evaluations / -0.37 kPa with 520 iterations | 0.05, 666 / 0.37, 520 |
| undrained OCR = 10, error in p': frozen K / exact, 10 increments | -28.74 / -1.70 kPa | 28.7 / 1.7 |
| undrained OCR = 10: frozen K with 1000 increments / exact with 50 | -0.38 / -0.11 kPa | 0.38 / 0.11 |
| drained NC, frozen K | no solution with 10 increments; with 20, eps_v 7.83% against 3.96% (+3.87 points) | the same |
| secant G against this work | errors changed by at most 2.5% in the NC tests, reduced by 1.9% to 20.1% in the OC tests | < 3%; 2 to 20% |
| secant G, drained OCR = 10, 10 increments | the first increment needs sub-stepping (two sub-steps; 11 sub-steps in the test) | the same |
| this work: single step per increment | 78 of 78 runs (`be_tensor` and library) | all the tests |
| drained tests (mixed control): equilibrium of the increments | all 90 runs with a solution: every increment solved by the Newton iterations of the driver, largest \|sigma_r + p'0\| 2.0e-8 kPa | (not reported) |

Single step in the first increment of the drained test with OCR = 10 (2001 lateral strains in [0.010, 0.030]):
the step of this work converges at all of them, in 5 to 8 local iterations; the single secant step fails (no
convergence in 50 iterations) at 1255, converges to the spurious state p' = 0 at 14 and converges at the other 732
(129 of them to another branch, more than 50 kPa away from the state of this work), with up to 50 iterations. The
single secant step converges for all the lateral strains up to 1.537%; above, it fails or ends at a spurious state
at 87% of the points (1269 of 1463). The solution of the increment is at eps_r = 1.52% with this work and at 1.65%
with the secant scheme, whose single step does not converge there: the sub-stepping of Bui et al. is needed.

Test B with this work and the tolerance of the local Newton iterations (relative error / local iterations):

| increments | 1e-2 | 1e-4 | 1e-8 | 1e-12 (default) | 1e-14 |
|---|---|---|---|---|---|
| 1 | 8.042e-2 / 7 | 8.041e-2 / 8 | 8.041e-2 / 9 | 8.041e-2 / 9 | 8.041e-2 / 9 |
| 4 | 1.353e-2 / 20 | 1.353e-2 / 22 | 1.353e-2 / 25 | 1.353e-2 / 28 | 1.353e-2 / 28 |
| 16 | 2.231e-3 / 48 | 2.237e-3 / 64 | 2.237e-3 / 80 | 2.237e-3 / 80 | 2.237e-3 / 84 |
| 64 | 3.835e-4 / 127 | 4.732e-4 / 191 | 4.734e-4 / 255 | 4.734e-4 / 256 | 4.734e-4 / 266 |

From 1e-8 to 1e-14 the error changes by less than 1e-9 in relative terms (at most 9.6e-10, with 64 increments;
with 1e-6 the change is at most 2.0e-6) (`schemes_newton_tolerance.csv`): a tighter tolerance only solves the
discrete equations more accurately, while the error is that of the first-order discretization of the increment.
With the loose tolerance 1e-2 the iterations stop before convergence and the error changes by chance (by 19% with
64 increments).

### The library against `be_tensor`

`TPZPlasticStepModifiedCamClay` (rotated Haigh-Westergaard space, spectral decomposition) and the tensor form
`be_tensor` give the same results in all the 39 BE runs of each variant (exact and frozen bulk modulus), with the
same numbers of local iterations in every run; the largest difference of p' or q along the paths is 5.5e-11 kPa
(relative 3e-13) with the exact integration and 3.7e-11 kPa with the frozen modulus.

### Comparison with the Python code

The program compares the 183 results of `data_rivais.pkl` (and the 5 closed-form references) with its own:

* the work counts (iterations or evaluations) and the convergence are identical in all of them;
* the closed-form references are identical (difference 0);
* the largest differences of the values are 5e-15 in the relative errors of test B, 6e-12 kPa in the undrained
  tests (p' about 830 kPa) and 2e-12 kPa in the drained tests (1.2e-11 kPa in the largest q), except in one run.

The exception is the drained test with OCR = 10 with the secant G in 10 increments: q - q_exact = 5.278 kPa and
largest q = 440.57 kPa here, against 5.316 kPa and 442.77 kPa in Python, with the same work (79) and two sub-steps
in the first increment in both codes. In that increment the iterations of the drained driver pass through lateral
strains where the single secant step fails or converges to spurious states (p' = 0) depending on round-off (see
the sweep above), so that sigma_r(eps_r) + p'0 of the sub-stepped update jumps between branches. Here the Newton
iterations of the driver converge, at eps_r = 0.016479, to an equilibrated state (|sigma_r + p'0| = 3e-14 kPa). In
the Python run (traced with `rivais.py` of v0.6) the Newton iterations reached a spurious p' = 0 state near
eps_r = 0.0164194 and diverged; the bisection safeguard then closed in on the jump at eps_r = 0.0164014 between the
sub-stepped branch (sigma_r + p'0 = -2.56 kPa) and a spurious single-step branch (+202 kPa), stopped after its 200
halvings without reaching the tolerance, and accepted the state with sigma_r + p'0 = -2.56 kPa, out of equilibrium
(1.3% of p'0); the largest q of the Python run, 442.77 kPa, is the q of that state. The value of this program is
therefore the solution of the drained increment and the Python value is an artefact of its driver. Both drivers
accept the last midpoint of a bisection that does not converge; this program counts such increments
(`unbalanced_increments` of `schemes_results.csv`) and reports the largest |sigma_r + p'0| of every drained run:
there are none, every drained increment is solved by the Newton iterations, and the largest |sigma_r + p'0| is
2.0e-8 kPa (the tolerance 1e-10 p'0). The run is not in Table 4 or Fig. 5; the statements of Sect. 6.2 that use it
hold with both values (the secant G reduces the error of this work by 16% here, 15% in Python).

Two details of the transcription make the comparison exact: the loading criterion of the Runge-Kutta schemes is
evaluated in closed form, a : (lambda tr(Delta eps) m + 2G Delta eps), which is exactly zero in the isochoric first
increment from the isotropic NC state, as the matrix product of the Python code (the increment then starts with an
elastic part of about 3e-7 of the increment); and Brent's method is the algorithm of `scipy.optimize.brentq`.

## Files

* `IntegrationSchemes.h`: the class `IntegrationSchemes`: materials, backward-Euler schemes (`BackwardEuler`,
  `BackwardEulerSubstep`), Runge-Kutta schemes (`TTableau`, `TRKModel`, `RungeKutta`), drivers (`UndrainedPath`,
  `DrainedPath`), closed forms (`UndrainedClosed`, `Brent`), update functions (`BEUpdate`, `RKUpdate`,
  `LibraryUpdate`), tests (`RunXieTestB`, `RunKLUndrained`, `RunKLDrained`, `RunSecantSweep`, `RunToleranceStudy`)
  and post-processing.
* `IntegrationSchemesReference.h`: the results of `data_rivais.pkl` (generated by `make_reference.py`; used only
  for the comparison).
* `make_reference.py`: writes `IntegrationSchemesReference.h` from `data_rivais.pkl` of the Python code.
* `main.cpp`: creates the example and calls `RunAll()`.
* `plot_figures.py`: Fig. 5 and the supplementary figure.
