# FrozenBulkModulus: exact versus frozen integration of the porous law (Sect. 6.1, Table 3)

This example reproduces Table 3 and the related results in the text of Sect. 6.1 of the article
*"Return mapping for Modified Cam-Clay plasticity in rotated Haigh-Westergaard space with consistent
tangent operator and coupled u-p consolidation"* (D. Lira Cecilio). The Python reference is the pair of
functions `frozen()` and `undrained_point()` in `gen_data.py`.

## What is computed

During the plastic correction, the volumetric equation of the local system (18) is solved in one of two
forms:

    exact (this work):        R1 = (xi - xi_tr exp(-v0 dal / kappa)) / a_n
    frozen bulk modulus:      R1 = (xi - xi_tr (1 - v0 dal / kappa)) / a_n      (Sanei et al. 2020)

The frozen form reverses the sign of `xi_tr` when `v0 dal / kappa > 1`. For `pt = 0` this limits the
growth of the preconsolidation pressure in one increment to `p'c / p'c,n < exp(kappa / (lambda - kappa)) = 1.098`.
The form is selected with `TPZPlasticStepModifiedCamClay::SetPorousIntegration(TPZYCModifiedCamClayRHW::EExact | EFrozen)`;
the example uses the library enumeration `TPZYCModifiedCamClayRHW::EPorousIntegration` directly. The rest of
the algorithm is the same for both forms.

* **Material.** RS2 clay: M = 1.2, lambda = 0.077, kappa = 0.0066, v0 = 1.70. Elasticity is porous, with a
  constant G = 20 MPa.
* **Initial states.** Two isotropic states: normally consolidated (NC, p'0 = p'c0 = 200 kPa) and OCR = 5
  (p'0 = 100 kPa, p'c0 = 500 kPa).
* **Increments.** n = 10, 15, 20, 25, 50, 100, 200, 400, 800 and 1600 increments, up to eps_a = 20%.
* **Drained test** (`mcc::TriaxialDrained`, constant cell pressure). The lateral strains come from a
  Newton iteration with the xx-yy block of the tangent. The program reports:
  * the error at the end, q_end - q_closed(eps_a,end);
  * the largest error along the path, max |q - q_closed(eps_a)|. Here q_closed is the closed form of
    Appendix B.1 (`mcc::TriaxialDrainedClosed` with 20000 points), interpolated linearly at the computed
    eps_a;
  * the mean and largest numbers of local Newton iterations per plastic projection (`mcc::TLocalStats`).

  A failed local projection means that the test has no solution, shown as a dash in Table 3. In that
  case the program also checks whether the first increment alone converges.
* **Undrained test** (`mcc::TriaxialUndrained`, with d eps_xx = d eps_yy = -d eps_zz / 2, so eps_v = 0).
  The program reports:
  * the largest |p' - p'(eta)| at the plastic states (|p' - p'0| > 1e-9 p'0), where p'(eta) is the
    closed-form path (B.7) at the same eta = q/p' (`mcc::UndrainedClosedP`, R = p'c0/p'0);
  * the error at the end, p'_end - p'0 (R/2)^Lambda, with Lambda = (lambda - kappa)/lambda;
  * the mean and largest numbers of local Newton iterations.
* **First increment of the drained NC test with 20 increments.** The program computes q, eps_v,
  p'c/p'c,n from the consistency condition (p'c = p' + q^2/(M^2 p')) and
  v0 dal/kappa = (lambda - kappa)/kappa ln(p'c/p'c,n) for both forms, then compares them with the
  closed form at eps_a = 1% and with the bound 1.098.
* **First increment of the frozen form with 15 and 10 increments** (`first_frozen` in `frozen()`).
  The program sweeps the lateral strain eps_r over the 4001 points of `numpy.linspace(-0.03, 0.01, 4001)`
  and applies the strain (eps_r, eps_r, -0.2/n) to the initial NC state in a single stress update. The
  change of sign of sigma_r + p'0 between two consecutive converged states brackets the solution of the
  increment. The program also interpolates linearly between the two states at sigma_r = -p'0, which gives
  the values quoted in the article for 15 increments (q = 41.9 kPa, eps_v = 3.5%).

This is a material-point example, so there is no finite element mesh and no VTK output. The class
`FrozenBulkModulus` follows the same sequence as the finite element examples:

1. set up the material (`CreateMaterial`) and the reference solutions (`ClosedForm`, `FirstIncrementClosed`);
2. solve (`RunDrained`, `RunUndrained`, gathered in `RunTests`, then `FirstIncrement` and `FirstIncrementSweep`);
3. post-process (`PostProcess` for the CSV files; `Print`, split into `PrintTable3`, `PrintOtherErrors`,
   `PrintIterations`, `PrintFirstIncrement` and `PrintAgreement`, for the comparison tables).

## How to run

The example is built together with the other examples of the article. Configure the top-level CMake with
`-DBUILD_PLASTICITY_MATERIALS=ON -DBUILD_PROJECTS=ON`, then build the target `FrozenBulkModulus`:

    cmake -G Ninja -B build -DBUILD_PLASTICITY_MATERIALS=ON -DBUILD_PROJECTS=ON -DCMAKE_BUILD_TYPE=Release
    ninja -C build FrozenBulkModulus
    mkdir run_frozen && cd run_frozen && ../build/Projects/FrozenBulkModulus/FrozenBulkModulus

The program takes no arguments and writes its files to the current directory. It runs in about 0.6 s
(Release build).

## Output files

| file | content |
|---|---|
| `frozen_table3.csv` | Table 3 for every n: `n`, largest error in q of the drained NC test (exact, frozen), mean local iterations (exact, frozen), largest error in p' of the undrained tests with the frozen modulus (NC, OCR5) and with the exact integration (NC, OCR5). `nan` means no solution. |
| `frozen_runs.csv` | One row per test (80 tests: drained or undrained, NC or OCR5, exact or frozen, n): convergence, end error, largest error, mean and largest local iterations, and the same values from the Python transcription. |
| `frozen_drained_<state>_<integration>.csv` | Paths of the drained tests for all n (column `n`): `eps_a, p_eff, q, eps_v, eps_q, sigma_a`, plus the closed-form `q_closed` at the same eps_a and `q_error`. |
| `frozen_undrained_<state>_<integration>.csv` | Paths of the undrained tests: `eps_a, p_eff, q, eta`, plus the flag `plastic`, `p_closed_B7` (the path (B.7) at the same eta) and `p_error`. |
| `frozen_closed_drained_<state>.csv` | Closed form of the drained test (B.1)-(B.6) with 600 points, for plotting. The errors are computed with 20000 points. |
| `frozen_closed_undrained_<state>.csv` | Closed-form undrained path: `eta, p_eff, q`. The elastic part (p' = p'0, from the origin to the yield ratio eta_y = M sqrt(R - 1)) and then (B.7) from eta_y to M (NC: eta from 0 up to M; OCR = 5: eta from 2.4 down to M). |
| `frozen_first_increment.csv` | First increment of the drained NC test with 20 increments (exact, frozen, closed form): `q, eps_v, pc_over_pcn, v0_dal_over_kappa`. |
| `frozen_first_sweep_n15.csv`, `frozen_first_sweep_n10.csv` | Sweep of eps_r in the first increment of the frozen form: `eps_r, converged, sigma_r_plus_p0, q, pc_over_pcn`. |
| `frozen_first_sweep_bracket.csv` | For n = 15 and 10: the two converged states around the change of sign and the interpolated state at sigma_r = -p'0 (`eps_r, sigma_r_plus_p0, q, pc_over_pcn, eps_v`). This is the `first_frozen` output of `frozen()`. |

`<state>` is `NC` or `OCR5`, and `<integration>` is `exact` or `frozen`.

## Figures

```
python3 <neopz>/Projects/FrozenBulkModulus/plot_figures.py [run directory] [-o output directory]
```

Run it after the executable, with the directory of its CSV files (default: the current directory). The figures
are written as PDF and PNG to `<run directory>/figures` (or to the output directory); Python 3 with numpy and
matplotlib is needed.

The article has no figure for this example (Table 3 only); the script draws a supplementary figure, labelled as
not in the article:

| File | Content | Data |
|---|---|---|
| `supplementary_table03_frozen_bulk_modulus` | exact integration against frozen bulk modulus versus the number of increments n (10 … 1600), NC and OCR = 5: (a) drained, largest error in q; (b) drained, mean local iterations; (c) undrained, largest error in p'; (d) undrained, mean local iterations. Table 3 holds the drained NC curves of (a) and (b) and the frozen curves of (c); crosses mark the frozen NC drained tests without solution (n ≤ 15) | `frozen_runs.csv` |

## Results

The program prints each value of this work with the reference in brackets. The reference is the
article, or the Python transcription (`reference_numbers.txt`, `gen_data.py frozen()`) for rows and
digits that the article does not print.

### Table 3

| n | drained NC, largest error in q (kPa), exact | frozen | local iterations, exact | frozen | undrained frozen, largest error in p' (kPa), NC | OCR = 5 |
|---|---|---|---|---|---|---|
| 10 | 31.07 (art. 31.07) | - (art. -) | 9.05 (9.05) | - (-) | 1.340 (1.340) | 5.161 (5.161) |
| 15 | 22.82 (Python 22.82) | - (Python -) | 7.90 (7.90) | - (-) | 1.237 (1.237) | 4.090 (4.090) |
| 20 | 18.15 (art. 18.15) | 94.06 (94.06) | 7.20 (7.20) | 8.98 (8.98) | 1.135 (1.135) | 3.301 (3.301) |
| 25 | 15.18 (Python 15.18) | 70.89 (70.89) | 6.78 (6.78) | 7.83 (7.83) | 1.037 (1.037) | 2.750 (2.750) |
| 50 | 8.47 (art. 8.47) | 28.11 (28.11) | 5.72 (5.72) | 5.86 (5.86) | 0.674 (0.674) | 2.383 (2.383) |
| 100 | 4.59 (art. 4.59) | 12.73 (12.73) | 5.16 (5.16) | 5.15 (5.15) | 0.401 (0.401) | 1.400 (1.400) |
| 200 | 2.43 (art. 2.43) | 6.23 (6.23) | 4.35 (4.35) | 4.33 (4.33) | 0.234 (0.234) | 0.769 (0.769) |
| 400 | 1.26 (art. 1.26) | 3.11 (3.11) | 4.05 (4.05) | 4.04 (4.04) | 0.127 (0.127) | 0.405 (0.405) |
| 800 | 0.64 (art. 0.64) | 1.56 (1.56) | 3.69 (3.69) | 3.65 (3.65) | 0.066 (0.066) | 0.208 (0.208) |
| 1600 | 0.32 (art. 0.32) | 0.78 (0.78) | 3.19 (3.19) | 3.12 (3.12) | 0.034 (0.034) | 0.106 (0.106) |

Every value agrees with the article to all the printed digits. With the frozen modulus, the drained NC
tests with 10 and 15 increments have no solution: the local projection fails already in the first
increment, as in the Python code and the article.

### Other results of Sect. 6.1

| quantity | this work | article | Python |
|---|---|---|---|
| undrained, exact integration: largest \|p' - p'(B.7)\|, n = 10..1600, NC / OCR5 (kPa) | 2.0e-11 / 6.5e-11 | < 1e-10 | 2.0e-11 / 6.5e-11 |
| drained NC: largest error of the frozen form / exact form | 2.41 to 5.18 | 2.4 to 5 | 2.41 to 5.18 |
| drained NC, n >= 50: difference in the mean local iterations | at most 2.43% | less than 3% | 2.43% |
| drained OCR5: largest error of the frozen form, larger by | 9.52% to 20.00% | 10 to 20% | 9.52% to 20.00% |
| drained OCR5: error at eps_a = 20% of the frozen form, larger by | 3.95% to 6.19% | 4 to 6% | 3.95% to 6.19% |
| undrained frozen, largest error in p' with 100 increments, NC / OCR5 (kPa) | 0.401 / 1.400 | 0.40 / 1.40 | 0.4006 / 1.3995 |

The program also prints, for every n, the drained end errors (NC and OCR5), the largest drained errors of
OCR5, the undrained end errors of the frozen form and the mean/largest local iterations of all the tests.
All of them agree with the Python values to the printed digits; the iteration counts are identical.

### First increment of the drained NC test with 20 increments (Delta eps_a = 1%)

| | q (kPa) | eps_v | p'c / p'c,n | v0 dal / kappa |
|---|---|---|---|---|
| closed form | 100.94 (art. 100.9) | 0.01209 (art. 1.21%) | | |
| exact | 85.00 (art. 85.0; Python 84.9955) | 0.00981 (Python 0.009805) | 1.2515 (art. 1.25) | 2.3931 (art. 2.39) |
| frozen | 41.82 (art. 41.8; Python 41.8229) | 0.02467 (art. 2.47%) | 1.0981 (art. 1.098) | 0.9981 (art. 0.998) |
| bound exp(kappa/(lambda - kappa)) | | | 1.0983 (art. 1.098) | |

### First increment of the frozen form with 15 and 10 increments (sweep of eps_r)

| n | sigma_r + p'0 changes sign between eps_r | q (kPa) | p'c / p'c,n | interpolated at sigma_r = -p'0 | failed projections |
|---|---|---|---|---|---|
| 15 | -0.01067 and -0.01065 (Python: same) | 41.7352 and 42.0285 (Python: same) | 1.098271 and 1.098270 (Python: same) | q = 41.88 kPa, eps_v = 0.0347 (art. 41.9 kPa and 3.5%) | 1836 of 4001 (Python 1847) |
| 10 | -0.01936 and -0.01523 (Python: same) | 10.3333 and 69.7412 (Python: same) | 1.098285 (Python: same) | q = 36.46 kPa, eps_v = 0.0551 | 2207 of 4001 (Python 2212) |

The solution of the increment lies at the bound p'c/p'c,n = 1.098. Near the bound the local projection
fails in about half of the trial states (with 10 increments, in all the 412 states between the bracketing
states).

The numbers of failed projections differ slightly from the Python run. The outcome differs in 41 states
(n = 15) and 19 states (n = 10) of the 4001, all with eps_r between -0.0199 and -0.0056 and trial mean
stresses between 6e5 and 2e8 kPa (the porous law is exponential). There the local Newton iterations of the
frozen form wander (up to 27 iterations when they converge), and the outcome depends on the last bits
of the trial invariants xi_tr and rho_tr. These come from the spectral decomposition of the trial stress
(`TPZPlasticStepModifiedCamClay::EigenSystem` in C++, `numpy.linalg.eigh` in Python), which differ by
round-off relative to that magnitude. When the Python trial invariants are given to
`TPZYCModifiedCamClayRHW::ProjectHW`, the C++ local Newton gives the same outcome and the same number of
iterations as the Python `newton()` in all 4001 states, for both n = 15 and n = 10 (checked outside the
example). The bracketing states, which are what the article uses, agree to all the printed digits.

### Agreement with the Python transcription

All 80 tests (drained and undrained, NC and OCR5, exact and frozen, 10 values of n) converge or fail
exactly as in the Python code. Over all the tests:

* the largest relative difference is 9.4e-11 for the end errors and 1.7e-11 for the largest errors;
* the undrained exact errors are at round-off level (absolute difference 1.9e-12 kPa);
* the mean and maximum numbers of local iterations are identical in every test.

## Files

* `FrozenBulkModulus.h`: the class `FrozenBulkModulus`. It contains the material set-up, the closed
  forms, the drained and undrained tests, the first-increment analyses, the CSV post-processing, the
  printout, the reference values of Table 3 (article) and the 80 Python reference values.
* `main.cpp`: creates the example and calls `RunAll()`.
* `CMakeLists.txt`: `add_mcc_example(FrozenBulkModulus main.cpp FrozenBulkModulus.h)`.

The example uses `TPZPlasticStepModifiedCamClay` and `TPZYCModifiedCamClayRHW` from
`Material/Plasticity`, and the material point drivers and closed forms of `Projects/Common/MCCPaperTools.h`.
