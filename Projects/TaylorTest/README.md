# TaylorTest: Taylor test of the consistent tangent (Sect. 4.5, Fig. 3)

This example checks that the tangent returned by the Modified Cam-Clay return mapping
(`TPZPlasticStepModifiedCamClay`, which projects in rotated Haigh-Westergaard space with
`TPZYCModifiedCamClayRHW`) is the consistent linearisation of the stress update. It reproduces
Sect. 4.5 and Fig. 3 of the article *"Return mapping for Modified Cam-Clay plasticity in rotated
Haigh-Westergaard space with consistent tangent operator and coupled u-p consolidation"*
(D. Lira Cecilio). The Python reference is the function `taylor()` in `gen_data.py`, and the figure
is drawn by `fig_taylor()` in `figs.py`.

## What is computed

The test compares the stress update with its first-order Taylor expansion, eq. (23):

    E(alpha) = || sigma(eps0 + alpha*deps) - sigma(eps0) - alpha*D*deps || = C*alpha^2 + O(alpha^3)
    p = (log E(alpha2) - log E(alpha1)) / (log alpha2 - log alpha1)  ~ 2

* Material: the clay of the Abaqus benchmark (Table 1): M = 1, lambda = 0.174, kappa = 0.026,
  v0 = 2.08, porous elasticity, nu = 0.3 (the shear modulus comes from the Poisson ratio).
* Every evaluation applies a total strain in one increment from an isotropic converged state with
  eps_n = 0 and p'_c,n = 116.6 kPa. The call is `mcc::ApplyStrain` from `Projects/Common`, which calls
  `TPZPlasticStepModifiedCamClay::ApplyStrainComputeSigma` and returns the 6x6 tangent:

  | kind | state | p'_n (kPa) | components of eps0 drawn in |
  |---|---|---|---|
  | 0 | elastic | 100 | [-0.002, 0.002] |
  | 1 | subcritical (compaction) | 100 | [-0.005, 0.005] |
  | 2 | supercritical (dilation) | 50 | [-0.02, 0.02] |

* For each kind, eps0 is drawn (six engineering Voigt components XX, XY, XZ, YY, YZ, ZZ) until the
  response is of that kind. The program then computes sigma0, D0 and the asymmetry
  ||D0 - D0^T|| / ||D0|| (Frobenius norms).
* For D = D0 and D = D0^T, the program draws 300 pairs of amplitudes alpha1, alpha2 ~ U(1e-4, 1e-2)
  and directions uniform in [-1, 1]^6, scaled to ||deps|| = 1e-3. Pairs whose perturbed states
  change kind are discarded. The program reports the least-squares slope of log E1 against
  log alpha1 over the 300 points and the median of the 300 pairwise slopes p.

### Random numbers

The random numbers are drawn in exactly the order used by `taylor()`, and the test runs twice:

1. **`TaylorTest::TNumpyRandom`** is a transcription of numpy's `default_rng(2026)`: SeedSequence,
   then PCG64 with XSL-RR output, then `Generator.uniform`. It uses portable 64-bit arithmetic. The
   program checks it against the first three raw outputs of numpy at start-up. Because it uses the
   same draws as the Python script, this run gives the same states, asymmetries and slopes as Fig. 3.
2. **`std::mt19937_64`** with seed 2026 gives an independent sample, which you can compare with the
   article only statistically. Its states (p', q) differ from those of Fig. 3. In the plastic states
   the slopes must be close to 2 with D and close to 1 with D^T.

The fitted slopes of the elastic panels have no reference in the article. They differ from 2 (Python
gives 2.04 and 2.28) because the elastic error E is very small: it comes only from the curvature of
the porous law, so it reaches round-off level at the smallest amplitudes. Only the median slope, 2.000,
is meaningful there.

## How to run

The example is built with the other examples of the article. The top-level CMake needs
`-DBUILD_PLASTICITY_MATERIALS=ON -DBUILD_PROJECTS=ON`; then build the target `TaylorTest`:

    cmake -G Ninja -B build -DBUILD_PLASTICITY_MATERIALS=ON -DBUILD_PROJECTS=ON -DCMAKE_BUILD_TYPE=Release
    ninja -C build TaylorTest
    mkdir run_taylor && cd run_taylor && ../build/Projects/TaylorTest/TaylorTest

The program takes no arguments, runs in well under a second and writes its files to the current
directory.

## Output files

| file | content |
|---|---|
| `taylor_pcg64_<kind>_<op>.csv` | Points of one panel of Fig. 3, with `<kind>` = `elastic`, `subcritical`, `supercritical` and `<op>` = `D`, `DT`. Columns: `alpha1, E1, log_alpha1, log_E1, fit_log_E1` (the least-squares line at log alpha1), `alpha2, E2, pair_slope`. Fig. 3 shows panels (a) `subcritical_D`, (b) `supercritical_D`, (c) `subcritical_DT` and (d) `supercritical_DT`. |
| `taylor_pcg64_summary.csv` | One row per panel: `kind, transposed, p_eff, q, pc, asym, fit_slope, fit_intercept, median_slope, min_slope, max_slope, rejected, draws`, then the state `x0_*` (engineering strains) and `sig0_*`. |
| `taylor_mt19937_*.csv` | The same files for the std::mt19937_64 sample. |

To redraw Fig. 3, plot `log_E1` against `log_alpha1` and the line `fit_intercept + fit_slope * log alpha`.
These are the same quantities that `fig_taylor()` plots from the Python results.

## Results

The program prints this work, then the article value in brackets and the Python value in braces.
The table below shows the run with the numpy-compatible stream (Fig. 3); the Python values come from
`reference_numbers.txt` (`gen_data.py taylor()`).

| state | op | p' (kPa) | q (kPa) | asym | fitted slope | median slope |
|---|---|---|---|---|---|---|
| elastic | D | 105.255 (Python 105.255) | 16.836 (16.836) | 0 (0) | 2.043056 (Python 2.043056) | 1.999994 (art. 2.000; Python 1.999994) |
| elastic | D^T | 105.255 | 16.836 | 0 | 2.284137 (Python 2.284155) | 1.999990 (art. 2.000; Python 1.999992) |
| subcritical | D | 103.412 (art. 103.4) | 40.704 (art. 40.7) | 0.01584 (art. 1.6%; Python 0.01584) | **1.996015** (art. 1.996; Python 1.996015) | 1.999987 (art. 2.000; Python 1.999987) |
| subcritical | D^T | 103.412 | 40.704 | 0.01584 | **1.039380** (art. 1.039; Python 1.039380) | 1.000271 (art. 1.000; Python 1.000271) |
| supercritical | D | 23.671 (art. 23.7) | 45.631 (art. 45.6) | 0.03288 (art. 3.3%; Python 0.03288) | **1.976463** (art. 1.976; Python 1.976463) | 2.000000 (art. 2.000; Python 1.999998) |
| supercritical | D^T | 23.671 | 45.631 | 0.03288 | **1.012999** (art. 1.013; Python 1.012999) | 0.999809 (art. 1.000; Python 0.999809) |

* The states sigma0 agree with the Python transcription to 12 or more digits: p' = 103.411561037478
  and q = 40.7038106469357 at the subcritical state, p' = 23.6707002558066 and q = 45.6311078924772 at
  the supercritical state. The asymmetries agree to 12 digits.
* The values of log E agree with those of the Python script to about 1e-8 (median over the 300
  points). The largest differences, about 1e-5, occur at the smallest amplitudes, where E reaches
  round-off level. As a result, the fitted slopes agree to 6 digits. Running the Python script on
  another machine changes them at the 1e-7 level.
* No pair was discarded, matching the Python run.

Independent sample with std::mt19937_64 (seed 2026). The states differ from Fig. 3, so compare only
the orders:

| state | op | p' (kPa) | q (kPa) | asym | fitted slope | median slope |
|---|---|---|---|---|---|---|
| elastic | D / D^T | 91.765 | 14.250 | 0 | 2.184 / 1.972 | 2.000010 / 2.000008 |
| subcritical | D / D^T | 76.154 | 55.985 | 0.01096 | 1.989 / 0.984 | 1.999997 / 1.000548 |
| supercritical | D / D^T | 56.637 | 58.024 | 0.08830 | 2.002 / 1.007 | 1.999998 / 1.000001 |

The conclusion is the same as in the article. With the consistent operator D the test is of second
order at both plastic states (median 2.000). With the transpose D^T, which is what the column-wise
Voigt assembly returns without the correction of Sect. 3, it is only of first order (median 1.000).

## Files

* `TaylorTest.h`: the class `TaylorTest`. It contains the material set-up, the drawing of the
  states, the perturbations and fits, the CSV post-processing, the printout with the reference values,
  and the generators `TNumpyRandom` and `TMersenneRandom`.
* `main.cpp`: creates the example and calls `RunAll()`.
