# YieldSurfaceProjection: Figs. 1 and 2 (MCC yield surface and closest-point projection)

This example reproduces Figs. 1 and 2 of the article *Return mapping for Modified Cam-Clay plasticity in rotated
Haigh-Westergaard space with consistent tangent operator and coupled u-p consolidation* (D. Lira Cecilio). Its
reference is the script `fig_surface.py` of the Python transcription (functions `fig_surface`, `cpp_linear`,
`ellipse` and `fig_meridian`).

It is not a finite element problem. It runs the two constitutive classes of the article and checks them:

- `TPZYCModifiedCamClayRHW`: the yield function (12), the hardening law (13) and the local Newton projection (18).
- `TPZPlasticStepModifiedCamClay`: the full stress update of Algorithm 1 and the consistent tangent.

The example also uses these native NeoPZ classes:

- `TPZHWTools`: Haigh-Westergaard cylindrical coordinates to principal stresses and to RHW Cartesian coordinates.
- `TPZTensor`: the tensors.
- `TPZElasticResponse`: the strain increment that produces a given trial stress.
- `TPZGeoMesh` and `TPZVTKGeoMesh`: the surface mesh and its VTK output.

## What is computed

**Fig. 1** shows the MCC yield surface (12) for M = 1, p'c = 100 kPa, pt = 0 and omega = 1, so a = 50 kPa.

- The surface is sampled on the grid of `fig_surface.py`:
  - 121 values of xi, from sqrt3·pt to -sqrt3·p'c;
  - 97 values of the Lode angle beta, from 0 to 2·pi;
  - the radius is rho_s(xi) = sqrt(2/3)·M·sqrt(a² - pbar²/b²).
- Principal stresses and RHW coordinates sigma* = (xi, rho·cos beta, rho·sin beta) are computed at each point.
- `TPZYCModifiedCamClayRHW::YieldFunction` is evaluated at every grid point. Phi/a² is at most 1.3e-15.
- The program prints the characteristic values from the caption: the length sqrt3·p'c, the width 2·rho_max and
  rho_max.
- The surface is also stored as a closed `TPZGeoMesh` and written to VTK:
  - quadrilaterals, with triangles at the two apexes;
  - material 1 for the subcritical part (pbar < 0) and material 2 for the supercritical part;
  - line elements with material 3 on the critical state circle;
  - the cell value is the deviatoric stress q.

**Fig. 2** shows the closest-point projection in the meridian plane.

- Model: linear elasticity with K = 6 MPa and G = 3 MPa, M = 1, p_c,n = 200 kPa, lambda = 0.2, kappa = 0.05 and
  v0 = 2, so v0/(lambda - kappa) = 13.3.
- Two trial states are projected: (p', q) = (190, 150) kPa in the subcritical region and (55, 128) kPa in the
  supercritical region.
- For each trial state:
  1. **Local problem (18).** It is solved with `TPZYCModifiedCamClayRHW::ProjectHW`, using xi_tr = -sqrt3·p'_tr,
     rho_tr = sqrt(2/3)·q_tr and b = 1, exactly as `cpp_linear` does. The program prints p', q, a, a_n, Delta
     alpha, Delta gamma and the number of Newton iterations next to the Python values.
  2. **Tangency check.** It checks that the energy-norm contour (p'-p'_tr)²/K + (q-q_tr)²/(3G) = d² through the
     projected state is tangent to the end-of-step ellipse. The sine of the angle between the two normals is
     printed.
  3. **Full stress update.** The trial state is passed through `TPZPlasticStepModifiedCamClay::ApplyStrainComputeSigma`.
     It starts from the converged state sigma_n = -100·I kPa, eps_n = 0 and p_c,n = 200 kPa. The total strain is
     eps = C⁻¹(sigma_tr - sigma_n), so the elastic predictor gives the Fig. 2 trial stress. Two kinds of trial
     state are used:
     - **Triaxial trial stress** (sigma_zz axial). The program prints the projected p', q and p_c, and six entries
       of the non-symmetric consistent tangent. They are compared with `camclay_hw.apply_strain` for the same
       trial stress.
     - **12 Lode angles** (beta = k·pi/6) with the same trial invariants, in a rotated frame. The program prints the
       number of failed projections, the largest deviations of p', q and p_c from the reference, and the largest
       change of the unit deviatoric direction n. The return is radial, so n is preserved.

## Build and run

The example is one of the targets of `Projects/CMakeLists.txt`
(`add_mcc_example(YieldSurfaceProjection main.cpp YieldSurfaceProjection.h)`), built when NeoPZ is configured with
`-DBUILD_PLASTICITY_MATERIALS=ON -DBUILD_PROJECTS=ON`:

```
cmake --build <build dir> --target YieldSurfaceProjection
cd <run dir> && <build dir>/Projects/YieldSurfaceProjection/YieldSurfaceProjection
```

The program takes no arguments and writes its files to the current directory. It runs in about 0.4 s.

The class `YieldSurfaceProjection` keeps the structure of the other examples: the data of the problem as members,
`CreateGeoMesh` (the closed surface mesh of Fig. 1), `RunSurface` (Fig. 1), `MeridianCases`, `ProjectMeridian`,
`StressUpdate`, `LodeAngleCheck` and `WriteMeridianFiles` (Fig. 2), `RunProjection` (printout of Fig. 2) and
`RunAll`, called by `main.cpp`.

## Output files

| File | Content |
|---|---|
| `fig1_surface.csv` | 121 x 97 grid of Fig. 1, beta outer loop and xi inner loop (the numpy meshgrid order of `fig_surface.py`). Columns: `xi, beta, rho, p_eff, q, sigma1, sigma2, sigma3` (principal stresses, tension positive), `sstar1, sstar2, sstar3` (RHW coordinates), `phi_over_a2` (Phi/a² from `TPZYCModifiedCamClayRHW`). |
| `fig1_critical_state_circle.csv` | Critical state circle, pbar = 0 (p' = p'c/2). Columns: `beta, sigma1..3, sstar1..3`. |
| `fig1_surface_principal.vtk` | Surface mesh in principal stresses with compression axes (-sigma1, -sigma2, -sigma3), as in Fig. 1a. |
| `fig1_surface_rhw.vtk` | Surface mesh in RHW space (xi, rho·cos beta, rho·sin beta), as in Fig. 1b. Cell fields: `ElData` = q (kPa), `Material` (1 subcritical, 2 supercritical, 3 critical state circle). |
| `fig2a_subcritical_curves.csv`, `fig2b_supercritical_curves.csv` | 400 points per curve. Columns: `p_start, q_start` (ellipse Phi_n = 0 at the start of the step), `p_end, q_end` (ellipse Phi = 0 at the end of the step), `p_energy, q_energy` (energy-norm contour through the projected state). |
| `fig2a_subcritical_points.csv`, `fig2b_supercritical_points.csv` | `point, p, q`. Point 0 is the trial state, point 1 the projected state, and point 2 the tip of the flow-direction arrow dPhi/dsigma (length 28 kPa, as in the figure). |
| `fig2_critical_state_line.csv` | Critical state line q = M·p' (two points). |

All curves match the arrays of `fig_surface.py` within the 12 significant digits written to the CSV files. The
largest absolute difference is 5e-10 kPa.

## Figures

```
python3 <neopz>/Projects/YieldSurfaceProjection/plot_figures.py [run directory] [-o output directory]
```

Run it after the executable, with the directory of its CSV files (default: the current directory). The figures
are written as PDF and PNG to `<run directory>/figures` (or to the output directory); Python 3 with numpy and
matplotlib is needed.

| File | Article | Data |
|---|---|---|
| `fig01_mcc_surface` | Fig. 1: MCC yield surface with the critical state circle, (a) principal stresses (compression axes), (b) RHW space | `fig1_surface.csv`, `fig1_critical_state_circle.csv` |
| `fig02_meridian_projection` | Fig. 2: closest-point projection in the meridian plane, (a) subcritical, (b) supercritical | `fig2a_subcritical_*.csv`, `fig2b_supercritical_*.csv`, `fig2_critical_state_line.csv` |

The script is a port of `fig_surface()` and `fig_meridian()` of `fig_surface.py` that reads the CSV files.

## Results compared with the article and the Python code

### Fig. 1

| Quantity | This code | Article (caption) |
|---|---|---|
| a = (p'c + pt)/(1 + omega) | 50 kPa | 50 kPa |
| Length sqrt3·p'c | 173.205 kPa | 173 kPa |
| Width 2·rho_max | 81.6497 kPa | 82 kPa |
| rho_max = sqrt(2/3)·M·a | 40.8248 kPa | 40.8 kPa |
| p' of the critical state section | 50 kPa | p'c/2 |
| q on the critical state circle | 50 kPa | M·a = 50 kPa (the CSL crosses the apex) |
| max \|Phi\|/a² on the grid | 1.3e-15 | 0 |

### Fig. 2: local projection

Python values come from `cpp_linear` in `fig_surface.py`. The two codes agree to machine precision: the relative
differences are at most 3e-16.

| Quantity | (a) this code | (a) Python | (a) article | (b) this code | (b) Python | (b) article |
|---|---|---|---|---|---|---|
| p' (kPa) | 163.661807029208 | 163.661807029208 | | 63.9281670721794 | 63.9281670721794 | |
| q (kPa) | 88.9952382809623 | 88.9952382809623 | | 91.9111041468245 | 91.9111041468245 | |
| a (kPa) | 106.027607010832 | 106.027607010832 | 106.0 | 98.0355153663691 | 98.0355153663691 | 98.0 |
| a_n (kPa) | 100 | 100 | 100 | 100 | 100 | 100 |
| Delta alpha | 4.38969882846531e-3 | 4.38969882846531e-3 | > 0 (expands) | -1.48802784536323e-3 | -1.48802784536323e-3 | < 0 (contracts) |
| Delta gamma | 3.80824131077182e-5 | 3.80824131077182e-5 | | 2.1813889378447e-5 | 2.1813889378447e-5 | |
| Newton iterations | 5 | 5 | | 5 | 5 | |
| sin(angle) of the energy contour and Phi = 0 | 3.2e-16 | | tangent | -2.7e-16 | | tangent |

### Fig. 2: full stress update

Python values come from `camclay_hw.apply_strain` with the same trial stress (p' and q are printed against the
`cpp_linear` values, which agree with those of `apply_strain` to the last bit).

| Quantity | (a) this code | (a) Python | (b) this code | (b) Python |
|---|---|---|---|---|
| p' (kPa) | 163.661807029208 | 163.661807029208 | 63.9281670721794 | 63.9281670721794 |
| q (kPa) | 88.9952382809623 | 88.9952382809623 | 91.9111041468245 | 91.9111041468245 |
| p_c (kPa) | 212.055214021664 | 212.055214021664 | 196.071030732738 | 196.071030732738 |
| Type | 1 (subcritical) | 1 | 2 (supercritical) | 2 |
| Local Newton iterations | 5 | 5 | 5 | 5 |
| D(xx,xx) | 6602.24726500131 | 6602.24726500131 | 5364.93535860097 | 5364.93535860097 |
| D(xx,yy) | 3042.43773376281 | 3042.43773376281 | 1056.60235171857 | 1056.60235171857 |
| D(xx,zz) | 2574.29639922714 | 2574.29639922714 | 4690.63762352755 | 4690.63762352755 |
| D(zz,xx) | 2901.65252916347 | 2901.65252916347 | 4952.36669843129 | 4952.36669843129 |
| D(zz,zz) | 2439.92318569832 | 2439.92318569832 | 6782.81618365257 | 6782.81618365257 |
| D(xy,xy) | 1779.90476561925 | 1779.90476561925 | 2154.1665034412 | 2154.1665034412 |
| 12 Lode angles, rotated frame: failures | 0 | | 0 | |
| 12 Lode angles: max \|Δp'\|, max \|Δq\|, max \|Δp_c\| (kPa) | 5.7e-14, 7.1e-14, 2.8e-14 | | 3.6e-14, 4.3e-14, 2.8e-14 | |
| 12 Lode angles: max \|n - n_tr\| | 1.3e-15 | | 7.0e-16 | |

The tangent entries agree to a relative difference of at most 3.4e-15. D(xx,zz) and D(zz,xx) differ from each
other: the consistent tangent of the MCC model is not symmetric (Sect. 3 of the article).
