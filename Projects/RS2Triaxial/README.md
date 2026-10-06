# RS2Triaxial: drained triaxial tests of the RS2 manual (Sect. 6.1, Fig. 4, Table 2)

This example reproduces Sect. 6.1 of the article *Return mapping for Modified Cam-Clay plasticity in rotated
Haigh-Westergaard space with consistent tangent operator and coupled u-p consolidation* (D. Lira Cecilio):
Fig. 4 and Table 2. These are the drained triaxial compression tests of the RS2 manual (Rocscience, Sect. 8.7),
computed at a material point. The C++ code transcribes the function `rs2()` of `gen_data.py`, which calls
`camclay_hw.triaxial_point` and `camclay_hw.triaxial_closed` in the Python version.

## Problem

- Material (Table 1): M = 1.2, lambda = 0.077, kappa = 0.0066, v0 = 1.70. Elasticity is porous, with either a
  constant Poisson's ratio nu = 0.3 or a constant shear modulus G = 20 MPa. pt = 0 and omega = 1.
- Loading: the cell pressure stays equal to the initial isotropic effective stress p'0. The axial strain
  increases to eps_a = 20%.
- The four cases are named after the figures of the RS2 manual:

| case | tag | p'0 (kPa) | p'c0 (kPa) | shear law |
|---|---|---|---|---|
| Fig8.5, NC | `nc_nu` | 200 | 200 | constant nu |
| Fig8.6, NC | `nc_g` | 200 | 200 | constant G |
| Fig8.7, OCR = 2 | `ocr2` | 100 | 200 | constant nu |
| Fig8.8, OCR = 5 | `ocr5` | 100 | 500 | constant nu |

## What the program does

The class `RS2Triaxial` (`RS2Triaxial.h`) runs these steps for each case:

1. It computes the closed form of Appendix B.1 with 600 points per branch (`mcc::TriaxialDrainedClosed`).
2. It computes the material point solution (`mcc::TriaxialDrained`, the transcription of `TriaxialPointCC`)
   with 100, 200, 400, 800 and 1600 increments. eps_zz is prescribed in each increment. The lateral strains
   come from a Newton iteration on sigma'_xx = sigma'_yy = -p'0 that uses the xx-yy block of the consistent
   tangent of `TPZPlasticStepModifiedCamClay`.
3. It computes the error q_num - q_exact at eps_a = 20% (Table 2) with two definitions of q_exact:
   - **(a) interpolated:** `numpy.interp` on the 600-point closed-form table, as in `gen_data.py`. This is the
     definition used by the Python reference numbers.
   - **(b) exact:** bisection on the stress ratio eta of the parametric closed form, solving eps_a(eta) = 0.2
     (`RS2Triaxial::ExactDeviatoricStress`). The two values differ by less than 1e-3 kPa. Definition (b)
     reproduces the digits printed in Table 2 of the article.
4. It finds the largest difference along the path, max |q - q_closed(eps_a)|, between the 400-increment
   solution and the closed form, interpolated at the numerical axial strains.
5. For OCR = 5, it reports the peak of q, the closed-form peak, q at 20% and the maximum compression eps_v
   before dilation.
6. It runs a finite element check that is not in the article. The 400-increment test is repeated with one
   axisymmetric Q8-Q4 u-p element of 1 m x 1 m with 2 x 2 Gauss points. The radial displacement is fixed on the
   axis and the vertical displacement is fixed at the base. The cell pressure acts on the lateral face, the top
   displacement is controlled, and p = 0 is prescribed on the whole boundary (drained). The check follows the
   usual NeoPZ sequence, in the same style as `footing.h`:
   - `CreateGeoMesh` builds the geometric mesh.
   - `CreateCompMesh` builds the atomic displacement and pressure meshes and the multiphysics mesh with memory.
   - `RunFiniteElement` builds `TPZPoroElastoPlasticUPAnalysis` with `TPZSkylineNSymStructMatrix` and
     `TPZStepSolver` (ELU), solves the increments, and post-processes the results to VTK and CSV.

   The strain is homogeneous, so the element must reproduce the material point.

## Build and run

The example is built with the other projects (`add_mcc_example(RS2Triaxial main.cpp RS2Triaxial.h)` in
`Projects/CMakeLists.txt`). After building, run it from an empty directory:

```
mkdir run_rs2 && cd run_rs2
<build>/Projects/RS2Triaxial/RS2Triaxial
```

The run takes about 2 s. The program prints Table 2 and the values of Fig. 4 next to the Python and
article reference values, which are hard-coded in `RS2Triaxial::Cases()` and `RS2Triaxial::RunAll()`.

## Output files

| file | contents |
|---|---|
| `rs2_<case>_closed.csv` | closed form up to eps_a = 20%: `eps_a,p_eff,q,eps_v,eps_q,sigma_a` (dashed lines of Fig. 4) |
| `rs2_<case>_n400.csv` | material point, 400 increments: the same columns plus `q_closed`, `q_minus_q_closed` and `eps_v_closed` (closed form interpolated at eps_a); solid lines of Fig. 4 |
| `rs2_table2.csv` | Table 2: `err_interp_<case>` and `err_exact_<case>` for 100 to 1600 increments; the last row (increments = 0) holds q_exact |
| `rs2_<case>_fe.csv` | FE check: `eps_a,p_eff,q,eps_v,eps_q,q_point,q_minus_q_point` |
| `rs2_<case>_fe.scal_vec.0.vtk` | FE check: nodal displacement and pore pressure at eps_a = 20% |
| `rs2_<case>_fe_gauss.vtk` | FE check: integration points (stress, p', q, p'c, v0, type of response) |

Fig. 4 plots q against `eps_q` in panels (a-d) and `eps_v` against `eps_a` in panels (e-h).

## Results compared with the reference

The reference values come from the Python transcription (`gen_data.py rs2`) and from the article (Table 2 and
the text of Sect. 6.1).

**Table 2 (a), q_num - q_exact (kPa) with q_exact interpolated.** This code and Python agree to all 8 printed
decimals, so a single value is shown:

| increments | NC, constant nu | NC, constant G | OCR = 2 | OCR = 5 |
|---|---|---|---|---|
| 100 | -0.92711526 | -0.93228889 | -0.20723535 | 0.19603905 |
| 200 | -0.46534719 | -0.46766796 | -0.10393817 | 0.09867007 |
| 400 | -0.23319764 | -0.23427883 | -0.05187259 | 0.04912938 |
| 800 | -0.11660995 | -0.11715512 | -0.02585924 | 0.02455855 |
| 1600 | -0.05815191 | -0.05845798 | -0.01285719 | 0.01224430 |
| q_exact | 388.3450145 | 387.8793175 | 196.8657244 | 202.8776536 |

The 400-increment paths (all columns, 401 rows) and the closed-form tables agree with the Python arrays to the
12 significant digits written in the CSV files: the largest absolute difference is 5e-10 kPa.

**Table 2 (b), with q_exact solved exactly.** Each cell shows this code | article:

| increments | NC, constant nu | NC, constant G | OCR = 2 | OCR = 5 |
|---|---|---|---|---|
| 100 | -0.928 \| -0.928 | -0.933 \| -0.933 | -0.207 \| -0.207 | 0.196 \| 0.196 |
| 200 | -0.466 \| -0.466 | -0.468 \| -0.468 | -0.104 \| -0.104 | 0.099 \| 0.099 |
| 400 | -0.234 \| -0.234 | -0.235 \| -0.235 | -0.052 \| -0.052 | 0.049 \| 0.049 |
| 800 | -0.117 \| -0.117 | -0.117 \| -0.118 | -0.026 \| -0.026 | 0.025 \| 0.025 |
| 1600 | -0.059 \| -0.059 | -0.059 \| -0.059 | -0.013 \| -0.013 | 0.012 \| 0.012 |
| q_exact | 388.345 \| 388.345 | 387.880 \| 387.880 | 196.866 \| 196.866 | 202.878 \| 202.878 |

One entry differs: NC with constant G at 800 increments gives -0.11750, which sits on the rounding boundary.
With the interpolated q_exact, 8 entries and the q_exact of NC with constant G (387.879) would differ from the
article in the third decimal. The article's Table 2 therefore uses the exact closed-form value, while the Python
reference numbers use the interpolated one.

The largest error at eps_a = 20% with 100 increments is 0.24%, the same value as the article.

**Fig. 4, 400 increments.** Each cell shows this code | Python | article; `-` means the article gives no value:

| quantity | NC, constant nu | NC, constant G | OCR = 2 | OCR = 5 |
|---|---|---|---|---|
| q(20%) (kPa) | 388.111817 \| 388.111817 | 387.645039 \| 387.645039 | 196.813852 \| 196.813852 | 202.926783 \| 202.926783 \| 202.9 |
| eps_v(20%) | 0.05055284 \| 0.05055284 | 0.05050174 \| 0.05050174 | 0.02244948 \| 0.02244948 | -0.01418197 \| -0.01418197 |
| max \|q - q_closed\| (kPa) | 1.3834 \| 1.3834 \| 1.4 | 1.2560 \| 1.2560 \| - | 2.7353 \| 2.7353 \| - | 8.4665 \| 8.4665 \| 8.5 |
| at eps_a | 0.85% \| 0.85% | 1.10% \| 1.10% | 0.30% \| 0.30% | 0.65% \| 0.65% |

**OCR = 5.** Each cell shows this code | Python | article:

| quantity | value |
|---|---|
| peak q, 400 increments | 292.9522 \| 292.9522 \| 293.0 kPa, at eps_a = 0.700 \| 0.700 \| 0.70% |
| closed-form peak | 293.3863 \| 293.3863 \| 293.4 kPa, at eps_a = 0.662 \| 0.662 \| 0.66% |
| q(20%) | 202.9268 \| 202.9268 \| 202.9 kPa |
| maximum compression eps_v | 0.2577% \| 0.2577% \| 0.26%, at eps_a = 0.70% |

**Local Newton iterations, 400 increments.** This includes every projection made inside the lateral-strain
iterations. Each cell shows this code | Python:

| case | mean iterations per projection | max | projections |
|---|---|---|---|
| NC, constant nu | 4.0708 \| 4.0708 | 5 \| 5 | 1215 \| 1215 |
| NC, constant G | 4.0511 \| 4.0511 | 5 \| 5 | 1214 \| 1214 |
| OCR = 2 | 4.0000 \| 4.0000 | 4 \| 4 | 1182 \| 1182 |
| OCR = 5 | 4.0000 \| 4.0000 | 4 \| 4 | 1163 \| 1163 |

**Finite element check (one axisymmetric Q8-Q4 element, 400 increments).** The element reproduces the material
point to the tolerance of the global Newton iterations (normalized residual 1e-8). The mean number of global
residual evaluations per increment includes the predictor.

| case | q(20%) FE (kPa) | q(20%) material point (kPa) | max \|q_FE - q_point\| (kPa) | evaluations per increment |
|---|---|---|---|---|
| NC, constant nu | 388.111817677 | 388.111816880 | 2.2e-6 | 2.64 |
| NC, constant G | 387.645039512 | 387.645038686 | 9.8e-7 | 2.62 |
| OCR = 2 | 196.813852106 | 196.813851849 | 4.8e-7 | 2.45 |
| OCR = 5 | 202.926783192 | 202.926782950 | 5.3e-7 | 2.41 |
