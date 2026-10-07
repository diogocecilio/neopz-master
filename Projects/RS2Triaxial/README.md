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

1. It computes the closed form of Appendix B.1 (`mcc::TriaxialDrainedClosed`) with `npts = 600`, as `gen_data.py`
   does: 41 points on the elastic branch, then two sets of 600 stress ratios on the plastic branch, one uniform
   and one clustered near M.
2. It computes the material point solution (`mcc::TriaxialDrained`, the transcription of `TriaxialPointCC`)
   with 100, 200, 400, 800 and 1600 increments. eps_zz is prescribed in each increment. The lateral strains
   come from a Newton iteration on sigma'_xx = sigma'_yy = -p'0 that uses the xx-yy block of the consistent
   tangent of `TPZPlasticStepModifiedCamClay`.
3. It computes the error q_num - q_exact at eps_a = 20% (Table 2) with two definitions of q_exact:
   - **(a) interpolated:** `numpy.interp` on the 600-point closed-form table, as in `gen_data.py`. This is the
     definition used by the Python reference numbers.
   - **(b) exact:** bisection on the stress ratio eta of the parametric closed form, solving eps_a(eta) = 0.2
     (`RS2Triaxial::ExactDeviatoricStress`). The two values differ by less than 1e-3 kPa. Definition (b)
     agrees with the digits printed in Table 2 of the article (see below).
4. It finds the largest difference along the path, max |q - q_closed(eps_a)|, between the 400-increment
   solution and the closed form, interpolated at the numerical axial strains.
5. For OCR = 5, it reports the peak of q, the closed-form peak, q at 20% and the maximum compression eps_v
   before dilation.
6. It runs a finite element check. The 400-increment test is repeated with one Hex20-Hex8 u-p element
   (quadratic serendipity displacement, trilinear pore pressure) on the unit cube with 2 x 2 x 2 Gauss points,
   the same element as the FLAC3D tests (FLAC3DTriaxial). u_x = 0 on the face x = 0, u_y = 0 on y = 0 and
   u_z = 0 on z = 0 (symmetry planes); the cell pressure p'0 acts on the faces x = 1 and y = 1; the vertical
   displacement of the top z = 1 is controlled; p_w = 0 is prescribed at the eight vertices (drained). The check
   follows the usual NeoPZ sequence, in the same style as `footing.h`:
   - `CreateGeoMesh` builds the geometric mesh (`mcc::CreateUnitCubeMesh`, a boundary quadrilateral for the
     displacement condition and a coincident one for the pore pressure on each face).
   - `CreateCompMesh` builds the atomic displacement and pressure meshes and the multiphysics mesh with memory.
   - `RunFiniteElement` builds `TPZPoroElastoPlasticUPAnalysis` with `TPZSkylineNSymStructMatrix` and
     `TPZStepSolver` (ELU), solves the increments, and post-processes the results to VTK and CSV.

   The strain is homogeneous, so the element must reproduce the material point.

## Build and run

The example is built with the other examples of the article. `Projects/RS2Triaxial/CMakeLists.txt` declares
the target with `add_mcc_example(RS2Triaxial main.cpp RS2Triaxial.h)`, and the function `add_mcc_example` is
defined in `Projects/CMakeLists.txt`. Configure the top-level CMake with `-DBUILD_PLASTICITY_MATERIALS=ON
-DBUILD_PROJECTS=ON`, build the target `RS2Triaxial`, and run it from an empty directory:

```
cmake -G Ninja -B build -DBUILD_PLASTICITY_MATERIALS=ON -DBUILD_PROJECTS=ON -DCMAKE_BUILD_TYPE=Release
ninja -C build RS2Triaxial
mkdir run_rs2 && cd run_rs2 && ../build/Projects/RS2Triaxial/RS2Triaxial
```

The run takes about 5.5 s (Release, one core; 1.2 to 1.6 s for each finite element check). The program prints Table 2 and the values of Fig. 4 next to the Python and
article reference values, which are hard-coded in `RS2Triaxial::Cases()` and `RS2Triaxial::RunAll()`.

## Output files

| file | contents |
|---|---|
| `rs2_<case>_closed.csv` | closed form up to eps_a = 20%: `eps_a,p_eff,q,eps_v,eps_q,sigma_a` (dashed lines of Fig. 4) |
| `rs2_<case>_n400.csv` | material point, 400 increments: the same columns plus `q_closed`, `q_minus_q_closed` and `eps_v_closed` (closed form interpolated at eps_a); solid lines of Fig. 4 |
| `rs2_table2.csv` | Table 2: `err_interp_<case>` and `err_exact_<case>` for 100 to 1600 increments; the last row (increments = 0) holds q_exact |
| `rs2_<case>_fe.csv` | FE check: `eps_a,p_eff,q,eps_v,eps_q` of the element, `p_point,q_point,eps_v_point` of the material point with the same increments, `q_minus_q_point` and the residual `evaluations` of the increment |
| `rs2_fe_check.csv` | FE check summary, one row per case (`case` 0-3 in the order nc_nu, nc_g, ocr2, ocr5): q, p' and eps_v at 20% of the element and of the material point, largest differences along the path, spread of q over the 8 integration points, mean and largest evaluations per increment, global iterations, bisections, wall time, equations |
| `rs2_mesh_*.csv` | geometric mesh of the element (`mcc::WriteMeshCSV`: nodes, elements, faces with their boundary ids, edges), for the figure of the model |
| `rs2_<case>_fe.scal_vec.0.vtk` | FE check: nodal displacement and pore pressure at eps_a = 20% (NeoPZ's graphical mesh appends `.scal_vec.0` to the name `rs2_<case>_fe.vtk`) |
| `rs2_<case>_fe_gauss.vtk` | FE check: integration points (stress, p', q, p'c, v0, type of response) |

Fig. 4 plots q against `eps_q` in panels (a-d) and `eps_v` against `eps_a` in panels (e-h).

## Figures

```
python3 <neopz>/Projects/RS2Triaxial/plot_figures.py [run directory] [-o output directory]
```

Run it after the executable, with the directory of its CSV files (default: the current directory). The figures
are written as PDF and PNG to `<run directory>/figures` (or to the output directory); Python 3 with numpy and
matplotlib is needed.

| File | Article | Data |
|---|---|---|
| `fig04_rs2_triaxial` | Fig. 4: q against ε_q (a–d) and ε_v against ε_a (e–h) for NC with constant ν, NC with constant G, OCR = 2 and OCR = 5: this work with 400 increments, closed form and the digitized RS2 curves (Figs. 8.5–8.8 of the RS2 manual) | `rs2_<case>_n400.csv`, `rs2_<case>_closed.csv`, `reference/rs2_fig85_88_digitized.json` |
| `supplementary_rs2_element_check` | finite element check: (a) the Hex20–Hex8 element with the boundary conditions of the drained tests; (b) q against ε_a of the element (markers) and of the material point (lines); (c) \|q_FE − q_point\| along the path | `rs2_mesh_*.csv`, `rs2_<case>_fe.csv`, `rs2_fe_check.csv` |

The script also writes `rs2_digitized_comparison.csv` (output directory): the last points of the digitized RS2
curves against this work interpolated at the same abscissa. The end values of q differ from the RS2 analytical
curves by 0.03 to 0.20 %; the RS2 finite element curves of the normally consolidated cases are softer by 1.3
(constant G) and 2.2 % (constant ν) in q, with volumetric strains 6.5 and 7.9 % larger (text of Sect. 6.1).

The model of the element is drawn with the module `Common/mcc_hexmodel.py`; the figure of the model for the
article, with the boundary conditions of the drained (RS2 and FLAC3D) and undrained (FLAC3D) tests, is
`fig_single_element_model` of FLAC3DTriaxial.

## Results compared with the reference

The reference values come from the Python transcription (`gen_data.py rs2`) and from the article (Table 2, the
text of Sect. 6.1 and, for NC with constant G, Table 3).

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

One entry differs: NC with constant G at 800 increments gives -0.117497, which sits on the rounding boundary of
-0.118 (the program prints a note with the extra digits). With the interpolated q_exact, 8 entries and the q_exact
of NC with constant G (387.879) would differ from the article in the third decimal. The article's Table 2 is
therefore consistent with the exact closed-form value, while the Python reference numbers (`gen_data.py`) use the
interpolated one.

The largest error at eps_a = 20% with 100 increments is 0.24%, the same value as the article.

**Fig. 4, 400 increments.** Each cell shows this code | Python | article; `-` means the article gives no value.
The Python values are those of the arrays of `gen_data.py rs2` (`data_rs2.pkl`) with the same definitions. The
article gives the largest difference along the path in the text of Sect. 6.1 (1.4 kPa for the normally
consolidated clay, 8.5 kPa for OCR = 5) and, for NC with constant G, in Table 3 (1.26 kPa, closed form with
npts = 20000):

| quantity | NC, constant nu | NC, constant G | OCR = 2 | OCR = 5 |
|---|---|---|---|---|
| q(20%) (kPa) | 388.111817 \| 388.111817 | 387.645039 \| 387.645039 | 196.813852 \| 196.813852 | 202.926783 \| 202.926783 \| 202.9 |
| eps_v(20%) | 0.05055284 \| 0.05055284 | 0.05050174 \| 0.05050174 | 0.02244948 \| 0.02244948 | -0.01418197 \| -0.01418197 |
| max \|q - q_closed\| (kPa) | 1.3834 \| 1.3834 \| 1.4 | 1.2560 \| 1.2560 \| 1.26 | 2.7353 \| 2.7353 \| - | 8.4665 \| 8.4665 \| 8.5 |
| at eps_a | 0.85% \| 0.85% | 1.10% \| 1.10% | 0.30% \| 0.30% | 0.65% \| 0.65% |

**OCR = 5.** Each cell shows this code | Python | article:

| quantity | value |
|---|---|
| peak q, 400 increments | 292.9522 \| 292.9522 \| 293.0 kPa, at eps_a = 0.700 \| 0.700 \| 0.70% |
| closed-form peak | 293.3863 \| 293.3863 \| 293.4 kPa, at eps_a = 0.662 \| 0.662 \| 0.66% |
| q(20%) | 202.9268 \| 202.9268 \| 202.9 kPa |
| maximum compression eps_v | 0.2577% \| 0.2577% \| 0.26%, at eps_a = 0.70% |

**Local Newton iterations, 400 increments.** This includes every projection made inside the lateral-strain
iterations. Each cell shows this code | Python, plus the article (Table 3, NC with constant G) for the mean:

| case | mean iterations per projection | max | projections |
|---|---|---|---|
| NC, constant nu | 4.0708 \| 4.0708 | 5 \| 5 | 1215 \| 1215 |
| NC, constant G | 4.0511 \| 4.0511 \| 4.05 | 5 \| 5 | 1214 \| 1214 |
| OCR = 2 | 4.0000 \| 4.0000 | 4 \| 4 | 1182 \| 1182 |
| OCR = 5 | 4.0000 \| 4.0000 | 4 \| 4 | 1163 \| 1163 |

**Finite element check (one Hex20-Hex8 element, 2 x 2 x 2 points, 400 increments).** The element reproduces the
material point to the tolerance of the global Newton iterations (normalized residual 1e-8): the largest
differences along the path are 2.1e-6 kPa in q, 1.9e-6 kPa in p' and 1.2e-10 in eps_v, and q differs between the
eight integration points by less than 5e-7 kPa. The mean number of global residual evaluations per increment
includes the predictor; no increment was bisected. The values are those of the axisymmetric Q8-Q4 element of the
previous version to the 9 printed decimals.

| case | q(20%) FE (kPa) | q(20%) material point (kPa) | max \|q_FE - q_point\| (kPa) | evaluations per increment (largest) | total | time (s) |
|---|---|---|---|---|---|---|
| NC, constant nu | 388.111817677 | 388.111816880 | 2.1e-6 | 2.6425 (5) | 1057 | 1.6 |
| NC, constant G | 387.645039512 | 387.645038686 | 9.8e-7 | 2.6225 (5) | 1049 | 1.3 |
| OCR = 2 | 196.813852106 | 196.813851849 | 4.8e-7 | 2.4475 (5) | 979 | 1.2 |
| OCR = 5 | 202.926783192 | 202.926782950 | 5.3e-7 | 2.4075 (5) | 963 | 1.3 |
