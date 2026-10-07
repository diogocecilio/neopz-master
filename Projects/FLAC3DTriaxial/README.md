# FLAC3DTriaxial — Sect. 6.3, Fig. 6 and Table 5; FLAC3D column of Table 10 (Sect. 6.7)

Drained and undrained triaxial tests of the FLAC3D verification problem on a single three-dimensional u–p
element, from an isotropic effective stress p'0 = 5 kPa with overconsolidation ratios R = p'c0/p'0 = 1.6
(subcritical) and R = 8 (supercritical).

## Model

* One Hex20–Hex8 element (quadratic serendipity displacement, trilinear pore pressure) on the unit cube, the FLAC3D
  zone of 1 m, with 2 × 2 × 2 Gauss points (`SetIntegrationOrder(3)`).
* Boundary conditions: u_x = 0 on the face x = 0, u_y = 0 on y = 0 and u_z = 0 on z = 0 (symmetry planes); total
  cell pressure p'0 on the faces x = 1 and y = 1 (`ENeumannU`, traction −p'0 n); vertical displacement of the top
  z = 1 controlled (`SetControlledDisplacement(EZ1, 2)`, ε_a = −u_z).
* Drained tests: p_w = 0 prescribed at the eight vertices (a drained boundary element on each face), ε_a up to 50 %
  in 500 increments.
* Undrained tests: no flow (k = 0, Δt = 0, no drained face), pore fluid with K_w = 2·10⁴ kPa (1/M_B = n/K_w,
  n = (v0 − 1)/v0), ε_a up to 10 % in 400 increments.
* Material (Table 1): M = 1.02, λ = 0.2, κ = 0.05, v_λ = 3.32 (v0 on the NCL), G = 250 kPa, porous law.
* Sect. 6.7 (Table 10): drained test with R = 1.6 to ε_a = 5 % in 50 increments with the five tangent operators of
  `TPZPlasticStepModifiedCamClay::SetTangentMode`: consistent `D`, central differences `fd`, symmetric part `sym`,
  continuum tangent `cont` and transpose `DT`.

The state is homogeneous: the eight integration points have the same p' and q to 4·10⁻¹¹ kPa and the eight
vertices the same pore pressure to 2·10⁻¹⁰ kPa, and the element reproduces the axisymmetric Q8–Q4 element of v0.6
of the article (Python code `gen_data.py itasca`) to 10⁻¹² kPa.

The class `FLAC3DTriaxial` follows the structure of the NeoPZ examples: `CreateGeoMesh` (unit cube with
`mcc::CreateUnitCubeMesh` and `mcc::FaceOnPlane`), `CreateCompMesh` (displacement, pore pressure and
multiphysics meshes with the material `TPZMatPoroElastoPlasticUP` and its boundary conditions), `Run` (analysis
`TPZPoroElastoPlasticUPAnalysis` with `TPZSkylineNSymStructMatrix` and `TPZStepSolver` LU, increments, work
counters and post-processing), the closed forms of Appendix B (`DrainedClosed`, `UndrainedClosed`,
`CriticalState`, `ClosedPeak`, evaluated without interpolation by bisection on the stress ratio), `RunTests`
(Table 5, Fig. 6), `RunTangents` (Table 10) and `WriteMeshes` (model figure).

## Running

```
./FLAC3DTriaxial
```

Run time: 5.3 s (Release, one core; the incremental solutions take 1.4 s for each drained test, 0.6 s for each
undrained test and 0.16 to 0.25 s for each run of Table 10).

| File | Contents |
|---|---|
| `flac3d_<test>.csv` | history of the test (`drained_R1.6`, `drained_R8`, `undrained_R1.6`, `undrained_R8`): `eps_a, p_eff, q, v, u` (mean pore pressure of the vertices) and `evaluations` (residual evaluations of the increment) |
| `flac3d_<test>_closed.csv` | closed form up to the final ε_a (dashed lines of Fig. 6): `eps_a, p_eff, q, v` (drained) or `eps_a, p_eff, q, u` (undrained) |
| `flac3d_table5.csv` | Table 5: `test, quantity, this_work, closed_form, flac3d, critical_state` and the differences (%) from this work |
| `flac3d_summary.csv` | per test: final state, η at the end, peak of q and of the closed form, mean and largest evaluations per increment, global iterations, bisections, wall time, spread over the integration points, equations, differences from the Python final state, material |
| `flac3d_table10.csv` | Table 10, FLAC3D column: `tangent, nsteps, mean_evaluations, max_evaluations, global_iterations, bisections, wall_time_s, completed, p_end, q_end, v_end, python_mean_evaluations` |
| `flac3d_tangents_evaluations.csv` | evaluations of each of the 50 increments with each operator |
| `flac3d_drained_R1.6_50steps_D.csv` | history of the 50-increment test with D |
| `flac3d_mesh_<drained\|undrained>_*.csv` | geometric mesh of the element (`mcc::WriteMeshCSV`: nodes, elements, faces with their boundary ids, edges) |
| `flac3d_<test>.scal_vec.0.vtk`, `flac3d_<test>_gauss.vtk` | nodal fields and integration points at the end of the test |

## Figures

```
python3 <neopz>/Projects/FLAC3DTriaxial/plot_figures.py [run directory] [-o output directory]
```

Run it after the executable, with the directory of its CSV files (default: the current directory). The figures
are written as PDF and PNG to `<run directory>/figures` (or to the output directory); Python 3 with numpy and
matplotlib is needed. The module `Common/mcc_hexmodel.py` (drawing of the hexahedral models from the CSV files of
`mcc::WriteMeshCSV`) is also used by the scripts of RS2Triaxial and TerzaghiConsolidation.

| File | Article | Data |
|---|---|---|
| `fig_single_element_model` | model of the single-element tests: (a) drained tests of RS2 (Sect. 6.1) and FLAC3D (Sect. 6.3), p_w = 0 at the vertices; (b) undrained tests of FLAC3D (no drained face) | `flac3d_mesh_<drained\|undrained>_*.csv`, `flac3d_summary.csv` |
| `fig06_flac3d_triaxial` | Fig. 6: drained (a, b) and undrained (c, d) tests with R = 1.6 and R = 8, q–ε_a and stress paths p'–q with the CSL and the initial yield surface; dashed lines: closed-form solutions; squares: FLAC3D final states (Table 5) | `flac3d_<test>.csv`, `flac3d_<test>_closed.csv`, `flac3d_table5.csv`, `flac3d_summary.csv` |

The script also prints the final states, the peaks and the closed-form values at the same ε_a.

## Results

**Table 5** (final states, p', q and u in kPa; drained at ε_a = 50 %, undrained at 10 %). The values of this code
are those of v0.6 of the article (axisymmetric Q8–Q4 element, Python code) to 10⁻¹² kPa.

| Test | | this code | closed form | FLAC3D | critical state |
|---|---|---|---|---|---|
| Drained, R = 1.6 | p' / q / v | 7.573 / 7.720 / 2.811 | 7.573 / 7.720 / 2.811 | 7.573 / 7.718 / 2.811 | 7.576 / 7.727 / 2.811 |
| Drained, R = 8 | p' / q / v | 7.584 / 7.752 / 2.811 | 7.584 / 7.751 / 2.811 | 7.583 / 7.747 / 2.811 | 7.576 / 7.727 / 2.811 |
| Undrained, R = 1.6 | p' / q / u | 4.234 / 4.319 / 2.205 | 4.230 / 4.314 / 2.208 | 4.234 / 4.312 / 2.203 | 4.229 / 4.314 / 2.209 |
| Undrained, R = 8 | p' / q / u | 14.048 / 14.422 / −4.240 | 14.081 / 14.446 / −4.266 | 14.05 / 14.42 / −4.241 | 14.142 / 14.425 / −4.334 |

The closed forms agree with the element to 0.011 % in the drained tests and to 0.60 % in the undrained tests
(u with R = 8; the closed form assumes an incompressible fluid); FLAC3D agrees to 0.15 %. Peaks with R = 8:
q = 18.108 kPa at ε_a = 3.0 % (closed form 18.262 kPa at 2.93 %) in the drained test and 15.048 kPa at
ε_a = 4.0 % (closed form 15.051 kPa at 4.03 %) in the undrained test; η = 14.422/14.048 = 1.027 at the end of
the undrained test with R = 8. Residual evaluations per increment: 2.228 and 2.34 (drained), 2.0025 (undrained),
identical to the Python code increment by increment; no bisection.

**Table 10, FLAC3D column** (drained, R = 1.6, 50 increments to 5 %):

| Operator | evaluations per increment (largest) | total | wall time (s) | v0.6 (Python) |
|---|---|---|---|---|
| consistent D | 3.02 (4) | 151 | 0.16 | 3.02 |
| finite differences | 3.02 (4) | 151 | 0.19 | 3.02 |
| symmetric part (D + Dᵀ)/2 | 3.02 (4) | 151 | 0.16 | 3.02 |
| continuum tangent | 5.48 (7) | 274 | 0.25 | 5.48 |
| transpose Dᵀ | 3.02 (4) | 151 | 0.17 | 3.02 |

The final q coincide to 2.5·10⁻⁹ kPa. In the homogeneous test the Newton correction involves only the xx–yy block
of the operator, which is symmetric (σ_xx = σ_yy), so D, Dᵀ and (D + Dᵀ)/2 give the same iterations; the
evaluations of every increment are those of the Python code for the five operators.
