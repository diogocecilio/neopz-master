# TerzaghiConsolidation — Sect. 6.4, Fig. 8 and Table 6

Consolidation of a linear elastic column (H = 10 m, E = 10⁴ kPa, ν = 0.25, k = 10⁻⁶ m²/(kPa·s),
c_v = kE_oed = 0.012 m²/s) loaded with q = 10 kPa on the drained top; incompressible constituents
(α_B = 1, 1/M_B = 0).

## Model

* Column [0, 1] × [0, 1] × [0, H] with 1 × 1 × 10 Hex20–Hex8 elements (quadratic serendipity displacement,
  trilinear pore pressure) and 3 × 3 × 3 Gauss points (`SetIntegrationOrder(4)`): 10 elements, 270 integration
  points, 428 equations.
* Base z = 0: u_z = 0, impermeable. Lateral faces x = 0, x = 1, y = 0, y = 1: zero normal displacement,
  impermeable. Top z = H: load q (traction −q e_z, `ENeumannU`) and drained (p_w = 0, `EDirichletP`).
* The load is applied in an undrained step (Δt = 0), followed by 102 time steps: 101 with 20 steps per decade from
  T = c_v t/H² = 10⁻⁵ to 1, plus T = 0.5 (T = 0.001, 0.01 and 0.1 are on the grid). The elastic response is given by
  `TPZPlasticStepModifiedCamClay::SetLinearElastic`.
* Monitored: the settlement of the top vertex of the vertical edge x = y = 0 and the pore pressure at the eleven
  vertices of this edge (z = 0, 1, ..., 10 m).

The solution is one-dimensional: the four top vertices settle by the same amount (to 2·10⁻¹⁸ m) and the four
vertical edges have the same pore pressures (to 8·10⁻¹⁴ kPa). The nodal values are those of the plane strain
1 × 10 Q8–Q4 column of v0.6 of the article (Python code `gen_data.py terzaghi`) to 6·10⁻¹⁵ m in the settlement and
5·10⁻¹¹ kPa in the pore pressures, and those of the Hex20–Hex8 column of the Python code (`gen_data3d.py`) to
10⁻¹⁰ kPa (the precision of the CSV files is 12 digits).

The class `TerzaghiConsolidation` follows the structure of the NeoPZ examples: `CreateGeoMesh`
(`mcc::CreateBoxMesh` with the boundary faces marked by `mcc::FaceOnPlane`), `CreateCompMesh` (displacement,
pore pressure and multiphysics meshes, material and boundary conditions), `Times`, `Run` (analysis
`TPZPoroElastoPlasticUPAnalysis` with `TPZSkylineNSymStructMatrix` and `TPZStepSolver` LU, time steps and
post-processing) and `RunAll` (Table 6, Fig. 8 and the exact series solution of Appendix B.3).

## Running

```
./TerzaghiConsolidation
```

Run time: about 2 s (Release, one core; 1.9 to 2.1 s in repeated runs).

| File | Contents |
|---|---|
| `terzaghi_history.csv` | `t`, `settlement` of the top and the pore pressures `p_z0` ... `p_z10` at the vertices of the edge x = y = 0 (initial state, undrained step and 102 time steps) |
| `terzaghi_isochrones.csv` | `T, z, pw_over_q, exact`: p_w/q at the vertices of x = y = 0 for T = 0.001, 0.01, 0.1 and 0.5 and the series solution (markers of Fig. 8b) |
| `terzaghi_degree.csv` | `T, U_numerical, U_exact`: degree of consolidation at the times of the analysis (solid line of Fig. 8c) |
| `terzaghi_exact_isochrones.csv` | `T, z, pw_over_q`: series solution at 201 heights (lines of Fig. 8b) |
| `terzaghi_exact_degree.csv` | `T, U`: series solution at 301 values of T from 10⁻⁵ to 1 (dashed line of Fig. 8c) |
| `terzaghi_table6.csv` | Table 6: `T, increment, max_err_pw, z_max_err, settlement_mm, settlement_exact_mm, settlement_diff_percent` |
| `terzaghi_summary.csv` | data of the problem, c_v, E_oed, w_∞, elements, points, equations, increments, evaluations per increment, global iterations, bisections, wall time, spread of the solution over the cross-section, p_w/q at z = 9 and 8 m and settlement after the undrained step |
| `terzaghi_mesh_*.csv` | geometric mesh of the column (`mcc::WriteMeshCSV`) for the model of Fig. 8a |
| `terzaghi.scal_vec.<step>.vtk` | nodal displacement and pore pressure at T = 0.001, 0.01, 0.1 and 0.5 (steps 42, 62, 82, 96) |

## Figures

```
python3 <neopz>/Projects/TerzaghiConsolidation/plot_figures.py [run directory] [-o output directory]
```

Run it after the executable, with the directory of its CSV files (default: the current directory). The figures
are written as PDF and PNG to `<run directory>/figures` (or to the output directory); Python 3 with numpy and
matplotlib is needed. The model is drawn with the module `Common/mcc_hexmodel.py`.

| File | Article | Data |
|---|---|---|
| `fig08_terzaghi_consolidation` | Fig. 8 (Fig. 7 of v0.6 with the model added): (a) the model (column of 1 × 1 × 10 Hex20–Hex8 elements, loaded and drained top with p_w = 0 at its vertices, lateral faces with u_n = 0 and impermeable, fixed impermeable base, monitored vertices of x = y = 0 and the top vertex); (b) p_w/q along x = y = 0 for T = 0.001, 0.01, 0.1 and 0.5 with the exact series solution; (c) degree of consolidation U = w/w∞ against T = c_v t/H² with the exact solution | `terzaghi_mesh_*.csv`, `terzaghi_isochrones.csv`, `terzaghi_degree.csv`, `terzaghi_exact_isochrones.csv`, `terzaghi_exact_degree.csv`, `terzaghi_summary.csv` |

The script also prints the largest error of p_w at each time.

## Results (Table 6)

| T | 0.001 | 0.01 | 0.1 | 0.5 |
|---|---|---|---|---|
| max \|p_w/q − exact\| (this code / v0.6) | 0.075 / 0.075 | 0.011 / 0.011 | 0.007 / 0.007 | 0.015 / 0.015 |
| settlement, this code (mm) | 0.366 | 0.954 | 2.959 | 6.289 |
| settlement, v0.6 (mm) | 0.366 | 0.954 | 2.959 | 6.289 |
| settlement, exact (mm) | 0.297 | 0.940 | 2.974 | 6.366 |

The errors and settlements differ from the Python values of v0.6 by less than 2·10⁻¹² (mm). Undrained step:
p_w/q = 1.268 at z = 9 m and 0.928 at z = 8 m (v0.6: 1.27 and 0.93), settlement 0.241 mm. At T = 0.5 the
settlement is 1.22 % below the exact value. Two evaluations of the residual per increment (linear problem), 206
global iterations in all, no bisection.
