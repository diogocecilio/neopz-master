# AbaqusTriaxialConsolidation — Sects. 6.5 and 6.7, Figs. 9–11, Tables 7, 9 and 10

Abaqus benchmark 1.15.2: drained triaxial test with displacement control on a cylindrical clay specimen,
modelled in three dimensions with the mixed u–p element Hex20–Hex8 (serendipity quadratic displacements,
trilinear pore pressures). Modified Cam-Clay with porous elasticity: M = 1, λ = 0.174, κ = 0.026, v0 = 2.08,
ν = 0.3, p'c0 = 116.6 kPa, k = 2·10⁻¹⁰ m²/(kPa·s). Initial isotropic effective stress p'0 = 100 kPa equal to
the cell pressure; the platen moves down to δ/H = 0.6 in 400 days in 150 increments; free drainage at the top.

## Model

A quarter of the upper half of the specimen (H = 60 mm, R = 20 mm) is meshed with
`mcc::CreateQuarterCylinderMesh` (quadratic geometry, `TPZQuadraticCube`, the nodes of the lateral faces on the
cylinder): 2 × 2 × 4 divisions, i.e. 48 Hex20–Hex8 elements, 321 nodes, 95 pore pressure nodes and 1058
degrees of freedom, of which 811 (smooth platen) or 731 (rough platen) are free equations (Fig. 9a), and the
refined mesh 4 × 4 × 8 of the softening study, 384 elements, 2009 nodes and 6576 degrees of freedom, 5727 free
with the smooth platen (Fig. 9b). Boundary conditions: symmetry planes x = 0 (u_x = 0) and y = 0 (u_y = 0); mid-plane z = 0
impermeable with u_z = 0; cell pressure p'0 on the curved lateral face (radial traction −p'0 (x, y, 0)/r);
platen z = H with prescribed u_z and p_w = 0, smooth (only u_z) or rough (also u_x = u_y = 0). Integration:
2 × 2 × 2 (reduced) or 3 × 3 × 3 (full) Gauss points. The axisymmetric Q8–Q4 models of the previous version of
the article were removed.

Monitored quantities: stress at point A (x = 5 mm, y = 0, z = 7.5 mm, i.e. r = 5 mm, z = 7.5 mm as in the
benchmark) interpolated from the integration points of the element that contains it (`mcc::StressAtPoint`).
In the axisymmetric mesh of the benchmark A is the centroid of the element next to the axis; in the quarter of
cylinder it lies on the edge x = 5 mm, y = 0 of the element next to the axis (it is the mid-edge node of that
edge, Fig. 9a), and the first element that contains it is used, as in the Python code: its stress is
extrapolated from the integration points, linearly with 2 × 2 × 2 points and quadratically with 3 × 3 × 3 points.
The smallest, mean and largest q at the integration points of that element are also recorded. Further: average
axial stress on the platen σ_a (reaction divided by πR²/4) and the global deviatoric stress σ_a − p'0; largest
excess pore pressure; volumetric strain from the displacement of the outer top node (exact only for a
homogeneous deformation, i.e. the smooth platen from p'0 = 100 kPa); number of evaluations of the residual in
each increment.

## Parts and running

Parts (command line arguments, default `all`; `novtk` disables the VTK series):

| Part | Contents | Time |
|---|---|---|
| `mesh` | CSV files of the meshes 2 × 2 × 4 and 4 × 4 × 8 (`mcc::WriteMeshCSV`) for Fig. 9 | < 1 s |
| `mp` | material point: drained tests from p'0 = 100 and 20 kPa with 600 increments and the closed form of Appendix B.1 (Fig. 10); the benchmark state with 30, 150 and 3000 increments (column Δδ/H = 0.02 of Table 7) | 0.3 s |
| `fe` | 48 elements, 150 increments: smooth platen 2 × 2 × 2, rough platen 2 × 2 × 2 and 3 × 3 × 3 (Fig. 11, Table 7), comparison with the digitized Abaqus curves | 56 s (+ 98 s for the VTK series) |
| `tangents` | rough platen 2 × 2 × 2 with the operators D, central differences (`fd`), (D + Dᵀ)/2 (`sym`), continuum tangent (`cont`) and Dᵀ (`DT`): Tables 9 and 10 (no VTK series) | 3.2 min |
| `tolerance` | the five runs of `tangents` with the tolerance 10⁻⁹ instead of 10⁻⁸ (sensitivity of Table 10, see Results; no VTK series) | 4.0 min |
| `states` | 48 elements, smooth platen, 600 increments from p'0 = 100 and 20 kPa (Sect. 6.5) | 1.7 min (+ 55 s) |
| `softening` | p'0 = 20 kPa with the refined mesh 4 × 4 × 8, 600 increments (softening study) | 25 min (23.5 min + 1.2 min for the VTK series) |

```
./AbaqusTriaxialConsolidation                                             # everything (about 38 min)
./AbaqusTriaxialConsolidation mesh mp fe tangents tolerance states novtk  # all but the refined mesh, about 10 min
./AbaqusTriaxialConsolidation mesh mp fe novtk                            # Figs. 9 to 11 and Table 7, about 1 min
```

Times: Release build, one thread, wall times measured with other jobs on the machine (load average 2 to 4 on 4
cores; the parts `softening` and the others ran at the same time); the VTK series add about 0.2 s per state of
the 48-element mesh and 1.4 s per state of the refined mesh. The summary files also give the processor time of
each solution (`cpu_time_s`, without the monitored quantities and the VTK series), which does not depend on the
load of the machine and is the time used for Table 10.

The refined mesh has 5727 free equations, and the LU decomposition of the skyline matrix dominates its cost. With
the numbering of `TPZPoroElastoPlasticUPAnalysis` (pressure equations after the displacement equations, required
by undrained steps) one iteration takes 6.7 s; the part `softening` therefore renumbers the equations with the
bandwidth optimization of NeoPZ (`TConfig::fOptimizeBandwidth`, `TPZSloanRenumbering`) before the analysis is
built, which reduces the work of the factorization 4.4 times (1.1 to 1.9 s per iteration; the 600 increments,
1240 iterations, take 23 to 32 min instead of about 2.3 h). Every increment of the
benchmark has a positive time step, so the pressure block −Δt H is negative definite and the interleaved
pressure pivots are safe; the results are unchanged (largest difference 2·10⁻¹¹ kPa in a test with the rough
platen, same iterations). The other parts keep the default numbering, so that the times of Table 10 are those
of the standard set-up.

## Files

| File | Contents (part) |
|---|---|
| `abaqus_<run>.csv` | history of a finite element run: `delta_H`, `p_A`, `q_A`, `sigma_a_platen`, `q_platen` (σ_a − p'0), `max_pw`, `eps_v` (from the outer top node, meaningful only for a homogeneous deformation: smooth platen from p'0 = 100 kPa), `evaluations` (of the residual in the increment), `q_elA_min`, `q_elA_mean`, `q_elA_max` (q at the integration points of the element that contains A); runs `smooth_2x2x2`, `rough_2x2x2`, `rough_3x3x3` (`fe`), `rough_2x2x2_tangent_<D, fd, sym, cont, DT>` (`tangents`), `rough_2x2x2_tolerance_tangent_<D, fd, sym, cont, DT>` (`tolerance`), `smooth_600_p0_100`, `smooth_600_p0_20` (`states`), `smooth_600_p0_20_mesh4x4x8` (`softening`) |
| `abaqus_<run>_profile.csv` | final displacement of the outer generatrix (y = 0, r = R): `z`, `u_r`, `u_z` |
| `abaqus_<run>.scal_vec.0.vtk`, `abaqus_<run>_gauss.vtk` | final nodal fields and integration points |
| `abaqus_table7.csv` | Table 7: q at A at the abscissas of the digitized Abaqus curves (`delta_H_smooth`, `delta_H_rough`) and at 0.6 (`fe`) |
| `abaqus_fe_summary.csv`, `abaqus_fe_numbers.csv` | one row per run (platen, integration points, tangent operator, tolerance, sizes, degrees of freedom and free equations, evaluations, iterations, bisections, wall and processor times of the solution, wall times of the monitored quantities and of the VTK series, final and largest values); named numbers of Sect. 6.5 (RMS differences to Abaqus, locking of the full integration, ...) (`fe`) |
| `abaqus_table9.csv` | Table 9: `increment`, `delta_H`, `transposed` (0 D, 1 Dᵀ), `iteration`, normalized `residual` (`tangents`) |
| `abaqus_table10.csv` | Table 10: one row per operator (`tangent`, `tolerance`), with `mean_evaluations`, `max_evaluations`, `global_iterations` (`NGlobalIterations`), `bisections` (`NBisections`), `cpu_time_s` (processor time, the time of Table 10), `wall_time_s` and the final values (`tangents`) |
| `abaqus_tangents_evaluations.csv` | evaluations of each increment with each operator (`tangents`) |
| `abaqus_table10_tolerance.csv` | the rows of `abaqus_table10.csv` with the tolerance 10⁻⁹ (runs `rough_2x2x2_tolerance_tangent_<D, fd, sym, cont, DT>`, whose histories give the evaluations of each increment) (`tolerance`) |
| `abaqus_states_summary.csv`, `abaqus_states_numbers.csv`, `abaqus_softening_summary.csv`, `abaqus_softening_numbers.csv` | runs with 600 increments: peaks, final values, global q from the platen force, material point (`states`, `softening`) |
| `abaqus_material_point_*.csv`, `abaqus_closed_form_*.csv` | material point (`mp`) |
| `abaqus_mesh_<n>_{nodes,elements,faces,edges}.csv` | meshes 2 × 2 × 4 and 4 × 4 × 8 (`mesh`, format of `mcc::WriteMeshCSV`) |
| `vtk/<run>/` | VTK file series of the runs of `fe`, `states` and `softening` (see below) |

## Viewing the solution in ParaView

The runs of the parts `fe`, `states` and `softening` write their converged states with `mcc::TVTKSeries` in
`vtk/<run>/`: all the 151 states of the 150-increment runs, every fourth state of the 600-increment runs with 48
elements (151 states) and every twelfth state of the refined run (51 states). The file names start with
`abaqus_<run>`. The series take 324 MB in 2442 files (38 to 42 MB per run with 48 elements, 65 MB with 3 × 3 × 3 points and
100 MB for the refined mesh):

| File | Contents |
|---|---|
| `abaqus_<run>_nodal.vtk.series` → `abaqus_<run>_nodal.scal_vec.<k>.vtk` | `Displacement` (vector) and `PorePressure` (excess pore pressure, zero at the start) at the nodes, written by the NeoPZ graph mesh of the multiphysics mesh |
| `abaqus_<run>_intpoints.vtk.series` → `abaqus_<run>_intpoints.scal_vec.<k>.vtk` | the variables of the integration points of `TPZMatPoroElastoPlasticUP`, projected element by element on a discontinuous mesh by `TPZPostProcAnalysis`: `MeanEffectiveStress` p', `DeviatoricStress` q, `PreconsolidationPressure` p_c, `PlasticType` (0 elastic, 1 subcritical, 2 supercritical), `VolumetricStrain` (tr ε, negative in compression), `SpecificVolume`, the components `EffectiveStressXX/YY/ZZ/XY/XZ/YZ` and `TotalStressXX/YY/ZZ/XY/XZ/YZ`, `PrincipalEffectiveStress` (vector σ'1 ≥ σ'2 ≥ σ'3 at the points) and the tensors `EffectiveStress` and `TotalStress` |
| `abaqus_<run>_gausspoints.vtk.series` → `abaqus_<run>_gausspoints.<k>.vtk` | the integration points as a point cloud with their exact values: p', q, p_c, v0, type and σ' |
| `abaqus_<run>_states.csv` | index k of each state with δ/H, the time t (s) and the platen displacement u_c (m) |

The time of the series is δ/H (0 to 0.6). Stresses are positive in tension; p' and q follow the soil mechanics
convention (p' > 0 in compression). The projection has the order n − 1 for n × n × n Gauss points: trilinear
with the reduced 2 × 2 × 2 rule and triquadratic with the full 3 × 3 × 3 rule, so that its values at the element
vertices are the Lagrange extrapolation of the Gauss values (the interpolation used for the point A,
`mcc::StressAtPoint`). Being extrapolations, the vertex values can overshoot: the projected `PlasticType` is not
an integer (the point cloud has the type of each point), q can be slightly negative, and the projected principal
stresses, sorted at the points, can be out of order where two of them are close. For the principal values of
the extrapolated tensor apply *Filters → Alphabetical → Tensor Principal Invariants* (ParaView 5.10 or later)
to `EffectiveStress`.

In ParaView (5.5 or later, which reads the `.series` files):

1. *File → Open* `vtk/rough_2x2x2/abaqus_rough_2x2x2_nodal.vtk.series` (the `.series` file, not the numbered
   files) and press *Apply*.
2. Choose the field in the *Coloring* box of the toolbar and use *Rescale to data range over all
   timesteps* for a fixed color scale.
3. Play the increments with the VCR buttons of the *Time* toolbar (the time shown is δ/H).
4. Deformed shape: *Filters → Alphabetical → Warp By Vector* with *Vectors* = `Displacement` and
   *Scale Factor* 1 (the true shape: the platen moves down 36 mm on the 60 mm half specimen).
5. Integration point fields: open `vtk/rough_2x2x2/abaqus_rough_2x2x2_intpoints.vtk.series` and color by
   `DeviatoricStress`, `MeanEffectiveStress`, `PreconsolidationPressure`, ... To see them on the deformed mesh,
   apply *Filters → Resample With Dataset* (source: the nodal reader, destination: the integration point
   reader) and then *Warp By Vector*.
6. Point cloud: open `abaqus_rough_2x2x2_gausspoints.vtk.series`, set *Representation* to *Point Gaussian*
   (a *Gaussian Radius* of about 0.5 mm) and color by `PlasticType` or `DeviatoricStress`.
7. Whole specimen: *Reflect* with *Plane* X and then with *Plane* Y, keeping the input (*Copy Input*); a further
   *Reflect* with *Plane* Z gives the lower half.

The pore pressure is the excess pore pressure (the initial value is zero and gravity is not modelled). The test
is drained: the excess pore pressure is largest in the first increments and decreases to about 4·10⁻⁴ kPa at
the end.

## Figures

```
python3 <neopz>/Projects/AbaqusTriaxialConsolidation/plot_figures.py [run directory] [-o output directory]
```

Run it after the executable, with the directory of its CSV files (default: the current directory). The figures
are written as PDF and PNG to `<run directory>/figures` (or to the output directory); Python 3 with numpy and
matplotlib is needed. A figure whose CSV files are missing is skipped with a message. The figures use only the
CSV files (`./AbaqusTriaxialConsolidation mesh mp fe novtk` is enough for Figs. 9 to 11):

| File | Article | Data (part of the executable) |
|---|---|---|
| `fig09_abaqus_model` | Fig. 9: (a) quarter of the specimen with 48 Hex20–Hex8 elements, boundary conditions, vertex and mid-edge nodes and point A; (b) refined mesh 4 × 4 × 8 | `abaqus_mesh_2x2x4_*.csv`, `abaqus_mesh_4x4x8_*.csv` (`mesh`) |
| `fig10_abaqus_states` | Fig. 10: material point from p'0 = 100 and 20 kPa against the closed form: p'–q, q–ε1 and εv–ε1 | `abaqus_material_point_p0_*.csv`, `abaqus_closed_form_p0_*.csv` (`mp`) |
| `fig11_abaqus_results` | Fig. 11: q at A against δ/H with the smooth (a) and rough (b) platens and the stress paths at A (c), with the digitized Abaqus curves | `abaqus_smooth_2x2x2.csv`, `abaqus_rough_2x2x2.csv`, `abaqus_rough_3x3x3.csv` (`fe`), `abaqus_material_point_30.csv` (`mp`), `reference/abaqus_1_15_2_digitalizado.json` |
| `supplementary_abaqus_softening` | not in the article: p'0 = 20 kPa, global q (platen force) and q at A with the two meshes against the material point; deformed outer generatrix | `abaqus_smooth_600_p0_20*.csv` (`states`, `softening`) |
| `supplementary_abaqus_tangents` | not in the article: cumulative evaluations of the residual with the five operators of Table 10, (a) tolerance 10⁻⁸ and (b) 10⁻⁹ (panel (b) only if the part `tolerance` was run) | `abaqus_tangents_evaluations.csv`, `abaqus_table10.csv` (`tangents`); `abaqus_table10_tolerance.csv`, `abaqus_rough_2x2x2_tolerance_tangent_*.csv` (`tolerance`) |

## Results

All the values below are those of the files of the executable (`abaqus_table7.csv`, `abaqus_fe_numbers.csv`,
`abaqus_table9.csv`, `abaqus_table10.csv`, `abaqus_states_numbers.csv`, `abaqus_softening_numbers.csv`).

Table 7, q at A (kPa), at the abscissas of the digitized Abaqus curves (0.0297, 0.1301, ... for the smooth and
0.0304, 0.1294, ... for the rough platen); in parentheses the values of the previous version of the article
(axisymmetric 2 × 2 and 3 × 3 models; for the rough 2 × 2 × 2 column, the former 3D column):

| δ/H | Abaqus smooth | smooth 2 × 2 × 2 | material point Δδ/H = 0.02 | Abaqus rough | rough 2 × 2 × 2 | rough 3 × 3 × 3 |
|---|---|---|---|---|---|---|
| 0.03 | 59.7 | 60.2 (60.2) | 56.7 (56.7) | 59.7 | 61.5 (61.5) | 61.6 (61.6) |
| 0.13 | 109.9 | 114.0 (114.0) | 110.1 (110.1) | 111.8 | 116.4 (116.4) | 117.2 (117.4) |
| 0.22 | 131.1 | 133.8 (133.8) | 131.0 (131.0) | 134.7 | 136.9 (136.9) | 138.5 (138.7) |
| 0.32 | 141.1 | 143.5 (143.5) | 141.9 (141.9) | 145.0 | 146.9 (146.9) | 146.9 (147.4) |
| 0.41 | 145.8 | 147.1 (147.1) | 146.2 (146.2) | 150.1 | 151.6 (151.6) | 147.7 (149.0) |
| 0.51 | 148.1 | 148.8 (148.8) | 148.4 (148.4) | 153.1 | 154.5 (154.5) | 145.2 (147.7) |
| 0.60 | — | 149.5 (149.5) | 149.3 (149.3) | — | 155.9 (155.9) | 141.8 (145.5) |

| Quantity (150 increments, 48 elements) | this code | previous version / reference |
|---|---|---|
| End q_A: smooth 2 × 2 × 2 / rough 2 × 2 × 2 / rough 3 × 3 × 3 (kPa) | 149.502 / 155.938 / 141.801 | 149.50 (2D and 3D) / 155.94 (3D; 2D 155.56) / 2D 3 × 3: 145.5 |
| Smooth platen against the material point with 150 increments | largest difference 0.064 kPa (δ/H = 0.008, consolidation transient), 0.036 kPa for δ/H ≥ 0.1 | — |
| RMS difference to Abaqus: smooth FE / material point 30 increments / 150 increments | 2.27 / 1.08 / 2.29 kPa | 2.3 / 1.1 kPa |
| Largest difference to Abaqus, smooth platen | 4.08 kPa at δ/H = 0.13 | 4 kPa at 0.13 |
| RMS difference to Abaqus, rough 2 × 2 × 2 / 3 × 3 × 3 | 2.38 / 4.08 kPa | — |
| Full integration: largest q_A; first δ/H with \|q_A(3×3×3) − q_A(2×2×2)\| > 1 kPa; first δ/H with q_A(3×3×3) below q_A(2×2×2) | 147.93 kPa at δ/H = 0.38; 0.144; 0.32 | 2D: 149.0 at 0.40; 0.172; 0.332 |
| q_A at δ/H = 0.5096, rough 2 × 2 × 2 (Abaqus 153.09) | 154.50 | 2D 154.3 |
| Average axial stress on the platen at the end: smooth / rough 2 × 2 × 2 / rough 3 × 3 × 3 (kPa) | 249.49 / 250.78 / 250.93 | 249.5 / 250.8 / 2D 251.0 |
| Largest excess pore pressure at the end, rough 2 × 2 × 2 (kPa) | 3.8·10⁻⁴ (0.14 kPa at δ/H = 0.008) | 4·10⁻⁴ |
| Evaluations per increment: smooth / rough 2 × 2 × 2 / rough 3 × 3 × 3; bisections | 2.48 / 2.73 / 2.71; none | 3D 2.48 / 2.73; 2D 3 × 3: 2.97 |
| Material point, q(0.6) with 30 / 150 / 3000 increments | 149.262 / 149.503 / 149.553 | 149.262 / 149.503 / 149.553 |

Checks: the histories of the smooth and rough 2 × 2 × 2 runs (p', q at A, σ_a, largest p_w) coincide with those of
the Python code (`abaqus3d` of `gen_data3d.py`, `data_3d.pkl`) to 5·10⁻¹⁰ kPa (the precision of the CSV files),
with the same number of evaluations in every increment. With respect to the independent implementation of the
previous version (return mapping of de Souza Neto et al. with the unknowns Δγ and α), the largest differences of
q at A are 5.9·10⁻⁵ (smooth), 1.3·10⁻⁴ (rough 2 × 2 × 2) and 5.0·10⁻⁵ kPa (rough 3 × 3 × 3), and of σ_a 6·10⁻⁶ kPa
after the initial state. The smooth 3D values of q at A differ from the former axisymmetric ones by less than
3.2·10⁻⁴ kPa. With the rough platen the 3D model gives 0.38 kPa more than the axisymmetric 2 × 2 model at the end
(largest difference 0.47 kPa at δ/H = 0.19): the quarter of cylinder is not the axisymmetric mesh (different
discretization of the section), and A is not at the same position in the element (see below).

Full integration. The 3 × 3 × 3 rule locks as the 3 × 3 rule of the axisymmetric model did, and at about the same
stage. At A the full integration is first slightly stiffer than the reduced one (at most 1.60 kPa above it, at
δ/H = 0.22), falls below it from δ/H = 0.32 and by more than 1 kPa from 0.348, and q at A then reaches a maximum of
147.93 kPa at δ/H = 0.38 and decreases to 141.80 kPa (axisymmetric model of the previous version, same
definitions: 1.31 kPa at 0.228, 0.332, 0.364, maximum 149.0 kPa at 0.40, 145.5 kPa at the end; the first δ/H
with a difference larger than 1 kPa in either direction, 0.144 against 0.172, falls in the stiffer phase;
`abaqus_fe_numbers.csv`). Beyond δ/H ≈ 0.35 the stresses oscillate within the elements, while the average
axial stress on the platen is unaffected (250.93 against 250.78 kPa). At the end, q at the 27 points of the
element that contains A ranges from 142.9 to 167.6 kPa (mean 158.1; with 2 × 2 × 2 points 155.6 to 156.8, mean
155.9; columns `q_elA_min`, `q_elA_mean`, `q_elA_max`, band in Fig. 11b). The value at A depends on where A is
in the locked element: in the axisymmetric mesh A was the centroid of its element, a Gauss point of the 3 × 3
rule (145.5 kPa at the end, within a range of 137.7 to 168.8 kPa, mean 156.2, in that element), whereas in the
quarter of cylinder A lies on an edge of the element and its stress is extrapolated quadratically in two
directions (141.8 kPa, below the smallest value at the points). The lower end value of q at A in 3D (141.8
against 145.5 kPa) therefore does not mean a stronger locking: the oscillation within the element is of the same
size (25 kPa against 31 kPa).

Table 9, normalized residual r = ‖R‖/max(‖f_ext‖, 1) of each evaluation (rough platen, 2 × 2 × 2). ‖f_ext‖ < 1 kN
in both the quarter model and the axisymmetric model of the previous version, so that r is the norm of the
residual in kN and the tolerance 10⁻⁸ is absolute. The quarter model carries a quarter of the forces of the
axisymmetric model (which integrates over the whole ring, 2πr) spread over 321 instead of 37 nodes, and its
residuals are 6 to 10 times smaller in the same state (1.30·10⁻² against 1.30·10⁻¹ in the predictor of the first
increment, 8.3·10⁻⁵ against 5.4·10⁻⁴ in the 41st); the convergence pattern is the same:

| Iter. | D, 1st (δ/H = 0.004) | Dᵀ | D, 41st (0.164) | Dᵀ | D, 76th (0.304) | Dᵀ | D, 150th (0.600) | Dᵀ |
|---|---|---|---|---|---|---|---|---|
| 1 | 1.30e-02 | 1.30e-02 | 8.32e-05 | 8.32e-05 | 2.08e-05 | 2.08e-05 | 6.17e-06 | 6.17e-06 |
| 2 | 1.46e-03 | 2.22e-03 | 1.81e-07 | 3.26e-06 | 3.92e-08 | 5.11e-06 | 1.61e-09 | 1.19e-06 |
| 3 | 8.86e-05 | 4.72e-04 | 4.70e-12 | 1.79e-06 | 3.41e-13 | 3.39e-06 |  | 5.88e-07 |
| 4 | 9.84e-06 | 1.61e-05 |  | 8.54e-08 |  | 4.25e-07 |  | 2.04e-07 |
| 5 | 9.07e-09 | 1.96e-07 |  | 1.38e-07 |  | 9.70e-07 |  | 2.93e-07 |
| 6 |  | 2.68e-09 |  | 7.00e-09 |  | 1.53e-07 |  | 3.74e-08 |
| 7 |  |  |  |  |  | 2.12e-07 |  | 9.64e-08 |
| 8 |  |  |  |  |  | 9.31e-08 |  | 4.88e-08 |
| 9 |  |  |  |  |  | 3.10e-08 |  | 1.79e-08 |
| 10 |  |  |  |  |  | 3.19e-08 |  | 2.47e-08 |
| 11 |  |  |  |  |  | 1.35e-09 |  | 2.86e-09 |

Table 10, Abaqus columns (rough platen, 2 × 2 × 2, 150 increments, tolerance 10⁻⁸, no bisection with any operator;
time: processor time of the solution, `cpu_time_s`, one thread, Release build; the wall time, `wall_time_s`, was
0.2 to 0.6 s larger in this run and grows with the load of the machine; a previous run gave 14.9, 18.4, 34.8, 56.9
and 57.5 s):

| Operator | evaluations per increment (largest) | total (`NGlobalIterations`) | time (s) | final q_A (kPa) | previous version (axisymmetric, Python) |
|---|---|---|---|---|---|
| consistent D | 2.73 (5) | 409 | 16.0 | 155.937765 | 3.03 (5), 454, 19 s |
| central differences (`fd`) | 2.73 (5) | 409 | 18.1 | 155.937765 | 3.03 (5), 454, 234 s |
| symmetric part (D + Dᵀ)/2 | 5.77 (7) | 866 | 35.2 | 155.937729 | 7.11 (9), 1066, 43 s |
| continuum tangent | 9.22 (14) | 1383 | 57.4 | 155.937787 | 11.47 (16), 1720, 71 s |
| transpose Dᵀ | 9.18 (11) | 1377 | 58.1 | 155.937754 | 12.23 (16), 1835, 77 s |

The operators D and central differences need the same evaluations in every increment; the ratios to D are 2.12
(symmetric part), 3.38 (continuum) and 3.37 (transpose), against 2.35, 3.79 and 4.04 in the axisymmetric model.
The final values of q at A differ by at most 3.6·10⁻⁵ kPa, consistent with the tolerance of the global
iterations.

The smaller ratios are an effect of the tolerance, not of the three-dimensional model. The normalized residual
is the residual norm in kN (‖f_ext‖ < 1 kN, see Table 9), and the residuals of the quarter model are about ten
times smaller than those of the axisymmetric model, so that the tolerance 10⁻⁸ is about ten times looser here;
since the operators other than D converge linearly, they stop earlier. With the tolerance 10⁻⁹ (part
`tolerance`, `abaqus_table10_tolerance.csv`) the counts are those of the axisymmetric model with 10⁻⁸:

| Operator | evaluations per increment (largest) | total | ratio to D | time (s) | final q_A (kPa) |
|---|---|---|---|---|---|
| consistent D | 3.03 (6) | 455 | 1 | 17.1 | 155.937757 |
| central differences (`fd`) | 3.03 (6) | 455 | 1.00 | 20.6 | 155.937757 |
| symmetric part (D + Dᵀ)/2 | 7.38 (9) | 1107 | 2.43 | 45.1 | 155.937760 |
| continuum tangent | 11.68 (17) | 1752 | 3.85 | 73.5 | 155.937760 |
| transpose Dᵀ | 12.50 (16) | 1875 | 4.12 | 78.9 | 155.937757 |

Smooth platen with 600 increments (48 elements):

| Quantity | p'0 = 100 kPa | p'0 = 20 kPa | previous version (axisymmetric 2 × 4) |
|---|---|---|---|
| Peak q at A (δ/H) | 149.542 (0.6) | 54.614 (0.021) | 149.54 (0.6); 54.61 (0.021) |
| End p' and q at A | 149.848, 149.542 | 38.551, 38.706 | 149.85; 36.99, 37.18 |
| Global q = σ_a − p'0: peak, end | 149.529 | 54.610 (0.021), 33.702 | —; 54.61, 33.4 |
| Material point at the end | 149.543 | 30.091 | 149.54; 30.1 |
| Evaluations per increment (largest), bisections | 2.04 (6), none | 2.07 (5), none | 2.17; 2.09 |

Softening study, p'0 = 20 kPa, smooth platen, 600 increments (`abaqus_states_numbers.csv`,
`abaqus_softening_numbers.csv`):

| Quantity | 2 × 2 × 4 (48 elements) | 4 × 4 × 8 (384 elements) | previous version (axisymmetric 2 × 4 / 4 × 8) |
|---|---|---|---|
| Global q = σ_a − p'0 at δ/H = 0.6 (kPa) | 33.702 | 33.701 | 33.4 / 33.6 |
| Peak of the global q (δ/H) | 54.610 (0.021) | 54.613 (0.021) | 54.614 |
| q at A at 0.6 | 38.706 | 37.687 | 37.18 / 39.66 |
| Material point at 0.6 | 30.091 | 30.091 | 30.1 |
| Evaluations per increment (largest), bisections | 2.07 (5), none | 2.07 (6), none | 2.09 |
| Processor time of the solution (wall time) | 48.4 s (49.1 s) | 22.9 min (23.5 min; renumbered equations) | — |

The softening response is not homogeneous: the specimen bulges at the platen and at the mid-plane (outer radius
33.1 mm at both ends of the deformed half specimen, 24 mm high, with the coarse mesh and 33.4 mm with the refined
one, against 20.4 and 20.6 mm at mid-height; `abaqus_smooth_600_p0_20*_profile.csv`), and the two meshes give the
same global response (largest difference 0.45 kPa at δ/H = 0.34, 4·10⁻⁴ kPa at the end), 3.6 kPa above the
material point. The column `eps_v` of these runs is not a volumetric strain (the deformation is not
homogeneous).
