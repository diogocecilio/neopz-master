# AbaqusTriaxialConsolidation — Sects. 6.4 and 6.6, Figs. 7–9, Tables 6 and 8

Abaqus benchmark 1.15.2: drained triaxial test with displacement control on a cylindrical clay specimen
(upper half: H = 60 mm, R = 20 mm). Modified Cam-Clay with porous elasticity: M = 1, λ = 0.174, κ = 0.026,
v0 = 2.08, ν = 0.3, p'c0 = 116.6 kPa, k = 2·10⁻¹⁰ m²/(kPa·s). Initial isotropic effective stress
p'0 = 100 kPa equal to the cell pressure; the platen moves down to δ/H = 0.6 in 400 days in 150 increments;
free drainage at the top. Smooth platen: only u_z prescribed; rough platen: also u_r = 0.

Parts (command line arguments, default `all`):

* `mp` — material point: drained tests from p'0 = 100 and 20 kPa with 600 increments and the closed form of
  Appendix B.1 (Fig. 8); the benchmark state with 30, 150 and 3000 increments (column Δδ/H = 0.02 of Table 6).
* `axi` — axisymmetric 2 × 4 Q8–Q4 models: smooth platen, rough platen with 2 × 2 (reduced) and 3 × 3 (full)
  integration (Fig. 9, Table 6); rough platen with the transposed operator Dᵀ (Sect. 6.6, Table 8).
* `states` — finite element solutions of the two initial states with 600 increments and the softening
  state with the 2 × 4 and 4 × 8 meshes (Sect. 6.4).
* `3d` — quarter of the specimen with 48 Hex20–Hex8 elements (quadratic geometry, `TPZQuadraticCube`,
  lateral nodes on the cylinder) and 2 × 2 × 2 points, smooth and rough platens (Fig. 9b, Table 6).

Monitored quantities: stress at point A (r = 5 mm, z = 7.5 mm) interpolated from the integration points of the
element that contains it, average axial stress on the platen (reaction divided by the area), largest excess pore
pressure and volumetric strain.

## Running

```
./AbaqusTriaxialConsolidation              # all parts, with the VTK series of every increment (about 3 min)
./AbaqusTriaxialConsolidation mp axi       # material point and axisymmetric models
./AbaqusTriaxialConsolidation novtk        # all parts, without the VTK series (about 1.5 min)
./AbaqusTriaxialConsolidation axi novtk    # axisymmetric models without the VTK series, about 5 s
```

Files: `abaqus_<run>.csv` (δ/H, p'_A, q_A, σ_a, max p_w, ε_v), `abaqus_<run>.scal_vec.0.vtk` and
`abaqus_<run>_gauss.vtk` (final nodal fields and integration points), `abaqus_material_point_*.csv`,
`abaqus_closed_form_*.csv` and `abaqus_table8.csv`; unless `novtk` is given, the VTK file series of every
increment of each finite element run in `vtk/<run>/` (see below).

## Viewing the solution in ParaView

Every finite element run writes all its converged increments (151 states for 150 increments, 601 for the
runs with 600 increments) with `mcc::TVTKSeries` in its own directory `vtk/<run>/`, where `<run>` is
`smooth_2x2`, `rough_2x2`, `rough_3x3`, `rough_2x2_transposed` (part `axi`), `smooth_600_p0_100`,
`smooth_600_p0_20`, `smooth_600_p0_20_mesh2x4`, `smooth_600_p0_20_mesh4x8` (part `states`), `3d_smooth` and
`3d_rough` (part `3d`). The file names start with `abaqus_<run>`. All the series take about 160 MB in
10 000 files (about 40 MB for each 3D run and for the 4 × 8 mesh, 3 to 11 MB for the other 2D runs) and add
about 2 min to the run (0.01 s per state for the 2 × 4 mesh and 0.2 s for the 3D mesh, mostly the NeoPZ graph
mesh of the integration point fields):

| File | Contents |
|---|---|
| `abaqus_<run>_nodal.vtk.series` → `abaqus_<run>_nodal.scal_vec.<k>.vtk` | `Displacement` (vector) and `PorePressure` (excess pore pressure, zero at the start) at the nodes, written by the NeoPZ graph mesh of the multiphysics mesh |
| `abaqus_<run>_intpoints.vtk.series` → `abaqus_<run>_intpoints.scal_vec.<k>.vtk` | the variables of the integration points of `TPZMatPoroElastoPlasticUP`, projected element by element on a discontinuous mesh by `TPZPostProcAnalysis`: `MeanEffectiveStress` p', `DeviatoricStress` q, `PreconsolidationPressure` p_c, `PlasticType` (0 elastic, 1 subcritical, 2 supercritical), `VolumetricStrain` (tr ε, negative in compression), `SpecificVolume`, the components `EffectiveStressXX/YY/ZZ/XY` and `TotalStressXX/YY/ZZ/XY` (also `XZ`, `YZ` in 3D), `PrincipalEffectiveStress` (vector σ'1 ≥ σ'2 ≥ σ'3 at the points) and the tensors `EffectiveStress` and `TotalStress` |
| `abaqus_<run>_gausspoints.vtk.series` → `abaqus_<run>_gausspoints.<k>.vtk` | the integration points as a point cloud with their exact values: p', q, p_c, v0, type and σ' |
| `abaqus_<run>_states.csv` | index k of each increment with δ/H, the time t (s) and the platen displacement u_c (m) |

The time of the series is δ/H (0 to 0.6). In the axisymmetric runs x is the radius r and y the height z
(the components XX, YY, ZZ and XY of the stresses are rr, zz, θθ and rz). Stresses are positive in tension;
p' and q follow the soil mechanics convention (p' > 0 in compression). The projection has the order n − 1
for n × n (× n) Gauss points: bilinear (trilinear in 3D) with the reduced 2 × 2 (2 × 2 × 2) rule and
biquadratic with the full 3 × 3 rule, so that its values at the element vertices are the Lagrange
extrapolation of the Gauss values (the interpolation used for the point A, `mcc::StressAtPoint`). Being
extrapolations, the vertex values can overshoot: the projected `PlasticType` is not an integer (it ranges
from about −1.2 to 3.3 in `rough_3x3`; the point cloud has the type of each point), q can be slightly
negative, and the projected principal stresses, sorted at the points, can be out of order where two of them
are close (by up to 10 kPa near the axis in the softening runs). For the principal values of the
extrapolated tensor apply *Filters → Alphabetical → Tensor Principal Invariants* (ParaView 5.10 or later) to
`EffectiveStress`.

In ParaView (5.5 or later, which reads the `.series` files):

1. *File → Open* `vtk/rough_2x2/abaqus_rough_2x2_nodal.vtk.series` (the `.series` file, not the numbered
   files) and press *Apply*.
2. Choose the field in the *Coloring* box of the toolbar and use *Rescale to data range over all
   timesteps* for a fixed color scale.
3. Play the increments with the VCR buttons of the *Time* toolbar (the time shown is δ/H).
4. Deformed shape: *Filters → Alphabetical → Warp By Vector* with *Vectors* = `Displacement` and
   *Scale Factor* 1 (the true shape: the platen moves down 36 mm on the 60 mm half specimen).
5. Integration point fields: open `vtk/rough_2x2/abaqus_rough_2x2_intpoints.vtk.series` and color by
   `DeviatoricStress`, `MeanEffectiveStress`, `PreconsolidationPressure`, ... (the type of response is
   shown by the point cloud, item 6). To see them on the deformed mesh, apply *Filters → Resample With
   Dataset* (source: the nodal reader, destination: the integration point reader) and then *Warp By Vector*.
6. Point cloud: open `abaqus_rough_2x2_gausspoints.vtk.series`, set *Representation* to *Point Gaussian*
   (a *Gaussian Radius* of about 0.5 mm) and color by `PlasticType` or `DeviatoricStress`.
7. Whole specimen from the axisymmetric section: *Transform* (rotate 90° about X, so that the axis y
   becomes z), *Extract Surface* and *Rotational Extrusion* (angle 360°). From the 3D quarter: *Reflect*
   with *Plane* X and then with *Plane* Y, keeping the input (*Copy Input*); a further *Reflect* with
   *Plane* Z gives the lower half.

The pore pressure is the excess pore pressure (the initial value is zero and gravity is not modelled).
The test is drained: the excess pore pressure is largest in the first increments (0.14 to 0.15 kPa at
δ/H ≤ 0.01 with p'0 = 100 kPa, 0.07 kPa at δ/H ≈ 0.02 with p'0 = 20 kPa) and decreases to about
4·10⁻⁴ kPa at the end.

## Figures

```
python3 <neopz>/Projects/AbaqusTriaxialConsolidation/plot_figures.py [run directory] [-o output directory]
```

Run it after the executable, with the directory of its CSV files (default: the current directory). The figures
are written as PDF and PNG to `<run directory>/figures` (or to the output directory); Python 3 with numpy and
matplotlib is needed. A figure whose CSV files are missing is skipped with a message.

The figures use only the CSV files (`./AbaqusTriaxialConsolidation mp axi 3d novtk` is enough):

| File | Article | Data (part of the executable) |
|---|---|---|
| `fig07_abaqus_model` | Fig. 7: (a) axisymmetric 2 × 4 Q8–Q4 mesh with the boundary conditions and the point A; (b) quarter of the specimen with 48 Hex20–Hex8 elements | none: the meshes are rebuilt with the generators of the C++ code |
| `fig08_abaqus_states` | Fig. 8: material point from p'0 = 100 and 20 kPa against the closed form: p'–q, q–ε1 and εv–ε1 | `abaqus_material_point_p0_*.csv`, `abaqus_closed_form_p0_*.csv` (`mp`) |
| `fig09_abaqus_results` | Fig. 9: q at A against δ/H with the smooth (a) and rough (b) platens and the stress paths at A (c), with the digitized Abaqus curves | `abaqus_smooth_2x2.csv`, `abaqus_rough_2x2.csv`, `abaqus_rough_3x3.csv` (`axi`), `abaqus_3d_rough.csv` (`3d`), `abaqus_material_point_30.csv` (`mp`), `reference/abaqus_1_15_2_digitalizado.json` |

## Results

Table 6, q at A (kPa); this code = article ("this work") in all entries:

| δ/H | Abaqus smooth | smooth | material point Δδ/H = 0.02 | Abaqus rough | rough 2 × 2 | rough 3 × 3 | 3D rough |
|---|---|---|---|---|---|---|---|
| 0.03 | 59.7 | 60.2 | 56.7 | 59.7 | 61.6 | 61.6 | 61.5 |
| 0.13 | 109.9 | 114.0 | 110.1 | 111.8 | 116.8 | 117.4 | 116.4 |
| 0.22 | 131.1 | 133.8 | 131.0 | 134.7 | 137.4 | 138.7 | 136.9 |
| 0.32 | 141.1 | 143.5 | 141.9 | 145.0 | 147.2 | 147.4 | 146.9 |
| 0.41 | 145.8 | 147.1 | 146.2 | 150.1 | 151.6 | 149.0 | 151.6 |
| 0.51 | 148.1 | 148.8 | 148.4 | 153.1 | 154.3 | 147.7 | 154.5 |
| 0.60 | — | 149.5 | 149.3 | — | 155.6 | 145.5 | 155.9 |

| Quantity | this code | Python / article |
|---|---|---|
| End q_A smooth / rough 2 × 2 / rough 3 × 3 (kPa) | 149.502 / 155.559 / 145.546 | 149.50 / 155.56 / 145.5 |
| Maximum q_A, full integration | 149.0 at δ/H = 0.40 | 149.0 at 0.40 |
| Average axial stress on the platen, rough 2 × 2 / 3 × 3 (kPa) | 250.831 / 251.008 | 250.8 / 251.0 |
| Largest excess pore pressure at the end, rough 2 × 2 (kPa) | 3.8·10⁻⁴ | 4·10⁻⁴ |
| Evaluations per increment: smooth / rough 2 × 2 / rough 3 × 3 | 2.653 / 3.027 / 2.973 | 2.65 / 3.03 / 2.97 |
| Evaluations per increment with Dᵀ (rough 2 × 2) | 12.233 | 12.2 |
| 3D: nodes, pressure nodes | 321, 95 | 321, 95 |
| 3D end q_A smooth / rough (kPa), evaluations | 149.502 / 155.938, 2.48 / 2.727 | 149.50 / 155.94, 2.48 / 2.73 |
| Material point, q(0.6) with 30 / 150 / 3000 increments | 149.262 / 149.503 / 149.553 | 149.262 / 149.503 / 149.553 |
| FE, p'0 = 20 kPa: peak q_A, end p' / q / ε_v | 54.614 at 0.021, 36.994 / 37.179 / −0.681 | 54.61 at 0.021, 36.99 / 37.18 / −0.681 |
| Softening, global q from the platen: mesh 2 × 4 / 4 × 8 at the end | 33.425 / 33.635 | 33.4 / 33.6 (material point 30.1) |

Table 8 (normalized residual of the iterations, rough platen) is reproduced exactly, e.g. first increment
with D: 1.30e−1, 1.27e−2, 7.37e−4, 1.09e−5, 3.87e−9, and with Dᵀ: 1.30e−1, 1.85e−2, 2.52e−3, 9.02e−5,
5.39e−7, 1.83e−8, 1.47e−10.
