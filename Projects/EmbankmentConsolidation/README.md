# EmbankmentConsolidation — Sect. 6.6, Figs. 12–14, Table 8 and the embankment column of Table 10

Embankment loading on a Modified Cam-Clay foundation (FLAC3D example *Embankment loading on a Cam-Clay
foundation*): a 50 kPa strip load on a 10 m saturated clay layer, applied in undrained increments and
followed by consolidation up to t = 10⁸ s. This is the function `aterro()` of `gen_data.py`; the elastic
variant is `aterro_elastic.py` and the comparison of the tangent operators is the embankment part of
`tangentes()` of `gen_data.py`. Since v0.7 of the article the model is three-dimensional.

* Domain: the 1 m slice of half of the problem modelled by FLAC3D, the slab [0, 20] × [0, 10] × [0, 1] m
  (x horizontal, y vertical, z the thickness), with 20 × 10 × 1 Hex20–Hex8 elements of 1 m (serendipity
  quadratic displacement, trilinear pore pressure, i.e. pressure one order lower): 1553 displacement nodes
  (462 vertices and 1091 mid-edge nodes), 462 pore pressure nodes, 5121 equations, 3430 after the elimination
  of the Dirichlet conditions; 3 × 3 × 3 Gauss points (`SetIntegrationOrder(4)`), 5400 integration points
  (Fig. 12). The mesh is `mcc::CreateSlabMesh`; its CSV files are written with `mcc::WriteMeshCSV`.
* Material (Table 1): M = 0.888, λ = 0.161, κ = 0.062, v_λ = 2.858, porous elasticity with ν = 0.3,
  uniform p'c0 = 160 kPa. γ_sat = 23 kN/m³, γ_w = 10 kN/m³, K_f = 2·10⁵ kPa, n = 0.3 (α_B = 1,
  1/M_B = n/K_f), mobility k = 10⁻⁹ m²/(kPa s); body force (0, −23, 0) and fluid weight (0, −10, 0).
* Initial state at the depth d = 10 − y: p_w = 10 d, σ'_yy = −13 d, σ'_xx = σ'_zz = 0.7(−23 d) + 10 d = −6.1 d,
  p'c = 160 kPa and v0 = v_λ − λ ln p'c0 + κ ln(p'c0/p'0) at each integration point
  (`TPZMatPoroElastoPlasticUP::InitializeMemory`); hydrostatic nodal pore pressures (`mcc::SetInitialPressure`).
* Boundary conditions (boundary quadrilaterals of `CreateGeoMesh`): u_x = 0 at x = 0 and x = 20 m
  (`EDirichletUDirectional`), u_z = 0 on the faces z = 0 and z = 1 m (plane strain, `EDirichletUDirectional`),
  fixed and impermeable base (`EDirichletU`), drained top (`EDirichletP` on the whole top) and the strip load
  q = 50 kPa on 0 ≤ x ≤ 4 m (`ENeumannU` on coincident faces, scaled by the load factor).
* Loading: 10 undrained increments of the load factor (Δt = 0), then 25 consolidation steps,
  t = 10^(2 + j/4) s, j = 0…24; no predictor, no bisection needed.
* Monitoring: settlements of the top at x = 0, 2, 4, 6 m (vertices of the face z = 0) and the pore pressures
  pp1 and pp2 (means of the eight vertices of the elements centred at (0.5, 9.5, 0.5) and (1.5, 7.5, 0.5) m).

**Plane strain.** The solution of the 20 × 10 Q8–Q4 plane strain model of v0.6 does not depend on z, has
u_z = 0 and belongs to the Hex20–Hex8 space of the slab; with a z-independent stress and σ'_xz = σ'_yz = 0, the
integral over the thickness of each Hex20 test function is a combination of Q8 test functions with the same
boundary conditions (vertex: (N_c − Σ N_e)/6; mid-edge node of the edge along z: (2 N_c + Σ N_e)/3, the sums
over the Q8 mid-edge functions of the edges at the vertex c; mid-edge node of the faces z = 0, 1: N_e/2), and that
of each Hex8 test function is half of a Q4 one. The slab therefore has the 2D solution as its discrete solution:
the run checks it (differences between the faces z = 0 and z = 1 m below 3·10⁻¹⁶ m and 2·10⁻¹⁴ kPa,
|u_z| < 10⁻¹⁶ m at the mid-edge nodes z = 0.5 m, the three layers of integration points in the same state). What changes is the normalized residual ‖R‖/‖f_ext‖, whose
nodal components are those of another basis: for the same iterate the 3D value is 0.61 to 0.95 times the 2D one
(median 0.67), so that the third undrained increment, whose fourth residual is 1.04·10⁻⁸ in 2D and 6.8·10⁻⁹ in
3D, converges with one evaluation less (4 instead of 5), and the histories differ from the 2D ones at the level of
the Newton tolerance (5·10⁻¹⁰ m, 3·10⁻⁸ kPa).

Two models are solved (`EmbankmentConsolidation::EVariant`), and the Cam-Clay model is repeated with each tangent
operator of Table 10 (`TPZPlasticStepModifiedCamClay::SetTangentMode`):

| Run | Python | What it is |
|---|---|---|
| `ECamClay`, D | `aterro()` | the model above: undrained loading and consolidation, CSV and VTK files |
| `EElastic`, D | `aterro_elastic.py` | p_c = 10⁷ kPa in the state (no yielding), v0 still computed with p'c0 = 160 kPa; gives the plastic share of the final settlement |
| `ECamClay` with D, sym, cont, DT, fd | `tangentes()` (`'aterro'`) | Table 10: global iterations with the consistent tangent D, its symmetric part (D + Dᵀ)/2, the continuum tangent, the transpose Dᵀ and central differences of the stress update; the run with Dᵀ also gives the undrained loading with the transposed tangent of the end of `aterro()` (`logUt`) |

The class `EmbankmentConsolidation` follows the structure of the NeoPZ examples:

* `CreateGeoMesh` builds the slab and its boundary quadrilaterals with `mcc::CreateSlabMesh` and `mcc::FaceOnPlane`.
* `CreateCompMesh` builds the displacement (Hex20), pore pressure (Hex8) and multiphysics meshes. It also sets the
  `TPZMatPoroElastoPlasticUP` material (`EThreeDimensional`, tangent mode), the boundary conditions and the
  initial state.
* `Run` holds the analysis, `TPZPoroElastoPlasticUPAnalysis` with `TPZSkylineNSymStructMatrix` and
  `TPZStepSolver` LU. It applies the increments and does the post-processing (reactions, plane strain check,
  plastic points per layer, matrix profile, CSV and VTK files).
* `RunAll` runs the selected parts and prints the comparison with the 2D model of v0.6 and FLAC3D; the numbers
  are also written to `embankment_summary.csv`.

## Running

```
./EmbankmentConsolidation                     # all parts, CSV files and the VTK series of every converged state
./EmbankmentConsolidation novtk               # without the VTK series
./EmbankmentConsolidation model elastic       # only these parts (model, elastic, tangents)
./EmbankmentConsolidation tangents modes=D,DT # Table 10 with these operators only (D, sym, cont, DT, fd)
./EmbankmentConsolidation verbose             # one line per increment
```

Run times (Release build, one core, on a shared machine; `embankment_summary.csv`, `embankment_tangents.csv`): all
parts 16.4 min (983 s): the Cam-Clay model 177 s (45 s of which writing the CSV files and the VTK series), the
elastic model 158 s, and the five runs of the comparison of the tangent operators 121 s (D), 128 s (sym), 143 s
(cont), 131 s (DT) and 126 s (fd). Without any file output the Cam-Clay run takes 121 s (the run with D of the
comparison). The timings vary with the load of the machine (a first run of the same code took 1337 s in all,
36 % more). The 2D model of v0.6 (20 × 10 Q8–Q4, 1553 equations) took 16 s in NeoPZ and 20.5 s in Python.

Almost all of the time goes into the skyline LU decomposition: a global iteration takes 0.92–0.97 s, and a
separate measurement of one iteration gave 0.12 s for the assembly against 1.3 s for the LU and the substitutions
(about 90 %). The analysis does
not renumber the equations, so that the LU without pivoting is stable in the undrained steps: the pore pressure
equations (420 after the elimination) follow the 3010 displacement equations, and their columns reach back to the
first displacement equation of their elements. The upper profile has 1.82·10⁶ entries, 0.88·10⁶ in the
displacement columns (mean height 292, largest 336) and 0.94·10⁶ in the pore pressure columns (mean height 2249,
largest 3099); one LU costs about 1.1·10⁹ multiply-adds (`embankment_matrix_profile.csv`: height of each column of
the filtered system). The CSV files and the VTK series add about 45 s to the Cam-Clay run and as much to the elastic one.

Files written in the working directory:

| File | Contents |
|---|---|
| `embankment_mesh_{nodes,elements,faces,edges}.csv` | the geometric mesh (`mcc::WriteMeshCSV`): 462 vertices, 200 hexahedra, 464 boundary quadrilaterals with their material ids (−1 base, −2 x = 20 m, −4 x = 0, −5 loaded strip, −6 z = 0, −7 z = 1 m, −13 drained top) and 1091 edges (Fig. 12) |
| `embankment_monitor.csv` | settlement points (vertices) and the boxes of the elements of the zones pp1, pp2 (Fig. 12) |
| `embankment_history.csv` | t, λ, settlements at x = 0, 2, 4, 6 m, pp1, pp2 (Fig. 13, Table 8); 11 undrained states (t = 0) and 25 consolidation steps |
| `embankment_convergence.csv` | normalized residual of every evaluation of every increment (stage 0 undrained, 1 consolidation; Sect. 6.7) |
| `embankment_types.csv` | elastic, subcritical and supercritical integration points in each of the three layers z = const of points, end of the loading and t = 10⁸ s |
| `embankment_nodal.scal_vec.{0,1,2}.vtk` | nodal displacement and pore pressure (native NeoPZ VTK) at the end of the loading, t = 10⁶ s and t = 10⁸ s |
| `embankment_gauss_{undrained,t1e6,t1e8}.vtk` | integration points: p', q, p_c, v0, type (0 elastic, 1 subcritical, 2 supercritical), σ' |
| `embankment_nodal_{undrained,t1e6,t1e8}.csv` | vertices: x, y, z, u_x, u_y, u_z, p, excess pore pressure p − γ_w(10 − y) (Fig. 14a, b, d, face z = 0) |
| `embankment_gauss_{undrained,t1e6,t1e8}.csv` | integration points: x, y, z, p', q, p_c, v0, type (Fig. 14c) |
| `embankment_matrix_profile.csv` | height of each column of the skyline matrix of the filtered system and 1 for a pore pressure equation |
| `embankment_elastic_{history,convergence,types}.csv` | the same for the elastic variant |
| `embankment_tangents.csv` | Table 10, embankment column: per operator the mean and largest evaluations per undrained increment and per consolidation step, total evaluations, global iterations (`NGlobalIterations`, failed attempts included), failed attempts (`NBisections`), run time, time per iteration, final settlement and pore pressure, largest difference of the histories from D, largest relative change of the undrained residuals from D, and the 2D values of v0.6 |
| `embankment_tangents_iterations.csv` | evaluations of every increment of every operator |
| `embankment_tangents_convergence.csv` | normalized residual of every evaluation of every increment of every operator (mode, stage, increment, t, λ, evaluation, residual) |
| `embankment_summary.csv` | every number printed by the program (quantity, value, value of the 2D model of v0.6, unit) |

The FLAC3D markers of Fig. 13 are the digitized histories of the Python package
(`reference/flac_historicos_digitalizados.json`); they are not written here.

## Viewing the solution in ParaView

The Cam-Clay and the elastic runs write all their 36 converged states with `mcc::TVTKSeries`, in
`vtk/embankment/` and `vtk/embankment_elastic/` (prefix `embankment` or `embankment_elastic`; 112 files and
about 64 MB in each directory):

| File | Contents |
|---|---|
| `embankment_nodal.vtk.series` → `embankment_nodal.scal_vec.<k>.vtk` | `Displacement` (vector) and `PorePressure` at the vertices of the hexahedra, written by the NeoPZ graph mesh of the multiphysics mesh (`TPZAnalysis::DefineGraphMesh`/`PostProcess`) |
| `embankment_intpoints.vtk.series` → `embankment_intpoints.scal_vec.<k>.vtk` | the variables of the integration points of `TPZMatPoroElastoPlasticUP`, projected element by element on a discontinuous mesh by `TPZPostProcAnalysis` (as in the NeoPZ footing example): `MeanEffectiveStress` p', `DeviatoricStress` q, `PreconsolidationPressure` p_c, `PlasticType` (0 elastic, 1 subcritical, 2 supercritical), `VolumetricStrain` (tr ε, negative in compression), `SpecificVolume`, `EffectiveStressXX/YY/ZZ/XY/XZ/YZ`, `TotalStressXX/.../YZ` (σ' − p_w I), `PrincipalEffectiveStress` (vector σ'1 ≥ σ'2 ≥ σ'3 at the points) and the tensors `EffectiveStress` and `TotalStress` |
| `embankment_gausspoints.vtk.series` → `embankment_gausspoints.<k>.vtk` | the 5400 integration points as a point cloud with their exact values: p', q, p_c, v0, type and σ' |
| `embankment_states.csv` | index k of each state with its time t (s) and load factor λ |

The time of the series is the index k of the state: 0 is the geostatic state, 1 to 10 are the undrained
increments (λ = 0.1 … 1, t = 0) and 11 to 35 the consolidation steps (t = 10^(2 + (k − 11)/4) s, up to
10⁸ s); `embankment_states.csv` gives t and λ of each k. Stresses are positive in tension (as in NeoPZ);
p' and q follow the soil mechanics convention (p' > 0 in compression). With 3 × 3 × 3 Gauss points the
projection is triquadratic: its values at the element vertices are the Lagrange extrapolation of the 27
Gauss values (the interpolation of `mcc::StressAtPoint`), and the jumps between neighbouring elements show
the discretization error. Being extrapolations, the vertex values can overshoot: the projected `PlasticType`
is not an integer (it ranges from about −2.4 to 4.4 at t = 10⁸ s; the point cloud has the type of each point), q can
be slightly negative, and the projected principal stresses, sorted at the points, can be out of order where two
of them are close. For the principal values of the extrapolated tensor apply *Filters → Alphabetical →
Tensor Principal Invariants* (ParaView 5.10 or later) to `EffectiveStress`.

In ParaView (5.5 or later, which reads the `.series` files):

1. *File → Open* `vtk/embankment/embankment_nodal.vtk.series` (open the `.series` file, not the numbered
   files) and press *Apply*.
2. Choose the field in the *Coloring* box of the toolbar (`PorePressure`, or `Displacement` with its
   magnitude or a component) and use *Rescale to data range over all timesteps* for a fixed color scale.
3. Play the states with the VCR buttons of the *Time* toolbar (*Last Frame* is t = 10⁸ s).
4. Deformed shape: select the reader and apply *Filters → Alphabetical → Warp By Vector* with
   *Vectors* = `Displacement` and a *Scale Factor* of 10 to 20 (the settlements are below 0.3 m on a
   20 m model).
5. Excess pore pressure (Fig. 14a, b): select the reader (before the warp, so that the coordinates are the
   undeformed ones) and apply *Filters → Calculator* with *Result Array Name* `ExcessPorePressure` and the
   expression `PorePressure - 10*(10 - coordsY)`, i.e. p − γ_w (10 − y) kPa.
6. Integration point fields: open `vtk/embankment/embankment_intpoints.vtk.series` and color by
   `MeanEffectiveStress`, `DeviatoricStress`, `PreconsolidationPressure`, ... (the type of response of
   Fig. 14c is shown by the point cloud, item 7).
   To see them on the deformed mesh, apply *Filters → Resample With Dataset* (source: the nodal reader,
   destination: the integration point reader) to bring `Displacement` to its points, then *Warp By Vector*.
7. Point cloud: open `vtk/embankment/embankment_gausspoints.vtk.series`, set *Representation* to
   *Point Gaussian* (a *Gaussian Radius* of about 0.1 m) and color by `PlasticType`: at the last state the
   3 × 140 subcritical and 3 × 32 supercritical points of Fig. 14c appear, the same in the three layers of points.
8. The fields do not depend on z: the view along −z (*Set view direction to −Z*) shows the face z = 1 m, and
   *Filters → Alphabetical → Slice* with the normal (0, 0, 1) gives any section.

The elastic run (`vtk/embankment_elastic/`) has the same files. The VTK files of the three states of Fig. 14
listed above are still written in the working directory.

## Figures

```
python3 <neopz>/Projects/EmbankmentConsolidation/plot_figures.py [run directory] [-o output directory]
```

Run it after the executable, with the directory of its CSV files (default: the current directory). The figures
are written as PDF and PNG to `<run directory>/figures` (or to the output directory); Python 3 with numpy and
matplotlib is needed. A figure whose CSV files are missing is skipped with a message.

| File | Article | Data |
|---|---|---|
| `fig12_embankment_model` | Fig. 12: the slab in an oblique projection (x to the right, y up, z towards the reader: the faces z = 1 m, y = 10 m and x = 20 m are visible) with its 20 × 10 × 1 mesh, the load, the drained top, the boundary conditions, the settlement points and the zones pp1, pp2; on the right, the Hex20–Hex8 element with its 20 displacement nodes and 8 pore pressure nodes | `embankment_mesh_*.csv`, `embankment_monitor.csv` |
| `fig13_embankment_history` | Fig. 13: (a) settlements at x = 0, 2, 4, 6 m and (b) pore pressures pp1, pp2 from t = 0 to 10⁸ s, with the inset of pp2 against log t (Mandel–Cryer effect); markers: FLAC3D | `embankment_history.csv`, `reference/flac_historicos_digitalizados.json`; the numbers of the discussion read from the FLAC3D histories (share of the excess pore pressure of pp2 dissipated and of the consolidation settlement at x = 0 developed at t = 2.5·10⁵ s and 10⁶ s: 47 %, 71 %, 13 %), with those of this work, are written to `fig13_flac3d_numbers.csv` |
| `fig14_embankment_fields` | Fig. 14, on the face z = 0: excess pore pressure at the end of the undrained loading (a) and at t = 10⁶ s (b), plastic integration points of the layer z = 0.113 m (c) and settlement (d) at t = 10⁸ s | `embankment_nodal_{undrained,t1e6,t1e8}.csv`, `embankment_gauss_t1e8.csv` |

The script draws the three states of the article; the fields of all 36 converged states are in the VTK series
(see *Viewing the solution in ParaView*).

## Results

Table 8 (settlements in m, positive downwards; pore pressures in kPa). The 2D column is the 20 × 10 Q8–Q4 plane
strain model of v0.6 (Python `data_aterro.pkl`, reproduced to round-off by the 2D NeoPZ model); FLAC3D values read
from its histories:

| | End of loading: 3D | 2D v0.6 | FLAC3D | t = 10⁸ s: 3D | 2D v0.6 | FLAC3D |
|---|---|---|---|---|---|---|
| Settlement, x = 0 | 0.152818 | 0.152818 | 0.140 | 0.275113 | 0.275113 | 0.193 |
| Settlement, x = 2 m | 0.152329 | 0.152329 | 0.135 | 0.269305 | 0.269305 | 0.186 |
| Settlement, x = 4 m | 0.067338 | 0.067338 | 0.055 | 0.164231 | 0.164231 | 0.104 |
| Settlement, x = 6 m | −0.038550 | −0.038550 | −0.042 | 0.022668 | 0.022668 | 0.004 |
| Pore pressure, pp1 | 33.0561 | 33.0561 | 18.1 | 5.0224 | 5.0224 | 5.1 |
| Pore pressure, pp2 | 55.7250 | 55.7250 | 62.4 | 25.1075 | 25.1075 | 25.1 |

Over the 36 monitored states the histories differ from the 2D ones by at most 4.7·10⁻¹⁰ m (settlements) and
3.3·10⁻⁸ kPa (pore pressures), the level of the Newton tolerance (the third undrained increment stops one
iteration earlier in 3D, see the plane strain paragraph above); the vertex fields of Fig. 14 differ by at most
4.9·10⁻¹⁰ m and 6·10⁻⁸ kPa, the states of the integration points by 10⁻⁶ kPa, with the same type of response at
every point.

| Quantity | 3D (this code) | 2D v0.6 |
|---|---|---|
| Elements, displacement nodes, pore pressure nodes, equations | 200 Hex20–Hex8, 1553, 462, 5121 (3430 free) | 200 Q8–Q4, 661, 231, 1553 |
| Integration points | 5400 (3 layers of 1800) | 1800 |
| Initial residual, largest component at the free equations | 1.3·10⁻¹³ kN | 8.4·10⁻¹³ kN/m (article 8·10⁻¹³) |
| Initial vertical reaction of the base | 4600.000 kN (slab of 1 m) | 4600.000 kN/m |
| p'0, q0 at the base | 84.0, 69.0 kPa | 84, 69 |
| Difference from FLAC3D at the end of the loading: settlement x = 0, heave x = 6 m | 9.2 %, 7.3 % (of the FLAC3D history value −0.0416 m) | 9 %, 7 % |
| Reactions at t = 10⁸ s: x = 0, x = 20 m, base x, base y (kN) | 898.088, −795.210, −102.878, 4800.000 | 898.088, −795.210, −102.878, 4800.000 kN/m |
| Sum of the horizontal reactions | 9·10⁻⁷ kN | 9·10⁻⁷ kN/m (article < 10⁻⁵) |
| Out-of-plane reactions of the faces z = 0 and z = 1 m, all their nodes (plane strain force, −/+∫σ_zz dA, total stress) | 17165.275, −17165.275 kN; −∫σ_zz dV/T over the integration points 17165.275 kN | 17165.275 kN/m (−∫σ_zz dA from the stresses of `data_aterro.pkl`) |
| Plane strain check at t = 10⁸ s: \|u(z=0) − u(z=1)\|, \|p(z=0) − p(z=1)\|, \|u_z\| at z = 0.5 m | 2.7·10⁻¹⁶ m, 1.4·10⁻¹⁴ kPa, 8.4·10⁻¹⁷ m | — |
| Integration points [elastic, sub, super], end of loading | [5400, 0, 0] | [1800, 0, 0] |
| Integration points [elastic, sub, super], t = 10⁸ s | [4884, 420, 96]; [1628, 140, 32] in each layer z = 0.113, 0.5, 0.887 m | [1628, 140, 32] |
| Evaluations per increment, undrained | [5,5,4,4,4,4,4,4,4,4], mean 4.20 | [5,5,5,4,4,4,4,4,4,4], mean 4.30 |
| Evaluations per increment, consolidation | [2,2,2,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,5,5,5,5,4], mean 3.56 | identical |
| Global iterations, bisections | 131, 0 | 132, 0 |
| Residuals of the last undrained increment | 2.1e−2, 4.3e−4, 5.3e−7, 1.3e−12 | 2.2e−2, 6.6e−4, 8.1e−7, 2.0e−12 |
| Residuals of the step ending at t = 10⁶ s | 4.7e−5, 4.7e−3, 1.2e−5, 1.3e−10 | 7.7e−5, 5.6e−3, 1.7e−5, 2.0e−10 |
| Mandel–Cryer effect: largest pp2 | 58.840 kPa at t = 3.16·10⁵ s | 58.840 kPa at 3.16·10⁵ s |
| pp2 at t = 10⁶ s | 55.815 kPa | 55.815 |
| Settlement at x = 0, t = 10⁶ s, as a share of its consolidation part | 31.3 % | 31.3 % |
| Consolidation coefficient c_v = k(K + 4G/3) in zone pp2, end of loading and t = 10⁸ s | 1.24·10⁻⁶, 2.58·10⁻⁶ m²/s | 1.24·10⁻⁶, 2.58·10⁻⁶ |
| c_v t/d² at t = 2.5·10⁵ s, d = 2.5 m | 0.050–0.103 | 0.050–0.103 (article 0.04–0.1) |
| Largest excess pore pressure, end of loading | 56.14 kPa at (1, 9) on both faces | 56.14 at (1, 9) |
| Largest excess pore pressure, t = 10⁶ s | 33.42 kPa at (0, 8) | 33.42 at (0, 8) |
| Heave of the top at t = 10⁸ s | 5.11 mm, from x = 9 m | 5.11 mm from x = 9 m |
| Elastic variant, final settlements x = 0, 2, 4, 6 m | 0.267524, 0.262540, 0.159354, 0.020157 | identical (article 0.268 at x = 0) |
| Plastic share of the final settlement at x = 0 | 0.0076 m | 0.0076 (article 0.008) |

The residuals are the normalized residuals ‖R‖/max(‖f_ext‖, 1) of the free equations in the nodal basis; their
values change with the basis (3D), not their quadratic decrease. The c_v values are averages over the 27
integration points of the pp2 element (the nine of the 2D element in 2D), with K = v0 p'/κ and
G = 3K(1 − 2ν)/(2(1 + ν)).

Table 10, embankment column (`embankment_tangents.csv`): evaluations of the residual per consolidation step,
mean (largest) over the 25 steps; per undrained increment (all the points are elastic, and the elastic tangent is
symmetric, so that every operator gives the same undrained iterations; Dᵀ changes no residual of the undrained
stage); total evaluations of the run (undrained + consolidation; there are no failed attempts, so this is also
`NGlobalIterations`); run time of the C++ code (one core):

| Operator | Undrained | Consolidation | Total | Time (s) | 2D v0.6: consolidation, total, time (Python, s) |
|---|---|---|---|---|---|
| Consistent, D | 4.20 | 3.56 (5) | 131 | 121 | 3.56 (5), 132, 20.5 |
| Finite differences | 4.20 | 3.52 (5) | 130 | 126 | not run |
| Symmetric part, (D + Dᵀ)/2 | 4.20 | 3.76 (6) | 136 | 128 | 3.84 (6), 139, 23.5 |
| Continuum tangent | 4.20 | 4.36 (8) | 151 | 143 | 4.40 (8), 153, 26.4 |
| Transpose, Dᵀ | 4.20 | 3.88 (6) | 139 | 131 | 3.88 (6), 140, 23.6 |

The 3D counts are equal to or one less than the 2D ones in every step (sym: steps 18 and 22; cont: step 19), as
expected from the smaller normalized residual of the 3D model in steps that converge linearly. All the operators
reach the same solution: the histories differ from those with D by at most 7·10⁻¹⁰ m and
2·10⁻⁸ kPa.

The finite-difference operator (central differences with h = 10⁻⁷) gives the counts of D in every step but the
24th (4 against 5 evaluations), but, unlike in the Abaqus benchmark, it does not reproduce the residuals of D from
the 13th consolidation step on (`embankment_tangents_convergence.csv` against `embankment_convergence.csv`). The
consolidation steps have no predictor, so the first iteration of a step evaluates the operator at a zero strain
increment, and the integration points that yielded in the previous step (the first ones in step 12) are then
exactly on the yield surface, where the stress update is not differentiable. D is there the elastic operator (a
zero increment is an elastic step), whereas the central differences straddle the surface and return the mean of
the elastic and the elastoplastic operators (checked at a material point: D(0) = D_e to 4·10⁻⁶ kPa, and the
finite differences equal (D_e + D_ep)/2 to 10⁻⁶ relative). The first correction of the step is therefore different: the
second residual with the finite differences is 1–4 % smaller in steps 13 to 18, 13 % in step 19 and about 30 %
in steps 20 to 25, after which both converge quadratically. In step 24 this brings the fourth residual to
5.5·10⁻⁹, below the tolerance 10⁻⁸, against 5.9·10⁻⁸ with D. In steps 1 to 12 (no integration point on the yield surface at the
start of the step; the first 9 points yield in step 12) and in the undrained loading the residuals agree to
round-off: to 2·10⁻⁷ relative in every residual above 10⁻⁹; the larger relative differences, up to 6·10⁻⁵ in the
undrained stage (`maxreldiff_undrained_residuals_vs_D` of `embankment_tangents.csv`), are in residuals below 10⁻⁹.

The time per global iteration is about the same with every operator (column `seconds_per_iteration`: 0.92 s with
D, 0.94–0.95 s with sym, cont and DT, 0.97 s with fd): the LU of the skyline matrix dominates, and the twelve
extra stress updates per point and iteration of the finite differences add at most 6 % per iteration (−1 to +6 % in four runs on the
shared machine), so that the run times follow the iterations.

The program prints these comparisons; `embankment_summary.csv` lists every number with the value of the 2D model
of v0.6.
